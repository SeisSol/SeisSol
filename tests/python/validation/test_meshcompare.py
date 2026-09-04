# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

"""Tests for postprocessing/validation/meshcompare.py

meshcompare.compare() is monolithic — it:
  1. Opens two XDMF files via seissolxdmf
  2. Matches the cells of the two files geometrically, independently of the cell
     order and of the vertex order within a cell; falls back to a per-element
     comparison via global-id if they are not the same cells
  3. Computes L2 errors per quantity
  4. Has a hardcoded workaround for a known SeisSol bug: DS field zeros
  5. Calls sys.exit(1) on threshold violation

To test without building real HDF5/XDMF files, we monkey-patch
seissolxdmf.seissolxdmf with a test double backed by numpy arrays.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

SEISSOL_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(SEISSOL_ROOT / "postprocessing" / "validation"))

import meshcompare  # noqa: E402


class FakeSeissolXdmf:
    """Test double for seissolxdmf.seissolxdmf.

    Backed by a dict registry keyed by "file path" (any string). Supports
    the subset of methods meshcompare.compare actually calls.
    """

    _registry: dict = {}

    def __init__(self, path):
        data = self._registry[path]
        self.geom = np.asarray(data["geom"], dtype=float)
        self.connect = np.asarray(data["connect"], dtype=int)
        self.fields = {
            k: np.asarray(v, dtype=float) for k, v in data.get("fields", {}).items()
        }
        self.int_fields = {
            k: np.asarray(v, dtype=int) for k, v in data.get("int_fields", {}).items()
        }
        self.nElements = self.connect.shape[0]
        self.ndt = data.get("ndt", 1)

    # API used by meshcompare.compare()
    def ReadGeometry(self):
        return self.geom

    def ReadConnect(self):
        return self.connect

    def ReadAvailableDataFields(self):
        return list(self.fields.keys()) + list(self.int_fields.keys())

    def Read1dData(self, name, n, isInt=False):
        if isInt:
            return self.int_fields[name]
        return self.fields[name]

    def ReadData(self, name, index):
        # meshcompare always asks for the last time index
        arr = self.fields[name]
        if arr.ndim == 1:
            return arr
        return arr[index]


@pytest.fixture
def patch_seissolxdmf(monkeypatch):
    """Patch meshcompare's seissolxdmf symbol and reset the registry."""
    FakeSeissolXdmf._registry = {}
    monkeypatch.setattr(meshcompare.sx, "seissolxdmf", FakeSeissolXdmf)
    return FakeSeissolXdmf._registry


# Standard 2-tetrahedron mesh: one shared face, 5 distinct vertices
# Cell 0: (0, 1, 2, 3)  barycenter ≈ (0.25, 0.25, 0.25)
# Cell 1: (0, 1, 2, 4)  barycenter ≈ (0.25, 0.25, -0.25)
GEOM = np.array(
    [
        [0.0, 0.0, 0.0],  # 0
        [1.0, 0.0, 0.0],  # 1
        [0.0, 1.0, 0.0],  # 2
        [0.0, 0.0, 1.0],  # 3
        [0.0, 0.0, -1.0],  # 4
    ]
)
CONNECT = np.array(
    [
        [0, 1, 2, 3],
        [0, 1, 2, 4],
    ]
)


class TestMeshCompareExact:
    """Two identical meshes with identical data should pass."""

    def test_identical_meshes_and_data_pass(self, patch_seissolxdmf, capsys):
        data = {
            "geom": GEOM,
            "connect": CONNECT,
            "fields": {"v1": np.array([1.0, 2.0])},
            "int_fields": {"global-id": np.array([0, 1])},
        }
        patch_seissolxdmf["sim.xdmf"] = data
        patch_seissolxdmf["ref.xdmf"] = data

        # Should complete without sys.exit
        meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=0.01)

    def test_identical_without_global_id_still_matches(self, patch_seissolxdmf, capsys):
        data = {
            "geom": GEOM,
            "connect": CONNECT,
            "fields": {"v1": np.array([1.0, 2.0])},
        }
        patch_seissolxdmf["sim.xdmf"] = data
        patch_seissolxdmf["ref.xdmf"] = data

        meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=0.01)
        out = capsys.readouterr().out
        assert "Matched all 2 cells geometrically" in out


class TestMeshCompareFailures:
    """Violations should trigger sys.exit(1)."""

    def test_large_field_difference_exits(self, patch_seissolxdmf):
        patch_seissolxdmf["sim.xdmf"] = {
            "geom": GEOM,
            "connect": CONNECT,
            "fields": {"v1": np.array([1.0, 2.0])},
            "int_fields": {"global-id": np.array([0, 1])},
        }
        patch_seissolxdmf["ref.xdmf"] = {
            "geom": GEOM,
            "connect": CONNECT,
            "fields": {"v1": np.array([10.0, 20.0])},  # 10x too big
            "int_fields": {"global-id": np.array([0, 1])},
        }
        with pytest.raises(SystemExit) as exc_info:
            meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=0.01)
        assert exc_info.value.code == 1

    def test_geometry_mismatch_exits(self, patch_seissolxdmf):
        """Cells that are not in the same place must not be compared."""
        shifted_geom = GEOM + 1.0  # different geometry
        patch_seissolxdmf["sim.xdmf"] = {
            "geom": GEOM,
            "connect": CONNECT,
            "fields": {"v1": np.array([1.0, 2.0])},
            "int_fields": {"global-id": np.array([0, 1])},
        }
        patch_seissolxdmf["ref.xdmf"] = {
            "geom": shifted_geom,
            "connect": CONNECT,
            "fields": {"v1": np.array([1.0, 2.0])},
            "int_fields": {"global-id": np.array([0, 1])},
        }
        with pytest.raises(SystemExit) as exc_info:
            meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=0.01)
        assert exc_info.value.code == 1

    def test_mismatched_global_ids_fail(self, patch_seissolxdmf):
        """If the cells do not align geometrically and cannot be aggregated into
        matching elements either, the comparison gives up."""
        patch_seissolxdmf["sim.xdmf"] = {
            "geom": GEOM,
            "connect": CONNECT,
            "fields": {"v1": np.array([1.0, 2.0])},
            "int_fields": {"global-id": np.array([0, 1])},
        }
        # Reference has the SAME global-ids but genuinely different geometry
        # for the cell indexed 0 (one extra vertex at (99, 99, 99))
        ref_geom = np.array(
            [
                [99.0, 99.0, 99.0],  # 0 — completely different location
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, -1.0],
            ]
        )
        patch_seissolxdmf["ref.xdmf"] = {
            "geom": ref_geom,
            "connect": CONNECT,
            "fields": {"v1": np.array([1.0, 2.0])},
            "int_fields": {"global-id": np.array([0, 1])},
        }
        with pytest.raises(SystemExit) as exc_info:
            meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=0.01)
        assert exc_info.value.code == 1


class TestMeshCompareMetafieldsExcluded:
    """Fields like 'partition', 'fault-tag', 'clustering' must be skipped."""

    def test_partition_field_not_compared(self, patch_seissolxdmf, capsys):
        # Reference has a bogus partition field; sim has a different one
        # — must not cause failure since 'partition' is in the ignore list
        patch_seissolxdmf["sim.xdmf"] = {
            "geom": GEOM,
            "connect": CONNECT,
            "fields": {
                "v1": np.array([1.0, 2.0]),
                "partition": np.array([0.0, 0.0]),  # different
            },
            "int_fields": {"global-id": np.array([0, 1])},
        }
        patch_seissolxdmf["ref.xdmf"] = {
            "geom": GEOM,
            "connect": CONNECT,
            "fields": {
                "v1": np.array([1.0, 2.0]),
                "partition": np.array([99.0, 99.0]),  # very different
            },
            "int_fields": {"global-id": np.array([0, 1])},
        }
        meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=0.01)

    def test_fault_tag_field_not_compared(self, patch_seissolxdmf):
        patch_seissolxdmf["sim.xdmf"] = {
            "geom": GEOM,
            "connect": CONNECT,
            "fields": {
                "v1": np.array([1.0, 2.0]),
                "fault-tag": np.array([1.0, 2.0]),
            },
            "int_fields": {"global-id": np.array([0, 1])},
        }
        patch_seissolxdmf["ref.xdmf"] = {
            "geom": GEOM,
            "connect": CONNECT,
            "fields": {
                "v1": np.array([1.0, 2.0]),
                "fault-tag": np.array([77.0, 88.0]),
            },
            "int_fields": {"global-id": np.array([0, 1])},
        }
        meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=0.01)

    def test_clustering_field_not_compared(self, patch_seissolxdmf):
        patch_seissolxdmf["sim.xdmf"] = {
            "geom": GEOM,
            "connect": CONNECT,
            "fields": {
                "v1": np.array([1.0, 2.0]),
                "clustering": np.array([0.0, 1.0]),
            },
            "int_fields": {"global-id": np.array([0, 1])},
        }
        patch_seissolxdmf["ref.xdmf"] = {
            "geom": GEOM,
            "connect": CONNECT,
            "fields": {
                "v1": np.array([1.0, 2.0]),
                "clustering": np.array([7.0, 8.0]),
            },
            "int_fields": {"global-id": np.array([0, 1])},
        }
        meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=0.01)


class TestMeshCompareDSBugWorkaround:
    """meshcompare has an inline comment and workaround for the DS field:
        # There is a bug on the master branch, which sets DS output to zero
        # in wrong places.
    When the REFERENCE is ~zero, the simulation value is forced to zero
    before the comparison — so any simulation DS value passes.
    """

    def test_ds_zero_in_reference_masks_sim_value(self, patch_seissolxdmf):
        # Reference DS = 0 everywhere → sim DS gets zeroed → passes
        patch_seissolxdmf["sim.xdmf"] = {
            "geom": GEOM,
            "connect": CONNECT,
            "fields": {"DS": np.array([1e5, 1e5])},  # huge sim DS values
            "int_fields": {"global-id": np.array([0, 1])},
        }
        patch_seissolxdmf["ref.xdmf"] = {
            "geom": GEOM,
            "connect": CONNECT,
            "fields": {"DS": np.array([0.0, 0.0])},
            "int_fields": {"global-id": np.array([0, 1])},
        }
        # The workaround zeroes the sim where ref is ~0, so this passes
        # despite the gigantic raw difference
        meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=0.01)


class TestMeshCompareVelocityRenaming:
    """v1/v2/v3 should map to u/v/w in the reference if only the legacy
    names are present there."""

    def test_v1_in_sim_falls_back_to_u_in_ref(self, patch_seissolxdmf):
        patch_seissolxdmf["sim.xdmf"] = {
            "geom": GEOM,
            "connect": CONNECT,
            "fields": {"v1": np.array([1.0, 2.0])},
            "int_fields": {"global-id": np.array([0, 1])},
        }
        # ref has the LEGACY name "u" instead of "v1"
        patch_seissolxdmf["ref.xdmf"] = {
            "geom": GEOM,
            "connect": CONNECT,
            "fields": {"u": np.array([1.0, 2.0])},
            "int_fields": {"global-id": np.array([0, 1])},
        }
        # meshcompare should use u from ref when comparing v1
        meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=0.01)


# A single tetrahedron split into four subcells by its barycenter, as the volume
# output does for `refinement = 1`. All four carry the same global-id.
REFINED_GEOM = np.array(
    [
        [0.0, 0.0, 0.0],  # 0
        [1.0, 0.0, 0.0],  # 1
        [0.0, 1.0, 0.0],  # 2
        [0.0, 0.0, 1.0],  # 3
        [0.25, 0.25, 0.25],  # 4 — barycenter
    ]
)
REFINED_CONNECT = np.array(
    [
        [0, 1, 2, 4],
        [0, 1, 4, 3],
        [0, 2, 3, 4],
        [1, 2, 4, 3],
    ]
)
REFINED_VALUES = np.array([1.0, 2.0, 3.0, 4.0])
REFINED_IDS = np.array([7, 7, 7, 7])


class TestMeshCompareReordering:
    """The whole point of the geometric matching: neither the cell order nor the
    vertex order within a cell may influence the result."""

    def test_shuffled_cells_still_match(self, patch_seissolxdmf, capsys):
        order = [2, 0, 3, 1]
        patch_seissolxdmf["sim.xdmf"] = {
            "geom": REFINED_GEOM,
            "connect": REFINED_CONNECT[order],
            "fields": {"v1": REFINED_VALUES[order]},
            "int_fields": {"global-id": REFINED_IDS},
        }
        patch_seissolxdmf["ref.xdmf"] = {
            "geom": REFINED_GEOM,
            "connect": REFINED_CONNECT,
            "fields": {"v1": REFINED_VALUES},
            "int_fields": {"global-id": REFINED_IDS},
        }

        meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=1e-12)
        out = capsys.readouterr().out
        assert "Matched all 4 cells geometrically" in out
        assert "4 of them reordered" in out or "reordered" in out

    def test_reoriented_cells_still_match(self, patch_seissolxdmf, capsys):
        # swap the last two vertices of every cell, i.e. flip the orientation
        reoriented = REFINED_CONNECT[:, [0, 1, 3, 2]]
        patch_seissolxdmf["sim.xdmf"] = {
            "geom": REFINED_GEOM,
            "connect": reoriented,
            "fields": {"v1": REFINED_VALUES},
            "int_fields": {"global-id": REFINED_IDS},
        }
        patch_seissolxdmf["ref.xdmf"] = {
            "geom": REFINED_GEOM,
            "connect": REFINED_CONNECT,
            "fields": {"v1": REFINED_VALUES},
            "int_fields": {"global-id": REFINED_IDS},
        }

        meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=1e-12)
        assert "Matched all 4 cells geometrically" in capsys.readouterr().out

    def test_permuted_values_are_reported_not_hidden(self, patch_seissolxdmf, capsys):
        """Values attached to the wrong subcell must fail, even though the cells
        themselves match perfectly."""
        patch_seissolxdmf["sim.xdmf"] = {
            "geom": REFINED_GEOM,
            "connect": REFINED_CONNECT,
            "fields": {"v1": REFINED_VALUES[[1, 2, 0, 3]]},
            "int_fields": {"global-id": REFINED_IDS},
        }
        patch_seissolxdmf["ref.xdmf"] = {
            "geom": REFINED_GEOM,
            "connect": REFINED_CONNECT,
            "fields": {"v1": REFINED_VALUES},
            "int_fields": {"global-id": REFINED_IDS},
        }

        with pytest.raises(SystemExit) as exc_info:
            meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=0.01)
        assert exc_info.value.code == 1
        assert "Matched all 4 cells geometrically" in capsys.readouterr().out


class TestMeshCompareAggregation:
    """Different subdivisions of the same element fall back to a per-element
    comparison instead of failing outright."""

    def test_different_tiling_aggregates(self, patch_seissolxdmf, capsys):
        # the reference does not refine at all; the simulation splits by 4
        patch_seissolxdmf["sim.xdmf"] = {
            "geom": REFINED_GEOM,
            "connect": REFINED_CONNECT,
            # the four subcells have equal volume, so the mean is 2.5
            "fields": {"v1": REFINED_VALUES},
            "int_fields": {"global-id": REFINED_IDS},
        }
        patch_seissolxdmf["ref.xdmf"] = {
            "geom": REFINED_GEOM[:4],
            "connect": np.array([[0, 1, 2, 3]]),
            "fields": {"v1": np.array([2.5])},
            "int_fields": {"global-id": np.array([7])},
        }

        meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=1e-12)
        out = capsys.readouterr().out
        assert "Falling back to a per-element comparison" in out
        assert "Aggregated 4 cells into 1 elements" in out

    def test_different_tiling_with_wrong_mean_fails(self, patch_seissolxdmf):
        patch_seissolxdmf["sim.xdmf"] = {
            "geom": REFINED_GEOM,
            "connect": REFINED_CONNECT,
            "fields": {"v1": REFINED_VALUES},
            "int_fields": {"global-id": REFINED_IDS},
        }
        patch_seissolxdmf["ref.xdmf"] = {
            "geom": REFINED_GEOM[:4],
            "connect": np.array([[0, 1, 2, 3]]),
            "fields": {"v1": np.array([25.0])},  # 10x the actual mean
            "int_fields": {"global-id": np.array([7])},
        }

        with pytest.raises(SystemExit) as exc_info:
            meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=0.01)
        assert exc_info.value.code == 1
