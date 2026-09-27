# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

"""Tests for scripts/validate/meshcompare.py

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
sys.path.insert(0, str(SEISSOL_ROOT / "scripts" / "validate"))

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


class TestMeshCompareComponentGroups:
    """A component is measured against the largest norm of its vector or tensor."""

    @staticmethod
    def errors(registry, tmp_path, sim, ref, epsilon=1e-3):
        for name, fields in (("sim.xdmf", sim), ("ref.xdmf", ref)):
            registry[name] = {"geom": GEOM, "connect": CONNECT, "fields": fields}
        report = tmp_path / "report.json"
        meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=epsilon, report_json=report)
        import json

        return json.loads(report.read_text())["quantities"]

    def test_small_component_is_relative_to_its_vector(
        self, patch_seissolxdmf, tmp_path
    ):
        # v2 is nearly zero, as the fault-normal traction change of a symmetric setup; its
        # noise relative to itself would be 1e-2, relative to the velocity 1e-8 or so
        ref = {
            "v1": np.array([1.0, 2.0]),
            "v2": np.array([1e-3, -1e-3]),
            "v3": np.array([0.5, 0.5]),
        }
        sim = dict(ref, v2=ref["v2"] + 1e-4)
        errors = self.errors(patch_seissolxdmf, tmp_path, sim, ref)
        assert errors["v2"] == pytest.approx(2 * 1e-8 / (1.0 + 4.0))
        assert errors["v1"] == 0.0

    def test_other_vectors_do_not_scale_a_component(self, patch_seissolxdmf):
        ref = {
            "v2": np.array([1e-3, -1e-3]),
            "u1": np.array([1e3, 1e3]),
        }
        sim = dict(ref, v2=ref["v2"] + 1e-4)
        for name, fields in (("sim.xdmf", sim), ("ref.xdmf", ref)):
            patch_seissolxdmf[name] = {
                "geom": GEOM,
                "connect": CONNECT,
                "fields": fields,
            }
        with pytest.raises(SystemExit):
            meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=1e-3)

    def test_zero_component_is_relative_to_its_vector(
        self, patch_seissolxdmf, tmp_path
    ):
        # identically zero in the reference, as the strike slip of a pure dip-slip source;
        # compared on its own it would be an absolute error weighted by the cell volumes
        ref = {"Sls": np.array([0.0, 0.0]), "Sld": np.array([2.0, 2.0])}
        sim = dict(ref, Sls=np.array([1e-6, -1e-6]))
        errors = self.errors(patch_seissolxdmf, tmp_path, sim, ref)
        assert errors["Sls"] == pytest.approx(1e-12 / 4.0)

    def test_fused_simulations_are_apart(self, patch_seissolxdmf, tmp_path):
        ref = {
            "v1-1": np.array([1.0, 1.0]),
            "v2-1": np.array([0.0, 0.0]),
            "v1-2": np.array([1e-3, 1e-3]),
            "v2-2": np.array([0.0, 0.0]),
        }
        sim = dict(ref, **{"v2-2": np.array([1e-4, 1e-4])})
        errors = self.errors(patch_seissolxdmf, tmp_path, sim, ref, epsilon=1.0)
        # relative to the velocity of its own simulation, not to that of the first
        assert errors["v2-2"] == pytest.approx(1e-2)


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


class TestMeshCompareDuplicateGeometry:
    """The free-surface output writes a face of an elastic-acoustic interface once
    for either side: the same geometry twice, with values that may differ, told apart
    by the locationFlag and the global-id. The two sides must never be swapped."""

    GEOM = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]]
    )
    CONNECT = np.array([[0, 1, 2], [0, 1, 2], [1, 3, 2]])
    VALUES = np.array([1.0, 5.0, 2.0])
    FLAGS = np.array([0, 1, 3])
    IDS = np.array([10, 11, 12])

    @pytest.mark.parametrize("order", [[0, 1, 2], [1, 0, 2]])
    @pytest.mark.parametrize("with_flags", [True, False])
    def test_sides_are_not_swapped(self, patch_seissolxdmf, capsys, order, with_flags):
        def entry(permutation):
            int_fields = {"global-id": self.IDS[permutation]}
            if with_flags:
                int_fields["locationFlag"] = self.FLAGS[permutation]
            return {
                "geom": self.GEOM,
                "connect": self.CONNECT[permutation],
                "fields": {"v1": self.VALUES[permutation]},
                "int_fields": int_fields,
            }

        patch_seissolxdmf["sim.xdmf"] = entry(order)
        patch_seissolxdmf["ref.xdmf"] = entry([0, 1, 2])

        meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=1e-12)
        out = capsys.readouterr().out
        assert "Matched all 3 cells geometrically" in out
        assert "conformant: True" in out


class TestMeshCompareUnreadableFlags:
    """A locationFlag that does not come back with one value per cell, as happens when a
    reader misinterprets its width, has to stop the comparison instead of shortening it.
    """

    def test_short_flags_fail(self, patch_seissolxdmf, capsys):
        geom = TestMeshCompareDuplicateGeometry.GEOM
        connect = TestMeshCompareDuplicateGeometry.CONNECT
        values = TestMeshCompareDuplicateGeometry.VALUES
        ids = TestMeshCompareDuplicateGeometry.IDS
        patch_seissolxdmf["sim.xdmf"] = {
            "geom": geom,
            "connect": connect,
            "fields": {"v1": values},
            "int_fields": {"global-id": ids, "locationFlag": np.array([0])},
        }
        patch_seissolxdmf["ref.xdmf"] = {
            "geom": geom,
            "connect": connect,
            "fields": {"v1": values},
            "int_fields": {"global-id": ids, "locationFlag": np.array([0, 1, 3])},
        }

        with pytest.raises(SystemExit) as exc_info:
            meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=1e-12)
        assert exc_info.value.code == 1
        assert "cannot be read as described" in capsys.readouterr().out


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


# A 120 km wide free surface as in tpv5. The x coordinate -13002.9759 lies exactly half-way
# between two points of the grid scale * 1e-9 = 1.2e-4 the matching used to snap to, so a
# single ulp of round-off decided which way it was rounded.
SURFACE_GEOM = np.array(
    [
        [-60000.0, -60000.0, 0.0],
        [60000.0, -60000.0, 0.0],
        [60000.0, 60000.0, 0.0],
        [-60000.0, 60000.0, 0.0],
        [-13002.9759, 5184.76159817, 0.0],
    ]
)
SURFACE_CONNECT = np.array([[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]])


class TestMeshCompareRoundOff:
    @pytest.mark.parametrize("ulps", [-2, -1, 1, 2])
    def test_vertex_on_a_grid_boundary_still_matches(
        self, patch_seissolxdmf, capsys, ulps
    ):
        shifted = SURFACE_GEOM.copy()
        for _ in range(abs(ulps)):
            shifted[4, 0] = np.nextafter(shifted[4, 0], np.sign(ulps) * np.inf)
        values = np.array([1.0, 2.0, 3.0, 4.0])
        patch_seissolxdmf["sim.xdmf"] = {
            "geom": shifted,
            "connect": SURFACE_CONNECT[:, ::-1],
            "fields": {"v1": values},
            "int_fields": {"global-id": np.array([4, 9, 14, 19])},
        }
        patch_seissolxdmf["ref.xdmf"] = {
            "geom": SURFACE_GEOM,
            "connect": SURFACE_CONNECT,
            "fields": {"v1": values},
            "int_fields": {"global-id": np.array([4, 9, 14, 19])},
        }
        meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=1e-12)
        assert "Matched all 4 cells geometrically" in capsys.readouterr().out

    def test_labels_ignore_round_off_but_separate_vertices(self):
        x = -13002.9759
        points = np.array(
            [[x, 0.0, 0.0], [np.nextafter(x, np.inf), 0.0, 0.0], [x + 1e-3, 0.0, 0.0]]
        )
        labels = meshcompare.point_labels(points, 120000.0 * 1e-9)
        assert labels[0] == labels[1] != labels[2]


# One surface triangle per element, split differently in the two files (by its edge
# midpoints, and by its barycenter), so that no cell matches and the comparison falls back
# to per-element means. The global-id is 4 * element + local side, and the local side
# differs between the two files, as it does between runs on differently labelled meshes.
def split4(a, b, c):
    ab, bc, ca = (a + b) / 2, (b + c) / 2, (c + a) / 2
    return [(a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)]


def split3(a, b, c):
    m = (a + b + c) / 3
    return [(a, b, m), (b, c, m), (c, a, m)]


def surface(triangles):
    geom = np.array([p for tri in triangles for p in tri])
    return geom, np.arange(len(geom)).reshape(-1, 3)


PARENTS = [
    (np.array([0.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])),
    (np.array([1.0, 0.0, 0.0]), np.array([1.0, 1.0, 0.0]), np.array([0.0, 1.0, 0.0])),
]


class TestMeshCompareFallbackWithoutComparableIds:
    def entries(self, means_ref, ids_ref):
        geom, connect = surface([t for p in PARENTS for t in split4(*p)])
        geom_ref, connect_ref = surface([t for p in PARENTS for t in split3(*p)])
        return (
            {
                "geom": geom,
                "connect": connect,
                "fields": {"v1": np.repeat([1.0, 5.0], 4)},
                # element 7 on side 1, element 3 on side 0
                "int_fields": {"global-id": np.repeat([29, 12], 4)},
            },
            {
                "geom": geom_ref,
                "connect": connect_ref,
                "fields": {"v1": np.repeat(means_ref, 3)},
                "int_fields": {"global-id": np.repeat(ids_ref, 3)},
            },
        )

    # the same elements on other sides, as on a relabelled mesh
    def test_elements_are_paired_by_position(self, patch_seissolxdmf, capsys):
        sim, ref = self.entries([1.0, 5.0], [31, 14])
        patch_seissolxdmf["sim.xdmf"], patch_seissolxdmf["ref.xdmf"] = sim, ref
        meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=1e-12)
        out = capsys.readouterr().out
        assert "Falling back to a per-element comparison" in out
        assert "Aggregated 8 cells into 2 elements" in out
        assert "conformant: True" in out

    # other elements altogether, for which even the sorted order flips: the elements are
    # still paired by position, but their ids have to agree as for matched cells
    def test_other_elements_fail(self, patch_seissolxdmf, capsys):
        sim, ref = self.entries([1.0, 5.0], [9, 38])
        patch_seissolxdmf["sim.xdmf"], patch_seissolxdmf["ref.xdmf"] = sim, ref
        with pytest.raises(SystemExit):
            meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=1e-12)
        out = capsys.readouterr().out
        assert "Aggregated 8 cells into 2 elements" in out
        assert "Global IDs present, but did not match" in out

    def test_elements_of_another_size_fail(self, patch_seissolxdmf, capsys):
        """Paired by the same centroid, but the reference element is larger."""
        sim, ref = self.entries([1.0, 5.0], [31, 14])
        centroid = np.mean(PARENTS[0], axis=0)
        grown = [tuple(centroid + 1.5 * (p - centroid) for p in PARENTS[0]), PARENTS[1]]
        ref["geom"], ref["connect"] = surface([t for p in grown for t in split3(*p)])
        patch_seissolxdmf["sim.xdmf"], patch_seissolxdmf["ref.xdmf"] = sim, ref
        with pytest.raises(SystemExit):
            meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=0.01)
        assert "do not cover the same volume" in capsys.readouterr().out

    @pytest.mark.parametrize("ids_ref", [[31, 14], [9, 38]])
    def test_wrong_means_still_fail(self, patch_seissolxdmf, ids_ref):
        sim, ref = self.entries([5.0, 1.0], ids_ref)
        patch_seissolxdmf["sim.xdmf"], patch_seissolxdmf["ref.xdmf"] = sim, ref
        with pytest.raises(SystemExit):
            meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=0.01)

    def test_elements_in_different_places_fail(self, patch_seissolxdmf, capsys):
        sim, ref = self.entries([1.0, 5.0], [31, 14])
        ref["geom"] = ref["geom"] + np.array([0.0, 0.0, 1.0])
        patch_seissolxdmf["sim.xdmf"], patch_seissolxdmf["ref.xdmf"] = sim, ref
        with pytest.raises(SystemExit):
            meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=0.01)
        assert "not in the same place" in capsys.readouterr().out

    @pytest.mark.parametrize("swap", [False, True])
    def test_interface_sides_are_not_swapped(self, patch_seissolxdmf, capsys, swap):
        """An elastic-acoustic interface face is written once per side, with the same
        centroid; the locationFlag keeps the two apart."""
        parent = PARENTS[0]
        geom, connect = surface(split4(*parent) * 2)
        geom_ref, connect_ref = surface(split3(*parent) * 2)
        order = [1, 0] if swap else [0, 1]
        patch_seissolxdmf["sim.xdmf"] = {
            "geom": geom,
            "connect": connect,
            "fields": {"v1": np.repeat([1.0, 5.0], 4)},
            "int_fields": {
                "global-id": np.repeat([29, 42], 4),
                "locationFlag": np.repeat([0, 1], 4),
            },
        }
        patch_seissolxdmf["ref.xdmf"] = {
            "geom": geom_ref,
            "connect": connect_ref,
            "fields": {"v1": np.repeat(np.array([1.0, 5.0])[order], 3)},
            "int_fields": {
                "global-id": np.repeat(np.array([31, 40])[order], 3),
                "locationFlag": np.repeat(np.array([0, 1])[order], 3),
            },
        }
        meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=1e-12)
        assert "Aggregated 8 cells into 2 elements" in capsys.readouterr().out


# The edge-midpoint split of a tetrahedron cuts its inner octahedron along one of three
# diagonals, and which one depends on the vertex order: relabelled meshes give another
# tiling, whose corner cells match and whose inner ones do not.
def refine8(a, b, c, d):
    ab, ac, ad, bc, bd, cd = (
        (a + b) / 2,
        (a + c) / 2,
        (a + d) / 2,
        (b + c) / 2,
        (b + d) / 2,
        (c + d) / 2,
    )
    return [
        (a, ab, ac, ad),
        (b, ab, bd, bc),
        (c, ac, bc, cd),
        (d, ad, cd, bd),
        (ab, ac, ad, bd),
        (ab, ac, bd, bc),
        (ac, ad, bd, cd),
        (ac, bc, cd, bd),
    ]


class TestMeshCompareRefine8:
    def test_other_diagonal_aggregates(self, patch_seissolxdmf, capsys):
        a, b, c, d = REFINED_GEOM[:4]
        geom = np.array([p for t in refine8(a, b, c, d) for p in t])
        geom_ref = np.array([p for t in refine8(b, c, a, d) for p in t])
        connect = np.arange(len(geom)).reshape(-1, 4)
        connect_ref = np.arange(len(geom_ref)).reshape(-1, 4)
        patch_seissolxdmf["sim.xdmf"] = {
            "geom": geom,
            "connect": connect,
            "fields": {"v1": np.full(8, 2.0)},
            "int_fields": {"global-id": np.full(8, 7)},
        }
        patch_seissolxdmf["ref.xdmf"] = {
            "geom": geom_ref,
            "connect": connect_ref,
            "fields": {"v1": np.full(8, 2.0)},
            "int_fields": {"global-id": np.full(8, 7)},
        }
        meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=1e-12)
        out = capsys.readouterr().out
        assert "Falling back to a per-element comparison" in out
        assert "conformant: True" in out

    def test_other_element_fails(self, patch_seissolxdmf, capsys):
        """A volume output carries the element itself, also in the fallback."""
        a, b, c, d = REFINED_GEOM[:4]
        geom = np.array([p for t in refine8(a, b, c, d) for p in t])
        geom_ref = np.array([p for t in refine8(b, c, a, d) for p in t])
        connect = np.arange(len(geom)).reshape(-1, 4)
        patch_seissolxdmf["sim.xdmf"] = {
            "geom": geom,
            "connect": connect,
            "fields": {"v1": np.full(8, 2.0)},
            "int_fields": {"global-id": np.full(8, 7)},
        }
        patch_seissolxdmf["ref.xdmf"] = {
            "geom": geom_ref,
            "connect": connect,
            "fields": {"v1": np.full(8, 2.0)},
            "int_fields": {"global-id": np.full(8, 6)},
        }
        with pytest.raises(SystemExit):
            meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=1e-12)
        assert "Global IDs present, but did not match" in capsys.readouterr().out


class TestMeshCompareGlobalIdGrouping:
    """Matched cells whose global-ids differ only by a renumbering of the groups."""

    def entry(self, ids):
        geom, connect = surface([t for p in PARENTS for t in split4(*p)])
        return {
            "geom": geom,
            "connect": connect,
            "fields": {"v1": np.arange(8.0)},
            "int_fields": {"global-id": np.asarray(ids)},
        }

    def test_renumbered_sides_are_conformant(self, patch_seissolxdmf, capsys):
        patch_seissolxdmf["sim.xdmf"] = self.entry(np.repeat([29, 12], 4))
        patch_seissolxdmf["ref.xdmf"] = self.entry(np.repeat([31, 14], 4))
        meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=1e-12)
        assert "conformant: True" in capsys.readouterr().out

    def test_regrouped_cells_are_not(self, patch_seissolxdmf, capsys):
        patch_seissolxdmf["sim.xdmf"] = self.entry(np.repeat([29, 12], 4))
        patch_seissolxdmf["ref.xdmf"] = self.entry([29, 29, 29, 12, 12, 12, 12, 12])
        with pytest.raises(SystemExit):
            meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=1e-12)
        assert "Global IDs present, but did not match" in capsys.readouterr().out

    def test_regrouped_faces_of_the_same_element_are_not(
        self, patch_seissolxdmf, capsys
    ):
        """The same element parts, but one cell is put into another face."""
        patch_seissolxdmf["sim.xdmf"] = self.entry(np.repeat([29, 12], 4))
        patch_seissolxdmf["ref.xdmf"] = self.entry([29, 29, 29, 28, 12, 12, 12, 12])
        with pytest.raises(SystemExit):
            meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=1e-12)
        assert "Global IDs present, but did not match" in capsys.readouterr().out

    def test_other_elements_are_not(self, patch_seissolxdmf, capsys):
        """The same grouping, but on other elements: the element part has to agree."""
        patch_seissolxdmf["sim.xdmf"] = self.entry(np.repeat([29, 12], 4))
        patch_seissolxdmf["ref.xdmf"] = self.entry(np.repeat([37, 16], 4))
        with pytest.raises(SystemExit):
            meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=1e-12)
        assert "Global IDs present, but did not match" in capsys.readouterr().out

    def test_renumbered_elements_of_a_volume_are_not(self, patch_seissolxdmf, capsys):
        """A volume output carries the element itself, which has to agree as it is."""
        entry = {
            "geom": REFINED_GEOM,
            "connect": REFINED_CONNECT,
            "fields": {"v1": REFINED_VALUES},
            "int_fields": {"global-id": REFINED_IDS},
        }
        patch_seissolxdmf["sim.xdmf"] = entry
        # 7 and 6 share the element part 7 // 4 == 6 // 4 a face output would compare, so
        # only the rule for volume outputs rejects them
        patch_seissolxdmf["ref.xdmf"] = dict(
            entry, int_fields={"global-id": REFINED_IDS - 1}
        )
        with pytest.raises(SystemExit):
            meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=1e-12)
        assert "Global IDs present, but did not match" in capsys.readouterr().out
