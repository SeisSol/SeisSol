# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

"""Tests for postprocessing/validation/meshcompare.py

meshcompare.compare() is monolithic — it:
  1. Opens two XDMF files via seissolxdmf
  2. Reorders cells (by global-id field if present, else by barycenter heuristic)
  3. Asserts geometries match
  4. Computes L1/L2 errors per quantity
  5. Has a hardcoded workaround for a known SeisSol bug: DS field zeros
  6. Calls sys.exit(1) on threshold violation

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

    def test_identical_without_global_id_falls_back_to_barycenter(
        self, patch_seissolxdmf, capsys
    ):
        data = {
            "geom": GEOM,
            "connect": CONNECT,
            "fields": {"v1": np.array([1.0, 2.0])},
        }
        patch_seissolxdmf["sim.xdmf"] = data
        patch_seissolxdmf["ref.xdmf"] = data

        meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=0.01)
        out = capsys.readouterr().out
        assert "Order the cells by their barycenter" in out


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

    def test_geometry_mismatch_raises_assertion(self, patch_seissolxdmf):
        """meshcompare has an inline `assert np.all(... < 1e-10)` on geometry"""
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
        with pytest.raises(AssertionError):
            meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=0.01)

    def test_mismatched_global_ids_fail(self, patch_seissolxdmf):
        """If global-id field is present but cells don't align geometrically,
        the inline assertion on geometry equality fires (AssertionError)."""
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
        # The inline assert `np.all(np.abs(geom[connect]-geom_ref[connect_ref])<1e-10)`
        # fires on mismatched cell geometry
        with pytest.raises(AssertionError):
            meshcompare.compare("sim.xdmf", "ref.xdmf", epsilon=0.01)


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
