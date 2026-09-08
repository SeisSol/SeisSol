# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

"""Tests for scripts/validate/compare-energies.py

This script gates energy-conservation and seismic-moment regression
checks in the E2E CI. Handles two distinct CSV schemas (pre/post PR #773).
"""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SEISSOL_ROOT = Path(__file__).resolve().parents[3]
_spec = importlib.util.spec_from_file_location(
    "compare_energies",
    SEISSOL_ROOT / "scripts" / "validate" / "compare-energies.py",
)
ce = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ce)


@pytest.fixture(autouse=True)
def _inject_np_pd_into_module_namespace():
    """Workaround for Bug #1 below: perform_check uses np.any() but np is
    only imported in `if __name__ == "__main__"`. Without this fixture,
    every perform_check test would fail with NameError.
    """
    ce.np = np
    ce.pd = pd
    yield


# =============================================================================
# pivot_if_necessary — normalizes the two historical CSV formats
# =============================================================================


class TestPivotIfNecessary:

    def test_old_format_passthrough(self):
        df = pd.DataFrame(
            {
                "time": [0.0, 0.1, 0.2],
                "elastic_energy": [1.0, 2.0, 3.0],
                "elastic_kinetic_energy": [0.5, 0.6, 0.7],
            }
        )
        result = ce.pivot_if_necessary(df)
        pd.testing.assert_frame_equal(result, df)

    def test_new_long_format_is_pivoted(self):
        df = pd.DataFrame(
            {
                "time": [0, 0, 1, 1],
                "variable": ["elastic_energy", "elastic_kinetic_energy"] * 2,
                "measurement": [1.0, 0.5, 2.0, 0.6],
            }
        )
        result = ce.pivot_if_necessary(df)
        assert "elastic_energy" in result.columns
        assert "elastic_kinetic_energy" in result.columns
        assert len(result) == 2

    def test_new_format_fused_keeps_simulation_index(self):
        df = pd.DataFrame(
            {
                "time": [0, 0, 0, 0, 1, 1, 1, 1],
                "simulation_index": [0, 0, 1, 1, 0, 0, 1, 1],
                "variable": ["elastic_energy", "elastic_kinetic_energy"] * 4,
                "measurement": [1.0, 0.5, 1.1, 0.55, 2.0, 0.6, 2.1, 0.65],
            }
        )
        result = ce.pivot_if_necessary(df)
        assert "simulation_index" in result.columns
        assert len(result) == 4

    def test_pivot_preserves_int_types_on_index(self):
        df = pd.DataFrame(
            {
                "time": [0, 1, 2],
                "variable": ["A", "A", "A"],
                "measurement": [1.0, 2.0, 3.0],
            }
        )
        result = ce.pivot_if_necessary(df)
        assert result["time"].dtype == np.int64


# =============================================================================
# get_number_of_fused_sims — documents current (partially buggy) behavior
# =============================================================================


class TestGetNumberOfFusedSims:

    def test_fused_by_simulation_index(self):
        df = pd.DataFrame(
            {
                "time": [0, 0, 1, 1],
                "simulation_index": [0, 1, 0, 1],
                "elastic_energy": [1.0, 1.1, 2.0, 2.1],
            }
        )
        assert ce.get_number_of_fused_sims(df) == 2

    def test_higher_simulation_indices(self):
        df = pd.DataFrame(
            {
                "time": [0, 0, 0, 0],
                "simulation_index": [0, 1, 2, 3],
                "elastic_energy": [1.0, 2.0, 3.0, 4.0],
            }
        )
        assert ce.get_number_of_fused_sims(df) == 4

    def test_pure_digit_suffix_columns_detected(self):
        df = pd.DataFrame(
            {
                "elastic_energy0": [1.0],
                "elastic_energy1": [1.1],
                "elastic_energy2": [1.2],
            }
        )
        assert ce.get_number_of_fused_sims(df) == 3

    def test_nonfused_new_format_returns_neg1(self):
        # Non-digit last chars cause the int() call to raise -> -1
        df = pd.DataFrame({"time": [0, 1], "elastic_energy": [1.0, 2.0]})
        assert ce.get_number_of_fused_sims(df) == -1


# =============================================================================
# get_sub_simulation — slice out one of N fused sims
# =============================================================================


class TestGetSubSimulation:

    def test_nonfused_returns_full_dataframe(self):
        df = pd.DataFrame({"time": [0, 1, 2], "elastic_energy": [1.0, 2.0, 3.0]})
        result = ce.get_sub_simulation(df, fused_index=0)
        assert len(result) == 3
        assert "elastic_energy" in result.columns

    def test_fused_selects_only_that_index(self):
        df = pd.DataFrame(
            {
                "time": [0, 0, 1, 1],
                "simulation_index": [0, 1, 0, 1],
                "elastic_energy": [1.0, 1.1, 2.0, 2.1],
            }
        )
        sim0 = ce.get_sub_simulation(df, fused_index=0)
        sim1 = ce.get_sub_simulation(df, fused_index=1)
        assert list(sim0["elastic_energy"]) == [1.0, 2.0]
        assert list(sim1["elastic_energy"]) == [1.1, 2.1]

    def test_fused_nonexistent_index_returns_empty(self):
        df = pd.DataFrame(
            {
                "time": [0, 0],
                "simulation_index": [0, 1],
                "elastic_energy": [1.0, 1.1],
            }
        )
        result = ce.get_sub_simulation(df, fused_index=99)
        assert len(result) == 0


# =============================================================================
# perform_check — the green/red gate (uses np injection from autouse fixture)
# =============================================================================


class TestPerformCheck:
    """perform_check returns (exceeded, relative_difference).

    Unpack it. A bare `assert perform_check(...)` passes for any two-element
    tuple and therefore asserts nothing.
    """

    @staticmethod
    def exceeded(sim, ref, epsilon):
        result, _ = ce.perform_check(sim, ref, epsilon=epsilon)
        return result

    def test_identical_frames_pass(self):
        df = pd.DataFrame(
            {
                "elastic_energy": [1.0, 1.1, 1.2],
                "total_frictional_work": [0.0, 0.5, 1.0],
            }
        )
        assert self.exceeded(df, df, 0.01) is False

    def test_large_difference_fails(self):
        ref = pd.DataFrame({"elastic_energy": [1.0, 1.1, 1.2]})
        sim = pd.DataFrame({"elastic_energy": [1.0, 1.5, 1.8]})
        assert self.exceeded(sim, ref, 0.05) is True

    def test_exact_threshold_comparison(self):
        ref = pd.DataFrame({"elastic_energy": [1.0, 1.0, 1.0]})
        sim_high = pd.DataFrame({"elastic_energy": [1.0, 1.051, 1.051]})
        sim_low = pd.DataFrame({"elastic_energy": [1.0, 1.049, 1.049]})
        assert self.exceeded(sim_high, ref, 0.05) is True
        assert self.exceeded(sim_low, ref, 0.05) is False

    def test_first_row_excluded(self):
        # A huge difference in row 0 should not fail (iloc[1:] drops it)
        ref = pd.DataFrame({"elastic_energy": [1.0, 1.0, 1.0]})
        sim = pd.DataFrame({"elastic_energy": [999.0, 1.0, 1.0]})
        assert self.exceeded(sim, ref, 0.01) is False

    def test_multiple_quantities_any_failure_is_failure(self):
        ref = pd.DataFrame(
            {
                "elastic_energy": [1.0, 1.0, 1.0],
                "total_frictional_work": [1.0, 1.0, 1.0],
            }
        )
        sim = pd.DataFrame(
            {
                "elastic_energy": [1.0, 1.0, 1.0],
                "total_frictional_work": [1.0, 2.0, 2.0],
            }
        )
        assert self.exceeded(sim, ref, 0.01) is True

    def test_relative_difference_is_handed_back(self):
        """The caller builds the machine-readable summary from this."""
        ref = pd.DataFrame({"elastic_energy": [1.0, 1.0, 1.0]})
        sim = pd.DataFrame({"elastic_energy": [1.0, 1.5, 1.25]})
        _, rel_diff = ce.perform_check(sim, ref, epsilon=0.05)
        # t=0 is dropped, so two rows remain
        assert list(rel_diff.index) == [1, 2]
        assert rel_diff["elastic_energy"].tolist() == pytest.approx([0.5, 0.25])

    def test_relative_difference_is_returned_even_when_passing(self):
        ref = pd.DataFrame({"elastic_energy": [1.0, 1.0, 1.0]})
        exceeded, rel_diff = ce.perform_check(ref, ref, epsilon=0.01)
        assert exceeded is False
        assert rel_diff["elastic_energy"].tolist() == pytest.approx([0.0, 0.0])

    def test_prints_content_for_debugging(self, capsys):
        df = pd.DataFrame({"elastic_energy": [1.0, 1.0, 1.0]})
        ce.perform_check(df, df, epsilon=0.01)
        out = capsys.readouterr().out
        assert "Energies" in out
        assert "Relative difference" in out


class TestDocumentedFormerBugs:

    def test_bug2_get_number_of_fused_sims_fails_on_mixed_columns(self):
        """BUG: The old fused CSV format had columns like
            ['time', 'elastic_energy0', 'elastic_energy1', ...]
        The function iterates ALL columns including 'time' and does
        int(c[-1]). int('e') raises, the except swallows it, returns -1 —
        i.e. it silently reports "not fused" even when it IS fused.

        Consequence: old-format fused-sim runs in the precomputed repo
        might be slipping through as non-fused, masking per-sub-simulation
        regressions in CI.

        Fix: filter columns by re.match(r'.*(\\d+)$', c) before int-casting.
        """
        old_fused_format = pd.DataFrame(
            {
                "time": [0.0, 0.1, 0.2],
                "elastic_energy0": [1.0, 1.1, 1.2],
                "elastic_energy1": [2.0, 2.1, 2.2],
            }
        )
        result = ce.get_number_of_fused_sims(old_fused_format)

        assert result == 2
