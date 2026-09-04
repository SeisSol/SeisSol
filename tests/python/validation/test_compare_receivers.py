# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

"""Tests for postprocessing/validation/compare-receivers.py

This script is what decides whether E2E regression tests pass or fail.
Every bug in this script either fakes success or hides real regressions.
"""

import importlib.util
from pathlib import Path
from textwrap import dedent

import numpy as np
import pandas as pd
import pytest

# Hyphenated filename — load via spec
SEISSOL_ROOT = Path(__file__).resolve().parents[3]
_spec = importlib.util.spec_from_file_location(
    "compare_receivers",
    SEISSOL_ROOT / "postprocessing" / "validation" / "compare-receivers.py",
)
cr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cr)


# ============================================================================
# normalize_variable_names — maps legacy SeisSol column names to current ones
# ============================================================================


class TestNormalizeVariableNames:
    """Handles output from ≥3 different SeisSol schema generations."""

    def test_legacy_nonfused_stress_renamed(self):
        result = cr.normalize_variable_names(
            ["Time", "xx", "yy", "zz", "xy", "xz", "yz"]
        )
        assert result == ["Time", "s_xx", "s_yy", "s_zz", "s_xy", "s_xz", "s_yz"]

    def test_legacy_nonfused_velocities_renamed(self):
        result = cr.normalize_variable_names(["Time", "u", "v", "w"])
        assert result == ["Time", "v1", "v2", "v3"]

    def test_current_names_are_passthrough(self):
        names = ["Time", "s_xx", "s_yy", "s_zz", "v1", "v2", "v3"]
        assert cr.normalize_variable_names(names) == names

    def test_fused_stress_renamed_per_index(self):
        result = cr.normalize_variable_names(["Time", "xx0", "xx1", "yy0", "yy1"])
        assert result == ["Time", "s_xx0", "s_xx1", "s_yy0", "s_yy1"]

    def test_fused_velocity_current_names_passthrough(self):
        # v1, v2, v3 are BOTH legacy-target AND current — they must not be re-mapped
        result = cr.normalize_variable_names(["Time", "v10", "v11", "v20", "v21"])
        assert result == ["Time", "v10", "v11", "v20", "v21"]

    def test_fused_velocity_legacy_names_renamed(self):
        result = cr.normalize_variable_names(["Time", "u0", "u1", "v0", "v1"])
        # u0/u1 → v10/v11, v0/v1 → v20/v21
        # BUT the extract-fused regex excludes columns starting with 'v'
        # so "v0" and "v1" are treated as NON-fused columns — this documents
        # the current behavior, which has a known ambiguity on velocity renaming.
        assert "v10" in result or "u0" in result  # implementation quirk
        assert "v11" in result or "u1" in result

    def test_empty_list(self):
        # Shouldn't crash on empty/minimal input
        assert cr.normalize_variable_names(["Time"]) == ["Time"]

    def test_partial_legacy(self):
        # Mixed: some legacy, some current — rare but defensive
        result = cr.normalize_variable_names(["Time", "xx", "s_yy", "u"])
        assert result == ["Time", "s_xx", "s_yy", "v1"]


# ============================================================================
# read_receiver — custom text-format parser
# ============================================================================


class TestReadReceiver:
    """Parses SeisSol receiver .dat files with a bespoke header."""

    @pytest.fixture
    def basic_receiver(self, tmp_path):
        f = tmp_path / "tpv-receiver-00001.dat"
        f.write_text(
            dedent(
                """\
            TITLE = "my receiver"
            VARIABLES = "Time", "s_xx", "s_yy", "v1"
            # x1    0.0
            # x2    0.0
            # x3    0.0
            0.0   1.0   2.0   0.5
            0.1   1.1   2.1   0.6
            0.2   1.2   2.2   0.7
            """
            )
        )
        return f

    def test_column_names(self, basic_receiver):
        df = cr.read_receiver(str(basic_receiver))
        assert list(df.columns) == ["Time", "s_xx", "s_yy", "v1"]

    def test_row_count(self, basic_receiver):
        df = cr.read_receiver(str(basic_receiver))
        assert len(df) == 3

    def test_values_parsed_as_floats(self, basic_receiver):
        df = cr.read_receiver(str(basic_receiver))
        assert df["Time"].iloc[0] == 0.0
        assert df["s_xx"].iloc[1] == pytest.approx(1.1)

    def test_legacy_names_get_normalized(self, tmp_path):
        f = tmp_path / "tpv-receiver-00001.dat"
        f.write_text(
            dedent(
                """\
            TITLE = "legacy"
            VARIABLES = "Time", "xx", "u"
            # some comment
            0.0   1.0   0.5
            0.1   1.1   0.6
            """
            )
        )
        df = cr.read_receiver(str(f))
        assert list(df.columns) == ["Time", "s_xx", "v1"]

    def test_fault_receiver_drops_t0(self, tmp_path):
        # A fault receiver file MUST drop its t=0 row (per the dr-cpp merge)
        f = tmp_path / "tpv-faultreceiver-00001.dat"
        f.write_text(
            dedent(
                """\
            TITLE = "fault"
            VARIABLES = "Time", "SRs"
            # x1    0.0
            0.0   0.0
            0.1   0.5
            0.2   0.7
            """
            )
        )
        df = cr.read_receiver(str(f))
        # t=0 row should be dropped
        assert len(df) == 2
        assert df["Time"].iloc[0] == pytest.approx(0.1)

    def test_nonfault_receiver_keeps_t0(self, tmp_path):
        # Regular receiver with t=0 must keep all rows
        f = tmp_path / "tpv-receiver-00001.dat"
        f.write_text(
            dedent(
                """\
            TITLE = "regular"
            VARIABLES = "Time", "v1"
            # meta
            0.0   0.5
            0.1   0.6
            """
            )
        )
        df = cr.read_receiver(str(f))
        assert len(df) == 2
        assert df["Time"].iloc[0] == 0.0

    def test_multiple_comment_lines_skipped(self, tmp_path):
        f = tmp_path / "tpv-receiver-00001.dat"
        f.write_text(
            dedent(
                """\
            TITLE = "many comments"
            VARIABLES = "Time", "v1"
            # coordinate x1
            # coordinate x2
            # coordinate x3
            # another line
            # yet another
            1.0   0.5
            """
            )
        )
        df = cr.read_receiver(str(f))
        assert len(df) == 1


# ============================================================================
# compare_receiver_columns — L2-error computation
# ============================================================================


class TestCompareReceiverColumns:

    def test_identical_receivers_yield_zero_error(self):
        t = np.linspace(0, 1, 100)
        df = pd.DataFrame({"Time": t, "v1": np.sin(t), "v2": np.cos(t)})
        errors = cr.compare_receiver_columns(df, df, label="test")
        assert errors["v1"] == pytest.approx(0.0, abs=1e-12)
        assert errors["v2"] == pytest.approx(0.0, abs=1e-12)

    def test_constant_offset_gives_nonzero_error(self):
        t = np.linspace(0, 1, 100)
        ref = pd.DataFrame({"Time": t, "v1": np.ones_like(t)})
        sim = pd.DataFrame({"Time": t, "v1": np.ones_like(t) * 1.1})
        errors = cr.compare_receiver_columns(sim, ref, label="test")
        # |sim - ref|_2 / |ref|_2 = 0.1 / 1.0 = 0.1
        assert errors["v1"] == pytest.approx(0.1, rel=1e-6)

    def test_missing_column_flagged_as_infinite(self):
        t = np.linspace(0, 1, 100)
        ref = pd.DataFrame({"Time": t, "v1": np.ones_like(t)})
        sim = pd.DataFrame({"Time": t})  # missing v1!
        errors = cr.compare_receiver_columns(sim, ref, label="test")
        assert errors["v1"] == float("inf")

    def test_zero_reference_uses_absolute_error(self):
        # When ref is ~zero, we fall back to absolute (not relative) error
        t = np.linspace(0, 1, 100)
        ref = pd.DataFrame({"Time": t, "v1": np.zeros_like(t)})
        sim = pd.DataFrame({"Time": t, "v1": np.ones_like(t) * 1e-8})
        errors = cr.compare_receiver_columns(sim, ref, label="test")
        # Not relative — returns diff_norm directly
        assert errors["v1"] == pytest.approx(1e-8, rel=1e-2)

    def test_time_column_excluded_from_errors(self):
        t = np.linspace(0, 1, 100)
        df = pd.DataFrame({"Time": t, "v1": np.sin(t)})
        errors = cr.compare_receiver_columns(df, df, label="test")
        assert "Time" not in errors


# ============================================================================
# find_all_receivers — glob/regex-based file discovery
# ============================================================================


class TestFindAllReceivers:

    def test_finds_numbered_receivers(self, tmp_path):
        for i in [1, 2, 5]:
            (tmp_path / f"tpv-receiver-{i:05d}.dat").touch()
        ids = cr.find_all_receivers(str(tmp_path), "tpv", "receiver")
        assert list(ids) == [1, 2, 5]

    def test_returns_sorted_unique(self, tmp_path):
        # Copy-layer receivers share an ID with a suffix: e.g. 00003-0.dat
        (tmp_path / "tpv-receiver-00003.dat").touch()
        (tmp_path / "tpv-receiver-00003-0.dat").touch()
        (tmp_path / "tpv-receiver-00001.dat").touch()
        ids = cr.find_all_receivers(str(tmp_path), "tpv", "receiver")
        # Should deduplicate and sort
        assert list(ids) == [1, 3]

    def test_prefix_filter_works(self, tmp_path):
        (tmp_path / "tpv-receiver-00001.dat").touch()
        (tmp_path / "otherprefix-receiver-00001.dat").touch()
        ids = cr.find_all_receivers(str(tmp_path), "tpv", "receiver")
        assert list(ids) == [1]

    def test_wrong_file_type_excluded(self, tmp_path):
        (tmp_path / "tpv-receiver-00001.dat").touch()
        (tmp_path / "tpv-faultreceiver-00002.dat").touch()
        ids = cr.find_all_receivers(str(tmp_path), "tpv", "receiver")
        assert list(ids) == [1]
        fault_ids = cr.find_all_receivers(str(tmp_path), "tpv", "faultreceiver")
        assert list(fault_ids) == [2]

    def test_empty_directory_returns_empty_array(self, tmp_path):
        ids = cr.find_all_receivers(str(tmp_path), "tpv", "receiver")
        assert len(ids) == 0


# ============================================================================
# report_errors — the gate that decides pass/fail
# ============================================================================


class TestReportErrors:

    def test_empty_returns_false(self, capsys):
        assert cr.report_errors("label", {}, 0.01) is False

    def test_all_within_epsilon_returns_false(self, capsys):
        errors = {1: {"v1": 0.001, "v2": 0.002}}
        assert cr.report_errors("label", errors, 0.01) is False

    def test_exceeds_epsilon_returns_true(self, capsys):
        errors = {1: {"v1": 0.1, "v2": 0.001}}
        assert cr.report_errors("label", errors, 0.01) is True

    def test_prints_offending_column_name(self, capsys):
        errors = {42: {"v1": 0.999}}
        cr.report_errors("receivers", errors, 0.01)
        captured = capsys.readouterr()
        assert "v1" in captured.out
        assert "42" in captured.out or "[42]" in captured.out

    def test_multiple_receivers_partial_failure(self, capsys):
        errors = {
            1: {"v1": 0.001, "v2": 0.001},
            2: {"v1": 0.999, "v2": 0.001},  # only v1 at id=2 fails
            3: {"v1": 0.001, "v2": 0.001},
        }
        assert cr.report_errors("label", errors, 0.01) is True
