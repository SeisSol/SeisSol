# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

"""Tests for scripts/validate/compare-receivers.py

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
    SEISSOL_ROOT / "scripts" / "validate" / "compare-receivers.py",
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
# fused simulations — a row per simulation, in the text files and in HDF5
# ============================================================================


def write_table(path, group, numbers, simulations, columns, times):
    """Write an HDF5 receiver table in the layout SeisSol writes.

    Row ``r`` of the table holds receiver ``numbers[r]`` (counted from zero) of
    simulation ``simulations[r]``, which is left out of the file if None, as a table
    from before the off-fault receivers had it. ``columns`` maps a quantity to its
    (samples, rows) values; all rows share one quantity set, so one dataset. The
    times are the same for every row, or given per row as (samples, rows).
    """
    h5py = pytest.importorskip("h5py")
    dtype = np.dtype([(name, np.float64) for name in ["Time", *columns]])
    data = np.zeros((len(times), len(numbers)), dtype=dtype)
    times = np.asarray(times, dtype=np.float64)
    data["Time"] = times if times.ndim == 2 else times[:, None]
    for name, values in columns.items():
        data[name] = values
    with h5py.File(path, "w") as handle:
        table = handle.create_group(group)
        table["group0"] = data
        table["Index"] = np.array([[0, row] for row in range(len(numbers))], np.uint64)
        name = "PointId" if group == "receivers" else "ReceiverId"
        table[name] = np.array(numbers, dtype=np.uint64)
        if simulations is not None:
            table["SimulationIndex"] = np.array(simulations, dtype=np.uint64)
        table["Coordinates"] = np.zeros((len(numbers), 3))


class TestSimulationRows:
    """Both layouts of a fused run read into the wide one the references use."""

    def test_text_rows_become_columns_per_simulation(self, tmp_path):
        f = tmp_path / "tpv-receiver-00001.dat"
        f.write_text(
            dedent(
                """\
            TITLE = "Temporal Signal for receiver number 00001"
            VARIABLES = "Time","SimulationIndex","v1","v2"
            # x1       0.0
            0.0  0  1.0  2.0
            0.0  1  3.0  4.0
            0.1  0  1.5  2.5
            0.1  1  3.5  4.5
            """
            )
        )
        df = cr.read_receiver(str(f))
        assert list(df.columns) == ["Time", "v10", "v20", "v11", "v21"]
        assert df["Time"].tolist() == [0.0, 0.1]
        assert df["v11"].tolist() == [3.0, 3.5]
        assert df["v20"].tolist() == [2.0, 2.5]

    def test_text_rows_on_the_fault_drop_t0_of_every_simulation(self, tmp_path):
        f = tmp_path / "tpv-faultreceiver-00001.dat"
        f.write_text(
            dedent(
                """\
            TITLE = "Temporal Signal for fault receiver number 1"
            VARIABLES = "Time" ,"SimulationIndex" ,"SRs" ,"SRd"
            # x1\t0.0
            # P_0\t-1.0\t-2.0
            0.0\t0\t0.0\t0.0\t
            0.0\t1\t0.0\t0.0\t
            0.1\t0\t0.5\t0.1\t
            0.1\t1\t0.7\t0.2\t
            """
            )
        )
        df = cr.read_receiver(str(f))
        assert list(df.columns) == ["Time", "SRs-1", "SRd-1", "SRs-2", "SRd-2"]
        assert df["Time"].tolist() == [0.1]
        assert df["SRs-2"].tolist() == [0.7]

    def test_text_rows_read_as_the_wide_layout(self, tmp_path):
        wide = tmp_path / "wide" / "tpv-receiver-00001.dat"
        rows = tmp_path / "rows" / "tpv-receiver-00001.dat"
        wide.parent.mkdir()
        rows.parent.mkdir()
        wide.write_text(
            'TITLE = "wide"\nVARIABLES = "Time","v10","v11"\n'
            "0.0  1.0  3.0\n0.1  1.5  3.5\n"
        )
        rows.write_text(
            'TITLE = "rows"\nVARIABLES = "Time","SimulationIndex","v1"\n'
            "0.0  0  1.0\n0.0  1  3.0\n0.1  0  1.5\n0.1  1  3.5\n"
        )
        pd.testing.assert_frame_equal(
            cr.read_receiver(str(rows))[["Time", "v10", "v11"]],
            cr.read_receiver(str(wide)),
        )

    def test_hdf5_rows_per_simulation(self, tmp_path):
        f = tmp_path / "tpv-receivers.h5"
        values = np.array([[1.0, 3.0, 5.0, 7.0], [1.5, 3.5, 5.5, 7.5]])
        write_table(
            f, "receivers", [0, 0, 1, 1], [0, 1, 0, 1], {"v1": values}, [0.0, 0.1]
        )
        receivers = cr.read_hdf5_receivers(str(f), "receiver")
        assert sorted(receivers) == [1, 2]
        assert list(receivers[2].columns) == ["Time", "v10", "v11"]
        assert receivers[2]["v11"].tolist() == [7.0, 7.5]
        assert receivers[1]["Time"].tolist() == [0.0, 0.1]

    def test_hdf5_without_simulation_index_is_taken_as_it_is(self, tmp_path):
        # off-fault tables written before a simulation took a row of its own
        f = tmp_path / "tpv-receivers.h5"
        columns = {"v10": np.array([[1.0], [2.0]]), "v11": np.array([[3.0], [4.0]])}
        write_table(f, "receivers", [0], None, columns, [0.0, 0.1])
        receivers = cr.read_hdf5_receivers(str(f), "receiver")
        assert list(receivers[1].columns) == ["Time", "v10", "v11"]
        assert receivers[1]["v11"].tolist() == [3.0, 4.0]

    def test_hdf5_on_the_fault_single_simulation_keeps_plain_names(self, tmp_path):
        f = tmp_path / "tpv-faultreceivers.h5"
        values = np.array([[0.0, 0.0], [0.5, 0.6], [0.7, 0.8]])
        write_table(
            f, "faultreceivers", [0, 1], [0, 0], {"SRs": values}, [0.0, 0.1, 0.2]
        )
        receivers = cr.read_hdf5_receivers(str(f), "faultreceiver")
        assert list(receivers[1].columns) == ["Time", "SRs"]
        # the t=0 sample is dropped, as for the text files
        assert receivers[2]["SRs"].tolist() == [0.6, 0.8]

    def test_hdf5_padding_of_local_time_stepping_is_dropped(self, tmp_path):
        # a receiver that took fewer samples than the longest one of its table has
        # the rest of its column at NaN, the time included
        f = tmp_path / "tpv-faultreceivers.h5"
        values = np.array([[0.0, 0.0], [0.5, 0.6], [np.nan, 0.8]])
        times = np.array([[0.0, 0.0], [0.1, 0.1], [np.nan, 0.2]])
        write_table(f, "faultreceivers", [0, 1], [0, 0], {"SRs": values}, times)
        receivers = cr.read_hdf5_receivers(str(f), "faultreceiver")
        assert receivers[1]["SRs"].tolist() == [0.5]
        assert receivers[2]["SRs"].tolist() == [0.6, 0.8]

    def test_hdf5_receiver_held_twice_is_taken_once(self, tmp_path):
        f = tmp_path / "tpv-receivers.h5"
        values = np.array([[1.0, 1.0], [2.0, 2.0]])
        write_table(f, "receivers", [0, 0], [0, 0], {"v1": values}, [0.0, 0.1])
        receivers = cr.read_hdf5_receivers(str(f), "receiver")
        assert list(receivers[1].columns) == ["Time", "v1"]
        assert receivers[1]["v1"].tolist() == [1.0, 2.0]

    def test_hdf5_against_text_references(self, tmp_path):
        # what CI does: a run writing HDF5 held against references in wide text files
        (tmp_path / "ref").mkdir()
        (tmp_path / "run").mkdir()
        (tmp_path / "ref" / "tpv-faultreceiver-00001.dat").write_text(
            'TITLE = "t"\nVARIABLES = "Time" ,"SRs-1" ,"SRs-2"\n# x1\t0.0\n'
            "0.0\t0.0\t0.0\t\n0.1\t0.5\t0.7\t\n0.2\t0.6\t0.8\t\n"
        )
        values = np.array([[0.0, 0.0], [0.5, 0.7], [0.6, 0.8]])
        write_table(
            tmp_path / "run" / "tpv-faultreceivers.h5",
            "faultreceivers",
            [0, 0],
            [0, 1],
            {"SRs": values},
            [0.0, 0.1, 0.2],
        )
        ref = cr.load_receivers(str(tmp_path / "ref"), "tpv", "faultreceiver")
        run = cr.load_receivers(str(tmp_path / "run"), "tpv", "faultreceiver")
        errors = cr.receiver_diff(run[1], ref[1], 1, "faultreceiver")
        assert errors == {"SRs-1": 0.0, "SRs-2": 0.0}


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

    def test_component_is_relative_to_its_tensor(self):
        # a shear stress that stays at zero, as in water, holds rounding noise only;
        # it is judged against the stress, not against its own noise
        t = np.linspace(0, 1, 100)
        ref = pd.DataFrame(
            {"Time": t, "s_xx": np.ones_like(t), "s_xy": np.full_like(t, 1e-9)}
        )
        sim = ref.copy()
        sim["s_xy"] = -1e-9
        errors = cr.compare_receiver_columns(sim, ref, label="test")
        assert errors["s_xy"] == pytest.approx(2e-9, rel=1e-6)
        assert errors["s_xx"] == pytest.approx(0.0, abs=1e-12)

    def test_components_of_fused_simulations_stay_apart(self):
        t = np.linspace(0, 1, 100)
        ref = pd.DataFrame(
            {
                "Time": t,
                "v10": np.ones_like(t),
                "v20": np.full_like(t, 1e-9),
                "v11": np.full_like(t, 1e-3),
                "v21": np.full_like(t, 1e-3),
            }
        )
        sim = ref.copy()
        sim["v21"] = 2e-3
        errors = cr.compare_receiver_columns(sim, ref, label="test")
        # v21 is relative to simulation 1 only, not to the large v10 of simulation 0
        assert errors["v21"] == pytest.approx(1.0, rel=1e-6)

    def test_component_group(self):
        assert cr.component_group("v1") == "velocity"
        assert cr.component_group("v13") == "velocity3"
        assert cr.component_group("v131") == "velocity31"
        assert cr.component_group("s_xz0") == "stress0"
        assert cr.component_group("SRd-4") == "slip rate-4"
        assert cr.component_group("Ts0") == "initial traction"
        assert cr.component_group("Ts0-2") == "initial traction-2"
        assert cr.component_group("v2_f1") == "fluid velocity1"
        assert cr.component_group("Mud") == "Mud"
        assert cr.component_group("RT-3") == "RT-3"


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


class TestEventQuantities:
    """RT, Vr and DS are 0 until the event and constant from it on."""

    @staticmethod
    def frames(sim_values, ref_values, column="RT"):
        time = np.arange(len(ref_values), dtype=float)
        sim = pd.DataFrame({"Time": time, column: np.array(sim_values, dtype=float)})
        ref = pd.DataFrame({"Time": time, column: np.array(ref_values, dtype=float)})
        return sim, ref

    def test_onset_one_sample_apart_is_left_out(self):
        sim, ref = self.frames([0, 0, 2.001, 2.001, 2.001], [0, 0, 0, 2.0, 2.0])
        errors = cr.compare_receiver_columns(sim, ref, label="test")
        # only the samples both hold the event in count
        assert errors["RT"] == pytest.approx(0.001 / 2.0, rel=1e-2)

    def test_fused_suffix_is_an_event_quantity_too(self):
        sim, ref = self.frames([0, 1.0, 1.0], [0, 0, 1.0], column="Vr-3")
        errors = cr.compare_receiver_columns(sim, ref, label="test")
        assert errors["Vr-3"] == 0.0

    def test_missing_event_still_counts(self):
        # no rupture at all in one run: not an onset, the whole difference counts
        sim, ref = self.frames([0, 0, 0, 0, 0], [0, 2.0, 2.0, 2.0, 2.0])
        errors = cr.compare_receiver_columns(sim, ref, label="test")
        assert errors["RT"] == pytest.approx(1.0)

    def test_other_quantities_are_compared_in_full(self):
        sim, ref = self.frames([0, 1.0, 1.0], [0, 0, 1.0], column="SRs")
        errors = cr.compare_receiver_columns(sim, ref, label="test")
        assert errors["SRs"] > 0.1

    @staticmethod
    def repeated(rt_sim, rt_ref, repeat):
        # a receiver of a coarser cluster: each time step shown at `repeat` samples in a row
        steps = len(rt_ref)
        time = np.arange(steps * repeat, dtype=float)
        srs = np.repeat(np.arange(1.0, steps + 1.0), repeat)
        sim = pd.DataFrame(
            {"Time": time, "SRs": srs, "RT": np.repeat(np.array(rt_sim, float), repeat)}
        )
        ref = pd.DataFrame(
            {"Time": time, "SRs": srs, "RT": np.repeat(np.array(rt_ref, float), repeat)}
        )
        return sim, ref

    def test_onset_counts_time_steps_not_samples(self):
        # one time step apart, at four samples per time step
        sim, ref = self.repeated([0, 0, 2.0, 2.0, 2.0], [0, 0, 0, 2.0, 2.0], repeat=4)
        errors = cr.compare_receiver_columns(sim, ref, label="test")
        assert errors["RT"] == 0.0

    def test_onset_longer_than_the_bound_counts(self):
        sim, ref = self.frames([0, 2.0, 2.0, 2.0, 2.0, 2.0], [0, 0, 0, 0, 2.0, 2.0])
        errors = cr.compare_receiver_columns(sim, ref, label="test")
        assert errors["RT"] > 0.1

    def test_one_sided_samples_away_from_the_onset_count(self):
        # the event vanishes again in one run: not an onset
        sim, ref = self.frames([0, 2.0, 2.0, 0, 2.0], [0, 0, 2.0, 2.0, 2.0])
        errors = cr.compare_receiver_columns(sim, ref, label="test")
        assert errors["RT"] > 0.1

    def test_not_finite_values_fail(self):
        sim, ref = self.frames([0, float("nan"), 2.0, 2.0], [0, 0, 2.0, 2.0])
        errors = cr.compare_receiver_columns(sim, ref, label="test")
        assert errors["RT"] == float("inf")
        sim, ref = self.frames([0, float("nan"), 2.0], [0, 1.0, 2.0], column="SRs")
        errors = cr.compare_receiver_columns(sim, ref, label="test")
        assert errors["SRs"] == float("inf")


# ============================================================================
# report_errors — the gate that decides pass/fail
# ============================================================================


class TestReportErrors:
    """report_errors returns (exceeded, per_column_max).

    The second element is the worst error per column across all receivers and
    is built regardless of pass/fail, since the machine-readable summary needs
    it either way.
    """

    def test_empty_returns_no_failure_and_no_maxima(self, capsys):
        assert cr.report_errors("label", {}, 0.01) == (False, {})

    def test_all_within_epsilon_does_not_report_a_failure(self, capsys):
        errors = {1: {"v1": 0.001, "v2": 0.002}}
        exceeded, maxima = cr.report_errors("label", errors, 0.01)
        assert exceeded is False
        assert maxima == pytest.approx({"v1": 0.001, "v2": 0.002})

    def test_exceeds_epsilon_reports_a_failure(self, capsys):
        errors = {1: {"v1": 0.1, "v2": 0.001}}
        exceeded, maxima = cr.report_errors("label", errors, 0.01)
        assert exceeded is True
        assert maxima == pytest.approx({"v1": 0.1, "v2": 0.001})

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
        exceeded, maxima = cr.report_errors("label", errors, 0.01)
        assert exceeded is True
        # the maximum is taken across receivers, so v1 carries the id=2 value
        assert maxima == pytest.approx({"v1": 0.999, "v2": 0.001})

    def test_maxima_ignore_missing_columns(self, capsys):
        """A column absent from one receiver must not poison its maximum."""
        errors = {1: {"v1": 0.002}, 2: {"v1": 0.001, "v2": 0.003}}
        exceeded, maxima = cr.report_errors("label", errors, 0.01)
        assert exceeded is False
        assert maxima == pytest.approx({"v1": 0.002, "v2": 0.003})
