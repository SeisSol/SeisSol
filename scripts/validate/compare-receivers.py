#!/usr/bin/env python3

# SPDX-FileCopyrightText: 2022 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from validation_report import write_report_json

if hasattr(np, "trapezoid"):
    trapz_func = np.trapezoid
else:
    trapz_func = np.trapz


# Maps legacy variable names (written by older SeisSol versions) to current names.
_LEGACY_NAMES = {
    "xx": "s_xx",
    "yy": "s_yy",
    "zz": "s_zz",
    "xy": "s_xy",
    "xz": "s_xz",
    "yz": "s_yz",
    "u": "v1",
    "v": "v2",
    "w": "v3",
    "-p": "pprime",
}


def normalize_variable_names(variables: list[str]) -> list[str]:
    """Rename legacy column names to current ones, for both fused and non-fused files.

    Fused receiver files append a numeric simulation index to each column (e.g. "xx0",
    "u1"). Non-fused files use bare names (e.g. "xx", "u") or already-current names
    (e.g. "s_xx", "v1"). Velocity columns are excluded from fused-index detection to
    avoid ambiguity with the already-current names v1, v2, v3.
    """
    # Detect fused files via digit suffixes on non-velocity columns.
    extract_fused = re.compile(r"^[^v].*?(\d+)$")
    max_index = max(
        (int(m.group(1)) for col in variables[1:] if (m := extract_fused.search(col))),
        default=-2,
    )
    n_fused = max_index + 1  # negative means non-fused

    result = list(variables)
    for old, new in _LEGACY_NAMES.items():
        if n_fused < 1:
            # Non-fused: rename bare legacy names only (exact match).
            if old in result:
                result[result.index(old)] = new
        else:
            # Fused: rename "old{i}" -> "new{i}" for each simulation index.
            for i in range(n_fused):
                old_i, new_i = f"{old}{i}", f"{new}{i}"
                if old_i in result:
                    result[result.index(old_i)] = new_i
    return result


def simulation_suffix(simulation: int, file_type: str) -> str:
    """The suffix a wide text file gives the quantities of a fused simulation.

    Before a simulation took rows of its own, the text files of a fused run named
    the quantities after it: counted from zero and appended in the volume receivers
    (v10, v11, ...), counted from one and after a dash on the fault (SRs-1, ...).
    """
    return f"-{simulation + 1}" if file_type == "faultreceiver" else str(simulation)


def join_simulations(
    parts: list[tuple[int, pd.DataFrame]], file_type: str
) -> pd.DataFrame:
    """Put the simulations of one receiver side by side, a row per sample.

    Each part holds the samples of one simulation, with the time and the quantities
    under their plain names. A single simulation is the receiver as it is; fused ones
    get their quantities named as a wide text file names them (see simulation_suffix),
    which is the layout the references were recorded in.
    """
    parts = sorted(parts, key=lambda part: part[0])
    if len(parts) == 1 and parts[0][0] == 0:
        return parts[0][1].reset_index(drop=True)
    time = parts[0][1]["Time"].to_numpy()
    columns = [pd.DataFrame({"Time": time})]
    for simulation, frame in parts:
        assert np.array_equal(frame["Time"].to_numpy(), time), (
            f"the samples of simulation {simulation} are not taken at the times of "
            f"simulation {parts[0][0]}"
        )
        suffix = simulation_suffix(simulation, file_type)
        columns.append(
            frame.drop(columns="Time")
            .rename(columns=lambda name: name + suffix)
            .reset_index(drop=True)
        )
    return pd.concat(columns, axis=1)


def read_receiver(filename: str) -> pd.DataFrame:
    """
    Read the receiver using the receiver filename and return a pandas DataFrame
    """
    with open(filename) as receiver_file:
        # We expect the header to look like the following:
        # TITLE = "title"
        # VARIABLES = "variable 1", "variable 2", ..., "variable n"
        # # x1	x coordinate of receiver
        # # x2	y coordinate of receiver
        # # x3	z coordinate of receiver
        # # possibly more comments, starting with #
        lines = receiver_file.readlines()
        # remove the first 12 characters from the line ("VARIABLES = ")
        variable_line = lines[1][12:].split(",")
        variables = [s.strip().replace('"', "") for s in variable_line]
        # find first row without comments
        first_row = 2
        while first_row < len(lines) and lines[first_row][0] == "#":
            first_row += 1

        assert first_row < len(lines), f"Empty file: {filename}"
    receiver = pd.read_csv(filename, header=None, skiprows=first_row, sep=r"\s+")
    receiver.columns = variables
    name = os.path.basename(filename)
    file_type = "faultreceiver" if "faultreceiver" in name else "receiver"
    if "SimulationIndex" in variables:
        # a fused run writes a row per simulation, which that column names
        simulations = receiver["SimulationIndex"].to_numpy().astype(np.int64)
        receiver = join_simulations(
            [
                (
                    int(simulation),
                    receiver[simulations == simulation].drop(columns="SimulationIndex"),
                )
                for simulation in np.unique(simulations)
            ],
            file_type,
        )
    # since dr-cpp merge, fault receiver files start writing at Time=0
    # (before they were writing at Time=dt)
    # We then skip the first timestep written if Time = 0
    if (
        file_type == "faultreceiver"
        and len(receiver) > 0
        and receiver["Time"].iloc[0] == 0
    ):
        receiver = receiver.iloc[1:].reset_index(drop=True)
    receiver.columns = normalize_variable_names(list(receiver.columns))
    return receiver


# the HDF5 receiver files, one per kind, and the group everything sits under
_HDF5_NAMES = {"receiver": "receivers", "faultreceiver": "faultreceivers"}


def hdf5_receiver_file(directory: str, prefix: str, file_type: str) -> str | None:
    """The HDF5 file holding all receivers of one kind, if the run wrote one."""
    path = os.path.join(directory, f"{prefix}-{_HDF5_NAMES[file_type]}.h5")
    return path if os.path.isfile(path) else None


def read_hdf5_receivers(filename: str, file_type: str) -> dict[int, pd.DataFrame]:
    """Read every receiver of an HDF5 receiver file, by the number of its text file.

    A row of the table is one receiver of one simulation, and the receivers come in
    the layout of the wide text files, so that either format can be compared against
    the other (see join_simulations). A table of off-fault receivers written before
    they took a row per simulation has no SimulationIndex; its rows carry all
    simulations already, under the names of the text files. The samples a receiver
    did not take, which a table under local time stepping pads with NaN, are dropped,
    and so is the one at t = 0 of an on-fault receiver, as read_receiver does.
    """
    import h5py

    with h5py.File(filename, "r") as handle:
        group = handle[_HDF5_NAMES[file_type]]
        index = group["Index"][:]
        numbers = group["PointId" if file_type == "receiver" else "ReceiverId"][:] + 1
        if "SimulationIndex" in group:
            simulations = group["SimulationIndex"][:]
        else:
            simulations = np.zeros(len(numbers), dtype=np.int64)

        tables = {}
        parts: dict[int, dict[int, pd.DataFrame]] = {}
        for row, (table, column) in enumerate(index):
            if table not in tables:
                tables[table] = group[f"group{table}"][:]
            samples = tables[table][:, column]
            frame = pd.DataFrame(
                {name: samples[name].astype(np.float64) for name in samples.dtype.names}
            )
            frame = frame[np.isfinite(frame["Time"].to_numpy())]
            # a receiver more than one rank holds is taken once, as of its text files
            parts.setdefault(int(numbers[row]), {}).setdefault(
                int(simulations[row]), frame
            )

    receivers = {}
    for number, receiver_parts in parts.items():
        receiver = join_simulations(list(receiver_parts.items()), file_type)
        if (
            file_type == "faultreceiver"
            and len(receiver) > 0
            and receiver["Time"].iloc[0] == 0
        ):
            receiver = receiver.iloc[1:].reset_index(drop=True)
        receivers[number] = receiver
    return receivers


def load_receivers(
    directory: str, prefix: str, file_type: str = "receiver"
) -> dict[int, pd.DataFrame]:
    """Every receiver of one kind in a directory, by number, from either format.

    An HDF5 receiver file is read if the run wrote one, and the text files otherwise;
    of the text files of one receiver, which a receiver in a copy layer may have
    several of, the first is taken.
    """
    hdf5 = hdf5_receiver_file(directory, prefix, file_type)
    if hdf5 is not None:
        return read_hdf5_receivers(hdf5, file_type)
    receivers = {}
    for number in find_all_receivers(directory, prefix, file_type):
        files = sorted(glob.glob(f"{directory}/{prefix}-{file_type}-{number:05d}*.dat"))
        receiver = read_receiver(files[0])
        # the t=0 row of a fault receiver may have been dropped; count rows from zero
        receivers[int(number)] = receiver.reset_index(drop=True)
    return receivers


# receiver columns that are components of one vector or tensor
_COMPONENTS = {
    **{name: "stress" for name in ("s_xx", "s_yy", "s_zz", "s_xy", "s_yz", "s_xz")},
    **{
        name: "strain rate"
        for name in ("epsxx", "epsyy", "epszz", "epsxy", "epsyz", "epsxz")
    },
    **{name: "velocity" for name in ("v1", "v2", "v3")},
    **{name: "fluid velocity" for name in ("v1_f", "v2_f", "v3_f")},
    **{name: "traction" for name in ("T_s", "T_d", "P_n")},
    **{name: "initial traction" for name in ("Ts0", "Td0", "Pn0")},
    **{name: "slip" for name in ("Sls", "Sld")},
    **{name: "slip rate" for name in ("SRs", "SRd")},
}


def component_group(column: str) -> str:
    """The vector or tensor a receiver column is a component of, per simulation.

    A fused simulation appends its index to every column (s_xx3, v13) or, on the
    fault, a dash and its number (SRs-4); the components of one simulation form
    one group. A column that is no component is a group of its own.
    """
    for end in range(len(column), 0, -1):
        name, suffix = column[:end], column[end:]
        if name in _COMPONENTS and re.fullmatch(r"(-?\d+)?", suffix):
            return _COMPONENTS[name] + suffix
    return column


def compare_receiver_columns(
    sim_receiver: pd.DataFrame, ref_receiver: pd.DataFrame, label: str
) -> dict[str, float]:
    """Compare all columns present in the reference receiver against the simulated one.

    Returns a dict mapping column name -> relative L2 error (or absolute if ref is ~zero).
    The components of a vector or a tensor are relative to the largest reference norm
    among them: a component that stays at zero, like the shear stress in water, holds
    rounding noise only, and relative to that noise any other rounding would look
    like a change of order one.
    """
    time = ref_receiver["Time"].values
    columns = [col for col in ref_receiver.columns if col != "Time"]
    ref_norms = {
        col: np.sqrt(trapz_func(ref_receiver[col].values ** 2, x=time))
        for col in columns
    }
    scale = {}
    for col, norm in ref_norms.items():
        group = component_group(col)
        scale[group] = max(scale.get(group, 0.0), norm)
    errors = {}
    for col in columns:
        if col not in sim_receiver.columns:
            print(f"Warning: column '{col}' missing in simulated output for {label}")
            errors[col] = float("inf")
            continue
        diff_col = sim_receiver[col].values - ref_receiver[col].values
        diff_norm = np.sqrt(trapz_func(diff_col**2, x=time))
        ref_norm = scale[component_group(col)]
        errors[col] = (
            float(diff_norm / ref_norm) if ref_norm > 1e-10 else float(diff_norm)
        )
    return errors


def receiver_diff(
    sim_receiver: pd.DataFrame,
    ref_receiver: pd.DataFrame,
    index: int,
    file_type: str = "receiver",
) -> dict[str, float]:
    """
    Checks if the receivers have same time axis, and returns the relative L2 errors
    """
    assert len(sim_receiver) == len(ref_receiver), (
        f"Record count mismatch at {file_type} {index}: "
        f"{len(sim_receiver)} vs {len(ref_receiver)} samples"
    )
    max_time_diff = np.max(
        np.abs(sim_receiver["Time"].values - ref_receiver["Time"].values),
        initial=0.0,
    )
    assert (
        max_time_diff < 1e-6
    ), f"Record time mismatch at {file_type} {index}: max |Δt| = {max_time_diff:.3e}"

    return compare_receiver_columns(
        sim_receiver, ref_receiver, f"{file_type}-{index:05d}"
    )


def find_all_receivers(
    directory: str, prefix: str, file_type: str = "receiver"
) -> np.typing.NDArray[np.int_]:
    """
    Returns list of receivers in the directory with a particular prefix
    """
    file_candidates = glob.glob(f"{directory}/{prefix}-{file_type}-*.dat")

    extract_id = re.compile(r".+/\w+-\w+-(\d+)(?:-\d+)?\.dat$")
    receiver_ids = []
    for fn in file_candidates:
        extract_id_result = extract_id.search(fn)
        if extract_id_result:
            receiver_ids.append(int(extract_id_result.group(1)))
    return np.array(sorted(list(set(receiver_ids))))


def report_errors(
    label: str, all_errors: dict[int, dict[str, float]], epsilon: float
) -> tuple[bool, dict[str, float]]:
    """Print a per-receiver × per-column error table.

    Returns ``(exceeded, per_column_max)`` where ``per_column_max`` maps each
    column to its worst error across all receivers -- used to build the
    machine-readable summary regardless of pass/fail.
    """
    if not all_errors:
        return False, {}

    # Build a DataFrame: rows = receiver IDs, columns = all unique column names seen
    all_cols = sorted({col for errs in all_errors.values() for col in errs})
    df = pd.DataFrame(index=sorted(all_errors.keys()), columns=all_cols, dtype=float)
    for rec_id, errs in all_errors.items():
        for col, val in errs.items():
            df.loc[rec_id, col] = val

    print(f"\nRelative L2 error of quantities at {label}:")
    print(df.to_string())

    per_column_max = {col: float(np.nanmax(df[col].to_numpy())) for col in df.columns}

    exceeded = False
    for col in df.columns:
        broken = df.index[df[col] > epsilon].tolist()
        if broken:
            print(f"  '{col}' exceeds relative error of {epsilon} at {label} {broken}")
            exceeded = True
    return exceeded, per_column_max


def main():
    parser = argparse.ArgumentParser(description="Compare two sets of receivers.")
    parser.add_argument("output", type=str)
    parser.add_argument("output_ref", type=str, nargs="?", default=None)
    parser.add_argument(
        "--list-quantities",
        action="store_true",
        help="Print the output's quantity names as a JSON array and exit "
        "(output_ref is not needed).",
    )
    parser.add_argument("--epsilon", type=float, default=0.01)
    parser.add_argument("--prefix", type=str, default="tpv", required=False)
    parser.add_argument(
        "--report-json",
        type=str,
        default=None,
        help="Write a machine-readable summary of the achieved errors to this "
        "path (always written, regardless of pass/fail). Does not affect the "
        "exit code.",
    )
    args = parser.parse_args()

    if args.list_quantities:
        # Structural: the report keys are f"{file_type}:{col}" for every column
        # except Time (see compare_receiver_columns), read from one file per type.
        import json

        names = []
        for file_type in ("receiver", "faultreceiver"):
            receivers = load_receivers(args.output, args.prefix, file_type)
            if not receivers:
                continue
            for col in receivers[min(receivers)].columns:
                if col != "Time":
                    names.append(f"{file_type}:{col}")
        print(json.dumps(sorted(names)))
        return

    if args.output_ref is None:
        parser.error("output_ref is required unless --list-quantities is given")

    ANY_FAILURE = False
    quantities: dict[str, float] = {}
    for file_type in ("receiver", "faultreceiver"):
        label = f"{file_type}s"
        # either side may be text files or an HDF5 file; both read into the same layout
        sim = load_receivers(args.output, args.prefix, file_type)
        ref = load_receivers(args.output_ref, args.prefix, file_type)
        missing = sorted(set(ref) - set(sim))
        assert not missing, f"some {label} IDs are missing: {missing}"
        errors = {
            index: receiver_diff(sim[index], ref[index], index, file_type=file_type)
            for index in sorted(ref)
        }
        exceeded, per_column_max = report_errors(label, errors, args.epsilon)
        ANY_FAILURE |= exceeded
        # Prefix by receiver type so "receiver" and "faultreceiver" don't collide.
        for col, val in per_column_max.items():
            quantities[f"{file_type}:{col}"] = val

    if args.report_json is not None:
        write_report_json(
            args.report_json, "receiver", args.epsilon, not ANY_FAILURE, quantities
        )

    sys.exit(1 if ANY_FAILURE else 0)


if __name__ == "__main__":
    main()
