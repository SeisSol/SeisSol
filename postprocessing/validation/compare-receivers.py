#!/usr/bin/env python3

# SPDX-FileCopyrightText: 2022 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

import argparse
import glob
import re
import sys

import numpy as np
import pandas as pd

if hasattr(np, "trapezoid"):
    trapz_func = np.trapezoid
else:
    trapz_func = np.trapz


# Maps legacy variable names (written by older SeisSol versions) to current names.
_LEGACY_NAMES = {
    "xx": "s_xx", "yy": "s_yy", "zz": "s_zz",
    "xy": "s_xy", "xz": "s_xz", "yz": "s_yz",
    "u": "v1", "v": "v2", "w": "v3", "-p": "pprime",
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

        # since dr-cpp merge, fault receiver files start writing at Time=0
        # (before they were writing at Time=dt)
        # We then skip the first timestep written if Time = 0
        t0 = float(lines[first_row].split()[0])
        is_fault_receiver = "faultreceiver" in filename
        if t0 == 0 and is_fault_receiver:
            first_row += 1
    receiver = pd.read_csv(filename, header=None, skiprows=first_row, sep=r"\s+")
    receiver.columns = normalize_variable_names(variables)
    return receiver


def compare_receiver_columns(
    sim_receiver: pd.DataFrame, ref_receiver: pd.DataFrame, label: str
) -> dict[str, float]:
    """Compare all columns present in the reference receiver against the simulated one.

    Returns a dict mapping column name -> relative L2 error (or absolute if ref is ~zero).
    """
    time = ref_receiver["Time"].values
    errors = {}
    for col in ref_receiver.columns:
        if col == "Time":
            continue
        if col not in sim_receiver.columns:
            print(f"Warning: column '{col}' missing in simulated output for {label}")
            errors[col] = float("inf")
            continue
        ref_col = ref_receiver[col].values
        diff_col = sim_receiver[col].values - ref_col
        ref_norm = np.sqrt(trapz_func(ref_col**2, x=time))
        diff_norm = np.sqrt(trapz_func(diff_col**2, x=time))
        errors[col] = float(diff_norm / ref_norm) if ref_norm > 1e-10 else float(diff_norm)
    return errors


def receiver_diff(args: argparse.Namespace, index: int, file_type: str = "receiver") -> dict[str, float]:
    """
    Checks if the receivers have same time axis, and returns the relative L2 errors
    """
    sim_files = glob.glob(f"{args.output}/{args.prefix}-{file_type}-{index:05d}*.dat")
    ref_files = glob.glob(f"{args.output_ref}/{args.prefix}-{file_type}-{index:05d}*.dat")

    # allow copy layer receivers
    assert len(sim_files) >= 1
    assert len(ref_files) == 1

    sim_receiver = read_receiver(sim_files[0])
    ref_receiver = read_receiver(ref_files[0])

    # Fault receivers may have the t=0 row dropped; reset to ensure a clean 0-based
    # integer index before comparing. Non-fault receivers are never row-dropped so
    # their index is already clean.
    sim_receiver.reset_index(drop=True, inplace=True)
    ref_receiver.reset_index(drop=True, inplace=True)

    max_time_diff = np.max(np.abs(sim_receiver["Time"].values - ref_receiver["Time"].values))
    assert max_time_diff < 1e-6, (
        f"Record time mismatch at {file_type} {index}: max |Δt| = {max_time_diff:.3e}"
    )

    return compare_receiver_columns(sim_receiver, ref_receiver, f"{file_type}-{index:05d}")


def find_all_receivers(directory: str, prefix: str, file_type: str = "receiver") -> np.typing.NDArray[np.int_]:
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


def report_errors(label: str, all_errors: dict[int, dict[str, float]], epsilon: float) -> bool:
    """Print a per-receiver × per-column error table and return True if any exceed epsilon."""
    if not all_errors:
        return False

    # Build a DataFrame: rows = receiver IDs, columns = all unique column names seen
    all_cols = sorted({col for errs in all_errors.values() for col in errs})
    df = pd.DataFrame(index=sorted(all_errors.keys()), columns=all_cols, dtype=float)
    for rec_id, errs in all_errors.items():
        for col, val in errs.items():
            df.loc[rec_id, col] = val

    print(f"\nRelative L2 error of quantities at {label}:")
    print(df.to_string())

    exceeded = False
    for col in df.columns:
        broken = df.index[df[col] > epsilon].tolist()
        if broken:
            print(f"  '{col}' exceeds relative error of {epsilon} at {label} {broken}")
            exceeded = True
    return exceeded

def main():
    parser = argparse.ArgumentParser(description="Compare two sets of receivers.")
    parser.add_argument("output", type=str)
    parser.add_argument("output_ref", type=str)
    parser.add_argument("--epsilon", type=float, default=0.01)
    parser.add_argument("--prefix", type=str, default="tpv", required=False)
    args = parser.parse_args()

    ANY_FAILURE = False
    for file_type in ("receiver", "faultreceiver"):
        label = f"{file_type}s"
        sim_ids = find_all_receivers(args.output, args.prefix, file_type)
        ref_ids = find_all_receivers(args.output_ref, args.prefix, file_type)
        ids = np.intersect1d(sim_ids, ref_ids)
        assert len(ids) == len(ref_ids), (
            f"some {label} IDs are missing: {ids} vs {ref_ids}"
        )
        errors = {index: receiver_diff(args, index, file_type=file_type) for index in ref_ids}
        ANY_FAILURE |= report_errors(label, errors, args.epsilon)

    sys.exit(1 if ANY_FAILURE else 0)

if __name__ == "__main__":
    main()
