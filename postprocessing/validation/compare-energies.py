#!/usr/bin/env python3

# SPDX-FileCopyrightText: 2022 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

import argparse
import sys
import re

import numpy as np
import pandas as pd

from validation_report import write_report_json


def pivot_if_necessary(df):
    if "variable" in df:
        # the format of the energy output changed following PR #773 (02.2023), allowing
        # to compute volume energies less frequently
        indices = ["time"]
        if "simulation_index" in df.columns:
            indices += ["simulation_index"]
        return (
            df.pivot_table(index=indices, columns="variable", values="measurement")
            .reset_index()
            .astype({idx: "int64" for idx in indices})
        )
    else:
        return df


def get_number_of_fused_sims(df):
    try:
        if "simulation_index" not in df.columns:
            # for historic reasons we want a negative number (e.g. -1) to be returned for non-fused; thus start with -1 here
            max_index = -2
            for c in df.columns:
                idxres = re.match(r'.*(\d+)$', c)
                if idxres:
                    current_index = int(idxres.group(1))
                    max_index = current_index if current_index > max_index else max_index
        else:
            max_index = df["simulation_index"].max()
        return max_index + 1
    except:
        return -1


def get_sub_simulation(df, fused_index):
    if "simulation_index" not in df.columns:
        # if there's no simulation_index, we're not fused
        return df.loc[:, :].reset_index()
    else:
        is_subsim = df["simulation_index"] == fused_index
        return df.loc[is_subsim, :].reset_index()


def perform_check(energy, energy_ref, epsilon):
    print("Energies")
    print(energy.to_string())
    print("Energies reference")
    print(energy_ref.to_string())
    relative_difference = ((energy - energy_ref).abs() / energy_ref).iloc[1:, :]
    print("Relative difference")
    print(relative_difference.to_string())

    # NOTE: ``relative_difference`` already drops the t=0 row (``.iloc[1:]`` above).
    # It must NOT be sliced a second time here -- doing so silently excluded the
    # first post-t0 timestep from the pass/fail decision.
    relative_difference_larger_eps = (relative_difference > epsilon).values
    exceeded = bool(np.any(relative_difference_larger_eps))
    # Also hand back the raw matrix so the caller can build a machine-readable
    # summary of the *achieved* errors (independent of pass/fail).
    return exceeded, relative_difference


def _numeric_column_maxima(rel_diff):
    """Per-column max of |relative difference|, dropping non-numeric/index columns.

    Returns a plain ``dict[str, float]``. NaNs (e.g. from 0/0) are ignored;
    +inf (from x/0) is preserved so it still shows up as a violation.
    """
    maxima = {}
    for col in rel_diff.columns:
        series = pd.to_numeric(rel_diff[col], errors="coerce")
        if series.notna().any():
            maxima[str(col)] = float(np.nanmax(series.to_numpy()))
    # 'time' compares to itself and is ~0; it carries no signal, drop it.
    maxima.pop("time", None)
    return maxima


def main():

    parser = argparse.ArgumentParser(description="Compare energy output csv files.")
    parser.add_argument("energy", type=str)
    parser.add_argument("energy_ref", type=str, nargs="?", default=None)
    parser.add_argument(
        "--list-quantities",
        action="store_true",
        help="Print the output's quantity names as a JSON array and exit "
        "(energy_ref is not needed).",
    )
    parser.add_argument("--epsilon", type=float, default=0.01)
    parser.add_argument(
        "--report-json",
        type=str,
        default=None,
        help="Write a machine-readable summary of the achieved errors to this "
        "path (always written, regardless of pass/fail). Does not affect the "
        "exit code.",
    )

    args = parser.parse_args()

    # In list mode we compare the output against itself, purely to reuse the
    # exact key-generation path below (the errors come out zero and are ignored).
    # The comparison is chatty, so its stdout is swallowed until the JSON is ready.
    _real_stdout = sys.stdout
    if args.list_quantities:
        import io
        if args.energy_ref is None:
            args.energy_ref = args.energy
        sys.stdout = io.StringIO()
    if args.energy_ref is None:
        parser.error("energy_ref is required unless --list-quantities is given")

    relevant_quantities = [
        "elastic_energy",
        "elastic_kinetic_energy",
        "total_frictional_work",
        "static_frictional_work",
        "seismic_moment",
    ]
    energy = pd.read_csv(args.energy)
    energy = pivot_if_necessary(energy)
    energy_ref = pd.read_csv(args.energy_ref)
    energy_ref = pivot_if_necessary(energy_ref)

    sims = get_number_of_fused_sims(energy)
    sims_ref = get_number_of_fused_sims(energy_ref)
    assert sims >= sims_ref, f"Simulation count mismatch: {sims} vs. {sims_ref}"
    number_of_fused_sims = get_number_of_fused_sims(energy)

    # restrict to present quantities only (the rest is assumed to be zero or undefined)
    relevant_quantities = list(
        set(energy.columns) & set(energy_ref.columns) & set(relevant_quantities)
    )
    relevant_quantities = list(sorted(relevant_quantities))

    failed = False
    quantities: dict[str, float] = {}  # column -> worst error seen across sub-sims

    if number_of_fused_sims < 0:
        exceeded, rel_diff = perform_check(energy, energy_ref, args.epsilon)
        failed |= exceeded
        for col, val in _numeric_column_maxima(rel_diff).items():
            quantities[col] = max(quantities.get(col, 0.0), val)
    else:
        for fused_index in range(number_of_fused_sims):
            energy_f = get_sub_simulation(energy, fused_index)
            energy_f = energy_f[relevant_quantities]

            energy_ref_f = get_sub_simulation(energy_ref, fused_index)
            energy_ref_f = energy_ref_f[relevant_quantities]
            exceeded, rel_diff = perform_check(energy_f, energy_ref_f, args.epsilon)
            failed |= exceeded
            for col, val in _numeric_column_maxima(rel_diff).items():
                key = f"{col}[{fused_index}]" if number_of_fused_sims > 1 else col
                quantities[key] = max(quantities.get(key, 0.0), val)

    if args.list_quantities:
        import json
        sys.stdout = _real_stdout
        print(json.dumps(sorted(quantities)))
        return

    if args.report_json is not None:
        write_report_json(
            args.report_json, "energy", args.epsilon, not failed, quantities
        )

    if failed:
        sys.exit(1)

if __name__ == "__main__":
    main()
