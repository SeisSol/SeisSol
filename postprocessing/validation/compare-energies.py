#!/usr/bin/env python3


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
            max_index = 0
            for c in df.columns:
                current_index = int(c[-1])
                max_index = current_index if current_index > max_index else max_index
        else:
            max_index = df["simulation_index"].max()
        return max_index + 1
    except:
        return -1


def get_sub_simulation(df, fused_index):
    if "simulation_index" not in df.columns:
        # if there's no simulation_index, we're not fused
        return 0
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

    relative_difference_larger_eps = (relative_difference.iloc[1:, :] > epsilon).values
    return np.any(relative_difference_larger_eps)


if __name__ == "__main__":
    import argparse
    import numpy as np
    import sys
    import pandas as pd

    parser = argparse.ArgumentParser(description="Compare energy output csv files.")
    parser.add_argument("energy", type=str)
    parser.add_argument("energy_ref", type=str)
    parser.add_argument("--epsilon", type=float, default=0.01)

    args = parser.parse_args()

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

    if number_of_fused_sims < 0:
        result = perform_check(energy, energy_ref, args.epsilon)
        if result:
            sys.exit(1)

    else:
        for fused_index in range(number_of_fused_sims):
            energy_f = get_sub_simulation(energy, fused_index)
            energy_f = energy_f[relevant_quantities]

            energy_ref_f = get_sub_simulation(energy_ref, fused_index)
            energy_ref_f = energy_ref_f[relevant_quantities]
            result = perform_check(energy_f, energy_ref_f, args.epsilon)
            if result:
                sys.exit(1)
