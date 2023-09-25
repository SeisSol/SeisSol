#!/usr/bin/env python3


def pivot_if_necessary(df):
    if "variable" in df:
        # the format of the energy output changed following PR #773 (02.2023), allowing
        # to compute volume energies less frequently
        return df.pivot_table(index="time", columns="variable", values="measurement")
    else:
        return df


def get_number_of_fused_sims(df):
    if "simulation_index" not in df.columns:
        max_index = 0
        for c in df.columns:
            current_index = int(c[-1])
            max_index = current_index if current_index > max_index else max_index
    else:
        max_index = df["simulation_index"].max()
    return max_index + 1


def get_sub_simulation(df, fused_index):
    if "simulation_index" not in df.columns:
        sub_df = pd.DataFrame()
        for c in df.columns:
            if int(c[-1]) == fused_index:
                sub_df[c[:-1]] = df[c]
        return sub_df
    else:
        is_subsim = df["simulation_index"] == fused_index
        return df.loc[is_subsim, :].reset_index()


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

    assert get_number_of_fused_sims(energy) == get_number_of_fused_sims(energy_ref)
    number_of_fused_sims = get_number_of_fused_sims(energy)

    for fused_index in range(number_of_fused_sims):
        energy_f = get_sub_simulation(energy, fused_index)
        energy_f = energy_f[relevant_quantities]

        energy_ref_f = get_sub_simulation(energy_ref, fused_index)
        energy_ref_f = energy_ref_f[relevant_quantities]
        print(f"Fused sim = {fused_index}")
        print("Energies")
        print(energy_f)
        print("Energies reference")
        print(energy_ref_f)
        relative_difference = ((energy_f - energy_ref_f).abs() / energy_ref_f).iloc[
            1:, :
        ]
        print("Relative difference")
        print(relative_difference)

        relative_difference_larger_eps = (
            relative_difference.iloc[1:, :] > args.epsilon
        ).values
        if np.any(relative_difference_larger_eps):
            sys.exit(1)
