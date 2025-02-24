#!/usr/bin/env python3


def pivot_if_necessary(df):
    if "variable" in df:
        # the format of the energy output changed following PR #773 (02.2023), allowing
        # to compute volume energies less frequently
        return df.pivot_table(index="time", columns="variable", values="measurement")
    else:
        return df


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
    energy = energy[relevant_quantities]

    energy_ref = pd.read_csv(args.energy_ref)
    energy_ref = pivot_if_necessary(energy_ref)
    energy_ref = energy_ref[relevant_quantities]
    print("Energies")
    print(energy.to_string())
    print("Energies reference")
    print(energy_ref.to_string())
    relative_difference = ((energy - energy_ref).abs() / energy_ref).iloc[1:, :]
    print("Relative difference")
    print(relative_difference.to_string())

    relative_difference_larger_eps = (
        relative_difference.iloc[1:, :] > args.epsilon
    ).values
    if np.any(relative_difference_larger_eps):
        sys.exit(1)
