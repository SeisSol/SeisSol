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

if hasattr(np, "trapezoid"):
    trapz_func = np.trapezoid
else:
    trapz_func = np.trapz


def velocity_norm(receiver, fused_index=""):
    names = [f"v1{fused_index}", f"v2{fused_index}", f"v3{fused_index}"]

    assert (
        names[0] in receiver.columns
        and names[1] in receiver.columns
        and names[2] in receiver.columns
    )

    return np.sqrt(
        receiver[names[0]] ** 2 + receiver[names[1]] ** 2 + receiver[names[2]] ** 2
    )


def stress_norm(receiver, fused_index=""):
    names = [
        f"s_xx{fused_index}",
        f"s_yy{fused_index}",
        f"s_zz{fused_index}",
        f"s_xy{fused_index}",
        f"s_yz{fused_index}",
        f"s_xz{fused_index}",
    ]
    return np.sqrt(
        receiver[names[0]] ** 2
        + receiver[names[1]] ** 2
        + receiver[names[2]] ** 2
        + receiver[names[3]] ** 2
        + receiver[names[4]] ** 2
        + receiver[names[5]] ** 2
    )


def absolute_slip_norm(receiver, fused_suffix=""):
    return np.sqrt(receiver["ASl" + fused_suffix] ** 2)


def dynstress_norm(receiver, fused_suffix=""):
    return np.sqrt(receiver["DS" + fused_suffix] ** 2)


def friction_coefficient_norm(receiver, fused_suffix=""):
    return np.sqrt(receiver["Mud" + fused_suffix] ** 2)


def peak_sliprate_norm(receiver, fused_suffix=""):
    return np.sqrt(receiver["PSR" + fused_suffix] ** 2)


def pressure_norm(receiver, fused_suffix=""):
    return np.sqrt(receiver["P_f" + fused_suffix] ** 2)


def traction_norm(receiver, fused_suffix=""):
    return np.sqrt(
        receiver["Td0" + fused_suffix] ** 2
        + receiver["Ts0" + fused_suffix] ** 2
        + receiver["Pn0" + fused_suffix] ** 2
    )


def rupture_time_norm(receiver, fused_suffix=""):
    return np.sqrt(receiver["RT" + fused_suffix] ** 2)


def sliprate_norm(receiver, fused_suffix=""):
    return np.sqrt(
        receiver["SRs" + fused_suffix] ** 2 + receiver["SRd" + fused_suffix] ** 2
    )


def slip_norm(receiver, fused_suffix=""):
    return np.sqrt(
        receiver["Sls" + fused_suffix] ** 2 + receiver["Sld" + fused_suffix] ** 2
    )


def statevariable_norm(receiver, fused_suffix=""):
    return np.sqrt(receiver["StV" + fused_suffix] ** 2)


def temperature_norm(receiver, fused_suffix=""):
    return np.sqrt(receiver["Tmp" + fused_suffix] ** 2)


def rupture_velocity_norm(receiver, fused_suffix=""):
    return np.sqrt(receiver["Vr" + fused_suffix] ** 2)


def normal_velocity_norm(receiver, fused_suffix=""):
    return np.sqrt(receiver["u_n" + fused_suffix] ** 2)


def integrate_in_time(time, samples):
    return trapz_func(samples, x=time)


def integrate_quantity_in_time(
    receiver, quantity, fused_index=0, number_of_fused_sims=1
):
    quantity_to_norm = {
        "absolute slip": absolute_slip_norm,
        "friction coefficient": friction_coefficient_norm,
        "peak sliprate": peak_sliprate_norm,
        "traction": traction_norm,
        "rupture time": rupture_time_norm,
        "sliprate": sliprate_norm,
        "slip": slip_norm,
        "rupture velocity": rupture_velocity_norm,
        "normal velocity": normal_velocity_norm,
        "dynstress time": dynstress_norm,
        "state variable": statevariable_norm,
        "pressure": pressure_norm,
        "temperature": temperature_norm,
    }
    fused_suffix = ""
    if number_of_fused_sims > 1:
        fused_suffix += "-" + str(
            fused_index + 1
        )  # +1 because we want to use the fused index as per fault receiver numbering
    return integrate_in_time(
        receiver["Time"], quantity_to_norm[quantity](receiver, fused_suffix)
    )


def get_number_of_fused_sims(columns):
    # a _very_ hacky way to ignore the velocities here
    # (which are by now: v1,v2,v3 without fusing)
    # (with fusing: v1N, v2N, v3N. E.g. v14, v24, v34 for simulation 4)
    extract_fused = re.compile(r"^[^v].*?(\d+)$")
    # omit time
    relevant_columns = columns[1:]
    try:
        max_index = -2
        for c in relevant_columns:
            fused = extract_fused.search(c)
            if fused:
                current_index = int(fused.group(1))
                max_index = current_index if current_index > max_index else max_index
        return max_index + 1
    except Exception as e:
        print(e)
        return -1


def read_receiver(filename):
    with open(filename) as receiver_file:
        # find variable names
        # We expect the header to look like the following
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
        isFaultReceiver = "faultreceiver" in filename
        if t0 == 0 and isFaultReceiver:
            first_row += 1
    receiver = pd.read_csv(filename, header=None, skiprows=first_row, sep=r"\s+")

    def replace(x, y, l, max_fused=-1):
        if max_fused < 0:
            if x in l:
                x_index = l.index(x)
                l[x_index] = y
        else:
            for fused_index in range(max_fused):
                x_ = f"{x}{fused_index}"
                y_ = f"{y}{fused_index}"
                if x_ in l:
                    x_index = l.index(x_)
                    l[x_index] = y_
        return l

    # Accomodate variable name changes
    number_of_fused_sims = get_number_of_fused_sims(variables)
    variables = replace("xx", "s_xx", variables, number_of_fused_sims)
    variables = replace("yy", "s_yy", variables, number_of_fused_sims)
    variables = replace("zz", "s_zz", variables, number_of_fused_sims)
    variables = replace("xy", "s_xy", variables, number_of_fused_sims)
    variables = replace("xz", "s_xz", variables, number_of_fused_sims)
    variables = replace("yz", "s_yz", variables, number_of_fused_sims)
    variables = replace("u", "v1", variables, number_of_fused_sims)
    variables = replace("v", "v2", variables, number_of_fused_sims)
    variables = replace("w", "v3", variables, number_of_fused_sims)

    receiver.columns = variables
    return receiver


def receiver_diff(args, i):
    sim_files = glob.glob(f"{args.output}/{args.prefix}-receiver-{i:05d}*.dat")
    ref_files = glob.glob(f"{args.output_ref}/{args.prefix}-receiver-{i:05d}*.dat")

    # allow copy layer receivers
    assert len(sim_files) >= 1
    assert len(ref_files) == 1

    sim_receiver = read_receiver(sim_files[0])
    ref_receiver = read_receiver(ref_files[0])
    # both receivers must have the same time axis
    assert (
        np.max(np.abs(sim_receiver["Time"] - ref_receiver["Time"])) < 1e-6
    ), f'Record time mismatch: {sim_receiver["Time"]} vs. {ref_receiver["Time"]}'
    time = sim_receiver["Time"]
    difference = sim_receiver - ref_receiver

    number_of_fused_sims = get_number_of_fused_sims(sim_receiver.columns)
    if number_of_fused_sims < 1:
        print(
            "Setting the number of fused simulations to 1, because the receiver file does not contain any fused simulations."
        )
        number_of_fused_sims = 1

    max_velocity = 0
    max_stress = 0

    for fused_index in range(number_of_fused_sims):
        fused_suffix = f"{fused_index}" if number_of_fused_sims > 1 else ""
        ref_velocity_norm = integrate_in_time(
            time, velocity_norm(ref_receiver, fused_suffix)
        )
        diff_velocity_norm = integrate_in_time(
            time, velocity_norm(difference, fused_suffix)
        )
        rel_velocity_diff = diff_velocity_norm / ref_velocity_norm
        max_velocity = (
            rel_velocity_diff if rel_velocity_diff > max_velocity else max_velocity
        )

        ref_stress_norm = integrate_in_time(
            time, stress_norm(ref_receiver, fused_suffix)
        )
        diff_stress_norm = integrate_in_time(
            time, stress_norm(difference, fused_suffix)
        )
        rel_stress_diff = diff_stress_norm / ref_stress_norm
        max_stress = rel_stress_diff if rel_stress_diff > max_stress else max_stress

    return (
        max_velocity,
        max_stress,
    )


def faultreceiver_diff(args, i, quantities):
    sim_files = glob.glob(f"{args.output}/{args.prefix}-faultreceiver-{i:05d}*.dat")
    ref_files = glob.glob(f"{args.output_ref}/{args.prefix}-faultreceiver-{i:05d}*.dat")
    # allow copy layer receivers
    assert len(sim_files) >= 1
    assert len(ref_files) == 1
    sim_receiver = read_receiver(sim_files[0])
    ref_receiver = read_receiver(ref_files[0])

    sim_receiver.reset_index(drop=True, inplace=True)
    # both receivers must have the same time axis
    assert (
        np.max(np.abs(sim_receiver["Time"] - ref_receiver["Time"])) < 1e-6
    ), f'Record time mismatch: {sim_receiver["Time"]} vs. {ref_receiver["Time"]}'
    time = sim_receiver["Time"]
    difference = sim_receiver - ref_receiver
    number_of_fused_sims = (
        get_number_of_fused_sims(sim_receiver.columns) - 1
    )  # -1 because the numbering of fault receivers is different for fused when compared to off-fault receivers
    if number_of_fused_sims < 1:
        print(
            "Setting the number of fused simulations to 1, because the receiver file does not contain any fused simulations."
        )
        number_of_fused_sims = 1
    # We still want to use the same time and not the difference in time steps.
    difference["Time"] = ref_receiver["Time"]

    errors = pd.DataFrame(0.0, index=[i], columns=quantities)

    possible_quantity_names = [
        "absolute slip",
        "friction coefficient",
        "peak sliprate",
        "traction",
        "rupture time",
        "sliprate",
        "slip",
        "rupture velocity",
        "normal velocity",
        "dynstress time",
        "state variable",
        "pressure",
        "temperature",
    ]

    for fused_index in range(number_of_fused_sims):
        for quantity_name in possible_quantity_names:
            if quantity_name in quantities:
                ref_norm = integrate_quantity_in_time(
                    ref_receiver,
                    quantity_name,
                    fused_index=fused_index,
                    number_of_fused_sims=number_of_fused_sims,
                )
                diff_norm = integrate_quantity_in_time(
                    difference,
                    quantity_name,
                    fused_index=fused_index,
                    number_of_fused_sims=number_of_fused_sims,
                )
                errors.loc[i, quantity_name] = (
                    diff_norm / ref_norm
                    if diff_norm / ref_norm > errors.loc[i, quantity_name]
                    else errors.loc[i, quantity_name]
                )
    return errors


def find_all_receivers(directory, prefix, faultreceiver=False):
    if faultreceiver:
        file_candidates = glob.glob(f"{directory}/{prefix}-faultreceiver-*.dat")
    else:
        file_candidates = glob.glob(f"{directory}/{prefix}-receiver-*.dat")

    extract_id = re.compile(r".+/\w+-\w+-(\d+)(?:-\d+)?\.dat$")
    receiver_ids = []
    for fn in file_candidates:
        extract_id_result = extract_id.search(fn)
        if extract_id_result:
            receiver_ids.append(int(extract_id_result.group(1)))
    return np.array(sorted(list(set(receiver_ids))))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare two sets of receivers.")
    parser.add_argument("output", type=str)
    parser.add_argument("output_ref", type=str)
    parser.add_argument("--epsilon", type=float, default=0.01)
    parser.add_argument(
        "--mode", type=str, default="rs", required=False, choices=["rs", "lsw", "tp"]
    )
    parser.add_argument("--prefix", type=str, default="tpv", required=False)
    parser.add_argument(
        "--compare_fault_receivers", type=bool, default=True, required=False
    )
    args = parser.parse_args()

    if args.compare_fault_receivers:
        sim_faultreceiver_ids = find_all_receivers(args.output, args.prefix, True)
        ref_faultreceiver_ids = find_all_receivers(args.output_ref, args.prefix, True)
        faultreceiver_ids = np.intersect1d(sim_faultreceiver_ids, ref_faultreceiver_ids)
        # Make sure, we actually compare some faultreceivers
        assert len(faultreceiver_ids) == len(
            ref_faultreceiver_ids
        ), f"some on-fault receiver IDs are missing: {faultreceiver_ids} vs {ref_faultreceiver_ids}"
    else:
        sim_faultreceiver_ids = []
        ref_faultreceiver_ids = []
        faultreceiver_ids = []

    sim_receiver_ids = find_all_receivers(args.output, args.prefix, False)
    ref_receiver_ids = find_all_receivers(args.output_ref, args.prefix, False)
    receiver_ids = np.intersect1d(sim_receiver_ids, ref_receiver_ids)
    # Make sure, we actually compare some receivers
    assert len(receiver_ids) == len(
        ref_receiver_ids
    ), f"some off-fault receiver IDs are missing: {receiver_ids} vs {ref_receiver_ids}"

    receiver_errors = pd.DataFrame(index=receiver_ids, columns=["velocity", "stress"])
    for i in ref_receiver_ids:
        receiver_errors.loc[i, :] = receiver_diff(args, i)
    print("Relative L2 error of the different quantities at the different receivers")
    print(receiver_errors.to_string())

    for q in receiver_errors.columns:
        broken_receivers = receiver_errors.index[
            receiver_errors[q] > args.epsilon
        ].tolist()
        print(
            f"{q} exceeds relative error of {args.epsilon} at receivers {broken_receivers}"
        )
        if len(broken_receivers) > 0:
            sys.exit(1)

    quantities = [
        "absolute slip",
        "friction coefficient",
        "peak sliprate",
        "traction",
        "rupture time",
        "sliprate",
        "slip",
        "normal velocity",
    ]
    if args.mode == "lsw":
        quantities.append("dynstress time")
    if args.mode == "rs":
        quantities.append("state variable")
    if args.mode == "tp":
        quantities.append("pressure")
        quantities.append("temperature")

    faultreceiver_errors = pd.DataFrame(index=faultreceiver_ids, columns=quantities)
    for i in ref_faultreceiver_ids:
        local_errors = faultreceiver_diff(args, i, quantities)
        faultreceiver_errors.loc[i, :] = local_errors.loc[i, :]

    print("")
    print(
        "Relative L2 error of the different quantities at the different faultreceivers"
    )
    print(faultreceiver_errors.to_string())

    for q in faultreceiver_errors.columns:
        broken_faultreceivers = faultreceiver_errors.index[
            faultreceiver_errors[q] > args.epsilon
        ].tolist()
        print(
            f"{q} exceeds relative error of {args.epsilon} at faultreceivers {broken_faultreceivers}"
        )

    if (receiver_errors > args.epsilon).any().any() or (
        faultreceiver_errors > args.epsilon
    ).any().any():
        sys.exit(1)

    sys.exit(0)
