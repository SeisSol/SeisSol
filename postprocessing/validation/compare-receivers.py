#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import sys
import os
import re
import glob


def velocity_norm(receiver):
    assert (
        "v1" in receiver.columns
        and "v2" in receiver.columns
        and "v3" in receiver.columns
    )
    return np.sqrt(receiver["v1"] ** 2 + receiver["v2"] ** 2 + receiver["v3"] ** 2)


def stress_norm(receiver):
    return np.sqrt(
        receiver["s_xx"] ** 2
        + receiver["s_yy"] ** 2
        + receiver["s_zz"] ** 2
        + receiver["s_xy"] ** 2
        + receiver["s_yz"] ** 2
        + receiver["s_xz"] ** 2
    )


def absolute_slip_norm(receiver):
    return np.sqrt(receiver["ASl"] ** 2)


def dynstress_norm(receiver):
    return np.sqrt(receiver["DS"] ** 2)


def friction_coefficient_norm(receiver):
    return np.sqrt(receiver["Mud"] ** 2)


def peak_sliprate_norm(receiver):
    return np.sqrt(receiver["PSR"] ** 2)


def pressure_norm(receiver):
    return np.sqrt(receiver["P_f"] ** 2)


def traction_norm(receiver):
    return np.sqrt(receiver["Td0"] ** 2 + receiver["Ts0"] ** 2 + receiver["Pn0"] ** 2)


def rupture_time_norm(receiver):
    return np.sqrt(receiver["RT"] ** 2)


def sliprate_norm(receiver):
    return np.sqrt(receiver["SRs"] ** 2 + receiver["SRd"] ** 2)


def slip_norm(receiver):
    return np.sqrt(receiver["Sls"] ** 2 + receiver["Sld"] ** 2)


def statevariable_norm(receiver):
    return np.sqrt(receiver["StV"] ** 2)


def temperature_norm(receiver):
    return np.sqrt(receiver["Tmp"] ** 2)


def rupture_velocity_norm(receiver):
    return np.sqrt(receiver["Vr"] ** 2)


def normal_velocity_norm(receiver):
    return np.sqrt(receiver["u_n"] ** 2)


def integrate_in_time(time, samples):
    return np.trapz(samples, x=time)


def integrate_quantity_in_time(receiver, quantity):
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
    return integrate_in_time(receiver["Time"], quantity_to_norm[quantity](receiver))


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
        while lines[first_row][0] == "#":
            first_row += 1

        # since dr-cpp merge, fault receiver files start writing at Time=0
        # (before they were writing at Time=dt)
        # We then skip the first timestep written if Time = 0
        t0 = float(lines[first_row].split()[0])
        isFaultReceiver = "faultreceiver" in filename
        if t0 == 0 and isFaultReceiver:
            first_row += 1
    receiver = pd.read_csv(filename, header=None, skiprows=first_row, sep="\s+")

    def replace(x, y, l):
        if x in l:
            x_index = l.index(x)
            l[x_index] = y
        return l

    # Accomodate variable name changes
    variables = replace("xx", "s_xx", variables)
    variables = replace("yy", "s_yy", variables)
    variables = replace("zz", "s_zz", variables)
    variables = replace("xy", "s_xy", variables)
    variables = replace("xz", "s_xz", variables)
    variables = replace("yz", "s_yz", variables)
    variables = replace("u", "v1", variables)
    variables = replace("v", "v2", variables)
    variables = replace("w", "v3", variables)
    
    receiver.columns = variables
    return receiver


def receiver_diff(args, i):
    sim_files = glob.glob(f"{args.output}/{args.prefix}-receiver-{i:05d}-*.dat")
    ref_files = glob.glob(f"{args.output_ref}/{args.prefix}-receiver-{i:05d}-*.dat")

    # allow copy layer receivers
    assert len(sim_files) >= 1
    assert len(ref_files) == 1

    sim_receiver = read_receiver(sim_files[0])
    ref_receiver = read_receiver(ref_files[0])
    # both receivers must have the same time axis
    assert np.max(np.abs(sim_receiver["Time"] - ref_receiver["Time"])) < 1e-7, f'Record time mismatch: {sim_receiver["Time"]} vs. {ref_receiver["Time"]}'
    time = sim_receiver["Time"]
    difference = sim_receiver - ref_receiver

    ref_velocity_norm = integrate_in_time(time, velocity_norm(ref_receiver))
    diff_velocity_norm = integrate_in_time(time, velocity_norm(difference))

    ref_stress_norm = integrate_in_time(time, stress_norm(ref_receiver))
    diff_stress_norm = integrate_in_time(time, stress_norm(difference))

    return (
        diff_velocity_norm / ref_velocity_norm,
        diff_stress_norm / ref_stress_norm,
    )

def faultreceiver_diff(args, i, quantities):
    sim_files = glob.glob(f"{args.output}/{args.prefix}-faultreceiver-{i:05d}-*.dat")
    ref_files = glob.glob(
        f"{args.output_ref}/{args.prefix}-faultreceiver-{i:05d}-*.dat"
    )
    # allow copy layer receivers
    assert len(sim_files) >= 1
    assert len(ref_files) == 1
    sim_receiver = read_receiver(sim_files[0])
    ref_receiver = read_receiver(ref_files[0])

    sim_receiver.reset_index(drop=True, inplace=True)

    # both receivers must have the same time axis
    assert np.max(np.abs(sim_receiver["Time"] - ref_receiver["Time"])) < 1e-7, f'Record time mismatch: {sim_receiver["Time"]} vs. {ref_receiver["Time"]}'
    time = sim_receiver["Time"]
    difference = sim_receiver - ref_receiver
    # We still want to use the same time and not the difference in time steps.
    difference["Time"] = ref_receiver["Time"]

    errors = pd.DataFrame(index=[i], columns=quantities)

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

    for quantity_name in possible_quantity_names:
        if quantity_name in quantities:
            ref_norm = integrate_quantity_in_time(ref_receiver, quantity_name)
            diff_norm = integrate_quantity_in_time(difference, quantity_name)
            errors.loc[i, quantity_name] = diff_norm / ref_norm

    return errors


def find_all_receivers(directory, prefix, faultreceiver=False):
    if faultreceiver:
        file_candidates = glob.glob(f"{directory}/{prefix}-faultreceiver-*-*.dat")
    else:
        file_candidates = glob.glob(f"{directory}/{prefix}-receiver-*-*.dat")
    extract_id = re.compile(r".+/\w+-\w+-(\d+)-\d+.dat$")

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
        assert len(faultreceiver_ids) == len(ref_faultreceiver_ids), f'some on-fault receiver IDs are missing: {faultreceiver_ids} vs {ref_faultreceiver_ids}'
    else:
        sim_faultreceiver_ids = []
        ref_faultreceiver_ids = []
        faultreceiver_ids = []

    sim_receiver_ids = find_all_receivers(args.output, args.prefix, False)
    ref_receiver_ids = find_all_receivers(args.output_ref, args.prefix, False)
    receiver_ids = np.intersect1d(sim_receiver_ids, ref_receiver_ids)
    # Make sure, we actually compare some receivers
    assert len(receiver_ids) == len(ref_receiver_ids), f'some off-fault receiver IDs are missing: {receiver_ids} vs {ref_receiver_ids}'

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
