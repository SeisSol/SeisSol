#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import sys
import os
import re

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare two sets of receivers.")
    parser.add_argument("output", type=str)
    parser.add_argument("output_ref", type=str)
    parser.add_argument("--epsilon", type=float, default=0.01)
    parser.add_argument(
        "--mode", type=str, default="rs", required=False, choices=["rs", "lsw", "tp"]
    )
    args = parser.parse_args()

    def velocity_norm(receiver):
        if (
            "u" in receiver.columns
            and "v" in receiver.columns
            and "w" in receiver.columns
        ):
            return np.sqrt(receiver["u"] ** 2 + receiver["v"] ** 2 + receiver["w"] ** 2)
        elif (
            "v1" in receiver.columns
            and "v2" in receiver.columns
            and "v3" in receiver.columns
        ):
            return np.sqrt(
                receiver["v1"] ** 2 + receiver["v2"] ** 2 + receiver["v3"] ** 2
            )
        else:
            raise ValueError("Could not find velocities in off-fault receiver.")

    def stress_norm(receiver):
        return np.sqrt(
            receiver["xx"] ** 2
            + receiver["yy"] ** 2
            + receiver["zz"] ** 2
            + receiver["xy"] ** 2
            + receiver["yz"] ** 2
            + receiver["xz"] ** 2
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
        return np.sqrt(
            receiver["Td0"] ** 2 + receiver["Ts0"] ** 2 + receiver["Pn0"] ** 2
        )

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
        receiver = pd.read_csv(filename, header=None, skiprows=first_row, sep="\s+")
        receiver.columns = variables
        return receiver

    def receiver_diff(args, i):
        sim_receiver = read_receiver(f"{args.output}/tpv-receiver-{i:05d}-00000.dat")
        ref_receiver = read_receiver(
            f"{args.output_ref}/tpv-receiver-{i:05d}-00000.dat"
        )
        # both receivers must have the same time axis
        assert np.max(np.abs(sim_receiver["Time"] - ref_receiver["Time"])) < 1e-14
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
        ref_filename = f"{args.output_ref}/tpv-faultreceiver-{i:05d}-00000.dat"
        ref_receiver = read_receiver(ref_filename)

        sim_filename = f"{args.output}/tpv-faultreceiver-{i:05d}-00000.dat"
        sim_receiver = read_receiver(sim_filename).iloc[1:]
        sim_receiver.reset_index(drop=True, inplace=True)

        # both receivers must have the same time axis
        assert np.max(np.abs(sim_receiver["Time"] - ref_receiver["Time"])) < 1e-14
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

    def find_all_receivers(directory, faultreceiver=False):
        if faultreceiver:
            receiver_re = re.compile("tpv-faultreceiver-(\d+)-\d+\.dat")
        else:
            receiver_re = re.compile("tpv-receiver-(\d+)-\d+\.dat")
        receiver_ids = []
        for fn in os.listdir(directory):
            match = receiver_re.match(fn)
            if match:
                receiver_ids.append(int(match.groups()[0]))
        return np.array(sorted(receiver_ids))

    sim_faultreceiver_ids = find_all_receivers(args.output, True)
    ref_faultreceiver_ids = find_all_receivers(args.output_ref, True)
    faultreceiver_ids = np.intersect1d(sim_faultreceiver_ids, ref_faultreceiver_ids)
    # Make sure, we actually compare some faultreceivers
    assert len(faultreceiver_ids) == len(ref_faultreceiver_ids)

    sim_receiver_ids = find_all_receivers(args.output, False)
    ref_receiver_ids = find_all_receivers(args.output_ref, False)
    receiver_ids = np.intersect1d(sim_receiver_ids, ref_receiver_ids)
    # Make sure, we actually compare some receivers
    assert len(receiver_ids) == len(ref_receiver_ids)

    receiver_errors = pd.DataFrame(index=receiver_ids, columns=["velocity", "stress"])
    for i in receiver_ids:
        receiver_errors.loc[i, :] = receiver_diff(args, i)
    print("Relative L2 error of the different quantities at the different receivers")
    print(receiver_errors)

    for q in receiver_errors.columns:
        broken_receivers = receiver_errors.index[
            receiver_errors[q] > args.epsilon
        ].tolist()
        print(
            f"{q} exceeds relative error of {args.epsilon} at receveivers {broken_receivers}"
        )

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
    for i in faultreceiver_ids:
        local_errors = faultreceiver_diff(args, i, quantities)
        faultreceiver_errors.loc[i, :] = local_errors.loc[i, :]

    print("")
    print(
        "Relative L2 error of the different quantities at the different faultreceivers"
    )
    print(faultreceiver_errors)

    for q in faultreceiver_errors.columns:
        broken_faultreceivers = faultreceiver_errors.index[
            faultreceiver_errors[q] > args.epsilon
        ].tolist()
        print(
            f"{q} exceeds relative error of {args.epsilon} at faultreceveivers {broken_faultreceivers}"
        )

    if (receiver_errors > args.epsilon).any().any() or (
        faultreceiver_errors > args.epsilon
    ).any().any():
        sys.exit(1)
