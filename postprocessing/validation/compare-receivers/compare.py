import argparse
import numpy as np
import pandas as pd
import sys
import os 
import re

parser = argparse.ArgumentParser(description='Compare two sets of receivers.')
parser.add_argument('output', type=str)
parser.add_argument('output_ref', type=str)
parser.add_argument('--epsilon', type=float, default=0.01, required=False)
parser.add_argument('--mode', type=str, default="rs", required=False, choices=["rs", "lsw", "tp"])
args = parser.parse_args()

def velocity_norm(receiver):
    return np.sqrt(receiver[7]**2 + receiver[8]**2 + receiver[9]**2)

def stress_norm(receiver):
    return np.sqrt(receiver[1]**2 + receiver[2]**2 + receiver[3]**2
            + receiver[4]**2 + receiver[5]**2 + receiver[6]**2)

def absolute_slip_norm(receiver):
    return np.sqrt(receiver["ASl"]**2)

def dynstress_norm(receiver):
    return np.sqrt(receiver["DS"]**2)

def friction_coefficient_norm(receiver):
    return np.sqrt(receiver["Mud"]**2)

def peak_sliprate_norm(receiver):
    return np.sqrt(receiver["PSR"]**2)

def pressure_norm(receiver):
    return np.sqrt(receiver["P_f"]**2)

def traction_norm(receiver):
    return np.sqrt(receiver["Td0"]**2 + receiver["Ts0"]**2 + receiver["Pn0"]**2)

def rupture_time_norm(receiver):
    return np.sqrt(receiver["RT"]**2)

def sliprate_norm(receiver):
    return np.sqrt(receiver["SRs"]**2 + receiver["SRd"]**2)

def slip_norm(receiver):
    return np.sqrt(receiver["Sls"]**2 + receiver["Sld"]**2)

def statevariable_norm(receiver):
    return np.sqrt(receiver["StV"]**2)

def temperature_norm(receiver):
    return np.sqrt(receiver["Tmp"]**2)

def rupture_velocity_norm(receiver):
    return np.sqrt(receiver["Vr"]**2)

def normal_velocity_norm(receiver):
    return np.sqrt(receiver["u_n"]**2)

def integrate_in_time(time, samples):
    return np.trapz(samples, x=time)

def receiver_diff(args, i):
    sim_receiver = pd.read_csv(f"{args.output}/tpv-receiver-0000{i}-00000.dat", header=None, skiprows=5, sep="\s+")
    ref_receiver = pd.read_csv(f"{args.output_ref}/tpv-receiver-0000{i}-00000.dat", header=None, skiprows=5, sep="\s+")
    # both receivers must have the same time axis
    assert(np.max(np.abs(sim_receiver[0] - ref_receiver[0])) < 1e-14)
    time = sim_receiver[0]
    difference = sim_receiver - ref_receiver

    ref_velocity_norm = integrate_in_time(time, velocity_norm(ref_receiver))
    diff_velocity_norm = integrate_in_time(time, velocity_norm(difference))
    
    ref_stress_norm = integrate_in_time(time, stress_norm(ref_receiver))
    diff_stress_norm = integrate_in_time(time, stress_norm(difference))

    return diff_velocity_norm / ref_velocity_norm, diff_stress_norm / ref_stress_norm

def read_initial_stress(filename):
    with open(filename) as receiver_file:
        lines = receiver_file.readlines()
        P_0 = float(lines[5][5:])
        T_s = float(lines[6][5:])
        T_d = float(lines[7][5:])
    return [P_0, T_s, T_d]

def read_faultreceiver(filename):
    with open(filename) as receiver_file:
        # find variable names
        lines = receiver_file.readlines()
        variable_line = lines[1][12:].split(",")
        variables = [s.strip().replace('"', '') for s in variable_line]
        # find first row without comments
        first_row = 2
        while lines[first_row][0] == "#":
            first_row += 1
    receiver = pd.read_csv(filename, header=None, skiprows=first_row, sep="\s+")
    receiver.columns = variables
    return receiver

def faultreceiver_diff(args, i, quantities):
    ref_filename = f"{args.output_ref}/tpv-faultreceiver-{i:05d}-00000.dat"
    ref_receiver = read_faultreceiver(ref_filename)

    sim_filename = f"{args.output}/tpv-faultreceiver-{i:05d}-00000.dat"
    sim_receiver = read_faultreceiver(sim_filename).iloc[1:]
    sim_receiver.reset_index(drop=True, inplace=True)

    # both receivers must have the same time axis
    assert(np.max(np.abs(sim_receiver["Time"] - ref_receiver["Time"])) < 1e-14)
    time = sim_receiver["Time"]
    difference = sim_receiver - ref_receiver

    errors = pd.DataFrame(index=[i], columns=quantities)
    
    if "absolute slip" in quantities:   
        ref_absolute_slip_norm = integrate_in_time(time, absolute_slip_norm(ref_receiver))
        diff_absolute_slip_norm = integrate_in_time(time, absolute_slip_norm(difference))
        errors.loc[i, "absolute slip"] = diff_absolute_slip_norm / ref_absolute_slip_norm

    if "friction coefficient" in quantities:   
        ref_friction_coefficient_norm = integrate_in_time(time, friction_coefficient_norm(ref_receiver))
        diff_friction_coefficient_norm = integrate_in_time(time, friction_coefficient_norm(difference))
        errors.loc[i, "friction coefficient"] = diff_friction_coefficient_norm / ref_friction_coefficient_norm

    if "peak sliprate" in quantities:   
        ref_peak_sliprate_norm = integrate_in_time(time, peak_sliprate_norm(ref_receiver))
        diff_peak_sliprate_norm = integrate_in_time(time, peak_sliprate_norm(difference))
        errors.loc[i, "peak sliprate"] = diff_peak_sliprate_norm / ref_peak_sliprate_norm

    if "traction" in quantities:   
        ref_traction_norm = integrate_in_time(time, traction_norm(ref_receiver))
        diff_traction_norm = integrate_in_time(time, traction_norm(difference))
        errors.loc[i, "traction"] = diff_traction_norm / ref_traction_norm

    if "rupture time" in quantities:   
        ref_rupture_time_norm = integrate_in_time(time, rupture_time_norm(ref_receiver))
        diff_rupture_time_norm = integrate_in_time(time, rupture_time_norm(difference))
        errors.loc[i, "rupture time"] = diff_rupture_time_norm / ref_rupture_time_norm

    if "sliprate" in quantities:   
        ref_sliprate_norm = integrate_in_time(time, sliprate_norm(ref_receiver))
        diff_sliprate_norm = integrate_in_time(time, sliprate_norm(difference))
        errors.loc[i, "sliprate"] = diff_sliprate_norm / ref_sliprate_norm

    if "slip" in quantities:   
        ref_slip_norm = integrate_in_time(time, slip_norm(ref_receiver))
        diff_slip_norm = integrate_in_time(time, slip_norm(difference))
        errors.loc[i, "slip"] = diff_slip_norm / ref_slip_norm

    if "rupture velocity" in quantities:   
        ref_rupture_velocity_norm = integrate_in_time(time, rupture_velocity_norm(ref_receiver))
        diff_rupture_velocity_norm = integrate_in_time(time, rupture_velocity_norm(difference))
        errors.loc[i, "rupture velocity"] = diff_rupture_velocity_norm / ref_rupture_velocity_norm

    if "normal velocity" in quantities:   
        ref_normal_velocity_norm = integrate_in_time(time, normal_velocity_norm(ref_receiver))
        diff_normal_velocity_norm = integrate_in_time(time, normal_velocity_norm(difference))
        errors.loc[i, "normal velocity"] = diff_normal_velocity_norm / ref_normal_velocity_norm

    if "dynstress time" in quantities:   
        ref_dynstress_norm = integrate_in_time(time, dynstress_norm(ref_receiver))
        diff_dynstress_norm = integrate_in_time(time, dynstress_norm(difference))
        errors.loc[i, "dynstress time"] = diff_dynstress_norm / ref_dynstress_norm

    if "state variable" in quantities:   
        ref_statevariable_norm = integrate_in_time(time, statevariable_norm(ref_receiver))
        diff_statevariable_norm = integrate_in_time(time, statevariable_norm(difference))
        errors.loc[i, "state variable"] = diff_statevariable_norm / ref_statevariable_norm

    if "pressure" in quantities:   
        ref_pressure_norm = integrate_in_time(time, pressure_norm(ref_receiver))
        diff_pressure_norm = integrate_in_time(time, pressure_norm(difference))
        errors.loc[i, "pressure"] = diff_pressure_norm / ref_pressure_norm

    if "temperature" in quantities:   
        ref_temperature_norm = integrate_in_time(time, temperature_norm(ref_receiver))
        diff_temperature_norm = integrate_in_time(time, temperature_norm(difference))
        errors.loc[i, "temperature"] = diff_temperature_norm / ref_temperature_norm

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

sim_receiver_ids = find_all_receivers(args.output, False)
ref_receiver_ids = find_all_receivers(args.output_ref, False)
receiver_ids = np.intersect1d(sim_receiver_ids, ref_receiver_ids)    


receiver_errors = pd.DataFrame(index=receiver_ids, columns=["velocity", "stress"])
for i in receiver_ids:
    receiver_errors.loc[i,:] = receiver_diff(args, i)
print(receiver_errors)

for q in receiver_errors.columns:
    broken_receivers = receiver_errors.index[receiver_errors[q] > args.epsilon].tolist()
    print(f"{q} exceeds relative error of {args.epsilon} at receveivers {broken_receivers}")

quantities = ["absolute slip", "friction coefficient", "peak sliprate", "traction", "rupture time", "sliprate", "slip", "normal velocity"]
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
    faultreceiver_errors.loc[i,:] = local_errors.loc[i,:]
print(faultreceiver_errors)

for q in faultreceiver_errors.columns:
    broken_faultreceivers = faultreceiver_errors.index[faultreceiver_errors[q] > args.epsilon].tolist()
    print(f"{q} exceeds relative error of {args.epsilon} at faultreceveivers {broken_faultreceivers}")

if (receiver_errors > args.epsilon).any().any() or (faultreceiver_errors > args.epsilon).any().any():
    sys.exit(1)

