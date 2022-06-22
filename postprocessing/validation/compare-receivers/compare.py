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
args = parser.parse_args()

def velocity_norm(receiver):
    return np.sqrt(receiver[7]**2 + receiver[8]**2 + receiver[9]**2)

def stress_norm(receiver):
    return np.sqrt(receiver[1]**2 + receiver[2]**2 + receiver[3]**2
            + receiver[4]**2 + receiver[5]**2 + receiver[6]**2)

def sliprate_norm(receiver):
    return np.sqrt(receiver["SRs"]**2 + receiver["SRd"]**2)

def traction_norm(receiver):
    return np.sqrt(receiver["Td0"]**2 + receiver["Ts0"]**2 + receiver["Pn0"]**2)

def normal_velocity_norm(receiver):
    return np.sqrt(receiver["u_n"]**2)

def friction_coefficient_norm(receiver):
    return np.sqrt(receiver["Mud"]**2)

def statevariable_norm(receiver):
    return np.sqrt(receiver["StV"]**2)

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
    #initial_stress= read_initial_stress(filename)
    #receiver[3] += initial_stress[1]
    #receiver[4] += initial_stress[2]
    #receiver[5] += initial_stress[0]
    return receiver

def faultreceiver_diff(args, i):
    sim_filename = f"{args.output}/tpv-faultreceiver-0000{i-1}-00000.dat"
    ref_filename = f"{args.output_ref}/tpv-faultreceiver-0000{i}-00000.dat"
    sim_receiver = read_faultreceiver(sim_filename)
    ref_receiver = read_faultreceiver(ref_filename).iloc[:-1]

    # both receivers must have the same time axis
    assert(np.max(np.abs(sim_receiver["Time"] - ref_receiver["Time"])) < 1e-14)
    time = sim_receiver["Time"]
    difference = sim_receiver - ref_receiver

    ref_sliprate_norm = integrate_in_time(time, sliprate_norm(ref_receiver))
    diff_sliprate_norm = integrate_in_time(time, sliprate_norm(difference))

    ref_traction_norm = integrate_in_time(time, traction_norm(ref_receiver))
    diff_traction_norm = integrate_in_time(time, traction_norm(difference))

    ref_normal_velocity_norm = integrate_in_time(time, normal_velocity_norm(ref_receiver))
    diff_normal_velocity_norm = integrate_in_time(time, normal_velocity_norm(difference))

    ref_friction_coefficient_norm = integrate_in_time(time, friction_coefficient_norm(ref_receiver))
    diff_friction_coefficient_norm = integrate_in_time(time, friction_coefficient_norm(difference))

    ref_statevariable_norm = integrate_in_time(time, statevariable_norm(ref_receiver))
    diff_statevariable_norm = integrate_in_time(time, statevariable_norm(difference))

    return diff_sliprate_norm / ref_sliprate_norm, diff_traction_norm / ref_traction_norm, diff_normal_velocity_norm / ref_normal_velocity_norm, diff_friction_coefficient_norm / ref_friction_coefficient_norm, diff_statevariable_norm / ref_statevariable_norm

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

sim_faultreceiver_ids = find_all_receivers(args.output, True) + 1
ref_faultreceiver_ids = find_all_receivers(args.output_ref, True)
faultreceiver_ids = np.intersect1d(sim_faultreceiver_ids, ref_faultreceiver_ids)    

sim_receiver_ids = find_all_receivers(args.output, False)
ref_receiver_ids = find_all_receivers(args.output_ref, False)
receiver_ids = np.intersect1d(sim_receiver_ids, ref_receiver_ids)    


receiver_errors = pd.DataFrame(index=receiver_ids, columns=["velocity", "stress"])
for i in receiver_ids:
    receiver_errors.iloc[i-1,:] = receiver_diff(args, i)
print(receiver_errors)

for q in receiver_errors.columns:
    broken_receivers = receiver_errors.index[receiver_errors[q] > args.epsilon].tolist()
    print(f"{q} exceeds relative error of {args.epsilon} at receveivers {broken_receivers}")

faultreceiver_errors = pd.DataFrame(index=faultreceiver_ids, columns=["sliprate", "traction", "normal_velocity", "friction_coefficient", "statevariable"])
for i in faultreceiver_ids:
    faultreceiver_errors.iloc[i-1,:] = faultreceiver_diff(args, i)
print(faultreceiver_errors)

for q in faultreceiver_errors.columns:
    broken_faultreceivers = faultreceiver_errors.index[faultreceiver_errors[q] > args.epsilon].tolist()
    print(f"{q} exceeds relative error of {args.epsilon} at faultreceveivers {broken_faultreceivers}")

if (receiver_errors > args.epsilon).any().any() or (faultreceiver_errors > args.epsilon).any().any():
    sys.exit(1)

