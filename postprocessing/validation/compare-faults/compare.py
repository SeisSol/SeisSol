import argparse
import numpy as np
import seissolxdmf as sx
import sys

parser = argparse.ArgumentParser(description='Compare two faults.')
parser.add_argument('fault', type=str)
parser.add_argument('fault_ref', type=str)
parser.add_argument('--epsilon', type=float, default=0.01, required=False)

args = parser.parse_args()
fault = sx.seissolxdmf(args.fault)
fault_ref = sx.seissolxdmf(args.fault_ref)

# read mesh
geom = fault.ReadGeometry()
connect = fault.ReadConnect()
# assert both simulations were run on the same mesh
assert(np.all(np.abs(geom - fault_ref.ReadGeometry()) < 1e-10))
assert(np.all(connect == fault_ref.ReadConnect()))

def compute_integral(geom, connect, q):
    triangles = geom[connect,:]
    a = triangles[:,1,:] - triangles[:,0,:]
    b = triangles[:,2,:] - triangles[:,0,:]
    areas = 0.5*np.linalg.norm(np.cross(a, b), axis=1)

    return np.dot(areas, q)

def l1_norm(q):
    return compute_integral(geom, connect, np.abs(q))

def l1_difference(q_0, q_1):
    return l1_norm(q_0 - q_1)

def l2_norm(q):
    return compute_integral(geom, connect, np.power(q, 2))

def l2_difference(q_0, q_1):
    return l2_norm(q_0 - q_1)

quantity_names = ["ASl", "Mud", "PSR", "P_n", "Pn0", "RT", "SRd", "SRs", "Sld", "Sls", "T_d", "T_s", "Td0", "Ts0", "Vr", "u_n"]
errors = np.zeros((len(quantity_names)))

for i, q in enumerate(quantity_names):
    # extract quantity
    quantity = fault.ReadData(q, 2)
    quantity_ref = fault_ref.ReadData(q, 2)
    # compute error
    relative_error = l2_difference(quantity, quantity_ref) / l2_norm(quantity_ref)
    print(f"{q:3}: {relative_error}")
    errors[i] = relative_error

if np.any(errors > args.epsilon):
    print(f"Relative error {args.epsilon} exceeded for quantities")
    print([quantity_names[i] for i in np.where(errors > args.epsilon)[0]])
    sys.exit(1)

    
