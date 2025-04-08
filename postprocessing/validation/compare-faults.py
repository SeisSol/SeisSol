#!/usr/bin/env python3

import argparse
import numpy as np
import seissolxdmf as sx
import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare two faults.")
    parser.add_argument("fault", type=str)
    parser.add_argument("fault_ref", type=str)
    parser.add_argument("--epsilon", type=float, default=0.01)

    args = parser.parse_args()
    fault = sx.seissolxdmf(args.fault)
    fault_ref = sx.seissolxdmf(args.fault_ref)

    # read mesh
    geom = fault.ReadGeometry()
    connect = fault.ReadConnect()
    # assert both simulations were run on the same mesh
    assert np.all(np.abs(geom - fault_ref.ReadGeometry()) < 1e-10)
    assert np.all(connect == fault_ref.ReadConnect())

    def compute_integral(geom, connect, q):
        triangles = geom[connect, :]
        a = triangles[:, 1, :] - triangles[:, 0, :]
        b = triangles[:, 2, :] - triangles[:, 0, :]
        areas = 0.5 * np.linalg.norm(np.cross(a, b), axis=1)

        return np.dot(areas, q)

    def l1_norm(q):
        return compute_integral(geom, connect, np.abs(q))

    def l1_difference(q_0, q_1):
        return l1_norm(q_0 - q_1)

    def l2_norm(q):
        return compute_integral(geom, connect, np.power(q, 2))

    def l2_difference(q_0, q_1):
        return l2_norm(q_0 - q_1)

    quantity_names = sorted(fault.ReadAvailableDataFields())
    quantity_names.remove("partition")
    if "fault-tag" in quantity_names:
        quantity_names.remove("fault-tag")
    errors = np.zeros((len(quantity_names)))

    last_index = fault.ndt
    assert last_index == fault_ref.ndt
    for i, q in enumerate(quantity_names):
        # extract quantity
        quantity = fault.ReadData(q, last_index - 1)
        quantity_ref = fault_ref.ReadData(q, last_index - 1)
        if q == "DS":
            print(
                "There is a bug on the master branch, which sets DS output to zero in wrong places. In order to make a fair comparison, we only compare the parts of DS, where it is non-zero."
            )
            ds_equals_zero = np.where(quantity_ref < 1e-10)
            quantity[ds_equals_zero] = 0.0

        # compute error
        ref_norm = l2_norm(quantity_ref)
        if ref_norm < 1e-10:
            absolute_error = l2_difference(quantity, quantity_ref)
            print(f"{q:3}: {absolute_error}")
            errors[i] = absolute_error
        else:
            relative_error = l2_difference(quantity, quantity_ref) / ref_norm
            print(f"{q:3}: {relative_error}")
            errors[i] = relative_error

    if np.any(errors > args.epsilon):
        print(f"Relative error {args.epsilon} exceeded for quantities")
        print([quantity_names[i] for i in np.where(errors > args.epsilon)[0]])
        sys.exit(1)
