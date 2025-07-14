#!/usr/bin/env python3

import argparse
import sys

import numpy as np
import seissolxdmf as sx

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare two faults.")
    parser.add_argument("fault", type=str)
    parser.add_argument("fault_ref", type=str)
    parser.add_argument("--epsilon", type=float, default=0.01)

    args = parser.parse_args()
    fault = sx.seissolxdmf(args.fault)
    fault_ref = sx.seissolxdmf(args.fault_ref)

    # read mesh
    preIds = fault.Read1dData("global-id", fault.nElements, isInt=True)
    ids = np.argsort(preIds)
    geom = fault.ReadGeometry()
    connect = fault.ReadConnect()[ids]
    geom_ref = fault_ref.ReadGeometry()
    connect_ref = fault_ref.ReadConnect()

    # read reference mesh
    if "global-id" in fault_ref.ReadAvailableDataFields():
        # rely on stably-sorting the refined triangles here
        preIds_ref = fault_ref.Read1dData("global-id", fault.nElements, isInt=True)
        ids_ref = np.argsort(preIds)
    else:
        # if the reference solution has no global IDs, painstakingly compare cell by cell
        ids_ref = np.empty((fault_ref.nElements,), dtype=np.int64)
        fullcoords = geom[connect]
        for i in range(fault_ref.nElements):
            coords = geom_ref[connect_ref[i, :]]
            matchId = np.where(np.all(np.abs(fullcoords - coords) < 1e-10, axis=(1, 2)))
            assert len(matchId[0]) == 1
            ids_ref[matchId[0][0]] = i

    connect_ref = connect_ref[ids_ref]

    # assert both simulations were run on the same mesh
    assert np.all(np.abs(geom[connect] - geom_ref[connect_ref]) < 1e-10)

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
    if "global-id" in quantity_names:
        quantity_names.remove("global-id")
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

        # re-assign data
        quantity = quantity[ids]
        quantity_ref = quantity_ref[ids_ref]

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
