# SPDX-FileCopyrightText: 2022 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

import argparse
import sys

import numpy as np
import seissolxdmf as sx


def compare(file, file_ref, epsilon):
    mesh = sx.seissolxdmf(file)
    mesh_ref = sx.seissolxdmf(file_ref)

    # read mesh
    geom = mesh.ReadGeometry()
    connect = mesh.ReadConnect()
    geom_ref = mesh_ref.ReadGeometry()
    connect_ref = mesh_ref.ReadConnect()

    # read reference mesh
    if "global-id" in mesh_ref.ReadAvailableDataFields():
        # rely on stably-sorting the refined cells here
        preIds = mesh.Read1dData("global-id", mesh.nElements, isInt=True)
        ids = np.argsort(preIds)
        preIds_ref = mesh_ref.Read1dData("global-id", mesh.nElements, isInt=True)
        ids_ref = np.argsort(preIds)
    else:
        # if the reference solution has no global IDs, heuristically sort by barycenter, weighed by 3D cube position
        print(
            "No global-id field found in the reference. Order the cells by their barycenter in a 3D array."
        )
        rngs = [(np.max(geom[:, i]), np.min(geom[:, i])) for i in range(3)]
        offset = np.array([rng[1] for rng in rngs])
        scale = np.array(
            [
                1,
                rngs[0][1] - rngs[0][0],
                (rngs[0][1] - rngs[0][0]) * (rngs[1][1] - rngs[1][0]),
            ]
        )

        fullcoords = np.einsum(
            "s,Is->I", scale, np.mean(geom[connect], axis=1) - offset
        )
        ids = np.argsort(fullcoords, kind="stable")
        fullcoords_ref = np.einsum(
            "s,Is->I", scale, np.mean(geom_ref[connect_ref], axis=1) - offset
        )
        ids_ref = np.argsort(fullcoords_ref, kind="stable")

    connect = connect[ids]
    connect_ref = connect_ref[ids_ref]

    # assert both simulations were run on the same mesh
    assert np.all(np.abs(geom[connect] - geom_ref[connect_ref]) < 1e-10)

    def compute_integral(geom, connect, q):
        cells = geom[connect, :]
        if cells.shape[1] == 1:
            areas = 0
        elif cells.shape[1] == 2:
            areas = cells[:, 1, :] - cells[:, 0, :]
        elif cells.shape[1] == 3:
            a = cells[:, 1, :] - cells[:, 0, :]
            b = cells[:, 2, :] - cells[:, 0, :]
            areas = 0.5 * np.linalg.norm(np.cross(a, b), axis=1)
        elif cells.shape[1] == 4:
            a = cells[:, 1, :] - cells[:, 0, :]
            b = cells[:, 2, :] - cells[:, 0, :]
            c = cells[:, 3, :] - cells[:, 0, :]
            areas = np.abs(np.einsum("nI,nI->n", np.cross(a, b), c)) / 6.0
        else:
            raise NotImplementedError()
        return np.dot(areas, q)

    def l1_norm(q):
        return compute_integral(geom, connect, np.abs(q))

    def l1_difference(q_0, q_1):
        return l1_norm(q_0 - q_1)

    def l2_norm(q):
        return compute_integral(geom, connect, np.power(q, 2))

    def l2_difference(q_0, q_1):
        return l2_norm(q_0 - q_1)

    quantity_names = sorted(mesh.ReadAvailableDataFields())
    for metafield in [
        "partition",
        "fault-tag",
        "global-id",
        "clustering",
        "locationFlag",
    ]:
        if metafield in quantity_names:
            quantity_names.remove(metafield)
    errors = np.zeros((len(quantity_names)))

    last_index = mesh.ndt
    assert last_index == mesh_ref.ndt
    for i, q in enumerate(quantity_names):
        # extract quantity
        quantity = mesh.ReadData(q, last_index - 1)
        if q in mesh_ref.ReadAvailableDataFields():
            q_ref = q
        else:
            q_ref = {"v1": "u", "v2": "v", "v3": "w"}[q]
        quantity_ref = mesh_ref.ReadData(q_ref, last_index - 1)

        # re-assign data
        quantity = quantity[ids]
        quantity_ref = quantity_ref[ids_ref]

        # we can leave this one in. A field with the name "DS" only appears on the mesh
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
            print(f"{q:3}: {absolute_error} [abs.]")
            errors[i] = absolute_error
        else:
            relative_error = l2_difference(quantity, quantity_ref) / ref_norm
            print(f"{q:3}: {relative_error} [rel.]")
            errors[i] = relative_error

    if np.any(errors > epsilon):
        print(f"Relative/absolute error {epsilon} exceeded for quantities")
        print([quantity_names[i] for i in np.where(errors > epsilon)[0]])
        sys.exit(1)
