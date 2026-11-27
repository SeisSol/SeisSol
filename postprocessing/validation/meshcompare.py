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

    matching = False
    global_id_correct = None

    # read reference mesh
    if "global-id" in mesh_ref.ReadAvailableDataFields():
        # rely on stably-sorting the refined cells here
        preIds = mesh.Read1dData("global-id", mesh.nElements, isInt=True)
        ids = np.argsort(preIds)
        preIds_ref = mesh_ref.Read1dData("global-id", mesh_ref.nElements, isInt=True)
        ids_ref = np.argsort(preIds_ref)

        global_id_correct = np.all(np.abs(geom[connect[ids]] - geom_ref[connect_ref[ids_ref]]) < 1e-10)
        matching = global_id_correct
        print(f'Global IDs present; conformant: {global_id_correct}')

    # positions closer than this count as the same, in the orders below: a vertex or a barycenter
    # can carry the rounding of how its run computed it (e.g. the subcells of a refined output)
    quantum = 1e-8 * np.max(np.max(geom, axis=0) - np.min(geom, axis=0))

    def cell_vertices(geom, connect):
        # the vertices of each cell in lexicographic order: a partition may list the vertices of a
        # cell (e.g. of a fault or surface triangle) in another order than the reference does
        cells = geom[connect]
        snapped = np.round(cells / quantum)
        order = np.lexsort((snapped[:, :, 2], snapped[:, :, 1], snapped[:, :, 0]), axis=-1)
        return np.take_along_axis(cells, order[:, :, None], axis=1)

    if not matching:
        # if the reference solution has no global IDs, sort the cells by their barycenter,
        # coordinate by coordinate: snapped to the grid of `quantum`, so that rounding in the
        # barycenters does not reorder them, and once more on a grid shifted by half a step, for a
        # barycenter that lies on a boundary of the first one
        print(
            "No or incorrect global-id field found in the reference. Order the cells by their barycenter."
        )
        def element_part(mesh):
            # Two cells may lie in the same place: a surface output writes a face of an interface
            # once for either side. They belong to different elements, and that tells them apart:
            # the surface and fault outputs number a face as 4 * element + local side, where only
            # the element is independent of how SeisSol has read the vertices.
            if "global-id" not in mesh.ReadAvailableDataFields():
                return np.zeros(mesh.nElements, dtype=np.int64)
            gids = mesh.Read1dData("global-id", mesh.nElements, isInt=True)
            return gids // 4 if mesh.ReadConnect().shape[1] == 3 else gids

        element = element_part(mesh)
        element_ref = element_part(mesh_ref)

        def by_barycenter(geom, connect, element, shift):
            snapped = np.floor(np.mean(geom[connect], axis=1) / quantum + shift)
            return np.lexsort((element, snapped[:, 2], snapped[:, 1], snapped[:, 0]))

        for shift in (0.0, 0.5):
            ids = by_barycenter(geom, connect, element, shift)
            ids_ref = by_barycenter(geom_ref, connect_ref, element_ref, shift)
            if np.all(
                np.abs(cell_vertices(geom, connect[ids]) - cell_vertices(geom_ref, connect_ref[ids_ref]))
                < 1e-10
            ):
                break

    connect = connect[ids]
    connect_ref = connect_ref[ids_ref]

    # assert both simulations were run on the same mesh
    assert np.all(np.abs(cell_vertices(geom, connect) - cell_vertices(geom_ref, connect_ref)) < 1e-10)

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

    failure = False
    if global_id_correct == False:
        # Not an error by itself: the cells were then paired by their barycenters, and the
        # assertion above checked that they agree. The order of a refined output's subcells
        # within a cell, and with it the order the global IDs give, depends on the partitioning.
        print('Global IDs present, but they order the cells differently; paired by barycenter.')

    if np.any(errors > epsilon):
        print(f"Relative/absolute error {epsilon} exceeded for quantities")
        print([quantity_names[i] for i in np.where(errors > epsilon)[0]])
        failure = True

    if failure:
        sys.exit(1)
