# SPDX-FileCopyrightText: 2022 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Compare two SeisSol XDMF outputs cell by cell.

The two files may enumerate their cells in any order, and may enumerate the
vertices within a cell in any order: cells are identified by their geometry, not
by their position in the file. That matters because the cell order of a refined
output depends on how the refinement is composed, and the vertex order within a
cell depends on how its subcells are oriented -- neither of which says anything
about whether the data agree.

If the two files do not consist of the same cells at all (a genuinely different
subdivision of the same elements), the cells are aggregated per ``global-id`` and
the volume-weighted means are compared instead.
"""

import argparse
import sys

import numpy as np
import seissolxdmf as sx

# Cell data that describes the run rather than the solution.
METAFIELDS = [
    "partition",
    "fault-tag",
    "global-id",
    "clustering",
    "locationFlag",
]

# Quantities that were renamed when the output modules were unified.
RENAMED = {"v1": "u", "v2": "v", "v3": "w", "pprime": "-p"}


def cell_volumes(geom, connect):
    """Length / area / volume of each cell, depending on its vertex count."""
    cells = geom[connect, :]
    if cells.shape[1] == 1:
        return np.zeros(cells.shape[0])
    if cells.shape[1] == 2:
        return np.linalg.norm(cells[:, 1, :] - cells[:, 0, :], axis=1)
    if cells.shape[1] == 3:
        a = cells[:, 1, :] - cells[:, 0, :]
        b = cells[:, 2, :] - cells[:, 0, :]
        return 0.5 * np.linalg.norm(np.cross(a, b), axis=1)
    if cells.shape[1] == 4:
        a = cells[:, 1, :] - cells[:, 0, :]
        b = cells[:, 2, :] - cells[:, 0, :]
        c = cells[:, 3, :] - cells[:, 0, :]
        return np.abs(np.einsum("nI,nI->n", np.cross(a, b), c)) / 6.0
    raise NotImplementedError(f"cells with {cells.shape[1]} vertices")


def cell_keys(geom, connect, quantum):
    """A hashable fingerprint per cell: its vertices, snapped to a grid and sorted.

    Sorting makes the key independent of the vertex order within the cell, which
    differs between outputs whose subcells are oriented differently.
    """
    snapped = np.rint(geom[connect, :] / quantum).astype(np.int64)
    order = np.lexsort((snapped[:, :, 2], snapped[:, :, 1], snapped[:, :, 0]), axis=1)
    snapped = np.take_along_axis(snapped, order[:, :, None], axis=1)
    return [tuple(cell.ravel()) for cell in snapped]


def match_cells(geom, connect, geom_ref, connect_ref):
    """Permutations ``ids``, ``ids_ref`` bringing the two cell lists into the same order.

    Returns ``None`` if the two files do not consist of the same cells. Vertices are
    snapped to a grid before hashing, so coordinates differing only by round-off still
    match; any ambiguity introduced by the snapping is ruled out afterwards by checking
    the matched cells against each other with a real tolerance.
    """
    if len(connect) != len(connect_ref):
        print(
            f"The two files have {len(connect)} and {len(connect_ref)} cells; "
            "they cannot consist of the same cells."
        )
        return None

    scale = np.max(np.max(geom, axis=0) - np.min(geom, axis=0))
    if scale <= 0:
        return None
    # coarse enough to absorb round-off, fine enough to separate distinct vertices
    quantum = scale * 1e-9

    keys = cell_keys(geom, connect, quantum)
    keys_ref = cell_keys(geom_ref, connect_ref, quantum)

    lookup = {}
    for index, key in enumerate(keys_ref):
        lookup.setdefault(key, []).append(index)

    ids = np.arange(len(keys))
    ids_ref = np.zeros(len(keys), dtype=np.int64)
    unmatched = 0
    for index, key in enumerate(keys):
        candidates = lookup.get(key)
        if not candidates:
            unmatched += 1
            if unmatched <= 3:
                print(f"  cell {index} has no counterpart in the reference")
            continue
        ids_ref[index] = candidates.pop()

    if unmatched > 0:
        print(f"{unmatched} of {len(keys)} cells could not be matched geometrically.")
        return None

    residual = np.max(
        np.abs(
            np.sort(geom[connect[ids]], axis=1)
            - np.sort(geom_ref[connect_ref[ids_ref]], axis=1)
        )
    )
    if residual > scale * 1e-8:
        print(
            f"Matched cells differ by up to {residual:.3e}; the match is not trustworthy."
        )
        return None

    reordered = int(np.count_nonzero(ids_ref != ids))
    print(
        f"Matched all {len(keys)} cells geometrically "
        f"(residual {residual:.3e}, {reordered} of them reordered)."
    )
    return ids, ids_ref


def report_permutation(ids_ref, group_ref, limit=8):
    """Print how the matched cells are permuted inside each parent element.

    One and the same permutation in every parent is the signature of a systematic
    mismatch -- a differently enumerated or differently sampled subdivision -- rather
    than of noise, and is worth seeing explicitly.
    """
    if group_ref is None or len(ids_ref) == 0:
        return
    grouped = group_ref[ids_ref]
    _, first, counts = np.unique(grouped, return_index=True, return_counts=True)
    per_parent = int(counts[0])
    if per_parent <= 1 or not np.all(counts == per_parent):
        return
    # the blocks have to be contiguous for the pattern below to mean anything
    if not all(
        np.all(grouped[start : start + per_parent] == grouped[start]) for start in first
    ):
        return

    patterns = {}
    for start in first:
        block = ids_ref[start : start + per_parent]
        pattern = tuple(int(i) for i in np.argsort(np.argsort(block)))
        patterns[pattern] = patterns.get(pattern, 0) + 1
    if len(patterns) == 1 and next(iter(patterns)) == tuple(range(per_parent)):
        return

    print(f"Subcell order within a parent element ({per_parent} subcells each):")
    for pattern, count in sorted(patterns.items(), key=lambda item: -item[1])[:limit]:
        print(f"  {list(pattern)}: {count} elements")
    if len(patterns) > limit:
        print(f"  ... and {len(patterns) - limit} further patterns")


def diagnose_subcell_permutation(groups, quantity, quantity_ref, sample=200, limit=6):
    """Check whether the values of a quantity are permuted *within* each parent element.

    Cells that match geometrically can still carry each other's values -- that is what a
    subdivision sampled in a different vertex labelling looks like. One and the same
    permutation across many elements is a strong hint at such a systematic mismatch, as
    opposed to a genuine numerical difference.
    """
    _, first, counts = np.unique(groups, return_index=True, return_counts=True)
    if counts.size == 0:
        return
    per_parent = int(counts[0])
    if per_parent <= 1 or not np.all(counts == per_parent):
        return
    if not all(
        np.all(groups[start : start + per_parent] == groups[start]) for start in first
    ):
        return

    patterns = {}
    considered = 0
    for start in first[:: max(1, len(first) // sample)]:
        block = quantity[start : start + per_parent]
        block_ref = quantity_ref[start : start + per_parent]
        spread = np.max(block_ref) - np.min(block_ref)
        if spread <= 0:
            continue
        # assign every value to the closest reference value inside the same element
        pattern = tuple(
            int(i)
            for i in np.argmin(np.abs(block[:, None] - block_ref[None, :]), axis=1)
        )
        if sorted(pattern) != list(range(per_parent)):
            continue
        residual = np.max(np.abs(block - block_ref[list(pattern)]))
        if residual > 1e-6 * spread:
            continue
        patterns[pattern] = patterns.get(pattern, 0) + 1
        considered += 1

    identity = tuple(range(per_parent))
    patterns.pop(identity, None)
    if considered == 0 or not patterns:
        return

    ranked = sorted(patterns.items(), key=lambda item: -item[1])
    share = 100.0 * ranked[0][1] / considered
    if share < 50.0:
        return
    print(
        f"  the values look permuted within each element: {list(ranked[0][0])} "
        f"in {share:.0f}% of the {considered} elements sampled "
        f"({per_parent} subcells each)"
    )
    for pattern, count in ranked[1:limit]:
        print(f"    {list(pattern)}: {100.0 * count / considered:.0f}%")


def aggregate_by_parent(groups, volumes, quantity):
    """Volume-weighted mean of ``quantity`` per parent element, plus its total volume."""
    order = np.argsort(groups, kind="stable")
    unique, first = np.unique(groups[order], return_index=True)
    weight = np.add.reduceat(volumes[order], first)
    weighted = np.add.reduceat((volumes * quantity)[order], first)
    return unique, weighted / np.where(weight > 0, weight, 1.0), weight


def compare(file, file_ref, epsilon):
    mesh = sx.seissolxdmf(file)
    mesh_ref = sx.seissolxdmf(file_ref)

    geom = mesh.ReadGeometry()
    connect = mesh.ReadConnect()
    geom_ref = mesh_ref.ReadGeometry()
    connect_ref = mesh_ref.ReadConnect()

    fields = mesh.ReadAvailableDataFields()
    fields_ref = mesh_ref.ReadAvailableDataFields()

    ids_global = (
        mesh.Read1dData("global-id", mesh.nElements, isInt=True)
        if "global-id" in fields
        else None
    )
    ids_global_ref = (
        mesh_ref.Read1dData("global-id", mesh_ref.nElements, isInt=True)
        if "global-id" in fields_ref
        else None
    )

    matched = match_cells(geom, connect, geom_ref, connect_ref)
    aggregated = matched is None

    if aggregated:
        if ids_global is None or ids_global_ref is None:
            print(
                "The two files do not consist of the same cells and carry no global-id, "
                "so there is nothing left to compare them by."
            )
            sys.exit(1)
        print(
            "Falling back to a per-element comparison: cells are aggregated by global-id "
            "and their volume-weighted means are compared."
        )
        ids = np.arange(len(connect))
        ids_ref = np.arange(len(connect_ref))
    else:
        ids, ids_ref = matched

    volumes = cell_volumes(geom, connect)[ids]
    volumes_ref = cell_volumes(geom_ref, connect_ref)[ids_ref]

    global_id_correct = None
    if not aggregated and ids_global is not None and ids_global_ref is not None:
        global_id_correct = bool(np.all(ids_global[ids] == ids_global_ref[ids_ref]))
        print(f"Global IDs present; conformant: {global_id_correct}")
        if not global_id_correct:
            report_permutation(ids_ref, ids_global_ref)

    if aggregated:
        scale = np.max(np.max(geom, axis=0) - np.min(geom, axis=0))
        barycenters = np.mean(geom[connect], axis=1)
        barycenters_ref = np.mean(geom_ref[connect_ref], axis=1)

        parents, _, weights = aggregate_by_parent(
            ids_global, volumes, np.zeros(len(volumes))
        )
        parents_ref, _, weights_ref = aggregate_by_parent(
            ids_global_ref, volumes_ref, np.zeros(len(volumes_ref))
        )
        if not np.array_equal(parents, parents_ref):
            print("The two files do not even cover the same elements.")
            sys.exit(1)
        if np.max(np.abs(weights - weights_ref)) > 1e-8 * np.max(weights):
            print("The subdivisions of an element do not cover the same volume.")
            sys.exit(1)
        for axis in range(barycenters.shape[1]):
            _, centroid, _ = aggregate_by_parent(
                ids_global, volumes, barycenters[:, axis]
            )
            _, centroid_ref, _ = aggregate_by_parent(
                ids_global_ref, volumes_ref, barycenters_ref[:, axis]
            )
            if np.max(np.abs(centroid - centroid_ref)) > 1e-8 * scale:
                print("The elements of the two files are not in the same place.")
                sys.exit(1)
        print(f"Aggregated {len(volumes)} cells into {len(parents)} elements.")

    def l2_error(quantity, quantity_ref):
        if aggregated:
            _, mean, weight = aggregate_by_parent(ids_global, volumes, quantity)
            _, mean_ref, _ = aggregate_by_parent(
                ids_global_ref, volumes_ref, quantity_ref
            )
            return (
                np.dot(weight, np.power(mean - mean_ref, 2)),
                np.dot(weight, np.power(mean_ref, 2)),
            )
        return (
            np.dot(volumes, np.power(quantity - quantity_ref, 2)),
            np.dot(volumes, np.power(quantity_ref, 2)),
        )

    quantity_names = sorted(name for name in fields if name not in METAFIELDS)
    errors = np.zeros(len(quantity_names))

    last_index = mesh.ndt
    assert last_index == mesh_ref.ndt
    for i, q in enumerate(quantity_names):
        quantity = mesh.ReadData(q, last_index - 1)[ids]
        q_ref = q if q in fields_ref else RENAMED[q]
        quantity_ref = mesh_ref.ReadData(q_ref, last_index - 1)[ids_ref]

        # we can leave this one in. A field with the name "DS" only appears on the mesh
        if q == "DS":
            print(
                "There is a bug on the master branch, which sets DS output to zero in wrong "
                "places. In order to make a fair comparison, we only compare the parts of DS, "
                "where it is non-zero."
            )
            quantity = np.where(quantity_ref < 1e-10, 0.0, quantity)

        difference, reference = l2_error(quantity, quantity_ref)
        if reference < 1e-10:
            print(f"{q:3}: {difference} [abs.]")
            errors[i] = difference
        else:
            print(f"{q:3}: {difference / reference} [rel.]")
            errors[i] = difference / reference

        if errors[i] > epsilon and not aggregated and ids_global is not None:
            diagnose_subcell_permutation(ids_global[ids], quantity, quantity_ref)

    failure = False
    if global_id_correct is False:
        print("Global IDs present, but did not match.")
        failure = True

    if np.any(errors > epsilon):
        print(f"Relative/absolute error {epsilon} exceeded for quantities")
        print([quantity_names[i] for i in np.where(errors > epsilon)[0]])
        failure = True

    if failure:
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare two SeisSol XDMF outputs.")
    parser.add_argument("mesh", type=str)
    parser.add_argument("mesh_ref", type=str)
    parser.add_argument("--epsilon", type=float, default=0.01)
    args = parser.parse_args()

    compare(args.mesh, args.mesh_ref, args.epsilon)
