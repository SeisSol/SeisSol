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
the volume-weighted means are compared instead; the elements are paired by where
they are, since their ``global-id`` need not agree between the two files (see
``FACES_PER_ELEMENT``).
"""

import argparse
import sys

import numpy as np
import seissolxdmf as sx
from validation_report import write_report_json

# Cell data that describes the run rather than the solution.
METAFIELDS = [
    "partition",
    "fault-tag",
    "global-id",
    "clustering",
    "locationFlag",
]

# The surface and fault outputs number a face as FACES_PER_ELEMENT * element + local side.
# The local side follows the order in which SeisSol has read the vertices of the element,
# which may differ between two runs (a relabelled mesh, a different SeisSol version), so for
# those outputs only the element part of the global-id has to agree.
FACES_PER_ELEMENT = 4

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


def extent(*geometries):
    """The largest extent of the given point sets taken together, along any axis."""
    points = np.concatenate(geometries)
    return np.max(np.max(points, axis=0) - np.min(points, axis=0))


def point_labels(points, tolerance):
    """An integer per point, the same for points that agree to within ``tolerance``.

    Per axis, the sorted coordinates are split wherever two neighbours are more than
    ``tolerance`` apart. Unlike snapping to a grid there is no rounding boundary two copies
    of the same coordinate could end up on either side of, however close they are.
    """
    labels = np.empty(points.shape, dtype=np.int64)
    for axis in range(points.shape[1]):
        order = np.argsort(points[:, axis], kind="stable")
        gaps = np.diff(points[order, axis]) > tolerance
        labels[order, axis] = np.concatenate(([0], np.cumsum(gaps)))
    order = np.lexsort(labels.T[::-1])
    distinct = np.any(np.diff(labels[order], axis=0) != 0, axis=1)
    label = np.empty(len(points), dtype=np.int64)
    label[order] = np.concatenate(([0], np.cumsum(distinct)))
    return label


def cell_keys(labels, connect):
    """A hashable fingerprint per cell: the labels of its vertices, sorted.

    Sorting makes the key independent of the vertex order within the cell, which
    differs between outputs whose subcells are oriented differently.
    """
    return [tuple(cell) for cell in np.sort(labels[connect], axis=1).tolist()]


def match_cells(
    geom,
    connect,
    geom_ref,
    connect_ref,
    tags=None,
    tags_ref=None,
    owners=None,
    owners_ref=None,
    scale=None,
    what="cells",
):
    """Permutations ``ids``, ``ids_ref`` bringing the two cell lists into the same order.

    Returns ``None`` if the two files do not consist of the same cells. Vertices of both
    files that agree to within ``scale * 1e-9`` are given a common label (see
    ``point_labels``), so coordinates differing only by round-off match however close to
    each other they are; the matched cells are then checked against each other with a real
    tolerance. ``scale`` defaults to the extent of both files, ``what`` names the cells in
    the messages.

    Some outputs hold the same geometry more than once: the free-surface output writes a
    face of an elastic-acoustic interface once for either side. ``tags`` (the
    ``locationFlag``) become part of the key, so the two sides are never swapped, and
    ``owners`` (the ``global-id``) decide between cells that are alike even then.
    """
    if len(connect) != len(connect_ref):
        print(
            f"The two files have {len(connect)} and {len(connect_ref)} {what}; "
            f"they cannot consist of the same {what}."
        )
        return None

    if scale is None:
        scale = extent(geom, geom_ref)
    if scale <= 0:
        return None
    # coarse enough to absorb round-off, fine enough to separate distinct vertices
    labels = point_labels(np.concatenate((geom, geom_ref)), scale * 1e-9)

    keys = cell_keys(labels[: len(geom)], connect)
    keys_ref = cell_keys(labels[len(geom) :], connect_ref)
    if tags is not None and tags_ref is not None:
        keys = [key + (int(tag),) for key, tag in zip(keys, tags)]
        keys_ref = [key + (int(tag),) for key, tag in zip(keys_ref, tags_ref)]

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
                print(f"  {what[:-1]} {index} has no counterpart in the reference")
            continue
        choice = len(candidates) - 1
        if owners is not None and owners_ref is not None and len(candidates) > 1:
            alike = [
                position
                for position, candidate in enumerate(candidates)
                if owners_ref[candidate] == owners[index]
            ]
            if alike:
                choice = alike[-1]
        ids_ref[index] = candidates.pop(choice)

    if unmatched > 0:
        print(f"{unmatched} of {len(keys)} {what} could not be matched geometrically.")
        return None

    residual = np.max(
        np.abs(
            np.sort(geom[connect[ids]], axis=1)
            - np.sort(geom_ref[connect_ref[ids_ref]], axis=1)
        )
    )
    if residual > scale * 1e-8:
        print(
            f"Matched {what} differ by up to {residual:.3e}; the match is not trustworthy."
        )
        return None

    reordered = int(np.count_nonzero(ids_ref != ids))
    print(
        f"Matched all {len(keys)} {what} geometrically "
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
    subdivision sampled in a different vertex labeling looks like. One and the same
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


def global_ids_agree(ids, ids_ref, faces):
    """Whether the global-ids of matched cells, or of paired elements, say the same.

    For a volume output they have to be identical. For a face output (``faces``) only the
    element part has to be (see FACES_PER_ELEMENT), and the cells have to fall into the
    same faces in both files.
    """
    if np.array_equal(ids, ids_ref):
        return True
    if not faces:
        return False
    pairs = np.unique(np.column_stack((ids, ids_ref)), axis=0)
    same_faces = len(pairs) == len(np.unique(ids)) == len(np.unique(ids_ref))
    same_elements = np.array_equal(
        ids // FACES_PER_ELEMENT, ids_ref // FACES_PER_ELEMENT
    )
    return bool(same_faces and same_elements)


def list_quantities(file):
    """Physical data-field names in a mesh output (bookkeeping fields removed).

    Matches exactly the keys ``compare()`` reports, so a generated data file's
    per-quantity entries line up with what a real comparison produces.
    """
    mesh = sx.seissolxdmf(file)
    return [q for q in sorted(mesh.ReadAvailableDataFields()) if q not in METAFIELDS]


def compare(file, file_ref, epsilon, report_json=None, category="mesh"):
    mesh = sx.seissolxdmf(file)
    mesh_ref = sx.seissolxdmf(file_ref)

    geom = mesh.ReadGeometry()
    connect = mesh.ReadConnect()
    geom_ref = mesh_ref.ReadGeometry()
    connect_ref = mesh_ref.ReadConnect()
    for name, points in (("", geom), ("reference ", geom_ref)):
        if not np.all(np.isfinite(points)):
            print(f"The {name}file has vertex coordinates that are not finite.")
            sys.exit(1)

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

    tags = (
        mesh.Read1dData("locationFlag", mesh.nElements, isInt=True)
        if "locationFlag" in fields and "locationFlag" in fields_ref
        else None
    )
    tags_ref = (
        mesh_ref.Read1dData("locationFlag", mesh_ref.nElements, isInt=True)
        if tags is not None
        else None
    )
    # a flag per cell, or the file could not be read as it says; matching on a shorter list would
    # silently leave most of the cells out of the comparison
    for name, flags, cells in (
        ("", tags, connect),
        ("reference ", tags_ref, connect_ref),
    ):
        if flags is not None and len(flags) != len(cells):
            print(
                f"The {name}file holds {len(flags)} values of locationFlag for {len(cells)} "
                "cells; it cannot be read as described."
            )
            sys.exit(1)

    scale = extent(geom, geom_ref)
    matched = match_cells(
        geom,
        connect,
        geom_ref,
        connect_ref,
        tags=tags,
        tags_ref=tags_ref,
        owners=ids_global,
        owners_ref=ids_global_ref,
        scale=scale,
    )
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

    faces = connect.shape[1] == 3
    global_id_correct = None
    if not aggregated and ids_global is not None and ids_global_ref is not None:
        global_id_correct = global_ids_agree(
            ids_global[ids], ids_global_ref[ids_ref], faces
        )
        print(f"Global IDs present; conformant: {global_id_correct}")
        if not np.array_equal(ids_global[ids], ids_global_ref[ids_ref]):
            if global_id_correct:
                print("The global IDs differ in the local side of the faces only.")
            report_permutation(ids_ref, ids_global_ref)

    if aggregated:
        # global-id groups the cells of one file into elements, but its values need not agree
        # between the two files (see FACES_PER_ELEMENT); the elements are paired by their
        # centroid instead
        barycenters = np.mean(geom[connect], axis=1)
        barycenters_ref = np.mean(geom_ref[connect_ref], axis=1)

        parents, _, weights = aggregate_by_parent(
            ids_global, volumes, np.zeros(len(volumes))
        )
        parents_ref, _, weights_ref = aggregate_by_parent(
            ids_global_ref, volumes_ref, np.zeros(len(volumes_ref))
        )
        centroids = np.column_stack(
            [
                aggregate_by_parent(ids_global, volumes, barycenters[:, axis])[1]
                for axis in range(barycenters.shape[1])
            ]
        )
        centroids_ref = np.column_stack(
            [
                aggregate_by_parent(
                    ids_global_ref, volumes_ref, barycenters_ref[:, axis]
                )[1]
                for axis in range(barycenters_ref.shape[1])
            ]
        )
        parent_tags = parent_tags_ref = None
        if tags is not None:
            parent_tags = tags[np.unique(ids_global, return_index=True)[1]]
            parent_tags_ref = tags_ref[np.unique(ids_global_ref, return_index=True)[1]]
        paired = match_cells(
            centroids,
            np.arange(len(parents))[:, None],
            centroids_ref,
            np.arange(len(parents_ref))[:, None],
            tags=parent_tags,
            tags_ref=parent_tags_ref,
            owners=parents,
            owners_ref=parents_ref,
            scale=scale,
            what="elements",
        )
        if paired is None:
            print("The elements of the two files are not in the same place.")
            sys.exit(1)
        order, order_ref = paired
        if np.max(np.abs(weights[order] - weights_ref[order_ref])) > 1e-8 * np.max(
            weights
        ):
            print("The subdivisions of an element do not cover the same volume.")
            sys.exit(1)
        print(f"Aggregated {len(volumes)} cells into {len(parents)} elements.")
        global_id_correct = global_ids_agree(
            parents[order], parents_ref[order_ref], faces
        )
        print(f"Global IDs of the paired elements; conformant: {global_id_correct}")

    def l2_error(quantity, quantity_ref):
        if aggregated:
            _, mean, weight = aggregate_by_parent(ids_global, volumes, quantity)
            _, mean_ref, _ = aggregate_by_parent(
                ids_global_ref, volumes_ref, quantity_ref
            )
            mean, weight, mean_ref = mean[order], weight[order], mean_ref[order_ref]
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

    if report_json is not None:
        quantities = {q: float(errors[i]) for i, q in enumerate(quantity_names)}
        checks = None if global_id_correct is None else {"global-id": global_id_correct}
        write_report_json(
            report_json, category, epsilon, not failure, quantities, checks=checks
        )

    if failure:
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare two SeisSol XDMF outputs.")
    parser.add_argument("mesh", type=str)
    parser.add_argument("mesh_ref", type=str)
    parser.add_argument("--epsilon", type=float, default=0.01)
    args = parser.parse_args()

    compare(args.mesh, args.mesh_ref, args.epsilon)
