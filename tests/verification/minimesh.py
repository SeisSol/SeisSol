#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Generate and validate the tiny PUML meshes used by the verification pipeline.

The meshes here have between 6 and roughly 200 cells. They are not meant to
resolve anything; they exist so that every code path (material, boundary
condition, friction law, plasticity) can be exercised in seconds, and so that a
drift in the numerics is visible in a handful of numbers.

The meshes are generated rather than checked in: the specification is text, the
HDF5 file is a build artifact. A generated mesh can be parameterised (skew,
fault plane, periodicity, cell ordering), which a checked-in binary cannot.

Cell ordering matters: PUML's ``none`` partitioner keeps the contiguous
distribution implied by the file order, so the order chosen here fully
determines the MPI decomposition. That is what makes an MPI-versus-serial
comparison reproducible; a graph partitioner would decide differently for every
rank count and library version.
"""

import argparse
import hashlib
import json
import sys
import tempfile
from dataclasses import asdict, dataclass, replace
from itertools import permutations
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------
# PUML conventions
#
# The tables below mirror what SeisSol's reader expects; they are not free
# choices. A face slot in /boundary is a *PUML* face index, i.e. it is defined
# by the local vertex triple listed here. SeisSol maps that to its own side
# numbering internally.
# --------------------------------------------------------------------------

FACE_VERTICES = ((1, 0, 2), (0, 1, 3), (1, 2, 3), (2, 0, 3))

FACE_TYPES = {
    "regular": 0,
    "free-surface": 1,
    "free-surface-gravity": 2,
    "dynamic-rupture": 3,
    "dirichlet": 4,
    "outflow": 5,
    "analytical": 7,
}

# a face carrying one of these has to have a neighbour; every other tag must
# not have one. SeisSol warns (and then aborts) if that is violated.
INTERNAL_FACE_TYPES = frozenset({0, 3})

# name -> (dtype, bits per face; None for one entry per face)
BOUNDARY_FORMATS = {
    "i32": (np.int32, 8),
    "i64": (np.int64, 16),
    "i32x4": (np.int32, None),
}

AXES = ("x", "y", "z")
PLANES = ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")

MASK64 = (1 << 64) - 1


# --------------------------------------------------------------------------
# specification
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Spec:
    """A fully resolved mesh specification.

    ``cells`` counts cubes per axis; every cube is split into six tetrahedra,
    so the cell count is ``6 * nx * ny * nz``.
    """

    name: str
    cells: tuple = (2, 2, 2)
    size: tuple = (100.0, 100.0, 100.0)
    origin: tuple = None
    boundary: tuple = (("all", 1),)
    faults: tuple = ()
    periodic: tuple = (False, False, False)
    group_splits: tuple = ()
    group_base: int = 1
    skew: float = 0.0
    seed: int = 0
    order: str = "lexicographic"
    boundary_format: str = "i32"

    def resolved_origin(self):
        if self.origin is not None:
            return tuple(self.origin)
        return tuple(-0.5 * length for length in self.size)

    def boundary_tags(self):
        """Expand the boundary specification into one tag per outer plane."""
        tags = {}
        for key, tag in self.boundary:
            if key == "all":
                tags.update({plane: tag for plane in PLANES})
            else:
                tags[key] = tag
        for axis, periodic in enumerate(self.periodic):
            if periodic:
                # a periodic face is an interior face; it has a neighbour
                tags[f"{AXES[axis]}min"] = FACE_TYPES["regular"]
                tags[f"{AXES[axis]}max"] = FACE_TYPES["regular"]
        missing = [plane for plane in PLANES if plane not in tags]
        if missing:
            raise ValueError(f"no boundary condition given for {', '.join(missing)}")
        return tags

    def cell_count(self):
        nx, ny, nz = self.cells
        return 6 * nx * ny * nz


PRESETS = {
    # the smallest thing that is still a mesh: one cube. Every outer plane can
    # carry a different condition, which makes this the natural host for the
    # boundary-condition matrix and for the negative tests.
    "box6": Spec(name="box6", cells=(1, 1, 1)),
    # 2x2x2 cubes: has interior faces, so it exercises the neighbour flux and
    # can be decomposed over up to eight ranks.
    "box48": Spec(name="box48"),
    # same topology, interior vertices displaced. Axis-aligned cubes hide
    # errors in face normals, Jacobians and tensor rotation; this one does not.
    "skew48": Spec(name="skew48", skew=0.2),
    # a fault plane through the centre, giving eight dynamic-rupture faces.
    "fault48": Spec(name="fault48", faults=(("x", None),)),
    "faultskew48": Spec(name="faultskew48", faults=(("x", None),), skew=0.2),
    # two material groups separated by a plane, for material interfaces and
    # (later) multiple configurations.
    "mat48": Spec(name="mat48", group_splits=(("x", None),)),
    # periodic along the propagation direction only: a wave leaves on one side
    # and comes back on the other, so after one period the state must equal the
    # initial state. That is a sharp test which needs no reference data.
    "periodicx72": Spec(
        name="periodicx72", cells=(3, 2, 2), periodic=(True, False, False)
    ),
    # periodic in all directions, for a plane wave without any lateral
    # condition to satisfy. Three cells per axis is the minimum, hence the size.
    "periodic162": Spec(
        name="periodic162", cells=(3, 3, 3), periodic=(True, True, True)
    ),
    # elongated and open on every side: energy has to decay monotonically.
    "open96": Spec(
        name="open96",
        cells=(4, 2, 2),
        size=(200.0, 100.0, 100.0),
        boundary=(("all", 5),),
    ),
}


# --------------------------------------------------------------------------
# construction
# --------------------------------------------------------------------------


@dataclass
class Mesh:
    spec: Spec
    coords: np.ndarray  # (nv, 3) float64
    lattice: np.ndarray  # (nv, 3) int, the structured index of each vertex
    connect: np.ndarray  # (nc, 4) int
    tags: np.ndarray  # (nc, 4) int, in PUML face order
    groups: np.ndarray  # (nc,) int
    identify: np.ndarray = None  # (nv,) int, only for periodic meshes


def _jitter(index, axis, seed):
    """A deterministic value in [-1, 1) from an integer lattice position.

    Hashing rather than seeding a PRNG keeps the displacement of a vertex
    independent of iteration order and of the mesh it appears in, so the same
    vertex is displaced identically in a 2x2x2 and a 4x4x4 mesh.
    """
    value = (
        int(index[0]) * 73856093
        ^ int(index[1]) * 19349663
        ^ int(index[2]) * 83492791
        ^ (axis + 1) * 2654435761
        ^ seed * 40503
    ) & MASK64
    value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & MASK64
    value ^= value >> 31
    return (value / float(1 << 64)) * 2.0 - 1.0


def _preserved_planes(spec):
    """Grid indices per axis whose plane must stay flat.

    The outer planes stay flat so that a free surface really is a plane, and
    fault and group-split planes stay flat so that a half-space description in
    easi matches the mesh. Everything else may move.
    """
    preserved = [set() for _ in AXES]
    for axis, count in enumerate(spec.cells):
        preserved[axis].update({0, count})
    for axis, index in _normalise_planes(spec, spec.faults):
        preserved[axis].add(index)
    for axis, index in _normalise_planes(spec, spec.group_splits):
        preserved[axis].add(index)
    return preserved


def _normalise_planes(spec, planes):
    """Resolve ``(axis, index)`` entries, filling in the centre for ``None``."""
    resolved = []
    for axis_name, index in planes:
        axis = AXES.index(axis_name)
        count = spec.cells[axis]
        if index is None:
            if count % 2 != 0:
                raise ValueError(
                    f"cannot centre a plane along {axis_name}: {count} cells is odd"
                )
            index = count // 2
        if not 0 < index < count:
            raise ValueError(
                f"plane index {index} along {axis_name} is not interior "
                f"(expected 0 < index < {count})"
            )
        resolved.append((axis, index))
    return tuple(resolved)


def _vertices(spec):
    nx, ny, nz = spec.cells
    size = np.asarray(spec.size, dtype=np.float64)
    origin = np.asarray(spec.resolved_origin(), dtype=np.float64)
    counts = np.asarray(spec.cells, dtype=np.int64)
    spacing = size / counts

    lattice = np.array(
        [
            (i, j, k)
            for i in range(nx + 1)
            for j in range(ny + 1)
            for k in range(nz + 1)
        ],
        dtype=np.int64,
    )
    coords = origin + lattice * spacing

    if spec.skew != 0.0:
        preserved = _preserved_planes(spec)
        # a periodic pair of vertices has to be displaced identically, so the
        # jitter is taken from the reduced index
        reduced = lattice.copy()
        for axis, periodic in enumerate(spec.periodic):
            if periodic:
                reduced[:, axis] %= spec.cells[axis]
        for vertex in range(lattice.shape[0]):
            for axis in range(3):
                if lattice[vertex, axis] in preserved[axis]:
                    continue
                offset = _jitter(reduced[vertex], axis, spec.seed)
                coords[vertex, axis] += spec.skew * spacing[axis] * offset

    return coords, lattice


def _vertex_ids(spec):
    nx, ny, nz = spec.cells

    def vid(i, j, k):
        return (i * (ny + 1) + j) * (nz + 1) + k

    return vid


def _kuhn_cells(spec):
    """Split every cube into six tetrahedra (Kuhn/Freudenthal).

    Each tetrahedron follows one of the six monotone paths from the lowest to
    the highest cube corner. The resulting triangulation is conforming across
    cube faces because the diagonal on a shared face is determined by the
    global lattice order, which is translation invariant. It also puts a face
    on every cube face and nowhere else inside a coordinate plane, so a fault
    plane on the grid catches exactly the faces meant for it.
    """
    nx, ny, nz = spec.cells
    vid = _vertex_ids(spec)
    cells = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                base = np.array((i, j, k), dtype=np.int64)
                for path in permutations(range(3)):
                    corners = [base.copy()]
                    walk = base.copy()
                    for axis in path[:2]:
                        walk = walk.copy()
                        walk[axis] += 1
                        corners.append(walk)
                    corners.append(base + 1)
                    cells.append(((i, j, k), tuple(vid(*c) for c in corners)))
    return cells


def _morton(index):
    code = 0
    for bit in range(21):
        for axis in range(3):
            code |= ((int(index[axis]) >> bit) & 1) << (3 * bit + axis)
    return code


def _sort_cells(cells, order):
    if order == "lexicographic":
        key = lambda entry: entry[0]  # noqa: E731
    elif order == "morton":
        key = lambda entry: _morton(entry[0])  # noqa: E731
    elif order.startswith("slab-") and order[5:] in AXES:
        axis = AXES.index(order[5:])
        key = lambda entry: (  # noqa: E731
            entry[0][axis],
            entry[0],
        )
    else:
        raise ValueError(f"unknown cell order: {order}")
    return sorted(range(len(cells)), key=lambda n: (key(cells[n]), n))


def _fix_orientation(coords, connect):
    """Make every tetrahedron right-handed, as SeisSol requires."""
    edges = coords[connect[:, 1:]] - coords[connect[:, 0]][:, None, :]
    volume = np.linalg.det(edges)
    flipped = volume < 0.0
    connect[flipped, 2], connect[flipped, 3] = (
        connect[flipped, 3],
        connect[flipped, 2].copy(),
    )
    return connect


def _assign_tags(spec, connect, lattice):
    """Derive one boundary tag per face from the structured index.

    Working on integer lattice indices instead of coordinates keeps this exact:
    a face lies on a plane if and only if its three vertices share the index
    along that axis, no tolerance involved. Since the tags are derived after
    the orientation fix, a flipped tetrahedron cannot silently move a tag to
    the wrong slot.
    """
    plane_tags = spec.boundary_tags()
    faults = dict(_normalise_planes(spec, spec.faults))
    tags = np.zeros((connect.shape[0], 4), dtype=np.int64)

    for slot, local in enumerate(FACE_VERTICES):
        face = lattice[connect[:, local]]  # (nc, 3 vertices, 3 axes)
        for axis in range(3):
            column = face[:, :, axis]
            flat = (column[:, 0] == column[:, 1]) & (column[:, 1] == column[:, 2])
            plane = column[:, 0]
            on_min = flat & (plane == 0)
            on_max = flat & (plane == spec.cells[axis])
            tags[on_min, slot] = plane_tags[f"{AXES[axis]}min"]
            tags[on_max, slot] = plane_tags[f"{AXES[axis]}max"]
            if axis in faults:
                on_fault = flat & (plane == faults[axis])
                tags[on_fault, slot] = FACE_TYPES["dynamic-rupture"]
    return tags


def _assign_groups(spec, cubes):
    splits = _normalise_planes(spec, spec.group_splits)
    groups = np.full(len(cubes), spec.group_base, dtype=np.int64)
    for offset, (axis, index) in enumerate(splits):
        beyond = np.array([cube[axis] >= index for cube in cubes])
        groups += beyond.astype(np.int64) << offset
    return groups


def _assign_identify(spec, lattice):
    if not any(spec.periodic):
        return None
    dims = [
        spec.cells[axis] if periodic else spec.cells[axis] + 1
        for axis, periodic in enumerate(spec.periodic)
    ]
    reduced = lattice.copy()
    for axis, periodic in enumerate(spec.periodic):
        if periodic:
            reduced[:, axis] %= spec.cells[axis]
    return (reduced[:, 0] * dims[1] + reduced[:, 1]) * dims[2] + reduced[:, 2]


def build(spec):
    """Build a mesh from a specification."""
    coords, lattice = _vertices(spec)
    cells = _kuhn_cells(spec)
    permutation = _sort_cells(cells, spec.order)
    cubes = [cells[n][0] for n in permutation]
    connect = np.array([cells[n][1] for n in permutation], dtype=np.int64)
    connect = _fix_orientation(coords, connect)
    return Mesh(
        spec=spec,
        coords=coords,
        lattice=lattice,
        connect=connect,
        tags=_assign_tags(spec, connect, lattice),
        groups=_assign_groups(spec, cubes),
        identify=_assign_identify(spec, lattice),
    )


# --------------------------------------------------------------------------
# encoding and writing
# --------------------------------------------------------------------------


def encode_boundary(tags, boundary_format):
    dtype, bits = BOUNDARY_FORMATS[boundary_format]
    if bits is None:
        return tags.astype(dtype)
    limit = 1 << bits
    if tags.max(initial=0) >= limit:
        raise ValueError(f"boundary tag does not fit into {bits} bits")
    packed = np.zeros(tags.shape[0], dtype=np.int64)
    for slot in range(4):
        packed |= tags[:, slot] << (bits * slot)
    return packed.astype(dtype)


def decode_boundary(raw, boundary_format):
    """Undo :func:`encode_boundary`, the way SeisSol's reader does."""
    _, bits = BOUNDARY_FORMATS[boundary_format]
    if bits is None:
        return np.asarray(raw, dtype=np.int64)
    values = np.asarray(raw)
    # the reader interprets the packed word as unsigned
    unsigned = values.astype(np.uint64 if bits == 16 else np.uint32)
    mask = (1 << bits) - 1
    return np.stack(
        [(unsigned >> np.uint32(bits * slot)) & mask for slot in range(4)],
        axis=1,
    ).astype(np.int64)


def write(mesh, path):
    """Write the mesh as PUML/HDF5, matching what PUMgen produces."""
    import h5py

    spec = mesh.spec
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    ascii_str = h5py.string_dtype(encoding="ascii")

    with h5py.File(path, "w") as handle:
        handle.create_dataset("connect", data=mesh.connect.astype(np.uint64))
        handle.create_dataset("geometry", data=mesh.coords.astype(np.float64))
        handle.create_dataset("group", data=mesh.groups.astype(np.int32))
        handle.create_dataset(
            "boundary", data=encode_boundary(mesh.tags, spec.boundary_format)
        )
        handle.attrs.create("boundary-format", spec.boundary_format, dtype=ascii_str)
        if mesh.identify is None:
            handle.attrs.create("topology-format", "geometric", dtype=ascii_str)
        else:
            handle.create_dataset("identify", data=mesh.identify.astype(np.uint64))
            handle.attrs.create("topology-format", "identify-vertex", dtype=ascii_str)
    return path


def write_xdmf(mesh, mesh_path, path):
    """Write a minimal XDMF next to the mesh, for looking at it in ParaView."""
    cells, vertices = mesh.connect.shape[0], mesh.coords.shape[0]
    name = Path(mesh_path).name
    path = Path(path)
    path.write_text(
        '<?xml version="1.0" ?>\n'
        '<Xdmf Version="2.0">\n'
        " <Domain>\n"
        f'  <Grid Name="{mesh.spec.name}">\n'
        f'   <Topology TopologyType="Tetrahedron" NumberOfElements="{cells}">\n'
        f'    <DataItem NumberType="UInt" Precision="8" Format="HDF" Dimensions="{cells} 4">'
        f"{name}:/connect</DataItem>\n"
        "   </Topology>\n"
        '   <Geometry GeometryType="XYZ">\n'
        f'    <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="{vertices} 3">'
        f"{name}:/geometry</DataItem>\n"
        "   </Geometry>\n"
        '   <Attribute Name="group" Center="Cell">\n'
        f'    <DataItem NumberType="Int" Precision="4" Format="HDF" Dimensions="{cells}">'
        f"{name}:/group</DataItem>\n"
        "   </Attribute>\n"
        "  </Grid>\n"
        " </Domain>\n"
        "</Xdmf>\n",
        encoding="utf-8",
    )
    return path


def manifest(mesh, mesh_path):
    """Describe the mesh, including a hash of every dataset.

    A reference snapshot records this identifier. If the mesh changes, the
    comparison fails with a clear reason instead of looking like a drift in the
    numerics.
    """
    digests = {}
    for name, array in (
        ("connect", mesh.connect.astype(np.uint64)),
        ("geometry", mesh.coords.astype(np.float64)),
        ("group", mesh.groups.astype(np.int32)),
        ("boundary", encode_boundary(mesh.tags, mesh.spec.boundary_format)),
    ):
        digests[name] = hashlib.sha256(
            np.ascontiguousarray(array).tobytes()
        ).hexdigest()
    combined = hashlib.sha256(
        "".join(digests[key] for key in sorted(digests)).encode("ascii")
    ).hexdigest()
    tags, counts = np.unique(mesh.tags, return_counts=True)
    return {
        "name": mesh.spec.name,
        "mesh-id": combined[:16],
        "spec": asdict(mesh.spec),
        "cells": int(mesh.connect.shape[0]),
        "vertices": int(mesh.coords.shape[0]),
        "file": Path(mesh_path).name,
        "face-slot-tags": {
            str(int(tag)): int(count) for tag, count in zip(tags, counts)
        },
        "datasets": digests,
    }


# --------------------------------------------------------------------------
# validation
# --------------------------------------------------------------------------


class Report:
    def __init__(self, label):
        self.label = label
        self.failures = []
        self.notes = []

    def require(self, condition, message):
        if not condition:
            self.failures.append(message)
        return condition

    def note(self, message):
        self.notes.append(message)

    def ok(self):
        return not self.failures

    def show(self, stream=sys.stdout):
        status = "ok" if self.ok() else "FAILED"
        print(f"[{status}] {self.label}", file=stream)
        for note in self.notes:
            print(f"    {note}", file=stream)
        for failure in self.failures:
            print(f"    error: {failure}", file=stream)


def face_map(connect):
    """Map each face, identified by its sorted vertex triple, to its owners."""
    owners = {}
    for cell in range(connect.shape[0]):
        for slot, local in enumerate(FACE_VERTICES):
            key = tuple(sorted(int(connect[cell, index]) for index in local))
            owners.setdefault(key, []).append((cell, slot))
    return owners


def block_distribution(entities, ranks):
    """The contiguous distribution PUML applies when reading a mesh.

    The first ``entities % ranks`` ranks take one entity more than the rest,
    so the sizes are not simply an even split; a test that reasons about which
    cell ends up on which rank has to use this and not something close to it.
    """
    per_rank, missing = divmod(entities, ranks)
    return [per_rank + 1] * missing + [per_rank] * (ranks - missing)


def validate(coords, connect, tags, groups, identify, label, rank_counts=(2, 3, 4)):
    """Check everything SeisSol will later assume about the mesh."""
    report = Report(label)
    cells, vertices = connect.shape[0], coords.shape[0]
    report.note(f"{cells} cells, {vertices} vertices")

    report.require(
        connect.min() >= 0 and connect.max() < vertices,
        "connectivity references a vertex outside /geometry",
    )
    used = np.unique(connect)
    if used.size != vertices:
        report.note(f"{vertices - used.size} vertices are unused")

    edges = coords[connect[:, 1:]] - coords[connect[:, 0]][:, None, :]
    volume = np.linalg.det(edges) / 6.0
    report.require(
        bool(np.all(volume > 0.0)),
        f"{int(np.sum(volume <= 0.0))} cells are degenerate or left-handed",
    )
    if np.all(volume > 0.0):
        report.note(f"volume ratio max/min = {volume.max() / volume.min():.3f}")

    # face topology. On a periodic mesh the neighbour relation only exists
    # after the vertices have been identified, which is what SeisSol builds its
    # topology mesh from, so the checks below have to see the identified mesh.
    topology = connect if identify is None else identify[connect]
    owners = face_map(topology)
    overfull = {key: owner for key, owner in owners.items() if len(owner) > 2}
    report.require(
        not overfull, f"{len(overfull)} faces are shared by more than two cells"
    )

    # the boundary surface has to be closed: every edge of it is used twice
    boundary_edges = {}
    for key, owner in owners.items():
        if len(owner) != 1:
            continue
        for a in range(3):
            edge = tuple(sorted((key[a], key[(a + 1) % 3])))
            boundary_edges[edge] = boundary_edges.get(edge, 0) + 1
    dangling = [edge for edge, count in boundary_edges.items() if count != 2]
    report.require(not dangling, f"{len(dangling)} edges leave the surface open")

    # the check SeisSol performs itself
    interior_mismatch = 0
    exterior_mismatch = 0
    unknown = set()
    for key, owner in owners.items():
        shared = len(owner) == 2
        for cell, slot in owner:
            tag = int(tags[cell, slot])
            if tag not in FACE_TYPES.values():
                unknown.add(tag)
            elif tag in INTERNAL_FACE_TYPES and not shared:
                interior_mismatch += 1
            elif tag not in INTERNAL_FACE_TYPES and shared:
                exterior_mismatch += 1
    report.require(
        interior_mismatch == 0,
        f"{interior_mismatch} interior faces (regular/rupture) have no neighbour",
    )
    report.require(
        exterior_mismatch == 0,
        f"{exterior_mismatch} exterior faces have a neighbour",
    )
    report.require(not unknown, f"unknown boundary tags: {sorted(unknown)}")

    # count physical faces, not slots: an interior face is seen by two cells
    names = {value: key for key, value in FACE_TYPES.items()}
    per_tag = {}
    for owner in owners.values():
        tag = int(tags[owner[0][0], owner[0][1]])
        per_tag[tag] = per_tag.get(tag, 0) + 1
    report.note(
        "faces: "
        + ", ".join(
            f"{names.get(tag, tag)}={count}" for tag, count in sorted(per_tag.items())
        )
    )

    if identify is not None:
        report.require(
            identify.shape[0] == vertices,
            "/identify does not have one entry per vertex",
        )
        merged = vertices - int(np.unique(identify).size)
        report.require(merged > 0, "/identify is present but merges no vertices")
        report.note(f"{merged} vertices merged by periodicity")

    # contiguous decomposition, which is what PUML's `none` partitioner uses
    for ranks in rank_counts:
        if ranks > cells:
            continue
        sizes = block_distribution(cells, ranks)
        report.require(
            all(size > 0 for size in sizes), f"a rank would be empty at {ranks} ranks"
        )
        rank_of = np.repeat(np.arange(ranks), sizes)
        cut = sum(
            1
            for owner in owners.values()
            if len(owner) == 2 and rank_of[owner[0][0]] != rank_of[owner[1][0]]
        )
        report.note(f"{ranks} ranks: {sizes} cells, {cut} faces cut")

    if cells > 200:
        report.note(f"{cells} cells is beyond what this tooling is meant for")
    return report


def infer_formats(handle, boundary_format=None):
    """Resolve the boundary and topology format the way SeisSol does.

    Meshes written before the attributes existed carry neither, and SeisSol
    then guesses: rank 2 means one entry per face, rank 1 means a packed
    32-bit word. The guess cannot distinguish a packed 64-bit word, and the
    topology guess does not work at all (it looks for an attribute named
    ``identify`` where the dataset is), so a mesh without the attributes is
    reported here as the ambiguity it is.
    """
    warnings = []

    given = handle.attrs.get("boundary-format")
    if isinstance(given, bytes):
        given = given.decode("ascii")
    if boundary_format is not None:
        resolved = boundary_format
        if given is not None and given != boundary_format:
            warnings.append(
                f"boundary-format attribute says {given!r}, overridden with "
                f"{boundary_format!r}"
            )
    elif given is not None:
        resolved = given
    else:
        rank = handle["boundary"].ndim
        resolved = "i32x4" if rank == 2 else "i32"
        warnings.append(
            f"no boundary-format attribute; assuming {resolved!r} from rank {rank}, "
            "as SeisSol does"
        )
        if rank == 1 and handle["boundary"].dtype.itemsize == 8:
            warnings.append(
                "the packed word is 64 bit wide, so the assumption above is "
                "probably wrong; pass --boundary-format i64"
            )

    topology = handle.attrs.get("topology-format")
    if isinstance(topology, bytes):
        topology = topology.decode("ascii")

    return resolved, topology, warnings


def validate_file(path, label=None, boundary_format=None):
    """Validate a written mesh by reading it back."""
    import h5py

    path = Path(path)
    report = Report(label or str(path))
    with h5py.File(path, "r") as handle:
        for name in ("connect", "geometry", "group", "boundary"):
            if name not in handle:
                report.require(False, f"/{name} is missing")
                return report
        resolved, topology, warnings = infer_formats(handle, boundary_format)
        report.require(
            resolved in BOUNDARY_FORMATS, f"unknown boundary format {resolved!r}"
        )
        if not report.ok():
            return report
        connect = handle["connect"][:].astype(np.int64)
        coords = handle["geometry"][:]
        groups = handle["group"][:]
        raw = handle["boundary"][:]
        identify = handle["identify"][:] if "identify" in handle else None

    tags = decode_boundary(raw, resolved)
    inner = validate(coords, connect, tags, groups, identify, label or str(path))
    inner.notes.insert(0, f"boundary format {resolved}")
    for warning in warnings:
        inner.notes.append(f"warning: {warning}")

    if identify is not None and topology != "identify-vertex":
        # SeisSol's fallback looks for an *attribute* named identify, not the
        # dataset, so without the attribute it silently treats the mesh as
        # non-periodic and computes something else without complaining
        inner.failures.append(
            "/identify is present but topology-format is "
            f"{topology!r}; SeisSol would ignore the periodicity"
        )
    elif identify is None and topology == "identify-vertex":
        inner.failures.append(
            "topology-format claims identify-vertex but /identify is missing"
        )

    inner.failures = report.failures + inner.failures
    return inner


# --------------------------------------------------------------------------
# command line
# --------------------------------------------------------------------------


def _parse_triple(text, cast):
    parts = text.replace(",", " ").split()
    if len(parts) != 3:
        raise argparse.ArgumentTypeError(f"expected three values, got {text!r}")
    return tuple(cast(part) for part in parts)


def _parse_tag(text):
    if text in FACE_TYPES:
        return FACE_TYPES[text]
    try:
        return int(text)
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            f"unknown boundary condition {text!r}; use a number or one of "
            + ", ".join(sorted(FACE_TYPES))
        ) from error


def _parse_boundary(entries):
    parsed = []
    for entry in entries:
        for item in entry.split(","):
            if not item:
                continue
            if "=" not in item:
                raise argparse.ArgumentTypeError(
                    f"expected plane=condition, got {item!r}"
                )
            plane, tag = item.split("=", 1)
            if plane not in PLANES and plane != "all":
                raise argparse.ArgumentTypeError(
                    f"unknown plane {plane!r}; use all or one of {', '.join(PLANES)}"
                )
            parsed.append((plane, _parse_tag(tag)))
    return tuple(parsed)


def _parse_plane(entries):
    parsed = []
    for entry in entries:
        for item in entry.split(","):
            if not item:
                continue
            axis, _, index = item.partition("=")
            if axis not in AXES:
                raise argparse.ArgumentTypeError(f"unknown axis {axis!r}")
            parsed.append((axis, int(index) if index else None))
    return tuple(parsed)


def _apply_overrides(spec, args):
    changes = {}
    if args.cells is not None:
        changes["cells"] = args.cells
    if args.size is not None:
        changes["size"] = args.size
    if args.origin is not None:
        changes["origin"] = args.origin
    if args.bc:
        changes["boundary"] = spec.boundary + _parse_boundary(args.bc)
    if args.fault:
        changes["faults"] = _parse_plane(args.fault)
    if args.group_split:
        changes["group_splits"] = _parse_plane(args.group_split)
    if args.group_base is not None:
        changes["group_base"] = args.group_base
    if args.periodic:
        axes = _parse_plane(args.periodic)
        changes["periodic"] = tuple(
            any(AXES[axis] == name for name, _ in axes) for axis in range(3)
        )
    if args.skew is not None:
        changes["skew"] = args.skew
    if args.seed is not None:
        changes["seed"] = args.seed
    if args.order is not None:
        changes["order"] = args.order
    if args.boundary_format is not None:
        changes["boundary_format"] = args.boundary_format
    return replace(spec, **changes) if changes else spec


def _add_override_arguments(parser):
    parser.add_argument("--cells", type=lambda t: _parse_triple(t, int))
    parser.add_argument("--size", type=lambda t: _parse_triple(t, float))
    parser.add_argument("--origin", type=lambda t: _parse_triple(t, float))
    parser.add_argument(
        "--bc",
        action="append",
        default=[],
        metavar="PLANE=COND",
        help="boundary condition per plane, e.g. xmin=free-surface or all=5",
    )
    parser.add_argument(
        "--fault",
        action="append",
        default=[],
        metavar="AXIS[=INDEX]",
        help="rupture plane at a grid index, centred if the index is omitted",
    )
    parser.add_argument(
        "--group-split", action="append", default=[], metavar="AXIS[=INDEX]"
    )
    parser.add_argument("--group-base", type=int)
    parser.add_argument(
        "--periodic", action="append", default=[], metavar="AXIS", help="identify axis"
    )
    parser.add_argument("--skew", type=float, help="displacement as a fraction of h")
    parser.add_argument("--seed", type=int)
    parser.add_argument(
        "--order",
        choices=("lexicographic", "morton", "slab-x", "slab-y", "slab-z"),
        help="cell order, which fixes the contiguous MPI decomposition",
    )
    parser.add_argument("--boundary-format", choices=sorted(BOUNDARY_FORMATS))


def _command_list(args):
    for name, spec in PRESETS.items():
        details = [f"{spec.cell_count()} cells"]
        if spec.skew:
            details.append(f"skew {spec.skew}")
        if spec.faults:
            details.append("fault")
        if any(spec.periodic):
            details.append("periodic")
        if spec.group_splits:
            details.append("groups")
        print(f"{name:16s} {', '.join(details)}")
    return 0


def _command_generate(args):
    if args.name not in PRESETS:
        print(
            f"unknown mesh {args.name!r}; known: {', '.join(PRESETS)}", file=sys.stderr
        )
        return 2
    spec = _apply_overrides(PRESETS[args.name], args)
    mesh = build(spec)
    path = write(mesh, args.output)
    if args.xdmf:
        write_xdmf(mesh, path, path.with_suffix(".xdmf"))
    info = manifest(mesh, path)
    if args.manifest:
        Path(args.manifest).write_text(json.dumps(info, indent=2) + "\n", "utf-8")
    if args.check:
        report = validate_file(path, label=f"{spec.name} -> {path}")
        report.show()
        if not report.ok():
            return 1
    print(f"{path} ({info['cells']} cells, mesh-id {info['mesh-id']})")
    return 0


def _command_check(args):
    failed = False
    for path in args.files:
        report = validate_file(path, boundary_format=args.boundary_format)
        report.show()
        failed |= not report.ok()
    return 1 if failed else 0


def _command_self_test(args):
    failed = False
    with tempfile.TemporaryDirectory() as directory:
        for name, spec in PRESETS.items():
            reference = None
            for boundary_format in sorted(BOUNDARY_FORMATS):
                variant = replace(spec, boundary_format=boundary_format)
                mesh = build(variant)
                path = write(mesh, Path(directory) / f"{name}-{boundary_format}.h5")
                report = validate_file(path, label=f"{name} [{boundary_format}]")
                # the packing must not change what the mesh means
                decoded = decode_boundary(
                    encode_boundary(mesh.tags, boundary_format), boundary_format
                )
                report.require(
                    bool(np.array_equal(decoded, mesh.tags)),
                    "boundary packing does not survive a round trip",
                )
                if reference is None:
                    reference = manifest(mesh, path)["mesh-id"]
                report.show()
                failed |= not report.ok()
            # rebuilding has to give the same bytes, or snapshots are worthless
            rebuilt = manifest(build(spec), "x")["mesh-id"]
            if rebuilt != reference:
                print(f"[FAILED] {name}: not reproducible", file=sys.stderr)
                failed = True
    return 1 if failed else 0


def locate(coords, connect, point):
    """Return the smallest barycentric coordinate of a point in its cell.

    The value is negative if no cell contains the point, zero if it lies on a
    face, an edge or a vertex, and positive if it is strictly interior.
    """
    best = -np.inf
    for cell in connect:
        x = coords[cell]
        frame = np.column_stack([x[1] - x[0], x[2] - x[0], x[3] - x[0]])
        local = np.linalg.solve(frame, np.asarray(point) - x[0])
        smallest = min(1.0 - local.sum(), *local)
        best = max(best, smallest)
    return best


def locate_on_fault(coords, connect, tags, point, tolerance=1e-9):
    """Return the smallest barycentric coordinate of a point on the fault.

    Only faces tagged as dynamic rupture are considered, and only those whose
    plane contains the point. The value is negative if no such face contains
    it, zero if it lies on an edge or a vertex of the fault triangulation.
    """
    point = np.asarray(point, dtype=np.float64)
    best = -np.inf
    rupture = FACE_TYPES["dynamic-rupture"]
    for cell, slot in zip(*np.nonzero(tags == rupture)):
        a, b, c = coords[connect[cell, list(FACE_VERTICES[slot])]]
        normal = np.cross(b - a, c - a)
        area = np.linalg.norm(normal)
        scale = max(np.linalg.norm(b - a), np.linalg.norm(c - a))
        if abs(np.dot(point - a, normal)) / area > tolerance * scale:
            continue
        # barycentric coordinates from signed sub-areas
        weights = [
            np.dot(np.cross(c - b, point - b), normal),
            np.dot(np.cross(a - c, point - c), normal),
            np.dot(np.cross(b - a, point - a), normal),
        ]
        best = max(best, min(weights) / area**2)
    return best


def _command_locate(args):
    """Check that every receiver is strictly inside a cell of every mesh.

    A point on a face, an edge or a vertex has no well-defined cell; which one
    it is assigned to can depend on the decomposition, and since a DG
    solution is discontinuous across faces, the same run then reports
    different values on different rank counts. On these meshes that is easy
    to hit: every cube diagonal is an edge shared by all six tetrahedra.
    """
    import h5py

    points = []
    for line in Path(args.receivers).read_text(encoding="utf-8").splitlines():
        if line.strip():
            points.append([float(value) for value in line.split()])
    failed = False
    for path in args.files:
        with h5py.File(path, "r") as handle:
            coords = handle["geometry"][:]
            connect = handle["connect"][:].astype(np.int64)
            if args.fault:
                boundary_format, _, _ = infer_formats(handle)
                tags = decode_boundary(handle["boundary"][:], boundary_format)
        if args.fault:
            margins = [locate_on_fault(coords, connect, tags, p) for p in points]
        else:
            margins = [locate(coords, connect, point) for point in points]
        worst = min(margins)
        ok = worst >= args.margin
        failed |= not ok
        status = "ok" if ok else "FAILED"
        print(
            f"[{status}] {Path(path).name}: smallest barycentric coordinate "
            + ", ".join(f"{margin:.4f}" for margin in margins)
        )
    return 1 if failed else 0


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="generate the small PUML meshes used for verification"
    )
    commands = parser.add_subparsers(dest="command", required=True)

    generate = commands.add_parser("generate", help="write a mesh")
    generate.add_argument("name")
    generate.add_argument("-o", "--output", required=True)
    generate.add_argument("--manifest", help="write a JSON description next to it")
    generate.add_argument("--xdmf", action="store_true", help="also write an XDMF")
    generate.add_argument(
        "--no-check",
        dest="check",
        action="store_false",
        help="skip validation of the result",
    )
    generate.set_defaults(check=True, func=_command_generate)
    _add_override_arguments(generate)

    check = commands.add_parser("check", help="validate existing meshes")
    check.add_argument("files", nargs="+")
    check.add_argument(
        "--boundary-format",
        choices=sorted(BOUNDARY_FORMATS),
        help="override the format instead of inferring it, for older meshes "
        "that carry no attributes",
    )
    check.set_defaults(func=_command_check)

    listing = commands.add_parser("list", help="show the known meshes")
    listing.set_defaults(func=_command_list)

    locate_parser = commands.add_parser(
        "locate", help="check that receivers lie strictly inside a cell"
    )
    locate_parser.add_argument("receivers")
    locate_parser.add_argument("files", nargs="+")
    locate_parser.add_argument(
        "--fault",
        action="store_true",
        help="the points are on-fault receivers; locate them on the rupture faces",
    )
    locate_parser.add_argument(
        "--margin",
        type=float,
        default=0.01,
        help="required distance from the cell boundary, in barycentric units",
    )
    locate_parser.set_defaults(func=_command_locate)

    self_test = commands.add_parser("self-test", help="build and validate all meshes")
    self_test.set_defaults(func=_command_self_test)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
