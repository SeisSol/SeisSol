# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Declarative description of a material's quantity layout.

A material states its layout as two lists of groups: the primary quantities,
and -- for materials with relaxation -- one mechanism block. How often the
mechanism block is instantiated along the quantity axis is a property of the
solver, not of the material: the fused layout repeats it once per mechanism,
the split layout carries a single copy and puts the mechanism index into a
separate tensor dimension. The two lists are therefore combined in exactly one
place, :func:`layout`, and everything that follows from the layout -- rotation
sparsity, the face selectors, the offset of the velocity components -- is
derived here rather than spelled out per material.
"""

from dataclasses import dataclass
from enum import Enum
from typing import Optional

import numpy as np


class QuantityKind(Enum):
    """How a group of quantities transforms under the face rotation."""

    SCALAR = "scalar"
    VECTOR = "vector"
    SYM_TENSOR2 = "sym_tensor2"


#: Number of components each kind occupies along the quantity axis.
EXTENT = {
    QuantityKind.SCALAR: 1,
    QuantityKind.VECTOR: 3,
    QuantityKind.SYM_TENSOR2: 6,
}

#: Component order of a symmetric second-order tensor.
VOIGT_ORDER = ("xx", "yy", "zz", "xy", "yz", "xz")

#: Components of a symmetric second-order tensor that carry the traction of a
#: face whose normal is the local x axis, i.e. sigma_xx, sigma_xy, sigma_xz.
SYM_TENSOR2_TRACTION = (0, 3, 5)

#: Number of rows in the face-local traction and velocity vectors.
FACE_VECTOR_ROWS = 3


class FaceRole(Enum):
    """Whether and how a group takes part in the face-local quantities.

    ``TRACTION`` and ``VELOCITY`` mark the groups that carry the mechanical
    traction and velocity vectors. ``EXTRA_*`` groups append one further row
    each, below those vectors; only their normal component participates.
    """

    NONE = "none"
    TRACTION = "traction"
    VELOCITY = "velocity"
    EXTRA_TRACTION = "extra_traction"
    EXTRA_VELOCITY = "extra_velocity"


@dataclass(frozen=True)
class QuantityGroup:
    """A contiguous run of quantities that share a transformation behaviour."""

    name: str
    kind: QuantityKind
    role: FaceRole = FaceRole.NONE

    @property
    def extent(self) -> int:
        return EXTENT[self.kind]


@dataclass(frozen=True)
class Block:
    """One instantiation of a group, placed at a concrete offset."""

    group: QuantityGroup
    offset: int
    #: Index of the mechanism this block belongs to, or None for a primary group.
    mechanism: Optional[int] = None

    @property
    def extent(self) -> int:
        return self.group.extent

    @property
    def slice(self) -> slice:
        return slice(self.offset, self.offset + self.extent)


def layout(primary, mechanism=(), repetitions=0):
    """Places the primary groups, then ``repetitions`` copies of the mechanism
    block, and returns the resulting blocks in quantity order."""
    if repetitions < 0:
        raise ValueError("repetitions must not be negative")
    blocks = []
    offset = 0
    for group in primary:
        blocks.append(Block(group, offset))
        offset += group.extent
    for rep in range(repetitions):
        for group in mechanism:
            blocks.append(Block(group, offset, rep))
            offset += group.extent
    return blocks


def total_extent(blocks) -> int:
    """Number of quantities the blocks occupy."""
    return sum(block.extent for block in blocks)


def rotation_spp(blocks):
    """Sparsity of the face rotation matrix: one dense block per group."""
    size = total_extent(blocks)
    spp = np.zeros((size, size), dtype=bool)
    for block in blocks:
        spp[block.slice, block.slice] = True
    return spp


def _face_vector_rows(block):
    """(row, quantity) pairs feeding the face-local vector of ``block``."""
    kind = block.group.kind
    if kind is QuantityKind.SYM_TENSOR2:
        return [
            (row, block.offset + component)
            for row, component in enumerate(SYM_TENSOR2_TRACTION)
        ]
    if kind is QuantityKind.VECTOR:
        return [(row, block.offset + row) for row in range(FACE_VECTOR_ROWS)]
    if kind is QuantityKind.SCALAR:
        # An isotropic stress carries a normal traction only, so the shear rows
        # of the face-local vector stay structurally empty.
        return [(0, block.offset)]
    raise ValueError(f"no face vector defined for {kind}")


def _selector(blocks, role, extra_role):
    main = [block for block in blocks if block.group.role is role]
    extra = [block for block in blocks if block.group.role is extra_role]
    if len(main) != 1:
        raise ValueError(f"expected exactly one {role.value} group, got {len(main)}")
    selector = np.zeros((FACE_VECTOR_ROWS + len(extra), total_extent(blocks)))
    for row, quantity in _face_vector_rows(main[0]):
        selector[row, quantity] = 1
    for index, block in enumerate(extra):
        # Only the normal component of an extra group takes part.
        selector[FACE_VECTOR_ROWS + index, block.offset] = 1
    return selector


def traction_selector(blocks):
    """Maps the quantities onto the face-local traction vector."""
    return _selector(blocks, FaceRole.TRACTION, FaceRole.EXTRA_TRACTION)


def velocity_selector(blocks):
    """Maps the quantities onto the face-local velocity vector."""
    return _selector(blocks, FaceRole.VELOCITY, FaceRole.EXTRA_VELOCITY)


def role_offset(blocks, role) -> int:
    """Quantity index at which the group with ``role`` starts."""
    for block in blocks:
        if block.group.role is role:
            return block.offset
    raise ValueError(f"no {role.value} group in this layout")


def well_formed(blocks, num_quantities=None):
    """Raises if the blocks cannot describe a usable layout."""
    names = [block.group.name for block in blocks if block.mechanism is None]
    if any(not name for name in names):
        raise ValueError("every group needs a name")
    if len(set(names)) != len(names):
        raise ValueError(f"duplicate group names: {names}")
    for role in (FaceRole.TRACTION, FaceRole.VELOCITY):
        count = sum(1 for block in blocks if block.group.role is role)
        if count > 1:
            raise ValueError(f"at most one {role.value} group, got {count}")
    if num_quantities is not None and total_extent(blocks) != num_quantities:
        raise ValueError(
            f"layout covers {total_extent(blocks)} quantities, expected {num_quantities}"
        )
    return True


_CXX_KIND = {
    QuantityKind.SCALAR: "Scalar",
    QuantityKind.VECTOR: "Vector",
    QuantityKind.SYM_TENSOR2: "SymTensor2",
}


def emit_header(aderdg, output_dir):
    """Writes the layout the codegen used, for the C++ side to check against.

    The two sides declare their groups independently -- the codegen cannot read
    C++ and the material headers are not generated. What this file buys is that
    a disagreement becomes a compile error instead of a wrong rotation matrix.
    """
    import os

    def render(name, blocks):
        kinds = ", ".join(f"QuantityKind::{_CXX_KIND[b.group.kind]}" for b in blocks)
        return (
            f"inline constexpr std::array<QuantityKind, {len(blocks)}> {name}"
            "{" + kinds + "};\n"
        )

    lines = [
        "// SPDX-FileCopyrightText: 2026 SeisSol Group",
        "//",
        "// SPDX-License-Identifier: BSD-3-Clause",
        "// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/",
        "//",
        "// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff",
        "",
        "// Generated by codegen/kernels/quantities.py. Do not edit.",
        "",
        "#pragma once",
        "",
        '#include "Model/Quantities.h"',
        "",
        "#include <array>",
        "",
        "namespace seissol::generated {",
        "",
        "using seissol::model::QuantityKind;",
        "",
        render("RotationGroupKinds", aderdg.extendedBlocks()),
        render("InverseRotationGroupKinds", aderdg.inverseRotationBlocks()),
        "} // namespace seissol::generated",
        "",
    ]
    with open(os.path.join(output_dir, "quantities.h"), "w") as file:
        file.write("\n".join(lines))
