# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""How a material's transposed coefficient matrices decompose.

Each matrix is linear in a handful of scalars the material supplies, and its
entries are those scalars times a constant. Stating which scalar sits where
lets a cell carry the coefficients and the rows of its Jacobian apart instead
of folded together, and lets the generator write the constants into the kernel
rather than read them from memory.

The scalars themselves stay on the C++ side, where the material lives; only
their placement is described here. The generated header carries the placement
of every material, not just the one a build is configured for, so that the
declarations can be checked against what each material writes itself.
"""

from dataclasses import dataclass, field
from typing import List, Tuple

from kernels import quantities

#: The two shapes the operator of a cell can be applied in where the material
#: varies inside it. Factored keeps the coefficients apart and scales the fixed
#: structures at every application; assembled folds them into one operator per
#: sample point beforehand. Which one is cheaper depends on how often a kernel
#: applies the operator and on what the machine does with the temporary the
#: fold needs, so a build chooses.
OPERATOR_FORMS = ("factored", "assembled")


@dataclass(frozen=True)
class Entry:
    """One entry of the transposed coefficient matrix of one direction."""

    coefficient: int
    dim: int
    row: int
    column: int
    factor: float = -1.0


@dataclass(frozen=True)
class AnelasticEntry:
    """One entry of the coupling block of a single relaxation mechanism.

    The column is relative to that mechanism's block, and there is no
    coefficient: the whole block is weighted by one scalar, and which one is
    the solver's decision.
    """

    dim: int
    row: int
    column_offset: int
    factor: float


@dataclass(frozen=True)
class SourceEntry:
    """One entry a single relaxation mechanism contributes to the source."""

    coefficient: int
    row: int
    column: int
    factor: float = 1.0


@dataclass(frozen=True)
class Decomposition:
    name: str
    coefficients: List[str]
    entries: List[Entry] = field(default_factory=list)
    anelastic: List[AnelasticEntry] = field(default_factory=list)
    source_coefficients: List[str] = field(default_factory=list)
    source: List[SourceEntry] = field(default_factory=list)


# ---------------------------------------------------------------- elastic

# lambda + 2 mu, lambda, mu, 1/rho, and the 1/rho that couples the shear
# stresses. The latter is a coefficient of its own because it vanishes for
# acoustic material while the first one does not.
_L2M, _LAM, _MU, _RHO, _RHO_SHEAR = range(5)

ELASTIC = Decomposition(
    name="Elastic",
    coefficients=["lambda + 2 mu", "lambda", "mu", "1/rho", "1/rho (shear)"],
    # fmt: off
    entries=[
        Entry(_L2M, 0, 6, 0), Entry(_LAM, 0, 6, 1), Entry(_LAM, 0, 6, 2),
        Entry(_MU, 0, 7, 3), Entry(_MU, 0, 8, 5),
        Entry(_RHO, 0, 0, 6), Entry(_RHO_SHEAR, 0, 3, 7), Entry(_RHO_SHEAR, 0, 5, 8),

        Entry(_LAM, 1, 7, 0), Entry(_L2M, 1, 7, 1), Entry(_LAM, 1, 7, 2),
        Entry(_MU, 1, 6, 3), Entry(_MU, 1, 8, 4),
        Entry(_RHO, 1, 1, 7), Entry(_RHO_SHEAR, 1, 3, 6), Entry(_RHO_SHEAR, 1, 4, 8),

        Entry(_LAM, 2, 8, 0), Entry(_LAM, 2, 8, 1), Entry(_L2M, 2, 8, 2),
        Entry(_MU, 2, 7, 4), Entry(_MU, 2, 6, 5),
        Entry(_RHO, 2, 2, 8), Entry(_RHO_SHEAR, 2, 5, 6), Entry(_RHO_SHEAR, 2, 4, 7),
    ],
    # fmt: on
)

# ---------------------------------------------------------------- acoustic

ACOUSTIC = Decomposition(
    name="Acoustic",
    coefficients=["lambda", "1/rho"],
    # fmt: off
    entries=[
        Entry(0, 0, 1, 0), Entry(1, 0, 0, 1),
        Entry(0, 1, 2, 0), Entry(1, 1, 0, 2),
        Entry(0, 2, 3, 0), Entry(1, 2, 0, 3),
    ],
    # fmt: on
)

# ------------------------------------------------------------- anisotropic

# fmt: off
_VOIGT = [
    "c11", "c12", "c13", "c14", "c15", "c16",
    "c22", "c23", "c24", "c25", "c26",
    "c33", "c34", "c35", "c36",
    "c44", "c45", "c46",
    "c55", "c56",
    "c66",
]
# fmt: on
_ANISO_RHO = len(_VOIGT)


@dataclass(frozen=True)
class FluxEntry:
    """One entry of a flux operator, as a scalar of the face times a factor."""

    coefficient: int
    row: int
    column: int
    factor: float


def flux_decomposition(blocks):
    """The scalars the flux operator of a face is linear in, and where they sit.

    That the operator is a handful of scalars at all is measured, not derived:
    over six hundred random isotropic material pairs the nine-quantity operator
    occupies sixteen of its eighty-one entries, and those sixteen span a space
    of exactly ten dimensions, the singular values falling from 1e-3 to 1e-17.
    Entries outside the ten are equal to one of them, not merely a multiple,
    because the two directions in the face enter the same way. The scalars
    themselves are whatever the Riemann solver makes of the two materials at
    the point; they are not linear in either, which is why they are read off a
    computed operator rather than assembled from material parameters.

    Which scalars a layout has, and which entries each fills, does follow from
    that layout. The traction group of a face supplies the component along its
    normal, the two that go with the directions in the face, and the two normal
    to those; the velocity group supplies the first three of those. Pairing them
    gives the ``p`` block the normal picks out, the ``s`` block the two face
    directions share, and -- once per relaxation mechanism -- the ``a`` block of
    columns that mechanism couples into both.

    A layout without a shear traction has no shear block at all: there is no
    shear wave to carry across the face, so the two tangential velocities are
    left uncoupled, which is what the operator of an acoustic medium shows.
    Coefficients that end up with no entry are dropped, so a layout gets exactly
    the scalars it has positions for.

    Returns the coefficient names, where each is read off a computed operator,
    and the entries each one fills.
    """
    extra = quantities.extra_face_blocks(blocks)
    if extra:
        names = ", ".join(block.group.name for block in extra)
        raise ValueError(
            "the flux operator of a layout with further face-local rows "
            f"({names}) couples more than the traction and the velocity of one "
            "medium, and is not these scalars"
        )

    traction = quantities.face_block(blocks, quantities.FaceRole.TRACTION)
    velocity = quantities.face_block(blocks, quantities.FaceRole.VELOCITY)
    tractionNormal, tractionShears, tractionTransverse = quantities.face_components(
        traction
    )
    velocityNormal, velocityShears, _ = quantities.face_components(velocity)

    groups = [
        ("pNormalNormal", ((tractionNormal, tractionNormal),)),
        # the two transverse normal stresses enter alike
        (
            "pNormalTransverse",
            tuple((tractionNormal, column) for column in tractionTransverse),
        ),
        ("pNormalVelocity", ((tractionNormal, velocityNormal),)),
        ("pVelocityNormal", ((velocityNormal, tractionNormal),)),
        (
            "pVelocityTransverse",
            tuple((velocityNormal, column) for column in tractionTransverse),
        ),
        ("pVelocityVelocity", ((velocityNormal, velocityNormal),)),
    ]
    if tractionShears:
        # and so do the two directions in the face
        groups += [
            ("sShearShear", tuple(zip(tractionShears, tractionShears))),
            ("sShearVelocity", tuple(zip(tractionShears, velocityShears))),
            ("sVelocityShear", tuple(zip(velocityShears, tractionShears))),
            ("sVelocityVelocity", tuple(zip(velocityShears, velocityShears))),
        ]

    mechanisms = [block for block in blocks if block.mechanism is not None]
    for index, block in enumerate(mechanisms):
        mechanismNormal, mechanismShears, _ = quantities.face_components(block)
        # one set per block, since what a mechanism carries is its own number
        suffix = f"[{index}]" if len(mechanisms) > 1 else ""
        groups += [
            (f"aNormalNormal{suffix}", ((tractionNormal, mechanismNormal),)),
            (
                f"aShearShear{suffix}",
                tuple(zip(tractionShears, mechanismShears)),
            ),
            (f"aVelocityNormal{suffix}", ((velocityNormal, mechanismNormal),)),
            (
                f"aVelocityShear{suffix}",
                tuple(zip(velocityShears, mechanismShears)),
            ),
        ]

    groups = [(name, positions) for name, positions in groups if positions]
    names = tuple(name for name, _ in groups)
    # the first position of a coefficient is where it is read off; the others
    # are equal to it, which is what makes them one coefficient
    sources = {name: positions[0] for name, positions in groups}
    entries = tuple(
        FluxEntry(index, row, column, 1.0)
        for index, (_, positions) in enumerate(groups)
        for row, column in positions
    )
    return names, sources, entries


def _c(name: str) -> int:
    return _VOIGT.index(name)


def _anisotropic_entries() -> List[Entry]:
    # rows 6..8 of direction d hold the stress response to a velocity
    # gradient; the constant at (row, column) is the one the Voigt pair of the
    # two picks out. The three directions share this block entirely, unlike
    # the isotropic case, where they are disjointly occupied.
    # fmt: off
    stress = {
        0: [("c11", "c16", "c15"), ("c12", "c26", "c25"), ("c13", "c36", "c35"),
            ("c16", "c66", "c56"), ("c14", "c46", "c45"), ("c15", "c56", "c55")],
        1: [("c16", "c12", "c14"), ("c26", "c22", "c24"), ("c36", "c23", "c34"),
            ("c66", "c26", "c46"), ("c46", "c24", "c44"), ("c56", "c25", "c45")],
        2: [("c15", "c14", "c13"), ("c25", "c24", "c23"), ("c35", "c34", "c33"),
            ("c56", "c46", "c36"), ("c45", "c44", "c34"), ("c55", "c45", "c35")],
    }
    # fmt: on
    # the velocity rows take 1/rho, at the one column of each direction
    velocity = {
        0: [(0, 6), (3, 7), (5, 8)],
        1: [(3, 6), (1, 7), (4, 8)],
        2: [(5, 6), (4, 7), (2, 8)],
    }

    entries = []
    for dim in (0, 1, 2):
        for column, names in enumerate(stress[dim]):
            for offset, name in enumerate(names):
                entries.append(Entry(_c(name), dim, 6 + offset, column))
        for row, column in velocity[dim]:
            entries.append(Entry(_ANISO_RHO, dim, row, column))
    return entries


ANISOTROPIC = Decomposition(
    name="Anisotropic",
    coefficients=_VOIGT + ["1/rho"],
    entries=_anisotropic_entries(),
)

# ------------------------------------------------------------- poroelastic

# The frame is isotropic, so the 6x6 cBar holds three distinct values and the
# six-vector alpha one. Most of what the matrices index is therefore
# structurally zero and does not appear here at all.
_CBAR_DIAG, _CBAR_OFF, _CBAR_SHEAR, _M_ALPHA, _M, _RHO1, _RHO2, _B1, _B2 = range(9)

# AT(row, column) of the stress rows reads cBar at this row index
_PORO_COLUMN_TO_CBAR = [0, 1, 2, 5, 3, 4]


def _cbar(i: int, j: int):
    if i < 3 and j < 3:
        return _CBAR_DIAG if i == j else _CBAR_OFF
    if i >= 3 and i == j:
        return _CBAR_SHEAR
    return None


def _poroelastic_entries() -> List[Entry]:
    # fmt: off
    spec = {
        0: dict(fluid=[(0, 6, _RHO1), (0, 10, _RHO2), (3, 7, _RHO1),
                       (3, 11, _RHO2), (5, 8, _RHO1), (5, 12, _RHO2)],
                stress=[(6, 0), (7, 5), (8, 4)],
                beta=[(9, 6, _B1), (9, 10, _B2)], pressure=10),
        1: dict(fluid=[(1, 7, _RHO1), (1, 11, _RHO2), (3, 6, _RHO1),
                       (3, 10, _RHO2), (4, 8, _RHO1), (4, 12, _RHO2)],
                stress=[(6, 5), (7, 1), (8, 3)],
                beta=[(9, 7, _B1), (9, 11, _B2)], pressure=11),
        2: dict(fluid=[(2, 8, _RHO1), (2, 12, _RHO2), (4, 7, _RHO1),
                       (4, 11, _RHO2), (5, 6, _RHO1), (5, 10, _RHO2)],
                stress=[(6, 4), (7, 3), (8, 2)],
                beta=[(9, 8, _B1), (9, 12, _B2)], pressure=12),
    }
    # fmt: on

    entries = []
    for dim in (0, 1, 2):
        block = spec[dim]
        for row, column, coefficient in block["fluid"]:
            entries.append(Entry(coefficient, dim, row, column))
        for row, cbar_column in block["stress"]:
            for column in range(6):
                coefficient = _cbar(_PORO_COLUMN_TO_CBAR[column], cbar_column)
                if coefficient is not None:
                    entries.append(Entry(coefficient, dim, row, column))
            if cbar_column < 3:
                entries.append(Entry(_M_ALPHA, dim, row, 9, 1.0))
        for row, column, coefficient in block["beta"]:
            entries.append(Entry(coefficient, dim, row, column))
        pressure = block["pressure"]
        for column in range(6):
            if _PORO_COLUMN_TO_CBAR[column] < 3:
                entries.append(Entry(_M_ALPHA, dim, pressure, column))
        entries.append(Entry(_M, dim, pressure, 9, 1.0))
    return entries


POROELASTIC = Decomposition(
    name="PoroElastic",
    coefficients=[
        "cBar(0,0)",
        "cBar(1,0)",
        "cBar(3,3)",
        "M alpha",
        "M",
        "1/rho1",
        "1/rho2",
        "beta1/rho1",
        "beta2/rho2",
    ],
    entries=_poroelastic_entries(),
    # the Biot drag: the relative motion of fluid and solid relaxes against the
    # two densities, one scalar each
    source_coefficients=[
        "beta1 eta / (rho1 kappa)",
        "beta2 eta / (rho2 kappa)",
    ],
    # fmt: off
    source=[
        SourceEntry(0, 10, 6), SourceEntry(0, 11, 7), SourceEntry(0, 12, 8),
        SourceEntry(1, 10, 10), SourceEntry(1, 11, 11), SourceEntry(1, 12, 12),
    ],
    # fmt: on
)

# -------------------------------------------------------- viscous materials

VISCOELASTIC = Decomposition(
    name="ViscoElastic",
    coefficients=[],  # the flux is the base material's
    anelastic=[
        AnelasticEntry(0, 6, 0, -1.0),
        AnelasticEntry(0, 7, 3, -0.5),
        AnelasticEntry(0, 8, 5, -0.5),
        AnelasticEntry(1, 7, 1, -1.0),
        AnelasticEntry(1, 6, 3, -0.5),
        AnelasticEntry(1, 8, 4, -0.5),
        AnelasticEntry(2, 8, 2, -1.0),
        AnelasticEntry(2, 7, 4, -0.5),
        AnelasticEntry(2, 6, 5, -0.5),
    ],
    source_coefficients=["theta[0]", "theta[1]", "theta[2]"],
    source=[
        SourceEntry(0, 0, 0),
        SourceEntry(1, 1, 0),
        SourceEntry(1, 2, 0),
        SourceEntry(1, 0, 1),
        SourceEntry(0, 1, 1),
        SourceEntry(1, 2, 1),
        SourceEntry(1, 0, 2),
        SourceEntry(1, 1, 2),
        SourceEntry(0, 2, 2),
        SourceEntry(2, 3, 3),
        SourceEntry(2, 4, 4),
        SourceEntry(2, 5, 5),
    ],
)

VISCOACOUSTIC = Decomposition(
    name="ViscoAcoustic",
    coefficients=[],
    anelastic=[
        AnelasticEntry(0, 1, 0, -1.0),
        AnelasticEntry(1, 2, 0, -1.0),
        AnelasticEntry(2, 3, 0, -1.0),
    ],
    source_coefficients=["theta[0]"],
    source=[SourceEntry(0, 0, 0)],
)


ALL = [ELASTIC, ACOUSTIC, ANISOTROPIC, POROELASTIC, VISCOELASTIC, VISCOACOUSTIC]


BY_EQUATION = {
    "elastic": ELASTIC,
    "acoustic": ACOUSTIC,
    "anisotropic": ANISOTROPIC,
    "poroelastic": POROELASTIC,
    "viscoelastic": VISCOELASTIC,
    "viscoacoustic": VISCOACOUSTIC,
}

#: Which material a viscous one takes its flux, and so its decomposition, from.
BASE_OF = {"viscoelastic": ELASTIC, "viscoacoustic": ACOUSTIC}


#: A coefficient read off the material. It varies from cell to cell, and
#: within a cell wherever the material is sampled at the nodal points.
MATERIAL = "Material"
#: A coefficient that is one number for the whole domain, fixed once the run is
#: set up -- the relaxation frequencies, which follow only the frequency band.
#: It is not a build constant, so it cannot be written into the kernel, but it
#: does not belong in every cell either.
GLOBAL = "Global"


def composed(
    equation: str,
    solver: str,
    mechanisms: int,
    elastic_quantities: int,
    per_mechanism: int,
) -> Tuple[int, List[Entry]]:
    """The decomposition of the operator a solver applies, coefficients first.

    Mirrors what SolverSetup does on the C++ side: the material's own entries,
    then the coupling block of every relaxation mechanism at its own columns.
    Which scalar weights a block is the solver's decision -- one relaxation
    frequency per block, or a single block of unit weight where the solver
    holds the frequencies elsewhere.
    """
    decomposition = BY_EQUATION[equation]
    base = BASE_OF.get(equation, decomposition)
    count = len(base.coefficients)
    entries = list(base.entries)
    # every coefficient a material declares is read off that material, so it is
    # a field; the weights a solver adds for its relaxation blocks are not
    origins = [MATERIAL] * count

    if mechanisms > 0 and decomposition.anelastic:
        if solver == "linearckanelastic":
            blocks, weights = 1, [count]
            count += 1
            origins.append(GLOBAL)
        else:
            blocks, weights = mechanisms, [count + m for m in range(mechanisms)]
            count += mechanisms
            origins += [GLOBAL] * mechanisms
        for block in range(blocks):
            column = elastic_quantities + block * per_mechanism
            for entry in decomposition.anelastic:
                entries.append(
                    Entry(
                        weights[block],
                        entry.dim,
                        entry.row,
                        column + entry.column_offset,
                        entry.factor,
                    )
                )

    return count, entries, origins


def source_composed(
    equation: str,
    solver: str,
    mechanisms: int,
    shape: Tuple[int, ...],
    elastic_quantities: int = 0,
    per_mechanism: int = 0,
):
    """The decomposition of the source term a solver applies.

    The same idea as `composed`, over the tensor each solver states its source
    in. A material without relaxation states it as a matrix and its scalars are
    read off the material. With relaxation there is one block per mechanism,
    and where the solver keeps the mechanism in a dimension of its own the
    relaxation frequency sits outside the source altogether; where it folds the
    blocks into one matrix, the frequency is on that block's diagonal and
    becomes a scalar of the decomposition -- one that follows the frequency
    band rather than the material.

    Returns the number of scalars, where each one's entries go, and where each
    one comes from.
    """
    decomposition = BY_EQUATION[equation]
    if not decomposition.source_coefficients:
        return 0, {}, []

    perBlock = len(decomposition.source_coefficients)
    split = len(shape) == 3
    blocks = mechanisms if mechanisms > 0 else 1
    # the folded form carries the relaxation frequency of every block with it
    perMechanism = perBlock if split or mechanisms == 0 else perBlock + 1

    values = {}
    origins = []
    for block in range(blocks):
        offset = elastic_quantities + block * per_mechanism
        for entry in decomposition.source:
            if split:
                index = (entry.row, block, entry.column)
            elif mechanisms > 0:
                index = (offset + entry.row, entry.column)
            else:
                index = (entry.row, entry.column)
            key = (block * perMechanism + entry.coefficient,) + index
            values[key] = values.get(key, 0.0) + entry.factor
        origins += [MATERIAL] * perBlock

        if not split and mechanisms > 0:
            relaxation = block * perMechanism + perBlock
            for i in range(per_mechanism):
                values[(relaxation, offset + i, offset + i)] = -1.0
            origins.append(GLOBAL)

    return perMechanism * blocks, values, origins


def structure_values(
    coefficient_count: int, entries: List[Entry], star_shape: Tuple[int, int]
) -> Tuple[Tuple[int, ...], dict]:
    """The decomposition as a tensor the generator can write into a kernel.

    Shape is (coefficients, 3, quantities, quantities); the values are the
    constants, which are plus or minus one almost everywhere. Handing these to
    a tensor with the immediate addressing mode is what turns the assembly
    into one product per entry: a factor of one is not a multiplication, and
    the zeros never become operations at all.
    """
    shape = (coefficient_count, 3) + tuple(star_shape)
    values = {
        (entry.coefficient, entry.dim, entry.row, entry.column): entry.factor
        for entry in entries
    }
    return shape, values


def _format(value: float) -> str:
    return f"{value:.1f}" if value == int(value) else repr(value)


def _flux_table(kind: str, name: str, rows: List[str]) -> List[str]:
    """Like :func:`_table`, but declared even when there is nothing in it: the
    code that reads these tables names them whether or not the build it is
    compiled for has a decomposition to read."""
    if not rows:
        return [f"inline constexpr std::array<{kind}, 0> {name}{{}};\n\n"]
    return _table(kind, name, rows)


def _table(kind: str, name: str, rows: List[str]) -> List[str]:
    if not rows:
        return []
    return (
        [f"inline constexpr std::array<{kind}, {len(rows)}> {name}{{{{\n"]
        + [f"    {row},\n" for row in rows]
        + ["}};\n\n"]
    )


def generate(
    path: str,
    solver_count: int = None,
    solver_origins: List[str] = None,
    solver_source_count: int = None,
    solver_source_deviations: int = None,
    material_samples: int = 1,
    face_permutations=(),
    flux_blocks=(),
    flux_decomposes: bool = True,
) -> None:
    """Write the declarations of every material into a C++ header.

    With a count for the configured build, the header also states how many
    coefficients its operator has, so that a cell can be sized without
    instantiating the solver's declaration.
    """
    lines = [
        "// SPDX-FileCopyrightText: 2026 SeisSol Group\n",
        "//\n",
        "// SPDX-License-Identifier: BSD-3-Clause\n",
        "\n",
        "#ifndef SEISSOL_GENERATEDCODE_COEFFICIENTS_H_\n",
        "#define SEISSOL_GENERATEDCODE_COEFFICIENTS_H_\n",
        "\n",
        '#include "Model/CommonDatastructures.h"\n',
        "\n",
        "#include <array>\n",
        "#include <cstddef>\n",
        "\n",
        "namespace seissol::generated {\n",
        "\n",
    ]

    for decomposition in ALL:
        name = decomposition.name
        if decomposition.coefficients:
            lines.append(f"// {name}: " + ", ".join(decomposition.coefficients) + "\n")
            lines.append(
                f"inline constexpr std::size_t {name}NumCoefficients = "
                f"{len(decomposition.coefficients)};\n"
            )
            lines += _table(
                "model::CoefficientEntry",
                f"{name}CoefficientEntries",
                [
                    f"{{{e.coefficient}, {e.dim}, {e.row}, {e.column}, {_format(e.factor)}}}"
                    for e in decomposition.entries
                ],
            )
        lines += _table(
            "model::AnelasticCoefficientEntry",
            f"{name}AnelasticEntries",
            [
                f"{{{e.dim}, {e.row}, {e.column_offset}, {_format(e.factor)}}}"
                for e in decomposition.anelastic
            ],
        )
        if decomposition.source_coefficients:
            lines.append(
                f"// {name} source: "
                + ", ".join(decomposition.source_coefficients)
                + "\n"
            )
            lines.append(
                f"inline constexpr std::size_t {name}NumSourceCoefficients = "
                f"{len(decomposition.source_coefficients)};\n"
            )
            lines += _table(
                "model::SourceCoefficientEntry",
                f"{name}SourceEntries",
                [
                    f"{{{e.coefficient}, {e.row}, {e.column}, {_format(e.factor)}}}"
                    for e in decomposition.source
                ],
            )

    if face_permutations:
        lines += [
            "// how the nodes of a face are renumbered between the two cells\n",
            "// sharing it, one row per reparametrisation\n",
            f"inline constexpr std::size_t FaceOrientations = {len(face_permutations)};\n",
            f"inline constexpr std::size_t FaceNodes = {len(face_permutations[0])};\n",
            "inline constexpr std::array<std::array<std::size_t, "
            f"{len(face_permutations[0])}>, {len(face_permutations)}> "
            "FaceOrientationPermutations{{\n",
        ]
        for perm in face_permutations:
            lines.append("    {{" + ", ".join(str(i) for i in perm) + "}},\n")
        lines.append("}};\n\n")

    if flux_decomposes:
        flux_names, flux_sources, flux_entries = flux_decomposition(flux_blocks)
    else:
        # nothing to state, and a build that would need it is refused where the
        # kernels are generated
        flux_names, flux_sources, flux_entries = (), {}, ()
    lines.append("// the flux operator of a face, as scalars of that face times\n")
    lines.append("// fixed entries: " + (", ".join(flux_names) or "none") + "\n")
    lines.append(
        f"inline constexpr std::size_t FluxNumCoefficients = {len(flux_names)};\n"
    )
    lines += _flux_table(
        "model::FluxCoefficientEntry",
        "FluxCoefficientEntries",
        [
            f"{{{e.coefficient}, {e.row}, {e.column}, {_format(e.factor)}}}"
            for e in flux_entries
        ],
    )
    lines += _flux_table(
        "model::FluxCoefficientSource",
        "FluxCoefficientSources",
        [
            f"{{{flux_sources[name][0]}, {flux_sources[name][1]}}}"
            for name in flux_names
        ],
    )

    lines += [
        "// how many samples of the material a cell carries. One, where the\n",
        "// material does not vary inside a cell.\n",
        f"inline constexpr std::size_t MaterialSampleCount = {material_samples};\n",
        "\n",
    ]

    if solver_count is not None:
        lines += [
            "// the operator the configured solver applies, material and any\n",
            "// relaxation blocks together\n",
            f"inline constexpr std::size_t SolverNumCoefficients = {solver_count};\n",
            "\n",
            "// the same for its source term, zero where it has none\n",
            "inline constexpr std::size_t SolverNumSourceCoefficients = "
            f"{solver_source_count or 0};\n",
            "\n",
            "// how many of those a cell carries as the deviation of a sample\n",
            "// point from the cell, for a solver that factorises the source\n",
            "// term into a solve it does once per cell\n",
            "inline constexpr std::size_t SolverNumSourceDeviations = "
            f"{solver_source_deviations or 0};\n",
            "\n",
        ]
    if solver_count is not None:
        lines += [
            "// where each of them comes from, which decides whether a cell has\n",
            "// to carry it. Empty where the build does not factor the star, and\n",
            "// so carries no coefficients at all.\n",
        ]
        rows = [
            f"model::CoefficientOrigin::{origin}" for origin in (solver_origins or [])
        ]
        if rows:
            lines += _table(
                "model::CoefficientOrigin", "SolverCoefficientOrigins", rows
            )
        else:
            # the name has to exist even then: it is looked up unconditionally
            lines.append(
                "inline constexpr std::array<model::CoefficientOrigin, 0> "
                "SolverCoefficientOrigins{};\n\n"
            )

    lines += [
        "} // namespace seissol::generated\n",
        "\n",
        "#endif // SEISSOL_GENERATEDCODE_COEFFICIENTS_H_\n",
    ]

    with open(path, "w") as file:
        file.writelines(lines)
