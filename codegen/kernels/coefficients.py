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
)

# ---------------------------------------------------------------- acoustic

ACOUSTIC = Decomposition(
    name="Acoustic",
    coefficients=["lambda", "1/rho"],
    entries=[
        Entry(0, 0, 1, 0), Entry(1, 0, 0, 1),
        Entry(0, 1, 2, 0), Entry(1, 1, 0, 2),
        Entry(0, 2, 3, 0), Entry(1, 2, 0, 3),
    ],
)

# ------------------------------------------------------------- anisotropic

_VOIGT = [
    "c11", "c12", "c13", "c14", "c15", "c16",
    "c22", "c23", "c24", "c25", "c26",
    "c33", "c34", "c35", "c36",
    "c44", "c45", "c46",
    "c55", "c56",
    "c66",
]
_ANISO_RHO = len(_VOIGT)


def _c(name: str) -> int:
    return _VOIGT.index(name)


def _anisotropic_entries() -> List[Entry]:
    # rows 6..8 of direction d hold the stress response to a velocity
    # gradient; the constant at (row, column) is the one the Voigt pair of the
    # two picks out. The three directions share this block entirely, unlike
    # the isotropic case, where they are disjointly occupied.
    stress = {
        0: [("c11", "c16", "c15"), ("c12", "c26", "c25"), ("c13", "c36", "c35"),
            ("c16", "c66", "c56"), ("c14", "c46", "c45"), ("c15", "c56", "c55")],
        1: [("c16", "c12", "c14"), ("c26", "c22", "c24"), ("c36", "c23", "c34"),
            ("c66", "c26", "c46"), ("c46", "c24", "c44"), ("c56", "c25", "c45")],
        2: [("c15", "c14", "c13"), ("c25", "c24", "c23"), ("c35", "c34", "c33"),
            ("c56", "c46", "c36"), ("c45", "c44", "c34"), ("c55", "c45", "c35")],
    }
    # the velocity rows take 1/rho, at the one column of each direction
    velocity = {0: [(0, 6), (3, 7), (5, 8)],
                1: [(3, 6), (1, 7), (4, 8)],
                2: [(5, 6), (4, 7), (2, 8)]}

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
    coefficients=["cBar(0,0)", "cBar(1,0)", "cBar(3,3)", "M alpha", "M",
                  "1/rho1", "1/rho2", "beta1/rho1", "beta2/rho2"],
    entries=_poroelastic_entries(),
)

# -------------------------------------------------------- viscous materials

VISCOELASTIC = Decomposition(
    name="ViscoElastic",
    coefficients=[],  # the flux is the base material's
    anelastic=[
        AnelasticEntry(0, 6, 0, -1.0), AnelasticEntry(0, 7, 3, -0.5),
        AnelasticEntry(0, 8, 5, -0.5),
        AnelasticEntry(1, 7, 1, -1.0), AnelasticEntry(1, 6, 3, -0.5),
        AnelasticEntry(1, 8, 4, -0.5),
        AnelasticEntry(2, 8, 2, -1.0), AnelasticEntry(2, 7, 4, -0.5),
        AnelasticEntry(2, 6, 5, -0.5),
    ],
    source_coefficients=["theta[0]", "theta[1]", "theta[2]"],
    source=[
        SourceEntry(0, 0, 0), SourceEntry(1, 1, 0), SourceEntry(1, 2, 0),
        SourceEntry(1, 0, 1), SourceEntry(0, 1, 1), SourceEntry(1, 2, 1),
        SourceEntry(1, 0, 2), SourceEntry(1, 1, 2), SourceEntry(0, 2, 2),
        SourceEntry(2, 3, 3), SourceEntry(2, 4, 4), SourceEntry(2, 5, 5),
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


def composed(equation: str,
             solver: str,
             mechanisms: int,
             elastic_quantities: int,
             per_mechanism: int) -> Tuple[int, List[Entry]]:
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
                entries.append(Entry(weights[block], entry.dim, entry.row,
                                     column + entry.column_offset, entry.factor))

    return count, entries, origins


def structure_values(coefficient_count: int,
                     entries: List[Entry],
                     star_shape: Tuple[int, int]) -> Tuple[Tuple[int, ...], dict]:
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


def _table(kind: str, name: str, rows: List[str]) -> List[str]:
    if not rows:
        return []
    return ([f"inline constexpr std::array<{kind}, {len(rows)}> {name}{{{{\n"]
            + [f"    {row},\n" for row in rows]
            + ["}};\n\n"])


def generate(path: str, solver_count: int = None,
             solver_origins: List[str] = None) -> None:
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
            lines.append(f"// {name}: "
                         + ", ".join(decomposition.coefficients) + "\n")
            lines.append(f"inline constexpr std::size_t {name}NumCoefficients = "
                         f"{len(decomposition.coefficients)};\n")
            lines += _table(
                "model::CoefficientEntry",
                f"{name}CoefficientEntries",
                [f"{{{e.coefficient}, {e.dim}, {e.row}, {e.column}, {_format(e.factor)}}}"
                 for e in decomposition.entries],
            )
        lines += _table(
            "model::AnelasticCoefficientEntry",
            f"{name}AnelasticEntries",
            [f"{{{e.dim}, {e.row}, {e.column_offset}, {_format(e.factor)}}}"
             for e in decomposition.anelastic],
        )
        if decomposition.source_coefficients:
            lines.append(f"// {name} source: "
                         + ", ".join(decomposition.source_coefficients) + "\n")
            lines.append(f"inline constexpr std::size_t {name}NumSourceCoefficients = "
                         f"{len(decomposition.source_coefficients)};\n")
            lines += _table(
                "model::SourceCoefficientEntry",
                f"{name}SourceEntries",
                [f"{{{e.coefficient}, {e.row}, {e.column}, {_format(e.factor)}}}"
                 for e in decomposition.source],
            )

    if solver_count is not None:
        lines += [
            "// the operator the configured solver applies, material and any\n",
            "// relaxation blocks together\n",
            f"inline constexpr std::size_t SolverNumCoefficients = {solver_count};\n",
            "\n",
        ]
    if solver_count is not None:
        lines += [
            "// where each of them comes from, which decides whether a cell has\n",
            "// to carry it. Empty where the build does not factor the star, and\n",
            "// so carries no coefficients at all.\n",
        ]
        rows = [f"model::CoefficientOrigin::{origin}" for origin in (solver_origins or [])]
        if rows:
            lines += _table("model::CoefficientOrigin", "SolverCoefficientOrigins", rows)
        else:
            # the name has to exist even then: it is looked up unconditionally
            lines.append("inline constexpr std::array<model::CoefficientOrigin, 0> "
                         "SolverCoefficientOrigins{};\n\n")

    lines += [
        "} // namespace seissol::generated\n",
        "\n",
        "#endif // SEISSOL_GENERATEDCODE_COEFFICIENTS_H_\n",
    ]

    with open(path, "w") as file:
        file.writelines(lines)
