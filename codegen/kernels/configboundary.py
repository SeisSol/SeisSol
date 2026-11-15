# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""The kernels for the faces between cells of different configurations.

The configurations whose materials pose the Riemann problem in the same
quantities, and that fuse the same number of simulations, form a family (cf.
model::SameFamily). An anisotropic solid poses its Riemann problem in a material
of its own, but in the quantities of an elastic one: it is of the family of the
isotropic solids. The cells of two configurations can be face neighbors if they
fuse the same number of simulations. The
neighbor flux of a cell reads the time integral of a neighbor only through its
trace on the shared face, tested with the face basis of the cell. The time
integral of a neighbor of another configuration of the family is converted in
two steps, each one generated for a single configuration:

- toCanonical: from the configuration into the canonical form of its family:
  the largest order of the family in the build, the quantities of the Riemann
  problem, and double precision if a configuration of the build computes in
  double precision, single precision otherwise. The modal basis is
  hierarchical, so a lower order is padded with zeros. A configuration in
  single precision pads and selects in single precision and widens last: both
  only copy values, which is exact in any precision, and products in double
  precision run at a fraction of the rate on many GPUs.
- fromCanonical(side): from the canonical form into the configuration, for a
  neighbor that touches the cell with its side `side`. For a configuration of a
  lower order than the canonical one, the trace on that side is projected to the
  face basis of the configuration and lifted back into its volume basis, with
  the least norm; the result agrees with the neighbor on that face only. The
  quantities outside the canonical form, i.e. memory variables, are zero. For a
  configuration in single precision and a canonical form in double precision,
  the canonical form is narrowed first, and lifted and selected in single
  precision: TensorForge computes a kernel in the precision of its result,
  intermediate results included.

Hence the number of kernels grows linearly with the number of configurations,
and the kernels of a configuration only depend on the canonical form of its
family.

The families (FAMILIES: the solids, the fluids and the poroelastic solids) pose
the Riemann problem in different quantities, but the cells of any two of them
can be face neighbors as well: their families are coupled. A cell of the one
converts the time integral of a neighbor of the other from the canonical form
of the family of the neighbor:

- fromCoupledCanonical(family, side): as fromCanonical, from the canonical form
  of the coupled family with the index `family`, with the quantities mapped
  across by the matrix of the face (the runtime tensor coupling(family)): the
  contribution of the neighbor to the Godunov state at the face, in the
  quantities of the cell, which the code sets up from the impedances of both
  sides (model::coupledNeighborMap). The flux matrix of the cell at such a face
  is then its own flux alone. If the order of the configuration is the larger
  one, the basis is padded instead of lifted. As in fromCanonical, a
  configuration in single precision narrows first.

For that, toCanonical exists for the configurations of a family that another
family in the build is coupled to, too.
"""

import os
from typing import Dict, NamedTuple, Optional

import numpy as np
from kernels.common import generate_kernel_name_prefix
from kernels.multsim import OptionalDimTensor
from kernels.output import write_if_changed
from kernels.quantities import (
    FaceRole,
    QuantityGroup,
    QuantityKind,
    layout,
    total_extent,
)
from yateto import Tensor, parameterSpaceFromRanges, simpleParameterSpace
from yateto.functions import cast
from yateto.input import parseXMLMatrixFile
from yateto.memory import CSCMemoryLayout
from yateto.type import Datatype

#: The material whose quantities the Riemann problem of each equation has (cf.
#: RiemannQuantitiesOf in the C++ material models): an anisotropic solid has
#: those of an elastic one.
RIEMANN_MATERIAL = {
    "elastic": "elastic",
    "viscoelastic": "elastic",
    "acoustic": "acoustic",
    "viscoacoustic": "acoustic",
    "anisotropic": "elastic",
    "poroelastic": "poroelastic",
}

#: The families, by the materials whose quantities their Riemann problems have:
#: the solids, the fluids and the poroelastic solids. The index of a family is
#: its position here, as model::FamilyIndex has it. The cells of any two
#: families can be face neighbors: their families are coupled.
FAMILIES = ("elastic", "acoustic", "poroelastic")

#: The quantities of the Riemann problem of each family, as its equations lay
#: them out.
RIEMANN_GROUPS = {
    "elastic": (
        QuantityGroup("s", QuantityKind.SYM_TENSOR2, FaceRole.TRACTION),
        QuantityGroup("v", QuantityKind.VECTOR, FaceRole.VELOCITY),
    ),
    "acoustic": (
        QuantityGroup("pprime", QuantityKind.SCALAR, FaceRole.TRACTION),
        QuantityGroup("v", QuantityKind.VECTOR, FaceRole.VELOCITY),
    ),
    "poroelastic": (
        QuantityGroup("s", QuantityKind.SYM_TENSOR2, FaceRole.TRACTION),
        QuantityGroup("v", QuantityKind.VECTOR, FaceRole.VELOCITY),
        QuantityGroup("p", QuantityKind.SCALAR, FaceRole.EXTRA_TRACTION),
        QuantityGroup("vf", QuantityKind.VECTOR, FaceRole.EXTRA_VELOCITY),
    ),
}


def family(args):
    """What the configurations whose cells can be face neighbors share."""
    return (RIEMANN_MATERIAL[args.equations], args.multipleSimulations)


class Plan(NamedTuple):
    """The kernels of a configuration for its faces towards others."""

    #: the order of the canonical form of the family, which toCanonical
    #: converts to; None if toCanonical is not generated
    canonical_order: Optional[int]
    #: whether the canonical forms are in double precision: if a configuration
    #: of the build computes in double precision (the same for all families, so
    #: that the C++ side has one type for them)
    canonical_double: bool
    #: whether fromCanonical, from the canonical form of the family, is
    #: generated
    from_family: bool
    #: per coupled family with a configuration in the build, by its name: the
    #: order of its canonical form, which fromCoupledCanonical converts from
    coupled: Dict[str, int]


def is_double(args):
    """Whether the configuration of `args` computes in double precision."""
    return args.precision in ["d", "f64"]


def plans(config_args):
    """Per configuration: its Plan. The canonical form of a family has the
    largest order of the family in the build, and double precision if any
    configuration of the build has it. toCanonical and fromCanonical exist for
    the configurations of a family with more than one configuration in the
    build; toCanonical and fromCoupledCanonical for those of a family that
    another family in the build is coupled to, for each of them."""
    canonical_double = any(is_double(args) for args in config_args)
    orders = {}
    counts = {}
    for args in config_args:
        key = family(args)
        orders[key] = max(orders.get(key, 0), args.order)
        counts[key] = counts.get(key, 0) + 1

    result = []
    for args in config_args:
        material, simulations = family(args)
        coupled = {
            other: orders[(other, simulations)]
            for other in FAMILIES
            if other != material and (other, simulations) in orders
        }
        with_family = counts[(material, simulations)] > 1
        result.append(
            Plan(
                canonical_order=(
                    orders[(material, simulations)] if with_family or coupled else None
                ),
                canonical_double=canonical_double,
                from_family=with_family,
                coupled=coupled,
            )
        )
    return result


def num_bases(order):
    return order * (order + 1) * (order + 2) // 6


def num_face_bases(order):
    return order * (order + 1) // 2


def _dense(tensor):
    matrix = np.zeros(tensor.shape())
    for index, value in tensor.values().items():
        matrix[index] = float(value)
    return matrix


def face_traces(matrices_dir, order):
    """Per side of a cell: the trace of its volume basis on that side, in the
    face basis as a face neighbor sees it, i.e. fP(0) rT(side), as
    (face basis) x (volume basis)."""
    db = parseXMLMatrixFile(f"{matrices_dir}/aderdg-{order}.xml")
    orientation = _dense(db.fP[0])
    return [orientation @ _dense(db.rT[side]) for side in range(4)]


def lifts(matrices_dir, order, canonical_order):
    """Per side of the neighbor: the map from the volume basis of
    `canonical_order` to the one of `order` that keeps the trace on that side,
    projected to the face basis of `order`, with the least norm."""
    traces = face_traces(matrices_dir, order)
    canonical_traces = face_traces(matrices_dir, canonical_order)
    face_bases = num_face_bases(order)
    result = []
    for trace, canonical_trace in zip(traces, canonical_traces):
        right_inverse = trace.T @ np.linalg.inv(trace @ trace.T)
        lift = right_inverse @ canonical_trace[:face_bases, :]
        lift[np.abs(lift) < 1e-14] = 0.0
        result.append(lift)
    return result


def canonical_quantities(aderdg):
    """The quantities of the Riemann problem: the primary ones, which every
    configuration of a family has first, in the same order."""
    return total_extent(layout(aderdg.primaryGroups()))


def _fused(integral, name, shape, datatype, **kwargs):
    """A tensor with the fused simulations of `integral`, of `datatype`."""
    return OptionalDimTensor(
        name,
        integral.optName(),
        integral.optSize(),
        integral.optPos(),
        shape,
        datatype=datatype,
        **kwargs,
    )


def add_kernels(generator, aderdg, matrices_dir, plan, material, precision, targets):
    """Adds the kernels that `plan` names (see Plan) of the configuration of
    `aderdg`, whose equations pose the Riemann problem in `material`, for
    `targets`."""
    if plan.canonical_order is not None:
        _add_family_kernels(
            generator,
            aderdg,
            matrices_dir,
            plan.canonical_order,
            plan.from_family,
            canonical_type(plan),
            precision,
            targets,
        )
    if plan.coupled:
        _add_coupled_kernels(
            generator,
            aderdg,
            matrices_dir,
            plan.coupled,
            material,
            canonical_type(plan),
            precision,
            targets,
        )


def canonical_type(plan):
    """The datatype of the canonical forms (see Plan)."""
    return Datatype.F64 if plan.canonical_double else Datatype.F32


def _add_family_kernels(
    generator,
    aderdg,
    matrices_dir,
    canonical_order,
    from_family,
    wide,
    precision,
    targets,
):
    """Adds toCanonical and, if `from_family`, the family fromCanonical (over
    the side of the neighbor), with the canonical form of datatype `wide`."""
    real = Datatype.F64 if precision == "double" else Datatype.F32
    assert real == wide or wide == Datatype.F64

    integral = aderdg.I
    bases = aderdg.num3DBasisFunctions()
    quantities = aderdg.numQuantities()
    canonical_bases = num_bases(canonical_order)
    canonical_count = canonical_quantities(aderdg)
    assert bases <= canonical_bases and canonical_count <= quantities

    def fused(name, shape, **kwargs):
        return _fused(integral, name, shape, **kwargs)

    # unpadded, so that all configurations of a family lay it out alike
    canonical = fused("canonicalI", (canonical_bases, canonical_count), datatype=wide)

    padded = bases < canonical_bases
    selected = canonical_count < quantities

    # in the precision of the configuration, as toCanonical and fromCanonical
    # compute
    embed = Tensor(
        "canonicalEmbed",
        (canonical_bases, bases),
        np.eye(canonical_bases, bases),
        CSCMemoryLayout,
        datatype=real,
    )
    select = np.eye(quantities, canonical_count)
    selectTo = Tensor(
        "canonicalSelect", select.shape, select, CSCMemoryLayout, datatype=real
    )
    selectFrom = Tensor(
        "canonicalSelectT", select.T.shape, select.T, CSCMemoryLayout, datatype=real
    )
    lift = None
    if padded and from_family:
        lift = [
            Tensor(
                f"canonicalLift({side})",
                (bases, canonical_bases),
                values,
                datatype=real,
            )
            for side, values in enumerate(
                lifts(matrices_dir, aderdg.order, canonical_order)
            )
        ]

    for target in targets:
        prefix = generate_kernel_name_prefix(target)

        # to the canonical form: pad and select first, then widen
        if not padded and not selected:
            source = integral["kq"]
            toCanonical = canonical["kq"] <= (
                cast(source, wide) if real != wide else source
            )
        else:
            if padded and selected:
                product = embed["kl"] * integral["lp"] * selectTo["pq"]
            elif padded:
                product = embed["kl"] * integral["lq"]
            else:
                product = integral["kp"] * selectTo["pq"]
            if real != wide:
                narrow = fused(
                    f"{prefix}narrowCanonicalI",
                    (canonical_bases, canonical_count),
                    datatype=real,
                    temporary=True,
                )
                toCanonical = [
                    narrow["kq"] <= product,
                    canonical["kq"] <= cast(narrow["kq"], wide),
                ]
            else:
                toCanonical = canonical["kq"] <= product
        generator.add(f"{prefix}toCanonical", toCanonical, target=target)

        # from the canonical form: narrow first, then lift and select
        def fromCanonical(side):
            if not padded and not selected:
                source = canonical["kp"]
                return integral["kp"] <= (
                    cast(source, real) if real != wide else source
                )
            if real != wide:
                narrow = fused(
                    f"{prefix}narrowCanonicalI",
                    (canonical_bases, canonical_count),
                    datatype=real,
                    temporary=True,
                )
                narrowing = [narrow["lq"] <= cast(canonical["lq"], real)]
            else:
                narrow = canonical
                narrowing = []
            if padded and selected:
                product = lift[side]["kl"] * narrow["lq"] * selectFrom["qp"]
            elif padded:
                product = lift[side]["kl"] * narrow["lp"]
            else:
                product = narrow["kq"] * selectFrom["qp"]
            return narrowing + [integral["kp"] <= product]

        if from_family:
            generator.addFamily(
                f"{prefix}fromCanonical",
                simpleParameterSpace(4),
                fromCanonical,
                target=target,
            )


def _add_coupled_kernels(
    generator, aderdg, matrices_dir, coupled, material, wide, precision, targets
):
    """Adds the family fromCoupledCanonical, over the index of the coupled family
    of the neighbor (see FAMILIES) and the side of the neighbor, for the coupled
    families in `coupled` (see Plan), whose canonical forms are of datatype
    `wide`. Its tensors are groups over the index of the coupled family."""
    real = Datatype.F64 if precision == "double" else Datatype.F32
    assert real == wide or wide == Datatype.F64

    integral = aderdg.I
    order = aderdg.order
    bases = aderdg.num3DBasisFunctions()
    quantities = aderdg.numQuantities()

    assert tuple(aderdg.primaryGroups()) == RIEMANN_GROUPS[material]

    shape = {}
    canonical = {}
    coupling = {}
    lift = {}
    for name, coupled_order in coupled.items():
        index = FAMILIES.index(name)
        coupled_count = total_extent(layout(RIEMANN_GROUPS[name]))
        coupled_bases = num_bases(coupled_order)

        shape[index] = (coupled_bases, coupled_count)
        canonical[index] = _fused(
            integral, f"coupledCanonicalI({index})", shape[index], datatype=wide
        )

        # per face: the quantities of the neighbor into the contribution of the
        # neighbor to the Godunov state, in the quantities of the cell and in the
        # precision of the configuration
        coupling[index] = Tensor(
            f"coupling({index})", (coupled_count, quantities), datatype=real
        )

        # the basis: lifted from a larger order as in fromCanonical, padded from
        # a smaller one
        if coupled_order > order:
            lift[index] = [
                Tensor(
                    f"coupledLift({index},{side})",
                    (bases, coupled_bases),
                    values,
                    datatype=real,
                )
                for side, values in enumerate(lifts(matrices_dir, order, coupled_order))
            ]
        elif coupled_order < order:
            embed = Tensor(
                f"coupledEmbed({index})",
                (bases, coupled_bases),
                np.eye(bases, coupled_bases),
                CSCMemoryLayout,
                datatype=real,
            )
            lift[index] = [embed] * 4
        else:
            lift[index] = None

    for target in targets:
        prefix = generate_kernel_name_prefix(target)

        def fromCoupledCanonical(index, side):
            source = canonical[index]
            if real != wide:
                narrow = _fused(
                    integral,
                    f"{prefix}narrowCoupledCanonicalI({index})",
                    shape[index],
                    datatype=real,
                    temporary=True,
                )
                narrowing = [narrow["lq"] <= cast(source["lq"], real)]
            else:
                narrow = source
                narrowing = []

            if lift[index] is None:
                product = narrow["kq"] * coupling[index]["qp"]
            else:
                product = lift[index][side]["kl"] * narrow["lq"] * coupling[index]["qp"]
            return narrowing + [integral["kp"] <= product]

        generator.addFamily(
            f"{prefix}fromCoupledCanonical",
            parameterSpaceFromRanges(sorted(canonical), range(4)),
            fromCoupledCanonical,
            target=target,
        )


def emit_header(aderdg, output_dir, key, plan, targets):
    """Tells the C++ side which kernels of the configuration `key` exist, and
    for which canonical forms."""
    to_canonical = plan.canonical_order is not None
    count = canonical_quantities(aderdg) if to_canonical else 0

    def boolean(value):
        return "true" if value else "false"

    def per_family(value):
        return "{" + ", ".join(value(name) for name in FAMILIES) + "}"

    def coupled_kernels(target):
        return per_family(
            lambda name: boolean(name in plan.coupled and target in targets)
        )

    coupled_orders = per_family(lambda name: str(plan.coupled.get(name, 0)))
    coupled_counts = per_family(
        lambda name: (
            str(total_extent(layout(RIEMANN_GROUPS[name])))
            if name in plan.coupled
            else "0"
        )
    )

    lines = [
        "// Generated by codegen/kernels/configboundary.py. Do not edit.",
        "",
        "#pragma once",
        "",
        '#include "Config.h"',
        "",
        "#include <array>",
        "#include <cstddef>",
        "",
        "namespace seissol::generated {",
        "",
        "template <typename Cfg>",
        "struct ConfigBoundaryKernels;",
        "",
        "/// Which kernels of the faces towards other configurations exist on the host and on the device,",
        "/// and the canonical forms they convert to and from: their order, their number of quantities",
        "/// and their precision.",
        "template <>",
        f"struct ConfigBoundaryKernels<{key}> {{",
        "  // toCanonical and fromCanonical, into and from the canonical form of the family",
        f"  static constexpr bool Host = {boolean(plan.from_family and 'cpu' in targets)};",
        f"  static constexpr bool Device = {boolean(plan.from_family and 'gpu' in targets)};",
        "  // toCanonical alone, which the configurations of the coupled families read as well",
        f"  static constexpr bool ToCanonicalHost = {boolean(to_canonical and 'cpu' in targets)};",
        f"  static constexpr bool ToCanonicalDevice = {boolean(to_canonical and 'gpu' in targets)};",
        f"  static constexpr std::size_t CanonicalOrder = {plan.canonical_order or 0};",
        "  // the precision of the canonical forms, of all families: double if a configuration of the",
        "  // build computes in double precision",
        f"  using CanonicalReal = {'double' if plan.canonical_double else 'float'};",
        f"  static constexpr std::size_t CanonicalQuantities = {count};",
        "  // fromCoupledCanonical, by the family of the neighbor (the solids, the fluids, the",
        "  // poroelastic solids): whether it converts from the canonical form of that family, and the",
        "  // order and the quantities of the form",
        f"  static constexpr std::array<bool, 3> CoupledHost = {coupled_kernels('cpu')};",
        f"  static constexpr std::array<bool, 3> CoupledDevice = {coupled_kernels('gpu')};",
        f"  static constexpr std::array<std::size_t, 3> CoupledCanonicalOrder = {coupled_orders};",
        f"  static constexpr std::array<std::size_t, 3> CoupledCanonicalQuantities = {coupled_counts};",
        "};",
        "",
        "} // namespace seissol::generated",
        "",
    ]
    write_if_changed(os.path.join(output_dir, "configboundary.h"), "\n".join(lines))
