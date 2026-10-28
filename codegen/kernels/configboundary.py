# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""The kernels for the faces between cells of different configurations.

The cells of two configurations can be face neighbors if their materials pose
the Riemann problem in the same material and if they fuse the same number of
simulations (cf. model::CanNeighbor); such configurations form a family. The
neighbor flux of a cell reads the time integral of a neighbor only through its
trace on the shared face, tested with the face basis of the cell. The time
integral of a neighbor of another configuration of the family is converted in
two steps, each one generated for a single configuration:

- toCanonical: from the configuration into the canonical form of its family:
  the largest order of the family in the build, the quantities of the Riemann
  problem, double precision. The modal basis is hierarchical, so a lower order
  is padded with zeros.
- fromCanonical(side): from the canonical form into the configuration, for a
  neighbor that touches the cell with its side `side`. For a configuration of a
  lower order than the canonical one, the trace on that side is projected to the
  face basis of the configuration and lifted back into its volume basis, with
  the least norm; the result agrees with the neighbor on that face only. The
  quantities outside the canonical form, i.e. memory variables, are zero. For a
  configuration in single precision, the canonical form is narrowed first, and
  lifted and selected in single precision: TensorForge computes a kernel in the
  precision of its result, intermediate results included.

Hence the number of kernels grows linearly with the number of configurations,
and the kernels of a configuration only depend on the canonical form of its
family.

A solid and a fluid pose the Riemann problem in different materials, but their
cells can be face neighbors as well: their families are coupled. A cell of the
one converts the time integral of a neighbor of the other from the canonical
form of the family of the neighbor:

- fromCoupledCanonical(side): as fromCanonical, from the canonical form of the
  coupled family, with the quantities mapped across. A solid takes the pressure
  of a fluid as an isotropic stress. A fluid takes the normal stress of a solid
  on the shared face as its pressure, with weights that the face gives (the
  runtime tensor normalStress). The velocities carry over. If the order of the
  configuration is the larger one, the basis is padded instead of lifted. As
  in fromCanonical, a configuration in single precision narrows first.

For that, toCanonical exists for the configurations of a family whose coupled
family is in the build, too.
"""

import os
from typing import NamedTuple, Optional

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
from yateto import Tensor, simpleParameterSpace
from yateto.functions import cast
from yateto.input import parseXMLMatrixFile
from yateto.memory import CSCMemoryLayout
from yateto.type import Datatype

#: The material each equation poses the Riemann problem in (cf. RiemannMaterial
#: in the C++ material models).
RIEMANN_MATERIAL = {
    "elastic": "elastic",
    "viscoelastic": "elastic",
    "acoustic": "acoustic",
    "viscoacoustic": "acoustic",
    "anisotropic": "anisotropic",
    "poroelastic": "poroelastic",
}


#: The Riemann materials whose cells are face neighbors although they pose the
#: Riemann problem in different materials: a solid and a fluid (cf.
#: model::CanNeighbor).
COUPLED = {"elastic": "acoustic", "acoustic": "elastic"}

#: The quantities of the Riemann problem of the coupled Riemann materials, as
#: their equations lay them out.
RIEMANN_GROUPS = {
    "elastic": (
        QuantityGroup("s", QuantityKind.SYM_TENSOR2, FaceRole.TRACTION),
        QuantityGroup("v", QuantityKind.VECTOR, FaceRole.VELOCITY),
    ),
    "acoustic": (
        QuantityGroup("pprime", QuantityKind.SCALAR, FaceRole.TRACTION),
        QuantityGroup("v", QuantityKind.VECTOR, FaceRole.VELOCITY),
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
    #: whether fromCanonical, from the canonical form of the family, is
    #: generated
    from_family: bool
    #: the order of the canonical form of the coupled family, which
    #: fromCoupledCanonical converts from; None if that is not generated
    coupled_order: Optional[int]


def plans(config_args):
    """Per configuration: its Plan. The canonical form of a family has the
    largest order of the family in the build. toCanonical and fromCanonical
    exist for the configurations of a family with more than one configuration in
    the build; toCanonical and fromCoupledCanonical for those of a family whose
    coupled family is in the build as well."""
    orders = {}
    counts = {}
    for args in config_args:
        key = family(args)
        orders[key] = max(orders.get(key, 0), args.order)
        counts[key] = counts.get(key, 0) + 1

    result = []
    for args in config_args:
        material, simulations = family(args)
        coupled = (COUPLED.get(material), simulations)
        with_family = counts[(material, simulations)] > 1
        with_coupled = coupled in counts
        result.append(
            Plan(
                canonical_order=(
                    orders[(material, simulations)]
                    if with_family or with_coupled
                    else None
                ),
                from_family=with_family,
                coupled_order=orders[coupled] if with_coupled else None,
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


def _fused(integral, name, shape, datatype=Datatype.F64, **kwargs):
    """A tensor with the fused simulations of `integral`, in double precision
    unless `datatype` says otherwise."""
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
            precision,
            targets,
        )
    if plan.coupled_order is not None:
        _add_coupled_kernels(
            generator,
            aderdg,
            matrices_dir,
            plan.coupled_order,
            material,
            precision,
            targets,
        )


def _add_family_kernels(
    generator, aderdg, matrices_dir, canonical_order, from_family, precision, targets
):
    """Adds toCanonical and, if `from_family`, the family fromCanonical (over
    the side of the neighbor)."""
    real = Datatype.F64 if precision == "double" else Datatype.F32
    f64 = Datatype.F64

    integral = aderdg.I
    bases = aderdg.num3DBasisFunctions()
    quantities = aderdg.numQuantities()
    canonical_bases = num_bases(canonical_order)
    canonical_count = canonical_quantities(aderdg)
    assert bases <= canonical_bases and canonical_count <= quantities

    def fused(name, shape, **kwargs):
        return _fused(integral, name, shape, **kwargs)

    # unpadded, so that all configurations of a family lay it out alike
    canonical = fused("canonicalI", (canonical_bases, canonical_count))

    padded = bases < canonical_bases
    selected = canonical_count < quantities

    embed = Tensor(
        "canonicalEmbed",
        (canonical_bases, bases),
        np.eye(canonical_bases, bases),
        CSCMemoryLayout,
        datatype=f64,
    )
    select = np.eye(quantities, canonical_count)
    selectTo = Tensor(
        "canonicalSelect", select.shape, select, CSCMemoryLayout, datatype=f64
    )
    # in the precision of the configuration, as fromCanonical computes
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

        # to the canonical form: widen first, then pad and select
        if not padded and not selected:
            source = integral["kq"]
            toCanonical = canonical["kq"] <= (
                cast(source, f64) if real != f64 else source
            )
        else:
            if real != f64:
                wide = fused(f"{prefix}wideI", (bases, quantities), temporary=True)
                widen = [wide["kp"] <= cast(integral["kp"], f64)]
            else:
                wide = integral
                widen = []
            if padded and selected:
                product = embed["kl"] * wide["lp"] * selectTo["pq"]
            elif padded:
                product = embed["kl"] * wide["lq"]
            else:
                product = wide["kp"] * selectTo["pq"]
            toCanonical = widen + [canonical["kq"] <= product]
        generator.add(f"{prefix}toCanonical", toCanonical, target=target)

        # from the canonical form: narrow first, then lift and select
        def fromCanonical(side):
            if not padded and not selected:
                source = canonical["kp"]
                return integral["kp"] <= (cast(source, real) if real != f64 else source)
            if real != f64:
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
    generator, aderdg, matrices_dir, coupled_order, material, precision, targets
):
    """Adds the family fromCoupledCanonical (over the side of the neighbor)."""
    real = Datatype.F64 if precision == "double" else Datatype.F32
    f64 = Datatype.F64

    integral = aderdg.I
    order = aderdg.order
    bases = aderdg.num3DBasisFunctions()
    quantities = aderdg.numQuantities()

    groups = tuple(aderdg.primaryGroups())
    assert groups == RIEMANN_GROUPS[material]
    blocks = layout(groups)
    coupled_blocks = layout(RIEMANN_GROUPS[COUPLED[material]])
    coupled_count = total_extent(coupled_blocks)
    coupled_bases = num_bases(coupled_order)

    coupled = _fused(integral, "coupledCanonicalI", (coupled_bases, coupled_count))

    # the quantities: what carries over as it is, and the traction of a fluid as
    # the isotropic stress of a solid; the traction of a solid becomes the
    # pressure of a fluid with the weights of the face
    mapping = np.zeros((coupled_count, quantities))
    pressure = None
    for block in blocks:
        (source,) = [
            other for other in coupled_blocks if other.group.role is block.group.role
        ]
        if block.group.kind is source.group.kind:
            for component in range(block.extent):
                mapping[source.offset + component, block.offset + component] = 1
        elif source.group.kind is QuantityKind.SCALAR:
            assert block.group.kind is QuantityKind.SYM_TENSOR2
            # the normal components, xx, yy and zz
            for component in range(3):
                mapping[source.offset, block.offset + component] = 1
        else:
            assert source.group.kind is QuantityKind.SYM_TENSOR2
            assert block.group.kind is QuantityKind.SCALAR
            pressure = np.zeros((1, quantities))
            pressure[0, block.offset] = 1

    # in the precision of the configuration, as fromCoupledCanonical computes
    mappingTensor = Tensor(
        "coupledMap", mapping.shape, mapping, CSCMemoryLayout, datatype=real
    )
    weights = None
    pressureTensor = None
    if pressure is not None:
        # per face: the weights of the stress components in the normal stress
        weights = Tensor("normalStress", (coupled_count, 1), datatype=real)
        pressureTensor = Tensor(
            "coupledPressure",
            pressure.shape,
            pressure,
            CSCMemoryLayout,
            datatype=real,
        )

    # the basis: lifted from a larger order as in fromCanonical, padded from a
    # smaller one
    lift = None
    if coupled_order > order:
        lift = [
            Tensor(
                f"coupledLift({side})",
                (bases, coupled_bases),
                values,
                datatype=real,
            )
            for side, values in enumerate(lifts(matrices_dir, order, coupled_order))
        ]
    elif coupled_order < order:
        embed = Tensor(
            "coupledEmbed",
            (bases, coupled_bases),
            np.eye(bases, coupled_bases),
            CSCMemoryLayout,
            datatype=real,
        )
        lift = [embed] * 4

    for target in targets:
        prefix = generate_kernel_name_prefix(target)

        def fromCoupledCanonical(side):
            if real != f64:
                narrow = _fused(
                    integral,
                    f"{prefix}narrowCoupledCanonicalI",
                    (coupled_bases, coupled_count),
                    datatype=real,
                    temporary=True,
                )
                narrowing = [narrow["lq"] <= cast(coupled["lq"], real)]
            else:
                narrow = coupled
                narrowing = []

            def basis(right):
                if lift is None:
                    return narrow["kq"] * right
                return lift[side]["kl"] * narrow["lq"] * right

            product = basis(mappingTensor["qp"])
            if weights is not None:
                product = product + basis(weights["qa"] * pressureTensor["ap"])
            return narrowing + [integral["kp"] <= product]

        generator.addFamily(
            f"{prefix}fromCoupledCanonical",
            simpleParameterSpace(4),
            fromCoupledCanonical,
            target=target,
        )


def emit_header(aderdg, output_dir, key, plan, material, targets):
    """Tells the C++ side which kernels of the configuration `key` exist, and
    for which canonical forms."""
    to_canonical = plan.canonical_order is not None
    coupled = plan.coupled_order is not None
    count = canonical_quantities(aderdg) if to_canonical else 0
    coupled_count = (
        total_extent(layout(RIEMANN_GROUPS[COUPLED[material]])) if coupled else 0
    )
    # a fluid weighs the stress of a solid with the normal of the face
    normal_stress = (
        coupled
        and RIEMANN_GROUPS[COUPLED[material]][0].kind is QuantityKind.SYM_TENSOR2
    )

    def boolean(value):
        return "true" if value else "false"

    lines = [
        "// Generated by codegen/kernels/configboundary.py. Do not edit.",
        "",
        "#pragma once",
        "",
        '#include "Config.h"',
        "",
        "#include <cstddef>",
        "",
        "namespace seissol::generated {",
        "",
        "template <typename Cfg>",
        "struct ConfigBoundaryKernels;",
        "",
        "/// Which kernels of the faces towards other configurations exist on the host and on the device,",
        "/// and the canonical forms they convert to and from: their order and their number of quantities.",
        "template <>",
        f"struct ConfigBoundaryKernels<{key}> {{",
        "  // toCanonical and fromCanonical, into and from the canonical form of the family",
        f"  static constexpr bool Host = {boolean(plan.from_family and 'cpu' in targets)};",
        f"  static constexpr bool Device = {boolean(plan.from_family and 'gpu' in targets)};",
        "  // toCanonical alone, which the configurations of the coupled family read as well",
        f"  static constexpr bool ToCanonicalHost = {boolean(to_canonical and 'cpu' in targets)};",
        f"  static constexpr bool ToCanonicalDevice = {boolean(to_canonical and 'gpu' in targets)};",
        f"  static constexpr std::size_t CanonicalOrder = {plan.canonical_order or 0};",
        f"  static constexpr std::size_t CanonicalQuantities = {count};",
        "  // fromCoupledCanonical, from the canonical form of the coupled family (a solid and a fluid)",
        f"  static constexpr bool CoupledHost = {boolean(coupled and 'cpu' in targets)};",
        f"  static constexpr bool CoupledDevice = {boolean(coupled and 'gpu' in targets)};",
        f"  static constexpr std::size_t CoupledCanonicalOrder = {plan.coupled_order or 0};",
        f"  static constexpr std::size_t CoupledCanonicalQuantities = {coupled_count};",
        "  // whether fromCoupledCanonical reads the weights of the stress components in the normal",
        "  // stress of the face (normalStress)",
        f"  static constexpr bool NormalStress = {boolean(normal_stress)};",
        "};",
        "",
        "} // namespace seissol::generated",
        "",
    ]
    write_if_changed(os.path.join(output_dir, "configboundary.h"), "\n".join(lines))
