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
"""

import os

import numpy as np
from kernels.common import generate_kernel_name_prefix
from kernels.multsim import OptionalDimTensor
from kernels.quantities import layout, total_extent
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


def family(args):
    """What the configurations whose cells can be face neighbors share."""
    return (RIEMANN_MATERIAL[args.equations], args.multipleSimulations)


def canonical_orders(config_args):
    """Per configuration: the order of the canonical form of its family, i.e.
    the largest order of the family in the build; None for a configuration
    whose family has no other configuration in the build."""
    orders = {}
    counts = {}
    for args in config_args:
        key = family(args)
        orders[key] = max(orders.get(key, 0), args.order)
        counts[key] = counts.get(key, 0) + 1
    return [
        orders[family(args)] if counts[family(args)] > 1 else None
        for args in config_args
    ]


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


def add_kernels(generator, aderdg, matrices_dir, canonical_order, precision, targets):
    """Adds toCanonical and the family fromCanonical (over the side of the
    neighbor) of the configuration of `aderdg` for `targets`."""
    real = Datatype.F64 if precision == "double" else Datatype.F32
    f64 = Datatype.F64

    integral = aderdg.I
    bases = aderdg.num3DBasisFunctions()
    quantities = aderdg.numQuantities()
    canonical_bases = num_bases(canonical_order)
    canonical_count = canonical_quantities(aderdg)
    assert bases <= canonical_bases and canonical_count <= quantities

    def fused(name, shape, datatype=f64, **kwargs):
        return OptionalDimTensor(
            name,
            integral.optName(),
            integral.optSize(),
            integral.optPos(),
            shape,
            datatype=datatype,
            **kwargs,
        )

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
    if padded:
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

        generator.addFamily(
            f"{prefix}fromCanonical",
            simpleParameterSpace(4),
            fromCanonical,
            target=target,
        )


def emit_header(aderdg, output_dir, key, canonical_order, targets):
    """Tells the C++ side whether the kernels of the configuration `key` exist,
    and for which canonical form."""
    host = canonical_order is not None and "cpu" in targets
    device = canonical_order is not None and "gpu" in targets
    count = canonical_quantities(aderdg) if canonical_order is not None else 0

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
        "/// Whether the kernels toCanonical and fromCanonical exist on the host and on the device, and",
        "/// the canonical form they convert to and from: its order and its number of quantities.",
        "template <>",
        f"struct ConfigBoundaryKernels<{key}> {{",
        f"  static constexpr bool Host = {boolean(host)};",
        f"  static constexpr bool Device = {boolean(device)};",
        f"  static constexpr std::size_t CanonicalOrder = {canonical_order or 0};",
        f"  static constexpr std::size_t CanonicalQuantities = {count};",
        "};",
        "",
        "} // namespace seissol::generated",
        "",
    ]
    with open(os.path.join(output_dir, "configboundary.h"), "w") as file:
        file.write("\n".join(lines))
