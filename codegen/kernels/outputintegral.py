# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""The time integral since the start that the outputs read (LTS::Integrals)."""

from kernels.common import generate_kernel_name_prefix
from kernels.multsim import OptionalDimTensor
from yateto.functions import cast
from yateto.type import Datatype

# the kernels, which the configurations of a build may not all have: declared for all of them
KERNELS = ["integrateOutput", "gpu_integrateOutput"]


def addKernels(generator, aderdg, precision, targets):
    """Adds integrateOutput, for `targets`: the time integral of the predictor over a step
    (aderdg.I) added to the time integral since the start, in binary64 in every configuration. In
    binary32, the integral of a step falls below the rounding of a sum that has grown over many
    steps. Laid out as the unknowns: the alignment of a layout counts elements of the precision
    of the configuration, whatever the type of the tensor."""
    integral = OptionalDimTensor(
        "QIntegral",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (aderdg.num3DBasisFunctions(), aderdg.numQuantities()),
        alignStride=True,
        datatype=Datatype.F64,
    )
    step = aderdg.I["kp"]
    widened = step if precision == "double" else cast(step, Datatype.F64)
    for target in targets:
        prefix = generate_kernel_name_prefix(target)
        generator.add(
            f"{prefix}integrateOutput",
            integral["kp"] <= integral["kp"] + widened,
            target=target,
        )
