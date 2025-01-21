# SPDX-FileCopyrightText: 2017-2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff

from kernels.common import generate_kernel_name_prefix
from kernels.multsim import OptionalDimTensor
from yateto import Tensor, simpleParameterSpace


def addKernels(generator, aderdg, include_tensors, targets):
    maxDepth = 3

    numberOf3DBasisFunctions = aderdg.numberOf3DBasisFunctions()
    numberOf2DBasisFunctions = aderdg.numberOf2DBasisFunctions()

    faceDisplacement = OptionalDimTensor(
        "faceDisplacement",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (numberOf2DBasisFunctions, 3),
        alignStride=True,
    )
    averageNormalDisplacement = OptionalDimTensor(
        "averageNormalDisplacement",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (numberOf2DBasisFunctions,),
        alignStride=True,
    )

    include_tensors.add(averageNormalDisplacement)

    subTriangleDofs = [
        OptionalDimTensor(
            "subTriangleDofs({})".format(depth),
            aderdg.Q.optName(),
            aderdg.Q.optSize(),
            aderdg.Q.optPos(),
            (4**depth, 3),
            alignStride=True,
        )
        for depth in range(maxDepth + 1)
    ]
    subTriangleProjection = [
        Tensor(
            "subTriangleProjection({})".format(depth),
            (4**depth, numberOf3DBasisFunctions),
            alignStride=True,
        )
        for depth in range(maxDepth + 1)
    ]
    subTriangleProjectionFromFace = [
        Tensor(
            "subTriangleProjectionFromFace({})".format(depth),
            (4**depth, numberOf2DBasisFunctions),
            alignStride=True,
        )
        for depth in range(maxDepth + 1)
    ]

    displacementRotationMatrix = Tensor(
        "displacementRotationMatrix", (3, 3), alignStride=True
    )
    subTriangleDisplacement = (
        lambda depth: subTriangleDofs[depth]["kp"]
        <= subTriangleProjectionFromFace[depth]["kl"]
        * aderdg.db.MV2nTo2m["lm"]
        * faceDisplacement["mp"]
    )
    subTriangleVelocity = (
        lambda depth: subTriangleDofs[depth]["kp"]
        <= subTriangleProjection[depth]["kl"]
        * aderdg.Q["lq"]
        * aderdg.selectVelocity["qp"]
    )

    generator.addFamily(
        "subTriangleDisplacement",
        simpleParameterSpace(maxDepth + 1),
        subTriangleDisplacement,
    )
    generator.addFamily(
        "subTriangleVelocity",
        simpleParameterSpace(maxDepth + 1),
        subTriangleVelocity,
    )

    rotatedFaceDisplacement = OptionalDimTensor(
        "rotatedFaceDisplacement",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (numberOf2DBasisFunctions, 3),
        alignStride=True,
    )
    for target in targets:
        name_prefix = generate_kernel_name_prefix(target)
        generator.add(
            f"{name_prefix}rotateFaceDisplacement",
            rotatedFaceDisplacement["mp"]
            <= faceDisplacement["mn"] * displacementRotationMatrix["pn"],
            target=target,
        )

    addVelocity = (
        lambda f: faceDisplacement["kp"]
        <= faceDisplacement["kp"]
        + aderdg.db.V3mTo2nFace[f]["kl"] * aderdg.I["lq"] * aderdg.selectVelocity["qp"]
    )
    generator.addFamily("addVelocity", simpleParameterSpace(4), addVelocity)

    numberOfQuadratureNodes = (aderdg.order + 1) ** 2
    rotatedFaceDisplacementAtQuadratureNodes = OptionalDimTensor(
        "rotatedFaceDisplacementAtQuadratureNodes",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (numberOfQuadratureNodes, 3),
        alignStride=True,
    )
    generator.add(
        "rotateFaceDisplacementsAndEvaluateAtQuadratureNodes",
        rotatedFaceDisplacementAtQuadratureNodes["in"]
        <= aderdg.V2nTo2JacobiQuad["ij"]
        * rotatedFaceDisplacement["jp"]
        * displacementRotationMatrix["np"],
    )

    if "gpu" in targets:
        name_prefix = generate_kernel_name_prefix(target="gpu")

        integratedVelocities = OptionalDimTensor(
            "integratedVelocities",
            aderdg.I.optName(),
            aderdg.I.optSize(),
            aderdg.I.optPos(),
            (numberOf3DBasisFunctions, 3),
            alignStride=True,
        )

        addVelocity = (
            lambda f: faceDisplacement["kp"]
            <= faceDisplacement["kp"]
            + aderdg.db.V3mTo2nFace[f]["kl"] * integratedVelocities["lp"]
        )

        generator.addFamily(
            f"{name_prefix}addVelocity",
            simpleParameterSpace(4),
            addVelocity,
            target="gpu",
        )
