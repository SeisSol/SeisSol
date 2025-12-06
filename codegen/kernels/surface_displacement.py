# SPDX-FileCopyrightText: 2017 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff

from kernels.common import generate_kernel_name_prefix
from kernels.multsim import OptionalDimTensor
from yateto import Scalar, Tensor, simpleParameterSpace


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
        + aderdg.db.V3mTo2nFace[f][aderdg.t("kl")]
        * aderdg.I["lq"]
        * aderdg.selectVelocity["qp"]
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
            + aderdg.db.V3mTo2nFace[f][aderdg.t("kl")] * integratedVelocities["lp"]
        )

        generator.addFamily(
            f"{name_prefix}addVelocity",
            simpleParameterSpace(4),
            addVelocity,
            target="gpu",
        )

    Iprev = OptionalDimTensor(
        "Iprev",
        aderdg.INodal.optName(),
        aderdg.INodal.optSize(),
        aderdg.INodal.optPos(),
        (aderdg.numberOf2DBasisFunctions(), 1),
        alignStride=True,
        temporary=True,
    )
    averageNormalDisplacement = OptionalDimTensor(
        "Iint",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (numberOf2DBasisFunctions, 1),
        alignStride=True,
    )
    faceDisplacementTmp = OptionalDimTensor(
        "IAcc",
        aderdg.INodal.optName(),
        aderdg.INodal.optSize(),
        aderdg.INodal.optPos(),
        (aderdg.numberOf2DBasisFunctions(), 3),
        alignStride=True,
        temporary=True,
    )
    INodalTmp = OptionalDimTensor(
        "INodalTmp",
        aderdg.INodal.optName(),
        aderdg.INodal.optSize(),
        aderdg.INodal.optPos(),
        (aderdg.numberOf2DBasisFunctions(), aderdg.numberOfQuantities()),
        alignStride=True,
        temporary=True,
    )

    coeffs = [Scalar(f"coeff({i})") for i in range(aderdg.order)]
    factor = Tensor("invImpFactor", ())

    vidx = aderdg.velocityOffset()

    for target in targets:
        name_prefix = generate_kernel_name_prefix(target)

        def kernelPerFace(f):
            kernel = [
                faceDisplacementTmp["mp"]
                <= faceDisplacement["mn"]
                * aderdg.Tinv["pn"].subslice("p", vidx, 9).subslice("n", vidx, 9),
                Iprev["mp"] <= faceDisplacementTmp["mp"].subslice("p", 0, 1),
                averageNormalDisplacement["mp"]
                <= aderdg.powers[0] * faceDisplacementTmp["mp"].subslice("p", 0, 1),
            ]

            for i in range(1, aderdg.order):

                kernel += [
                    INodalTmp["nq"]
                    <= aderdg.db.V3mTo2nFace[f]["nk"] * aderdg.dQs[i]["kq"]
                ]
                velocitiesVW = INodalTmp["nq"].subslice("q", vidx + 1, vidx + 3)
                velocitiesU = INodalTmp["nq"].subslice("q", vidx, vidx + 1)
                pressure = INodalTmp["nq"].subslice("q", 0, 1)

                kernel += [
                    Iprev["nq"] <= velocitiesU + factor[""] * Iprev["nq"] + pressure,
                    faceDisplacementTmp["nq"].subslice("q", 0, 1)
                    <= coeffs[i] * Iprev["nq"],
                    faceDisplacementTmp["nq"].subslice("q", 1, 3)
                    <= coeffs[i] * velocitiesVW,
                    averageNormalDisplacement["nq"]
                    <= averageNormalDisplacement["nq"] + aderdg.powers[i] * Iprev["nq"],
                ]

            kernel += [
                faceDisplacement["mp"]
                <= faceDisplacementTmp["mn"]
                * aderdg.T["pn"]
                .subslice("p", vidx, vidx + 3)
                .subslice("n", vidx, vidx + 3),
            ]

            return kernel

        generator.addFamily(
            f"{name_prefix}fsgKernel",
            simpleParameterSpace(4),
            kernelPerFace,
            target=target,
        )
