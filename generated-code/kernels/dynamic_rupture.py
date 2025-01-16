# SPDX-FileCopyrightText: 2016-2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff

from kernels.common import generate_kernel_name_prefix
from kernels.multsim import OptionalDimTensor
from yateto import Scalar, Tensor, simpleParameterSpace
from yateto.input import parseJSONMatrixFile


def addKernels(generator, aderdg, matricesDir, drQuadRule, targets):

    clones = dict()

    # Load matrices
    db = parseJSONMatrixFile(
        f"{matricesDir}/dr_{drQuadRule}_matrices_{aderdg.order}.json",
        clones,
        alignStride=aderdg.alignStride,
        transpose=aderdg.transpose,
    )
    numberOfPoints = db.resample.shape()[0]

    # Determine matrices
    # Note: This does only work because the flux does not depend
    # on the mechanisms in the case of viscoelastic attenuation
    trans_inv_spp_T = aderdg.transformation_inv_spp().transpose()
    TinvT = Tensor("TinvT", trans_inv_spp_T.shape, spp=trans_inv_spp_T)
    flux_solver_spp = aderdg.flux_solver_spp()
    fluxSolver = Tensor("fluxSolver", flux_solver_spp.shape, spp=flux_solver_spp)

    gShape = (numberOfPoints, aderdg.numberOfQuantities())
    QInterpolated = OptionalDimTensor(
        "QInterpolated",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        gShape,
        alignStride=True,
    )

    stressRotationMatrix = Tensor("stressRotationMatrix", (6, 6))
    initialStress = Tensor("initialStress", (6,))
    rotatedStress = Tensor("rotatedStress", (6,))
    rotationKernel = (
        rotatedStress["i"] <= stressRotationMatrix["ij"] * initialStress["j"]
    )
    generator.add("rotateStress", rotationKernel)

    reducedFaceAlignedMatrix = Tensor("reducedFaceAlignedMatrix", (6, 6))
    generator.add(
        "rotateInitStress",
        rotatedStress["k"]
        <= stressRotationMatrix["ki"]
        * reducedFaceAlignedMatrix["ij"]
        * initialStress["j"],
    )

    originalQ = Tensor("originalQ", (numberOfPoints,))
    resampledQ = Tensor("resampledQ", (numberOfPoints,))
    resampleKernel = resampledQ["i"] <= db.resample["ij"] * originalQ["j"]
    generator.add("resampleParameter", resampleKernel)

    generator.add("transposeTinv", TinvT["ij"] <= aderdg.Tinv["ji"])

    fluxScale = Scalar("fluxScaleDR")
    generator.add(
        "rotateFluxMatrix",
        fluxSolver["qp"] <= fluxScale * aderdg.starMatrix(0)["qk"] * aderdg.T["pk"],
    )

    numberOf3DBasisFunctions = aderdg.numberOf3DBasisFunctions()
    numberOfQuantities = aderdg.numberOfQuantities()
    basisFunctionsAtPoint = Tensor("basisFunctionsAtPoint", (numberOf3DBasisFunctions,))
    QAtPoint = OptionalDimTensor(
        "QAtPoint",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (numberOfQuantities,),
    )

    generator.add(
        "evaluateFaceAlignedDOFSAtPoint",
        QAtPoint["q"]
        <= aderdg.Tinv["qp"] * aderdg.Q["lp"] * basisFunctionsAtPoint["l"],
    )

    def interpolateQGenerator(i, h):
        return (
            QInterpolated["kp"]
            <= db.V3mTo2n[i, h][aderdg.t("kl")] * aderdg.Q["lq"] * TinvT["qp"]
        )

    interpolateQPrefetch = lambda i, h: QInterpolated
    for target in targets:
        name_prefix = generate_kernel_name_prefix(target)
        generator.addFamily(
            f"{name_prefix}evaluateAndRotateQAtInterpolationPoints",
            simpleParameterSpace(4, 4),
            interpolateQGenerator,
            interpolateQPrefetch if target == "cpu" else None,
            target=target,
        )

    nodalFluxGenerator = (
        lambda i, h: aderdg.extendedQTensor()["kp"]
        <= aderdg.extendedQTensor()["kp"]
        + db.V3mTo2nTWDivM[i, h][aderdg.t("kl")]
        * QInterpolated["lq"]
        * fluxSolver["qp"]
    )
    nodalFluxPrefetch = lambda i, h: aderdg.I

    for target in targets:
        name_prefix = generate_kernel_name_prefix(target)
        generator.addFamily(
            f"{name_prefix}nodalFlux",
            simpleParameterSpace(4, 4),
            nodalFluxGenerator,
            nodalFluxPrefetch if target == "cpu" else None,
            target=target,
        )

    # Energy output
    # Minus and plus refer to the original implementation of Christian Pelties,
    # where the normal points from the plus side to the minus side
    QInterpolatedPlus = OptionalDimTensor(
        "QInterpolatedPlus",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        gShape,
        alignStride=True,
    )
    QInterpolatedMinus = OptionalDimTensor(
        "QInterpolatedMinus",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        gShape,
        alignStride=True,
    )
    slipInterpolated = Tensor("slipInterpolated", (numberOfPoints, 3), alignStride=True)
    tractionInterpolated = Tensor(
        "tractionInterpolated", (numberOfPoints, 3), alignStride=True
    )
    staticFrictionalWork = Tensor("staticFrictionalWork", ())
    minusSurfaceArea = Scalar("minusSurfaceArea")
    spaceWeights = Tensor("spaceWeights", (numberOfPoints,), alignStride=True)

    computeTractionInterpolated = (
        tractionInterpolated["kp"]
        <= QInterpolatedMinus["kq"] * aderdg.tractionMinusMatrix["qp"]
        + QInterpolatedPlus["kq"] * aderdg.tractionPlusMatrix["qp"]
    )
    generator.add("computeTractionInterpolated", computeTractionInterpolated)

    accumulateStaticFrictionalWork = (
        staticFrictionalWork[""]
        <= staticFrictionalWork[""]
        + minusSurfaceArea
        * tractionInterpolated["kp"]
        * slipInterpolated["kp"]
        * spaceWeights["k"]
    )
    generator.add("accumulateStaticFrictionalWork", accumulateStaticFrictionalWork)

    # Dynamic Rupture Precompute
    qPlus = OptionalDimTensor(
        "Qplus",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        gShape,
        alignStride=True,
    )
    qMinus = OptionalDimTensor(
        "Qminus",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        gShape,
        alignStride=True,
    )

    extractVelocitiesSPP = aderdg.extractVelocities()
    extractVelocities = Tensor(
        "extractVelocities",
        extractVelocitiesSPP.shape,
        spp=extractVelocitiesSPP,
    )
    extractTractionsSPP = aderdg.extractTractions()
    extractTractions = Tensor(
        "extractTractions", extractTractionsSPP.shape, spp=extractTractionsSPP
    )

    N = extractTractionsSPP.shape[0]
    eta = Tensor("eta", (N, N))
    zPlus = Tensor("Zplus", (N, N))
    zMinus = Tensor("Zminus", (N, N))
    theta = OptionalDimTensor(
        "theta",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos,
        (numberOfPoints, N),
        alignStride=True,
    )

    velocityJump = (
        extractVelocities["lj"] * qMinus["ij"] - extractVelocities["lj"] * qPlus["ij"]
    )
    tractionsPlus = extractTractions["mn"] * qPlus["in"]
    tractionsMinus = extractTractions["mn"] * qMinus["in"]
    computeTheta = (
        theta["ik"]
        <= eta["kl"] * velocityJump
        + eta["kl"] * zPlus["lm"] * tractionsPlus
        + eta["kl"] * zMinus["lm"] * tractionsMinus
    )
    generator.add("computeTheta", computeTheta)

    mapToVelocitiesSPP = aderdg.mapToVelocities()
    mapToVelocities = Tensor(
        "mapToVelocities", mapToVelocitiesSPP.shape, spp=mapToVelocitiesSPP
    )
    mapToTractionsSPP = aderdg.mapToTractions()
    mapToTractions = Tensor(
        "mapToTractions", mapToTractionsSPP.shape, spp=mapToTractionsSPP
    )
    imposedState = OptionalDimTensor(
        "imposedState",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        gShape,
        alignStride=True,
    )
    weight = Scalar("weight")
    computeImposedStateM = (
        imposedState["ik"]
        <= imposedState["ik"]
        + weight
        * mapToVelocities["kl"]
        * (
            extractVelocities["lm"] * qMinus["im"]
            - zMinus["lm"] * theta["im"]
            + zMinus["lm"] * tractionsMinus
        )
        + weight * mapToTractions["kl"] * theta["il"]
    )
    computeImposedStateP = (
        imposedState["ik"]
        <= imposedState["ik"]
        + weight
        * mapToVelocities["kl"]
        * (
            extractVelocities["lm"] * qPlus["im"]
            - zPlus["lm"] * tractionsPlus
            + zPlus["lm"] * theta["im"]
        )
        + weight * mapToTractions["kl"] * theta["il"]
    )
    generator.add("computeImposedStateM", computeImposedStateM)
    generator.add("computeImposedStateP", computeImposedStateP)

    return {db.resample, db.quadpoints, db.quadweights}
