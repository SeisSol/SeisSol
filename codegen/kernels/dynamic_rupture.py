# SPDX-FileCopyrightText: 2016 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff

import numpy as np
from kernels.common import generate_kernel_name_prefix
from kernels.multsim import OptionalDimTensor
from yateto import Scalar, Tensor, simpleParameterSpace
from yateto.input import parseJSONMatrixFile
from yateto.memory import CSCMemoryLayout, PatternMemoryLayout
from yateto.type import AddressingMode

# The face relation index of a dynamic rupture face: 0 selects the plus side, 1 the minus side.
# The minus side carries the face orientation index of the shared face, which the canonical
# vertex numbering pins to zero.
NumFaceRelations = 2


def addKernels(
    generator,
    aderdg,
    matricesDir,
    drQuadRule,
    targets,
    isOldGpuInterface,
):

    clones = dict()

    # Load matrices
    db = parseJSONMatrixFile(
        f"{matricesDir}/dr_{drQuadRule}_matrices_{aderdg.order}.json",
        clones,
        alignStride=aderdg.alignStride,
        transpose=aderdg.transpose,
    )
    numPoints = aderdg.t(db.resample.shape())[0]

    # Determine matrices
    # Note: This does only work because the flux does not depend
    # on the mechanisms in the case of viscoelastic attenuation
    trans_inv_spp_T = aderdg.transformation_inv_spp().transpose()
    TinvT = Tensor("TinvT", trans_inv_spp_T.shape, spp=trans_inv_spp_T)
    # The face rotation is block diagonal -- one block per quantity group -- so
    # most of TinvT is structurally zero, and it is stored once per fault face.
    # Storing only the pattern shrinks that and lets the two projections below
    # skip the empty blocks. The old GPU interface (gemmforge/chainforge) reads
    # its operands as dense, so it keeps the dense layout.
    if not (isOldGpuInterface and "gpu" in targets):
        TinvT.setMemoryLayout(CSCMemoryLayout)

    # A curved fault face turns its frame along it, so it carries its rotation
    # at each of its points: into the coordinates of the fault, and back out of
    # them in the lift. Every point has the pattern of the face, packed by it.
    curvilinear = getattr(aderdg, "curvilinear", False)
    TinvTPoints = Tensor(
        "TinvTPoints",
        (numPoints,) + tuple(trans_inv_spp_T.shape),
        spp=np.broadcast_to(
            trans_inv_spp_T, (numPoints,) + tuple(trans_inv_spp_T.shape)
        ).copy(),
    )
    trans_spp = aderdg.transformation_spp()
    TPoints = Tensor(
        "TPoints",
        (numPoints,) + tuple(trans_spp.shape),
        spp=np.broadcast_to(trans_spp, (numPoints,) + tuple(trans_spp.shape)).copy(),
    )
    if not (isOldGpuInterface and "gpu" in targets):
        TinvTPoints.setMemoryLayout(PatternMemoryLayout, alignStride=False)
        TPoints.setMemoryLayout(PatternMemoryLayout, alignStride=False)
    # the rotation into the coordinates of the fault, contracted with the
    # values of a side at the points (index k) and its quantities (index q)
    rotationIntoFault = TinvTPoints["kqp"] if curvilinear else TinvT["qp"]

    flux_solver_spp = aderdg.flux_solver_spp()
    fluxSolver = Tensor("fluxSolver", flux_solver_spp.shape, spp=flux_solver_spp)

    gShape = (numPoints, aderdg.numQuantities())
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

    originalQ = OptionalDimTensor(
        "originalQ",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (numPoints,),
        alignStride=True,
    )
    resampledQ = OptionalDimTensor(
        "resampledQ",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (numPoints,),
        alignStride=True,
    )
    resampleKernel = resampledQ["i"] <= db.resample[aderdg.t("ij")] * originalQ["j"]
    generator.add("resampleParameter", resampleKernel)

    fluxScale = Scalar("fluxScaleDR")
    generator.add(
        "rotateFluxMatrix",
        fluxSolver["qp"]
        <= fluxScale * aderdg.starMatrixSetup(0)["qk"] * aderdg.T["pk"],
    )

    num3DBasisFunctions = aderdg.num3DBasisFunctions()
    numQuantities = aderdg.numQuantities()
    basisFunctionsAtPoint = Tensor("basisFunctionsAtPoint", (num3DBasisFunctions,))
    QAtPoint = OptionalDimTensor(
        "QAtPoint",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (numQuantities,),
    )

    generator.add(
        "evaluateFaceAlignedDOFSAtPoint",
        QAtPoint["q"]
        <= aderdg.Tinv["qp"] * aderdg.Q["lp"] * basisFunctionsAtPoint["l"],
    )

    # the same from the trace of a side on the face, at a point of the face
    faceBasisAtPoint = Tensor("faceBasisAtPoint", (aderdg.num2DBasisFunctions(),))
    generator.add(
        "evaluateFaceAlignedTraceAtPoint",
        QAtPoint["q"]
        <= aderdg.Tinv["qp"] * aderdg.faceTraceNeighbor["lp"] * faceBasisAtPoint["l"],
    )

    # Either side at the points of the fault, from the traces it writes on the
    # fault face, both in the parametrisation of the plus side (the minus side
    # writes its traces reflected, see _configureTraces): one matrix for both
    # sides and every face of a cell, the face basis at the points of the fault.
    traceToDR = faultTraceEvaluation(aderdg, db)
    generator.add(
        "evaluateAndRotateTraceAtInterpolationPoints",
        QInterpolated["kp"]
        <= traceToDR[aderdg.t("km")]
        * aderdg.faceTraceNeighbor["mq"]
        * rotationIntoFault,
        target="cpu",
    )

    steps = aderdg.order
    scalars = [
        [Scalar(f"coeffDR({i * aderdg.order + p})") for p in range(aderdg.order)]
        for i in range(steps)
    ]
    QDR = [
        OptionalDimTensor(
            f"QDR({i})",
            aderdg.Q.optName(),
            aderdg.Q.optSize(),
            aderdg.Q.optPos(),
            gShape,
            alignStride=True,
        )
        for i in range(steps)
    ]

    def projectTraces():
        # all time points of a side at once, from the traces of its time
        # derivatives (faceTraceAtTime)
        calc = []
        for c in range(steps):
            interm = aderdg.faceTraceAtTime(scalars[c])

            if isOldGpuInterface:
                # the "old" GPU implementation (gemmforge/chainforge) needs an explicit intermediate
                calc += [aderdg.faceTraceNeighbor["mq"] <= interm]
                interm = aderdg.faceTraceNeighbor["mq"]

            calc += [
                QDR[c]["kp"] <= traceToDR[aderdg.t("km")] * interm * rotationIntoFault
            ]
        return calc

    for target in targets:
        name_prefix = generate_kernel_name_prefix(target)
        generator.add(
            f"{name_prefix}projectTracesToDR",
            projectTraces(),
            target=target,
        )

    faultFlux = faultFluxTensors(aderdg, numPoints)
    if faultFlux is None:
        nodalFluxGenerator = (
            lambda i, h: aderdg.extendedQTensor()["kp"]
            <= aderdg.extendedQTensor()["kp"]
            + db.V3mTo2nTWDivM[i, h][aderdg.t("kl")]
            * QInterpolated["lq"]
            * fluxSolver["qp"]
        )
    else:
        # the rotation back out of the coordinates of the fault, contracted
        # with the lifted state at the points (index l)
        rotationOutOfFault = TPoints["lpq"] if curvilinear else aderdg.T["pq"]
        nodalFluxGenerator = lambda i, h: pointwiseLift(
            aderdg,
            faultFlux,
            QInterpolated,
            db.V3mTo2nTWDivM[i, h][aderdg.t("kl")],
            rotationOutOfFault,
        )
    for target in targets:
        name_prefix = generate_kernel_name_prefix(target)
        # a device kernel writes each temporary once, see singleDefinitions
        with aderdg.singleDefinitions(target == "gpu"):
            generator.addFamily(
                f"{name_prefix}nodalFlux",
                simpleParameterSpace(4, NumFaceRelations),
                nodalFluxGenerator,
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
    slipInterpolated = OptionalDimTensor(
        "slipInterpolated",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (numPoints, 3),
        alignStride=True,
    )

    tractionInterpolated = OptionalDimTensor(
        "tractionInterpolated",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (numPoints, 3),
        alignStride=True,
    )
    staticFrictionalWork = OptionalDimTensor(
        "staticFrictionalWork",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (1,),
        alignStride=True,
    )
    minusSurfaceArea = Scalar("minusSurfaceArea")

    computeTractionInterpolated = (
        tractionInterpolated["kp"]
        <= QInterpolatedMinus["kq"] * aderdg.tractionMinusMatrix["qp"]
        + QInterpolatedPlus["kq"] * aderdg.tractionPlusMatrix["qp"]
    )
    generator.add("computeTractionInterpolated", computeTractionInterpolated)

    accumulateStaticFrictionalWork = (
        staticFrictionalWork["l"]
        <= staticFrictionalWork["l"]
        + minusSurfaceArea
        * tractionInterpolated["kp"]
        * slipInterpolated["kp"]
        * db.quadweights["k"]
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
        aderdg.Q.optPos(),
        (numPoints, N),
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

    # the layout of a fault face names both rotations, whichever the kernels use
    return {
        db.resample,
        db.quadpoints,
        db.quadweights,
        TinvT,
        TinvTPoints,
        TPoints,
    }


def faultTraceEvaluation(aderdg, db):
    """The face basis at the points of the fault, in the parametrisation of the
    plus side: what a trace written in it (see _configureTraces) gives there.

    It is what reads the field of a side at those points with the trace taken
    out, V3mTo2n(i, 0) = W rT(i) for every face i of the plus side, and for the
    minus side, which writes its traces reflected, V3mTo2n(i, 1) = W R rT(i).
    Checked here for every face, against the matrices the field is read with.
    """

    def matrix(tensor):
        values = np.zeros(tensor.shape())
        for index, value in tensor.values().items():
            values[index] = float(value)
        return values.T if aderdg.transpose(tensor.name()) else values

    def plain(tensor):
        values = np.zeros(tensor.shape())
        for index, value in tensor.values().items():
            values[index] = float(value)
        return values

    reflection = plain(aderdg.db.M2inv) @ matrix(aderdg.db.fP[0])
    evaluation = matrix(db.V3mTo2n[0, 0]) @ np.linalg.pinv(matrix(aderdg.db.rT[0]))
    for face in range(4):
        trace = matrix(aderdg.db.rT[face])
        for relation, sideTrace in ((0, trace), (1, reflection @ trace)):
            expected = matrix(db.V3mTo2n[face, relation])
            deviation = np.abs(evaluation @ sideTrace - expected).max()
            if deviation > 1e-10 * max(1.0, np.abs(expected).max()):
                raise RuntimeError(
                    "the points of the fault do not read a side through its face "
                    f"trace (face {face}, relation {relation}: {deviation})"
                )
    evaluation[np.abs(evaluation) < 1e-14] = 0.0
    stored = evaluation.T if aderdg.transpose("traceToDR") else evaluation
    return Tensor(
        "traceToDR",
        stored.shape,
        spp={
            index: repr(float(value))
            for index, value in np.ndenumerate(stored)
            if value != 0.0
        },
        alignStride=aderdg.alignStride("traceToDR"),
    )


def faultFluxTensors(aderdg, numPoints):
    """The operands of the lift where a fault face carries it per point.

    The lift turns the imposed state of a side, given in the coordinates of the
    face at its quadrature points, into what that side's cell receives: the
    coefficient matrix of the fault normal applied to it, rotated back, and
    projected into the cell. In face coordinates the matrix is the star of the
    first direction, which is a handful of scalars times fixed entries -- the
    decomposition the volume operator already reads. So a face carries those
    scalars at each point, together with its rotation, instead of one matrix
    per side that the rotation is folded into. Stored as scalars, a material
    that varies along the face costs a few numbers per point, where a matrix
    per point would cost the whole operator.

    None where the face keeps the one matrix per side, which is wherever the
    material does not vary inside a cell.
    """
    indices = aderdg.faultFluxCoefficients()
    if not indices:
        return None
    starShape = aderdg.starMatrixSetup(0).shape()
    structures = []
    for position, coefficient in enumerate(indices):
        values = {}
        for entry in aderdg.solverCoefficientEntries():
            if entry.dim == 0 and entry.coefficient == coefficient:
                index = (entry.row, entry.column)
                values[index] = values.get(index, 0.0) + entry.factor
        structures.append(
            Tensor(
                f"faultFluxStructure({position})",
                starShape,
                spp={index: repr(float(value)) for index, value in values.items()},
                addressing=AddressingMode.IMMEDIATE,
            )
        )
    coefficients = [
        Tensor(f"faultFluxCoefficients({position})", (numPoints,))
        for position in range(len(indices))
    ]
    product = aderdg.nodalTemporary("faultFluxProduct", (numPoints, starShape[1]))
    return {
        "structures": structures,
        "coefficients": coefficients,
        "product": product,
    }


def pointwiseLift(aderdg, faultFlux, imposedState, lift, rotation):
    """The statements of the lift where a face carries its scalars per point:
    the operator at each point, in face coordinates, then the rotation back
    -- one for the face, or one per point where the face may be curved -- and
    the projection into the cell."""
    product = aderdg.definedOnce(faultFlux["product"])
    statements = []
    for position, coefficient in enumerate(faultFlux["coefficients"]):
        term = (
            coefficient["l"]
            * imposedState["lq"]
            * faultFlux["structures"][position]["qk"]
        )
        statements.append(
            product["lk"] <= (term if position == 0 else product["lk"] + term)
        )
    target = aderdg.extendedQTensor()
    statements.append(target["kp"] <= target["kp"] + lift * product["lq"] * rotation)
    return statements


def addKernelsGeneral(generator):
    stressRotationMatrix = Tensor("stressRotationMatrix", (6, 6))
    initialStress = Tensor("initialStress", (6,))
    rotatedStress = Tensor("rotatedStress", (6,))
    rotationKernel = (
        rotatedStress["i"] <= stressRotationMatrix["ij"] * initialStress["j"]
    )
    generator.add("rotateStress", rotationKernel)
