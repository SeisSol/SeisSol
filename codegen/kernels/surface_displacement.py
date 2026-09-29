# SPDX-FileCopyrightText: 2017 SeisSol Group
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


def addKernels(generator, aderdg, include_tensors, targets):
    maxDepth = 3

    num3DBasisFunctions = aderdg.num3DBasisFunctions()
    num2DBasisFunctions = aderdg.num2DBasisFunctions()

    faceDisplacement = OptionalDimTensor(
        "faceDisplacement",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (num2DBasisFunctions, 3),
        alignStride=True,
    )
    averageNormalDisplacement = OptionalDimTensor(
        "averageNormalDisplacement",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (num2DBasisFunctions,),
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
            (4**depth, num3DBasisFunctions),
            alignStride=True,
        )
        for depth in range(maxDepth + 1)
    ]
    subTriangleProjectionFromFace = [
        Tensor(
            "subTriangleProjectionFromFace({})".format(depth),
            (4**depth, num2DBasisFunctions),
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
        (num2DBasisFunctions, 3),
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

    numQuadratureNodes = (aderdg.order + 1) ** 2
    rotatedFaceDisplacementAtQuadratureNodes = OptionalDimTensor(
        "rotatedFaceDisplacementAtQuadratureNodes",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (numQuadratureNodes, 3),
        alignStride=True,
    )
    generator.add(
        "rotateFaceDisplacementsAndEvaluateAtQuadratureNodes",
        rotatedFaceDisplacementAtQuadratureNodes["in"]
        <= aderdg.V2nTo2JacobiQuad["ij"]
        * rotatedFaceDisplacement["jp"]
        * displacementRotationMatrix["np"],
    )

    # Explicit temporary: without it the two rotated-displacement factors carry
    # different contraction indices and yateto cannot see that they are the same
    # subexpression, so it recomputes the rotation for both.
    faceDisplacementModal = OptionalDimTensor(
        "faceDisplacementModal",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (num2DBasisFunctions, 3),
        alignStride=True,
        temporary=True,
    )
    faceDisplacementSquared = OptionalDimTensor(
        "faceDisplacementSquared",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (3,),
    )
    generator.add(
        "faceDisplacementSquaredCompute",
        [
            faceDisplacementModal["mn"]
            <= aderdg.db.MV2nTo2m["mI"]
            * rotatedFaceDisplacement["Ip"]
            * displacementRotationMatrix["np"],
            faceDisplacementSquared["n"]
            <= aderdg.db.M2["ij"]
            * faceDisplacementModal["in"]
            * faceDisplacementModal["jn"],
        ],
    )

    if "gpu" in targets:
        name_prefix = generate_kernel_name_prefix(target="gpu")

        integratedVelocities = OptionalDimTensor(
            "integratedVelocities",
            aderdg.I.optName(),
            aderdg.I.optSize(),
            aderdg.I.optPos(),
            (num3DBasisFunctions, 3),
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
        (aderdg.num2DBasisFunctions(), 1),
        alignStride=True,
        temporary=True,
    )
    averageNormalDisplacement = OptionalDimTensor(
        "Iint",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (num2DBasisFunctions, 1),
        alignStride=True,
    )
    faceDisplacementTmp = OptionalDimTensor(
        "IAcc",
        aderdg.INodal.optName(),
        aderdg.INodal.optSize(),
        aderdg.INodal.optPos(),
        (aderdg.num2DBasisFunctions(), 3),
        alignStride=True,
        temporary=True,
    )
    # Number of displacement components that the Taylor series contributes to:
    # the normal one always, the two tangential ones only if the material
    # carries shear.
    numDisplacementComponents = 3 if aderdg.velocityOffset() > 1 else 1

    # Modal accumulators. The Taylor coefficients of eta are assembled in modal
    # space and evaluated at the face nodes once, at the end of the kernel.
    MPrev = OptionalDimTensor(
        "MPrev",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (num3DBasisFunctions, 1),
        alignStride=True,
        temporary=True,
    )
    MDisp = OptionalDimTensor(
        "MDisp",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (num3DBasisFunctions, numDisplacementComponents),
        alignStride=True,
        temporary=True,
    )
    MInt = OptionalDimTensor(
        "MInt",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (num3DBasisFunctions, 1),
        alignStride=True,
        temporary=True,
    )

    coeffs = [Scalar(f"coeff({i})") for i in range(aderdg.order + 1)]
    powers = [Scalar(f"fsgpower({i})") for i in range(aderdg.order + 1)]
    invImp = Tensor("invImp", ())
    rhoG = Tensor("rhoG", ())

    vidx = aderdg.velocityOffset()

    # The series reads the normal stress row and the velocity rows of Tinv,
    # and the displacement block of T. The face rotation may be stored by its
    # pattern (CSC), which yateto can neither slice at a row offset nor check
    # against a read that leaves out its first row. So each kernel gathers the
    # first row and the velocity rows of both matrices into dense temporaries,
    # through constant selection matrices, and slices those instead; the first
    # row of T is gathered only to keep that check satisfied. The temporaries
    # keep the pattern of the rows they hold, so that the products reading
    # them skip the zeros the rotation has there.
    def gatherRows(rotation, name):
        values = np.zeros((4, rotation.shape()[0]))
        values[0, 0] = 1.0
        for j in range(3):
            values[1 + j, vidx + j] = 1.0
        select = Tensor(f"fsgSelect{name}", values.shape, values)
        pattern = (values != 0).astype(int) @ rotation.spp().as_ndarray().astype(int)
        rows = Tensor(f"fsg{name}", pattern.shape, spp=pattern != 0, temporary=True)
        return rows, lambda: rows["qm"] <= select["qp"] * rotation["pm"]

    TinvRows, gatherTinvRows = gatherRows(aderdg.Tinv, "TinvRows")
    TRows, gatherTRows = gatherRows(aderdg.T, "TRows")

    for target in targets:
        name_prefix = generate_kernel_name_prefix(target)

        # This kernel does two things:
        # 1: Compute eta (for all three dimensions) at the end of the timestep
        # 2: Compute the integral of eta in normal direction over the timestep
        # We do this by building up the Taylor series of eta.
        # Eta is defined by the ODE eta_t = u^R - 1/Z * (rho g eta - p^R)
        # Compute coefficients by differentiating ODE recursively, e.g.:
        # eta_tt = u^R_t - 1/Z * (rho eta_t g - p^R_t)
        # and substituting the previous coefficient eta_t
        # This implementation sums up the Taylor series directly without storing
        # all coefficients.

        def kernelPerFace(f):
            # The recursion is affine with node-independent coefficients, so it
            # splits into a part driven by eta at the beginning of the timestep
            # and a part driven by the interior derivatives. The former only
            # picks up powers of gamma = -1/Z rho g and stays at the face nodes;
            # the latter is assembled in modal space and evaluated at the nodes
            # once, after the series is complete.
            kernel = [
                gatherTinvRows(),
                gatherTRows(),
                faceDisplacementTmp["mp"]
                <= faceDisplacement["mn"]
                * TinvRows["pn"].subslice("p", 1, 4).subslice("n", vidx, vidx + 3),
                Iprev["mp"] <= faceDisplacementTmp["mp"].subslice("p", 0, 1),
                averageNormalDisplacement["mp"]
                <= powers[0] * faceDisplacementTmp["mp"].subslice("p", 0, 1),
            ]

            for i in range(1, aderdg.order + 1):
                velocitiesU = aderdg.dQs[i - 1]["lm"] * TinvRows["pm"].subslice(
                    "p", 1, 2
                )
                pressure = aderdg.dQs[i - 1]["lm"] * TinvRows["pm"].subslice("p", 0, 1)

                if i == 1:
                    kernel += [MPrev["lp"] <= velocitiesU - invImp[""] * pressure]
                else:
                    kernel += [
                        MPrev["lp"]
                        <= velocitiesU
                        - invImp[""] * (rhoG[""] * MPrev["lp"] + pressure)
                    ]

                kernel += [
                    Iprev["nq"] <= -invImp[""] * (rhoG[""] * Iprev["nq"]),
                ]

                if i == 1:
                    kernel += [
                        MDisp["lp"].subslice("p", 0, 1) <= coeffs[i] * MPrev["lp"],
                        MInt["lp"] <= powers[i] * MPrev["lp"],
                    ]
                else:
                    kernel += [
                        MDisp["lp"].subslice("p", 0, 1)
                        <= MDisp["lp"].subslice("p", 0, 1) + coeffs[i] * MPrev["lp"],
                        MInt["lp"] <= MInt["lp"] + powers[i] * MPrev["lp"],
                    ]

                kernel += [
                    faceDisplacementTmp["nq"].subslice("q", 0, 1)
                    <= faceDisplacementTmp["nq"].subslice("q", 0, 1)
                    + coeffs[i] * Iprev["nq"],
                    averageNormalDisplacement["nq"]
                    <= averageNormalDisplacement["nq"] + powers[i] * Iprev["nq"],
                ]

                if aderdg.velocityOffset() > 1:
                    velocitiesVW = aderdg.dQs[i - 1]["lm"] * TinvRows["pm"].subslice(
                        "p", 2, 4
                    )
                    if i == 1:
                        kernel += [
                            MDisp["lp"].subslice("p", 1, 3) <= coeffs[i] * velocitiesVW
                        ]
                    else:
                        kernel += [
                            MDisp["lp"].subslice("p", 1, 3)
                            <= MDisp["lp"].subslice("p", 1, 3)
                            + coeffs[i] * velocitiesVW
                        ]

            displacementTarget = faceDisplacementTmp["nq"]
            if numDisplacementComponents < 3:
                displacementTarget = displacementTarget.subslice(
                    "q", 0, numDisplacementComponents
                )

            kernel += [
                displacementTarget
                <= displacementTarget
                + aderdg.db.V3mTo2nFace[f][aderdg.t("nl")] * MDisp["lq"],
                averageNormalDisplacement["nq"]
                <= averageNormalDisplacement["nq"]
                + aderdg.db.V3mTo2nFace[f][aderdg.t("nl")] * MInt["lq"],
            ]

            kernel += [
                faceDisplacement["mp"]
                <= faceDisplacementTmp["mn"]
                * TRows["pn"].subslice("p", 1, 4).subslice("n", vidx, vidx + 3),
            ]

            return kernel

        generator.addFamily(
            f"{name_prefix}fsgKernel",
            simpleParameterSpace(4),
            kernelPerFace,
            target=target,
        )
