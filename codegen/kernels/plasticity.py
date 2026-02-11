# SPDX-FileCopyrightText: 2017 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff

from kernels.common import generate_kernel_name_prefix
from kernels.multsim import OptionalDimTensor
from yateto import Tensor
from yateto.input import parseJSONMatrixFile


def addKernels(generator, aderdg, matricesDir, PlasticityMethod, targets):
    # Load matrices
    db = parseJSONMatrixFile(
        f"{matricesDir}/plasticity-{PlasticityMethod}-matrices-{aderdg.order}.json",
        clones=dict(),
        alignStride=aderdg.alignStride,
    )
    numberOfNodes = db.v.shape()[0]

    numberOf3DBasisFunctions = aderdg.numberOf3DBasisFunctions()
    sShape = (numberOf3DBasisFunctions, 6)
    QStress = OptionalDimTensor(
        "QStress",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        sShape,
        alignStride=True,
    )

    iShape = (numberOfNodes, 6)

    initialLoading = OptionalDimTensor(
        "initialLoading",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        iShape,
    )

    QStressNodal = OptionalDimTensor(
        "QStressNodal",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        iShape,
        alignStride=True,
    )
    meanStress = OptionalDimTensor(
        "meanStress",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (numberOfNodes,),
        alignStride=True,
    )
    secondInvariant = OptionalDimTensor(
        "secondInvariant",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (numberOfNodes,),
        alignStride=True,
    )

    QEtaNodal = OptionalDimTensor(
        "QEtaNodal",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (numberOfNodes,),
        alignStride=True,
    )
    QEtaNodalProject = OptionalDimTensor(
        "QEtaNodalProject",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (aderdg.numberOf3DQuadraturePoints(),),
        alignStride=True,
    )

    selectBulkAverage = Tensor(
        "selectBulkAverage", (6,), spp={(i,): str(1.0 / 3.0) for i in range(3)}
    )
    selectBulkNegative = Tensor(
        "selectBulkNegative", (6,), spp={(i,): "-1.0" for i in range(3)}
    )
    weightSecondInvariant = Tensor(
        "weightSecondInvariant",
        (6,),
        spp={(i,): str(1.0 / 2.0) if i < 3 else "1.0" for i in range(6)},
    )

    generator.add(
        "plComputeMean",
        meanStress["k"] <= QStressNodal["kq"] * selectBulkAverage["q"],
    )
    generator.add(
        "plSubtractMean",
        QStressNodal["kp"]
        <= QStressNodal["kp"] + meanStress["k"] * selectBulkNegative["p"],
    )
    generator.add(
        "plComputeSecondInvariant",
        secondInvariant["k"]
        <= QStressNodal["kq"] * QStressNodal["kq"] * weightSecondInvariant["q"],
    )

    # for the "old" output
    nodalVar = OptionalDimTensor(
        "nodalVar",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (numberOfNodes,),
        alignStride=True,
    )
    modalVar = OptionalDimTensor(
        "modalVar",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (numberOf3DBasisFunctions,),
        alignStride=True,
    )
    generator.add(
        "plOutput",
        modalVar["k"] <= db.vInv["kl"] * nodalVar["l"],
    )
    # end for the "old" output

    if QEtaNodal.shape() == QEtaNodalProject.shape():
        generator.add("plProject", QEtaNodalProject["l"] <= QEtaNodal["l"])
    else:
        generator.add(
            "plProject",
            QEtaNodalProject["p"]
            <= aderdg.db.evalAtQP[aderdg.t("pk")] * db.vInv["kl"] * QEtaNodal["l"],
        )

    for target in targets:
        name_prefix = generate_kernel_name_prefix(target)

        generator.add(
            f"{name_prefix}plConvertToNodal",
            QStressNodal["kp"] <= db.v["kl"] * QStress["lp"] + initialLoading["kp"],
            target=target,
        )

        generator.add(
            f"{name_prefix}plConvertToModal",
            QStress["kp"] <= QStress["kp"] + db.vInv["kl"] * QStressNodal["lp"],
            target=target,
        )


def includeTensors(matricesDir, aderdg, PlasticityMethod, includeTensors):
    db = parseJSONMatrixFile(
        f"{matricesDir}/plasticity-{PlasticityMethod}-matrices-{aderdg.order}.json",
        clones=dict(),
        alignStride=aderdg.alignStride,
    )
    includeTensors.add(db.vNodes)

    numberOfNodes = db.v.shape()[0]

    yieldFactor = OptionalDimTensor(
        "yieldFactor",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (numberOfNodes,),
        alignStride=True,
    )
    includeTensors.add(yieldFactor)
