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
from yateto import Tensor
from yateto.input import parseXMLMatrixFile


def addKernels(generator, aderdg, matricesDir, PlasticityMethod, targets):
    # Load matrices
    db = parseXMLMatrixFile(
        "{}/plasticity_{}_matrices_{}.xml".format(
            matricesDir, PlasticityMethod, aderdg.order
        ),
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

    initialLoading = OptionalDimTensor(
        "initialLoading",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (6,),
    )

    replicateIniLShape = (numberOfNodes,)
    replicateIniLSpp = np.ones(replicateIniLShape)

    replicateInitialLoading = Tensor(
        "replicateInitialLoading",
        replicateIniLShape,
        spp=replicateIniLSpp,
        alignStride=aderdg.multipleSimulations == 1,
    )

    iShape = (numberOfNodes, 6)
    QStressNodal = OptionalDimTensor(
        "QStressNodal",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        iShape,
        alignStride=True,
    )
    QStressNodalPrev = OptionalDimTensor(
        "QStressNodalPrev",
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

    yieldFactor = OptionalDimTensor(
        "yieldFactor",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (numberOfNodes,),
    )

    generator.add(
        "plConvertToNodal",
        [
            QStressNodalPrev["kp"] <= db.v["kl"] * QStress["lp"],
            QStressNodal["kp"]
            <= QStressNodalPrev["kp"]
            + replicateInitialLoading["k"] * initialLoading["p"],
        ],
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

    generator.add(
        "plConvertToModal",
        QStress["kp"] <= QStress["kp"] + db.vInv["kl"] * QStressNodal["lp"],
    )

    generator.add(
        "plKeep",
        [QEtaNodal["k"] <= QEtaNodal["k"], yieldFactor["k"] <= yieldFactor["k"]],
    )

    gpu_target = "gpu"
    if gpu_target in targets:
        name_prefix = generate_kernel_name_prefix(gpu_target)

        if aderdg.multipleSimulations > 1:
            # for now, copy the tensors into here and rename them; until gemmforge/chainforge is deprecated
            initialLoadingM = OptionalDimTensor(
                "initialLoadingM",
                aderdg.Q.optName(),
                aderdg.Q.optSize(),
                aderdg.Q.optPos(),
                (6,),
            )
            replicateInitialLoadingM = Tensor(
                "replicateInitialLoadingM",
                replicateIniLShape,
                spp=replicateIniLSpp,
                alignStride=False,
            )

            matreplace = replicateInitialLoadingM["k"] * initialLoadingM["p"]
        else:
            # suffix `M` stands for `Matrix`
            replicateInitialLoadingM = Tensor(
                name="replicateInitialLoadingM",
                shape=(numberOfNodes, 1),
                spp=np.ones((numberOfNodes, 1)),
            )
            initialLoadingM = Tensor("initialLoadingM", (1, 6))

            matreplace = replicateInitialLoadingM["km"] * initialLoadingM["mp"]

        # Note: the last term was changed on purpose because
        # GemmForge doesn't currently support tensor product operation
        convert_to_nodal = [
            QStressNodalPrev["kp"] <= db.v["kl"] * QStress["lp"],
            QStressNodal["kp"] <= QStressNodalPrev["kp"] + matreplace,
        ]

        generator.add(
            name=f"{name_prefix}plConvertToNodal",
            ast=convert_to_nodal,
            target=gpu_target,
        )

        generator.add(
            f"{name_prefix}plConvertToModal",
            QStress["kp"] <= QStress["kp"] + db.vInv["kl"] * QStressNodal["lp"],
            target=gpu_target,
        )
