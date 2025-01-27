# SPDX-FileCopyrightText: 2017-2024 SeisSol Group
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
    numberOfNodes = aderdg.t(db.v.shape())[0]

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

    sShape_eta = (numberOf3DBasisFunctions,)
    QEtaModal = OptionalDimTensor(
        "QEtaModal",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        sShape_eta,
        alignStride=True,
    )

    initialLoading = Tensor("initialLoading", (6,))

    replicateIniLShape = (numberOfNodes,)
    replicateIniLSpp = np.ones(
        aderdg.Q.insertOptDim(replicateIniLShape, (aderdg.Q.optSize(),))
    )
    replicateInitialLoading = OptionalDimTensor(
        "replicateInitialLoading",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        replicateIniLShape,
        spp=replicateIniLSpp,
        alignStride=True,
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
    yieldFactor = Tensor("yieldFactor", (numberOfNodes,))

    generator.add(
        "plConvertToNodal",
        QStressNodal["kp"]
        <= db.v[aderdg.t("kl")] * QStress["lp"]
        + replicateInitialLoading["k"] * initialLoading["p"],
    )

    for target in targets:
        name_prefix = generate_kernel_name_prefix(target)
        generator.add(
            name=f"{name_prefix}plConvertToNodalNoLoading",
            ast=QStressNodal["kp"] <= db.v[aderdg.t("kl")] * QStress["lp"],
            target=target,
        )

        generator.add(
            name=f"{name_prefix}plConvertEtaModal2Nodal",
            ast=QEtaNodal["k"] <= db.v[aderdg.t("kl")] * QEtaModal["l"],
            target=target,
        )

        generator.add(
            name=f"{name_prefix}plConvertEtaNodal2Modal",
            ast=QEtaModal["k"] <= db.vInv[aderdg.t("kl")] * QEtaNodal["l"],
            target=target,
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
        "plAdjustStresses",
        QStress["kp"]
        <= QStress["kp"]
        + db.vInv[aderdg.t("kl")] * QStressNodal["lp"] * yieldFactor["l"],
    )

    gpu_target = "gpu"
    if gpu_target in targets:
        name_prefix = generate_kernel_name_prefix(gpu_target)

        # suffix `M` stands for `Matrix`
        replicateInitialLoadingM = Tensor(
            name="replicateInitialLoadingM",
            shape=(numberOfNodes, 1),
            spp=np.ones((numberOfNodes, 1)),
        )
        initialLoadingM = Tensor("initialLoadingM", (1, 6))

        # Note: the last term was change on purpose because
        # GemmForge doesn't currently support tensor product operation
        convert_to_nodal = (
            QStressNodal["kp"]
            <= db.v[aderdg.t("kl")] * QStress["lp"]
            + replicateInitialLoadingM["km"] * initialLoadingM["mp"]
        )

        generator.add(
            name=f"{name_prefix}plConvertToNodal",
            ast=convert_to_nodal,
            target=gpu_target,
        )

        generator.add(
            f"{name_prefix}plConvertToModal",
            QStress["kp"]
            <= QStress["kp"] + db.vInv[aderdg.t("kl")] * QStressNodal["lp"],
            target=gpu_target,
        )
