# SPDX-FileCopyrightText: 2019 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff
# SPDX-FileContributor: Sebastian Wolf

import kernels.equations.acoustic as acoustic
import numpy as np
from kernels.multsim import OptionalDimTensor
from yateto import Scalar, Tensor


def addKernels(generator, aderdg):
    numberOf3DBasisFunctions = aderdg.numberOf3DBasisFunctions()
    numberOfQuantities = aderdg.numberOfQuantities()
    order = aderdg.order
    # Point sources
    mStiffnessTensor = Tensor("stiffnessTensor", (3, 3, 3, 3))
    mNormal = Tensor("mNormal", (3,))
    mArea = Scalar("mArea")
    basisFunctionsAtPoint = Tensor("basisFunctionsAtPoint", (numberOf3DBasisFunctions,))
    basisFunctionDerivativesAtPoint = Tensor(
        "basisFunctionDerivativesAtPoint", (numberOf3DBasisFunctions, 3)
    )
    timeBasisFunctionsAtPoint = Tensor("timeBasisFunctionsAtPoint", (order,))
    mInvJInvPhisAtSources = Tensor(
        "mInvJInvPhisAtSources", (numberOf3DBasisFunctions,), alignStride=True
    )
    JInv = Scalar("JInv")

    generator.add(
        "computeMInvJInvPhisAtSources",
        mInvJInvPhisAtSources["k"]
        <= JInv * aderdg.db.M3inv["kl"] * basisFunctionsAtPoint["l"],
    )

    # extract the moment tensors entries in SeisSol ordering
    # i.e.: (xx, yy, zz, xy, yz, xz)
    if not isinstance(aderdg, acoustic.AcousticADERDG):
        assert numberOfQuantities >= 6
        momentToNRF_spp = np.zeros((numberOfQuantities, 3, 3))
        momentToNRF_spp[0, 0, 0] = 1
        momentToNRF_spp[1, 1, 1] = 1
        momentToNRF_spp[2, 2, 2] = 1
        momentToNRF_spp[3, 0, 1] = 1
        momentToNRF_spp[4, 1, 2] = 1
        momentToNRF_spp[5, 0, 2] = 1
    else:
        momentToNRF_spp = np.zeros((numberOfQuantities, 3, 3))
        momentToNRF_spp[0, 0, 0] = 1
    momentToNRF = Tensor("momentToNRF", (numberOfQuantities, 3, 3), spp=momentToNRF_spp)

    rotateNRF = Tensor("rotateNRF", (3, 3))
    momentNRFKernel = (
        momentToNRF["tpq"]
        * mArea
        * mStiffnessTensor["pqIj"]
        * mNormal["j"]
        * rotateNRF["Ii"]
    )

    tensorNRF = Tensor("tensorNRF", (numberOfQuantities, 3))

    generator.add("transformNRF", tensorNRF["ti"] <= momentNRFKernel)

    update = Tensor("update", (numberOfQuantities,))

    if aderdg.Q.hasOptDim():
        generator.add(
            "addPointSource",
            aderdg.Q["kt"]
            <= aderdg.Q["kt"]
            + mInvJInvPhisAtSources["k"] * update["t"] * aderdg.oneSimToMultSim["s"],
        )
    else:
        generator.add(
            "addPointSource",
            aderdg.Q["kt"] <= aderdg.Q["kt"] + mInvJInvPhisAtSources["k"] * update["t"],
        )

    # Receiver output
    QAtPoint = OptionalDimTensor(
        "QAtPoint",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (numberOfQuantities,),
    )
    evaluateDOFSAtPoint = QAtPoint["p"] <= aderdg.Q["kp"] * basisFunctionsAtPoint["k"]
    generator.add("evaluateDOFSAtPoint", evaluateDOFSAtPoint)
    QDerivativeAtPoint = OptionalDimTensor(
        "QDerivativeAtPoint",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (numberOfQuantities, 3),
    )
    evaluateDerivativeDOFSAtPoint = (
        QDerivativeAtPoint["pd"]
        <= aderdg.Q["kp"] * basisFunctionDerivativesAtPoint["kd"]
    )
    generator.add("evaluateDerivativeDOFSAtPoint", evaluateDerivativeDOFSAtPoint)

    stpShape = (numberOf3DBasisFunctions, numberOfQuantities, order)
    spaceTimePredictor = OptionalDimTensor(
        "spaceTimePredictor",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        stpShape,
        alignStride=True,
    )
    evaluateDOFSAtPointSTP = (
        QAtPoint["p"]
        <= spaceTimePredictor["kpt"]
        * basisFunctionsAtPoint["k"]
        * timeBasisFunctionsAtPoint["t"]
    )
    generator.add("evaluateDOFSAtPointSTP", evaluateDOFSAtPointSTP)
    spaceTimePredictor = OptionalDimTensor(
        "spaceTimePredictor",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        stpShape,
        alignStride=True,
    )
    evaluateDerivativeDOFSAtPointSTP = (
        QDerivativeAtPoint["pd"]
        <= spaceTimePredictor["kpt"]
        * basisFunctionDerivativesAtPoint["kd"]
        * timeBasisFunctionsAtPoint["t"]
    )
    generator.add("evaluateDerivativeDOFSAtPointSTP", evaluateDerivativeDOFSAtPointSTP)
