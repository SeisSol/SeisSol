# SPDX-FileCopyrightText: 2019 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff
# SPDX-FileContributor: Sebastian Wolf

from kernels.multsim import OptionalDimTensor
from kernels.quantities import moment_tensor_selector
from yateto import Scalar, Tensor


def addKernels(generator, aderdg):
    num3DBasisFunctions = aderdg.num3DBasisFunctions()
    numQuantities = aderdg.numQuantities()
    # Point sources
    mStiffnessTensor = Tensor("stiffnessTensor", (3, 3, 3, 3))
    mNormal = Tensor("mNormal", (3,))
    mArea = Scalar("mArea")
    basisFunctionsAtPoint = Tensor("basisFunctionsAtPoint", (num3DBasisFunctions,))
    basisFunctionDerivativesAtPoint = Tensor(
        "basisFunctionDerivativesAtPoint", (num3DBasisFunctions, 3)
    )
    mInvJInvPhisAtSources = Tensor(
        "mInvJInvPhisAtSources", (num3DBasisFunctions,), alignStride=True
    )
    JInv = Scalar("JInv")

    generator.add(
        "computeMInvJInvPhisAtSources",
        mInvJInvPhisAtSources["k"]
        <= JInv * aderdg.db.M3inv["kl"] * basisFunctionsAtPoint["l"],
    )

    # extract the moment tensor entries in SeisSol ordering
    # i.e.: (xx, yy, zz, xy, yz, xz); a scalar stress (acoustic) takes xx only
    momentToNRF_spp = moment_tensor_selector(aderdg.quantityBlocks())
    momentToNRF = Tensor("momentToNRF", (numQuantities, 3, 3), spp=momentToNRF_spp)

    rotateNRF = Tensor("rotateNRF", (3, 3))
    momentNRFKernel = (
        momentToNRF["tpq"]
        * mArea
        * mStiffnessTensor["pqIj"]
        * mNormal["j"]
        * rotateNRF["Ii"]
    )

    tensorNRF = Tensor("tensorNRF", (numQuantities, 3))

    generator.add("transformNRF", tensorNRF["ti"] <= momentNRFKernel)

    update = Tensor("update", (numQuantities,))

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
        (numQuantities,),
    )
    evaluateDOFSAtPoint = QAtPoint["p"] <= aderdg.Q["kp"] * basisFunctionsAtPoint["k"]
    generator.add("evaluateDOFSAtPoint", evaluateDOFSAtPoint)
    QDerivativeAtPoint = OptionalDimTensor(
        "QDerivativeAtPoint",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (numQuantities, 3),
    )
    evaluateDerivativeDOFSAtPoint = (
        QDerivativeAtPoint["pd"]
        <= aderdg.Q["kp"] * basisFunctionDerivativesAtPoint["kd"]
    )
    generator.add("evaluateDerivativeDOFSAtPoint", evaluateDerivativeDOFSAtPoint)
