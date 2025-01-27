# SPDX-FileCopyrightText: 2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

import numpy as np
from kernels.aderdg import LinearADERDG
from kernels.multsim import OptionalDimTensor
from yateto import Scalar, Tensor
from yateto.input import (memoryLayoutFromFile, parseJSONMatrixFile,
                          parseXMLMatrixFile)


def choose(n, k):
    num = np.prod(np.arange(n, n - k, -1))
    denom = np.prod(np.arange(1, k + 1))
    return num // denom


class PoroelasticADERDG(LinearADERDG):
    def __init__(
        self,
        order,
        multipleSimulations,
        matricesDir,
        memLayout,
        numberOfMechanisms,
        **kwargs
    ):

        super().__init__(order, multipleSimulations, matricesDir)
        clones = {
            "star": ["star(0)", "star(1)", "star(2)"],
        }
        self.db.update(
            parseXMLMatrixFile(
                "{}/matrices_poroelastic.xml".format(matricesDir), clones
            )
        )
        self.db.update(
            parseJSONMatrixFile("{}/stp_{}.json".format(matricesDir, order), clones)
        )

        memoryLayoutFromFile(memLayout, self.db, clones)

        self.kwargs = kwargs

    def numberOfQuantities(self):
        return 13

    def numberOfExtendedQuantities(self):
        return self.numberOfQuantities()

    def starMatrix(self, dim):
        return self.db.star[dim]

    def sourceMatrix(self):
        return self.db.ET

    def extractVelocities(self):
        extractVelocitiesSPP = np.zeros((4, self.numberOfQuantities()))
        extractVelocitiesSPP[0, 6] = 1
        extractVelocitiesSPP[1, 7] = 1
        extractVelocitiesSPP[2, 8] = 1
        extractVelocitiesSPP[3, 10] = 1
        return extractVelocitiesSPP

    def extractTractions(self):
        extractTractionsSPP = np.zeros((4, self.numberOfQuantities()))
        extractTractionsSPP[0, 0] = 1
        extractTractionsSPP[1, 3] = 1
        extractTractionsSPP[2, 5] = 1
        extractTractionsSPP[3, 9] = 1
        return extractTractionsSPP

    def name(self):
        return "poroelastic"

    def transformationSpp(self):
        spp = np.zeros(
            (self.numberOfQuantities(), self.numberOfQuantities()), dtype=bool
        )
        spp[0:6, 0:6] = 1
        spp[6:9, 6:9] = 1
        spp[9, 9] = 1
        spp[10:13, 10:13] = 1
        return spp

    def transformationInvSpp(self):
        return self.transformationSpp()

    def addTime(self, generator, targets):
        super().addTime(generator, targets)

        stiffnessValues = [self.db.kDivMT[d].values_as_ndarray() for d in range(3)]
        fullShape = (
            self.numberOf3DBasisFunctions(),
            self.numberOf3DBasisFunctions(),
        )

        stpShape = (
            self.numberOf3DBasisFunctions(),
            self.numberOfQuantities(),
            self.order,
        )
        spaceTimePredictorRhs = OptionalDimTensor(
            "spaceTimePredictorRhs",
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            stpShape,
            alignStride=True,
        )
        spaceTimePredictor = OptionalDimTensor(
            "spaceTimePredictor",
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            stpShape,
            alignStride=True,
        )
        testRhs = OptionalDimTensor(
            "testRhs",
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            stpShape,
            alignStride=True,
        )
        testLhs = OptionalDimTensor(
            "testLhs",
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            stpShape,
            alignStride=True,
        )
        timestep = Scalar("timestep")
        G = {10: Scalar("Gk"), 11: Scalar("Gl"), 12: Scalar("Gm")}

        # Compute the index range for basis functions of a certain degree
        #
        # The basis functions are ordered with increasing degree, i.e. the first
        # basis function has degree 0, the next three basis functions have degree
        # 1, the next six basis functions have degree 2 and so forth.
        # This method computes the indices Bn_lower, Bn_upper, such that
        # forall Bn_lower =< i < Bn_upper: degree(phi_i) == n
        #
        # @param n The desired polynomial degree
        def modeRange(n):
            Bn_lower = choose(n - 1 + 3, 3)
            Bn_upper = choose(n + 3, 3)
            return (Bn_lower, Bn_upper)

        # Compute a matrix, which filters out all basis function of degree n
        #
        # Compute a matrix M, such that for any DOF vector Q, M*Q only contains
        # the parts of Q, which correspond to basis functions of degree n
        #
        # @param n The desired polynomial degree
        def selectModes(n):
            Bn_1, Bn = modeRange(n)
            selectModesSpp = np.zeros(fullShape)
            selectModesSpp[Bn_1:Bn, Bn_1:Bn] = np.eye(Bn - Bn_1)
            return Tensor("selectModes({})".format(n), fullShape, spp=selectModesSpp)

        # Compute a matrix, which slices out one quantity
        #
        # @param quantityNumber The number of the quantity, which is sliced out
        def selectQuantity(quantityNumber):
            selectSpp = np.zeros((self.numberOfQuantities(), self.numberOfQuantities()))
            selectSpp[quantityNumber, quantityNumber] = 1
            return Tensor(
                "selectQuantity({})".format(quantityNumber),
                selectSpp.shape,
                spp=selectSpp,
            )

        # Compute a matrix, which slides out one quantity from the upper triangular
        #  part of the source matrix
        #
        # Note: G = E - diag(E)
        #
        # @param quantityNumber The number of the quantity, which is sliced out
        def selectQuantityG(quantityNumber):
            selectSpp = np.zeros((self.numberOfQuantities(), self.numberOfQuantities()))
            # The source matrix G only contains values at (o-4, o)
            selectSpp[quantityNumber - 4, quantityNumber] = 1
            return Tensor(
                "selectQuantityG({})".format(quantityNumber),
                selectSpp.shape,
                spp=selectSpp,
            )

        # Zinv(o) = $(Z - E^*_{oo} * I)^{-1}$
        #
        # @param o Index as described above
        def Zinv(o):
            return Tensor("Zinv({})".format(o), (self.order, self.order))

        # Submatrix of the stiffness matrix
        #
        # The stiffness matrix for derivatives in direction d, where only the
        # contribution from basis functions of degree n is considered
        #
        # @param d The direction in which derivatives are taken 0 =< d < 3
        # @param n The polynomial degree n, which is taken into account
        def kSub(d, n):
            Bn_1, Bn = modeRange(n)
            stiffnessSpp = np.zeros(fullShape)
            stiffnessSpp[:, Bn_1:Bn] = -stiffnessValues[d][:, Bn_1:Bn]
            return Tensor("kDivMTSub({},{})".format(d, n), fullShape, spp=stiffnessSpp)

        kernels = list()

        kernels.append(spaceTimePredictorRhs["kpt"] <= self.Q["kp"] * self.db.wHat["t"])
        for n in range(self.order - 1, -1, -1):
            for o in range(self.numberOfQuantities() - 1, -1, -1):
                kernels.append(
                    spaceTimePredictor["kpt"]
                    <= spaceTimePredictor["kpt"]
                    + selectModes(n)["kl"]
                    * selectQuantity(o)["pq"]
                    * spaceTimePredictorRhs["lqu"]
                    * Zinv(o)["ut"]
                )
                # G only has one relevant non-zero entry in each iteration, so we make it a scalar
                # G[o] = E[o-4, o] * timestep
                # In addition E only has non-zero entries, if o > 10
                if o >= 10:
                    kernels.append(
                        spaceTimePredictorRhs["kpt"]
                        <= spaceTimePredictorRhs["kpt"]
                        + G[o]
                        * selectQuantityG(o)["pv"]
                        * selectQuantity(o)["vq"]
                        * selectModes(n)["kl"]
                        * spaceTimePredictor["lqt"]
                    )
            if n > 0:
                derivativeSum = spaceTimePredictorRhs["kpt"]
                for d in range(3):
                    derivativeSum += (
                        kSub(d, n)["kl"]
                        * spaceTimePredictor["lqt"]
                        * self.starMatrix(d)["qp"]
                    )
            kernels.append(spaceTimePredictorRhs["kpt"] <= derivativeSum)
        kernels.append(
            self.I["kp"] <= timestep * spaceTimePredictor["kpt"] * self.db.timeInt["t"]
        )

        generator.add("spaceTimePredictor", kernels)

        # Test to see if the kernel actually solves the system of equations
        # This part is not used in the time kernel, but for unit testing
        deltaSppLarge = np.eye(self.numberOfQuantities())
        deltaLarge = Tensor("deltaLarge", deltaSppLarge.shape, spp=deltaSppLarge)
        deltaSppSmall = np.eye(self.order)
        deltaSmall = Tensor("deltaSmall", deltaSppSmall.shape, spp=deltaSppSmall)
        minus = Scalar("minus")

        lhs = deltaLarge["oq"] * self.db.Z["uk"] * spaceTimePredictor["lqk"]
        lhs += (
            minus
            * self.sourceMatrix()["qo"]
            * deltaSmall["uk"]
            * spaceTimePredictor["lqk"]
        )
        generator.add("stpTestLhs", testLhs["lou"] <= lhs)

        rhs = self.Q["lo"] * self.db.wHat["u"]
        for d in range(3):
            rhs += (
                minus
                * self.starMatrix(d)["qo"]
                * self.db.kDivMT[d]["lm"]
                * spaceTimePredictor["mqu"]
            )
        generator.add("stpTestRhs", testRhs["lou"] <= rhs)

        QAtTimeSTP = OptionalDimTensor(
            "QAtTimeSTP",
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            self.Q.shape(),
            alignStride=True,
        )
        timeBasisFunctionsAtPoint = Tensor("timeBasisFunctionsAtPoint", (self.order,))
        evaluateDOFSAtTimeSTP = (
            QAtTimeSTP["kp"]
            <= spaceTimePredictor["kpt"] * timeBasisFunctionsAtPoint["t"]
        )
        generator.add("evaluateDOFSAtTimeSTP", evaluateDOFSAtTimeSTP)

    def add_include_tensors(self, include_tensors):
        super().add_include_tensors(include_tensors)
        include_tensors.add(self.db.Z)


EQUATION_CLASS = PoroelasticADERDG
