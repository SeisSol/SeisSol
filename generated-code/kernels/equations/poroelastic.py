# SPDX-FileCopyrightText: 2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

import numpy as np
from kernels.aderdg import LinearADERDG
from kernels.common import generate_kernel_name_prefix
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
        **kwargs,
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
            return Tensor("kDivMTSub({},{})".format(d, n), fullShape, spp=stiffnessSpp, alignStride=True)

        QAtTimeSTP = OptionalDimTensor(
            "QAtTimeSTP",
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            self.Q.shape(),
            alignStride=True,
        )
        timeBasisFunctionsAtPoint = Tensor("timeBasisFunctionsAtPoint", (self.order,))

        for target in targets:
            name_prefix = generate_kernel_name_prefix(target)

            if target == "cpu":
                G = {10: Scalar("Gk"), 11: Scalar("Gl"), 12: Scalar("Gm")}
            else:
                G = {
                    10: Tensor("Gkt", ())[""] * timestep,
                    11: Tensor("Glt", ())[""] * timestep,
                    12: Tensor("Gmt", ())[""] * timestep,
                }

            kernels = list()

            kernels.append(
                spaceTimePredictorRhs["kpt"] <= self.Q["kp"] * self.db.wHat["t"]
            )
            for n in range(self.order - 1, -1, -1):
                for o in range(self.numberOfQuantities() - 1, -1, -1):
                    kernels.append(
                        spaceTimePredictor["kpt"].subslice('k', *modeRange(n)).subslice('p', o, o+1)
                        <= spaceTimePredictor["kpt"].subslice('k', *modeRange(n)).subslice('p', o, o+1)
                        + spaceTimePredictorRhs["kpu"].subslice('k', *modeRange(n)).subslice('p', o, o+1)
                        * Zinv(o)["ut"]
                    )
                    # G only has one relevant non-zero entry in each iteration, so we make it a scalar
                    # G[o] = E[o-4, o] * timestep
                    # In addition E only has non-zero entries, if o > 10
                    if o >= 10:
                        o2 = o-4
                        kernels.append(
                            spaceTimePredictorRhs["kpt"].subslice('k', *modeRange(n)).subslice('p', o2, o2+1)
                            <= spaceTimePredictorRhs["kpt"].subslice('k', *modeRange(n)).subslice('p', o2, o2+1)
                            + G[o]
                            * spaceTimePredictor["kpt"].subslice('k', *modeRange(n)).subslice('p', o, o+1)
                        )
                if n > 0:
                    derivativeSum = spaceTimePredictorRhs["kpt"]
                    star = (
                        (lambda d: self.starMatrix(d)["qp"] * timestep)
                        if target == "gpu"
                        else lambda d: self.starMatrix(d)["qp"]
                    )
                    for d in range(3):
                        derivativeSum += (
                            kSub(d, n)["kl"] * spaceTimePredictor["lqt"] * star(d)
                        )
                kernels.append(spaceTimePredictorRhs["kpt"] <= derivativeSum)
            kernels.append(
                self.I["kp"]
                <= timestep * spaceTimePredictor["kpt"] * self.db.timeInt["t"]
            )

            generator.add(f"{name_prefix}spaceTimePredictor", kernels, target=target)

            evaluateDOFSAtTimeSTP = (
                QAtTimeSTP["kp"]
                <= spaceTimePredictor["kpt"] * timeBasisFunctionsAtPoint["t"]
            )
            generator.add(
                f"{name_prefix}evaluateDOFSAtTimeSTP",
                evaluateDOFSAtTimeSTP,
                target=target,
            )

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

    def add_include_tensors(self, include_tensors):
        super().add_include_tensors(include_tensors)
        include_tensors.add(self.db.Z)


EQUATION_CLASS = PoroelasticADERDG
