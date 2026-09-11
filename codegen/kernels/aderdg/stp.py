# SPDX-FileCopyrightText: 2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

import numpy as np
from kernels.common import generate_kernel_name_prefix
from kernels.multsim import OptionalDimTensor
from yateto import Scalar, Tensor

from .linearck import LinearCK


def choose(n, k):
    num = np.prod(np.arange(n, n - k, -1))
    denom = np.prod(np.arange(1, k + 1))
    return num // denom


class STP(LinearCK):
    """
    Space-time predictor for ADER-DG. The volume and flux kernels
    are the same as in the LinearCK case.

    The stiff source rows are factorised separately and substituted back
    through G; which rows those are comes from the equation, not from here.
    """

    def __init__(
        self,
        order,
        multipleSimulations,
        matricesDir,
        memLayout,
        numMechanisms,
        **kwargs,
    ):

        super().__init__(order, multipleSimulations, matricesDir)
        self.configure(
            matricesDir, memLayout, kwargs, extra=[f"{matricesDir}/stp_{order}.json"]
        )

    def numExtendedQuantities(self):
        return self.numQuantities()

    def starMatrix(self, dim):
        return self.db.star[dim]

    def sourceMatrix(self):
        return None

    def name(self):
        return "stp"

    def addTime(self, generator, targets):
        super().addTime(generator, targets)

        stpShape = (
            self.num3DBasisFunctions(),
            self.numQuantities(),
            self.order,
        )
        spaceTimePredictorRhs = OptionalDimTensor(
            "spaceTimePredictorRhs",
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            stpShape,
            alignStride=True,
            temporary=True,
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

            stiffRows = {
                q: (target_q, name) for q, target_q, name in self.stiffSourceRows()
            }
            if target == "cpu":
                G = {q: Scalar(name) for q, (_, name) in stiffRows.items()}
                OptTimestep = lambda x: x
            else:
                G = {
                    q: Tensor(f"{name}t", ())[""] for q, (_, name) in stiffRows.items()
                }

                # needed due to a current Yateto bug not allowing e.g. (Gkt * timestep)
                OptTimestep = lambda x: x * timestep

            kernels = list()

            kernels.append(
                spaceTimePredictorRhs["kpt"] <= self.Q["kp"] * self.db.wHat["t"]
            )
            for n in range(self.order - 1, -1, -1):
                for o in range(self.numQuantities() - 1, -1, -1):
                    kernels.append(
                        spaceTimePredictor["kpt"]
                        .subslice("k", *modeRange(n))
                        .subslice("p", o, o + 1)
                        <= spaceTimePredictor["kpt"]
                        .subslice("k", *modeRange(n))
                        .subslice("p", o, o + 1)
                        + spaceTimePredictorRhs["kpu"]
                        .subslice("k", *modeRange(n))
                        .subslice("p", o, o + 1)
                        * Zinv(o)["ut"]
                    )
                    # G has one relevant non-zero entry per stiff row, so it is a
                    # scalar: G[o] = E[target, o] * timestep. Rows that are not
                    # stiff contribute nothing.
                    if o in stiffRows:
                        o2 = stiffRows[o][0]
                        kernels.append(
                            spaceTimePredictorRhs["kpt"]
                            .subslice("k", *modeRange(n))
                            .subslice("p", o2, o2 + 1)
                            <= spaceTimePredictorRhs["kpt"]
                            .subslice("k", *modeRange(n))
                            .subslice("p", o2, o2 + 1)
                            + OptTimestep(
                                G[o]
                                * spaceTimePredictor["kpt"]
                                .subslice("k", *modeRange(n))
                                .subslice("p", o, o + 1)
                            )
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
                            self.db.kDivMT[d]["kl"].subslice("l", *modeRange(n))
                            * spaceTimePredictor["lqt"].subslice("l", *modeRange(n))
                            * star(d)
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
        deltaSppLarge = np.eye(self.numQuantities())
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
