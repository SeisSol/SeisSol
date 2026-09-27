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

        # The predictor carries a time index the operator does not touch, so
        # the two places a nodal operator is formed at are that much wider.
        nodalValuesInTime = None
        nodalProductInTime = None
        if self.nodalMaterial:
            nodalValuesInTime = self.nodalTemporary(
                "nodalValuesInTime", self.nodalValuesShape + (self.order,)
            )
            nodalProductInTime = self.nodalTemporary(
                "nodalProductInTime", self.nodalProductShape + (self.order,)
            )

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
            (self.num3DBasisFunctions(), self.numQuantities()),
            alignStride=True,
        )
        timeBasisFunctionsAtPoint = Tensor("timeBasisFunctionsAtPoint", (self.order,))

        for target in targets:
            name_prefix = generate_kernel_name_prefix(target)

            # One entry per stiff row, as a family indexed by position: the
            # kernels then do not care how many rows a material declares.
            stiffRows = {
                quantity: (targetQuantity, index)
                for index, (quantity, targetQuantity) in enumerate(
                    self.stiffSourceRows()
                )
            }
            if target == "cpu":
                G = {q: Scalar(f"G({i})") for q, (_, i) in stiffRows.items()}
                OptTimestep = lambda x: x
            else:
                G = {q: Tensor(f"Gt({i})", ())[""] for q, (_, i) in stiffRows.items()}

                # needed due to a current Yateto bug not allowing e.g. (Gkt * timestep)
                OptTimestep = lambda x: x * timestep

            kernels = list()

            def quantitySolve(modes, accumulate):
                """The time system of one block of space modes, solved quantity
                by quantity.

                What the source puts on the diagonal is already in Zinv, and
                the rows it couples are substituted back through G, which is
                why the quantities run downwards: a row is solved before
                anything writes into it.
                """

                def block(expr):
                    return expr if modes is None else expr.subslice("k", *modes)

                statements = []
                for o in range(self.numQuantities() - 1, -1, -1):
                    solved = (
                        block(spaceTimePredictorRhs["kpu"]).subslice("p", o, o + 1)
                        * Zinv(o)["ut"]
                    )
                    if accumulate:
                        solved = (
                            block(spaceTimePredictor["kpt"]).subslice("p", o, o + 1)
                            + solved
                        )
                    statements.append(
                        block(spaceTimePredictor["kpt"]).subslice("p", o, o + 1)
                        <= solved
                    )
                    # G has one relevant non-zero entry per stiff row, so it is a
                    # scalar: G[o] = E[target, o] * timestep. Rows that are not
                    # stiff contribute nothing.
                    if o in stiffRows:
                        o2 = stiffRows[o][0]
                        statements.append(
                            block(spaceTimePredictorRhs["kpt"]).subslice(
                                "p", o2, o2 + 1
                            )
                            <= block(spaceTimePredictorRhs["kpt"]).subslice(
                                "p", o2, o2 + 1
                            )
                            + OptTimestep(
                                G[o]
                                * block(spaceTimePredictor["kpt"]).subslice(
                                    "p", o, o + 1
                                )
                            )
                        )
                return statements

            if self.nodalMaterial:
                # An operator that is constant over the cell lowers the degree,
                # and that is what lets the blocks be taken one degree at a
                # time: the derivative of the block just solved reaches only
                # blocks still to come. One that varies inside the cell raises
                # the degree by as much as it carries itself, so it reaches
                # every block, and the sweep has to become a fixed point.
                #
                # The two agree where they overlap: with a constant operator
                # each step below makes one more degree exact, so after as many
                # steps as there are degrees the iteration is the sweep. Where
                # the material varies, what is left after those steps is the
                # scheme's own truncation error, since each step carries a
                # factor of the timestep.
                for iteration in range(self.order):
                    kernels.append(
                        spaceTimePredictorRhs["kpt"] <= self.Q["kp"] * self.db.wHat["t"]
                    )
                    if iteration > 0:
                        kernels += self.nodalApply(
                            spaceTimePredictor,
                            spaceTimePredictorRhs,
                            self.db.kDivMT,
                            spectator="t",
                            temporaries=(nodalValuesInTime, nodalProductInTime),
                            accumulate=True,
                            scalar=timestep,
                        )
                    kernels += quantitySolve(None, accumulate=False)
            else:
                kernels.append(
                    spaceTimePredictorRhs["kpt"] <= self.Q["kp"] * self.db.wHat["t"]
                )
                for n in range(self.order - 1, -1, -1):
                    kernels += quantitySolve(modeRange(n), accumulate=True)
                    if n > 0:
                        derivativeSum = spaceTimePredictorRhs["kpt"]
                        for d in range(3):
                            derivativeSum += (
                                self.db.kDivMT[d]["kl"].subslice("l", *modeRange(n))
                                * spaceTimePredictor["lqt"].subslice("l", *modeRange(n))
                                * self.starMatrix(d)["qp"]
                                * timestep
                            )
                        kernels.append(spaceTimePredictorRhs["kpt"] <= derivativeSum)
            kernels.append(
                self.I["kp"]
                <= timestep * spaceTimePredictor["kpt"] * self.db.timeInt["t"]
            )

            generator.add(
                f"{name_prefix}spaceTimePredictor",
                self.starAssembly() + kernels,
                target=target,
            )

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
        # This part is not used in the time kernel, but for unit testing.
        # The matrices are operands here, not something the kernel assembles:
        # the point of the check is that whatever shape a cell carries its
        # operator in, the predictor solves the system those matrices state.
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

        # the derivative matrix carries the sign of the term it stands for, so
        # the scalar here is the timestep the predictor scales its operator by
        rhs = self.Q["lo"] * self.db.wHat["u"]
        for d in range(3):
            rhs += (
                timestep
                * self.starMatrixSetup(d)["qo"]
                * self.db.kDivMT[d]["lm"]
                * spaceTimePredictor["mqu"]
            )
        generator.add("stpTestRhs", testRhs["lou"] <= rhs)

    def add_include_tensors(self, include_tensors):
        super().add_include_tensors(include_tensors)
        include_tensors.add(self.db.Z)
