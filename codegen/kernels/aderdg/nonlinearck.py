# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Cauchy-Kovalevskaya predictor for materials whose flux is nonlinear in Q.

The predictor keeps the derivative recursion of :class:`LinearCK`, run on star
matrices that linearise the material about the cell average. What differs is
everything downstream of it: the flux is evaluated pointwise at the time
quadrature nodes rather than applied as a constant operator, and the faces are
coupled through a Rusanov flux instead of a Riemann solver. Those parts are not
in this class yet; they arrive with the constitutive law that defines them.
"""

from kernels.common import generate_kernel_name_prefix
from kernels.multsim import OptionalDimTensor
from yateto import Scalar, ops, simpleParameterSpace
from yateto.ast.node import Accumulate
from yateto.ast.transformer import DeduceIndices, EquivalentSparsityPattern

from .aderdg import ADERDGBase


class NonLinearCK(ADERDGBase):
    def name(self):
        return "nonlinearck"

    def sourceMatrix(self):
        return None

    def extendedQTensor(self):
        return self.Q

    def numExtendedQuantities(self):
        return self.numQuantities()

    def addInit(self, generator):
        super().addInit(generator)

        iniShape = (
            self.num3DQuadraturePoints(),
            self.numQuantities(),
        )
        iniCond = OptionalDimTensor(
            "iniCond",
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            iniShape,
            alignStride=True,
        )
        dofsQP = OptionalDimTensor(
            "dofsQP",
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            iniShape,
            alignStride=True,
        )

        generator.add(
            "projectIniCond",
            self.Q["kp"] <= self.db.projectQP[self.t("kl")] * iniCond["lp"],
        )
        generator.add(
            "evalAtQP",
            dofsQP["kp"] <= self.db.evalAtQP[self.t("kl")] * self.Q["lp"],
        )

        nodalShape = (
            self.num3DQuadraturePoints(),
            self.numQuantities(),
        )
        # The time-integrated stress, carried alongside the integrated state.
        # It is what lets a neighbour evaluate its half of the flux without the
        # other cell's material, and it is modal because that is the form it
        # travels and is stored in.
        self.sigmaI = OptionalDimTensor(
            "sigmaI",
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            (self.num3DBasisFunctions(), 6),
            alignStride=True,
        )

        self.QNodal = OptionalDimTensor(
            "QNodal",
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            nodalShape,
            alignStride=True,
        )

    def addConstitutive(self, generator, target, prefix):
        """Pointwise material response at the nodes of :attr:`QNodal`.

        Empty here: what the flux and the source terms look like is a
        property of the constitutive law, so the material fills this in.
        """

    def nodalTensor(self, name, columns=None):
        """A tensor over the nodes of :attr:`QNodal`, optionally with a
        second axis of ``columns`` entries."""
        shape = (self.num3DQuadraturePoints(),)
        if columns is not None:
            shape = shape + (columns,)
        return OptionalDimTensor(
            name,
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            shape,
            alignStride=True,
        )

    def faceTensor(self, name, columns=None):
        """A tensor over the nodes of one face."""
        shape = (self.num2DBasisFunctions(),)
        if columns is not None:
            shape = shape + (columns,)
        return OptionalDimTensor(
            name,
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            shape,
            alignStride=True,
        )

    def addFaceProjection(self, generator, target, prefix):
        """Face-nodal values of the time-integrated state, and the way back.

        The local side reads its own face directly. The neighbour is restricted
        to its face and re-parameterised to ours before being evaluated at the
        nodes, which is the chain LinearCK walks for the neighbouring flux, cut
        short of the lift. Evaluating nodally decouples the two halves: the
        neighbour projection no longer depends on which of our faces it lands
        on, so the family is twelve rather than forty-eight, and the lift that
        does depend on it is four on its own.
        """
        atFace = self.faceTensor("QAtFace", self.numQuantities())
        fromNeighbor = self.faceTensor("QAtFaceNeighbor", self.numQuantities())
        flux = self.faceTensor("fluxAtFace", self.numQuantities())

        generator.addFamily(
            f"{prefix}projectToFace",
            simpleParameterSpace(4),
            lambda i: atFace["kp"]
            <= self.db.V3mTo2nFace[i][self.t("kl")] * self.I["lp"],
            target=target,
        )
        generator.addFamily(
            f"{prefix}projectNeighborToFace",
            simpleParameterSpace(3, 4),
            lambda h, j: fromNeighbor["kp"]
            <= self.db.V2mTo2n[self.t("km")]
            * self.db.fP[h][self.t("mn")]
            * self.db.rT[j][self.t("nl")]
            * self.I["lp"],
            target=target,
        )
        generator.addFamily(
            f"{prefix}faceIntegral",
            simpleParameterSpace(4),
            lambda i: self.Q["kp"]
            <= self.Q["kp"] + self.db.project2nFaceTo3m[i]["kn"] * flux["np"],
            target=target,
        )

        self.QAtFace = atFace
        self.QAtFaceNeighbor = fromNeighbor
        self.fluxAtFace = flux

    def addFaceFlux(self, generator, target, prefix):
        """The numerical flux at the face nodes.

        Empty here: which flux couples two cells follows from the constitutive
        law, so the material fills this in. It runs after
        :meth:`addFaceProjection`, whose tensors it reads.
        """

    def addLocal(self, generator, targets):
        for target in targets:
            prefix = generate_kernel_name_prefix(target)
            generator.add(
                f"{prefix}convertToNodal",
                self.QNodal["lp"] <= self.db.evalAtQP[self.t("lk")] * self.Q["kp"],
                target=target,
            )
            generator.add(
                f"{prefix}convertToModal",
                self.Q["kp"] <= self.db.projectQP[self.t("kl")] * self.QNodal["lp"],
                target=target,
            )
            self.addConstitutive(generator, target, prefix)
            self.addFaceProjection(generator, target, prefix)
            self.addFaceFlux(generator, target, prefix)

    def addNeighbor(self, generator, targets):
        pass

    def addTime(self, generator, targets):
        powers = [Scalar(f"power({i})") for i in range(self.order)]
        for target in targets:
            name_prefix = generate_kernel_name_prefix(target)

            qShape = (
                self.num3DBasisFunctions(),
                self.numQuantities(),
            )
            dQ0 = OptionalDimTensor(
                "dQ(0)",
                self.Q.optName(),
                self.Q.optSize(),
                self.Q.optPos(),
                qShape,
                alignStride=True,
            )
            power = powers[0]

            dQ0True = self.Q if target == "gpu" else dQ0

            derivatives = [dQ0True]

            derivativeExpr = [self.I["kp"] <= power * dQ0True["kp"]]
            derivativeTaylorExpansion = power * dQ0["kp"]

            if target == "gpu":
                derivativeExpr += [dQ0["kp"] <= self.Q["kp"]]

            self.dQs = [dQ0]

            for i in range(1, self.order):
                power = powers[i]
                derivativeSum = Accumulate(ops.Add())
                if self.sourceMatrix():
                    derivativeSum += derivatives[-1]["kq"] * self.sourceMatrix()["qp"]
                for j in range(3):
                    derivativeSum += (
                        self.db.kDivMT[j][self.t("kl")]
                        * derivatives[-1]["lq"]
                        * self.starMatrix(j)["qp"]
                    )

                derivativeSum = DeduceIndices(self.Q["kp"].indices).visit(derivativeSum)
                derivativeSum = EquivalentSparsityPattern().visit(derivativeSum)
                dQ = OptionalDimTensor(
                    "dQ({})".format(i),
                    self.Q.optName(),
                    self.Q.optSize(),
                    self.Q.optPos(),
                    qShape,
                    spp=derivativeSum.eqspp(),
                    alignStride=True,
                )
                self.dQs.append(dQ)

                derivativeExpr += [
                    dQ["kp"] <= derivativeSum,
                    self.I["kp"] <= self.I["kp"] + power * dQ["kp"],
                ]
                derivativeTaylorExpansion += power * dQ["kp"]

                derivatives.append(dQ)

            derivativeTaylorExpansionExpr = self.I["kp"] <= derivativeTaylorExpansion
            generator.add(f"{name_prefix}derivative", derivativeExpr, target=target)
            generator.add(
                f"{name_prefix}derivativeTaylorExpansion",
                derivativeTaylorExpansionExpr,
                target=target,
            )

    def add_include_tensors(self, include_tensors):
        super().add_include_tensors(include_tensors)
        include_tensors.add(self.db.nodes2D)
        include_tensors.add(self.sigmaI)
