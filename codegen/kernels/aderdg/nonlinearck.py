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

from dataclasses import replace

import numpy as np
from kernels.common import generate_kernel_name_prefix
from kernels.multsim import OptionalDimTensor
from kernels.quantities import (
    FaceRole,
    QuantityGroup,
    QuantityKind,
    layout,
    total_extent,
)
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

    def transportGroups(self):
        """Groups of the transported integrals.

        The face-coupled groups of the state, then the stress and the
        dissipation coefficient. The stress takes the traction role here: it is
        the mechanical traction, and it is present rather than derived, so the
        group that stands in for it in the state layout gives the role up.

        Groups without a face role are left out. They evolve through source
        terms local to the cell, so no neighbour ever reads them, and their
        own time integral stays where it is computed.
        """
        coupled = [
            (
                replace(group, role=FaceRole.NONE)
                if group.role is FaceRole.TRACTION
                else group
            )
            for group in self.primaryGroups()
            if group.role is not FaceRole.NONE
        ]
        return coupled + [
            QuantityGroup("sigma", QuantityKind.SYM_TENSOR2, FaceRole.TRACTION),
            QuantityGroup("lambdaMax", QuantityKind.INVARIANT),
        ]

    def transportBlocks(self):
        return layout(self.transportGroups())

    def transportStateExtent(self):
        """Quantities the transported tensor shares with the state, and in the
        same order: the Taylor expansion writes exactly these."""
        return total_extent(
            [
                block
                for block in self.transportBlocks()
                if block.group.name != "sigma" and block.group.name != "lambdaMax"
            ]
        )

    def nodalMeanWeights(self):
        """Weights that average a nodal field over the cell.

        The constant row of the modal projection is that average already:
        it is the quadrature weights, normalised. Reading it off the matrix
        that is there anyway keeps one definition of where the nodes are.
        """
        projection = self.db.projectQP
        shape = projection.shape()
        values = np.zeros(shape)
        for entry, value in projection.values().items():
            values[entry] = float(value)
        row = values[:, 0] if shape[0] == self.num3DQuadraturePoints() else values[0, :]
        return row / row.sum()

    def transportGroupSlice(self, name):
        """``(start, stop)`` of a transport group along the quantity axis."""
        for block in self.transportBlocks():
            if block.group.name == name:
                return block.offset, block.offset + block.extent
        raise ValueError(f"no transport group named {name}")

    def transportSpp(self):
        """The dissipation coefficient is one number per cell, so its column
        carries the constant mode alone."""
        spp = np.ones(
            (self.num3DBasisFunctions(), self.numTransportQuantities()), dtype=bool
        )
        spp[1:, -1] = False
        return spp

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

        self.addStepTensors()

        nodalShape = (
            self.num3DQuadraturePoints(),
            self.numQuantities(),
        )

        # The time-integrated source of the internal variables. They carry no
        # flux, so this never leaves the cell; it is the only thing the
        # corrector needs besides the transported tensor.
        self.sourceI = OptionalDimTensor(
            "sourceI",
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            (self.num3DBasisFunctions(), self.numInternalVariables()),
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

    def numInternalVariables(self):
        """Quantities of the state that no face reads: the difference between
        the state and what the cell hands over."""
        return self.numQuantities() - self.transportStateExtent()

    def addStepTensors(self):
        """Tensors the step kernel needs. Filled in by the material."""

    def finishStatements(self):
        """What the step leaves behind once all nodes are done. Filled in by
        the material."""
        return []

    def addStep(self, generator, target, prefix):
        """One kernel for the whole nonlinear part of a timestep.

        The time nodes are written out rather than looped over: how many there
        are is fixed at generation time, and writing them out lets one set of
        temporaries serve all of them, so the scratch the kernel needs does not
        grow with the number of nodes. What the nodes do not share is the
        coefficients, and those are scalars the launch code sets.

        The state at a node comes from the stored expansion, the material
        response from the constitutive law, and the integrals accumulate here
        rather than in a pass of their own -- there is nothing to be gained
        from writing a nodal stress out only to read it back and weigh it.
        """
        nodes = self.order
        evaluate = [
            [Scalar(f"evaluate({q},{i})") for i in range(self.order)]
            for q in range(nodes)
        ]
        weights = [Scalar(f"weight({q})") for q in range(nodes)]
        march = [Scalar(f"march({q})") for q in range(nodes)]

        # The nodal state at a time node lives inside the kernel: the step is
        # the only thing that looks at it.
        self.nodalState = self.nodalTensor(
            "QNodalAtTime", self.numQuantities(), temporary=True
        )

        state = OptionalDimTensor(
            "QAtTime",
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            (self.num3DBasisFunctions(), self.numQuantities()),
            alignStride=True,
            temporary=True,
        )

        statements = []
        for q in range(nodes):
            expansion = evaluate[q][0] * self.dQs[0]["kp"]
            for i in range(1, self.order):
                expansion = expansion + evaluate[q][i] * self.dQs[i]["kp"]
            statements += [
                state["kp"] <= expansion,
                self.nodalState["lp"] <= self.db.evalAtQP[self.t("lk")] * state["kp"],
            ]
            statements += self.stepStatements(q, weights[q], march[q])

        statements += self.finishStatements()
        generator.add(f"{prefix}damageStep", statements, target=target)

    def addConstitutive(self, generator, target, prefix):
        """Pointwise material response at the nodes of :attr:`QNodal`.

        Empty here: what the flux and the source terms look like is a
        property of the constitutive law, so the material fills this in.
        """

    def nodalTensor(self, name, columns=None, temporary=False):
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
            temporary=temporary,
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
        atFace = self.faceTensor("QAtFace", self.numTransportQuantities())
        fromNeighbor = self.faceTensor("QAtFaceNeighbor", self.numTransportQuantities())
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

            # The transported tensor is wider than the state and shares its
            # first columns; the stress and the dissipation coefficient are
            # filled where they are computed.
            shared = self.transportStateExtent()

            def state(tensor):
                return tensor["kp"].subslice("p", 0, shared)

            derivativeExpr = [state(self.I) <= power * state(dQ0True)]
            derivativeTaylorExpansion = power * state(dQ0)

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
                    state(self.I) <= state(self.I) + power * state(dQ),
                ]
                derivativeTaylorExpansion += power * state(dQ)

                derivatives.append(dQ)

            derivativeTaylorExpansionExpr = state(self.I) <= derivativeTaylorExpansion
            generator.add(f"{name_prefix}derivative", derivativeExpr, target=target)
            generator.add(
                f"{name_prefix}derivativeTaylorExpansion",
                derivativeTaylorExpansionExpr,
                target=target,
            )
            self.addStep(generator, target, name_prefix)

    def add_include_tensors(self, include_tensors):
        super().add_include_tensors(include_tensors)
        include_tensors.add(self.db.nodes2D)
