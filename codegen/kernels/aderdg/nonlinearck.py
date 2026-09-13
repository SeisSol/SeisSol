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
from yateto import Scalar, Tensor, ops, simpleParameterSpace
from yateto.ast.node import Accumulate
from yateto.ast.transformer import DeduceIndices, EquivalentSparsityPattern
from yateto.memory import CSCMemoryLayout

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
            # The dissipation is scaled with a wave speed, and a wave speed is
            # not something that may be accumulated: over two timesteps of a
            # neighbour, the sum of two speeds is not a speed. Its square
            # integrated over the step is, and so is the length of the step --
            # both are time integrals, so the accumulation of a coarser
            # cluster sums them the way it sums every other column, and the
            # face that reads them divides and takes the root.
            QuantityGroup("waveIntegral", QuantityKind.INVARIANT),
            QuantityGroup("interval", QuantityKind.INVARIANT),
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
                if block.group.name in {group.name for group in self.primaryGroups()}
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
        """The two scalars of the step are one number per cell each, so their
        columns carry the constant mode alone."""
        spp = np.ones(
            (self.num3DBasisFunctions(), self.numTransportQuantities()), dtype=bool
        )
        for name in ("waveIntegral", "interval"):
            start, _ = self.transportGroupSlice(name)
            spp[1:, start] = False
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

    def flux_solver_spp(self):
        """Sparsity of the face flux solvers.

        A flux solver here maps the transported quantities onto the rows of
        the state, so it is as wide as the transport layout. Its pattern is
        the flux's, plus the diagonal of the coupled quantities: that is where
        the dissipation sits, and it is not where the flux is -- a strain row
        is fed by a velocity and a velocity row by a stress, never by itself.
        """
        shape = (self.numTransportQuantities(), self.numQuantities())
        spp = np.zeros(shape, dtype=bool)
        for source, target in self.fluxPattern():
            spp[source, target] = True
        for column in range(self.transportStateExtent()):
            spp[column, column] = True
        return spp

    def fluxPattern(self):
        """(transport quantity, state row) pairs the face flux connects.

        Empty here: which quantity feeds which equation is the constitutive
        law's business.
        """
        return []

    def numTimeNodes(self):
        """Time nodes the step kernel samples at.

        One more than the basis has functions, because the rule includes the
        ends of the timestep: a predictor that carries an internal variable
        across the nodes cannot account for an interval that no node sits in,
        and a rule whose nodes are all interior leaves the first one out. With
        the ends, the extra node costs its own arithmetic and nothing in
        accuracy -- Lobatto with n nodes integrates as far as Gauss-Legendre
        with n-1.
        """
        return self.order + 1

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
        nodes = self.numTimeNodes()
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

    def addCellIntegral(self, generator, target, prefix):
        """What the cell adds to its own state from its own integrals.

        Empty here: which map that is follows from the flux, so the material
        fills it in.
        """

    def addFluxSolver(self, generator, target, prefix):
        """The face flux, as two matrices per face and a scalar per step.

        The numerical flux is the average of the two sides plus a jump scaled
        with a wave speed, and both halves are linear in what a cell
        transports. So each half is a matrix, and the two differ only in the
        sign of the dissipation:

            Aplus  = C + lambda * D,   Aminus = C - lambda * D

        C carries the face geometry and the material, and is built once. D is
        the dissipation, and it is neither: it is the identity on the
        quantities the two cells couple through, so it is a constant of the
        layout. What is left per timestep is a scalar and an addition over the
        pattern of D -- nine entries, against rebuilding a flux solver.
        """
        fluxScale = Scalar("fluxScale")
        lambdaMax = Scalar("lambdaMax")
        normal = Tensor("faceNormal", (3,))

        dissipation = np.zeros(self.flux_solver_spp().shape)
        for column in range(self.transportStateExtent()):
            dissipation[column, column] = 0.5
        self.fluxDissipation = Tensor(
            "fluxDissipation",
            dissipation.shape,
            dissipation,
            CSCMemoryLayout,
        )

        self.fluxConstant = Tensor(
            "fluxConstant", self.flux_solver_spp().shape, spp=self.flux_solver_spp()
        )

        generator.add(
            f"{prefix}damageFluxSolver",
            self.fluxSolverStatements(fluxScale, normal),
            target=target,
        )
        generator.add(
            f"{prefix}damageFluxDissipation",
            [
                self.AplusT["qp"]
                <= self.fluxConstant["qp"] + lambdaMax * self.fluxDissipation["qp"],
                self.AminusT["qp"]
                <= self.fluxConstant["qp"] - lambdaMax * self.fluxDissipation["qp"],
            ],
            target=target,
        )

        generator.addFamily(
            f"{prefix}damageLocalFlux",
            simpleParameterSpace(4),
            lambda i: self.Q["kp"]
            <= self.Q["kp"]
            + self.db.rDivM[i][self.t("km")]
            * self.db.fMrT[i][self.t("ml")]
            * self.I["lq"]
            * self.AplusT["qp"],
            target=target,
        )
        generator.addFamily(
            f"{prefix}damageNeighborFlux",
            simpleParameterSpace(3, 4, 4),
            lambda h, j, i: self.Q["kp"]
            <= self.Q["kp"]
            + self.db.rDivM[i][self.t("km")]
            * self.db.fP[h][self.t("mn")]
            * self.db.rT[j][self.t("nl")]
            * self.I["lq"]
            * self.AminusT["qp"],
            target=target,
        )

    def fluxSolverStatements(self, fluxScale, normal):
        """How the constant half of the flux solver is built from a face
        normal. Empty here: it is the flux, so the material writes it."""
        return []

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
            self.addCellIntegral(generator, target, prefix)
            self.addFluxSolver(generator, target, prefix)

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
