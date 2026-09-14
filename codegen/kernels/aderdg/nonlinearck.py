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
import yateto.functions as yf
from kernels.common import generate_kernel_name_prefix
from kernels.multsim import OptionalDimTensor
from kernels.quantities import (
    FaceRole,
    QuantityGroup,
    QuantityKind,
    layout,
    rotation_spp,
    total_extent,
)
from yateto import Scalar, Tensor, ops, simpleParameterSpace
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
        return (
            coupled
            + [
                QuantityGroup("sigma", QuantityKind.SYM_TENSOR2, FaceRole.TRACTION),
                # What the state keeps to itself is transported after all: a face
                # needs the moduli a wave sees there, and those follow the
                # internal variables. They are carried rather than shared with the
                # state, because their own expansion is the projected one -- they
                # march through the step from a source, which no Taylor sum of the
                # state describes.
            ]
            + [
                replace(group, role=FaceRole.NONE)
                for group in self.primaryGroups()
                if group.role is FaceRole.NONE
            ]
            + [
                # The dissipation is scaled with a wave speed, and Rusanov wants
                # an upper bound on the instantaneous one. So what is carried is
                # the largest square of it over the step: the square, because
                # that is what the moduli are affine in, and the largest rather
                # than a mean, because a mean is below the bound at every
                # instant where the speed is above it. A maximum also survives
                # accumulation, which a mean does not -- a coarser cluster takes
                # the larger of what it has and what it reads.
                QuantityGroup("waveIntegral", QuantityKind.INVARIANT),
                QuantityGroup("shearIntegral", QuantityKind.INVARIANT),
            ]
        )

    def transportBlocks(self):
        return layout(self.transportGroups())

    def accumulateStatements(self):
        """Every column is a time integral and therefore a sum -- except the
        one that carries the bound a face scales its dissipation with, which
        is a maximum over the step. The larger of two bounds is a bound; their
        sum is not one, and it would grow with the cluster ratio."""
        bound = self.transportBoundColumn()
        last = bound + self.transportBoundCount()
        total = self.numTransportQuantities()
        statements = []
        if bound > 0:
            statements += [
                self.IAccumulated["kc"].subslice("c", 0, bound)
                <= self.IAccumulated["kc"].subslice("c", 0, bound)
                + self.I["kc"].subslice("c", 0, bound)
            ]
        statements += [
            self.IAccumulated["kc"].subslice("c", bound, last)
            <= yf.maximum(
                self.IAccumulated["kc"].subslice("c", bound, last),
                self.I["kc"].subslice("c", bound, last),
            )
        ]
        if last < total:
            statements += [
                self.IAccumulated["kc"].subslice("c", last, total)
                <= self.IAccumulated["kc"].subslice("c", last, total)
                + self.I["kc"].subslice("c", last, total)
            ]
        return statements

    def transportBoundColumn(self):
        """The first of the columns that carry a bound. They are adjacent, so
        whoever treats them differently from a time integral treats a run of
        columns rather than one and then another."""
        return self.transportGroupSlice("waveIntegral")[0]

    def transportBoundCount(self):
        """How many columns carry a bound: one per wave family the dissipation
        is scaled with separately."""
        return 2

    def transportStateExtent(self):
        """Quantities the transported tensor shares with the state, and in the
        same order: the Taylor expansion writes exactly these."""
        # The leading blocks only. A group of the state may appear further
        # back as well -- the internal variables are carried there, with an
        # expansion of their own -- and the Taylor sum writes a prefix, not a
        # selection.
        names = [group.name for group in self.primaryGroups()]
        prefix = []
        for block, name in zip(self.transportBlocks(), names):
            if block.group.name != name:
                break
            prefix.append(block)
        return total_extent(prefix)

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
        for name in ("waveIntegral", "shearIntegral"):
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

        # The far side of a face hands over a tensor of the same shape. It is
        # an argument of its own rather than a second use of I, because a face
        # reads both at once.
        self.INeighbor = OptionalDimTensor(
            "INeighbor",
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            (self.num3DBasisFunctions(), self.numTransportQuantities()),
            spp=self.transportSpp(),
            alignStride=True,
        )

        # The expansion in time of everything the transported tensor carries
        # beyond the state: the stress and the two scalars. The state has the
        # derivative family; this is its counterpart, and it is what a
        # neighbour on a coarser cluster reconstructs a subinterval from.
        carried = self.numTransportQuantities() - self.transportStateExtent()

        # The bounds a face scales its dissipation with are a maximum over the
        # step, and a maximum has no derivative: the first member carries them
        # and the rest have nowhere to put one. Said as a layout rather than
        # written as a zero -- a zero has to be assigned by somebody, and an
        # assignment of zero is the kind of statement that survives reading and
        # not optimisation.
        def carriedSpp(member):
            spp = np.ones((self.num3DBasisFunctions(), carried), dtype=bool)
            if member > 0:
                bound = self.transportBoundColumn() - self.transportStateExtent()
                spp[:, bound : bound + self.transportBoundCount()] = False
            return spp

        self.transportDer = [
            OptionalDimTensor(
                f"transportDer({i})",
                self.Q.optName(),
                self.Q.optSize(),
                self.Q.optPos(),
                (self.num3DBasisFunctions(), carried),
                spp=carriedSpp(i),
                alignStride=True,
            )
            for i in range(self.order)
        ]

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

        # A face without a neighbour folds its ghost rule into the pair, and a
        # ghost rule is a rotation: it mixes the rows of a group among
        # themselves. The pattern has to be closed under that, or the folding
        # would write where the layout says there is nothing.
        rotation = rotation_spp(self.transportBlocks())
        return spp | (rotation.T.astype(int) @ spp.astype(int)).astype(bool)

    def fluxPattern(self):
        """(transport quantity, state row) pairs the face flux connects.

        Empty here: which quantity feeds which equation is the constitutive
        law's business.
        """
        return []

    def timeNodes(self):
        """The nodes the step samples at, in time scaled onto [0, 1].

        Gauss-Lobatto with one node more than the basis has functions. The
        same rule the launch code asks its time basis for -- the two are
        derived independently and agree to 4e-16 at order six, which is the
        one coupling in this construction that no compiler checks.
        """
        points = self.numTimeNodes()
        if points == 2:
            reference = np.array([-1.0, 1.0])
        else:
            inner = np.polynomial.legendre.legroots(
                np.polynomial.legendre.legder([0] * (points - 1) + [1])
            )
            reference = np.concatenate(([-1.0], inner, [1.0]))
        return 0.5 * (reference + 1.0)

    def timeWeights(self):
        """Quadrature weights of :meth:`timeNodes`, in time scaled onto
        [0, 1].

        Gauss-Lobatto, so the weight of a node follows from the Legendre
        polynomial of the rule's degree there. They sum to one, because the
        interval they cover is the scaled step -- which is the same rule the
        launch code asks its time basis for, scaled by the timestep.
        """
        points = self.numTimeNodes()
        reference = 2.0 * self.timeNodes() - 1.0
        legendre = np.polynomial.legendre.legval(reference, [0] * (points - 1) + [1])
        return 1.0 / (points * (points - 1) * legendre * legendre)

    def timeProjection(self):
        """Coefficients of the shifted Legendre expansion, from the values at
        the time nodes.

        What a cell transports beyond its state has no recursion to come out
        of, so its expansion in time is won from the samples -- and won in a
        Legendre basis, because the monomial one loses five digits at order
        six and leaves two of them in single precision.

        There is one sample more than there are coefficients, so the fit is a
        least-squares one, and which inner product it minimises in decides
        what it answers. The one the nodes carry is the quadrature's: under it
        the coefficient of the constant is the rule's own integral, so
        reconstructing the whole step reproduces the integral the step
        accumulates. An unweighted fit minimises in an inner product the nodes
        do not have and lands on a different number -- at order two, the plain
        mean of three samples where the rule is Simpson's.

        The result is constant: the nodes are in scaled time, and so is the
        expansion, which is also the convention LegendreBasis::integrate
        reads.
        """
        nodes = self.timeNodes()
        vandermonde = np.array(
            [
                [
                    np.polynomial.legendre.legval(2.0 * node - 1.0, [0] * i + [1])
                    for i in range(self.order)
                ]
                for node in nodes
            ]
        )
        weighted = vandermonde.T * self.timeWeights()
        return np.linalg.solve(weighted @ vandermonde, weighted)

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

    def nodalTensor(self, name, columns=None, temporary=False, datatype=None):
        """A tensor over the nodes of :attr:`QNodal`, optionally with a
        second axis of ``columns`` entries.

        ``datatype`` for a temporary that does not hold a number: a tensor
        without one is the working precision, and a truth value written into
        it is then one type where it is written and another where it is read.
        """
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
            datatype=datatype,
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
        normal = Tensor("faceNormal", (3,))

        # Selectors for the two scalars a transported tensor carries. A wave
        # speed is a property of an element, and a scalar argument is uniform
        # over a batch, so the speed of a face is read out of the tensors and
        # formed inside the kernel that needs it.
        shape = (self.num3DBasisFunctions(), self.numTransportQuantities())

        def pick(name):
            values = np.zeros(shape)
            values[0, self.transportGroupSlice(name)[0]] = 1.0
            return values

        self.pickWave = Tensor("pickWaveIntegral", shape, pick("waveIntegral"))
        self.pickShear = Tensor("pickShearIntegral", shape, pick("shearIntegral"))

        # Both halves of the pair are per cell and face. The dissipation is
        # the identity on the coupled quantities for a face with a neighbour,
        # but a face without one has a ghost rule, and folding that rule into
        # the pair is what makes a boundary condition the same two matrices
        # and the same scalar as everything else.
        self.fluxConstant = Tensor(
            "fluxConstant", self.flux_solver_spp().shape, spp=self.flux_solver_spp()
        )
        # One per wave family. What is in them is what decides which flux this
        # is: the identity on the coupled quantities in the first and nothing
        # in the second is Rusanov, because the first bound is the larger; the
        # two rotated projectors of the face frame is an upwind flux. The
        # kernel never learns which, so there is no case to distinguish here.
        self.fluxDissipation = Tensor(
            "fluxDissipation", self.flux_solver_spp().shape, spp=self.flux_solver_spp()
        )
        self.fluxDissipationShear = Tensor(
            "fluxDissipationShear",
            self.flux_solver_spp().shape,
            spp=self.flux_solver_spp(),
        )

        self.ghostMap = Tensor(
            "ghostMap",
            (self.numTransportQuantities(), self.numTransportQuantities()),
            spp=rotation_spp(self.transportBlocks()),
        )
        self.fluxSource = Tensor(
            "fluxSource", self.flux_solver_spp().shape, spp=self.flux_solver_spp()
        )
        self.fluxFolded = Tensor(
            "fluxFolded", self.flux_solver_spp().shape, spp=self.flux_solver_spp()
        )
        ghostSign = Scalar("ghostSign")

        generator.add(
            f"{prefix}damageFluxGhost",
            self.fluxFolded["rp"]
            <= self.fluxSource["rp"]
            + ghostSign * self.ghostMap["qr"] * self.fluxSource["qp"],
            target=target,
        )

        generator.add(
            f"{prefix}damageFluxSolver",
            self.fluxSolverStatements(fluxScale, normal),
            target=target,
        )
        # The pair is assembled inside each half rather than between them: a
        # matrix that crosses a kernel boundary is a matrix the batched path
        # has to hold per element, and building it twice costs a hundred flops
        # against the thousands the half itself costs.
        # Temporaries of their own rather than the declared pair: the pair is
        # what a cell stores per face, and the C++ sizes those arrays with it,
        # so it has to stay a tensor even though no kernel takes it as an
        # argument any more.
        self.fluxPlus = Tensor(
            "fluxPlus",
            self.flux_solver_spp().shape,
            spp=self.flux_solver_spp(),
            temporary=True,
        )
        self.fluxMinus = Tensor(
            "fluxMinus",
            self.flux_solver_spp().shape,
            spp=self.flux_solver_spp(),
            temporary=True,
        )

        def dissipation(own, other):
            # The larger of the two sides' bounds, per wave family, on the
            # matrix that family is scaled with. Read out of the tensors
            # rather than passed in: a scalar argument is uniform over a batch
            # and a wave speed is not.
            def speed(selector):
                bound = lambda tensor: tensor["kc"] * selector["kc"]
                return yf.sqrt(yf.maximum(bound(own), bound(other)))

            return (
                speed(self.pickWave) * self.fluxDissipation["qp"]
                + speed(self.pickShear) * self.fluxDissipationShear["qp"]
            )

        generator.addFamily(
            f"{prefix}damageLocalFlux",
            simpleParameterSpace(4),
            lambda i: [
                self.fluxPlus["qp"]
                <= self.fluxConstant["qp"] + dissipation(self.I, self.INeighbor),
                self.Q["kp"]
                <= self.Q["kp"]
                + self.db.rDivM[i][self.t("km")]
                * self.db.fMrT[i][self.t("ml")]
                * self.I["lq"]
                * self.fluxPlus["qp"],
            ],
            target=target,
        )
        generator.addFamily(
            f"{prefix}damageNeighborFlux",
            simpleParameterSpace(3, 4, 4),
            lambda h, j, i: [
                self.fluxMinus["qp"]
                <= self.fluxConstant["qp"] - dissipation(self.I, self.INeighbor),
                self.Q["kp"]
                <= self.Q["kp"]
                + self.db.rDivM[i][self.t("km")]
                * self.db.fP[h][self.t("mn")]
                * self.db.rT[j][self.t("nl")]
                * self.INeighbor["lq"]
                * self.fluxMinus["qp"],
            ],
            target=target,
        )

    def addTransportToState(self, generator, targets):
        """The state of a cell, read back out of what it transports.

        Two blocks rather than one: the quantities a face couples through sit
        at the front of both layouts, and the rest of the state sits where the
        transport layout keeps the groups that carry no flux. Reading the
        transported tensor through the state's view would take the first
        component of the stress for the first internal variable.
        """
        shared = self.transportStateExtent()
        internal = self.transportInternalOffset()
        carried = self.numQuantities() - shared

        for target in targets:
            prefix = generate_kernel_name_prefix(target)
            generator.add(
                f"{prefix}transportToState",
                [
                    self.Q["kp"].subslice("p", 0, shared)
                    <= self.I["kp"].subslice("p", 0, shared),
                    self.Q["kp"].subslice("p", shared, self.numQuantities())
                    <= self.I["kp"].subslice("p", internal, internal + carried),
                ],
                target=target,
            )

    def fusedInterpolationStatements(self, coeffs, extraCoeffs):
        """Both sums, into the tensor a face then reads.

        The columns a face shares with the state come out of the state's
        expansion and its coefficients; the rest out of the expansion
        projected for them and the coefficients of that basis. Which is the
        same pair of sums the serial evaluation performs -- here they land in
        one tensor first, because what follows contracts it as a whole.
        """
        shared = self.transportStateExtent()
        carried = (shared, self.numTransportQuantities())

        state = Accumulate(ops.Add())
        for i, coefficient in enumerate(coeffs):
            state = state + coefficient * self.dQs[i]["lq"].subslice("q", 0, shared)

        rest = Accumulate(ops.Add())
        for i, coefficient in enumerate(extraCoeffs):
            rest = rest + coefficient * self.transportDer[i]["lq"]

        return [
            self.I["lq"].subslice("q", 0, shared) <= state,
            self.I["lq"].subslice("q", *carried) <= rest,
        ]

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

            # The other half of a reconstruction: the columns that are not the
            # state, out of the expansion that was projected for them. The
            # coefficients are of the other basis, and they are scalars
            # because an interval is a property of a cluster and not of a
            # cell.
            extraPowers = [Scalar(f"extraPower({i})") for i in range(self.order)]
            carried = (self.transportStateExtent(), self.numTransportQuantities())

            def rest(tensor):
                return tensor["kp"].subslice("p", *carried)

            carriedExpansion = Accumulate(ops.Add())
            for i in range(self.order):
                carriedExpansion += extraPowers[i] * self.transportDer[i]["kp"]
            generator.add(
                f"{name_prefix}carriedTaylorExpansion",
                rest(self.I) <= carriedExpansion,
                target=target,
            )
            self.addStep(generator, target, name_prefix)

    def add_include_tensors(self, include_tensors):
        super().add_include_tensors(include_tensors)
        include_tensors.add(self.db.nodes2D)
        # A cell stores the pair of a face and the C++ sizes those arrays with
        # it. No kernel of this solver takes it as an argument -- each half
        # assembles what it applies -- so nothing would pull it in by use.
        include_tensors.add(self.AplusT)
        include_tensors.add(self.AminusT)
