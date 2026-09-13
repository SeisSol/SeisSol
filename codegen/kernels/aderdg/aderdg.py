# SPDX-FileCopyrightText: 2019 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff

from abc import ABC, abstractmethod

import numpy as np
from kernels.common import generate_kernel_name_prefix
from kernels.multsim import OptionalDimTensor
from kernels.quantities import (
    FaceRole,
    layout,
    role_offset,
    rotation_spp,
    total_extent,
    traction_selector,
    velocity_selector,
    well_formed,
)
from yateto import Scalar, Tensor, simpleParameterSpace
from yateto.input import (
    memoryLayoutFromFile,
    parseJSONMatrixFile,
    parseXMLMatrixFile,
)
from yateto.memory import CSCMemoryLayout
from yateto.util import (
    tensor_collection_from_constant_expression,
    tensor_from_constant_expression,
)


class ADERDGBase(ABC):
    def __init__(self, order, multipleSimulations, matricesDir):
        self.order = order

        self.alignStride = lambda name: True
        self.multipleSimulations = multipleSimulations
        if multipleSimulations > 1:
            self.alignStride = lambda name: name.startswith("fP")
        transpose = multipleSimulations > 1
        self.transpose = lambda name: transpose
        self.t = (lambda x: x[::-1]) if transpose else (lambda x: x)

        self.db = parseXMLMatrixFile(
            f"{matricesDir}/aderdg-{order}.xml",
            transpose=self.transpose,
            alignStride=self.alignStride,
        )
        clonesQP = {"v": ["evalAtQP"], "vInv": ["projectQP"]}
        self.db.update(
            parseJSONMatrixFile(
                f"{matricesDir}/plasticity-ip-matrices-{order}.json",
                clonesQP,
                transpose=self.transpose,
                alignStride=self.alignStride,
            )
        )
        self.db.update(parseJSONMatrixFile(f"{matricesDir}/sampling_directions.json"))
        self.db.update(parseJSONMatrixFile(f"{matricesDir}/mass-{order}.json"))

        # mass matrices are diagonal; treat them as sparse for now
        self.db.M2.setMemoryLayout(CSCMemoryLayout)
        self.db.M3.setMemoryLayout(CSCMemoryLayout)

        qShape = (self.num3DBasisFunctions(), self.numQuantities())
        self.Q = OptionalDimTensor(
            "Q", "s", multipleSimulations, 0, qShape, alignStride=True
        )

        self.I = OptionalDimTensor(
            "I",
            "s",
            multipleSimulations,
            0,
            (self.num3DBasisFunctions(), self.numTransportQuantities()),
            spp=self.transportSpp(),
            alignStride=True,
        )

        # What a cell hands a coarser neighbour, which is the same tensor one
        # cluster step further out. A second name for the same shape, because
        # the accumulation reads one and writes the other.
        self.IAccumulated = OptionalDimTensor(
            "IAccumulated",
            "s",
            multipleSimulations,
            0,
            (self.num3DBasisFunctions(), self.numTransportQuantities()),
            spp=self.transportSpp(),
            alignStride=True,
        )

        Aplusminus_spp = self.flux_solver_spp()
        self.AplusT = Tensor("AplusT", Aplusminus_spp.shape, spp=Aplusminus_spp)
        self.AplusTAll = [
            Tensor(f"AplusTAll({i})", Aplusminus_spp.shape, spp=Aplusminus_spp)
            for i in range(4)
        ]
        self.AminusT = Tensor("AminusT", Aplusminus_spp.shape, spp=Aplusminus_spp)
        trans_spp = self.transformation_spp()
        self.T = Tensor("T", trans_spp.shape, spp=trans_spp)
        trans_inv_spp = self.transformation_inv_spp()
        self.Tinv = Tensor("Tinv", trans_inv_spp.shape, spp=trans_inv_spp)
        godunov_spp = self.godunov_spp()
        self.QgodLocal = Tensor("QgodLocal", godunov_spp.shape, spp=godunov_spp)
        self.QgodNeighbor = Tensor("QgodNeighbor", godunov_spp.shape, spp=godunov_spp)

        self.oneSimToMultSim = Tensor(
            "oneSimToMultSim",
            (self.Q.optSize(),),
            spp={(i,): "1.0" for i in range(self.Q.optSize())},
        )

        self.db.update(
            parseJSONMatrixFile(
                "{}/nodal/nodalBoundary_matrices_{}.json".format(
                    matricesDir, self.order
                ),
                {},
                alignStride=self.alignStride,
                transpose=self.transpose,
                namespace="nodal",
            )
        )
        self.db.update(
            parseXMLMatrixFile(
                f"{matricesDir}/nodal/gravitational_energy_matrices_{self.order}.xml",
                alignStride=self.alignStride,
            )
        )

        # Note: MV2nTo2m is Vandermonde matrix from nodal
        # to modal representation WITHOUT mass matrix factor
        self.V2nTo2JacobiQuad = tensor_from_constant_expression(
            "V2nTo2JacobiQuad",
            self.db.V2mTo2JacobiQuad["ik"] * self.db.MV2nTo2m["kj"],
            target_indices="ij",
        )

        self.INodal = OptionalDimTensor(
            "INodal",
            "s",
            multipleSimulations,
            0,
            (self.num2DBasisFunctions(), self.numQuantities()),
            alignStride=True,
        )

        project2nFaceTo3m = tensor_collection_from_constant_expression(
            base_name="project2nFaceTo3m",
            expressions=lambda i: self.db.rDivM[i][self.t("jk")]
            * self.db.V2nTo2m["kl"],
            group_indices=simpleParameterSpace(4),
            target_indices="jl",
        )

        self.db.update(project2nFaceTo3m)

        # Into and out of the state: the free-surface displacement integrates
        # the velocity of the state itself, so this one stays with the
        # quantity layout.
        selectVelocitySpp = self.mapToVelocities()[:, :3]
        self.selectVelocity = Tensor(
            "selectVelocity",
            selectVelocitySpp.shape,
            selectVelocitySpp,
            CSCMemoryLayout,
        )

        # The same selection out of what crossed a face rather than out of the
        # state. One object where the two layouts agree, so that nothing about
        # a solver with a linear flux changes -- not even a name in the pool.
        if self.transportMatchesQuantities():
            self.selectVelocityTransported = self.selectVelocity
        else:
            transportedSpp = self.extractVelocities().T[:, :3]
            self.selectVelocityTransported = Tensor(
                "selectVelocityTransported",
                transportedSpp.shape,
                transportedSpp,
                CSCMemoryLayout,
            )

        self.selectTractionSpp = self.extractTractions().T[:, :3]
        self.tractionPlusMatrix = Tensor(
            "tractionPlusMatrix",
            self.selectTractionSpp.shape,
            self.selectTractionSpp,
            CSCMemoryLayout,
        )
        self.tractionMinusMatrix = Tensor(
            "tractionMinusMatrix",
            self.selectTractionSpp.shape,
            self.selectTractionSpp,
            CSCMemoryLayout,
        )

        # add an empty source matrix so that `ET` as name is defined
        if not self.db.containsName("ET"):
            self.db.ET = Tensor("ET", self.godunov_spp().shape)

    def name(self):
        return ""

    def num2DBasisFunctions(self):
        return self.order * (self.order + 1) // 2

    def num3DBasisFunctions(self):
        return self.order * (self.order + 1) * (self.order + 2) // 6

    def num3DQuadraturePoints(self):
        return (self.order + 1) ** 3

    def godunov_spp(self):
        shape = (self.numQuantities(), self.numQuantities())
        return np.ones(shape, dtype=bool)

    def flux_solver_spp(self):
        shape = (self.numQuantities(), self.numExtendedQuantities())
        return np.ones(shape, dtype=bool)

    def transformation_spp(self):
        return rotation_spp(self.extendedBlocks())

    def transformation_inv_spp(self):
        return rotation_spp(self.inverseRotationBlocks())

    def fusedInterpolationStatements(self, coeffs, extraCoeffs):
        """Time evaluation of everything a face reads, as statements leaving
        the result in `I`.

        Empty here: what a face reads is the state, so the sum over the state's
        own expansion is the whole of it and the rupture module writes it
        inline, without a tensor in between.
        """
        return []

    def drFluxSolverStatements(self, fluxScale, fluxSolver):
        """How a fault's flux solver is built, where the star matrix does not
        build it. Empty here: the Godunov flux of a state is the star matrix
        rotated, and that is what the rupture module writes."""
        return []

    def accumulateStatements(self):
        """How one step's integrals fold into what a coarser cluster has.

        Every column is a time integral, so every column is a sum.
        """
        return [self.IAccumulated["kp"] <= self.IAccumulated["kp"] + self.I["kp"]]

    def addAccumulate(self, generator, targets):
        """The accumulation of a step into a coarser cluster's buffer.

        A kernel rather than a loop so that the two paths run the same
        arithmetic, and so that a solver whose columns are not all sums says
        so once, here, instead of in each of them.
        """
        for target in targets:
            prefix = generate_kernel_name_prefix(target)
            generator.add(
                f"{prefix}accumulateIntegrals",
                self.accumulateStatements(),
                target=target,
            )

    def addStateToTransport(self, generator, targets):
        """What a cell transports, at one instant, from its state at that
        instant.

        Not a timestep: nothing marches and nothing is integrated. Where the
        two layouts coincide it is a copy, and that is every solver whose flux
        is linear. Where they do not, the transported tensor holds quantities
        that are functions of the state -- a stress, most of all -- and the
        material says how they are formed.
        """
        for target in targets:
            prefix = generate_kernel_name_prefix(target)
            generator.add(
                f"{prefix}stateToTransport",
                self.I["kp"] <= self.Q["kp"],
                target=target,
            )

    def transportTinv(self):
        """The inverse rotation a face applies to what crosses it.

        Always a tensor of its own, so that the kernels and the code binding
        them name one thing. Its pattern is the rotation of the state
        wherever the two layouts coincide, which is every solver whose flux
        is linear.
        """
        if not hasattr(self, "_transportTinv"):
            spp = self.transportTransformationInvSpp()
            self._transportTinv = Tensor("transportTinv", spp.shape, spp=spp)
        return self._transportTinv

    def transportTransformationInvSpp(self):
        """Inverse rotation over what a cell transports.

        A face rotates what crosses it, and what crosses it is the transported
        tensor rather than the state. Where the two coincide this is the
        rotation of the state, which is every solver whose flux is linear.
        """
        return rotation_spp(self.transportBlocks())

    #: The three directional star matrices share one sparsity pattern.
    StarClones = {"star": ["star(0)", "star(1)", "star(2)"]}

    def readMatrices(self, matricesDir, clones):
        """Reads this equation's matrix file."""
        return parseJSONMatrixFile(f"{matricesDir}/equation-{self.name()}.json", clones)

    def finishConfigure(self, memLayout, clones, kwargs):
        """Resolves the memory layout and stores the generator arguments. Split
        out so that an equation can reshape its matrices in between."""
        memoryLayoutFromFile(memLayout, self.db, clones)
        self.kwargs = kwargs

    def configure(self, matricesDir, memLayout, kwargs, extra=()):
        """Reads this equation's matrix file, plus any the solver needs, and
        resolves the memory layout across all of them."""
        clones = dict(self.StarClones)
        self.db.update(self.readMatrices(matricesDir, clones))
        for path in extra:
            self.db.update(parseJSONMatrixFile(path, clones))
        self.finishConfigure(memLayout, clones, kwargs)
        return clones

    def starMatrix(self, dim):
        return self.db.star[dim]

    def stiffSourceRows(self):
        """Source rows a space-time predictor has to factorise separately, as
        (quantity, target, scalar name). Empty unless the source term is
        stiff."""
        return []

    def mapToVelocities(self):
        """Into the state: where a fault's imposed velocity is added."""
        return velocity_selector(self.quantityBlocks()).T

    def mapToTractions(self):
        """Into the state: where a fault's imposed traction is added."""
        return traction_selector(self.quantityBlocks()).T

    @abstractmethod
    def primaryGroups(self):
        """Quantity groups of the underlying material, in quantity order."""

    def mechanismGroups(self):
        """Quantity groups of a single relaxation mechanism, if any."""
        return []

    def mechanismRepetitions(self):
        """How often the mechanism block sits on the quantity axis of Q."""
        return 0

    def quantityBlocks(self):
        """Layout of Q."""
        return layout(
            self.primaryGroups(), self.mechanismGroups(), self.mechanismRepetitions()
        )

    def extendedBlocks(self):
        """Layout the face rotation operates on."""
        return self.quantityBlocks()

    def transportBoundColumn(self):
        """Which column carries the bound a face scales its dissipation with.

        The number of transported quantities where there is none, which is
        every solver whose face reads a wave speed out of the material.
        """
        return self.numTransportQuantities()

    def transportInternalOffset(self):
        """Where the groups of the state that carry no flux begin, for a solver
        that transports them anyway. Zero where there are none."""
        leading = self.transportStateExtent()
        names = {group.name for group in self.primaryGroups()}
        for block in self.transportBlocks():
            if block.offset >= leading and block.group.name in names:
                return block.offset
        return 0

    def transportBlocks(self):
        """Layout of the time-integrated quantities, i.e. of :attr:`I`.

        A solver whose flux is nonlinear in the state carries more across a
        face than the state itself, because the integral of a nonlinear flux is
        not the flux of the integrated state. What it carries is described
        here; for everyone else the two coincide.
        """
        return self.quantityBlocks()

    def numTransportQuantities(self):
        return total_extent(self.transportBlocks())

    def materialParameterNames(self):
        """Material parameters the kernels read per element, in the order they
        sit in the tensor that carries them. Empty where a material's
        parameters reach the kernels as scalars."""
        return []

    def transportStateExtent(self):
        """Quantities the transported tensor shares with the state, and in the
        same order. Where the two layouts coincide, that is all of them."""
        return self.numQuantities()

    def transportSpp(self):
        """Sparsity of :attr:`I`. Dense unless a solver says otherwise."""
        return None

    def transportMatchesQuantities(self):
        """Whether :attr:`I` is laid out like :attr:`Q`.

        Where it is not, a tensor shaped by the quantity layout cannot be
        contracted with the integrals, and the kernels that do so are not
        available.
        """
        return self.numTransportQuantities() == self.numQuantities()

    def inverseRotationBlocks(self):
        """Layout the inverse face rotation operates on. It need not match the
        forward one: a solver keeping the mechanism index in its own tensor
        dimension rotates one anelastic block forwards and none back."""
        return self.extendedBlocks()

    def numQuantities(self):
        return total_extent(self.quantityBlocks())

    def velocityOffset(self):
        return role_offset(self.quantityBlocks(), FaceRole.VELOCITY)

    def extractVelocities(self):
        """Out of what crossed a face, which is the transported tensor -- the
        state wherever the two coincide."""
        return velocity_selector(self.transportBlocks())

    def extractTractions(self):
        """Out of what crossed a face. See extractVelocities."""
        return traction_selector(self.transportBlocks())

    @abstractmethod
    def numExtendedQuantities(self):
        pass

    @abstractmethod
    def extendedQTensor(self):
        pass

    def addInit(self, generator):
        well_formed(self.quantityBlocks())
        well_formed(self.extendedBlocks(), self.numExtendedQuantities())

        flux_solver_spp = self.flux_solver_spp()
        # The correction shares the flux solver's sparsity: it is added to the
        # same operator, so it cannot be populated where AplusT/AminusT are
        # structurally zero.
        self.QcorrLocal = Tensor(
            "QcorrLocal", flux_solver_spp.shape, spp=flux_solver_spp
        )
        self.QcorrNeighbor = Tensor(
            "QcorrNeighbor", flux_solver_spp.shape, spp=flux_solver_spp
        )

        # A Godunov state is a state, so the solver it builds maps the state
        # onto itself. Where a cell transports more than that, the flux solver
        # is wider than the rotation and is assembled from the flux instead.
        if self.transportMatchesQuantities():
            fluxScale = Scalar("fluxScale")
            computeFluxSolverLocal = (
                self.AplusT["ij"]
                <= fluxScale
                * self.Tinv["ki"]
                * (
                    self.QgodLocal["kq"] * self.starMatrix(0)["ql"]
                    + self.QcorrLocal["kl"]
                )
                * self.T["jl"]
            )
            generator.add("computeFluxSolverLocal", computeFluxSolverLocal)

            computeFluxSolverNeighbor = (
                self.AminusT["ij"]
                <= fluxScale
                * self.Tinv["ki"]
                * (
                    self.QgodNeighbor["kq"] * self.starMatrix(0)["ql"]
                    + self.QcorrNeighbor["kl"]
                )
                * self.T["jl"]
            )
            generator.add("computeFluxSolverNeighbor", computeFluxSolverNeighbor)

        stiffnessTensor = Tensor("stiffnessTensor", (3, 3, 3, 3))
        direction = Tensor("direction", (3,))
        christoffel = Tensor("christoffel", (3, 3))

        computeChristoffel = (
            christoffel["ik"]
            <= stiffnessTensor["ijkl"] * direction["j"] * direction["l"]
        )
        generator.add("computeChristoffel", computeChristoffel)

        self.addEnergyProducts(generator)

    def addEnergyProducts(self, generator):
        """Mass-matrix moments of Q, used by the volume energy output.

        momentQ[0,J]  == \\int_{T_ref} Q_J
        momentQQ[I,J] == \\int_{T_ref} Q_I Q_J

        Multiply by the Jacobi determinant to obtain the physical integral. Both
        are exact, as opposed to evaluating at quadrature points.

        Note: this lives in ADERDGBase (not LinearCK), because the
        viscoelastic2 generator derives directly from ADERDGBase and would
        otherwise not get the kernels at all.
        """
        # Only the cell integral is needed, so M3 is narrowed to its first row.
        # subselect keeps the rank and sets the extent to 1, which turns the
        # kernel from an nb x nq product into an nq one.
        momentQ = OptionalDimTensor(
            "momentQ",
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            (1, self.numQuantities()),
            alignStride=True,
        )
        generator.add(
            "momentQCompute",
            momentQ["IJ"] <= self.db.M3["Ij"].subselect("I", 0) * self.Q["jJ"],
        )

        # The fused-simulation index 's' occurs in the result and in both Q
        # factors, i.e. it is a batch index. yateto handles that as of
        # <yateto batch-index fix>; without it this asserts in the GEMM factory.
        momentQQ = OptionalDimTensor(
            "momentQQ",
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            (self.numQuantities(), self.numQuantities()),
        )
        generator.add(
            "momentQQCompute",
            momentQQ["IJ"] <= self.db.M3["ij"] * self.Q["iI"] * self.Q["jJ"],
        )

    @abstractmethod
    def addLocal(self, generator, targets):
        pass

    @abstractmethod
    def addNeighbor(self, generator, targets):
        pass

    @abstractmethod
    def addTime(self, generator, targets):
        pass

    def add_include_tensors(self, include_tensors):
        if not self.transportMatchesQuantities():
            # The Godunov state of a solver that transports more than its
            # state is not assembled, so nothing pulls these in by use -- and
            # the setup of every material is compiled whatever the material.
            include_tensors.add(self.QgodLocal)
            include_tensors.add(self.QgodNeighbor)
            include_tensors.add(self.QcorrLocal)
            include_tensors.add(self.QcorrNeighbor)
        include_tensors.add(self.db.samplingDirections)
        include_tensors.add(self.db.M2inv)
        include_tensors.add(self.db.ET)
