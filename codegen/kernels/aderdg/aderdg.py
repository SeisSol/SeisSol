# SPDX-FileCopyrightText: 2019 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff

from abc import ABC, abstractmethod

import numpy as np
from kernels.multsim import OptionalDimTensor
from kernels.quantities import (
    FaceRole,
    layout,
    role_offset,
    total_extent,
    traction_selector,
    velocity_selector,
    well_formed,
)
from yateto import Scalar, Tensor, simpleParameterSpace
from yateto.input import parseJSONMatrixFile, parseXMLMatrixFile
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
            "I", "s", multipleSimulations, 0, qShape, alignStride=True
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

        selectVelocitySpp = self.mapToVelocities()[:, :3]
        self.selectVelocity = Tensor(
            "selectVelocity",
            selectVelocitySpp.shape,
            selectVelocitySpp,
            CSCMemoryLayout,
        )

        self.selectTractionSpp = self.mapToTractions()[:, :3]
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
        shape = (
            self.numExtendedQuantities(),
            self.numExtendedQuantities(),
        )
        return np.ones(shape, dtype=bool)

    def transformation_inv_spp(self):
        return self.godunov_spp()

    def mapToVelocities(self):
        return self.extractVelocities().T

    def mapToTractions(self):
        return self.extractTractions().T

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
        return velocity_selector(self.quantityBlocks())

    def extractTractions(self):
        return traction_selector(self.quantityBlocks())

    @abstractmethod
    def numExtendedQuantities(self):
        pass

    @abstractmethod
    def extendedQTensor(self):
        pass

    @abstractmethod
    def starMatrix(self, dim):
        pass

    def addInit(self, generator):
        well_formed(self.quantityBlocks())
        well_formed(self.extendedBlocks(), self.numExtendedQuantities())

        flux_solver_spp = self.flux_solver_spp()
        # The correction shares the flux solver's sparsity: it is added to the
        # same operator, so it cannot be populated where AplusT/AminusT are
        # structurally zero.
        self.QcorrLocal = Tensor("QcorrLocal", flux_solver_spp.shape, spp=flux_solver_spp)
        self.QcorrNeighbor = Tensor(
            "QcorrNeighbor", flux_solver_spp.shape, spp=flux_solver_spp
        )

        fluxScale = Scalar("fluxScale")
        computeFluxSolverLocal = (
            self.AplusT["ij"]
            <= fluxScale
            * self.Tinv["ki"]
            * (self.QgodLocal["kq"] * self.starMatrix(0)["ql"] + self.QcorrLocal["kl"])
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
        include_tensors.add(self.db.samplingDirections)
        include_tensors.add(self.db.M2inv)
        include_tensors.add(self.db.ET)
