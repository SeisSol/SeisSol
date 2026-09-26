# SPDX-FileCopyrightText: 2019 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff

from abc import ABC, abstractmethod

import numpy as np
from kernels import coefficients, material
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
from yateto.type import AddressingMode
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


def negateFamily(db, baseName, alignStride):
    """Replaces a tensor family in `db` with its negation.

    The sparsity pattern is unaffected by a sign, so the rebuilt tensors carry
    the same pattern, shape and stride alignment and can take the place of the
    originals without anything downstream noticing. Values stay decimal text,
    the form the matrix files deliver them in and the form the emitter prints.
    """
    db[baseName] = {
        group: Tensor(
            tensor.name(),
            tensor.shape(),
            {index: repr(-float(value)) for index, value in tensor.values().items()},
            alignStride=alignStride(tensor.name()),
        )
        for group, tensor in db[baseName].items()
    }


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
        # The derivative kernels contract against -kDivMT throughout. Folding
        # the sign into the matrix here rather than scaling the operand keeps
        # the global matrices read-only, which is what lets them be shared and
        # placed wherever a target wants them. It has to happen before the
        # memory layout configuration is applied, since that mutates whatever
        # tensors are in the database at the time.
        negateFamily(self.db, "kDivMT", self.alignStride)

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

        # Which simulation a point source acts on is set per source at run
        # time, so this is a pattern, not numbers: numbers would make it a
        # constant, held in the pool and bound by bindGlobals.
        self.oneSimToMultSim = Tensor(
            "oneSimToMultSim",
            (self.Q.optSize(),),
            spp={(i,): True for i in range(self.Q.optSize())},
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

        # The traction weights are computed per fault face from the impedances
        # (DynamicRuptureMatrices), so only their pattern is known here. Passed
        # as booleans: a float array would be taken for the values, which
        # would make them constants held in the pool and bound by bindGlobals.
        self.selectTractionSpp = self.tractionMatrixSpp() != 0
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

        # The canonical vertex numbering forces the face orientation index to
        # zero on every interior face, so the neighbouring flux matrix is the
        # constant fP(0) and folds into the neighbour change of basis. The
        # stride alignment is given explicitly, since the name based rule would
        # key off the "fP" prefix, while this tensor takes the role, and hence
        # the alignment, of a change of basis matrix.
        self.db.update(
            tensor_collection_from_constant_expression(
                "fPrT",
                lambda j: self.db.fP[0][self.t("mn")] * self.db.rT[j][self.t("nl")],
                simpleParameterSpace(4),
                target_indices=self.t("ml"),
                tensor_args={"alignStride": self.multipleSimulations == 1},
                zero_tolerance=1e-14,
            )
        )

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
        self._configureStarAssembly(kwargs)

    def _configureStarAssembly(self, kwargs):
        """Sets up the tensors a cell carries where it holds the coefficients
        of its operator rather than the matrices they fold into.

        The structure the two fold into is a signed permutation, so it is
        stated as an immediate operand: the generator writes it into the
        kernel, where a factor of one is not a multiplication and the zeros
        never become operations."""
        # the space-time predictor scales the star matrices by the timestep
        # outside the kernel, which a cell that does not carry them cannot do
        self.factoredStar = bool(kwargs.get("factored_star", False)) and kwargs.get(
            "solver"
        ) not in ("stp",)
        if not self.factoredStar:
            return

        mechanisms = getattr(self, "numMechanisms", 0)
        elastic = total_extent(self.primaryGroups())
        perMechanism = total_extent(self.mechanismGroups()) if mechanisms > 0 else 0

        count, entries, origins = coefficients.composed(
            self.name(), kwargs.get("solver"), mechanisms, elastic, perMechanism
        )
        self._solverCoefficientCount = count
        self._solverCoefficientOrigins = origins
        # a solver that keeps the mechanism index in a dimension of its own
        # carries a narrower star than the quantity count suggests, so take the
        # extents from the star itself
        starSpp = self.db.star[0].spp()
        shape, values = coefficients.structure_values(count, entries, starSpp.shape)

        self.starStructure = Tensor(
            "starStructure", shape, spp=values, addressing=AddressingMode.IMMEDIATE
        )
        self.materialCoefficients = Tensor("materialCoefficients", (count,))
        self.referenceGradients = [
            Tensor(f"referenceGradients({dim})", (3,)) for dim in range(3)
        ]
        self.starAssembled = [
            Tensor(f"starAssembled({dim})", starSpp.shape, spp=starSpp, temporary=True)
            for dim in range(3)
        ]

        self._configureNodalMaterial(kwargs, count, entries, starSpp)

    def _configureNodalMaterial(self, kwargs, count, entries, starSpp):
        """Sets up the tensors for a material that varies inside a cell.

        The coefficients then carry a point index, and the product of material
        and derivative has to be formed where the samples are and projected
        back. Which points those are the build decides; nothing here depends on
        the choice beyond the two matrices that read and project.

        The structure is split per coefficient. Written as one product of
        coefficients, nodal values and structure, the generator materializes an
        intermediate over (coefficient, point, quantity), which at order six
        and the conical-product set is larger than everything else in the
        kernel together.
        """
        self.nodalMaterial = self.factoredStar and bool(kwargs.get("material_nodal", False))
        if not self.nodalMaterial:
            return

        points = material.tensors(
            self._matricesDir, self, kwargs["material_points"]
        )
        self.materialEval = points["materialEval"]
        self.materialProject = points["materialProject"]
        npoints = self.materialEval.shape()[0]

        # one structure per coefficient, written into the kernel as before
        perCoefficient = [
            [e for e in entries if e.coefficient == a] for a in range(count)
        ]
        shape = (3,) + tuple(starSpp.shape)
        self.coefficientStructure = []
        for a in range(count):
            values = {}
            for entry in perCoefficient[a]:
                idx = (entry.dim, entry.row, entry.column)
                values[idx] = repr(float(values.get(idx, 0.0)) + entry.factor)
            self.coefficientStructure.append(
                Tensor(f"coefficientStructure({a})", shape, spp=values,
                       addressing=AddressingMode.IMMEDIATE)
            )

        # the Jacobian rows are a per-cell constant, so the fold happens once
        # and is reused by every step of the chain
        self.structureFolded = [
            [Tensor(f"structureFolded({dim},{a})", tuple(starSpp.shape), temporary=True)
             for a in range(count)]
            for dim in range(3)
        ]
        self.nodalCoefficients = [
            Tensor(f"nodalCoefficients({a})", (npoints,)) for a in range(count)
        ]
        quantities = starSpp.shape[0]
        self.nodalValues = Tensor("nodalValues", (npoints, quantities), temporary=True)
        self.nodalProduct = Tensor("nodalProduct", (npoints, starSpp.shape[1]),
                                   temporary=True)

    def solverCoefficientCount(self):
        """How many scalars the operator this solver applies is linear in."""
        return getattr(self, "_solverCoefficientCount", 0)

    def solverCoefficientOrigins(self):
        """Where each of those scalars comes from -- the material, or the run."""
        return getattr(self, "_solverCoefficientOrigins", [])

    def nodalAssembly(self):
        """Folds the Jacobian rows into the structure, once per kernel."""
        if not getattr(self, "nodalMaterial", False):
            return []
        return [
            self.structureFolded[dim][a]["qp"]
            <= self.referenceGradients[dim]["j"] * self.coefficientStructure[a]["jqp"]
            for dim in range(3)
            for a in range(len(self.coefficientStructure))
        ]

    def nodalApply(self, source, target, operators):
        """One application of the operator where the material varies inside the
        cell: read the derivative at the sample points, multiply by the
        material there, and project the result back.

        `operators` gives the modal operator per direction -- the stiffness for
        a derivative step, whatever the caller needs otherwise.
        """
        statements = []
        first = True
        for dim in range(3):
            statements.append(
                self.nodalValues["nq"]
                <= self.materialEval["nk"] * operators[dim][self.t("kl")] * source["lq"]
            )
            for a, coefficient in enumerate(self.nodalCoefficients):
                term = (
                    coefficient["n"]
                    * self.nodalValues["nq"]
                    * self.structureFolded[dim][a]["qp"]
                )
                statements.append(
                    self.nodalProduct["np"]
                    <= (term if first else self.nodalProduct["np"] + term)
                )
                first = False
        statements.append(target["kp"] <= self.materialProject["kn"] * self.nodalProduct["np"])
        return statements

    def starAssembly(self):
        """The statements that put the star matrices together, or none where a
        cell carries them assembled already."""
        if getattr(self, "nodalMaterial", False):
            return self.nodalAssembly()
        if not self.factoredStar:
            return []
        return [
            self.starAssembled[dim]["qp"]
            <= self.referenceGradients[dim]["j"]
            * self.materialCoefficients["a"]
            * self.starStructure["ajqp"]
            for dim in range(3)
        ]

    def configure(self, matricesDir, memLayout, kwargs, extra=()):
        """Reads this equation's matrix file, plus any the solver needs, and
        resolves the memory layout across all of them."""
        self._matricesDir = matricesDir
        clones = dict(self.StarClones)
        self.db.update(self.readMatrices(matricesDir, clones))
        for path in extra:
            self.db.update(parseJSONMatrixFile(path, clones))
        self.finishConfigure(memLayout, clones, kwargs)
        return clones

    def starMatrix(self, dim):
        """The star matrix a time-stepping kernel applies."""
        return self.starAssembled[dim] if self.factoredStar else self.db.star[dim]

    def starMatrixSetup(self, dim):
        """The star matrix an initialization kernel is handed.

        Always the assembled one: these run once on the host, where the cell's
        coefficients are at hand and putting the matrix together costs nothing
        worth generating a kernel for."""
        return self.db.star[dim]

    def stiffSourceRows(self):
        """Source rows a space-time predictor has to factorise separately, as
        (quantity, target, scalar name). Empty unless the source term is
        stiff."""
        return []

    def mapToVelocities(self):
        return self.extractVelocities().T

    def mapToTractions(self):
        return self.extractTractions().T

    def tractionMatrixSpp(self):
        """Sparsity pattern of the traction averaging matrices b+ and b-.

        Entry (q, p) is the weight with which quantity q of one side enters component p of the
        interface traction, so the pattern is the transpose of extractTractions, restricted to
        the three components the frictional work is computed with. Materials whose impedance is
        not diagonal reach more rows and override this.
        """
        return self.mapToTractions()[:, :3]

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

        fluxScale = Scalar("fluxScale")
        computeFluxSolverLocal = (
            self.AplusT["ij"]
            <= fluxScale
            * self.Tinv["ki"]
            * (self.QgodLocal["kq"] * self.starMatrixSetup(0)["ql"] + self.QcorrLocal["kl"])
            * self.T["jl"]
        )
        generator.add("computeFluxSolverLocal", computeFluxSolverLocal)

        computeFluxSolverNeighbor = (
            self.AminusT["ij"]
            <= fluxScale
            * self.Tinv["ki"]
            * (
                self.QgodNeighbor["kq"] * self.starMatrixSetup(0)["ql"]
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
