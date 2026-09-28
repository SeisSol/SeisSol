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
    extra_face_blocks,
    layout,
    role_offset,
    rotation_spp,
    total_extent,
    traction_selector,
    velocity_selector,
    voigt_weights,
    well_formed,
)
from yateto import Scalar, Tensor, simpleParameterSpace
from yateto.input import (
    memoryLayoutFromFile,
    parseJSONMatrixFile,
    parseXMLMatrixFile,
)
from yateto.memory import CSCMemoryLayout, PatternMemoryLayout
from yateto.type import AddressingMode
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

        # every equation reads further matrices from here, whether or not it
        # goes through configure() to do so
        self._matricesDir = matricesDir

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
            self.db.V2mTo2JacobiQuad["ik"] * self.db.MV2nTo2m[self.t("kj")],
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
            * self.db.V2nTo2m[self.t("kl")],
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
        self._configureRotationLayout(kwargs)
        self._configureStarAssembly(kwargs)

    def _configureRotationLayout(self, kwargs):
        """Stores the face rotation by its pattern rather than as a full square.

        The rotation is block diagonal -- one block per quantity group, and no
        group mixes with another -- so a dense square carries a majority of
        structural zeros. The pattern is already declared; only the layout was
        dense. Both matrices sit in the boundary face data and both feed the
        nodal boundary projections every timestep, so the empty blocks cost
        memory and operations there.

        The old GPU interface (gemmforge/chainforge) reads its operands as
        dense, so a GPU build served by it keeps the dense layout -- the same
        reservation the dynamic rupture rotation makes."""
        if kwargs.get("old_gpu_interface", True) and "gpu" in (
            kwargs.get("targets") or []
        ):
            return
        self.T.setMemoryLayout(CSCMemoryLayout)
        self.Tinv.setMemoryLayout(CSCMemoryLayout)

    def _configureStarAssembly(self, kwargs):
        """Sets up the tensors a cell carries where it holds the coefficients
        of its operator rather than the matrices they fold into.

        The structure the two fold into is a signed permutation, so it is
        stated as an immediate operand: the generator writes it into the
        kernel, where a factor of one is not a multiplication and the zeros
        never become operations."""
        self.factoredStar = bool(kwargs.get("factored_star", False))
        # set here as well, so that every solver can ask without knowing whether
        # the build got as far as the nodal configuration
        self.nodalMaterial = False
        self.nodalFaceFlux = False
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
        self.nodalMaterial = self.factoredStar and bool(
            kwargs.get("material_nodal", False)
        )
        if not self.nodalMaterial:
            return
        # The nodal chain contracts the derivative matrices over the modes
        # that carry a derivative at all. At the lowest orders that range
        # starts past the first stored row, where a CSC layout cannot be
        # sliced; a layout by pattern can. It is not aligned: the kernels
        # unroll the pattern, and the view the C++ side reads it through
        # strides by the tensor's own rows, not by an aligned extent.
        for derivative in self.db.kDivM.values():
            derivative.setMemoryLayout(PatternMemoryLayout, alignStride=False)
        # A face carries its operator as scalars only where that operator is
        # those scalars. Where it is not, the material still varies inside the
        # cell and the face keeps the one operator per side that is assembled
        # from the material of the two cells sharing it.
        self.nodalFaceFlux = self.fluxDecomposes()

        points = material.tensors(self._matricesDir, self, kwargs["material_points"])
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
                Tensor(
                    f"coefficientStructure({a})",
                    shape,
                    spp=values,
                    addressing=AddressingMode.IMMEDIATE,
                )
            )

        # the Jacobian rows are a per-cell constant, so the fold happens once
        # and is reused by every step of the chain
        self.structureFolded = [
            [
                Tensor(
                    f"structureFolded({dim},{a})", tuple(starSpp.shape), temporary=True
                )
                for a in range(count)
            ]
            for dim in range(3)
        ]
        self.nodalCoefficients = [
            Tensor(f"nodalCoefficients({a})", (npoints,)) for a in range(count)
        ]
        quantities = starSpp.shape[0]
        self.nodalOperatorAssembled = (
            kwargs.get("material_operator", coefficients.OPERATOR_FORMS[0])
            == "assembled"
        )
        self.starAtPoint = [
            Tensor(
                f"starAtPoint({dim})",
                (npoints,) + tuple(starSpp.shape),
                temporary=True,
            )
            for dim in range(3)
        ]

        # kept as well, since a solver whose field carries more than modes and
        # quantities needs the same two shapes one index wider
        self.nodalValuesShape = (npoints, quantities)
        self.nodalProductShape = (npoints, starSpp.shape[1])
        self.nodalValues = self.nodalTemporary("nodalValues", self.nodalValuesShape)
        self.nodalProduct = self.nodalTemporary("nodalProduct", self.nodalProductShape)

        if self.nodalFaceFlux:
            self._configureNodalFlux()
        self._configureNodalSource(kwargs)

    def _configureNodalSource(self, kwargs):
        """The tensors a source term needs where the material varies inside a
        cell.

        The same idea as the flux: the source is a handful of scalars the
        material supplies at fixed entries, so a cell carries those at the
        sample points and the kernel puts the term together where they are.
        A solver without a source term declares none and gets none.
        """
        self._sourceCoefficientCount = 0
        prototype = self.sourceStructurePrototype()
        if prototype is None:
            return

        shape = tuple(prototype.shape())
        count, values, origins = coefficients.source_composed(
            self.name(),
            kwargs.get("solver"),
            getattr(self, "numMechanisms", 0),
            shape,
            total_extent(self.primaryGroups()),
            total_extent(self.mechanismGroups()),
        )
        if count == 0:
            return

        self._sourceCoefficientCount = count
        self._sourceCoefficientOrigins = origins
        self.sourceStructure = [
            Tensor(
                f"sourceStructure({a})",
                shape,
                spp={
                    key[1:]: repr(float(factor))
                    for key, factor in values.items()
                    if key[0] == a
                },
                addressing=AddressingMode.IMMEDIATE,
            )
            for a in range(count)
        ]
        npoints = self.materialEval.shape()[0]
        self.sourceCoefficients = [
            Tensor(f"sourceCoefficients({a})", (npoints,)) for a in range(count)
        ]
        # the field at the sample points, in whatever indices the source term
        # sums over -- the quantities, and the mechanisms where there are any --
        # and what the source makes of it. The second one is the source's own:
        # a solver may write fewer quantities here than its operator does, and
        # what it does not write has to stay out of the projection.
        self.nodalSourceValues = self.nodalTemporary(
            "nodalSourceValues", (npoints,) + shape[:-1]
        )
        self.nodalSourceProduct = self.nodalTemporary(
            "nodalSourceProduct", (npoints, shape[-1])
        )
        self.sourceDeviation = [
            Tensor(f"sourceDeviation({a})", (npoints,))
            for a in range(self.sourceDeviationCount())
        ]

    def sourceStructurePrototype(self):
        """The tensor this solver states its source term in, or none where it
        has no source term to state."""
        matrix = getattr(self, "sourceMatrix", None)
        return matrix() if matrix is not None else None

    def sourceCoefficientCount(self):
        """How many scalars the source term of this solver is linear in."""
        return getattr(self, "_sourceCoefficientCount", 0)

    def sourceCoefficientOrigins(self):
        """Where each of those scalars comes from -- the material, or the run."""
        return getattr(self, "_sourceCoefficientOrigins", [])

    def sourceDeviationCount(self):
        """How many scalars a cell carries as the difference between its source
        term at a sample point and the one it carries for itself.

        Only a solver that puts the source term inside a solve needs them: it
        factorises that solve once for the cell, and what a sample point
        deviates from it has to be carried separately. Zero for a solver that
        applies the source as a product, which reads the samples directly.
        """
        return 0

    def _configureNodalFlux(self):
        """The tensors a face carries where the material varies along it.

        In face coordinates the flux operator is a few scalars times fixed
        entries -- ten for the Godunov flux of an elastic medium, one more for
        the Rusanov penalty, and a set per relaxation mechanism; measured, and
        stated in the generated tables -- so a face
        holds those scalars per node instead of a matrix. In global coordinates
        it is not: rotated, the same operator occupies every entry and spans
        far more than ten dimensions, so the rotation belongs in the kernel and
        not in what a face stores.

        Only the forward rotation is stored. Its inverse follows from it by the
        Voigt weights, and neither weight costs a multiplication: one half is
        folded into the structure the coefficients scale, the other is a
        constant diagonal the field passes through on its way to the face.
        """
        faceNodes = material.addNeighborFaceMatrices(self, self._matricesDir)
        quantities = self.numQuantities()
        extended = self.numExtendedQuantities()
        weights = voigt_weights(self.quantityBlocks())

        # the Rusanov diagonal spans the square part of the Godunov state
        names, self.fluxSources, entries = coefficients.flux_decomposition(
            self.extendedBlocks(), self.QgodLocal.shape()[0]
        )
        count = len(names)
        self._fluxCoefficientCount = count
        self._fluxEntries = entries

        self.fluxStructure = [
            Tensor(
                f"fluxStructure({a})",
                (extended, extended),
                spp={
                    (e.row, e.column): repr(float(e.factor) * weights[e.row])
                    for e in entries
                    if e.coefficient == a
                },
                addressing=AddressingMode.IMMEDIATE,
            )
            for a in range(count)
        ]
        # the field arrives with the quantities the cell carries and leaves with
        # the ones the operator writes, so the weights inject as well as scale
        self.inverseVoigtWeights = Tensor(
            "inverseVoigtWeights",
            (quantities, extended),
            spp={(q, q): repr(1.0 / weights[q]) for q in range(quantities)},
            addressing=AddressingMode.IMMEDIATE,
        )
        self.fluxCoefficientsLocal = [
            Tensor(f"fluxCoefficientsLocal({a})", (faceNodes,)) for a in range(count)
        ]
        self.fluxCoefficientsNeighbor = [
            Tensor(f"fluxCoefficientsNeighbor({a})", (faceNodes,)) for a in range(count)
        ]
        # A kernel that applies the flux of all four faces at once needs the
        # operands of each face apart: the scalars and the rotation both belong
        # to one face. The rotations alias what a face stores, so they take the
        # layout of the one rotation the per-face kernels read.
        self.fluxCoefficientsLocalAll = [
            [
                Tensor(f"fluxCoefficientsLocalAll({face},{a})", (faceNodes,))
                for a in range(count)
            ]
            for face in range(4)
        ]
        rotationLayout = self.T.memoryLayout()
        self.TAll = []
        for face in range(4):
            rotation = Tensor(f"TAll({face})", self.T.shape(), spp=self.T.spp())
            rotation.setMemoryLayout(
                rotationLayout.__class__,
                alignStride=rotationLayout.alignedStride(),
                alignmentArch=rotationLayout.alignmentArch(),
            )
            self.TAll.append(rotation)
        shape = (faceNodes, extended)
        self.faceValues = self.nodalTemporary("faceValues", shape)
        self.faceRotated = self.nodalTemporary("faceRotated", shape)
        self.faceProduct = self.nodalTemporary("faceProduct", shape)
        self.faceBack = self.nodalTemporary("faceBack", shape)

    def nodalFlux(
        self, source, target, toFace, lift, coefficientsOfFace, rotation=None
    ):
        """One face contribution where the operator varies along the face.

        The field is read at the nodes of the face and turned into the face
        coordinates the scalars are stated in, the operator is applied
        there, and the result is turned back and lifted into the cell with the
        operator the nodal boundary conditions already use. The rotation is the
        same matrix both ways, once transposed against the quantity the field
        carries and once against the quantity the result is written in.

        `toFace` reads the field at the nodes of the face, with the node index
        first and the mode index second; `lift` goes the other way. Both come
        indexed, because how a matrix is laid out is the caller's to state.
        `rotation` is the face rotation the kernel reads, T unless a kernel
        applies more than one face and needs one per face.
        """
        rotation = self.T if rotation is None else rotation
        statements = [
            self.faceValues["nq"]
            <= toFace * source["lk"] * self.inverseVoigtWeights["kq"],
            self.faceRotated["nk"] <= self.faceValues["nq"] * rotation["qk"],
        ]
        first = True
        for a, coefficient in enumerate(coefficientsOfFace):
            term = (
                coefficient["n"] * self.faceRotated["nk"] * self.fluxStructure[a]["kl"]
            )
            statements.append(
                self.faceProduct["nl"]
                <= (term if first else self.faceProduct["nl"] + term)
            )
            first = False
        statements.append(
            self.faceBack["np"] <= self.faceProduct["nl"] * rotation["pl"]
        )
        statements.append(target["kp"] <= target["kp"] + lift * self.faceBack["np"])
        return statements

    def nodalLocalFluxAll(self, source, target):
        """The local flux of all four faces in one kernel, where the operator
        varies along a face.

        A face whose flux the cell does not take -- a fault face -- carries
        zero scalars, the same way the matrix form carries a zero matrix, so
        the four faces need no case distinction here.
        """
        statements = []
        for face in range(4):
            statements += self.nodalFlux(
                source,
                target,
                self.db.V3mTo2nFace[face][self.t("nl")],
                self.db.project2nFaceTo3m[face]["kn"],
                self.fluxCoefficientsLocalAll[face],
                rotation=self.TAll[face],
            )
        return statements

    def solverCoefficientCount(self):
        """How many scalars the operator this solver applies is linear in."""
        return getattr(self, "_solverCoefficientCount", 0)

    def solverCoefficientOrigins(self):
        """Where each of those scalars comes from -- the material, or the run."""
        return getattr(self, "_solverCoefficientOrigins", [])

    def nodalTemporary(self, name, shape):
        """A temporary of the nodal path.

        It carries the same field a kernel's operands do, so a build that fuses
        simulations gives it that index as well; everything else about the
        nodal path is per cell and shared across them."""
        return OptionalDimTensor(
            name,
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            shape,
            temporary=True,
        )

    def nodalAssembly(self):
        """Folds the Jacobian rows into the structure, once per kernel.

        Where the build asks for the assembled form, the coefficients go in as
        well and what comes out is one operator per sample point. That trades
        the products a kernel does at every application for a temporary over
        the points, so which one is cheaper depends on how often the kernel
        applies the operator and on the machine.
        """
        if not self.nodalMaterial:
            return []
        statements = [
            self.structureFolded[dim][a]["qp"]
            <= self.referenceGradients[dim]["j"] * self.coefficientStructure[a]["jqp"]
            for dim in range(3)
            for a in range(len(self.coefficientStructure))
        ]
        if self.nodalOperatorAssembled:
            for dim in range(3):
                folded = None
                for a, coefficient in enumerate(self.nodalCoefficients):
                    term = coefficient["n"] * self.structureFolded[dim][a]["qp"]
                    folded = term if folded is None else folded + term
                statements.append(self.starAtPoint[dim]["nqp"] <= folded)
        return statements

    def nodalApply(
        self,
        source,
        target,
        operators,
        spectator="",
        temporaries=None,
        accumulate=False,
        scalar=None,
    ):
        """One application of the operator where the material varies inside the
        cell: read the derivative at the sample points, multiply by the
        material there, and project the result back.

        `operators` gives the modal operator per direction -- the stiffness for
        a derivative step, whatever the caller needs otherwise. `spectator`
        names indices the operator leaves alone, for a field that carries more
        than modes and quantities; the caller then hands over the two
        temporaries those indices widen. `scalar` scales the result, and
        `accumulate` adds it to what the target holds instead of replacing it.
        """
        values, product = (
            temporaries
            if temporaries is not None
            else (self.nodalValues, self.nodalProduct)
        )
        statements = []
        first = True
        for dim in range(3):
            statements.append(
                values["nq" + spectator]
                <= self.materialEval["nk"]
                * operators[dim][self.t("kl")]
                * source["lq" + spectator]
            )
            if self.nodalOperatorAssembled:
                terms = [values["nq" + spectator] * self.starAtPoint[dim]["nqp"]]
            else:
                terms = [
                    coefficient["n"]
                    * values["nq" + spectator]
                    * self.structureFolded[dim][a]["qp"]
                    for a, coefficient in enumerate(self.nodalCoefficients)
                ]
            for term in terms:
                statements.append(
                    product["np" + spectator]
                    <= (term if first else product["np" + spectator] + term)
                )
                first = False
        projected = self.materialProject["kn"] * product["np" + spectator]
        if scalar is not None:
            projected = scalar * projected
        statements.append(
            target["kp" + spectator]
            <= (target["kp" + spectator] + projected if accumulate else projected)
        )
        return statements

    def sourceTerm(self, source, target):
        """The source term added to a target, in whichever shape this build
        forms it: from the matrix a cell carries, or from the scalars it
        carries at the sample points. Nothing at all where the solver has no
        source term."""
        if self.sourceMatrix() is None:
            return []
        if self.sourceCoefficientCount() == 0:
            return [
                target["kp"] <= target["kp"] + source["kq"] * self.sourceMatrix()["qp"]
            ]
        return self.nodalSource(source, target, "nq")

    def nodalSource(
        self,
        source,
        target,
        contract,
        spectator="",
        temporaries=None,
        coefficients=None,
        scalar=None,
    ):
        """The source term where the material varies inside the cell.

        No derivative is taken here, so the field goes straight to the sample
        points, is multiplied by the source the material has there, and comes
        back. `contract` names the indices the source term sums over -- the
        quantities, and the mechanisms where a solver keeps them in a dimension
        of their own. `coefficients` is what the material says at those points,
        or what it says beyond what the cell carries.
        """
        coefficients = self.sourceCoefficients if coefficients is None else coefficients
        values, product = (
            temporaries
            if temporaries is not None
            else (self.nodalSourceValues, self.nodalSourceProduct)
        )
        statements = [
            values[contract + spectator]
            <= self.materialEval["nk"] * source["k" + contract[1:] + spectator]
        ]
        first = True
        for a, coefficient in enumerate(coefficients):
            term = (
                coefficient["n"]
                * values[contract + spectator]
                * self.sourceStructure[a][contract[1:] + "p"]
            )
            statements.append(
                product["np" + spectator]
                <= (term if first else product["np" + spectator] + term)
            )
            first = False
        projected = self.materialProject["kn"] * product["np" + spectator]
        if scalar is not None:
            projected = scalar * projected
        statements.append(
            target["kp" + spectator] <= target["kp" + spectator] + projected
        )
        return statements

    def starAssembly(self):
        """The statements that put the star matrices together, or none where a
        cell carries them assembled already."""
        if self.nodalMaterial:
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

    def fluxDecomposes(self):
        """Whether the flux operator of a face is the handful of scalars
        :func:`kernels.coefficients.flux_decomposition` states it as.

        A layout whose face-local vectors carry rows beyond the traction and
        the velocity of one medium couples across a face in ways those scalars
        do not name: a second, fluid, medium reaches the whole operator, which
        then occupies forty-three of its one hundred and sixty-nine entries
        instead of thirteen, and the scalars leave a third of it behind.
        """
        return not extra_face_blocks(self.extendedBlocks())

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
            * (
                self.QgodLocal["kq"] * self.starMatrixSetup(0)["ql"]
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
        # the reparametrisation of a shared face. The neighbour flux reads it
        # folded into fPrT and the nodal flux as a renumbering of the face
        # nodes, so no kernel names it; it is what the renumbering is checked
        # against, in every build, so it has to reach the generated code.
        for orientation in self.db.fP.values():
            include_tensors.add(orientation)
        if self.nodalFaceFlux:
            include_tensors.add(self.db.M2)
            # the nodal flux is checked against the matrix form, which is
            # built from these. No flux kernel of such a build names them, and
            # a device build premultiplies them even where the matrix form
            # remains, so they reach the generated code only from here.
            for family in (self.db.rDivM, self.db.fMrT):
                for member in family.values():
                    include_tensors.add(member)
