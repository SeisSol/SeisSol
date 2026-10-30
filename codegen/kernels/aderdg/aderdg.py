# SPDX-FileCopyrightText: 2019 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff

from abc import ABC, abstractmethod
from contextlib import contextmanager

import numpy as np
from kernels import coefficients, material
from kernels.basis import Basis, count, dense, shared
from kernels.common import generate_kernel_name_prefix
from kernels.multsim import OptionalDimTensor
from kernels.quantities import (
    FaceRole,
    extra_face_blocks,
    face_block,
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
from yateto.ast.indices import Indices
from yateto.ast.transformer import DeduceIndices, EquivalentSparsityPattern
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


def cloneOperand(tensor, name):
    """A tensor of its own with the shape, the pattern and the memory layout of
    `tensor`, for a kernel that reads several of it at once: the copies are
    laid out as the one a cell stores, so each can be bound to where the cell
    stores one of them."""
    clone = Tensor(name, tensor.shape(), spp=tensor.spp().as_ndarray())
    layout = tensor.memoryLayout()
    clone.setMemoryLayout(type(layout), alignStride=layout.alignedStride())
    if clone.memoryLayout().requiredReals() != layout.requiredReals():
        raise ValueError(f"{name} is not laid out as {tensor.name()}")
    return clone


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


def embedFamily(db, baseName, size, alignStride):
    """Replaces a tensor family in `db` with its members placed in the leading
    block of a size x size matrix: a matrix between the bases of two degrees
    then applies to the leading rows of a tensor of the cell, which hold a
    field of that degree. The entries keep their indices, so this holds for a
    matrix stored transposed as well."""
    db[baseName] = {
        group: Tensor(
            tensor.name(),
            (size, size),
            dict(tensor.values().items()),
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
        self.basis = Basis(matricesDir)

        self.db = parseXMLMatrixFile(
            f"{matricesDir}/aderdg-{order}.xml",
            transpose=self.transpose,
            alignStride=self.alignStride,
        )
        if not (self.basis.hierarchical and self.basis.orthogonal):
            # the derivatives through the basis of one degree less, and the bases
            # of the lower degrees (modalDerivative,
            # material.projectionOntoDegree), read like the derivatives
            self.db.update(
                parseJSONMatrixFile(
                    f"{matricesDir}/lower-{order}.json",
                    transpose=self.transpose,
                    alignStride=self.alignStride,
                )
            )
        # The derivative kernels contract against -kDivMT throughout. Folding
        # the sign into the matrix here rather than scaling the operand keeps
        # the global matrices read-only, which is what lets them be shared and
        # placed wherever a target wants them. It has to happen before the
        # memory layout configuration is applied, since that mutates whatever
        # tensors are in the database at the time.
        negateFamily(self.db, "kDivMT", self.alignStride)
        if "kDivMTDegree" in self.db:
            negateFamily(self.db, "kDivMTDegree", self.alignStride)
            for family in ("kDivMTDegree", "elevateDegree"):
                embedFamily(
                    self.db, family, self.num3DBasisFunctions(), self.alignStride
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
        self.db.update(
            parseJSONMatrixFile(shared(matricesDir, "sampling_directions.json"))
        )
        self.db.update(parseJSONMatrixFile(f"{matricesDir}/mass-{order}.json"))

        # A diagonal mass matrix is sparse, a dense one is not; a memory layout
        # file that names them decides over this.
        for mass in (self.db.M2, self.db.M3):
            pattern = mass.spp().as_ndarray()
            if not pattern[~np.eye(pattern.shape[0], dtype=bool)].any():
                mass.setMemoryLayout(CSCMemoryLayout)

        qShape = (self.num3DBasisFunctions(), self.numQuantities())
        self.Q = OptionalDimTensor(
            "Q", "s", multipleSimulations, 0, qShape, alignStride=True
        )

        self.I = OptionalDimTensor(
            "I", "s", multipleSimulations, 0, qShape, alignStride=True
        )

        Aplusminus_spp = self.flux_solver_spp()
        self.AplusT = Tensor("AplusT", Aplusminus_spp.shape, spp=Aplusminus_spp)
        self.AminusT = Tensor("AminusT", Aplusminus_spp.shape, spp=Aplusminus_spp)
        trans_spp = self.transformation_spp()
        self.T = Tensor("T", trans_spp.shape, spp=trans_spp)
        trans_inv_spp = self.transformation_inv_spp()
        self.Tinv = Tensor("Tinv", trans_inv_spp.shape, spp=trans_inv_spp)
        # The points of a face a flux that varies along it is integrated at:
        # those of the symmetric rule exact for the product of two face basis
        # functions, with what reads a field there and integrates against it.
        # The points themselves are read by the host, which queries the material
        # there, so they keep the orientation of the file in every build.
        self.db.update(
            parseJSONMatrixFile(
                f"{matricesDir}/face-sq-matrices-{order}.json",
                {},
                alignStride=self.alignStride,
                transpose=lambda name: self.transpose(name)
                and name not in ("faceQuadPoints", "faceQuadWeights"),
            )
        )
        # The rotation once per point of a face, for a face whose normal turns
        # along it. It is declared in every build, used or not, so that the host
        # can size what a face stores without asking which build it is.
        rotationNodesSpp = np.broadcast_to(
            trans_spp, (self.numFaceFluxPoints(),) + tuple(trans_spp.shape)
        ).copy()
        self.TNodes = Tensor("TNodes", rotationNodesSpp.shape, spp=rotationNodesSpp)
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

        # the lift of what a flux gives at the points of a face, integrated
        # against the face basis with the weights of the rule
        projectFaceQuadTo3m = tensor_collection_from_constant_expression(
            base_name="projectFaceQuadTo3m",
            expressions=lambda i: self.db.rDivM[i][self.t("jk")]
            * self.db.VFaceQuadTo2m[self.t("kl")],
            group_indices=simpleParameterSpace(4),
            target_indices="jl",
        )
        self.db.update(projectFaceQuadTo3m)

        selectVelocitySpp = self.mapToVelocities()[:, :3]
        self.selectVelocity = Tensor(
            "selectVelocity",
            selectVelocitySpp.shape,
            selectVelocitySpp,
            CSCMemoryLayout,
        )

        # The traction weights are computed per fault face from the impedances
        # (DynamicRuptureMatrices), or per point of it where the material varies
        # inside a cell, so only their pattern is known here. Passed
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

    def numFaceFluxPoints(self):
        """How many points of a face a flux that varies along it is integrated
        at."""
        return self.db.faceQuadPoints.shape()[0]

    def faceLift(self, face):
        """Lifts what a flux gives at the points of a face into the cell: the
        integral against the face basis, by the rule the points are of, and the
        lift of the face basis into the cell."""
        return self.db.projectFaceQuadTo3m[face]["kn"]

    def num3DBasisFunctions(self):
        return self.order * (self.order + 1) * (self.order + 2) // 6

    def num3DQuadraturePoints(self):
        return (self.order + 1) ** 3

    def godunov_spp(self):
        shape = (self.numQuantities(), self.numQuantities())
        return np.ones(shape, dtype=bool)

    def rusanovDiagonal(self):
        """The quantities the Rusanov form puts its penalty on: the diagonal of
        the Godunov state. A solver that folds the relaxation into its
        quantities states the Godunov state over the elastic rows alone."""
        spp = self.godunov_spp()
        return [q for q in range(min(spp.shape)) if spp[q, q]]

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
        return parseJSONMatrixFile(
            shared(matricesDir, f"equation-{self.name()}.json"), clones
        )

    def finishConfigure(self, memLayout, clones, kwargs):
        """Resolves the memory layout and stores the generator arguments. Split
        out so that an equation can reshape its matrices in between."""
        memoryLayoutFromFile(memLayout, self.db, clones)
        self.kwargs = kwargs
        self._configureRotationLayout(kwargs)
        self._configureStarAssembly(kwargs)
        # after the memory layout: the trace matrices take the layout that
        # traceMatrixArgs gives them, whatever a layout file says about the
        # matrices they are formed from
        self._configureTraces()

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
        # one block-diagonal rotation per point, packed by the pattern they share
        self.TNodes.setMemoryLayout(PatternMemoryLayout, alignStride=False)

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
        self.narrowDerivatives = False
        # the points the operator is formed at: the barycentre, where it is one
        # per cell; the nodal configuration replaces it with its point set
        self.operatorNodes = Tensor(
            "operatorNodes",
            (1, 3),
            spp={(0, d): repr(0.25) for d in range(3)},
        )
        # the rows of the inverse Jacobian a cell carries apart from its
        # material; a curved cell has them at every point the operator is
        # formed at, see _configureCurvilinear
        self.referenceGradients = [
            Tensor(f"referenceGradients({dim})", (3,)) for dim in range(3)
        ]
        self.curvilinear = bool(kwargs.get("curvilinear", False))
        if self.curvilinear and not (
            self.factoredStar and kwargs.get("material_nodal", False)
        ):
            raise ValueError(
                "a curved cell carries its metric at the points the operator is "
                "formed at, which only a build with MATERIAL_NODAL has"
            )
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
        self._solverCoefficientEntries = entries
        # a solver that keeps the mechanism index in a dimension of its own
        # carries a narrower star than the quantity count suggests, so take the
        # extents from the star itself
        starSpp = self.db.star[0].spp()
        shape, values = coefficients.structure_values(count, entries, starSpp.shape)

        self.starStructure = Tensor(
            "starStructure", shape, spp=values, addressing=AddressingMode.IMMEDIATE
        )
        self.materialCoefficients = Tensor("materialCoefficients", (count,))
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
        # whether the build lets the time derivatives of the predictor keep
        # fewer modes with every derivative; the solver decides whether its
        # chain can (narrowsDerivatives)
        self.narrowDerivatives = bool(kwargs.get("material_narrow_derivatives", False))
        # The nodal chain contracts the derivative matrices over the modes
        # that carry a derivative at all. At the lowest orders that range
        # starts past the first stored row, or ends before the rows an aligned
        # CSC layout pads to, and a CSC layout can be sliced at neither; a
        # layout by pattern can. It is not aligned: the kernels unroll the
        # pattern, and the view the C++ side reads it through strides by the
        # tensor's own rows, not by an aligned extent. Both the derivative
        # matrices and their transposes go that way, since the chain reads
        # either depending on the solver.
        for derivative in list(self.db.kDivM.values()) + list(self.db.kDivMT.values()):
            derivative.setMemoryLayout(PatternMemoryLayout, alignStride=False)
        # A face carries its operator as scalars only where that operator is
        # those scalars. Where it is not, the material still varies inside the
        # cell and the face keeps the one operator per side that is assembled
        # from the material of the two cells sharing it.
        self.nodalFaceFlux = self.fluxDecomposes()

        self.materialPoints = kwargs["material_points"]
        points = material.tensors(self._matricesDir, self, self.materialPoints)
        self.materialEval = points["materialEval"]
        self.materialProject = points["materialProject"]
        # The samples are the points of a quadrature rule, so the operator is
        # formed at them and projected back with the rule's weights: reading a
        # field there is materialEval, projecting what is formed there back to
        # the modes is materialProject.
        self.operatorEval = self.materialEval
        self.operatorProject = self.materialProject
        samples = self.materialEval.shape()[0]
        npoints = samples
        # the points the host evaluates the metric of a curved cell at
        self.operatorNodes = Tensor(
            "operatorNodes",
            points["materialNodes"].shape(),
            spp=dict(points["materialNodes"].values()),
            alignStride=self.alignStride("operatorNodes"),
        )
        if self.curvilinear:
            self._configureCurvilinear(kwargs, npoints)

        self.nodalCoefficients = [
            Tensor(f"nodalCoefficients({a})", (samples,)) for a in range(count)
        ]
        quantities = starSpp.shape[0]
        self.nodalOperatorAssembled = (
            kwargs.get("material_operator", coefficients.OPERATOR_FORMS[0])
            == "assembled"
        )

        # kept as well, since a solver whose field carries more than modes and
        # quantities needs the same two shapes one index wider
        self.nodalValuesShape = (npoints, quantities)
        self.nodalProductShape = (npoints, starSpp.shape[1])
        self.nodalValues = self.nodalTemporary("nodalValues", self.nodalValuesShape)
        self.nodalProduct = self.nodalTemporary("nodalProduct", self.nodalProductShape)

        # Where the metric is one per cell, the operator splits into a constant
        # part and the material (_configureLeftFactor). In a curved cell the
        # metric varies from point to point as the material does, and the two
        # form the operator at every point together (_configurePointOperator).
        self.leftFactor = not self.curvilinear
        if self.leftFactor:
            self._configureLeftFactor(entries, count, starSpp, npoints)
        else:
            self._configurePointOperator(entries, count, starSpp, npoints)
        if self.curvilinear:
            self._configureCurvedSplit(entries, count, starSpp, npoints)

        if self.nodalFaceFlux:
            self._configureNodalFlux()
        self._configureWeightAdjust()
        self._configureNodalSource(kwargs)

    def _configurePointOperator(self, entries, count, starSpp, npoints):
        """Sets up the operator formed at every sample point of a curved cell:
        one structure per coefficient, immediates in the kernel, which the
        metric and the coefficients at a point fold into the operator of each
        direction there."""
        shape = (3,) + tuple(starSpp.shape)
        operatorPattern = np.zeros(tuple(starSpp.shape), dtype=bool)
        self.coefficientStructure = []
        for a in range(count):
            values = {}
            for entry in entries:
                if entry.coefficient != a:
                    continue
                idx = (entry.dim, entry.row, entry.column)
                values[idx] = repr(float(values.get(idx, 0.0)) + entry.factor)
                operatorPattern[entry.row, entry.column] = True
            self.coefficientStructure.append(
                Tensor(
                    f"coefficientStructure({a})",
                    shape,
                    spp=values,
                    addressing=AddressingMode.IMMEDIATE,
                )
            )
        # The operator at a point has the entries some coefficient reaches in
        # some direction. Declared without them, it is a dense block, and every
        # product with it runs over all of its entries, most of them zero.
        pointPattern = np.broadcast_to(
            operatorPattern, (npoints,) + tuple(starSpp.shape)
        ).copy()
        self.starAtPoint = [
            Tensor(
                f"starAtPoint({dim})",
                (npoints,) + tuple(starSpp.shape),
                spp=pointPattern,
                temporary=True,
            )
            for dim in range(3)
        ]

    def _configureLeftFactor(self, entries, count, starSpp, npoints):
        """Sets up the operator of a straight-sided cell as a constant operator
        followed by the material.

        The kernels apply the operator of direction j as the field times its
        transpose, and a decomposition gives that transpose entry by entry: the
        entry of input quantity q and output column p is a sum of coefficients
        times factors. Two input rows (j, q) that enter every column with the
        same coefficients and factors are the same function of the material,
        so the operator is U_j L(c): U_j sends each row to its group, ones and
        zeros that do not depend on the material, and L(c) holds what each
        group gives in each column, linear in the coefficients.

        Where the metric is one per cell it goes into U alone. The derivatives
        of the three directions then combine into the groups in modal space,
        the field is read at the sample points once instead of once per
        direction, and the material is applied to as many columns as there are
        groups.
        """
        usage = {}
        for entry in entries:
            row = usage.setdefault((entry.dim, entry.row), {})
            key = (entry.coefficient, entry.column)
            row[key] = row.get(key, 0.0) + entry.factor
        groups = {}
        for row, use in usage.items():
            signature = tuple(
                sorted((key, factor) for key, factor in use.items() if factor)
            )
            if signature:
                groups.setdefault(signature, []).append(row)
        signatures = list(groups)
        self.leftGroups = len(signatures)

        quantities, columns = starSpp.shape
        unit = {}
        for group, signature in enumerate(signatures):
            for dim, row in groups[signature]:
                unit[(dim, row, group)] = "1.0"
        self.leftUnit = Tensor(
            "leftUnit",
            (3, quantities, self.leftGroups),
            spp=unit,
            addressing=AddressingMode.IMMEDIATE,
        )
        unitPattern = np.zeros((quantities, self.leftGroups), dtype=bool)
        for _, row, group in unit:
            unitPattern[row, group] = True
        # the unit of each direction with the metric of the cell folded in,
        # once per kernel
        self.leftUnitFolded = [
            Tensor(
                f"leftUnitFolded({dim})",
                unitPattern.shape,
                spp=unitPattern,
                temporary=True,
            )
            for dim in range(3)
        ]

        groupPattern = np.zeros((self.leftGroups, columns), dtype=bool)
        self.leftStructure = []
        for a in range(count):
            values = {}
            for group, signature in enumerate(signatures):
                for (coefficient, column), factor in signature:
                    if coefficient == a:
                        values[(group, column)] = repr(factor)
                        groupPattern[group, column] = True
            self.leftStructure.append(
                Tensor(
                    f"leftStructure({a})",
                    (self.leftGroups, columns),
                    spp=values,
                    addressing=AddressingMode.IMMEDIATE,
                )
            )
        # what the groups give at every point, where the build asks for the
        # assembled form
        self.leftAtPoint = Tensor(
            "leftAtPoint",
            (npoints,) + groupPattern.shape,
            spp=np.broadcast_to(groupPattern, (npoints,) + groupPattern.shape).copy(),
            temporary=True,
        )
        # the temporaries of the groups, one per shape and pattern
        self._leftModalTemporaries = {}
        self._leftValuesTemporaries = {}

    def _configureCurvilinear(self, kwargs, npoints):
        """Sets up what a curved cell carries beyond a straight-sided one.

        The rows of the inverse Jacobian vary inside a curved cell, so a cell
        carries them at the points the operator is formed at, where they go
        into the operator together with the material; the predictor applies
        the operator there in the strong form, as on a straight-sided cell.

        The corrector of a curved cell is weight-adjusted with the weight J W,
        J the Jacobian determinant: the cell lifts its face terms without
        material and with the surface Jacobian alone into its flux buffer, and
        the inverse weight applied there carries 1 / J at every point
        (weightScale, see weightAdjust). Its volume term goes through the
        weight as well, in the split form (_configureCurvedSplit).

        A face carries its rotation per point, since its normal turns along it;
        that makes the flux one that has to be applied at the points of the
        face, which is where only a flux decomposed into scalars per point is.
        """
        if (
            kwargs.get("material_operator", coefficients.OPERATOR_FORMS[0])
            != "assembled"
        ):
            raise ValueError(
                "a curved cell folds its metric into the operator at every point, "
                "which is the assembled form; MATERIAL_OPERATOR=factored would "
                "fold it at every application"
            )
        if not self.fluxDecomposes():
            raise ValueError(
                "a curved face needs its flux at its points, and this layout's "
                "flux is not the scalars per point that form takes"
            )
        self.referenceGradients = [
            Tensor(f"referenceGradients({dim})", (npoints, 3)) for dim in range(3)
        ]
        # 1 / J at the points the operator is formed at, which the inverse
        # weight of a curved cell carries on top of the one of its medium
        self.weightScale = Tensor("weightScale", (npoints,))

    def momentumColumns(self):
        """The columns of the operator that give the rate of the momentum: the
        ones of the velocity, which the divergence of the stress drives."""
        velocity = face_block(self.quantityBlocks(), FaceRole.VELOCITY)
        return set(range(velocity.offset, velocity.offset + velocity.extent))

    def _configureCurvedSplit(self, entries, count, starSpp, npoints):
        """Sets up the split form of the corrector of a curved cell.

        The scheme tests with the basis and weighs with J W (see
        _configureCurvilinear). In the strong form its energy balance has a
        remainder that no rule of the volume and of the faces integrates away,
        the error of the rules on the divergence theorem of the momentum rows
        on a curved cell: the strong form grows with the rules of some orders.
        The split form takes the rows of the momentum in the weak form and the
        others in the strong form. In the energy the two volume terms then
        cancel at every point the rule has, whatever the rule, and the faces
        leave the dissipation of the upwind flux.

        So a curved cell forms its volume term itself, into its flux buffer and
        through the weight, as moments without material (splitVolume): with J
        times the operator without material at every point,

        - the columns other than the momentum ones in the strong form, the
          derivative of the time integral at the points against the operator,
          projected back;
        - the momentum columns in the weak form, the time integral at the
          points against the operator, projected back, against the derivatives
          of the test functions (kDivM).

        The faces lift the strong form of the flux as on any face. For the
        momentum rows the weak form wants the flux alone, so a curved cell adds
        the normal flux of its own trace for those rows back (ownTraction), at
        the points of each face and with the scale of the face there, except on
        a fault, whose lift leaves the traction of the cell's own trace in
        place (see FaultSubtractsOwnTrace on the host).

        A curved cell holds J times the coefficients of the operator without
        material at its points (splitCoefficients) and those times the scale of
        each point of each face (ownTractionCoefficients). A constant state
        gives nothing, and on a straight-sided cell the split form is the
        strong form: the rules integrate both exactly there.
        """
        if self.sourceStructurePrototype() is not None:
            raise ValueError(
                "a curved cell takes the momentum rows in the weak form, and its "
                "volume term together with a source term is not done yet"
            )
        momentum = self.momentumColumns()

        def structures(name, keep, dims=(0, 1, 2), shape=None, factor=None):
            # one immediate per coefficient that has entries in the columns
            # `keep` takes, as _configurePointOperator states them
            result = {}
            for a in range(count):
                values = {}
                for entry in entries:
                    if (
                        entry.coefficient != a
                        or entry.dim not in dims
                        or not keep(entry.column)
                    ):
                        continue
                    idx = (
                        (entry.dim, entry.row, entry.column)
                        if shape is None
                        else (entry.row, entry.column)
                    )
                    scale = 1.0 if factor is None else factor(entry)
                    values[idx] = repr(
                        float(values.get(idx, 0.0)) + entry.factor * scale
                    )
                if values:
                    result[a] = Tensor(
                        f"{name}({a})",
                        (3,) + tuple(starSpp.shape) if shape is None else shape,
                        spp=values,
                        addressing=AddressingMode.IMMEDIATE,
                    )
            return result

        self.splitCoefficients = [
            Tensor(f"splitCoefficients({a})", (npoints,)) for a in range(count)
        ]
        # the structures without the metric, which the kernel applies to the
        # field in physical coordinates (_configurePointOperator states them
        # with it): the columns of the strong form, read from the derivative,
        # and the momentum columns, read from the field
        self.splitStrongStructure = structures(
            "splitStrongStructure", lambda column: column not in momentum
        )
        self.splitWeakStructure = structures(
            "splitWeakStructure", lambda column: column in momentum
        )
        quantities, columns = starSpp.shape

        def rowsOf(parts):
            rows = np.zeros(quantities, dtype=bool)
            for structure in parts.values():
                rows |= structure.spp().as_ndarray().any(axis=(0, 2))
            return rows

        def columnsOf(parts):
            used = np.zeros(columns, dtype=bool)
            for structure in parts.values():
                used |= structure.spp().as_ndarray().any(axis=(0, 1))
            return used

        def inputs(name, rows):
            # the quantities a part reads, so that the field is read at the
            # points for those alone
            return Tensor(
                name,
                (quantities, quantities),
                spp={(q, q): "1.0" for q in np.nonzero(rows)[0]},
                addressing=AddressingMode.IMMEDIATE,
            )

        def over(shape, last):
            # every entry of the leading indices, the last one where `last` is
            return np.broadcast_to(last, tuple(shape)).copy()

        strongRows = rowsOf(self.splitStrongStructure)
        weakRows = rowsOf(self.splitWeakStructure)
        strongColumns = columnsOf(self.splitStrongStructure)
        weakColumns = columnsOf(self.splitWeakStructure)
        self.splitStrongInputs = inputs("splitStrongInputs", strongRows)
        self.splitWeakInputs = inputs("splitWeakInputs", weakRows)
        self.splitGradient = self.nodalTemporary(
            "splitGradient",
            (npoints, quantities),
            over((npoints, quantities), strongRows),
        )
        self.splitPhysical = self.nodalTemporary(
            "splitPhysical",
            (3, npoints, quantities),
            over((3, npoints, quantities), strongRows),
        )
        self.splitRate = self.nodalTemporary(
            "splitRate", (npoints, columns), over((npoints, columns), strongColumns)
        )
        self.splitValues = self.nodalTemporary(
            "splitValues", (npoints, quantities), over((npoints, quantities), weakRows)
        )
        self.splitFlux = self.nodalTemporary(
            "splitFlux", (3, npoints, columns), over((3, npoints, columns), weakColumns)
        )
        self.splitFluxRef = self.nodalTemporary(
            "splitFluxRef", (npoints, columns), over((npoints, columns), weakColumns)
        )

        # The normal flux of the cell's own trace in the momentum rows, in the
        # coordinates of a face, where the normal is the first direction: the
        # entries of that direction, with the Voigt weights of the flux
        # (_configureNodalFlux) folded in.
        weights = voigt_weights(self.quantityBlocks())
        extended = self.numExtendedQuantities()
        self.ownTractionStructure = structures(
            "ownTractionStructure",
            lambda column: column in momentum,
            dims=(0,),
            shape=(extended, extended),
            factor=lambda entry: weights[entry.row],
        )
        self._ownTractionIndices = sorted(self.ownTractionStructure)
        faceNodes = self.numFaceFluxPoints()
        self.ownTractionCoefficients = {
            a: Tensor(f"ownTractionCoefficients({a})", (faceNodes,))
            for a in self._ownTractionIndices
        }
        # the field reaches the face with the rotation of every quantity of a
        # block that the own traction reads, and with no other one: the inverse
        # Voigt weights of the flux restricted to those blocks
        read = set()
        for structure in self.ownTractionStructure.values():
            read |= set(np.nonzero(structure.spp().as_ndarray().any(axis=1))[0])
        blocks = [
            block
            for block in self.quantityBlocks()
            if read & set(range(block.offset, block.offset + block.extent))
        ]
        self.ownTractionWeights = Tensor(
            "ownTractionWeights",
            (self.numQuantities(), extended),
            spp={
                (q, q): repr(1.0 / weights[q])
                for block in blocks
                for q in range(block.offset, block.offset + block.extent)
            },
            addressing=AddressingMode.IMMEDIATE,
        )
        # the temporaries of nodalFlux, with the entries the own traction has
        readRows = np.zeros(extended, dtype=bool)
        for block in blocks:
            readRows[block.offset : block.offset + block.extent] = True
        writes = np.zeros(extended, dtype=bool)
        for column in momentum:
            writes[column] = True
        shape = (faceNodes, extended)
        self.ownTractionTemporaries = (
            self.nodalTemporary("ownTractionValues", shape, over(shape, readRows)),
            self.nodalTemporary("ownTractionRotated", shape, over(shape, readRows)),
            self.nodalTemporary("ownTractionProduct", shape, over(shape, writes)),
            self.nodalTemporary("ownTractionBack", shape, over(shape, writes)),
        )

    def ownTractionCoefficientIndices(self):
        """The scalars of the operator without material the own traction of a
        face reads at each of its points (_configureCurvedSplit), as indices
        into the scalars the solver's operator is linear in. Empty where no
        cell may be curved."""
        return list(getattr(self, "_ownTractionIndices", []))

    def splitVolumeStatements(self, source, target):
        """The volume term of a curved cell in the split form, added to
        `target`, from the time integral `source` (_configureCurvedSplit).

        The strong columns: the derivative of the field at the points in each
        reference direction, turned into the physical one with the rows of the
        inverse Jacobian, against the structure and the coefficients (J times
        those of the operator without material), projected back. The momentum
        columns: the field at the points against the structure and the
        coefficients, the flux of each physical direction, turned into the one
        of each reference direction, projected back and tested with the
        derivative of the test functions of that direction (kDivM)."""
        statements = []
        physical = self.definedOnce(self.splitPhysical)
        for dim in range(3):
            gradient = self.definedOnce(self.splitGradient)
            statements.append(
                gradient["nq"]
                <= self.operatorEval["nk"]
                * self.db.kDivMT[dim][self.t("kl")]
                * source["lr"]
                * self.splitStrongInputs["rq"]
            )
            term = self.referenceGradients[dim]["nj"] * gradient["nq"]
            statements.append(
                physical["jnq"] <= (term if dim == 0 else physical["jnq"] + term)
            )
        rate = self.definedOnce(self.splitRate)
        terms = [
            self.splitCoefficients[a]["n"] * physical["jnq"] * structure["jqp"]
            for a, structure in self.splitStrongStructure.items()
        ]
        for position, term in enumerate(terms):
            statements.append(
                rate["np"] <= (term if position == 0 else rate["np"] + term)
            )
        statements.append(
            target["kp"] <= target["kp"] + self.operatorProject["kn"] * rate["np"]
        )

        values = self.definedOnce(self.splitValues)
        statements.append(
            values["nq"]
            <= self.operatorEval["nk"] * source["kr"] * self.splitWeakInputs["rq"]
        )
        flux = self.definedOnce(self.splitFlux)
        terms = [
            self.splitCoefficients[a]["n"] * values["nq"] * structure["jqp"]
            for a, structure in self.splitWeakStructure.items()
        ]
        for position, term in enumerate(terms):
            statements.append(
                flux["jnp"] <= (term if position == 0 else flux["jnp"] + term)
            )
        for dim in range(3):
            fluxOfDirection = self.definedOnce(self.splitFluxRef)
            statements.append(
                fluxOfDirection["np"]
                <= self.referenceGradients[dim]["nj"] * flux["jnp"]
            )
            statements.append(
                target["kp"]
                <= target["kp"]
                + self.db.kDivM[dim][self.t("kl")]
                * self.operatorProject["ln"]
                * fluxOfDirection["np"]
            )
        return statements

    def ownTractionStatements(self, face, side, source, target):
        """The normal flux of the cell's own trace `source` on the face `face`,
        in the momentum rows, added to `target` (_configureCurvedSplit): read
        at the points of the face as the corrector reads the jump there, from
        the side of the face the cell is on."""
        return self.nodalFlux(
            source,
            target,
            self.traceToFace[side][self.t("nl")],
            self.faceLift(face),
            [self.ownTractionCoefficients[a] for a in self._ownTractionIndices],
            structures=[self.ownTractionStructure[a] for a in self._ownTractionIndices],
            voigt=self.ownTractionWeights,
            temporaries=self.ownTractionTemporaries,
        )

    def faceRotation(self):
        """The rotation a face's flux is applied with: one per face, or one per
        point where a face may be curved."""
        return self.TNodes if getattr(self, "curvilinear", False) else self.T

    def faceRotationIndices(self):
        """How the rotation is contracted into the face coordinates and back:
        a curved face has a rotation per point, which is an index of its own."""
        return ("nqk", "npl") if getattr(self, "curvilinear", False) else ("qk", "pl")

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
        samples = self.materialEval.shape()[0]
        npoints = self.operatorEval.shape()[0]
        self.sourceCoefficients = [
            Tensor(f"sourceCoefficients({a})", (samples,)) for a in range(count)
        ]
        # the field at the operator's points, in whatever indices the source term
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
            Tensor(f"sourceDeviation({a})", (samples,))
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
        holds those scalars per point instead of a matrix. In global coordinates
        it is not: rotated, the same operator occupies every entry and spans
        far more than ten dimensions, so the rotation belongs in the kernel and
        not in what a face stores.

        Only the forward rotation is stored. Its inverse follows from it by the
        Voigt weights, and neither weight costs a multiplication: one half is
        folded into the structure the coefficients scale, the other is a
        constant diagonal the field passes through on its way to the face.
        """
        faceNodes = self.numFaceFluxPoints()
        quantities = self.numQuantities()
        extended = self.numExtendedQuantities()
        weights = voigt_weights(self.quantityBlocks())

        # the Rusanov penalty sits on the diagonal the Godunov state stores
        names, self.fluxSources, entries = coefficients.flux_decomposition(
            self.extendedBlocks(), self.rusanovDiagonal()
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
        shape = (faceNodes, extended)
        self.faceValues = self.nodalTemporary("faceValues", shape)
        self.faceRotated = self.nodalTemporary("faceRotated", shape)
        self.faceProduct = self.nodalTemporary("faceProduct", shape)
        self.faceBack = self.nodalTemporary("faceBack", shape)

    def _configureWeightAdjust(self):
        """What a cell needs to lift its fluxes without material and apply the
        inverse weight of its medium afterwards (weight-adjusted DG).

        The operator of a medium is the inverse of its weight times an operator
        without material, A_d = W^-1 B_d (coefficients.WeightEntry). A cell
        that takes the two apart lifts the flux of every face with B into a
        buffer of its own and applies W^-1 once for all of them: it reads the
        buffer at the sample points, multiplies by the inverse weight there and
        projects back, the three steps the operator takes inside the cell. Only
        a medium that states its weight has this, a medium with relaxation the
        one of its base material.

        The buffer has the shape of what the faces lift into, the field a
        solver adds its flux to before anything else: Q, or the extended field
        of a solver that keeps the memory variables apart. The weight reaches
        the quantities of the medium itself. The rows the faces lift beyond
        those -- the memory variables, or the strain rate a solver hands to
        them -- carry no material, so they pass through as they are.
        """
        self._weightStructure = {}
        decomposition = coefficients.weighted(self.name())
        if decomposition is None or not self.nodalMaterial:
            return
        weighted = total_extent(layout(self.primaryGroups()))
        lifted = self.numExtendedQuantities()
        self._weightedQuantities = weighted
        for a in sorted({entry.coefficient for entry in decomposition.weight}):
            self._weightStructure[a] = Tensor(
                f"weightStructure({a})",
                (weighted, weighted),
                spp={
                    (entry.row, entry.column): repr(float(entry.factor))
                    for entry in decomposition.weight
                    if entry.coefficient == a
                },
                addressing=AddressingMode.IMMEDIATE,
            )
        # what the faces of a cell lift without material, before the weight
        self.fluxBuffer = OptionalDimTensor(
            "fluxBuffer",
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            (self.num3DBasisFunctions(), lifted),
            alignStride=True,
        )
        npoints = self.operatorEval.shape()[0]
        self.weightValues = self.nodalTemporary("weightValues", (npoints, weighted))
        self.weightProduct = self.nodalTemporary("weightProduct", (npoints, weighted))
        self._weightInject = None
        self._weightPassThrough = None
        if lifted > weighted:
            self._weightInject = Tensor(
                "weightInject",
                (weighted, lifted),
                spp={(q, q): "1.0" for q in range(weighted)},
                addressing=AddressingMode.IMMEDIATE,
            )
            self._weightPassThrough = Tensor(
                "weightPassThrough",
                (lifted, lifted),
                spp={(q, q): "1.0" for q in range(weighted, lifted)},
                addressing=AddressingMode.IMMEDIATE,
            )

    #: Below this magnitude relative to its largest entry, an entry of a trace
    #: matrix is rounding of the matrix files it is formed from. The rounding
    #: reaches 5e-12 of the largest entry in the Bernstein basis at order 8,
    #: where every true entry of a trace is one; the smallest true entry in the
    #: Dubiner basis is above 1e-6 of the largest up to order 8.
    TraceZeroTolerance = 1e-10

    def tracesSelectModes(self):
        """Whether the trace of every function of the basis on a face is a
        function of the face basis or zero, as in the Bernstein basis: then
        every row of rT(i) has a single entry, and every mode of a face trace
        is one mode of the cell."""
        for i in range(4):
            tensor = self.db.rT[i]
            trace = dense(tensor, self.transpose(tensor.name()))
            if np.any(np.count_nonzero(trace, axis=1) != 1):
                return False
        return True

    def traceMatrixArgs(self, name):
        """The tensor arguments of a trace matrix: where the traces select modes
        (tracesSelectModes), a sparse layout of exactly its entries -- an
        aligned one would pad every column to a block of rows, which for a
        matrix of one entry per row is most of the bounding box again -- and a
        dense one otherwise, where a trace fills its bounding box. The old GPU
        interface (gemmforge and chainforge) reads its operands as dense, so a
        GPU build served by it keeps the dense layout, as for the face
        rotation."""
        oldGpu = self.kwargs.get("old_gpu_interface", True) and "gpu" in (
            self.kwargs.get("targets") or []
        )
        if self.tracesSelectModes() and not oldGpu:
            return {"memoryLayoutClass": CSCMemoryLayout, "alignStride": False}
        return {"alignStride": self.alignStride(name)}

    def _configureTraces(self):
        """The face traces the corrector of a cell works from.

        After the predictor, a cell writes the trace of its time integral on
        each of its faces, in the face basis: rT(i), the restriction of the
        field to face i, which is exact since the face basis spans the
        polynomials of the degree of the cell. The corrector of a face reads
        that trace and the one the neighbour wrote on the same face, applies the
        flux to the two and lifts the result into the cell; it reads no field of
        the volume. Where a neighbour takes another time step than the cell, it
        writes the trace of its time derivatives instead, and the cell evaluates
        that over its own interval (addTraces).

        The two cells sharing a face parametrise it differently: each sees the
        other's parametrisation reflected, by R = M2^-1 fP(0), which is its own
        inverse (the face orientation index is zero on every interior face).
        Each face therefore has one parametrisation both of its traces are
        written in, the one of the cell the setup calls side 0 of the face; the
        other side writes its trace reflected, R rT(i). The corrector forms the
        jump of the two traces as they are and reads it at its own points, which
        for side 1 means through the reflection: V2m R, a renumbering of the
        rows of V2m. Where the flux is one matrix per side, side 1 lifts through
        it, rDivM(i) fP(0) instead of rDivM(i) M2.

        The traces carry no padding: a face has as many modes as its basis,
        which is what the cells exchange.

        Where the trace selects modes (tracesSelectModes), so does R rT(i), R
        renumbering the face modes, and the trace matrices take a sparse layout
        (traceMatrixArgs): writing a trace gathers modes of the cell rather than
        multiplying the bounding box of a matrix of zeros and ones. Their
        entries below TraceZeroTolerance of the largest are rounding of the
        matrix files they are formed from, and are dropped.
        """
        modes = self.num2DBasisFunctions()
        self.traceMatrix = tensor_collection_from_constant_expression(
            "faceTraceMatrix",
            lambda i, side: (
                self.db.rT[i][self.t("ml")]
                if side == 0
                else self.db.M2inv["mn"] * self.db.fPrT[i][self.t("nl")]
            ),
            simpleParameterSpace(4, 2),
            target_indices=self.t("ml"),
            tensor_args=self.traceMatrixArgs("faceTraceMatrix"),
            zero_tolerance=self.TraceZeroTolerance,
        )["faceTraceMatrix"]
        shape = (modes, self.numQuantities())
        self.faceTrace = OptionalDimTensor(
            "faceTrace", self.Q.optName(), self.Q.optSize(), self.Q.optPos(), shape
        )
        self.faceTraceNeighbor = OptionalDimTensor(
            "faceTraceNeighbor",
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            shape,
        )
        if self.nodalFaceFlux:
            # the jump read at the points of the face, from either side
            self.traceToFace = tensor_collection_from_constant_expression(
                "traceToFace",
                lambda side: (
                    self.db.V2mToFaceQuad[self.t("nm")]
                    if side == 0
                    else self.db.V2mToFaceQuad[self.t("nk")]
                    * self.db.M2inv["kj"]
                    * self.db.fP[0][self.t("jm")]
                ),
                simpleParameterSpace(2),
                target_indices=self.t("nm"),
                tensor_args={"alignStride": self.alignStride("traceToFace")},
                zero_tolerance=1e-14,
            )["traceToFace"]
            self.traceJump = self.nodalTemporary("traceJump", shape)
            # one set of scalars per face: the corrector applies them to the jump
            faceNodes = self.numFaceFluxPoints()
            self.fluxCoefficients = [
                Tensor(f"fluxCoefficients({a})", (faceNodes,))
                for a in range(len(self.fluxStructure))
            ]
        else:
            self.traceLift = tensor_collection_from_constant_expression(
                "traceLift",
                lambda i, side: self.db.rDivM[i][self.t("km")]
                * (self.db.M2["mn"] if side == 0 else self.db.fP[0][self.t("mn")]),
                simpleParameterSpace(4, 2),
                target_indices=self.t("kn"),
                tensor_args={"alignStride": self.alignStride("traceLift")},
                zero_tolerance=1e-14,
            )["traceLift"]
            self.traceProduct = self.nodalTemporary(
                "traceProduct", (modes, self.numExtendedQuantities())
            )
            # the jump of the two traces, which the corrector applies the flux to
            self.traceJump = self.nodalTemporary("traceJump", shape)
        if any(self.fusesFaces(target) for target in self.kwargs.get("targets") or []):
            self._configureTraceSides()

    def fusesFaces(self, target):
        """Whether the face kernels of a target take the four faces of a cell
        in one kernel, with the side of each face as data
        (_configureTraceSides): on the device, where a launch costs more than
        applying the reflection of a side does. The host keeps a kernel per
        face and side, where the reflection is folded into the operators."""
        return target == "gpu"

    def traceSideArgs(self, name):
        """The tensor arguments of the side operators (_configureTraceSides): a
        sparse layout of the pattern they share, which the old GPU interface
        reads as dense (see traceMatrixArgs)."""
        oldGpu = self.kwargs.get("old_gpu_interface", True) and "gpu" in (
            self.kwargs.get("targets") or []
        )
        if oldGpu:
            return {"alignStride": self.alignStride(name)}
        return {"memoryLayoutClass": CSCMemoryLayout, "alignStride": False}

    def _configureTraceSides(self):
        """The side of a face as data, for the kernels that take the four faces
        of a cell at once (fusesFaces).

        The operators of side 1 of a face are those of side 0 with the
        reflection R = M2^-1 fP(0) applied (see _configureTraces): the trace
        R rT(i), the lift L(i, 0) R, the points V(0) R. A kernel over the four
        faces of a cell applies the operators of side 0 and, between the trace
        and what reads it, the operator of the side of each face, I or R, which
        it reads where the cell points it (traceSide). The two are stored by
        the pattern they share (traceSideOperator), so that the kernel is the
        same for either: in the Bernstein basis R renumbers the face modes and
        a row of the union has two entries at most; in a hierarchical basis R
        keeps the degree of a mode, and the union is the blocks of the
        degrees, which keeps a trace that reaches the modes up to a degree
        there.

        The four faces of a kernel read their operators at once, so each has
        tensors of its own: its two traces (faceTraces, faceTracesNeighbor),
        and the operator of the corrector in the shape the cell carries it,
        laid out as the one of a single face (faceAplusT, or the rotation and
        the scalars per point).
        """
        name = "traceSideOperator"
        modes = self.num2DBasisFunctions()
        reflection = dense(self.db.M2inv, self.transpose(self.db.M2inv.name())) @ dense(
            self.db.fP[0], self.transpose(self.db.fP[0].name())
        )
        scale = np.max(np.abs(reflection))
        reflection[np.abs(reflection) <= self.TraceZeroTolerance * scale] = 0.0
        operators = [np.eye(modes), reflection]
        pattern = (operators[0] != 0.0) | (operators[1] != 0.0)
        # face mode by face mode, what checkSideKeeps reads
        self._traceSidePattern = pattern
        if self.transpose(name):
            operators = [operator.T for operator in operators]
            pattern = pattern.T
        entries = list(zip(*np.nonzero(pattern)))
        args = self.traceSideArgs(name)
        self.traceSideOperator = [
            Tensor(
                f"{name}({side})",
                pattern.shape,
                spp={entry: repr(float(operator[entry])) for entry in entries},
                **args,
            )
            for side, operator in enumerate(operators)
        ]
        self.traceSide = [
            Tensor(f"traceSide({face})", pattern.shape, spp=pattern, **args)
            for face in range(4)
        ]
        shape = (modes, self.numQuantities())
        self.faceTraces, self.faceTracesNeighbor = (
            [
                OptionalDimTensor(
                    f"{base}({face})",
                    self.Q.optName(),
                    self.Q.optSize(),
                    self.Q.optPos(),
                    shape,
                )
                for face in range(4)
            ]
            for base in ("faceTraces", "faceTracesNeighbor")
        )
        if self.nodalFaceFlux:
            rotation = self.faceRotation()
            self.faceRotations = [
                cloneOperand(rotation, f"face{rotation.name()}({face})")
                for face in range(4)
            ]
            self.faceFluxCoefficients = [
                [
                    cloneOperand(coefficient, f"faceFluxCoefficients({face},{a})")
                    for a, coefficient in enumerate(self.fluxCoefficients)
                ]
                for face in range(4)
            ]
        else:
            self.faceAplusT = [
                cloneOperand(self.AplusT, f"faceAplusT({face})") for face in range(4)
            ]

    def faceOperands(self, face):
        """What the corrector of a face reads of the cell in a kernel that takes
        the four faces of a cell at once (_configureTraceSides)."""
        if self.nodalFaceFlux:
            return {
                "rotation": self.faceRotations[face],
                "coefficients": self.faceFluxCoefficients[face],
            }
        return {"fluxMatrix": self.faceAplusT[face]}

    def sideStatement(self, target, sideOperator, source, indices="mq"):
        """target = the side operator applied to source (_configureTraceSides):
        `source` in the parametrisation of side 0 of a face, `target` in the
        one of the side the operator stands for. `indices` names the face modes
        first."""
        inner = "z"
        return (
            target[indices]
            <= sideOperator[self.t(indices[0] + inner)] * source[inner + indices[1:]]
        )

    def traceCorrector(
        self,
        face,
        side,
        ownTrace=None,
        neighborTrace=None,
        sideOperator=None,
        operands=None,
    ):
        """What the corrector of one face adds to the field the faces lift into:
        the flux of the two traces of the face, lifted into the cell, less the
        normal flux of the cell's own trace, since the volume term is in the
        strong form (see volumeDerivatives).

        Where the flux is the scalars per point, it is the jump of the two
        traces that goes to the points, and the one set of scalars the face
        carries is applied to it: for a regular face the neighbour's part of
        the Godunov state, which is the local part less the identity, negated.
        A face without neighbour reads a zero trace there and carries the local
        scalars negated, so that the same statements apply its own flux.

        Where the flux is one matrix per side, the matrix of the cell's own side
        is its flux less the normal flux (toCorrectorForm); on a regular face
        that is the neighbour's matrix negated, since the two add up to the
        normal flux, and any other face reads a zero trace for the neighbour. So
        the one matrix, applied to the jump of the traces (the cell's own less
        the neighbour's), is the flux of both on every face.

        Given a `sideOperator`, the corrector takes the side as data: it
        applies the operator to the jump and goes on as side 0 does
        (_configureTraceSides). `operands` are the tensors of the face it
        reads the operator of the corrector from (faceOperands); the ones of
        the single face a kernel applies otherwise.
        """
        own = self.faceTrace if ownTrace is None else ownTrace
        neighbor = self.faceTraceNeighbor if neighborTrace is None else neighborTrace
        operands = operands or {}
        target = self.extendedQTensor()
        jump = self.definedOnce(self.traceJump)
        if self.nodalFaceFlux:
            statements = [jump["mq"] <= neighbor["mq"] - own["mq"]]
        else:
            statements = [jump["mq"] <= own["mq"] - neighbor["mq"]]
        if sideOperator is not None:
            sided = self.definedOnce(self.traceJump)
            statements.append(self.sideStatement(sided, sideOperator, jump))
            jump, side = sided, 0
        if self.nodalFaceFlux:
            return statements + self.nodalFlux(
                jump,
                target,
                self.traceToFace[side][self.t("nl")],
                self.faceLift(face),
                operands.get("coefficients", self.fluxCoefficients),
                rotation=operands.get("rotation"),
            )
        product = self.definedOnce(self.traceProduct)
        fluxMatrix = operands.get("fluxMatrix", self.AplusT)
        return statements + [
            product["mp"] <= jump["mq"] * fluxMatrix["qp"],
            target["kp"]
            <= target["kp"] + self.traceLift[face, side][self.t("km")] * product["mp"],
        ]

    def faceDerivativeTraces(self):
        """What a cell writes on a face whose neighbour evaluates it over an
        interval of its own: the traces of the time derivatives of the
        predictor, one per derivative and each with the modes its derivative
        reaches on the face. A derivative the chain narrows keeps fewer modes
        with every order, and so does its trace; together they hold as many
        modes as the cell has. A derivative kept in the basis of its degree
        (derivativesByDegree) is raised into the face basis of the cell by its
        trace (faceTraceOfDegree), and its trace has every face mode."""
        if not hasattr(self, "_faceDerivatives"):
            modes = self.num2DBasisFunctions()
            shape = (modes, self.numQuantities())
            fullShape = self.faceTrace.shape()
            self._faceDerivatives = []
            for k, derivative in enumerate(self.dQs):
                # the trace of a derivative in the basis of its degree is raised
                # into the face basis of the cell, which it fills
                reached = (
                    modes
                    if self.derivativesByDegree() and k > 0
                    else self._traceModesOf(derivative)
                )
                rows = self.faceTrace.insertOptDim(
                    (slice(0, reached), slice(None)), (slice(None),)
                )
                pattern = np.zeros(fullShape, dtype=bool)
                pattern[rows] = True
                self._faceDerivatives.append(
                    OptionalDimTensor(
                        f"faceDerivatives({k})",
                        self.Q.optName(),
                        self.Q.optSize(),
                        self.Q.optPos(),
                        shape,
                        spp=None if reached == modes else pattern,
                    )
                )
        return self._faceDerivatives

    def _traceModesOf(self, field):
        """How many face modes the trace of a field with the pattern of `field`
        reaches: the face basis is ordered by degree, and the trace of the
        modes of a degree reaches the face modes up to that degree."""
        pattern = field.spp().as_ndarray()
        if isinstance(field, OptionalDimTensor) and field.hasOptDim():
            pattern = pattern.any(axis=field.optPos())
        volumeModes = int(np.nonzero(pattern.any(axis=1))[0].max()) + 1
        reached = 0
        for matrix in self.traceMatrix.values():
            values = matrix.spp().as_ndarray()
            values = values if not self.transpose(matrix.name()) else values.T
            rows = np.nonzero(values[:, :volumeModes].any(axis=1))[0]
            if rows.size > 0:
                reached = max(reached, int(rows.max()) + 1)
        return reached

    def faceDerivativeStatements(self, face, side, sideOperator=None):
        """The traces of the time derivatives on one face, see
        faceDerivativeTraces. Given a `sideOperator`, the traces of side 0
        with the side of the face as data (_configureTraceSides)."""
        statements = []
        for k, (trace, derivative) in enumerate(
            zip(self.faceDerivativeTraces(), self.dQs)
        ):
            matrix = (
                self.faceTraceOfDegree(face, side, self.order - 1 - k)
                if self.derivativesByDegree() and k > 0
                else self.traceMatrix[face, side]
            )
            if sideOperator is None:
                statements.append(
                    trace["mq"] <= matrix[self.t("ml")] * derivative["lq"]
                )
                continue
            self.checkSideKeeps(trace)
            gathered = self.definedOnce(trace)
            statements.append(gathered["mq"] <= matrix[self.t("ml")] * derivative["lq"])
            statements.append(self.sideStatement(trace, sideOperator, gathered))
        return statements

    def checkSideKeeps(self, trace):
        """Refuses a trace that reaches the face modes up to a degree where the
        side operators do not keep it there (_configureTraceSides): the trace
        of the other side has to have the pattern of the trace of side 0."""
        pattern = trace.spp().as_ndarray()
        if isinstance(trace, OptionalDimTensor) and trace.hasOptDim():
            pattern = pattern.any(axis=trace.optPos())
        reached = (
            int(np.nonzero(pattern.any(axis=tuple(range(1, pattern.ndim))))[0].max())
            + 1
        )
        if self._traceSidePattern[reached:, :reached].any():
            raise ValueError(
                f"the side operators take the first {reached} face modes of "
                f"{trace.name()} beyond them"
            )

    def faceTracesStatements(self):
        """The traces of the time integral of a cell on its four faces, each in
        the parametrisation of the side its face takes (_configureTraceSides):
        the trace of side 0, then the side operator of the face."""
        statements = []
        for face in range(4):
            gathered = self.definedOnce(self.faceTraces[face])
            statements.append(
                gathered["mq"] <= self.traceMatrix[face, 0][self.t("ml")] * self.I["lq"]
            )
            statements.append(
                self.sideStatement(
                    self.faceTraces[face], self.traceSide[face], gathered
                )
            )
        return statements

    def faceCorrectorsStatements(self):
        """The correctors of the four faces of a cell (traceCorrector), each
        with the side of its face as data and the operands of its face
        (_configureTraceSides)."""
        statements = []
        for face in range(4):
            statements += self.traceCorrector(
                face,
                0,
                ownTrace=self.faceTraces[face],
                neighborTrace=self.faceTracesNeighbor[face],
                sideOperator=self.traceSide[face],
                operands=self.faceOperands(face),
            )
        return statements

    def faceTraceAtTime(self, coefficients):
        """The trace of the field over an interval, or at an instant, from the
        traces of its time derivatives: the Taylor expansion the predictor
        evaluates in the volume, one scalar of `coefficients` per derivative,
        the same the volume is evaluated with."""
        expansion = None
        for coefficient, trace in zip(coefficients, self.faceDerivativeTraces()):
            term = coefficient * trace["mq"]
            expansion = term if expansion is None else expansion + term
        return expansion

    def faceDerivativesAtTimeStatements(self):
        """The trace over an interval, or at an instant, into the trace a
        corrector reads as the neighbour's (faceTraceAtTime), with the scalars
        power(k)."""
        powers = [Scalar(f"power({k})") for k in range(self.order)]
        return [self.faceTraceNeighbor["mq"] <= self.faceTraceAtTime(powers)]

    def addTraces(self, generator, targets):
        """The kernels of the face traces (see _configureTraces): the trace of
        the time integral of a cell on one of its faces, the traces of its time
        derivatives on one face, their evaluation over an interval, and the
        corrector of one face. A family is indexed by the face and by the side
        of the face the cell is on.

        A target that takes the side as data (fusesFaces) has a kernel for
        the traces of the time integral on the four faces of a cell
        (faceTraces), one for the correctors of the four (faceCorrectors), and
        a family of the traces of the time derivatives indexed by the face
        alone."""
        for target in targets:
            prefix = generate_kernel_name_prefix(target)
            with self.singleDefinitions(target == "gpu"):
                if self.fusesFaces(target):
                    generator.add(
                        f"{prefix}faceTraces",
                        self.faceTracesStatements(),
                        target=target,
                    )
                    generator.addFamily(
                        f"{prefix}faceDerivatives",
                        simpleParameterSpace(4),
                        lambda face: self.faceDerivativeStatements(
                            face, 0, self.traceSide[face]
                        ),
                        target=target,
                    )
                    generator.add(
                        f"{prefix}faceCorrectors",
                        self.faceCorrectorsStatements(),
                        target=target,
                    )
                    generator.add(
                        f"{prefix}faceDerivativesAtTime",
                        self.faceDerivativesAtTimeStatements(),
                        target=target,
                    )
                    continue
                generator.addFamily(
                    f"{prefix}faceTrace",
                    simpleParameterSpace(4, 2),
                    lambda i, side: self.faceTrace["mq"]
                    <= self.traceMatrix[i, side][self.t("ml")] * self.I["lq"],
                    target=target,
                )
                generator.addFamily(
                    f"{prefix}faceDerivatives",
                    simpleParameterSpace(4, 2),
                    self.faceDerivativeStatements,
                    target=target,
                )
                generator.add(
                    f"{prefix}faceDerivativesAtTime",
                    self.faceDerivativesAtTimeStatements(),
                    target=target,
                )
                generator.addFamily(
                    f"{prefix}faceCorrector",
                    simpleParameterSpace(4, 2),
                    self.traceCorrector,
                    target=target,
                )

    def weightAdjustable(self):
        """Whether a cell of this build can apply the inverse weight of its
        medium after the lift of its fluxes, see _configureWeightAdjust."""
        return bool(getattr(self, "_weightStructure", None))

    def weightAdjust(self, source, target):
        """Adds the inverse weight applied to `source` to `target`: `source` read
        at the sample points, multiplied by W^-1 there and projected back, and
        the rows beyond the medium's own as they are."""
        values = self.definedOnce(self.weightValues)
        product = self.definedOnce(self.weightProduct)
        weighted = self._weightedQuantities
        read = source["kq"]
        if self._weightInject is not None:
            read = read.subslice("q", 0, weighted)
        # a curved cell weighs with J W, so 1 / J comes on top at every point
        # (_configureCurvilinear)
        scale = getattr(self, "weightScale", None)
        read = self.operatorEval["nk"] * read
        if scale is not None:
            read = scale["n"] * read
        statements = [values["nq"] <= read]
        for position, (a, structure) in enumerate(self._weightStructure.items()):
            term = self.nodalCoefficients[a]["n"] * values["nq"] * structure["qp"]
            statements.append(
                product["np"] <= (term if position == 0 else product["np"] + term)
            )
        if self._weightInject is None:
            statements.append(
                target["kp"]
                <= target["kp"] + self.operatorProject["kn"] * product["np"]
            )
        else:
            statements.append(
                target["kp"]
                <= target["kp"]
                + self.operatorProject["kn"] * product["nq"] * self._weightInject["qp"]
                + source["kq"] * self._weightPassThrough["qp"]
            )
        return statements

    def nodalFlux(
        self,
        source,
        target,
        toFace,
        lift,
        coefficientsOfFace,
        rotation=None,
        structures=None,
        voigt=None,
        temporaries=None,
    ):
        """One face contribution where the operator varies along the face.

        The field is read at the points of the face and turned into the face
        coordinates the scalars are stated in, the operator is applied there,
        and the result is turned back, integrated against the face basis with
        the weights of the rule the points are of, and lifted into the cell.
        The rotation is the same matrix both ways, once transposed against the
        quantity the field carries and once against the quantity the result is
        written in.

        `toFace` reads the field at the points of the face, with the point index
        first and the mode index second; `lift` goes the other way. Both come
        indexed, because how a matrix is laid out is the caller's to state.
        `rotation` is the tensor the rotation is read from, the one of the
        single face a kernel applies where not given (faceRotation).
        `structures` are the fixed entries the scalars scale, one per scalar,
        those of the flux where not given, and `voigt` the inverse Voigt
        weights the field reaches the face with, all of them where not given.
        `temporaries` are the four the field passes through, the ones of the
        flux where not given.
        """
        rotation = self.faceRotation() if rotation is None else rotation
        structures = self.fluxStructure if structures is None else structures
        voigt = self.inverseVoigtWeights if voigt is None else voigt
        into, back = self.faceRotationIndices()
        if temporaries is None:
            temporaries = (
                self.faceValues,
                self.faceRotated,
                self.faceProduct,
                self.faceBack,
            )
        faceValues, faceRotated, faceProduct, faceBack = (
            self.definedOnce(temporary) for temporary in temporaries
        )
        statements = [
            faceValues["nq"] <= toFace * source["lk"] * voigt["kq"],
            faceRotated["nk"] <= faceValues["nq"] * rotation[into],
        ]
        first = True
        for coefficient, structure in zip(coefficientsOfFace, structures):
            term = coefficient["n"] * faceRotated["nk"] * structure["kl"]
            statements.append(
                faceProduct["nl"] <= (term if first else faceProduct["nl"] + term)
            )
            first = False
        statements.append(faceBack["np"] <= faceProduct["nl"] * rotation[back])
        statements.append(target["kp"] <= target["kp"] + lift * faceBack["np"])
        return statements

    def solverCoefficientCount(self):
        """How many scalars the operator this solver applies is linear in."""
        return getattr(self, "_solverCoefficientCount", 0)

    def solverCoefficientOrigins(self):
        """Where each of those scalars comes from -- the material, or the run."""
        return getattr(self, "_solverCoefficientOrigins", [])

    def solverCoefficientEntries(self):
        """Where each of those scalars sits, per direction."""
        return getattr(self, "_solverCoefficientEntries", [])

    def faultFluxCoefficients(self):
        """The scalars the lift of a fault face reads at each of its points.

        The lift applies the coefficient matrix of the fault normal to the
        imposed state, which the fault states in its own coordinates; there the
        normal is the first direction, so the matrix is the star of that
        direction, and the scalars are the ones with an entry in it. Where the
        material varies inside a cell it varies along a fault too, and a face
        then carries these scalars at its points in place of the one matrix per
        side. Empty where the material does not vary, and the face keeps the
        matrix. Given as indices into the scalars the solver's operator is
        linear in.
        """
        if not getattr(self, "nodalMaterial", False):
            return []
        return sorted(
            {
                entry.coefficient
                for entry in self.solverCoefficientEntries()
                if entry.dim == 0
            }
        )

    def nodalTemporary(self, name, shape, pattern=None):
        """A temporary of the nodal path.

        It carries the same field a kernel's operands do, so a build that fuses
        simulations gives it that index as well; everything else about the
        nodal path is per cell and shared across them. `pattern`, over `shape`,
        names the entries it has, where those are fewer than all."""
        spp = None
        if pattern is not None:
            spp = np.asarray(pattern, dtype=bool)
            if self.Q.optSize() > 1:
                spp = np.repeat(
                    np.expand_dims(spp, self.Q.optPos()),
                    self.Q.optSize(),
                    axis=self.Q.optPos(),
                )
        return OptionalDimTensor(
            name,
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            shape,
            spp=spp,
            temporary=True,
        )

    @contextmanager
    def singleDefinitions(self, enabled=True):
        """Lets every nodal temporary a kernel writes be written once.

        The nodal forms write their temporaries again and again inside one
        kernel: the field at the samples once per direction, the product once
        per application of the operator, the values at a face once per face.
        The write that starts each round covers only what that round has in
        its sparsity pattern -- the samples a direction reaches, the columns
        of the first coefficient of a face -- and leaves the rest to be zero.
        The host generator clears that rest at every such write. The device
        generator clears what no operation of the kernel writes, judged over
        the whole kernel, so a round reads what an earlier round, or another
        temporary sharing the buffer, left there. Inside this block every
        round writes a temporary of its own, once, which the device generator
        clears for what that one write leaves out.
        """
        previous = getattr(self, "_singleDefinitions", False)
        self._singleDefinitions = enabled
        try:
            yield
        finally:
            self._singleDefinitions = previous

    def definedOnce(self, prototype):
        """`prototype`, or a temporary of its shape and sparsity pattern that
        nothing else writes, where the kernel being built defines each
        temporary once."""
        if not getattr(self, "_singleDefinitions", False):
            return prototype
        self._definitionCount = getattr(self, "_definitionCount", 0) + 1
        # a member of a family carries its group in its name, which a name of
        # its own cannot; the group goes into the base name instead
        name = "".join(
            [prototype.baseName()]
            + [f"_{index}" for index in prototype.group()]
            + [f"Once{self._definitionCount}"]
        )
        if not isinstance(prototype, OptionalDimTensor):
            # per cell, so shared by the simulations a build bundles
            return Tensor(name, prototype.shape(), spp=prototype.spp(), temporary=True)
        return OptionalDimTensor(
            name,
            prototype.optName(),
            prototype.optSize(),
            prototype.optPos(),
            self._fieldShape(prototype),
            spp=prototype.spp(),
            temporary=True,
        )

    def nodalAssembly(self):
        """Folds the Jacobian rows into the operator, once per kernel.

        On a straight-sided cell they go into the constant part
        (_configureLeftFactor). Where the build asks for the assembled form, the
        coefficients go in as well and what comes out is what the groups give at
        every sample point. That trades the products a kernel does at every
        application for a temporary over the points, so which one is cheaper
        depends on how often the kernel applies the operator and on the
        machine. In a curved cell the metric varies from point to point as the
        material does, and the two go into the operator at every point
        together, which is the assembled form.
        """
        if not self.nodalMaterial:
            return []
        if self.leftFactor:
            statements = [
                self.leftUnitFolded[dim]["qu"]
                <= self.referenceGradients[dim]["j"] * self.leftUnit["jqu"]
                for dim in range(3)
            ]
            if self.nodalOperatorAssembled:
                folded = None
                for a, coefficient in enumerate(self.nodalCoefficients):
                    term = coefficient["n"] * self.leftStructure[a]["up"]
                    folded = term if folded is None else folded + term
                statements.append(self.leftAtPoint["nup"] <= folded)
            return statements
        statements = []
        for dim in range(3):
            folded = None
            for a, coefficient in enumerate(self.nodalCoefficients):
                term = (
                    coefficient["n"]
                    * self.referenceGradients[dim]["nj"]
                    * self.coefficientStructure[a]["jqp"]
                )
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
        project=None,
        pointSum=None,
    ):
        """One application of the operator where the material varies inside the
        cell: read the derivative at the sample points, multiply by the
        material there, and project the result back.

        `operators` gives the modal operator per direction -- the stiffness for
        a derivative step, whatever the caller needs otherwise. `spectator`
        names indices the operator leaves alone, for a field that carries more
        than modes and quantities; the caller then hands over the two
        temporaries those indices widen. On a straight-sided cell the field is
        read at the points after the constant part of the operator, through
        temporaries over the groups (_configureLeftFactor), and of the two only
        the product is used. `scalar` scales the result, and `accumulate` adds
        it to what the target holds instead of replacing it. `project` is the
        projection the result comes back with, for a target that keeps fewer
        modes than the field; the operator's own where not given. `pointSum`,
        a temporary over the points with a weight and whether it starts the
        sum, collects the product at the points as well, for a caller that
        projects a weighted sum of several applications at once.
        """
        valuesPrototype, product = (
            temporaries
            if temporaries is not None
            else (self.nodalValues, self.nodalProduct)
        )
        product = self.definedOnce(product)
        if self.leftFactor:
            statements = self._leftProduct(source, operators, spectator, product)
        else:
            statements = []
            for dim in range(3):
                values = self.definedOnce(valuesPrototype)
                statements.append(
                    values["nq" + spectator]
                    <= self.operatorEval["nk"]
                    * operators[dim][self.t("kl")]
                    * source["lq" + spectator]
                )
                term = values["nq" + spectator] * self.starAtPoint[dim]["nqp"]
                statements.append(
                    product["np" + spectator]
                    <= (term if dim == 0 else product["np" + spectator] + term)
                )
        if pointSum is not None:
            total, weight, starts = pointSum
            weighted = weight * product["np" + spectator]
            statements.append(
                total["np" + spectator]
                <= (weighted if starts else total["np" + spectator] + weighted)
            )
        projector = self.operatorProject if project is None else project
        projected = projector["kn"] * product["np" + spectator]
        if scalar is not None:
            projected = scalar * projected
        statements.append(
            target["kp" + spectator]
            <= (target["kp" + spectator] + projected if accumulate else projected)
        )
        return statements

    def _leftProduct(self, source, operators, spectator, product):
        """The operator times the derivative at the sample points of a
        straight-sided cell (_configureLeftFactor): the derivatives of the
        three directions combined into the groups in modal space, the groups
        read at the points, and the material there."""
        extents = self._fieldShape(product)[2:]
        modalSum = None
        for dim in range(3):
            term = (
                operators[dim][self.t("kl")]
                * source["lq" + spectator]
                * self.leftUnitFolded[dim]["qu"]
            )
            modalSum = term if modalSum is None else modalSum + term
        modal, modalSum = self._leftModal(modalSum, spectator, extents)
        modal = self.definedOnce(modal)
        values = self.definedOnce(self._leftValues(spectator, extents))
        statements = [
            modal["ku" + spectator] <= modalSum,
            values["nu" + spectator]
            <= self.operatorEval["nk"] * modal["ku" + spectator],
        ]
        if self.nodalOperatorAssembled:
            terms = [values["nu" + spectator] * self.leftAtPoint["nup"]]
        else:
            terms = [
                coefficient["n"]
                * values["nu" + spectator]
                * self.leftStructure[a]["up"]
                for a, coefficient in enumerate(self.nodalCoefficients)
            ]
        for position, term in enumerate(terms):
            statements.append(
                product["np" + spectator]
                <= (term if position == 0 else product["np" + spectator] + term)
            )
        return statements

    def _leftModal(self, modalSum, spectator, extents):
        """The temporary the constant part of the operator writes, together with
        the sum it is written from. It has the modes the derivative matrices
        reach from the modes the field has, which is fewer for a derivative of
        the chain that keeps fewer (narrowsDerivatives) and never the highest
        degree; the field is read at the points from those alone."""
        shape = (self.num3DBasisFunctions(), self.leftGroups) + tuple(extents)
        indices = Indices(
            self.Q.insertOptDim("ku" + spectator, self.Q.optName()),
            self.Q.insertOptDim(shape, (self.Q.optSize(),)),
        )
        modalSum = DeduceIndices(indices).visit(modalSum)
        modalSum = EquivalentSparsityPattern().visit(modalSum)
        pattern = modalSum.eqspp().as_ndarray()
        temporaries = self._leftModalTemporaries
        key = (spectator, tuple(extents), pattern.tobytes())
        if key not in temporaries:
            suffix = f"_{spectator}" if spectator else ""
            index = sum(1 for other in temporaries if other[0] == spectator)
            temporaries[key] = OptionalDimTensor(
                f"leftModal{suffix}({index})",
                self.Q.optName(),
                self.Q.optSize(),
                self.Q.optPos(),
                shape,
                spp=pattern,
                temporary=True,
            )
        return temporaries[key], modalSum

    def _leftValues(self, spectator, extents):
        """The temporary the groups are read at the sample points into."""
        temporaries = self._leftValuesTemporaries
        key = (spectator, tuple(extents))
        if key not in temporaries:
            suffix = f"_{spectator}" if spectator else ""
            temporaries[key] = self.nodalTemporary(
                f"leftValues{suffix}",
                (self.materialEval.shape()[0], self.leftGroups) + tuple(extents),
            )
        return temporaries[key]

    def _fieldShape(self, tensor):
        """The shape of a tensor without the index of the simulations a build
        bundles."""
        if isinstance(tensor, OptionalDimTensor) and tensor.hasOptDim():
            return tuple(
                extent
                for position, extent in enumerate(tensor.shape())
                if position != tensor.optPos()
            )
        return tuple(tensor.shape())

    def narrowsDerivatives(self):
        """Whether the k-th time derivative of the predictor keeps only the
        modes up to the degree order - 1 - k where the operator varies inside
        the cell. Not here; a solver whose chain can says so."""
        return False

    def modalDerivative(self, source):
        """One step of a derivative chain where the operator is one per cell:
        the derivative of `source` in each direction times the star matrix of
        that direction, as statements to run first and the terms of the sum
        that is the step.

        In a hierarchical basis the step is the derivative matrices kDivMT
        against the field, and the pattern of a derivative shows which modes
        it has. In another basis, of a chain that keeps the degree (a source
        term keeps it), the derivative of the degree of the cell is of one
        degree less, and kDivMT is that lower degree raised back into the
        basis of the cell, which has every function of the cell in it. The
        step then goes through the basis of the lower degree, where the
        derivative matrices (kDivMTDegree) have a few entries per row and the
        star matrices see as many functions as that basis has, and is raised
        back once, after the star matrices of all three directions are summed
        (elevateDegree). A chain that drops the degree keeps its derivatives
        in the bases of their degrees instead (derivativesByDegree).
        """
        if self.basis.hierarchical:
            return [], [
                self.db.kDivMT[j][self.t("kl")]
                * source["lq"]
                * self.starMatrix(j)["qp"]
                for j in range(3)
            ]
        degree = self.order - 1
        if not hasattr(self, "_derivativeLower"):
            # what the star matrices give in the basis of the lower degree, in
            # the leading rows
            columns = self.starMatrix(0).shape()[1]
            self._derivativeLower = OptionalDimTensor(
                "derivativeLower",
                self.Q.optName(),
                self.Q.optSize(),
                self.Q.optPos(),
                (self.num3DBasisFunctions(), columns),
                spp=self.leadingRows(count(degree - 1), columns),
                alignStride=True,
                temporary=True,
            )
        lower = self.definedOnce(self._derivativeLower)
        lowerSum = None
        for j in range(3):
            term = (
                self.derivativeOfDegree(degree, j)[self.t("kl")]
                * source["lq"]
                * self.starMatrix(j)["qp"]
            )
            lowerSum = term if lowerSum is None else lowerSum + term
        raised = self.elevationOfDegree(degree - 1)[self.t("kl")] * lower["lp"]
        return [lower["kp"] <= lowerSum], [raised]

    def derivativesByDegree(self):
        """Whether the time derivatives of a chain that drops a degree with every
        derivative are kept in the bases of their degrees: derivative k in the
        basis of degree order - 1 - k, in the leading rows of its tensor. A
        hierarchical basis has those bases as its leading functions and keeps
        them there by the patterns of its derivatives; another basis states the
        derivatives and elevations between its degrees (kDivMTDegree,
        elevateDegree), and whatever reads a derivative goes through them."""
        return not self.basis.hierarchical and self.chainDropsDegree()

    def chainDropsDegree(self):
        """Whether every time derivative of the predictor has one degree less
        than the one before. Not here; a solver whose chain does says so."""
        return False

    def embeddingOfDegree(self, degree):
        """The basis of a degree in the basis of the cell, as a dense matrix
        (embedDegree); the identity for the degree of the cell."""
        size = self.num3DBasisFunctions()
        if degree == self.order - 1:
            return np.eye(size)
        tensor = self.db.embedDegree[degree]
        return dense(tensor, self.transpose(tensor.name()))

    def embeddingOperator(self, degree):
        """The embedding of a degree (embedDegree) as a matrix of the size of
        the cell basis, its columns in the leading ones, so that it raises a
        field kept in the leading rows of a tensor of the cell."""
        name = f"embeddingOfDegree({degree})"
        if not hasattr(self, "_embeddingOperators"):
            self._embeddingOperators = {}
        if name not in self._embeddingOperators:
            size = self.num3DBasisFunctions()
            logical = self.embeddingOfDegree(degree)
            full = np.zeros((size, size))
            full[:, : logical.shape[1]] = logical
            stored = full.T if self.transpose(name) else full
            self._embeddingOperators[name] = Tensor(
                name,
                stored.shape,
                spp={
                    index: repr(float(value))
                    for index, value in np.ndenumerate(stored)
                    if value != 0.0
                },
                alignStride=self.alignStride(name),
            )
        return self._embeddingOperators[name]

    def derivativeInCellBasis(self, k, indices, inner):
        """Time derivative k in the basis of the cell, as a term indexed by
        `indices` (the index of the basis first): the derivative itself where it
        is kept there, else the embedding of its degree times the derivative
        (derivativesByDegree), contracted over `inner`, an index the caller
        leaves free. As a term the embedding is a factor of the product it
        enters, which can apply it after the factors that leave fewer columns."""
        derivative = self.dQs[k]
        if not self.derivativesByDegree() or k == 0:
            return derivative[indices]
        embedding = self.embeddingOperator(self.order - 1 - k)
        return embedding[self.t(indices[0] + inner)] * derivative[inner + indices[1:]]

    def faceTraceOfDegree(self, face, side, degree):
        """The face trace (traceMatrix) of a field given in the basis of a degree,
        in the leading rows of a tensor of the cell: the trace times the
        embedding of that degree, into the face basis of the cell."""
        name = f"faceTraceOfDegree({face},{side},{degree})"
        if not hasattr(self, "_faceTracesOfDegree"):
            self._faceTracesOfDegree = {}
        if name not in self._faceTracesOfDegree:
            tensor = self.traceMatrix[face, side]
            trace = dense(tensor, self.transpose(tensor.name()))
            folded = trace @ self.embeddingOfDegree(degree)
            scale = np.max(np.abs(folded))
            folded[np.abs(folded) <= self.TraceZeroTolerance * scale] = 0.0
            full = np.zeros(trace.shape)
            full[:, : folded.shape[1]] = folded
            stored = full.T if self.transpose(name) else full
            self._faceTracesOfDegree[name] = Tensor(
                name,
                stored.shape,
                spp={
                    index: repr(float(value))
                    for index, value in np.ndenumerate(stored)
                    if value != 0.0
                },
                **self.traceMatrixArgs(name),
            )
        return self._faceTracesOfDegree[name]

    def derivativeOfDegree(self, degree, direction):
        """The derivative in a direction from the basis of a degree into the one
        of one degree less (kDivMTDegree), negated as kDivMT is (negateFamily)
        and in the leading block of a matrix of the cell (embedFamily)."""
        return self.db.kDivMTDegree[degree, direction]

    def elevationOfDegree(self, degree):
        """The basis of a degree in the one of one degree more (elevateDegree),
        in the leading block of a matrix of the cell (embedFamily)."""
        return self.db.elevateDegree[degree]

    def leadingRows(self, rows, columns=None):
        """The pattern of a tensor of the cell, of `columns` quantities (those of
        Q by default), with entries in its leading rows."""
        columns = self.numQuantities() if columns is None else columns
        shape = self.Q.insertOptDim(
            (self.num3DBasisFunctions(), columns), (self.Q.optSize(),)
        )
        index = self.Q.insertOptDim((slice(0, rows), slice(None)), (slice(None),))
        pattern = np.zeros(shape, dtype=bool)
        pattern[index] = True
        return pattern

    def hornerTaylor(self, target, first, derivatives, scalars, name):
        """target = first (where there is one) + sum over k >= 1 of scalars[k]
        derivatives[k], with derivative k in the basis of degree order - 1 - k
        (derivativesByDegree): from the last derivative up, each partial sum
        raised by one degree before the next derivative is added (Horner's
        scheme), and the last one into the basis of the cell."""
        statements = []
        partial = None
        for k in range(len(derivatives) - 1, 0, -1):
            degree = self.order - 1 - k
            term = scalars[k] * derivatives[k]["kp"]
            if partial is not None:
                term = (
                    term
                    + self.elevationOfDegree(degree - 1)[self.t("kl")] * partial["lp"]
                )
            current = OptionalDimTensor(
                f"{name}({k})",
                self.Q.optName(),
                self.Q.optSize(),
                self.Q.optPos(),
                (self.num3DBasisFunctions(), self.numQuantities()),
                spp=self.leadingRows(count(degree)),
                alignStride=True,
                temporary=True,
            )
            current = self.definedOnce(current)
            statements.append(current["kp"] <= term)
            partial = current
        raised = self.elevationOfDegree(self.order - 2)[self.t("kl")] * partial["lp"]
        statements.append(target["kp"] <= (raised if first is None else first + raised))
        return statements

    def volumeDerivatives(self):
        """The derivative matrices the volume term of the corrector applies:
        the strong ones, the same the derivative chain applies, so that the
        volume term differentiates the field rather than the test function.

        The weak ones, which differentiate the test function, would return the
        derivative of the field less the lift of its trace on the four faces.
        With one operator per cell that lift cancels against the flux of the
        cell's own trace, and the two forms are the same scheme. With an
        operator that varies inside the cell it does not: the lift is multiplied
        by the operator at the points inside the cell, while the flux states the
        operator at the face, and a constant state would not stay constant, by
        an amount that does not shrink with the mesh -- for a material that
        varies smoothly inside the cell as much as for a cell whose metric
        varies.

        In the strong form the trace is not in the volume term, and the flux of
        every face subtracts the normal flux of the cell's own trace instead,
        where the operator is stated at the face (traceCorrector), a fault face
        included, whose flux is that subtraction alone. The host does that
        subtraction (initializer::toCorrectorForm). With one operator per cell
        the volume term in this form is what the chain forms anyway, a step of
        it applied to the time integral, or the time integral of the time
        derivative (see fusesVolume)."""
        return self.db.kDivMT

    def predictorAppliesVolume(self, target):
        """Whether the predictor of a target applies the volume term of the
        corrector itself, so that the target has no volume kernel: on the
        device, which keeps nothing per cell between the kernels of a batch
        but what the cells hold.

        The volume term is a step of the chain applied to the time integral
        (volumeDerivatives), through the operators the predictor has put
        together for the chain already (starAssembly), so the predictor applies
        it once the time integral is complete. It reads the field at the start
        of the step only, and nothing on the device reads the field after the
        predictor and before the corrector adds the fluxes of the faces: the
        receivers gather it before the predictor runs, the traces and the
        conditions of the boundary read the time integral, the displacement of
        a free surface the derivatives."""
        return target == "gpu"

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
        values = self.definedOnce(values)
        product = self.definedOnce(product)
        statements = [
            values[contract + spectator]
            <= self.operatorEval["nk"] * source["k" + contract[1:] + spectator]
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
        projected = self.operatorProject["kn"] * product["np" + spectator]
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
        # Only the cell integral is needed: the integrals of the basis functions,
        # a row of their own where the basis states them. A basis whose first
        # function is the constant one has them as the first row of M3, which
        # subselect narrows it to; it keeps the rank and sets the extent to 1,
        # which turns the kernel from an nb x nq product into an nq one.
        momentQ = OptionalDimTensor(
            "momentQ",
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            (1, self.numQuantities()),
            alignStride=True,
        )
        if "basisIntegrals" in self.db:
            integrals = self.db.basisIntegrals["Ij"]
        else:
            integrals = self.db.M3["Ij"].subselect("I", 0)
        generator.add("momentQCompute", momentQ["IJ"] <= integrals * self.Q["jJ"])

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
        # the reparametrisation of a shared face. The corrector reads it folded
        # into the trace and the lift of the side that writes its traces
        # reflected, and the setup of a face as a renumbering of the points of
        # the face, so no kernel names it; it is what both are checked against,
        # in every build, so it has to reach the generated code.
        for orientation in self.db.fP.values():
            include_tensors.add(orientation)
        # what a face stores where its rotation is per point, where those points
        # are, and where the operator is formed, in every build so that the host
        # can read them without asking which build it is
        include_tensors.add(self.TNodes)
        include_tensors.add(self.db.faceQuadPoints)
        include_tensors.add(self.db.faceQuadWeights)
        include_tensors.add(self.db.V2mToFaceQuad)
        include_tensors.add(self.operatorNodes)
        # and the metric, which a build that carries assembled star matrices
        # never reads but the host's setup of a cell still names
        for gradients in self.referenceGradients:
            include_tensors.add(gradients)
        # the weak derivative matrices, which no kernel applies; the strong form
        # of the volume term is checked against them (volumeDerivatives)
        for member in self.db.kDivM.values():
            include_tensors.add(member)
        # the strong ones, which the chain of a basis that is not hierarchical
        # takes through the lower degrees instead (modalDerivative,
        # derivativesByDegree); every chain is checked against them
        for member in self.db.kDivMT.values():
            include_tensors.add(member)
        # the flux of a face in the volume, which the corrector is checked
        # against in every build. No kernel names them: the corrector reads the
        # traces of the face and lifts what it forms with traceLift, or with
        # faceLift where the flux is scalars at the points of the face.
        for family in (self.db.rDivM, self.db.fMrT):
            for member in family.values():
                include_tensors.add(member)
        # the operators of the two sides of a face, which the host hands to the
        # kernels that take the side as data (_configureTraceSides); no kernel
        # names them
        for operator in getattr(self, "traceSideOperator", []):
            include_tensors.add(operator)
