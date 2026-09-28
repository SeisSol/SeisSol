# SPDX-FileCopyrightText: 2016 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff

import numpy as np
from kernels.common import generate_kernel_name_prefix
from kernels.multsim import OptionalDimTensor
from kernels.quantities import layout, total_extent
from yateto import Scalar, Tensor, ops, simpleParameterSpace
from yateto.ast.node import Accumulate
from yateto.input import parseJSONMatrixFile
from yateto.memory import CSCMemoryLayout
from yateto.util import tensor_collection_from_constant_expression

from .aderdg import ADERDGBase


class LinearCKAnelastic(ADERDGBase):
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

        self.numMechanisms = numMechanisms

        self._qShapeExtended = (
            self.num3DBasisFunctions(),
            self.numExtendedQuantities(),
        )
        self._qShapeAnelastic = (
            self.num3DBasisFunctions(),
            self.numAnelasticQuantities(),
            self.numMechanisms,
        )
        self.Qext = OptionalDimTensor(
            "Qext",
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            self._qShapeExtended,
            alignStride=True,
        )
        self.Qane = OptionalDimTensor(
            "Qane",
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            self._qShapeAnelastic,
            alignStride=True,
        )
        self.Iane = OptionalDimTensor(
            "Iane",
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            self._qShapeAnelastic,
            alignStride=True,
        )

        self.E = Tensor(
            "E",
            (
                self.numAnelasticQuantities(),
                self.numMechanisms,
                self.numQuantities(),
            ),
        )
        self.w = Tensor("w", (self.numMechanisms,))

        self.W = Tensor(
            "W",
            (self.numMechanisms, self.numMechanisms),
            spp=np.eye(self.numMechanisms, dtype=bool),
            memoryLayoutClass=CSCMemoryLayout,
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

        self.kwargs = kwargs

    def extendedBlocks(self):
        """The face rotation sees one anelastic block, not one per mechanism:
        this solver keeps the mechanism index as a separate tensor dimension."""
        return layout(self.primaryGroups(), self.mechanismGroups(), 1)

    def inverseRotationBlocks(self):
        return layout(self.primaryGroups())

    def numAnelasticQuantities(self):
        return total_extent(layout(self.mechanismGroups()))

    def anelasticSource(self, source, target):
        """The relaxation added to a target, in whichever shape this build
        forms it: from the tensor a cell carries, or from the scalars it
        carries at the sample points."""
        if self.sourceCoefficientCount() == 0:
            return [target["kp"] <= target["kp"] + source["kqm"] * self.E["qmp"]]
        return self.nodalSource(source, target, "nqm")

    def sourceStructurePrototype(self):
        """The relaxation, which this solver keeps in a tensor of its own with
        the mechanism as a dimension."""
        return self.E

    def numExtendedQuantities(self):
        """Return the number of quantities for fused computation of elastic and anelastic update."""
        return self.numQuantities() + self.numAnelasticQuantities()

    def numFullQuantities(self):
        """Return the number of quantities when unrolling anelastic tensor into a matrix."""
        return self.numQuantities() + self.numMechanisms * self.numAnelasticQuantities()

    def extendedQTensor(self):
        return self.Qext

    def name(self):
        return "linearckanelastic"

    def addInit(self, generator):
        super().addInit(generator)

        selectElaFullSpp = np.zeros((self.numFullQuantities(), self.numQuantities()))
        selectElaFullSpp[0 : self.numQuantities(), 0 : self.numQuantities()] = np.eye(
            self.numQuantities()
        )
        selectElaFull = Tensor(
            "selectElaFull",
            (self.numFullQuantities(), self.numQuantities()),
            selectElaFullSpp,
            CSCMemoryLayout,
        )

        selectAneFullSpp = np.zeros(
            (
                self.numFullQuantities(),
                self.numAnelasticQuantities(),
                self.numMechanisms,
            )
        )
        for mech in range(self.numMechanisms):
            q1 = self.numQuantities() + mech * self.numAnelasticQuantities()
            q2 = q1 + self.numAnelasticQuantities()
            selectAneFullSpp[q1:q2, :, mech] = np.eye(self.numAnelasticQuantities())
        selectAneFull = Tensor(
            "selectAneFull",
            (
                self.numFullQuantities(),
                self.numAnelasticQuantities(),
                self.numMechanisms,
            ),
            selectAneFullSpp,
        )

        iniShape = (
            self.num3DQuadraturePoints(),
            self.numFullQuantities(),
        )
        iniCond = OptionalDimTensor(
            "iniCond",
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            iniShape,
            alignStride=True,
        )
        dofsShape = (
            self.num3DQuadraturePoints(),
            self.numQuantities(),
        )
        dofsQP = OptionalDimTensor(
            "dofsQP",
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            dofsShape,
            alignStride=True,
        )

        projectIniCondEla = (
            self.Q["kp"]
            <= self.db.projectQP[self.t("kl")] * iniCond["lq"] * selectElaFull["qp"]
        )
        projectIniCondAne = (
            self.Qane["kpm"]
            <= self.db.projectQP[self.t("kl")] * iniCond["lq"] * selectAneFull["qpm"]
        )
        generator.add("projectIniCond", [projectIniCondEla, projectIniCondAne])
        generator.add(
            "evalAtQP",
            dofsQP["kp"] <= self.db.evalAtQP[self.t("kl")] * self.Q["lp"],
        )

        self.addAnelasticEnergyProducts(generator)

    def addAnelasticEnergyProducts(self, generator):
        """Mass-matrix moments involving the anelastic variables.

        The viscoelastic energies are quadratic forms in the combined state
        z = (sigma, vartheta^(1..L)). ADERDGBase.addEnergyProducts already covers
        the sigma-sigma block via momentQQCompute; these two add the remaining ones:

          momentQaneQane[I,J,m,n]  == \\int_{T_ref} Qane_{I,m} Qane_{J,n}
          momentQQane[I,J,m]  == \\int_{T_ref} Q_I Qane_{J,m}

        Cross-mechanism terms (m != n) are needed: the relaxed strain is
        eps = C_r^-1 (sigma - sum_l D^(l) Qane^(l) / omega_l), so squaring it
        couples every pair of mechanisms.

        Both are quadratic in an OptionalDimTensor, i.e. they rely on yateto
        handling the fused-simulation index as a batch index.
        """
        momentQaneQane = OptionalDimTensor(
            "momentQaneQane",
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            (
                self.numAnelasticQuantities(),
                self.numAnelasticQuantities(),
                self.numMechanisms,
                self.numMechanisms,
            ),
        )
        generator.add(
            "momentQaneQaneCompute",
            momentQaneQane["IJmn"]
            <= self.db.M3["ij"] * self.Qane["iIm"] * self.Qane["jJn"],
        )

        momentQQane = OptionalDimTensor(
            "momentQQane",
            self.Q.optName(),
            self.Q.optSize(),
            self.Q.optPos(),
            (
                self.numQuantities(),
                self.numAnelasticQuantities(),
                self.numMechanisms,
            ),
        )
        generator.add(
            "momentQQaneCompute",
            momentQQane["IJm"] <= self.db.M3["ij"] * self.Q["iI"] * self.Qane["jJm"],
        )

    def addLocal(self, generator, targets):
        for target in targets:
            name_prefix = generate_kernel_name_prefix(target)

            if self.nodalMaterial:
                volumeExpr = self.nodalApply(self.I, self.Qext, self.db.kDivM)
            else:
                volumeSum = Accumulate(ops.Add())
                for i in range(3):
                    volumeSum += (
                        self.db.kDivM[i][self.t("kl")]
                        * self.I["lq"]
                        * self.starMatrix(i)["qp"]
                    )
                volumeExpr = [self.Qext["kp"] <= volumeSum]
            generator.add(
                f"{name_prefix}volumeExt",
                self.starAssembly() + volumeExpr,
                target=target,
            )

            plusFluxMatrixAccessor = (
                lambda i: self.db.rDivM[i][self.t("km")] * self.db.fMrT[i][self.t("ml")]
            )
            if self.kwargs["enable_premultiply_flux"] and target == "gpu":
                contractionResult = tensor_collection_from_constant_expression(
                    "plusFluxMatrices",
                    plusFluxMatrixAccessor,
                    simpleParameterSpace(4),
                    target_indices="kl",
                )
                self.db.update(contractionResult)
                plusFluxMatrixAccessor = lambda i: self.db.plusFluxMatrices[i]["kl"]

            if self.nodalFaceFlux:
                localFluxExt = lambda i: self.nodalFlux(
                    self.I,
                    self.Qext,
                    self.db.V3mTo2nFace[i][self.t("nl")],
                    self.db.project2nFaceTo3m[i]["kn"],
                    self.fluxCoefficientsLocal,
                )
            else:
                localFluxExt = (
                    lambda i: self.Qext["kp"]
                    <= self.Qext["kp"]
                    + plusFluxMatrixAccessor(i) * self.I["lq"] * self.AplusT["qp"]
                )
            localFluxExtPrefetch = lambda i: (
                self.I if i == 0 else (self.Q if i == 1 else None)
            )
            generator.addFamily(
                f"{name_prefix}localFluxExt",
                simpleParameterSpace(4),
                localFluxExt,
                localFluxExtPrefetch,
                target=target,
            )

            local_ops = [
                self.Qane["kpm"]
                <= self.Qane["kpm"]
                + self.w["m"]
                * self.Qext["kp"].subslice(
                    "p",
                    self.numQuantities(),
                    self.numExtendedQuantities(),
                )
                + self.Iane["kpl"] * self.W["lm"],
                self.Q["kp"]
                <= self.Q["kp"]
                + self.Qext["kp"].subslice("p", 0, self.numQuantities()),
            ] + self.anelasticSource(self.Iane, self.Q)
            generator.add(
                f"{name_prefix}local",
                local_ops,
                target=target,
            )

            if self.nodalFaceFlux:
                flux_ops = self.nodalLocalFluxAll(self.I, self.Qext)
            else:
                flux_ops = [
                    self.Qext["kp"]
                    <= sum(
                        [
                            plusFluxMatrixAccessor(i)
                            * self.I["lq"]
                            * self.AplusTAll[i]["qp"]
                            for i in range(4)
                        ],
                        start=self.Qext["kp"],
                    )
                ]
            generator.add(
                f"{name_prefix}fluxLocalAll",
                flux_ops + local_ops,
                target=target,
            )

    def addNeighbor(self, generator, targets):
        for target in targets:
            name_prefix = generate_kernel_name_prefix(target)

            minusFluxMatrixAccessor = (
                lambda j, i: self.db.rDivM[i][self.t("km")]
                * self.db.fPrT[j][self.t("ml")]
            )
            if self.kwargs["enable_premultiply_flux"] and target == "gpu":
                contractionResult = tensor_collection_from_constant_expression(
                    "minusFluxMatrices",
                    minusFluxMatrixAccessor,
                    simpleParameterSpace(4, 4),
                    target_indices="kl",
                )
                self.db.update(contractionResult)
                minusFluxMatrixAccessor = lambda j, i: self.db.minusFluxMatrices[j, i][
                    "kl"
                ]

            if self.nodalFaceFlux:
                # every regular face has the face orientation index zero, see
                # LinearCK.addNeighbor
                neighborFluxExt = lambda j, i: self.nodalFlux(
                    self.I,
                    self.Qext,
                    self.db.neighborToFace[0, j]["nl"],
                    self.db.project2nFaceTo3m[i]["kn"],
                    self.fluxCoefficientsNeighbor,
                )
            else:
                neighborFluxExt = (
                    lambda j, i: self.Qext["kp"]
                    <= self.Qext["kp"]
                    + minusFluxMatrixAccessor(j, i) * self.I["lq"] * self.AminusT["qp"]
                )
            neighborFluxExtPrefetch = lambda j, i: self.I
            generator.addFamily(
                f"{name_prefix}neighborFluxExt",
                simpleParameterSpace(4, 4),
                neighborFluxExt,
                neighborFluxExtPrefetch,
                target=target,
            )

            generator.add(
                f"{name_prefix}neighbor",
                [
                    self.Qane["kpm"]
                    <= self.Qane["kpm"]
                    + self.w["m"]
                    * self.Qext["kp"].subslice(
                        "p",
                        self.numQuantities(),
                        self.numExtendedQuantities(),
                    ),
                    self.Q["kp"]
                    <= self.Q["kp"]
                    + self.Qext["kp"].subslice("p", 0, self.numQuantities()),
                ],
                target=target,
            )

    def addTime(self, generator, targets):
        qShape = (self.num3DBasisFunctions(), self.numQuantities())
        dQ = [
            OptionalDimTensor(
                "dQ({})".format(d),
                self.Q.optName(),
                self.Q.optSize(),
                self.Q.optPos(),
                qShape,
                alignStride=True,
            )
            for d in range(self.order)
        ]
        self.dQs = dQ
        dQext = [
            OptionalDimTensor(
                "dQext({})".format(d),
                self.Q.optName(),
                self.Q.optSize(),
                self.Q.optPos(),
                self._qShapeExtended,
                alignStride=True,
            )
            for d in range(self.order)
        ]
        dQane = [
            OptionalDimTensor(
                "dQane({})".format(d),
                self.Q.optName(),
                self.Q.optSize(),
                self.Q.optPos(),
                self._qShapeAnelastic,
                alignStride=True,
            )
            for d in range(self.order)
        ]

        powers = [Scalar(f"power({i})") for i in range(self.order)]

        for target in targets:
            name_prefix = generate_kernel_name_prefix(target)

            derivativeTaylorExpansionEla = Accumulate(ops.Add())
            # derivativeTaylorExpansionAne = Accumulate(ops.Add())
            for d in range(0, self.order):
                derivativeTaylorExpansionEla += powers[d] * dQ[d]["kp"]
                # derivativeTaylorExpansionAne += powers[d] * dQane[d]['kpm']
            derivativeTaylorExpansionElaExpr = (
                self.I["kp"] <= derivativeTaylorExpansionEla
            )
            # derivativeTaylorExpansionAneExpr = self.Iane['kpm'] <= derivativeTaylorExpansionAne

            def derivative(kthDer):
                derivativeSum = Accumulate(ops.Add())
                for j in range(3):
                    derivativeSum += (
                        self.db.kDivMT[j][self.t("kl")]
                        * dQ[kthDer - 1]["lq"]
                        * self.starMatrix(j)["qp"]
                    )
                return derivativeSum

            def derivativeStep(kthDer):
                """One step of the chain, in whichever shape the operator has.

                A constant operator drops a degree with every derivative and the
                generator narrows the matrices accordingly; one that varies
                inside the cell is read where its samples are, so nothing
                narrows and every derivative stays full.
                """
                if self.nodalMaterial:
                    return self.nodalApply(
                        dQ[kthDer - 1], dQext[kthDer], self.db.kDivMT
                    )
                return [dQext[kthDer]["kp"] <= derivative(kthDer)]

            # WARNING: the following kernel may produce incorrect results,
            # if not executed in the order as specified here
            # the reason for that is that dQext, dQane (except dQane(0))
            # and potentially dQ (except dQ(0)) are allocated in temporary arrays
            # which are smaller than the whole tensor families
            # (even indices share the same buffer,
            # and odd indices share the same buffer)

            if target == "gpu":
                derivativeExpr = self.starAssembly() + [
                    dQ[0]["kp"] <= self.Q["kp"],
                    self.I["kp"] <= powers[0] * self.Q["kp"],  # == dQ[0]
                ]
            else:
                derivativeExpr = self.starAssembly() + [
                    self.I["kp"] <= powers[0] * dQ[0]["kp"],
                ]

            derivativeExpr += [
                self.Iane["kpm"] <= powers[0] * dQane[0]["kpm"],
            ]

            for d in range(1, self.order):
                derivativeExpr += derivativeStep(d)
                derivativeExpr += [
                    dQ[d]["kp"]
                    <= dQext[d]["kp"].subslice("p", 0, self.numQuantities()),
                ]
                derivativeExpr += self.anelasticSource(dQane[d - 1], dQ[d])
                derivativeExpr += [
                    dQane[d]["kpm"]
                    <= self.w["m"]
                    * dQext[d]["kp"].subslice(
                        "p",
                        self.numQuantities(),
                        self.numExtendedQuantities(),
                    )
                    + dQane[d - 1]["kpl"] * self.W["lm"],
                    self.I["kp"] <= self.I["kp"] + powers[d] * dQ[d]["kp"],
                    self.Iane["kpm"] <= self.Iane["kpm"] + powers[d] * dQane[d]["kpm"],
                ]
            # TODO(David): we'll need to add intermediate results to Yateto,
            # then the temporary storage needed can be reduced.
            # for now, we'll interleave the Taylor
            # expansion with the derivative computation
            # derivativeExpr += [
            #   derivativeTaylorExpansionElaExpr,
            #   derivativeTaylorExpansionAneExpr
            # ]

            generator.add(f"{name_prefix}derivative", derivativeExpr, target=target)
            generator.add(
                f"{name_prefix}derivativeTaylorExpansionEla",
                derivativeTaylorExpansionElaExpr,
                target=target,
            )

    def add_include_tensors(self, include_tensors):
        super().add_include_tensors(include_tensors)
        include_tensors.add(self.db.nodes2D)
        # the relaxation a cell carries, which no kernel names where the
        # material varies inside the cell and the term is formed at its samples
        include_tensors.add(self.E)
        # Nodal flux kernel uses this matrix but is not supported by visco2
        include_tensors.update([self.db.project2nFaceTo3m[i] for i in range(4)])
