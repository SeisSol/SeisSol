# SPDX-FileCopyrightText: 2016-2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff

import numpy as np
from kernels.aderdg import ADERDGBase
from kernels.multsim import OptionalDimTensor
from yateto import Scalar, Tensor, simpleParameterSpace
from yateto.ast.node import Add
from yateto.input import (memoryLayoutFromFile, parseJSONMatrixFile,
                          parseXMLMatrixFile)
from yateto.memory import CSCMemoryLayout


class Viscoelastic2ADERDG(ADERDGBase):
    def __init__(
        self,
        order,
        multipleSimulations,
        matricesDir,
        memLayout,
        numberOfMechanisms,
        **kwargs,
    ):
        super().__init__(order, multipleSimulations, matricesDir)

        self.numberOfMechanisms = numberOfMechanisms

        clones = {
            "star": ["star(0)", "star(1)", "star(2)"],
        }
        self.db.update(
            parseXMLMatrixFile(
                "{}/matrices_viscoelastic.xml".format(matricesDir), clones
            )
        )
        memoryLayoutFromFile(memLayout, self.db, clones)

        self._qShapeExtended = (
            self.numberOf3DBasisFunctions(),
            self.numberOfExtendedQuantities(),
        )
        self._qShapeAnelastic = (
            self.numberOf3DBasisFunctions(),
            self.numberOfAnelasticQuantities(),
            self.numberOfMechanisms,
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
                self.numberOfAnelasticQuantities(),
                self.numberOfMechanisms,
                self.numberOfQuantities(),
            ),
        )
        self.w = Tensor("w", (self.numberOfMechanisms,))
        self.W = Tensor(
            "W",
            (self.numberOfMechanisms, self.numberOfMechanisms),
            np.eye(self.numberOfMechanisms, dtype=bool),
            CSCMemoryLayout,
        )

        selectElaSpp = np.zeros(
            (self.numberOfExtendedQuantities(), self.numberOfQuantities())
        )
        selectElaSpp[0 : self.numberOfQuantities(), 0 : self.numberOfQuantities()] = (
            np.eye(self.numberOfQuantities())
        )
        self.selectEla = Tensor(
            "selectEla",
            (self.numberOfExtendedQuantities(), self.numberOfQuantities()),
            selectElaSpp,
            CSCMemoryLayout,
        )

        selectAneSpp = np.zeros(
            (
                self.numberOfExtendedQuantities(),
                self.numberOfAnelasticQuantities(),
            )
        )
        selectAneSpp[
            self.numberOfQuantities() : self.numberOfExtendedQuantities(),
            0 : self.numberOfAnelasticQuantities(),
        ] = np.eye(self.numberOfAnelasticQuantities())
        self.selectAne = Tensor(
            "selectAne",
            (
                self.numberOfExtendedQuantities(),
                self.numberOfAnelasticQuantities(),
            ),
            selectAneSpp,
            CSCMemoryLayout,
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

    def numberOfQuantities(self):
        return 9

    def numberOfAnelasticQuantities(self):
        return 6

    def numberOfExtendedQuantities(self):
        """Number of quantities for fused computation of elastic and anelastic update."""
        return self.numberOfQuantities() + self.numberOfAnelasticQuantities()

    def numberOfFullQuantities(self):
        """Number of quantities when unrolling anelastic tensor into a matrix."""
        return (
            self.numberOfQuantities()
            + self.numberOfMechanisms * self.numberOfAnelasticQuantities()
        )

    def extendedQTensor(self):
        return self.Qext

    def starMatrix(self, dim):
        return self.db.star[dim]

    def name(self):
        return "viscoelastic2"

    def addInit(self, generator):
        super().addInit(generator)

        selectElaFullSpp = np.zeros(
            (self.numberOfFullQuantities(), self.numberOfQuantities())
        )
        selectElaFullSpp[
            0 : self.numberOfQuantities(), 0 : self.numberOfQuantities()
        ] = np.eye(self.numberOfQuantities())
        selectElaFull = Tensor(
            "selectElaFull",
            (self.numberOfFullQuantities(), self.numberOfQuantities()),
            selectElaFullSpp,
            CSCMemoryLayout,
        )

        selectAneFullSpp = np.zeros(
            (
                self.numberOfFullQuantities(),
                self.numberOfAnelasticQuantities(),
                self.numberOfMechanisms,
            )
        )
        for mech in range(self.numberOfMechanisms):
            q1 = self.numberOfQuantities() + mech * self.numberOfAnelasticQuantities()
            q2 = q1 + self.numberOfAnelasticQuantities()
            selectAneFullSpp[q1:q2, :, mech] = np.eye(
                self.numberOfAnelasticQuantities()
            )
        selectAneFull = Tensor(
            "selectAneFull",
            (
                self.numberOfFullQuantities(),
                self.numberOfAnelasticQuantities(),
                self.numberOfMechanisms,
            ),
            selectAneFullSpp,
        )

        iniShape = (
            self.numberOf3DQuadraturePoints(),
            self.numberOfFullQuantities(),
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
            self.numberOf3DQuadraturePoints(),
            self.numberOfQuantities(),
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

    def addLocal(self, generator, targets):
        volumeSum = Add()
        for i in range(3):
            volumeSum += (
                self.db.kDivM[i][self.t("kl")] * self.I["lq"] * self.db.star[i]["qp"]
            )
        volumeExt = self.Qext["kp"] <= volumeSum
        generator.add("volumeExt", volumeExt)

        localFluxExt = (
            lambda i: self.Qext["kp"]
            <= self.Qext["kp"]
            + self.db.rDivM[i][self.t("km")]
            * self.db.fMrT[i][self.t("ml")]
            * self.I["lq"]
            * self.AplusT["qp"]
        )
        localFluxExtPrefetch = lambda i: (
            self.I if i == 0 else (self.Q if i == 1 else None)
        )
        generator.addFamily(
            "localFluxExt",
            simpleParameterSpace(4),
            localFluxExt,
            localFluxExtPrefetch,
        )

        generator.add(
            "local",
            [
                self.Qane["kpm"]
                <= self.Qane["kpm"]
                + self.w["m"] * self.Qext["kq"] * self.selectAne["qp"]
                + self.Iane["kpl"] * self.W["lm"],
                self.Q["kp"]
                <= self.Q["kp"]
                + self.Qext["kq"] * self.selectEla["qp"]
                + self.Iane["kqm"] * self.E["qmp"],
            ],
        )

    def addNeighbor(self, generator, targets):
        neighbourFluxExt = (
            lambda h, j, i: self.Qext["kp"]
            <= self.Qext["kp"]
            + self.db.rDivM[i][self.t("km")]
            * self.db.fP[h][self.t("mn")]
            * self.db.rT[j][self.t("nl")]
            * self.I["lq"]
            * self.AminusT["qp"]
        )
        neighbourFluxExtPrefetch = lambda h, j, i: self.I
        generator.addFamily(
            "neighbourFluxExt",
            simpleParameterSpace(3, 4, 4),
            neighbourFluxExt,
            neighbourFluxExtPrefetch,
        )

        generator.add(
            "neighbour",
            [
                self.Qane["kpm"]
                <= self.Qane["kpm"]
                + self.w["m"] * self.Qext["kq"] * self.selectAne["qp"],
                self.Q["kp"] <= self.Q["kp"] + self.Qext["kq"] * self.selectEla["qp"],
            ],
        )

    def addTime(self, generator, targets):
        qShape = (self.numberOf3DBasisFunctions(), self.numberOfQuantities())
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

        derivativeTaylorExpansionEla = Add()
        # derivativeTaylorExpansionAne = Add()
        for d in range(0, self.order):
            derivativeTaylorExpansionEla += powers[d] * dQ[d]["kp"]
            # derivativeTaylorExpansionAne += powers[d] * dQane[d]['kpm']
        derivativeTaylorExpansionElaExpr = self.I["kp"] <= derivativeTaylorExpansionEla
        # derivativeTaylorExpansionAneExpr = self.Iane['kpm'] <= derivativeTaylorExpansionAne

        def derivative(kthDer):
            derivativeSum = Add()
            for j in range(3):
                derivativeSum += (
                    self.db.kDivMT[j][self.t("kl")]
                    * dQ[kthDer - 1]["lq"]
                    * self.db.star[j]["qp"]
                )
            return derivativeSum

        # WARNING: the following kernel may produce incorrect results,
        # if not executed in the order as specified here
        # the reason for that is that dQext, dQane (except dQane(0))
        # and potentially dQ (except dQ(0)) are allocated in temporary arrays
        # which are smaller than the whole tensor families
        # (even indices share the same buffer,
        # and odd indices share the same buffer)
        derivativeExpr = [
            self.I["kp"] <= powers[0] * dQ[0]["kp"],
            self.Iane["kpm"] <= powers[0] * dQane[0]["kpm"],
        ]
        for d in range(1, self.order):
            derivativeExpr += [
                dQext[d]["kp"] <= derivative(d),
                dQ[d]["kp"]
                <= dQext[d]["kq"] * self.selectEla["qp"]
                + dQane[d - 1]["kqm"] * self.E["qmp"],
                dQane[d]["kpm"]
                <= self.w["m"] * dQext[d]["kq"] * self.selectAne["qp"]
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

        generator.add("derivative", derivativeExpr)
        generator.add("derivativeTaylorExpansionEla", derivativeTaylorExpansionElaExpr)

    def add_include_tensors(self, include_tensors):
        super().add_include_tensors(include_tensors)
        include_tensors.add(self.db.nodes2D)
        # Nodal flux kernel uses this matrix but is not supported by visco2
        include_tensors.update([self.db.project2nFaceTo3m[i] for i in range(4)])


EQUATION_CLASS = Viscoelastic2ADERDG
