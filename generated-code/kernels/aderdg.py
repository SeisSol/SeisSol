# SPDX-FileCopyrightText: 2019-2024 SeisSol Group
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
from yateto import Scalar, Tensor, simpleParameterSpace
from yateto.ast.node import Add
from yateto.ast.transformer import DeduceIndices, EquivalentSparsityPattern
from yateto.input import parseJSONMatrixFile, parseXMLMatrixFile
from yateto.memory import CSCMemoryLayout
from yateto.util import (tensor_collection_from_constant_expression,
                         tensor_from_constant_expression)


class ADERDGBase(ABC):
    def __init__(self, order, multipleSimulations, matricesDir):
        self.order = order

        self.alignStride = lambda name: True
        if multipleSimulations > 1:
            self.alignStride = lambda name: name.startswith("fP")
        transpose = multipleSimulations > 1
        self.transpose = lambda name: transpose
        self.t = (lambda x: x[::-1]) if transpose else (lambda x: x)

        self.db = parseXMLMatrixFile(
            "{}/matrices_{}.xml".format(matricesDir, self.numberOf3DBasisFunctions()),
            transpose=self.transpose,
            alignStride=self.alignStride,
        )
        clonesQP = {"v": ["evalAtQP"], "vInv": ["projectQP"]}
        self.db.update(
            parseXMLMatrixFile(
                "{}/plasticity_ip_matrices_{}.xml".format(matricesDir, order),
                clonesQP,
                transpose=self.transpose,
                alignStride=self.alignStride,
            )
        )
        self.db.update(
            parseJSONMatrixFile("{}/sampling_directions.json".format(matricesDir))
        )
        self.db.update(
            parseJSONMatrixFile("{}/mass_{}.json".format(matricesDir, order))
        )

        qShape = (self.numberOf3DBasisFunctions(), self.numberOfQuantities())
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
            False,  # multipleSimulations,
            0,
            (self.numberOf2DBasisFunctions(), self.numberOfQuantities()),
            alignStride=True,
        )

        project2nFaceTo3m = tensor_collection_from_constant_expression(
            base_name="project2nFaceTo3m",
            expressions=lambda i: self.db.rDivM[i]["jk"] * self.db.V2nTo2m["kl"],
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

    def name(self):
        return ""

    def numberOf2DBasisFunctions(self):
        return self.order * (self.order + 1) // 2

    def numberOf3DBasisFunctions(self):
        return self.order * (self.order + 1) * (self.order + 2) // 6

    def numberOf3DQuadraturePoints(self):
        return (self.order + 1) ** 3

    def godunov_spp(self):
        shape = (self.numberOfQuantities(), self.numberOfQuantities())
        return np.ones(shape, dtype=bool)

    def flux_solver_spp(self):
        shape = (self.numberOfQuantities(), self.numberOfExtendedQuantities())
        return np.ones(shape, dtype=bool)

    def transformation_spp(self):
        shape = (
            self.numberOfExtendedQuantities(),
            self.numberOfExtendedQuantities(),
        )
        return np.ones(shape, dtype=bool)

    def transformation_inv_spp(self):
        return self.godunov_spp()

    def extractVelocities(self):
        extractVelocitiesSPP = np.zeros((3, self.numberOfQuantities()))
        extractVelocitiesSPP[0, 6] = 1
        extractVelocitiesSPP[1, 7] = 1
        extractVelocitiesSPP[2, 8] = 1
        return extractVelocitiesSPP

    def mapToVelocities(self):
        return self.extractVelocities().T

    def extractTractions(self):
        extractTractionsSPP = np.zeros((3, self.numberOfQuantities()))
        extractTractionsSPP[0, 0] = 1
        extractTractionsSPP[1, 3] = 1
        extractTractionsSPP[2, 5] = 1
        return extractTractionsSPP

    def mapToTractions(self):
        return self.extractTractions().T

    @abstractmethod
    def numberOfQuantities(self):
        pass

    @abstractmethod
    def numberOfExtendedQuantities(self):
        pass

    @abstractmethod
    def extendedQTensor(self):
        pass

    @abstractmethod
    def starMatrix(self, dim):
        pass

    def addInit(self, generator):
        flux_solver_spp = self.flux_solver_spp()
        self.QcorrLocal = Tensor("QcorrLocal", flux_solver_spp.shape)
        self.QcorrNeighbor = Tensor("QcorrNeighbor", flux_solver_spp.shape)

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

        QFortran = Tensor(
            "QFortran",
            (self.numberOf3DBasisFunctions(), self.numberOfQuantities()),
        )
        multSimToFirstSim = Tensor(
            "multSimToFirstSim", (self.Q.optSize(),), spp={(0,): "1.0"}
        )
        if self.Q.hasOptDim():
            copyQToQFortran = QFortran["kp"] <= self.Q["kp"] * multSimToFirstSim["s"]
        else:
            copyQToQFortran = QFortran["kp"] <= self.Q["kp"]

        generator.add("copyQToQFortran", copyQToQFortran)

        stiffnessTensor = Tensor("stiffnessTensor", (3, 3, 3, 3))
        direction = Tensor("direction", (3,))
        christoffel = Tensor("christoffel", (3, 3))

        computeChristoffel = (
            christoffel["ik"]
            <= stiffnessTensor["ijkl"] * direction["j"] * direction["l"]
        )
        generator.add("computeChristoffel", computeChristoffel)

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


class LinearADERDG(ADERDGBase):
    def name(self):
        return "linear"

    def sourceMatrix(self):
        return None

    def extendedQTensor(self):
        return self.Q

    def numberOfExtendedQuantities(self):
        return self.numberOfQuantities()

    def addInit(self, generator):
        super().addInit(generator)

        iniShape = (
            self.numberOf3DQuadraturePoints(),
            self.numberOfQuantities(),
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

    def addLocal(self, generator, targets):
        for target in targets:
            name_prefix = generate_kernel_name_prefix(target)
            volumeSum = self.Q["kp"]
            for i in range(3):
                volumeSum += (
                    self.db.kDivM[i][self.t("kl")]
                    * self.I["lq"]
                    * self.starMatrix(i)["qp"]
                )
            if self.sourceMatrix():
                volumeSum += self.I["kq"] * self.sourceMatrix()["qp"]
            volume = self.Q["kp"] <= volumeSum
            generator.add(f"{name_prefix}volume", volume, target=target)

            localFluxNodal = (
                lambda i: self.Q["kp"]
                <= self.Q["kp"]
                + self.db.project2nFaceTo3m[i]["kn"]
                * self.INodal["no"]
                * self.AminusT["op"]
            )
            localFluxNodalPrefetch = lambda i: (
                self.I if i == 0 else (self.Q if i == 1 else None)
            )
            generator.addFamily(
                f"{name_prefix}localFluxNodal",
                simpleParameterSpace(4),
                localFluxNodal,
                localFluxNodalPrefetch,
                target=target,
            )

        localFlux = (
            lambda i: self.Q["kp"]
            <= self.Q["kp"]
            + self.db.rDivM[i][self.t("km")]
            * self.db.fMrT[i][self.t("ml")]
            * self.I["lq"]
            * self.AplusT["qp"]
        )
        localFluxPrefetch = lambda i: (
            self.I if i == 0 else (self.Q if i == 1 else None)
        )
        generator.addFamily(
            "localFlux",
            simpleParameterSpace(4),
            localFlux,
            localFluxPrefetch,
            target="cpu",
        )

        if "gpu" in targets:
            plusFluxMatrixAccessor = (
                lambda i: self.db.rDivM[i][self.t("km")] * self.db.fMrT[i][self.t("ml")]
            )
            if self.kwargs["enable_premultiply_flux"]:
                contractionResult = tensor_collection_from_constant_expression(
                    "plusFluxMatrices",
                    plusFluxMatrixAccessor,
                    simpleParameterSpace(4),
                    target_indices="kl",
                )
                self.db.update(contractionResult)
                plusFluxMatrixAccessor = lambda i: self.db.plusFluxMatrices[i]["kl"]

            localFlux = (
                lambda i: self.Q["kp"]
                <= self.Q["kp"]
                + plusFluxMatrixAccessor(i) * self.I["lq"] * self.AplusT["qp"]
            )
            generator.addFamily(
                "gpu_localFlux",
                simpleParameterSpace(4),
                localFlux,
                target="gpu",
            )

    def addNeighbor(self, generator, targets):
        neighborFlux = (
            lambda h, j, i: self.Q["kp"]
            <= self.Q["kp"]
            + self.db.rDivM[i][self.t("km")]
            * self.db.fP[h][self.t("mn")]
            * self.db.rT[j][self.t("nl")]
            * self.I["lq"]
            * self.AminusT["qp"]
        )
        neighborFluxPrefetch = lambda h, j, i: self.I
        generator.addFamily(
            "neighboringFlux",
            simpleParameterSpace(3, 4, 4),
            neighborFlux,
            neighborFluxPrefetch,
            target="cpu",
        )

        if "gpu" in targets:
            minusFluxMatrixAccessor = (
                lambda h, j, i: self.db.rDivM[i][self.t("km")]
                * self.db.fP[h][self.t("mn")]
                * self.db.rT[j][self.t("nl")]
            )
            if self.kwargs["enable_premultiply_flux"]:
                contractionResult = tensor_collection_from_constant_expression(
                    "minusFluxMatrices",
                    minusFluxMatrixAccessor,
                    simpleParameterSpace(3, 4, 4),
                    target_indices="kl",
                )
                self.db.update(contractionResult)
                minusFluxMatrixAccessor = lambda h, j, i: self.db.minusFluxMatrices[
                    h, j, i
                ]["kl"]

            neighborFlux = (
                lambda h, j, i: self.Q["kp"]
                <= self.Q["kp"]
                + minusFluxMatrixAccessor(h, j, i) * self.I["lq"] * self.AminusT["qp"]
            )
            generator.addFamily(
                "gpu_neighboringFlux",
                simpleParameterSpace(3, 4, 4),
                neighborFlux,
                target="gpu",
            )

    def addTime(self, generator, targets):
        powers = [Scalar(f"power({i})") for i in range(self.order)]
        for target in targets:
            name_prefix = generate_kernel_name_prefix(target)

            qShape = (
                self.numberOf3DBasisFunctions(),
                self.numberOfQuantities(),
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
            derivatives = [dQ0]

            # for now, interleave Taylor expansion and derivative computation
            derivativeExpr = [self.I["kp"] <= power * dQ0["kp"]]
            derivativeTaylorExpansion = power * dQ0["kp"]

            self.dQs = [dQ0]

            for i in range(1, self.order):
                power = powers[i]
                derivativeSum = Add()
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

                # for now, we interleave derivative
                # and derivativeTaylorExpansion kernels.
                derivativeExpr += [
                    dQ["kp"] <= derivativeSum,
                    self.I["kp"] <= self.I["kp"] + power * dQ["kp"],
                ]
                derivativeTaylorExpansion += power * dQ["kp"]

                derivatives.append(dQ)

            derivativeTaylorExpansionExpr = self.I["kp"] <= derivativeTaylorExpansion
            # derivativeExpr += [derivativeTaylorExpansionExpr]
            generator.add(f"{name_prefix}derivative", derivativeExpr, target=target)
            generator.add(
                f"{name_prefix}derivativeTaylorExpansion",
                derivativeTaylorExpansionExpr,
                target=target,
            )

    def add_include_tensors(self, include_tensors):
        super().add_include_tensors(include_tensors)
        include_tensors.add(self.db.nodes2D)
