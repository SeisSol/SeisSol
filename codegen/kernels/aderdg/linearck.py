# SPDX-FileCopyrightText: 2019 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff

from kernels.common import generate_kernel_name_prefix
from kernels.multsim import OptionalDimTensor
from yateto import Scalar, ops, simpleParameterSpace
from yateto.ast.node import Accumulate
from yateto.ast.transformer import DeduceIndices, EquivalentSparsityPattern
from yateto.util import (
    tensor_collection_from_constant_expression,
)

from .aderdg import ADERDGBase


class LinearCK(ADERDGBase):
    def name(self):
        return "linearck"

    def sourceMatrix(self):
        return None

    def extendedQTensor(self):
        return self.Q

    def numExtendedQuantities(self):
        return self.numQuantities()

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

    def addLocal(self, generator, targets):
        for target in targets:
            name_prefix = generate_kernel_name_prefix(target)
            if self.nodalMaterial:
                volumeUpdate = OptionalDimTensor(
                    "volumeUpdate",
                    self.Q.optName(),
                    self.Q.optSize(),
                    self.Q.optPos(),
                    (self.num3DBasisFunctions(), self.numQuantities()),
                    alignStride=True,
                    temporary=True,
                )
                volumeExpr = self.nodalApply(self.I, volumeUpdate, self.db.kDivM)
                volumeExpr += [self.Q["kp"] <= self.Q["kp"] + volumeUpdate["kp"]]
                volumeExpr += self.sourceTerm(self.I, self.Q)
            else:
                volumeSum = self.Q["kp"]
                for i in range(3):
                    volumeSum += (
                        self.db.kDivM[i][self.t("kl")]
                        * self.I["lq"]
                        * self.starMatrix(i)["qp"]
                    )
                if self.sourceMatrix():
                    volumeSum += self.I["kq"] * self.sourceMatrix()["qp"]
                volumeExpr = [self.Q["kp"] <= volumeSum]
            generator.add(
                f"{name_prefix}volume", self.starAssembly() + volumeExpr, target=target
            )

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

        if self.nodalFaceFlux:
            localFlux = lambda i: self.nodalFlux(
                self.I,
                self.Q,
                self.db.V3mTo2nFace[i][self.t("nl")],
                self.db.project2nFaceTo3m[i]["kn"],
                self.fluxCoefficientsLocal,
            )
        else:
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

        for target in targets:
            name_prefix = generate_kernel_name_prefix(target)
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

            if self.nodalFaceFlux:
                localFlux = lambda i: self.nodalFlux(
                    self.I,
                    self.Q,
                    self.db.V3mTo2nFace[i][self.t("nl")],
                    self.db.project2nFaceTo3m[i]["kn"],
                    self.fluxCoefficientsLocal,
                )
            else:
                localFlux = (
                    lambda i: self.Q["kp"]
                    <= self.Q["kp"]
                    + plusFluxMatrixAccessor(i) * self.I["lq"] * self.AplusT["qp"]
                )
            if target == "gpu":
                generator.addFamily(
                    f"{name_prefix}localFlux",
                    simpleParameterSpace(4),
                    localFlux,
                    target=target,
                )

            localFluxAll = self.Q["kp"] <= sum(
                [
                    plusFluxMatrixAccessor(i) * self.I["lq"] * self.AplusTAll[i]["qp"]
                    for i in range(4)
                ],
                start=self.Q["kp"],
            )
            generator.add(
                f"{name_prefix}localFluxAll",
                localFluxAll,
                target=target,
            )

    def addNeighbor(self, generator, targets):
        if self.nodalFaceFlux:
            # The canonical vertex order pins the face orientation of every
            # regular face to zero, so the neighbour's face values reach this
            # cell's ordering through the first renumbering only.
            neighborFlux = lambda j, i: self.nodalFlux(
                self.I,
                self.Q,
                self.db.neighborToFace[0, j]["nl"],
                self.db.project2nFaceTo3m[i]["kn"],
                self.fluxCoefficientsNeighbor,
            )
        else:
            neighborFlux = (
                lambda j, i: self.Q["kp"]
                <= self.Q["kp"]
                + self.db.rDivM[i][self.t("km")]
                * self.db.fPrT[j][self.t("ml")]
                * self.I["lq"]
                * self.AminusT["qp"]
            )
        neighborFluxPrefetch = lambda j, i: self.I
        generator.addFamily(
            "neighboringFlux",
            simpleParameterSpace(4, 4),
            neighborFlux,
            neighborFluxPrefetch,
            target="cpu",
        )

        if "gpu" in targets:
            minusFluxMatrixAccessor = (
                lambda j, i: self.db.rDivM[i][self.t("km")]
                * self.db.fPrT[j][self.t("ml")]
            )
            if self.kwargs["enable_premultiply_flux"]:
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

            neighborFlux = (
                lambda j, i: self.Q["kp"]
                <= self.Q["kp"]
                + minusFluxMatrixAccessor(j, i) * self.I["lq"] * self.AminusT["qp"]
            )
            generator.addFamily(
                "gpu_neighboringFlux",
                simpleParameterSpace(4, 4),
                neighborFlux,
                target="gpu",
            )

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

            # for now, interleave Taylor expansion and derivative computation
            derivativeExpr = self.starAssembly() + [
                self.I["kp"] <= power * dQ0True["kp"]
            ]
            derivativeTaylorExpansion = power * dQ0["kp"]

            if target == "gpu":
                derivativeExpr += [dQ0["kp"] <= self.Q["kp"]]

            self.dQs = [dQ0]

            for i in range(1, self.order):
                power = powers[i]
                if self.nodalMaterial:
                    # A constant operator drops a degree with every derivative,
                    # so each one occupies fewer modes than the last and the
                    # chain narrows. A material that varies inside the cell
                    # raises the degree again by as much as it carries itself,
                    # so nothing narrows and every derivative stays full.
                    dQ = OptionalDimTensor(
                        "dQ({})".format(i),
                        self.Q.optName(),
                        self.Q.optSize(),
                        self.Q.optPos(),
                        qShape,
                        alignStride=True,
                    )
                    self.dQs.append(dQ)
                    derivativeExpr += self.nodalApply(
                        derivatives[-1], dQ, self.db.kDivMT
                    )
                    derivativeExpr += self.sourceTerm(derivatives[-1], dQ)
                    derivativeExpr += [self.I["kp"] <= self.I["kp"] + power * dQ["kp"]]
                    derivativeTaylorExpansion += power * dQ["kp"]
                    derivatives.append(dQ)
                    continue

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
