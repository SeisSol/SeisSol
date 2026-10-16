# SPDX-FileCopyrightText: 2019 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff

from kernels.common import cold_kernel_attrs, generate_kernel_name_prefix
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
            attrs=cold_kernel_attrs(),
        )
        generator.add(
            "evalAtQP",
            dofsQP["kp"] <= self.db.evalAtQP[self.t("kl")] * self.Q["lp"],
            attrs=cold_kernel_attrs(),
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
                lambda i: self.extendedQTensor()["kp"]
                <= self.extendedQTensor()["kp"]
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
        qShape = (
            self.num3DBasisFunctions(),
            self.numQuantities(),
        )

        def derivative(target, writeDerivatives):
            """The predictor: the derivatives dQ and, interleaved with them, the
            time integral I of their Taylor expansion. Without
            writeDerivatives, the derivatives other than the DOFs in dQ(0) are
            temporaries, i.e. only I leaves the kernel."""
            dQ0 = OptionalDimTensor(
                "dQ(0)",
                self.Q.optName(),
                self.Q.optSize(),
                self.Q.optPos(),
                qShape,
                alignStride=True,
            )

            # on GPUs, the kernel reads the DOFs as they are, and copies them
            # into dQ(0) only if it writes the derivatives
            dQ0True = self.Q if target == "gpu" else dQ0

            derivatives = [dQ0True]
            dQs = [dQ0]

            derivativeExpr = [self.I["kp"] <= powers[0] * dQ0True["kp"]]

            if target == "gpu" and writeDerivatives:
                derivativeExpr += [dQ0["kp"] <= self.Q["kp"]]

            for i in range(1, self.order):
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
                    temporary=not writeDerivatives,
                )
                dQs.append(dQ)

                # for now, we interleave derivative
                # and derivativeTaylorExpansion kernels.
                derivativeExpr += [
                    dQ["kp"] <= derivativeSum,
                    self.I["kp"] <= self.I["kp"] + powers[i] * dQ["kp"],
                ]

                derivatives.append(dQ)

            return derivativeExpr, dQs

        # the derivatives as written out, for the kernels that read them
        _, self.dQs = derivative("cpu", True)

        for target in targets:
            name_prefix = generate_kernel_name_prefix(target)

            # member 0 writes the derivatives, member 1 only the time integral;
            # where the kernel would keep them in memory, it has member 0 alone
            # (whether a device uses member 1 is up to Spacetime)
            generator.addFamily(
                f"{name_prefix}derivative",
                simpleParameterSpace(2 if self.keepsTemporaries(target) else 1),
                lambda i: derivative(target, i == 0)[0],
                target=target,
            )

            derivativeTaylorExpansion = powers[0] * self.dQs[0]["kp"]
            for i in range(1, self.order):
                derivativeTaylorExpansion += powers[i] * self.dQs[i]["kp"]

            derivativeTaylorExpansionExpr = self.I["kp"] <= derivativeTaylorExpansion
            generator.add(
                f"{name_prefix}derivativeTaylorExpansion",
                derivativeTaylorExpansionExpr,
                target=target,
            )

    def add_include_tensors(self, include_tensors):
        super().add_include_tensors(include_tensors)
        include_tensors.add(self.db.nodes2D)
