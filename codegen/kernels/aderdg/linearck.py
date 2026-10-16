# SPDX-FileCopyrightText: 2019 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff

import numpy as np
from kernels.common import generate_kernel_name_prefix
from kernels.multsim import OptionalDimTensor
from yateto import Scalar, Tensor, ops, simpleParameterSpace
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

    def narrowsDerivatives(self):
        """A constant operator drops a degree with every derivative, so its
        k-th time derivative has the degree order - 1 - k. One that varies
        inside the cell raises the degree again, by as much as it varies. The
        field is accurate to the order of the scheme and its k-th derivative
        to k orders less, which the degree order - 1 - k resolves already: the
        modes above it are of the size of that error, and in the predictor
        they come with the k-th power of the time step, which brings them to
        the order of the scheme. Where the build asks for it
        (MATERIAL_NARROW_DERIVATIVES), the k-th derivative keeps the modes up
        to that degree, projected onto those alone.

        A source term keeps the degree of what it acts on, and its derivatives
        then have every mode even where the operator is constant, so a solver
        with one keeps them all."""
        return (
            self.nodalMaterial
            and self.narrowDerivatives
            and self.sourceMatrix() is None
        )

    def fusesVolume(self, target):
        """Whether the derivative kernel of a target forms the volume term of
        the corrector as well.

        The volume term is the operator applied to the time integral of the
        predictor, the sum of power(k) A dQ(k), projected onto every mode. The
        chain forms A dQ(k) at the sample points for every k but the last, and
        projects it onto the modes the next derivative keeps. Where the chain
        narrows, the last derivative is constant and A dQ(last) vanishes, so the
        products the chain forms are all the volume term needs: summed at the
        points with the weights of the time integral and projected once onto
        every mode, they are the volume term, and the time integral is not read
        at the points again. The host hands the result from the predictor to
        the corrector of the same cell (LocalTmp::volumeFused).

        The device runs the predictor and the corrector of a batch as two
        kernels without a buffer per cell between them, so its corrector forms
        the volume term itself."""
        return target == "cpu" and self.narrowsDerivatives()

    def volumeFused(self):
        """The volume term the derivative kernel hands to the corrector, and the
        sum at the points it is projected from; see fusesVolume."""
        if not hasattr(self, "_volumeFused"):
            self._volumeFused = OptionalDimTensor(
                "volumeFused",
                self.Q.optName(),
                self.Q.optSize(),
                self.Q.optPos(),
                (self.num3DBasisFunctions(), self.numQuantities()),
                alignStride=True,
            )
            self._volumeAtPoints = self.nodalTemporary(
                "volumeAtPoints", self.nodalProductShape
            )
        return self._volumeFused, self._volumeAtPoints

    def derivativeModes(self, k):
        """How many modes the k-th time derivative of the predictor keeps where
        the operator varies inside the cell, see narrowsDerivatives."""
        if not self.narrowsDerivatives():
            return self.num3DBasisFunctions()
        degrees = self.order - k
        return degrees * (degrees + 1) * (degrees + 2) // 6

    def derivativeProjection(self, k):
        """The projection the k-th time derivative of the predictor comes back
        with: onto the modes it keeps. The basis is orthogonal and ordered by
        degree, so that is the operator's projection cut to its leading rows."""
        modes = self.derivativeModes(k)
        if modes == self.num3DBasisFunctions():
            return self.operatorProject
        # one tensor per member, however many targets ask for it
        if not hasattr(self, "_derivativeProjections"):
            self._derivativeProjections = {}
        if k not in self._derivativeProjections:
            self._derivativeProjections[k] = Tensor(
                f"derivativeProject({k})",
                self.operatorProject.shape(),
                spp={
                    index: value
                    for index, value in self.operatorProject.values().items()
                    if index[0] < modes
                },
                alignStride=self.alignStride("derivativeProject"),
            )
        return self._derivativeProjections[k]

    def derivativePattern(self, k):
        """Where the k-th time derivative of the predictor has entries where the
        operator varies inside the cell: the modes it keeps, every quantity.
        None where it keeps every mode."""
        modes = self.derivativeModes(k)
        if modes == self.num3DBasisFunctions():
            return None
        rows = self.Q.insertOptDim((slice(0, modes), slice(None)), (slice(None),))
        pattern = np.zeros(self.Q.shape(), dtype=bool)
        pattern[rows] = True
        return pattern

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
            # a device kernel writes each temporary once, see singleDefinitions
            with self.singleDefinitions(target == "gpu"):
                if self.fusesVolume(target):
                    # the derivative kernel formed the volume term already
                    volumeFused, _ = self.volumeFused()
                    volumeExpr = [self.Q["kp"] <= self.Q["kp"] + volumeFused["kp"]]
                elif self.nodalMaterial:
                    volumeUpdate = OptionalDimTensor(
                        "volumeUpdate",
                        self.Q.optName(),
                        self.Q.optSize(),
                        self.Q.optPos(),
                        (self.num3DBasisFunctions(), self.numQuantities()),
                        alignStride=True,
                        temporary=True,
                    )
                    volumeExpr = self.nodalApply(
                        self.I, volumeUpdate, self.nodalVolumeDerivatives()
                    )
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
                if not self.fusesVolume(target):
                    volumeExpr = self.starAssembly() + volumeExpr
                generator.add(f"{name_prefix}volume", volumeExpr, target=target)

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

        if self.nodalFaceFlux:
            localFlux = lambda i: self.nodalFlux(
                self.I,
                self.Q,
                self.faceEval(i),
                self.faceLift(i),
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
                    self.faceEval(i),
                    self.faceLift(i),
                    self.fluxCoefficientsLocal,
                )
            else:
                localFlux = (
                    lambda i: self.Q["kp"]
                    <= self.Q["kp"]
                    + plusFluxMatrixAccessor(i) * self.I["lq"] * self.AplusT["qp"]
                )
            with self.singleDefinitions(target == "gpu"):
                if target == "gpu":
                    generator.addFamily(
                        f"{name_prefix}localFlux",
                        simpleParameterSpace(4),
                        localFlux,
                        target=target,
                    )

                if self.nodalFaceFlux:
                    localFluxAll = self.nodalLocalFluxAll(self.I, self.Q)
                else:
                    localFluxAll = self.Q["kp"] <= sum(
                        [
                            plusFluxMatrixAccessor(i)
                            * self.I["lq"]
                            * self.AplusTAll[i]["qp"]
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
                self.faceLift(i),
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

        if self.weightAdjustable():
            # the end of the neighbour phase, where a cell that lifts its fluxes
            # without material has them all (see _configureWeightAdjust)
            generator.add(
                "weightAdjust",
                self.weightAdjust(self.fluxBuffer, self.Q),
                target="cpu",
            )

        if "gpu" in targets and self.nodalFaceFlux:
            # the same statements as on the host; only the premultiplied
            # matrices below are a device matter, and the nodal form has none
            with self.singleDefinitions():
                generator.addFamily(
                    "gpu_neighboringFlux",
                    simpleParameterSpace(4, 4),
                    neighborFlux,
                    target="gpu",
                )
        elif "gpu" in targets:
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
            with self.singleDefinitions(target == "gpu"):
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

                fused = self.fusesVolume(target)
                if fused:
                    volumeFused, volumeAtPoints = self.volumeFused()

                for i in range(1, self.order):
                    power = powers[i]
                    if self.nodalMaterial:
                        # The operator varies inside the cell, so the degree of
                        # a derivative is what the build lets it keep, not what
                        # the generator can deduce (narrowsDerivatives).
                        dQ = OptionalDimTensor(
                            "dQ({})".format(i),
                            self.Q.optName(),
                            self.Q.optSize(),
                            self.Q.optPos(),
                            qShape,
                            spp=self.derivativePattern(i),
                            alignStride=True,
                        )
                        self.dQs.append(dQ)
                        # the product at the points enters the volume term
                        # with the weight of the derivative it is formed from
                        derivativeExpr += self.nodalApply(
                            derivatives[-1],
                            dQ,
                            self.db.kDivMT,
                            project=self.derivativeProjection(i),
                            pointSum=(
                                (volumeAtPoints, powers[i - 1], i == 1)
                                if fused
                                else None
                            ),
                        )
                        derivativeExpr += self.sourceTerm(derivatives[-1], dQ)
                        derivativeExpr += [
                            self.I["kp"] <= self.I["kp"] + power * dQ["kp"]
                        ]
                        derivativeTaylorExpansion += power * dQ["kp"]
                        derivatives.append(dQ)
                        continue

                    derivativeSum = Accumulate(ops.Add())
                    if self.sourceMatrix():
                        derivativeSum += (
                            derivatives[-1]["kq"] * self.sourceMatrix()["qp"]
                        )
                    for j in range(3):
                        derivativeSum += (
                            self.db.kDivMT[j][self.t("kl")]
                            * derivatives[-1]["lq"]
                            * self.starMatrix(j)["qp"]
                        )

                    derivativeSum = DeduceIndices(self.Q["kp"].indices).visit(
                        derivativeSum
                    )
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

                if fused:
                    derivativeExpr += [
                        volumeFused["kp"]
                        <= self.operatorProject["kn"] * volumeAtPoints["np"]
                    ]

                derivativeTaylorExpansionExpr = (
                    self.I["kp"] <= derivativeTaylorExpansion
                )
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
