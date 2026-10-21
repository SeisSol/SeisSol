# SPDX-FileCopyrightText: 2019 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff

import numpy as np
from kernels import material
from kernels.common import generate_kernel_name_prefix
from kernels.multsim import OptionalDimTensor
from yateto import Scalar, ops, simpleParameterSpace
from yateto.ast.node import Accumulate
from yateto.ast.transformer import DeduceIndices, EquivalentSparsityPattern

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
        the operator varies inside the cell, see narrowsDerivatives: in a
        hierarchical basis the ones up to its degree, in another one all of
        them, since the polynomials of a lower degree have entries in every
        function of the cell there."""
        if not self.narrowsDerivatives() or not self.basis.hierarchical:
            return self.num3DBasisFunctions()
        degrees = self.order - k
        return degrees * (degrees + 1) * (degrees + 2) // 6

    def derivativeProjection(self, k):
        """The projection the k-th time derivative of the predictor comes back
        with: onto the polynomials of the degree it keeps, in the basis of the
        cell (material.projectionOntoDegree)."""
        if not self.narrowsDerivatives():
            return self.operatorProject
        # one tensor per member, however many targets ask for it
        if not hasattr(self, "_derivativeProjections"):
            self._derivativeProjections = {}
        if k not in self._derivativeProjections:
            self._derivativeProjections[k] = material.projectionOntoDegree(
                self._matricesDir,
                self,
                self.materialPoints,
                self.order - 1 - k,
                f"derivativeProject({k})",
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

    def addNeighbor(self, generator, targets):
        if self.weightAdjustable():
            # the end of the neighbour phase, where a cell that lifts its fluxes
            # without material has them all (see _configureWeightAdjust)
            generator.add(
                "weightAdjust",
                self.weightAdjust(self.fluxBuffer, self.Q),
                target="cpu",
            )
            if "gpu" in targets:
                # a device kernel writes each temporary once, see
                # singleDefinitions
                with self.singleDefinitions():
                    generator.add(
                        "gpu_weightAdjust",
                        self.weightAdjust(self.fluxBuffer, self.Q),
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

                    lowerStatements, terms = self.modalDerivative(derivatives[-1])
                    derivativeSum = Accumulate(ops.Add())
                    if self.sourceMatrix():
                        derivativeSum += (
                            derivatives[-1]["kq"] * self.sourceMatrix()["qp"]
                        )
                    for term in terms:
                        derivativeSum += term

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
                    derivativeExpr += lowerStatements + [
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
