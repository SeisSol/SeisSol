# SPDX-FileCopyrightText: 2019 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff

import numpy as np
from kernels import material
from kernels.basis import count
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

    def chainDropsDegree(self):
        """A constant operator without a source term drops a degree with every
        derivative."""
        return not self.nodalMaterial and self.sourceMatrix() is None

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
        """Whether the predictor of a target hands the corrector its volume term
        (LocalTmp::volumeFused), so that the volume kernel only adds it: on the
        host, where the derivative chain forms it (chainFormsVolume). The
        device has no volume kernel (predictorAppliesVolume)."""
        return target == "cpu" and self.chainFormsVolume()

    def chainFormsVolume(self):
        """Whether the derivative chain forms the volume term of the corrector
        as well.

        The volume term is the operator applied to the time integral of the
        predictor, the sum of power(k) A dQ(k), projected onto every mode. The
        chain forms A dQ(k) at the sample points for every k but the last, and
        projects it onto the modes the next derivative keeps. Where the chain
        narrows, the last derivative is constant and A dQ(last) vanishes, so the
        products the chain forms are all the volume term needs: summed at the
        points with the weights of the time integral and projected once onto
        every mode, they are the volume term, and the time integral is not read
        at the points again. The host hands the result from the predictor to
        the corrector of the same cell (LocalTmp::volumeFused); the device adds
        it to the field in the predictor (predictorAppliesVolume).

        Where the operator is one per cell and there is no source term
        (chainDropsDegree), every derivative the chain forms is exact: A dQ(k)
        is dQ(k + 1), and the volume term is the sum of power(k) dQ(k + 1), the
        time integral of the time derivative -- the Taylor expansion at the end
        of the step less the field at its start."""
        return self.narrowsDerivatives() or self.chainDropsDegree()

    def volumeStatements(self):
        """The volume term, applied to the time integral and added to the
        field: a step of the chain (volumeDerivatives) and the source term,
        which a cell carries as one matrix where its material does not vary and
        puts together at its sample points where it does. The operators are the
        ones the chain applies, assembled by starAssembly. A predictor that
        applies the volume term (predictorAppliesVolume) adds these where its
        chain does not form the term (chainFormsVolume)."""
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
            statements = self.nodalApply(self.I, volumeUpdate, self.volumeDerivatives())
            statements += [self.Q["kp"] <= self.Q["kp"] + volumeUpdate["kp"]]
            return statements + self.sourceTerm(self.I, self.Q)
        statements, terms = self.modalDerivative(self.I)
        volumeSum = self.Q["kp"]
        for term in terms:
            volumeSum += term
        if self.sourceMatrix() is not None:
            volumeSum += self.I["kq"] * self.sourceMatrix()["qp"]
        return statements + [self.Q["kp"] <= volumeSum]

    def chainVolume(self, target, first, derivatives, powers, volumeAtPoints):
        """target = first (where there is one) + the volume term the chain
        forms (chainFormsVolume): the sum at the points projected onto every
        mode where the chain narrows, and the time integral of the time
        derivative, power(k - 1) dQ(k) summed, where it drops a degree with
        every derivative -- by Horner's scheme where the derivatives are kept in
        the bases of their degrees."""
        if self.nodalMaterial:
            projected = self.operatorProject["kn"] * volumeAtPoints["np"]
            return [target["kp"] <= (projected if first is None else first + projected)]
        shifted = [None] + powers[:-1]
        if self.derivativesByDegree():
            return self.hornerTaylor(target, first, derivatives, shifted, "volumeSum")
        volumeSum = first
        for k in range(1, len(derivatives)):
            term = shifted[k] * derivatives[k]["kp"]
            volumeSum = term if volumeSum is None else volumeSum + term
        return [target["kp"] <= volumeSum]

    def volumeFused(self):
        """The volume term the predictor hands to the corrector, and the sum at
        the points the chain projects it from; see fusesVolume."""
        if not hasattr(self, "_volumeFused"):
            self._volumeFused = OptionalDimTensor(
                "volumeFused",
                self.Q.optName(),
                self.Q.optSize(),
                self.Q.optPos(),
                (self.num3DBasisFunctions(), self.numQuantities()),
                alignStride=True,
            )
            # the chain of a varying operator sums its products at the points
            self._volumeAtPoints = (
                self.nodalTemporary("volumeAtPoints", self.nodalProductShape)
                if self.nodalMaterial
                else None
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
            if not self.predictorAppliesVolume(target):
                with self.singleDefinitions(target == "gpu"):
                    if self.fusesVolume(target):
                        # the derivative kernel formed the volume term already
                        volumeFused, _ = self.volumeFused()
                        volumeExpr = [self.Q["kp"] <= self.Q["kp"] + volumeFused["kp"]]
                    else:
                        volumeExpr = self.starAssembly() + self.volumeStatements()
                    generator.add(f"{name_prefix}volume", volumeExpr, target=target)

            if self.curvilinear and target == "cpu":
                # the corrector of a curved cell forms its volume term through
                # its weight, in the split form, and adds the normal flux of its
                # own trace back for the momentum rows (_configureCurvedSplit)
                generator.add(
                    "splitVolume",
                    self.splitVolumeStatements(self.I, self.fluxBuffer),
                    target=target,
                )
                generator.addFamily(
                    "ownTraction",
                    simpleParameterSpace(4, 2),
                    lambda face, side: self.ownTractionStatements(
                        face, side, self.faceTrace, self.fluxBuffer
                    ),
                    target=target,
                )

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

                # derivatives in the bases of their degrees, summed by Horner's
                # scheme once the chain is done (derivativesByDegree); else
                # the Taylor expansion is interleaved with the chain
                byDegree = self.derivativesByDegree()

                derivativeExpr = self.starAssembly()
                if not byDegree:
                    derivativeExpr += [self.I["kp"] <= power * dQ0True["kp"]]
                derivativeTaylorExpansion = power * dQ0["kp"]

                if target == "gpu":
                    derivativeExpr += [dQ0["kp"] <= self.Q["kp"]]

                self.dQs = [dQ0]

                fused = self.chainFormsVolume()
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

                    if byDegree:
                        # derivative i in the basis of its degree, from the one
                        # before, of one degree more (derivativesByDegree)
                        degree = self.order - i
                        dQ = OptionalDimTensor(
                            "dQ({})".format(i),
                            self.Q.optName(),
                            self.Q.optSize(),
                            self.Q.optPos(),
                            qShape,
                            spp=self.leadingRows(count(degree - 1)),
                            alignStride=True,
                        )
                        self.dQs.append(dQ)
                        derivativeSum = None
                        for j in range(3):
                            term = (
                                self.derivativeOfDegree(degree, j)[self.t("kl")]
                                * derivatives[-1]["lq"]
                                * self.starMatrix(j)["qp"]
                            )
                            derivativeSum = (
                                term if derivativeSum is None else derivativeSum + term
                            )
                        derivativeExpr += [dQ["kp"] <= derivativeSum]
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

                if fused and not self.predictorAppliesVolume(target):
                    # into the buffer the host hands its corrector
                    derivativeExpr += self.chainVolume(
                        volumeFused, None, derivatives, powers, volumeAtPoints
                    )

                if byDegree:
                    derivativeExpr += self.hornerTaylor(
                        self.I,
                        powers[0] * dQ0True["kp"],
                        derivatives,
                        powers,
                        "taylorSum",
                    )
                    derivativeTaylorExpansionExpr = self.hornerTaylor(
                        self.I,
                        powers[0] * dQ0["kp"],
                        self.dQs,
                        powers,
                        "taylorExpansion",
                    )
                else:
                    derivativeTaylorExpansionExpr = (
                        self.I["kp"] <= derivativeTaylorExpansion
                    )

                if self.predictorAppliesVolume(target):
                    # into the field, last: the chain reads the field as its
                    # first derivative
                    derivativeExpr += (
                        self.chainVolume(
                            self.Q, self.Q["kp"], derivatives, powers, volumeAtPoints
                        )
                        if fused
                        else self.volumeStatements()
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
