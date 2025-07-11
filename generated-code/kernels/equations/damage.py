# SPDX-FileCopyrightText: 2016 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff

from kernels.equations.elastic import ElasticADERDG as ADERDGBase
from kernels.common import generate_kernel_name_prefix

from yateto.input import memoryLayoutFromFile, parseXMLMatrixFile, parseJSONMatrixFile

from kernels.multsim import OptionalDimTensor
from yateto import Scalar, simpleParameterSpace
from yateto.ast.node import Add
from yateto.ast.transformer import DeduceIndices, EquivalentSparsityPattern


class DamageADERDG(ADERDGBase):
    def __init__(self, order, multipleSimulations, matricesDir, memLayout, **kwargs):
        super().__init__(order, multipleSimulations, matricesDir, memLayout)
        clones = {
            "star": ["star(0)", "star(1)", "star(2)"],
        }
        self.db.update(parseXMLMatrixFile("{}/matrices_damage.xml".format(matricesDir), clones))
        self.db.update(parseXMLMatrixFile('{}/plasticity_{}_matrices_{}.xml'.format(matricesDir, 'nb', self.order), clones=dict(), alignStride=self.alignStride))
        numberOfNodes = self.t(self.db.v.shape())[0]

        # For cell-average
        self.db.update(parseXMLMatrixFile('{}/phi_ave_{}.xml'.format(matricesDir, order), transpose=self.transpose, alignStride=self.alignStride))
        qAveShape = (self.numberOfQuantities(),)
        self.QAve = OptionalDimTensor('QAve', 's', multipleSimulations, 0, qAveShape, alignStride=True)
        self.db.update(parseJSONMatrixFile(f'{matricesDir}/dr_stroud_matrices_{self.order}.json', clones, alignStride=self.alignStride, transpose=self.transpose))

        # For storing nodal values
        qNShape = (numberOfNodes, self.numberOfQuantities())
        self.QNodal = OptionalDimTensor('QNodal', self.Q.optName(), self.Q.optSize(), self.Q.optPos(), qNShape, alignStride=True)
        dqMShape = (numberOfNodes, self.numberOfQuantities())
        self.dQModal = OptionalDimTensor('dQModal', self.Q.optName(), self.Q.optSize(), self.Q.optPos(), dqMShape, alignStride=True)
        FNShape = (numberOfNodes, self.numberOfQuantities())
        self.FNodal = OptionalDimTensor('FNodal', self.Q.optName(), self.Q.optSize(), self.Q.optPos(), FNShape, alignStride=True)

        memoryLayoutFromFile(memLayout, self.db, clones)
        self.kwargs = kwargs

    def numberOfQuantities(self):
        return 10

    def name(self):
        return "damage"

    def addLocal(self, generator, targets):
        super().addLocal(generator, targets)

        # For cell average
        generator.add('cellAve', self.QAve['p'] <= self.db.phiAve[self.t('l')] * self.Q['lp'] * 6.0)

        # Volumetric flux integration
        qNShape = (self.t(self.db.v.shape())[0], self.numberOfQuantities())
        Tweight = Scalar('Tweight')
        for i in range(0, self.order):
            QTNodal = OptionalDimTensor('QTNodal({})'.format(i), self.Q.optName(), self.Q.optSize(), self.Q.optPos(), qNShape, alignStride=True)
            generator.add(f'damageIntegration({i})', self.Q['kp'] <= self.Q['kp'] + self.db.vInv[self.t('kl')] * QTNodal['lp'] * Tweight)

        # Kernels for volumetric flux integration
        gradXiEtaZetaX0 = Scalar('gradXiEtaZetaX0')
        gradXiEtaZetaY0 = Scalar('gradXiEtaZetaY0')
        gradXiEtaZetaZ0 = Scalar('gradXiEtaZetaZ0')
        gradXiEtaZetaX1 = Scalar('gradXiEtaZetaX1')
        gradXiEtaZetaY1 = Scalar('gradXiEtaZetaY1')
        gradXiEtaZetaZ1 = Scalar('gradXiEtaZetaZ1')
        gradXiEtaZetaX2 = Scalar('gradXiEtaZetaX2')
        gradXiEtaZetaY2 = Scalar('gradXiEtaZetaY2')
        gradXiEtaZetaZ2 = Scalar('gradXiEtaZetaZ2')
        # TweightN = Scalar('TweightN')

        numberOfNodes = self.t(self.db.v.shape())[0]
        fluxVolShape = (numberOfNodes, self.numberOfQuantities())

        for i in range(0, self.order):
            FluxVolX = OptionalDimTensor('FluxVolX({})'.format(i), self.Q.optName(), self.Q.optSize(), self.Q.optPos(), fluxVolShape, alignStride=True)
            FluxVolY = OptionalDimTensor('FluxVolY({})'.format(i), self.Q.optName(), self.Q.optSize(), self.Q.optPos(), fluxVolShape, alignStride=True)
            FluxVolZ = OptionalDimTensor('FluxVolZ({})'.format(i), self.Q.optName(), self.Q.optSize(), self.Q.optPos(), fluxVolShape, alignStride=True)

            volumeNonl = (self.Q['kp'] <= self.Q['kp'] + (self.db.kDivM[0][self.t('kl')]
                          * self.db.vInv[self.t('lm')]
                          * (gradXiEtaZetaX0 * FluxVolX['mp']
                          + gradXiEtaZetaY0 * FluxVolY['mp']
                          + gradXiEtaZetaZ0 * FluxVolZ['mp'])
                          + self.db.kDivM[1][self.t('kl')]
                          * self.db.vInv[self.t('lm')]
                          * (gradXiEtaZetaX1 * FluxVolX['mp']
                          + gradXiEtaZetaY1 * FluxVolY['mp']
                          + gradXiEtaZetaZ1 * FluxVolZ['mp'])
                          + self.db.kDivM[2][self.t('kl')]
                          * self.db.vInv[self.t('lm')]
                          * (gradXiEtaZetaX2 * FluxVolX['mp']
                          + gradXiEtaZetaY2 * FluxVolY['mp']
                          + gradXiEtaZetaZ2 * FluxVolZ['mp']))  # *TweightN
                          )
            generator.add(f'nonlinearVolumeIntegration({i})', volumeNonl)

        # Kernel for Interpolating on the surface quadrature points
        numberOfPoints = self.db.resample.shape()[0]
        gShape = (numberOfPoints, self.numberOfQuantities())
        QInterpolated = OptionalDimTensor('QInterpolated', self.Q.optName(), self.Q.optSize(), self.Q.optPos(), gShape, alignStride=True)

        def interpolateQGenerator(i, h):
            return QInterpolated['kp'] <= self.db.V3mTo2n[i, h][self.t('kl')] * self.Q['lp']  # *self.TinvT['qp']
        # interpolateQPrefetch = lambda i,h: QInterpolated
        generator.addFamily('nonlEvaluateAndRotateQAtInterpolationPoints',
                            simpleParameterSpace(4, 4),
                            interpolateQGenerator)

        # Surface integration based on the Rusanov flux values on surface quadrature points.
        fluxScale = Scalar('fluxScale')
        Flux = OptionalDimTensor('Flux', self.Q.optName(), self.Q.optSize(), self.Q.optPos(), gShape, alignStride=True)

        nodalFluxGenerator = lambda i, h: self.extendedQTensor()['kp'] <= self.extendedQTensor()['kp'] + fluxScale * self.db.V3mTo2nTWDivM[i, h][self.t('kl')] * Flux['lp']  # *self.TT['qp']
        # nodalFluxPrefetch = lambda i,h: self.I

        generator.addFamily('nonlinearSurfaceIntegral',
                            simpleParameterSpace(4, 4),
                            nodalFluxGenerator)

    def addNeighbor(self, generator, targets):
        super().addNeighbor(generator, targets)

    def addTime(self, generator, targets):
        super().addTime(generator, targets)
        # Kernel 1
        generator.add('damageConvertToNodal', self.QNodal['kp'] <= self.db.v[self.t('kl')] * self.Q['lp'])
        # Kernel 2
        generator.add('damageAssignFToDQ', self.dQModal['kp'] <= self.db.vInv[self.t('kl')] * self.FNodal['lp'])
        # Kernel 3
        powers = [Scalar(f"nlPower({i})") for i in range(self.order)]
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
            generator.add(f"{name_prefix}derivativeDamage", derivativeExpr, target=target)
            generator.add(
                f"{name_prefix}derivativeDamageTaylorExpansion",
                derivativeTaylorExpansionExpr,
                target=target,
            )

    def add_include_tensors(self, include_tensors):
        super().add_include_tensors(include_tensors)


EQUATION_CLASS = DamageADERDG
