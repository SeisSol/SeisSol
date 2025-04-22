# SPDX-FileCopyrightText: 2016 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff

from kernels.equations.elastic import ElasticADERDG as ADERDGBase
from kernels.common import generate_kernel_name_prefix

from yateto.input import memoryLayoutFromFile, parseXMLMatrixFile

from kernels.multsim import OptionalDimTensor
from yateto import Scalar, Tensor
from yateto.ast.node import Add
from yateto.ast.transformer import DeduceIndices, EquivalentSparsityPattern


class DamageADERDG(ADERDGBase):
    def __init__(self, order, multipleSimulations, matricesDir, memLayout, **kwargs):
        super().__init__(order, multipleSimulations, matricesDir, memLayout)
        clones = {
            "star": ["star(0)", "star(1)", "star(2)"],
        }
        self.db.update(parseXMLMatrixFile("{}/matrices_damage.xml".format(matricesDir), clones))
        self.db.update( parseXMLMatrixFile('{}/plasticity_{}_matrices_{}.xml'.format(matricesDir, 'nb', self.order), clones=dict(), alignStride=self.alignStride) )
        numberOfNodes = self.t(self.db.v.shape())[0]

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

    def addNeighbor(self, generator, targets):
        super().addNeighbor(generator, targets)
    
    def addTime(self, generator, targets):
        super().addTime(generator, targets)
        ## Kernel 1
        generator.add('damageConvertToNodal', self.QNodal['kp'] <= self.db.v[self.t('kl')] * self.Q['lp'] )
        ## Kernel 2
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
