#!/usr/bin/env python3
##
# @file
# This file is part of SeisSol.
#
# @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
#
# @section LICENSE
# Copyright (c) 2016-2018, SeisSol Group
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# @section DESCRIPTION
#
  
import numpy as np
from yateto import Tensor, Scalar, simpleParameterSpace, parameterSpaceFromRanges
from yateto.input import parseXMLMatrixFile, parseJSONMatrixFile, memoryLayoutFromFile
from yateto.ast.node import Add
from yateto.ast.transformer import DeduceIndices, EquivalentSparsityPattern
from yateto.memory import CSCMemoryLayout

from aderdg import ADERDGBase
from multSim import OptionalDimTensor

class Viscoelastic2ADERDG(ADERDGBase):
  def __init__(self, order, multipleSimulations, matricesDir, memLayout, numberOfMechanisms, **kwargs):
    super().__init__(order, multipleSimulations, matricesDir)

    self.numberOfMechanisms = numberOfMechanisms

    clones = {
      'star': ['star(0)', 'star(1)', 'star(2)'],
    }
    self.db.update( parseXMLMatrixFile('{}/matrices_viscoelastic.xml'.format(matricesDir), clones) )
    memoryLayoutFromFile(memLayout, self.db, clones)

    self._qShapeExtended = (self.numberOf3DBasisFunctions(), self.numberOfExtendedQuantities())
    self._qShapeAnelastic = (self.numberOf3DBasisFunctions(), self.numberOfAnelasticQuantities(), self.numberOfMechanisms)
    self.Qext = OptionalDimTensor('Qext', self.Q.optName(), self.Q.optSize(), self.Q.optPos(), self._qShapeExtended, alignStride=True)
    self.Qane = OptionalDimTensor('Qane', self.Q.optName(), self.Q.optSize(), self.Q.optPos(), self._qShapeAnelastic, alignStride=True)
    self.Iane = OptionalDimTensor('Iane', self.Q.optName(), self.Q.optSize(), self.Q.optPos(), self._qShapeAnelastic, alignStride=True)

    self.E = Tensor('E', (self.numberOfAnelasticQuantities(), self.numberOfMechanisms, self.numberOfQuantities()))
    self.w = Tensor('w', (self.numberOfMechanisms,))
    self.W = Tensor('W', (self.numberOfMechanisms, self.numberOfMechanisms), np.eye(self.numberOfMechanisms, dtype=bool), CSCMemoryLayout)

    selectElaSpp = np.zeros((self.numberOfExtendedQuantities(), self.numberOfQuantities()))
    selectElaSpp[0:self.numberOfQuantities(),0:self.numberOfQuantities()] = np.eye(self.numberOfQuantities())
    self.selectEla = Tensor('selectEla', (self.numberOfExtendedQuantities(), self.numberOfQuantities()), selectElaSpp, CSCMemoryLayout)

    selectAneSpp = np.zeros((self.numberOfExtendedQuantities(), self.numberOfAnelasticQuantities()))
    selectAneSpp[self.numberOfQuantities():self.numberOfExtendedQuantities(),0:self.numberOfAnelasticQuantities()] = np.eye(self.numberOfAnelasticQuantities())
    self.selectAne = Tensor('selectAne', (self.numberOfExtendedQuantities(), self.numberOfAnelasticQuantities()), selectAneSpp, CSCMemoryLayout)

    self.db.update(
      parseJSONMatrixFile('{}/nodal/nodalBoundary_matrices_{}.json'.format(matricesDir,
                                                                           self.order),
                          {},
                          alignStride=self.alignStride,
                          transpose=self.transpose,
                          namespace='nodal')
    )

  def numberOfQuantities(self):
    return 9

  def numberOfAnelasticQuantities(self):
    return 6

  def numberOfExtendedQuantities(self):
    """Number of quantities for fused computation of elastic and anelastic update."""
    return self.numberOfQuantities() + self.numberOfAnelasticQuantities()

  def numberOfFullQuantities(self):
    """Number of quantities when unrolling anelastic tensor into a matrix."""
    return self.numberOfQuantities() + self.numberOfMechanisms * self.numberOfAnelasticQuantities()

  def extendedQTensor(self):
    return self.Qext

  def starMatrix(self, dim):
    return self.db.star[dim]

  def addInit(self, generator):
    super().addInit(generator)

    selectElaFullSpp = np.zeros((self.numberOfFullQuantities(), self.numberOfQuantities()))
    selectElaFullSpp[0:self.numberOfQuantities(),0:self.numberOfQuantities()] = np.eye(self.numberOfQuantities())
    selectElaFull = Tensor('selectElaFull', (self.numberOfFullQuantities(), self.numberOfQuantities()), selectElaFullSpp, CSCMemoryLayout)

    selectAneFullSpp = np.zeros((self.numberOfFullQuantities(), self.numberOfAnelasticQuantities(), self.numberOfMechanisms))
    for mech in range(self.numberOfMechanisms):
      q1 = self.numberOfQuantities()+mech*self.numberOfAnelasticQuantities()
      q2 = q1 + self.numberOfAnelasticQuantities()
      selectAneFullSpp[q1:q2,:,mech] = np.eye(self.numberOfAnelasticQuantities())
    selectAneFull = Tensor('selectAneFull', (self.numberOfFullQuantities(), self.numberOfAnelasticQuantities(), self.numberOfMechanisms), selectAneFullSpp)

    iniShape = (self.numberOf3DQuadraturePoints(), self.numberOfFullQuantities())
    iniCond = OptionalDimTensor('iniCond', self.Q.optName(), self.Q.optSize(), self.Q.optPos(), iniShape, alignStride=True)
    dofsShape = (self.numberOf3DQuadraturePoints(), self.numberOfQuantities())
    dofsQP = OptionalDimTensor('dofsQP', self.Q.optName(), self.Q.optSize(), self.Q.optPos(), dofsShape, alignStride=True)

    projectIniCondEla = self.Q['kp'] <= self.db.projectQP[self.t('kl')] * iniCond['lq'] * selectElaFull['qp']
    projectIniCondAne = self.Qane['kpm'] <= self.db.projectQP[self.t('kl')] * iniCond['lq'] * selectAneFull['qpm']
    generator.add('projectIniCond', [projectIniCondEla, projectIniCondAne])
    generator.add('evalAtQP', dofsQP['kp'] <= self.db.evalAtQP[self.t('kl')] * self.Q['lp'])

  def addLocal(self, generator, targets):
    volumeSum = Add()
    for i in range(3):
      volumeSum += self.db.kDivM[i][self.t('kl')] * self.I['lq'] * self.db.star[i]['qp']
    volumeExt = (self.Qext['kp'] <= volumeSum)
    generator.add('volumeExt', volumeExt)

    localFluxExt = lambda i: self.Qext['kp'] <= self.Qext['kp'] + self.db.rDivM[i][self.t('km')] * self.db.fMrT[i][self.t('ml')] * self.I['lq'] * self.AplusT['qp']
    localFluxExtPrefetch = lambda i: self.I if i == 0 else (self.Q if i == 1 else None)
    generator.addFamily('localFluxExt', simpleParameterSpace(4), localFluxExt, localFluxExtPrefetch)

    generator.add('local', [
      self.Qane['kpm'] <= self.Qane['kpm'] + self.w['m'] * self.Qext['kq'] * self.selectAne['qp'] + self.Iane['kpl'] * self.W['lm'],
      self.Q['kp'] <= self.Q['kp'] + self.Qext['kq'] * self.selectEla['qp'] + self.Iane['kqm'] * self.E['qmp']
    ])

  def addNeighbor(self, generator, targets):
    neighbourFluxExt = lambda h,j,i: self.Qext['kp'] <= self.Qext['kp'] + self.db.rDivM[i][self.t('km')] * self.db.fP[h][self.t('mn')] * self.db.rT[j][self.t('nl')] * self.I['lq'] * self.AminusT['qp']
    neighbourFluxExtPrefetch = lambda h,j,i: self.I
    generator.addFamily('neighbourFluxExt', simpleParameterSpace(3,4,4), neighbourFluxExt, neighbourFluxExtPrefetch)

    generator.add('neighbour', [
      self.Qane['kpm'] <= self.Qane['kpm'] + self.w['m'] * self.Qext['kq'] * self.selectAne['qp'],
      self.Q['kp'] <= self.Q['kp'] + self.Qext['kq'] * self.selectEla['qp']
    ])

  def addTime(self, generator, targets):
    qShape = (self.numberOf3DBasisFunctions(), self.numberOfQuantities())
    dQ = [OptionalDimTensor('dQ({})'.format(d), self.Q.optName(), self.Q.optSize(), self.Q.optPos(), qShape, alignStride=True) for d in range(self.order)]
    self.dQs = dQ
    dQext = [OptionalDimTensor('dQext({})'.format(d), self.Q.optName(), self.Q.optSize(), self.Q.optPos(), self._qShapeExtended, alignStride=True) for d in range(self.order)]
    dQane = [OptionalDimTensor('dQane({})'.format(d), self.Q.optName(), self.Q.optSize(), self.Q.optPos(), self._qShapeAnelastic, alignStride=True) for d in range(self.order)]

    power = Scalar('power')

    derivativeTaylorExpansionEla = lambda d: (self.I['kp'] <= self.I['kp'] + power * dQ[d]['kp']) if d > 0 else (self.I['kp'] <= power * dQ[0]['kp'])
    derivativeTaylorExpansionAne = lambda d: (self.Iane['kpm'] <= self.Iane['kpm'] + power * dQane[d]['kpm']) if d > 0 else (self.Iane['kpm'] <= power * dQane[0]['kpm'])

    def derivative(kthDer):
      derivativeSum = Add()
      for j in range(3):
        derivativeSum += self.db.kDivMT[j][self.t('kl')] * dQ[kthDer-1]['lq'] * self.db.star[j]['qp']
      return derivativeSum

    generator.addFamily('derivative', parameterSpaceFromRanges(range(1,self.order)), lambda d: [
      dQext[d]['kp'] <= derivative(d),
      dQ[d]['kp'] <= dQext[d]['kq'] * self.selectEla['qp'] + dQane[d-1]['kqm'] * self.E['qmp'],
      dQane[d]['kpm'] <= self.w['m'] * dQext[d]['kq'] * self.selectAne['qp'] + dQane[d-1]['kpl'] * self.W['lm']
    ])
    generator.addFamily('derivativeTaylorExpansion', simpleParameterSpace(self.order), lambda d: [
      derivativeTaylorExpansionEla(d),
      derivativeTaylorExpansionAne(d)
    ])
    generator.addFamily('derivativeTaylorExpansionEla', simpleParameterSpace(self.order), derivativeTaylorExpansionEla)

  def add_include_tensors(self, include_tensors):
    super().add_include_tensors(include_tensors)
    include_tensors.add(self.db.nodes2D)
    # Nodal flux kernel uses this matrix but is not supported by visco2
    include_tensors.update([
      self.db.project2nFaceTo3m[i] for i in range(4)
    ])