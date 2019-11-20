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
  
from yateto import Tensor, Scalar, simpleParameterSpace
from yateto.input import parseXMLMatrixFile, parseJSONMatrixFile, memoryLayoutFromFile
from yateto.ast.node import Add
from yateto.ast.transformer import DeduceIndices, EquivalentSparsityPattern

from aderdg import LinearADERDG
from multSim import OptionalDimTensor

class ElasticADERDG(LinearADERDG):
  def __init__(self, order, multipleSimulations, matricesDir, memLayout, **kwargs):
    super().__init__(order, multipleSimulations, matricesDir)
    clones = {
      'star': ['star(0)', 'star(1)', 'star(2)'],
    }
    self.db.update(
      parseXMLMatrixFile('{}/star.xml'.format(matricesDir), clones)
    )

    memoryLayoutFromFile(memLayout, self.db, clones)

  def numberOfQuantities(self):
    return 9

  def starMatrix(self, dim):
    return self.db.star[dim]

  def addLocal(self, generator):
    super().addLocal(generator)

  def addTime(self, generator):
    qShape = (self.numberOf3DBasisFunctions(), self.numberOfQuantities())
    dQ0 = OptionalDimTensor('dQ(0)', self.Q.optName(), self.Q.optSize(), self.Q.optPos(), qShape, alignStride=True)

    qNodalShape = (self.numberOf2DBasisFunctions(), self.numberOfQuantities())
    dQ0Nodal = OptionalDimTensor('dQNodal(0)',
                                 self.INodal.optName(),
                                 self.INodal.optSize(),
                                 self.INodal.optPos(),
                                 qNodalShape,
                                 alignStride=True)

    power = Scalar('power')
    derivatives = [dQ0]
    generator.add('derivativeTaylorExpansion(0)', self.I['kp'] <= power * dQ0['kp'])
    for i in range(1,self.order):
      derivativeSum = Add()
      for j in range(3):
        derivativeSum += self.db.kDivMT[j][self.t('kl')] * derivatives[-1]['lq'] * self.db.star[j]['qp']
      derivativeSum = DeduceIndices( self.Q['kp'].indices ).visit(derivativeSum)
      derivativeSum = EquivalentSparsityPattern().visit(derivativeSum)
      dQ = OptionalDimTensor('dQ({})'.format(i), self.Q.optName(), self.Q.optSize(), self.Q.optPos(), qShape, spp=derivativeSum.eqspp(), alignStride=True)
      generator.add('derivative({})'.format(i), dQ['kp'] <= derivativeSum)
      generator.add('derivativeTaylorExpansion({})'.format(i), self.I['kp'] <= self.I['kp'] + power * dQ['kp'])
      derivatives.append(dQ)

  def add_include_tensors(self, include_tensors):
    super().add_include_tensors(include_tensors)
