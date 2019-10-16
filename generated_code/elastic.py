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

import operator
from functools import reduce
import numpy as np

from yateto import Tensor, Scalar
from yateto.input import parseXMLMatrixFile, parseJSONMatrixFile, memoryLayoutFromFile
from yateto.ast.transformer import DeduceIndices, EquivalentSparsityPattern

from aderdg import ADERDGStandard
from multSim import OptionalDimTensor


def choose(n, k):
  num = reduce(operator.mul, range(n, n-k, -1), 1)
  denom = reduce(operator.mul, range(1, k+1), 1)
  return num // denom

class ADERDG(ADERDGStandard):
  def __init__(self, order, multipleSimulations, matricesDir, memLayout):
    super().__init__(order, multipleSimulations, matricesDir)
    clones = {
      'star': ['star(0)', 'star(1)', 'star(2)'],
    }
    self.db.update( parseXMLMatrixFile('{}/star.xml'.format(matricesDir), clones) )
    self.db.update( parseJSONMatrixFile('{}/stp_{}.json'.format(matricesDir, order)) )
    memoryLayoutFromFile(memLayout, self.db, clones)

  def numberOfQuantities(self):
    return 9

  def numberOfExtendedQuantities(self):
    return self.numberOfQuantities()

  def starMatrix(self, dim):
    return self.db.star[dim]

  def addTime(self, generator):
    super().addTime(generator)

    stiffnessValues = [self.db.kDivMT[d].values_as_ndarray() for d in range(3)]

    qShape = (self.numberOf3DBasisFunctions(), self.numberOfQuantities())
    stpShape = (self.numberOf3DBasisFunctions(), self.numberOfQuantities(), self.order)
    stpRhs = OptionalDimTensor('stpRhs', self.Q.optName(), self.Q.optSize(), self.Q.optPos(), stpShape, alignStride=True)
    stp = OptionalDimTensor('stp', self.Q.optName(), self.Q.optSize(), self.Q.optPos(), stpShape, alignStride=True)
    timestep = Scalar('timestep')

    kernels = [stpRhs['kpt'] <= self.Q['kp'] * self.db.wHat['t']]

    for n in range(self.order-1,-1,-1):
      Bn_1 = choose(n-1+3,3)
      Bn = choose(n+3,3)
      fullShape = (self.numberOf3DBasisFunctions(), self.numberOf3DBasisFunctions())
      selectModesSpp = np.zeros(fullShape)
      selectModesSpp[Bn_1:Bn,Bn_1:Bn] = np.eye(Bn-Bn_1) 
      selectModes = Tensor('selectModes({})'.format(n), fullShape, spp=selectModesSpp)
      kSub = [None] * 3
      if n > 0: 
        for d in range(3):
          stiffnessSpp = np.zeros(fullShape)
          stiffnessSpp[:,Bn_1:Bn] = -stiffnessValues[d][:,Bn_1:Bn]
          kSub[d] = Tensor('kDivMTSub({},{})'.format(d,n-1), fullShape, spp=stiffnessSpp)

      kernel = stp['kpt'] <= stp['kpt'] + selectModes['kl'] * stpRhs['lpu'] * self.db.Zinv['tu']
      kernels.append(kernel)
      if n > 0:
        derivativeSum = stpRhs['kpt']
        for d in range(3):
          derivativeSum += timestep * kSub[d]['kl'] * stp['lqt'] * self.starMatrix(d)['qp']
        updateRhs = stpRhs['kpt'] <= derivativeSum
        kernels.append(updateRhs)

    timeInt = self.I['kp'] <= timestep * stp['kpt'] * self.db.timeInt['t']
    kernels.append(timeInt)

    generator.add('stp', kernels)

