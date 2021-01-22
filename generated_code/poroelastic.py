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
import numpy as np
from functools import reduce

from yateto import Tensor, Scalar
from yateto.input import parseXMLMatrixFile, parseJSONMatrixFile, memoryLayoutFromFile
from yateto.ast.transformer import DeduceIndices, EquivalentSparsityPattern

from aderdg import LinearADERDG
from multSim import OptionalDimTensor

def choose(n, k):
  num = reduce(operator.mul, range(n, n-k, -1), 1)
  denom = reduce(operator.mul, range(1, k+1), 1)
  return num // denom

class PoroelasticADERDG(LinearADERDG):
  def __init__(self, order, multipleSimulations, matricesDir, memLayout, numberOfMechanisms, **kwargs):

    super().__init__(order, multipleSimulations, matricesDir)
    clones = {
      'star': ['star(0)', 'star(1)', 'star(2)'],
    }
    self.db.update( parseXMLMatrixFile('{}/matrices_poroelastic.xml'.format(matricesDir), clones) )
    self.db.update( parseJSONMatrixFile('{}/stp_{}.json'.format(matricesDir, order), clones) )

    memoryLayoutFromFile(memLayout, self.db, clones)

  def numberOfQuantities(self):
    return 13 

  def numberOfExtendedQuantities(self):
    return self.numberOfQuantities()

  def starMatrix(self, dim):
    return self.db.star[dim]

  def sourceMatrix(self):
    return self.db.ET

  def transformation_spp(self):
    spp = np.zeros((self.numberOfQuantities(), self.numberOfQuantities()), dtype=bool)
    spp[0:6,0:6] = 1
    spp[6:9,6:9] = 1
    spp[9,9] = 1
    spp[10:13,10:13] = 1
    return spp

  def transformation_inv_spp(self):
    return self.transformation_spp()

  def addTime(self, generator, targets):
    super().addTime(generator, targets)

    stiffnessValues = [self.db.kDivMT[d].values_as_ndarray() for d in range(3)]
    fullShape = (self.numberOf3DBasisFunctions(), self.numberOf3DBasisFunctions())

    qShape = (self.numberOf3DBasisFunctions(), self.numberOfQuantities())
    stpShape = (self.numberOf3DBasisFunctions(), self.numberOfQuantities(), self.order)
    quantityShape = (self.numberOfQuantities(), self.numberOfQuantities())
    stpRhs = OptionalDimTensor('stpRhs', self.Q.optName(), self.Q.optSize(), self.Q.optPos(), stpShape, alignStride=True)
    stp = OptionalDimTensor('stp', self.Q.optName(), self.Q.optSize(), self.Q.optPos(), stpShape, alignStride=True)
    G = Tensor('G', quantityShape, spp=self.db.ET.spp().as_ndarray())
    Zinv = Tensor('Zinv', (self.numberOfQuantities(), self.order, self.order))
    timestep = Scalar('timestep')

    def modeRange(n):
      Bn_1 = choose(n-1+3,3)
      Bn = choose(n+3,3)
      return (Bn_1,Bn)

    def selectModes(n):
      Bn_1, Bn = modeRange(n)
      selectModesSpp = np.zeros(fullShape)
      selectModesSpp[Bn_1:Bn,Bn_1:Bn] = np.eye(Bn-Bn_1) 
      return Tensor('selectModes({})'.format(n), fullShape, spp=selectModesSpp)

    def selectQuantity(o):
      selectSpp = np.zeros((self.numberOfQuantities(),))
      selectSpp[o] = 1
      return Tensor('selectQuantity({})'.format(o), selectSpp.shape, spp = selectSpp)

    def selectQuantity_G(o):
      selectSpp = np.zeros((self.numberOfQuantities(),))
      selectSpp[o-4] = 1
      return Tensor('selectQuantity_G({})'.format(o), selectSpp.shape, spp = selectSpp)

    def Zinv(o):
        return Tensor('Zinv({})'.format(o), (self.order, self.order))

    def selectQuantity_Z(o):
      selectSpp = np.zeros((self.numberOfQuantities(),))
      selectSpp[o] = 1
      return Tensor('selectQuantity_Z({})'.format(o), selectSpp.shape, spp = selectSpp)

    def kSub(d,n):
      Bn_1, Bn = modeRange(n)
      stiffnessSpp = np.zeros(fullShape)
      stiffnessSpp[:,Bn_1:Bn] = -stiffnessValues[d][:,Bn_1:Bn]
      return Tensor('kDivMTSub({},{})'.format(d,n), fullShape, spp=stiffnessSpp)

    kernels = list()

    kernels.append( stpRhs['kpt'] <= self.Q['kp'] * self.db.wHat['t'] )
    for n in range(self.order-1,-1,-1):
      for o in range(self.numberOfQuantities()-1,-1,-1):
        kernels.append( stp['kpt'] <= stp['kpt'] + selectModes(n)['kl'] * selectQuantity(o)['p'] * stpRhs['lpu'] * Zinv(o)['tu'] )
        if o >= 10:
          kernels.append( stpRhs['kpt'] <= stpRhs['kpt'] + timestep * selectQuantity_G(o)['p'] * selectQuantity(o)['u'] * G['up'] * selectModes(n)['kl'] * stp['lut'] )
      if n > 0:
        derivativeSum = stpRhs['kpt']
        for d in range(3):
          derivativeSum += timestep * kSub(d,n)['kl'] * stp['lqt'] * self.starMatrix(d)['qp']
      kernels.append( stpRhs['kpt'] <=  derivativeSum )
    kernels.append( self.I['kp'] <= timestep * stp['kpt'] * self.db.timeInt['t'] )

    generator.add('stp', kernels)

  def add_include_tensors(self, include_tensors):
    super().add_include_tensors(include_tensors)
    include_tensors.add(self.db.Z)
