#!/usr/bin/env python3
##
# @file
# This file is part of SeisSol.
#
# @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
#
# @section LICENSE
# Copyright (c) 2019, SeisSol Group
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

from abc import ABC, abstractmethod

from yateto import *
from yateto.input import parseXMLMatrixFile
from multSim import OptionalDimTensor

class ADERDGBase(ABC):
  def __init__(self, order, multipleSimulations, matricesDir):
    self.order = order

    self.alignStride = lambda name: True
    if multipleSimulations > 1:
      self.alignStride = lambda name: name.startswith('fP')
    transpose = multipleSimulations > 1
    self.transpose = lambda name: transpose
    self.t = (lambda x: x[::-1]) if transpose else (lambda x: x)

    self.db = parseXMLMatrixFile('{}/matrices_{}.xml'.format(matricesDir, self.numberOf3DBasisFunctions()), transpose=self.transpose, alignStride=self.alignStride)
    clonesQP = {
      'v': [ 'evalAtQP' ],
      'vInv': [ 'projectQP' ]
    }
    self.db.update( parseXMLMatrixFile('{}/plasticity_ip_matrices_{}.xml'.format(matricesDir, order), clonesQP, transpose=self.transpose, alignStride=self.alignStride))

    qShape = (self.numberOf3DBasisFunctions(), self.numberOfQuantities())
    self.Q = OptionalDimTensor('Q', 's', multipleSimulations, 0, qShape, alignStride=True)
    self.I = OptionalDimTensor('I', 's', multipleSimulations, 0, qShape, alignStride=True)

    Ashape = (self.numberOfQuantities(), self.numberOfExtendedQuantities())
    self.AplusT = Tensor('AplusT', Ashape)
    self.AminusT = Tensor('AminusT', Ashape)
    Tshape = (self.numberOfExtendedQuantities(), self.numberOfExtendedQuantities())
    self.T = Tensor('T', Tshape)
    QgodShape = (self.numberOfQuantities(), self.numberOfQuantities())
    self.Tinv = Tensor('Tinv', QgodShape)
    self.QgodLocal = Tensor('QgodLocal', QgodShape)
    self.QgodNeighbor = Tensor('QgodNeighbor', QgodShape)

    self.oneSimToMultSim = Tensor('oneSimToMultSim', (self.Q.optSize(),), spp={(i,): '1.0' for i in range(self.Q.optSize())})

  def numberOf2DBasisFunctions(self):
    return self.order*(self.order+1)//2

  def numberOf3DBasisFunctions(self):
    return self.order*(self.order+1)*(self.order+2)//6

  def numberOf3DQuadraturePoints(self):
    return (self.order+1)**3

  @abstractmethod
  def numberOfQuantities(self):
    pass

  @abstractmethod
  def numberOfExtendedQuantities(self):
    pass

  @abstractmethod
  def extendedQTensor(self):
    pass

  @abstractmethod
  def starMatrix(self, dim):
    pass

  def addInit(self, generator):
    fluxScale = Scalar('fluxScale')
    computeFluxSolverLocal = self.AplusT['ij'] <= fluxScale * self.Tinv['ki'] * self.QgodLocal['kq'] * self.db.star[0]['ql'] * self.T['jl']
    generator.add('computeFluxSolverLocal', computeFluxSolverLocal)

    computeFluxSolverNeighbor = self.AminusT['ij'] <= fluxScale * self.Tinv['ki'] * self.QgodNeighbor['kq'] * self.db.star[0]['ql'] * self.T['jl']
    generator.add('computeFluxSolverNeighbor', computeFluxSolverNeighbor)

    QFortran = Tensor('QFortran', (self.numberOf3DBasisFunctions(), self.numberOfQuantities()))
    multSimToFirstSim = Tensor('multSimToFirstSim', (self.Q.optSize(),), spp={(0,): '1.0'})
    if self.Q.hasOptDim():
      copyQToQFortran = QFortran['kp'] <= self.Q['kp'] * multSimToFirstSim['s']
    else:
      copyQToQFortran = QFortran['kp'] <= self.Q['kp']

    generator.add('copyQToQFortran', copyQToQFortran)

  @abstractmethod
  def addLocal(self, generator):
    pass

  @abstractmethod
  def addNeighbor(self, generator):
    pass

  @abstractmethod
  def addTime(self, generator):
    pass
