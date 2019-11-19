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
from yateto import Tensor, Scalar, simpleParameterSpace
from yateto.util import tensor_collection_from_constant_expression
from yateto.input import parseXMLMatrixFile, parseJSONMatrixFile, memoryLayoutFromFile
from yateto.memory import CSCMemoryLayout
from yateto.ast.node import Add
from yateto.ast.transformer import DeduceIndices, EquivalentSparsityPattern

from aderdg import ADERDGBase
import aderdg
from multSim import OptionalDimTensor

class ADERDG(ADERDGBase):
  def __init__(self, order, multipleSimulations, matricesDir, memLayout):
    super().__init__(order, multipleSimulations, matricesDir)
    clones = {
      'star': ['star(0)', 'star(1)', 'star(2)'],
    }
    self.db.update(
      parseXMLMatrixFile('{}/star.xml'.format(matricesDir), clones)
    )

    memoryLayoutFromFile(memLayout, self.db, clones)
    self.INodal = OptionalDimTensor('INodal',
                                    's',
                                    False, #multipleSimulations,
                                    0,
                                    (self.numberOf2DBasisFunctions(), self.numberOfQuantities()),
                                    alignStride=True)



  def numberOfQuantities(self):
    return 9

  def numberOfExtendedQuantities(self):
    return self.numberOfQuantities()

  def extendedQTensor(self):
    return self.Q

  def starMatrix(self, dim):
    return self.db.star[dim]

  def addInit(self, generator):
    super().addInit(generator)

    iniShape = (self.numberOf3DQuadraturePoints(), self.numberOfQuantities())
    iniCond = OptionalDimTensor('iniCond', self.Q.optName(), self.Q.optSize(), self.Q.optPos(), iniShape, alignStride=True)
    dofsQP = OptionalDimTensor('dofsQP', self.Q.optName(), self.Q.optSize(), self.Q.optPos(), iniShape, alignStride=True)

    generator.add('projectIniCond', self.Q['kp'] <= self.db.projectQP[self.t('kl')] * iniCond['lp'])
    generator.add('evalAtQP', dofsQP['kp'] <= self.db.evalAtQP[self.t('kl')] * self.Q['lp'])

  def addLocal(self, generator):
    volumeSum = self.Q['kp']
    for i in range(3):
      volumeSum += self.db.kDivM[i][self.t('kl')] * self.I['lq'] * self.db.star[i]['qp']
    volume = (self.Q['kp'] <= volumeSum)
    generator.add('volume', volume)

    localFlux = lambda i: self.Q['kp'] <= self.Q['kp'] + self.db.rDivM[i][self.t('km')] * self.db.fMrT[i][self.t('ml')] * self.I['lq'] * self.AplusT['qp']
    localFluxPrefetch = lambda i: self.I if i == 0 else (self.Q if i == 1 else None)
    generator.addFamily('localFlux', simpleParameterSpace(4), localFlux, localFluxPrefetch)

    rDivM_mult_V2nTo2m = tensor_collection_from_constant_expression(
      base_name='rDivMMultV2nTo2m',
      expressions=lambda i: self.db.rDivM[i]['jk'] * self.db.V2nTo2m['kl'],
      group_indices=range(4),
      target_indices='jl')
    self.db.update(rDivM_mult_V2nTo2m)

    localFluxNodal = lambda i: self.Q['kp'] <= self.Q['kp'] + self.db.rDivMMultV2nTo2m[i]['kn'] * self.INodal['no'] * self.AminusT['op']
    localFluxNodalPrefetch = localFluxPrefetch
    generator.addFamily('localFluxNodal', simpleParameterSpace(4), localFluxNodal, localFluxNodalPrefetch)

    easi_ident_map = np.stack([np.eye(self.numberOfQuantities())] * self.numberOf2DBasisFunctions(), axis=2)
    assert(easi_ident_map.shape ==
           (self.numberOfQuantities(), self.numberOfQuantities(), self.numberOf2DBasisFunctions()))
    easi_ident_map = Tensor('easiIdentMap',
                            easi_ident_map.shape,
                            easi_ident_map,
                            alignStride=False)
    easi_boundary_constant = Tensor('easiBoundaryConstant',
                                    (self.numberOfQuantities(), self.numberOf2DBasisFunctions()),
                                    alignStride=False)
    easi_boundary_map = Tensor('easiBoundaryMap',
                               (self.numberOfQuantities(), self.numberOfQuantities(), self.numberOf2DBasisFunctions(),),
                               alignStride=False)
    create_easi_boundary_ghost_cells = (
            self.INodal['la'] <= easi_boundary_map['abl'] * self.INodal['lb'] + easi_ident_map['abl'] * easi_boundary_constant['bl']
    )
    generator.add('createEasiBoundaryGhostCells', create_easi_boundary_ghost_cells)

  def addNeighbor(self, generator):
    neighbourFlux = lambda h,j,i: self.Q['kp'] <= self.Q['kp'] + self.db.rDivM[i][self.t('km')] * self.db.fP[h][self.t('mn')] * self.db.rT[j][self.t('nl')] * self.I['lq'] * self.AminusT['qp']
    neighbourFluxPrefetch = lambda h,j,i: self.I
    generator.addFamily('neighboringFlux', simpleParameterSpace(3,4,4), neighbourFlux, neighbourFluxPrefetch)

    projectToNodalBoundary = lambda j: self.INodal['kp'] <= self.db.V3mTo2nFace[j]['km'] * self.I['mp']
    generator.addFamily('projectToNodalBoundary',
                        simpleParameterSpace(4),
                        projectToNodalBoundary)

    projectToNodalBoundaryRotated = lambda j: self.INodal['kp'] <= self.db.V3mTo2nFace[j]['kl'] \
                                              * self.I['lm'] \
                                              * self.T['mp']

    generator.addFamily('projectToNodalBoundaryRotated',
                        simpleParameterSpace(4),
                        projectToNodalBoundaryRotated)

    rotateBoundaryDofsBack = self.INodal['kp'] <= self.INodal['kl'] * self.Tinv['lp']
    generator.add('rotateBoundaryDofsBack', rotateBoundaryDofsBack)
    
    selectZDisplacement = np.zeros((self.numberOfQuantities(), 1))
    selectZDisplacement[8, 0] = 1
    selectZDisplacement = Tensor('selectZDisplacement',
                                 selectZDisplacement.shape,
                                 selectZDisplacement,
                                 CSCMemoryLayout)

    selectZDisplacementFromDisplacements = np.zeros((3, 1))
    selectZDisplacementFromDisplacements[2, 0] = 1
    selectZDisplacementFromDisplacements = Tensor('selectZDisplacementFromDisplacements',
                                                   selectZDisplacementFromDisplacements.shape,
                                                   selectZDisplacementFromDisplacements,
                                                   CSCMemoryLayout)

    self.INodalDisplacement = OptionalDimTensor('INodalDisplacement',
                                                self.Q.optName(),
                                                self.Q.optSize(),
                                                self.Q.optPos(),
                                                (self.numberOf2DBasisFunctions(), 1),
                                                alignStride=True)

    displacement = OptionalDimTensor('displacement',
                                     self.Q.optName(),
                                     self.Q.optSize(),
                                     self.Q.optPos(),
                                     (self.numberOf3DBasisFunctions(), 3),
                                     alignStride=True)

    dt = Scalar('dt')
    displacementAvgNodal = lambda side: self.INodalDisplacement['ip'] <= self.db.V3mTo2nFace[side]['ij'] * self.I['jk'] * selectZDisplacement['kp'] \
        + dt * self.db.V3mTo2nFace[side]['ij'] * displacement['jk'] * selectZDisplacementFromDisplacements['kp']

    generator.addFamily('displacementAvgNodal',
                        simpleParameterSpace(4),
                        displacementAvgNodal)

    self.INodalUpdate = OptionalDimTensor('INodalUpdate',
                                          self.INodal.optName(),
                                          self.INodal.optSize(),
                                          self.INodal.optPos(),
                                          (self.numberOf2DBasisFunctions(), self.numberOfQuantities()),
                                          alignStride=True)

    factor = Scalar('factor')
    updateINodal = self.INodal['kp'] <= self.INodal['kp'] + factor * self.INodalUpdate['kp']
    generator.add('updateINodal', updateINodal)


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
    include_tensors.add(self.db.V2nTo2m)
