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
    # todo change path
    self.db.update(
      parseJSONMatrixFile('{}/nodal/nodalBoundary_matrices_{}.json'.format(matricesDir,
                                                                           self.order - 1),
                          {},
                          alignStride=self.alignStride,
                          transpose=self.transpose)
    )

    memoryLayoutFromFile(memLayout, self.db, clones)

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

    self.INodal = OptionalDimTensor('INodal',
                                    's',
                                    False, #multipleSimulations,
                                    0,
                                    (self.numberOf2DBasisFunctions(), self.numberOfQuantities()),
                                    alignStride=True)

    # Add nodal 2d weights, needed in this generator
    self.nodes2D = Tensor('nodes2DTmp', (self.numberOf2DBasisFunctions(),2)) 

  def addLocal(self, generator):
    volumeSum = self.Q['kp']
    for i in range(3):
      volumeSum += self.db.kDivM[i][self.t('kl')] * self.I['lq'] * self.db.star[i]['qp']
    volume = (self.Q['kp'] <= volumeSum)
    generator.add('volume', volume)

    localFlux = lambda i: self.Q['kp'] <= self.Q['kp'] + self.db.rDivM[i][self.t('km')] * self.db.fMrT[i][self.t('ml')] * self.I['lq'] * self.AplusT['qp']
    localFluxPrefetch = lambda i: self.I if i == 0 else (self.Q if i == 1 else None)
    generator.addFamily('localFlux', simpleParameterSpace(4), localFlux, localFluxPrefetch)

    # TODO(Lukas) rDivM[i] * V2nTo2m can be simplified
    localFluxNodal = lambda i: self.Q['kp'] <= self.Q['kp'] + self.db.rDivM[i][self.t('km')] * self.db.V2nTo2m[self.t('mn')] * self.INodal['no'] * self.AplusT['op']
    localFluxNodalPrefetch = localFluxPrefetch
    generator.addFamily('localFluxNodal', simpleParameterSpace(4), localFluxNodal, localFluxNodalPrefetch)

  def addNeighbor(self, generator):
    neighbourFlux = lambda h,j,i: self.Q['kp'] <= self.Q['kp'] + self.db.rDivM[i][self.t('km')] * self.db.fP[h][self.t('mn')] * self.db.rT[j][self.t('nl')] * self.I['lq'] * self.AminusT['qp']
    neighbourFluxPrefetch = lambda h,j,i: self.I
    generator.addFamily('neighboringFlux', simpleParameterSpace(3,4,4), neighbourFlux, neighbourFluxPrefetch)

    # todo: transpose? not local flux?
    # neighbourFluxNodal = lambda h,i: self.Q['kp'] <= self.Q['kp'] + self.db.rDivM[i][self.t('km')] * self.db.fP[h][self.t('mn')] * self.db.V2nTo2m[self.t('nl')] * self.INodal['lq'] * self.AminusT['qp']
    # neighbourFluxNodalPrefetch = lambda h,i: None #self.I ?
    # generator.addFamily('neighboringFluxNodal',
    #                     simpleParameterSpace(3,4),
    #                     neighbourFluxNodal,
    #                     neighbourFluxNodalPrefetch)

    projectToNodalBoundary = lambda j: self.INodal['kp'] <= self.db.V3mTo2nFace[j]['km'] * self.I['mp']
    generator.addFamily('projectToNodalBoundary',
                        simpleParameterSpace(4),
                        projectToNodalBoundary)

    # todo maybe integrate this better in the flux 
    # todo tinv? t transpose?
    projectToNodalBoundaryRotated = lambda j: self.INodal['kp'] <= self.db.V3mTo2nFace[j]['kl'] \
                                    * self.I['lm'] \
                                    * self.T['mp']

    generator.addFamily('projectToNodalBoundaryRotated',
                        simpleParameterSpace(4),
                        projectToNodalBoundaryRotated)

    selectZDisplacement = np.zeros((self.numberOfQuantities(), 1))
    selectZDisplacement[8,0] = 1 # todo is correct?
    selectZDisplacement = Tensor('selectZDisplacement',
                            selectZDisplacement.shape,
                                 selectZDisplacement,
                                 CSCMemoryLayout)

    selectZDisplacementFromDisplacements = np.zeros((3, 1))
    selectZDisplacementFromDisplacements[2,0] = 1 # todo is correct?
    selectZDisplacementFromDisplacements = Tensor('selectZDisplacementFromDisplacements',
                                                   selectZDisplacementFromDisplacements.shape,
                                                   selectZDisplacementFromDisplacements,
                                                   CSCMemoryLayout)

    # displacement = OptionalDimTensor('displacement',
    #                                  aderdg.Q.optName(),
    #                                  aderdg.Q.optSize(),
    #                                  aderdg.Q.optPos(),
    #                                  (numberOf3DBasisFunctions, 3), alignStride=True)
    
    self.INodalDisplacement = OptionalDimTensor('INodalDisplacement',
                                                self.Q.optName(),
                                                self.Q.optSize(),
                                                self.Q.optPos(),
                                                (self.numberOf2DBasisFunctions(), 1),
                                                alignStride=True)

    #copied from SurfaceDisplacement.py
    displacement = OptionalDimTensor('displacement',
                                     self.Q.optName(),
                                     self.Q.optSize(),
                                     self.Q.optPos(),
                                     (self.numberOf3DBasisFunctions(), 3), alignStride=True)

    dt = Scalar('dt')
    displacementAvgNodal = lambda side: self.INodalDisplacement['ip'] <= self.db.V3mTo2nFace[side]['ij'] * self.I['jk'] * selectZDisplacement['kp'] \
        + dt * self.db.V3mTo2nFace[side]['ij'] * displacement['jk'] * selectZDisplacementFromDisplacements['kp']
    generator.addFamily('displacementAvgNodal',
                        simpleParameterSpace(4),
                        displacementAvgNodal)

    nodes2DGen = self.nodes2D['ij'] <= self.db.nodes2D['ij']
    generator.add('nodes2DGen', nodes2DGen)

  def addTime(self, generator):
    qShape = (self.numberOf3DBasisFunctions(), self.numberOfQuantities())
    dQ0 = OptionalDimTensor('dQ(0)', self.Q.optName(), self.Q.optSize(), self.Q.optPos(), qShape, alignStride=True)
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
