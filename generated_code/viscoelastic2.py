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
  
import argparse
import numpy as np
from yateto import *
from yateto.input import parseXMLMatrixFile, memoryLayoutFromFile
from yateto.gemm_configuration import *
from yateto.ast.node import Add
from yateto.ast.transformer import DeduceIndices, EquivalentSparsityPattern
from yateto.memory import CSCMemoryLayout

from multSim import OptionalDimTensor
import init
import DynamicRupture
import Plasticity
import point

cmdLineParser = argparse.ArgumentParser()
cmdLineParser.add_argument('--matricesDir')
cmdLineParser.add_argument('--outputDir')
cmdLineParser.add_argument('--arch')
cmdLineParser.add_argument('--order')
cmdLineParser.add_argument('--numberOfMechanisms')
cmdLineParser.add_argument('--memLayout')
cmdLineParser.add_argument('--multipleSimulations')
cmdLineParser.add_argument('--dynamicRuptureMethod')
cmdLineParser.add_argument('--PlasticityMethod')
cmdLineArgs = cmdLineParser.parse_args()

arch = useArchitectureIdentifiedBy(cmdLineArgs.arch)

order = int(cmdLineArgs.order)
numberOf2DBasisFunctions = order*(order+1)//2
numberOf3DBasisFunctions = order*(order+1)*(order+2)//6
numberOfQuantities = 9
numberOfAnelasticQuantities = 6
numberOfExtendedQuantities = numberOfQuantities + numberOfAnelasticQuantities
numberOfMechanisms = int(cmdLineArgs.numberOfMechanisms)
multipleSimulations = int(cmdLineArgs.multipleSimulations)

qShape = (numberOf3DBasisFunctions, numberOfQuantities)
qShapeExtended = (numberOf3DBasisFunctions, numberOfExtendedQuantities)
qShapeAnelastic = (numberOf3DBasisFunctions, numberOfAnelasticQuantities, numberOfMechanisms)
  
# Quantities
if multipleSimulations > 1:
  alignStride=set(['fP({})'.format(i) for i in range(3)])
  transpose=True
else:
  alignStride=True
  transpose=False

if transpose:
  t = lambda x: x[::-1]
else:
  t = lambda x: x

clones = {
  'star': ['star(0)', 'star(1)', 'star(2)'],
}
db = parseXMLMatrixFile('{}/matrices_{}.xml'.format(cmdLineArgs.matricesDir, numberOf3DBasisFunctions), transpose=transpose, alignStride=alignStride)
db.update( parseXMLMatrixFile('{}/matrices_viscoelastic.xml'.format(cmdLineArgs.matricesDir, numberOf3DBasisFunctions), clones) )
memoryLayoutFromFile(cmdLineArgs.memLayout, db, clones)

msName = 's'
msPos = 0

Q = OptionalDimTensor('Q', msName, multipleSimulations, msPos, qShape, alignStride=alignStride)
Qext = OptionalDimTensor('Qext', msName, multipleSimulations, msPos, qShapeExtended, alignStride=alignStride)
Qane = OptionalDimTensor('Qane', msName, multipleSimulations, msPos, qShapeAnelastic, alignStride=alignStride)

dQ = [OptionalDimTensor('dQ({})'.format(d), msName, multipleSimulations, msPos, qShape, alignStride=alignStride) for d in range(order)]
dQext = [OptionalDimTensor('dQext({})'.format(d), msName, multipleSimulations, msPos, qShapeExtended, alignStride=alignStride) for d in range(order)]
dQane = [OptionalDimTensor('dQane({})'.format(d), msName, multipleSimulations, msPos, qShapeAnelastic, alignStride=alignStride) for d in range(order)]
I = OptionalDimTensor('I', msName, multipleSimulations, msPos, qShape, alignStride=alignStride)
Iane = OptionalDimTensor('Iane', msName, multipleSimulations, msPos, qShapeAnelastic, alignStride=alignStride)

selectElaSpp = np.zeros((numberOfExtendedQuantities, numberOfQuantities))
selectElaSpp[0:numberOfQuantities,0:numberOfQuantities] = np.eye(numberOfQuantities)
selectEla = Tensor('selectEla', (numberOfExtendedQuantities, numberOfQuantities), selectElaSpp, CSCMemoryLayout)

selectAneSpp = np.zeros((numberOfExtendedQuantities, numberOfAnelasticQuantities))
selectAneSpp[numberOfQuantities:numberOfExtendedQuantities,0:numberOfAnelasticQuantities] = np.eye(numberOfAnelasticQuantities)
selectAne = Tensor('selectAne', (numberOfExtendedQuantities, numberOfAnelasticQuantities), selectAneSpp, CSCMemoryLayout)

E = Tensor('E', (numberOfAnelasticQuantities, numberOfMechanisms, numberOfQuantities))
w = Tensor('w', (numberOfMechanisms,))
W = Tensor('W', (numberOfMechanisms, numberOfMechanisms), np.eye(numberOfMechanisms, dtype=bool), CSCMemoryLayout)

# Kernels
g = Generator(arch)

## Initialization
ti = init.addKernels(g, db, Q, order, numberOfQuantities, numberOfExtendedQuantities)

## Local
volumeSum = Add()
for i in range(3):
  volumeSum += db.kDivM[i][t('kl')] * I['lq'] * db.star[i]['qp']
volumeExt = (Qext['kp'] <= volumeSum)
g.add('volumeExt', volumeExt)

localFluxExt = lambda i: Qext['kp'] <= Qext['kp'] + db.rDivM[i][t('km')] * db.fMrT[i][t('ml')] * I['lq'] * ti.AplusT['qp']
localFluxExtPrefetch = lambda i: I if i == 0 else (Q if i == 1 else None)
g.addFamily('localFluxExt', simpleParameterSpace(4), localFluxExt, localFluxExtPrefetch)

g.add('local', [
  Qane['kpm'] <= Qane['kpm'] + w['m'] * Qext['kq'] * selectAne['qp'] + Iane['kpl'] * W['lm'],
  Q['kp'] <= Qext['kq'] * selectEla['qp'] + Iane['kqm'] * E['qmp']
])

## Neighbour
neighbourFluxExt = lambda h,j,i: Qext['kp'] <= Qext['kp'] + db.rDivM[i][t('km')] * db.fP[h][t('mn')] * db.rT[j][t('nl')] * I['lq'] * ti.AminusT['qp']
neighbourFluxExtPrefetch = lambda h,j,i: I
g.addFamily('neighbourFluxExt', simpleParameterSpace(3,4,4), neighbourFluxExt, neighbourFluxExtPrefetch)

g.add('neighbour', [
  Qane['kpm'] <= Qane['kpm'] + w['m'] * Qext['kq'] * selectAne['qp'],
  Q['kp'] <= Qext['kq'] * selectEla['qp']
])

## ADER
power = Scalar('power')

derivativeTaylorExpansionEla = lambda d: (I['kp'] <= I['kp'] + power * dQ[d]['kp']) if d > 0 else (I['kp'] <= power * dQ[0]['kp'])
derivativeTaylorExpansionAne = lambda d: (Iane['kpm'] <= Iane['kpm'] + power * dQane[d]['kpm']) if d > 0 else (Iane['kpm'] <= power * dQane[0]['kpm'])

def derivative(d):
  derivativeSum = Add()
  for j in range(3):
    derivativeSum += db.kDivMT[j][t('kl')] * dQ[d-1]['lq'] * db.star[j]['qp']
  return derivativeSum

g.addFamily('derivative', parameterSpaceFromRanges(range(1,order)), lambda d: [
  dQext[d]['kp'] <= derivative(d),
  dQ[d]['kp'] <= dQext[d]['kq'] * selectEla['qp'] + dQane[d-1]['kqm'] * E['qmp'],
  dQane[d]['kpm'] <= w['m'] * dQext[d]['kq'] * selectAne['qp'] + dQane[d-1]['kpl'] * W['lm']
])
g.addFamily('derivativeTaylorExpansion', simpleParameterSpace(order), lambda d: [
  derivativeTaylorExpansionEla(d),
  derivativeTaylorExpansionAne(d)
])
g.addFamily('derivativeTaylorExpansionEla', simpleParameterSpace(order), derivativeTaylorExpansionEla)

## Other
DynamicRupture.addKernels(g, Q, Qext, I, alignStride, cmdLineArgs.matricesDir, order, cmdLineArgs.dynamicRuptureMethod, numberOfQuantities, numberOfExtendedQuantities)
Plasticity.addKernels(g, Q, alignStride, cmdLineArgs.matricesDir, order, cmdLineArgs.PlasticityMethod)
point.addKernels(g, Q, ti.oneSimToMultSim, numberOf3DBasisFunctions, numberOfQuantities)

# Generate code
gemmTools = GeneratorCollection([LIBXSMM(arch), PSpaMM(arch)])
g.generate(cmdLineArgs.outputDir, 'seissol', gemmTools)

