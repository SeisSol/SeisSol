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
numberOf3DQuadraturePoints = (order+1)**3
numberOfQuantities = 9
multipleSimulations = int(cmdLineArgs.multipleSimulations)

qShape = (numberOf3DBasisFunctions, numberOfQuantities)

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
db.update( parseXMLMatrixFile('{}/star.xml'.format(cmdLineArgs.matricesDir, numberOf3DBasisFunctions), clones) )
clonesQP = {
  'v': [ 'evalAtQP' ],
  'vInv': [ 'projectQP' ]
}
db.update( parseXMLMatrixFile('{}/plasticity_ip_matrices_{}.xml'.format(cmdLineArgs.matricesDir, order), clonesQP, transpose=transpose))
memoryLayoutFromFile(cmdLineArgs.memLayout, db, clones)

Q = OptionalDimTensor('Q', 's', multipleSimulations, 0, qShape, alignStride=alignStride)
dQ0 = OptionalDimTensor('dQ(0)', 's', multipleSimulations, 0, qShape, alignStride=alignStride)
I = OptionalDimTensor('I', 's', multipleSimulations, 0, qShape, alignStride=alignStride)
iniCond = OptionalDimTensor('iniCond', 's', multipleSimulations, 0, (numberOf3DQuadraturePoints, numberOfQuantities))
dofsQP = OptionalDimTensor('dofsQP', 's', multipleSimulations, 0, (numberOf3DQuadraturePoints, numberOfQuantities))


# Kernels
g = Generator(arch)

## Initialization
ti = init.addKernels(g, db, Q, order, numberOfQuantities, numberOfQuantities)
g.add('projectIniCond', Q['kp'] <= db.projectQP[t('kl')] * iniCond['lp'])
g.add('evalAtQP', dofsQP['kp'] <= db.evalAtQP[t('kl')] * Q['lp'])


## Local
volumeSum = Q['kp']
for i in range(3):
  volumeSum += db.kDivM[i][t('kl')] * I['lq'] * db.star[i]['qp']
volume = (Q['kp'] <= volumeSum)
g.add('volume', volume)

localFlux = lambda i: Q['kp'] <= Q['kp'] + db.rDivM[i][t('km')] * db.fMrT[i][t('ml')] * I['lq'] * ti.AplusT['qp']
localFluxPrefetch = lambda i: I if i == 0 else (Q if i == 1 else None)
g.addFamily('localFlux', simpleParameterSpace(4), localFlux, localFluxPrefetch)

## Neighbour
neighbourFlux = lambda h,j,i: Q['kp'] <= Q['kp'] + db.rDivM[i][t('km')] * db.fP[h][t('mn')] * db.rT[j][t('nl')] * I['lq'] * ti.AminusT['qp']
neighbourFluxPrefetch = lambda h,j,i: I
g.addFamily('neighboringFlux', simpleParameterSpace(3,4,4), neighbourFlux, neighbourFluxPrefetch)

## ADER
power = Scalar('power')
derivatives = [dQ0]
g.add('derivativeTaylorExpansion(0)', I['kp'] <= power * dQ0['kp'])
for i in range(1,order):
  derivativeSum = Add()
  for j in range(3):
    derivativeSum += db.kDivMT[j][t('kl')] * derivatives[-1]['lq'] * db.star[j]['qp']
  derivativeSum = DeduceIndices( Q['kp'].indices ).visit(derivativeSum)
  derivativeSum = EquivalentSparsityPattern().visit(derivativeSum)
  dQ = OptionalDimTensor('dQ({})'.format(i), 's', multipleSimulations, 0, qShape, spp=derivativeSum.eqspp(), alignStride=True)
  g.add('derivative({})'.format(i), dQ['kp'] <= derivativeSum)
  g.add('derivativeTaylorExpansion({})'.format(i), I['kp'] <= I['kp'] + power * dQ['kp'])
  derivatives.append(dQ)

## Other
DynamicRupture.addKernels(g, db, ti, Q, Q, I, alignStride, cmdLineArgs.matricesDir, order, cmdLineArgs.dynamicRuptureMethod, numberOfQuantities, numberOfQuantities)
Plasticity.addKernels(g, Q, alignStride, cmdLineArgs.matricesDir, order, cmdLineArgs.PlasticityMethod)
point.addKernels(g, Q, ti.oneSimToMultSim, numberOf3DBasisFunctions, numberOfQuantities)

# Generate code
gemmTools = GeneratorCollection([LIBXSMM(arch), PSpaMM(arch)])
g.generate(cmdLineArgs.outputDir, 'seissol', gemmTools)

