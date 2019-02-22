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
from yateto import *
from yateto.input import parseXMLMatrixFile, memoryLayoutFromFile
from yateto.gemm_configuration import *
from yateto.ast.node import Add
from yateto.ast.transformer import DeduceIndices, EquivalentSparsityPattern

import DynamicRupture
import Plasticity

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
multipleSimulations = int(cmdLineArgs.multipleSimulations)

# Quantities
if multipleSimulations > 1:
  qShape = (multipleSimulations, numberOf3DBasisFunctions, numberOfQuantities)
  qi = lambda x: 's' + x
  alignStride=set(['fP({})'.format(i) for i in range(3)])
  transpose=True
else:
  qShape = (numberOf3DBasisFunctions, numberOfQuantities)
  qi = lambda x: x
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
memoryLayoutFromFile(cmdLineArgs.memLayout, db, clones)

Q = Tensor('Q', qShape, alignStride=alignStride)
QFortran = Tensor('QFortran', (numberOf3DBasisFunctions, numberOfQuantities))
dQ0 = Tensor('dQ(0)', qShape, alignStride=alignStride)
I = Tensor('I', qShape, alignStride=alignStride)

# Flux solver
AplusT = Tensor('AplusT', (numberOfQuantities, numberOfQuantities))
AminusT = Tensor('AminusT', (numberOfQuantities, numberOfQuantities))
T = Tensor('T', (numberOfQuantities, numberOfQuantities))
Tinv = Tensor('Tinv', (numberOfQuantities, numberOfQuantities))
QgodLocal = Tensor('QgodLocal', (numberOfQuantities, numberOfQuantities))
QgodNeighbor = Tensor('QgodNeighbor', (numberOfQuantities, numberOfQuantities))

# Kernels
g = Generator(arch)

## Main kernels
volumeSum = Q[qi('kp')]
for i in range(3):
  volumeSum += db.kDivM[i][t('kl')] * I[qi('lq')] * db.star[i]['qp']
volume = (Q[qi('kp')] <= volumeSum)
g.add('volume', volume)

localFlux = lambda i: Q[qi('kp')] <= Q[qi('kp')] + db.rDivM[i][t('km')] * db.fMrT[i][t('ml')] * I[qi('lq')] * AplusT['qp']
localFluxPrefetch = lambda i: I if i == 0 else (Q if i == 1 else None)
g.addFamily('localFlux', simpleParameterSpace(4), localFlux, localFluxPrefetch)

neighbourFlux = lambda h,j,i: Q[qi('kp')] <= Q[qi('kp')] + db.rDivM[i][t('km')] * db.fP[h][t('mn')] * db.rT[j][t('nl')] * I[qi('lq')] * AminusT['qp']
neighbourFluxPrefetch = lambda h,j,i: I
g.addFamily('neighboringFlux', simpleParameterSpace(3,4,4), neighbourFlux, neighbourFluxPrefetch)

power = Scalar('power')
derivatives = [dQ0]
g.add('derivativeTaylorExpansion(0)', I[qi('kp')] <= power * dQ0[qi('kp')])
for i in range(1,order):
  derivativeSum = Add()
  for j in range(3):
    derivativeSum += db.kDivMT[j][t('kl')] * derivatives[-1][qi('lq')] * db.star[j]['qp']
  derivativeSum = DeduceIndices( Q[qi('kp')].indices ).visit(derivativeSum)
  derivativeSum = EquivalentSparsityPattern().visit(derivativeSum)
  dQ = Tensor('dQ({})'.format(i), qShape, spp=derivativeSum.eqspp(), alignStride=True)
  g.add('derivative({})'.format(i), dQ[qi('kp')] <= derivativeSum)
  g.add('derivativeTaylorExpansion({})'.format(i), I[qi('kp')] <= I[qi('kp')] + power * dQ[qi('kp')])
  derivatives.append(dQ)

## Initialization kernels
fluxScale = Scalar('fluxScale')
computeFluxSolverLocal = AplusT['ij'] <= fluxScale * Tinv['ki'] * QgodLocal['kq'] * db.star[0]['ql'] * T['jl']
g.add('computeFluxSolverLocal', computeFluxSolverLocal)

computeFluxSolverNeighbor = AminusT['ij'] <= fluxScale * Tinv['ki'] * QgodNeighbor['kq'] * db.star[0]['ql'] * T['jl']
g.add('computeFluxSolverNeighbor', computeFluxSolverNeighbor)

if multipleSimulations > 1:
  oneSimToMultSim = Tensor('oneSimToMultSim', (multipleSimulations,), spp={(i,): '1.0' for i in range(multipleSimulations)})
  addQFortran = Q[qi('kp')] <= Q[qi('kp')] + QFortran['kp'] * oneSimToMultSim['s']

  multSimToFirstSim = Tensor('multSimToFirstSim', (multipleSimulations,), spp={(0,): '1.0'})
  copyQToQFortran = QFortran['kp'] <= Q[qi('kp')] * multSimToFirstSim['s']

  copyDQToDQFortran = lambda degree: QFortran['kp'] <= derivatives[degree][qi('kp')] * multSimToFirstSim['s']
  g.addFamily('copyDQToDQFortran', simpleParameterSpace(order), copyDQToDQFortran)
else:
  addQFortran = Q[qi('kp')] <= Q[qi('kp')] + QFortran['kp']
  copyQToQFortran = QFortran['kp'] <= Q[qi('kp')]

  copyDQToDQFortran = lambda degree: QFortran['kp'] <= derivatives[degree][qi('kp')]
  g.addFamily('copyDQToDQFortran', simpleParameterSpace(order), copyDQToDQFortran)

g.add('addQFortran', addQFortran)
g.add('copyQToQFortran', copyQToQFortran)

DynamicRupture.addKernels(g, Q, I, qi, qShape, alignStride, cmdLineArgs.matricesDir, order, cmdLineArgs.dynamicRuptureMethod, numberOfQuantities, numberOfQuantities)
Plasticity.addKernels(g, qi, qShape, alignStride, cmdLineArgs.matricesDir, order, cmdLineArgs.PlasticityMethod)

# Generate code
gemmTools = GeneratorCollection([LIBXSMM(arch), PSpaMM(arch)])
g.generate(cmdLineArgs.outputDir, 'seissol', gemmTools)

