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
from yateto.input import parseXMLMatrixFile
from yateto.ast.node import Add
from yateto.ast.transformer import DeduceIndices, EquivalentSparsityPattern

#import DynamicRupture

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
multipleSimulations = bool(cmdLineArgs.multipleSimulations)

# Quantities
if multipleSimulations:
  qShape = (arch.alignedReals, numberOf3DBasisFunctions, numberOfQuantities)
  qi = lambda x: 's' + x
  alignStride=False
  transpose=True
  t = lambda x: x[::-1]
else:
  qShape = (numberOf3DBasisFunctions, numberOfQuantities)
  qi = lambda x: x
  alignStride=True
  transpose=False
  t = lambda x: x

clones = {
  'star': ['star[0]', 'star[1]', 'star[2]'],
}
db = parseXMLMatrixFile('{}/matrices_{}.xml'.format(cmdLineArgs.matricesDir, numberOf3DBasisFunctions), transpose=transpose, alignStride=alignStride)
db.update( parseXMLMatrixFile('{}/star.xml'.format(cmdLineArgs.matricesDir, numberOf3DBasisFunctions), clones) )

Q = Tensor('Q', qShape, alignStride=alignStride)
I = Tensor('I', qShape, alignStride=alignStride)
Ineigh = [Tensor('Ineigh[{}]'.format(i), qShape, alignStride=alignStride) for i in range(4)]

# Flux solver
AplusT = [Tensor('AplusT[{}]'.format(dim), (numberOfQuantities, numberOfQuantities)) for dim in range(4)]
AminusT = [Tensor('AminusT[{}]'.format(dim), (numberOfQuantities, numberOfQuantities)) for dim in range(4)]

# Kernels
g = Generator(arch)

volumeSum = Q[qi('kp')]
for i in range(3):
  volumeSum += db.kDivM[i][t('kl')] * I[qi('lq')] * db.star[i]['qp']
volume = (Q[qi('kp')] <= volumeSum)
g.add('volume', volume)

localFlux = lambda i: (Q[qi('kp')] <= db.rDivM[i][t('km')] * db.fMrT[i][t('ml')] * I[qi('lq')] * AplusT[i]['qp'])
g.addFamily('localFlux', simpleParameterSpace(4), localFlux)

neighbourFlux = lambda h,j,i: Q[qi('kp')] <= Q[qi('kp')] + db.rDivM[i][t('km')] * db.fP[h][t('mn')] * db.rT[j][t('nl')] * Ineigh[i][qi('lq')] * AminusT[i]['qp']
g.addFamily('neighboringFlux', simpleParameterSpace(3,4,4), neighbourFlux)

lastDQ = Q
for i in range(order-1):
  derivativeSum = Add()
  for j in range(3):
    derivativeSum += db.kDivMT[j][t('kl')] * lastDQ[qi('lq')] * db.star[j]['qp']
  derivativeSum = DeduceIndices( Q[qi('kp')].indices ).visit(derivativeSum)
  derivativeSum = EquivalentSparsityPattern().visit(derivativeSum)
  dQ = Tensor('dQ[{0}]'.format(i), qShape, spp=derivativeSum.eqspp(), alignStride=True)
  g.add('derivative[{}]'.format(i), dQ[qi('kp')] <= derivativeSum)
  lastDQ = dQ

#DynamicRupture.addKernels(g, Q, cmdLineArgs.matricesDir, order, cmdLineArgs.dynamicRuptureMethod, numberOfQuantities, numberOfQuantities)

# Generate code
g.generate(cmdLineArgs.outputDir, 'seissol')
