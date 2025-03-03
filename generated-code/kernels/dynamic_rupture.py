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
from kernels.common import *
from yateto import Tensor, Scalar, simpleParameterSpace
from yateto.ast.node import Add
from yateto.ast.transformer import DeduceIndices, EquivalentSparsityPattern
from yateto.input import parseJSONMatrixFile
from kernels.multsim import OptionalDimTensor
from copy import deepcopy

def addKernels(generator, aderdg, matricesDir, drQuadRule, targets):

  clones = dict()

  # Load matrices
  db = parseJSONMatrixFile(f'{matricesDir}/dr_{drQuadRule}_matrices_{aderdg.order}.json', clones, alignStride=aderdg.alignStride, transpose=aderdg.transpose)  
  numberOfPoints = db.resample.shape()[0]

  # Determine matrices
  # Note: This does only work because the flux does not depend on the mechanisms in the case of viscoelastic attenuation
  trans_inv_spp_T = aderdg.transformation_inv_spp().transpose()
  TinvT = Tensor('TinvT', trans_inv_spp_T.shape, spp=trans_inv_spp_T)
  flux_solver_spp = aderdg.flux_solver_spp()
  fluxSolver    = Tensor('fluxSolver', flux_solver_spp.shape, spp=flux_solver_spp)
  
  gShape = (numberOfPoints, aderdg.numberOfQuantities())
  # QInterpolated = OptionalDimTensor('QInterpolated', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), gShape, alignStride=True)
  # TODO: (VK) Make this work with the original tensors
  QInterpolated = Tensor('QInterpolated', gShape, alignStride=False)

  stressRotationMatrix = Tensor("stressRotationMatrix", (6, 6))
  initialStress = Tensor("initialStress", (6, ))
  rotatedStress = Tensor("rotatedStress", (6, ))
  rotationKernel = rotatedStress['i'] <= stressRotationMatrix['ij'] * initialStress['j']
  generator.add('rotateStress', rotationKernel)

  reducedFaceAlignedMatrix = Tensor("reducedFaceAlignedMatrix", (6, 6))
  generator.add('rotateInitStress',
                rotatedStress['k'] <= stressRotationMatrix['ki'] * reducedFaceAlignedMatrix['ij'] * initialStress['j'])

  originalQ = Tensor('originalQ', (numberOfPoints,))
  resampledQ = Tensor('resampledQ', (numberOfPoints,))
  if aderdg.multipleSimulations > 1:
    resampleKernel = resampledQ['i'] <= db.resample['ji'] * originalQ['j']
  else:
    resampleKernel = resampledQ['i'] <= db.resample['ij'] * originalQ['j']
  generator.add('resampleParameter', resampleKernel )

  generator.add('transposeTinv', TinvT['ij'] <= aderdg.Tinv['ji'])

  fluxScale = Scalar('fluxScaleDR')
  generator.add('rotateFluxMatrix', fluxSolver['qp'] <= fluxScale * aderdg.starMatrix(0)['qk'] * aderdg.T['pk'])

  numberOf3DBasisFunctions = aderdg.numberOf3DBasisFunctions()
  numberOfQuantities = aderdg.numberOfQuantities()
  basisFunctionsAtPoint = Tensor('basisFunctionsAtPoint', (numberOf3DBasisFunctions,))
  qShape = (aderdg.numberOf3DBasisFunctions(), aderdg.numberOfQuantities())
  QuantitiesSingleSim = Tensor('singleSimQ', qShape, alignStride=False)
  QAtPoint = Tensor('QAtPoint', (aderdg.numberOfQuantities(), ))
  # QAtPoint = OptionalDimTensor('QAtPoint', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), (numberOfQuantities,))

  # generator.add('evaluateFaceAlignedDOFSAtPoint',
  #               QAtPoint['q'] <= aderdg.Tinv['qp'] * aderdg.Q['lp'] * basisFunctionsAtPoint['l'])

  generator.add('evaluateFaceAlignedDOFSAtPoint', QAtPoint['q'] <= aderdg.Tinv['qp']*QuantitiesSingleSim['lp']*basisFunctionsAtPoint['l'])
  
  dQ0_DR = Tensor('dQ_DR(0)', qShape, alignStride=False)
  I_DR = Tensor('I_DR', qShape, alignStride=False)
  powers_DR = [Scalar(f'power_DR({i})') for i in range(aderdg.order)]
  power_DR = powers_DR[0]
  derivatives_DR = [dQ0_DR]
  derivativeExpr_DR = [I_DR['kp'] <= power_DR * dQ0_DR['kp']]
  derivativeTaylorExpansion_DR = power_DR*dQ0_DR['kp']

  for i in range(1, aderdg.order):
    power_DR = powers_DR[i]
    derivativeSum_DR = Add()
    for j in range(3):
      derivativeSum_DR += aderdg.db.kDivMT[j][aderdg.t('kl')] * derivatives_DR[-1]['lq'] * aderdg.starMatrix(j)['qp']

    derivativeSum_DR = DeduceIndices(QuantitiesSingleSim['kp'].indices).visit(derivativeSum_DR)
    derivativeSum_DR = EquivalentSparsityPattern().visit(derivativeSum_DR)
    dQ_DR = Tensor('dQ_DR({})'.format(i), qShape, spp=derivativeSum_DR.eqspp(), alignStride=False)
    derivativeTaylorExpansion_DR += power_DR * dQ_DR['kp']

    derivatives_DR.append(dQ_DR)
  

  
  derivativeTaylorExpansion_DR = I_DR['kp'] <= derivativeTaylorExpansion_DR
  for target in targets:
    name_prefix = generate_kernel_name_prefix(target)
    generator.add(f'{name_prefix}derivativeTaylorExpansion_DR', derivativeTaylorExpansion_DR, target=target)
  
  # TODO : (VK) Make this work with the original tensors 
  def interpolateQGenerator(i,h):
    return QInterpolated['kp'] <= db.V3mTo2n[i,h][aderdg.t('kl')] * QuantitiesSingleSim['lq'] * TinvT['qp']

  interpolateQPrefetch = lambda i,h: QInterpolated
  for target in targets:
    name_prefix = generate_kernel_name_prefix(target)
    generator.addFamily(f'{name_prefix}evaluateAndRotateQAtInterpolationPoints',
                        simpleParameterSpace(4,4),
                        interpolateQGenerator,
                        interpolateQPrefetch if target == 'cpu' else None,
                        target=target)

  QInterpolatedSingleSim = Tensor('QInterpolatedSingleSim', gShape, alignStride=False)
#  nodalFluxGenerator = lambda i,h: aderdg.extendedQTensor()['kp'] <= aderdg.extendedQTensor()['kp'] + db.V3mTo2nTWDivM[i,h][aderdg.t('kl')] * QInterpolated['lq'] * fluxSolver['qp']
  # TODO (VK): make this work in the original tensor format by making the other variables in tensor format later
  nodalFluxGenerator = lambda i,h: QuantitiesSingleSim['kp'] <= QuantitiesSingleSim['kp'] + db.V3mTo2nTWDivM[i,h][aderdg.t('kl')] * QInterpolatedSingleSim['lq'] * fluxSolver['qp']
  nodalFluxPrefetch = lambda i,h: aderdg.I

  for target in targets:
    name_prefix = generate_kernel_name_prefix(target)
    generator.addFamily(f'{name_prefix}nodalFlux',
                        simpleParameterSpace(4,4),
                        nodalFluxGenerator,
                        nodalFluxPrefetch if target =='cpu' else None,
                        target=target)

  # Energy output
  # Minus and plus refer to the original implementation of Christian Pelties,
  # where the normal points from the plus side to the minus side
  # QInterpolatedPlus = OptionalDimTensor('QInterpolatedPlus', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), gShape, alignStride=True)
  QInterpolatedPlus = Tensor('QInterpolatedPlus', gShape, alignStride=False)
  # QInterpolatedMinus = OptionalDimTensor('QInterpolatedMinus', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), gShape, alignStride=True)
  QInterpolatedMinus = Tensor('QInterpolatedMinus', gShape, alignStride=False)
  # slipInterpolated = OptionalDimTensor('slipInterpolated', 's', aderdg.multipleSimulations, 0, (numberOfPoints,3), alignStride=True)
  slipInterpolated = Tensor('slipInterpolated', (numberOfPoints, 3), alignStride=False)
  # tractionInterpolated = OptionalDimTensor('tractionInterpolated', 's', aderdg.multipleSimulations, 0, (numberOfPoints,3), alignStride=True)
  tractionInterpolated = Tensor('tractionInterpolated', (numberOfPoints, 3), alignStride=False)
  # staticFrictionalWork = OptionalDimTensor('staticFrictionalWork', 's', aderdg.multipleSimulations, 0, (1,), alignStride=True)
  staticFrictionalWork = Tensor('staticFrictionalWork', (1,), alignStride=False)
  minusSurfaceArea = Scalar('minusSurfaceArea')
  spaceWeights = Tensor('spaceWeights', (numberOfPoints, 1), alignStride=False)

  computeTractionInterpolated = tractionInterpolated['kp'] <= QInterpolatedMinus['kq'] * aderdg.tractionMinusMatrix['qp'] + QInterpolatedPlus['kq'] * aderdg.tractionPlusMatrix['qp']
  generator.add('computeTractionInterpolated', computeTractionInterpolated)

  accumulateStaticFrictionalWork = staticFrictionalWork['l'] <= staticFrictionalWork['l'] + minusSurfaceArea * tractionInterpolated['kp'] * slipInterpolated['kp'] * spaceWeights['kl']
  generator.add('accumulateStaticFrictionalWork', accumulateStaticFrictionalWork)

  ## Dynamic Rupture Precompute
  # qPlus = OptionalDimTensor('Qplus', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), gShape, alignStride=True)
  qPlus = Tensor('Qplus', gShape, alignStride=False)
  # qMinus = OptionalDimTensor('Qminus', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), gShape, alignStride=True)
  qMinus = Tensor('Qminus', gShape, alignStride=False)

  extractVelocitiesSPP = aderdg.extractVelocities()
  extractVelocities = Tensor('extractVelocities', extractVelocitiesSPP.shape, spp=extractVelocitiesSPP)
  extractTractionsSPP = aderdg.extractTractions()
  extractTractions = Tensor('extractTractions', extractTractionsSPP.shape, spp=extractTractionsSPP)

  N = extractTractionsSPP.shape[0]
  eta = Tensor('eta', (N,N))
  zPlus = Tensor('Zplus', (N,N))
  zMinus = Tensor('Zminus', (N,N))
  # theta = OptionalDimTensor('theta', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), (numberOfPoints, N), alignStride=True)
  theta = Tensor('theta', (numberOfPoints, N), alignStride=False)
  velocityJump = extractVelocities['lj'] * qMinus['ij'] - extractVelocities['lj'] * qPlus['ij']
  tractionsPlus = extractTractions['mn'] * qPlus['in']
  tractionsMinus = extractTractions['mn'] * qMinus['in']
  computeTheta = theta['ik'] <= eta['kl'] * velocityJump + eta['kl'] * zPlus['lm'] * tractionsPlus + eta['kl'] * zMinus['lm'] * tractionsMinus
  generator.add('computeTheta', computeTheta)

  mapToVelocitiesSPP = aderdg.mapToVelocities()
  mapToVelocities = Tensor('mapToVelocities', mapToVelocitiesSPP.shape, spp=mapToVelocitiesSPP)
  mapToTractionsSPP = aderdg.mapToTractions()
  mapToTractions = Tensor('mapToTractions', mapToTractionsSPP.shape, spp=mapToTractionsSPP)
  # imposedState = OptionalDimTensor('imposedState', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), gShape, alignStride=True)
  imposedState = Tensor('imposedState', gShape, alignStride=False)
  weight = Scalar('weight')
  computeImposedStateM = imposedState['ik'] <= imposedState['ik'] + weight * mapToVelocities['kl'] * (extractVelocities['lm'] * qMinus['im'] - zMinus['lm'] * theta['im'] + zMinus['lm'] * tractionsMinus) + weight * mapToTractions['kl'] * theta['il']
  computeImposedStateP = imposedState['ik'] <= imposedState['ik'] + weight * mapToVelocities['kl'] * (extractVelocities['lm'] * qPlus['im'] - zPlus['lm'] * tractionsPlus + zPlus['lm'] * theta['im']) + weight * mapToTractions['kl'] * theta['il']
  generator.add('computeImposedStateM', computeImposedStateM)
  generator.add('computeImposedStateP', computeImposedStateP)

  return {db.resample, db.quadpoints, db.quadweights}
