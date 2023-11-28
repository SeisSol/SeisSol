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
from common import *
from yateto import Tensor, Scalar, simpleParameterSpace
from yateto.input import parseJSONMatrixFile
from multSim import OptionalDimTensor
from copy import deepcopy
import numpy as np

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
  QInterpolated = OptionalDimTensor('QInterpolated', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), gShape, alignStride=True)

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
  resampleKernel = resampledQ['i'] <= db.resample['ij'] * originalQ['j']
  generator.add('resampleParameter', resampleKernel )

  # DR filter
  numberOf2DBasisFunctions = aderdg.numberOf2DBasisFunctions()
  # Filter matrix
  drFilter = Tensor("drFilter", (numberOfPoints, numberOfPoints))
  # Filter weights matrix
  filterWeights = Tensor("filterWeights", (numberOf2DBasisFunctions, numberOf2DBasisFunctions))
  # Compute filter matrix kernel
  generator.add('computeFilterMatrix', drFilter['qp'] <= db.V2mTo2Quad['ql'] * filterWeights['lk'] * db.V2QuadTo2m['kp'])
  filteredQ = Tensor('filteredQ', (numberOfPoints,))
  generator.add('filterParameter', filteredQ['i'] <= drFilter['ij'] * originalQ['j'])

  generator.add('transposeTinv', TinvT['ij'] <= aderdg.Tinv['ji'])

  fluxScale = Scalar('fluxScale')
  generator.add('rotateFluxMatrix', fluxSolver['qp'] <= fluxScale * aderdg.starMatrix(0)['qk'] * aderdg.T['pk'])

  numberOf3DBasisFunctions = aderdg.numberOf3DBasisFunctions()
  numberOfQuantities = aderdg.numberOfQuantities()
  basisFunctionsAtPoint = Tensor('basisFunctionsAtPoint', (numberOf3DBasisFunctions,))
  QAtPoint = OptionalDimTensor('QAtPoint', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), (numberOfQuantities,))

  generator.add('evaluateFaceAlignedDOFSAtPoint',
                QAtPoint['q'] <= aderdg.Tinv['qp'] * aderdg.Q['lp'] * basisFunctionsAtPoint['l'])

  def interpolateQGenerator(i,h):
    return QInterpolated['kp'] <= db.V3mTo2n[i,h][aderdg.t('kl')] * aderdg.Q['lq'] * TinvT['qp']

  interpolateQPrefetch = lambda i,h: QInterpolated
  for target in targets:
    name_prefix = generate_kernel_name_prefix(target)
    generator.addFamily(f'{name_prefix}evaluateAndRotateQAtInterpolationPoints',
                        simpleParameterSpace(4,4),
                        interpolateQGenerator,
                        interpolateQPrefetch if target == 'cpu' else None,
                        target=target)

  nodalFluxGenerator = lambda i,h: aderdg.extendedQTensor()['kp'] <= aderdg.extendedQTensor()['kp'] + db.V3mTo2nTWDivM[i,h][aderdg.t('kl')] * QInterpolated['lq'] * fluxSolver['qp']
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
  QInterpolatedPlus = OptionalDimTensor('QInterpolatedPlus', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), gShape, alignStride=True)
  QInterpolatedMinus = OptionalDimTensor('QInterpolatedMinus', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), gShape, alignStride=True)
  slipRateInterpolated = Tensor('slipRateInterpolated', (numberOfPoints,3), alignStride=True)
  tractionInterpolated = Tensor('tractionInterpolated', (numberOfPoints,3), alignStride=True)
  frictionalEnergy = Tensor('frictionalEnergy', ())
  timeWeight = Scalar('timeWeight')
  spaceWeights = Tensor('spaceWeights', (numberOfPoints,), alignStride=True)

  computeTractionInterpolated = tractionInterpolated['kp'] <= QInterpolatedMinus['kq'] * aderdg.tractionMinusMatrix['qp'] + QInterpolatedPlus['kq'] * aderdg.tractionPlusMatrix['qp']
  generator.add('computeTractionInterpolated', computeTractionInterpolated)

  accumulateFrictionalEnergy = frictionalEnergy[''] <= frictionalEnergy[''] + timeWeight * tractionInterpolated['kp'] * slipRateInterpolated['kp'] * spaceWeights['k']
  generator.add('accumulateFrictionalEnergy', accumulateFrictionalEnergy)

  ## Dynamic Rupture Precompute
  qPlus = OptionalDimTensor('Qplus', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), gShape, alignStride=True)
  qMinus = OptionalDimTensor('Qminus', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), gShape, alignStride=True)

  extractVelocitiesSPP = aderdg.extractVelocities()
  extractVelocities = Tensor('extractVelocities', extractVelocitiesSPP.shape, spp=extractVelocitiesSPP)
  extractTractionsSPP = aderdg.extractTractions()
  extractTractions = Tensor('extractTractions', extractTractionsSPP.shape, spp=extractTractionsSPP)

  N = extractTractionsSPP.shape[0]
  eta = Tensor('eta', (N,N))
  zPlus = Tensor('Zplus', (N,N))
  zMinus = Tensor('Zminus', (N,N))
  theta = OptionalDimTensor('theta', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos, (numberOfPoints, N), alignStride=True)

  velocityJump = extractVelocities['lj'] * qMinus['ij'] - extractVelocities['lj'] * qPlus['ij']
  tractionsPlus = extractTractions['mn'] * qPlus['in']
  tractionsMinus = extractTractions['mn'] * qMinus['in']
  computeTheta = theta['ik'] <= eta['kl'] * velocityJump + eta['kl'] * zPlus['lm'] * tractionsPlus + eta['kl'] * zMinus['lm'] * tractionsMinus
  generator.add('computeTheta', computeTheta)

  mapToVelocitiesSPP = aderdg.mapToVelocities()
  mapToVelocities = Tensor('mapToVelocities', mapToVelocitiesSPP.shape, spp=mapToVelocitiesSPP)
  mapToTractionsSPP = aderdg.mapToTractions()
  mapToTractions = Tensor('mapToTractions', mapToTractionsSPP.shape, spp=mapToTractionsSPP)
  imposedState = OptionalDimTensor('imposedState', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), gShape, alignStride=True)
  weight = Scalar('weight')
  computeImposedStateM = imposedState['ik'] <= imposedState['ik'] + weight * mapToVelocities['kl'] * (extractVelocities['lm'] * qMinus['im'] - zMinus['lm'] * theta['im'] + zMinus['lm'] * tractionsMinus) + weight * mapToTractions['kl'] * theta['il']
  computeImposedStateP = imposedState['ik'] <= imposedState['ik'] + weight * mapToVelocities['kl'] * (extractVelocities['lm'] * qPlus['im'] - zPlus['lm'] * tractionsPlus + zPlus['lm'] * theta['im']) + weight * mapToTractions['kl'] * theta['il']
  generator.add('computeImposedStateM', computeImposedStateM)
  generator.add('computeImposedStateP', computeImposedStateP)


  return {db.resample, db.quadpoints, db.quadweights}
