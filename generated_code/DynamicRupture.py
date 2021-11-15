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



from common import *
from yateto import Tensor, Scalar, simpleParameterSpace
from yateto.input import parseJSONMatrixFile
from multSim import OptionalDimTensor
#adrian numpy added:
import numpy as np

def addKernels(generator, aderdg, matricesDir, dynamicRuptureMethod, targets):
  if dynamicRuptureMethod == 'quadrature':
    numberOfPoints = (aderdg.order+1)**2
  else:
    raise ValueError('Unknown dynamic rupture method.')

  clones = dict()

  # Load matrices
  db = parseJSONMatrixFile('{}/dr_{}_matrices_{}.json'.format(matricesDir, dynamicRuptureMethod, aderdg.order), clones, alignStride=aderdg.alignStride, transpose=aderdg.transpose)
  db.update( parseJSONMatrixFile('{}/resample_{}.json'.format(matricesDir, aderdg.order)) )

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
  generator.add('rotateStressToFaultCS', rotationKernel )


  resamplePar = Tensor('resamplePar', (numberOfPoints,))
  resampledPar = Tensor('resampledPar', (numberOfPoints,))
  resampleM = Tensor('resampleM', (numberOfPoints, numberOfPoints) )
  resampleKernel = resampledPar['i'] <= resampleM['ij'] * resamplePar['j']
  generator.add('resampleParameter', resampleKernel )


  QInterpolatedPlus = OptionalDimTensor('QInterpolatedPlus', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), gShape, alignStride=True)
  QInterpolatedMinus = OptionalDimTensor('QInterpolatedMinus', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), gShape, alignStride=True)
  NorStressGP = Tensor('NorStressGP', (numberOfPoints,))
  XYStressGP = Tensor('XYStressGP', (numberOfPoints,))
  XZStressGP = Tensor('XZStressGP', (numberOfPoints,))
  eta_p = Scalar('eta_p')
  eta_s = Scalar('eta_s')
  inv_Zp = Scalar('inv_Zp')
  inv_Zp_neig = Scalar('inv_Zp_neig')
  inv_Zs = Scalar('inv_Zs')
  inv_Zs_neig = Scalar('inv_Zs_neig')

  select0Spp = np.zeros(aderdg.numberOfQuantities())
  select0Spp[0] = 1
  select0 = Tensor('select0', select0Spp.shape, select0Spp,)

  select6Spp = np.zeros(aderdg.numberOfQuantities())
  select6Spp[6] = 1
  select6 = Tensor('select6', select6Spp.shape, select6Spp,)

  select7Spp = np.zeros(aderdg.numberOfQuantities())
  select7Spp[7] = 1
  select7 = Tensor('select7', select7Spp.shape, select7Spp,)

  select8Spp = np.zeros(aderdg.numberOfQuantities())
  select8Spp[8] = 1
  select8 = Tensor('select8', select8Spp.shape, select8Spp,)

  select3Spp = np.zeros(aderdg.numberOfQuantities())
  select3Spp[3] = 1
  select3 = Tensor('select3', select3Spp.shape, select3Spp,)

  select5Spp = np.zeros(aderdg.numberOfQuantities())
  select5Spp[5] = 1
  select5 = Tensor('select5', select5Spp.shape, select5Spp,)

  NorStressFromQInterpolatedKernel = NorStressGP['i'] <= eta_p * \
                                     ( \
                                               select6['k'] * QInterpolatedMinus['ik'] - select6['k'] * QInterpolatedPlus['ik'] + \
                                               inv_Zp * select0['k'] * QInterpolatedPlus['ik'] + inv_Zp_neig * select0['k'] * QInterpolatedMinus['ik'] \
                                       )

  XYStressFromQInterpolatedKernel = XYStressGP['i'] <= eta_s * \
                                    ( \
                                              select7['k'] * QInterpolatedMinus['ik'] - select7['k'] * QInterpolatedPlus['ik'] + \
                                              inv_Zs * select3['k'] * QInterpolatedPlus['ik'] +  inv_Zs_neig * select3['k']  * QInterpolatedMinus['ik'] \
                                      )

  XZStressFromQInterpolatedKernel = XZStressGP['i'] <= eta_s * \
                                    ( \
                                              select8['k'] * QInterpolatedMinus['ik'] - select8['k'] * QInterpolatedPlus['ik'] + \
                                              inv_Zs * select5['k'] * QInterpolatedPlus['ik'] +  inv_Zs_neig * select5['k']  * QInterpolatedMinus['ik'] \
                                      )

  generator.add('StressFromQInterpolated', [NorStressFromQInterpolatedKernel, XYStressFromQInterpolatedKernel, XZStressFromQInterpolatedKernel] )

  timeWeights = Scalar('timeWeights')
  TractionGP_XY = Tensor('TractionGP_XY', (numberOfPoints,))
  TractionGP_XZ = Tensor('TractionGP_XZ', (numberOfPoints,))

  imposedStatePlus = OptionalDimTensor('imposedStatePlus', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), gShape, alignStride=True)
  imposedStateMinus = OptionalDimTensor('imposedStateMinus', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), gShape, alignStride=True)

  imposedStatePlus0 = imposedStatePlus['ik'] <= imposedStatePlus['ik'] + select0['k'] * timeWeights * NorStressGP['i']
  imposedStatePlus3 = imposedStatePlus['ik'] <= imposedStatePlus['ik'] + select3['k'] * timeWeights * TractionGP_XY['i']
  imposedStatePlus5 = imposedStatePlus['ik'] <= imposedStatePlus['ik'] + select5['k'] * timeWeights * TractionGP_XZ['i']
  imposedStatePlus6 = imposedStatePlus['ik'] <= imposedStatePlus['ik'] + select6['k'] * timeWeights * ( \
            select6['l'] * QInterpolatedPlus['il'] + ( NorStressGP['i'] - select0['l']*QInterpolatedPlus['il'] ) * inv_Zp \
    )
  imposedStatePlus7 = imposedStatePlus['ik'] <= imposedStatePlus['ik'] + select7['k'] * timeWeights * ( \
            select7['l'] * QInterpolatedPlus['il'] +  ( TractionGP_XY['i'] - select3['l'] * QInterpolatedPlus['il'] ) * inv_Zs \
    )
  imposedStatePlus8 = imposedStatePlus['ik'] <= imposedStatePlus['ik'] + select8['k'] * timeWeights * ( \
            select8['l'] * QInterpolatedPlus['il'] +  ( TractionGP_XZ['i'] - select5['l']*QInterpolatedPlus['il'] ) * inv_Zs \
    )

  imposedStateMinus0 = imposedStateMinus['ik'] <= imposedStateMinus['ik'] + select0['k'] * timeWeights * NorStressGP['i']
  imposedStateMinus3 = imposedStateMinus['ik'] <= imposedStateMinus['ik'] + select3['k'] * timeWeights * TractionGP_XY['i']
  imposedStateMinus5 = imposedStateMinus['ik'] <= imposedStateMinus['ik'] + select5['k'] * timeWeights * TractionGP_XZ['i']
  imposedStateMinus6 = imposedStateMinus['ik'] <= imposedStateMinus['ik'] + select6['k'] * timeWeights * ( \
            select6['l'] * QInterpolatedMinus['il'] + ( select0['l']*QInterpolatedMinus['il'] - NorStressGP['i'] ) * inv_Zp_neig \
    )
  imposedStateMinus7 = imposedStateMinus['ik'] <= imposedStateMinus['ik'] + select7['k'] * timeWeights * ( \
            select7['l'] * QInterpolatedMinus['il'] + ( select3['l']*QInterpolatedMinus['il'] - TractionGP_XY['i'] ) * inv_Zs_neig \
    )
  imposedStateMinus8 = imposedStateMinus['ik'] <= imposedStateMinus['ik'] + select8['k'] * timeWeights * ( \
            select8['l'] * QInterpolatedMinus['il'] + ( select5['l']*QInterpolatedMinus['il'] - TractionGP_XZ['i'] ) * inv_Zs_neig \
    )

  generator.add('ImposedStateFromNewStress', \
                [imposedStatePlus0, imposedStatePlus3, imposedStatePlus5, imposedStatePlus6, imposedStatePlus7, imposedStatePlus8, \
                 imposedStateMinus0, imposedStateMinus3, imposedStateMinus5, imposedStateMinus6, imposedStateMinus7, imposedStateMinus8] )

  generator.add('transposeTinv', TinvT['ij'] <= aderdg.Tinv['ji'])

  fluxScale = Scalar('fluxScale')
  generator.add('rotateFluxMatrix', fluxSolver['qp'] <= fluxScale * aderdg.starMatrix(0)['qk'] * aderdg.T['pk'])

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

  return {db.resample}
