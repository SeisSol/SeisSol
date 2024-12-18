#!/usr/bin/env python
##
# @file
# This file is part of SeisSol.
#
# @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
#
# @section LICENSE
# Copyright (c) 2017, SeisSol Group
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
from yateto import Tensor
from yateto.input import parseXMLMatrixFile
from multSim import OptionalDimTensor
from common import generate_kernel_name_prefix

def addKernels(generator, aderdg, matricesDir, PlasticityMethod, targets):
  # Load matrices
  db = parseXMLMatrixFile('{}/plasticity_{}_matrices_{}.xml'.format(matricesDir, PlasticityMethod, aderdg.order), clones=dict(), alignStride=aderdg.alignStride)
  numberOfNodes = aderdg.t(db.v.shape())[0]

  numberOf3DBasisFunctions = aderdg.numberOf3DBasisFunctions()
  sShape = (numberOf3DBasisFunctions, 6)

  sShape_eta = (numberOf3DBasisFunctions,)
  iShape = (numberOfNodes, 6)
  replicateIniLShape = (numberOfNodes,)
  replicateIniLSpp = np.ones(aderdg.Q.insertOptDim(replicateIniLShape, (aderdg.Q.optSize(),)))
  
  pstrainShape = (numberOfNodes, 7)

  pstrain = OptionalDimTensor('pstrain', 's', aderdg.multipleSimulations, 0, pstrainShape, alignStride=True)
  pstrain_ijs = Tensor('pstrain_ijs', (pstrainShape[0], pstrainShape[1], aderdg.multipleSimulations))
  pstrainModified = pstrain_ijs['ijs'] <= pstrain['ij']
  generator.add('pstrainModified', pstrainModified)
  pstrainModifiedReversed = pstrain['ij'] <= pstrain_ijs['ijs']
  generator.add('pstrainModifiedReversed', pstrainModifiedReversed)

  if aderdg.multipleSimulations > 1:
    QEtaModal = Tensor('QEtaModal', sShape_eta, alignStride=False)
    QStressNodal = Tensor('QStressNodal', iShape, alignStride=False)
    meanStress = Tensor('meanStress', (numberOfNodes, ), alignStride=False)
    QEtaNodal = Tensor('QEtaNodal', (numberOfNodes, ), alignStride=False)
    secondInvariant = Tensor('secondInvariant', (numberOfNodes, ), alignStride=False)
    QStress = Tensor('QStress', sShape, alignStride=False)
    replicateInitialLoading = Tensor('replicateInitialLoading', replicateIniLShape, spp=np.ones(replicateIniLShape, ), alignStride=False)
  else:
    QEtaModal = Tensor('QEtaModal', sShape_eta, alignStride=True)
    QStressNodal = Tensor('QStressNodal', iShape, alignStride=True)
    meanStress = Tensor('meanStress', (numberOfNodes, ), alignStride=True)
    QEtaNodal = Tensor('QEtaNodal', (numberOfNodes, ), alignStride=True)
    secondInvariant = Tensor('secondInvariant', (numberOfNodes, ), alignStride=True)
    QStress = Tensor('QStress', sShape, alignStride=True)
    replicateInitialLoading = Tensor('replicateInitialLoading', replicateIniLShape, spp=np.ones(replicateIniLShape, ), alignStride=True)
  
  initialLoading = Tensor('initialLoading', (6,))

  selectBulkAverage = Tensor('selectBulkAverage', (6,), spp={(i,): str(1.0/3.0) for i in range(3)})
  selectBulkNegative = Tensor('selectBulkNegative', (6,), spp={(i,): '-1.0' for i in range(3)})
  weightSecondInvariant = Tensor('weightSecondInvariant', (6,), spp={(i,): str(1.0/2.0) if i < 3 else '1.0' for i in range(6)})
  yieldFactor = Tensor('yieldFactor', (numberOfNodes,))

  generator.add('plConvertToNodal', QStressNodal['kp'] <= db.v[aderdg.t('kl')] * QStress['lp'] + replicateInitialLoading['k'] * initialLoading['p'])

  for target in targets:
    name_prefix = generate_kernel_name_prefix(target)
    generator.add(name=f'{name_prefix}plConvertToNodalNoLoading',
                  ast=QStressNodal['kp'] <= db.v[aderdg.t('kl')] * QStress['lp'],
                  target=target)

    generator.add(name=f'{name_prefix}plConvertEtaModal2Nodal',
                  ast=QEtaNodal['k'] <= db.v[aderdg.t('kl')] * QEtaModal['l'],
                  target=target)

    generator.add(name=f'{name_prefix}plConvertEtaNodal2Modal',
                  ast=QEtaModal['k'] <= db.vInv[aderdg.t('kl')] * QEtaNodal['l'],
                  target=target)

  generator.add('plComputeMean', meanStress['k'] <= QStressNodal['kq'] * selectBulkAverage['q'])
  generator.add('plSubtractMean', QStressNodal['kp'] <= QStressNodal['kp'] + meanStress['k'] * selectBulkNegative['p'])
  generator.add('plComputeSecondInvariant', secondInvariant['k'] <= QStressNodal['kq'] * QStressNodal['kq'] * weightSecondInvariant['q'])
  generator.add('plAdjustStresses', QStress['kp'] <= QStress['kp'] + db.vInv[aderdg.t('kl')] * QStressNodal['lp'] * yieldFactor['l'])

  gpu_target = 'gpu'
  if gpu_target in targets:
    name_prefix = generate_kernel_name_prefix(gpu_target)

    # suffix `M` stands for `Matrix`
    replicateInitialLoadingM = Tensor(name='replicateInitialLoadingM',
                                      shape=(numberOfNodes, 1),
                                      spp=np.ones((numberOfNodes, 1)))
    initialLoadingM = Tensor('initialLoadingM', (1, 6))
    # initialLoadingM = OptionalDimTensor('initialLoadingM', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), (1,6), alignStride=True)

    # Note: the last term was change on purpose because
    # GemmForge doesn't currently support tensor product operation
    convert_to_nodal = QStressNodal['kp'] <= \
                       db.v[aderdg.t('kl')] * QStress['lp'] + \
                       replicateInitialLoadingM['km'] * initialLoadingM['mp']

    generator.add(name=f'{name_prefix}plConvertToNodal',
                  ast=convert_to_nodal,
                  target=gpu_target)


    generator.add(f'{name_prefix}plConvertToModal',
                  QStress['kp'] <= QStress['kp'] + db.vInv[aderdg.t('kl')] * QStressNodal['lp'],
                  target=gpu_target)

