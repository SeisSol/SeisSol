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

def addKernels(generator, aderdg, matricesDir, PlasticityMethod):
  # Load matrices
  db = parseXMLMatrixFile('{}/plasticity_{}_matrices_{}.xml'.format(matricesDir, PlasticityMethod, aderdg.order), clones=dict(), alignStride=aderdg.alignStride)
  numberOfNodes = aderdg.t(db.v.shape())[0]

  sShape = (aderdg.numberOf3DBasisFunctions(), 6)
  QStress = OptionalDimTensor('QStress', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), sShape, alignStride=True)
  initialLoading = Tensor('initialLoading', (6,))

  replicateIniLShape = (numberOfNodes,)
  replicateIniLSpp = np.ones(aderdg.Q.insertOptDim(replicateIniLShape, (aderdg.Q.optSize(),)))
  replicateInitialLoading = OptionalDimTensor('replicateInitialLoading', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), replicateIniLShape, spp=replicateIniLSpp, alignStride=True)

  iShape = (numberOfNodes, 6)
  QStressNodal = OptionalDimTensor('QStressNodal', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), iShape, alignStride=True)
  meanStress = OptionalDimTensor('meanStress', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), (numberOfNodes,), alignStride=True)
  secondInvariant = OptionalDimTensor('secondInvariant', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), (numberOfNodes,), alignStride=True)

  selectBulkAverage = Tensor('selectBulkAverage', (6,), spp={(i,): str(1.0/3.0) for i in range(3)})
  selectBulkNegative = Tensor('selectBulkNegative', (6,), spp={(i,): '-1.0' for i in range(3)})
  weightSecondInvariant = Tensor('weightSecondInvariant', (6,), spp={(i,): str(1.0/2.0) if i < 3 else '1.0' for i in range(6)})
  yieldFactor = Tensor('yieldFactor', (numberOfNodes,))

  generator.add('plConvertToNodal', QStressNodal['kp'] <= db.v[aderdg.t('kl')] * QStress['lp'] + replicateInitialLoading['k'] * initialLoading['p'])
  generator.add('plComputeMean', meanStress['k'] <= QStressNodal['kq'] * selectBulkAverage['q'])
  generator.add('plSubtractMean', QStressNodal['kp'] <= QStressNodal['kp'] + meanStress['k'] * selectBulkNegative['p'])
  generator.add('plComputeSecondInvariant', secondInvariant['k'] <= QStressNodal['kq'] * QStressNodal['kq'] * weightSecondInvariant['q'])
  generator.add('plAdjustStresses', QStress['kp'] <= QStress['kp'] + db.vInv[aderdg.t('kl')] * QStressNodal['lp'] * yieldFactor['l'])
