#!/usr/bin/env python3
##
# @file
# This file is part of SeisSol.
#
# @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
# @author Sebastian Wolf (wolf.sebastian AT tum.de, https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
#
# @section LICENSE
# Copyright (c) 2019, SeisSol Group
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
import acoustic
from yateto import Tensor, Scalar
from multSim import OptionalDimTensor

def addKernels(generator, aderdg):
  numberOf3DBasisFunctions = aderdg.numberOf3DBasisFunctions()
  numberOfQuantities = aderdg.numberOfQuantities()
  order = aderdg.order
  ## Point sources
  mStiffnessTensor = Tensor('stiffnessTensor', (3,3,3,3))
  mSlip = Tensor('mSlip', (3,))
  mNormal = Tensor('mNormal', (3,))
  mArea = Scalar('mArea')
  basisFunctionsAtPoint = Tensor('basisFunctionsAtPoint', (numberOf3DBasisFunctions,))
  basisFunctionDerivativesAtPoint = Tensor('basisFunctionDerivativesAtPoint', (numberOf3DBasisFunctions, 3))
  timeBasisFunctionsAtPoint = Tensor('timeBasisFunctionsAtPoint', (order,))
  mInvJInvPhisAtSources = Tensor('mInvJInvPhisAtSources', (numberOf3DBasisFunctions,))
  JInv = Scalar('JInv')

  generator.add('computeMInvJInvPhisAtSources',
    mInvJInvPhisAtSources['k'] <= JInv * aderdg.db.M3inv['kl'] * basisFunctionsAtPoint['l'])

  #extract the moment tensors entries in SeisSol ordering (xx, yy, zz, xy, yz, xz)
  if not isinstance(aderdg, acoustic.AcousticADERDG):
    assert(numberOfQuantities >= 6)
    momentToNRF_spp = np.zeros((numberOfQuantities, 3, 3))
    momentToNRF_spp[0, 0, 0] = 1
    momentToNRF_spp[1, 1, 1] = 1
    momentToNRF_spp[2, 2, 2] = 1
    momentToNRF_spp[3, 0, 1] = 1
    momentToNRF_spp[4, 1, 2] = 1
    momentToNRF_spp[5, 0, 2] = 1
  else:
    momentToNRF_spp = np.zeros((numberOfQuantities, 3, 3))
    momentToNRF_spp[0, 0, 0] = 1
  momentToNRF = Tensor('momentToNRF', (numberOfQuantities, 3, 3), spp=momentToNRF_spp) 

  momentNRFKernel = momentToNRF['tpq'] * mArea * mStiffnessTensor['pqij'] * mSlip['i'] * mNormal['j'] 

  if aderdg.Q.hasOptDim():
    sourceNRF = aderdg.Q['kt'] <= aderdg.Q['kt'] + mInvJInvPhisAtSources['k'] * momentNRFKernel * aderdg.oneSimToMultSim['s'] 
  else:
    sourceNRF = aderdg.Q['kt'] <= aderdg.Q['kt'] + mInvJInvPhisAtSources['k'] * momentNRFKernel 
  generator.add('sourceNRF', sourceNRF)

  momentFSRM = Tensor('momentFSRM', (numberOfQuantities,))
  stfIntegral = Scalar('stfIntegral')
  if aderdg.Q.hasOptDim():
    sourceFSRM = aderdg.Q['kp'] <= aderdg.Q['kp'] + stfIntegral * mInvJInvPhisAtSources['k'] * momentFSRM['p'] * aderdg.oneSimToMultSim['s']
  else:
    sourceFSRM = aderdg.Q['kp'] <= aderdg.Q['kp'] + stfIntegral * mInvJInvPhisAtSources['k'] * momentFSRM['p']
  generator.add('sourceFSRM', sourceFSRM)

  ## Receiver output
  QAtPoint = OptionalDimTensor('QAtPoint', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), (numberOfQuantities,))
  evaluateDOFSAtPoint = QAtPoint['p'] <= aderdg.Q['kp'] * basisFunctionsAtPoint['k']
  generator.add('evaluateDOFSAtPoint', evaluateDOFSAtPoint)
  QDerivativeAtPoint = OptionalDimTensor('QDerivativeAtPoint', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), (numberOfQuantities, 3))
  evaluateDerivativeDOFSAtPoint = QDerivativeAtPoint['pd'] <= aderdg.Q['kp'] * basisFunctionDerivativesAtPoint['kd']
  generator.add('evaluateDerivativeDOFSAtPoint', evaluateDerivativeDOFSAtPoint)

  stpShape = (numberOf3DBasisFunctions, numberOfQuantities, order)
  spaceTimePredictor = OptionalDimTensor('spaceTimePredictor', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), stpShape, alignStride=True)
  evaluateDOFSAtPointSTP = QAtPoint['p'] <= spaceTimePredictor['kpt'] * basisFunctionsAtPoint['k'] * timeBasisFunctionsAtPoint['t']
  generator.add('evaluateDOFSAtPointSTP', evaluateDOFSAtPointSTP)
  spaceTimePredictor = OptionalDimTensor('spaceTimePredictor', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), stpShape, alignStride=True)
  evaluateDerivativeDOFSAtPointSTP = QDerivativeAtPoint['pd'] <= spaceTimePredictor['kpt'] * basisFunctionDerivativesAtPoint['kd'] * timeBasisFunctionsAtPoint['t']
  generator.add('evaluateDerivativeDOFSAtPointSTP', evaluateDerivativeDOFSAtPointSTP)
