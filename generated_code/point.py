#!/usr/bin/env python3
##
# @file
# This file is part of SeisSol.
#
# @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
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
from yateto import *
from multSim import OptionalDimTensor

def addKernels(g, Q, oneSimToMultSim, numberOf3DBasisFunctions, numberOfQuantities):
  ## Point sources
  mInvJInvPhisAtSources = Tensor('mInvJInvPhisAtSources', (numberOf3DBasisFunctions,))

  momentNRF = Tensor('momentNRF', (numberOfQuantities,), spp=np.array([1]*6 + [0]*(numberOfQuantities-6), dtype=bool))
  if Q.hasOptDim():
    sourceNRF = Q['kp'] <= Q['kp'] - mInvJInvPhisAtSources['k'] * momentNRF['p'] * oneSimToMultSim['s']
  else:
    sourceNRF = Q['kp'] <= Q['kp'] - mInvJInvPhisAtSources['k'] * momentNRF['p']
  g.add('sourceNRF', sourceNRF)

  momentFSRM = Tensor('momentFSRM', (numberOfQuantities,))
  stfIntegral = Scalar('stfIntegral')
  if Q.hasOptDim():
    sourceFSRM = Q['kp'] <= Q['kp'] + stfIntegral * mInvJInvPhisAtSources['k'] * momentFSRM['p'] * oneSimToMultSim['s']
  else:
    sourceFSRM = Q['kp'] <= Q['kp'] + stfIntegral * mInvJInvPhisAtSources['k'] * momentFSRM['p']
  g.add('sourceFSRM', sourceFSRM)

  ## Receiver output
  basisFunctionsAtPoint = Tensor('basisFunctions', (numberOf3DBasisFunctions,))
  QAtPoint = OptionalDimTensor('QAtPoint', 's', Q.optSize(), 0, (numberOfQuantities,))
  evaluateDOFSAtPoint = QAtPoint['p'] <= Q['kp'] * basisFunctionsAtPoint['k']
  g.add('evaluateDOFSAtPoint', evaluateDOFSAtPoint)
