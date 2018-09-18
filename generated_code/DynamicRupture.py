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

from yateto import *
from yateto.input import parseJSONMatrixFile

def addKernels(generator, Q, matricesDir, order, dynamicRuptureMethod, numberOfElasticQuantities, numberOfQuantities):
  numberOfBasisFunctions = order*(order+1)*(order+2)//6

  if dynamicRuptureMethod == 'quadrature':
    numberOfPoints = (order+1)**2
  elif dynamicRuptureMethod == 'cellaverage':
    numberOfPoints = int(4**math.ceil(math.log(order*(order+1)/2,4)))
  else:
    raise ValueError('Unknown dynamic rupture method.')

  clones = dict()

  # Load matrices
  db = parseJSONMatrixFile('{}/dr_{}_matrices_{}.json'.format(matricesDir, dynamicRuptureMethod, order), clones)

  # Determine matrices  
  # Note: This does only work because the flux does not depend on the mechanisms in the case of viscoelastic attenuation
  godunovMatrix = Tensor('godunovMatrix', (numberOfElasticQuantities, numberOfElasticQuantities))
  fluxSolver    = Tensor('fluxSolver', (numberOfElasticQuantities, numberOfQuantities))
  godunovState  = Tensor('godunovState', (numberOfPoints, numberOfElasticQuantities), alignStride=True)

  def godunovStateGenerator(i, h):
    if h == 0:
      return Q['kp'] <= godunovState['kp'] <= db.nP[i]['kl'] * Q['lq'] * godunovMatrix['qp']
    return Q['kp'] <= godunovState['kp'] <= godunovState['kp'] + db.nM[i+4*(h-1)]['kl'] * Q['lq'] * godunovMatrix['qp']

  generator.addFamily('godunovState', simpleParameterSpace(4,4), godunovStateGenerator)

  def nodalFluxGenerator(i, h):
    if h == 0:
      return Q['kp'] <= db.pP[i]['kl'] * godunovState['lq'] * fluxSolver['qp']
    return Q['kp'] <= db.pM[i+4*(h-1)]['kl'] * godunovState['lq'] * fluxSolver['qp']
  
  generator.addFamily('nodalFlux', simpleParameterSpace(4,4), nodalFluxGenerator)
