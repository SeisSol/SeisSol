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

def addKernels(g, db, Q, order, numberOfQuantities, numberOfExtendedQuantities):
  numberOf3DBasisFunctions = order*(order+1)*(order+2)//6

  ti = Collection()
  ti.AplusT = Tensor('AplusT', (numberOfQuantities, numberOfExtendedQuantities))
  ti.AminusT = Tensor('AminusT', (numberOfQuantities, numberOfExtendedQuantities))
  T = Tensor('T', (numberOfExtendedQuantities, numberOfExtendedQuantities))
  Tinv = Tensor('Tinv', (numberOfQuantities, numberOfQuantities))
  QgodLocal = Tensor('QgodLocal', (numberOfQuantities, numberOfQuantities))
  QgodNeighbor = Tensor('QgodNeighbor', (numberOfQuantities, numberOfQuantities))
  QFortran = Tensor('QFortran', (numberOf3DBasisFunctions, numberOfQuantities))

  ti.oneSimToMultSim = Tensor('oneSimToMultSim', (Q.optSize(),), spp={(i,): '1.0' for i in range(Q.optSize())})
  multSimToFirstSim = Tensor('multSimToFirstSim', (Q.optSize(),), spp={(0,): '1.0'})

  fluxScale = Scalar('fluxScale')
  computeFluxSolverLocal = ti.AplusT['ij'] <= fluxScale * Tinv['ki'] * QgodLocal['kq'] * db.star[0]['ql'] * T['jl']
  g.add('computeFluxSolverLocal', computeFluxSolverLocal)

  computeFluxSolverNeighbor = ti.AminusT['ij'] <= fluxScale * Tinv['ki'] * QgodNeighbor['kq'] * db.star[0]['ql'] * T['jl']
  g.add('computeFluxSolverNeighbor', computeFluxSolverNeighbor)

  if Q.hasOptDim():
    addQFortran = Q['kp'] <= Q['kp'] + QFortran['kp'] * ti.oneSimToMultSim['s']
    copyQToQFortran = QFortran['kp'] <= Q['kp'] * multSimToFirstSim['s']
  else:
    addQFortran = Q['kp'] <= Q['kp'] + QFortran['kp']
    copyQToQFortran = QFortran['kp'] <= Q['kp']

  g.add('addQFortran', addQFortran)
  g.add('copyQToQFortran', copyQToQFortran)

  return ti
