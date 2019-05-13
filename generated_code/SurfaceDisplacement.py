#!/usr/bin/env python
##
# @file
# This file is part of SeisSol.
#
# @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
#
# @section LICENSE
# Copyright (c) 2017-2019, SeisSol Group
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
from yateto import Tensor, simpleParameterSpace
from yateto.memory import CSCMemoryLayout
from multSim import OptionalDimTensor

def addKernels(generator, aderdg):
  maxDepth = 3

  numberOf3DBasisFunctions = aderdg.numberOf3DBasisFunctions()
  numberOfQuantities = aderdg.numberOfQuantities()
  selectVelocitySpp = np.zeros((numberOfQuantities, 3))
  selectVelocitySpp[6:9,0:3] = np.eye(3)
  selectVelocity = Tensor('selectVelocity', selectVelocitySpp.shape, selectVelocitySpp, CSCMemoryLayout)

  displacement = OptionalDimTensor('displacement', aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), (numberOf3DBasisFunctions, 3), alignStride=True)
  generator.add('addVelocity', displacement['kp'] <= displacement['kp'] + aderdg.I['kq'] * selectVelocity['qp'])

  subTriangleDofs = [OptionalDimTensor('subTriangleDofs({})'.format(depth), aderdg.Q.optName(), aderdg.Q.optSize(), aderdg.Q.optPos(), (4**depth, 3), alignStride=True) for depth in range(maxDepth+1)]
  subTriangleProjection = [Tensor('subTriangleProjection({})'.format(depth), (4**depth, numberOf3DBasisFunctions), alignStride=True) for depth in range(maxDepth+1)]

  subTriangleDisplacement = lambda depth: subTriangleDofs[depth]['kp'] <= subTriangleProjection[depth]['kl'] * displacement['lp']
  subTriangleVelocity     = lambda depth: subTriangleDofs[depth]['kp'] <= subTriangleProjection[depth]['kl'] * aderdg.Q['lq'] * selectVelocity['qp']

  generator.addFamily('subTriangleDisplacement', simpleParameterSpace(maxDepth+1), subTriangleDisplacement)
  generator.addFamily('subTriangleVelocity', simpleParameterSpace(maxDepth+1), subTriangleVelocity)
