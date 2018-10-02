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

from yateto import *
from yateto.input import parseXMLMatrixFile

def addKernels(generator, qi, qShape, alignStride, matricesDir, order, PlasticityMethod):
  # Load matrices
  db = parseXMLMatrixFile('{}/plasticity_{}_matrices_{}.xml'.format(matricesDir, PlasticityMethod, order), clones=dict(), alignStride=alignStride)
  numberOfNodes = db.v.shape()[0]

  qPos = qi('lq').find('q')
  sShape = tuple(6 if i == qPos else q for i,q in enumerate(qShape))
  stressDOFS  = Tensor('stressDOFS', sShape, alignStride=alignStride)

  bPos = qi('lq').find('l')
  iShape = tuple(numberOfNodes if i == bPos else q for i,q in enumerate(sShape))
  interpolationDOFS = Tensor('interpolationDOFS', iShape, alignStride=alignStride)
  
  generator.add('interpolationDOFS', interpolationDOFS[qi('kp')] <= db.v['kl'] * stressDOFS[qi('lp')])
  generator.add('convertToModal', stressDOFS[qi('kp')] <= db.vInv['kl'] * interpolationDOFS[qi('lp')])
