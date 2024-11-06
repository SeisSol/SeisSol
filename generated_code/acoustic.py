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
from yateto.input import parseXMLMatrixFile, memoryLayoutFromFile

from aderdg import LinearADERDG

class AcousticADERDG(LinearADERDG):
  def __init__(self, order, multipleSimulations, matricesDir, memLayout, **kwargs):
    super().__init__(order, multipleSimulations, matricesDir)
    clones = {
      'star': ['star(0)', 'star(1)', 'star(2)'],
    }
    self.db.update(
      parseXMLMatrixFile('{}/star_acoustic.xml'.format(matricesDir), clones)
    )

    memoryLayoutFromFile(memLayout, self.db, clones)
    self.kwargs = kwargs

  def numberOfQuantities(self):
    return 4

  def starMatrix(self, dim):
    return self.db.star[dim]

  def addLocal(self, generator, targets):
    super().addLocal(generator, targets)
  
  def extractVelocities(self):
    extractVelocitiesSPP = np.zeros((3, self.numberOfQuantities()))
    extractVelocitiesSPP[0, 1] = 1
    extractVelocitiesSPP[1, 2] = 1
    extractVelocitiesSPP[2, 3] = 1
    return extractVelocitiesSPP

  def extractTractions(self):
    extractTractionsSPP = np.zeros((3, self.numberOfQuantities()))
    extractTractionsSPP[0, 0] = 1
    return extractTractionsSPP

  def setBottomToIdentity(self, selectVelocitySpp):
    selectVelocitySpp[1:4,0:3] = np.eye(3)

  def setSelectTractionSppToTrue(self, selectTractionSpp):
    selectTractionSpp[0,0] = True
