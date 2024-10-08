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
  
from yateto.type import Tensor
from yateto.ast.node import IndexedTensor
from yateto.memory import DenseMemoryLayout

class OptionalDimTensor(Tensor):
  # dimSize = 1 is considered optional
  def __init__(self, name, optName, optSize, optPos, shape, spp=None, memoryLayoutClass=DenseMemoryLayout, alignStride=False):
    self._optName = optName
    self._optSize = optSize
    self._optPos = optPos
    shape = self.insertOptDim(shape, (self._optSize,))
    super().__init__(name, shape, spp, memoryLayoutClass, alignStride)

  def hasOptDim(self):
    return self._optSize > 1

  def insertOptDim(self, sliceable, item):
    if self.hasOptDim():
      return sliceable[0:self._optPos] + item + sliceable[self._optPos:]
    return sliceable

  def __getitem__(self, indexNames):
    indexNames = self.insertOptDim(indexNames, self._optName)
    return IndexedTensor(self, indexNames)

  def optName(self):
    return self._optName

  def optSize(self):
    return self._optSize

  def optPos(self):
    return self._optPos
