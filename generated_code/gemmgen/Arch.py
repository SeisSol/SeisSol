##
# @file
# This file is part of SeisSol.
#
# @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
#
# @section LICENSE
# Copyright (c) 2015, SeisSol Group
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

alignments = {
  'noarch': 16,
  'wsm':    16,
  'snb':    32,
  'hsw':    32,
  'skx':    64,
  'knc':    64,
  'knl':    64
}

# Libxsmm currently supports prefetch only for KNL kernels
enablePrefetch = {
  'noarch': False,
  'wsm':    False,
  'snb':    False,
  'hsw':    False,
  'skx':    False,
  'knc':    False,
  'knl':    True
}

class Architecture(object):
  def __init__(self, name, precision, alignment, enablePrefetch=False):
    self.name = name
    self.precision = precision.upper()
    if self.precision == 'D':
      self.bytesPerReal = 8
      self.typename = 'double'
      self.epsilon = 2.22e-16
    elif self.precision == 'S':
      self.bytesPerReal = 4
      self.typename = 'float'
      self.epsilon = 1.19e-7
    else:
      raise ValueError('Unknown precision type ' + self.precision)
    self.alignment = alignment
    self.alignedReals = self.alignment / self.bytesPerReal
    self.enablePrefetch = enablePrefetch
    
  def getAlignedIndex(self, index):
    return index - index % self.alignedReals

  def getAlignedDim(self, dim):
    return dim + (self.alignedReals - dim % self.alignedReals) % self.alignedReals

  def checkAlignment(self, offset):
    return offset % self.alignedReals == 0

def getArchitectureByIdentifier(ident):
  precision = ident[0].upper()
  name = ident[1:]
  return Architecture(name, precision, alignments[name], enablePrefetch[name])
