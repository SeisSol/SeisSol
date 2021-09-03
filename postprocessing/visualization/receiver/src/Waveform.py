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

import numpy

class Waveform:
  def __init__(self, names, data, coordinates):
    data = numpy.array(data)
    
    self.waveforms = dict()
    self.norms = dict()
    self.printIndicator = dict()
    for i in range(0, len(names)):
      if names[i] == 'Time':
        self.time = data[:,i]
      else:
        self.waveforms[ names[i] ] = data[:,i]
        self.norms[ names[i] ] = numpy.max(numpy.abs(data[:,i]))
        self.printIndicator[ names[i] ] = True
    
    self.coordinates = numpy.array(coordinates)

  def subtract(self, other):
    newTime = numpy.union1d(self.time, other.time)
    self.waveforms = { key: value for key,value in self.waveforms.items() if key in other.waveforms }
    for name, wf in self.waveforms.items():
      wf0 = numpy.interp(newTime, other.time, other.waveforms[name])
      wf1 = numpy.interp(newTime, self.time, self.waveforms[name])
      self.waveforms[name] = wf1 - wf0
    self.time = newTime

  def normalize(self):
    for name, wf in self.waveforms.items():
      if self.norms[name] > numpy.finfo(float).eps:
        self.waveforms[name] = wf / self.norms[name]

