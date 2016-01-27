#!/usr/bin/env python
##
# @file
# This file is part of SeisSol.
#
# @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
#
# @section LICENSE
# Copyright (c) 2016, SeisSol Group
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

import Proxy
import os
import re

def writeTimes(denseTime, times):
  with open('times.txt', 'w') as f:
    f.write('dense {}\n'.format(denseTime))
    for name, value in sorted(times.iteritems()):
      for idx, time in sorted(value.iteritems()):
        f.write('{}{} {}\n'.format(name, idx, time))

def analyse():
  denseTime = float('nan')
  times = dict()
  timePattern = re.compile(r'^time\D+([0-9\.]+)', re.MULTILINE)
  for resultFile in os.listdir(Proxy.OutputDir):
    matrix, extension = os.path.splitext(resultFile)
    if extension == '.run':
      content = open(os.path.join(Proxy.OutputDir, resultFile), 'r').read()
      time = timePattern.search(content).group(1)
      if matrix == 'dense':
        denseTime = time
      else:
        matrixBase = matrix[:-1]
        matrixVariant = matrix[-1]
        if not times.has_key(matrixBase):
          times[matrixBase] = dict()
        times[matrixBase][matrixVariant] = time
        
  writeTimes(denseTime, times)
  
  matrices = list()
  for key, value in times.iteritems():
    minTimeKey = min(value, key=value.get)
    if value[minTimeKey] < denseTime:
      matrices.append((key, int(minTimeKey)))
      
  return matrices
  

