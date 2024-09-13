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
import math

def median(values):
  n = len(values)
  x = sorted(values)
  if n%2 == 0:
    return 0.5 * (x[n/2-1] + x[n/2])
  return x[n/2]

def mean(values):
  n = float(len(values))
  return sum(values) / n
  
def var(values):
  if len(values) < 2:
    return float('nan')
  m = mean(values)
  n = float(len(values))
  return sum([(v-m)**2. for v in values]) / (n-1.)
  
def stdev(values):
  return math.sqrt(var(values))

def writeTimes(times):
  with open('times.txt', 'w') as f:
    f.write('{:16} {:10} {:10} {:10}\n'.format('Name', 'Median', 'Mean', 'Std. dev.'))
    for name, value in sorted(times.iteritems()):
      for idx, time in sorted(value.iteritems()):
        f.write('{:16} {:>10f} {:>10f} {:>10f}\n'.format(name + str(idx), median(time), mean(time), stdev(time)))

def analyse():
  times = dict()
  timePattern = re.compile(r'^time\D+([0-9\.]+)', re.MULTILINE)
  for resultFile in os.listdir(Proxy.OutputDir):
    resultFileName, extension = os.path.splitext(resultFile)
    if extension == '.run':
      matrix = resultFileName.split('_')[0]
      content = open(os.path.join(Proxy.OutputDir, resultFile), 'r').read()
      timeSearch = timePattern.search(content)
      if timeSearch:
        time = float(timeSearch.group(1))
        matrixBase = matrix[:-1]
        matrixVariant = int(matrix[-1])
	if matrixVariant > 0:
		continue
        if not times.has_key(matrixBase):
          times[matrixBase] = dict()
        if not times[matrixBase].has_key(matrixVariant):
          times[matrixBase][matrixVariant] = list()
        times[matrixBase][matrixVariant].append(time)
      else:
        print('Warning: Invalid result file {}.'.format(resultFileName))
        
  writeTimes(times)
  
  matrices = list()
  dense = times.pop('dense').pop(0)
  denseTime = median(dense)
  timeSaveVsDense = 0.
  for key, variants in times.iteritems():
    matrixTimes = {variant: median(timeSeries) for variant, timeSeries in variants.iteritems()}
    minTimeKey = min(matrixTimes, key=matrixTimes.get)
    if matrixTimes[minTimeKey] < denseTime:
      matrices.append((key, minTimeKey))
      timeSaveVsDense += denseTime - matrixTimes[minTimeKey]
      
  print('Estimated tuned time: {}'.format(denseTime - timeSaveVsDense))
  return matrices
  

