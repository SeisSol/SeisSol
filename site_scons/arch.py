#! /usr/bin/env python
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

def getArchitectures():
  # wsm = Westmere
  # snb = Sandy Bridge
  # knc = Knights Corner (Xeon Phi)
  # hsw = Haswell
  # knl = Knight Landing (Xeon Phi)
  cpus = ['noarch', 'wsm', 'snb', 'knc', 'hsw', 'knl']
  precisions = ['s', 'd']
  return [p + c for c in cpus for p in precisions]
  
def getCpu(architecture):
  return architecture[1:]
  
def getAlignment(architecture):
  alignments = {
      'noarch': 16,
      'wsm': 16,
      'snb': 32,
      'hsw': 32,
      'knc': 64,
      'knl': 64
  }
  return alignments[ getCpu(architecture) ]
  
def getFlags(architecture, compiler):
  if architecture not in getArchitectures():
    raise ValueError('Unknown architecture.')
  
  cpu = getCpu(architecture)
  
  if cpu == 'wsm':
    flags = ['-msse3']
  elif cpu == 'snb':
    flags =  ['-mavx']
  elif cpu == 'hsw':
    if compiler == 'intel':
      flags = ['-xCORE-AVX2', '-fma']
    else:
      flags = ['-mavx2', '-mfma']
  elif cpu == 'knc':
    flags = ['-mmic', '-fma']
  elif cpu == 'knl':
    if compiler == 'intel':
      flags = ['-xMIC-AVX512', '-fma', '-DENABLE_MATRIX_PREFETCH', '-DENABLE_STREAM_MATRIX_PREFETCH', '-DNUMBER_OF_THREADS_PER_GLOBALDATA_COPY=4']
    else:
      flags = ['-mavx512f', '-mavx512cd', '-mavx512pf', '-mavx512er', '-mfma', '-DENABLE_MATRIX_PREFETCH', '-DENABLE_STREAM_MATRIX_PREFETCH', '-DNUMBER_OF_THREADS_PER_GLOBALDATA_COPY=4']
  else:
    flags = []
  
  # enable interproc. opts for small cores
  if cpu in ['knc', 'knl']:
    flags.extend(['-ip', '-ipo'])
              
  return flags
  
