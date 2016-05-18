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

import numpy as np
import copy
import sys

def equivalentMultiplicationPatterns(sparsityPatterns):
  layers = list()
  for pattern in sparsityPatterns:
    layers.append(np.zeros(pattern.shape[0], dtype=np.int))
  layers.append(np.zeros(sparsityPatterns[-1].shape[1], dtype=np.int))
  
  sweepedPatterns = copy.deepcopy(sparsityPatterns)

  for r, pattern in enumerate(sparsityPatterns):
    for i in range(0, pattern.shape[0]):
      if r == 0 or layers[r][i] != 0:
        for j in range(0, pattern.shape[1]):
          if pattern[i,j] != 0.0:
            layers[r+1][j] += 1 # incoming
      else:
        sweepedPatterns[r][i,:] = 0.0
  
  for layer in layers:
    layer[:] = 0
    
  for r, pattern in enumerate(list(reversed(sparsityPatterns))):
    R = len(layers)-1-r
    for i in range(0, pattern.shape[1]):
      if r == 0 or layers[R][i] != 0:
        for j in range(0, pattern.shape[0]):
          if pattern[j,i] != 0.0:
            layers[R-1][j] += 1 # incoming
      else:
        sweepedPatterns[R-1][:,i] = 0.0

  return sweepedPatterns


def calculateOptimalSparseFlops(equivalentSparsityPatterns):
  s, m = sparseMatrixChainOrder(equivalentSparsityPatterns)
  return 2 * int(m[0, -1]) # additions and multiplications

  
def sparseMatrixChainOrder(sparsityPatterns):
  p = list()
  p.append(sparsityPatterns[0].shape[0])
  for pattern in sparsityPatterns:
    p.append(pattern.shape[1])
  
  n = len(sparsityPatterns)
  m = np.matlib.zeros((n, n))
  s = np.matlib.zeros((n, n), dtype=int)
  products = [[0 for j in range(0,n)] for i in range(0,n)]
  for i,matrix in enumerate(sparsityPatterns):
    products[i][i] = matrix
  for L in range(1,n):
    for i in range(0, n-L):
      j = i + L
      m[i, j] = sys.maxint
      for k in range(i, j):
        product = products[i][k]*products[k+1][j]
        q = m[i,k] + m[k+1,j] + np.sum(product)
        if q <= m[i,j]:
          m[i,j] = q
          s[i,j] = k
          product[np.abs(product) > 0] = 1.0
          products[i][j] = product

  return (s, m)
