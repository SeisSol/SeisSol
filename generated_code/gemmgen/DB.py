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
import numpy.matlib
import scipy.sparse
import sympy
import copy
import re

def getImplementationPattern(fittedBlocks, pattern):
  implementationPattern = copy.deepcopy(pattern)
  for block in fittedBlocks:
    implementationPattern[block.slice()] = 1.0
  return implementationPattern

class MatrixBlock:
  def __init__(self, startrow, stoprow, startcol, stopcol):
    self.startrow = startrow
    self.stoprow = stoprow
    self.startcol = startcol
    self.stopcol = stopcol
    self.ld = -1
    self.offset = -1

  def rows(self):
    return self.stoprow - self.startrow

  def cols(self):
    return self.stopcol - self.startcol
    
  def slice(self):
    return (slice(self.startrow, self.stoprow), slice(self.startcol, self.stopcol))

  def __repr__(self):
    return '{{slice: [{}:{},{}:{}], ld: {}, offset: {}}}'.format(self.startrow, self.stoprow, self.startcol, self.stopcol, self.ld, self.offset)

def findLargestNonzeroBlock(pattern):
  nnz = pattern.nonzero()
  if nnz[0].shape[0] == 0 or nnz[0].shape[1] == 0:
    return MatrixBlock(0, 0, 0, 0)
  else:
    return MatrixBlock(nnz[0].min(), nnz[0].max()+1, nnz[1].min(), nnz[1].max()+1)

def patternFittedBlocks(blocks, pattern):
  fittedBlocks = list()
  for block in blocks:
    bounds = findLargestNonzeroBlock(pattern[block.slice()])
    bounds.startrow += block.startrow
    bounds.stoprow += block.startrow
    bounds.startcol += block.startcol
    bounds.stopcol += block.startcol
    fittedBlocks.append(bounds)
  return fittedBlocks
  
def determineGlobalMatrixIds(globalMatrixIdRules, db):
  for key, value in db.iteritems():
    for rule in globalMatrixIdRules:
      match = re.search(rule[0], key)
      if match != None:
        value.globalMatrixId = rule[1](match.groups())

class MatrixInfo:
  def __init__(self, name, rows = 0, cols = 0, sparsityPattern = None, values = None):
    self.name = name
    self.rows = rows
    self.cols = cols
    self.blocks = [MatrixBlock(0, self.rows, 0, self.cols)]
    self.requiredReals = -1
    self.leftMultiplication = False
    self.rightMultiplication = False
    self.symbol = sympy.MatrixSymbol(name, self.rows, self.cols)
    self.globalMatrixId = -1

    if isinstance(sparsityPattern, tuple):
      self.spp = scipy.sparse.coo_matrix((np.ones(len(sparsityPattern[0])), sparsityPattern), shape=(self.rows, self.cols), dtype=np.float64).todense()
    elif sparsityPattern is None:
      self.spp = np.matlib.ones((self.rows, self.cols), dtype=np.float64)
    else:
      self.spp = sparsityPattern
    
    # ensure that spp has only zeros and ones
    self.spp[np.abs(self.spp) > 0] = 1.0
      
    if self.spp.shape[0] != self.rows or self.spp.shape[1] != self.cols:
      raise ValueError('Matrix dimensions are different to the dimensions of the sparsity pattern.')
      
    self.values = values
    if self.values != None:
      self.values.sort(key=lambda entry: (entry[1], entry[0]))
    
  def __mul__(self, other):
    spp = self.spp * other.spp
    result = MatrixInfo('{}*{}'.format(self.name, other.name), self.rows, other.cols, spp)
    result.symbol = self.symbol * other.symbol
    return result
    
  def __add__(self, other):
    spp = self.spp + other.spp
    result = MatrixInfo('{}+{}'.format(self.name, other.name), self.rows, self.cols, spp)
    result.symbol = self.symbol + other.symbol
    return result
    
  def setBlocks(self, blocks):
    self.blocks = blocks
    self.__checkBlocks()
    
  def fitBlocksToSparsityPattern(self):
    self.setBlocks(patternFittedBlocks(self.blocks, self.spp))

  def __checkBlocks(self):
    # look for overlapping blocks
    # naive algorithm should do as I do not expect a large amount of blocks
    for i in range(0, len(self.blocks)):
      bi = self.blocks[i]
      if bi.startrow >= bi.stoprow or bi.startcol >= bi.stopcol or bi.startrow < 0 or bi.startcol < 0 or bi.stoprow > self.rows or bi.stopcol > self.cols:
        raise ValueError('{} ({}x{}): The memory block {} is invalid.'.format(self.name, self.rows, self.cols, bi))
      for j in range(i+1, len(self.blocks)):
        bj = self.blocks[j]
        if not (bj.startrow >= bi.stoprow or bj.stoprow <= bi.startrow or bj.startcol >= bi.stopcol or bj.stopcol <= bi.startcol):
          raise ValueError('{}: The memory blocks {} and {} overlap.'.format(self.name, bi, bj))
    
    # check if all nnzs are contained in the blocks
    nnz = np.sum(self.spp)
    for block in self.blocks:
      nnz -= np.sum(self.spp[block.slice()])
    
    if nnz != 0:
      raise ValueError('{}: The memory blocks do not cover all nonzeros.'.format(self.name))

  def generateMemoryLayout(self, architecture, alignStartrow=False):
    self.requiredReals = 0
    if alignStartrow == True:
      for block in self.blocks:
        if architecture.checkAlignment(block.startrow) == False:
          block.startrow = architecture.getAlignedIndex(block.startrow)
    for block in self.blocks:
      block.ld = architecture.getAlignedDim(block.rows()) if self.leftMultiplication or not self.rightMultiplication else block.rows()
      block.offset = self.requiredReals
      self.requiredReals += block.ld * block.cols()
    
  def flat(self, name):
    return MatrixInfo(name, self.rows, self.cols, sparsityPattern=self.spp)
    
  def getValuesAsStoredInMemory(self):
    if self.values != None:
      values = list()
      for block in self.blocks:
        blockValues = ['0.'] * (block.ld * block.cols())
        for entry in self.values:
          if entry[0] >= block.startrow and entry[0] < block.stoprow and entry[1] >= block.startcol and entry[1] < block.stopcol:
            blockValues[(entry[0] - block.startrow) + (entry[1] - block.startcol) * block.ld] = entry[2]
        values.extend(blockValues)
      return values
    else:
      return None
      
  def getIndexLUT(self):
    lut = [-1] * (self.rows * self.cols)
    for block in self.blocks:
      for col in range(block.startcol, block.stopcol):
        for row in range(block.startrow, block.stoprow):
          if self.spp[row, col] != 0.0:
            lut[row + col * self.rows] = block.offset + row + col * block.ld
    return lut
    
