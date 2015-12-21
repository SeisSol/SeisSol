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

import DB
import Sparse
import MDS
import itertools
import operator
import numpy

class Operation:
  MEMSET = 1,
  GEMM = 2
  
class Kernel(object):
  def __init__(self, kernel, db, architecture):
    self.db = db
    self.arch = architecture
    self.involvedMatrices = set()
    self.temps = list()
    self.operations = list()
    self.tempBaseName = 'temporaryResult'
    self.resultName = 'result'
    self.kernel = kernel

    for mul in kernel.symbol:
      self.gemms([self.db[name] for name in mul])

    self.temps.sort(key=lambda temp: temp.name)
    self.involvedMatrices = set([name for mul in kernel.symbol for name in mul])
    self.involvedMatrices = sorted(list(self.involvedMatrices))
    self.involvedMatrices.append(self.resultName)
      
  def gemms(self, matrices):
    raise NotImplementedError()

class GeneratedKernel(Kernel):
  def __init__(self, kernel, db, architecture):
    self.gemmlist = list()
    self.nonZeroFlops = 0
    self.hardwareFlops = 0
    
    super(GeneratedKernel, self).__init__(kernel, db, architecture)

  def __intersect(self, blockA, blockB):
    return (max(blockA.startcol, blockB.startrow), min(blockA.stopcol, blockB.stoprow))

  def gemms(self, matrices):
    nameToIndex = dict()
    for i,matrix in enumerate(matrices):
      nameToIndex[matrix.name] = i

    sparsityPatterns = [matrix.spp for matrix in matrices]
    # Eliminate irrelevant entries in the matrix multiplication
    equivalentSparsityPatterns = Sparse.equivalentMultiplicationPatterns(sparsityPatterns)
    # Fit the equivalent sparsity pattern tightly for each memory block
    fittedBlocks = [DB.patternFittedBlocks(matrix.blocks, equivalentSparsityPatterns[i]) for i, matrix in enumerate(matrices)]
    # Find the actual implementation pattern, i.e. dense matrix -> dense pattern
    implementationPatterns = [matrix.getImplementationPattern(fittedBlocks[i], equivalentSparsityPatterns[i]) for i, matrix in enumerate(matrices)]
    # Determine the matrix multiplication order based on the implementation pattern
    chainOrder, dummy = Sparse.sparseMatrixChainOrder(implementationPatterns)
    
    self.nonZeroFlops += Sparse.calculateOptimalSparseFlops(equivalentSparsityPatterns)

    # convert matrix chain order to postfix
    stack = list()
    output = list()
    stack.append((0, len(matrices)-1))
    while len(stack) > 0:
      current = stack[-1]
      i = current[0]
      j = current[1]
      stack.pop()
      if (i == j): # matrix chain is a single matrix
        output.append(matrices[i]) # post the matrix A_i
      else: # subproblem A_i * ... * A_j
        output.append('*') # post a multiplication
        k = chainOrder[current[0], current[1]] # split position A_i..k * A_(k+1)j
        stack.append((current[0], k))   # post A_i..k
        stack.append((k+1, current[1])) # post A_(k+1)j

    # parse postfix
    operands = list()
    tempCounter = len(self.temps)
    while len(output) > 0:
      top = output.pop()
      if top != '*':
        operands.append(top)
      else:
        # In the following we generate instructions for op1 * op2.
        op2 = operands.pop()
        op1 = operands.pop()
        
        blocks1 = fittedBlocks[nameToIndex[op1.name]] if nameToIndex.has_key(op1.name) else op1.blocks
        blocks2 = fittedBlocks[nameToIndex[op2.name]] if nameToIndex.has_key(op2.name) else op2.blocks
        
        # Manage temporary variables
        if len(output) > 0:
          if len(self.temps) == 0:
            tempCounter += 1
            self.temps.append(DB.MatrixInfo(self.tempBaseName + str(tempCounter)))
          result = self.temps.pop()
          resultRequiredReals = result.requiredReals
          spp1 = implementationPatterns[nameToIndex[op1.name]] if nameToIndex.has_key(op1.name) else op1.spp
          spp2 = implementationPatterns[nameToIndex[op2.name]] if nameToIndex.has_key(op2.name) else op2.spp
          result = DB.MatrixInfo(result.name, op1.rows, op2.cols, sparsityPattern = spp1 * spp2)
          result.fitBlocksToSparsityPattern()
          result.generateMemoryLayout(self.arch, alignStartrow=True)
          resultName = result.name
          operands.append(result)
          beta = 0
          result.requiredReals = max(resultRequiredReals, result.requiredReals)
        else:
          beta = 1
          result = self.kernel
          resultName = self.resultName
        
        ops = []
        writes = []
        # op1 and op2 may be partitioned in several blocks.
        # Here we split the blocks of op1 and op2 in sums, i.e.
        # op1 * op2 = (op11 + op12 + ... + op1m) * (op21 + op22 + ... + op2n)
        #           = op11*op21 + op11*op22 + ...
        # opij is a matrix where the j-th block of opi is nonzero.
        # E.g. the first block of op1 (e.g. a 3x3 matrix) is given by (1, 2, 0, 2), then
        #        0 0 0
        # op11 = x x 0
        #        0 0 0
        for i1, block1 in enumerate(blocks1):
          for i2, block2 in enumerate(blocks2):
            # op1k * op2l is only nonzero if the columns of op1k and
            # the rows of op2l intersect.
            self.__gemm(op1.name, block1, op1.blocks[i1], op2.name, block2, op2.blocks[i2], resultName, result.blocks[0], beta, ops, writes)

        # Reorder ops in order to find betas
        if len(writes) > 0 and beta == 0:
          targetCard = result.blocks[0].ld * result.blocks[0].cols()
          mdsIn = MDS.maxDisjointSet(writes, targetCard)
          mdsOut = list( set(range(len(writes))).difference(set(mdsIn)) )
          order = mdsIn + mdsOut          
          memsetInterval = set(range(targetCard))
          for m in mdsIn:
            memsetInterval.difference_update(set( [i + j*result.blocks[0].ld for j in range(writes[m].startcol, writes[m].stopcol) for i in range(writes[m].startrow, writes[m].stoprow)] ))
          memsetInterval = list(memsetInterval)
          ranges = []
          for key, group in itertools.groupby(enumerate(memsetInterval), lambda(index, value): value - index):
            group = map(operator.itemgetter(1), group)
            start = group[0]
            end = group[-1] if len(group) > 1 else start
            self.operations.append(dict(
              type=Operation.MEMSET,
              pointer=resultName,
              offset=start,
              numberOfReals=1+end-start,
              dataType=self.arch.typename
            ))
          for m in mdsOut:
            ops[m]['gemm']['beta'] = 1
          ops = [ops[o] for o in order]

        for op in ops:
          self.gemmlist.append(op['gemm'])
          self.operations.append(op)
          if op['gemm']['spp'] is not None:
            NNZ = int(numpy.sum(op['gemm']['spp']))
            if op['gemm']['LDA'] < 1:
              self.hardwareFlops += 2 * NNZ * op['gemm']['N']
            else:
              self.hardwareFlops += 2 * op['gemm']['M'] * NNZ
          else:
            self.hardwareFlops += 2 * op['gemm']['M'] * op['gemm']['N'] * op['gemm']['K']
            
        # if op is temporary
        if not nameToIndex.has_key(op1.name):
          self.temps.append(op1)
        if not nameToIndex.has_key(op2.name):
          self.temps.append(op2)

  def __gemm(self, nameA, blockA, memoryBlockA, nameB, blockB, memoryBlockB, nameC, memoryBlockC, beta, ops, writes):
    if memoryBlockA.sparse and memoryBlockB.sparse:
      raise NotImplementedError('The generator does not support sparse * sparse multiplications.')

    k1, k2 = self.__intersect(blockA, blockB)
    if k2 > k1:
      if memoryBlockA.sparse:
        m1 = memoryBlockA.startrow
        m2 = memoryBlockA.stoprow
        sparsityPattern = memoryBlockA.sparsityPattern(k1, k2)
        mmType = 'sparse'
        spMtxName = nameA
      elif memoryBlockB.sparse:
        m1 = self.arch.getAlignedIndex(blockA.startrow)
        m2 = m1 + self.arch.getAlignedDim(blockA.stoprow - m1)
        k1 = memoryBlockB.startrow
        k2 = memoryBlockB.stoprow
        sparsityPattern = memoryBlockB.sparsityPattern(blockB.startcol, blockB.stopcol)
        mmType = 'sparse'
        spMtxName = nameB
      else:
        m1 = self.arch.getAlignedIndex(blockA.startrow)
        m2 = m1 + self.arch.getAlignedDim(blockA.stoprow - m1)
        sparsityPattern = None
        mmType = 'dense'
        spMtxName = ''
      
      offsetA = memoryBlockA.calcOffset(m1, k1)
      offsetB = memoryBlockB.calcOffset(k1, blockB.startcol)

      M = m2-m1
      N = blockB.cols()
      K = k2-k1
      offsetC = memoryBlockC.calcOffset(m1, blockB.startcol)

      startrowC = m1 - memoryBlockC.startrow
      startcolC = blockB.startcol - memoryBlockC.startcol
      writes.append(DB.MatrixBlock(startrowC, startrowC + M, startcolC, startcolC + N))

      alignedA = self.arch.checkAlignment(offsetA)
      alignedC = self.arch.checkAlignment(offsetC)
      
      if not memoryBlockA.sparse and (not alignedA or not alignedC):
        print('WARNING for {}*{}: Either A or C cannot be aligned.'.format(nameA, nameB))
    
      gemm = {
        'type':         mmType,
        'M':            M,
        'N':            N,
        'K':            K,
        'LDA':          memoryBlockA.ld,
        'LDB':          memoryBlockB.ld,
        'LDC':          memoryBlockC.ld,
        'alpha':        1,
        'beta':         beta,
        'alignedA':     int(alignedA),
        'alignedC':     int(alignedC),
        'spp':          sparsityPattern,
        'spMtxName':    spMtxName
      }
      ops.append(dict(
        type=Operation.GEMM,
        gemm=gemm,
        nameA=nameA,
        nameB=nameB,
        nameC=nameC,
        offsetA=offsetA,
        offsetB=offsetB,
        offsetC=offsetC        
      ))
    

class ReferenceKernel(Kernel):
  def __init__(self, kernel, db, architecture):
    super(ReferenceKernel, self).__init__(kernel, db, architecture)

  def gemms(self, matrices):
    resultSize = 0
    leftName = matrices[0].name
    leftRows = matrices[0].rows
    leftCols = matrices[0].cols
    tempCounter = len(self.temps)
    for i in range(1, len(matrices)):
      if i < len(matrices)-1:
        tempCounter += 1
        tempName = self.tempBaseName + str(tempCounter)
        temp = DB.MatrixInfo(tempName, leftRows, matrices[i].cols)
        self.temps.append(temp)
        beta = 0.0
      else:
        tempName = 'result'
        beta = 1.0
      self.operations.append(dict(
        type=Operation.GEMM,
        m=leftRows,
        n=matrices[i].cols,
        k=leftCols,
        A=leftName,
        B=matrices[i].name,
        beta=beta,
        C=tempName
      ))
      leftName = tempName
      leftCols = matrices[i].cols
