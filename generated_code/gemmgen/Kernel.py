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
    
    self.__traverseExprTree(kernel.symbol)

    self.temps.sort(key=lambda temp: temp.name)
    self.involvedMatrices = sorted(list(self.involvedMatrices))
    self.involvedMatrices.append(self.resultName)
    
  def __traverseExprTree(self, expr):
    if expr.is_MatAdd:
      for arg in expr.args:
          self.__traverseExprTree(arg)
    elif expr.is_MatMul:
      matrices = list()
      for symbol in expr.args:
        matrices.append(self.db[symbol.name])
        self.involvedMatrices.update([symbol.name])
      self.gemms(matrices)
    else:
      raise ValueError('Unexpected subexpression {}.'.format(expr))
      
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
    implementationPatterns = [DB.getImplementationPattern(fittedBlocks[i], equivalentSparsityPatterns[i]) for i in range(0, len(matrices))]
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
          result = DB.MatrixInfo(result.name, op1.rows, op2.cols, sparsityPattern = op1.spp * op2.spp)
          result.fitBlocksToSparsityPattern()
          result.generateMemoryLayout(self.arch, alignStartrow=True)
          resultName = result.name
          operands.append(result)
          # check how many gemm-calls will be generated
          gemmCallCount = 0
          for i1, block1 in enumerate(blocks1):
            for i2, block2 in enumerate(blocks2):
              n1, n2 = self.__intersect(block1, block2)
              if n2 > n1:
                gemmCallCount += 1              
          if gemmCallCount > 1:
            beta = 1
            self.operations.append(dict(
              type=Operation.MEMSET,
              pointer=resultName,
              numberOfReals=result.requiredReals,
              dataType=self.arch.typename
            ))
          else:
            beta = 0
          result.requiredReals = max(resultRequiredReals, result.requiredReals)
        else:
          beta = 1
          result = self.kernel
          resultName = self.resultName
        
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
            self.__gemm(op1.name, block1, op1.blocks[i1], op2.name, block2, op2.blocks[i2], resultName, result.blocks[0], beta)
            
        # if op is temporary
        if not nameToIndex.has_key(op1.name):
          self.temps.append(op1)
        if not nameToIndex.has_key(op2.name):
          self.temps.append(op2)

  def __gemm(self, nameA, blockA, memoryBlockA, nameB, blockB, memoryBlockB, nameC, memoryBlockC, beta):
    n1, n2 = self.__intersect(blockA, blockB)
    if n2 > n1:
      m1 = self.arch.getAlignedIndex(blockA.startrow)
      m2 = m1 + self.arch.getAlignedDim(blockA.stoprow - m1)
      offsetA = memoryBlockA.offset + (m1 - memoryBlockA.startrow) + memoryBlockA.ld * (n1 - memoryBlockA.startcol)
      offsetB = memoryBlockB.offset + (n1 - memoryBlockB.startrow) + memoryBlockB.ld * (blockB.startcol - memoryBlockB.startcol)
      offsetC = memoryBlockC.offset + (m1 - memoryBlockC.startrow) + memoryBlockC.ld * (blockB.startcol - memoryBlockC.startcol)
      alignedA = self.arch.checkAlignment(offsetA)
      alignedC = self.arch.checkAlignment(offsetC)
      
      if not alignedA or not alignedC:
        print('WARNING for {}*{}: Either A or C cannot be aligned.'.format(nameA, nameB))
    
      gemm = {
        'type':         'dense',
        'M':            m2-m1,
        'N':            blockB.cols(),
        'K':            n2-n1,
        'LDA':          memoryBlockA.ld,
        'LDB':          memoryBlockB.ld,
        'LDC':          memoryBlockC.ld,
        'alpha':        1,
        'beta':         beta,
        'alignedA':     int(alignedA),
        'alignedC':     int(alignedC)
      }
      self.gemmlist.append(gemm)
      self.operations.append(dict(
        type=Operation.GEMM,
        gemm=gemm,
        nameA=nameA,
        nameB=nameB,
        nameC=nameC,
        offsetA=offsetA,
        offsetB=offsetB,
        offsetC=offsetC        
      ))
      self.hardwareFlops += 2 * gemm['M'] * gemm['N'] * gemm['K']
    

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
