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

from DB import *
from Sparse import equivalentMultiplicationPatterns
from sys import maxint
from os import system, path
from numpy import matlib, sum, abs
from itertools import groupby
from copy import deepcopy
import re
  
def sparseMatrixChainOrder(sparsityPatterns):
  p = list()
  p.append(sparsityPatterns[0].shape[0])
  for pattern in sparsityPatterns:
    p.append(pattern.shape[1])
  
  n = len(sparsityPatterns)
  m = matlib.zeros((n, n))
  s = matlib.zeros((n, n), dtype=int)
  products = [[0 for j in range(0,n)] for i in range(0,n)]
  for i,matrix in enumerate(sparsityPatterns):
    products[i][i] = matrix
  for L in range(1,n):
    for i in range(0, n-L):
      j = i + L
      m[i, j] = maxint
      for k in range(i, j):
        product = products[i][k]*products[k+1][j]
        q = m[i,k] + m[k+1,j] + sum(product)
        if q < m[i,j]:
          m[i,j] = q
          s[i,j] = k
          product[abs(product) > 0] = 1.0
          products[i][j] = product

  return s
    
def generateRoutineName(gemm):
  return 'dgemm_m{}_n{}_k{}_ldA{}_ldB{}_ldC{}_beta{}_alignedA{}_alignedC{}_pfsigonly'.format(
    gemm['M'],
    gemm['N'],
    gemm['K'],
    gemm['LDA'],
    gemm['LDB'],
    gemm['LDC'],
    gemm['beta'],
    gemm['alignedA'],
    gemm['alignedC']
  )

def formatOffset(name, offset):
  if offset == 0:
    return name
  else:
    return '&{}[{}]'.format(name, offset)

class Generator:
  def __init__(self, db, libxsmmGenerator, architecture):
    self.db = db
    self.gemmlist = list()
    self.tempBaseName = 'temporaryResult'
    self.libxsmmGenerator = libxsmmGenerator
    self.architecture = architecture
    self.__result = None
    self.resultName = 'result'
    
  def generateGemms(self, outputDir):
    cppFilename = outputDir + '/gemms.cpp'
    hFilename = outputDir + '/gemms.h'
    with open(cppFilename, 'w+') as cpp: # delete existing file and create an empty one
      cpp.write("""#ifndef NDEBUG
extern long long libxsmm_num_total_flops;
#endif\n\n""")
    with open(hFilename, 'w+') as header:
      self.__headerGuardBegin(header, 'GEMMS')
      uniquegemms = [a for a,b in groupby(sorted(self.gemmlist))]
      for gemm in uniquegemms:
        header.write('void {name}({tn} const* A, {tn} const* B, {tn}* C, {tn} const* A_prefetch, {tn} const* B_prefetch, {tn} const* C_prefetch);\n'.format(
          name=generateRoutineName(gemm),
          tn=self.architecture.typename
        ))
        system(self.generateLibxsmmGeneratorCall(cppFilename, gemm))
      self.__headerGuardEnd(header)
  
  def generateLibxsmmGeneratorCall(self, filename, gemm):
    return '{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} pfsigonly {}P'.format(
      self.libxsmmGenerator,
      gemm['type'],
      filename,
      generateRoutineName(gemm),
      gemm['M'],
      gemm['N'],
      gemm['K'],
      gemm['LDA'],
      gemm['LDB'],
      gemm['LDC'],
      gemm['alpha'],
      gemm['beta'],
      gemm['alignedA'],
      gemm['alignedC'],
      self.architecture.name,
      self.architecture.precision
    )
    
  def __functionSignature(self, names, writeNames=True):
    if writeNames:
      signature = ', '.join(['{} const* {}'.format(self.architecture.typename, name) for name in names[0:-1]])
      signature += ', {}* {}'.format(self.architecture.typename, names[-1])
    else:
      signature = ''
      for i in range(0, len(names)-1):
        signature += self.architecture.typename + ' const*, '
      signature += self.architecture.typename + '*'
    return signature
    
  def __localArray(self, name, reals, aligned=True):
    aligned = ' __attribute__((aligned(PAGESIZE_STACK)))' if aligned else ''
    return '  {} {}[{}]{};\n'.format(self.architecture.typename, name, reals, aligned)
    

  def generateKernels(self, outputDir, kernels):
    with open(outputDir + '/kernels.h', 'w+') as header, open(outputDir + '/kernels.cpp', 'w+') as cpp:
      self.__headerGuardBegin(header, 'KERNELS')
      namespaceBegin = """namespace seissol {
namespace generatedKernels {
"""
      header.write(namespaceBegin)
      cpp.write("""#include <cstring>
#include <Initializer/preProcessorMacros.fpp>
#include "gemms.h"
""")
      cpp.write(namespaceBegin)
      luts = dict()
      signatures = dict()
      for functionName, kernel in kernels:
        self.__result = kernel
        temps = list()
        involvedMatrices = set()
        inner = self.__traverseForGemms(kernel.symbol, involvedMatrices, temps)
        temps.sort(key=lambda temp: temp.name)
        involvedMatrices = sorted(list(involvedMatrices))
        involvedMatrices.append(self.resultName)
          
        brackets = re.match(r'^(\w+)\[(\d+)\]$', functionName)
        if brackets != None:
          functionName = brackets.expand(r'\1\2')
          base = brackets.group(1)
          index = int(brackets.group(2))
          if not luts.has_key(base):
            luts[base] = dict()
            signatures[base] = self.__functionSignature(involvedMatrices, writeNames=False)
          luts[base].update({index: base})

        functionSignature = 'void {}({})'.format(functionName, self.__functionSignature(involvedMatrices))
        code = 'void {}({}) {{\n'.format(functionName, self.__functionSignature(involvedMatrices))
        localVariables = ''
        for temp in temps:
          localVariables += self.__localArray(temp.name, temp.requiredReals)
        header.write(functionSignature + ';\n')
        cpp.write('{} {{\n{}{}}}\n'.format(functionSignature, localVariables, inner))
      for key, value in luts.iteritems():
        maxkey = max(value.keys())
        header.write('static void (*{}[])({}) = {{ '.format(key, signatures[key]))
        for i in range(0, maxkey+1):
          header.write(value[i] + str(i) if value.has_key(i) else 'NULL')
          if i != maxkey:
            header.write(', ')
        header.write(' };\n')
      namespaceEnd = '}\n}\n'
      header.write(namespaceEnd)
      cpp.write(namespaceEnd)
      self.__headerGuardEnd(header)

  def __traverseForGemms(self, expr, involvedMatrices, temps):
    code = ''    
    if expr.is_MatAdd:
      for arg in expr.args:
          code += self.__traverseForGemms(arg, involvedMatrices, temps)
    elif expr.is_MatMul:
      matrices = list()
      for symbol in expr.args:
        matrices.append(self.db[symbol.name])
        involvedMatrices.update([symbol.name])
      code += self.__gemms(matrices, temps)
    else:
      raise ValueError('Unexpected subexpression {}.'.format(expr))
    return code
    
  def __gemms(self, matrices, temps):
    code = ''
    nameToIndex = dict()
    for i,matrix in enumerate(matrices):
      nameToIndex[matrix.name] = i
    
    sparsityPatterns = [matrix.spp for matrix in matrices]
    equivalentSparsityPatterns = equivalentMultiplicationPatterns(sparsityPatterns)
    fittedBlocks = [patternFittedBlocks(matrix.blocks, equivalentSparsityPatterns[i]) for i, matrix in enumerate(matrices)]
    implementationPatterns = [getImplementationPattern(fittedBlocks[i], equivalentSparsityPatterns[i]) for i in range(0, len(matrices))]
    chainOrder = sparseMatrixChainOrder(implementationPatterns)

    stack = list()
    output = list()
    stack.append((0, len(matrices)-1))
    # convert matrix chain order to postfix
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
    tempCounter = len(temps)
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
          if len(temps) == 0:
            tempCounter += 1
            temps.append(MatrixInfo(self.tempBaseName + str(tempCounter)))
          result = temps.pop()
          resultRequiredReals = result.requiredReals
          result = MatrixInfo(result.name, op1.rows, op2.cols, sparsityPattern = op1.spp * op2.spp)
          result.fitBlocksToSparsityPattern()
          result.generateMemoryLayout(self.architecture, alignStartrow=True)
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
            code += self.__memset(resultName, result.requiredReals)
          else:
            beta = 0
          result.requiredReals = max(resultRequiredReals, result.requiredReals)
        else:
          beta = 1
          result = self.__result
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
            code += self.__gemm(op1.name, block1, op1.blocks[i1], op2.name, block2, op2.blocks[i2], resultName, result.blocks[0], beta)
            
        # if op is temporary
        if not nameToIndex.has_key(op1.name):
          temps.append(op1)
        if not nameToIndex.has_key(op2.name):
          temps.append(op2)
    return code
  
  def __intersect(self, blockA, blockB):
    return (max(blockA.startcol, blockB.startrow), min(blockA.stopcol, blockB.stoprow))
  
  def __memset(self, name, reals):
    return '  memset({}, 0, {} * sizeof({}));\n'.format(name, reals, self.architecture.typename)
    
  def __headerGuardBegin(self, header, name):
    header.write('#ifndef GENERATED_{}_H_\n'.format(name))
    header.write('#define GENERATED_{}_H_\n'.format(name))

  def __headerGuardEnd(self, header):
    header.write('#endif\n')

  def __gemm(self, nameA, blockA, memoryBlockA, nameB, blockB, memoryBlockB, nameC, memoryBlockC, beta):
    code = ''
    n1, n2 = self.__intersect(blockA, blockB)
    if n2 > n1:
      m1 = getAlignedIndex(blockA.startrow, self.architecture)
      m2 = m1 + getAlignedDim(blockA.stoprow - m1, self.architecture)
      offsetA = memoryBlockA.offset + (m1 - memoryBlockA.startrow) + memoryBlockA.ld * (n1 - memoryBlockA.startcol)
      offsetB = memoryBlockB.offset + (n1 - memoryBlockB.startrow) + memoryBlockB.ld * (blockB.startcol - memoryBlockB.startcol)
      offsetC = memoryBlockC.offset + (m1 - memoryBlockC.startrow) + memoryBlockC.ld * (blockB.startcol - memoryBlockC.startcol)
      alignedA = checkAlignment(offsetA, self.architecture)
      alignedC = checkAlignment(offsetC, self.architecture)
      
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
      code += '  {}({}, {}, {}, NULL, NULL, NULL);\n'.format(
                        generateRoutineName(gemm),
                        formatOffset(nameA, offsetA),
                        formatOffset(nameB, offsetB),
                        formatOffset(nameC, offsetC)
                      )
      
    return code
    
  def generateInitializer(self, outputDir):
    hFilename = 'init.h'
    with open(outputDir + '/' + hFilename, 'w+') as header, open(outputDir + '/init.cpp', 'w+') as cpp:
      self.__headerGuardBegin(header, 'INIT')
      header.write('#include <cstring>\n')
      cpp.write('#include "{}"\n'.format(path.basename(hFilename)))
      namespaceBegin = """
namespace seissol {
namespace model {
"""
      header.write(namespaceBegin)
      cpp.write(namespaceBegin)
      globalMatrixValues = dict()
      maxGlobalMatrixId = -1
      for matrixInfo in self.db.itervalues():
        cpp.write('void {namespace}::convertToDense({typename}* const matrix, {typename}* denseMatrix) {{\n'.format(
          namespace=matrixInfo.name,
          typename=self.architecture.typename
        ))
        cpp.write(self.__memset('denseMatrix', matrixInfo.rows * matrixInfo.cols))
        for block in matrixInfo.blocks:
          cpp.write("""
  for (unsigned col = {block.startcol}; col < {block.stopcol}; ++col) {{
    for (unsigned row = {block.startrow}; row < {block.stoprow}; ++row) {{
      denseMatrix[col * {matrixInfo.rows} + row] = matrix[{block.offset} + (col - {block.startcol}) * {block.ld} + (row - {block.startrow})];
    }}
  }}
""".format(matrixInfo=matrixInfo, block=block, typename=self.architecture.typename))
        cpp.write('}\n')

        if matrixInfo.values != None:
          init = '  extern {} const values[];'.format(self.architecture.typename)
          cpp.write('{} const {}::values[] = {{ {} }};\n'.format(
            self.architecture.typename,
            matrixInfo.name,
            ', '.join(matrixInfo.getValuesAsStoredInMemory())
          ))
          globalMatrixValues[matrixInfo.globalMatrixId] = matrixInfo.name
          maxGlobalMatrixId = max(maxGlobalMatrixId, matrixInfo.globalMatrixId)
        else:
          init = """  static unsigned index(unsigned row, unsigned column) {{
    static unsigned const lut[] = {{ {lut} }};
    return lut[row + column*{matrixInfo.rows}];
  }}
  static void setZero({typename}* data) {{
    memset(data, 0, {matrixInfo.requiredReals}*sizeof({typename}));
  }}""".format(lut=', '.join(map(str, matrixInfo.getIndexLUT())), matrixInfo=matrixInfo, typename=self.architecture.typename)
  #~ void setRandom({typename}* data) {{
    #~ for (unsigned col = 0; col < {matrixInfo.cols}; ++col) {{
      #~ for (unsigned row = 0; row < {matrixInfo.rows}; ++row) {{
        #~ unsigned idx = index(row, col);
        #~ if (idx != -1) {{
          #~ data[idx] = drand48();
        #~ }}
      #~ }}
    #~ }}
  #~ }}""".format(lut=', '.join(map(str, matrixInfo.getIndexLUT())), matrixInfo=matrixInfo, typename=self.architecture.typename)

        code = """
namespace {matrixInfo.name} {{
  unsigned const rows = {matrixInfo.rows};
  unsigned const cols = {matrixInfo.cols};
  unsigned const reals = {matrixInfo.requiredReals};
{init}
  void convertToDense({typename}* const matrix, {typename}* denseMatrix);
}}\n""".format(typename=self.architecture.typename, init=init, matrixInfo=matrixInfo)
        header.write(code)
        
      # Pointer table for global matrices
      cpp.write('{} const*const globalMatrixValues[] = {{ {} }};\n'.format(
        self.architecture.typename,
        ', '.join(['&' + globalMatrixValues[i] + '::values[0]' if globalMatrixValues.has_key(i) else 'NULL' for i in range(0, maxGlobalMatrixId+1)])))
      cpp.write('unsigned const globalMatrixOffsets[] = { 0')
      globalMatrixOffset = 0
      for i in range(0, maxGlobalMatrixId+1):
        globalMatrixOffset += self.db[globalMatrixValues[i]].requiredReals if globalMatrixValues.has_key(i) else 0
        cpp.write(', {}'.format(globalMatrixOffset))
      cpp.write('};\n')
      header.write('extern {} const*const globalMatrixValues[];\n'.format(self.architecture.typename))
      header.write('extern unsigned const globalMatrixOffsets[];\n'.format(self.architecture.typename))
      header.write('unsigned const numGlobalMatrices = {};\n'.format(maxGlobalMatrixId+1))
        
      namespaceEnd = '}\n}\n'
      header.write(namespaceEnd)
      cpp.write(namespaceEnd)
      self.__headerGuardEnd(header)

  def generateUnitTests(self, functionName, kernel):
    involvedMatrices = set()
    inner = self.__traverseForUnitTests(kernel.symbol, involvedMatrices)
    involvedMatrices = sorted(list(involvedMatrices))
    code = 'void {}_reference({}) {{\n'.format(functionName, self.__functionSignature(involvedMatrices + ['result_dense']))
    declarations = ''
    initializations = ''
    for matrix in involvedMatrices:
      matrixInfo = self.db[matrix]
      declarations += self.__localArray(matrix + '_dense', matrixInfo.rows * matrixInfo.cols, aligned=False)
      initializations += '  {name}::convertToDense({name}, {name}_dense);\n'.format(name=matrix)
    code += declarations
    code += initializations
    code += inner
    code += '}\n'
    code += 'void test_{}() {{\n'.format(functionName)
    declarations = ''
    initializations = ''
    for matrix in involvedMatrices:
      matrixInfo = self.db[matrix]
      declarations += self.__localArray(matrix, matrixInfo.requiredReals)
      if matrixInfo.values != None:
        initializations += '  memcpy({name}, {name}::values, {reals} * sizeof({typename}));\n'.format(name=matrix, reals=matrixInfo.requiredReals, typename=self.architecture.typename)
      else:
        initializations += '  {name}::setRandom({name});\n'.format(name=matrix)        
    declarations += self.__localArray('reference', kernel.rows * kernel.cols)
    initializations += self.__memset('reference', kernel.rows * kernel.cols)
    declarations += self.__localArray('result', kernel.requiredReals)
    initializations += self.__memset('result', kernel.requiredReals)
    code += declarations
    code += initializations
    signature = ', '.join(involvedMatrices)
    code += '  {}_reference({}, {});\n'.format(functionName, signature, 'reference')
    code += '  {}({}, {});\n'.format(functionName, signature, 'result')
    code += """  double error = 0.0;
  double refNorm = 0.0;
  for (unsigned col = {block.startcol}; col < {block.stopcol}; ++col) {{
    for (unsigned row = {block.startrow}; row < {block.stoprow}; ++row) {{
      double ref = reference[col*{rows} + row];
      double diff = ref - result[col*{block.ld} + row];
      error += diff * diff;
      refNorm += ref * ref;
    }}
  }}
  //printf("{name} error (Frobenius, absolute): %e\\n", sqrt(error));
  printf("{name} error (Frobenius, relative): %e\\n", sqrt(error) / sqrt(refNorm));\n""".format(rows=kernel.rows, block=kernel.blocks[0], name=functionName)
    code += '}\n'
    return code

  def __traverseForUnitTests(self, expr, involvedMatrices):
    code = ''    
    if expr.is_MatAdd:
      for arg in expr.args:
          code += self.__traverseForUnitTests(arg, involvedMatrices)
    elif expr.is_MatMul:
      matrices = list()
      for symbol in expr.args:
        matrices.append(self.db[symbol.name])
        involvedMatrices.update([symbol.name])
      resultSize = 0
      leftName = matrices[0].name
      leftRows = matrices[0].rows
      leftCols = matrices[0].cols
      for i in range(1, len(matrices)):
        if i < len(matrices)-1:
          tempName = '{}_temp{}'.format('_'.join([symbol.name for symbol in expr.args]), i)
          code += '  {} {}_dense[{}];\n'.format(self.architecture.typename, tempName, leftRows * matrices[i].cols)
          beta = 0.0
        else:
          tempName = self.resultName
          beta = 1.0
        code += '  cblas_{}gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, {m}, {n}, {k}, 1.0, {A}_dense, {ldA}, {B}_dense, {ldB}, {beta}, {C}_dense, {ldC});\n'.format(
          self.architecture.precision.lower(),
          m=leftRows,
          n=matrices[i].cols,
          k=leftCols,
          A=leftName,
          ldA=leftRows,
          B=matrices[i].name,
          ldB=matrices[i].rows,
          beta=beta,
          C=tempName,
          ldC=leftRows
        )
        leftName = tempName
        leftCols = matrices[i].cols
    else:
      raise ValueError('Unexpected subexpression {}.'.format(expr))
    return code
