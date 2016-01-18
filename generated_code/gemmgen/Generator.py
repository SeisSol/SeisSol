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

import Code
import Kernel
import os
import itertools
import re
import tempfile
import scipy.io
import scipy.sparse
    
def generateRoutineName(gemm):
  name = 'sparse_' + gemm['spMtxName'] if gemm['spp'] is not None else 'gemm'
  lda = 'Asparse' if gemm['LDA'] <= 0 else 'ldA{}'.format(gemm['LDA'])
  ldb = 'Bsparse' if gemm['LDB'] <= 0 else 'ldB{}'.format(gemm['LDB'])
  return '{}_m{}_n{}_k{}_{}_{}_ldC{}_beta{}_alignedA{}_alignedC{}_pfsigonly'.format(
    name,
    gemm['M'],
    gemm['N'],
    gemm['K'],
    lda,
    ldb,
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
    
def functionName(name):
  brackets = re.match(r'^(\w+)\[(\d+)\]$', name)
  if brackets != None:
    functionName = brackets.expand(r'\1\2')
    base = brackets.group(1)
    index = int(brackets.group(2))
    if index < 0:
      raise ValueError('Invalid name {}. Index must be greater than 0.'.format(name))
  else:
    functionName = name
    base = name
    index = -1
  return (functionName, base, index)

class Generator:
  def __init__(self, db, libxsmmGenerator, architecture):
    self.db = db
    self.libxsmmGenerator = libxsmmGenerator
    self.architecture = architecture
    
  def __generateGemms(self, outputDir, gemmlist):
    cppFilename = outputDir + '/gemms.cpp'
    hFilename = outputDir + '/gemms.h'
    
    with Code.Cpp(cppFilename) as cpp:
      cpp('#ifndef NDEBUG')
      cpp('extern long long libxsmm_num_total_flops;')
      cpp('#endif')
      cpp('#if defined( __SSE3__) || defined(__MIC__)')
      cpp('#include <immintrin.h>')
      cpp('#endif')

    with Code.Cpp(hFilename) as header:
      with header.HeaderGuard('GEMMS'):
        indexnamelist = [(i, generateRoutineName(gemm)) for i, gemm in enumerate(gemmlist)]
        keyFunc = lambda x: x[1]
        indexnamelist.sort(key=keyFunc)
        uniqueindexnames = [group.next() for key, group in itertools.groupby(indexnamelist, key=keyFunc)]
        for index, name in uniqueindexnames:
          header('void {name}({tn} const* A, {tn} const* B, {tn}* C, {tn} const* A_prefetch, {tn} const* B_prefetch, {tn} const* C_prefetch);'.format(
            name=name,
            tn=self.architecture.typename
          ))
          spp = gemmlist[index]['spp']
          if spp is not None:
            temp = tempfile.NamedTemporaryFile()
            scipy.io.mmwrite(temp, scipy.sparse.coo_matrix(spp).asformat('csc'))
            sppFile = temp.name
          else:
            sppFile = ''
          os.system(self.__generateLibxsmmGeneratorCall(cppFilename, gemmlist[index], sppFile))

  def __generateLibxsmmGeneratorCall(self, filename, gemm, sppFile):
    return '{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} pfsigonly {}P {}'.format(
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
      self.architecture.precision,
      sppFile
    )
    
  def __gemmSignature(self, names, writeNames=True):
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
    return '{} {}[{}]{};\n'.format(self.architecture.typename, name, reals, aligned)

  def generateKernels(self, outputDir, kernels):
    luts = dict()
    signatures = dict()
    flops = dict()
    generatedKernels = list()
    gemmlist = list()
    for name, kernel in kernels:
      gk = Kernel.GeneratedKernel(kernel, self.db, self.architecture)
      flop = (gk.nonZeroFlops, gk.hardwareFlops)
      funName, base, index = functionName(name)
      if index >= 0:
        if not luts.has_key(base):
          luts[base] = dict()
          signatures[base] = self.__gemmSignature(gk.involvedMatrices, writeNames=False)
        luts[base].update({index: base})        
        if not flops.has_key(base):
          flops[base] = dict()
        flops[base].update({index: flop})
      else:
        flops[funName] = flop
      generatedKernels.append( (funName, gk) )
      gemmlist.extend(gk.gemmlist)
      
    self.__generateGemms(outputDir, gemmlist)

    with Code.Cpp(outputDir + '/kernels.h') as header:
      with header.HeaderGuard('KERNELS'):
        with header.Namespace('seissol'):
          with header.Namespace('generatedKernels'):
            for name, gk in generatedKernels:
              header('void {}({});'.format(name, self.__gemmSignature(gk.involvedMatrices)))
            for key, value in luts.iteritems():
              maxkey = max(value.keys())
              pointers = [value[i] + str(i) if value.has_key(i) else '0' for i in range(0, maxkey+1)]
              header('static void (*{}[])({}) = {{ {} }};'.format(key, signatures[key], ', '.join(pointers)))
            
    with Code.Cpp(outputDir + '/kernels.cpp') as cpp:
      cpp.includeSys('cstring')
      cpp.includeSys('Initializer/preProcessorMacros.fpp')
      cpp.include('gemms.h')
      with cpp.Namespace('seissol'):
        with cpp.Namespace('generatedKernels'):
          for name, gk in generatedKernels:
            with cpp.Function('void {}({})'.format(name, self.__gemmSignature(gk.involvedMatrices))):
              for temp in gk.temps:
                cpp(self.__localArray(temp.name, temp.requiredReals))
              for operation in gk.operations:
                if operation['type'] == Kernel.Operation.MEMSET:
                  cpp.memset(operation['pointer'], operation['numberOfReals'], operation['dataType'], operation['offset'])
                elif operation['type'] == Kernel.Operation.GEMM:
                  cpp('{}({}, {}, {}, NULL, NULL, NULL);'.format(
                    generateRoutineName(operation['gemm']),
                    formatOffset(operation['nameA'], operation['offsetA']),
                    formatOffset(operation['nameB'], operation['offsetB']),
                    formatOffset(operation['nameC'], operation['offsetC'])
                  ))
                  
    with Code.Cpp(outputDir + '/flops.h') as header:
      with header.HeaderGuard('FLOPS'):
        with header.Namespace('seissol'):
          with header.Namespace('flops'):
            for key, value in flops.iteritems():
              if isinstance(value, dict):
                maxkey = max(value.keys())
                nonZeroFlops = [str(value[i][0]) if value.has_key(i) else '0' for i in range(0, maxkey+1)]
                header('unsigned const {}_nonZero[] = {{ {} }};'.format(key, ', '.join(nonZeroFlops)))
                hardwareFlops = [str(value[i][1]) if value.has_key(i) else '0' for i in range(0, maxkey+1)]
                header('unsigned const {}_hardware[] = {{ {} }};'.format(key, ', '.join(hardwareFlops)))
              else:
                header('unsigned const {}_nonZero = {};'.format(key, value[0]))
                header('unsigned const {}_hardware = {};'.format(key, value[1]))
    
  def generateInitializer(self, outputDir):
    globalMatrixValues = dict()
    maxGlobalMatrixId = -1
    for matrixInfo in self.db.itervalues():
      if matrixInfo.values != None:
        globalMatrixValues[matrixInfo.globalMatrixId] = matrixInfo.name
        maxGlobalMatrixId = max(maxGlobalMatrixId, matrixInfo.globalMatrixId)
        
    globalMatrixOffsets = [0]
    for i in range(0, maxGlobalMatrixId+1):
      offset = self.db[globalMatrixValues[i]].requiredReals if globalMatrixValues.has_key(i) else 0
      globalMatrixOffsets.append(globalMatrixOffsets[-1] + offset)
      
    with Code.Cpp(outputDir + '/sizes.h') as header:
      with header.HeaderGuard('SIZES'):
        with header.Namespace('seissol'):
          with header.Namespace('model'):
            for matrixInfo in self.db.itervalues():
              with header.Namespace(matrixInfo.name):
                header('unsigned const rows = {};'.format(matrixInfo.rows))
                header('unsigned const cols = {};'.format(matrixInfo.cols))
                header('unsigned const reals = {};'.format(matrixInfo.requiredReals))
                if len(matrixInfo.blocks) == 1 and matrixInfo.blocks[0].ld > 0:
                  header('unsigned const ld = {};'.format(matrixInfo.blocks[0].ld))
    
    hFilename = 'init.h'
    with Code.Cpp(outputDir + '/' + hFilename) as header:
      with header.HeaderGuard('INIT'):
        header.includeSys('cstring')
        with header.Namespace('seissol'):
          with header.Namespace('model'):
            for matrixInfo in self.db.itervalues():
              with header.Namespace(matrixInfo.name):
                if matrixInfo.values == None:
                  header('void convertToDense({typename} const* matrix, {typename}* denseMatrix);'.format(typename=self.architecture.typename))
                  with header.Function('static int index(unsigned row, unsigned column)'):
                    header('static int const lut[] = {{ {} }};'.format(
                      ', '.join(map(str, matrixInfo.getIndexLUT()))
                    ))
                    header('return lut[row + column*{}];'.format(matrixInfo.rows))
                  with header.Function('static void setZero({}* data)'.format(self.architecture.typename)):
                    header.memset('data', matrixInfo.requiredReals, self.architecture.typename)
                else:
                  header('extern {} const values[];'.format(self.architecture.typename))
                  with header.Function('static void convertToDense({typename}* denseMatrix)'.format(typename=self.architecture.typename)):
                    header('static {} const denseValues[] = {{ {} }};'.format(
                      self.architecture.typename,
                      ', '.join(matrixInfo.getValuesDense())
                    ))
                    header('memcpy(denseMatrix, denseValues, {reals} * sizeof({typename}));'.format(reals=matrixInfo.rows*matrixInfo.cols, typename=self.architecture.typename))
                  
            header('extern {} const*const globalMatrixValues[];'.format(self.architecture.typename))
            header('extern unsigned const globalMatrixOffsets[];')
            header('unsigned const numGlobalMatrices = {};'.format(maxGlobalMatrixId+1))
                  

    with Code.Cpp(outputDir + '/init.cpp') as cpp:
      cpp.include(hFilename)      
      with cpp.Namespace('seissol'):
        with cpp.Namespace('model'):
          for matrixInfo in self.db.itervalues():
            if matrixInfo.values == None:
              with cpp.Function('void {namespace}::convertToDense({typename} const* matrix, {typename}* denseMatrix)'.format(namespace=matrixInfo.name, typename=self.architecture.typename)):
                cpp.memset('denseMatrix', matrixInfo.rows * matrixInfo.cols, self.architecture.typename)
                for block in matrixInfo.blocks:
                  with cpp.For('unsigned col = {block.startcol}; col < {block.stopcol}; ++col'.format(block=block)):
                    with cpp.For('unsigned row = {block.startrow} + {block.startpaddingrows}; row < {block.stoprow}; ++row'.format(block=block)):
                      if block.sparse:
                        cpp('unsigned idx = index(row, col);')
                        with cpp.If('idx != -1'):
                          cpp('denseMatrix[col * {matrixInfo.rows} + row] = matrix[idx];'.format(matrixInfo=matrixInfo))
                      else:
                        cpp('denseMatrix[col * {matrixInfo.rows} + row] = matrix[{block.offset} + (col - {block.startcol}) * {block.ld} + (row - {block.startrow})];'.format(
                          matrixInfo=matrixInfo,
                          block=block
                        ))            
            else:
              cpp('{} const {}::values[] = {{ {} }};'.format(
                self.architecture.typename,
                matrixInfo.name,
                ', '.join(matrixInfo.getValuesAsStoredInMemory())
              ))

          cpp('{} const*const globalMatrixValues[] = {{ {} }};'.format(
            self.architecture.typename,
            ', '.join(['&' + globalMatrixValues[i] + '::values[0]' if globalMatrixValues.has_key(i) else 'NULL' for i in range(0, maxGlobalMatrixId+1)])
          ))            
          cpp('unsigned const globalMatrixOffsets[] = {{ {} }};'.format(
            ', '.join(map(str, globalMatrixOffsets))
          ))
          
  def generateUnitTests(self, outputDir, kernels):
    referenceKernels = list()
    for name, kernel in kernels:
      rk = Kernel.ReferenceKernel(kernel, self.db, self.architecture)
      funName, base, index = functionName(name)  
      referenceKernels.append( (funName, rk) )
    
    with Code.Cpp(outputDir + '/KernelTests.t.h') as test:
      with test.HeaderGuard('TEST'):
        test.includeSys('cstdlib')
        test.includeSys('cstring')
        test.includeSys('ctime')
        test.includeSys('cxxtest/TestSuite.h')
        test.includeSys('Initializer/preProcessorMacros.fpp')
        test.include('init.h')
        test.include('kernels.h')
        with test.Namespace('seissol'):
          with test.Namespace('unit_test'):
            with test.Function('void gemm(unsigned m, unsigned n, unsigned k, {tn}* A, {tn}* B, {tn} beta, {tn}* C)'.format(tn=self.architecture.typename)):
              with test.If('beta == 0.0'):
                test.memset('C', 'm*n', self.architecture.typename)
              with test.For('unsigned col = 0; col < n; ++col'):
                with test.For('unsigned row = 0; row < m; ++row'):
                  with test.For('unsigned s = 0; s < k; ++s'):
                    test('C[row + col * m] += A[row + s * m] * B[s + col * k];')
            for name, rk in referenceKernels:
              with test.Function('void {}({})'.format(name, self.__gemmSignature(rk.involvedMatrices))):
                for matrix in rk.involvedMatrices:
                  if self.db.has_key(matrix):
                    matrixInfo = self.db[matrix]
                    test(self.__localArray(matrix + '_dense', matrixInfo.rows * matrixInfo.cols, aligned=False))
                    if matrixInfo.values == None:
                      test('seissol::model::{name}::convertToDense({name}, {name}_dense);'.format(name=matrix))
                    else:
                      test('seissol::model::{name}::convertToDense({name}_dense);'.format(name=matrix))
                for temp in rk.temps:
                  test(self.__localArray(temp.name, temp.rows * temp.cols, aligned=False))
                for operation in rk.operations:
                  if operation['type'] == Kernel.Operation.GEMM:
                    test('gemm({}, {}, {}, {}, {}, {}, {});'.format(
                      operation['m'],
                      operation['n'],
                      operation['k'],
                      operation['A'] + ('_dense' if self.db.has_key(operation['A']) else ''),
                      operation['B'] + ('_dense' if self.db.has_key(operation['B']) else ''),
                      operation['beta'],                        
                      operation['C']
                    ))
            test('class KernelTestSuite;')

        with test.Class('seissol::unit_test::KernelTestSuite : public CxxTest::TestSuite'):
          test.label('public')
          with test.Function('KernelTestSuite* createSuite()'):
            test('srand48(time(NULL));')
            test('return new KernelTestSuite();')
          for name, rk in referenceKernels:              
            with test.Function('void test{}()'.format(name.capitalize())):
              for matrix in rk.involvedMatrices:
                if self.db.has_key(matrix):
                  matrixInfo = self.db[matrix]
                  test(self.__localArray(matrix, matrixInfo.requiredReals))
                  if matrixInfo.values != None:
                    test('memcpy({name}, seissol::model::{name}::values, {reals} * sizeof({typename}));'.format(name=matrix, reals=matrixInfo.requiredReals, typename=self.architecture.typename))
                  else:
                    test.memset(matrix, matrixInfo.requiredReals, self.architecture.typename)
                    with test.For('unsigned col = 0; col < {}; ++col'.format(matrixInfo.cols)):
                      with test.For('unsigned row = 0; row < {}; ++row'.format(matrixInfo.rows)):
                        test('unsigned idx = seissol::model::{}::index(row, col);'.format(matrix))
                        with test.If('idx != -1'):
                          test('{name}[idx] = drand48();'.format(name=matrix))
              referenceReals = rk.kernel.rows * rk.kernel.cols
              test(self.__localArray('reference', referenceReals))
              test.memset('reference', referenceReals, self.architecture.typename)
              test(self.__localArray('result', rk.kernel.requiredReals))
              test.memset('result', rk.kernel.requiredReals, self.architecture.typename)
              
              test('seissol::unit_test::{}({});'.format(name, ', '.join(rk.involvedMatrices[0:-1] + ['reference'])))
              test('seissol::generatedKernels::{}({});'.format(name, ', '.join(rk.involvedMatrices)))

              test('double error = 0.0;')
              test('double refNorm = 0.0;')
              resultBlock = rk.kernel.blocks[0]
              with test.For('unsigned col = {block.startcol}; col < {block.stopcol}; ++col'.format(block=resultBlock)):
                with test.For('unsigned row = {block.startrow}; row < {block.stoprow}; ++row'.format(block=resultBlock)):
                  test('double ref = reference[col*{} + row];'.format(rk.kernel.rows))
                  test('double diff = ref - result[col*{} + row];'.format(resultBlock.ld))
                  test('error += diff * diff;')
                  test('refNorm += ref * ref;')
              test('TS_ASSERT_LESS_THAN(sqrt(error/refNorm), 1e-15);')
