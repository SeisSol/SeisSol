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

import lxml.etree
import DB
import Expr
import Generator
import numpy
import scipy.sparse

def __complain(child):
  raise ValueError('Unknown tag ' + child.tag)

def __parseMatrix(node, clones):
  name = node.get('name')
  rows = int( node.get('rows') )
  columns = int( node.get('columns') )
  if len(node) != 0:
    indices = (list(), list())
    data = list()
    for child in node:
      if child.tag == 'entry':
        row = int(child.get('row'))-1
        col = int(child.get('column'))-1
        indices[0].append(row)
        indices[1].append(col)
        if data != None and 'value' in child.keys():
          data.append(child.get('value'))
        else:
          data = None
      else:
        self.__complain(child)
  else:
    indices = None
    data = None

  spp = None
  isConstantGlobalMatrix = False
  if indices is not None:
    if data is None:
      dtype = numpy.float64
      data = numpy.ones(len(indices[0]))
    else:
      dtype = numpy.str_
      isConstantGlobalMatrix = True
    spp = scipy.sparse.coo_matrix((data, indices), shape=(rows, columns), dtype=dtype)

  dbUpdate = DB.DB()
  if clones.has_key(name):
    for clone in clones[name]:
      dbUpdate[clone] = DB.MatrixInfo(clone, rows, columns, spp, isConstantGlobalMatrix)
  else:
    dbUpdate[name] = DB.MatrixInfo(name, rows, columns, spp, isConstantGlobalMatrix)
  
  return dbUpdate

def parseMatrixFile(xmlFile, clones):
  tree = lxml.etree.parse(xmlFile)
  root = tree.getroot()
  
  matrices = DB.DB()
  
  for child in root:
    if child.tag == 'matrix':
      matrices.update( __parseMatrix(child, clones) )
    else:
      __complain(child)
  
  return matrices
  
def memoryLayoutFromFile(xmlFile, db, clones):
  tree = lxml.etree.parse(xmlFile)
  root = tree.getroot()
  strtobool = ['yes', 'true', '1']
  nofits = dict()
  groups = dict()
  
  for group in root.findall('group'):
    groupName = group.get('name')
    noMutualSparsityPattern = group.get('noMutualSparsityPattern', '').lower() in strtobool
    groups[groupName] = list()
    for matrix in group:
      if matrix.tag == 'matrix':
        matrixName = matrix.get('name')
        if not db.has_key(matrixName):
          raise ValueError('Unrecognized matrix name ' + matrixName)
        if len(groups[groupName]) > 0:
          lastMatrixInGroup = groups[groupName][-1]
          if db[lastMatrixInGroup].rows != db[matrixName].rows or db[lastMatrixInGroup].cols != db[matrixName].cols:
            raise ValueError('Matrix {} cannot be in the same group as matrix {} due to different shapes.'.format(matrixName, lastMatrixInGroup))
        groups[groupName].append( matrixName )
      else:
        __complain(group)
    # equalize sparsity pattern
    if not noMutualSparsityPattern:
      spp = None
      for matrix in groups[groupName]:
        spp = spp + db[matrix].spp if spp is not None else db[matrix].spp
      spp[numpy.abs(spp) > 0] = 1.0
      for matrix in groups[groupName]:
        db[matrix].updateSparsityPattern(spp)

  for matrix in root.findall('matrix'):
    group = matrix.get('group')
    name = matrix.get('name')
    nofit = matrix.get('nofit', '').lower() in strtobool
    sparse = matrix.get('sparse', '').lower() in strtobool
    
    if groups.has_key(group) or clones.has_key(name) or db.has_key(name):
      blocks = []
      for block in matrix:
        if block.tag == 'block':
          startrow = int(block.get('startrow'))
          stoprow = int(block.get('stoprow'))
          startcol = int(block.get('startcol'))
          stopcol = int(block.get('stopcol'))
          blksparse = (block.get('sparse') == None and sparse) or block.get('sparse', '').lower() in strtobool
          blocks.append(DB.MatrixBlock(startrow, stoprow, startcol, stopcol, blksparse))
        else:
          __complain(block)
      names = groups[group] if groups.has_key(group) else (clones[name] if clones.has_key(name) else [name])
      for n in names:
        nofits[n] = nofit
        if len(blocks) == 0:
          db[n].setSingleBlock(sparse)
        else:
          db[n].setBlocks(blocks)
        if not nofit:
          db[n].fitBlocksToSparsityPattern()
    else:
      raise ValueError('Unrecognized matrix name ' + name)
  for name, matrix in db.iteritems():
    if not nofits.has_key(name) or nofits[name] == False:
      matrix.fitBlocksToSparsityPattern()

def numberOfBasisFunctions2D(order):
  return order * (order + 1) / 2

def alignedNumberOfBasisFunctions2D(order, architecture):
  return architecture.getAlignedDim(numberOfBasisFunctions2D(order))
  
def numberOfBasisFunctions(order):
  return order * (order + 1) * (order + 2) / 6
  
def alignedNumberOfBasisFunctions(order, architecture):
  return architecture.getAlignedDim(numberOfBasisFunctions(order))

def generate(outputDir, db, kernels, libxsmmGenerator, architecture, prefix=''):
  Expr.analyseKernels(db, kernels)

  for matrixInfo in db.itervalues():
    matrixInfo.generateMemoryLayout(architecture, alignStartrow=matrixInfo.leftMultiplication)    
  for prototype in kernels:
    prototype.kernel.generateMemoryLayout(architecture, alignStartrow=True)

  print('\nKernels')
  print('-------')
  for prototype in kernels:
    print(u'{}: {}'.format(prototype.name, prototype.kernel.symbol))

  print('\nMemory layout')
  print('-------------')
  keys = db.keys()
  keys.sort(key=lambda s: s.lower())
  for key in keys:
    name = db[key].name
    for block in db[key].blocks:
      print('{:16} {}'.format(name, block))
      name = ''

  generator = Generator.Generator(db, libxsmmGenerator, architecture, prefix)
  generator.generateKernels(outputDir, kernels)
  generator.generateInitializer(outputDir)
  generator.generateUnitTests(outputDir, kernels)
  
