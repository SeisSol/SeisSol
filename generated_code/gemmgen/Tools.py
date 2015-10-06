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

from lxml import etree
from DB import MatrixInfo, getAlignedDim
from Expr import analyseKernels
from Generator import Generator
from sympy import pprint, pretty

def __complain(child):
  raise ValueError('Unknown tag ' + child.tag)

def __parseMatrix(node, clones):
  name = node.get('name')
  rows = int( node.get('rows') )
  columns = int( node.get('columns') )
  if len(node) != 0:
    spp = (list(), list())
    values = list()
    for child in node:
      if child.tag == 'entry':
        row = int(child.get('row'))-1
        col = int(child.get('column'))-1
        spp[0].append(row)
        spp[1].append(col)
        if values != None and 'value' in child.keys():
          values.append((row, col, child.get('value')))
        else:
          values = None
      else:
        self.__complain(child)
  else:
    spp = None
    values = None

  dbUpdate = dict()
  if clones.has_key(name):
    for clone in clones[name]:
      dbUpdate[clone] = MatrixInfo(clone, rows, columns, spp, values)
  else:
    dbUpdate[name] = MatrixInfo(name, rows, columns, spp, values)
  
  return dbUpdate  

def parseMatrixFile(xmlFile, clones):
  tree = etree.parse(xmlFile)
  root = tree.getroot()
  
  matrices = dict()
  
  for child in root:
    if child.tag == 'matrix':
      matrices.update( __parseMatrix(child, clones) )
    else:
      __complain(child)
  
  return matrices
  
def numberOfBasisFunctions(order):
  return order * (order + 1) * (order + 2) / 6
  
def alignedNumberOfBasisFunctions(order, architecture):
  return getAlignedDim(numberOfBasisFunctions(order), architecture)

def generate(outputDir, db, kernels, libxsmmGenerator, architecture):
  analyseKernels(db, kernels)

  for matrixInfo in db.itervalues():
    matrixInfo.generateMemoryLayout(architecture, alignStartrow=matrixInfo.leftMultiplication)    
  for name, kernel in kernels:
    kernel.generateMemoryLayout(architecture, alignStartrow=True)

  print('\nKernels')
  print('-------')
  for name, kernel in kernels:
    print(u'{}: {}'.format(name, pretty(kernel.symbol)))

  print('\nMemory layout')
  print('-------------')
  keys = db.keys()
  keys.sort(key=lambda s: s.lower())
  for key in keys:
    name = db[key].name
    for block in db[key].blocks:
      print('{:16} {}'.format(name, block))
      name = ''

  generator = Generator(db, libxsmmGenerator, architecture)
  generator.generateKernels(outputDir, kernels)
  generator.generateInitializer(outputDir)
  generator.generateGemms(outputDir)
  
