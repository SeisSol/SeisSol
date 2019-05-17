#!/usr/bin/env python
##
# @file
# This file is part of SeisSol.
#
# @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
#
# @section LICENSE
# Copyright (c) 2015-2016, SeisSol Group
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

from gemmgen import Arch, Tools
from lxml import etree
import numpy
import argparse
import copy

OutputDir = 'layouts'

def writeConfig(fileName, groups, configs):
  root = etree.Element('memory_layouts')
  
  for name, matrices in groups.iteritems():
    group = etree.SubElement(root, 'group', {'name': name})
    for matrix in matrices:
      etree.SubElement(group, 'matrix', {'name': matrix})
  
  for name, layout in configs.iteritems():
    if groups.has_key(name):
      attribs = {'group': name}
    else:
      attribs = {'name': name}
    if type(layout) == list:
      matrix = etree.SubElement(root, 'matrix', attribs)
      for l in layout:
        sparse = l[4] if len(l) >= 5 else False
        etree.SubElement(matrix, 'block', startrow=str(l[0]), stoprow=str(l[1]), startcol=str(l[2]), stopcol=str(l[3]), sparse=str(sparse))
    elif type(layout) == bool:
      attribs.update({'sparse': str(layout)})
      matrix = etree.SubElement(root, 'matrix', attribs)
    else:
      raise ValueError('Unsupported type.')
  
  xml = etree.ElementTree(root)
  with open(fileName, 'w') as f:
    xml.write(f, pretty_print=True)

def generateTuningLayoutFiles(memoryLayouts):
  prefix = OutputDir + '/'
  writeConfig(prefix + 'dense0.xml', memoryLayouts[0], {})

  for name, layouts in memoryLayouts[1].iteritems():
    for idx, layout in enumerate(layouts):
      writeConfig('{}{}{}.xml'.format(prefix, name, idx), memoryLayouts[0], { name: layout })

def generateLayoutFile(configNames, memoryLayouts):
  selection = dict()
  for name, idx in configNames:
    selection[name] = memoryLayouts[1][name][idx]
  writeConfig('tuned_layout.xml', memoryLayouts[0], selection)
  
def mergeBlock(block1, block2):
  return  ( min(block1[0], block2[0]),
            max(block1[1], block2[1]),
            min(block1[2], block2[2]),
            max(block1[3], block2[3]) )

def getGlobalMatrices(order, arch):
  architecture = Arch.getArchitectureByIdentifier(arch)
  stiffnessMatrices = ['kXiDivM', 'kEtaDivM', 'kZetaDivM']
  groups = {
    'stiffnessTransposed': map(lambda x: x + 'T', stiffnessMatrices)
  }

  if architecture.name in ['knc', 'knl']:
    blockMerge = 2
    configs = {
      'kXiDivM': [],
      'kEtaDivM': [],
      'kZetaDivM': [],
      'stiffnessTransposed': [],
      'fP1': [ ],
      'fP2': [ ],
      'fP3': [ ],
      'r1DivM': [ ],
      'r2DivM': [ ],
      'r3DivM': [ ],
      'r4DivM': [ ],
      'rT1': [ ],
      'rT2': [ ],
      'rT3': [ ],
      'rT4': [ ]
    }
  else:
    blockMerge = 1
    configs = {
      'kXiDivM': [ True ],
      'kEtaDivM': [ True ],
      'kZetaDivM': [ True ],
      'stiffnessTransposed': [ True ],
      'fP1': [ True ],
      'fP2': [ True ],
      'fP3': [ True ],
      'r1DivM': [ True ],
      'r2DivM': [ True ],
      'r3DivM': [ True ],
      'r4DivM': [ True ],
      'rT1': [ True ],
      'rT2': [ True ],
      'rT3': [ True ],
      'rT4': [ True ]
    }

  """
  transposedStiffnessBlocks = list()
  for o in range(2, order+1):
    stoprow = Tools.numberOfBasisFunctions(o-1)
    startcol = Tools.numberOfBasisFunctions(o-1)
    stopcol = Tools.numberOfBasisFunctions(o)
    transposedStiffnessBlocks.append((0, stoprow, startcol, stopcol))
  for i in range(blockMerge):
    if len(transposedStiffnessBlocks) >= 2:
      # merge first blocks
      block1 = transposedStiffnessBlocks.pop(0)
      block2 = transposedStiffnessBlocks.pop(0)
      transposedStiffnessBlocks.insert(0, mergeBlock(block1, block2))

  stiffnessBlocks = [(block[2], block[3], block[0], block[1]) for block in transposedStiffnessBlocks]
  noMemsetStiffnessBlocks = list()
  for i, block in enumerate(stiffnessBlocks):
    startrow = noMemsetStiffnessBlocks[i-1][1] if i > 0 else block[0]
    stoprow = architecture.getAlignedIndex(block[1]) if i != len(stiffnessBlocks)-1 else block[1]
    noMemsetStiffnessBlocks.append( (startrow, stoprow, block[2], block[3]) )

  for matrix in stiffnessMatrices:
    configs[matrix].append(stiffnessBlocks)
    configs[matrix].append(noMemsetStiffnessBlocks)

  if groups.has_key('stiffnessTransposed'):
    configs['stiffnessTransposed'].append(transposedStiffnessBlocks)
  else:
    for matrix in stiffnessMatrices:
      configs[matrix + 'T'].append(transposedStiffnessBlocks)
    
  # fP matrices
  fPBlocks = list()
  for o in range(1, order+1):
    start = Tools.numberOfBasisFunctions2D(o-1)
    stop = Tools.numberOfBasisFunctions2D(o)
    fPBlocks.append((start, stop, start, stop))
  # merge first three blocks
  for i in range(blockMerge+1):
    if len(fPBlocks) >= 2:
      block1 = fPBlocks.pop(0)
      block2 = fPBlocks.pop(0)
      fPBlocks.insert(0, mergeBlock(block1, block2))
  for i in range(1,4):
    configs['fP{}'.format(i)].append(fPBlocks)
  
  # rT matrices
  rTBlocks = list()
  for o in range(1, order+1):
    stoprow = Tools.numberOfBasisFunctions2D(o)
    startcol = Tools.numberOfBasisFunctions(o-1)
    stopcol = Tools.numberOfBasisFunctions(o)
    rTBlocks.append((0, stoprow, startcol, stopcol))    
  # merge first three blocks
  for i in range(blockMerge+1):
    if len(rTBlocks) >= 2:
      block1 = rTBlocks.pop(0)
      block2 = rTBlocks.pop(0)
      rTBlocks.insert(0, mergeBlock(block1, block2))
  rBlocks = [(block[2], block[3], block[0], block[1]) for block in rTBlocks]
  noMemsetRBlocks = list()
  for i, block in enumerate(rBlocks):
    startrow = noMemsetRBlocks[i-1][1] if i > 0 else block[0]
    stoprow = architecture.getAlignedIndex(block[1]) if i != len(rBlocks)-1 else block[1]
    noMemsetRBlocks.append( (startrow, stoprow, block[2], block[3]) )
  for i in range(1,5):
    configs['r{}DivM'.format(i)].append(rBlocks)
    configs['r{}DivM'.format(i)].append(noMemsetRBlocks)
    configs['rT{}'.format(i)].append(rTBlocks)
  """

  # fMrT and rT have the same sparsity pattern
  for i in range(1,5):
    configs['fMrT{}'.format(i)] = copy.deepcopy(configs['rT{}'.format(i)])
    
  return (groups, configs)
  
def getStarMatrices(Q):
  star1 = [(0, 6, 6, 9), (6, 9, 0, Q)]
  star2 = [(0, 6, 6, 9), (6, 9, 0, 6)]
  if Q > 9:
    star2.append((6, 9, 9, Q))
  return { 'star': [ True, star1, star2 ] }
  
def getElasticMemoryLayouts(order, arch):
  groups, configs = getGlobalMatrices(order, arch)
  configs.update( getStarMatrices(9) )
  
  return (groups, configs)
    
def getViscoelasticMemoryLayouts(order, numberOfMechanisms, arch):  
  Q = 9 + 6 * numberOfMechanisms
  
  groups, configs = getGlobalMatrices(order, arch)
  configs.update( getStarMatrices(Q) )
  configs.update( { 'source': [ True, [(9, Q, 0, 9), (9, Q, 9, Q, True)] ] } )
  
  return (groups, configs)
  
def getViscoelastic2MemoryLayouts(order, arch):  
  groups, configs = getGlobalMatrices(order, arch)
  configs.update( getStarMatrices(15) )
  configs.update( { 'ET': [ True ] } )
  
  return (groups, configs)
  
  
