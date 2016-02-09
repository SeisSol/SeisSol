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

OutputDir = 'layouts'

def writeConfig(fileName, config):
  root = etree.Element('memory_layouts')
  
  for name, layout in config.iteritems():
    if type(layout) == list:
      matrix = etree.SubElement(root, 'matrix', name=name)
      for l in layout:
        sparse = l[4] if len(l) >= 5 else False
        etree.SubElement(matrix, 'block', startrow=str(l[0]), stoprow=str(l[1]), startcol=str(l[2]), stopcol=str(l[3]), sparse=str(sparse))
    elif type(layout) == bool:
      matrix = etree.SubElement(root, 'matrix', name=name, sparse=str(layout))
    else:
      raise ValueError('Unsupported type.')
  
  xml = etree.ElementTree(root)
  with open(fileName, 'w') as f:
    xml.write(f, pretty_print=True)

def generateTuningLayoutFiles(configs):
  prefix = OutputDir + '/'
  writeConfig(prefix + 'dense0.xml', {})

  for name, layouts in configs.iteritems():
    for idx, layout in enumerate(layouts):
      writeConfig('{}{}{}.xml'.format(prefix, name, idx), { name: layout })

def generateLayoutFile(matrices, configs):
  selection = dict()
  for name, idx in matrices:
    selection[name] = configs[name][idx]
  writeConfig('tuned_layout.xml', selection)


def getGlobalMatrices(order, arch):
  architecture = Arch.getArchitectureByIdentifier(arch)

  configs = {
    'kXiDivM': [ True ],
    'kEtaDivM': [ True ],
    'kZetaDivM': [ True ],
    'kXiDivMT': [ True ],
    'kEtaDivMT': [ True ],
    'kZetaDivMT': [ True ],
    'fM1': [ True ],
    'fM2': [ True ],
    'fP111': [ True ],
    'fP112': [ True ],
    'fP113': [ True ],
    'fP121': [ True ],
    'fP211': [ True ],
    'fP222': [ True ]
  }

  stiffnessMatrices = ['kXiDivM', 'kEtaDivM', 'kZetaDivM']
  transposedStiffnessBlocks = list()
  for o in range(2, order+1):
    stoprow = Tools.alignedNumberOfBasisFunctions(o-1, architecture)
    startcol = Tools.numberOfBasisFunctions(o-1)
    stopcol = Tools.numberOfBasisFunctions(o)
    transposedStiffnessBlocks.append((0, stoprow, startcol, stopcol))
  if len(transposedStiffnessBlocks) >= 2:
    # merge first two blocks
    block1 = transposedStiffnessBlocks.pop(0)
    block2 = transposedStiffnessBlocks.pop(0)
    mergedBlock = ( min(block1[0], block2[0]),
                    max(block1[1], block2[1]),
                    min(block1[2], block2[2]),
                    max(block1[3], block2[3]) )
    transposedStiffnessBlocks.insert(0, mergedBlock)

  stiffnessBlocks = [(block[2], block[3], block[0], block[1]) for block in transposedStiffnessBlocks]
  noMemsetStiffnessBlocks = list()
  for i, block in enumerate(stiffnessBlocks):
    startrow = noMemsetStiffnessBlocks[i-1][1] if i > 0 else block[0]
    stoprow = architecture.getAlignedIndex(block[1])
    noMemsetStiffnessBlocks.append( (startrow, stoprow, block[2], block[3]) )

  for matrix in stiffnessMatrices:
    configs[matrix].append(stiffnessBlocks)
    configs[matrix].append(noMemsetStiffnessBlocks)
    configs[matrix + 'T'].append(transposedStiffnessBlocks)
    
  return configs
  
def getStarMatrices(Q):
  return { 'star': [ True, [(0, 6, 6, 9), (6, 9, 0, Q)], [(0, 6, 6, 9), (6, 9, 0, 6), (6, 9, 9, Q)] ] }
    
def getViscoelasticMemoryLayouts(order, numberOfMechanisms, arch):  
  Q = 9 + 6 * numberOfMechanisms
  
  configs = getGlobalMatrices(order, arch)
  configs.update( getStarMatrices(Q) )
  configs.update( { 'source': [ True, [(9, Q, 0, 9), (9, Q, 9, Q, True)] ] } )
  
  return configs
  
def getViscoelastic2MemoryLayouts(order, arch):  
  configs = getGlobalMatrices(order, arch)
  configs.update( getStarMatrices(15) )
  configs.update( { 'ET': [ True ] } )
  
  return configs
  
  
