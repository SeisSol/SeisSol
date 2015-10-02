#!/usr/bin/env python
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

from gemmgen import DB, Tools, Arch
from numpy import matlib, bmat, float64
import argparse

cmdLineParser = argparse.ArgumentParser()
cmdLineParser.add_argument('--matricesDir')
cmdLineParser.add_argument('--outputDir')
cmdLineParser.add_argument('--arch')
cmdLineParser.add_argument('--order')
cmdLineParser.add_argument('--numberOfMechanisms')
cmdLineParser.add_argument('--generator')
cmdLineArgs = cmdLineParser.parse_args()

architecture = Arch.getArchitectureByIdentifier(cmdLineArgs.arch)
libxsmmGenerator = cmdLineArgs.generator
order = int(cmdLineArgs.order)
numberOfMechanisms = int(cmdLineArgs.numberOfMechanisms)
numberOfBasisFunctions = Tools.numberOfBasisFunctions(order)
numberOfQuantities = 9 + 6*numberOfMechanisms

clones = {
  'star': [ 'AstarT', 'BstarT', 'CstarT' ]
}

db = Tools.parseMatrixFile('{}/matrices_{}.xml'.format(cmdLineArgs.matricesDir, numberOfBasisFunctions), clones)
db.update(Tools.parseMatrixFile('{}/matrices_viscoelastic.xml'.format(cmdLineArgs.matricesDir), clones))
db['timeIntegrated'] = DB.MatrixInfo('timeIntegrated', numberOfBasisFunctions, numberOfQuantities)
db['timeDerivative0'] = DB.MatrixInfo('timeDerivative0', numberOfBasisFunctions, numberOfQuantities)
riemannSolverSpp = bmat([[matlib.ones((9, numberOfQuantities), dtype=float64)], [matlib.zeros((numberOfQuantities-9, numberOfQuantities), dtype=float64)]])
db['AplusT'] = DB.MatrixInfo('AplusT', numberOfQuantities, numberOfQuantities, sparsityPattern=riemannSolverSpp)
#~ db['Aplus'].fitBlocksToSparsityPattern()
db['AminusT'] = DB.MatrixInfo('AminusT', numberOfQuantities, numberOfQuantities, sparsityPattern=riemannSolverSpp)
#~ db['Aminus'].fitBlocksToSparsityPattern()

mechMatrix = matlib.zeros((15, numberOfQuantities))
mechMatrix[0:9,0:9] = matlib.identity(9)
for m in range(0, numberOfMechanisms):
  mechMatrix[9:15,9+6*m:9+6*(m+1)] = matlib.identity(6)
tallMatrix = matlib.zeros((numberOfQuantities, 15))
tallMatrix[0:15,0:15] = db[clones['star'][0]].spp
starMatrix = tallMatrix * mechMatrix
for clone in clones['star']:
  db[clone] = DB.MatrixInfo(clone, starMatrix.shape[0], starMatrix.shape[1], starMatrix)
  #~ db[clone].setBlocks([DB.MatrixBlock(0, 9, 0, 9), DB.MatrixBlock(0, 9, 9, numberOfQuantities)])
  #~ db[clone].fitBlocksToSparsityPattern()
  
source = matlib.zeros((numberOfQuantities, numberOfQuantities))
for m in range(0, numberOfMechanisms):
  r = slice(9+6*m, 9+6*(m+1))
  source[r,0:6] = db['YT'].spp
  source[r,r] = matlib.identity(6)
db['source'] = DB.MatrixInfo('source', numberOfQuantities, numberOfQuantities, source)
#~ db['source'].setBlocks([DB.MatrixBlock(9, numberOfQuantities, 0, 9), DB.MatrixBlock(9, numberOfQuantities, 9, numberOfQuantities)])
#~ db['source'].fitBlocksToSparsityPattern()
  

stiffnessMatrices = ['kXiDivM', 'kEtaDivM', 'kZetaDivM']
transposedStiffnessBlocks = list()
for o in range(2, order+1):
  mb = Tools.alignedNumberOfBasisFunctions(o-1, architecture)
  col = Tools.numberOfBasisFunctions(o-1)
  nb = Tools.numberOfBasisFunctions(o)
  transposedStiffnessBlocks.append(DB.MatrixBlock(0, mb, col, nb))
for matrixName in stiffnessMatrices:
  db[matrixName].fitBlocksToSparsityPattern()
  db[matrixName + 'T'].setBlocks(transposedStiffnessBlocks)
  db[matrixName + 'T'].fitBlocksToSparsityPattern()

stiffnessOrder = { 'Xi': 0, 'Eta': 1, 'Zeta': 2 }
globalMatrixIdRules = [
  (r'^fM(\d{1})$', lambda x: int(x[0])-1),
  (r'^fP(\d{1})(\d{1})(\d{1})$', lambda x: 4 + (int(x[0])-1)*12 + (int(x[1])-1)*3 + (int(x[2])-1)),
  (r'^k(Xi|Eta|Zeta)DivMT$', lambda x: 52 + stiffnessOrder[x[0]]),
  (r'^k(Xi|Eta|Zeta)DivM$', lambda x: 55 + stiffnessOrder[x[0]])
]
DB.determineGlobalMatrixIds(globalMatrixIdRules, db)

kernels = list()

volume = db['timeIntegrated'] * db['source'] \
       + db['kXiDivM'] * db['timeIntegrated'] * db['AstarT'] \
       + db['kEtaDivM'] * db['timeIntegrated'] * db['BstarT'] \
       + db['kZetaDivM'] * db['timeIntegrated'] * db['CstarT']
kernels.append(('volume', volume))

for i in range(0, 4):
  localFlux = db['fM{}'.format(i+1)] * db['timeIntegrated'] * db['AplusT']
  kernels.append(('localFlux[{}]'.format(i), localFlux))

for i in range(0, 4):
  for j in range(0, 4):
    for h in range(0, 3):
      neighboringFlux = db['fP{}{}{}'.format(i+1, j+1, h+1)] * db['timeIntegrated'] * db['AminusT']
      kernels.append(('neighboringFlux[{}]'.format(i*12+j*3+h), neighboringFlux))

for i in range(1, order):
  lastD = 'timeDerivative{}'.format(str(i-1))
  newD  = 'timeDerivative{}'.format(str(i))
  derivative = db['kXiDivMT'] * db[lastD] * db['AstarT'] \
             + db['kEtaDivMT'] * db[lastD] * db['BstarT'] \
             + db['kZetaDivMT'] * db[lastD] * db['CstarT'] \
             + db[lastD] * db['source']
  derivative.fitBlocksToSparsityPattern()
  kernels.append(('derivative[{}]'.format(i), derivative))
  db[newD] = derivative.flat(newD)
  db[newD].fitBlocksToSparsityPattern()

Tools.generate(cmdLineArgs.outputDir, db, kernels, libxsmmGenerator, architecture)

