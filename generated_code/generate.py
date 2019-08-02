#!/usr/bin/env python3
##
# @file
# This file is part of SeisSol.
#
# @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
#
# @section LICENSE
# Copyright (c) 2019, SeisSol Group
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

import argparse
import importlib.util
import inspect

from yateto import useArchitectureIdentifiedBy, Generator
from yateto.gemm_configuration import GeneratorCollection, LIBXSMM, PSpaMM 

import DynamicRupture
import Plasticity
import SurfaceDisplacement
import Point
import memlayout

cmdLineParser = argparse.ArgumentParser()
cmdLineParser.add_argument('--equations')
cmdLineParser.add_argument('--matricesDir')
cmdLineParser.add_argument('--outputDir')
cmdLineParser.add_argument('--arch')
cmdLineParser.add_argument('--order', type=int)
cmdLineParser.add_argument('--numberOfMechanisms', type=int)
cmdLineParser.add_argument('--memLayout')
cmdLineParser.add_argument('--multipleSimulations', type=int)
cmdLineParser.add_argument('--dynamicRuptureMethod')
cmdLineParser.add_argument('--PlasticityMethod')
cmdLineArgs = cmdLineParser.parse_args()

if cmdLineArgs.memLayout == 'auto':
  # TODO(Lukas) Don't hardcode this
  env = {
    'equations': cmdLineArgs.equations,
    'order': cmdLineArgs.order,
    'arch': cmdLineArgs.arch,
    'multipleSimulations': cmdLineArgs.multipleSimulations
  }
  mem_layout = memlayout.guessMemoryLayout(env)
else:
  mem_layout = cmdLineArgs.memLayout
  

arch = useArchitectureIdentifiedBy(cmdLineArgs.arch)

equationsSpec = importlib.util.find_spec(cmdLineArgs.equations)
try:
  equations = equationsSpec.loader.load_module()
except:
  raise RuntimeError('Could not find kernels for ' + cmdLineArgs.equations)

adgArgs = inspect.getargspec(equations.ADERDG.__init__).args[1:]
cmdArgsDict = vars(cmdLineArgs)
cmdArgsDict['memLayout'] = mem_layout
args = [cmdArgsDict[key] for key in adgArgs]
adg = equations.ADERDG(*args)

g = Generator(arch)

# Equation-specific kernels
adg.addInit(g)
adg.addLocal(g)
adg.addNeighbor(g)
adg.addTime(g)

# Common kernels
DynamicRupture.addKernels(g, adg, cmdLineArgs.matricesDir, cmdLineArgs.dynamicRuptureMethod) 
Plasticity.addKernels(g, adg, cmdLineArgs.matricesDir, cmdLineArgs.PlasticityMethod)
SurfaceDisplacement.addKernels(g, adg)
Point.addKernels(g, adg)

# Generate code
gemmTools = GeneratorCollection([LIBXSMM(arch), PSpaMM(arch)])
g.generate(cmdLineArgs.outputDir, 'seissol', gemmTools)
