#!/usr/bin/env python3
##
# @file
# This file is part of SeisSol.
#
# @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
# @author Sebastian Wolf (wolf.sebastian AT tum.de, https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
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
import sys

from yateto import useArchitectureIdentifiedBy, Generator, NamespacedGenerator
from yateto import gemm_configuration
from yateto.gemm_configuration import GeneratorCollection, LIBXSMM, PSpaMM, MKL, BLIS, OpenBLAS, GemmForge

import DynamicRupture
import Plasticity
import SurfaceDisplacement
import Point
import NodalBoundaryConditions
import memlayout

cmdLineParser = argparse.ArgumentParser()
cmdLineParser.add_argument('--equations')
cmdLineParser.add_argument('--matricesDir')
cmdLineParser.add_argument('--outputDir')
cmdLineParser.add_argument('--host_arch')
cmdLineParser.add_argument('--device_arch', default=None)
cmdLineParser.add_argument('--device_sub_arch', default=None)
cmdLineParser.add_argument('--order', type=int)
cmdLineParser.add_argument('--numberOfMechanisms', type=int)
cmdLineParser.add_argument('--memLayout')
cmdLineParser.add_argument('--multipleSimulations', type=int)
cmdLineParser.add_argument('--dynamicRuptureMethod')
cmdLineParser.add_argument('--PlasticityMethod')
cmdLineParser.add_argument('--gemm_tools')
cmdLineArgs = cmdLineParser.parse_args()

# derive the compute platform
gpu_platforms = ['nvidia', 'amd-gpu']
targets = ['gpu', 'cpu'] if cmdLineArgs.device_arch[1:] in gpu_platforms else ['cpu']

if cmdLineArgs.memLayout == 'auto':
  # TODO(Lukas) Don't hardcode this
  env = {
    'equations': cmdLineArgs.equations,
    'order': cmdLineArgs.order,
    'arch': cmdLineArgs.host_arch,
    'device': cmdLineArgs.device_arch,
    'multipleSimulations': cmdLineArgs.multipleSimulations,
    'targets': targets
  }
  mem_layout = memlayout.guessMemoryLayout(env)
else:
  mem_layout = cmdLineArgs.memLayout


if cmdLineArgs.device_arch == 'none':
    arch = useArchitectureIdentifiedBy(cmdLineArgs.host_arch)
else:
    arch = useArchitectureIdentifiedBy(cmdLineArgs.device_arch, cmdLineArgs.device_sub_arch, cmdLineArgs.host_arch)


equationsSpec = importlib.util.find_spec(cmdLineArgs.equations)
try:
  equations = equationsSpec.loader.load_module()
except:
  raise RuntimeError('Could not find kernels for ' + cmdLineArgs.equations)

cmdArgsDict = vars(cmdLineArgs)
cmdArgsDict['memLayout'] = mem_layout

if cmdLineArgs.equations == 'anisotropic':
    adg = equations.AnisotropicADERDG(**cmdArgsDict)
elif cmdLineArgs.equations == 'elastic':
    adg = equations.ElasticADERDG(**cmdArgsDict)
elif cmdLineArgs.equations == 'viscoelastic':
    adg = equations.ViscoelasticADERDG(**cmdArgsDict)
else:
    adg = equations.Viscoelastic2ADERDG(**cmdArgsDict)

include_tensors = set()
generator = Generator(arch)

# Equation-specific kernels
adg.addInit(generator)
adg.addLocal(generator, targets)
adg.addNeighbor(generator, targets)
adg.addTime(generator, targets)
adg.add_include_tensors(include_tensors)

# Common kernels
include_tensors.update(DynamicRupture.addKernels(NamespacedGenerator(generator, namespace="dynamicRupture"),
                                                 adg,
                                                 cmdLineArgs.matricesDir,
                                                 cmdLineArgs.dynamicRuptureMethod,
                                                 targets))

Plasticity.addKernels(generator, 
                      adg,
                      cmdLineArgs.matricesDir,
                      cmdLineArgs.PlasticityMethod,
                      targets)
NodalBoundaryConditions.addKernels(generator, adg, include_tensors, cmdLineArgs.matricesDir, cmdLineArgs)
SurfaceDisplacement.addKernels(generator, adg)
Point.addKernels(generator, adg)

# pick up the user's defined gemm tools
gemm_tool_list = cmdLineArgs.gemm_tools.replace(" ", "").split(",")
gemm_generators = []

for tool in gemm_tool_list:
  if hasattr(gemm_configuration, tool):
    specific_gemm_class = getattr(gemm_configuration, tool)
    gemm_generators.append(specific_gemm_class(arch))
  else:
    print("YATETO::ERROR: unknown \"{}\" GEMM tool. "
          "Please, refer to the documentation".format(tool))
    sys.exit("failure")

# Generate code
gemmTools = GeneratorCollection(gemm_generators)
generator.generate(outputDir=cmdLineArgs.outputDir,
                   namespace='seissol',
                   gemm_cfg=gemmTools,
                   include_tensors=include_tensors)
