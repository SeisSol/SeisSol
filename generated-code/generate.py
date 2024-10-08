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
import os

from yateto import useArchitectureIdentifiedBy, Generator, NamespacedGenerator
from yateto import gemm_configuration
from yateto.gemm_configuration import GeneratorCollection, Eigen, LIBXSMM_JIT, PSpaMM, MKL, BLIS, OpenBLAS, GemmForge
from yateto.ast.cost import BoundingBoxCostEstimator, FusedGemmsBoundingBoxCostEstimator

import kernels.DynamicRupture as DynamicRupture
import kernels.Plasticity as Plasticity
import kernels.SurfaceDisplacement as SurfaceDisplacement
import kernels.Point as Point
import kernels.NodalBoundaryConditions as NodalBoundaryConditions
import kernels.memlayout as memlayout
import kernels.vtkproject as vtkproject
import kernels.general as general

def main():
  cmdLineParser = argparse.ArgumentParser()
  cmdLineParser.add_argument('--equations')
  cmdLineParser.add_argument('--matricesDir')
  cmdLineParser.add_argument('--outputDir')
  cmdLineParser.add_argument('--host_arch')
  cmdLineParser.add_argument('--device_backend', default=None)
  cmdLineParser.add_argument('--device_arch', default=None)
  cmdLineParser.add_argument('--order', type=int)
  cmdLineParser.add_argument('--numberOfMechanisms', type=int)
  cmdLineParser.add_argument('--memLayout')
  cmdLineParser.add_argument('--multipleSimulations', type=int)
  cmdLineParser.add_argument('--PlasticityMethod')
  cmdLineParser.add_argument('--gemm_tools')
  cmdLineParser.add_argument('--drQuadRule')
  cmdLineParser.add_argument('--enable_premultiply_flux', action='store_true')
  cmdLineParser.add_argument('--disable_premultiply_flux', dest='enable_premultiply_flux', action='store_false')
  cmdLineParser.add_argument('--executable_libxsmm', default='')
  cmdLineParser.add_argument('--executable_pspamm', default='')
  cmdLineParser.set_defaults(enable_premultiply_flux=False)
  cmdLineArgs = cmdLineParser.parse_args()

  # derive the compute platform
  gpu_platforms = ['cuda', 'hip', 'hipsycl', 'oneapi']
  targets = ['gpu', 'cpu'] if cmdLineArgs.device_backend in gpu_platforms else ['cpu']

  subfolders = []

  if cmdLineArgs.memLayout == 'auto':
    # TODO(Lukas) Don't hardcode this
    env = {
      'equations': cmdLineArgs.equations,
      'order': cmdLineArgs.order,
      'arch': cmdLineArgs.host_arch,
      'device_arch': cmdLineArgs.device_arch,
      'multipleSimulations': cmdLineArgs.multipleSimulations,
      'targets': targets
    }
    mem_layout = memlayout.guessMemoryLayout(env)
  else:
    mem_layout = cmdLineArgs.memLayout


  if cmdLineArgs.device_backend == 'none':
      arch = useArchitectureIdentifiedBy(cmdLineArgs.host_arch)
  else:
      arch = useArchitectureIdentifiedBy(cmdLineArgs.host_arch, cmdLineArgs.device_arch, cmdLineArgs.device_backend)


  equationsSpec = importlib.util.find_spec(f'kernels.{cmdLineArgs.equations}')
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
  elif cmdLineArgs.equations == 'viscoelastic2':
      adg = equations.Viscoelastic2ADERDG(**cmdArgsDict)
  else:
      adg = equations.PoroelasticADERDG(**cmdArgsDict)

  include_tensors = set()
  generator = Generator(arch)

  # Equation-specific kernels
  adg.addInit(generator)
  adg.addLocal(generator, targets)
  adg.addNeighbor(generator, targets)
  adg.addTime(generator, targets)
  adg.add_include_tensors(include_tensors)

  vtkproject.addKernels(generator, adg, cmdLineArgs.matricesDir, targets)
  vtkproject.includeTensors(cmdLineArgs.matricesDir, include_tensors)

  # Common kernels
  include_tensors.update(DynamicRupture.addKernels(NamespacedGenerator(generator, namespace="dynamicRupture"),
                                                  adg,
                                                  cmdLineArgs.matricesDir,
                                                  cmdLineArgs.drQuadRule,
                                                  targets))

  Plasticity.addKernels(generator, 
                        adg,
                        cmdLineArgs.matricesDir,
                        cmdLineArgs.PlasticityMethod,
                        targets)
  NodalBoundaryConditions.addKernels(generator, adg, include_tensors, cmdLineArgs.matricesDir, cmdLineArgs, targets)
  SurfaceDisplacement.addKernels(generator, adg, include_tensors, targets)
  Point.addKernels(generator, adg)

  # pick up the user's defined gemm tools
  gemm_tool_list = cmdLineArgs.gemm_tools.replace(" ", "").split(",")
  gemm_generators = []

  for tool in gemm_tool_list:
    if hasattr(gemm_configuration, tool):
      specific_gemm_class = getattr(gemm_configuration, tool)
      # take executable arguments, but only if they are not empty
      if specific_gemm_class is gemm_configuration.LIBXSMM and cmdLineArgs.executable_libxsmm != '':
        gemm_generators.append(specific_gemm_class(arch, cmdLineArgs.executable_libxsmm))
      elif specific_gemm_class is gemm_configuration.PSpaMM and cmdLineArgs.executable_pspamm != '':
        gemm_generators.append(specific_gemm_class(arch, cmdLineArgs.executable_pspamm))
      else:
        gemm_generators.append(specific_gemm_class(arch))
    elif tool.strip().lower() != 'none':
      print("YATETO::ERROR: unknown \"{}\" GEMM tool. "
            "Please, refer to the documentation".format(tool))
      sys.exit("failure")


  cost_estimators = BoundingBoxCostEstimator
  if 'gpu' in targets and cmdLineArgs.equations == 'elastic':
    try:
      chainforge_spec = importlib.util.find_spec('chainforge')
      chainforge_spec.loader.load_module()
      cost_estimators = FusedGemmsBoundingBoxCostEstimator
    except:
      print('WARNING: ChainForge was not found. Falling back to GemmForge.')

  precision = 'double' if cmdLineArgs.host_arch[0] == 'd' else 'single'
  outputDirName = f'equation-{cmdLineArgs.equations}-{cmdLineArgs.order}-{precision}'
  trueOutputDir = os.path.join(cmdLineArgs.outputDir, outputDirName)
  if not os.path.exists(trueOutputDir):
    os.mkdir(trueOutputDir)

  subfolders += [outputDirName]

  # Generate code
  gemmTools = GeneratorCollection(gemm_generators)
  generator.generate(outputDir=trueOutputDir,
                    namespace='seissol',
                    gemm_cfg=gemmTools,
                    cost_estimator=cost_estimators,
                    include_tensors=include_tensors)

  def generate_general(subfolders):
    # we use always use double here, since these kernels are only used in the initialization
    arch = useArchitectureIdentifiedBy('d' + cmdLineArgs.host_arch[1:])

    outputDir = os.path.join(cmdLineArgs.outputDir, 'general')
    if not os.path.exists(outputDir):
      os.mkdir(outputDir)

    subfolders += [f'general']

    # for now, enforce Eigen as a code generator here... Until we have a shared subroutine cache
    generator = Generator(arch)
    general.addStiffnessTensor(generator)
    generator.generate(outputDir=outputDir,
                      namespace='seissol_general',
                      gemm_cfg=GeneratorCollection([Eigen(arch)]),
                      cost_estimator=cost_estimators,
                      include_tensors=general.includeMatrices(cmdLineArgs.matricesDir))

  generate_general(subfolders)

  def forward_files(filename):
    with open(os.path.join(cmdLineArgs.outputDir, filename), 'w') as file:
      file.writelines(['// IWYU pragma: begin_exports\n'])
      file.writelines([f'#include "{os.path.join(folder, filename)}"\n' for folder in subfolders])
      file.writelines(['// IWYU pragma: end_exports\n'])

  forward_files('init.h')
  forward_files('kernel.h')
  forward_files('subroutine.h')
  forward_files('tensor.h')
  forward_files('init.cpp')
  forward_files('kernel.cpp')
  forward_files('subroutine.cpp')
  forward_files('tensor.cpp')
  forward_files('gpulike_subroutine.cpp')

if __name__ == '__main__':
  main()
