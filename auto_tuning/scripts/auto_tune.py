#!/usr/bin/env python
##
# @file
# This file is part of SeisSol.
#
# @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
#
# @section LICENSE
# Copyright (c) 2016, SeisSol Group
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
import os
from tune import MemoryLayout, Proxy, Analysis

cmdLineParser = argparse.ArgumentParser()
cmdLineParser.add_argument('--equations', required=True)
cmdLineParser.add_argument('--order', required=True, type=int)
cmdLineParser.add_argument('--numberOfMechanisms', type=int, default=0)
cmdLineParser.add_argument('--arch', required=True)
cmdLineParser.add_argument('--workingDir', required=True)
cmdLineParser.add_argument('--nelem', default=10000, type=int)
cmdLineParser.add_argument('--ntimesteps', default=100, type=int)
cmdLineParser.add_argument('--ncompileJobs', default=4, type=int)
cmdLineParser.add_argument('--host', default='')
cmdLineParser.add_argument('action', choices=['build', 'test', 'run', 'analyse', 'all'])
args = cmdLineParser.parse_args()

build = args.action == 'build' or args.action == 'all'
test = args.action == 'test' or args.action == 'all'
run = args.action == 'run' or args.action == 'all'
analyse = args.action == 'analyse' or args.action == 'all'

if args.numberOfMechanisms == 0 and args.equations.startswith('viscoelastic'):
  raise ValueError('The number of mechanisms must be greater than 0 for equations=viscoelastic.')


# Generate working directory
try:
  os.mkdir(args.workingDir)
except OSError:
  print('Warning: Working directory does already exist. Previous results may be overwritten.')
  
os.chdir(args.workingDir)

def tryMkdir(dirName):
  try:
    os.mkdir(dirName)
  except OSError:
    pass

tryMkdir(MemoryLayout.OutputDir)
tryMkdir(Proxy.BuildDir)
tryMkdir(Proxy.OutputDir)

if build or analyse:
  # Generate memory layouts  
  if args.equations == 'elastic':
    memoryLayouts = MemoryLayout.getElasticMemoryLayouts(args.order, args.arch)
  if args.equations == 'viscoelastic':
    memoryLayouts = MemoryLayout.getViscoelasticMemoryLayouts(args.order, args.numberOfMechanisms, args.arch)
  elif args.equations == 'viscoelastic2':
    memoryLayouts = MemoryLayout.getViscoelastic2MemoryLayouts(args.order, args.arch)

if build:
  MemoryLayout.generateTuningLayoutFiles(memoryLayouts)
  
# Build versions, run tests and proxy
options = {
  'equations':          args.equations,
  'order':              args.order,
  'numberOfMechanisms': args.numberOfMechanisms,
  'arch':               args.arch,
  'compileMode':        'release',
  '--jobs':             args.ncompileJobs
}
Proxy.buildOrRun(options, args.nelem, args.ntimesteps, build=build, test=test, run=run, host=args.host)

if analyse:
  # Analysis
  matrices = Analysis.analyse()

  # Generate best memory layout
  MemoryLayout.generateLayoutFile(matrices, memoryLayouts)
