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

import MemoryLayout
import os
import time
import random
import string

BuildDir = 'builds'
OutputDir = 'results'

def _run(cmd, host=''):
  syscmd = 'ssh {} "{}"'.format(host, cmd) if host else cmd
  print(syscmd)
  if os.system(syscmd) != 0:
    raise Exception('Last system command failed.')

def buildOrRun(proxyOptions, nElem, nTimesteps, build=True, test=True, run=True, host=''):
  layoutDir = MemoryLayout.OutputDir
  buildDir = os.path.abspath(BuildDir)
  outputDir = os.path.abspath(OutputDir)
  configFiles = [os.path.abspath(os.path.join(layoutDir, filename)) for filename in os.listdir(layoutDir)]
  
  cwd = os.getcwd()
  os.chdir(os.path.dirname(os.path.realpath(__file__)) + '/../../proxy/')

  for configFile in configFiles:
    configName = os.path.basename(configFile).split('.')[0]
    options = ' '.join(['{}={}'.format(key, value) for key, value in proxyOptions.items()])
    configBuildDir = '{}/build_{}'.format(buildDir, configName)
    if build:
      _run( 'scons {} memLayout={} buildDir={} > {}/{}.build'.format(
                    options,
                    configFile,
                    configBuildDir,
                    outputDir,
                    configName ), host)
    if test:
      _run( '{}/bin/generated_kernels_test_suite > {}/{}.test'.format(configBuildDir, outputDir, configName), host)
    if run:
      # Create a random string in order to prevent file name collision
      outputName = '{}_{}_{}'.format(configName, time.strftime('%Y-%m-%d-%H-%M-%S'), ''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(5)))
      _run('{}/bin/seissol_proxy {} {} all > {}/{}.run'.format(configBuildDir, nElem, nTimesteps, outputDir, outputName), host)
  
  os.chdir(cwd)

