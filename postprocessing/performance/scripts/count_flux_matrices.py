#!/bin/env python
##
# @file
# This file is part of SeisSol.
#
# @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
#
# @section LICENSE
# Copyright (c) 2013, SeisSol Group
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
# Sums up the the usage of flux matrices in SeisSol.
#

import argparse
import os
import re

# set up command line parser and get job log
l_commandLineParser = argparse.ArgumentParser( description=
  'This simple script sums up the calls of flux matrices in SeisSol.'
)

l_commandLineParser.add_argument( 'pathToJobLog',
                                  type=str, nargs='+',
                                  help='path to the job log')

l_commandLineArguments = l_commandLineParser.parse_args()

l_pathToLogFile = l_commandLineArguments.pathToJobLog[0]

print '#given job log: ' + l_pathToLogFile

# open log gile
l_logFile = open(l_pathToLogFile, 'r')

#! our flux matrix counter, intialized to zero
l_fluxMatrixCounter = [0]*52

# iterate over lines in the log file
for l_line in l_logFile:
  # search for matrix counter lines
  if re.search('Debug(\s)+\|(\s)+(\d)+.*', l_line):
    # get measurements
    l_line = l_line.split()
    
    # get matrix index and #calls for this rank
    #   in case of ignored boundary conditions increase the index by four: We don't have local flux matrices
    l_index = int(l_line[-2])
    l_numberOfCalls = int(l_line[-1])

    # add to the counter
    l_fluxMatrixCounter[l_index] = l_fluxMatrixCounter[l_index] + l_numberOfCalls
l_logFile.close()

# return the result
print '#printing flux matrix calls'
for l_index in range(len(l_fluxMatrixCounter)):
  print str(l_index) + ',' + str(l_fluxMatrixCounter[l_index])
