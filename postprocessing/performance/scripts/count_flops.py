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
# Sums up the total number of FLOPs over multiple MPI-ranks.
#

import argparse
import os
import re

# set up command line parser and get job log
l_commandLineParser = argparse.ArgumentParser( description=
  'This simple script sums of the number of FLOPs over multiple MPI-ranks measured with kernel generation.'
)

l_commandLineParser.add_argument( 'pathToJobLog',
                                  #metavar='pathToJobLog',
                                  type=str, nargs='+',
                                  help='path to the job log')

l_commandLineArguments = l_commandLineParser.parse_args()

l_pathToLogFile = l_commandLineArguments.pathToJobLog[0]

print 'given job log: ' + l_pathToLogFile

# open log gile
l_logFile = open(l_pathToLogFile, 'r')

# define dictionary for the FLOP measurements
l_measurements = {'time integration': 0, 'volume integration': 0, 'boundary integration': 0}

# iterate over lines in the log file
for l_line in l_logFile:
  # search for FLOP measurements
  if re.search('.*Info.*FLOPS', l_line):
    # get measurements
    l_line = l_line.split('-')[1]

    # get the measurements of this core
    l_coreMeasurements = [token for token in l_line.split() if token.isdigit()]

    # add measurements of this core to the overall FLOP count
    l_measurements['time integration']     = l_measurements['time integration']     +  int(l_coreMeasurements[0])
    l_measurements['volume integration']   = l_measurements['volume integration']   +  int(l_coreMeasurements[1])
    l_measurements['boundary integration'] = l_measurements['boundary integration'] +  int(l_coreMeasurements[2])
l_logFile.close()

# return the result
print 'printing total flop numbers over all ranks'
print l_measurements
