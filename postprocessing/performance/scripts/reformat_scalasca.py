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
# Reformats cube_stat output to be readable by spreadsheets.
#
import argparse

# set up command line parser
l_commandLineParser = argparse.ArgumentParser( description=
  'This simple script reformats cube_stat output for spreadsheets.'
)

l_commandLineParser.add_argument( 'mergedCubeFile',
                                  type=str, nargs=1,
                                  help='file, which holds merged output of Cube')

l_commandLineParser.add_argument( 'outputFile',
                                  type=str, nargs=1,
                                  help='file, where the output will be written to')

l_commandLineArguments = l_commandLineParser.parse_args()

# get the file containing cube output
l_csvFile = l_commandLineArguments.mergedCubeFile[0]
l_outputFile = l_commandLineArguments.outputFile[0]

print 'processing ' + l_csvFile + '..'

# open the cube csv file
l_csvFile = open( l_csvFile, 'r' )

# open output file
l_outputFile = open( l_outputFile, 'w' )

# single measurement
l_measurement = ''

# iterate over rows
for l_row in l_csvFile:
  if( l_row == '\n' ): # ignore empty rows
    continue
  else: # append new routine of the current measurement
    # get row description
    l_rowDescription = l_row.split(',')[0]

    if l_rowDescription == 'directory': # next measurement
      # write to output file and remove last colon
      l_outputFile.write(l_measurement[:-1]+'\n')
      l_measurement = l_row.split(',')[1].replace('\n', '') + ','
    elif l_rowDescription == 'time':
      l_measurement = l_measurement + l_row.split(',')[4].replace('\n', '') + ','
