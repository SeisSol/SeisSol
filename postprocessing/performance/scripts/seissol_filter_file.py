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
# Summarizes the performance output of multiple Jobs.
#

import argparse
import os
import re
import itertools

l_commandLineParser = argparse.ArgumentParser( description=
  'This simple script parse a couple of job logs.'
)

l_commandLineParser.add_argument( 'jobLogDirectory',
                                  #metavar='jobLogDirectory',
                                  type=str, nargs='+',
                                  help='directory, which contains the job logs (ending with *.out).')

l_commandLineArguments = l_commandLineParser.parse_args()

print 'given directory: ' + l_commandLineArguments.jobLogDirectory[0]

l_mergedFilePath = l_commandLineArguments.jobLogDirectory[0] + '/merged.csv'
print 'creating ' + l_mergedFilePath

l_mergedFile = open(l_mergedFilePath, 'w+')

l_logFiles = os.listdir(l_commandLineArguments.jobLogDirectory[0])
l_logFiles.sort()

# define column names in the job scripts
l_columnNames = ['MPI rank', 'code version', 'kernel', '#basis functions', '#kernel calls', 'CPU cycles (mean), execution time (mean)']

# write file header for the merged file
l_mergedFile.write( 'job' )
for l_columnName in l_columnNames:
  l_mergedFile.write( ',' + l_columnName )
l_mergedFile.write( '\n' )

# columns which have to be identical to be part of the same measurement
l_csvFileLayout = [True,False,True,True,True,False,False,False]

# setup measurement layout
l_measurementLayout = [not x for x in l_csvFileLayout]
# ignore MPI rank
l_measurementLayout[1] = False

for l_file in l_logFiles:

  if l_file.endswith('.out'):
    print '  processing ', l_file

    # create job sum
    l_jobSum = {}

    l_pathToLogFile = l_commandLineArguments.jobLogDirectory[0] + l_file

    l_logFile = open(l_pathToLogFile, 'r')

    for l_line in l_logFile:
      # search for performance measurements
      if re.search('^[0-9]', l_line):
        # print with original file name in front
        l_mergedFile.write( l_file + ',' + l_line )

        # split this measurement into its component and add the job name in front
        l_splittedMeasurement = [l_file] +l_line[0:-1].split(',')

        # generate a unique key for this measurement
        l_multiKey = ()
        for l_key in itertools.compress(l_splittedMeasurement, l_csvFileLayout):
          l_multiKey = l_multiKey + (l_key,)

        # get the measurements
        l_measurements = []

        for l_measurement in itertools.compress(l_splittedMeasurement, l_measurementLayout):
          l_measurements = l_measurements + [float(l_measurement)]

        # Add the measurement to the dictionary
        if l_multiKey in l_jobSum:
          # update the current measurements
          for l_index in xrange(len(l_measurements)):
            l_jobSum[l_multiKey][l_index] = l_jobSum[l_multiKey][l_index] + l_measurements[l_index]
        else:
          # add the key
          l_jobSum[l_multiKey] = l_measurements
      
    # write job sum to merged file
    for l_keys in l_jobSum.keys():
      # write keys
      l_mergedFile.write( l_keys[0] + ',sum_all_ranks' )
      for l_key in l_keys[1:]:
        l_mergedFile.write( ',' + l_key )
      
      # write values
      for l_value in l_jobSum[l_keys]:
        l_mergedFile.write( ',' + str(l_value) )
      l_mergedFile.write( '\n' )

    l_logFile.close()

l_mergedFile.close()
