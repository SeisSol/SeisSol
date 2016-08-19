#!/usr/bin/python
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
# @section DESCRIPTION
#
# Compares receivers in two different directories against each other.
#

import argparse
import os
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
from matplotlib.backends.backend_pdf import PdfPages

# set up command line parser and get job log
l_commandLineParser = argparse.ArgumentParser( description=
  'This simple script sums up the calls of flux matrices in SeisSol.'
)

l_commandLineParser.add_argument( 'receiverDirectories',
                                  type=str, nargs=2,
                                  help='directories which hold the reveiver outputs')

l_commandLineParser.add_argument( 'outputFiles',
                                  type=str, nargs=2,
                                  help='PDF-file where the plots are stores (will be overwritten) and csv-file for the misfits (will be appended)')

l_commandLineArguments = l_commandLineParser.parse_args()


#! dictionaries containing information about the date
l_reference = {}
l_simulation = {}

# get the directories containing the receiver outputs
l_reference['directory'] = l_commandLineArguments.receiverDirectories[0]
l_simulation['directory'] = l_commandLineArguments.receiverDirectories[1]

# get the output files
l_output = {}
l_output['pdf'] = l_commandLineArguments.outputFiles[0]
l_output['csv'] = l_commandLineArguments.outputFiles[1]

print 'given input directories: ' + l_reference['directory']  + ' and ' + l_simulation['directory']
print 'given output file: ' + l_output['pdf']

# get and sort the files
l_reference['files'] = os.listdir(l_reference['directory'])
l_reference['files'].sort()
l_simulation['files'] = os.listdir(l_simulation['directory'])
l_simulation['files'].sort()

# remove non-receiver data
for l_files in l_reference['files'], l_simulation['files']:
  for l_file in l_files:
    if( 'receiver' not in  l_file ):
      l_files.remove( l_file )

# assert each file of the reference is covered by the simulation
for l_file in l_reference['files']:
  assert( l_file in l_simulation['files'] )

# open the output files
l_output['pdf'] = PdfPages( l_output['pdf'] )

if( os.path.exists(l_output['csv']) ):
  l_output['csv'] = open( l_output['csv'], 'a' )
else:
  # create a new file with a header if the file does not exist already
  l_output['csv'] = open( l_output['csv'], 'w' )
  l_output['csv'].write( 'reference_directory,simulation_directory,receiver,variable,seismogram_misfit\n' )

# iterate over the files and plot them
for l_file in l_reference['files']:
  # get the reference and simulation's solution of this receiver
  for l_solution in l_reference, l_simulation:
    # reset the current solution
    l_solution['current'] = {'time': [],  'u': [],  'v': [],  'w': [],
                                         'xx': [], 'yy': [], 'zz': [],
                                         'xy': [], 'yz': [], 'xz': [],
                                         'rx': [], 'ry': [], 'rz': [],
                                         'srs': [], 'srd': [],
                                         't_s': [], 't_d': [],
                                         'p_n': [], 'u_n': []          }

    # get solution
    for l_row in open(l_solution['directory'] + '/' + l_file ):
      l_row = l_row.split()
      
      # jump over the title
      if( l_row[0] == "TITLE" ):
        continue
      # found the variables
      elif( l_row[0] == "VARIABLES" ):
        # read and format the variables
        l_solution['current']['variables'] = ''.join(l_row[2:])
        l_solution['current']['variables'] =  l_solution['current']['variables'].replace('"', '').lower().split(',')
        continue
      # continue in case of an header row
      elif (len(l_row) != len( l_solution['current']['variables'] ) ):
        continue

      for l_variableIndex in xrange( len(l_solution['current']['variables']) ):
        l_variable = l_solution['current']['variables'][l_variableIndex]      
        l_solution['current'][l_variable] = l_solution['current'][l_variable] + [float(l_row[l_variableIndex])]

  # iterate over variables
  for l_variable in l_solution['current']['variables'][1:]:
    #
    # compute seismogram misfit
    #
    l_firstSum = 0
    l_secondSum = 0

    for l_index in xrange( min( len(l_reference['current']['time']), len(l_simulation['current']['time']) ) ):
      l_firstSum = l_firstSum + (l_simulation['current'][l_variable][l_index] - l_reference['current'][l_variable][l_index])**2
      l_secondSum = l_secondSum + (l_reference['current'][l_variable][l_index])**2
    
    if( l_secondSum > 0 ):
      absoluteError = False
      l_seismogramMisfit = l_firstSum / l_secondSum
    else:
      absoluteError = True
      l_seismogramMisfit = l_firstSum

    l_output['csv'].write( l_reference['directory']  + ',' +
                           l_simulation['directory'] + ',' +
                           l_file                    + ',' +
                           l_variable                + ',' +
                           str(l_seismogramMisfit)   + '\n' )

    #
    # plot the reference against the simulation
    #

    # define a new figure
    l_figure = matplotlib.pyplot.figure( figsize=(20, 10) )

    # setup the plot
    matplotlib.pyplot.title('{}, {}seismogram misfit={}'.format(l_file, '(absolute) ' if absoluteError else '', l_seismogramMisfit))
    matplotlib.pyplot.suptitle('reference: ' + l_reference['directory'] + '\n simulation: ' + l_simulation['directory'] )
    matplotlib.pyplot.xlabel('time')
    matplotlib.pyplot.ylabel(l_variable)

    # plot both solutions
    matplotlib.pyplot.plot( l_reference['current']['time'],  l_reference['current'][l_variable], 's', mfc='none', label='reference' )
    matplotlib.pyplot.plot( l_simulation['current']['time'], l_simulation['current'][l_variable], 'x', label='test' )

    # set legend
    matplotlib.pyplot.legend( ['reference', 'simulation'] )

    l_output['pdf'].savefig( l_figure )

# close files
l_output['pdf'].close()
l_output['csv'].close()
