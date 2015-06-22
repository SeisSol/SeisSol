#!/bin/env python
##
# @file
# This file is part of SeisSol.
#
# @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
#
# @section LICENSE
# Copyright (c) 2014, SeisSol Group
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
# Performs a convergence analysis.
#

import argparse
import logging
import glob
import os
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
from matplotlib.backends.backend_pdf import PdfPages
import csv

# Gets the values of the norms defined in the given log file.
# @param i_pathToLog path to the log file.
def getValues( i_pathToLog ):
  # directory holding the results of this file
  l_results = {}

  l_fileContents = open( i_pathToLog, "r" )

  # regular expressions for the norms
  l_normExp = { 'float'    : '\s+[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?',
                'var'      : 'Rank:        0 \| Info    \| MPI Error analysis of variable  "',
                'l1'       : 'Rank:        0 \| Info    \| MPI L1_norm   :',
                'l2'       : 'Rank:        0 \| Info    \| MPI L2_norm   :',
                'linf'     : 'Rank:        0 \| Info    \| MPI Linf_norm :',
                'time'     : 'Rank:        0 \| Info    \| CPU-Time:',
                'timesteps': 'Rank:        0 \| Info    \| total number of performed time steps:' }

  # iterate over log file
  for l_line in l_fileContents:
    # check if this is a line describing a variable
    l_varSearch = re.search('(?<='+l_normExp['var']+')\w+', l_line )
    if l_varSearch:
      # store the variable
      l_variable = l_varSearch.group(0)

      # create a place in the result dictionary
      l_results[l_variable] = {}

    # check if this line describes a norm
    for l_norm in ['l1', 'l2', 'linf']:
      # check if this line describes a norm
      l_normSearch = re.search('(?<='+l_normExp[l_norm]+')'+l_normExp['float'], l_line )
      if l_normSearch:
        # store the norm
        l_results[l_variable][l_norm] = l_normSearch.group(0).strip()

    # check if this line describes the exceution tim
    l_timeSearch = re.search('(?<='+l_normExp['time']+')'+l_normExp['float'], l_line )
    if l_timeSearch:
      l_timeSearch = l_timeSearch.group(0)
      l_results['execution_time'] = float( l_timeSearch )

    # check if this line describes the number of time steps
    l_timeStepsSearch = re.search('(?<='+l_normExp['timesteps']+')'+l_normExp['float'], l_line )
    if l_timeStepsSearch:
      l_timeStepsSearch = l_timeStepsSearch.group(0)
      l_results['time_steps'] = int(l_timeStepsSearch)

  l_fileContents.close()

  return l_results

# set up logging level
logging.basicConfig( level=logging.DEBUG,
                     format='%(asctime)s - %(name)s - %(levelname)s - %(message)s' )

# command line interface
logging.info( "parsing command line arguments" )

l_parser    = argparse.ArgumentParser( description='Performs a convergence analysis.' )
l_parser.add_argument( '--log_dir',
                       dest     = "log_dir",
                       required = True,
                       help     = "Directory where the log files are located.",
                       metavar  = "LOG_DIR" )

l_parser.add_argument( '--log_regexp',
                       dest     = "log_regexp",
                       required = True,
                       help     = "Regular expression describing the format of the log files.\n\
                                   The three tags PRECISION_TAG, ORDER_TAG and CUBES_PER_DIM_TAG are required in\n\
                                   order to distinct between the runs.\n\
                                   Example \"release_generated_PRECISION_TAGknc_hybrid_none_Netcdf_9_ORDER_TAG_auto_periodic_cube_CUBES_PER_DIM_TAG_CUBES_PER_DIM_TAG_CUBES_PER_DIM_TAG_1_1_1_wsm_1_1_1_16.out\".",
                       metavar  = "LOG_REGEXP" )

l_parser.add_argument( '--bytes_per_cell',
                       dest     = "bytes_per_cell",
                       required = False,
                       help     = "Path to a CSV file, which specifies the bytes per cell.\n\
                                   Format:\n\
                                   oder,bytes_per_cell\\n1,9576\\n2,11088",
                       metavar  = "PATH_TO_CSV" )

l_parser.add_argument( '--output_dir',
                       dest     = "output_dir",
                       required = True,
                       help     = "Path to the output directory where the PDFs go.",
                       metavar  = "OUTPUT_DIR" )

l_arguments = vars(l_parser.parse_args())

#
# Parse files
#
logging.info( "parsing files" )

# get all files matching the regexp
l_files = glob.glob( l_arguments['log_dir']+'/'+l_arguments['log_regexp'].replace("PRECISION_TAG", "*").replace("ORDER_TAG", "*").replace("CUBES_PER_DIM_TAG", "*") )

# cut off directories
for l_file in range(len(l_files)):
  l_files[l_file] = os.path.basename(l_files[l_file])

# set up regular epressions
l_regularExp = { 'precision':     l_arguments['log_regexp'].replace( "PRECISION_TAG", "([sSdD]?)" ).replace( "CUBES_PER_DIM_TAG", "[0-9]+"   ).replace( "ORDER_TAG", "[0-9]+"   ),
                 'order':         l_arguments['log_regexp'].replace( "PRECISION_TAG", "[sSdD]"   ).replace( "CUBES_PER_DIM_TAG", "[0-9]+"   ).replace( "ORDER_TAG", "([0-9]+)" ),
                 'cubes_per_dim': l_arguments['log_regexp'].replace( "PRECISION_TAG", "[sSdD]"   ).replace( "CUBES_PER_DIM_TAG", "([0-9]+)" ).replace( "ORDER_TAG", "[0-9]+"   ) }

# directory containing all fieles sorted by order and cubes per dim
l_results = {}

# iterate over files and get order dimension and norms
for l_file in range(len(l_files)):
  # get order and number of cubes per dim
  l_precision = re.search(l_regularExp['precision'], l_files[l_file]).group(1)

  l_order     = re.search(l_regularExp['order'], l_files[l_file]).group(1)
  l_cubes     = re.search(l_regularExp['cubes_per_dim'], l_files[l_file]).group(1)
  l_order = int(l_order)
  l_cubes = int(l_cubes)

  # create a new directory for the precision
  if not l_precision in l_results:
    l_results[l_precision] = {}

  # create a new directory for new orders
  if not l_order in l_results[l_precision]:
    l_results[l_precision][l_order] = {}

  # get values
  l_values = getValues( l_arguments['log_dir']+'/'+l_files[l_file] )

  # add directory for this number of cubes per dimension (if log file is complete)
  if( 'execution_time' in l_values and 'w' in l_values ):
    l_results[l_precision][l_order][l_cubes] = l_values

# parse csv specifying bytes per cells
# TODO: add support for single precision
l_definitions = {}
if l_arguments['bytes_per_cell'] != None:
  l_definitions['bytes_per_cell'] = {}

  with open(l_arguments['bytes_per_cell'], 'r') as l_file:
    # jump over header
    l_file.next()

    # add all rows to the dictionary
    for l_row in l_file:
      l_row = l_row[:-1].split(",")
      l_definitions['bytes_per_cell'][int(l_row[0])] = int(l_row[1])
#
# Write csv files
#
logging.info( "writing csvs" )

# setup header
l_header = [ 'order', 'precision', 'cubes_per_dim', 'execution_time', 'time_steps',
             'u_l1', 'v_l1',   'w_l1',
             'sigma_xx_l1',   'sigma_yy_l1',   'sigma_zz_l1',
             'sigma_xy_l1',   'sigma_xz_l1',   'sigma_yz_l1',
             'u_l2',          'v_l2',          'w_l2',
             'sigma_xx_l2',   'sigma_yy_l2',   'sigma_zz_l2',
             'sigma_xy_l2',   'sigma_xz_l2',   'sigma_yz_l2',
             'u_linf',        'v_linf',        'w_linf',
             'sigma_xx_linf', 'sigma_yy_linf', 'sigma_zz_linf',
             'sigma_xy_linf', 'sigma_xz_linf', 'sigma_yz_linf'
           ]

for l_precision in l_results.keys():
  for l_order in l_results[l_precision]:
    # save to csv file
    with open( l_arguments['output_dir']+'/' + l_precision + str(l_order) + '.csv', 'w' ) as l_csvFile:
      l_writer = csv.DictWriter(l_csvFile, fieldnames=l_header)
      l_writer.writeheader()

      for l_cubes in l_results[l_precision][l_order].keys():
        # generate a flat temporary dictionary
        l_tempDict = {}
        l_tempDict['order']          = str(l_order)
        l_tempDict['precision']      = l_precision
        l_tempDict['cubes_per_dim']  = l_cubes
        l_tempDict['execution_time'] = l_results[l_precision][l_order][l_cubes]['execution_time']

        for l_var in ['u', 'v', 'w', 'sigma_xx', 'sigma_yy', 'sigma_zz', 'sigma_xy', 'sigma_xz', 'sigma_yz']:
          for l_norm in ['l1', 'l2', 'linf']:
            l_tempDict[l_var+'_'+l_norm] = l_results[l_precision][l_order][l_cubes][l_var][l_norm]

        l_writer.writerow( l_tempDict )

#
# Plot
#
logging.info( "plotting results" )

l_style = { 'markers': {  2: '^',
                          3: 'o',
                          4: '*',
                          5: 'p', 
                          6: 's',
                          7: 'v'       },
           'line':     { 'd': '-',
                         's': '--'     },
           'color':    { 'd': '#E37222',
                         's': '#0065BD'     },
           'alpha':    { 'd': 1,
                         's': 1
                                       }
          }

for l_variable in ['execution_time', 'GiB/s']:
  # jump over gibs, if no configuration is given
  if l_arguments['bytes_per_cell'] == None and l_variable == 'GiB/s':
    continue

  # define a new figure
  l_figure = matplotlib.pyplot.figure( figsize=(10, 7) )

  # generate plot showing time per time step
  if l_variable is 'execution_time':
    l_outputFile = PdfPages( l_arguments['output_dir']+'/time_per_update.pdf' )
  else:
    l_outputFile = PdfPages( l_arguments['output_dir']+'/gibs.pdf' )

  l_legend = []

  for l_precision in sorted(l_results.keys(), reverse=True):
    for l_order in sorted(l_results[l_precision].keys()):
      # mesh width and time
      l_h    = []
      l_values = []

      for l_cubesPerDim in sorted(l_results[l_precision][l_order].keys()):
        # only consider runs with positive CPU times
        if( l_results[l_precision][l_order][l_cubesPerDim]['execution_time'] < 0 ):
          continue

        l_h.append( l_cubesPerDim**3 * 5 )
        l_values.append(   l_results[l_precision][l_order][l_cubesPerDim]['execution_time'] / \
                           ( l_results[l_precision][l_order][l_cubesPerDim]['time_steps'] * l_cubesPerDim**3 * 5 )   )
        if l_variable is 'GiB/s':
          l_values[-1] = ( l_definitions['bytes_per_cell'][l_order] / l_values[-1] ) / 1024**3

      # plot the results
      matplotlib.pyplot.plot( l_h,
                              l_values,
                              linestyle = l_style['line'][l_precision],
                              color = l_style['color'][l_precision],
                              alpha = l_style['alpha'][l_precision],
                              label='reference',
                              marker=l_style['markers'][l_order] )

      # add labels to last point
      matplotlib.pyplot.annotate( "{:10.2e}".format( l_values[-1] ),
                                  xy = (l_h[-1], l_values[-1]),
                                  xytext = (-5, 2),
                                  ha='right',
                                  textcoords = 'offset points',
                                  fontsize = 6 )

      # add item to legend
      l_legend = l_legend + [l_precision + str(l_order)]

  # add legend
  matplotlib.pyplot.legend( l_legend )

  # use logarithmic scale
  #matplotlib.pyplot.xscale( 'log' )

  # use scientific notation on the y-axis
  matplotlib.pyplot.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

  # set labels
  matplotlib.pyplot.xlabel('number of tetrahedrons')

  if l_variable is 'execution_time':
    matplotlib.pyplot.ylabel('time per element update')
  else:
    matplotlib.pyplot.ylabel('GiB/s')

  # add a grid to the plot
  matplotlib.pyplot.grid( which='both' )

  # save the figure
  l_outputFile.savefig( l_figure )

  # close file
  l_outputFile.close()

# generate error plots
for l_parameter in ['mesh_width', 'execution_time']:
  # open the output file
  l_outputFile = PdfPages( l_arguments['output_dir']+'/'+l_parameter+'.pdf' )
    

  # iterate over norms
  for l_norm in ['l1', 'l2', 'linf']:
    # iterate over variables
    for l_variable in ['sigma_xx', 'sigma_yy', 'sigma_zz',
                       'sigma_xy', 'sigma_yz', 'sigma_xz',
                       'u',        'v',        'w' ]:
      # define a new figure
      l_figure = matplotlib.pyplot.figure( figsize=(10, 7) )

      # plot the resutls
      for l_precision in sorted(l_results.keys(), reverse=True):
        for l_order in sorted(l_results[l_precision].keys()):
          # values of the parameters for the x-axis
          l_xValues = []
          l_errors = []

          for l_cubesPerDim in sorted(l_results[l_precision][l_order].keys()):
            if( l_parameter == 'mesh_width' ):
              l_xValues.append( 1.0 / l_cubesPerDim )
            elif( l_parameter == 'execution_time' ):
              # continue for invalid execution times
              if( l_results[l_precision][l_order][l_cubesPerDim]['execution_time'] < 0 ):
                continue

              l_xValues.append( l_results[l_precision][l_order][l_cubesPerDim][l_parameter] )

            l_errors.append( l_results[l_precision][l_order][l_cubesPerDim][l_variable][l_norm] )

          # add item to legend
          l_legend = l_legend + [l_precision + str(l_order)]

          # plot the results
          matplotlib.pyplot.plot( l_xValues,
                                  l_errors,
                                  linestyle = l_style['line'][l_precision],
                                  color = l_style['color'][l_precision],
                                  alpha = l_style['alpha'][l_precision],
                                  label='reference',
                                  marker=l_style['markers'][l_order] )

        # set title
        matplotlib.pyplot.title( l_variable )

      # add legend
      matplotlib.pyplot.legend( l_legend )

      # use logarithmic scale
      matplotlib.pyplot.xscale( 'log' )
      matplotlib.pyplot.yscale( 'log' )

      # set labels
      if( l_parameter == 'mesh_width' ):
        matplotlib.pyplot.xlabel('mesh width (cubes)')
      elif( l_parameter == 'execution_time' ):
        matplotlib.pyplot.xlabel('execution time (s)')
      matplotlib.pyplot.ylabel('error '+ l_norm)

      # add a grid to the plot
      matplotlib.pyplot.grid( which='both' )

      # save the figure
      l_outputFile.savefig( l_figure )

      # close figure
      matplotlib.pyplot.close( l_figure )

  # close file
  l_outputFile.close()
