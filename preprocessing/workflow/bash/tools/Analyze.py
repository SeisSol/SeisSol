#!/bin/env python
##
# @file
# This file is part of SeisSol.
#
# @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
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
# Sets up the benchmark runs.
#
import logging

## Analyzation phase of the benchmark runs.
class Analyze:
  def __init__( self ):
    self.m_logger = logging.getLogger('Setup')

  ## Prepares the directories for the analyzed output.
  #
  # @param i_benchmarkRuns benchmarks.
  # @param i_builds builds configurations.
  def setupOutputDirectories( self,
                              i_benchmarkRuns,
                              i_builds ):
    l_bashCommands = '\n# preparing output directories\n'
    l_bashCommands = l_bashCommands + "echo 'preparing output directories'\n"

    l_bashCommands = l_bashCommands + 'mkdir $WORKING_DIRECTORY/output/\n'

    for l_build in i_builds:
      l_bashCommands = l_bashCommands+\
                       'mkdir $WORKING_DIRECTORY/output/' + l_build['id'] + '\n'
      for l_benchmark in i_benchmarkRuns:
        l_bashCommands = l_bashCommands+\
                         'mkdir $WORKING_DIRECTORY/output/' + l_build['id'] + '/' + l_benchmark['id'] + '\n'
        l_bashCommands = l_bashCommands+\
                         'mkdir $WORKING_DIRECTORY/output/' + l_build['id'] + '/' + l_benchmark['id'] + '/receivers' + '\n'
    return l_bashCommands

  ## Prepares the receivers for comparisons
  #
  # @param i_type type of the workflow.
  # @param i_benchmarkRuns benchmarks.
  # @param i_builds build configurations.
  def prepareReceivers( self,
                        i_type,
                        i_benchmarkRuns,
                        i_builds ):
    self.m_logger.info( "preparing receivers" )

    l_bashCommands = '\n# preparing receivers\n'
    l_bashCommands = l_bashCommands + "echo 'preparing receivers'\n"

    for l_build in i_builds:
      for l_benchmark in i_benchmarkRuns:
        l_bashArguments = " -i " + "$WORKING_DIRECTORY/runs/"   + l_build['id'] + "/" + l_benchmark['id'] + "/output"+\
                          " -o " + "$WORKING_DIRECTORY/output/" + l_build['id'] + "/" + l_benchmark['id'] + "/receivers"

        l_bashCommands = l_bashCommands+\
                         "sh ${SCRIPTS_DIRECTORY}/analyze/remove_ranks.sh" + l_bashArguments +"\n"
    return l_bashCommands

  ## Compares the receivers with a reference solution.
  #
  # @param i_benchmarkRuns benchmark runs to compare.
  def compareReceivers( self,
                        i_benchmarkRuns,
                        i_builds ):
    self.m_logger.info( "comparing receivers" )
    l_bashCommands = "echo 'comparing receivers'\n"

    for l_build in i_builds:
      for l_benchmark in i_benchmarkRuns:
        l_pythonArguments = "$INPUT_DIRECTORY/benchmarks/" + l_benchmark['name'] + "/references/" + l_build['order']+ " "+\
                            "$WORKING_DIRECTORY/output/"   + l_build['id'] + "/" + l_benchmark['id'] + "/receivers" + " "+\
                            "$WORKING_DIRECTORY/output/"   + l_build['id'] + "/" + l_benchmark['id'] + "/plots.pdf" + " "+\
                            "$WORKING_DIRECTORY/output/"   + l_build['id'] + "/" + l_benchmark['id'] + "/misfits.csv"
        l_bashCommands    = l_bashCommands + 'python ${SCRIPTS_DIRECTORY}/analyze/compare_receivers.py '  + l_pythonArguments + " & \n"
    return l_bashCommands

  # Extracts the data of the log files from convergence runs.
  #
  # @param i_regularExpression regular expression for the log files.
  def extractLogData( self,
                      i_regularExpressions ):
    l_bashCommands = ''

    # iterate over the regular expressions and genrate plots for each
    for l_regularExpression in i_regularExpressions:
      l_bashCommands = l_bashCommands+\
                       'mkdir -p $WORKING_DIRECTORY/output/' + l_regularExpression + '\n'
      l_pythonArguments = "--log_dir=$WORKING_DIRECTORY/logs/ "+\
                          "--log_regexp="+l_regularExpression+".out "+\
                          "--output_dir=$WORKING_DIRECTORY/output/" + l_regularExpression
      l_bashCommands = l_bashCommands + 'python ${SCRIPTS_DIRECTORY}/analyze/convergence.py '  + l_pythonArguments + " & \n"
    return l_bashCommands

  ## Preprocessing of the benchmark runs.
  #
  # @param i_type type of the workflow.
  # @param i_benchmarks benchmarks.
  # @param i_builds build configurations.
  def postProcess( self,
                   i_type,
                   i_benchmarks,
                   i_builds,
                   i_regularExpressions = set() ):
    self.m_logger.info( "postprocessing benchmarks" )

    #
    # bash interface
    #
    if( i_type == 'bash' ): 
      l_bashCommands ='''
#
# Postprocessing benchmarks
#
if [ ${ANALYZE} != "0" ]
then

echo "$(date) postprocessing benchmarks"
'''
      # process runs with receivers
      l_receiverRuns = [l_entry for l_entry in i_benchmarks['runs'] if l_entry['name'] != 'periodic']
      if( len(l_receiverRuns) != 0 ):
        l_bashCommands = l_bashCommands+\
                         self.setupOutputDirectories( i_benchmarkRuns = i_benchmarks['runs'],
                                                      i_builds        = i_builds )

        l_bashCommands = l_bashCommands+\
                         self.prepareReceivers( i_type           = i_type,
                                                i_benchmarkRuns = i_benchmarks['runs'],
                                                i_builds         = i_builds )

        l_bashCommands = l_bashCommands+\
                         self.compareReceivers( i_benchmarkRuns = i_benchmarks['runs'],
                                                i_builds        = i_builds )

      # process convergence runs
      if( len(i_regularExpressions) > 0 ):
        l_bashCommands = l_bashCommands + self.extractLogData( i_regularExpressions )

      # close analysis
      l_bashCommands = l_bashCommands+\
                       "\n echo \'waiting for postprocessing jobs to finish\' \n"+\
                       "wait\n"+\
                       "\nfi\n"+\
                       "# finished postprocessing benchmarks\n"

      return l_bashCommands
