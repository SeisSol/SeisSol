#!/bin/env python
##
# @file
# This file is part of SeisSol.
#
# @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
#
# @author Leonhard Rannabauer (lrannabauer AT mytum.de)
#
# @section LICENSE
# Copyright (c) 2014-2015, SeisSol Group
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
# Executes the benchmark runs.

import logging
import os

## Execution phase of the benchmark runs.
class Run:
  def __init__( self ):
    self.m_logger = logging.getLogger('Run')



  def submit( self,
              i_type,
              i_resource_management,
              i_benchmarks,
              i_builds ):
    if i_resource_management == "loadleveler":
      return self.submitLoadLeveler(i_type,i_benchmarks,i_builds)
    elif i_resource_management == "slurm":
      return self.submitSlurm(i_type,i_benchmarks,i_builds)
    else:
      raise ValueError("Invalid resource management: "+i_resource_management+". Must be loadleveler or slurm")


  ## Submits the runs to the batch system.
  #
  # @param i_type type of the workflow.
  # @param i_benchmarks benchmarks.
  # @param i_builds build configurations.
  def submitSlurm( self,
                   i_type,
                   i_benchmarks,
                   i_builds ):
    self.m_logger.info( "submitting benchmarks" )

    if( i_type == 'bash' ): 
      l_bashCommands ='''
#
# Submitting benchmarks
#
if [ ${BENCHMARKS} != "0" ]
then

echo "$(date) submitting benchmarks"
'''

      for l_build in i_builds:
        for l_benchmark in i_benchmarks['runs']:
          l_slurmArguments  =  " -o "          + '${WORKING_DIRECTORY}/logs/' + l_build['id'] + '_' + l_benchmark['id'] + '.out'+\
                               " -J "          + l_build['id'] + '_' + l_benchmark['id']+\
                               " --partition=" + l_benchmark['queue']+\
                               " --nodes="    + l_benchmark['number_of_nodes']+\
                               " --time="      + l_benchmark['maximum_runtime']+\
                               " --ntasks-per-node=" + l_benchmark['ranks_per_node']

          l_submitArguments =  " -d " + l_build['compile_mode']+\
                               " -c " + l_build['code_version']+\
                               " -a " + l_build['architecture']+\
                               " -p " + l_build['parallelization']+\
                               " -u " + l_build['communication_thread']+\
                               " -s " + l_build['scalasca']+\
                               " -q " + "9"+\
                               " -o " + l_build['order']+\
                               " -n " + "${WORKING_DIRECTORY}/runs/" + l_build['id'] + '/' + l_benchmark['id'] +\
                               " -m " + l_benchmark['number_of_mpi_ranks']+\
                               " -r " + l_benchmark['ranks_per_node']+\
                               " -t " + l_benchmark['threads_per_rank']+\
                               " -b " + l_benchmark['thread_binding']+\
                               " -e " + l_build['id']+\
                               " -f " + " ${SCRIPTS_DIRECTORY}/run/"+\
                               " -i " + "${INPUT_DIRECTORY}/env_vars.sh"
          l_bashCommands = l_bashCommands + "sbatch" + l_slurmArguments + " ${SCRIPTS_DIRECTORY}/run/submit_benchmark.slurm" + l_submitArguments + "\n"

      l_bashCommands = l_bashCommands+\
                       "\nfi\n"+\
                       "# finished submitting benchmarks\n"
      return l_bashCommands

  def submitLoadLeveler( self,
                         i_type,
                         i_benchmarks,
                         i_builds ):
    self.m_logger.info( "submitting benchmarks" )

    if( i_type == 'bash' ): 
      l_bashCommands ='''
#
# Submitting benchmarks
#
if [ ${BENCHMARKS} != "0" ]
then

echo "$(date) submitting benchmarks"
'''
      for l_build in i_builds:
        for l_benchmark in i_benchmarks['runs']:
          l_benchmarkDir = '${WORKING_DIRECTORY}/runs/' + l_build['id'] + '/' + l_benchmark['id'] + '/'
          
          arguments =         '--output '+'${WORKING_DIRECTORY}/logs/' + l_build['id'] + '_' + l_benchmark['id'] + '.out'
          arguments=arguments+' --queue '+l_benchmark['queue']
          arguments=arguments+' --job_name '+l_build['id'] + '_' + l_benchmark['id']
          arguments=arguments+' --numberofnodes '+l_benchmark['number_of_nodes']
          arguments=arguments+' --time '+l_benchmark['maximum_runtime']
          arguments=arguments+' --compile_mode '+ l_build['compile_mode']
          arguments=arguments+' --code_version '+ l_build['code_version']
          arguments=arguments+' --architecture '+l_build['architecture']
          arguments=arguments+' --parallelization '+l_build['parallelization']
          arguments=arguments+' --communicationThread '+l_build['communication_thread']
          arguments=arguments+' --scalasca '+l_build['scalasca']
          arguments=arguments+' --number_of_quantities '+'9'
          arguments=arguments+' --order '+l_build['order']
          arguments=arguments+' --benchmark_directory '+l_benchmarkDir
          arguments=arguments+' --number_of_mpi_ranks '+l_benchmark['number_of_mpi_ranks']
          arguments=arguments+' --ranks_per_node '+l_benchmark['ranks_per_node']
          arguments=arguments+' --threads_per_rank '+l_benchmark['threads_per_rank']
          arguments=arguments+' --executable '+l_build['id']
          arguments=arguments+' --scripts_directory '+' ${SCRIPTS_DIRECTORY}/run/'
          arguments=arguments+' --filename '+l_benchmarkDir+'submit_benchmark.ll'
          arguments=arguments+' --genll '+"${SCRIPTS_DIRECTORY}/run/submit_benchmark.ll"


          l_bashCommands = l_bashCommands + "python ${SCRIPTS_DIRECTORY}/run/generateLoadLevelerScript.py " + arguments + "\n"
          
          l_bashCommands = l_bashCommands + "llsubmit " + l_benchmarkDir + "submit_benchmark.ll" +"\n"

      l_bashCommands = l_bashCommands+\
         "\nfi\n"+\
         "# finished submitting benchmarks\n"

      return l_bashCommands
