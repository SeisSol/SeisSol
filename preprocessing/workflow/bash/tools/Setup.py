#!/bin/env python
##
# @file
# This file is part of SeisSol.
#
# @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
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
# Sets up the benchmark runs.
#
import logging
import tools.Workflow

## Setup phase of the benchmark runs.
class Setup:
  def __init__( self ):
    self.m_logger = logging.getLogger('Setup')

  ## Copys the data to the working directory
  #
  # @param i_type type of the workflow.
  # @param i_benchmarks benchmarks.
  # @param i_builds build configurations.
  # @param i_type workflow type.
  def copyData( self,
                i_benchmarks,
                i_type ):
    #
    # bash interface
    #
    if( i_type == 'bash' ):
      l_bashCommands = '''
# create directory for logs
mkdir $WORKING_DIRECTORY/logs

# copy and extract the source code
mkdir $WORKING_DIRECTORY/code
cd $WORKING_DIRECTORY/code
cp $INPUT_DIRECTORY/code.tar.gz .
tar -xzf code.tar.gz
rm code.tar.gz

# copy and extract the benchmarks
mkdir $WORKING_DIRECTORY/benchmarks
cd $WORKING_DIRECTORY/benchmarks
'''
      for l_benchmark in i_benchmarks['setups']:
        l_bashCommands = l_bashCommands \
          + "cp -r $INPUT_DIRECTORY/benchmarks/" + l_benchmark + " .\n"

      l_bashCommands = l_bashCommands \
        + "\n# create run directories and link benchmark setups\n" \
        + "mkdir $WORKING_DIRECTORY/runs\n" \
        + "cd $WORKING_DIRECTORY/runs\n"

    return l_bashCommands

  ## Builds the code versions.
  #
  # @param i_type type of the workflow.
  # @param i_builds build configurations.
  def buildCodeVersions( self,
                         i_type,
                         i_builds ):
    l_bashCommands = "\n# build the code versions\n" \
                   + "cd $WORKING_DIRECTORY/code\n"
    for l_build in i_builds:
      l_buildArguments =  " -e " + "${INPUT_DIRECTORY}/env_vars.sh"+\
                          " -d " + l_build['compile_mode']+\
                          " -c " + l_build['code_version']+\
                          " -a " + l_build['architecture']+\
                          " -p " + l_build['parallelization']+\
                          " -r " + l_build['communication_thread']+\
                          " -s " + l_build['scalasca']+\
                          " -m " + l_build['mesh_format']+\
                          " -u " + l_build['output_format']+\
                          " -q " + l_build['number_of_quantities']+\
                          " -o " + l_build['order']+\
                          " -t " + l_build['number_of_temporal_integration_points']+\
                          " -n " + l_build['id']
      l_bashCommands = l_bashCommands + "sh $SCRIPTS_DIRECTORY/setup/build.sh"+\
                       l_buildArguments+\
                       " > ${WORKING_DIRECTORY}/logs/build_" +  l_build['id'] +".log &\n"

    return l_bashCommands+"wait\n"

  ## Set up the run directories.
  #
  # @param i_type type of the workflow.
  # @param i_benchmarks benchmarks.
  # @param i_builds build configurations.
  def setupRunDirectories( self,
                           i_type,
                           i_benchmarks,
                           i_builds ):
    if( i_type == 'bash' ): 
      l_bashCommands = "\n# create run directories\n"
      for l_build in i_builds:
        l_bashCommands = l_bashCommands \
          + "mkdir $WORKING_DIRECTORY/runs/" + l_build['id'] + "\n" \
          + "cd $WORKING_DIRECTORY/runs/" + l_build['id'] + "\n"

        for l_benchmark in i_benchmarks['runs']:
          l_bashCommands = l_bashCommands \
            + "mkdir " + l_benchmark['id'] + "\n" \
            + "mkdir " + l_benchmark['id'] + "/output\n"\
            + "ln -s ${WORKING_DIRECTORY}/benchmarks/" + l_benchmark['name'] + "/setup/* " + l_benchmark['id'] + "\n"\
            + "ln -s ${WORKING_DIRECTORY}/code/preprocessing/matrices/*" + " " + l_benchmark['id'] + "\n"\
            + "echo  ${WORKING_DIRECTORY}/code/Maple/ > " + l_benchmark['id'] + '/' + 'DGPATH' + "\n"\
            + "ln -s ${WORKING_DIRECTORY}/code/build/" + l_build['id'] + " " + l_benchmark['id'] + '/' + l_build['id'] + "\n"\
            + "ln -s ${WORKING_DIRECTORY}/code/Maple " + l_benchmark['id'] + "\n"

          l_parameterFileArguments = " -e " + "${INPUT_DIRECTORY}/env_vars.sh"+\
                                     " -p " + "cpp"+\
                                     " -t " + l_benchmark['id'] + "/parameters.template"+\
                                     " -o " + l_benchmark['id'] + "/parameters.par"+\
                                     " -r " + l_build['order']+\
                                     " -f " + l_build['mesh_format']+\
                                     " -q " + ("6" if (l_build['output_format'] == 'hdf5') else "10")+\
                                     " -l " + l_benchmark['clustered_lts']+\
                                     " -b " + l_benchmark['mesh_base_name']

          if( l_build['mesh_format'] == 'Netcdf' ):
            l_parameterFileArguments = l_parameterFileArguments + '_' + l_benchmark['number_of_mpi_ranks']

          l_bashCommands = l_bashCommands +  "sh ${SCRIPTS_DIRECTORY}/setup/generate_parameter_file.sh" + l_parameterFileArguments + "\n"

      return l_bashCommands
  

  ## Prepares the working directory.
  #
  # @param i_type type of the workflow.
  # @param i_benchmarks benchmarks.
  # @param i_builds build configurations.
  def prepareWorkingDirectory( self,
                               i_type,
                               i_benchmarks,
                               i_builds ):
    self.m_logger.info( "preparing working directory" )

    #
    # bash interface
    #
    if( i_type == 'bash' ): 
      l_bashCommands ='''
#
# Preparing benchmarks
#
if [ ${SETUP} != "0" ]
then

echo "$(date) peparing the working directory"

mkdir $WORKING_DIRECTORY
'''

      l_bashCommands = l_bashCommands+\
                       "echo \"$(date) copying data\"\n"+\
                       self.copyData( i_benchmarks = i_benchmarks,
                                      i_type       = i_type )

      l_bashCommands = l_bashCommands+\
                       "echo \"$(date) building code versions\"\n"+\
                       self.buildCodeVersions( i_type   = i_type,
                                               i_builds = i_builds )

      l_bashCommands = l_bashCommands+\
                       "echo \"$(date) setting up run directories\"\n"+\
                       self.setupRunDirectories( i_type       = i_type,
                                                 i_builds     = i_builds,
                                                 i_benchmarks = i_benchmarks )

      l_bashCommands = l_bashCommands+\
                       "\nfi\n"+\
                       "# finished preparing benchmarks\n"

      return l_bashCommands
