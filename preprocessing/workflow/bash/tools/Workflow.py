#!/bin/env python
##
# @file
# This file is part of SeisSol.
#
# @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
#
# @author Leonhard Rannabauer (lrannabauer AT mytum.de)
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
# Generates the workflow itself.
#

import logging
import sys
import datetime

import tools.Setup
import tools.Run
import tools.Analyze

class Workflow:
  m_management = dict()
  m_setup   = tools.Setup.Setup()
  m_run     = tools.Run.Run()
  m_analyze = tools.Analyze.Analyze()

  ## Initializes the workflow.
  #
  # @param i_management underlying workflow management.
  def __init__( self, i_management ):
    self.m_logger = logging.getLogger('Workflow')

    self.m_management = i_management

    # clean output location
    if self.m_management['type'] == 'bash':
      open( self.m_management['output_location'], 'w')

  ## Writes the string to the output location.
  #
  # @param i_string string to write
  def write( self, i_string ):
    #
    # bash interface
    #
    if( self.m_management['type'] == 'bash' ):
      if( self.m_management['output_location'] == 'stdout' ):
        l_location = sys.stdout
      else:
        l_location = open( self.m_management['output_location'], 'a' )
      

      l_location.write( i_string )

  ## Writes the line to the output locations.
  #
  # @param i_string string to write.
  def writeLine( self, i_string ):
    self.write( i_string+"\n" )

  ## Generates the file header.
  #
  # @param i_license license of the file header.
  def generateHeader( self, i_license ):
    self.m_logger.info( 'generating workflow header' )
    #
    # bash interface
    #
    if( self.m_management['type'] == 'bash' ):
      self.m_logger.debug( 'bash header' )
      self.writeLine( "#!/bin/sh" )
      self.write( i_license )
      self.writeLine( " -- time of generation: " + str(datetime.datetime.now()) )

  ## Generates the user interface to the workflows
  def generateUserInterface( self ):
    self.m_logger.info( 'generating workflow interface' )

    #
    # bash interface
    #
    if( self.m_management['type'] == 'bash' ):
      self.m_logger.debug( 'bash interface' )

      l_userInterface='''
#
# User interface
#
show_help() {
cat << EOF
Usage: ${0##*/} [-h -s -b -a][-i INPUT_DIRECTORY -e SCRIPTS_DIRECTORY -w WORKING_DIRECTORY]
Workflow management. 
     -h display this help and exit
     -s runs the setup phase of the workflow
     -b runs the benchmarks of the workflow.
     -a runs the analyzation phase of the workflow.
     -i INPUT_DIRECTORY
     -e SCRIPTS_DIRECTORY
     -w WORKING_DIRECTORY
EOF
}

# parse command line arguments
SETUP=0
BENCHMARKS=0
ANALYZE=0
INPUT_DIRECTORY=NOT_SET
WORKING_DIRECTORY=NOT_SET

OPTIND=1
while getopts "hsbai:e:w:" opt; do
    case "$opt" in
        h)
            show_help
            exit 0
            ;;
        s) SETUP=1
            ;;
        b) BENCHMARKS=1
            ;;
        a) ANALYZE=1
            ;;
        i) if [ $(echo $OPTARG | head -c 1) != "/" ]
           then
             # extend relative path with current path
             INPUT_DIRECTORY="$(pwd)/${OPTARG}"
           else
             INPUT_DIRECTORY=$OPTARG
           fi
            ;;
        e) if [ $(echo $OPTARG | head -c 1) != "/" ]
           then
             # extend relative path with current path
             SCRIPTS_DIRECTORY="$(pwd)/${OPTARG}"
           else
             SCRIPTS_DIRECTORY=$OPTARG
           fi
            ;;
        w) if [ $(echo $OPTARG | head -c 1) != "/" ]
           then
             # extend relative path with current path
             WORKING_DIRECTORY="$(pwd)/${OPTARG}"
           else
             WORKING_DIRECTORY=$OPTARG
           fi
            ;;
        '?')
            show_help >&2
            exit 1
            ;;
    esac
done
shift "$((OPTIND-1))" # Shift off the options and optional --.

# setting up environment
. ${INPUT_DIRECTORY}/env_vars.sh'''

      if self.m_management['resource_management'] == 'loadleveler':
         l_userInterface = l_userInterface + '''
rm $HOME/tmp_env.sh
touch $HOME/tmp_env.sh
echo 'export WORKING_DIRECTORY='${WORKING_DIRECTORY} >> $HOME/tmp_env.sh
echo 'export INPUT_DIRECTORY='${INPUT_DIRECTORY} >> $HOME/tmp_env.sh
echo 'export SCRIPTS_DIRECTORY='${SCRIPTS_DIRECTORY} >> $HOME/tmp_env.sh'''

      self.writeLine( l_userInterface )

  ## Generates the directory structure
  def prepareWorkingDirectory( self, i_builds, i_benchmarks ):
    self.write( self.m_setup.prepareWorkingDirectory( i_type        = self.m_management['type'],
                                                      i_builds      = i_builds,
                                                      i_benchmarks  = i_benchmarks ) )

  ## Submits the runs
  def submit( self, i_builds, i_benchmarks ):
    self.write( self.m_run.submit( i_type        = self.m_management['type'],
                                   i_resource_management = self.m_management['resource_management'],
                                   i_builds      = i_builds,
                                   i_benchmarks  = i_benchmarks ) )

  def postProcess( self, i_builds, i_benchmarks ):
    self.write( self.m_analyze.postProcess( i_type               = self.m_management['type'],
                                            i_builds             = i_builds,
                                            i_benchmarks         = i_benchmarks,
                                            i_regularExpressions = self.m_management['log_reg_exps'] ) )
