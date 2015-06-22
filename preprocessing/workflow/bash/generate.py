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
# Interface to the workflow generation.
#

import argparse
import tools.Configuration
import tools.Workflow
import logging
import xmltodict
import os

# set up logging level
logging.basicConfig( level=logging.INFO,
                     format='%(asctime)s - %(name)s - %(levelname)s - %(message)s' )

# command line interface
logging.info( "parsing command line arguments" )

l_parser    = argparse.ArgumentParser( description='Generates the workflows based on the given XML-configuration.' )
l_parser.add_argument( '--workflows_xml',
                       dest     = "workflows_xml",
                       required = True,
                       help     = "Configurations of the workflows.",
                       metavar  = "WORKFLOWS_XML" )

l_arguments = vars(l_parser.parse_args())

# generate the workflows
l_workflows = xmltodict.parse( open( l_arguments['workflows_xml'] ) )

for l_workflow in l_workflows['workflows']:
  l_configuration = tools.Configuration.Configuration( i_workflowDirectory = os.path.dirname(l_arguments['workflows_xml']),
                                                       i_workflow = l_workflows['workflows'][l_workflow] )

  l_workflow = tools.Workflow.Workflow( i_management = l_configuration.m_workflowManagement )
  l_workflow.generateHeader( l_configuration.m_license )
  l_workflow.generateUserInterface()

  l_workflow.prepareWorkingDirectory( i_builds =      l_configuration.m_builds,
                                      i_benchmarks  = l_configuration.m_benchmarks )

  l_workflow.submit( i_builds =      l_configuration.m_builds,
                     i_benchmarks  = l_configuration.m_benchmarks )

  l_workflow.postProcess( i_builds            = l_configuration.m_builds,
                          i_benchmarks        = l_configuration.m_benchmarks )
