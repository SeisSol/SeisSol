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
# Workflow configuration.
#

import itertools
import xmltodict
import logging
from collections import OrderedDict

## Configuration of the workflow.
#
class Configuration:
  
  # member variables
  m_workflowManagement                = dict()
  m_builds                            = dict()
  m_benchmarks                        = dict()
  m_license = '''##
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
# Generated workflow.'''

  def cleanKeys( i_dictionary ):
    l_tmpDict = i_dictionary
    for k in l_tmpDict:
            if hasattr(i_dictionary[k], '__getitem__'):
                    change_keys(i_dictionary[k])
            if '.' in k:
                    i_dictionary[k.replace('.', '$')] = i_dictionary[k]
                    del i_dictionary[k]

  ## Converts an dictionary to a list
  #
  # @param i_dictionary dictionary which is converted
  def convertToList( self, i_dictionary ):
    l_list = list()

    for l_key in i_dictionary.keys():
      l_list.append(i_dictionary[l_key])

    return l_list

  ## Removes the attribute flag "@" in the keys of a dictionary and all nested dictionaries.
  #
  # @param i_dictionary dictionary which attribubtes are removed.
  def cleanKeys( self, io_dictionary ):
   for l_key in io_dictionary.keys():
     if isinstance(io_dictionary[l_key], dict):
       self.cleanKeys(io_dictionary[l_key])

     io_dictionary[l_key.replace("@", "")] = io_dictionary.pop(l_key)

  ## Splits the white space separated lists of a dictionary.
  #
  # @param i_dictionary dictionary which lists are split
  def splitLists( self, i_dictionary ):
    for l_key in i_dictionary.keys():
      i_dictionary[l_key] = i_dictionary[l_key].split()

  ## generates the cartesian product of the dictionaries including the keys 
  # @param i_dictionary dictionary which cartesian product is build
  def cartesianProduct( self, i_dictionary):
    return [dict(zip(i_dictionary, v)) for v in itertools.product(*i_dictionary.values())]

  ## converts benchmark runs to a list with one benchmark run per entry.
  #
  # @param i_dictionary benchmark runs as dictionary
  def simplifyBenchmarkRuns( self, i_dictionary ):
    l_list = list()

    for l_key in i_dictionary:
      if isinstance(i_dictionary[l_key], list):
        for l_subEntry in range( len(i_dictionary[l_key]) ):
          l_list.append( dict( i_dictionary[l_key][l_subEntry] ) )
          self.cleanKeys( l_list[-1] )
          l_list[-1]['name'] = l_key
      else:
        l_list.append( dict( i_dictionary[l_key] ) ) 
        l_list[-1]['name'] = l_key
    return l_list

  ## Constructor of the configuration
  #
  # @param i_pathToWorkflowXml configuration of the workflow and the underlying management system.
  # @param i_pathToBuildXml configuration of the code versions.
  # @param i_pathToBenchmarkXml path of the benchmark configurations.
  #
  def __init__( self,
                i_workflowDirectory,
                i_workflow ):
    self.m_logger = logging.getLogger('Configuration')

    # parse xml files
    self.m_workflowManagement  = i_workflow
    self.cleanKeys(  self.m_workflowManagement )

    self.m_builds              = xmltodict.parse( open( i_workflowDirectory+'/'+self.m_workflowManagement['builds'] )    )['builds']
    self.m_benchmarks['runs']  = xmltodict.parse( open( i_workflowDirectory+'/'+self.m_workflowManagement['benchmarks'] ))['benchmarks']

    self.cleanKeys(  self.m_builds )
    self.splitLists( self.m_builds )
    self.cleanKeys(  self.m_benchmarks['runs'] )
    self.m_benchmarks['setups'] = self.m_benchmarks['runs'].keys()

    # convert benchmark settings to simple list
    self.m_benchmarks['runs'] = self.simplifyBenchmarkRuns( self.m_benchmarks['runs'] )

    # asssemble unique benchmark run ids
    for l_i in range( len(self.m_benchmarks['runs']) ):
      self.m_benchmarks['runs'][l_i]['id'] = self.m_benchmarks['runs'][l_i]['name']+\
                                       "_" + self.m_benchmarks['runs'][l_i]['mesh_base_name']+\
                                       "_" + self.m_benchmarks['runs'][l_i]['clustered_lts']+\
                                       "_" + self.m_benchmarks['runs'][l_i]['queue']+\
                                       "_" + self.m_benchmarks['runs'][l_i]['number_of_nodes']+\
                                       "_" + self.m_benchmarks['runs'][l_i]['number_of_mpi_ranks']+\
                                       "_" + self.m_benchmarks['runs'][l_i]['ranks_per_node']+\
                                       "_" + self.m_benchmarks['runs'][l_i]['threads_per_rank']

    # build cartesian products, which represent the true builds and runs
    self.m_builds = self.cartesianProduct( self.m_builds )

    # assemble unqiue build ids
    for l_i in range(len(self.m_builds)):
     self.m_builds[l_i]['id'] = self.m_builds[l_i]['compile_mode']+\
                          "_" + self.m_builds[l_i]['code_version']+\
                          "_" + self.m_builds[l_i]['architecture']+\
                          "_" + self.m_builds[l_i]['parallelization']+\
                          "_" + self.m_builds[l_i]['communication_thread']+\
                          "_" + self.m_builds[l_i]['scalasca']+\
                          "_" + self.m_builds[l_i]['mesh_format']+\
                          "_" + self.m_builds[l_i]['output_format']+\
                          "_" + self.m_builds[l_i]['number_of_quantities']+\
                          "_" + self.m_builds[l_i]['order']+\
                          "_" + self.m_builds[l_i]['number_of_temporal_integration_points']

    self.m_logger.info( "Workflow configuration complete." )
    self.m_logger.info( str(len(self.m_builds))     + " build configurations." )
    self.m_logger.info( str(len(self.m_benchmarks['runs'])) + " benchmarks runs." )
    self.m_logger.info( str(len(self.m_builds) * len(self.m_benchmarks['runs'])) + " run configurations." )

    # assemble regular expressions for analysis of log files
    self.m_workflowManagement['log_reg_exps'] = set()

    if( 'mesh_regular_expression' in self.m_workflowManagement ):
      for l_build in self.m_builds:
        for l_benchmark in self.m_benchmarks['runs']:
          self.m_workflowManagement['log_reg_exps'].add( l_build['compile_mode']+\
                                                   "_" + l_build['code_version']+\
                                                   "_" + 'PRECISION_TAG'+\
                                                         l_build['architecture'][1:]+\
                                                   "_" + l_build['parallelization']+\
                                                   "_" + l_build['communication_thread']+\
                                                   "_" + l_build['scalasca']+\
                                                   "_" + l_build['mesh_format']+\
                                                   "_" + l_build['output_format']+\
                                                   "_" + l_build['number_of_quantities']+\
                                                   "_" + 'ORDER_TAG'+\
                                                   "_" + l_build['number_of_temporal_integration_points']+\
                                                   "_" + l_benchmark['name']+\
                                                   "_" + self.m_workflowManagement['mesh_regular_expression']+\
                                                   "_" + l_benchmark['queue']+\
                                                   "_" + l_benchmark['number_of_nodes']+\
                                                   "_" + l_benchmark['number_of_mpi_ranks']+\
                                                   "_" + l_benchmark['ranks_per_node']+\
                                                   "_" + l_benchmark['threads_per_rank'] )

    # debug information
    self.m_logger.debug( "workflow management:" )
    self.m_logger.debug( self.m_workflowManagement )

    self.m_logger.debug( "builds:" )
    self.m_logger.debug( self.m_builds )

    self.m_logger.debug( "benchmarks:" )
    self.m_logger.debug( self.m_benchmarks )
