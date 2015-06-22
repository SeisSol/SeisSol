#!/usr/bin/python
##
# @file
# This file is part of SeisSol.
#
# @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
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

import metis

import subprocess

METIS_MESH = 'metis.mesh'
METIS_GRAPH = 'metis.graph'

class Partitioner:
    """Converts a mesh into graph and partitions it using metis"""
    
    def __init__(self, mesh, partitions, tmpdir):
        metisMesh = tmpdir.path(METIS_MESH)
        
        # Write metis mesh
        metis.MeshWriter(metisMesh, mesh.elements())
        
        # Convert to graph
        metisGraph = tmpdir.path(METIS_GRAPH)
        p = subprocess.Popen(['m2gmetis', '-ncommon=3', metisMesh, metisGraph],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        _, errmsg = p.communicate()
        if p.returncode:
            raise Exception(errmsg.strip())
        
        # Run metis
        p = subprocess.Popen(['gpmetis', '-ptype=rb', metisGraph, str(partitions)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        _, errmsg = p.communicate()
        if p.returncode:
            raise Exception(errmsg.strip())
        
        # Read partitions
        self.__partition = metis.PartitionReader(metisGraph+'.part.'+str(partitions),
            partitions, len(mesh.elements()))
        
        if self.__partition.size() != len(mesh.elements()):
            raise Exception('Mesh size and partition size do not match: mesh size = '
                +str(len(mesh.elements()))+' != partition size = '+str(self.__partition.size()))
            
    def partition(self):
        return self.__partition
