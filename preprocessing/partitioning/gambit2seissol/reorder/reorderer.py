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

from mesh import ReorderedMesh
from partition import ReorderedPartition

import copy
import os
import subprocess
import threading

class Reorderer:
    """Reorders a partitioned mesh using Zoltan"""
    
    def __init__(self, mesh, partition):
        # First reorder the mesh, so partitions are no longer distributed
        self.__partitionStart = []
        sum = 0
        for i in range(len(partition)):
            self.__partitionStart.append(sum)
            sum += partition.size(i)
        
        nextPos = copy.copy(self.__partitionStart)
        old2new = []
        for p in partition:
            old2new.append(nextPos[p])
            nextPos[p] += 1
            
        self.__mesh = ReorderedMesh(mesh, old2new)
        self.__partition = ReorderedPartition(partition, old2new)
        
        self.__old2new = [None] * len(self.__mesh.elements())
        
        # Run Zoltan for every partition
        for p in range(len(partition)):
            self.__zoltanError = ''
            try:
                zoltan = subprocess.Popen(['mpiexec', '-n', '1',
                        os.path.join(os.path.dirname(__file__), 'zoltan', 'zoltan'),
                        str(partition.size(p))],
                    stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            except OSError:
                raise Exception('Could not run mpiexec. Please provide correct $PATH')
        
            threading.Thread(target=self.__feedElements, args=(zoltan.stdin, p)).start()
            outThread = threading.Thread(target=self.__readZoltanOrder, args=(zoltan.stdout, p))
            outThread.start()
            errorThread = threading.Thread(target=self.__readZoltanError, args=(zoltan.stderr,))
            errorThread.start()
        
            if zoltan.wait():
                errorThread.join()
                raise Exception(self.__zoltanError.strip())
        
            outThread.join()
            errorThread.join()
        
        # Create the reordered mesh
        self.__mesh = ReorderedMesh(self.__mesh, self.__old2new)
        
    def __feedElements(self, stdin, partition):
        """Writes the element coordinates to stdin of the zoltan process"""
        
        for element in self.__mesh.elements()[self.__partitionStart[partition]:
                self.__partitionStart[partition]+self.__partition.size(partition)]:
            coords = [sum(c) / 4 for c in zip(self.__mesh.coords()[element[0]],
                self.__mesh.coords()[element[1]], self.__mesh.coords()[element[2]],
                self.__mesh.coords()[element[3]])]
            try:
                print >> stdin, coords[0], coords[1], coords[2]
            except IOError:
                # Something went wrong, zoltan should report an error to stderr
                break
            
    def __readZoltanOrder(self, stdout, partition):
        """Reads the Zoltan stdout stream and creates the ordering"""
        
        for line in stdout:
            try:
                new, old = map(int, line.split())
                self.__old2new[old+self.__partitionStart[partition]] = new+self.__partitionStart[partition]
            except Exception:
                # Something went wrong, we check for inconsistent values later
                pass
            
    def __readZoltanError(self, stderr):
        """Reads the Zoltan error stream"""
        
        for line in stderr:
            self.__zoltanError += line + os.linesep
        
    def mesh(self):
        return self.__mesh
    
    def partition(self):
        return self.__partition
        
