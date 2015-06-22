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

import collections

class MeshWriter:
    """Writes a mesh in metis format"""
    
    def __init__(self, file, elements):
        if isinstance(file, basestring):
            file = open(file, 'w')
        
        # Write the elements
        print >> file, len(elements)
        for element in elements:
            print >> file, ' '.join(map(lambda x: str(x+1), element))
            
        file.close()
        
class PartitionReader(collections.Iterable):
    """Reads a metis partition file"""
    
    def __init__(self, file, length, size):
        if isinstance(file, basestring):
            file = open(file, 'rU')
            
        self.__file = file
        
        # Length = number of partitions
        self.__len = length
        
        # Size0 = number of elements
        self.__size0 = size
        
        self.__partition = None
    
    def __del__(self):
        self.__file.close()
        
    def __len__(self):
        return self.__len
    
    def __readFile(self):
        if self.__partition:
            # File was already parsed
            return
        
        # Read the partition for each element
        self.__file.seek(0)
        self.__partition = [int(line) for line in self.__file]
        
        # Get the size of each partition
        self.__size = [0] * self.__len
        for p in self.__partition:
            self.__size[p] += 1
    
    def size(self, partition = -1):
        if partition < 0:
            return self.__size0
        
        self.__readFile()
        return self.__size[partition]
    
    def __getitem__(self, key):
        self.__readFile()
        return self.__partition[key]
    
    def __iter__(self, default=None):
        if self.__partition:
            for p in self.__partition:
                yield p
        else:
            for line in self.__file:
                yield int(line)
            
