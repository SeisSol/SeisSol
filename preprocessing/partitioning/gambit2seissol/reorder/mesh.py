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

import copy

class ReorderedMesh:
    """A mesh that reorders elements, groups and boundaries"""
    
    def __init__(self, mesh, old2new):
        self.version = mesh.version
        self.name = mesh.name
        self.date = mesh.date
        self.problemSize = mesh.problemSize
        
        self.__coords = mesh.coords()
        
        self.__elements = [None] * len(mesh.elements())
        for i, element in enumerate(mesh.elements()):
            self.__elements[old2new[i]] = element
            
        self.__groups = []
        for group in mesh.groups():
            newGroup = group.copy()
            newGroup['cells'] = []
            for cell in group['cells']:
                newGroup['cells'].append(old2new[cell])
                
            self.__groups.append(newGroup)
                
        self.__boundaries = []
        for boundary in mesh.boundaries():
            newBoundary = copy.deepcopy(boundary)
            for i in range(len(newBoundary['sides'])):
                newBoundary['sides'][i][0] = old2new[newBoundary['sides'][i][0]]
                
            self.__boundaries.append(newBoundary)
            
    def coords(self):
        return self.__coords
    
    def elements(self):
        return self.__elements
    
    def groups(self):
        return self.__groups
    
    def boundaries(self):
        return self.__boundaries
