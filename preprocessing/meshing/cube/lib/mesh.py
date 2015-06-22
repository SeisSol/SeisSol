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
import datetime
import numpy

class Coords(collections.Sequence):
    def __init__(self, mesh):
        self.__mesh = mesh
        self.__len1 = [2*i + 1 for i in self.__mesh.size()] # number of coords in each dimension
        self.__len = reduce(lambda x, y: x*y, self.__len1)
        
    def len1(self):
        return self.__len1
        
    def __len__(self):
        return self.__len
    
    def __getitem__(self, index):
        if index >= self.__len:
            raise StopIteration
        
        def scale(x, dim):
            return (x - (self.__len1[dim]-1)/2.) * 5. / self.__mesh.size()[dim]
        
        x = scale(index % self.__len1[0], 0)
        y = scale((index / self.__len1[0]) % self.__len1[1], 1) 
        z = scale(index / (self.__len1[0] * self.__len1[1]), 2)
        
        return (x, y, z)
    
class Elements(collections.Sequence):
    def __init__(self, mesh):
        self.__mesh = mesh
        self.__len = 40 * reduce(lambda x, y: x*y, self.__mesh.size())
        self.__len1 = [2 * i for i in self.__mesh.size()] # number of cubes in one dimension
        
        self.__checkOrientation()
        
    def __checkOrientation(self):
        """Checks the correct orientation for the different tetrahedra"""
        
        for i in range(10):
            vertices = self[i]
            
            coords = [None]*4
            for j in range(4):
                coords[j] = numpy.array(self.__mesh.coords()[vertices[j]])
                
            n = numpy.cross(-coords[0]+coords[1], -coords[0]+coords[2])
            orientation = numpy.dot(n, -coords[0]+coords[3])
            
            if orientation <= 0:
                raise Exception('Wrong orientation in element '+str(i))
        
    def __len__(self):
        return self.__len
    
    def tetUnitCube(self, index, even):
        """Returns the coordinates of one of the 5 tetrahedron in an even or odd unit cube"""
        
        if even:
            if index == 0:
                return ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
            if index == 1:
                return ((1, 0, 0), (0, 1, 0), (1, 1, 1), (1, 1, 0))
            if index == 2:
                return ((1, 0, 0), (1, 1, 1), (0, 0, 1), (1, 0, 1))
            if index == 3:
                return ((0, 1, 0), (0, 1, 1), (0, 0, 1), (1, 1, 1))
            # index == 4
            return ((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1))
    
        if index == 0:
            return ((0, 0, 0), (0, 1, 0), (0, 1, 1), (1, 1, 0))
        if index == 1:
            return ((0, 0, 0), (1, 1, 0), (1, 0, 1), (1, 0, 0))
        if index == 2:
            return ((0, 0, 0), (1, 0, 1), (0, 1, 1), (0, 0, 1))
        if index == 3:
            return ((1, 1, 0), (1, 0, 1), (1, 1, 1), (0, 1, 1))
        # index == 4
        return ((0, 0, 0), (1, 1, 0), (0, 1, 1), (1, 0, 1))
    
    def __getitem__(self, index):
        if index >= self.__len:
            raise StopIteration
        
        # Index of the cube we are currently working on
        cIndex = index/5
        
        cx = cIndex % self.__len1[0]
        cy = (cIndex / self.__len1[0]) % self.__len1[1]
        cz = cIndex / (self.__len1[0] * self.__len1[1])
        
        # The index inside the cube 
        i = index % 5
        
        # Odd cube?
        odd = (cx+cy+cz) % 2
        
        # Number of coords in one dimension
        coordLength = self.__mesh.coords().len1()
        
        def crd((x, y, z)):
            return x+cx + ((y+cy) + (z+cz) * coordLength[1]) * coordLength[0]
        
        return map(crd, self.tetUnitCube(i, not odd))

class Group:
    def __init__(self, name, material, mesh):
        self._name = name
        self._size = len(mesh.elements())
        self._material = material
        self._cells = xrange(self._size)
        
    def __getitem__(self, key):
        return getattr(self, '_'+str(key))
        
class Boundary:
    class Sides(collections.Sequence):
        """Generates the boundary condition for all face elements"""
        
        def __init__(self, mesh):
            self.__len1 = [8 * i for i in  # Number of faces on one side
                [mesh.size()[0]*mesh.size()[1], # Top / Bottom
                 mesh.size()[0]*mesh.size()[2], # Left / Right behind
                 mesh.size()[1]*mesh.size()[2]]] # Right / Left behind
            self.__len = 2 * reduce(lambda x, y: x+y, self.__len1)
            self.__mesh = mesh
            
            def face(points):
                """Returns the number of the face or None if its not a face"""
        
                if points == [0, 0, 0, 1]:
                    return 1
                if points == [0, 0, 1, 0]:
                    return 2
                if points == [0, 1, 0, 0]:
                    return 4
                if points == [1, 0, 0, 0]:
                    return 3
                
                return None
        
            self.__top = []
            self.__left = []
            self.__right = []
            self.__bot = []
            self.__rightb = []
            self.__leftb = []
            
            for i in range(5):
                element1 = mesh.elements().tetUnitCube(i, 1)
                
                f = face(map(lambda e: e[0], element1))
                if f:
                    self.__right.append((i, f))
                f = face(map(lambda e: e[1], element1))
                if f:
                    self.__left.append((i, f))
                f = face(map(lambda e: e[2], element1))
                if f:
                    self.__top.append((i, f))
                
            for i in range(5):
                element0 = mesh.elements().tetUnitCube(i, 0)
                
                f = face(map(lambda e: e[0], element0))
                if f:
                    self.__right.append((i, f))
                f = face(map(lambda e: e[1], element0))
                if f:
                    self.__left.append((i, f))
                f = face(map(lambda e: e[2], element0))
                if f:
                    self.__top.append((i, f))
                    
                f = face(map(lambda e: 1-e[0], element0))
                if f:
                    self.__leftb.append((i, f))
                f = face(map(lambda e: 1-e[1], element0))
                if f:
                    self.__rightb.append((i, f))
                f = face(map(lambda e: 1-e[2], element0))
                if f:
                    self.__bot.append((i, f))
                        
            for i in range(5):
                element1 = mesh.elements().tetUnitCube(i, 1)
                
                f = face(map(lambda e: 1-e[0], element1))
                if f:
                    self.__leftb.append((i, f))
                f = face(map(lambda e: 1-e[1], element1))
                if f:
                    self.__rightb.append((i, f))
                f = face(map(lambda e: 1-e[2], element1))
                if f:
                    self.__bot.append((i, f))
                    
            assert(len(self.__top) == 4)
            assert(len(self.__left) == 4)
            assert(len(self.__right) == 4)
            assert(len(self.__bot) == 4)
            assert(len(self.__rightb) == 4)
            assert(len(self.__leftb) == 4)
            
        def __len__(self):
            return self.__len
        
        def __getitem__(self, index):
            if index >= self.__len:
                raise StopIteration
            
            len1All = self.__len1 + self.__len1
            for side in range(6):
                if index < len1All[side]:
                    break
                index -= len1All[side]
            
            def getFace(index, sizeX, sizeY):
                x = (index / 2) % sizeX
                y = (index / (2 * sizeX)) % sizeY
                
                face = index % 4
            
                if y%2 == 1:
                    # odd rows ...
                    face = (face+2)%4
                    
                return (x, y, face)
            
            size = [2 * i for i in self.__mesh.size()] # number of cubes in each dimension
            
            if side == 0: # top (x = x; y = y)
                (x, y, face) = getFace(index, size[0], size[1])
                
                offset = (x + y*size[0]) * 5
                f = self.__top[face]
            elif side == 1: # left (x = x; y = z)
                (x, y, face) = getFace(index, size[0], size[2])
                
                offset = (x + y*size[1]*size[0]) * 5
                f = self.__left[face]
            elif side == 2: # right (x = y; y = z)
                (x, y, face) = getFace(index, size[1], size[2])
                
                offset = (x + y*size[1]) * size[0] * 5
                f = self.__right[face]
            elif side == 3: # bottom (x = x; y = y)
                (x, y, face) = getFace(index, size[0], size[1])
                
                offset = (x + (y + (size[2]-1)*size[1])*size[0]) * 5
                f = self.__bot[face]
            elif side == 4: # right behind (x = x; y = z)
                (x, y, face) = getFace(index, size[0], size[2])
                
                offset = (x + (y*size[1] + size[1]-1)*size[0]) * 5
                f = self.__rightb[face]
            else: # side == 5 # left behind (x = y; y = z)
                (x, y, face) = getFace(index, size[1], size[2])
                
                offset = ((x + y*size[1]) * size[0] + size[0]-1) * 5
                f = self.__leftb[face]
        
            return (f[0]+offset, f[1])
    
    def __init__(self, name, mesh):
        self._name = name
        self._sides = Boundary.Sides(mesh)
        self._size = len(self._sides)
        
    def __getitem__(self, key):
        return getattr(self, '_'+str(key))

class Mesh:
    version = '2.0.0'
    name = 'cube'
    date = datetime.datetime.now().strftime('%b %Y')

    def __init__(self, size, boundaryCond):
        # The smallest cube for this mesh has size 2*2*2
        self.__size = [i / 2 for i in size]

        self.name += '_'.join(map(str, size))
        
        self.__coords = Coords(self)
        self.__elements = Elements(self)
        self.__groups = [Group('fluid', 2, self)]
        self.__boundaries = [Boundary(boundaryCond, self)]
        
    def size(self):
        return self.__size
        
    def coords(self):
        return self.__coords
    
    def elements(self):
        return self.__elements

    def groups(self):
        return self.__groups
    
    def boundaries(self):
        return self.__boundaries
