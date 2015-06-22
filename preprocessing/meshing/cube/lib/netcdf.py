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
#
# @section DESCRIPTION
# TODO Currently all partitions must have the same size !!!
#

import operator
import functools
from netCDF4 import Dataset

def unique(seq, idfun=None):  
    # order preserving 
    if idfun is None: 
        def idfun(x): return x 
    seen = {} 
    result = [] 
    for item in seq: 
        marker = idfun(item) 
        # in old Python versions: 
        # if seen.has_key(marker) 
        # but in new ones: 
        if marker in seen: continue 
        seen[marker] = 1 
        result.append(item) 
    return result

class NetcdfWriter:
    def __init__(self, mesh, partitions, file):
        file.close() # Close the file object because argsparse opens it automatically
        
        numPartitions = reduce(operator.mul, partitions, 1)
        partitionCubes = [(mesh.size()[i] + partitions[i] - 1) / partitions[i] * 2 for i in range(3)]
        partitionSize = reduce(operator.mul, partitionCubes, 5) # 5*8 elements per double cube
        for i in range(3):
            if mesh.size()[i] % partitions[i] != 0:
                raise Exception('Different partition sizes currently not supported')
        
        file = Dataset(file.name, 'w', format='NETCDF4')
        
        # Dimensions
        dimDimension = file.createDimension('dimension', 3)
        dimPartitions = file.createDimension('partitions', numPartitions)
        
        dimElem = file.createDimension('elements', partitionSize)
        dimElemSides = file.createDimension('element_sides', 4)
        dimElemVertices = file.createDimension('element_vertices_dim', 4) # Python bindings are to stupid so we need another name for the dimension
        
        vrtxSize = reduce(operator.mul, map(lambda i: i+1, partitionCubes), 1)
        dimVrtx = file.createDimension('vertices', vrtxSize)
        
        boundarySizes = [2*partitionCubes[1]*partitionCubes[2], 2*partitionCubes[0]*partitionCubes[2], 2*partitionCubes[0]*partitionCubes[1]]
        boundarySizes = boundarySizes[::-1] + boundarySizes
        dimBnd = file.createDimension('boundaries', 6)
        dimBndElem = file.createDimension('boundary_elements', max(boundarySizes))
        
        # Variables
        varElemSize = file.createVariable('element_size', 'i4', ('partitions',))
        varElemVertices = file.createVariable('element_vertices', 'i4', ('partitions', 'elements', 'element_vertices_dim'))
        varElemNeighbors = file.createVariable('element_neighbors', 'i4', ('partitions', 'elements', 'element_sides'))
        varElemBoundaries = file.createVariable('element_boundaries', 'i4', ('partitions', 'elements', 'element_sides'))
        varElemNeighborSides = file.createVariable('element_neighbor_sides', 'i4', ('partitions', 'elements', 'element_sides'))
        varElemSideOrientations = file.createVariable('element_side_orientations', 'i4', ('partitions', 'elements', 'element_sides'))
        varElemNeighborRanks = file.createVariable('element_neighbor_ranks', 'i4', ('partitions', 'elements', 'element_sides'))
        varElemMPIIndices = file.createVariable('element_mpi_indices', 'i4', ('partitions', 'elements', 'element_sides'))
        
        varVrtxSize = file.createVariable('vertex_size', 'i4', ('partitions',))
        varVrtxCoords = file.createVariable('vertex_coordinates', 'f8', ('partitions', 'vertices', 'dimension'))
        
        varBndSize = file.createVariable('boundary_size', 'i4', ('partitions',))
        varBndElemSize = file.createVariable('boundary_element_size', 'i4', ('partitions', 'boundaries'))
        varBndElemRank = file.createVariable('boundary_element_rank', 'i4', ('partitions', 'boundaries'))
        varBndElemLocalIds = file.createVariable('boundary_element_localids', 'i4', ('partitions', 'boundaries', 'boundary_elements'))
        
        # Elements
        varElemSize[:] = numPartitions*[partitionSize]
        
#         for z in range(partitions[2]):
#             for y in range(partitions[1]):
#                 for x in range(partitions[0]):
        vertices = []
                    
        for zz in range(partitionCubes[2]):
            for yy in range(partitionCubes[1]):
                for xx in range(partitionCubes[0]):
                    for i in range(5):
                        #vertices.append(mesh.elements()[(((z*partitionCubes[2]+zz)*mesh.size()[1]*2 + (y*partitionCubes[1]+yy))*mesh.size()[0]*2 + (x*partitionCubes[0]+xx))*5 + i])
                        vertices.append(mesh.elements()[((zz*mesh.size()[1]*2 + yy)*mesh.size()[0]*2 + xx)*5 + i])
                    
        verticesMap = {}
        for i in range(len(vertices)):
            for j in range(len(vertices[i])):
                if vertices[i][j] in verticesMap:
                    vertices[i][j] = verticesMap[vertices[i][j]]
                else:
                    n = len(verticesMap)
                    verticesMap[vertices[i][j]] = n
                    vertices[i][j] = n
        
        for i in range(numPartitions):
            varElemVertices[i,:,:] = vertices                                        
        #varElemVertices[(z*partitions[1] + y)*partitions[0] + x,:,:] = vertices
                                        
        print 'element_vertices done'
            
        neighbors = []
        for z in range(partitionCubes[2]):
            for y in range(partitionCubes[1]):
                for x in range(partitionCubes[0]):
                    odd = (x + y + z) % 2
                    
                    if odd:
                        n = [[partitionSize if x == 0 else -4,
                              partitionSize if z == 0 else -partitionCubes[1]*partitionCubes[0]*5+3,
                              4,
                              partitionSize if y == partitionCubes[1]-1 else partitionCubes[0]*5],
                             [4,
                              partitionSize if z == 0 else -partitionCubes[1]*partitionCubes[0]*5+2,
                              partitionSize if y == 0 else -partitionCubes[0]*5+1,
                              partitionSize if x == partitionCubes[0]-1 else 5],
                             [4,
                              partitionSize if y == 0 else -partitionCubes[0]*5+3,
                              partitionSize if x == 0 else -3,
                              partitionSize if z == partitionCubes[2]-1 else partitionCubes[1]*partitionCubes[0]*5],
                             [partitionSize if x == partitionCubes[0]-1 else 8,
                              4,
                              partitionSize if y == partitionCubes[1]-1 else partitionCubes[0]*5+2,
                              partitionSize if z == partitionCubes[2]-1 else partitionCubes[1]*partitionCubes[0]*5+1],
                             [0, 1, 2, 3]]
                    else:
                        n = [[partitionSize if z == 0 else -partitionCubes[1]*partitionCubes[0]*5+2, 
                              partitionSize if y == 0 else -partitionCubes[0]*5,
                              partitionSize if x == 0 else -4,
                              4],
                             [4,
                              partitionSize if z == 0 else -partitionCubes[1]*partitionCubes[0]*5+3,
                              partitionSize if x == partitionCubes[0]-1 else 5,
                              partitionSize if y == partitionCubes[1]-1 else partitionCubes[0]*5+1],
                             [4,
                              partitionSize if x == partitionCubes[0]-1 else 7,
                              partitionSize if y == 0 else -partitionCubes[0]*5+3,
                              partitionSize if z == partitionCubes[2]-1 else partitionCubes[1]*partitionCubes[0]*5+1],
                             [partitionSize if x == 0 else -2,
                              partitionSize if y == partitionCubes[1]-1 else partitionCubes[0]*5+2,
                              4,
                              partitionSize if z == partitionCubes[2]-1 else partitionCubes[1]*partitionCubes[0]*5],
                             [0, 1, 2, 3]]
                        
                    offset = ((z*partitionCubes[1] + y)*partitionCubes[0] + x)*5
                    n = map(lambda a: map(lambda i: i + offset if i < partitionSize else i, a), n)
                    neighbors.extend(n)
                    
        for i in range(numPartitions):
            varElemNeighbors[i,:,:] = neighbors
            
        print 'element_neighbors done'
            
        for z in range(partitions[2]):
            for y in range(partitions[1]):
                for x in range(partitions[0]):
                    boundaries = []
                    
                    for zz in range(partitionCubes[2]):
                        for yy in range(partitionCubes[1]):
                            for xx in range(partitionCubes[0]):
                                odd = (xx + yy + zz) % 2

				# TODO support different boundaries
                                
                                if odd:
                                    b = [[1 if x == 0 and xx == 0 else 0,
                                          1 if z == 0 and zz == 0 else 0,
                                          0,
                                          1 if y == partitions[1]-1 and yy == partitionCubes[1]-1 else 0],
                                         [0,
                                          1 if z == 0 and zz == 0 else 0,
                                          1 if y == 0 and yy == 0 else 0,
                                          1 if x == partitions[0]-1 and xx == partitionCubes[0]-1 else 0],
                                         [0,
                                          1 if y == 0 and yy == 0 else 0,
                                          1 if x == 0 and xx == 0 else 0,
                                          1 if z == partitions[2]-1 and zz == partitionCubes[2]-1 else 0],
                                         [1 if x == partitions[0]-1 and xx == partitionCubes[0]-1 else 0,
                                          0,
                                          1 if y == partitions[1]-1 and yy == partitionCubes[1]-1 else 0,
                                          1 if z == partitions[2]-1 and zz == partitionCubes[2]-1 else 0],
                                         [0, 0, 0, 0]]
                                else:
                                    b = [[1 if z == 0 and zz == 0 else 0, 
                                          1 if y == 0 and yy == 0 else 0,
                                          1 if x == 0 and xx == 0 else 0,
                                          0],
                                         [0,
                                          1 if z == 0 and zz == 0 else 0,
                                          1 if x == partitions[0]-1 and xx == partitionCubes[0]-1 else 0,
                                          1 if y == partitions[1]-1 and yy == partitionCubes[1]-1 else 0],
                                         [0,
                                          1 if x == partitions[0]-1 and xx == partitionCubes[0]-1 else 0,
                                          1 if y == 0 and yy == 0 else 0,
                                          1 if z == partitions[2]-1 and zz == partitionCubes[2]-1 else 0],
                                         [1 if x == 0 and xx == 0 else 0,
                                          1 if y == partitions[1]-1 and yy == partitionCubes[1]-1 else 0,
                                          0,
                                          1 if z == partitions[2]-1 and zz == partitionCubes[2]-1 else 0],
                                         [0, 0, 0, 0]]
                                
                                boundaries.extend(b)
                                
                    varElemBoundaries[(z*partitions[1] + y)*partitions[0] + x,:,:] = boundaries
                    
        print 'element_boundaries done'
                    
        for z in range(partitions[2]):
            for y in range(partitions[1]):
                for x in range(partitions[0]):
                    sides = []
                    
                    for zz in range(partitionCubes[2]):
                        for yy in range(partitionCubes[1]):
                            for xx in range(partitionCubes[0]):
                                odd = (xx + yy + zz) % 2
                                
                                if odd:
                                    s = [[0 if x == 0 and xx == 0 else 2,
                                          0 if z == 0 and zz == 0 else 3,
                                          0,
                                          0 if y == partitions[1]-1 and yy == partitionCubes[1]-1 else 1],
                                         [1,
                                          0 if z == 0 and zz == 0 else 3,
                                          0 if y == 0 and yy == 0 else 3,
                                          0 if x == partitions[0]-1 and xx == partitionCubes[0]-1 else 2],
                                         [2,
                                          0 if y == 0 and yy == 0 else 1,
                                          0 if x == 0 and xx == 0 else 1,
                                          0],# if z == partitions[2]-1 and zz == partitionCubes[2]-1 else 0],
                                         [0,# if x == partitions[0]-1 and xx == partitionCubes[0]-1 else 0,
                                          3,
                                          0 if y == partitions[1]-1 and yy == partitionCubes[1]-1 else 2,
                                          0 if z == partitions[2]-1 and zz == partitionCubes[2]-1 else 1],
                                         [2, 0, 0, 1]]
                                else:
                                    s = [[0 if z == 0 and zz == 0 else 3, 
                                          0 if y == 0 and yy == 0 else 3,
                                          0 if x == 0 and xx == 0 else 3,
                                          0],
                                         [1,
                                          0 if z == 0 and zz == 0 else 3,
                                          0,# if x == partitions[0]-1 and xx == partitionCubes[0]-1 else 0,
                                          0 if y == partitions[1]-1 and yy == partitionCubes[1]-1 else 2],
                                         [2,
                                          0 if x == partitions[0]-1 and xx == partitionCubes[0]-1 else 2,
                                          0 if y == 0 and yy == 0 else 2,
                                          0 if z == partitions[2]-1 and zz == partitionCubes[2]-1 else 1],
                                         [0 if x == 0 and xx == 0 else 0,
                                          0 if y == partitions[1]-1 and yy == partitionCubes[1]-1 else 1,
                                          3,
                                          0 if z == partitions[2]-1 and zz == partitionCubes[2]-1 else 1],
                                         [3, 0, 0, 2]]
                                
                                sides.extend(s)
                                
                    varElemNeighborSides[(z*partitions[1] + y)*partitions[0] + x,:,:] = sides
                    
        print 'element_neighbor_sides done'
        
        for z in range(partitions[2]):
            for y in range(partitions[1]):
                for x in range(partitions[0]):
                    orientations = []
                    
                    for zz in range(partitionCubes[2]):
                        for yy in range(partitionCubes[1]):
                            for xx in range(partitionCubes[0]):
                                odd = (xx + yy + zz) % 2
                                
                                if odd:
                                    o = [[0,# if x == 0 and xx == 0 else 0,
                                          0 if z == 0 and zz == 0 else 1,
                                          0,
                                          0],# if y == partitions[1]-1 and yy == partitionCubes[1]-1 else 0],
                                         [0,
                                          0 if z == 0 and zz == 0 else 1,
                                          0,# if y == 0 and yy == 0 else 0,
                                          0 if x == partitions[0]-1 and xx == partitionCubes[0]-1 else 2],
                                         [0,
                                          0,# if y == 0 and yy == 0 else 0,
                                          0,# if x == 0 and xx == 0 else 0,
                                          0 if z == partitions[2]-1 and zz == partitionCubes[2]-1 else 2],
                                         [0,# if x == partitions[0]-1 and xx == partitionCubes[0]-1 else 0,
                                          0,
                                          0,# if y == partitions[1]-1 and yy == partitionCubes[1]-1 else 0,
                                          0],# if z == partitions[2]-1 and zz == partitionCubes[2]-1 else 0],
                                         [0, 0, 0, 0]]
                                else:
                                    o = [[0 if z == 0 and zz == 0 else 2, 
                                          0,# if y == 0 and yy == 0 else 0,
                                          0 if x == 0 and xx == 0 else 2,
                                          0],
                                         [0,
                                          0,# if z == 0 and zz == 0 else 0,
                                          0,# if x == partitions[0]-1 and xx == partitionCubes[0]-1 else 0,
                                          0],# if y == partitions[1]-1 and yy == partitionCubes[1]-1 else 0],
                                         [0,
                                          0,# if x == partitions[0]-1 and xx == partitionCubes[0]-1 else 0,
                                          0,# if y == 0 and yy == 0 else 0,
                                          0 if z == partitions[2]-1 and zz == partitionCubes[2]-1 else 1],
                                         [0,# if x == 0 and xx == 0 else 0,
                                          0,# if y == partitions[1]-1 and yy == partitionCubes[1]-1 else 0,
                                          0,
                                          0 if z == partitions[2]-1 and zz == partitionCubes[2]-1 else 1],
                                         [0, 0, 0, 0]]
                                
                                orientations.extend(o)
                                
                    varElemSideOrientations[(z*partitions[1] + y)*partitions[0] + x,:,:] = orientations
                    
        print 'element_side_orientations done'
        
        for z in range(partitions[2]):
            for y in range(partitions[1]):
                for x in range(partitions[0]):
                    myrank = (z*partitions[1] + y)*partitions[0] + x
                    ranks = []
                    
                    for zz in range(partitionCubes[2]):
                        for yy in range(partitionCubes[1]):
                            for xx in range(partitionCubes[0]):
                                odd = (xx + yy + zz) % 2
                                
                                if odd:
                                    r = [[myrank if xx != 0 else myrank if x == 0 else myrank-1,
                                          myrank if zz != 0 else myrank if z == 0 else myrank-partitions[1]*partitions[0],
                                          myrank,
                                          myrank if yy != partitionCubes[1]-1 else myrank if y == partitions[1]-1 else myrank+partitions[0]],
                                         [myrank,
                                          myrank if zz != 0 else myrank if z == 0 else myrank-partitions[1]*partitions[0],
                                          myrank if yy != 0 else myrank if y == 0 else myrank-partitions[0],
                                          myrank if xx != partitionCubes[0]-1 else myrank if x == partitions[0]-1 else myrank+1],
                                         [myrank,
                                          myrank if yy != 0 else myrank if y == 0 else myrank-partitions[0],
                                          myrank if xx != 0 else myrank if x == 0 else myrank-1,
                                          myrank if zz != partitionCubes[2]-1 else myrank if z == partitions[2]-1 else myrank+partitions[1]*partitions[0]],
                                         [myrank if xx != partitionCubes[0]-1 else myrank if x == partitions[0]-1 else myrank+1,
                                          myrank,
                                          myrank if yy != partitionCubes[1]-1 else myrank if y == partitions[1]-1 else myrank+partitions[0],
                                          myrank if zz != partitionCubes[2]-1 else myrank if z == partitions[2]-1 else myrank+partitions[1]*partitions[0]],
                                         [myrank, myrank, myrank, myrank]]
                                else:
                                    r = [[myrank if zz != 0 else myrank if z == 0 else myrank-partitions[1]*partitions[0], 
                                          myrank if yy != 0 else myrank if y == 0 else myrank-partitions[0],
                                          myrank if xx != 0 else myrank if x == 0 else myrank-1,
                                          myrank],
                                         [myrank,
                                          myrank if zz != 0 else myrank if z == 0 else myrank-partitions[1]*partitions[0],
                                          myrank if xx != partitionCubes[0]-1 else myrank if x == partitions[0]-1 else myrank+1,
                                          myrank if yy != partitionCubes[1]-1 else myrank if y == partitions[1]-1 else myrank+partitions[0]],
                                         [myrank,
                                          myrank if xx != partitionCubes[0]-1 else myrank if x == partitions[0]-1 else myrank+1,
                                          myrank if yy != 0 else myrank if y == 0 else myrank-partitions[0],
                                          myrank if zz != partitionCubes[2]-1 else myrank if z == partitions[2]-1 else myrank+partitions[1]*partitions[0]],
                                         [myrank if xx != 0 else myrank if x == 0 else myrank-1,
                                          myrank if yy != partitionCubes[1]-1 else myrank if y == partitions[1]-1 else myrank+partitions[0],
                                          myrank,
                                          myrank if zz != partitionCubes[2]-1 else myrank if z == partitions[2]-1 else myrank+partitions[1]*partitions[0]],
                                         [myrank, myrank, myrank, myrank]]
                                
                                ranks.extend(r)
                       
                    lastIndices = [-1]*6
                    bndLocalIds = [None]*6
                    for i in range(6):
                        bndLocalIds[i] = []
                    def mpiIndex(e, r):
                        if r == myrank:
                            return 0
                        if r == myrank-partitions[0]*partitions[1]:
                            lastIndices[0] += 1
                            if ((lastIndices[0] % (partitionCubes[0]*4)) / (partitionCubes[0]*2) == 0 and lastIndices[0] % 4 == 2) \
                                    or ((lastIndices[0] % (partitionCubes[0]*4)) / (partitionCubes[0]*2) == 1 and (lastIndices[0] % 4) == 0):
                                bndLocalIds[0].extend([None, e])
                                return lastIndices[0]+1
                            elif ((lastIndices[0] % (partitionCubes[0]*4)) / (partitionCubes[0]*2) == 0 and lastIndices[0] % 4 == 3) \
                                    or ((lastIndices[0] % (partitionCubes[0]*4)) / (partitionCubes[0]*2) == 1 and (lastIndices[0] % 4) == 1):
                                bndLocalIds[0][len(bndLocalIds[0])-2] = e
                                return lastIndices[0]-1
                            bndLocalIds[0].append(e)
                            return lastIndices[0]
                        if r == myrank-partitions[0]:
                            lastIndices[1] += 1
                            bndLocalIds[1].append(e)
                            return lastIndices[1]
                        if r == myrank-1:
                            lastIndices[2] += 1
                            bndLocalIds[2].append(e)
                            return lastIndices[2]
                        if r == myrank+1:
                            lastIndices[3] += 1
                            bndLocalIds[3].append(e)
                            return lastIndices[3]
                        if r == myrank+partitions[0]:
                            lastIndices[4] += 1
                            bndLocalIds[4].append(e)
                            return lastIndices[4]
                        if r == myrank+partitions[0]*partitions[1]:
                            lastIndices[5] += 1
                            bndLocalIds[5].append(e)
                            return lastIndices[5]
                        raise Exception('Invalid neighbor rank found')
                    indices = map(lambda (e, r): map(functools.partial(mpiIndex, e), r), enumerate(ranks))
                    varElemNeighborRanks[(z*partitions[1] + y)*partitions[0] + x,:,:] = ranks
                    varElemMPIIndices[(z*partitions[1] + y)*partitions[0] + x,:,:] = indices
                            
                    k = 0
                    for i in range(6):
                        if bndLocalIds[i]:
                            varBndElemLocalIds[(z*partitions[1] + y)*partitions[0] + x,k,0:len(bndLocalIds[i])] = bndLocalIds[i]
                            k += 1
                            
        print 'element_neighbor_ranks, element_mpi_indices and boundary_element_localids done'
                    
        varVrtxSize[:] = vrtxSize
        
        for z in range(partitions[2]):
            for y in range(partitions[1]):
                for x in range(partitions[0]):
                    globalVrtxIds = []
                    
                    for zz in range(partitionCubes[2]):
                        for yy in range(partitionCubes[1]):
                            for xx in range(partitionCubes[0]):
                                for i in range(5):
                                    globalVrtxIds += mesh.elements()[(((z*partitionCubes[2]+zz)*mesh.size()[1]*2 + (y*partitionCubes[1]+yy))*mesh.size()[0]*2 + (x*partitionCubes[0]+xx))*5 + i]
                                    
                    coords = map(lambda v: mesh.coords()[v], unique(globalVrtxIds))
                    varVrtxCoords[(z*partitions[1] + y)*partitions[0] + x,0:vrtxSize,:] = coords
                    
        print 'vertex_coordinates done'
                    
        boundaries = [None]*numPartitions
        for z in range(partitions[2]):
            for y in range(partitions[1]):
                for x in range(partitions[0]):
                    b = [True]*6
                    if x == 0:
                        b[2] = False
                    if y == 0:
                        b[1] = False
                    if z == 0:
                        b[0] = False
                    if x == partitions[0]-1:
                        b[3] = False
                    if y == partitions[1]-1:
                        b[4] = False
                    if z == partitions[2]-1:
                        b[5] = False
                    boundaries[(z*partitions[1] + y)*partitions[0] + x] = b
        
        varBndSize[:] = map(sum, boundaries)
        
        boundaryRankDiff = [1, partitions[0], partitions[0]*partitions[1]]
        boundaryRankDiff = map(lambda i: -i, boundaryRankDiff[::-1]) + boundaryRankDiff
        
        for i in range(numPartitions):
            bndElemSize = []
            bndElemRank = []
            
            for j in range(6):
                if boundaries[i][j]:
                    bndElemSize.append(boundarySizes[j])
                    bndElemRank.append(boundaryRankDiff[j]+i)
                    
            varBndElemSize[i,0:len(bndElemSize)] = bndElemSize
            varBndElemRank[i,0:len(bndElemRank)] = bndElemRank
        
        file.close()
