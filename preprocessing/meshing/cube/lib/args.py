#! /usr/bin/python
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

import argparse
import operator
import sys

class Args:
    """Parses the command line arguments using
    argparse.ArgumentParser"""
    
    def __init__(self):
        parser = argparse.ArgumentParser(prog='cubeGenerator')
        parser.add_argument('-s', '--size', type=int,
            help='number of cubes in each dimensions')
        parser.add_argument('-x', '--size-x', type=int,
            help='number of cubes in x dimension')
        parser.add_argument('-y', '--size-y', type=int,
            help='number of cubes in y dimension')
        parser.add_argument('-z', '--size-z', type=int,
            help='number of cubes in z dimension')
        parser.add_argument('--px', type=int,
            help='number of partitions x dimension')
        parser.add_argument('--py', type=int,
            help='number of partitions in y dimension')
        parser.add_argument('--pz', type=int,
            help='number of partitions in z dimension')
        parser.add_argument('-o', '--output', type=argparse.FileType('w'),
            required=True, help='output file for resulting Gambit or netCDF mesh')
        parser.add_argument('-n', '--netcdf', action='store_true',
            help='Create netCDF files'),
        parser.add_argument('-b', '--boundary', type=int,
            help='boundary condition (default 106 = periodic)', default=106)
        parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
        
        try:
            # Parse cmd line options
            self.__options = parser.parse_args(sys.argv[1:])
            if self.__options.size == None \
                and (self.__options.size_x == None \
                     or self.__options.size_y == None \
                     or self.__options.size_z == None):
                raise IOError('Either size or x/y/z is required')
            
            if self.__options.netcdf:                
                partitions = [getattr(self.__options, 'p'+i) for i in ['x', 'y', 'z']]
                p = sum(p != None for p in partitions)
                if p != 3:
                    raise IOError('Number of partitions is required for netCDF files')
            
        except IOError, e:
            parser.error(str(e))
            
        # Set correct size options
        self.__options.size = [self.__options.size for _ in range(3)]
        
        if self.__options.size_x != None:
            self.__options.size[0] = self.__options.size_x
        if self.__options.size_y != None:
            self.__options.size[1] = self.__options.size_y
        if self.__options.size_z != None:
            self.__options.size[2] = self.__options.size_z
    
    def size(self):
        return self.__options.size
    
    def outputFile(self):
        return self.__options.output
    
    def boundary(self):
        return self.__options.boundary
    
    def netcdf(self):
        return self.__options.netcdf
    
    def partitions(self):
        return (self.__options.px, self.__options.py, self.__options.pz)
