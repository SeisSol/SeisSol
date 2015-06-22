#!/usr/bin/python
##
# @file
# This file is part of SeisSol.
#
# @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
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
# Gathers information about matrix distributions.

import argparse
import sys
from netCDF4 import Dataset

class Args:
    
    def __init__(self):
        parser = argparse.ArgumentParser(prog='netcdfMatrixCounter')
        parser.add_argument("file", type=file, help='netCDF file which should be analysed')
    
        try:
            # Parse cmd line options
            self.__options = parser.parse_args(sys.argv[1:])
            self.__options.file.close() # Let netcdf handle this ...
        except IOError, e:
            parser.error(str(e))
            
    def file(self):
        return self.__options.file.name
    
def progress(val, end_val, bar_length=10):
    percent = float(val) / end_val
    hashes = '#' * int(round(percent * bar_length))
    spaces = ' ' * (bar_length - len(hashes))
    sys.stderr.write("Done: [{0}] {1}%\r".format(hashes + spaces, int(round(percent * 100))))
    sys.stderr.flush()
    
def main():
    # Parse command line arguements
    args = Args()
    
    progress(0, 1)
    
    dataset = Dataset(args.file(), 'r')
    
    dimPartitions = dataset.dimensions['partitions']
    dimElementSides = dataset.dimensions['element_sides']
    
    varElemSize = dataset.variables['element_size']
    varBoundaries = dataset.variables['element_boundaries']
    varNeighborSide = dataset.variables['element_neighbor_sides']
    varOrientation = dataset.variables['element_side_orientations']
    varNeighborRanks = dataset.variables['element_neighbor_ranks']
    
    # Known boundries
    boundaries = {1: 0, 3: 0, 5: 0}
    mpiDrBoundaries = 0
    local = [0]*4
    fluxes = [[[0]*3 for j in range(4)] for i in range(4)]
    
    for i in range(len(dimPartitions)):
        size = varElemSize[i]
        b = varBoundaries[i]
        n = varNeighborSide[i]
        o = varOrientation[i]
        r = varNeighborRanks[i]
        for j in xrange(size):
            for k in range(len(dimElementSides)):
                if b[j][k] == 0:
                    fluxes[k][n[j][k]][o[j][k]] += 1
                    local[k] += 1
                else:
                    boundaries[b[j][k]] += 1
                    if b[j][k] == 1:
                        local[k] += 2
                    elif b[j][k] == 3:
                        if r[j][k] != i:
                            mpiDrBoundaries += 1
                    elif b[j][k]:
                        local[k] += 1
               
        progress(i+1, len(dimPartitions))
        
    sys.stderr.write(' '*25+"\r") # Clear status bar
                    
    print '# Boundary faces:'
    for b,n in boundaries.items():
        print str(b)+':', n
    print 'DR MPI boundary faces:', mpiDrBoundaries    
    print
    
    print 'Local matrices'
    for i in range(4):
        print str(i)+','+str(local[i])
    
    print 'local face,neighbor face,orientation'
    for i in range(4):
        for j in range(4):
            for k in range(3):
                print ','.join(map(str, [i, j, k]))+','+str(fluxes[i][j][k])
    
    dataset.close()

if __name__ == '__main__':
    main()
