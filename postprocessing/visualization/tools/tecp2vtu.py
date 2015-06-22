#!/usr/bin/env python
##
# @file
# This file is part of SeisSol.
#
# @author Alice Gabriel (gabriel AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/gabriel)
#
# @section LICENSE
# Copyright (c) SeisSol Group
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
# usage: 
# python tecp2vtu <tecpfile>
#
# creates vtu-file 'tecpfile.vtu' for one tecplot output file
#
# for tetrahedral meshes only!

import sys
import numpy as np

fi = open(sys.argv[1], 'r')
fo = open(sys.argv[1] + '.vtu', 'w')

title = fi.readline().split('"')[1].strip()

vars = fi.readline().split()[2:]

for i, var in enumerate(vars):
    vars[i] = var.strip('"')

nvar = len(vars)
s = fi.readline().split()
npts = int(s[2])
nelem = int(s[4])

data = np.empty((nvar, npts))

for i in range(npts):
    data[:,i] = np.array(fi.readline().split(), dtype=float)

elems = np.zeros((4, nelem), dtype=int)
for i in range(nelem):
    elems[:,i] = np.array(fi.readline().split(), dtype=int) - 1

fi.close()

fo.write('<VTKFile type=\'UnstructuredGrid\' version=\'0.1\' byte_order=\'BigEndian\'>\n')
fo.write('<UnstructuredGrid>\n')
fo.write('<Piece NumberOfPoints=\'' + str(npts) + '\' NumberOfCells=\'' +  str(nelem) + '\'>\n')
fo.write('<Points>\n')
fo.write('<DataArray type=\'Float64\' NumberOfComponents=\'3\' Format=\'ascii\'>\n')

for i in range(npts):
    fo.write('   ' + str('{: .16E}'.format(data[0,i])) + '   ' + str('{: .16E}'.format(data[1,i])) + '   ' + str('{: .16E}'.format(data[2,i])) + '\n')

fo.write('</DataArray>\n')
fo.write('</Points>\n')
fo.write('<Cells>\n')
fo.write('<DataArray type=\'Int32\' Name=\'connectivity\' Format=\'ascii\'>\n')

for i in range(nelem):
     fo.write(' ' + str('{:6d}'.format(elems[0,i])) + ' ' + str('{:6d}'.format(elems[1,i])) + ' ' + str('{:6d}'.format(elems[2,i])) + ' ' + str('{:6d}'.format(elems[3,i])) + '\n')
fo.write('</DataArray>\n')
fo.write('<DataArray type=\'Int32\' Name=\'offsets\' Format=\'ascii\'>\n')

for i in range(nelem):
    fo.write(' ' + str('{:6d}'.format((i+1)*4)) + '\n')

fo.write('</DataArray>\n')
fo.write('<DataArray type=\'UInt8\' Name=\'types\' Format=\'ascii\'>\n')

for i in range(nelem):
    fo.write(' ' + str('{:6d}'.format(10)) + '\n')

fo.write('</DataArray>\n')
fo.write('</Cells>\n')

fo.write('<PointData Scalars=\'scalars\'>\n')

for j, var in enumerate(vars[3:]):
  
    fo.write('<DataArray type=\'Float64\' Name=\'' + var + '\' Format=\'ascii\'>\n') 	

    for i in range(npts):
        fo.write('   ' + str('{: .16E}'.format(data[j+3,i])) + '\n')
    
    fo.write('</DataArray>\n')


fo.write('</PointData>\n')
fo.write('</Piece>\n')
fo.write('</UnstructuredGrid>\n')
fo.write('</VTKFile>\n')
fo.close()
