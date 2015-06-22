#! /usr/bin/python
##
# @file
# This file is part of SeisSol.
#
# @author Stefan Wenk (wenk AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/wenk)
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
# Creates vtk file for receiver coordinates from SeisSol parameter file
# 
# usage: 
# python generate_rec_vtk <SeisSol parameter file>
# i.e. python generate_rec_vtk.py PARAMETER.par
#
# output:
# rec.vtp receiver coordinates in vtk format

import sys

try:
    fi = open(sys.argv[1], 'r')
except:
    sys.exit('please give SeisSol parameter file and call: python generate_rec_vtk.py PARAMETER.par')

count = 0
for line in fi:
    if line[0] == '<':
        count = count+1
        if count == 27:
            break

for i,line in enumerate(fi):
    if i == 10:
        nrec = int(line.split()[0])
        break

a = []
for i,line in enumerate(fi):
    a = a+[line]
    if i == nrec-1:
        break

fo = open('rec.vtp', 'w')

fo.write('<VTKFile type="PolyData" version="0.1">\n')
fo.write('<PolyData>\n')
fo.write('<Piece NumberOfPoints="'+str(nrec)+'" NumberOfVerts="'+str(nrec)+'">\n')
fo.write('<Points>\n')
fo.write('<DataArray NumberOfComponents="3" format="ascii" type="Float32">\n')
for i in a:
    fo.write(i)

fo.write('</DataArray>\n')
fo.write('</Points>\n')
fo.write('<Verts>\n')
fo.write('<DataArray type="Int32" Name="connectivity" format="ascii">\n')
for i in range(nrec):
    fo.write(str(i)+' ')

fo.write('\n')
fo.write('</DataArray>')
fo.write('<DataArray type="Int32" Name="offsets" format="ascii">\n')
for i in range(1,nrec+1):
    fo.write(str(i)+' ')

fo.write('\n')
fo.write('</DataArray>\n')
fo.write('</Verts>\n')
fo.write('<PointData>\n')
fo.write('<DataArray format="ascii" type="Float32" Name="val">\n')
for i in range(1,nrec+1):
    fo.write(str(i)+'\n')

fo.write('</DataArray>\n')
fo.write('</PointData>\n')
fo.write('</Piece>\n')
fo.write('</PolyData>\n')
fo.write('</VTKFile>\n')
