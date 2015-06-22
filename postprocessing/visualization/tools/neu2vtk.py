#!/usr/bin/env python
##
# @file
# This file is part of SeisSol.
#
# @author Martin van Driel (driel AT geophysik.uni-muenchen.de, http://www.geophysik.lmu.de/Members/driel)
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
# neu2vtk <neufile>
#
# creates vtk-file 'neufile.vtk'
#
# for tetrahedral meshes only!!!!

import sys

fi = open(sys.argv[1], 'r')
fo = open(sys.argv[1] + '.vtk', 'w')

for i in range(7):
    s = fi.readline()
s = s.split()

npts = int(s[0])
nelem = int(s[1])

fo.write('# vtk DataFile Version 3.0\n')
fo.write('vtk ' + sys.argv[1] + '\n')
fo.write('ASCII\n')
fo.write('DATASET UNSTRUCTURED_GRID\n')

fo.write('POINTS ' + str(npts) + ' float\n')

fi.readline()
fi.readline()

for i in range(npts):
    s = fi.readline().split()
    fo.write(s[1] + ' ' + s[2] + ' ' + s[3] + '\n')

fi.readline()
fi.readline()

fo.write('CELLS ' + str(nelem) + ' ' + str(nelem*5) + '\n')

for i in range(nelem):
    s = fi.readline().split()
    s3 = str(int(s[3]) - 1)
    s4 = str(int(s[4]) - 1)
    s5 = str(int(s[5]) - 1)
    s6 = str(int(s[6]) - 1)
    fo.write('4 ' + s3 + ' ' + s4 + ' ' + s5 + ' ' + s6 + '\n')

fo.write('CELL_TYPES ' + str(nelem) + '\n')

for i in range(nelem):
    fo.write('10\n')

fi.close()
fo.close()
