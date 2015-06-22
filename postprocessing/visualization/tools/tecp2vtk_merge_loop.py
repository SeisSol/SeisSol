#!/usr/bin/env python
##
# @file
# This file is part of SeisSol.
#
# @author Martin van Driel (driel AT geophysik.uni-muenchen.de, http://www.geophysik.lmu.de/Members/driel)
# @author Simon Kremers (kremers AT geophysik.uni-muenchen.de)
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
# python tecp2vtk_merge [-b] <tecpfile-rootname>
# i.e. python tecp2vtk_merge_loop.py data
#
# creates vtk-files '<timestep>.vtk' for all timesteps in <tecpfile-rootname> subfolders
#
# options:
#   -b: binary output
#
# for tetrahedral meshes only!!! 
# 
# requires pyvtk (http://pypi.python.org/pypi/PyVTK)

import sys
import numpy as np
import glob
import os

try:
    import pyvtk
except:
    sys.exit('please install pyvtk (binary = Falsehttp://pypi.python.org/pypi/PyVTK)')
    
binary = False
tetfiles = []
ts = []

if len(sys.argv) == 2:
    files = sys.argv[1]
elif len(sys.argv) == 3:
    if sys.argv[1] == '-b':
        binary = True
    files = sys.argv[2]
else:
    sys.exit('usage: python tecp2vtk_merge [-b] <tecpfile-rootname> (e.g. "data" )')

for path, dirs, fls in os.walk('./'):
	for d in dirs:
	    if d.startswith(files):
		for f in glob.iglob(os.path.join(path, d, '*.tet.dat')):
			   tetfiles.append(f)
                

if len(tetfiles) == 0:
    sys.exit('no such file: ./' + files + '*/*.tet.dat')

# get number and name of timesteps
for path, dirs, fls in os.walk('./'):
	for d in dirs:
	    if d.startswith(files):
		for f in glob.iglob(os.path.join(path, d, '*.tet.dat')):
		        ts.append(f)
		print ts
	        break

# loop over timesteps
for time in ts:
    t_step = (time.split('-')[2]).split('.')[0] 
    elems = []
    datal = []
    ntelem = 0
    ntpts = 0

    for file in tetfiles:
     if file.endswith(t_step + '.tet.dat'):
        print 'reading ' + file
        fi = open(file, 'r')
    
        title = fi.readline().split('"')[1].strip()
    
        vars = fi.readline().split()[2:]
    
        for i, var in enumerate(vars):
            vars[i] = var.strip('"')
    
        nvar = len(vars)
        s = fi.readline().split()
        npts = int(s[2])
        nelem = int(s[4])
    
        for i in range(npts):
            datal.append(fi.readline().split())
    
        for i in range(nelem):
            elems.append((np.array(fi.readline().split(), dtype=int) - 1 + ntpts).tolist())
    
        ntelem += nelem
        ntpts += npts
    
        fi.close()
    
    data = np.array(datal, dtype=float).T
    elems = np.array(elems, dtype=int).T
    
    pl = []
    for i in np.arange(data.shape[1]):
        pl.append(tuple(data[:3,i]))
    
    el = []
    for i in np.arange(elems.shape[1]):
        el.append(tuple(elems[:,i]))
    
    grid = pyvtk.UnstructuredGrid(pl, tetra=el)
    
    pointdata = pyvtk.PointData()
    for j, var in enumerate(vars[3:]):
        pointdata.append(pyvtk.Scalars(data[j+3,:].tolist(), name=var, lookup_table='default'))
    
    vtk = pyvtk.VtkData(grid, pointdata, title)
    if binary:
        vtk.tofile(t_step + '.vtk', 'binary')
    else:
        vtk.tofile(t_step + '.vtk', 'ascii')
    
    print t_step + '.vtk written'
