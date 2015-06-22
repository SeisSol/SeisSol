#!/usr/bin/env python
##
# @file
# This file is part of SeisSol.
#
# @author Martin van Driel (driel AT geophysik.uni-muenchen.de)
#
# @section LICENSE
# Copyright (c) 2011, Martin van Driel
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
#---------------------------------------------------------------------
#  Purpose: Transform station coordinates to reference tetrahedron
#---------------------------------------------------------------------
#
# usage:
#   python stations2reference.py <mesh.neu>
#
# input:
#   assumes to find a file 'stations.txt' containing the station coordinates in
#   the working directory
#
# output:
#   produces a file 'xi_eta_zeta.txt' containing the reference coordinates of
#   the stations

import numpy as np
import sys

fi = open(sys.argv[1], 'r')

stats = np.loadtxt('stations.txt')

xP = stats[:,0]
yP = stats[:,1]
zP = stats[:,2]

EP = np.zeros(len(xP), dtype=int)
XI = np.zeros(len(xP))
ETA = np.zeros(len(xP))
ZETA = np.zeros(len(xP))
found = np.zeros(len(xP), dtype=bool)

# read mesh

for i in range(7):
    s = fi.readline()
s = s.split()

npts = int(s[0])
nelem = int(s[1])

elems = np.empty((nelem, 4), dtype=int)
points = np.empty((npts, 3), dtype=float)

fi.readline()
fi.readline()

for i in range(npts):
    s = fi.readline().split()
    points[i,0] = float(s[1])
    points[i,1] = float(s[2])
    points[i,2] = float(s[3])

fi.readline()
fi.readline()

for i in range(nelem):
    s = fi.readline().split()
    elems[i,0] = int(s[3]) - 1
    elems[i,1] = int(s[4]) - 1
    elems[i,2] = int(s[5]) - 1
    elems[i,3] = int(s[6]) - 1

fi.close()

x = np.empty(4)
y = np.empty(4)
z = np.empty(4)

# loop over all elements in mesh

percent = np.arange(1,11) * (nelem - 1) / 10
for i in range(nelem):
    if i in percent:
        print '%3d %% done' % (((percent == i)* np.arange(1, 11)).sum() * 10)

    x[0] = points[elems[i,0],0]    
    x[1] = points[elems[i,1],0]    
    x[2] = points[elems[i,2],0]    
    x[3] = points[elems[i,3],0]    
    
    y[0] = points[elems[i,0],1]    
    y[1] = points[elems[i,1],1]    
    y[2] = points[elems[i,2],1]    
    y[3] = points[elems[i,3],1]    
    
    z[0] = points[elems[i,0],2]    
    z[1] = points[elems[i,1],2]    
    z[2] = points[elems[i,2],2]    
    z[3] = points[elems[i,3],2]    

    # transform to reference tetrahedron

    J = -z[0]*x[1]*y[2]+x[0]*z[2]*y[3]+x[0]*z[1]*y[2]+z[3]*x[1]*y[2]+             \
         z[2]*x[3]*y[1]+z[1]*y[3]*x[2]+z[3]*y[0]*x[2]-z[2]*y[0]*x[3]+z[2]*        \
         y[0]*x[1]-x[0]*z[1]*y[3]+x[0]*z[3]*y[1]+z[0]*x[1]*y[3]+z[1]*y[0]*        \
         x[3]-z[3]*y[0]*x[1]-x[0]*z[3]*y[2]-z[2]*y[3]*x[1]-z[0]*x[2]*y[3]-        \
         z[3]*x[2]*y[1]-z[1]*x[3]*y[2]+z[0]*x[2]*y[1]-z[1]*y[0]*x[2]-z[0]*x[3]*   \
         y[1]-x[0]*z[2]*y[1]+z[0]*x[3]*y[2]


    xi = ( z[3]*y[2]+z[0]*y[3]-z[2]*y[3]-z[0]*y[2]-z[3]*y[0]+z[2]*y[0])/J*xP +    \
         (-x[0]*z[2]+z[3]*x[0]-z[3]*x[2]+z[0]*x[2]-z[0]*x[3]+x[3]*z[2])/J*yP +    \
         (-y[0]*x[2]+y[2]*x[0]+y[0]*x[3]+y[3]*x[2]-y[2]*x[3]-y[3]*x[0])/J*zP +    \
         ( x[0]*z[2]*y[3]-z[2]*y[0]*x[3]+z[0]*x[3]*y[2]-x[0]*z[3]*y[2] -          \
           z[0]*x[2]*y[3]+z[3]*y[0]*x[2])/J
    
    eta = -( z[3]*y[1]-z[0]*y[1]+z[0]*y[3]+z[1]*y[0]-z[1]*y[3]-z[3]*y[0])/J*xP -  \
           ( z[1]*x[3]-z[1]*x[0]-z[0]*x[3]-z[3]*x[1]+z[3]*x[0]+z[0]*x[1])/J*yP -  \
           (-y[1]*x[3]+y[0]*x[3]+y[3]*x[1]-y[0]*x[1]-y[3]*x[0]+y[1]*x[0])/J*zP -  \
           (-z[1]*y[0]*x[3]+x[0]*z[1]*y[3]-x[0]*z[3]*y[1]-z[0]*x[1]*y[3] +        \
             z[0]*x[3]*y[1]+z[3]*y[0]*x[1])/J
    
    zeta = (z[0]*y[2]-z[2]*y[0]-z[0]*y[1]+z[2]*y[1]+z[1]*y[0]-z[1]*y[2])/J*xP +   \
           (x[0]*z[2]-z[1]*x[0]-z[2]*x[1]+z[0]*x[1]+z[1]*x[2]-z[0]*x[2])/J*yP +   \
           (y[1]*x[0]-y[2]*x[0]+x[1]*y[2]-y[0]*x[1]+y[0]*x[2]-x[2]*y[1])/J*zP +   \
           (-x[0]*z[2]*y[3]-z[3]*x[1]*y[2]+z[0]*x[2]*y[3]-z[2]*x[3]*y[1] -        \
             z[1]*y[3]*x[2]-z[3]*y[0]*x[2]+z[2]*y[0]*x[3]-z[1]*y[0]*x[3] +        \
             x[0]*z[1]*y[3]-x[0]*z[3]*y[1]-z[0]*x[1]*y[3]+z[0]*x[3]*y[1] +        \
             z[3]*y[0]*x[1]+x[0]*z[3]*y[2]+z[2]*y[3]*x[1]+z[3]*x[2]*y[1] -        \
             z[0]*x[3]*y[2]+z[1]*x[3]*y[2]+J)/J
    
    # test if station is in this element
    
    inelem = ((xi <= 0.) | (eta <= 0.) | (zeta <= 0.) | (zeta >= (1.-xi-eta)))
    inelem = (inelem == False)

    EP += inelem * i
    XI += inelem * xi
    ETA += inelem * eta
    ZETA += inelem * zeta
    found += inelem * True

if False in found:
    print 'some station are exactly on a boundary or outside the mesh'

XI[found == False] = np.nan
ETA[found == False] = np.nan
ZETA[found == False] = np.nan

# write to file
np.savetxt('stations_xi_eta_zeta.txt', np.array((XI, ETA, ZETA)).T)
