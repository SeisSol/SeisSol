#!/usr/bin/env python
##
# @file
# This file is part of SeisSol.
#
# @author Martin van Driel (driel@geophysik.uni-muenchen.de)
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
#  Purpose: compute distance to faces and nodes in the reference
#           tetrahedron
#---------------------------------------------------------------------
#
# usage:
#   python distances.py
#
# input:
#   assumes to find a file 'xi_eta_zeta.txt' in the working directory
#   containing the reference coordinates of the stations
#
# output:
#   produces to files 'dist_to_nodes.txt' and 'dist_to_faces.txt' containg the
#   distance to the closest node an face of each station. Thresholds for
#   warnings on the command line can be defined below.

import numpy as np

# define warning threshold:
# for nodes:
en = 1e-2
# for faces:
ef = 1e-5

m = np.loadtxt('stations_xi_eta_zeta.txt')
dnodes = np.empty((m.shape[0], 4))
dfaces = np.empty((m.shape[0], 4))

# compute distance to closest node
for i in np.arange(m.shape[0]):
    dist1 = (m[i,0]**2 + m[i,1]**2 + m[i,2]**2)**.5
    dist2 = ((m[i,0]-1)**2 + m[i,1]**2 + m[i,2]**2)**.5
    dist3 = (m[i,0]**2 + (m[i,1]-1)**2 + m[i,2]**2)**.5
    dist4 = (m[i,0]**2 + m[i,1]**2 + (m[i,2]-1)**2)**.5
    dnodes[i] = np.sort((dist1, dist2, dist3, dist4))
    if dnodes[i,0] < en:
        print 'station %d is close to a node, distance is: %1.2e' % (i,dnodes[i,0])

# write to file
np.savetxt('dist_to_nodes.txt', dnodes)

nv4 = np.matrix([1.,1.,1.]) / np.sqrt(3)

# compute distance to closest face 
for i in np.arange(m.shape[0]):
    vec = np.matrix([m[i,0]-1,m[i,1],m[i,2]])
    dist4 = (nv4 * vec.T)[0,0]
    dfaces[i] = np.sort((m[i,0], m[i,1], m[i,2], np.abs(dist4)))
    if dfaces[i,0] < en:
        print 'station %d is close to a face, distance is: %1.2e' % (i,dfaces[i,0])

# write to file
np.savetxt('dist_to_faces.txt', dfaces)
