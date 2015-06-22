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
# python pickpoint2mseed <root-filename>
#
# creates two files:
# <root-filename>.coords contains the station coordinates
# <root-filename>.mseed contains the seismograms in mini-seed format
#
# requires obspy (http://www.obspy.org)

import numpy as np
import sys
import glob
try:
    from obspy.core import Trace, Stream, UTCDateTime
except:
    sys.exit('please install obspy (http://www.obspy.org)')

t = UTCDateTime(0)

try:
    files = sys.argv[1]
except:
    sys.exit('usage: python pickpoint2mseed <root-filename>')

traces = []

f = open(files + '.coords', 'w')

for file in glob.iglob(files + '*.dat'):
    print file
    header = [x.replace("\"","").replace(","," ").replace("="," ").split() for x in open(file).readlines()[:5]]
    stationID = header[0][-1]
    channels = header[1][2:]
    x1 = header[2][2]
    x2 = header[3][2]
    x3 = header[4][2]

    if x1[0] != '-':
        x1 = ' ' + x1
    if x2[0] != '-':
        x2 = ' ' + x2
    if x3[0] != '-':
        x3 = ' ' + x3
    
    f.write(stationID[-4:] + ' ' + x1 + ' ' + x2 + ' ' + x3 + '\n')
    
    dat = np.loadtxt(file, skiprows=5)
    
    for i, chan in enumerate(channels):
        stats = {'network': 'SG', 
                 'station': stationID[-4:], 
                 'location': '',
                 'channel': chan, 
                 'npts': len(dat[:,i+1]), 
                 'sampling_rate': 1./(dat[1,0] - dat[0,0]),
                 'starttime': t,
                 'mseed' : {'dataquality': 'D'}}
        traces.append(Trace(data=dat[:,i+1], header=stats))


f.close()

st = Stream(traces)
st.sort()
print st
fname =  files + '.mseed'
print fname
st.write(fname, format='MSEED')
