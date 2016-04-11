##
# @file
# This file is part of SeisSol.
#
# @author Thomas Ulrich (ulrich AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/ulrich)
#
# @section LICENSE
# Copyright (c) 2005-2016, SeisSol Group
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

import glob
import numpy as np
import argparse
from math import log10

parser = argparse.ArgumentParser(description='concatenate fault energy rate files, write the result in a file and if asked plot it')
parser.add_argument('prefix', help='folder/prefix')
parser.add_argument('--plot', dest='display_plot', action='store_true', help='display a plot of the moment rate')
args = parser.parse_args()

filelist = glob.glob(args.prefix+'-EnF_t*')
print '%d files found' %len(filelist)

En_conc = np.loadtxt(filelist[0], skiprows=1)
ndt0 = np.shape(En_conc)[0]

for fid in filelist[1:]:
   En = np.loadtxt(fid, skiprows=1)
   ndt = np.shape(En)[0]
   ndt0 = min(ndt ,ndt0) 
   En_conc[0:ndt0,1] = En_conc[0:ndt0,1] + En[0:ndt0,1]
   En_conc[0:ndt0,2] = En_conc[0:ndt0,2] + En[0:ndt0,2]

En_conc = En_conc[0:ndt0,:]
print 'Moment rate:'
print En_conc[:,1]

print 'Frictional energy rate'
print En_conc[:,2]

#calculate moment magnitude
M0 = np.trapz(En_conc[:,1], x=En_conc[:,0])
Mw = 2.*log10(M0)/3.-6.07
print 'Moment magnitude: %f (M0 = %e)' % (Mw, M0)

#calculate total frictional energy release
FricEn_total = np.trapz(En_conc[:,2], x=En_conc[:,0])
print 'Frictional energy release: %f ' % (FricEn_total)

FricEn = np.zeros(ndt0)
dt = En_conc[2,0]-En_conc[1,0]
for idt in range(ndt0):
    FricEn[idt] = FricEn[idt-1] + dt*En_conc[idt,2]

print 'Frictional energy over time:  ' 
print FricEn

#save data
#time moment-rate frictional-energy-rate frictional-energy
output = np.c_[En_conc[:,0], En_conc[:,1], En_conc[:,2],  FricEn[:]]
np.savetxt(args.prefix+'-EnF_t-all.dat', output)

if args.display_plot:
   import matplotlib.pyplot as plt
   #fig = plt.figure()
   #ax = fig.add_subplot(111)
   f, axarr = plt.subplots(2)
   axarr[0].plot(En_conc[:,0], En_conc[:,1])
   axarr[0].set_ylabel('Nm/s')
   axarr[0].set_title('Seismic moment rate')
   axarr[1].plot(En_conc[:,0], FricEn[:])
   axarr[1].set_ylabel('Joule')
   axarr[1].set_xlabel('time (s)')
   axarr[1].set_title('Frictional energy')
   plt.show()


