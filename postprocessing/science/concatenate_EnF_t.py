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
import subprocess
import sys
from multiprocessing import Pool,cpu_count,Manager
import time

def loadEnF(args):
#####function used for loading the data using pool_async
   i,q  = args
   En = np.loadtxt(filelist[i], skiprows=1)
   myres = np.zeros((2*ndt,))
   myres[0:ndt] =  np.repeat(En[:,1], int(round(dt[i]/dt0)))[0:ndt]
   myres[ndt:2*ndt] =  np.repeat(En[:,2], int(round(dt[i]/dt0)))[0:ndt]
   if q!=0:
      q.put(i)
   return myres


parser = argparse.ArgumentParser(description='concatenate fault energy rate files, write the result in a file and if asked plot it')
parser.add_argument('prefix', help='folder/prefix')
parser.add_argument('--plot', dest='display_plot', action='store_true', help='display a plot of the moment rate')
parser.add_argument('--MP', nargs=1, metavar=('ncpu'), default=([1]), help='use np.pool to speed-up calculations' ,type=int)
args = parser.parse_args()

filelist = glob.glob(args.prefix+'-EnF_t*')
print('%d files found' %len(filelist))

tmax = sys.float_info.max
dt=np.zeros(len(filelist))
for i, fname in enumerate(filelist):
   #check for tmax
   line = subprocess.check_output(['tail', '-1', fname])
   tmax = min(tmax, float(line.split()[0]))
   #check for dt
   fid = open(fname)
   fid.readline()
   line = fid.readline()
   line = fid.readline()
   fid.close()
   dt[i] = float(line.split()[0])

dt0 = min(dt)
print('tfmin=%f, dtmin=%f' %(tmax,dt0))
print('dt/dt0', dt/dt0)

ndt = int(round(tmax/dt0+1))
En_conc = np.zeros((ndt,3))

En_conc[:,0]=np.linspace(0, tmax, ndt)

nprocs  = args.MP[0]
assert(nprocs<=cpu_count())
pool = Pool(processes=nprocs)
m = Manager()
q = m.Queue()
N=len(filelist)
inputs = list(range(0,N))
args2 = [(i, q) for i in inputs]
Result = pool.map_async(loadEnF, args2)
pool.close()
while (True):
  if (Result.ready()): break
  remaining = N+1 - q.qsize()
  print("Waiting for", remaining, "tasks to complete...")
  time.sleep(2.0)
a = np.sum(np.array(Result.get()),axis=0)
print(a)
print(ndt, a.shape)
a = a.reshape((ndt,2),order='F')
En_conc[:,1:3] = En_conc[:,1:3] + a

print('Moment rate:')
print(En_conc[:,1])

print('Frictional energy rate')
print(En_conc[:,2])

#calculate moment magnitude
M0 = np.trapz(En_conc[:,1], x=En_conc[:,0])
Mw = 2.*log10(M0)/3.-6.07
print('Moment magnitude: %f (M0 = %e)' % (Mw, M0))

#calculate total frictional energy release
FricEn_total = np.trapz(En_conc[:,2], x=En_conc[:,0])
print('Frictional energy release: %f ' % (FricEn_total))

FricEn = np.zeros(ndt)
dt = dt0
for idt in range(ndt):
    FricEn[idt] = FricEn[idt-1] + dt*En_conc[idt,2]

print('Frictional energy over time:  ') 
print(FricEn)

#save data
#time moment-rate frictional-energy-rate frictional-energy
output = np.c_[En_conc[:,0], En_conc[:,1], En_conc[:,2],  FricEn[:]]
np.savetxt(args.prefix+'-EnF_0t-all.dat', output)

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


