##
# @file
# This file is part of SeisSol.
#
# @author Thomas Ulrich  
#
# @section LICENSE
# Copyright (c) 2016, SeisSol Group
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
#

#Author: Thomas Ulrich
#Date: 29.09.17
#aim: 
#1 Read time history from free surface output in either hdf5 or posix format
#2 compute ground motion parameter (PGA,PGV,PGD, SA(T))
#3 store output in a hdf5 file readable by paraview

import sys
import os
import h5py
import numpy as np
import argparse
from multiprocessing import Pool,cpu_count,Manager
import time
import lxml.etree as ET
import seissolxdmf
from scipy import signal
from scipy.integrate import cumtrapz

sys.path.append("%s/gmpe-smtk/" %(os.path.dirname(sys.argv[0])))
try:
   from smtk.intensity_measures import gmrotipp
except ImportError:
   print('module smtk not found: please follow instruction on the readme (README.md)')
   raise

def chunk(xs, n):
    '''Split the list, xs, into n evenly sized chunks'''
    L = len(xs)
    assert 0 < n <= L
    s, r = divmod(L, n)
    t = s + 1
    return ([xs[p:p+t] for p in range(0, r*t, t)] +
            [xs[p:p+s] for p in range(r*t, L, s)])


def low_pass_filter(waveform, fs, cutoff_freq):
    " applied 2 pass zero-phase order 2 low pass filter on data"
    order=2
    b,a = signal.butter(order, cutoff_freq, 'low', fs=fs)
    return signal.filtfilt(b, a, waveform)

def compute_cav_gmrot(acceleration_x, time_step_x, acceleration_y, time_step_y, angles, percentile):
    """ compute the cumulative velocity using gmrot """
    from smtk.intensity_measures import get_cav, rotate_horizontal
    cav_theta = np.zeros(len(angles), dtype=float)
    for iloc, theta in enumerate(angles):
        if iloc == 0:
            cav_theta[iloc] = np.sqrt(get_cav(acceleration_x, time_step_x) * 
                    get_cav(acceleration_y, time_step_y))
        else:
            rot_x, rot_y = rotate_horizontal(acceleration_x, acceleration_y, theta)
            cav_theta[iloc] = np.sqrt(get_cav(rot_x, time_step_x) * 
                    get_cav(rot_y, time_step_y))
    return np.percentile(cav_theta, percentile)

def gmrotdpp_withPG(acceleration_x, time_step_x, acceleration_y, time_step_y, periods,
        percentile, damping=0.05, units="cm/s/s", method="Nigam-Jennings"):
    """
    modified from gmrotdpp to also return gmrotdpp(PGA, PGV and PGD)
    This is much faster than gmrotdpp_slow
    """
    from smtk.intensity_measures import get_response_spectrum,equalise_series,rotate_horizontal
    if (percentile > 100. + 1E-9) or (percentile < 0.):
        raise ValueError("Percentile for GMRotDpp must be between 0. and 100.")
    # Get the time-series corresponding to the SDOF
    sax, _, x_a, _, _ = get_response_spectrum(acceleration_x,
                                              time_step_x,
                                              periods, damping,
                                              units, method)
    say, _, y_a, _, _ = get_response_spectrum(acceleration_y,
                                              time_step_y,
                                              periods, damping,
                                              units, method)
    x_a, y_a = equalise_series(x_a, y_a)

    #TU: this is the part I m adding
    #compute vel and disp from acceleration and
    #add to the spectral acceleration time series
    velocity_x = time_step_x * cumtrapz(acceleration_x[0:-1], initial=0.)
    displacement_x = time_step_x * cumtrapz(velocity_x, initial=0.)
    x_a = np.column_stack((acceleration_x[0:-1], velocity_x, displacement_x, x_a))

    velocity_y = time_step_y * cumtrapz(acceleration_y[0:-1], initial=0.)
    displacement_y = time_step_y * cumtrapz(velocity_y, initial=0.)
    y_a = np.column_stack((acceleration_y[0:-1], velocity_y, displacement_y, y_a))

    angles = np.arange(0., 90., 1.)
    max_a_theta = np.zeros([len(angles), len(periods)+3], dtype=float)
    max_a_theta[0, :] = np.sqrt(np.max(np.fabs(x_a), axis=0) *
                                np.max(np.fabs(y_a), axis=0))
    for iloc, theta in enumerate(angles):
        if iloc == 0:
            max_a_theta[iloc, :] = np.sqrt(np.max(np.fabs(x_a), axis=0) *
                                           np.max(np.fabs(y_a), axis=0))
        else:
            rot_x, rot_y = rotate_horizontal(x_a, y_a, theta)
            max_a_theta[iloc, :] = np.sqrt(np.max(np.fabs(rot_x), axis=0) *
                                           np.max(np.fabs(rot_y), axis=0))

    gmrotd = np.percentile(max_a_theta, percentile, axis=0)

    res =  {"PGA": gmrotd[0],
            "PGV": gmrotd[1],
            "PGD": gmrotd[2],
            "Acceleration": gmrotd[3:]}

    if args.CAV:
        cav = compute_cav_gmrot(acceleration_x, time_step_x, acceleration_y, time_step_y, angles, percentile)
        res['CAV']=cav

    return res

def ComputeGroundMotionParameters(args2):
   # This function encapsulate the ground motion estimates calculation
   # q is the queue required for keeping track of the remaining tasks of pool.map_async
   ir,q  = args2

   #Xdirection
   velocity = u[ir,:]
   if args.lowpass:
       velocity = low_pass_filter(velocity, fs=1./dt, cutoff_freq=args.lowpass[0])
   acceleration_x = np.gradient(velocity, dt)

   #Ydirection
   velocity = v[ir,:]
   if args.lowpass:
       velocity = low_pass_filter(velocity, fs=1./dt, cutoff_freq=args.lowpass[0])
   acceleration_y = np.gradient(velocity, dt)


   #compute SA(T) and PGA,...
   res = myfunc(acceleration_x, dt, acceleration_y, dt, periods,
            percentile=50, damping=0.05, units="cm/s/s", method="Nigam-Jennings")
   

   myres = [res['PGA'],res['PGV'],res['PGD']]
   if args.CAV:
        myres.append(res['CAV'])
   myres += list(res['Acceleration'])

   if q!=0:
      q.put(ir)
   return myres

start = time.time()
# parsing python arguments
parser = argparse.ArgumentParser(description='compute ground motions parameters maps (e.g. PGA, PGV, etc) from free surface output (this script uses mpi4py and np.pool to speed-up calculations)')
parser.add_argument('filename', help='surface output (xdmf file)')
parser.add_argument('--MP', nargs=1, metavar=('nthreads'), default=([1]), help='use np.pool to speed-up calculations' ,type=int)
parser.add_argument('--noMPI',  dest='noMPI', action='store_true' , default = False, help='run on one node (skip mpi4py)')
parser.add_argument('--periods', nargs='+', help='list of periods for computing SA(T)' ,type=float)
parser.add_argument('--ipp',  dest='ipp', action='store_true' , default = False, help='use gmrotipp rather than the period dependant estimate (~30 times slower)')
parser.add_argument('--CAV',  dest='CAV', action='store_true' , default = False, help='compute cumulated absolute velocity (slow down the code by a factor ~3)')
parser.add_argument('--lowpass', nargs=1, metavar=('cutoff_freq'), help='low pass filter the velocity waveforms prior to computing the GME. cutoff_freq in Hz' ,type=float)
args = parser.parse_args()

if args.noMPI:
   nranks = 1
   irank = 0
else:
   from mpi4py import MPI
   comm = MPI.COMM_WORLD
   nranks = comm.size
   irank = comm.rank

if args.ipp:
   myfunc = gmrotipp
else:
   myfunc = gmrotdpp_withPG

sx = seissolxdmf.seissolxdmf(args.filename)
dt = sx.ReadTimeStep()
nElements = sx.nElements

#split the input array in nranks inputs
inputs0 = np.arange(0, nElements)
inputsi = chunk(inputs0,nranks)
inputs = inputsi[irank]
nElements_i = len(inputs)

#This reads only the chunk of horizontal velocity data required by the rank
u = sx.ReadDataChunk('u', firstElement=inputsi[irank][0], nchunk=nElements_i, idt=-1)
v = sx.ReadDataChunk('v', firstElement=inputsi[irank][0], nchunk=nElements_i, idt=-1)

u=np.transpose(u)
v=np.transpose(v)
ndt = u.shape[1]
if irank==0:
   print("done reading surface data: %f" % (time.time()-start))

#parameters for response spectrum
damping=0.05
units="cm/s/s"

if args.periods is None:
   #periods = np.array([0.100,0.106,0.113,0.120,0.128,0.136,0.145,0.154,0.164,0.175,0.186,0.198,0.211,0.224,0.239,0.254,0.270,0.288,0.306,0.326,0.347,0.369,0.392,0.418,0.444,0.473,0.503,0.535,0.570,0.606,0.645,0.687,0.731,0.777,0.827,0.880,0.937,0.997,1.061,1.129,1.201,1.278,1.360,1.447,1.540,1.639,1.744,1.856,1.975,2.101,2.236,2.379,2.532,2.694,2.867,3.051,3.247,3.455,3.676,3.912,4.163,4.430,4.714,5.016,5.338,5.680,6.044,6.431,6.844,7.283,7.75,8.25,8.78,9.34,9.94,10.57,11.25,11.97,12.74,13.56,14.43,15.35,16.34,17.38,18.50,19.68,20.95,22.29,23.72,25.24,26.86,28.58,30.41,32.36,34.44,36.65,39.00,41.50,44.16])
   periods = np.array([0.100,0.125,0.25,0.4,0.5,0.75,1,1.5,2,2.5,5])
   if irank==0:
      print("no periods specified: using default periods", periods)
else:
   periods = np.array([float(per) for per in args.periods])  
   if irank==0:
      print("periods", periods)

nTasksPerRank  = args.MP[0]
assert(nTasksPerRank<=cpu_count())
pool = Pool(processes=nTasksPerRank)
m = Manager()
q = m.Queue()


args2 = [(i-inputs[0], q) for i in inputs]
Result = pool.map_async(ComputeGroundMotionParameters, args2)
pool.close()

while (True):
  if (Result.ready()): break
  remaining = nElements_i+1 - q.qsize()
  print(irank, ": Waiting for", remaining, "tasks to complete...")
  time.sleep(10.0)
myres = np.ascontiguousarray(np.array(Result.get()))

if irank==0:
   print("done computing ground motion estimates: %f" % (time.time()-start))

dataName = ['PGA','PGV','PGD']

if args.CAV:
    dataName.append('CAV')

for per in periods:
   dataName.append('SA%06.3fs' %per)

#write the output in a temp file
h5f = h5py.File('temp'+str(irank)+'.h5','w')
for i,dname in enumerate(dataName):
     hdname = "mesh0/"+dname
     h5f.create_dataset(hdname, (nElements_i,), dtype='d')
     h5f[hdname][:] = myres[:,i]
h5f.close()


mypath, fn = os.path.split(args.filename)
prefix=fn.split('-')[-2]

sLowPass=f'_lp{args.lowpass[0]:.1f}' if args.lowpass else ''
prefixGME = f'{prefix}{sLowPass}-GME'

if not args.noMPI:
   comm.Barrier()

if irank==0:
  #Write GME data to hdf5
  h5f = h5py.File(prefixGME + '-surface_cell.h5','w')

  for i,dname in enumerate(dataName):
     hdname = "mesh0/"+dname
     h5f.create_dataset(hdname, (nElements,), dtype='d')
     for j in range(nranks):
        fname = 'temp'+str(j)+'.h5'
        h5fj = h5py.File(fname,'r')
        h5f[hdname][inputsi[j][0]:(1+inputsi[j][-1])] = h5fj[hdname][:]
        h5fj.close()
  
  for j in range(nranks):
     fname = 'temp'+str(j)+'.h5'
     os.remove(fname)

  # Copy connect and geometry: 
  h5f.create_dataset('mesh0/connect', (nElements,3), dtype='i8')
  connect = sx.ReadConnect()
  h5f['mesh0/connect'][:,:] = connect[:,:]
  del connect
  h5f.close()
  h5f = h5py.File(prefixGME + '-surface_vertex.h5','w')
  xyz = sx.ReadGeometry()
  Nvertex=xyz.shape[0]
  h5f.create_dataset('mesh0/geometry', xyz.shape, dtype='d')
  h5f['mesh0/geometry'][:,:] = xyz[:,:]
  del xyz
  h5f.close()

  #This write a new xdmf file
  a0 = ET.Element('Xdmf', attrib={'Version':'2.0'})
  a = ET.SubElement(a0,'Domain')
  b = ET.SubElement(a, 'Grid', attrib={'Name':"TimeSeries", 'GridType':'Collection','CollectionType':'Temporal'})
  c = ET.SubElement(b, 'Grid', attrib={'Name':'step_000000000000', 'GridType':'Uniform'})
  d = ET.SubElement(c, 'Topology', attrib={'TopologyType':'Triangle', 'NumberOfElements':str(nElements)})
  dimension = '%d 3' %nElements
  e = ET.SubElement(d, 'DataItem', attrib={'NumberType':'Int', 'Precision':'8', 'Format':'HDF', 'Dimensions': dimension})
  e.text = "%s-surface_cell.h5:/mesh0/connect" %prefixGME

  d1 = ET.SubElement(c, 'Geometry', attrib={'Name':'geo', 'GeometryType':'XYZ', 'NumberOfElements':str(Nvertex)})
  dimension = '%d 3' %Nvertex
  e1 = ET.SubElement(d1, 'DataItem', attrib={'NumberType':'Float', 'Precision':'8', 'Format':'HDF', 'Dimensions': dimension})
  e1.text = "%s-surface_vertex.h5:/mesh0/geometry" %prefixGME
  d2 = ET.SubElement(c, 'Time', attrib={'Value':'0'})

  for i,dname in enumerate(dataName):
     d3 = ET.SubElement(c, 'Attribute', attrib={'Name':dname, 'Center':'Cell'})
     e3 = ET.SubElement(d3, 'DataItem', attrib={'ItemType':'HyperSlab', 'Dimensions':str(nElements)})
     f1 = ET.SubElement(e3, 'DataItem', attrib={'NumberType':'UInt', 'Precision':'4', 'Format':'XML', 'Dimensions':'3 2'})
     f1.text = '0 0 1 1 1 %d' %nElements
     dimension = '1 %d' %nElements
     f2 = ET.SubElement(e3, 'DataItem', attrib={'NumberType':'Float', 'Precision':'8', 'Format':'HDF', 'Dimensions':dimension})
     f2.text = "%s-surface_cell.h5:/mesh0/%s" %(prefixGME, dname)

  fout=open(prefixGME+'.xdmf','w')
  fout.write('<?xml version="1.0" ?>\n')
  fout.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
  fout.write(ET.tostring(a0, pretty_print=True, encoding="unicode"))
  fout.close()
  print("done writing the file: %f" % (time.time()-start))
