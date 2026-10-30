#!/usr/bin/env python3
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
import seissolxdmfwriter as sxw
from scipy import signal
try:
    from scipy.integrate import cumulative_trapezoid
except ImportError:
    from scipy.integrate import cumtrapz as cumulative_trapezoid

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

def butter_bandpass_filter(waveform, fs, freqmin, freqmax, order=2):
    " applied 2 pass zero-phase order 2 band pass filter on data"
    nyquist = 0.5 * fs
    if freqmin <= 0.0 or freqmin >= freqmax or freqmax >= nyquist:
        raise ValueError("Band-pass frequencies must satisfy 0 < fmin < fmax < Nyquist.")
    b, a = signal.butter(order, [freqmin, freqmax], btype='band', fs=fs)
    return signal.filtfilt(b, a, waveform)

def _rotation_penalty_index(max_a_theta, percentile):
    """ Select the RotI/GMRotI angle closest to the requested target spectrum. """
    target = np.percentile(max_a_theta, percentile, axis=0)
    spectra_by_angle = max_a_theta[:, 3:]
    target_spectrum = target[3:]

    if not np.all(np.isfinite(target_spectrum)):
        raise ValueError("The target response spectrum contains non-finite values.")
    valid_periods = target_spectrum > 0.0
    if not np.any(valid_periods):
        raise ValueError("Cannot determine RotI angle: target response spectrum is zero at all periods.")

    relative_difference = (
        spectra_by_angle[:, valid_periods] - target_spectrum[valid_periods]
    ) / target_spectrum[valid_periods]
    penalty = np.mean(relative_difference ** 2, axis=1)
    if not np.all(np.isfinite(penalty)):
        raise ValueError("The RotI angle penalty contains non-finite values.")
    return int(np.argmin(penalty)), target

def compute_cav_gmrot(acceleration_x, time_step_x, acceleration_y, time_step_y, angles, percentile, idx_angle=None):
    """ compute the cumulative velocity using gmrot """
    from smtk.intensity_measures import get_cav, rotate_horizontal
    if idx_angle is not None:
        if idx_angle == 0:
            return np.sqrt(get_cav(acceleration_x, time_step_x) * 
                           get_cav(acceleration_y, time_step_y))
        rot_x, rot_y = rotate_horizontal(acceleration_x, acceleration_y, angles[idx_angle])
        return np.sqrt(get_cav(rot_x, time_step_x) * get_cav(rot_y, time_step_y))

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

def compute_cav_rotd(acceleration_x, time_step_x, acceleration_y, time_step_y, angles, percentile, idx_angle=None):
    """ compute the cumulative velocity using rotd """
    from smtk.intensity_measures import get_cav
    if idx_angle is not None:
        rad = angles[idx_angle]
        rotated = acceleration_x * np.cos(rad) + acceleration_y * np.sin(rad)
        return get_cav(rotated, time_step_x)

    cav_theta = np.zeros(len(angles), dtype=float)
    for iloc, rad in enumerate(angles):
        rotated = acceleration_x * np.cos(rad) + acceleration_y * np.sin(rad)
        cav_theta[iloc] = get_cav(rotated, time_step_x)
    return np.percentile(cav_theta, percentile)

def gmrotdpp_withPG(acceleration_x, time_step_x, acceleration_y, time_step_y, periods,
        percentile=50.0, damping=0.05, units="cm/s/s", method="Nigam-Jennings", ipp=False, include_cav=False):
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
    velocity_x = time_step_x * cumulative_trapezoid(acceleration_x[0:-1], initial=0.)
    displacement_x = time_step_x * cumulative_trapezoid(velocity_x, initial=0.)
    x_a = np.column_stack((acceleration_x[0:-1], velocity_x, displacement_x, x_a))

    velocity_y = time_step_y * cumulative_trapezoid(acceleration_y[0:-1], initial=0.)
    displacement_y = time_step_y * cumulative_trapezoid(velocity_y, initial=0.)
    y_a = np.column_stack((acceleration_y[0:-1], velocity_y, displacement_y, y_a))

    angles = np.arange(0., 90., 1.)
    max_a_theta = np.zeros([len(angles), len(periods)+3], dtype=float)
    for iloc, theta in enumerate(angles):
        if iloc == 0:
            max_a_theta[iloc, :] = np.sqrt(np.max(np.fabs(x_a), axis=0) *
                                           np.max(np.fabs(y_a), axis=0))
        else:
            rot_x, rot_y = rotate_horizontal(x_a, y_a, theta)
            max_a_theta[iloc, :] = np.sqrt(np.max(np.fabs(rot_x), axis=0) *
                                           np.max(np.fabs(rot_y), axis=0))

    if not ipp:
        gmrotd = np.percentile(max_a_theta, percentile, axis=0)
        res =  {"PGA": gmrotd[0],
                "PGV": gmrotd[1],
                "PGD": gmrotd[2],
                "Acceleration": gmrotd[3:]}
        if include_cav or (args and getattr(args, 'CAV', False)):
            cav = compute_cav_gmrot(acceleration_x, time_step_x, acceleration_y, time_step_y, angles, percentile)
            res['CAV'] = cav
    else:
        idx_i_perc, _ = _rotation_penalty_index(max_a_theta, percentile)
        gmroti = max_a_theta[idx_i_perc, :]
        res =  {"PGA": gmroti[0],
                "PGV": gmroti[1],
                "PGD": gmroti[2],
                "Acceleration": gmroti[3:]}
        if include_cav or (args and getattr(args, 'CAV', False)):
            cav = compute_cav_gmrot(acceleration_x, time_step_x, acceleration_y, time_step_y, angles, percentile, idx_angle=idx_i_perc)
            res['CAV'] = cav

    return res

def rotdpp_withPG(acceleration_x, time_step_x, acceleration_y, time_step_y, periods,
        percentile=50.0, damping=0.05, units="cm/s/s", method="Nigam-Jennings", ipp=False, include_cav=False):
    """
    Compute RotD/RotI ground motion parameters (PGA, PGV, PGD, SA(T), optionally CAV)
    """
    from smtk.intensity_measures import get_response_spectrum, equalise_series
    if (percentile > 100. + 1E-9) or (percentile < 0.):
        raise ValueError("Percentile for RotDpp must be between 0. and 100.")
    sax, _, x_a, _, _ = get_response_spectrum(acceleration_x,
                                              time_step_x,
                                              periods, damping,
                                              units, method)
    say, _, y_a, _, _ = get_response_spectrum(acceleration_y,
                                              time_step_y,
                                              periods, damping,
                                              units, method)
    x_a, y_a = equalise_series(x_a, y_a)

    velocity_x = time_step_x * cumulative_trapezoid(acceleration_x[0:-1], initial=0.)
    displacement_x = time_step_x * cumulative_trapezoid(velocity_x, initial=0.)
    x_a = np.column_stack((acceleration_x[0:-1], velocity_x, displacement_x, x_a))

    velocity_y = time_step_y * cumulative_trapezoid(acceleration_y[0:-1], initial=0.)
    displacement_y = time_step_y * cumulative_trapezoid(velocity_y, initial=0.)
    y_a = np.column_stack((acceleration_y[0:-1], velocity_y, displacement_y, y_a))

    angles = np.radians(np.arange(0., 180., 1.))
    max_a_theta = np.zeros([len(angles), len(periods)+3], dtype=float)

    for iloc, rad in enumerate(angles):
        rotated = x_a * np.cos(rad) + y_a * np.sin(rad)
        max_a_theta[iloc, :] = np.max(np.fabs(rotated), axis=0)

    if not ipp:
        rotd = np.percentile(max_a_theta, percentile, axis=0)
        res = {"PGA": rotd[0],
               "PGV": rotd[1],
               "PGD": rotd[2],
               "Acceleration": rotd[3:]}
        if include_cav or (args and getattr(args, 'CAV', False)):
            cav = compute_cav_rotd(acceleration_x, time_step_x, acceleration_y, time_step_y, angles, percentile)
            res['CAV'] = cav
    else:
        idx_i_perc, _ = _rotation_penalty_index(max_a_theta, percentile)
        roti = max_a_theta[idx_i_perc, :]
        res = {"PGA": roti[0],
               "PGV": roti[1],
               "PGD": roti[2],
               "Acceleration": roti[3:]}
        if include_cav or (args and getattr(args, 'CAV', False)):
            cav = compute_cav_rotd(acceleration_x, time_step_x, acceleration_y, time_step_y, angles, percentile, idx_angle=idx_i_perc)
            res['CAV'] = cav

    return res

def ComputeGroundMotionParameters(args2):
   # This function encapsulate the ground motion estimates calculation
   # q is the queue required for keeping track of the remaining tasks of pool.map_async
   ir,q  = args2

   #Xdirection
   velocity = u[ir,:]
   if args.lowpass:
       velocity = low_pass_filter(velocity, fs=1./dt, cutoff_freq=args.lowpass[0])
   elif args.bandpass:
       velocity = butter_bandpass_filter(velocity, fs=1./dt, freqmin=args.bandpass[0], freqmax=args.bandpass[1])
   acceleration_x = np.gradient(velocity, dt)

   #Ydirection
   velocity = v[ir,:]
   if args.lowpass:
       velocity = low_pass_filter(velocity, fs=1./dt, cutoff_freq=args.lowpass[0])
   elif args.bandpass:
       velocity = butter_bandpass_filter(velocity, fs=1./dt, freqmin=args.bandpass[0], freqmax=args.bandpass[1])
   acceleration_y = np.gradient(velocity, dt)


   #compute SA(T) and PGA,...
   res = myfunc(acceleration_x, dt, acceleration_y, dt, periods,
            percentile=args.percentile, damping=0.05, units="cm/s/s", method="Nigam-Jennings", ipp=args.ipp, include_cav=args.CAV)
   

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
parser.add_argument('--metric', choices=['GMRot', 'Rot'], default='GMRot', help='rotation metric family: GMRot or Rot (default: GMRot)')
parser.add_argument('--percentile', type=float, default=50.0, help='percentile across rotation angles (0 to 100, default: 50)')

filter_group = parser.add_mutually_exclusive_group()
filter_group.add_argument('--lowpass', nargs=1, metavar=('cutoff_freq'), help='low pass filter the velocity waveforms prior to computing the GME. cutoff_freq in Hz' ,type=float)
filter_group.add_argument('--bandpass', nargs=2, metavar=('fmin', 'fmax'), help='band pass filter the velocity waveforms prior to computing the GME. fmin and fmax in Hz', type=float)

args = parser.parse_args()

if not 0.0 <= args.percentile <= 100.0:
    parser.error("--percentile must be between 0 and 100.")

if args.noMPI:
   nranks = 1
   irank = 0
else:
   from mpi4py import MPI
   comm = MPI.COMM_WORLD
   nranks = comm.size
   irank = comm.rank

if args.metric == 'Rot':
    myfunc = rotdpp_withPG
else:
    myfunc = gmrotdpp_withPG

sx = seissolxdmf.seissolxdmf(args.filename)
dt = sx.ReadTimeStep()
if args.lowpass and args.lowpass[0] >= 0.5 / dt:
    parser.error("--lowpass cutoff frequency must be below the Nyquist frequency (0.5 / dt).")
if args.bandpass and args.bandpass[1] >= 0.5 / dt:
    parser.error("--bandpass fmax must be below the Nyquist frequency (0.5 / dt).")

nElements = sx.nElements
connect = sx.ReadConnect()

isVolumeData = True if connect.shape[1] == 4 else False

if isVolumeData:
    prefixType = "-volume"
else:
    prefixType = "-surface"

#split the input array in nranks inputs
inputs0 = np.arange(0, nElements)
inputsi = chunk(inputs0,nranks)
inputs = inputsi[irank]
nElements_i = len(inputs)

#This reads only the chunk of horizontal velocity data required by the rank
componentVariables = ['u', 'v', 'w'] if 'u' in sx.ReadAvailableDataFields() else ['v1', 'v2', 'v3']
u = sx.ReadDataChunk(componentVariables[0], firstElement=inputsi[irank][0], nchunk=nElements_i, idt=-1)
v = sx.ReadDataChunk(componentVariables[1], firstElement=inputsi[irank][0], nchunk=nElements_i, idt=-1)

u=np.transpose(u)
v=np.transpose(v)
ndt = u.shape[1]
if irank==0:
    print(f"done reading {prefixType[1:]} data: {time.time()-start}")

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

mypath, fn = os.path.split(args.filename)
prefix = fn.split('-')[-2] if not isVolumeData else fn.split('.')[-2]

perc_str = f"{int(args.percentile)}" if args.percentile.is_integer() else f"{args.percentile:g}"
type_letter = "I" if args.ipp else "D"
effective_metric = f"{args.metric}{type_letter}{perc_str}"

sFilter = ''
if args.lowpass:
    sFilter = f"_lp{args.lowpass[0]:.1f}"
elif args.bandpass:
    sFilter = f"_bp{args.bandpass[0]:.2f}-{args.bandpass[1]:.2f}"

prefixGME = f"{prefix}_{effective_metric}{sFilter}-GME"

if not args.noMPI:
    comm.Barrier()
    myres = comm.gather(myres, root=0)

if irank==0:
    if nranks > 1:
        myres = np.concatenate(myres, axis=0)
    geom = sx.ReadGeometry()
    sxw.write(
        f"{prefixGME}{prefixType}",
        geom,
        connect,
        {name: myres[:,i] for i, name in enumerate(dataName)},
        dictTime=None,
    )
