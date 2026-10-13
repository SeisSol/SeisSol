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
#2 compute ground motion parameters (PGA, PGV, PGD, SA(T), CAV) using GMRotD50, RotD50, GMRotI50, or RotI50
#3 store output in a hdf5 file readable by paraview

import sys
import os
import h5py
import numpy as np
import argparse
from multiprocessing import Pool, cpu_count, Manager
import time
import seissolxdmf
import seissolxdmfwriter as sxw
from scipy import signal

try:
    from scipy.integrate import cumulative_trapezoid
except ImportError:
    from scipy.integrate import cumtrapz as cumulative_trapezoid

sys.path.append("%s/gmpe-smtk/" % (os.path.dirname(sys.argv[0])))


def chunk(xs, n):
    '''Split the list, xs, into n evenly sized chunks'''
    L = len(xs)
    assert 0 < n <= L
    s, r = divmod(L, n)
    t = s + 1
    return ([xs[p:p+t] for p in range(0, r*t, t)] +
            [xs[p:p+s] for p in range(r*t, L, s)])


def low_pass_filter(waveform, fs, cutoff_freq):
    "Applied 2-pass zero-phase order 2 low pass filter on data"
    order = 2
    b, a = signal.butter(order, cutoff_freq, 'low', fs=fs)
    return signal.filtfilt(b, a, waveform)


def butter_bandpass_filter(waveform, fs, freqmin, freqmax, order=2):
    "Applied 2-pass zero-phase order 2 bandpass filter on data"
    nyq = 0.5 * fs
    low = freqmin / nyq
    high = freqmax / nyq
    if high >= 1.0:
        high = 0.99
    b, a = signal.butter(order, [low, high], btype='band')
    return signal.filtfilt(b, a, waveform)


def compute_cav_gmrot(acceleration_x, time_step_x, acceleration_y, time_step_y, angles, percentile=50, idx_angle=None):
    """Compute CAV using GMRot (Boore et al., 2006) for percentile or single optimal angle"""
    from smtk.intensity_measures import get_cav, rotate_horizontal
    if idx_angle is not None:
        theta = angles[idx_angle]
        if idx_angle == 0:
            return np.sqrt(get_cav(acceleration_x, time_step_x) * get_cav(acceleration_y, time_step_y))
        else:
            rot_x, rot_y = rotate_horizontal(acceleration_x, acceleration_y, theta)
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


def compute_cav_rotd(acceleration_x, time_step_x, acceleration_y, time_step_y, angles, percentile=50, idx_angle=None):
    """Compute CAV using single vector RotD (Boore, 2010) for percentile or single optimal angle"""
    from smtk.intensity_measures import get_cav
    if idx_angle is not None:
        rad = angles[idx_angle]
        rot_single = acceleration_x * np.cos(rad) + acceleration_y * np.sin(rad)
        return get_cav(rot_single, time_step_x)

    cav_theta = np.zeros(len(angles), dtype=float)
    for iloc, rad in enumerate(angles):
        rot_single = acceleration_x * np.cos(rad) + acceleration_y * np.sin(rad)
        cav_theta[iloc] = get_cav(rot_single, time_step_x)
    return np.percentile(cav_theta, percentile)


def gmrotdpp_withPG(acceleration_x, time_step_x, acceleration_y, time_step_y, periods,
                    percentile=50, damping=0.05, units="cm/s/s", method="Nigam-Jennings", ipp=False):
    """
    Computes GMRotD50 or GMRotI50 ground motion parameters (Boore et al., 2006).
    """
    from smtk.intensity_measures import get_response_spectrum, equalise_series, rotate_horizontal
    if (percentile > 100. + 1E-9) or (percentile < 0.):
        raise ValueError("Percentile for GMRotDpp/GMRotIpp must be between 0. and 100.")

    sax, _, x_a, _, _ = get_response_spectrum(acceleration_x, time_step_x, periods, damping, units, method)
    say, _, y_a, _, _ = get_response_spectrum(acceleration_y, time_step_y, periods, damping, units, method)
    x_a, y_a = equalise_series(x_a, y_a)

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
            max_a_theta[iloc, :] = np.sqrt(np.max(np.fabs(x_a), axis=0) * np.max(np.fabs(y_a), axis=0))
        else:
            rot_x, rot_y = rotate_horizontal(x_a, y_a, theta)
            max_a_theta[iloc, :] = np.sqrt(np.max(np.fabs(rot_x), axis=0) * np.max(np.fabs(rot_y), axis=0))

    if not ipp:
        gmrotd = np.percentile(max_a_theta, percentile, axis=0)
        res = {"PGA": gmrotd[0],
               "PGV": gmrotd[1],
               "PGD": gmrotd[2],
               "Acceleration": gmrotd[3:]}
        if args.CAV:
            res['CAV'] = compute_cav_gmrot(acceleration_x, time_step_x, acceleration_y, time_step_y, angles, percentile)
    else:
        # GMRotI50 implementation (Boore et al., 2006)
        d50 = np.percentile(max_a_theta, percentile, axis=0)
        sa_theta = max_a_theta[:, 3:]
        sa_d50 = d50[3:]
        penalty = np.mean(((sa_theta - sa_d50) / sa_d50) ** 2, axis=1)
        idx_i50 = np.argmin(penalty)
        gmroti = max_a_theta[idx_i50, :]

        res = {"PGA": gmroti[0],
               "PGV": gmroti[1],
               "PGD": gmroti[2],
               "Acceleration": gmroti[3:]}
        if args.CAV:
            res['CAV'] = compute_cav_gmrot(acceleration_x, time_step_x, acceleration_y, time_step_y, angles, idx_angle=idx_i50)

    return res


def rotdpp_withPG(acceleration_x, time_step_x, acceleration_y, time_step_y, periods,
                  percentile=50, damping=0.05, units="cm/s/s", method="Nigam-Jennings", ipp=False):
    """
    Computes RotD50 or RotI50 ground motion parameters (Boore, 2010).
    Rotates single 1D component over angles 0 to 180 degrees.
    """
    from smtk.intensity_measures import get_response_spectrum, equalise_series
    if (percentile > 100. + 1E-9) or (percentile < 0.):
        raise ValueError("Percentile for RotDpp/RotIpp must be between 0. and 100.")

    sax, _, x_a, _, _ = get_response_spectrum(acceleration_x, time_step_x, periods, damping, units, method)
    say, _, y_a, _, _ = get_response_spectrum(acceleration_y, time_step_y, periods, damping, units, method)
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
        rot_single = x_a * np.cos(rad) + y_a * np.sin(rad)
        max_a_theta[iloc, :] = np.max(np.abs(rot_single), axis=0)

    if not ipp:
        rotd = np.percentile(max_a_theta, percentile, axis=0)
        res = {"PGA": rotd[0],
               "PGV": rotd[1],
               "PGD": rotd[2],
               "Acceleration": rotd[3:]}
        if args.CAV:
            res['CAV'] = compute_cav_rotd(acceleration_x, time_step_x, acceleration_y, time_step_y, angles, percentile)
    else:
        # RotI50 implementation (Boore, 2010)
        d50 = np.percentile(max_a_theta, percentile, axis=0)
        sa_theta = max_a_theta[:, 3:]
        sa_d50 = d50[3:]
        penalty = np.mean(((sa_theta - sa_d50) / sa_d50) ** 2, axis=1)
        idx_i50 = np.argmin(penalty)
        roti = max_a_theta[idx_i50, :]

        res = {"PGA": roti[0],
               "PGV": roti[1],
               "PGD": roti[2],
               "Acceleration": roti[3:]}
        if args.CAV:
            res['CAV'] = compute_cav_rotd(acceleration_x, time_step_x, acceleration_y, time_step_y, angles, idx_angle=idx_i50)

    return res


def ComputeGroundMotionParameters(args2):
    ir, q = args2

    # Xdirection
    velocity_x = u[ir, :]
    if args.lowpass:
        velocity_x = low_pass_filter(velocity_x, fs=1./dt, cutoff_freq=args.lowpass[0])
    elif args.bandpass:
        velocity_x = butter_bandpass_filter(velocity_x, fs=1./dt, freqmin=args.bandpass[0], freqmax=args.bandpass[1])
    acceleration_x = np.gradient(velocity_x, dt)

    # Ydirection
    velocity_y = v[ir, :]
    if args.lowpass:
        velocity_y = low_pass_filter(velocity_y, fs=1./dt, cutoff_freq=args.lowpass[0])
    elif args.bandpass:
        velocity_y = butter_bandpass_filter(velocity_y, fs=1./dt, freqmin=args.bandpass[0], freqmax=args.bandpass[1])
    acceleration_y = np.gradient(velocity_y, dt)

    # Compute SA(T), PGA, PGV, PGD, CAV...
    res = myfunc(acceleration_x, dt, acceleration_y, dt, periods,
                 percentile=50, damping=0.05, units="cm/s/s", method="Nigam-Jennings", ipp=args.ipp)

    myres = [res['PGA'], res['PGV'], res['PGD']]
    if args.CAV:
        myres.append(res['CAV'])
    myres += list(res['Acceleration'])

    if q != 0:
        q.put(ir)
    return myres


start = time.time()

# Parsing python arguments
parser = argparse.ArgumentParser(description='Compute ground motion parameters maps (e.g. PGA, PGV, PGD, SA) from free surface output.')
parser.add_argument('filename', help='Surface output (xdmf file)')
parser.add_argument('--MP', nargs=1, metavar=('nthreads'), default=([1]), help='Use multiprocessing Pool to speed up calculations', type=int)
parser.add_argument('--noMPI', dest='noMPI', action='store_true', default=False, help='Run on one node (skip mpi4py)')
parser.add_argument('--periods', nargs='+', help='List of periods for computing SA(T)', type=float)
parser.add_argument('--ipp', dest='ipp', action='store_true', default=False,
                    help='Use period-independent estimate (GMRotI50 or RotI50) based on Boore et al. (2006) / Boore (2010)')
parser.add_argument('--CAV', dest='CAV', action='store_true', default=False, help='Compute cumulated absolute velocity')
parser.add_argument('--metric', choices=['GMRotD50', 'RotD50'], default='GMRotD50',
                    help='Horizontal rotation metric type (GMRotD50 or RotD50). Default: GMRotD50')

# Filter arguments (mutually exclusive)
filter_group = parser.add_mutually_exclusive_group()
filter_group.add_argument('--lowpass', nargs=1, metavar=('cutoff_freq'), type=float,
                          help='Low-pass filter velocity waveforms prior to computing ground motion estimates. Cutoff frequency in Hz.')
filter_group.add_argument('--bandpass', nargs=2, metavar=('fmin', 'fmax'), type=float,
                          help='Band-pass filter velocity waveforms prior to computing ground motion estimates. fmin and fmax in Hz.')

args = parser.parse_args()

if args.noMPI:
    nranks = 1
    irank = 0
else:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    nranks = comm.size
    irank = comm.rank

if args.metric == 'RotD50':
    myfunc = rotdpp_withPG
else:
    myfunc = gmrotdpp_withPG

sx = seissolxdmf.seissolxdmf(args.filename)
dt = sx.ReadTimeStep()
nElements = sx.nElements
connect = sx.ReadConnect()

isVolumeData = True if connect.shape[1] == 4 else False
prefixType = "-volume" if isVolumeData else "-surface"

# Split input array across MPI ranks
inputs0 = np.arange(0, nElements)
inputsi = chunk(inputs0, nranks)
inputs = inputsi[irank]
nElements_i = len(inputs)

# Read horizontal velocity data chunks required by rank
componentVariables = ['u', 'v', 'w'] if 'u' in sx.ReadAvailableDataFields() else ['v1', 'v2', 'v3']
u = sx.ReadDataChunk(componentVariables[0], firstElement=inputsi[irank][0], nchunk=nElements_i, idt=-1)
v = sx.ReadDataChunk(componentVariables[1], firstElement=inputsi[irank][0], nchunk=nElements_i, idt=-1)

u = np.transpose(u)
v = np.transpose(v)
ndt = u.shape[1]

if irank == 0:
    print(f"Done reading {prefixType[1:]} data: {time.time()-start:.2f}s")

# Parameters for response spectrum
damping = 0.05
units = "cm/s/s"

if args.periods is None:
    periods = np.array([0.100, 0.125, 0.25, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 5.0])
    if irank == 0:
        print("No periods specified: using default periods", periods)
else:
    periods = np.array([float(per) for per in args.periods])  
    if irank == 0:
        print("Periods", periods)

nTasksPerRank = args.MP[0]
assert nTasksPerRank <= cpu_count()
pool = Pool(processes=nTasksPerRank)
m = Manager()
q = m.Queue()

args2 = [(i - inputs[0], q) for i in inputs]
Result = pool.map_async(ComputeGroundMotionParameters, args2)
pool.close()

while True:
    if Result.ready():
        break
    remaining = nElements_i + 1 - q.qsize()
    print(irank, ": Waiting for", remaining, "tasks to complete...")
    time.sleep(10.0)

myres = np.ascontiguousarray(np.array(Result.get()))

if irank == 0:
    print("Done computing ground motion estimates: %f" % (time.time() - start))

dataName = ['PGA', 'PGV', 'PGD']

if args.CAV:
    dataName.append('CAV')

for per in periods:
    dataName.append('SA%06.3fs' % per)

mypath, fn = os.path.split(args.filename)
prefix = fn.split('-')[-2] if not isVolumeData else fn.split('.')[-2]

effective_metric = f"{args.metric[:-3]}I50" if args.ipp else args.metric

sFilter = ''
if args.lowpass:
    sFilter = f'_lp{args.lowpass[0]:.1f}'
elif args.bandpass:
    sFilter = f'_bp{args.bandpass[0]:.2f}-{args.bandpass[1]:.2f}'

prefixGME = f'{prefix}_{effective_metric}{sFilter}-GME'

if not args.noMPI:
    comm.Barrier()
    myres = comm.gather(myres, root=0)

if irank == 0:
    if nranks > 1:
        myres = np.concatenate(myres, axis=0)
    geom = sx.ReadGeometry()
    sxw.write_seissol_output(f"{prefixGME}{prefixType}", geom, connect, dataName,
                             [myres[:, i] for i in range(myres.shape[1])], 0, [0])
