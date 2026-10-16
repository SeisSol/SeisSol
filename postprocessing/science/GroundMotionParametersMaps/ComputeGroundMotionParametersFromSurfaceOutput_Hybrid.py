#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
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
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN,
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# @section DESCRIPTION
# aim:
# 1 Read time histories from SeisSol surface or volume output in either HDF5 or POSIX format
# 2 Compute ground motion parameters (PGA, PGV, PGD, SA(T), optionally CAV) using GMRotDxx, RotDxx, GMRotIxx, or RotIxx
# 3 Store output in SeisSol-compatible XDMF/HDF5 files readable by ParaView
# Author(s): Thomas Ulrich, Fabian Kutschera
# Date(s): 29.09.17, 14.10.26

import sys
import os
import argparse
import time
from multiprocessing import Pool, cpu_count, Manager

import numpy as np
import seissolxdmf
import seissolxdmfwriter as sxw
from scipy import signal

try:
    from scipy.integrate import cumulative_trapezoid
except ImportError:
    from scipy.integrate import cumtrapz as cumulative_trapezoid

sys.path.append("%s/gmpe-smtk/" % os.path.dirname(sys.argv[0]))


def chunk(xs, n):
    """Split list/array xs into n evenly sized, non-empty chunks."""
    length = len(xs)
    if not 0 < n <= length:
        raise ValueError("The number of chunks must be between 1 and len(xs).")
    chunk_size, remainder = divmod(length, n)
    larger_size = chunk_size + 1
    return (
        [xs[p:p + larger_size] for p in range(0, remainder * larger_size, larger_size)]
        + [xs[p:p + chunk_size] for p in range(remainder * larger_size, length, chunk_size)]
    )


def validate_percentile(percentile):
    """Raise ValueError unless percentile is finite and lies in [0, 100]."""
    if not np.isfinite(percentile) or not 0.0 <= percentile <= 100.0:
        raise ValueError("Percentile must be a finite number between 0 and 100.")


def low_pass_filter(waveform, fs, cutoff_freq):
    """Apply a two-pass, zero-phase, second-order low-pass filter."""
    order = 2
    b, a = signal.butter(order, cutoff_freq, "low", fs=fs)
    return signal.filtfilt(b, a, waveform)


def butter_bandpass_filter(waveform, fs, freqmin, freqmax, order=2):
    """Apply a two-pass, zero-phase Butterworth band-pass filter."""
    nyquist = 0.5 * fs
    low = freqmin / nyquist
    high = freqmax / nyquist
    if low <= 0.0 or low >= high or high >= 1.0:
        raise ValueError("Band-pass frequencies must satisfy 0 < fmin < fmax < Nyquist.")
    b, a = signal.butter(order, [low, high], btype="band")
    return signal.filtfilt(b, a, waveform)


def compute_cav_gmrot(
    acceleration_x, time_step_x, acceleration_y, time_step_y,
    angles, percentile, idx_angle=None
):
    """Calculate GMRot CAV at a percentile or at one specified angle.

    The angles are in degrees, matching smtk.intensity_measures.rotate_horizontal.
    """
    from smtk.intensity_measures import get_cav, rotate_horizontal

    validate_percentile(percentile)
    if idx_angle is not None:
        if idx_angle == 0:
            return np.sqrt(
                get_cav(acceleration_x, time_step_x)
                * get_cav(acceleration_y, time_step_y)
            )
        theta = angles[idx_angle]
        rot_x, rot_y = rotate_horizontal(acceleration_x, acceleration_y, theta)
        return np.sqrt(get_cav(rot_x, time_step_x) * get_cav(rot_y, time_step_y))

    cav_theta = np.zeros(len(angles), dtype=float)
    for iloc, theta in enumerate(angles):
        if iloc == 0:
            cav_theta[iloc] = np.sqrt(
                get_cav(acceleration_x, time_step_x)
                * get_cav(acceleration_y, time_step_y)
            )
        else:
            rot_x, rot_y = rotate_horizontal(acceleration_x, acceleration_y, theta)
            cav_theta[iloc] = np.sqrt(
                get_cav(rot_x, time_step_x) * get_cav(rot_y, time_step_y)
            )
    return np.percentile(cav_theta, percentile)


def compute_cav_rotd(
    acceleration_x, time_step_x, acceleration_y, time_step_y,
    angles, percentile, idx_angle=None
):
    """Calculate RotD CAV at a percentile or at one specified angle.

    The angles are in radians, as used in the RotD calculation.
    """
    from smtk.intensity_measures import get_cav

    validate_percentile(percentile)
    if idx_angle is not None:
        rad = angles[idx_angle]
        rotated = acceleration_x * np.cos(rad) + acceleration_y * np.sin(rad)
        return get_cav(rotated, time_step_x)

    cav_theta = np.zeros(len(angles), dtype=float)
    for iloc, rad in enumerate(angles):
        rotated = acceleration_x * np.cos(rad) + acceleration_y * np.sin(rad)
        cav_theta[iloc] = get_cav(rotated, time_step_x)
    return np.percentile(cav_theta, percentile)


def _rotation_penalty_index(max_a_theta, percentile):
    """Select the RotI/GMRotI angle closest to the requested target spectrum.

    The first three columns of max_a_theta are PGA, PGV and PGD; remaining
    columns are spectral acceleration ordinates. The penalty is the mean
    squared relative difference between each angle's spectrum and the target
    percentile spectrum. Zero target ordinates are excluded because their
    relative error is undefined.
    """
    target = np.percentile(max_a_theta, percentile, axis=0)
    spectra_by_angle = max_a_theta[:, 3:]
    target_spectrum = target[3:]

    if not np.all(np.isfinite(target_spectrum)):
        raise ValueError("The target response spectrum contains non-finite values.")
    valid_periods = target_spectrum > 0.0
    if not np.any(valid_periods):
        raise ValueError(
            "Cannot determine RotI angle: the target response spectrum is zero "
            "at every requested period."
        )

    relative_difference = (
        spectra_by_angle[:, valid_periods] - target_spectrum[valid_periods]
    ) / target_spectrum[valid_periods]
    penalty = np.mean(relative_difference ** 2, axis=1)
    if not np.all(np.isfinite(penalty)):
        raise ValueError("The RotI angle penalty contains non-finite values.")
    return int(np.argmin(penalty)), target


def gmrotdpp_withPG(
    acceleration_x, time_step_x, acceleration_y, time_step_y, periods,
    percentile, damping=0.05, units="cm/s/s", method="Nigam-Jennings",
    ipp=False, include_cav=False
):
    """Compute GMRotD/GMRotI using the legacy ground-motion calculation.

    For backward compatibility, acceleration is truncated to all but the last
    sample for PGA, then integrated trapezoidally to obtain the velocity and
    displacement histories used for PGV and PGD. Those histories are assembled
    with the SMTK oscillator-acceleration histories before rotation, exactly as
    in the original script. The acceleration keeps the original script's
    numerical scaling; do not convert the SeisSol velocity input to cm/s.
    """
    from smtk.intensity_measures import (
        get_response_spectrum, equalise_series, rotate_horizontal
    )

    validate_percentile(percentile)
    if len(acceleration_x) != len(acceleration_y):
        raise ValueError("The two horizontal acceleration histories must have equal lengths.")
    if not np.isclose(time_step_x, time_step_y, rtol=0.0, atol=1.0e-10):
        raise ValueError("The two horizontal components must have the same time step.")

    _, _, response_acc_x, _, _ = get_response_spectrum(
        acceleration_x, time_step_x, periods, damping, units, method
    )
    _, _, response_acc_y, _, _ = get_response_spectrum(
        acceleration_y, time_step_y, periods, damping, units, method
    )
    response_acc_x, response_acc_y = equalise_series(response_acc_x, response_acc_y)

    # Preserve the original algorithm and array lengths: acceleration[0:-1]
    # is integrated twice, and each series has the same length as the SMTK
    # Nigam-Jennings oscillator histories (typically n_samples - 1).
    velocity_x = time_step_x * cumulative_trapezoid(acceleration_x[0:-1], initial=0.0)
    displacement_x = time_step_x * cumulative_trapezoid(velocity_x, initial=0.0)
    x_a = np.column_stack((acceleration_x[0:-1], velocity_x, displacement_x, response_acc_x))

    velocity_y = time_step_y * cumulative_trapezoid(acceleration_y[0:-1], initial=0.0)
    displacement_y = time_step_y * cumulative_trapezoid(velocity_y, initial=0.0)
    y_a = np.column_stack((acceleration_y[0:-1], velocity_y, displacement_y, response_acc_y))

    angles = np.arange(0.0, 90.0, 1.0)  # degrees; 90 degrees is redundant for GMRot
    max_a_theta = np.zeros((len(angles), len(periods) + 3), dtype=float)

    for iloc, theta in enumerate(angles):
        if iloc == 0:
            max_a_theta[iloc, :] = np.sqrt(
                np.max(np.abs(x_a), axis=0) * np.max(np.abs(y_a), axis=0)
            )
        else:
            rot_x, rot_y = rotate_horizontal(x_a, y_a, theta)
            max_a_theta[iloc, :] = np.sqrt(
                np.max(np.abs(rot_x), axis=0) * np.max(np.abs(rot_y), axis=0)
            )

    if not ipp:
        gmrotd = np.percentile(max_a_theta, percentile, axis=0)
        result = {
            "PGA": gmrotd[0],
            "PGV": gmrotd[1],
            "PGD": gmrotd[2],
            "Acceleration": gmrotd[3:],
        }
        if include_cav:
            result["CAV"] = compute_cav_gmrot(
                acceleration_x, time_step_x, acceleration_y, time_step_y,
                angles, percentile
            )
    else:
        idx_i_perc, _ = _rotation_penalty_index(max_a_theta, percentile)
        gmroti = max_a_theta[idx_i_perc, :]
        result = {
            "PGA": gmroti[0],
            "PGV": gmroti[1],
            "PGD": gmroti[2],
            "Acceleration": gmroti[3:],
        }
        if include_cav:
            result["CAV"] = compute_cav_gmrot(
                acceleration_x, time_step_x, acceleration_y, time_step_y,
                angles, percentile, idx_angle=idx_i_perc
            )

    return result


def rotdpp_withPG(
    acceleration_x, time_step_x, acceleration_y, time_step_y, periods,
    percentile, damping=0.05, units="cm/s/s", method="Nigam-Jennings",
    ipp=False, include_cav=False
):
    """Compute RotD/RotI using the legacy ground-motion calculation.

    The acceleration histories are truncated to all but the last sample for
    PGA and integrated twice for PGV and PGD, then combined with the oscillator
    acceleration histories before rotation, preserving the original script's
    calculation and array lengths. The acceleration keeps the original script's
    numerical scaling; do not convert the SeisSol velocity input to cm/s.
    """
    from smtk.intensity_measures import get_response_spectrum, equalise_series

    validate_percentile(percentile)
    if len(acceleration_x) != len(acceleration_y):
        raise ValueError("The two horizontal acceleration histories must have equal lengths.")
    if not np.isclose(time_step_x, time_step_y, rtol=0.0, atol=1.0e-10):
        raise ValueError("The two horizontal components must have the same time step.")

    _, _, response_acc_x, _, _ = get_response_spectrum(
        acceleration_x, time_step_x, periods, damping, units, method
    )
    _, _, response_acc_y, _, _ = get_response_spectrum(
        acceleration_y, time_step_y, periods, damping, units, method
    )
    response_acc_x, response_acc_y = equalise_series(response_acc_x, response_acc_y)

    velocity_x = time_step_x * cumulative_trapezoid(acceleration_x[0:-1], initial=0.0)
    displacement_x = time_step_x * cumulative_trapezoid(velocity_x, initial=0.0)
    x_a = np.column_stack((acceleration_x[0:-1], velocity_x, displacement_x, response_acc_x))

    velocity_y = time_step_y * cumulative_trapezoid(acceleration_y[0:-1], initial=0.0)
    displacement_y = time_step_y * cumulative_trapezoid(velocity_y, initial=0.0)
    y_a = np.column_stack((acceleration_y[0:-1], velocity_y, displacement_y, response_acc_y))

    angles = np.radians(np.arange(0.0, 180.0, 1.0))
    max_a_theta = np.zeros((len(angles), len(periods) + 3), dtype=float)

    for iloc, rad in enumerate(angles):
        rotated = x_a * np.cos(rad) + y_a * np.sin(rad)
        max_a_theta[iloc, :] = np.max(np.abs(rotated), axis=0)

    if not ipp:
        rotd = np.percentile(max_a_theta, percentile, axis=0)
        result = {
            "PGA": rotd[0],
            "PGV": rotd[1],
            "PGD": rotd[2],
            "Acceleration": rotd[3:],
        }
        if include_cav:
            result["CAV"] = compute_cav_rotd(
                acceleration_x, time_step_x, acceleration_y, time_step_y,
                angles, percentile
            )
    else:
        idx_i_perc, _ = _rotation_penalty_index(max_a_theta, percentile)
        roti = max_a_theta[idx_i_perc, :]
        result = {
            "PGA": roti[0],
            "PGV": roti[1],
            "PGD": roti[2],
            "Acceleration": roti[3:],
        }
        if include_cav:
            result["CAV"] = compute_cav_rotd(
                acceleration_x, time_step_x, acceleration_y, time_step_y,
                angles, percentile, idx_angle=idx_i_perc
            )

    return result


def ComputeGroundMotionParameters(args2):
    ir, q = args2

    # Preserve the original script's output scaling: use SeisSol velocities
    # directly, without converting m/s to cm/s before differentiation.
    velocity_x = u[ir, :].copy()
    if args.lowpass:
        velocity_x = low_pass_filter(
            velocity_x, fs=1.0 / dt, cutoff_freq=args.lowpass[0]
        )
    elif args.bandpass:
        velocity_x = butter_bandpass_filter(
            velocity_x, fs=1.0 / dt,
            freqmin=args.bandpass[0], freqmax=args.bandpass[1]
        )
    acceleration_x = np.gradient(velocity_x, dt)

    velocity_y = v[ir, :].copy()
    if args.lowpass:
        velocity_y = low_pass_filter(
            velocity_y, fs=1.0 / dt, cutoff_freq=args.lowpass[0]
        )
    elif args.bandpass:
        velocity_y = butter_bandpass_filter(
            velocity_y, fs=1.0 / dt,
            freqmin=args.bandpass[0], freqmax=args.bandpass[1]
        )
    acceleration_y = np.gradient(velocity_y, dt)

    # Keep the original call pattern: PGA comes from differentiated velocity;
    # PGV and PGD are obtained by integrating that acceleration in the metric
    # routines, not from the original SeisSol velocity field directly.
    result = myfunc(
        acceleration_x, dt, acceleration_y, dt, periods,
        percentile=percentile, damping=0.05, units="cm/s/s",
        method="Nigam-Jennings", ipp=args.ipp, include_cav=args.CAV
    )

    row = [result["PGA"], result["PGV"], result["PGD"]]
    if args.CAV:
        row.append(result["CAV"])
    row.extend(result["Acceleration"])

    if q != 0:
        q.put(ir)
    return row


start = time.time()

# Command-line arguments
parser = argparse.ArgumentParser(
    description=(
        "Compute ground-motion parameter maps (PGA, PGV, PGD, SA(T), optionally CAV) "
        "from SeisSol free-surface or volume output."
    )
)
parser.add_argument("filename", help="Surface/volume output XDMF file")
parser.add_argument(
    "--MP", nargs=1, metavar=("nthreads"), default=[1], type=int,
    help="Number of multiprocessing workers for each MPI rank (default: 1)"
)
parser.add_argument(
    "--noMPI", dest="noMPI", action="store_true", default=False,
    help="Run without MPI (single rank)"
)
parser.add_argument(
    "--periods", nargs="+", type=float,
    help="Periods in seconds for computing spectral acceleration"
)
parser.add_argument(
    "--ipp", dest="ipp", action="store_true", default=False,
    help="Use period-independent RotI/GMRotI instead of period-dependent RotD/GMRotD"
)
parser.add_argument(
    "--CAV", dest="CAV", action="store_true", default=False,
    help="Also compute cumulative absolute velocity"
)
parser.add_argument(
    "--metric", choices=["GMRot", "Rot"], default="GMRot",
    help=(
        "Rotation metric family: GMRot takes the geometric mean of the peaks of "
        "the two rotated components; Rot uses one projected horizontal component. "
        "Without --ipp, these produce GMRotD/RotD; with --ipp, GMRotI/RotI. "
        "Default: GMRot."
    )
)
parser.add_argument(
    "--percentile", type=float, default=50.0,
    help="Percentile across rotation angles (0 to 100; default: 50)"
)
filter_group = parser.add_mutually_exclusive_group()
filter_group.add_argument(
    "--lowpass", nargs=1, metavar=("cutoff_freq"), type=float,
    help="Low-pass filter velocity waveforms before computing ground-motion measures; cutoff frequency in Hz"
)
filter_group.add_argument(
    "--bandpass", nargs=2, metavar=("fmin", "fmax"), type=float,
    help="Band-pass filter velocity waveforms before computing ground-motion measures; fmin and fmax in Hz"
)

args = parser.parse_args()

if not np.isfinite(args.percentile) or not 0.0 <= args.percentile <= 100.0:
    parser.error("--percentile must be a finite number between 0 and 100.")
if args.MP[0] < 1:
    parser.error("--MP must be at least 1.")
if args.periods is not None and (
    not args.periods or not np.all(np.isfinite(args.periods)) or np.any(np.asarray(args.periods) <= 0.0)
):
    parser.error("--periods must contain finite, positive period values.")
if args.lowpass and (not np.isfinite(args.lowpass[0]) or args.lowpass[0] <= 0.0):
    parser.error("--lowpass cutoff frequency must be finite and positive.")
if args.bandpass and (
    not np.all(np.isfinite(args.bandpass))
    or args.bandpass[0] <= 0.0
    or args.bandpass[1] <= args.bandpass[0]
):
    parser.error("--bandpass requires finite frequencies satisfying 0 < fmin < fmax.")

if args.noMPI:
    nranks = 1
    irank = 0
else:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    nranks = comm.size
    irank = comm.rank

sx = seissolxdmf.seissolxdmf(args.filename)
dt = sx.ReadTimeStep()
if args.lowpass and args.lowpass[0] >= 0.5 / dt:
    parser.error("--lowpass cutoff frequency must be below the Nyquist frequency (0.5 / dt).")
if args.bandpass and args.bandpass[1] >= 0.5 / dt:
    parser.error("--bandpass fmax must be below the Nyquist frequency (0.5 / dt).")
nElements = sx.nElements
connect = sx.ReadConnect()

isVolumeData = connect.shape[1] == 4
prefixType = "-volume" if isVolumeData else "-surface"

# Split element indices across MPI ranks.
inputs0 = np.arange(0, nElements)
inputsi = chunk(inputs0, nranks)
inputs = inputsi[irank]
nElements_i = len(inputs)

# Read horizontal velocity data for this rank.
available_fields = sx.ReadAvailableDataFields()
componentVariables = ["u", "v", "w"] if "u" in available_fields else ["v1", "v2", "v3"]
u = sx.ReadDataChunk(
    componentVariables[0], firstElement=inputsi[irank][0], nchunk=nElements_i, idt=-1
)
v = sx.ReadDataChunk(
    componentVariables[1], firstElement=inputsi[irank][0], nchunk=nElements_i, idt=-1
)
u = np.transpose(u)
v = np.transpose(v)

if irank == 0:
    print(f"Done reading {prefixType[1:]} data: {time.time() - start:.2f}s")

if args.periods is None:
    periods = np.array([0.100, 0.125, 0.25, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 5.0])
    if irank == 0:
        print("No periods specified: using default periods", periods)
else:
    periods = np.asarray(args.periods, dtype=float)
    if irank == 0:
        print("Periods", periods)

percentile = args.percentile
myfunc = rotdpp_withPG if args.metric == "Rot" else gmrotdpp_withPG

nTasksPerRank = args.MP[0]
if nTasksPerRank > cpu_count():
    parser.error(f"--MP ({nTasksPerRank}) cannot exceed available CPUs ({cpu_count()}).")
if nElements_i < 1:
    raise RuntimeError("This MPI rank has no elements to process; reduce the number of MPI ranks.")

pool = Pool(processes=nTasksPerRank)
manager = Manager()
q = manager.Queue()
args2 = [(i - inputs[0], q) for i in inputs]
result_async = pool.map_async(ComputeGroundMotionParameters, args2)
pool.close()

while not result_async.ready():
    remaining = nElements_i - q.qsize()
    print(irank, ": Waiting for", remaining, "tasks to complete...")
    time.sleep(10.0)

myres = np.ascontiguousarray(np.asarray(result_async.get()))
pool.join()
manager.shutdown()

if irank == 0:
    print("Done computing ground-motion estimates: %f" % (time.time() - start))

dataName = ["PGA", "PGV", "PGD"]
if args.CAV:
    dataName.append("CAV")
for per in periods:
    dataName.append("SA%06.3fs" % per)

mypath, fn = os.path.split(args.filename)
prefix = fn.split("-")[-2] if not isVolumeData else fn.split(".")[-2]
perc_str = f"{int(percentile)}" if percentile.is_integer() else f"{percentile:g}"
base_metric = args.metric
type_letter = "I" if args.ipp else "D"
effective_metric = f"{base_metric}{type_letter}{perc_str}"

filter_suffix = ""
if args.lowpass:
    filter_suffix = f"_lp{args.lowpass[0]:.1f}"
elif args.bandpass:
    filter_suffix = f"_bp{args.bandpass[0]:.2f}-{args.bandpass[1]:.2f}"

prefixGME = f"{prefix}_{effective_metric}{filter_suffix}-GME"

if not args.noMPI:
    comm.Barrier()
    myres = comm.gather(myres, root=0)

if irank == 0:
    if nranks > 1:
        myres = np.concatenate(myres, axis=0)
    geom = sx.ReadGeometry()
    sxw.write_seissol_output(
        f"{prefixGME}{prefixType}", geom, connect, dataName,
        [myres[:, i] for i in range(myres.shape[1])], 0, [0]
    )
