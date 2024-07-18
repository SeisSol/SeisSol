# Script to determine the rise time (time the absolute slip rate is above a certain threshold) from seissol fault.xdmf output

import numpy as np
import seissolxdmf as sx
import seissolxdmfwriter as sxw
import argparse
import os
import timeit
from numba import jit
from multiprocessing import Pool, cpu_count

parser = argparse.ArgumentParser(
    description="calculate rise times from fault xdmf file"
)
parser.add_argument("filename", help="path+filename-fault.xdmf")
parser.add_argument(
    "--threshold",
    nargs=1,
    metavar=("threshold"),
    default=([0.1]),
    help="specify slip rate threshold",
    type=float,
)
parser.add_argument(
    "--events",
    nargs=1,
    metavar=("events"),
    choices=[1, 2],
    default=([1]),
    help="number of events [1,2] in xdmf file",
    type=int,
)
parser.add_argument(
    "--median_filter_radius",
    nargs=1,
    default=([0.0]),
    metavar=("radius of the median filter in m"),
    help="the rise time value of a fault cell is set to the median value of all surrounding fault cells within \
                    a sphere with this radius; 0 turns the filter off",
    type=float,
)
parser.add_argument(
    "--nprocs",
    nargs=1,
    metavar=("number of processes"),
    default=([4]),
    help="use multiprocessing to speed up the median filtering",
    type=int,
)
parser.add_argument(
    "--complicated_fault",
    dest="complicated_fault",
    action="store_true",
    default=False,
    help="Improves median filter for branching faults and removes isolated slipping elements",
)
args = parser.parse_args()


@jit(nopython=True)
def counting_loop(row, array, threshold=0.1):
    """counting the time steps for which a certain point on the fault has an ASR above the threshold"""
    count = 0
    for j in range(array.shape[1]):
        if array[row, j] < threshold and count > 0:
            break
        if array[row, j] < threshold and count == 0:
            continue
        if array[row, j] >= threshold:
            count += 1
    return count


@jit(nopython=True)
def get_rise_time(ASR_array, threshold=0.1, dt=0.1, events=1):

    RT = np.zeros_like(ASR_array[0:2, :])
    mid = int(ASR_array.shape[0] / 2)
    ASR_array = np.transpose(ASR_array)  # to speed up

    if events == 1:
        for i in range(ASR_array.shape[0]):
            count = counting_loop(i, ASR_array, threshold=threshold)
            RT[0, i] = count * dt
        return RT

    if events == 2:
        for i in range(ASR_array.shape[0]):
            count = counting_loop(i, ASR_array[:, 0:mid], threshold=threshold)
            RT[0, i] = count * dt
            count = counting_loop(i, ASR_array[:, mid:], threshold=threshold)
            RT[1, i] = count * dt
        return RT


def get_absolute_slip_rate(sx_object):
    SRs = sx_object.ReadData("SRs")
    SRd = sx_object.ReadData("SRd")
    ASR = np.sqrt(SRs**2 + SRd**2)
    return ASR


@jit(nopython=True)
def MedianSmoothing(input_tuple):
    # input_tuple: (variable, xyz, indices, radius)
    variable, xyz, indices, radius = input_tuple
    variableSmooth = np.zeros(indices.shape[0])
    for i, j in enumerate(indices):
        xyz0 = xyz[j, :]
        distances = np.sqrt(np.sum((xyz - xyz0) ** 2, axis=1))
        if complicated_fault:
            tmp_ind1 = distances <= radius
            tmp_ind2 = variable != 0
            tmp_ind3 = tmp_ind1 & tmp_ind2
            non_zero_portion = np.sum(tmp_ind3) / np.sum(tmp_ind1)
            if (
                non_zero_portion < 0.75 and variable[j] == 0
            ):  # improve filter for branching faults
                variableSmooth[i] = 0
            elif non_zero_portion < 0.35:  # remove isolated slipping elements
                variableSmooth[i] = 0
            else:
                variableSmooth[i] = np.median(variable[tmp_ind3])
        else:
            variableSmooth[i] = np.median(variable[distances <= radius])
    return variableSmooth


def MedianSmoothingParallel(variable, xyz, radius=250.0, nprocs=4):

    chunks = [
        (variable, xyz, i, radius)
        for i in np.array_split(np.arange(variable.shape[0]), nprocs)
    ]
    with Pool(nprocs) as pool:
        result = pool.map(MedianSmoothing, chunks)
        pool.close()
        pool.join()
    return np.concatenate(result)


def get_xyz_from_connect(geom, connect):
    # Genrate array with coordinates of triangle midpoints (same dimension as connect)
    xyz = np.zeros_like(connect, dtype=float)
    xyz = (1.0 / 3.0) * (
        geom[connect[:, 0], :] + geom[connect[:, 1], :] + geom[connect[:, 2], :]
    )
    return xyz


print("Loading file and calculating absolute slip rates...")
sx = sx.seissolxdmf(args.filename)
dt = sx.ReadTimeStep()
ASR = get_absolute_slip_rate(sx)

print("Determining rise time...")
start = timeit.default_timer()
RT = get_rise_time(ASR, threshold=args.threshold[0], dt=dt, events=args.events[0])
stop = timeit.default_timer()
print(f"Time to compute rise times: {stop - start:.2f}")

if args.median_filter_radius[0]:

    connect = sx.ReadConnect()
    geom = sx.ReadGeometry()
    xyz = get_xyz_from_connect(geom, connect)
    complicated_fault = args.complicated_fault

    print(f"Median filter radius: {args.median_filter_radius[0]}")
    assert args.nprocs[0] >= 1 and args.nprocs[0] <= cpu_count()
    print(f"Median smoothing using {args.nprocs[0]} processes...")

    start = timeit.default_timer()
    for i in range(args.events[0]):
        RT[i] = MedianSmoothingParallel(
            RT[i], xyz, radius=args.median_filter_radius[0], nprocs=args.nprocs[0]
        )

    stop = timeit.default_timer()
    print(f"Time to apply median filter: {stop - start:.2f}")


print("Writing output file...")
connect = sx.ReadConnect()
xyz = sx.ReadGeometry()
aDataName = ["risetime"] if args.events[0] == 1 else ["risetime1", "risetime2"]
fn = os.path.basename(args.filename).rsplit(".", 1)[0] + "_risetime"
lRT = [RT[0, :]] if args.events[0] == 1 else [RT[0, :], RT[1, :]]
sxw.write_seissol_output(fn, xyz, connect, aDataName, lRT, dt, [0])
