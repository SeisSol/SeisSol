# Script to compute isochrones (for S-wave radiation) of a single receiver from SeisSol's surface/fault.xdmf output
# Isochrones are calculated from RT output and peak slip rate times

import numpy as np
import seissolxdmf
import argparse
import seissolxdmfwriter as sxw
import timeit
from multiprocessing import Pool, cpu_count

parser = argparse.ArgumentParser(description='get isochrones from SeisSol\'s xdmf output')

parser.add_argument('filename', help='path+filename-surface.xdmf')

parser.add_argument("x", nargs=1, type=float, help="receiver x")
parser.add_argument("y", nargs=1, type=float, help="receiver y")
parser.add_argument("--z", nargs=1, type=float, default=([0.]), help="receiver z")

parser.add_argument('--faultXdmf', type=str, metavar=('-fault.xdmf'), 
                    help='provide path+filename-fault.xdmf; only needed when the path differs from the -surface.xdmf file')

parser.add_argument('--output', choices=["numpy","xdmf","both"], default="xdmf", 
                    help='choose the output format')

parser.add_argument('--sliprateThreshold', nargs=1, default=([0.05]), metavar=('peak sliprate threshold in m/s'),
                    help='fault cells below the treshold are assumed to not radiate waves' , type=float)

parser.add_argument('--sliprateParameters', choices=["absolute","components","both"], default="both", 
                    help="absolute: compute isochrones from peak times of absolute slip rates;"+ 
                    " components: peak slip rate times are computed for SRs and SRd separately")

parser.add_argument('--medianFilterRadius', nargs=1, default=([0.]), metavar=('radius of the median filter in m'),
                    help='peak slip rate time of a fault cell is set to the median value of all surrounding fault cells within \
                    a sphere with this radius' , type=float)
# This prevents strong oscillations of peak slip rate times where different rupture fronts have similar peak slip rates 
# and picks the time of the dominant rupture front.
# Default = 0. turns the filter off
# Radii between 150 m and 300 m are advised (also depending on the resolution of the fault output)

parser.add_argument('--threads', nargs=1, metavar=('threads'), default=([4]),
                    help='use multiprocessing to speed up the median filtering' ,type=int)

parser.add_argument('--events', nargs=1, choices=[0,1,2], default=([0]), 
                    help='0: xdmf file contains just one event; 1: process the first of two events; 2: process the second of two events' ,type=int) 
# Needed for Ridgecrest output, where the xdmf file contains two events 

args = parser.parse_args()


def ComputeTriangleMidpoints(geom, connect):
    """Generates an array with coordinates of triangle midpoints (same dimension as connect)"""
    xyz = np.zeros_like(connect, dtype=float)
    xyz = (1./3.)*(geom[connect[:,0],:]+geom[connect[:,1],:]+geom[connect[:,2],:])   
    return xyz


def GetTimeIndices(xdmfFile):
    """Selects the correct time steps for output that contains two events"""    
    if args.events[0] == 0:
        return [0, xdmfFile.ReadNdt()]
    elif args.events[0] == 1:
        return [0, int(xdmfFile.ReadNdt()/2)]
    else:
        return [int(xdmfFile.ReadNdt()/2), xdmfFile.ReadNdt()]


def LoadSliprates(sliprateParameters):
    sliprates = []
    if "SRs" in sliprateParameters:
        sliprates.append(faultXdmf.ReadData("SRs")[timeIndicesFault[0]:timeIndicesFault[1]].T)
    if "SRd" in sliprateParameters:
        sliprates.append(faultXdmf.ReadData("SRd")[timeIndicesFault[0]:timeIndicesFault[1]].T)
    if "ASR" in sliprateParameters:
        if "SRs" in sliprateParameters and "SRd" in sliprateParameters:
            sliprates.append(np.sqrt(sliprates[0]**2 + sliprates[1]**2))
        else:
            sliprates.append(np.sqrt(faultXdmf.ReadData("SRs")[timeIndicesFault[0]:timeIndicesFault[1]].T**2 +
                                     faultXdmf.ReadData("SRd")[timeIndicesFault[0]:timeIndicesFault[1]].T**2))
    return np.array(sliprates)


def ClosestReceiverIndex(coords, receiver_array):
    dist = np.sum((receiver_array-coords)**2, axis=1)
    return np.argmin(dist)


def PTraveltime(receiverIndex, trigger=0.001):
    """Picks the P-wave arrival. This function is only suited for synthetic data without noise"""
    
    timeIndices = GetTimeIndices(surfaceXdmf)
    absWaveform = np.sqrt(surfaceXdmf.ReadDataChunk(surfaceVariables[0], receiverIndex, 1)**2 +
                      surfaceXdmf.ReadDataChunk(surfaceVariables[1], receiverIndex, 1)**2 +
                      surfaceXdmf.ReadDataChunk(surfaceVariables[2], receiverIndex, 1)**2)[timeIndices[0]:timeIndices[1]]
    
    threshold = np.amax(absWaveform)*trigger           # trigger dictates the portion of the maximum ground velocity,
    t = np.argmax(absWaveform>threshold)*dtSurface     # which needs to be reached to pick the p-phase
    return t                                           # p-phase arrival time (equal to traveltime if nucleation at t=0)


def ApproximateHypocenter(absoluteSliprates, faultCells, slipRateThreshold=args.sliprateThreshold[0]):
    
    slippingFault = np.where(absoluteSliprates > slipRateThreshold, 1, 0)
    slippingFault = np.where(slippingFault > 0, slippingFault - slippingFault[:,0].reshape(-1,1), 0) 
    # ignore elements that already slip at t=0 (important when the output contains more than one event)
    
    indHypocenter = np.argwhere(slippingFault[:,np.argmax(np.sum(slippingFault, axis=0) > 100)]==1)
    hypocenter = np.mean(faultCells[indHypocenter,:], axis=0)[0]
    return hypocenter


def GetTraveltimesToAllFaultCells(receiverCoords, faultCells, swaveVelocity):
    
    distances = np.sqrt(np.sum((receiverCoords - faultCells) ** 2, axis=1))
    return distances * (1 / swaveVelocity)


def MedianSmoothing(input_tuple):
    # input_tuple: (parameter, xyz, indices, radius)
    parameter = input_tuple[0]
    xyz = input_tuple[1]
    indices = input_tuple[2]
    radius = input_tuple[3]
    parameter_new = np.zeros(indices.shape[0])
    for i, j in enumerate(indices):
        parameter_new[i] = np.median(parameter[(xyz[:,0] < xyz[j,0]+radius) & (xyz[:,0] > xyz[j,0]-radius) & 
                            (xyz[:,1] < xyz[j,1]+radius) & (xyz[:,1] > xyz[j,1]-radius) &                    
                            (xyz[:,2] < xyz[j,2]+radius) & (xyz[:,2] > xyz[j,2]-radius)])
    return parameter_new


def MedianSmoothingParallel(parameter, xyz, radius=250., nprocs=4):
    
    chunks = [tuple((parameter, xyz, i, radius)) for i in np.array_split(np.arange(parameter.shape[0]), nprocs)]
        
    p = Pool(nprocs)
    result = p.map(MedianSmoothing, chunks)
    p.close()
    p.join()
    parameter_new = np.concatenate(result)
    return parameter_new


start = timeit.default_timer()

print("Loading data...")
surfaceXdmf = seissolxdmf.seissolxdmf(args.filename)
surfaceCells = ComputeTriangleMidpoints(surfaceXdmf.ReadGeometry(), surfaceXdmf.ReadConnect())
dtSurface = surfaceXdmf.ReadTimeStep()

if not args.faultXdmf:
    args.faultXdmf = args.filename[:-12]+"fault.xdmf"
    
faultXdmf = seissolxdmf.seissolxdmf(args.faultXdmf) 
faultGeom = faultXdmf.ReadGeometry()
faultConnect = faultXdmf.ReadConnect()
faultCells = ComputeTriangleMidpoints(faultGeom, faultConnect)
dtFault = faultXdmf.ReadTimeStep()
timeIndicesFault = GetTimeIndices(faultXdmf)
ruptureTimes = faultXdmf.ReadData("RT", timeIndicesFault[1]-1)

if args.events[0] == 2:
    ruptureTimes = ruptureTimes - faultXdmf.ReadData("RT", timeIndicesFault[0]-1)
    ruptureTimes = np.where(ruptureTimes > 0., ruptureTimes - dtFault * timeIndicesFault[0], 0.)

if args.sliprateParameters == "absolute":
    sliprateParameters = ["ASR"]
elif args.sliprateParameters == "components":
    sliprateParameters = ["SRs", "SRd"]
else:
    sliprateParameters = ["SRs", "SRd", "ASR"] 
    
sliprates = LoadSliprates(sliprateParameters)
sliprates = np.abs(sliprates)

print(f"surfaceCells.shape: {surfaceCells.shape}")
print(f"sliprates.shape: {sliprates.shape}")

stop1 = timeit.default_timer()
print(f"Time to load data: {np.round(stop1 - start, 2)}")

print("Computing travel times and peak slip rate times...")

receiverIndex = ClosestReceiverIndex(np.array([args.x[0], args.y[0], args.z[0]]), surfaceCells)
receiverCoords = surfaceCells[receiverIndex]

print(f"Closest receiver: {np.int_(receiverCoords)}")

surfaceVariables = ['u', 'v', 'w'] if 'u' in surfaceXdmf.ReadAvailableDataFields() else ['v1', 'v2', 'v3']

receiverPwaveTraveltime = PTraveltime(receiverIndex, trigger = 0.01 if args.events[0]==2 else 0.001)
print(f"P-wave traveltime to receiver: {np.round(receiverPwaveTraveltime, 2)}")

absoluteSliprates = sliprates[-1,:,:] if "ASR" in sliprateParameters else np.linalg.norm(sliprates, axis=0)
hypocenter = ApproximateHypocenter(absoluteSliprates, faultCells, slipRateThreshold=args.sliprateThreshold[0])
print(f"Hypocenter approximated at: {np.round(hypocenter)}")

avgSwaveVelocity = np.linalg.norm(receiverCoords - hypocenter) / receiverPwaveTraveltime / np.sqrt(3) # Assumes a Poisson solid
print(f"Average S-wave velocity between fault and receiver: {np.round(avgSwaveVelocity,2)}")

allFaultTraveltimes = GetTraveltimesToAllFaultCells(receiverCoords, faultCells, avgSwaveVelocity)
print(f"Mean traveltime between receiver and fault cells: {np.round(np.mean(allFaultTraveltimes), 2)}")

print(f"Sliprate threshold: {args.sliprateThreshold[0]}")
sliprates = np.where(sliprates >= args.sliprateThreshold[0], sliprates, 0.)  # Set slip rates below the threshold to zero
peakSliprateTimes = np.argmax(sliprates, axis=2) * dtFault

stop2 = timeit.default_timer()
print(f"Time to compute values: {np.round(stop2 - stop1, 2)}")

if args.medianFilterRadius[0] > 0.:
    
    print(f"Median filter radius: {args.medianFilterRadius[0]}")
    assert(args.threads[0] >= 1 and args.threads[0] <= cpu_count())
    print(f"Median smoothing using {args.threads[0]} threads...")
    
    for i in range(peakSliprateTimes.shape[0]):
        peakSliprateTimes[i] = MedianSmoothingParallel(peakSliprateTimes[i], faultCells, radius=args.medianFilterRadius[0], 
                                                       nprocs=args.threads[0])
    stop3 = timeit.default_timer()
    print(f"Time to apply median filter: {np.round(stop3 - stop2, 2)}")

isochroneTimes = np.zeros((peakSliprateTimes.shape[0]+1, peakSliprateTimes.shape[1]))
isochroneTimes[0] = np.where(ruptureTimes > 0., ruptureTimes+allFaultTraveltimes, 0.)
isochroneTimes[1:] = np.where(peakSliprateTimes > 0., peakSliprateTimes+allFaultTraveltimes, 0.)

stop3 = timeit.default_timer()
print("Saving output...")

prefix = f"isochrones_{np.int(receiverCoords[0])}_{np.int(receiverCoords[1])}_{np.int(receiverCoords[2])}"
if args.medianFilterRadius[0] > 0.:
    prefix += f"_r{np.int(args.medianFilterRadius[0])}"
sliprateParameters.insert(0, "RT")

if args.output == "xdmf" or args.output == "both":

    outDataNames = ["iso_time_"+ i for i in sliprateParameters]
    outData = [isochroneTimes[i] for i in range(len(outDataNames))]

    sxw.write_seissol_output(prefix, faultGeom, faultConnect, outDataNames, outData, dtFault, [0])

if args.output == "numpy" or args.output == "both":
    for i in sliprateParameters:
        prefix = prefix+"_"+i

    np.save(prefix+"_xyz", np.concatenate((isochroneTimes.T, faultCells), axis=1))

stop4 = timeit.default_timer()
print(f"Time to save output: {np.round(stop4 - stop3, 2)}")
print(f"Total time: {np.round(stop4 - start, 2)}")
print("Done. Goodbye.")