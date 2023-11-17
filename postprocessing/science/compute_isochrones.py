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
                    help='output format')

parser.add_argument('--slipRateThreshold', nargs=1, default=([0.05]), metavar=('peak slip rate threshold in m/s'),
                    help='fault cells below the threshold are assumed to not radiate waves' , type=float)

parser.add_argument('--slipRateParameters', choices=["absolute","components","both"], default="both", 
                    help="absolute: compute isochrones from peak times of absolute slip rates;"+ 
                    " components: peak slip rate times are computed for SRs and SRd separately")

parser.add_argument('--medianFilterRadius', nargs=1, required=True, metavar=('radius of the median filter in m'),
                    help='peak slip rate time of a fault cell is set to the median value of all surrounding fault cells within \
                    a sphere with this radius; 0 turns the filter off' , type=float)
# This prevents strong oscillations of peak slip rate times where different rupture fronts have similar peak slip rates 
# and picks the time of the dominant rupture front.
# Radii between 150 m and 300 m are advised (also depending on the resolution of the fault output)

parser.add_argument('--threads', nargs=1, metavar=('threads'), default=([4]),
                    help='use multiprocessing to speed up the median filtering' ,type=int)

parser.add_argument('--timeStepRange', nargs=2, default=([0,-1]), metavar=("minIndex", "maxIndex"),
                    help='parse the fault output time step index range that includes the desired part of the rupture; e.g., \
                    if the output contains more than one event; first time step has index 0', type=int) 

args = parser.parse_args()


def ComputeTriangleMidpoints(geom, connect):
    """Generates an array with coordinates of triangle midpoints (same dimension as connect)"""
    xyz = (1./3.)*(geom[connect[:,0],:]+geom[connect[:,1],:]+geom[connect[:,2],:])   
    return xyz


def LoadSliprates(sliprateParameters):
    sliprates = []
    if "SRs" in sliprateParameters:
        sliprates.append(faultXdmf.ReadData("SRs")[indexBoundsFault[0]:indexBoundsFault[1]].T)
    if "SRd" in sliprateParameters:
        sliprates.append(faultXdmf.ReadData("SRd")[indexBoundsFault[0]:indexBoundsFault[1]].T)
    if "ASR" in sliprateParameters:
        if "SRs" in sliprateParameters and "SRd" in sliprateParameters:
            sliprates.append(np.sqrt(sliprates[0]**2 + sliprates[1]**2))
        else:
            sliprates.append(np.sqrt(faultXdmf.ReadData("SRs")[indexBoundsFault[0]:indexBoundsFault[1]].T**2 +
                                     faultXdmf.ReadData("SRd")[indexBoundsFault[0]:indexBoundsFault[1]].T**2))
    return np.array(sliprates)


def ClosestReceiverIndex(coords, receiver_array):
    dist = np.sum((receiver_array-coords)**2, axis=1)
    return np.argmin(dist)


def PTraveltime(receiverIndex, trigger=0.01):
    """Picks the P-wave arrival. This function is only suited for synthetic data with negligible noise"""
    
    indexBoundsSurface = np.int_(indexBoundsFault * (dtFault / dtSurface))
    absWaveform = np.sqrt(surfaceXdmf.ReadDataChunk(surfaceVariables[0], receiverIndex, 1)**2 +
                      surfaceXdmf.ReadDataChunk(surfaceVariables[1], receiverIndex, 1)**2 +
                      surfaceXdmf.ReadDataChunk(surfaceVariables[2],
                                                receiverIndex, 1)**2)[indexBoundsSurface [0]:indexBoundsSurface [1]]
    
    threshold = np.amax(absWaveform)*trigger           # trigger dictates the portion of the maximum ground velocity,
                                                       # which needs to be reached to pick the p-phase
    return np.argmax(absWaveform>threshold)*dtSurface  # p-phase arrival time (equal to traveltime if nucleation at t=0)


def ApproximateHypocenter(ruptureTimes, faultMidPoints):
    
    nucleationTime = dtFault
    indHypocenter = (ruptureTimes > 0) & (ruptureTimes <= nucleationTime)
    
    # Take the mean of the first 100+ slipping elements to reduce impact from degenerated elements
    while indHypocenter.sum() < 100:
        nucleationTime += dtFault
        indHypocenter = (ruptureTimes > 0) & (ruptureTimes <= nucleationTime)
    
    hypocenter = np.mean(faultMidPoints[indHypocenter,:], axis=0)
    return hypocenter


def GetTraveltimesToAllFaultCells(receiverCoords, faultMidPoints, swaveVelocity):
    
    distances = np.sqrt(np.sum((receiverCoords - faultMidPoints) ** 2, axis=1))
    return distances * (1 / swaveVelocity)


def MedianSmoothing(input_tuple):
    # input_tuple: (peakSliprateTimesLocal, xyz, indices, radius)
    peakSliprateTimesLocal, xyz, indices, radius = input_tuple
    peakSliprateTimesSmooth = np.zeros(indices.shape[0])
    for i, j in enumerate(indices):
        xyz0=  xyz[j,:]
        distances = np.linalg.norm(xyz - xyz0, axis=1)
        indices_within_sphere = np.argwhere(distances <= radius)
        peakSliprateTimesSmooth[i] =np.median(peakSliprateTimesLocal[indices_within_sphere])
    return peakSliprateTimesSmooth


def MedianSmoothingParallel(peakSliprateTimesLocal, xyz, radius=250., nprocs=4):
    
    chunks = [(peakSliprateTimesLocal, xyz, i, radius) for i in np.array_split(np.arange(peakSliprateTimesLocal.shape[0]),
                                                                               nprocs)]
        
    p = Pool(nprocs)
    result = p.map(MedianSmoothing, chunks)
    p.close()
    p.join()
    return np.concatenate(result)


start = timeit.default_timer()

print("Loading data...")
surfaceXdmf = seissolxdmf.seissolxdmf(args.filename)
surfaceCellBarycenters = ComputeTriangleMidpoints(surfaceXdmf.ReadGeometry(), surfaceXdmf.ReadConnect())
dtSurface = surfaceXdmf.ReadTimeStep()

if not args.faultXdmf:
    args.faultXdmf = args.filename[:-12]+"fault.xdmf"
    
faultXdmf = seissolxdmf.seissolxdmf(args.faultXdmf) 
faultGeom = faultXdmf.ReadGeometry()
faultConnect = faultXdmf.ReadConnect()
faultMidPoints = ComputeTriangleMidpoints(faultGeom, faultConnect)
dtFault = faultXdmf.ReadTimeStep()
indexBoundsFault = np.array(args.timeStepRange)
indexBoundsFault[1] = faultXdmf.ReadNdt()-1 if indexBoundsFault[1] == -1 else indexBoundsFault[1]

# set rupture times before indexBoundsFault[0] to 0:
ruptureTimes = faultXdmf.ReadData("RT", indexBoundsFault[1]) - faultXdmf.ReadData("RT", indexBoundsFault[0])
# shift rupture times so that indexBoundsFault[0] corresponds to 0:
ruptureTimes = np.where(ruptureTimes > 0., ruptureTimes - dtFault * indexBoundsFault[0], 0.)


if args.slipRateParameters == "absolute":
    sliprateParameters = ["ASR"]
elif args.slipRateParameters == "components":
    sliprateParameters = ["SRs", "SRd"]
else:
    sliprateParameters = ["SRs", "SRd", "ASR"] 
    
sliprates = LoadSliprates(sliprateParameters)
sliprates = np.abs(sliprates)

print(f"surfaceCellBarycenters.shape: {surfaceCellBarycenters.shape}")
print(f"sliprates.shape: {sliprates.shape}")

stop1 = timeit.default_timer()
print(f"Time to load data: {stop1 - start:.2f}")

print("Computing travel times and peak slip rate times...")

receiverIndex = ClosestReceiverIndex(np.array([args.x[0], args.y[0], args.z[0]]), surfaceCellBarycenters)
receiverCoords = surfaceCellBarycenters[receiverIndex]

print(f"Closest receiver: {np.int_(receiverCoords)}")

surfaceVariables = ['u', 'v', 'w'] if 'u' in surfaceXdmf.ReadAvailableDataFields() else ['v1', 'v2', 'v3']

receiverPwaveTraveltime = PTraveltime(receiverIndex, trigger = 0.01)
print(f"P-wave traveltime to receiver: {receiverPwaveTraveltime:.2f}")

absoluteSliprates = sliprates[-1,:,:] if "ASR" in sliprateParameters else np.linalg.norm(sliprates, axis=0)
hypocenter = ApproximateHypocenter(ruptureTimes, faultMidPoints)
print(f"Hypocenter approximated at: {np.round(hypocenter)}")

avgSwaveVelocity = np.linalg.norm(receiverCoords - hypocenter) / receiverPwaveTraveltime / np.sqrt(3) # Assumes a Poisson solid
print(f"Average S-wave velocity between hypocenter and receiver: {avgSwaveVelocity:.2f}")

allFaultTraveltimes = GetTraveltimesToAllFaultCells(receiverCoords, faultMidPoints, avgSwaveVelocity)
print(f"Mean traveltime between receiver and fault cells: {np.mean(allFaultTraveltimes):.2f}")

print(f"Sliprate threshold: {args.slipRateThreshold[0]}")
sliprates = np.where(sliprates >= args.slipRateThreshold[0], sliprates, 0.)  # Set slip rates below the threshold to zero
peakSliprateTimes = np.argmax(sliprates, axis=2) * dtFault

stop2 = timeit.default_timer()
print(f"Time to compute values: {stop2 - stop1:.2f}")

if args.medianFilterRadius[0]:
    
    print(f"Median filter radius: {args.medianFilterRadius[0]}")
    assert(args.threads[0] >= 1 and args.threads[0] <= cpu_count())
    print(f"Median smoothing using {args.threads[0]} threads...")
    
    for i in range(peakSliprateTimes.shape[0]):
        peakSliprateTimes[i] = MedianSmoothingParallel(peakSliprateTimes[i], faultMidPoints, radius=args.medianFilterRadius[0], 
                                                       nprocs=args.threads[0])
    stop3 = timeit.default_timer()
    print(f"Time to apply median filter: {stop3 - stop2:.2f}")

isochroneTimes = np.zeros((peakSliprateTimes.shape[0]+1, peakSliprateTimes.shape[1]))
isochroneTimes[0] = np.where(ruptureTimes > 0., ruptureTimes+allFaultTraveltimes, 0.)
isochroneTimes[1:] = np.where(peakSliprateTimes > 0., peakSliprateTimes+allFaultTraveltimes, 0.)

stop3 = timeit.default_timer()
print("Saving output...")

prefix = f"isochrones_{receiverCoords[0]:.0f}_{receiverCoords[1]:.0f}_{receiverCoords[2]:.0f}"
if args.medianFilterRadius[0] > 0.:
    prefix += f"_r{int(args.medianFilterRadius[0])}"
sliprateParameters.insert(0, "RT")

if args.output == "xdmf" or args.output == "both":

    outDataNames = ["iso_time_"+ i for i in sliprateParameters]
    outData = [isochroneTimes[i] for i in range(len(outDataNames))]

    sxw.write_seissol_output(prefix, faultGeom, faultConnect, outDataNames, outData, dtFault, [0])

if args.output == "numpy" or args.output == "both":
    for i in sliprateParameters:
        prefix = prefix+"_"+i

    np.save(prefix+"_xyz", np.concatenate((isochroneTimes.T, faultMidPoints), axis=1))

stop4 = timeit.default_timer()
print(f"Time to save output: {stop4 - stop3:.2f}")
print(f"Total time: {stop4 - start:.2f}")
print("Done. Goodbye.")
