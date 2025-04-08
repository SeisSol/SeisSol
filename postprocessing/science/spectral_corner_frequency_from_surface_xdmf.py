# Script for computing spectral corner frequencies of surface waveforms from seissol surface.xdmf output
# Output of this script is either a numpy or xdmf file

import numpy as np
import seissolxdmf as sx
import seissolxdmfwriter as sxw
import argparse
import timeit
import os
import sys
import scipy as sp
from scipy import signal
from multiprocessing import Pool,cpu_count

parser = argparse.ArgumentParser(description='calculate spectral corner frequencies from surface xdmf file')
parser.add_argument('filename', help='path+filename-surface.xdmf')
parser.add_argument('--step', nargs=1, default=([1]), metavar=('step'), type=int,
                    help='choose step size to skip receivers for faster computation (only works for numpy output)')

parser.add_argument('--MP', nargs=1, metavar=('ncpu'), default=([1]),
                    help='use np.pool to speed-up calculations' ,type=int)
        
parser.add_argument('--maxFreq', nargs=1, default=([1.0]), metavar=('maxFreq'),
                    help='maximum frequency of the spectrum used for the corner frequency computation' ,type=float)
                    # maxFreq that should be sufficiently resolved by the simulation
    
parser.add_argument('--maxCornerFreq', nargs=1, metavar=('maxCornerFreq'),
                    help='maximum possible corner frequency (default = maxFreq)' ,type=float)
                    # Should be specified to reduce computational costs, if the maximum expected corner frequency is
                    # significantly less than maxFreq
        
parser.add_argument('--output', choices=["numpy","xdmf","both"], default="both", 
                    help='choose the output format')

parser.add_argument('--bodyWaveWindow',  dest='bodyWaveWindow', action='store_true' , default = False, 
                    help='cuts out a body wave window before performing the fft to mitigate the impact of surface waves')
                    # -fault.xdmf file needed (within the same directory as -surface.xdmf or arg --faultxdmf needed)
    
parser.add_argument('--slipRateThreshold', nargs=1, default=([0.1]), metavar=('slipRateThreshold'),
                    help='slip rate threshold that is used to approximate event duration and hypocenter location' ,type=float)
                    # only needed when --bodyWaveWindow is active
    
parser.add_argument('--events', nargs=1, choices=[0,1,2], default=([0]), 
                    help='0: xdmf file contains just one event; 1: process the first of two events; 2: process the second of two events' ,type=int) 
                    # Needed for Ridgecrest output, where the xdmf file contains two events
    
parser.add_argument('--rotate',  dest='rotate', action='store_true' , default = False, 
                    help='horizontal components are rotated w.r.t. the slip centroid') 
                    # -fault.xdmf file needed
    
parser.add_argument('--outputprefix', default="data", type=str, metavar=('outputprefix'), 
                    help='provide an outputprefix')

parser.add_argument('--faultXdmf', type=str, metavar=('-fault.xdmf'), 
                    help='provide path+filename-fault.xdmf; only needed when the path differs from the -surface.xdmf file')

parser.add_argument('--parallelLoading',  dest='parallelLoading', action='store_true' , default = False, 
                    help='load components in parallel (weak scaling)')
                    # This can lead to an enconding error of multiprocessing, when the input arrays are too large

parser.add_argument('--avgSWaveVelocity', nargs=1, metavar=('avgSWaveVelocity'),
                    help='provide an average S-wave velocity near the surface (otherwise it is approximated from P-arrivals)',
                    type=float)
args = parser.parse_args()

def LoadSingleComponentParallel(component, surfaceData=True):
    chunks = np.int_(np.linspace(0, nElements, nprocs+1)) if surfaceData else np.int_(np.linspace(0, nElementsFault, nprocs+1))
    chunks = np.array([chunks[:-1], chunks[1:]-chunks[:-1], np.full(nprocs, component)]).T
    chunks = [tuple(i) for i in chunks]
    
    if surfaceData:
        p = Pool(nprocs)
        result = p.map(ReadSurfaceWithPool, chunks)
        p.close()
        p.join()
        return np.concatenate(result, axis=1).T[::stepsize]
    else:
        p = Pool(nprocs)
        result = p.map(ReadFaultWithPool, chunks)
        p.close()
        p.join()
        return np.concatenate(result, axis=1).T[::stepsize]
    
def ReadSurfaceWithPool(chunk):
    return surfacexdmf.ReadDataChunk(chunk[2], int(chunk[0]), int(chunk[1]))[timeIndices[0]:timeIndices[1]]

def ReadFaultWithPool(chunk):
    return faultxdmf.ReadDataChunk(chunk[2], int(chunk[0]), int(chunk[1]))[timeIndicesFault[0]:timeIndicesFault[1]]

def ComputeTriangleMidpoints(geom, connect):
    """Generates an array with coordinates of triangle midpoints (same dimension as connect)"""
    xyz = np.zeros_like(connect, dtype=float)
    xyz = (1./3.)*(geom[connect[:,0],:]+geom[connect[:,1],:]+geom[connect[:,2],:])   
    return xyz

def CalculateSlipCentroid(faultxdmf, faultxyz):
    ASl = faultxdmf.ReadData("ASl", idt=timeIndicesFault[1]-1).T
    if timeIndicesFault[0]!=0:
        ASl -= faultxdmf.ReadData("ASl", idt=timeIndicesFault[0]).T   
    return np.average(faultxyz, axis=0, weights=ASl)

def ComputeBackazimuth(xyz, centroid):
    ba = np.arctan2(xyz[:,0]-centroid[0], xyz[:,1]-centroid[1]) + np.pi
    return ba

def RotateToRadialTransversal(u, v, xyz, centroid):
    ba = ComputeBackazimuth(xyz, centroid)
    ba = ba.reshape(-1,1)
    r = -u * np.sin(ba) - v * np.cos(ba)
    t = -u * np.cos(ba) + v * np.sin(ba)
    return r, t

def ComputeCornerFrequencyParallel(component, pphases, windowduration):
    chunk_inds = np.linspace(0, component.shape[0], nprocs+1, dtype=int)
    chunks = []
    for i in range(nprocs):
        chunks.append((component[chunk_inds[i]:chunk_inds[i+1]], pphases[chunk_inds[i]:chunk_inds[i+1]], 
                     windowduration[chunk_inds[i]:chunk_inds[i+1]]))
       
    p = Pool(nprocs)
    result = p.map(ComputeCornerFrequency, chunks)
    p.close()
    p.join()
    return np.concatenate(result)        


def ComputeCornerFrequency(tuples):
    """Preproccesing of a single waveform before calling fc_Qzero_gridsearch to determine the corner frequency"""
    
    fc_arr = np.zeros(tuples[0].shape[0])
    endfreq = args.maxFreq[0]
    startfreq = np.max([minFreq, 0.01])
    
    for i in range(0, fc_arr.size):
        
        velocity = tuples[0][i]
        windowduration = tuples[2][i] * 1.05     # add 5% to not lose information due to the taper
        pphase = tuples[1][i]
    
        # preprocessing (taper and padding)
        if windowduration != 0:
            data = cut_window_and_pad(velocity, windowduration, dt, pphase=pphase, alpha=0.05)
        else:
            taper = sp.signal.tukey(velocity.shape[0], alpha=0.05, sym=True)
            data = velocity * taper

        # fourier transform and down trimming to maxFreq
        data, freqs = complete_fft(data, dt)
        maxFreqIndices = int(np.argwhere(freqs > endfreq)[1])
        data = data[1:maxFreqIndices]
        freqs = freqs[1:maxFreqIndices]

        # integrate within the frequency domain and convert to amplitude spectrum
        data = np.abs(data*(1/(1j*2*np.pi*freqs)))

        newfreqs = np.arange(startfreq, endfreq+0.005, 0.005)
        data = np.interp(newfreqs, freqs, data)

        fc_arr[i] = fc_Qzero_gridsearch(data, newfreqs, misfit="SSM", n=2, 
                                 fc_end = 0 if not args.maxCornerFreq else args.maxCornerFreq[0])[0]
    return fc_arr

def cut_window_and_pad(data, window_duration, dt, pphase=0, alpha=0.05):
    """Cuts a tapered body wave window and pads remaining timesteps with zeros"""
    
    window_len = int(window_duration / dt)
    ind1 = np.max([int(pphase / dt - 0.025 * window_len), 0])
  
    if window_len + ind1 >= data.shape[0]:
        window_len = data.shape[0] - ind1
    
    taper = sp.signal.tukey(window_len, alpha=alpha, sym=True)
    data2 = data[ind1:ind1+window_len] * taper
    
    result = np.zeros_like(data)
    result[ind1:ind1+window_len] = data2
    
    return result

def fc_Qzero_gridsearch(amplitude_tr, freq_tr, misfit="SSM", n=2, fc_end=0., q_dim=0):
    """Determines the corner frequency of a displacement spectrum by performing a gridsearch 
    to find the best fitting Brune-type spectrum"""
    
    misfit_types = ["RMS", "SSM"]
    if misfit not in misfit_types:
        raise ValueError("Invalid misfit type. Expected one of: %s" % misfit_types)
    
    # Prepare tensors for broadcasting
    fc_ind = -2 if fc_end == 0. else int(np.argwhere(freq_tr>fc_end)[0])
         
    fc_tensor = np.empty((1,freq_tr[:fc_ind].shape[0], 1))
    fc_tensor[0,:,0] = freq_tr[:fc_ind]
    
    freq_tensor = np.empty((freq_tr.shape[0],1,1))
    freq_tensor[:,0,0] = freq_tr
    
    amplitude_tensor = np.empty((freq_tr.shape[0],1,1))
    amplitude_tensor[:,0,0] = amplitude_tr
    
    # Initialise tensor of possible Q_zero values depending on amplitude_tr[0] or take mean amplitude below corner frequency
    if q_dim != 0:
        q_zero = np.array([[np.geomspace(0.1*amplitude_tr[0], 10*amplitude_tr[0], num=q_dim)]])
    else:
        q_zero = np.empty((1,freq_tr[:fc_ind].shape[0], 1))
        for i in range(q_zero.shape[1]):
            q_zero[0,i,0] = np.mean(amplitude_tr[:i+1])
    
    # Calculate Brune spectra of every (Q, fc) combination by broadcasting
    brune_function = q_zero*(1/(1+(freq_tensor/fc_tensor)**n))
    
    # Calculate the misfit of every (Q, fc) combination
    if misfit == "RMS":
        misfit_tensor = (brune_function - amplitude_tensor)**2
        misfit_tensor = np.sqrt(np.mean(misfit_tensor, axis=0))
    if misfit == "SSM":
        misfit_tensor = np.abs(np.log10(amplitude_tensor/brune_function))
        misfit_tensor = np.mean(misfit_tensor, axis=0)
    
    # Get indices of minimum misfit (Q, fc) combination
    index = np.unravel_index(np.argmin(misfit_tensor), misfit_tensor.shape)
    
    fc = fc_tensor[0,index[0],0]
    if q_dim != 0:
        q = q_zero[0,0,index[1]]
    else:
        q = q_zero[0,index[0],0]
        
    return fc, q

def complete_fft(tr, dt):
    tr_fft=np.fft.rfft(tr.data)
    times_comp = np.fft.rfftfreq(tr.shape[0],dt)
    return tr_fft, times_comp

def ApproximateEventDurationAndHypocenter(faultxdmf, faultxyz):
    """Calculates the derivative of ASl (probably faster than loading SR1 and SR2)
     and approximates event duration by multiplying dt with the number of timesteps,
     where maximum on-fault slip rate is above a threshold (slipRateThreshold)"""
    
    slipRateThreshold=args.slipRateThreshold[0]
    if nprocs == 1 or not args.parallelLoading:
        ASl = faultxdmf.ReadData("ASl").T[::stepsize,timeIndicesFault[0]:timeIndicesFault[1]]
    else:
        ASl = LoadSingleComponentParallel("ASl", surfaceData=False)
    ASR = np.abs(np.gradient(ASl, dtFault, axis=1))
    slippingElement = np.where(np.amax(ASR, axis=0) > slipRateThreshold, 1, 0)
    if np.all(slippingElement == 1):
        print("Warning: at least one element is always slipping")
    elif np.argmax(np.flip(slippingElement)) == 0:
        print("Warning: at least one element is still slipping at the last time step")
    duration = np.trim_zeros(slippingElement).size * dtFault
    slippingFault = np.where(ASR > slipRateThreshold, 1, 0)
    indHypocenter = np.argmax(slippingFault[:,np.argmax(np.amax(slippingFault, axis=0))])
    hypocenter = faultxyz[indHypocenter,:]
    return duration, hypocenter

def PickPPhases(trigger=0.001):
    """Picks the P-wave arrival time for each receiver. This function is only suited for synthetic data without noise"""
    
    data = np.sqrt(u**2 + v**2 + w**2)
    thresholds = np.amax(data, axis=1)*trigger  # trigger dictates the portion of the maximum ground velocity,
    thresholds = thresholds.reshape(-1,1)       # which needs to be reached to pick the p-phase
    ind = np.argmax(data>thresholds, axis=1)
    t = ind*dt
    return np.transpose(np.vstack((ind,t))) # p-phase output: first column index and second column simulation time

def compute_P2S_window_duration(xyz, pphase, hypocenter, centroid):
    """Approximates the difference between P- and S-wave arrivals for each receiver"""
    
    if not args.avgSWaveVelocity:
        distance1 = np.linalg.norm((xyz - hypocenter), axis=1)
        swaveVelocity = np.mean(distance1 / pphase) / np.sqrt(3)
    else:
        swaveVelocity = args.avgSWaveVelocity[0]
    print("Average S-wave velocity near the surface: "+str(swaveVelocity))
    distance2 = np.linalg.norm((xyz - centroid), axis=1)
    swaveDelay = distance2 / swaveVelocity * (1-1/np.sqrt(3))
    return swaveDelay

def GetTimeIndices(xdmfFile):
    if args.events[0] == 0:
        return [0, xdmfFile.ReadNdt()]
    elif args.events[0] == 1:
        return [0, int(xdmfFile.ReadNdt()/2)]
    else:
        return [int(xdmfFile.ReadNdt()/2), xdmfFile.ReadNdt()]

start = timeit.default_timer()

print("Loading data...")
surfacexdmf = sx.seissolxdmf(args.filename)
dt = surfacexdmf.ReadTimeStep()
nElements = surfacexdmf.ReadNElements()
stepsize = args.step[0]
surfacexyz = ComputeTriangleMidpoints(surfacexdmf.ReadGeometry(), surfacexdmf.ReadConnect())[::stepsize]

if stepsize != 1 and args.output != "numpy":
    sys.exit("Xdmf output is not compatible with stepsize != 1")

timeIndices = GetTimeIndices(surfacexdmf)
minFreq = 1. / ((timeIndices[1]-timeIndices[0]-1.) * dt)
print(f"Minimum frequency: {minFreq}")

nprocs  = args.MP[0]
assert(nprocs<=cpu_count())
print(f"Using {nprocs} ranks")

if not args.parallelLoading:
    print("Serial loading of components...")
    
surfaceVariables = ['u', 'v', 'w'] if 'u' in surfacexdmf.ReadAvailableDataFields() else ['v1', 'v2', 'v3']
    
if nprocs == 1 or not args.parallelLoading:
    u = surfacexdmf.ReadData(surfaceVariables[0]).T[::stepsize,timeIndices[0]:timeIndices[1]]
    v = surfacexdmf.ReadData(surfaceVariables[1]).T[::stepsize,timeIndices[0]:timeIndices[1]]
    w = surfacexdmf.ReadData(surfaceVariables[2]).T[::stepsize,timeIndices[0]:timeIndices[1]]
else:
    u = LoadSingleComponentParallel(surfaceVariables[0])
    v = LoadSingleComponentParallel(surfaceVariables[1])
    w = LoadSingleComponentParallel(surfaceVariables[2])
print("Shape of components: "+str(u.shape))    

stop1 = timeit.default_timer()
print('Time to load data: ', stop1 - start)

print("Preparing data...")

if args.rotate or args.bodyWaveWindow:
    if not args.faultXdmf:
        args.faultXdmf = args.filename[:-12]+"fault.xdmf"
    faultxdmf = sx.seissolxdmf(args.faultXdmf) 
    faultxyz = ComputeTriangleMidpoints(faultxdmf.ReadGeometry(), faultxdmf.ReadConnect())
    timeIndicesFault = GetTimeIndices(faultxdmf)
    slipCentroid = CalculateSlipCentroid(faultxdmf, faultxyz)
    print("Calculated centroid: "+str(slipCentroid))   

if args.rotate:
    u, v = RotateToRadialTransversal(u, v, surfacexyz, slipCentroid)
    
if args.bodyWaveWindow:
    pPhases = PickPPhases(trigger=0.001)[:,1]
    dtFault = faultxdmf.ReadTimeStep()
    nElementsFault = faultxdmf.ReadNElements()
    eventDuration, hypocenter = ApproximateEventDurationAndHypocenter(faultxdmf, faultxyz)
    print("Slip rate threshold: "+ str(args.slipRateThreshold[0]))
    print("Approximated event duration: " + str(eventDuration))
    print("Approximated hypocenter: " + str(hypocenter))
    windowDuration = compute_P2S_window_duration(surfacexyz, pPhases, hypocenter, slipCentroid) + eventDuration
    print("Mean body wave window duration: "+str(np.mean(windowDuration)))
else:
    pPhases = np.zeros(nElements)[::stepsize]
    windowDuration = np.zeros(nElements)[::stepsize]
    
    
stop2 = timeit.default_timer()
print('Time to prepare data: ', stop2 - stop1)


print("Computing corner frequencies...")
fc_1 = ComputeCornerFrequencyParallel(u, pPhases, windowDuration)
fc_2 = ComputeCornerFrequencyParallel(v, pPhases, windowDuration)
fc_3 = ComputeCornerFrequencyParallel(w, pPhases, windowDuration)
print("Shape of corner frequency arrays: "+str(fc_3.shape)) 

stop3 = timeit.default_timer()
print('Time to compute corner frequencies: ', stop3 - stop2)


print("Saving output...")
output_strings = ["fc_radial", "fc_transversal", "fc_vertical"] if args.rotate else ["fc_v1", "fc_v2", "fc_v3"]
prefix = args.outputprefix
    
if args.output == "numpy" or args.output == "both":
    np.save(prefix+'_'+output_strings[0], fc_1)
    np.save(prefix+'_'+output_strings[1], fc_2)
    np.save(prefix+'_'+output_strings[2], fc_3)
    np.save(prefix+"_fc_xyz", surfacexyz)
    
if args.output == "xdmf" or args.output == "both":
    connect = surfacexdmf.ReadConnect()[::stepsize]
    geom = surfacexdmf.ReadGeometry()
    aDataName = output_strings
    fn = prefix+'_fc'
    fc_list = [fc_1, fc_2, fc_3]
    sxw.write_seissol_output(fn, geom, connect, aDataName, fc_list, dt, [0])
    
stop4 = timeit.default_timer()
print('Time to save output: ', stop4 - stop3)
print('Total time: ', stop4 - start)
print("Done.")
