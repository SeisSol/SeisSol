# Script to determine the rise time (time the absolute slip rate is above a certain threshold) from seissol fault.xdmf output
# Output of this script is a binary file. If events = 1, the binary file can be exchanged with the rupture time fault output
# and opened in paraview without creating a new xdmf file.

import numpy as np
import seissolxdmf
import argparse
import sys

parser = argparse.ArgumentParser(description='calculate rise times from fault xdmf file')
parser.add_argument('filename', help='path+filename.xdmf')
parser.add_argument('--threshold', nargs=1, metavar=('threshold'), default=([0.01]), help='specify slip rate threshold' ,type=float)
parser.add_argument('--dt', nargs=1, metavar=('dt'), default=([0.1]), help='insert time step of xdmf file' ,type=float)
parser.add_argument('--events', nargs=1, metavar=('events'), default=([1]), help='number of events [1,2] in xdmf file' ,type=int)
args = parser.parse_args()

if args.events[0] != 1 and args.events[0] != 2:
    sys.exit("Script only compatible with events = [1,2]")


def get_rise_time(ASR_array, threshold=0.01, dt=0.1, events=1):
    
    RT = np.zeros_like(ASR_array)
    
    if events == 1:
        for i in range(ASR_array.shape[1]):
            count=0
            for j in range(ASR_array.shape[0]):
                if ASR_array[j,i] < threshold and count > 0: break
                if ASR_array[j,i] <  threshold and count == 0: continue
                if ASR_array[j,i] >= threshold: count += 1        
            RT[:,i] = (count * dt)
            
    if events == 2:
        mid = int(ASR_array.shape[0]/2)
        for i in range(ASR_array.shape[1]):
            count = 0
            for j in range(0, mid):
                if ASR_array[j,i] < threshold and count > 0: break
                if ASR_array[j,i] < threshold and count == 0: continue
                if ASR_array[j,i] >= threshold: count += 1        
            RT[:mid,i] = (count * dt)            
            
            count = 0
            for j in range(mid, ASR_array.shape[0]):
                if ASR_array[j,i] < threshold and count > 0: break
                if ASR_array[j,i] < threshold and count == 0: continue
                if ASR_array[j,i] >= threshold: count += 1       
            RT[mid:,i] = (count * dt)
            
    return RT

print("Loading file...")          
sx = seissolxdmf.seissolxdmf(args.filename)

print("Calculating absolute slip rate...")          
SRs = sx.ReadData('SRs')
SRd = sx.ReadData('SRd')
ASR = np.sqrt(SRs**2 + SRd**2)
          
print("Determining rise time...")
RT = get_rise_time(ASR, threshold=args.threshold[0], dt=args.dt[0], events=args.events[0])

print("Writing output file...")
fname = "RT.bin"
output_file = open(fname, "wb")
RT.tofile(output_file)
output_file.close()
print("Done writing RT.bin")
