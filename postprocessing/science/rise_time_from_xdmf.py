# Script to determine the rise time (time the absolute slip rate is above a certain threshold) from seissol fault.xdmf output
# Output of this script is a new xdmf file

import numpy as np
import seissolxdmf as sx
import seissolxdmfwriter as sxw
import argparse
import os

parser = argparse.ArgumentParser(description='calculate rise times from fault xdmf file')
parser.add_argument('filename', help='path+filename-fault.xdmf')
parser.add_argument('--threshold', nargs=1, metavar=('threshold'), default=([0.01]), help='specify slip rate threshold' ,type=float)
parser.add_argument('--events', nargs=1, metavar=('events'), choices=[1,2], default=([1]), help='number of events [1,2] in xdmf file' ,type=int)
args = parser.parse_args()

def counting_loop(row, array, threshold=0.01):
    """counting the time steps for which a certain point on the fault has an ASR above the threshold"""
    count = 0
    for j in range(array.shape[1]):
        if array[row,j] < threshold and count > 0: break
        if array[row,j] <  threshold and count == 0: continue
        if array[row,j] >= threshold: count += 1
    return count

def get_rise_time(ASR_array, threshold=0.01, dt=0.1, events=1):
    if events == 1:
        RT = np.zeros(ASR_array.shape[1])
        ASR_array=np.transpose(ASR_array) # to speed up
        for i in range(ASR_array.shape[0]):
            count = counting_loop(i, ASR_array, threshold=threshold) 
            RT[i] = (count * dt)
            
    if events == 2:
        RT = np.zeros_like(ASR_array[0:2,:])
        mid = int(ASR_array.shape[0]/2)
        ASR_array=np.transpose(ASR_array)
        for i in range(ASR_array.shape[0]):
            count = counting_loop(i, ASR_array[:,0:mid], threshold=threshold)      
            RT[0,i] = (count * dt) 
            count = counting_loop(i, ASR_array[:,mid:], threshold=threshold)    
            RT[1,i] = (count * dt)          
    return RT

def get_absolute_slip_rate(sx_object):
    SRs = sx_object.ReadData('SRs')
    SRd = sx_object.ReadData('SRd')
    ASR = np.sqrt(SRs**2 + SRd**2)
    return ASR

print("Loading file and calculating absolute slip rates...") 
sx = sx.seissolxdmf(args.filename)
dt = sx.ReadTimeStep()
ASR = get_absolute_slip_rate(sx)
          
print("Determining rise time...")
RT = get_rise_time(ASR, threshold=args.threshold[0], dt=dt, events=args.events[0])

print("Writing output file...")
connect = sx.ReadConnect()
xyz = sx.ReadGeometry()
aDataName = ['RT'] if args.events[0] == 1 else ['RT1', 'RT2']
prefix = os.path.splitext(args.filename)[0]
fn = prefix+'_RT'
lRT = [RT] if args.events[0] == 1 else [RT[0,:], RT[1,:]]
sxw.write_seissol_output(fn, xyz, connect, aDataName, lRT, dt, [0])
