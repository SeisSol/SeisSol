import argparse
import seissolxdmf as sx
import numpy as np
parser = argparse.ArgumentParser(description='estimate minimum cohesive length from fault output')
parser.add_argument('filename', help='fault output filename (xdmf)')
args = parser.parse_args()


xdmfFilename = args.filename
ndt = sx.ReadNdt(xdmfFilename)
nElements = sx.ReadNElements(xdmfFilename)
DS = sx.ReadData(xdmfFilename, 'DS', ndt-1)
Vr = sx.ReadData(xdmfFilename, 'Vr', ndt-1)
RT = sx.ReadData(xdmfFilename, 'RT', ndt-1)
ASl = sx.ReadData(xdmfFilename, 'ASl', ndt-1)
cohesiveZoneWidth = (DS-RT)*Vr
cohesiveZoneWidth = cohesiveZoneWidth[DS>0]

ndsp = (len(np.where(DS>0)[0]))
naslp = (len(np.where(ASl>0.1)[0]))

print('DS defined over %d pc of elements that sliped more that 10cm' %(ndsp/naslp*100))

for i in range(5,100,5):
   print('%dth percentile: %.2f' %(i, np.percentile(cohesiveZoneWidth, i)))
