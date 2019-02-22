import argparse
from submodules.pythonXdmfReader.pythonXdmfReader import *
import numpy as np
parser = argparse.ArgumentParser(description='estimate minimum cohesive length from fault output')
parser.add_argument('filename', help='fault output filename (xdmf)')
args = parser.parse_args()


xdmfFilename = args.filename
ndt = ReadNdt(xdmfFilename)
nElements = ReadNElements(xdmfFilename)
DS, data_prec = LoadData(xdmfFilename, 'DS', nElements, ndt-1, oneDtMem=True)
Vr, data_prec = LoadData(xdmfFilename, 'Vr', nElements, ndt-1, oneDtMem=True)
RT, data_prec = LoadData(xdmfFilename, 'RT', nElements, ndt-1, oneDtMem=True)
ASl, data_prec = LoadData(xdmfFilename, 'RT', nElements, ndt-1, oneDtMem=True)
cohesiveZoneWidth = (DS-RT)*Vr
cohesiveZoneWidth = cohesiveZoneWidth[DS>0]

ndsp = (len(np.where(DS>0)[0]))
naslp = (len(np.where(ASl>0.1)[0]))

print('DS defined over %d pc of elements that sliped more that 10cm' %(ndsp/naslp*100))

for i in range(5,100,5):
   print('%dth percentile: %.2f' %(i, np.percentile(cohesiveZoneWidth, i)))
