#Author Thomas Ulrich, LMU
import numpy as np

# parsing python arguments
import argparse
import os
parser = argparse.ArgumentParser(description='create surfaces linking 2 close curves'+ 
'(GOCAD has this function but it does not work perfectly: overlapping triangles).' +
'This script works only if both curves are rectangles (side = constant coordinates).' + 
'This  script can be useful for limiting the amount of data to be processed in gocad:' + 
'example coarse mesh on the whole domain exept close to the fault')
parser.add_argument('input_file', help='pl file with 2 close curves')
parser.add_argument('output_file', help='gocad output file with the 4 surfaces written')
args = parser.parse_args()

inputfn = args.input_file
outputfn = args.output_file

########Read curves
fid=open(inputfn,'r')

lines = fid.readlines()
nline = len(lines)
readVTX = False
xyz1 = []
for i in range(nline):
   if lines[i].startswith('SEG'):
      icurrent = i
      break
   if lines[i].startswith('VRTX 1 '):
      readVTX = True
   if readVTX:
      vtxdata = lines[i].split()
      xyz1.append([float(xxx) for xxx in vtxdata[2:5]])
xyz1 = np.array(xyz1)

readVTX = False
xyz2 = []
for i in range(icurrent,nline):
   if readVTX & lines[i].startswith('SEG 1'):
      break
   if lines[i].startswith('VRTX 1 '):
      readVTX = True
   if readVTX:
      vtxdata = lines[i].split()
      xyz2.append([float(xxx) for xxx in vtxdata[2:5]])
xyz2 = np.array(xyz2)

##############create the surface inside the curves (4 parts)
xyz1 = xyz1[np.argsort(xyz1[:, 0])]
xyz2 = xyz2[np.argsort(xyz2[:, 0])]

min1 = min(xyz1[:,1])
i1 = np.where(xyz1[:,1]==min1)
min2 = min(xyz2[:,1])
i2 = np.where(xyz2[:,1]==min2)

fout = open('trash.dat','w')
np.savetxt(fout, xyz1[i1], fmt='%25.20f')
np.savetxt(fout, xyz2[i2], fmt='%25.20f')
fout.close()

os.system("python createGOCADTSurf.py trash.dat %s --axis 1 --large_precision --objectname trans_B " %args.output_file)

min1 = max(xyz1[:,1])
i1 = np.where(xyz1[:,1]==min1)
min2 = max(xyz2[:,1])
i2 = np.where(xyz2[:,1]==min2)

fout = open('trash.dat','w')
np.savetxt(fout, xyz1[i1], fmt='%25.20f')
np.savetxt(fout, xyz2[i2], fmt='%25.20f')
fout.close()

os.system("python createGOCADTSurf.py trash.dat %s --axis 1 --large_precision --objectname trans_U --append" %args.output_file)

xyz1 = xyz1[np.argsort(xyz1[:, 1])]
xyz2 = xyz2[np.argsort(xyz2[:, 1])]

min1 = min(xyz1[:,0])
i1 = np.where(xyz1[:,0]==min1)
min2 = min(xyz2[:,0])
i2 = np.where(xyz2[:,0]==min2)

fout = open('trash.dat','w')
np.savetxt(fout, xyz1[i1], fmt='%25.20f')
np.savetxt(fout, xyz2[i2], fmt='%25.20f')
fout.close()

os.system("python createGOCADTSurf.py trash.dat %s --axis 0 --large_precision --objectname trans_L --append" %args.output_file)

min1 = max(xyz1[:,0])
i1 = np.where(xyz1[:,0]==min1)
min2 = max(xyz2[:,0])
i2 = np.where(xyz2[:,0]==min2)

fout = open('trash.dat','w')
np.savetxt(fout, xyz1[i1], fmt='%25.20f')
np.savetxt(fout, xyz2[i2], fmt='%25.20f')
fout.close()

os.system("python createGOCADTSurf.py trash.dat %s --axis 0 --large_precision --objectname trans_R --append" %args.output_file)
os.remove('trash.dat')
