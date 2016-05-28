import numpy as np

#Author: Thomas Ulrich, LMU 
#(inspired from a script by J. Klicpera)
#create surface from a structured grid of nodes

#Mesh example:

#3 - 3
#  /
#2 - 2
#  /
#1 - 1

import sys

if len(sys.argv)!=5:
	print 'Number of arguments (', len(sys.argv), ') not egal to 5'
	print 'usage: python *.py inputfilename objectname NX NY'
	print 'inputfile: x y z, one node coordinate by line'
	print 'objectname: name of the surface in gocad'
	print 'NX: number of nodes in the x direction'
        exit(-1)

inputfn = sys.argv[1]
objectname = sys.argv[2]
NX = int(sys.argv[3])
NY = int(sys.argv[4])

dataxyz = np.loadtxt(inputfn)
nvertex = np.shape(dataxyz)[0]

assert (NX*NY==nvertex), "NX*NY != nvertex %d x %d = %d != %d" %(NX,NY,NX*NY,nvertex)

triangles=[]

for j in range(NY-1):
   for i in range(1,NX):
      triangles.append([i+j*NX,i+1+j*NX,i+1+(j+1)*NX])
      triangles.append([i+j*NX,i+(j+1)*NX,i+1+(j+1)*NX])

print("GOCAD TSURF 1\nHEADER {\nname:"+objectname+"\n}\nTRIANGLES")

for i in range(0,nvertex):
   print("VRTX %d %f %f %f" %(i+1, dataxyz[i,0], dataxyz[i,1], dataxyz[i,2]))

for tr in triangles:
   print("TRGL %d %d %d" %(tr[0],tr[1],tr[2]))


print("END")

