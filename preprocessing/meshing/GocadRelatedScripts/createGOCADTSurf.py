import numpy as np
#Author: Thomas Ulrich, LMU 
#(inspired from a script by J. Klicpera)
#create surface from a structured lat/lon/depth points set
#structured in the sense:  the point set should be serveral lines of constant latitudes (i.e. Y coordinate)

#Mesh example:

#    5   
#  /
# /  4
#  /
#3 - 3
#  /
#2 - 2
#  /
#1 - 1

import sys

if len(sys.argv)!=3:
	print 'Number of arguments (', len(sys.argv), ') not egal to 3'
	print 'usage: python *.py inputfilename objectname > outputfilename'
	print 'inputfile: x y z, one node coordinate by line'
	print 'objectname: name of the surface in gocad'
        exit(-1)

inputfn = sys.argv[1]
objectname = sys.argv[2]

dataxyz = np.loadtxt(inputfn)
nvertex = np.shape(dataxyz)[0]


lats= set(dataxyz[:,1])
lats=sorted(lats)
#precision for selection the nodes with a given latitude
dx=1e-4
triangles=[]

for i, lat in enumerate(lats):
   if i==0:
      p1 = np.where(abs(dataxyz[:,1] -lat)<dx)[0]
      np1 = np.size(p1)
      continue
   p0 = p1
   np0 = np1
   p1 = np.where(abs(dataxyz[:,1] -lat)<dx)[0]
   np1 = np.size(p1)
   if (np0<np1):
   	for i0 in range(np0-1):
           triangles.append([ p0[i0],p0[i0+1],p1[i0] ])
           triangles.append([ p0[i0+1],p1[i0+1],p1[i0] ])
   	for i0 in range(np0-1, np1-1):
           triangles.append([ p0[np0-1],p1[i0+1],p1[i0] ])
   else:
   	for i1 in range(np1-1):
           triangles.append([ p1[i1],p1[i1+1],p0[i1] ])
           triangles.append([ p1[i1+1],p0[i1+1],p0[i1] ])
   	for i0 in range(np1-1, np0-1):
           triangles.append([ p1[np1-1],p0[i1+1],p0[i1] ])


print("GOCAD TSURF 1\nHEADER {\nname:"+objectname+"\n}\nTRIANGLES")

for i in range(0,nvertex):
   print("VRTX %d %f %f %f" %(i+1, dataxyz[i,0], dataxyz[i,1], dataxyz[i,2]))

for tr in triangles:
   print("TRGL %d %d %d" %(tr[0]+1,tr[1]+1,tr[2]+1))


print("END")

