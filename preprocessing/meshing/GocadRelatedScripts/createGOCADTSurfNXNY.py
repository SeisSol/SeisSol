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

# parsing python arguments
import argparse
import os
parser = argparse.ArgumentParser(description='create surface from a structured grid of nodes')
parser.add_argument('input_file', help='x y z, one node coordinate by line')
parser.add_argument('output_file', help='gocad output file')
parser.add_argument('--NXNY', nargs=2, metavar=(('NX'),('NY')), default = (''), help='NX: number of nodes in the x direction, NY idem y')
parser.add_argument('--subsample', nargs=1, metavar=('onesample_every'), default = (''), help='use only one value every onesample_every in both direction')
parser.add_argument('--objectname', nargs=1, metavar=('objectname'), default = (''), help='name of the surface in gocad')
parser.add_argument('--hole', nargs=4, metavar=(('x0'),('x1'),('y0'),('y1')), default = (''), help='create a hole in surface defined by x0<=x<=x1 and y0<=y<=y1')
parser.add_argument('--crop', nargs=4, metavar=(('x0'),('x1'),('y0'),('y1')), default = (''), help='select only surfaces in x0<=x<=x1 and y0<=y<=y1')
parser.add_argument('--proj', nargs=1, metavar=('projname'), default = (''), help='name of the projection (ex EPSG:32646 (UTM46N), or geocent (cartesian global)) if a projection is considered')
args = parser.parse_args()

if args.objectname == '':
   base = os.path.basename(args.input_file)
   args.objectname = os.path.splitext(base)[0]
else:
   args.objectname = args.objectname[0]

dataxyz = np.loadtxt(args.input_file)
nvertex = np.shape(dataxyz)[0]

if args.crop!='':
   if args.NXNY!='':
      print("uncompatible inputs")
      exit()
   print "croping the surface"
   x0c = float(args.crop[0])
   x1c = float(args.crop[1])
   y0c = float(args.crop[2])
   y1c = float(args.crop[3])
   indexes = np.where((dataxyz[:,0] >= x0c) & (dataxyz[:,0] <= x1c) & (dataxyz[:,1] >= y0c) & (dataxyz[:,1] <= y1c))
   dataxyz = dataxyz[indexes[0],:]
   nvertex = np.shape(dataxyz)[0]
   #print test
   #indexes = np.where((test >= x0c))
   # & (dataxyz(:,1) >= y0) & (dataxyz(:,1) <= y1))
   #print indexes


if args.NXNY=='':
   print "NX and NY not defined: trying to guess them..."
   if dataxyz[1,1]==dataxyz[0,1]:
      #we can try to get NX
      for i in range(1,nvertex):
         if dataxyz[i,1]!=dataxyz[i-1,1]:
            NX = i
            NY = nvertex/NX
            break
      print ((NX,NY))
   elif dataxyz[1,0]==dataxyz[0,0]:
      #we can try to get NY
      for i in range(1,nvertex):
         if dataxyz[i,0]!=dataxyz[i-1,0]:
            NX = i
            NY = nvertex/NX
            break
      print ((NX,NY))
   else:
      print "unable to guess NX and NY"
      exit()
else:
   print "using user defined NX and NY"
   NX=int(args.NXNY[0])
   NY=int(args.NXNY[1])
   assert (NX*NY==nvertex), "NX*NY != nvertex %d x %d = %d != %d" %(NX,NY,NX*NY,nvertex)

if args.hole!='':
   print "a hole will be left in the surface"
   x0hole = float(args.hole[0])
   x1hole = float(args.hole[1])
   y0hole = float(args.hole[2])
   y1hole = float(args.hole[3])
   print "hole coordinates %f %f %f %f" %(x0hole,x1hole,y0hole,y1hole)

if args.subsample!='':
   onesample_every = int(args.subsample[0])
   print "subsampling : 1/%d" %onesample_every
else:
   onesample_every = 1

dataxyz = dataxyz.reshape((NY,NX,3))

triangles=[]

dataxyz = dataxyz[::onesample_every,::onesample_every,:]
NX = np.shape(dataxyz)[0]
NY = np.shape(dataxyz)[1]
nvertex = NX*NY


for j in range(NY-1):
   for i in range(1,NX):
      write_triangle=True
      if args.hole!='':
         for ij in [[i-1,j-1], [i-1,j],[i,j-1],[i,j]]:
            if  ((dataxyz[ij[0],ij[1],0]>x0hole) & (dataxyz[ij[0],ij[1],0]<x1hole))&((dataxyz[ij[0],ij[1],1]>y0hole) & (dataxyz[ij[0],ij[1],1]<y1hole)):
               write_triangle=False
      if write_triangle:
         triangles.append([i+j*NX,i+1+j*NX,i+1+(j+1)*NX])
         triangles.append([i+j*NX,i+(j+1)*NX,i+1+(j+1)*NX])

if args.proj!='':
   print "Projecting the nodes coordinates"
   import mpl_toolkits.basemap.pyproj as pyproj
   lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
   if args.proj[0]!='geocent':
      sProj = "+init=%s" %args.proj[0]
      myproj=pyproj.Proj(sProj)
   else:
      myproj = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
else:
   print "no projection carried out"



fout = open(args.output_file,'w')
### WRITE THE GOCAD TS FILE
fout.write("GOCAD TSURF 1\nHEADER {\nname:"+args.objectname+"\n}\nTRIANGLES\n")
for j in range(0,NY):
    for i in range(0,NX):
        if args.proj!='':
            xyz = pyproj.transform(lla, myproj, dataxyz[i,j,0],dataxyz[i,j,1], dataxyz[i,j,2], radians=False)
            fout.write('VRTX '+str(i+j*NX+1)+' %.10e %.10e %.10e\n' %tuple(xyz))
        else:
            fout.write("VRTX %d %f %f %f\n" %(i+j*NX+1, dataxyz[i,j,0], dataxyz[i,j,1], dataxyz[i,j,2]))
for tr in triangles:
   fout.write("TRGL %d %d %d\n" %(tr[0],tr[1],tr[2]))
fout.write("END")

fout.close()
