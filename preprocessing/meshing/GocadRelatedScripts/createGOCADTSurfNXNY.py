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
parser.add_argument('--NX', nargs=1, metavar=('NX'), default = (''), help='NX: number of nodes in the first structured dimension')
parser.add_argument('--subsample', nargs=1, metavar=('onesample_every'), default = (''), help='use only one value every onesample_every in both direction')
parser.add_argument('--objectname', nargs=1, metavar=('objectname'), default = (''), help='name of the surface in gocad')
parser.add_argument('--hole', nargs=4, metavar=(('x0'),('x1'),('y0'),('y1')), default = (''), help='create a hole in surface defined by x0<=x<=x1 and y0<=y<=y1')
parser.add_argument('--crop', nargs=4, metavar=(('x0'),('x1'),('y0'),('y1')), default = (''), help='select only surfaces in x0<=x<=x1 and y0<=y<=y1')
parser.add_argument('--proj', nargs=1, metavar=('projname'), default = (''), help='string describing its projection (ex: +init=EPSG:32646 (UTM46N), or geocent (cartesian global)) if a projection is considered')
args = parser.parse_args()

if args.objectname == '':
   base = os.path.basename(args.input_file)
   args.objectname = os.path.splitext(base)[0]
else:
   args.objectname = args.objectname[0]

dataxyz = np.loadtxt(args.input_file)
nvertex = np.shape(dataxyz)[0]

if args.crop!='':
   if args.NX!='':
      print("uncompatible inputs")
      exit()
   print("croping the surface")
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


if args.NX=='':
   print("NX not defined: trying to guess it...")
   rowdiff = dataxyz[1,:]-dataxyz[0,:]
   ix = -1
   ids = np.where(abs(rowdiff)<1e-16)[0]
   if len(ids)==1:
      #only one column starts with constant values
      ix = ids[0]
      for i in range(1,nvertex):
         if abs(dataxyz[i,ix]-dataxyz[i-1,ix])>1e-16:
            NX = i
            assert (nvertex%NX==0), "nvertex%%NX!=0 nvertex/NX = %f" %(float(nvertex)/NX)
            NY = nvertex/NX
            print("NX,NY = %d,%d" %(NX,NY))
            break
   elif len(ids)>1:
      print("2 columns starts with constant values")
      nx =[]
      #find other dimension
      for ix in range(0,3):
         if ix not in ids:
            iy=ix
      for ix in ids:
         for i in range(1,nvertex):
            if abs(dataxyz[i,ix]-dataxyz[i-1,ix])>1e-16:
               break
         if abs(dataxyz[0,iy]-dataxyz[0+i,iy])>1e-16:
            nx.append(1e20)
         else:
            nx.append(i)
      NX = min(nx)

      if NX==1e10:
         print("unable to guess NX and NY")
         exit()

      assert (nvertex%NX==0), "nvertex%%NX!=0 nvertex/NX = %f" %(float(nvertex)/NX)
      NY = nvertex/NX
      print("NX,NY = %d,%d" %(NX,NY))
   else:
      print("unable to guess NX and NY")
      exit()
else:
   print("using user defined NX")
   NX=int(args.NX[0])
   assert (nvertex%NX==0), "nvertex%%NX!=0 nvertex/NX = %f" %(float(nvertex)/NX)
   NY = int(nvertex/NX)

if args.hole!='':
   print("a hole will be left in the surface")
   x0hole = float(args.hole[0])
   x1hole = float(args.hole[1])
   y0hole = float(args.hole[2])
   y1hole = float(args.hole[3])
   print("hole coordinates %f %f %f %f" %(x0hole,x1hole,y0hole,y1hole))

if args.subsample!='':
   onesample_every = int(args.subsample[0])
   print("subsampling : 1/%d" %onesample_every)
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
         triangles.append([i+j*NX,i+1+(j+1)*NX,i+(j+1)*NX])

if args.proj!='':
   print("Projecting the nodes coordinates")
   import mpl_toolkits.basemap.pyproj as pyproj
   lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
   if args.proj[0]!='geocent':
      sProj = args.proj[0]
      myproj=pyproj.Proj(sProj)
   else:
      myproj = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
else:
   print("no projection carried out")



fout = open(args.output_file,'w')
### WRITE THE GOCAD TS FILE
fout.write("GOCAD TSURF 1\nHEADER {\nname:"+args.objectname+"\n}\nTRIANGLES\n")
for j in range(0,NY):
    for i in range(0,NX):
        if args.proj!='':
            xyz = pyproj.transform(lla, myproj, dataxyz[i,j,0],dataxyz[i,j,1], 1e3*dataxyz[i,j,2], radians=False)
            fout.write('VRTX '+str(i+j*NX+1)+' %.10e %.10e %.10e\n' %tuple(xyz))
        else:
            fout.write("VRTX %d %f %f %f\n" %(i+j*NX+1, dataxyz[i,j,0], dataxyz[i,j,1], dataxyz[i,j,2]))
for tr in triangles:
   fout.write("TRGL %d %d %d\n" %(tr[0],tr[1],tr[2]))
fout.write("END")

fout.close()
