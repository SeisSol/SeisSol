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

# parsing python arguments
import argparse
import os
parser = argparse.ArgumentParser(description='create surface from a ascii Esri grid file')
parser.add_argument('input_file', help='ascii Esri grid file')
parser.add_argument('output_file', help='gocad or stl output file')
parser.add_argument('--subsample', nargs=1, type=int, metavar=('onesample_every'), default = (1), help='use only one value every onesample_every in both direction')
parser.add_argument('--objectname', nargs=1, metavar=('objectname'), default = (''), help='name of the surface in gocad')
parser.add_argument('--hole', nargs=4, metavar=(('x0'),('x1'),('y0'),('y1')), default = (''), help='isolate a hole in surface defined by x0<=x<=x1 and y0<=y<=y1 (stl output only)')
parser.add_argument('--crop', nargs=4, metavar=(('x0'),('x1'),('y0'),('y1')), default = (''), help='select only surfaces in x0<=x<=x1 and y0<=y<=y1')
parser.add_argument('--proj', nargs=1, metavar=('projname'), default = (''), help='string describing the projection acting on the hole and crop coordinates (ex: +init=EPSG:32646 (UTM46N))')
args = parser.parse_args()

if args.objectname == '':
   base = os.path.basename(args.input_file)
   args.objectname = os.path.splitext(base)[0]
else:
   args.objectname = args.objectname[0]

if args.proj!='':
   print("Projecting the nodes coordinates")
   import mpl_toolkits.basemap.pyproj as pyproj
   lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
   if args.proj[0]!='geocent':
      sProj = args.proj[0]
      myproj=pyproj.Proj(sProj)
   else:
      myproj = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')

if args.hole!='':
   print("a hole will be isolated in the surface (stl only)")
   x0hole = float(args.hole[0])
   x1hole = float(args.hole[1])
   y0hole = float(args.hole[2])
   y1hole = float(args.hole[3])
   print("hole coordinates %f %f %f %f" %(x0hole,x1hole,y0hole,y1hole))
   if args.proj!='':
      xy = pyproj.transform(lla, myproj, x0hole, y0hole, 0, radians=False)
      x0hole = xy[0]
      y0hole = xy[1]
      xy = pyproj.transform(lla, myproj, x1hole, y1hole, 0, radians=False)
      x1hole = xy[0]
      y1hole = xy[1]
      print("hole coordinates (projected) %f %f %f %f" %(x0hole,x1hole,y0hole,y1hole))

fh = open(args.input_file)
NX = int(fh.readline().split()[1])
NY = int(fh.readline().split()[1])
xbotleft = float(fh.readline().split()[1])
ybotleft = float(fh.readline().split()[1])
dx = float(fh.readline().split()[1])
fh.readline()

lon =  np.arange(xbotleft, xbotleft + NX*dx, dx)
lat =  np.arange(ybotleft, ybotleft + NY*dx, dx)
#reverse array
lat = lat[::-1]
print(lon, lat)
print("reading elevation")
#using pandas rather than loadtxt because much faster
#elevation = np.loadtxt(fh)
import pandas as pd
elevation = pd.read_csv(fh, delimiter = " ", dtype=np.float64, header=None).values[:,:-1]
#elevation = pd.read_csv(fh, delimiter = " ", nrows=20, dtype=np.float64, header=None).values[:,:-1]
print(elevation)
print(elevation.shape)
print("done reading")
fh.close()


if args.crop!='':
   print("croping the surface")
   x0c = float(args.crop[0])
   x1c = float(args.crop[1])
   y0c = float(args.crop[2])
   y1c = float(args.crop[3])
   print("crop coordinates %f %f %f %f" %(x0c,x1c,y0c,y1c))
   if args.proj!='':
      xy = pyproj.transform(lla, myproj, x0c, y0c, 0, radians=False)
      x0c = xy[0]
      y0c = xy[1]
      xy = pyproj.transform(lla, myproj, x1c, y1c, 0, radians=False)
      x1c = xy[0]
      y1c = xy[1]
      print("crop coordinates (projected) %f %f %f %f" %(x0c,x1c,y0c,y1c))
   indexesXc = np.where((lon >= x0c) & (lon <= x1c))[0]
   indexesYc = np.where((lat >= y0c) & (lat <= y1c))[0]
   print(indexesXc)
   print(indexesXc.shape)
   print(indexesYc)
   print(indexesYc.shape)
   lon = lon[indexesXc]
   lat = lat[indexesYc]
   #elevation = elevation[indexesXc,indexesYc]
   elevation = elevation[:,indexesXc]
   elevation = elevation[indexesYc,:]

#subsampling
lat = lat[0::args.subsample[0]]
lon = lon[0::args.subsample[0]]
elevation =  elevation[0::args.subsample[0],0::args.subsample[0]]
print(elevation)

NY,NX = elevation.shape
print(NX,NY, lat.shape, lon.shape)

nnodes = NX*NY
ntriangles=2*(NX-1)*(NY-1)

nodes=np.zeros((nnodes+1,3))
triangles=np.zeros((ntriangles,3))

k=1
for j in range(0,NY):
    for i in range(0,NX):
       a = lon[i]
       b = lat[j]
       c = elevation[j,i]
       nodes[k,:]=  [lon[i], lat[j], elevation[j,i]]
       k=k+1
k=0
for j in range(NY-1):
   for i in range(1,NX):
      triangles[k,:] = [i+j*NX,i+1+j*NX,i+1+(j+1)*NX]
      triangles[k+1,:] = [i+j*NX,i+(j+1)*NX,i+1+(j+1)*NX]
      k=k+2

solid_id = np.zeros(ntriangles)
if args.hole!='':
   for k in range(ntriangles):
      xmin = nodes[triangles[k,0],0]
      xmax = nodes[triangles[k,2],0]
      ymin = nodes[triangles[k,0],1]
      ymax = nodes[triangles[k,2],1]
      if  ((xmin>x0hole) & (xmax<x1hole))&((ymin>y0hole) & (ymax<y1hole)):
         solid_id[k]=1
      else:
         solid_id[k]=0
nsolid=int(max(solid_id))

_, ext = os.path.splitext(args.output_file)


if ext=='.ts':
   fout = open(args.output_file,'w')
   fout.write("GOCAD TSURF 1\nHEADER {\nname:"+args.objectname+"\n}\nTRIANGLES\n")
   for k in range(1,nnodes+1):
      fout.write("VRTX %d %f %f %f\n" %(k, nodes[k,0], nodes[k,1], nodes[k,2]))

   for k in range(ntriangles):
      fout.write("TRGL %d %d %d\n" %(triangles[k,0],triangles[k,1],triangles[k,2]))
   fout.write("END")
elif ext=='.stl':
   fout = open(args.output_file,'w')
   for sid in range(nsolid+1):
      fout.write("solid %s%d\n" %(args.objectname, sid))
      idtr = np.where(solid_id==sid)[0]
      for k in idtr:
         normal = np.cross(nodes[triangles[k,1],:]-nodes[triangles[k,0],:],nodes[triangles[k,2],:]-nodes[triangles[k,0],:])
      norm=np.linalg.norm(normal)
      fout.write('facet normal %e %e %e\n' %tuple(normal/norm))
      fout.write('outer loop\n')
      for i in range(0,3):
         fout.write('vertex %.10e %.10e %.10e\n' % tuple(nodes[triangles[k,i],:]))
      fout.write("endloop\n")
      fout.write("endfacet\n")
      fout.write("endsolid %s%d\n" %(args.objectname, sid))
elif ext=='.bstl':
   import struct
   fout = open(args.output_file,'wb')
   fout.seek(80)
   fout.write(struct.pack('<L', ntriangles))
   for sid in range(nsolid+1):
      idtr = np.where(solid_id==sid)[0]
      for k in idtr:
         normal = np.cross(nodes[triangles[k,1],:]-nodes[triangles[k,0],:],nodes[triangles[k,2],:]-nodes[triangles[k,0],:])
         norm=np.linalg.norm(normal)
         normal = normal/norm
         fout.write(struct.pack('<3f', *normal))
         for i in range(0,3):
            fout.write(struct.pack('<3f', *nodes[triangles[k,i],:]))
         fout.write(struct.pack('<H', sid))  
else:
   print("only bsl, stl and ts are valid output formats")

fout.close()
