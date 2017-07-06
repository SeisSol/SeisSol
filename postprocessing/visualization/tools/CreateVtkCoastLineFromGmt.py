###
#Author: Thomas Ulrich, LMU, 24.09.2015
#Read a coastline file from GMT
#Create a vtk file from that input
#aim: diplaying the coastline with the simulation results for instance

# parsing python arguments
import argparse
import os
parser = argparse.ArgumentParser(description='create surface from a structured grid of nodes')
parser.add_argument('--lon', nargs=2, metavar=(('lonmin'),('lonmax')), default = (''), help='lonmin: minimum longitude, lonmax: maximum longitude')
parser.add_argument('--lat', nargs=2, metavar=(('latmin'),('latmax')), default = (''), help='latmin: minimum latitude, lonmax: maximum latitude')
parser.add_argument('--proj', nargs=1, metavar=('projname'), default = (''), help='name of the projection (ex +init=EPSG:32646 (UTM46N), or geocent (cartesian global)) if a projection is considered')
parser.add_argument('--resolution', nargs=1, metavar=('resolution'), default = ('i'), help='resolution of the coastline,  (f)ull, (h)igh, (i)ntermediate, (l)ow, and (c)rude')
args = parser.parse_args()

if args.proj!='':
   import mpl_toolkits.basemap.pyproj as pyproj
   lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
   sProj = args.proj[0]
   myproj=pyproj.Proj(sProj)
   print("using pyproj to project the coordinates...Please check that the projection used corresponds with your lat/lon range") 

if args.lon=='':
   print "no longitude range specified"
   exit()
else:
   lonmin = float(args.lon[0])
   lonmax = float(args.lon[1])

if args.lat=='':
   print "no longitude range specified"
   exit()
else:
   latmin = float(args.lat[0])
   latmax = float(args.lat[1])

#export cordinates from GMT
command = "module load gmt;gmt pscoast -R%f/%f/%f/%f -D%s -M -W > coastline.dat" %(lonmin, lonmax, latmin, latmax, args.resolution[0])
os.system(command)

#Read GMT file
vertices=[]
segments=[]
nvert=0
newPolyLine=True
fid=open('coastline.dat')
for line in fid:
   if (line.startswith('#')):
      continue
   if (line.startswith('>')):
      newPolyLine=True
   else:
      vertices.append(line.split())
      nvert=nvert+1
      if newPolyLine==False:
         segments.append([nvert-1,nvert])
      newPolyLine=False
fid.close()


#Now write vtk file
fout=open('CoastLine.vtk','w')
nlines=0

fout.write('# vtk DataFile Version 2.0\n\
parabola - polyline\n\
ASCII\n\n\
DATASET POLYDATA\n')


fout.write('POINTS      %d float\n' %(len(vertices)+1))
#dummy point for avoiding changing numbering scheme
fout.write('0 0 0\n')
for vert in vertices:
   if args.proj!='':
      latlon=[float(v) for v in vert[0:2]]
      xyz = pyproj.transform(lla, myproj,latlon[0],latlon[1],0, radians=False)
      fout.write('%e %e %e\n' %tuple(xyz))
      #fout.write('%e %e 0.\n' %(xyz[0],xyz[2]))
   else:
      fout.write('%s %s 0.\n' %tuple(vert[0:2]))


fout.write('\nLINES %d %d\n' %(len(segments),3*len(segments)))
for seg in segments:
   fout.write('2 %s %s\n' %tuple(seg))

fout.write('\n')
fout.close()

print("CoastLine.vtk successfully created")
