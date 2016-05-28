###
#Author: Thomas Ulrich, LMU, 24.09.2015
#Read a coastline file from GMT
#Create a vtk file from that input
#aim: diplaying the coastline with the simulation results for instance
projectlatlon = True
lonmin=90
lonmax=100
latmin=-5
latmax=15
###

if projectlatlon:
   import mpl_toolkits.basemap.pyproj as pyproj
   lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
   sProj = "+init=%s" %"EPSG:32646"
   myproj=pyproj.Proj(sProj)
   #myproj = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
   print("using pyproj to project the coordinates...Please check that the projection used corresponds with your lat/lon range") 

#export cordinates from GMT
import os
command = "module load gmt;gmt pscoast -R%f/%f/%f/%f -Di -M -W > coastline.dat" %(lonmin, lonmax, latmin, latmax)
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
   if projectlatlon:
      latlon=[float(v) for v in vert[0:2]]
      xyz = pyproj.transform(lla, myproj,latlon[0],latlon[1],0, radians=False)
      #fout.write('%e %e %e\n' %tuple(xyz))
      fout.write('%e %e 0.\n' %(xyz[0],xyz[2]))
   else:
      fout.write('%s %s 0.\n' %tuple(vert[0:2]))


fout.write('\nLINES %d %d\n' %(len(segments),3*len(segments)))
for seg in segments:
   fout.write('2 %s %s\n' %tuple(seg))

fout.write('\n')
fout.close()

print("CoastLine.vtk successfully created")
