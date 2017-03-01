#Author Thomas Ulrich, LMU
#Convert dxf file (example GOCAD) to binary stl file (for SimModeler)
import numpy as np
import sys
import os
import struct

#list all the surface names
ntriangle=0
lsurfname=[]
###########READ DXF FUNCTIONS
def RW_SECTION():
	for i in range(0,2):
		fid.readline()
        if not args.isolate:
	   basename = os.path.basename(args.dxf_filename)

def RW_POINT():
	for i in range(0,8):
		fid.readline()

def RW_3DFACE():
	fid.readline()
	surfname = fid.readline().strip()
        if not surfname in lsurfname:
                 print "Reading surface : %s" %surfname
                 lsurfname.append(surfname)

	x=np.zeros((4,3))
	for i in range(0,4):
		for j in range(0,3):
			fid.readline()
			x[i,j]=fid.readline()

	normal = np.cross(x[1,:]-x[0,:],x[2,:]-x[0,:])
	norm=np.linalg.norm(normal)
        normal = normal/norm
        fout.write(struct.pack('<3f', *normal))
	for i in range(0,3):
		if args.proj!='':
                   xyz = pyproj.transform(lla, myproj, x[i,0],x[i,1], 1e3*x[i,2], radians=False)
                   fout.write(struct.pack('<3f', *xyz))

		else:
                   fout.write(struct.pack('<3f', *x[i,:]))
        if args.isolate:
           fout.write(struct.pack('<H', len(lsurfname)-1))
        else:
           fout.write(struct.pack('<H', 0))

def RW_ENDSEC():
        if not args.isolate:
	   basename = os.path.basename(args.dxf_filename)
	fid.close()
        fout.seek(80)
        print ntriangle
        fout.write(struct.pack('<L', ntriangle))
        fout.close()
	exit()

def RW_POLYLINE():
        fid.readline()
        print "Reading polyline: %s" %fid.readline().strip()
        for i in range(0,10):
                fid.readline()

def RW_VERTEX():
        for i in range(0,10):
                fid.readline()

def RW_SEQEND():
        return

options = {"SECTION": RW_SECTION, "POINT" : RW_POINT, "3DFACE": RW_3DFACE, "ENDSEC": RW_ENDSEC, "POLYLINE": RW_POLYLINE, "VERTEX": RW_VERTEX, "SEQEND": RW_SEQEND}
#########################

# parsing python arguments
import argparse
parser = argparse.ArgumentParser(description='Convert dxf file to bstl')
parser.add_argument('dxf_filename', help='dxf filename')
parser.add_argument('stl_filename', nargs='?',help='stl filname (if not used, stl_filename = dfxbasename.stl)',default='')
parser.add_argument('--proj', nargs=1, metavar=('projname'), default = (''), help='name of the projection (ex EPSG:32646 (UTM46N), or geocent (cartesian global)) if a projection is considered')
parser.add_argument('--isolate', dest='isolate', action='store_true', help='isolate every Gocad surface in a different stl solid (feature currently not supported by simModeler)')
args = parser.parse_args()

if args.stl_filename == '':
   args.stl_filename = args.dxf_filename[0:-4]+'.bstl'

fid = open(args.dxf_filename)
fout = open(args.stl_filename,'wb')
fout.seek(80)
#number of triangle. will be filled afterwards
fout.write(struct.pack('<L', 0))

if args.isolate:
   print "Isolating every Gocad surface in a different stl solid"

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

while fid.readline():  #0
	sId = fid.readline().strip() #eg 3DFACE
	
	if not sId in ["SECTION","POINT","3DFACE","ENDSEC","POLYLINE","VERTEX","SEQEND"]:
		print sId
		print("unknown format file")
		break
        if sId=='3DFACE':
           ntriangle=ntriangle+1
	options[sId]();

