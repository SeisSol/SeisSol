#Author Thomas Ulrich, LMU
#Convert dxf file (example GOCAD) to stl file (for SimModeler)
import numpy as np
import sys
import os

###########READ DXF FUNCTIONS
def RW_SECTION():
	for i in range(0,2):
		fid.readline()
	basename = os.path.basename(args.dxf_filename)
	fout.write("solid %s\n" %(basename))

def RW_POINT():
	for i in range(0,8):
		fid.readline()

def RW_3DFACE():
	fid.readline()
	fid.readline()
	x=np.zeros((4,3))
	for i in range(0,4):
		for j in range(0,3):
			fid.readline()
			x[i,j]=fid.readline()

	normal = np.cross(x[1,:]-x[0,:],x[2,:]-x[0,:])
	norm=np.linalg.norm(normal)
	fout.write('facet normal %e %e %e\n' %tuple(normal/norm))
	fout.write('outer loop\n')
	for i in range(0,3):
		if args.proj!='':
                   xy=myproj(x[i,0],x[i,1])
                   fout.write('vertex %.10e %.10e %.10e\n' %(xy[0],xy[1],1e3*x[i,2]))
		else:
		   fout.write('vertex %.10e %.10e %.10e\n' %tuple(x[i,:]))
	fout.write("endloop\n")
	fout.write("endfacet\n")

def RW_ENDSEC():
	basename = os.path.basename(args.dxf_filename)
	fout.write("endsolid %s\n" %(basename))
	fout.close()
	fid.close()
	exit()

options = {"SECTION": RW_SECTION, "POINT" : RW_POINT, "3DFACE": RW_3DFACE, "ENDSEC": RW_ENDSEC}
#########################

# parsing python arguments
import argparse
parser = argparse.ArgumentParser(description='Convert dxf file to stl')
parser.add_argument('dxf_filename', help='dxf filename')
parser.add_argument('stl_filename', nargs='?',help='stl filname (if not used, stl_filename = dfxbasename.stl)',default='')
parser.add_argument('--proj', nargs=1, metavar=('projname'), default = (''), help='name of the projection (ex EPSG:32646 (UTM46N)) if a projection is considered')

args = parser.parse_args()

if args.stl_filename == '':
   args.stl_filename = args.dxf_filename[0:-4]+'.stl'

fid = open(args.dxf_filename)
fout = open(args.stl_filename,'w')

if args.proj!='':
   print "projecting the nodes coordinates"
   import mpl_toolkits.basemap.pyproj as pyproj
   sProj = "+init=%s" %args.proj[0]
   myproj=pyproj.Proj(sProj)
else:
   print "no projection carried out"

while fid.readline():  #0
	sId = fid.readline().strip() #eg 3DFACE
	
	if not sId in ["SECTION","POINT","3DFACE","ENDSEC"]:
		print sId
		print("unknown format file")
		break

	options[sId]();


