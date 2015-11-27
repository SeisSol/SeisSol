#Author Thomas Ulrich, LMU
#Convert dxf file (example GOCAD) to stl file (for SimModeler)
import numpy as np
import sys

UtmProjection = False

if UtmProjection==True:
   print "projecting the nodes coordinates"
   import mpl_toolkits.basemap.pyproj as pyproj
   UTM46N=pyproj.Proj("+init=EPSG:32646")
else:
   print "no projection carried out"

###########READ DXF FUNCTIONS
def RW_SECTION():
	for i in range(0,2):
		fid.readline()
	fout.write("solid AssimpScene\n")

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
		if UtmProjection==True:
                   xy=UTM46N(x[i,0]/1e4,x[i,1]/1e4)
                   fout.write('vertex %.10e %.10e %.10e\n' %(xy[0],xy[1],x[i,2]))
		else:
		   fout.write('vertex %.10e %.10e %.10e\n' %tuple(x[i,:]))
	fout.write("endloop\n")
	fout.write("endfacet\n")

def RW_ENDSEC():
	fout.write("endsolid AssimpScene\n")
	fout.close()
	fid.close()
	exit()

options = {"SECTION": RW_SECTION, "POINT" : RW_POINT, "3DFACE": RW_3DFACE, "ENDSEC": RW_ENDSEC}
#########################

narg = len(sys.argv)
if (narg <2) |  (narg>4):
	print("usage (1): python convertDxf2Stl.py input.dxf output.stl")
	print("usage (2): python convertDxf2Stl.py input.dxf")
	print("usage (2): then outputfn = input.stl")
	print("UtmProjection (bool) : hardcoded")
	exit()
elif narg==2:
   fid=open(sys.argv[1])
   fnout = sys.argv[1][0:-4]+'.stl'
   fout=open(fnout,"w")
elif narg==3:
   fid=open(sys.argv[1])
   fout=open(sys.argv[2],"w")


while fid.readline():  #0
	sId = fid.readline().strip() #eg 3DFACE
	
	if not sId in ["SECTION","POINT","3DFACE","ENDSEC"]:
		print sId
		print("unknown format file")
		break

	options[sId]();


