###
# Author: Thomas Ulrich, LMU, 24.09.2015
# Read a coastline file from GMT
# Create a vtk file from that input
# aim: diplaying the coastline with the simulation results for instance

# parsing python arguments
import argparse
import os
import numpy as np

parser = argparse.ArgumentParser(description="create surface from a structured grid of nodes")
parser.add_argument("--lon", nargs=2, metavar=(("lonmin"), ("lonmax")), default=(""), help="lonmin: minimum longitude, lonmax: maximum longitude", type=float)
parser.add_argument("--lat", nargs=2, metavar=(("latmin"), ("latmax")), default=(""), help="latmin: minimum latitude, lonmax: maximum latitude", type=float)
parser.add_argument("--proj", nargs=1, metavar=("projname"), default=(""), help="name of the projection (ex +init=EPSG:32646 (UTM46N), or geocent (cartesian global)) if a projection is considered")
parser.add_argument("--resolution", nargs=1, metavar=("resolution"), default=("i"), help="resolution of the coastline,  (f)ull, (h)igh, (i)ntermediate, (l)ow, and (c)rude")
parser.add_argument("--coords", nargs=3, metavar=("x", "y", "z"), default=([0, 0, 0]), help="translation vector", type=float)
args = parser.parse_args()

if args.proj != "":
    import pyproj

    lla = pyproj.Proj(proj="latlong", ellps="WGS84", datum="WGS84")
    myproj = pyproj.Proj(args.proj[0])
    print("using pyproj to project the coordinates...Please check that the projection used corresponds with your lat/lon range")

if args.lon == "":
    print("no longitude range specified")
    exit()
else:
    lonmin, lonmax = args.lon[0], args.lon[1]

if args.lat == "":
    print("no longitude range specified")
    exit()
else:
    latmin, latmax = args.lat[0], args.lat[1]

# export cordinates from GMT
command = "module load gmt;gmt pscoast -R%f/%f/%f/%f -D%s -M -W > coastline.dat" % (lonmin, lonmax, latmin, latmax, args.resolution[0])
os.system(command)

# Read GMT file
vertices = []
segments = []
nvert = 0
newPolyLine = True
fid = open("coastline.dat")
for line in fid:
    if line.startswith("#"):
        continue
    if line.startswith(">"):
        newPolyLine = True
    else:
        vertices.append([float(val) for val in line.split()])
        nvert = nvert + 1
        if newPolyLine == False:
            segments.append([nvert - 1, nvert])
        newPolyLine = False
fid.close()

vertices = np.asarray(vertices)
segments = np.asarray(segments)

# Now write vtk file
fout = open("CoastLine.vtk", "w")
nlines = 0

fout.write(
    "# vtk DataFile Version 2.0\n\
parabola - polyline\n\
ASCII\n\n\
DATASET POLYDATA\n"
)
nvert = vertices.shape[0]
fout.write("POINTS      %d float\n" % (nvert))

if args.proj != "":
    x, y = pyproj.transform(lla, myproj, vertices[:, 0], vertices[:, 1], radians=False)
    xy = np.vstack((x, y)).T
else:
    xy = vertices
xy = xy - args.coords[0:2]
xyz = np.zeros((nvert, 3))
xyz[:, 0:2] = xy
xyz[:, 2] = args.coords[2]

np.savetxt(fout, xyz, "%e %e %e")

nseg = segments.shape[0]
fout.write("\nLINES %d %d\n" % (nseg, 3 * nseg))
np.savetxt(fout, segments - 1, "2 %d %d")

fout.write("\n")
fout.close()

print("CoastLine.vtk successfully created")
