#!/usr/bin/env python3
import argparse
import subprocess
import numpy as np

parser = argparse.ArgumentParser(description="create coastline vtk file from gmt")
parser.add_argument("--lon", nargs=2, metavar=(("lonmin"), ("lonmax")), help="lonmin: minimum longitude, lonmax: maximum longitude", type=float)
parser.add_argument("--lat", nargs=2, metavar=(("latmin"), ("latmax")), help="latmin: minimum latitude, lonmax: maximum latitude", type=float)
parser.add_argument("--proj", nargs=1, metavar=("projname"), help="project the data. projname: PROJ string describing the projection (ex epsg:32646 for UTM46N). Use geocent for geocentric cartesian")
parser.add_argument("--resolution", nargs=1, metavar=("resolution"), default=("i"), help="Coastline resolution,  (f)ull, (h)igh, (i)ntermediate, (l)ow, and (c)rude")
parser.add_argument("--recenter", nargs=2, metavar=("x", "y"), default=([0, 0]), help="translate coordinate array (e.g. x_new = x_old - x)", type=float)
parser.add_argument("--z", nargs=1, metavar=("elevation"), default=([0]), help="z coordinate of coastline", type=float)
parser.add_argument("--filter_by_area", nargs=1, metavar=("min_area"), help="remove closed feature of area smaller than min_area (-A option of pscoast)", type=float)
args = parser.parse_args()

command = f"gmt pscoast -R{args.lon[0]}/{args.lon[1]}/{args.lat[0]}/{args.lat[1]} -D{args.resolution[0]} -M -W"
if args.filter_by_area:
    command += f" -A{args.filter_by_area[0]}"

# export cordinates from GMT
result = subprocess.run(command.split(), capture_output=True, text=True)

# Read GMT file
xyz = []
segments = []
nvert = 0
newPolyLine = True
for line in result.stdout.split('\n'):
    if not line:
        continue
    if line.startswith("#"):
        continue
    if line.startswith(">"):
        newPolyLine = True
    else:
        xyz.append([float(val) for val in line.split()])
        nvert = nvert + 1
        if not newPolyLine:
            segments.append([nvert - 1, nvert])
        newPolyLine = False

xyz = np.asarray(xyz)
# add extra column for z coordinates
xyz = np.insert(xyz, xyz.shape[1], 0, axis=1)

segments = np.asarray(segments)

if args.proj:
    from pyproj import Transformer

    f = lambda x: {"proj": "geocent", "ellps": "WGS84", "datum": "WGS84"} if x == "geocent" else x
    transformer = Transformer.from_crs("epsg:4326", f(args.proj[0]), always_xy=True)
    xyz[:, 0], xyz[:, 1], xyz[:, 2] = transformer.transform(xyz[:, 0], xyz[:, 1], xyz[:, 2])
    if args.proj[0] != "geocent":
        xyz[:, 2] = args.z[0]
else:
    xyz[:, 2] = args.z[0]

xyz[:, 0] -= args.recenter[0]
xyz[:, 1] -= args.recenter[1]

nvert = xyz.shape[0]
nseg = segments.shape[0]


# Now write vtk file
with open("CoastLine.vtk", "w") as fout:
    fout.write(
        f"""# vtk DataFile Version 2.0
parabola - polyline
ASCII
DATASET POLYDATA
POINTS {nvert} float
"""
    )
    np.savetxt(fout, xyz, "%e %e %e")
    fout.write(f"\nLINES {nseg} {3*nseg}\n")
    np.savetxt(fout, segments - 1, "2 %d %d")
    fout.write("\n")

print("CoastLine.vtk successfully created")
