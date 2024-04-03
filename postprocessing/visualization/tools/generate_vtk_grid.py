#!/usr/bin/env python3
import pyproj
import numpy as np
import argparse
import seissolxdmf as sx
from pyproj import Transformer

parser = argparse.ArgumentParser(
    description="generate lat/lon grid for paraview vizualisation"
)
parser.add_argument(
    "--dx", nargs=1, help="sampling between grid lines", type=float, required=True
)
parser.add_argument(
    "--lon",
    nargs=2,
    metavar=(("lonmin"), ("lonmax")),
    help="minimum (lonmin) and maximum (lonmax) longitude",
    type=float,
)
parser.add_argument(
    "--lat",
    nargs=2,
    metavar=(("latmin"), ("latmax")),
    help="minimum (latmin) and maximum (latmax) latitude",
    type=float,
    type=float,
)

parser.add_argument(
    "--from_fault_output",
    nargs=1,
    metavar=("filename added_degree"),
    help="infer domain from range of fault xdmf + added extension in degree",
)
parser.add_argument(
    "--output_prefix", nargs=1, help="prefix if output vtk", default=["grid"]
)

parser.add_argument(
    "--lineout",
    nargs=2,
    metavar=(("lon2"), ("lat2")),
    default=("0.35", "0.35"),
    help="draw line on out from the grid of length lon2 (or lat2)",
)
parser.add_argument(
    "--proj",
    nargs=1,
    metavar=("projname"),
    help="project the data. projname: proj4 string describing the projection",
)
args = parser.parse_args()

if args.from_fault_output:
    assert (not args.lon) and (not args.lat)
    fn, added_deg = args.from_fault_output[0].split()
    s = sx.seissolxdmf(fn)
    added_deg = float(added_deg)
    xyz = s.ReadGeometry()
    assert args.proj
    transformer = Transformer.from_crs(args.proj[0], "epsg:4326", always_xy=True)
    xyz[:, 0], xyz[:, 1], xyz[:, 2] = transformer.transform(
        xyz[:, 0], xyz[:, 1], xyz[:, 2]
    )
    args.lon = [0, 0]
    args.lat = [0, 0]
    args.lon[1], args.lat[1], _ = np.amax(xyz, axis=0) + added_deg
    args.lon[0], args.lat[0], _ = np.amin(xyz, axis=0) - added_deg

    def my_round(x, dx, floor=True):
        func = np.floor if floor else np.ceil
        return func(x / dx) * dx

    args.lon[1], args.lat[1] = [
        my_round(x, args.dx[0], False) for x in [args.lon[1], args.lat[1]]
    ]
    args.lon[0], args.lat[0] = [
        my_round(x, args.dx[0], True) for x in [args.lon[0], args.lat[0]]
    ]

lons = np.arange(args.lon[0], args.lon[1], args.dx[0])
lats = np.arange(args.lat[0], args.lat[1], args.dx[0])
print("lon", lons)
print("lat", lats)

nx = lons.shape[0]
ny = lats.shape[0]

nodes = []
seg = []

k = 0
dx0 = float(args.lineout[0])
dy0 = float(args.lineout[1])

for i in range(nx):
    nodes.append([lons[i], lats.min() - dy0])
    nodes.append([lons[i], lats.max() + dy0])
    k = len(nodes)
    seg.append([k - 1, k])

for j in range(ny):
    nodes.append([lons.min() - dx0, lats[j]])
    nodes.append([lons.max() + dx0, lats[j]])
    k = len(nodes)
    seg.append([k - 1, k])

nodes = np.array(nodes)
seg = np.array(seg)

transformer = Transformer.from_crs("epsg:4326", args.proj[0], always_xy=True)
nodes[:, 0], nodes[:, 1] = transformer.transform(nodes[:, 0], nodes[:, 1])

extra_column = np.zeros((nodes.shape[0], 1))
xyz = np.hstack((nodes, extra_column))
fn = args.output_prefix[0] + ".vtk"
with open(fn, "w") as fout:
    fout.write(
        f"""# vtk DataFile Version 2.0
parabola - polyline
ASCII
DATASET POLYDATA
POINTS {xyz.shape[0]} float
"""
    )
    np.savetxt(fout, xyz, "%e %e %e")
    fout.write(f"\nLINES {seg.shape[0]} {3*seg.shape[0]}\n")
    np.savetxt(fout, seg - 1, "2 %d %d")
    fout.write("\n")
print(f"done writing {fn}")
