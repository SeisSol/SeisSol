#!/usr/bin/env python3
import os
import argparse
from FaultPlane import FaultPlane, MultiFaultPlane
import os.path
from sklearn.decomposition import PCA
import numpy as np

parser = argparse.ArgumentParser(
    description="infer spatial zoom to match fault mesh size"
)
parser.add_argument("filename", help="filename of the srf file")
parser.add_argument("--mesh_size", help="fault mesh size", type=float)

parser.add_argument(
    "--proj",
    metavar="proj",
    nargs=1,
    help="proj4 string describing the projection",
    required=True,
)
args = parser.parse_args()

prefix, ext = os.path.splitext(args.filename)
prefix = os.path.basename(prefix)

if ext == ".srf":
    mfp = MultiFaultPlane.from_srf(args.filename)
elif ext == ".param":
    mfp = MultiFaultPlane.from_usgs_param_file(args.filename)
elif ext == ".fsp":
    mfp = MultiFaultPlane.from_usgs_fsp_file(args.filename)
elif ext == ".txt":
    mfp = MultiFaultPlane.from_slipnear_param_file(args.filename)
else:
    raise NotImplementedError(f" unknown extension: {ext}")

dx = float("inf")
dy = float("inf")
min_fault_plane_area = float("inf")
total_area = 0
for p, p1 in enumerate(mfp.fault_planes):
    p1.compute_xy_from_latlon(args.proj[0])
    xyz = np.column_stack((p1.x.flatten(), p1.y.flatten(), -p1.depth.flatten() * 1e3))
    # Perform PCA to get principal axes
    pca = PCA(n_components=2)
    points = pca.fit_transform(xyz) / 1e3
    la, lb = np.amax(points, axis=0) - np.amin(points, axis=0)
    min_fault_plane_area = min(min_fault_plane_area, la * lb)
    total_area += la * lb
    print("inferred fault dimensions (km)", la, lb)
    points = points.reshape((p1.ny, p1.nx, 2)) * 1e3
    iy = p1.ny // 2
    ix = p1.nx // 2
    dx = min(dx, abs(points[iy, ix, 0] - points[iy, ix - 1, 0]))
    dy = min(dy, abs(points[iy, ix, 1] - points[iy - 1, ix, 1]))
    print(dx, dy)


def next_odd_integer(x):
    if (x == int(x)) & (int(x) % 2 == 1):
        return int(x)
    nearest_integer = int(x)
    # If nearest integer is even, add 1, else add 2
    next_odd = nearest_integer + (1 if nearest_integer % 2 == 0 else 2)
    return next_odd


def get_fault_mesh_size(min_plane_area, total_area):
    if total_area < 40 * 80:
        return 500
    elif total_area > 100 * 200:
        return 1000
    else:
        return 700


if not os.path.exists("tmp"):
    os.makedirs("tmp")
if args.mesh_size:
    mesh_size = args.mesh_size[0]
else:
    mesh_size = get_fault_mesh_size(min_fault_plane_area, total_area)
print(f"using a mesh size of {mesh_size}")

inferred_spatial_zoom = next_odd_integer(min(dx, dy) / mesh_size)

print(f"inferred spatial zoom {inferred_spatial_zoom}")
fns = {
    "inferred_spatial_zoom.txt": inferred_spatial_zoom,
    "inferred_fault_mesh_size.txt": mesh_size,
}

for fn, value in fns.items():
    with open(f"tmp/{fn}", "w") as f:
        f.write(f"{value}\n")
    print(f"done writing tmp/{fn}")
