#!/usr/bin/env python3
import gmsh
import numpy as np
import glob
import re
import argparse
import os

parser = argparse.ArgumentParser(
    description="generate yaml and netcdf input to be used with friction law 33/34 based on a (here"
    + "upsampled) kinematic model in the standard rupture format srf file."
)

parser.add_argument(
    "--domain_mesh_size",
    help="mesh size in the domain",
    nargs=1,
    type=float,
    default=[20000],
)
parser.add_argument(
    "--fault_mesh_size",
    help="mesh size on the faults",
    nargs=1,
    type=float,
    default=[1000],
)
parser.add_argument(
    "--interactive",
    dest="interactive",
    action="store_true",
    help="open the gui of gmsh once the mesh is generated",
)

args = parser.parse_args()


# mesh sizes
h_domain = args.domain_mesh_size[0]
h_fault = args.fault_mesh_size[0]

# domain dimensions
length_added = 120e3
z0, z1 = -length_added, 0

gmsh.initialize()
gmsh.model.add("finite-fault")

# We can log all messages for further processing with:
gmsh.logger.start()

# Regular expression patterns to match vertex lines
vertex_pattern = re.compile(r"VRTX (\d+) ([\d.-]+) ([\d.-]+) ([\d.-]+)")
allv = []
faults = []
ts_files = sorted(glob.glob(f"tmp/*.ts"))

for i, fn in enumerate(ts_files):
    vertices = []
    with open(fn, "r") as file:
        for line in file:
            match = vertex_pattern.match(line)
            if match:
                vertex_id, x, y, z = map(float, match.groups())
                if z > -h_fault:
                    z = 0.0
                vertices.append([x, y, z])
    allv.extend(vertices)
    vertices = np.array(vertices)
    # Define the points
    point1 = gmsh.model.occ.addPoint(*vertices[0, :])
    point2 = gmsh.model.occ.addPoint(*vertices[1, :])
    # Define the B-spline curve
    curve = gmsh.model.occ.addBSpline([point1, point2])
    # Extrude the curve
    extrude_dir = vertices[2, :] - vertices[1, :]
    extrusion_results = gmsh.model.occ.extrude([(1, curve)], *extrude_dir)
    # Extract the tuple with 2 as the first item
    fault = next((t for t in extrusion_results if t[0] == 2), None)
    faults.append(fault)

# compute fault normals
gmsh.model.occ.synchronize()
fault_normals = []
for i, fn in enumerate(ts_files):
    fault_normal = gmsh.model.getNormal(fault[1], [0.5, 0.5])
    fault_normals.append(fault_normal)


allv = np.array(allv)
min_x, max_x = np.min(allv[:, 0]), np.max(allv[:, 0])
min_y, max_y = np.min(allv[:, 1]), np.max(allv[:, 1])

x0 = min_x - length_added
x1 = max_x + length_added
y0 = min_y - length_added
y1 = max_y + length_added
box = gmsh.model.occ.addBox(x0, y0, z0, x1 - x0, y1 - y0, z1 - z0)

gmsh.model.occ.synchronize()

ov, ovv = gmsh.model.occ.fragment([(3, box)], faults)
gmsh.model.occ.synchronize()

# now retrieves which surface is where

all_points = gmsh.model.getEntities(dim=0)
coords = np.zeros((len(all_points), 3))
pt2vertexid = {}
for i, pt in enumerate(all_points):
    coords[i, :] = gmsh.model.getValue(*pt, [])
    pt2vertexid[pt[1]] = i
print(coords)

all_surfaces = gmsh.model.getEntities(dim=2)

points_of_surface = {}
tags = {}
for i in [1, 5]:
    tags[i] = []
for i, fn in enumerate(ts_files):
    fault_id = 3 if i == 0 else 64 + i
    tags[fault_id] = []


for surface in all_surfaces:
    curves = gmsh.model.getBoundary([surface])
    pts = set()
    for cu in curves:
        points = gmsh.model.getBoundary([cu])
        points = [pt[1] for pt in points]
        pts.update(points)
    points_of_surface[surface] = pts
    vids = [pt2vertexid[pt] for pt in pts]
    surf_coords = coords[list(vids), :]
    zmin = np.min(surf_coords[:, 2])
    stag = surface[1]
    # print('s',stag, pts)
    if zmin > -0.01:
        # free surface
        tags[1].append(stag)
    elif abs(zmin - z0) < 0.01:
        # absorbing
        tags[5].append(stag)
    else:
        tagged = False
        for i, fn in enumerate(ts_files):
            fault_id = 3 if i == 0 else 64 + i
            j = i * 4
            nb_match = 0
            for k in range(4):
                v1 = allv[j + k, :]
                dist = np.min(np.linalg.norm(surf_coords - v1, axis=1))
                # print(i,k, dist)
                if dist < 200:
                    nb_match += 1
            if nb_match > 1:
                tags[fault_id].append(stag)
                tagged = True
                break
        if not tagged:
            raise ValueError(
                f"surface {stag} could not be tagged, {nb_match}, {surf_coords}"
            )
print(tags)


for key in tags.keys():
    h = h_domain if key in [1, 5] else h_fault
    pairs = [(2, tag) for tag in tags[key]]
    gmsh.model.mesh.setSize(gmsh.model.getBoundary(pairs, False, False, True), h)
    gmsh.model.addPhysicalGroup(2, tags[key], key)


fault_faces = []
for key in tags.keys():
    if key not in [1, 5]:
        fault_faces.extend(tags[key])
print(fault_faces)

# Set mesh size based on a distance from faults field
gmsh.model.mesh.field.add("Distance", 1)
gmsh.model.mesh.field.setNumbers(1, "SurfacesList", fault_faces)
gmsh.model.mesh.field.setNumber(1, "Sampling", 100)
gmsh.model.mesh.field.add("MathEval", 2)
gmsh.model.mesh.field.setString(2, "F", f"20*F1^(0.5) + 0.1*F1 + {h_fault}")
gmsh.model.mesh.field.setAsBackgroundMesh(2)
gmsh.model.geo.synchronize()

gmsh.model.addPhysicalGroup(3, [1], 1)
gmsh.option.setNumber("Mesh.Algorithm3D", 10)
gmsh.model.mesh.generate(3)

if not os.path.exists("tmp"):
    os.makedirs("tmp")
gmsh.write("tmp/mesh.msh")

if args.interactive:
    gmsh.fltk.run()
gmsh.finalize()
