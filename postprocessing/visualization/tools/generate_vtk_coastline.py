#!/usr/bin/env python3
import argparse
import os
import numpy as np
import seissolxdmf as sx
from pyproj import Transformer

def parse_args():
    parser = argparse.ArgumentParser(description="Create coastline VTK file from GMT")
    parser.add_argument(
        "--lon",
        nargs=2,
        metavar=("lonmin", "lonmax"),
        help="Minimum and maximum longitude",
        type=float,
    )
    parser.add_argument(
        "--lat",
        nargs=2,
        metavar=("latmin", "latmax"),
        help="Minimum and maximum latitude",
        type=float,
    )
    parser.add_argument(
        "--from_fault_output",
        nargs=1,
        metavar=("filename added_degree"),
        help="Infer domain from range of fault xdmf + added extension in degree",
    )
    parser.add_argument(
        "--proj",
        nargs=1,
        metavar=("projname"),
        help="Project the data. projname: proj4 string describing the projection",
    )
    parser.add_argument(
        "--resolution",
        nargs=1,
        metavar=("resolution"),
        default=("i"),
        help="Coastline resolution, (f)ull, (h)igh, (i)ntermediate, (l)ow, and (c)rude",
    )
    parser.add_argument(
        "--translate",
        nargs=2,
        metavar=("x", "y"),
        help="Translate coordinate array (e.g. x_new = x_old - x)",
        type=float,
    )
    parser.add_argument(
        "--z",
        nargs=1,
        metavar=("elevation"),
        default=([0]),
        help="Z coordinate of coastline",
        type=float,
    )
    return parser.parse_args()

def infer_domain_from_fault_output(args):
    fn, added_deg = args.from_fault_output[0].split()
    s = sx.seissolxdmf(fn)
    added_deg = float(added_deg)
    xyz = s.ReadGeometry()
    if args.proj:
        transformer = Transformer.from_crs(args.proj[0], "epsg:4326", always_xy=True)
        xyz[:, 0], xyz[:, 1], xyz[:, 2] = transformer.transform(xyz[:, 0], xyz[:, 1], xyz[:, 2])
        args.lon = [0, 0]
        args.lat = [0, 0]
        args.lon[1], args.lat[1], _ = np.amax(xyz, axis=0) + added_deg
        args.lon[0], args.lat[0], _ = np.amin(xyz, axis=0) - added_deg
        print(f"Inferred domain dims, lon: {args.lon}, lat: {args.lat}")

def run_gmt_pscoast(args):
    os.system(
        f"gmt pscoast -R{args.lon[0]}/{args.lon[1]}/{args.lat[0]}/{args.lat[1]} -D{args.resolution[0]} -M -W > coastline.dat"
    )

def read_gmt_file():
    xyz = []
    segments = []
    nvert = 0
    new_polyline = True
    with open("coastline.dat") as fid:
        for line in fid:
            if line.startswith("#"):
                continue
            if line.startswith(">"):
                new_polyline = True
            else:
                xyz.append([float(val) for val in line.split()])
                nvert = nvert + 1
                if not new_polyline:
                    segments.append([nvert - 1, nvert])
                new_polyline = False
    return xyz, segments

def project_coordinates(args, xyz):
    transformer = Transformer.from_crs("epsg:4326", args.proj[0], always_xy=True)
    xyz[:, 0], xyz[:, 1], xyz[:, 2] = transformer.transform(xyz[:, 0], xyz[:, 1], xyz[:, 2])

def translate_coordinates(args, xyz):
    xyz[:, 0] -= args.translate[0]
    xyz[:, 1] -= args.translate[1]

def write_vtk_file(xyz, segments):
    nvert = xyz.shape[0]
    nseg = segments.shape[0]
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

def main():
    args = parse_args()
    if args.from_fault_output:
        infer_domain_from_fault_output(args)
    run_gmt_pscoast(args)
    xyz, segments = read_gmt_file()
    xyz = np.asarray(xyz)
    xyz = np.insert(xyz, xyz.shape[1], 0, axis=1)
    segments = np.asarray(segments)
    if args.proj:
        project_coordinates(args, xyz)
    xyz[:, 2] = args.z[0]
    if args.translate:
        translate_coordinates(args, xyz)
    write_vtk_file(xyz, segments)

if __name__ == "__main__":
    main()

