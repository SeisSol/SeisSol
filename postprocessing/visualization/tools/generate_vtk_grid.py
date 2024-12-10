#!/usr/bin/env python3
import pyproj
import numpy as np
import argparse
import seissolxdmf as sx
from pyproj import Transformer


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate lat/lon grid for Paraview visualization"
    )
    parser.add_argument(
        "--dx", nargs=1, help="Sampling between grid lines", type=float, required=True
    )
    parser.add_argument(
        "--lon",
        nargs=2,
        metavar=("lonmin", "lonmax"),
        help="Minimum (lonmin) and maximum (lonmax) longitude",
        type=float,
    )
    parser.add_argument(
        "--lat",
        nargs=2,
        metavar=("latmin", "latmax"),
        help="Minimum (latmin) and maximum (latmax) latitude",
        type=float,
    )
    parser.add_argument(
        "--from_fault_output",
        nargs=1,
        metavar=("filename added_degree"),
        help="Infer domain from range of fault xdmf + added extension in degree",
    )
    parser.add_argument(
        "--output_prefix", nargs=1, help="Prefix for output VTK file", default=["grid"]
    )
    parser.add_argument(
        "--lineout",
        nargs=2,
        metavar=("lon2", "lat2"),
        default=("0.35", "0.35"),
        help="Draw line out from the grid of length lon2 (or lat2)",
    )
    parser.add_argument(
        "--proj",
        nargs=1,
        metavar=("projname"),
        help="Project the data. projname: proj4 string describing the projection",
    )
    return parser.parse_args()


def my_round(x: float, dx: float, floor: bool = True) -> float:
    func = np.floor if floor else np.ceil
    return func(x / dx) * dx


def generate_grid(args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray]:
    if args.from_fault_output:
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
        args.lon[1], args.lat[1] = [
            my_round(x, args.dx[0], False) for x in [args.lon[1], args.lat[1]]
        ]
        args.lon[0], args.lat[0] = [
            my_round(x, args.dx[0], True) for x in [args.lon[0], args.lat[0]]
        ]
    lons = np.arange(args.lon[0], args.lon[1], args.dx[0])
    lats = np.arange(args.lat[0], args.lat[1], args.dx[0])
    return lons, lats


def create_nodes_and_segments(
    lons: np.ndarray, lats: np.ndarray, args: argparse.Namespace
) -> tuple[np.ndarray, np.ndarray]:
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
    return nodes, seg


def project_nodes(args: argparse.Namespace, nodes: np.ndarray) -> np.ndarray:
    transformer = Transformer.from_crs("epsg:4326", args.proj[0], always_xy=True)
    nodes[:, 0], nodes[:, 1] = transformer.transform(nodes[:, 0], nodes[:, 1])
    return nodes


def write_vtk_file(
    nodes: np.ndarray, seg: np.ndarray, args: argparse.Namespace
) -> None:
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


def main() -> None:
    args = parse_args()
    lons, lats = generate_grid(args)
    nodes, seg = create_nodes_and_segments(lons, lats, args)
    nodes = project_nodes(args, nodes)
    write_vtk_file(nodes, seg, args)


if __name__ == "__main__":
    main()
