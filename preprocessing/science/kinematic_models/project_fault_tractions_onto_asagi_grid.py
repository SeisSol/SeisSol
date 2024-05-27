#!/usr/bin/env python3
import argparse
from kinmodmodules.ugrid_data_projector import generate_input_files

if __name__ == "__main__":
    # parsing python arguments
    parser = argparse.ArgumentParser(
        description=(
            "project 3d fault output onto 2d grids to be read with Asagi. One grid per"
            " fault tag"
        )
    )
    parser.add_argument("fault_filename", help="fault.xdmf filename")
    parser.add_argument(
        "--dx",
        nargs=1,
        help="grid smapling",
        type=float,
        default=([100]),
    )
    parser.add_argument(
        "--gaussian_kernel",
        metavar="sigma_m",
        nargs=1,
        help="apply a gaussian kernel to smooth out intput stresses",
        type=float,
    )

    parser.add_argument(
        "--taper",
        nargs=1,
        help="tapper stress value (MPa)",
        type=float,
    )
    parser.add_argument(
        "--paraview_readable",
        dest="paraview_readable",
        action="store_true",
        help="write also netcdf readable by paraview",
        default=False,
    )

    args = parser.parse_args()
    generate_input_files(
        args.fault_filename,
        args.dx[0],
        args.gaussian_kernel,
        args.taper,
        args.paraview_readable,
    )
