#!/usr/bin/env python3

import argparse
from kinmodmodules.ugrid_data_projector import generate_input_files


def main() -> None:
    """
    Main function to parse arguments and generate input files
    by projecting 3D fault output onto 2D grids for ASAGI.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Project 3D fault output onto 2D grids to be read with ASAGI. "
            "One grid per fault tag."
        )
    )
    parser.add_argument("fault_filename", help="Fault.xdmf filename")
    parser.add_argument(
        "--dx",
        nargs=1,
        help="Grid sampling",
        type=float,
        default=[100.0],
    )
    parser.add_argument(
        "--gaussian_kernel",
        metavar="sigma_m",
        nargs=1,
        help="Apply a Gaussian kernel to smooth out input stresses",
        type=float,
    )
    parser.add_argument(
        "--taper",
        nargs=1,
        help="Taper stress value (MPa)",
        type=float,
    )
    parser.add_argument(
        "--paraview_readable",
        dest="paraview_readable",
        action="store_true",
        help="Write NetCDF files readable by ParaView",
        default=False,
    )

    args = parser.parse_args()
    generate_input_files(
        args.fault_filename,
        args.dx[0],
        args.gaussian_kernel[0] if args.gaussian_kernel else None,
        args.taper[0] if args.taper else None,
        args.paraview_readable,
    )


if __name__ == "__main__":
    main()
