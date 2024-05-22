#!/usr/bin/env python3
import argparse
from kinmodmodules.moment_rate_calculator import compute

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "compute moment rate release from kinematic model. Assumptions: Gaussian source"
            " time function, no bimaterial conditions"
        )
    )
    parser.add_argument(
        "--dt",
        nargs=1,
        metavar="dt",
        default=[0.25],
        help="sampling time of the output file",
        type=float,
    )

    parser.add_argument("filename", help="filename of the srf file")
    parser.add_argument("yaml_filename", help="fault easi/yaml filename")


    parser.add_argument(
        "--proj",
        metavar="proj",
        nargs=1,
        help="proj4 string describing the projection",
        required=True,
    )
    args = parser.parse_args()
    compute(args.filename, args.yaml_file, args.dt[0], args.proj[0])
