#!/usr/bin/env python3
import argparse
from kinmodmodules.fault_output_generator import generate

if __name__ == "__main__":
    # parsing python arguments
    parser = argparse.ArgumentParser(
        description="generate a fault output from FL33 input files"
    )
    parser.add_argument("fault_filename", help="fault.xdmf filename")
    parser.add_argument("yaml_filename", help="fault easi/yaml filename")

    parser.add_argument(
        "--output_file",
        help="path and prefix of the output file",
        nargs=1,
        default=["fault_from_fl33_input"],
    )
    parser.add_argument(
        "--stf",
        type=str,
        choices=["Yoffe", "Gaussian", "AsymmetricCosine"],
        default="Gaussian",
        help="the source time function to use",
    )
    parser.add_argument(
        "--dt",
        nargs=1,
        metavar="dt",
        default=[0.5],
        help="sampling time of the output file",
        type=float,
    )

    args = parser.parse_args()
    generate(
        args.fault_filename,
        args.yaml_filename,
        args.output_file[0],
        args.stf,
        args.dt[0],
    )
