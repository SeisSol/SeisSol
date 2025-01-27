# SPDX-FileCopyrightText: 2021-2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

import argparse
import os
import sys

import seissol_proxy_bindings as pb


def main():
    currentdir = os.path.dirname(os.path.abspath(__file__))
    parentdir = os.path.dirname(currentdir)
    sys.path.insert(0, parentdir)

    parser = argparse.ArgumentParser()
    kernels_options = pb.Aux.get_allowed_kernels()
    parser.add_argument(
        "-c",
        "--cells",
        default=100000,
        type=int,
        help="num cells in a time cluster",
    )
    parser.add_argument(
        "-t",
        "--timesteps",
        default=20,
        type=int,
        help="num time steps/repeats",
    )
    parser.add_argument(
        "-k",
        "--kernel",
        default="all",
        choices=kernels_options,
        type=str,
        help="kernel types",
    )
    args = parser.parse_args()

    config = pb.ProxyConfig()
    config.cells = args.cells
    config.timesteps = args.timesteps
    config.kernel = pb.Aux.str_to_kernel(args.kernel)

    output = pb.run_proxy(config)
    pb.Aux.display_output(output, args.kernel)


if __name__ == "__main__":
    main()
