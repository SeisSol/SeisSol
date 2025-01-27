# SPDX-FileCopyrightText: 2021-2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

import argparse
import os
import sys
import time

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seissol_proxy_bindings as pb
from progress.bar import Bar


def main():
    currentdir = os.path.dirname(os.path.abspath(__file__))
    parentdir = os.path.dirname(currentdir)
    sys.path.insert(0, parentdir)
    parser = argparse.ArgumentParser()
    plotting_options = ["csv", "json", "matplotlib"]
    parser.add_argument(
        "-c",
        "--cells",
        default=50000,
        type=int,
        help="num cells in a time cluster",
    )
    parser.add_argument(
        "-t",
        "--timesteps",
        default=30,
        type=int,
        help="num time steps/repeats",
    )
    parser.add_argument(
        "--output_type",
        choices=plotting_options,
        default="csv",
        help="format to save data",
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        default=".",
        type=str,
        help="relative output directory",
    )
    args = parser.parse_args()

    config = pb.ProxyConfig()
    config.cells = args.cells
    config.timesteps = args.timesteps
    config.verbose = False

    df = pd.DataFrame(columns=["kernel type", "time", "HW GFLOPS", "NZ GFLOPS"])
    kernels = pb.Aux.get_allowed_kernels()
    with Bar("proxy...    ", max=len(kernels)) as bar:
        for index, kernel in enumerate(kernels):
            config.kernel = pb.Aux.str_to_kernel(kernel)
            result = pb.run_proxy(config)
            df.loc[index] = [
                kernel,
                result.time,
                result.hardware_gflops,
                result.non_zero_gflops,
            ]
            bar.next()

    # prepare a unique file suffix
    time_obj = time.localtime()
    now = time.strftime("%B-%H-%M-%S_%Y", time_obj)

    # prepare a directory where to write data
    output_dir = os.path.join(os.getcwd(), args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # write data to files
    if args.output_type == "json":
        with open(f"{output_dir}/proxy_{now}.json", "w") as file:
            file.write(df.to_json())
    elif args.output_type == "matplotlib":
        mpl.style.use("ggplot")
        df.plot(kind="bar", x="kernel type", y="HW GFLOPS", rot=45)
        plt.savefig(
            f"{output_dir}/proxy_hw-gflops_{now}.png",
            dpi=300,
            bbox_inches="tight",
        )

        df.plot(kind="bar", x="kernel type", y="NZ GFLOPS", rot=45)
        plt.savefig(
            f"{output_dir}/proxy_nz-gflops_{now}.png",
            dpi=300,
            bbox_inches="tight",
        )

        df.plot(kind="bar", x="kernel type", y="time", rot=45)
        plt.savefig(f"{output_dir}/proxy_time_{now}.png", dpi=300, bbox_inches="tight")
    elif args.output_type == "csv":
        with open(f"{output_dir}/proxy_{now}.pd", "w") as file:
            file.write(df.to_csv())

    # print out to console
    print(df)


if __name__ == "__main__":
    main()
