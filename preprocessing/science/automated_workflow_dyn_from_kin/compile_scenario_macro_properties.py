#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import argparse
import matplotlib
import os
import glob
import re
from obspy.signal.cross_correlation import correlate, correlate_template, xcorr_max
from scipy import integrate


def infer_duration(time, moment_rate):
    moment = integrate.cumulative_trapezoid(moment_rate, time, initial=0)
    M0 = np.trapz(moment_rate[:], x=time[:])
    return np.amax(time[moment < 0.99 * M0])


def computeMw(label, time, moment_rate):
    M0 = np.trapz(moment_rate[:], x=time[:])
    Mw = 2.0 * np.log10(M0) / 3.0 - 6.07
    # print(f"{label} moment magnitude: {Mw:.2} (M0 = {M0:.4e})")
    return M0, Mw


def extractBCR(fname):
  # Define a regular expression pattern to extract B, C, and R values
  pattern = r'B([\d.]+)_C([\d.]+)_R([\d._]+)-energy.csv'

  # Match the pattern against the filename
  match = re.search(pattern, fname)

  if match:
      # Extract B and C values
      B_value = float(match.group(1))
      C_value = float(match.group(2))

      # Extract R values and convert to a list
      R_value = list(map(float, match.group(3).split('_')))
      return B_value, C_value, R_value
  else:
        # Raise a custom exception if no match is found
        raise ValueError(f"No match found in the file name: {fname}")


def generate_XY_panel(
    name_arr1, arr1, name_arr2, arr2, name3, val3, name_col, ax, cmap
):
    "generate a 2D plot with name_arr1 and name_arr2 in X and Y axis"
    "colored by name_col"
    unique_arr1 = np.unique(arr1)
    unique_arr2 = np.unique(arr2)
    X, Y = np.meshgrid(unique_arr1, unique_arr2)
    gof_array = np.zeros_like(X) + np.nan
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            result = result_df[
                (result_df[name_arr1] == X[i, j])
                & (result_df[name_arr2] == Y[i, j])
                & (result_df[name3] == val3)
            ]
            if not result.empty:
                specific_index = result.index[0]
                gof_array[i, j] = result.loc[specific_index, name_col]
    im = ax.pcolormesh(X, Y, gof_array, cmap=cmap)
    im.set_clim(0, result_df[name_col].max())
    ax.set_xlabel(name_arr1)
    ax.set_ylabel(name_arr2)
    ax.set_xticks(unique_arr1)
    ax.set_yticks(unique_arr2)
    ax.set_title(f"{name3}={val3}")
    if name_col == "ccmax":
        label = "gof usgs moment rate release"
    elif name_col == "M0mis":
        label = "gof usgs moment release"
    else:
        label = name_col

    fig.colorbar(im, label=label, ax=ax)


def generate_BCR_plots(B, C, R):
    unique_R = np.unique(R)
    n_div = 2
    nrow, ncol = len(unique_R) // n_div, 2 * n_div
    # nrow, ncol = 2, 2
    fig, axarr = plt.subplots(
        nrow,
        ncol,
        figsize=(ncol * 4, nrow * 4),
        dpi=160,
        sharex=False,
        sharey=True,
        squeeze=False,
    )

    for k, Rk in enumerate(unique_R):
        row = k % nrow
        col = k // nrow * n_div
        generate_XY_panel("B", B, "C", C, "R0", Rk, "ccmax", axarr[row, col], "viridis")
        generate_XY_panel(
            "B", B, "C", C, "R0", Rk, "M0mis", axarr[row, col + 1], "plasma"
        )
    fname = "BC_constant_R.pdf"
    plt.savefig(fname)
    print(f"done writing {fname}")

    unique_B = np.unique(B)
    n_div = 1
    nrow, ncol = len(unique_B) // n_div, 2 * n_div
    fig, axarr = plt.subplots(
        nrow,
        ncol,
        figsize=(ncol * 4, nrow * 4),
        dpi=160,
        sharex=False,
        sharey=True,
        squeeze=False,
    )

    for k, Bk in enumerate(unique_B):
        row = k % nrow
        col = k // nrow * n_div
        generate_XY_panel("R0", R, "C", C, "B", Bk, "ccmax", axarr[row, col], "viridis")
        generate_XY_panel(
            "R0", R, "C", C, "B", Bk, "M0mis", axarr[row, col + 1], "plasma"
        )
    fname = "R0C_constant_B.pdf"
    plt.savefig(fname)
    print(f"done writing {fname}")


if __name__ == "__main__":
    ps = 8
    matplotlib.rcParams.update({"font.size": ps})
    plt.rcParams["font.family"] = "sans"
    matplotlib.rc("xtick", labelsize=ps)
    matplotlib.rc("ytick", labelsize=ps)
    matplotlib.rcParams["lines.linewidth"] = 1.0

    parser = argparse.ArgumentParser(description="compute scenario properties")
    parser.add_argument(
        "output_folder",
        help="path to output folder or full path to specific energy.csv file",
    )
    parser.add_argument(
        "--extension",
        help="figure extension (without the .)",
        nargs=1,
        default=["pdf"],
    )
    parser.add_argument(
        "--gof_threshold",
        help="gof threshold from which results are selected",
        nargs=1,
        type=float,
        default=[0.6],
    )
    parser.add_argument(
        "--nmin",
        help="minimum number of synthetic moment rates drawn in color",
        nargs=1,
        type=int,
        default=[3],
    )
    parser.add_argument(
        "--nmax",
        help="maximum number of synthetic moment rates drawn in color",
        nargs=1,
        type=int,
        default=[10],
    )

    args = parser.parse_args()

    if not os.path.exists("plots"):
        os.makedirs("plots")

    if args.output_folder.endswith("energy.csv"):
        energy_files = [args.output_folder]
    else:
        if os.path.exists(args.output_folder):
            args.output_folder += "/"
        energy_files = sorted(glob.glob(f"{args.output_folder}*-energy.csv"))
    one_model_shown = args.nmax[0] == 1

    cm = 1 / 2.54
    figsize = (6.5 * cm, 3.5 * cm) if one_model_shown else (8, 4)
    fig = plt.figure(figsize=figsize, dpi=80)
    ax = fig.add_subplot(111)

    results = {
        "B": [],
        "C": [],
        "R0": [],
        "Mw": [],
        "ccmax": [],
        "M0mis": [],
        "faultfn": [],
        "shift_syn_ref_sec": [],
    }

    def read_usgs_moment_rate():
        mr_ref = np.loadtxt("tmp/moment_rate.mr", skiprows=2)
        ref_name = "usgs"
        # Conversion factor from dyne-cm/sec to Nm/sec (for older usgs files)
        scaling_factor = 1.0 if np.amax(mr_ref[:, 1]) < 1e23 else 1e-7
        mr_ref[:, 1] *= scaling_factor
        return mr_ref

    def trim_trailing_zero(mr_ref):
        last_index_non_zero = np.nonzero(mr_ref[:, 1])[0][-1]
        return mr_ref[:last_index_non_zero, :]

    if os.path.exists("tmp/moment_rate_slipnear.txt"):
        print("loading slipnear moment rate")
        mr_ref = np.loadtxt("tmp/moment_rate_slipnear.txt")
        ref_name = "slipnear"
    else:
        mr_ref = read_usgs_moment_rate()
        ref_name = "usgs"

    mr_ref = trim_trailing_zero(mr_ref)

    M0ref, Mwref = computeMw(ref_name, mr_ref[:, 0], mr_ref[:, 1])

    inferred_duration = infer_duration(mr_ref[:, 0], mr_ref[:, 1])

    # Create a new time array with the desired time step
    dt = 0.25
    new_time = np.arange(mr_ref[:, 0].min(), mr_ref[:, 0].max() + dt, dt)
    # Use numpy's interp function to interpolate values at the new time points
    mr_ref_interp = np.interp(new_time, mr_ref[:, 0], mr_ref[:, 1])

    if 2 * inferred_duration > new_time[-1]:
        added = int((2 * inferred_duration - new_time.max()) / dt)
        mr_ref_interp = np.pad(
            mr_ref_interp, (0, added), "constant", constant_values=(0, 0)
        )

    for i, fn in enumerate(energy_files):
        if "fl33" in fn:
            continue
        df = pd.read_csv(fn)
        df = df.pivot_table(index="time", columns="variable", values="measurement")
        if len(df) < 2:
            print(f"skipping empty {fn}")
            continue
        dt = df.index[1] - df.index[0]
        assert dt == 0.25
        df["seismic_moment_rate"] = np.gradient(df["seismic_moment"], dt)
        label = os.path.basename(fn)
        prefix = fn.split("-energy.csv")[0]
        faultfn = glob.glob(f"{prefix}*-fault.xdmf")[0]
        B, C, R = extractBCR(fn)
        M0, Mw = computeMw(label, df.index.values, df["seismic_moment_rate"])
        results["Mw"].append(Mw)
        results["B"].append(B)
        results["C"].append(C)
        results["R0"].append(R)
        results["faultfn"].append(faultfn)
        if len(mr_ref_interp) > len(df["seismic_moment_rate"]):
            cc = correlate_template(
                mr_ref_interp,
                df["seismic_moment_rate"],
                mode="same",
                normalize="naive",
            )
        else:
            cc = correlate_template(
                df["seismic_moment_rate"],
                mr_ref_interp,
                mode="same",
                normalize="naive",
            )
        # cc = correlate(df["seismic_moment_rate"], mr_ref_interp, shift=int(0.25 * inferred_duration/dt))
        shift, ccmax = xcorr_max(cc, abs_max=False)
        d1 = max(len(df["seismic_moment_rate"]), len(mr_ref_interp))
        d2 = min(len(df["seismic_moment_rate"]), len(mr_ref_interp))
        results["shift_syn_ref_sec"].append((shift + 0.5 * (d1 - d2)) * dt)
        if abs(results["shift_syn_ref_sec"][-1]) > 0.5 * inferred_duration:
            ccmax = 0.0
        results["ccmax"].append(ccmax)
        M0_gof = 1 - abs(M0 - M0ref) / M0ref
        results["M0mis"].append(M0_gof)

    result_df = pd.DataFrame(results)
    result_df["overall_gof"] = np.sqrt(result_df["M0mis"] * result_df["ccmax"])
    print(result_df)
    B = result_df["B"].values
    C = result_df["C"].values
    R = result_df["R0"].values

    generate_BCR_plots(B, C, R)

    overall_gof = result_df["overall_gof"].values
    Mw = result_df["Mw"].values
    indices_of_nlargest_values = result_df["overall_gof"].nlargest(args.nmin[0]).index
    indices_of_nmax_largest_values = (
        result_df["overall_gof"].nlargest(args.nmax[0]).index
    )
    indices_greater_than_threshold = result_df[
        result_df["overall_gof"] > args.gof_threshold[0]
    ].index
    if len(indices_greater_than_threshold) > args.nmax[0]:
        selected_indices = indices_of_nmax_largest_values
    else:
        selected_indices = indices_greater_than_threshold
    i = 0
    for fn in energy_files:
        prefix_to_match = fn.split("-energy.csv")[0]
        row_with_prefix1 = result_df[
            result_df["faultfn"].str.startswith(prefix_to_match + "_")
        ]
        row_with_prefix2 = result_df[
            result_df["faultfn"].str.startswith(prefix_to_match + "-")
        ]
        if row_with_prefix1.empty and row_with_prefix2.empty:
            continue
        df = pd.read_csv(fn)
        df = df.pivot_table(index="time", columns="variable", values="measurement")
        dt = df.index[1] - df.index[0]
        assert dt == 0.25
        df["seismic_moment_rate"] = np.gradient(df["seismic_moment"], dt)

        if one_model_shown:
            label = f"simulation"
        else:
            label = f"B={B[i]}, C={C[i]}, R={R[i]}"
        if i in selected_indices or i in indices_of_nlargest_values:
            if one_model_shown:
                labelargs = {"label": f"{label} (Mw={Mw[i]:.2f})"}
            else:
                labelargs = {
                    "label": f"{label} (Mw={Mw[i]:.2f}, gof={overall_gof[i]:.2})"
                }
            alpha = 1.0
        else:
            labelargs = {"color": "lightgrey", "zorder": 1}
            alpha = 0.5

        ax.plot(
            df.index.values,
            df["seismic_moment_rate"] / 1e19,
            alpha=alpha,
            **labelargs,
        )
        i += 1

    ax.plot(
        mr_ref[:, 0],
        mr_ref[:, 1] / 1e19,
        label=f"{ref_name} (Mw={Mwref:.2f})",
        color="black",
    )
    if ref_name != "usgs":
        mr_usgs = read_usgs_moment_rate()
        mr_usgs = trim_trailing_zero(mr_usgs)
        M0usgs, Mwusgs = computeMw("usgs", mr_usgs[:, 0], mr_usgs[:, 1])
        ax.plot(
            mr_usgs[:, 0],
            mr_usgs[:, 1] / 1e19,
            label=f"usgs (Mw={Mwusgs:.2f})",
            color="darkblue",
        )

    selected_rows = result_df[result_df["Mw"] > 6]
    selected_rows = selected_rows.sort_values(
        by="overall_gof", ascending=False
    ).reset_index(drop=True)
    print(selected_rows)

    fname = "tmp/selected_output.txt"
    with open(fname, "w") as fid:
        for index, row in selected_rows.iterrows():
            fid.write(f"{row['faultfn']}\n")
            if index == 3:
                break
        fid.write("output/dyn-usgs-fault.xdmf\n")
    print(f"done writing {fname}")

    col = 1 if one_model_shown else 2
    ax.legend(frameon=False, loc="upper right", ncol=col, fontsize=8)
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    ax.set_ylabel(r"moment rate (e19 $\times$ Nm/s)")
    ax.set_xlabel("time (s)")

    fn = f"plots/moment_rate.{args.extension[0]}"
    fig.savefig(fn, bbox_inches="tight", transparent=True)
    print(f"done write {fn}")
    full_path = os.path.abspath(fn)
    print(f"full path: {full_path}")
