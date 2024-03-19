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
    # Define a regular expression pattern to extract numeric values
    pattern = re.compile(r"B([\d.]+)_C([\d.]+)_R([\d.]+)-energy.csv")

    # Use the pattern to find matches in the file name
    matches = pattern.search(fname)

    if matches:
        # Extract the values of B, C, and R from the matched groups
        B_value = float(matches.group(1))
        C_value = float(matches.group(2))
        R_value = float(matches.group(3))
        return B_value, C_value, R_value
    else:
        # Raise a custom exception if no match is found
        raise ValueError(f"No match found in the file name: {fname}")


if __name__ == "__main__":
    ps = 8
    matplotlib.rcParams.update({"font.size": ps})
    plt.rcParams["font.family"] = "sans"
    matplotlib.rc("xtick", labelsize=ps)
    matplotlib.rc("ytick", labelsize=ps)
    matplotlib.rcParams["lines.linewidth"] = 1.0

    parser = argparse.ArgumentParser(description="compute scenario properties")
    parser.add_argument("output_folder", help="path to output folder")
    args = parser.parse_args()

    if not os.path.exists("plots"):
        os.makedirs("plots")

    fig = plt.figure(figsize=(8, 4), dpi=80)
    ax = fig.add_subplot(111)

    energy_files = sorted(glob.glob(f"{args.output_folder}/*-energy.csv"))
    results = {
        "B": [],
        "C": [],
        "R0": [],
        "Mw": [],
        "ccmax": [],
        "M0mis": [],
        "faultfn": [],
    }

    mr_usgs = np.loadtxt("tmp/moment_rate.mr", skiprows=2)
    last_index_non_zero = np.nonzero(mr_usgs[:, 1])[0][-1]
    mr_usgs = mr_usgs[:last_index_non_zero, :]
    # Conversion factor from dyne-cm/sec to Nm/sec (for older usgs files)
    scaling_factor = 1.0 if np.amax(mr_usgs[:, 1]) < 1e23 else 1e-7
    mr_usgs[:, 1] *= scaling_factor

    M0usgs, Mwusgs = computeMw("usgs", mr_usgs[:, 0], mr_usgs[:, 1])

    usgs_duration = infer_duration(mr_usgs[:, 0], mr_usgs[:, 1])
    print(usgs_duration)
    # Create a new time array with the desired time step
    dt = 0.25
    new_time = np.arange(mr_usgs[:, 0].min(), mr_usgs[:, 0].max() + dt, dt)
    # Use numpy's interp function to interpolate values at the new time points
    mr_usgs_interp = np.interp(new_time, mr_usgs[:, 0], mr_usgs[:, 1])

    for i, fn in enumerate(energy_files):
        if "fl33" in fn:
            continue
        df = pd.read_csv(fn)
        df = df.pivot_table(index="time", columns="variable", values="measurement")
        dt = df.index[1] - df.index[0]
        assert dt == 0.25
        df["seismic_moment_rate"] = np.gradient(df["seismic_moment"], dt)
        label = os.path.basename(fn)
        faultfn = fn.split("-energy.csv")[0] + "-fault.xdmf"
        B, C, R = extractBCR(fn)
        M0, Mw = computeMw(label, df.index.values, df["seismic_moment_rate"])
        results["Mw"].append(Mw)
        results["B"].append(B)
        results["C"].append(C)
        results["R0"].append(R)
        results["faultfn"].append(faultfn)
        if len(mr_usgs_interp) > len(df["seismic_moment_rate"]):
            cc = correlate_template(
                mr_usgs_interp,
                df["seismic_moment_rate"],
                mode="same",
                normalize="naive",
            )
        else:
            cc = correlate_template(
                df["seismic_moment_rate"],
                mr_usgs_interp,
                mode="same",
                normalize="naive",
            )
        # cc = correlate(stcc[i], stcc[j], shift=100)
        shift, ccmax = xcorr_max(cc, abs_max=False)
        if abs(shift * dt) > 0.75 * usgs_duration:
            ccmax = 0.0
        results["ccmax"].append(ccmax)
        M0_gof = 1 - abs(M0 - M0usgs) / M0usgs
        results["M0mis"].append(M0_gof)

        label = f"B={B}, C={C}, R={R}"
        if Mw < 6.0:
            continue
        overall_gof = M0_gof * ccmax
        if overall_gof > 0.6:
            labelargs = {"label": f"{label} (Mw={Mw:.2f}, gof={overall_gof:.2})"}
            alpha = 1.0
        else:
            labelargs = {"color": "lightgrey", "zorder": 1}
            alpha = 0.5

        ax.plot(
            df.index.values,
            df["seismic_moment_rate"] / 1e19,
            alpha=alpha,
            linestyle="--" if C == 0.3 else "-",
            **labelargs,
        )
    ax.plot(
        mr_usgs[:, 0],
        mr_usgs[:, 1] / 1e19,
        label=f"usgs (Mw={Mwusgs:.2f})",
        color="black",
    )
    result_df = pd.DataFrame(results)
    result_df["overall_gof"] = result_df["M0mis"] * result_df["ccmax"]
    print(result_df)

    selected_rows = result_df[result_df["Mw"] > 6]
    selected_rows = selected_rows.sort_values(by="overall_gof", ascending=False).reset_index(drop=True)
    print(selected_rows)

    fname = "tmp/selected_output.txt"
    with open(fname, "w") as fid:
        for index, row in selected_rows.iterrows():
            fid.write(f"{row['faultfn']}\n")
            if index == 3:
                break
        fid.write("output/dyn-usgs-fault.xdmf\n")
    print(f"done writing {fname}")

    ax.legend(frameon=False, loc="upper right", ncol=2, fontsize=6)
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    ax.set_ylabel(r"moment rate (e19 $\times$ Nm/s)")
    ax.set_xlabel("time (s)")

    fn = f"plots/moment_rate.pdf"
    fig.savefig(fn, bbox_inches="tight", transparent=True)
    print(f"done write {fn}")
