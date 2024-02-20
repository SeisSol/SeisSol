#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import argparse
import matplotlib
import os
import glob
import re


def computeMw(label, time, moment_rate):
    M0 = np.trapz(moment_rate[:], x=time[:])
    Mw = 2.0 * np.log10(M0) / 3.0 - 6.07
    # print(f"{label} moment magnitude: {Mw:.2} (M0 = {M0:.4e})")
    return Mw


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
results = {"B": [], "C": [], "R0": [], "Mw": []}

for i, fn in enumerate(energy_files):
    df = pd.read_csv(fn)
    df = df.pivot_table(index="time", columns="variable", values="measurement")
    df["seismic_moment_rate"] = np.gradient(
        df["seismic_moment"], df.index[1] - df.index[0]
    )
    label = os.path.basename(fn)
    B, C, R = extractBCR(fn)
    Mw = computeMw(label, df.index.values, df["seismic_moment_rate"])
    results["Mw"].append(Mw)
    results["B"].append(B)
    results["C"].append(C)
    results["R0"].append(R)
    label = f"B={B}, C={C}, R={R}"
    ax.plot(
        df.index.values,
        df["seismic_moment_rate"] / 1e19,
        label=f"{label} (Mw={Mw:.2f})",
        linestyle="--" if C == 0.3 else "-",
    )
mr_usgs = np.loadtxt("tmp/moment_rate.mr", skiprows=2)
ax.plot(mr_usgs[:, 0], mr_usgs[:, 1] / 1e19, label="usgs", color="black")

result_df = pd.DataFrame(results)
print(result_df)

selected_rows = result_df[result_df["Mw"] > 6]
print(selected_rows)

ax.legend(frameon=False, loc="upper right", ncol=2, fontsize=6)
ax.set_ylim(bottom=0)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

ax.set_ylabel(r"moment rate (e19 $\times$ Nm/s)")
ax.set_xlabel("time (s)")

fn = f"plots/moment_rate.pdf"
fig.savefig(fn, bbox_inches="tight", transparent=True)
print(f"done write {fn}")
