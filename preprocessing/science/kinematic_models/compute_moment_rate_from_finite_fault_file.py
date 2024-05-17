#!/usr/bin/env python3
import easi
import os
import argparse
from FaultPlane import FaultPlane, MultiFaultPlane
import os.path
import numpy as np
from GaussianSTF import GaussianSTF, SmoothStep

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

prefix, ext = os.path.splitext(args.filename)
prefix = os.path.basename(prefix)

if ext == ".srf":
    mfp = MultiFaultPlane.from_srf(args.filename)
elif ext == ".param":
    mfp = MultiFaultPlane.from_usgs_param_file(args.filename)
elif ext == ".fsp":
    mfp = MultiFaultPlane.from_usgs_fsp_file(args.filename)
elif ext == ".txt":
    mfp = MultiFaultPlane.from_slipnear_param_file(args.filename)
else:
    raise NotImplementedError(f" unknown extension: {ext}")

dt = args.dt[0]
duration = 0
for p, fp in enumerate(mfp.fault_planes):
    duration = max(duration, fp.t0.max() + fp.rise_time.max())

time = np.arange(0, duration, dt)
moment_rate = np.zeros_like(time)

for p, fp in enumerate(mfp.fault_planes):
    fp.compute_xy_from_latlon(args.proj[0])
    centers = np.column_stack(
        (fp.x.flatten(), fp.y.flatten(), -fp.depth.flatten() * 1e3)
    )
    tags = np.zeros_like(centers[:, 0]) + 1
    out = easi.evaluate_model(
        centers,
        tags,
        ["mu"],
        args.yaml_filename,
    )
    mu = out["mu"].reshape(fp.x.shape)
    for k, tk in enumerate(time):
        STF = GaussianSTF(tk - fp.t0[:, :], fp.rise_time[:, :], dt)
        for j in range(fp.ny):
            for i in range(fp.nx):
                moment_rate[k] += (
                    mu[j, i] * fp.dx * fp.dy * 1e6 * STF[j, i] * fp.slip1[j, i] * 0.01
                )
M0 = np.trapz(moment_rate[:], x=time[:])
Mw = 2.0 * np.log10(M0) / 3.0 - 6.07
print(f"inferred Mw {Mw} and duration {duration}:")

if not os.path.exists("tmp"):
    os.makedirs("tmp")
fname = "tmp/moment_rate_from_finite_source_file.txt"
with open(fname, "w") as f:
    np.savetxt(f, np.column_stack((time, moment_rate)), fmt="%g")
print(f"done writing {fname}")
