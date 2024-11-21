#!/usr/bin/env python3
import easi
import seissolxdmf
import seissolxdmfwriter as sxw
import argparse
import numpy as np
from stf import regularizedYoffe, gaussianSTF, smoothStep, asymmetric_cosine
from tqdm import tqdm
import argparse


class seissolxdmfExtended(seissolxdmf.seissolxdmf):
    def __init__(self, xdmfFilename):
        super().__init__(xdmfFilename)
        self.xyz = self.ReadGeometry()
        self.connect = self.ReadConnect()

    def ReadFaultTags(self):
        """Read partition array"""
        return self.Read1dData("fault-tag", self.nElements, isInt=True).T

    def ComputeCellCenters(self):
        """compute cell center array"""
        return (
            self.xyz[self.connect[:, 0]]
            + self.xyz[self.connect[:, 1]]
            + self.xyz[self.connect[:, 2]]
        ) / 3.0


def generate(fault_filename, yaml_filename, output_file, stf, dt_output):
    sx = seissolxdmfExtended(fault_filename)
    centers = sx.ComputeCellCenters()
    tags = sx.ReadFaultTags()

    print(f"using {stf}")
    if stf == "Yoffe":
        out = easi.evaluate_model(
            centers,
            tags,
            ["strike_slip", "dip_slip", "tau_S", "tau_R", "rupture_onset"],
            yaml_filename,
        )
        ts = out["tau_S"]
        tr = out["tau_R"]
        ts = np.maximum(ts, 0)
        tr = np.maximum(tr, ts)
        rise_time = tr + 2.0 * ts
        dt_required = 1.27 * np.amin(ts) / 3.0
        n_time_sub = max(1, round(dt_output / dt_required))
        dt = dt_output / n_time_sub
        print(f"STF numerically integrated (dt used {dt})")
    elif stf == "AsymmetricCosine":
        out = easi.evaluate_model(
            centers,
            tags,
            ["strike_slip", "dip_slip", "tau_S", "rupture_rise_time", "rupture_onset"],
            yaml_filename,
        )
        acc_time = np.maximum(out["tau_S"], 0) * 1.27
        rise_time = out["rupture_rise_time"]
        dt_required = np.amin(acc_time) / 3.0
        n_time_sub = max(1, round(dt_output / dt_required))
        dt = dt_output / n_time_sub
        print(f"STF numerically integrated (dt used {dt})")
    elif stf == "Gaussian":
        n_time_sub = 1
        dt = dt_output
        out = easi.evaluate_model(
            centers,
            tags,
            ["strike_slip", "dip_slip", "rupture_rise_time", "rupture_onset"],
            yaml_filename,
        )
        rise_time = out["rupture_rise_time"]

    onset = out["rupture_onset"]
    onset[onset > 1e50] = np.nan
    rise_time[rise_time > 1e50] = np.nan

    tmax = np.nanmax(rise_time + onset) + 2.0
    print(tmax)

    time = np.arange(0, tmax, dt)
    time_out = time[::n_time_sub]

    ndt_output = time_out.shape[0]
    nel = sx.ReadNElements()
    STF = np.zeros((nel,))
    intSTF = np.zeros((nel,))

    ASl = np.zeros((ndt_output, nel))
    SR = np.zeros((ndt_output, nel))
    Sls = np.zeros((ndt_output, nel))
    Sld = np.zeros((ndt_output, nel))

    sls = out["strike_slip"]
    sld = out["dip_slip"]
    slip = np.sqrt(sls**2 + sld**2)
    print(time)
    for k, ti in enumerate(tqdm(time)):
        if stf == "Gaussian":
            STF = gaussianSTF(ti - onset, rise_time, dt_output)
            intSTF = smoothStep(ti - onset, rise_time)
        elif stf == "AsymmetricCosine":
            STF = asymmetric_cosine(ti - onset, acc_time, rise_time - acc_time)
            intSTF += dt * STF
        elif stf == "Yoffe":
            for i in range(nel):
                STF[i] = regularizedYoffe(ti - onset[i], ts[i], tr[i])
            intSTF += dt * STF
        if k % n_time_sub == 0:
            p = k // n_time_sub
            SR[p, :] = STF * slip
            ASl[p, :] = intSTF * slip
            Sls[p, :] = intSTF * sls
            Sld[p, :] = intSTF * sld
    assert ASl.max() > 0
    if not stf == "Gaussian":
        id_where_slip = np.where(slip > 0.02 * slip)[0]
        error_slip = np.abs(ASl[-1, id_where_slip] - slip[id_where_slip])
        error_slip_rel = error_slip / slip[id_where_slip]
        print(
            "max error/relative on slip due to numerical integration"
            f" {np.amax(error_slip)} {np.amax(error_slip_rel)}"
        )
    dictTime = {time_out[i]: i for i in range(time_out.shape[0])}

    sxw.write(
        output_file,
        sx.xyz,
        sx.connect,
        {"ASl": ASl, "SR": SR, "Sls": Sls, "Sld": Sld},
        dictTime,
        reduce_precision=True,
        backend="raw",
    )


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
