#!/usr/bin/env python3
import os
import argparse
from .FaultPlane import FaultPlane, MultiFaultPlane
import os.path

def main(
    filename,
    interpolation_method,
    spatial_zoom,
    projection,
    write_paraview,
    PSRthreshold,
):
    prefix, ext = os.path.splitext(filename)
    prefix = os.path.basename(prefix)

    if ext == ".srf":
        mfp = MultiFaultPlane.from_srf(filename)
    elif ext == ".param":
        mfp = MultiFaultPlane.from_usgs_param_file(filename)
    elif ext == ".fsp":
        mfp = MultiFaultPlane.from_usgs_fsp_file(filename)
    elif ext == ".txt":
        mfp = MultiFaultPlane.from_slipnear_param_file(filename)
    else:
        raise NotImplementedError(f" unknown extension: {ext}")

    for p, p1 in enumerate(mfp.fault_planes):
        p1.compute_time_array()
        if ext == ".srf":
            p1.assess_STF_parameters(PSRthreshold)
        p1.generate_netcdf_fl33(
            f"{prefix}{p+1}",
            method=interpolation_method,
            spatial_zoom=spatial_zoom,
            proj=projection,
            write_paraview=write_paraview,
            slip_cutoff=0.0,
        )

    mfp.generate_fault_ts_yaml_fl33(
        prefix,
        method=interpolation_method,
        spatial_zoom=spatial_zoom,
        proj=projection,
    )
