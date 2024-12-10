#!/usr/bin/env python3
import os
import argparse
from FaultPlane import FaultPlane, MultiFaultPlane
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "generate yaml and netcdf input to be used with friction law 33/34 based on a"
            " (here upsampled) kinematic model in the standard rupture format srf file."
        )
    )
    parser.add_argument("filename", help="filename of the srf file")

    parser.add_argument(
        "--interpolation_method",
        nargs=1,
        metavar="interpolation_method",
        default=["linear"],
        help="interpolation method",
        choices=["linear", "nearest", "slinear", "cubic", "quintic"],
    )
    parser.add_argument(
        "--proj",
        metavar=("proj"),
        nargs=1,
        help=("proj4 string describing the projection"),
        required=True,
    )
    parser.add_argument(
        "--PSRthreshold",
        help="peak slip rate threshold (0-1) to determine onset time and duration of STF",
        nargs=1,
        metavar="PSRthreshold",
        type=float,
        default=[0.0],
    )
    parser.add_argument(
        "--spatial_zoom",
        nargs=1,
        metavar="spatial_zoom",
        required=True,
        help="level of spatial upsampling",
        type=int,
    )
    parser.add_argument(
        "--write_paraview",
        dest="write_paraview",
        action="store_true",
        help="write also netcdf readable by paraview",
        default=False,
    )
    args = parser.parse_args()
    main(
        args.filename,
        args.interpolation_method[0],
        args.spatial_zoom[0],
        args.proj[0],
        args.write_paraview,
        args.PSRthreshold[0],
    )
