#!/usr/bin/env python3
import os
import argparse
from FaultPlane import FaultPlane, MultiFaultPlane
import os.path

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
parser.add_argument(
    "--generate_ts_yaml",
    metavar=("proj", "instantaneous"),
    nargs=2,
    help=(
        "generate fault yaml file and write fault geometry (ts file). proj: proj4"
        " string describing the projection. instantaneous: (0 or 1) to apply the slip"
        " over 0.1s everywhere or not"
    ),
)

parser.add_argument(
    "--PSRthreshold",
    help="peak slip rate threshold (0-1) to determine onset time and duration of STF",
    nargs=1,
    metavar="PSRthreshold",
    type=float,
    default=[0.0],
)

args = parser.parse_args()

prefix, ext = os.path.splitext(args.filename)
prefix = os.path.basename(prefix)

if ext == ".srf":
    mfp = MultiFaultPlane.from_srf(args.filename)
elif ext == ".param":
    mfp = MultiFaultPlane.from_usgs_param_file(args.filename)
elif ext == ".txt":
    mfp = MultiFaultPlane.from_slipnear_param_file(args.filename)
else:
    raise NotImplementedError(f" unknown extension: {ext}")

for p, p1 in enumerate(mfp.fault_planes):
    p1.compute_time_array()
    if ext == ".srf":
        p1.assess_STF_parameters(args.PSRthreshold[0])
    p1.generate_netcdf_fl33(
        f"{prefix}{p+1}",
        method=args.interpolation_method[0],
        spatial_zoom=args.spatial_zoom[0],
        proj=args.generate_ts_yaml[0],
        write_paraview=args.write_paraview,
    )

if args.generate_ts_yaml:
    mfp.generate_fault_ts_yaml_fl33(
        prefix,
        method=args.interpolation_method[0],
        spatial_zoom=args.spatial_zoom[0],
        proj=args.generate_ts_yaml[0],
        instantaneous=True if args.generate_ts_yaml[1] else False,
    )
