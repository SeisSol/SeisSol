import os
import argparse
from FaultPlane import FaultPlane
import os.path

parser = argparse.ArgumentParser(
    description="generate yaml and netcdf input to be used with friction law 33/34 based on a (here"
    + "upsampled) kinematic model in the standard rupture format srf file."
)
parser.add_argument("filename", help="filename of the srf file")
parser.add_argument(
    "--interpolation_method",
    nargs=1,
    metavar=("interpolation_method"),
    default=["linear"],
    help="interpolation method",
    choices=["linear", "nearest", "slinear", "cubic", "quintic"],
)
parser.add_argument(
    "--spatial_zoom",
    nargs=1,
    metavar=("spatial_zoom"),
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
    metavar=("proj"),
    nargs=1,
    help="generate fault yaml file and write fault geometry (ts file). proj: proj4 string describing the projection.",
)

parser.add_argument(
    "--PSRthreshold",
    help="peak slip rate threshold (0-1) to determine onset time and duration of STF",
    nargs=1,
    metavar=("PSRthreshold"),
    type=float,
    default=[0.0],
)

args = parser.parse_args()

p1 = FaultPlane()
p1.init_from_srf(args.filename)
p1.compute_time_array()
p1.assess_STF_parameters(args.PSRthreshold[0])
prefix, ext = os.path.splitext(args.filename)
prefix = os.path.basename(prefix)

p1.generate_netcdf_fl33(
    prefix,
    method=args.interpolation_method[0],
    spatial_zoom=args.spatial_zoom[0],
    proj=args.generate_ts_yaml,
    write_paraview=args.write_paraview,
)
if args.generate_ts_yaml:
    p1.generate_fault_ts_yaml_fl33(
        prefix,
        method=args.interpolation_method[0],
        spatial_zoom=args.spatial_zoom[0],
        proj=args.generate_ts_yaml,
    )
