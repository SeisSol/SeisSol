import os
import argparse
from FaultPlane import FaultPlane
import os.path

parser = argparse.ArgumentParser(
    description="generate yaml and netcdf input to be used with friction law 33 based on a (here"
    + "upsampled) kinematic model in the standard rupture format srf file."
)
parser.add_argument("filename", help="filename of the srf file")
parser.add_argument(
    "--spatial_order",
    nargs=1,
    metavar=("spatial_order"),
    default=[3],
    help="spatial order of the interpolation",
    type=int,
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
    "--generate_yaml",
    metavar=("proj"),
    nargs=1,
    help="generate fault yaml file. In this case the proj4 string describing the projection is required",
)
args = parser.parse_args()

p1 = FaultPlane()
p1.init_from_srf(args.filename)
p1.compute_time_array()
p1.assess_Yoffe_parameters()
prefix, ext = os.path.splitext(args.filename)
prefix = os.path.basename(prefix)

p1.generate_netcdf_fl33(
    prefix,
    spatial_order=args.spatial_order[0],
    spatial_zoom=args.spatial_zoom[0],
    proj=args.generate_yaml,
    write_paraview=args.write_paraview,
)
if args.generate_yaml:
    p1.generate_fault_yaml_fl33(
        prefix,
        spatial_order=args.spatial_order[0],
        spatial_zoom=args.spatial_zoom[0],
        proj=args.generate_yaml,
    )
