import os
import argparse
from FaultPlane import FaultPlane

parser = argparse.ArgumentParser(
    description="upsample temporally and spatially a kinematic model (should be a planar model) in the standard rupture format"
)
parser.add_argument("filename", help="filename of the srf file")
parser.add_argument(
    "--proj",
    help="proj4 string (might be better to upsample the geometry in the local coordinate system)",
)
parser.add_argument(
    "--spatial_order",
    nargs=1,
    metavar=("spatial_order"),
    default=([3]),
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
    "--temporal_zoom",
    nargs=1,
    metavar=("temporal_zoom"),
    required=True,
    help="level of temporal upsampling",
    type=int,
)

parser.add_argument(
    "--use_Yoffe",
    help="replace the discretized STF with a Yoffe function (e.g. for comparison with FL33)",
    dest="use_Yoffe",
    action="store_true",
)

args = parser.parse_args()

p1 = FaultPlane()
p1.init_from_srf(args.filename)
p1.compute_xy_from_latlon(args.proj)
p1.compute_time_array()

p2 = p1.upsample_fault(
    spatial_order=args.spatial_order[0],
    spatial_zoom=args.spatial_zoom[0],
    temporal_zoom=args.temporal_zoom[0],
    proj=args.proj,
    use_Yoffe=args.use_Yoffe,
)
prefix, ext = os.path.splitext(args.filename)
fnout = prefix + "_resampled" + ".srf"
p2.write_srf(fnout)
