import os
import argparse
from FaultPlane import FaultPlane

parser = argparse.ArgumentParser(
    description="upsample temporally and spatially a kinematic model (should consist of only one segment) in the standard rupture format"
)
parser.add_argument("filename", help="filename of the srf file")
parser.add_argument(
    "--proj",
    help="transform geometry given proj4 string (as it might be better to upsample the geometry in the local coordinate system)",
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
    "--time_smoothing_kernel_as_dt_fraction",
    nargs=1,
    metavar=("alpha_dt"),
    default=([0.5]),
    help="sigma, expressed as a portion of dt, of the gaussian kernel used to smooth SR",
    type=float,
)

parser.add_argument(
    "--use_Yoffe",
    help="replace the discretized STF with a Yoffe function (e.g. for comparison with FL33).\
       Requires peak slip rate threshold (0-1) to determine onset time and duration of STF",
    dest="use_Yoffe",
    nargs=1,
    metavar=("PSRthreshold"),
    type=float,
)

args = parser.parse_args()

p1 = FaultPlane()
p1.init_from_srf(args.filename)
p1.compute_xy_from_latlon(args.proj)
p1.compute_time_array()

use_Yoffe = True if args.use_Yoffe else False
if use_Yoffe:
    p1.assess_STF_parameters(args.use_Yoffe[0])

p2 = p1.upsample_fault(
    spatial_order=args.spatial_order[0],
    spatial_zoom=args.spatial_zoom[0],
    temporal_zoom=args.temporal_zoom[0],
    proj=args.proj,
    use_Yoffe=use_Yoffe,
    time_smoothing_kernel_as_dt_fraction=args.time_smoothing_kernel_as_dt_fraction[0],
)
prefix, ext = os.path.splitext(args.filename)
fnout = prefix + "_resampled" + ".srf"
p2.write_srf(fnout)
