import os
import argparse
from FaultPlane3D import FaultPlane


parser = argparse.ArgumentParser(
    description="Upsample temporally and spatially a kinematic model in the SRF")

parser.add_argument("filename", help="filename of the srf file")
parser.add_argument("--proj", help="transform geometry given proj4 string \
                    as it might be better to upsample the geometry in the local coordinate system")

parser.add_argument("--spatial_order",
                   nargs=1,
                   metavar=("spatial_order"),
                   default=([3]),
                   help="spatial order of the interpolation",
                   type=int)

parser.add_argument("--spatial_zoom",
                   nargs=1,
                   metavar=("spatial_zoom"),
                   default=([3]),
                   help="level of spatial upsampling",
                   type=int)

parser.add_argument("--temporal_zoom",
                   nargs=1,
                   metavar=("temproal_zoom"),
                   default=([3]),
                   help="level of temporal upsampling",
                   type=int)

parser.add_argument("--time_smoothing_kernal_as_dt_fraction",
                   nargs=1,
                   metavar=("alpha_dt"),
                   default=([0.5]),
                   help="sigma, expressed as a portion of dt, of hte gaussian keranl to smooth SR",
                   type=float)

parser.add_argument("--use_Yoffee",
                   help="replace the discretized STF with a Ypffe function (e.g. for comparison with FL33)",
                   dest="use_Yoffee",
                   action="store_true")

args = parser.args()

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
    time_smoothing_kernal_as_dt_fraction=args.time_smoothing_kernal_as_dt_fraction[0])

prefix, ext = os.path.splitext(args.filename)
fnout = prefix + "_resampled" + ".srf"
p2.write_srf(fnout)