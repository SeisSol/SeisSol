import os
import argparse
from FaultPlane import FaultPlane

def refine(
    filename,
    projection,
    yoffePSRthreshold,
    spatial_order,
    spatial_zoom,
    temporal_zoom,
    time_smoothing_kernel_as_dt_fraction,
):

    p1 = FaultPlane()
    p1.init_from_srf(filename)
    p1.compute_xy_from_latlon(projection)
    p1.compute_time_array()

    use_Yoffe = True if yoffePSRthreshold else False
    if use_Yoffe:
        p1.assess_STF_parameters(yoffePSRthreshold)

    p2 = p1.upsample_fault(
        spatial_order=spatial_order,
        spatial_zoom=spatial_zoom,
        temporal_zoom=temporal_zoom,
        proj=projection,
        use_Yoffe=use_Yoffe,
        time_smoothing_kernel_as_dt_fraction=time_smoothing_kernel_as_dt_fraction,
    )
    prefix, ext = os.path.splitext(args.filename)
    fnout = prefix + "_resampled" + ".srf"
    p2.write_srf(fnout)


if __name__ == "__main__":
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
        default=[None]
    )

    args = parser.parse_args()


    yoffePSRthreshold= args.use_Yoffe[0]

    refine(
        ags.filename,
        args.proj,
        yoffePSRthreshold,
        args.spatial_order[0],
        args.spatial_zoom[0],
        args.temporal_zoom[0],
        args.time_smoothing_kernel_as_dt_fraction[0],
    )
