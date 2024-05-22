import os
import argparse
from .FaultPlane import FaultPlane


def refine(
    filename,
    projection,
    YoffePSRthreshold,
    spatial_order,
    spatial_zoom,
    temporal_zoom,
    time_smoothing_kernel_as_dt_fraction,
):

    p1 = FaultPlane()
    p1.init_from_srf(filename)
    p1.compute_xy_from_latlon(projection)
    p1.compute_time_array()

    use_Yoffe = True if YoffePSRthreshold else False
    if use_Yoffe:
        p1.assess_STF_parameters(YoffePSRthreshold)

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
