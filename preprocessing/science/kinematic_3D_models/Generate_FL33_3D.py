"""
Generate FL33 files
With reference to Thomas Ulrich generate_FL33_input_files.py
"""
import os 
import argparse
import numpy as np 
from FaultPlane3D import FaultPlane
import os.path

parser = argparse.ArgumentParser(
    description="Upscale and generate yaml and netCDF input file for fiction law 33/34" +\
    "based on a 3D kinematic model in the standard rupture " +\
    "format .srf file"
)

parser.add_argument("filename", help="filename of the srf file")
parser.add_argument("nx", type=int, help="no. of x nodes in the srf file")
parser.add_argument("ny", type=int, help="no. of y nodes in the srf file")

parser.add_argument(
    "--spatial_order",
    nargs=1,
    metavar=("spatial_order"),
    default=3,
    help="spatial order of the interpolation",
    type=int)

parser.add_argument(
    "--spatial_zoom",
    metavar=("spatial_zoom"),
    required=True,
    help="level of spatial upsampling",
    type=int)

parser.add_argument(
    "--write_paraview",
    dest="write_paraview",
    action="store_true",
    help="write also netcdf readable by paraview",
    default=False)

parser.add_argument(
    "--generate_yaml",
    metavar=("proj"),
    nargs=1,
    help="generate fault yaml file. in this case the proj4 string describing the projection is required.")

parser.add_argument(
    "--Netcdf_nx_ny",
    nargs=2,
    metavar=("nx","ny"),
    required=True,
    type=int,
    help="regular array across the fault as required by Asagi")

parser.add_argument(
    "--area_correction",
    dest="area_correction",
    action="store_true",
    help="Compute area correction for the 3D fault geometry",
    default=True)

parser.add_argument(
    "--plot",
    type=int,
    help="0: no plotting (default) | 1: plot results | 2: plot comparision",
    default=0)

args = parser.parse_args()
proj = args.generate_yaml



###
# 1 Read file and project to XY plane
prefix, ext = os.path.splitext(args.filename)
prefix = os.path.basename(prefix)

fault = FaultPlane()
fault.init_from_srf(args.filename, args.nx, args.ny)
fault.compute_xy_from_latlon(proj)


# 2 Area correction on the fault
if args.area_correction:
    slip1_, slip2_ = fault.slip1.copy(), fault.slip2.copy()

    fault.slip1, fault.slip2, area_factor = fault.area_correction(proj)
    area_sum = sum(sum(area_factor))/(area_factor.shape[0]*area_factor.shape[1])
    area_factor = area_factor/area_sum
    fault.slip1, fault.slip2 = fault.slip1/area_sum, fault.slip2/area_sum


# 3 upscale fault
fault_us = fault.upsample_fault_slip(
    spatial_order=args.spatial_order,
    spatial_zoom=args.spatial_zoom,
    proj=proj,
    verbose=False)

print("####")
ratio_potency1 = np.sum(np.abs(fault_us.slip1)) * fault_us.PSarea_cm2/ (np.sum(np.abs(fault.slip1)) * fault.PSarea_cm2)
ratio_potency2 = np.sum(np.abs(fault_us.slip2)) * fault_us.PSarea_cm2/ (np.sum(np.abs(fault.slip2)) * fault.PSarea_cm2)
print(f"seismic potency ratio - slip direction 1(upscaled and area correction over initial): {ratio_potency1}" )
print(f"seismic potency ratio - slip direction 2(upscaled and area correction over initial): {ratio_potency2}" )


# 4 output yaml and netcdf
fault_us.generate_netcdf_fl33(prefix, args.Netcdf_nx_ny, write_paraview=args.write_paraview)
fault_us.generate_fault_yaml_fl33("Tohoku", prefix)


# 5 Plotting 
if args.plot == 1:
    from FaultPlane_plotting import plot_slip, plot_slip2
    if not os.path.isdir('plot'):
        os.mkdir('plot')
        
    plot_slip(fault, fault_us, savefig=f'plot/{prefix}_upsample.png')
    
elif args.plot == 2:
    print("#####\n start plotting")
    import os.path 
    from FaultPlane_plotting import plot_slip, plot_slip2, plot_netcdf
    from matplotlib import cm
    import matplotlib.pyplot as plt
    
    if not os.path.isdir('plot'):
        os.mkdir('plot')
        
    # Plot 
    cmap_red = cm.get_cmap('Reds')
    a = plt.contourf(fault.x, fault.y, area_factor,cmap=cmap_red,levels=6)
    plt.colorbar(a, label='area factor')
    plt.savefig(f"plot/{prefix}_area_fac.png")

    plt.figure()
    cmap_red = cm.get_cmap('Blues')
    a = plt.contourf(fault.x, fault.y, fault.dip,cmap=cmap_red,levels=6)
    plt.colorbar(a, label='dip')
    plt.savefig(f"plot/{prefix}_dip.png")

    # Plot
    plot_slip2(fault, slip1_, slip2_, fault, savefig=f'plot/{prefix}_area_adjustment')
    plot_slip(fault, fault_us, savefig=f'plot/{prefix}_upsample')
    plot_netcdf(fault_us, savefig=f'plot/{prefix}_NETCDF.png')
