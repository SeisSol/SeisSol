#!/bin/sh

proj="+proj=tmerc +datum=WGS84 +k=0.9996 +lon_0=143 +lat_0=39"

# List of srf and corresponding Nx Ny nodes are in SRF list
tail -n +2 srf/srf_list | while read line
do
    file=`echo $line| awk '{print $1}'`
    nx=`echo $line| awk '{print $2}'` 
    ny=`echo $line| awk '{print $3}'`

   echo $file $nx $ny 

   python Generate_FL33_3D.py srf/$file $nx $ny \
        --spatial_zoom 4 \
        --generate_yaml "+proj=tmerc +datum=WGS84 +k=0.9996 +lon_0=143 +lat_0=39"  \
        --Netcdf_nx_ny 100 100 \
        --write_paraview \
        --area_correction  
done

