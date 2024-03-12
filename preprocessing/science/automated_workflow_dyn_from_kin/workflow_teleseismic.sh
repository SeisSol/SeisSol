#!/bin/bash
set -euov pipefail

if [ $# -ne 2 ]
then
    echo "illegal number of parameters, usage:"
    echo "$0 path2seissol path2teleseismicScripts"
    exit 1
fi

read lon lat _ < tmp/hypocenter.txt
proj="+proj=tmerc +datum=WGS84 +k=0.9996 +lon_0=${lon} +lat_0=${lat}"
for faultoutput in output/dyn*-fault.xdmf; do
   echo $faultoutput
   #todo change to easi
   $2/compute_multi_cmt.py spatial $faultoutput  0 30e9 --DH 10 --proj "${proj}"
done
mv PointSource* tmp
# todo: adapt duration Tmin and Tmax
$1/preprocessing/science/automated_workflow_dyn_from_kin/generate_teleseismic_config_from_usgs.py
$2/select_stations_azimuthal.py teleseismic_config.ini 6
$2/generate_figure_teleseismic_synthetics.py teleseismic_config.ini
