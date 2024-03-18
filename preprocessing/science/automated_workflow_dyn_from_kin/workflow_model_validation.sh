#!/bin/bash
set -euov pipefail

if [[ ! -v seissolpath ]]; then
    default=$HOME/SeisSol/SeisSol
    export seissolpath=$default
fi

if [[ ! -v teleseismicspath ]]; then
    default=$HOME/TuSeisSolScripts/TeleseismicDataRelated/
    export teleseismicspath=$default
fi


$seissolpath/preprocessing/science/automated_workflow_dyn_from_kin/compile_scenario_macro_properties.py output

read lon lat _ < tmp/hypocenter.txt
proj="+proj=tmerc +datum=WGS84 +k=0.9996 +lon_0=${lon} +lat_0=${lat}"

for faultoutput in output/dyn*-fault.xdmf; do
   echo $faultoutput
   $teleseismicspath/compute_multi_cmt.py spatial $faultoutput  1 $seissolpath/preprocessing/science/automated_workflow_dyn_from_kin/input_files/smooth_PREM.txt --DH 10 --proj "${proj}" --NZ 3
done
mv PointSource* tmp
$seissolpath/preprocessing/science/automated_workflow_dyn_from_kin/generate_teleseismic_config_from_usgs.py
$teleseismicspath/select_stations_azimuthal.py teleseismic_config.ini 5
$teleseismicspath/generate_figure_teleseismic_synthetics.py teleseismic_config.ini
