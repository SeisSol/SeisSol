#!/bin/bash
set -euov pipefail

if [[ ! -v seissolpath ]]; then
    default=$HOME/SeisSol/SeisSol
    export seissolpath=$default
fi

usgs_event_id_or_dtgeo_file=$1

if [ $# -ne 1 ]
then
    echo "illegal number of parameters, usage:"
    echo "./workflow_usgs_finite_model.sh usgs_event_id_or_dtgeo_file"
    exit 1
fi

folder_name=$($seissolpath/preprocessing/science/automated_workflow_dyn_from_kin/get_usgs_finite_fault_data.py $usgs_event_id_or_dtgeo_file)
cd $folder_name
read lon lat _ < tmp/hypocenter.txt
proj="+proj=tmerc +datum=WGS84 +k=0.9996 +lon_0=${lon} +lat_0=${lat}"
echo $proj
$seissolpath/preprocessing/science/kinematic_models/generate_FL33_input_files.py tmp/basic_inversion.param --spatial_zoom 5 --proj "${proj}"  --interpolation_method cubic
$seissolpath/preprocessing/science/kinematic_models/modify_FL33_34_fault_instantaneous_slip.py

$seissolpath/preprocessing/science/automated_workflow_dyn_from_kin/generate_usgs_finite_fault_mesh.py --fault_mesh_size 500
module load pumgen
pumgen -s msh4 tmp/mesh.msh
$seissolpath/preprocessing/science/automated_workflow_dyn_from_kin/generate_input_seissol_fl33.py
$seissolpath/preprocessing/science/automated_workflow_dyn_from_kin/extract_velocity_model_from_usgs_fsp.py
echo "now run seissol"
