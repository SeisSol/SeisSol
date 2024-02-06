#!/bin/bash
set -euov pipefail

event_id=$1
if [ $# -ne 2 ]
then
    echo "illegal number of parameters, usage:"
    echo "./workflow_usgs_finite_model.sh event_id path2seissol"
    exit 1
fi

if ! which pumgen > /dev/null; then
    echo "pumgen not in the PATH."
    exit 1
fi

$2/preprocessing/science/automated_workflow_dyn_from_kin/get_usgs_finite_fault_data.py $event_id
cd $event_id
read lon lat _ < tmp/hypocenter.txt
proj="+proj=tmerc +datum=WGS84 +k=0.9996 +lon_0=${lon} +lat_0=${lat}"
echo $proj
$2/preprocessing/science/kinematic_models/generate_FL33_input_files.py tmp/basic_inversion.param --spatial_zoom 3 --generate_ts_yaml "${proj}" 1 --interpolation_method cubic --write_paraview
$2/preprocessing/science/automated_workflow_dyn_from_kin/generate_usgs_finite_fault_mesh.py --fault_mesh_size 500
pumgen -s msh4 tmp/mesh.msh
$2/preprocessing/science/automated_workflow_dyn_from_kin/generate_input_seissol_fl33.py
echo "now run seissol"
