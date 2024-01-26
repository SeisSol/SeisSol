#!/bin/bash
set -euo pipefail

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

$2/preprocessing/science/kinematic_models/get_usgs_finite_fault_data.py $event_id
cd $event_id
read lon lat _ < hypocenter.txt
proj="+proj=tmerc +datum=WGS84 +k=0.9996 +lon_0=${lon} +lat_0=${lat}"
echo $proj
$2/preprocessing/science/kinematic_models/generate_FL33_input_files.py basic_inversion.param --spatial_zoom 3 --generate_ts_yaml "${proj}"
$2/preprocessing/science/kinematic_models/generate_usgs_finite_fault_mesh.py
pumgen -s msh4 mesh.msh
