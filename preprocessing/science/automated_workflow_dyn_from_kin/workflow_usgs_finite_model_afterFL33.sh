#!/bin/bash
set -euov pipefail

if [ $# -ne 1 ]
then
    echo "illegal number of parameters, usage:"
    echo "$0 path2seissol"
    exit 1
fi

$1/preprocessing/science/automated_workflow_dyn_from_kin/project_fault_tractions_asagi_grid.py --dx 200 output/fl33-fault.xdmf
$1/preprocessing/science/automated_workflow_dyn_from_kin/infer_segment_average_rake.py output/fl33-fault.xdmf
$1/preprocessing/science/automated_workflow_dyn_from_kin/generate_input_seissol_dr.py
read lon lat _ < tmp/hypocenter.txt
proj="+proj=tmerc +datum=WGS84 +k=0.9996 +lon_0=${lon} +lat_0=${lat}"
$1/preprocessing/science/kinematic_models/generate_fault_output_from_fl33_input_files.py output/fl33-fault.xdmf yaml_files/FL33_34_fault.yaml --stf Gaussian --output output/dyn-usgs-fault
