#!/bin/bash
set -euov pipefail

if [[ ! -v seissolpath ]]; then
    default=$HOME/SeisSol/SeisSol
    export seissolpath=$default
fi

$seissolpath/preprocessing/science/kinematic_models/generate_fault_output_from_fl33_input_files.py output/fl33-fault.xdmf yaml_files/FL33_34_fault.yaml --stf Gaussian --output output/dyn-usgs-fault
$seissolpath/preprocessing/science/automated_workflow_dyn_from_kin/project_fault_tractions_asagi_grid.py --dx 200 output/fl33-fault.xdmf
$seissolpath/preprocessing/science/automated_workflow_dyn_from_kin/infer_segment_average_rake.py output/fl33-fault.xdmf
module load seissol
$seissolpath/preprocessing/science/automated_workflow_dyn_from_kin/generate_input_seissol_dr.py
read lon lat _ < tmp/hypocenter.txt
proj="+proj=tmerc +datum=WGS84 +k=0.9996 +lon_0=${lon} +lat_0=${lat}"
