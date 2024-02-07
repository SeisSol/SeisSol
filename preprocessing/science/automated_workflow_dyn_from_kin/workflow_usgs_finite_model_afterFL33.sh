#!/bin/bash
set -euov pipefail

event_id=$1
if [ $# -ne 1 ]
then
    echo "illegal number of parameters, usage:"
    echo "$0 path2seissol"
    exit 1
fi

$1/preprocessing/science/automated_workflow_dyn_from_kin/project_fault_tractions_asagi_grid.py --dx 100 output/fl33-fault.xdmf
$1/preprocessing/science/automated_workflow_dyn_from_kin/infer_segment_average_rake.py output/fl33-fault.xdmf
$1/preprocessing/science/automated_workflow_dyn_from_kin/generate_input_seissol_dr.py
