#!/bin/bash
set -euo pipefail
module unload seissol
module load seissol/master-gcc-10.2.1-o4-elas-dunav-single-cuda-7x4s4mk
mkdir -p logs
nvidia-smi
export SEISSOL_PREFERRED_MPI_DATA_TRANSFER_MODE=host
export SEISSOL_MPI_PERSISTENT=1

export OMP_NUM_THREADS=10
export OMP_PLACES="cores"
export OMP_PROC_BIND=spread
counter=0
total_params=$(ls parameters_dyn* | wc -l)

for par_file in parameters_dyn*; do
    counter=$((counter+1))
    echo "Processing file $counter of $total_params: $par_file"
    mpirun -n 2 --map-by ppr:2:numa:pe=$OMP_NUM_THREADS --report-bindings seissol-launch SeisSol_Release_ssm_86_cuda_4_elastic $par_file  2>&1 | tee logs/$par_file.out
done
