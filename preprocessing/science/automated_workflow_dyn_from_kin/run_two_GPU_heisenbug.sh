#!/bin/bash
set -euov pipefail
mkdir -p logs
nvidia-smi
mem_usage=$(nvidia-smi --query-gpu=memory.used --format=csv | awk ' {s+=$1} END {printf "%.0f\n", s}')
export SEISSOL_PREFERRED_MPI_DATA_TRANSFER_MODE=host
export SEISSOL_MPI_PERSISTENT=1

if [ $mem_usage == '0' ]; then
   export OMP_NUM_THREADS=10
   export OMP_PLACES="cores"
   export OMP_PROC_BIND=spread

   for par_file in parameters_dyn*; do
       echo "Processing file: $par_file"
       mpirun -n 2 --map-by ppr:1:numa:pe=10 --report-bindings seissol-launch SeisSol_Release_ssm_86_cuda_4_elastic $par_file  2>&1 | tee logs/$par_file.out
   done
else
    echo "GPU already in use!"
fi
