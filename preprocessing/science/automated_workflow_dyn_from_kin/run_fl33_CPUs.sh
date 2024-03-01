#!/bin/bash
set -euo pipefail
module unload seissol
module load seissol/master-gcc-12.2.0-o4-elas-dunav-single-sggt32q
mkdir -p logs

export SEISSOL_PREFERRED_MPI_DATA_TRANSFER_MODE=host
export SEISSOL_MPI_PERSISTENT=1
export SEISSOL_ASAGI_MPI_MODE=OFF

export OMP_NUM_THREADS=31
export OMP_PLACES="cores"
export OMP_PROC_BIND=spread
export par_file=parameters_fl34.par  

mpirun -n 2 --map-by ppr:2:numa:pe=$OMP_NUM_THREADS SeisSol_Release_srome_4_elastic $par_file >&1 | tee logs/$par_file.out
