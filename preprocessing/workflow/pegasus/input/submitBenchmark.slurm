#!/bin/sh
##
# @file
# This file is part of SeisSol.
#
# @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
#
# @section LICENSE
# Copyright (c) 2014, SeisSol Group
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
##SBATCH -o tpv16.%j.%N.out                   # command line parameter
##SBATCH -J tpv16                             # command line parameter
#SBATCH --get-user-env 
##SBATCH --partition=snb                      # command line parameter
##SBATCH --constraint=turbo_off               # command line parameter
##SBATCH --ntasks=16                          # command line parameter
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=end                                             
#SBATCH --mail-user=breuera@in.tum.de
#SBATCH --export=NONE                                               
##SBATCH --time=10:00:00                      # command line parameter
#SBATCH --nice=500

. /etc/profile 
. /etc/profile.d/modules.sh

# Usage info
show_help() {
cat << EOF
Usage: ${0##*/} [-h] [-d COMPILE_MODE -c CODE_VERSION -v VECTOR_INSTRUCTION_SET -p PARALLELIZATION -s SCALASCA -q NUMBER_OF_QUANTITIES -b NUMBER_OF_BASIS_FUNCTIONS -n BENCHMARK_NAME -m NUMBER_OF_MPI_RANKS -r RANKS_PER_NODE -t THREADS_PER_RANK]
Submits a Benchmark-job to Slurm.
     -h                           display this help and exit
     -d COMPILE_MODE              release mode.
     -c CODE_VERSION              code version used for execution.
     -v VECTOR_INSTRUCTION_SET    vector instruction set.
     -p PARALLELIZATION           parallelization.
     -s SCALASCA                  Scalasca.
     -q NUMBER_OF_QUANTITIES      number of quantities.
     -b NUMBER_OF_BASIS_FUNCTIONS number of basis functions used.
     -n BENCHMARK_NAME            name of the benchmark.
     -m NUMBER_OF_MPI_RANKS       total number of MPI-ranks.
     -r RANKS_PER_NODE            MPI-ranks used per node.
     -t THREADS_PER_RANK          number of OMP-threads per rank.
EOF
}

#
# parse command line arguments
#
COMPILE_MODE=NOT_SET
CODE_VERSION=NOT_SET
VECTOR_INSTRUCTION_SET=NOT_SET
PARALLELIZATION=NOT_SET
SCALASCA=NOT_SET
NUMBER_OF_QUANTITIES=NOT_SET
NUMBER_OF_BASIS_FUNCTIONS=NOT_SET
BENCHMARK_NAME=NOT_SET
NUMBER_OF_MPI_RANKS=NOT_SET
RANKS_PER_NODE=NOT_SET
THREADS_PER_RANK=NOT_SET

OPTIND=1
while getopts "hd:c:v:p:s:q:b:n:m:r:t:" opt; do
    case "$opt" in
        h)
            show_help
            exit 0
            ;;
        d) COMPILE_MODE=$OPTARG
            ;;
        c) CODE_VERSION=$OPTARG
            ;;
        v) VECTOR_INSTRUCTION_SET=$OPTARG
            ;;
        p) PARALLELIZATION=$OPTARG
            ;;
        s) SCALASCA=$OPTARG
            ;;
        q) NUMBER_OF_QUANTITIES=$OPTARG
            ;;
        b) NUMBER_OF_BASIS_FUNCTIONS=$OPTARG
            ;;
        n) BENCHMARK_NAME=$OPTARG
            echo $BENCHMARK_NAME
            ;;
        m) NUMBER_OF_MPI_RANKS=$OPTARG
            ;;
        r) RANKS_PER_NODE=$OPTARG
            ;;
        t) THREADS_PER_RANK=$OPTARG
            ;;
        '?')
            show_help >&2
            exit 1
            ;;
    esac
done
shift "$((OPTIND-1))" # Shift off the options and optional --.
 
#
# Set up Phi-speicfic environment
#
if [ ${VECTOR_INSTRUCTION_SET} == "mic" ]
then
  echo 'loading custom toolchain for the phi-version'
  module unload mpi.intel
  module unload ccomp
  module unload fortran
  module unload mkl
  module load gcc/4.6
  source /lrz/sys/intel/compiler140_106/composer_xe_2013_sp1.1.106/bin/compilervars.sh intel64
  source /lrz/sys/intel/mpi_41_1_036/intel64/bin/mpivars.sh
fi

#
# Setup variables
#
if [[ $SCALASCA != "none" ]]
then
SCALASCA_WRAPPER="scalasca -analyze"
fi

# Should be parent of trunk
STARTDIR=$(pwd)
export OMP_NUM_THREADS=${THREADS_PER_RANK}
BENCHMARK_DIR="Benchmark_${BENCHMARK_NAME}_${NUMBER_OF_BASIS_FUNCTIONS}"
TIME=$(date +"%y_%m_%d-%H_%M_%S")

EXECUTABLE="SeisSol_${COMPILE_MODE}_${CODE_VERSION}_${VECTOR_INSTRUCTION_SET}_${PARALLELIZATION}_${SCALASCA}_${NUMBER_OF_QUANTITIES}_${NUMBER_OF_BASIS_FUNCTIONS}"

echo 'start time: '$(date)
echo "Running ${BENCHMARK_NAME}: ${BENCHMARK_DIR}"

# print information about intel compiler
mpiicc --version
mpiifort --version

#
# Set up the benchmark
#
echo 'setting up benchmark'
cd ${BENCHMARK_DIR} # move to directory
echo $STARTDIR/trunk/Maple/ > DGPATH # create DGPATH-FILE
ln -s ${STARTDIR}/trunk/preprocessing/matrices/matrices_*.xml . # link matrices
ln -s ${STARTDIR}/trunk/build/${EXECUTABLE} # link executable
mkdir output # create output folder

echo 'ready, waiting 30 seconds for slow file systems'
date
sleep 30

# upload data to phis
if [ ${VECTOR_INSTRUCTION_SET} == "mic" ]
then
  ln -s ${STARTDIR}/trunk/Maple
  echo /home/${USER}/benchmark/Maple/ > DGPATH
  ${STARTDIR}/phi_communication -m upload -c ${RANKS_PER_NODE} -b ${STARTDIR}/${BENCHMARK_DIR}
fi

#
# Execute the code
#
echo "executing SeisSol: ${date}"
if [ ${VECTOR_INSTRUCTION_SET} != "mic" ]
then
  $SCALASCA_WRAPPER mpiexec.hydra -ppn ${RANKS_PER_NODE} -n ${NUMBER_OF_MPI_RANKS} ./${EXECUTABLE} parameters.par
else
  # load Phi toolchain
  module unload mpi.intel
  module unload ccomp
  module unload fortran
  module unload mkl
  module load gcc/4.6
  source /lrz/sys/intel/compiler140_106/composer_xe_2013_sp1.1.106/bin/compilervars.sh intel64
  source /lrz/sys/intel/mpi_41_1_036/intel64/bin/mpivars.sh
  export I_MPI_MIC=1

  export I_MPI_FABRICS=shm:dapl
  export I_MPI_HYDRA_BOOTSTRAP=ssh

  # wait for MIC setup
  srun uptime

  mpiexec.hydra -genv I_MPI_DEBUG 4 -genv OMP_NUM_THREADS=240 -genv KMP_AFFINITY proclist=[1-240],granularity=thread,explicit -genv I_MPI_DAPL_PROVIDER ofa-v2-mlx4_0-1 -genv I_MPI_DEVICE rdma -f mic.${SLURM_JOB_ID}.hosts -ppn 1 -n ${NUMBER_OF_MPI_RANKS} -path /home/${USER}/benchmark -wdir /home/${USER}/benchmark  ./${EXECUTABLE} parameters.par
fi

echo "finished execution, getting the data: $(date)"
# download data from phis
if [ ${VECTOR_INSTRUCTION_SET} == "mic" ]
then
  ${STARTDIR}/phi_communication -m download -c ${RANKS_PER_NODE} -b ${STARTDIR}/${BENCHMARK_DIR}
fi

echo "that's it; let's get out of here!"
