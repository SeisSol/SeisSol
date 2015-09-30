#!/bin/sh
##
# @file
# This file is part of SeisSol.
#
# @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
#
# @section LICENSE
# Copyright (c) 2014-2015, SeisSol Group
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
#SBATCH --get-user-env 
# this is for sorting out nodes
##SBATCH --exclude=mac-bdz02,mac-bdz03,mac-bdz09,mac-bdz10,mac-bdz11,mac-bdz14,mac-snb01,mac-snb07,mac-snb08,mac-snb09,mac-snb11,mac-snb13,mac-snb16,mac-snb20,mac-snb26 # weak nodes
#SBATCH --nice=500

# Usage info
show_help() {
cat << EOF
Usage: ${0##*/} [-h] [-d COMPILE_MODE -c CODE_VERSION -a ARCHITECTURE -p PARALLELIZATION -u COMMUNICATION_THREAD -s SCALASCA -q NUMBER_OF_QUANTITIES -o ORDER -n BENCHMARK_DIRECTORY -m NUMBER_OF_MPI_RANKS -r RANKS_PER_NODE -t THREADS_PER_RANK -e EXECUTABLE -f SCRIPTS_DIRECTORY -i ENVIRONMENT]
Submits a Benchmark-job to Slurm.
     -h                      display this help and exit
     -d COMPILE_MODE         release mode.
     -c CODE_VERSION         code version used for execution.
     -a ARCHITECTURE         architecture.
     -p PARALLELIZATION      parallelization.
     -u COMMUNICATION_THREAD communication thread.
     -s SCALASCA             Scalasca.
     -q NUMBER_OF_QUANTITIES number of quantities.
     -o ORDER                order of convergence of the ADER-DG method.
     -n BENCHMARK_DIRECTORY  name of the benchmark.
     -m NUMBER_OF_MPI_RANKS  total number of MPI-ranks.
     -r RANKS_PER_NODE       MPI-ranks used per node.
     -t THREADS_PER_RANK     number of OMP-threads per rank.
     -b THREAD_BINDING       thread pinning
     -e EXECUTABLE           name of the executable.
     -f SCRIPTS_DIRECTORY    directory containing further bash scripts.
     -i ENVIRONMENT          script executed to setup bash-variables
EOF
}

#
# parse command line arguments
#
COMPILE_MODE=NOT_SET
CODE_VERSION=NOT_SET
ARCHITECTURE=NOT_SET
PARALLELIZATION=NOT_SET
COMMUNICATION_THREAD=NOT_SET
SCALASCA=NOT_SET
NUMBER_OF_QUANTITIES=NOT_SET
ORDER=NOT_SET
BENCHMARK_DIRECTORY=NOT_SET
NUMBER_OF_MPI_RANKS=NOT_SET
RANKS_PER_NODE=NOT_SET
THREADS_PER_RANK=NOT_SET
THREAD_BINDING=NOT_SET
EXECUTABLE=NOT_SET
SCRIPTS_DIRECTORY=NOT_SET
ENVIRONMENT=NOT_SET

OPTIND=1
while getopts "hd:c:a:p:u:s:q:o:n:m:r:t:b:e:f:i:" opt; do
    case "$opt" in
        h)
            show_help
            exit 0
            ;;
        d) COMPILE_MODE=$OPTARG
            ;;
        c) CODE_VERSION=$OPTARG
            ;;
        a) ARCHITECTURE=$OPTARG
            ;;
        p) PARALLELIZATION=$OPTARG
            ;;
        u) COMMUNICATION_THREAD=$OPTARG
            ;;
        s) SCALASCA=$OPTARG
            ;;
        q) NUMBER_OF_QUANTITIES=$OPTARG
            ;;
        o) ORDER=$OPTARG
            ;;
        n) BENCHMARK_DIRECTORY=$OPTARG
            echo $BENCHMARK_DIRECTORY
            ;;
        m) NUMBER_OF_MPI_RANKS=$OPTARG
            ;;
        r) RANKS_PER_NODE=$OPTARG
            ;;
        t) THREADS_PER_RANK=$OPTARG
            ;;
        b) THREAD_BINDING=$OPTARG
            ;;
        e) EXECUTABLE=$OPTARG
            ;;
        f) SCRIPTS_DIRECTORY=$OPTARG
            ;;
        i) ENVIRONMENT=$OPTARG
            ;;
        '?')
            show_help >&2
            exit 1
            ;;
    esac
done
shift "$((OPTIND-1))" # Shift off the options and optional --.

if [ "${ENVIRONMENT}" != "NOT_SET" ]
then
  source ${ENVIRONMENT}
fi

export OMP_NUM_THREADS=${THREADS_PER_RANK}
TIME=$(date +"%y_%m_%d-%H_%M_%S")

if [[ $EXECUTABLE == "NOT_SET" ]]
then
EXECUTABLE="SeisSol_${COMPILE_MODE}_${CODE_VERSION}_${ARCHITECTURE}_${PARALLELIZATION}_${COMMUNICATION_THREAD}_${SCALASCA}_${NUMBER_OF_QUANTITIES}_${ORDER}"
fi

echo 'start time: '$(date)
echo 'this are our nodes: '${SLURM_JOB_NODELIST}
echo "Running ${BENCHMARK_DIRECTORY}"

# print information about intel compiler
mpiicc --version
mpiifort --version

#
# Set up the benchmark
#
echo 'moving to benchmark directory'
cd ${BENCHMARK_DIRECTORY} # move to directory
date

# upload data to phis
if [[ ${ARCHITECTURE} == *"knc" ]]
then
  echo /home/${USER}/benchmark/Maple/ > DGPATH

  ${SCRIPTS_DIRECTORY}/phi_communication.sh -m upload -c ${RANKS_PER_NODE} -b ${BENCHMARK_DIRECTORY}
fi

#
# Execute the code
#
echo "executing SeisSol: ${date}"
if [[ ${ARCHITECTURE} != *"knc" ]]
then
  MPIEXEC_ARGUMENTS="-genv I_MPI_DEBUG 4 -genv OMP_NUM_THREADS ${THREADS_PER_RANK} -genv KMP_AFFINITY ${THREAD_BINDING} -ppn ${RANKS_PER_NODE} -n ${NUMBER_OF_MPI_RANKS}"
  if [[ $SCALASCA != "none" ]]
  then
    scalasca -analyze mpiexec.hydra "${MPIEXEC_ARGUMENTS}" ./${EXECUTABLE} parameters.par
  else
    mpiexec.hydra ${MPIEXEC_ARGUMENTS} ./${EXECUTABLE} parameters.par
  fi
else
  # wait for MIC setup
  srun uptime

  mpiexec.hydra -genv I_MPI_DEBUG 4 -genv OMP_NUM_THREADS=240 -genv KMP_AFFINITY proclist=[1-240],granularity=thread,explicit -genv I_MPI_DAPL_PROVIDER ofa-v2-mlx4_0-1 -genv I_MPI_DEVICE rdma -f ${BENCHMARK_DIRECTORY}/mic.${SLURM_JOB_ID}.hosts -ppn 1 -n ${NUMBER_OF_MPI_RANKS} -path /home/${USER}/benchmark -wdir /home/${USER}/benchmark  ./${EXECUTABLE} parameters.par
fi

echo "finished execution, getting the data: $(date)"
# download data from phis
if [[ ${ARCHITECTURE} == *"knc" ]]
then
  ${SCRIPTS_DIRECTORY}/phi_communication.sh -m download -c ${RANKS_PER_NODE} -b ${BENCHMARK_DIRECTORY}
fi

echo "that's it; let's get out of here!"
