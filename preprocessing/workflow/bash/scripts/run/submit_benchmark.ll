#!/bin/sh
##
#
# This file is part of SeisSol.
#
# author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
#
# author Leonhard Rannabauer (lrannabauer AT mytum.de)
#
# section LICENSE
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
#@ job_type = parallel
#@ job_name = {job_name}
#@ output = {output}
#@ error = {output}
#@ class = {queue}
#@ node = {numberofnodes}
#@ tasks_per_node = {ranks_per_node}
#@ wall_clock_limit = {time}
#@ queue

. /etc/profile 
. /etc/profile.d/modules.sh
source ${HOME}/tmp_env.sh

echo ${WORKING_DIRECTORY}
echo ${INPUT_DIRECTORY}
echo ${SCRIPTS_DIRECTORY}

MIC={mic}

#
# parse command line arguments
#
COMPILE_MODE={compile_mode}
CODE_VERSION={code_version}
ARCHITECTURE={architecture}
PARALLELIZATION={parallelization}
COMMUNICATION_THREAD={communication_thread}
SCALASCA={scalasca}
NUMBER_OF_QUANTITIES={number_of_quantities}
ORDER={order}
BENCHMARK_DIRECTORY={benchmark_directory}
NUMBER_OF_MPI_RANKS={number_of_mpi_ranks}
RANKS_PER_NODE={ranks_per_node}
THREADS_PER_RANK={threads_per_rank}
EXECUTABLE={executable}
SCRIPTS_DIRECTORY={scripts_directory}

#
# Setup variables
#
if [[ $SCALASCA != "none" ]]
then
SCALASCA_WRAPPER="scalasca -analyze"
fi

# setup omp
export MP_PE_AFFINITY=no
export OMP_NUM_THREADS=${THREADS_PER_RANK}

if [[ ${ARCHITECTURE} == *"hsw" ]]
then
  export KMP_AFFINITY=compact,granularity=thread
fi

TIME=$(date +"%y_%m_%d-%H_%M_%S")

if [[ $EXECUTABLE == "NOT_SET" ]]
then
EXECUTABLE="SeisSol_${COMPILE_MODE}_${CODE_VERSION}_${ARCHITECTURE}_${PARALLELIZATION}_${COMMUNICATION_THREAD}_${SCALASCA}_${NUMBER_OF_QUANTITIES}_${ORDER}"
fi

echo 'start time: '$(date)

# List nodes
echo "listing nodes"
cat ${LOADL_HOSTFILE}
echo

echo "Running ${BENCHMARK_DIRECTORY}"

echo "getting some rpm in those fans"
/lrz/sys/supermic/bin/fan-max

# print information about intel compiler
icc --version
ifort --version

#
# Set up the benchmark
#
echo 'moving to benchmark directory'
cd ${BENCHMARK_DIRECTORY} # move to directory
echo ${BENCHMARK_DIRECTORY}
date

#
# Execute the code
#
echo "executing SeisSol: ${date}"
if [[ ${ARCHITECTURE} == *"knc" ]]
then
    echo "executin SeisSol on: $HOSTNAME-$MIC"
    ssh $HOSTNAME-$MIC "cd $BENCHMARK_DIRECTORY ; mpiexec -genv I_MPI_DAPL_PROVIDER_LIST ofa-v2-mlx4_0-1 -genv I_MPI_DEBUG 4 -genv I_MPI_PIN_DOMAIN socket -genv MIC_ENV_PREFIX MIC -genv MIC_OMP_NUM_THREADS 236 -genv KMP_AFFINITY proclist=[1-236],granularity=thread,explicit,verbose -genv MIC_USE_2MB_BUFFERS 64K -ppn ${RANKS_PER_NODE} -n ${NUMBER_OF_MPI_RANKS} ./${EXECUTABLE} parameters.par"
else
    poe ./${EXECUTABLE} parameters.par
fi

echo "that's it; let's get out of here!"
