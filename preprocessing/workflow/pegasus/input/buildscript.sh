#!/bin/bash
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

# Usage info
show_help() {
cat << EOF
Usage: ${0##*/} [-h] [-d COMPILE_MODE -e ENVIRONMENT -d COMPILE_MODE -c CODE_VERSION-v VECTOR_INSTRUCTION_SET -p PARALLELIZATION -s SCALASCA -m MESH_FORMAT -q NUMBER_OF_QUANTITIES -b NUMBER_OF_BASIS_FUNCTIONS]
Builds the code.
     -h                                       display this help and exit
     -e ENVIRONMENT                           script executed to setup bash-variables
     -d COMPILE_MODE                          release mode.
     -c CODE_VERSION                          code version used for execution.
     -v VECTOR_INSTRUCTION_SET                vector instruction set.
     -p PARALLELIZATION                       parallelization.
     -s SCALASCA                              Scalasca.
     -m INPUT_FORMAT                          mesh format (ascii, Netcdf).
     -q NUMBER_OF_QUANTITIES                  number of quantities.
     -b NUMBER_OF_BASIS_FUNCTIONS             number of basis functions used.
     -t NUMBER_OF_TEMPORAL_INTEGRATION_POINTS number of temporal integration points.
EOF
}

#
# parse command line arguments
#
ENVIRONMENT=NOT_SET
COMPILE_MODE=NOT_SET
CODE_VERSION=NOT_SET
VECTOR_INSTRUCTION_SET=NOT_SET
PARALLELIZATION=NOT_SET
SCALASCA=NOT_SET
MESH_FORMAT=NOT_SET
NUMBER_OF_QUANTITIES=NOT_SET
NUMBER_OF_BASIS_FUNCTIONS=NOT_SET
NUMBER_OF_TEMPORAL_INTEGRATION_POINTS=auto

OPTIND=1
while getopts "he:d:c:v:p:s:m:q:b:t:" opt; do
    case "$opt" in
        h)
            show_help
            exit 0
            ;;
        e) ENVIRONMENT=$OPTARG
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
        m) MESH_FORMAT=$OPTARG
            ;;
        q) NUMBER_OF_QUANTITIES=$OPTARG
            ;;
        b) NUMBER_OF_BASIS_FUNCTIONS=$OPTARG
            ;;
        t) NUMBER_OF_TEMPORAL_INTEGRATION_POINTS=$OPTARG
            ;;
        '?')
            show_help >&2
            exit 1
            ;;
    esac
done
shift "$((OPTIND-1))" # Shift off the options and optional --.

#
# Set up environment
#
echo "loading modules"
source /etc/profile
source /etc/profile.d/modules.sh

module load fortran/intel
module load ccomp/intel

if [ "${ENVIRONMENT}" != "NOT_SET" ]
then
  source ${ENVIRONMENT}
fi

# load custom toolchain for the Phis
# TODO: This is a fix for the broken environment
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

# list modules
module list
mpiicpc --version
mpiifort --version

TIME=$(date +"%y_%m_%d-%H_%M_%S")
MAIN_DIRECTORY=$(pwd)/trunk
TMP_BUILD_DIRECTORY=${TMPDIR}/build_${COMPILE_MODE}_${CODE_VERSION}_${VECTOR_INSTRUCTION_SET}_${PARALLELIZATION}_${SCALASCA}_${NUMBER_OF_QUANTITIES}_${NUMBER_OF_BASIS_FUNCTIONS}_${TIME}

#
# copy files to tmp dir
#
echo "copying trunk and moving to ${TMP_BUILD_DIRECTORY}"
mkdir ${TMP_BUILD_DIRECTORY}
cp -r ${MAIN_DIRECTORY}/* ${TMP_BUILD_DIRECTORY}
cd ${TMP_BUILD_DIRECTORY}

# setup netcdf format switch
if [ "${MESH_FORMAT}" == "Netcdf" ]
then
  NETCDF_SWITCH=yes
else
  NETCDF_SWITCH=no
fi

#
# build the code
#
echo "building the code"
date
if [ "${CODE_VERSION}" == "generatedKernels" ]
then
  # building for SuperMUC / MAC Cluster
  scons -f SConstruct_generatedKernels logLevel=warning logLevel0=info compileMode=${COMPILE_MODE} generatedKernels=yes vectorInstructionSet=${VECTOR_INSTRUCTION_SET} parallelization=${PARALLELIZATION} scalasca=${SCALASCA} numberOfBasisFunctions=${NUMBER_OF_BASIS_FUNCTIONS} cppCompiler=mpiicpc fortranCompiler=mpiifort useExecutionEnvironment=yes netcdf=${NETCDF_SWITCH} netcdfDir=${NETCDF4DIR} numberOfTemporalIntegrationPoints=${NUMBER_OF_TEMPORAL_INTEGRATION_POINTS}
elif [ "${CODE_VERSION}" == "classic" ]
then
  # build old classic SeisSol code
  scons -f SConstruct_generatedKernels logLevel=warning logLevel0=info compileMode=${COMPILE_MODE} generatedKernels=no vectorInstructionSet=${VECTOR_INSTRUCTION_SET} parallelization=mpi scalasca=${SCALASCA} numberOfBasisFunctions=${NUMBER_OF_BASIS_FUNCTIONS} cppCompiler=mpiicpc fortranCompiler=mpiifort useExecutionEnvironment=yes netcdf=${NETCDF_SWITCH} netcdfDir=${NETCDF4DIR} numberOfTemporalIntegrationPoints=${NUMBER_OF_TEMPORAL_INTEGRATION_POINTS}
else
  echo "unknown code version, aborting"
  exit 1
fi
date

#
# Copy executables back and remove temporary directory
#
echo "copying executables back.."
cp ${TMP_BUILD_DIRECTORY}/build/SeisSol* ${MAIN_DIRECTORY}/build
rm -r ${TMP_BUILD_DIRECTORY}

echo "Finished! at: " $(date)
