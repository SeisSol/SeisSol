#!/bin/bash
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
# @section DESCRIPTION
# Builds SeisSol with several options in a temporary build directory.

# Usage info
show_help() {
cat << EOF
Usage: ${0##*/} [-h] [-d COMPILE_MODE -e ENVIRONMENT -d COMPILE_MODE -c CODE_VERSION -a ARCHITECTURE -p PARALLELIZATION -r COMMUNICATION_THREAD -s SCALASCA -m MESH_FORMAT -u OUTPUT_FORMAT -q NUMBER_OF_QUANTITIES -o CONVERGENCE_ORDER -n PROGRAM_NAME]
Builds the code.
     -h                                       display this help and exit
     -e ENVIRONMENT                           script executed to setup bash-variables
     -d COMPILE_MODE                          release mode.
     -c CODE_VERSION                          code version used for execution.
     -a ARCHITECTURE                          architecture.
     -p PARALLELIZATION                       parallelization.
     -r COMMUNICATION_THREAD                  communication thread.
     -s SCALASCA                              Scalasca.
     -m INPUT_FORMAT                          mesh format (ascii, Netcdf).
     -u OUTPUT_FORMAT                         output format (hdf5, default: none).
     -q NUMBER_OF_QUANTITIES                  number of quantities.
     -o ORDER                                 order of convergence of the used ADER-DG method.
     -t NUMBER_OF_TEMPORAL_INTEGRATION_POINTS number of temporal integration points.
     -n PROGRAM_NAME                          name of the executable.
EOF
}

#
# parse command line arguments
#
ENVIRONMENT=NOT_SET
COMPILE_MODE=NOT_SET
CODE_VERSION=NOT_SET
ARCHITECTURE=NOT_SET
PARALLELIZATION=NOT_SET
COMMUNICATION_THREAD=NOT_SET
SCALASCA=NOT_SET
OUTPUT_FORMAT=none
MESH_FORMAT=NOT_SET
NUMBER_OF_QUANTITIES=NOT_SET
ORDER=NOT_SET
NUMBER_OF_TEMPORAL_INTEGRATION_POINTS=auto
PROGRAM_NAME=none

OPTIND=1
while getopts "he:d:c:a:p:r:s:m:u:q:o:t:n:" opt; do
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
        a) ARCHITECTURE=$OPTARG
            ;;
        p) PARALLELIZATION=$OPTARG
            ;;
        r) COMMUNICATION_THREAD=$OPTARG
            ;;
        s) SCALASCA=$OPTARG
            ;;
        m) MESH_FORMAT=$OPTARG
            ;;
        u) OUTPUT_FORMAT=$OPTARG
            ;;
        q) NUMBER_OF_QUANTITIES=$OPTARG
            ;;
        o) ORDER=$OPTARG
            ;;
        t) NUMBER_OF_TEMPORAL_INTEGRATION_POINTS=$OPTARG
            ;;
        n) PROGRAM_NAME=$OPTARG
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

# take some guesses and print our environment
module list
echo $PATH
echo $LIBPATH
echo $LD_LIBRARY_PATH
mpiicpc --version
mpiifort --version
mpiCC --version
mpicxx --version
mpif90 --version

TIME=$(date +"%y_%m_%d-%H_%M_%S")
MAIN_DIRECTORY=$(pwd)
TMP_BUILD_DIRECTORY=${TMPDIR}/build_${PROGRAM_NAME}_${TIME}

#
# copy files to tmp dir
#
echo "copying trunk and moving to ${TMP_BUILD_DIRECTORY}"
mkdir ${TMP_BUILD_DIRECTORY}
cp -r ${MAIN_DIRECTORY}/* ${TMP_BUILD_DIRECTORY}
cd ${TMP_BUILD_DIRECTORY}

# setup hdf5 format switch
if [ "${OUTPUT_FORMAT}" == "hdf5" ]
then
  HDF5_SWITCH=yes
else
  HDF5_SWITCH=no
fi

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
if [ "${CODE_VERSION}" == "generated" ]
then
    # building generated kernels version
#    scons -f SConstruct_generatedKernels logLevel=warning logLevel0=info compileMode=${COMPILE_MODE} generatedKernels=yes arch=${ARCHITECTURE} parallelization=${PARALLELIZATION} commThread=${COMMUNICATION_THREAD} scalasca=${SCALASCA} order=${ORDER} cppCompiler=mpiicpc fortranCompiler=mpiifort useExecutionEnvironment=yes netcdf=${NETCDF_SWITCH} hdf5=${HDF5_SWITCH} zlibDir=${ZLIB_DIR} hdf5Dir=${HDF5_DIR} netcdfDir=${NETCDF4DIR} numberOfTemporalIntegrationPoints=${NUMBER_OF_TEMPORAL_INTEGRATION_POINTS} programName=${PROGRAM_NAME}
    scons logLevel=warning logLevel0=info compileMode=${COMPILE_MODE} generatedKernels=yes arch=${ARCHITECTURE} parallelization=${PARALLELIZATION} commThread=${COMMUNICATION_THREAD} scalasca=${SCALASCA} order=${ORDER} compiler=intel useExecutionEnvironment=yes netcdf=${NETCDF_SWITCH} hdf5=${HDF5_SWITCH} zlibDir=${ZLIB_DIR} hdf5Dir=${HDF5_DIR} netcdfDir=${NETCDF4DIR} numberOfTemporalIntegrationPoints=${NUMBER_OF_TEMPORAL_INTEGRATION_POINTS} programName=${PROGRAM_NAME}
elif [ "${CODE_VERSION}" == "classic" ]
then
  # build old classic SeisSol code
#  scons -f SConstruct_generatedKernels logLevel=warning logLevel0=info compileMode=${COMPILE_MODE} generatedKernels=no arch=${ARCHITECTURE} parallelization=mpi scalasca=${SCALASCA} order=${ORDER} cppCompiler=mpiicpc fortranCompiler=mpiifort useExecutionEnvironment=yes netcdf=${NETCDF_SWITCH} hdf5=${HDF5_SWITCH} zlibDir=${ZLIB_DIR} hdf5Dir=${HDF5_DIR} netcdfDir=${NETCDF4DIR} numberOfTemporalIntegrationPoints=${NUMBER_OF_TEMPORAL_INTEGRATION_POINTS} programName=${PROGRAM_NAME}
  scons logLevel=warning logLevel0=info compileMode=${COMPILE_MODE} generatedKernels=no arch=${ARCHITECTURE} parallelization=mpi scalasca=${SCALASCA} order=${ORDER} compiler=intel useExecutionEnvironment=yes netcdf=${NETCDF_SWITCH} hdf5=${HDF5_SWITCH} zlibDir=${ZLIB_DIR} hdf5Dir=${HDF5_DIR} netcdfDir=${NETCDF4DIR} numberOfTemporalIntegrationPoints=${NUMBER_OF_TEMPORAL_INTEGRATION_POINTS} programName=${PROGRAM_NAME}
else
  echo "unknown code version, aborting"
  exit 1
fi
date

#
# Copy executables back and remove temporary directory
#
echo "copying executables back.."
if [ "${PROGRAM_NAME}" == "none" ]
then
  cp ${TMP_BUILD_DIRECTORY}/build/SeisSol* ${MAIN_DIRECTORY}/build
else
  cp ${TMP_BUILD_DIRECTORY}/build/${PROGRAM_NAME} ${MAIN_DIRECTORY}/build
fi
rm -r ${TMP_BUILD_DIRECTORY}

echo "Finished! at: " $(date)
