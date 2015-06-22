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
Usage: ${0##*/} [-h] [-e ENVIRONMENT -p PREPROCESSOR -t TEMPLATE -o OUTPUT -r ORDER -f MESH_FORMAT -b MESH_BASE_NAME]
Submits a Benchmark-job to Slurm.
     -h              display this help and exit
     -e ENVIRONMENT               script executed to setup bash-variables
     -p PREPROCESSOR   path to executable of the used preprocessor.
     -t TEMPLATE       path of the template file.
     -o OUTPUT         path of the output file.
     -r ORDER          order of the ADER-DG methods (2, 3, 4, 5, 6)
     -f MESH_FORMAT    mesh format ('Gambit3D-Mixed', 'Gambit3D-fast', 'Netcdf').
     -b MESH_BASE_NAME base name of the mesh (without file ending).
EOF
}

#
# parse command line arguments
#
ENVIRONMENT=NOT_SET
PREPROCESSOR=NOT_SET
TEMPLATE=NOT_SET
OUTPUT=NOT_SET
ORDER=NOT_SET
MESH_FORMAT=NOT_SET
MESH_BASE_NAME=NOT_SET

OPTIND=1
while getopts "he:p:t:o:r:f:b:" opt; do
    case "$opt" in
        h)
            show_help
            exit 0
            ;;
        e) ENVIRONMENT=$OPTARG
            ;;
        p) PREPROCESSOR=$OPTARG
            ;;
        t) TEMPLATE=$OPTARG
            ;;
        o) OUTPUT=$OPTARG
            ;;
        r) ORDER=$OPTARG
            ;;
        f) MESH_FORMAT=$OPTARG
            ;;
        b) MESH_BASE_NAME=$OPTARG
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

echo ${PREPROCESSOR} -P ${TEMPLATE} -DORDER_TAG=${ORDER} -DMESH_FORMAT_TAG=\'${MESH_FORMAT}\' -DMESH_BASE_NAME_TAG=\'${MESH_BASE_NAME}\'> ${OUTPUT}
