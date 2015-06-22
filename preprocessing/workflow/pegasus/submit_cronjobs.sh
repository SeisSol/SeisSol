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

#
# Read command line arguments
#
show_help() {
cat << EOF
Usage: ${0##*/} [-h] [-d CRONJOB_DIRECTORY]
Submits a Benchmark-job to Slurm.
     -h                           display this help and exit.
     -d CRONJOB_DIRECTORY         directory where the dax-files of the cronjobs are located.
     -s SUBMIT_DIRECTORY          directory where the executable workflows are generated to.
EOF
}

CRONJOB_DIRECTORY=NOT_SET
SUBMIT_DIRECTORY=NOT_SET

OPTIND=1
while getopts "hd:s:" opt; do
    case "$opt" in
        h)
            show_help
            exit 0
            ;;
        d) CRONJOB_DIRECTORY=$OPTARG
            ;;
        s) SUBMIT_DIRECTORY=$OPTARG
            ;;
        '?')
            show_help >&2
            exit 1
            ;;
    esac
done
shift "$((OPTIND-1))" # Shift off the options and optional --.

TIME=$(date +"%y_%m_%d-%H_%M_%S")

# Execute the workflows
for FILE in ${CRONJOB_DIRECTORY}/*;
do
  if [[ ${FILE} == *.dax* ]]
  then
    # TODO: Remove nocleanup, when our catalogues are save
    pegasus-plan --conf pegasus.conf --dax ${FILE} --dir ${SUBMIT_DIRECTORY} --sites local --output-site local --output-dir ${HOME}/${FILE}_${TIME} --submit --nocleanup
    # give the svn repo some space, TODO: we should use a single checkout
    sleep 10m
  fi
done
