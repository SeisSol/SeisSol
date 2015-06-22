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
# @section DESCRIPTION
# Removes ranks from SeisSol receiver naming scheme.

# Usage info
show_help() {
cat << EOF
Usage: ${0##*/} [-h] [-m MODE -i INPUT_DIRECTORY -o OUTPUT_DIRECTORY]
Removes the rank information in the receiver files of the input directory and stores them in the output directory
     -h                  display this help and exit
     -m MODE             mode: "copy" or "move"
     -i INPUT_DIRECTORY  input directory with full receiver names
     -o OUTPUT_DIRECTORY output directory where the receivers will be copied to with shortened names.
EOF
}

#
# parse command line arguments
#
MODE=copy
INPUT_DIRECTORY=NOT_SET
OUTPUT_DIRECTORY=NOT_SET

OPTIND=1
while getopts "hm:i:o:" opt; do
    case "$opt" in
        h)
            show_help
            exit 0
            ;;
        m) MODE=$OPTARG
            ;;
        i) INPUT_DIRECTORY=$OPTARG
            ;;
        o) OUTPUT_DIRECTORY=$OPTARG
            ;;
        '?')
            show_help >&2
            exit 1
            ;;
    esac
done
shift "$((OPTIND-1))" # Shift off the options and optional --.

# move to input directory
cd $INPUT_DIRECTORY

# iterate over receivers
for RECEIVER in *receiver*.dat;
do
  if [[ $RECEIVER =~ .+-[0-9]+-[0-9]+\.dat ]]
  then
    # MPI rank is attached, remove
    SHORTENED=${RECEIVER%-*.dat}
  else
    # just remove the ending
    SHORTENED=${RECEIVER%.dat}
  fi

  echo "renaming ${RECEIVER} to ${SHORTENED}"
  if [[ $MODE == "copy" ]]
  then
    cp ${RECEIVER} ${OUTPUT_DIRECTORY}/${SHORTENED}
  elif [[ $MODE == "move" ]]
  then
    mv ${RECEIVER} ${OUTPUT_DIRECTORY}/${SHORTENED}
  fi
done
