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
Usage: ${0##*/} [-h] [-m MODE -c NUMBER_OF_CARDS_PER_HOST -b BENCHMARK_DIR]
Sets up the Phis and upload the benchmark directory to the local homes, if -m upload is used.
Downloads the receiver data from the Phis ${BENCHMARK_DIR}/output, if -m download is used. 
     -h                          display this help and exit
     -m MODE                     upload or download
     -c NUMBER_OF_CARDS_PER_HOST number of cards per host
     -b BENCHMARK_DIR            path to the benchmark directory, which will be uploaded/downloaded to/from all phis.
EOF
}

#
# parse command line arguments
#
MODE=NOT_SET
BENCHMARK_DIR=NOT_SET
NUMBER_OF_CARDS_PER_HOST=1

OPTIND=1
while getopts "hm:b:c:" opt; do
    case "$opt" in
        h)
            show_help
            exit 0
            ;;
        m) MODE=$OPTARG
            ;;
        b) BENCHMARK_DIR=$OPTARG
            ;;
        c) NUMBER_OF_CARDS_PER_HOST=$OPTARG
            ;;
        '?')
            show_help >&2
            exit 1
            ;;
    esac
done
shift "$((OPTIND-1))" # Shift off the options and optional --.

echo "running: $(date)"

#
# upload mode
#
if [ "${MODE}" == "upload" ]
then
  # ensure cards are booted
  srun uptime

  # create host files
  scontrol show hostnames ${SLURM_JOB_NODELIST} > ${SLURM_JOB_ID}.hosts
  for (( CARD=0; CARD<${NUMBER_OF_CARDS_PER_HOST}; CARD++ ))
  do
    scontrol show hostnames ${SLURM_JOB_NODELIST} > mic_${CARD}.${SLURM_JOB_ID}.hosts
    sed -e "s/$/-mic${CARD}/" -i mic_${CARD}.${SLURM_JOB_ID}.hosts
  done

  # create mic host file
  cat mic_*.${SLURM_JOB_ID}.hosts | sort > mic.${SLURM_JOB_ID}.hosts

  rm mic_*.${SLURM_JOB_ID}.hosts

  # print mic hosts
  echo "Phi cards:"
  cat mic.${SLURM_JOB_ID}.hosts

  for CARD in `cat mic.${SLURM_JOB_ID}.hosts`
  do
    ssh ${CARD} -o StrictHostKeyChecking=no 'rm -rf benchmark'
    ssh ${CARD} 'mkdir benchmark'
    scp -r ${BENCHMARK_DIR}/* ${CARD}:~/benchmark &
    ssh ${CARD} 'echo ~/benchmark/Maple/ > ~/benchmark/DGPATH' &
  done

  echo "Waiting for all copies to Phis: "$(date)

  # wait for all child processes
  wait
  echo "Done: "$(date)
#
# download mode
#
elif [ "${MODE}" == "download" ]
then
  for CARD in `cat mic.${SLURM_JOB_ID}.hosts`
  do
    scp -r ${CARD}:~/benchmark/output/* ${BENCHMARK_DIR}/output/ &
  done
  
  echo "Waiting for all downloads from Phis: $(date)"

  # wait for all child processes
  wait
  echo "Done: $(date)"
fi
