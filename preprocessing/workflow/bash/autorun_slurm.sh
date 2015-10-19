#!/bin/bash
#
# Copyright (c) 2015, Intel Corporation
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
#     * Redistributions of source code must retain the above copyright notice,
#       this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Intel Corporation nor the names of its contributors
#       may be used to endorse or promote products derived from this software
#       without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# check if we got 4 arguments
if [[ $# -ne 4 ]]
then
  echo "you need to specify 4 argumenrs: WORKFLOW_SCRIPT WORKFLOW_INPUT_DIR SCRIPTS_DIR WORKFLOW_OUTPUT_DIR"
  exit
fi

# read arguments
if [[ $# -eq 4 ]]
then
  WORKFLOW_SCRIPT=$1
  WORKFLOW_INPUT_DIR=$2
  SCRIPTS_DIR=$3
  WORKFLOW_OUTPUT_DIR=$4
fi

# remove output dir
rm -rf ${WORKFLOW_OUTPUT_DIR}

#setup everything
bash ${WORKFLOW_SCRIPT} -s -i ${WORKFLOW_INPUT_DIR} -e ${SCRIPTS_DIR} -w ${WORKFLOW_OUTPUT_DIR}

# submitting jobs to SLURM
bash ${WORKFLOW_SCRIPT} -b -i ${WORKFLOW_INPUT_DIR} -e ${SCRIPTS_DIR} -w ${WORKFLOW_OUTPUT_DIR} | grep "Submitted" | awk '{print $4}' > ${WORKFLOW_OUTPUT_DIR}/jobs

# querying for completion of SLURM jobs
SLURM_WORKFLOW_JOBS=`cat ${WORKFLOW_OUTPUT_DIR}/jobs | wc | awk '{print $1}'`
SLURM_WORKFLOW_DONE=`cat ${WORKFLOW_OUTPUT_DIR}/jobs | wc | awk '{print $1}'`
while [[  ${SLURM_WORKFLOW_DONE} -gt 0 ]] 
do
  SLURM_WORKFLOW_DONE=0
  # loop over submitted jobs
  for i in `cat ${WORKFLOW_OUTPUT_DIR}/jobs`
  do 
    # check job
    JOB_STATUS=`squeue -j $i 2> /dev/null | wc | awk '{print $1}'`
    if [[ ${JOB_STATUS} -eq 2 ]]
    then 
      SLURM_WORKFLOW_DONE=$((SLURM_WORKFLOW_DONE + 1))
    fi
  done
  # print status
  CUR_TIMESTAMP=`date`
  echo "${CUR_TIMESTAMP} - running ${WORKFLOW_SCRIPT}, still to go: ${SLURM_WORKFLOW_DONE} of ${SLURM_WORKFLOW_JOBS}"
  # let's sleep 60s to avoid DOS attack of SLURM daemon
  if [[ ${SLURM_WORKFLOW_DONE} -ne 0 ]]
  then 
    sleep 60
  fi
done

#analyse everything
bash ${WORKFLOW_SCRIPT} -a -i ${WORKFLOW_INPUT_DIR} -e ${SCRIPTS_DIR} -w ${WORKFLOW_OUTPUT_DIR}
