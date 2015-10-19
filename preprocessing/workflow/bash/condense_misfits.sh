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
################################################################################
# Alexander Heinecke (Intel Corp.)
################################################################################

# check if we got 1 arguments
if [[ $# -ne 1 ]]
then
  echo "you need to specify 4 argumenrs: WORKFLOW_OUTPUT_DIR"
  exit
fi

# read arguments
if [[ $# -eq 1 ]]
then
  WORKFLOW_OUTPUT_DIR=$1
fi

CONDENSED_MISFITS=${WORKFLOW_OUTPUT_DIR}/maxmisfits
rm -rf ${CONDENSED_MISFITS}
CUR_DATE=`date +"%Y-%m-%d"`
SUMMARY=${WORKFLOW_OUTPUT_DIR}/../${CUR_DATE}_maxmisfits.txt
touch ${SUMMARY}

#@TODO we might need some special handling for tpv16 since p_n is broken
for MISFITS in `find ${WORKFLOW_OUTPUT_DIR}/ -name "misfits.csv" | sort`
do
  MAX_MISFIT=`cat ${MISFITS} | awk -F"," '{print $5}' | sort -g | tail -n 1`
  cat ${MISFITS} | grep ${MAX_MISFIT} | awk -F"," '{print $2 ":\n" $3 " " $4 " " $5}' >> ${CONDENSED_MISFITS}
  echo "" >> ${CONDENSED_MISFITS}
done

cat ${CONDENSED_MISFITS} >> ${SUMMARY}
