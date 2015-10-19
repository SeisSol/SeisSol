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
  echo "you need to specify 4 argumenrs: WORKFLOW_INPUT_DIR"
  exit
fi

# read arguments
if [[ $# -eq 1 ]]
then
  WORKFLOW_INPUT_DIR=$1
fi

CODE_ARCHIVE_OLD=${WORKFLOW_INPUT_DIR}/code.tar.gz
CODE_ARCHIVE_NEW=${WORKFLOW_INPUT_DIR}/code_new.tar.gz
HASH_OLD=00000000000000000000000000000000
HASH_NEW=00000000000000000000000000000000
#return 1 means no new code, nothing to prepare
RETURN=1

# generate md5 hase of current code
if [[ -a ${CODE_ARCHIVE_OLD} ]]
then
  HASH_OLD=`md5sum ${CODE_ARCHIVE_OLD} | awk '{print $1}'`
fi

# clone current github
SAVED_DIR=`pwd`
cd /tmp/
git clone --recursive https://github.com/SeisSol/SeisSol.git
cd SeisSol
tar -czf $CODE_ARCHIVE_NEW *
cd ..
rm -rf SeisSol
cd $SAVE_DIR

# create hash of new code
HASH_NEW=`md5sum ${CODE_ARCHIVE_NEW} | awk '{print $1}'`

# compare and take action
if [[ ${HASH_NEW} != ${HASH_OLD} ]]
then
  mv ${CODE_ARCHIVE_NEW} ${CODE_ARCHIVE_OLD}
  RETURN=0
else
  rm ${CODE_ARCHIVE_NEW}
fi

exit ${RETURN}
