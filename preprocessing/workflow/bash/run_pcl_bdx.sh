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

# setup log
CUR_DATE_LOG=`date +"%Y-%m-%d"`
CRON_LOG=/scratch/aheineck/${CUR_DATE_LOG}_seissol_cornjob.log
touch ${CRON_LOG}

# source env
source /nfs_home/aheineck/.bashrc

# set some exports
export TMPDIR=/scratch/aheineck
export NETCDF4DIR=/swtools/netcdf/netcdf-4.3.0/
export HDF5_DIR=/swtools/hdf5/hdf5-1.8.11/
export ZLIB_DIR=/lib64

# dump env
env > ${CRON_LOG}

# make sure we are in the right directory
cd /nfs_home/aheineck/Projects/SeisSol_workspace/SeisSol/preprocessing/workflow/bash
pwd >> ${CRON_LOG}

# re-build workflows
cd configurations/pcl_bdx/
rm *.sh
python ./../../generate.py --workflows_xml ./workflows.xml >> ${CRON_LOG}
cd ./../../

# prepare code
bash ./prepare_code.sh /scratch/aheineck/input_workflow >> ${CRON_LOG}

# run and check dynamic rupture
bash ./autorun_slurm.sh /nfs_home/aheineck/Projects/SeisSol_workspace/SeisSol/preprocessing/workflow/bash/configurations/pcl_bdx/dynamic_rupture_generated_hybrid.sh /scratch/aheineck/input_workflow /nfs_home/aheineck/Projects/SeisSol_workspace/SeisSol/preprocessing/workflow/bash/scripts/ /scratch/aheineck/current_workflow_dr >> ${CRON_LOG}
bash ./condense_misfits.sh /scratch/aheineck/current_workflow_dr/ >> ${CRON_LOG}

# run and check wave propagation
bash ./autorun_slurm.sh /nfs_home/aheineck/Projects/SeisSol_workspace/SeisSol/preprocessing/workflow/bash/configurations/pcl_bdx/wave_propagation_generated_hybrid.sh /scratch/aheineck/input_workflow /nfs_home/aheineck/Projects/SeisSol_workspace/SeisSol/preprocessing/workflow/bash/scripts/ /scratch/aheineck/current_workflow_wp >> ${CRON_LOG}
bash ./condense_misfits.sh /scratch/aheineck/current_workflow_wp/ >> ${CRON_LOG}

unset TMPDIR
unset NETCDF4DIR
unset HDF5_DIR
unset ZLIB_DIR
