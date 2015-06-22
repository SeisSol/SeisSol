#!/bin/bash
##
# @file
# This file is part of SeisSol.
#
# @author Fabio Gratl (f.gratl AT in.tum.de)
# @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
#
# @section LICENSE
# Copyright (c) 2013-2014, SeisSol Group
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
Usage: ${0##*/} [-h] [-i INPUT_DIRECTORY -r RC_FILE -s]
Generates a replica catalogue based on the given setting.
     -h                 display this help and exit
     -i INPUT_DIRECTORY location of the input directory.
     -r RC_FILE         output file of the replica catalogue
     -s                 usingle temporal integration point in the dynamic rupture boundary conditions if set (optional).
EOF
}

#
# parse command line arguments
#
INPUT_DIRECTORY=NOT_SET
RC_FILE=NOT_SET
TEMPORAL_INTEGRATION_POINT_SWITCH=auto_temporal_int_points

OPTIND=1
while getopts "hi:r:s" opt; do
    case "$opt" in
        h)
            show_help
            exit 0
            ;;
        i) if [ $(echo $OPTARG | head -c 1) != "/" ]
           then
             # extend relative path with current path
             INPUT_DIRECTORY="$(pwd)/${OPTARG}"
           else
             # use absolute path
             INPUT_DIRECTORY=$OPTARG
           fi
            ;;
        r) RC_FILE=$OPTARG
            ;;
        s) TEMPORAL_INTEGRATION_POINT_SWITCH=single_temporal_int_point
            ;;
        '?')
            show_help >&2
            exit 1
            ;;
    esac
done
shift "$((OPTIND-1))" # Shift off the options and optional --.


echo "Input directory: ${INPUT_DIRECTORY}"
echo "Output file:  ${RC_FILE}"

echo -e "# This is the replica catalog. It lists information about each of the
# input files used by the workflow.

# The format is:
# LFN     PFN    pool="SITE"

env_vars.sh              file://${INPUT_DIRECTORY}/env_vars.sh    							pool=\"local\"
trunk.tar.gz             file://${INPUT_DIRECTORY}/trunk.tar.gz    							pool=\"local\"
mesh_TetraElastic.tar.gz file://${INPUT_DIRECTORY}/TetraElasticVerification_parameters/mesh_TetraElastic.tar.gz	pool=\"local\"
mesh_TetraLOH4.tar.gz    file://${INPUT_DIRECTORY}/TetraLOH4_parameters/mesh_TetraLOH4.tar.gz			pool=\"local\"
mesh_LOH1Scaling.tar.gz  file://${INPUT_DIRECTORY}/LOH1Scaling_parameters/mesh_LOH1Scaling.tar.gz			pool=\"local\"
mesh_Merapi.tar.gz       file://${INPUT_DIRECTORY}/Merapi_parameters/mesh_Merapi.tar.gz				pool=\"local\"
mesh_tpv16.tar.gz        file://${INPUT_DIRECTORY}/tpv16_parameters/mesh_tpv16.tar.gz				pool=\"local\"
mesh_landers.tar.gz      file://${INPUT_DIRECTORY}/landers/input/mesh_landers.tar.gz				pool=\"local\"

data_output.f90			file://${INPUT_DIRECTORY}/data_output_NoOutput.f90	pool=\"local\"

gnuplot_cube3_stat		file://${INPUT_DIRECTORY}/gnuplot_cube3_stat		pool=\"local\"
gnuplot_TetraElastic	file://${INPUT_DIRECTORY}/gnuplot_TetraElastic		pool=\"local\"

submitBenchmark			file://${INPUT_DIRECTORY}/submitBenchmark.slurm		pool=\"local\"
runTetraElastic    		file://${INPUT_DIRECTORY}/runTetraElastic.slurm	pool=\"local\"
runTetraLOH4			file://${INPUT_DIRECTORY}/runTetraLOH4.slurm	pool=\"local\"
runLOH1Scaling 			file://${INPUT_DIRECTORY}/runLOH1Scaling.slurm	pool=\"local\"
runMerapi			file://${INPUT_DIRECTORY}/runMerapi.slurm		pool=\"local\"
runtpv16			file://${INPUT_DIRECTORY}/runTpv16.slurm		pool=\"local\"
submitBenchmark  file://${INPUT_DIRECTORY}/submitBenchmark.slurm		pool=\"local\"
phi_communication  file://${INPUT_DIRECTORY}/phi_communication.sh		pool=\"local\"

landers_parameters.template        file://${INPUT_DIRECTORY}/landers/input/parameters.template pool=\"local\"
TetraElastic_parameters.template   file://${INPUT_DIRECTORY}/TetraElasticVerification_parameters/parameters.template pool=\"local\"
TetraLOH4_parameters.template      file://${INPUT_DIRECTORY}/TetraLOH4_parameters/parameters.template pool=\"local\"
LOH1Scaling_parameters.template    file://${INPUT_DIRECTORY}/LOH1Scaling_parameters/parameters.template pool=\"local\"
Merapi_parameters.template         file://${INPUT_DIRECTORY}/Merapi_parameters/parameters.template pool=\"local\"
tpv16_parameters.template          file://${INPUT_DIRECTORY}/tpv16_parameters/parameters.template pool=\"local\"

TetraElastic_referenceFile_4	file://${INPUT_DIRECTORY}/TetraElasticVerification_reference/Tetra_Elastic_generatedKernels_4	pool=\"local\"
TetraElastic_referenceFile_10	file://${INPUT_DIRECTORY}/TetraElasticVerification_reference/Tetra_Elastic_generatedKernels_10	pool=\"local\"
TetraElastic_referenceFile_20	file://${INPUT_DIRECTORY}/TetraElasticVerification_reference/Tetra_Elastic_generatedKernels_20	pool=\"local\"
TetraElastic_referenceFile_35	file://${INPUT_DIRECTORY}/TetraElasticVerification_reference/Tetra_Elastic_generatedKernels_35	pool=\"local\"
TetraElastic_referenceFile_56	file://${INPUT_DIRECTORY}/TetraElasticVerification_reference/Tetra_Elastic_generatedKernels_56	pool=\"local\"

TetraLOH4_referenceFiles_4	file://${INPUT_DIRECTORY}/TetraLOH4_reference/TetraLOH4_referenceFiles_4.tar.gz 	pool=\"local\"
TetraLOH4_referenceFiles_10	file://${INPUT_DIRECTORY}/TetraLOH4_reference/TetraLOH4_referenceFiles_10.tar.gz 	pool=\"local\"
TetraLOH4_referenceFiles_20	file://${INPUT_DIRECTORY}/TetraLOH4_reference/TetraLOH4_referenceFiles_20.tar.gz 	pool=\"local\"
TetraLOH4_referenceFiles_35	file://${INPUT_DIRECTORY}/TetraLOH4_reference/TetraLOH4_referenceFiles_35.tar.gz 	pool=\"local\"
TetraLOH4_referenceFiles_56	file://${INPUT_DIRECTORY}/TetraLOH4_reference/TetraLOH4_referenceFiles_56.tar.gz 	pool=\"local\"

LOH1Scaling_referenceFiles_4	file://${INPUT_DIRECTORY}/LOH1Scaling_reference/LOH1Scaling_referenceFiles_4.tar.gz		pool=\"local\"
LOH1Scaling_referenceFiles_10	file://${INPUT_DIRECTORY}/LOH1Scaling_reference/LOH1Scaling_referenceFiles_10.tar.gz 	pool=\"local\"
LOH1Scaling_referenceFiles_20	file://${INPUT_DIRECTORY}/LOH1Scaling_reference/LOH1Scaling_referenceFiles_20.tar.gz 	pool=\"local\"
LOH1Scaling_referenceFiles_35	file://${INPUT_DIRECTORY}/LOH1Scaling_reference/LOH1Scaling_referenceFiles_35.tar.gz 	pool=\"local\"
LOH1Scaling_referenceFiles_56	file://${INPUT_DIRECTORY}/LOH1Scaling_reference/LOH1Scaling_referenceFiles_56.tar.gz 	pool=\"local\"

tpv16_referenceFiles_4		file://${INPUT_DIRECTORY}/tpv16_reference/${TEMPORAL_INTEGRATION_POINT_SWITCH}/tpv16_referenceFiles_4.tar.gz	pool=\"local\"
tpv16_referenceFiles_10		file://${INPUT_DIRECTORY}/tpv16_reference/${TEMPORAL_INTEGRATION_POINT_SWITCH}/tpv16_referenceFiles_10.tar.gz 	pool=\"local\"
tpv16_referenceFiles_20		file://${INPUT_DIRECTORY}/tpv16_reference/${TEMPORAL_INTEGRATION_POINT_SWITCH}/tpv16_referenceFiles_20.tar.gz 	pool=\"local\"
tpv16_referenceFiles_35		file://${INPUT_DIRECTORY}/tpv16_reference/${TEMPORAL_INTEGRATION_POINT_SWITCH}/tpv16_referenceFiles_35.tar.gz 	pool=\"local\"
tpv16_referenceFiles_56		file://${INPUT_DIRECTORY}/tpv16_reference/${TEMPORAL_INTEGRATION_POINT_SWITCH}/tpv16_referenceFiles_56.tar.gz 	pool=\"local\"

landers_referenceFiles_4		file://${INPUT_DIRECTORY}/landers/references/${TEMPORAL_INTEGRATION_POINT_SWITCH}/landers_referenceFiles_4.tar.gz	pool=\"local\"
landers_referenceFiles_10		file://${INPUT_DIRECTORY}/landers/references/${TEMPORAL_INTEGRATION_POINT_SWITCH}/landers_referenceFiles_10.tar.gz 	pool=\"local\"
landers_referenceFiles_20		file://${INPUT_DIRECTORY}/landers/references/${TEMPORAL_INTEGRATION_POINT_SWITCH}/landers_referenceFiles_20.tar.gz 	pool=\"local\"
landers_referenceFiles_35		file://${INPUT_DIRECTORY}/landers/references/${TEMPORAL_INTEGRATION_POINT_SWITCH}/landers_referenceFiles_35.tar.gz 	pool=\"local\"
landers_referenceFiles_56		file://${INPUT_DIRECTORY}/landers/references/${TEMPORAL_INTEGRATION_POINT_SWITCH}/landers_referenceFiles_56.tar.gz 	pool=\"local\"
" > ${RC_FILE}
