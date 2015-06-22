#!/bin/bash
##
# @file
# This file is part of SeisSol.
#
# @author Fabio Gratl (f.gratl AT in.tum.de)
#
# @section LICENSE
# Copyright (c) 2013, SeisSol Group
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
# reads jobID from file
# monitors job and reacts on exit state
# renames output file for further processing

runOutput=$1
outputFile=$2
jobstate=INCOMPLETE
completed=" COMPLETED "
f=" FAILED "
ca=" CANCELLED "
nf=" NODE_FAIL "
pr=" PREEMPTED "
to=" TIMEOUT "

echo read from ${runOutput}
while read line
do
		jobID=$line
		echo "read        $jobID"
done < ${runOutput}

jobID=${jobID##S* }
echo "modified to $jobID"


#monitor job until state is completed
while [[ $jobstate != *$completed* ]]
do
	jobstate=$(sacct -aX | grep $jobID)
	echo -e $jobstate
	sleep 30 
	
	#exit if job stops without completion
	if [[ $jobstate == *$ca* ]] || [[ $jobstate == *$f* ]] || [[ $jobstate == *$nf* ]] || [[ $jobstate == *$pr* ]] || [[ $jobstate == *$to* ]]
	then
		echo "Job didn't finish normally"
		exit 1
	fi
done
echo -e "Job is done\n"

#renaming output
echo "renaming output"
#mv *$jobID* outputSlurm_${bench}_${numberOfBasisFunctions}${slurmScript}${ranks}${ID}
mv *$jobID* $outputFile
