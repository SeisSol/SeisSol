#!/bin/env python
##
# @file
# This file is part of 
#
# @author Leonhard Rannabauer (lrannabauer AT mytum.de)
# @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
#
# @section LICENSE
# Copyright (c) 2014-2015, SeisSol Group
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
# Generates a loadlever script from a generic example script

import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--output", required = True, help="" , dest="output", action="store")
parser.add_argument("--queue", required = True, help="" , dest="queue", action="store")
parser.add_argument("--job_name", required = True, help="" , dest="job_name", action="store")
parser.add_argument("--numberofnodes", required = True, help="" , dest="numberofnodes", action="store")
parser.add_argument("--time", required = True, help="" , dest="time", action="store")
parser.add_argument("--ranks_per_node", required = True, help="" , dest="ranks_per_node", action="store")
parser.add_argument("--compile_mode", required = True, help="" , dest="compile_mode", action="store")
parser.add_argument("--code_version", required = True, help="" , dest="code_version", action="store")
parser.add_argument("--architecture", required = True, help="" , dest="architecture", action="store")
parser.add_argument("--parallelization", required = True, help="" , dest="parallelization", action="store")
parser.add_argument("--communicationThread", required = True, help="" , dest="communicationThread", action="store")
parser.add_argument("--scalasca", required = True, help="" , dest="scalasca", action="store")
parser.add_argument("--number_of_quantities", required = True, help="" , dest="number_of_quantities", action="store")
parser.add_argument("--order", required = True, help="" , dest="order", action="store")
parser.add_argument("--benchmark_directory", required = True, help="" , dest="benchmark_directory", action="store")
parser.add_argument("--number_of_mpi_ranks", required = True, help="" , dest="number_of_mpi_ranks", action="store")
parser.add_argument("--threads_per_rank", required = True, help="" , dest="threads_per_rank", action="store")
parser.add_argument("--executable", required = True, help="" , dest="executable", action="store")
parser.add_argument("--scripts_directory", required = True, help="" , dest="scripts_directory", action="store")
parser.add_argument("--filename", required = True, help="" , dest='filename', action="store")
parser.add_argument("--genll", required = True, help="" , dest='genll', action="store")

args = parser.parse_args()

replacements = [    ("{output}",args.output),\
                    ("{queue}", args.queue),\
                    ("{job_name}",args.job_name),\
                    ("{numberofnodes}",args.numberofnodes),\
                    ("{time}",args.time),\
                    ("{ranks_per_node}",args.ranks_per_node),\
                    ("{compile_mode}",args.compile_mode),\
                    ("{code_version}",args.code_version),\
                    ("{architecture}",args.architecture),\
                    ("{parallelization}",args.parallelization),\
                    ("{communicationThread}",args.communicationThread),\
                    ("{scalasca}",args.scalasca),\
                    ("{number_of_quantities}",args.number_of_quantities),\
                    ("{order}",args.order),\
                    ("{benchmark_directory}",args.benchmark_directory),\
                    ("{number_of_mpi_ranks}",args.number_of_mpi_ranks),\
                    ("{ranks_per_node}",args.ranks_per_node),\
                    ("{threads_per_rank}",args.threads_per_rank),\
                    ("{executable}",args.executable),\
                    ("{scripts_directory}",args.scripts_directory) ]

if 'knc' in args.architecture:
  replacements = replacements + [ ('{mic}', 'mic0') ]
else:
  replacements = replacements + [ ('{mic}', '') ]
    
llfile = open(args.filename,'w+')
genfile= open(args.genll)
for line in genfile:
    for tag, replacement in replacements:
        if tag in line:
            line = line.replace(tag,replacement)
    llfile.write(line)
                
llfile.close()
genfile.close()

