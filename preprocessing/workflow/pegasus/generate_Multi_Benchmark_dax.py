#!/usr/bin/env python
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

import os
import subprocess
import time
import xmltodict
import argparse
# Import the Python DAX library
from Pegasus.DAX3 import *


m_basisFunctionsToOrder = {  '4': '2',
                            '10': '3',
                            '20': '4',
                            '35': '5',
                            '56': '6' }

l_parser    = argparse.ArgumentParser( description='Setup up a multiple verification benchmarks' )
l_parser.add_argument( '--xml',
                       dest="xml_file",
                       required=True,
                       help="XML configuration for the benchmark runs.",
                       metavar="XML_FILE" )
l_parser.add_argument( '--dax',
                       dest="dax",
                       required=True,
                       help="Name of the dax and the corresponding dax file.",
                       metavar="DAX" )
l_parser.add_argument( '--svn',
                       dest='svn',
                       help="Use svn to checkout the repository.",
                       action='store_true' )

l_arguments = vars(l_parser.parse_args())


# parse the XML-setup
l_setup = xmltodict.parse( open( l_arguments['xml_file'] ) )['setup']

# create a abstract dag
print "creating ADAG..."
dax = ADAG(l_arguments['dax'])


if l_arguments['svn'] == False:
  print "using trunk.tar.gz from input folder"
  #Add files for unpackTrunkJob
  trunkPacked = File("trunk.tar.gz")
  trunk = File("trunk")
  errorFileTar = File("tarTrunk.err")
  outputFileTar = File("tarTrunk.out")
  # Copy packed Trunk to site and unpack it there
  print "Adding unpackTrunk job ..."
  print 'unpacked folder must be named "trunk"'
  unpackTrunkJob = Job(name="tar")
  unpackTrunkJob.addArguments("-x","-z","-v", "-f", trunkPacked)
  unpackTrunkJob.setStdout(outputFileTar)
  unpackTrunkJob.setStderr(errorFileTar)
  unpackTrunkJob.uses(trunkPacked, link=Link.INPUT)
  unpackTrunkJob.uses(trunk, link=Link.OUTPUT, transfer=False, register=False)
  dax.addJob(unpackTrunkJob)

  #Remove packed Trunk
  print "Adding removeTrunk job ..."
  removeTrunkJob = Job(name="rm")
  removeTrunkJob.addArguments(trunkPacked)
  dax.addJob(removeTrunkJob)
  #Trunk should be removed only after its unpacked
  dax.depends(parent=unpackTrunkJob, child=removeTrunkJob)
  
else:
  print "using svn to download SeisSol trunk"
  #Add files for git job
  trunk = File("trunk")
  errorFileSvn = File("svn.err")
  outputFileSvn = File("svn.out")
  #using git svn to download trunk 
  #job has same name as tar job because of dependencies
  unpackTrunkJob = Job(name="svn")
  unpackTrunkJob.addArguments("co", "https://svn.geophysik.uni-muenchen.de/svn/seissol/trunk")
  unpackTrunkJob.setStdout(outputFileSvn)
  unpackTrunkJob.setStderr(errorFileSvn)
  dax.addJob(unpackTrunkJob)

# dictionary for saving build jobs
l_buildJobs = {}

#build an run SeisSol for various number of basis functions parallel
for i in  l_setup['global']['@number_of_basis_functions'].split():
  
  # build SeisSol with script
  outputFileBuild_i = File("build_%s.out" %(i))
  errorFileBuild_i = File("build_%s.err" %(i))
  env_vars = File("env_vars.sh")
  
  l_buildArguments =   " -d " + l_setup['global']['@compile_mode']+\
                       " -c " + l_setup['global']['@code_version']+\
                       " -v " + l_setup['global']['@vector_instruction_set']+\
                       " -p " + l_setup['global']['@parallelization']+\
                       " -s " + l_setup['global']['@scalasca']+\
                       " -m " + l_setup['global']['@mesh_format']+\
                       " -q " + l_setup['global']['@number_of_quantities']+\
                       " -b " + str(i)
  if( '@number_of_temporal_integration_points' in l_setup['global'] ):
    l_buildArguments = l_buildArguments + " -t " + l_setup['global']['@number_of_temporal_integration_points']

  print "Adding buildscript.sh job for numberOfBasisFunctions=%s ..." %(i)
  l_buildJobs[i] = Job(name="buildscript")
  l_buildJobs[i].addArguments( "-e", env_vars, l_buildArguments )
  l_buildJobs[i].uses(env_vars, link=Link.INPUT)
  l_buildJobs[i].setStdout(outputFileBuild_i)
  l_buildJobs[i].setStderr(errorFileBuild_i)
  dax.addJob(l_buildJobs[i])

  #build needs unpacked trunk
  dax.depends(parent=unpackTrunkJob, child=l_buildJobs[i] )

#copy, unpack and remove mesh for every Benchmark
for benchmark in l_setup['benchmarks']:
  #files for unpackMeshJob
  meshPacked_bench = File("mesh_%s.tar.gz" %(benchmark))
  errorFileTar2_bench = File("tarMesh_%s.err" %(benchmark))
  outputFileTar2_bench = File("tarMesh_%s.out" %(benchmark))

  #copy packed mesh folder to site and unpack it there
  print "Adding unpackMesh job for %s ..." %(benchmark)
  print 'unpacked file must be named "mesh_%s"' %(benchmark)
  unpackMeshJob_bench = Job(name="tar")
  unpackMeshJob_bench.addArguments("-x","-z","-v", "-f", meshPacked_bench)
  unpackMeshJob_bench.setStdout(outputFileTar2_bench)
  unpackMeshJob_bench.setStderr(errorFileTar2_bench)
  unpackMeshJob_bench.uses(meshPacked_bench, link=Link.INPUT)
  dax.addJob(unpackMeshJob_bench)

  #Remove packed Mesh
  print "Adding removeMesh job for %s ..." %(benchmark)
  removeMeshJob_bench= Job(name="rm")
  removeMeshJob_bench.addArguments(meshPacked_bench)
  dax.addJob(removeMeshJob_bench)
  #Mesh should be removed only after its unpacked
  dax.depends(parent=unpackMeshJob_bench, child=removeMeshJob_bench)
  
  #counter for identifying the job in buildIDs list
  n=0
  
  l_benchmarkSettings = l_setup['benchmarks'][benchmark]

  for i in  l_setup['global']['@number_of_basis_functions'].split():
    #create Benchmark dir and copy mesh Files to it using script
    outputCopyMesh_i_benchmark = File("copyMeshToBench_%s_%s.out" %(benchmark, i))
    errorCopyMesh_i_benchmark = File("copyMeshToBench_%s_%s.err" %(benchmark, i))
    print "Adding copyMeshToBench.sh job for %s and numberOfBasisFunctions=%s ..." %(benchmark, i)
    copyMeshJob_i_benchmark = Job(name="copyMeshToBench")
    copyMeshJob_i_benchmark.addArguments("%s" %(i), "%s" %(benchmark))
    copyMeshJob_i_benchmark.setStdout(outputCopyMesh_i_benchmark)
    copyMeshJob_i_benchmark.setStderr(errorCopyMesh_i_benchmark)
    dax.addJob(copyMeshJob_i_benchmark)

    #copyMesh needs unpacked Mesh
    dax.depends(parent=unpackMeshJob_bench, child=copyMeshJob_i_benchmark)

    # generate parameter file
    l_generateParameterFile = Job(name="generate_parameter_file")
    l_generateParameterFile.uses( File("%s_parameters.template"%(benchmark)), link=Link.INPUT )
    l_generateParameterFile.addArguments( " -p " + Job(name="cpp").name+\
                                          " -t " + File( "%s_parameters.template"%(benchmark)).name+\
                                          " -o " + "Benchmark_%s_%s/parameters.par"%(benchmark, i)+\
                                          " -r " + m_basisFunctionsToOrder[i]+\
                                          " -f " + l_setup['global']['@mesh_format']+\
                                          " -b " + l_benchmarkSettings['@mesh_base_name'],
                                          "-e", env_vars )
    dax.addJob(l_generateParameterFile)
    l_generateParameterFile.uses(env_vars, link=Link.INPUT)
    dax.depends( parent=copyMeshJob_i_benchmark, child=l_generateParameterFile )

    #Adding new files for runJob
    outputFileRun_i_benchmark = File("run_%s_%s.out" %(benchmark, i))
    errorFileRun_i_benchmark = File("run_%s_%s.err" %(benchmark, i))

    # set up slurm job
    if "snb" in l_benchmarkSettings['@queue']:
      submitBenchmark = "--constraint=turbo_off "
    else:
      submitBenchmark = ""
    submitBenchmark =  submitBenchmark+\
                       "-o " + benchmark + ".%j.%N.out"+\
                       " -J " + l_arguments['dax'] + "_" + benchmark + "_" + str(i) + "_" + datetime.datetime.now().isoformat()+\
                       " --ntasks=" + l_benchmarkSettings['@number_of_nodes']+\
                       " --partition=" + l_benchmarkSettings['@queue']+\
                       " --time=" +  l_benchmarkSettings['@maximum_runtime']+\
                       " submitBenchmark"+\
                       " -d " + l_setup['global']['@compile_mode']+\
                       " -c " + l_setup['global']['@code_version']+\
                       " -v " + l_setup['global']['@vector_instruction_set']+\
                       " -p " + l_setup['global']['@parallelization']+\
                       " -s " + l_setup['global']['@scalasca']+\
                       " -q " + l_setup['global']['@number_of_quantities']+\
                       " -b " + str(i)+\
                       " -n " + benchmark+\
                       " -m " + l_benchmarkSettings['@number_of_mpi_ranks']+\
                       " -r " + l_benchmarkSettings['@ranks_per_node']+\
                       " -t " + l_benchmarkSettings['@threads_per_rank']
 
    #submit slurm run script with sbatch
    print "Adding sbatch run%s.slurm job for numberOfBasisFunctions=%s ..." %(benchmark, i)
    runJob_i_benchmark = Job(name="sbatch")
    runJob_i_benchmark.addArguments(submitBenchmark)
    runJob_i_benchmark.setStdout(outputFileRun_i_benchmark)
    runJob_i_benchmark.setStderr(errorFileRun_i_benchmark)
    runJob_i_benchmark.uses("submitBenchmark", link=Link.INPUT)
    if ( l_setup['global']['@vector_instruction_set'] == 'mic' ):
      runJob_i_benchmark.uses("phi_communication", link=Link.INPUT)
    dax.addJob(runJob_i_benchmark)  
    dax.depends(parent=l_generateParameterFile, child=runJob_i_benchmark)
    dax.depends(parent=l_buildJobs[i],              child=runJob_i_benchmark)
    dax.depends(parent=copyMeshJob_i_benchmark, child=runJob_i_benchmark)
    
    #Adding new files for monitor Job
    outputFileMonitor_i_benchmark = File("monitor_%s_%s.out" %(benchmark, i))
    errorFileMonitor_i_benchmark = File("monitor_%s_%s.err" %(benchmark, i))
    actualOutput_i_benchmark = File("outputSlurm_%s_%s" %(benchmark, i))

    #retrieving job id from run_i_benchmark.out
    #monitor slurm state and react on exit state
    #rename output file of slurm job to outputSlurm_benchmark_i
    monitorJob_i_benchmark =Job(name="monitorSlurmState")
    monitorJob_i_benchmark.addArguments(outputFileRun_i_benchmark, actualOutput_i_benchmark)
    monitorJob_i_benchmark.setStdout(outputFileMonitor_i_benchmark)
    monitorJob_i_benchmark.setStderr(errorFileMonitor_i_benchmark)
    dax.addJob(monitorJob_i_benchmark)

    #monitor needs run job to monitor s.th.
    dax.depends(parent=runJob_i_benchmark, child=monitorJob_i_benchmark)
   
    #if cases for different postprocessing of Benchmark output
    if benchmark == "TetraElastic":
      #Add new files for compareJob
      referenceFile_i_benchmark = File("%s_referenceFile_%s" %(benchmark, i))
      comparison_i_benchmark = File("compare_%s_%s.txt" %(benchmark, i))
      errorCompare_i_benchmark = File("compare_%s_%s.err" %(benchmark, i))

      #Add compare Job
      print "Adding compare_%s job for numberOfBasisFunctions=%s ..." %(benchmark, i)
      compareJob_i_benchmark = Job(name="compare_%s" %(benchmark))
      compareJob_i_benchmark.addArguments(referenceFile_i_benchmark, actualOutput_i_benchmark)
      compareJob_i_benchmark.setStdout(comparison_i_benchmark)
      compareJob_i_benchmark.setStderr(errorCompare_i_benchmark)
      compareJob_i_benchmark.uses(referenceFile_i_benchmark, link=Link.INPUT)
      compareJob_i_benchmark.uses(comparison_i_benchmark, link=Link.OUTPUT, transfer=True, register=False)
      dax.addJob(compareJob_i_benchmark)

      #when monitor finishes the slurm job is finished and output can be compared
      dax.depends(parent=monitorJob_i_benchmark, child=compareJob_i_benchmark)
            
    elif benchmark in ["TetraLOH4", "LOH1Scaling", "tpv16", "landers"]:
      #copy and extract reference receiver files to site
      receiversPacked = File("%s_referenceFiles_%s" %(benchmark, i))
      errorFileTarReceivers_i_benchmark = File("tarReceivers_%s_%s.err" %(benchmark, i))
      outputFileTarReceivers_i_benchmark = File("tarReceivers_%s_%s.out" %(benchmark, i))

      print "Adding unpackReceivers job for %s for numberOfBasisFunctions=%s ..." %(benchmark, i)
      print 'unpacked folder must be named "%s_referenceFiles_%s"' %(benchmark, i)
      #reference receiver files are expected to already have their rank info and ending removed
      unpackReceivers_i_benchmark = Job(name="tar")
      unpackReceivers_i_benchmark.addArguments("-x","-z","-v", "-f", receiversPacked)
      unpackReceivers_i_benchmark.setStdout(outputFileTarReceivers_i_benchmark)
      unpackReceivers_i_benchmark.setStderr(errorFileTarReceivers_i_benchmark)
      unpackReceivers_i_benchmark.uses(receiversPacked, link=Link.INPUT)
      dax.addJob(unpackReceivers_i_benchmark)
          
      #create new directory for new receiver files and prepare them for compare_receivers.py script
      remove_ranks_out_i_benchmark = File("remove_ranks_%s_%s.out" %(benchmark, i))
      remove_ranks_err_i_benchmark = File("remove_ranks_%s_%s.err" %(benchmark, i))
          
      remove_ranks_i_benchmark = Job(name="remove_ranks")
      remove_ranks_i_benchmark.addArguments("Benchmark_%s_%s/output"  %(benchmark, i))
      remove_ranks_i_benchmark.setStdout(remove_ranks_out_i_benchmark)
      remove_ranks_i_benchmark.setStderr(remove_ranks_err_i_benchmark)
      dax.addJob(remove_ranks_i_benchmark)
          
      #remove ranks job needs the slurm job to have finished
      dax.depends(parent=monitorJob_i_benchmark, child=remove_ranks_i_benchmark)
          
      #Add new files for compareJob
      comparison_i_benchmark = File("compare_%s_%s.out" %(benchmark, i))
      errorCompare_i_benchmark = File("compare_%s_%s.err" %(benchmark, i))
      plot_i_benchmark = File("plot_%s_%s.pdf" %(benchmark, i))
      csv_i_benchmark = File("values_%s_%s.csv" %(benchmark, i))
  
      #Add compare Job
      print "Adding compare_%s job for numberOfBasisFunctions=%s ..." %(benchmark, i)
      compareJob_i_benchmark = Job(name="compare_receivers")
      compareJob_i_benchmark.addArguments( "%s_referenceFiles_%s/" %(benchmark, i), 
                                           "Benchmark_%s_%s/output/receivers/" %(benchmark, i), 
                                           "plot_%s_%s.pdf" %(benchmark, i), 
                                           "values_%s_%s.csv" %(benchmark, i))
      compareJob_i_benchmark.setStdout(comparison_i_benchmark)
      compareJob_i_benchmark.setStderr(errorCompare_i_benchmark)
      compareJob_i_benchmark.uses(plot_i_benchmark, link=Link.OUTPUT, transfer=True, register=False)
      compareJob_i_benchmark.uses(csv_i_benchmark, link=Link.OUTPUT, transfer=True, register=False)
      dax.addJob(compareJob_i_benchmark)
          
      #Compare Job needs remove ranks job to have finished
      dax.depends(parent=remove_ranks_i_benchmark, child=compareJob_i_benchmark)
      dax.depends(parent=unpackReceivers_i_benchmark, child=compareJob_i_benchmark)

    else:
      print "no postprocessing defined for %s!" %(benchmark)

#Write the DAX to stdout
print "Writing %s" % l_arguments['dax']
f = open(l_arguments['dax'], "w")
dax.writeXML(f)
f.close()
