#! /usr/bin/python
##
# @file
# This file is part of SeisSol.
#
# @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
#
# @section LICENSE
# Copyright (c) 2012, SeisSol Group
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
# Build parameters for MAC Research Cluster (TUM).
# partition "nvd" features:
#   - 4 nodes: dual socket Intel SandyBridge-EP Xeon E5-2670, 128 GB RAM,
#     two NVIDIA M2090 GPUs and FDR infiniband
#    
# partition "ati" features:
#   - 4 nodes: dual socket Intel SandyBridge-EP Xeon E5-2670, 128 GB RAM,
#     two AMD FirePro W8000 GPUs and FDR infiniband
#
# partition "wsm" features:
#   - 2 nodes: quad socket Intel Westmere-EX Xeon E7-4830, 512 GB RAM and
#     FDR infiniband (QDR speed due to PCIe 2.0)
#  
# partition "snb" features:
#   - 28 nodes: dual socket Intel SandyBridge-EP Xeon E5-2670, 128 GB RAM
#     and QDR infiniband
#
# partition "bdz" features:
#   - 19 nodes: quad socket AMD Bulldozer Opteron 6274, 256 GB RAM
#     and QDR infiniband
#

compileMode='release'
parallelization='hybrid'
arch='dnoarch'
order='4'
generatedKernels = 'yes'
compiler = 'gcc'
logLevel = 'info'

aux_libs_folder='/ss_aux_libs'

netcdf='yes'
netcdfDir=aux_libs_folder
hdf5='yes'
hdf5Dir=aux_libs_folder

##  additionally for puml mesh format
metis = 'yes'
metisDir=aux_libs_folder

##  optional for ASAGI
asagi='yes'
zlibDir=aux_libs_folder #e.g. <path_to_ASAGI>/build/lib/
