#! /usr/bin/python
##
# @file
# This file is part of SeisSol.
#
# @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
#
# @section LICENSE
# Copyright (c) 2013-2017, SeisSol Group
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
# Finds netCDF and add inlcude pathes, libs and library pathes
#
import sys
sys.path.append('../../../site_scons')


import utils.checks

netcdf_prog_src_serial = """
#include <netcdf.h>

int main() {
    int ncFile;
    nc_open("", 0, &ncFile);
    
    return 0;
}
"""

netcdf_prog_src_parallel = """
#include <netcdf_par.h>

int main() {
    int ncFile;
    nc_open_par("", 0, MPI_COMM_WORLD, MPI_INFO_NULL, &ncFile);
    
    return 0;
}
"""

def CheckNetcdfLinking(context, parallel, message):
    context.Message(message+"... ")
    if parallel:
        src = netcdf_prog_src_parallel
    else:
        src = netcdf_prog_src_serial
    ret = context.TryLink(src, '.c')
    context.Result(ret)
    
    return ret

def generate(env, **kw):
    conf = env.Configure(custom_tests = {'CheckLib': utils.checks.CheckLib,
                                         'CheckLibWithHeader': utils.checks.CheckLibWithHeader,
                                         'CheckNetcdfLinking' : CheckNetcdfLinking})
    
    if 'parallel' in kw and kw['parallel']:
        parallel = True
        header = 'netcdf_par.h'
    else:
        parallel = False
        header = 'netcdf.h'
        
    if 'required' in kw:
        required = kw['required']
    else:
        required = False
        
    if not conf.CheckLibWithHeader('netcdf', header, 'c'):
        if required:
            print('Could not find netCDF!')
            env.Exit(1)
        else:
            conf.Finish()
            return
        
    ret = conf.CheckNetcdfLinking(parallel, "Checking whether shared netCDF library is used")
    if not ret:
        # Static library, link with HDF5 and zlib as well
        ret = [conf.CheckLib(lib) for lib in ['hdf5_hl', 'hdf5', 'z']]
        
        if not all(ret):
            if required:
                print('Could not find HDF5 or zlib!')
                env.Exit(1)
            else:
                conf.Finish()
                return

        # Try to find all additional libraries
        conf.CheckLib('curl')
        conf.CheckLib('gpfs')
        conf.CheckLib('dl')

        ret = conf.CheckNetcdfLinking(parallel, "Checking whether all netCDF dependencies are found")
        if not ret:
            if required:
                print('Could not find all netCDF dependencies!')
                env.Exit(1)
            else:
                conf.Finish()
                return
            
    conf.Finish()

def exists(env):
    return True
