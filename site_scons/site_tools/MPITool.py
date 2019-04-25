#! /usr/bin/python
##
# @file
# This file is part of SeisSol.
#
# @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
#
# @section LICENSE
# Copyright (c) 2015, SeisSol Group
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
# Handles MPI environments
#

import os
from SCons.Variables import PathVariable

# Add more variables here for other MPI libraries
MPI_ENV_VARIABLES = {
    'CC': ['OMPI_CC', 'MPICH_CC', 'I_MPI_CC'],
    'CXX': ['OMPI_CXX', 'MPICH_CXX', 'I_MPI_CXX'],
    'F90': ['OMPI_FC', 'MPICH_F90', 'I_MPI_F90']}

def CheckCompiler(context, command):
    context.Message('Checking for '+command+'... ')
    # Use the version key word since this should be available
    # for all compilers
    ret = context.TryAction(command+' --version')[0]
    context.Result(ret)
    
    return ret

def generate(env, **kw):
    # Do we need to set the scons variables?
    if 'vars' in kw:
        kw['vars'].AddVariables(
            PathVariable( 'mpicc',
                          'MPI C compiler wrapper (default: mpicc)',
                          None,
                          PathVariable.PathAccept ),
                  
            PathVariable( 'mpicxx',
                          'MPI C++ compiler wrapper (default: mpiCC)',
                          None,
                          PathVariable.PathAccept ),
                        
            PathVariable( 'mpif90',
                          'MPI Fortran compiler wrapper (default: mpif90)',
                          None,
                          PathVariable.PathAccept ))
        
        return
    
    conf = env.Configure(custom_tests = {'CheckCompiler' : CheckCompiler})
    
    # Set default compilers
    # Use I_MPI_ROOT to detect special Intel MPI wrappers
    if not 'mpicc' in env:
        if env['compiler'].startswith('cray_'):
            env['mpicc'] = 'cc'
        elif 'I_MPI_ROOT' in env['ENV']:
            if env['compiler'] == 'intel':
                env['mpicc'] = 'mpiicc'
            else:
                env['mpicc'] = 'mpicc'
        else:
            env['mpicc'] = 'mpicc'
    if not 'mpicxx' in env:
        if env['compiler'].startswith('cray_'):
            env['mpicxx'] = 'CC'
        elif 'I_MPI_ROOT' in env['ENV']:
            if env['compiler'] == 'intel':
                env['mpicxx'] = 'mpiicpc'
            else:
                env['mpicxx'] = 'mpicxx'
        else:
            if conf.CheckCompiler('mpicxx'):
                env['mpicxx'] = 'mpicxx'
            else:
                env['mpicxx'] = 'mpiCC'
    if not 'mpif90' in env:
        if env['compiler'].startswith('cray_'):
            env['mpif90'] = 'ftn'
        elif 'I_MPI_ROOT' in env['ENV']:
            if env['compiler'] == 'intel':
                env['mpif90'] = 'mpiifort'
            else:
                env['mpif90'] = 'mpif90'
        else:
            env['mpif90'] = 'mpif90'
            
    # Check all MPI wrappers
    ret = [conf.CheckCompiler(env[cc]) for cc in ['mpicc', 'mpicxx', 'mpif90']]
    if not all(ret):
        print('Could not find all MPI wrappers!')
        env.Exit(1)

    # Update the environment and build environment
    for (compiler, wrapper) in [('CC', 'mpicc'), ('CXX', 'mpicxx'), ('F90', 'mpif90')]:
        # Set all known env variables
        for envVar in MPI_ENV_VARIABLES[compiler]:
            env['ENV'][envVar] = env[compiler]
            
        # Update the build environment
        env[compiler] = env[wrapper]
        
    conf.Finish()

def exists(env):
    return True
