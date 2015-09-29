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
    'CC': ['OMPI_CC', 'MPICH_CC'],
    'CXX': ['OMPI_CXX', 'MPICH_CXX'],
    'F90': ['OMPI_FC', 'MPICH_F90']}

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
    
    # Set default compilers
    if not 'mpicc' in env:
        env['mpicc'] = 'mpicc'
    if not 'mpicxx' in env:
        env['mpicxx'] = 'mpiCC'
    if not 'mpif90' in env:
        env['mpif90'] = 'mpif90'
    
    # Update the environment and build environment
    for (compiler, wrapper) in [('CC', 'mpicc'), ('CXX', 'mpicxx'), ('F90', 'mpif90')]:
        # Set all known env variables
        for envVar in MPI_ENV_VARIABLES[compiler]:
            env['ENV'][envVar] = env[compiler]
            
        # Update the build environment
        env[compiler] = env[wrapper]

def exists(env):
    return True
