#! /usr/bin/python
##
# @file
# This file is part of SeisSol.
#
# @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
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
# Checks if PAPI is installed and sets up include path and library path
#

import os
import sys

import SCons

def find_file(file, pathes):
    """Search for a file in a list of given pathes"""
    
    for path in pathes:
        if os.path.exists(os.path.join(path, file)):
            return path
        
    return None

def find_lib(env, name, pathes):
    """Search for a library in a list of given pathes"""
    
    prefixes = SCons.Util.flatten(env.subst_list('$LIBPREFIXES'))
    suffixes = SCons.Util.flatten(env.subst_list('$LIBSUFFIXES'))
    
    for prefix in prefixes:
        for suffix in suffixes:
            path = find_file(str(prefix)+name+str(suffix), pathes)
            if path:
                return path
            
    return None

def generate(env):
    papi_base = None
    
    try:
        papi_base = env['PAPI_BASE']
    except KeyError:
        # PAPI_BASE not set
        # check os.environ, PATH
        if 'PAPI_BASE' in os.environ:
            papi_base = os.environ['PAPI_BASE']
        elif env.WhereIs('papi_version'):
            papi_base = os.path.join(os.path.dirname(env.WhereIs('papi_version')), os.pardir)
            papi_base = os.path.normpath(papi_base)
            
    # Set up possible library and include pathes
    include_dirs = []
    library_dirs = []
    
    if papi_base:
        include_dirs.append(os.path.join(papi_base, 'include'))
        if os.path.exists(os.path.join(papi_base, 'lib64')):
            library_dirs.append(os.path.join(papi_base, 'lib64'))
        library_dirs.append(os.path.join(papi_base, 'lib'))
        
    if 'CPLUS_INCLUDE_PATH' in os.environ:
        include_dirs.extend(os.environ['CPLUS_INCLUDE_PATH'].split(os.path.pathsep))
        
    if 'LIBRARY_PATH' in os.environ:
        library_dirs.extend(os.environ['LIBRARY_PATH'].split(os.path.pathsep))
    
    # Set the include and library path we found in the environment, so
    # they are accessible form the SConstruct file.
    # The values can be None. In this case PAPI is either installed in
    # a default path (e.g. /usr/include, /usr/lib) or not installed at
    # all.
    env['PAPI_INCLUDE_PATH'] = find_file('papi.h', include_dirs)
    env['PAPI_LIBRARY_PATH'] = find_lib(env, 'papi', library_dirs)
    
    # Add include/library path if found and link with PAPI
    if env['PAPI_INCLUDE_PATH']:
        env.Append(CPPPATH=[env['PAPI_INCLUDE_PATH']])
    if env['PAPI_LIBRARY_PATH']:
        env.Append(LIBPATH=[env['PAPI_LIBRARY_PATH']])
        env.Append(RPATH=[env['PAPI_LIBRARY_PATH']])
    env.Append(LIBS=['papi'])

def exists(env):
    """Always returns true. This tool does not do any checks, it just
     tries to add the correct include and library path"""
    return True
