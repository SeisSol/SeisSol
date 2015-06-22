#! /usr/bin/python
##
# @file
# This file is part of SeisSol.
#
# @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
#
# @section LICENSE
# Copyright (c) SeisSol Group
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

import SCons

def _pathListExists(key, value, env):
    for path in env[key]:
        SCons.Script.PathVariable.PathExists(key, path, env)
        
def _pathToList(value):
    if not value or value == 'none':
        return []
    
    return value.split(os.path.pathsep)

def generate(env, *kw):
    if 'prefixPath' in env:
        
        if 'rpath' in kw:
            rpath = kw['rpath']
        else:
            rpath = True
        
        # Append include/lib and add them to the list if they exist
        incPathes = [p for p in map(lambda p: os.path.join(p, 'include'), env['prefixPath']) if os.path.exists(p)]
        libPathes = [p for p in map(lambda p: os.path.join(p, 'lib'), env['prefixPath']) if os.path.exists(p)]
                    
        env.AppendUnique(CPPPATH=incPathes)
        env.AppendUnique(LIBPATH=libPathes)
        if rpath:
            env.AppendUnique(RPATH=libPathes)
    
    env['PREFIX_PATH_VARIABLE'] = ('prefixPath',
        'Used when searching for include files, binaries or libraries ( /prefix/path1'+os.path.pathsep+'/prefix/path2 )',
        None,
        _pathListExists,
        _pathToList)

def exists(env):
    return True
