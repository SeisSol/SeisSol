#! /usr/bin/python
##
# @file
# This file is part of SeisSol.
#
# @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
#
# @section LICENSE
# Copyright (c) 2016, SeisSol Group
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
# Finds ASAGI and add inlcude pathes, libs and library pathes
#

asagi_prog_src = """
#include <asagi.h>

int main() {
    asagi::Grid* grid = asagi::Grid::create();
    
    return 0;
}
"""

def CheckPKG(context, lib):
    context.Message('Checking for ASAGI... ')
    ret = context.TryAction('pkg-config --exists \'%s\'' % lib)[0]
    context.Result(ret)
    return ret

def CheckASAGILinking(context, message):
    context.Message(message+'... ')
    ret = context.TryLink(asagi_prog_src, '.cpp')
    context.Result(ret)
    
    return ret

def generate(env, **kw):
    conf = env.Configure(custom_tests = {'CheckPKG': CheckPKG,
                                         'CheckASAGILinking': CheckASAGILinking})
    
    if 'required' in kw:
        required = kw['required']
    else:
        required = False
        
    if 'parallel' in kw:
        parallel = kw['parallel']
    else:
        # Use parallel as default, since ASAGI makes not too much
        # sense in serial code 
        parallel = True
        
    if not parallel:
        env.Append(CPPDEFINES=['ASAGI_NOMPI'])
        
    lib = 'asagi' if parallel else 'asagi_nompi'
        
    ret = conf.CheckPKG(lib)
    if not ret:
        if required:
            print('Could not find ASAGI!')
            env.Exit(1)
        else:
            conf.Finish()
            return
        
    # Try shared first
    conf.env.ParseConfig('pkg-config --cflags --libs '+lib)
        
    ret = conf.CheckASAGILinking('Checking for shared ASAGI library')
    if not ret:
        conf.env.ParseConfig('pkg-config --cflags --libs --static '+lib)
        
        ret = conf.CheckASAGILinking('Checking for static ASAGI library')
        if not ret:
            if required:
                print('Could not find ASAGI!')
                env.Exit(1)
            else:
                conf.Finish()
                return
            
    conf.Finish()

def exists(env):
    return True
