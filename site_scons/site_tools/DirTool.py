#! /usr/bin/python
##
# @file
# This file is part of SeisSol.
#
# @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
#
# @section LICENSE
# Copyright (c) 2013-2016, SeisSol Group
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
# Adds all xxxDir pathes in env to lib pathes and include pathes
#

import os

def generate(env, **kw):
    if 'rpath' in kw:
        rpath = kw['rpath']
    else:
        rpath = True

    if 'fortran' in kw:
        fortran = kw['fortran']
    else:
        fortran = False

    # Collect all pathes here, otherwise we change the directory during iteration
    binPathes = []
    incPathes = []
    libPathes = []
    pkgPathes = []

    for var in env.Dictionary():
        if var.endswith('Dir') and var != 'Dir':
            value = env[var]

            binPath = os.path.join(value, 'bin')
            incPath = os.path.join(value, 'include')
            libPath = os.path.join(value, 'lib')
            pkgPath = [os.path.join(value, 'lib', 'pkgconfig'),
                       os.path.join(value, 'share', 'pkgconfig')]

            if os.path.exists(binPath):
                binPathes.append(binPath)

            if os.path.exists(incPath):
                incPathes.append(incPath)

            if os.path.exists(libPath):
                libPathes.append(libPath)

            for p in pkgPath:
                if os.path.exists(p):
                    pkgPathes.append(p)

    env.PrependENVPath('PATH', binPathes)

    env.AppendUnique(CPPPATH=incPathes)
    if fortran:
        env.AppendUnique(F90PATH=incPathes)
    env.AppendUnique(LIBPATH=libPathes)
    if rpath:
        env.AppendUnique(RPATH=libPathes)

    env.PrependENVPath('PKG_CONFIG_PATH', pkgPathes)

def exists(env):
    return True
