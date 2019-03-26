#! /usr/bin/python
##
# @file
# This file is part of SeisSol.
#
# @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
#
# @section LICENSE
# Copyright (c) 2014-2016, SeisSol Group
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
# Finds HDF5 and add inlcude pathes, libs and library pathes
#

import utils.checks
import utils.compiler
import utils.pkgconfig

hdf5_fortran_prog_src = """
program HDF5_Test
    use hdf5
end program HDF5_Test
"""

def CheckHDF5FortranInclude(context):
    context.Message("Checking for Fortran HDF5 module... ")
    ret = context.TryCompile(hdf5_fortran_prog_src, '.f90')
    context.Result(ret)

    return ret

def CheckV18API(context, h5cc):
    context.Message('Checking for HDF5 v18 API... ')
    if h5cc:
        ret = context.TryAction(h5cc + ' -showconfig | fgrep v18')[0]
    else:
        # FIXME currently assuming the default if the h5cc does not exist
        ret = True
    context.Result(ret)

    return ret

def generate(env, required = False, parallel = False, fortran = False, **kw):
    if env.GetOption('help') or env.GetOption('clean'):
        return

    conf = env.Configure(custom_tests = {'CheckProg': utils.checks.CheckProg,
                                         'CheckLibWithHeader': utils.checks.CheckLibWithHeader,
                                         'CheckHDF5FortranInclude' : CheckHDF5FortranInclude,
                                         'CheckV18API' : CheckV18API})

    # Find h5cc or h5pcc
    h5ccs = ['h5cc', 'h5pcc']
    if parallel:
        h5ccs[0], h5ccs[1] = h5ccs[1], h5ccs[0]
    for h5cc in h5ccs:
        h5cc = conf.CheckProg(h5cc)
        if h5cc:
            break

    if h5cc:
        # Parse the output from the h5cc compiler wrapper
        def parse_func(env, cmd):
            # remove the compiler
            cmd = cmd.partition(' ')[2]
            # remove unknown arguments
            cmd = utils.compiler.removeUnknownOptions(cmd)
            return env.ParseFlags(cmd)
        flags = env.ParseConfig([h5cc, '-show', '-shlib'], parse_func)
        if flags['LIBS']==[] or flags['LIBPATH']==[]:
           h5cc=None
    if not h5cc:
        # Try pkg-config
        hdf5s = ['hdf5_hl', 'hdf5_hl_parallel']
        if parallel:
            hdf5s[0], hdf5s[1] = hdf5s[1], hdf5s[0]
        for hdf5 in hdf5s:
            flags = utils.pkgconfig.parse(conf, hdf5)
            if flags:
                break

    if not flags:
            if required:
                print 'Could not find h5cc or h5pcc. Make sure the path to the HDF5 library is correct!'
                env.Exit(1)
            else:
                conf.Finish()
                return

    # Add the lib path
    env.AppendUnique(LIBPATH=flags['LIBPATH'])

    if fortran:
        # Fortran module file
        ret = conf.CheckHDF5FortranInclude()
        if not ret:
            if required:
                print 'Could not find HDF5 for Fortran!'
                env.Exit(1)
            else:
                conf.Finish()
                return

        # Fortran library
        ret = conf.CheckLib('hdf5_fortran')
        if not ret:
            if required:
                print 'Could not find HDF5 for Fortran!'
                env.Exit(1)
            else:
                conf.Finish()
                return

    ret = conf.CheckLibWithHeader(flags['LIBS'][0], 'hdf5.h', 'c', extra_libs=flags['LIBS'][1:])
    if not ret:
        if required:
            print 'Could not find the HDF5 library!'
            env.Exit(1)
        else:
            conf.Finish()
            return

    # Check API Mapping
    if not conf.CheckV18API(h5cc):
        # TODO We might need to extent this list
        conf.env.Append(CPPDEFINES=['H5Dcreate_vers=2',
                                    'H5Dopen_vers=2',
                                    'H5Gcreate_vers=2',
                                    'H5Gopen_vers=2',
                                    'H5Acreate_vers=2',
                                    'H5Eget_auto_vers=2',
                                    'H5Eset_auto_vers=2'])

    conf.Finish()

def exists(env):
    return True
