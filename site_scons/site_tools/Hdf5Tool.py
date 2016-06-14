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

import SCons.SConf
import utils.compiler

hdf5_fortran_prog_src = """
program HDF5_Test
    use hdf5
end program HDF5_Test
"""

def CheckProg(context, prog_name):
    """
    This function is from the latest version of SCons to support
    older SCons version.
    Configure check for a specific program.
    Check whether program prog_name exists in path.  If it is found,
    returns the path for it, otherwise returns None.
    """

    context.Message("Checking whether %s program exists..." % prog_name)
    path = context.env.WhereIs(prog_name)
    context.Result(bool(path))

    return path

def CheckLibWithHeader(context, libs, header, language,
                       call = None, extra_libs = None, autoadd = 1):
     """
     This function is from SCons but extended with additional flags, e.g.
     the extra_libs.
     Another (more sophisticated) test for a library.
     Checks, if library and header is available for language (may be 'C'
     or 'CXX'). Call maybe be a valid expression _with_ a trailing ';'.
     As in CheckLib, we support library=None, to test if the call compiles
     without extra link flags.
     """
     prog_prefix, dummy = \
                  SCons.SConf.createIncludesFromHeaders(header, 0)
     if libs == []:
         libs = [None]

     if not SCons.Util.is_List(libs):
         libs = [libs]

     res = SCons.Conftest.CheckLib(context, libs, None, prog_prefix,
             call = call, language = language, extra_libs = extra_libs,
             autoadd = autoadd)
     context.did_show_result = 1
     return not res

def CheckHDF5FortranInclude(context):
    context.Message("Checking for Fortran HDF5 module... ")
    ret = context.TryCompile(hdf5_fortran_prog_src, '.f90')
    context.Result(ret)
    
    return ret

def CheckV18API(context, h5cc):
    context.Message('Checking for HDF5 v18 API... ')
    ret = context.TryAction(h5cc + ' -showconfig | fgrep v18')[0]
    context.Result(ret)
    
    return ret

def generate(env, required = False, parallel = False, fortran = False, **kw):
    if env.GetOption('help') or env.GetOption('clean'):
        return

    conf = env.Configure(custom_tests = {'CheckProg': CheckProg,
                                         'CheckLibWithHeader': CheckLibWithHeader,
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

    if not h5cc:
        if required:
            print 'Could not find h5cc or h5pcc. Make sure the path to the HDF5 library is correct!'
            env.Exit(1)
        else:
            conf.Finish()
            return

    # Parse the output from the h5cc compiler wrapper
    def parse_func(env, cmd):
        # remove the compiler
        cmd = cmd.partition(' ')[2]
	# remove unknown arguments
	cmd = utils.compiler.removeUnknownOptions(cmd)
        return env.ParseFlags(cmd)
    flags = env.ParseConfig([h5cc, '-show', '-shlib'], parse_func)

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
                                    'H5Acreate_vers=2',
                                    'H5Eget_auto_vers=2',
                                    'H5Eset_auto_vers=2'])
            
    conf.Finish()

def exists(env):
    return True
