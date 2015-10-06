#! /usr/bin/python
##
# @file
# This file is part of SeisSol.
#
# @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
#
# @section LICENSE
# Copyright (c) 2014, SeisSol Group
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

hdf5_fortran_prog_src = """
program HDF5_Test
    use hdf5
end program HDF5_Test
"""

hdf5_prog_src_serial = """
#include <hdf5.h>

int main() {
    H5open();
    
    return 0;
}
"""

def CheckHDF5FortranInclude(context):
    context.Message("Checking for Fortran HDF5 module... ")
    ret = context.TryCompile(hdf5_fortran_prog_src, '.f90')
    context.Result(ret)
    
    return ret

def CheckHDF5Linking(context, message):
    context.Message(message+"... ")
    # TODO serial/parallel switch
    ret = context.TryLink(hdf5_prog_src_serial, '.c')
    context.Result(ret)
    
    return ret

def CheckV18API(context):
    context.Message('Checking for HDF5 v18 API... ')
    ret = context.TryAction('h5cc -showconfig | fgrep v18')[0]
    context.Result(ret)
    
    return ret

def generate(env, **kw):
    conf = env.Configure(custom_tests = {'CheckHDF5FortranInclude' : CheckHDF5FortranInclude,
                                         'CheckHDF5Linking' : CheckHDF5Linking,
                                         'CheckV18API' : CheckV18API})
    
    if 'required' in kw:
        required = kw['required']
    else:
        required = False
        
    if 'fortran' in kw:
        fortran = kw['fortran']
    else:
        fortran = False
    
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
        
    ret = [conf.CheckLibWithHeader('hdf5_hl', 'hdf5.h', 'c'), conf.CheckLib('hdf5')]
    if not all(ret):
        if required:
            print 'Could not find HDF5 or zlib!'
            env.Exit(1)
        else:
            conf.Finish()
            return
        
    ret = conf.CheckHDF5Linking("Checking whether shared HDF5 library is used")
    if not ret:
        # Static library, link zlib as well
        ret = conf.CheckLib('z')
        if not ret:
            if required:
                print 'Could not find zlib!'
                env.Exit(1)
            else:
                conf.Finish()
                return
 
        # Try to find all additional libraries
        conf.CheckLib('curl')
        conf.CheckLib('gpfs')
        conf.CheckLib('dl')
 
        ret = conf.CheckHDF5Linking("Checking whether all HDF5 dependencies are found")
        if not ret:
            if required:
                print 'Could not find all HDF5 dependencies!'
                env.Exit(1)
            else:
                conf.Finish()            
                return
    
    # Check API Mapping
    if not conf.CheckV18API():
        # TODO We might need to extent this list
        conf.env.Append(CPPDEFINES=['H5Dcreate_vers=2',
                                    'H5Dopen_vers=2',
                                    'H5Acreate_vers=2',
                                    'H5Eget_auto_vers=2',
                                    'H5Eset_auto_vers=2'])
            
    conf.Finish()

def exists(env):
    return True
