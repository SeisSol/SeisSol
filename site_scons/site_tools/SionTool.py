#! /usr/bin/python
##
# @file
# This file is part of SeisSol.
#
# @author Gilbert Brietzke (brietzke AT lrz.de, http://www.lrz.de)
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
# uses sionconfig to add include pathes, libs and library pathes for sionlib


def generate(env, **kw):
    import commands
    import re
    import os 
    conf = env.Configure()

    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    
    print 'Checking for sionlib ...',

    if not('parallel' in kw and kw['parallel']):
        print (FAIL+'Failed!')
        print (WARNING+'Warning: sionlib: non mpi or hybrid build but sionlib requires mpi:'+ENDC)
        print (WARNING+'Warning: sionlib: ... will continue without sionlib'+ENDC)
        conf.env['sionlib']=False

    if conf.env['sionlib']:
        try:
            sionlibDir=conf.env['sionlibDir']
        except:
            pass
        if 'sionlibDir' in locals():
            if type(sionlibDir) == str:
                if os.path.isdir(sionlibDir):
                    if os.path.isfile(os.path.join(sionlibDir, 'bin','sionconfig')):
                        print ' using sionconfig in '+os.path.join(sionlibDir, 'bin')+' ... ',
                        sionconfig = os.path.join(sionlibDir, 'bin','sionconfig')
                    elif os.path.isfile(os.path.join(sionlibDir, '..','bin','sionconfig')):
                        print ' using sionconfig in '+os.path.join(sionlibDir, '..','bin')+' ...',
                        sionconfig = os.path.join(sionlibDir, '..','bin','sionconfig')
        else:
            print ' using sionconfig from $PATH ...',
            sionconfig = commands.getstatusoutput('which sionconfig')[1]
            
        if commands.getstatusoutput(sionconfig + ' --cxx --mpi --libs --64')[0] == 0: 
            out=commands.getstatusoutput(sionconfig + ' --cxx --mpi --libs --64')[1].split()
            for item in out:
                if re.match("-L",item):
                    conf.env.Append(LINKFLAGS=[item])
                if re.match("-l",item):
                    conf.env.Append(LIBS=[item[2:]])
                    
            out=commands.getstatusoutput(sionconfig + ' --cxx --mpi --cflags --64')[1].split()
            for item in out:
                if re.match("-I/",item):
                    conf.env.Append(CXXFLAGS=[item])

            conf.env.Append(F90FLAGS=['-DUSE_SIONLIB'])
            conf.env.Append(CPPDEFINES=['USE_SIONLIB'])
            print (OKGREEN+'yes'+ENDC)
        else:
            print (FAIL+'no!')
            print (WARNING+'Warning: sionlib: could not execute sionconfig'+ENDC)
            print (WARNING+'Warning: sionlib: please add location of sionconfig to your search path'+ENDC)
            print (WARNING+'Warning: sionlib: ... will continue without sionlib'+ENDC)
            conf.env['sionlib']=False

    conf.Finish()

def exists(env):
    return True
