#!/usr/bin/python
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

import argparse
import os
import sys

def writableDir(dir):
    if not os.path.isdir(dir):
        raise IOError("{0} is not a valid directory".format(dir))
    
    if os.access(dir, os.W_OK):
        return dir
    
    raise IOError("{0} is not a writable directory".format(dir))

class Args:
    """Parses the command line arguments using
    argparse.ArgumentParser"""
    
    def __init__(self):
        parser = argparse.ArgumentParser(prog='gambit2seissol')
        parser.add_argument('-i', '--input', type=argparse.FileType('rU'),
            required=True, help='gambit mesh file', metavar='FILE')
        parser.add_argument('-p', '--partitions', type=int,
            required=True, help='number of partitions, set to 1 to skip partitioning')
        parser.add_argument('-o', '--output', type=writableDir,
            required=True, help='output directory for resulting Gambit and partition file',
            metavar='DIR')
        parser.add_argument('-r', '--no-reorder', action='store_true',
            help='do not reorder the mesh using Zoltan')
        parser.add_argument('-t', '--tmp', type=writableDir, default=None,
            help='tmp directory that should be used instead of the default', metavar='DIR')
        parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
        
        try:
            # Parse cmd line options
            self.__options = parser.parse_args(sys.argv[1:])
        except IOError, e:
            parser.error(str(e))
            
    def inputFile(self):
        return self.__options.input
    
    def partitions(self):
        return self.__options.partitions
    
    def outputDir(self):
        return self.__options.output
    
    def noReorder(self):
        return self.__options.no_reorder
    
    def tmpDir(self):
        return self.__options.tmp
