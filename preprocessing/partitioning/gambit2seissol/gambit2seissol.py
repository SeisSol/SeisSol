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

from lib.args import Args
from lib.gambit import GambitReader
from lib.gambit import GambitWriter
from lib.partition import DummyPartition
from lib.partition import PartitionWriter
from lib.tmp import TmpDir

from partition.partitioner import Partitioner
from reorder.reorderer import Reorderer

import os
import sys
import shutil

def main():
    # Parse command line arguements
    args = Args()
    
    # Create a temporary working directory
    tmpdir = TmpDir('gambit2seissol', args.tmpDir())
    
    # Read the Gambit file
    try:
        mesh = GambitReader(args.inputFile(), args.noReorder())
    except IOError, e:
        print >> sys.stderr, 'Could not parse GAMBIT file'
        print >> sys.stderr, str(e)
        sys.exit(1)
    
    if args.partitions() > 1:
        print 'Starting partitioner'
        try:
            partitioner = Partitioner(mesh, args.partitions(), tmpdir)
            partition = partitioner.partition()
        except Exception, e:
            print >> sys.stderr, 'Could not partition mesh'
            print >> sys.stderr, str(e)
            sys.exit(1)
    else:
        partition = DummyPartition(len(mesh.elements()))
        
    if not args.noReorder():
        print 'Starting reorderer'
        try:
            reorderer = Reorderer(mesh, partition)
            mesh = reorderer.mesh()
            partition = reorderer.partition()
        except Exception, e:
            print >> sys.stderr, 'Could not reorder mesh'
            print >> sys.stderr, str(e)
            sys.exit(2)
        
    # Filenames
    basename = os.path.basename(args.inputFile().name)
    if basename[-4:] == '.neu':
        basename = basename[:-4]
        
    if not args.noReorder():
        basename += '.'+str(args.partitions())
    
    # Write the new mesh and partition
    outputFilename = os.path.join(args.outputDir(), basename+'.neu')
    if args.noReorder():
        if not os.path.exists(outputFilename) or not os.path.samefile(outputFilename, args.inputFile().name):
            # No reorder but file does not exist -> copy it
            shutil.copyfile(args.inputFile().name, outputFilename)
    else:
        # Write Gambit mesh only if it is not the the input file
        GambitWriter(mesh, os.path.join(args.outputDir(), basename+'.neu'))
    PartitionWriter(partition, os.path.join(args.outputDir(), basename
        +'.met.epart.'+str(args.partitions())))

if __name__ == '__main__':
    main()
