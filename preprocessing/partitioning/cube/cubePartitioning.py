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

import collections
import math
import operator
import os

def main():
    # Parse command line arguements
    args = Args()
    
    size = args.size()
    partitions = args.partitions()
    
    if isinstance(partitions, collections.Iterable):
        pd = partitions
    else:
        # try to find good partitions
        pd = [0] * 3
    
        pd[0] = int(math.ceil(math.pow(partitions, 1./3.)))
        while partitions % pd[0] != 0:
            pd[0] -= 1
        
        p = partitions/pd[0]
        
        pd[1] = int(math.ceil(math.sqrt(p)))
        while p % pd[1] != 0:
            pd[1] -=  1
        
        pd[2] = p/pd[1]
            
        pd.sort()
    
        # Get sorted indices of the sizes and sort partitions in the same order
        sortedSizes = [i[0] for i in sorted(enumerate(size),  key=operator.itemgetter(1))]
        pd = [pd[sortedSizes[i]] for i in range(3)]
    
    maxp = pd[0] * pd[1] * pd[2]
    
    print 'Number of partitions in each direction: '+', '.join(map(str, pd))+'; total partitions: '+str(maxp)
    
    # Write partition
    file = args.output()
    
    for z in range(size[2]):
        pz = z / ((size[2] + pd[2] - 1) / pd[2])
        for y in range(size[1]):
            py = y / ((size[1] + pd[1] - 1) / pd[1])
            for x in range(size[0]):
                px = x / ((size[0] + pd[0] - 1) / pd[0])
                
                p = px + (py + pz * pd[1]) * pd[0]
                
                if p < 0 or p >= maxp:
                    raise IOError("Wrong partition number computed: "+str(p))
                
                for i in range(5):
                    print >> file, p
                
    file.close()

if __name__ == '__main__':
    main()
