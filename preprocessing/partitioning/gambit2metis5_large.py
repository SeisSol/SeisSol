##
# @file
# This file is part of SeisSol.
#
# @author Marek Simon (http://www.geophysik.uni-muenchen.de/Members/msimon)
# @author Martin Kaeser (martin.kaeser AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/kaeser)
#
# @section LICENSE
# Copyright (c) 2008, SeisSol Group
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

import numpy as np
import os
import timeit
import sys
#import ipdb; ipdb.set_trace()

print '     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
print '     %%                                                             %%'
print '     %% Gambit2Metis converts a GAMBIT-meshfile of                  %%'
print '     %% tetrahedrons or hexahedrons                                 %%'
print '     %% stored as "filename.neu"                                    %%'
print '     %% into a METIS-file and calls METIS for mesh partitioning.    %%'
print '     %% METIS mesh partitioner:                                     %%'
print '     %% http://www-users.cs.umn.edu/~karypis/metis/metis/index.html %%'
print '     %%                                                             %%'
print '     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
print '\n \n'
filename = raw_input('     Specify GAMBIT-filename:   ')
proc = raw_input('     Specify number of processor of the partition:   ')
procs = proc.split()

t=timeit.time.time()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%        Read Gambit Data 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=open(filename+'.neu');
print '\n'
print '-----------------------------------------------------------------------------------'
print ' Reading data from: %s' % filename+'.neu'
for i in range(6):
    junk  = f.readline()
tmp=f.readline()
tmp=[int(s) for s in tmp.split() if s.isdigit()]
NUMNP  = tmp[0] #                     %Number of Vertices
NELEM  = tmp[1] #                     %Number of Elements 
NGRPS  = tmp[2] #                     %Number of Element Groups
NBSETS = tmp[3] #                     %Number of Boundary Condition Sets
NDFCD  = tmp[4] #                     %Number of Dimensions (2 or 3)
try:
    if 3 == NDFCD:
        f.readline()
        f.readline()
        vertices=np.zeros((4, NUMNP))
        for v in range(NUMNP):
            line = f.readline()
            numbers = line.split()
            vertices[:,v] = [int(numbers[0]),float(numbers[1]),float(numbers[2]),float(numbers[3])]
        print 'Vertices read'
        f.readline()
        f.readline()
        position=f.tell()
        first_line=f.readline()
        first_line=[int(s) for s in first_line.split() if s.isdigit()]
        if first_line[1] is 6:
            n=7
        if first_line[1] is 4:
            n=11
        f.seek(position, 0)
        elements = np.zeros((7, NELEM))
        for e in range(NELEM):
            line = f.readline()
            elements[:,e] = [int(num) for num in line.split()]
        print 'Read the elements'
        if elements[1][0] == 6:
            print 'Read all tetrahedrons successfully!'
            ncommon = '3'
        if elements[1][0] == 4:
            print 'Read all hexahedrons successfully!'
            ncommon = '4'
    
    elif 2 == NDFCD:
        f.readline()
        f.readline()
        vertices=np.fromfile(file=f, count=NUMNP*3, sep=" ")
        vertices=vertices.reshape(NUMNP, 3)
        vertices=vertices.transpose()
        f.readline()
        f.readline()
        position=f.tell()
        first_line=f.readline()
        first_line=[int(s) for s in first_line.split() if s.isdigit()]
        if first_line[1] is 3:
            n=6
        if first_line[1] is 2:
            n=7
        f.seek(position, 0)
        elements=np.fromfile(file=f, count=NELEM*n, sep=" ").reshape(NELEM, n).transpose()
        if elements[1][0] == 3:
            print 'Read all triangles successfully!'
            ncommon = '2'
        if elements[1][0] == 2:
            print 'Read all quadrilaterals successfully!'
            ncommon = '2'
except:
    print '##ERROR: The GAMBIT input file shows a problem\n\n ##', sys.exc_info()[0]
    print 'Line: ', line
#   return
f.close()

print '-----------------------------------------------------------------------------------'

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%        Writing Metis Data 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print 'Write METIS-file: %s' % filename+'.met'
f=open(filename+'.met','w')
if 3 == NDFCD:
    if int(elements[1][0]) is 6:
        f.write('%12i\n' % NELEM)
        #np.savetxt(f,  elements[3:n].transpose())
        for [a, b, c, d] in elements[3:n].transpose():
            f.write('%12i %12i %12i %12i \n' % (a, b, c, d))
    if int(elements[1][0]) is 4:
        f.write('%12i\n' % NELEM)
        #np.savetxt(f,  elements[3:n].transpose())    <------- is somehow slower ??
        for [a,b,c,d,e,ff,g,h] in elements[3:n].transpose():
            f.write('%12i %12i %12i %12i %12i %12i %12i %12i\n' % (a, b, c, d,e,ff,g,h))
elif 2 == NDFCD:
    if int(elements[1][0]) is 3:
        f.write('%12i\n' % NELEM)
        for [a, b, c] in elements[3:n].transpose():
            f.write('%12i %12i %12i\n' % (a, b, c))
    if int(elements[1][0]) is 2:
        f.write('%12i\n' % NELEM)
        for [a,b,c,d] in elements[3:n].transpose():
            f.write('%12i %12i %12i %12i\n' % (a, b, c, d))
f.close();
t=t-timeit.time.time()
print '-----------------------------------------------------------------------------------'
print ' Conversion finished successfully!  (%s CPU sec)\n' % t
print '-----------------------------------------------------------------------------------'
print '\n'

print '-----------------------------------------------------------------------------------'
print 'METIS partition starts!\n'
print '-----------------------------------------------------------------------------------'
print '\n'
os.system('m2gmetis -ncommon='+ncommon+' '+filename+'.met '+filename+'.met.dgraph')
for proc in procs :
  os.system('gpmetis -ptype=rb '+filename+'.met.dgraph '+proc)
  os.system('mv '+filename+'.met.dgraph.part.'+proc+'  '+filename+'.met.epart.'+proc)
  print 'Wrote final METIS-file: %s' % filename+'.met.epart.'+proc

#%clean up of files
os.system('rm '+filename+'.met')
os.system('rm '+filename+'.met.dgraph')    

print '-----------------------------------------------------------------------------------'
print 'METIS partition done!'
print '-----------------------------------------------------------------------------------'
print '\n'
print '-----------------------------------------------------------------------------------'
print 'Wrote final METIS-file: %s' % filename+'.met.epart.'+proc
print '-----------------------------------------------------------------------------------'
print '\n'

