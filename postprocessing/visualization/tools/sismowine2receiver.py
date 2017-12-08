#!/usr/bin/env python
##
# @file
# This file is part of SeisSol.
#
# @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
#
# @section LICENSE
# Copyright (c) SeisSol Group
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
# usage: 
# sismowine2receiver <directory-containing-x.dat-y.dat-and-z.dat> <output-file>

import sys
import numpy

class SismowineFile:
  def __init__(self, filename):
    f = open(filename, 'r')
    header = f.readline()
    header = header.split()
    self.nsamples = int(header[0])
    self.dt = float(header[1])
    self.nrec = int(header[2])
    self.data = numpy.empty([self.nsamples, self.nrec], order='F', dtype=float)
    for rec in range(self.nrec):
      for sample in range(self.nsamples):
        self.data[sample, rec] = float(f.readline())
  
  

xdat = SismowineFile(sys.argv[1] + '/x.dat')
ydat = SismowineFile(sys.argv[1] + '/y.dat')
zdat = SismowineFile(sys.argv[1] + '/z.dat')

for rec in range(xdat.nrec):
  fo = open('{0}-{1:03d}.dat'.format(sys.argv[2], rec+1), 'w')
  fo.write('TITLE = "Reference for receiver number {0:03d}"\n'.format(rec))
  fo.write('VARIABLES = "Time","u","v","w"\n')
  fo.write('# x1 NaN\n')
  fo.write('# x2 NaN\n')
  fo.write('# x3 NaN\n')
  for sample in range(xdat.nsamples):
    fo.write('{} {} {} {}\n'.format(sample * xdat.dt, xdat.data[sample, rec], ydat.data[sample, rec], zdat.data[sample, rec]))
