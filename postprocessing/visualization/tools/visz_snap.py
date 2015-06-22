#!/usr/bin/env python
##
# @file
# This file is part of SeisSol.
#
# @author Simon Kremers (kremers AT geophysik.uni-muenchen.de)
# @author Christian Pelties (pelties AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/pelties)
# @author Tobias Megies (tobias.megies AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/megies)
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
# visualize Seissol output

import numpy as np
import sys
import matplotlib.pyplot as plt
from StringIO import StringIO

try:
    infile = sys.argv[1]
except:
    raise NameError('no input file!')

# read file header
lines = open(infile).readlines()
line = lines.pop(0)
info = line[line.find('"')+1:line.rfind('"')]

line = lines.pop(0)
variables = line.replace('"', '').split()[2:]
num_var = len(variables)
titl = variables
titl.remove("x")
titl.remove("y")
num_titl = len(titl)

line = lines.pop(0)
num_cells = int(line.split()[2])

msg = ["%s: %s" % (i+1, item) for i, item in enumerate(titl)]
msg = "What component do you want to look at?\n" + "\n".join(msg) + "\n"

comp = -1
while not 0 <= comp < num_titl:
    try:
        comp = int(raw_input(msg)) - 1
    except:
        pass

# read in data matrix
data = "".join(lines[:num_cells])
data = StringIO(data)
data = np.loadtxt(data, dtype="float").T

x = data[0]
y = data[1]
data = data[2:]

# read in connectivity matrix
triangles = "".join(lines[num_cells:])
triangles = StringIO(triangles)
triangles = np.loadtxt(triangles, dtype="int")
# matplotlib expects the first point to have index 0,
# in the file it has index 1!
triangles -= 1

plt.figure()
ax = plt.gca()
ax.set_aspect('equal')
plt.tripcolor(x, y, triangles, data[comp, :], shading='faceted')
cb = plt.colorbar()
cb.set_label(titl[comp])
#plt.triplot(x, y, triangles)
plt.title('%s -- %s component of wavefield' % (info, titl[comp]))
plt.xlabel('X [m]')
plt.ylabel('Y [m]')

plt.show()
