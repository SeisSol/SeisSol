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
#
# visualize Seissol output

import numpy as np
import sys
import matplotlib.pyplot as plt

try:
    infile = sys.argv[1]
except:
    raise NameError('no input file!')

# read file header
f = open(infile)
line1 = f.readline()
line2 = f.readline()
f.close()
variables = line1.replace('"', '').split()[2:]
num_x = int(line2.split()[2])
num_y = int(line2.split()[4])
num_var = len(variables)
titl = variables
titl.remove("x")
titl.remove("y")
num_titl = len(titl)

msg = ["%s: %s" % (i+1, item) for i, item in enumerate(titl)]
msg = "What component do you want to look at?\n" + "\n".join(msg) + "\n"

comp = -1
while not 0 <= comp < num_titl:
    try:
        comp = int(raw_input(msg)) - 1
    except:
        pass

# read in data matrix
data = np.loadtxt(infile, skiprows=2).T
x = data[0].reshape((num_y, num_x))
y = data[1].reshape((num_y, num_x))
data = data[2:]
z = data[comp].reshape((num_y, num_x))
z = np.flipud(z)

plt.figure()
ax = plt.gca()
ax.set_aspect('equal')
plt.imshow(z, extent=(x.min(), x.max(), y.min(), y.max()), interpolation="nearest")
cb = plt.colorbar()
cb.set_label(titl[comp])
#plt.triplot(x, y, triangles)
plt.title('%s component of wavefield' % titl[comp])
plt.xlabel('X [m]')
plt.ylabel('Y [m]')

plt.show()
