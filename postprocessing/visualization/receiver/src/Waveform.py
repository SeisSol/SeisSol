# SPDX-License-Identifier: BSD-3-Clause
##
# @file
# This file is part of SeisSol.
#
# @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
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
#

import numpy
import re

# the lines of the header of a fault receiver giving the stress the fault starts out
# under, and the quantity each of them is the initial value of
INITIAL_STRESS = {'P_0': 'P_n', 'T_s': 'T_s', 'T_d': 'T_d'}

class Waveform:
  def __init__(self, names, data, coordinates, simulation = None):
    data = numpy.array(data)

    self.waveforms = dict()
    self.norm = dict()
    self.show = dict()
    # which of the fused simulations this is, or None for a run of a single one
    self.simulation = simulation

    for i in range(0, len(names)):
      if names[i] == 'Time':
        self.time = data[:,i]
      else:
        self.waveforms[ names[i] ] = data[:,i]
        self.norm[ names[i] ] = numpy.max(numpy.abs(data[:,i]))
        self.show[ names[i] ] = True

    self.coordinates = numpy.array(coordinates)

  def subtract(self, other):
    newTime = numpy.union1d(self.time, other.time)
    self.waveforms = { key: value for key,value in self.waveforms.items() if key in other.waveforms }
    for name, wf in self.waveforms.items():
      wf0 = numpy.interp(newTime, other.time, other.waveforms[name])
      wf1 = numpy.interp(newTime, self.time, self.waveforms[name])
      self.waveforms[name] = wf1 - wf0
    self.time = newTime

  def normalize(self):
    for name, wf in self.waveforms.items():
      if self.norm[name] > numpy.finfo(float).eps:
        self.waveforms[name] = wf / self.norm[name]

  def differentiate(self):
    for name, wf in self.waveforms.items():
      self.waveforms[name] = numpy.gradient(wf, self.time)

  def integrate(self):
    dt = self.time[1] - self.time[0]
    for name, wf in self.waveforms.items():
      self.waveforms[name] = numpy.cumsum(wf) * dt

def splitSimulations(names, data):
  """Splits the columns of a receiver by simulation.

  Returns a dict from the simulation, counted from zero, to its names and columns,
  the time first; a run of a single simulation is the one entry None. A fused run
  writes a row per sample and simulation, which a column SimulationIndex names.
  Before, it wrote a column per quantity and simulation, with the simulation in
  the name: counted from zero and appended in the volume (v10, v11, ...), counted
  from one after a dash on the fault (SRs-1, SRs-2, ...).
  """
  data = numpy.asarray(data, dtype=float)
  if 'SimulationIndex' in names:
    column = names.index('SimulationIndex')
    kept = [i for i in range(len(names)) if i != column]
    simulations = data[:, column].astype(int)
    return { int(simulation): ([names[i] for i in kept], data[simulations == simulation][:, kept])
             for simulation in numpy.unique(simulations) }

  quantities = names[1:]
  dashed = [re.fullmatch(r'(.+)-(\d+)', name) for name in quantities]
  if quantities and all(dashed):
    columns = dict()
    for i, match in enumerate(dashed):
      columns.setdefault(int(match.group(2)) - 1, []).append((match.group(1), i + 1))
    return { simulation: (['Time'] + [name for name, _ in entries],
                          data[:, [0] + [i for _, i in entries]])
             for simulation, entries in sorted(columns.items()) }

  # in the volume, a block of columns per simulation, all of them ending in its index
  for count in range(2, len(quantities) + 1):
    if len(quantities) % count != 0:
      continue
    width = len(quantities) // count
    blocks = [quantities[k * width:(k + 1) * width] for k in range(count)]
    if not all(name.endswith(str(k)) for k, block in enumerate(blocks) for name in block):
      continue
    bases = [[name[:len(name) - len(str(k))] for name in block] for k, block in enumerate(blocks)]
    if all(base == bases[0] for base in bases):
      return { k: (['Time'] + bases[0], data[:, [0] + list(range(1 + k * width, 1 + (k + 1) * width))])
               for k in range(count) }

  return { None: (names, data) }

def addInitialStress(names, data, stresses):
  """Adds the stress a fault starts out under to the tractions, which hold its change.

  stresses maps the name of a header line (P_0, T_s, T_d) to its value.
  """
  for header, quantity in INITIAL_STRESS.items():
    if quantity in names and header in stresses:
      data[:, names.index(quantity)] += stresses[header]
