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

import Waveform
import re

def read(fileName):
  """Reads a text receiver file, as a list of waveforms, one per simulation."""
  data = []
  coordinates = [float('nan')] * 3
  coordComment = re.compile(r'#\s*x(\d)\s+([0-9\.eE\+\-]+)')
  # "P_0 <value> <value> ...": a value per simulation. Files of a fused run of an older
  # version have a line per simulation instead, counted from one, e.g. "P_03 <value>".
  stressComment = re.compile(r'#\s*(P_0|T_s|T_d)(\d*)\s+(.*)$')

  stresses = {}  # e.g. {('P_0', 0): value, ('T_s', 1): value, ...}
  variables = []

  with open(fileName) as f:
    for row in f:
      row = row.strip()
      if not row:
        continue
      if row[0] in '-+0123456789.':
        values = row.split()
        data.append([float(x) for x in values])
      elif row.startswith('#'):
        match_coord = coordComment.match(row)
        if match_coord:
          coordinates[int(match_coord.group(1))-1] = float(match_coord.group(2))
          continue
        match_stress = stressComment.match(row)
        if match_stress:
          key, number, values = match_stress.groups()
          values = [float(x) for x in values.split()]
          if number:
            stresses[(key, int(number) - 1)] = values[0]
          else:
            for simulation, value in enumerate(values):
              stresses[(key, simulation)] = value
      elif row.startswith('VARIABLES'):
        var_line = row.split('=', 1)[1]
        variables = [v.strip().strip('"') for v in var_line.split(',')]

  if not data:
    return None

  n_cols = len(data[0])
  # Fallback: if variable list not matching, auto-generate unnamed
  if len(variables) < n_cols:
    variables += [f"unnamed_{i}" for i in range(len(variables), n_cols)]
  elif len(variables) > n_cols:
    variables = variables[:n_cols]

  waveforms = []
  for simulation, (names, values) in Waveform.splitSimulations(variables, data).items():
    index = 0 if simulation is None else simulation
    Waveform.addInitialStress(names, values, { key: value for (key, sim), value in stresses.items() if sim == index })
    waveforms.append(Waveform.Waveform(names, values, coordinates, simulation))
  return waveforms
