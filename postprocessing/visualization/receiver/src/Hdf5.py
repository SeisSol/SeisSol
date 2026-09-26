# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

import numpy
import os
import re

import Waveform

def isReceiverFile(fileName):
  return re.fullmatch(r'.+-(fault)?receivers\.h5', os.path.basename(fileName)) is not None

def read(fileName):
  """Reads an HDF5 receiver file, as (label, waveforms) per receiver, one waveform per simulation.

  The label is the name the receiver has as a text file. A row of the table is one
  receiver of one simulation, which SimulationIndex names; a table of off-fault
  receivers written before it had that column holds all simulations in a row, as
  the text files did (see Waveform.splitSimulations).
  """
  try:
    import h5py
  except ImportError:
    print('Reading {} needs h5py.'.format(fileName))
    return []

  match = re.fullmatch(r'(.+)-((fault)?receiver)s\.h5', os.path.basename(fileName))
  if match is None:
    return []
  prefix, kind = match.group(1), match.group(2)

  receivers = dict()
  with h5py.File(fileName, 'r') as handle:
    group = handle.get(kind + 's')
    if group is None:
      return []
    index = group['Index'][:]
    numbers = group['PointId' if kind == 'receiver' else 'ReceiverId'][:]
    simulations = group['SimulationIndex'][:] if 'SimulationIndex' in group else None
    coordinates = group['Coordinates'][:]
    # P_0, T_s, T_d of the fault receivers, as the header of a text file gives them
    stresses = group['InitialStress'][:] if 'InitialStress' in group else None

    tables = dict()
    for row, (table, column) in enumerate(index):
      if table not in tables:
        tables[table] = group['group{}'.format(table)][:]
      samples = tables[table][:, column]
      # a receiver that took fewer samples than the longest one of its table is padded with NaN
      samples = samples[numpy.isfinite(samples['Time'])]
      names = list(samples.dtype.names)
      data = numpy.column_stack([samples[name].astype(float) for name in names])

      if simulations is None:
        parts = Waveform.splitSimulations(names, data)
      else:
        parts = { int(simulations[row]): (names, data) }
      for simulation, (partNames, partData) in parts.items():
        if stresses is not None:
          Waveform.addInitialStress(partNames, partData, dict(zip(Waveform.INITIAL_STRESS, stresses[row])))
        receivers.setdefault(int(numbers[row]), []).append(
          Waveform.Waveform(partNames, partData, coordinates[row], simulation))

  result = []
  for number in sorted(receivers):
    waveforms = sorted(receivers[number], key=lambda wf: -1 if wf.simulation is None else wf.simulation)
    # a run of a single simulation has none to tell apart
    if len(waveforms) == 1:
      waveforms[0].simulation = None
    result.append(('{}-{}-{:05d}'.format(prefix, kind, number + 1), waveforms))
  return result
