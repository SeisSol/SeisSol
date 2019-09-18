#! /usr/bin/env python
##
# @file
# This file is part of SeisSol.
#
# @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
#
# @section LICENSE
# Copyright (c) 2016, SeisSol Group
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

import os
import arch
import re

class Candidate(object):
  """A Candidate measures if a memory layout is suitable
     for a build configuration. If a build configuration
     shares an attribute with a Candidate, the Candidate
     gets a higher score.
     The scoring system is chosen such that the best
     Candidate is unique (i.e. 2**Importance).
  """

  IMPORTANCE = {'precision': 1, 'equations': 2, 'order': 3, 'pe': 4, 'multipleSimulations': 5}

  def __init__(self, atts):
    self.atts = atts

  def score(self, reqs):
    sc = 0
    for key,val in reqs.items():
      if key in self.atts and val == self.atts[key]:
        sc += 2**self.IMPORTANCE[key]
    return sc

  def __repr__(self):
    return repr(self.atts)

def findCandidates(search_path):
  """Determine Candidate attributes from file name."""

  archs = arch.getArchitectures()
  pes = [arch.getCpu(a) for a in archs]

  candidates = dict()
  for c in os.listdir(search_path):
    name, ext = os.path.splitext(c)
    atts = dict()
    for att in name.split('_'):
      multipleSimulations = re.match('ms([0-9]+)', att)
      order = re.match('O([0-9]+)', att)
      if multipleSimulations:
        atts['multipleSimulations'] = int(multipleSimulations.group(1))
      elif order:
        atts['order'] = int(order.group(1))
      elif att.lower() in ['s', 'd']:
        atts['precision'] = att.lower()
      elif att.lower() in pes:
        atts['pe'] = att.lower()
      else:
        atts['equations'] = att
    candidates[c] = Candidate(atts)
  return candidates

def guessMemoryLayout(env):    
  script_dir = os.path.dirname(os.path.abspath(__file__))
  path = os.path.join(script_dir, '..', 'auto_tuning', 'config')

  # from least to most
  importance = ['precision', 'equations', 'order', 'pe', 'multipleSimulations']
  values = {
    'precision': env['arch'][0].lower(),
    'equations': env['equations'].lower(),
    'order': int(env['order']),
    'pe': arch.getCpu(env['arch']),
    'multipleSimulations': int(env['multipleSimulations'])
  }

  candidates = findCandidates(search_path=path)
  bestFit = max(candidates.keys(), key=lambda key: candidates[key].score(values))
  bestScore = candidates[bestFit].score(values)

  if bestScore == 0:
    print('WARNING: No suitable memory layout found. (Will fall back to all dense.)')
    bestFit = 'dense.xml'
  print('Using memory layout {}'.format(bestFit))
  return os.path.join(path, bestFit)
