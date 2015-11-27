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

def cardinalityDS(blocks):
  return sum([block.rows()*block.cols() for block in blocks])
  
def intersect(block1, block2):
  return min(block1.stoprow, block2.stoprow) > max(block1.startrow, block2.startrow) and min(block1.stopcol, block2.stopcol) > max(block1.startcol, block2.startcol)

def pairwiseDisjoint(blocks):
  n = len(blocks)
  for i in range(1, n):
    for j in range(i):
      if intersect(blocks[j], blocks[i]):
        return False
  return True

def maxDisjointSet(blocks, targetCard):
  n = len(blocks)
  if n > 24:
    raise NotImplementedError("This implementation is the exact solution of an NP problem. It sucks pretty much and should not be used for large len(blocks).")
  maxnum = 0
  maxcard = 0
  for i in range(1 << n):
    subset = [blocks[j] for j in range(n) if (i & (1 << j))]
    if pairwiseDisjoint(subset):
      card = cardinalityDS(subset)
      if card == targetCard:
        maxnum = i
        break
      if card > maxcard:
        maxcard = card
        maxnum = i
  return [j for j in range(n) if (maxnum & (1 << j))]
