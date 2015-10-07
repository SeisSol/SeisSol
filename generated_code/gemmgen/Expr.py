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

import sympy

def checkExprTree(expr):
  valid = True
  if expr.is_MatAdd:
    for arg in expr.args:
      valid &= checkExprTree(arg)
  elif expr.is_MatMul:
    for arg in expr.args:
      valid &= isinstance(arg, sympy.MatrixSymbol)
  else:
    valid = False
  return valid

def analyseLeftAndRightMultiplication(db, expr):
  if expr.is_MatAdd:
    for arg in expr.args:
      analyseLeftAndRightMultiplication(db, arg)
  elif expr.is_MatMul:
    if len(expr.args) > 1:
      for arg in expr.args[0:-1]:
        db[arg.name].leftMultiplication = True
      db[expr.args[-1].name].rightMultiplication = True

def analyseKernels(db, kernels):
  for name, kernel in kernels:
    if checkExprTree(kernel.symbol) is False:
      raise ValueError('The {} kernel is not given in the canonical form \sum_{i=1}^s \prod_{j=1}^p(i) M_{ij}, where s >= 1 and p(i) >= 2 ({}).'.format(name, kernel.symbol))
    analyseLeftAndRightMultiplication(db, kernel.symbol)
