#!/usr/bin/env python3
##
# @file
# This file is part of SeisSol.
#
# @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
#
# @section LICENSE
# Copyright (c) 2016-2019, SeisSol Group
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
  
import numpy as np
from yateto import Tensor
from yateto.input import parseXMLMatrixFile, memoryLayoutFromFile

from aderdg import LinearADERDG

class ViscoelasticADERDG(LinearADERDG):
  def __init__(self, order, multipleSimulations, matricesDir, memLayout, numberOfMechanisms, **kwargs):
    self.numberOfMechanisms = numberOfMechanisms
    self.numberOfElasticQuantities = 9

    super().__init__(order, multipleSimulations, matricesDir)
    clones = {
      'star': ['star(0)', 'star(1)', 'star(2)'],
    }
    self.db.update( parseXMLMatrixFile('{}/matrices_viscoelastic.xml'.format(matricesDir), clones) )

    star_spp = self.db.star[0].spp().as_ndarray()
    star_rows, star_cols = star_spp.shape
    aniso_cols = star_cols - self.numberOfElasticQuantities
    star_spp_new = np.zeros((self.numberOfQuantities(), self.numberOfQuantities()), dtype=bool)
    star_spp_new[0:star_rows,0:star_cols] = star_spp
    ''' The last 6 columns of star_spp contain the prototype sparsity pattern for
        a mechanism. Therefore, the spp is repeated for every mechanism. '''
    for mech in range(1,numberOfMechanisms):
      offset0 = self.numberOfElasticQuantities
      offsetm = self.numberOfElasticQuantities + mech*aniso_cols
      star_spp_new[0:star_rows,offsetm:offsetm+aniso_cols] = star_spp[0:star_rows,offset0:offset0+aniso_cols]
    for dim in range(3):
      self.db.star[dim] = Tensor(self.db.star[dim].name(), star_spp_new.shape, spp=star_spp_new)

    source_spp = np.zeros((self.numberOfQuantities(), self.numberOfQuantities()), dtype=bool)
    ET_spp = self.db['ET'].spp().as_ndarray()
    ''' ET is a prototype sparsity pattern for a mechanism. Therefore, repeated for every
        mechanism. See Kaeser and Dumbser 2006, III. Viscoelastic attenuation.
    '''
    for mech in range(numberOfMechanisms):
      offset = self.numberOfElasticQuantities + mech*aniso_cols
      r = slice(offset, offset+aniso_cols)
      source_spp[r,0:aniso_cols] = ET_spp
      source_spp[r,r] = np.identity(aniso_cols, dtype=bool)
    self.db.ET = Tensor('ET', source_spp.shape, spp=source_spp)

    memoryLayoutFromFile(memLayout, self.db, clones)

    self.kwargs = kwargs

  def numberOfQuantities(self):
    return 9 + 6*self.numberOfMechanisms

  def starMatrix(self, dim):
    return self.db.star[dim]

  def sourceMatrix(self):
    return self.db.ET

  def godunov_spp(self):
    spp = np.zeros((self.numberOfQuantities(), self.numberOfQuantities()), dtype=bool)
    spp[0:self.numberOfElasticQuantities,:] = True
    return spp

  def flux_solver_spp(self):
    return self.godunov_spp()

  def transformation_spp(self):
    spp = np.zeros((self.numberOfQuantities(), self.numberOfQuantities()), dtype=bool)
    spp[0:6,0:6] = 1
    spp[6:9,6:9] = 1
    for mechs in range(self.numberOfMechanisms):
      offset = 9 + mechs*6
      spp[offset:offset+6,offset:offset+6] = 1
    return spp

  def transformation_inv_spp(self):
    return self.transformation_spp()
