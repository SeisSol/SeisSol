/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Wolf (wolf.sebastian AT tum.de,
 *https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2015 - 2020, SeisSol Group
 * Copyright (c) 2023, Intel corporation
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Point source computation.
 **/

#include "PointSource.h"
#include <cmath>
#include <algorithm>

void seissol::sourceterm::transformMomentTensor(
    real const i_localMomentTensor[3][3],
    real const i_localSolidVelocityComponent[3],
    real i_localPressureComponent,
    real const i_localFluidVelocityComponent[3],
    real strike,
    real dip,
    real rake,
    AlignedArray<real, PointSources::TensorSize>& o_forceComponents) {
  real cstrike = cos(strike);
  real sstrike = sin(strike);
  real cdip = cos(dip);
  real sdip = sin(dip);
  real crake = cos(rake);
  real srake = sin(rake);

  // Note, that R[j][i] = R_{ij} here.
  real R[3][3] = {{crake * cstrike + cdip * srake * sstrike,
                   cdip * crake * sstrike - cstrike * srake,
                   sdip * sstrike},
                  {cdip * cstrike * srake - crake * sstrike,
                   srake * sstrike + cdip * crake * cstrike,
                   cstrike * sdip},
                  {-sdip * srake, -crake * sdip, cdip}};

  real M[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

  // Calculate M_{ij} = R_{ki} * LM_{kl} * R_{lj}.
  // Note, again, that X[j][i] = X_{ij} here.
  // As M is symmetric, it is sufficient to calculate
  // (i,j) = (0,0), (1,0), (2,0), (1,1), (2,1), (2,2)
  for (unsigned j = 0; j < 3; ++j) {
    for (unsigned i = j; i < 3; ++i) {
      for (unsigned k = 0; k < 3; ++k) {
        for (unsigned l = 0; l < 3; ++l) {
          M[j][i] += R[i][k] * i_localMomentTensor[l][k] * R[j][l];
        }
      }
    }
  }
  real f[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  for (unsigned j = 0; j < 3; ++j) {
    for (unsigned k = 0; k < 3; ++k) {
      f[k] += R[k][j] * i_localSolidVelocityComponent[j];
      f[k + 3] += R[k][j] * i_localFluidVelocityComponent[j];
    }
  }

#if NUMBER_OF_QUANTITIES < 6
#error You cannot use PointSource with less than 6 quantities.
#endif

  std::fill(o_forceComponents.data(), o_forceComponents.data() + o_forceComponents.size(), 0);
  // Save in order (\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \sigma_{xy}, \sigma_{yz}, \sigma_{xz}, u,
  // v, w, p, u_f, v_f, w_f)
  o_forceComponents[0] = M[0][0];
  o_forceComponents[1] = M[1][1];
  o_forceComponents[2] = M[2][2];
  o_forceComponents[3] = M[0][1];
  o_forceComponents[4] = M[1][2];
  o_forceComponents[5] = M[0][2];
  o_forceComponents[6] = f[0];
  o_forceComponents[7] = f[1];
  o_forceComponents[8] = f[2];
#ifdef USE_POROELASTIC
  o_forceComponents[9] = i_localPressureComponent;
  o_forceComponents[10] = f[3];
  o_forceComponents[11] = f[4];
  o_forceComponents[12] = f[5];
#endif
}
