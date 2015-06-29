/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2015, SeisSol Group
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
 
void seissol::physics::transformMomentTensor(real const i_localMomentTensor[3][3],
                                             real strike,
                                             real dip,
                                             real rake,
                                             real o_momentTensor[NUMBER_OF_QUANTITIES])
{
  real cstrike = cos(strike);
  real sstrike = sin(strike);
  real cdip = cos(dip);
  real sdip = sin(dip);
  real crake = cos(rake);
  real srake = sin(rake);
  
  // Note, that R[j][i] = R_{ij} here.
  real R[3][3] = {
    { crake*cstrike + cdip*srake*sstrike, cdip*crake*sstrike - cstrike*srake, sdip*sstrike },
    { cdip*cstrike*srake - crake*sstrike, srake*sstrike + cdip*crake*cstrike, cstrike*sdip },
    {                        -sdip*srake,                        -crake*sdip,         cdip }
  };
  
  real M[3][3] = { 0.0 };

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
  
#if NUMBER_OF_QUANTITIES < 6
  #error You cannot use PointSource with less than 6 quantities.
#endif
  
  // Save in order (\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \sigma_{xy}, \sigma_{yz}, \sigma_{xz}, u, v, w)
  o_momentTensor[0] = M[0][0];
  o_momentTensor[1] = M[1][1];
  o_momentTensor[2] = M[2][2];
  o_momentTensor[3] = M[0][1];
  o_momentTensor[4] = M[1][2];
  o_momentTensor[5] = M[0][2];
  for (unsigned m = 6; m < NUMBER_OF_QUANTITIES; ++m) {
    o_momentTensor[m] = 0.0;
  }
}

real seissol::physics::computePwLFTimeIntegral(PiecewiseLinearFunction1D const* i_pwLF,
                                               double i_fromTime,
                                               double i_toTime)
{   
   real l_integral;
   // j_{from} := \argmax_j s.t. t_{from} >= t_{onset} + j*dt   =   floor[(t_{from} - t_{onset}) / dt]
   int l_fromIndex = (i_fromTime - i_pwLF->onsetTime) / i_pwLF->samplingInterval;
   // j_{to}   := \argmin_j s.t. t_{to}   >= t_{onset} + j*dt   =   floor[(t_{to} - t_{onset}) / dt]
   int l_toIndex = (i_toTime - i_pwLF->onsetTime) / i_pwLF->samplingInterval;
   
   l_fromIndex = std::max(0, l_fromIndex);
   l_toIndex = std::min(static_cast<int>(i_pwLF->numberOfPieces)-1, l_toIndex);
   
  /* The indefinite integral of the j-th linear function is
   * int m_j * t + n_j dt = 1 / 2 * m_j * t^2 + n_j * t
   */
   real l_time = i_pwLF->onsetTime + l_fromIndex * i_pwLF->samplingInterval;
   l_integral = 0.0;
   for (int j = l_fromIndex; j <= l_toIndex; ++j) {
     real tFrom = std::max(i_fromTime, l_time);
     l_time += i_pwLF->samplingInterval;
     real tTo = std::min(i_toTime, l_time);
     l_integral += 0.5 * i_pwLF->slopes[j] * (tTo * tTo - tFrom * tFrom) + i_pwLF->intercepts[j] * (tTo - tFrom);
   }
   
   return l_integral;
}

void seissol::physics::addTimeIntegratedPointSource(real const i_mInvJInvPhisAtSources[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
                                                    real const i_momentTensor[NUMBER_OF_QUANTITIES],
                                                    PiecewiseLinearFunction1D const* i_pwLF,
                                                    double i_fromTime,
                                                    double i_toTime,
                                                    real o_dofUpdate[NUMBER_OF_ALIGNED_DOFS])
{
  real l_integral = computePwLFTimeIntegral(i_pwLF, i_fromTime, i_toTime);  

  for (unsigned quantity = 0; quantity < NUMBER_OF_QUANTITIES; ++quantity) {
    for (unsigned basisFunction = 0; basisFunction < NUMBER_OF_ALIGNED_BASIS_FUNCTIONS; ++basisFunction) {
      o_dofUpdate[quantity * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + basisFunction] += l_integral * i_mInvJInvPhisAtSources[basisFunction] * i_momentTensor[quantity];
    }
  }
}
