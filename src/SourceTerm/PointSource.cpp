/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Wolf (wolf.sebastian AT tum.de, https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2015 - 2020, SeisSol Group
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
#include <generated_code/kernel.h>
#include <generated_code/init.h>
#include <iostream>
 
void seissol::sourceterm::transformMomentTensor(real const i_localMomentTensor[3][3],
                                                real const i_localVelocityComponent[3],
                                                real strike,
                                                real dip,
                                                real rake,
                                                real o_forceComponents[seissol::sourceterm::PointSources::TensorSize])
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
  
  real M[3][3] = {
		  {0.0, 0.0, 0.0},
		  {0.0, 0.0, 0.0},
		  {0.0, 0.0, 0.0}
  };

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
  real f[3] = {0.0, 0.0, 0.0};
  for (unsigned j = 0; j < 3; ++j) {
    for (unsigned k = 0; k < 3; ++k) {
        f[k] += R[k][j] * i_localVelocityComponent[j];
    }
  }
  
#if NUMBER_OF_QUANTITIES < 6
  #error You cannot use PointSource with less than 6 quantities.
#endif
  
  // Save in order (\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \sigma_{xy}, \sigma_{yz}, \sigma_{xz}, u, v, w)
  o_forceComponents[0] = M[0][0];
  o_forceComponents[1] = M[1][1];
  o_forceComponents[2] = M[2][2];
  o_forceComponents[3] = M[0][1];
  o_forceComponents[4] = M[1][2];
  o_forceComponents[5] = M[0][2];
  o_forceComponents[6] = f[0];
  o_forceComponents[7] = f[1];
  o_forceComponents[8] = f[2];

  for (unsigned m = 9; m < seissol::sourceterm::PointSources::TensorSize; ++m) {
    o_forceComponents[m] = 0.0;
  }
}

real seissol::sourceterm::computePwLFTimeIntegral(PiecewiseLinearFunction1D const& i_pwLF,
                                               double i_fromTime,
                                               double i_toTime)
{
   real l_integral;
   // j_{from} := \argmax_j s.t. t_{from} >= t_{onset} + j*dt   =   floor[(t_{from} - t_{onset}) / dt]
   int l_fromIndex = (i_fromTime - i_pwLF.onsetTime) / i_pwLF.samplingInterval;
   // j_{to}   := \argmin_j s.t. t_{to}   >= t_{onset} + j*dt   =   floor[(t_{to} - t_{onset}) / dt]
   int l_toIndex = (i_toTime - i_pwLF.onsetTime) / i_pwLF.samplingInterval;
   
   l_fromIndex = std::max(0, l_fromIndex);
   l_toIndex = std::min(static_cast<int>(i_pwLF.numberOfPieces)-1, l_toIndex);
   
  /* The indefinite integral of the j-th linear function is
   * int m_j * t + n_j dt = 1 / 2 * m_j * t^2 + n_j * t
   */
   real l_time = i_pwLF.onsetTime + l_fromIndex * i_pwLF.samplingInterval;
   l_integral = 0.0;
   for (int j = l_fromIndex; j <= l_toIndex; ++j) {
     real tFrom = std::max((real)i_fromTime, l_time);
     l_time += i_pwLF.samplingInterval;
     real tTo = std::min((real)i_toTime, l_time);
     l_integral += 0.5 * i_pwLF.slopes[j] * (tTo * tTo - tFrom * tFrom) + i_pwLF.intercepts[j] * (tTo - tFrom);
   }
   
   return l_integral;
}

void seissol::sourceterm::addTimeIntegratedPointSourceNRF( real const i_mInvJInvPhisAtSources[tensor::mInvJInvPhisAtSources::size()],
                                                           real const faultBasis[9],
                                                           real A,
                                                           std::array<real, 81> const &stiffnessTensor,
                                                           std::array<PiecewiseLinearFunction1D, 3> const &slipRates,
                                                           double i_fromTime,
                                                           double i_toTime,
                                                           real o_dofUpdate[tensor::Q::size()],
                                                           unsigned int sourceNumber)
{  
  real slip[] = { 0.0, 0.0, 0.0};
  for (unsigned i = 0; i < 3; ++i) {
    if (slipRates[i].numberOfPieces > 0) {
      slip[i] = computePwLFTimeIntegral(slipRates[i], i_fromTime, i_toTime);
    }
  }
  
  real rotatedSlip[] = { 0.0, 0.0, 0.0 };
  for (unsigned i = 0; i < 3; ++i) {
    for (unsigned j = 0; j < 3; ++j) {
      rotatedSlip[j] += faultBasis[j + i*3] * slip[i];
    }
  }
  
  kernel::sourceNRF krnl;
  krnl.Q = o_dofUpdate;
  krnl.mInvJInvPhisAtSources = i_mInvJInvPhisAtSources;
  krnl.stiffnessTensor = stiffnessTensor.data();
  krnl.mSlip = rotatedSlip;
  krnl.mNormal = faultBasis + 6;
  krnl.mArea = -A;
  krnl.momentToNRF = init::momentToNRF::Values;
#ifdef MULTIPLE_SIMULATIONS
  std::array<double, MULTIPLE_SIMULATIONS> sourceToMultSim{};
  sourceToMultSim[sourceNumber % MULTIPLE_SIMULATIONS] = 1.0;
  krnl.oneSimToMultSim = sourceToMultSim.data();
#endif
  krnl.execute();
}

void seissol::sourceterm::addTimeIntegratedPointSourceFSRM( real const i_mInvJInvPhisAtSources[tensor::mInvJInvPhisAtSources::size()],
                                                            real const i_forceComponents[tensor::momentFSRM::size()],
                                                            PiecewiseLinearFunction1D const& i_pwLF,
                                                            double i_fromTime,
                                                            double i_toTime,
                                                            real o_dofUpdate[tensor::Q::size()],
                                                            unsigned int sourceNumber)
{
  kernel::sourceFSRM krnl;
  krnl.Q = o_dofUpdate;
  krnl.mInvJInvPhisAtSources = i_mInvJInvPhisAtSources;
  krnl.momentFSRM = i_forceComponents;
  krnl.stfIntegral = computePwLFTimeIntegral(i_pwLF, i_fromTime, i_toTime);
#ifdef MULTIPLE_SIMULATIONS
  std::array<double, MULTIPLE_SIMULATIONS> sourceToMultSim{};
  sourceToMultSim[sourceNumber % MULTIPLE_SIMULATIONS] = 1.0;
  krnl.oneSimToMultSim = sourceToMultSim.data();
#endif
  krnl.execute();
}
