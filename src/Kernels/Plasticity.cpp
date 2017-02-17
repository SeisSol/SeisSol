/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Stephanie Wollherr (wollherr AT geophysik.uni-muenchen.de, https://www.geophysik.uni-muenchen.de/Members/wollherr)
 *
 * @section LICENSE
 * Copyright (c) 2017, SeisSol Group
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
 * IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Plasticity kernel of SeisSol.
 **/

#include "Plasticity.h"

#include <cstring>
#include <algorithm>
#include <cmath>
#include <generated_code/kernels.h>

void seissol::kernels::Plasticity::computePlasticity( double                      relaxTime,
                                                      GlobalData const*           global,
                                                      PlasticityData const*       plasticityData,
                                                      real                        degreesOfFreedom[ NUMBER_OF_ALIGNED_DOFS ] )
{
  real interpolationDofs[seissol::model::interpolationDOFS::reals] __attribute__((aligned(ALIGNMENT)));
  real meanStress[seissol::model::interpolationDOFS::ld] __attribute__((aligned(ALIGNMENT)));
  real tau[seissol::model::interpolationDOFS::ld] __attribute__((aligned(ALIGNMENT)));
  real taulim[seissol::model::interpolationDOFS::ld] __attribute__((aligned(ALIGNMENT)));
  real yieldFactor[seissol::model::interpolationDOFS::ld] __attribute__((aligned(ALIGNMENT)));
  
  seissol::generatedKernels::evaluateAtNodes(degreesOfFreedom, global->vandermondeMatrix, interpolationDofs);
  
  for (unsigned q = 0; q < 6; ++q) {
    real initialLoading = plasticityData->initialLoading[q][0];
    for (unsigned ip = 0; ip < seissol::model::interpolationDOFS::ld; ++ip) {
      interpolationDofs[q * seissol::model::interpolationDOFS::ld + ip] += initialLoading;
    }
  }
  
  for (unsigned ip = 0; ip < seissol::model::interpolationDOFS::ld; ++ip) {
    meanStress[ip] = ( interpolationDofs[0 * seissol::model::interpolationDOFS::ld + ip]
                      + interpolationDofs[1 * seissol::model::interpolationDOFS::ld + ip]
                      + interpolationDofs[2 * seissol::model::interpolationDOFS::ld + ip] ) * (1.0 / 3.0);
  }

  for (unsigned q = 0; q < 3; ++q) {
    for (unsigned ip = 0; ip < seissol::model::interpolationDOFS::ld; ++ip) {
      interpolationDofs[q * seissol::model::interpolationDOFS::ld + ip] -= meanStress[ip];
    }
  }
  
  for (unsigned ip = 0; ip < seissol::model::interpolationDOFS::ld; ++ip) {
    tau[ip] = sqrt(0.5 * (interpolationDofs[0 * seissol::model::interpolationDOFS::ld + ip] * interpolationDofs[0 * seissol::model::interpolationDOFS::ld + ip]
            + interpolationDofs[1 * seissol::model::interpolationDOFS::ld + ip] * interpolationDofs[1 * seissol::model::interpolationDOFS::ld + ip]
            + interpolationDofs[2 * seissol::model::interpolationDOFS::ld + ip] * interpolationDofs[2 * seissol::model::interpolationDOFS::ld + ip])
            + interpolationDofs[3 * seissol::model::interpolationDOFS::ld + ip] * interpolationDofs[3 * seissol::model::interpolationDOFS::ld + ip]
            + interpolationDofs[4 * seissol::model::interpolationDOFS::ld + ip] * interpolationDofs[4 * seissol::model::interpolationDOFS::ld + ip]
            + interpolationDofs[5 * seissol::model::interpolationDOFS::ld + ip] * interpolationDofs[5 * seissol::model::interpolationDOFS::ld + ip]);
  }
  
  for (unsigned ip = 0; ip < seissol::model::interpolationDOFS::ld; ++ip) {
    taulim[ip] = std::max((real) 0.0, plasticityData->cohesionTimesCosAngularFriction - meanStress[ip] * plasticityData->sinAngularFriction);
  }
  
  bool adjust = false;
  for (unsigned ip = 0; ip < seissol::model::interpolationDOFS::rows; ++ip) {
    if (tau[ip] > taulim[ip]) {
      adjust = true;
      yieldFactor[ip] = 1.0 - (1.0 - taulim[ip] / tau[ip]) * relaxTime;
    } else {
      yieldFactor[ip] = 1.0;
    }
  }
  
  if (adjust) {
    for (unsigned q = 0; q < 3; ++q) {
      real initialLoading = plasticityData->initialLoading[q][0];
      for (unsigned ip = 0; ip < seissol::model::interpolationDOFS::ld; ++ip) {
        interpolationDofs[q * seissol::model::interpolationDOFS::ld + ip] = yieldFactor[ip] * interpolationDofs[q * seissol::model::interpolationDOFS::ld + ip]
                                                                          + meanStress[ip]
                                                                          - initialLoading;
      }
    }
    for (unsigned q = 3; q < 6; ++q) {
      real initialLoading = plasticityData->initialLoading[q][0];
      for (unsigned ip = 0; ip < seissol::model::interpolationDOFS::ld; ++ip) {
        interpolationDofs[q * seissol::model::interpolationDOFS::ld + ip] = yieldFactor[ip] * interpolationDofs[q * seissol::model::interpolationDOFS::ld + ip]
                                                                          - initialLoading;
      }
    }
    generatedKernels::convertToModal(interpolationDofs, global->vandermondeMatrixInverse, degreesOfFreedom);
  }
  
}
