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
#include <generated_code/flops.h>

unsigned seissol::kernels::Plasticity::computePlasticity( double                      relaxTime,
                                                      double                      timeStepWidth,
                                                      GlobalData const*           global,
                                                      PlasticityData const*       plasticityData,
                                                      real                        degreesOfFreedom[ NUMBER_OF_ALIGNED_DOFS ],
                                                      double*                     pstrain)
{
  real interpolationDofs[seissol::model::interpolationDOFS::reals] __attribute__((aligned(ALIGNMENT)));
  real meanStress[seissol::model::interpolationDOFS::ld] __attribute__((aligned(ALIGNMENT)));
  real tau[seissol::model::interpolationDOFS::ld] __attribute__((aligned(ALIGNMENT)));
  real taulim[seissol::model::interpolationDOFS::ld] __attribute__((aligned(ALIGNMENT)));
  real yieldFactor[seissol::model::interpolationDOFS::ld] __attribute__((aligned(ALIGNMENT)));
  real dudt_pstrain[7];
  
  seissol::generatedKernels::evaluateAtNodes(degreesOfFreedom, global->vandermondeMatrix, interpolationDofs);
  
  //copy dofs for later comparison, only first dof of stresses required
  real prev_degreesOfFreedom[6];
  for (unsigned q = 0; q < 6; ++q) {
	  prev_degreesOfFreedom[q] = degreesOfFreedom[q * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];
  }

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

    // calculate plastic strain with first dof only (for now)
    for (unsigned q = 0; q < 6; ++q) {
        real mufactor = plasticityData->mufactor;
        dudt_pstrain[q] = mufactor*(prev_degreesOfFreedom[q] - degreesOfFreedom[q * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS]);
        pstrain[q] += dudt_pstrain[q];
    }

    //accumulated plastic strain
    pstrain[6] += timeStepWidth * sqrt(0.5 * (dudt_pstrain[0]*dudt_pstrain[0] + dudt_pstrain[1]*dudt_pstrain[1]
            + dudt_pstrain[2]*dudt_pstrain[2])+ dudt_pstrain[3]*dudt_pstrain[3]
			+ dudt_pstrain[4]*dudt_pstrain[4] + dudt_pstrain[5]*dudt_pstrain[5]);
      
    return 1;
  }
  
  return 0;
}

void seissol::kernels::Plasticity::flopsPlasticity( long long&  o_nonZeroFlopsCheck,
                                                    long long&  o_hardwareFlopsCheck,
                                                    long long&  o_nonZeroFlopsYield,
                                                    long long&  o_hardwareFlopsYield )
{
  // reset flops
  o_nonZeroFlopsCheck = 0; o_hardwareFlopsCheck = 0;
  o_nonZeroFlopsYield = 0; o_hardwareFlopsYield = 0;
  
  // flops from checking, i.e. outside if (adjust) {}
  o_nonZeroFlopsCheck  += seissol::flops::evaluateAtNodes_nonZero;
  o_hardwareFlopsCheck += seissol::flops::evaluateAtNodes_hardware;
  
  // add initial loading
  o_nonZeroFlopsCheck  += 6 * seissol::model::interpolationDOFS::rows;
  o_hardwareFlopsCheck += 6 * seissol::model::interpolationDOFS::ld;
  
  // compute mean stress (2 adds, 1 mul)
  o_nonZeroFlopsCheck  += 3 * seissol::model::interpolationDOFS::rows;
  o_hardwareFlopsCheck += 3 * seissol::model::interpolationDOFS::ld;
  
  // subtract mean stress
  o_nonZeroFlopsCheck  += 3 * seissol::model::interpolationDOFS::rows;
  o_hardwareFlopsCheck += 3 * seissol::model::interpolationDOFS::ld;
  
  // compute tau (5 adds, 7 muls, sqrt NOT counted)
  o_nonZeroFlopsCheck  += 12 * seissol::model::interpolationDOFS::rows;
  o_hardwareFlopsCheck += 12 * seissol::model::interpolationDOFS::ld;
  
  // compute taulim (1 add, 1 mul, max NOT counted)
  o_nonZeroFlopsCheck  += 2 * seissol::model::interpolationDOFS::rows;
  o_hardwareFlopsCheck += 2 * seissol::model::interpolationDOFS::ld;
  
  // check for yield (NOT counted, as it would require counting the number of yielding points) 

  // flops from plastic yielding, i.e. inside if (adjust) {}
  o_nonZeroFlopsYield  += seissol::flops::convertToModal_nonZero;
  o_hardwareFlopsYield += seissol::flops::convertToModal_hardware;
    
  // adjust 3 bulk stresses (2 adds, 1 mul)
  o_nonZeroFlopsYield  += 3 * 3 * seissol::model::interpolationDOFS::rows;
  o_hardwareFlopsYield += 3 * 3 * seissol::model::interpolationDOFS::ld;

  // adjust 3 shear stresses (1 adds, 1 mul)
  o_nonZeroFlopsYield  += 3 * 2 * seissol::model::interpolationDOFS::rows;
  o_hardwareFlopsYield += 3 * 2 * seissol::model::interpolationDOFS::ld;
}
