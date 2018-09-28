/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2016, SeisSol Group
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
 * Dynamic Rupture kernel of SeisSol.
 **/

#include "DynamicRupture.h"

#ifndef NDEBUG
#pragma message "compiling dynamic rupture kernel with assertions"
#endif

#include <cassert>
#include <cstring>
#include <stdint.h>

#include <generated_code/kernel.h>
#include <Kernels/common.hpp>
#include <Kernels/denseMatrixOps.hpp>
#include <Numerical_aux/Quadrature.h>
#include <yateto.h>

seissol::kernels::DynamicRupture::DynamicRupture() {
  m_derivativesOffsets[0] = 0;
  for (int order = 0; order < CONVERGENCE_ORDER; ++order) {
    if (order > 0) {
      m_derivativesOffsets[order] = tensor::dQ::size(order-1) + m_derivativesOffsets[order-1];
    }
  }
}

void seissol::kernels::DynamicRupture::setTimeStepWidth(double timestep)
{
#ifdef USE_DR_CELLAVERAGE
  double subIntervalWidth = timestep / CONVERGENCE_ORDER;
  for (unsigned timeInterval = 0; timeInterval < CONVERGENCE_ORDER; ++timeInterval) {
    double t1 = timeInterval * subIntervalWidth;
    double t2 = t1 + subIntervalWidth;
    /// Compute time-integrated Taylor expansion (at t0=0) weights for interval [t1,t2].
    unsigned factorial = 1;
    for (unsigned derivative = 0; derivative < CONVERGENCE_ORDER; ++derivative) {
      m_timeFactors[timeInterval][derivative] = (t2-t1) / (factorial * subIntervalWidth);
      t1 *= t1;
      t2 *= t2;
      factorial *= (derivative+2);
    }
    /// We define the time "point" of the interval as the centre of the interval in order
    /// to be somewhat compatible to legacy code.
    timePoints[timeInterval] = timeInterval * subIntervalWidth + subIntervalWidth / 2.;
    timeWeights[timeInterval] = subIntervalWidth;
  }
#else
  seissol::quadrature::GaussLegendre(timePoints, timeWeights, CONVERGENCE_ORDER);
  for (unsigned point = 0; point < CONVERGENCE_ORDER; ++point) {
    timePoints[point] = 0.5 * (timestep * timePoints[point] + timestep);
    timeWeights[point] = 0.5 * timestep * timeWeights[point];
  }

  for (unsigned point = 0; point < CONVERGENCE_ORDER; ++point) {
    double time = 1.0;
    unsigned factorial = 1;
    for (unsigned derivative = 0; derivative < CONVERGENCE_ORDER; ++derivative) {
      m_timeFactors[point][derivative] = time / factorial;
      time *= timePoints[point];
      factorial *= (derivative+1);
    }
  }
#endif
}

void seissol::kernels::DynamicRupture::computeGodunovState( DRFaceInformation const&    faceInfo,
                                                            GlobalData const*           global,
                                                            DRGodunovData const*        godunovData,
                                                            real const*                 timeDerivativePlus,
                                                            real const*                 timeDerivativeMinus,
                                                            real                        godunov[CONVERGENCE_ORDER][seissol::tensor::godunovState::size()],
                                                            real const*                 timeDerivativePlus_prefetch,
                                                            real const*                 timeDerivativeMinus_prefetch ) {
  // assert alignments
#ifndef NDEBUG
  for (unsigned face = 0; face < 4; ++face) {
    for (unsigned h = 0; h < 4; ++h) {
      assert( ((uintptr_t)global->faceToNodalMatrices[face][h]) % ALIGNMENT == 0 );
    }
  }
  assert( ((uintptr_t)timeDerivativePlus) % ALIGNMENT == 0 );
  assert( ((uintptr_t)timeDerivativeMinus) % ALIGNMENT == 0 );
  assert( ((uintptr_t)&godunov[0])         % ALIGNMENT == 0 );
#endif

  real degreesOfFreedomPlus[NUMBER_OF_ALIGNED_DOFS] __attribute__((aligned(PAGESIZE_STACK)));
  real degreesOfFreedomMinus[NUMBER_OF_ALIGNED_DOFS] __attribute__((aligned(PAGESIZE_STACK)));

  kernel::derivativeTaylorExpansion taylorKrnlPlus;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    taylorKrnlPlus.dQ(i) = timeDerivativePlus + m_derivativesOffsets[i];
  }
  taylorKrnlPlus.I = degreesOfFreedomPlus;

  kernel::derivativeTaylorExpansion taylorKrnlMinus;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    taylorKrnlMinus.dQ(i) = timeDerivativeMinus + m_derivativesOffsets[i];
  }
  taylorKrnlMinus.I = degreesOfFreedomMinus;

  kernel::godunovState krnl;
  krnl.V3mTo2n = global->faceToNodalMatrices;

  for (unsigned timeInterval = 0; timeInterval < CONVERGENCE_ORDER; ++timeInterval) {    
    for (unsigned der = 0; der < yateto::numFamilyMembers<tensor::dQ>(); ++der) {
      taylorKrnlPlus.power = m_timeFactors[timeInterval][der];
      taylorKrnlPlus.execute(der);
    }
    for (unsigned der = 0; der < yateto::numFamilyMembers<tensor::dQ>(); ++der) {
      taylorKrnlMinus.power = m_timeFactors[timeInterval][der];
      taylorKrnlMinus.execute(der);
    }

    real const* plusPrefetch = (timeInterval < CONVERGENCE_ORDER-1) ? &godunov[timeInterval+1][0] : timeDerivativePlus_prefetch;
    real const* minusPrefetch = (timeInterval < CONVERGENCE_ORDER-1) ? &godunov[timeInterval+1][0] : timeDerivativeMinus_prefetch;
    
    krnl.godunovState = &godunov[timeInterval][0];
    
    krnl.Q = degreesOfFreedomPlus;
    krnl.godunovMatrix = godunovData->godunovMatrixPlus;
    krnl._prefetch.godunovState = plusPrefetch;
    krnl.execute(faceInfo.plusSide, 0);
    
    krnl.Q = degreesOfFreedomMinus;
    krnl.godunovMatrix = godunovData->godunovMatrixMinus;
    krnl._prefetch.godunovState = minusPrefetch;
    krnl.execute(faceInfo.minusSide, faceInfo.faceRelation);
  }
}

void seissol::kernels::DynamicRupture::flopsGodunovState( DRFaceInformation const&  faceInfo,
                                                          long long&                o_nonZeroFlops,
                                                          long long&                o_hardwareFlops )
{
  // reset flops
  o_nonZeroFlops = 0; o_hardwareFlops = 0;

  // evaluateTaylorExpansion flops
  for (unsigned der = 0; der < CONVERGENCE_ORDER; ++der) {
    o_nonZeroFlops += kernel::derivativeTaylorExpansion::nonZeroFlops(der);
    o_hardwareFlops += kernel::derivativeTaylorExpansion::hardwareFlops(der);
  }
  
  // 2x evaluateTaylorExpansion
  o_nonZeroFlops *= 2;
  o_hardwareFlops *= 2;
  
  o_nonZeroFlops += kernel::godunovState::nonZeroFlops(faceInfo.plusSide, 0);
  o_hardwareFlops += kernel::godunovState::hardwareFlops(faceInfo.plusSide, 0);
  
  o_nonZeroFlops += kernel::godunovState::nonZeroFlops(faceInfo.minusSide, faceInfo.faceRelation);
  o_hardwareFlops += kernel::godunovState::hardwareFlops(faceInfo.minusSide, faceInfo.faceRelation);
  
  o_nonZeroFlops *= CONVERGENCE_ORDER;
  o_hardwareFlops *= CONVERGENCE_ORDER;
}
