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

#include <generated_code/dr_kernels.h>
#include <Kernels/common.hpp>
#include <Kernels/denseMatrixOps.hpp>

/// \todo Make this information available in the code generator
/// \todo FIXME Will not work for viscoelasticity
seissol::kernels::DynamicRupture::DynamicRupture() {
  m_derivativesOffsets[0] = 0;
  for( int l_order = 0; l_order < CONVERGENCE_ORDER; l_order++ ) {
    m_numberOfAlignedBasisFunctions[l_order] = getNumberOfAlignedBasisFunctions( CONVERGENCE_ORDER-l_order, ALIGNMENT );

    if( l_order > 0 ) {
      m_derivativesOffsets[l_order]  =  m_numberOfAlignedBasisFunctions[l_order-1] * NUMBER_OF_QUANTITIES;
      m_derivativesOffsets[l_order] +=  m_derivativesOffsets[l_order-1];
    }
  }
}

void seissol::kernels::DynamicRupture::setTimeStepWidth(double timestep) {
  double subIntervalWidth = timestep / CONVERGENCE_ORDER;
  for (unsigned timeInterval = 0; timeInterval < CONVERGENCE_ORDER; ++timeInterval) {
    double t1 = timeInterval * subIntervalWidth;
    double t2 = t1 + subIntervalWidth;
    /// Compute time-integrated Taylor expansion (at t0=0) weights for interval [t1,t2].
    unsigned factorial = 1;
    for (unsigned derivative = 0; derivative < CONVERGENCE_ORDER; ++derivative) {
      m_timeAverageFactors[timeInterval][derivative] = (t2-t1) / (factorial * subIntervalWidth);
      t1 *= t1;
      t2 *= t2;
      factorial *= (derivative+2);
    }
  }
}

void seissol::kernels::DynamicRupture::computeTimeAverage( unsigned timeInterval,
                                                           real const* timeDerivatives,
                                                           real timeAverage[NUMBER_OF_ALIGNED_DOFS] ) {
  SXt(  m_timeAverageFactors[timeInterval][0],
        m_numberOfAlignedBasisFunctions[0],
        NUMBER_OF_QUANTITIES,
        timeDerivatives,
        m_numberOfAlignedBasisFunctions[0],
        timeAverage,
        NUMBER_OF_ALIGNED_BASIS_FUNCTIONS );

  for (unsigned derivative = 1; derivative < CONVERGENCE_ORDER; ++derivative) {
    SXtYp(  m_timeAverageFactors[timeInterval][derivative],
            m_numberOfAlignedBasisFunctions[derivative],
            NUMBER_OF_QUANTITIES,
            timeDerivatives + m_derivativesOffsets[derivative],
            m_numberOfAlignedBasisFunctions[derivative],
            timeAverage,
            NUMBER_OF_ALIGNED_BASIS_FUNCTIONS );
  }
}

void seissol::kernels::DynamicRupture::computeGodunovState( DRFaceInformation const&    faceInfo,
                                                            GlobalData const*           global,
                                                            DRGodunovData const*        godunovData,
                                                            real const*                 timeDerivativePlus,
                                                            real const*                 timeDerivativeMinus,
                                                            real                        godunov[CONVERGENCE_ORDER][seissol::model::godunovState::reals] ) {
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

  memset(godunov, 0, CONVERGENCE_ORDER * seissol::model::godunovState::reals * sizeof(real));
  
  for (unsigned timeInterval = 0; timeInterval < CONVERGENCE_ORDER; ++timeInterval) {
    real timeAveragePlus[NUMBER_OF_ALIGNED_DOFS] __attribute__((aligned(PAGESIZE_STACK)));
    real timeAverageMinus[NUMBER_OF_ALIGNED_DOFS] __attribute__((aligned(PAGESIZE_STACK)));
    computeTimeAverage(timeInterval, timeDerivativePlus, timeAveragePlus);
    computeTimeAverage(timeInterval, timeDerivativeMinus, timeAverageMinus);
    
    seissol::generatedKernels::godunovState[4*faceInfo.plusSide](
      godunovData->godunovMatrixPlus,
      global->faceToNodalMatrices[faceInfo.plusSide][0],
      timeAveragePlus,
      &godunov[timeInterval][0]
    );
  
    seissol::generatedKernels::godunovState[4*faceInfo.minusSide + faceInfo.faceRelation](
      godunovData->godunovMatrixMinus,
      global->faceToNodalMatrices[faceInfo.minusSide][faceInfo.faceRelation],
      timeAverageMinus,
      &godunov[timeInterval][0]
    ); 
  }
}
