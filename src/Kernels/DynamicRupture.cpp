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

#include <generated_code/kernels.h>
#include <generated_code/flops.h>
#include <Kernels/common.hpp>
#include <Kernels/denseMatrixOps.hpp>
#include <Numerical_aux/Quadrature.h>

seissol::kernels::DynamicRupture::DynamicRupture() {
  m_derivativesOffsets[0] = 0;
  for( int l_order = 0; l_order < CONVERGENCE_ORDER; l_order++ ) {

#ifdef NUMBER_OF_ALIGNED_DER_BASIS_FUNCTIONS /// Derivatives are compressed
    m_numberOfAlignedBasisFunctions[l_order] = getNumberOfAlignedBasisFunctions( CONVERGENCE_ORDER-l_order, ALIGNMENT );
#else /// Derivatives are NOT compressed
    m_numberOfAlignedBasisFunctions[l_order] = getNumberOfAlignedBasisFunctions( CONVERGENCE_ORDER, ALIGNMENT );
#endif

    if( l_order > 0 ) {
      m_derivativesOffsets[l_order]  =  m_numberOfAlignedBasisFunctions[l_order-1] * NUMBER_OF_QUANTITIES;
      m_derivativesOffsets[l_order] +=  m_derivativesOffsets[l_order-1];
    }
  }
  
  double points[NUMBER_OF_SPACE_QUADRATURE_POINTS][2];  
  seissol::quadrature::TriangleQuadrature(points, spaceWeights, CONVERGENCE_ORDER+1);
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

void seissol::kernels::DynamicRupture::evaluateTaylorExpansion( unsigned timeInterval,
                                                                real const* timeDerivatives,
                                                                real degreesOfFreedom[NUMBER_OF_ALIGNED_DOFS] ) {
  SXt(  m_timeFactors[timeInterval][0],
        m_numberOfAlignedBasisFunctions[0],
        NUMBER_OF_QUANTITIES,
        timeDerivatives,
        m_numberOfAlignedBasisFunctions[0],
        degreesOfFreedom,
        NUMBER_OF_ALIGNED_BASIS_FUNCTIONS );

  for (unsigned derivative = 1; derivative < CONVERGENCE_ORDER; ++derivative) {
    SXtYp(  m_timeFactors[timeInterval][derivative],
            m_numberOfAlignedBasisFunctions[derivative],
            NUMBER_OF_QUANTITIES,
            timeDerivatives + m_derivativesOffsets[derivative],
            m_numberOfAlignedBasisFunctions[derivative],
            degreesOfFreedom,
            NUMBER_OF_ALIGNED_BASIS_FUNCTIONS );
  }
}

void seissol::kernels::DynamicRupture::computeGodunovState( DRFaceInformation const&    faceInfo,
                                                            GlobalData const*           global,
                                                            DRGodunovData const&        godunovData,
                                                            real const*                 timeDerivativePlus,
                                                            real const*                 timeDerivativeMinus,
                                                            real                        godunov[CONVERGENCE_ORDER][seissol::model::godunovState::reals],
                                                            real const*                 timeDerivativePlus_prefetch,
                                                            real const*                 timeDerivativeMinus_prefetch,
                                                            DROutput&                   drOutput ) {
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

  for (unsigned timeInterval = 0; timeInterval < CONVERGENCE_ORDER; ++timeInterval) {
    real degreesOfFreedomPlus[NUMBER_OF_ALIGNED_DOFS] __attribute__((aligned(PAGESIZE_STACK)));
    real degreesOfFreedomMinus[NUMBER_OF_ALIGNED_DOFS] __attribute__((aligned(PAGESIZE_STACK)));
    evaluateTaylorExpansion(timeInterval, timeDerivativePlus, degreesOfFreedomPlus);
    evaluateTaylorExpansion(timeInterval, timeDerivativeMinus, degreesOfFreedomMinus);

    real const* plusPrefetch = (timeInterval < CONVERGENCE_ORDER-1) ? &godunov[timeInterval+1][0] : timeDerivativePlus_prefetch;
    real const* minusPrefetch = (timeInterval < CONVERGENCE_ORDER-1) ? &godunov[timeInterval+1][0] : timeDerivativeMinus_prefetch;        

    real atQPPointsPlus[seissol::model::godunovState::reals] __attribute__((aligned(PAGESIZE_STACK)));
    real atQPPointsMinus[seissol::model::godunovState::reals] __attribute__((aligned(PAGESIZE_STACK)));

    seissol::generatedKernels::evaluateAtQuadraturePoints[4*faceInfo.plusSide](
      global->faceToNodalMatrices[faceInfo.plusSide][0],
      degreesOfFreedomPlus,
      atQPPointsPlus,
      plusPrefetch
    );

    seissol::generatedKernels::evaluateAtQuadraturePoints[4*faceInfo.minusSide + faceInfo.faceRelation](
      global->faceToNodalMatrices[faceInfo.minusSide][faceInfo.faceRelation],
      degreesOfFreedomMinus,
      atQPPointsMinus,
      minusPrefetch
    );

    seissol::generatedKernels::godunovStatePlus(
      godunovData.godunovMatrixPlus,
      atQPPointsPlus,
      &godunov[timeInterval][0]
    );

    seissol::generatedKernels::godunovStateMinus(
      godunovData.godunovMatrixMinus,
      atQPPointsMinus,
      &godunov[timeInterval][0]
    );
    
    real tractionAndSlipRate[seissol::model::godunovState::ld*seissol::model::tractionAndSlipRateMatrix::cols] __attribute__((aligned(ALIGNMENT)));
    computeTractionAndSlipRate(godunovData, atQPPointsPlus, atQPPointsMinus, tractionAndSlipRate);    

    double norm[seissol::model::godunovState::ld] __attribute__((aligned(ALIGNMENT))) = {};
    for (unsigned d = 0; d < 3; ++d) {
      for (unsigned i = 0; i < seissol::model::godunovState::rows; ++i) {
        auto sr = tractionAndSlipRate[(3+d)*seissol::model::godunovState::ld + i];
        drOutput.slip[d*seissol::model::godunovState::ld + i] += timeWeights[timeInterval] * sr;
        norm[i] += sr*sr;
      }
    }
    for (unsigned i = 0; i < seissol::model::godunovState::rows; ++i) {
      drOutput.absoluteSlip[i] += timeWeights[timeInterval] * sqrt(norm[i]);
    }

    drOutput.frictionalEnergy += - timeWeights[timeInterval] * energySpaceIntegral(godunovData, tractionAndSlipRate);
  }
}

void seissol::kernels::DynamicRupture::computeTractionAndSlipRate(  DRGodunovData const&  godunovData,
                                                                    real const            atQPPointsPlus[seissol::model::godunovState::reals],
                                                                    real const            atQPPointsMinus[seissol::model::godunovState::reals],
                                                                    real                  tractionAndSlipRate[seissol::model::godunovState::ld*seissol::model::tractionAndSlipRateMatrix::cols] )
{
  real stressAndSlipRate[seissol::model::godunovState::reals] __attribute__((aligned(ALIGNMENT)));
  for (unsigned d = 0; d < 6; ++d) {
    for (unsigned i = 0; i < seissol::model::godunovState::rows; ++i) {
      stressAndSlipRate[d*seissol::model::godunovState::ld + i] = 0.5 * (atQPPointsMinus[d*seissol::model::godunovState::ld + i] + atQPPointsPlus[d*seissol::model::godunovState::ld + i]);
    }
  }
  for (unsigned d = 0; d < 3; ++d) {
    for (unsigned i = 0; i < seissol::model::godunovState::rows; ++i) {
      stressAndSlipRate[(6+d)*seissol::model::godunovState::ld + i] = atQPPointsMinus[(6+d)*seissol::model::godunovState::ld + i] - atQPPointsPlus[(6+d)*seissol::model::godunovState::ld + i];
    }
  }

  seissol::generatedKernels::computeTractionAndRotateSlipRate(stressAndSlipRate, godunovData.tractionAndSlipRateMatrix, tractionAndSlipRate);
}

double seissol::kernels::DynamicRupture::energySpaceIntegral( DRGodunovData const&  godunovData,
                                                              real                  tractionAndSlipRate[seissol::model::godunovState::ld*seissol::model::tractionAndSlipRateMatrix::cols] )
{
  double energy = 0.0;
  for (unsigned d = 0; d < 3; ++d) {
    for (unsigned i = 0; i < NUMBER_OF_SPACE_QUADRATURE_POINTS; ++i) {
      energy += spaceWeights[i] * tractionAndSlipRate[d*seissol::model::godunovState::ld + i] * tractionAndSlipRate[(3+d)*seissol::model::godunovState::ld + i];
    }
  }
  return godunovData.doubledSurfaceArea * energy;
}

void seissol::kernels::DynamicRupture::flopsGodunovState( DRFaceInformation const&  faceInfo,
                                                          long long&                o_nonZeroFlops,
                                                          long long&                o_hardwareFlops )
{
  // reset flops
  o_nonZeroFlops = 0; o_hardwareFlops = 0;
  
  // SXt flops
  o_nonZeroFlops += m_numberOfAlignedBasisFunctions[0] * NUMBER_OF_QUANTITIES;
  o_hardwareFlops += m_numberOfAlignedBasisFunctions[0] * NUMBER_OF_QUANTITIES;
  
  // SXtYp flops
  for (unsigned derivative = 1; derivative < CONVERGENCE_ORDER; ++derivative) {
    o_nonZeroFlops += 2 * m_numberOfAlignedBasisFunctions[derivative] * NUMBER_OF_QUANTITIES;
    o_hardwareFlops += 2 * m_numberOfAlignedBasisFunctions[derivative] * NUMBER_OF_QUANTITIES;
  }
  
  // 2x evaluateTaylorExpansion
  o_nonZeroFlops *= 2;
  o_hardwareFlops *= 2;

  o_nonZeroFlops += seissol::flops::evaluateAtQuadraturePoints_nonZero[4*faceInfo.plusSide];
  o_hardwareFlops += seissol::flops::evaluateAtQuadraturePoints_hardware[4*faceInfo.plusSide];
  
  o_nonZeroFlops += seissol::flops::evaluateAtQuadraturePoints_nonZero[4*faceInfo.minusSide + faceInfo.faceRelation];
  o_hardwareFlops += seissol::flops::evaluateAtQuadraturePoints_hardware[4*faceInfo.minusSide + faceInfo.faceRelation];
  
  o_nonZeroFlops += seissol::flops::godunovStatePlus_nonZero;
  o_hardwareFlops += seissol::flops::godunovStatePlus_hardware;

  o_nonZeroFlops += seissol::flops::godunovStateMinus_nonZero;
  o_hardwareFlops += seissol::flops::godunovStateMinus_hardware;
  
  o_nonZeroFlops *= CONVERGENCE_ORDER;
  o_hardwareFlops *= CONVERGENCE_ORDER;
}
