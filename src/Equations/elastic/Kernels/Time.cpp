/******************************************************************************
** Copyright (c) 2014-2015, Intel Corporation                                **
** All rights reserved.                                                      **
**                                                                           **
** Redistribution and use in source and binary forms, with or without        **
** modification, are permitted provided that the following conditions        **
** are met:                                                                  **
** 1. Redistributions of source code must retain the above copyright         **
**    notice, this list of conditions and the following disclaimer.          **
** 2. Redistributions in binary form must reproduce the above copyright      **
**    notice, this list of conditions and the following disclaimer in the    **
**    documentation and/or other materials provided with the distribution.   **
** 3. Neither the name of the copyright holder nor the names of its          **
**    contributors may be used to endorse or promote products derived        **
**    from this software without specific prior written permission.          **
**                                                                           **
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       **
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT         **
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR     **
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT      **
** HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,    **
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED  **
** TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR    **
** PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    **
** LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      **
** NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        **
** SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              **
******************************************************************************/
/* Alexander Heinecke (Intel Corp.)
******************************************************************************/
/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2013-2015, SeisSol Group
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
 * Time kernel of SeisSol.
 **/

#include "Time.h"

#ifndef NDEBUG
#pragma message "compiling time kernel with assertions"
extern long long libxsmm_num_total_flops;
#endif

#include <Kernels/denseMatrixOps.hpp>
#include <generated_code/kernels.h>
#include <generated_code/flops.h>

#include <cstring>
#include <cassert>
#include <stdint.h>

seissol::kernels::Time::Time() {
  // compute the aligned number of basis functions and offsets of the derivatives
  m_derivativesOffsets[0] = 0;
  for( int l_order = 0; l_order < CONVERGENCE_ORDER; l_order++ ) {
    m_numberOfAlignedBasisFunctions[l_order] = getNumberOfAlignedBasisFunctions( CONVERGENCE_ORDER-l_order, ALIGNMENT );

    if( l_order > 0 ) {
      m_derivativesOffsets[l_order]  =  m_numberOfAlignedBasisFunctions[l_order-1] * NUMBER_OF_QUANTITIES;
      m_derivativesOffsets[l_order] +=  m_derivativesOffsets[l_order-1];
    }
  }
}

void seissol::kernels::Time::computeAder( double                      i_timeStepWidth,
                                          GlobalData const*           global,
                                          LocalIntegrationData const* local,
                                          real const*                 i_degreesOfFreedom,
                                          real*                       o_timeIntegrated,
                                          real*                       o_timeDerivatives )
{
  /*
   * assert alignments.
   */
  assert( ((uintptr_t)i_degreesOfFreedom)     % ALIGNMENT == 0 );
  assert( ((uintptr_t)global->stiffnessMatricesTransposed[0]) % ALIGNMENT == 0 );
  assert( ((uintptr_t)global->stiffnessMatricesTransposed[1]) % ALIGNMENT == 0 );
  assert( ((uintptr_t)global->stiffnessMatricesTransposed[2]) % ALIGNMENT == 0 );
  assert( ((uintptr_t)o_timeIntegrated )      % ALIGNMENT == 0 );
  assert( ((uintptr_t)o_timeDerivatives)      % ALIGNMENT == 0 || o_timeDerivatives == NULL );

  /*
   * compute ADER scheme.
   */
  // scalars in the taylor-series expansion
  real l_scalar = i_timeStepWidth;

  // temporary result
  real l_derivativesBuffer[NUMBER_OF_ALIGNED_DERS] __attribute__((aligned(PAGESIZE_STACK)));

  // initialize time integrated DOFs
  SXt(  l_scalar,
        NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
        NUMBER_OF_QUANTITIES,
        i_degreesOfFreedom,
        NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
        o_timeIntegrated,
        NUMBER_OF_ALIGNED_BASIS_FUNCTIONS  );

  // stream out frist derivative (order 0)
  if ( o_timeDerivatives != NULL ) {
    streamstore(NUMBER_OF_ALIGNED_DOFS, i_degreesOfFreedom, o_timeDerivatives);
  }
  
  real const* lastDerivative = i_degreesOfFreedom;
  for( unsigned l_derivative = 1; l_derivative < CONVERGENCE_ORDER; l_derivative++ ) {
    real* currentDerivative = l_derivativesBuffer + m_derivativesOffsets[l_derivative];
    seissol::generatedKernels::derivative[l_derivative](
      local->starMatrices[0],
      local->starMatrices[1],
      local->starMatrices[2],
      global->stiffnessMatricesTransposed[1],
      global->stiffnessMatricesTransposed[0],
      global->stiffnessMatricesTransposed[2],
      lastDerivative,
      currentDerivative
    );
    lastDerivative = currentDerivative;

    // update scalar for this derivative
    l_scalar *= i_timeStepWidth / real(l_derivative+1);

    // stream out derivatives if not NULL
    real* derivativesStore = NULL;
    if (o_timeDerivatives != NULL) {
      derivativesStore = o_timeDerivatives + m_derivativesOffsets[l_derivative];
    }

    // update time integrated DOFs
    SXtYp(  l_scalar,
            m_numberOfAlignedBasisFunctions[l_derivative],
            NUMBER_OF_QUANTITIES,
            currentDerivative,
            m_numberOfAlignedBasisFunctions[l_derivative],
            o_timeIntegrated,
            NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
            derivativesStore );
  }
}

void seissol::kernels::Time::flopsAder( unsigned int        &o_nonZeroFlops,
                                        unsigned int        &o_hardwareFlops ) {
  // reset flops
  o_nonZeroFlops = 0; o_hardwareFlops =0;

  // initialization
  o_nonZeroFlops  += NUMBER_OF_DOFS;
  o_hardwareFlops += NUMBER_OF_ALIGNED_DOFS;

  // interate over derivatives
  for( unsigned l_derivative = 1; l_derivative < CONVERGENCE_ORDER; l_derivative++ ) {
    o_nonZeroFlops  += seissol::flops::derivative_nonZero[l_derivative];
    o_hardwareFlops += seissol::flops::derivative_hardware[l_derivative];

    // update of time integrated DOFs
    o_nonZeroFlops  += seissol::kernels::getNumberOfBasisFunctions(        CONVERGENCE_ORDER - l_derivative ) * NUMBER_OF_QUANTITIES * 2;
    o_hardwareFlops += seissol::kernels::getNumberOfAlignedBasisFunctions( CONVERGENCE_ORDER - l_derivative ) * NUMBER_OF_QUANTITIES * 2;
  }

}

unsigned seissol::kernels::Time::bytesAder()
{
  unsigned reals = 0;
  
  // DOFs load, tDOFs load, tDOFs write
  reals += 3 * NUMBER_OF_ALIGNED_DOFS;
  // star matrices, source matrix
  reals += seissol::model::AstarT::reals
           + seissol::model::BstarT::reals
           + seissol::model::CstarT::reals;
           
  /// \todo incorporate derivatives

  return reals * sizeof(real);
}

void seissol::kernels::Time::computeIntegral( double                            i_expansionPoint,
                                              double                            i_integrationStart,
                                              double                            i_integrationEnd,
                                              const real*                       i_timeDerivatives,
                                              real                              o_timeIntegrated[NUMBER_OF_ALIGNED_DOFS] )
{
  /*
   * assert alignments.
   */
  assert( ((uintptr_t)i_timeDerivatives)  % ALIGNMENT == 0 );
  assert( ((uintptr_t)o_timeIntegrated)   % ALIGNMENT == 0 );

  // assert that this is a forwared integration in time
  assert( i_integrationStart + (real) 1.E-10 > i_expansionPoint   );
  assert( i_integrationEnd                   > i_integrationStart );

  /*
   * compute time integral.
   */
  // reset time integrated degrees of freedom
  memset( o_timeIntegrated, 0, NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES*sizeof(real) );

  // compute lengths of integration intervals
  real l_deltaTLower = i_integrationStart - i_expansionPoint;
  real l_deltaTUpper = i_integrationEnd   - i_expansionPoint;

  // initialization of scalars in the taylor series expansion (0th term)
  real l_firstTerm  = (real) 1;
  real l_secondTerm = (real) 1;
  real l_factorial  = (real) 1;
  real l_scalar;
 
  // iterate over time derivatives
  for(int l_derivative = 0; l_derivative < CONVERGENCE_ORDER; l_derivative++ ) {
    l_firstTerm  *= l_deltaTUpper;
    l_secondTerm *= l_deltaTLower;
    l_factorial  *= (real)(l_derivative+1);

    l_scalar  = l_firstTerm - l_secondTerm;
    l_scalar /= l_factorial;

    SXtYp(  l_scalar,
            m_numberOfAlignedBasisFunctions[l_derivative],
            NUMBER_OF_QUANTITIES,
            i_timeDerivatives + m_derivativesOffsets[l_derivative],
            m_numberOfAlignedBasisFunctions[l_derivative],
            o_timeIntegrated,
            NUMBER_OF_ALIGNED_BASIS_FUNCTIONS );
  }
}
