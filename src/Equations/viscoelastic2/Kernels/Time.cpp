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

  // initialize derivatives
  for( unsigned int l_dof = 0; l_dof < NUMBER_OF_ALIGNED_DOFS; ++l_dof ) {
    l_derivativesBuffer[l_dof] = i_degreesOfFreedom[l_dof];
    o_timeIntegrated[l_dof]  = i_degreesOfFreedom[l_dof] * l_scalar;
  }
  for( unsigned int l_dof = NUMBER_OF_ALIGNED_DOFS; l_dof < NUMBER_OF_ALIGNED_DERS; l_dof++ ) {
    l_derivativesBuffer[l_dof] = 0.0;
  }
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
  libxsmm_num_total_flops += NUMBER_OF_ALIGNED_DOFS;
#endif

  // stream out frist derivative (order 0)
  if ( o_timeDerivatives != NULL ) {
    streamstore(NUMBER_OF_ALIGNED_DOFS, i_degreesOfFreedom, o_timeDerivatives);
  }

  // compute all derivatives and contributions to the time integrated DOFs
  for( unsigned l_derivative = 1; l_derivative < CONVERGENCE_ORDER; l_derivative++ ) {
    real const* lastDerivative = l_derivativesBuffer + (l_derivative-1) * NUMBER_OF_ALIGNED_DOFS;
    real* currentDerivative = l_derivativesBuffer + l_derivative * NUMBER_OF_ALIGNED_DOFS;
    seissol::generatedKernels::derivative(
      local->starMatrices[0],
      local->starMatrices[1],
      local->starMatrices[2],
      global->stiffnessMatricesTransposed[1],
      global->stiffnessMatricesTransposed[0],
      global->stiffnessMatricesTransposed[2],
      lastDerivative,
      currentDerivative
    );

    for (int mech = NUMBER_OF_RELAXATION_MECHANISMS-1; mech >= 0; --mech) {
      unsigned mechOffset = NUMBER_OF_ALIGNED_ELASTIC_DOFS + mech * NUMBER_OF_ALIGNED_MECHANISM_DOFS;
      
      seissol::generatedKernels::source(  local->specific.ET + mech * seissol::model::ET::reals,
                                          lastDerivative + mechOffset,
                                          currentDerivative );

      XYmSt(  local->specific.omega[mech],
              NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
              NUMBER_OF_MECHANISM_QUANTITIES,
              currentDerivative + NUMBER_OF_ALIGNED_ELASTIC_DOFS,
              NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
              lastDerivative + mechOffset,
              NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
              currentDerivative + mechOffset,
              NUMBER_OF_ALIGNED_BASIS_FUNCTIONS );
    }

    l_scalar *= i_timeStepWidth / real(l_derivative+1);
    
    real* derivativesStore = NULL;
    if (o_timeDerivatives != NULL) {
      derivativesStore = o_timeDerivatives + l_derivative * NUMBER_OF_ALIGNED_DOFS;
    }

    SXtYp(  l_scalar,
            NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
            NUMBER_OF_QUANTITIES,
            currentDerivative,
            NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
            o_timeIntegrated,
            NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
            derivativesStore );
  }
}

void seissol::kernels::Time::flopsAder( unsigned int        &o_nonZeroFlops,
                                        unsigned int        &o_hardwareFlops ) {
  // reset flops
  o_nonZeroFlops = 0; o_hardwareFlops = 0;
  
  // initialization
  o_nonZeroFlops += NUMBER_OF_DOFS;
  o_hardwareFlops += NUMBER_OF_ALIGNED_DOFS;
  
  for( unsigned l_derivative = 1; l_derivative < CONVERGENCE_ORDER; ++l_derivative ) {
    o_nonZeroFlops  += seissol::flops::derivative_nonZero;
    o_hardwareFlops += seissol::flops::derivative_hardware;
    
    o_nonZeroFlops  += NUMBER_OF_RELAXATION_MECHANISMS * seissol::flops::source_nonZero;
    o_hardwareFlops += NUMBER_OF_RELAXATION_MECHANISMS * seissol::flops::source_hardware;
    
    // XYmSt
    o_nonZeroFlops  += NUMBER_OF_RELAXATION_MECHANISMS * NUMBER_OF_BASIS_FUNCTIONS * NUMBER_OF_MECHANISM_QUANTITIES * 2;
    o_hardwareFlops += NUMBER_OF_RELAXATION_MECHANISMS * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS * NUMBER_OF_MECHANISM_QUANTITIES * 2;

    // update of time integrated DOFs
    o_nonZeroFlops  += NUMBER_OF_DOFS * 2;
    o_hardwareFlops += NUMBER_OF_ALIGNED_DOFS * 2;
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
           + seissol::model::CstarT::reals
           + NUMBER_OF_RELAXATION_MECHANISMS * seissol::model::ET::reals;
           
  /// \todo incorporate derivatives

  return reals * sizeof(real);
}

void seissol::kernels::Time::computeIntegral( double                                      i_expansionPoint,
                                              double                                      i_integrationStart,
                                              double                                      i_integrationEnd,
                                              real const*                                 i_timeDerivatives,
                                              real                                        o_timeIntegrated[NUMBER_OF_ALIGNED_DOFS] )
{

  /*
   * assert alignments.
   */
  assert( ((uintptr_t)i_timeDerivatives)  % ALIGNMENT == 0 );
  assert( ((uintptr_t)o_timeIntegrated)   % ALIGNMENT == 0 );

  // compute lengths of integration intervals
  real l_deltaTLower = i_integrationStart - i_expansionPoint;
  real l_deltaTUpper = i_integrationEnd   - i_expansionPoint;
  
  // initialization of scalars in the taylor series expansion (0th term)
  real l_firstTerm  = l_deltaTUpper;
  real l_secondTerm = l_deltaTLower;
  real l_factorial  = 1.0;
  real l_scalar = l_deltaTUpper - l_deltaTLower;

  for( unsigned int l_dof = 0; l_dof < NUMBER_OF_ALIGNED_DOFS; ++l_dof ) {
    o_timeIntegrated[l_dof]  = i_timeDerivatives[l_dof] * l_scalar;
  }
 
  // iterate over time derivatives
  for(int l_derivative = 1; l_derivative < CONVERGENCE_ORDER; l_derivative++ ) {
    l_firstTerm  *= l_deltaTUpper;
    l_secondTerm *= l_deltaTLower;
    l_factorial  *= (real)(l_derivative+1);
    l_scalar  = (l_firstTerm - l_secondTerm) / l_factorial;

    SXtYp(  l_scalar,
            NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
            NUMBER_OF_QUANTITIES,
            i_timeDerivatives + l_derivative * NUMBER_OF_ALIGNED_DOFS,
            NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
            o_timeIntegrated,
            NUMBER_OF_ALIGNED_BASIS_FUNCTIONS );
  }
}
