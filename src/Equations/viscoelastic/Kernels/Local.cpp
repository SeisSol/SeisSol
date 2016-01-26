/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2013-2014, SeisSol Group
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
 * Local kernel of SeisSol.
 **/

#include "Local.h"

#ifndef NDEBUG
#pragma message "compiling local kernel with assertions"
#endif

#include <cassert>
#include <stdint.h>

#include <generated_code/kernels.h>
#include <generated_code/flops.h>

void seissol::kernels::Local::computeIntegral(  enum faceType const         i_faceTypes[4],
                                                GlobalData const*           global,
                                                LocalIntegrationData const* local,
                                                real*                       i_timeIntegratedDegreesOfFreedom,
                                                real*                       io_degreesOfFreedom ) {
  // assert alignments
#ifndef NDEBUG
  for (unsigned stiffness = 0; stiffness < 3; ++stiffness) {
    assert( ((uintptr_t)global->stiffnessMatrices[stiffness]) % ALIGNMENT == 0 );
  }
  for (unsigned flux = 0; flux < 4; ++flux) {
    assert( ((uintptr_t)global->fluxMatrices[flux]) % ALIGNMENT == 0 );
  }
  assert( ((uintptr_t)i_timeIntegratedDegreesOfFreedom) % ALIGNMENT == 0 );
  assert( ((uintptr_t)io_degreesOfFreedom)              % ALIGNMENT == 0 );
#endif
  
  seissol::generatedKernels::volume(
    local->starMatrices[0],
    local->starMatrices[1],
    local->starMatrices[2],
    global->stiffnessMatrices[1],
    global->stiffnessMatrices[0],
    global->stiffnessMatrices[2],
    local->specific.sourceMatrix,
    i_timeIntegratedDegreesOfFreedom,
    io_degreesOfFreedom
  );

  for( unsigned int face = 0; face < 4; ++face ) {
    // no element local contribution in the case of dynamic rupture boundary conditions
    if( i_faceTypes[face] != dynamicRupture ) {
      seissol::generatedKernels::localFlux[face](
        local->nApNm1[face],
        global->fluxMatrices[face],
        i_timeIntegratedDegreesOfFreedom,
        io_degreesOfFreedom
      );
    }
  }
}

void seissol::kernels::Local::flopsIntegral(  enum faceType const i_faceTypes[4],
                                              unsigned int        &o_nonZeroFlops,
                                              unsigned int        &o_hardwareFlops )
{
  o_nonZeroFlops = seissol::flops::volume_nonZero;
  o_hardwareFlops = seissol::flops::volume_hardware;

  for( unsigned int face = 0; face < 4; ++face ) {
    if( i_faceTypes[face] != dynamicRupture ) {
      o_nonZeroFlops  += seissol::flops::localFlux_nonZero[face];
      o_hardwareFlops += seissol::flops::localFlux_hardware[face];
    }
  }
}

unsigned seissol::kernels::Local::bytesIntegral()
{
  unsigned reals = 0;

  // star matrices load, source matrix load
  reals += seissol::model::AstarT::reals
           + seissol::model::BstarT::reals
           + seissol::model::CstarT::reals
           + seissol::model::source::reals;
  // flux solvers
  reals += 4 * seissol::model::AplusT::reals;

  // DOFs write
  reals += NUMBER_OF_ALIGNED_DOFS;
  
  return reals * sizeof(real);
}
