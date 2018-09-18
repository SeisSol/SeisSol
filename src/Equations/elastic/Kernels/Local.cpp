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

#include <generated_code/kernel.h>
#include <yateto.h>


#include <cassert>
#include <stdint.h>
#include <cstring>

void seissol::kernels::Local::computeIntegral(  enum faceType const         i_faceTypes[4],
                                                GlobalData const*           global,
                                                LocalIntegrationData const* local,
                                                real*                       i_timeIntegratedDegreesOfFreedom,
                                                real*                       io_degreesOfFreedom )
{
  // assert alignments
#ifndef NDEBUG
  //~ for (unsigned stiffness = 0; stiffness < 3; ++stiffness) {
    //~ assert( ((uintptr_t)global->stiffnessMatrices[stiffness]) % ALIGNMENT == 0 );
  //~ }
  //~ for (unsigned flux = 0; flux < 4; ++flux) {
    //~ assert( ((uintptr_t)global->localChangeOfBasisMatricesTransposed[flux]) % ALIGNMENT == 0 );
    //~ assert( ((uintptr_t)global->changeOfBasisMatrices[flux]) % ALIGNMENT == 0 );
  //~ }
  assert( ((uintptr_t)i_timeIntegratedDegreesOfFreedom) % ALIGNMENT == 0 );
  assert( ((uintptr_t)io_degreesOfFreedom)              % ALIGNMENT == 0 );
#endif

  kernel::volume volKrnl;
  volKrnl.Q = io_degreesOfFreedom;
  volKrnl.I = i_timeIntegratedDegreesOfFreedom;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::kDivM>(); ++i) {
    volKrnl.kDivM[i] = global->stiffnessMatrices[i];
  }
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
    volKrnl.star[i] = local->starMatrices[i];
  }
  
  kernel::localFlux lfKrnl;
  lfKrnl.Q = io_degreesOfFreedom;
  lfKrnl.I = i_timeIntegratedDegreesOfFreedom;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::AplusT>(); ++i) {
    lfKrnl.AplusT[i] = local->nApNm1[i];
  }
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::rDivM>(); ++i) {
    lfKrnl.rDivM[i] = global->changeOfBasisMatrices[i];
  }
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::rDivM>(); ++i) {
    lfKrnl.fMrT[i] = global->localChangeOfBasisMatricesTransposed[i];
  }

  
  volKrnl.execute();
  
  for( unsigned int face = 0; face < 4; face++ ) {
    real const* prefetch = NULL;
    if (face == 0) {
      prefetch = i_timeIntegratedDegreesOfFreedom + NUMBER_OF_ALIGNED_DOFS;
    } else if (face == 1) {
      prefetch = io_degreesOfFreedom + NUMBER_OF_ALIGNED_DOFS;
    }
    // no element local contribution in the case of dynamic rupture boundary conditions
    if( i_faceTypes[face] != dynamicRupture ) {
      (lfKrnl.*lfKrnl.findExecute(face))();
    }
  }
}

void seissol::kernels::Local::flopsIntegral(  enum faceType const i_faceTypes[4],
                                              unsigned int        &o_nonZeroFlops,
                                              unsigned int        &o_hardwareFlops )
{
  o_nonZeroFlops = seissol::kernel::volume::NonZeroFlops;
  o_hardwareFlops = seissol::kernel::volume::HardwareFlops;

  for( unsigned int face = 0; face < 4; ++face ) {
    if( i_faceTypes[face] != dynamicRupture ) {
      o_nonZeroFlops  += seissol::kernel::localFlux::nonZeroFlops(face);
      o_hardwareFlops += seissol::kernel::localFlux::hardwareFlops(face);
    }
  }
}

unsigned seissol::kernels::Local::bytesIntegral()
{
  unsigned reals = 0;

  // star matrices load, source matrix load
  reals += yateto::computeFamilySize<tensor::star>();
  // flux solvers
  reals += yateto::computeFamilySize<tensor::AplusT>();

  // DOFs write
  reals += NUMBER_OF_ALIGNED_DOFS;
  
  return reals * sizeof(real);
}
