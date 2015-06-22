/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
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
 * Volume kernel of SeisSol.
 **/

#include "Volume.h"

#include <matrix_kernels/sparse.h>
#include <matrix_kernels/dense.h>

#ifndef NDEBUG
#pragma message "compiling volume kernel with assertions"
#endif

#include <cassert>
#include <stdint.h>

seissol::kernels::Volume::Volume() {
  // intialize the function pointers to the matrix kernels
#define VOLUME_KERNEL
#include <initialization/bind.h>
#undef VOLUME_KERNEL
}

void seissol::kernels::Volume::computeIntegral( real** i_stiffnessMatrices,
                                                real*  i_timeIntegratedDegreesOfFreedom,
                                                real   i_starMatrices[3][STAR_NNZ],
                                                real*  io_degreesOfFreedom ) {
  // assert alignments
  assert( ((uintptr_t)i_stiffnessMatrices[0])           % ALIGNMENT == 0 );
  assert( ((uintptr_t)i_stiffnessMatrices[1])           % ALIGNMENT == 0 );
  assert( ((uintptr_t)i_stiffnessMatrices[2])           % ALIGNMENT == 0 );
  assert( ((uintptr_t)i_timeIntegratedDegreesOfFreedom) % ALIGNMENT == 0 );
  assert( ((uintptr_t)io_degreesOfFreedom)              % ALIGNMENT == 0 );

  // temporary result
  real l_temporaryResult[NUMBER_OF_ALIGNED_DOFS] __attribute__((aligned(4096)));

  // iterate over dimensions 
  for( unsigned int l_c = 0; l_c < 3; l_c++ ) {
    m_matrixKernels[l_c] ( i_stiffnessMatrices[l_c], i_timeIntegratedDegreesOfFreedom, l_temporaryResult,
                           NULL,                     NULL,                             NULL                 ); // These will be be ignored
    m_matrixKernels[3]   ( l_temporaryResult,        i_starMatrices[l_c],              io_degreesOfFreedom,
                           NULL,                     NULL,                             NULL                 ); // These will be be ignored
  }
}

void seissol::kernels::Volume::flopsIntegral( unsigned int &o_nonZeroFlops,
                                              unsigned int &o_hardwareFlops ) {
  // reset flops
  o_nonZeroFlops = 0; o_hardwareFlops = 0;

  // iterate over dimensions
  for( unsigned int l_c = 0; l_c < 3; l_c++ ) {
    o_nonZeroFlops  += m_nonZeroFlops[ l_c];
    o_hardwareFlops += m_hardwareFlops[l_c];

    o_nonZeroFlops  += m_nonZeroFlops[3];
    o_hardwareFlops += m_hardwareFlops[3];
  }
}
