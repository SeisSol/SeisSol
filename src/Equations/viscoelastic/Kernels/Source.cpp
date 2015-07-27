/**
 * @file
 * This file is part of SeisSol.
 *
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
 * Source kernel of SeisSol.
 **/

#include "Source.h"

#include <matrix_kernels/dense.h>

#ifndef NDEBUG
#pragma message "compiling source kernel with assertions"
#endif

#include <cassert>
#include <stdint.h>

seissol::kernels::Source::Source() {
#define VOLUME_KERNEL
#include <initialization/bind.h>
#undef VOLUME_KERNEL
}
#include <iostream>
void seissol::kernels::Source::computeIntegral( real*  timeIntegratedDegreesOfFreedom,
                                                real   sourceMatrix[NUMBER_OF_QUANTITIES * NUMBER_OF_QUANTITIES],
                                                real*  degreesOfFreedom ) {
  // assert alignments
  assert( ((uintptr_t)timeIntegratedDegreesOfFreedom) % ALIGNMENT == 0 );
  assert( ((uintptr_t)degreesOfFreedom)               % ALIGNMENT == 0 );
  
      // DEBUG
    /*for (unsigned i = 0; i < NUMBER_OF_QUANTITIES * NUMBER_OF_QUANTITIES; ++i) {
      
        std::cout << sourceMatrix[i] << " ";
    }
    std::cout << std::endl << std::endl;*/
    // END DEBUG
  
  m_matrixKernels[3]   ( timeIntegratedDegreesOfFreedom, sourceMatrix, degreesOfFreedom, NULL, NULL, NULL );
}
