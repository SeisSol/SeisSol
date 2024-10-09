/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
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
 * Common kernel-level functions
 **/

#ifndef KERNElS_DENSEMATRIXOPS_HPP_
#define KERNElS_DENSEMATRIXOPS_HPP_

#if defined(__SSE3__) || defined(__MIC__)
#include <immintrin.h>
#endif

#if defined(__AVX512F__)
#include "DenseMatrixOpsAVX512.h"
#elif defined(__MIC__)
#include "DenseMatrixOpsMIC.h"
#elif defined(__AVX2__)
#include "DenseMatrixOpsAVX2.h"
#elif defined(__AVX__)
#include "DenseMatrixOpsAVX.h"
#elif defined(__SSE3__)
#include "DenseMatrixOpsSSE3.h"
#elif defined(__ARM_FEATURE_SVE)
#include "DenseMatrixOpsSVE.h"
#elif defined(__aarch64__)
#include "DenseMatrixOpsAARCH64.h"
#else
#include "DenseMatrixOpsNoarch.h"
#endif

#include <cassert>

namespace seissol::kernels {
/** Stores X in Y with non-temporal hint.
 *
 * @param numberOfReals The size of X and Y.
 * @param X
 * @param Y
 */
inline void streamstore(std::size_t numberOfReals, const real* x, real* y) {
  assert(numberOfReals % DMO_INCREMENT == 0);

  for (std::size_t i = 0; i < numberOfReals; i += DMO_INCREMENT) {
    DMO_STREAM(&x[i], &y[i])
  }
}
} // namespace seissol::kernels

#endif
