// SPDX-FileCopyrightText: 2016-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_KERNELS_DENSEMATRIXOPS_H_
#define SEISSOL_SRC_KERNELS_DENSEMATRIXOPS_H_

#include <Kernels/Precision.h>

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

#endif // SEISSOL_SRC_KERNELS_DENSEMATRIXOPS_H_
