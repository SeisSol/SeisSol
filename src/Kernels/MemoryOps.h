// SPDX-FileCopyrightText: 2016 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_KERNELS_MEMORYOPS_H_
#define SEISSOL_SRC_KERNELS_MEMORYOPS_H_

#include "Kernels/Precision.h"

#ifdef __AVX512F__

#include <immintrin.h>

#define DMO_INCREMENT64 8
#define DMO_STREAM64(IN, OUT) _mm512_stream_pd(OUT, _mm512_load_pd(IN));

#define DMO_INCREMENT32 16
#define DMO_STREAM32(IN, OUT) _mm512_stream_ps(OUT, _mm512_load_ps(IN));

#elif defined(__MIC__)

#include <immintrin.h>

#define DMO_INCREMENT64 8
#define DMO_STREAM64(IN, OUT) _mm512_storenrngo_pd(OUT, _mm512_load_pd(IN));

#define DMO_INCREMENT32 16
#define DMO_STREAM32(IN, OUT) _mm512_storenrngo_ps(OUT, _mm512_load_ps(IN));

#elif defined(__AVX__)

#include <immintrin.h>

#define DMO_INCREMENT64 4
#define DMO_STREAM64(IN, OUT) _mm256_stream_pd(OUT, _mm256_load_pd(IN));

#define DMO_INCREMENT32 8
#define DMO_STREAM32(IN, OUT) _mm256_stream_ps(OUT, _mm256_load_ps(IN));

#elif defined(__SSE3__)

#include <immintrin.h>

#define DMO_INCREMENT64 2
#define DMO_STREAM64(IN, OUT) _mm_stream_pd(OUT, _mm_load_pd(IN));

#define DMO_INCREMENT32 4
#define DMO_STREAM32(IN, OUT) _mm_stream_ps(OUT, _mm_load_ps(IN));

#elif defined(__ARM_FEATURE_SVE)

#include <arm_sve.h>

#define DMO_INCREMENT64 svcntd()
#define DMO_STREAM64(IN, OUT) svstnt1_f64(svptrue_b64(), OUT, svld1_f64(svptrue_b64(), IN));

#define DMO_INCREMENT32 svcntw()
#define DMO_STREAM32(IN, OUT) svstnt1_f32(svptrue_b32(), OUT, svld1_f32(svptrue_b32(), IN));

#elif defined(__aarch64__)
// cf. https://stackoverflow.com/a/61248308

#define DMO_INCREMENT64 2
#define DMO_STREAM64(IN, OUT)                                                                      \
  uint64_t v1 = *(reinterpret_cast<const uint64_t*>(IN));                                          \
  uint64_t v2 = *(reinterpret_cast<const uint64_t*>(IN) + 1);                                      \
  asm volatile(                                                                                    \
      "stnp %[IN1],%[IN2],[%[OUTaddr]]" ::[IN1] "r"(v1), [IN2] "r"(v2), [OUTaddr] "r"(OUT)         \
      :);

#define DMO_INCREMENT32 4
#define DMO_STREAM32(IN, OUT)                                                                      \
  uint64_t v1 = *(reinterpret_cast<const uint64_t*>(IN));                                          \
  uint64_t v2 = *(reinterpret_cast<const uint64_t*>(IN) + 1);                                      \
  asm volatile(                                                                                    \
      "stnp %[IN1],%[IN2],[%[OUTaddr]]" ::[IN1] "r"(v1), [IN2] "r"(v2), [OUTaddr] "r"(OUT)         \
      :);

#else

#define DMO_INCREMENT64 1
#define DMO_STREAM64(IN, OUT) *(OUT) = *(IN);
#define DMO_INCREMENT32 1
#define DMO_STREAM32(IN, OUT) *(OUT) = *(IN);

#endif

#include <cassert>

namespace seissol::kernels {
/** Stores X in Y with non-temporal hint.
 *
 * @param numberOfReals The size of X and Y.
 * @param X
 * @param Y
 */
template <typename T>
inline void streamstore(std::size_t numberOfReals, const T* x, T* y);

template <>
inline void streamstore<float>(std::size_t numberOfReals, const float* x, float* y) {
  assert(numberOfReals % DMO_INCREMENT32 == 0);

  for (std::size_t i = 0; i < numberOfReals; i += DMO_INCREMENT32) {
    DMO_STREAM32(&x[i], &y[i])
  }
}

template <>
inline void streamstore<double>(std::size_t numberOfReals, const double* x, double* y) {
  assert(numberOfReals % DMO_INCREMENT64 == 0);

  for (std::size_t i = 0; i < numberOfReals; i += DMO_INCREMENT64) {
    DMO_STREAM64(&x[i], &y[i])
  }
}
} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_MEMORYOPS_H_
