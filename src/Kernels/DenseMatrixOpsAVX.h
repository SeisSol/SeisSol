// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 */
 
#ifndef SEISSOL_SRC_KERNELS_DENSEMATRIXOPSAVX_H_
#define SEISSOL_SRC_KERNELS_DENSEMATRIXOPSAVX_H_

#if defined(DOUBLE_PRECISION)

#define DMO_INCREMENT 4
#define DMO_BROADCAST(IN, OUT) __m256d OUT = _mm256_broadcast_sd(IN);
#define DMO_STREAM(IN, OUT) _mm256_stream_pd(OUT, _mm256_load_pd(IN));
#define DMO_SXT(S, X, Y) __m256d x = _mm256_load_pd(X); _mm256_store_pd(Y, _mm256_mul_pd(S, x));
#define DMO_SXTYP(S, X, Y) __m256d x = _mm256_load_pd(X); _mm256_store_pd(Y, _mm256_add_pd(_mm256_mul_pd(S, x), _mm256_load_pd(Y)));
#define DMO_XYMST(S, X, Y, Z) __m256d x = _mm256_load_pd(X); __m256d y = _mm256_load_pd(Y); _mm256_store_pd(Z, _mm256_mul_pd(S, _mm256_sub_pd(x, y)));
#define DMO_XYMSTZP(S, X, Y, Z) __m256d x = _mm256_load_pd(X); __m256d y = _mm256_load_pd(Y); _mm256_store_pd(Z, _mm256_add_pd(_mm256_mul_pd(S, _mm256_sub_pd(x, y)), _mm256_load_pd(Z)));

#elif defined(SINGLE_PRECISION)

#define DMO_INCREMENT 8
#define DMO_BROADCAST(IN, OUT) __m256 OUT = _mm256_broadcast_ss(IN);
#define DMO_STREAM(IN, OUT) _mm256_stream_ps(OUT, _mm256_load_ps(IN));
#define DMO_SXT(S, X, Y) __m256 x = _mm256_load_ps(X); _mm256_store_ps(Y, _mm256_mul_ps(S, x));
#define DMO_SXTYP(S, X, Y) __m256 x = _mm256_load_ps(X); _mm256_store_ps(Y, _mm256_add_ps(_mm256_mul_ps(S, x), _mm256_load_ps(Y)));
#define DMO_XYMST(S, X, Y, Z) __m256 x = _mm256_load_ps(X); __m256 y = _mm256_load_ps(Y); _mm256_store_ps(Z, _mm256_mul_ps(S, _mm256_sub_ps(x, y)));
#define DMO_XYMSTZP(S, X, Y, Z) __m256 x = _mm256_load_ps(X); __m256 y = _mm256_load_ps(Y); _mm256_store_ps(Z, _mm256_add_ps(_mm256_mul_ps(S, _mm256_sub_ps(x, y)), _mm256_load_ps(Z)));

#else
#error no precision was defined
#endif

#endif // SEISSOL_SRC_KERNELS_DENSEMATRIXOPSAVX_H_

