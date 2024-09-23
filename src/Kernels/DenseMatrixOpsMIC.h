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

#ifndef SEISSOL_SRC_KERNELS_DENSEMATRIXOPSMIC_H_
#define SEISSOL_SRC_KERNELS_DENSEMATRIXOPSMIC_H_

#if defined(DOUBLE_PRECISION)

#define DMO_INCREMENT 8
#define DMO_BROADCAST(IN, OUT) __m512d OUT = _mm512_extload_pd(IN, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
#define DMO_STREAM(IN, OUT) _mm512_storenrngo_pd(OUT, _mm512_load_pd(IN));
#define DMO_SXT(S, X, Y) __m512d x = _mm512_load_pd(X); _mm512_store_pd(Y, _mm512_mul_pd(S, x));
#define DMO_SXTYP(S, X, Y) __m512d x = _mm512_load_pd(X); _mm512_store_pd(Y, _mm512_fmadd_pd(S, x, _mm512_load_pd(Y)));
#define DMO_XYMST(S, X, Y, Z) __m512d x = _mm512_load_pd(X); __m512d y = _mm512_load_pd(Y); _mm512_storenrngo_pd(Z, _mm512_mul_pd(S, _mm512_sub_pd(x, y)));
#define DMO_XYMSTZP(S, X, Y, Z) __m512d x = _mm512_load_pd(X); __m512d y = _mm512_load_pd(Y); _mm512_store_pd(Z, _mm512_fmadd_pd(S, _mm512_sub_pd(x, y), _mm512_load_pd(Z)));

#elif defined(SINGLE_PRECISION)

#define DMO_INCREMENT 16
#define DMO_BROADCAST(IN, OUT) __m512 OUT = _mm512_extload_ps(IN, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);
#define DMO_STREAM(IN, OUT) _mm512_storenrngo_ps(OUT, _mm512_load_ps(IN));
#define DMO_SXT(S, X, Y) __m512 x = _mm512_load_ps(X); _mm512_store_ps(Y, _mm512_mul_ps(S, x));
#define DMO_SXTYP(S, X, Y) __m512 x = _mm512_load_ps(X); _mm512_store_ps(Y, _mm512_fmadd_ps(S, x, _mm512_load_ps(Y)));
#define DMO_XYMST(S, X, Y, Z) __m512 x = _mm512_load_ps(X); __m512 y = _mm512_load_ps(Y); _mm512_storenrngo_ps(Z, _mm512_mul_ps(S, _mm512_sub_ps(x, y)));
#define DMO_XYMSTZP(S, X, Y, Z) __m512 x = _mm512_load_ps(X); __m512 y = _mm512_load_ps(Y); _mm512_store_ps(Z, _mm512_fmadd_ps(S, _mm512_sub_ps(x, y), _mm512_load_ps(Z)));

#else
#error no precision was defined
#endif

#endif // SEISSOL_SRC_KERNELS_DENSEMATRIXOPSMIC_H_

