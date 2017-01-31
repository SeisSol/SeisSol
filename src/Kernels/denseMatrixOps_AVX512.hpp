/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
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

#include <Kernels/precision.hpp>

#if defined(DOUBLE_PRECISION)

#define DMO_INCREMENT 8
#define DMO_BROADCAST(IN, OUT) __m512d OUT = _mm512_broadcastsd_pd(_mm_load_sd(IN));
#define DMO_STREAM(IN, OUT) _mm512_stream_pd(OUT, _mm512_load_pd(IN));
#define DMO_SXT(S, X, Y) __m512d x = _mm512_load_pd(X); _mm512_store_pd(Y, _mm512_mul_pd(S, x));
#define DMO_SXTYP(S, X, Y) __m512d x = _mm512_load_pd(X); _mm512_store_pd(Y, _mm512_fmadd_pd(S, x, _mm512_load_pd(Y)));
#define DMO_XYMST(S, X, Y, Z) __m512d x = _mm512_load_pd(X); __m512d y = _mm512_load_pd(Y); _mm512_store_pd(Z, _mm512_mul_pd(S, _mm512_sub_pd(x, y)));
#define DMO_XYMSTZP(S, X, Y, Z) __m512d x = _mm512_load_pd(X); __m512d y = _mm512_load_pd(Y); _mm512_store_pd(Z, _mm512_fmadd_pd(S, _mm512_sub_pd(x, y), _mm512_load_pd(Z)));

#elif defined(SINGLE_PRECISION)

#define DMO_INCREMENT 16
#define DMO_BROADCAST(IN, OUT) __m512 OUT = _mm512_broadcastss_ps(_mm_load_ss(IN));
#define DMO_STREAM(IN, OUT) _mm512_stream_ps(OUT, _mm512_load_ps(IN));
#define DMO_SXT(S, X, Y) __m512 x = _mm512_load_ps(X); _mm512_store_ps(Y, _mm512_mul_ps(S, x));
#define DMO_SXTYP(S, X, Y) __m512 x = _mm512_load_ps(X); _mm512_store_ps(Y, _mm512_fmadd_ps(S, x, _mm512_load_ps(Y)));
#define DMO_XYMST(S, X, Y, Z) __m512 x = _mm512_load_ps(X); __m512 y = _mm512_load_ps(Y); _mm512_store_ps(Z, _mm512_mul_ps(S, _mm512_sub_ps(x, y)));
#define DMO_XYMSTZP(S, X, Y, Z) __m512 x = _mm512_load_ps(X); __m512 y = _mm512_load_ps(Y); _mm512_store_ps(Z, _mm512_fmadd_ps(S, _mm512_sub_ps(x, y), _mm512_load_ps(Z)));

#else
#error no precision was defined
#endif
