// Copyright (c) 2015, Intel Corporation
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 
//     * Redistributions of source code must retain the above copyright notice,
//       this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the copyright holder nor the names of its contributors
//       may be used to endorse or promote products derived from this software
//       without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// @file
// This file is part of SeisSol.
// 
// @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
// @author Alexander Heinecke (alexander.heinecke AT mytum.de, http://www5.in.tum.de/wiki/index.php/Alexander_Heinecke,_M.Sc.,_M.Sc._with_honors)
// 
// @date 2015-11-21 13:46:07.999565
// 
// @section LICENSE
// Copyright (c) 2012-2015, SeisSol Group
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
// 
// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from this
//    software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
// 
// @section DESCRIPTION
// Remark: This file was generated.
#ifndef SPARSESSNBCPP
#define SPARSESSNBCPP

#if defined( __SSE3__) || defined(__MIC__)
#include <immintrin.h>
#endif

#include <cstddef>
#ifndef NDEBUG
extern long long libxsmm_num_total_flops;
#endif

void ssparse_kXiDivMT_m1_n9_k4_ldAna2_ldB8_ldC8_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_n = 0;
  #pragma nounroll_and_jam
  for ( l_n = 0; l_n < 9; l_n++) {
    unsigned int l_m = 0;
    for ( l_m = 0; l_m < 1; l_m++) {
      C[(l_n*8)+l_m] = 0.0f;
    }
#if defined(__SSE3__) || defined(__AVX__)
#else
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b1 = _mm_broadcast_ss(&B[(l_n*8)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b1 = _mm_load_ss(&B[(l_n*8)+1]);    b1 = _mm_shuffle_ps(b1, b1, 0x00);
#endif
    __m128 c1_0 = _mm_load_ss(&C[(l_n*8)+0]);
    __m128 a1_0 = _mm_load_ss(&A[0]);
    c1_0 = _mm_add_ss(c1_0, _mm_mul_ss(a1_0, b1));
    _mm_store_ss(&C[(l_n*8)+0], c1_0);
#else
    C[(l_n*8)+0] += A[0] * B[(l_n*8)+1];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#else
#endif

#if defined(__SSE3__) || defined(__AVX__)
#else
#endif

  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 18;
#endif
}

void ssparse_starMatrix_m1_n9_k9_ldA8_ldBna2_ldC8_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  for ( l_m = 0; l_m < 1; l_m++) {
    C[0+l_m] += A[48+l_m] * B[0];
    C[0+l_m] += A[56+l_m] * B[1];
    C[0+l_m] += A[64+l_m] * B[2];
    C[8+l_m] += A[48+l_m] * B[3];
    C[8+l_m] += A[56+l_m] * B[4];
    C[8+l_m] += A[64+l_m] * B[5];
    C[16+l_m] += A[48+l_m] * B[6];
    C[16+l_m] += A[56+l_m] * B[7];
    C[16+l_m] += A[64+l_m] * B[8];
    C[24+l_m] += A[48+l_m] * B[9];
    C[24+l_m] += A[56+l_m] * B[10];
    C[32+l_m] += A[56+l_m] * B[11];
    C[32+l_m] += A[64+l_m] * B[12];
    C[40+l_m] += A[48+l_m] * B[13];
    C[40+l_m] += A[64+l_m] * B[14];
    C[48+l_m] += A[0+l_m] * B[15];
    C[48+l_m] += A[24+l_m] * B[16];
    C[48+l_m] += A[40+l_m] * B[17];
    C[56+l_m] += A[8+l_m] * B[18];
    C[56+l_m] += A[24+l_m] * B[19];
    C[56+l_m] += A[32+l_m] * B[20];
    C[64+l_m] += A[16+l_m] * B[21];
    C[64+l_m] += A[32+l_m] * B[22];
    C[64+l_m] += A[40+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 48;
#endif
}

void ssparse_kXiDivMT_m4_n9_k10_ldAna3_ldB16_ldC8_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_n = 0;
  #pragma nounroll_and_jam
  for ( l_n = 0; l_n < 9; l_n++) {
    unsigned int l_m = 0;
   #pragma simd
    for ( l_m = 0; l_m < 4; l_m++) {
      C[(l_n*8)+l_m] = 0.0f;
    }
#if defined(__SSE3__) || defined(__AVX__)
#else
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b1 = _mm_broadcast_ss(&B[(l_n*16)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b1 = _mm_load_ss(&B[(l_n*16)+1]);    b1 = _mm_shuffle_ps(b1, b1, 0x00);
#endif
    __m128 c1_0 = _mm_load_ss(&C[(l_n*8)+0]);
    __m128 a1_0 = _mm_load_ss(&A[0]);
    c1_0 = _mm_add_ss(c1_0, _mm_mul_ss(a1_0, b1));
    _mm_store_ss(&C[(l_n*8)+0], c1_0);
#else
    C[(l_n*8)+0] += A[0] * B[(l_n*16)+1];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#else
#endif

#if defined(__SSE3__) || defined(__AVX__)
#else
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b4 = _mm_broadcast_ss(&B[(l_n*16)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b4 = _mm_load_ss(&B[(l_n*16)+4]);    b4 = _mm_shuffle_ps(b4, b4, 0x00);
#endif
    __m128 c4_0 = _mm_load_ss(&C[(l_n*8)+1]);
    __m128 a4_0 = _mm_load_ss(&A[1]);
    c4_0 = _mm_add_ss(c4_0, _mm_mul_ss(a4_0, b4));
    _mm_store_ss(&C[(l_n*8)+1], c4_0);
#else
    C[(l_n*8)+1] += A[1] * B[(l_n*16)+4];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b5 = _mm_broadcast_ss(&B[(l_n*16)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b5 = _mm_load_ss(&B[(l_n*16)+5]);    b5 = _mm_shuffle_ps(b5, b5, 0x00);
#endif
    __m128 c5_0 = _mm_load_ss(&C[(l_n*8)+0]);
    __m128 a5_0 = _mm_load_ss(&A[2]);
    c5_0 = _mm_add_ss(c5_0, _mm_mul_ss(a5_0, b5));
    _mm_store_ss(&C[(l_n*8)+0], c5_0);
    __m128 c5_1 = _mm_castpd_ps(_mm_load_sd((const double*)&C[(l_n*8)+2]));
    __m128 a5_1 = _mm_castpd_ps(_mm_load_sd((const double*)&A[3]));
    c5_1 = _mm_add_ps(c5_1, _mm_mul_ps(a5_1, b5));
    _mm_store_sd((double*)&C[(l_n*8)+2], _mm_castps_pd(c5_1));
#else
    C[(l_n*8)+0] += A[2] * B[(l_n*16)+5];
    C[(l_n*8)+2] += A[3] * B[(l_n*16)+5];
    C[(l_n*8)+3] += A[4] * B[(l_n*16)+5];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#else
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b7 = _mm_broadcast_ss(&B[(l_n*16)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b7 = _mm_load_ss(&B[(l_n*16)+7]);    b7 = _mm_shuffle_ps(b7, b7, 0x00);
#endif
    __m128 c7_0 = _mm_load_ss(&C[(l_n*8)+0]);
    __m128 a7_0 = _mm_load_ss(&A[5]);
    c7_0 = _mm_add_ss(c7_0, _mm_mul_ss(a7_0, b7));
    _mm_store_ss(&C[(l_n*8)+0], c7_0);
    __m128 c7_1 = _mm_load_ss(&C[(l_n*8)+3]);
    __m128 a7_1 = _mm_load_ss(&A[6]);
    c7_1 = _mm_add_ss(c7_1, _mm_mul_ss(a7_1, b7));
    _mm_store_ss(&C[(l_n*8)+3], c7_1);
#else
    C[(l_n*8)+0] += A[5] * B[(l_n*16)+7];
    C[(l_n*8)+3] += A[6] * B[(l_n*16)+7];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#else
#endif

#if defined(__SSE3__) || defined(__AVX__)
#else
#endif

  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 126;
#endif
}

void ssparse_starMatrix_m4_n9_k9_ldA8_ldBna3_ldC8_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(4)
  #pragma vector aligned
  for ( l_m = 0; l_m < 4; l_m++) {
    C[0+l_m] += A[48+l_m] * B[0];
    C[0+l_m] += A[56+l_m] * B[1];
    C[0+l_m] += A[64+l_m] * B[2];
    C[8+l_m] += A[48+l_m] * B[3];
    C[8+l_m] += A[56+l_m] * B[4];
    C[8+l_m] += A[64+l_m] * B[5];
    C[16+l_m] += A[48+l_m] * B[6];
    C[16+l_m] += A[56+l_m] * B[7];
    C[16+l_m] += A[64+l_m] * B[8];
    C[24+l_m] += A[48+l_m] * B[9];
    C[24+l_m] += A[56+l_m] * B[10];
    C[32+l_m] += A[56+l_m] * B[11];
    C[32+l_m] += A[64+l_m] * B[12];
    C[40+l_m] += A[48+l_m] * B[13];
    C[40+l_m] += A[64+l_m] * B[14];
    C[48+l_m] += A[0+l_m] * B[15];
    C[48+l_m] += A[24+l_m] * B[16];
    C[48+l_m] += A[40+l_m] * B[17];
    C[56+l_m] += A[8+l_m] * B[18];
    C[56+l_m] += A[24+l_m] * B[19];
    C[56+l_m] += A[32+l_m] * B[20];
    C[64+l_m] += A[16+l_m] * B[21];
    C[64+l_m] += A[32+l_m] * B[22];
    C[64+l_m] += A[40+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 192;
#endif
}

void ssparse_kXiDivMT_m1_n9_k4_ldAna3_ldB8_ldC8_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_n = 0;
  #pragma nounroll_and_jam
  for ( l_n = 0; l_n < 9; l_n++) {
    unsigned int l_m = 0;
    for ( l_m = 0; l_m < 1; l_m++) {
      C[(l_n*8)+l_m] = 0.0f;
    }
#if defined(__SSE3__) || defined(__AVX__)
#else
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b1 = _mm_broadcast_ss(&B[(l_n*8)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b1 = _mm_load_ss(&B[(l_n*8)+1]);    b1 = _mm_shuffle_ps(b1, b1, 0x00);
#endif
    __m128 c1_0 = _mm_load_ss(&C[(l_n*8)+0]);
    __m128 a1_0 = _mm_load_ss(&A[0]);
    c1_0 = _mm_add_ss(c1_0, _mm_mul_ss(a1_0, b1));
    _mm_store_ss(&C[(l_n*8)+0], c1_0);
#else
    C[(l_n*8)+0] += A[0] * B[(l_n*8)+1];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#else
#endif

#if defined(__SSE3__) || defined(__AVX__)
#else
#endif

  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 18;
#endif
}

void ssparse_starMatrix_m1_n9_k9_ldA8_ldBna3_ldC8_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  for ( l_m = 0; l_m < 1; l_m++) {
    C[0+l_m] += A[48+l_m] * B[0];
    C[0+l_m] += A[56+l_m] * B[1];
    C[0+l_m] += A[64+l_m] * B[2];
    C[8+l_m] += A[48+l_m] * B[3];
    C[8+l_m] += A[56+l_m] * B[4];
    C[8+l_m] += A[64+l_m] * B[5];
    C[16+l_m] += A[48+l_m] * B[6];
    C[16+l_m] += A[56+l_m] * B[7];
    C[16+l_m] += A[64+l_m] * B[8];
    C[24+l_m] += A[48+l_m] * B[9];
    C[24+l_m] += A[56+l_m] * B[10];
    C[32+l_m] += A[56+l_m] * B[11];
    C[32+l_m] += A[64+l_m] * B[12];
    C[40+l_m] += A[48+l_m] * B[13];
    C[40+l_m] += A[64+l_m] * B[14];
    C[48+l_m] += A[0+l_m] * B[15];
    C[48+l_m] += A[24+l_m] * B[16];
    C[48+l_m] += A[40+l_m] * B[17];
    C[56+l_m] += A[8+l_m] * B[18];
    C[56+l_m] += A[24+l_m] * B[19];
    C[56+l_m] += A[32+l_m] * B[20];
    C[64+l_m] += A[16+l_m] * B[21];
    C[64+l_m] += A[32+l_m] * B[22];
    C[64+l_m] += A[40+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 48;
#endif
}

void ssparse_starMatrix_m10_n9_k9_ldA16_ldBna4_ldC16_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 10; l_m++) {
    C[0+l_m] += A[96+l_m] * B[0];
    C[0+l_m] += A[112+l_m] * B[1];
    C[0+l_m] += A[128+l_m] * B[2];
    C[16+l_m] += A[96+l_m] * B[3];
    C[16+l_m] += A[112+l_m] * B[4];
    C[16+l_m] += A[128+l_m] * B[5];
    C[32+l_m] += A[96+l_m] * B[6];
    C[32+l_m] += A[112+l_m] * B[7];
    C[32+l_m] += A[128+l_m] * B[8];
    C[48+l_m] += A[96+l_m] * B[9];
    C[48+l_m] += A[112+l_m] * B[10];
    C[64+l_m] += A[112+l_m] * B[11];
    C[64+l_m] += A[128+l_m] * B[12];
    C[80+l_m] += A[96+l_m] * B[13];
    C[80+l_m] += A[128+l_m] * B[14];
    C[96+l_m] += A[0+l_m] * B[15];
    C[96+l_m] += A[48+l_m] * B[16];
    C[96+l_m] += A[80+l_m] * B[17];
    C[112+l_m] += A[16+l_m] * B[18];
    C[112+l_m] += A[48+l_m] * B[19];
    C[112+l_m] += A[64+l_m] * B[20];
    C[128+l_m] += A[32+l_m] * B[21];
    C[128+l_m] += A[64+l_m] * B[22];
    C[128+l_m] += A[80+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 480;
#endif
}

void ssparse_starMatrix_m4_n9_k9_ldA8_ldBna4_ldC8_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(4)
  #pragma vector aligned
  for ( l_m = 0; l_m < 4; l_m++) {
    C[0+l_m] += A[48+l_m] * B[0];
    C[0+l_m] += A[56+l_m] * B[1];
    C[0+l_m] += A[64+l_m] * B[2];
    C[8+l_m] += A[48+l_m] * B[3];
    C[8+l_m] += A[56+l_m] * B[4];
    C[8+l_m] += A[64+l_m] * B[5];
    C[16+l_m] += A[48+l_m] * B[6];
    C[16+l_m] += A[56+l_m] * B[7];
    C[16+l_m] += A[64+l_m] * B[8];
    C[24+l_m] += A[48+l_m] * B[9];
    C[24+l_m] += A[56+l_m] * B[10];
    C[32+l_m] += A[56+l_m] * B[11];
    C[32+l_m] += A[64+l_m] * B[12];
    C[40+l_m] += A[48+l_m] * B[13];
    C[40+l_m] += A[64+l_m] * B[14];
    C[48+l_m] += A[0+l_m] * B[15];
    C[48+l_m] += A[24+l_m] * B[16];
    C[48+l_m] += A[40+l_m] * B[17];
    C[56+l_m] += A[8+l_m] * B[18];
    C[56+l_m] += A[24+l_m] * B[19];
    C[56+l_m] += A[32+l_m] * B[20];
    C[64+l_m] += A[16+l_m] * B[21];
    C[64+l_m] += A[32+l_m] * B[22];
    C[64+l_m] += A[40+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 192;
#endif
}

void ssparse_starMatrix_m1_n9_k9_ldA8_ldBna4_ldC8_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  for ( l_m = 0; l_m < 1; l_m++) {
    C[0+l_m] += A[48+l_m] * B[0];
    C[0+l_m] += A[56+l_m] * B[1];
    C[0+l_m] += A[64+l_m] * B[2];
    C[8+l_m] += A[48+l_m] * B[3];
    C[8+l_m] += A[56+l_m] * B[4];
    C[8+l_m] += A[64+l_m] * B[5];
    C[16+l_m] += A[48+l_m] * B[6];
    C[16+l_m] += A[56+l_m] * B[7];
    C[16+l_m] += A[64+l_m] * B[8];
    C[24+l_m] += A[48+l_m] * B[9];
    C[24+l_m] += A[56+l_m] * B[10];
    C[32+l_m] += A[56+l_m] * B[11];
    C[32+l_m] += A[64+l_m] * B[12];
    C[40+l_m] += A[48+l_m] * B[13];
    C[40+l_m] += A[64+l_m] * B[14];
    C[48+l_m] += A[0+l_m] * B[15];
    C[48+l_m] += A[24+l_m] * B[16];
    C[48+l_m] += A[40+l_m] * B[17];
    C[56+l_m] += A[8+l_m] * B[18];
    C[56+l_m] += A[24+l_m] * B[19];
    C[56+l_m] += A[32+l_m] * B[20];
    C[64+l_m] += A[16+l_m] * B[21];
    C[64+l_m] += A[32+l_m] * B[22];
    C[64+l_m] += A[40+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 48;
#endif
}

void ssparse_starMatrix_m20_n9_k9_ldA24_ldBna5_ldC24_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 20; l_m++) {
    C[0+l_m] += A[144+l_m] * B[0];
    C[0+l_m] += A[168+l_m] * B[1];
    C[0+l_m] += A[192+l_m] * B[2];
    C[24+l_m] += A[144+l_m] * B[3];
    C[24+l_m] += A[168+l_m] * B[4];
    C[24+l_m] += A[192+l_m] * B[5];
    C[48+l_m] += A[144+l_m] * B[6];
    C[48+l_m] += A[168+l_m] * B[7];
    C[48+l_m] += A[192+l_m] * B[8];
    C[72+l_m] += A[144+l_m] * B[9];
    C[72+l_m] += A[168+l_m] * B[10];
    C[96+l_m] += A[168+l_m] * B[11];
    C[96+l_m] += A[192+l_m] * B[12];
    C[120+l_m] += A[144+l_m] * B[13];
    C[120+l_m] += A[192+l_m] * B[14];
    C[144+l_m] += A[0+l_m] * B[15];
    C[144+l_m] += A[72+l_m] * B[16];
    C[144+l_m] += A[120+l_m] * B[17];
    C[168+l_m] += A[24+l_m] * B[18];
    C[168+l_m] += A[72+l_m] * B[19];
    C[168+l_m] += A[96+l_m] * B[20];
    C[192+l_m] += A[48+l_m] * B[21];
    C[192+l_m] += A[96+l_m] * B[22];
    C[192+l_m] += A[120+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 960;
#endif
}

void ssparse_starMatrix_m10_n9_k9_ldA16_ldBna5_ldC16_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 10; l_m++) {
    C[0+l_m] += A[96+l_m] * B[0];
    C[0+l_m] += A[112+l_m] * B[1];
    C[0+l_m] += A[128+l_m] * B[2];
    C[16+l_m] += A[96+l_m] * B[3];
    C[16+l_m] += A[112+l_m] * B[4];
    C[16+l_m] += A[128+l_m] * B[5];
    C[32+l_m] += A[96+l_m] * B[6];
    C[32+l_m] += A[112+l_m] * B[7];
    C[32+l_m] += A[128+l_m] * B[8];
    C[48+l_m] += A[96+l_m] * B[9];
    C[48+l_m] += A[112+l_m] * B[10];
    C[64+l_m] += A[112+l_m] * B[11];
    C[64+l_m] += A[128+l_m] * B[12];
    C[80+l_m] += A[96+l_m] * B[13];
    C[80+l_m] += A[128+l_m] * B[14];
    C[96+l_m] += A[0+l_m] * B[15];
    C[96+l_m] += A[48+l_m] * B[16];
    C[96+l_m] += A[80+l_m] * B[17];
    C[112+l_m] += A[16+l_m] * B[18];
    C[112+l_m] += A[48+l_m] * B[19];
    C[112+l_m] += A[64+l_m] * B[20];
    C[128+l_m] += A[32+l_m] * B[21];
    C[128+l_m] += A[64+l_m] * B[22];
    C[128+l_m] += A[80+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 480;
#endif
}

void ssparse_starMatrix_m4_n9_k9_ldA8_ldBna5_ldC8_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(4)
  #pragma vector aligned
  for ( l_m = 0; l_m < 4; l_m++) {
    C[0+l_m] += A[48+l_m] * B[0];
    C[0+l_m] += A[56+l_m] * B[1];
    C[0+l_m] += A[64+l_m] * B[2];
    C[8+l_m] += A[48+l_m] * B[3];
    C[8+l_m] += A[56+l_m] * B[4];
    C[8+l_m] += A[64+l_m] * B[5];
    C[16+l_m] += A[48+l_m] * B[6];
    C[16+l_m] += A[56+l_m] * B[7];
    C[16+l_m] += A[64+l_m] * B[8];
    C[24+l_m] += A[48+l_m] * B[9];
    C[24+l_m] += A[56+l_m] * B[10];
    C[32+l_m] += A[56+l_m] * B[11];
    C[32+l_m] += A[64+l_m] * B[12];
    C[40+l_m] += A[48+l_m] * B[13];
    C[40+l_m] += A[64+l_m] * B[14];
    C[48+l_m] += A[0+l_m] * B[15];
    C[48+l_m] += A[24+l_m] * B[16];
    C[48+l_m] += A[40+l_m] * B[17];
    C[56+l_m] += A[8+l_m] * B[18];
    C[56+l_m] += A[24+l_m] * B[19];
    C[56+l_m] += A[32+l_m] * B[20];
    C[64+l_m] += A[16+l_m] * B[21];
    C[64+l_m] += A[32+l_m] * B[22];
    C[64+l_m] += A[40+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 192;
#endif
}

void ssparse_starMatrix_m1_n9_k9_ldA8_ldBna5_ldC8_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  for ( l_m = 0; l_m < 1; l_m++) {
    C[0+l_m] += A[48+l_m] * B[0];
    C[0+l_m] += A[56+l_m] * B[1];
    C[0+l_m] += A[64+l_m] * B[2];
    C[8+l_m] += A[48+l_m] * B[3];
    C[8+l_m] += A[56+l_m] * B[4];
    C[8+l_m] += A[64+l_m] * B[5];
    C[16+l_m] += A[48+l_m] * B[6];
    C[16+l_m] += A[56+l_m] * B[7];
    C[16+l_m] += A[64+l_m] * B[8];
    C[24+l_m] += A[48+l_m] * B[9];
    C[24+l_m] += A[56+l_m] * B[10];
    C[32+l_m] += A[56+l_m] * B[11];
    C[32+l_m] += A[64+l_m] * B[12];
    C[40+l_m] += A[48+l_m] * B[13];
    C[40+l_m] += A[64+l_m] * B[14];
    C[48+l_m] += A[0+l_m] * B[15];
    C[48+l_m] += A[24+l_m] * B[16];
    C[48+l_m] += A[40+l_m] * B[17];
    C[56+l_m] += A[8+l_m] * B[18];
    C[56+l_m] += A[24+l_m] * B[19];
    C[56+l_m] += A[32+l_m] * B[20];
    C[64+l_m] += A[16+l_m] * B[21];
    C[64+l_m] += A[32+l_m] * B[22];
    C[64+l_m] += A[40+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 48;
#endif
}

void ssparse_starMatrix_m35_n9_k9_ldA40_ldBna6_ldC40_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 35; l_m++) {
    C[0+l_m] += A[240+l_m] * B[0];
    C[0+l_m] += A[280+l_m] * B[1];
    C[0+l_m] += A[320+l_m] * B[2];
    C[40+l_m] += A[240+l_m] * B[3];
    C[40+l_m] += A[280+l_m] * B[4];
    C[40+l_m] += A[320+l_m] * B[5];
    C[80+l_m] += A[240+l_m] * B[6];
    C[80+l_m] += A[280+l_m] * B[7];
    C[80+l_m] += A[320+l_m] * B[8];
    C[120+l_m] += A[240+l_m] * B[9];
    C[120+l_m] += A[280+l_m] * B[10];
    C[160+l_m] += A[280+l_m] * B[11];
    C[160+l_m] += A[320+l_m] * B[12];
    C[200+l_m] += A[240+l_m] * B[13];
    C[200+l_m] += A[320+l_m] * B[14];
    C[240+l_m] += A[0+l_m] * B[15];
    C[240+l_m] += A[120+l_m] * B[16];
    C[240+l_m] += A[200+l_m] * B[17];
    C[280+l_m] += A[40+l_m] * B[18];
    C[280+l_m] += A[120+l_m] * B[19];
    C[280+l_m] += A[160+l_m] * B[20];
    C[320+l_m] += A[80+l_m] * B[21];
    C[320+l_m] += A[160+l_m] * B[22];
    C[320+l_m] += A[200+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 1680;
#endif
}

void ssparse_starMatrix_m20_n9_k9_ldA24_ldBna6_ldC24_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 20; l_m++) {
    C[0+l_m] += A[144+l_m] * B[0];
    C[0+l_m] += A[168+l_m] * B[1];
    C[0+l_m] += A[192+l_m] * B[2];
    C[24+l_m] += A[144+l_m] * B[3];
    C[24+l_m] += A[168+l_m] * B[4];
    C[24+l_m] += A[192+l_m] * B[5];
    C[48+l_m] += A[144+l_m] * B[6];
    C[48+l_m] += A[168+l_m] * B[7];
    C[48+l_m] += A[192+l_m] * B[8];
    C[72+l_m] += A[144+l_m] * B[9];
    C[72+l_m] += A[168+l_m] * B[10];
    C[96+l_m] += A[168+l_m] * B[11];
    C[96+l_m] += A[192+l_m] * B[12];
    C[120+l_m] += A[144+l_m] * B[13];
    C[120+l_m] += A[192+l_m] * B[14];
    C[144+l_m] += A[0+l_m] * B[15];
    C[144+l_m] += A[72+l_m] * B[16];
    C[144+l_m] += A[120+l_m] * B[17];
    C[168+l_m] += A[24+l_m] * B[18];
    C[168+l_m] += A[72+l_m] * B[19];
    C[168+l_m] += A[96+l_m] * B[20];
    C[192+l_m] += A[48+l_m] * B[21];
    C[192+l_m] += A[96+l_m] * B[22];
    C[192+l_m] += A[120+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 960;
#endif
}

void ssparse_starMatrix_m10_n9_k9_ldA16_ldBna6_ldC16_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 10; l_m++) {
    C[0+l_m] += A[96+l_m] * B[0];
    C[0+l_m] += A[112+l_m] * B[1];
    C[0+l_m] += A[128+l_m] * B[2];
    C[16+l_m] += A[96+l_m] * B[3];
    C[16+l_m] += A[112+l_m] * B[4];
    C[16+l_m] += A[128+l_m] * B[5];
    C[32+l_m] += A[96+l_m] * B[6];
    C[32+l_m] += A[112+l_m] * B[7];
    C[32+l_m] += A[128+l_m] * B[8];
    C[48+l_m] += A[96+l_m] * B[9];
    C[48+l_m] += A[112+l_m] * B[10];
    C[64+l_m] += A[112+l_m] * B[11];
    C[64+l_m] += A[128+l_m] * B[12];
    C[80+l_m] += A[96+l_m] * B[13];
    C[80+l_m] += A[128+l_m] * B[14];
    C[96+l_m] += A[0+l_m] * B[15];
    C[96+l_m] += A[48+l_m] * B[16];
    C[96+l_m] += A[80+l_m] * B[17];
    C[112+l_m] += A[16+l_m] * B[18];
    C[112+l_m] += A[48+l_m] * B[19];
    C[112+l_m] += A[64+l_m] * B[20];
    C[128+l_m] += A[32+l_m] * B[21];
    C[128+l_m] += A[64+l_m] * B[22];
    C[128+l_m] += A[80+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 480;
#endif
}

void ssparse_starMatrix_m4_n9_k9_ldA8_ldBna6_ldC8_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(4)
  #pragma vector aligned
  for ( l_m = 0; l_m < 4; l_m++) {
    C[0+l_m] += A[48+l_m] * B[0];
    C[0+l_m] += A[56+l_m] * B[1];
    C[0+l_m] += A[64+l_m] * B[2];
    C[8+l_m] += A[48+l_m] * B[3];
    C[8+l_m] += A[56+l_m] * B[4];
    C[8+l_m] += A[64+l_m] * B[5];
    C[16+l_m] += A[48+l_m] * B[6];
    C[16+l_m] += A[56+l_m] * B[7];
    C[16+l_m] += A[64+l_m] * B[8];
    C[24+l_m] += A[48+l_m] * B[9];
    C[24+l_m] += A[56+l_m] * B[10];
    C[32+l_m] += A[56+l_m] * B[11];
    C[32+l_m] += A[64+l_m] * B[12];
    C[40+l_m] += A[48+l_m] * B[13];
    C[40+l_m] += A[64+l_m] * B[14];
    C[48+l_m] += A[0+l_m] * B[15];
    C[48+l_m] += A[24+l_m] * B[16];
    C[48+l_m] += A[40+l_m] * B[17];
    C[56+l_m] += A[8+l_m] * B[18];
    C[56+l_m] += A[24+l_m] * B[19];
    C[56+l_m] += A[32+l_m] * B[20];
    C[64+l_m] += A[16+l_m] * B[21];
    C[64+l_m] += A[32+l_m] * B[22];
    C[64+l_m] += A[40+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 192;
#endif
}

void ssparse_starMatrix_m1_n9_k9_ldA8_ldBna6_ldC8_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  for ( l_m = 0; l_m < 1; l_m++) {
    C[0+l_m] += A[48+l_m] * B[0];
    C[0+l_m] += A[56+l_m] * B[1];
    C[0+l_m] += A[64+l_m] * B[2];
    C[8+l_m] += A[48+l_m] * B[3];
    C[8+l_m] += A[56+l_m] * B[4];
    C[8+l_m] += A[64+l_m] * B[5];
    C[16+l_m] += A[48+l_m] * B[6];
    C[16+l_m] += A[56+l_m] * B[7];
    C[16+l_m] += A[64+l_m] * B[8];
    C[24+l_m] += A[48+l_m] * B[9];
    C[24+l_m] += A[56+l_m] * B[10];
    C[32+l_m] += A[56+l_m] * B[11];
    C[32+l_m] += A[64+l_m] * B[12];
    C[40+l_m] += A[48+l_m] * B[13];
    C[40+l_m] += A[64+l_m] * B[14];
    C[48+l_m] += A[0+l_m] * B[15];
    C[48+l_m] += A[24+l_m] * B[16];
    C[48+l_m] += A[40+l_m] * B[17];
    C[56+l_m] += A[8+l_m] * B[18];
    C[56+l_m] += A[24+l_m] * B[19];
    C[56+l_m] += A[32+l_m] * B[20];
    C[64+l_m] += A[16+l_m] * B[21];
    C[64+l_m] += A[32+l_m] * B[22];
    C[64+l_m] += A[40+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 48;
#endif
}

void ssparse_starMatrix_m56_n9_k9_ldA56_ldBna7_ldC56_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 56; l_m++) {
    C[0+l_m] += A[336+l_m] * B[0];
    C[0+l_m] += A[392+l_m] * B[1];
    C[0+l_m] += A[448+l_m] * B[2];
    C[56+l_m] += A[336+l_m] * B[3];
    C[56+l_m] += A[392+l_m] * B[4];
    C[56+l_m] += A[448+l_m] * B[5];
    C[112+l_m] += A[336+l_m] * B[6];
    C[112+l_m] += A[392+l_m] * B[7];
    C[112+l_m] += A[448+l_m] * B[8];
    C[168+l_m] += A[336+l_m] * B[9];
    C[168+l_m] += A[392+l_m] * B[10];
    C[224+l_m] += A[392+l_m] * B[11];
    C[224+l_m] += A[448+l_m] * B[12];
    C[280+l_m] += A[336+l_m] * B[13];
    C[280+l_m] += A[448+l_m] * B[14];
    C[336+l_m] += A[0+l_m] * B[15];
    C[336+l_m] += A[168+l_m] * B[16];
    C[336+l_m] += A[280+l_m] * B[17];
    C[392+l_m] += A[56+l_m] * B[18];
    C[392+l_m] += A[168+l_m] * B[19];
    C[392+l_m] += A[224+l_m] * B[20];
    C[448+l_m] += A[112+l_m] * B[21];
    C[448+l_m] += A[224+l_m] * B[22];
    C[448+l_m] += A[280+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 2688;
#endif
}

void ssparse_starMatrix_m35_n9_k9_ldA40_ldBna7_ldC40_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 35; l_m++) {
    C[0+l_m] += A[240+l_m] * B[0];
    C[0+l_m] += A[280+l_m] * B[1];
    C[0+l_m] += A[320+l_m] * B[2];
    C[40+l_m] += A[240+l_m] * B[3];
    C[40+l_m] += A[280+l_m] * B[4];
    C[40+l_m] += A[320+l_m] * B[5];
    C[80+l_m] += A[240+l_m] * B[6];
    C[80+l_m] += A[280+l_m] * B[7];
    C[80+l_m] += A[320+l_m] * B[8];
    C[120+l_m] += A[240+l_m] * B[9];
    C[120+l_m] += A[280+l_m] * B[10];
    C[160+l_m] += A[280+l_m] * B[11];
    C[160+l_m] += A[320+l_m] * B[12];
    C[200+l_m] += A[240+l_m] * B[13];
    C[200+l_m] += A[320+l_m] * B[14];
    C[240+l_m] += A[0+l_m] * B[15];
    C[240+l_m] += A[120+l_m] * B[16];
    C[240+l_m] += A[200+l_m] * B[17];
    C[280+l_m] += A[40+l_m] * B[18];
    C[280+l_m] += A[120+l_m] * B[19];
    C[280+l_m] += A[160+l_m] * B[20];
    C[320+l_m] += A[80+l_m] * B[21];
    C[320+l_m] += A[160+l_m] * B[22];
    C[320+l_m] += A[200+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 1680;
#endif
}

void ssparse_starMatrix_m20_n9_k9_ldA24_ldBna7_ldC24_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 20; l_m++) {
    C[0+l_m] += A[144+l_m] * B[0];
    C[0+l_m] += A[168+l_m] * B[1];
    C[0+l_m] += A[192+l_m] * B[2];
    C[24+l_m] += A[144+l_m] * B[3];
    C[24+l_m] += A[168+l_m] * B[4];
    C[24+l_m] += A[192+l_m] * B[5];
    C[48+l_m] += A[144+l_m] * B[6];
    C[48+l_m] += A[168+l_m] * B[7];
    C[48+l_m] += A[192+l_m] * B[8];
    C[72+l_m] += A[144+l_m] * B[9];
    C[72+l_m] += A[168+l_m] * B[10];
    C[96+l_m] += A[168+l_m] * B[11];
    C[96+l_m] += A[192+l_m] * B[12];
    C[120+l_m] += A[144+l_m] * B[13];
    C[120+l_m] += A[192+l_m] * B[14];
    C[144+l_m] += A[0+l_m] * B[15];
    C[144+l_m] += A[72+l_m] * B[16];
    C[144+l_m] += A[120+l_m] * B[17];
    C[168+l_m] += A[24+l_m] * B[18];
    C[168+l_m] += A[72+l_m] * B[19];
    C[168+l_m] += A[96+l_m] * B[20];
    C[192+l_m] += A[48+l_m] * B[21];
    C[192+l_m] += A[96+l_m] * B[22];
    C[192+l_m] += A[120+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 960;
#endif
}

void ssparse_starMatrix_m10_n9_k9_ldA16_ldBna7_ldC16_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 10; l_m++) {
    C[0+l_m] += A[96+l_m] * B[0];
    C[0+l_m] += A[112+l_m] * B[1];
    C[0+l_m] += A[128+l_m] * B[2];
    C[16+l_m] += A[96+l_m] * B[3];
    C[16+l_m] += A[112+l_m] * B[4];
    C[16+l_m] += A[128+l_m] * B[5];
    C[32+l_m] += A[96+l_m] * B[6];
    C[32+l_m] += A[112+l_m] * B[7];
    C[32+l_m] += A[128+l_m] * B[8];
    C[48+l_m] += A[96+l_m] * B[9];
    C[48+l_m] += A[112+l_m] * B[10];
    C[64+l_m] += A[112+l_m] * B[11];
    C[64+l_m] += A[128+l_m] * B[12];
    C[80+l_m] += A[96+l_m] * B[13];
    C[80+l_m] += A[128+l_m] * B[14];
    C[96+l_m] += A[0+l_m] * B[15];
    C[96+l_m] += A[48+l_m] * B[16];
    C[96+l_m] += A[80+l_m] * B[17];
    C[112+l_m] += A[16+l_m] * B[18];
    C[112+l_m] += A[48+l_m] * B[19];
    C[112+l_m] += A[64+l_m] * B[20];
    C[128+l_m] += A[32+l_m] * B[21];
    C[128+l_m] += A[64+l_m] * B[22];
    C[128+l_m] += A[80+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 480;
#endif
}

void ssparse_starMatrix_m4_n9_k9_ldA8_ldBna7_ldC8_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(4)
  #pragma vector aligned
  for ( l_m = 0; l_m < 4; l_m++) {
    C[0+l_m] += A[48+l_m] * B[0];
    C[0+l_m] += A[56+l_m] * B[1];
    C[0+l_m] += A[64+l_m] * B[2];
    C[8+l_m] += A[48+l_m] * B[3];
    C[8+l_m] += A[56+l_m] * B[4];
    C[8+l_m] += A[64+l_m] * B[5];
    C[16+l_m] += A[48+l_m] * B[6];
    C[16+l_m] += A[56+l_m] * B[7];
    C[16+l_m] += A[64+l_m] * B[8];
    C[24+l_m] += A[48+l_m] * B[9];
    C[24+l_m] += A[56+l_m] * B[10];
    C[32+l_m] += A[56+l_m] * B[11];
    C[32+l_m] += A[64+l_m] * B[12];
    C[40+l_m] += A[48+l_m] * B[13];
    C[40+l_m] += A[64+l_m] * B[14];
    C[48+l_m] += A[0+l_m] * B[15];
    C[48+l_m] += A[24+l_m] * B[16];
    C[48+l_m] += A[40+l_m] * B[17];
    C[56+l_m] += A[8+l_m] * B[18];
    C[56+l_m] += A[24+l_m] * B[19];
    C[56+l_m] += A[32+l_m] * B[20];
    C[64+l_m] += A[16+l_m] * B[21];
    C[64+l_m] += A[32+l_m] * B[22];
    C[64+l_m] += A[40+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 192;
#endif
}

void ssparse_starMatrix_m1_n9_k9_ldA8_ldBna7_ldC8_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  for ( l_m = 0; l_m < 1; l_m++) {
    C[0+l_m] += A[48+l_m] * B[0];
    C[0+l_m] += A[56+l_m] * B[1];
    C[0+l_m] += A[64+l_m] * B[2];
    C[8+l_m] += A[48+l_m] * B[3];
    C[8+l_m] += A[56+l_m] * B[4];
    C[8+l_m] += A[64+l_m] * B[5];
    C[16+l_m] += A[48+l_m] * B[6];
    C[16+l_m] += A[56+l_m] * B[7];
    C[16+l_m] += A[64+l_m] * B[8];
    C[24+l_m] += A[48+l_m] * B[9];
    C[24+l_m] += A[56+l_m] * B[10];
    C[32+l_m] += A[56+l_m] * B[11];
    C[32+l_m] += A[64+l_m] * B[12];
    C[40+l_m] += A[48+l_m] * B[13];
    C[40+l_m] += A[64+l_m] * B[14];
    C[48+l_m] += A[0+l_m] * B[15];
    C[48+l_m] += A[24+l_m] * B[16];
    C[48+l_m] += A[40+l_m] * B[17];
    C[56+l_m] += A[8+l_m] * B[18];
    C[56+l_m] += A[24+l_m] * B[19];
    C[56+l_m] += A[32+l_m] * B[20];
    C[64+l_m] += A[16+l_m] * B[21];
    C[64+l_m] += A[32+l_m] * B[22];
    C[64+l_m] += A[40+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 48;
#endif
}

void ssparse_starMatrix_m84_n9_k9_ldA88_ldBna8_ldC88_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 84; l_m++) {
    C[0+l_m] += A[528+l_m] * B[0];
    C[0+l_m] += A[616+l_m] * B[1];
    C[0+l_m] += A[704+l_m] * B[2];
    C[88+l_m] += A[528+l_m] * B[3];
    C[88+l_m] += A[616+l_m] * B[4];
    C[88+l_m] += A[704+l_m] * B[5];
    C[176+l_m] += A[528+l_m] * B[6];
    C[176+l_m] += A[616+l_m] * B[7];
    C[176+l_m] += A[704+l_m] * B[8];
    C[264+l_m] += A[528+l_m] * B[9];
    C[264+l_m] += A[616+l_m] * B[10];
    C[352+l_m] += A[616+l_m] * B[11];
    C[352+l_m] += A[704+l_m] * B[12];
    C[440+l_m] += A[528+l_m] * B[13];
    C[440+l_m] += A[704+l_m] * B[14];
    C[528+l_m] += A[0+l_m] * B[15];
    C[528+l_m] += A[264+l_m] * B[16];
    C[528+l_m] += A[440+l_m] * B[17];
    C[616+l_m] += A[88+l_m] * B[18];
    C[616+l_m] += A[264+l_m] * B[19];
    C[616+l_m] += A[352+l_m] * B[20];
    C[704+l_m] += A[176+l_m] * B[21];
    C[704+l_m] += A[352+l_m] * B[22];
    C[704+l_m] += A[440+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 4032;
#endif
}

void ssparse_starMatrix_m56_n9_k9_ldA56_ldBna8_ldC56_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 56; l_m++) {
    C[0+l_m] += A[336+l_m] * B[0];
    C[0+l_m] += A[392+l_m] * B[1];
    C[0+l_m] += A[448+l_m] * B[2];
    C[56+l_m] += A[336+l_m] * B[3];
    C[56+l_m] += A[392+l_m] * B[4];
    C[56+l_m] += A[448+l_m] * B[5];
    C[112+l_m] += A[336+l_m] * B[6];
    C[112+l_m] += A[392+l_m] * B[7];
    C[112+l_m] += A[448+l_m] * B[8];
    C[168+l_m] += A[336+l_m] * B[9];
    C[168+l_m] += A[392+l_m] * B[10];
    C[224+l_m] += A[392+l_m] * B[11];
    C[224+l_m] += A[448+l_m] * B[12];
    C[280+l_m] += A[336+l_m] * B[13];
    C[280+l_m] += A[448+l_m] * B[14];
    C[336+l_m] += A[0+l_m] * B[15];
    C[336+l_m] += A[168+l_m] * B[16];
    C[336+l_m] += A[280+l_m] * B[17];
    C[392+l_m] += A[56+l_m] * B[18];
    C[392+l_m] += A[168+l_m] * B[19];
    C[392+l_m] += A[224+l_m] * B[20];
    C[448+l_m] += A[112+l_m] * B[21];
    C[448+l_m] += A[224+l_m] * B[22];
    C[448+l_m] += A[280+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 2688;
#endif
}

void ssparse_starMatrix_m35_n9_k9_ldA40_ldBna8_ldC40_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 35; l_m++) {
    C[0+l_m] += A[240+l_m] * B[0];
    C[0+l_m] += A[280+l_m] * B[1];
    C[0+l_m] += A[320+l_m] * B[2];
    C[40+l_m] += A[240+l_m] * B[3];
    C[40+l_m] += A[280+l_m] * B[4];
    C[40+l_m] += A[320+l_m] * B[5];
    C[80+l_m] += A[240+l_m] * B[6];
    C[80+l_m] += A[280+l_m] * B[7];
    C[80+l_m] += A[320+l_m] * B[8];
    C[120+l_m] += A[240+l_m] * B[9];
    C[120+l_m] += A[280+l_m] * B[10];
    C[160+l_m] += A[280+l_m] * B[11];
    C[160+l_m] += A[320+l_m] * B[12];
    C[200+l_m] += A[240+l_m] * B[13];
    C[200+l_m] += A[320+l_m] * B[14];
    C[240+l_m] += A[0+l_m] * B[15];
    C[240+l_m] += A[120+l_m] * B[16];
    C[240+l_m] += A[200+l_m] * B[17];
    C[280+l_m] += A[40+l_m] * B[18];
    C[280+l_m] += A[120+l_m] * B[19];
    C[280+l_m] += A[160+l_m] * B[20];
    C[320+l_m] += A[80+l_m] * B[21];
    C[320+l_m] += A[160+l_m] * B[22];
    C[320+l_m] += A[200+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 1680;
#endif
}

void ssparse_starMatrix_m20_n9_k9_ldA24_ldBna8_ldC24_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 20; l_m++) {
    C[0+l_m] += A[144+l_m] * B[0];
    C[0+l_m] += A[168+l_m] * B[1];
    C[0+l_m] += A[192+l_m] * B[2];
    C[24+l_m] += A[144+l_m] * B[3];
    C[24+l_m] += A[168+l_m] * B[4];
    C[24+l_m] += A[192+l_m] * B[5];
    C[48+l_m] += A[144+l_m] * B[6];
    C[48+l_m] += A[168+l_m] * B[7];
    C[48+l_m] += A[192+l_m] * B[8];
    C[72+l_m] += A[144+l_m] * B[9];
    C[72+l_m] += A[168+l_m] * B[10];
    C[96+l_m] += A[168+l_m] * B[11];
    C[96+l_m] += A[192+l_m] * B[12];
    C[120+l_m] += A[144+l_m] * B[13];
    C[120+l_m] += A[192+l_m] * B[14];
    C[144+l_m] += A[0+l_m] * B[15];
    C[144+l_m] += A[72+l_m] * B[16];
    C[144+l_m] += A[120+l_m] * B[17];
    C[168+l_m] += A[24+l_m] * B[18];
    C[168+l_m] += A[72+l_m] * B[19];
    C[168+l_m] += A[96+l_m] * B[20];
    C[192+l_m] += A[48+l_m] * B[21];
    C[192+l_m] += A[96+l_m] * B[22];
    C[192+l_m] += A[120+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 960;
#endif
}

void ssparse_starMatrix_m10_n9_k9_ldA16_ldBna8_ldC16_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 10; l_m++) {
    C[0+l_m] += A[96+l_m] * B[0];
    C[0+l_m] += A[112+l_m] * B[1];
    C[0+l_m] += A[128+l_m] * B[2];
    C[16+l_m] += A[96+l_m] * B[3];
    C[16+l_m] += A[112+l_m] * B[4];
    C[16+l_m] += A[128+l_m] * B[5];
    C[32+l_m] += A[96+l_m] * B[6];
    C[32+l_m] += A[112+l_m] * B[7];
    C[32+l_m] += A[128+l_m] * B[8];
    C[48+l_m] += A[96+l_m] * B[9];
    C[48+l_m] += A[112+l_m] * B[10];
    C[64+l_m] += A[112+l_m] * B[11];
    C[64+l_m] += A[128+l_m] * B[12];
    C[80+l_m] += A[96+l_m] * B[13];
    C[80+l_m] += A[128+l_m] * B[14];
    C[96+l_m] += A[0+l_m] * B[15];
    C[96+l_m] += A[48+l_m] * B[16];
    C[96+l_m] += A[80+l_m] * B[17];
    C[112+l_m] += A[16+l_m] * B[18];
    C[112+l_m] += A[48+l_m] * B[19];
    C[112+l_m] += A[64+l_m] * B[20];
    C[128+l_m] += A[32+l_m] * B[21];
    C[128+l_m] += A[64+l_m] * B[22];
    C[128+l_m] += A[80+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 480;
#endif
}

void ssparse_starMatrix_m4_n9_k9_ldA8_ldBna8_ldC8_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(4)
  #pragma vector aligned
  for ( l_m = 0; l_m < 4; l_m++) {
    C[0+l_m] += A[48+l_m] * B[0];
    C[0+l_m] += A[56+l_m] * B[1];
    C[0+l_m] += A[64+l_m] * B[2];
    C[8+l_m] += A[48+l_m] * B[3];
    C[8+l_m] += A[56+l_m] * B[4];
    C[8+l_m] += A[64+l_m] * B[5];
    C[16+l_m] += A[48+l_m] * B[6];
    C[16+l_m] += A[56+l_m] * B[7];
    C[16+l_m] += A[64+l_m] * B[8];
    C[24+l_m] += A[48+l_m] * B[9];
    C[24+l_m] += A[56+l_m] * B[10];
    C[32+l_m] += A[56+l_m] * B[11];
    C[32+l_m] += A[64+l_m] * B[12];
    C[40+l_m] += A[48+l_m] * B[13];
    C[40+l_m] += A[64+l_m] * B[14];
    C[48+l_m] += A[0+l_m] * B[15];
    C[48+l_m] += A[24+l_m] * B[16];
    C[48+l_m] += A[40+l_m] * B[17];
    C[56+l_m] += A[8+l_m] * B[18];
    C[56+l_m] += A[24+l_m] * B[19];
    C[56+l_m] += A[32+l_m] * B[20];
    C[64+l_m] += A[16+l_m] * B[21];
    C[64+l_m] += A[32+l_m] * B[22];
    C[64+l_m] += A[40+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 192;
#endif
}

void ssparse_starMatrix_m1_n9_k9_ldA8_ldBna8_ldC8_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  for ( l_m = 0; l_m < 1; l_m++) {
    C[0+l_m] += A[48+l_m] * B[0];
    C[0+l_m] += A[56+l_m] * B[1];
    C[0+l_m] += A[64+l_m] * B[2];
    C[8+l_m] += A[48+l_m] * B[3];
    C[8+l_m] += A[56+l_m] * B[4];
    C[8+l_m] += A[64+l_m] * B[5];
    C[16+l_m] += A[48+l_m] * B[6];
    C[16+l_m] += A[56+l_m] * B[7];
    C[16+l_m] += A[64+l_m] * B[8];
    C[24+l_m] += A[48+l_m] * B[9];
    C[24+l_m] += A[56+l_m] * B[10];
    C[32+l_m] += A[56+l_m] * B[11];
    C[32+l_m] += A[64+l_m] * B[12];
    C[40+l_m] += A[48+l_m] * B[13];
    C[40+l_m] += A[64+l_m] * B[14];
    C[48+l_m] += A[0+l_m] * B[15];
    C[48+l_m] += A[24+l_m] * B[16];
    C[48+l_m] += A[40+l_m] * B[17];
    C[56+l_m] += A[8+l_m] * B[18];
    C[56+l_m] += A[24+l_m] * B[19];
    C[56+l_m] += A[32+l_m] * B[20];
    C[64+l_m] += A[16+l_m] * B[21];
    C[64+l_m] += A[32+l_m] * B[22];
    C[64+l_m] += A[40+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 48;
#endif
}

void ssparse_kXiDivM_m4_n9_k4_ldAna2_ldB8_ldC8_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_n = 0;
  #pragma nounroll_and_jam
  for ( l_n = 0; l_n < 9; l_n++) {
    unsigned int l_m = 0;
   #pragma simd
    for ( l_m = 0; l_m < 4; l_m++) {
      C[(l_n*8)+l_m] = 0.0f;
    }
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b0 = _mm_broadcast_ss(&B[(l_n*8)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b0 = _mm_load_ss(&B[(l_n*8)+0]);    b0 = _mm_shuffle_ps(b0, b0, 0x00);
#endif
    __m128 c0_0 = _mm_load_ss(&C[(l_n*8)+1]);
    __m128 a0_0 = _mm_load_ss(&A[0]);
    c0_0 = _mm_add_ss(c0_0, _mm_mul_ss(a0_0, b0));
    _mm_store_ss(&C[(l_n*8)+1], c0_0);
#else
    C[(l_n*8)+1] += A[0] * B[(l_n*8)+0];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#else
#endif

#if defined(__SSE3__) || defined(__AVX__)
#else
#endif

#if defined(__SSE3__) || defined(__AVX__)
#else
#endif

  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 18;
#endif
}

void ssparse_starMatrix_m4_n9_k9_ldA8_ldBna2_ldC8_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(4)
  #pragma vector aligned
  for ( l_m = 0; l_m < 4; l_m++) {
    C[0+l_m] += A[48+l_m] * B[0];
    C[0+l_m] += A[56+l_m] * B[1];
    C[0+l_m] += A[64+l_m] * B[2];
    C[8+l_m] += A[48+l_m] * B[3];
    C[8+l_m] += A[56+l_m] * B[4];
    C[8+l_m] += A[64+l_m] * B[5];
    C[16+l_m] += A[48+l_m] * B[6];
    C[16+l_m] += A[56+l_m] * B[7];
    C[16+l_m] += A[64+l_m] * B[8];
    C[24+l_m] += A[48+l_m] * B[9];
    C[24+l_m] += A[56+l_m] * B[10];
    C[32+l_m] += A[56+l_m] * B[11];
    C[32+l_m] += A[64+l_m] * B[12];
    C[40+l_m] += A[48+l_m] * B[13];
    C[40+l_m] += A[64+l_m] * B[14];
    C[48+l_m] += A[0+l_m] * B[15];
    C[48+l_m] += A[24+l_m] * B[16];
    C[48+l_m] += A[40+l_m] * B[17];
    C[56+l_m] += A[8+l_m] * B[18];
    C[56+l_m] += A[24+l_m] * B[19];
    C[56+l_m] += A[32+l_m] * B[20];
    C[64+l_m] += A[16+l_m] * B[21];
    C[64+l_m] += A[32+l_m] * B[22];
    C[64+l_m] += A[40+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 192;
#endif
}

void ssparse_fP113DivM_m4_n9_k4_ldAna2_ldB8_ldC8_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_n = 0;
  #pragma nounroll_and_jam
  for ( l_n = 0; l_n < 9; l_n++) {
    unsigned int l_m = 0;
   #pragma simd
    for ( l_m = 0; l_m < 4; l_m++) {
      C[(l_n*8)+l_m] = 0.0f;
    }
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b0 = _mm_broadcast_ss(&B[(l_n*8)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b0 = _mm_load_ss(&B[(l_n*8)+0]);    b0 = _mm_shuffle_ps(b0, b0, 0x00);
#endif
    __m128 c0_0 = _mm_load_ss(&C[(l_n*8)+0]);
    __m128 a0_0 = _mm_load_ss(&A[0]);
    c0_0 = _mm_add_ss(c0_0, _mm_mul_ss(a0_0, b0));
    _mm_store_ss(&C[(l_n*8)+0], c0_0);
    __m128 c0_1 = _mm_load_ss(&C[(l_n*8)+3]);
    __m128 a0_1 = _mm_load_ss(&A[1]);
    c0_1 = _mm_add_ss(c0_1, _mm_mul_ss(a0_1, b0));
    _mm_store_ss(&C[(l_n*8)+3], c0_1);
#else
    C[(l_n*8)+0] += A[0] * B[(l_n*8)+0];
    C[(l_n*8)+3] += A[1] * B[(l_n*8)+0];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b1 = _mm_broadcast_ss(&B[(l_n*8)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b1 = _mm_load_ss(&B[(l_n*8)+1]);    b1 = _mm_shuffle_ps(b1, b1, 0x00);
#endif
    __m128 c1_0 = _mm_load_ss(&C[(l_n*8)+1]);
    __m128 a1_0 = _mm_load_ss(&A[2]);
    c1_0 = _mm_add_ss(c1_0, _mm_mul_ss(a1_0, b1));
    _mm_store_ss(&C[(l_n*8)+1], c1_0);
#else
    C[(l_n*8)+1] += A[2] * B[(l_n*8)+1];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b2 = _mm_broadcast_ss(&B[(l_n*8)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b2 = _mm_load_ss(&B[(l_n*8)+2]);    b2 = _mm_shuffle_ps(b2, b2, 0x00);
#endif
    __m128 c2_0 = _mm_load_ss(&C[(l_n*8)+2]);
    __m128 a2_0 = _mm_load_ss(&A[3]);
    c2_0 = _mm_add_ss(c2_0, _mm_mul_ss(a2_0, b2));
    _mm_store_ss(&C[(l_n*8)+2], c2_0);
#else
    C[(l_n*8)+2] += A[3] * B[(l_n*8)+2];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b3 = _mm_broadcast_ss(&B[(l_n*8)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b3 = _mm_load_ss(&B[(l_n*8)+3]);    b3 = _mm_shuffle_ps(b3, b3, 0x00);
#endif
    __m128 c3_0 = _mm_load_ss(&C[(l_n*8)+0]);
    __m128 a3_0 = _mm_load_ss(&A[4]);
    c3_0 = _mm_add_ss(c3_0, _mm_mul_ss(a3_0, b3));
    _mm_store_ss(&C[(l_n*8)+0], c3_0);
    __m128 c3_1 = _mm_load_ss(&C[(l_n*8)+3]);
    __m128 a3_1 = _mm_load_ss(&A[5]);
    c3_1 = _mm_add_ss(c3_1, _mm_mul_ss(a3_1, b3));
    _mm_store_ss(&C[(l_n*8)+3], c3_1);
#else
    C[(l_n*8)+0] += A[4] * B[(l_n*8)+3];
    C[(l_n*8)+3] += A[5] * B[(l_n*8)+3];
#endif

  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 108;
#endif
}

void ssparse_fP111DivM_m4_n9_k4_ldAna2_ldB8_ldC8_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_n = 0;
  #pragma nounroll_and_jam
  for ( l_n = 0; l_n < 9; l_n++) {
    unsigned int l_m = 0;
   #pragma simd
    for ( l_m = 0; l_m < 4; l_m++) {
      C[(l_n*8)+l_m] = 0.0f;
    }
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b0 = _mm_broadcast_ss(&B[(l_n*8)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b0 = _mm_load_ss(&B[(l_n*8)+0]);    b0 = _mm_shuffle_ps(b0, b0, 0x00);
#endif
    __m128 c0_0 = _mm_load_ss(&C[(l_n*8)+0]);
    __m128 a0_0 = _mm_load_ss(&A[0]);
    c0_0 = _mm_add_ss(c0_0, _mm_mul_ss(a0_0, b0));
    _mm_store_ss(&C[(l_n*8)+0], c0_0);
    __m128 c0_1 = _mm_load_ss(&C[(l_n*8)+3]);
    __m128 a0_1 = _mm_load_ss(&A[1]);
    c0_1 = _mm_add_ss(c0_1, _mm_mul_ss(a0_1, b0));
    _mm_store_ss(&C[(l_n*8)+3], c0_1);
#else
    C[(l_n*8)+0] += A[0] * B[(l_n*8)+0];
    C[(l_n*8)+3] += A[1] * B[(l_n*8)+0];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b1 = _mm_broadcast_ss(&B[(l_n*8)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b1 = _mm_load_ss(&B[(l_n*8)+1]);    b1 = _mm_shuffle_ps(b1, b1, 0x00);
#endif
    __m128 c1_0 = _mm_castpd_ps(_mm_load_sd((const double*)&C[(l_n*8)+1]));
    __m128 a1_0 = _mm_castpd_ps(_mm_load_sd((const double*)&A[2]));
    c1_0 = _mm_add_ps(c1_0, _mm_mul_ps(a1_0, b1));
    _mm_store_sd((double*)&C[(l_n*8)+1], _mm_castps_pd(c1_0));
#else
    C[(l_n*8)+1] += A[2] * B[(l_n*8)+1];
    C[(l_n*8)+2] += A[3] * B[(l_n*8)+1];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b2 = _mm_broadcast_ss(&B[(l_n*8)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b2 = _mm_load_ss(&B[(l_n*8)+2]);    b2 = _mm_shuffle_ps(b2, b2, 0x00);
#endif
    __m128 c2_0 = _mm_castpd_ps(_mm_load_sd((const double*)&C[(l_n*8)+1]));
    __m128 a2_0 = _mm_castpd_ps(_mm_load_sd((const double*)&A[4]));
    c2_0 = _mm_add_ps(c2_0, _mm_mul_ps(a2_0, b2));
    _mm_store_sd((double*)&C[(l_n*8)+1], _mm_castps_pd(c2_0));
#else
    C[(l_n*8)+1] += A[4] * B[(l_n*8)+2];
    C[(l_n*8)+2] += A[5] * B[(l_n*8)+2];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b3 = _mm_broadcast_ss(&B[(l_n*8)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b3 = _mm_load_ss(&B[(l_n*8)+3]);    b3 = _mm_shuffle_ps(b3, b3, 0x00);
#endif
    __m128 c3_0 = _mm_load_ss(&C[(l_n*8)+0]);
    __m128 a3_0 = _mm_load_ss(&A[6]);
    c3_0 = _mm_add_ss(c3_0, _mm_mul_ss(a3_0, b3));
    _mm_store_ss(&C[(l_n*8)+0], c3_0);
    __m128 c3_1 = _mm_load_ss(&C[(l_n*8)+3]);
    __m128 a3_1 = _mm_load_ss(&A[7]);
    c3_1 = _mm_add_ss(c3_1, _mm_mul_ss(a3_1, b3));
    _mm_store_ss(&C[(l_n*8)+3], c3_1);
#else
    C[(l_n*8)+0] += A[6] * B[(l_n*8)+3];
    C[(l_n*8)+3] += A[7] * B[(l_n*8)+3];
#endif

  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 144;
#endif
}

void ssparse_starMatrix_m10_n9_k9_ldA16_ldBna3_ldC16_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 10; l_m++) {
    C[0+l_m] += A[96+l_m] * B[0];
    C[0+l_m] += A[112+l_m] * B[1];
    C[0+l_m] += A[128+l_m] * B[2];
    C[16+l_m] += A[96+l_m] * B[3];
    C[16+l_m] += A[112+l_m] * B[4];
    C[16+l_m] += A[128+l_m] * B[5];
    C[32+l_m] += A[96+l_m] * B[6];
    C[32+l_m] += A[112+l_m] * B[7];
    C[32+l_m] += A[128+l_m] * B[8];
    C[48+l_m] += A[96+l_m] * B[9];
    C[48+l_m] += A[112+l_m] * B[10];
    C[64+l_m] += A[112+l_m] * B[11];
    C[64+l_m] += A[128+l_m] * B[12];
    C[80+l_m] += A[96+l_m] * B[13];
    C[80+l_m] += A[128+l_m] * B[14];
    C[96+l_m] += A[0+l_m] * B[15];
    C[96+l_m] += A[48+l_m] * B[16];
    C[96+l_m] += A[80+l_m] * B[17];
    C[112+l_m] += A[16+l_m] * B[18];
    C[112+l_m] += A[48+l_m] * B[19];
    C[112+l_m] += A[64+l_m] * B[20];
    C[128+l_m] += A[32+l_m] * B[21];
    C[128+l_m] += A[64+l_m] * B[22];
    C[128+l_m] += A[80+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 480;
#endif
}

void ssparse_fP113DivM_m10_n9_k10_ldAna3_ldB16_ldC16_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_n = 0;
  #pragma nounroll_and_jam
  for ( l_n = 0; l_n < 9; l_n++) {
    unsigned int l_m = 0;
   #pragma simd
    for ( l_m = 0; l_m < 10; l_m++) {
      C[(l_n*16)+l_m] = 0.0f;
    }
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b0 = _mm_broadcast_ss(&B[(l_n*16)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b0 = _mm_load_ss(&B[(l_n*16)+0]);    b0 = _mm_shuffle_ps(b0, b0, 0x00);
#endif
    __m128 c0_0 = _mm_load_ss(&C[(l_n*16)+0]);
    __m128 a0_0 = _mm_load_ss(&A[0]);
    c0_0 = _mm_add_ss(c0_0, _mm_mul_ss(a0_0, b0));
    _mm_store_ss(&C[(l_n*16)+0], c0_0);
    __m128 c0_1 = _mm_load_ss(&C[(l_n*16)+3]);
    __m128 a0_1 = _mm_load_ss(&A[1]);
    c0_1 = _mm_add_ss(c0_1, _mm_mul_ss(a0_1, b0));
    _mm_store_ss(&C[(l_n*16)+3], c0_1);
    __m128 c0_2 = _mm_load_ss(&C[(l_n*16)+9]);
    __m128 a0_2 = _mm_load_ss(&A[2]);
    c0_2 = _mm_add_ss(c0_2, _mm_mul_ss(a0_2, b0));
    _mm_store_ss(&C[(l_n*16)+9], c0_2);
#else
    C[(l_n*16)+0] += A[0] * B[(l_n*16)+0];
    C[(l_n*16)+3] += A[1] * B[(l_n*16)+0];
    C[(l_n*16)+9] += A[2] * B[(l_n*16)+0];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b1 = _mm_broadcast_ss(&B[(l_n*16)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b1 = _mm_load_ss(&B[(l_n*16)+1]);    b1 = _mm_shuffle_ps(b1, b1, 0x00);
#endif
    __m128 c1_0 = _mm_load_ss(&C[(l_n*16)+1]);
    __m128 a1_0 = _mm_load_ss(&A[3]);
    c1_0 = _mm_add_ss(c1_0, _mm_mul_ss(a1_0, b1));
    _mm_store_ss(&C[(l_n*16)+1], c1_0);
    __m128 c1_1 = _mm_load_ss(&C[(l_n*16)+7]);
    __m128 a1_1 = _mm_load_ss(&A[4]);
    c1_1 = _mm_add_ss(c1_1, _mm_mul_ss(a1_1, b1));
    _mm_store_ss(&C[(l_n*16)+7], c1_1);
#else
    C[(l_n*16)+1] += A[3] * B[(l_n*16)+1];
    C[(l_n*16)+7] += A[4] * B[(l_n*16)+1];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b2 = _mm_broadcast_ss(&B[(l_n*16)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b2 = _mm_load_ss(&B[(l_n*16)+2]);    b2 = _mm_shuffle_ps(b2, b2, 0x00);
#endif
    __m128 c2_0 = _mm_load_ss(&C[(l_n*16)+2]);
    __m128 a2_0 = _mm_load_ss(&A[5]);
    c2_0 = _mm_add_ss(c2_0, _mm_mul_ss(a2_0, b2));
    _mm_store_ss(&C[(l_n*16)+2], c2_0);
    __m128 c2_1 = _mm_load_ss(&C[(l_n*16)+8]);
    __m128 a2_1 = _mm_load_ss(&A[6]);
    c2_1 = _mm_add_ss(c2_1, _mm_mul_ss(a2_1, b2));
    _mm_store_ss(&C[(l_n*16)+8], c2_1);
#else
    C[(l_n*16)+2] += A[5] * B[(l_n*16)+2];
    C[(l_n*16)+8] += A[6] * B[(l_n*16)+2];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b3 = _mm_broadcast_ss(&B[(l_n*16)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b3 = _mm_load_ss(&B[(l_n*16)+3]);    b3 = _mm_shuffle_ps(b3, b3, 0x00);
#endif
    __m128 c3_0 = _mm_load_ss(&C[(l_n*16)+0]);
    __m128 a3_0 = _mm_load_ss(&A[7]);
    c3_0 = _mm_add_ss(c3_0, _mm_mul_ss(a3_0, b3));
    _mm_store_ss(&C[(l_n*16)+0], c3_0);
    __m128 c3_1 = _mm_load_ss(&C[(l_n*16)+3]);
    __m128 a3_1 = _mm_load_ss(&A[8]);
    c3_1 = _mm_add_ss(c3_1, _mm_mul_ss(a3_1, b3));
    _mm_store_ss(&C[(l_n*16)+3], c3_1);
    __m128 c3_2 = _mm_load_ss(&C[(l_n*16)+9]);
    __m128 a3_2 = _mm_load_ss(&A[9]);
    c3_2 = _mm_add_ss(c3_2, _mm_mul_ss(a3_2, b3));
    _mm_store_ss(&C[(l_n*16)+9], c3_2);
#else
    C[(l_n*16)+0] += A[7] * B[(l_n*16)+3];
    C[(l_n*16)+3] += A[8] * B[(l_n*16)+3];
    C[(l_n*16)+9] += A[9] * B[(l_n*16)+3];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b4 = _mm_broadcast_ss(&B[(l_n*16)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b4 = _mm_load_ss(&B[(l_n*16)+4]);    b4 = _mm_shuffle_ps(b4, b4, 0x00);
#endif
    __m128 c4_0 = _mm_load_ss(&C[(l_n*16)+4]);
    __m128 a4_0 = _mm_load_ss(&A[10]);
    c4_0 = _mm_add_ss(c4_0, _mm_mul_ss(a4_0, b4));
    _mm_store_ss(&C[(l_n*16)+4], c4_0);
#else
    C[(l_n*16)+4] += A[10] * B[(l_n*16)+4];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b5 = _mm_broadcast_ss(&B[(l_n*16)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b5 = _mm_load_ss(&B[(l_n*16)+5]);    b5 = _mm_shuffle_ps(b5, b5, 0x00);
#endif
    __m128 c5_0 = _mm_load_ss(&C[(l_n*16)+5]);
    __m128 a5_0 = _mm_load_ss(&A[11]);
    c5_0 = _mm_add_ss(c5_0, _mm_mul_ss(a5_0, b5));
    _mm_store_ss(&C[(l_n*16)+5], c5_0);
#else
    C[(l_n*16)+5] += A[11] * B[(l_n*16)+5];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b6 = _mm_broadcast_ss(&B[(l_n*16)+6]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b6 = _mm_load_ss(&B[(l_n*16)+6]);    b6 = _mm_shuffle_ps(b6, b6, 0x00);
#endif
    __m128 c6_0 = _mm_load_ss(&C[(l_n*16)+6]);
    __m128 a6_0 = _mm_load_ss(&A[12]);
    c6_0 = _mm_add_ss(c6_0, _mm_mul_ss(a6_0, b6));
    _mm_store_ss(&C[(l_n*16)+6], c6_0);
#else
    C[(l_n*16)+6] += A[12] * B[(l_n*16)+6];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b7 = _mm_broadcast_ss(&B[(l_n*16)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b7 = _mm_load_ss(&B[(l_n*16)+7]);    b7 = _mm_shuffle_ps(b7, b7, 0x00);
#endif
    __m128 c7_0 = _mm_load_ss(&C[(l_n*16)+1]);
    __m128 a7_0 = _mm_load_ss(&A[13]);
    c7_0 = _mm_add_ss(c7_0, _mm_mul_ss(a7_0, b7));
    _mm_store_ss(&C[(l_n*16)+1], c7_0);
    __m128 c7_1 = _mm_load_ss(&C[(l_n*16)+7]);
    __m128 a7_1 = _mm_load_ss(&A[14]);
    c7_1 = _mm_add_ss(c7_1, _mm_mul_ss(a7_1, b7));
    _mm_store_ss(&C[(l_n*16)+7], c7_1);
#else
    C[(l_n*16)+1] += A[13] * B[(l_n*16)+7];
    C[(l_n*16)+7] += A[14] * B[(l_n*16)+7];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b8 = _mm_broadcast_ss(&B[(l_n*16)+8]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b8 = _mm_load_ss(&B[(l_n*16)+8]);    b8 = _mm_shuffle_ps(b8, b8, 0x00);
#endif
    __m128 c8_0 = _mm_load_ss(&C[(l_n*16)+2]);
    __m128 a8_0 = _mm_load_ss(&A[15]);
    c8_0 = _mm_add_ss(c8_0, _mm_mul_ss(a8_0, b8));
    _mm_store_ss(&C[(l_n*16)+2], c8_0);
    __m128 c8_1 = _mm_load_ss(&C[(l_n*16)+8]);
    __m128 a8_1 = _mm_load_ss(&A[16]);
    c8_1 = _mm_add_ss(c8_1, _mm_mul_ss(a8_1, b8));
    _mm_store_ss(&C[(l_n*16)+8], c8_1);
#else
    C[(l_n*16)+2] += A[15] * B[(l_n*16)+8];
    C[(l_n*16)+8] += A[16] * B[(l_n*16)+8];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b9 = _mm_broadcast_ss(&B[(l_n*16)+9]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b9 = _mm_load_ss(&B[(l_n*16)+9]);    b9 = _mm_shuffle_ps(b9, b9, 0x00);
#endif
    __m128 c9_0 = _mm_load_ss(&C[(l_n*16)+0]);
    __m128 a9_0 = _mm_load_ss(&A[17]);
    c9_0 = _mm_add_ss(c9_0, _mm_mul_ss(a9_0, b9));
    _mm_store_ss(&C[(l_n*16)+0], c9_0);
    __m128 c9_1 = _mm_load_ss(&C[(l_n*16)+3]);
    __m128 a9_1 = _mm_load_ss(&A[18]);
    c9_1 = _mm_add_ss(c9_1, _mm_mul_ss(a9_1, b9));
    _mm_store_ss(&C[(l_n*16)+3], c9_1);
    __m128 c9_2 = _mm_load_ss(&C[(l_n*16)+9]);
    __m128 a9_2 = _mm_load_ss(&A[19]);
    c9_2 = _mm_add_ss(c9_2, _mm_mul_ss(a9_2, b9));
    _mm_store_ss(&C[(l_n*16)+9], c9_2);
#else
    C[(l_n*16)+0] += A[17] * B[(l_n*16)+9];
    C[(l_n*16)+3] += A[18] * B[(l_n*16)+9];
    C[(l_n*16)+9] += A[19] * B[(l_n*16)+9];
#endif

  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 360;
#endif
}

void ssparse_fP111DivM_m10_n9_k10_ldAna3_ldB16_ldC16_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_n = 0;
  #pragma nounroll_and_jam
  for ( l_n = 0; l_n < 9; l_n++) {
    unsigned int l_m = 0;
   #pragma simd
    for ( l_m = 0; l_m < 10; l_m++) {
      C[(l_n*16)+l_m] = 0.0f;
    }
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b0 = _mm_broadcast_ss(&B[(l_n*16)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b0 = _mm_load_ss(&B[(l_n*16)+0]);    b0 = _mm_shuffle_ps(b0, b0, 0x00);
#endif
    __m128 c0_0 = _mm_load_ss(&C[(l_n*16)+0]);
    __m128 a0_0 = _mm_load_ss(&A[0]);
    c0_0 = _mm_add_ss(c0_0, _mm_mul_ss(a0_0, b0));
    _mm_store_ss(&C[(l_n*16)+0], c0_0);
    __m128 c0_1 = _mm_load_ss(&C[(l_n*16)+3]);
    __m128 a0_1 = _mm_load_ss(&A[1]);
    c0_1 = _mm_add_ss(c0_1, _mm_mul_ss(a0_1, b0));
    _mm_store_ss(&C[(l_n*16)+3], c0_1);
    __m128 c0_2 = _mm_load_ss(&C[(l_n*16)+9]);
    __m128 a0_2 = _mm_load_ss(&A[2]);
    c0_2 = _mm_add_ss(c0_2, _mm_mul_ss(a0_2, b0));
    _mm_store_ss(&C[(l_n*16)+9], c0_2);
#else
    C[(l_n*16)+0] += A[0] * B[(l_n*16)+0];
    C[(l_n*16)+3] += A[1] * B[(l_n*16)+0];
    C[(l_n*16)+9] += A[2] * B[(l_n*16)+0];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b1 = _mm_broadcast_ss(&B[(l_n*16)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b1 = _mm_load_ss(&B[(l_n*16)+1]);    b1 = _mm_shuffle_ps(b1, b1, 0x00);
#endif
    __m128 c1_0 = _mm_castpd_ps(_mm_load_sd((const double*)&C[(l_n*16)+1]));
    __m128 a1_0 = _mm_castpd_ps(_mm_load_sd((const double*)&A[3]));
    c1_0 = _mm_add_ps(c1_0, _mm_mul_ps(a1_0, b1));
    _mm_store_sd((double*)&C[(l_n*16)+1], _mm_castps_pd(c1_0));
    __m128 c1_2 = _mm_castpd_ps(_mm_load_sd((const double*)&C[(l_n*16)+7]));
    __m128 a1_2 = _mm_castpd_ps(_mm_load_sd((const double*)&A[5]));
    c1_2 = _mm_add_ps(c1_2, _mm_mul_ps(a1_2, b1));
    _mm_store_sd((double*)&C[(l_n*16)+7], _mm_castps_pd(c1_2));
#else
    C[(l_n*16)+1] += A[3] * B[(l_n*16)+1];
    C[(l_n*16)+2] += A[4] * B[(l_n*16)+1];
    C[(l_n*16)+7] += A[5] * B[(l_n*16)+1];
    C[(l_n*16)+8] += A[6] * B[(l_n*16)+1];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b2 = _mm_broadcast_ss(&B[(l_n*16)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b2 = _mm_load_ss(&B[(l_n*16)+2]);    b2 = _mm_shuffle_ps(b2, b2, 0x00);
#endif
    __m128 c2_0 = _mm_castpd_ps(_mm_load_sd((const double*)&C[(l_n*16)+1]));
    __m128 a2_0 = _mm_castpd_ps(_mm_load_sd((const double*)&A[7]));
    c2_0 = _mm_add_ps(c2_0, _mm_mul_ps(a2_0, b2));
    _mm_store_sd((double*)&C[(l_n*16)+1], _mm_castps_pd(c2_0));
    __m128 c2_2 = _mm_castpd_ps(_mm_load_sd((const double*)&C[(l_n*16)+7]));
    __m128 a2_2 = _mm_castpd_ps(_mm_load_sd((const double*)&A[9]));
    c2_2 = _mm_add_ps(c2_2, _mm_mul_ps(a2_2, b2));
    _mm_store_sd((double*)&C[(l_n*16)+7], _mm_castps_pd(c2_2));
#else
    C[(l_n*16)+1] += A[7] * B[(l_n*16)+2];
    C[(l_n*16)+2] += A[8] * B[(l_n*16)+2];
    C[(l_n*16)+7] += A[9] * B[(l_n*16)+2];
    C[(l_n*16)+8] += A[10] * B[(l_n*16)+2];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b3 = _mm_broadcast_ss(&B[(l_n*16)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b3 = _mm_load_ss(&B[(l_n*16)+3]);    b3 = _mm_shuffle_ps(b3, b3, 0x00);
#endif
    __m128 c3_0 = _mm_load_ss(&C[(l_n*16)+0]);
    __m128 a3_0 = _mm_load_ss(&A[11]);
    c3_0 = _mm_add_ss(c3_0, _mm_mul_ss(a3_0, b3));
    _mm_store_ss(&C[(l_n*16)+0], c3_0);
    __m128 c3_1 = _mm_load_ss(&C[(l_n*16)+3]);
    __m128 a3_1 = _mm_load_ss(&A[12]);
    c3_1 = _mm_add_ss(c3_1, _mm_mul_ss(a3_1, b3));
    _mm_store_ss(&C[(l_n*16)+3], c3_1);
    __m128 c3_2 = _mm_load_ss(&C[(l_n*16)+9]);
    __m128 a3_2 = _mm_load_ss(&A[13]);
    c3_2 = _mm_add_ss(c3_2, _mm_mul_ss(a3_2, b3));
    _mm_store_ss(&C[(l_n*16)+9], c3_2);
#else
    C[(l_n*16)+0] += A[11] * B[(l_n*16)+3];
    C[(l_n*16)+3] += A[12] * B[(l_n*16)+3];
    C[(l_n*16)+9] += A[13] * B[(l_n*16)+3];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b4 = _mm_broadcast_ss(&B[(l_n*16)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b4 = _mm_load_ss(&B[(l_n*16)+4]);    b4 = _mm_shuffle_ps(b4, b4, 0x00);
#endif
    __m128 c4_0 = _mm_castpd_ps(_mm_load_sd((const double*)&C[(l_n*16)+4]));
    __m128 a4_0 = _mm_castpd_ps(_mm_load_sd((const double*)&A[14]));
    c4_0 = _mm_add_ps(c4_0, _mm_mul_ps(a4_0, b4));
    _mm_store_sd((double*)&C[(l_n*16)+4], _mm_castps_pd(c4_0));
    __m128 c4_2 = _mm_load_ss(&C[(l_n*16)+6]);
    __m128 a4_2 = _mm_load_ss(&A[16]);
    c4_2 = _mm_add_ss(c4_2, _mm_mul_ss(a4_2, b4));
    _mm_store_ss(&C[(l_n*16)+6], c4_2);
#else
    C[(l_n*16)+4] += A[14] * B[(l_n*16)+4];
    C[(l_n*16)+5] += A[15] * B[(l_n*16)+4];
    C[(l_n*16)+6] += A[16] * B[(l_n*16)+4];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b5 = _mm_broadcast_ss(&B[(l_n*16)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b5 = _mm_load_ss(&B[(l_n*16)+5]);    b5 = _mm_shuffle_ps(b5, b5, 0x00);
#endif
    __m128 c5_0 = _mm_castpd_ps(_mm_load_sd((const double*)&C[(l_n*16)+4]));
    __m128 a5_0 = _mm_castpd_ps(_mm_load_sd((const double*)&A[17]));
    c5_0 = _mm_add_ps(c5_0, _mm_mul_ps(a5_0, b5));
    _mm_store_sd((double*)&C[(l_n*16)+4], _mm_castps_pd(c5_0));
    __m128 c5_2 = _mm_load_ss(&C[(l_n*16)+6]);
    __m128 a5_2 = _mm_load_ss(&A[19]);
    c5_2 = _mm_add_ss(c5_2, _mm_mul_ss(a5_2, b5));
    _mm_store_ss(&C[(l_n*16)+6], c5_2);
#else
    C[(l_n*16)+4] += A[17] * B[(l_n*16)+5];
    C[(l_n*16)+5] += A[18] * B[(l_n*16)+5];
    C[(l_n*16)+6] += A[19] * B[(l_n*16)+5];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b6 = _mm_broadcast_ss(&B[(l_n*16)+6]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b6 = _mm_load_ss(&B[(l_n*16)+6]);    b6 = _mm_shuffle_ps(b6, b6, 0x00);
#endif
    __m128 c6_0 = _mm_castpd_ps(_mm_load_sd((const double*)&C[(l_n*16)+4]));
    __m128 a6_0 = _mm_castpd_ps(_mm_load_sd((const double*)&A[20]));
    c6_0 = _mm_add_ps(c6_0, _mm_mul_ps(a6_0, b6));
    _mm_store_sd((double*)&C[(l_n*16)+4], _mm_castps_pd(c6_0));
    __m128 c6_2 = _mm_load_ss(&C[(l_n*16)+6]);
    __m128 a6_2 = _mm_load_ss(&A[22]);
    c6_2 = _mm_add_ss(c6_2, _mm_mul_ss(a6_2, b6));
    _mm_store_ss(&C[(l_n*16)+6], c6_2);
#else
    C[(l_n*16)+4] += A[20] * B[(l_n*16)+6];
    C[(l_n*16)+5] += A[21] * B[(l_n*16)+6];
    C[(l_n*16)+6] += A[22] * B[(l_n*16)+6];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b7 = _mm_broadcast_ss(&B[(l_n*16)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b7 = _mm_load_ss(&B[(l_n*16)+7]);    b7 = _mm_shuffle_ps(b7, b7, 0x00);
#endif
    __m128 c7_0 = _mm_castpd_ps(_mm_load_sd((const double*)&C[(l_n*16)+1]));
    __m128 a7_0 = _mm_castpd_ps(_mm_load_sd((const double*)&A[23]));
    c7_0 = _mm_add_ps(c7_0, _mm_mul_ps(a7_0, b7));
    _mm_store_sd((double*)&C[(l_n*16)+1], _mm_castps_pd(c7_0));
    __m128 c7_2 = _mm_castpd_ps(_mm_load_sd((const double*)&C[(l_n*16)+7]));
    __m128 a7_2 = _mm_castpd_ps(_mm_load_sd((const double*)&A[25]));
    c7_2 = _mm_add_ps(c7_2, _mm_mul_ps(a7_2, b7));
    _mm_store_sd((double*)&C[(l_n*16)+7], _mm_castps_pd(c7_2));
#else
    C[(l_n*16)+1] += A[23] * B[(l_n*16)+7];
    C[(l_n*16)+2] += A[24] * B[(l_n*16)+7];
    C[(l_n*16)+7] += A[25] * B[(l_n*16)+7];
    C[(l_n*16)+8] += A[26] * B[(l_n*16)+7];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b8 = _mm_broadcast_ss(&B[(l_n*16)+8]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b8 = _mm_load_ss(&B[(l_n*16)+8]);    b8 = _mm_shuffle_ps(b8, b8, 0x00);
#endif
    __m128 c8_0 = _mm_castpd_ps(_mm_load_sd((const double*)&C[(l_n*16)+1]));
    __m128 a8_0 = _mm_castpd_ps(_mm_load_sd((const double*)&A[27]));
    c8_0 = _mm_add_ps(c8_0, _mm_mul_ps(a8_0, b8));
    _mm_store_sd((double*)&C[(l_n*16)+1], _mm_castps_pd(c8_0));
    __m128 c8_2 = _mm_castpd_ps(_mm_load_sd((const double*)&C[(l_n*16)+7]));
    __m128 a8_2 = _mm_castpd_ps(_mm_load_sd((const double*)&A[29]));
    c8_2 = _mm_add_ps(c8_2, _mm_mul_ps(a8_2, b8));
    _mm_store_sd((double*)&C[(l_n*16)+7], _mm_castps_pd(c8_2));
#else
    C[(l_n*16)+1] += A[27] * B[(l_n*16)+8];
    C[(l_n*16)+2] += A[28] * B[(l_n*16)+8];
    C[(l_n*16)+7] += A[29] * B[(l_n*16)+8];
    C[(l_n*16)+8] += A[30] * B[(l_n*16)+8];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b9 = _mm_broadcast_ss(&B[(l_n*16)+9]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b9 = _mm_load_ss(&B[(l_n*16)+9]);    b9 = _mm_shuffle_ps(b9, b9, 0x00);
#endif
    __m128 c9_0 = _mm_load_ss(&C[(l_n*16)+0]);
    __m128 a9_0 = _mm_load_ss(&A[31]);
    c9_0 = _mm_add_ss(c9_0, _mm_mul_ss(a9_0, b9));
    _mm_store_ss(&C[(l_n*16)+0], c9_0);
    __m128 c9_1 = _mm_load_ss(&C[(l_n*16)+3]);
    __m128 a9_1 = _mm_load_ss(&A[32]);
    c9_1 = _mm_add_ss(c9_1, _mm_mul_ss(a9_1, b9));
    _mm_store_ss(&C[(l_n*16)+3], c9_1);
    __m128 c9_2 = _mm_load_ss(&C[(l_n*16)+9]);
    __m128 a9_2 = _mm_load_ss(&A[33]);
    c9_2 = _mm_add_ss(c9_2, _mm_mul_ss(a9_2, b9));
    _mm_store_ss(&C[(l_n*16)+9], c9_2);
#else
    C[(l_n*16)+0] += A[31] * B[(l_n*16)+9];
    C[(l_n*16)+3] += A[32] * B[(l_n*16)+9];
    C[(l_n*16)+9] += A[33] * B[(l_n*16)+9];
#endif

  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 612;
#endif
}

void ssparse_starMatrix_m20_n9_k9_ldA24_ldBna4_ldC24_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 20; l_m++) {
    C[0+l_m] += A[144+l_m] * B[0];
    C[0+l_m] += A[168+l_m] * B[1];
    C[0+l_m] += A[192+l_m] * B[2];
    C[24+l_m] += A[144+l_m] * B[3];
    C[24+l_m] += A[168+l_m] * B[4];
    C[24+l_m] += A[192+l_m] * B[5];
    C[48+l_m] += A[144+l_m] * B[6];
    C[48+l_m] += A[168+l_m] * B[7];
    C[48+l_m] += A[192+l_m] * B[8];
    C[72+l_m] += A[144+l_m] * B[9];
    C[72+l_m] += A[168+l_m] * B[10];
    C[96+l_m] += A[168+l_m] * B[11];
    C[96+l_m] += A[192+l_m] * B[12];
    C[120+l_m] += A[144+l_m] * B[13];
    C[120+l_m] += A[192+l_m] * B[14];
    C[144+l_m] += A[0+l_m] * B[15];
    C[144+l_m] += A[72+l_m] * B[16];
    C[144+l_m] += A[120+l_m] * B[17];
    C[168+l_m] += A[24+l_m] * B[18];
    C[168+l_m] += A[72+l_m] * B[19];
    C[168+l_m] += A[96+l_m] * B[20];
    C[192+l_m] += A[48+l_m] * B[21];
    C[192+l_m] += A[96+l_m] * B[22];
    C[192+l_m] += A[120+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 960;
#endif
}

void ssparse_starMatrix_m35_n9_k9_ldA40_ldBna5_ldC40_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 35; l_m++) {
    C[0+l_m] += A[240+l_m] * B[0];
    C[0+l_m] += A[280+l_m] * B[1];
    C[0+l_m] += A[320+l_m] * B[2];
    C[40+l_m] += A[240+l_m] * B[3];
    C[40+l_m] += A[280+l_m] * B[4];
    C[40+l_m] += A[320+l_m] * B[5];
    C[80+l_m] += A[240+l_m] * B[6];
    C[80+l_m] += A[280+l_m] * B[7];
    C[80+l_m] += A[320+l_m] * B[8];
    C[120+l_m] += A[240+l_m] * B[9];
    C[120+l_m] += A[280+l_m] * B[10];
    C[160+l_m] += A[280+l_m] * B[11];
    C[160+l_m] += A[320+l_m] * B[12];
    C[200+l_m] += A[240+l_m] * B[13];
    C[200+l_m] += A[320+l_m] * B[14];
    C[240+l_m] += A[0+l_m] * B[15];
    C[240+l_m] += A[120+l_m] * B[16];
    C[240+l_m] += A[200+l_m] * B[17];
    C[280+l_m] += A[40+l_m] * B[18];
    C[280+l_m] += A[120+l_m] * B[19];
    C[280+l_m] += A[160+l_m] * B[20];
    C[320+l_m] += A[80+l_m] * B[21];
    C[320+l_m] += A[160+l_m] * B[22];
    C[320+l_m] += A[200+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 1680;
#endif
}

void ssparse_starMatrix_m56_n9_k9_ldA56_ldBna6_ldC56_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 56; l_m++) {
    C[0+l_m] += A[336+l_m] * B[0];
    C[0+l_m] += A[392+l_m] * B[1];
    C[0+l_m] += A[448+l_m] * B[2];
    C[56+l_m] += A[336+l_m] * B[3];
    C[56+l_m] += A[392+l_m] * B[4];
    C[56+l_m] += A[448+l_m] * B[5];
    C[112+l_m] += A[336+l_m] * B[6];
    C[112+l_m] += A[392+l_m] * B[7];
    C[112+l_m] += A[448+l_m] * B[8];
    C[168+l_m] += A[336+l_m] * B[9];
    C[168+l_m] += A[392+l_m] * B[10];
    C[224+l_m] += A[392+l_m] * B[11];
    C[224+l_m] += A[448+l_m] * B[12];
    C[280+l_m] += A[336+l_m] * B[13];
    C[280+l_m] += A[448+l_m] * B[14];
    C[336+l_m] += A[0+l_m] * B[15];
    C[336+l_m] += A[168+l_m] * B[16];
    C[336+l_m] += A[280+l_m] * B[17];
    C[392+l_m] += A[56+l_m] * B[18];
    C[392+l_m] += A[168+l_m] * B[19];
    C[392+l_m] += A[224+l_m] * B[20];
    C[448+l_m] += A[112+l_m] * B[21];
    C[448+l_m] += A[224+l_m] * B[22];
    C[448+l_m] += A[280+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 2688;
#endif
}

void ssparse_starMatrix_m84_n9_k9_ldA88_ldBna7_ldC88_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 84; l_m++) {
    C[0+l_m] += A[528+l_m] * B[0];
    C[0+l_m] += A[616+l_m] * B[1];
    C[0+l_m] += A[704+l_m] * B[2];
    C[88+l_m] += A[528+l_m] * B[3];
    C[88+l_m] += A[616+l_m] * B[4];
    C[88+l_m] += A[704+l_m] * B[5];
    C[176+l_m] += A[528+l_m] * B[6];
    C[176+l_m] += A[616+l_m] * B[7];
    C[176+l_m] += A[704+l_m] * B[8];
    C[264+l_m] += A[528+l_m] * B[9];
    C[264+l_m] += A[616+l_m] * B[10];
    C[352+l_m] += A[616+l_m] * B[11];
    C[352+l_m] += A[704+l_m] * B[12];
    C[440+l_m] += A[528+l_m] * B[13];
    C[440+l_m] += A[704+l_m] * B[14];
    C[528+l_m] += A[0+l_m] * B[15];
    C[528+l_m] += A[264+l_m] * B[16];
    C[528+l_m] += A[440+l_m] * B[17];
    C[616+l_m] += A[88+l_m] * B[18];
    C[616+l_m] += A[264+l_m] * B[19];
    C[616+l_m] += A[352+l_m] * B[20];
    C[704+l_m] += A[176+l_m] * B[21];
    C[704+l_m] += A[352+l_m] * B[22];
    C[704+l_m] += A[440+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 4032;
#endif
}

void ssparse_fM1DivM_m84_n9_k84_ldAna7_ldB88_ldC88_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_n = 0;
  #pragma nounroll_and_jam
  for ( l_n = 0; l_n < 9; l_n++) {
    unsigned int l_m = 0;
   #pragma simd
    for ( l_m = 0; l_m < 84; l_m++) {
      C[(l_n*88)+l_m] = 0.0f;
    }
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b0 = _mm_broadcast_ss(&B[(l_n*88)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b0 = _mm_load_ss(&B[(l_n*88)+0]);    b0 = _mm_shuffle_ps(b0, b0, 0x00);
#endif
    __m128 c0_0 = _mm_load_ss(&C[(l_n*88)+0]);
    __m128 a0_0 = _mm_load_ss(&A[0]);
    c0_0 = _mm_add_ss(c0_0, _mm_mul_ss(a0_0, b0));
    _mm_store_ss(&C[(l_n*88)+0], c0_0);
    __m128 c0_1 = _mm_load_ss(&C[(l_n*88)+3]);
    __m128 a0_1 = _mm_load_ss(&A[1]);
    c0_1 = _mm_add_ss(c0_1, _mm_mul_ss(a0_1, b0));
    _mm_store_ss(&C[(l_n*88)+3], c0_1);
    __m128 c0_2 = _mm_load_ss(&C[(l_n*88)+9]);
    __m128 a0_2 = _mm_load_ss(&A[2]);
    c0_2 = _mm_add_ss(c0_2, _mm_mul_ss(a0_2, b0));
    _mm_store_ss(&C[(l_n*88)+9], c0_2);
    __m128 c0_3 = _mm_load_ss(&C[(l_n*88)+19]);
    __m128 a0_3 = _mm_load_ss(&A[3]);
    c0_3 = _mm_add_ss(c0_3, _mm_mul_ss(a0_3, b0));
    _mm_store_ss(&C[(l_n*88)+19], c0_3);
    __m128 c0_4 = _mm_load_ss(&C[(l_n*88)+34]);
    __m128 a0_4 = _mm_load_ss(&A[4]);
    c0_4 = _mm_add_ss(c0_4, _mm_mul_ss(a0_4, b0));
    _mm_store_ss(&C[(l_n*88)+34], c0_4);
    __m128 c0_5 = _mm_load_ss(&C[(l_n*88)+55]);
    __m128 a0_5 = _mm_load_ss(&A[5]);
    c0_5 = _mm_add_ss(c0_5, _mm_mul_ss(a0_5, b0));
    _mm_store_ss(&C[(l_n*88)+55], c0_5);
    __m128 c0_6 = _mm_load_ss(&C[(l_n*88)+83]);
    __m128 a0_6 = _mm_load_ss(&A[6]);
    c0_6 = _mm_add_ss(c0_6, _mm_mul_ss(a0_6, b0));
    _mm_store_ss(&C[(l_n*88)+83], c0_6);
#else
    C[(l_n*88)+0] += A[0] * B[(l_n*88)+0];
    C[(l_n*88)+3] += A[1] * B[(l_n*88)+0];
    C[(l_n*88)+9] += A[2] * B[(l_n*88)+0];
    C[(l_n*88)+19] += A[3] * B[(l_n*88)+0];
    C[(l_n*88)+34] += A[4] * B[(l_n*88)+0];
    C[(l_n*88)+55] += A[5] * B[(l_n*88)+0];
    C[(l_n*88)+83] += A[6] * B[(l_n*88)+0];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b1 = _mm_broadcast_ss(&B[(l_n*88)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b1 = _mm_load_ss(&B[(l_n*88)+1]);    b1 = _mm_shuffle_ps(b1, b1, 0x00);
#endif
    __m128 c1_0 = _mm_load_ss(&C[(l_n*88)+1]);
    __m128 a1_0 = _mm_load_ss(&A[7]);
    c1_0 = _mm_add_ss(c1_0, _mm_mul_ss(a1_0, b1));
    _mm_store_ss(&C[(l_n*88)+1], c1_0);
    __m128 c1_1 = _mm_load_ss(&C[(l_n*88)+7]);
    __m128 a1_1 = _mm_load_ss(&A[8]);
    c1_1 = _mm_add_ss(c1_1, _mm_mul_ss(a1_1, b1));
    _mm_store_ss(&C[(l_n*88)+7], c1_1);
    __m128 c1_2 = _mm_load_ss(&C[(l_n*88)+17]);
    __m128 a1_2 = _mm_load_ss(&A[9]);
    c1_2 = _mm_add_ss(c1_2, _mm_mul_ss(a1_2, b1));
    _mm_store_ss(&C[(l_n*88)+17], c1_2);
    __m128 c1_3 = _mm_load_ss(&C[(l_n*88)+32]);
    __m128 a1_3 = _mm_load_ss(&A[10]);
    c1_3 = _mm_add_ss(c1_3, _mm_mul_ss(a1_3, b1));
    _mm_store_ss(&C[(l_n*88)+32], c1_3);
    __m128 c1_4 = _mm_load_ss(&C[(l_n*88)+53]);
    __m128 a1_4 = _mm_load_ss(&A[11]);
    c1_4 = _mm_add_ss(c1_4, _mm_mul_ss(a1_4, b1));
    _mm_store_ss(&C[(l_n*88)+53], c1_4);
    __m128 c1_5 = _mm_load_ss(&C[(l_n*88)+81]);
    __m128 a1_5 = _mm_load_ss(&A[12]);
    c1_5 = _mm_add_ss(c1_5, _mm_mul_ss(a1_5, b1));
    _mm_store_ss(&C[(l_n*88)+81], c1_5);
#else
    C[(l_n*88)+1] += A[7] * B[(l_n*88)+1];
    C[(l_n*88)+7] += A[8] * B[(l_n*88)+1];
    C[(l_n*88)+17] += A[9] * B[(l_n*88)+1];
    C[(l_n*88)+32] += A[10] * B[(l_n*88)+1];
    C[(l_n*88)+53] += A[11] * B[(l_n*88)+1];
    C[(l_n*88)+81] += A[12] * B[(l_n*88)+1];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b2 = _mm_broadcast_ss(&B[(l_n*88)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b2 = _mm_load_ss(&B[(l_n*88)+2]);    b2 = _mm_shuffle_ps(b2, b2, 0x00);
#endif
    __m128 c2_0 = _mm_load_ss(&C[(l_n*88)+2]);
    __m128 a2_0 = _mm_load_ss(&A[13]);
    c2_0 = _mm_add_ss(c2_0, _mm_mul_ss(a2_0, b2));
    _mm_store_ss(&C[(l_n*88)+2], c2_0);
    __m128 c2_1 = _mm_load_ss(&C[(l_n*88)+8]);
    __m128 a2_1 = _mm_load_ss(&A[14]);
    c2_1 = _mm_add_ss(c2_1, _mm_mul_ss(a2_1, b2));
    _mm_store_ss(&C[(l_n*88)+8], c2_1);
    __m128 c2_2 = _mm_load_ss(&C[(l_n*88)+18]);
    __m128 a2_2 = _mm_load_ss(&A[15]);
    c2_2 = _mm_add_ss(c2_2, _mm_mul_ss(a2_2, b2));
    _mm_store_ss(&C[(l_n*88)+18], c2_2);
    __m128 c2_3 = _mm_load_ss(&C[(l_n*88)+33]);
    __m128 a2_3 = _mm_load_ss(&A[16]);
    c2_3 = _mm_add_ss(c2_3, _mm_mul_ss(a2_3, b2));
    _mm_store_ss(&C[(l_n*88)+33], c2_3);
    __m128 c2_4 = _mm_load_ss(&C[(l_n*88)+54]);
    __m128 a2_4 = _mm_load_ss(&A[17]);
    c2_4 = _mm_add_ss(c2_4, _mm_mul_ss(a2_4, b2));
    _mm_store_ss(&C[(l_n*88)+54], c2_4);
    __m128 c2_5 = _mm_load_ss(&C[(l_n*88)+82]);
    __m128 a2_5 = _mm_load_ss(&A[18]);
    c2_5 = _mm_add_ss(c2_5, _mm_mul_ss(a2_5, b2));
    _mm_store_ss(&C[(l_n*88)+82], c2_5);
#else
    C[(l_n*88)+2] += A[13] * B[(l_n*88)+2];
    C[(l_n*88)+8] += A[14] * B[(l_n*88)+2];
    C[(l_n*88)+18] += A[15] * B[(l_n*88)+2];
    C[(l_n*88)+33] += A[16] * B[(l_n*88)+2];
    C[(l_n*88)+54] += A[17] * B[(l_n*88)+2];
    C[(l_n*88)+82] += A[18] * B[(l_n*88)+2];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b3 = _mm_broadcast_ss(&B[(l_n*88)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b3 = _mm_load_ss(&B[(l_n*88)+3]);    b3 = _mm_shuffle_ps(b3, b3, 0x00);
#endif
    __m128 c3_0 = _mm_load_ss(&C[(l_n*88)+0]);
    __m128 a3_0 = _mm_load_ss(&A[19]);
    c3_0 = _mm_add_ss(c3_0, _mm_mul_ss(a3_0, b3));
    _mm_store_ss(&C[(l_n*88)+0], c3_0);
    __m128 c3_1 = _mm_load_ss(&C[(l_n*88)+3]);
    __m128 a3_1 = _mm_load_ss(&A[20]);
    c3_1 = _mm_add_ss(c3_1, _mm_mul_ss(a3_1, b3));
    _mm_store_ss(&C[(l_n*88)+3], c3_1);
    __m128 c3_2 = _mm_load_ss(&C[(l_n*88)+9]);
    __m128 a3_2 = _mm_load_ss(&A[21]);
    c3_2 = _mm_add_ss(c3_2, _mm_mul_ss(a3_2, b3));
    _mm_store_ss(&C[(l_n*88)+9], c3_2);
    __m128 c3_3 = _mm_load_ss(&C[(l_n*88)+19]);
    __m128 a3_3 = _mm_load_ss(&A[22]);
    c3_3 = _mm_add_ss(c3_3, _mm_mul_ss(a3_3, b3));
    _mm_store_ss(&C[(l_n*88)+19], c3_3);
    __m128 c3_4 = _mm_load_ss(&C[(l_n*88)+34]);
    __m128 a3_4 = _mm_load_ss(&A[23]);
    c3_4 = _mm_add_ss(c3_4, _mm_mul_ss(a3_4, b3));
    _mm_store_ss(&C[(l_n*88)+34], c3_4);
    __m128 c3_5 = _mm_load_ss(&C[(l_n*88)+55]);
    __m128 a3_5 = _mm_load_ss(&A[24]);
    c3_5 = _mm_add_ss(c3_5, _mm_mul_ss(a3_5, b3));
    _mm_store_ss(&C[(l_n*88)+55], c3_5);
    __m128 c3_6 = _mm_load_ss(&C[(l_n*88)+83]);
    __m128 a3_6 = _mm_load_ss(&A[25]);
    c3_6 = _mm_add_ss(c3_6, _mm_mul_ss(a3_6, b3));
    _mm_store_ss(&C[(l_n*88)+83], c3_6);
#else
    C[(l_n*88)+0] += A[19] * B[(l_n*88)+3];
    C[(l_n*88)+3] += A[20] * B[(l_n*88)+3];
    C[(l_n*88)+9] += A[21] * B[(l_n*88)+3];
    C[(l_n*88)+19] += A[22] * B[(l_n*88)+3];
    C[(l_n*88)+34] += A[23] * B[(l_n*88)+3];
    C[(l_n*88)+55] += A[24] * B[(l_n*88)+3];
    C[(l_n*88)+83] += A[25] * B[(l_n*88)+3];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b4 = _mm_broadcast_ss(&B[(l_n*88)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b4 = _mm_load_ss(&B[(l_n*88)+4]);    b4 = _mm_shuffle_ps(b4, b4, 0x00);
#endif
    __m128 c4_0 = _mm_load_ss(&C[(l_n*88)+4]);
    __m128 a4_0 = _mm_load_ss(&A[26]);
    c4_0 = _mm_add_ss(c4_0, _mm_mul_ss(a4_0, b4));
    _mm_store_ss(&C[(l_n*88)+4], c4_0);
    __m128 c4_1 = _mm_load_ss(&C[(l_n*88)+14]);
    __m128 a4_1 = _mm_load_ss(&A[27]);
    c4_1 = _mm_add_ss(c4_1, _mm_mul_ss(a4_1, b4));
    _mm_store_ss(&C[(l_n*88)+14], c4_1);
    __m128 c4_2 = _mm_load_ss(&C[(l_n*88)+29]);
    __m128 a4_2 = _mm_load_ss(&A[28]);
    c4_2 = _mm_add_ss(c4_2, _mm_mul_ss(a4_2, b4));
    _mm_store_ss(&C[(l_n*88)+29], c4_2);
    __m128 c4_3 = _mm_load_ss(&C[(l_n*88)+50]);
    __m128 a4_3 = _mm_load_ss(&A[29]);
    c4_3 = _mm_add_ss(c4_3, _mm_mul_ss(a4_3, b4));
    _mm_store_ss(&C[(l_n*88)+50], c4_3);
    __m128 c4_4 = _mm_load_ss(&C[(l_n*88)+78]);
    __m128 a4_4 = _mm_load_ss(&A[30]);
    c4_4 = _mm_add_ss(c4_4, _mm_mul_ss(a4_4, b4));
    _mm_store_ss(&C[(l_n*88)+78], c4_4);
#else
    C[(l_n*88)+4] += A[26] * B[(l_n*88)+4];
    C[(l_n*88)+14] += A[27] * B[(l_n*88)+4];
    C[(l_n*88)+29] += A[28] * B[(l_n*88)+4];
    C[(l_n*88)+50] += A[29] * B[(l_n*88)+4];
    C[(l_n*88)+78] += A[30] * B[(l_n*88)+4];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b5 = _mm_broadcast_ss(&B[(l_n*88)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b5 = _mm_load_ss(&B[(l_n*88)+5]);    b5 = _mm_shuffle_ps(b5, b5, 0x00);
#endif
    __m128 c5_0 = _mm_load_ss(&C[(l_n*88)+5]);
    __m128 a5_0 = _mm_load_ss(&A[31]);
    c5_0 = _mm_add_ss(c5_0, _mm_mul_ss(a5_0, b5));
    _mm_store_ss(&C[(l_n*88)+5], c5_0);
    __m128 c5_1 = _mm_load_ss(&C[(l_n*88)+15]);
    __m128 a5_1 = _mm_load_ss(&A[32]);
    c5_1 = _mm_add_ss(c5_1, _mm_mul_ss(a5_1, b5));
    _mm_store_ss(&C[(l_n*88)+15], c5_1);
    __m128 c5_2 = _mm_load_ss(&C[(l_n*88)+30]);
    __m128 a5_2 = _mm_load_ss(&A[33]);
    c5_2 = _mm_add_ss(c5_2, _mm_mul_ss(a5_2, b5));
    _mm_store_ss(&C[(l_n*88)+30], c5_2);
    __m128 c5_3 = _mm_load_ss(&C[(l_n*88)+51]);
    __m128 a5_3 = _mm_load_ss(&A[34]);
    c5_3 = _mm_add_ss(c5_3, _mm_mul_ss(a5_3, b5));
    _mm_store_ss(&C[(l_n*88)+51], c5_3);
    __m128 c5_4 = _mm_load_ss(&C[(l_n*88)+79]);
    __m128 a5_4 = _mm_load_ss(&A[35]);
    c5_4 = _mm_add_ss(c5_4, _mm_mul_ss(a5_4, b5));
    _mm_store_ss(&C[(l_n*88)+79], c5_4);
#else
    C[(l_n*88)+5] += A[31] * B[(l_n*88)+5];
    C[(l_n*88)+15] += A[32] * B[(l_n*88)+5];
    C[(l_n*88)+30] += A[33] * B[(l_n*88)+5];
    C[(l_n*88)+51] += A[34] * B[(l_n*88)+5];
    C[(l_n*88)+79] += A[35] * B[(l_n*88)+5];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b6 = _mm_broadcast_ss(&B[(l_n*88)+6]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b6 = _mm_load_ss(&B[(l_n*88)+6]);    b6 = _mm_shuffle_ps(b6, b6, 0x00);
#endif
    __m128 c6_0 = _mm_load_ss(&C[(l_n*88)+6]);
    __m128 a6_0 = _mm_load_ss(&A[36]);
    c6_0 = _mm_add_ss(c6_0, _mm_mul_ss(a6_0, b6));
    _mm_store_ss(&C[(l_n*88)+6], c6_0);
    __m128 c6_1 = _mm_load_ss(&C[(l_n*88)+16]);
    __m128 a6_1 = _mm_load_ss(&A[37]);
    c6_1 = _mm_add_ss(c6_1, _mm_mul_ss(a6_1, b6));
    _mm_store_ss(&C[(l_n*88)+16], c6_1);
    __m128 c6_2 = _mm_load_ss(&C[(l_n*88)+31]);
    __m128 a6_2 = _mm_load_ss(&A[38]);
    c6_2 = _mm_add_ss(c6_2, _mm_mul_ss(a6_2, b6));
    _mm_store_ss(&C[(l_n*88)+31], c6_2);
    __m128 c6_3 = _mm_load_ss(&C[(l_n*88)+52]);
    __m128 a6_3 = _mm_load_ss(&A[39]);
    c6_3 = _mm_add_ss(c6_3, _mm_mul_ss(a6_3, b6));
    _mm_store_ss(&C[(l_n*88)+52], c6_3);
    __m128 c6_4 = _mm_load_ss(&C[(l_n*88)+80]);
    __m128 a6_4 = _mm_load_ss(&A[40]);
    c6_4 = _mm_add_ss(c6_4, _mm_mul_ss(a6_4, b6));
    _mm_store_ss(&C[(l_n*88)+80], c6_4);
#else
    C[(l_n*88)+6] += A[36] * B[(l_n*88)+6];
    C[(l_n*88)+16] += A[37] * B[(l_n*88)+6];
    C[(l_n*88)+31] += A[38] * B[(l_n*88)+6];
    C[(l_n*88)+52] += A[39] * B[(l_n*88)+6];
    C[(l_n*88)+80] += A[40] * B[(l_n*88)+6];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b7 = _mm_broadcast_ss(&B[(l_n*88)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b7 = _mm_load_ss(&B[(l_n*88)+7]);    b7 = _mm_shuffle_ps(b7, b7, 0x00);
#endif
    __m128 c7_0 = _mm_load_ss(&C[(l_n*88)+1]);
    __m128 a7_0 = _mm_load_ss(&A[41]);
    c7_0 = _mm_add_ss(c7_0, _mm_mul_ss(a7_0, b7));
    _mm_store_ss(&C[(l_n*88)+1], c7_0);
    __m128 c7_1 = _mm_load_ss(&C[(l_n*88)+7]);
    __m128 a7_1 = _mm_load_ss(&A[42]);
    c7_1 = _mm_add_ss(c7_1, _mm_mul_ss(a7_1, b7));
    _mm_store_ss(&C[(l_n*88)+7], c7_1);
    __m128 c7_2 = _mm_load_ss(&C[(l_n*88)+17]);
    __m128 a7_2 = _mm_load_ss(&A[43]);
    c7_2 = _mm_add_ss(c7_2, _mm_mul_ss(a7_2, b7));
    _mm_store_ss(&C[(l_n*88)+17], c7_2);
    __m128 c7_3 = _mm_load_ss(&C[(l_n*88)+32]);
    __m128 a7_3 = _mm_load_ss(&A[44]);
    c7_3 = _mm_add_ss(c7_3, _mm_mul_ss(a7_3, b7));
    _mm_store_ss(&C[(l_n*88)+32], c7_3);
    __m128 c7_4 = _mm_load_ss(&C[(l_n*88)+53]);
    __m128 a7_4 = _mm_load_ss(&A[45]);
    c7_4 = _mm_add_ss(c7_4, _mm_mul_ss(a7_4, b7));
    _mm_store_ss(&C[(l_n*88)+53], c7_4);
    __m128 c7_5 = _mm_load_ss(&C[(l_n*88)+81]);
    __m128 a7_5 = _mm_load_ss(&A[46]);
    c7_5 = _mm_add_ss(c7_5, _mm_mul_ss(a7_5, b7));
    _mm_store_ss(&C[(l_n*88)+81], c7_5);
#else
    C[(l_n*88)+1] += A[41] * B[(l_n*88)+7];
    C[(l_n*88)+7] += A[42] * B[(l_n*88)+7];
    C[(l_n*88)+17] += A[43] * B[(l_n*88)+7];
    C[(l_n*88)+32] += A[44] * B[(l_n*88)+7];
    C[(l_n*88)+53] += A[45] * B[(l_n*88)+7];
    C[(l_n*88)+81] += A[46] * B[(l_n*88)+7];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b8 = _mm_broadcast_ss(&B[(l_n*88)+8]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b8 = _mm_load_ss(&B[(l_n*88)+8]);    b8 = _mm_shuffle_ps(b8, b8, 0x00);
#endif
    __m128 c8_0 = _mm_load_ss(&C[(l_n*88)+2]);
    __m128 a8_0 = _mm_load_ss(&A[47]);
    c8_0 = _mm_add_ss(c8_0, _mm_mul_ss(a8_0, b8));
    _mm_store_ss(&C[(l_n*88)+2], c8_0);
    __m128 c8_1 = _mm_load_ss(&C[(l_n*88)+8]);
    __m128 a8_1 = _mm_load_ss(&A[48]);
    c8_1 = _mm_add_ss(c8_1, _mm_mul_ss(a8_1, b8));
    _mm_store_ss(&C[(l_n*88)+8], c8_1);
    __m128 c8_2 = _mm_load_ss(&C[(l_n*88)+18]);
    __m128 a8_2 = _mm_load_ss(&A[49]);
    c8_2 = _mm_add_ss(c8_2, _mm_mul_ss(a8_2, b8));
    _mm_store_ss(&C[(l_n*88)+18], c8_2);
    __m128 c8_3 = _mm_load_ss(&C[(l_n*88)+33]);
    __m128 a8_3 = _mm_load_ss(&A[50]);
    c8_3 = _mm_add_ss(c8_3, _mm_mul_ss(a8_3, b8));
    _mm_store_ss(&C[(l_n*88)+33], c8_3);
    __m128 c8_4 = _mm_load_ss(&C[(l_n*88)+54]);
    __m128 a8_4 = _mm_load_ss(&A[51]);
    c8_4 = _mm_add_ss(c8_4, _mm_mul_ss(a8_4, b8));
    _mm_store_ss(&C[(l_n*88)+54], c8_4);
    __m128 c8_5 = _mm_load_ss(&C[(l_n*88)+82]);
    __m128 a8_5 = _mm_load_ss(&A[52]);
    c8_5 = _mm_add_ss(c8_5, _mm_mul_ss(a8_5, b8));
    _mm_store_ss(&C[(l_n*88)+82], c8_5);
#else
    C[(l_n*88)+2] += A[47] * B[(l_n*88)+8];
    C[(l_n*88)+8] += A[48] * B[(l_n*88)+8];
    C[(l_n*88)+18] += A[49] * B[(l_n*88)+8];
    C[(l_n*88)+33] += A[50] * B[(l_n*88)+8];
    C[(l_n*88)+54] += A[51] * B[(l_n*88)+8];
    C[(l_n*88)+82] += A[52] * B[(l_n*88)+8];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b9 = _mm_broadcast_ss(&B[(l_n*88)+9]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b9 = _mm_load_ss(&B[(l_n*88)+9]);    b9 = _mm_shuffle_ps(b9, b9, 0x00);
#endif
    __m128 c9_0 = _mm_load_ss(&C[(l_n*88)+0]);
    __m128 a9_0 = _mm_load_ss(&A[53]);
    c9_0 = _mm_add_ss(c9_0, _mm_mul_ss(a9_0, b9));
    _mm_store_ss(&C[(l_n*88)+0], c9_0);
    __m128 c9_1 = _mm_load_ss(&C[(l_n*88)+3]);
    __m128 a9_1 = _mm_load_ss(&A[54]);
    c9_1 = _mm_add_ss(c9_1, _mm_mul_ss(a9_1, b9));
    _mm_store_ss(&C[(l_n*88)+3], c9_1);
    __m128 c9_2 = _mm_load_ss(&C[(l_n*88)+9]);
    __m128 a9_2 = _mm_load_ss(&A[55]);
    c9_2 = _mm_add_ss(c9_2, _mm_mul_ss(a9_2, b9));
    _mm_store_ss(&C[(l_n*88)+9], c9_2);
    __m128 c9_3 = _mm_load_ss(&C[(l_n*88)+19]);
    __m128 a9_3 = _mm_load_ss(&A[56]);
    c9_3 = _mm_add_ss(c9_3, _mm_mul_ss(a9_3, b9));
    _mm_store_ss(&C[(l_n*88)+19], c9_3);
    __m128 c9_4 = _mm_load_ss(&C[(l_n*88)+34]);
    __m128 a9_4 = _mm_load_ss(&A[57]);
    c9_4 = _mm_add_ss(c9_4, _mm_mul_ss(a9_4, b9));
    _mm_store_ss(&C[(l_n*88)+34], c9_4);
    __m128 c9_5 = _mm_load_ss(&C[(l_n*88)+55]);
    __m128 a9_5 = _mm_load_ss(&A[58]);
    c9_5 = _mm_add_ss(c9_5, _mm_mul_ss(a9_5, b9));
    _mm_store_ss(&C[(l_n*88)+55], c9_5);
    __m128 c9_6 = _mm_load_ss(&C[(l_n*88)+83]);
    __m128 a9_6 = _mm_load_ss(&A[59]);
    c9_6 = _mm_add_ss(c9_6, _mm_mul_ss(a9_6, b9));
    _mm_store_ss(&C[(l_n*88)+83], c9_6);
#else
    C[(l_n*88)+0] += A[53] * B[(l_n*88)+9];
    C[(l_n*88)+3] += A[54] * B[(l_n*88)+9];
    C[(l_n*88)+9] += A[55] * B[(l_n*88)+9];
    C[(l_n*88)+19] += A[56] * B[(l_n*88)+9];
    C[(l_n*88)+34] += A[57] * B[(l_n*88)+9];
    C[(l_n*88)+55] += A[58] * B[(l_n*88)+9];
    C[(l_n*88)+83] += A[59] * B[(l_n*88)+9];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b10 = _mm_broadcast_ss(&B[(l_n*88)+10]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b10 = _mm_load_ss(&B[(l_n*88)+10]);    b10 = _mm_shuffle_ps(b10, b10, 0x00);
#endif
    __m128 c10_0 = _mm_load_ss(&C[(l_n*88)+10]);
    __m128 a10_0 = _mm_load_ss(&A[60]);
    c10_0 = _mm_add_ss(c10_0, _mm_mul_ss(a10_0, b10));
    _mm_store_ss(&C[(l_n*88)+10], c10_0);
    __m128 c10_1 = _mm_load_ss(&C[(l_n*88)+25]);
    __m128 a10_1 = _mm_load_ss(&A[61]);
    c10_1 = _mm_add_ss(c10_1, _mm_mul_ss(a10_1, b10));
    _mm_store_ss(&C[(l_n*88)+25], c10_1);
    __m128 c10_2 = _mm_load_ss(&C[(l_n*88)+46]);
    __m128 a10_2 = _mm_load_ss(&A[62]);
    c10_2 = _mm_add_ss(c10_2, _mm_mul_ss(a10_2, b10));
    _mm_store_ss(&C[(l_n*88)+46], c10_2);
    __m128 c10_3 = _mm_load_ss(&C[(l_n*88)+74]);
    __m128 a10_3 = _mm_load_ss(&A[63]);
    c10_3 = _mm_add_ss(c10_3, _mm_mul_ss(a10_3, b10));
    _mm_store_ss(&C[(l_n*88)+74], c10_3);
#else
    C[(l_n*88)+10] += A[60] * B[(l_n*88)+10];
    C[(l_n*88)+25] += A[61] * B[(l_n*88)+10];
    C[(l_n*88)+46] += A[62] * B[(l_n*88)+10];
    C[(l_n*88)+74] += A[63] * B[(l_n*88)+10];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b11 = _mm_broadcast_ss(&B[(l_n*88)+11]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b11 = _mm_load_ss(&B[(l_n*88)+11]);    b11 = _mm_shuffle_ps(b11, b11, 0x00);
#endif
    __m128 c11_0 = _mm_load_ss(&C[(l_n*88)+11]);
    __m128 a11_0 = _mm_load_ss(&A[64]);
    c11_0 = _mm_add_ss(c11_0, _mm_mul_ss(a11_0, b11));
    _mm_store_ss(&C[(l_n*88)+11], c11_0);
    __m128 c11_1 = _mm_load_ss(&C[(l_n*88)+26]);
    __m128 a11_1 = _mm_load_ss(&A[65]);
    c11_1 = _mm_add_ss(c11_1, _mm_mul_ss(a11_1, b11));
    _mm_store_ss(&C[(l_n*88)+26], c11_1);
    __m128 c11_2 = _mm_load_ss(&C[(l_n*88)+47]);
    __m128 a11_2 = _mm_load_ss(&A[66]);
    c11_2 = _mm_add_ss(c11_2, _mm_mul_ss(a11_2, b11));
    _mm_store_ss(&C[(l_n*88)+47], c11_2);
    __m128 c11_3 = _mm_load_ss(&C[(l_n*88)+75]);
    __m128 a11_3 = _mm_load_ss(&A[67]);
    c11_3 = _mm_add_ss(c11_3, _mm_mul_ss(a11_3, b11));
    _mm_store_ss(&C[(l_n*88)+75], c11_3);
#else
    C[(l_n*88)+11] += A[64] * B[(l_n*88)+11];
    C[(l_n*88)+26] += A[65] * B[(l_n*88)+11];
    C[(l_n*88)+47] += A[66] * B[(l_n*88)+11];
    C[(l_n*88)+75] += A[67] * B[(l_n*88)+11];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b12 = _mm_broadcast_ss(&B[(l_n*88)+12]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b12 = _mm_load_ss(&B[(l_n*88)+12]);    b12 = _mm_shuffle_ps(b12, b12, 0x00);
#endif
    __m128 c12_0 = _mm_load_ss(&C[(l_n*88)+12]);
    __m128 a12_0 = _mm_load_ss(&A[68]);
    c12_0 = _mm_add_ss(c12_0, _mm_mul_ss(a12_0, b12));
    _mm_store_ss(&C[(l_n*88)+12], c12_0);
    __m128 c12_1 = _mm_load_ss(&C[(l_n*88)+27]);
    __m128 a12_1 = _mm_load_ss(&A[69]);
    c12_1 = _mm_add_ss(c12_1, _mm_mul_ss(a12_1, b12));
    _mm_store_ss(&C[(l_n*88)+27], c12_1);
    __m128 c12_2 = _mm_load_ss(&C[(l_n*88)+48]);
    __m128 a12_2 = _mm_load_ss(&A[70]);
    c12_2 = _mm_add_ss(c12_2, _mm_mul_ss(a12_2, b12));
    _mm_store_ss(&C[(l_n*88)+48], c12_2);
    __m128 c12_3 = _mm_load_ss(&C[(l_n*88)+76]);
    __m128 a12_3 = _mm_load_ss(&A[71]);
    c12_3 = _mm_add_ss(c12_3, _mm_mul_ss(a12_3, b12));
    _mm_store_ss(&C[(l_n*88)+76], c12_3);
#else
    C[(l_n*88)+12] += A[68] * B[(l_n*88)+12];
    C[(l_n*88)+27] += A[69] * B[(l_n*88)+12];
    C[(l_n*88)+48] += A[70] * B[(l_n*88)+12];
    C[(l_n*88)+76] += A[71] * B[(l_n*88)+12];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b13 = _mm_broadcast_ss(&B[(l_n*88)+13]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b13 = _mm_load_ss(&B[(l_n*88)+13]);    b13 = _mm_shuffle_ps(b13, b13, 0x00);
#endif
    __m128 c13_0 = _mm_load_ss(&C[(l_n*88)+13]);
    __m128 a13_0 = _mm_load_ss(&A[72]);
    c13_0 = _mm_add_ss(c13_0, _mm_mul_ss(a13_0, b13));
    _mm_store_ss(&C[(l_n*88)+13], c13_0);
    __m128 c13_1 = _mm_load_ss(&C[(l_n*88)+28]);
    __m128 a13_1 = _mm_load_ss(&A[73]);
    c13_1 = _mm_add_ss(c13_1, _mm_mul_ss(a13_1, b13));
    _mm_store_ss(&C[(l_n*88)+28], c13_1);
    __m128 c13_2 = _mm_load_ss(&C[(l_n*88)+49]);
    __m128 a13_2 = _mm_load_ss(&A[74]);
    c13_2 = _mm_add_ss(c13_2, _mm_mul_ss(a13_2, b13));
    _mm_store_ss(&C[(l_n*88)+49], c13_2);
    __m128 c13_3 = _mm_load_ss(&C[(l_n*88)+77]);
    __m128 a13_3 = _mm_load_ss(&A[75]);
    c13_3 = _mm_add_ss(c13_3, _mm_mul_ss(a13_3, b13));
    _mm_store_ss(&C[(l_n*88)+77], c13_3);
#else
    C[(l_n*88)+13] += A[72] * B[(l_n*88)+13];
    C[(l_n*88)+28] += A[73] * B[(l_n*88)+13];
    C[(l_n*88)+49] += A[74] * B[(l_n*88)+13];
    C[(l_n*88)+77] += A[75] * B[(l_n*88)+13];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b14 = _mm_broadcast_ss(&B[(l_n*88)+14]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b14 = _mm_load_ss(&B[(l_n*88)+14]);    b14 = _mm_shuffle_ps(b14, b14, 0x00);
#endif
    __m128 c14_0 = _mm_load_ss(&C[(l_n*88)+4]);
    __m128 a14_0 = _mm_load_ss(&A[76]);
    c14_0 = _mm_add_ss(c14_0, _mm_mul_ss(a14_0, b14));
    _mm_store_ss(&C[(l_n*88)+4], c14_0);
    __m128 c14_1 = _mm_load_ss(&C[(l_n*88)+14]);
    __m128 a14_1 = _mm_load_ss(&A[77]);
    c14_1 = _mm_add_ss(c14_1, _mm_mul_ss(a14_1, b14));
    _mm_store_ss(&C[(l_n*88)+14], c14_1);
    __m128 c14_2 = _mm_load_ss(&C[(l_n*88)+29]);
    __m128 a14_2 = _mm_load_ss(&A[78]);
    c14_2 = _mm_add_ss(c14_2, _mm_mul_ss(a14_2, b14));
    _mm_store_ss(&C[(l_n*88)+29], c14_2);
    __m128 c14_3 = _mm_load_ss(&C[(l_n*88)+50]);
    __m128 a14_3 = _mm_load_ss(&A[79]);
    c14_3 = _mm_add_ss(c14_3, _mm_mul_ss(a14_3, b14));
    _mm_store_ss(&C[(l_n*88)+50], c14_3);
    __m128 c14_4 = _mm_load_ss(&C[(l_n*88)+78]);
    __m128 a14_4 = _mm_load_ss(&A[80]);
    c14_4 = _mm_add_ss(c14_4, _mm_mul_ss(a14_4, b14));
    _mm_store_ss(&C[(l_n*88)+78], c14_4);
#else
    C[(l_n*88)+4] += A[76] * B[(l_n*88)+14];
    C[(l_n*88)+14] += A[77] * B[(l_n*88)+14];
    C[(l_n*88)+29] += A[78] * B[(l_n*88)+14];
    C[(l_n*88)+50] += A[79] * B[(l_n*88)+14];
    C[(l_n*88)+78] += A[80] * B[(l_n*88)+14];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b15 = _mm_broadcast_ss(&B[(l_n*88)+15]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b15 = _mm_load_ss(&B[(l_n*88)+15]);    b15 = _mm_shuffle_ps(b15, b15, 0x00);
#endif
    __m128 c15_0 = _mm_load_ss(&C[(l_n*88)+5]);
    __m128 a15_0 = _mm_load_ss(&A[81]);
    c15_0 = _mm_add_ss(c15_0, _mm_mul_ss(a15_0, b15));
    _mm_store_ss(&C[(l_n*88)+5], c15_0);
    __m128 c15_1 = _mm_load_ss(&C[(l_n*88)+15]);
    __m128 a15_1 = _mm_load_ss(&A[82]);
    c15_1 = _mm_add_ss(c15_1, _mm_mul_ss(a15_1, b15));
    _mm_store_ss(&C[(l_n*88)+15], c15_1);
    __m128 c15_2 = _mm_load_ss(&C[(l_n*88)+30]);
    __m128 a15_2 = _mm_load_ss(&A[83]);
    c15_2 = _mm_add_ss(c15_2, _mm_mul_ss(a15_2, b15));
    _mm_store_ss(&C[(l_n*88)+30], c15_2);
    __m128 c15_3 = _mm_load_ss(&C[(l_n*88)+51]);
    __m128 a15_3 = _mm_load_ss(&A[84]);
    c15_3 = _mm_add_ss(c15_3, _mm_mul_ss(a15_3, b15));
    _mm_store_ss(&C[(l_n*88)+51], c15_3);
    __m128 c15_4 = _mm_load_ss(&C[(l_n*88)+79]);
    __m128 a15_4 = _mm_load_ss(&A[85]);
    c15_4 = _mm_add_ss(c15_4, _mm_mul_ss(a15_4, b15));
    _mm_store_ss(&C[(l_n*88)+79], c15_4);
#else
    C[(l_n*88)+5] += A[81] * B[(l_n*88)+15];
    C[(l_n*88)+15] += A[82] * B[(l_n*88)+15];
    C[(l_n*88)+30] += A[83] * B[(l_n*88)+15];
    C[(l_n*88)+51] += A[84] * B[(l_n*88)+15];
    C[(l_n*88)+79] += A[85] * B[(l_n*88)+15];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b16 = _mm_broadcast_ss(&B[(l_n*88)+16]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b16 = _mm_load_ss(&B[(l_n*88)+16]);    b16 = _mm_shuffle_ps(b16, b16, 0x00);
#endif
    __m128 c16_0 = _mm_load_ss(&C[(l_n*88)+6]);
    __m128 a16_0 = _mm_load_ss(&A[86]);
    c16_0 = _mm_add_ss(c16_0, _mm_mul_ss(a16_0, b16));
    _mm_store_ss(&C[(l_n*88)+6], c16_0);
    __m128 c16_1 = _mm_load_ss(&C[(l_n*88)+16]);
    __m128 a16_1 = _mm_load_ss(&A[87]);
    c16_1 = _mm_add_ss(c16_1, _mm_mul_ss(a16_1, b16));
    _mm_store_ss(&C[(l_n*88)+16], c16_1);
    __m128 c16_2 = _mm_load_ss(&C[(l_n*88)+31]);
    __m128 a16_2 = _mm_load_ss(&A[88]);
    c16_2 = _mm_add_ss(c16_2, _mm_mul_ss(a16_2, b16));
    _mm_store_ss(&C[(l_n*88)+31], c16_2);
    __m128 c16_3 = _mm_load_ss(&C[(l_n*88)+52]);
    __m128 a16_3 = _mm_load_ss(&A[89]);
    c16_3 = _mm_add_ss(c16_3, _mm_mul_ss(a16_3, b16));
    _mm_store_ss(&C[(l_n*88)+52], c16_3);
    __m128 c16_4 = _mm_load_ss(&C[(l_n*88)+80]);
    __m128 a16_4 = _mm_load_ss(&A[90]);
    c16_4 = _mm_add_ss(c16_4, _mm_mul_ss(a16_4, b16));
    _mm_store_ss(&C[(l_n*88)+80], c16_4);
#else
    C[(l_n*88)+6] += A[86] * B[(l_n*88)+16];
    C[(l_n*88)+16] += A[87] * B[(l_n*88)+16];
    C[(l_n*88)+31] += A[88] * B[(l_n*88)+16];
    C[(l_n*88)+52] += A[89] * B[(l_n*88)+16];
    C[(l_n*88)+80] += A[90] * B[(l_n*88)+16];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b17 = _mm_broadcast_ss(&B[(l_n*88)+17]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b17 = _mm_load_ss(&B[(l_n*88)+17]);    b17 = _mm_shuffle_ps(b17, b17, 0x00);
#endif
    __m128 c17_0 = _mm_load_ss(&C[(l_n*88)+1]);
    __m128 a17_0 = _mm_load_ss(&A[91]);
    c17_0 = _mm_add_ss(c17_0, _mm_mul_ss(a17_0, b17));
    _mm_store_ss(&C[(l_n*88)+1], c17_0);
    __m128 c17_1 = _mm_load_ss(&C[(l_n*88)+7]);
    __m128 a17_1 = _mm_load_ss(&A[92]);
    c17_1 = _mm_add_ss(c17_1, _mm_mul_ss(a17_1, b17));
    _mm_store_ss(&C[(l_n*88)+7], c17_1);
    __m128 c17_2 = _mm_load_ss(&C[(l_n*88)+17]);
    __m128 a17_2 = _mm_load_ss(&A[93]);
    c17_2 = _mm_add_ss(c17_2, _mm_mul_ss(a17_2, b17));
    _mm_store_ss(&C[(l_n*88)+17], c17_2);
    __m128 c17_3 = _mm_load_ss(&C[(l_n*88)+32]);
    __m128 a17_3 = _mm_load_ss(&A[94]);
    c17_3 = _mm_add_ss(c17_3, _mm_mul_ss(a17_3, b17));
    _mm_store_ss(&C[(l_n*88)+32], c17_3);
    __m128 c17_4 = _mm_load_ss(&C[(l_n*88)+53]);
    __m128 a17_4 = _mm_load_ss(&A[95]);
    c17_4 = _mm_add_ss(c17_4, _mm_mul_ss(a17_4, b17));
    _mm_store_ss(&C[(l_n*88)+53], c17_4);
    __m128 c17_5 = _mm_load_ss(&C[(l_n*88)+81]);
    __m128 a17_5 = _mm_load_ss(&A[96]);
    c17_5 = _mm_add_ss(c17_5, _mm_mul_ss(a17_5, b17));
    _mm_store_ss(&C[(l_n*88)+81], c17_5);
#else
    C[(l_n*88)+1] += A[91] * B[(l_n*88)+17];
    C[(l_n*88)+7] += A[92] * B[(l_n*88)+17];
    C[(l_n*88)+17] += A[93] * B[(l_n*88)+17];
    C[(l_n*88)+32] += A[94] * B[(l_n*88)+17];
    C[(l_n*88)+53] += A[95] * B[(l_n*88)+17];
    C[(l_n*88)+81] += A[96] * B[(l_n*88)+17];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b18 = _mm_broadcast_ss(&B[(l_n*88)+18]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b18 = _mm_load_ss(&B[(l_n*88)+18]);    b18 = _mm_shuffle_ps(b18, b18, 0x00);
#endif
    __m128 c18_0 = _mm_load_ss(&C[(l_n*88)+2]);
    __m128 a18_0 = _mm_load_ss(&A[97]);
    c18_0 = _mm_add_ss(c18_0, _mm_mul_ss(a18_0, b18));
    _mm_store_ss(&C[(l_n*88)+2], c18_0);
    __m128 c18_1 = _mm_load_ss(&C[(l_n*88)+8]);
    __m128 a18_1 = _mm_load_ss(&A[98]);
    c18_1 = _mm_add_ss(c18_1, _mm_mul_ss(a18_1, b18));
    _mm_store_ss(&C[(l_n*88)+8], c18_1);
    __m128 c18_2 = _mm_load_ss(&C[(l_n*88)+18]);
    __m128 a18_2 = _mm_load_ss(&A[99]);
    c18_2 = _mm_add_ss(c18_2, _mm_mul_ss(a18_2, b18));
    _mm_store_ss(&C[(l_n*88)+18], c18_2);
    __m128 c18_3 = _mm_load_ss(&C[(l_n*88)+33]);
    __m128 a18_3 = _mm_load_ss(&A[100]);
    c18_3 = _mm_add_ss(c18_3, _mm_mul_ss(a18_3, b18));
    _mm_store_ss(&C[(l_n*88)+33], c18_3);
    __m128 c18_4 = _mm_load_ss(&C[(l_n*88)+54]);
    __m128 a18_4 = _mm_load_ss(&A[101]);
    c18_4 = _mm_add_ss(c18_4, _mm_mul_ss(a18_4, b18));
    _mm_store_ss(&C[(l_n*88)+54], c18_4);
    __m128 c18_5 = _mm_load_ss(&C[(l_n*88)+82]);
    __m128 a18_5 = _mm_load_ss(&A[102]);
    c18_5 = _mm_add_ss(c18_5, _mm_mul_ss(a18_5, b18));
    _mm_store_ss(&C[(l_n*88)+82], c18_5);
#else
    C[(l_n*88)+2] += A[97] * B[(l_n*88)+18];
    C[(l_n*88)+8] += A[98] * B[(l_n*88)+18];
    C[(l_n*88)+18] += A[99] * B[(l_n*88)+18];
    C[(l_n*88)+33] += A[100] * B[(l_n*88)+18];
    C[(l_n*88)+54] += A[101] * B[(l_n*88)+18];
    C[(l_n*88)+82] += A[102] * B[(l_n*88)+18];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b19 = _mm_broadcast_ss(&B[(l_n*88)+19]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b19 = _mm_load_ss(&B[(l_n*88)+19]);    b19 = _mm_shuffle_ps(b19, b19, 0x00);
#endif
    __m128 c19_0 = _mm_load_ss(&C[(l_n*88)+0]);
    __m128 a19_0 = _mm_load_ss(&A[103]);
    c19_0 = _mm_add_ss(c19_0, _mm_mul_ss(a19_0, b19));
    _mm_store_ss(&C[(l_n*88)+0], c19_0);
    __m128 c19_1 = _mm_load_ss(&C[(l_n*88)+3]);
    __m128 a19_1 = _mm_load_ss(&A[104]);
    c19_1 = _mm_add_ss(c19_1, _mm_mul_ss(a19_1, b19));
    _mm_store_ss(&C[(l_n*88)+3], c19_1);
    __m128 c19_2 = _mm_load_ss(&C[(l_n*88)+9]);
    __m128 a19_2 = _mm_load_ss(&A[105]);
    c19_2 = _mm_add_ss(c19_2, _mm_mul_ss(a19_2, b19));
    _mm_store_ss(&C[(l_n*88)+9], c19_2);
    __m128 c19_3 = _mm_load_ss(&C[(l_n*88)+19]);
    __m128 a19_3 = _mm_load_ss(&A[106]);
    c19_3 = _mm_add_ss(c19_3, _mm_mul_ss(a19_3, b19));
    _mm_store_ss(&C[(l_n*88)+19], c19_3);
    __m128 c19_4 = _mm_load_ss(&C[(l_n*88)+34]);
    __m128 a19_4 = _mm_load_ss(&A[107]);
    c19_4 = _mm_add_ss(c19_4, _mm_mul_ss(a19_4, b19));
    _mm_store_ss(&C[(l_n*88)+34], c19_4);
    __m128 c19_5 = _mm_load_ss(&C[(l_n*88)+55]);
    __m128 a19_5 = _mm_load_ss(&A[108]);
    c19_5 = _mm_add_ss(c19_5, _mm_mul_ss(a19_5, b19));
    _mm_store_ss(&C[(l_n*88)+55], c19_5);
    __m128 c19_6 = _mm_load_ss(&C[(l_n*88)+83]);
    __m128 a19_6 = _mm_load_ss(&A[109]);
    c19_6 = _mm_add_ss(c19_6, _mm_mul_ss(a19_6, b19));
    _mm_store_ss(&C[(l_n*88)+83], c19_6);
#else
    C[(l_n*88)+0] += A[103] * B[(l_n*88)+19];
    C[(l_n*88)+3] += A[104] * B[(l_n*88)+19];
    C[(l_n*88)+9] += A[105] * B[(l_n*88)+19];
    C[(l_n*88)+19] += A[106] * B[(l_n*88)+19];
    C[(l_n*88)+34] += A[107] * B[(l_n*88)+19];
    C[(l_n*88)+55] += A[108] * B[(l_n*88)+19];
    C[(l_n*88)+83] += A[109] * B[(l_n*88)+19];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b20 = _mm_broadcast_ss(&B[(l_n*88)+20]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b20 = _mm_load_ss(&B[(l_n*88)+20]);    b20 = _mm_shuffle_ps(b20, b20, 0x00);
#endif
    __m128 c20_0 = _mm_load_ss(&C[(l_n*88)+20]);
    __m128 a20_0 = _mm_load_ss(&A[110]);
    c20_0 = _mm_add_ss(c20_0, _mm_mul_ss(a20_0, b20));
    _mm_store_ss(&C[(l_n*88)+20], c20_0);
    __m128 c20_1 = _mm_load_ss(&C[(l_n*88)+41]);
    __m128 a20_1 = _mm_load_ss(&A[111]);
    c20_1 = _mm_add_ss(c20_1, _mm_mul_ss(a20_1, b20));
    _mm_store_ss(&C[(l_n*88)+41], c20_1);
    __m128 c20_2 = _mm_load_ss(&C[(l_n*88)+69]);
    __m128 a20_2 = _mm_load_ss(&A[112]);
    c20_2 = _mm_add_ss(c20_2, _mm_mul_ss(a20_2, b20));
    _mm_store_ss(&C[(l_n*88)+69], c20_2);
#else
    C[(l_n*88)+20] += A[110] * B[(l_n*88)+20];
    C[(l_n*88)+41] += A[111] * B[(l_n*88)+20];
    C[(l_n*88)+69] += A[112] * B[(l_n*88)+20];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b21 = _mm_broadcast_ss(&B[(l_n*88)+21]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b21 = _mm_load_ss(&B[(l_n*88)+21]);    b21 = _mm_shuffle_ps(b21, b21, 0x00);
#endif
    __m128 c21_0 = _mm_load_ss(&C[(l_n*88)+21]);
    __m128 a21_0 = _mm_load_ss(&A[113]);
    c21_0 = _mm_add_ss(c21_0, _mm_mul_ss(a21_0, b21));
    _mm_store_ss(&C[(l_n*88)+21], c21_0);
    __m128 c21_1 = _mm_load_ss(&C[(l_n*88)+42]);
    __m128 a21_1 = _mm_load_ss(&A[114]);
    c21_1 = _mm_add_ss(c21_1, _mm_mul_ss(a21_1, b21));
    _mm_store_ss(&C[(l_n*88)+42], c21_1);
    __m128 c21_2 = _mm_load_ss(&C[(l_n*88)+70]);
    __m128 a21_2 = _mm_load_ss(&A[115]);
    c21_2 = _mm_add_ss(c21_2, _mm_mul_ss(a21_2, b21));
    _mm_store_ss(&C[(l_n*88)+70], c21_2);
#else
    C[(l_n*88)+21] += A[113] * B[(l_n*88)+21];
    C[(l_n*88)+42] += A[114] * B[(l_n*88)+21];
    C[(l_n*88)+70] += A[115] * B[(l_n*88)+21];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b22 = _mm_broadcast_ss(&B[(l_n*88)+22]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b22 = _mm_load_ss(&B[(l_n*88)+22]);    b22 = _mm_shuffle_ps(b22, b22, 0x00);
#endif
    __m128 c22_0 = _mm_load_ss(&C[(l_n*88)+22]);
    __m128 a22_0 = _mm_load_ss(&A[116]);
    c22_0 = _mm_add_ss(c22_0, _mm_mul_ss(a22_0, b22));
    _mm_store_ss(&C[(l_n*88)+22], c22_0);
    __m128 c22_1 = _mm_load_ss(&C[(l_n*88)+43]);
    __m128 a22_1 = _mm_load_ss(&A[117]);
    c22_1 = _mm_add_ss(c22_1, _mm_mul_ss(a22_1, b22));
    _mm_store_ss(&C[(l_n*88)+43], c22_1);
    __m128 c22_2 = _mm_load_ss(&C[(l_n*88)+71]);
    __m128 a22_2 = _mm_load_ss(&A[118]);
    c22_2 = _mm_add_ss(c22_2, _mm_mul_ss(a22_2, b22));
    _mm_store_ss(&C[(l_n*88)+71], c22_2);
#else
    C[(l_n*88)+22] += A[116] * B[(l_n*88)+22];
    C[(l_n*88)+43] += A[117] * B[(l_n*88)+22];
    C[(l_n*88)+71] += A[118] * B[(l_n*88)+22];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b23 = _mm_broadcast_ss(&B[(l_n*88)+23]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b23 = _mm_load_ss(&B[(l_n*88)+23]);    b23 = _mm_shuffle_ps(b23, b23, 0x00);
#endif
    __m128 c23_0 = _mm_load_ss(&C[(l_n*88)+23]);
    __m128 a23_0 = _mm_load_ss(&A[119]);
    c23_0 = _mm_add_ss(c23_0, _mm_mul_ss(a23_0, b23));
    _mm_store_ss(&C[(l_n*88)+23], c23_0);
    __m128 c23_1 = _mm_load_ss(&C[(l_n*88)+44]);
    __m128 a23_1 = _mm_load_ss(&A[120]);
    c23_1 = _mm_add_ss(c23_1, _mm_mul_ss(a23_1, b23));
    _mm_store_ss(&C[(l_n*88)+44], c23_1);
    __m128 c23_2 = _mm_load_ss(&C[(l_n*88)+72]);
    __m128 a23_2 = _mm_load_ss(&A[121]);
    c23_2 = _mm_add_ss(c23_2, _mm_mul_ss(a23_2, b23));
    _mm_store_ss(&C[(l_n*88)+72], c23_2);
#else
    C[(l_n*88)+23] += A[119] * B[(l_n*88)+23];
    C[(l_n*88)+44] += A[120] * B[(l_n*88)+23];
    C[(l_n*88)+72] += A[121] * B[(l_n*88)+23];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b24 = _mm_broadcast_ss(&B[(l_n*88)+24]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b24 = _mm_load_ss(&B[(l_n*88)+24]);    b24 = _mm_shuffle_ps(b24, b24, 0x00);
#endif
    __m128 c24_0 = _mm_load_ss(&C[(l_n*88)+24]);
    __m128 a24_0 = _mm_load_ss(&A[122]);
    c24_0 = _mm_add_ss(c24_0, _mm_mul_ss(a24_0, b24));
    _mm_store_ss(&C[(l_n*88)+24], c24_0);
    __m128 c24_1 = _mm_load_ss(&C[(l_n*88)+45]);
    __m128 a24_1 = _mm_load_ss(&A[123]);
    c24_1 = _mm_add_ss(c24_1, _mm_mul_ss(a24_1, b24));
    _mm_store_ss(&C[(l_n*88)+45], c24_1);
    __m128 c24_2 = _mm_load_ss(&C[(l_n*88)+73]);
    __m128 a24_2 = _mm_load_ss(&A[124]);
    c24_2 = _mm_add_ss(c24_2, _mm_mul_ss(a24_2, b24));
    _mm_store_ss(&C[(l_n*88)+73], c24_2);
#else
    C[(l_n*88)+24] += A[122] * B[(l_n*88)+24];
    C[(l_n*88)+45] += A[123] * B[(l_n*88)+24];
    C[(l_n*88)+73] += A[124] * B[(l_n*88)+24];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b25 = _mm_broadcast_ss(&B[(l_n*88)+25]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b25 = _mm_load_ss(&B[(l_n*88)+25]);    b25 = _mm_shuffle_ps(b25, b25, 0x00);
#endif
    __m128 c25_0 = _mm_load_ss(&C[(l_n*88)+10]);
    __m128 a25_0 = _mm_load_ss(&A[125]);
    c25_0 = _mm_add_ss(c25_0, _mm_mul_ss(a25_0, b25));
    _mm_store_ss(&C[(l_n*88)+10], c25_0);
    __m128 c25_1 = _mm_load_ss(&C[(l_n*88)+25]);
    __m128 a25_1 = _mm_load_ss(&A[126]);
    c25_1 = _mm_add_ss(c25_1, _mm_mul_ss(a25_1, b25));
    _mm_store_ss(&C[(l_n*88)+25], c25_1);
    __m128 c25_2 = _mm_load_ss(&C[(l_n*88)+46]);
    __m128 a25_2 = _mm_load_ss(&A[127]);
    c25_2 = _mm_add_ss(c25_2, _mm_mul_ss(a25_2, b25));
    _mm_store_ss(&C[(l_n*88)+46], c25_2);
    __m128 c25_3 = _mm_load_ss(&C[(l_n*88)+74]);
    __m128 a25_3 = _mm_load_ss(&A[128]);
    c25_3 = _mm_add_ss(c25_3, _mm_mul_ss(a25_3, b25));
    _mm_store_ss(&C[(l_n*88)+74], c25_3);
#else
    C[(l_n*88)+10] += A[125] * B[(l_n*88)+25];
    C[(l_n*88)+25] += A[126] * B[(l_n*88)+25];
    C[(l_n*88)+46] += A[127] * B[(l_n*88)+25];
    C[(l_n*88)+74] += A[128] * B[(l_n*88)+25];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b26 = _mm_broadcast_ss(&B[(l_n*88)+26]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b26 = _mm_load_ss(&B[(l_n*88)+26]);    b26 = _mm_shuffle_ps(b26, b26, 0x00);
#endif
    __m128 c26_0 = _mm_load_ss(&C[(l_n*88)+11]);
    __m128 a26_0 = _mm_load_ss(&A[129]);
    c26_0 = _mm_add_ss(c26_0, _mm_mul_ss(a26_0, b26));
    _mm_store_ss(&C[(l_n*88)+11], c26_0);
    __m128 c26_1 = _mm_load_ss(&C[(l_n*88)+26]);
    __m128 a26_1 = _mm_load_ss(&A[130]);
    c26_1 = _mm_add_ss(c26_1, _mm_mul_ss(a26_1, b26));
    _mm_store_ss(&C[(l_n*88)+26], c26_1);
    __m128 c26_2 = _mm_load_ss(&C[(l_n*88)+47]);
    __m128 a26_2 = _mm_load_ss(&A[131]);
    c26_2 = _mm_add_ss(c26_2, _mm_mul_ss(a26_2, b26));
    _mm_store_ss(&C[(l_n*88)+47], c26_2);
    __m128 c26_3 = _mm_load_ss(&C[(l_n*88)+75]);
    __m128 a26_3 = _mm_load_ss(&A[132]);
    c26_3 = _mm_add_ss(c26_3, _mm_mul_ss(a26_3, b26));
    _mm_store_ss(&C[(l_n*88)+75], c26_3);
#else
    C[(l_n*88)+11] += A[129] * B[(l_n*88)+26];
    C[(l_n*88)+26] += A[130] * B[(l_n*88)+26];
    C[(l_n*88)+47] += A[131] * B[(l_n*88)+26];
    C[(l_n*88)+75] += A[132] * B[(l_n*88)+26];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b27 = _mm_broadcast_ss(&B[(l_n*88)+27]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b27 = _mm_load_ss(&B[(l_n*88)+27]);    b27 = _mm_shuffle_ps(b27, b27, 0x00);
#endif
    __m128 c27_0 = _mm_load_ss(&C[(l_n*88)+12]);
    __m128 a27_0 = _mm_load_ss(&A[133]);
    c27_0 = _mm_add_ss(c27_0, _mm_mul_ss(a27_0, b27));
    _mm_store_ss(&C[(l_n*88)+12], c27_0);
    __m128 c27_1 = _mm_load_ss(&C[(l_n*88)+27]);
    __m128 a27_1 = _mm_load_ss(&A[134]);
    c27_1 = _mm_add_ss(c27_1, _mm_mul_ss(a27_1, b27));
    _mm_store_ss(&C[(l_n*88)+27], c27_1);
    __m128 c27_2 = _mm_load_ss(&C[(l_n*88)+48]);
    __m128 a27_2 = _mm_load_ss(&A[135]);
    c27_2 = _mm_add_ss(c27_2, _mm_mul_ss(a27_2, b27));
    _mm_store_ss(&C[(l_n*88)+48], c27_2);
    __m128 c27_3 = _mm_load_ss(&C[(l_n*88)+76]);
    __m128 a27_3 = _mm_load_ss(&A[136]);
    c27_3 = _mm_add_ss(c27_3, _mm_mul_ss(a27_3, b27));
    _mm_store_ss(&C[(l_n*88)+76], c27_3);
#else
    C[(l_n*88)+12] += A[133] * B[(l_n*88)+27];
    C[(l_n*88)+27] += A[134] * B[(l_n*88)+27];
    C[(l_n*88)+48] += A[135] * B[(l_n*88)+27];
    C[(l_n*88)+76] += A[136] * B[(l_n*88)+27];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b28 = _mm_broadcast_ss(&B[(l_n*88)+28]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b28 = _mm_load_ss(&B[(l_n*88)+28]);    b28 = _mm_shuffle_ps(b28, b28, 0x00);
#endif
    __m128 c28_0 = _mm_load_ss(&C[(l_n*88)+13]);
    __m128 a28_0 = _mm_load_ss(&A[137]);
    c28_0 = _mm_add_ss(c28_0, _mm_mul_ss(a28_0, b28));
    _mm_store_ss(&C[(l_n*88)+13], c28_0);
    __m128 c28_1 = _mm_load_ss(&C[(l_n*88)+28]);
    __m128 a28_1 = _mm_load_ss(&A[138]);
    c28_1 = _mm_add_ss(c28_1, _mm_mul_ss(a28_1, b28));
    _mm_store_ss(&C[(l_n*88)+28], c28_1);
    __m128 c28_2 = _mm_load_ss(&C[(l_n*88)+49]);
    __m128 a28_2 = _mm_load_ss(&A[139]);
    c28_2 = _mm_add_ss(c28_2, _mm_mul_ss(a28_2, b28));
    _mm_store_ss(&C[(l_n*88)+49], c28_2);
    __m128 c28_3 = _mm_load_ss(&C[(l_n*88)+77]);
    __m128 a28_3 = _mm_load_ss(&A[140]);
    c28_3 = _mm_add_ss(c28_3, _mm_mul_ss(a28_3, b28));
    _mm_store_ss(&C[(l_n*88)+77], c28_3);
#else
    C[(l_n*88)+13] += A[137] * B[(l_n*88)+28];
    C[(l_n*88)+28] += A[138] * B[(l_n*88)+28];
    C[(l_n*88)+49] += A[139] * B[(l_n*88)+28];
    C[(l_n*88)+77] += A[140] * B[(l_n*88)+28];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b29 = _mm_broadcast_ss(&B[(l_n*88)+29]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b29 = _mm_load_ss(&B[(l_n*88)+29]);    b29 = _mm_shuffle_ps(b29, b29, 0x00);
#endif
    __m128 c29_0 = _mm_load_ss(&C[(l_n*88)+4]);
    __m128 a29_0 = _mm_load_ss(&A[141]);
    c29_0 = _mm_add_ss(c29_0, _mm_mul_ss(a29_0, b29));
    _mm_store_ss(&C[(l_n*88)+4], c29_0);
    __m128 c29_1 = _mm_load_ss(&C[(l_n*88)+14]);
    __m128 a29_1 = _mm_load_ss(&A[142]);
    c29_1 = _mm_add_ss(c29_1, _mm_mul_ss(a29_1, b29));
    _mm_store_ss(&C[(l_n*88)+14], c29_1);
    __m128 c29_2 = _mm_load_ss(&C[(l_n*88)+29]);
    __m128 a29_2 = _mm_load_ss(&A[143]);
    c29_2 = _mm_add_ss(c29_2, _mm_mul_ss(a29_2, b29));
    _mm_store_ss(&C[(l_n*88)+29], c29_2);
    __m128 c29_3 = _mm_load_ss(&C[(l_n*88)+50]);
    __m128 a29_3 = _mm_load_ss(&A[144]);
    c29_3 = _mm_add_ss(c29_3, _mm_mul_ss(a29_3, b29));
    _mm_store_ss(&C[(l_n*88)+50], c29_3);
    __m128 c29_4 = _mm_load_ss(&C[(l_n*88)+78]);
    __m128 a29_4 = _mm_load_ss(&A[145]);
    c29_4 = _mm_add_ss(c29_4, _mm_mul_ss(a29_4, b29));
    _mm_store_ss(&C[(l_n*88)+78], c29_4);
#else
    C[(l_n*88)+4] += A[141] * B[(l_n*88)+29];
    C[(l_n*88)+14] += A[142] * B[(l_n*88)+29];
    C[(l_n*88)+29] += A[143] * B[(l_n*88)+29];
    C[(l_n*88)+50] += A[144] * B[(l_n*88)+29];
    C[(l_n*88)+78] += A[145] * B[(l_n*88)+29];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b30 = _mm_broadcast_ss(&B[(l_n*88)+30]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b30 = _mm_load_ss(&B[(l_n*88)+30]);    b30 = _mm_shuffle_ps(b30, b30, 0x00);
#endif
    __m128 c30_0 = _mm_load_ss(&C[(l_n*88)+5]);
    __m128 a30_0 = _mm_load_ss(&A[146]);
    c30_0 = _mm_add_ss(c30_0, _mm_mul_ss(a30_0, b30));
    _mm_store_ss(&C[(l_n*88)+5], c30_0);
    __m128 c30_1 = _mm_load_ss(&C[(l_n*88)+15]);
    __m128 a30_1 = _mm_load_ss(&A[147]);
    c30_1 = _mm_add_ss(c30_1, _mm_mul_ss(a30_1, b30));
    _mm_store_ss(&C[(l_n*88)+15], c30_1);
    __m128 c30_2 = _mm_load_ss(&C[(l_n*88)+30]);
    __m128 a30_2 = _mm_load_ss(&A[148]);
    c30_2 = _mm_add_ss(c30_2, _mm_mul_ss(a30_2, b30));
    _mm_store_ss(&C[(l_n*88)+30], c30_2);
    __m128 c30_3 = _mm_load_ss(&C[(l_n*88)+51]);
    __m128 a30_3 = _mm_load_ss(&A[149]);
    c30_3 = _mm_add_ss(c30_3, _mm_mul_ss(a30_3, b30));
    _mm_store_ss(&C[(l_n*88)+51], c30_3);
    __m128 c30_4 = _mm_load_ss(&C[(l_n*88)+79]);
    __m128 a30_4 = _mm_load_ss(&A[150]);
    c30_4 = _mm_add_ss(c30_4, _mm_mul_ss(a30_4, b30));
    _mm_store_ss(&C[(l_n*88)+79], c30_4);
#else
    C[(l_n*88)+5] += A[146] * B[(l_n*88)+30];
    C[(l_n*88)+15] += A[147] * B[(l_n*88)+30];
    C[(l_n*88)+30] += A[148] * B[(l_n*88)+30];
    C[(l_n*88)+51] += A[149] * B[(l_n*88)+30];
    C[(l_n*88)+79] += A[150] * B[(l_n*88)+30];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b31 = _mm_broadcast_ss(&B[(l_n*88)+31]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b31 = _mm_load_ss(&B[(l_n*88)+31]);    b31 = _mm_shuffle_ps(b31, b31, 0x00);
#endif
    __m128 c31_0 = _mm_load_ss(&C[(l_n*88)+6]);
    __m128 a31_0 = _mm_load_ss(&A[151]);
    c31_0 = _mm_add_ss(c31_0, _mm_mul_ss(a31_0, b31));
    _mm_store_ss(&C[(l_n*88)+6], c31_0);
    __m128 c31_1 = _mm_load_ss(&C[(l_n*88)+16]);
    __m128 a31_1 = _mm_load_ss(&A[152]);
    c31_1 = _mm_add_ss(c31_1, _mm_mul_ss(a31_1, b31));
    _mm_store_ss(&C[(l_n*88)+16], c31_1);
    __m128 c31_2 = _mm_load_ss(&C[(l_n*88)+31]);
    __m128 a31_2 = _mm_load_ss(&A[153]);
    c31_2 = _mm_add_ss(c31_2, _mm_mul_ss(a31_2, b31));
    _mm_store_ss(&C[(l_n*88)+31], c31_2);
    __m128 c31_3 = _mm_load_ss(&C[(l_n*88)+52]);
    __m128 a31_3 = _mm_load_ss(&A[154]);
    c31_3 = _mm_add_ss(c31_3, _mm_mul_ss(a31_3, b31));
    _mm_store_ss(&C[(l_n*88)+52], c31_3);
    __m128 c31_4 = _mm_load_ss(&C[(l_n*88)+80]);
    __m128 a31_4 = _mm_load_ss(&A[155]);
    c31_4 = _mm_add_ss(c31_4, _mm_mul_ss(a31_4, b31));
    _mm_store_ss(&C[(l_n*88)+80], c31_4);
#else
    C[(l_n*88)+6] += A[151] * B[(l_n*88)+31];
    C[(l_n*88)+16] += A[152] * B[(l_n*88)+31];
    C[(l_n*88)+31] += A[153] * B[(l_n*88)+31];
    C[(l_n*88)+52] += A[154] * B[(l_n*88)+31];
    C[(l_n*88)+80] += A[155] * B[(l_n*88)+31];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b32 = _mm_broadcast_ss(&B[(l_n*88)+32]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b32 = _mm_load_ss(&B[(l_n*88)+32]);    b32 = _mm_shuffle_ps(b32, b32, 0x00);
#endif
    __m128 c32_0 = _mm_load_ss(&C[(l_n*88)+1]);
    __m128 a32_0 = _mm_load_ss(&A[156]);
    c32_0 = _mm_add_ss(c32_0, _mm_mul_ss(a32_0, b32));
    _mm_store_ss(&C[(l_n*88)+1], c32_0);
    __m128 c32_1 = _mm_load_ss(&C[(l_n*88)+7]);
    __m128 a32_1 = _mm_load_ss(&A[157]);
    c32_1 = _mm_add_ss(c32_1, _mm_mul_ss(a32_1, b32));
    _mm_store_ss(&C[(l_n*88)+7], c32_1);
    __m128 c32_2 = _mm_load_ss(&C[(l_n*88)+17]);
    __m128 a32_2 = _mm_load_ss(&A[158]);
    c32_2 = _mm_add_ss(c32_2, _mm_mul_ss(a32_2, b32));
    _mm_store_ss(&C[(l_n*88)+17], c32_2);
    __m128 c32_3 = _mm_load_ss(&C[(l_n*88)+32]);
    __m128 a32_3 = _mm_load_ss(&A[159]);
    c32_3 = _mm_add_ss(c32_3, _mm_mul_ss(a32_3, b32));
    _mm_store_ss(&C[(l_n*88)+32], c32_3);
    __m128 c32_4 = _mm_load_ss(&C[(l_n*88)+53]);
    __m128 a32_4 = _mm_load_ss(&A[160]);
    c32_4 = _mm_add_ss(c32_4, _mm_mul_ss(a32_4, b32));
    _mm_store_ss(&C[(l_n*88)+53], c32_4);
    __m128 c32_5 = _mm_load_ss(&C[(l_n*88)+81]);
    __m128 a32_5 = _mm_load_ss(&A[161]);
    c32_5 = _mm_add_ss(c32_5, _mm_mul_ss(a32_5, b32));
    _mm_store_ss(&C[(l_n*88)+81], c32_5);
#else
    C[(l_n*88)+1] += A[156] * B[(l_n*88)+32];
    C[(l_n*88)+7] += A[157] * B[(l_n*88)+32];
    C[(l_n*88)+17] += A[158] * B[(l_n*88)+32];
    C[(l_n*88)+32] += A[159] * B[(l_n*88)+32];
    C[(l_n*88)+53] += A[160] * B[(l_n*88)+32];
    C[(l_n*88)+81] += A[161] * B[(l_n*88)+32];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b33 = _mm_broadcast_ss(&B[(l_n*88)+33]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b33 = _mm_load_ss(&B[(l_n*88)+33]);    b33 = _mm_shuffle_ps(b33, b33, 0x00);
#endif
    __m128 c33_0 = _mm_load_ss(&C[(l_n*88)+2]);
    __m128 a33_0 = _mm_load_ss(&A[162]);
    c33_0 = _mm_add_ss(c33_0, _mm_mul_ss(a33_0, b33));
    _mm_store_ss(&C[(l_n*88)+2], c33_0);
    __m128 c33_1 = _mm_load_ss(&C[(l_n*88)+8]);
    __m128 a33_1 = _mm_load_ss(&A[163]);
    c33_1 = _mm_add_ss(c33_1, _mm_mul_ss(a33_1, b33));
    _mm_store_ss(&C[(l_n*88)+8], c33_1);
    __m128 c33_2 = _mm_load_ss(&C[(l_n*88)+18]);
    __m128 a33_2 = _mm_load_ss(&A[164]);
    c33_2 = _mm_add_ss(c33_2, _mm_mul_ss(a33_2, b33));
    _mm_store_ss(&C[(l_n*88)+18], c33_2);
    __m128 c33_3 = _mm_load_ss(&C[(l_n*88)+33]);
    __m128 a33_3 = _mm_load_ss(&A[165]);
    c33_3 = _mm_add_ss(c33_3, _mm_mul_ss(a33_3, b33));
    _mm_store_ss(&C[(l_n*88)+33], c33_3);
    __m128 c33_4 = _mm_load_ss(&C[(l_n*88)+54]);
    __m128 a33_4 = _mm_load_ss(&A[166]);
    c33_4 = _mm_add_ss(c33_4, _mm_mul_ss(a33_4, b33));
    _mm_store_ss(&C[(l_n*88)+54], c33_4);
    __m128 c33_5 = _mm_load_ss(&C[(l_n*88)+82]);
    __m128 a33_5 = _mm_load_ss(&A[167]);
    c33_5 = _mm_add_ss(c33_5, _mm_mul_ss(a33_5, b33));
    _mm_store_ss(&C[(l_n*88)+82], c33_5);
#else
    C[(l_n*88)+2] += A[162] * B[(l_n*88)+33];
    C[(l_n*88)+8] += A[163] * B[(l_n*88)+33];
    C[(l_n*88)+18] += A[164] * B[(l_n*88)+33];
    C[(l_n*88)+33] += A[165] * B[(l_n*88)+33];
    C[(l_n*88)+54] += A[166] * B[(l_n*88)+33];
    C[(l_n*88)+82] += A[167] * B[(l_n*88)+33];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b34 = _mm_broadcast_ss(&B[(l_n*88)+34]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b34 = _mm_load_ss(&B[(l_n*88)+34]);    b34 = _mm_shuffle_ps(b34, b34, 0x00);
#endif
    __m128 c34_0 = _mm_load_ss(&C[(l_n*88)+0]);
    __m128 a34_0 = _mm_load_ss(&A[168]);
    c34_0 = _mm_add_ss(c34_0, _mm_mul_ss(a34_0, b34));
    _mm_store_ss(&C[(l_n*88)+0], c34_0);
    __m128 c34_1 = _mm_load_ss(&C[(l_n*88)+3]);
    __m128 a34_1 = _mm_load_ss(&A[169]);
    c34_1 = _mm_add_ss(c34_1, _mm_mul_ss(a34_1, b34));
    _mm_store_ss(&C[(l_n*88)+3], c34_1);
    __m128 c34_2 = _mm_load_ss(&C[(l_n*88)+9]);
    __m128 a34_2 = _mm_load_ss(&A[170]);
    c34_2 = _mm_add_ss(c34_2, _mm_mul_ss(a34_2, b34));
    _mm_store_ss(&C[(l_n*88)+9], c34_2);
    __m128 c34_3 = _mm_load_ss(&C[(l_n*88)+19]);
    __m128 a34_3 = _mm_load_ss(&A[171]);
    c34_3 = _mm_add_ss(c34_3, _mm_mul_ss(a34_3, b34));
    _mm_store_ss(&C[(l_n*88)+19], c34_3);
    __m128 c34_4 = _mm_load_ss(&C[(l_n*88)+34]);
    __m128 a34_4 = _mm_load_ss(&A[172]);
    c34_4 = _mm_add_ss(c34_4, _mm_mul_ss(a34_4, b34));
    _mm_store_ss(&C[(l_n*88)+34], c34_4);
    __m128 c34_5 = _mm_load_ss(&C[(l_n*88)+55]);
    __m128 a34_5 = _mm_load_ss(&A[173]);
    c34_5 = _mm_add_ss(c34_5, _mm_mul_ss(a34_5, b34));
    _mm_store_ss(&C[(l_n*88)+55], c34_5);
    __m128 c34_6 = _mm_load_ss(&C[(l_n*88)+83]);
    __m128 a34_6 = _mm_load_ss(&A[174]);
    c34_6 = _mm_add_ss(c34_6, _mm_mul_ss(a34_6, b34));
    _mm_store_ss(&C[(l_n*88)+83], c34_6);
#else
    C[(l_n*88)+0] += A[168] * B[(l_n*88)+34];
    C[(l_n*88)+3] += A[169] * B[(l_n*88)+34];
    C[(l_n*88)+9] += A[170] * B[(l_n*88)+34];
    C[(l_n*88)+19] += A[171] * B[(l_n*88)+34];
    C[(l_n*88)+34] += A[172] * B[(l_n*88)+34];
    C[(l_n*88)+55] += A[173] * B[(l_n*88)+34];
    C[(l_n*88)+83] += A[174] * B[(l_n*88)+34];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b35 = _mm_broadcast_ss(&B[(l_n*88)+35]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b35 = _mm_load_ss(&B[(l_n*88)+35]);    b35 = _mm_shuffle_ps(b35, b35, 0x00);
#endif
    __m128 c35_0 = _mm_load_ss(&C[(l_n*88)+35]);
    __m128 a35_0 = _mm_load_ss(&A[175]);
    c35_0 = _mm_add_ss(c35_0, _mm_mul_ss(a35_0, b35));
    _mm_store_ss(&C[(l_n*88)+35], c35_0);
    __m128 c35_1 = _mm_load_ss(&C[(l_n*88)+63]);
    __m128 a35_1 = _mm_load_ss(&A[176]);
    c35_1 = _mm_add_ss(c35_1, _mm_mul_ss(a35_1, b35));
    _mm_store_ss(&C[(l_n*88)+63], c35_1);
#else
    C[(l_n*88)+35] += A[175] * B[(l_n*88)+35];
    C[(l_n*88)+63] += A[176] * B[(l_n*88)+35];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b36 = _mm_broadcast_ss(&B[(l_n*88)+36]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b36 = _mm_load_ss(&B[(l_n*88)+36]);    b36 = _mm_shuffle_ps(b36, b36, 0x00);
#endif
    __m128 c36_0 = _mm_load_ss(&C[(l_n*88)+36]);
    __m128 a36_0 = _mm_load_ss(&A[177]);
    c36_0 = _mm_add_ss(c36_0, _mm_mul_ss(a36_0, b36));
    _mm_store_ss(&C[(l_n*88)+36], c36_0);
    __m128 c36_1 = _mm_load_ss(&C[(l_n*88)+64]);
    __m128 a36_1 = _mm_load_ss(&A[178]);
    c36_1 = _mm_add_ss(c36_1, _mm_mul_ss(a36_1, b36));
    _mm_store_ss(&C[(l_n*88)+64], c36_1);
#else
    C[(l_n*88)+36] += A[177] * B[(l_n*88)+36];
    C[(l_n*88)+64] += A[178] * B[(l_n*88)+36];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b37 = _mm_broadcast_ss(&B[(l_n*88)+37]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b37 = _mm_load_ss(&B[(l_n*88)+37]);    b37 = _mm_shuffle_ps(b37, b37, 0x00);
#endif
    __m128 c37_0 = _mm_load_ss(&C[(l_n*88)+37]);
    __m128 a37_0 = _mm_load_ss(&A[179]);
    c37_0 = _mm_add_ss(c37_0, _mm_mul_ss(a37_0, b37));
    _mm_store_ss(&C[(l_n*88)+37], c37_0);
    __m128 c37_1 = _mm_load_ss(&C[(l_n*88)+65]);
    __m128 a37_1 = _mm_load_ss(&A[180]);
    c37_1 = _mm_add_ss(c37_1, _mm_mul_ss(a37_1, b37));
    _mm_store_ss(&C[(l_n*88)+65], c37_1);
#else
    C[(l_n*88)+37] += A[179] * B[(l_n*88)+37];
    C[(l_n*88)+65] += A[180] * B[(l_n*88)+37];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b38 = _mm_broadcast_ss(&B[(l_n*88)+38]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b38 = _mm_load_ss(&B[(l_n*88)+38]);    b38 = _mm_shuffle_ps(b38, b38, 0x00);
#endif
    __m128 c38_0 = _mm_load_ss(&C[(l_n*88)+38]);
    __m128 a38_0 = _mm_load_ss(&A[181]);
    c38_0 = _mm_add_ss(c38_0, _mm_mul_ss(a38_0, b38));
    _mm_store_ss(&C[(l_n*88)+38], c38_0);
    __m128 c38_1 = _mm_load_ss(&C[(l_n*88)+66]);
    __m128 a38_1 = _mm_load_ss(&A[182]);
    c38_1 = _mm_add_ss(c38_1, _mm_mul_ss(a38_1, b38));
    _mm_store_ss(&C[(l_n*88)+66], c38_1);
#else
    C[(l_n*88)+38] += A[181] * B[(l_n*88)+38];
    C[(l_n*88)+66] += A[182] * B[(l_n*88)+38];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b39 = _mm_broadcast_ss(&B[(l_n*88)+39]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b39 = _mm_load_ss(&B[(l_n*88)+39]);    b39 = _mm_shuffle_ps(b39, b39, 0x00);
#endif
    __m128 c39_0 = _mm_load_ss(&C[(l_n*88)+39]);
    __m128 a39_0 = _mm_load_ss(&A[183]);
    c39_0 = _mm_add_ss(c39_0, _mm_mul_ss(a39_0, b39));
    _mm_store_ss(&C[(l_n*88)+39], c39_0);
    __m128 c39_1 = _mm_load_ss(&C[(l_n*88)+67]);
    __m128 a39_1 = _mm_load_ss(&A[184]);
    c39_1 = _mm_add_ss(c39_1, _mm_mul_ss(a39_1, b39));
    _mm_store_ss(&C[(l_n*88)+67], c39_1);
#else
    C[(l_n*88)+39] += A[183] * B[(l_n*88)+39];
    C[(l_n*88)+67] += A[184] * B[(l_n*88)+39];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b40 = _mm_broadcast_ss(&B[(l_n*88)+40]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b40 = _mm_load_ss(&B[(l_n*88)+40]);    b40 = _mm_shuffle_ps(b40, b40, 0x00);
#endif
    __m128 c40_0 = _mm_load_ss(&C[(l_n*88)+40]);
    __m128 a40_0 = _mm_load_ss(&A[185]);
    c40_0 = _mm_add_ss(c40_0, _mm_mul_ss(a40_0, b40));
    _mm_store_ss(&C[(l_n*88)+40], c40_0);
    __m128 c40_1 = _mm_load_ss(&C[(l_n*88)+68]);
    __m128 a40_1 = _mm_load_ss(&A[186]);
    c40_1 = _mm_add_ss(c40_1, _mm_mul_ss(a40_1, b40));
    _mm_store_ss(&C[(l_n*88)+68], c40_1);
#else
    C[(l_n*88)+40] += A[185] * B[(l_n*88)+40];
    C[(l_n*88)+68] += A[186] * B[(l_n*88)+40];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b41 = _mm_broadcast_ss(&B[(l_n*88)+41]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b41 = _mm_load_ss(&B[(l_n*88)+41]);    b41 = _mm_shuffle_ps(b41, b41, 0x00);
#endif
    __m128 c41_0 = _mm_load_ss(&C[(l_n*88)+20]);
    __m128 a41_0 = _mm_load_ss(&A[187]);
    c41_0 = _mm_add_ss(c41_0, _mm_mul_ss(a41_0, b41));
    _mm_store_ss(&C[(l_n*88)+20], c41_0);
    __m128 c41_1 = _mm_load_ss(&C[(l_n*88)+41]);
    __m128 a41_1 = _mm_load_ss(&A[188]);
    c41_1 = _mm_add_ss(c41_1, _mm_mul_ss(a41_1, b41));
    _mm_store_ss(&C[(l_n*88)+41], c41_1);
    __m128 c41_2 = _mm_load_ss(&C[(l_n*88)+69]);
    __m128 a41_2 = _mm_load_ss(&A[189]);
    c41_2 = _mm_add_ss(c41_2, _mm_mul_ss(a41_2, b41));
    _mm_store_ss(&C[(l_n*88)+69], c41_2);
#else
    C[(l_n*88)+20] += A[187] * B[(l_n*88)+41];
    C[(l_n*88)+41] += A[188] * B[(l_n*88)+41];
    C[(l_n*88)+69] += A[189] * B[(l_n*88)+41];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b42 = _mm_broadcast_ss(&B[(l_n*88)+42]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b42 = _mm_load_ss(&B[(l_n*88)+42]);    b42 = _mm_shuffle_ps(b42, b42, 0x00);
#endif
    __m128 c42_0 = _mm_load_ss(&C[(l_n*88)+21]);
    __m128 a42_0 = _mm_load_ss(&A[190]);
    c42_0 = _mm_add_ss(c42_0, _mm_mul_ss(a42_0, b42));
    _mm_store_ss(&C[(l_n*88)+21], c42_0);
    __m128 c42_1 = _mm_load_ss(&C[(l_n*88)+42]);
    __m128 a42_1 = _mm_load_ss(&A[191]);
    c42_1 = _mm_add_ss(c42_1, _mm_mul_ss(a42_1, b42));
    _mm_store_ss(&C[(l_n*88)+42], c42_1);
    __m128 c42_2 = _mm_load_ss(&C[(l_n*88)+70]);
    __m128 a42_2 = _mm_load_ss(&A[192]);
    c42_2 = _mm_add_ss(c42_2, _mm_mul_ss(a42_2, b42));
    _mm_store_ss(&C[(l_n*88)+70], c42_2);
#else
    C[(l_n*88)+21] += A[190] * B[(l_n*88)+42];
    C[(l_n*88)+42] += A[191] * B[(l_n*88)+42];
    C[(l_n*88)+70] += A[192] * B[(l_n*88)+42];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b43 = _mm_broadcast_ss(&B[(l_n*88)+43]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b43 = _mm_load_ss(&B[(l_n*88)+43]);    b43 = _mm_shuffle_ps(b43, b43, 0x00);
#endif
    __m128 c43_0 = _mm_load_ss(&C[(l_n*88)+22]);
    __m128 a43_0 = _mm_load_ss(&A[193]);
    c43_0 = _mm_add_ss(c43_0, _mm_mul_ss(a43_0, b43));
    _mm_store_ss(&C[(l_n*88)+22], c43_0);
    __m128 c43_1 = _mm_load_ss(&C[(l_n*88)+43]);
    __m128 a43_1 = _mm_load_ss(&A[194]);
    c43_1 = _mm_add_ss(c43_1, _mm_mul_ss(a43_1, b43));
    _mm_store_ss(&C[(l_n*88)+43], c43_1);
    __m128 c43_2 = _mm_load_ss(&C[(l_n*88)+71]);
    __m128 a43_2 = _mm_load_ss(&A[195]);
    c43_2 = _mm_add_ss(c43_2, _mm_mul_ss(a43_2, b43));
    _mm_store_ss(&C[(l_n*88)+71], c43_2);
#else
    C[(l_n*88)+22] += A[193] * B[(l_n*88)+43];
    C[(l_n*88)+43] += A[194] * B[(l_n*88)+43];
    C[(l_n*88)+71] += A[195] * B[(l_n*88)+43];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b44 = _mm_broadcast_ss(&B[(l_n*88)+44]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b44 = _mm_load_ss(&B[(l_n*88)+44]);    b44 = _mm_shuffle_ps(b44, b44, 0x00);
#endif
    __m128 c44_0 = _mm_load_ss(&C[(l_n*88)+23]);
    __m128 a44_0 = _mm_load_ss(&A[196]);
    c44_0 = _mm_add_ss(c44_0, _mm_mul_ss(a44_0, b44));
    _mm_store_ss(&C[(l_n*88)+23], c44_0);
    __m128 c44_1 = _mm_load_ss(&C[(l_n*88)+44]);
    __m128 a44_1 = _mm_load_ss(&A[197]);
    c44_1 = _mm_add_ss(c44_1, _mm_mul_ss(a44_1, b44));
    _mm_store_ss(&C[(l_n*88)+44], c44_1);
    __m128 c44_2 = _mm_load_ss(&C[(l_n*88)+72]);
    __m128 a44_2 = _mm_load_ss(&A[198]);
    c44_2 = _mm_add_ss(c44_2, _mm_mul_ss(a44_2, b44));
    _mm_store_ss(&C[(l_n*88)+72], c44_2);
#else
    C[(l_n*88)+23] += A[196] * B[(l_n*88)+44];
    C[(l_n*88)+44] += A[197] * B[(l_n*88)+44];
    C[(l_n*88)+72] += A[198] * B[(l_n*88)+44];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b45 = _mm_broadcast_ss(&B[(l_n*88)+45]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b45 = _mm_load_ss(&B[(l_n*88)+45]);    b45 = _mm_shuffle_ps(b45, b45, 0x00);
#endif
    __m128 c45_0 = _mm_load_ss(&C[(l_n*88)+24]);
    __m128 a45_0 = _mm_load_ss(&A[199]);
    c45_0 = _mm_add_ss(c45_0, _mm_mul_ss(a45_0, b45));
    _mm_store_ss(&C[(l_n*88)+24], c45_0);
    __m128 c45_1 = _mm_load_ss(&C[(l_n*88)+45]);
    __m128 a45_1 = _mm_load_ss(&A[200]);
    c45_1 = _mm_add_ss(c45_1, _mm_mul_ss(a45_1, b45));
    _mm_store_ss(&C[(l_n*88)+45], c45_1);
    __m128 c45_2 = _mm_load_ss(&C[(l_n*88)+73]);
    __m128 a45_2 = _mm_load_ss(&A[201]);
    c45_2 = _mm_add_ss(c45_2, _mm_mul_ss(a45_2, b45));
    _mm_store_ss(&C[(l_n*88)+73], c45_2);
#else
    C[(l_n*88)+24] += A[199] * B[(l_n*88)+45];
    C[(l_n*88)+45] += A[200] * B[(l_n*88)+45];
    C[(l_n*88)+73] += A[201] * B[(l_n*88)+45];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b46 = _mm_broadcast_ss(&B[(l_n*88)+46]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b46 = _mm_load_ss(&B[(l_n*88)+46]);    b46 = _mm_shuffle_ps(b46, b46, 0x00);
#endif
    __m128 c46_0 = _mm_load_ss(&C[(l_n*88)+10]);
    __m128 a46_0 = _mm_load_ss(&A[202]);
    c46_0 = _mm_add_ss(c46_0, _mm_mul_ss(a46_0, b46));
    _mm_store_ss(&C[(l_n*88)+10], c46_0);
    __m128 c46_1 = _mm_load_ss(&C[(l_n*88)+25]);
    __m128 a46_1 = _mm_load_ss(&A[203]);
    c46_1 = _mm_add_ss(c46_1, _mm_mul_ss(a46_1, b46));
    _mm_store_ss(&C[(l_n*88)+25], c46_1);
    __m128 c46_2 = _mm_load_ss(&C[(l_n*88)+46]);
    __m128 a46_2 = _mm_load_ss(&A[204]);
    c46_2 = _mm_add_ss(c46_2, _mm_mul_ss(a46_2, b46));
    _mm_store_ss(&C[(l_n*88)+46], c46_2);
    __m128 c46_3 = _mm_load_ss(&C[(l_n*88)+74]);
    __m128 a46_3 = _mm_load_ss(&A[205]);
    c46_3 = _mm_add_ss(c46_3, _mm_mul_ss(a46_3, b46));
    _mm_store_ss(&C[(l_n*88)+74], c46_3);
#else
    C[(l_n*88)+10] += A[202] * B[(l_n*88)+46];
    C[(l_n*88)+25] += A[203] * B[(l_n*88)+46];
    C[(l_n*88)+46] += A[204] * B[(l_n*88)+46];
    C[(l_n*88)+74] += A[205] * B[(l_n*88)+46];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b47 = _mm_broadcast_ss(&B[(l_n*88)+47]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b47 = _mm_load_ss(&B[(l_n*88)+47]);    b47 = _mm_shuffle_ps(b47, b47, 0x00);
#endif
    __m128 c47_0 = _mm_load_ss(&C[(l_n*88)+11]);
    __m128 a47_0 = _mm_load_ss(&A[206]);
    c47_0 = _mm_add_ss(c47_0, _mm_mul_ss(a47_0, b47));
    _mm_store_ss(&C[(l_n*88)+11], c47_0);
    __m128 c47_1 = _mm_load_ss(&C[(l_n*88)+26]);
    __m128 a47_1 = _mm_load_ss(&A[207]);
    c47_1 = _mm_add_ss(c47_1, _mm_mul_ss(a47_1, b47));
    _mm_store_ss(&C[(l_n*88)+26], c47_1);
    __m128 c47_2 = _mm_load_ss(&C[(l_n*88)+47]);
    __m128 a47_2 = _mm_load_ss(&A[208]);
    c47_2 = _mm_add_ss(c47_2, _mm_mul_ss(a47_2, b47));
    _mm_store_ss(&C[(l_n*88)+47], c47_2);
    __m128 c47_3 = _mm_load_ss(&C[(l_n*88)+75]);
    __m128 a47_3 = _mm_load_ss(&A[209]);
    c47_3 = _mm_add_ss(c47_3, _mm_mul_ss(a47_3, b47));
    _mm_store_ss(&C[(l_n*88)+75], c47_3);
#else
    C[(l_n*88)+11] += A[206] * B[(l_n*88)+47];
    C[(l_n*88)+26] += A[207] * B[(l_n*88)+47];
    C[(l_n*88)+47] += A[208] * B[(l_n*88)+47];
    C[(l_n*88)+75] += A[209] * B[(l_n*88)+47];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b48 = _mm_broadcast_ss(&B[(l_n*88)+48]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b48 = _mm_load_ss(&B[(l_n*88)+48]);    b48 = _mm_shuffle_ps(b48, b48, 0x00);
#endif
    __m128 c48_0 = _mm_load_ss(&C[(l_n*88)+12]);
    __m128 a48_0 = _mm_load_ss(&A[210]);
    c48_0 = _mm_add_ss(c48_0, _mm_mul_ss(a48_0, b48));
    _mm_store_ss(&C[(l_n*88)+12], c48_0);
    __m128 c48_1 = _mm_load_ss(&C[(l_n*88)+27]);
    __m128 a48_1 = _mm_load_ss(&A[211]);
    c48_1 = _mm_add_ss(c48_1, _mm_mul_ss(a48_1, b48));
    _mm_store_ss(&C[(l_n*88)+27], c48_1);
    __m128 c48_2 = _mm_load_ss(&C[(l_n*88)+48]);
    __m128 a48_2 = _mm_load_ss(&A[212]);
    c48_2 = _mm_add_ss(c48_2, _mm_mul_ss(a48_2, b48));
    _mm_store_ss(&C[(l_n*88)+48], c48_2);
    __m128 c48_3 = _mm_load_ss(&C[(l_n*88)+76]);
    __m128 a48_3 = _mm_load_ss(&A[213]);
    c48_3 = _mm_add_ss(c48_3, _mm_mul_ss(a48_3, b48));
    _mm_store_ss(&C[(l_n*88)+76], c48_3);
#else
    C[(l_n*88)+12] += A[210] * B[(l_n*88)+48];
    C[(l_n*88)+27] += A[211] * B[(l_n*88)+48];
    C[(l_n*88)+48] += A[212] * B[(l_n*88)+48];
    C[(l_n*88)+76] += A[213] * B[(l_n*88)+48];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b49 = _mm_broadcast_ss(&B[(l_n*88)+49]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b49 = _mm_load_ss(&B[(l_n*88)+49]);    b49 = _mm_shuffle_ps(b49, b49, 0x00);
#endif
    __m128 c49_0 = _mm_load_ss(&C[(l_n*88)+13]);
    __m128 a49_0 = _mm_load_ss(&A[214]);
    c49_0 = _mm_add_ss(c49_0, _mm_mul_ss(a49_0, b49));
    _mm_store_ss(&C[(l_n*88)+13], c49_0);
    __m128 c49_1 = _mm_load_ss(&C[(l_n*88)+28]);
    __m128 a49_1 = _mm_load_ss(&A[215]);
    c49_1 = _mm_add_ss(c49_1, _mm_mul_ss(a49_1, b49));
    _mm_store_ss(&C[(l_n*88)+28], c49_1);
    __m128 c49_2 = _mm_load_ss(&C[(l_n*88)+49]);
    __m128 a49_2 = _mm_load_ss(&A[216]);
    c49_2 = _mm_add_ss(c49_2, _mm_mul_ss(a49_2, b49));
    _mm_store_ss(&C[(l_n*88)+49], c49_2);
    __m128 c49_3 = _mm_load_ss(&C[(l_n*88)+77]);
    __m128 a49_3 = _mm_load_ss(&A[217]);
    c49_3 = _mm_add_ss(c49_3, _mm_mul_ss(a49_3, b49));
    _mm_store_ss(&C[(l_n*88)+77], c49_3);
#else
    C[(l_n*88)+13] += A[214] * B[(l_n*88)+49];
    C[(l_n*88)+28] += A[215] * B[(l_n*88)+49];
    C[(l_n*88)+49] += A[216] * B[(l_n*88)+49];
    C[(l_n*88)+77] += A[217] * B[(l_n*88)+49];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b50 = _mm_broadcast_ss(&B[(l_n*88)+50]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b50 = _mm_load_ss(&B[(l_n*88)+50]);    b50 = _mm_shuffle_ps(b50, b50, 0x00);
#endif
    __m128 c50_0 = _mm_load_ss(&C[(l_n*88)+4]);
    __m128 a50_0 = _mm_load_ss(&A[218]);
    c50_0 = _mm_add_ss(c50_0, _mm_mul_ss(a50_0, b50));
    _mm_store_ss(&C[(l_n*88)+4], c50_0);
    __m128 c50_1 = _mm_load_ss(&C[(l_n*88)+14]);
    __m128 a50_1 = _mm_load_ss(&A[219]);
    c50_1 = _mm_add_ss(c50_1, _mm_mul_ss(a50_1, b50));
    _mm_store_ss(&C[(l_n*88)+14], c50_1);
    __m128 c50_2 = _mm_load_ss(&C[(l_n*88)+29]);
    __m128 a50_2 = _mm_load_ss(&A[220]);
    c50_2 = _mm_add_ss(c50_2, _mm_mul_ss(a50_2, b50));
    _mm_store_ss(&C[(l_n*88)+29], c50_2);
    __m128 c50_3 = _mm_load_ss(&C[(l_n*88)+50]);
    __m128 a50_3 = _mm_load_ss(&A[221]);
    c50_3 = _mm_add_ss(c50_3, _mm_mul_ss(a50_3, b50));
    _mm_store_ss(&C[(l_n*88)+50], c50_3);
    __m128 c50_4 = _mm_load_ss(&C[(l_n*88)+78]);
    __m128 a50_4 = _mm_load_ss(&A[222]);
    c50_4 = _mm_add_ss(c50_4, _mm_mul_ss(a50_4, b50));
    _mm_store_ss(&C[(l_n*88)+78], c50_4);
#else
    C[(l_n*88)+4] += A[218] * B[(l_n*88)+50];
    C[(l_n*88)+14] += A[219] * B[(l_n*88)+50];
    C[(l_n*88)+29] += A[220] * B[(l_n*88)+50];
    C[(l_n*88)+50] += A[221] * B[(l_n*88)+50];
    C[(l_n*88)+78] += A[222] * B[(l_n*88)+50];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b51 = _mm_broadcast_ss(&B[(l_n*88)+51]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b51 = _mm_load_ss(&B[(l_n*88)+51]);    b51 = _mm_shuffle_ps(b51, b51, 0x00);
#endif
    __m128 c51_0 = _mm_load_ss(&C[(l_n*88)+5]);
    __m128 a51_0 = _mm_load_ss(&A[223]);
    c51_0 = _mm_add_ss(c51_0, _mm_mul_ss(a51_0, b51));
    _mm_store_ss(&C[(l_n*88)+5], c51_0);
    __m128 c51_1 = _mm_load_ss(&C[(l_n*88)+15]);
    __m128 a51_1 = _mm_load_ss(&A[224]);
    c51_1 = _mm_add_ss(c51_1, _mm_mul_ss(a51_1, b51));
    _mm_store_ss(&C[(l_n*88)+15], c51_1);
    __m128 c51_2 = _mm_load_ss(&C[(l_n*88)+30]);
    __m128 a51_2 = _mm_load_ss(&A[225]);
    c51_2 = _mm_add_ss(c51_2, _mm_mul_ss(a51_2, b51));
    _mm_store_ss(&C[(l_n*88)+30], c51_2);
    __m128 c51_3 = _mm_load_ss(&C[(l_n*88)+51]);
    __m128 a51_3 = _mm_load_ss(&A[226]);
    c51_3 = _mm_add_ss(c51_3, _mm_mul_ss(a51_3, b51));
    _mm_store_ss(&C[(l_n*88)+51], c51_3);
    __m128 c51_4 = _mm_load_ss(&C[(l_n*88)+79]);
    __m128 a51_4 = _mm_load_ss(&A[227]);
    c51_4 = _mm_add_ss(c51_4, _mm_mul_ss(a51_4, b51));
    _mm_store_ss(&C[(l_n*88)+79], c51_4);
#else
    C[(l_n*88)+5] += A[223] * B[(l_n*88)+51];
    C[(l_n*88)+15] += A[224] * B[(l_n*88)+51];
    C[(l_n*88)+30] += A[225] * B[(l_n*88)+51];
    C[(l_n*88)+51] += A[226] * B[(l_n*88)+51];
    C[(l_n*88)+79] += A[227] * B[(l_n*88)+51];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b52 = _mm_broadcast_ss(&B[(l_n*88)+52]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b52 = _mm_load_ss(&B[(l_n*88)+52]);    b52 = _mm_shuffle_ps(b52, b52, 0x00);
#endif
    __m128 c52_0 = _mm_load_ss(&C[(l_n*88)+6]);
    __m128 a52_0 = _mm_load_ss(&A[228]);
    c52_0 = _mm_add_ss(c52_0, _mm_mul_ss(a52_0, b52));
    _mm_store_ss(&C[(l_n*88)+6], c52_0);
    __m128 c52_1 = _mm_load_ss(&C[(l_n*88)+16]);
    __m128 a52_1 = _mm_load_ss(&A[229]);
    c52_1 = _mm_add_ss(c52_1, _mm_mul_ss(a52_1, b52));
    _mm_store_ss(&C[(l_n*88)+16], c52_1);
    __m128 c52_2 = _mm_load_ss(&C[(l_n*88)+31]);
    __m128 a52_2 = _mm_load_ss(&A[230]);
    c52_2 = _mm_add_ss(c52_2, _mm_mul_ss(a52_2, b52));
    _mm_store_ss(&C[(l_n*88)+31], c52_2);
    __m128 c52_3 = _mm_load_ss(&C[(l_n*88)+52]);
    __m128 a52_3 = _mm_load_ss(&A[231]);
    c52_3 = _mm_add_ss(c52_3, _mm_mul_ss(a52_3, b52));
    _mm_store_ss(&C[(l_n*88)+52], c52_3);
    __m128 c52_4 = _mm_load_ss(&C[(l_n*88)+80]);
    __m128 a52_4 = _mm_load_ss(&A[232]);
    c52_4 = _mm_add_ss(c52_4, _mm_mul_ss(a52_4, b52));
    _mm_store_ss(&C[(l_n*88)+80], c52_4);
#else
    C[(l_n*88)+6] += A[228] * B[(l_n*88)+52];
    C[(l_n*88)+16] += A[229] * B[(l_n*88)+52];
    C[(l_n*88)+31] += A[230] * B[(l_n*88)+52];
    C[(l_n*88)+52] += A[231] * B[(l_n*88)+52];
    C[(l_n*88)+80] += A[232] * B[(l_n*88)+52];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b53 = _mm_broadcast_ss(&B[(l_n*88)+53]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b53 = _mm_load_ss(&B[(l_n*88)+53]);    b53 = _mm_shuffle_ps(b53, b53, 0x00);
#endif
    __m128 c53_0 = _mm_load_ss(&C[(l_n*88)+1]);
    __m128 a53_0 = _mm_load_ss(&A[233]);
    c53_0 = _mm_add_ss(c53_0, _mm_mul_ss(a53_0, b53));
    _mm_store_ss(&C[(l_n*88)+1], c53_0);
    __m128 c53_1 = _mm_load_ss(&C[(l_n*88)+7]);
    __m128 a53_1 = _mm_load_ss(&A[234]);
    c53_1 = _mm_add_ss(c53_1, _mm_mul_ss(a53_1, b53));
    _mm_store_ss(&C[(l_n*88)+7], c53_1);
    __m128 c53_2 = _mm_load_ss(&C[(l_n*88)+17]);
    __m128 a53_2 = _mm_load_ss(&A[235]);
    c53_2 = _mm_add_ss(c53_2, _mm_mul_ss(a53_2, b53));
    _mm_store_ss(&C[(l_n*88)+17], c53_2);
    __m128 c53_3 = _mm_load_ss(&C[(l_n*88)+32]);
    __m128 a53_3 = _mm_load_ss(&A[236]);
    c53_3 = _mm_add_ss(c53_3, _mm_mul_ss(a53_3, b53));
    _mm_store_ss(&C[(l_n*88)+32], c53_3);
    __m128 c53_4 = _mm_load_ss(&C[(l_n*88)+53]);
    __m128 a53_4 = _mm_load_ss(&A[237]);
    c53_4 = _mm_add_ss(c53_4, _mm_mul_ss(a53_4, b53));
    _mm_store_ss(&C[(l_n*88)+53], c53_4);
    __m128 c53_5 = _mm_load_ss(&C[(l_n*88)+81]);
    __m128 a53_5 = _mm_load_ss(&A[238]);
    c53_5 = _mm_add_ss(c53_5, _mm_mul_ss(a53_5, b53));
    _mm_store_ss(&C[(l_n*88)+81], c53_5);
#else
    C[(l_n*88)+1] += A[233] * B[(l_n*88)+53];
    C[(l_n*88)+7] += A[234] * B[(l_n*88)+53];
    C[(l_n*88)+17] += A[235] * B[(l_n*88)+53];
    C[(l_n*88)+32] += A[236] * B[(l_n*88)+53];
    C[(l_n*88)+53] += A[237] * B[(l_n*88)+53];
    C[(l_n*88)+81] += A[238] * B[(l_n*88)+53];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b54 = _mm_broadcast_ss(&B[(l_n*88)+54]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b54 = _mm_load_ss(&B[(l_n*88)+54]);    b54 = _mm_shuffle_ps(b54, b54, 0x00);
#endif
    __m128 c54_0 = _mm_load_ss(&C[(l_n*88)+2]);
    __m128 a54_0 = _mm_load_ss(&A[239]);
    c54_0 = _mm_add_ss(c54_0, _mm_mul_ss(a54_0, b54));
    _mm_store_ss(&C[(l_n*88)+2], c54_0);
    __m128 c54_1 = _mm_load_ss(&C[(l_n*88)+8]);
    __m128 a54_1 = _mm_load_ss(&A[240]);
    c54_1 = _mm_add_ss(c54_1, _mm_mul_ss(a54_1, b54));
    _mm_store_ss(&C[(l_n*88)+8], c54_1);
    __m128 c54_2 = _mm_load_ss(&C[(l_n*88)+18]);
    __m128 a54_2 = _mm_load_ss(&A[241]);
    c54_2 = _mm_add_ss(c54_2, _mm_mul_ss(a54_2, b54));
    _mm_store_ss(&C[(l_n*88)+18], c54_2);
    __m128 c54_3 = _mm_load_ss(&C[(l_n*88)+33]);
    __m128 a54_3 = _mm_load_ss(&A[242]);
    c54_3 = _mm_add_ss(c54_3, _mm_mul_ss(a54_3, b54));
    _mm_store_ss(&C[(l_n*88)+33], c54_3);
    __m128 c54_4 = _mm_load_ss(&C[(l_n*88)+54]);
    __m128 a54_4 = _mm_load_ss(&A[243]);
    c54_4 = _mm_add_ss(c54_4, _mm_mul_ss(a54_4, b54));
    _mm_store_ss(&C[(l_n*88)+54], c54_4);
    __m128 c54_5 = _mm_load_ss(&C[(l_n*88)+82]);
    __m128 a54_5 = _mm_load_ss(&A[244]);
    c54_5 = _mm_add_ss(c54_5, _mm_mul_ss(a54_5, b54));
    _mm_store_ss(&C[(l_n*88)+82], c54_5);
#else
    C[(l_n*88)+2] += A[239] * B[(l_n*88)+54];
    C[(l_n*88)+8] += A[240] * B[(l_n*88)+54];
    C[(l_n*88)+18] += A[241] * B[(l_n*88)+54];
    C[(l_n*88)+33] += A[242] * B[(l_n*88)+54];
    C[(l_n*88)+54] += A[243] * B[(l_n*88)+54];
    C[(l_n*88)+82] += A[244] * B[(l_n*88)+54];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b55 = _mm_broadcast_ss(&B[(l_n*88)+55]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b55 = _mm_load_ss(&B[(l_n*88)+55]);    b55 = _mm_shuffle_ps(b55, b55, 0x00);
#endif
    __m128 c55_0 = _mm_load_ss(&C[(l_n*88)+0]);
    __m128 a55_0 = _mm_load_ss(&A[245]);
    c55_0 = _mm_add_ss(c55_0, _mm_mul_ss(a55_0, b55));
    _mm_store_ss(&C[(l_n*88)+0], c55_0);
    __m128 c55_1 = _mm_load_ss(&C[(l_n*88)+3]);
    __m128 a55_1 = _mm_load_ss(&A[246]);
    c55_1 = _mm_add_ss(c55_1, _mm_mul_ss(a55_1, b55));
    _mm_store_ss(&C[(l_n*88)+3], c55_1);
    __m128 c55_2 = _mm_load_ss(&C[(l_n*88)+9]);
    __m128 a55_2 = _mm_load_ss(&A[247]);
    c55_2 = _mm_add_ss(c55_2, _mm_mul_ss(a55_2, b55));
    _mm_store_ss(&C[(l_n*88)+9], c55_2);
    __m128 c55_3 = _mm_load_ss(&C[(l_n*88)+19]);
    __m128 a55_3 = _mm_load_ss(&A[248]);
    c55_3 = _mm_add_ss(c55_3, _mm_mul_ss(a55_3, b55));
    _mm_store_ss(&C[(l_n*88)+19], c55_3);
    __m128 c55_4 = _mm_load_ss(&C[(l_n*88)+34]);
    __m128 a55_4 = _mm_load_ss(&A[249]);
    c55_4 = _mm_add_ss(c55_4, _mm_mul_ss(a55_4, b55));
    _mm_store_ss(&C[(l_n*88)+34], c55_4);
    __m128 c55_5 = _mm_load_ss(&C[(l_n*88)+55]);
    __m128 a55_5 = _mm_load_ss(&A[250]);
    c55_5 = _mm_add_ss(c55_5, _mm_mul_ss(a55_5, b55));
    _mm_store_ss(&C[(l_n*88)+55], c55_5);
    __m128 c55_6 = _mm_load_ss(&C[(l_n*88)+83]);
    __m128 a55_6 = _mm_load_ss(&A[251]);
    c55_6 = _mm_add_ss(c55_6, _mm_mul_ss(a55_6, b55));
    _mm_store_ss(&C[(l_n*88)+83], c55_6);
#else
    C[(l_n*88)+0] += A[245] * B[(l_n*88)+55];
    C[(l_n*88)+3] += A[246] * B[(l_n*88)+55];
    C[(l_n*88)+9] += A[247] * B[(l_n*88)+55];
    C[(l_n*88)+19] += A[248] * B[(l_n*88)+55];
    C[(l_n*88)+34] += A[249] * B[(l_n*88)+55];
    C[(l_n*88)+55] += A[250] * B[(l_n*88)+55];
    C[(l_n*88)+83] += A[251] * B[(l_n*88)+55];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b56 = _mm_broadcast_ss(&B[(l_n*88)+56]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b56 = _mm_load_ss(&B[(l_n*88)+56]);    b56 = _mm_shuffle_ps(b56, b56, 0x00);
#endif
    __m128 c56_0 = _mm_load_ss(&C[(l_n*88)+56]);
    __m128 a56_0 = _mm_load_ss(&A[252]);
    c56_0 = _mm_add_ss(c56_0, _mm_mul_ss(a56_0, b56));
    _mm_store_ss(&C[(l_n*88)+56], c56_0);
#else
    C[(l_n*88)+56] += A[252] * B[(l_n*88)+56];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b57 = _mm_broadcast_ss(&B[(l_n*88)+57]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b57 = _mm_load_ss(&B[(l_n*88)+57]);    b57 = _mm_shuffle_ps(b57, b57, 0x00);
#endif
    __m128 c57_0 = _mm_load_ss(&C[(l_n*88)+57]);
    __m128 a57_0 = _mm_load_ss(&A[253]);
    c57_0 = _mm_add_ss(c57_0, _mm_mul_ss(a57_0, b57));
    _mm_store_ss(&C[(l_n*88)+57], c57_0);
#else
    C[(l_n*88)+57] += A[253] * B[(l_n*88)+57];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b58 = _mm_broadcast_ss(&B[(l_n*88)+58]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b58 = _mm_load_ss(&B[(l_n*88)+58]);    b58 = _mm_shuffle_ps(b58, b58, 0x00);
#endif
    __m128 c58_0 = _mm_load_ss(&C[(l_n*88)+58]);
    __m128 a58_0 = _mm_load_ss(&A[254]);
    c58_0 = _mm_add_ss(c58_0, _mm_mul_ss(a58_0, b58));
    _mm_store_ss(&C[(l_n*88)+58], c58_0);
#else
    C[(l_n*88)+58] += A[254] * B[(l_n*88)+58];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b59 = _mm_broadcast_ss(&B[(l_n*88)+59]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b59 = _mm_load_ss(&B[(l_n*88)+59]);    b59 = _mm_shuffle_ps(b59, b59, 0x00);
#endif
    __m128 c59_0 = _mm_load_ss(&C[(l_n*88)+59]);
    __m128 a59_0 = _mm_load_ss(&A[255]);
    c59_0 = _mm_add_ss(c59_0, _mm_mul_ss(a59_0, b59));
    _mm_store_ss(&C[(l_n*88)+59], c59_0);
#else
    C[(l_n*88)+59] += A[255] * B[(l_n*88)+59];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b60 = _mm_broadcast_ss(&B[(l_n*88)+60]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b60 = _mm_load_ss(&B[(l_n*88)+60]);    b60 = _mm_shuffle_ps(b60, b60, 0x00);
#endif
    __m128 c60_0 = _mm_load_ss(&C[(l_n*88)+60]);
    __m128 a60_0 = _mm_load_ss(&A[256]);
    c60_0 = _mm_add_ss(c60_0, _mm_mul_ss(a60_0, b60));
    _mm_store_ss(&C[(l_n*88)+60], c60_0);
#else
    C[(l_n*88)+60] += A[256] * B[(l_n*88)+60];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b61 = _mm_broadcast_ss(&B[(l_n*88)+61]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b61 = _mm_load_ss(&B[(l_n*88)+61]);    b61 = _mm_shuffle_ps(b61, b61, 0x00);
#endif
    __m128 c61_0 = _mm_load_ss(&C[(l_n*88)+61]);
    __m128 a61_0 = _mm_load_ss(&A[257]);
    c61_0 = _mm_add_ss(c61_0, _mm_mul_ss(a61_0, b61));
    _mm_store_ss(&C[(l_n*88)+61], c61_0);
#else
    C[(l_n*88)+61] += A[257] * B[(l_n*88)+61];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b62 = _mm_broadcast_ss(&B[(l_n*88)+62]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b62 = _mm_load_ss(&B[(l_n*88)+62]);    b62 = _mm_shuffle_ps(b62, b62, 0x00);
#endif
    __m128 c62_0 = _mm_load_ss(&C[(l_n*88)+62]);
    __m128 a62_0 = _mm_load_ss(&A[258]);
    c62_0 = _mm_add_ss(c62_0, _mm_mul_ss(a62_0, b62));
    _mm_store_ss(&C[(l_n*88)+62], c62_0);
#else
    C[(l_n*88)+62] += A[258] * B[(l_n*88)+62];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b63 = _mm_broadcast_ss(&B[(l_n*88)+63]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b63 = _mm_load_ss(&B[(l_n*88)+63]);    b63 = _mm_shuffle_ps(b63, b63, 0x00);
#endif
    __m128 c63_0 = _mm_load_ss(&C[(l_n*88)+35]);
    __m128 a63_0 = _mm_load_ss(&A[259]);
    c63_0 = _mm_add_ss(c63_0, _mm_mul_ss(a63_0, b63));
    _mm_store_ss(&C[(l_n*88)+35], c63_0);
    __m128 c63_1 = _mm_load_ss(&C[(l_n*88)+63]);
    __m128 a63_1 = _mm_load_ss(&A[260]);
    c63_1 = _mm_add_ss(c63_1, _mm_mul_ss(a63_1, b63));
    _mm_store_ss(&C[(l_n*88)+63], c63_1);
#else
    C[(l_n*88)+35] += A[259] * B[(l_n*88)+63];
    C[(l_n*88)+63] += A[260] * B[(l_n*88)+63];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b64 = _mm_broadcast_ss(&B[(l_n*88)+64]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b64 = _mm_load_ss(&B[(l_n*88)+64]);    b64 = _mm_shuffle_ps(b64, b64, 0x00);
#endif
    __m128 c64_0 = _mm_load_ss(&C[(l_n*88)+36]);
    __m128 a64_0 = _mm_load_ss(&A[261]);
    c64_0 = _mm_add_ss(c64_0, _mm_mul_ss(a64_0, b64));
    _mm_store_ss(&C[(l_n*88)+36], c64_0);
    __m128 c64_1 = _mm_load_ss(&C[(l_n*88)+64]);
    __m128 a64_1 = _mm_load_ss(&A[262]);
    c64_1 = _mm_add_ss(c64_1, _mm_mul_ss(a64_1, b64));
    _mm_store_ss(&C[(l_n*88)+64], c64_1);
#else
    C[(l_n*88)+36] += A[261] * B[(l_n*88)+64];
    C[(l_n*88)+64] += A[262] * B[(l_n*88)+64];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b65 = _mm_broadcast_ss(&B[(l_n*88)+65]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b65 = _mm_load_ss(&B[(l_n*88)+65]);    b65 = _mm_shuffle_ps(b65, b65, 0x00);
#endif
    __m128 c65_0 = _mm_load_ss(&C[(l_n*88)+37]);
    __m128 a65_0 = _mm_load_ss(&A[263]);
    c65_0 = _mm_add_ss(c65_0, _mm_mul_ss(a65_0, b65));
    _mm_store_ss(&C[(l_n*88)+37], c65_0);
    __m128 c65_1 = _mm_load_ss(&C[(l_n*88)+65]);
    __m128 a65_1 = _mm_load_ss(&A[264]);
    c65_1 = _mm_add_ss(c65_1, _mm_mul_ss(a65_1, b65));
    _mm_store_ss(&C[(l_n*88)+65], c65_1);
#else
    C[(l_n*88)+37] += A[263] * B[(l_n*88)+65];
    C[(l_n*88)+65] += A[264] * B[(l_n*88)+65];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b66 = _mm_broadcast_ss(&B[(l_n*88)+66]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b66 = _mm_load_ss(&B[(l_n*88)+66]);    b66 = _mm_shuffle_ps(b66, b66, 0x00);
#endif
    __m128 c66_0 = _mm_load_ss(&C[(l_n*88)+38]);
    __m128 a66_0 = _mm_load_ss(&A[265]);
    c66_0 = _mm_add_ss(c66_0, _mm_mul_ss(a66_0, b66));
    _mm_store_ss(&C[(l_n*88)+38], c66_0);
    __m128 c66_1 = _mm_load_ss(&C[(l_n*88)+66]);
    __m128 a66_1 = _mm_load_ss(&A[266]);
    c66_1 = _mm_add_ss(c66_1, _mm_mul_ss(a66_1, b66));
    _mm_store_ss(&C[(l_n*88)+66], c66_1);
#else
    C[(l_n*88)+38] += A[265] * B[(l_n*88)+66];
    C[(l_n*88)+66] += A[266] * B[(l_n*88)+66];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b67 = _mm_broadcast_ss(&B[(l_n*88)+67]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b67 = _mm_load_ss(&B[(l_n*88)+67]);    b67 = _mm_shuffle_ps(b67, b67, 0x00);
#endif
    __m128 c67_0 = _mm_load_ss(&C[(l_n*88)+39]);
    __m128 a67_0 = _mm_load_ss(&A[267]);
    c67_0 = _mm_add_ss(c67_0, _mm_mul_ss(a67_0, b67));
    _mm_store_ss(&C[(l_n*88)+39], c67_0);
    __m128 c67_1 = _mm_load_ss(&C[(l_n*88)+67]);
    __m128 a67_1 = _mm_load_ss(&A[268]);
    c67_1 = _mm_add_ss(c67_1, _mm_mul_ss(a67_1, b67));
    _mm_store_ss(&C[(l_n*88)+67], c67_1);
#else
    C[(l_n*88)+39] += A[267] * B[(l_n*88)+67];
    C[(l_n*88)+67] += A[268] * B[(l_n*88)+67];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b68 = _mm_broadcast_ss(&B[(l_n*88)+68]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b68 = _mm_load_ss(&B[(l_n*88)+68]);    b68 = _mm_shuffle_ps(b68, b68, 0x00);
#endif
    __m128 c68_0 = _mm_load_ss(&C[(l_n*88)+40]);
    __m128 a68_0 = _mm_load_ss(&A[269]);
    c68_0 = _mm_add_ss(c68_0, _mm_mul_ss(a68_0, b68));
    _mm_store_ss(&C[(l_n*88)+40], c68_0);
    __m128 c68_1 = _mm_load_ss(&C[(l_n*88)+68]);
    __m128 a68_1 = _mm_load_ss(&A[270]);
    c68_1 = _mm_add_ss(c68_1, _mm_mul_ss(a68_1, b68));
    _mm_store_ss(&C[(l_n*88)+68], c68_1);
#else
    C[(l_n*88)+40] += A[269] * B[(l_n*88)+68];
    C[(l_n*88)+68] += A[270] * B[(l_n*88)+68];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b69 = _mm_broadcast_ss(&B[(l_n*88)+69]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b69 = _mm_load_ss(&B[(l_n*88)+69]);    b69 = _mm_shuffle_ps(b69, b69, 0x00);
#endif
    __m128 c69_0 = _mm_load_ss(&C[(l_n*88)+20]);
    __m128 a69_0 = _mm_load_ss(&A[271]);
    c69_0 = _mm_add_ss(c69_0, _mm_mul_ss(a69_0, b69));
    _mm_store_ss(&C[(l_n*88)+20], c69_0);
    __m128 c69_1 = _mm_load_ss(&C[(l_n*88)+41]);
    __m128 a69_1 = _mm_load_ss(&A[272]);
    c69_1 = _mm_add_ss(c69_1, _mm_mul_ss(a69_1, b69));
    _mm_store_ss(&C[(l_n*88)+41], c69_1);
    __m128 c69_2 = _mm_load_ss(&C[(l_n*88)+69]);
    __m128 a69_2 = _mm_load_ss(&A[273]);
    c69_2 = _mm_add_ss(c69_2, _mm_mul_ss(a69_2, b69));
    _mm_store_ss(&C[(l_n*88)+69], c69_2);
#else
    C[(l_n*88)+20] += A[271] * B[(l_n*88)+69];
    C[(l_n*88)+41] += A[272] * B[(l_n*88)+69];
    C[(l_n*88)+69] += A[273] * B[(l_n*88)+69];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b70 = _mm_broadcast_ss(&B[(l_n*88)+70]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b70 = _mm_load_ss(&B[(l_n*88)+70]);    b70 = _mm_shuffle_ps(b70, b70, 0x00);
#endif
    __m128 c70_0 = _mm_load_ss(&C[(l_n*88)+21]);
    __m128 a70_0 = _mm_load_ss(&A[274]);
    c70_0 = _mm_add_ss(c70_0, _mm_mul_ss(a70_0, b70));
    _mm_store_ss(&C[(l_n*88)+21], c70_0);
    __m128 c70_1 = _mm_load_ss(&C[(l_n*88)+42]);
    __m128 a70_1 = _mm_load_ss(&A[275]);
    c70_1 = _mm_add_ss(c70_1, _mm_mul_ss(a70_1, b70));
    _mm_store_ss(&C[(l_n*88)+42], c70_1);
    __m128 c70_2 = _mm_load_ss(&C[(l_n*88)+70]);
    __m128 a70_2 = _mm_load_ss(&A[276]);
    c70_2 = _mm_add_ss(c70_2, _mm_mul_ss(a70_2, b70));
    _mm_store_ss(&C[(l_n*88)+70], c70_2);
#else
    C[(l_n*88)+21] += A[274] * B[(l_n*88)+70];
    C[(l_n*88)+42] += A[275] * B[(l_n*88)+70];
    C[(l_n*88)+70] += A[276] * B[(l_n*88)+70];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b71 = _mm_broadcast_ss(&B[(l_n*88)+71]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b71 = _mm_load_ss(&B[(l_n*88)+71]);    b71 = _mm_shuffle_ps(b71, b71, 0x00);
#endif
    __m128 c71_0 = _mm_load_ss(&C[(l_n*88)+22]);
    __m128 a71_0 = _mm_load_ss(&A[277]);
    c71_0 = _mm_add_ss(c71_0, _mm_mul_ss(a71_0, b71));
    _mm_store_ss(&C[(l_n*88)+22], c71_0);
    __m128 c71_1 = _mm_load_ss(&C[(l_n*88)+43]);
    __m128 a71_1 = _mm_load_ss(&A[278]);
    c71_1 = _mm_add_ss(c71_1, _mm_mul_ss(a71_1, b71));
    _mm_store_ss(&C[(l_n*88)+43], c71_1);
    __m128 c71_2 = _mm_load_ss(&C[(l_n*88)+71]);
    __m128 a71_2 = _mm_load_ss(&A[279]);
    c71_2 = _mm_add_ss(c71_2, _mm_mul_ss(a71_2, b71));
    _mm_store_ss(&C[(l_n*88)+71], c71_2);
#else
    C[(l_n*88)+22] += A[277] * B[(l_n*88)+71];
    C[(l_n*88)+43] += A[278] * B[(l_n*88)+71];
    C[(l_n*88)+71] += A[279] * B[(l_n*88)+71];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b72 = _mm_broadcast_ss(&B[(l_n*88)+72]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b72 = _mm_load_ss(&B[(l_n*88)+72]);    b72 = _mm_shuffle_ps(b72, b72, 0x00);
#endif
    __m128 c72_0 = _mm_load_ss(&C[(l_n*88)+23]);
    __m128 a72_0 = _mm_load_ss(&A[280]);
    c72_0 = _mm_add_ss(c72_0, _mm_mul_ss(a72_0, b72));
    _mm_store_ss(&C[(l_n*88)+23], c72_0);
    __m128 c72_1 = _mm_load_ss(&C[(l_n*88)+44]);
    __m128 a72_1 = _mm_load_ss(&A[281]);
    c72_1 = _mm_add_ss(c72_1, _mm_mul_ss(a72_1, b72));
    _mm_store_ss(&C[(l_n*88)+44], c72_1);
    __m128 c72_2 = _mm_load_ss(&C[(l_n*88)+72]);
    __m128 a72_2 = _mm_load_ss(&A[282]);
    c72_2 = _mm_add_ss(c72_2, _mm_mul_ss(a72_2, b72));
    _mm_store_ss(&C[(l_n*88)+72], c72_2);
#else
    C[(l_n*88)+23] += A[280] * B[(l_n*88)+72];
    C[(l_n*88)+44] += A[281] * B[(l_n*88)+72];
    C[(l_n*88)+72] += A[282] * B[(l_n*88)+72];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b73 = _mm_broadcast_ss(&B[(l_n*88)+73]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b73 = _mm_load_ss(&B[(l_n*88)+73]);    b73 = _mm_shuffle_ps(b73, b73, 0x00);
#endif
    __m128 c73_0 = _mm_load_ss(&C[(l_n*88)+24]);
    __m128 a73_0 = _mm_load_ss(&A[283]);
    c73_0 = _mm_add_ss(c73_0, _mm_mul_ss(a73_0, b73));
    _mm_store_ss(&C[(l_n*88)+24], c73_0);
    __m128 c73_1 = _mm_load_ss(&C[(l_n*88)+45]);
    __m128 a73_1 = _mm_load_ss(&A[284]);
    c73_1 = _mm_add_ss(c73_1, _mm_mul_ss(a73_1, b73));
    _mm_store_ss(&C[(l_n*88)+45], c73_1);
    __m128 c73_2 = _mm_load_ss(&C[(l_n*88)+73]);
    __m128 a73_2 = _mm_load_ss(&A[285]);
    c73_2 = _mm_add_ss(c73_2, _mm_mul_ss(a73_2, b73));
    _mm_store_ss(&C[(l_n*88)+73], c73_2);
#else
    C[(l_n*88)+24] += A[283] * B[(l_n*88)+73];
    C[(l_n*88)+45] += A[284] * B[(l_n*88)+73];
    C[(l_n*88)+73] += A[285] * B[(l_n*88)+73];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b74 = _mm_broadcast_ss(&B[(l_n*88)+74]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b74 = _mm_load_ss(&B[(l_n*88)+74]);    b74 = _mm_shuffle_ps(b74, b74, 0x00);
#endif
    __m128 c74_0 = _mm_load_ss(&C[(l_n*88)+10]);
    __m128 a74_0 = _mm_load_ss(&A[286]);
    c74_0 = _mm_add_ss(c74_0, _mm_mul_ss(a74_0, b74));
    _mm_store_ss(&C[(l_n*88)+10], c74_0);
    __m128 c74_1 = _mm_load_ss(&C[(l_n*88)+25]);
    __m128 a74_1 = _mm_load_ss(&A[287]);
    c74_1 = _mm_add_ss(c74_1, _mm_mul_ss(a74_1, b74));
    _mm_store_ss(&C[(l_n*88)+25], c74_1);
    __m128 c74_2 = _mm_load_ss(&C[(l_n*88)+46]);
    __m128 a74_2 = _mm_load_ss(&A[288]);
    c74_2 = _mm_add_ss(c74_2, _mm_mul_ss(a74_2, b74));
    _mm_store_ss(&C[(l_n*88)+46], c74_2);
    __m128 c74_3 = _mm_load_ss(&C[(l_n*88)+74]);
    __m128 a74_3 = _mm_load_ss(&A[289]);
    c74_3 = _mm_add_ss(c74_3, _mm_mul_ss(a74_3, b74));
    _mm_store_ss(&C[(l_n*88)+74], c74_3);
#else
    C[(l_n*88)+10] += A[286] * B[(l_n*88)+74];
    C[(l_n*88)+25] += A[287] * B[(l_n*88)+74];
    C[(l_n*88)+46] += A[288] * B[(l_n*88)+74];
    C[(l_n*88)+74] += A[289] * B[(l_n*88)+74];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b75 = _mm_broadcast_ss(&B[(l_n*88)+75]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b75 = _mm_load_ss(&B[(l_n*88)+75]);    b75 = _mm_shuffle_ps(b75, b75, 0x00);
#endif
    __m128 c75_0 = _mm_load_ss(&C[(l_n*88)+11]);
    __m128 a75_0 = _mm_load_ss(&A[290]);
    c75_0 = _mm_add_ss(c75_0, _mm_mul_ss(a75_0, b75));
    _mm_store_ss(&C[(l_n*88)+11], c75_0);
    __m128 c75_1 = _mm_load_ss(&C[(l_n*88)+26]);
    __m128 a75_1 = _mm_load_ss(&A[291]);
    c75_1 = _mm_add_ss(c75_1, _mm_mul_ss(a75_1, b75));
    _mm_store_ss(&C[(l_n*88)+26], c75_1);
    __m128 c75_2 = _mm_load_ss(&C[(l_n*88)+47]);
    __m128 a75_2 = _mm_load_ss(&A[292]);
    c75_2 = _mm_add_ss(c75_2, _mm_mul_ss(a75_2, b75));
    _mm_store_ss(&C[(l_n*88)+47], c75_2);
    __m128 c75_3 = _mm_load_ss(&C[(l_n*88)+75]);
    __m128 a75_3 = _mm_load_ss(&A[293]);
    c75_3 = _mm_add_ss(c75_3, _mm_mul_ss(a75_3, b75));
    _mm_store_ss(&C[(l_n*88)+75], c75_3);
#else
    C[(l_n*88)+11] += A[290] * B[(l_n*88)+75];
    C[(l_n*88)+26] += A[291] * B[(l_n*88)+75];
    C[(l_n*88)+47] += A[292] * B[(l_n*88)+75];
    C[(l_n*88)+75] += A[293] * B[(l_n*88)+75];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b76 = _mm_broadcast_ss(&B[(l_n*88)+76]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b76 = _mm_load_ss(&B[(l_n*88)+76]);    b76 = _mm_shuffle_ps(b76, b76, 0x00);
#endif
    __m128 c76_0 = _mm_load_ss(&C[(l_n*88)+12]);
    __m128 a76_0 = _mm_load_ss(&A[294]);
    c76_0 = _mm_add_ss(c76_0, _mm_mul_ss(a76_0, b76));
    _mm_store_ss(&C[(l_n*88)+12], c76_0);
    __m128 c76_1 = _mm_load_ss(&C[(l_n*88)+27]);
    __m128 a76_1 = _mm_load_ss(&A[295]);
    c76_1 = _mm_add_ss(c76_1, _mm_mul_ss(a76_1, b76));
    _mm_store_ss(&C[(l_n*88)+27], c76_1);
    __m128 c76_2 = _mm_load_ss(&C[(l_n*88)+48]);
    __m128 a76_2 = _mm_load_ss(&A[296]);
    c76_2 = _mm_add_ss(c76_2, _mm_mul_ss(a76_2, b76));
    _mm_store_ss(&C[(l_n*88)+48], c76_2);
    __m128 c76_3 = _mm_load_ss(&C[(l_n*88)+76]);
    __m128 a76_3 = _mm_load_ss(&A[297]);
    c76_3 = _mm_add_ss(c76_3, _mm_mul_ss(a76_3, b76));
    _mm_store_ss(&C[(l_n*88)+76], c76_3);
#else
    C[(l_n*88)+12] += A[294] * B[(l_n*88)+76];
    C[(l_n*88)+27] += A[295] * B[(l_n*88)+76];
    C[(l_n*88)+48] += A[296] * B[(l_n*88)+76];
    C[(l_n*88)+76] += A[297] * B[(l_n*88)+76];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b77 = _mm_broadcast_ss(&B[(l_n*88)+77]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b77 = _mm_load_ss(&B[(l_n*88)+77]);    b77 = _mm_shuffle_ps(b77, b77, 0x00);
#endif
    __m128 c77_0 = _mm_load_ss(&C[(l_n*88)+13]);
    __m128 a77_0 = _mm_load_ss(&A[298]);
    c77_0 = _mm_add_ss(c77_0, _mm_mul_ss(a77_0, b77));
    _mm_store_ss(&C[(l_n*88)+13], c77_0);
    __m128 c77_1 = _mm_load_ss(&C[(l_n*88)+28]);
    __m128 a77_1 = _mm_load_ss(&A[299]);
    c77_1 = _mm_add_ss(c77_1, _mm_mul_ss(a77_1, b77));
    _mm_store_ss(&C[(l_n*88)+28], c77_1);
    __m128 c77_2 = _mm_load_ss(&C[(l_n*88)+49]);
    __m128 a77_2 = _mm_load_ss(&A[300]);
    c77_2 = _mm_add_ss(c77_2, _mm_mul_ss(a77_2, b77));
    _mm_store_ss(&C[(l_n*88)+49], c77_2);
    __m128 c77_3 = _mm_load_ss(&C[(l_n*88)+77]);
    __m128 a77_3 = _mm_load_ss(&A[301]);
    c77_3 = _mm_add_ss(c77_3, _mm_mul_ss(a77_3, b77));
    _mm_store_ss(&C[(l_n*88)+77], c77_3);
#else
    C[(l_n*88)+13] += A[298] * B[(l_n*88)+77];
    C[(l_n*88)+28] += A[299] * B[(l_n*88)+77];
    C[(l_n*88)+49] += A[300] * B[(l_n*88)+77];
    C[(l_n*88)+77] += A[301] * B[(l_n*88)+77];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b78 = _mm_broadcast_ss(&B[(l_n*88)+78]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b78 = _mm_load_ss(&B[(l_n*88)+78]);    b78 = _mm_shuffle_ps(b78, b78, 0x00);
#endif
    __m128 c78_0 = _mm_load_ss(&C[(l_n*88)+4]);
    __m128 a78_0 = _mm_load_ss(&A[302]);
    c78_0 = _mm_add_ss(c78_0, _mm_mul_ss(a78_0, b78));
    _mm_store_ss(&C[(l_n*88)+4], c78_0);
    __m128 c78_1 = _mm_load_ss(&C[(l_n*88)+14]);
    __m128 a78_1 = _mm_load_ss(&A[303]);
    c78_1 = _mm_add_ss(c78_1, _mm_mul_ss(a78_1, b78));
    _mm_store_ss(&C[(l_n*88)+14], c78_1);
    __m128 c78_2 = _mm_load_ss(&C[(l_n*88)+29]);
    __m128 a78_2 = _mm_load_ss(&A[304]);
    c78_2 = _mm_add_ss(c78_2, _mm_mul_ss(a78_2, b78));
    _mm_store_ss(&C[(l_n*88)+29], c78_2);
    __m128 c78_3 = _mm_load_ss(&C[(l_n*88)+50]);
    __m128 a78_3 = _mm_load_ss(&A[305]);
    c78_3 = _mm_add_ss(c78_3, _mm_mul_ss(a78_3, b78));
    _mm_store_ss(&C[(l_n*88)+50], c78_3);
    __m128 c78_4 = _mm_load_ss(&C[(l_n*88)+78]);
    __m128 a78_4 = _mm_load_ss(&A[306]);
    c78_4 = _mm_add_ss(c78_4, _mm_mul_ss(a78_4, b78));
    _mm_store_ss(&C[(l_n*88)+78], c78_4);
#else
    C[(l_n*88)+4] += A[302] * B[(l_n*88)+78];
    C[(l_n*88)+14] += A[303] * B[(l_n*88)+78];
    C[(l_n*88)+29] += A[304] * B[(l_n*88)+78];
    C[(l_n*88)+50] += A[305] * B[(l_n*88)+78];
    C[(l_n*88)+78] += A[306] * B[(l_n*88)+78];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b79 = _mm_broadcast_ss(&B[(l_n*88)+79]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b79 = _mm_load_ss(&B[(l_n*88)+79]);    b79 = _mm_shuffle_ps(b79, b79, 0x00);
#endif
    __m128 c79_0 = _mm_load_ss(&C[(l_n*88)+5]);
    __m128 a79_0 = _mm_load_ss(&A[307]);
    c79_0 = _mm_add_ss(c79_0, _mm_mul_ss(a79_0, b79));
    _mm_store_ss(&C[(l_n*88)+5], c79_0);
    __m128 c79_1 = _mm_load_ss(&C[(l_n*88)+15]);
    __m128 a79_1 = _mm_load_ss(&A[308]);
    c79_1 = _mm_add_ss(c79_1, _mm_mul_ss(a79_1, b79));
    _mm_store_ss(&C[(l_n*88)+15], c79_1);
    __m128 c79_2 = _mm_load_ss(&C[(l_n*88)+30]);
    __m128 a79_2 = _mm_load_ss(&A[309]);
    c79_2 = _mm_add_ss(c79_2, _mm_mul_ss(a79_2, b79));
    _mm_store_ss(&C[(l_n*88)+30], c79_2);
    __m128 c79_3 = _mm_load_ss(&C[(l_n*88)+51]);
    __m128 a79_3 = _mm_load_ss(&A[310]);
    c79_3 = _mm_add_ss(c79_3, _mm_mul_ss(a79_3, b79));
    _mm_store_ss(&C[(l_n*88)+51], c79_3);
    __m128 c79_4 = _mm_load_ss(&C[(l_n*88)+79]);
    __m128 a79_4 = _mm_load_ss(&A[311]);
    c79_4 = _mm_add_ss(c79_4, _mm_mul_ss(a79_4, b79));
    _mm_store_ss(&C[(l_n*88)+79], c79_4);
#else
    C[(l_n*88)+5] += A[307] * B[(l_n*88)+79];
    C[(l_n*88)+15] += A[308] * B[(l_n*88)+79];
    C[(l_n*88)+30] += A[309] * B[(l_n*88)+79];
    C[(l_n*88)+51] += A[310] * B[(l_n*88)+79];
    C[(l_n*88)+79] += A[311] * B[(l_n*88)+79];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b80 = _mm_broadcast_ss(&B[(l_n*88)+80]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b80 = _mm_load_ss(&B[(l_n*88)+80]);    b80 = _mm_shuffle_ps(b80, b80, 0x00);
#endif
    __m128 c80_0 = _mm_load_ss(&C[(l_n*88)+6]);
    __m128 a80_0 = _mm_load_ss(&A[312]);
    c80_0 = _mm_add_ss(c80_0, _mm_mul_ss(a80_0, b80));
    _mm_store_ss(&C[(l_n*88)+6], c80_0);
    __m128 c80_1 = _mm_load_ss(&C[(l_n*88)+16]);
    __m128 a80_1 = _mm_load_ss(&A[313]);
    c80_1 = _mm_add_ss(c80_1, _mm_mul_ss(a80_1, b80));
    _mm_store_ss(&C[(l_n*88)+16], c80_1);
    __m128 c80_2 = _mm_load_ss(&C[(l_n*88)+31]);
    __m128 a80_2 = _mm_load_ss(&A[314]);
    c80_2 = _mm_add_ss(c80_2, _mm_mul_ss(a80_2, b80));
    _mm_store_ss(&C[(l_n*88)+31], c80_2);
    __m128 c80_3 = _mm_load_ss(&C[(l_n*88)+52]);
    __m128 a80_3 = _mm_load_ss(&A[315]);
    c80_3 = _mm_add_ss(c80_3, _mm_mul_ss(a80_3, b80));
    _mm_store_ss(&C[(l_n*88)+52], c80_3);
    __m128 c80_4 = _mm_load_ss(&C[(l_n*88)+80]);
    __m128 a80_4 = _mm_load_ss(&A[316]);
    c80_4 = _mm_add_ss(c80_4, _mm_mul_ss(a80_4, b80));
    _mm_store_ss(&C[(l_n*88)+80], c80_4);
#else
    C[(l_n*88)+6] += A[312] * B[(l_n*88)+80];
    C[(l_n*88)+16] += A[313] * B[(l_n*88)+80];
    C[(l_n*88)+31] += A[314] * B[(l_n*88)+80];
    C[(l_n*88)+52] += A[315] * B[(l_n*88)+80];
    C[(l_n*88)+80] += A[316] * B[(l_n*88)+80];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b81 = _mm_broadcast_ss(&B[(l_n*88)+81]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b81 = _mm_load_ss(&B[(l_n*88)+81]);    b81 = _mm_shuffle_ps(b81, b81, 0x00);
#endif
    __m128 c81_0 = _mm_load_ss(&C[(l_n*88)+1]);
    __m128 a81_0 = _mm_load_ss(&A[317]);
    c81_0 = _mm_add_ss(c81_0, _mm_mul_ss(a81_0, b81));
    _mm_store_ss(&C[(l_n*88)+1], c81_0);
    __m128 c81_1 = _mm_load_ss(&C[(l_n*88)+7]);
    __m128 a81_1 = _mm_load_ss(&A[318]);
    c81_1 = _mm_add_ss(c81_1, _mm_mul_ss(a81_1, b81));
    _mm_store_ss(&C[(l_n*88)+7], c81_1);
    __m128 c81_2 = _mm_load_ss(&C[(l_n*88)+17]);
    __m128 a81_2 = _mm_load_ss(&A[319]);
    c81_2 = _mm_add_ss(c81_2, _mm_mul_ss(a81_2, b81));
    _mm_store_ss(&C[(l_n*88)+17], c81_2);
    __m128 c81_3 = _mm_load_ss(&C[(l_n*88)+32]);
    __m128 a81_3 = _mm_load_ss(&A[320]);
    c81_3 = _mm_add_ss(c81_3, _mm_mul_ss(a81_3, b81));
    _mm_store_ss(&C[(l_n*88)+32], c81_3);
    __m128 c81_4 = _mm_load_ss(&C[(l_n*88)+53]);
    __m128 a81_4 = _mm_load_ss(&A[321]);
    c81_4 = _mm_add_ss(c81_4, _mm_mul_ss(a81_4, b81));
    _mm_store_ss(&C[(l_n*88)+53], c81_4);
    __m128 c81_5 = _mm_load_ss(&C[(l_n*88)+81]);
    __m128 a81_5 = _mm_load_ss(&A[322]);
    c81_5 = _mm_add_ss(c81_5, _mm_mul_ss(a81_5, b81));
    _mm_store_ss(&C[(l_n*88)+81], c81_5);
#else
    C[(l_n*88)+1] += A[317] * B[(l_n*88)+81];
    C[(l_n*88)+7] += A[318] * B[(l_n*88)+81];
    C[(l_n*88)+17] += A[319] * B[(l_n*88)+81];
    C[(l_n*88)+32] += A[320] * B[(l_n*88)+81];
    C[(l_n*88)+53] += A[321] * B[(l_n*88)+81];
    C[(l_n*88)+81] += A[322] * B[(l_n*88)+81];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b82 = _mm_broadcast_ss(&B[(l_n*88)+82]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b82 = _mm_load_ss(&B[(l_n*88)+82]);    b82 = _mm_shuffle_ps(b82, b82, 0x00);
#endif
    __m128 c82_0 = _mm_load_ss(&C[(l_n*88)+2]);
    __m128 a82_0 = _mm_load_ss(&A[323]);
    c82_0 = _mm_add_ss(c82_0, _mm_mul_ss(a82_0, b82));
    _mm_store_ss(&C[(l_n*88)+2], c82_0);
    __m128 c82_1 = _mm_load_ss(&C[(l_n*88)+8]);
    __m128 a82_1 = _mm_load_ss(&A[324]);
    c82_1 = _mm_add_ss(c82_1, _mm_mul_ss(a82_1, b82));
    _mm_store_ss(&C[(l_n*88)+8], c82_1);
    __m128 c82_2 = _mm_load_ss(&C[(l_n*88)+18]);
    __m128 a82_2 = _mm_load_ss(&A[325]);
    c82_2 = _mm_add_ss(c82_2, _mm_mul_ss(a82_2, b82));
    _mm_store_ss(&C[(l_n*88)+18], c82_2);
    __m128 c82_3 = _mm_load_ss(&C[(l_n*88)+33]);
    __m128 a82_3 = _mm_load_ss(&A[326]);
    c82_3 = _mm_add_ss(c82_3, _mm_mul_ss(a82_3, b82));
    _mm_store_ss(&C[(l_n*88)+33], c82_3);
    __m128 c82_4 = _mm_load_ss(&C[(l_n*88)+54]);
    __m128 a82_4 = _mm_load_ss(&A[327]);
    c82_4 = _mm_add_ss(c82_4, _mm_mul_ss(a82_4, b82));
    _mm_store_ss(&C[(l_n*88)+54], c82_4);
    __m128 c82_5 = _mm_load_ss(&C[(l_n*88)+82]);
    __m128 a82_5 = _mm_load_ss(&A[328]);
    c82_5 = _mm_add_ss(c82_5, _mm_mul_ss(a82_5, b82));
    _mm_store_ss(&C[(l_n*88)+82], c82_5);
#else
    C[(l_n*88)+2] += A[323] * B[(l_n*88)+82];
    C[(l_n*88)+8] += A[324] * B[(l_n*88)+82];
    C[(l_n*88)+18] += A[325] * B[(l_n*88)+82];
    C[(l_n*88)+33] += A[326] * B[(l_n*88)+82];
    C[(l_n*88)+54] += A[327] * B[(l_n*88)+82];
    C[(l_n*88)+82] += A[328] * B[(l_n*88)+82];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b83 = _mm_broadcast_ss(&B[(l_n*88)+83]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b83 = _mm_load_ss(&B[(l_n*88)+83]);    b83 = _mm_shuffle_ps(b83, b83, 0x00);
#endif
    __m128 c83_0 = _mm_load_ss(&C[(l_n*88)+0]);
    __m128 a83_0 = _mm_load_ss(&A[329]);
    c83_0 = _mm_add_ss(c83_0, _mm_mul_ss(a83_0, b83));
    _mm_store_ss(&C[(l_n*88)+0], c83_0);
    __m128 c83_1 = _mm_load_ss(&C[(l_n*88)+3]);
    __m128 a83_1 = _mm_load_ss(&A[330]);
    c83_1 = _mm_add_ss(c83_1, _mm_mul_ss(a83_1, b83));
    _mm_store_ss(&C[(l_n*88)+3], c83_1);
    __m128 c83_2 = _mm_load_ss(&C[(l_n*88)+9]);
    __m128 a83_2 = _mm_load_ss(&A[331]);
    c83_2 = _mm_add_ss(c83_2, _mm_mul_ss(a83_2, b83));
    _mm_store_ss(&C[(l_n*88)+9], c83_2);
    __m128 c83_3 = _mm_load_ss(&C[(l_n*88)+19]);
    __m128 a83_3 = _mm_load_ss(&A[332]);
    c83_3 = _mm_add_ss(c83_3, _mm_mul_ss(a83_3, b83));
    _mm_store_ss(&C[(l_n*88)+19], c83_3);
    __m128 c83_4 = _mm_load_ss(&C[(l_n*88)+34]);
    __m128 a83_4 = _mm_load_ss(&A[333]);
    c83_4 = _mm_add_ss(c83_4, _mm_mul_ss(a83_4, b83));
    _mm_store_ss(&C[(l_n*88)+34], c83_4);
    __m128 c83_5 = _mm_load_ss(&C[(l_n*88)+55]);
    __m128 a83_5 = _mm_load_ss(&A[334]);
    c83_5 = _mm_add_ss(c83_5, _mm_mul_ss(a83_5, b83));
    _mm_store_ss(&C[(l_n*88)+55], c83_5);
    __m128 c83_6 = _mm_load_ss(&C[(l_n*88)+83]);
    __m128 a83_6 = _mm_load_ss(&A[335]);
    c83_6 = _mm_add_ss(c83_6, _mm_mul_ss(a83_6, b83));
    _mm_store_ss(&C[(l_n*88)+83], c83_6);
#else
    C[(l_n*88)+0] += A[329] * B[(l_n*88)+83];
    C[(l_n*88)+3] += A[330] * B[(l_n*88)+83];
    C[(l_n*88)+9] += A[331] * B[(l_n*88)+83];
    C[(l_n*88)+19] += A[332] * B[(l_n*88)+83];
    C[(l_n*88)+34] += A[333] * B[(l_n*88)+83];
    C[(l_n*88)+55] += A[334] * B[(l_n*88)+83];
    C[(l_n*88)+83] += A[335] * B[(l_n*88)+83];
#endif

  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 6048;
#endif
}

void ssparse_fP113DivM_m84_n9_k84_ldAna7_ldB88_ldC88_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_n = 0;
  #pragma nounroll_and_jam
  for ( l_n = 0; l_n < 9; l_n++) {
    unsigned int l_m = 0;
   #pragma simd
    for ( l_m = 0; l_m < 84; l_m++) {
      C[(l_n*88)+l_m] = 0.0f;
    }
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b0 = _mm_broadcast_ss(&B[(l_n*88)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b0 = _mm_load_ss(&B[(l_n*88)+0]);    b0 = _mm_shuffle_ps(b0, b0, 0x00);
#endif
    __m128 c0_0 = _mm_load_ss(&C[(l_n*88)+0]);
    __m128 a0_0 = _mm_load_ss(&A[0]);
    c0_0 = _mm_add_ss(c0_0, _mm_mul_ss(a0_0, b0));
    _mm_store_ss(&C[(l_n*88)+0], c0_0);
    __m128 c0_1 = _mm_load_ss(&C[(l_n*88)+3]);
    __m128 a0_1 = _mm_load_ss(&A[1]);
    c0_1 = _mm_add_ss(c0_1, _mm_mul_ss(a0_1, b0));
    _mm_store_ss(&C[(l_n*88)+3], c0_1);
    __m128 c0_2 = _mm_load_ss(&C[(l_n*88)+9]);
    __m128 a0_2 = _mm_load_ss(&A[2]);
    c0_2 = _mm_add_ss(c0_2, _mm_mul_ss(a0_2, b0));
    _mm_store_ss(&C[(l_n*88)+9], c0_2);
    __m128 c0_3 = _mm_load_ss(&C[(l_n*88)+19]);
    __m128 a0_3 = _mm_load_ss(&A[3]);
    c0_3 = _mm_add_ss(c0_3, _mm_mul_ss(a0_3, b0));
    _mm_store_ss(&C[(l_n*88)+19], c0_3);
    __m128 c0_4 = _mm_load_ss(&C[(l_n*88)+34]);
    __m128 a0_4 = _mm_load_ss(&A[4]);
    c0_4 = _mm_add_ss(c0_4, _mm_mul_ss(a0_4, b0));
    _mm_store_ss(&C[(l_n*88)+34], c0_4);
    __m128 c0_5 = _mm_load_ss(&C[(l_n*88)+55]);
    __m128 a0_5 = _mm_load_ss(&A[5]);
    c0_5 = _mm_add_ss(c0_5, _mm_mul_ss(a0_5, b0));
    _mm_store_ss(&C[(l_n*88)+55], c0_5);
    __m128 c0_6 = _mm_load_ss(&C[(l_n*88)+83]);
    __m128 a0_6 = _mm_load_ss(&A[6]);
    c0_6 = _mm_add_ss(c0_6, _mm_mul_ss(a0_6, b0));
    _mm_store_ss(&C[(l_n*88)+83], c0_6);
#else
    C[(l_n*88)+0] += A[0] * B[(l_n*88)+0];
    C[(l_n*88)+3] += A[1] * B[(l_n*88)+0];
    C[(l_n*88)+9] += A[2] * B[(l_n*88)+0];
    C[(l_n*88)+19] += A[3] * B[(l_n*88)+0];
    C[(l_n*88)+34] += A[4] * B[(l_n*88)+0];
    C[(l_n*88)+55] += A[5] * B[(l_n*88)+0];
    C[(l_n*88)+83] += A[6] * B[(l_n*88)+0];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b1 = _mm_broadcast_ss(&B[(l_n*88)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b1 = _mm_load_ss(&B[(l_n*88)+1]);    b1 = _mm_shuffle_ps(b1, b1, 0x00);
#endif
    __m128 c1_0 = _mm_load_ss(&C[(l_n*88)+1]);
    __m128 a1_0 = _mm_load_ss(&A[7]);
    c1_0 = _mm_add_ss(c1_0, _mm_mul_ss(a1_0, b1));
    _mm_store_ss(&C[(l_n*88)+1], c1_0);
    __m128 c1_1 = _mm_load_ss(&C[(l_n*88)+7]);
    __m128 a1_1 = _mm_load_ss(&A[8]);
    c1_1 = _mm_add_ss(c1_1, _mm_mul_ss(a1_1, b1));
    _mm_store_ss(&C[(l_n*88)+7], c1_1);
    __m128 c1_2 = _mm_load_ss(&C[(l_n*88)+17]);
    __m128 a1_2 = _mm_load_ss(&A[9]);
    c1_2 = _mm_add_ss(c1_2, _mm_mul_ss(a1_2, b1));
    _mm_store_ss(&C[(l_n*88)+17], c1_2);
    __m128 c1_3 = _mm_load_ss(&C[(l_n*88)+32]);
    __m128 a1_3 = _mm_load_ss(&A[10]);
    c1_3 = _mm_add_ss(c1_3, _mm_mul_ss(a1_3, b1));
    _mm_store_ss(&C[(l_n*88)+32], c1_3);
    __m128 c1_4 = _mm_load_ss(&C[(l_n*88)+53]);
    __m128 a1_4 = _mm_load_ss(&A[11]);
    c1_4 = _mm_add_ss(c1_4, _mm_mul_ss(a1_4, b1));
    _mm_store_ss(&C[(l_n*88)+53], c1_4);
    __m128 c1_5 = _mm_load_ss(&C[(l_n*88)+81]);
    __m128 a1_5 = _mm_load_ss(&A[12]);
    c1_5 = _mm_add_ss(c1_5, _mm_mul_ss(a1_5, b1));
    _mm_store_ss(&C[(l_n*88)+81], c1_5);
#else
    C[(l_n*88)+1] += A[7] * B[(l_n*88)+1];
    C[(l_n*88)+7] += A[8] * B[(l_n*88)+1];
    C[(l_n*88)+17] += A[9] * B[(l_n*88)+1];
    C[(l_n*88)+32] += A[10] * B[(l_n*88)+1];
    C[(l_n*88)+53] += A[11] * B[(l_n*88)+1];
    C[(l_n*88)+81] += A[12] * B[(l_n*88)+1];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b2 = _mm_broadcast_ss(&B[(l_n*88)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b2 = _mm_load_ss(&B[(l_n*88)+2]);    b2 = _mm_shuffle_ps(b2, b2, 0x00);
#endif
    __m128 c2_0 = _mm_load_ss(&C[(l_n*88)+2]);
    __m128 a2_0 = _mm_load_ss(&A[13]);
    c2_0 = _mm_add_ss(c2_0, _mm_mul_ss(a2_0, b2));
    _mm_store_ss(&C[(l_n*88)+2], c2_0);
    __m128 c2_1 = _mm_load_ss(&C[(l_n*88)+8]);
    __m128 a2_1 = _mm_load_ss(&A[14]);
    c2_1 = _mm_add_ss(c2_1, _mm_mul_ss(a2_1, b2));
    _mm_store_ss(&C[(l_n*88)+8], c2_1);
    __m128 c2_2 = _mm_load_ss(&C[(l_n*88)+18]);
    __m128 a2_2 = _mm_load_ss(&A[15]);
    c2_2 = _mm_add_ss(c2_2, _mm_mul_ss(a2_2, b2));
    _mm_store_ss(&C[(l_n*88)+18], c2_2);
    __m128 c2_3 = _mm_load_ss(&C[(l_n*88)+33]);
    __m128 a2_3 = _mm_load_ss(&A[16]);
    c2_3 = _mm_add_ss(c2_3, _mm_mul_ss(a2_3, b2));
    _mm_store_ss(&C[(l_n*88)+33], c2_3);
    __m128 c2_4 = _mm_load_ss(&C[(l_n*88)+54]);
    __m128 a2_4 = _mm_load_ss(&A[17]);
    c2_4 = _mm_add_ss(c2_4, _mm_mul_ss(a2_4, b2));
    _mm_store_ss(&C[(l_n*88)+54], c2_4);
    __m128 c2_5 = _mm_load_ss(&C[(l_n*88)+82]);
    __m128 a2_5 = _mm_load_ss(&A[18]);
    c2_5 = _mm_add_ss(c2_5, _mm_mul_ss(a2_5, b2));
    _mm_store_ss(&C[(l_n*88)+82], c2_5);
#else
    C[(l_n*88)+2] += A[13] * B[(l_n*88)+2];
    C[(l_n*88)+8] += A[14] * B[(l_n*88)+2];
    C[(l_n*88)+18] += A[15] * B[(l_n*88)+2];
    C[(l_n*88)+33] += A[16] * B[(l_n*88)+2];
    C[(l_n*88)+54] += A[17] * B[(l_n*88)+2];
    C[(l_n*88)+82] += A[18] * B[(l_n*88)+2];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b3 = _mm_broadcast_ss(&B[(l_n*88)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b3 = _mm_load_ss(&B[(l_n*88)+3]);    b3 = _mm_shuffle_ps(b3, b3, 0x00);
#endif
    __m128 c3_0 = _mm_load_ss(&C[(l_n*88)+0]);
    __m128 a3_0 = _mm_load_ss(&A[19]);
    c3_0 = _mm_add_ss(c3_0, _mm_mul_ss(a3_0, b3));
    _mm_store_ss(&C[(l_n*88)+0], c3_0);
    __m128 c3_1 = _mm_load_ss(&C[(l_n*88)+3]);
    __m128 a3_1 = _mm_load_ss(&A[20]);
    c3_1 = _mm_add_ss(c3_1, _mm_mul_ss(a3_1, b3));
    _mm_store_ss(&C[(l_n*88)+3], c3_1);
    __m128 c3_2 = _mm_load_ss(&C[(l_n*88)+9]);
    __m128 a3_2 = _mm_load_ss(&A[21]);
    c3_2 = _mm_add_ss(c3_2, _mm_mul_ss(a3_2, b3));
    _mm_store_ss(&C[(l_n*88)+9], c3_2);
    __m128 c3_3 = _mm_load_ss(&C[(l_n*88)+19]);
    __m128 a3_3 = _mm_load_ss(&A[22]);
    c3_3 = _mm_add_ss(c3_3, _mm_mul_ss(a3_3, b3));
    _mm_store_ss(&C[(l_n*88)+19], c3_3);
    __m128 c3_4 = _mm_load_ss(&C[(l_n*88)+34]);
    __m128 a3_4 = _mm_load_ss(&A[23]);
    c3_4 = _mm_add_ss(c3_4, _mm_mul_ss(a3_4, b3));
    _mm_store_ss(&C[(l_n*88)+34], c3_4);
    __m128 c3_5 = _mm_load_ss(&C[(l_n*88)+55]);
    __m128 a3_5 = _mm_load_ss(&A[24]);
    c3_5 = _mm_add_ss(c3_5, _mm_mul_ss(a3_5, b3));
    _mm_store_ss(&C[(l_n*88)+55], c3_5);
    __m128 c3_6 = _mm_load_ss(&C[(l_n*88)+83]);
    __m128 a3_6 = _mm_load_ss(&A[25]);
    c3_6 = _mm_add_ss(c3_6, _mm_mul_ss(a3_6, b3));
    _mm_store_ss(&C[(l_n*88)+83], c3_6);
#else
    C[(l_n*88)+0] += A[19] * B[(l_n*88)+3];
    C[(l_n*88)+3] += A[20] * B[(l_n*88)+3];
    C[(l_n*88)+9] += A[21] * B[(l_n*88)+3];
    C[(l_n*88)+19] += A[22] * B[(l_n*88)+3];
    C[(l_n*88)+34] += A[23] * B[(l_n*88)+3];
    C[(l_n*88)+55] += A[24] * B[(l_n*88)+3];
    C[(l_n*88)+83] += A[25] * B[(l_n*88)+3];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b4 = _mm_broadcast_ss(&B[(l_n*88)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b4 = _mm_load_ss(&B[(l_n*88)+4]);    b4 = _mm_shuffle_ps(b4, b4, 0x00);
#endif
    __m128 c4_0 = _mm_load_ss(&C[(l_n*88)+4]);
    __m128 a4_0 = _mm_load_ss(&A[26]);
    c4_0 = _mm_add_ss(c4_0, _mm_mul_ss(a4_0, b4));
    _mm_store_ss(&C[(l_n*88)+4], c4_0);
    __m128 c4_1 = _mm_load_ss(&C[(l_n*88)+14]);
    __m128 a4_1 = _mm_load_ss(&A[27]);
    c4_1 = _mm_add_ss(c4_1, _mm_mul_ss(a4_1, b4));
    _mm_store_ss(&C[(l_n*88)+14], c4_1);
    __m128 c4_2 = _mm_load_ss(&C[(l_n*88)+29]);
    __m128 a4_2 = _mm_load_ss(&A[28]);
    c4_2 = _mm_add_ss(c4_2, _mm_mul_ss(a4_2, b4));
    _mm_store_ss(&C[(l_n*88)+29], c4_2);
    __m128 c4_3 = _mm_load_ss(&C[(l_n*88)+50]);
    __m128 a4_3 = _mm_load_ss(&A[29]);
    c4_3 = _mm_add_ss(c4_3, _mm_mul_ss(a4_3, b4));
    _mm_store_ss(&C[(l_n*88)+50], c4_3);
    __m128 c4_4 = _mm_load_ss(&C[(l_n*88)+78]);
    __m128 a4_4 = _mm_load_ss(&A[30]);
    c4_4 = _mm_add_ss(c4_4, _mm_mul_ss(a4_4, b4));
    _mm_store_ss(&C[(l_n*88)+78], c4_4);
#else
    C[(l_n*88)+4] += A[26] * B[(l_n*88)+4];
    C[(l_n*88)+14] += A[27] * B[(l_n*88)+4];
    C[(l_n*88)+29] += A[28] * B[(l_n*88)+4];
    C[(l_n*88)+50] += A[29] * B[(l_n*88)+4];
    C[(l_n*88)+78] += A[30] * B[(l_n*88)+4];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b5 = _mm_broadcast_ss(&B[(l_n*88)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b5 = _mm_load_ss(&B[(l_n*88)+5]);    b5 = _mm_shuffle_ps(b5, b5, 0x00);
#endif
    __m128 c5_0 = _mm_load_ss(&C[(l_n*88)+5]);
    __m128 a5_0 = _mm_load_ss(&A[31]);
    c5_0 = _mm_add_ss(c5_0, _mm_mul_ss(a5_0, b5));
    _mm_store_ss(&C[(l_n*88)+5], c5_0);
    __m128 c5_1 = _mm_load_ss(&C[(l_n*88)+15]);
    __m128 a5_1 = _mm_load_ss(&A[32]);
    c5_1 = _mm_add_ss(c5_1, _mm_mul_ss(a5_1, b5));
    _mm_store_ss(&C[(l_n*88)+15], c5_1);
    __m128 c5_2 = _mm_load_ss(&C[(l_n*88)+30]);
    __m128 a5_2 = _mm_load_ss(&A[33]);
    c5_2 = _mm_add_ss(c5_2, _mm_mul_ss(a5_2, b5));
    _mm_store_ss(&C[(l_n*88)+30], c5_2);
    __m128 c5_3 = _mm_load_ss(&C[(l_n*88)+51]);
    __m128 a5_3 = _mm_load_ss(&A[34]);
    c5_3 = _mm_add_ss(c5_3, _mm_mul_ss(a5_3, b5));
    _mm_store_ss(&C[(l_n*88)+51], c5_3);
    __m128 c5_4 = _mm_load_ss(&C[(l_n*88)+79]);
    __m128 a5_4 = _mm_load_ss(&A[35]);
    c5_4 = _mm_add_ss(c5_4, _mm_mul_ss(a5_4, b5));
    _mm_store_ss(&C[(l_n*88)+79], c5_4);
#else
    C[(l_n*88)+5] += A[31] * B[(l_n*88)+5];
    C[(l_n*88)+15] += A[32] * B[(l_n*88)+5];
    C[(l_n*88)+30] += A[33] * B[(l_n*88)+5];
    C[(l_n*88)+51] += A[34] * B[(l_n*88)+5];
    C[(l_n*88)+79] += A[35] * B[(l_n*88)+5];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b6 = _mm_broadcast_ss(&B[(l_n*88)+6]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b6 = _mm_load_ss(&B[(l_n*88)+6]);    b6 = _mm_shuffle_ps(b6, b6, 0x00);
#endif
    __m128 c6_0 = _mm_load_ss(&C[(l_n*88)+6]);
    __m128 a6_0 = _mm_load_ss(&A[36]);
    c6_0 = _mm_add_ss(c6_0, _mm_mul_ss(a6_0, b6));
    _mm_store_ss(&C[(l_n*88)+6], c6_0);
    __m128 c6_1 = _mm_load_ss(&C[(l_n*88)+16]);
    __m128 a6_1 = _mm_load_ss(&A[37]);
    c6_1 = _mm_add_ss(c6_1, _mm_mul_ss(a6_1, b6));
    _mm_store_ss(&C[(l_n*88)+16], c6_1);
    __m128 c6_2 = _mm_load_ss(&C[(l_n*88)+31]);
    __m128 a6_2 = _mm_load_ss(&A[38]);
    c6_2 = _mm_add_ss(c6_2, _mm_mul_ss(a6_2, b6));
    _mm_store_ss(&C[(l_n*88)+31], c6_2);
    __m128 c6_3 = _mm_load_ss(&C[(l_n*88)+52]);
    __m128 a6_3 = _mm_load_ss(&A[39]);
    c6_3 = _mm_add_ss(c6_3, _mm_mul_ss(a6_3, b6));
    _mm_store_ss(&C[(l_n*88)+52], c6_3);
    __m128 c6_4 = _mm_load_ss(&C[(l_n*88)+80]);
    __m128 a6_4 = _mm_load_ss(&A[40]);
    c6_4 = _mm_add_ss(c6_4, _mm_mul_ss(a6_4, b6));
    _mm_store_ss(&C[(l_n*88)+80], c6_4);
#else
    C[(l_n*88)+6] += A[36] * B[(l_n*88)+6];
    C[(l_n*88)+16] += A[37] * B[(l_n*88)+6];
    C[(l_n*88)+31] += A[38] * B[(l_n*88)+6];
    C[(l_n*88)+52] += A[39] * B[(l_n*88)+6];
    C[(l_n*88)+80] += A[40] * B[(l_n*88)+6];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b7 = _mm_broadcast_ss(&B[(l_n*88)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b7 = _mm_load_ss(&B[(l_n*88)+7]);    b7 = _mm_shuffle_ps(b7, b7, 0x00);
#endif
    __m128 c7_0 = _mm_load_ss(&C[(l_n*88)+1]);
    __m128 a7_0 = _mm_load_ss(&A[41]);
    c7_0 = _mm_add_ss(c7_0, _mm_mul_ss(a7_0, b7));
    _mm_store_ss(&C[(l_n*88)+1], c7_0);
    __m128 c7_1 = _mm_load_ss(&C[(l_n*88)+7]);
    __m128 a7_1 = _mm_load_ss(&A[42]);
    c7_1 = _mm_add_ss(c7_1, _mm_mul_ss(a7_1, b7));
    _mm_store_ss(&C[(l_n*88)+7], c7_1);
    __m128 c7_2 = _mm_load_ss(&C[(l_n*88)+17]);
    __m128 a7_2 = _mm_load_ss(&A[43]);
    c7_2 = _mm_add_ss(c7_2, _mm_mul_ss(a7_2, b7));
    _mm_store_ss(&C[(l_n*88)+17], c7_2);
    __m128 c7_3 = _mm_load_ss(&C[(l_n*88)+32]);
    __m128 a7_3 = _mm_load_ss(&A[44]);
    c7_3 = _mm_add_ss(c7_3, _mm_mul_ss(a7_3, b7));
    _mm_store_ss(&C[(l_n*88)+32], c7_3);
    __m128 c7_4 = _mm_load_ss(&C[(l_n*88)+53]);
    __m128 a7_4 = _mm_load_ss(&A[45]);
    c7_4 = _mm_add_ss(c7_4, _mm_mul_ss(a7_4, b7));
    _mm_store_ss(&C[(l_n*88)+53], c7_4);
    __m128 c7_5 = _mm_load_ss(&C[(l_n*88)+81]);
    __m128 a7_5 = _mm_load_ss(&A[46]);
    c7_5 = _mm_add_ss(c7_5, _mm_mul_ss(a7_5, b7));
    _mm_store_ss(&C[(l_n*88)+81], c7_5);
#else
    C[(l_n*88)+1] += A[41] * B[(l_n*88)+7];
    C[(l_n*88)+7] += A[42] * B[(l_n*88)+7];
    C[(l_n*88)+17] += A[43] * B[(l_n*88)+7];
    C[(l_n*88)+32] += A[44] * B[(l_n*88)+7];
    C[(l_n*88)+53] += A[45] * B[(l_n*88)+7];
    C[(l_n*88)+81] += A[46] * B[(l_n*88)+7];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b8 = _mm_broadcast_ss(&B[(l_n*88)+8]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b8 = _mm_load_ss(&B[(l_n*88)+8]);    b8 = _mm_shuffle_ps(b8, b8, 0x00);
#endif
    __m128 c8_0 = _mm_load_ss(&C[(l_n*88)+2]);
    __m128 a8_0 = _mm_load_ss(&A[47]);
    c8_0 = _mm_add_ss(c8_0, _mm_mul_ss(a8_0, b8));
    _mm_store_ss(&C[(l_n*88)+2], c8_0);
    __m128 c8_1 = _mm_load_ss(&C[(l_n*88)+8]);
    __m128 a8_1 = _mm_load_ss(&A[48]);
    c8_1 = _mm_add_ss(c8_1, _mm_mul_ss(a8_1, b8));
    _mm_store_ss(&C[(l_n*88)+8], c8_1);
    __m128 c8_2 = _mm_load_ss(&C[(l_n*88)+18]);
    __m128 a8_2 = _mm_load_ss(&A[49]);
    c8_2 = _mm_add_ss(c8_2, _mm_mul_ss(a8_2, b8));
    _mm_store_ss(&C[(l_n*88)+18], c8_2);
    __m128 c8_3 = _mm_load_ss(&C[(l_n*88)+33]);
    __m128 a8_3 = _mm_load_ss(&A[50]);
    c8_3 = _mm_add_ss(c8_3, _mm_mul_ss(a8_3, b8));
    _mm_store_ss(&C[(l_n*88)+33], c8_3);
    __m128 c8_4 = _mm_load_ss(&C[(l_n*88)+54]);
    __m128 a8_4 = _mm_load_ss(&A[51]);
    c8_4 = _mm_add_ss(c8_4, _mm_mul_ss(a8_4, b8));
    _mm_store_ss(&C[(l_n*88)+54], c8_4);
    __m128 c8_5 = _mm_load_ss(&C[(l_n*88)+82]);
    __m128 a8_5 = _mm_load_ss(&A[52]);
    c8_5 = _mm_add_ss(c8_5, _mm_mul_ss(a8_5, b8));
    _mm_store_ss(&C[(l_n*88)+82], c8_5);
#else
    C[(l_n*88)+2] += A[47] * B[(l_n*88)+8];
    C[(l_n*88)+8] += A[48] * B[(l_n*88)+8];
    C[(l_n*88)+18] += A[49] * B[(l_n*88)+8];
    C[(l_n*88)+33] += A[50] * B[(l_n*88)+8];
    C[(l_n*88)+54] += A[51] * B[(l_n*88)+8];
    C[(l_n*88)+82] += A[52] * B[(l_n*88)+8];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b9 = _mm_broadcast_ss(&B[(l_n*88)+9]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b9 = _mm_load_ss(&B[(l_n*88)+9]);    b9 = _mm_shuffle_ps(b9, b9, 0x00);
#endif
    __m128 c9_0 = _mm_load_ss(&C[(l_n*88)+0]);
    __m128 a9_0 = _mm_load_ss(&A[53]);
    c9_0 = _mm_add_ss(c9_0, _mm_mul_ss(a9_0, b9));
    _mm_store_ss(&C[(l_n*88)+0], c9_0);
    __m128 c9_1 = _mm_load_ss(&C[(l_n*88)+3]);
    __m128 a9_1 = _mm_load_ss(&A[54]);
    c9_1 = _mm_add_ss(c9_1, _mm_mul_ss(a9_1, b9));
    _mm_store_ss(&C[(l_n*88)+3], c9_1);
    __m128 c9_2 = _mm_load_ss(&C[(l_n*88)+9]);
    __m128 a9_2 = _mm_load_ss(&A[55]);
    c9_2 = _mm_add_ss(c9_2, _mm_mul_ss(a9_2, b9));
    _mm_store_ss(&C[(l_n*88)+9], c9_2);
    __m128 c9_3 = _mm_load_ss(&C[(l_n*88)+19]);
    __m128 a9_3 = _mm_load_ss(&A[56]);
    c9_3 = _mm_add_ss(c9_3, _mm_mul_ss(a9_3, b9));
    _mm_store_ss(&C[(l_n*88)+19], c9_3);
    __m128 c9_4 = _mm_load_ss(&C[(l_n*88)+34]);
    __m128 a9_4 = _mm_load_ss(&A[57]);
    c9_4 = _mm_add_ss(c9_4, _mm_mul_ss(a9_4, b9));
    _mm_store_ss(&C[(l_n*88)+34], c9_4);
    __m128 c9_5 = _mm_load_ss(&C[(l_n*88)+55]);
    __m128 a9_5 = _mm_load_ss(&A[58]);
    c9_5 = _mm_add_ss(c9_5, _mm_mul_ss(a9_5, b9));
    _mm_store_ss(&C[(l_n*88)+55], c9_5);
    __m128 c9_6 = _mm_load_ss(&C[(l_n*88)+83]);
    __m128 a9_6 = _mm_load_ss(&A[59]);
    c9_6 = _mm_add_ss(c9_6, _mm_mul_ss(a9_6, b9));
    _mm_store_ss(&C[(l_n*88)+83], c9_6);
#else
    C[(l_n*88)+0] += A[53] * B[(l_n*88)+9];
    C[(l_n*88)+3] += A[54] * B[(l_n*88)+9];
    C[(l_n*88)+9] += A[55] * B[(l_n*88)+9];
    C[(l_n*88)+19] += A[56] * B[(l_n*88)+9];
    C[(l_n*88)+34] += A[57] * B[(l_n*88)+9];
    C[(l_n*88)+55] += A[58] * B[(l_n*88)+9];
    C[(l_n*88)+83] += A[59] * B[(l_n*88)+9];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b10 = _mm_broadcast_ss(&B[(l_n*88)+10]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b10 = _mm_load_ss(&B[(l_n*88)+10]);    b10 = _mm_shuffle_ps(b10, b10, 0x00);
#endif
    __m128 c10_0 = _mm_load_ss(&C[(l_n*88)+10]);
    __m128 a10_0 = _mm_load_ss(&A[60]);
    c10_0 = _mm_add_ss(c10_0, _mm_mul_ss(a10_0, b10));
    _mm_store_ss(&C[(l_n*88)+10], c10_0);
    __m128 c10_1 = _mm_load_ss(&C[(l_n*88)+25]);
    __m128 a10_1 = _mm_load_ss(&A[61]);
    c10_1 = _mm_add_ss(c10_1, _mm_mul_ss(a10_1, b10));
    _mm_store_ss(&C[(l_n*88)+25], c10_1);
    __m128 c10_2 = _mm_load_ss(&C[(l_n*88)+46]);
    __m128 a10_2 = _mm_load_ss(&A[62]);
    c10_2 = _mm_add_ss(c10_2, _mm_mul_ss(a10_2, b10));
    _mm_store_ss(&C[(l_n*88)+46], c10_2);
    __m128 c10_3 = _mm_load_ss(&C[(l_n*88)+74]);
    __m128 a10_3 = _mm_load_ss(&A[63]);
    c10_3 = _mm_add_ss(c10_3, _mm_mul_ss(a10_3, b10));
    _mm_store_ss(&C[(l_n*88)+74], c10_3);
#else
    C[(l_n*88)+10] += A[60] * B[(l_n*88)+10];
    C[(l_n*88)+25] += A[61] * B[(l_n*88)+10];
    C[(l_n*88)+46] += A[62] * B[(l_n*88)+10];
    C[(l_n*88)+74] += A[63] * B[(l_n*88)+10];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b11 = _mm_broadcast_ss(&B[(l_n*88)+11]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b11 = _mm_load_ss(&B[(l_n*88)+11]);    b11 = _mm_shuffle_ps(b11, b11, 0x00);
#endif
    __m128 c11_0 = _mm_load_ss(&C[(l_n*88)+11]);
    __m128 a11_0 = _mm_load_ss(&A[64]);
    c11_0 = _mm_add_ss(c11_0, _mm_mul_ss(a11_0, b11));
    _mm_store_ss(&C[(l_n*88)+11], c11_0);
    __m128 c11_1 = _mm_load_ss(&C[(l_n*88)+26]);
    __m128 a11_1 = _mm_load_ss(&A[65]);
    c11_1 = _mm_add_ss(c11_1, _mm_mul_ss(a11_1, b11));
    _mm_store_ss(&C[(l_n*88)+26], c11_1);
    __m128 c11_2 = _mm_load_ss(&C[(l_n*88)+47]);
    __m128 a11_2 = _mm_load_ss(&A[66]);
    c11_2 = _mm_add_ss(c11_2, _mm_mul_ss(a11_2, b11));
    _mm_store_ss(&C[(l_n*88)+47], c11_2);
    __m128 c11_3 = _mm_load_ss(&C[(l_n*88)+75]);
    __m128 a11_3 = _mm_load_ss(&A[67]);
    c11_3 = _mm_add_ss(c11_3, _mm_mul_ss(a11_3, b11));
    _mm_store_ss(&C[(l_n*88)+75], c11_3);
#else
    C[(l_n*88)+11] += A[64] * B[(l_n*88)+11];
    C[(l_n*88)+26] += A[65] * B[(l_n*88)+11];
    C[(l_n*88)+47] += A[66] * B[(l_n*88)+11];
    C[(l_n*88)+75] += A[67] * B[(l_n*88)+11];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b12 = _mm_broadcast_ss(&B[(l_n*88)+12]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b12 = _mm_load_ss(&B[(l_n*88)+12]);    b12 = _mm_shuffle_ps(b12, b12, 0x00);
#endif
    __m128 c12_0 = _mm_load_ss(&C[(l_n*88)+12]);
    __m128 a12_0 = _mm_load_ss(&A[68]);
    c12_0 = _mm_add_ss(c12_0, _mm_mul_ss(a12_0, b12));
    _mm_store_ss(&C[(l_n*88)+12], c12_0);
    __m128 c12_1 = _mm_load_ss(&C[(l_n*88)+27]);
    __m128 a12_1 = _mm_load_ss(&A[69]);
    c12_1 = _mm_add_ss(c12_1, _mm_mul_ss(a12_1, b12));
    _mm_store_ss(&C[(l_n*88)+27], c12_1);
    __m128 c12_2 = _mm_load_ss(&C[(l_n*88)+48]);
    __m128 a12_2 = _mm_load_ss(&A[70]);
    c12_2 = _mm_add_ss(c12_2, _mm_mul_ss(a12_2, b12));
    _mm_store_ss(&C[(l_n*88)+48], c12_2);
    __m128 c12_3 = _mm_load_ss(&C[(l_n*88)+76]);
    __m128 a12_3 = _mm_load_ss(&A[71]);
    c12_3 = _mm_add_ss(c12_3, _mm_mul_ss(a12_3, b12));
    _mm_store_ss(&C[(l_n*88)+76], c12_3);
#else
    C[(l_n*88)+12] += A[68] * B[(l_n*88)+12];
    C[(l_n*88)+27] += A[69] * B[(l_n*88)+12];
    C[(l_n*88)+48] += A[70] * B[(l_n*88)+12];
    C[(l_n*88)+76] += A[71] * B[(l_n*88)+12];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b13 = _mm_broadcast_ss(&B[(l_n*88)+13]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b13 = _mm_load_ss(&B[(l_n*88)+13]);    b13 = _mm_shuffle_ps(b13, b13, 0x00);
#endif
    __m128 c13_0 = _mm_load_ss(&C[(l_n*88)+13]);
    __m128 a13_0 = _mm_load_ss(&A[72]);
    c13_0 = _mm_add_ss(c13_0, _mm_mul_ss(a13_0, b13));
    _mm_store_ss(&C[(l_n*88)+13], c13_0);
    __m128 c13_1 = _mm_load_ss(&C[(l_n*88)+28]);
    __m128 a13_1 = _mm_load_ss(&A[73]);
    c13_1 = _mm_add_ss(c13_1, _mm_mul_ss(a13_1, b13));
    _mm_store_ss(&C[(l_n*88)+28], c13_1);
    __m128 c13_2 = _mm_load_ss(&C[(l_n*88)+49]);
    __m128 a13_2 = _mm_load_ss(&A[74]);
    c13_2 = _mm_add_ss(c13_2, _mm_mul_ss(a13_2, b13));
    _mm_store_ss(&C[(l_n*88)+49], c13_2);
    __m128 c13_3 = _mm_load_ss(&C[(l_n*88)+77]);
    __m128 a13_3 = _mm_load_ss(&A[75]);
    c13_3 = _mm_add_ss(c13_3, _mm_mul_ss(a13_3, b13));
    _mm_store_ss(&C[(l_n*88)+77], c13_3);
#else
    C[(l_n*88)+13] += A[72] * B[(l_n*88)+13];
    C[(l_n*88)+28] += A[73] * B[(l_n*88)+13];
    C[(l_n*88)+49] += A[74] * B[(l_n*88)+13];
    C[(l_n*88)+77] += A[75] * B[(l_n*88)+13];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b14 = _mm_broadcast_ss(&B[(l_n*88)+14]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b14 = _mm_load_ss(&B[(l_n*88)+14]);    b14 = _mm_shuffle_ps(b14, b14, 0x00);
#endif
    __m128 c14_0 = _mm_load_ss(&C[(l_n*88)+4]);
    __m128 a14_0 = _mm_load_ss(&A[76]);
    c14_0 = _mm_add_ss(c14_0, _mm_mul_ss(a14_0, b14));
    _mm_store_ss(&C[(l_n*88)+4], c14_0);
    __m128 c14_1 = _mm_load_ss(&C[(l_n*88)+14]);
    __m128 a14_1 = _mm_load_ss(&A[77]);
    c14_1 = _mm_add_ss(c14_1, _mm_mul_ss(a14_1, b14));
    _mm_store_ss(&C[(l_n*88)+14], c14_1);
    __m128 c14_2 = _mm_load_ss(&C[(l_n*88)+29]);
    __m128 a14_2 = _mm_load_ss(&A[78]);
    c14_2 = _mm_add_ss(c14_2, _mm_mul_ss(a14_2, b14));
    _mm_store_ss(&C[(l_n*88)+29], c14_2);
    __m128 c14_3 = _mm_load_ss(&C[(l_n*88)+50]);
    __m128 a14_3 = _mm_load_ss(&A[79]);
    c14_3 = _mm_add_ss(c14_3, _mm_mul_ss(a14_3, b14));
    _mm_store_ss(&C[(l_n*88)+50], c14_3);
    __m128 c14_4 = _mm_load_ss(&C[(l_n*88)+78]);
    __m128 a14_4 = _mm_load_ss(&A[80]);
    c14_4 = _mm_add_ss(c14_4, _mm_mul_ss(a14_4, b14));
    _mm_store_ss(&C[(l_n*88)+78], c14_4);
#else
    C[(l_n*88)+4] += A[76] * B[(l_n*88)+14];
    C[(l_n*88)+14] += A[77] * B[(l_n*88)+14];
    C[(l_n*88)+29] += A[78] * B[(l_n*88)+14];
    C[(l_n*88)+50] += A[79] * B[(l_n*88)+14];
    C[(l_n*88)+78] += A[80] * B[(l_n*88)+14];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b15 = _mm_broadcast_ss(&B[(l_n*88)+15]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b15 = _mm_load_ss(&B[(l_n*88)+15]);    b15 = _mm_shuffle_ps(b15, b15, 0x00);
#endif
    __m128 c15_0 = _mm_load_ss(&C[(l_n*88)+5]);
    __m128 a15_0 = _mm_load_ss(&A[81]);
    c15_0 = _mm_add_ss(c15_0, _mm_mul_ss(a15_0, b15));
    _mm_store_ss(&C[(l_n*88)+5], c15_0);
    __m128 c15_1 = _mm_load_ss(&C[(l_n*88)+15]);
    __m128 a15_1 = _mm_load_ss(&A[82]);
    c15_1 = _mm_add_ss(c15_1, _mm_mul_ss(a15_1, b15));
    _mm_store_ss(&C[(l_n*88)+15], c15_1);
    __m128 c15_2 = _mm_load_ss(&C[(l_n*88)+30]);
    __m128 a15_2 = _mm_load_ss(&A[83]);
    c15_2 = _mm_add_ss(c15_2, _mm_mul_ss(a15_2, b15));
    _mm_store_ss(&C[(l_n*88)+30], c15_2);
    __m128 c15_3 = _mm_load_ss(&C[(l_n*88)+51]);
    __m128 a15_3 = _mm_load_ss(&A[84]);
    c15_3 = _mm_add_ss(c15_3, _mm_mul_ss(a15_3, b15));
    _mm_store_ss(&C[(l_n*88)+51], c15_3);
    __m128 c15_4 = _mm_load_ss(&C[(l_n*88)+79]);
    __m128 a15_4 = _mm_load_ss(&A[85]);
    c15_4 = _mm_add_ss(c15_4, _mm_mul_ss(a15_4, b15));
    _mm_store_ss(&C[(l_n*88)+79], c15_4);
#else
    C[(l_n*88)+5] += A[81] * B[(l_n*88)+15];
    C[(l_n*88)+15] += A[82] * B[(l_n*88)+15];
    C[(l_n*88)+30] += A[83] * B[(l_n*88)+15];
    C[(l_n*88)+51] += A[84] * B[(l_n*88)+15];
    C[(l_n*88)+79] += A[85] * B[(l_n*88)+15];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b16 = _mm_broadcast_ss(&B[(l_n*88)+16]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b16 = _mm_load_ss(&B[(l_n*88)+16]);    b16 = _mm_shuffle_ps(b16, b16, 0x00);
#endif
    __m128 c16_0 = _mm_load_ss(&C[(l_n*88)+6]);
    __m128 a16_0 = _mm_load_ss(&A[86]);
    c16_0 = _mm_add_ss(c16_0, _mm_mul_ss(a16_0, b16));
    _mm_store_ss(&C[(l_n*88)+6], c16_0);
    __m128 c16_1 = _mm_load_ss(&C[(l_n*88)+16]);
    __m128 a16_1 = _mm_load_ss(&A[87]);
    c16_1 = _mm_add_ss(c16_1, _mm_mul_ss(a16_1, b16));
    _mm_store_ss(&C[(l_n*88)+16], c16_1);
    __m128 c16_2 = _mm_load_ss(&C[(l_n*88)+31]);
    __m128 a16_2 = _mm_load_ss(&A[88]);
    c16_2 = _mm_add_ss(c16_2, _mm_mul_ss(a16_2, b16));
    _mm_store_ss(&C[(l_n*88)+31], c16_2);
    __m128 c16_3 = _mm_load_ss(&C[(l_n*88)+52]);
    __m128 a16_3 = _mm_load_ss(&A[89]);
    c16_3 = _mm_add_ss(c16_3, _mm_mul_ss(a16_3, b16));
    _mm_store_ss(&C[(l_n*88)+52], c16_3);
    __m128 c16_4 = _mm_load_ss(&C[(l_n*88)+80]);
    __m128 a16_4 = _mm_load_ss(&A[90]);
    c16_4 = _mm_add_ss(c16_4, _mm_mul_ss(a16_4, b16));
    _mm_store_ss(&C[(l_n*88)+80], c16_4);
#else
    C[(l_n*88)+6] += A[86] * B[(l_n*88)+16];
    C[(l_n*88)+16] += A[87] * B[(l_n*88)+16];
    C[(l_n*88)+31] += A[88] * B[(l_n*88)+16];
    C[(l_n*88)+52] += A[89] * B[(l_n*88)+16];
    C[(l_n*88)+80] += A[90] * B[(l_n*88)+16];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b17 = _mm_broadcast_ss(&B[(l_n*88)+17]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b17 = _mm_load_ss(&B[(l_n*88)+17]);    b17 = _mm_shuffle_ps(b17, b17, 0x00);
#endif
    __m128 c17_0 = _mm_load_ss(&C[(l_n*88)+1]);
    __m128 a17_0 = _mm_load_ss(&A[91]);
    c17_0 = _mm_add_ss(c17_0, _mm_mul_ss(a17_0, b17));
    _mm_store_ss(&C[(l_n*88)+1], c17_0);
    __m128 c17_1 = _mm_load_ss(&C[(l_n*88)+7]);
    __m128 a17_1 = _mm_load_ss(&A[92]);
    c17_1 = _mm_add_ss(c17_1, _mm_mul_ss(a17_1, b17));
    _mm_store_ss(&C[(l_n*88)+7], c17_1);
    __m128 c17_2 = _mm_load_ss(&C[(l_n*88)+17]);
    __m128 a17_2 = _mm_load_ss(&A[93]);
    c17_2 = _mm_add_ss(c17_2, _mm_mul_ss(a17_2, b17));
    _mm_store_ss(&C[(l_n*88)+17], c17_2);
    __m128 c17_3 = _mm_load_ss(&C[(l_n*88)+32]);
    __m128 a17_3 = _mm_load_ss(&A[94]);
    c17_3 = _mm_add_ss(c17_3, _mm_mul_ss(a17_3, b17));
    _mm_store_ss(&C[(l_n*88)+32], c17_3);
    __m128 c17_4 = _mm_load_ss(&C[(l_n*88)+53]);
    __m128 a17_4 = _mm_load_ss(&A[95]);
    c17_4 = _mm_add_ss(c17_4, _mm_mul_ss(a17_4, b17));
    _mm_store_ss(&C[(l_n*88)+53], c17_4);
    __m128 c17_5 = _mm_load_ss(&C[(l_n*88)+81]);
    __m128 a17_5 = _mm_load_ss(&A[96]);
    c17_5 = _mm_add_ss(c17_5, _mm_mul_ss(a17_5, b17));
    _mm_store_ss(&C[(l_n*88)+81], c17_5);
#else
    C[(l_n*88)+1] += A[91] * B[(l_n*88)+17];
    C[(l_n*88)+7] += A[92] * B[(l_n*88)+17];
    C[(l_n*88)+17] += A[93] * B[(l_n*88)+17];
    C[(l_n*88)+32] += A[94] * B[(l_n*88)+17];
    C[(l_n*88)+53] += A[95] * B[(l_n*88)+17];
    C[(l_n*88)+81] += A[96] * B[(l_n*88)+17];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b18 = _mm_broadcast_ss(&B[(l_n*88)+18]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b18 = _mm_load_ss(&B[(l_n*88)+18]);    b18 = _mm_shuffle_ps(b18, b18, 0x00);
#endif
    __m128 c18_0 = _mm_load_ss(&C[(l_n*88)+2]);
    __m128 a18_0 = _mm_load_ss(&A[97]);
    c18_0 = _mm_add_ss(c18_0, _mm_mul_ss(a18_0, b18));
    _mm_store_ss(&C[(l_n*88)+2], c18_0);
    __m128 c18_1 = _mm_load_ss(&C[(l_n*88)+8]);
    __m128 a18_1 = _mm_load_ss(&A[98]);
    c18_1 = _mm_add_ss(c18_1, _mm_mul_ss(a18_1, b18));
    _mm_store_ss(&C[(l_n*88)+8], c18_1);
    __m128 c18_2 = _mm_load_ss(&C[(l_n*88)+18]);
    __m128 a18_2 = _mm_load_ss(&A[99]);
    c18_2 = _mm_add_ss(c18_2, _mm_mul_ss(a18_2, b18));
    _mm_store_ss(&C[(l_n*88)+18], c18_2);
    __m128 c18_3 = _mm_load_ss(&C[(l_n*88)+33]);
    __m128 a18_3 = _mm_load_ss(&A[100]);
    c18_3 = _mm_add_ss(c18_3, _mm_mul_ss(a18_3, b18));
    _mm_store_ss(&C[(l_n*88)+33], c18_3);
    __m128 c18_4 = _mm_load_ss(&C[(l_n*88)+54]);
    __m128 a18_4 = _mm_load_ss(&A[101]);
    c18_4 = _mm_add_ss(c18_4, _mm_mul_ss(a18_4, b18));
    _mm_store_ss(&C[(l_n*88)+54], c18_4);
    __m128 c18_5 = _mm_load_ss(&C[(l_n*88)+82]);
    __m128 a18_5 = _mm_load_ss(&A[102]);
    c18_5 = _mm_add_ss(c18_5, _mm_mul_ss(a18_5, b18));
    _mm_store_ss(&C[(l_n*88)+82], c18_5);
#else
    C[(l_n*88)+2] += A[97] * B[(l_n*88)+18];
    C[(l_n*88)+8] += A[98] * B[(l_n*88)+18];
    C[(l_n*88)+18] += A[99] * B[(l_n*88)+18];
    C[(l_n*88)+33] += A[100] * B[(l_n*88)+18];
    C[(l_n*88)+54] += A[101] * B[(l_n*88)+18];
    C[(l_n*88)+82] += A[102] * B[(l_n*88)+18];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b19 = _mm_broadcast_ss(&B[(l_n*88)+19]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b19 = _mm_load_ss(&B[(l_n*88)+19]);    b19 = _mm_shuffle_ps(b19, b19, 0x00);
#endif
    __m128 c19_0 = _mm_load_ss(&C[(l_n*88)+0]);
    __m128 a19_0 = _mm_load_ss(&A[103]);
    c19_0 = _mm_add_ss(c19_0, _mm_mul_ss(a19_0, b19));
    _mm_store_ss(&C[(l_n*88)+0], c19_0);
    __m128 c19_1 = _mm_load_ss(&C[(l_n*88)+3]);
    __m128 a19_1 = _mm_load_ss(&A[104]);
    c19_1 = _mm_add_ss(c19_1, _mm_mul_ss(a19_1, b19));
    _mm_store_ss(&C[(l_n*88)+3], c19_1);
    __m128 c19_2 = _mm_load_ss(&C[(l_n*88)+9]);
    __m128 a19_2 = _mm_load_ss(&A[105]);
    c19_2 = _mm_add_ss(c19_2, _mm_mul_ss(a19_2, b19));
    _mm_store_ss(&C[(l_n*88)+9], c19_2);
    __m128 c19_3 = _mm_load_ss(&C[(l_n*88)+19]);
    __m128 a19_3 = _mm_load_ss(&A[106]);
    c19_3 = _mm_add_ss(c19_3, _mm_mul_ss(a19_3, b19));
    _mm_store_ss(&C[(l_n*88)+19], c19_3);
    __m128 c19_4 = _mm_load_ss(&C[(l_n*88)+34]);
    __m128 a19_4 = _mm_load_ss(&A[107]);
    c19_4 = _mm_add_ss(c19_4, _mm_mul_ss(a19_4, b19));
    _mm_store_ss(&C[(l_n*88)+34], c19_4);
    __m128 c19_5 = _mm_load_ss(&C[(l_n*88)+55]);
    __m128 a19_5 = _mm_load_ss(&A[108]);
    c19_5 = _mm_add_ss(c19_5, _mm_mul_ss(a19_5, b19));
    _mm_store_ss(&C[(l_n*88)+55], c19_5);
    __m128 c19_6 = _mm_load_ss(&C[(l_n*88)+83]);
    __m128 a19_6 = _mm_load_ss(&A[109]);
    c19_6 = _mm_add_ss(c19_6, _mm_mul_ss(a19_6, b19));
    _mm_store_ss(&C[(l_n*88)+83], c19_6);
#else
    C[(l_n*88)+0] += A[103] * B[(l_n*88)+19];
    C[(l_n*88)+3] += A[104] * B[(l_n*88)+19];
    C[(l_n*88)+9] += A[105] * B[(l_n*88)+19];
    C[(l_n*88)+19] += A[106] * B[(l_n*88)+19];
    C[(l_n*88)+34] += A[107] * B[(l_n*88)+19];
    C[(l_n*88)+55] += A[108] * B[(l_n*88)+19];
    C[(l_n*88)+83] += A[109] * B[(l_n*88)+19];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b20 = _mm_broadcast_ss(&B[(l_n*88)+20]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b20 = _mm_load_ss(&B[(l_n*88)+20]);    b20 = _mm_shuffle_ps(b20, b20, 0x00);
#endif
    __m128 c20_0 = _mm_load_ss(&C[(l_n*88)+20]);
    __m128 a20_0 = _mm_load_ss(&A[110]);
    c20_0 = _mm_add_ss(c20_0, _mm_mul_ss(a20_0, b20));
    _mm_store_ss(&C[(l_n*88)+20], c20_0);
    __m128 c20_1 = _mm_load_ss(&C[(l_n*88)+41]);
    __m128 a20_1 = _mm_load_ss(&A[111]);
    c20_1 = _mm_add_ss(c20_1, _mm_mul_ss(a20_1, b20));
    _mm_store_ss(&C[(l_n*88)+41], c20_1);
    __m128 c20_2 = _mm_load_ss(&C[(l_n*88)+69]);
    __m128 a20_2 = _mm_load_ss(&A[112]);
    c20_2 = _mm_add_ss(c20_2, _mm_mul_ss(a20_2, b20));
    _mm_store_ss(&C[(l_n*88)+69], c20_2);
#else
    C[(l_n*88)+20] += A[110] * B[(l_n*88)+20];
    C[(l_n*88)+41] += A[111] * B[(l_n*88)+20];
    C[(l_n*88)+69] += A[112] * B[(l_n*88)+20];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b21 = _mm_broadcast_ss(&B[(l_n*88)+21]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b21 = _mm_load_ss(&B[(l_n*88)+21]);    b21 = _mm_shuffle_ps(b21, b21, 0x00);
#endif
    __m128 c21_0 = _mm_load_ss(&C[(l_n*88)+21]);
    __m128 a21_0 = _mm_load_ss(&A[113]);
    c21_0 = _mm_add_ss(c21_0, _mm_mul_ss(a21_0, b21));
    _mm_store_ss(&C[(l_n*88)+21], c21_0);
    __m128 c21_1 = _mm_load_ss(&C[(l_n*88)+42]);
    __m128 a21_1 = _mm_load_ss(&A[114]);
    c21_1 = _mm_add_ss(c21_1, _mm_mul_ss(a21_1, b21));
    _mm_store_ss(&C[(l_n*88)+42], c21_1);
    __m128 c21_2 = _mm_load_ss(&C[(l_n*88)+70]);
    __m128 a21_2 = _mm_load_ss(&A[115]);
    c21_2 = _mm_add_ss(c21_2, _mm_mul_ss(a21_2, b21));
    _mm_store_ss(&C[(l_n*88)+70], c21_2);
#else
    C[(l_n*88)+21] += A[113] * B[(l_n*88)+21];
    C[(l_n*88)+42] += A[114] * B[(l_n*88)+21];
    C[(l_n*88)+70] += A[115] * B[(l_n*88)+21];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b22 = _mm_broadcast_ss(&B[(l_n*88)+22]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b22 = _mm_load_ss(&B[(l_n*88)+22]);    b22 = _mm_shuffle_ps(b22, b22, 0x00);
#endif
    __m128 c22_0 = _mm_load_ss(&C[(l_n*88)+22]);
    __m128 a22_0 = _mm_load_ss(&A[116]);
    c22_0 = _mm_add_ss(c22_0, _mm_mul_ss(a22_0, b22));
    _mm_store_ss(&C[(l_n*88)+22], c22_0);
    __m128 c22_1 = _mm_load_ss(&C[(l_n*88)+43]);
    __m128 a22_1 = _mm_load_ss(&A[117]);
    c22_1 = _mm_add_ss(c22_1, _mm_mul_ss(a22_1, b22));
    _mm_store_ss(&C[(l_n*88)+43], c22_1);
    __m128 c22_2 = _mm_load_ss(&C[(l_n*88)+71]);
    __m128 a22_2 = _mm_load_ss(&A[118]);
    c22_2 = _mm_add_ss(c22_2, _mm_mul_ss(a22_2, b22));
    _mm_store_ss(&C[(l_n*88)+71], c22_2);
#else
    C[(l_n*88)+22] += A[116] * B[(l_n*88)+22];
    C[(l_n*88)+43] += A[117] * B[(l_n*88)+22];
    C[(l_n*88)+71] += A[118] * B[(l_n*88)+22];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b23 = _mm_broadcast_ss(&B[(l_n*88)+23]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b23 = _mm_load_ss(&B[(l_n*88)+23]);    b23 = _mm_shuffle_ps(b23, b23, 0x00);
#endif
    __m128 c23_0 = _mm_load_ss(&C[(l_n*88)+23]);
    __m128 a23_0 = _mm_load_ss(&A[119]);
    c23_0 = _mm_add_ss(c23_0, _mm_mul_ss(a23_0, b23));
    _mm_store_ss(&C[(l_n*88)+23], c23_0);
    __m128 c23_1 = _mm_load_ss(&C[(l_n*88)+44]);
    __m128 a23_1 = _mm_load_ss(&A[120]);
    c23_1 = _mm_add_ss(c23_1, _mm_mul_ss(a23_1, b23));
    _mm_store_ss(&C[(l_n*88)+44], c23_1);
    __m128 c23_2 = _mm_load_ss(&C[(l_n*88)+72]);
    __m128 a23_2 = _mm_load_ss(&A[121]);
    c23_2 = _mm_add_ss(c23_2, _mm_mul_ss(a23_2, b23));
    _mm_store_ss(&C[(l_n*88)+72], c23_2);
#else
    C[(l_n*88)+23] += A[119] * B[(l_n*88)+23];
    C[(l_n*88)+44] += A[120] * B[(l_n*88)+23];
    C[(l_n*88)+72] += A[121] * B[(l_n*88)+23];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b24 = _mm_broadcast_ss(&B[(l_n*88)+24]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b24 = _mm_load_ss(&B[(l_n*88)+24]);    b24 = _mm_shuffle_ps(b24, b24, 0x00);
#endif
    __m128 c24_0 = _mm_load_ss(&C[(l_n*88)+24]);
    __m128 a24_0 = _mm_load_ss(&A[122]);
    c24_0 = _mm_add_ss(c24_0, _mm_mul_ss(a24_0, b24));
    _mm_store_ss(&C[(l_n*88)+24], c24_0);
    __m128 c24_1 = _mm_load_ss(&C[(l_n*88)+45]);
    __m128 a24_1 = _mm_load_ss(&A[123]);
    c24_1 = _mm_add_ss(c24_1, _mm_mul_ss(a24_1, b24));
    _mm_store_ss(&C[(l_n*88)+45], c24_1);
    __m128 c24_2 = _mm_load_ss(&C[(l_n*88)+73]);
    __m128 a24_2 = _mm_load_ss(&A[124]);
    c24_2 = _mm_add_ss(c24_2, _mm_mul_ss(a24_2, b24));
    _mm_store_ss(&C[(l_n*88)+73], c24_2);
#else
    C[(l_n*88)+24] += A[122] * B[(l_n*88)+24];
    C[(l_n*88)+45] += A[123] * B[(l_n*88)+24];
    C[(l_n*88)+73] += A[124] * B[(l_n*88)+24];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b25 = _mm_broadcast_ss(&B[(l_n*88)+25]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b25 = _mm_load_ss(&B[(l_n*88)+25]);    b25 = _mm_shuffle_ps(b25, b25, 0x00);
#endif
    __m128 c25_0 = _mm_load_ss(&C[(l_n*88)+10]);
    __m128 a25_0 = _mm_load_ss(&A[125]);
    c25_0 = _mm_add_ss(c25_0, _mm_mul_ss(a25_0, b25));
    _mm_store_ss(&C[(l_n*88)+10], c25_0);
    __m128 c25_1 = _mm_load_ss(&C[(l_n*88)+25]);
    __m128 a25_1 = _mm_load_ss(&A[126]);
    c25_1 = _mm_add_ss(c25_1, _mm_mul_ss(a25_1, b25));
    _mm_store_ss(&C[(l_n*88)+25], c25_1);
    __m128 c25_2 = _mm_load_ss(&C[(l_n*88)+46]);
    __m128 a25_2 = _mm_load_ss(&A[127]);
    c25_2 = _mm_add_ss(c25_2, _mm_mul_ss(a25_2, b25));
    _mm_store_ss(&C[(l_n*88)+46], c25_2);
    __m128 c25_3 = _mm_load_ss(&C[(l_n*88)+74]);
    __m128 a25_3 = _mm_load_ss(&A[128]);
    c25_3 = _mm_add_ss(c25_3, _mm_mul_ss(a25_3, b25));
    _mm_store_ss(&C[(l_n*88)+74], c25_3);
#else
    C[(l_n*88)+10] += A[125] * B[(l_n*88)+25];
    C[(l_n*88)+25] += A[126] * B[(l_n*88)+25];
    C[(l_n*88)+46] += A[127] * B[(l_n*88)+25];
    C[(l_n*88)+74] += A[128] * B[(l_n*88)+25];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b26 = _mm_broadcast_ss(&B[(l_n*88)+26]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b26 = _mm_load_ss(&B[(l_n*88)+26]);    b26 = _mm_shuffle_ps(b26, b26, 0x00);
#endif
    __m128 c26_0 = _mm_load_ss(&C[(l_n*88)+11]);
    __m128 a26_0 = _mm_load_ss(&A[129]);
    c26_0 = _mm_add_ss(c26_0, _mm_mul_ss(a26_0, b26));
    _mm_store_ss(&C[(l_n*88)+11], c26_0);
    __m128 c26_1 = _mm_load_ss(&C[(l_n*88)+26]);
    __m128 a26_1 = _mm_load_ss(&A[130]);
    c26_1 = _mm_add_ss(c26_1, _mm_mul_ss(a26_1, b26));
    _mm_store_ss(&C[(l_n*88)+26], c26_1);
    __m128 c26_2 = _mm_load_ss(&C[(l_n*88)+47]);
    __m128 a26_2 = _mm_load_ss(&A[131]);
    c26_2 = _mm_add_ss(c26_2, _mm_mul_ss(a26_2, b26));
    _mm_store_ss(&C[(l_n*88)+47], c26_2);
    __m128 c26_3 = _mm_load_ss(&C[(l_n*88)+75]);
    __m128 a26_3 = _mm_load_ss(&A[132]);
    c26_3 = _mm_add_ss(c26_3, _mm_mul_ss(a26_3, b26));
    _mm_store_ss(&C[(l_n*88)+75], c26_3);
#else
    C[(l_n*88)+11] += A[129] * B[(l_n*88)+26];
    C[(l_n*88)+26] += A[130] * B[(l_n*88)+26];
    C[(l_n*88)+47] += A[131] * B[(l_n*88)+26];
    C[(l_n*88)+75] += A[132] * B[(l_n*88)+26];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b27 = _mm_broadcast_ss(&B[(l_n*88)+27]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b27 = _mm_load_ss(&B[(l_n*88)+27]);    b27 = _mm_shuffle_ps(b27, b27, 0x00);
#endif
    __m128 c27_0 = _mm_load_ss(&C[(l_n*88)+12]);
    __m128 a27_0 = _mm_load_ss(&A[133]);
    c27_0 = _mm_add_ss(c27_0, _mm_mul_ss(a27_0, b27));
    _mm_store_ss(&C[(l_n*88)+12], c27_0);
    __m128 c27_1 = _mm_load_ss(&C[(l_n*88)+27]);
    __m128 a27_1 = _mm_load_ss(&A[134]);
    c27_1 = _mm_add_ss(c27_1, _mm_mul_ss(a27_1, b27));
    _mm_store_ss(&C[(l_n*88)+27], c27_1);
    __m128 c27_2 = _mm_load_ss(&C[(l_n*88)+48]);
    __m128 a27_2 = _mm_load_ss(&A[135]);
    c27_2 = _mm_add_ss(c27_2, _mm_mul_ss(a27_2, b27));
    _mm_store_ss(&C[(l_n*88)+48], c27_2);
    __m128 c27_3 = _mm_load_ss(&C[(l_n*88)+76]);
    __m128 a27_3 = _mm_load_ss(&A[136]);
    c27_3 = _mm_add_ss(c27_3, _mm_mul_ss(a27_3, b27));
    _mm_store_ss(&C[(l_n*88)+76], c27_3);
#else
    C[(l_n*88)+12] += A[133] * B[(l_n*88)+27];
    C[(l_n*88)+27] += A[134] * B[(l_n*88)+27];
    C[(l_n*88)+48] += A[135] * B[(l_n*88)+27];
    C[(l_n*88)+76] += A[136] * B[(l_n*88)+27];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b28 = _mm_broadcast_ss(&B[(l_n*88)+28]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b28 = _mm_load_ss(&B[(l_n*88)+28]);    b28 = _mm_shuffle_ps(b28, b28, 0x00);
#endif
    __m128 c28_0 = _mm_load_ss(&C[(l_n*88)+13]);
    __m128 a28_0 = _mm_load_ss(&A[137]);
    c28_0 = _mm_add_ss(c28_0, _mm_mul_ss(a28_0, b28));
    _mm_store_ss(&C[(l_n*88)+13], c28_0);
    __m128 c28_1 = _mm_load_ss(&C[(l_n*88)+28]);
    __m128 a28_1 = _mm_load_ss(&A[138]);
    c28_1 = _mm_add_ss(c28_1, _mm_mul_ss(a28_1, b28));
    _mm_store_ss(&C[(l_n*88)+28], c28_1);
    __m128 c28_2 = _mm_load_ss(&C[(l_n*88)+49]);
    __m128 a28_2 = _mm_load_ss(&A[139]);
    c28_2 = _mm_add_ss(c28_2, _mm_mul_ss(a28_2, b28));
    _mm_store_ss(&C[(l_n*88)+49], c28_2);
    __m128 c28_3 = _mm_load_ss(&C[(l_n*88)+77]);
    __m128 a28_3 = _mm_load_ss(&A[140]);
    c28_3 = _mm_add_ss(c28_3, _mm_mul_ss(a28_3, b28));
    _mm_store_ss(&C[(l_n*88)+77], c28_3);
#else
    C[(l_n*88)+13] += A[137] * B[(l_n*88)+28];
    C[(l_n*88)+28] += A[138] * B[(l_n*88)+28];
    C[(l_n*88)+49] += A[139] * B[(l_n*88)+28];
    C[(l_n*88)+77] += A[140] * B[(l_n*88)+28];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b29 = _mm_broadcast_ss(&B[(l_n*88)+29]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b29 = _mm_load_ss(&B[(l_n*88)+29]);    b29 = _mm_shuffle_ps(b29, b29, 0x00);
#endif
    __m128 c29_0 = _mm_load_ss(&C[(l_n*88)+4]);
    __m128 a29_0 = _mm_load_ss(&A[141]);
    c29_0 = _mm_add_ss(c29_0, _mm_mul_ss(a29_0, b29));
    _mm_store_ss(&C[(l_n*88)+4], c29_0);
    __m128 c29_1 = _mm_load_ss(&C[(l_n*88)+14]);
    __m128 a29_1 = _mm_load_ss(&A[142]);
    c29_1 = _mm_add_ss(c29_1, _mm_mul_ss(a29_1, b29));
    _mm_store_ss(&C[(l_n*88)+14], c29_1);
    __m128 c29_2 = _mm_load_ss(&C[(l_n*88)+29]);
    __m128 a29_2 = _mm_load_ss(&A[143]);
    c29_2 = _mm_add_ss(c29_2, _mm_mul_ss(a29_2, b29));
    _mm_store_ss(&C[(l_n*88)+29], c29_2);
    __m128 c29_3 = _mm_load_ss(&C[(l_n*88)+50]);
    __m128 a29_3 = _mm_load_ss(&A[144]);
    c29_3 = _mm_add_ss(c29_3, _mm_mul_ss(a29_3, b29));
    _mm_store_ss(&C[(l_n*88)+50], c29_3);
    __m128 c29_4 = _mm_load_ss(&C[(l_n*88)+78]);
    __m128 a29_4 = _mm_load_ss(&A[145]);
    c29_4 = _mm_add_ss(c29_4, _mm_mul_ss(a29_4, b29));
    _mm_store_ss(&C[(l_n*88)+78], c29_4);
#else
    C[(l_n*88)+4] += A[141] * B[(l_n*88)+29];
    C[(l_n*88)+14] += A[142] * B[(l_n*88)+29];
    C[(l_n*88)+29] += A[143] * B[(l_n*88)+29];
    C[(l_n*88)+50] += A[144] * B[(l_n*88)+29];
    C[(l_n*88)+78] += A[145] * B[(l_n*88)+29];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b30 = _mm_broadcast_ss(&B[(l_n*88)+30]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b30 = _mm_load_ss(&B[(l_n*88)+30]);    b30 = _mm_shuffle_ps(b30, b30, 0x00);
#endif
    __m128 c30_0 = _mm_load_ss(&C[(l_n*88)+5]);
    __m128 a30_0 = _mm_load_ss(&A[146]);
    c30_0 = _mm_add_ss(c30_0, _mm_mul_ss(a30_0, b30));
    _mm_store_ss(&C[(l_n*88)+5], c30_0);
    __m128 c30_1 = _mm_load_ss(&C[(l_n*88)+15]);
    __m128 a30_1 = _mm_load_ss(&A[147]);
    c30_1 = _mm_add_ss(c30_1, _mm_mul_ss(a30_1, b30));
    _mm_store_ss(&C[(l_n*88)+15], c30_1);
    __m128 c30_2 = _mm_load_ss(&C[(l_n*88)+30]);
    __m128 a30_2 = _mm_load_ss(&A[148]);
    c30_2 = _mm_add_ss(c30_2, _mm_mul_ss(a30_2, b30));
    _mm_store_ss(&C[(l_n*88)+30], c30_2);
    __m128 c30_3 = _mm_load_ss(&C[(l_n*88)+51]);
    __m128 a30_3 = _mm_load_ss(&A[149]);
    c30_3 = _mm_add_ss(c30_3, _mm_mul_ss(a30_3, b30));
    _mm_store_ss(&C[(l_n*88)+51], c30_3);
    __m128 c30_4 = _mm_load_ss(&C[(l_n*88)+79]);
    __m128 a30_4 = _mm_load_ss(&A[150]);
    c30_4 = _mm_add_ss(c30_4, _mm_mul_ss(a30_4, b30));
    _mm_store_ss(&C[(l_n*88)+79], c30_4);
#else
    C[(l_n*88)+5] += A[146] * B[(l_n*88)+30];
    C[(l_n*88)+15] += A[147] * B[(l_n*88)+30];
    C[(l_n*88)+30] += A[148] * B[(l_n*88)+30];
    C[(l_n*88)+51] += A[149] * B[(l_n*88)+30];
    C[(l_n*88)+79] += A[150] * B[(l_n*88)+30];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b31 = _mm_broadcast_ss(&B[(l_n*88)+31]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b31 = _mm_load_ss(&B[(l_n*88)+31]);    b31 = _mm_shuffle_ps(b31, b31, 0x00);
#endif
    __m128 c31_0 = _mm_load_ss(&C[(l_n*88)+6]);
    __m128 a31_0 = _mm_load_ss(&A[151]);
    c31_0 = _mm_add_ss(c31_0, _mm_mul_ss(a31_0, b31));
    _mm_store_ss(&C[(l_n*88)+6], c31_0);
    __m128 c31_1 = _mm_load_ss(&C[(l_n*88)+16]);
    __m128 a31_1 = _mm_load_ss(&A[152]);
    c31_1 = _mm_add_ss(c31_1, _mm_mul_ss(a31_1, b31));
    _mm_store_ss(&C[(l_n*88)+16], c31_1);
    __m128 c31_2 = _mm_load_ss(&C[(l_n*88)+31]);
    __m128 a31_2 = _mm_load_ss(&A[153]);
    c31_2 = _mm_add_ss(c31_2, _mm_mul_ss(a31_2, b31));
    _mm_store_ss(&C[(l_n*88)+31], c31_2);
    __m128 c31_3 = _mm_load_ss(&C[(l_n*88)+52]);
    __m128 a31_3 = _mm_load_ss(&A[154]);
    c31_3 = _mm_add_ss(c31_3, _mm_mul_ss(a31_3, b31));
    _mm_store_ss(&C[(l_n*88)+52], c31_3);
    __m128 c31_4 = _mm_load_ss(&C[(l_n*88)+80]);
    __m128 a31_4 = _mm_load_ss(&A[155]);
    c31_4 = _mm_add_ss(c31_4, _mm_mul_ss(a31_4, b31));
    _mm_store_ss(&C[(l_n*88)+80], c31_4);
#else
    C[(l_n*88)+6] += A[151] * B[(l_n*88)+31];
    C[(l_n*88)+16] += A[152] * B[(l_n*88)+31];
    C[(l_n*88)+31] += A[153] * B[(l_n*88)+31];
    C[(l_n*88)+52] += A[154] * B[(l_n*88)+31];
    C[(l_n*88)+80] += A[155] * B[(l_n*88)+31];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b32 = _mm_broadcast_ss(&B[(l_n*88)+32]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b32 = _mm_load_ss(&B[(l_n*88)+32]);    b32 = _mm_shuffle_ps(b32, b32, 0x00);
#endif
    __m128 c32_0 = _mm_load_ss(&C[(l_n*88)+1]);
    __m128 a32_0 = _mm_load_ss(&A[156]);
    c32_0 = _mm_add_ss(c32_0, _mm_mul_ss(a32_0, b32));
    _mm_store_ss(&C[(l_n*88)+1], c32_0);
    __m128 c32_1 = _mm_load_ss(&C[(l_n*88)+7]);
    __m128 a32_1 = _mm_load_ss(&A[157]);
    c32_1 = _mm_add_ss(c32_1, _mm_mul_ss(a32_1, b32));
    _mm_store_ss(&C[(l_n*88)+7], c32_1);
    __m128 c32_2 = _mm_load_ss(&C[(l_n*88)+17]);
    __m128 a32_2 = _mm_load_ss(&A[158]);
    c32_2 = _mm_add_ss(c32_2, _mm_mul_ss(a32_2, b32));
    _mm_store_ss(&C[(l_n*88)+17], c32_2);
    __m128 c32_3 = _mm_load_ss(&C[(l_n*88)+32]);
    __m128 a32_3 = _mm_load_ss(&A[159]);
    c32_3 = _mm_add_ss(c32_3, _mm_mul_ss(a32_3, b32));
    _mm_store_ss(&C[(l_n*88)+32], c32_3);
    __m128 c32_4 = _mm_load_ss(&C[(l_n*88)+53]);
    __m128 a32_4 = _mm_load_ss(&A[160]);
    c32_4 = _mm_add_ss(c32_4, _mm_mul_ss(a32_4, b32));
    _mm_store_ss(&C[(l_n*88)+53], c32_4);
    __m128 c32_5 = _mm_load_ss(&C[(l_n*88)+81]);
    __m128 a32_5 = _mm_load_ss(&A[161]);
    c32_5 = _mm_add_ss(c32_5, _mm_mul_ss(a32_5, b32));
    _mm_store_ss(&C[(l_n*88)+81], c32_5);
#else
    C[(l_n*88)+1] += A[156] * B[(l_n*88)+32];
    C[(l_n*88)+7] += A[157] * B[(l_n*88)+32];
    C[(l_n*88)+17] += A[158] * B[(l_n*88)+32];
    C[(l_n*88)+32] += A[159] * B[(l_n*88)+32];
    C[(l_n*88)+53] += A[160] * B[(l_n*88)+32];
    C[(l_n*88)+81] += A[161] * B[(l_n*88)+32];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b33 = _mm_broadcast_ss(&B[(l_n*88)+33]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b33 = _mm_load_ss(&B[(l_n*88)+33]);    b33 = _mm_shuffle_ps(b33, b33, 0x00);
#endif
    __m128 c33_0 = _mm_load_ss(&C[(l_n*88)+2]);
    __m128 a33_0 = _mm_load_ss(&A[162]);
    c33_0 = _mm_add_ss(c33_0, _mm_mul_ss(a33_0, b33));
    _mm_store_ss(&C[(l_n*88)+2], c33_0);
    __m128 c33_1 = _mm_load_ss(&C[(l_n*88)+8]);
    __m128 a33_1 = _mm_load_ss(&A[163]);
    c33_1 = _mm_add_ss(c33_1, _mm_mul_ss(a33_1, b33));
    _mm_store_ss(&C[(l_n*88)+8], c33_1);
    __m128 c33_2 = _mm_load_ss(&C[(l_n*88)+18]);
    __m128 a33_2 = _mm_load_ss(&A[164]);
    c33_2 = _mm_add_ss(c33_2, _mm_mul_ss(a33_2, b33));
    _mm_store_ss(&C[(l_n*88)+18], c33_2);
    __m128 c33_3 = _mm_load_ss(&C[(l_n*88)+33]);
    __m128 a33_3 = _mm_load_ss(&A[165]);
    c33_3 = _mm_add_ss(c33_3, _mm_mul_ss(a33_3, b33));
    _mm_store_ss(&C[(l_n*88)+33], c33_3);
    __m128 c33_4 = _mm_load_ss(&C[(l_n*88)+54]);
    __m128 a33_4 = _mm_load_ss(&A[166]);
    c33_4 = _mm_add_ss(c33_4, _mm_mul_ss(a33_4, b33));
    _mm_store_ss(&C[(l_n*88)+54], c33_4);
    __m128 c33_5 = _mm_load_ss(&C[(l_n*88)+82]);
    __m128 a33_5 = _mm_load_ss(&A[167]);
    c33_5 = _mm_add_ss(c33_5, _mm_mul_ss(a33_5, b33));
    _mm_store_ss(&C[(l_n*88)+82], c33_5);
#else
    C[(l_n*88)+2] += A[162] * B[(l_n*88)+33];
    C[(l_n*88)+8] += A[163] * B[(l_n*88)+33];
    C[(l_n*88)+18] += A[164] * B[(l_n*88)+33];
    C[(l_n*88)+33] += A[165] * B[(l_n*88)+33];
    C[(l_n*88)+54] += A[166] * B[(l_n*88)+33];
    C[(l_n*88)+82] += A[167] * B[(l_n*88)+33];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b34 = _mm_broadcast_ss(&B[(l_n*88)+34]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b34 = _mm_load_ss(&B[(l_n*88)+34]);    b34 = _mm_shuffle_ps(b34, b34, 0x00);
#endif
    __m128 c34_0 = _mm_load_ss(&C[(l_n*88)+0]);
    __m128 a34_0 = _mm_load_ss(&A[168]);
    c34_0 = _mm_add_ss(c34_0, _mm_mul_ss(a34_0, b34));
    _mm_store_ss(&C[(l_n*88)+0], c34_0);
    __m128 c34_1 = _mm_load_ss(&C[(l_n*88)+3]);
    __m128 a34_1 = _mm_load_ss(&A[169]);
    c34_1 = _mm_add_ss(c34_1, _mm_mul_ss(a34_1, b34));
    _mm_store_ss(&C[(l_n*88)+3], c34_1);
    __m128 c34_2 = _mm_load_ss(&C[(l_n*88)+9]);
    __m128 a34_2 = _mm_load_ss(&A[170]);
    c34_2 = _mm_add_ss(c34_2, _mm_mul_ss(a34_2, b34));
    _mm_store_ss(&C[(l_n*88)+9], c34_2);
    __m128 c34_3 = _mm_load_ss(&C[(l_n*88)+19]);
    __m128 a34_3 = _mm_load_ss(&A[171]);
    c34_3 = _mm_add_ss(c34_3, _mm_mul_ss(a34_3, b34));
    _mm_store_ss(&C[(l_n*88)+19], c34_3);
    __m128 c34_4 = _mm_load_ss(&C[(l_n*88)+34]);
    __m128 a34_4 = _mm_load_ss(&A[172]);
    c34_4 = _mm_add_ss(c34_4, _mm_mul_ss(a34_4, b34));
    _mm_store_ss(&C[(l_n*88)+34], c34_4);
    __m128 c34_5 = _mm_load_ss(&C[(l_n*88)+55]);
    __m128 a34_5 = _mm_load_ss(&A[173]);
    c34_5 = _mm_add_ss(c34_5, _mm_mul_ss(a34_5, b34));
    _mm_store_ss(&C[(l_n*88)+55], c34_5);
    __m128 c34_6 = _mm_load_ss(&C[(l_n*88)+83]);
    __m128 a34_6 = _mm_load_ss(&A[174]);
    c34_6 = _mm_add_ss(c34_6, _mm_mul_ss(a34_6, b34));
    _mm_store_ss(&C[(l_n*88)+83], c34_6);
#else
    C[(l_n*88)+0] += A[168] * B[(l_n*88)+34];
    C[(l_n*88)+3] += A[169] * B[(l_n*88)+34];
    C[(l_n*88)+9] += A[170] * B[(l_n*88)+34];
    C[(l_n*88)+19] += A[171] * B[(l_n*88)+34];
    C[(l_n*88)+34] += A[172] * B[(l_n*88)+34];
    C[(l_n*88)+55] += A[173] * B[(l_n*88)+34];
    C[(l_n*88)+83] += A[174] * B[(l_n*88)+34];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b35 = _mm_broadcast_ss(&B[(l_n*88)+35]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b35 = _mm_load_ss(&B[(l_n*88)+35]);    b35 = _mm_shuffle_ps(b35, b35, 0x00);
#endif
    __m128 c35_0 = _mm_load_ss(&C[(l_n*88)+35]);
    __m128 a35_0 = _mm_load_ss(&A[175]);
    c35_0 = _mm_add_ss(c35_0, _mm_mul_ss(a35_0, b35));
    _mm_store_ss(&C[(l_n*88)+35], c35_0);
    __m128 c35_1 = _mm_load_ss(&C[(l_n*88)+63]);
    __m128 a35_1 = _mm_load_ss(&A[176]);
    c35_1 = _mm_add_ss(c35_1, _mm_mul_ss(a35_1, b35));
    _mm_store_ss(&C[(l_n*88)+63], c35_1);
#else
    C[(l_n*88)+35] += A[175] * B[(l_n*88)+35];
    C[(l_n*88)+63] += A[176] * B[(l_n*88)+35];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b36 = _mm_broadcast_ss(&B[(l_n*88)+36]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b36 = _mm_load_ss(&B[(l_n*88)+36]);    b36 = _mm_shuffle_ps(b36, b36, 0x00);
#endif
    __m128 c36_0 = _mm_load_ss(&C[(l_n*88)+36]);
    __m128 a36_0 = _mm_load_ss(&A[177]);
    c36_0 = _mm_add_ss(c36_0, _mm_mul_ss(a36_0, b36));
    _mm_store_ss(&C[(l_n*88)+36], c36_0);
    __m128 c36_1 = _mm_load_ss(&C[(l_n*88)+64]);
    __m128 a36_1 = _mm_load_ss(&A[178]);
    c36_1 = _mm_add_ss(c36_1, _mm_mul_ss(a36_1, b36));
    _mm_store_ss(&C[(l_n*88)+64], c36_1);
#else
    C[(l_n*88)+36] += A[177] * B[(l_n*88)+36];
    C[(l_n*88)+64] += A[178] * B[(l_n*88)+36];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b37 = _mm_broadcast_ss(&B[(l_n*88)+37]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b37 = _mm_load_ss(&B[(l_n*88)+37]);    b37 = _mm_shuffle_ps(b37, b37, 0x00);
#endif
    __m128 c37_0 = _mm_load_ss(&C[(l_n*88)+37]);
    __m128 a37_0 = _mm_load_ss(&A[179]);
    c37_0 = _mm_add_ss(c37_0, _mm_mul_ss(a37_0, b37));
    _mm_store_ss(&C[(l_n*88)+37], c37_0);
    __m128 c37_1 = _mm_load_ss(&C[(l_n*88)+65]);
    __m128 a37_1 = _mm_load_ss(&A[180]);
    c37_1 = _mm_add_ss(c37_1, _mm_mul_ss(a37_1, b37));
    _mm_store_ss(&C[(l_n*88)+65], c37_1);
#else
    C[(l_n*88)+37] += A[179] * B[(l_n*88)+37];
    C[(l_n*88)+65] += A[180] * B[(l_n*88)+37];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b38 = _mm_broadcast_ss(&B[(l_n*88)+38]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b38 = _mm_load_ss(&B[(l_n*88)+38]);    b38 = _mm_shuffle_ps(b38, b38, 0x00);
#endif
    __m128 c38_0 = _mm_load_ss(&C[(l_n*88)+38]);
    __m128 a38_0 = _mm_load_ss(&A[181]);
    c38_0 = _mm_add_ss(c38_0, _mm_mul_ss(a38_0, b38));
    _mm_store_ss(&C[(l_n*88)+38], c38_0);
    __m128 c38_1 = _mm_load_ss(&C[(l_n*88)+66]);
    __m128 a38_1 = _mm_load_ss(&A[182]);
    c38_1 = _mm_add_ss(c38_1, _mm_mul_ss(a38_1, b38));
    _mm_store_ss(&C[(l_n*88)+66], c38_1);
#else
    C[(l_n*88)+38] += A[181] * B[(l_n*88)+38];
    C[(l_n*88)+66] += A[182] * B[(l_n*88)+38];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b39 = _mm_broadcast_ss(&B[(l_n*88)+39]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b39 = _mm_load_ss(&B[(l_n*88)+39]);    b39 = _mm_shuffle_ps(b39, b39, 0x00);
#endif
    __m128 c39_0 = _mm_load_ss(&C[(l_n*88)+39]);
    __m128 a39_0 = _mm_load_ss(&A[183]);
    c39_0 = _mm_add_ss(c39_0, _mm_mul_ss(a39_0, b39));
    _mm_store_ss(&C[(l_n*88)+39], c39_0);
    __m128 c39_1 = _mm_load_ss(&C[(l_n*88)+67]);
    __m128 a39_1 = _mm_load_ss(&A[184]);
    c39_1 = _mm_add_ss(c39_1, _mm_mul_ss(a39_1, b39));
    _mm_store_ss(&C[(l_n*88)+67], c39_1);
#else
    C[(l_n*88)+39] += A[183] * B[(l_n*88)+39];
    C[(l_n*88)+67] += A[184] * B[(l_n*88)+39];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b40 = _mm_broadcast_ss(&B[(l_n*88)+40]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b40 = _mm_load_ss(&B[(l_n*88)+40]);    b40 = _mm_shuffle_ps(b40, b40, 0x00);
#endif
    __m128 c40_0 = _mm_load_ss(&C[(l_n*88)+40]);
    __m128 a40_0 = _mm_load_ss(&A[185]);
    c40_0 = _mm_add_ss(c40_0, _mm_mul_ss(a40_0, b40));
    _mm_store_ss(&C[(l_n*88)+40], c40_0);
    __m128 c40_1 = _mm_load_ss(&C[(l_n*88)+68]);
    __m128 a40_1 = _mm_load_ss(&A[186]);
    c40_1 = _mm_add_ss(c40_1, _mm_mul_ss(a40_1, b40));
    _mm_store_ss(&C[(l_n*88)+68], c40_1);
#else
    C[(l_n*88)+40] += A[185] * B[(l_n*88)+40];
    C[(l_n*88)+68] += A[186] * B[(l_n*88)+40];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b41 = _mm_broadcast_ss(&B[(l_n*88)+41]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b41 = _mm_load_ss(&B[(l_n*88)+41]);    b41 = _mm_shuffle_ps(b41, b41, 0x00);
#endif
    __m128 c41_0 = _mm_load_ss(&C[(l_n*88)+20]);
    __m128 a41_0 = _mm_load_ss(&A[187]);
    c41_0 = _mm_add_ss(c41_0, _mm_mul_ss(a41_0, b41));
    _mm_store_ss(&C[(l_n*88)+20], c41_0);
    __m128 c41_1 = _mm_load_ss(&C[(l_n*88)+41]);
    __m128 a41_1 = _mm_load_ss(&A[188]);
    c41_1 = _mm_add_ss(c41_1, _mm_mul_ss(a41_1, b41));
    _mm_store_ss(&C[(l_n*88)+41], c41_1);
    __m128 c41_2 = _mm_load_ss(&C[(l_n*88)+69]);
    __m128 a41_2 = _mm_load_ss(&A[189]);
    c41_2 = _mm_add_ss(c41_2, _mm_mul_ss(a41_2, b41));
    _mm_store_ss(&C[(l_n*88)+69], c41_2);
#else
    C[(l_n*88)+20] += A[187] * B[(l_n*88)+41];
    C[(l_n*88)+41] += A[188] * B[(l_n*88)+41];
    C[(l_n*88)+69] += A[189] * B[(l_n*88)+41];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b42 = _mm_broadcast_ss(&B[(l_n*88)+42]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b42 = _mm_load_ss(&B[(l_n*88)+42]);    b42 = _mm_shuffle_ps(b42, b42, 0x00);
#endif
    __m128 c42_0 = _mm_load_ss(&C[(l_n*88)+21]);
    __m128 a42_0 = _mm_load_ss(&A[190]);
    c42_0 = _mm_add_ss(c42_0, _mm_mul_ss(a42_0, b42));
    _mm_store_ss(&C[(l_n*88)+21], c42_0);
    __m128 c42_1 = _mm_load_ss(&C[(l_n*88)+42]);
    __m128 a42_1 = _mm_load_ss(&A[191]);
    c42_1 = _mm_add_ss(c42_1, _mm_mul_ss(a42_1, b42));
    _mm_store_ss(&C[(l_n*88)+42], c42_1);
    __m128 c42_2 = _mm_load_ss(&C[(l_n*88)+70]);
    __m128 a42_2 = _mm_load_ss(&A[192]);
    c42_2 = _mm_add_ss(c42_2, _mm_mul_ss(a42_2, b42));
    _mm_store_ss(&C[(l_n*88)+70], c42_2);
#else
    C[(l_n*88)+21] += A[190] * B[(l_n*88)+42];
    C[(l_n*88)+42] += A[191] * B[(l_n*88)+42];
    C[(l_n*88)+70] += A[192] * B[(l_n*88)+42];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b43 = _mm_broadcast_ss(&B[(l_n*88)+43]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b43 = _mm_load_ss(&B[(l_n*88)+43]);    b43 = _mm_shuffle_ps(b43, b43, 0x00);
#endif
    __m128 c43_0 = _mm_load_ss(&C[(l_n*88)+22]);
    __m128 a43_0 = _mm_load_ss(&A[193]);
    c43_0 = _mm_add_ss(c43_0, _mm_mul_ss(a43_0, b43));
    _mm_store_ss(&C[(l_n*88)+22], c43_0);
    __m128 c43_1 = _mm_load_ss(&C[(l_n*88)+43]);
    __m128 a43_1 = _mm_load_ss(&A[194]);
    c43_1 = _mm_add_ss(c43_1, _mm_mul_ss(a43_1, b43));
    _mm_store_ss(&C[(l_n*88)+43], c43_1);
    __m128 c43_2 = _mm_load_ss(&C[(l_n*88)+71]);
    __m128 a43_2 = _mm_load_ss(&A[195]);
    c43_2 = _mm_add_ss(c43_2, _mm_mul_ss(a43_2, b43));
    _mm_store_ss(&C[(l_n*88)+71], c43_2);
#else
    C[(l_n*88)+22] += A[193] * B[(l_n*88)+43];
    C[(l_n*88)+43] += A[194] * B[(l_n*88)+43];
    C[(l_n*88)+71] += A[195] * B[(l_n*88)+43];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b44 = _mm_broadcast_ss(&B[(l_n*88)+44]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b44 = _mm_load_ss(&B[(l_n*88)+44]);    b44 = _mm_shuffle_ps(b44, b44, 0x00);
#endif
    __m128 c44_0 = _mm_load_ss(&C[(l_n*88)+23]);
    __m128 a44_0 = _mm_load_ss(&A[196]);
    c44_0 = _mm_add_ss(c44_0, _mm_mul_ss(a44_0, b44));
    _mm_store_ss(&C[(l_n*88)+23], c44_0);
    __m128 c44_1 = _mm_load_ss(&C[(l_n*88)+44]);
    __m128 a44_1 = _mm_load_ss(&A[197]);
    c44_1 = _mm_add_ss(c44_1, _mm_mul_ss(a44_1, b44));
    _mm_store_ss(&C[(l_n*88)+44], c44_1);
    __m128 c44_2 = _mm_load_ss(&C[(l_n*88)+72]);
    __m128 a44_2 = _mm_load_ss(&A[198]);
    c44_2 = _mm_add_ss(c44_2, _mm_mul_ss(a44_2, b44));
    _mm_store_ss(&C[(l_n*88)+72], c44_2);
#else
    C[(l_n*88)+23] += A[196] * B[(l_n*88)+44];
    C[(l_n*88)+44] += A[197] * B[(l_n*88)+44];
    C[(l_n*88)+72] += A[198] * B[(l_n*88)+44];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b45 = _mm_broadcast_ss(&B[(l_n*88)+45]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b45 = _mm_load_ss(&B[(l_n*88)+45]);    b45 = _mm_shuffle_ps(b45, b45, 0x00);
#endif
    __m128 c45_0 = _mm_load_ss(&C[(l_n*88)+24]);
    __m128 a45_0 = _mm_load_ss(&A[199]);
    c45_0 = _mm_add_ss(c45_0, _mm_mul_ss(a45_0, b45));
    _mm_store_ss(&C[(l_n*88)+24], c45_0);
    __m128 c45_1 = _mm_load_ss(&C[(l_n*88)+45]);
    __m128 a45_1 = _mm_load_ss(&A[200]);
    c45_1 = _mm_add_ss(c45_1, _mm_mul_ss(a45_1, b45));
    _mm_store_ss(&C[(l_n*88)+45], c45_1);
    __m128 c45_2 = _mm_load_ss(&C[(l_n*88)+73]);
    __m128 a45_2 = _mm_load_ss(&A[201]);
    c45_2 = _mm_add_ss(c45_2, _mm_mul_ss(a45_2, b45));
    _mm_store_ss(&C[(l_n*88)+73], c45_2);
#else
    C[(l_n*88)+24] += A[199] * B[(l_n*88)+45];
    C[(l_n*88)+45] += A[200] * B[(l_n*88)+45];
    C[(l_n*88)+73] += A[201] * B[(l_n*88)+45];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b46 = _mm_broadcast_ss(&B[(l_n*88)+46]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b46 = _mm_load_ss(&B[(l_n*88)+46]);    b46 = _mm_shuffle_ps(b46, b46, 0x00);
#endif
    __m128 c46_0 = _mm_load_ss(&C[(l_n*88)+10]);
    __m128 a46_0 = _mm_load_ss(&A[202]);
    c46_0 = _mm_add_ss(c46_0, _mm_mul_ss(a46_0, b46));
    _mm_store_ss(&C[(l_n*88)+10], c46_0);
    __m128 c46_1 = _mm_load_ss(&C[(l_n*88)+25]);
    __m128 a46_1 = _mm_load_ss(&A[203]);
    c46_1 = _mm_add_ss(c46_1, _mm_mul_ss(a46_1, b46));
    _mm_store_ss(&C[(l_n*88)+25], c46_1);
    __m128 c46_2 = _mm_load_ss(&C[(l_n*88)+46]);
    __m128 a46_2 = _mm_load_ss(&A[204]);
    c46_2 = _mm_add_ss(c46_2, _mm_mul_ss(a46_2, b46));
    _mm_store_ss(&C[(l_n*88)+46], c46_2);
    __m128 c46_3 = _mm_load_ss(&C[(l_n*88)+74]);
    __m128 a46_3 = _mm_load_ss(&A[205]);
    c46_3 = _mm_add_ss(c46_3, _mm_mul_ss(a46_3, b46));
    _mm_store_ss(&C[(l_n*88)+74], c46_3);
#else
    C[(l_n*88)+10] += A[202] * B[(l_n*88)+46];
    C[(l_n*88)+25] += A[203] * B[(l_n*88)+46];
    C[(l_n*88)+46] += A[204] * B[(l_n*88)+46];
    C[(l_n*88)+74] += A[205] * B[(l_n*88)+46];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b47 = _mm_broadcast_ss(&B[(l_n*88)+47]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b47 = _mm_load_ss(&B[(l_n*88)+47]);    b47 = _mm_shuffle_ps(b47, b47, 0x00);
#endif
    __m128 c47_0 = _mm_load_ss(&C[(l_n*88)+11]);
    __m128 a47_0 = _mm_load_ss(&A[206]);
    c47_0 = _mm_add_ss(c47_0, _mm_mul_ss(a47_0, b47));
    _mm_store_ss(&C[(l_n*88)+11], c47_0);
    __m128 c47_1 = _mm_load_ss(&C[(l_n*88)+26]);
    __m128 a47_1 = _mm_load_ss(&A[207]);
    c47_1 = _mm_add_ss(c47_1, _mm_mul_ss(a47_1, b47));
    _mm_store_ss(&C[(l_n*88)+26], c47_1);
    __m128 c47_2 = _mm_load_ss(&C[(l_n*88)+47]);
    __m128 a47_2 = _mm_load_ss(&A[208]);
    c47_2 = _mm_add_ss(c47_2, _mm_mul_ss(a47_2, b47));
    _mm_store_ss(&C[(l_n*88)+47], c47_2);
    __m128 c47_3 = _mm_load_ss(&C[(l_n*88)+75]);
    __m128 a47_3 = _mm_load_ss(&A[209]);
    c47_3 = _mm_add_ss(c47_3, _mm_mul_ss(a47_3, b47));
    _mm_store_ss(&C[(l_n*88)+75], c47_3);
#else
    C[(l_n*88)+11] += A[206] * B[(l_n*88)+47];
    C[(l_n*88)+26] += A[207] * B[(l_n*88)+47];
    C[(l_n*88)+47] += A[208] * B[(l_n*88)+47];
    C[(l_n*88)+75] += A[209] * B[(l_n*88)+47];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b48 = _mm_broadcast_ss(&B[(l_n*88)+48]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b48 = _mm_load_ss(&B[(l_n*88)+48]);    b48 = _mm_shuffle_ps(b48, b48, 0x00);
#endif
    __m128 c48_0 = _mm_load_ss(&C[(l_n*88)+12]);
    __m128 a48_0 = _mm_load_ss(&A[210]);
    c48_0 = _mm_add_ss(c48_0, _mm_mul_ss(a48_0, b48));
    _mm_store_ss(&C[(l_n*88)+12], c48_0);
    __m128 c48_1 = _mm_load_ss(&C[(l_n*88)+27]);
    __m128 a48_1 = _mm_load_ss(&A[211]);
    c48_1 = _mm_add_ss(c48_1, _mm_mul_ss(a48_1, b48));
    _mm_store_ss(&C[(l_n*88)+27], c48_1);
    __m128 c48_2 = _mm_load_ss(&C[(l_n*88)+48]);
    __m128 a48_2 = _mm_load_ss(&A[212]);
    c48_2 = _mm_add_ss(c48_2, _mm_mul_ss(a48_2, b48));
    _mm_store_ss(&C[(l_n*88)+48], c48_2);
    __m128 c48_3 = _mm_load_ss(&C[(l_n*88)+76]);
    __m128 a48_3 = _mm_load_ss(&A[213]);
    c48_3 = _mm_add_ss(c48_3, _mm_mul_ss(a48_3, b48));
    _mm_store_ss(&C[(l_n*88)+76], c48_3);
#else
    C[(l_n*88)+12] += A[210] * B[(l_n*88)+48];
    C[(l_n*88)+27] += A[211] * B[(l_n*88)+48];
    C[(l_n*88)+48] += A[212] * B[(l_n*88)+48];
    C[(l_n*88)+76] += A[213] * B[(l_n*88)+48];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b49 = _mm_broadcast_ss(&B[(l_n*88)+49]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b49 = _mm_load_ss(&B[(l_n*88)+49]);    b49 = _mm_shuffle_ps(b49, b49, 0x00);
#endif
    __m128 c49_0 = _mm_load_ss(&C[(l_n*88)+13]);
    __m128 a49_0 = _mm_load_ss(&A[214]);
    c49_0 = _mm_add_ss(c49_0, _mm_mul_ss(a49_0, b49));
    _mm_store_ss(&C[(l_n*88)+13], c49_0);
    __m128 c49_1 = _mm_load_ss(&C[(l_n*88)+28]);
    __m128 a49_1 = _mm_load_ss(&A[215]);
    c49_1 = _mm_add_ss(c49_1, _mm_mul_ss(a49_1, b49));
    _mm_store_ss(&C[(l_n*88)+28], c49_1);
    __m128 c49_2 = _mm_load_ss(&C[(l_n*88)+49]);
    __m128 a49_2 = _mm_load_ss(&A[216]);
    c49_2 = _mm_add_ss(c49_2, _mm_mul_ss(a49_2, b49));
    _mm_store_ss(&C[(l_n*88)+49], c49_2);
    __m128 c49_3 = _mm_load_ss(&C[(l_n*88)+77]);
    __m128 a49_3 = _mm_load_ss(&A[217]);
    c49_3 = _mm_add_ss(c49_3, _mm_mul_ss(a49_3, b49));
    _mm_store_ss(&C[(l_n*88)+77], c49_3);
#else
    C[(l_n*88)+13] += A[214] * B[(l_n*88)+49];
    C[(l_n*88)+28] += A[215] * B[(l_n*88)+49];
    C[(l_n*88)+49] += A[216] * B[(l_n*88)+49];
    C[(l_n*88)+77] += A[217] * B[(l_n*88)+49];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b50 = _mm_broadcast_ss(&B[(l_n*88)+50]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b50 = _mm_load_ss(&B[(l_n*88)+50]);    b50 = _mm_shuffle_ps(b50, b50, 0x00);
#endif
    __m128 c50_0 = _mm_load_ss(&C[(l_n*88)+4]);
    __m128 a50_0 = _mm_load_ss(&A[218]);
    c50_0 = _mm_add_ss(c50_0, _mm_mul_ss(a50_0, b50));
    _mm_store_ss(&C[(l_n*88)+4], c50_0);
    __m128 c50_1 = _mm_load_ss(&C[(l_n*88)+14]);
    __m128 a50_1 = _mm_load_ss(&A[219]);
    c50_1 = _mm_add_ss(c50_1, _mm_mul_ss(a50_1, b50));
    _mm_store_ss(&C[(l_n*88)+14], c50_1);
    __m128 c50_2 = _mm_load_ss(&C[(l_n*88)+29]);
    __m128 a50_2 = _mm_load_ss(&A[220]);
    c50_2 = _mm_add_ss(c50_2, _mm_mul_ss(a50_2, b50));
    _mm_store_ss(&C[(l_n*88)+29], c50_2);
    __m128 c50_3 = _mm_load_ss(&C[(l_n*88)+50]);
    __m128 a50_3 = _mm_load_ss(&A[221]);
    c50_3 = _mm_add_ss(c50_3, _mm_mul_ss(a50_3, b50));
    _mm_store_ss(&C[(l_n*88)+50], c50_3);
    __m128 c50_4 = _mm_load_ss(&C[(l_n*88)+78]);
    __m128 a50_4 = _mm_load_ss(&A[222]);
    c50_4 = _mm_add_ss(c50_4, _mm_mul_ss(a50_4, b50));
    _mm_store_ss(&C[(l_n*88)+78], c50_4);
#else
    C[(l_n*88)+4] += A[218] * B[(l_n*88)+50];
    C[(l_n*88)+14] += A[219] * B[(l_n*88)+50];
    C[(l_n*88)+29] += A[220] * B[(l_n*88)+50];
    C[(l_n*88)+50] += A[221] * B[(l_n*88)+50];
    C[(l_n*88)+78] += A[222] * B[(l_n*88)+50];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b51 = _mm_broadcast_ss(&B[(l_n*88)+51]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b51 = _mm_load_ss(&B[(l_n*88)+51]);    b51 = _mm_shuffle_ps(b51, b51, 0x00);
#endif
    __m128 c51_0 = _mm_load_ss(&C[(l_n*88)+5]);
    __m128 a51_0 = _mm_load_ss(&A[223]);
    c51_0 = _mm_add_ss(c51_0, _mm_mul_ss(a51_0, b51));
    _mm_store_ss(&C[(l_n*88)+5], c51_0);
    __m128 c51_1 = _mm_load_ss(&C[(l_n*88)+15]);
    __m128 a51_1 = _mm_load_ss(&A[224]);
    c51_1 = _mm_add_ss(c51_1, _mm_mul_ss(a51_1, b51));
    _mm_store_ss(&C[(l_n*88)+15], c51_1);
    __m128 c51_2 = _mm_load_ss(&C[(l_n*88)+30]);
    __m128 a51_2 = _mm_load_ss(&A[225]);
    c51_2 = _mm_add_ss(c51_2, _mm_mul_ss(a51_2, b51));
    _mm_store_ss(&C[(l_n*88)+30], c51_2);
    __m128 c51_3 = _mm_load_ss(&C[(l_n*88)+51]);
    __m128 a51_3 = _mm_load_ss(&A[226]);
    c51_3 = _mm_add_ss(c51_3, _mm_mul_ss(a51_3, b51));
    _mm_store_ss(&C[(l_n*88)+51], c51_3);
    __m128 c51_4 = _mm_load_ss(&C[(l_n*88)+79]);
    __m128 a51_4 = _mm_load_ss(&A[227]);
    c51_4 = _mm_add_ss(c51_4, _mm_mul_ss(a51_4, b51));
    _mm_store_ss(&C[(l_n*88)+79], c51_4);
#else
    C[(l_n*88)+5] += A[223] * B[(l_n*88)+51];
    C[(l_n*88)+15] += A[224] * B[(l_n*88)+51];
    C[(l_n*88)+30] += A[225] * B[(l_n*88)+51];
    C[(l_n*88)+51] += A[226] * B[(l_n*88)+51];
    C[(l_n*88)+79] += A[227] * B[(l_n*88)+51];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b52 = _mm_broadcast_ss(&B[(l_n*88)+52]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b52 = _mm_load_ss(&B[(l_n*88)+52]);    b52 = _mm_shuffle_ps(b52, b52, 0x00);
#endif
    __m128 c52_0 = _mm_load_ss(&C[(l_n*88)+6]);
    __m128 a52_0 = _mm_load_ss(&A[228]);
    c52_0 = _mm_add_ss(c52_0, _mm_mul_ss(a52_0, b52));
    _mm_store_ss(&C[(l_n*88)+6], c52_0);
    __m128 c52_1 = _mm_load_ss(&C[(l_n*88)+16]);
    __m128 a52_1 = _mm_load_ss(&A[229]);
    c52_1 = _mm_add_ss(c52_1, _mm_mul_ss(a52_1, b52));
    _mm_store_ss(&C[(l_n*88)+16], c52_1);
    __m128 c52_2 = _mm_load_ss(&C[(l_n*88)+31]);
    __m128 a52_2 = _mm_load_ss(&A[230]);
    c52_2 = _mm_add_ss(c52_2, _mm_mul_ss(a52_2, b52));
    _mm_store_ss(&C[(l_n*88)+31], c52_2);
    __m128 c52_3 = _mm_load_ss(&C[(l_n*88)+52]);
    __m128 a52_3 = _mm_load_ss(&A[231]);
    c52_3 = _mm_add_ss(c52_3, _mm_mul_ss(a52_3, b52));
    _mm_store_ss(&C[(l_n*88)+52], c52_3);
    __m128 c52_4 = _mm_load_ss(&C[(l_n*88)+80]);
    __m128 a52_4 = _mm_load_ss(&A[232]);
    c52_4 = _mm_add_ss(c52_4, _mm_mul_ss(a52_4, b52));
    _mm_store_ss(&C[(l_n*88)+80], c52_4);
#else
    C[(l_n*88)+6] += A[228] * B[(l_n*88)+52];
    C[(l_n*88)+16] += A[229] * B[(l_n*88)+52];
    C[(l_n*88)+31] += A[230] * B[(l_n*88)+52];
    C[(l_n*88)+52] += A[231] * B[(l_n*88)+52];
    C[(l_n*88)+80] += A[232] * B[(l_n*88)+52];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b53 = _mm_broadcast_ss(&B[(l_n*88)+53]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b53 = _mm_load_ss(&B[(l_n*88)+53]);    b53 = _mm_shuffle_ps(b53, b53, 0x00);
#endif
    __m128 c53_0 = _mm_load_ss(&C[(l_n*88)+1]);
    __m128 a53_0 = _mm_load_ss(&A[233]);
    c53_0 = _mm_add_ss(c53_0, _mm_mul_ss(a53_0, b53));
    _mm_store_ss(&C[(l_n*88)+1], c53_0);
    __m128 c53_1 = _mm_load_ss(&C[(l_n*88)+7]);
    __m128 a53_1 = _mm_load_ss(&A[234]);
    c53_1 = _mm_add_ss(c53_1, _mm_mul_ss(a53_1, b53));
    _mm_store_ss(&C[(l_n*88)+7], c53_1);
    __m128 c53_2 = _mm_load_ss(&C[(l_n*88)+17]);
    __m128 a53_2 = _mm_load_ss(&A[235]);
    c53_2 = _mm_add_ss(c53_2, _mm_mul_ss(a53_2, b53));
    _mm_store_ss(&C[(l_n*88)+17], c53_2);
    __m128 c53_3 = _mm_load_ss(&C[(l_n*88)+32]);
    __m128 a53_3 = _mm_load_ss(&A[236]);
    c53_3 = _mm_add_ss(c53_3, _mm_mul_ss(a53_3, b53));
    _mm_store_ss(&C[(l_n*88)+32], c53_3);
    __m128 c53_4 = _mm_load_ss(&C[(l_n*88)+53]);
    __m128 a53_4 = _mm_load_ss(&A[237]);
    c53_4 = _mm_add_ss(c53_4, _mm_mul_ss(a53_4, b53));
    _mm_store_ss(&C[(l_n*88)+53], c53_4);
    __m128 c53_5 = _mm_load_ss(&C[(l_n*88)+81]);
    __m128 a53_5 = _mm_load_ss(&A[238]);
    c53_5 = _mm_add_ss(c53_5, _mm_mul_ss(a53_5, b53));
    _mm_store_ss(&C[(l_n*88)+81], c53_5);
#else
    C[(l_n*88)+1] += A[233] * B[(l_n*88)+53];
    C[(l_n*88)+7] += A[234] * B[(l_n*88)+53];
    C[(l_n*88)+17] += A[235] * B[(l_n*88)+53];
    C[(l_n*88)+32] += A[236] * B[(l_n*88)+53];
    C[(l_n*88)+53] += A[237] * B[(l_n*88)+53];
    C[(l_n*88)+81] += A[238] * B[(l_n*88)+53];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b54 = _mm_broadcast_ss(&B[(l_n*88)+54]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b54 = _mm_load_ss(&B[(l_n*88)+54]);    b54 = _mm_shuffle_ps(b54, b54, 0x00);
#endif
    __m128 c54_0 = _mm_load_ss(&C[(l_n*88)+2]);
    __m128 a54_0 = _mm_load_ss(&A[239]);
    c54_0 = _mm_add_ss(c54_0, _mm_mul_ss(a54_0, b54));
    _mm_store_ss(&C[(l_n*88)+2], c54_0);
    __m128 c54_1 = _mm_load_ss(&C[(l_n*88)+8]);
    __m128 a54_1 = _mm_load_ss(&A[240]);
    c54_1 = _mm_add_ss(c54_1, _mm_mul_ss(a54_1, b54));
    _mm_store_ss(&C[(l_n*88)+8], c54_1);
    __m128 c54_2 = _mm_load_ss(&C[(l_n*88)+18]);
    __m128 a54_2 = _mm_load_ss(&A[241]);
    c54_2 = _mm_add_ss(c54_2, _mm_mul_ss(a54_2, b54));
    _mm_store_ss(&C[(l_n*88)+18], c54_2);
    __m128 c54_3 = _mm_load_ss(&C[(l_n*88)+33]);
    __m128 a54_3 = _mm_load_ss(&A[242]);
    c54_3 = _mm_add_ss(c54_3, _mm_mul_ss(a54_3, b54));
    _mm_store_ss(&C[(l_n*88)+33], c54_3);
    __m128 c54_4 = _mm_load_ss(&C[(l_n*88)+54]);
    __m128 a54_4 = _mm_load_ss(&A[243]);
    c54_4 = _mm_add_ss(c54_4, _mm_mul_ss(a54_4, b54));
    _mm_store_ss(&C[(l_n*88)+54], c54_4);
    __m128 c54_5 = _mm_load_ss(&C[(l_n*88)+82]);
    __m128 a54_5 = _mm_load_ss(&A[244]);
    c54_5 = _mm_add_ss(c54_5, _mm_mul_ss(a54_5, b54));
    _mm_store_ss(&C[(l_n*88)+82], c54_5);
#else
    C[(l_n*88)+2] += A[239] * B[(l_n*88)+54];
    C[(l_n*88)+8] += A[240] * B[(l_n*88)+54];
    C[(l_n*88)+18] += A[241] * B[(l_n*88)+54];
    C[(l_n*88)+33] += A[242] * B[(l_n*88)+54];
    C[(l_n*88)+54] += A[243] * B[(l_n*88)+54];
    C[(l_n*88)+82] += A[244] * B[(l_n*88)+54];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b55 = _mm_broadcast_ss(&B[(l_n*88)+55]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b55 = _mm_load_ss(&B[(l_n*88)+55]);    b55 = _mm_shuffle_ps(b55, b55, 0x00);
#endif
    __m128 c55_0 = _mm_load_ss(&C[(l_n*88)+0]);
    __m128 a55_0 = _mm_load_ss(&A[245]);
    c55_0 = _mm_add_ss(c55_0, _mm_mul_ss(a55_0, b55));
    _mm_store_ss(&C[(l_n*88)+0], c55_0);
    __m128 c55_1 = _mm_load_ss(&C[(l_n*88)+3]);
    __m128 a55_1 = _mm_load_ss(&A[246]);
    c55_1 = _mm_add_ss(c55_1, _mm_mul_ss(a55_1, b55));
    _mm_store_ss(&C[(l_n*88)+3], c55_1);
    __m128 c55_2 = _mm_load_ss(&C[(l_n*88)+9]);
    __m128 a55_2 = _mm_load_ss(&A[247]);
    c55_2 = _mm_add_ss(c55_2, _mm_mul_ss(a55_2, b55));
    _mm_store_ss(&C[(l_n*88)+9], c55_2);
    __m128 c55_3 = _mm_load_ss(&C[(l_n*88)+19]);
    __m128 a55_3 = _mm_load_ss(&A[248]);
    c55_3 = _mm_add_ss(c55_3, _mm_mul_ss(a55_3, b55));
    _mm_store_ss(&C[(l_n*88)+19], c55_3);
    __m128 c55_4 = _mm_load_ss(&C[(l_n*88)+34]);
    __m128 a55_4 = _mm_load_ss(&A[249]);
    c55_4 = _mm_add_ss(c55_4, _mm_mul_ss(a55_4, b55));
    _mm_store_ss(&C[(l_n*88)+34], c55_4);
    __m128 c55_5 = _mm_load_ss(&C[(l_n*88)+55]);
    __m128 a55_5 = _mm_load_ss(&A[250]);
    c55_5 = _mm_add_ss(c55_5, _mm_mul_ss(a55_5, b55));
    _mm_store_ss(&C[(l_n*88)+55], c55_5);
    __m128 c55_6 = _mm_load_ss(&C[(l_n*88)+83]);
    __m128 a55_6 = _mm_load_ss(&A[251]);
    c55_6 = _mm_add_ss(c55_6, _mm_mul_ss(a55_6, b55));
    _mm_store_ss(&C[(l_n*88)+83], c55_6);
#else
    C[(l_n*88)+0] += A[245] * B[(l_n*88)+55];
    C[(l_n*88)+3] += A[246] * B[(l_n*88)+55];
    C[(l_n*88)+9] += A[247] * B[(l_n*88)+55];
    C[(l_n*88)+19] += A[248] * B[(l_n*88)+55];
    C[(l_n*88)+34] += A[249] * B[(l_n*88)+55];
    C[(l_n*88)+55] += A[250] * B[(l_n*88)+55];
    C[(l_n*88)+83] += A[251] * B[(l_n*88)+55];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b56 = _mm_broadcast_ss(&B[(l_n*88)+56]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b56 = _mm_load_ss(&B[(l_n*88)+56]);    b56 = _mm_shuffle_ps(b56, b56, 0x00);
#endif
    __m128 c56_0 = _mm_load_ss(&C[(l_n*88)+56]);
    __m128 a56_0 = _mm_load_ss(&A[252]);
    c56_0 = _mm_add_ss(c56_0, _mm_mul_ss(a56_0, b56));
    _mm_store_ss(&C[(l_n*88)+56], c56_0);
#else
    C[(l_n*88)+56] += A[252] * B[(l_n*88)+56];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b57 = _mm_broadcast_ss(&B[(l_n*88)+57]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b57 = _mm_load_ss(&B[(l_n*88)+57]);    b57 = _mm_shuffle_ps(b57, b57, 0x00);
#endif
    __m128 c57_0 = _mm_load_ss(&C[(l_n*88)+57]);
    __m128 a57_0 = _mm_load_ss(&A[253]);
    c57_0 = _mm_add_ss(c57_0, _mm_mul_ss(a57_0, b57));
    _mm_store_ss(&C[(l_n*88)+57], c57_0);
#else
    C[(l_n*88)+57] += A[253] * B[(l_n*88)+57];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b58 = _mm_broadcast_ss(&B[(l_n*88)+58]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b58 = _mm_load_ss(&B[(l_n*88)+58]);    b58 = _mm_shuffle_ps(b58, b58, 0x00);
#endif
    __m128 c58_0 = _mm_load_ss(&C[(l_n*88)+58]);
    __m128 a58_0 = _mm_load_ss(&A[254]);
    c58_0 = _mm_add_ss(c58_0, _mm_mul_ss(a58_0, b58));
    _mm_store_ss(&C[(l_n*88)+58], c58_0);
#else
    C[(l_n*88)+58] += A[254] * B[(l_n*88)+58];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b59 = _mm_broadcast_ss(&B[(l_n*88)+59]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b59 = _mm_load_ss(&B[(l_n*88)+59]);    b59 = _mm_shuffle_ps(b59, b59, 0x00);
#endif
    __m128 c59_0 = _mm_load_ss(&C[(l_n*88)+59]);
    __m128 a59_0 = _mm_load_ss(&A[255]);
    c59_0 = _mm_add_ss(c59_0, _mm_mul_ss(a59_0, b59));
    _mm_store_ss(&C[(l_n*88)+59], c59_0);
#else
    C[(l_n*88)+59] += A[255] * B[(l_n*88)+59];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b60 = _mm_broadcast_ss(&B[(l_n*88)+60]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b60 = _mm_load_ss(&B[(l_n*88)+60]);    b60 = _mm_shuffle_ps(b60, b60, 0x00);
#endif
    __m128 c60_0 = _mm_load_ss(&C[(l_n*88)+60]);
    __m128 a60_0 = _mm_load_ss(&A[256]);
    c60_0 = _mm_add_ss(c60_0, _mm_mul_ss(a60_0, b60));
    _mm_store_ss(&C[(l_n*88)+60], c60_0);
#else
    C[(l_n*88)+60] += A[256] * B[(l_n*88)+60];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b61 = _mm_broadcast_ss(&B[(l_n*88)+61]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b61 = _mm_load_ss(&B[(l_n*88)+61]);    b61 = _mm_shuffle_ps(b61, b61, 0x00);
#endif
    __m128 c61_0 = _mm_load_ss(&C[(l_n*88)+61]);
    __m128 a61_0 = _mm_load_ss(&A[257]);
    c61_0 = _mm_add_ss(c61_0, _mm_mul_ss(a61_0, b61));
    _mm_store_ss(&C[(l_n*88)+61], c61_0);
#else
    C[(l_n*88)+61] += A[257] * B[(l_n*88)+61];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b62 = _mm_broadcast_ss(&B[(l_n*88)+62]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b62 = _mm_load_ss(&B[(l_n*88)+62]);    b62 = _mm_shuffle_ps(b62, b62, 0x00);
#endif
    __m128 c62_0 = _mm_load_ss(&C[(l_n*88)+62]);
    __m128 a62_0 = _mm_load_ss(&A[258]);
    c62_0 = _mm_add_ss(c62_0, _mm_mul_ss(a62_0, b62));
    _mm_store_ss(&C[(l_n*88)+62], c62_0);
#else
    C[(l_n*88)+62] += A[258] * B[(l_n*88)+62];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b63 = _mm_broadcast_ss(&B[(l_n*88)+63]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b63 = _mm_load_ss(&B[(l_n*88)+63]);    b63 = _mm_shuffle_ps(b63, b63, 0x00);
#endif
    __m128 c63_0 = _mm_load_ss(&C[(l_n*88)+35]);
    __m128 a63_0 = _mm_load_ss(&A[259]);
    c63_0 = _mm_add_ss(c63_0, _mm_mul_ss(a63_0, b63));
    _mm_store_ss(&C[(l_n*88)+35], c63_0);
    __m128 c63_1 = _mm_load_ss(&C[(l_n*88)+63]);
    __m128 a63_1 = _mm_load_ss(&A[260]);
    c63_1 = _mm_add_ss(c63_1, _mm_mul_ss(a63_1, b63));
    _mm_store_ss(&C[(l_n*88)+63], c63_1);
#else
    C[(l_n*88)+35] += A[259] * B[(l_n*88)+63];
    C[(l_n*88)+63] += A[260] * B[(l_n*88)+63];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b64 = _mm_broadcast_ss(&B[(l_n*88)+64]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b64 = _mm_load_ss(&B[(l_n*88)+64]);    b64 = _mm_shuffle_ps(b64, b64, 0x00);
#endif
    __m128 c64_0 = _mm_load_ss(&C[(l_n*88)+36]);
    __m128 a64_0 = _mm_load_ss(&A[261]);
    c64_0 = _mm_add_ss(c64_0, _mm_mul_ss(a64_0, b64));
    _mm_store_ss(&C[(l_n*88)+36], c64_0);
    __m128 c64_1 = _mm_load_ss(&C[(l_n*88)+64]);
    __m128 a64_1 = _mm_load_ss(&A[262]);
    c64_1 = _mm_add_ss(c64_1, _mm_mul_ss(a64_1, b64));
    _mm_store_ss(&C[(l_n*88)+64], c64_1);
#else
    C[(l_n*88)+36] += A[261] * B[(l_n*88)+64];
    C[(l_n*88)+64] += A[262] * B[(l_n*88)+64];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b65 = _mm_broadcast_ss(&B[(l_n*88)+65]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b65 = _mm_load_ss(&B[(l_n*88)+65]);    b65 = _mm_shuffle_ps(b65, b65, 0x00);
#endif
    __m128 c65_0 = _mm_load_ss(&C[(l_n*88)+37]);
    __m128 a65_0 = _mm_load_ss(&A[263]);
    c65_0 = _mm_add_ss(c65_0, _mm_mul_ss(a65_0, b65));
    _mm_store_ss(&C[(l_n*88)+37], c65_0);
    __m128 c65_1 = _mm_load_ss(&C[(l_n*88)+65]);
    __m128 a65_1 = _mm_load_ss(&A[264]);
    c65_1 = _mm_add_ss(c65_1, _mm_mul_ss(a65_1, b65));
    _mm_store_ss(&C[(l_n*88)+65], c65_1);
#else
    C[(l_n*88)+37] += A[263] * B[(l_n*88)+65];
    C[(l_n*88)+65] += A[264] * B[(l_n*88)+65];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b66 = _mm_broadcast_ss(&B[(l_n*88)+66]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b66 = _mm_load_ss(&B[(l_n*88)+66]);    b66 = _mm_shuffle_ps(b66, b66, 0x00);
#endif
    __m128 c66_0 = _mm_load_ss(&C[(l_n*88)+38]);
    __m128 a66_0 = _mm_load_ss(&A[265]);
    c66_0 = _mm_add_ss(c66_0, _mm_mul_ss(a66_0, b66));
    _mm_store_ss(&C[(l_n*88)+38], c66_0);
    __m128 c66_1 = _mm_load_ss(&C[(l_n*88)+66]);
    __m128 a66_1 = _mm_load_ss(&A[266]);
    c66_1 = _mm_add_ss(c66_1, _mm_mul_ss(a66_1, b66));
    _mm_store_ss(&C[(l_n*88)+66], c66_1);
#else
    C[(l_n*88)+38] += A[265] * B[(l_n*88)+66];
    C[(l_n*88)+66] += A[266] * B[(l_n*88)+66];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b67 = _mm_broadcast_ss(&B[(l_n*88)+67]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b67 = _mm_load_ss(&B[(l_n*88)+67]);    b67 = _mm_shuffle_ps(b67, b67, 0x00);
#endif
    __m128 c67_0 = _mm_load_ss(&C[(l_n*88)+39]);
    __m128 a67_0 = _mm_load_ss(&A[267]);
    c67_0 = _mm_add_ss(c67_0, _mm_mul_ss(a67_0, b67));
    _mm_store_ss(&C[(l_n*88)+39], c67_0);
    __m128 c67_1 = _mm_load_ss(&C[(l_n*88)+67]);
    __m128 a67_1 = _mm_load_ss(&A[268]);
    c67_1 = _mm_add_ss(c67_1, _mm_mul_ss(a67_1, b67));
    _mm_store_ss(&C[(l_n*88)+67], c67_1);
#else
    C[(l_n*88)+39] += A[267] * B[(l_n*88)+67];
    C[(l_n*88)+67] += A[268] * B[(l_n*88)+67];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b68 = _mm_broadcast_ss(&B[(l_n*88)+68]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b68 = _mm_load_ss(&B[(l_n*88)+68]);    b68 = _mm_shuffle_ps(b68, b68, 0x00);
#endif
    __m128 c68_0 = _mm_load_ss(&C[(l_n*88)+40]);
    __m128 a68_0 = _mm_load_ss(&A[269]);
    c68_0 = _mm_add_ss(c68_0, _mm_mul_ss(a68_0, b68));
    _mm_store_ss(&C[(l_n*88)+40], c68_0);
    __m128 c68_1 = _mm_load_ss(&C[(l_n*88)+68]);
    __m128 a68_1 = _mm_load_ss(&A[270]);
    c68_1 = _mm_add_ss(c68_1, _mm_mul_ss(a68_1, b68));
    _mm_store_ss(&C[(l_n*88)+68], c68_1);
#else
    C[(l_n*88)+40] += A[269] * B[(l_n*88)+68];
    C[(l_n*88)+68] += A[270] * B[(l_n*88)+68];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b69 = _mm_broadcast_ss(&B[(l_n*88)+69]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b69 = _mm_load_ss(&B[(l_n*88)+69]);    b69 = _mm_shuffle_ps(b69, b69, 0x00);
#endif
    __m128 c69_0 = _mm_load_ss(&C[(l_n*88)+20]);
    __m128 a69_0 = _mm_load_ss(&A[271]);
    c69_0 = _mm_add_ss(c69_0, _mm_mul_ss(a69_0, b69));
    _mm_store_ss(&C[(l_n*88)+20], c69_0);
    __m128 c69_1 = _mm_load_ss(&C[(l_n*88)+41]);
    __m128 a69_1 = _mm_load_ss(&A[272]);
    c69_1 = _mm_add_ss(c69_1, _mm_mul_ss(a69_1, b69));
    _mm_store_ss(&C[(l_n*88)+41], c69_1);
    __m128 c69_2 = _mm_load_ss(&C[(l_n*88)+69]);
    __m128 a69_2 = _mm_load_ss(&A[273]);
    c69_2 = _mm_add_ss(c69_2, _mm_mul_ss(a69_2, b69));
    _mm_store_ss(&C[(l_n*88)+69], c69_2);
#else
    C[(l_n*88)+20] += A[271] * B[(l_n*88)+69];
    C[(l_n*88)+41] += A[272] * B[(l_n*88)+69];
    C[(l_n*88)+69] += A[273] * B[(l_n*88)+69];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b70 = _mm_broadcast_ss(&B[(l_n*88)+70]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b70 = _mm_load_ss(&B[(l_n*88)+70]);    b70 = _mm_shuffle_ps(b70, b70, 0x00);
#endif
    __m128 c70_0 = _mm_load_ss(&C[(l_n*88)+21]);
    __m128 a70_0 = _mm_load_ss(&A[274]);
    c70_0 = _mm_add_ss(c70_0, _mm_mul_ss(a70_0, b70));
    _mm_store_ss(&C[(l_n*88)+21], c70_0);
    __m128 c70_1 = _mm_load_ss(&C[(l_n*88)+42]);
    __m128 a70_1 = _mm_load_ss(&A[275]);
    c70_1 = _mm_add_ss(c70_1, _mm_mul_ss(a70_1, b70));
    _mm_store_ss(&C[(l_n*88)+42], c70_1);
    __m128 c70_2 = _mm_load_ss(&C[(l_n*88)+70]);
    __m128 a70_2 = _mm_load_ss(&A[276]);
    c70_2 = _mm_add_ss(c70_2, _mm_mul_ss(a70_2, b70));
    _mm_store_ss(&C[(l_n*88)+70], c70_2);
#else
    C[(l_n*88)+21] += A[274] * B[(l_n*88)+70];
    C[(l_n*88)+42] += A[275] * B[(l_n*88)+70];
    C[(l_n*88)+70] += A[276] * B[(l_n*88)+70];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b71 = _mm_broadcast_ss(&B[(l_n*88)+71]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b71 = _mm_load_ss(&B[(l_n*88)+71]);    b71 = _mm_shuffle_ps(b71, b71, 0x00);
#endif
    __m128 c71_0 = _mm_load_ss(&C[(l_n*88)+22]);
    __m128 a71_0 = _mm_load_ss(&A[277]);
    c71_0 = _mm_add_ss(c71_0, _mm_mul_ss(a71_0, b71));
    _mm_store_ss(&C[(l_n*88)+22], c71_0);
    __m128 c71_1 = _mm_load_ss(&C[(l_n*88)+43]);
    __m128 a71_1 = _mm_load_ss(&A[278]);
    c71_1 = _mm_add_ss(c71_1, _mm_mul_ss(a71_1, b71));
    _mm_store_ss(&C[(l_n*88)+43], c71_1);
    __m128 c71_2 = _mm_load_ss(&C[(l_n*88)+71]);
    __m128 a71_2 = _mm_load_ss(&A[279]);
    c71_2 = _mm_add_ss(c71_2, _mm_mul_ss(a71_2, b71));
    _mm_store_ss(&C[(l_n*88)+71], c71_2);
#else
    C[(l_n*88)+22] += A[277] * B[(l_n*88)+71];
    C[(l_n*88)+43] += A[278] * B[(l_n*88)+71];
    C[(l_n*88)+71] += A[279] * B[(l_n*88)+71];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b72 = _mm_broadcast_ss(&B[(l_n*88)+72]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b72 = _mm_load_ss(&B[(l_n*88)+72]);    b72 = _mm_shuffle_ps(b72, b72, 0x00);
#endif
    __m128 c72_0 = _mm_load_ss(&C[(l_n*88)+23]);
    __m128 a72_0 = _mm_load_ss(&A[280]);
    c72_0 = _mm_add_ss(c72_0, _mm_mul_ss(a72_0, b72));
    _mm_store_ss(&C[(l_n*88)+23], c72_0);
    __m128 c72_1 = _mm_load_ss(&C[(l_n*88)+44]);
    __m128 a72_1 = _mm_load_ss(&A[281]);
    c72_1 = _mm_add_ss(c72_1, _mm_mul_ss(a72_1, b72));
    _mm_store_ss(&C[(l_n*88)+44], c72_1);
    __m128 c72_2 = _mm_load_ss(&C[(l_n*88)+72]);
    __m128 a72_2 = _mm_load_ss(&A[282]);
    c72_2 = _mm_add_ss(c72_2, _mm_mul_ss(a72_2, b72));
    _mm_store_ss(&C[(l_n*88)+72], c72_2);
#else
    C[(l_n*88)+23] += A[280] * B[(l_n*88)+72];
    C[(l_n*88)+44] += A[281] * B[(l_n*88)+72];
    C[(l_n*88)+72] += A[282] * B[(l_n*88)+72];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b73 = _mm_broadcast_ss(&B[(l_n*88)+73]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b73 = _mm_load_ss(&B[(l_n*88)+73]);    b73 = _mm_shuffle_ps(b73, b73, 0x00);
#endif
    __m128 c73_0 = _mm_load_ss(&C[(l_n*88)+24]);
    __m128 a73_0 = _mm_load_ss(&A[283]);
    c73_0 = _mm_add_ss(c73_0, _mm_mul_ss(a73_0, b73));
    _mm_store_ss(&C[(l_n*88)+24], c73_0);
    __m128 c73_1 = _mm_load_ss(&C[(l_n*88)+45]);
    __m128 a73_1 = _mm_load_ss(&A[284]);
    c73_1 = _mm_add_ss(c73_1, _mm_mul_ss(a73_1, b73));
    _mm_store_ss(&C[(l_n*88)+45], c73_1);
    __m128 c73_2 = _mm_load_ss(&C[(l_n*88)+73]);
    __m128 a73_2 = _mm_load_ss(&A[285]);
    c73_2 = _mm_add_ss(c73_2, _mm_mul_ss(a73_2, b73));
    _mm_store_ss(&C[(l_n*88)+73], c73_2);
#else
    C[(l_n*88)+24] += A[283] * B[(l_n*88)+73];
    C[(l_n*88)+45] += A[284] * B[(l_n*88)+73];
    C[(l_n*88)+73] += A[285] * B[(l_n*88)+73];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b74 = _mm_broadcast_ss(&B[(l_n*88)+74]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b74 = _mm_load_ss(&B[(l_n*88)+74]);    b74 = _mm_shuffle_ps(b74, b74, 0x00);
#endif
    __m128 c74_0 = _mm_load_ss(&C[(l_n*88)+10]);
    __m128 a74_0 = _mm_load_ss(&A[286]);
    c74_0 = _mm_add_ss(c74_0, _mm_mul_ss(a74_0, b74));
    _mm_store_ss(&C[(l_n*88)+10], c74_0);
    __m128 c74_1 = _mm_load_ss(&C[(l_n*88)+25]);
    __m128 a74_1 = _mm_load_ss(&A[287]);
    c74_1 = _mm_add_ss(c74_1, _mm_mul_ss(a74_1, b74));
    _mm_store_ss(&C[(l_n*88)+25], c74_1);
    __m128 c74_2 = _mm_load_ss(&C[(l_n*88)+46]);
    __m128 a74_2 = _mm_load_ss(&A[288]);
    c74_2 = _mm_add_ss(c74_2, _mm_mul_ss(a74_2, b74));
    _mm_store_ss(&C[(l_n*88)+46], c74_2);
    __m128 c74_3 = _mm_load_ss(&C[(l_n*88)+74]);
    __m128 a74_3 = _mm_load_ss(&A[289]);
    c74_3 = _mm_add_ss(c74_3, _mm_mul_ss(a74_3, b74));
    _mm_store_ss(&C[(l_n*88)+74], c74_3);
#else
    C[(l_n*88)+10] += A[286] * B[(l_n*88)+74];
    C[(l_n*88)+25] += A[287] * B[(l_n*88)+74];
    C[(l_n*88)+46] += A[288] * B[(l_n*88)+74];
    C[(l_n*88)+74] += A[289] * B[(l_n*88)+74];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b75 = _mm_broadcast_ss(&B[(l_n*88)+75]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b75 = _mm_load_ss(&B[(l_n*88)+75]);    b75 = _mm_shuffle_ps(b75, b75, 0x00);
#endif
    __m128 c75_0 = _mm_load_ss(&C[(l_n*88)+11]);
    __m128 a75_0 = _mm_load_ss(&A[290]);
    c75_0 = _mm_add_ss(c75_0, _mm_mul_ss(a75_0, b75));
    _mm_store_ss(&C[(l_n*88)+11], c75_0);
    __m128 c75_1 = _mm_load_ss(&C[(l_n*88)+26]);
    __m128 a75_1 = _mm_load_ss(&A[291]);
    c75_1 = _mm_add_ss(c75_1, _mm_mul_ss(a75_1, b75));
    _mm_store_ss(&C[(l_n*88)+26], c75_1);
    __m128 c75_2 = _mm_load_ss(&C[(l_n*88)+47]);
    __m128 a75_2 = _mm_load_ss(&A[292]);
    c75_2 = _mm_add_ss(c75_2, _mm_mul_ss(a75_2, b75));
    _mm_store_ss(&C[(l_n*88)+47], c75_2);
    __m128 c75_3 = _mm_load_ss(&C[(l_n*88)+75]);
    __m128 a75_3 = _mm_load_ss(&A[293]);
    c75_3 = _mm_add_ss(c75_3, _mm_mul_ss(a75_3, b75));
    _mm_store_ss(&C[(l_n*88)+75], c75_3);
#else
    C[(l_n*88)+11] += A[290] * B[(l_n*88)+75];
    C[(l_n*88)+26] += A[291] * B[(l_n*88)+75];
    C[(l_n*88)+47] += A[292] * B[(l_n*88)+75];
    C[(l_n*88)+75] += A[293] * B[(l_n*88)+75];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b76 = _mm_broadcast_ss(&B[(l_n*88)+76]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b76 = _mm_load_ss(&B[(l_n*88)+76]);    b76 = _mm_shuffle_ps(b76, b76, 0x00);
#endif
    __m128 c76_0 = _mm_load_ss(&C[(l_n*88)+12]);
    __m128 a76_0 = _mm_load_ss(&A[294]);
    c76_0 = _mm_add_ss(c76_0, _mm_mul_ss(a76_0, b76));
    _mm_store_ss(&C[(l_n*88)+12], c76_0);
    __m128 c76_1 = _mm_load_ss(&C[(l_n*88)+27]);
    __m128 a76_1 = _mm_load_ss(&A[295]);
    c76_1 = _mm_add_ss(c76_1, _mm_mul_ss(a76_1, b76));
    _mm_store_ss(&C[(l_n*88)+27], c76_1);
    __m128 c76_2 = _mm_load_ss(&C[(l_n*88)+48]);
    __m128 a76_2 = _mm_load_ss(&A[296]);
    c76_2 = _mm_add_ss(c76_2, _mm_mul_ss(a76_2, b76));
    _mm_store_ss(&C[(l_n*88)+48], c76_2);
    __m128 c76_3 = _mm_load_ss(&C[(l_n*88)+76]);
    __m128 a76_3 = _mm_load_ss(&A[297]);
    c76_3 = _mm_add_ss(c76_3, _mm_mul_ss(a76_3, b76));
    _mm_store_ss(&C[(l_n*88)+76], c76_3);
#else
    C[(l_n*88)+12] += A[294] * B[(l_n*88)+76];
    C[(l_n*88)+27] += A[295] * B[(l_n*88)+76];
    C[(l_n*88)+48] += A[296] * B[(l_n*88)+76];
    C[(l_n*88)+76] += A[297] * B[(l_n*88)+76];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b77 = _mm_broadcast_ss(&B[(l_n*88)+77]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b77 = _mm_load_ss(&B[(l_n*88)+77]);    b77 = _mm_shuffle_ps(b77, b77, 0x00);
#endif
    __m128 c77_0 = _mm_load_ss(&C[(l_n*88)+13]);
    __m128 a77_0 = _mm_load_ss(&A[298]);
    c77_0 = _mm_add_ss(c77_0, _mm_mul_ss(a77_0, b77));
    _mm_store_ss(&C[(l_n*88)+13], c77_0);
    __m128 c77_1 = _mm_load_ss(&C[(l_n*88)+28]);
    __m128 a77_1 = _mm_load_ss(&A[299]);
    c77_1 = _mm_add_ss(c77_1, _mm_mul_ss(a77_1, b77));
    _mm_store_ss(&C[(l_n*88)+28], c77_1);
    __m128 c77_2 = _mm_load_ss(&C[(l_n*88)+49]);
    __m128 a77_2 = _mm_load_ss(&A[300]);
    c77_2 = _mm_add_ss(c77_2, _mm_mul_ss(a77_2, b77));
    _mm_store_ss(&C[(l_n*88)+49], c77_2);
    __m128 c77_3 = _mm_load_ss(&C[(l_n*88)+77]);
    __m128 a77_3 = _mm_load_ss(&A[301]);
    c77_3 = _mm_add_ss(c77_3, _mm_mul_ss(a77_3, b77));
    _mm_store_ss(&C[(l_n*88)+77], c77_3);
#else
    C[(l_n*88)+13] += A[298] * B[(l_n*88)+77];
    C[(l_n*88)+28] += A[299] * B[(l_n*88)+77];
    C[(l_n*88)+49] += A[300] * B[(l_n*88)+77];
    C[(l_n*88)+77] += A[301] * B[(l_n*88)+77];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b78 = _mm_broadcast_ss(&B[(l_n*88)+78]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b78 = _mm_load_ss(&B[(l_n*88)+78]);    b78 = _mm_shuffle_ps(b78, b78, 0x00);
#endif
    __m128 c78_0 = _mm_load_ss(&C[(l_n*88)+4]);
    __m128 a78_0 = _mm_load_ss(&A[302]);
    c78_0 = _mm_add_ss(c78_0, _mm_mul_ss(a78_0, b78));
    _mm_store_ss(&C[(l_n*88)+4], c78_0);
    __m128 c78_1 = _mm_load_ss(&C[(l_n*88)+14]);
    __m128 a78_1 = _mm_load_ss(&A[303]);
    c78_1 = _mm_add_ss(c78_1, _mm_mul_ss(a78_1, b78));
    _mm_store_ss(&C[(l_n*88)+14], c78_1);
    __m128 c78_2 = _mm_load_ss(&C[(l_n*88)+29]);
    __m128 a78_2 = _mm_load_ss(&A[304]);
    c78_2 = _mm_add_ss(c78_2, _mm_mul_ss(a78_2, b78));
    _mm_store_ss(&C[(l_n*88)+29], c78_2);
    __m128 c78_3 = _mm_load_ss(&C[(l_n*88)+50]);
    __m128 a78_3 = _mm_load_ss(&A[305]);
    c78_3 = _mm_add_ss(c78_3, _mm_mul_ss(a78_3, b78));
    _mm_store_ss(&C[(l_n*88)+50], c78_3);
    __m128 c78_4 = _mm_load_ss(&C[(l_n*88)+78]);
    __m128 a78_4 = _mm_load_ss(&A[306]);
    c78_4 = _mm_add_ss(c78_4, _mm_mul_ss(a78_4, b78));
    _mm_store_ss(&C[(l_n*88)+78], c78_4);
#else
    C[(l_n*88)+4] += A[302] * B[(l_n*88)+78];
    C[(l_n*88)+14] += A[303] * B[(l_n*88)+78];
    C[(l_n*88)+29] += A[304] * B[(l_n*88)+78];
    C[(l_n*88)+50] += A[305] * B[(l_n*88)+78];
    C[(l_n*88)+78] += A[306] * B[(l_n*88)+78];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b79 = _mm_broadcast_ss(&B[(l_n*88)+79]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b79 = _mm_load_ss(&B[(l_n*88)+79]);    b79 = _mm_shuffle_ps(b79, b79, 0x00);
#endif
    __m128 c79_0 = _mm_load_ss(&C[(l_n*88)+5]);
    __m128 a79_0 = _mm_load_ss(&A[307]);
    c79_0 = _mm_add_ss(c79_0, _mm_mul_ss(a79_0, b79));
    _mm_store_ss(&C[(l_n*88)+5], c79_0);
    __m128 c79_1 = _mm_load_ss(&C[(l_n*88)+15]);
    __m128 a79_1 = _mm_load_ss(&A[308]);
    c79_1 = _mm_add_ss(c79_1, _mm_mul_ss(a79_1, b79));
    _mm_store_ss(&C[(l_n*88)+15], c79_1);
    __m128 c79_2 = _mm_load_ss(&C[(l_n*88)+30]);
    __m128 a79_2 = _mm_load_ss(&A[309]);
    c79_2 = _mm_add_ss(c79_2, _mm_mul_ss(a79_2, b79));
    _mm_store_ss(&C[(l_n*88)+30], c79_2);
    __m128 c79_3 = _mm_load_ss(&C[(l_n*88)+51]);
    __m128 a79_3 = _mm_load_ss(&A[310]);
    c79_3 = _mm_add_ss(c79_3, _mm_mul_ss(a79_3, b79));
    _mm_store_ss(&C[(l_n*88)+51], c79_3);
    __m128 c79_4 = _mm_load_ss(&C[(l_n*88)+79]);
    __m128 a79_4 = _mm_load_ss(&A[311]);
    c79_4 = _mm_add_ss(c79_4, _mm_mul_ss(a79_4, b79));
    _mm_store_ss(&C[(l_n*88)+79], c79_4);
#else
    C[(l_n*88)+5] += A[307] * B[(l_n*88)+79];
    C[(l_n*88)+15] += A[308] * B[(l_n*88)+79];
    C[(l_n*88)+30] += A[309] * B[(l_n*88)+79];
    C[(l_n*88)+51] += A[310] * B[(l_n*88)+79];
    C[(l_n*88)+79] += A[311] * B[(l_n*88)+79];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b80 = _mm_broadcast_ss(&B[(l_n*88)+80]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b80 = _mm_load_ss(&B[(l_n*88)+80]);    b80 = _mm_shuffle_ps(b80, b80, 0x00);
#endif
    __m128 c80_0 = _mm_load_ss(&C[(l_n*88)+6]);
    __m128 a80_0 = _mm_load_ss(&A[312]);
    c80_0 = _mm_add_ss(c80_0, _mm_mul_ss(a80_0, b80));
    _mm_store_ss(&C[(l_n*88)+6], c80_0);
    __m128 c80_1 = _mm_load_ss(&C[(l_n*88)+16]);
    __m128 a80_1 = _mm_load_ss(&A[313]);
    c80_1 = _mm_add_ss(c80_1, _mm_mul_ss(a80_1, b80));
    _mm_store_ss(&C[(l_n*88)+16], c80_1);
    __m128 c80_2 = _mm_load_ss(&C[(l_n*88)+31]);
    __m128 a80_2 = _mm_load_ss(&A[314]);
    c80_2 = _mm_add_ss(c80_2, _mm_mul_ss(a80_2, b80));
    _mm_store_ss(&C[(l_n*88)+31], c80_2);
    __m128 c80_3 = _mm_load_ss(&C[(l_n*88)+52]);
    __m128 a80_3 = _mm_load_ss(&A[315]);
    c80_3 = _mm_add_ss(c80_3, _mm_mul_ss(a80_3, b80));
    _mm_store_ss(&C[(l_n*88)+52], c80_3);
    __m128 c80_4 = _mm_load_ss(&C[(l_n*88)+80]);
    __m128 a80_4 = _mm_load_ss(&A[316]);
    c80_4 = _mm_add_ss(c80_4, _mm_mul_ss(a80_4, b80));
    _mm_store_ss(&C[(l_n*88)+80], c80_4);
#else
    C[(l_n*88)+6] += A[312] * B[(l_n*88)+80];
    C[(l_n*88)+16] += A[313] * B[(l_n*88)+80];
    C[(l_n*88)+31] += A[314] * B[(l_n*88)+80];
    C[(l_n*88)+52] += A[315] * B[(l_n*88)+80];
    C[(l_n*88)+80] += A[316] * B[(l_n*88)+80];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b81 = _mm_broadcast_ss(&B[(l_n*88)+81]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b81 = _mm_load_ss(&B[(l_n*88)+81]);    b81 = _mm_shuffle_ps(b81, b81, 0x00);
#endif
    __m128 c81_0 = _mm_load_ss(&C[(l_n*88)+1]);
    __m128 a81_0 = _mm_load_ss(&A[317]);
    c81_0 = _mm_add_ss(c81_0, _mm_mul_ss(a81_0, b81));
    _mm_store_ss(&C[(l_n*88)+1], c81_0);
    __m128 c81_1 = _mm_load_ss(&C[(l_n*88)+7]);
    __m128 a81_1 = _mm_load_ss(&A[318]);
    c81_1 = _mm_add_ss(c81_1, _mm_mul_ss(a81_1, b81));
    _mm_store_ss(&C[(l_n*88)+7], c81_1);
    __m128 c81_2 = _mm_load_ss(&C[(l_n*88)+17]);
    __m128 a81_2 = _mm_load_ss(&A[319]);
    c81_2 = _mm_add_ss(c81_2, _mm_mul_ss(a81_2, b81));
    _mm_store_ss(&C[(l_n*88)+17], c81_2);
    __m128 c81_3 = _mm_load_ss(&C[(l_n*88)+32]);
    __m128 a81_3 = _mm_load_ss(&A[320]);
    c81_3 = _mm_add_ss(c81_3, _mm_mul_ss(a81_3, b81));
    _mm_store_ss(&C[(l_n*88)+32], c81_3);
    __m128 c81_4 = _mm_load_ss(&C[(l_n*88)+53]);
    __m128 a81_4 = _mm_load_ss(&A[321]);
    c81_4 = _mm_add_ss(c81_4, _mm_mul_ss(a81_4, b81));
    _mm_store_ss(&C[(l_n*88)+53], c81_4);
    __m128 c81_5 = _mm_load_ss(&C[(l_n*88)+81]);
    __m128 a81_5 = _mm_load_ss(&A[322]);
    c81_5 = _mm_add_ss(c81_5, _mm_mul_ss(a81_5, b81));
    _mm_store_ss(&C[(l_n*88)+81], c81_5);
#else
    C[(l_n*88)+1] += A[317] * B[(l_n*88)+81];
    C[(l_n*88)+7] += A[318] * B[(l_n*88)+81];
    C[(l_n*88)+17] += A[319] * B[(l_n*88)+81];
    C[(l_n*88)+32] += A[320] * B[(l_n*88)+81];
    C[(l_n*88)+53] += A[321] * B[(l_n*88)+81];
    C[(l_n*88)+81] += A[322] * B[(l_n*88)+81];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b82 = _mm_broadcast_ss(&B[(l_n*88)+82]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b82 = _mm_load_ss(&B[(l_n*88)+82]);    b82 = _mm_shuffle_ps(b82, b82, 0x00);
#endif
    __m128 c82_0 = _mm_load_ss(&C[(l_n*88)+2]);
    __m128 a82_0 = _mm_load_ss(&A[323]);
    c82_0 = _mm_add_ss(c82_0, _mm_mul_ss(a82_0, b82));
    _mm_store_ss(&C[(l_n*88)+2], c82_0);
    __m128 c82_1 = _mm_load_ss(&C[(l_n*88)+8]);
    __m128 a82_1 = _mm_load_ss(&A[324]);
    c82_1 = _mm_add_ss(c82_1, _mm_mul_ss(a82_1, b82));
    _mm_store_ss(&C[(l_n*88)+8], c82_1);
    __m128 c82_2 = _mm_load_ss(&C[(l_n*88)+18]);
    __m128 a82_2 = _mm_load_ss(&A[325]);
    c82_2 = _mm_add_ss(c82_2, _mm_mul_ss(a82_2, b82));
    _mm_store_ss(&C[(l_n*88)+18], c82_2);
    __m128 c82_3 = _mm_load_ss(&C[(l_n*88)+33]);
    __m128 a82_3 = _mm_load_ss(&A[326]);
    c82_3 = _mm_add_ss(c82_3, _mm_mul_ss(a82_3, b82));
    _mm_store_ss(&C[(l_n*88)+33], c82_3);
    __m128 c82_4 = _mm_load_ss(&C[(l_n*88)+54]);
    __m128 a82_4 = _mm_load_ss(&A[327]);
    c82_4 = _mm_add_ss(c82_4, _mm_mul_ss(a82_4, b82));
    _mm_store_ss(&C[(l_n*88)+54], c82_4);
    __m128 c82_5 = _mm_load_ss(&C[(l_n*88)+82]);
    __m128 a82_5 = _mm_load_ss(&A[328]);
    c82_5 = _mm_add_ss(c82_5, _mm_mul_ss(a82_5, b82));
    _mm_store_ss(&C[(l_n*88)+82], c82_5);
#else
    C[(l_n*88)+2] += A[323] * B[(l_n*88)+82];
    C[(l_n*88)+8] += A[324] * B[(l_n*88)+82];
    C[(l_n*88)+18] += A[325] * B[(l_n*88)+82];
    C[(l_n*88)+33] += A[326] * B[(l_n*88)+82];
    C[(l_n*88)+54] += A[327] * B[(l_n*88)+82];
    C[(l_n*88)+82] += A[328] * B[(l_n*88)+82];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b83 = _mm_broadcast_ss(&B[(l_n*88)+83]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b83 = _mm_load_ss(&B[(l_n*88)+83]);    b83 = _mm_shuffle_ps(b83, b83, 0x00);
#endif
    __m128 c83_0 = _mm_load_ss(&C[(l_n*88)+0]);
    __m128 a83_0 = _mm_load_ss(&A[329]);
    c83_0 = _mm_add_ss(c83_0, _mm_mul_ss(a83_0, b83));
    _mm_store_ss(&C[(l_n*88)+0], c83_0);
    __m128 c83_1 = _mm_load_ss(&C[(l_n*88)+3]);
    __m128 a83_1 = _mm_load_ss(&A[330]);
    c83_1 = _mm_add_ss(c83_1, _mm_mul_ss(a83_1, b83));
    _mm_store_ss(&C[(l_n*88)+3], c83_1);
    __m128 c83_2 = _mm_load_ss(&C[(l_n*88)+9]);
    __m128 a83_2 = _mm_load_ss(&A[331]);
    c83_2 = _mm_add_ss(c83_2, _mm_mul_ss(a83_2, b83));
    _mm_store_ss(&C[(l_n*88)+9], c83_2);
    __m128 c83_3 = _mm_load_ss(&C[(l_n*88)+19]);
    __m128 a83_3 = _mm_load_ss(&A[332]);
    c83_3 = _mm_add_ss(c83_3, _mm_mul_ss(a83_3, b83));
    _mm_store_ss(&C[(l_n*88)+19], c83_3);
    __m128 c83_4 = _mm_load_ss(&C[(l_n*88)+34]);
    __m128 a83_4 = _mm_load_ss(&A[333]);
    c83_4 = _mm_add_ss(c83_4, _mm_mul_ss(a83_4, b83));
    _mm_store_ss(&C[(l_n*88)+34], c83_4);
    __m128 c83_5 = _mm_load_ss(&C[(l_n*88)+55]);
    __m128 a83_5 = _mm_load_ss(&A[334]);
    c83_5 = _mm_add_ss(c83_5, _mm_mul_ss(a83_5, b83));
    _mm_store_ss(&C[(l_n*88)+55], c83_5);
    __m128 c83_6 = _mm_load_ss(&C[(l_n*88)+83]);
    __m128 a83_6 = _mm_load_ss(&A[335]);
    c83_6 = _mm_add_ss(c83_6, _mm_mul_ss(a83_6, b83));
    _mm_store_ss(&C[(l_n*88)+83], c83_6);
#else
    C[(l_n*88)+0] += A[329] * B[(l_n*88)+83];
    C[(l_n*88)+3] += A[330] * B[(l_n*88)+83];
    C[(l_n*88)+9] += A[331] * B[(l_n*88)+83];
    C[(l_n*88)+19] += A[332] * B[(l_n*88)+83];
    C[(l_n*88)+34] += A[333] * B[(l_n*88)+83];
    C[(l_n*88)+55] += A[334] * B[(l_n*88)+83];
    C[(l_n*88)+83] += A[335] * B[(l_n*88)+83];
#endif

  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 6048;
#endif
}

void ssparse_starMatrix_m120_n9_k9_ldA120_ldBna8_ldC120_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 120; l_m++) {
    C[0+l_m] += A[720+l_m] * B[0];
    C[0+l_m] += A[840+l_m] * B[1];
    C[0+l_m] += A[960+l_m] * B[2];
    C[120+l_m] += A[720+l_m] * B[3];
    C[120+l_m] += A[840+l_m] * B[4];
    C[120+l_m] += A[960+l_m] * B[5];
    C[240+l_m] += A[720+l_m] * B[6];
    C[240+l_m] += A[840+l_m] * B[7];
    C[240+l_m] += A[960+l_m] * B[8];
    C[360+l_m] += A[720+l_m] * B[9];
    C[360+l_m] += A[840+l_m] * B[10];
    C[480+l_m] += A[840+l_m] * B[11];
    C[480+l_m] += A[960+l_m] * B[12];
    C[600+l_m] += A[720+l_m] * B[13];
    C[600+l_m] += A[960+l_m] * B[14];
    C[720+l_m] += A[0+l_m] * B[15];
    C[720+l_m] += A[360+l_m] * B[16];
    C[720+l_m] += A[600+l_m] * B[17];
    C[840+l_m] += A[120+l_m] * B[18];
    C[840+l_m] += A[360+l_m] * B[19];
    C[840+l_m] += A[480+l_m] * B[20];
    C[960+l_m] += A[240+l_m] * B[21];
    C[960+l_m] += A[480+l_m] * B[22];
    C[960+l_m] += A[600+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 5760;
#endif
}

#endif
