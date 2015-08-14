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
// @date 2015-08-14 10:02:43.271255
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
#ifndef SPARSEDKNLCPP
#define SPARSEDKNLCPP

#if defined( __SSE3__) || defined(__MIC__)
#include <immintrin.h>
#endif

#include <cstddef>
#ifndef NDEBUG
extern long long libxsmm_num_total_flops;
#endif

void dsparse_kXiDivMT_m1_n9_k4_ldAna2_ldB8_ldC8_beta0_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_n = 0;
  #pragma nounroll_and_jam
  for ( l_n = 0; l_n < 9; l_n++) {
    unsigned int l_m = 0;
    for ( l_m = 0; l_m < 1; l_m++) {
      C[(l_n*8)+l_m] = 0.0;
    }
#if defined(__SSE3__) || defined(__AVX__)
#else
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b1 = _mm256_broadcast_sd(&B[(l_n*8)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b1 = _mm_loaddup_pd(&B[(l_n*8)+1]);
#endif
    __m128d c1_0 = _mm_load_sd(&C[(l_n*8)+0]);
    __m128d a1_0 = _mm_load_sd(&A[0]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, b1));
#endif
    _mm_store_sd(&C[(l_n*8)+0], c1_0);
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

void dsparse_starMatrix_m1_n9_k9_ldA8_ldBna2_ldC8_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m4_n9_k9_ldA8_ldBna3_ldC8_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
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

void dsparse_starMatrix_m1_n9_k9_ldA8_ldBna3_ldC8_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m10_n9_k9_ldA16_ldBna4_ldC16_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
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

void dsparse_starMatrix_m4_n9_k9_ldA8_ldBna4_ldC8_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
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

void dsparse_starMatrix_m1_n9_k9_ldA8_ldBna4_ldC8_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m20_n9_k9_ldA24_ldBna5_ldC24_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
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

void dsparse_starMatrix_m10_n9_k9_ldA16_ldBna5_ldC16_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
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

void dsparse_starMatrix_m4_n9_k9_ldA8_ldBna5_ldC8_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
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

void dsparse_starMatrix_m1_n9_k9_ldA8_ldBna5_ldC8_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m35_n9_k9_ldA40_ldBna6_ldC40_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
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

void dsparse_starMatrix_m20_n9_k9_ldA24_ldBna6_ldC24_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
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

void dsparse_starMatrix_m10_n9_k9_ldA16_ldBna6_ldC16_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
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

void dsparse_starMatrix_m4_n9_k9_ldA8_ldBna6_ldC8_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
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

void dsparse_starMatrix_m1_n9_k9_ldA8_ldBna6_ldC8_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m56_n9_k9_ldA56_ldBna7_ldC56_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
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

void dsparse_starMatrix_m35_n9_k9_ldA40_ldBna7_ldC40_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
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

void dsparse_starMatrix_m20_n9_k9_ldA24_ldBna7_ldC24_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
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

void dsparse_starMatrix_m10_n9_k9_ldA16_ldBna7_ldC16_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
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

void dsparse_starMatrix_m4_n9_k9_ldA8_ldBna7_ldC8_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
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

void dsparse_starMatrix_m1_n9_k9_ldA8_ldBna7_ldC8_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m4_n9_k9_ldA8_ldBna2_ldC8_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
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

void dsparse_starMatrix_m10_n9_k9_ldA16_ldBna3_ldC16_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
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

void dsparse_starMatrix_m20_n9_k9_ldA24_ldBna4_ldC24_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
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

void dsparse_starMatrix_m35_n9_k9_ldA40_ldBna5_ldC40_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
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

void dsparse_starMatrix_m56_n9_k9_ldA56_ldBna6_ldC56_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
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

void dsparse_starMatrix_m84_n9_k9_ldA88_ldBna7_ldC88_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
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

#endif
