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
// @date 2015-09-27 13:24:48.701664
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
#ifndef SPARSESKNLCPP
#define SPARSESKNLCPP

#if defined( __SSE3__) || defined(__MIC__)
#include <immintrin.h>
#endif

#include <cstddef>
#ifndef NDEBUG
extern long long libxsmm_num_total_flops;
#endif

void ssparse_starMatrix_m1_n9_k9_ldA16_ldBna2_ldC16_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  for ( l_m = 0; l_m < 1; l_m++) {
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
libxsmm_num_total_flops += 48;
#endif
}

void ssparse_starMatrix_m4_n9_k9_ldA16_ldBna3_ldC16_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
  #pragma vector aligned
  for ( l_m = 0; l_m < 4; l_m++) {
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
libxsmm_num_total_flops += 192;
#endif
}

void ssparse_starMatrix_m1_n9_k9_ldA16_ldBna3_ldC16_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  for ( l_m = 0; l_m < 1; l_m++) {
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
libxsmm_num_total_flops += 48;
#endif
}

void ssparse_starMatrix_m10_n9_k9_ldA16_ldBna4_ldC16_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
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

void ssparse_starMatrix_m4_n9_k9_ldA16_ldBna4_ldC16_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
  #pragma vector aligned
  for ( l_m = 0; l_m < 4; l_m++) {
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
libxsmm_num_total_flops += 192;
#endif
}

void ssparse_starMatrix_m1_n9_k9_ldA16_ldBna4_ldC16_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  for ( l_m = 0; l_m < 1; l_m++) {
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
libxsmm_num_total_flops += 48;
#endif
}

void ssparse_starMatrix_m20_n9_k9_ldA32_ldBna5_ldC32_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
  #pragma vector aligned
  for ( l_m = 0; l_m < 20; l_m++) {
    C[0+l_m] += A[192+l_m] * B[0];
    C[0+l_m] += A[224+l_m] * B[1];
    C[0+l_m] += A[256+l_m] * B[2];
    C[32+l_m] += A[192+l_m] * B[3];
    C[32+l_m] += A[224+l_m] * B[4];
    C[32+l_m] += A[256+l_m] * B[5];
    C[64+l_m] += A[192+l_m] * B[6];
    C[64+l_m] += A[224+l_m] * B[7];
    C[64+l_m] += A[256+l_m] * B[8];
    C[96+l_m] += A[192+l_m] * B[9];
    C[96+l_m] += A[224+l_m] * B[10];
    C[128+l_m] += A[224+l_m] * B[11];
    C[128+l_m] += A[256+l_m] * B[12];
    C[160+l_m] += A[192+l_m] * B[13];
    C[160+l_m] += A[256+l_m] * B[14];
    C[192+l_m] += A[0+l_m] * B[15];
    C[192+l_m] += A[96+l_m] * B[16];
    C[192+l_m] += A[160+l_m] * B[17];
    C[224+l_m] += A[32+l_m] * B[18];
    C[224+l_m] += A[96+l_m] * B[19];
    C[224+l_m] += A[128+l_m] * B[20];
    C[256+l_m] += A[64+l_m] * B[21];
    C[256+l_m] += A[128+l_m] * B[22];
    C[256+l_m] += A[160+l_m] * B[23];
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

void ssparse_starMatrix_m4_n9_k9_ldA16_ldBna5_ldC16_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
  #pragma vector aligned
  for ( l_m = 0; l_m < 4; l_m++) {
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
libxsmm_num_total_flops += 192;
#endif
}

void ssparse_starMatrix_m1_n9_k9_ldA16_ldBna5_ldC16_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  for ( l_m = 0; l_m < 1; l_m++) {
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
libxsmm_num_total_flops += 48;
#endif
}

void ssparse_starMatrix_m35_n9_k9_ldA48_ldBna6_ldC48_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
  #pragma vector aligned
  for ( l_m = 0; l_m < 35; l_m++) {
    C[0+l_m] += A[288+l_m] * B[0];
    C[0+l_m] += A[336+l_m] * B[1];
    C[0+l_m] += A[384+l_m] * B[2];
    C[48+l_m] += A[288+l_m] * B[3];
    C[48+l_m] += A[336+l_m] * B[4];
    C[48+l_m] += A[384+l_m] * B[5];
    C[96+l_m] += A[288+l_m] * B[6];
    C[96+l_m] += A[336+l_m] * B[7];
    C[96+l_m] += A[384+l_m] * B[8];
    C[144+l_m] += A[288+l_m] * B[9];
    C[144+l_m] += A[336+l_m] * B[10];
    C[192+l_m] += A[336+l_m] * B[11];
    C[192+l_m] += A[384+l_m] * B[12];
    C[240+l_m] += A[288+l_m] * B[13];
    C[240+l_m] += A[384+l_m] * B[14];
    C[288+l_m] += A[0+l_m] * B[15];
    C[288+l_m] += A[144+l_m] * B[16];
    C[288+l_m] += A[240+l_m] * B[17];
    C[336+l_m] += A[48+l_m] * B[18];
    C[336+l_m] += A[144+l_m] * B[19];
    C[336+l_m] += A[192+l_m] * B[20];
    C[384+l_m] += A[96+l_m] * B[21];
    C[384+l_m] += A[192+l_m] * B[22];
    C[384+l_m] += A[240+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 1680;
#endif
}

void ssparse_starMatrix_m20_n9_k9_ldA32_ldBna6_ldC32_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
  #pragma vector aligned
  for ( l_m = 0; l_m < 20; l_m++) {
    C[0+l_m] += A[192+l_m] * B[0];
    C[0+l_m] += A[224+l_m] * B[1];
    C[0+l_m] += A[256+l_m] * B[2];
    C[32+l_m] += A[192+l_m] * B[3];
    C[32+l_m] += A[224+l_m] * B[4];
    C[32+l_m] += A[256+l_m] * B[5];
    C[64+l_m] += A[192+l_m] * B[6];
    C[64+l_m] += A[224+l_m] * B[7];
    C[64+l_m] += A[256+l_m] * B[8];
    C[96+l_m] += A[192+l_m] * B[9];
    C[96+l_m] += A[224+l_m] * B[10];
    C[128+l_m] += A[224+l_m] * B[11];
    C[128+l_m] += A[256+l_m] * B[12];
    C[160+l_m] += A[192+l_m] * B[13];
    C[160+l_m] += A[256+l_m] * B[14];
    C[192+l_m] += A[0+l_m] * B[15];
    C[192+l_m] += A[96+l_m] * B[16];
    C[192+l_m] += A[160+l_m] * B[17];
    C[224+l_m] += A[32+l_m] * B[18];
    C[224+l_m] += A[96+l_m] * B[19];
    C[224+l_m] += A[128+l_m] * B[20];
    C[256+l_m] += A[64+l_m] * B[21];
    C[256+l_m] += A[128+l_m] * B[22];
    C[256+l_m] += A[160+l_m] * B[23];
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

void ssparse_starMatrix_m4_n9_k9_ldA16_ldBna6_ldC16_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
  #pragma vector aligned
  for ( l_m = 0; l_m < 4; l_m++) {
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
libxsmm_num_total_flops += 192;
#endif
}

void ssparse_starMatrix_m1_n9_k9_ldA16_ldBna6_ldC16_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  for ( l_m = 0; l_m < 1; l_m++) {
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
libxsmm_num_total_flops += 48;
#endif
}

void ssparse_starMatrix_m56_n9_k9_ldA64_ldBna7_ldC64_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
  #pragma vector aligned
  for ( l_m = 0; l_m < 56; l_m++) {
    C[0+l_m] += A[384+l_m] * B[0];
    C[0+l_m] += A[448+l_m] * B[1];
    C[0+l_m] += A[512+l_m] * B[2];
    C[64+l_m] += A[384+l_m] * B[3];
    C[64+l_m] += A[448+l_m] * B[4];
    C[64+l_m] += A[512+l_m] * B[5];
    C[128+l_m] += A[384+l_m] * B[6];
    C[128+l_m] += A[448+l_m] * B[7];
    C[128+l_m] += A[512+l_m] * B[8];
    C[192+l_m] += A[384+l_m] * B[9];
    C[192+l_m] += A[448+l_m] * B[10];
    C[256+l_m] += A[448+l_m] * B[11];
    C[256+l_m] += A[512+l_m] * B[12];
    C[320+l_m] += A[384+l_m] * B[13];
    C[320+l_m] += A[512+l_m] * B[14];
    C[384+l_m] += A[0+l_m] * B[15];
    C[384+l_m] += A[192+l_m] * B[16];
    C[384+l_m] += A[320+l_m] * B[17];
    C[448+l_m] += A[64+l_m] * B[18];
    C[448+l_m] += A[192+l_m] * B[19];
    C[448+l_m] += A[256+l_m] * B[20];
    C[512+l_m] += A[128+l_m] * B[21];
    C[512+l_m] += A[256+l_m] * B[22];
    C[512+l_m] += A[320+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 2688;
#endif
}

void ssparse_starMatrix_m35_n9_k9_ldA48_ldBna7_ldC48_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
  #pragma vector aligned
  for ( l_m = 0; l_m < 35; l_m++) {
    C[0+l_m] += A[288+l_m] * B[0];
    C[0+l_m] += A[336+l_m] * B[1];
    C[0+l_m] += A[384+l_m] * B[2];
    C[48+l_m] += A[288+l_m] * B[3];
    C[48+l_m] += A[336+l_m] * B[4];
    C[48+l_m] += A[384+l_m] * B[5];
    C[96+l_m] += A[288+l_m] * B[6];
    C[96+l_m] += A[336+l_m] * B[7];
    C[96+l_m] += A[384+l_m] * B[8];
    C[144+l_m] += A[288+l_m] * B[9];
    C[144+l_m] += A[336+l_m] * B[10];
    C[192+l_m] += A[336+l_m] * B[11];
    C[192+l_m] += A[384+l_m] * B[12];
    C[240+l_m] += A[288+l_m] * B[13];
    C[240+l_m] += A[384+l_m] * B[14];
    C[288+l_m] += A[0+l_m] * B[15];
    C[288+l_m] += A[144+l_m] * B[16];
    C[288+l_m] += A[240+l_m] * B[17];
    C[336+l_m] += A[48+l_m] * B[18];
    C[336+l_m] += A[144+l_m] * B[19];
    C[336+l_m] += A[192+l_m] * B[20];
    C[384+l_m] += A[96+l_m] * B[21];
    C[384+l_m] += A[192+l_m] * B[22];
    C[384+l_m] += A[240+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 1680;
#endif
}

void ssparse_starMatrix_m20_n9_k9_ldA32_ldBna7_ldC32_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
  #pragma vector aligned
  for ( l_m = 0; l_m < 20; l_m++) {
    C[0+l_m] += A[192+l_m] * B[0];
    C[0+l_m] += A[224+l_m] * B[1];
    C[0+l_m] += A[256+l_m] * B[2];
    C[32+l_m] += A[192+l_m] * B[3];
    C[32+l_m] += A[224+l_m] * B[4];
    C[32+l_m] += A[256+l_m] * B[5];
    C[64+l_m] += A[192+l_m] * B[6];
    C[64+l_m] += A[224+l_m] * B[7];
    C[64+l_m] += A[256+l_m] * B[8];
    C[96+l_m] += A[192+l_m] * B[9];
    C[96+l_m] += A[224+l_m] * B[10];
    C[128+l_m] += A[224+l_m] * B[11];
    C[128+l_m] += A[256+l_m] * B[12];
    C[160+l_m] += A[192+l_m] * B[13];
    C[160+l_m] += A[256+l_m] * B[14];
    C[192+l_m] += A[0+l_m] * B[15];
    C[192+l_m] += A[96+l_m] * B[16];
    C[192+l_m] += A[160+l_m] * B[17];
    C[224+l_m] += A[32+l_m] * B[18];
    C[224+l_m] += A[96+l_m] * B[19];
    C[224+l_m] += A[128+l_m] * B[20];
    C[256+l_m] += A[64+l_m] * B[21];
    C[256+l_m] += A[128+l_m] * B[22];
    C[256+l_m] += A[160+l_m] * B[23];
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

void ssparse_starMatrix_m4_n9_k9_ldA16_ldBna7_ldC16_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
  #pragma vector aligned
  for ( l_m = 0; l_m < 4; l_m++) {
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
libxsmm_num_total_flops += 192;
#endif
}

void ssparse_starMatrix_m1_n9_k9_ldA16_ldBna7_ldC16_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  for ( l_m = 0; l_m < 1; l_m++) {
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
libxsmm_num_total_flops += 48;
#endif
}

void ssparse_starMatrix_m84_n9_k9_ldA96_ldBna8_ldC96_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
  #pragma vector aligned
  for ( l_m = 0; l_m < 84; l_m++) {
    C[0+l_m] += A[576+l_m] * B[0];
    C[0+l_m] += A[672+l_m] * B[1];
    C[0+l_m] += A[768+l_m] * B[2];
    C[96+l_m] += A[576+l_m] * B[3];
    C[96+l_m] += A[672+l_m] * B[4];
    C[96+l_m] += A[768+l_m] * B[5];
    C[192+l_m] += A[576+l_m] * B[6];
    C[192+l_m] += A[672+l_m] * B[7];
    C[192+l_m] += A[768+l_m] * B[8];
    C[288+l_m] += A[576+l_m] * B[9];
    C[288+l_m] += A[672+l_m] * B[10];
    C[384+l_m] += A[672+l_m] * B[11];
    C[384+l_m] += A[768+l_m] * B[12];
    C[480+l_m] += A[576+l_m] * B[13];
    C[480+l_m] += A[768+l_m] * B[14];
    C[576+l_m] += A[0+l_m] * B[15];
    C[576+l_m] += A[288+l_m] * B[16];
    C[576+l_m] += A[480+l_m] * B[17];
    C[672+l_m] += A[96+l_m] * B[18];
    C[672+l_m] += A[288+l_m] * B[19];
    C[672+l_m] += A[384+l_m] * B[20];
    C[768+l_m] += A[192+l_m] * B[21];
    C[768+l_m] += A[384+l_m] * B[22];
    C[768+l_m] += A[480+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 4032;
#endif
}

void ssparse_starMatrix_m56_n9_k9_ldA64_ldBna8_ldC64_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
  #pragma vector aligned
  for ( l_m = 0; l_m < 56; l_m++) {
    C[0+l_m] += A[384+l_m] * B[0];
    C[0+l_m] += A[448+l_m] * B[1];
    C[0+l_m] += A[512+l_m] * B[2];
    C[64+l_m] += A[384+l_m] * B[3];
    C[64+l_m] += A[448+l_m] * B[4];
    C[64+l_m] += A[512+l_m] * B[5];
    C[128+l_m] += A[384+l_m] * B[6];
    C[128+l_m] += A[448+l_m] * B[7];
    C[128+l_m] += A[512+l_m] * B[8];
    C[192+l_m] += A[384+l_m] * B[9];
    C[192+l_m] += A[448+l_m] * B[10];
    C[256+l_m] += A[448+l_m] * B[11];
    C[256+l_m] += A[512+l_m] * B[12];
    C[320+l_m] += A[384+l_m] * B[13];
    C[320+l_m] += A[512+l_m] * B[14];
    C[384+l_m] += A[0+l_m] * B[15];
    C[384+l_m] += A[192+l_m] * B[16];
    C[384+l_m] += A[320+l_m] * B[17];
    C[448+l_m] += A[64+l_m] * B[18];
    C[448+l_m] += A[192+l_m] * B[19];
    C[448+l_m] += A[256+l_m] * B[20];
    C[512+l_m] += A[128+l_m] * B[21];
    C[512+l_m] += A[256+l_m] * B[22];
    C[512+l_m] += A[320+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 2688;
#endif
}

void ssparse_starMatrix_m35_n9_k9_ldA48_ldBna8_ldC48_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
  #pragma vector aligned
  for ( l_m = 0; l_m < 35; l_m++) {
    C[0+l_m] += A[288+l_m] * B[0];
    C[0+l_m] += A[336+l_m] * B[1];
    C[0+l_m] += A[384+l_m] * B[2];
    C[48+l_m] += A[288+l_m] * B[3];
    C[48+l_m] += A[336+l_m] * B[4];
    C[48+l_m] += A[384+l_m] * B[5];
    C[96+l_m] += A[288+l_m] * B[6];
    C[96+l_m] += A[336+l_m] * B[7];
    C[96+l_m] += A[384+l_m] * B[8];
    C[144+l_m] += A[288+l_m] * B[9];
    C[144+l_m] += A[336+l_m] * B[10];
    C[192+l_m] += A[336+l_m] * B[11];
    C[192+l_m] += A[384+l_m] * B[12];
    C[240+l_m] += A[288+l_m] * B[13];
    C[240+l_m] += A[384+l_m] * B[14];
    C[288+l_m] += A[0+l_m] * B[15];
    C[288+l_m] += A[144+l_m] * B[16];
    C[288+l_m] += A[240+l_m] * B[17];
    C[336+l_m] += A[48+l_m] * B[18];
    C[336+l_m] += A[144+l_m] * B[19];
    C[336+l_m] += A[192+l_m] * B[20];
    C[384+l_m] += A[96+l_m] * B[21];
    C[384+l_m] += A[192+l_m] * B[22];
    C[384+l_m] += A[240+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 1680;
#endif
}

void ssparse_starMatrix_m20_n9_k9_ldA32_ldBna8_ldC32_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
  #pragma vector aligned
  for ( l_m = 0; l_m < 20; l_m++) {
    C[0+l_m] += A[192+l_m] * B[0];
    C[0+l_m] += A[224+l_m] * B[1];
    C[0+l_m] += A[256+l_m] * B[2];
    C[32+l_m] += A[192+l_m] * B[3];
    C[32+l_m] += A[224+l_m] * B[4];
    C[32+l_m] += A[256+l_m] * B[5];
    C[64+l_m] += A[192+l_m] * B[6];
    C[64+l_m] += A[224+l_m] * B[7];
    C[64+l_m] += A[256+l_m] * B[8];
    C[96+l_m] += A[192+l_m] * B[9];
    C[96+l_m] += A[224+l_m] * B[10];
    C[128+l_m] += A[224+l_m] * B[11];
    C[128+l_m] += A[256+l_m] * B[12];
    C[160+l_m] += A[192+l_m] * B[13];
    C[160+l_m] += A[256+l_m] * B[14];
    C[192+l_m] += A[0+l_m] * B[15];
    C[192+l_m] += A[96+l_m] * B[16];
    C[192+l_m] += A[160+l_m] * B[17];
    C[224+l_m] += A[32+l_m] * B[18];
    C[224+l_m] += A[96+l_m] * B[19];
    C[224+l_m] += A[128+l_m] * B[20];
    C[256+l_m] += A[64+l_m] * B[21];
    C[256+l_m] += A[128+l_m] * B[22];
    C[256+l_m] += A[160+l_m] * B[23];
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

void ssparse_starMatrix_m4_n9_k9_ldA16_ldBna8_ldC16_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
  #pragma vector aligned
  for ( l_m = 0; l_m < 4; l_m++) {
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
libxsmm_num_total_flops += 192;
#endif
}

void ssparse_starMatrix_m1_n9_k9_ldA16_ldBna8_ldC16_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  for ( l_m = 0; l_m < 1; l_m++) {
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
libxsmm_num_total_flops += 48;
#endif
}

void ssparse_starMatrix_m4_n9_k9_ldA16_ldBna2_ldC16_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
  #pragma vector aligned
  for ( l_m = 0; l_m < 4; l_m++) {
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
libxsmm_num_total_flops += 192;
#endif
}

void ssparse_starMatrix_m10_n9_k9_ldA16_ldBna3_ldC16_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
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

void ssparse_starMatrix_m20_n9_k9_ldA32_ldBna4_ldC32_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
  #pragma vector aligned
  for ( l_m = 0; l_m < 20; l_m++) {
    C[0+l_m] += A[192+l_m] * B[0];
    C[0+l_m] += A[224+l_m] * B[1];
    C[0+l_m] += A[256+l_m] * B[2];
    C[32+l_m] += A[192+l_m] * B[3];
    C[32+l_m] += A[224+l_m] * B[4];
    C[32+l_m] += A[256+l_m] * B[5];
    C[64+l_m] += A[192+l_m] * B[6];
    C[64+l_m] += A[224+l_m] * B[7];
    C[64+l_m] += A[256+l_m] * B[8];
    C[96+l_m] += A[192+l_m] * B[9];
    C[96+l_m] += A[224+l_m] * B[10];
    C[128+l_m] += A[224+l_m] * B[11];
    C[128+l_m] += A[256+l_m] * B[12];
    C[160+l_m] += A[192+l_m] * B[13];
    C[160+l_m] += A[256+l_m] * B[14];
    C[192+l_m] += A[0+l_m] * B[15];
    C[192+l_m] += A[96+l_m] * B[16];
    C[192+l_m] += A[160+l_m] * B[17];
    C[224+l_m] += A[32+l_m] * B[18];
    C[224+l_m] += A[96+l_m] * B[19];
    C[224+l_m] += A[128+l_m] * B[20];
    C[256+l_m] += A[64+l_m] * B[21];
    C[256+l_m] += A[128+l_m] * B[22];
    C[256+l_m] += A[160+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 960;
#endif
}

void ssparse_starMatrix_m35_n9_k9_ldA48_ldBna5_ldC48_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
  #pragma vector aligned
  for ( l_m = 0; l_m < 35; l_m++) {
    C[0+l_m] += A[288+l_m] * B[0];
    C[0+l_m] += A[336+l_m] * B[1];
    C[0+l_m] += A[384+l_m] * B[2];
    C[48+l_m] += A[288+l_m] * B[3];
    C[48+l_m] += A[336+l_m] * B[4];
    C[48+l_m] += A[384+l_m] * B[5];
    C[96+l_m] += A[288+l_m] * B[6];
    C[96+l_m] += A[336+l_m] * B[7];
    C[96+l_m] += A[384+l_m] * B[8];
    C[144+l_m] += A[288+l_m] * B[9];
    C[144+l_m] += A[336+l_m] * B[10];
    C[192+l_m] += A[336+l_m] * B[11];
    C[192+l_m] += A[384+l_m] * B[12];
    C[240+l_m] += A[288+l_m] * B[13];
    C[240+l_m] += A[384+l_m] * B[14];
    C[288+l_m] += A[0+l_m] * B[15];
    C[288+l_m] += A[144+l_m] * B[16];
    C[288+l_m] += A[240+l_m] * B[17];
    C[336+l_m] += A[48+l_m] * B[18];
    C[336+l_m] += A[144+l_m] * B[19];
    C[336+l_m] += A[192+l_m] * B[20];
    C[384+l_m] += A[96+l_m] * B[21];
    C[384+l_m] += A[192+l_m] * B[22];
    C[384+l_m] += A[240+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 1680;
#endif
}

void ssparse_starMatrix_m56_n9_k9_ldA64_ldBna6_ldC64_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
  #pragma vector aligned
  for ( l_m = 0; l_m < 56; l_m++) {
    C[0+l_m] += A[384+l_m] * B[0];
    C[0+l_m] += A[448+l_m] * B[1];
    C[0+l_m] += A[512+l_m] * B[2];
    C[64+l_m] += A[384+l_m] * B[3];
    C[64+l_m] += A[448+l_m] * B[4];
    C[64+l_m] += A[512+l_m] * B[5];
    C[128+l_m] += A[384+l_m] * B[6];
    C[128+l_m] += A[448+l_m] * B[7];
    C[128+l_m] += A[512+l_m] * B[8];
    C[192+l_m] += A[384+l_m] * B[9];
    C[192+l_m] += A[448+l_m] * B[10];
    C[256+l_m] += A[448+l_m] * B[11];
    C[256+l_m] += A[512+l_m] * B[12];
    C[320+l_m] += A[384+l_m] * B[13];
    C[320+l_m] += A[512+l_m] * B[14];
    C[384+l_m] += A[0+l_m] * B[15];
    C[384+l_m] += A[192+l_m] * B[16];
    C[384+l_m] += A[320+l_m] * B[17];
    C[448+l_m] += A[64+l_m] * B[18];
    C[448+l_m] += A[192+l_m] * B[19];
    C[448+l_m] += A[256+l_m] * B[20];
    C[512+l_m] += A[128+l_m] * B[21];
    C[512+l_m] += A[256+l_m] * B[22];
    C[512+l_m] += A[320+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 2688;
#endif
}

void ssparse_starMatrix_m84_n9_k9_ldA96_ldBna7_ldC96_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
  #pragma vector aligned
  for ( l_m = 0; l_m < 84; l_m++) {
    C[0+l_m] += A[576+l_m] * B[0];
    C[0+l_m] += A[672+l_m] * B[1];
    C[0+l_m] += A[768+l_m] * B[2];
    C[96+l_m] += A[576+l_m] * B[3];
    C[96+l_m] += A[672+l_m] * B[4];
    C[96+l_m] += A[768+l_m] * B[5];
    C[192+l_m] += A[576+l_m] * B[6];
    C[192+l_m] += A[672+l_m] * B[7];
    C[192+l_m] += A[768+l_m] * B[8];
    C[288+l_m] += A[576+l_m] * B[9];
    C[288+l_m] += A[672+l_m] * B[10];
    C[384+l_m] += A[672+l_m] * B[11];
    C[384+l_m] += A[768+l_m] * B[12];
    C[480+l_m] += A[576+l_m] * B[13];
    C[480+l_m] += A[768+l_m] * B[14];
    C[576+l_m] += A[0+l_m] * B[15];
    C[576+l_m] += A[288+l_m] * B[16];
    C[576+l_m] += A[480+l_m] * B[17];
    C[672+l_m] += A[96+l_m] * B[18];
    C[672+l_m] += A[288+l_m] * B[19];
    C[672+l_m] += A[384+l_m] * B[20];
    C[768+l_m] += A[192+l_m] * B[21];
    C[768+l_m] += A[384+l_m] * B[22];
    C[768+l_m] += A[480+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 4032;
#endif
}

void ssparse_starMatrix_m120_n9_k9_ldA128_ldBna8_ldC128_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(32)
  #pragma vector aligned
  for ( l_m = 0; l_m < 120; l_m++) {
    C[0+l_m] += A[768+l_m] * B[0];
    C[0+l_m] += A[896+l_m] * B[1];
    C[0+l_m] += A[1024+l_m] * B[2];
    C[128+l_m] += A[768+l_m] * B[3];
    C[128+l_m] += A[896+l_m] * B[4];
    C[128+l_m] += A[1024+l_m] * B[5];
    C[256+l_m] += A[768+l_m] * B[6];
    C[256+l_m] += A[896+l_m] * B[7];
    C[256+l_m] += A[1024+l_m] * B[8];
    C[384+l_m] += A[768+l_m] * B[9];
    C[384+l_m] += A[896+l_m] * B[10];
    C[512+l_m] += A[896+l_m] * B[11];
    C[512+l_m] += A[1024+l_m] * B[12];
    C[640+l_m] += A[768+l_m] * B[13];
    C[640+l_m] += A[1024+l_m] * B[14];
    C[768+l_m] += A[0+l_m] * B[15];
    C[768+l_m] += A[384+l_m] * B[16];
    C[768+l_m] += A[640+l_m] * B[17];
    C[896+l_m] += A[128+l_m] * B[18];
    C[896+l_m] += A[384+l_m] * B[19];
    C[896+l_m] += A[512+l_m] * B[20];
    C[1024+l_m] += A[256+l_m] * B[21];
    C[1024+l_m] += A[512+l_m] * B[22];
    C[1024+l_m] += A[640+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 5760;
#endif
}

#endif
