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
// @date 2015-09-27 13:24:48.675790
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
#ifndef SPARSESSKXCPP
#define SPARSESSKXCPP

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

void ssparse_fP113DivM_m4_n9_k4_ldAna2_ldB16_ldC16_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_n = 0;
  #pragma nounroll_and_jam
  for ( l_n = 0; l_n < 9; l_n++) {
    unsigned int l_m = 0;
   #pragma simd
    for ( l_m = 0; l_m < 4; l_m++) {
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
#else
    C[(l_n*16)+0] += A[0] * B[(l_n*16)+0];
    C[(l_n*16)+3] += A[1] * B[(l_n*16)+0];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b1 = _mm_broadcast_ss(&B[(l_n*16)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b1 = _mm_load_ss(&B[(l_n*16)+1]);    b1 = _mm_shuffle_ps(b1, b1, 0x00);
#endif
    __m128 c1_0 = _mm_load_ss(&C[(l_n*16)+1]);
    __m128 a1_0 = _mm_load_ss(&A[2]);
    c1_0 = _mm_add_ss(c1_0, _mm_mul_ss(a1_0, b1));
    _mm_store_ss(&C[(l_n*16)+1], c1_0);
#else
    C[(l_n*16)+1] += A[2] * B[(l_n*16)+1];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b2 = _mm_broadcast_ss(&B[(l_n*16)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b2 = _mm_load_ss(&B[(l_n*16)+2]);    b2 = _mm_shuffle_ps(b2, b2, 0x00);
#endif
    __m128 c2_0 = _mm_load_ss(&C[(l_n*16)+2]);
    __m128 a2_0 = _mm_load_ss(&A[3]);
    c2_0 = _mm_add_ss(c2_0, _mm_mul_ss(a2_0, b2));
    _mm_store_ss(&C[(l_n*16)+2], c2_0);
#else
    C[(l_n*16)+2] += A[3] * B[(l_n*16)+2];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b3 = _mm_broadcast_ss(&B[(l_n*16)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b3 = _mm_load_ss(&B[(l_n*16)+3]);    b3 = _mm_shuffle_ps(b3, b3, 0x00);
#endif
    __m128 c3_0 = _mm_load_ss(&C[(l_n*16)+0]);
    __m128 a3_0 = _mm_load_ss(&A[4]);
    c3_0 = _mm_add_ss(c3_0, _mm_mul_ss(a3_0, b3));
    _mm_store_ss(&C[(l_n*16)+0], c3_0);
    __m128 c3_1 = _mm_load_ss(&C[(l_n*16)+3]);
    __m128 a3_1 = _mm_load_ss(&A[5]);
    c3_1 = _mm_add_ss(c3_1, _mm_mul_ss(a3_1, b3));
    _mm_store_ss(&C[(l_n*16)+3], c3_1);
#else
    C[(l_n*16)+0] += A[4] * B[(l_n*16)+3];
    C[(l_n*16)+3] += A[5] * B[(l_n*16)+3];
#endif

  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 108;
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

void ssparse_fP113DivM_m56_n9_k56_ldAna6_ldB64_ldC64_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_n = 0;
  #pragma nounroll_and_jam
  for ( l_n = 0; l_n < 9; l_n++) {
    unsigned int l_m = 0;
   #pragma simd
    for ( l_m = 0; l_m < 56; l_m++) {
      C[(l_n*64)+l_m] = 0.0f;
    }
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b0 = _mm_broadcast_ss(&B[(l_n*64)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b0 = _mm_load_ss(&B[(l_n*64)+0]);    b0 = _mm_shuffle_ps(b0, b0, 0x00);
#endif
    __m128 c0_0 = _mm_load_ss(&C[(l_n*64)+0]);
    __m128 a0_0 = _mm_load_ss(&A[0]);
    c0_0 = _mm_add_ss(c0_0, _mm_mul_ss(a0_0, b0));
    _mm_store_ss(&C[(l_n*64)+0], c0_0);
    __m128 c0_1 = _mm_load_ss(&C[(l_n*64)+3]);
    __m128 a0_1 = _mm_load_ss(&A[1]);
    c0_1 = _mm_add_ss(c0_1, _mm_mul_ss(a0_1, b0));
    _mm_store_ss(&C[(l_n*64)+3], c0_1);
    __m128 c0_2 = _mm_load_ss(&C[(l_n*64)+9]);
    __m128 a0_2 = _mm_load_ss(&A[2]);
    c0_2 = _mm_add_ss(c0_2, _mm_mul_ss(a0_2, b0));
    _mm_store_ss(&C[(l_n*64)+9], c0_2);
    __m128 c0_3 = _mm_load_ss(&C[(l_n*64)+19]);
    __m128 a0_3 = _mm_load_ss(&A[3]);
    c0_3 = _mm_add_ss(c0_3, _mm_mul_ss(a0_3, b0));
    _mm_store_ss(&C[(l_n*64)+19], c0_3);
    __m128 c0_4 = _mm_load_ss(&C[(l_n*64)+34]);
    __m128 a0_4 = _mm_load_ss(&A[4]);
    c0_4 = _mm_add_ss(c0_4, _mm_mul_ss(a0_4, b0));
    _mm_store_ss(&C[(l_n*64)+34], c0_4);
    __m128 c0_5 = _mm_load_ss(&C[(l_n*64)+55]);
    __m128 a0_5 = _mm_load_ss(&A[5]);
    c0_5 = _mm_add_ss(c0_5, _mm_mul_ss(a0_5, b0));
    _mm_store_ss(&C[(l_n*64)+55], c0_5);
#else
    C[(l_n*64)+0] += A[0] * B[(l_n*64)+0];
    C[(l_n*64)+3] += A[1] * B[(l_n*64)+0];
    C[(l_n*64)+9] += A[2] * B[(l_n*64)+0];
    C[(l_n*64)+19] += A[3] * B[(l_n*64)+0];
    C[(l_n*64)+34] += A[4] * B[(l_n*64)+0];
    C[(l_n*64)+55] += A[5] * B[(l_n*64)+0];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b1 = _mm_broadcast_ss(&B[(l_n*64)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b1 = _mm_load_ss(&B[(l_n*64)+1]);    b1 = _mm_shuffle_ps(b1, b1, 0x00);
#endif
    __m128 c1_0 = _mm_load_ss(&C[(l_n*64)+1]);
    __m128 a1_0 = _mm_load_ss(&A[6]);
    c1_0 = _mm_add_ss(c1_0, _mm_mul_ss(a1_0, b1));
    _mm_store_ss(&C[(l_n*64)+1], c1_0);
    __m128 c1_1 = _mm_load_ss(&C[(l_n*64)+7]);
    __m128 a1_1 = _mm_load_ss(&A[7]);
    c1_1 = _mm_add_ss(c1_1, _mm_mul_ss(a1_1, b1));
    _mm_store_ss(&C[(l_n*64)+7], c1_1);
    __m128 c1_2 = _mm_load_ss(&C[(l_n*64)+17]);
    __m128 a1_2 = _mm_load_ss(&A[8]);
    c1_2 = _mm_add_ss(c1_2, _mm_mul_ss(a1_2, b1));
    _mm_store_ss(&C[(l_n*64)+17], c1_2);
    __m128 c1_3 = _mm_load_ss(&C[(l_n*64)+32]);
    __m128 a1_3 = _mm_load_ss(&A[9]);
    c1_3 = _mm_add_ss(c1_3, _mm_mul_ss(a1_3, b1));
    _mm_store_ss(&C[(l_n*64)+32], c1_3);
    __m128 c1_4 = _mm_load_ss(&C[(l_n*64)+53]);
    __m128 a1_4 = _mm_load_ss(&A[10]);
    c1_4 = _mm_add_ss(c1_4, _mm_mul_ss(a1_4, b1));
    _mm_store_ss(&C[(l_n*64)+53], c1_4);
#else
    C[(l_n*64)+1] += A[6] * B[(l_n*64)+1];
    C[(l_n*64)+7] += A[7] * B[(l_n*64)+1];
    C[(l_n*64)+17] += A[8] * B[(l_n*64)+1];
    C[(l_n*64)+32] += A[9] * B[(l_n*64)+1];
    C[(l_n*64)+53] += A[10] * B[(l_n*64)+1];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b2 = _mm_broadcast_ss(&B[(l_n*64)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b2 = _mm_load_ss(&B[(l_n*64)+2]);    b2 = _mm_shuffle_ps(b2, b2, 0x00);
#endif
    __m128 c2_0 = _mm_load_ss(&C[(l_n*64)+2]);
    __m128 a2_0 = _mm_load_ss(&A[11]);
    c2_0 = _mm_add_ss(c2_0, _mm_mul_ss(a2_0, b2));
    _mm_store_ss(&C[(l_n*64)+2], c2_0);
    __m128 c2_1 = _mm_load_ss(&C[(l_n*64)+8]);
    __m128 a2_1 = _mm_load_ss(&A[12]);
    c2_1 = _mm_add_ss(c2_1, _mm_mul_ss(a2_1, b2));
    _mm_store_ss(&C[(l_n*64)+8], c2_1);
    __m128 c2_2 = _mm_load_ss(&C[(l_n*64)+18]);
    __m128 a2_2 = _mm_load_ss(&A[13]);
    c2_2 = _mm_add_ss(c2_2, _mm_mul_ss(a2_2, b2));
    _mm_store_ss(&C[(l_n*64)+18], c2_2);
    __m128 c2_3 = _mm_load_ss(&C[(l_n*64)+33]);
    __m128 a2_3 = _mm_load_ss(&A[14]);
    c2_3 = _mm_add_ss(c2_3, _mm_mul_ss(a2_3, b2));
    _mm_store_ss(&C[(l_n*64)+33], c2_3);
    __m128 c2_4 = _mm_load_ss(&C[(l_n*64)+54]);
    __m128 a2_4 = _mm_load_ss(&A[15]);
    c2_4 = _mm_add_ss(c2_4, _mm_mul_ss(a2_4, b2));
    _mm_store_ss(&C[(l_n*64)+54], c2_4);
#else
    C[(l_n*64)+2] += A[11] * B[(l_n*64)+2];
    C[(l_n*64)+8] += A[12] * B[(l_n*64)+2];
    C[(l_n*64)+18] += A[13] * B[(l_n*64)+2];
    C[(l_n*64)+33] += A[14] * B[(l_n*64)+2];
    C[(l_n*64)+54] += A[15] * B[(l_n*64)+2];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b3 = _mm_broadcast_ss(&B[(l_n*64)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b3 = _mm_load_ss(&B[(l_n*64)+3]);    b3 = _mm_shuffle_ps(b3, b3, 0x00);
#endif
    __m128 c3_0 = _mm_load_ss(&C[(l_n*64)+0]);
    __m128 a3_0 = _mm_load_ss(&A[16]);
    c3_0 = _mm_add_ss(c3_0, _mm_mul_ss(a3_0, b3));
    _mm_store_ss(&C[(l_n*64)+0], c3_0);
    __m128 c3_1 = _mm_load_ss(&C[(l_n*64)+3]);
    __m128 a3_1 = _mm_load_ss(&A[17]);
    c3_1 = _mm_add_ss(c3_1, _mm_mul_ss(a3_1, b3));
    _mm_store_ss(&C[(l_n*64)+3], c3_1);
    __m128 c3_2 = _mm_load_ss(&C[(l_n*64)+9]);
    __m128 a3_2 = _mm_load_ss(&A[18]);
    c3_2 = _mm_add_ss(c3_2, _mm_mul_ss(a3_2, b3));
    _mm_store_ss(&C[(l_n*64)+9], c3_2);
    __m128 c3_3 = _mm_load_ss(&C[(l_n*64)+19]);
    __m128 a3_3 = _mm_load_ss(&A[19]);
    c3_3 = _mm_add_ss(c3_3, _mm_mul_ss(a3_3, b3));
    _mm_store_ss(&C[(l_n*64)+19], c3_3);
    __m128 c3_4 = _mm_load_ss(&C[(l_n*64)+34]);
    __m128 a3_4 = _mm_load_ss(&A[20]);
    c3_4 = _mm_add_ss(c3_4, _mm_mul_ss(a3_4, b3));
    _mm_store_ss(&C[(l_n*64)+34], c3_4);
    __m128 c3_5 = _mm_load_ss(&C[(l_n*64)+55]);
    __m128 a3_5 = _mm_load_ss(&A[21]);
    c3_5 = _mm_add_ss(c3_5, _mm_mul_ss(a3_5, b3));
    _mm_store_ss(&C[(l_n*64)+55], c3_5);
#else
    C[(l_n*64)+0] += A[16] * B[(l_n*64)+3];
    C[(l_n*64)+3] += A[17] * B[(l_n*64)+3];
    C[(l_n*64)+9] += A[18] * B[(l_n*64)+3];
    C[(l_n*64)+19] += A[19] * B[(l_n*64)+3];
    C[(l_n*64)+34] += A[20] * B[(l_n*64)+3];
    C[(l_n*64)+55] += A[21] * B[(l_n*64)+3];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b4 = _mm_broadcast_ss(&B[(l_n*64)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b4 = _mm_load_ss(&B[(l_n*64)+4]);    b4 = _mm_shuffle_ps(b4, b4, 0x00);
#endif
    __m128 c4_0 = _mm_load_ss(&C[(l_n*64)+4]);
    __m128 a4_0 = _mm_load_ss(&A[22]);
    c4_0 = _mm_add_ss(c4_0, _mm_mul_ss(a4_0, b4));
    _mm_store_ss(&C[(l_n*64)+4], c4_0);
    __m128 c4_1 = _mm_load_ss(&C[(l_n*64)+14]);
    __m128 a4_1 = _mm_load_ss(&A[23]);
    c4_1 = _mm_add_ss(c4_1, _mm_mul_ss(a4_1, b4));
    _mm_store_ss(&C[(l_n*64)+14], c4_1);
    __m128 c4_2 = _mm_load_ss(&C[(l_n*64)+29]);
    __m128 a4_2 = _mm_load_ss(&A[24]);
    c4_2 = _mm_add_ss(c4_2, _mm_mul_ss(a4_2, b4));
    _mm_store_ss(&C[(l_n*64)+29], c4_2);
    __m128 c4_3 = _mm_load_ss(&C[(l_n*64)+50]);
    __m128 a4_3 = _mm_load_ss(&A[25]);
    c4_3 = _mm_add_ss(c4_3, _mm_mul_ss(a4_3, b4));
    _mm_store_ss(&C[(l_n*64)+50], c4_3);
#else
    C[(l_n*64)+4] += A[22] * B[(l_n*64)+4];
    C[(l_n*64)+14] += A[23] * B[(l_n*64)+4];
    C[(l_n*64)+29] += A[24] * B[(l_n*64)+4];
    C[(l_n*64)+50] += A[25] * B[(l_n*64)+4];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b5 = _mm_broadcast_ss(&B[(l_n*64)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b5 = _mm_load_ss(&B[(l_n*64)+5]);    b5 = _mm_shuffle_ps(b5, b5, 0x00);
#endif
    __m128 c5_0 = _mm_load_ss(&C[(l_n*64)+5]);
    __m128 a5_0 = _mm_load_ss(&A[26]);
    c5_0 = _mm_add_ss(c5_0, _mm_mul_ss(a5_0, b5));
    _mm_store_ss(&C[(l_n*64)+5], c5_0);
    __m128 c5_1 = _mm_load_ss(&C[(l_n*64)+15]);
    __m128 a5_1 = _mm_load_ss(&A[27]);
    c5_1 = _mm_add_ss(c5_1, _mm_mul_ss(a5_1, b5));
    _mm_store_ss(&C[(l_n*64)+15], c5_1);
    __m128 c5_2 = _mm_load_ss(&C[(l_n*64)+30]);
    __m128 a5_2 = _mm_load_ss(&A[28]);
    c5_2 = _mm_add_ss(c5_2, _mm_mul_ss(a5_2, b5));
    _mm_store_ss(&C[(l_n*64)+30], c5_2);
    __m128 c5_3 = _mm_load_ss(&C[(l_n*64)+51]);
    __m128 a5_3 = _mm_load_ss(&A[29]);
    c5_3 = _mm_add_ss(c5_3, _mm_mul_ss(a5_3, b5));
    _mm_store_ss(&C[(l_n*64)+51], c5_3);
#else
    C[(l_n*64)+5] += A[26] * B[(l_n*64)+5];
    C[(l_n*64)+15] += A[27] * B[(l_n*64)+5];
    C[(l_n*64)+30] += A[28] * B[(l_n*64)+5];
    C[(l_n*64)+51] += A[29] * B[(l_n*64)+5];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b6 = _mm_broadcast_ss(&B[(l_n*64)+6]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b6 = _mm_load_ss(&B[(l_n*64)+6]);    b6 = _mm_shuffle_ps(b6, b6, 0x00);
#endif
    __m128 c6_0 = _mm_load_ss(&C[(l_n*64)+6]);
    __m128 a6_0 = _mm_load_ss(&A[30]);
    c6_0 = _mm_add_ss(c6_0, _mm_mul_ss(a6_0, b6));
    _mm_store_ss(&C[(l_n*64)+6], c6_0);
    __m128 c6_1 = _mm_load_ss(&C[(l_n*64)+16]);
    __m128 a6_1 = _mm_load_ss(&A[31]);
    c6_1 = _mm_add_ss(c6_1, _mm_mul_ss(a6_1, b6));
    _mm_store_ss(&C[(l_n*64)+16], c6_1);
    __m128 c6_2 = _mm_load_ss(&C[(l_n*64)+31]);
    __m128 a6_2 = _mm_load_ss(&A[32]);
    c6_2 = _mm_add_ss(c6_2, _mm_mul_ss(a6_2, b6));
    _mm_store_ss(&C[(l_n*64)+31], c6_2);
    __m128 c6_3 = _mm_load_ss(&C[(l_n*64)+52]);
    __m128 a6_3 = _mm_load_ss(&A[33]);
    c6_3 = _mm_add_ss(c6_3, _mm_mul_ss(a6_3, b6));
    _mm_store_ss(&C[(l_n*64)+52], c6_3);
#else
    C[(l_n*64)+6] += A[30] * B[(l_n*64)+6];
    C[(l_n*64)+16] += A[31] * B[(l_n*64)+6];
    C[(l_n*64)+31] += A[32] * B[(l_n*64)+6];
    C[(l_n*64)+52] += A[33] * B[(l_n*64)+6];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b7 = _mm_broadcast_ss(&B[(l_n*64)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b7 = _mm_load_ss(&B[(l_n*64)+7]);    b7 = _mm_shuffle_ps(b7, b7, 0x00);
#endif
    __m128 c7_0 = _mm_load_ss(&C[(l_n*64)+1]);
    __m128 a7_0 = _mm_load_ss(&A[34]);
    c7_0 = _mm_add_ss(c7_0, _mm_mul_ss(a7_0, b7));
    _mm_store_ss(&C[(l_n*64)+1], c7_0);
    __m128 c7_1 = _mm_load_ss(&C[(l_n*64)+7]);
    __m128 a7_1 = _mm_load_ss(&A[35]);
    c7_1 = _mm_add_ss(c7_1, _mm_mul_ss(a7_1, b7));
    _mm_store_ss(&C[(l_n*64)+7], c7_1);
    __m128 c7_2 = _mm_load_ss(&C[(l_n*64)+17]);
    __m128 a7_2 = _mm_load_ss(&A[36]);
    c7_2 = _mm_add_ss(c7_2, _mm_mul_ss(a7_2, b7));
    _mm_store_ss(&C[(l_n*64)+17], c7_2);
    __m128 c7_3 = _mm_load_ss(&C[(l_n*64)+32]);
    __m128 a7_3 = _mm_load_ss(&A[37]);
    c7_3 = _mm_add_ss(c7_3, _mm_mul_ss(a7_3, b7));
    _mm_store_ss(&C[(l_n*64)+32], c7_3);
    __m128 c7_4 = _mm_load_ss(&C[(l_n*64)+53]);
    __m128 a7_4 = _mm_load_ss(&A[38]);
    c7_4 = _mm_add_ss(c7_4, _mm_mul_ss(a7_4, b7));
    _mm_store_ss(&C[(l_n*64)+53], c7_4);
#else
    C[(l_n*64)+1] += A[34] * B[(l_n*64)+7];
    C[(l_n*64)+7] += A[35] * B[(l_n*64)+7];
    C[(l_n*64)+17] += A[36] * B[(l_n*64)+7];
    C[(l_n*64)+32] += A[37] * B[(l_n*64)+7];
    C[(l_n*64)+53] += A[38] * B[(l_n*64)+7];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b8 = _mm_broadcast_ss(&B[(l_n*64)+8]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b8 = _mm_load_ss(&B[(l_n*64)+8]);    b8 = _mm_shuffle_ps(b8, b8, 0x00);
#endif
    __m128 c8_0 = _mm_load_ss(&C[(l_n*64)+2]);
    __m128 a8_0 = _mm_load_ss(&A[39]);
    c8_0 = _mm_add_ss(c8_0, _mm_mul_ss(a8_0, b8));
    _mm_store_ss(&C[(l_n*64)+2], c8_0);
    __m128 c8_1 = _mm_load_ss(&C[(l_n*64)+8]);
    __m128 a8_1 = _mm_load_ss(&A[40]);
    c8_1 = _mm_add_ss(c8_1, _mm_mul_ss(a8_1, b8));
    _mm_store_ss(&C[(l_n*64)+8], c8_1);
    __m128 c8_2 = _mm_load_ss(&C[(l_n*64)+18]);
    __m128 a8_2 = _mm_load_ss(&A[41]);
    c8_2 = _mm_add_ss(c8_2, _mm_mul_ss(a8_2, b8));
    _mm_store_ss(&C[(l_n*64)+18], c8_2);
    __m128 c8_3 = _mm_load_ss(&C[(l_n*64)+33]);
    __m128 a8_3 = _mm_load_ss(&A[42]);
    c8_3 = _mm_add_ss(c8_3, _mm_mul_ss(a8_3, b8));
    _mm_store_ss(&C[(l_n*64)+33], c8_3);
    __m128 c8_4 = _mm_load_ss(&C[(l_n*64)+54]);
    __m128 a8_4 = _mm_load_ss(&A[43]);
    c8_4 = _mm_add_ss(c8_4, _mm_mul_ss(a8_4, b8));
    _mm_store_ss(&C[(l_n*64)+54], c8_4);
#else
    C[(l_n*64)+2] += A[39] * B[(l_n*64)+8];
    C[(l_n*64)+8] += A[40] * B[(l_n*64)+8];
    C[(l_n*64)+18] += A[41] * B[(l_n*64)+8];
    C[(l_n*64)+33] += A[42] * B[(l_n*64)+8];
    C[(l_n*64)+54] += A[43] * B[(l_n*64)+8];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b9 = _mm_broadcast_ss(&B[(l_n*64)+9]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b9 = _mm_load_ss(&B[(l_n*64)+9]);    b9 = _mm_shuffle_ps(b9, b9, 0x00);
#endif
    __m128 c9_0 = _mm_load_ss(&C[(l_n*64)+0]);
    __m128 a9_0 = _mm_load_ss(&A[44]);
    c9_0 = _mm_add_ss(c9_0, _mm_mul_ss(a9_0, b9));
    _mm_store_ss(&C[(l_n*64)+0], c9_0);
    __m128 c9_1 = _mm_load_ss(&C[(l_n*64)+3]);
    __m128 a9_1 = _mm_load_ss(&A[45]);
    c9_1 = _mm_add_ss(c9_1, _mm_mul_ss(a9_1, b9));
    _mm_store_ss(&C[(l_n*64)+3], c9_1);
    __m128 c9_2 = _mm_load_ss(&C[(l_n*64)+9]);
    __m128 a9_2 = _mm_load_ss(&A[46]);
    c9_2 = _mm_add_ss(c9_2, _mm_mul_ss(a9_2, b9));
    _mm_store_ss(&C[(l_n*64)+9], c9_2);
    __m128 c9_3 = _mm_load_ss(&C[(l_n*64)+19]);
    __m128 a9_3 = _mm_load_ss(&A[47]);
    c9_3 = _mm_add_ss(c9_3, _mm_mul_ss(a9_3, b9));
    _mm_store_ss(&C[(l_n*64)+19], c9_3);
    __m128 c9_4 = _mm_load_ss(&C[(l_n*64)+34]);
    __m128 a9_4 = _mm_load_ss(&A[48]);
    c9_4 = _mm_add_ss(c9_4, _mm_mul_ss(a9_4, b9));
    _mm_store_ss(&C[(l_n*64)+34], c9_4);
    __m128 c9_5 = _mm_load_ss(&C[(l_n*64)+55]);
    __m128 a9_5 = _mm_load_ss(&A[49]);
    c9_5 = _mm_add_ss(c9_5, _mm_mul_ss(a9_5, b9));
    _mm_store_ss(&C[(l_n*64)+55], c9_5);
#else
    C[(l_n*64)+0] += A[44] * B[(l_n*64)+9];
    C[(l_n*64)+3] += A[45] * B[(l_n*64)+9];
    C[(l_n*64)+9] += A[46] * B[(l_n*64)+9];
    C[(l_n*64)+19] += A[47] * B[(l_n*64)+9];
    C[(l_n*64)+34] += A[48] * B[(l_n*64)+9];
    C[(l_n*64)+55] += A[49] * B[(l_n*64)+9];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b10 = _mm_broadcast_ss(&B[(l_n*64)+10]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b10 = _mm_load_ss(&B[(l_n*64)+10]);    b10 = _mm_shuffle_ps(b10, b10, 0x00);
#endif
    __m128 c10_0 = _mm_load_ss(&C[(l_n*64)+10]);
    __m128 a10_0 = _mm_load_ss(&A[50]);
    c10_0 = _mm_add_ss(c10_0, _mm_mul_ss(a10_0, b10));
    _mm_store_ss(&C[(l_n*64)+10], c10_0);
    __m128 c10_1 = _mm_load_ss(&C[(l_n*64)+25]);
    __m128 a10_1 = _mm_load_ss(&A[51]);
    c10_1 = _mm_add_ss(c10_1, _mm_mul_ss(a10_1, b10));
    _mm_store_ss(&C[(l_n*64)+25], c10_1);
    __m128 c10_2 = _mm_load_ss(&C[(l_n*64)+46]);
    __m128 a10_2 = _mm_load_ss(&A[52]);
    c10_2 = _mm_add_ss(c10_2, _mm_mul_ss(a10_2, b10));
    _mm_store_ss(&C[(l_n*64)+46], c10_2);
#else
    C[(l_n*64)+10] += A[50] * B[(l_n*64)+10];
    C[(l_n*64)+25] += A[51] * B[(l_n*64)+10];
    C[(l_n*64)+46] += A[52] * B[(l_n*64)+10];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b11 = _mm_broadcast_ss(&B[(l_n*64)+11]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b11 = _mm_load_ss(&B[(l_n*64)+11]);    b11 = _mm_shuffle_ps(b11, b11, 0x00);
#endif
    __m128 c11_0 = _mm_load_ss(&C[(l_n*64)+11]);
    __m128 a11_0 = _mm_load_ss(&A[53]);
    c11_0 = _mm_add_ss(c11_0, _mm_mul_ss(a11_0, b11));
    _mm_store_ss(&C[(l_n*64)+11], c11_0);
    __m128 c11_1 = _mm_load_ss(&C[(l_n*64)+26]);
    __m128 a11_1 = _mm_load_ss(&A[54]);
    c11_1 = _mm_add_ss(c11_1, _mm_mul_ss(a11_1, b11));
    _mm_store_ss(&C[(l_n*64)+26], c11_1);
    __m128 c11_2 = _mm_load_ss(&C[(l_n*64)+47]);
    __m128 a11_2 = _mm_load_ss(&A[55]);
    c11_2 = _mm_add_ss(c11_2, _mm_mul_ss(a11_2, b11));
    _mm_store_ss(&C[(l_n*64)+47], c11_2);
#else
    C[(l_n*64)+11] += A[53] * B[(l_n*64)+11];
    C[(l_n*64)+26] += A[54] * B[(l_n*64)+11];
    C[(l_n*64)+47] += A[55] * B[(l_n*64)+11];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b12 = _mm_broadcast_ss(&B[(l_n*64)+12]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b12 = _mm_load_ss(&B[(l_n*64)+12]);    b12 = _mm_shuffle_ps(b12, b12, 0x00);
#endif
    __m128 c12_0 = _mm_load_ss(&C[(l_n*64)+12]);
    __m128 a12_0 = _mm_load_ss(&A[56]);
    c12_0 = _mm_add_ss(c12_0, _mm_mul_ss(a12_0, b12));
    _mm_store_ss(&C[(l_n*64)+12], c12_0);
    __m128 c12_1 = _mm_load_ss(&C[(l_n*64)+27]);
    __m128 a12_1 = _mm_load_ss(&A[57]);
    c12_1 = _mm_add_ss(c12_1, _mm_mul_ss(a12_1, b12));
    _mm_store_ss(&C[(l_n*64)+27], c12_1);
    __m128 c12_2 = _mm_load_ss(&C[(l_n*64)+48]);
    __m128 a12_2 = _mm_load_ss(&A[58]);
    c12_2 = _mm_add_ss(c12_2, _mm_mul_ss(a12_2, b12));
    _mm_store_ss(&C[(l_n*64)+48], c12_2);
#else
    C[(l_n*64)+12] += A[56] * B[(l_n*64)+12];
    C[(l_n*64)+27] += A[57] * B[(l_n*64)+12];
    C[(l_n*64)+48] += A[58] * B[(l_n*64)+12];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b13 = _mm_broadcast_ss(&B[(l_n*64)+13]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b13 = _mm_load_ss(&B[(l_n*64)+13]);    b13 = _mm_shuffle_ps(b13, b13, 0x00);
#endif
    __m128 c13_0 = _mm_load_ss(&C[(l_n*64)+13]);
    __m128 a13_0 = _mm_load_ss(&A[59]);
    c13_0 = _mm_add_ss(c13_0, _mm_mul_ss(a13_0, b13));
    _mm_store_ss(&C[(l_n*64)+13], c13_0);
    __m128 c13_1 = _mm_load_ss(&C[(l_n*64)+28]);
    __m128 a13_1 = _mm_load_ss(&A[60]);
    c13_1 = _mm_add_ss(c13_1, _mm_mul_ss(a13_1, b13));
    _mm_store_ss(&C[(l_n*64)+28], c13_1);
    __m128 c13_2 = _mm_load_ss(&C[(l_n*64)+49]);
    __m128 a13_2 = _mm_load_ss(&A[61]);
    c13_2 = _mm_add_ss(c13_2, _mm_mul_ss(a13_2, b13));
    _mm_store_ss(&C[(l_n*64)+49], c13_2);
#else
    C[(l_n*64)+13] += A[59] * B[(l_n*64)+13];
    C[(l_n*64)+28] += A[60] * B[(l_n*64)+13];
    C[(l_n*64)+49] += A[61] * B[(l_n*64)+13];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b14 = _mm_broadcast_ss(&B[(l_n*64)+14]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b14 = _mm_load_ss(&B[(l_n*64)+14]);    b14 = _mm_shuffle_ps(b14, b14, 0x00);
#endif
    __m128 c14_0 = _mm_load_ss(&C[(l_n*64)+4]);
    __m128 a14_0 = _mm_load_ss(&A[62]);
    c14_0 = _mm_add_ss(c14_0, _mm_mul_ss(a14_0, b14));
    _mm_store_ss(&C[(l_n*64)+4], c14_0);
    __m128 c14_1 = _mm_load_ss(&C[(l_n*64)+14]);
    __m128 a14_1 = _mm_load_ss(&A[63]);
    c14_1 = _mm_add_ss(c14_1, _mm_mul_ss(a14_1, b14));
    _mm_store_ss(&C[(l_n*64)+14], c14_1);
    __m128 c14_2 = _mm_load_ss(&C[(l_n*64)+29]);
    __m128 a14_2 = _mm_load_ss(&A[64]);
    c14_2 = _mm_add_ss(c14_2, _mm_mul_ss(a14_2, b14));
    _mm_store_ss(&C[(l_n*64)+29], c14_2);
    __m128 c14_3 = _mm_load_ss(&C[(l_n*64)+50]);
    __m128 a14_3 = _mm_load_ss(&A[65]);
    c14_3 = _mm_add_ss(c14_3, _mm_mul_ss(a14_3, b14));
    _mm_store_ss(&C[(l_n*64)+50], c14_3);
#else
    C[(l_n*64)+4] += A[62] * B[(l_n*64)+14];
    C[(l_n*64)+14] += A[63] * B[(l_n*64)+14];
    C[(l_n*64)+29] += A[64] * B[(l_n*64)+14];
    C[(l_n*64)+50] += A[65] * B[(l_n*64)+14];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b15 = _mm_broadcast_ss(&B[(l_n*64)+15]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b15 = _mm_load_ss(&B[(l_n*64)+15]);    b15 = _mm_shuffle_ps(b15, b15, 0x00);
#endif
    __m128 c15_0 = _mm_load_ss(&C[(l_n*64)+5]);
    __m128 a15_0 = _mm_load_ss(&A[66]);
    c15_0 = _mm_add_ss(c15_0, _mm_mul_ss(a15_0, b15));
    _mm_store_ss(&C[(l_n*64)+5], c15_0);
    __m128 c15_1 = _mm_load_ss(&C[(l_n*64)+15]);
    __m128 a15_1 = _mm_load_ss(&A[67]);
    c15_1 = _mm_add_ss(c15_1, _mm_mul_ss(a15_1, b15));
    _mm_store_ss(&C[(l_n*64)+15], c15_1);
    __m128 c15_2 = _mm_load_ss(&C[(l_n*64)+30]);
    __m128 a15_2 = _mm_load_ss(&A[68]);
    c15_2 = _mm_add_ss(c15_2, _mm_mul_ss(a15_2, b15));
    _mm_store_ss(&C[(l_n*64)+30], c15_2);
    __m128 c15_3 = _mm_load_ss(&C[(l_n*64)+51]);
    __m128 a15_3 = _mm_load_ss(&A[69]);
    c15_3 = _mm_add_ss(c15_3, _mm_mul_ss(a15_3, b15));
    _mm_store_ss(&C[(l_n*64)+51], c15_3);
#else
    C[(l_n*64)+5] += A[66] * B[(l_n*64)+15];
    C[(l_n*64)+15] += A[67] * B[(l_n*64)+15];
    C[(l_n*64)+30] += A[68] * B[(l_n*64)+15];
    C[(l_n*64)+51] += A[69] * B[(l_n*64)+15];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b16 = _mm_broadcast_ss(&B[(l_n*64)+16]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b16 = _mm_load_ss(&B[(l_n*64)+16]);    b16 = _mm_shuffle_ps(b16, b16, 0x00);
#endif
    __m128 c16_0 = _mm_load_ss(&C[(l_n*64)+6]);
    __m128 a16_0 = _mm_load_ss(&A[70]);
    c16_0 = _mm_add_ss(c16_0, _mm_mul_ss(a16_0, b16));
    _mm_store_ss(&C[(l_n*64)+6], c16_0);
    __m128 c16_1 = _mm_load_ss(&C[(l_n*64)+16]);
    __m128 a16_1 = _mm_load_ss(&A[71]);
    c16_1 = _mm_add_ss(c16_1, _mm_mul_ss(a16_1, b16));
    _mm_store_ss(&C[(l_n*64)+16], c16_1);
    __m128 c16_2 = _mm_load_ss(&C[(l_n*64)+31]);
    __m128 a16_2 = _mm_load_ss(&A[72]);
    c16_2 = _mm_add_ss(c16_2, _mm_mul_ss(a16_2, b16));
    _mm_store_ss(&C[(l_n*64)+31], c16_2);
    __m128 c16_3 = _mm_load_ss(&C[(l_n*64)+52]);
    __m128 a16_3 = _mm_load_ss(&A[73]);
    c16_3 = _mm_add_ss(c16_3, _mm_mul_ss(a16_3, b16));
    _mm_store_ss(&C[(l_n*64)+52], c16_3);
#else
    C[(l_n*64)+6] += A[70] * B[(l_n*64)+16];
    C[(l_n*64)+16] += A[71] * B[(l_n*64)+16];
    C[(l_n*64)+31] += A[72] * B[(l_n*64)+16];
    C[(l_n*64)+52] += A[73] * B[(l_n*64)+16];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b17 = _mm_broadcast_ss(&B[(l_n*64)+17]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b17 = _mm_load_ss(&B[(l_n*64)+17]);    b17 = _mm_shuffle_ps(b17, b17, 0x00);
#endif
    __m128 c17_0 = _mm_load_ss(&C[(l_n*64)+1]);
    __m128 a17_0 = _mm_load_ss(&A[74]);
    c17_0 = _mm_add_ss(c17_0, _mm_mul_ss(a17_0, b17));
    _mm_store_ss(&C[(l_n*64)+1], c17_0);
    __m128 c17_1 = _mm_load_ss(&C[(l_n*64)+7]);
    __m128 a17_1 = _mm_load_ss(&A[75]);
    c17_1 = _mm_add_ss(c17_1, _mm_mul_ss(a17_1, b17));
    _mm_store_ss(&C[(l_n*64)+7], c17_1);
    __m128 c17_2 = _mm_load_ss(&C[(l_n*64)+17]);
    __m128 a17_2 = _mm_load_ss(&A[76]);
    c17_2 = _mm_add_ss(c17_2, _mm_mul_ss(a17_2, b17));
    _mm_store_ss(&C[(l_n*64)+17], c17_2);
    __m128 c17_3 = _mm_load_ss(&C[(l_n*64)+32]);
    __m128 a17_3 = _mm_load_ss(&A[77]);
    c17_3 = _mm_add_ss(c17_3, _mm_mul_ss(a17_3, b17));
    _mm_store_ss(&C[(l_n*64)+32], c17_3);
    __m128 c17_4 = _mm_load_ss(&C[(l_n*64)+53]);
    __m128 a17_4 = _mm_load_ss(&A[78]);
    c17_4 = _mm_add_ss(c17_4, _mm_mul_ss(a17_4, b17));
    _mm_store_ss(&C[(l_n*64)+53], c17_4);
#else
    C[(l_n*64)+1] += A[74] * B[(l_n*64)+17];
    C[(l_n*64)+7] += A[75] * B[(l_n*64)+17];
    C[(l_n*64)+17] += A[76] * B[(l_n*64)+17];
    C[(l_n*64)+32] += A[77] * B[(l_n*64)+17];
    C[(l_n*64)+53] += A[78] * B[(l_n*64)+17];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b18 = _mm_broadcast_ss(&B[(l_n*64)+18]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b18 = _mm_load_ss(&B[(l_n*64)+18]);    b18 = _mm_shuffle_ps(b18, b18, 0x00);
#endif
    __m128 c18_0 = _mm_load_ss(&C[(l_n*64)+2]);
    __m128 a18_0 = _mm_load_ss(&A[79]);
    c18_0 = _mm_add_ss(c18_0, _mm_mul_ss(a18_0, b18));
    _mm_store_ss(&C[(l_n*64)+2], c18_0);
    __m128 c18_1 = _mm_load_ss(&C[(l_n*64)+8]);
    __m128 a18_1 = _mm_load_ss(&A[80]);
    c18_1 = _mm_add_ss(c18_1, _mm_mul_ss(a18_1, b18));
    _mm_store_ss(&C[(l_n*64)+8], c18_1);
    __m128 c18_2 = _mm_load_ss(&C[(l_n*64)+18]);
    __m128 a18_2 = _mm_load_ss(&A[81]);
    c18_2 = _mm_add_ss(c18_2, _mm_mul_ss(a18_2, b18));
    _mm_store_ss(&C[(l_n*64)+18], c18_2);
    __m128 c18_3 = _mm_load_ss(&C[(l_n*64)+33]);
    __m128 a18_3 = _mm_load_ss(&A[82]);
    c18_3 = _mm_add_ss(c18_3, _mm_mul_ss(a18_3, b18));
    _mm_store_ss(&C[(l_n*64)+33], c18_3);
    __m128 c18_4 = _mm_load_ss(&C[(l_n*64)+54]);
    __m128 a18_4 = _mm_load_ss(&A[83]);
    c18_4 = _mm_add_ss(c18_4, _mm_mul_ss(a18_4, b18));
    _mm_store_ss(&C[(l_n*64)+54], c18_4);
#else
    C[(l_n*64)+2] += A[79] * B[(l_n*64)+18];
    C[(l_n*64)+8] += A[80] * B[(l_n*64)+18];
    C[(l_n*64)+18] += A[81] * B[(l_n*64)+18];
    C[(l_n*64)+33] += A[82] * B[(l_n*64)+18];
    C[(l_n*64)+54] += A[83] * B[(l_n*64)+18];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b19 = _mm_broadcast_ss(&B[(l_n*64)+19]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b19 = _mm_load_ss(&B[(l_n*64)+19]);    b19 = _mm_shuffle_ps(b19, b19, 0x00);
#endif
    __m128 c19_0 = _mm_load_ss(&C[(l_n*64)+0]);
    __m128 a19_0 = _mm_load_ss(&A[84]);
    c19_0 = _mm_add_ss(c19_0, _mm_mul_ss(a19_0, b19));
    _mm_store_ss(&C[(l_n*64)+0], c19_0);
    __m128 c19_1 = _mm_load_ss(&C[(l_n*64)+3]);
    __m128 a19_1 = _mm_load_ss(&A[85]);
    c19_1 = _mm_add_ss(c19_1, _mm_mul_ss(a19_1, b19));
    _mm_store_ss(&C[(l_n*64)+3], c19_1);
    __m128 c19_2 = _mm_load_ss(&C[(l_n*64)+9]);
    __m128 a19_2 = _mm_load_ss(&A[86]);
    c19_2 = _mm_add_ss(c19_2, _mm_mul_ss(a19_2, b19));
    _mm_store_ss(&C[(l_n*64)+9], c19_2);
    __m128 c19_3 = _mm_load_ss(&C[(l_n*64)+19]);
    __m128 a19_3 = _mm_load_ss(&A[87]);
    c19_3 = _mm_add_ss(c19_3, _mm_mul_ss(a19_3, b19));
    _mm_store_ss(&C[(l_n*64)+19], c19_3);
    __m128 c19_4 = _mm_load_ss(&C[(l_n*64)+34]);
    __m128 a19_4 = _mm_load_ss(&A[88]);
    c19_4 = _mm_add_ss(c19_4, _mm_mul_ss(a19_4, b19));
    _mm_store_ss(&C[(l_n*64)+34], c19_4);
    __m128 c19_5 = _mm_load_ss(&C[(l_n*64)+55]);
    __m128 a19_5 = _mm_load_ss(&A[89]);
    c19_5 = _mm_add_ss(c19_5, _mm_mul_ss(a19_5, b19));
    _mm_store_ss(&C[(l_n*64)+55], c19_5);
#else
    C[(l_n*64)+0] += A[84] * B[(l_n*64)+19];
    C[(l_n*64)+3] += A[85] * B[(l_n*64)+19];
    C[(l_n*64)+9] += A[86] * B[(l_n*64)+19];
    C[(l_n*64)+19] += A[87] * B[(l_n*64)+19];
    C[(l_n*64)+34] += A[88] * B[(l_n*64)+19];
    C[(l_n*64)+55] += A[89] * B[(l_n*64)+19];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b20 = _mm_broadcast_ss(&B[(l_n*64)+20]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b20 = _mm_load_ss(&B[(l_n*64)+20]);    b20 = _mm_shuffle_ps(b20, b20, 0x00);
#endif
    __m128 c20_0 = _mm_load_ss(&C[(l_n*64)+20]);
    __m128 a20_0 = _mm_load_ss(&A[90]);
    c20_0 = _mm_add_ss(c20_0, _mm_mul_ss(a20_0, b20));
    _mm_store_ss(&C[(l_n*64)+20], c20_0);
    __m128 c20_1 = _mm_load_ss(&C[(l_n*64)+41]);
    __m128 a20_1 = _mm_load_ss(&A[91]);
    c20_1 = _mm_add_ss(c20_1, _mm_mul_ss(a20_1, b20));
    _mm_store_ss(&C[(l_n*64)+41], c20_1);
#else
    C[(l_n*64)+20] += A[90] * B[(l_n*64)+20];
    C[(l_n*64)+41] += A[91] * B[(l_n*64)+20];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b21 = _mm_broadcast_ss(&B[(l_n*64)+21]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b21 = _mm_load_ss(&B[(l_n*64)+21]);    b21 = _mm_shuffle_ps(b21, b21, 0x00);
#endif
    __m128 c21_0 = _mm_load_ss(&C[(l_n*64)+21]);
    __m128 a21_0 = _mm_load_ss(&A[92]);
    c21_0 = _mm_add_ss(c21_0, _mm_mul_ss(a21_0, b21));
    _mm_store_ss(&C[(l_n*64)+21], c21_0);
    __m128 c21_1 = _mm_load_ss(&C[(l_n*64)+42]);
    __m128 a21_1 = _mm_load_ss(&A[93]);
    c21_1 = _mm_add_ss(c21_1, _mm_mul_ss(a21_1, b21));
    _mm_store_ss(&C[(l_n*64)+42], c21_1);
#else
    C[(l_n*64)+21] += A[92] * B[(l_n*64)+21];
    C[(l_n*64)+42] += A[93] * B[(l_n*64)+21];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b22 = _mm_broadcast_ss(&B[(l_n*64)+22]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b22 = _mm_load_ss(&B[(l_n*64)+22]);    b22 = _mm_shuffle_ps(b22, b22, 0x00);
#endif
    __m128 c22_0 = _mm_load_ss(&C[(l_n*64)+22]);
    __m128 a22_0 = _mm_load_ss(&A[94]);
    c22_0 = _mm_add_ss(c22_0, _mm_mul_ss(a22_0, b22));
    _mm_store_ss(&C[(l_n*64)+22], c22_0);
    __m128 c22_1 = _mm_load_ss(&C[(l_n*64)+43]);
    __m128 a22_1 = _mm_load_ss(&A[95]);
    c22_1 = _mm_add_ss(c22_1, _mm_mul_ss(a22_1, b22));
    _mm_store_ss(&C[(l_n*64)+43], c22_1);
#else
    C[(l_n*64)+22] += A[94] * B[(l_n*64)+22];
    C[(l_n*64)+43] += A[95] * B[(l_n*64)+22];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b23 = _mm_broadcast_ss(&B[(l_n*64)+23]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b23 = _mm_load_ss(&B[(l_n*64)+23]);    b23 = _mm_shuffle_ps(b23, b23, 0x00);
#endif
    __m128 c23_0 = _mm_load_ss(&C[(l_n*64)+23]);
    __m128 a23_0 = _mm_load_ss(&A[96]);
    c23_0 = _mm_add_ss(c23_0, _mm_mul_ss(a23_0, b23));
    _mm_store_ss(&C[(l_n*64)+23], c23_0);
    __m128 c23_1 = _mm_load_ss(&C[(l_n*64)+44]);
    __m128 a23_1 = _mm_load_ss(&A[97]);
    c23_1 = _mm_add_ss(c23_1, _mm_mul_ss(a23_1, b23));
    _mm_store_ss(&C[(l_n*64)+44], c23_1);
#else
    C[(l_n*64)+23] += A[96] * B[(l_n*64)+23];
    C[(l_n*64)+44] += A[97] * B[(l_n*64)+23];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b24 = _mm_broadcast_ss(&B[(l_n*64)+24]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b24 = _mm_load_ss(&B[(l_n*64)+24]);    b24 = _mm_shuffle_ps(b24, b24, 0x00);
#endif
    __m128 c24_0 = _mm_load_ss(&C[(l_n*64)+24]);
    __m128 a24_0 = _mm_load_ss(&A[98]);
    c24_0 = _mm_add_ss(c24_0, _mm_mul_ss(a24_0, b24));
    _mm_store_ss(&C[(l_n*64)+24], c24_0);
    __m128 c24_1 = _mm_load_ss(&C[(l_n*64)+45]);
    __m128 a24_1 = _mm_load_ss(&A[99]);
    c24_1 = _mm_add_ss(c24_1, _mm_mul_ss(a24_1, b24));
    _mm_store_ss(&C[(l_n*64)+45], c24_1);
#else
    C[(l_n*64)+24] += A[98] * B[(l_n*64)+24];
    C[(l_n*64)+45] += A[99] * B[(l_n*64)+24];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b25 = _mm_broadcast_ss(&B[(l_n*64)+25]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b25 = _mm_load_ss(&B[(l_n*64)+25]);    b25 = _mm_shuffle_ps(b25, b25, 0x00);
#endif
    __m128 c25_0 = _mm_load_ss(&C[(l_n*64)+10]);
    __m128 a25_0 = _mm_load_ss(&A[100]);
    c25_0 = _mm_add_ss(c25_0, _mm_mul_ss(a25_0, b25));
    _mm_store_ss(&C[(l_n*64)+10], c25_0);
    __m128 c25_1 = _mm_load_ss(&C[(l_n*64)+25]);
    __m128 a25_1 = _mm_load_ss(&A[101]);
    c25_1 = _mm_add_ss(c25_1, _mm_mul_ss(a25_1, b25));
    _mm_store_ss(&C[(l_n*64)+25], c25_1);
    __m128 c25_2 = _mm_load_ss(&C[(l_n*64)+46]);
    __m128 a25_2 = _mm_load_ss(&A[102]);
    c25_2 = _mm_add_ss(c25_2, _mm_mul_ss(a25_2, b25));
    _mm_store_ss(&C[(l_n*64)+46], c25_2);
#else
    C[(l_n*64)+10] += A[100] * B[(l_n*64)+25];
    C[(l_n*64)+25] += A[101] * B[(l_n*64)+25];
    C[(l_n*64)+46] += A[102] * B[(l_n*64)+25];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b26 = _mm_broadcast_ss(&B[(l_n*64)+26]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b26 = _mm_load_ss(&B[(l_n*64)+26]);    b26 = _mm_shuffle_ps(b26, b26, 0x00);
#endif
    __m128 c26_0 = _mm_load_ss(&C[(l_n*64)+11]);
    __m128 a26_0 = _mm_load_ss(&A[103]);
    c26_0 = _mm_add_ss(c26_0, _mm_mul_ss(a26_0, b26));
    _mm_store_ss(&C[(l_n*64)+11], c26_0);
    __m128 c26_1 = _mm_load_ss(&C[(l_n*64)+26]);
    __m128 a26_1 = _mm_load_ss(&A[104]);
    c26_1 = _mm_add_ss(c26_1, _mm_mul_ss(a26_1, b26));
    _mm_store_ss(&C[(l_n*64)+26], c26_1);
    __m128 c26_2 = _mm_load_ss(&C[(l_n*64)+47]);
    __m128 a26_2 = _mm_load_ss(&A[105]);
    c26_2 = _mm_add_ss(c26_2, _mm_mul_ss(a26_2, b26));
    _mm_store_ss(&C[(l_n*64)+47], c26_2);
#else
    C[(l_n*64)+11] += A[103] * B[(l_n*64)+26];
    C[(l_n*64)+26] += A[104] * B[(l_n*64)+26];
    C[(l_n*64)+47] += A[105] * B[(l_n*64)+26];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b27 = _mm_broadcast_ss(&B[(l_n*64)+27]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b27 = _mm_load_ss(&B[(l_n*64)+27]);    b27 = _mm_shuffle_ps(b27, b27, 0x00);
#endif
    __m128 c27_0 = _mm_load_ss(&C[(l_n*64)+12]);
    __m128 a27_0 = _mm_load_ss(&A[106]);
    c27_0 = _mm_add_ss(c27_0, _mm_mul_ss(a27_0, b27));
    _mm_store_ss(&C[(l_n*64)+12], c27_0);
    __m128 c27_1 = _mm_load_ss(&C[(l_n*64)+27]);
    __m128 a27_1 = _mm_load_ss(&A[107]);
    c27_1 = _mm_add_ss(c27_1, _mm_mul_ss(a27_1, b27));
    _mm_store_ss(&C[(l_n*64)+27], c27_1);
    __m128 c27_2 = _mm_load_ss(&C[(l_n*64)+48]);
    __m128 a27_2 = _mm_load_ss(&A[108]);
    c27_2 = _mm_add_ss(c27_2, _mm_mul_ss(a27_2, b27));
    _mm_store_ss(&C[(l_n*64)+48], c27_2);
#else
    C[(l_n*64)+12] += A[106] * B[(l_n*64)+27];
    C[(l_n*64)+27] += A[107] * B[(l_n*64)+27];
    C[(l_n*64)+48] += A[108] * B[(l_n*64)+27];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b28 = _mm_broadcast_ss(&B[(l_n*64)+28]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b28 = _mm_load_ss(&B[(l_n*64)+28]);    b28 = _mm_shuffle_ps(b28, b28, 0x00);
#endif
    __m128 c28_0 = _mm_load_ss(&C[(l_n*64)+13]);
    __m128 a28_0 = _mm_load_ss(&A[109]);
    c28_0 = _mm_add_ss(c28_0, _mm_mul_ss(a28_0, b28));
    _mm_store_ss(&C[(l_n*64)+13], c28_0);
    __m128 c28_1 = _mm_load_ss(&C[(l_n*64)+28]);
    __m128 a28_1 = _mm_load_ss(&A[110]);
    c28_1 = _mm_add_ss(c28_1, _mm_mul_ss(a28_1, b28));
    _mm_store_ss(&C[(l_n*64)+28], c28_1);
    __m128 c28_2 = _mm_load_ss(&C[(l_n*64)+49]);
    __m128 a28_2 = _mm_load_ss(&A[111]);
    c28_2 = _mm_add_ss(c28_2, _mm_mul_ss(a28_2, b28));
    _mm_store_ss(&C[(l_n*64)+49], c28_2);
#else
    C[(l_n*64)+13] += A[109] * B[(l_n*64)+28];
    C[(l_n*64)+28] += A[110] * B[(l_n*64)+28];
    C[(l_n*64)+49] += A[111] * B[(l_n*64)+28];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b29 = _mm_broadcast_ss(&B[(l_n*64)+29]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b29 = _mm_load_ss(&B[(l_n*64)+29]);    b29 = _mm_shuffle_ps(b29, b29, 0x00);
#endif
    __m128 c29_0 = _mm_load_ss(&C[(l_n*64)+4]);
    __m128 a29_0 = _mm_load_ss(&A[112]);
    c29_0 = _mm_add_ss(c29_0, _mm_mul_ss(a29_0, b29));
    _mm_store_ss(&C[(l_n*64)+4], c29_0);
    __m128 c29_1 = _mm_load_ss(&C[(l_n*64)+14]);
    __m128 a29_1 = _mm_load_ss(&A[113]);
    c29_1 = _mm_add_ss(c29_1, _mm_mul_ss(a29_1, b29));
    _mm_store_ss(&C[(l_n*64)+14], c29_1);
    __m128 c29_2 = _mm_load_ss(&C[(l_n*64)+29]);
    __m128 a29_2 = _mm_load_ss(&A[114]);
    c29_2 = _mm_add_ss(c29_2, _mm_mul_ss(a29_2, b29));
    _mm_store_ss(&C[(l_n*64)+29], c29_2);
    __m128 c29_3 = _mm_load_ss(&C[(l_n*64)+50]);
    __m128 a29_3 = _mm_load_ss(&A[115]);
    c29_3 = _mm_add_ss(c29_3, _mm_mul_ss(a29_3, b29));
    _mm_store_ss(&C[(l_n*64)+50], c29_3);
#else
    C[(l_n*64)+4] += A[112] * B[(l_n*64)+29];
    C[(l_n*64)+14] += A[113] * B[(l_n*64)+29];
    C[(l_n*64)+29] += A[114] * B[(l_n*64)+29];
    C[(l_n*64)+50] += A[115] * B[(l_n*64)+29];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b30 = _mm_broadcast_ss(&B[(l_n*64)+30]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b30 = _mm_load_ss(&B[(l_n*64)+30]);    b30 = _mm_shuffle_ps(b30, b30, 0x00);
#endif
    __m128 c30_0 = _mm_load_ss(&C[(l_n*64)+5]);
    __m128 a30_0 = _mm_load_ss(&A[116]);
    c30_0 = _mm_add_ss(c30_0, _mm_mul_ss(a30_0, b30));
    _mm_store_ss(&C[(l_n*64)+5], c30_0);
    __m128 c30_1 = _mm_load_ss(&C[(l_n*64)+15]);
    __m128 a30_1 = _mm_load_ss(&A[117]);
    c30_1 = _mm_add_ss(c30_1, _mm_mul_ss(a30_1, b30));
    _mm_store_ss(&C[(l_n*64)+15], c30_1);
    __m128 c30_2 = _mm_load_ss(&C[(l_n*64)+30]);
    __m128 a30_2 = _mm_load_ss(&A[118]);
    c30_2 = _mm_add_ss(c30_2, _mm_mul_ss(a30_2, b30));
    _mm_store_ss(&C[(l_n*64)+30], c30_2);
    __m128 c30_3 = _mm_load_ss(&C[(l_n*64)+51]);
    __m128 a30_3 = _mm_load_ss(&A[119]);
    c30_3 = _mm_add_ss(c30_3, _mm_mul_ss(a30_3, b30));
    _mm_store_ss(&C[(l_n*64)+51], c30_3);
#else
    C[(l_n*64)+5] += A[116] * B[(l_n*64)+30];
    C[(l_n*64)+15] += A[117] * B[(l_n*64)+30];
    C[(l_n*64)+30] += A[118] * B[(l_n*64)+30];
    C[(l_n*64)+51] += A[119] * B[(l_n*64)+30];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b31 = _mm_broadcast_ss(&B[(l_n*64)+31]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b31 = _mm_load_ss(&B[(l_n*64)+31]);    b31 = _mm_shuffle_ps(b31, b31, 0x00);
#endif
    __m128 c31_0 = _mm_load_ss(&C[(l_n*64)+6]);
    __m128 a31_0 = _mm_load_ss(&A[120]);
    c31_0 = _mm_add_ss(c31_0, _mm_mul_ss(a31_0, b31));
    _mm_store_ss(&C[(l_n*64)+6], c31_0);
    __m128 c31_1 = _mm_load_ss(&C[(l_n*64)+16]);
    __m128 a31_1 = _mm_load_ss(&A[121]);
    c31_1 = _mm_add_ss(c31_1, _mm_mul_ss(a31_1, b31));
    _mm_store_ss(&C[(l_n*64)+16], c31_1);
    __m128 c31_2 = _mm_load_ss(&C[(l_n*64)+31]);
    __m128 a31_2 = _mm_load_ss(&A[122]);
    c31_2 = _mm_add_ss(c31_2, _mm_mul_ss(a31_2, b31));
    _mm_store_ss(&C[(l_n*64)+31], c31_2);
    __m128 c31_3 = _mm_load_ss(&C[(l_n*64)+52]);
    __m128 a31_3 = _mm_load_ss(&A[123]);
    c31_3 = _mm_add_ss(c31_3, _mm_mul_ss(a31_3, b31));
    _mm_store_ss(&C[(l_n*64)+52], c31_3);
#else
    C[(l_n*64)+6] += A[120] * B[(l_n*64)+31];
    C[(l_n*64)+16] += A[121] * B[(l_n*64)+31];
    C[(l_n*64)+31] += A[122] * B[(l_n*64)+31];
    C[(l_n*64)+52] += A[123] * B[(l_n*64)+31];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b32 = _mm_broadcast_ss(&B[(l_n*64)+32]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b32 = _mm_load_ss(&B[(l_n*64)+32]);    b32 = _mm_shuffle_ps(b32, b32, 0x00);
#endif
    __m128 c32_0 = _mm_load_ss(&C[(l_n*64)+1]);
    __m128 a32_0 = _mm_load_ss(&A[124]);
    c32_0 = _mm_add_ss(c32_0, _mm_mul_ss(a32_0, b32));
    _mm_store_ss(&C[(l_n*64)+1], c32_0);
    __m128 c32_1 = _mm_load_ss(&C[(l_n*64)+7]);
    __m128 a32_1 = _mm_load_ss(&A[125]);
    c32_1 = _mm_add_ss(c32_1, _mm_mul_ss(a32_1, b32));
    _mm_store_ss(&C[(l_n*64)+7], c32_1);
    __m128 c32_2 = _mm_load_ss(&C[(l_n*64)+17]);
    __m128 a32_2 = _mm_load_ss(&A[126]);
    c32_2 = _mm_add_ss(c32_2, _mm_mul_ss(a32_2, b32));
    _mm_store_ss(&C[(l_n*64)+17], c32_2);
    __m128 c32_3 = _mm_load_ss(&C[(l_n*64)+32]);
    __m128 a32_3 = _mm_load_ss(&A[127]);
    c32_3 = _mm_add_ss(c32_3, _mm_mul_ss(a32_3, b32));
    _mm_store_ss(&C[(l_n*64)+32], c32_3);
    __m128 c32_4 = _mm_load_ss(&C[(l_n*64)+53]);
    __m128 a32_4 = _mm_load_ss(&A[128]);
    c32_4 = _mm_add_ss(c32_4, _mm_mul_ss(a32_4, b32));
    _mm_store_ss(&C[(l_n*64)+53], c32_4);
#else
    C[(l_n*64)+1] += A[124] * B[(l_n*64)+32];
    C[(l_n*64)+7] += A[125] * B[(l_n*64)+32];
    C[(l_n*64)+17] += A[126] * B[(l_n*64)+32];
    C[(l_n*64)+32] += A[127] * B[(l_n*64)+32];
    C[(l_n*64)+53] += A[128] * B[(l_n*64)+32];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b33 = _mm_broadcast_ss(&B[(l_n*64)+33]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b33 = _mm_load_ss(&B[(l_n*64)+33]);    b33 = _mm_shuffle_ps(b33, b33, 0x00);
#endif
    __m128 c33_0 = _mm_load_ss(&C[(l_n*64)+2]);
    __m128 a33_0 = _mm_load_ss(&A[129]);
    c33_0 = _mm_add_ss(c33_0, _mm_mul_ss(a33_0, b33));
    _mm_store_ss(&C[(l_n*64)+2], c33_0);
    __m128 c33_1 = _mm_load_ss(&C[(l_n*64)+8]);
    __m128 a33_1 = _mm_load_ss(&A[130]);
    c33_1 = _mm_add_ss(c33_1, _mm_mul_ss(a33_1, b33));
    _mm_store_ss(&C[(l_n*64)+8], c33_1);
    __m128 c33_2 = _mm_load_ss(&C[(l_n*64)+18]);
    __m128 a33_2 = _mm_load_ss(&A[131]);
    c33_2 = _mm_add_ss(c33_2, _mm_mul_ss(a33_2, b33));
    _mm_store_ss(&C[(l_n*64)+18], c33_2);
    __m128 c33_3 = _mm_load_ss(&C[(l_n*64)+33]);
    __m128 a33_3 = _mm_load_ss(&A[132]);
    c33_3 = _mm_add_ss(c33_3, _mm_mul_ss(a33_3, b33));
    _mm_store_ss(&C[(l_n*64)+33], c33_3);
    __m128 c33_4 = _mm_load_ss(&C[(l_n*64)+54]);
    __m128 a33_4 = _mm_load_ss(&A[133]);
    c33_4 = _mm_add_ss(c33_4, _mm_mul_ss(a33_4, b33));
    _mm_store_ss(&C[(l_n*64)+54], c33_4);
#else
    C[(l_n*64)+2] += A[129] * B[(l_n*64)+33];
    C[(l_n*64)+8] += A[130] * B[(l_n*64)+33];
    C[(l_n*64)+18] += A[131] * B[(l_n*64)+33];
    C[(l_n*64)+33] += A[132] * B[(l_n*64)+33];
    C[(l_n*64)+54] += A[133] * B[(l_n*64)+33];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b34 = _mm_broadcast_ss(&B[(l_n*64)+34]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b34 = _mm_load_ss(&B[(l_n*64)+34]);    b34 = _mm_shuffle_ps(b34, b34, 0x00);
#endif
    __m128 c34_0 = _mm_load_ss(&C[(l_n*64)+0]);
    __m128 a34_0 = _mm_load_ss(&A[134]);
    c34_0 = _mm_add_ss(c34_0, _mm_mul_ss(a34_0, b34));
    _mm_store_ss(&C[(l_n*64)+0], c34_0);
    __m128 c34_1 = _mm_load_ss(&C[(l_n*64)+3]);
    __m128 a34_1 = _mm_load_ss(&A[135]);
    c34_1 = _mm_add_ss(c34_1, _mm_mul_ss(a34_1, b34));
    _mm_store_ss(&C[(l_n*64)+3], c34_1);
    __m128 c34_2 = _mm_load_ss(&C[(l_n*64)+9]);
    __m128 a34_2 = _mm_load_ss(&A[136]);
    c34_2 = _mm_add_ss(c34_2, _mm_mul_ss(a34_2, b34));
    _mm_store_ss(&C[(l_n*64)+9], c34_2);
    __m128 c34_3 = _mm_load_ss(&C[(l_n*64)+19]);
    __m128 a34_3 = _mm_load_ss(&A[137]);
    c34_3 = _mm_add_ss(c34_3, _mm_mul_ss(a34_3, b34));
    _mm_store_ss(&C[(l_n*64)+19], c34_3);
    __m128 c34_4 = _mm_load_ss(&C[(l_n*64)+34]);
    __m128 a34_4 = _mm_load_ss(&A[138]);
    c34_4 = _mm_add_ss(c34_4, _mm_mul_ss(a34_4, b34));
    _mm_store_ss(&C[(l_n*64)+34], c34_4);
    __m128 c34_5 = _mm_load_ss(&C[(l_n*64)+55]);
    __m128 a34_5 = _mm_load_ss(&A[139]);
    c34_5 = _mm_add_ss(c34_5, _mm_mul_ss(a34_5, b34));
    _mm_store_ss(&C[(l_n*64)+55], c34_5);
#else
    C[(l_n*64)+0] += A[134] * B[(l_n*64)+34];
    C[(l_n*64)+3] += A[135] * B[(l_n*64)+34];
    C[(l_n*64)+9] += A[136] * B[(l_n*64)+34];
    C[(l_n*64)+19] += A[137] * B[(l_n*64)+34];
    C[(l_n*64)+34] += A[138] * B[(l_n*64)+34];
    C[(l_n*64)+55] += A[139] * B[(l_n*64)+34];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b35 = _mm_broadcast_ss(&B[(l_n*64)+35]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b35 = _mm_load_ss(&B[(l_n*64)+35]);    b35 = _mm_shuffle_ps(b35, b35, 0x00);
#endif
    __m128 c35_0 = _mm_load_ss(&C[(l_n*64)+35]);
    __m128 a35_0 = _mm_load_ss(&A[140]);
    c35_0 = _mm_add_ss(c35_0, _mm_mul_ss(a35_0, b35));
    _mm_store_ss(&C[(l_n*64)+35], c35_0);
#else
    C[(l_n*64)+35] += A[140] * B[(l_n*64)+35];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b36 = _mm_broadcast_ss(&B[(l_n*64)+36]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b36 = _mm_load_ss(&B[(l_n*64)+36]);    b36 = _mm_shuffle_ps(b36, b36, 0x00);
#endif
    __m128 c36_0 = _mm_load_ss(&C[(l_n*64)+36]);
    __m128 a36_0 = _mm_load_ss(&A[141]);
    c36_0 = _mm_add_ss(c36_0, _mm_mul_ss(a36_0, b36));
    _mm_store_ss(&C[(l_n*64)+36], c36_0);
#else
    C[(l_n*64)+36] += A[141] * B[(l_n*64)+36];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b37 = _mm_broadcast_ss(&B[(l_n*64)+37]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b37 = _mm_load_ss(&B[(l_n*64)+37]);    b37 = _mm_shuffle_ps(b37, b37, 0x00);
#endif
    __m128 c37_0 = _mm_load_ss(&C[(l_n*64)+37]);
    __m128 a37_0 = _mm_load_ss(&A[142]);
    c37_0 = _mm_add_ss(c37_0, _mm_mul_ss(a37_0, b37));
    _mm_store_ss(&C[(l_n*64)+37], c37_0);
#else
    C[(l_n*64)+37] += A[142] * B[(l_n*64)+37];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b38 = _mm_broadcast_ss(&B[(l_n*64)+38]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b38 = _mm_load_ss(&B[(l_n*64)+38]);    b38 = _mm_shuffle_ps(b38, b38, 0x00);
#endif
    __m128 c38_0 = _mm_load_ss(&C[(l_n*64)+38]);
    __m128 a38_0 = _mm_load_ss(&A[143]);
    c38_0 = _mm_add_ss(c38_0, _mm_mul_ss(a38_0, b38));
    _mm_store_ss(&C[(l_n*64)+38], c38_0);
#else
    C[(l_n*64)+38] += A[143] * B[(l_n*64)+38];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b39 = _mm_broadcast_ss(&B[(l_n*64)+39]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b39 = _mm_load_ss(&B[(l_n*64)+39]);    b39 = _mm_shuffle_ps(b39, b39, 0x00);
#endif
    __m128 c39_0 = _mm_load_ss(&C[(l_n*64)+39]);
    __m128 a39_0 = _mm_load_ss(&A[144]);
    c39_0 = _mm_add_ss(c39_0, _mm_mul_ss(a39_0, b39));
    _mm_store_ss(&C[(l_n*64)+39], c39_0);
#else
    C[(l_n*64)+39] += A[144] * B[(l_n*64)+39];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b40 = _mm_broadcast_ss(&B[(l_n*64)+40]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b40 = _mm_load_ss(&B[(l_n*64)+40]);    b40 = _mm_shuffle_ps(b40, b40, 0x00);
#endif
    __m128 c40_0 = _mm_load_ss(&C[(l_n*64)+40]);
    __m128 a40_0 = _mm_load_ss(&A[145]);
    c40_0 = _mm_add_ss(c40_0, _mm_mul_ss(a40_0, b40));
    _mm_store_ss(&C[(l_n*64)+40], c40_0);
#else
    C[(l_n*64)+40] += A[145] * B[(l_n*64)+40];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b41 = _mm_broadcast_ss(&B[(l_n*64)+41]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b41 = _mm_load_ss(&B[(l_n*64)+41]);    b41 = _mm_shuffle_ps(b41, b41, 0x00);
#endif
    __m128 c41_0 = _mm_load_ss(&C[(l_n*64)+20]);
    __m128 a41_0 = _mm_load_ss(&A[146]);
    c41_0 = _mm_add_ss(c41_0, _mm_mul_ss(a41_0, b41));
    _mm_store_ss(&C[(l_n*64)+20], c41_0);
    __m128 c41_1 = _mm_load_ss(&C[(l_n*64)+41]);
    __m128 a41_1 = _mm_load_ss(&A[147]);
    c41_1 = _mm_add_ss(c41_1, _mm_mul_ss(a41_1, b41));
    _mm_store_ss(&C[(l_n*64)+41], c41_1);
#else
    C[(l_n*64)+20] += A[146] * B[(l_n*64)+41];
    C[(l_n*64)+41] += A[147] * B[(l_n*64)+41];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b42 = _mm_broadcast_ss(&B[(l_n*64)+42]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b42 = _mm_load_ss(&B[(l_n*64)+42]);    b42 = _mm_shuffle_ps(b42, b42, 0x00);
#endif
    __m128 c42_0 = _mm_load_ss(&C[(l_n*64)+21]);
    __m128 a42_0 = _mm_load_ss(&A[148]);
    c42_0 = _mm_add_ss(c42_0, _mm_mul_ss(a42_0, b42));
    _mm_store_ss(&C[(l_n*64)+21], c42_0);
    __m128 c42_1 = _mm_load_ss(&C[(l_n*64)+42]);
    __m128 a42_1 = _mm_load_ss(&A[149]);
    c42_1 = _mm_add_ss(c42_1, _mm_mul_ss(a42_1, b42));
    _mm_store_ss(&C[(l_n*64)+42], c42_1);
#else
    C[(l_n*64)+21] += A[148] * B[(l_n*64)+42];
    C[(l_n*64)+42] += A[149] * B[(l_n*64)+42];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b43 = _mm_broadcast_ss(&B[(l_n*64)+43]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b43 = _mm_load_ss(&B[(l_n*64)+43]);    b43 = _mm_shuffle_ps(b43, b43, 0x00);
#endif
    __m128 c43_0 = _mm_load_ss(&C[(l_n*64)+22]);
    __m128 a43_0 = _mm_load_ss(&A[150]);
    c43_0 = _mm_add_ss(c43_0, _mm_mul_ss(a43_0, b43));
    _mm_store_ss(&C[(l_n*64)+22], c43_0);
    __m128 c43_1 = _mm_load_ss(&C[(l_n*64)+43]);
    __m128 a43_1 = _mm_load_ss(&A[151]);
    c43_1 = _mm_add_ss(c43_1, _mm_mul_ss(a43_1, b43));
    _mm_store_ss(&C[(l_n*64)+43], c43_1);
#else
    C[(l_n*64)+22] += A[150] * B[(l_n*64)+43];
    C[(l_n*64)+43] += A[151] * B[(l_n*64)+43];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b44 = _mm_broadcast_ss(&B[(l_n*64)+44]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b44 = _mm_load_ss(&B[(l_n*64)+44]);    b44 = _mm_shuffle_ps(b44, b44, 0x00);
#endif
    __m128 c44_0 = _mm_load_ss(&C[(l_n*64)+23]);
    __m128 a44_0 = _mm_load_ss(&A[152]);
    c44_0 = _mm_add_ss(c44_0, _mm_mul_ss(a44_0, b44));
    _mm_store_ss(&C[(l_n*64)+23], c44_0);
    __m128 c44_1 = _mm_load_ss(&C[(l_n*64)+44]);
    __m128 a44_1 = _mm_load_ss(&A[153]);
    c44_1 = _mm_add_ss(c44_1, _mm_mul_ss(a44_1, b44));
    _mm_store_ss(&C[(l_n*64)+44], c44_1);
#else
    C[(l_n*64)+23] += A[152] * B[(l_n*64)+44];
    C[(l_n*64)+44] += A[153] * B[(l_n*64)+44];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b45 = _mm_broadcast_ss(&B[(l_n*64)+45]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b45 = _mm_load_ss(&B[(l_n*64)+45]);    b45 = _mm_shuffle_ps(b45, b45, 0x00);
#endif
    __m128 c45_0 = _mm_load_ss(&C[(l_n*64)+24]);
    __m128 a45_0 = _mm_load_ss(&A[154]);
    c45_0 = _mm_add_ss(c45_0, _mm_mul_ss(a45_0, b45));
    _mm_store_ss(&C[(l_n*64)+24], c45_0);
    __m128 c45_1 = _mm_load_ss(&C[(l_n*64)+45]);
    __m128 a45_1 = _mm_load_ss(&A[155]);
    c45_1 = _mm_add_ss(c45_1, _mm_mul_ss(a45_1, b45));
    _mm_store_ss(&C[(l_n*64)+45], c45_1);
#else
    C[(l_n*64)+24] += A[154] * B[(l_n*64)+45];
    C[(l_n*64)+45] += A[155] * B[(l_n*64)+45];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b46 = _mm_broadcast_ss(&B[(l_n*64)+46]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b46 = _mm_load_ss(&B[(l_n*64)+46]);    b46 = _mm_shuffle_ps(b46, b46, 0x00);
#endif
    __m128 c46_0 = _mm_load_ss(&C[(l_n*64)+10]);
    __m128 a46_0 = _mm_load_ss(&A[156]);
    c46_0 = _mm_add_ss(c46_0, _mm_mul_ss(a46_0, b46));
    _mm_store_ss(&C[(l_n*64)+10], c46_0);
    __m128 c46_1 = _mm_load_ss(&C[(l_n*64)+25]);
    __m128 a46_1 = _mm_load_ss(&A[157]);
    c46_1 = _mm_add_ss(c46_1, _mm_mul_ss(a46_1, b46));
    _mm_store_ss(&C[(l_n*64)+25], c46_1);
    __m128 c46_2 = _mm_load_ss(&C[(l_n*64)+46]);
    __m128 a46_2 = _mm_load_ss(&A[158]);
    c46_2 = _mm_add_ss(c46_2, _mm_mul_ss(a46_2, b46));
    _mm_store_ss(&C[(l_n*64)+46], c46_2);
#else
    C[(l_n*64)+10] += A[156] * B[(l_n*64)+46];
    C[(l_n*64)+25] += A[157] * B[(l_n*64)+46];
    C[(l_n*64)+46] += A[158] * B[(l_n*64)+46];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b47 = _mm_broadcast_ss(&B[(l_n*64)+47]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b47 = _mm_load_ss(&B[(l_n*64)+47]);    b47 = _mm_shuffle_ps(b47, b47, 0x00);
#endif
    __m128 c47_0 = _mm_load_ss(&C[(l_n*64)+11]);
    __m128 a47_0 = _mm_load_ss(&A[159]);
    c47_0 = _mm_add_ss(c47_0, _mm_mul_ss(a47_0, b47));
    _mm_store_ss(&C[(l_n*64)+11], c47_0);
    __m128 c47_1 = _mm_load_ss(&C[(l_n*64)+26]);
    __m128 a47_1 = _mm_load_ss(&A[160]);
    c47_1 = _mm_add_ss(c47_1, _mm_mul_ss(a47_1, b47));
    _mm_store_ss(&C[(l_n*64)+26], c47_1);
    __m128 c47_2 = _mm_load_ss(&C[(l_n*64)+47]);
    __m128 a47_2 = _mm_load_ss(&A[161]);
    c47_2 = _mm_add_ss(c47_2, _mm_mul_ss(a47_2, b47));
    _mm_store_ss(&C[(l_n*64)+47], c47_2);
#else
    C[(l_n*64)+11] += A[159] * B[(l_n*64)+47];
    C[(l_n*64)+26] += A[160] * B[(l_n*64)+47];
    C[(l_n*64)+47] += A[161] * B[(l_n*64)+47];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b48 = _mm_broadcast_ss(&B[(l_n*64)+48]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b48 = _mm_load_ss(&B[(l_n*64)+48]);    b48 = _mm_shuffle_ps(b48, b48, 0x00);
#endif
    __m128 c48_0 = _mm_load_ss(&C[(l_n*64)+12]);
    __m128 a48_0 = _mm_load_ss(&A[162]);
    c48_0 = _mm_add_ss(c48_0, _mm_mul_ss(a48_0, b48));
    _mm_store_ss(&C[(l_n*64)+12], c48_0);
    __m128 c48_1 = _mm_load_ss(&C[(l_n*64)+27]);
    __m128 a48_1 = _mm_load_ss(&A[163]);
    c48_1 = _mm_add_ss(c48_1, _mm_mul_ss(a48_1, b48));
    _mm_store_ss(&C[(l_n*64)+27], c48_1);
    __m128 c48_2 = _mm_load_ss(&C[(l_n*64)+48]);
    __m128 a48_2 = _mm_load_ss(&A[164]);
    c48_2 = _mm_add_ss(c48_2, _mm_mul_ss(a48_2, b48));
    _mm_store_ss(&C[(l_n*64)+48], c48_2);
#else
    C[(l_n*64)+12] += A[162] * B[(l_n*64)+48];
    C[(l_n*64)+27] += A[163] * B[(l_n*64)+48];
    C[(l_n*64)+48] += A[164] * B[(l_n*64)+48];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b49 = _mm_broadcast_ss(&B[(l_n*64)+49]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b49 = _mm_load_ss(&B[(l_n*64)+49]);    b49 = _mm_shuffle_ps(b49, b49, 0x00);
#endif
    __m128 c49_0 = _mm_load_ss(&C[(l_n*64)+13]);
    __m128 a49_0 = _mm_load_ss(&A[165]);
    c49_0 = _mm_add_ss(c49_0, _mm_mul_ss(a49_0, b49));
    _mm_store_ss(&C[(l_n*64)+13], c49_0);
    __m128 c49_1 = _mm_load_ss(&C[(l_n*64)+28]);
    __m128 a49_1 = _mm_load_ss(&A[166]);
    c49_1 = _mm_add_ss(c49_1, _mm_mul_ss(a49_1, b49));
    _mm_store_ss(&C[(l_n*64)+28], c49_1);
    __m128 c49_2 = _mm_load_ss(&C[(l_n*64)+49]);
    __m128 a49_2 = _mm_load_ss(&A[167]);
    c49_2 = _mm_add_ss(c49_2, _mm_mul_ss(a49_2, b49));
    _mm_store_ss(&C[(l_n*64)+49], c49_2);
#else
    C[(l_n*64)+13] += A[165] * B[(l_n*64)+49];
    C[(l_n*64)+28] += A[166] * B[(l_n*64)+49];
    C[(l_n*64)+49] += A[167] * B[(l_n*64)+49];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b50 = _mm_broadcast_ss(&B[(l_n*64)+50]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b50 = _mm_load_ss(&B[(l_n*64)+50]);    b50 = _mm_shuffle_ps(b50, b50, 0x00);
#endif
    __m128 c50_0 = _mm_load_ss(&C[(l_n*64)+4]);
    __m128 a50_0 = _mm_load_ss(&A[168]);
    c50_0 = _mm_add_ss(c50_0, _mm_mul_ss(a50_0, b50));
    _mm_store_ss(&C[(l_n*64)+4], c50_0);
    __m128 c50_1 = _mm_load_ss(&C[(l_n*64)+14]);
    __m128 a50_1 = _mm_load_ss(&A[169]);
    c50_1 = _mm_add_ss(c50_1, _mm_mul_ss(a50_1, b50));
    _mm_store_ss(&C[(l_n*64)+14], c50_1);
    __m128 c50_2 = _mm_load_ss(&C[(l_n*64)+29]);
    __m128 a50_2 = _mm_load_ss(&A[170]);
    c50_2 = _mm_add_ss(c50_2, _mm_mul_ss(a50_2, b50));
    _mm_store_ss(&C[(l_n*64)+29], c50_2);
    __m128 c50_3 = _mm_load_ss(&C[(l_n*64)+50]);
    __m128 a50_3 = _mm_load_ss(&A[171]);
    c50_3 = _mm_add_ss(c50_3, _mm_mul_ss(a50_3, b50));
    _mm_store_ss(&C[(l_n*64)+50], c50_3);
#else
    C[(l_n*64)+4] += A[168] * B[(l_n*64)+50];
    C[(l_n*64)+14] += A[169] * B[(l_n*64)+50];
    C[(l_n*64)+29] += A[170] * B[(l_n*64)+50];
    C[(l_n*64)+50] += A[171] * B[(l_n*64)+50];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b51 = _mm_broadcast_ss(&B[(l_n*64)+51]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b51 = _mm_load_ss(&B[(l_n*64)+51]);    b51 = _mm_shuffle_ps(b51, b51, 0x00);
#endif
    __m128 c51_0 = _mm_load_ss(&C[(l_n*64)+5]);
    __m128 a51_0 = _mm_load_ss(&A[172]);
    c51_0 = _mm_add_ss(c51_0, _mm_mul_ss(a51_0, b51));
    _mm_store_ss(&C[(l_n*64)+5], c51_0);
    __m128 c51_1 = _mm_load_ss(&C[(l_n*64)+15]);
    __m128 a51_1 = _mm_load_ss(&A[173]);
    c51_1 = _mm_add_ss(c51_1, _mm_mul_ss(a51_1, b51));
    _mm_store_ss(&C[(l_n*64)+15], c51_1);
    __m128 c51_2 = _mm_load_ss(&C[(l_n*64)+30]);
    __m128 a51_2 = _mm_load_ss(&A[174]);
    c51_2 = _mm_add_ss(c51_2, _mm_mul_ss(a51_2, b51));
    _mm_store_ss(&C[(l_n*64)+30], c51_2);
    __m128 c51_3 = _mm_load_ss(&C[(l_n*64)+51]);
    __m128 a51_3 = _mm_load_ss(&A[175]);
    c51_3 = _mm_add_ss(c51_3, _mm_mul_ss(a51_3, b51));
    _mm_store_ss(&C[(l_n*64)+51], c51_3);
#else
    C[(l_n*64)+5] += A[172] * B[(l_n*64)+51];
    C[(l_n*64)+15] += A[173] * B[(l_n*64)+51];
    C[(l_n*64)+30] += A[174] * B[(l_n*64)+51];
    C[(l_n*64)+51] += A[175] * B[(l_n*64)+51];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b52 = _mm_broadcast_ss(&B[(l_n*64)+52]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b52 = _mm_load_ss(&B[(l_n*64)+52]);    b52 = _mm_shuffle_ps(b52, b52, 0x00);
#endif
    __m128 c52_0 = _mm_load_ss(&C[(l_n*64)+6]);
    __m128 a52_0 = _mm_load_ss(&A[176]);
    c52_0 = _mm_add_ss(c52_0, _mm_mul_ss(a52_0, b52));
    _mm_store_ss(&C[(l_n*64)+6], c52_0);
    __m128 c52_1 = _mm_load_ss(&C[(l_n*64)+16]);
    __m128 a52_1 = _mm_load_ss(&A[177]);
    c52_1 = _mm_add_ss(c52_1, _mm_mul_ss(a52_1, b52));
    _mm_store_ss(&C[(l_n*64)+16], c52_1);
    __m128 c52_2 = _mm_load_ss(&C[(l_n*64)+31]);
    __m128 a52_2 = _mm_load_ss(&A[178]);
    c52_2 = _mm_add_ss(c52_2, _mm_mul_ss(a52_2, b52));
    _mm_store_ss(&C[(l_n*64)+31], c52_2);
    __m128 c52_3 = _mm_load_ss(&C[(l_n*64)+52]);
    __m128 a52_3 = _mm_load_ss(&A[179]);
    c52_3 = _mm_add_ss(c52_3, _mm_mul_ss(a52_3, b52));
    _mm_store_ss(&C[(l_n*64)+52], c52_3);
#else
    C[(l_n*64)+6] += A[176] * B[(l_n*64)+52];
    C[(l_n*64)+16] += A[177] * B[(l_n*64)+52];
    C[(l_n*64)+31] += A[178] * B[(l_n*64)+52];
    C[(l_n*64)+52] += A[179] * B[(l_n*64)+52];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b53 = _mm_broadcast_ss(&B[(l_n*64)+53]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b53 = _mm_load_ss(&B[(l_n*64)+53]);    b53 = _mm_shuffle_ps(b53, b53, 0x00);
#endif
    __m128 c53_0 = _mm_load_ss(&C[(l_n*64)+1]);
    __m128 a53_0 = _mm_load_ss(&A[180]);
    c53_0 = _mm_add_ss(c53_0, _mm_mul_ss(a53_0, b53));
    _mm_store_ss(&C[(l_n*64)+1], c53_0);
    __m128 c53_1 = _mm_load_ss(&C[(l_n*64)+7]);
    __m128 a53_1 = _mm_load_ss(&A[181]);
    c53_1 = _mm_add_ss(c53_1, _mm_mul_ss(a53_1, b53));
    _mm_store_ss(&C[(l_n*64)+7], c53_1);
    __m128 c53_2 = _mm_load_ss(&C[(l_n*64)+17]);
    __m128 a53_2 = _mm_load_ss(&A[182]);
    c53_2 = _mm_add_ss(c53_2, _mm_mul_ss(a53_2, b53));
    _mm_store_ss(&C[(l_n*64)+17], c53_2);
    __m128 c53_3 = _mm_load_ss(&C[(l_n*64)+32]);
    __m128 a53_3 = _mm_load_ss(&A[183]);
    c53_3 = _mm_add_ss(c53_3, _mm_mul_ss(a53_3, b53));
    _mm_store_ss(&C[(l_n*64)+32], c53_3);
    __m128 c53_4 = _mm_load_ss(&C[(l_n*64)+53]);
    __m128 a53_4 = _mm_load_ss(&A[184]);
    c53_4 = _mm_add_ss(c53_4, _mm_mul_ss(a53_4, b53));
    _mm_store_ss(&C[(l_n*64)+53], c53_4);
#else
    C[(l_n*64)+1] += A[180] * B[(l_n*64)+53];
    C[(l_n*64)+7] += A[181] * B[(l_n*64)+53];
    C[(l_n*64)+17] += A[182] * B[(l_n*64)+53];
    C[(l_n*64)+32] += A[183] * B[(l_n*64)+53];
    C[(l_n*64)+53] += A[184] * B[(l_n*64)+53];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b54 = _mm_broadcast_ss(&B[(l_n*64)+54]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b54 = _mm_load_ss(&B[(l_n*64)+54]);    b54 = _mm_shuffle_ps(b54, b54, 0x00);
#endif
    __m128 c54_0 = _mm_load_ss(&C[(l_n*64)+2]);
    __m128 a54_0 = _mm_load_ss(&A[185]);
    c54_0 = _mm_add_ss(c54_0, _mm_mul_ss(a54_0, b54));
    _mm_store_ss(&C[(l_n*64)+2], c54_0);
    __m128 c54_1 = _mm_load_ss(&C[(l_n*64)+8]);
    __m128 a54_1 = _mm_load_ss(&A[186]);
    c54_1 = _mm_add_ss(c54_1, _mm_mul_ss(a54_1, b54));
    _mm_store_ss(&C[(l_n*64)+8], c54_1);
    __m128 c54_2 = _mm_load_ss(&C[(l_n*64)+18]);
    __m128 a54_2 = _mm_load_ss(&A[187]);
    c54_2 = _mm_add_ss(c54_2, _mm_mul_ss(a54_2, b54));
    _mm_store_ss(&C[(l_n*64)+18], c54_2);
    __m128 c54_3 = _mm_load_ss(&C[(l_n*64)+33]);
    __m128 a54_3 = _mm_load_ss(&A[188]);
    c54_3 = _mm_add_ss(c54_3, _mm_mul_ss(a54_3, b54));
    _mm_store_ss(&C[(l_n*64)+33], c54_3);
    __m128 c54_4 = _mm_load_ss(&C[(l_n*64)+54]);
    __m128 a54_4 = _mm_load_ss(&A[189]);
    c54_4 = _mm_add_ss(c54_4, _mm_mul_ss(a54_4, b54));
    _mm_store_ss(&C[(l_n*64)+54], c54_4);
#else
    C[(l_n*64)+2] += A[185] * B[(l_n*64)+54];
    C[(l_n*64)+8] += A[186] * B[(l_n*64)+54];
    C[(l_n*64)+18] += A[187] * B[(l_n*64)+54];
    C[(l_n*64)+33] += A[188] * B[(l_n*64)+54];
    C[(l_n*64)+54] += A[189] * B[(l_n*64)+54];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m128 b55 = _mm_broadcast_ss(&B[(l_n*64)+55]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128 b55 = _mm_load_ss(&B[(l_n*64)+55]);    b55 = _mm_shuffle_ps(b55, b55, 0x00);
#endif
    __m128 c55_0 = _mm_load_ss(&C[(l_n*64)+0]);
    __m128 a55_0 = _mm_load_ss(&A[190]);
    c55_0 = _mm_add_ss(c55_0, _mm_mul_ss(a55_0, b55));
    _mm_store_ss(&C[(l_n*64)+0], c55_0);
    __m128 c55_1 = _mm_load_ss(&C[(l_n*64)+3]);
    __m128 a55_1 = _mm_load_ss(&A[191]);
    c55_1 = _mm_add_ss(c55_1, _mm_mul_ss(a55_1, b55));
    _mm_store_ss(&C[(l_n*64)+3], c55_1);
    __m128 c55_2 = _mm_load_ss(&C[(l_n*64)+9]);
    __m128 a55_2 = _mm_load_ss(&A[192]);
    c55_2 = _mm_add_ss(c55_2, _mm_mul_ss(a55_2, b55));
    _mm_store_ss(&C[(l_n*64)+9], c55_2);
    __m128 c55_3 = _mm_load_ss(&C[(l_n*64)+19]);
    __m128 a55_3 = _mm_load_ss(&A[193]);
    c55_3 = _mm_add_ss(c55_3, _mm_mul_ss(a55_3, b55));
    _mm_store_ss(&C[(l_n*64)+19], c55_3);
    __m128 c55_4 = _mm_load_ss(&C[(l_n*64)+34]);
    __m128 a55_4 = _mm_load_ss(&A[194]);
    c55_4 = _mm_add_ss(c55_4, _mm_mul_ss(a55_4, b55));
    _mm_store_ss(&C[(l_n*64)+34], c55_4);
    __m128 c55_5 = _mm_load_ss(&C[(l_n*64)+55]);
    __m128 a55_5 = _mm_load_ss(&A[195]);
    c55_5 = _mm_add_ss(c55_5, _mm_mul_ss(a55_5, b55));
    _mm_store_ss(&C[(l_n*64)+55], c55_5);
#else
    C[(l_n*64)+0] += A[190] * B[(l_n*64)+55];
    C[(l_n*64)+3] += A[191] * B[(l_n*64)+55];
    C[(l_n*64)+9] += A[192] * B[(l_n*64)+55];
    C[(l_n*64)+19] += A[193] * B[(l_n*64)+55];
    C[(l_n*64)+34] += A[194] * B[(l_n*64)+55];
    C[(l_n*64)+55] += A[195] * B[(l_n*64)+55];
#endif

  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 3528;
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
