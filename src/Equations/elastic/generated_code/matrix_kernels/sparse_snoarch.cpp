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
// @date 2015-11-22 19:13:49.766038
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
#ifndef SPARSESNOARCHCPP
#define SPARSESNOARCHCPP

#if defined( __SSE3__) || defined(__MIC__)
#include <immintrin.h>
#endif

#include <cstddef>
#ifndef NDEBUG
extern long long libxsmm_num_total_flops;
#endif

void ssparse_starMatrix_m1_n9_k9_ldA4_ldBna2_ldC4_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  for ( l_m = 0; l_m < 1; l_m++) {
    C[0+l_m] += A[24+l_m] * B[0];
    C[0+l_m] += A[28+l_m] * B[1];
    C[0+l_m] += A[32+l_m] * B[2];
    C[4+l_m] += A[24+l_m] * B[3];
    C[4+l_m] += A[28+l_m] * B[4];
    C[4+l_m] += A[32+l_m] * B[5];
    C[8+l_m] += A[24+l_m] * B[6];
    C[8+l_m] += A[28+l_m] * B[7];
    C[8+l_m] += A[32+l_m] * B[8];
    C[12+l_m] += A[24+l_m] * B[9];
    C[12+l_m] += A[28+l_m] * B[10];
    C[16+l_m] += A[28+l_m] * B[11];
    C[16+l_m] += A[32+l_m] * B[12];
    C[20+l_m] += A[24+l_m] * B[13];
    C[20+l_m] += A[32+l_m] * B[14];
    C[24+l_m] += A[0+l_m] * B[15];
    C[24+l_m] += A[12+l_m] * B[16];
    C[24+l_m] += A[20+l_m] * B[17];
    C[28+l_m] += A[4+l_m] * B[18];
    C[28+l_m] += A[12+l_m] * B[19];
    C[28+l_m] += A[16+l_m] * B[20];
    C[32+l_m] += A[8+l_m] * B[21];
    C[32+l_m] += A[16+l_m] * B[22];
    C[32+l_m] += A[20+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 48;
#endif
}

void ssparse_starMatrix_m4_n9_k9_ldA4_ldBna3_ldC4_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(4)
  #pragma vector aligned
  for ( l_m = 0; l_m < 4; l_m++) {
    C[0+l_m] += A[24+l_m] * B[0];
    C[0+l_m] += A[28+l_m] * B[1];
    C[0+l_m] += A[32+l_m] * B[2];
    C[4+l_m] += A[24+l_m] * B[3];
    C[4+l_m] += A[28+l_m] * B[4];
    C[4+l_m] += A[32+l_m] * B[5];
    C[8+l_m] += A[24+l_m] * B[6];
    C[8+l_m] += A[28+l_m] * B[7];
    C[8+l_m] += A[32+l_m] * B[8];
    C[12+l_m] += A[24+l_m] * B[9];
    C[12+l_m] += A[28+l_m] * B[10];
    C[16+l_m] += A[28+l_m] * B[11];
    C[16+l_m] += A[32+l_m] * B[12];
    C[20+l_m] += A[24+l_m] * B[13];
    C[20+l_m] += A[32+l_m] * B[14];
    C[24+l_m] += A[0+l_m] * B[15];
    C[24+l_m] += A[12+l_m] * B[16];
    C[24+l_m] += A[20+l_m] * B[17];
    C[28+l_m] += A[4+l_m] * B[18];
    C[28+l_m] += A[12+l_m] * B[19];
    C[28+l_m] += A[16+l_m] * B[20];
    C[32+l_m] += A[8+l_m] * B[21];
    C[32+l_m] += A[16+l_m] * B[22];
    C[32+l_m] += A[20+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 192;
#endif
}

void ssparse_starMatrix_m1_n9_k9_ldA4_ldBna3_ldC4_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  for ( l_m = 0; l_m < 1; l_m++) {
    C[0+l_m] += A[24+l_m] * B[0];
    C[0+l_m] += A[28+l_m] * B[1];
    C[0+l_m] += A[32+l_m] * B[2];
    C[4+l_m] += A[24+l_m] * B[3];
    C[4+l_m] += A[28+l_m] * B[4];
    C[4+l_m] += A[32+l_m] * B[5];
    C[8+l_m] += A[24+l_m] * B[6];
    C[8+l_m] += A[28+l_m] * B[7];
    C[8+l_m] += A[32+l_m] * B[8];
    C[12+l_m] += A[24+l_m] * B[9];
    C[12+l_m] += A[28+l_m] * B[10];
    C[16+l_m] += A[28+l_m] * B[11];
    C[16+l_m] += A[32+l_m] * B[12];
    C[20+l_m] += A[24+l_m] * B[13];
    C[20+l_m] += A[32+l_m] * B[14];
    C[24+l_m] += A[0+l_m] * B[15];
    C[24+l_m] += A[12+l_m] * B[16];
    C[24+l_m] += A[20+l_m] * B[17];
    C[28+l_m] += A[4+l_m] * B[18];
    C[28+l_m] += A[12+l_m] * B[19];
    C[28+l_m] += A[16+l_m] * B[20];
    C[32+l_m] += A[8+l_m] * B[21];
    C[32+l_m] += A[16+l_m] * B[22];
    C[32+l_m] += A[20+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 48;
#endif
}

void ssparse_starMatrix_m10_n9_k9_ldA12_ldBna4_ldC12_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 10; l_m++) {
    C[0+l_m] += A[72+l_m] * B[0];
    C[0+l_m] += A[84+l_m] * B[1];
    C[0+l_m] += A[96+l_m] * B[2];
    C[12+l_m] += A[72+l_m] * B[3];
    C[12+l_m] += A[84+l_m] * B[4];
    C[12+l_m] += A[96+l_m] * B[5];
    C[24+l_m] += A[72+l_m] * B[6];
    C[24+l_m] += A[84+l_m] * B[7];
    C[24+l_m] += A[96+l_m] * B[8];
    C[36+l_m] += A[72+l_m] * B[9];
    C[36+l_m] += A[84+l_m] * B[10];
    C[48+l_m] += A[84+l_m] * B[11];
    C[48+l_m] += A[96+l_m] * B[12];
    C[60+l_m] += A[72+l_m] * B[13];
    C[60+l_m] += A[96+l_m] * B[14];
    C[72+l_m] += A[0+l_m] * B[15];
    C[72+l_m] += A[36+l_m] * B[16];
    C[72+l_m] += A[60+l_m] * B[17];
    C[84+l_m] += A[12+l_m] * B[18];
    C[84+l_m] += A[36+l_m] * B[19];
    C[84+l_m] += A[48+l_m] * B[20];
    C[96+l_m] += A[24+l_m] * B[21];
    C[96+l_m] += A[48+l_m] * B[22];
    C[96+l_m] += A[60+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 480;
#endif
}

void ssparse_starMatrix_m4_n9_k9_ldA4_ldBna4_ldC4_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(4)
  #pragma vector aligned
  for ( l_m = 0; l_m < 4; l_m++) {
    C[0+l_m] += A[24+l_m] * B[0];
    C[0+l_m] += A[28+l_m] * B[1];
    C[0+l_m] += A[32+l_m] * B[2];
    C[4+l_m] += A[24+l_m] * B[3];
    C[4+l_m] += A[28+l_m] * B[4];
    C[4+l_m] += A[32+l_m] * B[5];
    C[8+l_m] += A[24+l_m] * B[6];
    C[8+l_m] += A[28+l_m] * B[7];
    C[8+l_m] += A[32+l_m] * B[8];
    C[12+l_m] += A[24+l_m] * B[9];
    C[12+l_m] += A[28+l_m] * B[10];
    C[16+l_m] += A[28+l_m] * B[11];
    C[16+l_m] += A[32+l_m] * B[12];
    C[20+l_m] += A[24+l_m] * B[13];
    C[20+l_m] += A[32+l_m] * B[14];
    C[24+l_m] += A[0+l_m] * B[15];
    C[24+l_m] += A[12+l_m] * B[16];
    C[24+l_m] += A[20+l_m] * B[17];
    C[28+l_m] += A[4+l_m] * B[18];
    C[28+l_m] += A[12+l_m] * B[19];
    C[28+l_m] += A[16+l_m] * B[20];
    C[32+l_m] += A[8+l_m] * B[21];
    C[32+l_m] += A[16+l_m] * B[22];
    C[32+l_m] += A[20+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 192;
#endif
}

void ssparse_starMatrix_m1_n9_k9_ldA4_ldBna4_ldC4_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  for ( l_m = 0; l_m < 1; l_m++) {
    C[0+l_m] += A[24+l_m] * B[0];
    C[0+l_m] += A[28+l_m] * B[1];
    C[0+l_m] += A[32+l_m] * B[2];
    C[4+l_m] += A[24+l_m] * B[3];
    C[4+l_m] += A[28+l_m] * B[4];
    C[4+l_m] += A[32+l_m] * B[5];
    C[8+l_m] += A[24+l_m] * B[6];
    C[8+l_m] += A[28+l_m] * B[7];
    C[8+l_m] += A[32+l_m] * B[8];
    C[12+l_m] += A[24+l_m] * B[9];
    C[12+l_m] += A[28+l_m] * B[10];
    C[16+l_m] += A[28+l_m] * B[11];
    C[16+l_m] += A[32+l_m] * B[12];
    C[20+l_m] += A[24+l_m] * B[13];
    C[20+l_m] += A[32+l_m] * B[14];
    C[24+l_m] += A[0+l_m] * B[15];
    C[24+l_m] += A[12+l_m] * B[16];
    C[24+l_m] += A[20+l_m] * B[17];
    C[28+l_m] += A[4+l_m] * B[18];
    C[28+l_m] += A[12+l_m] * B[19];
    C[28+l_m] += A[16+l_m] * B[20];
    C[32+l_m] += A[8+l_m] * B[21];
    C[32+l_m] += A[16+l_m] * B[22];
    C[32+l_m] += A[20+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 48;
#endif
}

void ssparse_starMatrix_m20_n9_k9_ldA20_ldBna5_ldC20_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 20; l_m++) {
    C[0+l_m] += A[120+l_m] * B[0];
    C[0+l_m] += A[140+l_m] * B[1];
    C[0+l_m] += A[160+l_m] * B[2];
    C[20+l_m] += A[120+l_m] * B[3];
    C[20+l_m] += A[140+l_m] * B[4];
    C[20+l_m] += A[160+l_m] * B[5];
    C[40+l_m] += A[120+l_m] * B[6];
    C[40+l_m] += A[140+l_m] * B[7];
    C[40+l_m] += A[160+l_m] * B[8];
    C[60+l_m] += A[120+l_m] * B[9];
    C[60+l_m] += A[140+l_m] * B[10];
    C[80+l_m] += A[140+l_m] * B[11];
    C[80+l_m] += A[160+l_m] * B[12];
    C[100+l_m] += A[120+l_m] * B[13];
    C[100+l_m] += A[160+l_m] * B[14];
    C[120+l_m] += A[0+l_m] * B[15];
    C[120+l_m] += A[60+l_m] * B[16];
    C[120+l_m] += A[100+l_m] * B[17];
    C[140+l_m] += A[20+l_m] * B[18];
    C[140+l_m] += A[60+l_m] * B[19];
    C[140+l_m] += A[80+l_m] * B[20];
    C[160+l_m] += A[40+l_m] * B[21];
    C[160+l_m] += A[80+l_m] * B[22];
    C[160+l_m] += A[100+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 960;
#endif
}

void ssparse_starMatrix_m10_n9_k9_ldA12_ldBna5_ldC12_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 10; l_m++) {
    C[0+l_m] += A[72+l_m] * B[0];
    C[0+l_m] += A[84+l_m] * B[1];
    C[0+l_m] += A[96+l_m] * B[2];
    C[12+l_m] += A[72+l_m] * B[3];
    C[12+l_m] += A[84+l_m] * B[4];
    C[12+l_m] += A[96+l_m] * B[5];
    C[24+l_m] += A[72+l_m] * B[6];
    C[24+l_m] += A[84+l_m] * B[7];
    C[24+l_m] += A[96+l_m] * B[8];
    C[36+l_m] += A[72+l_m] * B[9];
    C[36+l_m] += A[84+l_m] * B[10];
    C[48+l_m] += A[84+l_m] * B[11];
    C[48+l_m] += A[96+l_m] * B[12];
    C[60+l_m] += A[72+l_m] * B[13];
    C[60+l_m] += A[96+l_m] * B[14];
    C[72+l_m] += A[0+l_m] * B[15];
    C[72+l_m] += A[36+l_m] * B[16];
    C[72+l_m] += A[60+l_m] * B[17];
    C[84+l_m] += A[12+l_m] * B[18];
    C[84+l_m] += A[36+l_m] * B[19];
    C[84+l_m] += A[48+l_m] * B[20];
    C[96+l_m] += A[24+l_m] * B[21];
    C[96+l_m] += A[48+l_m] * B[22];
    C[96+l_m] += A[60+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 480;
#endif
}

void ssparse_starMatrix_m4_n9_k9_ldA4_ldBna5_ldC4_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(4)
  #pragma vector aligned
  for ( l_m = 0; l_m < 4; l_m++) {
    C[0+l_m] += A[24+l_m] * B[0];
    C[0+l_m] += A[28+l_m] * B[1];
    C[0+l_m] += A[32+l_m] * B[2];
    C[4+l_m] += A[24+l_m] * B[3];
    C[4+l_m] += A[28+l_m] * B[4];
    C[4+l_m] += A[32+l_m] * B[5];
    C[8+l_m] += A[24+l_m] * B[6];
    C[8+l_m] += A[28+l_m] * B[7];
    C[8+l_m] += A[32+l_m] * B[8];
    C[12+l_m] += A[24+l_m] * B[9];
    C[12+l_m] += A[28+l_m] * B[10];
    C[16+l_m] += A[28+l_m] * B[11];
    C[16+l_m] += A[32+l_m] * B[12];
    C[20+l_m] += A[24+l_m] * B[13];
    C[20+l_m] += A[32+l_m] * B[14];
    C[24+l_m] += A[0+l_m] * B[15];
    C[24+l_m] += A[12+l_m] * B[16];
    C[24+l_m] += A[20+l_m] * B[17];
    C[28+l_m] += A[4+l_m] * B[18];
    C[28+l_m] += A[12+l_m] * B[19];
    C[28+l_m] += A[16+l_m] * B[20];
    C[32+l_m] += A[8+l_m] * B[21];
    C[32+l_m] += A[16+l_m] * B[22];
    C[32+l_m] += A[20+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 192;
#endif
}

void ssparse_starMatrix_m1_n9_k9_ldA4_ldBna5_ldC4_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  for ( l_m = 0; l_m < 1; l_m++) {
    C[0+l_m] += A[24+l_m] * B[0];
    C[0+l_m] += A[28+l_m] * B[1];
    C[0+l_m] += A[32+l_m] * B[2];
    C[4+l_m] += A[24+l_m] * B[3];
    C[4+l_m] += A[28+l_m] * B[4];
    C[4+l_m] += A[32+l_m] * B[5];
    C[8+l_m] += A[24+l_m] * B[6];
    C[8+l_m] += A[28+l_m] * B[7];
    C[8+l_m] += A[32+l_m] * B[8];
    C[12+l_m] += A[24+l_m] * B[9];
    C[12+l_m] += A[28+l_m] * B[10];
    C[16+l_m] += A[28+l_m] * B[11];
    C[16+l_m] += A[32+l_m] * B[12];
    C[20+l_m] += A[24+l_m] * B[13];
    C[20+l_m] += A[32+l_m] * B[14];
    C[24+l_m] += A[0+l_m] * B[15];
    C[24+l_m] += A[12+l_m] * B[16];
    C[24+l_m] += A[20+l_m] * B[17];
    C[28+l_m] += A[4+l_m] * B[18];
    C[28+l_m] += A[12+l_m] * B[19];
    C[28+l_m] += A[16+l_m] * B[20];
    C[32+l_m] += A[8+l_m] * B[21];
    C[32+l_m] += A[16+l_m] * B[22];
    C[32+l_m] += A[20+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 48;
#endif
}

void ssparse_starMatrix_m35_n9_k9_ldA36_ldBna6_ldC36_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 35; l_m++) {
    C[0+l_m] += A[216+l_m] * B[0];
    C[0+l_m] += A[252+l_m] * B[1];
    C[0+l_m] += A[288+l_m] * B[2];
    C[36+l_m] += A[216+l_m] * B[3];
    C[36+l_m] += A[252+l_m] * B[4];
    C[36+l_m] += A[288+l_m] * B[5];
    C[72+l_m] += A[216+l_m] * B[6];
    C[72+l_m] += A[252+l_m] * B[7];
    C[72+l_m] += A[288+l_m] * B[8];
    C[108+l_m] += A[216+l_m] * B[9];
    C[108+l_m] += A[252+l_m] * B[10];
    C[144+l_m] += A[252+l_m] * B[11];
    C[144+l_m] += A[288+l_m] * B[12];
    C[180+l_m] += A[216+l_m] * B[13];
    C[180+l_m] += A[288+l_m] * B[14];
    C[216+l_m] += A[0+l_m] * B[15];
    C[216+l_m] += A[108+l_m] * B[16];
    C[216+l_m] += A[180+l_m] * B[17];
    C[252+l_m] += A[36+l_m] * B[18];
    C[252+l_m] += A[108+l_m] * B[19];
    C[252+l_m] += A[144+l_m] * B[20];
    C[288+l_m] += A[72+l_m] * B[21];
    C[288+l_m] += A[144+l_m] * B[22];
    C[288+l_m] += A[180+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 1680;
#endif
}

void ssparse_starMatrix_m20_n9_k9_ldA20_ldBna6_ldC20_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 20; l_m++) {
    C[0+l_m] += A[120+l_m] * B[0];
    C[0+l_m] += A[140+l_m] * B[1];
    C[0+l_m] += A[160+l_m] * B[2];
    C[20+l_m] += A[120+l_m] * B[3];
    C[20+l_m] += A[140+l_m] * B[4];
    C[20+l_m] += A[160+l_m] * B[5];
    C[40+l_m] += A[120+l_m] * B[6];
    C[40+l_m] += A[140+l_m] * B[7];
    C[40+l_m] += A[160+l_m] * B[8];
    C[60+l_m] += A[120+l_m] * B[9];
    C[60+l_m] += A[140+l_m] * B[10];
    C[80+l_m] += A[140+l_m] * B[11];
    C[80+l_m] += A[160+l_m] * B[12];
    C[100+l_m] += A[120+l_m] * B[13];
    C[100+l_m] += A[160+l_m] * B[14];
    C[120+l_m] += A[0+l_m] * B[15];
    C[120+l_m] += A[60+l_m] * B[16];
    C[120+l_m] += A[100+l_m] * B[17];
    C[140+l_m] += A[20+l_m] * B[18];
    C[140+l_m] += A[60+l_m] * B[19];
    C[140+l_m] += A[80+l_m] * B[20];
    C[160+l_m] += A[40+l_m] * B[21];
    C[160+l_m] += A[80+l_m] * B[22];
    C[160+l_m] += A[100+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 960;
#endif
}

void ssparse_starMatrix_m10_n9_k9_ldA12_ldBna6_ldC12_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 10; l_m++) {
    C[0+l_m] += A[72+l_m] * B[0];
    C[0+l_m] += A[84+l_m] * B[1];
    C[0+l_m] += A[96+l_m] * B[2];
    C[12+l_m] += A[72+l_m] * B[3];
    C[12+l_m] += A[84+l_m] * B[4];
    C[12+l_m] += A[96+l_m] * B[5];
    C[24+l_m] += A[72+l_m] * B[6];
    C[24+l_m] += A[84+l_m] * B[7];
    C[24+l_m] += A[96+l_m] * B[8];
    C[36+l_m] += A[72+l_m] * B[9];
    C[36+l_m] += A[84+l_m] * B[10];
    C[48+l_m] += A[84+l_m] * B[11];
    C[48+l_m] += A[96+l_m] * B[12];
    C[60+l_m] += A[72+l_m] * B[13];
    C[60+l_m] += A[96+l_m] * B[14];
    C[72+l_m] += A[0+l_m] * B[15];
    C[72+l_m] += A[36+l_m] * B[16];
    C[72+l_m] += A[60+l_m] * B[17];
    C[84+l_m] += A[12+l_m] * B[18];
    C[84+l_m] += A[36+l_m] * B[19];
    C[84+l_m] += A[48+l_m] * B[20];
    C[96+l_m] += A[24+l_m] * B[21];
    C[96+l_m] += A[48+l_m] * B[22];
    C[96+l_m] += A[60+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 480;
#endif
}

void ssparse_starMatrix_m4_n9_k9_ldA4_ldBna6_ldC4_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(4)
  #pragma vector aligned
  for ( l_m = 0; l_m < 4; l_m++) {
    C[0+l_m] += A[24+l_m] * B[0];
    C[0+l_m] += A[28+l_m] * B[1];
    C[0+l_m] += A[32+l_m] * B[2];
    C[4+l_m] += A[24+l_m] * B[3];
    C[4+l_m] += A[28+l_m] * B[4];
    C[4+l_m] += A[32+l_m] * B[5];
    C[8+l_m] += A[24+l_m] * B[6];
    C[8+l_m] += A[28+l_m] * B[7];
    C[8+l_m] += A[32+l_m] * B[8];
    C[12+l_m] += A[24+l_m] * B[9];
    C[12+l_m] += A[28+l_m] * B[10];
    C[16+l_m] += A[28+l_m] * B[11];
    C[16+l_m] += A[32+l_m] * B[12];
    C[20+l_m] += A[24+l_m] * B[13];
    C[20+l_m] += A[32+l_m] * B[14];
    C[24+l_m] += A[0+l_m] * B[15];
    C[24+l_m] += A[12+l_m] * B[16];
    C[24+l_m] += A[20+l_m] * B[17];
    C[28+l_m] += A[4+l_m] * B[18];
    C[28+l_m] += A[12+l_m] * B[19];
    C[28+l_m] += A[16+l_m] * B[20];
    C[32+l_m] += A[8+l_m] * B[21];
    C[32+l_m] += A[16+l_m] * B[22];
    C[32+l_m] += A[20+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 192;
#endif
}

void ssparse_starMatrix_m1_n9_k9_ldA4_ldBna6_ldC4_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  for ( l_m = 0; l_m < 1; l_m++) {
    C[0+l_m] += A[24+l_m] * B[0];
    C[0+l_m] += A[28+l_m] * B[1];
    C[0+l_m] += A[32+l_m] * B[2];
    C[4+l_m] += A[24+l_m] * B[3];
    C[4+l_m] += A[28+l_m] * B[4];
    C[4+l_m] += A[32+l_m] * B[5];
    C[8+l_m] += A[24+l_m] * B[6];
    C[8+l_m] += A[28+l_m] * B[7];
    C[8+l_m] += A[32+l_m] * B[8];
    C[12+l_m] += A[24+l_m] * B[9];
    C[12+l_m] += A[28+l_m] * B[10];
    C[16+l_m] += A[28+l_m] * B[11];
    C[16+l_m] += A[32+l_m] * B[12];
    C[20+l_m] += A[24+l_m] * B[13];
    C[20+l_m] += A[32+l_m] * B[14];
    C[24+l_m] += A[0+l_m] * B[15];
    C[24+l_m] += A[12+l_m] * B[16];
    C[24+l_m] += A[20+l_m] * B[17];
    C[28+l_m] += A[4+l_m] * B[18];
    C[28+l_m] += A[12+l_m] * B[19];
    C[28+l_m] += A[16+l_m] * B[20];
    C[32+l_m] += A[8+l_m] * B[21];
    C[32+l_m] += A[16+l_m] * B[22];
    C[32+l_m] += A[20+l_m] * B[23];
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

void ssparse_starMatrix_m35_n9_k9_ldA36_ldBna7_ldC36_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 35; l_m++) {
    C[0+l_m] += A[216+l_m] * B[0];
    C[0+l_m] += A[252+l_m] * B[1];
    C[0+l_m] += A[288+l_m] * B[2];
    C[36+l_m] += A[216+l_m] * B[3];
    C[36+l_m] += A[252+l_m] * B[4];
    C[36+l_m] += A[288+l_m] * B[5];
    C[72+l_m] += A[216+l_m] * B[6];
    C[72+l_m] += A[252+l_m] * B[7];
    C[72+l_m] += A[288+l_m] * B[8];
    C[108+l_m] += A[216+l_m] * B[9];
    C[108+l_m] += A[252+l_m] * B[10];
    C[144+l_m] += A[252+l_m] * B[11];
    C[144+l_m] += A[288+l_m] * B[12];
    C[180+l_m] += A[216+l_m] * B[13];
    C[180+l_m] += A[288+l_m] * B[14];
    C[216+l_m] += A[0+l_m] * B[15];
    C[216+l_m] += A[108+l_m] * B[16];
    C[216+l_m] += A[180+l_m] * B[17];
    C[252+l_m] += A[36+l_m] * B[18];
    C[252+l_m] += A[108+l_m] * B[19];
    C[252+l_m] += A[144+l_m] * B[20];
    C[288+l_m] += A[72+l_m] * B[21];
    C[288+l_m] += A[144+l_m] * B[22];
    C[288+l_m] += A[180+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 1680;
#endif
}

void ssparse_starMatrix_m20_n9_k9_ldA20_ldBna7_ldC20_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 20; l_m++) {
    C[0+l_m] += A[120+l_m] * B[0];
    C[0+l_m] += A[140+l_m] * B[1];
    C[0+l_m] += A[160+l_m] * B[2];
    C[20+l_m] += A[120+l_m] * B[3];
    C[20+l_m] += A[140+l_m] * B[4];
    C[20+l_m] += A[160+l_m] * B[5];
    C[40+l_m] += A[120+l_m] * B[6];
    C[40+l_m] += A[140+l_m] * B[7];
    C[40+l_m] += A[160+l_m] * B[8];
    C[60+l_m] += A[120+l_m] * B[9];
    C[60+l_m] += A[140+l_m] * B[10];
    C[80+l_m] += A[140+l_m] * B[11];
    C[80+l_m] += A[160+l_m] * B[12];
    C[100+l_m] += A[120+l_m] * B[13];
    C[100+l_m] += A[160+l_m] * B[14];
    C[120+l_m] += A[0+l_m] * B[15];
    C[120+l_m] += A[60+l_m] * B[16];
    C[120+l_m] += A[100+l_m] * B[17];
    C[140+l_m] += A[20+l_m] * B[18];
    C[140+l_m] += A[60+l_m] * B[19];
    C[140+l_m] += A[80+l_m] * B[20];
    C[160+l_m] += A[40+l_m] * B[21];
    C[160+l_m] += A[80+l_m] * B[22];
    C[160+l_m] += A[100+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 960;
#endif
}

void ssparse_starMatrix_m10_n9_k9_ldA12_ldBna7_ldC12_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 10; l_m++) {
    C[0+l_m] += A[72+l_m] * B[0];
    C[0+l_m] += A[84+l_m] * B[1];
    C[0+l_m] += A[96+l_m] * B[2];
    C[12+l_m] += A[72+l_m] * B[3];
    C[12+l_m] += A[84+l_m] * B[4];
    C[12+l_m] += A[96+l_m] * B[5];
    C[24+l_m] += A[72+l_m] * B[6];
    C[24+l_m] += A[84+l_m] * B[7];
    C[24+l_m] += A[96+l_m] * B[8];
    C[36+l_m] += A[72+l_m] * B[9];
    C[36+l_m] += A[84+l_m] * B[10];
    C[48+l_m] += A[84+l_m] * B[11];
    C[48+l_m] += A[96+l_m] * B[12];
    C[60+l_m] += A[72+l_m] * B[13];
    C[60+l_m] += A[96+l_m] * B[14];
    C[72+l_m] += A[0+l_m] * B[15];
    C[72+l_m] += A[36+l_m] * B[16];
    C[72+l_m] += A[60+l_m] * B[17];
    C[84+l_m] += A[12+l_m] * B[18];
    C[84+l_m] += A[36+l_m] * B[19];
    C[84+l_m] += A[48+l_m] * B[20];
    C[96+l_m] += A[24+l_m] * B[21];
    C[96+l_m] += A[48+l_m] * B[22];
    C[96+l_m] += A[60+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 480;
#endif
}

void ssparse_starMatrix_m4_n9_k9_ldA4_ldBna7_ldC4_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(4)
  #pragma vector aligned
  for ( l_m = 0; l_m < 4; l_m++) {
    C[0+l_m] += A[24+l_m] * B[0];
    C[0+l_m] += A[28+l_m] * B[1];
    C[0+l_m] += A[32+l_m] * B[2];
    C[4+l_m] += A[24+l_m] * B[3];
    C[4+l_m] += A[28+l_m] * B[4];
    C[4+l_m] += A[32+l_m] * B[5];
    C[8+l_m] += A[24+l_m] * B[6];
    C[8+l_m] += A[28+l_m] * B[7];
    C[8+l_m] += A[32+l_m] * B[8];
    C[12+l_m] += A[24+l_m] * B[9];
    C[12+l_m] += A[28+l_m] * B[10];
    C[16+l_m] += A[28+l_m] * B[11];
    C[16+l_m] += A[32+l_m] * B[12];
    C[20+l_m] += A[24+l_m] * B[13];
    C[20+l_m] += A[32+l_m] * B[14];
    C[24+l_m] += A[0+l_m] * B[15];
    C[24+l_m] += A[12+l_m] * B[16];
    C[24+l_m] += A[20+l_m] * B[17];
    C[28+l_m] += A[4+l_m] * B[18];
    C[28+l_m] += A[12+l_m] * B[19];
    C[28+l_m] += A[16+l_m] * B[20];
    C[32+l_m] += A[8+l_m] * B[21];
    C[32+l_m] += A[16+l_m] * B[22];
    C[32+l_m] += A[20+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 192;
#endif
}

void ssparse_starMatrix_m1_n9_k9_ldA4_ldBna7_ldC4_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  for ( l_m = 0; l_m < 1; l_m++) {
    C[0+l_m] += A[24+l_m] * B[0];
    C[0+l_m] += A[28+l_m] * B[1];
    C[0+l_m] += A[32+l_m] * B[2];
    C[4+l_m] += A[24+l_m] * B[3];
    C[4+l_m] += A[28+l_m] * B[4];
    C[4+l_m] += A[32+l_m] * B[5];
    C[8+l_m] += A[24+l_m] * B[6];
    C[8+l_m] += A[28+l_m] * B[7];
    C[8+l_m] += A[32+l_m] * B[8];
    C[12+l_m] += A[24+l_m] * B[9];
    C[12+l_m] += A[28+l_m] * B[10];
    C[16+l_m] += A[28+l_m] * B[11];
    C[16+l_m] += A[32+l_m] * B[12];
    C[20+l_m] += A[24+l_m] * B[13];
    C[20+l_m] += A[32+l_m] * B[14];
    C[24+l_m] += A[0+l_m] * B[15];
    C[24+l_m] += A[12+l_m] * B[16];
    C[24+l_m] += A[20+l_m] * B[17];
    C[28+l_m] += A[4+l_m] * B[18];
    C[28+l_m] += A[12+l_m] * B[19];
    C[28+l_m] += A[16+l_m] * B[20];
    C[32+l_m] += A[8+l_m] * B[21];
    C[32+l_m] += A[16+l_m] * B[22];
    C[32+l_m] += A[20+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 48;
#endif
}

void ssparse_starMatrix_m84_n9_k9_ldA84_ldBna8_ldC84_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 84; l_m++) {
    C[0+l_m] += A[504+l_m] * B[0];
    C[0+l_m] += A[588+l_m] * B[1];
    C[0+l_m] += A[672+l_m] * B[2];
    C[84+l_m] += A[504+l_m] * B[3];
    C[84+l_m] += A[588+l_m] * B[4];
    C[84+l_m] += A[672+l_m] * B[5];
    C[168+l_m] += A[504+l_m] * B[6];
    C[168+l_m] += A[588+l_m] * B[7];
    C[168+l_m] += A[672+l_m] * B[8];
    C[252+l_m] += A[504+l_m] * B[9];
    C[252+l_m] += A[588+l_m] * B[10];
    C[336+l_m] += A[588+l_m] * B[11];
    C[336+l_m] += A[672+l_m] * B[12];
    C[420+l_m] += A[504+l_m] * B[13];
    C[420+l_m] += A[672+l_m] * B[14];
    C[504+l_m] += A[0+l_m] * B[15];
    C[504+l_m] += A[252+l_m] * B[16];
    C[504+l_m] += A[420+l_m] * B[17];
    C[588+l_m] += A[84+l_m] * B[18];
    C[588+l_m] += A[252+l_m] * B[19];
    C[588+l_m] += A[336+l_m] * B[20];
    C[672+l_m] += A[168+l_m] * B[21];
    C[672+l_m] += A[336+l_m] * B[22];
    C[672+l_m] += A[420+l_m] * B[23];
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

void ssparse_starMatrix_m35_n9_k9_ldA36_ldBna8_ldC36_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 35; l_m++) {
    C[0+l_m] += A[216+l_m] * B[0];
    C[0+l_m] += A[252+l_m] * B[1];
    C[0+l_m] += A[288+l_m] * B[2];
    C[36+l_m] += A[216+l_m] * B[3];
    C[36+l_m] += A[252+l_m] * B[4];
    C[36+l_m] += A[288+l_m] * B[5];
    C[72+l_m] += A[216+l_m] * B[6];
    C[72+l_m] += A[252+l_m] * B[7];
    C[72+l_m] += A[288+l_m] * B[8];
    C[108+l_m] += A[216+l_m] * B[9];
    C[108+l_m] += A[252+l_m] * B[10];
    C[144+l_m] += A[252+l_m] * B[11];
    C[144+l_m] += A[288+l_m] * B[12];
    C[180+l_m] += A[216+l_m] * B[13];
    C[180+l_m] += A[288+l_m] * B[14];
    C[216+l_m] += A[0+l_m] * B[15];
    C[216+l_m] += A[108+l_m] * B[16];
    C[216+l_m] += A[180+l_m] * B[17];
    C[252+l_m] += A[36+l_m] * B[18];
    C[252+l_m] += A[108+l_m] * B[19];
    C[252+l_m] += A[144+l_m] * B[20];
    C[288+l_m] += A[72+l_m] * B[21];
    C[288+l_m] += A[144+l_m] * B[22];
    C[288+l_m] += A[180+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 1680;
#endif
}

void ssparse_starMatrix_m20_n9_k9_ldA20_ldBna8_ldC20_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 20; l_m++) {
    C[0+l_m] += A[120+l_m] * B[0];
    C[0+l_m] += A[140+l_m] * B[1];
    C[0+l_m] += A[160+l_m] * B[2];
    C[20+l_m] += A[120+l_m] * B[3];
    C[20+l_m] += A[140+l_m] * B[4];
    C[20+l_m] += A[160+l_m] * B[5];
    C[40+l_m] += A[120+l_m] * B[6];
    C[40+l_m] += A[140+l_m] * B[7];
    C[40+l_m] += A[160+l_m] * B[8];
    C[60+l_m] += A[120+l_m] * B[9];
    C[60+l_m] += A[140+l_m] * B[10];
    C[80+l_m] += A[140+l_m] * B[11];
    C[80+l_m] += A[160+l_m] * B[12];
    C[100+l_m] += A[120+l_m] * B[13];
    C[100+l_m] += A[160+l_m] * B[14];
    C[120+l_m] += A[0+l_m] * B[15];
    C[120+l_m] += A[60+l_m] * B[16];
    C[120+l_m] += A[100+l_m] * B[17];
    C[140+l_m] += A[20+l_m] * B[18];
    C[140+l_m] += A[60+l_m] * B[19];
    C[140+l_m] += A[80+l_m] * B[20];
    C[160+l_m] += A[40+l_m] * B[21];
    C[160+l_m] += A[80+l_m] * B[22];
    C[160+l_m] += A[100+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 960;
#endif
}

void ssparse_starMatrix_m10_n9_k9_ldA12_ldBna8_ldC12_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 10; l_m++) {
    C[0+l_m] += A[72+l_m] * B[0];
    C[0+l_m] += A[84+l_m] * B[1];
    C[0+l_m] += A[96+l_m] * B[2];
    C[12+l_m] += A[72+l_m] * B[3];
    C[12+l_m] += A[84+l_m] * B[4];
    C[12+l_m] += A[96+l_m] * B[5];
    C[24+l_m] += A[72+l_m] * B[6];
    C[24+l_m] += A[84+l_m] * B[7];
    C[24+l_m] += A[96+l_m] * B[8];
    C[36+l_m] += A[72+l_m] * B[9];
    C[36+l_m] += A[84+l_m] * B[10];
    C[48+l_m] += A[84+l_m] * B[11];
    C[48+l_m] += A[96+l_m] * B[12];
    C[60+l_m] += A[72+l_m] * B[13];
    C[60+l_m] += A[96+l_m] * B[14];
    C[72+l_m] += A[0+l_m] * B[15];
    C[72+l_m] += A[36+l_m] * B[16];
    C[72+l_m] += A[60+l_m] * B[17];
    C[84+l_m] += A[12+l_m] * B[18];
    C[84+l_m] += A[36+l_m] * B[19];
    C[84+l_m] += A[48+l_m] * B[20];
    C[96+l_m] += A[24+l_m] * B[21];
    C[96+l_m] += A[48+l_m] * B[22];
    C[96+l_m] += A[60+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 480;
#endif
}

void ssparse_starMatrix_m4_n9_k9_ldA4_ldBna8_ldC4_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(4)
  #pragma vector aligned
  for ( l_m = 0; l_m < 4; l_m++) {
    C[0+l_m] += A[24+l_m] * B[0];
    C[0+l_m] += A[28+l_m] * B[1];
    C[0+l_m] += A[32+l_m] * B[2];
    C[4+l_m] += A[24+l_m] * B[3];
    C[4+l_m] += A[28+l_m] * B[4];
    C[4+l_m] += A[32+l_m] * B[5];
    C[8+l_m] += A[24+l_m] * B[6];
    C[8+l_m] += A[28+l_m] * B[7];
    C[8+l_m] += A[32+l_m] * B[8];
    C[12+l_m] += A[24+l_m] * B[9];
    C[12+l_m] += A[28+l_m] * B[10];
    C[16+l_m] += A[28+l_m] * B[11];
    C[16+l_m] += A[32+l_m] * B[12];
    C[20+l_m] += A[24+l_m] * B[13];
    C[20+l_m] += A[32+l_m] * B[14];
    C[24+l_m] += A[0+l_m] * B[15];
    C[24+l_m] += A[12+l_m] * B[16];
    C[24+l_m] += A[20+l_m] * B[17];
    C[28+l_m] += A[4+l_m] * B[18];
    C[28+l_m] += A[12+l_m] * B[19];
    C[28+l_m] += A[16+l_m] * B[20];
    C[32+l_m] += A[8+l_m] * B[21];
    C[32+l_m] += A[16+l_m] * B[22];
    C[32+l_m] += A[20+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 192;
#endif
}

void ssparse_starMatrix_m1_n9_k9_ldA4_ldBna8_ldC4_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  for ( l_m = 0; l_m < 1; l_m++) {
    C[0+l_m] += A[24+l_m] * B[0];
    C[0+l_m] += A[28+l_m] * B[1];
    C[0+l_m] += A[32+l_m] * B[2];
    C[4+l_m] += A[24+l_m] * B[3];
    C[4+l_m] += A[28+l_m] * B[4];
    C[4+l_m] += A[32+l_m] * B[5];
    C[8+l_m] += A[24+l_m] * B[6];
    C[8+l_m] += A[28+l_m] * B[7];
    C[8+l_m] += A[32+l_m] * B[8];
    C[12+l_m] += A[24+l_m] * B[9];
    C[12+l_m] += A[28+l_m] * B[10];
    C[16+l_m] += A[28+l_m] * B[11];
    C[16+l_m] += A[32+l_m] * B[12];
    C[20+l_m] += A[24+l_m] * B[13];
    C[20+l_m] += A[32+l_m] * B[14];
    C[24+l_m] += A[0+l_m] * B[15];
    C[24+l_m] += A[12+l_m] * B[16];
    C[24+l_m] += A[20+l_m] * B[17];
    C[28+l_m] += A[4+l_m] * B[18];
    C[28+l_m] += A[12+l_m] * B[19];
    C[28+l_m] += A[16+l_m] * B[20];
    C[32+l_m] += A[8+l_m] * B[21];
    C[32+l_m] += A[16+l_m] * B[22];
    C[32+l_m] += A[20+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 48;
#endif
}

void ssparse_starMatrix_m4_n9_k9_ldA4_ldBna2_ldC4_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(4)
  #pragma vector aligned
  for ( l_m = 0; l_m < 4; l_m++) {
    C[0+l_m] += A[24+l_m] * B[0];
    C[0+l_m] += A[28+l_m] * B[1];
    C[0+l_m] += A[32+l_m] * B[2];
    C[4+l_m] += A[24+l_m] * B[3];
    C[4+l_m] += A[28+l_m] * B[4];
    C[4+l_m] += A[32+l_m] * B[5];
    C[8+l_m] += A[24+l_m] * B[6];
    C[8+l_m] += A[28+l_m] * B[7];
    C[8+l_m] += A[32+l_m] * B[8];
    C[12+l_m] += A[24+l_m] * B[9];
    C[12+l_m] += A[28+l_m] * B[10];
    C[16+l_m] += A[28+l_m] * B[11];
    C[16+l_m] += A[32+l_m] * B[12];
    C[20+l_m] += A[24+l_m] * B[13];
    C[20+l_m] += A[32+l_m] * B[14];
    C[24+l_m] += A[0+l_m] * B[15];
    C[24+l_m] += A[12+l_m] * B[16];
    C[24+l_m] += A[20+l_m] * B[17];
    C[28+l_m] += A[4+l_m] * B[18];
    C[28+l_m] += A[12+l_m] * B[19];
    C[28+l_m] += A[16+l_m] * B[20];
    C[32+l_m] += A[8+l_m] * B[21];
    C[32+l_m] += A[16+l_m] * B[22];
    C[32+l_m] += A[20+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 192;
#endif
}

void ssparse_starMatrix_m10_n9_k9_ldA12_ldBna3_ldC12_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 10; l_m++) {
    C[0+l_m] += A[72+l_m] * B[0];
    C[0+l_m] += A[84+l_m] * B[1];
    C[0+l_m] += A[96+l_m] * B[2];
    C[12+l_m] += A[72+l_m] * B[3];
    C[12+l_m] += A[84+l_m] * B[4];
    C[12+l_m] += A[96+l_m] * B[5];
    C[24+l_m] += A[72+l_m] * B[6];
    C[24+l_m] += A[84+l_m] * B[7];
    C[24+l_m] += A[96+l_m] * B[8];
    C[36+l_m] += A[72+l_m] * B[9];
    C[36+l_m] += A[84+l_m] * B[10];
    C[48+l_m] += A[84+l_m] * B[11];
    C[48+l_m] += A[96+l_m] * B[12];
    C[60+l_m] += A[72+l_m] * B[13];
    C[60+l_m] += A[96+l_m] * B[14];
    C[72+l_m] += A[0+l_m] * B[15];
    C[72+l_m] += A[36+l_m] * B[16];
    C[72+l_m] += A[60+l_m] * B[17];
    C[84+l_m] += A[12+l_m] * B[18];
    C[84+l_m] += A[36+l_m] * B[19];
    C[84+l_m] += A[48+l_m] * B[20];
    C[96+l_m] += A[24+l_m] * B[21];
    C[96+l_m] += A[48+l_m] * B[22];
    C[96+l_m] += A[60+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 480;
#endif
}

void ssparse_starMatrix_m20_n9_k9_ldA20_ldBna4_ldC20_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 20; l_m++) {
    C[0+l_m] += A[120+l_m] * B[0];
    C[0+l_m] += A[140+l_m] * B[1];
    C[0+l_m] += A[160+l_m] * B[2];
    C[20+l_m] += A[120+l_m] * B[3];
    C[20+l_m] += A[140+l_m] * B[4];
    C[20+l_m] += A[160+l_m] * B[5];
    C[40+l_m] += A[120+l_m] * B[6];
    C[40+l_m] += A[140+l_m] * B[7];
    C[40+l_m] += A[160+l_m] * B[8];
    C[60+l_m] += A[120+l_m] * B[9];
    C[60+l_m] += A[140+l_m] * B[10];
    C[80+l_m] += A[140+l_m] * B[11];
    C[80+l_m] += A[160+l_m] * B[12];
    C[100+l_m] += A[120+l_m] * B[13];
    C[100+l_m] += A[160+l_m] * B[14];
    C[120+l_m] += A[0+l_m] * B[15];
    C[120+l_m] += A[60+l_m] * B[16];
    C[120+l_m] += A[100+l_m] * B[17];
    C[140+l_m] += A[20+l_m] * B[18];
    C[140+l_m] += A[60+l_m] * B[19];
    C[140+l_m] += A[80+l_m] * B[20];
    C[160+l_m] += A[40+l_m] * B[21];
    C[160+l_m] += A[80+l_m] * B[22];
    C[160+l_m] += A[100+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 960;
#endif
}

void ssparse_starMatrix_m35_n9_k9_ldA36_ldBna5_ldC36_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 35; l_m++) {
    C[0+l_m] += A[216+l_m] * B[0];
    C[0+l_m] += A[252+l_m] * B[1];
    C[0+l_m] += A[288+l_m] * B[2];
    C[36+l_m] += A[216+l_m] * B[3];
    C[36+l_m] += A[252+l_m] * B[4];
    C[36+l_m] += A[288+l_m] * B[5];
    C[72+l_m] += A[216+l_m] * B[6];
    C[72+l_m] += A[252+l_m] * B[7];
    C[72+l_m] += A[288+l_m] * B[8];
    C[108+l_m] += A[216+l_m] * B[9];
    C[108+l_m] += A[252+l_m] * B[10];
    C[144+l_m] += A[252+l_m] * B[11];
    C[144+l_m] += A[288+l_m] * B[12];
    C[180+l_m] += A[216+l_m] * B[13];
    C[180+l_m] += A[288+l_m] * B[14];
    C[216+l_m] += A[0+l_m] * B[15];
    C[216+l_m] += A[108+l_m] * B[16];
    C[216+l_m] += A[180+l_m] * B[17];
    C[252+l_m] += A[36+l_m] * B[18];
    C[252+l_m] += A[108+l_m] * B[19];
    C[252+l_m] += A[144+l_m] * B[20];
    C[288+l_m] += A[72+l_m] * B[21];
    C[288+l_m] += A[144+l_m] * B[22];
    C[288+l_m] += A[180+l_m] * B[23];
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

void ssparse_starMatrix_m84_n9_k9_ldA84_ldBna7_ldC84_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch, const float* B_prefetch, const float* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 84; l_m++) {
    C[0+l_m] += A[504+l_m] * B[0];
    C[0+l_m] += A[588+l_m] * B[1];
    C[0+l_m] += A[672+l_m] * B[2];
    C[84+l_m] += A[504+l_m] * B[3];
    C[84+l_m] += A[588+l_m] * B[4];
    C[84+l_m] += A[672+l_m] * B[5];
    C[168+l_m] += A[504+l_m] * B[6];
    C[168+l_m] += A[588+l_m] * B[7];
    C[168+l_m] += A[672+l_m] * B[8];
    C[252+l_m] += A[504+l_m] * B[9];
    C[252+l_m] += A[588+l_m] * B[10];
    C[336+l_m] += A[588+l_m] * B[11];
    C[336+l_m] += A[672+l_m] * B[12];
    C[420+l_m] += A[504+l_m] * B[13];
    C[420+l_m] += A[672+l_m] * B[14];
    C[504+l_m] += A[0+l_m] * B[15];
    C[504+l_m] += A[252+l_m] * B[16];
    C[504+l_m] += A[420+l_m] * B[17];
    C[588+l_m] += A[84+l_m] * B[18];
    C[588+l_m] += A[252+l_m] * B[19];
    C[588+l_m] += A[336+l_m] * B[20];
    C[672+l_m] += A[168+l_m] * B[21];
    C[672+l_m] += A[336+l_m] * B[22];
    C[672+l_m] += A[420+l_m] * B[23];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 4032;
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
