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
// @date 2015-09-27 13:25:38.148859
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
#ifndef SPARSEDHSWCPP
#define SPARSEDHSWCPP

#if defined( __SSE3__) || defined(__MIC__)
#include <immintrin.h>
#endif

#include <cstddef>
#ifndef NDEBUG
extern long long libxsmm_num_total_flops;
#endif

void dsparse_starMatrix_m1_n9_k9_ldA4_ldBna2_ldC4_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m4_n9_k9_ldA4_ldBna3_ldC4_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m1_n9_k9_ldA4_ldBna3_ldC4_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m10_n9_k9_ldA12_ldBna4_ldC12_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m4_n9_k9_ldA4_ldBna4_ldC4_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m1_n9_k9_ldA4_ldBna4_ldC4_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m20_n9_k9_ldA20_ldBna5_ldC20_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m10_n9_k9_ldA12_ldBna5_ldC12_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m4_n9_k9_ldA4_ldBna5_ldC4_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m1_n9_k9_ldA4_ldBna5_ldC4_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m35_n9_k9_ldA36_ldBna6_ldC36_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m20_n9_k9_ldA20_ldBna6_ldC20_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m10_n9_k9_ldA12_ldBna6_ldC12_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m4_n9_k9_ldA4_ldBna6_ldC4_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m1_n9_k9_ldA4_ldBna6_ldC4_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m56_n9_k9_ldA56_ldBna7_ldC56_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m35_n9_k9_ldA36_ldBna7_ldC36_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m20_n9_k9_ldA20_ldBna7_ldC20_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m10_n9_k9_ldA12_ldBna7_ldC12_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m4_n9_k9_ldA4_ldBna7_ldC4_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m1_n9_k9_ldA4_ldBna7_ldC4_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m84_n9_k9_ldA84_ldBna8_ldC84_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m56_n9_k9_ldA56_ldBna8_ldC56_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m35_n9_k9_ldA36_ldBna8_ldC36_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m20_n9_k9_ldA20_ldBna8_ldC20_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m10_n9_k9_ldA12_ldBna8_ldC12_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m4_n9_k9_ldA4_ldBna8_ldC4_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m1_n9_k9_ldA4_ldBna8_ldC4_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m4_n9_k9_ldA4_ldBna2_ldC4_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m10_n9_k9_ldA12_ldBna3_ldC12_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_fP113DivM_m10_n9_k10_ldAna3_ldB12_ldC12_beta0_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_n = 0;
  #pragma nounroll_and_jam
  for ( l_n = 0; l_n < 9; l_n++) {
    unsigned int l_m = 0;
   #pragma simd
    for ( l_m = 0; l_m < 10; l_m++) {
      C[(l_n*12)+l_m] = 0.0;
    }
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b0 = _mm256_broadcast_sd(&B[(l_n*12)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b0 = _mm_loaddup_pd(&B[(l_n*12)+0]);
#endif
    __m128d c0_0 = _mm_load_sd(&C[(l_n*12)+0]);
    __m128d a0_0 = _mm_load_sd(&A[0]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, b0));
#endif
    _mm_store_sd(&C[(l_n*12)+0], c0_0);
    __m128d c0_1 = _mm_load_sd(&C[(l_n*12)+3]);
    __m128d a0_1 = _mm_load_sd(&A[1]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, b0));
#endif
    _mm_store_sd(&C[(l_n*12)+3], c0_1);
    __m128d c0_2 = _mm_load_sd(&C[(l_n*12)+9]);
    __m128d a0_2 = _mm_load_sd(&A[2]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, b0));
#endif
    _mm_store_sd(&C[(l_n*12)+9], c0_2);
#else
    C[(l_n*12)+0] += A[0] * B[(l_n*12)+0];
    C[(l_n*12)+3] += A[1] * B[(l_n*12)+0];
    C[(l_n*12)+9] += A[2] * B[(l_n*12)+0];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b1 = _mm256_broadcast_sd(&B[(l_n*12)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b1 = _mm_loaddup_pd(&B[(l_n*12)+1]);
#endif
    __m128d c1_0 = _mm_load_sd(&C[(l_n*12)+1]);
    __m128d a1_0 = _mm_load_sd(&A[3]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, b1));
#endif
    _mm_store_sd(&C[(l_n*12)+1], c1_0);
    __m128d c1_1 = _mm_load_sd(&C[(l_n*12)+7]);
    __m128d a1_1 = _mm_load_sd(&A[4]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_1 = _mm_add_sd(c1_1, _mm_mul_sd(a1_1, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_1 = _mm_add_sd(c1_1, _mm_mul_sd(a1_1, b1));
#endif
    _mm_store_sd(&C[(l_n*12)+7], c1_1);
#else
    C[(l_n*12)+1] += A[3] * B[(l_n*12)+1];
    C[(l_n*12)+7] += A[4] * B[(l_n*12)+1];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b2 = _mm256_broadcast_sd(&B[(l_n*12)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b2 = _mm_loaddup_pd(&B[(l_n*12)+2]);
#endif
    __m128d c2_0 = _mm_load_sd(&C[(l_n*12)+2]);
    __m128d a2_0 = _mm_load_sd(&A[5]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_0 = _mm_add_sd(c2_0, _mm_mul_sd(a2_0, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_0 = _mm_add_sd(c2_0, _mm_mul_sd(a2_0, b2));
#endif
    _mm_store_sd(&C[(l_n*12)+2], c2_0);
    __m128d c2_1 = _mm_load_sd(&C[(l_n*12)+8]);
    __m128d a2_1 = _mm_load_sd(&A[6]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_1 = _mm_add_sd(c2_1, _mm_mul_sd(a2_1, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_1 = _mm_add_sd(c2_1, _mm_mul_sd(a2_1, b2));
#endif
    _mm_store_sd(&C[(l_n*12)+8], c2_1);
#else
    C[(l_n*12)+2] += A[5] * B[(l_n*12)+2];
    C[(l_n*12)+8] += A[6] * B[(l_n*12)+2];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b3 = _mm256_broadcast_sd(&B[(l_n*12)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b3 = _mm_loaddup_pd(&B[(l_n*12)+3]);
#endif
    __m128d c3_0 = _mm_load_sd(&C[(l_n*12)+0]);
    __m128d a3_0 = _mm_load_sd(&A[7]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, b3));
#endif
    _mm_store_sd(&C[(l_n*12)+0], c3_0);
    __m128d c3_1 = _mm_load_sd(&C[(l_n*12)+3]);
    __m128d a3_1 = _mm_load_sd(&A[8]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, b3));
#endif
    _mm_store_sd(&C[(l_n*12)+3], c3_1);
    __m128d c3_2 = _mm_load_sd(&C[(l_n*12)+9]);
    __m128d a3_2 = _mm_load_sd(&A[9]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, b3));
#endif
    _mm_store_sd(&C[(l_n*12)+9], c3_2);
#else
    C[(l_n*12)+0] += A[7] * B[(l_n*12)+3];
    C[(l_n*12)+3] += A[8] * B[(l_n*12)+3];
    C[(l_n*12)+9] += A[9] * B[(l_n*12)+3];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b4 = _mm256_broadcast_sd(&B[(l_n*12)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b4 = _mm_loaddup_pd(&B[(l_n*12)+4]);
#endif
    __m128d c4_0 = _mm_load_sd(&C[(l_n*12)+4]);
    __m128d a4_0 = _mm_load_sd(&A[10]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_0 = _mm_add_sd(c4_0, _mm_mul_sd(a4_0, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_0 = _mm_add_sd(c4_0, _mm_mul_sd(a4_0, b4));
#endif
    _mm_store_sd(&C[(l_n*12)+4], c4_0);
#else
    C[(l_n*12)+4] += A[10] * B[(l_n*12)+4];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b5 = _mm256_broadcast_sd(&B[(l_n*12)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b5 = _mm_loaddup_pd(&B[(l_n*12)+5]);
#endif
    __m128d c5_0 = _mm_load_sd(&C[(l_n*12)+5]);
    __m128d a5_0 = _mm_load_sd(&A[11]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_0 = _mm_add_sd(c5_0, _mm_mul_sd(a5_0, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_0 = _mm_add_sd(c5_0, _mm_mul_sd(a5_0, b5));
#endif
    _mm_store_sd(&C[(l_n*12)+5], c5_0);
#else
    C[(l_n*12)+5] += A[11] * B[(l_n*12)+5];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b6 = _mm256_broadcast_sd(&B[(l_n*12)+6]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b6 = _mm_loaddup_pd(&B[(l_n*12)+6]);
#endif
    __m128d c6_0 = _mm_load_sd(&C[(l_n*12)+6]);
    __m128d a6_0 = _mm_load_sd(&A[12]);
#if defined(__SSE3__) && defined(__AVX__)
    c6_0 = _mm_add_sd(c6_0, _mm_mul_sd(a6_0, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c6_0 = _mm_add_sd(c6_0, _mm_mul_sd(a6_0, b6));
#endif
    _mm_store_sd(&C[(l_n*12)+6], c6_0);
#else
    C[(l_n*12)+6] += A[12] * B[(l_n*12)+6];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b7 = _mm256_broadcast_sd(&B[(l_n*12)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b7 = _mm_loaddup_pd(&B[(l_n*12)+7]);
#endif
    __m128d c7_0 = _mm_load_sd(&C[(l_n*12)+1]);
    __m128d a7_0 = _mm_load_sd(&A[13]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_0 = _mm_add_sd(c7_0, _mm_mul_sd(a7_0, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_0 = _mm_add_sd(c7_0, _mm_mul_sd(a7_0, b7));
#endif
    _mm_store_sd(&C[(l_n*12)+1], c7_0);
    __m128d c7_1 = _mm_load_sd(&C[(l_n*12)+7]);
    __m128d a7_1 = _mm_load_sd(&A[14]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_1 = _mm_add_sd(c7_1, _mm_mul_sd(a7_1, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_1 = _mm_add_sd(c7_1, _mm_mul_sd(a7_1, b7));
#endif
    _mm_store_sd(&C[(l_n*12)+7], c7_1);
#else
    C[(l_n*12)+1] += A[13] * B[(l_n*12)+7];
    C[(l_n*12)+7] += A[14] * B[(l_n*12)+7];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b8 = _mm256_broadcast_sd(&B[(l_n*12)+8]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b8 = _mm_loaddup_pd(&B[(l_n*12)+8]);
#endif
    __m128d c8_0 = _mm_load_sd(&C[(l_n*12)+2]);
    __m128d a8_0 = _mm_load_sd(&A[15]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_0 = _mm_add_sd(c8_0, _mm_mul_sd(a8_0, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_0 = _mm_add_sd(c8_0, _mm_mul_sd(a8_0, b8));
#endif
    _mm_store_sd(&C[(l_n*12)+2], c8_0);
    __m128d c8_1 = _mm_load_sd(&C[(l_n*12)+8]);
    __m128d a8_1 = _mm_load_sd(&A[16]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_1 = _mm_add_sd(c8_1, _mm_mul_sd(a8_1, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_1 = _mm_add_sd(c8_1, _mm_mul_sd(a8_1, b8));
#endif
    _mm_store_sd(&C[(l_n*12)+8], c8_1);
#else
    C[(l_n*12)+2] += A[15] * B[(l_n*12)+8];
    C[(l_n*12)+8] += A[16] * B[(l_n*12)+8];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b9 = _mm256_broadcast_sd(&B[(l_n*12)+9]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b9 = _mm_loaddup_pd(&B[(l_n*12)+9]);
#endif
    __m128d c9_0 = _mm_load_sd(&C[(l_n*12)+0]);
    __m128d a9_0 = _mm_load_sd(&A[17]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, b9));
#endif
    _mm_store_sd(&C[(l_n*12)+0], c9_0);
    __m128d c9_1 = _mm_load_sd(&C[(l_n*12)+3]);
    __m128d a9_1 = _mm_load_sd(&A[18]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, b9));
#endif
    _mm_store_sd(&C[(l_n*12)+3], c9_1);
    __m128d c9_2 = _mm_load_sd(&C[(l_n*12)+9]);
    __m128d a9_2 = _mm_load_sd(&A[19]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, b9));
#endif
    _mm_store_sd(&C[(l_n*12)+9], c9_2);
#else
    C[(l_n*12)+0] += A[17] * B[(l_n*12)+9];
    C[(l_n*12)+3] += A[18] * B[(l_n*12)+9];
    C[(l_n*12)+9] += A[19] * B[(l_n*12)+9];
#endif

  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 360;
#endif
}

void dsparse_fP111DivM_m10_n9_k10_ldAna3_ldB12_ldC12_beta0_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_n = 0;
  #pragma nounroll_and_jam
  for ( l_n = 0; l_n < 9; l_n++) {
    unsigned int l_m = 0;
   #pragma simd
    for ( l_m = 0; l_m < 10; l_m++) {
      C[(l_n*12)+l_m] = 0.0;
    }
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b0 = _mm256_broadcast_sd(&B[(l_n*12)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b0 = _mm_loaddup_pd(&B[(l_n*12)+0]);
#endif
    __m128d c0_0 = _mm_load_sd(&C[(l_n*12)+0]);
    __m128d a0_0 = _mm_load_sd(&A[0]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, b0));
#endif
    _mm_store_sd(&C[(l_n*12)+0], c0_0);
    __m128d c0_1 = _mm_load_sd(&C[(l_n*12)+3]);
    __m128d a0_1 = _mm_load_sd(&A[1]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, b0));
#endif
    _mm_store_sd(&C[(l_n*12)+3], c0_1);
    __m128d c0_2 = _mm_load_sd(&C[(l_n*12)+9]);
    __m128d a0_2 = _mm_load_sd(&A[2]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, b0));
#endif
    _mm_store_sd(&C[(l_n*12)+9], c0_2);
#else
    C[(l_n*12)+0] += A[0] * B[(l_n*12)+0];
    C[(l_n*12)+3] += A[1] * B[(l_n*12)+0];
    C[(l_n*12)+9] += A[2] * B[(l_n*12)+0];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b1 = _mm256_broadcast_sd(&B[(l_n*12)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b1 = _mm_loaddup_pd(&B[(l_n*12)+1]);
#endif
    __m128d c1_0 = _mm_loadu_pd(&C[(l_n*12)+1]);
    __m128d a1_0 = _mm_loadu_pd(&A[3]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_0 = _mm_add_pd(c1_0, _mm_mul_pd(a1_0, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_0 = _mm_add_pd(c1_0, _mm_mul_pd(a1_0, b1));
#endif
    _mm_storeu_pd(&C[(l_n*12)+1], c1_0);
    __m128d c1_2 = _mm_loadu_pd(&C[(l_n*12)+7]);
    __m128d a1_2 = _mm_loadu_pd(&A[5]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_2 = _mm_add_pd(c1_2, _mm_mul_pd(a1_2, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_2 = _mm_add_pd(c1_2, _mm_mul_pd(a1_2, b1));
#endif
    _mm_storeu_pd(&C[(l_n*12)+7], c1_2);
#else
    C[(l_n*12)+1] += A[3] * B[(l_n*12)+1];
    C[(l_n*12)+2] += A[4] * B[(l_n*12)+1];
    C[(l_n*12)+7] += A[5] * B[(l_n*12)+1];
    C[(l_n*12)+8] += A[6] * B[(l_n*12)+1];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b2 = _mm256_broadcast_sd(&B[(l_n*12)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b2 = _mm_loaddup_pd(&B[(l_n*12)+2]);
#endif
    __m128d c2_0 = _mm_loadu_pd(&C[(l_n*12)+1]);
    __m128d a2_0 = _mm_loadu_pd(&A[7]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_0 = _mm_add_pd(c2_0, _mm_mul_pd(a2_0, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_0 = _mm_add_pd(c2_0, _mm_mul_pd(a2_0, b2));
#endif
    _mm_storeu_pd(&C[(l_n*12)+1], c2_0);
    __m128d c2_2 = _mm_loadu_pd(&C[(l_n*12)+7]);
    __m128d a2_2 = _mm_loadu_pd(&A[9]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_2 = _mm_add_pd(c2_2, _mm_mul_pd(a2_2, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_2 = _mm_add_pd(c2_2, _mm_mul_pd(a2_2, b2));
#endif
    _mm_storeu_pd(&C[(l_n*12)+7], c2_2);
#else
    C[(l_n*12)+1] += A[7] * B[(l_n*12)+2];
    C[(l_n*12)+2] += A[8] * B[(l_n*12)+2];
    C[(l_n*12)+7] += A[9] * B[(l_n*12)+2];
    C[(l_n*12)+8] += A[10] * B[(l_n*12)+2];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b3 = _mm256_broadcast_sd(&B[(l_n*12)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b3 = _mm_loaddup_pd(&B[(l_n*12)+3]);
#endif
    __m128d c3_0 = _mm_load_sd(&C[(l_n*12)+0]);
    __m128d a3_0 = _mm_load_sd(&A[11]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, b3));
#endif
    _mm_store_sd(&C[(l_n*12)+0], c3_0);
    __m128d c3_1 = _mm_load_sd(&C[(l_n*12)+3]);
    __m128d a3_1 = _mm_load_sd(&A[12]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, b3));
#endif
    _mm_store_sd(&C[(l_n*12)+3], c3_1);
    __m128d c3_2 = _mm_load_sd(&C[(l_n*12)+9]);
    __m128d a3_2 = _mm_load_sd(&A[13]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, b3));
#endif
    _mm_store_sd(&C[(l_n*12)+9], c3_2);
#else
    C[(l_n*12)+0] += A[11] * B[(l_n*12)+3];
    C[(l_n*12)+3] += A[12] * B[(l_n*12)+3];
    C[(l_n*12)+9] += A[13] * B[(l_n*12)+3];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b4 = _mm256_broadcast_sd(&B[(l_n*12)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b4 = _mm_loaddup_pd(&B[(l_n*12)+4]);
#endif
    __m128d c4_0 = _mm_loadu_pd(&C[(l_n*12)+4]);
    __m128d a4_0 = _mm_loadu_pd(&A[14]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_0 = _mm_add_pd(c4_0, _mm_mul_pd(a4_0, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_0 = _mm_add_pd(c4_0, _mm_mul_pd(a4_0, b4));
#endif
    _mm_storeu_pd(&C[(l_n*12)+4], c4_0);
    __m128d c4_2 = _mm_load_sd(&C[(l_n*12)+6]);
    __m128d a4_2 = _mm_load_sd(&A[16]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, b4));
#endif
    _mm_store_sd(&C[(l_n*12)+6], c4_2);
#else
    C[(l_n*12)+4] += A[14] * B[(l_n*12)+4];
    C[(l_n*12)+5] += A[15] * B[(l_n*12)+4];
    C[(l_n*12)+6] += A[16] * B[(l_n*12)+4];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b5 = _mm256_broadcast_sd(&B[(l_n*12)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b5 = _mm_loaddup_pd(&B[(l_n*12)+5]);
#endif
    __m128d c5_0 = _mm_loadu_pd(&C[(l_n*12)+4]);
    __m128d a5_0 = _mm_loadu_pd(&A[17]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_0 = _mm_add_pd(c5_0, _mm_mul_pd(a5_0, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_0 = _mm_add_pd(c5_0, _mm_mul_pd(a5_0, b5));
#endif
    _mm_storeu_pd(&C[(l_n*12)+4], c5_0);
    __m128d c5_2 = _mm_load_sd(&C[(l_n*12)+6]);
    __m128d a5_2 = _mm_load_sd(&A[19]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, b5));
#endif
    _mm_store_sd(&C[(l_n*12)+6], c5_2);
#else
    C[(l_n*12)+4] += A[17] * B[(l_n*12)+5];
    C[(l_n*12)+5] += A[18] * B[(l_n*12)+5];
    C[(l_n*12)+6] += A[19] * B[(l_n*12)+5];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b6 = _mm256_broadcast_sd(&B[(l_n*12)+6]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b6 = _mm_loaddup_pd(&B[(l_n*12)+6]);
#endif
    __m128d c6_0 = _mm_loadu_pd(&C[(l_n*12)+4]);
    __m128d a6_0 = _mm_loadu_pd(&A[20]);
#if defined(__SSE3__) && defined(__AVX__)
    c6_0 = _mm_add_pd(c6_0, _mm_mul_pd(a6_0, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c6_0 = _mm_add_pd(c6_0, _mm_mul_pd(a6_0, b6));
#endif
    _mm_storeu_pd(&C[(l_n*12)+4], c6_0);
    __m128d c6_2 = _mm_load_sd(&C[(l_n*12)+6]);
    __m128d a6_2 = _mm_load_sd(&A[22]);
#if defined(__SSE3__) && defined(__AVX__)
    c6_2 = _mm_add_sd(c6_2, _mm_mul_sd(a6_2, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c6_2 = _mm_add_sd(c6_2, _mm_mul_sd(a6_2, b6));
#endif
    _mm_store_sd(&C[(l_n*12)+6], c6_2);
#else
    C[(l_n*12)+4] += A[20] * B[(l_n*12)+6];
    C[(l_n*12)+5] += A[21] * B[(l_n*12)+6];
    C[(l_n*12)+6] += A[22] * B[(l_n*12)+6];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b7 = _mm256_broadcast_sd(&B[(l_n*12)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b7 = _mm_loaddup_pd(&B[(l_n*12)+7]);
#endif
    __m128d c7_0 = _mm_loadu_pd(&C[(l_n*12)+1]);
    __m128d a7_0 = _mm_loadu_pd(&A[23]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_0 = _mm_add_pd(c7_0, _mm_mul_pd(a7_0, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_0 = _mm_add_pd(c7_0, _mm_mul_pd(a7_0, b7));
#endif
    _mm_storeu_pd(&C[(l_n*12)+1], c7_0);
    __m128d c7_2 = _mm_loadu_pd(&C[(l_n*12)+7]);
    __m128d a7_2 = _mm_loadu_pd(&A[25]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_2 = _mm_add_pd(c7_2, _mm_mul_pd(a7_2, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_2 = _mm_add_pd(c7_2, _mm_mul_pd(a7_2, b7));
#endif
    _mm_storeu_pd(&C[(l_n*12)+7], c7_2);
#else
    C[(l_n*12)+1] += A[23] * B[(l_n*12)+7];
    C[(l_n*12)+2] += A[24] * B[(l_n*12)+7];
    C[(l_n*12)+7] += A[25] * B[(l_n*12)+7];
    C[(l_n*12)+8] += A[26] * B[(l_n*12)+7];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b8 = _mm256_broadcast_sd(&B[(l_n*12)+8]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b8 = _mm_loaddup_pd(&B[(l_n*12)+8]);
#endif
    __m128d c8_0 = _mm_loadu_pd(&C[(l_n*12)+1]);
    __m128d a8_0 = _mm_loadu_pd(&A[27]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_0 = _mm_add_pd(c8_0, _mm_mul_pd(a8_0, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_0 = _mm_add_pd(c8_0, _mm_mul_pd(a8_0, b8));
#endif
    _mm_storeu_pd(&C[(l_n*12)+1], c8_0);
    __m128d c8_2 = _mm_loadu_pd(&C[(l_n*12)+7]);
    __m128d a8_2 = _mm_loadu_pd(&A[29]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_2 = _mm_add_pd(c8_2, _mm_mul_pd(a8_2, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_2 = _mm_add_pd(c8_2, _mm_mul_pd(a8_2, b8));
#endif
    _mm_storeu_pd(&C[(l_n*12)+7], c8_2);
#else
    C[(l_n*12)+1] += A[27] * B[(l_n*12)+8];
    C[(l_n*12)+2] += A[28] * B[(l_n*12)+8];
    C[(l_n*12)+7] += A[29] * B[(l_n*12)+8];
    C[(l_n*12)+8] += A[30] * B[(l_n*12)+8];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b9 = _mm256_broadcast_sd(&B[(l_n*12)+9]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b9 = _mm_loaddup_pd(&B[(l_n*12)+9]);
#endif
    __m128d c9_0 = _mm_load_sd(&C[(l_n*12)+0]);
    __m128d a9_0 = _mm_load_sd(&A[31]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, b9));
#endif
    _mm_store_sd(&C[(l_n*12)+0], c9_0);
    __m128d c9_1 = _mm_load_sd(&C[(l_n*12)+3]);
    __m128d a9_1 = _mm_load_sd(&A[32]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, b9));
#endif
    _mm_store_sd(&C[(l_n*12)+3], c9_1);
    __m128d c9_2 = _mm_load_sd(&C[(l_n*12)+9]);
    __m128d a9_2 = _mm_load_sd(&A[33]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, b9));
#endif
    _mm_store_sd(&C[(l_n*12)+9], c9_2);
#else
    C[(l_n*12)+0] += A[31] * B[(l_n*12)+9];
    C[(l_n*12)+3] += A[32] * B[(l_n*12)+9];
    C[(l_n*12)+9] += A[33] * B[(l_n*12)+9];
#endif

  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 612;
#endif
}

void dsparse_starMatrix_m20_n9_k9_ldA20_ldBna4_ldC20_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_starMatrix_m35_n9_k9_ldA36_ldBna5_ldC36_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_fM1DivM_m35_n9_k35_ldAna5_ldB36_ldC36_beta0_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_n = 0;
  #pragma nounroll_and_jam
  for ( l_n = 0; l_n < 9; l_n++) {
    unsigned int l_m = 0;
   #pragma simd
    for ( l_m = 0; l_m < 35; l_m++) {
      C[(l_n*36)+l_m] = 0.0;
    }
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b0 = _mm256_broadcast_sd(&B[(l_n*36)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b0 = _mm_loaddup_pd(&B[(l_n*36)+0]);
#endif
    __m128d c0_0 = _mm_load_sd(&C[(l_n*36)+0]);
    __m128d a0_0 = _mm_load_sd(&A[0]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, b0));
#endif
    _mm_store_sd(&C[(l_n*36)+0], c0_0);
    __m128d c0_1 = _mm_load_sd(&C[(l_n*36)+3]);
    __m128d a0_1 = _mm_load_sd(&A[1]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, b0));
#endif
    _mm_store_sd(&C[(l_n*36)+3], c0_1);
    __m128d c0_2 = _mm_load_sd(&C[(l_n*36)+9]);
    __m128d a0_2 = _mm_load_sd(&A[2]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, b0));
#endif
    _mm_store_sd(&C[(l_n*36)+9], c0_2);
    __m128d c0_3 = _mm_load_sd(&C[(l_n*36)+19]);
    __m128d a0_3 = _mm_load_sd(&A[3]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_3 = _mm_add_sd(c0_3, _mm_mul_sd(a0_3, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_3 = _mm_add_sd(c0_3, _mm_mul_sd(a0_3, b0));
#endif
    _mm_store_sd(&C[(l_n*36)+19], c0_3);
    __m128d c0_4 = _mm_load_sd(&C[(l_n*36)+34]);
    __m128d a0_4 = _mm_load_sd(&A[4]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_4 = _mm_add_sd(c0_4, _mm_mul_sd(a0_4, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_4 = _mm_add_sd(c0_4, _mm_mul_sd(a0_4, b0));
#endif
    _mm_store_sd(&C[(l_n*36)+34], c0_4);
#else
    C[(l_n*36)+0] += A[0] * B[(l_n*36)+0];
    C[(l_n*36)+3] += A[1] * B[(l_n*36)+0];
    C[(l_n*36)+9] += A[2] * B[(l_n*36)+0];
    C[(l_n*36)+19] += A[3] * B[(l_n*36)+0];
    C[(l_n*36)+34] += A[4] * B[(l_n*36)+0];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b1 = _mm256_broadcast_sd(&B[(l_n*36)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b1 = _mm_loaddup_pd(&B[(l_n*36)+1]);
#endif
    __m128d c1_0 = _mm_load_sd(&C[(l_n*36)+1]);
    __m128d a1_0 = _mm_load_sd(&A[5]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, b1));
#endif
    _mm_store_sd(&C[(l_n*36)+1], c1_0);
    __m128d c1_1 = _mm_load_sd(&C[(l_n*36)+7]);
    __m128d a1_1 = _mm_load_sd(&A[6]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_1 = _mm_add_sd(c1_1, _mm_mul_sd(a1_1, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_1 = _mm_add_sd(c1_1, _mm_mul_sd(a1_1, b1));
#endif
    _mm_store_sd(&C[(l_n*36)+7], c1_1);
    __m128d c1_2 = _mm_load_sd(&C[(l_n*36)+17]);
    __m128d a1_2 = _mm_load_sd(&A[7]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_2 = _mm_add_sd(c1_2, _mm_mul_sd(a1_2, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_2 = _mm_add_sd(c1_2, _mm_mul_sd(a1_2, b1));
#endif
    _mm_store_sd(&C[(l_n*36)+17], c1_2);
    __m128d c1_3 = _mm_load_sd(&C[(l_n*36)+32]);
    __m128d a1_3 = _mm_load_sd(&A[8]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_3 = _mm_add_sd(c1_3, _mm_mul_sd(a1_3, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_3 = _mm_add_sd(c1_3, _mm_mul_sd(a1_3, b1));
#endif
    _mm_store_sd(&C[(l_n*36)+32], c1_3);
#else
    C[(l_n*36)+1] += A[5] * B[(l_n*36)+1];
    C[(l_n*36)+7] += A[6] * B[(l_n*36)+1];
    C[(l_n*36)+17] += A[7] * B[(l_n*36)+1];
    C[(l_n*36)+32] += A[8] * B[(l_n*36)+1];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b2 = _mm256_broadcast_sd(&B[(l_n*36)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b2 = _mm_loaddup_pd(&B[(l_n*36)+2]);
#endif
    __m128d c2_0 = _mm_load_sd(&C[(l_n*36)+2]);
    __m128d a2_0 = _mm_load_sd(&A[9]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_0 = _mm_add_sd(c2_0, _mm_mul_sd(a2_0, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_0 = _mm_add_sd(c2_0, _mm_mul_sd(a2_0, b2));
#endif
    _mm_store_sd(&C[(l_n*36)+2], c2_0);
    __m128d c2_1 = _mm_load_sd(&C[(l_n*36)+8]);
    __m128d a2_1 = _mm_load_sd(&A[10]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_1 = _mm_add_sd(c2_1, _mm_mul_sd(a2_1, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_1 = _mm_add_sd(c2_1, _mm_mul_sd(a2_1, b2));
#endif
    _mm_store_sd(&C[(l_n*36)+8], c2_1);
    __m128d c2_2 = _mm_load_sd(&C[(l_n*36)+18]);
    __m128d a2_2 = _mm_load_sd(&A[11]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_2 = _mm_add_sd(c2_2, _mm_mul_sd(a2_2, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_2 = _mm_add_sd(c2_2, _mm_mul_sd(a2_2, b2));
#endif
    _mm_store_sd(&C[(l_n*36)+18], c2_2);
    __m128d c2_3 = _mm_load_sd(&C[(l_n*36)+33]);
    __m128d a2_3 = _mm_load_sd(&A[12]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_3 = _mm_add_sd(c2_3, _mm_mul_sd(a2_3, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_3 = _mm_add_sd(c2_3, _mm_mul_sd(a2_3, b2));
#endif
    _mm_store_sd(&C[(l_n*36)+33], c2_3);
#else
    C[(l_n*36)+2] += A[9] * B[(l_n*36)+2];
    C[(l_n*36)+8] += A[10] * B[(l_n*36)+2];
    C[(l_n*36)+18] += A[11] * B[(l_n*36)+2];
    C[(l_n*36)+33] += A[12] * B[(l_n*36)+2];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b3 = _mm256_broadcast_sd(&B[(l_n*36)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b3 = _mm_loaddup_pd(&B[(l_n*36)+3]);
#endif
    __m128d c3_0 = _mm_load_sd(&C[(l_n*36)+0]);
    __m128d a3_0 = _mm_load_sd(&A[13]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, b3));
#endif
    _mm_store_sd(&C[(l_n*36)+0], c3_0);
    __m128d c3_1 = _mm_load_sd(&C[(l_n*36)+3]);
    __m128d a3_1 = _mm_load_sd(&A[14]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, b3));
#endif
    _mm_store_sd(&C[(l_n*36)+3], c3_1);
    __m128d c3_2 = _mm_load_sd(&C[(l_n*36)+9]);
    __m128d a3_2 = _mm_load_sd(&A[15]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, b3));
#endif
    _mm_store_sd(&C[(l_n*36)+9], c3_2);
    __m128d c3_3 = _mm_load_sd(&C[(l_n*36)+19]);
    __m128d a3_3 = _mm_load_sd(&A[16]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_3 = _mm_add_sd(c3_3, _mm_mul_sd(a3_3, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_3 = _mm_add_sd(c3_3, _mm_mul_sd(a3_3, b3));
#endif
    _mm_store_sd(&C[(l_n*36)+19], c3_3);
    __m128d c3_4 = _mm_load_sd(&C[(l_n*36)+34]);
    __m128d a3_4 = _mm_load_sd(&A[17]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_4 = _mm_add_sd(c3_4, _mm_mul_sd(a3_4, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_4 = _mm_add_sd(c3_4, _mm_mul_sd(a3_4, b3));
#endif
    _mm_store_sd(&C[(l_n*36)+34], c3_4);
#else
    C[(l_n*36)+0] += A[13] * B[(l_n*36)+3];
    C[(l_n*36)+3] += A[14] * B[(l_n*36)+3];
    C[(l_n*36)+9] += A[15] * B[(l_n*36)+3];
    C[(l_n*36)+19] += A[16] * B[(l_n*36)+3];
    C[(l_n*36)+34] += A[17] * B[(l_n*36)+3];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b4 = _mm256_broadcast_sd(&B[(l_n*36)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b4 = _mm_loaddup_pd(&B[(l_n*36)+4]);
#endif
    __m128d c4_0 = _mm_load_sd(&C[(l_n*36)+4]);
    __m128d a4_0 = _mm_load_sd(&A[18]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_0 = _mm_add_sd(c4_0, _mm_mul_sd(a4_0, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_0 = _mm_add_sd(c4_0, _mm_mul_sd(a4_0, b4));
#endif
    _mm_store_sd(&C[(l_n*36)+4], c4_0);
    __m128d c4_1 = _mm_load_sd(&C[(l_n*36)+14]);
    __m128d a4_1 = _mm_load_sd(&A[19]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_1 = _mm_add_sd(c4_1, _mm_mul_sd(a4_1, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_1 = _mm_add_sd(c4_1, _mm_mul_sd(a4_1, b4));
#endif
    _mm_store_sd(&C[(l_n*36)+14], c4_1);
    __m128d c4_2 = _mm_load_sd(&C[(l_n*36)+29]);
    __m128d a4_2 = _mm_load_sd(&A[20]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, b4));
#endif
    _mm_store_sd(&C[(l_n*36)+29], c4_2);
#else
    C[(l_n*36)+4] += A[18] * B[(l_n*36)+4];
    C[(l_n*36)+14] += A[19] * B[(l_n*36)+4];
    C[(l_n*36)+29] += A[20] * B[(l_n*36)+4];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b5 = _mm256_broadcast_sd(&B[(l_n*36)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b5 = _mm_loaddup_pd(&B[(l_n*36)+5]);
#endif
    __m128d c5_0 = _mm_load_sd(&C[(l_n*36)+5]);
    __m128d a5_0 = _mm_load_sd(&A[21]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_0 = _mm_add_sd(c5_0, _mm_mul_sd(a5_0, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_0 = _mm_add_sd(c5_0, _mm_mul_sd(a5_0, b5));
#endif
    _mm_store_sd(&C[(l_n*36)+5], c5_0);
    __m128d c5_1 = _mm_load_sd(&C[(l_n*36)+15]);
    __m128d a5_1 = _mm_load_sd(&A[22]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_1 = _mm_add_sd(c5_1, _mm_mul_sd(a5_1, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_1 = _mm_add_sd(c5_1, _mm_mul_sd(a5_1, b5));
#endif
    _mm_store_sd(&C[(l_n*36)+15], c5_1);
    __m128d c5_2 = _mm_load_sd(&C[(l_n*36)+30]);
    __m128d a5_2 = _mm_load_sd(&A[23]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, b5));
#endif
    _mm_store_sd(&C[(l_n*36)+30], c5_2);
#else
    C[(l_n*36)+5] += A[21] * B[(l_n*36)+5];
    C[(l_n*36)+15] += A[22] * B[(l_n*36)+5];
    C[(l_n*36)+30] += A[23] * B[(l_n*36)+5];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b6 = _mm256_broadcast_sd(&B[(l_n*36)+6]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b6 = _mm_loaddup_pd(&B[(l_n*36)+6]);
#endif
    __m128d c6_0 = _mm_load_sd(&C[(l_n*36)+6]);
    __m128d a6_0 = _mm_load_sd(&A[24]);
#if defined(__SSE3__) && defined(__AVX__)
    c6_0 = _mm_add_sd(c6_0, _mm_mul_sd(a6_0, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c6_0 = _mm_add_sd(c6_0, _mm_mul_sd(a6_0, b6));
#endif
    _mm_store_sd(&C[(l_n*36)+6], c6_0);
    __m128d c6_1 = _mm_load_sd(&C[(l_n*36)+16]);
    __m128d a6_1 = _mm_load_sd(&A[25]);
#if defined(__SSE3__) && defined(__AVX__)
    c6_1 = _mm_add_sd(c6_1, _mm_mul_sd(a6_1, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c6_1 = _mm_add_sd(c6_1, _mm_mul_sd(a6_1, b6));
#endif
    _mm_store_sd(&C[(l_n*36)+16], c6_1);
    __m128d c6_2 = _mm_load_sd(&C[(l_n*36)+31]);
    __m128d a6_2 = _mm_load_sd(&A[26]);
#if defined(__SSE3__) && defined(__AVX__)
    c6_2 = _mm_add_sd(c6_2, _mm_mul_sd(a6_2, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c6_2 = _mm_add_sd(c6_2, _mm_mul_sd(a6_2, b6));
#endif
    _mm_store_sd(&C[(l_n*36)+31], c6_2);
#else
    C[(l_n*36)+6] += A[24] * B[(l_n*36)+6];
    C[(l_n*36)+16] += A[25] * B[(l_n*36)+6];
    C[(l_n*36)+31] += A[26] * B[(l_n*36)+6];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b7 = _mm256_broadcast_sd(&B[(l_n*36)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b7 = _mm_loaddup_pd(&B[(l_n*36)+7]);
#endif
    __m128d c7_0 = _mm_load_sd(&C[(l_n*36)+1]);
    __m128d a7_0 = _mm_load_sd(&A[27]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_0 = _mm_add_sd(c7_0, _mm_mul_sd(a7_0, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_0 = _mm_add_sd(c7_0, _mm_mul_sd(a7_0, b7));
#endif
    _mm_store_sd(&C[(l_n*36)+1], c7_0);
    __m128d c7_1 = _mm_load_sd(&C[(l_n*36)+7]);
    __m128d a7_1 = _mm_load_sd(&A[28]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_1 = _mm_add_sd(c7_1, _mm_mul_sd(a7_1, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_1 = _mm_add_sd(c7_1, _mm_mul_sd(a7_1, b7));
#endif
    _mm_store_sd(&C[(l_n*36)+7], c7_1);
    __m128d c7_2 = _mm_load_sd(&C[(l_n*36)+17]);
    __m128d a7_2 = _mm_load_sd(&A[29]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_2 = _mm_add_sd(c7_2, _mm_mul_sd(a7_2, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_2 = _mm_add_sd(c7_2, _mm_mul_sd(a7_2, b7));
#endif
    _mm_store_sd(&C[(l_n*36)+17], c7_2);
    __m128d c7_3 = _mm_load_sd(&C[(l_n*36)+32]);
    __m128d a7_3 = _mm_load_sd(&A[30]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_3 = _mm_add_sd(c7_3, _mm_mul_sd(a7_3, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_3 = _mm_add_sd(c7_3, _mm_mul_sd(a7_3, b7));
#endif
    _mm_store_sd(&C[(l_n*36)+32], c7_3);
#else
    C[(l_n*36)+1] += A[27] * B[(l_n*36)+7];
    C[(l_n*36)+7] += A[28] * B[(l_n*36)+7];
    C[(l_n*36)+17] += A[29] * B[(l_n*36)+7];
    C[(l_n*36)+32] += A[30] * B[(l_n*36)+7];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b8 = _mm256_broadcast_sd(&B[(l_n*36)+8]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b8 = _mm_loaddup_pd(&B[(l_n*36)+8]);
#endif
    __m128d c8_0 = _mm_load_sd(&C[(l_n*36)+2]);
    __m128d a8_0 = _mm_load_sd(&A[31]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_0 = _mm_add_sd(c8_0, _mm_mul_sd(a8_0, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_0 = _mm_add_sd(c8_0, _mm_mul_sd(a8_0, b8));
#endif
    _mm_store_sd(&C[(l_n*36)+2], c8_0);
    __m128d c8_1 = _mm_load_sd(&C[(l_n*36)+8]);
    __m128d a8_1 = _mm_load_sd(&A[32]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_1 = _mm_add_sd(c8_1, _mm_mul_sd(a8_1, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_1 = _mm_add_sd(c8_1, _mm_mul_sd(a8_1, b8));
#endif
    _mm_store_sd(&C[(l_n*36)+8], c8_1);
    __m128d c8_2 = _mm_load_sd(&C[(l_n*36)+18]);
    __m128d a8_2 = _mm_load_sd(&A[33]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_2 = _mm_add_sd(c8_2, _mm_mul_sd(a8_2, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_2 = _mm_add_sd(c8_2, _mm_mul_sd(a8_2, b8));
#endif
    _mm_store_sd(&C[(l_n*36)+18], c8_2);
    __m128d c8_3 = _mm_load_sd(&C[(l_n*36)+33]);
    __m128d a8_3 = _mm_load_sd(&A[34]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_3 = _mm_add_sd(c8_3, _mm_mul_sd(a8_3, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_3 = _mm_add_sd(c8_3, _mm_mul_sd(a8_3, b8));
#endif
    _mm_store_sd(&C[(l_n*36)+33], c8_3);
#else
    C[(l_n*36)+2] += A[31] * B[(l_n*36)+8];
    C[(l_n*36)+8] += A[32] * B[(l_n*36)+8];
    C[(l_n*36)+18] += A[33] * B[(l_n*36)+8];
    C[(l_n*36)+33] += A[34] * B[(l_n*36)+8];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b9 = _mm256_broadcast_sd(&B[(l_n*36)+9]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b9 = _mm_loaddup_pd(&B[(l_n*36)+9]);
#endif
    __m128d c9_0 = _mm_load_sd(&C[(l_n*36)+0]);
    __m128d a9_0 = _mm_load_sd(&A[35]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, b9));
#endif
    _mm_store_sd(&C[(l_n*36)+0], c9_0);
    __m128d c9_1 = _mm_load_sd(&C[(l_n*36)+3]);
    __m128d a9_1 = _mm_load_sd(&A[36]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, b9));
#endif
    _mm_store_sd(&C[(l_n*36)+3], c9_1);
    __m128d c9_2 = _mm_load_sd(&C[(l_n*36)+9]);
    __m128d a9_2 = _mm_load_sd(&A[37]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, b9));
#endif
    _mm_store_sd(&C[(l_n*36)+9], c9_2);
    __m128d c9_3 = _mm_load_sd(&C[(l_n*36)+19]);
    __m128d a9_3 = _mm_load_sd(&A[38]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_3 = _mm_add_sd(c9_3, _mm_mul_sd(a9_3, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_3 = _mm_add_sd(c9_3, _mm_mul_sd(a9_3, b9));
#endif
    _mm_store_sd(&C[(l_n*36)+19], c9_3);
    __m128d c9_4 = _mm_load_sd(&C[(l_n*36)+34]);
    __m128d a9_4 = _mm_load_sd(&A[39]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_4 = _mm_add_sd(c9_4, _mm_mul_sd(a9_4, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_4 = _mm_add_sd(c9_4, _mm_mul_sd(a9_4, b9));
#endif
    _mm_store_sd(&C[(l_n*36)+34], c9_4);
#else
    C[(l_n*36)+0] += A[35] * B[(l_n*36)+9];
    C[(l_n*36)+3] += A[36] * B[(l_n*36)+9];
    C[(l_n*36)+9] += A[37] * B[(l_n*36)+9];
    C[(l_n*36)+19] += A[38] * B[(l_n*36)+9];
    C[(l_n*36)+34] += A[39] * B[(l_n*36)+9];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b10 = _mm256_broadcast_sd(&B[(l_n*36)+10]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b10 = _mm_loaddup_pd(&B[(l_n*36)+10]);
#endif
    __m128d c10_0 = _mm_load_sd(&C[(l_n*36)+10]);
    __m128d a10_0 = _mm_load_sd(&A[40]);
#if defined(__SSE3__) && defined(__AVX__)
    c10_0 = _mm_add_sd(c10_0, _mm_mul_sd(a10_0, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c10_0 = _mm_add_sd(c10_0, _mm_mul_sd(a10_0, b10));
#endif
    _mm_store_sd(&C[(l_n*36)+10], c10_0);
    __m128d c10_1 = _mm_load_sd(&C[(l_n*36)+25]);
    __m128d a10_1 = _mm_load_sd(&A[41]);
#if defined(__SSE3__) && defined(__AVX__)
    c10_1 = _mm_add_sd(c10_1, _mm_mul_sd(a10_1, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c10_1 = _mm_add_sd(c10_1, _mm_mul_sd(a10_1, b10));
#endif
    _mm_store_sd(&C[(l_n*36)+25], c10_1);
#else
    C[(l_n*36)+10] += A[40] * B[(l_n*36)+10];
    C[(l_n*36)+25] += A[41] * B[(l_n*36)+10];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b11 = _mm256_broadcast_sd(&B[(l_n*36)+11]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b11 = _mm_loaddup_pd(&B[(l_n*36)+11]);
#endif
    __m128d c11_0 = _mm_load_sd(&C[(l_n*36)+11]);
    __m128d a11_0 = _mm_load_sd(&A[42]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_0 = _mm_add_sd(c11_0, _mm_mul_sd(a11_0, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_0 = _mm_add_sd(c11_0, _mm_mul_sd(a11_0, b11));
#endif
    _mm_store_sd(&C[(l_n*36)+11], c11_0);
    __m128d c11_1 = _mm_load_sd(&C[(l_n*36)+26]);
    __m128d a11_1 = _mm_load_sd(&A[43]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_1 = _mm_add_sd(c11_1, _mm_mul_sd(a11_1, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_1 = _mm_add_sd(c11_1, _mm_mul_sd(a11_1, b11));
#endif
    _mm_store_sd(&C[(l_n*36)+26], c11_1);
#else
    C[(l_n*36)+11] += A[42] * B[(l_n*36)+11];
    C[(l_n*36)+26] += A[43] * B[(l_n*36)+11];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b12 = _mm256_broadcast_sd(&B[(l_n*36)+12]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b12 = _mm_loaddup_pd(&B[(l_n*36)+12]);
#endif
    __m128d c12_0 = _mm_load_sd(&C[(l_n*36)+12]);
    __m128d a12_0 = _mm_load_sd(&A[44]);
#if defined(__SSE3__) && defined(__AVX__)
    c12_0 = _mm_add_sd(c12_0, _mm_mul_sd(a12_0, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c12_0 = _mm_add_sd(c12_0, _mm_mul_sd(a12_0, b12));
#endif
    _mm_store_sd(&C[(l_n*36)+12], c12_0);
    __m128d c12_1 = _mm_load_sd(&C[(l_n*36)+27]);
    __m128d a12_1 = _mm_load_sd(&A[45]);
#if defined(__SSE3__) && defined(__AVX__)
    c12_1 = _mm_add_sd(c12_1, _mm_mul_sd(a12_1, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c12_1 = _mm_add_sd(c12_1, _mm_mul_sd(a12_1, b12));
#endif
    _mm_store_sd(&C[(l_n*36)+27], c12_1);
#else
    C[(l_n*36)+12] += A[44] * B[(l_n*36)+12];
    C[(l_n*36)+27] += A[45] * B[(l_n*36)+12];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b13 = _mm256_broadcast_sd(&B[(l_n*36)+13]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b13 = _mm_loaddup_pd(&B[(l_n*36)+13]);
#endif
    __m128d c13_0 = _mm_load_sd(&C[(l_n*36)+13]);
    __m128d a13_0 = _mm_load_sd(&A[46]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_0 = _mm_add_sd(c13_0, _mm_mul_sd(a13_0, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_0 = _mm_add_sd(c13_0, _mm_mul_sd(a13_0, b13));
#endif
    _mm_store_sd(&C[(l_n*36)+13], c13_0);
    __m128d c13_1 = _mm_load_sd(&C[(l_n*36)+28]);
    __m128d a13_1 = _mm_load_sd(&A[47]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_1 = _mm_add_sd(c13_1, _mm_mul_sd(a13_1, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_1 = _mm_add_sd(c13_1, _mm_mul_sd(a13_1, b13));
#endif
    _mm_store_sd(&C[(l_n*36)+28], c13_1);
#else
    C[(l_n*36)+13] += A[46] * B[(l_n*36)+13];
    C[(l_n*36)+28] += A[47] * B[(l_n*36)+13];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b14 = _mm256_broadcast_sd(&B[(l_n*36)+14]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b14 = _mm_loaddup_pd(&B[(l_n*36)+14]);
#endif
    __m128d c14_0 = _mm_load_sd(&C[(l_n*36)+4]);
    __m128d a14_0 = _mm_load_sd(&A[48]);
#if defined(__SSE3__) && defined(__AVX__)
    c14_0 = _mm_add_sd(c14_0, _mm_mul_sd(a14_0, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c14_0 = _mm_add_sd(c14_0, _mm_mul_sd(a14_0, b14));
#endif
    _mm_store_sd(&C[(l_n*36)+4], c14_0);
    __m128d c14_1 = _mm_load_sd(&C[(l_n*36)+14]);
    __m128d a14_1 = _mm_load_sd(&A[49]);
#if defined(__SSE3__) && defined(__AVX__)
    c14_1 = _mm_add_sd(c14_1, _mm_mul_sd(a14_1, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c14_1 = _mm_add_sd(c14_1, _mm_mul_sd(a14_1, b14));
#endif
    _mm_store_sd(&C[(l_n*36)+14], c14_1);
    __m128d c14_2 = _mm_load_sd(&C[(l_n*36)+29]);
    __m128d a14_2 = _mm_load_sd(&A[50]);
#if defined(__SSE3__) && defined(__AVX__)
    c14_2 = _mm_add_sd(c14_2, _mm_mul_sd(a14_2, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c14_2 = _mm_add_sd(c14_2, _mm_mul_sd(a14_2, b14));
#endif
    _mm_store_sd(&C[(l_n*36)+29], c14_2);
#else
    C[(l_n*36)+4] += A[48] * B[(l_n*36)+14];
    C[(l_n*36)+14] += A[49] * B[(l_n*36)+14];
    C[(l_n*36)+29] += A[50] * B[(l_n*36)+14];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b15 = _mm256_broadcast_sd(&B[(l_n*36)+15]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b15 = _mm_loaddup_pd(&B[(l_n*36)+15]);
#endif
    __m128d c15_0 = _mm_load_sd(&C[(l_n*36)+5]);
    __m128d a15_0 = _mm_load_sd(&A[51]);
#if defined(__SSE3__) && defined(__AVX__)
    c15_0 = _mm_add_sd(c15_0, _mm_mul_sd(a15_0, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c15_0 = _mm_add_sd(c15_0, _mm_mul_sd(a15_0, b15));
#endif
    _mm_store_sd(&C[(l_n*36)+5], c15_0);
    __m128d c15_1 = _mm_load_sd(&C[(l_n*36)+15]);
    __m128d a15_1 = _mm_load_sd(&A[52]);
#if defined(__SSE3__) && defined(__AVX__)
    c15_1 = _mm_add_sd(c15_1, _mm_mul_sd(a15_1, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c15_1 = _mm_add_sd(c15_1, _mm_mul_sd(a15_1, b15));
#endif
    _mm_store_sd(&C[(l_n*36)+15], c15_1);
    __m128d c15_2 = _mm_load_sd(&C[(l_n*36)+30]);
    __m128d a15_2 = _mm_load_sd(&A[53]);
#if defined(__SSE3__) && defined(__AVX__)
    c15_2 = _mm_add_sd(c15_2, _mm_mul_sd(a15_2, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c15_2 = _mm_add_sd(c15_2, _mm_mul_sd(a15_2, b15));
#endif
    _mm_store_sd(&C[(l_n*36)+30], c15_2);
#else
    C[(l_n*36)+5] += A[51] * B[(l_n*36)+15];
    C[(l_n*36)+15] += A[52] * B[(l_n*36)+15];
    C[(l_n*36)+30] += A[53] * B[(l_n*36)+15];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b16 = _mm256_broadcast_sd(&B[(l_n*36)+16]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b16 = _mm_loaddup_pd(&B[(l_n*36)+16]);
#endif
    __m128d c16_0 = _mm_load_sd(&C[(l_n*36)+6]);
    __m128d a16_0 = _mm_load_sd(&A[54]);
#if defined(__SSE3__) && defined(__AVX__)
    c16_0 = _mm_add_sd(c16_0, _mm_mul_sd(a16_0, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c16_0 = _mm_add_sd(c16_0, _mm_mul_sd(a16_0, b16));
#endif
    _mm_store_sd(&C[(l_n*36)+6], c16_0);
    __m128d c16_1 = _mm_load_sd(&C[(l_n*36)+16]);
    __m128d a16_1 = _mm_load_sd(&A[55]);
#if defined(__SSE3__) && defined(__AVX__)
    c16_1 = _mm_add_sd(c16_1, _mm_mul_sd(a16_1, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c16_1 = _mm_add_sd(c16_1, _mm_mul_sd(a16_1, b16));
#endif
    _mm_store_sd(&C[(l_n*36)+16], c16_1);
    __m128d c16_2 = _mm_load_sd(&C[(l_n*36)+31]);
    __m128d a16_2 = _mm_load_sd(&A[56]);
#if defined(__SSE3__) && defined(__AVX__)
    c16_2 = _mm_add_sd(c16_2, _mm_mul_sd(a16_2, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c16_2 = _mm_add_sd(c16_2, _mm_mul_sd(a16_2, b16));
#endif
    _mm_store_sd(&C[(l_n*36)+31], c16_2);
#else
    C[(l_n*36)+6] += A[54] * B[(l_n*36)+16];
    C[(l_n*36)+16] += A[55] * B[(l_n*36)+16];
    C[(l_n*36)+31] += A[56] * B[(l_n*36)+16];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b17 = _mm256_broadcast_sd(&B[(l_n*36)+17]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b17 = _mm_loaddup_pd(&B[(l_n*36)+17]);
#endif
    __m128d c17_0 = _mm_load_sd(&C[(l_n*36)+1]);
    __m128d a17_0 = _mm_load_sd(&A[57]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_0 = _mm_add_sd(c17_0, _mm_mul_sd(a17_0, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_0 = _mm_add_sd(c17_0, _mm_mul_sd(a17_0, b17));
#endif
    _mm_store_sd(&C[(l_n*36)+1], c17_0);
    __m128d c17_1 = _mm_load_sd(&C[(l_n*36)+7]);
    __m128d a17_1 = _mm_load_sd(&A[58]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_1 = _mm_add_sd(c17_1, _mm_mul_sd(a17_1, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_1 = _mm_add_sd(c17_1, _mm_mul_sd(a17_1, b17));
#endif
    _mm_store_sd(&C[(l_n*36)+7], c17_1);
    __m128d c17_2 = _mm_load_sd(&C[(l_n*36)+17]);
    __m128d a17_2 = _mm_load_sd(&A[59]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_2 = _mm_add_sd(c17_2, _mm_mul_sd(a17_2, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_2 = _mm_add_sd(c17_2, _mm_mul_sd(a17_2, b17));
#endif
    _mm_store_sd(&C[(l_n*36)+17], c17_2);
    __m128d c17_3 = _mm_load_sd(&C[(l_n*36)+32]);
    __m128d a17_3 = _mm_load_sd(&A[60]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_3 = _mm_add_sd(c17_3, _mm_mul_sd(a17_3, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_3 = _mm_add_sd(c17_3, _mm_mul_sd(a17_3, b17));
#endif
    _mm_store_sd(&C[(l_n*36)+32], c17_3);
#else
    C[(l_n*36)+1] += A[57] * B[(l_n*36)+17];
    C[(l_n*36)+7] += A[58] * B[(l_n*36)+17];
    C[(l_n*36)+17] += A[59] * B[(l_n*36)+17];
    C[(l_n*36)+32] += A[60] * B[(l_n*36)+17];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b18 = _mm256_broadcast_sd(&B[(l_n*36)+18]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b18 = _mm_loaddup_pd(&B[(l_n*36)+18]);
#endif
    __m128d c18_0 = _mm_load_sd(&C[(l_n*36)+2]);
    __m128d a18_0 = _mm_load_sd(&A[61]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_0 = _mm_add_sd(c18_0, _mm_mul_sd(a18_0, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_0 = _mm_add_sd(c18_0, _mm_mul_sd(a18_0, b18));
#endif
    _mm_store_sd(&C[(l_n*36)+2], c18_0);
    __m128d c18_1 = _mm_load_sd(&C[(l_n*36)+8]);
    __m128d a18_1 = _mm_load_sd(&A[62]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_1 = _mm_add_sd(c18_1, _mm_mul_sd(a18_1, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_1 = _mm_add_sd(c18_1, _mm_mul_sd(a18_1, b18));
#endif
    _mm_store_sd(&C[(l_n*36)+8], c18_1);
    __m128d c18_2 = _mm_load_sd(&C[(l_n*36)+18]);
    __m128d a18_2 = _mm_load_sd(&A[63]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_2 = _mm_add_sd(c18_2, _mm_mul_sd(a18_2, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_2 = _mm_add_sd(c18_2, _mm_mul_sd(a18_2, b18));
#endif
    _mm_store_sd(&C[(l_n*36)+18], c18_2);
    __m128d c18_3 = _mm_load_sd(&C[(l_n*36)+33]);
    __m128d a18_3 = _mm_load_sd(&A[64]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_3 = _mm_add_sd(c18_3, _mm_mul_sd(a18_3, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_3 = _mm_add_sd(c18_3, _mm_mul_sd(a18_3, b18));
#endif
    _mm_store_sd(&C[(l_n*36)+33], c18_3);
#else
    C[(l_n*36)+2] += A[61] * B[(l_n*36)+18];
    C[(l_n*36)+8] += A[62] * B[(l_n*36)+18];
    C[(l_n*36)+18] += A[63] * B[(l_n*36)+18];
    C[(l_n*36)+33] += A[64] * B[(l_n*36)+18];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b19 = _mm256_broadcast_sd(&B[(l_n*36)+19]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b19 = _mm_loaddup_pd(&B[(l_n*36)+19]);
#endif
    __m128d c19_0 = _mm_load_sd(&C[(l_n*36)+0]);
    __m128d a19_0 = _mm_load_sd(&A[65]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_0 = _mm_add_sd(c19_0, _mm_mul_sd(a19_0, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_0 = _mm_add_sd(c19_0, _mm_mul_sd(a19_0, b19));
#endif
    _mm_store_sd(&C[(l_n*36)+0], c19_0);
    __m128d c19_1 = _mm_load_sd(&C[(l_n*36)+3]);
    __m128d a19_1 = _mm_load_sd(&A[66]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_1 = _mm_add_sd(c19_1, _mm_mul_sd(a19_1, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_1 = _mm_add_sd(c19_1, _mm_mul_sd(a19_1, b19));
#endif
    _mm_store_sd(&C[(l_n*36)+3], c19_1);
    __m128d c19_2 = _mm_load_sd(&C[(l_n*36)+9]);
    __m128d a19_2 = _mm_load_sd(&A[67]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_2 = _mm_add_sd(c19_2, _mm_mul_sd(a19_2, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_2 = _mm_add_sd(c19_2, _mm_mul_sd(a19_2, b19));
#endif
    _mm_store_sd(&C[(l_n*36)+9], c19_2);
    __m128d c19_3 = _mm_load_sd(&C[(l_n*36)+19]);
    __m128d a19_3 = _mm_load_sd(&A[68]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_3 = _mm_add_sd(c19_3, _mm_mul_sd(a19_3, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_3 = _mm_add_sd(c19_3, _mm_mul_sd(a19_3, b19));
#endif
    _mm_store_sd(&C[(l_n*36)+19], c19_3);
    __m128d c19_4 = _mm_load_sd(&C[(l_n*36)+34]);
    __m128d a19_4 = _mm_load_sd(&A[69]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_4 = _mm_add_sd(c19_4, _mm_mul_sd(a19_4, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_4 = _mm_add_sd(c19_4, _mm_mul_sd(a19_4, b19));
#endif
    _mm_store_sd(&C[(l_n*36)+34], c19_4);
#else
    C[(l_n*36)+0] += A[65] * B[(l_n*36)+19];
    C[(l_n*36)+3] += A[66] * B[(l_n*36)+19];
    C[(l_n*36)+9] += A[67] * B[(l_n*36)+19];
    C[(l_n*36)+19] += A[68] * B[(l_n*36)+19];
    C[(l_n*36)+34] += A[69] * B[(l_n*36)+19];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b20 = _mm256_broadcast_sd(&B[(l_n*36)+20]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b20 = _mm_loaddup_pd(&B[(l_n*36)+20]);
#endif
    __m128d c20_0 = _mm_load_sd(&C[(l_n*36)+20]);
    __m128d a20_0 = _mm_load_sd(&A[70]);
#if defined(__SSE3__) && defined(__AVX__)
    c20_0 = _mm_add_sd(c20_0, _mm_mul_sd(a20_0, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c20_0 = _mm_add_sd(c20_0, _mm_mul_sd(a20_0, b20));
#endif
    _mm_store_sd(&C[(l_n*36)+20], c20_0);
#else
    C[(l_n*36)+20] += A[70] * B[(l_n*36)+20];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b21 = _mm256_broadcast_sd(&B[(l_n*36)+21]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b21 = _mm_loaddup_pd(&B[(l_n*36)+21]);
#endif
    __m128d c21_0 = _mm_load_sd(&C[(l_n*36)+21]);
    __m128d a21_0 = _mm_load_sd(&A[71]);
#if defined(__SSE3__) && defined(__AVX__)
    c21_0 = _mm_add_sd(c21_0, _mm_mul_sd(a21_0, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c21_0 = _mm_add_sd(c21_0, _mm_mul_sd(a21_0, b21));
#endif
    _mm_store_sd(&C[(l_n*36)+21], c21_0);
#else
    C[(l_n*36)+21] += A[71] * B[(l_n*36)+21];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b22 = _mm256_broadcast_sd(&B[(l_n*36)+22]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b22 = _mm_loaddup_pd(&B[(l_n*36)+22]);
#endif
    __m128d c22_0 = _mm_load_sd(&C[(l_n*36)+22]);
    __m128d a22_0 = _mm_load_sd(&A[72]);
#if defined(__SSE3__) && defined(__AVX__)
    c22_0 = _mm_add_sd(c22_0, _mm_mul_sd(a22_0, _mm256_castpd256_pd128(b22)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c22_0 = _mm_add_sd(c22_0, _mm_mul_sd(a22_0, b22));
#endif
    _mm_store_sd(&C[(l_n*36)+22], c22_0);
#else
    C[(l_n*36)+22] += A[72] * B[(l_n*36)+22];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b23 = _mm256_broadcast_sd(&B[(l_n*36)+23]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b23 = _mm_loaddup_pd(&B[(l_n*36)+23]);
#endif
    __m128d c23_0 = _mm_load_sd(&C[(l_n*36)+23]);
    __m128d a23_0 = _mm_load_sd(&A[73]);
#if defined(__SSE3__) && defined(__AVX__)
    c23_0 = _mm_add_sd(c23_0, _mm_mul_sd(a23_0, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c23_0 = _mm_add_sd(c23_0, _mm_mul_sd(a23_0, b23));
#endif
    _mm_store_sd(&C[(l_n*36)+23], c23_0);
#else
    C[(l_n*36)+23] += A[73] * B[(l_n*36)+23];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b24 = _mm256_broadcast_sd(&B[(l_n*36)+24]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b24 = _mm_loaddup_pd(&B[(l_n*36)+24]);
#endif
    __m128d c24_0 = _mm_load_sd(&C[(l_n*36)+24]);
    __m128d a24_0 = _mm_load_sd(&A[74]);
#if defined(__SSE3__) && defined(__AVX__)
    c24_0 = _mm_add_sd(c24_0, _mm_mul_sd(a24_0, _mm256_castpd256_pd128(b24)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c24_0 = _mm_add_sd(c24_0, _mm_mul_sd(a24_0, b24));
#endif
    _mm_store_sd(&C[(l_n*36)+24], c24_0);
#else
    C[(l_n*36)+24] += A[74] * B[(l_n*36)+24];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b25 = _mm256_broadcast_sd(&B[(l_n*36)+25]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b25 = _mm_loaddup_pd(&B[(l_n*36)+25]);
#endif
    __m128d c25_0 = _mm_load_sd(&C[(l_n*36)+10]);
    __m128d a25_0 = _mm_load_sd(&A[75]);
#if defined(__SSE3__) && defined(__AVX__)
    c25_0 = _mm_add_sd(c25_0, _mm_mul_sd(a25_0, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c25_0 = _mm_add_sd(c25_0, _mm_mul_sd(a25_0, b25));
#endif
    _mm_store_sd(&C[(l_n*36)+10], c25_0);
    __m128d c25_1 = _mm_load_sd(&C[(l_n*36)+25]);
    __m128d a25_1 = _mm_load_sd(&A[76]);
#if defined(__SSE3__) && defined(__AVX__)
    c25_1 = _mm_add_sd(c25_1, _mm_mul_sd(a25_1, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c25_1 = _mm_add_sd(c25_1, _mm_mul_sd(a25_1, b25));
#endif
    _mm_store_sd(&C[(l_n*36)+25], c25_1);
#else
    C[(l_n*36)+10] += A[75] * B[(l_n*36)+25];
    C[(l_n*36)+25] += A[76] * B[(l_n*36)+25];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b26 = _mm256_broadcast_sd(&B[(l_n*36)+26]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b26 = _mm_loaddup_pd(&B[(l_n*36)+26]);
#endif
    __m128d c26_0 = _mm_load_sd(&C[(l_n*36)+11]);
    __m128d a26_0 = _mm_load_sd(&A[77]);
#if defined(__SSE3__) && defined(__AVX__)
    c26_0 = _mm_add_sd(c26_0, _mm_mul_sd(a26_0, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c26_0 = _mm_add_sd(c26_0, _mm_mul_sd(a26_0, b26));
#endif
    _mm_store_sd(&C[(l_n*36)+11], c26_0);
    __m128d c26_1 = _mm_load_sd(&C[(l_n*36)+26]);
    __m128d a26_1 = _mm_load_sd(&A[78]);
#if defined(__SSE3__) && defined(__AVX__)
    c26_1 = _mm_add_sd(c26_1, _mm_mul_sd(a26_1, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c26_1 = _mm_add_sd(c26_1, _mm_mul_sd(a26_1, b26));
#endif
    _mm_store_sd(&C[(l_n*36)+26], c26_1);
#else
    C[(l_n*36)+11] += A[77] * B[(l_n*36)+26];
    C[(l_n*36)+26] += A[78] * B[(l_n*36)+26];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b27 = _mm256_broadcast_sd(&B[(l_n*36)+27]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b27 = _mm_loaddup_pd(&B[(l_n*36)+27]);
#endif
    __m128d c27_0 = _mm_load_sd(&C[(l_n*36)+12]);
    __m128d a27_0 = _mm_load_sd(&A[79]);
#if defined(__SSE3__) && defined(__AVX__)
    c27_0 = _mm_add_sd(c27_0, _mm_mul_sd(a27_0, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c27_0 = _mm_add_sd(c27_0, _mm_mul_sd(a27_0, b27));
#endif
    _mm_store_sd(&C[(l_n*36)+12], c27_0);
    __m128d c27_1 = _mm_load_sd(&C[(l_n*36)+27]);
    __m128d a27_1 = _mm_load_sd(&A[80]);
#if defined(__SSE3__) && defined(__AVX__)
    c27_1 = _mm_add_sd(c27_1, _mm_mul_sd(a27_1, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c27_1 = _mm_add_sd(c27_1, _mm_mul_sd(a27_1, b27));
#endif
    _mm_store_sd(&C[(l_n*36)+27], c27_1);
#else
    C[(l_n*36)+12] += A[79] * B[(l_n*36)+27];
    C[(l_n*36)+27] += A[80] * B[(l_n*36)+27];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b28 = _mm256_broadcast_sd(&B[(l_n*36)+28]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b28 = _mm_loaddup_pd(&B[(l_n*36)+28]);
#endif
    __m128d c28_0 = _mm_load_sd(&C[(l_n*36)+13]);
    __m128d a28_0 = _mm_load_sd(&A[81]);
#if defined(__SSE3__) && defined(__AVX__)
    c28_0 = _mm_add_sd(c28_0, _mm_mul_sd(a28_0, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c28_0 = _mm_add_sd(c28_0, _mm_mul_sd(a28_0, b28));
#endif
    _mm_store_sd(&C[(l_n*36)+13], c28_0);
    __m128d c28_1 = _mm_load_sd(&C[(l_n*36)+28]);
    __m128d a28_1 = _mm_load_sd(&A[82]);
#if defined(__SSE3__) && defined(__AVX__)
    c28_1 = _mm_add_sd(c28_1, _mm_mul_sd(a28_1, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c28_1 = _mm_add_sd(c28_1, _mm_mul_sd(a28_1, b28));
#endif
    _mm_store_sd(&C[(l_n*36)+28], c28_1);
#else
    C[(l_n*36)+13] += A[81] * B[(l_n*36)+28];
    C[(l_n*36)+28] += A[82] * B[(l_n*36)+28];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b29 = _mm256_broadcast_sd(&B[(l_n*36)+29]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b29 = _mm_loaddup_pd(&B[(l_n*36)+29]);
#endif
    __m128d c29_0 = _mm_load_sd(&C[(l_n*36)+4]);
    __m128d a29_0 = _mm_load_sd(&A[83]);
#if defined(__SSE3__) && defined(__AVX__)
    c29_0 = _mm_add_sd(c29_0, _mm_mul_sd(a29_0, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c29_0 = _mm_add_sd(c29_0, _mm_mul_sd(a29_0, b29));
#endif
    _mm_store_sd(&C[(l_n*36)+4], c29_0);
    __m128d c29_1 = _mm_load_sd(&C[(l_n*36)+14]);
    __m128d a29_1 = _mm_load_sd(&A[84]);
#if defined(__SSE3__) && defined(__AVX__)
    c29_1 = _mm_add_sd(c29_1, _mm_mul_sd(a29_1, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c29_1 = _mm_add_sd(c29_1, _mm_mul_sd(a29_1, b29));
#endif
    _mm_store_sd(&C[(l_n*36)+14], c29_1);
    __m128d c29_2 = _mm_load_sd(&C[(l_n*36)+29]);
    __m128d a29_2 = _mm_load_sd(&A[85]);
#if defined(__SSE3__) && defined(__AVX__)
    c29_2 = _mm_add_sd(c29_2, _mm_mul_sd(a29_2, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c29_2 = _mm_add_sd(c29_2, _mm_mul_sd(a29_2, b29));
#endif
    _mm_store_sd(&C[(l_n*36)+29], c29_2);
#else
    C[(l_n*36)+4] += A[83] * B[(l_n*36)+29];
    C[(l_n*36)+14] += A[84] * B[(l_n*36)+29];
    C[(l_n*36)+29] += A[85] * B[(l_n*36)+29];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b30 = _mm256_broadcast_sd(&B[(l_n*36)+30]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b30 = _mm_loaddup_pd(&B[(l_n*36)+30]);
#endif
    __m128d c30_0 = _mm_load_sd(&C[(l_n*36)+5]);
    __m128d a30_0 = _mm_load_sd(&A[86]);
#if defined(__SSE3__) && defined(__AVX__)
    c30_0 = _mm_add_sd(c30_0, _mm_mul_sd(a30_0, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c30_0 = _mm_add_sd(c30_0, _mm_mul_sd(a30_0, b30));
#endif
    _mm_store_sd(&C[(l_n*36)+5], c30_0);
    __m128d c30_1 = _mm_load_sd(&C[(l_n*36)+15]);
    __m128d a30_1 = _mm_load_sd(&A[87]);
#if defined(__SSE3__) && defined(__AVX__)
    c30_1 = _mm_add_sd(c30_1, _mm_mul_sd(a30_1, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c30_1 = _mm_add_sd(c30_1, _mm_mul_sd(a30_1, b30));
#endif
    _mm_store_sd(&C[(l_n*36)+15], c30_1);
    __m128d c30_2 = _mm_load_sd(&C[(l_n*36)+30]);
    __m128d a30_2 = _mm_load_sd(&A[88]);
#if defined(__SSE3__) && defined(__AVX__)
    c30_2 = _mm_add_sd(c30_2, _mm_mul_sd(a30_2, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c30_2 = _mm_add_sd(c30_2, _mm_mul_sd(a30_2, b30));
#endif
    _mm_store_sd(&C[(l_n*36)+30], c30_2);
#else
    C[(l_n*36)+5] += A[86] * B[(l_n*36)+30];
    C[(l_n*36)+15] += A[87] * B[(l_n*36)+30];
    C[(l_n*36)+30] += A[88] * B[(l_n*36)+30];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b31 = _mm256_broadcast_sd(&B[(l_n*36)+31]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b31 = _mm_loaddup_pd(&B[(l_n*36)+31]);
#endif
    __m128d c31_0 = _mm_load_sd(&C[(l_n*36)+6]);
    __m128d a31_0 = _mm_load_sd(&A[89]);
#if defined(__SSE3__) && defined(__AVX__)
    c31_0 = _mm_add_sd(c31_0, _mm_mul_sd(a31_0, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c31_0 = _mm_add_sd(c31_0, _mm_mul_sd(a31_0, b31));
#endif
    _mm_store_sd(&C[(l_n*36)+6], c31_0);
    __m128d c31_1 = _mm_load_sd(&C[(l_n*36)+16]);
    __m128d a31_1 = _mm_load_sd(&A[90]);
#if defined(__SSE3__) && defined(__AVX__)
    c31_1 = _mm_add_sd(c31_1, _mm_mul_sd(a31_1, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c31_1 = _mm_add_sd(c31_1, _mm_mul_sd(a31_1, b31));
#endif
    _mm_store_sd(&C[(l_n*36)+16], c31_1);
    __m128d c31_2 = _mm_load_sd(&C[(l_n*36)+31]);
    __m128d a31_2 = _mm_load_sd(&A[91]);
#if defined(__SSE3__) && defined(__AVX__)
    c31_2 = _mm_add_sd(c31_2, _mm_mul_sd(a31_2, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c31_2 = _mm_add_sd(c31_2, _mm_mul_sd(a31_2, b31));
#endif
    _mm_store_sd(&C[(l_n*36)+31], c31_2);
#else
    C[(l_n*36)+6] += A[89] * B[(l_n*36)+31];
    C[(l_n*36)+16] += A[90] * B[(l_n*36)+31];
    C[(l_n*36)+31] += A[91] * B[(l_n*36)+31];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b32 = _mm256_broadcast_sd(&B[(l_n*36)+32]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b32 = _mm_loaddup_pd(&B[(l_n*36)+32]);
#endif
    __m128d c32_0 = _mm_load_sd(&C[(l_n*36)+1]);
    __m128d a32_0 = _mm_load_sd(&A[92]);
#if defined(__SSE3__) && defined(__AVX__)
    c32_0 = _mm_add_sd(c32_0, _mm_mul_sd(a32_0, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c32_0 = _mm_add_sd(c32_0, _mm_mul_sd(a32_0, b32));
#endif
    _mm_store_sd(&C[(l_n*36)+1], c32_0);
    __m128d c32_1 = _mm_load_sd(&C[(l_n*36)+7]);
    __m128d a32_1 = _mm_load_sd(&A[93]);
#if defined(__SSE3__) && defined(__AVX__)
    c32_1 = _mm_add_sd(c32_1, _mm_mul_sd(a32_1, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c32_1 = _mm_add_sd(c32_1, _mm_mul_sd(a32_1, b32));
#endif
    _mm_store_sd(&C[(l_n*36)+7], c32_1);
    __m128d c32_2 = _mm_load_sd(&C[(l_n*36)+17]);
    __m128d a32_2 = _mm_load_sd(&A[94]);
#if defined(__SSE3__) && defined(__AVX__)
    c32_2 = _mm_add_sd(c32_2, _mm_mul_sd(a32_2, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c32_2 = _mm_add_sd(c32_2, _mm_mul_sd(a32_2, b32));
#endif
    _mm_store_sd(&C[(l_n*36)+17], c32_2);
    __m128d c32_3 = _mm_load_sd(&C[(l_n*36)+32]);
    __m128d a32_3 = _mm_load_sd(&A[95]);
#if defined(__SSE3__) && defined(__AVX__)
    c32_3 = _mm_add_sd(c32_3, _mm_mul_sd(a32_3, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c32_3 = _mm_add_sd(c32_3, _mm_mul_sd(a32_3, b32));
#endif
    _mm_store_sd(&C[(l_n*36)+32], c32_3);
#else
    C[(l_n*36)+1] += A[92] * B[(l_n*36)+32];
    C[(l_n*36)+7] += A[93] * B[(l_n*36)+32];
    C[(l_n*36)+17] += A[94] * B[(l_n*36)+32];
    C[(l_n*36)+32] += A[95] * B[(l_n*36)+32];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b33 = _mm256_broadcast_sd(&B[(l_n*36)+33]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b33 = _mm_loaddup_pd(&B[(l_n*36)+33]);
#endif
    __m128d c33_0 = _mm_load_sd(&C[(l_n*36)+2]);
    __m128d a33_0 = _mm_load_sd(&A[96]);
#if defined(__SSE3__) && defined(__AVX__)
    c33_0 = _mm_add_sd(c33_0, _mm_mul_sd(a33_0, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c33_0 = _mm_add_sd(c33_0, _mm_mul_sd(a33_0, b33));
#endif
    _mm_store_sd(&C[(l_n*36)+2], c33_0);
    __m128d c33_1 = _mm_load_sd(&C[(l_n*36)+8]);
    __m128d a33_1 = _mm_load_sd(&A[97]);
#if defined(__SSE3__) && defined(__AVX__)
    c33_1 = _mm_add_sd(c33_1, _mm_mul_sd(a33_1, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c33_1 = _mm_add_sd(c33_1, _mm_mul_sd(a33_1, b33));
#endif
    _mm_store_sd(&C[(l_n*36)+8], c33_1);
    __m128d c33_2 = _mm_load_sd(&C[(l_n*36)+18]);
    __m128d a33_2 = _mm_load_sd(&A[98]);
#if defined(__SSE3__) && defined(__AVX__)
    c33_2 = _mm_add_sd(c33_2, _mm_mul_sd(a33_2, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c33_2 = _mm_add_sd(c33_2, _mm_mul_sd(a33_2, b33));
#endif
    _mm_store_sd(&C[(l_n*36)+18], c33_2);
    __m128d c33_3 = _mm_load_sd(&C[(l_n*36)+33]);
    __m128d a33_3 = _mm_load_sd(&A[99]);
#if defined(__SSE3__) && defined(__AVX__)
    c33_3 = _mm_add_sd(c33_3, _mm_mul_sd(a33_3, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c33_3 = _mm_add_sd(c33_3, _mm_mul_sd(a33_3, b33));
#endif
    _mm_store_sd(&C[(l_n*36)+33], c33_3);
#else
    C[(l_n*36)+2] += A[96] * B[(l_n*36)+33];
    C[(l_n*36)+8] += A[97] * B[(l_n*36)+33];
    C[(l_n*36)+18] += A[98] * B[(l_n*36)+33];
    C[(l_n*36)+33] += A[99] * B[(l_n*36)+33];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b34 = _mm256_broadcast_sd(&B[(l_n*36)+34]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b34 = _mm_loaddup_pd(&B[(l_n*36)+34]);
#endif
    __m128d c34_0 = _mm_load_sd(&C[(l_n*36)+0]);
    __m128d a34_0 = _mm_load_sd(&A[100]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_0 = _mm_add_sd(c34_0, _mm_mul_sd(a34_0, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_0 = _mm_add_sd(c34_0, _mm_mul_sd(a34_0, b34));
#endif
    _mm_store_sd(&C[(l_n*36)+0], c34_0);
    __m128d c34_1 = _mm_load_sd(&C[(l_n*36)+3]);
    __m128d a34_1 = _mm_load_sd(&A[101]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_1 = _mm_add_sd(c34_1, _mm_mul_sd(a34_1, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_1 = _mm_add_sd(c34_1, _mm_mul_sd(a34_1, b34));
#endif
    _mm_store_sd(&C[(l_n*36)+3], c34_1);
    __m128d c34_2 = _mm_load_sd(&C[(l_n*36)+9]);
    __m128d a34_2 = _mm_load_sd(&A[102]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_2 = _mm_add_sd(c34_2, _mm_mul_sd(a34_2, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_2 = _mm_add_sd(c34_2, _mm_mul_sd(a34_2, b34));
#endif
    _mm_store_sd(&C[(l_n*36)+9], c34_2);
    __m128d c34_3 = _mm_load_sd(&C[(l_n*36)+19]);
    __m128d a34_3 = _mm_load_sd(&A[103]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_3 = _mm_add_sd(c34_3, _mm_mul_sd(a34_3, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_3 = _mm_add_sd(c34_3, _mm_mul_sd(a34_3, b34));
#endif
    _mm_store_sd(&C[(l_n*36)+19], c34_3);
    __m128d c34_4 = _mm_load_sd(&C[(l_n*36)+34]);
    __m128d a34_4 = _mm_load_sd(&A[104]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_4 = _mm_add_sd(c34_4, _mm_mul_sd(a34_4, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_4 = _mm_add_sd(c34_4, _mm_mul_sd(a34_4, b34));
#endif
    _mm_store_sd(&C[(l_n*36)+34], c34_4);
#else
    C[(l_n*36)+0] += A[100] * B[(l_n*36)+34];
    C[(l_n*36)+3] += A[101] * B[(l_n*36)+34];
    C[(l_n*36)+9] += A[102] * B[(l_n*36)+34];
    C[(l_n*36)+19] += A[103] * B[(l_n*36)+34];
    C[(l_n*36)+34] += A[104] * B[(l_n*36)+34];
#endif

  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 1890;
#endif
}

void dsparse_fP113DivM_m35_n9_k35_ldAna5_ldB36_ldC36_beta0_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_n = 0;
  #pragma nounroll_and_jam
  for ( l_n = 0; l_n < 9; l_n++) {
    unsigned int l_m = 0;
   #pragma simd
    for ( l_m = 0; l_m < 35; l_m++) {
      C[(l_n*36)+l_m] = 0.0;
    }
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b0 = _mm256_broadcast_sd(&B[(l_n*36)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b0 = _mm_loaddup_pd(&B[(l_n*36)+0]);
#endif
    __m128d c0_0 = _mm_load_sd(&C[(l_n*36)+0]);
    __m128d a0_0 = _mm_load_sd(&A[0]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, b0));
#endif
    _mm_store_sd(&C[(l_n*36)+0], c0_0);
    __m128d c0_1 = _mm_load_sd(&C[(l_n*36)+3]);
    __m128d a0_1 = _mm_load_sd(&A[1]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, b0));
#endif
    _mm_store_sd(&C[(l_n*36)+3], c0_1);
    __m128d c0_2 = _mm_load_sd(&C[(l_n*36)+9]);
    __m128d a0_2 = _mm_load_sd(&A[2]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, b0));
#endif
    _mm_store_sd(&C[(l_n*36)+9], c0_2);
    __m128d c0_3 = _mm_load_sd(&C[(l_n*36)+19]);
    __m128d a0_3 = _mm_load_sd(&A[3]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_3 = _mm_add_sd(c0_3, _mm_mul_sd(a0_3, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_3 = _mm_add_sd(c0_3, _mm_mul_sd(a0_3, b0));
#endif
    _mm_store_sd(&C[(l_n*36)+19], c0_3);
    __m128d c0_4 = _mm_load_sd(&C[(l_n*36)+34]);
    __m128d a0_4 = _mm_load_sd(&A[4]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_4 = _mm_add_sd(c0_4, _mm_mul_sd(a0_4, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_4 = _mm_add_sd(c0_4, _mm_mul_sd(a0_4, b0));
#endif
    _mm_store_sd(&C[(l_n*36)+34], c0_4);
#else
    C[(l_n*36)+0] += A[0] * B[(l_n*36)+0];
    C[(l_n*36)+3] += A[1] * B[(l_n*36)+0];
    C[(l_n*36)+9] += A[2] * B[(l_n*36)+0];
    C[(l_n*36)+19] += A[3] * B[(l_n*36)+0];
    C[(l_n*36)+34] += A[4] * B[(l_n*36)+0];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b1 = _mm256_broadcast_sd(&B[(l_n*36)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b1 = _mm_loaddup_pd(&B[(l_n*36)+1]);
#endif
    __m128d c1_0 = _mm_load_sd(&C[(l_n*36)+1]);
    __m128d a1_0 = _mm_load_sd(&A[5]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, b1));
#endif
    _mm_store_sd(&C[(l_n*36)+1], c1_0);
    __m128d c1_1 = _mm_load_sd(&C[(l_n*36)+7]);
    __m128d a1_1 = _mm_load_sd(&A[6]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_1 = _mm_add_sd(c1_1, _mm_mul_sd(a1_1, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_1 = _mm_add_sd(c1_1, _mm_mul_sd(a1_1, b1));
#endif
    _mm_store_sd(&C[(l_n*36)+7], c1_1);
    __m128d c1_2 = _mm_load_sd(&C[(l_n*36)+17]);
    __m128d a1_2 = _mm_load_sd(&A[7]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_2 = _mm_add_sd(c1_2, _mm_mul_sd(a1_2, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_2 = _mm_add_sd(c1_2, _mm_mul_sd(a1_2, b1));
#endif
    _mm_store_sd(&C[(l_n*36)+17], c1_2);
    __m128d c1_3 = _mm_load_sd(&C[(l_n*36)+32]);
    __m128d a1_3 = _mm_load_sd(&A[8]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_3 = _mm_add_sd(c1_3, _mm_mul_sd(a1_3, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_3 = _mm_add_sd(c1_3, _mm_mul_sd(a1_3, b1));
#endif
    _mm_store_sd(&C[(l_n*36)+32], c1_3);
#else
    C[(l_n*36)+1] += A[5] * B[(l_n*36)+1];
    C[(l_n*36)+7] += A[6] * B[(l_n*36)+1];
    C[(l_n*36)+17] += A[7] * B[(l_n*36)+1];
    C[(l_n*36)+32] += A[8] * B[(l_n*36)+1];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b2 = _mm256_broadcast_sd(&B[(l_n*36)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b2 = _mm_loaddup_pd(&B[(l_n*36)+2]);
#endif
    __m128d c2_0 = _mm_load_sd(&C[(l_n*36)+2]);
    __m128d a2_0 = _mm_load_sd(&A[9]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_0 = _mm_add_sd(c2_0, _mm_mul_sd(a2_0, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_0 = _mm_add_sd(c2_0, _mm_mul_sd(a2_0, b2));
#endif
    _mm_store_sd(&C[(l_n*36)+2], c2_0);
    __m128d c2_1 = _mm_load_sd(&C[(l_n*36)+8]);
    __m128d a2_1 = _mm_load_sd(&A[10]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_1 = _mm_add_sd(c2_1, _mm_mul_sd(a2_1, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_1 = _mm_add_sd(c2_1, _mm_mul_sd(a2_1, b2));
#endif
    _mm_store_sd(&C[(l_n*36)+8], c2_1);
    __m128d c2_2 = _mm_load_sd(&C[(l_n*36)+18]);
    __m128d a2_2 = _mm_load_sd(&A[11]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_2 = _mm_add_sd(c2_2, _mm_mul_sd(a2_2, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_2 = _mm_add_sd(c2_2, _mm_mul_sd(a2_2, b2));
#endif
    _mm_store_sd(&C[(l_n*36)+18], c2_2);
    __m128d c2_3 = _mm_load_sd(&C[(l_n*36)+33]);
    __m128d a2_3 = _mm_load_sd(&A[12]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_3 = _mm_add_sd(c2_3, _mm_mul_sd(a2_3, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_3 = _mm_add_sd(c2_3, _mm_mul_sd(a2_3, b2));
#endif
    _mm_store_sd(&C[(l_n*36)+33], c2_3);
#else
    C[(l_n*36)+2] += A[9] * B[(l_n*36)+2];
    C[(l_n*36)+8] += A[10] * B[(l_n*36)+2];
    C[(l_n*36)+18] += A[11] * B[(l_n*36)+2];
    C[(l_n*36)+33] += A[12] * B[(l_n*36)+2];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b3 = _mm256_broadcast_sd(&B[(l_n*36)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b3 = _mm_loaddup_pd(&B[(l_n*36)+3]);
#endif
    __m128d c3_0 = _mm_load_sd(&C[(l_n*36)+0]);
    __m128d a3_0 = _mm_load_sd(&A[13]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, b3));
#endif
    _mm_store_sd(&C[(l_n*36)+0], c3_0);
    __m128d c3_1 = _mm_load_sd(&C[(l_n*36)+3]);
    __m128d a3_1 = _mm_load_sd(&A[14]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, b3));
#endif
    _mm_store_sd(&C[(l_n*36)+3], c3_1);
    __m128d c3_2 = _mm_load_sd(&C[(l_n*36)+9]);
    __m128d a3_2 = _mm_load_sd(&A[15]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, b3));
#endif
    _mm_store_sd(&C[(l_n*36)+9], c3_2);
    __m128d c3_3 = _mm_load_sd(&C[(l_n*36)+19]);
    __m128d a3_3 = _mm_load_sd(&A[16]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_3 = _mm_add_sd(c3_3, _mm_mul_sd(a3_3, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_3 = _mm_add_sd(c3_3, _mm_mul_sd(a3_3, b3));
#endif
    _mm_store_sd(&C[(l_n*36)+19], c3_3);
    __m128d c3_4 = _mm_load_sd(&C[(l_n*36)+34]);
    __m128d a3_4 = _mm_load_sd(&A[17]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_4 = _mm_add_sd(c3_4, _mm_mul_sd(a3_4, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_4 = _mm_add_sd(c3_4, _mm_mul_sd(a3_4, b3));
#endif
    _mm_store_sd(&C[(l_n*36)+34], c3_4);
#else
    C[(l_n*36)+0] += A[13] * B[(l_n*36)+3];
    C[(l_n*36)+3] += A[14] * B[(l_n*36)+3];
    C[(l_n*36)+9] += A[15] * B[(l_n*36)+3];
    C[(l_n*36)+19] += A[16] * B[(l_n*36)+3];
    C[(l_n*36)+34] += A[17] * B[(l_n*36)+3];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b4 = _mm256_broadcast_sd(&B[(l_n*36)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b4 = _mm_loaddup_pd(&B[(l_n*36)+4]);
#endif
    __m128d c4_0 = _mm_load_sd(&C[(l_n*36)+4]);
    __m128d a4_0 = _mm_load_sd(&A[18]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_0 = _mm_add_sd(c4_0, _mm_mul_sd(a4_0, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_0 = _mm_add_sd(c4_0, _mm_mul_sd(a4_0, b4));
#endif
    _mm_store_sd(&C[(l_n*36)+4], c4_0);
    __m128d c4_1 = _mm_load_sd(&C[(l_n*36)+14]);
    __m128d a4_1 = _mm_load_sd(&A[19]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_1 = _mm_add_sd(c4_1, _mm_mul_sd(a4_1, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_1 = _mm_add_sd(c4_1, _mm_mul_sd(a4_1, b4));
#endif
    _mm_store_sd(&C[(l_n*36)+14], c4_1);
    __m128d c4_2 = _mm_load_sd(&C[(l_n*36)+29]);
    __m128d a4_2 = _mm_load_sd(&A[20]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, b4));
#endif
    _mm_store_sd(&C[(l_n*36)+29], c4_2);
#else
    C[(l_n*36)+4] += A[18] * B[(l_n*36)+4];
    C[(l_n*36)+14] += A[19] * B[(l_n*36)+4];
    C[(l_n*36)+29] += A[20] * B[(l_n*36)+4];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b5 = _mm256_broadcast_sd(&B[(l_n*36)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b5 = _mm_loaddup_pd(&B[(l_n*36)+5]);
#endif
    __m128d c5_0 = _mm_load_sd(&C[(l_n*36)+5]);
    __m128d a5_0 = _mm_load_sd(&A[21]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_0 = _mm_add_sd(c5_0, _mm_mul_sd(a5_0, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_0 = _mm_add_sd(c5_0, _mm_mul_sd(a5_0, b5));
#endif
    _mm_store_sd(&C[(l_n*36)+5], c5_0);
    __m128d c5_1 = _mm_load_sd(&C[(l_n*36)+15]);
    __m128d a5_1 = _mm_load_sd(&A[22]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_1 = _mm_add_sd(c5_1, _mm_mul_sd(a5_1, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_1 = _mm_add_sd(c5_1, _mm_mul_sd(a5_1, b5));
#endif
    _mm_store_sd(&C[(l_n*36)+15], c5_1);
    __m128d c5_2 = _mm_load_sd(&C[(l_n*36)+30]);
    __m128d a5_2 = _mm_load_sd(&A[23]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, b5));
#endif
    _mm_store_sd(&C[(l_n*36)+30], c5_2);
#else
    C[(l_n*36)+5] += A[21] * B[(l_n*36)+5];
    C[(l_n*36)+15] += A[22] * B[(l_n*36)+5];
    C[(l_n*36)+30] += A[23] * B[(l_n*36)+5];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b6 = _mm256_broadcast_sd(&B[(l_n*36)+6]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b6 = _mm_loaddup_pd(&B[(l_n*36)+6]);
#endif
    __m128d c6_0 = _mm_load_sd(&C[(l_n*36)+6]);
    __m128d a6_0 = _mm_load_sd(&A[24]);
#if defined(__SSE3__) && defined(__AVX__)
    c6_0 = _mm_add_sd(c6_0, _mm_mul_sd(a6_0, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c6_0 = _mm_add_sd(c6_0, _mm_mul_sd(a6_0, b6));
#endif
    _mm_store_sd(&C[(l_n*36)+6], c6_0);
    __m128d c6_1 = _mm_load_sd(&C[(l_n*36)+16]);
    __m128d a6_1 = _mm_load_sd(&A[25]);
#if defined(__SSE3__) && defined(__AVX__)
    c6_1 = _mm_add_sd(c6_1, _mm_mul_sd(a6_1, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c6_1 = _mm_add_sd(c6_1, _mm_mul_sd(a6_1, b6));
#endif
    _mm_store_sd(&C[(l_n*36)+16], c6_1);
    __m128d c6_2 = _mm_load_sd(&C[(l_n*36)+31]);
    __m128d a6_2 = _mm_load_sd(&A[26]);
#if defined(__SSE3__) && defined(__AVX__)
    c6_2 = _mm_add_sd(c6_2, _mm_mul_sd(a6_2, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c6_2 = _mm_add_sd(c6_2, _mm_mul_sd(a6_2, b6));
#endif
    _mm_store_sd(&C[(l_n*36)+31], c6_2);
#else
    C[(l_n*36)+6] += A[24] * B[(l_n*36)+6];
    C[(l_n*36)+16] += A[25] * B[(l_n*36)+6];
    C[(l_n*36)+31] += A[26] * B[(l_n*36)+6];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b7 = _mm256_broadcast_sd(&B[(l_n*36)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b7 = _mm_loaddup_pd(&B[(l_n*36)+7]);
#endif
    __m128d c7_0 = _mm_load_sd(&C[(l_n*36)+1]);
    __m128d a7_0 = _mm_load_sd(&A[27]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_0 = _mm_add_sd(c7_0, _mm_mul_sd(a7_0, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_0 = _mm_add_sd(c7_0, _mm_mul_sd(a7_0, b7));
#endif
    _mm_store_sd(&C[(l_n*36)+1], c7_0);
    __m128d c7_1 = _mm_load_sd(&C[(l_n*36)+7]);
    __m128d a7_1 = _mm_load_sd(&A[28]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_1 = _mm_add_sd(c7_1, _mm_mul_sd(a7_1, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_1 = _mm_add_sd(c7_1, _mm_mul_sd(a7_1, b7));
#endif
    _mm_store_sd(&C[(l_n*36)+7], c7_1);
    __m128d c7_2 = _mm_load_sd(&C[(l_n*36)+17]);
    __m128d a7_2 = _mm_load_sd(&A[29]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_2 = _mm_add_sd(c7_2, _mm_mul_sd(a7_2, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_2 = _mm_add_sd(c7_2, _mm_mul_sd(a7_2, b7));
#endif
    _mm_store_sd(&C[(l_n*36)+17], c7_2);
    __m128d c7_3 = _mm_load_sd(&C[(l_n*36)+32]);
    __m128d a7_3 = _mm_load_sd(&A[30]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_3 = _mm_add_sd(c7_3, _mm_mul_sd(a7_3, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_3 = _mm_add_sd(c7_3, _mm_mul_sd(a7_3, b7));
#endif
    _mm_store_sd(&C[(l_n*36)+32], c7_3);
#else
    C[(l_n*36)+1] += A[27] * B[(l_n*36)+7];
    C[(l_n*36)+7] += A[28] * B[(l_n*36)+7];
    C[(l_n*36)+17] += A[29] * B[(l_n*36)+7];
    C[(l_n*36)+32] += A[30] * B[(l_n*36)+7];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b8 = _mm256_broadcast_sd(&B[(l_n*36)+8]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b8 = _mm_loaddup_pd(&B[(l_n*36)+8]);
#endif
    __m128d c8_0 = _mm_load_sd(&C[(l_n*36)+2]);
    __m128d a8_0 = _mm_load_sd(&A[31]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_0 = _mm_add_sd(c8_0, _mm_mul_sd(a8_0, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_0 = _mm_add_sd(c8_0, _mm_mul_sd(a8_0, b8));
#endif
    _mm_store_sd(&C[(l_n*36)+2], c8_0);
    __m128d c8_1 = _mm_load_sd(&C[(l_n*36)+8]);
    __m128d a8_1 = _mm_load_sd(&A[32]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_1 = _mm_add_sd(c8_1, _mm_mul_sd(a8_1, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_1 = _mm_add_sd(c8_1, _mm_mul_sd(a8_1, b8));
#endif
    _mm_store_sd(&C[(l_n*36)+8], c8_1);
    __m128d c8_2 = _mm_load_sd(&C[(l_n*36)+18]);
    __m128d a8_2 = _mm_load_sd(&A[33]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_2 = _mm_add_sd(c8_2, _mm_mul_sd(a8_2, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_2 = _mm_add_sd(c8_2, _mm_mul_sd(a8_2, b8));
#endif
    _mm_store_sd(&C[(l_n*36)+18], c8_2);
    __m128d c8_3 = _mm_load_sd(&C[(l_n*36)+33]);
    __m128d a8_3 = _mm_load_sd(&A[34]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_3 = _mm_add_sd(c8_3, _mm_mul_sd(a8_3, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_3 = _mm_add_sd(c8_3, _mm_mul_sd(a8_3, b8));
#endif
    _mm_store_sd(&C[(l_n*36)+33], c8_3);
#else
    C[(l_n*36)+2] += A[31] * B[(l_n*36)+8];
    C[(l_n*36)+8] += A[32] * B[(l_n*36)+8];
    C[(l_n*36)+18] += A[33] * B[(l_n*36)+8];
    C[(l_n*36)+33] += A[34] * B[(l_n*36)+8];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b9 = _mm256_broadcast_sd(&B[(l_n*36)+9]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b9 = _mm_loaddup_pd(&B[(l_n*36)+9]);
#endif
    __m128d c9_0 = _mm_load_sd(&C[(l_n*36)+0]);
    __m128d a9_0 = _mm_load_sd(&A[35]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, b9));
#endif
    _mm_store_sd(&C[(l_n*36)+0], c9_0);
    __m128d c9_1 = _mm_load_sd(&C[(l_n*36)+3]);
    __m128d a9_1 = _mm_load_sd(&A[36]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, b9));
#endif
    _mm_store_sd(&C[(l_n*36)+3], c9_1);
    __m128d c9_2 = _mm_load_sd(&C[(l_n*36)+9]);
    __m128d a9_2 = _mm_load_sd(&A[37]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, b9));
#endif
    _mm_store_sd(&C[(l_n*36)+9], c9_2);
    __m128d c9_3 = _mm_load_sd(&C[(l_n*36)+19]);
    __m128d a9_3 = _mm_load_sd(&A[38]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_3 = _mm_add_sd(c9_3, _mm_mul_sd(a9_3, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_3 = _mm_add_sd(c9_3, _mm_mul_sd(a9_3, b9));
#endif
    _mm_store_sd(&C[(l_n*36)+19], c9_3);
    __m128d c9_4 = _mm_load_sd(&C[(l_n*36)+34]);
    __m128d a9_4 = _mm_load_sd(&A[39]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_4 = _mm_add_sd(c9_4, _mm_mul_sd(a9_4, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_4 = _mm_add_sd(c9_4, _mm_mul_sd(a9_4, b9));
#endif
    _mm_store_sd(&C[(l_n*36)+34], c9_4);
#else
    C[(l_n*36)+0] += A[35] * B[(l_n*36)+9];
    C[(l_n*36)+3] += A[36] * B[(l_n*36)+9];
    C[(l_n*36)+9] += A[37] * B[(l_n*36)+9];
    C[(l_n*36)+19] += A[38] * B[(l_n*36)+9];
    C[(l_n*36)+34] += A[39] * B[(l_n*36)+9];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b10 = _mm256_broadcast_sd(&B[(l_n*36)+10]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b10 = _mm_loaddup_pd(&B[(l_n*36)+10]);
#endif
    __m128d c10_0 = _mm_load_sd(&C[(l_n*36)+10]);
    __m128d a10_0 = _mm_load_sd(&A[40]);
#if defined(__SSE3__) && defined(__AVX__)
    c10_0 = _mm_add_sd(c10_0, _mm_mul_sd(a10_0, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c10_0 = _mm_add_sd(c10_0, _mm_mul_sd(a10_0, b10));
#endif
    _mm_store_sd(&C[(l_n*36)+10], c10_0);
    __m128d c10_1 = _mm_load_sd(&C[(l_n*36)+25]);
    __m128d a10_1 = _mm_load_sd(&A[41]);
#if defined(__SSE3__) && defined(__AVX__)
    c10_1 = _mm_add_sd(c10_1, _mm_mul_sd(a10_1, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c10_1 = _mm_add_sd(c10_1, _mm_mul_sd(a10_1, b10));
#endif
    _mm_store_sd(&C[(l_n*36)+25], c10_1);
#else
    C[(l_n*36)+10] += A[40] * B[(l_n*36)+10];
    C[(l_n*36)+25] += A[41] * B[(l_n*36)+10];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b11 = _mm256_broadcast_sd(&B[(l_n*36)+11]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b11 = _mm_loaddup_pd(&B[(l_n*36)+11]);
#endif
    __m128d c11_0 = _mm_load_sd(&C[(l_n*36)+11]);
    __m128d a11_0 = _mm_load_sd(&A[42]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_0 = _mm_add_sd(c11_0, _mm_mul_sd(a11_0, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_0 = _mm_add_sd(c11_0, _mm_mul_sd(a11_0, b11));
#endif
    _mm_store_sd(&C[(l_n*36)+11], c11_0);
    __m128d c11_1 = _mm_load_sd(&C[(l_n*36)+26]);
    __m128d a11_1 = _mm_load_sd(&A[43]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_1 = _mm_add_sd(c11_1, _mm_mul_sd(a11_1, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_1 = _mm_add_sd(c11_1, _mm_mul_sd(a11_1, b11));
#endif
    _mm_store_sd(&C[(l_n*36)+26], c11_1);
#else
    C[(l_n*36)+11] += A[42] * B[(l_n*36)+11];
    C[(l_n*36)+26] += A[43] * B[(l_n*36)+11];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b12 = _mm256_broadcast_sd(&B[(l_n*36)+12]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b12 = _mm_loaddup_pd(&B[(l_n*36)+12]);
#endif
    __m128d c12_0 = _mm_load_sd(&C[(l_n*36)+12]);
    __m128d a12_0 = _mm_load_sd(&A[44]);
#if defined(__SSE3__) && defined(__AVX__)
    c12_0 = _mm_add_sd(c12_0, _mm_mul_sd(a12_0, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c12_0 = _mm_add_sd(c12_0, _mm_mul_sd(a12_0, b12));
#endif
    _mm_store_sd(&C[(l_n*36)+12], c12_0);
    __m128d c12_1 = _mm_load_sd(&C[(l_n*36)+27]);
    __m128d a12_1 = _mm_load_sd(&A[45]);
#if defined(__SSE3__) && defined(__AVX__)
    c12_1 = _mm_add_sd(c12_1, _mm_mul_sd(a12_1, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c12_1 = _mm_add_sd(c12_1, _mm_mul_sd(a12_1, b12));
#endif
    _mm_store_sd(&C[(l_n*36)+27], c12_1);
#else
    C[(l_n*36)+12] += A[44] * B[(l_n*36)+12];
    C[(l_n*36)+27] += A[45] * B[(l_n*36)+12];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b13 = _mm256_broadcast_sd(&B[(l_n*36)+13]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b13 = _mm_loaddup_pd(&B[(l_n*36)+13]);
#endif
    __m128d c13_0 = _mm_load_sd(&C[(l_n*36)+13]);
    __m128d a13_0 = _mm_load_sd(&A[46]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_0 = _mm_add_sd(c13_0, _mm_mul_sd(a13_0, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_0 = _mm_add_sd(c13_0, _mm_mul_sd(a13_0, b13));
#endif
    _mm_store_sd(&C[(l_n*36)+13], c13_0);
    __m128d c13_1 = _mm_load_sd(&C[(l_n*36)+28]);
    __m128d a13_1 = _mm_load_sd(&A[47]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_1 = _mm_add_sd(c13_1, _mm_mul_sd(a13_1, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_1 = _mm_add_sd(c13_1, _mm_mul_sd(a13_1, b13));
#endif
    _mm_store_sd(&C[(l_n*36)+28], c13_1);
#else
    C[(l_n*36)+13] += A[46] * B[(l_n*36)+13];
    C[(l_n*36)+28] += A[47] * B[(l_n*36)+13];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b14 = _mm256_broadcast_sd(&B[(l_n*36)+14]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b14 = _mm_loaddup_pd(&B[(l_n*36)+14]);
#endif
    __m128d c14_0 = _mm_load_sd(&C[(l_n*36)+4]);
    __m128d a14_0 = _mm_load_sd(&A[48]);
#if defined(__SSE3__) && defined(__AVX__)
    c14_0 = _mm_add_sd(c14_0, _mm_mul_sd(a14_0, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c14_0 = _mm_add_sd(c14_0, _mm_mul_sd(a14_0, b14));
#endif
    _mm_store_sd(&C[(l_n*36)+4], c14_0);
    __m128d c14_1 = _mm_load_sd(&C[(l_n*36)+14]);
    __m128d a14_1 = _mm_load_sd(&A[49]);
#if defined(__SSE3__) && defined(__AVX__)
    c14_1 = _mm_add_sd(c14_1, _mm_mul_sd(a14_1, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c14_1 = _mm_add_sd(c14_1, _mm_mul_sd(a14_1, b14));
#endif
    _mm_store_sd(&C[(l_n*36)+14], c14_1);
    __m128d c14_2 = _mm_load_sd(&C[(l_n*36)+29]);
    __m128d a14_2 = _mm_load_sd(&A[50]);
#if defined(__SSE3__) && defined(__AVX__)
    c14_2 = _mm_add_sd(c14_2, _mm_mul_sd(a14_2, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c14_2 = _mm_add_sd(c14_2, _mm_mul_sd(a14_2, b14));
#endif
    _mm_store_sd(&C[(l_n*36)+29], c14_2);
#else
    C[(l_n*36)+4] += A[48] * B[(l_n*36)+14];
    C[(l_n*36)+14] += A[49] * B[(l_n*36)+14];
    C[(l_n*36)+29] += A[50] * B[(l_n*36)+14];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b15 = _mm256_broadcast_sd(&B[(l_n*36)+15]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b15 = _mm_loaddup_pd(&B[(l_n*36)+15]);
#endif
    __m128d c15_0 = _mm_load_sd(&C[(l_n*36)+5]);
    __m128d a15_0 = _mm_load_sd(&A[51]);
#if defined(__SSE3__) && defined(__AVX__)
    c15_0 = _mm_add_sd(c15_0, _mm_mul_sd(a15_0, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c15_0 = _mm_add_sd(c15_0, _mm_mul_sd(a15_0, b15));
#endif
    _mm_store_sd(&C[(l_n*36)+5], c15_0);
    __m128d c15_1 = _mm_load_sd(&C[(l_n*36)+15]);
    __m128d a15_1 = _mm_load_sd(&A[52]);
#if defined(__SSE3__) && defined(__AVX__)
    c15_1 = _mm_add_sd(c15_1, _mm_mul_sd(a15_1, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c15_1 = _mm_add_sd(c15_1, _mm_mul_sd(a15_1, b15));
#endif
    _mm_store_sd(&C[(l_n*36)+15], c15_1);
    __m128d c15_2 = _mm_load_sd(&C[(l_n*36)+30]);
    __m128d a15_2 = _mm_load_sd(&A[53]);
#if defined(__SSE3__) && defined(__AVX__)
    c15_2 = _mm_add_sd(c15_2, _mm_mul_sd(a15_2, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c15_2 = _mm_add_sd(c15_2, _mm_mul_sd(a15_2, b15));
#endif
    _mm_store_sd(&C[(l_n*36)+30], c15_2);
#else
    C[(l_n*36)+5] += A[51] * B[(l_n*36)+15];
    C[(l_n*36)+15] += A[52] * B[(l_n*36)+15];
    C[(l_n*36)+30] += A[53] * B[(l_n*36)+15];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b16 = _mm256_broadcast_sd(&B[(l_n*36)+16]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b16 = _mm_loaddup_pd(&B[(l_n*36)+16]);
#endif
    __m128d c16_0 = _mm_load_sd(&C[(l_n*36)+6]);
    __m128d a16_0 = _mm_load_sd(&A[54]);
#if defined(__SSE3__) && defined(__AVX__)
    c16_0 = _mm_add_sd(c16_0, _mm_mul_sd(a16_0, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c16_0 = _mm_add_sd(c16_0, _mm_mul_sd(a16_0, b16));
#endif
    _mm_store_sd(&C[(l_n*36)+6], c16_0);
    __m128d c16_1 = _mm_load_sd(&C[(l_n*36)+16]);
    __m128d a16_1 = _mm_load_sd(&A[55]);
#if defined(__SSE3__) && defined(__AVX__)
    c16_1 = _mm_add_sd(c16_1, _mm_mul_sd(a16_1, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c16_1 = _mm_add_sd(c16_1, _mm_mul_sd(a16_1, b16));
#endif
    _mm_store_sd(&C[(l_n*36)+16], c16_1);
    __m128d c16_2 = _mm_load_sd(&C[(l_n*36)+31]);
    __m128d a16_2 = _mm_load_sd(&A[56]);
#if defined(__SSE3__) && defined(__AVX__)
    c16_2 = _mm_add_sd(c16_2, _mm_mul_sd(a16_2, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c16_2 = _mm_add_sd(c16_2, _mm_mul_sd(a16_2, b16));
#endif
    _mm_store_sd(&C[(l_n*36)+31], c16_2);
#else
    C[(l_n*36)+6] += A[54] * B[(l_n*36)+16];
    C[(l_n*36)+16] += A[55] * B[(l_n*36)+16];
    C[(l_n*36)+31] += A[56] * B[(l_n*36)+16];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b17 = _mm256_broadcast_sd(&B[(l_n*36)+17]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b17 = _mm_loaddup_pd(&B[(l_n*36)+17]);
#endif
    __m128d c17_0 = _mm_load_sd(&C[(l_n*36)+1]);
    __m128d a17_0 = _mm_load_sd(&A[57]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_0 = _mm_add_sd(c17_0, _mm_mul_sd(a17_0, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_0 = _mm_add_sd(c17_0, _mm_mul_sd(a17_0, b17));
#endif
    _mm_store_sd(&C[(l_n*36)+1], c17_0);
    __m128d c17_1 = _mm_load_sd(&C[(l_n*36)+7]);
    __m128d a17_1 = _mm_load_sd(&A[58]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_1 = _mm_add_sd(c17_1, _mm_mul_sd(a17_1, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_1 = _mm_add_sd(c17_1, _mm_mul_sd(a17_1, b17));
#endif
    _mm_store_sd(&C[(l_n*36)+7], c17_1);
    __m128d c17_2 = _mm_load_sd(&C[(l_n*36)+17]);
    __m128d a17_2 = _mm_load_sd(&A[59]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_2 = _mm_add_sd(c17_2, _mm_mul_sd(a17_2, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_2 = _mm_add_sd(c17_2, _mm_mul_sd(a17_2, b17));
#endif
    _mm_store_sd(&C[(l_n*36)+17], c17_2);
    __m128d c17_3 = _mm_load_sd(&C[(l_n*36)+32]);
    __m128d a17_3 = _mm_load_sd(&A[60]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_3 = _mm_add_sd(c17_3, _mm_mul_sd(a17_3, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_3 = _mm_add_sd(c17_3, _mm_mul_sd(a17_3, b17));
#endif
    _mm_store_sd(&C[(l_n*36)+32], c17_3);
#else
    C[(l_n*36)+1] += A[57] * B[(l_n*36)+17];
    C[(l_n*36)+7] += A[58] * B[(l_n*36)+17];
    C[(l_n*36)+17] += A[59] * B[(l_n*36)+17];
    C[(l_n*36)+32] += A[60] * B[(l_n*36)+17];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b18 = _mm256_broadcast_sd(&B[(l_n*36)+18]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b18 = _mm_loaddup_pd(&B[(l_n*36)+18]);
#endif
    __m128d c18_0 = _mm_load_sd(&C[(l_n*36)+2]);
    __m128d a18_0 = _mm_load_sd(&A[61]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_0 = _mm_add_sd(c18_0, _mm_mul_sd(a18_0, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_0 = _mm_add_sd(c18_0, _mm_mul_sd(a18_0, b18));
#endif
    _mm_store_sd(&C[(l_n*36)+2], c18_0);
    __m128d c18_1 = _mm_load_sd(&C[(l_n*36)+8]);
    __m128d a18_1 = _mm_load_sd(&A[62]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_1 = _mm_add_sd(c18_1, _mm_mul_sd(a18_1, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_1 = _mm_add_sd(c18_1, _mm_mul_sd(a18_1, b18));
#endif
    _mm_store_sd(&C[(l_n*36)+8], c18_1);
    __m128d c18_2 = _mm_load_sd(&C[(l_n*36)+18]);
    __m128d a18_2 = _mm_load_sd(&A[63]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_2 = _mm_add_sd(c18_2, _mm_mul_sd(a18_2, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_2 = _mm_add_sd(c18_2, _mm_mul_sd(a18_2, b18));
#endif
    _mm_store_sd(&C[(l_n*36)+18], c18_2);
    __m128d c18_3 = _mm_load_sd(&C[(l_n*36)+33]);
    __m128d a18_3 = _mm_load_sd(&A[64]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_3 = _mm_add_sd(c18_3, _mm_mul_sd(a18_3, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_3 = _mm_add_sd(c18_3, _mm_mul_sd(a18_3, b18));
#endif
    _mm_store_sd(&C[(l_n*36)+33], c18_3);
#else
    C[(l_n*36)+2] += A[61] * B[(l_n*36)+18];
    C[(l_n*36)+8] += A[62] * B[(l_n*36)+18];
    C[(l_n*36)+18] += A[63] * B[(l_n*36)+18];
    C[(l_n*36)+33] += A[64] * B[(l_n*36)+18];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b19 = _mm256_broadcast_sd(&B[(l_n*36)+19]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b19 = _mm_loaddup_pd(&B[(l_n*36)+19]);
#endif
    __m128d c19_0 = _mm_load_sd(&C[(l_n*36)+0]);
    __m128d a19_0 = _mm_load_sd(&A[65]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_0 = _mm_add_sd(c19_0, _mm_mul_sd(a19_0, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_0 = _mm_add_sd(c19_0, _mm_mul_sd(a19_0, b19));
#endif
    _mm_store_sd(&C[(l_n*36)+0], c19_0);
    __m128d c19_1 = _mm_load_sd(&C[(l_n*36)+3]);
    __m128d a19_1 = _mm_load_sd(&A[66]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_1 = _mm_add_sd(c19_1, _mm_mul_sd(a19_1, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_1 = _mm_add_sd(c19_1, _mm_mul_sd(a19_1, b19));
#endif
    _mm_store_sd(&C[(l_n*36)+3], c19_1);
    __m128d c19_2 = _mm_load_sd(&C[(l_n*36)+9]);
    __m128d a19_2 = _mm_load_sd(&A[67]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_2 = _mm_add_sd(c19_2, _mm_mul_sd(a19_2, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_2 = _mm_add_sd(c19_2, _mm_mul_sd(a19_2, b19));
#endif
    _mm_store_sd(&C[(l_n*36)+9], c19_2);
    __m128d c19_3 = _mm_load_sd(&C[(l_n*36)+19]);
    __m128d a19_3 = _mm_load_sd(&A[68]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_3 = _mm_add_sd(c19_3, _mm_mul_sd(a19_3, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_3 = _mm_add_sd(c19_3, _mm_mul_sd(a19_3, b19));
#endif
    _mm_store_sd(&C[(l_n*36)+19], c19_3);
    __m128d c19_4 = _mm_load_sd(&C[(l_n*36)+34]);
    __m128d a19_4 = _mm_load_sd(&A[69]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_4 = _mm_add_sd(c19_4, _mm_mul_sd(a19_4, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_4 = _mm_add_sd(c19_4, _mm_mul_sd(a19_4, b19));
#endif
    _mm_store_sd(&C[(l_n*36)+34], c19_4);
#else
    C[(l_n*36)+0] += A[65] * B[(l_n*36)+19];
    C[(l_n*36)+3] += A[66] * B[(l_n*36)+19];
    C[(l_n*36)+9] += A[67] * B[(l_n*36)+19];
    C[(l_n*36)+19] += A[68] * B[(l_n*36)+19];
    C[(l_n*36)+34] += A[69] * B[(l_n*36)+19];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b20 = _mm256_broadcast_sd(&B[(l_n*36)+20]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b20 = _mm_loaddup_pd(&B[(l_n*36)+20]);
#endif
    __m128d c20_0 = _mm_load_sd(&C[(l_n*36)+20]);
    __m128d a20_0 = _mm_load_sd(&A[70]);
#if defined(__SSE3__) && defined(__AVX__)
    c20_0 = _mm_add_sd(c20_0, _mm_mul_sd(a20_0, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c20_0 = _mm_add_sd(c20_0, _mm_mul_sd(a20_0, b20));
#endif
    _mm_store_sd(&C[(l_n*36)+20], c20_0);
#else
    C[(l_n*36)+20] += A[70] * B[(l_n*36)+20];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b21 = _mm256_broadcast_sd(&B[(l_n*36)+21]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b21 = _mm_loaddup_pd(&B[(l_n*36)+21]);
#endif
    __m128d c21_0 = _mm_load_sd(&C[(l_n*36)+21]);
    __m128d a21_0 = _mm_load_sd(&A[71]);
#if defined(__SSE3__) && defined(__AVX__)
    c21_0 = _mm_add_sd(c21_0, _mm_mul_sd(a21_0, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c21_0 = _mm_add_sd(c21_0, _mm_mul_sd(a21_0, b21));
#endif
    _mm_store_sd(&C[(l_n*36)+21], c21_0);
#else
    C[(l_n*36)+21] += A[71] * B[(l_n*36)+21];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b22 = _mm256_broadcast_sd(&B[(l_n*36)+22]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b22 = _mm_loaddup_pd(&B[(l_n*36)+22]);
#endif
    __m128d c22_0 = _mm_load_sd(&C[(l_n*36)+22]);
    __m128d a22_0 = _mm_load_sd(&A[72]);
#if defined(__SSE3__) && defined(__AVX__)
    c22_0 = _mm_add_sd(c22_0, _mm_mul_sd(a22_0, _mm256_castpd256_pd128(b22)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c22_0 = _mm_add_sd(c22_0, _mm_mul_sd(a22_0, b22));
#endif
    _mm_store_sd(&C[(l_n*36)+22], c22_0);
#else
    C[(l_n*36)+22] += A[72] * B[(l_n*36)+22];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b23 = _mm256_broadcast_sd(&B[(l_n*36)+23]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b23 = _mm_loaddup_pd(&B[(l_n*36)+23]);
#endif
    __m128d c23_0 = _mm_load_sd(&C[(l_n*36)+23]);
    __m128d a23_0 = _mm_load_sd(&A[73]);
#if defined(__SSE3__) && defined(__AVX__)
    c23_0 = _mm_add_sd(c23_0, _mm_mul_sd(a23_0, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c23_0 = _mm_add_sd(c23_0, _mm_mul_sd(a23_0, b23));
#endif
    _mm_store_sd(&C[(l_n*36)+23], c23_0);
#else
    C[(l_n*36)+23] += A[73] * B[(l_n*36)+23];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b24 = _mm256_broadcast_sd(&B[(l_n*36)+24]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b24 = _mm_loaddup_pd(&B[(l_n*36)+24]);
#endif
    __m128d c24_0 = _mm_load_sd(&C[(l_n*36)+24]);
    __m128d a24_0 = _mm_load_sd(&A[74]);
#if defined(__SSE3__) && defined(__AVX__)
    c24_0 = _mm_add_sd(c24_0, _mm_mul_sd(a24_0, _mm256_castpd256_pd128(b24)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c24_0 = _mm_add_sd(c24_0, _mm_mul_sd(a24_0, b24));
#endif
    _mm_store_sd(&C[(l_n*36)+24], c24_0);
#else
    C[(l_n*36)+24] += A[74] * B[(l_n*36)+24];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b25 = _mm256_broadcast_sd(&B[(l_n*36)+25]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b25 = _mm_loaddup_pd(&B[(l_n*36)+25]);
#endif
    __m128d c25_0 = _mm_load_sd(&C[(l_n*36)+10]);
    __m128d a25_0 = _mm_load_sd(&A[75]);
#if defined(__SSE3__) && defined(__AVX__)
    c25_0 = _mm_add_sd(c25_0, _mm_mul_sd(a25_0, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c25_0 = _mm_add_sd(c25_0, _mm_mul_sd(a25_0, b25));
#endif
    _mm_store_sd(&C[(l_n*36)+10], c25_0);
    __m128d c25_1 = _mm_load_sd(&C[(l_n*36)+25]);
    __m128d a25_1 = _mm_load_sd(&A[76]);
#if defined(__SSE3__) && defined(__AVX__)
    c25_1 = _mm_add_sd(c25_1, _mm_mul_sd(a25_1, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c25_1 = _mm_add_sd(c25_1, _mm_mul_sd(a25_1, b25));
#endif
    _mm_store_sd(&C[(l_n*36)+25], c25_1);
#else
    C[(l_n*36)+10] += A[75] * B[(l_n*36)+25];
    C[(l_n*36)+25] += A[76] * B[(l_n*36)+25];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b26 = _mm256_broadcast_sd(&B[(l_n*36)+26]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b26 = _mm_loaddup_pd(&B[(l_n*36)+26]);
#endif
    __m128d c26_0 = _mm_load_sd(&C[(l_n*36)+11]);
    __m128d a26_0 = _mm_load_sd(&A[77]);
#if defined(__SSE3__) && defined(__AVX__)
    c26_0 = _mm_add_sd(c26_0, _mm_mul_sd(a26_0, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c26_0 = _mm_add_sd(c26_0, _mm_mul_sd(a26_0, b26));
#endif
    _mm_store_sd(&C[(l_n*36)+11], c26_0);
    __m128d c26_1 = _mm_load_sd(&C[(l_n*36)+26]);
    __m128d a26_1 = _mm_load_sd(&A[78]);
#if defined(__SSE3__) && defined(__AVX__)
    c26_1 = _mm_add_sd(c26_1, _mm_mul_sd(a26_1, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c26_1 = _mm_add_sd(c26_1, _mm_mul_sd(a26_1, b26));
#endif
    _mm_store_sd(&C[(l_n*36)+26], c26_1);
#else
    C[(l_n*36)+11] += A[77] * B[(l_n*36)+26];
    C[(l_n*36)+26] += A[78] * B[(l_n*36)+26];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b27 = _mm256_broadcast_sd(&B[(l_n*36)+27]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b27 = _mm_loaddup_pd(&B[(l_n*36)+27]);
#endif
    __m128d c27_0 = _mm_load_sd(&C[(l_n*36)+12]);
    __m128d a27_0 = _mm_load_sd(&A[79]);
#if defined(__SSE3__) && defined(__AVX__)
    c27_0 = _mm_add_sd(c27_0, _mm_mul_sd(a27_0, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c27_0 = _mm_add_sd(c27_0, _mm_mul_sd(a27_0, b27));
#endif
    _mm_store_sd(&C[(l_n*36)+12], c27_0);
    __m128d c27_1 = _mm_load_sd(&C[(l_n*36)+27]);
    __m128d a27_1 = _mm_load_sd(&A[80]);
#if defined(__SSE3__) && defined(__AVX__)
    c27_1 = _mm_add_sd(c27_1, _mm_mul_sd(a27_1, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c27_1 = _mm_add_sd(c27_1, _mm_mul_sd(a27_1, b27));
#endif
    _mm_store_sd(&C[(l_n*36)+27], c27_1);
#else
    C[(l_n*36)+12] += A[79] * B[(l_n*36)+27];
    C[(l_n*36)+27] += A[80] * B[(l_n*36)+27];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b28 = _mm256_broadcast_sd(&B[(l_n*36)+28]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b28 = _mm_loaddup_pd(&B[(l_n*36)+28]);
#endif
    __m128d c28_0 = _mm_load_sd(&C[(l_n*36)+13]);
    __m128d a28_0 = _mm_load_sd(&A[81]);
#if defined(__SSE3__) && defined(__AVX__)
    c28_0 = _mm_add_sd(c28_0, _mm_mul_sd(a28_0, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c28_0 = _mm_add_sd(c28_0, _mm_mul_sd(a28_0, b28));
#endif
    _mm_store_sd(&C[(l_n*36)+13], c28_0);
    __m128d c28_1 = _mm_load_sd(&C[(l_n*36)+28]);
    __m128d a28_1 = _mm_load_sd(&A[82]);
#if defined(__SSE3__) && defined(__AVX__)
    c28_1 = _mm_add_sd(c28_1, _mm_mul_sd(a28_1, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c28_1 = _mm_add_sd(c28_1, _mm_mul_sd(a28_1, b28));
#endif
    _mm_store_sd(&C[(l_n*36)+28], c28_1);
#else
    C[(l_n*36)+13] += A[81] * B[(l_n*36)+28];
    C[(l_n*36)+28] += A[82] * B[(l_n*36)+28];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b29 = _mm256_broadcast_sd(&B[(l_n*36)+29]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b29 = _mm_loaddup_pd(&B[(l_n*36)+29]);
#endif
    __m128d c29_0 = _mm_load_sd(&C[(l_n*36)+4]);
    __m128d a29_0 = _mm_load_sd(&A[83]);
#if defined(__SSE3__) && defined(__AVX__)
    c29_0 = _mm_add_sd(c29_0, _mm_mul_sd(a29_0, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c29_0 = _mm_add_sd(c29_0, _mm_mul_sd(a29_0, b29));
#endif
    _mm_store_sd(&C[(l_n*36)+4], c29_0);
    __m128d c29_1 = _mm_load_sd(&C[(l_n*36)+14]);
    __m128d a29_1 = _mm_load_sd(&A[84]);
#if defined(__SSE3__) && defined(__AVX__)
    c29_1 = _mm_add_sd(c29_1, _mm_mul_sd(a29_1, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c29_1 = _mm_add_sd(c29_1, _mm_mul_sd(a29_1, b29));
#endif
    _mm_store_sd(&C[(l_n*36)+14], c29_1);
    __m128d c29_2 = _mm_load_sd(&C[(l_n*36)+29]);
    __m128d a29_2 = _mm_load_sd(&A[85]);
#if defined(__SSE3__) && defined(__AVX__)
    c29_2 = _mm_add_sd(c29_2, _mm_mul_sd(a29_2, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c29_2 = _mm_add_sd(c29_2, _mm_mul_sd(a29_2, b29));
#endif
    _mm_store_sd(&C[(l_n*36)+29], c29_2);
#else
    C[(l_n*36)+4] += A[83] * B[(l_n*36)+29];
    C[(l_n*36)+14] += A[84] * B[(l_n*36)+29];
    C[(l_n*36)+29] += A[85] * B[(l_n*36)+29];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b30 = _mm256_broadcast_sd(&B[(l_n*36)+30]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b30 = _mm_loaddup_pd(&B[(l_n*36)+30]);
#endif
    __m128d c30_0 = _mm_load_sd(&C[(l_n*36)+5]);
    __m128d a30_0 = _mm_load_sd(&A[86]);
#if defined(__SSE3__) && defined(__AVX__)
    c30_0 = _mm_add_sd(c30_0, _mm_mul_sd(a30_0, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c30_0 = _mm_add_sd(c30_0, _mm_mul_sd(a30_0, b30));
#endif
    _mm_store_sd(&C[(l_n*36)+5], c30_0);
    __m128d c30_1 = _mm_load_sd(&C[(l_n*36)+15]);
    __m128d a30_1 = _mm_load_sd(&A[87]);
#if defined(__SSE3__) && defined(__AVX__)
    c30_1 = _mm_add_sd(c30_1, _mm_mul_sd(a30_1, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c30_1 = _mm_add_sd(c30_1, _mm_mul_sd(a30_1, b30));
#endif
    _mm_store_sd(&C[(l_n*36)+15], c30_1);
    __m128d c30_2 = _mm_load_sd(&C[(l_n*36)+30]);
    __m128d a30_2 = _mm_load_sd(&A[88]);
#if defined(__SSE3__) && defined(__AVX__)
    c30_2 = _mm_add_sd(c30_2, _mm_mul_sd(a30_2, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c30_2 = _mm_add_sd(c30_2, _mm_mul_sd(a30_2, b30));
#endif
    _mm_store_sd(&C[(l_n*36)+30], c30_2);
#else
    C[(l_n*36)+5] += A[86] * B[(l_n*36)+30];
    C[(l_n*36)+15] += A[87] * B[(l_n*36)+30];
    C[(l_n*36)+30] += A[88] * B[(l_n*36)+30];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b31 = _mm256_broadcast_sd(&B[(l_n*36)+31]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b31 = _mm_loaddup_pd(&B[(l_n*36)+31]);
#endif
    __m128d c31_0 = _mm_load_sd(&C[(l_n*36)+6]);
    __m128d a31_0 = _mm_load_sd(&A[89]);
#if defined(__SSE3__) && defined(__AVX__)
    c31_0 = _mm_add_sd(c31_0, _mm_mul_sd(a31_0, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c31_0 = _mm_add_sd(c31_0, _mm_mul_sd(a31_0, b31));
#endif
    _mm_store_sd(&C[(l_n*36)+6], c31_0);
    __m128d c31_1 = _mm_load_sd(&C[(l_n*36)+16]);
    __m128d a31_1 = _mm_load_sd(&A[90]);
#if defined(__SSE3__) && defined(__AVX__)
    c31_1 = _mm_add_sd(c31_1, _mm_mul_sd(a31_1, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c31_1 = _mm_add_sd(c31_1, _mm_mul_sd(a31_1, b31));
#endif
    _mm_store_sd(&C[(l_n*36)+16], c31_1);
    __m128d c31_2 = _mm_load_sd(&C[(l_n*36)+31]);
    __m128d a31_2 = _mm_load_sd(&A[91]);
#if defined(__SSE3__) && defined(__AVX__)
    c31_2 = _mm_add_sd(c31_2, _mm_mul_sd(a31_2, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c31_2 = _mm_add_sd(c31_2, _mm_mul_sd(a31_2, b31));
#endif
    _mm_store_sd(&C[(l_n*36)+31], c31_2);
#else
    C[(l_n*36)+6] += A[89] * B[(l_n*36)+31];
    C[(l_n*36)+16] += A[90] * B[(l_n*36)+31];
    C[(l_n*36)+31] += A[91] * B[(l_n*36)+31];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b32 = _mm256_broadcast_sd(&B[(l_n*36)+32]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b32 = _mm_loaddup_pd(&B[(l_n*36)+32]);
#endif
    __m128d c32_0 = _mm_load_sd(&C[(l_n*36)+1]);
    __m128d a32_0 = _mm_load_sd(&A[92]);
#if defined(__SSE3__) && defined(__AVX__)
    c32_0 = _mm_add_sd(c32_0, _mm_mul_sd(a32_0, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c32_0 = _mm_add_sd(c32_0, _mm_mul_sd(a32_0, b32));
#endif
    _mm_store_sd(&C[(l_n*36)+1], c32_0);
    __m128d c32_1 = _mm_load_sd(&C[(l_n*36)+7]);
    __m128d a32_1 = _mm_load_sd(&A[93]);
#if defined(__SSE3__) && defined(__AVX__)
    c32_1 = _mm_add_sd(c32_1, _mm_mul_sd(a32_1, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c32_1 = _mm_add_sd(c32_1, _mm_mul_sd(a32_1, b32));
#endif
    _mm_store_sd(&C[(l_n*36)+7], c32_1);
    __m128d c32_2 = _mm_load_sd(&C[(l_n*36)+17]);
    __m128d a32_2 = _mm_load_sd(&A[94]);
#if defined(__SSE3__) && defined(__AVX__)
    c32_2 = _mm_add_sd(c32_2, _mm_mul_sd(a32_2, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c32_2 = _mm_add_sd(c32_2, _mm_mul_sd(a32_2, b32));
#endif
    _mm_store_sd(&C[(l_n*36)+17], c32_2);
    __m128d c32_3 = _mm_load_sd(&C[(l_n*36)+32]);
    __m128d a32_3 = _mm_load_sd(&A[95]);
#if defined(__SSE3__) && defined(__AVX__)
    c32_3 = _mm_add_sd(c32_3, _mm_mul_sd(a32_3, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c32_3 = _mm_add_sd(c32_3, _mm_mul_sd(a32_3, b32));
#endif
    _mm_store_sd(&C[(l_n*36)+32], c32_3);
#else
    C[(l_n*36)+1] += A[92] * B[(l_n*36)+32];
    C[(l_n*36)+7] += A[93] * B[(l_n*36)+32];
    C[(l_n*36)+17] += A[94] * B[(l_n*36)+32];
    C[(l_n*36)+32] += A[95] * B[(l_n*36)+32];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b33 = _mm256_broadcast_sd(&B[(l_n*36)+33]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b33 = _mm_loaddup_pd(&B[(l_n*36)+33]);
#endif
    __m128d c33_0 = _mm_load_sd(&C[(l_n*36)+2]);
    __m128d a33_0 = _mm_load_sd(&A[96]);
#if defined(__SSE3__) && defined(__AVX__)
    c33_0 = _mm_add_sd(c33_0, _mm_mul_sd(a33_0, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c33_0 = _mm_add_sd(c33_0, _mm_mul_sd(a33_0, b33));
#endif
    _mm_store_sd(&C[(l_n*36)+2], c33_0);
    __m128d c33_1 = _mm_load_sd(&C[(l_n*36)+8]);
    __m128d a33_1 = _mm_load_sd(&A[97]);
#if defined(__SSE3__) && defined(__AVX__)
    c33_1 = _mm_add_sd(c33_1, _mm_mul_sd(a33_1, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c33_1 = _mm_add_sd(c33_1, _mm_mul_sd(a33_1, b33));
#endif
    _mm_store_sd(&C[(l_n*36)+8], c33_1);
    __m128d c33_2 = _mm_load_sd(&C[(l_n*36)+18]);
    __m128d a33_2 = _mm_load_sd(&A[98]);
#if defined(__SSE3__) && defined(__AVX__)
    c33_2 = _mm_add_sd(c33_2, _mm_mul_sd(a33_2, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c33_2 = _mm_add_sd(c33_2, _mm_mul_sd(a33_2, b33));
#endif
    _mm_store_sd(&C[(l_n*36)+18], c33_2);
    __m128d c33_3 = _mm_load_sd(&C[(l_n*36)+33]);
    __m128d a33_3 = _mm_load_sd(&A[99]);
#if defined(__SSE3__) && defined(__AVX__)
    c33_3 = _mm_add_sd(c33_3, _mm_mul_sd(a33_3, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c33_3 = _mm_add_sd(c33_3, _mm_mul_sd(a33_3, b33));
#endif
    _mm_store_sd(&C[(l_n*36)+33], c33_3);
#else
    C[(l_n*36)+2] += A[96] * B[(l_n*36)+33];
    C[(l_n*36)+8] += A[97] * B[(l_n*36)+33];
    C[(l_n*36)+18] += A[98] * B[(l_n*36)+33];
    C[(l_n*36)+33] += A[99] * B[(l_n*36)+33];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b34 = _mm256_broadcast_sd(&B[(l_n*36)+34]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b34 = _mm_loaddup_pd(&B[(l_n*36)+34]);
#endif
    __m128d c34_0 = _mm_load_sd(&C[(l_n*36)+0]);
    __m128d a34_0 = _mm_load_sd(&A[100]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_0 = _mm_add_sd(c34_0, _mm_mul_sd(a34_0, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_0 = _mm_add_sd(c34_0, _mm_mul_sd(a34_0, b34));
#endif
    _mm_store_sd(&C[(l_n*36)+0], c34_0);
    __m128d c34_1 = _mm_load_sd(&C[(l_n*36)+3]);
    __m128d a34_1 = _mm_load_sd(&A[101]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_1 = _mm_add_sd(c34_1, _mm_mul_sd(a34_1, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_1 = _mm_add_sd(c34_1, _mm_mul_sd(a34_1, b34));
#endif
    _mm_store_sd(&C[(l_n*36)+3], c34_1);
    __m128d c34_2 = _mm_load_sd(&C[(l_n*36)+9]);
    __m128d a34_2 = _mm_load_sd(&A[102]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_2 = _mm_add_sd(c34_2, _mm_mul_sd(a34_2, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_2 = _mm_add_sd(c34_2, _mm_mul_sd(a34_2, b34));
#endif
    _mm_store_sd(&C[(l_n*36)+9], c34_2);
    __m128d c34_3 = _mm_load_sd(&C[(l_n*36)+19]);
    __m128d a34_3 = _mm_load_sd(&A[103]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_3 = _mm_add_sd(c34_3, _mm_mul_sd(a34_3, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_3 = _mm_add_sd(c34_3, _mm_mul_sd(a34_3, b34));
#endif
    _mm_store_sd(&C[(l_n*36)+19], c34_3);
    __m128d c34_4 = _mm_load_sd(&C[(l_n*36)+34]);
    __m128d a34_4 = _mm_load_sd(&A[104]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_4 = _mm_add_sd(c34_4, _mm_mul_sd(a34_4, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_4 = _mm_add_sd(c34_4, _mm_mul_sd(a34_4, b34));
#endif
    _mm_store_sd(&C[(l_n*36)+34], c34_4);
#else
    C[(l_n*36)+0] += A[100] * B[(l_n*36)+34];
    C[(l_n*36)+3] += A[101] * B[(l_n*36)+34];
    C[(l_n*36)+9] += A[102] * B[(l_n*36)+34];
    C[(l_n*36)+19] += A[103] * B[(l_n*36)+34];
    C[(l_n*36)+34] += A[104] * B[(l_n*36)+34];
#endif

  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 1890;
#endif
}

void dsparse_starMatrix_m56_n9_k9_ldA56_ldBna6_ldC56_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_fM1DivM_m56_n9_k56_ldAna6_ldB56_ldC56_beta0_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_n = 0;
  #pragma nounroll_and_jam
  for ( l_n = 0; l_n < 9; l_n++) {
    unsigned int l_m = 0;
   #pragma simd
    for ( l_m = 0; l_m < 56; l_m++) {
      C[(l_n*56)+l_m] = 0.0;
    }
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b0 = _mm256_broadcast_sd(&B[(l_n*56)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b0 = _mm_loaddup_pd(&B[(l_n*56)+0]);
#endif
    __m128d c0_0 = _mm_load_sd(&C[(l_n*56)+0]);
    __m128d a0_0 = _mm_load_sd(&A[0]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, b0));
#endif
    _mm_store_sd(&C[(l_n*56)+0], c0_0);
    __m128d c0_1 = _mm_load_sd(&C[(l_n*56)+3]);
    __m128d a0_1 = _mm_load_sd(&A[1]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, b0));
#endif
    _mm_store_sd(&C[(l_n*56)+3], c0_1);
    __m128d c0_2 = _mm_load_sd(&C[(l_n*56)+9]);
    __m128d a0_2 = _mm_load_sd(&A[2]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, b0));
#endif
    _mm_store_sd(&C[(l_n*56)+9], c0_2);
    __m128d c0_3 = _mm_load_sd(&C[(l_n*56)+19]);
    __m128d a0_3 = _mm_load_sd(&A[3]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_3 = _mm_add_sd(c0_3, _mm_mul_sd(a0_3, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_3 = _mm_add_sd(c0_3, _mm_mul_sd(a0_3, b0));
#endif
    _mm_store_sd(&C[(l_n*56)+19], c0_3);
    __m128d c0_4 = _mm_load_sd(&C[(l_n*56)+34]);
    __m128d a0_4 = _mm_load_sd(&A[4]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_4 = _mm_add_sd(c0_4, _mm_mul_sd(a0_4, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_4 = _mm_add_sd(c0_4, _mm_mul_sd(a0_4, b0));
#endif
    _mm_store_sd(&C[(l_n*56)+34], c0_4);
    __m128d c0_5 = _mm_load_sd(&C[(l_n*56)+55]);
    __m128d a0_5 = _mm_load_sd(&A[5]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_5 = _mm_add_sd(c0_5, _mm_mul_sd(a0_5, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_5 = _mm_add_sd(c0_5, _mm_mul_sd(a0_5, b0));
#endif
    _mm_store_sd(&C[(l_n*56)+55], c0_5);
#else
    C[(l_n*56)+0] += A[0] * B[(l_n*56)+0];
    C[(l_n*56)+3] += A[1] * B[(l_n*56)+0];
    C[(l_n*56)+9] += A[2] * B[(l_n*56)+0];
    C[(l_n*56)+19] += A[3] * B[(l_n*56)+0];
    C[(l_n*56)+34] += A[4] * B[(l_n*56)+0];
    C[(l_n*56)+55] += A[5] * B[(l_n*56)+0];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b1 = _mm256_broadcast_sd(&B[(l_n*56)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b1 = _mm_loaddup_pd(&B[(l_n*56)+1]);
#endif
    __m128d c1_0 = _mm_load_sd(&C[(l_n*56)+1]);
    __m128d a1_0 = _mm_load_sd(&A[6]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, b1));
#endif
    _mm_store_sd(&C[(l_n*56)+1], c1_0);
    __m128d c1_1 = _mm_load_sd(&C[(l_n*56)+7]);
    __m128d a1_1 = _mm_load_sd(&A[7]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_1 = _mm_add_sd(c1_1, _mm_mul_sd(a1_1, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_1 = _mm_add_sd(c1_1, _mm_mul_sd(a1_1, b1));
#endif
    _mm_store_sd(&C[(l_n*56)+7], c1_1);
    __m128d c1_2 = _mm_load_sd(&C[(l_n*56)+17]);
    __m128d a1_2 = _mm_load_sd(&A[8]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_2 = _mm_add_sd(c1_2, _mm_mul_sd(a1_2, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_2 = _mm_add_sd(c1_2, _mm_mul_sd(a1_2, b1));
#endif
    _mm_store_sd(&C[(l_n*56)+17], c1_2);
    __m128d c1_3 = _mm_load_sd(&C[(l_n*56)+32]);
    __m128d a1_3 = _mm_load_sd(&A[9]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_3 = _mm_add_sd(c1_3, _mm_mul_sd(a1_3, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_3 = _mm_add_sd(c1_3, _mm_mul_sd(a1_3, b1));
#endif
    _mm_store_sd(&C[(l_n*56)+32], c1_3);
    __m128d c1_4 = _mm_load_sd(&C[(l_n*56)+53]);
    __m128d a1_4 = _mm_load_sd(&A[10]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_4 = _mm_add_sd(c1_4, _mm_mul_sd(a1_4, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_4 = _mm_add_sd(c1_4, _mm_mul_sd(a1_4, b1));
#endif
    _mm_store_sd(&C[(l_n*56)+53], c1_4);
#else
    C[(l_n*56)+1] += A[6] * B[(l_n*56)+1];
    C[(l_n*56)+7] += A[7] * B[(l_n*56)+1];
    C[(l_n*56)+17] += A[8] * B[(l_n*56)+1];
    C[(l_n*56)+32] += A[9] * B[(l_n*56)+1];
    C[(l_n*56)+53] += A[10] * B[(l_n*56)+1];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b2 = _mm256_broadcast_sd(&B[(l_n*56)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b2 = _mm_loaddup_pd(&B[(l_n*56)+2]);
#endif
    __m128d c2_0 = _mm_load_sd(&C[(l_n*56)+2]);
    __m128d a2_0 = _mm_load_sd(&A[11]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_0 = _mm_add_sd(c2_0, _mm_mul_sd(a2_0, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_0 = _mm_add_sd(c2_0, _mm_mul_sd(a2_0, b2));
#endif
    _mm_store_sd(&C[(l_n*56)+2], c2_0);
    __m128d c2_1 = _mm_load_sd(&C[(l_n*56)+8]);
    __m128d a2_1 = _mm_load_sd(&A[12]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_1 = _mm_add_sd(c2_1, _mm_mul_sd(a2_1, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_1 = _mm_add_sd(c2_1, _mm_mul_sd(a2_1, b2));
#endif
    _mm_store_sd(&C[(l_n*56)+8], c2_1);
    __m128d c2_2 = _mm_load_sd(&C[(l_n*56)+18]);
    __m128d a2_2 = _mm_load_sd(&A[13]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_2 = _mm_add_sd(c2_2, _mm_mul_sd(a2_2, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_2 = _mm_add_sd(c2_2, _mm_mul_sd(a2_2, b2));
#endif
    _mm_store_sd(&C[(l_n*56)+18], c2_2);
    __m128d c2_3 = _mm_load_sd(&C[(l_n*56)+33]);
    __m128d a2_3 = _mm_load_sd(&A[14]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_3 = _mm_add_sd(c2_3, _mm_mul_sd(a2_3, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_3 = _mm_add_sd(c2_3, _mm_mul_sd(a2_3, b2));
#endif
    _mm_store_sd(&C[(l_n*56)+33], c2_3);
    __m128d c2_4 = _mm_load_sd(&C[(l_n*56)+54]);
    __m128d a2_4 = _mm_load_sd(&A[15]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_4 = _mm_add_sd(c2_4, _mm_mul_sd(a2_4, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_4 = _mm_add_sd(c2_4, _mm_mul_sd(a2_4, b2));
#endif
    _mm_store_sd(&C[(l_n*56)+54], c2_4);
#else
    C[(l_n*56)+2] += A[11] * B[(l_n*56)+2];
    C[(l_n*56)+8] += A[12] * B[(l_n*56)+2];
    C[(l_n*56)+18] += A[13] * B[(l_n*56)+2];
    C[(l_n*56)+33] += A[14] * B[(l_n*56)+2];
    C[(l_n*56)+54] += A[15] * B[(l_n*56)+2];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b3 = _mm256_broadcast_sd(&B[(l_n*56)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b3 = _mm_loaddup_pd(&B[(l_n*56)+3]);
#endif
    __m128d c3_0 = _mm_load_sd(&C[(l_n*56)+0]);
    __m128d a3_0 = _mm_load_sd(&A[16]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, b3));
#endif
    _mm_store_sd(&C[(l_n*56)+0], c3_0);
    __m128d c3_1 = _mm_load_sd(&C[(l_n*56)+3]);
    __m128d a3_1 = _mm_load_sd(&A[17]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, b3));
#endif
    _mm_store_sd(&C[(l_n*56)+3], c3_1);
    __m128d c3_2 = _mm_load_sd(&C[(l_n*56)+9]);
    __m128d a3_2 = _mm_load_sd(&A[18]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, b3));
#endif
    _mm_store_sd(&C[(l_n*56)+9], c3_2);
    __m128d c3_3 = _mm_load_sd(&C[(l_n*56)+19]);
    __m128d a3_3 = _mm_load_sd(&A[19]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_3 = _mm_add_sd(c3_3, _mm_mul_sd(a3_3, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_3 = _mm_add_sd(c3_3, _mm_mul_sd(a3_3, b3));
#endif
    _mm_store_sd(&C[(l_n*56)+19], c3_3);
    __m128d c3_4 = _mm_load_sd(&C[(l_n*56)+34]);
    __m128d a3_4 = _mm_load_sd(&A[20]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_4 = _mm_add_sd(c3_4, _mm_mul_sd(a3_4, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_4 = _mm_add_sd(c3_4, _mm_mul_sd(a3_4, b3));
#endif
    _mm_store_sd(&C[(l_n*56)+34], c3_4);
    __m128d c3_5 = _mm_load_sd(&C[(l_n*56)+55]);
    __m128d a3_5 = _mm_load_sd(&A[21]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_5 = _mm_add_sd(c3_5, _mm_mul_sd(a3_5, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_5 = _mm_add_sd(c3_5, _mm_mul_sd(a3_5, b3));
#endif
    _mm_store_sd(&C[(l_n*56)+55], c3_5);
#else
    C[(l_n*56)+0] += A[16] * B[(l_n*56)+3];
    C[(l_n*56)+3] += A[17] * B[(l_n*56)+3];
    C[(l_n*56)+9] += A[18] * B[(l_n*56)+3];
    C[(l_n*56)+19] += A[19] * B[(l_n*56)+3];
    C[(l_n*56)+34] += A[20] * B[(l_n*56)+3];
    C[(l_n*56)+55] += A[21] * B[(l_n*56)+3];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b4 = _mm256_broadcast_sd(&B[(l_n*56)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b4 = _mm_loaddup_pd(&B[(l_n*56)+4]);
#endif
    __m128d c4_0 = _mm_load_sd(&C[(l_n*56)+4]);
    __m128d a4_0 = _mm_load_sd(&A[22]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_0 = _mm_add_sd(c4_0, _mm_mul_sd(a4_0, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_0 = _mm_add_sd(c4_0, _mm_mul_sd(a4_0, b4));
#endif
    _mm_store_sd(&C[(l_n*56)+4], c4_0);
    __m128d c4_1 = _mm_load_sd(&C[(l_n*56)+14]);
    __m128d a4_1 = _mm_load_sd(&A[23]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_1 = _mm_add_sd(c4_1, _mm_mul_sd(a4_1, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_1 = _mm_add_sd(c4_1, _mm_mul_sd(a4_1, b4));
#endif
    _mm_store_sd(&C[(l_n*56)+14], c4_1);
    __m128d c4_2 = _mm_load_sd(&C[(l_n*56)+29]);
    __m128d a4_2 = _mm_load_sd(&A[24]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, b4));
#endif
    _mm_store_sd(&C[(l_n*56)+29], c4_2);
    __m128d c4_3 = _mm_load_sd(&C[(l_n*56)+50]);
    __m128d a4_3 = _mm_load_sd(&A[25]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_3 = _mm_add_sd(c4_3, _mm_mul_sd(a4_3, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_3 = _mm_add_sd(c4_3, _mm_mul_sd(a4_3, b4));
#endif
    _mm_store_sd(&C[(l_n*56)+50], c4_3);
#else
    C[(l_n*56)+4] += A[22] * B[(l_n*56)+4];
    C[(l_n*56)+14] += A[23] * B[(l_n*56)+4];
    C[(l_n*56)+29] += A[24] * B[(l_n*56)+4];
    C[(l_n*56)+50] += A[25] * B[(l_n*56)+4];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b5 = _mm256_broadcast_sd(&B[(l_n*56)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b5 = _mm_loaddup_pd(&B[(l_n*56)+5]);
#endif
    __m128d c5_0 = _mm_load_sd(&C[(l_n*56)+5]);
    __m128d a5_0 = _mm_load_sd(&A[26]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_0 = _mm_add_sd(c5_0, _mm_mul_sd(a5_0, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_0 = _mm_add_sd(c5_0, _mm_mul_sd(a5_0, b5));
#endif
    _mm_store_sd(&C[(l_n*56)+5], c5_0);
    __m128d c5_1 = _mm_load_sd(&C[(l_n*56)+15]);
    __m128d a5_1 = _mm_load_sd(&A[27]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_1 = _mm_add_sd(c5_1, _mm_mul_sd(a5_1, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_1 = _mm_add_sd(c5_1, _mm_mul_sd(a5_1, b5));
#endif
    _mm_store_sd(&C[(l_n*56)+15], c5_1);
    __m128d c5_2 = _mm_load_sd(&C[(l_n*56)+30]);
    __m128d a5_2 = _mm_load_sd(&A[28]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, b5));
#endif
    _mm_store_sd(&C[(l_n*56)+30], c5_2);
    __m128d c5_3 = _mm_load_sd(&C[(l_n*56)+51]);
    __m128d a5_3 = _mm_load_sd(&A[29]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_3 = _mm_add_sd(c5_3, _mm_mul_sd(a5_3, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_3 = _mm_add_sd(c5_3, _mm_mul_sd(a5_3, b5));
#endif
    _mm_store_sd(&C[(l_n*56)+51], c5_3);
#else
    C[(l_n*56)+5] += A[26] * B[(l_n*56)+5];
    C[(l_n*56)+15] += A[27] * B[(l_n*56)+5];
    C[(l_n*56)+30] += A[28] * B[(l_n*56)+5];
    C[(l_n*56)+51] += A[29] * B[(l_n*56)+5];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b6 = _mm256_broadcast_sd(&B[(l_n*56)+6]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b6 = _mm_loaddup_pd(&B[(l_n*56)+6]);
#endif
    __m128d c6_0 = _mm_load_sd(&C[(l_n*56)+6]);
    __m128d a6_0 = _mm_load_sd(&A[30]);
#if defined(__SSE3__) && defined(__AVX__)
    c6_0 = _mm_add_sd(c6_0, _mm_mul_sd(a6_0, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c6_0 = _mm_add_sd(c6_0, _mm_mul_sd(a6_0, b6));
#endif
    _mm_store_sd(&C[(l_n*56)+6], c6_0);
    __m128d c6_1 = _mm_load_sd(&C[(l_n*56)+16]);
    __m128d a6_1 = _mm_load_sd(&A[31]);
#if defined(__SSE3__) && defined(__AVX__)
    c6_1 = _mm_add_sd(c6_1, _mm_mul_sd(a6_1, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c6_1 = _mm_add_sd(c6_1, _mm_mul_sd(a6_1, b6));
#endif
    _mm_store_sd(&C[(l_n*56)+16], c6_1);
    __m128d c6_2 = _mm_load_sd(&C[(l_n*56)+31]);
    __m128d a6_2 = _mm_load_sd(&A[32]);
#if defined(__SSE3__) && defined(__AVX__)
    c6_2 = _mm_add_sd(c6_2, _mm_mul_sd(a6_2, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c6_2 = _mm_add_sd(c6_2, _mm_mul_sd(a6_2, b6));
#endif
    _mm_store_sd(&C[(l_n*56)+31], c6_2);
    __m128d c6_3 = _mm_load_sd(&C[(l_n*56)+52]);
    __m128d a6_3 = _mm_load_sd(&A[33]);
#if defined(__SSE3__) && defined(__AVX__)
    c6_3 = _mm_add_sd(c6_3, _mm_mul_sd(a6_3, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c6_3 = _mm_add_sd(c6_3, _mm_mul_sd(a6_3, b6));
#endif
    _mm_store_sd(&C[(l_n*56)+52], c6_3);
#else
    C[(l_n*56)+6] += A[30] * B[(l_n*56)+6];
    C[(l_n*56)+16] += A[31] * B[(l_n*56)+6];
    C[(l_n*56)+31] += A[32] * B[(l_n*56)+6];
    C[(l_n*56)+52] += A[33] * B[(l_n*56)+6];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b7 = _mm256_broadcast_sd(&B[(l_n*56)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b7 = _mm_loaddup_pd(&B[(l_n*56)+7]);
#endif
    __m128d c7_0 = _mm_load_sd(&C[(l_n*56)+1]);
    __m128d a7_0 = _mm_load_sd(&A[34]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_0 = _mm_add_sd(c7_0, _mm_mul_sd(a7_0, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_0 = _mm_add_sd(c7_0, _mm_mul_sd(a7_0, b7));
#endif
    _mm_store_sd(&C[(l_n*56)+1], c7_0);
    __m128d c7_1 = _mm_load_sd(&C[(l_n*56)+7]);
    __m128d a7_1 = _mm_load_sd(&A[35]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_1 = _mm_add_sd(c7_1, _mm_mul_sd(a7_1, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_1 = _mm_add_sd(c7_1, _mm_mul_sd(a7_1, b7));
#endif
    _mm_store_sd(&C[(l_n*56)+7], c7_1);
    __m128d c7_2 = _mm_load_sd(&C[(l_n*56)+17]);
    __m128d a7_2 = _mm_load_sd(&A[36]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_2 = _mm_add_sd(c7_2, _mm_mul_sd(a7_2, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_2 = _mm_add_sd(c7_2, _mm_mul_sd(a7_2, b7));
#endif
    _mm_store_sd(&C[(l_n*56)+17], c7_2);
    __m128d c7_3 = _mm_load_sd(&C[(l_n*56)+32]);
    __m128d a7_3 = _mm_load_sd(&A[37]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_3 = _mm_add_sd(c7_3, _mm_mul_sd(a7_3, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_3 = _mm_add_sd(c7_3, _mm_mul_sd(a7_3, b7));
#endif
    _mm_store_sd(&C[(l_n*56)+32], c7_3);
    __m128d c7_4 = _mm_load_sd(&C[(l_n*56)+53]);
    __m128d a7_4 = _mm_load_sd(&A[38]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_4 = _mm_add_sd(c7_4, _mm_mul_sd(a7_4, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_4 = _mm_add_sd(c7_4, _mm_mul_sd(a7_4, b7));
#endif
    _mm_store_sd(&C[(l_n*56)+53], c7_4);
#else
    C[(l_n*56)+1] += A[34] * B[(l_n*56)+7];
    C[(l_n*56)+7] += A[35] * B[(l_n*56)+7];
    C[(l_n*56)+17] += A[36] * B[(l_n*56)+7];
    C[(l_n*56)+32] += A[37] * B[(l_n*56)+7];
    C[(l_n*56)+53] += A[38] * B[(l_n*56)+7];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b8 = _mm256_broadcast_sd(&B[(l_n*56)+8]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b8 = _mm_loaddup_pd(&B[(l_n*56)+8]);
#endif
    __m128d c8_0 = _mm_load_sd(&C[(l_n*56)+2]);
    __m128d a8_0 = _mm_load_sd(&A[39]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_0 = _mm_add_sd(c8_0, _mm_mul_sd(a8_0, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_0 = _mm_add_sd(c8_0, _mm_mul_sd(a8_0, b8));
#endif
    _mm_store_sd(&C[(l_n*56)+2], c8_0);
    __m128d c8_1 = _mm_load_sd(&C[(l_n*56)+8]);
    __m128d a8_1 = _mm_load_sd(&A[40]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_1 = _mm_add_sd(c8_1, _mm_mul_sd(a8_1, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_1 = _mm_add_sd(c8_1, _mm_mul_sd(a8_1, b8));
#endif
    _mm_store_sd(&C[(l_n*56)+8], c8_1);
    __m128d c8_2 = _mm_load_sd(&C[(l_n*56)+18]);
    __m128d a8_2 = _mm_load_sd(&A[41]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_2 = _mm_add_sd(c8_2, _mm_mul_sd(a8_2, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_2 = _mm_add_sd(c8_2, _mm_mul_sd(a8_2, b8));
#endif
    _mm_store_sd(&C[(l_n*56)+18], c8_2);
    __m128d c8_3 = _mm_load_sd(&C[(l_n*56)+33]);
    __m128d a8_3 = _mm_load_sd(&A[42]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_3 = _mm_add_sd(c8_3, _mm_mul_sd(a8_3, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_3 = _mm_add_sd(c8_3, _mm_mul_sd(a8_3, b8));
#endif
    _mm_store_sd(&C[(l_n*56)+33], c8_3);
    __m128d c8_4 = _mm_load_sd(&C[(l_n*56)+54]);
    __m128d a8_4 = _mm_load_sd(&A[43]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_4 = _mm_add_sd(c8_4, _mm_mul_sd(a8_4, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_4 = _mm_add_sd(c8_4, _mm_mul_sd(a8_4, b8));
#endif
    _mm_store_sd(&C[(l_n*56)+54], c8_4);
#else
    C[(l_n*56)+2] += A[39] * B[(l_n*56)+8];
    C[(l_n*56)+8] += A[40] * B[(l_n*56)+8];
    C[(l_n*56)+18] += A[41] * B[(l_n*56)+8];
    C[(l_n*56)+33] += A[42] * B[(l_n*56)+8];
    C[(l_n*56)+54] += A[43] * B[(l_n*56)+8];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b9 = _mm256_broadcast_sd(&B[(l_n*56)+9]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b9 = _mm_loaddup_pd(&B[(l_n*56)+9]);
#endif
    __m128d c9_0 = _mm_load_sd(&C[(l_n*56)+0]);
    __m128d a9_0 = _mm_load_sd(&A[44]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, b9));
#endif
    _mm_store_sd(&C[(l_n*56)+0], c9_0);
    __m128d c9_1 = _mm_load_sd(&C[(l_n*56)+3]);
    __m128d a9_1 = _mm_load_sd(&A[45]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, b9));
#endif
    _mm_store_sd(&C[(l_n*56)+3], c9_1);
    __m128d c9_2 = _mm_load_sd(&C[(l_n*56)+9]);
    __m128d a9_2 = _mm_load_sd(&A[46]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, b9));
#endif
    _mm_store_sd(&C[(l_n*56)+9], c9_2);
    __m128d c9_3 = _mm_load_sd(&C[(l_n*56)+19]);
    __m128d a9_3 = _mm_load_sd(&A[47]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_3 = _mm_add_sd(c9_3, _mm_mul_sd(a9_3, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_3 = _mm_add_sd(c9_3, _mm_mul_sd(a9_3, b9));
#endif
    _mm_store_sd(&C[(l_n*56)+19], c9_3);
    __m128d c9_4 = _mm_load_sd(&C[(l_n*56)+34]);
    __m128d a9_4 = _mm_load_sd(&A[48]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_4 = _mm_add_sd(c9_4, _mm_mul_sd(a9_4, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_4 = _mm_add_sd(c9_4, _mm_mul_sd(a9_4, b9));
#endif
    _mm_store_sd(&C[(l_n*56)+34], c9_4);
    __m128d c9_5 = _mm_load_sd(&C[(l_n*56)+55]);
    __m128d a9_5 = _mm_load_sd(&A[49]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_5 = _mm_add_sd(c9_5, _mm_mul_sd(a9_5, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_5 = _mm_add_sd(c9_5, _mm_mul_sd(a9_5, b9));
#endif
    _mm_store_sd(&C[(l_n*56)+55], c9_5);
#else
    C[(l_n*56)+0] += A[44] * B[(l_n*56)+9];
    C[(l_n*56)+3] += A[45] * B[(l_n*56)+9];
    C[(l_n*56)+9] += A[46] * B[(l_n*56)+9];
    C[(l_n*56)+19] += A[47] * B[(l_n*56)+9];
    C[(l_n*56)+34] += A[48] * B[(l_n*56)+9];
    C[(l_n*56)+55] += A[49] * B[(l_n*56)+9];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b10 = _mm256_broadcast_sd(&B[(l_n*56)+10]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b10 = _mm_loaddup_pd(&B[(l_n*56)+10]);
#endif
    __m128d c10_0 = _mm_load_sd(&C[(l_n*56)+10]);
    __m128d a10_0 = _mm_load_sd(&A[50]);
#if defined(__SSE3__) && defined(__AVX__)
    c10_0 = _mm_add_sd(c10_0, _mm_mul_sd(a10_0, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c10_0 = _mm_add_sd(c10_0, _mm_mul_sd(a10_0, b10));
#endif
    _mm_store_sd(&C[(l_n*56)+10], c10_0);
    __m128d c10_1 = _mm_load_sd(&C[(l_n*56)+25]);
    __m128d a10_1 = _mm_load_sd(&A[51]);
#if defined(__SSE3__) && defined(__AVX__)
    c10_1 = _mm_add_sd(c10_1, _mm_mul_sd(a10_1, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c10_1 = _mm_add_sd(c10_1, _mm_mul_sd(a10_1, b10));
#endif
    _mm_store_sd(&C[(l_n*56)+25], c10_1);
    __m128d c10_2 = _mm_load_sd(&C[(l_n*56)+46]);
    __m128d a10_2 = _mm_load_sd(&A[52]);
#if defined(__SSE3__) && defined(__AVX__)
    c10_2 = _mm_add_sd(c10_2, _mm_mul_sd(a10_2, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c10_2 = _mm_add_sd(c10_2, _mm_mul_sd(a10_2, b10));
#endif
    _mm_store_sd(&C[(l_n*56)+46], c10_2);
#else
    C[(l_n*56)+10] += A[50] * B[(l_n*56)+10];
    C[(l_n*56)+25] += A[51] * B[(l_n*56)+10];
    C[(l_n*56)+46] += A[52] * B[(l_n*56)+10];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b11 = _mm256_broadcast_sd(&B[(l_n*56)+11]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b11 = _mm_loaddup_pd(&B[(l_n*56)+11]);
#endif
    __m128d c11_0 = _mm_load_sd(&C[(l_n*56)+11]);
    __m128d a11_0 = _mm_load_sd(&A[53]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_0 = _mm_add_sd(c11_0, _mm_mul_sd(a11_0, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_0 = _mm_add_sd(c11_0, _mm_mul_sd(a11_0, b11));
#endif
    _mm_store_sd(&C[(l_n*56)+11], c11_0);
    __m128d c11_1 = _mm_load_sd(&C[(l_n*56)+26]);
    __m128d a11_1 = _mm_load_sd(&A[54]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_1 = _mm_add_sd(c11_1, _mm_mul_sd(a11_1, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_1 = _mm_add_sd(c11_1, _mm_mul_sd(a11_1, b11));
#endif
    _mm_store_sd(&C[(l_n*56)+26], c11_1);
    __m128d c11_2 = _mm_load_sd(&C[(l_n*56)+47]);
    __m128d a11_2 = _mm_load_sd(&A[55]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_2 = _mm_add_sd(c11_2, _mm_mul_sd(a11_2, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_2 = _mm_add_sd(c11_2, _mm_mul_sd(a11_2, b11));
#endif
    _mm_store_sd(&C[(l_n*56)+47], c11_2);
#else
    C[(l_n*56)+11] += A[53] * B[(l_n*56)+11];
    C[(l_n*56)+26] += A[54] * B[(l_n*56)+11];
    C[(l_n*56)+47] += A[55] * B[(l_n*56)+11];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b12 = _mm256_broadcast_sd(&B[(l_n*56)+12]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b12 = _mm_loaddup_pd(&B[(l_n*56)+12]);
#endif
    __m128d c12_0 = _mm_load_sd(&C[(l_n*56)+12]);
    __m128d a12_0 = _mm_load_sd(&A[56]);
#if defined(__SSE3__) && defined(__AVX__)
    c12_0 = _mm_add_sd(c12_0, _mm_mul_sd(a12_0, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c12_0 = _mm_add_sd(c12_0, _mm_mul_sd(a12_0, b12));
#endif
    _mm_store_sd(&C[(l_n*56)+12], c12_0);
    __m128d c12_1 = _mm_load_sd(&C[(l_n*56)+27]);
    __m128d a12_1 = _mm_load_sd(&A[57]);
#if defined(__SSE3__) && defined(__AVX__)
    c12_1 = _mm_add_sd(c12_1, _mm_mul_sd(a12_1, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c12_1 = _mm_add_sd(c12_1, _mm_mul_sd(a12_1, b12));
#endif
    _mm_store_sd(&C[(l_n*56)+27], c12_1);
    __m128d c12_2 = _mm_load_sd(&C[(l_n*56)+48]);
    __m128d a12_2 = _mm_load_sd(&A[58]);
#if defined(__SSE3__) && defined(__AVX__)
    c12_2 = _mm_add_sd(c12_2, _mm_mul_sd(a12_2, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c12_2 = _mm_add_sd(c12_2, _mm_mul_sd(a12_2, b12));
#endif
    _mm_store_sd(&C[(l_n*56)+48], c12_2);
#else
    C[(l_n*56)+12] += A[56] * B[(l_n*56)+12];
    C[(l_n*56)+27] += A[57] * B[(l_n*56)+12];
    C[(l_n*56)+48] += A[58] * B[(l_n*56)+12];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b13 = _mm256_broadcast_sd(&B[(l_n*56)+13]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b13 = _mm_loaddup_pd(&B[(l_n*56)+13]);
#endif
    __m128d c13_0 = _mm_load_sd(&C[(l_n*56)+13]);
    __m128d a13_0 = _mm_load_sd(&A[59]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_0 = _mm_add_sd(c13_0, _mm_mul_sd(a13_0, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_0 = _mm_add_sd(c13_0, _mm_mul_sd(a13_0, b13));
#endif
    _mm_store_sd(&C[(l_n*56)+13], c13_0);
    __m128d c13_1 = _mm_load_sd(&C[(l_n*56)+28]);
    __m128d a13_1 = _mm_load_sd(&A[60]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_1 = _mm_add_sd(c13_1, _mm_mul_sd(a13_1, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_1 = _mm_add_sd(c13_1, _mm_mul_sd(a13_1, b13));
#endif
    _mm_store_sd(&C[(l_n*56)+28], c13_1);
    __m128d c13_2 = _mm_load_sd(&C[(l_n*56)+49]);
    __m128d a13_2 = _mm_load_sd(&A[61]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_2 = _mm_add_sd(c13_2, _mm_mul_sd(a13_2, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_2 = _mm_add_sd(c13_2, _mm_mul_sd(a13_2, b13));
#endif
    _mm_store_sd(&C[(l_n*56)+49], c13_2);
#else
    C[(l_n*56)+13] += A[59] * B[(l_n*56)+13];
    C[(l_n*56)+28] += A[60] * B[(l_n*56)+13];
    C[(l_n*56)+49] += A[61] * B[(l_n*56)+13];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b14 = _mm256_broadcast_sd(&B[(l_n*56)+14]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b14 = _mm_loaddup_pd(&B[(l_n*56)+14]);
#endif
    __m128d c14_0 = _mm_load_sd(&C[(l_n*56)+4]);
    __m128d a14_0 = _mm_load_sd(&A[62]);
#if defined(__SSE3__) && defined(__AVX__)
    c14_0 = _mm_add_sd(c14_0, _mm_mul_sd(a14_0, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c14_0 = _mm_add_sd(c14_0, _mm_mul_sd(a14_0, b14));
#endif
    _mm_store_sd(&C[(l_n*56)+4], c14_0);
    __m128d c14_1 = _mm_load_sd(&C[(l_n*56)+14]);
    __m128d a14_1 = _mm_load_sd(&A[63]);
#if defined(__SSE3__) && defined(__AVX__)
    c14_1 = _mm_add_sd(c14_1, _mm_mul_sd(a14_1, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c14_1 = _mm_add_sd(c14_1, _mm_mul_sd(a14_1, b14));
#endif
    _mm_store_sd(&C[(l_n*56)+14], c14_1);
    __m128d c14_2 = _mm_load_sd(&C[(l_n*56)+29]);
    __m128d a14_2 = _mm_load_sd(&A[64]);
#if defined(__SSE3__) && defined(__AVX__)
    c14_2 = _mm_add_sd(c14_2, _mm_mul_sd(a14_2, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c14_2 = _mm_add_sd(c14_2, _mm_mul_sd(a14_2, b14));
#endif
    _mm_store_sd(&C[(l_n*56)+29], c14_2);
    __m128d c14_3 = _mm_load_sd(&C[(l_n*56)+50]);
    __m128d a14_3 = _mm_load_sd(&A[65]);
#if defined(__SSE3__) && defined(__AVX__)
    c14_3 = _mm_add_sd(c14_3, _mm_mul_sd(a14_3, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c14_3 = _mm_add_sd(c14_3, _mm_mul_sd(a14_3, b14));
#endif
    _mm_store_sd(&C[(l_n*56)+50], c14_3);
#else
    C[(l_n*56)+4] += A[62] * B[(l_n*56)+14];
    C[(l_n*56)+14] += A[63] * B[(l_n*56)+14];
    C[(l_n*56)+29] += A[64] * B[(l_n*56)+14];
    C[(l_n*56)+50] += A[65] * B[(l_n*56)+14];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b15 = _mm256_broadcast_sd(&B[(l_n*56)+15]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b15 = _mm_loaddup_pd(&B[(l_n*56)+15]);
#endif
    __m128d c15_0 = _mm_load_sd(&C[(l_n*56)+5]);
    __m128d a15_0 = _mm_load_sd(&A[66]);
#if defined(__SSE3__) && defined(__AVX__)
    c15_0 = _mm_add_sd(c15_0, _mm_mul_sd(a15_0, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c15_0 = _mm_add_sd(c15_0, _mm_mul_sd(a15_0, b15));
#endif
    _mm_store_sd(&C[(l_n*56)+5], c15_0);
    __m128d c15_1 = _mm_load_sd(&C[(l_n*56)+15]);
    __m128d a15_1 = _mm_load_sd(&A[67]);
#if defined(__SSE3__) && defined(__AVX__)
    c15_1 = _mm_add_sd(c15_1, _mm_mul_sd(a15_1, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c15_1 = _mm_add_sd(c15_1, _mm_mul_sd(a15_1, b15));
#endif
    _mm_store_sd(&C[(l_n*56)+15], c15_1);
    __m128d c15_2 = _mm_load_sd(&C[(l_n*56)+30]);
    __m128d a15_2 = _mm_load_sd(&A[68]);
#if defined(__SSE3__) && defined(__AVX__)
    c15_2 = _mm_add_sd(c15_2, _mm_mul_sd(a15_2, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c15_2 = _mm_add_sd(c15_2, _mm_mul_sd(a15_2, b15));
#endif
    _mm_store_sd(&C[(l_n*56)+30], c15_2);
    __m128d c15_3 = _mm_load_sd(&C[(l_n*56)+51]);
    __m128d a15_3 = _mm_load_sd(&A[69]);
#if defined(__SSE3__) && defined(__AVX__)
    c15_3 = _mm_add_sd(c15_3, _mm_mul_sd(a15_3, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c15_3 = _mm_add_sd(c15_3, _mm_mul_sd(a15_3, b15));
#endif
    _mm_store_sd(&C[(l_n*56)+51], c15_3);
#else
    C[(l_n*56)+5] += A[66] * B[(l_n*56)+15];
    C[(l_n*56)+15] += A[67] * B[(l_n*56)+15];
    C[(l_n*56)+30] += A[68] * B[(l_n*56)+15];
    C[(l_n*56)+51] += A[69] * B[(l_n*56)+15];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b16 = _mm256_broadcast_sd(&B[(l_n*56)+16]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b16 = _mm_loaddup_pd(&B[(l_n*56)+16]);
#endif
    __m128d c16_0 = _mm_load_sd(&C[(l_n*56)+6]);
    __m128d a16_0 = _mm_load_sd(&A[70]);
#if defined(__SSE3__) && defined(__AVX__)
    c16_0 = _mm_add_sd(c16_0, _mm_mul_sd(a16_0, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c16_0 = _mm_add_sd(c16_0, _mm_mul_sd(a16_0, b16));
#endif
    _mm_store_sd(&C[(l_n*56)+6], c16_0);
    __m128d c16_1 = _mm_load_sd(&C[(l_n*56)+16]);
    __m128d a16_1 = _mm_load_sd(&A[71]);
#if defined(__SSE3__) && defined(__AVX__)
    c16_1 = _mm_add_sd(c16_1, _mm_mul_sd(a16_1, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c16_1 = _mm_add_sd(c16_1, _mm_mul_sd(a16_1, b16));
#endif
    _mm_store_sd(&C[(l_n*56)+16], c16_1);
    __m128d c16_2 = _mm_load_sd(&C[(l_n*56)+31]);
    __m128d a16_2 = _mm_load_sd(&A[72]);
#if defined(__SSE3__) && defined(__AVX__)
    c16_2 = _mm_add_sd(c16_2, _mm_mul_sd(a16_2, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c16_2 = _mm_add_sd(c16_2, _mm_mul_sd(a16_2, b16));
#endif
    _mm_store_sd(&C[(l_n*56)+31], c16_2);
    __m128d c16_3 = _mm_load_sd(&C[(l_n*56)+52]);
    __m128d a16_3 = _mm_load_sd(&A[73]);
#if defined(__SSE3__) && defined(__AVX__)
    c16_3 = _mm_add_sd(c16_3, _mm_mul_sd(a16_3, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c16_3 = _mm_add_sd(c16_3, _mm_mul_sd(a16_3, b16));
#endif
    _mm_store_sd(&C[(l_n*56)+52], c16_3);
#else
    C[(l_n*56)+6] += A[70] * B[(l_n*56)+16];
    C[(l_n*56)+16] += A[71] * B[(l_n*56)+16];
    C[(l_n*56)+31] += A[72] * B[(l_n*56)+16];
    C[(l_n*56)+52] += A[73] * B[(l_n*56)+16];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b17 = _mm256_broadcast_sd(&B[(l_n*56)+17]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b17 = _mm_loaddup_pd(&B[(l_n*56)+17]);
#endif
    __m128d c17_0 = _mm_load_sd(&C[(l_n*56)+1]);
    __m128d a17_0 = _mm_load_sd(&A[74]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_0 = _mm_add_sd(c17_0, _mm_mul_sd(a17_0, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_0 = _mm_add_sd(c17_0, _mm_mul_sd(a17_0, b17));
#endif
    _mm_store_sd(&C[(l_n*56)+1], c17_0);
    __m128d c17_1 = _mm_load_sd(&C[(l_n*56)+7]);
    __m128d a17_1 = _mm_load_sd(&A[75]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_1 = _mm_add_sd(c17_1, _mm_mul_sd(a17_1, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_1 = _mm_add_sd(c17_1, _mm_mul_sd(a17_1, b17));
#endif
    _mm_store_sd(&C[(l_n*56)+7], c17_1);
    __m128d c17_2 = _mm_load_sd(&C[(l_n*56)+17]);
    __m128d a17_2 = _mm_load_sd(&A[76]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_2 = _mm_add_sd(c17_2, _mm_mul_sd(a17_2, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_2 = _mm_add_sd(c17_2, _mm_mul_sd(a17_2, b17));
#endif
    _mm_store_sd(&C[(l_n*56)+17], c17_2);
    __m128d c17_3 = _mm_load_sd(&C[(l_n*56)+32]);
    __m128d a17_3 = _mm_load_sd(&A[77]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_3 = _mm_add_sd(c17_3, _mm_mul_sd(a17_3, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_3 = _mm_add_sd(c17_3, _mm_mul_sd(a17_3, b17));
#endif
    _mm_store_sd(&C[(l_n*56)+32], c17_3);
    __m128d c17_4 = _mm_load_sd(&C[(l_n*56)+53]);
    __m128d a17_4 = _mm_load_sd(&A[78]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_4 = _mm_add_sd(c17_4, _mm_mul_sd(a17_4, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_4 = _mm_add_sd(c17_4, _mm_mul_sd(a17_4, b17));
#endif
    _mm_store_sd(&C[(l_n*56)+53], c17_4);
#else
    C[(l_n*56)+1] += A[74] * B[(l_n*56)+17];
    C[(l_n*56)+7] += A[75] * B[(l_n*56)+17];
    C[(l_n*56)+17] += A[76] * B[(l_n*56)+17];
    C[(l_n*56)+32] += A[77] * B[(l_n*56)+17];
    C[(l_n*56)+53] += A[78] * B[(l_n*56)+17];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b18 = _mm256_broadcast_sd(&B[(l_n*56)+18]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b18 = _mm_loaddup_pd(&B[(l_n*56)+18]);
#endif
    __m128d c18_0 = _mm_load_sd(&C[(l_n*56)+2]);
    __m128d a18_0 = _mm_load_sd(&A[79]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_0 = _mm_add_sd(c18_0, _mm_mul_sd(a18_0, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_0 = _mm_add_sd(c18_0, _mm_mul_sd(a18_0, b18));
#endif
    _mm_store_sd(&C[(l_n*56)+2], c18_0);
    __m128d c18_1 = _mm_load_sd(&C[(l_n*56)+8]);
    __m128d a18_1 = _mm_load_sd(&A[80]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_1 = _mm_add_sd(c18_1, _mm_mul_sd(a18_1, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_1 = _mm_add_sd(c18_1, _mm_mul_sd(a18_1, b18));
#endif
    _mm_store_sd(&C[(l_n*56)+8], c18_1);
    __m128d c18_2 = _mm_load_sd(&C[(l_n*56)+18]);
    __m128d a18_2 = _mm_load_sd(&A[81]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_2 = _mm_add_sd(c18_2, _mm_mul_sd(a18_2, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_2 = _mm_add_sd(c18_2, _mm_mul_sd(a18_2, b18));
#endif
    _mm_store_sd(&C[(l_n*56)+18], c18_2);
    __m128d c18_3 = _mm_load_sd(&C[(l_n*56)+33]);
    __m128d a18_3 = _mm_load_sd(&A[82]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_3 = _mm_add_sd(c18_3, _mm_mul_sd(a18_3, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_3 = _mm_add_sd(c18_3, _mm_mul_sd(a18_3, b18));
#endif
    _mm_store_sd(&C[(l_n*56)+33], c18_3);
    __m128d c18_4 = _mm_load_sd(&C[(l_n*56)+54]);
    __m128d a18_4 = _mm_load_sd(&A[83]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_4 = _mm_add_sd(c18_4, _mm_mul_sd(a18_4, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_4 = _mm_add_sd(c18_4, _mm_mul_sd(a18_4, b18));
#endif
    _mm_store_sd(&C[(l_n*56)+54], c18_4);
#else
    C[(l_n*56)+2] += A[79] * B[(l_n*56)+18];
    C[(l_n*56)+8] += A[80] * B[(l_n*56)+18];
    C[(l_n*56)+18] += A[81] * B[(l_n*56)+18];
    C[(l_n*56)+33] += A[82] * B[(l_n*56)+18];
    C[(l_n*56)+54] += A[83] * B[(l_n*56)+18];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b19 = _mm256_broadcast_sd(&B[(l_n*56)+19]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b19 = _mm_loaddup_pd(&B[(l_n*56)+19]);
#endif
    __m128d c19_0 = _mm_load_sd(&C[(l_n*56)+0]);
    __m128d a19_0 = _mm_load_sd(&A[84]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_0 = _mm_add_sd(c19_0, _mm_mul_sd(a19_0, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_0 = _mm_add_sd(c19_0, _mm_mul_sd(a19_0, b19));
#endif
    _mm_store_sd(&C[(l_n*56)+0], c19_0);
    __m128d c19_1 = _mm_load_sd(&C[(l_n*56)+3]);
    __m128d a19_1 = _mm_load_sd(&A[85]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_1 = _mm_add_sd(c19_1, _mm_mul_sd(a19_1, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_1 = _mm_add_sd(c19_1, _mm_mul_sd(a19_1, b19));
#endif
    _mm_store_sd(&C[(l_n*56)+3], c19_1);
    __m128d c19_2 = _mm_load_sd(&C[(l_n*56)+9]);
    __m128d a19_2 = _mm_load_sd(&A[86]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_2 = _mm_add_sd(c19_2, _mm_mul_sd(a19_2, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_2 = _mm_add_sd(c19_2, _mm_mul_sd(a19_2, b19));
#endif
    _mm_store_sd(&C[(l_n*56)+9], c19_2);
    __m128d c19_3 = _mm_load_sd(&C[(l_n*56)+19]);
    __m128d a19_3 = _mm_load_sd(&A[87]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_3 = _mm_add_sd(c19_3, _mm_mul_sd(a19_3, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_3 = _mm_add_sd(c19_3, _mm_mul_sd(a19_3, b19));
#endif
    _mm_store_sd(&C[(l_n*56)+19], c19_3);
    __m128d c19_4 = _mm_load_sd(&C[(l_n*56)+34]);
    __m128d a19_4 = _mm_load_sd(&A[88]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_4 = _mm_add_sd(c19_4, _mm_mul_sd(a19_4, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_4 = _mm_add_sd(c19_4, _mm_mul_sd(a19_4, b19));
#endif
    _mm_store_sd(&C[(l_n*56)+34], c19_4);
    __m128d c19_5 = _mm_load_sd(&C[(l_n*56)+55]);
    __m128d a19_5 = _mm_load_sd(&A[89]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_5 = _mm_add_sd(c19_5, _mm_mul_sd(a19_5, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_5 = _mm_add_sd(c19_5, _mm_mul_sd(a19_5, b19));
#endif
    _mm_store_sd(&C[(l_n*56)+55], c19_5);
#else
    C[(l_n*56)+0] += A[84] * B[(l_n*56)+19];
    C[(l_n*56)+3] += A[85] * B[(l_n*56)+19];
    C[(l_n*56)+9] += A[86] * B[(l_n*56)+19];
    C[(l_n*56)+19] += A[87] * B[(l_n*56)+19];
    C[(l_n*56)+34] += A[88] * B[(l_n*56)+19];
    C[(l_n*56)+55] += A[89] * B[(l_n*56)+19];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b20 = _mm256_broadcast_sd(&B[(l_n*56)+20]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b20 = _mm_loaddup_pd(&B[(l_n*56)+20]);
#endif
    __m128d c20_0 = _mm_load_sd(&C[(l_n*56)+20]);
    __m128d a20_0 = _mm_load_sd(&A[90]);
#if defined(__SSE3__) && defined(__AVX__)
    c20_0 = _mm_add_sd(c20_0, _mm_mul_sd(a20_0, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c20_0 = _mm_add_sd(c20_0, _mm_mul_sd(a20_0, b20));
#endif
    _mm_store_sd(&C[(l_n*56)+20], c20_0);
    __m128d c20_1 = _mm_load_sd(&C[(l_n*56)+41]);
    __m128d a20_1 = _mm_load_sd(&A[91]);
#if defined(__SSE3__) && defined(__AVX__)
    c20_1 = _mm_add_sd(c20_1, _mm_mul_sd(a20_1, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c20_1 = _mm_add_sd(c20_1, _mm_mul_sd(a20_1, b20));
#endif
    _mm_store_sd(&C[(l_n*56)+41], c20_1);
#else
    C[(l_n*56)+20] += A[90] * B[(l_n*56)+20];
    C[(l_n*56)+41] += A[91] * B[(l_n*56)+20];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b21 = _mm256_broadcast_sd(&B[(l_n*56)+21]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b21 = _mm_loaddup_pd(&B[(l_n*56)+21]);
#endif
    __m128d c21_0 = _mm_load_sd(&C[(l_n*56)+21]);
    __m128d a21_0 = _mm_load_sd(&A[92]);
#if defined(__SSE3__) && defined(__AVX__)
    c21_0 = _mm_add_sd(c21_0, _mm_mul_sd(a21_0, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c21_0 = _mm_add_sd(c21_0, _mm_mul_sd(a21_0, b21));
#endif
    _mm_store_sd(&C[(l_n*56)+21], c21_0);
    __m128d c21_1 = _mm_load_sd(&C[(l_n*56)+42]);
    __m128d a21_1 = _mm_load_sd(&A[93]);
#if defined(__SSE3__) && defined(__AVX__)
    c21_1 = _mm_add_sd(c21_1, _mm_mul_sd(a21_1, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c21_1 = _mm_add_sd(c21_1, _mm_mul_sd(a21_1, b21));
#endif
    _mm_store_sd(&C[(l_n*56)+42], c21_1);
#else
    C[(l_n*56)+21] += A[92] * B[(l_n*56)+21];
    C[(l_n*56)+42] += A[93] * B[(l_n*56)+21];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b22 = _mm256_broadcast_sd(&B[(l_n*56)+22]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b22 = _mm_loaddup_pd(&B[(l_n*56)+22]);
#endif
    __m128d c22_0 = _mm_load_sd(&C[(l_n*56)+22]);
    __m128d a22_0 = _mm_load_sd(&A[94]);
#if defined(__SSE3__) && defined(__AVX__)
    c22_0 = _mm_add_sd(c22_0, _mm_mul_sd(a22_0, _mm256_castpd256_pd128(b22)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c22_0 = _mm_add_sd(c22_0, _mm_mul_sd(a22_0, b22));
#endif
    _mm_store_sd(&C[(l_n*56)+22], c22_0);
    __m128d c22_1 = _mm_load_sd(&C[(l_n*56)+43]);
    __m128d a22_1 = _mm_load_sd(&A[95]);
#if defined(__SSE3__) && defined(__AVX__)
    c22_1 = _mm_add_sd(c22_1, _mm_mul_sd(a22_1, _mm256_castpd256_pd128(b22)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c22_1 = _mm_add_sd(c22_1, _mm_mul_sd(a22_1, b22));
#endif
    _mm_store_sd(&C[(l_n*56)+43], c22_1);
#else
    C[(l_n*56)+22] += A[94] * B[(l_n*56)+22];
    C[(l_n*56)+43] += A[95] * B[(l_n*56)+22];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b23 = _mm256_broadcast_sd(&B[(l_n*56)+23]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b23 = _mm_loaddup_pd(&B[(l_n*56)+23]);
#endif
    __m128d c23_0 = _mm_load_sd(&C[(l_n*56)+23]);
    __m128d a23_0 = _mm_load_sd(&A[96]);
#if defined(__SSE3__) && defined(__AVX__)
    c23_0 = _mm_add_sd(c23_0, _mm_mul_sd(a23_0, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c23_0 = _mm_add_sd(c23_0, _mm_mul_sd(a23_0, b23));
#endif
    _mm_store_sd(&C[(l_n*56)+23], c23_0);
    __m128d c23_1 = _mm_load_sd(&C[(l_n*56)+44]);
    __m128d a23_1 = _mm_load_sd(&A[97]);
#if defined(__SSE3__) && defined(__AVX__)
    c23_1 = _mm_add_sd(c23_1, _mm_mul_sd(a23_1, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c23_1 = _mm_add_sd(c23_1, _mm_mul_sd(a23_1, b23));
#endif
    _mm_store_sd(&C[(l_n*56)+44], c23_1);
#else
    C[(l_n*56)+23] += A[96] * B[(l_n*56)+23];
    C[(l_n*56)+44] += A[97] * B[(l_n*56)+23];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b24 = _mm256_broadcast_sd(&B[(l_n*56)+24]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b24 = _mm_loaddup_pd(&B[(l_n*56)+24]);
#endif
    __m128d c24_0 = _mm_load_sd(&C[(l_n*56)+24]);
    __m128d a24_0 = _mm_load_sd(&A[98]);
#if defined(__SSE3__) && defined(__AVX__)
    c24_0 = _mm_add_sd(c24_0, _mm_mul_sd(a24_0, _mm256_castpd256_pd128(b24)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c24_0 = _mm_add_sd(c24_0, _mm_mul_sd(a24_0, b24));
#endif
    _mm_store_sd(&C[(l_n*56)+24], c24_0);
    __m128d c24_1 = _mm_load_sd(&C[(l_n*56)+45]);
    __m128d a24_1 = _mm_load_sd(&A[99]);
#if defined(__SSE3__) && defined(__AVX__)
    c24_1 = _mm_add_sd(c24_1, _mm_mul_sd(a24_1, _mm256_castpd256_pd128(b24)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c24_1 = _mm_add_sd(c24_1, _mm_mul_sd(a24_1, b24));
#endif
    _mm_store_sd(&C[(l_n*56)+45], c24_1);
#else
    C[(l_n*56)+24] += A[98] * B[(l_n*56)+24];
    C[(l_n*56)+45] += A[99] * B[(l_n*56)+24];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b25 = _mm256_broadcast_sd(&B[(l_n*56)+25]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b25 = _mm_loaddup_pd(&B[(l_n*56)+25]);
#endif
    __m128d c25_0 = _mm_load_sd(&C[(l_n*56)+10]);
    __m128d a25_0 = _mm_load_sd(&A[100]);
#if defined(__SSE3__) && defined(__AVX__)
    c25_0 = _mm_add_sd(c25_0, _mm_mul_sd(a25_0, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c25_0 = _mm_add_sd(c25_0, _mm_mul_sd(a25_0, b25));
#endif
    _mm_store_sd(&C[(l_n*56)+10], c25_0);
    __m128d c25_1 = _mm_load_sd(&C[(l_n*56)+25]);
    __m128d a25_1 = _mm_load_sd(&A[101]);
#if defined(__SSE3__) && defined(__AVX__)
    c25_1 = _mm_add_sd(c25_1, _mm_mul_sd(a25_1, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c25_1 = _mm_add_sd(c25_1, _mm_mul_sd(a25_1, b25));
#endif
    _mm_store_sd(&C[(l_n*56)+25], c25_1);
    __m128d c25_2 = _mm_load_sd(&C[(l_n*56)+46]);
    __m128d a25_2 = _mm_load_sd(&A[102]);
#if defined(__SSE3__) && defined(__AVX__)
    c25_2 = _mm_add_sd(c25_2, _mm_mul_sd(a25_2, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c25_2 = _mm_add_sd(c25_2, _mm_mul_sd(a25_2, b25));
#endif
    _mm_store_sd(&C[(l_n*56)+46], c25_2);
#else
    C[(l_n*56)+10] += A[100] * B[(l_n*56)+25];
    C[(l_n*56)+25] += A[101] * B[(l_n*56)+25];
    C[(l_n*56)+46] += A[102] * B[(l_n*56)+25];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b26 = _mm256_broadcast_sd(&B[(l_n*56)+26]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b26 = _mm_loaddup_pd(&B[(l_n*56)+26]);
#endif
    __m128d c26_0 = _mm_load_sd(&C[(l_n*56)+11]);
    __m128d a26_0 = _mm_load_sd(&A[103]);
#if defined(__SSE3__) && defined(__AVX__)
    c26_0 = _mm_add_sd(c26_0, _mm_mul_sd(a26_0, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c26_0 = _mm_add_sd(c26_0, _mm_mul_sd(a26_0, b26));
#endif
    _mm_store_sd(&C[(l_n*56)+11], c26_0);
    __m128d c26_1 = _mm_load_sd(&C[(l_n*56)+26]);
    __m128d a26_1 = _mm_load_sd(&A[104]);
#if defined(__SSE3__) && defined(__AVX__)
    c26_1 = _mm_add_sd(c26_1, _mm_mul_sd(a26_1, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c26_1 = _mm_add_sd(c26_1, _mm_mul_sd(a26_1, b26));
#endif
    _mm_store_sd(&C[(l_n*56)+26], c26_1);
    __m128d c26_2 = _mm_load_sd(&C[(l_n*56)+47]);
    __m128d a26_2 = _mm_load_sd(&A[105]);
#if defined(__SSE3__) && defined(__AVX__)
    c26_2 = _mm_add_sd(c26_2, _mm_mul_sd(a26_2, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c26_2 = _mm_add_sd(c26_2, _mm_mul_sd(a26_2, b26));
#endif
    _mm_store_sd(&C[(l_n*56)+47], c26_2);
#else
    C[(l_n*56)+11] += A[103] * B[(l_n*56)+26];
    C[(l_n*56)+26] += A[104] * B[(l_n*56)+26];
    C[(l_n*56)+47] += A[105] * B[(l_n*56)+26];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b27 = _mm256_broadcast_sd(&B[(l_n*56)+27]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b27 = _mm_loaddup_pd(&B[(l_n*56)+27]);
#endif
    __m128d c27_0 = _mm_load_sd(&C[(l_n*56)+12]);
    __m128d a27_0 = _mm_load_sd(&A[106]);
#if defined(__SSE3__) && defined(__AVX__)
    c27_0 = _mm_add_sd(c27_0, _mm_mul_sd(a27_0, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c27_0 = _mm_add_sd(c27_0, _mm_mul_sd(a27_0, b27));
#endif
    _mm_store_sd(&C[(l_n*56)+12], c27_0);
    __m128d c27_1 = _mm_load_sd(&C[(l_n*56)+27]);
    __m128d a27_1 = _mm_load_sd(&A[107]);
#if defined(__SSE3__) && defined(__AVX__)
    c27_1 = _mm_add_sd(c27_1, _mm_mul_sd(a27_1, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c27_1 = _mm_add_sd(c27_1, _mm_mul_sd(a27_1, b27));
#endif
    _mm_store_sd(&C[(l_n*56)+27], c27_1);
    __m128d c27_2 = _mm_load_sd(&C[(l_n*56)+48]);
    __m128d a27_2 = _mm_load_sd(&A[108]);
#if defined(__SSE3__) && defined(__AVX__)
    c27_2 = _mm_add_sd(c27_2, _mm_mul_sd(a27_2, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c27_2 = _mm_add_sd(c27_2, _mm_mul_sd(a27_2, b27));
#endif
    _mm_store_sd(&C[(l_n*56)+48], c27_2);
#else
    C[(l_n*56)+12] += A[106] * B[(l_n*56)+27];
    C[(l_n*56)+27] += A[107] * B[(l_n*56)+27];
    C[(l_n*56)+48] += A[108] * B[(l_n*56)+27];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b28 = _mm256_broadcast_sd(&B[(l_n*56)+28]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b28 = _mm_loaddup_pd(&B[(l_n*56)+28]);
#endif
    __m128d c28_0 = _mm_load_sd(&C[(l_n*56)+13]);
    __m128d a28_0 = _mm_load_sd(&A[109]);
#if defined(__SSE3__) && defined(__AVX__)
    c28_0 = _mm_add_sd(c28_0, _mm_mul_sd(a28_0, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c28_0 = _mm_add_sd(c28_0, _mm_mul_sd(a28_0, b28));
#endif
    _mm_store_sd(&C[(l_n*56)+13], c28_0);
    __m128d c28_1 = _mm_load_sd(&C[(l_n*56)+28]);
    __m128d a28_1 = _mm_load_sd(&A[110]);
#if defined(__SSE3__) && defined(__AVX__)
    c28_1 = _mm_add_sd(c28_1, _mm_mul_sd(a28_1, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c28_1 = _mm_add_sd(c28_1, _mm_mul_sd(a28_1, b28));
#endif
    _mm_store_sd(&C[(l_n*56)+28], c28_1);
    __m128d c28_2 = _mm_load_sd(&C[(l_n*56)+49]);
    __m128d a28_2 = _mm_load_sd(&A[111]);
#if defined(__SSE3__) && defined(__AVX__)
    c28_2 = _mm_add_sd(c28_2, _mm_mul_sd(a28_2, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c28_2 = _mm_add_sd(c28_2, _mm_mul_sd(a28_2, b28));
#endif
    _mm_store_sd(&C[(l_n*56)+49], c28_2);
#else
    C[(l_n*56)+13] += A[109] * B[(l_n*56)+28];
    C[(l_n*56)+28] += A[110] * B[(l_n*56)+28];
    C[(l_n*56)+49] += A[111] * B[(l_n*56)+28];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b29 = _mm256_broadcast_sd(&B[(l_n*56)+29]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b29 = _mm_loaddup_pd(&B[(l_n*56)+29]);
#endif
    __m128d c29_0 = _mm_load_sd(&C[(l_n*56)+4]);
    __m128d a29_0 = _mm_load_sd(&A[112]);
#if defined(__SSE3__) && defined(__AVX__)
    c29_0 = _mm_add_sd(c29_0, _mm_mul_sd(a29_0, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c29_0 = _mm_add_sd(c29_0, _mm_mul_sd(a29_0, b29));
#endif
    _mm_store_sd(&C[(l_n*56)+4], c29_0);
    __m128d c29_1 = _mm_load_sd(&C[(l_n*56)+14]);
    __m128d a29_1 = _mm_load_sd(&A[113]);
#if defined(__SSE3__) && defined(__AVX__)
    c29_1 = _mm_add_sd(c29_1, _mm_mul_sd(a29_1, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c29_1 = _mm_add_sd(c29_1, _mm_mul_sd(a29_1, b29));
#endif
    _mm_store_sd(&C[(l_n*56)+14], c29_1);
    __m128d c29_2 = _mm_load_sd(&C[(l_n*56)+29]);
    __m128d a29_2 = _mm_load_sd(&A[114]);
#if defined(__SSE3__) && defined(__AVX__)
    c29_2 = _mm_add_sd(c29_2, _mm_mul_sd(a29_2, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c29_2 = _mm_add_sd(c29_2, _mm_mul_sd(a29_2, b29));
#endif
    _mm_store_sd(&C[(l_n*56)+29], c29_2);
    __m128d c29_3 = _mm_load_sd(&C[(l_n*56)+50]);
    __m128d a29_3 = _mm_load_sd(&A[115]);
#if defined(__SSE3__) && defined(__AVX__)
    c29_3 = _mm_add_sd(c29_3, _mm_mul_sd(a29_3, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c29_3 = _mm_add_sd(c29_3, _mm_mul_sd(a29_3, b29));
#endif
    _mm_store_sd(&C[(l_n*56)+50], c29_3);
#else
    C[(l_n*56)+4] += A[112] * B[(l_n*56)+29];
    C[(l_n*56)+14] += A[113] * B[(l_n*56)+29];
    C[(l_n*56)+29] += A[114] * B[(l_n*56)+29];
    C[(l_n*56)+50] += A[115] * B[(l_n*56)+29];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b30 = _mm256_broadcast_sd(&B[(l_n*56)+30]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b30 = _mm_loaddup_pd(&B[(l_n*56)+30]);
#endif
    __m128d c30_0 = _mm_load_sd(&C[(l_n*56)+5]);
    __m128d a30_0 = _mm_load_sd(&A[116]);
#if defined(__SSE3__) && defined(__AVX__)
    c30_0 = _mm_add_sd(c30_0, _mm_mul_sd(a30_0, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c30_0 = _mm_add_sd(c30_0, _mm_mul_sd(a30_0, b30));
#endif
    _mm_store_sd(&C[(l_n*56)+5], c30_0);
    __m128d c30_1 = _mm_load_sd(&C[(l_n*56)+15]);
    __m128d a30_1 = _mm_load_sd(&A[117]);
#if defined(__SSE3__) && defined(__AVX__)
    c30_1 = _mm_add_sd(c30_1, _mm_mul_sd(a30_1, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c30_1 = _mm_add_sd(c30_1, _mm_mul_sd(a30_1, b30));
#endif
    _mm_store_sd(&C[(l_n*56)+15], c30_1);
    __m128d c30_2 = _mm_load_sd(&C[(l_n*56)+30]);
    __m128d a30_2 = _mm_load_sd(&A[118]);
#if defined(__SSE3__) && defined(__AVX__)
    c30_2 = _mm_add_sd(c30_2, _mm_mul_sd(a30_2, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c30_2 = _mm_add_sd(c30_2, _mm_mul_sd(a30_2, b30));
#endif
    _mm_store_sd(&C[(l_n*56)+30], c30_2);
    __m128d c30_3 = _mm_load_sd(&C[(l_n*56)+51]);
    __m128d a30_3 = _mm_load_sd(&A[119]);
#if defined(__SSE3__) && defined(__AVX__)
    c30_3 = _mm_add_sd(c30_3, _mm_mul_sd(a30_3, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c30_3 = _mm_add_sd(c30_3, _mm_mul_sd(a30_3, b30));
#endif
    _mm_store_sd(&C[(l_n*56)+51], c30_3);
#else
    C[(l_n*56)+5] += A[116] * B[(l_n*56)+30];
    C[(l_n*56)+15] += A[117] * B[(l_n*56)+30];
    C[(l_n*56)+30] += A[118] * B[(l_n*56)+30];
    C[(l_n*56)+51] += A[119] * B[(l_n*56)+30];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b31 = _mm256_broadcast_sd(&B[(l_n*56)+31]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b31 = _mm_loaddup_pd(&B[(l_n*56)+31]);
#endif
    __m128d c31_0 = _mm_load_sd(&C[(l_n*56)+6]);
    __m128d a31_0 = _mm_load_sd(&A[120]);
#if defined(__SSE3__) && defined(__AVX__)
    c31_0 = _mm_add_sd(c31_0, _mm_mul_sd(a31_0, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c31_0 = _mm_add_sd(c31_0, _mm_mul_sd(a31_0, b31));
#endif
    _mm_store_sd(&C[(l_n*56)+6], c31_0);
    __m128d c31_1 = _mm_load_sd(&C[(l_n*56)+16]);
    __m128d a31_1 = _mm_load_sd(&A[121]);
#if defined(__SSE3__) && defined(__AVX__)
    c31_1 = _mm_add_sd(c31_1, _mm_mul_sd(a31_1, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c31_1 = _mm_add_sd(c31_1, _mm_mul_sd(a31_1, b31));
#endif
    _mm_store_sd(&C[(l_n*56)+16], c31_1);
    __m128d c31_2 = _mm_load_sd(&C[(l_n*56)+31]);
    __m128d a31_2 = _mm_load_sd(&A[122]);
#if defined(__SSE3__) && defined(__AVX__)
    c31_2 = _mm_add_sd(c31_2, _mm_mul_sd(a31_2, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c31_2 = _mm_add_sd(c31_2, _mm_mul_sd(a31_2, b31));
#endif
    _mm_store_sd(&C[(l_n*56)+31], c31_2);
    __m128d c31_3 = _mm_load_sd(&C[(l_n*56)+52]);
    __m128d a31_3 = _mm_load_sd(&A[123]);
#if defined(__SSE3__) && defined(__AVX__)
    c31_3 = _mm_add_sd(c31_3, _mm_mul_sd(a31_3, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c31_3 = _mm_add_sd(c31_3, _mm_mul_sd(a31_3, b31));
#endif
    _mm_store_sd(&C[(l_n*56)+52], c31_3);
#else
    C[(l_n*56)+6] += A[120] * B[(l_n*56)+31];
    C[(l_n*56)+16] += A[121] * B[(l_n*56)+31];
    C[(l_n*56)+31] += A[122] * B[(l_n*56)+31];
    C[(l_n*56)+52] += A[123] * B[(l_n*56)+31];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b32 = _mm256_broadcast_sd(&B[(l_n*56)+32]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b32 = _mm_loaddup_pd(&B[(l_n*56)+32]);
#endif
    __m128d c32_0 = _mm_load_sd(&C[(l_n*56)+1]);
    __m128d a32_0 = _mm_load_sd(&A[124]);
#if defined(__SSE3__) && defined(__AVX__)
    c32_0 = _mm_add_sd(c32_0, _mm_mul_sd(a32_0, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c32_0 = _mm_add_sd(c32_0, _mm_mul_sd(a32_0, b32));
#endif
    _mm_store_sd(&C[(l_n*56)+1], c32_0);
    __m128d c32_1 = _mm_load_sd(&C[(l_n*56)+7]);
    __m128d a32_1 = _mm_load_sd(&A[125]);
#if defined(__SSE3__) && defined(__AVX__)
    c32_1 = _mm_add_sd(c32_1, _mm_mul_sd(a32_1, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c32_1 = _mm_add_sd(c32_1, _mm_mul_sd(a32_1, b32));
#endif
    _mm_store_sd(&C[(l_n*56)+7], c32_1);
    __m128d c32_2 = _mm_load_sd(&C[(l_n*56)+17]);
    __m128d a32_2 = _mm_load_sd(&A[126]);
#if defined(__SSE3__) && defined(__AVX__)
    c32_2 = _mm_add_sd(c32_2, _mm_mul_sd(a32_2, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c32_2 = _mm_add_sd(c32_2, _mm_mul_sd(a32_2, b32));
#endif
    _mm_store_sd(&C[(l_n*56)+17], c32_2);
    __m128d c32_3 = _mm_load_sd(&C[(l_n*56)+32]);
    __m128d a32_3 = _mm_load_sd(&A[127]);
#if defined(__SSE3__) && defined(__AVX__)
    c32_3 = _mm_add_sd(c32_3, _mm_mul_sd(a32_3, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c32_3 = _mm_add_sd(c32_3, _mm_mul_sd(a32_3, b32));
#endif
    _mm_store_sd(&C[(l_n*56)+32], c32_3);
    __m128d c32_4 = _mm_load_sd(&C[(l_n*56)+53]);
    __m128d a32_4 = _mm_load_sd(&A[128]);
#if defined(__SSE3__) && defined(__AVX__)
    c32_4 = _mm_add_sd(c32_4, _mm_mul_sd(a32_4, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c32_4 = _mm_add_sd(c32_4, _mm_mul_sd(a32_4, b32));
#endif
    _mm_store_sd(&C[(l_n*56)+53], c32_4);
#else
    C[(l_n*56)+1] += A[124] * B[(l_n*56)+32];
    C[(l_n*56)+7] += A[125] * B[(l_n*56)+32];
    C[(l_n*56)+17] += A[126] * B[(l_n*56)+32];
    C[(l_n*56)+32] += A[127] * B[(l_n*56)+32];
    C[(l_n*56)+53] += A[128] * B[(l_n*56)+32];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b33 = _mm256_broadcast_sd(&B[(l_n*56)+33]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b33 = _mm_loaddup_pd(&B[(l_n*56)+33]);
#endif
    __m128d c33_0 = _mm_load_sd(&C[(l_n*56)+2]);
    __m128d a33_0 = _mm_load_sd(&A[129]);
#if defined(__SSE3__) && defined(__AVX__)
    c33_0 = _mm_add_sd(c33_0, _mm_mul_sd(a33_0, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c33_0 = _mm_add_sd(c33_0, _mm_mul_sd(a33_0, b33));
#endif
    _mm_store_sd(&C[(l_n*56)+2], c33_0);
    __m128d c33_1 = _mm_load_sd(&C[(l_n*56)+8]);
    __m128d a33_1 = _mm_load_sd(&A[130]);
#if defined(__SSE3__) && defined(__AVX__)
    c33_1 = _mm_add_sd(c33_1, _mm_mul_sd(a33_1, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c33_1 = _mm_add_sd(c33_1, _mm_mul_sd(a33_1, b33));
#endif
    _mm_store_sd(&C[(l_n*56)+8], c33_1);
    __m128d c33_2 = _mm_load_sd(&C[(l_n*56)+18]);
    __m128d a33_2 = _mm_load_sd(&A[131]);
#if defined(__SSE3__) && defined(__AVX__)
    c33_2 = _mm_add_sd(c33_2, _mm_mul_sd(a33_2, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c33_2 = _mm_add_sd(c33_2, _mm_mul_sd(a33_2, b33));
#endif
    _mm_store_sd(&C[(l_n*56)+18], c33_2);
    __m128d c33_3 = _mm_load_sd(&C[(l_n*56)+33]);
    __m128d a33_3 = _mm_load_sd(&A[132]);
#if defined(__SSE3__) && defined(__AVX__)
    c33_3 = _mm_add_sd(c33_3, _mm_mul_sd(a33_3, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c33_3 = _mm_add_sd(c33_3, _mm_mul_sd(a33_3, b33));
#endif
    _mm_store_sd(&C[(l_n*56)+33], c33_3);
    __m128d c33_4 = _mm_load_sd(&C[(l_n*56)+54]);
    __m128d a33_4 = _mm_load_sd(&A[133]);
#if defined(__SSE3__) && defined(__AVX__)
    c33_4 = _mm_add_sd(c33_4, _mm_mul_sd(a33_4, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c33_4 = _mm_add_sd(c33_4, _mm_mul_sd(a33_4, b33));
#endif
    _mm_store_sd(&C[(l_n*56)+54], c33_4);
#else
    C[(l_n*56)+2] += A[129] * B[(l_n*56)+33];
    C[(l_n*56)+8] += A[130] * B[(l_n*56)+33];
    C[(l_n*56)+18] += A[131] * B[(l_n*56)+33];
    C[(l_n*56)+33] += A[132] * B[(l_n*56)+33];
    C[(l_n*56)+54] += A[133] * B[(l_n*56)+33];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b34 = _mm256_broadcast_sd(&B[(l_n*56)+34]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b34 = _mm_loaddup_pd(&B[(l_n*56)+34]);
#endif
    __m128d c34_0 = _mm_load_sd(&C[(l_n*56)+0]);
    __m128d a34_0 = _mm_load_sd(&A[134]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_0 = _mm_add_sd(c34_0, _mm_mul_sd(a34_0, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_0 = _mm_add_sd(c34_0, _mm_mul_sd(a34_0, b34));
#endif
    _mm_store_sd(&C[(l_n*56)+0], c34_0);
    __m128d c34_1 = _mm_load_sd(&C[(l_n*56)+3]);
    __m128d a34_1 = _mm_load_sd(&A[135]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_1 = _mm_add_sd(c34_1, _mm_mul_sd(a34_1, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_1 = _mm_add_sd(c34_1, _mm_mul_sd(a34_1, b34));
#endif
    _mm_store_sd(&C[(l_n*56)+3], c34_1);
    __m128d c34_2 = _mm_load_sd(&C[(l_n*56)+9]);
    __m128d a34_2 = _mm_load_sd(&A[136]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_2 = _mm_add_sd(c34_2, _mm_mul_sd(a34_2, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_2 = _mm_add_sd(c34_2, _mm_mul_sd(a34_2, b34));
#endif
    _mm_store_sd(&C[(l_n*56)+9], c34_2);
    __m128d c34_3 = _mm_load_sd(&C[(l_n*56)+19]);
    __m128d a34_3 = _mm_load_sd(&A[137]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_3 = _mm_add_sd(c34_3, _mm_mul_sd(a34_3, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_3 = _mm_add_sd(c34_3, _mm_mul_sd(a34_3, b34));
#endif
    _mm_store_sd(&C[(l_n*56)+19], c34_3);
    __m128d c34_4 = _mm_load_sd(&C[(l_n*56)+34]);
    __m128d a34_4 = _mm_load_sd(&A[138]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_4 = _mm_add_sd(c34_4, _mm_mul_sd(a34_4, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_4 = _mm_add_sd(c34_4, _mm_mul_sd(a34_4, b34));
#endif
    _mm_store_sd(&C[(l_n*56)+34], c34_4);
    __m128d c34_5 = _mm_load_sd(&C[(l_n*56)+55]);
    __m128d a34_5 = _mm_load_sd(&A[139]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_5 = _mm_add_sd(c34_5, _mm_mul_sd(a34_5, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_5 = _mm_add_sd(c34_5, _mm_mul_sd(a34_5, b34));
#endif
    _mm_store_sd(&C[(l_n*56)+55], c34_5);
#else
    C[(l_n*56)+0] += A[134] * B[(l_n*56)+34];
    C[(l_n*56)+3] += A[135] * B[(l_n*56)+34];
    C[(l_n*56)+9] += A[136] * B[(l_n*56)+34];
    C[(l_n*56)+19] += A[137] * B[(l_n*56)+34];
    C[(l_n*56)+34] += A[138] * B[(l_n*56)+34];
    C[(l_n*56)+55] += A[139] * B[(l_n*56)+34];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b35 = _mm256_broadcast_sd(&B[(l_n*56)+35]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b35 = _mm_loaddup_pd(&B[(l_n*56)+35]);
#endif
    __m128d c35_0 = _mm_load_sd(&C[(l_n*56)+35]);
    __m128d a35_0 = _mm_load_sd(&A[140]);
#if defined(__SSE3__) && defined(__AVX__)
    c35_0 = _mm_add_sd(c35_0, _mm_mul_sd(a35_0, _mm256_castpd256_pd128(b35)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c35_0 = _mm_add_sd(c35_0, _mm_mul_sd(a35_0, b35));
#endif
    _mm_store_sd(&C[(l_n*56)+35], c35_0);
#else
    C[(l_n*56)+35] += A[140] * B[(l_n*56)+35];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b36 = _mm256_broadcast_sd(&B[(l_n*56)+36]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b36 = _mm_loaddup_pd(&B[(l_n*56)+36]);
#endif
    __m128d c36_0 = _mm_load_sd(&C[(l_n*56)+36]);
    __m128d a36_0 = _mm_load_sd(&A[141]);
#if defined(__SSE3__) && defined(__AVX__)
    c36_0 = _mm_add_sd(c36_0, _mm_mul_sd(a36_0, _mm256_castpd256_pd128(b36)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c36_0 = _mm_add_sd(c36_0, _mm_mul_sd(a36_0, b36));
#endif
    _mm_store_sd(&C[(l_n*56)+36], c36_0);
#else
    C[(l_n*56)+36] += A[141] * B[(l_n*56)+36];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b37 = _mm256_broadcast_sd(&B[(l_n*56)+37]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b37 = _mm_loaddup_pd(&B[(l_n*56)+37]);
#endif
    __m128d c37_0 = _mm_load_sd(&C[(l_n*56)+37]);
    __m128d a37_0 = _mm_load_sd(&A[142]);
#if defined(__SSE3__) && defined(__AVX__)
    c37_0 = _mm_add_sd(c37_0, _mm_mul_sd(a37_0, _mm256_castpd256_pd128(b37)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c37_0 = _mm_add_sd(c37_0, _mm_mul_sd(a37_0, b37));
#endif
    _mm_store_sd(&C[(l_n*56)+37], c37_0);
#else
    C[(l_n*56)+37] += A[142] * B[(l_n*56)+37];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b38 = _mm256_broadcast_sd(&B[(l_n*56)+38]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b38 = _mm_loaddup_pd(&B[(l_n*56)+38]);
#endif
    __m128d c38_0 = _mm_load_sd(&C[(l_n*56)+38]);
    __m128d a38_0 = _mm_load_sd(&A[143]);
#if defined(__SSE3__) && defined(__AVX__)
    c38_0 = _mm_add_sd(c38_0, _mm_mul_sd(a38_0, _mm256_castpd256_pd128(b38)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c38_0 = _mm_add_sd(c38_0, _mm_mul_sd(a38_0, b38));
#endif
    _mm_store_sd(&C[(l_n*56)+38], c38_0);
#else
    C[(l_n*56)+38] += A[143] * B[(l_n*56)+38];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b39 = _mm256_broadcast_sd(&B[(l_n*56)+39]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b39 = _mm_loaddup_pd(&B[(l_n*56)+39]);
#endif
    __m128d c39_0 = _mm_load_sd(&C[(l_n*56)+39]);
    __m128d a39_0 = _mm_load_sd(&A[144]);
#if defined(__SSE3__) && defined(__AVX__)
    c39_0 = _mm_add_sd(c39_0, _mm_mul_sd(a39_0, _mm256_castpd256_pd128(b39)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c39_0 = _mm_add_sd(c39_0, _mm_mul_sd(a39_0, b39));
#endif
    _mm_store_sd(&C[(l_n*56)+39], c39_0);
#else
    C[(l_n*56)+39] += A[144] * B[(l_n*56)+39];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b40 = _mm256_broadcast_sd(&B[(l_n*56)+40]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b40 = _mm_loaddup_pd(&B[(l_n*56)+40]);
#endif
    __m128d c40_0 = _mm_load_sd(&C[(l_n*56)+40]);
    __m128d a40_0 = _mm_load_sd(&A[145]);
#if defined(__SSE3__) && defined(__AVX__)
    c40_0 = _mm_add_sd(c40_0, _mm_mul_sd(a40_0, _mm256_castpd256_pd128(b40)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c40_0 = _mm_add_sd(c40_0, _mm_mul_sd(a40_0, b40));
#endif
    _mm_store_sd(&C[(l_n*56)+40], c40_0);
#else
    C[(l_n*56)+40] += A[145] * B[(l_n*56)+40];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b41 = _mm256_broadcast_sd(&B[(l_n*56)+41]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b41 = _mm_loaddup_pd(&B[(l_n*56)+41]);
#endif
    __m128d c41_0 = _mm_load_sd(&C[(l_n*56)+20]);
    __m128d a41_0 = _mm_load_sd(&A[146]);
#if defined(__SSE3__) && defined(__AVX__)
    c41_0 = _mm_add_sd(c41_0, _mm_mul_sd(a41_0, _mm256_castpd256_pd128(b41)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c41_0 = _mm_add_sd(c41_0, _mm_mul_sd(a41_0, b41));
#endif
    _mm_store_sd(&C[(l_n*56)+20], c41_0);
    __m128d c41_1 = _mm_load_sd(&C[(l_n*56)+41]);
    __m128d a41_1 = _mm_load_sd(&A[147]);
#if defined(__SSE3__) && defined(__AVX__)
    c41_1 = _mm_add_sd(c41_1, _mm_mul_sd(a41_1, _mm256_castpd256_pd128(b41)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c41_1 = _mm_add_sd(c41_1, _mm_mul_sd(a41_1, b41));
#endif
    _mm_store_sd(&C[(l_n*56)+41], c41_1);
#else
    C[(l_n*56)+20] += A[146] * B[(l_n*56)+41];
    C[(l_n*56)+41] += A[147] * B[(l_n*56)+41];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b42 = _mm256_broadcast_sd(&B[(l_n*56)+42]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b42 = _mm_loaddup_pd(&B[(l_n*56)+42]);
#endif
    __m128d c42_0 = _mm_load_sd(&C[(l_n*56)+21]);
    __m128d a42_0 = _mm_load_sd(&A[148]);
#if defined(__SSE3__) && defined(__AVX__)
    c42_0 = _mm_add_sd(c42_0, _mm_mul_sd(a42_0, _mm256_castpd256_pd128(b42)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c42_0 = _mm_add_sd(c42_0, _mm_mul_sd(a42_0, b42));
#endif
    _mm_store_sd(&C[(l_n*56)+21], c42_0);
    __m128d c42_1 = _mm_load_sd(&C[(l_n*56)+42]);
    __m128d a42_1 = _mm_load_sd(&A[149]);
#if defined(__SSE3__) && defined(__AVX__)
    c42_1 = _mm_add_sd(c42_1, _mm_mul_sd(a42_1, _mm256_castpd256_pd128(b42)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c42_1 = _mm_add_sd(c42_1, _mm_mul_sd(a42_1, b42));
#endif
    _mm_store_sd(&C[(l_n*56)+42], c42_1);
#else
    C[(l_n*56)+21] += A[148] * B[(l_n*56)+42];
    C[(l_n*56)+42] += A[149] * B[(l_n*56)+42];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b43 = _mm256_broadcast_sd(&B[(l_n*56)+43]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b43 = _mm_loaddup_pd(&B[(l_n*56)+43]);
#endif
    __m128d c43_0 = _mm_load_sd(&C[(l_n*56)+22]);
    __m128d a43_0 = _mm_load_sd(&A[150]);
#if defined(__SSE3__) && defined(__AVX__)
    c43_0 = _mm_add_sd(c43_0, _mm_mul_sd(a43_0, _mm256_castpd256_pd128(b43)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c43_0 = _mm_add_sd(c43_0, _mm_mul_sd(a43_0, b43));
#endif
    _mm_store_sd(&C[(l_n*56)+22], c43_0);
    __m128d c43_1 = _mm_load_sd(&C[(l_n*56)+43]);
    __m128d a43_1 = _mm_load_sd(&A[151]);
#if defined(__SSE3__) && defined(__AVX__)
    c43_1 = _mm_add_sd(c43_1, _mm_mul_sd(a43_1, _mm256_castpd256_pd128(b43)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c43_1 = _mm_add_sd(c43_1, _mm_mul_sd(a43_1, b43));
#endif
    _mm_store_sd(&C[(l_n*56)+43], c43_1);
#else
    C[(l_n*56)+22] += A[150] * B[(l_n*56)+43];
    C[(l_n*56)+43] += A[151] * B[(l_n*56)+43];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b44 = _mm256_broadcast_sd(&B[(l_n*56)+44]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b44 = _mm_loaddup_pd(&B[(l_n*56)+44]);
#endif
    __m128d c44_0 = _mm_load_sd(&C[(l_n*56)+23]);
    __m128d a44_0 = _mm_load_sd(&A[152]);
#if defined(__SSE3__) && defined(__AVX__)
    c44_0 = _mm_add_sd(c44_0, _mm_mul_sd(a44_0, _mm256_castpd256_pd128(b44)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c44_0 = _mm_add_sd(c44_0, _mm_mul_sd(a44_0, b44));
#endif
    _mm_store_sd(&C[(l_n*56)+23], c44_0);
    __m128d c44_1 = _mm_load_sd(&C[(l_n*56)+44]);
    __m128d a44_1 = _mm_load_sd(&A[153]);
#if defined(__SSE3__) && defined(__AVX__)
    c44_1 = _mm_add_sd(c44_1, _mm_mul_sd(a44_1, _mm256_castpd256_pd128(b44)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c44_1 = _mm_add_sd(c44_1, _mm_mul_sd(a44_1, b44));
#endif
    _mm_store_sd(&C[(l_n*56)+44], c44_1);
#else
    C[(l_n*56)+23] += A[152] * B[(l_n*56)+44];
    C[(l_n*56)+44] += A[153] * B[(l_n*56)+44];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b45 = _mm256_broadcast_sd(&B[(l_n*56)+45]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b45 = _mm_loaddup_pd(&B[(l_n*56)+45]);
#endif
    __m128d c45_0 = _mm_load_sd(&C[(l_n*56)+24]);
    __m128d a45_0 = _mm_load_sd(&A[154]);
#if defined(__SSE3__) && defined(__AVX__)
    c45_0 = _mm_add_sd(c45_0, _mm_mul_sd(a45_0, _mm256_castpd256_pd128(b45)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c45_0 = _mm_add_sd(c45_0, _mm_mul_sd(a45_0, b45));
#endif
    _mm_store_sd(&C[(l_n*56)+24], c45_0);
    __m128d c45_1 = _mm_load_sd(&C[(l_n*56)+45]);
    __m128d a45_1 = _mm_load_sd(&A[155]);
#if defined(__SSE3__) && defined(__AVX__)
    c45_1 = _mm_add_sd(c45_1, _mm_mul_sd(a45_1, _mm256_castpd256_pd128(b45)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c45_1 = _mm_add_sd(c45_1, _mm_mul_sd(a45_1, b45));
#endif
    _mm_store_sd(&C[(l_n*56)+45], c45_1);
#else
    C[(l_n*56)+24] += A[154] * B[(l_n*56)+45];
    C[(l_n*56)+45] += A[155] * B[(l_n*56)+45];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b46 = _mm256_broadcast_sd(&B[(l_n*56)+46]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b46 = _mm_loaddup_pd(&B[(l_n*56)+46]);
#endif
    __m128d c46_0 = _mm_load_sd(&C[(l_n*56)+10]);
    __m128d a46_0 = _mm_load_sd(&A[156]);
#if defined(__SSE3__) && defined(__AVX__)
    c46_0 = _mm_add_sd(c46_0, _mm_mul_sd(a46_0, _mm256_castpd256_pd128(b46)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c46_0 = _mm_add_sd(c46_0, _mm_mul_sd(a46_0, b46));
#endif
    _mm_store_sd(&C[(l_n*56)+10], c46_0);
    __m128d c46_1 = _mm_load_sd(&C[(l_n*56)+25]);
    __m128d a46_1 = _mm_load_sd(&A[157]);
#if defined(__SSE3__) && defined(__AVX__)
    c46_1 = _mm_add_sd(c46_1, _mm_mul_sd(a46_1, _mm256_castpd256_pd128(b46)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c46_1 = _mm_add_sd(c46_1, _mm_mul_sd(a46_1, b46));
#endif
    _mm_store_sd(&C[(l_n*56)+25], c46_1);
    __m128d c46_2 = _mm_load_sd(&C[(l_n*56)+46]);
    __m128d a46_2 = _mm_load_sd(&A[158]);
#if defined(__SSE3__) && defined(__AVX__)
    c46_2 = _mm_add_sd(c46_2, _mm_mul_sd(a46_2, _mm256_castpd256_pd128(b46)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c46_2 = _mm_add_sd(c46_2, _mm_mul_sd(a46_2, b46));
#endif
    _mm_store_sd(&C[(l_n*56)+46], c46_2);
#else
    C[(l_n*56)+10] += A[156] * B[(l_n*56)+46];
    C[(l_n*56)+25] += A[157] * B[(l_n*56)+46];
    C[(l_n*56)+46] += A[158] * B[(l_n*56)+46];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b47 = _mm256_broadcast_sd(&B[(l_n*56)+47]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b47 = _mm_loaddup_pd(&B[(l_n*56)+47]);
#endif
    __m128d c47_0 = _mm_load_sd(&C[(l_n*56)+11]);
    __m128d a47_0 = _mm_load_sd(&A[159]);
#if defined(__SSE3__) && defined(__AVX__)
    c47_0 = _mm_add_sd(c47_0, _mm_mul_sd(a47_0, _mm256_castpd256_pd128(b47)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c47_0 = _mm_add_sd(c47_0, _mm_mul_sd(a47_0, b47));
#endif
    _mm_store_sd(&C[(l_n*56)+11], c47_0);
    __m128d c47_1 = _mm_load_sd(&C[(l_n*56)+26]);
    __m128d a47_1 = _mm_load_sd(&A[160]);
#if defined(__SSE3__) && defined(__AVX__)
    c47_1 = _mm_add_sd(c47_1, _mm_mul_sd(a47_1, _mm256_castpd256_pd128(b47)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c47_1 = _mm_add_sd(c47_1, _mm_mul_sd(a47_1, b47));
#endif
    _mm_store_sd(&C[(l_n*56)+26], c47_1);
    __m128d c47_2 = _mm_load_sd(&C[(l_n*56)+47]);
    __m128d a47_2 = _mm_load_sd(&A[161]);
#if defined(__SSE3__) && defined(__AVX__)
    c47_2 = _mm_add_sd(c47_2, _mm_mul_sd(a47_2, _mm256_castpd256_pd128(b47)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c47_2 = _mm_add_sd(c47_2, _mm_mul_sd(a47_2, b47));
#endif
    _mm_store_sd(&C[(l_n*56)+47], c47_2);
#else
    C[(l_n*56)+11] += A[159] * B[(l_n*56)+47];
    C[(l_n*56)+26] += A[160] * B[(l_n*56)+47];
    C[(l_n*56)+47] += A[161] * B[(l_n*56)+47];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b48 = _mm256_broadcast_sd(&B[(l_n*56)+48]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b48 = _mm_loaddup_pd(&B[(l_n*56)+48]);
#endif
    __m128d c48_0 = _mm_load_sd(&C[(l_n*56)+12]);
    __m128d a48_0 = _mm_load_sd(&A[162]);
#if defined(__SSE3__) && defined(__AVX__)
    c48_0 = _mm_add_sd(c48_0, _mm_mul_sd(a48_0, _mm256_castpd256_pd128(b48)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c48_0 = _mm_add_sd(c48_0, _mm_mul_sd(a48_0, b48));
#endif
    _mm_store_sd(&C[(l_n*56)+12], c48_0);
    __m128d c48_1 = _mm_load_sd(&C[(l_n*56)+27]);
    __m128d a48_1 = _mm_load_sd(&A[163]);
#if defined(__SSE3__) && defined(__AVX__)
    c48_1 = _mm_add_sd(c48_1, _mm_mul_sd(a48_1, _mm256_castpd256_pd128(b48)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c48_1 = _mm_add_sd(c48_1, _mm_mul_sd(a48_1, b48));
#endif
    _mm_store_sd(&C[(l_n*56)+27], c48_1);
    __m128d c48_2 = _mm_load_sd(&C[(l_n*56)+48]);
    __m128d a48_2 = _mm_load_sd(&A[164]);
#if defined(__SSE3__) && defined(__AVX__)
    c48_2 = _mm_add_sd(c48_2, _mm_mul_sd(a48_2, _mm256_castpd256_pd128(b48)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c48_2 = _mm_add_sd(c48_2, _mm_mul_sd(a48_2, b48));
#endif
    _mm_store_sd(&C[(l_n*56)+48], c48_2);
#else
    C[(l_n*56)+12] += A[162] * B[(l_n*56)+48];
    C[(l_n*56)+27] += A[163] * B[(l_n*56)+48];
    C[(l_n*56)+48] += A[164] * B[(l_n*56)+48];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b49 = _mm256_broadcast_sd(&B[(l_n*56)+49]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b49 = _mm_loaddup_pd(&B[(l_n*56)+49]);
#endif
    __m128d c49_0 = _mm_load_sd(&C[(l_n*56)+13]);
    __m128d a49_0 = _mm_load_sd(&A[165]);
#if defined(__SSE3__) && defined(__AVX__)
    c49_0 = _mm_add_sd(c49_0, _mm_mul_sd(a49_0, _mm256_castpd256_pd128(b49)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c49_0 = _mm_add_sd(c49_0, _mm_mul_sd(a49_0, b49));
#endif
    _mm_store_sd(&C[(l_n*56)+13], c49_0);
    __m128d c49_1 = _mm_load_sd(&C[(l_n*56)+28]);
    __m128d a49_1 = _mm_load_sd(&A[166]);
#if defined(__SSE3__) && defined(__AVX__)
    c49_1 = _mm_add_sd(c49_1, _mm_mul_sd(a49_1, _mm256_castpd256_pd128(b49)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c49_1 = _mm_add_sd(c49_1, _mm_mul_sd(a49_1, b49));
#endif
    _mm_store_sd(&C[(l_n*56)+28], c49_1);
    __m128d c49_2 = _mm_load_sd(&C[(l_n*56)+49]);
    __m128d a49_2 = _mm_load_sd(&A[167]);
#if defined(__SSE3__) && defined(__AVX__)
    c49_2 = _mm_add_sd(c49_2, _mm_mul_sd(a49_2, _mm256_castpd256_pd128(b49)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c49_2 = _mm_add_sd(c49_2, _mm_mul_sd(a49_2, b49));
#endif
    _mm_store_sd(&C[(l_n*56)+49], c49_2);
#else
    C[(l_n*56)+13] += A[165] * B[(l_n*56)+49];
    C[(l_n*56)+28] += A[166] * B[(l_n*56)+49];
    C[(l_n*56)+49] += A[167] * B[(l_n*56)+49];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b50 = _mm256_broadcast_sd(&B[(l_n*56)+50]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b50 = _mm_loaddup_pd(&B[(l_n*56)+50]);
#endif
    __m128d c50_0 = _mm_load_sd(&C[(l_n*56)+4]);
    __m128d a50_0 = _mm_load_sd(&A[168]);
#if defined(__SSE3__) && defined(__AVX__)
    c50_0 = _mm_add_sd(c50_0, _mm_mul_sd(a50_0, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c50_0 = _mm_add_sd(c50_0, _mm_mul_sd(a50_0, b50));
#endif
    _mm_store_sd(&C[(l_n*56)+4], c50_0);
    __m128d c50_1 = _mm_load_sd(&C[(l_n*56)+14]);
    __m128d a50_1 = _mm_load_sd(&A[169]);
#if defined(__SSE3__) && defined(__AVX__)
    c50_1 = _mm_add_sd(c50_1, _mm_mul_sd(a50_1, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c50_1 = _mm_add_sd(c50_1, _mm_mul_sd(a50_1, b50));
#endif
    _mm_store_sd(&C[(l_n*56)+14], c50_1);
    __m128d c50_2 = _mm_load_sd(&C[(l_n*56)+29]);
    __m128d a50_2 = _mm_load_sd(&A[170]);
#if defined(__SSE3__) && defined(__AVX__)
    c50_2 = _mm_add_sd(c50_2, _mm_mul_sd(a50_2, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c50_2 = _mm_add_sd(c50_2, _mm_mul_sd(a50_2, b50));
#endif
    _mm_store_sd(&C[(l_n*56)+29], c50_2);
    __m128d c50_3 = _mm_load_sd(&C[(l_n*56)+50]);
    __m128d a50_3 = _mm_load_sd(&A[171]);
#if defined(__SSE3__) && defined(__AVX__)
    c50_3 = _mm_add_sd(c50_3, _mm_mul_sd(a50_3, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c50_3 = _mm_add_sd(c50_3, _mm_mul_sd(a50_3, b50));
#endif
    _mm_store_sd(&C[(l_n*56)+50], c50_3);
#else
    C[(l_n*56)+4] += A[168] * B[(l_n*56)+50];
    C[(l_n*56)+14] += A[169] * B[(l_n*56)+50];
    C[(l_n*56)+29] += A[170] * B[(l_n*56)+50];
    C[(l_n*56)+50] += A[171] * B[(l_n*56)+50];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b51 = _mm256_broadcast_sd(&B[(l_n*56)+51]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b51 = _mm_loaddup_pd(&B[(l_n*56)+51]);
#endif
    __m128d c51_0 = _mm_load_sd(&C[(l_n*56)+5]);
    __m128d a51_0 = _mm_load_sd(&A[172]);
#if defined(__SSE3__) && defined(__AVX__)
    c51_0 = _mm_add_sd(c51_0, _mm_mul_sd(a51_0, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c51_0 = _mm_add_sd(c51_0, _mm_mul_sd(a51_0, b51));
#endif
    _mm_store_sd(&C[(l_n*56)+5], c51_0);
    __m128d c51_1 = _mm_load_sd(&C[(l_n*56)+15]);
    __m128d a51_1 = _mm_load_sd(&A[173]);
#if defined(__SSE3__) && defined(__AVX__)
    c51_1 = _mm_add_sd(c51_1, _mm_mul_sd(a51_1, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c51_1 = _mm_add_sd(c51_1, _mm_mul_sd(a51_1, b51));
#endif
    _mm_store_sd(&C[(l_n*56)+15], c51_1);
    __m128d c51_2 = _mm_load_sd(&C[(l_n*56)+30]);
    __m128d a51_2 = _mm_load_sd(&A[174]);
#if defined(__SSE3__) && defined(__AVX__)
    c51_2 = _mm_add_sd(c51_2, _mm_mul_sd(a51_2, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c51_2 = _mm_add_sd(c51_2, _mm_mul_sd(a51_2, b51));
#endif
    _mm_store_sd(&C[(l_n*56)+30], c51_2);
    __m128d c51_3 = _mm_load_sd(&C[(l_n*56)+51]);
    __m128d a51_3 = _mm_load_sd(&A[175]);
#if defined(__SSE3__) && defined(__AVX__)
    c51_3 = _mm_add_sd(c51_3, _mm_mul_sd(a51_3, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c51_3 = _mm_add_sd(c51_3, _mm_mul_sd(a51_3, b51));
#endif
    _mm_store_sd(&C[(l_n*56)+51], c51_3);
#else
    C[(l_n*56)+5] += A[172] * B[(l_n*56)+51];
    C[(l_n*56)+15] += A[173] * B[(l_n*56)+51];
    C[(l_n*56)+30] += A[174] * B[(l_n*56)+51];
    C[(l_n*56)+51] += A[175] * B[(l_n*56)+51];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b52 = _mm256_broadcast_sd(&B[(l_n*56)+52]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b52 = _mm_loaddup_pd(&B[(l_n*56)+52]);
#endif
    __m128d c52_0 = _mm_load_sd(&C[(l_n*56)+6]);
    __m128d a52_0 = _mm_load_sd(&A[176]);
#if defined(__SSE3__) && defined(__AVX__)
    c52_0 = _mm_add_sd(c52_0, _mm_mul_sd(a52_0, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c52_0 = _mm_add_sd(c52_0, _mm_mul_sd(a52_0, b52));
#endif
    _mm_store_sd(&C[(l_n*56)+6], c52_0);
    __m128d c52_1 = _mm_load_sd(&C[(l_n*56)+16]);
    __m128d a52_1 = _mm_load_sd(&A[177]);
#if defined(__SSE3__) && defined(__AVX__)
    c52_1 = _mm_add_sd(c52_1, _mm_mul_sd(a52_1, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c52_1 = _mm_add_sd(c52_1, _mm_mul_sd(a52_1, b52));
#endif
    _mm_store_sd(&C[(l_n*56)+16], c52_1);
    __m128d c52_2 = _mm_load_sd(&C[(l_n*56)+31]);
    __m128d a52_2 = _mm_load_sd(&A[178]);
#if defined(__SSE3__) && defined(__AVX__)
    c52_2 = _mm_add_sd(c52_2, _mm_mul_sd(a52_2, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c52_2 = _mm_add_sd(c52_2, _mm_mul_sd(a52_2, b52));
#endif
    _mm_store_sd(&C[(l_n*56)+31], c52_2);
    __m128d c52_3 = _mm_load_sd(&C[(l_n*56)+52]);
    __m128d a52_3 = _mm_load_sd(&A[179]);
#if defined(__SSE3__) && defined(__AVX__)
    c52_3 = _mm_add_sd(c52_3, _mm_mul_sd(a52_3, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c52_3 = _mm_add_sd(c52_3, _mm_mul_sd(a52_3, b52));
#endif
    _mm_store_sd(&C[(l_n*56)+52], c52_3);
#else
    C[(l_n*56)+6] += A[176] * B[(l_n*56)+52];
    C[(l_n*56)+16] += A[177] * B[(l_n*56)+52];
    C[(l_n*56)+31] += A[178] * B[(l_n*56)+52];
    C[(l_n*56)+52] += A[179] * B[(l_n*56)+52];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b53 = _mm256_broadcast_sd(&B[(l_n*56)+53]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b53 = _mm_loaddup_pd(&B[(l_n*56)+53]);
#endif
    __m128d c53_0 = _mm_load_sd(&C[(l_n*56)+1]);
    __m128d a53_0 = _mm_load_sd(&A[180]);
#if defined(__SSE3__) && defined(__AVX__)
    c53_0 = _mm_add_sd(c53_0, _mm_mul_sd(a53_0, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c53_0 = _mm_add_sd(c53_0, _mm_mul_sd(a53_0, b53));
#endif
    _mm_store_sd(&C[(l_n*56)+1], c53_0);
    __m128d c53_1 = _mm_load_sd(&C[(l_n*56)+7]);
    __m128d a53_1 = _mm_load_sd(&A[181]);
#if defined(__SSE3__) && defined(__AVX__)
    c53_1 = _mm_add_sd(c53_1, _mm_mul_sd(a53_1, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c53_1 = _mm_add_sd(c53_1, _mm_mul_sd(a53_1, b53));
#endif
    _mm_store_sd(&C[(l_n*56)+7], c53_1);
    __m128d c53_2 = _mm_load_sd(&C[(l_n*56)+17]);
    __m128d a53_2 = _mm_load_sd(&A[182]);
#if defined(__SSE3__) && defined(__AVX__)
    c53_2 = _mm_add_sd(c53_2, _mm_mul_sd(a53_2, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c53_2 = _mm_add_sd(c53_2, _mm_mul_sd(a53_2, b53));
#endif
    _mm_store_sd(&C[(l_n*56)+17], c53_2);
    __m128d c53_3 = _mm_load_sd(&C[(l_n*56)+32]);
    __m128d a53_3 = _mm_load_sd(&A[183]);
#if defined(__SSE3__) && defined(__AVX__)
    c53_3 = _mm_add_sd(c53_3, _mm_mul_sd(a53_3, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c53_3 = _mm_add_sd(c53_3, _mm_mul_sd(a53_3, b53));
#endif
    _mm_store_sd(&C[(l_n*56)+32], c53_3);
    __m128d c53_4 = _mm_load_sd(&C[(l_n*56)+53]);
    __m128d a53_4 = _mm_load_sd(&A[184]);
#if defined(__SSE3__) && defined(__AVX__)
    c53_4 = _mm_add_sd(c53_4, _mm_mul_sd(a53_4, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c53_4 = _mm_add_sd(c53_4, _mm_mul_sd(a53_4, b53));
#endif
    _mm_store_sd(&C[(l_n*56)+53], c53_4);
#else
    C[(l_n*56)+1] += A[180] * B[(l_n*56)+53];
    C[(l_n*56)+7] += A[181] * B[(l_n*56)+53];
    C[(l_n*56)+17] += A[182] * B[(l_n*56)+53];
    C[(l_n*56)+32] += A[183] * B[(l_n*56)+53];
    C[(l_n*56)+53] += A[184] * B[(l_n*56)+53];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b54 = _mm256_broadcast_sd(&B[(l_n*56)+54]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b54 = _mm_loaddup_pd(&B[(l_n*56)+54]);
#endif
    __m128d c54_0 = _mm_load_sd(&C[(l_n*56)+2]);
    __m128d a54_0 = _mm_load_sd(&A[185]);
#if defined(__SSE3__) && defined(__AVX__)
    c54_0 = _mm_add_sd(c54_0, _mm_mul_sd(a54_0, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c54_0 = _mm_add_sd(c54_0, _mm_mul_sd(a54_0, b54));
#endif
    _mm_store_sd(&C[(l_n*56)+2], c54_0);
    __m128d c54_1 = _mm_load_sd(&C[(l_n*56)+8]);
    __m128d a54_1 = _mm_load_sd(&A[186]);
#if defined(__SSE3__) && defined(__AVX__)
    c54_1 = _mm_add_sd(c54_1, _mm_mul_sd(a54_1, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c54_1 = _mm_add_sd(c54_1, _mm_mul_sd(a54_1, b54));
#endif
    _mm_store_sd(&C[(l_n*56)+8], c54_1);
    __m128d c54_2 = _mm_load_sd(&C[(l_n*56)+18]);
    __m128d a54_2 = _mm_load_sd(&A[187]);
#if defined(__SSE3__) && defined(__AVX__)
    c54_2 = _mm_add_sd(c54_2, _mm_mul_sd(a54_2, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c54_2 = _mm_add_sd(c54_2, _mm_mul_sd(a54_2, b54));
#endif
    _mm_store_sd(&C[(l_n*56)+18], c54_2);
    __m128d c54_3 = _mm_load_sd(&C[(l_n*56)+33]);
    __m128d a54_3 = _mm_load_sd(&A[188]);
#if defined(__SSE3__) && defined(__AVX__)
    c54_3 = _mm_add_sd(c54_3, _mm_mul_sd(a54_3, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c54_3 = _mm_add_sd(c54_3, _mm_mul_sd(a54_3, b54));
#endif
    _mm_store_sd(&C[(l_n*56)+33], c54_3);
    __m128d c54_4 = _mm_load_sd(&C[(l_n*56)+54]);
    __m128d a54_4 = _mm_load_sd(&A[189]);
#if defined(__SSE3__) && defined(__AVX__)
    c54_4 = _mm_add_sd(c54_4, _mm_mul_sd(a54_4, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c54_4 = _mm_add_sd(c54_4, _mm_mul_sd(a54_4, b54));
#endif
    _mm_store_sd(&C[(l_n*56)+54], c54_4);
#else
    C[(l_n*56)+2] += A[185] * B[(l_n*56)+54];
    C[(l_n*56)+8] += A[186] * B[(l_n*56)+54];
    C[(l_n*56)+18] += A[187] * B[(l_n*56)+54];
    C[(l_n*56)+33] += A[188] * B[(l_n*56)+54];
    C[(l_n*56)+54] += A[189] * B[(l_n*56)+54];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b55 = _mm256_broadcast_sd(&B[(l_n*56)+55]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b55 = _mm_loaddup_pd(&B[(l_n*56)+55]);
#endif
    __m128d c55_0 = _mm_load_sd(&C[(l_n*56)+0]);
    __m128d a55_0 = _mm_load_sd(&A[190]);
#if defined(__SSE3__) && defined(__AVX__)
    c55_0 = _mm_add_sd(c55_0, _mm_mul_sd(a55_0, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c55_0 = _mm_add_sd(c55_0, _mm_mul_sd(a55_0, b55));
#endif
    _mm_store_sd(&C[(l_n*56)+0], c55_0);
    __m128d c55_1 = _mm_load_sd(&C[(l_n*56)+3]);
    __m128d a55_1 = _mm_load_sd(&A[191]);
#if defined(__SSE3__) && defined(__AVX__)
    c55_1 = _mm_add_sd(c55_1, _mm_mul_sd(a55_1, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c55_1 = _mm_add_sd(c55_1, _mm_mul_sd(a55_1, b55));
#endif
    _mm_store_sd(&C[(l_n*56)+3], c55_1);
    __m128d c55_2 = _mm_load_sd(&C[(l_n*56)+9]);
    __m128d a55_2 = _mm_load_sd(&A[192]);
#if defined(__SSE3__) && defined(__AVX__)
    c55_2 = _mm_add_sd(c55_2, _mm_mul_sd(a55_2, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c55_2 = _mm_add_sd(c55_2, _mm_mul_sd(a55_2, b55));
#endif
    _mm_store_sd(&C[(l_n*56)+9], c55_2);
    __m128d c55_3 = _mm_load_sd(&C[(l_n*56)+19]);
    __m128d a55_3 = _mm_load_sd(&A[193]);
#if defined(__SSE3__) && defined(__AVX__)
    c55_3 = _mm_add_sd(c55_3, _mm_mul_sd(a55_3, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c55_3 = _mm_add_sd(c55_3, _mm_mul_sd(a55_3, b55));
#endif
    _mm_store_sd(&C[(l_n*56)+19], c55_3);
    __m128d c55_4 = _mm_load_sd(&C[(l_n*56)+34]);
    __m128d a55_4 = _mm_load_sd(&A[194]);
#if defined(__SSE3__) && defined(__AVX__)
    c55_4 = _mm_add_sd(c55_4, _mm_mul_sd(a55_4, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c55_4 = _mm_add_sd(c55_4, _mm_mul_sd(a55_4, b55));
#endif
    _mm_store_sd(&C[(l_n*56)+34], c55_4);
    __m128d c55_5 = _mm_load_sd(&C[(l_n*56)+55]);
    __m128d a55_5 = _mm_load_sd(&A[195]);
#if defined(__SSE3__) && defined(__AVX__)
    c55_5 = _mm_add_sd(c55_5, _mm_mul_sd(a55_5, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c55_5 = _mm_add_sd(c55_5, _mm_mul_sd(a55_5, b55));
#endif
    _mm_store_sd(&C[(l_n*56)+55], c55_5);
#else
    C[(l_n*56)+0] += A[190] * B[(l_n*56)+55];
    C[(l_n*56)+3] += A[191] * B[(l_n*56)+55];
    C[(l_n*56)+9] += A[192] * B[(l_n*56)+55];
    C[(l_n*56)+19] += A[193] * B[(l_n*56)+55];
    C[(l_n*56)+34] += A[194] * B[(l_n*56)+55];
    C[(l_n*56)+55] += A[195] * B[(l_n*56)+55];
#endif

  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 3528;
#endif
}

void dsparse_fP113DivM_m56_n9_k56_ldAna6_ldB56_ldC56_beta0_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_n = 0;
  #pragma nounroll_and_jam
  for ( l_n = 0; l_n < 9; l_n++) {
    unsigned int l_m = 0;
   #pragma simd
    for ( l_m = 0; l_m < 56; l_m++) {
      C[(l_n*56)+l_m] = 0.0;
    }
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b0 = _mm256_broadcast_sd(&B[(l_n*56)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b0 = _mm_loaddup_pd(&B[(l_n*56)+0]);
#endif
    __m128d c0_0 = _mm_load_sd(&C[(l_n*56)+0]);
    __m128d a0_0 = _mm_load_sd(&A[0]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, b0));
#endif
    _mm_store_sd(&C[(l_n*56)+0], c0_0);
    __m128d c0_1 = _mm_load_sd(&C[(l_n*56)+3]);
    __m128d a0_1 = _mm_load_sd(&A[1]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, b0));
#endif
    _mm_store_sd(&C[(l_n*56)+3], c0_1);
    __m128d c0_2 = _mm_load_sd(&C[(l_n*56)+9]);
    __m128d a0_2 = _mm_load_sd(&A[2]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, b0));
#endif
    _mm_store_sd(&C[(l_n*56)+9], c0_2);
    __m128d c0_3 = _mm_load_sd(&C[(l_n*56)+19]);
    __m128d a0_3 = _mm_load_sd(&A[3]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_3 = _mm_add_sd(c0_3, _mm_mul_sd(a0_3, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_3 = _mm_add_sd(c0_3, _mm_mul_sd(a0_3, b0));
#endif
    _mm_store_sd(&C[(l_n*56)+19], c0_3);
    __m128d c0_4 = _mm_load_sd(&C[(l_n*56)+34]);
    __m128d a0_4 = _mm_load_sd(&A[4]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_4 = _mm_add_sd(c0_4, _mm_mul_sd(a0_4, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_4 = _mm_add_sd(c0_4, _mm_mul_sd(a0_4, b0));
#endif
    _mm_store_sd(&C[(l_n*56)+34], c0_4);
    __m128d c0_5 = _mm_load_sd(&C[(l_n*56)+55]);
    __m128d a0_5 = _mm_load_sd(&A[5]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_5 = _mm_add_sd(c0_5, _mm_mul_sd(a0_5, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_5 = _mm_add_sd(c0_5, _mm_mul_sd(a0_5, b0));
#endif
    _mm_store_sd(&C[(l_n*56)+55], c0_5);
#else
    C[(l_n*56)+0] += A[0] * B[(l_n*56)+0];
    C[(l_n*56)+3] += A[1] * B[(l_n*56)+0];
    C[(l_n*56)+9] += A[2] * B[(l_n*56)+0];
    C[(l_n*56)+19] += A[3] * B[(l_n*56)+0];
    C[(l_n*56)+34] += A[4] * B[(l_n*56)+0];
    C[(l_n*56)+55] += A[5] * B[(l_n*56)+0];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b1 = _mm256_broadcast_sd(&B[(l_n*56)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b1 = _mm_loaddup_pd(&B[(l_n*56)+1]);
#endif
    __m128d c1_0 = _mm_load_sd(&C[(l_n*56)+1]);
    __m128d a1_0 = _mm_load_sd(&A[6]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, b1));
#endif
    _mm_store_sd(&C[(l_n*56)+1], c1_0);
    __m128d c1_1 = _mm_load_sd(&C[(l_n*56)+7]);
    __m128d a1_1 = _mm_load_sd(&A[7]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_1 = _mm_add_sd(c1_1, _mm_mul_sd(a1_1, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_1 = _mm_add_sd(c1_1, _mm_mul_sd(a1_1, b1));
#endif
    _mm_store_sd(&C[(l_n*56)+7], c1_1);
    __m128d c1_2 = _mm_load_sd(&C[(l_n*56)+17]);
    __m128d a1_2 = _mm_load_sd(&A[8]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_2 = _mm_add_sd(c1_2, _mm_mul_sd(a1_2, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_2 = _mm_add_sd(c1_2, _mm_mul_sd(a1_2, b1));
#endif
    _mm_store_sd(&C[(l_n*56)+17], c1_2);
    __m128d c1_3 = _mm_load_sd(&C[(l_n*56)+32]);
    __m128d a1_3 = _mm_load_sd(&A[9]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_3 = _mm_add_sd(c1_3, _mm_mul_sd(a1_3, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_3 = _mm_add_sd(c1_3, _mm_mul_sd(a1_3, b1));
#endif
    _mm_store_sd(&C[(l_n*56)+32], c1_3);
    __m128d c1_4 = _mm_load_sd(&C[(l_n*56)+53]);
    __m128d a1_4 = _mm_load_sd(&A[10]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_4 = _mm_add_sd(c1_4, _mm_mul_sd(a1_4, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_4 = _mm_add_sd(c1_4, _mm_mul_sd(a1_4, b1));
#endif
    _mm_store_sd(&C[(l_n*56)+53], c1_4);
#else
    C[(l_n*56)+1] += A[6] * B[(l_n*56)+1];
    C[(l_n*56)+7] += A[7] * B[(l_n*56)+1];
    C[(l_n*56)+17] += A[8] * B[(l_n*56)+1];
    C[(l_n*56)+32] += A[9] * B[(l_n*56)+1];
    C[(l_n*56)+53] += A[10] * B[(l_n*56)+1];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b2 = _mm256_broadcast_sd(&B[(l_n*56)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b2 = _mm_loaddup_pd(&B[(l_n*56)+2]);
#endif
    __m128d c2_0 = _mm_load_sd(&C[(l_n*56)+2]);
    __m128d a2_0 = _mm_load_sd(&A[11]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_0 = _mm_add_sd(c2_0, _mm_mul_sd(a2_0, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_0 = _mm_add_sd(c2_0, _mm_mul_sd(a2_0, b2));
#endif
    _mm_store_sd(&C[(l_n*56)+2], c2_0);
    __m128d c2_1 = _mm_load_sd(&C[(l_n*56)+8]);
    __m128d a2_1 = _mm_load_sd(&A[12]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_1 = _mm_add_sd(c2_1, _mm_mul_sd(a2_1, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_1 = _mm_add_sd(c2_1, _mm_mul_sd(a2_1, b2));
#endif
    _mm_store_sd(&C[(l_n*56)+8], c2_1);
    __m128d c2_2 = _mm_load_sd(&C[(l_n*56)+18]);
    __m128d a2_2 = _mm_load_sd(&A[13]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_2 = _mm_add_sd(c2_2, _mm_mul_sd(a2_2, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_2 = _mm_add_sd(c2_2, _mm_mul_sd(a2_2, b2));
#endif
    _mm_store_sd(&C[(l_n*56)+18], c2_2);
    __m128d c2_3 = _mm_load_sd(&C[(l_n*56)+33]);
    __m128d a2_3 = _mm_load_sd(&A[14]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_3 = _mm_add_sd(c2_3, _mm_mul_sd(a2_3, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_3 = _mm_add_sd(c2_3, _mm_mul_sd(a2_3, b2));
#endif
    _mm_store_sd(&C[(l_n*56)+33], c2_3);
    __m128d c2_4 = _mm_load_sd(&C[(l_n*56)+54]);
    __m128d a2_4 = _mm_load_sd(&A[15]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_4 = _mm_add_sd(c2_4, _mm_mul_sd(a2_4, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_4 = _mm_add_sd(c2_4, _mm_mul_sd(a2_4, b2));
#endif
    _mm_store_sd(&C[(l_n*56)+54], c2_4);
#else
    C[(l_n*56)+2] += A[11] * B[(l_n*56)+2];
    C[(l_n*56)+8] += A[12] * B[(l_n*56)+2];
    C[(l_n*56)+18] += A[13] * B[(l_n*56)+2];
    C[(l_n*56)+33] += A[14] * B[(l_n*56)+2];
    C[(l_n*56)+54] += A[15] * B[(l_n*56)+2];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b3 = _mm256_broadcast_sd(&B[(l_n*56)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b3 = _mm_loaddup_pd(&B[(l_n*56)+3]);
#endif
    __m128d c3_0 = _mm_load_sd(&C[(l_n*56)+0]);
    __m128d a3_0 = _mm_load_sd(&A[16]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, b3));
#endif
    _mm_store_sd(&C[(l_n*56)+0], c3_0);
    __m128d c3_1 = _mm_load_sd(&C[(l_n*56)+3]);
    __m128d a3_1 = _mm_load_sd(&A[17]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, b3));
#endif
    _mm_store_sd(&C[(l_n*56)+3], c3_1);
    __m128d c3_2 = _mm_load_sd(&C[(l_n*56)+9]);
    __m128d a3_2 = _mm_load_sd(&A[18]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, b3));
#endif
    _mm_store_sd(&C[(l_n*56)+9], c3_2);
    __m128d c3_3 = _mm_load_sd(&C[(l_n*56)+19]);
    __m128d a3_3 = _mm_load_sd(&A[19]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_3 = _mm_add_sd(c3_3, _mm_mul_sd(a3_3, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_3 = _mm_add_sd(c3_3, _mm_mul_sd(a3_3, b3));
#endif
    _mm_store_sd(&C[(l_n*56)+19], c3_3);
    __m128d c3_4 = _mm_load_sd(&C[(l_n*56)+34]);
    __m128d a3_4 = _mm_load_sd(&A[20]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_4 = _mm_add_sd(c3_4, _mm_mul_sd(a3_4, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_4 = _mm_add_sd(c3_4, _mm_mul_sd(a3_4, b3));
#endif
    _mm_store_sd(&C[(l_n*56)+34], c3_4);
    __m128d c3_5 = _mm_load_sd(&C[(l_n*56)+55]);
    __m128d a3_5 = _mm_load_sd(&A[21]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_5 = _mm_add_sd(c3_5, _mm_mul_sd(a3_5, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_5 = _mm_add_sd(c3_5, _mm_mul_sd(a3_5, b3));
#endif
    _mm_store_sd(&C[(l_n*56)+55], c3_5);
#else
    C[(l_n*56)+0] += A[16] * B[(l_n*56)+3];
    C[(l_n*56)+3] += A[17] * B[(l_n*56)+3];
    C[(l_n*56)+9] += A[18] * B[(l_n*56)+3];
    C[(l_n*56)+19] += A[19] * B[(l_n*56)+3];
    C[(l_n*56)+34] += A[20] * B[(l_n*56)+3];
    C[(l_n*56)+55] += A[21] * B[(l_n*56)+3];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b4 = _mm256_broadcast_sd(&B[(l_n*56)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b4 = _mm_loaddup_pd(&B[(l_n*56)+4]);
#endif
    __m128d c4_0 = _mm_load_sd(&C[(l_n*56)+4]);
    __m128d a4_0 = _mm_load_sd(&A[22]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_0 = _mm_add_sd(c4_0, _mm_mul_sd(a4_0, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_0 = _mm_add_sd(c4_0, _mm_mul_sd(a4_0, b4));
#endif
    _mm_store_sd(&C[(l_n*56)+4], c4_0);
    __m128d c4_1 = _mm_load_sd(&C[(l_n*56)+14]);
    __m128d a4_1 = _mm_load_sd(&A[23]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_1 = _mm_add_sd(c4_1, _mm_mul_sd(a4_1, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_1 = _mm_add_sd(c4_1, _mm_mul_sd(a4_1, b4));
#endif
    _mm_store_sd(&C[(l_n*56)+14], c4_1);
    __m128d c4_2 = _mm_load_sd(&C[(l_n*56)+29]);
    __m128d a4_2 = _mm_load_sd(&A[24]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, b4));
#endif
    _mm_store_sd(&C[(l_n*56)+29], c4_2);
    __m128d c4_3 = _mm_load_sd(&C[(l_n*56)+50]);
    __m128d a4_3 = _mm_load_sd(&A[25]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_3 = _mm_add_sd(c4_3, _mm_mul_sd(a4_3, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_3 = _mm_add_sd(c4_3, _mm_mul_sd(a4_3, b4));
#endif
    _mm_store_sd(&C[(l_n*56)+50], c4_3);
#else
    C[(l_n*56)+4] += A[22] * B[(l_n*56)+4];
    C[(l_n*56)+14] += A[23] * B[(l_n*56)+4];
    C[(l_n*56)+29] += A[24] * B[(l_n*56)+4];
    C[(l_n*56)+50] += A[25] * B[(l_n*56)+4];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b5 = _mm256_broadcast_sd(&B[(l_n*56)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b5 = _mm_loaddup_pd(&B[(l_n*56)+5]);
#endif
    __m128d c5_0 = _mm_load_sd(&C[(l_n*56)+5]);
    __m128d a5_0 = _mm_load_sd(&A[26]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_0 = _mm_add_sd(c5_0, _mm_mul_sd(a5_0, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_0 = _mm_add_sd(c5_0, _mm_mul_sd(a5_0, b5));
#endif
    _mm_store_sd(&C[(l_n*56)+5], c5_0);
    __m128d c5_1 = _mm_load_sd(&C[(l_n*56)+15]);
    __m128d a5_1 = _mm_load_sd(&A[27]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_1 = _mm_add_sd(c5_1, _mm_mul_sd(a5_1, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_1 = _mm_add_sd(c5_1, _mm_mul_sd(a5_1, b5));
#endif
    _mm_store_sd(&C[(l_n*56)+15], c5_1);
    __m128d c5_2 = _mm_load_sd(&C[(l_n*56)+30]);
    __m128d a5_2 = _mm_load_sd(&A[28]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, b5));
#endif
    _mm_store_sd(&C[(l_n*56)+30], c5_2);
    __m128d c5_3 = _mm_load_sd(&C[(l_n*56)+51]);
    __m128d a5_3 = _mm_load_sd(&A[29]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_3 = _mm_add_sd(c5_3, _mm_mul_sd(a5_3, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_3 = _mm_add_sd(c5_3, _mm_mul_sd(a5_3, b5));
#endif
    _mm_store_sd(&C[(l_n*56)+51], c5_3);
#else
    C[(l_n*56)+5] += A[26] * B[(l_n*56)+5];
    C[(l_n*56)+15] += A[27] * B[(l_n*56)+5];
    C[(l_n*56)+30] += A[28] * B[(l_n*56)+5];
    C[(l_n*56)+51] += A[29] * B[(l_n*56)+5];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b6 = _mm256_broadcast_sd(&B[(l_n*56)+6]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b6 = _mm_loaddup_pd(&B[(l_n*56)+6]);
#endif
    __m128d c6_0 = _mm_load_sd(&C[(l_n*56)+6]);
    __m128d a6_0 = _mm_load_sd(&A[30]);
#if defined(__SSE3__) && defined(__AVX__)
    c6_0 = _mm_add_sd(c6_0, _mm_mul_sd(a6_0, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c6_0 = _mm_add_sd(c6_0, _mm_mul_sd(a6_0, b6));
#endif
    _mm_store_sd(&C[(l_n*56)+6], c6_0);
    __m128d c6_1 = _mm_load_sd(&C[(l_n*56)+16]);
    __m128d a6_1 = _mm_load_sd(&A[31]);
#if defined(__SSE3__) && defined(__AVX__)
    c6_1 = _mm_add_sd(c6_1, _mm_mul_sd(a6_1, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c6_1 = _mm_add_sd(c6_1, _mm_mul_sd(a6_1, b6));
#endif
    _mm_store_sd(&C[(l_n*56)+16], c6_1);
    __m128d c6_2 = _mm_load_sd(&C[(l_n*56)+31]);
    __m128d a6_2 = _mm_load_sd(&A[32]);
#if defined(__SSE3__) && defined(__AVX__)
    c6_2 = _mm_add_sd(c6_2, _mm_mul_sd(a6_2, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c6_2 = _mm_add_sd(c6_2, _mm_mul_sd(a6_2, b6));
#endif
    _mm_store_sd(&C[(l_n*56)+31], c6_2);
    __m128d c6_3 = _mm_load_sd(&C[(l_n*56)+52]);
    __m128d a6_3 = _mm_load_sd(&A[33]);
#if defined(__SSE3__) && defined(__AVX__)
    c6_3 = _mm_add_sd(c6_3, _mm_mul_sd(a6_3, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c6_3 = _mm_add_sd(c6_3, _mm_mul_sd(a6_3, b6));
#endif
    _mm_store_sd(&C[(l_n*56)+52], c6_3);
#else
    C[(l_n*56)+6] += A[30] * B[(l_n*56)+6];
    C[(l_n*56)+16] += A[31] * B[(l_n*56)+6];
    C[(l_n*56)+31] += A[32] * B[(l_n*56)+6];
    C[(l_n*56)+52] += A[33] * B[(l_n*56)+6];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b7 = _mm256_broadcast_sd(&B[(l_n*56)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b7 = _mm_loaddup_pd(&B[(l_n*56)+7]);
#endif
    __m128d c7_0 = _mm_load_sd(&C[(l_n*56)+1]);
    __m128d a7_0 = _mm_load_sd(&A[34]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_0 = _mm_add_sd(c7_0, _mm_mul_sd(a7_0, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_0 = _mm_add_sd(c7_0, _mm_mul_sd(a7_0, b7));
#endif
    _mm_store_sd(&C[(l_n*56)+1], c7_0);
    __m128d c7_1 = _mm_load_sd(&C[(l_n*56)+7]);
    __m128d a7_1 = _mm_load_sd(&A[35]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_1 = _mm_add_sd(c7_1, _mm_mul_sd(a7_1, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_1 = _mm_add_sd(c7_1, _mm_mul_sd(a7_1, b7));
#endif
    _mm_store_sd(&C[(l_n*56)+7], c7_1);
    __m128d c7_2 = _mm_load_sd(&C[(l_n*56)+17]);
    __m128d a7_2 = _mm_load_sd(&A[36]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_2 = _mm_add_sd(c7_2, _mm_mul_sd(a7_2, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_2 = _mm_add_sd(c7_2, _mm_mul_sd(a7_2, b7));
#endif
    _mm_store_sd(&C[(l_n*56)+17], c7_2);
    __m128d c7_3 = _mm_load_sd(&C[(l_n*56)+32]);
    __m128d a7_3 = _mm_load_sd(&A[37]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_3 = _mm_add_sd(c7_3, _mm_mul_sd(a7_3, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_3 = _mm_add_sd(c7_3, _mm_mul_sd(a7_3, b7));
#endif
    _mm_store_sd(&C[(l_n*56)+32], c7_3);
    __m128d c7_4 = _mm_load_sd(&C[(l_n*56)+53]);
    __m128d a7_4 = _mm_load_sd(&A[38]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_4 = _mm_add_sd(c7_4, _mm_mul_sd(a7_4, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_4 = _mm_add_sd(c7_4, _mm_mul_sd(a7_4, b7));
#endif
    _mm_store_sd(&C[(l_n*56)+53], c7_4);
#else
    C[(l_n*56)+1] += A[34] * B[(l_n*56)+7];
    C[(l_n*56)+7] += A[35] * B[(l_n*56)+7];
    C[(l_n*56)+17] += A[36] * B[(l_n*56)+7];
    C[(l_n*56)+32] += A[37] * B[(l_n*56)+7];
    C[(l_n*56)+53] += A[38] * B[(l_n*56)+7];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b8 = _mm256_broadcast_sd(&B[(l_n*56)+8]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b8 = _mm_loaddup_pd(&B[(l_n*56)+8]);
#endif
    __m128d c8_0 = _mm_load_sd(&C[(l_n*56)+2]);
    __m128d a8_0 = _mm_load_sd(&A[39]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_0 = _mm_add_sd(c8_0, _mm_mul_sd(a8_0, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_0 = _mm_add_sd(c8_0, _mm_mul_sd(a8_0, b8));
#endif
    _mm_store_sd(&C[(l_n*56)+2], c8_0);
    __m128d c8_1 = _mm_load_sd(&C[(l_n*56)+8]);
    __m128d a8_1 = _mm_load_sd(&A[40]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_1 = _mm_add_sd(c8_1, _mm_mul_sd(a8_1, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_1 = _mm_add_sd(c8_1, _mm_mul_sd(a8_1, b8));
#endif
    _mm_store_sd(&C[(l_n*56)+8], c8_1);
    __m128d c8_2 = _mm_load_sd(&C[(l_n*56)+18]);
    __m128d a8_2 = _mm_load_sd(&A[41]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_2 = _mm_add_sd(c8_2, _mm_mul_sd(a8_2, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_2 = _mm_add_sd(c8_2, _mm_mul_sd(a8_2, b8));
#endif
    _mm_store_sd(&C[(l_n*56)+18], c8_2);
    __m128d c8_3 = _mm_load_sd(&C[(l_n*56)+33]);
    __m128d a8_3 = _mm_load_sd(&A[42]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_3 = _mm_add_sd(c8_3, _mm_mul_sd(a8_3, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_3 = _mm_add_sd(c8_3, _mm_mul_sd(a8_3, b8));
#endif
    _mm_store_sd(&C[(l_n*56)+33], c8_3);
    __m128d c8_4 = _mm_load_sd(&C[(l_n*56)+54]);
    __m128d a8_4 = _mm_load_sd(&A[43]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_4 = _mm_add_sd(c8_4, _mm_mul_sd(a8_4, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_4 = _mm_add_sd(c8_4, _mm_mul_sd(a8_4, b8));
#endif
    _mm_store_sd(&C[(l_n*56)+54], c8_4);
#else
    C[(l_n*56)+2] += A[39] * B[(l_n*56)+8];
    C[(l_n*56)+8] += A[40] * B[(l_n*56)+8];
    C[(l_n*56)+18] += A[41] * B[(l_n*56)+8];
    C[(l_n*56)+33] += A[42] * B[(l_n*56)+8];
    C[(l_n*56)+54] += A[43] * B[(l_n*56)+8];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b9 = _mm256_broadcast_sd(&B[(l_n*56)+9]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b9 = _mm_loaddup_pd(&B[(l_n*56)+9]);
#endif
    __m128d c9_0 = _mm_load_sd(&C[(l_n*56)+0]);
    __m128d a9_0 = _mm_load_sd(&A[44]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, b9));
#endif
    _mm_store_sd(&C[(l_n*56)+0], c9_0);
    __m128d c9_1 = _mm_load_sd(&C[(l_n*56)+3]);
    __m128d a9_1 = _mm_load_sd(&A[45]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, b9));
#endif
    _mm_store_sd(&C[(l_n*56)+3], c9_1);
    __m128d c9_2 = _mm_load_sd(&C[(l_n*56)+9]);
    __m128d a9_2 = _mm_load_sd(&A[46]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, b9));
#endif
    _mm_store_sd(&C[(l_n*56)+9], c9_2);
    __m128d c9_3 = _mm_load_sd(&C[(l_n*56)+19]);
    __m128d a9_3 = _mm_load_sd(&A[47]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_3 = _mm_add_sd(c9_3, _mm_mul_sd(a9_3, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_3 = _mm_add_sd(c9_3, _mm_mul_sd(a9_3, b9));
#endif
    _mm_store_sd(&C[(l_n*56)+19], c9_3);
    __m128d c9_4 = _mm_load_sd(&C[(l_n*56)+34]);
    __m128d a9_4 = _mm_load_sd(&A[48]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_4 = _mm_add_sd(c9_4, _mm_mul_sd(a9_4, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_4 = _mm_add_sd(c9_4, _mm_mul_sd(a9_4, b9));
#endif
    _mm_store_sd(&C[(l_n*56)+34], c9_4);
    __m128d c9_5 = _mm_load_sd(&C[(l_n*56)+55]);
    __m128d a9_5 = _mm_load_sd(&A[49]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_5 = _mm_add_sd(c9_5, _mm_mul_sd(a9_5, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_5 = _mm_add_sd(c9_5, _mm_mul_sd(a9_5, b9));
#endif
    _mm_store_sd(&C[(l_n*56)+55], c9_5);
#else
    C[(l_n*56)+0] += A[44] * B[(l_n*56)+9];
    C[(l_n*56)+3] += A[45] * B[(l_n*56)+9];
    C[(l_n*56)+9] += A[46] * B[(l_n*56)+9];
    C[(l_n*56)+19] += A[47] * B[(l_n*56)+9];
    C[(l_n*56)+34] += A[48] * B[(l_n*56)+9];
    C[(l_n*56)+55] += A[49] * B[(l_n*56)+9];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b10 = _mm256_broadcast_sd(&B[(l_n*56)+10]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b10 = _mm_loaddup_pd(&B[(l_n*56)+10]);
#endif
    __m128d c10_0 = _mm_load_sd(&C[(l_n*56)+10]);
    __m128d a10_0 = _mm_load_sd(&A[50]);
#if defined(__SSE3__) && defined(__AVX__)
    c10_0 = _mm_add_sd(c10_0, _mm_mul_sd(a10_0, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c10_0 = _mm_add_sd(c10_0, _mm_mul_sd(a10_0, b10));
#endif
    _mm_store_sd(&C[(l_n*56)+10], c10_0);
    __m128d c10_1 = _mm_load_sd(&C[(l_n*56)+25]);
    __m128d a10_1 = _mm_load_sd(&A[51]);
#if defined(__SSE3__) && defined(__AVX__)
    c10_1 = _mm_add_sd(c10_1, _mm_mul_sd(a10_1, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c10_1 = _mm_add_sd(c10_1, _mm_mul_sd(a10_1, b10));
#endif
    _mm_store_sd(&C[(l_n*56)+25], c10_1);
    __m128d c10_2 = _mm_load_sd(&C[(l_n*56)+46]);
    __m128d a10_2 = _mm_load_sd(&A[52]);
#if defined(__SSE3__) && defined(__AVX__)
    c10_2 = _mm_add_sd(c10_2, _mm_mul_sd(a10_2, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c10_2 = _mm_add_sd(c10_2, _mm_mul_sd(a10_2, b10));
#endif
    _mm_store_sd(&C[(l_n*56)+46], c10_2);
#else
    C[(l_n*56)+10] += A[50] * B[(l_n*56)+10];
    C[(l_n*56)+25] += A[51] * B[(l_n*56)+10];
    C[(l_n*56)+46] += A[52] * B[(l_n*56)+10];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b11 = _mm256_broadcast_sd(&B[(l_n*56)+11]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b11 = _mm_loaddup_pd(&B[(l_n*56)+11]);
#endif
    __m128d c11_0 = _mm_load_sd(&C[(l_n*56)+11]);
    __m128d a11_0 = _mm_load_sd(&A[53]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_0 = _mm_add_sd(c11_0, _mm_mul_sd(a11_0, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_0 = _mm_add_sd(c11_0, _mm_mul_sd(a11_0, b11));
#endif
    _mm_store_sd(&C[(l_n*56)+11], c11_0);
    __m128d c11_1 = _mm_load_sd(&C[(l_n*56)+26]);
    __m128d a11_1 = _mm_load_sd(&A[54]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_1 = _mm_add_sd(c11_1, _mm_mul_sd(a11_1, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_1 = _mm_add_sd(c11_1, _mm_mul_sd(a11_1, b11));
#endif
    _mm_store_sd(&C[(l_n*56)+26], c11_1);
    __m128d c11_2 = _mm_load_sd(&C[(l_n*56)+47]);
    __m128d a11_2 = _mm_load_sd(&A[55]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_2 = _mm_add_sd(c11_2, _mm_mul_sd(a11_2, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_2 = _mm_add_sd(c11_2, _mm_mul_sd(a11_2, b11));
#endif
    _mm_store_sd(&C[(l_n*56)+47], c11_2);
#else
    C[(l_n*56)+11] += A[53] * B[(l_n*56)+11];
    C[(l_n*56)+26] += A[54] * B[(l_n*56)+11];
    C[(l_n*56)+47] += A[55] * B[(l_n*56)+11];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b12 = _mm256_broadcast_sd(&B[(l_n*56)+12]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b12 = _mm_loaddup_pd(&B[(l_n*56)+12]);
#endif
    __m128d c12_0 = _mm_load_sd(&C[(l_n*56)+12]);
    __m128d a12_0 = _mm_load_sd(&A[56]);
#if defined(__SSE3__) && defined(__AVX__)
    c12_0 = _mm_add_sd(c12_0, _mm_mul_sd(a12_0, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c12_0 = _mm_add_sd(c12_0, _mm_mul_sd(a12_0, b12));
#endif
    _mm_store_sd(&C[(l_n*56)+12], c12_0);
    __m128d c12_1 = _mm_load_sd(&C[(l_n*56)+27]);
    __m128d a12_1 = _mm_load_sd(&A[57]);
#if defined(__SSE3__) && defined(__AVX__)
    c12_1 = _mm_add_sd(c12_1, _mm_mul_sd(a12_1, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c12_1 = _mm_add_sd(c12_1, _mm_mul_sd(a12_1, b12));
#endif
    _mm_store_sd(&C[(l_n*56)+27], c12_1);
    __m128d c12_2 = _mm_load_sd(&C[(l_n*56)+48]);
    __m128d a12_2 = _mm_load_sd(&A[58]);
#if defined(__SSE3__) && defined(__AVX__)
    c12_2 = _mm_add_sd(c12_2, _mm_mul_sd(a12_2, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c12_2 = _mm_add_sd(c12_2, _mm_mul_sd(a12_2, b12));
#endif
    _mm_store_sd(&C[(l_n*56)+48], c12_2);
#else
    C[(l_n*56)+12] += A[56] * B[(l_n*56)+12];
    C[(l_n*56)+27] += A[57] * B[(l_n*56)+12];
    C[(l_n*56)+48] += A[58] * B[(l_n*56)+12];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b13 = _mm256_broadcast_sd(&B[(l_n*56)+13]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b13 = _mm_loaddup_pd(&B[(l_n*56)+13]);
#endif
    __m128d c13_0 = _mm_load_sd(&C[(l_n*56)+13]);
    __m128d a13_0 = _mm_load_sd(&A[59]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_0 = _mm_add_sd(c13_0, _mm_mul_sd(a13_0, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_0 = _mm_add_sd(c13_0, _mm_mul_sd(a13_0, b13));
#endif
    _mm_store_sd(&C[(l_n*56)+13], c13_0);
    __m128d c13_1 = _mm_load_sd(&C[(l_n*56)+28]);
    __m128d a13_1 = _mm_load_sd(&A[60]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_1 = _mm_add_sd(c13_1, _mm_mul_sd(a13_1, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_1 = _mm_add_sd(c13_1, _mm_mul_sd(a13_1, b13));
#endif
    _mm_store_sd(&C[(l_n*56)+28], c13_1);
    __m128d c13_2 = _mm_load_sd(&C[(l_n*56)+49]);
    __m128d a13_2 = _mm_load_sd(&A[61]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_2 = _mm_add_sd(c13_2, _mm_mul_sd(a13_2, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_2 = _mm_add_sd(c13_2, _mm_mul_sd(a13_2, b13));
#endif
    _mm_store_sd(&C[(l_n*56)+49], c13_2);
#else
    C[(l_n*56)+13] += A[59] * B[(l_n*56)+13];
    C[(l_n*56)+28] += A[60] * B[(l_n*56)+13];
    C[(l_n*56)+49] += A[61] * B[(l_n*56)+13];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b14 = _mm256_broadcast_sd(&B[(l_n*56)+14]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b14 = _mm_loaddup_pd(&B[(l_n*56)+14]);
#endif
    __m128d c14_0 = _mm_load_sd(&C[(l_n*56)+4]);
    __m128d a14_0 = _mm_load_sd(&A[62]);
#if defined(__SSE3__) && defined(__AVX__)
    c14_0 = _mm_add_sd(c14_0, _mm_mul_sd(a14_0, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c14_0 = _mm_add_sd(c14_0, _mm_mul_sd(a14_0, b14));
#endif
    _mm_store_sd(&C[(l_n*56)+4], c14_0);
    __m128d c14_1 = _mm_load_sd(&C[(l_n*56)+14]);
    __m128d a14_1 = _mm_load_sd(&A[63]);
#if defined(__SSE3__) && defined(__AVX__)
    c14_1 = _mm_add_sd(c14_1, _mm_mul_sd(a14_1, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c14_1 = _mm_add_sd(c14_1, _mm_mul_sd(a14_1, b14));
#endif
    _mm_store_sd(&C[(l_n*56)+14], c14_1);
    __m128d c14_2 = _mm_load_sd(&C[(l_n*56)+29]);
    __m128d a14_2 = _mm_load_sd(&A[64]);
#if defined(__SSE3__) && defined(__AVX__)
    c14_2 = _mm_add_sd(c14_2, _mm_mul_sd(a14_2, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c14_2 = _mm_add_sd(c14_2, _mm_mul_sd(a14_2, b14));
#endif
    _mm_store_sd(&C[(l_n*56)+29], c14_2);
    __m128d c14_3 = _mm_load_sd(&C[(l_n*56)+50]);
    __m128d a14_3 = _mm_load_sd(&A[65]);
#if defined(__SSE3__) && defined(__AVX__)
    c14_3 = _mm_add_sd(c14_3, _mm_mul_sd(a14_3, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c14_3 = _mm_add_sd(c14_3, _mm_mul_sd(a14_3, b14));
#endif
    _mm_store_sd(&C[(l_n*56)+50], c14_3);
#else
    C[(l_n*56)+4] += A[62] * B[(l_n*56)+14];
    C[(l_n*56)+14] += A[63] * B[(l_n*56)+14];
    C[(l_n*56)+29] += A[64] * B[(l_n*56)+14];
    C[(l_n*56)+50] += A[65] * B[(l_n*56)+14];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b15 = _mm256_broadcast_sd(&B[(l_n*56)+15]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b15 = _mm_loaddup_pd(&B[(l_n*56)+15]);
#endif
    __m128d c15_0 = _mm_load_sd(&C[(l_n*56)+5]);
    __m128d a15_0 = _mm_load_sd(&A[66]);
#if defined(__SSE3__) && defined(__AVX__)
    c15_0 = _mm_add_sd(c15_0, _mm_mul_sd(a15_0, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c15_0 = _mm_add_sd(c15_0, _mm_mul_sd(a15_0, b15));
#endif
    _mm_store_sd(&C[(l_n*56)+5], c15_0);
    __m128d c15_1 = _mm_load_sd(&C[(l_n*56)+15]);
    __m128d a15_1 = _mm_load_sd(&A[67]);
#if defined(__SSE3__) && defined(__AVX__)
    c15_1 = _mm_add_sd(c15_1, _mm_mul_sd(a15_1, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c15_1 = _mm_add_sd(c15_1, _mm_mul_sd(a15_1, b15));
#endif
    _mm_store_sd(&C[(l_n*56)+15], c15_1);
    __m128d c15_2 = _mm_load_sd(&C[(l_n*56)+30]);
    __m128d a15_2 = _mm_load_sd(&A[68]);
#if defined(__SSE3__) && defined(__AVX__)
    c15_2 = _mm_add_sd(c15_2, _mm_mul_sd(a15_2, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c15_2 = _mm_add_sd(c15_2, _mm_mul_sd(a15_2, b15));
#endif
    _mm_store_sd(&C[(l_n*56)+30], c15_2);
    __m128d c15_3 = _mm_load_sd(&C[(l_n*56)+51]);
    __m128d a15_3 = _mm_load_sd(&A[69]);
#if defined(__SSE3__) && defined(__AVX__)
    c15_3 = _mm_add_sd(c15_3, _mm_mul_sd(a15_3, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c15_3 = _mm_add_sd(c15_3, _mm_mul_sd(a15_3, b15));
#endif
    _mm_store_sd(&C[(l_n*56)+51], c15_3);
#else
    C[(l_n*56)+5] += A[66] * B[(l_n*56)+15];
    C[(l_n*56)+15] += A[67] * B[(l_n*56)+15];
    C[(l_n*56)+30] += A[68] * B[(l_n*56)+15];
    C[(l_n*56)+51] += A[69] * B[(l_n*56)+15];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b16 = _mm256_broadcast_sd(&B[(l_n*56)+16]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b16 = _mm_loaddup_pd(&B[(l_n*56)+16]);
#endif
    __m128d c16_0 = _mm_load_sd(&C[(l_n*56)+6]);
    __m128d a16_0 = _mm_load_sd(&A[70]);
#if defined(__SSE3__) && defined(__AVX__)
    c16_0 = _mm_add_sd(c16_0, _mm_mul_sd(a16_0, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c16_0 = _mm_add_sd(c16_0, _mm_mul_sd(a16_0, b16));
#endif
    _mm_store_sd(&C[(l_n*56)+6], c16_0);
    __m128d c16_1 = _mm_load_sd(&C[(l_n*56)+16]);
    __m128d a16_1 = _mm_load_sd(&A[71]);
#if defined(__SSE3__) && defined(__AVX__)
    c16_1 = _mm_add_sd(c16_1, _mm_mul_sd(a16_1, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c16_1 = _mm_add_sd(c16_1, _mm_mul_sd(a16_1, b16));
#endif
    _mm_store_sd(&C[(l_n*56)+16], c16_1);
    __m128d c16_2 = _mm_load_sd(&C[(l_n*56)+31]);
    __m128d a16_2 = _mm_load_sd(&A[72]);
#if defined(__SSE3__) && defined(__AVX__)
    c16_2 = _mm_add_sd(c16_2, _mm_mul_sd(a16_2, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c16_2 = _mm_add_sd(c16_2, _mm_mul_sd(a16_2, b16));
#endif
    _mm_store_sd(&C[(l_n*56)+31], c16_2);
    __m128d c16_3 = _mm_load_sd(&C[(l_n*56)+52]);
    __m128d a16_3 = _mm_load_sd(&A[73]);
#if defined(__SSE3__) && defined(__AVX__)
    c16_3 = _mm_add_sd(c16_3, _mm_mul_sd(a16_3, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c16_3 = _mm_add_sd(c16_3, _mm_mul_sd(a16_3, b16));
#endif
    _mm_store_sd(&C[(l_n*56)+52], c16_3);
#else
    C[(l_n*56)+6] += A[70] * B[(l_n*56)+16];
    C[(l_n*56)+16] += A[71] * B[(l_n*56)+16];
    C[(l_n*56)+31] += A[72] * B[(l_n*56)+16];
    C[(l_n*56)+52] += A[73] * B[(l_n*56)+16];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b17 = _mm256_broadcast_sd(&B[(l_n*56)+17]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b17 = _mm_loaddup_pd(&B[(l_n*56)+17]);
#endif
    __m128d c17_0 = _mm_load_sd(&C[(l_n*56)+1]);
    __m128d a17_0 = _mm_load_sd(&A[74]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_0 = _mm_add_sd(c17_0, _mm_mul_sd(a17_0, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_0 = _mm_add_sd(c17_0, _mm_mul_sd(a17_0, b17));
#endif
    _mm_store_sd(&C[(l_n*56)+1], c17_0);
    __m128d c17_1 = _mm_load_sd(&C[(l_n*56)+7]);
    __m128d a17_1 = _mm_load_sd(&A[75]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_1 = _mm_add_sd(c17_1, _mm_mul_sd(a17_1, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_1 = _mm_add_sd(c17_1, _mm_mul_sd(a17_1, b17));
#endif
    _mm_store_sd(&C[(l_n*56)+7], c17_1);
    __m128d c17_2 = _mm_load_sd(&C[(l_n*56)+17]);
    __m128d a17_2 = _mm_load_sd(&A[76]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_2 = _mm_add_sd(c17_2, _mm_mul_sd(a17_2, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_2 = _mm_add_sd(c17_2, _mm_mul_sd(a17_2, b17));
#endif
    _mm_store_sd(&C[(l_n*56)+17], c17_2);
    __m128d c17_3 = _mm_load_sd(&C[(l_n*56)+32]);
    __m128d a17_3 = _mm_load_sd(&A[77]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_3 = _mm_add_sd(c17_3, _mm_mul_sd(a17_3, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_3 = _mm_add_sd(c17_3, _mm_mul_sd(a17_3, b17));
#endif
    _mm_store_sd(&C[(l_n*56)+32], c17_3);
    __m128d c17_4 = _mm_load_sd(&C[(l_n*56)+53]);
    __m128d a17_4 = _mm_load_sd(&A[78]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_4 = _mm_add_sd(c17_4, _mm_mul_sd(a17_4, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_4 = _mm_add_sd(c17_4, _mm_mul_sd(a17_4, b17));
#endif
    _mm_store_sd(&C[(l_n*56)+53], c17_4);
#else
    C[(l_n*56)+1] += A[74] * B[(l_n*56)+17];
    C[(l_n*56)+7] += A[75] * B[(l_n*56)+17];
    C[(l_n*56)+17] += A[76] * B[(l_n*56)+17];
    C[(l_n*56)+32] += A[77] * B[(l_n*56)+17];
    C[(l_n*56)+53] += A[78] * B[(l_n*56)+17];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b18 = _mm256_broadcast_sd(&B[(l_n*56)+18]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b18 = _mm_loaddup_pd(&B[(l_n*56)+18]);
#endif
    __m128d c18_0 = _mm_load_sd(&C[(l_n*56)+2]);
    __m128d a18_0 = _mm_load_sd(&A[79]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_0 = _mm_add_sd(c18_0, _mm_mul_sd(a18_0, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_0 = _mm_add_sd(c18_0, _mm_mul_sd(a18_0, b18));
#endif
    _mm_store_sd(&C[(l_n*56)+2], c18_0);
    __m128d c18_1 = _mm_load_sd(&C[(l_n*56)+8]);
    __m128d a18_1 = _mm_load_sd(&A[80]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_1 = _mm_add_sd(c18_1, _mm_mul_sd(a18_1, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_1 = _mm_add_sd(c18_1, _mm_mul_sd(a18_1, b18));
#endif
    _mm_store_sd(&C[(l_n*56)+8], c18_1);
    __m128d c18_2 = _mm_load_sd(&C[(l_n*56)+18]);
    __m128d a18_2 = _mm_load_sd(&A[81]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_2 = _mm_add_sd(c18_2, _mm_mul_sd(a18_2, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_2 = _mm_add_sd(c18_2, _mm_mul_sd(a18_2, b18));
#endif
    _mm_store_sd(&C[(l_n*56)+18], c18_2);
    __m128d c18_3 = _mm_load_sd(&C[(l_n*56)+33]);
    __m128d a18_3 = _mm_load_sd(&A[82]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_3 = _mm_add_sd(c18_3, _mm_mul_sd(a18_3, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_3 = _mm_add_sd(c18_3, _mm_mul_sd(a18_3, b18));
#endif
    _mm_store_sd(&C[(l_n*56)+33], c18_3);
    __m128d c18_4 = _mm_load_sd(&C[(l_n*56)+54]);
    __m128d a18_4 = _mm_load_sd(&A[83]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_4 = _mm_add_sd(c18_4, _mm_mul_sd(a18_4, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_4 = _mm_add_sd(c18_4, _mm_mul_sd(a18_4, b18));
#endif
    _mm_store_sd(&C[(l_n*56)+54], c18_4);
#else
    C[(l_n*56)+2] += A[79] * B[(l_n*56)+18];
    C[(l_n*56)+8] += A[80] * B[(l_n*56)+18];
    C[(l_n*56)+18] += A[81] * B[(l_n*56)+18];
    C[(l_n*56)+33] += A[82] * B[(l_n*56)+18];
    C[(l_n*56)+54] += A[83] * B[(l_n*56)+18];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b19 = _mm256_broadcast_sd(&B[(l_n*56)+19]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b19 = _mm_loaddup_pd(&B[(l_n*56)+19]);
#endif
    __m128d c19_0 = _mm_load_sd(&C[(l_n*56)+0]);
    __m128d a19_0 = _mm_load_sd(&A[84]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_0 = _mm_add_sd(c19_0, _mm_mul_sd(a19_0, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_0 = _mm_add_sd(c19_0, _mm_mul_sd(a19_0, b19));
#endif
    _mm_store_sd(&C[(l_n*56)+0], c19_0);
    __m128d c19_1 = _mm_load_sd(&C[(l_n*56)+3]);
    __m128d a19_1 = _mm_load_sd(&A[85]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_1 = _mm_add_sd(c19_1, _mm_mul_sd(a19_1, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_1 = _mm_add_sd(c19_1, _mm_mul_sd(a19_1, b19));
#endif
    _mm_store_sd(&C[(l_n*56)+3], c19_1);
    __m128d c19_2 = _mm_load_sd(&C[(l_n*56)+9]);
    __m128d a19_2 = _mm_load_sd(&A[86]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_2 = _mm_add_sd(c19_2, _mm_mul_sd(a19_2, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_2 = _mm_add_sd(c19_2, _mm_mul_sd(a19_2, b19));
#endif
    _mm_store_sd(&C[(l_n*56)+9], c19_2);
    __m128d c19_3 = _mm_load_sd(&C[(l_n*56)+19]);
    __m128d a19_3 = _mm_load_sd(&A[87]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_3 = _mm_add_sd(c19_3, _mm_mul_sd(a19_3, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_3 = _mm_add_sd(c19_3, _mm_mul_sd(a19_3, b19));
#endif
    _mm_store_sd(&C[(l_n*56)+19], c19_3);
    __m128d c19_4 = _mm_load_sd(&C[(l_n*56)+34]);
    __m128d a19_4 = _mm_load_sd(&A[88]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_4 = _mm_add_sd(c19_4, _mm_mul_sd(a19_4, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_4 = _mm_add_sd(c19_4, _mm_mul_sd(a19_4, b19));
#endif
    _mm_store_sd(&C[(l_n*56)+34], c19_4);
    __m128d c19_5 = _mm_load_sd(&C[(l_n*56)+55]);
    __m128d a19_5 = _mm_load_sd(&A[89]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_5 = _mm_add_sd(c19_5, _mm_mul_sd(a19_5, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_5 = _mm_add_sd(c19_5, _mm_mul_sd(a19_5, b19));
#endif
    _mm_store_sd(&C[(l_n*56)+55], c19_5);
#else
    C[(l_n*56)+0] += A[84] * B[(l_n*56)+19];
    C[(l_n*56)+3] += A[85] * B[(l_n*56)+19];
    C[(l_n*56)+9] += A[86] * B[(l_n*56)+19];
    C[(l_n*56)+19] += A[87] * B[(l_n*56)+19];
    C[(l_n*56)+34] += A[88] * B[(l_n*56)+19];
    C[(l_n*56)+55] += A[89] * B[(l_n*56)+19];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b20 = _mm256_broadcast_sd(&B[(l_n*56)+20]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b20 = _mm_loaddup_pd(&B[(l_n*56)+20]);
#endif
    __m128d c20_0 = _mm_load_sd(&C[(l_n*56)+20]);
    __m128d a20_0 = _mm_load_sd(&A[90]);
#if defined(__SSE3__) && defined(__AVX__)
    c20_0 = _mm_add_sd(c20_0, _mm_mul_sd(a20_0, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c20_0 = _mm_add_sd(c20_0, _mm_mul_sd(a20_0, b20));
#endif
    _mm_store_sd(&C[(l_n*56)+20], c20_0);
    __m128d c20_1 = _mm_load_sd(&C[(l_n*56)+41]);
    __m128d a20_1 = _mm_load_sd(&A[91]);
#if defined(__SSE3__) && defined(__AVX__)
    c20_1 = _mm_add_sd(c20_1, _mm_mul_sd(a20_1, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c20_1 = _mm_add_sd(c20_1, _mm_mul_sd(a20_1, b20));
#endif
    _mm_store_sd(&C[(l_n*56)+41], c20_1);
#else
    C[(l_n*56)+20] += A[90] * B[(l_n*56)+20];
    C[(l_n*56)+41] += A[91] * B[(l_n*56)+20];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b21 = _mm256_broadcast_sd(&B[(l_n*56)+21]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b21 = _mm_loaddup_pd(&B[(l_n*56)+21]);
#endif
    __m128d c21_0 = _mm_load_sd(&C[(l_n*56)+21]);
    __m128d a21_0 = _mm_load_sd(&A[92]);
#if defined(__SSE3__) && defined(__AVX__)
    c21_0 = _mm_add_sd(c21_0, _mm_mul_sd(a21_0, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c21_0 = _mm_add_sd(c21_0, _mm_mul_sd(a21_0, b21));
#endif
    _mm_store_sd(&C[(l_n*56)+21], c21_0);
    __m128d c21_1 = _mm_load_sd(&C[(l_n*56)+42]);
    __m128d a21_1 = _mm_load_sd(&A[93]);
#if defined(__SSE3__) && defined(__AVX__)
    c21_1 = _mm_add_sd(c21_1, _mm_mul_sd(a21_1, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c21_1 = _mm_add_sd(c21_1, _mm_mul_sd(a21_1, b21));
#endif
    _mm_store_sd(&C[(l_n*56)+42], c21_1);
#else
    C[(l_n*56)+21] += A[92] * B[(l_n*56)+21];
    C[(l_n*56)+42] += A[93] * B[(l_n*56)+21];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b22 = _mm256_broadcast_sd(&B[(l_n*56)+22]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b22 = _mm_loaddup_pd(&B[(l_n*56)+22]);
#endif
    __m128d c22_0 = _mm_load_sd(&C[(l_n*56)+22]);
    __m128d a22_0 = _mm_load_sd(&A[94]);
#if defined(__SSE3__) && defined(__AVX__)
    c22_0 = _mm_add_sd(c22_0, _mm_mul_sd(a22_0, _mm256_castpd256_pd128(b22)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c22_0 = _mm_add_sd(c22_0, _mm_mul_sd(a22_0, b22));
#endif
    _mm_store_sd(&C[(l_n*56)+22], c22_0);
    __m128d c22_1 = _mm_load_sd(&C[(l_n*56)+43]);
    __m128d a22_1 = _mm_load_sd(&A[95]);
#if defined(__SSE3__) && defined(__AVX__)
    c22_1 = _mm_add_sd(c22_1, _mm_mul_sd(a22_1, _mm256_castpd256_pd128(b22)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c22_1 = _mm_add_sd(c22_1, _mm_mul_sd(a22_1, b22));
#endif
    _mm_store_sd(&C[(l_n*56)+43], c22_1);
#else
    C[(l_n*56)+22] += A[94] * B[(l_n*56)+22];
    C[(l_n*56)+43] += A[95] * B[(l_n*56)+22];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b23 = _mm256_broadcast_sd(&B[(l_n*56)+23]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b23 = _mm_loaddup_pd(&B[(l_n*56)+23]);
#endif
    __m128d c23_0 = _mm_load_sd(&C[(l_n*56)+23]);
    __m128d a23_0 = _mm_load_sd(&A[96]);
#if defined(__SSE3__) && defined(__AVX__)
    c23_0 = _mm_add_sd(c23_0, _mm_mul_sd(a23_0, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c23_0 = _mm_add_sd(c23_0, _mm_mul_sd(a23_0, b23));
#endif
    _mm_store_sd(&C[(l_n*56)+23], c23_0);
    __m128d c23_1 = _mm_load_sd(&C[(l_n*56)+44]);
    __m128d a23_1 = _mm_load_sd(&A[97]);
#if defined(__SSE3__) && defined(__AVX__)
    c23_1 = _mm_add_sd(c23_1, _mm_mul_sd(a23_1, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c23_1 = _mm_add_sd(c23_1, _mm_mul_sd(a23_1, b23));
#endif
    _mm_store_sd(&C[(l_n*56)+44], c23_1);
#else
    C[(l_n*56)+23] += A[96] * B[(l_n*56)+23];
    C[(l_n*56)+44] += A[97] * B[(l_n*56)+23];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b24 = _mm256_broadcast_sd(&B[(l_n*56)+24]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b24 = _mm_loaddup_pd(&B[(l_n*56)+24]);
#endif
    __m128d c24_0 = _mm_load_sd(&C[(l_n*56)+24]);
    __m128d a24_0 = _mm_load_sd(&A[98]);
#if defined(__SSE3__) && defined(__AVX__)
    c24_0 = _mm_add_sd(c24_0, _mm_mul_sd(a24_0, _mm256_castpd256_pd128(b24)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c24_0 = _mm_add_sd(c24_0, _mm_mul_sd(a24_0, b24));
#endif
    _mm_store_sd(&C[(l_n*56)+24], c24_0);
    __m128d c24_1 = _mm_load_sd(&C[(l_n*56)+45]);
    __m128d a24_1 = _mm_load_sd(&A[99]);
#if defined(__SSE3__) && defined(__AVX__)
    c24_1 = _mm_add_sd(c24_1, _mm_mul_sd(a24_1, _mm256_castpd256_pd128(b24)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c24_1 = _mm_add_sd(c24_1, _mm_mul_sd(a24_1, b24));
#endif
    _mm_store_sd(&C[(l_n*56)+45], c24_1);
#else
    C[(l_n*56)+24] += A[98] * B[(l_n*56)+24];
    C[(l_n*56)+45] += A[99] * B[(l_n*56)+24];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b25 = _mm256_broadcast_sd(&B[(l_n*56)+25]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b25 = _mm_loaddup_pd(&B[(l_n*56)+25]);
#endif
    __m128d c25_0 = _mm_load_sd(&C[(l_n*56)+10]);
    __m128d a25_0 = _mm_load_sd(&A[100]);
#if defined(__SSE3__) && defined(__AVX__)
    c25_0 = _mm_add_sd(c25_0, _mm_mul_sd(a25_0, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c25_0 = _mm_add_sd(c25_0, _mm_mul_sd(a25_0, b25));
#endif
    _mm_store_sd(&C[(l_n*56)+10], c25_0);
    __m128d c25_1 = _mm_load_sd(&C[(l_n*56)+25]);
    __m128d a25_1 = _mm_load_sd(&A[101]);
#if defined(__SSE3__) && defined(__AVX__)
    c25_1 = _mm_add_sd(c25_1, _mm_mul_sd(a25_1, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c25_1 = _mm_add_sd(c25_1, _mm_mul_sd(a25_1, b25));
#endif
    _mm_store_sd(&C[(l_n*56)+25], c25_1);
    __m128d c25_2 = _mm_load_sd(&C[(l_n*56)+46]);
    __m128d a25_2 = _mm_load_sd(&A[102]);
#if defined(__SSE3__) && defined(__AVX__)
    c25_2 = _mm_add_sd(c25_2, _mm_mul_sd(a25_2, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c25_2 = _mm_add_sd(c25_2, _mm_mul_sd(a25_2, b25));
#endif
    _mm_store_sd(&C[(l_n*56)+46], c25_2);
#else
    C[(l_n*56)+10] += A[100] * B[(l_n*56)+25];
    C[(l_n*56)+25] += A[101] * B[(l_n*56)+25];
    C[(l_n*56)+46] += A[102] * B[(l_n*56)+25];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b26 = _mm256_broadcast_sd(&B[(l_n*56)+26]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b26 = _mm_loaddup_pd(&B[(l_n*56)+26]);
#endif
    __m128d c26_0 = _mm_load_sd(&C[(l_n*56)+11]);
    __m128d a26_0 = _mm_load_sd(&A[103]);
#if defined(__SSE3__) && defined(__AVX__)
    c26_0 = _mm_add_sd(c26_0, _mm_mul_sd(a26_0, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c26_0 = _mm_add_sd(c26_0, _mm_mul_sd(a26_0, b26));
#endif
    _mm_store_sd(&C[(l_n*56)+11], c26_0);
    __m128d c26_1 = _mm_load_sd(&C[(l_n*56)+26]);
    __m128d a26_1 = _mm_load_sd(&A[104]);
#if defined(__SSE3__) && defined(__AVX__)
    c26_1 = _mm_add_sd(c26_1, _mm_mul_sd(a26_1, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c26_1 = _mm_add_sd(c26_1, _mm_mul_sd(a26_1, b26));
#endif
    _mm_store_sd(&C[(l_n*56)+26], c26_1);
    __m128d c26_2 = _mm_load_sd(&C[(l_n*56)+47]);
    __m128d a26_2 = _mm_load_sd(&A[105]);
#if defined(__SSE3__) && defined(__AVX__)
    c26_2 = _mm_add_sd(c26_2, _mm_mul_sd(a26_2, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c26_2 = _mm_add_sd(c26_2, _mm_mul_sd(a26_2, b26));
#endif
    _mm_store_sd(&C[(l_n*56)+47], c26_2);
#else
    C[(l_n*56)+11] += A[103] * B[(l_n*56)+26];
    C[(l_n*56)+26] += A[104] * B[(l_n*56)+26];
    C[(l_n*56)+47] += A[105] * B[(l_n*56)+26];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b27 = _mm256_broadcast_sd(&B[(l_n*56)+27]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b27 = _mm_loaddup_pd(&B[(l_n*56)+27]);
#endif
    __m128d c27_0 = _mm_load_sd(&C[(l_n*56)+12]);
    __m128d a27_0 = _mm_load_sd(&A[106]);
#if defined(__SSE3__) && defined(__AVX__)
    c27_0 = _mm_add_sd(c27_0, _mm_mul_sd(a27_0, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c27_0 = _mm_add_sd(c27_0, _mm_mul_sd(a27_0, b27));
#endif
    _mm_store_sd(&C[(l_n*56)+12], c27_0);
    __m128d c27_1 = _mm_load_sd(&C[(l_n*56)+27]);
    __m128d a27_1 = _mm_load_sd(&A[107]);
#if defined(__SSE3__) && defined(__AVX__)
    c27_1 = _mm_add_sd(c27_1, _mm_mul_sd(a27_1, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c27_1 = _mm_add_sd(c27_1, _mm_mul_sd(a27_1, b27));
#endif
    _mm_store_sd(&C[(l_n*56)+27], c27_1);
    __m128d c27_2 = _mm_load_sd(&C[(l_n*56)+48]);
    __m128d a27_2 = _mm_load_sd(&A[108]);
#if defined(__SSE3__) && defined(__AVX__)
    c27_2 = _mm_add_sd(c27_2, _mm_mul_sd(a27_2, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c27_2 = _mm_add_sd(c27_2, _mm_mul_sd(a27_2, b27));
#endif
    _mm_store_sd(&C[(l_n*56)+48], c27_2);
#else
    C[(l_n*56)+12] += A[106] * B[(l_n*56)+27];
    C[(l_n*56)+27] += A[107] * B[(l_n*56)+27];
    C[(l_n*56)+48] += A[108] * B[(l_n*56)+27];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b28 = _mm256_broadcast_sd(&B[(l_n*56)+28]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b28 = _mm_loaddup_pd(&B[(l_n*56)+28]);
#endif
    __m128d c28_0 = _mm_load_sd(&C[(l_n*56)+13]);
    __m128d a28_0 = _mm_load_sd(&A[109]);
#if defined(__SSE3__) && defined(__AVX__)
    c28_0 = _mm_add_sd(c28_0, _mm_mul_sd(a28_0, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c28_0 = _mm_add_sd(c28_0, _mm_mul_sd(a28_0, b28));
#endif
    _mm_store_sd(&C[(l_n*56)+13], c28_0);
    __m128d c28_1 = _mm_load_sd(&C[(l_n*56)+28]);
    __m128d a28_1 = _mm_load_sd(&A[110]);
#if defined(__SSE3__) && defined(__AVX__)
    c28_1 = _mm_add_sd(c28_1, _mm_mul_sd(a28_1, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c28_1 = _mm_add_sd(c28_1, _mm_mul_sd(a28_1, b28));
#endif
    _mm_store_sd(&C[(l_n*56)+28], c28_1);
    __m128d c28_2 = _mm_load_sd(&C[(l_n*56)+49]);
    __m128d a28_2 = _mm_load_sd(&A[111]);
#if defined(__SSE3__) && defined(__AVX__)
    c28_2 = _mm_add_sd(c28_2, _mm_mul_sd(a28_2, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c28_2 = _mm_add_sd(c28_2, _mm_mul_sd(a28_2, b28));
#endif
    _mm_store_sd(&C[(l_n*56)+49], c28_2);
#else
    C[(l_n*56)+13] += A[109] * B[(l_n*56)+28];
    C[(l_n*56)+28] += A[110] * B[(l_n*56)+28];
    C[(l_n*56)+49] += A[111] * B[(l_n*56)+28];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b29 = _mm256_broadcast_sd(&B[(l_n*56)+29]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b29 = _mm_loaddup_pd(&B[(l_n*56)+29]);
#endif
    __m128d c29_0 = _mm_load_sd(&C[(l_n*56)+4]);
    __m128d a29_0 = _mm_load_sd(&A[112]);
#if defined(__SSE3__) && defined(__AVX__)
    c29_0 = _mm_add_sd(c29_0, _mm_mul_sd(a29_0, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c29_0 = _mm_add_sd(c29_0, _mm_mul_sd(a29_0, b29));
#endif
    _mm_store_sd(&C[(l_n*56)+4], c29_0);
    __m128d c29_1 = _mm_load_sd(&C[(l_n*56)+14]);
    __m128d a29_1 = _mm_load_sd(&A[113]);
#if defined(__SSE3__) && defined(__AVX__)
    c29_1 = _mm_add_sd(c29_1, _mm_mul_sd(a29_1, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c29_1 = _mm_add_sd(c29_1, _mm_mul_sd(a29_1, b29));
#endif
    _mm_store_sd(&C[(l_n*56)+14], c29_1);
    __m128d c29_2 = _mm_load_sd(&C[(l_n*56)+29]);
    __m128d a29_2 = _mm_load_sd(&A[114]);
#if defined(__SSE3__) && defined(__AVX__)
    c29_2 = _mm_add_sd(c29_2, _mm_mul_sd(a29_2, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c29_2 = _mm_add_sd(c29_2, _mm_mul_sd(a29_2, b29));
#endif
    _mm_store_sd(&C[(l_n*56)+29], c29_2);
    __m128d c29_3 = _mm_load_sd(&C[(l_n*56)+50]);
    __m128d a29_3 = _mm_load_sd(&A[115]);
#if defined(__SSE3__) && defined(__AVX__)
    c29_3 = _mm_add_sd(c29_3, _mm_mul_sd(a29_3, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c29_3 = _mm_add_sd(c29_3, _mm_mul_sd(a29_3, b29));
#endif
    _mm_store_sd(&C[(l_n*56)+50], c29_3);
#else
    C[(l_n*56)+4] += A[112] * B[(l_n*56)+29];
    C[(l_n*56)+14] += A[113] * B[(l_n*56)+29];
    C[(l_n*56)+29] += A[114] * B[(l_n*56)+29];
    C[(l_n*56)+50] += A[115] * B[(l_n*56)+29];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b30 = _mm256_broadcast_sd(&B[(l_n*56)+30]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b30 = _mm_loaddup_pd(&B[(l_n*56)+30]);
#endif
    __m128d c30_0 = _mm_load_sd(&C[(l_n*56)+5]);
    __m128d a30_0 = _mm_load_sd(&A[116]);
#if defined(__SSE3__) && defined(__AVX__)
    c30_0 = _mm_add_sd(c30_0, _mm_mul_sd(a30_0, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c30_0 = _mm_add_sd(c30_0, _mm_mul_sd(a30_0, b30));
#endif
    _mm_store_sd(&C[(l_n*56)+5], c30_0);
    __m128d c30_1 = _mm_load_sd(&C[(l_n*56)+15]);
    __m128d a30_1 = _mm_load_sd(&A[117]);
#if defined(__SSE3__) && defined(__AVX__)
    c30_1 = _mm_add_sd(c30_1, _mm_mul_sd(a30_1, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c30_1 = _mm_add_sd(c30_1, _mm_mul_sd(a30_1, b30));
#endif
    _mm_store_sd(&C[(l_n*56)+15], c30_1);
    __m128d c30_2 = _mm_load_sd(&C[(l_n*56)+30]);
    __m128d a30_2 = _mm_load_sd(&A[118]);
#if defined(__SSE3__) && defined(__AVX__)
    c30_2 = _mm_add_sd(c30_2, _mm_mul_sd(a30_2, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c30_2 = _mm_add_sd(c30_2, _mm_mul_sd(a30_2, b30));
#endif
    _mm_store_sd(&C[(l_n*56)+30], c30_2);
    __m128d c30_3 = _mm_load_sd(&C[(l_n*56)+51]);
    __m128d a30_3 = _mm_load_sd(&A[119]);
#if defined(__SSE3__) && defined(__AVX__)
    c30_3 = _mm_add_sd(c30_3, _mm_mul_sd(a30_3, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c30_3 = _mm_add_sd(c30_3, _mm_mul_sd(a30_3, b30));
#endif
    _mm_store_sd(&C[(l_n*56)+51], c30_3);
#else
    C[(l_n*56)+5] += A[116] * B[(l_n*56)+30];
    C[(l_n*56)+15] += A[117] * B[(l_n*56)+30];
    C[(l_n*56)+30] += A[118] * B[(l_n*56)+30];
    C[(l_n*56)+51] += A[119] * B[(l_n*56)+30];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b31 = _mm256_broadcast_sd(&B[(l_n*56)+31]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b31 = _mm_loaddup_pd(&B[(l_n*56)+31]);
#endif
    __m128d c31_0 = _mm_load_sd(&C[(l_n*56)+6]);
    __m128d a31_0 = _mm_load_sd(&A[120]);
#if defined(__SSE3__) && defined(__AVX__)
    c31_0 = _mm_add_sd(c31_0, _mm_mul_sd(a31_0, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c31_0 = _mm_add_sd(c31_0, _mm_mul_sd(a31_0, b31));
#endif
    _mm_store_sd(&C[(l_n*56)+6], c31_0);
    __m128d c31_1 = _mm_load_sd(&C[(l_n*56)+16]);
    __m128d a31_1 = _mm_load_sd(&A[121]);
#if defined(__SSE3__) && defined(__AVX__)
    c31_1 = _mm_add_sd(c31_1, _mm_mul_sd(a31_1, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c31_1 = _mm_add_sd(c31_1, _mm_mul_sd(a31_1, b31));
#endif
    _mm_store_sd(&C[(l_n*56)+16], c31_1);
    __m128d c31_2 = _mm_load_sd(&C[(l_n*56)+31]);
    __m128d a31_2 = _mm_load_sd(&A[122]);
#if defined(__SSE3__) && defined(__AVX__)
    c31_2 = _mm_add_sd(c31_2, _mm_mul_sd(a31_2, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c31_2 = _mm_add_sd(c31_2, _mm_mul_sd(a31_2, b31));
#endif
    _mm_store_sd(&C[(l_n*56)+31], c31_2);
    __m128d c31_3 = _mm_load_sd(&C[(l_n*56)+52]);
    __m128d a31_3 = _mm_load_sd(&A[123]);
#if defined(__SSE3__) && defined(__AVX__)
    c31_3 = _mm_add_sd(c31_3, _mm_mul_sd(a31_3, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c31_3 = _mm_add_sd(c31_3, _mm_mul_sd(a31_3, b31));
#endif
    _mm_store_sd(&C[(l_n*56)+52], c31_3);
#else
    C[(l_n*56)+6] += A[120] * B[(l_n*56)+31];
    C[(l_n*56)+16] += A[121] * B[(l_n*56)+31];
    C[(l_n*56)+31] += A[122] * B[(l_n*56)+31];
    C[(l_n*56)+52] += A[123] * B[(l_n*56)+31];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b32 = _mm256_broadcast_sd(&B[(l_n*56)+32]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b32 = _mm_loaddup_pd(&B[(l_n*56)+32]);
#endif
    __m128d c32_0 = _mm_load_sd(&C[(l_n*56)+1]);
    __m128d a32_0 = _mm_load_sd(&A[124]);
#if defined(__SSE3__) && defined(__AVX__)
    c32_0 = _mm_add_sd(c32_0, _mm_mul_sd(a32_0, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c32_0 = _mm_add_sd(c32_0, _mm_mul_sd(a32_0, b32));
#endif
    _mm_store_sd(&C[(l_n*56)+1], c32_0);
    __m128d c32_1 = _mm_load_sd(&C[(l_n*56)+7]);
    __m128d a32_1 = _mm_load_sd(&A[125]);
#if defined(__SSE3__) && defined(__AVX__)
    c32_1 = _mm_add_sd(c32_1, _mm_mul_sd(a32_1, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c32_1 = _mm_add_sd(c32_1, _mm_mul_sd(a32_1, b32));
#endif
    _mm_store_sd(&C[(l_n*56)+7], c32_1);
    __m128d c32_2 = _mm_load_sd(&C[(l_n*56)+17]);
    __m128d a32_2 = _mm_load_sd(&A[126]);
#if defined(__SSE3__) && defined(__AVX__)
    c32_2 = _mm_add_sd(c32_2, _mm_mul_sd(a32_2, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c32_2 = _mm_add_sd(c32_2, _mm_mul_sd(a32_2, b32));
#endif
    _mm_store_sd(&C[(l_n*56)+17], c32_2);
    __m128d c32_3 = _mm_load_sd(&C[(l_n*56)+32]);
    __m128d a32_3 = _mm_load_sd(&A[127]);
#if defined(__SSE3__) && defined(__AVX__)
    c32_3 = _mm_add_sd(c32_3, _mm_mul_sd(a32_3, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c32_3 = _mm_add_sd(c32_3, _mm_mul_sd(a32_3, b32));
#endif
    _mm_store_sd(&C[(l_n*56)+32], c32_3);
    __m128d c32_4 = _mm_load_sd(&C[(l_n*56)+53]);
    __m128d a32_4 = _mm_load_sd(&A[128]);
#if defined(__SSE3__) && defined(__AVX__)
    c32_4 = _mm_add_sd(c32_4, _mm_mul_sd(a32_4, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c32_4 = _mm_add_sd(c32_4, _mm_mul_sd(a32_4, b32));
#endif
    _mm_store_sd(&C[(l_n*56)+53], c32_4);
#else
    C[(l_n*56)+1] += A[124] * B[(l_n*56)+32];
    C[(l_n*56)+7] += A[125] * B[(l_n*56)+32];
    C[(l_n*56)+17] += A[126] * B[(l_n*56)+32];
    C[(l_n*56)+32] += A[127] * B[(l_n*56)+32];
    C[(l_n*56)+53] += A[128] * B[(l_n*56)+32];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b33 = _mm256_broadcast_sd(&B[(l_n*56)+33]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b33 = _mm_loaddup_pd(&B[(l_n*56)+33]);
#endif
    __m128d c33_0 = _mm_load_sd(&C[(l_n*56)+2]);
    __m128d a33_0 = _mm_load_sd(&A[129]);
#if defined(__SSE3__) && defined(__AVX__)
    c33_0 = _mm_add_sd(c33_0, _mm_mul_sd(a33_0, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c33_0 = _mm_add_sd(c33_0, _mm_mul_sd(a33_0, b33));
#endif
    _mm_store_sd(&C[(l_n*56)+2], c33_0);
    __m128d c33_1 = _mm_load_sd(&C[(l_n*56)+8]);
    __m128d a33_1 = _mm_load_sd(&A[130]);
#if defined(__SSE3__) && defined(__AVX__)
    c33_1 = _mm_add_sd(c33_1, _mm_mul_sd(a33_1, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c33_1 = _mm_add_sd(c33_1, _mm_mul_sd(a33_1, b33));
#endif
    _mm_store_sd(&C[(l_n*56)+8], c33_1);
    __m128d c33_2 = _mm_load_sd(&C[(l_n*56)+18]);
    __m128d a33_2 = _mm_load_sd(&A[131]);
#if defined(__SSE3__) && defined(__AVX__)
    c33_2 = _mm_add_sd(c33_2, _mm_mul_sd(a33_2, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c33_2 = _mm_add_sd(c33_2, _mm_mul_sd(a33_2, b33));
#endif
    _mm_store_sd(&C[(l_n*56)+18], c33_2);
    __m128d c33_3 = _mm_load_sd(&C[(l_n*56)+33]);
    __m128d a33_3 = _mm_load_sd(&A[132]);
#if defined(__SSE3__) && defined(__AVX__)
    c33_3 = _mm_add_sd(c33_3, _mm_mul_sd(a33_3, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c33_3 = _mm_add_sd(c33_3, _mm_mul_sd(a33_3, b33));
#endif
    _mm_store_sd(&C[(l_n*56)+33], c33_3);
    __m128d c33_4 = _mm_load_sd(&C[(l_n*56)+54]);
    __m128d a33_4 = _mm_load_sd(&A[133]);
#if defined(__SSE3__) && defined(__AVX__)
    c33_4 = _mm_add_sd(c33_4, _mm_mul_sd(a33_4, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c33_4 = _mm_add_sd(c33_4, _mm_mul_sd(a33_4, b33));
#endif
    _mm_store_sd(&C[(l_n*56)+54], c33_4);
#else
    C[(l_n*56)+2] += A[129] * B[(l_n*56)+33];
    C[(l_n*56)+8] += A[130] * B[(l_n*56)+33];
    C[(l_n*56)+18] += A[131] * B[(l_n*56)+33];
    C[(l_n*56)+33] += A[132] * B[(l_n*56)+33];
    C[(l_n*56)+54] += A[133] * B[(l_n*56)+33];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b34 = _mm256_broadcast_sd(&B[(l_n*56)+34]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b34 = _mm_loaddup_pd(&B[(l_n*56)+34]);
#endif
    __m128d c34_0 = _mm_load_sd(&C[(l_n*56)+0]);
    __m128d a34_0 = _mm_load_sd(&A[134]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_0 = _mm_add_sd(c34_0, _mm_mul_sd(a34_0, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_0 = _mm_add_sd(c34_0, _mm_mul_sd(a34_0, b34));
#endif
    _mm_store_sd(&C[(l_n*56)+0], c34_0);
    __m128d c34_1 = _mm_load_sd(&C[(l_n*56)+3]);
    __m128d a34_1 = _mm_load_sd(&A[135]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_1 = _mm_add_sd(c34_1, _mm_mul_sd(a34_1, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_1 = _mm_add_sd(c34_1, _mm_mul_sd(a34_1, b34));
#endif
    _mm_store_sd(&C[(l_n*56)+3], c34_1);
    __m128d c34_2 = _mm_load_sd(&C[(l_n*56)+9]);
    __m128d a34_2 = _mm_load_sd(&A[136]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_2 = _mm_add_sd(c34_2, _mm_mul_sd(a34_2, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_2 = _mm_add_sd(c34_2, _mm_mul_sd(a34_2, b34));
#endif
    _mm_store_sd(&C[(l_n*56)+9], c34_2);
    __m128d c34_3 = _mm_load_sd(&C[(l_n*56)+19]);
    __m128d a34_3 = _mm_load_sd(&A[137]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_3 = _mm_add_sd(c34_3, _mm_mul_sd(a34_3, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_3 = _mm_add_sd(c34_3, _mm_mul_sd(a34_3, b34));
#endif
    _mm_store_sd(&C[(l_n*56)+19], c34_3);
    __m128d c34_4 = _mm_load_sd(&C[(l_n*56)+34]);
    __m128d a34_4 = _mm_load_sd(&A[138]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_4 = _mm_add_sd(c34_4, _mm_mul_sd(a34_4, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_4 = _mm_add_sd(c34_4, _mm_mul_sd(a34_4, b34));
#endif
    _mm_store_sd(&C[(l_n*56)+34], c34_4);
    __m128d c34_5 = _mm_load_sd(&C[(l_n*56)+55]);
    __m128d a34_5 = _mm_load_sd(&A[139]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_5 = _mm_add_sd(c34_5, _mm_mul_sd(a34_5, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_5 = _mm_add_sd(c34_5, _mm_mul_sd(a34_5, b34));
#endif
    _mm_store_sd(&C[(l_n*56)+55], c34_5);
#else
    C[(l_n*56)+0] += A[134] * B[(l_n*56)+34];
    C[(l_n*56)+3] += A[135] * B[(l_n*56)+34];
    C[(l_n*56)+9] += A[136] * B[(l_n*56)+34];
    C[(l_n*56)+19] += A[137] * B[(l_n*56)+34];
    C[(l_n*56)+34] += A[138] * B[(l_n*56)+34];
    C[(l_n*56)+55] += A[139] * B[(l_n*56)+34];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b35 = _mm256_broadcast_sd(&B[(l_n*56)+35]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b35 = _mm_loaddup_pd(&B[(l_n*56)+35]);
#endif
    __m128d c35_0 = _mm_load_sd(&C[(l_n*56)+35]);
    __m128d a35_0 = _mm_load_sd(&A[140]);
#if defined(__SSE3__) && defined(__AVX__)
    c35_0 = _mm_add_sd(c35_0, _mm_mul_sd(a35_0, _mm256_castpd256_pd128(b35)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c35_0 = _mm_add_sd(c35_0, _mm_mul_sd(a35_0, b35));
#endif
    _mm_store_sd(&C[(l_n*56)+35], c35_0);
#else
    C[(l_n*56)+35] += A[140] * B[(l_n*56)+35];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b36 = _mm256_broadcast_sd(&B[(l_n*56)+36]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b36 = _mm_loaddup_pd(&B[(l_n*56)+36]);
#endif
    __m128d c36_0 = _mm_load_sd(&C[(l_n*56)+36]);
    __m128d a36_0 = _mm_load_sd(&A[141]);
#if defined(__SSE3__) && defined(__AVX__)
    c36_0 = _mm_add_sd(c36_0, _mm_mul_sd(a36_0, _mm256_castpd256_pd128(b36)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c36_0 = _mm_add_sd(c36_0, _mm_mul_sd(a36_0, b36));
#endif
    _mm_store_sd(&C[(l_n*56)+36], c36_0);
#else
    C[(l_n*56)+36] += A[141] * B[(l_n*56)+36];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b37 = _mm256_broadcast_sd(&B[(l_n*56)+37]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b37 = _mm_loaddup_pd(&B[(l_n*56)+37]);
#endif
    __m128d c37_0 = _mm_load_sd(&C[(l_n*56)+37]);
    __m128d a37_0 = _mm_load_sd(&A[142]);
#if defined(__SSE3__) && defined(__AVX__)
    c37_0 = _mm_add_sd(c37_0, _mm_mul_sd(a37_0, _mm256_castpd256_pd128(b37)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c37_0 = _mm_add_sd(c37_0, _mm_mul_sd(a37_0, b37));
#endif
    _mm_store_sd(&C[(l_n*56)+37], c37_0);
#else
    C[(l_n*56)+37] += A[142] * B[(l_n*56)+37];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b38 = _mm256_broadcast_sd(&B[(l_n*56)+38]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b38 = _mm_loaddup_pd(&B[(l_n*56)+38]);
#endif
    __m128d c38_0 = _mm_load_sd(&C[(l_n*56)+38]);
    __m128d a38_0 = _mm_load_sd(&A[143]);
#if defined(__SSE3__) && defined(__AVX__)
    c38_0 = _mm_add_sd(c38_0, _mm_mul_sd(a38_0, _mm256_castpd256_pd128(b38)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c38_0 = _mm_add_sd(c38_0, _mm_mul_sd(a38_0, b38));
#endif
    _mm_store_sd(&C[(l_n*56)+38], c38_0);
#else
    C[(l_n*56)+38] += A[143] * B[(l_n*56)+38];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b39 = _mm256_broadcast_sd(&B[(l_n*56)+39]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b39 = _mm_loaddup_pd(&B[(l_n*56)+39]);
#endif
    __m128d c39_0 = _mm_load_sd(&C[(l_n*56)+39]);
    __m128d a39_0 = _mm_load_sd(&A[144]);
#if defined(__SSE3__) && defined(__AVX__)
    c39_0 = _mm_add_sd(c39_0, _mm_mul_sd(a39_0, _mm256_castpd256_pd128(b39)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c39_0 = _mm_add_sd(c39_0, _mm_mul_sd(a39_0, b39));
#endif
    _mm_store_sd(&C[(l_n*56)+39], c39_0);
#else
    C[(l_n*56)+39] += A[144] * B[(l_n*56)+39];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b40 = _mm256_broadcast_sd(&B[(l_n*56)+40]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b40 = _mm_loaddup_pd(&B[(l_n*56)+40]);
#endif
    __m128d c40_0 = _mm_load_sd(&C[(l_n*56)+40]);
    __m128d a40_0 = _mm_load_sd(&A[145]);
#if defined(__SSE3__) && defined(__AVX__)
    c40_0 = _mm_add_sd(c40_0, _mm_mul_sd(a40_0, _mm256_castpd256_pd128(b40)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c40_0 = _mm_add_sd(c40_0, _mm_mul_sd(a40_0, b40));
#endif
    _mm_store_sd(&C[(l_n*56)+40], c40_0);
#else
    C[(l_n*56)+40] += A[145] * B[(l_n*56)+40];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b41 = _mm256_broadcast_sd(&B[(l_n*56)+41]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b41 = _mm_loaddup_pd(&B[(l_n*56)+41]);
#endif
    __m128d c41_0 = _mm_load_sd(&C[(l_n*56)+20]);
    __m128d a41_0 = _mm_load_sd(&A[146]);
#if defined(__SSE3__) && defined(__AVX__)
    c41_0 = _mm_add_sd(c41_0, _mm_mul_sd(a41_0, _mm256_castpd256_pd128(b41)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c41_0 = _mm_add_sd(c41_0, _mm_mul_sd(a41_0, b41));
#endif
    _mm_store_sd(&C[(l_n*56)+20], c41_0);
    __m128d c41_1 = _mm_load_sd(&C[(l_n*56)+41]);
    __m128d a41_1 = _mm_load_sd(&A[147]);
#if defined(__SSE3__) && defined(__AVX__)
    c41_1 = _mm_add_sd(c41_1, _mm_mul_sd(a41_1, _mm256_castpd256_pd128(b41)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c41_1 = _mm_add_sd(c41_1, _mm_mul_sd(a41_1, b41));
#endif
    _mm_store_sd(&C[(l_n*56)+41], c41_1);
#else
    C[(l_n*56)+20] += A[146] * B[(l_n*56)+41];
    C[(l_n*56)+41] += A[147] * B[(l_n*56)+41];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b42 = _mm256_broadcast_sd(&B[(l_n*56)+42]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b42 = _mm_loaddup_pd(&B[(l_n*56)+42]);
#endif
    __m128d c42_0 = _mm_load_sd(&C[(l_n*56)+21]);
    __m128d a42_0 = _mm_load_sd(&A[148]);
#if defined(__SSE3__) && defined(__AVX__)
    c42_0 = _mm_add_sd(c42_0, _mm_mul_sd(a42_0, _mm256_castpd256_pd128(b42)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c42_0 = _mm_add_sd(c42_0, _mm_mul_sd(a42_0, b42));
#endif
    _mm_store_sd(&C[(l_n*56)+21], c42_0);
    __m128d c42_1 = _mm_load_sd(&C[(l_n*56)+42]);
    __m128d a42_1 = _mm_load_sd(&A[149]);
#if defined(__SSE3__) && defined(__AVX__)
    c42_1 = _mm_add_sd(c42_1, _mm_mul_sd(a42_1, _mm256_castpd256_pd128(b42)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c42_1 = _mm_add_sd(c42_1, _mm_mul_sd(a42_1, b42));
#endif
    _mm_store_sd(&C[(l_n*56)+42], c42_1);
#else
    C[(l_n*56)+21] += A[148] * B[(l_n*56)+42];
    C[(l_n*56)+42] += A[149] * B[(l_n*56)+42];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b43 = _mm256_broadcast_sd(&B[(l_n*56)+43]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b43 = _mm_loaddup_pd(&B[(l_n*56)+43]);
#endif
    __m128d c43_0 = _mm_load_sd(&C[(l_n*56)+22]);
    __m128d a43_0 = _mm_load_sd(&A[150]);
#if defined(__SSE3__) && defined(__AVX__)
    c43_0 = _mm_add_sd(c43_0, _mm_mul_sd(a43_0, _mm256_castpd256_pd128(b43)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c43_0 = _mm_add_sd(c43_0, _mm_mul_sd(a43_0, b43));
#endif
    _mm_store_sd(&C[(l_n*56)+22], c43_0);
    __m128d c43_1 = _mm_load_sd(&C[(l_n*56)+43]);
    __m128d a43_1 = _mm_load_sd(&A[151]);
#if defined(__SSE3__) && defined(__AVX__)
    c43_1 = _mm_add_sd(c43_1, _mm_mul_sd(a43_1, _mm256_castpd256_pd128(b43)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c43_1 = _mm_add_sd(c43_1, _mm_mul_sd(a43_1, b43));
#endif
    _mm_store_sd(&C[(l_n*56)+43], c43_1);
#else
    C[(l_n*56)+22] += A[150] * B[(l_n*56)+43];
    C[(l_n*56)+43] += A[151] * B[(l_n*56)+43];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b44 = _mm256_broadcast_sd(&B[(l_n*56)+44]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b44 = _mm_loaddup_pd(&B[(l_n*56)+44]);
#endif
    __m128d c44_0 = _mm_load_sd(&C[(l_n*56)+23]);
    __m128d a44_0 = _mm_load_sd(&A[152]);
#if defined(__SSE3__) && defined(__AVX__)
    c44_0 = _mm_add_sd(c44_0, _mm_mul_sd(a44_0, _mm256_castpd256_pd128(b44)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c44_0 = _mm_add_sd(c44_0, _mm_mul_sd(a44_0, b44));
#endif
    _mm_store_sd(&C[(l_n*56)+23], c44_0);
    __m128d c44_1 = _mm_load_sd(&C[(l_n*56)+44]);
    __m128d a44_1 = _mm_load_sd(&A[153]);
#if defined(__SSE3__) && defined(__AVX__)
    c44_1 = _mm_add_sd(c44_1, _mm_mul_sd(a44_1, _mm256_castpd256_pd128(b44)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c44_1 = _mm_add_sd(c44_1, _mm_mul_sd(a44_1, b44));
#endif
    _mm_store_sd(&C[(l_n*56)+44], c44_1);
#else
    C[(l_n*56)+23] += A[152] * B[(l_n*56)+44];
    C[(l_n*56)+44] += A[153] * B[(l_n*56)+44];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b45 = _mm256_broadcast_sd(&B[(l_n*56)+45]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b45 = _mm_loaddup_pd(&B[(l_n*56)+45]);
#endif
    __m128d c45_0 = _mm_load_sd(&C[(l_n*56)+24]);
    __m128d a45_0 = _mm_load_sd(&A[154]);
#if defined(__SSE3__) && defined(__AVX__)
    c45_0 = _mm_add_sd(c45_0, _mm_mul_sd(a45_0, _mm256_castpd256_pd128(b45)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c45_0 = _mm_add_sd(c45_0, _mm_mul_sd(a45_0, b45));
#endif
    _mm_store_sd(&C[(l_n*56)+24], c45_0);
    __m128d c45_1 = _mm_load_sd(&C[(l_n*56)+45]);
    __m128d a45_1 = _mm_load_sd(&A[155]);
#if defined(__SSE3__) && defined(__AVX__)
    c45_1 = _mm_add_sd(c45_1, _mm_mul_sd(a45_1, _mm256_castpd256_pd128(b45)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c45_1 = _mm_add_sd(c45_1, _mm_mul_sd(a45_1, b45));
#endif
    _mm_store_sd(&C[(l_n*56)+45], c45_1);
#else
    C[(l_n*56)+24] += A[154] * B[(l_n*56)+45];
    C[(l_n*56)+45] += A[155] * B[(l_n*56)+45];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b46 = _mm256_broadcast_sd(&B[(l_n*56)+46]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b46 = _mm_loaddup_pd(&B[(l_n*56)+46]);
#endif
    __m128d c46_0 = _mm_load_sd(&C[(l_n*56)+10]);
    __m128d a46_0 = _mm_load_sd(&A[156]);
#if defined(__SSE3__) && defined(__AVX__)
    c46_0 = _mm_add_sd(c46_0, _mm_mul_sd(a46_0, _mm256_castpd256_pd128(b46)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c46_0 = _mm_add_sd(c46_0, _mm_mul_sd(a46_0, b46));
#endif
    _mm_store_sd(&C[(l_n*56)+10], c46_0);
    __m128d c46_1 = _mm_load_sd(&C[(l_n*56)+25]);
    __m128d a46_1 = _mm_load_sd(&A[157]);
#if defined(__SSE3__) && defined(__AVX__)
    c46_1 = _mm_add_sd(c46_1, _mm_mul_sd(a46_1, _mm256_castpd256_pd128(b46)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c46_1 = _mm_add_sd(c46_1, _mm_mul_sd(a46_1, b46));
#endif
    _mm_store_sd(&C[(l_n*56)+25], c46_1);
    __m128d c46_2 = _mm_load_sd(&C[(l_n*56)+46]);
    __m128d a46_2 = _mm_load_sd(&A[158]);
#if defined(__SSE3__) && defined(__AVX__)
    c46_2 = _mm_add_sd(c46_2, _mm_mul_sd(a46_2, _mm256_castpd256_pd128(b46)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c46_2 = _mm_add_sd(c46_2, _mm_mul_sd(a46_2, b46));
#endif
    _mm_store_sd(&C[(l_n*56)+46], c46_2);
#else
    C[(l_n*56)+10] += A[156] * B[(l_n*56)+46];
    C[(l_n*56)+25] += A[157] * B[(l_n*56)+46];
    C[(l_n*56)+46] += A[158] * B[(l_n*56)+46];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b47 = _mm256_broadcast_sd(&B[(l_n*56)+47]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b47 = _mm_loaddup_pd(&B[(l_n*56)+47]);
#endif
    __m128d c47_0 = _mm_load_sd(&C[(l_n*56)+11]);
    __m128d a47_0 = _mm_load_sd(&A[159]);
#if defined(__SSE3__) && defined(__AVX__)
    c47_0 = _mm_add_sd(c47_0, _mm_mul_sd(a47_0, _mm256_castpd256_pd128(b47)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c47_0 = _mm_add_sd(c47_0, _mm_mul_sd(a47_0, b47));
#endif
    _mm_store_sd(&C[(l_n*56)+11], c47_0);
    __m128d c47_1 = _mm_load_sd(&C[(l_n*56)+26]);
    __m128d a47_1 = _mm_load_sd(&A[160]);
#if defined(__SSE3__) && defined(__AVX__)
    c47_1 = _mm_add_sd(c47_1, _mm_mul_sd(a47_1, _mm256_castpd256_pd128(b47)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c47_1 = _mm_add_sd(c47_1, _mm_mul_sd(a47_1, b47));
#endif
    _mm_store_sd(&C[(l_n*56)+26], c47_1);
    __m128d c47_2 = _mm_load_sd(&C[(l_n*56)+47]);
    __m128d a47_2 = _mm_load_sd(&A[161]);
#if defined(__SSE3__) && defined(__AVX__)
    c47_2 = _mm_add_sd(c47_2, _mm_mul_sd(a47_2, _mm256_castpd256_pd128(b47)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c47_2 = _mm_add_sd(c47_2, _mm_mul_sd(a47_2, b47));
#endif
    _mm_store_sd(&C[(l_n*56)+47], c47_2);
#else
    C[(l_n*56)+11] += A[159] * B[(l_n*56)+47];
    C[(l_n*56)+26] += A[160] * B[(l_n*56)+47];
    C[(l_n*56)+47] += A[161] * B[(l_n*56)+47];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b48 = _mm256_broadcast_sd(&B[(l_n*56)+48]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b48 = _mm_loaddup_pd(&B[(l_n*56)+48]);
#endif
    __m128d c48_0 = _mm_load_sd(&C[(l_n*56)+12]);
    __m128d a48_0 = _mm_load_sd(&A[162]);
#if defined(__SSE3__) && defined(__AVX__)
    c48_0 = _mm_add_sd(c48_0, _mm_mul_sd(a48_0, _mm256_castpd256_pd128(b48)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c48_0 = _mm_add_sd(c48_0, _mm_mul_sd(a48_0, b48));
#endif
    _mm_store_sd(&C[(l_n*56)+12], c48_0);
    __m128d c48_1 = _mm_load_sd(&C[(l_n*56)+27]);
    __m128d a48_1 = _mm_load_sd(&A[163]);
#if defined(__SSE3__) && defined(__AVX__)
    c48_1 = _mm_add_sd(c48_1, _mm_mul_sd(a48_1, _mm256_castpd256_pd128(b48)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c48_1 = _mm_add_sd(c48_1, _mm_mul_sd(a48_1, b48));
#endif
    _mm_store_sd(&C[(l_n*56)+27], c48_1);
    __m128d c48_2 = _mm_load_sd(&C[(l_n*56)+48]);
    __m128d a48_2 = _mm_load_sd(&A[164]);
#if defined(__SSE3__) && defined(__AVX__)
    c48_2 = _mm_add_sd(c48_2, _mm_mul_sd(a48_2, _mm256_castpd256_pd128(b48)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c48_2 = _mm_add_sd(c48_2, _mm_mul_sd(a48_2, b48));
#endif
    _mm_store_sd(&C[(l_n*56)+48], c48_2);
#else
    C[(l_n*56)+12] += A[162] * B[(l_n*56)+48];
    C[(l_n*56)+27] += A[163] * B[(l_n*56)+48];
    C[(l_n*56)+48] += A[164] * B[(l_n*56)+48];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b49 = _mm256_broadcast_sd(&B[(l_n*56)+49]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b49 = _mm_loaddup_pd(&B[(l_n*56)+49]);
#endif
    __m128d c49_0 = _mm_load_sd(&C[(l_n*56)+13]);
    __m128d a49_0 = _mm_load_sd(&A[165]);
#if defined(__SSE3__) && defined(__AVX__)
    c49_0 = _mm_add_sd(c49_0, _mm_mul_sd(a49_0, _mm256_castpd256_pd128(b49)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c49_0 = _mm_add_sd(c49_0, _mm_mul_sd(a49_0, b49));
#endif
    _mm_store_sd(&C[(l_n*56)+13], c49_0);
    __m128d c49_1 = _mm_load_sd(&C[(l_n*56)+28]);
    __m128d a49_1 = _mm_load_sd(&A[166]);
#if defined(__SSE3__) && defined(__AVX__)
    c49_1 = _mm_add_sd(c49_1, _mm_mul_sd(a49_1, _mm256_castpd256_pd128(b49)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c49_1 = _mm_add_sd(c49_1, _mm_mul_sd(a49_1, b49));
#endif
    _mm_store_sd(&C[(l_n*56)+28], c49_1);
    __m128d c49_2 = _mm_load_sd(&C[(l_n*56)+49]);
    __m128d a49_2 = _mm_load_sd(&A[167]);
#if defined(__SSE3__) && defined(__AVX__)
    c49_2 = _mm_add_sd(c49_2, _mm_mul_sd(a49_2, _mm256_castpd256_pd128(b49)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c49_2 = _mm_add_sd(c49_2, _mm_mul_sd(a49_2, b49));
#endif
    _mm_store_sd(&C[(l_n*56)+49], c49_2);
#else
    C[(l_n*56)+13] += A[165] * B[(l_n*56)+49];
    C[(l_n*56)+28] += A[166] * B[(l_n*56)+49];
    C[(l_n*56)+49] += A[167] * B[(l_n*56)+49];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b50 = _mm256_broadcast_sd(&B[(l_n*56)+50]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b50 = _mm_loaddup_pd(&B[(l_n*56)+50]);
#endif
    __m128d c50_0 = _mm_load_sd(&C[(l_n*56)+4]);
    __m128d a50_0 = _mm_load_sd(&A[168]);
#if defined(__SSE3__) && defined(__AVX__)
    c50_0 = _mm_add_sd(c50_0, _mm_mul_sd(a50_0, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c50_0 = _mm_add_sd(c50_0, _mm_mul_sd(a50_0, b50));
#endif
    _mm_store_sd(&C[(l_n*56)+4], c50_0);
    __m128d c50_1 = _mm_load_sd(&C[(l_n*56)+14]);
    __m128d a50_1 = _mm_load_sd(&A[169]);
#if defined(__SSE3__) && defined(__AVX__)
    c50_1 = _mm_add_sd(c50_1, _mm_mul_sd(a50_1, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c50_1 = _mm_add_sd(c50_1, _mm_mul_sd(a50_1, b50));
#endif
    _mm_store_sd(&C[(l_n*56)+14], c50_1);
    __m128d c50_2 = _mm_load_sd(&C[(l_n*56)+29]);
    __m128d a50_2 = _mm_load_sd(&A[170]);
#if defined(__SSE3__) && defined(__AVX__)
    c50_2 = _mm_add_sd(c50_2, _mm_mul_sd(a50_2, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c50_2 = _mm_add_sd(c50_2, _mm_mul_sd(a50_2, b50));
#endif
    _mm_store_sd(&C[(l_n*56)+29], c50_2);
    __m128d c50_3 = _mm_load_sd(&C[(l_n*56)+50]);
    __m128d a50_3 = _mm_load_sd(&A[171]);
#if defined(__SSE3__) && defined(__AVX__)
    c50_3 = _mm_add_sd(c50_3, _mm_mul_sd(a50_3, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c50_3 = _mm_add_sd(c50_3, _mm_mul_sd(a50_3, b50));
#endif
    _mm_store_sd(&C[(l_n*56)+50], c50_3);
#else
    C[(l_n*56)+4] += A[168] * B[(l_n*56)+50];
    C[(l_n*56)+14] += A[169] * B[(l_n*56)+50];
    C[(l_n*56)+29] += A[170] * B[(l_n*56)+50];
    C[(l_n*56)+50] += A[171] * B[(l_n*56)+50];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b51 = _mm256_broadcast_sd(&B[(l_n*56)+51]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b51 = _mm_loaddup_pd(&B[(l_n*56)+51]);
#endif
    __m128d c51_0 = _mm_load_sd(&C[(l_n*56)+5]);
    __m128d a51_0 = _mm_load_sd(&A[172]);
#if defined(__SSE3__) && defined(__AVX__)
    c51_0 = _mm_add_sd(c51_0, _mm_mul_sd(a51_0, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c51_0 = _mm_add_sd(c51_0, _mm_mul_sd(a51_0, b51));
#endif
    _mm_store_sd(&C[(l_n*56)+5], c51_0);
    __m128d c51_1 = _mm_load_sd(&C[(l_n*56)+15]);
    __m128d a51_1 = _mm_load_sd(&A[173]);
#if defined(__SSE3__) && defined(__AVX__)
    c51_1 = _mm_add_sd(c51_1, _mm_mul_sd(a51_1, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c51_1 = _mm_add_sd(c51_1, _mm_mul_sd(a51_1, b51));
#endif
    _mm_store_sd(&C[(l_n*56)+15], c51_1);
    __m128d c51_2 = _mm_load_sd(&C[(l_n*56)+30]);
    __m128d a51_2 = _mm_load_sd(&A[174]);
#if defined(__SSE3__) && defined(__AVX__)
    c51_2 = _mm_add_sd(c51_2, _mm_mul_sd(a51_2, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c51_2 = _mm_add_sd(c51_2, _mm_mul_sd(a51_2, b51));
#endif
    _mm_store_sd(&C[(l_n*56)+30], c51_2);
    __m128d c51_3 = _mm_load_sd(&C[(l_n*56)+51]);
    __m128d a51_3 = _mm_load_sd(&A[175]);
#if defined(__SSE3__) && defined(__AVX__)
    c51_3 = _mm_add_sd(c51_3, _mm_mul_sd(a51_3, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c51_3 = _mm_add_sd(c51_3, _mm_mul_sd(a51_3, b51));
#endif
    _mm_store_sd(&C[(l_n*56)+51], c51_3);
#else
    C[(l_n*56)+5] += A[172] * B[(l_n*56)+51];
    C[(l_n*56)+15] += A[173] * B[(l_n*56)+51];
    C[(l_n*56)+30] += A[174] * B[(l_n*56)+51];
    C[(l_n*56)+51] += A[175] * B[(l_n*56)+51];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b52 = _mm256_broadcast_sd(&B[(l_n*56)+52]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b52 = _mm_loaddup_pd(&B[(l_n*56)+52]);
#endif
    __m128d c52_0 = _mm_load_sd(&C[(l_n*56)+6]);
    __m128d a52_0 = _mm_load_sd(&A[176]);
#if defined(__SSE3__) && defined(__AVX__)
    c52_0 = _mm_add_sd(c52_0, _mm_mul_sd(a52_0, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c52_0 = _mm_add_sd(c52_0, _mm_mul_sd(a52_0, b52));
#endif
    _mm_store_sd(&C[(l_n*56)+6], c52_0);
    __m128d c52_1 = _mm_load_sd(&C[(l_n*56)+16]);
    __m128d a52_1 = _mm_load_sd(&A[177]);
#if defined(__SSE3__) && defined(__AVX__)
    c52_1 = _mm_add_sd(c52_1, _mm_mul_sd(a52_1, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c52_1 = _mm_add_sd(c52_1, _mm_mul_sd(a52_1, b52));
#endif
    _mm_store_sd(&C[(l_n*56)+16], c52_1);
    __m128d c52_2 = _mm_load_sd(&C[(l_n*56)+31]);
    __m128d a52_2 = _mm_load_sd(&A[178]);
#if defined(__SSE3__) && defined(__AVX__)
    c52_2 = _mm_add_sd(c52_2, _mm_mul_sd(a52_2, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c52_2 = _mm_add_sd(c52_2, _mm_mul_sd(a52_2, b52));
#endif
    _mm_store_sd(&C[(l_n*56)+31], c52_2);
    __m128d c52_3 = _mm_load_sd(&C[(l_n*56)+52]);
    __m128d a52_3 = _mm_load_sd(&A[179]);
#if defined(__SSE3__) && defined(__AVX__)
    c52_3 = _mm_add_sd(c52_3, _mm_mul_sd(a52_3, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c52_3 = _mm_add_sd(c52_3, _mm_mul_sd(a52_3, b52));
#endif
    _mm_store_sd(&C[(l_n*56)+52], c52_3);
#else
    C[(l_n*56)+6] += A[176] * B[(l_n*56)+52];
    C[(l_n*56)+16] += A[177] * B[(l_n*56)+52];
    C[(l_n*56)+31] += A[178] * B[(l_n*56)+52];
    C[(l_n*56)+52] += A[179] * B[(l_n*56)+52];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b53 = _mm256_broadcast_sd(&B[(l_n*56)+53]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b53 = _mm_loaddup_pd(&B[(l_n*56)+53]);
#endif
    __m128d c53_0 = _mm_load_sd(&C[(l_n*56)+1]);
    __m128d a53_0 = _mm_load_sd(&A[180]);
#if defined(__SSE3__) && defined(__AVX__)
    c53_0 = _mm_add_sd(c53_0, _mm_mul_sd(a53_0, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c53_0 = _mm_add_sd(c53_0, _mm_mul_sd(a53_0, b53));
#endif
    _mm_store_sd(&C[(l_n*56)+1], c53_0);
    __m128d c53_1 = _mm_load_sd(&C[(l_n*56)+7]);
    __m128d a53_1 = _mm_load_sd(&A[181]);
#if defined(__SSE3__) && defined(__AVX__)
    c53_1 = _mm_add_sd(c53_1, _mm_mul_sd(a53_1, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c53_1 = _mm_add_sd(c53_1, _mm_mul_sd(a53_1, b53));
#endif
    _mm_store_sd(&C[(l_n*56)+7], c53_1);
    __m128d c53_2 = _mm_load_sd(&C[(l_n*56)+17]);
    __m128d a53_2 = _mm_load_sd(&A[182]);
#if defined(__SSE3__) && defined(__AVX__)
    c53_2 = _mm_add_sd(c53_2, _mm_mul_sd(a53_2, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c53_2 = _mm_add_sd(c53_2, _mm_mul_sd(a53_2, b53));
#endif
    _mm_store_sd(&C[(l_n*56)+17], c53_2);
    __m128d c53_3 = _mm_load_sd(&C[(l_n*56)+32]);
    __m128d a53_3 = _mm_load_sd(&A[183]);
#if defined(__SSE3__) && defined(__AVX__)
    c53_3 = _mm_add_sd(c53_3, _mm_mul_sd(a53_3, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c53_3 = _mm_add_sd(c53_3, _mm_mul_sd(a53_3, b53));
#endif
    _mm_store_sd(&C[(l_n*56)+32], c53_3);
    __m128d c53_4 = _mm_load_sd(&C[(l_n*56)+53]);
    __m128d a53_4 = _mm_load_sd(&A[184]);
#if defined(__SSE3__) && defined(__AVX__)
    c53_4 = _mm_add_sd(c53_4, _mm_mul_sd(a53_4, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c53_4 = _mm_add_sd(c53_4, _mm_mul_sd(a53_4, b53));
#endif
    _mm_store_sd(&C[(l_n*56)+53], c53_4);
#else
    C[(l_n*56)+1] += A[180] * B[(l_n*56)+53];
    C[(l_n*56)+7] += A[181] * B[(l_n*56)+53];
    C[(l_n*56)+17] += A[182] * B[(l_n*56)+53];
    C[(l_n*56)+32] += A[183] * B[(l_n*56)+53];
    C[(l_n*56)+53] += A[184] * B[(l_n*56)+53];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b54 = _mm256_broadcast_sd(&B[(l_n*56)+54]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b54 = _mm_loaddup_pd(&B[(l_n*56)+54]);
#endif
    __m128d c54_0 = _mm_load_sd(&C[(l_n*56)+2]);
    __m128d a54_0 = _mm_load_sd(&A[185]);
#if defined(__SSE3__) && defined(__AVX__)
    c54_0 = _mm_add_sd(c54_0, _mm_mul_sd(a54_0, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c54_0 = _mm_add_sd(c54_0, _mm_mul_sd(a54_0, b54));
#endif
    _mm_store_sd(&C[(l_n*56)+2], c54_0);
    __m128d c54_1 = _mm_load_sd(&C[(l_n*56)+8]);
    __m128d a54_1 = _mm_load_sd(&A[186]);
#if defined(__SSE3__) && defined(__AVX__)
    c54_1 = _mm_add_sd(c54_1, _mm_mul_sd(a54_1, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c54_1 = _mm_add_sd(c54_1, _mm_mul_sd(a54_1, b54));
#endif
    _mm_store_sd(&C[(l_n*56)+8], c54_1);
    __m128d c54_2 = _mm_load_sd(&C[(l_n*56)+18]);
    __m128d a54_2 = _mm_load_sd(&A[187]);
#if defined(__SSE3__) && defined(__AVX__)
    c54_2 = _mm_add_sd(c54_2, _mm_mul_sd(a54_2, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c54_2 = _mm_add_sd(c54_2, _mm_mul_sd(a54_2, b54));
#endif
    _mm_store_sd(&C[(l_n*56)+18], c54_2);
    __m128d c54_3 = _mm_load_sd(&C[(l_n*56)+33]);
    __m128d a54_3 = _mm_load_sd(&A[188]);
#if defined(__SSE3__) && defined(__AVX__)
    c54_3 = _mm_add_sd(c54_3, _mm_mul_sd(a54_3, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c54_3 = _mm_add_sd(c54_3, _mm_mul_sd(a54_3, b54));
#endif
    _mm_store_sd(&C[(l_n*56)+33], c54_3);
    __m128d c54_4 = _mm_load_sd(&C[(l_n*56)+54]);
    __m128d a54_4 = _mm_load_sd(&A[189]);
#if defined(__SSE3__) && defined(__AVX__)
    c54_4 = _mm_add_sd(c54_4, _mm_mul_sd(a54_4, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c54_4 = _mm_add_sd(c54_4, _mm_mul_sd(a54_4, b54));
#endif
    _mm_store_sd(&C[(l_n*56)+54], c54_4);
#else
    C[(l_n*56)+2] += A[185] * B[(l_n*56)+54];
    C[(l_n*56)+8] += A[186] * B[(l_n*56)+54];
    C[(l_n*56)+18] += A[187] * B[(l_n*56)+54];
    C[(l_n*56)+33] += A[188] * B[(l_n*56)+54];
    C[(l_n*56)+54] += A[189] * B[(l_n*56)+54];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b55 = _mm256_broadcast_sd(&B[(l_n*56)+55]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b55 = _mm_loaddup_pd(&B[(l_n*56)+55]);
#endif
    __m128d c55_0 = _mm_load_sd(&C[(l_n*56)+0]);
    __m128d a55_0 = _mm_load_sd(&A[190]);
#if defined(__SSE3__) && defined(__AVX__)
    c55_0 = _mm_add_sd(c55_0, _mm_mul_sd(a55_0, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c55_0 = _mm_add_sd(c55_0, _mm_mul_sd(a55_0, b55));
#endif
    _mm_store_sd(&C[(l_n*56)+0], c55_0);
    __m128d c55_1 = _mm_load_sd(&C[(l_n*56)+3]);
    __m128d a55_1 = _mm_load_sd(&A[191]);
#if defined(__SSE3__) && defined(__AVX__)
    c55_1 = _mm_add_sd(c55_1, _mm_mul_sd(a55_1, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c55_1 = _mm_add_sd(c55_1, _mm_mul_sd(a55_1, b55));
#endif
    _mm_store_sd(&C[(l_n*56)+3], c55_1);
    __m128d c55_2 = _mm_load_sd(&C[(l_n*56)+9]);
    __m128d a55_2 = _mm_load_sd(&A[192]);
#if defined(__SSE3__) && defined(__AVX__)
    c55_2 = _mm_add_sd(c55_2, _mm_mul_sd(a55_2, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c55_2 = _mm_add_sd(c55_2, _mm_mul_sd(a55_2, b55));
#endif
    _mm_store_sd(&C[(l_n*56)+9], c55_2);
    __m128d c55_3 = _mm_load_sd(&C[(l_n*56)+19]);
    __m128d a55_3 = _mm_load_sd(&A[193]);
#if defined(__SSE3__) && defined(__AVX__)
    c55_3 = _mm_add_sd(c55_3, _mm_mul_sd(a55_3, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c55_3 = _mm_add_sd(c55_3, _mm_mul_sd(a55_3, b55));
#endif
    _mm_store_sd(&C[(l_n*56)+19], c55_3);
    __m128d c55_4 = _mm_load_sd(&C[(l_n*56)+34]);
    __m128d a55_4 = _mm_load_sd(&A[194]);
#if defined(__SSE3__) && defined(__AVX__)
    c55_4 = _mm_add_sd(c55_4, _mm_mul_sd(a55_4, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c55_4 = _mm_add_sd(c55_4, _mm_mul_sd(a55_4, b55));
#endif
    _mm_store_sd(&C[(l_n*56)+34], c55_4);
    __m128d c55_5 = _mm_load_sd(&C[(l_n*56)+55]);
    __m128d a55_5 = _mm_load_sd(&A[195]);
#if defined(__SSE3__) && defined(__AVX__)
    c55_5 = _mm_add_sd(c55_5, _mm_mul_sd(a55_5, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c55_5 = _mm_add_sd(c55_5, _mm_mul_sd(a55_5, b55));
#endif
    _mm_store_sd(&C[(l_n*56)+55], c55_5);
#else
    C[(l_n*56)+0] += A[190] * B[(l_n*56)+55];
    C[(l_n*56)+3] += A[191] * B[(l_n*56)+55];
    C[(l_n*56)+9] += A[192] * B[(l_n*56)+55];
    C[(l_n*56)+19] += A[193] * B[(l_n*56)+55];
    C[(l_n*56)+34] += A[194] * B[(l_n*56)+55];
    C[(l_n*56)+55] += A[195] * B[(l_n*56)+55];
#endif

  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 3528;
#endif
}

void dsparse_starMatrix_m84_n9_k9_ldA84_ldBna7_ldC84_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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

void dsparse_fM1DivM_m84_n9_k84_ldAna7_ldB84_ldC84_beta0_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_n = 0;
  #pragma nounroll_and_jam
  for ( l_n = 0; l_n < 9; l_n++) {
    unsigned int l_m = 0;
   #pragma simd
    for ( l_m = 0; l_m < 84; l_m++) {
      C[(l_n*84)+l_m] = 0.0;
    }
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b0 = _mm256_broadcast_sd(&B[(l_n*84)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b0 = _mm_loaddup_pd(&B[(l_n*84)+0]);
#endif
    __m128d c0_0 = _mm_load_sd(&C[(l_n*84)+0]);
    __m128d a0_0 = _mm_load_sd(&A[0]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, b0));
#endif
    _mm_store_sd(&C[(l_n*84)+0], c0_0);
    __m128d c0_1 = _mm_load_sd(&C[(l_n*84)+3]);
    __m128d a0_1 = _mm_load_sd(&A[1]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, b0));
#endif
    _mm_store_sd(&C[(l_n*84)+3], c0_1);
    __m128d c0_2 = _mm_load_sd(&C[(l_n*84)+9]);
    __m128d a0_2 = _mm_load_sd(&A[2]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, b0));
#endif
    _mm_store_sd(&C[(l_n*84)+9], c0_2);
    __m128d c0_3 = _mm_load_sd(&C[(l_n*84)+19]);
    __m128d a0_3 = _mm_load_sd(&A[3]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_3 = _mm_add_sd(c0_3, _mm_mul_sd(a0_3, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_3 = _mm_add_sd(c0_3, _mm_mul_sd(a0_3, b0));
#endif
    _mm_store_sd(&C[(l_n*84)+19], c0_3);
    __m128d c0_4 = _mm_load_sd(&C[(l_n*84)+34]);
    __m128d a0_4 = _mm_load_sd(&A[4]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_4 = _mm_add_sd(c0_4, _mm_mul_sd(a0_4, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_4 = _mm_add_sd(c0_4, _mm_mul_sd(a0_4, b0));
#endif
    _mm_store_sd(&C[(l_n*84)+34], c0_4);
    __m128d c0_5 = _mm_load_sd(&C[(l_n*84)+55]);
    __m128d a0_5 = _mm_load_sd(&A[5]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_5 = _mm_add_sd(c0_5, _mm_mul_sd(a0_5, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_5 = _mm_add_sd(c0_5, _mm_mul_sd(a0_5, b0));
#endif
    _mm_store_sd(&C[(l_n*84)+55], c0_5);
    __m128d c0_6 = _mm_load_sd(&C[(l_n*84)+83]);
    __m128d a0_6 = _mm_load_sd(&A[6]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_6 = _mm_add_sd(c0_6, _mm_mul_sd(a0_6, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_6 = _mm_add_sd(c0_6, _mm_mul_sd(a0_6, b0));
#endif
    _mm_store_sd(&C[(l_n*84)+83], c0_6);
#else
    C[(l_n*84)+0] += A[0] * B[(l_n*84)+0];
    C[(l_n*84)+3] += A[1] * B[(l_n*84)+0];
    C[(l_n*84)+9] += A[2] * B[(l_n*84)+0];
    C[(l_n*84)+19] += A[3] * B[(l_n*84)+0];
    C[(l_n*84)+34] += A[4] * B[(l_n*84)+0];
    C[(l_n*84)+55] += A[5] * B[(l_n*84)+0];
    C[(l_n*84)+83] += A[6] * B[(l_n*84)+0];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b1 = _mm256_broadcast_sd(&B[(l_n*84)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b1 = _mm_loaddup_pd(&B[(l_n*84)+1]);
#endif
    __m128d c1_0 = _mm_load_sd(&C[(l_n*84)+1]);
    __m128d a1_0 = _mm_load_sd(&A[7]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, b1));
#endif
    _mm_store_sd(&C[(l_n*84)+1], c1_0);
    __m128d c1_1 = _mm_load_sd(&C[(l_n*84)+7]);
    __m128d a1_1 = _mm_load_sd(&A[8]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_1 = _mm_add_sd(c1_1, _mm_mul_sd(a1_1, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_1 = _mm_add_sd(c1_1, _mm_mul_sd(a1_1, b1));
#endif
    _mm_store_sd(&C[(l_n*84)+7], c1_1);
    __m128d c1_2 = _mm_load_sd(&C[(l_n*84)+17]);
    __m128d a1_2 = _mm_load_sd(&A[9]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_2 = _mm_add_sd(c1_2, _mm_mul_sd(a1_2, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_2 = _mm_add_sd(c1_2, _mm_mul_sd(a1_2, b1));
#endif
    _mm_store_sd(&C[(l_n*84)+17], c1_2);
    __m128d c1_3 = _mm_load_sd(&C[(l_n*84)+32]);
    __m128d a1_3 = _mm_load_sd(&A[10]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_3 = _mm_add_sd(c1_3, _mm_mul_sd(a1_3, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_3 = _mm_add_sd(c1_3, _mm_mul_sd(a1_3, b1));
#endif
    _mm_store_sd(&C[(l_n*84)+32], c1_3);
    __m128d c1_4 = _mm_load_sd(&C[(l_n*84)+53]);
    __m128d a1_4 = _mm_load_sd(&A[11]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_4 = _mm_add_sd(c1_4, _mm_mul_sd(a1_4, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_4 = _mm_add_sd(c1_4, _mm_mul_sd(a1_4, b1));
#endif
    _mm_store_sd(&C[(l_n*84)+53], c1_4);
    __m128d c1_5 = _mm_load_sd(&C[(l_n*84)+81]);
    __m128d a1_5 = _mm_load_sd(&A[12]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_5 = _mm_add_sd(c1_5, _mm_mul_sd(a1_5, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_5 = _mm_add_sd(c1_5, _mm_mul_sd(a1_5, b1));
#endif
    _mm_store_sd(&C[(l_n*84)+81], c1_5);
#else
    C[(l_n*84)+1] += A[7] * B[(l_n*84)+1];
    C[(l_n*84)+7] += A[8] * B[(l_n*84)+1];
    C[(l_n*84)+17] += A[9] * B[(l_n*84)+1];
    C[(l_n*84)+32] += A[10] * B[(l_n*84)+1];
    C[(l_n*84)+53] += A[11] * B[(l_n*84)+1];
    C[(l_n*84)+81] += A[12] * B[(l_n*84)+1];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b2 = _mm256_broadcast_sd(&B[(l_n*84)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b2 = _mm_loaddup_pd(&B[(l_n*84)+2]);
#endif
    __m128d c2_0 = _mm_load_sd(&C[(l_n*84)+2]);
    __m128d a2_0 = _mm_load_sd(&A[13]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_0 = _mm_add_sd(c2_0, _mm_mul_sd(a2_0, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_0 = _mm_add_sd(c2_0, _mm_mul_sd(a2_0, b2));
#endif
    _mm_store_sd(&C[(l_n*84)+2], c2_0);
    __m128d c2_1 = _mm_load_sd(&C[(l_n*84)+8]);
    __m128d a2_1 = _mm_load_sd(&A[14]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_1 = _mm_add_sd(c2_1, _mm_mul_sd(a2_1, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_1 = _mm_add_sd(c2_1, _mm_mul_sd(a2_1, b2));
#endif
    _mm_store_sd(&C[(l_n*84)+8], c2_1);
    __m128d c2_2 = _mm_load_sd(&C[(l_n*84)+18]);
    __m128d a2_2 = _mm_load_sd(&A[15]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_2 = _mm_add_sd(c2_2, _mm_mul_sd(a2_2, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_2 = _mm_add_sd(c2_2, _mm_mul_sd(a2_2, b2));
#endif
    _mm_store_sd(&C[(l_n*84)+18], c2_2);
    __m128d c2_3 = _mm_load_sd(&C[(l_n*84)+33]);
    __m128d a2_3 = _mm_load_sd(&A[16]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_3 = _mm_add_sd(c2_3, _mm_mul_sd(a2_3, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_3 = _mm_add_sd(c2_3, _mm_mul_sd(a2_3, b2));
#endif
    _mm_store_sd(&C[(l_n*84)+33], c2_3);
    __m128d c2_4 = _mm_load_sd(&C[(l_n*84)+54]);
    __m128d a2_4 = _mm_load_sd(&A[17]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_4 = _mm_add_sd(c2_4, _mm_mul_sd(a2_4, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_4 = _mm_add_sd(c2_4, _mm_mul_sd(a2_4, b2));
#endif
    _mm_store_sd(&C[(l_n*84)+54], c2_4);
    __m128d c2_5 = _mm_load_sd(&C[(l_n*84)+82]);
    __m128d a2_5 = _mm_load_sd(&A[18]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_5 = _mm_add_sd(c2_5, _mm_mul_sd(a2_5, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_5 = _mm_add_sd(c2_5, _mm_mul_sd(a2_5, b2));
#endif
    _mm_store_sd(&C[(l_n*84)+82], c2_5);
#else
    C[(l_n*84)+2] += A[13] * B[(l_n*84)+2];
    C[(l_n*84)+8] += A[14] * B[(l_n*84)+2];
    C[(l_n*84)+18] += A[15] * B[(l_n*84)+2];
    C[(l_n*84)+33] += A[16] * B[(l_n*84)+2];
    C[(l_n*84)+54] += A[17] * B[(l_n*84)+2];
    C[(l_n*84)+82] += A[18] * B[(l_n*84)+2];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b3 = _mm256_broadcast_sd(&B[(l_n*84)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b3 = _mm_loaddup_pd(&B[(l_n*84)+3]);
#endif
    __m128d c3_0 = _mm_load_sd(&C[(l_n*84)+0]);
    __m128d a3_0 = _mm_load_sd(&A[19]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, b3));
#endif
    _mm_store_sd(&C[(l_n*84)+0], c3_0);
    __m128d c3_1 = _mm_load_sd(&C[(l_n*84)+3]);
    __m128d a3_1 = _mm_load_sd(&A[20]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, b3));
#endif
    _mm_store_sd(&C[(l_n*84)+3], c3_1);
    __m128d c3_2 = _mm_load_sd(&C[(l_n*84)+9]);
    __m128d a3_2 = _mm_load_sd(&A[21]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, b3));
#endif
    _mm_store_sd(&C[(l_n*84)+9], c3_2);
    __m128d c3_3 = _mm_load_sd(&C[(l_n*84)+19]);
    __m128d a3_3 = _mm_load_sd(&A[22]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_3 = _mm_add_sd(c3_3, _mm_mul_sd(a3_3, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_3 = _mm_add_sd(c3_3, _mm_mul_sd(a3_3, b3));
#endif
    _mm_store_sd(&C[(l_n*84)+19], c3_3);
    __m128d c3_4 = _mm_load_sd(&C[(l_n*84)+34]);
    __m128d a3_4 = _mm_load_sd(&A[23]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_4 = _mm_add_sd(c3_4, _mm_mul_sd(a3_4, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_4 = _mm_add_sd(c3_4, _mm_mul_sd(a3_4, b3));
#endif
    _mm_store_sd(&C[(l_n*84)+34], c3_4);
    __m128d c3_5 = _mm_load_sd(&C[(l_n*84)+55]);
    __m128d a3_5 = _mm_load_sd(&A[24]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_5 = _mm_add_sd(c3_5, _mm_mul_sd(a3_5, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_5 = _mm_add_sd(c3_5, _mm_mul_sd(a3_5, b3));
#endif
    _mm_store_sd(&C[(l_n*84)+55], c3_5);
    __m128d c3_6 = _mm_load_sd(&C[(l_n*84)+83]);
    __m128d a3_6 = _mm_load_sd(&A[25]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_6 = _mm_add_sd(c3_6, _mm_mul_sd(a3_6, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_6 = _mm_add_sd(c3_6, _mm_mul_sd(a3_6, b3));
#endif
    _mm_store_sd(&C[(l_n*84)+83], c3_6);
#else
    C[(l_n*84)+0] += A[19] * B[(l_n*84)+3];
    C[(l_n*84)+3] += A[20] * B[(l_n*84)+3];
    C[(l_n*84)+9] += A[21] * B[(l_n*84)+3];
    C[(l_n*84)+19] += A[22] * B[(l_n*84)+3];
    C[(l_n*84)+34] += A[23] * B[(l_n*84)+3];
    C[(l_n*84)+55] += A[24] * B[(l_n*84)+3];
    C[(l_n*84)+83] += A[25] * B[(l_n*84)+3];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b4 = _mm256_broadcast_sd(&B[(l_n*84)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b4 = _mm_loaddup_pd(&B[(l_n*84)+4]);
#endif
    __m128d c4_0 = _mm_load_sd(&C[(l_n*84)+4]);
    __m128d a4_0 = _mm_load_sd(&A[26]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_0 = _mm_add_sd(c4_0, _mm_mul_sd(a4_0, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_0 = _mm_add_sd(c4_0, _mm_mul_sd(a4_0, b4));
#endif
    _mm_store_sd(&C[(l_n*84)+4], c4_0);
    __m128d c4_1 = _mm_load_sd(&C[(l_n*84)+14]);
    __m128d a4_1 = _mm_load_sd(&A[27]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_1 = _mm_add_sd(c4_1, _mm_mul_sd(a4_1, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_1 = _mm_add_sd(c4_1, _mm_mul_sd(a4_1, b4));
#endif
    _mm_store_sd(&C[(l_n*84)+14], c4_1);
    __m128d c4_2 = _mm_load_sd(&C[(l_n*84)+29]);
    __m128d a4_2 = _mm_load_sd(&A[28]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, b4));
#endif
    _mm_store_sd(&C[(l_n*84)+29], c4_2);
    __m128d c4_3 = _mm_load_sd(&C[(l_n*84)+50]);
    __m128d a4_3 = _mm_load_sd(&A[29]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_3 = _mm_add_sd(c4_3, _mm_mul_sd(a4_3, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_3 = _mm_add_sd(c4_3, _mm_mul_sd(a4_3, b4));
#endif
    _mm_store_sd(&C[(l_n*84)+50], c4_3);
    __m128d c4_4 = _mm_load_sd(&C[(l_n*84)+78]);
    __m128d a4_4 = _mm_load_sd(&A[30]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_4 = _mm_add_sd(c4_4, _mm_mul_sd(a4_4, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_4 = _mm_add_sd(c4_4, _mm_mul_sd(a4_4, b4));
#endif
    _mm_store_sd(&C[(l_n*84)+78], c4_4);
#else
    C[(l_n*84)+4] += A[26] * B[(l_n*84)+4];
    C[(l_n*84)+14] += A[27] * B[(l_n*84)+4];
    C[(l_n*84)+29] += A[28] * B[(l_n*84)+4];
    C[(l_n*84)+50] += A[29] * B[(l_n*84)+4];
    C[(l_n*84)+78] += A[30] * B[(l_n*84)+4];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b5 = _mm256_broadcast_sd(&B[(l_n*84)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b5 = _mm_loaddup_pd(&B[(l_n*84)+5]);
#endif
    __m128d c5_0 = _mm_load_sd(&C[(l_n*84)+5]);
    __m128d a5_0 = _mm_load_sd(&A[31]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_0 = _mm_add_sd(c5_0, _mm_mul_sd(a5_0, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_0 = _mm_add_sd(c5_0, _mm_mul_sd(a5_0, b5));
#endif
    _mm_store_sd(&C[(l_n*84)+5], c5_0);
    __m128d c5_1 = _mm_load_sd(&C[(l_n*84)+15]);
    __m128d a5_1 = _mm_load_sd(&A[32]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_1 = _mm_add_sd(c5_1, _mm_mul_sd(a5_1, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_1 = _mm_add_sd(c5_1, _mm_mul_sd(a5_1, b5));
#endif
    _mm_store_sd(&C[(l_n*84)+15], c5_1);
    __m128d c5_2 = _mm_load_sd(&C[(l_n*84)+30]);
    __m128d a5_2 = _mm_load_sd(&A[33]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, b5));
#endif
    _mm_store_sd(&C[(l_n*84)+30], c5_2);
    __m128d c5_3 = _mm_load_sd(&C[(l_n*84)+51]);
    __m128d a5_3 = _mm_load_sd(&A[34]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_3 = _mm_add_sd(c5_3, _mm_mul_sd(a5_3, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_3 = _mm_add_sd(c5_3, _mm_mul_sd(a5_3, b5));
#endif
    _mm_store_sd(&C[(l_n*84)+51], c5_3);
    __m128d c5_4 = _mm_load_sd(&C[(l_n*84)+79]);
    __m128d a5_4 = _mm_load_sd(&A[35]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_4 = _mm_add_sd(c5_4, _mm_mul_sd(a5_4, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_4 = _mm_add_sd(c5_4, _mm_mul_sd(a5_4, b5));
#endif
    _mm_store_sd(&C[(l_n*84)+79], c5_4);
#else
    C[(l_n*84)+5] += A[31] * B[(l_n*84)+5];
    C[(l_n*84)+15] += A[32] * B[(l_n*84)+5];
    C[(l_n*84)+30] += A[33] * B[(l_n*84)+5];
    C[(l_n*84)+51] += A[34] * B[(l_n*84)+5];
    C[(l_n*84)+79] += A[35] * B[(l_n*84)+5];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b6 = _mm256_broadcast_sd(&B[(l_n*84)+6]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b6 = _mm_loaddup_pd(&B[(l_n*84)+6]);
#endif
    __m128d c6_0 = _mm_load_sd(&C[(l_n*84)+6]);
    __m128d a6_0 = _mm_load_sd(&A[36]);
#if defined(__SSE3__) && defined(__AVX__)
    c6_0 = _mm_add_sd(c6_0, _mm_mul_sd(a6_0, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c6_0 = _mm_add_sd(c6_0, _mm_mul_sd(a6_0, b6));
#endif
    _mm_store_sd(&C[(l_n*84)+6], c6_0);
    __m128d c6_1 = _mm_load_sd(&C[(l_n*84)+16]);
    __m128d a6_1 = _mm_load_sd(&A[37]);
#if defined(__SSE3__) && defined(__AVX__)
    c6_1 = _mm_add_sd(c6_1, _mm_mul_sd(a6_1, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c6_1 = _mm_add_sd(c6_1, _mm_mul_sd(a6_1, b6));
#endif
    _mm_store_sd(&C[(l_n*84)+16], c6_1);
    __m128d c6_2 = _mm_load_sd(&C[(l_n*84)+31]);
    __m128d a6_2 = _mm_load_sd(&A[38]);
#if defined(__SSE3__) && defined(__AVX__)
    c6_2 = _mm_add_sd(c6_2, _mm_mul_sd(a6_2, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c6_2 = _mm_add_sd(c6_2, _mm_mul_sd(a6_2, b6));
#endif
    _mm_store_sd(&C[(l_n*84)+31], c6_2);
    __m128d c6_3 = _mm_load_sd(&C[(l_n*84)+52]);
    __m128d a6_3 = _mm_load_sd(&A[39]);
#if defined(__SSE3__) && defined(__AVX__)
    c6_3 = _mm_add_sd(c6_3, _mm_mul_sd(a6_3, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c6_3 = _mm_add_sd(c6_3, _mm_mul_sd(a6_3, b6));
#endif
    _mm_store_sd(&C[(l_n*84)+52], c6_3);
    __m128d c6_4 = _mm_load_sd(&C[(l_n*84)+80]);
    __m128d a6_4 = _mm_load_sd(&A[40]);
#if defined(__SSE3__) && defined(__AVX__)
    c6_4 = _mm_add_sd(c6_4, _mm_mul_sd(a6_4, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c6_4 = _mm_add_sd(c6_4, _mm_mul_sd(a6_4, b6));
#endif
    _mm_store_sd(&C[(l_n*84)+80], c6_4);
#else
    C[(l_n*84)+6] += A[36] * B[(l_n*84)+6];
    C[(l_n*84)+16] += A[37] * B[(l_n*84)+6];
    C[(l_n*84)+31] += A[38] * B[(l_n*84)+6];
    C[(l_n*84)+52] += A[39] * B[(l_n*84)+6];
    C[(l_n*84)+80] += A[40] * B[(l_n*84)+6];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b7 = _mm256_broadcast_sd(&B[(l_n*84)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b7 = _mm_loaddup_pd(&B[(l_n*84)+7]);
#endif
    __m128d c7_0 = _mm_load_sd(&C[(l_n*84)+1]);
    __m128d a7_0 = _mm_load_sd(&A[41]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_0 = _mm_add_sd(c7_0, _mm_mul_sd(a7_0, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_0 = _mm_add_sd(c7_0, _mm_mul_sd(a7_0, b7));
#endif
    _mm_store_sd(&C[(l_n*84)+1], c7_0);
    __m128d c7_1 = _mm_load_sd(&C[(l_n*84)+7]);
    __m128d a7_1 = _mm_load_sd(&A[42]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_1 = _mm_add_sd(c7_1, _mm_mul_sd(a7_1, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_1 = _mm_add_sd(c7_1, _mm_mul_sd(a7_1, b7));
#endif
    _mm_store_sd(&C[(l_n*84)+7], c7_1);
    __m128d c7_2 = _mm_load_sd(&C[(l_n*84)+17]);
    __m128d a7_2 = _mm_load_sd(&A[43]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_2 = _mm_add_sd(c7_2, _mm_mul_sd(a7_2, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_2 = _mm_add_sd(c7_2, _mm_mul_sd(a7_2, b7));
#endif
    _mm_store_sd(&C[(l_n*84)+17], c7_2);
    __m128d c7_3 = _mm_load_sd(&C[(l_n*84)+32]);
    __m128d a7_3 = _mm_load_sd(&A[44]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_3 = _mm_add_sd(c7_3, _mm_mul_sd(a7_3, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_3 = _mm_add_sd(c7_3, _mm_mul_sd(a7_3, b7));
#endif
    _mm_store_sd(&C[(l_n*84)+32], c7_3);
    __m128d c7_4 = _mm_load_sd(&C[(l_n*84)+53]);
    __m128d a7_4 = _mm_load_sd(&A[45]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_4 = _mm_add_sd(c7_4, _mm_mul_sd(a7_4, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_4 = _mm_add_sd(c7_4, _mm_mul_sd(a7_4, b7));
#endif
    _mm_store_sd(&C[(l_n*84)+53], c7_4);
    __m128d c7_5 = _mm_load_sd(&C[(l_n*84)+81]);
    __m128d a7_5 = _mm_load_sd(&A[46]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_5 = _mm_add_sd(c7_5, _mm_mul_sd(a7_5, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_5 = _mm_add_sd(c7_5, _mm_mul_sd(a7_5, b7));
#endif
    _mm_store_sd(&C[(l_n*84)+81], c7_5);
#else
    C[(l_n*84)+1] += A[41] * B[(l_n*84)+7];
    C[(l_n*84)+7] += A[42] * B[(l_n*84)+7];
    C[(l_n*84)+17] += A[43] * B[(l_n*84)+7];
    C[(l_n*84)+32] += A[44] * B[(l_n*84)+7];
    C[(l_n*84)+53] += A[45] * B[(l_n*84)+7];
    C[(l_n*84)+81] += A[46] * B[(l_n*84)+7];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b8 = _mm256_broadcast_sd(&B[(l_n*84)+8]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b8 = _mm_loaddup_pd(&B[(l_n*84)+8]);
#endif
    __m128d c8_0 = _mm_load_sd(&C[(l_n*84)+2]);
    __m128d a8_0 = _mm_load_sd(&A[47]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_0 = _mm_add_sd(c8_0, _mm_mul_sd(a8_0, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_0 = _mm_add_sd(c8_0, _mm_mul_sd(a8_0, b8));
#endif
    _mm_store_sd(&C[(l_n*84)+2], c8_0);
    __m128d c8_1 = _mm_load_sd(&C[(l_n*84)+8]);
    __m128d a8_1 = _mm_load_sd(&A[48]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_1 = _mm_add_sd(c8_1, _mm_mul_sd(a8_1, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_1 = _mm_add_sd(c8_1, _mm_mul_sd(a8_1, b8));
#endif
    _mm_store_sd(&C[(l_n*84)+8], c8_1);
    __m128d c8_2 = _mm_load_sd(&C[(l_n*84)+18]);
    __m128d a8_2 = _mm_load_sd(&A[49]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_2 = _mm_add_sd(c8_2, _mm_mul_sd(a8_2, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_2 = _mm_add_sd(c8_2, _mm_mul_sd(a8_2, b8));
#endif
    _mm_store_sd(&C[(l_n*84)+18], c8_2);
    __m128d c8_3 = _mm_load_sd(&C[(l_n*84)+33]);
    __m128d a8_3 = _mm_load_sd(&A[50]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_3 = _mm_add_sd(c8_3, _mm_mul_sd(a8_3, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_3 = _mm_add_sd(c8_3, _mm_mul_sd(a8_3, b8));
#endif
    _mm_store_sd(&C[(l_n*84)+33], c8_3);
    __m128d c8_4 = _mm_load_sd(&C[(l_n*84)+54]);
    __m128d a8_4 = _mm_load_sd(&A[51]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_4 = _mm_add_sd(c8_4, _mm_mul_sd(a8_4, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_4 = _mm_add_sd(c8_4, _mm_mul_sd(a8_4, b8));
#endif
    _mm_store_sd(&C[(l_n*84)+54], c8_4);
    __m128d c8_5 = _mm_load_sd(&C[(l_n*84)+82]);
    __m128d a8_5 = _mm_load_sd(&A[52]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_5 = _mm_add_sd(c8_5, _mm_mul_sd(a8_5, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_5 = _mm_add_sd(c8_5, _mm_mul_sd(a8_5, b8));
#endif
    _mm_store_sd(&C[(l_n*84)+82], c8_5);
#else
    C[(l_n*84)+2] += A[47] * B[(l_n*84)+8];
    C[(l_n*84)+8] += A[48] * B[(l_n*84)+8];
    C[(l_n*84)+18] += A[49] * B[(l_n*84)+8];
    C[(l_n*84)+33] += A[50] * B[(l_n*84)+8];
    C[(l_n*84)+54] += A[51] * B[(l_n*84)+8];
    C[(l_n*84)+82] += A[52] * B[(l_n*84)+8];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b9 = _mm256_broadcast_sd(&B[(l_n*84)+9]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b9 = _mm_loaddup_pd(&B[(l_n*84)+9]);
#endif
    __m128d c9_0 = _mm_load_sd(&C[(l_n*84)+0]);
    __m128d a9_0 = _mm_load_sd(&A[53]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, b9));
#endif
    _mm_store_sd(&C[(l_n*84)+0], c9_0);
    __m128d c9_1 = _mm_load_sd(&C[(l_n*84)+3]);
    __m128d a9_1 = _mm_load_sd(&A[54]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, b9));
#endif
    _mm_store_sd(&C[(l_n*84)+3], c9_1);
    __m128d c9_2 = _mm_load_sd(&C[(l_n*84)+9]);
    __m128d a9_2 = _mm_load_sd(&A[55]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, b9));
#endif
    _mm_store_sd(&C[(l_n*84)+9], c9_2);
    __m128d c9_3 = _mm_load_sd(&C[(l_n*84)+19]);
    __m128d a9_3 = _mm_load_sd(&A[56]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_3 = _mm_add_sd(c9_3, _mm_mul_sd(a9_3, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_3 = _mm_add_sd(c9_3, _mm_mul_sd(a9_3, b9));
#endif
    _mm_store_sd(&C[(l_n*84)+19], c9_3);
    __m128d c9_4 = _mm_load_sd(&C[(l_n*84)+34]);
    __m128d a9_4 = _mm_load_sd(&A[57]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_4 = _mm_add_sd(c9_4, _mm_mul_sd(a9_4, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_4 = _mm_add_sd(c9_4, _mm_mul_sd(a9_4, b9));
#endif
    _mm_store_sd(&C[(l_n*84)+34], c9_4);
    __m128d c9_5 = _mm_load_sd(&C[(l_n*84)+55]);
    __m128d a9_5 = _mm_load_sd(&A[58]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_5 = _mm_add_sd(c9_5, _mm_mul_sd(a9_5, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_5 = _mm_add_sd(c9_5, _mm_mul_sd(a9_5, b9));
#endif
    _mm_store_sd(&C[(l_n*84)+55], c9_5);
    __m128d c9_6 = _mm_load_sd(&C[(l_n*84)+83]);
    __m128d a9_6 = _mm_load_sd(&A[59]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_6 = _mm_add_sd(c9_6, _mm_mul_sd(a9_6, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_6 = _mm_add_sd(c9_6, _mm_mul_sd(a9_6, b9));
#endif
    _mm_store_sd(&C[(l_n*84)+83], c9_6);
#else
    C[(l_n*84)+0] += A[53] * B[(l_n*84)+9];
    C[(l_n*84)+3] += A[54] * B[(l_n*84)+9];
    C[(l_n*84)+9] += A[55] * B[(l_n*84)+9];
    C[(l_n*84)+19] += A[56] * B[(l_n*84)+9];
    C[(l_n*84)+34] += A[57] * B[(l_n*84)+9];
    C[(l_n*84)+55] += A[58] * B[(l_n*84)+9];
    C[(l_n*84)+83] += A[59] * B[(l_n*84)+9];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b10 = _mm256_broadcast_sd(&B[(l_n*84)+10]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b10 = _mm_loaddup_pd(&B[(l_n*84)+10]);
#endif
    __m128d c10_0 = _mm_load_sd(&C[(l_n*84)+10]);
    __m128d a10_0 = _mm_load_sd(&A[60]);
#if defined(__SSE3__) && defined(__AVX__)
    c10_0 = _mm_add_sd(c10_0, _mm_mul_sd(a10_0, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c10_0 = _mm_add_sd(c10_0, _mm_mul_sd(a10_0, b10));
#endif
    _mm_store_sd(&C[(l_n*84)+10], c10_0);
    __m128d c10_1 = _mm_load_sd(&C[(l_n*84)+25]);
    __m128d a10_1 = _mm_load_sd(&A[61]);
#if defined(__SSE3__) && defined(__AVX__)
    c10_1 = _mm_add_sd(c10_1, _mm_mul_sd(a10_1, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c10_1 = _mm_add_sd(c10_1, _mm_mul_sd(a10_1, b10));
#endif
    _mm_store_sd(&C[(l_n*84)+25], c10_1);
    __m128d c10_2 = _mm_load_sd(&C[(l_n*84)+46]);
    __m128d a10_2 = _mm_load_sd(&A[62]);
#if defined(__SSE3__) && defined(__AVX__)
    c10_2 = _mm_add_sd(c10_2, _mm_mul_sd(a10_2, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c10_2 = _mm_add_sd(c10_2, _mm_mul_sd(a10_2, b10));
#endif
    _mm_store_sd(&C[(l_n*84)+46], c10_2);
    __m128d c10_3 = _mm_load_sd(&C[(l_n*84)+74]);
    __m128d a10_3 = _mm_load_sd(&A[63]);
#if defined(__SSE3__) && defined(__AVX__)
    c10_3 = _mm_add_sd(c10_3, _mm_mul_sd(a10_3, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c10_3 = _mm_add_sd(c10_3, _mm_mul_sd(a10_3, b10));
#endif
    _mm_store_sd(&C[(l_n*84)+74], c10_3);
#else
    C[(l_n*84)+10] += A[60] * B[(l_n*84)+10];
    C[(l_n*84)+25] += A[61] * B[(l_n*84)+10];
    C[(l_n*84)+46] += A[62] * B[(l_n*84)+10];
    C[(l_n*84)+74] += A[63] * B[(l_n*84)+10];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b11 = _mm256_broadcast_sd(&B[(l_n*84)+11]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b11 = _mm_loaddup_pd(&B[(l_n*84)+11]);
#endif
    __m128d c11_0 = _mm_load_sd(&C[(l_n*84)+11]);
    __m128d a11_0 = _mm_load_sd(&A[64]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_0 = _mm_add_sd(c11_0, _mm_mul_sd(a11_0, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_0 = _mm_add_sd(c11_0, _mm_mul_sd(a11_0, b11));
#endif
    _mm_store_sd(&C[(l_n*84)+11], c11_0);
    __m128d c11_1 = _mm_load_sd(&C[(l_n*84)+26]);
    __m128d a11_1 = _mm_load_sd(&A[65]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_1 = _mm_add_sd(c11_1, _mm_mul_sd(a11_1, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_1 = _mm_add_sd(c11_1, _mm_mul_sd(a11_1, b11));
#endif
    _mm_store_sd(&C[(l_n*84)+26], c11_1);
    __m128d c11_2 = _mm_load_sd(&C[(l_n*84)+47]);
    __m128d a11_2 = _mm_load_sd(&A[66]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_2 = _mm_add_sd(c11_2, _mm_mul_sd(a11_2, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_2 = _mm_add_sd(c11_2, _mm_mul_sd(a11_2, b11));
#endif
    _mm_store_sd(&C[(l_n*84)+47], c11_2);
    __m128d c11_3 = _mm_load_sd(&C[(l_n*84)+75]);
    __m128d a11_3 = _mm_load_sd(&A[67]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_3 = _mm_add_sd(c11_3, _mm_mul_sd(a11_3, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_3 = _mm_add_sd(c11_3, _mm_mul_sd(a11_3, b11));
#endif
    _mm_store_sd(&C[(l_n*84)+75], c11_3);
#else
    C[(l_n*84)+11] += A[64] * B[(l_n*84)+11];
    C[(l_n*84)+26] += A[65] * B[(l_n*84)+11];
    C[(l_n*84)+47] += A[66] * B[(l_n*84)+11];
    C[(l_n*84)+75] += A[67] * B[(l_n*84)+11];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b12 = _mm256_broadcast_sd(&B[(l_n*84)+12]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b12 = _mm_loaddup_pd(&B[(l_n*84)+12]);
#endif
    __m128d c12_0 = _mm_load_sd(&C[(l_n*84)+12]);
    __m128d a12_0 = _mm_load_sd(&A[68]);
#if defined(__SSE3__) && defined(__AVX__)
    c12_0 = _mm_add_sd(c12_0, _mm_mul_sd(a12_0, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c12_0 = _mm_add_sd(c12_0, _mm_mul_sd(a12_0, b12));
#endif
    _mm_store_sd(&C[(l_n*84)+12], c12_0);
    __m128d c12_1 = _mm_load_sd(&C[(l_n*84)+27]);
    __m128d a12_1 = _mm_load_sd(&A[69]);
#if defined(__SSE3__) && defined(__AVX__)
    c12_1 = _mm_add_sd(c12_1, _mm_mul_sd(a12_1, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c12_1 = _mm_add_sd(c12_1, _mm_mul_sd(a12_1, b12));
#endif
    _mm_store_sd(&C[(l_n*84)+27], c12_1);
    __m128d c12_2 = _mm_load_sd(&C[(l_n*84)+48]);
    __m128d a12_2 = _mm_load_sd(&A[70]);
#if defined(__SSE3__) && defined(__AVX__)
    c12_2 = _mm_add_sd(c12_2, _mm_mul_sd(a12_2, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c12_2 = _mm_add_sd(c12_2, _mm_mul_sd(a12_2, b12));
#endif
    _mm_store_sd(&C[(l_n*84)+48], c12_2);
    __m128d c12_3 = _mm_load_sd(&C[(l_n*84)+76]);
    __m128d a12_3 = _mm_load_sd(&A[71]);
#if defined(__SSE3__) && defined(__AVX__)
    c12_3 = _mm_add_sd(c12_3, _mm_mul_sd(a12_3, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c12_3 = _mm_add_sd(c12_3, _mm_mul_sd(a12_3, b12));
#endif
    _mm_store_sd(&C[(l_n*84)+76], c12_3);
#else
    C[(l_n*84)+12] += A[68] * B[(l_n*84)+12];
    C[(l_n*84)+27] += A[69] * B[(l_n*84)+12];
    C[(l_n*84)+48] += A[70] * B[(l_n*84)+12];
    C[(l_n*84)+76] += A[71] * B[(l_n*84)+12];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b13 = _mm256_broadcast_sd(&B[(l_n*84)+13]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b13 = _mm_loaddup_pd(&B[(l_n*84)+13]);
#endif
    __m128d c13_0 = _mm_load_sd(&C[(l_n*84)+13]);
    __m128d a13_0 = _mm_load_sd(&A[72]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_0 = _mm_add_sd(c13_0, _mm_mul_sd(a13_0, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_0 = _mm_add_sd(c13_0, _mm_mul_sd(a13_0, b13));
#endif
    _mm_store_sd(&C[(l_n*84)+13], c13_0);
    __m128d c13_1 = _mm_load_sd(&C[(l_n*84)+28]);
    __m128d a13_1 = _mm_load_sd(&A[73]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_1 = _mm_add_sd(c13_1, _mm_mul_sd(a13_1, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_1 = _mm_add_sd(c13_1, _mm_mul_sd(a13_1, b13));
#endif
    _mm_store_sd(&C[(l_n*84)+28], c13_1);
    __m128d c13_2 = _mm_load_sd(&C[(l_n*84)+49]);
    __m128d a13_2 = _mm_load_sd(&A[74]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_2 = _mm_add_sd(c13_2, _mm_mul_sd(a13_2, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_2 = _mm_add_sd(c13_2, _mm_mul_sd(a13_2, b13));
#endif
    _mm_store_sd(&C[(l_n*84)+49], c13_2);
    __m128d c13_3 = _mm_load_sd(&C[(l_n*84)+77]);
    __m128d a13_3 = _mm_load_sd(&A[75]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_3 = _mm_add_sd(c13_3, _mm_mul_sd(a13_3, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_3 = _mm_add_sd(c13_3, _mm_mul_sd(a13_3, b13));
#endif
    _mm_store_sd(&C[(l_n*84)+77], c13_3);
#else
    C[(l_n*84)+13] += A[72] * B[(l_n*84)+13];
    C[(l_n*84)+28] += A[73] * B[(l_n*84)+13];
    C[(l_n*84)+49] += A[74] * B[(l_n*84)+13];
    C[(l_n*84)+77] += A[75] * B[(l_n*84)+13];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b14 = _mm256_broadcast_sd(&B[(l_n*84)+14]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b14 = _mm_loaddup_pd(&B[(l_n*84)+14]);
#endif
    __m128d c14_0 = _mm_load_sd(&C[(l_n*84)+4]);
    __m128d a14_0 = _mm_load_sd(&A[76]);
#if defined(__SSE3__) && defined(__AVX__)
    c14_0 = _mm_add_sd(c14_0, _mm_mul_sd(a14_0, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c14_0 = _mm_add_sd(c14_0, _mm_mul_sd(a14_0, b14));
#endif
    _mm_store_sd(&C[(l_n*84)+4], c14_0);
    __m128d c14_1 = _mm_load_sd(&C[(l_n*84)+14]);
    __m128d a14_1 = _mm_load_sd(&A[77]);
#if defined(__SSE3__) && defined(__AVX__)
    c14_1 = _mm_add_sd(c14_1, _mm_mul_sd(a14_1, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c14_1 = _mm_add_sd(c14_1, _mm_mul_sd(a14_1, b14));
#endif
    _mm_store_sd(&C[(l_n*84)+14], c14_1);
    __m128d c14_2 = _mm_load_sd(&C[(l_n*84)+29]);
    __m128d a14_2 = _mm_load_sd(&A[78]);
#if defined(__SSE3__) && defined(__AVX__)
    c14_2 = _mm_add_sd(c14_2, _mm_mul_sd(a14_2, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c14_2 = _mm_add_sd(c14_2, _mm_mul_sd(a14_2, b14));
#endif
    _mm_store_sd(&C[(l_n*84)+29], c14_2);
    __m128d c14_3 = _mm_load_sd(&C[(l_n*84)+50]);
    __m128d a14_3 = _mm_load_sd(&A[79]);
#if defined(__SSE3__) && defined(__AVX__)
    c14_3 = _mm_add_sd(c14_3, _mm_mul_sd(a14_3, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c14_3 = _mm_add_sd(c14_3, _mm_mul_sd(a14_3, b14));
#endif
    _mm_store_sd(&C[(l_n*84)+50], c14_3);
    __m128d c14_4 = _mm_load_sd(&C[(l_n*84)+78]);
    __m128d a14_4 = _mm_load_sd(&A[80]);
#if defined(__SSE3__) && defined(__AVX__)
    c14_4 = _mm_add_sd(c14_4, _mm_mul_sd(a14_4, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c14_4 = _mm_add_sd(c14_4, _mm_mul_sd(a14_4, b14));
#endif
    _mm_store_sd(&C[(l_n*84)+78], c14_4);
#else
    C[(l_n*84)+4] += A[76] * B[(l_n*84)+14];
    C[(l_n*84)+14] += A[77] * B[(l_n*84)+14];
    C[(l_n*84)+29] += A[78] * B[(l_n*84)+14];
    C[(l_n*84)+50] += A[79] * B[(l_n*84)+14];
    C[(l_n*84)+78] += A[80] * B[(l_n*84)+14];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b15 = _mm256_broadcast_sd(&B[(l_n*84)+15]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b15 = _mm_loaddup_pd(&B[(l_n*84)+15]);
#endif
    __m128d c15_0 = _mm_load_sd(&C[(l_n*84)+5]);
    __m128d a15_0 = _mm_load_sd(&A[81]);
#if defined(__SSE3__) && defined(__AVX__)
    c15_0 = _mm_add_sd(c15_0, _mm_mul_sd(a15_0, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c15_0 = _mm_add_sd(c15_0, _mm_mul_sd(a15_0, b15));
#endif
    _mm_store_sd(&C[(l_n*84)+5], c15_0);
    __m128d c15_1 = _mm_load_sd(&C[(l_n*84)+15]);
    __m128d a15_1 = _mm_load_sd(&A[82]);
#if defined(__SSE3__) && defined(__AVX__)
    c15_1 = _mm_add_sd(c15_1, _mm_mul_sd(a15_1, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c15_1 = _mm_add_sd(c15_1, _mm_mul_sd(a15_1, b15));
#endif
    _mm_store_sd(&C[(l_n*84)+15], c15_1);
    __m128d c15_2 = _mm_load_sd(&C[(l_n*84)+30]);
    __m128d a15_2 = _mm_load_sd(&A[83]);
#if defined(__SSE3__) && defined(__AVX__)
    c15_2 = _mm_add_sd(c15_2, _mm_mul_sd(a15_2, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c15_2 = _mm_add_sd(c15_2, _mm_mul_sd(a15_2, b15));
#endif
    _mm_store_sd(&C[(l_n*84)+30], c15_2);
    __m128d c15_3 = _mm_load_sd(&C[(l_n*84)+51]);
    __m128d a15_3 = _mm_load_sd(&A[84]);
#if defined(__SSE3__) && defined(__AVX__)
    c15_3 = _mm_add_sd(c15_3, _mm_mul_sd(a15_3, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c15_3 = _mm_add_sd(c15_3, _mm_mul_sd(a15_3, b15));
#endif
    _mm_store_sd(&C[(l_n*84)+51], c15_3);
    __m128d c15_4 = _mm_load_sd(&C[(l_n*84)+79]);
    __m128d a15_4 = _mm_load_sd(&A[85]);
#if defined(__SSE3__) && defined(__AVX__)
    c15_4 = _mm_add_sd(c15_4, _mm_mul_sd(a15_4, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c15_4 = _mm_add_sd(c15_4, _mm_mul_sd(a15_4, b15));
#endif
    _mm_store_sd(&C[(l_n*84)+79], c15_4);
#else
    C[(l_n*84)+5] += A[81] * B[(l_n*84)+15];
    C[(l_n*84)+15] += A[82] * B[(l_n*84)+15];
    C[(l_n*84)+30] += A[83] * B[(l_n*84)+15];
    C[(l_n*84)+51] += A[84] * B[(l_n*84)+15];
    C[(l_n*84)+79] += A[85] * B[(l_n*84)+15];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b16 = _mm256_broadcast_sd(&B[(l_n*84)+16]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b16 = _mm_loaddup_pd(&B[(l_n*84)+16]);
#endif
    __m128d c16_0 = _mm_load_sd(&C[(l_n*84)+6]);
    __m128d a16_0 = _mm_load_sd(&A[86]);
#if defined(__SSE3__) && defined(__AVX__)
    c16_0 = _mm_add_sd(c16_0, _mm_mul_sd(a16_0, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c16_0 = _mm_add_sd(c16_0, _mm_mul_sd(a16_0, b16));
#endif
    _mm_store_sd(&C[(l_n*84)+6], c16_0);
    __m128d c16_1 = _mm_load_sd(&C[(l_n*84)+16]);
    __m128d a16_1 = _mm_load_sd(&A[87]);
#if defined(__SSE3__) && defined(__AVX__)
    c16_1 = _mm_add_sd(c16_1, _mm_mul_sd(a16_1, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c16_1 = _mm_add_sd(c16_1, _mm_mul_sd(a16_1, b16));
#endif
    _mm_store_sd(&C[(l_n*84)+16], c16_1);
    __m128d c16_2 = _mm_load_sd(&C[(l_n*84)+31]);
    __m128d a16_2 = _mm_load_sd(&A[88]);
#if defined(__SSE3__) && defined(__AVX__)
    c16_2 = _mm_add_sd(c16_2, _mm_mul_sd(a16_2, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c16_2 = _mm_add_sd(c16_2, _mm_mul_sd(a16_2, b16));
#endif
    _mm_store_sd(&C[(l_n*84)+31], c16_2);
    __m128d c16_3 = _mm_load_sd(&C[(l_n*84)+52]);
    __m128d a16_3 = _mm_load_sd(&A[89]);
#if defined(__SSE3__) && defined(__AVX__)
    c16_3 = _mm_add_sd(c16_3, _mm_mul_sd(a16_3, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c16_3 = _mm_add_sd(c16_3, _mm_mul_sd(a16_3, b16));
#endif
    _mm_store_sd(&C[(l_n*84)+52], c16_3);
    __m128d c16_4 = _mm_load_sd(&C[(l_n*84)+80]);
    __m128d a16_4 = _mm_load_sd(&A[90]);
#if defined(__SSE3__) && defined(__AVX__)
    c16_4 = _mm_add_sd(c16_4, _mm_mul_sd(a16_4, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c16_4 = _mm_add_sd(c16_4, _mm_mul_sd(a16_4, b16));
#endif
    _mm_store_sd(&C[(l_n*84)+80], c16_4);
#else
    C[(l_n*84)+6] += A[86] * B[(l_n*84)+16];
    C[(l_n*84)+16] += A[87] * B[(l_n*84)+16];
    C[(l_n*84)+31] += A[88] * B[(l_n*84)+16];
    C[(l_n*84)+52] += A[89] * B[(l_n*84)+16];
    C[(l_n*84)+80] += A[90] * B[(l_n*84)+16];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b17 = _mm256_broadcast_sd(&B[(l_n*84)+17]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b17 = _mm_loaddup_pd(&B[(l_n*84)+17]);
#endif
    __m128d c17_0 = _mm_load_sd(&C[(l_n*84)+1]);
    __m128d a17_0 = _mm_load_sd(&A[91]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_0 = _mm_add_sd(c17_0, _mm_mul_sd(a17_0, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_0 = _mm_add_sd(c17_0, _mm_mul_sd(a17_0, b17));
#endif
    _mm_store_sd(&C[(l_n*84)+1], c17_0);
    __m128d c17_1 = _mm_load_sd(&C[(l_n*84)+7]);
    __m128d a17_1 = _mm_load_sd(&A[92]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_1 = _mm_add_sd(c17_1, _mm_mul_sd(a17_1, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_1 = _mm_add_sd(c17_1, _mm_mul_sd(a17_1, b17));
#endif
    _mm_store_sd(&C[(l_n*84)+7], c17_1);
    __m128d c17_2 = _mm_load_sd(&C[(l_n*84)+17]);
    __m128d a17_2 = _mm_load_sd(&A[93]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_2 = _mm_add_sd(c17_2, _mm_mul_sd(a17_2, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_2 = _mm_add_sd(c17_2, _mm_mul_sd(a17_2, b17));
#endif
    _mm_store_sd(&C[(l_n*84)+17], c17_2);
    __m128d c17_3 = _mm_load_sd(&C[(l_n*84)+32]);
    __m128d a17_3 = _mm_load_sd(&A[94]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_3 = _mm_add_sd(c17_3, _mm_mul_sd(a17_3, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_3 = _mm_add_sd(c17_3, _mm_mul_sd(a17_3, b17));
#endif
    _mm_store_sd(&C[(l_n*84)+32], c17_3);
    __m128d c17_4 = _mm_load_sd(&C[(l_n*84)+53]);
    __m128d a17_4 = _mm_load_sd(&A[95]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_4 = _mm_add_sd(c17_4, _mm_mul_sd(a17_4, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_4 = _mm_add_sd(c17_4, _mm_mul_sd(a17_4, b17));
#endif
    _mm_store_sd(&C[(l_n*84)+53], c17_4);
    __m128d c17_5 = _mm_load_sd(&C[(l_n*84)+81]);
    __m128d a17_5 = _mm_load_sd(&A[96]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_5 = _mm_add_sd(c17_5, _mm_mul_sd(a17_5, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_5 = _mm_add_sd(c17_5, _mm_mul_sd(a17_5, b17));
#endif
    _mm_store_sd(&C[(l_n*84)+81], c17_5);
#else
    C[(l_n*84)+1] += A[91] * B[(l_n*84)+17];
    C[(l_n*84)+7] += A[92] * B[(l_n*84)+17];
    C[(l_n*84)+17] += A[93] * B[(l_n*84)+17];
    C[(l_n*84)+32] += A[94] * B[(l_n*84)+17];
    C[(l_n*84)+53] += A[95] * B[(l_n*84)+17];
    C[(l_n*84)+81] += A[96] * B[(l_n*84)+17];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b18 = _mm256_broadcast_sd(&B[(l_n*84)+18]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b18 = _mm_loaddup_pd(&B[(l_n*84)+18]);
#endif
    __m128d c18_0 = _mm_load_sd(&C[(l_n*84)+2]);
    __m128d a18_0 = _mm_load_sd(&A[97]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_0 = _mm_add_sd(c18_0, _mm_mul_sd(a18_0, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_0 = _mm_add_sd(c18_0, _mm_mul_sd(a18_0, b18));
#endif
    _mm_store_sd(&C[(l_n*84)+2], c18_0);
    __m128d c18_1 = _mm_load_sd(&C[(l_n*84)+8]);
    __m128d a18_1 = _mm_load_sd(&A[98]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_1 = _mm_add_sd(c18_1, _mm_mul_sd(a18_1, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_1 = _mm_add_sd(c18_1, _mm_mul_sd(a18_1, b18));
#endif
    _mm_store_sd(&C[(l_n*84)+8], c18_1);
    __m128d c18_2 = _mm_load_sd(&C[(l_n*84)+18]);
    __m128d a18_2 = _mm_load_sd(&A[99]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_2 = _mm_add_sd(c18_2, _mm_mul_sd(a18_2, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_2 = _mm_add_sd(c18_2, _mm_mul_sd(a18_2, b18));
#endif
    _mm_store_sd(&C[(l_n*84)+18], c18_2);
    __m128d c18_3 = _mm_load_sd(&C[(l_n*84)+33]);
    __m128d a18_3 = _mm_load_sd(&A[100]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_3 = _mm_add_sd(c18_3, _mm_mul_sd(a18_3, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_3 = _mm_add_sd(c18_3, _mm_mul_sd(a18_3, b18));
#endif
    _mm_store_sd(&C[(l_n*84)+33], c18_3);
    __m128d c18_4 = _mm_load_sd(&C[(l_n*84)+54]);
    __m128d a18_4 = _mm_load_sd(&A[101]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_4 = _mm_add_sd(c18_4, _mm_mul_sd(a18_4, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_4 = _mm_add_sd(c18_4, _mm_mul_sd(a18_4, b18));
#endif
    _mm_store_sd(&C[(l_n*84)+54], c18_4);
    __m128d c18_5 = _mm_load_sd(&C[(l_n*84)+82]);
    __m128d a18_5 = _mm_load_sd(&A[102]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_5 = _mm_add_sd(c18_5, _mm_mul_sd(a18_5, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_5 = _mm_add_sd(c18_5, _mm_mul_sd(a18_5, b18));
#endif
    _mm_store_sd(&C[(l_n*84)+82], c18_5);
#else
    C[(l_n*84)+2] += A[97] * B[(l_n*84)+18];
    C[(l_n*84)+8] += A[98] * B[(l_n*84)+18];
    C[(l_n*84)+18] += A[99] * B[(l_n*84)+18];
    C[(l_n*84)+33] += A[100] * B[(l_n*84)+18];
    C[(l_n*84)+54] += A[101] * B[(l_n*84)+18];
    C[(l_n*84)+82] += A[102] * B[(l_n*84)+18];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b19 = _mm256_broadcast_sd(&B[(l_n*84)+19]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b19 = _mm_loaddup_pd(&B[(l_n*84)+19]);
#endif
    __m128d c19_0 = _mm_load_sd(&C[(l_n*84)+0]);
    __m128d a19_0 = _mm_load_sd(&A[103]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_0 = _mm_add_sd(c19_0, _mm_mul_sd(a19_0, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_0 = _mm_add_sd(c19_0, _mm_mul_sd(a19_0, b19));
#endif
    _mm_store_sd(&C[(l_n*84)+0], c19_0);
    __m128d c19_1 = _mm_load_sd(&C[(l_n*84)+3]);
    __m128d a19_1 = _mm_load_sd(&A[104]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_1 = _mm_add_sd(c19_1, _mm_mul_sd(a19_1, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_1 = _mm_add_sd(c19_1, _mm_mul_sd(a19_1, b19));
#endif
    _mm_store_sd(&C[(l_n*84)+3], c19_1);
    __m128d c19_2 = _mm_load_sd(&C[(l_n*84)+9]);
    __m128d a19_2 = _mm_load_sd(&A[105]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_2 = _mm_add_sd(c19_2, _mm_mul_sd(a19_2, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_2 = _mm_add_sd(c19_2, _mm_mul_sd(a19_2, b19));
#endif
    _mm_store_sd(&C[(l_n*84)+9], c19_2);
    __m128d c19_3 = _mm_load_sd(&C[(l_n*84)+19]);
    __m128d a19_3 = _mm_load_sd(&A[106]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_3 = _mm_add_sd(c19_3, _mm_mul_sd(a19_3, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_3 = _mm_add_sd(c19_3, _mm_mul_sd(a19_3, b19));
#endif
    _mm_store_sd(&C[(l_n*84)+19], c19_3);
    __m128d c19_4 = _mm_load_sd(&C[(l_n*84)+34]);
    __m128d a19_4 = _mm_load_sd(&A[107]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_4 = _mm_add_sd(c19_4, _mm_mul_sd(a19_4, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_4 = _mm_add_sd(c19_4, _mm_mul_sd(a19_4, b19));
#endif
    _mm_store_sd(&C[(l_n*84)+34], c19_4);
    __m128d c19_5 = _mm_load_sd(&C[(l_n*84)+55]);
    __m128d a19_5 = _mm_load_sd(&A[108]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_5 = _mm_add_sd(c19_5, _mm_mul_sd(a19_5, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_5 = _mm_add_sd(c19_5, _mm_mul_sd(a19_5, b19));
#endif
    _mm_store_sd(&C[(l_n*84)+55], c19_5);
    __m128d c19_6 = _mm_load_sd(&C[(l_n*84)+83]);
    __m128d a19_6 = _mm_load_sd(&A[109]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_6 = _mm_add_sd(c19_6, _mm_mul_sd(a19_6, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_6 = _mm_add_sd(c19_6, _mm_mul_sd(a19_6, b19));
#endif
    _mm_store_sd(&C[(l_n*84)+83], c19_6);
#else
    C[(l_n*84)+0] += A[103] * B[(l_n*84)+19];
    C[(l_n*84)+3] += A[104] * B[(l_n*84)+19];
    C[(l_n*84)+9] += A[105] * B[(l_n*84)+19];
    C[(l_n*84)+19] += A[106] * B[(l_n*84)+19];
    C[(l_n*84)+34] += A[107] * B[(l_n*84)+19];
    C[(l_n*84)+55] += A[108] * B[(l_n*84)+19];
    C[(l_n*84)+83] += A[109] * B[(l_n*84)+19];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b20 = _mm256_broadcast_sd(&B[(l_n*84)+20]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b20 = _mm_loaddup_pd(&B[(l_n*84)+20]);
#endif
    __m128d c20_0 = _mm_load_sd(&C[(l_n*84)+20]);
    __m128d a20_0 = _mm_load_sd(&A[110]);
#if defined(__SSE3__) && defined(__AVX__)
    c20_0 = _mm_add_sd(c20_0, _mm_mul_sd(a20_0, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c20_0 = _mm_add_sd(c20_0, _mm_mul_sd(a20_0, b20));
#endif
    _mm_store_sd(&C[(l_n*84)+20], c20_0);
    __m128d c20_1 = _mm_load_sd(&C[(l_n*84)+41]);
    __m128d a20_1 = _mm_load_sd(&A[111]);
#if defined(__SSE3__) && defined(__AVX__)
    c20_1 = _mm_add_sd(c20_1, _mm_mul_sd(a20_1, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c20_1 = _mm_add_sd(c20_1, _mm_mul_sd(a20_1, b20));
#endif
    _mm_store_sd(&C[(l_n*84)+41], c20_1);
    __m128d c20_2 = _mm_load_sd(&C[(l_n*84)+69]);
    __m128d a20_2 = _mm_load_sd(&A[112]);
#if defined(__SSE3__) && defined(__AVX__)
    c20_2 = _mm_add_sd(c20_2, _mm_mul_sd(a20_2, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c20_2 = _mm_add_sd(c20_2, _mm_mul_sd(a20_2, b20));
#endif
    _mm_store_sd(&C[(l_n*84)+69], c20_2);
#else
    C[(l_n*84)+20] += A[110] * B[(l_n*84)+20];
    C[(l_n*84)+41] += A[111] * B[(l_n*84)+20];
    C[(l_n*84)+69] += A[112] * B[(l_n*84)+20];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b21 = _mm256_broadcast_sd(&B[(l_n*84)+21]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b21 = _mm_loaddup_pd(&B[(l_n*84)+21]);
#endif
    __m128d c21_0 = _mm_load_sd(&C[(l_n*84)+21]);
    __m128d a21_0 = _mm_load_sd(&A[113]);
#if defined(__SSE3__) && defined(__AVX__)
    c21_0 = _mm_add_sd(c21_0, _mm_mul_sd(a21_0, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c21_0 = _mm_add_sd(c21_0, _mm_mul_sd(a21_0, b21));
#endif
    _mm_store_sd(&C[(l_n*84)+21], c21_0);
    __m128d c21_1 = _mm_load_sd(&C[(l_n*84)+42]);
    __m128d a21_1 = _mm_load_sd(&A[114]);
#if defined(__SSE3__) && defined(__AVX__)
    c21_1 = _mm_add_sd(c21_1, _mm_mul_sd(a21_1, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c21_1 = _mm_add_sd(c21_1, _mm_mul_sd(a21_1, b21));
#endif
    _mm_store_sd(&C[(l_n*84)+42], c21_1);
    __m128d c21_2 = _mm_load_sd(&C[(l_n*84)+70]);
    __m128d a21_2 = _mm_load_sd(&A[115]);
#if defined(__SSE3__) && defined(__AVX__)
    c21_2 = _mm_add_sd(c21_2, _mm_mul_sd(a21_2, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c21_2 = _mm_add_sd(c21_2, _mm_mul_sd(a21_2, b21));
#endif
    _mm_store_sd(&C[(l_n*84)+70], c21_2);
#else
    C[(l_n*84)+21] += A[113] * B[(l_n*84)+21];
    C[(l_n*84)+42] += A[114] * B[(l_n*84)+21];
    C[(l_n*84)+70] += A[115] * B[(l_n*84)+21];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b22 = _mm256_broadcast_sd(&B[(l_n*84)+22]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b22 = _mm_loaddup_pd(&B[(l_n*84)+22]);
#endif
    __m128d c22_0 = _mm_load_sd(&C[(l_n*84)+22]);
    __m128d a22_0 = _mm_load_sd(&A[116]);
#if defined(__SSE3__) && defined(__AVX__)
    c22_0 = _mm_add_sd(c22_0, _mm_mul_sd(a22_0, _mm256_castpd256_pd128(b22)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c22_0 = _mm_add_sd(c22_0, _mm_mul_sd(a22_0, b22));
#endif
    _mm_store_sd(&C[(l_n*84)+22], c22_0);
    __m128d c22_1 = _mm_load_sd(&C[(l_n*84)+43]);
    __m128d a22_1 = _mm_load_sd(&A[117]);
#if defined(__SSE3__) && defined(__AVX__)
    c22_1 = _mm_add_sd(c22_1, _mm_mul_sd(a22_1, _mm256_castpd256_pd128(b22)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c22_1 = _mm_add_sd(c22_1, _mm_mul_sd(a22_1, b22));
#endif
    _mm_store_sd(&C[(l_n*84)+43], c22_1);
    __m128d c22_2 = _mm_load_sd(&C[(l_n*84)+71]);
    __m128d a22_2 = _mm_load_sd(&A[118]);
#if defined(__SSE3__) && defined(__AVX__)
    c22_2 = _mm_add_sd(c22_2, _mm_mul_sd(a22_2, _mm256_castpd256_pd128(b22)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c22_2 = _mm_add_sd(c22_2, _mm_mul_sd(a22_2, b22));
#endif
    _mm_store_sd(&C[(l_n*84)+71], c22_2);
#else
    C[(l_n*84)+22] += A[116] * B[(l_n*84)+22];
    C[(l_n*84)+43] += A[117] * B[(l_n*84)+22];
    C[(l_n*84)+71] += A[118] * B[(l_n*84)+22];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b23 = _mm256_broadcast_sd(&B[(l_n*84)+23]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b23 = _mm_loaddup_pd(&B[(l_n*84)+23]);
#endif
    __m128d c23_0 = _mm_load_sd(&C[(l_n*84)+23]);
    __m128d a23_0 = _mm_load_sd(&A[119]);
#if defined(__SSE3__) && defined(__AVX__)
    c23_0 = _mm_add_sd(c23_0, _mm_mul_sd(a23_0, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c23_0 = _mm_add_sd(c23_0, _mm_mul_sd(a23_0, b23));
#endif
    _mm_store_sd(&C[(l_n*84)+23], c23_0);
    __m128d c23_1 = _mm_load_sd(&C[(l_n*84)+44]);
    __m128d a23_1 = _mm_load_sd(&A[120]);
#if defined(__SSE3__) && defined(__AVX__)
    c23_1 = _mm_add_sd(c23_1, _mm_mul_sd(a23_1, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c23_1 = _mm_add_sd(c23_1, _mm_mul_sd(a23_1, b23));
#endif
    _mm_store_sd(&C[(l_n*84)+44], c23_1);
    __m128d c23_2 = _mm_load_sd(&C[(l_n*84)+72]);
    __m128d a23_2 = _mm_load_sd(&A[121]);
#if defined(__SSE3__) && defined(__AVX__)
    c23_2 = _mm_add_sd(c23_2, _mm_mul_sd(a23_2, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c23_2 = _mm_add_sd(c23_2, _mm_mul_sd(a23_2, b23));
#endif
    _mm_store_sd(&C[(l_n*84)+72], c23_2);
#else
    C[(l_n*84)+23] += A[119] * B[(l_n*84)+23];
    C[(l_n*84)+44] += A[120] * B[(l_n*84)+23];
    C[(l_n*84)+72] += A[121] * B[(l_n*84)+23];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b24 = _mm256_broadcast_sd(&B[(l_n*84)+24]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b24 = _mm_loaddup_pd(&B[(l_n*84)+24]);
#endif
    __m128d c24_0 = _mm_load_sd(&C[(l_n*84)+24]);
    __m128d a24_0 = _mm_load_sd(&A[122]);
#if defined(__SSE3__) && defined(__AVX__)
    c24_0 = _mm_add_sd(c24_0, _mm_mul_sd(a24_0, _mm256_castpd256_pd128(b24)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c24_0 = _mm_add_sd(c24_0, _mm_mul_sd(a24_0, b24));
#endif
    _mm_store_sd(&C[(l_n*84)+24], c24_0);
    __m128d c24_1 = _mm_load_sd(&C[(l_n*84)+45]);
    __m128d a24_1 = _mm_load_sd(&A[123]);
#if defined(__SSE3__) && defined(__AVX__)
    c24_1 = _mm_add_sd(c24_1, _mm_mul_sd(a24_1, _mm256_castpd256_pd128(b24)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c24_1 = _mm_add_sd(c24_1, _mm_mul_sd(a24_1, b24));
#endif
    _mm_store_sd(&C[(l_n*84)+45], c24_1);
    __m128d c24_2 = _mm_load_sd(&C[(l_n*84)+73]);
    __m128d a24_2 = _mm_load_sd(&A[124]);
#if defined(__SSE3__) && defined(__AVX__)
    c24_2 = _mm_add_sd(c24_2, _mm_mul_sd(a24_2, _mm256_castpd256_pd128(b24)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c24_2 = _mm_add_sd(c24_2, _mm_mul_sd(a24_2, b24));
#endif
    _mm_store_sd(&C[(l_n*84)+73], c24_2);
#else
    C[(l_n*84)+24] += A[122] * B[(l_n*84)+24];
    C[(l_n*84)+45] += A[123] * B[(l_n*84)+24];
    C[(l_n*84)+73] += A[124] * B[(l_n*84)+24];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b25 = _mm256_broadcast_sd(&B[(l_n*84)+25]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b25 = _mm_loaddup_pd(&B[(l_n*84)+25]);
#endif
    __m128d c25_0 = _mm_load_sd(&C[(l_n*84)+10]);
    __m128d a25_0 = _mm_load_sd(&A[125]);
#if defined(__SSE3__) && defined(__AVX__)
    c25_0 = _mm_add_sd(c25_0, _mm_mul_sd(a25_0, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c25_0 = _mm_add_sd(c25_0, _mm_mul_sd(a25_0, b25));
#endif
    _mm_store_sd(&C[(l_n*84)+10], c25_0);
    __m128d c25_1 = _mm_load_sd(&C[(l_n*84)+25]);
    __m128d a25_1 = _mm_load_sd(&A[126]);
#if defined(__SSE3__) && defined(__AVX__)
    c25_1 = _mm_add_sd(c25_1, _mm_mul_sd(a25_1, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c25_1 = _mm_add_sd(c25_1, _mm_mul_sd(a25_1, b25));
#endif
    _mm_store_sd(&C[(l_n*84)+25], c25_1);
    __m128d c25_2 = _mm_load_sd(&C[(l_n*84)+46]);
    __m128d a25_2 = _mm_load_sd(&A[127]);
#if defined(__SSE3__) && defined(__AVX__)
    c25_2 = _mm_add_sd(c25_2, _mm_mul_sd(a25_2, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c25_2 = _mm_add_sd(c25_2, _mm_mul_sd(a25_2, b25));
#endif
    _mm_store_sd(&C[(l_n*84)+46], c25_2);
    __m128d c25_3 = _mm_load_sd(&C[(l_n*84)+74]);
    __m128d a25_3 = _mm_load_sd(&A[128]);
#if defined(__SSE3__) && defined(__AVX__)
    c25_3 = _mm_add_sd(c25_3, _mm_mul_sd(a25_3, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c25_3 = _mm_add_sd(c25_3, _mm_mul_sd(a25_3, b25));
#endif
    _mm_store_sd(&C[(l_n*84)+74], c25_3);
#else
    C[(l_n*84)+10] += A[125] * B[(l_n*84)+25];
    C[(l_n*84)+25] += A[126] * B[(l_n*84)+25];
    C[(l_n*84)+46] += A[127] * B[(l_n*84)+25];
    C[(l_n*84)+74] += A[128] * B[(l_n*84)+25];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b26 = _mm256_broadcast_sd(&B[(l_n*84)+26]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b26 = _mm_loaddup_pd(&B[(l_n*84)+26]);
#endif
    __m128d c26_0 = _mm_load_sd(&C[(l_n*84)+11]);
    __m128d a26_0 = _mm_load_sd(&A[129]);
#if defined(__SSE3__) && defined(__AVX__)
    c26_0 = _mm_add_sd(c26_0, _mm_mul_sd(a26_0, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c26_0 = _mm_add_sd(c26_0, _mm_mul_sd(a26_0, b26));
#endif
    _mm_store_sd(&C[(l_n*84)+11], c26_0);
    __m128d c26_1 = _mm_load_sd(&C[(l_n*84)+26]);
    __m128d a26_1 = _mm_load_sd(&A[130]);
#if defined(__SSE3__) && defined(__AVX__)
    c26_1 = _mm_add_sd(c26_1, _mm_mul_sd(a26_1, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c26_1 = _mm_add_sd(c26_1, _mm_mul_sd(a26_1, b26));
#endif
    _mm_store_sd(&C[(l_n*84)+26], c26_1);
    __m128d c26_2 = _mm_load_sd(&C[(l_n*84)+47]);
    __m128d a26_2 = _mm_load_sd(&A[131]);
#if defined(__SSE3__) && defined(__AVX__)
    c26_2 = _mm_add_sd(c26_2, _mm_mul_sd(a26_2, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c26_2 = _mm_add_sd(c26_2, _mm_mul_sd(a26_2, b26));
#endif
    _mm_store_sd(&C[(l_n*84)+47], c26_2);
    __m128d c26_3 = _mm_load_sd(&C[(l_n*84)+75]);
    __m128d a26_3 = _mm_load_sd(&A[132]);
#if defined(__SSE3__) && defined(__AVX__)
    c26_3 = _mm_add_sd(c26_3, _mm_mul_sd(a26_3, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c26_3 = _mm_add_sd(c26_3, _mm_mul_sd(a26_3, b26));
#endif
    _mm_store_sd(&C[(l_n*84)+75], c26_3);
#else
    C[(l_n*84)+11] += A[129] * B[(l_n*84)+26];
    C[(l_n*84)+26] += A[130] * B[(l_n*84)+26];
    C[(l_n*84)+47] += A[131] * B[(l_n*84)+26];
    C[(l_n*84)+75] += A[132] * B[(l_n*84)+26];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b27 = _mm256_broadcast_sd(&B[(l_n*84)+27]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b27 = _mm_loaddup_pd(&B[(l_n*84)+27]);
#endif
    __m128d c27_0 = _mm_load_sd(&C[(l_n*84)+12]);
    __m128d a27_0 = _mm_load_sd(&A[133]);
#if defined(__SSE3__) && defined(__AVX__)
    c27_0 = _mm_add_sd(c27_0, _mm_mul_sd(a27_0, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c27_0 = _mm_add_sd(c27_0, _mm_mul_sd(a27_0, b27));
#endif
    _mm_store_sd(&C[(l_n*84)+12], c27_0);
    __m128d c27_1 = _mm_load_sd(&C[(l_n*84)+27]);
    __m128d a27_1 = _mm_load_sd(&A[134]);
#if defined(__SSE3__) && defined(__AVX__)
    c27_1 = _mm_add_sd(c27_1, _mm_mul_sd(a27_1, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c27_1 = _mm_add_sd(c27_1, _mm_mul_sd(a27_1, b27));
#endif
    _mm_store_sd(&C[(l_n*84)+27], c27_1);
    __m128d c27_2 = _mm_load_sd(&C[(l_n*84)+48]);
    __m128d a27_2 = _mm_load_sd(&A[135]);
#if defined(__SSE3__) && defined(__AVX__)
    c27_2 = _mm_add_sd(c27_2, _mm_mul_sd(a27_2, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c27_2 = _mm_add_sd(c27_2, _mm_mul_sd(a27_2, b27));
#endif
    _mm_store_sd(&C[(l_n*84)+48], c27_2);
    __m128d c27_3 = _mm_load_sd(&C[(l_n*84)+76]);
    __m128d a27_3 = _mm_load_sd(&A[136]);
#if defined(__SSE3__) && defined(__AVX__)
    c27_3 = _mm_add_sd(c27_3, _mm_mul_sd(a27_3, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c27_3 = _mm_add_sd(c27_3, _mm_mul_sd(a27_3, b27));
#endif
    _mm_store_sd(&C[(l_n*84)+76], c27_3);
#else
    C[(l_n*84)+12] += A[133] * B[(l_n*84)+27];
    C[(l_n*84)+27] += A[134] * B[(l_n*84)+27];
    C[(l_n*84)+48] += A[135] * B[(l_n*84)+27];
    C[(l_n*84)+76] += A[136] * B[(l_n*84)+27];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b28 = _mm256_broadcast_sd(&B[(l_n*84)+28]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b28 = _mm_loaddup_pd(&B[(l_n*84)+28]);
#endif
    __m128d c28_0 = _mm_load_sd(&C[(l_n*84)+13]);
    __m128d a28_0 = _mm_load_sd(&A[137]);
#if defined(__SSE3__) && defined(__AVX__)
    c28_0 = _mm_add_sd(c28_0, _mm_mul_sd(a28_0, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c28_0 = _mm_add_sd(c28_0, _mm_mul_sd(a28_0, b28));
#endif
    _mm_store_sd(&C[(l_n*84)+13], c28_0);
    __m128d c28_1 = _mm_load_sd(&C[(l_n*84)+28]);
    __m128d a28_1 = _mm_load_sd(&A[138]);
#if defined(__SSE3__) && defined(__AVX__)
    c28_1 = _mm_add_sd(c28_1, _mm_mul_sd(a28_1, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c28_1 = _mm_add_sd(c28_1, _mm_mul_sd(a28_1, b28));
#endif
    _mm_store_sd(&C[(l_n*84)+28], c28_1);
    __m128d c28_2 = _mm_load_sd(&C[(l_n*84)+49]);
    __m128d a28_2 = _mm_load_sd(&A[139]);
#if defined(__SSE3__) && defined(__AVX__)
    c28_2 = _mm_add_sd(c28_2, _mm_mul_sd(a28_2, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c28_2 = _mm_add_sd(c28_2, _mm_mul_sd(a28_2, b28));
#endif
    _mm_store_sd(&C[(l_n*84)+49], c28_2);
    __m128d c28_3 = _mm_load_sd(&C[(l_n*84)+77]);
    __m128d a28_3 = _mm_load_sd(&A[140]);
#if defined(__SSE3__) && defined(__AVX__)
    c28_3 = _mm_add_sd(c28_3, _mm_mul_sd(a28_3, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c28_3 = _mm_add_sd(c28_3, _mm_mul_sd(a28_3, b28));
#endif
    _mm_store_sd(&C[(l_n*84)+77], c28_3);
#else
    C[(l_n*84)+13] += A[137] * B[(l_n*84)+28];
    C[(l_n*84)+28] += A[138] * B[(l_n*84)+28];
    C[(l_n*84)+49] += A[139] * B[(l_n*84)+28];
    C[(l_n*84)+77] += A[140] * B[(l_n*84)+28];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b29 = _mm256_broadcast_sd(&B[(l_n*84)+29]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b29 = _mm_loaddup_pd(&B[(l_n*84)+29]);
#endif
    __m128d c29_0 = _mm_load_sd(&C[(l_n*84)+4]);
    __m128d a29_0 = _mm_load_sd(&A[141]);
#if defined(__SSE3__) && defined(__AVX__)
    c29_0 = _mm_add_sd(c29_0, _mm_mul_sd(a29_0, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c29_0 = _mm_add_sd(c29_0, _mm_mul_sd(a29_0, b29));
#endif
    _mm_store_sd(&C[(l_n*84)+4], c29_0);
    __m128d c29_1 = _mm_load_sd(&C[(l_n*84)+14]);
    __m128d a29_1 = _mm_load_sd(&A[142]);
#if defined(__SSE3__) && defined(__AVX__)
    c29_1 = _mm_add_sd(c29_1, _mm_mul_sd(a29_1, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c29_1 = _mm_add_sd(c29_1, _mm_mul_sd(a29_1, b29));
#endif
    _mm_store_sd(&C[(l_n*84)+14], c29_1);
    __m128d c29_2 = _mm_load_sd(&C[(l_n*84)+29]);
    __m128d a29_2 = _mm_load_sd(&A[143]);
#if defined(__SSE3__) && defined(__AVX__)
    c29_2 = _mm_add_sd(c29_2, _mm_mul_sd(a29_2, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c29_2 = _mm_add_sd(c29_2, _mm_mul_sd(a29_2, b29));
#endif
    _mm_store_sd(&C[(l_n*84)+29], c29_2);
    __m128d c29_3 = _mm_load_sd(&C[(l_n*84)+50]);
    __m128d a29_3 = _mm_load_sd(&A[144]);
#if defined(__SSE3__) && defined(__AVX__)
    c29_3 = _mm_add_sd(c29_3, _mm_mul_sd(a29_3, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c29_3 = _mm_add_sd(c29_3, _mm_mul_sd(a29_3, b29));
#endif
    _mm_store_sd(&C[(l_n*84)+50], c29_3);
    __m128d c29_4 = _mm_load_sd(&C[(l_n*84)+78]);
    __m128d a29_4 = _mm_load_sd(&A[145]);
#if defined(__SSE3__) && defined(__AVX__)
    c29_4 = _mm_add_sd(c29_4, _mm_mul_sd(a29_4, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c29_4 = _mm_add_sd(c29_4, _mm_mul_sd(a29_4, b29));
#endif
    _mm_store_sd(&C[(l_n*84)+78], c29_4);
#else
    C[(l_n*84)+4] += A[141] * B[(l_n*84)+29];
    C[(l_n*84)+14] += A[142] * B[(l_n*84)+29];
    C[(l_n*84)+29] += A[143] * B[(l_n*84)+29];
    C[(l_n*84)+50] += A[144] * B[(l_n*84)+29];
    C[(l_n*84)+78] += A[145] * B[(l_n*84)+29];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b30 = _mm256_broadcast_sd(&B[(l_n*84)+30]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b30 = _mm_loaddup_pd(&B[(l_n*84)+30]);
#endif
    __m128d c30_0 = _mm_load_sd(&C[(l_n*84)+5]);
    __m128d a30_0 = _mm_load_sd(&A[146]);
#if defined(__SSE3__) && defined(__AVX__)
    c30_0 = _mm_add_sd(c30_0, _mm_mul_sd(a30_0, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c30_0 = _mm_add_sd(c30_0, _mm_mul_sd(a30_0, b30));
#endif
    _mm_store_sd(&C[(l_n*84)+5], c30_0);
    __m128d c30_1 = _mm_load_sd(&C[(l_n*84)+15]);
    __m128d a30_1 = _mm_load_sd(&A[147]);
#if defined(__SSE3__) && defined(__AVX__)
    c30_1 = _mm_add_sd(c30_1, _mm_mul_sd(a30_1, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c30_1 = _mm_add_sd(c30_1, _mm_mul_sd(a30_1, b30));
#endif
    _mm_store_sd(&C[(l_n*84)+15], c30_1);
    __m128d c30_2 = _mm_load_sd(&C[(l_n*84)+30]);
    __m128d a30_2 = _mm_load_sd(&A[148]);
#if defined(__SSE3__) && defined(__AVX__)
    c30_2 = _mm_add_sd(c30_2, _mm_mul_sd(a30_2, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c30_2 = _mm_add_sd(c30_2, _mm_mul_sd(a30_2, b30));
#endif
    _mm_store_sd(&C[(l_n*84)+30], c30_2);
    __m128d c30_3 = _mm_load_sd(&C[(l_n*84)+51]);
    __m128d a30_3 = _mm_load_sd(&A[149]);
#if defined(__SSE3__) && defined(__AVX__)
    c30_3 = _mm_add_sd(c30_3, _mm_mul_sd(a30_3, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c30_3 = _mm_add_sd(c30_3, _mm_mul_sd(a30_3, b30));
#endif
    _mm_store_sd(&C[(l_n*84)+51], c30_3);
    __m128d c30_4 = _mm_load_sd(&C[(l_n*84)+79]);
    __m128d a30_4 = _mm_load_sd(&A[150]);
#if defined(__SSE3__) && defined(__AVX__)
    c30_4 = _mm_add_sd(c30_4, _mm_mul_sd(a30_4, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c30_4 = _mm_add_sd(c30_4, _mm_mul_sd(a30_4, b30));
#endif
    _mm_store_sd(&C[(l_n*84)+79], c30_4);
#else
    C[(l_n*84)+5] += A[146] * B[(l_n*84)+30];
    C[(l_n*84)+15] += A[147] * B[(l_n*84)+30];
    C[(l_n*84)+30] += A[148] * B[(l_n*84)+30];
    C[(l_n*84)+51] += A[149] * B[(l_n*84)+30];
    C[(l_n*84)+79] += A[150] * B[(l_n*84)+30];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b31 = _mm256_broadcast_sd(&B[(l_n*84)+31]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b31 = _mm_loaddup_pd(&B[(l_n*84)+31]);
#endif
    __m128d c31_0 = _mm_load_sd(&C[(l_n*84)+6]);
    __m128d a31_0 = _mm_load_sd(&A[151]);
#if defined(__SSE3__) && defined(__AVX__)
    c31_0 = _mm_add_sd(c31_0, _mm_mul_sd(a31_0, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c31_0 = _mm_add_sd(c31_0, _mm_mul_sd(a31_0, b31));
#endif
    _mm_store_sd(&C[(l_n*84)+6], c31_0);
    __m128d c31_1 = _mm_load_sd(&C[(l_n*84)+16]);
    __m128d a31_1 = _mm_load_sd(&A[152]);
#if defined(__SSE3__) && defined(__AVX__)
    c31_1 = _mm_add_sd(c31_1, _mm_mul_sd(a31_1, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c31_1 = _mm_add_sd(c31_1, _mm_mul_sd(a31_1, b31));
#endif
    _mm_store_sd(&C[(l_n*84)+16], c31_1);
    __m128d c31_2 = _mm_load_sd(&C[(l_n*84)+31]);
    __m128d a31_2 = _mm_load_sd(&A[153]);
#if defined(__SSE3__) && defined(__AVX__)
    c31_2 = _mm_add_sd(c31_2, _mm_mul_sd(a31_2, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c31_2 = _mm_add_sd(c31_2, _mm_mul_sd(a31_2, b31));
#endif
    _mm_store_sd(&C[(l_n*84)+31], c31_2);
    __m128d c31_3 = _mm_load_sd(&C[(l_n*84)+52]);
    __m128d a31_3 = _mm_load_sd(&A[154]);
#if defined(__SSE3__) && defined(__AVX__)
    c31_3 = _mm_add_sd(c31_3, _mm_mul_sd(a31_3, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c31_3 = _mm_add_sd(c31_3, _mm_mul_sd(a31_3, b31));
#endif
    _mm_store_sd(&C[(l_n*84)+52], c31_3);
    __m128d c31_4 = _mm_load_sd(&C[(l_n*84)+80]);
    __m128d a31_4 = _mm_load_sd(&A[155]);
#if defined(__SSE3__) && defined(__AVX__)
    c31_4 = _mm_add_sd(c31_4, _mm_mul_sd(a31_4, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c31_4 = _mm_add_sd(c31_4, _mm_mul_sd(a31_4, b31));
#endif
    _mm_store_sd(&C[(l_n*84)+80], c31_4);
#else
    C[(l_n*84)+6] += A[151] * B[(l_n*84)+31];
    C[(l_n*84)+16] += A[152] * B[(l_n*84)+31];
    C[(l_n*84)+31] += A[153] * B[(l_n*84)+31];
    C[(l_n*84)+52] += A[154] * B[(l_n*84)+31];
    C[(l_n*84)+80] += A[155] * B[(l_n*84)+31];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b32 = _mm256_broadcast_sd(&B[(l_n*84)+32]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b32 = _mm_loaddup_pd(&B[(l_n*84)+32]);
#endif
    __m128d c32_0 = _mm_load_sd(&C[(l_n*84)+1]);
    __m128d a32_0 = _mm_load_sd(&A[156]);
#if defined(__SSE3__) && defined(__AVX__)
    c32_0 = _mm_add_sd(c32_0, _mm_mul_sd(a32_0, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c32_0 = _mm_add_sd(c32_0, _mm_mul_sd(a32_0, b32));
#endif
    _mm_store_sd(&C[(l_n*84)+1], c32_0);
    __m128d c32_1 = _mm_load_sd(&C[(l_n*84)+7]);
    __m128d a32_1 = _mm_load_sd(&A[157]);
#if defined(__SSE3__) && defined(__AVX__)
    c32_1 = _mm_add_sd(c32_1, _mm_mul_sd(a32_1, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c32_1 = _mm_add_sd(c32_1, _mm_mul_sd(a32_1, b32));
#endif
    _mm_store_sd(&C[(l_n*84)+7], c32_1);
    __m128d c32_2 = _mm_load_sd(&C[(l_n*84)+17]);
    __m128d a32_2 = _mm_load_sd(&A[158]);
#if defined(__SSE3__) && defined(__AVX__)
    c32_2 = _mm_add_sd(c32_2, _mm_mul_sd(a32_2, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c32_2 = _mm_add_sd(c32_2, _mm_mul_sd(a32_2, b32));
#endif
    _mm_store_sd(&C[(l_n*84)+17], c32_2);
    __m128d c32_3 = _mm_load_sd(&C[(l_n*84)+32]);
    __m128d a32_3 = _mm_load_sd(&A[159]);
#if defined(__SSE3__) && defined(__AVX__)
    c32_3 = _mm_add_sd(c32_3, _mm_mul_sd(a32_3, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c32_3 = _mm_add_sd(c32_3, _mm_mul_sd(a32_3, b32));
#endif
    _mm_store_sd(&C[(l_n*84)+32], c32_3);
    __m128d c32_4 = _mm_load_sd(&C[(l_n*84)+53]);
    __m128d a32_4 = _mm_load_sd(&A[160]);
#if defined(__SSE3__) && defined(__AVX__)
    c32_4 = _mm_add_sd(c32_4, _mm_mul_sd(a32_4, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c32_4 = _mm_add_sd(c32_4, _mm_mul_sd(a32_4, b32));
#endif
    _mm_store_sd(&C[(l_n*84)+53], c32_4);
    __m128d c32_5 = _mm_load_sd(&C[(l_n*84)+81]);
    __m128d a32_5 = _mm_load_sd(&A[161]);
#if defined(__SSE3__) && defined(__AVX__)
    c32_5 = _mm_add_sd(c32_5, _mm_mul_sd(a32_5, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c32_5 = _mm_add_sd(c32_5, _mm_mul_sd(a32_5, b32));
#endif
    _mm_store_sd(&C[(l_n*84)+81], c32_5);
#else
    C[(l_n*84)+1] += A[156] * B[(l_n*84)+32];
    C[(l_n*84)+7] += A[157] * B[(l_n*84)+32];
    C[(l_n*84)+17] += A[158] * B[(l_n*84)+32];
    C[(l_n*84)+32] += A[159] * B[(l_n*84)+32];
    C[(l_n*84)+53] += A[160] * B[(l_n*84)+32];
    C[(l_n*84)+81] += A[161] * B[(l_n*84)+32];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b33 = _mm256_broadcast_sd(&B[(l_n*84)+33]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b33 = _mm_loaddup_pd(&B[(l_n*84)+33]);
#endif
    __m128d c33_0 = _mm_load_sd(&C[(l_n*84)+2]);
    __m128d a33_0 = _mm_load_sd(&A[162]);
#if defined(__SSE3__) && defined(__AVX__)
    c33_0 = _mm_add_sd(c33_0, _mm_mul_sd(a33_0, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c33_0 = _mm_add_sd(c33_0, _mm_mul_sd(a33_0, b33));
#endif
    _mm_store_sd(&C[(l_n*84)+2], c33_0);
    __m128d c33_1 = _mm_load_sd(&C[(l_n*84)+8]);
    __m128d a33_1 = _mm_load_sd(&A[163]);
#if defined(__SSE3__) && defined(__AVX__)
    c33_1 = _mm_add_sd(c33_1, _mm_mul_sd(a33_1, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c33_1 = _mm_add_sd(c33_1, _mm_mul_sd(a33_1, b33));
#endif
    _mm_store_sd(&C[(l_n*84)+8], c33_1);
    __m128d c33_2 = _mm_load_sd(&C[(l_n*84)+18]);
    __m128d a33_2 = _mm_load_sd(&A[164]);
#if defined(__SSE3__) && defined(__AVX__)
    c33_2 = _mm_add_sd(c33_2, _mm_mul_sd(a33_2, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c33_2 = _mm_add_sd(c33_2, _mm_mul_sd(a33_2, b33));
#endif
    _mm_store_sd(&C[(l_n*84)+18], c33_2);
    __m128d c33_3 = _mm_load_sd(&C[(l_n*84)+33]);
    __m128d a33_3 = _mm_load_sd(&A[165]);
#if defined(__SSE3__) && defined(__AVX__)
    c33_3 = _mm_add_sd(c33_3, _mm_mul_sd(a33_3, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c33_3 = _mm_add_sd(c33_3, _mm_mul_sd(a33_3, b33));
#endif
    _mm_store_sd(&C[(l_n*84)+33], c33_3);
    __m128d c33_4 = _mm_load_sd(&C[(l_n*84)+54]);
    __m128d a33_4 = _mm_load_sd(&A[166]);
#if defined(__SSE3__) && defined(__AVX__)
    c33_4 = _mm_add_sd(c33_4, _mm_mul_sd(a33_4, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c33_4 = _mm_add_sd(c33_4, _mm_mul_sd(a33_4, b33));
#endif
    _mm_store_sd(&C[(l_n*84)+54], c33_4);
    __m128d c33_5 = _mm_load_sd(&C[(l_n*84)+82]);
    __m128d a33_5 = _mm_load_sd(&A[167]);
#if defined(__SSE3__) && defined(__AVX__)
    c33_5 = _mm_add_sd(c33_5, _mm_mul_sd(a33_5, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c33_5 = _mm_add_sd(c33_5, _mm_mul_sd(a33_5, b33));
#endif
    _mm_store_sd(&C[(l_n*84)+82], c33_5);
#else
    C[(l_n*84)+2] += A[162] * B[(l_n*84)+33];
    C[(l_n*84)+8] += A[163] * B[(l_n*84)+33];
    C[(l_n*84)+18] += A[164] * B[(l_n*84)+33];
    C[(l_n*84)+33] += A[165] * B[(l_n*84)+33];
    C[(l_n*84)+54] += A[166] * B[(l_n*84)+33];
    C[(l_n*84)+82] += A[167] * B[(l_n*84)+33];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b34 = _mm256_broadcast_sd(&B[(l_n*84)+34]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b34 = _mm_loaddup_pd(&B[(l_n*84)+34]);
#endif
    __m128d c34_0 = _mm_load_sd(&C[(l_n*84)+0]);
    __m128d a34_0 = _mm_load_sd(&A[168]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_0 = _mm_add_sd(c34_0, _mm_mul_sd(a34_0, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_0 = _mm_add_sd(c34_0, _mm_mul_sd(a34_0, b34));
#endif
    _mm_store_sd(&C[(l_n*84)+0], c34_0);
    __m128d c34_1 = _mm_load_sd(&C[(l_n*84)+3]);
    __m128d a34_1 = _mm_load_sd(&A[169]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_1 = _mm_add_sd(c34_1, _mm_mul_sd(a34_1, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_1 = _mm_add_sd(c34_1, _mm_mul_sd(a34_1, b34));
#endif
    _mm_store_sd(&C[(l_n*84)+3], c34_1);
    __m128d c34_2 = _mm_load_sd(&C[(l_n*84)+9]);
    __m128d a34_2 = _mm_load_sd(&A[170]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_2 = _mm_add_sd(c34_2, _mm_mul_sd(a34_2, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_2 = _mm_add_sd(c34_2, _mm_mul_sd(a34_2, b34));
#endif
    _mm_store_sd(&C[(l_n*84)+9], c34_2);
    __m128d c34_3 = _mm_load_sd(&C[(l_n*84)+19]);
    __m128d a34_3 = _mm_load_sd(&A[171]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_3 = _mm_add_sd(c34_3, _mm_mul_sd(a34_3, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_3 = _mm_add_sd(c34_3, _mm_mul_sd(a34_3, b34));
#endif
    _mm_store_sd(&C[(l_n*84)+19], c34_3);
    __m128d c34_4 = _mm_load_sd(&C[(l_n*84)+34]);
    __m128d a34_4 = _mm_load_sd(&A[172]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_4 = _mm_add_sd(c34_4, _mm_mul_sd(a34_4, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_4 = _mm_add_sd(c34_4, _mm_mul_sd(a34_4, b34));
#endif
    _mm_store_sd(&C[(l_n*84)+34], c34_4);
    __m128d c34_5 = _mm_load_sd(&C[(l_n*84)+55]);
    __m128d a34_5 = _mm_load_sd(&A[173]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_5 = _mm_add_sd(c34_5, _mm_mul_sd(a34_5, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_5 = _mm_add_sd(c34_5, _mm_mul_sd(a34_5, b34));
#endif
    _mm_store_sd(&C[(l_n*84)+55], c34_5);
    __m128d c34_6 = _mm_load_sd(&C[(l_n*84)+83]);
    __m128d a34_6 = _mm_load_sd(&A[174]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_6 = _mm_add_sd(c34_6, _mm_mul_sd(a34_6, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_6 = _mm_add_sd(c34_6, _mm_mul_sd(a34_6, b34));
#endif
    _mm_store_sd(&C[(l_n*84)+83], c34_6);
#else
    C[(l_n*84)+0] += A[168] * B[(l_n*84)+34];
    C[(l_n*84)+3] += A[169] * B[(l_n*84)+34];
    C[(l_n*84)+9] += A[170] * B[(l_n*84)+34];
    C[(l_n*84)+19] += A[171] * B[(l_n*84)+34];
    C[(l_n*84)+34] += A[172] * B[(l_n*84)+34];
    C[(l_n*84)+55] += A[173] * B[(l_n*84)+34];
    C[(l_n*84)+83] += A[174] * B[(l_n*84)+34];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b35 = _mm256_broadcast_sd(&B[(l_n*84)+35]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b35 = _mm_loaddup_pd(&B[(l_n*84)+35]);
#endif
    __m128d c35_0 = _mm_load_sd(&C[(l_n*84)+35]);
    __m128d a35_0 = _mm_load_sd(&A[175]);
#if defined(__SSE3__) && defined(__AVX__)
    c35_0 = _mm_add_sd(c35_0, _mm_mul_sd(a35_0, _mm256_castpd256_pd128(b35)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c35_0 = _mm_add_sd(c35_0, _mm_mul_sd(a35_0, b35));
#endif
    _mm_store_sd(&C[(l_n*84)+35], c35_0);
    __m128d c35_1 = _mm_load_sd(&C[(l_n*84)+63]);
    __m128d a35_1 = _mm_load_sd(&A[176]);
#if defined(__SSE3__) && defined(__AVX__)
    c35_1 = _mm_add_sd(c35_1, _mm_mul_sd(a35_1, _mm256_castpd256_pd128(b35)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c35_1 = _mm_add_sd(c35_1, _mm_mul_sd(a35_1, b35));
#endif
    _mm_store_sd(&C[(l_n*84)+63], c35_1);
#else
    C[(l_n*84)+35] += A[175] * B[(l_n*84)+35];
    C[(l_n*84)+63] += A[176] * B[(l_n*84)+35];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b36 = _mm256_broadcast_sd(&B[(l_n*84)+36]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b36 = _mm_loaddup_pd(&B[(l_n*84)+36]);
#endif
    __m128d c36_0 = _mm_load_sd(&C[(l_n*84)+36]);
    __m128d a36_0 = _mm_load_sd(&A[177]);
#if defined(__SSE3__) && defined(__AVX__)
    c36_0 = _mm_add_sd(c36_0, _mm_mul_sd(a36_0, _mm256_castpd256_pd128(b36)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c36_0 = _mm_add_sd(c36_0, _mm_mul_sd(a36_0, b36));
#endif
    _mm_store_sd(&C[(l_n*84)+36], c36_0);
    __m128d c36_1 = _mm_load_sd(&C[(l_n*84)+64]);
    __m128d a36_1 = _mm_load_sd(&A[178]);
#if defined(__SSE3__) && defined(__AVX__)
    c36_1 = _mm_add_sd(c36_1, _mm_mul_sd(a36_1, _mm256_castpd256_pd128(b36)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c36_1 = _mm_add_sd(c36_1, _mm_mul_sd(a36_1, b36));
#endif
    _mm_store_sd(&C[(l_n*84)+64], c36_1);
#else
    C[(l_n*84)+36] += A[177] * B[(l_n*84)+36];
    C[(l_n*84)+64] += A[178] * B[(l_n*84)+36];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b37 = _mm256_broadcast_sd(&B[(l_n*84)+37]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b37 = _mm_loaddup_pd(&B[(l_n*84)+37]);
#endif
    __m128d c37_0 = _mm_load_sd(&C[(l_n*84)+37]);
    __m128d a37_0 = _mm_load_sd(&A[179]);
#if defined(__SSE3__) && defined(__AVX__)
    c37_0 = _mm_add_sd(c37_0, _mm_mul_sd(a37_0, _mm256_castpd256_pd128(b37)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c37_0 = _mm_add_sd(c37_0, _mm_mul_sd(a37_0, b37));
#endif
    _mm_store_sd(&C[(l_n*84)+37], c37_0);
    __m128d c37_1 = _mm_load_sd(&C[(l_n*84)+65]);
    __m128d a37_1 = _mm_load_sd(&A[180]);
#if defined(__SSE3__) && defined(__AVX__)
    c37_1 = _mm_add_sd(c37_1, _mm_mul_sd(a37_1, _mm256_castpd256_pd128(b37)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c37_1 = _mm_add_sd(c37_1, _mm_mul_sd(a37_1, b37));
#endif
    _mm_store_sd(&C[(l_n*84)+65], c37_1);
#else
    C[(l_n*84)+37] += A[179] * B[(l_n*84)+37];
    C[(l_n*84)+65] += A[180] * B[(l_n*84)+37];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b38 = _mm256_broadcast_sd(&B[(l_n*84)+38]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b38 = _mm_loaddup_pd(&B[(l_n*84)+38]);
#endif
    __m128d c38_0 = _mm_load_sd(&C[(l_n*84)+38]);
    __m128d a38_0 = _mm_load_sd(&A[181]);
#if defined(__SSE3__) && defined(__AVX__)
    c38_0 = _mm_add_sd(c38_0, _mm_mul_sd(a38_0, _mm256_castpd256_pd128(b38)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c38_0 = _mm_add_sd(c38_0, _mm_mul_sd(a38_0, b38));
#endif
    _mm_store_sd(&C[(l_n*84)+38], c38_0);
    __m128d c38_1 = _mm_load_sd(&C[(l_n*84)+66]);
    __m128d a38_1 = _mm_load_sd(&A[182]);
#if defined(__SSE3__) && defined(__AVX__)
    c38_1 = _mm_add_sd(c38_1, _mm_mul_sd(a38_1, _mm256_castpd256_pd128(b38)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c38_1 = _mm_add_sd(c38_1, _mm_mul_sd(a38_1, b38));
#endif
    _mm_store_sd(&C[(l_n*84)+66], c38_1);
#else
    C[(l_n*84)+38] += A[181] * B[(l_n*84)+38];
    C[(l_n*84)+66] += A[182] * B[(l_n*84)+38];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b39 = _mm256_broadcast_sd(&B[(l_n*84)+39]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b39 = _mm_loaddup_pd(&B[(l_n*84)+39]);
#endif
    __m128d c39_0 = _mm_load_sd(&C[(l_n*84)+39]);
    __m128d a39_0 = _mm_load_sd(&A[183]);
#if defined(__SSE3__) && defined(__AVX__)
    c39_0 = _mm_add_sd(c39_0, _mm_mul_sd(a39_0, _mm256_castpd256_pd128(b39)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c39_0 = _mm_add_sd(c39_0, _mm_mul_sd(a39_0, b39));
#endif
    _mm_store_sd(&C[(l_n*84)+39], c39_0);
    __m128d c39_1 = _mm_load_sd(&C[(l_n*84)+67]);
    __m128d a39_1 = _mm_load_sd(&A[184]);
#if defined(__SSE3__) && defined(__AVX__)
    c39_1 = _mm_add_sd(c39_1, _mm_mul_sd(a39_1, _mm256_castpd256_pd128(b39)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c39_1 = _mm_add_sd(c39_1, _mm_mul_sd(a39_1, b39));
#endif
    _mm_store_sd(&C[(l_n*84)+67], c39_1);
#else
    C[(l_n*84)+39] += A[183] * B[(l_n*84)+39];
    C[(l_n*84)+67] += A[184] * B[(l_n*84)+39];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b40 = _mm256_broadcast_sd(&B[(l_n*84)+40]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b40 = _mm_loaddup_pd(&B[(l_n*84)+40]);
#endif
    __m128d c40_0 = _mm_load_sd(&C[(l_n*84)+40]);
    __m128d a40_0 = _mm_load_sd(&A[185]);
#if defined(__SSE3__) && defined(__AVX__)
    c40_0 = _mm_add_sd(c40_0, _mm_mul_sd(a40_0, _mm256_castpd256_pd128(b40)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c40_0 = _mm_add_sd(c40_0, _mm_mul_sd(a40_0, b40));
#endif
    _mm_store_sd(&C[(l_n*84)+40], c40_0);
    __m128d c40_1 = _mm_load_sd(&C[(l_n*84)+68]);
    __m128d a40_1 = _mm_load_sd(&A[186]);
#if defined(__SSE3__) && defined(__AVX__)
    c40_1 = _mm_add_sd(c40_1, _mm_mul_sd(a40_1, _mm256_castpd256_pd128(b40)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c40_1 = _mm_add_sd(c40_1, _mm_mul_sd(a40_1, b40));
#endif
    _mm_store_sd(&C[(l_n*84)+68], c40_1);
#else
    C[(l_n*84)+40] += A[185] * B[(l_n*84)+40];
    C[(l_n*84)+68] += A[186] * B[(l_n*84)+40];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b41 = _mm256_broadcast_sd(&B[(l_n*84)+41]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b41 = _mm_loaddup_pd(&B[(l_n*84)+41]);
#endif
    __m128d c41_0 = _mm_load_sd(&C[(l_n*84)+20]);
    __m128d a41_0 = _mm_load_sd(&A[187]);
#if defined(__SSE3__) && defined(__AVX__)
    c41_0 = _mm_add_sd(c41_0, _mm_mul_sd(a41_0, _mm256_castpd256_pd128(b41)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c41_0 = _mm_add_sd(c41_0, _mm_mul_sd(a41_0, b41));
#endif
    _mm_store_sd(&C[(l_n*84)+20], c41_0);
    __m128d c41_1 = _mm_load_sd(&C[(l_n*84)+41]);
    __m128d a41_1 = _mm_load_sd(&A[188]);
#if defined(__SSE3__) && defined(__AVX__)
    c41_1 = _mm_add_sd(c41_1, _mm_mul_sd(a41_1, _mm256_castpd256_pd128(b41)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c41_1 = _mm_add_sd(c41_1, _mm_mul_sd(a41_1, b41));
#endif
    _mm_store_sd(&C[(l_n*84)+41], c41_1);
    __m128d c41_2 = _mm_load_sd(&C[(l_n*84)+69]);
    __m128d a41_2 = _mm_load_sd(&A[189]);
#if defined(__SSE3__) && defined(__AVX__)
    c41_2 = _mm_add_sd(c41_2, _mm_mul_sd(a41_2, _mm256_castpd256_pd128(b41)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c41_2 = _mm_add_sd(c41_2, _mm_mul_sd(a41_2, b41));
#endif
    _mm_store_sd(&C[(l_n*84)+69], c41_2);
#else
    C[(l_n*84)+20] += A[187] * B[(l_n*84)+41];
    C[(l_n*84)+41] += A[188] * B[(l_n*84)+41];
    C[(l_n*84)+69] += A[189] * B[(l_n*84)+41];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b42 = _mm256_broadcast_sd(&B[(l_n*84)+42]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b42 = _mm_loaddup_pd(&B[(l_n*84)+42]);
#endif
    __m128d c42_0 = _mm_load_sd(&C[(l_n*84)+21]);
    __m128d a42_0 = _mm_load_sd(&A[190]);
#if defined(__SSE3__) && defined(__AVX__)
    c42_0 = _mm_add_sd(c42_0, _mm_mul_sd(a42_0, _mm256_castpd256_pd128(b42)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c42_0 = _mm_add_sd(c42_0, _mm_mul_sd(a42_0, b42));
#endif
    _mm_store_sd(&C[(l_n*84)+21], c42_0);
    __m128d c42_1 = _mm_load_sd(&C[(l_n*84)+42]);
    __m128d a42_1 = _mm_load_sd(&A[191]);
#if defined(__SSE3__) && defined(__AVX__)
    c42_1 = _mm_add_sd(c42_1, _mm_mul_sd(a42_1, _mm256_castpd256_pd128(b42)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c42_1 = _mm_add_sd(c42_1, _mm_mul_sd(a42_1, b42));
#endif
    _mm_store_sd(&C[(l_n*84)+42], c42_1);
    __m128d c42_2 = _mm_load_sd(&C[(l_n*84)+70]);
    __m128d a42_2 = _mm_load_sd(&A[192]);
#if defined(__SSE3__) && defined(__AVX__)
    c42_2 = _mm_add_sd(c42_2, _mm_mul_sd(a42_2, _mm256_castpd256_pd128(b42)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c42_2 = _mm_add_sd(c42_2, _mm_mul_sd(a42_2, b42));
#endif
    _mm_store_sd(&C[(l_n*84)+70], c42_2);
#else
    C[(l_n*84)+21] += A[190] * B[(l_n*84)+42];
    C[(l_n*84)+42] += A[191] * B[(l_n*84)+42];
    C[(l_n*84)+70] += A[192] * B[(l_n*84)+42];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b43 = _mm256_broadcast_sd(&B[(l_n*84)+43]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b43 = _mm_loaddup_pd(&B[(l_n*84)+43]);
#endif
    __m128d c43_0 = _mm_load_sd(&C[(l_n*84)+22]);
    __m128d a43_0 = _mm_load_sd(&A[193]);
#if defined(__SSE3__) && defined(__AVX__)
    c43_0 = _mm_add_sd(c43_0, _mm_mul_sd(a43_0, _mm256_castpd256_pd128(b43)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c43_0 = _mm_add_sd(c43_0, _mm_mul_sd(a43_0, b43));
#endif
    _mm_store_sd(&C[(l_n*84)+22], c43_0);
    __m128d c43_1 = _mm_load_sd(&C[(l_n*84)+43]);
    __m128d a43_1 = _mm_load_sd(&A[194]);
#if defined(__SSE3__) && defined(__AVX__)
    c43_1 = _mm_add_sd(c43_1, _mm_mul_sd(a43_1, _mm256_castpd256_pd128(b43)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c43_1 = _mm_add_sd(c43_1, _mm_mul_sd(a43_1, b43));
#endif
    _mm_store_sd(&C[(l_n*84)+43], c43_1);
    __m128d c43_2 = _mm_load_sd(&C[(l_n*84)+71]);
    __m128d a43_2 = _mm_load_sd(&A[195]);
#if defined(__SSE3__) && defined(__AVX__)
    c43_2 = _mm_add_sd(c43_2, _mm_mul_sd(a43_2, _mm256_castpd256_pd128(b43)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c43_2 = _mm_add_sd(c43_2, _mm_mul_sd(a43_2, b43));
#endif
    _mm_store_sd(&C[(l_n*84)+71], c43_2);
#else
    C[(l_n*84)+22] += A[193] * B[(l_n*84)+43];
    C[(l_n*84)+43] += A[194] * B[(l_n*84)+43];
    C[(l_n*84)+71] += A[195] * B[(l_n*84)+43];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b44 = _mm256_broadcast_sd(&B[(l_n*84)+44]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b44 = _mm_loaddup_pd(&B[(l_n*84)+44]);
#endif
    __m128d c44_0 = _mm_load_sd(&C[(l_n*84)+23]);
    __m128d a44_0 = _mm_load_sd(&A[196]);
#if defined(__SSE3__) && defined(__AVX__)
    c44_0 = _mm_add_sd(c44_0, _mm_mul_sd(a44_0, _mm256_castpd256_pd128(b44)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c44_0 = _mm_add_sd(c44_0, _mm_mul_sd(a44_0, b44));
#endif
    _mm_store_sd(&C[(l_n*84)+23], c44_0);
    __m128d c44_1 = _mm_load_sd(&C[(l_n*84)+44]);
    __m128d a44_1 = _mm_load_sd(&A[197]);
#if defined(__SSE3__) && defined(__AVX__)
    c44_1 = _mm_add_sd(c44_1, _mm_mul_sd(a44_1, _mm256_castpd256_pd128(b44)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c44_1 = _mm_add_sd(c44_1, _mm_mul_sd(a44_1, b44));
#endif
    _mm_store_sd(&C[(l_n*84)+44], c44_1);
    __m128d c44_2 = _mm_load_sd(&C[(l_n*84)+72]);
    __m128d a44_2 = _mm_load_sd(&A[198]);
#if defined(__SSE3__) && defined(__AVX__)
    c44_2 = _mm_add_sd(c44_2, _mm_mul_sd(a44_2, _mm256_castpd256_pd128(b44)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c44_2 = _mm_add_sd(c44_2, _mm_mul_sd(a44_2, b44));
#endif
    _mm_store_sd(&C[(l_n*84)+72], c44_2);
#else
    C[(l_n*84)+23] += A[196] * B[(l_n*84)+44];
    C[(l_n*84)+44] += A[197] * B[(l_n*84)+44];
    C[(l_n*84)+72] += A[198] * B[(l_n*84)+44];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b45 = _mm256_broadcast_sd(&B[(l_n*84)+45]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b45 = _mm_loaddup_pd(&B[(l_n*84)+45]);
#endif
    __m128d c45_0 = _mm_load_sd(&C[(l_n*84)+24]);
    __m128d a45_0 = _mm_load_sd(&A[199]);
#if defined(__SSE3__) && defined(__AVX__)
    c45_0 = _mm_add_sd(c45_0, _mm_mul_sd(a45_0, _mm256_castpd256_pd128(b45)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c45_0 = _mm_add_sd(c45_0, _mm_mul_sd(a45_0, b45));
#endif
    _mm_store_sd(&C[(l_n*84)+24], c45_0);
    __m128d c45_1 = _mm_load_sd(&C[(l_n*84)+45]);
    __m128d a45_1 = _mm_load_sd(&A[200]);
#if defined(__SSE3__) && defined(__AVX__)
    c45_1 = _mm_add_sd(c45_1, _mm_mul_sd(a45_1, _mm256_castpd256_pd128(b45)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c45_1 = _mm_add_sd(c45_1, _mm_mul_sd(a45_1, b45));
#endif
    _mm_store_sd(&C[(l_n*84)+45], c45_1);
    __m128d c45_2 = _mm_load_sd(&C[(l_n*84)+73]);
    __m128d a45_2 = _mm_load_sd(&A[201]);
#if defined(__SSE3__) && defined(__AVX__)
    c45_2 = _mm_add_sd(c45_2, _mm_mul_sd(a45_2, _mm256_castpd256_pd128(b45)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c45_2 = _mm_add_sd(c45_2, _mm_mul_sd(a45_2, b45));
#endif
    _mm_store_sd(&C[(l_n*84)+73], c45_2);
#else
    C[(l_n*84)+24] += A[199] * B[(l_n*84)+45];
    C[(l_n*84)+45] += A[200] * B[(l_n*84)+45];
    C[(l_n*84)+73] += A[201] * B[(l_n*84)+45];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b46 = _mm256_broadcast_sd(&B[(l_n*84)+46]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b46 = _mm_loaddup_pd(&B[(l_n*84)+46]);
#endif
    __m128d c46_0 = _mm_load_sd(&C[(l_n*84)+10]);
    __m128d a46_0 = _mm_load_sd(&A[202]);
#if defined(__SSE3__) && defined(__AVX__)
    c46_0 = _mm_add_sd(c46_0, _mm_mul_sd(a46_0, _mm256_castpd256_pd128(b46)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c46_0 = _mm_add_sd(c46_0, _mm_mul_sd(a46_0, b46));
#endif
    _mm_store_sd(&C[(l_n*84)+10], c46_0);
    __m128d c46_1 = _mm_load_sd(&C[(l_n*84)+25]);
    __m128d a46_1 = _mm_load_sd(&A[203]);
#if defined(__SSE3__) && defined(__AVX__)
    c46_1 = _mm_add_sd(c46_1, _mm_mul_sd(a46_1, _mm256_castpd256_pd128(b46)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c46_1 = _mm_add_sd(c46_1, _mm_mul_sd(a46_1, b46));
#endif
    _mm_store_sd(&C[(l_n*84)+25], c46_1);
    __m128d c46_2 = _mm_load_sd(&C[(l_n*84)+46]);
    __m128d a46_2 = _mm_load_sd(&A[204]);
#if defined(__SSE3__) && defined(__AVX__)
    c46_2 = _mm_add_sd(c46_2, _mm_mul_sd(a46_2, _mm256_castpd256_pd128(b46)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c46_2 = _mm_add_sd(c46_2, _mm_mul_sd(a46_2, b46));
#endif
    _mm_store_sd(&C[(l_n*84)+46], c46_2);
    __m128d c46_3 = _mm_load_sd(&C[(l_n*84)+74]);
    __m128d a46_3 = _mm_load_sd(&A[205]);
#if defined(__SSE3__) && defined(__AVX__)
    c46_3 = _mm_add_sd(c46_3, _mm_mul_sd(a46_3, _mm256_castpd256_pd128(b46)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c46_3 = _mm_add_sd(c46_3, _mm_mul_sd(a46_3, b46));
#endif
    _mm_store_sd(&C[(l_n*84)+74], c46_3);
#else
    C[(l_n*84)+10] += A[202] * B[(l_n*84)+46];
    C[(l_n*84)+25] += A[203] * B[(l_n*84)+46];
    C[(l_n*84)+46] += A[204] * B[(l_n*84)+46];
    C[(l_n*84)+74] += A[205] * B[(l_n*84)+46];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b47 = _mm256_broadcast_sd(&B[(l_n*84)+47]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b47 = _mm_loaddup_pd(&B[(l_n*84)+47]);
#endif
    __m128d c47_0 = _mm_load_sd(&C[(l_n*84)+11]);
    __m128d a47_0 = _mm_load_sd(&A[206]);
#if defined(__SSE3__) && defined(__AVX__)
    c47_0 = _mm_add_sd(c47_0, _mm_mul_sd(a47_0, _mm256_castpd256_pd128(b47)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c47_0 = _mm_add_sd(c47_0, _mm_mul_sd(a47_0, b47));
#endif
    _mm_store_sd(&C[(l_n*84)+11], c47_0);
    __m128d c47_1 = _mm_load_sd(&C[(l_n*84)+26]);
    __m128d a47_1 = _mm_load_sd(&A[207]);
#if defined(__SSE3__) && defined(__AVX__)
    c47_1 = _mm_add_sd(c47_1, _mm_mul_sd(a47_1, _mm256_castpd256_pd128(b47)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c47_1 = _mm_add_sd(c47_1, _mm_mul_sd(a47_1, b47));
#endif
    _mm_store_sd(&C[(l_n*84)+26], c47_1);
    __m128d c47_2 = _mm_load_sd(&C[(l_n*84)+47]);
    __m128d a47_2 = _mm_load_sd(&A[208]);
#if defined(__SSE3__) && defined(__AVX__)
    c47_2 = _mm_add_sd(c47_2, _mm_mul_sd(a47_2, _mm256_castpd256_pd128(b47)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c47_2 = _mm_add_sd(c47_2, _mm_mul_sd(a47_2, b47));
#endif
    _mm_store_sd(&C[(l_n*84)+47], c47_2);
    __m128d c47_3 = _mm_load_sd(&C[(l_n*84)+75]);
    __m128d a47_3 = _mm_load_sd(&A[209]);
#if defined(__SSE3__) && defined(__AVX__)
    c47_3 = _mm_add_sd(c47_3, _mm_mul_sd(a47_3, _mm256_castpd256_pd128(b47)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c47_3 = _mm_add_sd(c47_3, _mm_mul_sd(a47_3, b47));
#endif
    _mm_store_sd(&C[(l_n*84)+75], c47_3);
#else
    C[(l_n*84)+11] += A[206] * B[(l_n*84)+47];
    C[(l_n*84)+26] += A[207] * B[(l_n*84)+47];
    C[(l_n*84)+47] += A[208] * B[(l_n*84)+47];
    C[(l_n*84)+75] += A[209] * B[(l_n*84)+47];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b48 = _mm256_broadcast_sd(&B[(l_n*84)+48]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b48 = _mm_loaddup_pd(&B[(l_n*84)+48]);
#endif
    __m128d c48_0 = _mm_load_sd(&C[(l_n*84)+12]);
    __m128d a48_0 = _mm_load_sd(&A[210]);
#if defined(__SSE3__) && defined(__AVX__)
    c48_0 = _mm_add_sd(c48_0, _mm_mul_sd(a48_0, _mm256_castpd256_pd128(b48)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c48_0 = _mm_add_sd(c48_0, _mm_mul_sd(a48_0, b48));
#endif
    _mm_store_sd(&C[(l_n*84)+12], c48_0);
    __m128d c48_1 = _mm_load_sd(&C[(l_n*84)+27]);
    __m128d a48_1 = _mm_load_sd(&A[211]);
#if defined(__SSE3__) && defined(__AVX__)
    c48_1 = _mm_add_sd(c48_1, _mm_mul_sd(a48_1, _mm256_castpd256_pd128(b48)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c48_1 = _mm_add_sd(c48_1, _mm_mul_sd(a48_1, b48));
#endif
    _mm_store_sd(&C[(l_n*84)+27], c48_1);
    __m128d c48_2 = _mm_load_sd(&C[(l_n*84)+48]);
    __m128d a48_2 = _mm_load_sd(&A[212]);
#if defined(__SSE3__) && defined(__AVX__)
    c48_2 = _mm_add_sd(c48_2, _mm_mul_sd(a48_2, _mm256_castpd256_pd128(b48)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c48_2 = _mm_add_sd(c48_2, _mm_mul_sd(a48_2, b48));
#endif
    _mm_store_sd(&C[(l_n*84)+48], c48_2);
    __m128d c48_3 = _mm_load_sd(&C[(l_n*84)+76]);
    __m128d a48_3 = _mm_load_sd(&A[213]);
#if defined(__SSE3__) && defined(__AVX__)
    c48_3 = _mm_add_sd(c48_3, _mm_mul_sd(a48_3, _mm256_castpd256_pd128(b48)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c48_3 = _mm_add_sd(c48_3, _mm_mul_sd(a48_3, b48));
#endif
    _mm_store_sd(&C[(l_n*84)+76], c48_3);
#else
    C[(l_n*84)+12] += A[210] * B[(l_n*84)+48];
    C[(l_n*84)+27] += A[211] * B[(l_n*84)+48];
    C[(l_n*84)+48] += A[212] * B[(l_n*84)+48];
    C[(l_n*84)+76] += A[213] * B[(l_n*84)+48];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b49 = _mm256_broadcast_sd(&B[(l_n*84)+49]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b49 = _mm_loaddup_pd(&B[(l_n*84)+49]);
#endif
    __m128d c49_0 = _mm_load_sd(&C[(l_n*84)+13]);
    __m128d a49_0 = _mm_load_sd(&A[214]);
#if defined(__SSE3__) && defined(__AVX__)
    c49_0 = _mm_add_sd(c49_0, _mm_mul_sd(a49_0, _mm256_castpd256_pd128(b49)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c49_0 = _mm_add_sd(c49_0, _mm_mul_sd(a49_0, b49));
#endif
    _mm_store_sd(&C[(l_n*84)+13], c49_0);
    __m128d c49_1 = _mm_load_sd(&C[(l_n*84)+28]);
    __m128d a49_1 = _mm_load_sd(&A[215]);
#if defined(__SSE3__) && defined(__AVX__)
    c49_1 = _mm_add_sd(c49_1, _mm_mul_sd(a49_1, _mm256_castpd256_pd128(b49)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c49_1 = _mm_add_sd(c49_1, _mm_mul_sd(a49_1, b49));
#endif
    _mm_store_sd(&C[(l_n*84)+28], c49_1);
    __m128d c49_2 = _mm_load_sd(&C[(l_n*84)+49]);
    __m128d a49_2 = _mm_load_sd(&A[216]);
#if defined(__SSE3__) && defined(__AVX__)
    c49_2 = _mm_add_sd(c49_2, _mm_mul_sd(a49_2, _mm256_castpd256_pd128(b49)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c49_2 = _mm_add_sd(c49_2, _mm_mul_sd(a49_2, b49));
#endif
    _mm_store_sd(&C[(l_n*84)+49], c49_2);
    __m128d c49_3 = _mm_load_sd(&C[(l_n*84)+77]);
    __m128d a49_3 = _mm_load_sd(&A[217]);
#if defined(__SSE3__) && defined(__AVX__)
    c49_3 = _mm_add_sd(c49_3, _mm_mul_sd(a49_3, _mm256_castpd256_pd128(b49)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c49_3 = _mm_add_sd(c49_3, _mm_mul_sd(a49_3, b49));
#endif
    _mm_store_sd(&C[(l_n*84)+77], c49_3);
#else
    C[(l_n*84)+13] += A[214] * B[(l_n*84)+49];
    C[(l_n*84)+28] += A[215] * B[(l_n*84)+49];
    C[(l_n*84)+49] += A[216] * B[(l_n*84)+49];
    C[(l_n*84)+77] += A[217] * B[(l_n*84)+49];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b50 = _mm256_broadcast_sd(&B[(l_n*84)+50]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b50 = _mm_loaddup_pd(&B[(l_n*84)+50]);
#endif
    __m128d c50_0 = _mm_load_sd(&C[(l_n*84)+4]);
    __m128d a50_0 = _mm_load_sd(&A[218]);
#if defined(__SSE3__) && defined(__AVX__)
    c50_0 = _mm_add_sd(c50_0, _mm_mul_sd(a50_0, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c50_0 = _mm_add_sd(c50_0, _mm_mul_sd(a50_0, b50));
#endif
    _mm_store_sd(&C[(l_n*84)+4], c50_0);
    __m128d c50_1 = _mm_load_sd(&C[(l_n*84)+14]);
    __m128d a50_1 = _mm_load_sd(&A[219]);
#if defined(__SSE3__) && defined(__AVX__)
    c50_1 = _mm_add_sd(c50_1, _mm_mul_sd(a50_1, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c50_1 = _mm_add_sd(c50_1, _mm_mul_sd(a50_1, b50));
#endif
    _mm_store_sd(&C[(l_n*84)+14], c50_1);
    __m128d c50_2 = _mm_load_sd(&C[(l_n*84)+29]);
    __m128d a50_2 = _mm_load_sd(&A[220]);
#if defined(__SSE3__) && defined(__AVX__)
    c50_2 = _mm_add_sd(c50_2, _mm_mul_sd(a50_2, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c50_2 = _mm_add_sd(c50_2, _mm_mul_sd(a50_2, b50));
#endif
    _mm_store_sd(&C[(l_n*84)+29], c50_2);
    __m128d c50_3 = _mm_load_sd(&C[(l_n*84)+50]);
    __m128d a50_3 = _mm_load_sd(&A[221]);
#if defined(__SSE3__) && defined(__AVX__)
    c50_3 = _mm_add_sd(c50_3, _mm_mul_sd(a50_3, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c50_3 = _mm_add_sd(c50_3, _mm_mul_sd(a50_3, b50));
#endif
    _mm_store_sd(&C[(l_n*84)+50], c50_3);
    __m128d c50_4 = _mm_load_sd(&C[(l_n*84)+78]);
    __m128d a50_4 = _mm_load_sd(&A[222]);
#if defined(__SSE3__) && defined(__AVX__)
    c50_4 = _mm_add_sd(c50_4, _mm_mul_sd(a50_4, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c50_4 = _mm_add_sd(c50_4, _mm_mul_sd(a50_4, b50));
#endif
    _mm_store_sd(&C[(l_n*84)+78], c50_4);
#else
    C[(l_n*84)+4] += A[218] * B[(l_n*84)+50];
    C[(l_n*84)+14] += A[219] * B[(l_n*84)+50];
    C[(l_n*84)+29] += A[220] * B[(l_n*84)+50];
    C[(l_n*84)+50] += A[221] * B[(l_n*84)+50];
    C[(l_n*84)+78] += A[222] * B[(l_n*84)+50];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b51 = _mm256_broadcast_sd(&B[(l_n*84)+51]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b51 = _mm_loaddup_pd(&B[(l_n*84)+51]);
#endif
    __m128d c51_0 = _mm_load_sd(&C[(l_n*84)+5]);
    __m128d a51_0 = _mm_load_sd(&A[223]);
#if defined(__SSE3__) && defined(__AVX__)
    c51_0 = _mm_add_sd(c51_0, _mm_mul_sd(a51_0, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c51_0 = _mm_add_sd(c51_0, _mm_mul_sd(a51_0, b51));
#endif
    _mm_store_sd(&C[(l_n*84)+5], c51_0);
    __m128d c51_1 = _mm_load_sd(&C[(l_n*84)+15]);
    __m128d a51_1 = _mm_load_sd(&A[224]);
#if defined(__SSE3__) && defined(__AVX__)
    c51_1 = _mm_add_sd(c51_1, _mm_mul_sd(a51_1, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c51_1 = _mm_add_sd(c51_1, _mm_mul_sd(a51_1, b51));
#endif
    _mm_store_sd(&C[(l_n*84)+15], c51_1);
    __m128d c51_2 = _mm_load_sd(&C[(l_n*84)+30]);
    __m128d a51_2 = _mm_load_sd(&A[225]);
#if defined(__SSE3__) && defined(__AVX__)
    c51_2 = _mm_add_sd(c51_2, _mm_mul_sd(a51_2, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c51_2 = _mm_add_sd(c51_2, _mm_mul_sd(a51_2, b51));
#endif
    _mm_store_sd(&C[(l_n*84)+30], c51_2);
    __m128d c51_3 = _mm_load_sd(&C[(l_n*84)+51]);
    __m128d a51_3 = _mm_load_sd(&A[226]);
#if defined(__SSE3__) && defined(__AVX__)
    c51_3 = _mm_add_sd(c51_3, _mm_mul_sd(a51_3, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c51_3 = _mm_add_sd(c51_3, _mm_mul_sd(a51_3, b51));
#endif
    _mm_store_sd(&C[(l_n*84)+51], c51_3);
    __m128d c51_4 = _mm_load_sd(&C[(l_n*84)+79]);
    __m128d a51_4 = _mm_load_sd(&A[227]);
#if defined(__SSE3__) && defined(__AVX__)
    c51_4 = _mm_add_sd(c51_4, _mm_mul_sd(a51_4, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c51_4 = _mm_add_sd(c51_4, _mm_mul_sd(a51_4, b51));
#endif
    _mm_store_sd(&C[(l_n*84)+79], c51_4);
#else
    C[(l_n*84)+5] += A[223] * B[(l_n*84)+51];
    C[(l_n*84)+15] += A[224] * B[(l_n*84)+51];
    C[(l_n*84)+30] += A[225] * B[(l_n*84)+51];
    C[(l_n*84)+51] += A[226] * B[(l_n*84)+51];
    C[(l_n*84)+79] += A[227] * B[(l_n*84)+51];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b52 = _mm256_broadcast_sd(&B[(l_n*84)+52]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b52 = _mm_loaddup_pd(&B[(l_n*84)+52]);
#endif
    __m128d c52_0 = _mm_load_sd(&C[(l_n*84)+6]);
    __m128d a52_0 = _mm_load_sd(&A[228]);
#if defined(__SSE3__) && defined(__AVX__)
    c52_0 = _mm_add_sd(c52_0, _mm_mul_sd(a52_0, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c52_0 = _mm_add_sd(c52_0, _mm_mul_sd(a52_0, b52));
#endif
    _mm_store_sd(&C[(l_n*84)+6], c52_0);
    __m128d c52_1 = _mm_load_sd(&C[(l_n*84)+16]);
    __m128d a52_1 = _mm_load_sd(&A[229]);
#if defined(__SSE3__) && defined(__AVX__)
    c52_1 = _mm_add_sd(c52_1, _mm_mul_sd(a52_1, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c52_1 = _mm_add_sd(c52_1, _mm_mul_sd(a52_1, b52));
#endif
    _mm_store_sd(&C[(l_n*84)+16], c52_1);
    __m128d c52_2 = _mm_load_sd(&C[(l_n*84)+31]);
    __m128d a52_2 = _mm_load_sd(&A[230]);
#if defined(__SSE3__) && defined(__AVX__)
    c52_2 = _mm_add_sd(c52_2, _mm_mul_sd(a52_2, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c52_2 = _mm_add_sd(c52_2, _mm_mul_sd(a52_2, b52));
#endif
    _mm_store_sd(&C[(l_n*84)+31], c52_2);
    __m128d c52_3 = _mm_load_sd(&C[(l_n*84)+52]);
    __m128d a52_3 = _mm_load_sd(&A[231]);
#if defined(__SSE3__) && defined(__AVX__)
    c52_3 = _mm_add_sd(c52_3, _mm_mul_sd(a52_3, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c52_3 = _mm_add_sd(c52_3, _mm_mul_sd(a52_3, b52));
#endif
    _mm_store_sd(&C[(l_n*84)+52], c52_3);
    __m128d c52_4 = _mm_load_sd(&C[(l_n*84)+80]);
    __m128d a52_4 = _mm_load_sd(&A[232]);
#if defined(__SSE3__) && defined(__AVX__)
    c52_4 = _mm_add_sd(c52_4, _mm_mul_sd(a52_4, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c52_4 = _mm_add_sd(c52_4, _mm_mul_sd(a52_4, b52));
#endif
    _mm_store_sd(&C[(l_n*84)+80], c52_4);
#else
    C[(l_n*84)+6] += A[228] * B[(l_n*84)+52];
    C[(l_n*84)+16] += A[229] * B[(l_n*84)+52];
    C[(l_n*84)+31] += A[230] * B[(l_n*84)+52];
    C[(l_n*84)+52] += A[231] * B[(l_n*84)+52];
    C[(l_n*84)+80] += A[232] * B[(l_n*84)+52];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b53 = _mm256_broadcast_sd(&B[(l_n*84)+53]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b53 = _mm_loaddup_pd(&B[(l_n*84)+53]);
#endif
    __m128d c53_0 = _mm_load_sd(&C[(l_n*84)+1]);
    __m128d a53_0 = _mm_load_sd(&A[233]);
#if defined(__SSE3__) && defined(__AVX__)
    c53_0 = _mm_add_sd(c53_0, _mm_mul_sd(a53_0, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c53_0 = _mm_add_sd(c53_0, _mm_mul_sd(a53_0, b53));
#endif
    _mm_store_sd(&C[(l_n*84)+1], c53_0);
    __m128d c53_1 = _mm_load_sd(&C[(l_n*84)+7]);
    __m128d a53_1 = _mm_load_sd(&A[234]);
#if defined(__SSE3__) && defined(__AVX__)
    c53_1 = _mm_add_sd(c53_1, _mm_mul_sd(a53_1, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c53_1 = _mm_add_sd(c53_1, _mm_mul_sd(a53_1, b53));
#endif
    _mm_store_sd(&C[(l_n*84)+7], c53_1);
    __m128d c53_2 = _mm_load_sd(&C[(l_n*84)+17]);
    __m128d a53_2 = _mm_load_sd(&A[235]);
#if defined(__SSE3__) && defined(__AVX__)
    c53_2 = _mm_add_sd(c53_2, _mm_mul_sd(a53_2, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c53_2 = _mm_add_sd(c53_2, _mm_mul_sd(a53_2, b53));
#endif
    _mm_store_sd(&C[(l_n*84)+17], c53_2);
    __m128d c53_3 = _mm_load_sd(&C[(l_n*84)+32]);
    __m128d a53_3 = _mm_load_sd(&A[236]);
#if defined(__SSE3__) && defined(__AVX__)
    c53_3 = _mm_add_sd(c53_3, _mm_mul_sd(a53_3, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c53_3 = _mm_add_sd(c53_3, _mm_mul_sd(a53_3, b53));
#endif
    _mm_store_sd(&C[(l_n*84)+32], c53_3);
    __m128d c53_4 = _mm_load_sd(&C[(l_n*84)+53]);
    __m128d a53_4 = _mm_load_sd(&A[237]);
#if defined(__SSE3__) && defined(__AVX__)
    c53_4 = _mm_add_sd(c53_4, _mm_mul_sd(a53_4, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c53_4 = _mm_add_sd(c53_4, _mm_mul_sd(a53_4, b53));
#endif
    _mm_store_sd(&C[(l_n*84)+53], c53_4);
    __m128d c53_5 = _mm_load_sd(&C[(l_n*84)+81]);
    __m128d a53_5 = _mm_load_sd(&A[238]);
#if defined(__SSE3__) && defined(__AVX__)
    c53_5 = _mm_add_sd(c53_5, _mm_mul_sd(a53_5, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c53_5 = _mm_add_sd(c53_5, _mm_mul_sd(a53_5, b53));
#endif
    _mm_store_sd(&C[(l_n*84)+81], c53_5);
#else
    C[(l_n*84)+1] += A[233] * B[(l_n*84)+53];
    C[(l_n*84)+7] += A[234] * B[(l_n*84)+53];
    C[(l_n*84)+17] += A[235] * B[(l_n*84)+53];
    C[(l_n*84)+32] += A[236] * B[(l_n*84)+53];
    C[(l_n*84)+53] += A[237] * B[(l_n*84)+53];
    C[(l_n*84)+81] += A[238] * B[(l_n*84)+53];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b54 = _mm256_broadcast_sd(&B[(l_n*84)+54]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b54 = _mm_loaddup_pd(&B[(l_n*84)+54]);
#endif
    __m128d c54_0 = _mm_load_sd(&C[(l_n*84)+2]);
    __m128d a54_0 = _mm_load_sd(&A[239]);
#if defined(__SSE3__) && defined(__AVX__)
    c54_0 = _mm_add_sd(c54_0, _mm_mul_sd(a54_0, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c54_0 = _mm_add_sd(c54_0, _mm_mul_sd(a54_0, b54));
#endif
    _mm_store_sd(&C[(l_n*84)+2], c54_0);
    __m128d c54_1 = _mm_load_sd(&C[(l_n*84)+8]);
    __m128d a54_1 = _mm_load_sd(&A[240]);
#if defined(__SSE3__) && defined(__AVX__)
    c54_1 = _mm_add_sd(c54_1, _mm_mul_sd(a54_1, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c54_1 = _mm_add_sd(c54_1, _mm_mul_sd(a54_1, b54));
#endif
    _mm_store_sd(&C[(l_n*84)+8], c54_1);
    __m128d c54_2 = _mm_load_sd(&C[(l_n*84)+18]);
    __m128d a54_2 = _mm_load_sd(&A[241]);
#if defined(__SSE3__) && defined(__AVX__)
    c54_2 = _mm_add_sd(c54_2, _mm_mul_sd(a54_2, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c54_2 = _mm_add_sd(c54_2, _mm_mul_sd(a54_2, b54));
#endif
    _mm_store_sd(&C[(l_n*84)+18], c54_2);
    __m128d c54_3 = _mm_load_sd(&C[(l_n*84)+33]);
    __m128d a54_3 = _mm_load_sd(&A[242]);
#if defined(__SSE3__) && defined(__AVX__)
    c54_3 = _mm_add_sd(c54_3, _mm_mul_sd(a54_3, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c54_3 = _mm_add_sd(c54_3, _mm_mul_sd(a54_3, b54));
#endif
    _mm_store_sd(&C[(l_n*84)+33], c54_3);
    __m128d c54_4 = _mm_load_sd(&C[(l_n*84)+54]);
    __m128d a54_4 = _mm_load_sd(&A[243]);
#if defined(__SSE3__) && defined(__AVX__)
    c54_4 = _mm_add_sd(c54_4, _mm_mul_sd(a54_4, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c54_4 = _mm_add_sd(c54_4, _mm_mul_sd(a54_4, b54));
#endif
    _mm_store_sd(&C[(l_n*84)+54], c54_4);
    __m128d c54_5 = _mm_load_sd(&C[(l_n*84)+82]);
    __m128d a54_5 = _mm_load_sd(&A[244]);
#if defined(__SSE3__) && defined(__AVX__)
    c54_5 = _mm_add_sd(c54_5, _mm_mul_sd(a54_5, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c54_5 = _mm_add_sd(c54_5, _mm_mul_sd(a54_5, b54));
#endif
    _mm_store_sd(&C[(l_n*84)+82], c54_5);
#else
    C[(l_n*84)+2] += A[239] * B[(l_n*84)+54];
    C[(l_n*84)+8] += A[240] * B[(l_n*84)+54];
    C[(l_n*84)+18] += A[241] * B[(l_n*84)+54];
    C[(l_n*84)+33] += A[242] * B[(l_n*84)+54];
    C[(l_n*84)+54] += A[243] * B[(l_n*84)+54];
    C[(l_n*84)+82] += A[244] * B[(l_n*84)+54];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b55 = _mm256_broadcast_sd(&B[(l_n*84)+55]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b55 = _mm_loaddup_pd(&B[(l_n*84)+55]);
#endif
    __m128d c55_0 = _mm_load_sd(&C[(l_n*84)+0]);
    __m128d a55_0 = _mm_load_sd(&A[245]);
#if defined(__SSE3__) && defined(__AVX__)
    c55_0 = _mm_add_sd(c55_0, _mm_mul_sd(a55_0, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c55_0 = _mm_add_sd(c55_0, _mm_mul_sd(a55_0, b55));
#endif
    _mm_store_sd(&C[(l_n*84)+0], c55_0);
    __m128d c55_1 = _mm_load_sd(&C[(l_n*84)+3]);
    __m128d a55_1 = _mm_load_sd(&A[246]);
#if defined(__SSE3__) && defined(__AVX__)
    c55_1 = _mm_add_sd(c55_1, _mm_mul_sd(a55_1, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c55_1 = _mm_add_sd(c55_1, _mm_mul_sd(a55_1, b55));
#endif
    _mm_store_sd(&C[(l_n*84)+3], c55_1);
    __m128d c55_2 = _mm_load_sd(&C[(l_n*84)+9]);
    __m128d a55_2 = _mm_load_sd(&A[247]);
#if defined(__SSE3__) && defined(__AVX__)
    c55_2 = _mm_add_sd(c55_2, _mm_mul_sd(a55_2, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c55_2 = _mm_add_sd(c55_2, _mm_mul_sd(a55_2, b55));
#endif
    _mm_store_sd(&C[(l_n*84)+9], c55_2);
    __m128d c55_3 = _mm_load_sd(&C[(l_n*84)+19]);
    __m128d a55_3 = _mm_load_sd(&A[248]);
#if defined(__SSE3__) && defined(__AVX__)
    c55_3 = _mm_add_sd(c55_3, _mm_mul_sd(a55_3, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c55_3 = _mm_add_sd(c55_3, _mm_mul_sd(a55_3, b55));
#endif
    _mm_store_sd(&C[(l_n*84)+19], c55_3);
    __m128d c55_4 = _mm_load_sd(&C[(l_n*84)+34]);
    __m128d a55_4 = _mm_load_sd(&A[249]);
#if defined(__SSE3__) && defined(__AVX__)
    c55_4 = _mm_add_sd(c55_4, _mm_mul_sd(a55_4, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c55_4 = _mm_add_sd(c55_4, _mm_mul_sd(a55_4, b55));
#endif
    _mm_store_sd(&C[(l_n*84)+34], c55_4);
    __m128d c55_5 = _mm_load_sd(&C[(l_n*84)+55]);
    __m128d a55_5 = _mm_load_sd(&A[250]);
#if defined(__SSE3__) && defined(__AVX__)
    c55_5 = _mm_add_sd(c55_5, _mm_mul_sd(a55_5, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c55_5 = _mm_add_sd(c55_5, _mm_mul_sd(a55_5, b55));
#endif
    _mm_store_sd(&C[(l_n*84)+55], c55_5);
    __m128d c55_6 = _mm_load_sd(&C[(l_n*84)+83]);
    __m128d a55_6 = _mm_load_sd(&A[251]);
#if defined(__SSE3__) && defined(__AVX__)
    c55_6 = _mm_add_sd(c55_6, _mm_mul_sd(a55_6, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c55_6 = _mm_add_sd(c55_6, _mm_mul_sd(a55_6, b55));
#endif
    _mm_store_sd(&C[(l_n*84)+83], c55_6);
#else
    C[(l_n*84)+0] += A[245] * B[(l_n*84)+55];
    C[(l_n*84)+3] += A[246] * B[(l_n*84)+55];
    C[(l_n*84)+9] += A[247] * B[(l_n*84)+55];
    C[(l_n*84)+19] += A[248] * B[(l_n*84)+55];
    C[(l_n*84)+34] += A[249] * B[(l_n*84)+55];
    C[(l_n*84)+55] += A[250] * B[(l_n*84)+55];
    C[(l_n*84)+83] += A[251] * B[(l_n*84)+55];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b56 = _mm256_broadcast_sd(&B[(l_n*84)+56]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b56 = _mm_loaddup_pd(&B[(l_n*84)+56]);
#endif
    __m128d c56_0 = _mm_load_sd(&C[(l_n*84)+56]);
    __m128d a56_0 = _mm_load_sd(&A[252]);
#if defined(__SSE3__) && defined(__AVX__)
    c56_0 = _mm_add_sd(c56_0, _mm_mul_sd(a56_0, _mm256_castpd256_pd128(b56)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c56_0 = _mm_add_sd(c56_0, _mm_mul_sd(a56_0, b56));
#endif
    _mm_store_sd(&C[(l_n*84)+56], c56_0);
#else
    C[(l_n*84)+56] += A[252] * B[(l_n*84)+56];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b57 = _mm256_broadcast_sd(&B[(l_n*84)+57]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b57 = _mm_loaddup_pd(&B[(l_n*84)+57]);
#endif
    __m128d c57_0 = _mm_load_sd(&C[(l_n*84)+57]);
    __m128d a57_0 = _mm_load_sd(&A[253]);
#if defined(__SSE3__) && defined(__AVX__)
    c57_0 = _mm_add_sd(c57_0, _mm_mul_sd(a57_0, _mm256_castpd256_pd128(b57)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c57_0 = _mm_add_sd(c57_0, _mm_mul_sd(a57_0, b57));
#endif
    _mm_store_sd(&C[(l_n*84)+57], c57_0);
#else
    C[(l_n*84)+57] += A[253] * B[(l_n*84)+57];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b58 = _mm256_broadcast_sd(&B[(l_n*84)+58]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b58 = _mm_loaddup_pd(&B[(l_n*84)+58]);
#endif
    __m128d c58_0 = _mm_load_sd(&C[(l_n*84)+58]);
    __m128d a58_0 = _mm_load_sd(&A[254]);
#if defined(__SSE3__) && defined(__AVX__)
    c58_0 = _mm_add_sd(c58_0, _mm_mul_sd(a58_0, _mm256_castpd256_pd128(b58)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c58_0 = _mm_add_sd(c58_0, _mm_mul_sd(a58_0, b58));
#endif
    _mm_store_sd(&C[(l_n*84)+58], c58_0);
#else
    C[(l_n*84)+58] += A[254] * B[(l_n*84)+58];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b59 = _mm256_broadcast_sd(&B[(l_n*84)+59]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b59 = _mm_loaddup_pd(&B[(l_n*84)+59]);
#endif
    __m128d c59_0 = _mm_load_sd(&C[(l_n*84)+59]);
    __m128d a59_0 = _mm_load_sd(&A[255]);
#if defined(__SSE3__) && defined(__AVX__)
    c59_0 = _mm_add_sd(c59_0, _mm_mul_sd(a59_0, _mm256_castpd256_pd128(b59)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c59_0 = _mm_add_sd(c59_0, _mm_mul_sd(a59_0, b59));
#endif
    _mm_store_sd(&C[(l_n*84)+59], c59_0);
#else
    C[(l_n*84)+59] += A[255] * B[(l_n*84)+59];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b60 = _mm256_broadcast_sd(&B[(l_n*84)+60]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b60 = _mm_loaddup_pd(&B[(l_n*84)+60]);
#endif
    __m128d c60_0 = _mm_load_sd(&C[(l_n*84)+60]);
    __m128d a60_0 = _mm_load_sd(&A[256]);
#if defined(__SSE3__) && defined(__AVX__)
    c60_0 = _mm_add_sd(c60_0, _mm_mul_sd(a60_0, _mm256_castpd256_pd128(b60)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c60_0 = _mm_add_sd(c60_0, _mm_mul_sd(a60_0, b60));
#endif
    _mm_store_sd(&C[(l_n*84)+60], c60_0);
#else
    C[(l_n*84)+60] += A[256] * B[(l_n*84)+60];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b61 = _mm256_broadcast_sd(&B[(l_n*84)+61]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b61 = _mm_loaddup_pd(&B[(l_n*84)+61]);
#endif
    __m128d c61_0 = _mm_load_sd(&C[(l_n*84)+61]);
    __m128d a61_0 = _mm_load_sd(&A[257]);
#if defined(__SSE3__) && defined(__AVX__)
    c61_0 = _mm_add_sd(c61_0, _mm_mul_sd(a61_0, _mm256_castpd256_pd128(b61)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c61_0 = _mm_add_sd(c61_0, _mm_mul_sd(a61_0, b61));
#endif
    _mm_store_sd(&C[(l_n*84)+61], c61_0);
#else
    C[(l_n*84)+61] += A[257] * B[(l_n*84)+61];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b62 = _mm256_broadcast_sd(&B[(l_n*84)+62]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b62 = _mm_loaddup_pd(&B[(l_n*84)+62]);
#endif
    __m128d c62_0 = _mm_load_sd(&C[(l_n*84)+62]);
    __m128d a62_0 = _mm_load_sd(&A[258]);
#if defined(__SSE3__) && defined(__AVX__)
    c62_0 = _mm_add_sd(c62_0, _mm_mul_sd(a62_0, _mm256_castpd256_pd128(b62)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c62_0 = _mm_add_sd(c62_0, _mm_mul_sd(a62_0, b62));
#endif
    _mm_store_sd(&C[(l_n*84)+62], c62_0);
#else
    C[(l_n*84)+62] += A[258] * B[(l_n*84)+62];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b63 = _mm256_broadcast_sd(&B[(l_n*84)+63]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b63 = _mm_loaddup_pd(&B[(l_n*84)+63]);
#endif
    __m128d c63_0 = _mm_load_sd(&C[(l_n*84)+35]);
    __m128d a63_0 = _mm_load_sd(&A[259]);
#if defined(__SSE3__) && defined(__AVX__)
    c63_0 = _mm_add_sd(c63_0, _mm_mul_sd(a63_0, _mm256_castpd256_pd128(b63)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c63_0 = _mm_add_sd(c63_0, _mm_mul_sd(a63_0, b63));
#endif
    _mm_store_sd(&C[(l_n*84)+35], c63_0);
    __m128d c63_1 = _mm_load_sd(&C[(l_n*84)+63]);
    __m128d a63_1 = _mm_load_sd(&A[260]);
#if defined(__SSE3__) && defined(__AVX__)
    c63_1 = _mm_add_sd(c63_1, _mm_mul_sd(a63_1, _mm256_castpd256_pd128(b63)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c63_1 = _mm_add_sd(c63_1, _mm_mul_sd(a63_1, b63));
#endif
    _mm_store_sd(&C[(l_n*84)+63], c63_1);
#else
    C[(l_n*84)+35] += A[259] * B[(l_n*84)+63];
    C[(l_n*84)+63] += A[260] * B[(l_n*84)+63];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b64 = _mm256_broadcast_sd(&B[(l_n*84)+64]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b64 = _mm_loaddup_pd(&B[(l_n*84)+64]);
#endif
    __m128d c64_0 = _mm_load_sd(&C[(l_n*84)+36]);
    __m128d a64_0 = _mm_load_sd(&A[261]);
#if defined(__SSE3__) && defined(__AVX__)
    c64_0 = _mm_add_sd(c64_0, _mm_mul_sd(a64_0, _mm256_castpd256_pd128(b64)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c64_0 = _mm_add_sd(c64_0, _mm_mul_sd(a64_0, b64));
#endif
    _mm_store_sd(&C[(l_n*84)+36], c64_0);
    __m128d c64_1 = _mm_load_sd(&C[(l_n*84)+64]);
    __m128d a64_1 = _mm_load_sd(&A[262]);
#if defined(__SSE3__) && defined(__AVX__)
    c64_1 = _mm_add_sd(c64_1, _mm_mul_sd(a64_1, _mm256_castpd256_pd128(b64)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c64_1 = _mm_add_sd(c64_1, _mm_mul_sd(a64_1, b64));
#endif
    _mm_store_sd(&C[(l_n*84)+64], c64_1);
#else
    C[(l_n*84)+36] += A[261] * B[(l_n*84)+64];
    C[(l_n*84)+64] += A[262] * B[(l_n*84)+64];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b65 = _mm256_broadcast_sd(&B[(l_n*84)+65]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b65 = _mm_loaddup_pd(&B[(l_n*84)+65]);
#endif
    __m128d c65_0 = _mm_load_sd(&C[(l_n*84)+37]);
    __m128d a65_0 = _mm_load_sd(&A[263]);
#if defined(__SSE3__) && defined(__AVX__)
    c65_0 = _mm_add_sd(c65_0, _mm_mul_sd(a65_0, _mm256_castpd256_pd128(b65)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c65_0 = _mm_add_sd(c65_0, _mm_mul_sd(a65_0, b65));
#endif
    _mm_store_sd(&C[(l_n*84)+37], c65_0);
    __m128d c65_1 = _mm_load_sd(&C[(l_n*84)+65]);
    __m128d a65_1 = _mm_load_sd(&A[264]);
#if defined(__SSE3__) && defined(__AVX__)
    c65_1 = _mm_add_sd(c65_1, _mm_mul_sd(a65_1, _mm256_castpd256_pd128(b65)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c65_1 = _mm_add_sd(c65_1, _mm_mul_sd(a65_1, b65));
#endif
    _mm_store_sd(&C[(l_n*84)+65], c65_1);
#else
    C[(l_n*84)+37] += A[263] * B[(l_n*84)+65];
    C[(l_n*84)+65] += A[264] * B[(l_n*84)+65];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b66 = _mm256_broadcast_sd(&B[(l_n*84)+66]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b66 = _mm_loaddup_pd(&B[(l_n*84)+66]);
#endif
    __m128d c66_0 = _mm_load_sd(&C[(l_n*84)+38]);
    __m128d a66_0 = _mm_load_sd(&A[265]);
#if defined(__SSE3__) && defined(__AVX__)
    c66_0 = _mm_add_sd(c66_0, _mm_mul_sd(a66_0, _mm256_castpd256_pd128(b66)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c66_0 = _mm_add_sd(c66_0, _mm_mul_sd(a66_0, b66));
#endif
    _mm_store_sd(&C[(l_n*84)+38], c66_0);
    __m128d c66_1 = _mm_load_sd(&C[(l_n*84)+66]);
    __m128d a66_1 = _mm_load_sd(&A[266]);
#if defined(__SSE3__) && defined(__AVX__)
    c66_1 = _mm_add_sd(c66_1, _mm_mul_sd(a66_1, _mm256_castpd256_pd128(b66)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c66_1 = _mm_add_sd(c66_1, _mm_mul_sd(a66_1, b66));
#endif
    _mm_store_sd(&C[(l_n*84)+66], c66_1);
#else
    C[(l_n*84)+38] += A[265] * B[(l_n*84)+66];
    C[(l_n*84)+66] += A[266] * B[(l_n*84)+66];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b67 = _mm256_broadcast_sd(&B[(l_n*84)+67]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b67 = _mm_loaddup_pd(&B[(l_n*84)+67]);
#endif
    __m128d c67_0 = _mm_load_sd(&C[(l_n*84)+39]);
    __m128d a67_0 = _mm_load_sd(&A[267]);
#if defined(__SSE3__) && defined(__AVX__)
    c67_0 = _mm_add_sd(c67_0, _mm_mul_sd(a67_0, _mm256_castpd256_pd128(b67)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c67_0 = _mm_add_sd(c67_0, _mm_mul_sd(a67_0, b67));
#endif
    _mm_store_sd(&C[(l_n*84)+39], c67_0);
    __m128d c67_1 = _mm_load_sd(&C[(l_n*84)+67]);
    __m128d a67_1 = _mm_load_sd(&A[268]);
#if defined(__SSE3__) && defined(__AVX__)
    c67_1 = _mm_add_sd(c67_1, _mm_mul_sd(a67_1, _mm256_castpd256_pd128(b67)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c67_1 = _mm_add_sd(c67_1, _mm_mul_sd(a67_1, b67));
#endif
    _mm_store_sd(&C[(l_n*84)+67], c67_1);
#else
    C[(l_n*84)+39] += A[267] * B[(l_n*84)+67];
    C[(l_n*84)+67] += A[268] * B[(l_n*84)+67];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b68 = _mm256_broadcast_sd(&B[(l_n*84)+68]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b68 = _mm_loaddup_pd(&B[(l_n*84)+68]);
#endif
    __m128d c68_0 = _mm_load_sd(&C[(l_n*84)+40]);
    __m128d a68_0 = _mm_load_sd(&A[269]);
#if defined(__SSE3__) && defined(__AVX__)
    c68_0 = _mm_add_sd(c68_0, _mm_mul_sd(a68_0, _mm256_castpd256_pd128(b68)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c68_0 = _mm_add_sd(c68_0, _mm_mul_sd(a68_0, b68));
#endif
    _mm_store_sd(&C[(l_n*84)+40], c68_0);
    __m128d c68_1 = _mm_load_sd(&C[(l_n*84)+68]);
    __m128d a68_1 = _mm_load_sd(&A[270]);
#if defined(__SSE3__) && defined(__AVX__)
    c68_1 = _mm_add_sd(c68_1, _mm_mul_sd(a68_1, _mm256_castpd256_pd128(b68)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c68_1 = _mm_add_sd(c68_1, _mm_mul_sd(a68_1, b68));
#endif
    _mm_store_sd(&C[(l_n*84)+68], c68_1);
#else
    C[(l_n*84)+40] += A[269] * B[(l_n*84)+68];
    C[(l_n*84)+68] += A[270] * B[(l_n*84)+68];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b69 = _mm256_broadcast_sd(&B[(l_n*84)+69]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b69 = _mm_loaddup_pd(&B[(l_n*84)+69]);
#endif
    __m128d c69_0 = _mm_load_sd(&C[(l_n*84)+20]);
    __m128d a69_0 = _mm_load_sd(&A[271]);
#if defined(__SSE3__) && defined(__AVX__)
    c69_0 = _mm_add_sd(c69_0, _mm_mul_sd(a69_0, _mm256_castpd256_pd128(b69)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c69_0 = _mm_add_sd(c69_0, _mm_mul_sd(a69_0, b69));
#endif
    _mm_store_sd(&C[(l_n*84)+20], c69_0);
    __m128d c69_1 = _mm_load_sd(&C[(l_n*84)+41]);
    __m128d a69_1 = _mm_load_sd(&A[272]);
#if defined(__SSE3__) && defined(__AVX__)
    c69_1 = _mm_add_sd(c69_1, _mm_mul_sd(a69_1, _mm256_castpd256_pd128(b69)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c69_1 = _mm_add_sd(c69_1, _mm_mul_sd(a69_1, b69));
#endif
    _mm_store_sd(&C[(l_n*84)+41], c69_1);
    __m128d c69_2 = _mm_load_sd(&C[(l_n*84)+69]);
    __m128d a69_2 = _mm_load_sd(&A[273]);
#if defined(__SSE3__) && defined(__AVX__)
    c69_2 = _mm_add_sd(c69_2, _mm_mul_sd(a69_2, _mm256_castpd256_pd128(b69)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c69_2 = _mm_add_sd(c69_2, _mm_mul_sd(a69_2, b69));
#endif
    _mm_store_sd(&C[(l_n*84)+69], c69_2);
#else
    C[(l_n*84)+20] += A[271] * B[(l_n*84)+69];
    C[(l_n*84)+41] += A[272] * B[(l_n*84)+69];
    C[(l_n*84)+69] += A[273] * B[(l_n*84)+69];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b70 = _mm256_broadcast_sd(&B[(l_n*84)+70]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b70 = _mm_loaddup_pd(&B[(l_n*84)+70]);
#endif
    __m128d c70_0 = _mm_load_sd(&C[(l_n*84)+21]);
    __m128d a70_0 = _mm_load_sd(&A[274]);
#if defined(__SSE3__) && defined(__AVX__)
    c70_0 = _mm_add_sd(c70_0, _mm_mul_sd(a70_0, _mm256_castpd256_pd128(b70)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c70_0 = _mm_add_sd(c70_0, _mm_mul_sd(a70_0, b70));
#endif
    _mm_store_sd(&C[(l_n*84)+21], c70_0);
    __m128d c70_1 = _mm_load_sd(&C[(l_n*84)+42]);
    __m128d a70_1 = _mm_load_sd(&A[275]);
#if defined(__SSE3__) && defined(__AVX__)
    c70_1 = _mm_add_sd(c70_1, _mm_mul_sd(a70_1, _mm256_castpd256_pd128(b70)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c70_1 = _mm_add_sd(c70_1, _mm_mul_sd(a70_1, b70));
#endif
    _mm_store_sd(&C[(l_n*84)+42], c70_1);
    __m128d c70_2 = _mm_load_sd(&C[(l_n*84)+70]);
    __m128d a70_2 = _mm_load_sd(&A[276]);
#if defined(__SSE3__) && defined(__AVX__)
    c70_2 = _mm_add_sd(c70_2, _mm_mul_sd(a70_2, _mm256_castpd256_pd128(b70)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c70_2 = _mm_add_sd(c70_2, _mm_mul_sd(a70_2, b70));
#endif
    _mm_store_sd(&C[(l_n*84)+70], c70_2);
#else
    C[(l_n*84)+21] += A[274] * B[(l_n*84)+70];
    C[(l_n*84)+42] += A[275] * B[(l_n*84)+70];
    C[(l_n*84)+70] += A[276] * B[(l_n*84)+70];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b71 = _mm256_broadcast_sd(&B[(l_n*84)+71]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b71 = _mm_loaddup_pd(&B[(l_n*84)+71]);
#endif
    __m128d c71_0 = _mm_load_sd(&C[(l_n*84)+22]);
    __m128d a71_0 = _mm_load_sd(&A[277]);
#if defined(__SSE3__) && defined(__AVX__)
    c71_0 = _mm_add_sd(c71_0, _mm_mul_sd(a71_0, _mm256_castpd256_pd128(b71)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c71_0 = _mm_add_sd(c71_0, _mm_mul_sd(a71_0, b71));
#endif
    _mm_store_sd(&C[(l_n*84)+22], c71_0);
    __m128d c71_1 = _mm_load_sd(&C[(l_n*84)+43]);
    __m128d a71_1 = _mm_load_sd(&A[278]);
#if defined(__SSE3__) && defined(__AVX__)
    c71_1 = _mm_add_sd(c71_1, _mm_mul_sd(a71_1, _mm256_castpd256_pd128(b71)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c71_1 = _mm_add_sd(c71_1, _mm_mul_sd(a71_1, b71));
#endif
    _mm_store_sd(&C[(l_n*84)+43], c71_1);
    __m128d c71_2 = _mm_load_sd(&C[(l_n*84)+71]);
    __m128d a71_2 = _mm_load_sd(&A[279]);
#if defined(__SSE3__) && defined(__AVX__)
    c71_2 = _mm_add_sd(c71_2, _mm_mul_sd(a71_2, _mm256_castpd256_pd128(b71)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c71_2 = _mm_add_sd(c71_2, _mm_mul_sd(a71_2, b71));
#endif
    _mm_store_sd(&C[(l_n*84)+71], c71_2);
#else
    C[(l_n*84)+22] += A[277] * B[(l_n*84)+71];
    C[(l_n*84)+43] += A[278] * B[(l_n*84)+71];
    C[(l_n*84)+71] += A[279] * B[(l_n*84)+71];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b72 = _mm256_broadcast_sd(&B[(l_n*84)+72]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b72 = _mm_loaddup_pd(&B[(l_n*84)+72]);
#endif
    __m128d c72_0 = _mm_load_sd(&C[(l_n*84)+23]);
    __m128d a72_0 = _mm_load_sd(&A[280]);
#if defined(__SSE3__) && defined(__AVX__)
    c72_0 = _mm_add_sd(c72_0, _mm_mul_sd(a72_0, _mm256_castpd256_pd128(b72)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c72_0 = _mm_add_sd(c72_0, _mm_mul_sd(a72_0, b72));
#endif
    _mm_store_sd(&C[(l_n*84)+23], c72_0);
    __m128d c72_1 = _mm_load_sd(&C[(l_n*84)+44]);
    __m128d a72_1 = _mm_load_sd(&A[281]);
#if defined(__SSE3__) && defined(__AVX__)
    c72_1 = _mm_add_sd(c72_1, _mm_mul_sd(a72_1, _mm256_castpd256_pd128(b72)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c72_1 = _mm_add_sd(c72_1, _mm_mul_sd(a72_1, b72));
#endif
    _mm_store_sd(&C[(l_n*84)+44], c72_1);
    __m128d c72_2 = _mm_load_sd(&C[(l_n*84)+72]);
    __m128d a72_2 = _mm_load_sd(&A[282]);
#if defined(__SSE3__) && defined(__AVX__)
    c72_2 = _mm_add_sd(c72_2, _mm_mul_sd(a72_2, _mm256_castpd256_pd128(b72)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c72_2 = _mm_add_sd(c72_2, _mm_mul_sd(a72_2, b72));
#endif
    _mm_store_sd(&C[(l_n*84)+72], c72_2);
#else
    C[(l_n*84)+23] += A[280] * B[(l_n*84)+72];
    C[(l_n*84)+44] += A[281] * B[(l_n*84)+72];
    C[(l_n*84)+72] += A[282] * B[(l_n*84)+72];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b73 = _mm256_broadcast_sd(&B[(l_n*84)+73]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b73 = _mm_loaddup_pd(&B[(l_n*84)+73]);
#endif
    __m128d c73_0 = _mm_load_sd(&C[(l_n*84)+24]);
    __m128d a73_0 = _mm_load_sd(&A[283]);
#if defined(__SSE3__) && defined(__AVX__)
    c73_0 = _mm_add_sd(c73_0, _mm_mul_sd(a73_0, _mm256_castpd256_pd128(b73)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c73_0 = _mm_add_sd(c73_0, _mm_mul_sd(a73_0, b73));
#endif
    _mm_store_sd(&C[(l_n*84)+24], c73_0);
    __m128d c73_1 = _mm_load_sd(&C[(l_n*84)+45]);
    __m128d a73_1 = _mm_load_sd(&A[284]);
#if defined(__SSE3__) && defined(__AVX__)
    c73_1 = _mm_add_sd(c73_1, _mm_mul_sd(a73_1, _mm256_castpd256_pd128(b73)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c73_1 = _mm_add_sd(c73_1, _mm_mul_sd(a73_1, b73));
#endif
    _mm_store_sd(&C[(l_n*84)+45], c73_1);
    __m128d c73_2 = _mm_load_sd(&C[(l_n*84)+73]);
    __m128d a73_2 = _mm_load_sd(&A[285]);
#if defined(__SSE3__) && defined(__AVX__)
    c73_2 = _mm_add_sd(c73_2, _mm_mul_sd(a73_2, _mm256_castpd256_pd128(b73)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c73_2 = _mm_add_sd(c73_2, _mm_mul_sd(a73_2, b73));
#endif
    _mm_store_sd(&C[(l_n*84)+73], c73_2);
#else
    C[(l_n*84)+24] += A[283] * B[(l_n*84)+73];
    C[(l_n*84)+45] += A[284] * B[(l_n*84)+73];
    C[(l_n*84)+73] += A[285] * B[(l_n*84)+73];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b74 = _mm256_broadcast_sd(&B[(l_n*84)+74]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b74 = _mm_loaddup_pd(&B[(l_n*84)+74]);
#endif
    __m128d c74_0 = _mm_load_sd(&C[(l_n*84)+10]);
    __m128d a74_0 = _mm_load_sd(&A[286]);
#if defined(__SSE3__) && defined(__AVX__)
    c74_0 = _mm_add_sd(c74_0, _mm_mul_sd(a74_0, _mm256_castpd256_pd128(b74)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c74_0 = _mm_add_sd(c74_0, _mm_mul_sd(a74_0, b74));
#endif
    _mm_store_sd(&C[(l_n*84)+10], c74_0);
    __m128d c74_1 = _mm_load_sd(&C[(l_n*84)+25]);
    __m128d a74_1 = _mm_load_sd(&A[287]);
#if defined(__SSE3__) && defined(__AVX__)
    c74_1 = _mm_add_sd(c74_1, _mm_mul_sd(a74_1, _mm256_castpd256_pd128(b74)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c74_1 = _mm_add_sd(c74_1, _mm_mul_sd(a74_1, b74));
#endif
    _mm_store_sd(&C[(l_n*84)+25], c74_1);
    __m128d c74_2 = _mm_load_sd(&C[(l_n*84)+46]);
    __m128d a74_2 = _mm_load_sd(&A[288]);
#if defined(__SSE3__) && defined(__AVX__)
    c74_2 = _mm_add_sd(c74_2, _mm_mul_sd(a74_2, _mm256_castpd256_pd128(b74)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c74_2 = _mm_add_sd(c74_2, _mm_mul_sd(a74_2, b74));
#endif
    _mm_store_sd(&C[(l_n*84)+46], c74_2);
    __m128d c74_3 = _mm_load_sd(&C[(l_n*84)+74]);
    __m128d a74_3 = _mm_load_sd(&A[289]);
#if defined(__SSE3__) && defined(__AVX__)
    c74_3 = _mm_add_sd(c74_3, _mm_mul_sd(a74_3, _mm256_castpd256_pd128(b74)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c74_3 = _mm_add_sd(c74_3, _mm_mul_sd(a74_3, b74));
#endif
    _mm_store_sd(&C[(l_n*84)+74], c74_3);
#else
    C[(l_n*84)+10] += A[286] * B[(l_n*84)+74];
    C[(l_n*84)+25] += A[287] * B[(l_n*84)+74];
    C[(l_n*84)+46] += A[288] * B[(l_n*84)+74];
    C[(l_n*84)+74] += A[289] * B[(l_n*84)+74];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b75 = _mm256_broadcast_sd(&B[(l_n*84)+75]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b75 = _mm_loaddup_pd(&B[(l_n*84)+75]);
#endif
    __m128d c75_0 = _mm_load_sd(&C[(l_n*84)+11]);
    __m128d a75_0 = _mm_load_sd(&A[290]);
#if defined(__SSE3__) && defined(__AVX__)
    c75_0 = _mm_add_sd(c75_0, _mm_mul_sd(a75_0, _mm256_castpd256_pd128(b75)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c75_0 = _mm_add_sd(c75_0, _mm_mul_sd(a75_0, b75));
#endif
    _mm_store_sd(&C[(l_n*84)+11], c75_0);
    __m128d c75_1 = _mm_load_sd(&C[(l_n*84)+26]);
    __m128d a75_1 = _mm_load_sd(&A[291]);
#if defined(__SSE3__) && defined(__AVX__)
    c75_1 = _mm_add_sd(c75_1, _mm_mul_sd(a75_1, _mm256_castpd256_pd128(b75)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c75_1 = _mm_add_sd(c75_1, _mm_mul_sd(a75_1, b75));
#endif
    _mm_store_sd(&C[(l_n*84)+26], c75_1);
    __m128d c75_2 = _mm_load_sd(&C[(l_n*84)+47]);
    __m128d a75_2 = _mm_load_sd(&A[292]);
#if defined(__SSE3__) && defined(__AVX__)
    c75_2 = _mm_add_sd(c75_2, _mm_mul_sd(a75_2, _mm256_castpd256_pd128(b75)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c75_2 = _mm_add_sd(c75_2, _mm_mul_sd(a75_2, b75));
#endif
    _mm_store_sd(&C[(l_n*84)+47], c75_2);
    __m128d c75_3 = _mm_load_sd(&C[(l_n*84)+75]);
    __m128d a75_3 = _mm_load_sd(&A[293]);
#if defined(__SSE3__) && defined(__AVX__)
    c75_3 = _mm_add_sd(c75_3, _mm_mul_sd(a75_3, _mm256_castpd256_pd128(b75)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c75_3 = _mm_add_sd(c75_3, _mm_mul_sd(a75_3, b75));
#endif
    _mm_store_sd(&C[(l_n*84)+75], c75_3);
#else
    C[(l_n*84)+11] += A[290] * B[(l_n*84)+75];
    C[(l_n*84)+26] += A[291] * B[(l_n*84)+75];
    C[(l_n*84)+47] += A[292] * B[(l_n*84)+75];
    C[(l_n*84)+75] += A[293] * B[(l_n*84)+75];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b76 = _mm256_broadcast_sd(&B[(l_n*84)+76]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b76 = _mm_loaddup_pd(&B[(l_n*84)+76]);
#endif
    __m128d c76_0 = _mm_load_sd(&C[(l_n*84)+12]);
    __m128d a76_0 = _mm_load_sd(&A[294]);
#if defined(__SSE3__) && defined(__AVX__)
    c76_0 = _mm_add_sd(c76_0, _mm_mul_sd(a76_0, _mm256_castpd256_pd128(b76)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c76_0 = _mm_add_sd(c76_0, _mm_mul_sd(a76_0, b76));
#endif
    _mm_store_sd(&C[(l_n*84)+12], c76_0);
    __m128d c76_1 = _mm_load_sd(&C[(l_n*84)+27]);
    __m128d a76_1 = _mm_load_sd(&A[295]);
#if defined(__SSE3__) && defined(__AVX__)
    c76_1 = _mm_add_sd(c76_1, _mm_mul_sd(a76_1, _mm256_castpd256_pd128(b76)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c76_1 = _mm_add_sd(c76_1, _mm_mul_sd(a76_1, b76));
#endif
    _mm_store_sd(&C[(l_n*84)+27], c76_1);
    __m128d c76_2 = _mm_load_sd(&C[(l_n*84)+48]);
    __m128d a76_2 = _mm_load_sd(&A[296]);
#if defined(__SSE3__) && defined(__AVX__)
    c76_2 = _mm_add_sd(c76_2, _mm_mul_sd(a76_2, _mm256_castpd256_pd128(b76)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c76_2 = _mm_add_sd(c76_2, _mm_mul_sd(a76_2, b76));
#endif
    _mm_store_sd(&C[(l_n*84)+48], c76_2);
    __m128d c76_3 = _mm_load_sd(&C[(l_n*84)+76]);
    __m128d a76_3 = _mm_load_sd(&A[297]);
#if defined(__SSE3__) && defined(__AVX__)
    c76_3 = _mm_add_sd(c76_3, _mm_mul_sd(a76_3, _mm256_castpd256_pd128(b76)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c76_3 = _mm_add_sd(c76_3, _mm_mul_sd(a76_3, b76));
#endif
    _mm_store_sd(&C[(l_n*84)+76], c76_3);
#else
    C[(l_n*84)+12] += A[294] * B[(l_n*84)+76];
    C[(l_n*84)+27] += A[295] * B[(l_n*84)+76];
    C[(l_n*84)+48] += A[296] * B[(l_n*84)+76];
    C[(l_n*84)+76] += A[297] * B[(l_n*84)+76];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b77 = _mm256_broadcast_sd(&B[(l_n*84)+77]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b77 = _mm_loaddup_pd(&B[(l_n*84)+77]);
#endif
    __m128d c77_0 = _mm_load_sd(&C[(l_n*84)+13]);
    __m128d a77_0 = _mm_load_sd(&A[298]);
#if defined(__SSE3__) && defined(__AVX__)
    c77_0 = _mm_add_sd(c77_0, _mm_mul_sd(a77_0, _mm256_castpd256_pd128(b77)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c77_0 = _mm_add_sd(c77_0, _mm_mul_sd(a77_0, b77));
#endif
    _mm_store_sd(&C[(l_n*84)+13], c77_0);
    __m128d c77_1 = _mm_load_sd(&C[(l_n*84)+28]);
    __m128d a77_1 = _mm_load_sd(&A[299]);
#if defined(__SSE3__) && defined(__AVX__)
    c77_1 = _mm_add_sd(c77_1, _mm_mul_sd(a77_1, _mm256_castpd256_pd128(b77)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c77_1 = _mm_add_sd(c77_1, _mm_mul_sd(a77_1, b77));
#endif
    _mm_store_sd(&C[(l_n*84)+28], c77_1);
    __m128d c77_2 = _mm_load_sd(&C[(l_n*84)+49]);
    __m128d a77_2 = _mm_load_sd(&A[300]);
#if defined(__SSE3__) && defined(__AVX__)
    c77_2 = _mm_add_sd(c77_2, _mm_mul_sd(a77_2, _mm256_castpd256_pd128(b77)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c77_2 = _mm_add_sd(c77_2, _mm_mul_sd(a77_2, b77));
#endif
    _mm_store_sd(&C[(l_n*84)+49], c77_2);
    __m128d c77_3 = _mm_load_sd(&C[(l_n*84)+77]);
    __m128d a77_3 = _mm_load_sd(&A[301]);
#if defined(__SSE3__) && defined(__AVX__)
    c77_3 = _mm_add_sd(c77_3, _mm_mul_sd(a77_3, _mm256_castpd256_pd128(b77)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c77_3 = _mm_add_sd(c77_3, _mm_mul_sd(a77_3, b77));
#endif
    _mm_store_sd(&C[(l_n*84)+77], c77_3);
#else
    C[(l_n*84)+13] += A[298] * B[(l_n*84)+77];
    C[(l_n*84)+28] += A[299] * B[(l_n*84)+77];
    C[(l_n*84)+49] += A[300] * B[(l_n*84)+77];
    C[(l_n*84)+77] += A[301] * B[(l_n*84)+77];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b78 = _mm256_broadcast_sd(&B[(l_n*84)+78]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b78 = _mm_loaddup_pd(&B[(l_n*84)+78]);
#endif
    __m128d c78_0 = _mm_load_sd(&C[(l_n*84)+4]);
    __m128d a78_0 = _mm_load_sd(&A[302]);
#if defined(__SSE3__) && defined(__AVX__)
    c78_0 = _mm_add_sd(c78_0, _mm_mul_sd(a78_0, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c78_0 = _mm_add_sd(c78_0, _mm_mul_sd(a78_0, b78));
#endif
    _mm_store_sd(&C[(l_n*84)+4], c78_0);
    __m128d c78_1 = _mm_load_sd(&C[(l_n*84)+14]);
    __m128d a78_1 = _mm_load_sd(&A[303]);
#if defined(__SSE3__) && defined(__AVX__)
    c78_1 = _mm_add_sd(c78_1, _mm_mul_sd(a78_1, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c78_1 = _mm_add_sd(c78_1, _mm_mul_sd(a78_1, b78));
#endif
    _mm_store_sd(&C[(l_n*84)+14], c78_1);
    __m128d c78_2 = _mm_load_sd(&C[(l_n*84)+29]);
    __m128d a78_2 = _mm_load_sd(&A[304]);
#if defined(__SSE3__) && defined(__AVX__)
    c78_2 = _mm_add_sd(c78_2, _mm_mul_sd(a78_2, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c78_2 = _mm_add_sd(c78_2, _mm_mul_sd(a78_2, b78));
#endif
    _mm_store_sd(&C[(l_n*84)+29], c78_2);
    __m128d c78_3 = _mm_load_sd(&C[(l_n*84)+50]);
    __m128d a78_3 = _mm_load_sd(&A[305]);
#if defined(__SSE3__) && defined(__AVX__)
    c78_3 = _mm_add_sd(c78_3, _mm_mul_sd(a78_3, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c78_3 = _mm_add_sd(c78_3, _mm_mul_sd(a78_3, b78));
#endif
    _mm_store_sd(&C[(l_n*84)+50], c78_3);
    __m128d c78_4 = _mm_load_sd(&C[(l_n*84)+78]);
    __m128d a78_4 = _mm_load_sd(&A[306]);
#if defined(__SSE3__) && defined(__AVX__)
    c78_4 = _mm_add_sd(c78_4, _mm_mul_sd(a78_4, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c78_4 = _mm_add_sd(c78_4, _mm_mul_sd(a78_4, b78));
#endif
    _mm_store_sd(&C[(l_n*84)+78], c78_4);
#else
    C[(l_n*84)+4] += A[302] * B[(l_n*84)+78];
    C[(l_n*84)+14] += A[303] * B[(l_n*84)+78];
    C[(l_n*84)+29] += A[304] * B[(l_n*84)+78];
    C[(l_n*84)+50] += A[305] * B[(l_n*84)+78];
    C[(l_n*84)+78] += A[306] * B[(l_n*84)+78];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b79 = _mm256_broadcast_sd(&B[(l_n*84)+79]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b79 = _mm_loaddup_pd(&B[(l_n*84)+79]);
#endif
    __m128d c79_0 = _mm_load_sd(&C[(l_n*84)+5]);
    __m128d a79_0 = _mm_load_sd(&A[307]);
#if defined(__SSE3__) && defined(__AVX__)
    c79_0 = _mm_add_sd(c79_0, _mm_mul_sd(a79_0, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c79_0 = _mm_add_sd(c79_0, _mm_mul_sd(a79_0, b79));
#endif
    _mm_store_sd(&C[(l_n*84)+5], c79_0);
    __m128d c79_1 = _mm_load_sd(&C[(l_n*84)+15]);
    __m128d a79_1 = _mm_load_sd(&A[308]);
#if defined(__SSE3__) && defined(__AVX__)
    c79_1 = _mm_add_sd(c79_1, _mm_mul_sd(a79_1, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c79_1 = _mm_add_sd(c79_1, _mm_mul_sd(a79_1, b79));
#endif
    _mm_store_sd(&C[(l_n*84)+15], c79_1);
    __m128d c79_2 = _mm_load_sd(&C[(l_n*84)+30]);
    __m128d a79_2 = _mm_load_sd(&A[309]);
#if defined(__SSE3__) && defined(__AVX__)
    c79_2 = _mm_add_sd(c79_2, _mm_mul_sd(a79_2, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c79_2 = _mm_add_sd(c79_2, _mm_mul_sd(a79_2, b79));
#endif
    _mm_store_sd(&C[(l_n*84)+30], c79_2);
    __m128d c79_3 = _mm_load_sd(&C[(l_n*84)+51]);
    __m128d a79_3 = _mm_load_sd(&A[310]);
#if defined(__SSE3__) && defined(__AVX__)
    c79_3 = _mm_add_sd(c79_3, _mm_mul_sd(a79_3, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c79_3 = _mm_add_sd(c79_3, _mm_mul_sd(a79_3, b79));
#endif
    _mm_store_sd(&C[(l_n*84)+51], c79_3);
    __m128d c79_4 = _mm_load_sd(&C[(l_n*84)+79]);
    __m128d a79_4 = _mm_load_sd(&A[311]);
#if defined(__SSE3__) && defined(__AVX__)
    c79_4 = _mm_add_sd(c79_4, _mm_mul_sd(a79_4, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c79_4 = _mm_add_sd(c79_4, _mm_mul_sd(a79_4, b79));
#endif
    _mm_store_sd(&C[(l_n*84)+79], c79_4);
#else
    C[(l_n*84)+5] += A[307] * B[(l_n*84)+79];
    C[(l_n*84)+15] += A[308] * B[(l_n*84)+79];
    C[(l_n*84)+30] += A[309] * B[(l_n*84)+79];
    C[(l_n*84)+51] += A[310] * B[(l_n*84)+79];
    C[(l_n*84)+79] += A[311] * B[(l_n*84)+79];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b80 = _mm256_broadcast_sd(&B[(l_n*84)+80]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b80 = _mm_loaddup_pd(&B[(l_n*84)+80]);
#endif
    __m128d c80_0 = _mm_load_sd(&C[(l_n*84)+6]);
    __m128d a80_0 = _mm_load_sd(&A[312]);
#if defined(__SSE3__) && defined(__AVX__)
    c80_0 = _mm_add_sd(c80_0, _mm_mul_sd(a80_0, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c80_0 = _mm_add_sd(c80_0, _mm_mul_sd(a80_0, b80));
#endif
    _mm_store_sd(&C[(l_n*84)+6], c80_0);
    __m128d c80_1 = _mm_load_sd(&C[(l_n*84)+16]);
    __m128d a80_1 = _mm_load_sd(&A[313]);
#if defined(__SSE3__) && defined(__AVX__)
    c80_1 = _mm_add_sd(c80_1, _mm_mul_sd(a80_1, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c80_1 = _mm_add_sd(c80_1, _mm_mul_sd(a80_1, b80));
#endif
    _mm_store_sd(&C[(l_n*84)+16], c80_1);
    __m128d c80_2 = _mm_load_sd(&C[(l_n*84)+31]);
    __m128d a80_2 = _mm_load_sd(&A[314]);
#if defined(__SSE3__) && defined(__AVX__)
    c80_2 = _mm_add_sd(c80_2, _mm_mul_sd(a80_2, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c80_2 = _mm_add_sd(c80_2, _mm_mul_sd(a80_2, b80));
#endif
    _mm_store_sd(&C[(l_n*84)+31], c80_2);
    __m128d c80_3 = _mm_load_sd(&C[(l_n*84)+52]);
    __m128d a80_3 = _mm_load_sd(&A[315]);
#if defined(__SSE3__) && defined(__AVX__)
    c80_3 = _mm_add_sd(c80_3, _mm_mul_sd(a80_3, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c80_3 = _mm_add_sd(c80_3, _mm_mul_sd(a80_3, b80));
#endif
    _mm_store_sd(&C[(l_n*84)+52], c80_3);
    __m128d c80_4 = _mm_load_sd(&C[(l_n*84)+80]);
    __m128d a80_4 = _mm_load_sd(&A[316]);
#if defined(__SSE3__) && defined(__AVX__)
    c80_4 = _mm_add_sd(c80_4, _mm_mul_sd(a80_4, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c80_4 = _mm_add_sd(c80_4, _mm_mul_sd(a80_4, b80));
#endif
    _mm_store_sd(&C[(l_n*84)+80], c80_4);
#else
    C[(l_n*84)+6] += A[312] * B[(l_n*84)+80];
    C[(l_n*84)+16] += A[313] * B[(l_n*84)+80];
    C[(l_n*84)+31] += A[314] * B[(l_n*84)+80];
    C[(l_n*84)+52] += A[315] * B[(l_n*84)+80];
    C[(l_n*84)+80] += A[316] * B[(l_n*84)+80];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b81 = _mm256_broadcast_sd(&B[(l_n*84)+81]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b81 = _mm_loaddup_pd(&B[(l_n*84)+81]);
#endif
    __m128d c81_0 = _mm_load_sd(&C[(l_n*84)+1]);
    __m128d a81_0 = _mm_load_sd(&A[317]);
#if defined(__SSE3__) && defined(__AVX__)
    c81_0 = _mm_add_sd(c81_0, _mm_mul_sd(a81_0, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c81_0 = _mm_add_sd(c81_0, _mm_mul_sd(a81_0, b81));
#endif
    _mm_store_sd(&C[(l_n*84)+1], c81_0);
    __m128d c81_1 = _mm_load_sd(&C[(l_n*84)+7]);
    __m128d a81_1 = _mm_load_sd(&A[318]);
#if defined(__SSE3__) && defined(__AVX__)
    c81_1 = _mm_add_sd(c81_1, _mm_mul_sd(a81_1, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c81_1 = _mm_add_sd(c81_1, _mm_mul_sd(a81_1, b81));
#endif
    _mm_store_sd(&C[(l_n*84)+7], c81_1);
    __m128d c81_2 = _mm_load_sd(&C[(l_n*84)+17]);
    __m128d a81_2 = _mm_load_sd(&A[319]);
#if defined(__SSE3__) && defined(__AVX__)
    c81_2 = _mm_add_sd(c81_2, _mm_mul_sd(a81_2, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c81_2 = _mm_add_sd(c81_2, _mm_mul_sd(a81_2, b81));
#endif
    _mm_store_sd(&C[(l_n*84)+17], c81_2);
    __m128d c81_3 = _mm_load_sd(&C[(l_n*84)+32]);
    __m128d a81_3 = _mm_load_sd(&A[320]);
#if defined(__SSE3__) && defined(__AVX__)
    c81_3 = _mm_add_sd(c81_3, _mm_mul_sd(a81_3, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c81_3 = _mm_add_sd(c81_3, _mm_mul_sd(a81_3, b81));
#endif
    _mm_store_sd(&C[(l_n*84)+32], c81_3);
    __m128d c81_4 = _mm_load_sd(&C[(l_n*84)+53]);
    __m128d a81_4 = _mm_load_sd(&A[321]);
#if defined(__SSE3__) && defined(__AVX__)
    c81_4 = _mm_add_sd(c81_4, _mm_mul_sd(a81_4, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c81_4 = _mm_add_sd(c81_4, _mm_mul_sd(a81_4, b81));
#endif
    _mm_store_sd(&C[(l_n*84)+53], c81_4);
    __m128d c81_5 = _mm_load_sd(&C[(l_n*84)+81]);
    __m128d a81_5 = _mm_load_sd(&A[322]);
#if defined(__SSE3__) && defined(__AVX__)
    c81_5 = _mm_add_sd(c81_5, _mm_mul_sd(a81_5, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c81_5 = _mm_add_sd(c81_5, _mm_mul_sd(a81_5, b81));
#endif
    _mm_store_sd(&C[(l_n*84)+81], c81_5);
#else
    C[(l_n*84)+1] += A[317] * B[(l_n*84)+81];
    C[(l_n*84)+7] += A[318] * B[(l_n*84)+81];
    C[(l_n*84)+17] += A[319] * B[(l_n*84)+81];
    C[(l_n*84)+32] += A[320] * B[(l_n*84)+81];
    C[(l_n*84)+53] += A[321] * B[(l_n*84)+81];
    C[(l_n*84)+81] += A[322] * B[(l_n*84)+81];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b82 = _mm256_broadcast_sd(&B[(l_n*84)+82]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b82 = _mm_loaddup_pd(&B[(l_n*84)+82]);
#endif
    __m128d c82_0 = _mm_load_sd(&C[(l_n*84)+2]);
    __m128d a82_0 = _mm_load_sd(&A[323]);
#if defined(__SSE3__) && defined(__AVX__)
    c82_0 = _mm_add_sd(c82_0, _mm_mul_sd(a82_0, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c82_0 = _mm_add_sd(c82_0, _mm_mul_sd(a82_0, b82));
#endif
    _mm_store_sd(&C[(l_n*84)+2], c82_0);
    __m128d c82_1 = _mm_load_sd(&C[(l_n*84)+8]);
    __m128d a82_1 = _mm_load_sd(&A[324]);
#if defined(__SSE3__) && defined(__AVX__)
    c82_1 = _mm_add_sd(c82_1, _mm_mul_sd(a82_1, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c82_1 = _mm_add_sd(c82_1, _mm_mul_sd(a82_1, b82));
#endif
    _mm_store_sd(&C[(l_n*84)+8], c82_1);
    __m128d c82_2 = _mm_load_sd(&C[(l_n*84)+18]);
    __m128d a82_2 = _mm_load_sd(&A[325]);
#if defined(__SSE3__) && defined(__AVX__)
    c82_2 = _mm_add_sd(c82_2, _mm_mul_sd(a82_2, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c82_2 = _mm_add_sd(c82_2, _mm_mul_sd(a82_2, b82));
#endif
    _mm_store_sd(&C[(l_n*84)+18], c82_2);
    __m128d c82_3 = _mm_load_sd(&C[(l_n*84)+33]);
    __m128d a82_3 = _mm_load_sd(&A[326]);
#if defined(__SSE3__) && defined(__AVX__)
    c82_3 = _mm_add_sd(c82_3, _mm_mul_sd(a82_3, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c82_3 = _mm_add_sd(c82_3, _mm_mul_sd(a82_3, b82));
#endif
    _mm_store_sd(&C[(l_n*84)+33], c82_3);
    __m128d c82_4 = _mm_load_sd(&C[(l_n*84)+54]);
    __m128d a82_4 = _mm_load_sd(&A[327]);
#if defined(__SSE3__) && defined(__AVX__)
    c82_4 = _mm_add_sd(c82_4, _mm_mul_sd(a82_4, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c82_4 = _mm_add_sd(c82_4, _mm_mul_sd(a82_4, b82));
#endif
    _mm_store_sd(&C[(l_n*84)+54], c82_4);
    __m128d c82_5 = _mm_load_sd(&C[(l_n*84)+82]);
    __m128d a82_5 = _mm_load_sd(&A[328]);
#if defined(__SSE3__) && defined(__AVX__)
    c82_5 = _mm_add_sd(c82_5, _mm_mul_sd(a82_5, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c82_5 = _mm_add_sd(c82_5, _mm_mul_sd(a82_5, b82));
#endif
    _mm_store_sd(&C[(l_n*84)+82], c82_5);
#else
    C[(l_n*84)+2] += A[323] * B[(l_n*84)+82];
    C[(l_n*84)+8] += A[324] * B[(l_n*84)+82];
    C[(l_n*84)+18] += A[325] * B[(l_n*84)+82];
    C[(l_n*84)+33] += A[326] * B[(l_n*84)+82];
    C[(l_n*84)+54] += A[327] * B[(l_n*84)+82];
    C[(l_n*84)+82] += A[328] * B[(l_n*84)+82];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b83 = _mm256_broadcast_sd(&B[(l_n*84)+83]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b83 = _mm_loaddup_pd(&B[(l_n*84)+83]);
#endif
    __m128d c83_0 = _mm_load_sd(&C[(l_n*84)+0]);
    __m128d a83_0 = _mm_load_sd(&A[329]);
#if defined(__SSE3__) && defined(__AVX__)
    c83_0 = _mm_add_sd(c83_0, _mm_mul_sd(a83_0, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c83_0 = _mm_add_sd(c83_0, _mm_mul_sd(a83_0, b83));
#endif
    _mm_store_sd(&C[(l_n*84)+0], c83_0);
    __m128d c83_1 = _mm_load_sd(&C[(l_n*84)+3]);
    __m128d a83_1 = _mm_load_sd(&A[330]);
#if defined(__SSE3__) && defined(__AVX__)
    c83_1 = _mm_add_sd(c83_1, _mm_mul_sd(a83_1, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c83_1 = _mm_add_sd(c83_1, _mm_mul_sd(a83_1, b83));
#endif
    _mm_store_sd(&C[(l_n*84)+3], c83_1);
    __m128d c83_2 = _mm_load_sd(&C[(l_n*84)+9]);
    __m128d a83_2 = _mm_load_sd(&A[331]);
#if defined(__SSE3__) && defined(__AVX__)
    c83_2 = _mm_add_sd(c83_2, _mm_mul_sd(a83_2, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c83_2 = _mm_add_sd(c83_2, _mm_mul_sd(a83_2, b83));
#endif
    _mm_store_sd(&C[(l_n*84)+9], c83_2);
    __m128d c83_3 = _mm_load_sd(&C[(l_n*84)+19]);
    __m128d a83_3 = _mm_load_sd(&A[332]);
#if defined(__SSE3__) && defined(__AVX__)
    c83_3 = _mm_add_sd(c83_3, _mm_mul_sd(a83_3, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c83_3 = _mm_add_sd(c83_3, _mm_mul_sd(a83_3, b83));
#endif
    _mm_store_sd(&C[(l_n*84)+19], c83_3);
    __m128d c83_4 = _mm_load_sd(&C[(l_n*84)+34]);
    __m128d a83_4 = _mm_load_sd(&A[333]);
#if defined(__SSE3__) && defined(__AVX__)
    c83_4 = _mm_add_sd(c83_4, _mm_mul_sd(a83_4, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c83_4 = _mm_add_sd(c83_4, _mm_mul_sd(a83_4, b83));
#endif
    _mm_store_sd(&C[(l_n*84)+34], c83_4);
    __m128d c83_5 = _mm_load_sd(&C[(l_n*84)+55]);
    __m128d a83_5 = _mm_load_sd(&A[334]);
#if defined(__SSE3__) && defined(__AVX__)
    c83_5 = _mm_add_sd(c83_5, _mm_mul_sd(a83_5, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c83_5 = _mm_add_sd(c83_5, _mm_mul_sd(a83_5, b83));
#endif
    _mm_store_sd(&C[(l_n*84)+55], c83_5);
    __m128d c83_6 = _mm_load_sd(&C[(l_n*84)+83]);
    __m128d a83_6 = _mm_load_sd(&A[335]);
#if defined(__SSE3__) && defined(__AVX__)
    c83_6 = _mm_add_sd(c83_6, _mm_mul_sd(a83_6, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c83_6 = _mm_add_sd(c83_6, _mm_mul_sd(a83_6, b83));
#endif
    _mm_store_sd(&C[(l_n*84)+83], c83_6);
#else
    C[(l_n*84)+0] += A[329] * B[(l_n*84)+83];
    C[(l_n*84)+3] += A[330] * B[(l_n*84)+83];
    C[(l_n*84)+9] += A[331] * B[(l_n*84)+83];
    C[(l_n*84)+19] += A[332] * B[(l_n*84)+83];
    C[(l_n*84)+34] += A[333] * B[(l_n*84)+83];
    C[(l_n*84)+55] += A[334] * B[(l_n*84)+83];
    C[(l_n*84)+83] += A[335] * B[(l_n*84)+83];
#endif

  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 6048;
#endif
}

void dsparse_fP113DivM_m84_n9_k84_ldAna7_ldB84_ldC84_beta0_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_n = 0;
  #pragma nounroll_and_jam
  for ( l_n = 0; l_n < 9; l_n++) {
    unsigned int l_m = 0;
   #pragma simd
    for ( l_m = 0; l_m < 84; l_m++) {
      C[(l_n*84)+l_m] = 0.0;
    }
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b0 = _mm256_broadcast_sd(&B[(l_n*84)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b0 = _mm_loaddup_pd(&B[(l_n*84)+0]);
#endif
    __m128d c0_0 = _mm_load_sd(&C[(l_n*84)+0]);
    __m128d a0_0 = _mm_load_sd(&A[0]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, b0));
#endif
    _mm_store_sd(&C[(l_n*84)+0], c0_0);
    __m128d c0_1 = _mm_load_sd(&C[(l_n*84)+3]);
    __m128d a0_1 = _mm_load_sd(&A[1]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, b0));
#endif
    _mm_store_sd(&C[(l_n*84)+3], c0_1);
    __m128d c0_2 = _mm_load_sd(&C[(l_n*84)+9]);
    __m128d a0_2 = _mm_load_sd(&A[2]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, b0));
#endif
    _mm_store_sd(&C[(l_n*84)+9], c0_2);
    __m128d c0_3 = _mm_load_sd(&C[(l_n*84)+19]);
    __m128d a0_3 = _mm_load_sd(&A[3]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_3 = _mm_add_sd(c0_3, _mm_mul_sd(a0_3, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_3 = _mm_add_sd(c0_3, _mm_mul_sd(a0_3, b0));
#endif
    _mm_store_sd(&C[(l_n*84)+19], c0_3);
    __m128d c0_4 = _mm_load_sd(&C[(l_n*84)+34]);
    __m128d a0_4 = _mm_load_sd(&A[4]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_4 = _mm_add_sd(c0_4, _mm_mul_sd(a0_4, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_4 = _mm_add_sd(c0_4, _mm_mul_sd(a0_4, b0));
#endif
    _mm_store_sd(&C[(l_n*84)+34], c0_4);
    __m128d c0_5 = _mm_load_sd(&C[(l_n*84)+55]);
    __m128d a0_5 = _mm_load_sd(&A[5]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_5 = _mm_add_sd(c0_5, _mm_mul_sd(a0_5, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_5 = _mm_add_sd(c0_5, _mm_mul_sd(a0_5, b0));
#endif
    _mm_store_sd(&C[(l_n*84)+55], c0_5);
    __m128d c0_6 = _mm_load_sd(&C[(l_n*84)+83]);
    __m128d a0_6 = _mm_load_sd(&A[6]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_6 = _mm_add_sd(c0_6, _mm_mul_sd(a0_6, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_6 = _mm_add_sd(c0_6, _mm_mul_sd(a0_6, b0));
#endif
    _mm_store_sd(&C[(l_n*84)+83], c0_6);
#else
    C[(l_n*84)+0] += A[0] * B[(l_n*84)+0];
    C[(l_n*84)+3] += A[1] * B[(l_n*84)+0];
    C[(l_n*84)+9] += A[2] * B[(l_n*84)+0];
    C[(l_n*84)+19] += A[3] * B[(l_n*84)+0];
    C[(l_n*84)+34] += A[4] * B[(l_n*84)+0];
    C[(l_n*84)+55] += A[5] * B[(l_n*84)+0];
    C[(l_n*84)+83] += A[6] * B[(l_n*84)+0];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b1 = _mm256_broadcast_sd(&B[(l_n*84)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b1 = _mm_loaddup_pd(&B[(l_n*84)+1]);
#endif
    __m128d c1_0 = _mm_load_sd(&C[(l_n*84)+1]);
    __m128d a1_0 = _mm_load_sd(&A[7]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, b1));
#endif
    _mm_store_sd(&C[(l_n*84)+1], c1_0);
    __m128d c1_1 = _mm_load_sd(&C[(l_n*84)+7]);
    __m128d a1_1 = _mm_load_sd(&A[8]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_1 = _mm_add_sd(c1_1, _mm_mul_sd(a1_1, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_1 = _mm_add_sd(c1_1, _mm_mul_sd(a1_1, b1));
#endif
    _mm_store_sd(&C[(l_n*84)+7], c1_1);
    __m128d c1_2 = _mm_load_sd(&C[(l_n*84)+17]);
    __m128d a1_2 = _mm_load_sd(&A[9]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_2 = _mm_add_sd(c1_2, _mm_mul_sd(a1_2, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_2 = _mm_add_sd(c1_2, _mm_mul_sd(a1_2, b1));
#endif
    _mm_store_sd(&C[(l_n*84)+17], c1_2);
    __m128d c1_3 = _mm_load_sd(&C[(l_n*84)+32]);
    __m128d a1_3 = _mm_load_sd(&A[10]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_3 = _mm_add_sd(c1_3, _mm_mul_sd(a1_3, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_3 = _mm_add_sd(c1_3, _mm_mul_sd(a1_3, b1));
#endif
    _mm_store_sd(&C[(l_n*84)+32], c1_3);
    __m128d c1_4 = _mm_load_sd(&C[(l_n*84)+53]);
    __m128d a1_4 = _mm_load_sd(&A[11]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_4 = _mm_add_sd(c1_4, _mm_mul_sd(a1_4, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_4 = _mm_add_sd(c1_4, _mm_mul_sd(a1_4, b1));
#endif
    _mm_store_sd(&C[(l_n*84)+53], c1_4);
    __m128d c1_5 = _mm_load_sd(&C[(l_n*84)+81]);
    __m128d a1_5 = _mm_load_sd(&A[12]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_5 = _mm_add_sd(c1_5, _mm_mul_sd(a1_5, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_5 = _mm_add_sd(c1_5, _mm_mul_sd(a1_5, b1));
#endif
    _mm_store_sd(&C[(l_n*84)+81], c1_5);
#else
    C[(l_n*84)+1] += A[7] * B[(l_n*84)+1];
    C[(l_n*84)+7] += A[8] * B[(l_n*84)+1];
    C[(l_n*84)+17] += A[9] * B[(l_n*84)+1];
    C[(l_n*84)+32] += A[10] * B[(l_n*84)+1];
    C[(l_n*84)+53] += A[11] * B[(l_n*84)+1];
    C[(l_n*84)+81] += A[12] * B[(l_n*84)+1];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b2 = _mm256_broadcast_sd(&B[(l_n*84)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b2 = _mm_loaddup_pd(&B[(l_n*84)+2]);
#endif
    __m128d c2_0 = _mm_load_sd(&C[(l_n*84)+2]);
    __m128d a2_0 = _mm_load_sd(&A[13]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_0 = _mm_add_sd(c2_0, _mm_mul_sd(a2_0, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_0 = _mm_add_sd(c2_0, _mm_mul_sd(a2_0, b2));
#endif
    _mm_store_sd(&C[(l_n*84)+2], c2_0);
    __m128d c2_1 = _mm_load_sd(&C[(l_n*84)+8]);
    __m128d a2_1 = _mm_load_sd(&A[14]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_1 = _mm_add_sd(c2_1, _mm_mul_sd(a2_1, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_1 = _mm_add_sd(c2_1, _mm_mul_sd(a2_1, b2));
#endif
    _mm_store_sd(&C[(l_n*84)+8], c2_1);
    __m128d c2_2 = _mm_load_sd(&C[(l_n*84)+18]);
    __m128d a2_2 = _mm_load_sd(&A[15]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_2 = _mm_add_sd(c2_2, _mm_mul_sd(a2_2, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_2 = _mm_add_sd(c2_2, _mm_mul_sd(a2_2, b2));
#endif
    _mm_store_sd(&C[(l_n*84)+18], c2_2);
    __m128d c2_3 = _mm_load_sd(&C[(l_n*84)+33]);
    __m128d a2_3 = _mm_load_sd(&A[16]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_3 = _mm_add_sd(c2_3, _mm_mul_sd(a2_3, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_3 = _mm_add_sd(c2_3, _mm_mul_sd(a2_3, b2));
#endif
    _mm_store_sd(&C[(l_n*84)+33], c2_3);
    __m128d c2_4 = _mm_load_sd(&C[(l_n*84)+54]);
    __m128d a2_4 = _mm_load_sd(&A[17]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_4 = _mm_add_sd(c2_4, _mm_mul_sd(a2_4, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_4 = _mm_add_sd(c2_4, _mm_mul_sd(a2_4, b2));
#endif
    _mm_store_sd(&C[(l_n*84)+54], c2_4);
    __m128d c2_5 = _mm_load_sd(&C[(l_n*84)+82]);
    __m128d a2_5 = _mm_load_sd(&A[18]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_5 = _mm_add_sd(c2_5, _mm_mul_sd(a2_5, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_5 = _mm_add_sd(c2_5, _mm_mul_sd(a2_5, b2));
#endif
    _mm_store_sd(&C[(l_n*84)+82], c2_5);
#else
    C[(l_n*84)+2] += A[13] * B[(l_n*84)+2];
    C[(l_n*84)+8] += A[14] * B[(l_n*84)+2];
    C[(l_n*84)+18] += A[15] * B[(l_n*84)+2];
    C[(l_n*84)+33] += A[16] * B[(l_n*84)+2];
    C[(l_n*84)+54] += A[17] * B[(l_n*84)+2];
    C[(l_n*84)+82] += A[18] * B[(l_n*84)+2];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b3 = _mm256_broadcast_sd(&B[(l_n*84)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b3 = _mm_loaddup_pd(&B[(l_n*84)+3]);
#endif
    __m128d c3_0 = _mm_load_sd(&C[(l_n*84)+0]);
    __m128d a3_0 = _mm_load_sd(&A[19]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, b3));
#endif
    _mm_store_sd(&C[(l_n*84)+0], c3_0);
    __m128d c3_1 = _mm_load_sd(&C[(l_n*84)+3]);
    __m128d a3_1 = _mm_load_sd(&A[20]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, b3));
#endif
    _mm_store_sd(&C[(l_n*84)+3], c3_1);
    __m128d c3_2 = _mm_load_sd(&C[(l_n*84)+9]);
    __m128d a3_2 = _mm_load_sd(&A[21]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, b3));
#endif
    _mm_store_sd(&C[(l_n*84)+9], c3_2);
    __m128d c3_3 = _mm_load_sd(&C[(l_n*84)+19]);
    __m128d a3_3 = _mm_load_sd(&A[22]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_3 = _mm_add_sd(c3_3, _mm_mul_sd(a3_3, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_3 = _mm_add_sd(c3_3, _mm_mul_sd(a3_3, b3));
#endif
    _mm_store_sd(&C[(l_n*84)+19], c3_3);
    __m128d c3_4 = _mm_load_sd(&C[(l_n*84)+34]);
    __m128d a3_4 = _mm_load_sd(&A[23]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_4 = _mm_add_sd(c3_4, _mm_mul_sd(a3_4, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_4 = _mm_add_sd(c3_4, _mm_mul_sd(a3_4, b3));
#endif
    _mm_store_sd(&C[(l_n*84)+34], c3_4);
    __m128d c3_5 = _mm_load_sd(&C[(l_n*84)+55]);
    __m128d a3_5 = _mm_load_sd(&A[24]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_5 = _mm_add_sd(c3_5, _mm_mul_sd(a3_5, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_5 = _mm_add_sd(c3_5, _mm_mul_sd(a3_5, b3));
#endif
    _mm_store_sd(&C[(l_n*84)+55], c3_5);
    __m128d c3_6 = _mm_load_sd(&C[(l_n*84)+83]);
    __m128d a3_6 = _mm_load_sd(&A[25]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_6 = _mm_add_sd(c3_6, _mm_mul_sd(a3_6, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_6 = _mm_add_sd(c3_6, _mm_mul_sd(a3_6, b3));
#endif
    _mm_store_sd(&C[(l_n*84)+83], c3_6);
#else
    C[(l_n*84)+0] += A[19] * B[(l_n*84)+3];
    C[(l_n*84)+3] += A[20] * B[(l_n*84)+3];
    C[(l_n*84)+9] += A[21] * B[(l_n*84)+3];
    C[(l_n*84)+19] += A[22] * B[(l_n*84)+3];
    C[(l_n*84)+34] += A[23] * B[(l_n*84)+3];
    C[(l_n*84)+55] += A[24] * B[(l_n*84)+3];
    C[(l_n*84)+83] += A[25] * B[(l_n*84)+3];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b4 = _mm256_broadcast_sd(&B[(l_n*84)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b4 = _mm_loaddup_pd(&B[(l_n*84)+4]);
#endif
    __m128d c4_0 = _mm_load_sd(&C[(l_n*84)+4]);
    __m128d a4_0 = _mm_load_sd(&A[26]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_0 = _mm_add_sd(c4_0, _mm_mul_sd(a4_0, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_0 = _mm_add_sd(c4_0, _mm_mul_sd(a4_0, b4));
#endif
    _mm_store_sd(&C[(l_n*84)+4], c4_0);
    __m128d c4_1 = _mm_load_sd(&C[(l_n*84)+14]);
    __m128d a4_1 = _mm_load_sd(&A[27]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_1 = _mm_add_sd(c4_1, _mm_mul_sd(a4_1, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_1 = _mm_add_sd(c4_1, _mm_mul_sd(a4_1, b4));
#endif
    _mm_store_sd(&C[(l_n*84)+14], c4_1);
    __m128d c4_2 = _mm_load_sd(&C[(l_n*84)+29]);
    __m128d a4_2 = _mm_load_sd(&A[28]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, b4));
#endif
    _mm_store_sd(&C[(l_n*84)+29], c4_2);
    __m128d c4_3 = _mm_load_sd(&C[(l_n*84)+50]);
    __m128d a4_3 = _mm_load_sd(&A[29]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_3 = _mm_add_sd(c4_3, _mm_mul_sd(a4_3, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_3 = _mm_add_sd(c4_3, _mm_mul_sd(a4_3, b4));
#endif
    _mm_store_sd(&C[(l_n*84)+50], c4_3);
    __m128d c4_4 = _mm_load_sd(&C[(l_n*84)+78]);
    __m128d a4_4 = _mm_load_sd(&A[30]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_4 = _mm_add_sd(c4_4, _mm_mul_sd(a4_4, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_4 = _mm_add_sd(c4_4, _mm_mul_sd(a4_4, b4));
#endif
    _mm_store_sd(&C[(l_n*84)+78], c4_4);
#else
    C[(l_n*84)+4] += A[26] * B[(l_n*84)+4];
    C[(l_n*84)+14] += A[27] * B[(l_n*84)+4];
    C[(l_n*84)+29] += A[28] * B[(l_n*84)+4];
    C[(l_n*84)+50] += A[29] * B[(l_n*84)+4];
    C[(l_n*84)+78] += A[30] * B[(l_n*84)+4];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b5 = _mm256_broadcast_sd(&B[(l_n*84)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b5 = _mm_loaddup_pd(&B[(l_n*84)+5]);
#endif
    __m128d c5_0 = _mm_load_sd(&C[(l_n*84)+5]);
    __m128d a5_0 = _mm_load_sd(&A[31]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_0 = _mm_add_sd(c5_0, _mm_mul_sd(a5_0, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_0 = _mm_add_sd(c5_0, _mm_mul_sd(a5_0, b5));
#endif
    _mm_store_sd(&C[(l_n*84)+5], c5_0);
    __m128d c5_1 = _mm_load_sd(&C[(l_n*84)+15]);
    __m128d a5_1 = _mm_load_sd(&A[32]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_1 = _mm_add_sd(c5_1, _mm_mul_sd(a5_1, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_1 = _mm_add_sd(c5_1, _mm_mul_sd(a5_1, b5));
#endif
    _mm_store_sd(&C[(l_n*84)+15], c5_1);
    __m128d c5_2 = _mm_load_sd(&C[(l_n*84)+30]);
    __m128d a5_2 = _mm_load_sd(&A[33]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, b5));
#endif
    _mm_store_sd(&C[(l_n*84)+30], c5_2);
    __m128d c5_3 = _mm_load_sd(&C[(l_n*84)+51]);
    __m128d a5_3 = _mm_load_sd(&A[34]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_3 = _mm_add_sd(c5_3, _mm_mul_sd(a5_3, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_3 = _mm_add_sd(c5_3, _mm_mul_sd(a5_3, b5));
#endif
    _mm_store_sd(&C[(l_n*84)+51], c5_3);
    __m128d c5_4 = _mm_load_sd(&C[(l_n*84)+79]);
    __m128d a5_4 = _mm_load_sd(&A[35]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_4 = _mm_add_sd(c5_4, _mm_mul_sd(a5_4, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_4 = _mm_add_sd(c5_4, _mm_mul_sd(a5_4, b5));
#endif
    _mm_store_sd(&C[(l_n*84)+79], c5_4);
#else
    C[(l_n*84)+5] += A[31] * B[(l_n*84)+5];
    C[(l_n*84)+15] += A[32] * B[(l_n*84)+5];
    C[(l_n*84)+30] += A[33] * B[(l_n*84)+5];
    C[(l_n*84)+51] += A[34] * B[(l_n*84)+5];
    C[(l_n*84)+79] += A[35] * B[(l_n*84)+5];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b6 = _mm256_broadcast_sd(&B[(l_n*84)+6]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b6 = _mm_loaddup_pd(&B[(l_n*84)+6]);
#endif
    __m128d c6_0 = _mm_load_sd(&C[(l_n*84)+6]);
    __m128d a6_0 = _mm_load_sd(&A[36]);
#if defined(__SSE3__) && defined(__AVX__)
    c6_0 = _mm_add_sd(c6_0, _mm_mul_sd(a6_0, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c6_0 = _mm_add_sd(c6_0, _mm_mul_sd(a6_0, b6));
#endif
    _mm_store_sd(&C[(l_n*84)+6], c6_0);
    __m128d c6_1 = _mm_load_sd(&C[(l_n*84)+16]);
    __m128d a6_1 = _mm_load_sd(&A[37]);
#if defined(__SSE3__) && defined(__AVX__)
    c6_1 = _mm_add_sd(c6_1, _mm_mul_sd(a6_1, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c6_1 = _mm_add_sd(c6_1, _mm_mul_sd(a6_1, b6));
#endif
    _mm_store_sd(&C[(l_n*84)+16], c6_1);
    __m128d c6_2 = _mm_load_sd(&C[(l_n*84)+31]);
    __m128d a6_2 = _mm_load_sd(&A[38]);
#if defined(__SSE3__) && defined(__AVX__)
    c6_2 = _mm_add_sd(c6_2, _mm_mul_sd(a6_2, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c6_2 = _mm_add_sd(c6_2, _mm_mul_sd(a6_2, b6));
#endif
    _mm_store_sd(&C[(l_n*84)+31], c6_2);
    __m128d c6_3 = _mm_load_sd(&C[(l_n*84)+52]);
    __m128d a6_3 = _mm_load_sd(&A[39]);
#if defined(__SSE3__) && defined(__AVX__)
    c6_3 = _mm_add_sd(c6_3, _mm_mul_sd(a6_3, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c6_3 = _mm_add_sd(c6_3, _mm_mul_sd(a6_3, b6));
#endif
    _mm_store_sd(&C[(l_n*84)+52], c6_3);
    __m128d c6_4 = _mm_load_sd(&C[(l_n*84)+80]);
    __m128d a6_4 = _mm_load_sd(&A[40]);
#if defined(__SSE3__) && defined(__AVX__)
    c6_4 = _mm_add_sd(c6_4, _mm_mul_sd(a6_4, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c6_4 = _mm_add_sd(c6_4, _mm_mul_sd(a6_4, b6));
#endif
    _mm_store_sd(&C[(l_n*84)+80], c6_4);
#else
    C[(l_n*84)+6] += A[36] * B[(l_n*84)+6];
    C[(l_n*84)+16] += A[37] * B[(l_n*84)+6];
    C[(l_n*84)+31] += A[38] * B[(l_n*84)+6];
    C[(l_n*84)+52] += A[39] * B[(l_n*84)+6];
    C[(l_n*84)+80] += A[40] * B[(l_n*84)+6];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b7 = _mm256_broadcast_sd(&B[(l_n*84)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b7 = _mm_loaddup_pd(&B[(l_n*84)+7]);
#endif
    __m128d c7_0 = _mm_load_sd(&C[(l_n*84)+1]);
    __m128d a7_0 = _mm_load_sd(&A[41]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_0 = _mm_add_sd(c7_0, _mm_mul_sd(a7_0, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_0 = _mm_add_sd(c7_0, _mm_mul_sd(a7_0, b7));
#endif
    _mm_store_sd(&C[(l_n*84)+1], c7_0);
    __m128d c7_1 = _mm_load_sd(&C[(l_n*84)+7]);
    __m128d a7_1 = _mm_load_sd(&A[42]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_1 = _mm_add_sd(c7_1, _mm_mul_sd(a7_1, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_1 = _mm_add_sd(c7_1, _mm_mul_sd(a7_1, b7));
#endif
    _mm_store_sd(&C[(l_n*84)+7], c7_1);
    __m128d c7_2 = _mm_load_sd(&C[(l_n*84)+17]);
    __m128d a7_2 = _mm_load_sd(&A[43]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_2 = _mm_add_sd(c7_2, _mm_mul_sd(a7_2, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_2 = _mm_add_sd(c7_2, _mm_mul_sd(a7_2, b7));
#endif
    _mm_store_sd(&C[(l_n*84)+17], c7_2);
    __m128d c7_3 = _mm_load_sd(&C[(l_n*84)+32]);
    __m128d a7_3 = _mm_load_sd(&A[44]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_3 = _mm_add_sd(c7_3, _mm_mul_sd(a7_3, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_3 = _mm_add_sd(c7_3, _mm_mul_sd(a7_3, b7));
#endif
    _mm_store_sd(&C[(l_n*84)+32], c7_3);
    __m128d c7_4 = _mm_load_sd(&C[(l_n*84)+53]);
    __m128d a7_4 = _mm_load_sd(&A[45]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_4 = _mm_add_sd(c7_4, _mm_mul_sd(a7_4, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_4 = _mm_add_sd(c7_4, _mm_mul_sd(a7_4, b7));
#endif
    _mm_store_sd(&C[(l_n*84)+53], c7_4);
    __m128d c7_5 = _mm_load_sd(&C[(l_n*84)+81]);
    __m128d a7_5 = _mm_load_sd(&A[46]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_5 = _mm_add_sd(c7_5, _mm_mul_sd(a7_5, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_5 = _mm_add_sd(c7_5, _mm_mul_sd(a7_5, b7));
#endif
    _mm_store_sd(&C[(l_n*84)+81], c7_5);
#else
    C[(l_n*84)+1] += A[41] * B[(l_n*84)+7];
    C[(l_n*84)+7] += A[42] * B[(l_n*84)+7];
    C[(l_n*84)+17] += A[43] * B[(l_n*84)+7];
    C[(l_n*84)+32] += A[44] * B[(l_n*84)+7];
    C[(l_n*84)+53] += A[45] * B[(l_n*84)+7];
    C[(l_n*84)+81] += A[46] * B[(l_n*84)+7];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b8 = _mm256_broadcast_sd(&B[(l_n*84)+8]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b8 = _mm_loaddup_pd(&B[(l_n*84)+8]);
#endif
    __m128d c8_0 = _mm_load_sd(&C[(l_n*84)+2]);
    __m128d a8_0 = _mm_load_sd(&A[47]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_0 = _mm_add_sd(c8_0, _mm_mul_sd(a8_0, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_0 = _mm_add_sd(c8_0, _mm_mul_sd(a8_0, b8));
#endif
    _mm_store_sd(&C[(l_n*84)+2], c8_0);
    __m128d c8_1 = _mm_load_sd(&C[(l_n*84)+8]);
    __m128d a8_1 = _mm_load_sd(&A[48]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_1 = _mm_add_sd(c8_1, _mm_mul_sd(a8_1, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_1 = _mm_add_sd(c8_1, _mm_mul_sd(a8_1, b8));
#endif
    _mm_store_sd(&C[(l_n*84)+8], c8_1);
    __m128d c8_2 = _mm_load_sd(&C[(l_n*84)+18]);
    __m128d a8_2 = _mm_load_sd(&A[49]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_2 = _mm_add_sd(c8_2, _mm_mul_sd(a8_2, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_2 = _mm_add_sd(c8_2, _mm_mul_sd(a8_2, b8));
#endif
    _mm_store_sd(&C[(l_n*84)+18], c8_2);
    __m128d c8_3 = _mm_load_sd(&C[(l_n*84)+33]);
    __m128d a8_3 = _mm_load_sd(&A[50]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_3 = _mm_add_sd(c8_3, _mm_mul_sd(a8_3, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_3 = _mm_add_sd(c8_3, _mm_mul_sd(a8_3, b8));
#endif
    _mm_store_sd(&C[(l_n*84)+33], c8_3);
    __m128d c8_4 = _mm_load_sd(&C[(l_n*84)+54]);
    __m128d a8_4 = _mm_load_sd(&A[51]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_4 = _mm_add_sd(c8_4, _mm_mul_sd(a8_4, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_4 = _mm_add_sd(c8_4, _mm_mul_sd(a8_4, b8));
#endif
    _mm_store_sd(&C[(l_n*84)+54], c8_4);
    __m128d c8_5 = _mm_load_sd(&C[(l_n*84)+82]);
    __m128d a8_5 = _mm_load_sd(&A[52]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_5 = _mm_add_sd(c8_5, _mm_mul_sd(a8_5, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_5 = _mm_add_sd(c8_5, _mm_mul_sd(a8_5, b8));
#endif
    _mm_store_sd(&C[(l_n*84)+82], c8_5);
#else
    C[(l_n*84)+2] += A[47] * B[(l_n*84)+8];
    C[(l_n*84)+8] += A[48] * B[(l_n*84)+8];
    C[(l_n*84)+18] += A[49] * B[(l_n*84)+8];
    C[(l_n*84)+33] += A[50] * B[(l_n*84)+8];
    C[(l_n*84)+54] += A[51] * B[(l_n*84)+8];
    C[(l_n*84)+82] += A[52] * B[(l_n*84)+8];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b9 = _mm256_broadcast_sd(&B[(l_n*84)+9]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b9 = _mm_loaddup_pd(&B[(l_n*84)+9]);
#endif
    __m128d c9_0 = _mm_load_sd(&C[(l_n*84)+0]);
    __m128d a9_0 = _mm_load_sd(&A[53]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, b9));
#endif
    _mm_store_sd(&C[(l_n*84)+0], c9_0);
    __m128d c9_1 = _mm_load_sd(&C[(l_n*84)+3]);
    __m128d a9_1 = _mm_load_sd(&A[54]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, b9));
#endif
    _mm_store_sd(&C[(l_n*84)+3], c9_1);
    __m128d c9_2 = _mm_load_sd(&C[(l_n*84)+9]);
    __m128d a9_2 = _mm_load_sd(&A[55]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, b9));
#endif
    _mm_store_sd(&C[(l_n*84)+9], c9_2);
    __m128d c9_3 = _mm_load_sd(&C[(l_n*84)+19]);
    __m128d a9_3 = _mm_load_sd(&A[56]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_3 = _mm_add_sd(c9_3, _mm_mul_sd(a9_3, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_3 = _mm_add_sd(c9_3, _mm_mul_sd(a9_3, b9));
#endif
    _mm_store_sd(&C[(l_n*84)+19], c9_3);
    __m128d c9_4 = _mm_load_sd(&C[(l_n*84)+34]);
    __m128d a9_4 = _mm_load_sd(&A[57]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_4 = _mm_add_sd(c9_4, _mm_mul_sd(a9_4, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_4 = _mm_add_sd(c9_4, _mm_mul_sd(a9_4, b9));
#endif
    _mm_store_sd(&C[(l_n*84)+34], c9_4);
    __m128d c9_5 = _mm_load_sd(&C[(l_n*84)+55]);
    __m128d a9_5 = _mm_load_sd(&A[58]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_5 = _mm_add_sd(c9_5, _mm_mul_sd(a9_5, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_5 = _mm_add_sd(c9_5, _mm_mul_sd(a9_5, b9));
#endif
    _mm_store_sd(&C[(l_n*84)+55], c9_5);
    __m128d c9_6 = _mm_load_sd(&C[(l_n*84)+83]);
    __m128d a9_6 = _mm_load_sd(&A[59]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_6 = _mm_add_sd(c9_6, _mm_mul_sd(a9_6, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_6 = _mm_add_sd(c9_6, _mm_mul_sd(a9_6, b9));
#endif
    _mm_store_sd(&C[(l_n*84)+83], c9_6);
#else
    C[(l_n*84)+0] += A[53] * B[(l_n*84)+9];
    C[(l_n*84)+3] += A[54] * B[(l_n*84)+9];
    C[(l_n*84)+9] += A[55] * B[(l_n*84)+9];
    C[(l_n*84)+19] += A[56] * B[(l_n*84)+9];
    C[(l_n*84)+34] += A[57] * B[(l_n*84)+9];
    C[(l_n*84)+55] += A[58] * B[(l_n*84)+9];
    C[(l_n*84)+83] += A[59] * B[(l_n*84)+9];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b10 = _mm256_broadcast_sd(&B[(l_n*84)+10]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b10 = _mm_loaddup_pd(&B[(l_n*84)+10]);
#endif
    __m128d c10_0 = _mm_load_sd(&C[(l_n*84)+10]);
    __m128d a10_0 = _mm_load_sd(&A[60]);
#if defined(__SSE3__) && defined(__AVX__)
    c10_0 = _mm_add_sd(c10_0, _mm_mul_sd(a10_0, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c10_0 = _mm_add_sd(c10_0, _mm_mul_sd(a10_0, b10));
#endif
    _mm_store_sd(&C[(l_n*84)+10], c10_0);
    __m128d c10_1 = _mm_load_sd(&C[(l_n*84)+25]);
    __m128d a10_1 = _mm_load_sd(&A[61]);
#if defined(__SSE3__) && defined(__AVX__)
    c10_1 = _mm_add_sd(c10_1, _mm_mul_sd(a10_1, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c10_1 = _mm_add_sd(c10_1, _mm_mul_sd(a10_1, b10));
#endif
    _mm_store_sd(&C[(l_n*84)+25], c10_1);
    __m128d c10_2 = _mm_load_sd(&C[(l_n*84)+46]);
    __m128d a10_2 = _mm_load_sd(&A[62]);
#if defined(__SSE3__) && defined(__AVX__)
    c10_2 = _mm_add_sd(c10_2, _mm_mul_sd(a10_2, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c10_2 = _mm_add_sd(c10_2, _mm_mul_sd(a10_2, b10));
#endif
    _mm_store_sd(&C[(l_n*84)+46], c10_2);
    __m128d c10_3 = _mm_load_sd(&C[(l_n*84)+74]);
    __m128d a10_3 = _mm_load_sd(&A[63]);
#if defined(__SSE3__) && defined(__AVX__)
    c10_3 = _mm_add_sd(c10_3, _mm_mul_sd(a10_3, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c10_3 = _mm_add_sd(c10_3, _mm_mul_sd(a10_3, b10));
#endif
    _mm_store_sd(&C[(l_n*84)+74], c10_3);
#else
    C[(l_n*84)+10] += A[60] * B[(l_n*84)+10];
    C[(l_n*84)+25] += A[61] * B[(l_n*84)+10];
    C[(l_n*84)+46] += A[62] * B[(l_n*84)+10];
    C[(l_n*84)+74] += A[63] * B[(l_n*84)+10];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b11 = _mm256_broadcast_sd(&B[(l_n*84)+11]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b11 = _mm_loaddup_pd(&B[(l_n*84)+11]);
#endif
    __m128d c11_0 = _mm_load_sd(&C[(l_n*84)+11]);
    __m128d a11_0 = _mm_load_sd(&A[64]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_0 = _mm_add_sd(c11_0, _mm_mul_sd(a11_0, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_0 = _mm_add_sd(c11_0, _mm_mul_sd(a11_0, b11));
#endif
    _mm_store_sd(&C[(l_n*84)+11], c11_0);
    __m128d c11_1 = _mm_load_sd(&C[(l_n*84)+26]);
    __m128d a11_1 = _mm_load_sd(&A[65]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_1 = _mm_add_sd(c11_1, _mm_mul_sd(a11_1, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_1 = _mm_add_sd(c11_1, _mm_mul_sd(a11_1, b11));
#endif
    _mm_store_sd(&C[(l_n*84)+26], c11_1);
    __m128d c11_2 = _mm_load_sd(&C[(l_n*84)+47]);
    __m128d a11_2 = _mm_load_sd(&A[66]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_2 = _mm_add_sd(c11_2, _mm_mul_sd(a11_2, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_2 = _mm_add_sd(c11_2, _mm_mul_sd(a11_2, b11));
#endif
    _mm_store_sd(&C[(l_n*84)+47], c11_2);
    __m128d c11_3 = _mm_load_sd(&C[(l_n*84)+75]);
    __m128d a11_3 = _mm_load_sd(&A[67]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_3 = _mm_add_sd(c11_3, _mm_mul_sd(a11_3, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_3 = _mm_add_sd(c11_3, _mm_mul_sd(a11_3, b11));
#endif
    _mm_store_sd(&C[(l_n*84)+75], c11_3);
#else
    C[(l_n*84)+11] += A[64] * B[(l_n*84)+11];
    C[(l_n*84)+26] += A[65] * B[(l_n*84)+11];
    C[(l_n*84)+47] += A[66] * B[(l_n*84)+11];
    C[(l_n*84)+75] += A[67] * B[(l_n*84)+11];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b12 = _mm256_broadcast_sd(&B[(l_n*84)+12]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b12 = _mm_loaddup_pd(&B[(l_n*84)+12]);
#endif
    __m128d c12_0 = _mm_load_sd(&C[(l_n*84)+12]);
    __m128d a12_0 = _mm_load_sd(&A[68]);
#if defined(__SSE3__) && defined(__AVX__)
    c12_0 = _mm_add_sd(c12_0, _mm_mul_sd(a12_0, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c12_0 = _mm_add_sd(c12_0, _mm_mul_sd(a12_0, b12));
#endif
    _mm_store_sd(&C[(l_n*84)+12], c12_0);
    __m128d c12_1 = _mm_load_sd(&C[(l_n*84)+27]);
    __m128d a12_1 = _mm_load_sd(&A[69]);
#if defined(__SSE3__) && defined(__AVX__)
    c12_1 = _mm_add_sd(c12_1, _mm_mul_sd(a12_1, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c12_1 = _mm_add_sd(c12_1, _mm_mul_sd(a12_1, b12));
#endif
    _mm_store_sd(&C[(l_n*84)+27], c12_1);
    __m128d c12_2 = _mm_load_sd(&C[(l_n*84)+48]);
    __m128d a12_2 = _mm_load_sd(&A[70]);
#if defined(__SSE3__) && defined(__AVX__)
    c12_2 = _mm_add_sd(c12_2, _mm_mul_sd(a12_2, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c12_2 = _mm_add_sd(c12_2, _mm_mul_sd(a12_2, b12));
#endif
    _mm_store_sd(&C[(l_n*84)+48], c12_2);
    __m128d c12_3 = _mm_load_sd(&C[(l_n*84)+76]);
    __m128d a12_3 = _mm_load_sd(&A[71]);
#if defined(__SSE3__) && defined(__AVX__)
    c12_3 = _mm_add_sd(c12_3, _mm_mul_sd(a12_3, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c12_3 = _mm_add_sd(c12_3, _mm_mul_sd(a12_3, b12));
#endif
    _mm_store_sd(&C[(l_n*84)+76], c12_3);
#else
    C[(l_n*84)+12] += A[68] * B[(l_n*84)+12];
    C[(l_n*84)+27] += A[69] * B[(l_n*84)+12];
    C[(l_n*84)+48] += A[70] * B[(l_n*84)+12];
    C[(l_n*84)+76] += A[71] * B[(l_n*84)+12];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b13 = _mm256_broadcast_sd(&B[(l_n*84)+13]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b13 = _mm_loaddup_pd(&B[(l_n*84)+13]);
#endif
    __m128d c13_0 = _mm_load_sd(&C[(l_n*84)+13]);
    __m128d a13_0 = _mm_load_sd(&A[72]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_0 = _mm_add_sd(c13_0, _mm_mul_sd(a13_0, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_0 = _mm_add_sd(c13_0, _mm_mul_sd(a13_0, b13));
#endif
    _mm_store_sd(&C[(l_n*84)+13], c13_0);
    __m128d c13_1 = _mm_load_sd(&C[(l_n*84)+28]);
    __m128d a13_1 = _mm_load_sd(&A[73]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_1 = _mm_add_sd(c13_1, _mm_mul_sd(a13_1, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_1 = _mm_add_sd(c13_1, _mm_mul_sd(a13_1, b13));
#endif
    _mm_store_sd(&C[(l_n*84)+28], c13_1);
    __m128d c13_2 = _mm_load_sd(&C[(l_n*84)+49]);
    __m128d a13_2 = _mm_load_sd(&A[74]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_2 = _mm_add_sd(c13_2, _mm_mul_sd(a13_2, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_2 = _mm_add_sd(c13_2, _mm_mul_sd(a13_2, b13));
#endif
    _mm_store_sd(&C[(l_n*84)+49], c13_2);
    __m128d c13_3 = _mm_load_sd(&C[(l_n*84)+77]);
    __m128d a13_3 = _mm_load_sd(&A[75]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_3 = _mm_add_sd(c13_3, _mm_mul_sd(a13_3, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_3 = _mm_add_sd(c13_3, _mm_mul_sd(a13_3, b13));
#endif
    _mm_store_sd(&C[(l_n*84)+77], c13_3);
#else
    C[(l_n*84)+13] += A[72] * B[(l_n*84)+13];
    C[(l_n*84)+28] += A[73] * B[(l_n*84)+13];
    C[(l_n*84)+49] += A[74] * B[(l_n*84)+13];
    C[(l_n*84)+77] += A[75] * B[(l_n*84)+13];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b14 = _mm256_broadcast_sd(&B[(l_n*84)+14]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b14 = _mm_loaddup_pd(&B[(l_n*84)+14]);
#endif
    __m128d c14_0 = _mm_load_sd(&C[(l_n*84)+4]);
    __m128d a14_0 = _mm_load_sd(&A[76]);
#if defined(__SSE3__) && defined(__AVX__)
    c14_0 = _mm_add_sd(c14_0, _mm_mul_sd(a14_0, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c14_0 = _mm_add_sd(c14_0, _mm_mul_sd(a14_0, b14));
#endif
    _mm_store_sd(&C[(l_n*84)+4], c14_0);
    __m128d c14_1 = _mm_load_sd(&C[(l_n*84)+14]);
    __m128d a14_1 = _mm_load_sd(&A[77]);
#if defined(__SSE3__) && defined(__AVX__)
    c14_1 = _mm_add_sd(c14_1, _mm_mul_sd(a14_1, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c14_1 = _mm_add_sd(c14_1, _mm_mul_sd(a14_1, b14));
#endif
    _mm_store_sd(&C[(l_n*84)+14], c14_1);
    __m128d c14_2 = _mm_load_sd(&C[(l_n*84)+29]);
    __m128d a14_2 = _mm_load_sd(&A[78]);
#if defined(__SSE3__) && defined(__AVX__)
    c14_2 = _mm_add_sd(c14_2, _mm_mul_sd(a14_2, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c14_2 = _mm_add_sd(c14_2, _mm_mul_sd(a14_2, b14));
#endif
    _mm_store_sd(&C[(l_n*84)+29], c14_2);
    __m128d c14_3 = _mm_load_sd(&C[(l_n*84)+50]);
    __m128d a14_3 = _mm_load_sd(&A[79]);
#if defined(__SSE3__) && defined(__AVX__)
    c14_3 = _mm_add_sd(c14_3, _mm_mul_sd(a14_3, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c14_3 = _mm_add_sd(c14_3, _mm_mul_sd(a14_3, b14));
#endif
    _mm_store_sd(&C[(l_n*84)+50], c14_3);
    __m128d c14_4 = _mm_load_sd(&C[(l_n*84)+78]);
    __m128d a14_4 = _mm_load_sd(&A[80]);
#if defined(__SSE3__) && defined(__AVX__)
    c14_4 = _mm_add_sd(c14_4, _mm_mul_sd(a14_4, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c14_4 = _mm_add_sd(c14_4, _mm_mul_sd(a14_4, b14));
#endif
    _mm_store_sd(&C[(l_n*84)+78], c14_4);
#else
    C[(l_n*84)+4] += A[76] * B[(l_n*84)+14];
    C[(l_n*84)+14] += A[77] * B[(l_n*84)+14];
    C[(l_n*84)+29] += A[78] * B[(l_n*84)+14];
    C[(l_n*84)+50] += A[79] * B[(l_n*84)+14];
    C[(l_n*84)+78] += A[80] * B[(l_n*84)+14];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b15 = _mm256_broadcast_sd(&B[(l_n*84)+15]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b15 = _mm_loaddup_pd(&B[(l_n*84)+15]);
#endif
    __m128d c15_0 = _mm_load_sd(&C[(l_n*84)+5]);
    __m128d a15_0 = _mm_load_sd(&A[81]);
#if defined(__SSE3__) && defined(__AVX__)
    c15_0 = _mm_add_sd(c15_0, _mm_mul_sd(a15_0, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c15_0 = _mm_add_sd(c15_0, _mm_mul_sd(a15_0, b15));
#endif
    _mm_store_sd(&C[(l_n*84)+5], c15_0);
    __m128d c15_1 = _mm_load_sd(&C[(l_n*84)+15]);
    __m128d a15_1 = _mm_load_sd(&A[82]);
#if defined(__SSE3__) && defined(__AVX__)
    c15_1 = _mm_add_sd(c15_1, _mm_mul_sd(a15_1, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c15_1 = _mm_add_sd(c15_1, _mm_mul_sd(a15_1, b15));
#endif
    _mm_store_sd(&C[(l_n*84)+15], c15_1);
    __m128d c15_2 = _mm_load_sd(&C[(l_n*84)+30]);
    __m128d a15_2 = _mm_load_sd(&A[83]);
#if defined(__SSE3__) && defined(__AVX__)
    c15_2 = _mm_add_sd(c15_2, _mm_mul_sd(a15_2, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c15_2 = _mm_add_sd(c15_2, _mm_mul_sd(a15_2, b15));
#endif
    _mm_store_sd(&C[(l_n*84)+30], c15_2);
    __m128d c15_3 = _mm_load_sd(&C[(l_n*84)+51]);
    __m128d a15_3 = _mm_load_sd(&A[84]);
#if defined(__SSE3__) && defined(__AVX__)
    c15_3 = _mm_add_sd(c15_3, _mm_mul_sd(a15_3, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c15_3 = _mm_add_sd(c15_3, _mm_mul_sd(a15_3, b15));
#endif
    _mm_store_sd(&C[(l_n*84)+51], c15_3);
    __m128d c15_4 = _mm_load_sd(&C[(l_n*84)+79]);
    __m128d a15_4 = _mm_load_sd(&A[85]);
#if defined(__SSE3__) && defined(__AVX__)
    c15_4 = _mm_add_sd(c15_4, _mm_mul_sd(a15_4, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c15_4 = _mm_add_sd(c15_4, _mm_mul_sd(a15_4, b15));
#endif
    _mm_store_sd(&C[(l_n*84)+79], c15_4);
#else
    C[(l_n*84)+5] += A[81] * B[(l_n*84)+15];
    C[(l_n*84)+15] += A[82] * B[(l_n*84)+15];
    C[(l_n*84)+30] += A[83] * B[(l_n*84)+15];
    C[(l_n*84)+51] += A[84] * B[(l_n*84)+15];
    C[(l_n*84)+79] += A[85] * B[(l_n*84)+15];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b16 = _mm256_broadcast_sd(&B[(l_n*84)+16]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b16 = _mm_loaddup_pd(&B[(l_n*84)+16]);
#endif
    __m128d c16_0 = _mm_load_sd(&C[(l_n*84)+6]);
    __m128d a16_0 = _mm_load_sd(&A[86]);
#if defined(__SSE3__) && defined(__AVX__)
    c16_0 = _mm_add_sd(c16_0, _mm_mul_sd(a16_0, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c16_0 = _mm_add_sd(c16_0, _mm_mul_sd(a16_0, b16));
#endif
    _mm_store_sd(&C[(l_n*84)+6], c16_0);
    __m128d c16_1 = _mm_load_sd(&C[(l_n*84)+16]);
    __m128d a16_1 = _mm_load_sd(&A[87]);
#if defined(__SSE3__) && defined(__AVX__)
    c16_1 = _mm_add_sd(c16_1, _mm_mul_sd(a16_1, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c16_1 = _mm_add_sd(c16_1, _mm_mul_sd(a16_1, b16));
#endif
    _mm_store_sd(&C[(l_n*84)+16], c16_1);
    __m128d c16_2 = _mm_load_sd(&C[(l_n*84)+31]);
    __m128d a16_2 = _mm_load_sd(&A[88]);
#if defined(__SSE3__) && defined(__AVX__)
    c16_2 = _mm_add_sd(c16_2, _mm_mul_sd(a16_2, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c16_2 = _mm_add_sd(c16_2, _mm_mul_sd(a16_2, b16));
#endif
    _mm_store_sd(&C[(l_n*84)+31], c16_2);
    __m128d c16_3 = _mm_load_sd(&C[(l_n*84)+52]);
    __m128d a16_3 = _mm_load_sd(&A[89]);
#if defined(__SSE3__) && defined(__AVX__)
    c16_3 = _mm_add_sd(c16_3, _mm_mul_sd(a16_3, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c16_3 = _mm_add_sd(c16_3, _mm_mul_sd(a16_3, b16));
#endif
    _mm_store_sd(&C[(l_n*84)+52], c16_3);
    __m128d c16_4 = _mm_load_sd(&C[(l_n*84)+80]);
    __m128d a16_4 = _mm_load_sd(&A[90]);
#if defined(__SSE3__) && defined(__AVX__)
    c16_4 = _mm_add_sd(c16_4, _mm_mul_sd(a16_4, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c16_4 = _mm_add_sd(c16_4, _mm_mul_sd(a16_4, b16));
#endif
    _mm_store_sd(&C[(l_n*84)+80], c16_4);
#else
    C[(l_n*84)+6] += A[86] * B[(l_n*84)+16];
    C[(l_n*84)+16] += A[87] * B[(l_n*84)+16];
    C[(l_n*84)+31] += A[88] * B[(l_n*84)+16];
    C[(l_n*84)+52] += A[89] * B[(l_n*84)+16];
    C[(l_n*84)+80] += A[90] * B[(l_n*84)+16];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b17 = _mm256_broadcast_sd(&B[(l_n*84)+17]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b17 = _mm_loaddup_pd(&B[(l_n*84)+17]);
#endif
    __m128d c17_0 = _mm_load_sd(&C[(l_n*84)+1]);
    __m128d a17_0 = _mm_load_sd(&A[91]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_0 = _mm_add_sd(c17_0, _mm_mul_sd(a17_0, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_0 = _mm_add_sd(c17_0, _mm_mul_sd(a17_0, b17));
#endif
    _mm_store_sd(&C[(l_n*84)+1], c17_0);
    __m128d c17_1 = _mm_load_sd(&C[(l_n*84)+7]);
    __m128d a17_1 = _mm_load_sd(&A[92]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_1 = _mm_add_sd(c17_1, _mm_mul_sd(a17_1, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_1 = _mm_add_sd(c17_1, _mm_mul_sd(a17_1, b17));
#endif
    _mm_store_sd(&C[(l_n*84)+7], c17_1);
    __m128d c17_2 = _mm_load_sd(&C[(l_n*84)+17]);
    __m128d a17_2 = _mm_load_sd(&A[93]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_2 = _mm_add_sd(c17_2, _mm_mul_sd(a17_2, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_2 = _mm_add_sd(c17_2, _mm_mul_sd(a17_2, b17));
#endif
    _mm_store_sd(&C[(l_n*84)+17], c17_2);
    __m128d c17_3 = _mm_load_sd(&C[(l_n*84)+32]);
    __m128d a17_3 = _mm_load_sd(&A[94]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_3 = _mm_add_sd(c17_3, _mm_mul_sd(a17_3, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_3 = _mm_add_sd(c17_3, _mm_mul_sd(a17_3, b17));
#endif
    _mm_store_sd(&C[(l_n*84)+32], c17_3);
    __m128d c17_4 = _mm_load_sd(&C[(l_n*84)+53]);
    __m128d a17_4 = _mm_load_sd(&A[95]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_4 = _mm_add_sd(c17_4, _mm_mul_sd(a17_4, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_4 = _mm_add_sd(c17_4, _mm_mul_sd(a17_4, b17));
#endif
    _mm_store_sd(&C[(l_n*84)+53], c17_4);
    __m128d c17_5 = _mm_load_sd(&C[(l_n*84)+81]);
    __m128d a17_5 = _mm_load_sd(&A[96]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_5 = _mm_add_sd(c17_5, _mm_mul_sd(a17_5, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_5 = _mm_add_sd(c17_5, _mm_mul_sd(a17_5, b17));
#endif
    _mm_store_sd(&C[(l_n*84)+81], c17_5);
#else
    C[(l_n*84)+1] += A[91] * B[(l_n*84)+17];
    C[(l_n*84)+7] += A[92] * B[(l_n*84)+17];
    C[(l_n*84)+17] += A[93] * B[(l_n*84)+17];
    C[(l_n*84)+32] += A[94] * B[(l_n*84)+17];
    C[(l_n*84)+53] += A[95] * B[(l_n*84)+17];
    C[(l_n*84)+81] += A[96] * B[(l_n*84)+17];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b18 = _mm256_broadcast_sd(&B[(l_n*84)+18]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b18 = _mm_loaddup_pd(&B[(l_n*84)+18]);
#endif
    __m128d c18_0 = _mm_load_sd(&C[(l_n*84)+2]);
    __m128d a18_0 = _mm_load_sd(&A[97]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_0 = _mm_add_sd(c18_0, _mm_mul_sd(a18_0, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_0 = _mm_add_sd(c18_0, _mm_mul_sd(a18_0, b18));
#endif
    _mm_store_sd(&C[(l_n*84)+2], c18_0);
    __m128d c18_1 = _mm_load_sd(&C[(l_n*84)+8]);
    __m128d a18_1 = _mm_load_sd(&A[98]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_1 = _mm_add_sd(c18_1, _mm_mul_sd(a18_1, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_1 = _mm_add_sd(c18_1, _mm_mul_sd(a18_1, b18));
#endif
    _mm_store_sd(&C[(l_n*84)+8], c18_1);
    __m128d c18_2 = _mm_load_sd(&C[(l_n*84)+18]);
    __m128d a18_2 = _mm_load_sd(&A[99]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_2 = _mm_add_sd(c18_2, _mm_mul_sd(a18_2, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_2 = _mm_add_sd(c18_2, _mm_mul_sd(a18_2, b18));
#endif
    _mm_store_sd(&C[(l_n*84)+18], c18_2);
    __m128d c18_3 = _mm_load_sd(&C[(l_n*84)+33]);
    __m128d a18_3 = _mm_load_sd(&A[100]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_3 = _mm_add_sd(c18_3, _mm_mul_sd(a18_3, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_3 = _mm_add_sd(c18_3, _mm_mul_sd(a18_3, b18));
#endif
    _mm_store_sd(&C[(l_n*84)+33], c18_3);
    __m128d c18_4 = _mm_load_sd(&C[(l_n*84)+54]);
    __m128d a18_4 = _mm_load_sd(&A[101]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_4 = _mm_add_sd(c18_4, _mm_mul_sd(a18_4, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_4 = _mm_add_sd(c18_4, _mm_mul_sd(a18_4, b18));
#endif
    _mm_store_sd(&C[(l_n*84)+54], c18_4);
    __m128d c18_5 = _mm_load_sd(&C[(l_n*84)+82]);
    __m128d a18_5 = _mm_load_sd(&A[102]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_5 = _mm_add_sd(c18_5, _mm_mul_sd(a18_5, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_5 = _mm_add_sd(c18_5, _mm_mul_sd(a18_5, b18));
#endif
    _mm_store_sd(&C[(l_n*84)+82], c18_5);
#else
    C[(l_n*84)+2] += A[97] * B[(l_n*84)+18];
    C[(l_n*84)+8] += A[98] * B[(l_n*84)+18];
    C[(l_n*84)+18] += A[99] * B[(l_n*84)+18];
    C[(l_n*84)+33] += A[100] * B[(l_n*84)+18];
    C[(l_n*84)+54] += A[101] * B[(l_n*84)+18];
    C[(l_n*84)+82] += A[102] * B[(l_n*84)+18];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b19 = _mm256_broadcast_sd(&B[(l_n*84)+19]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b19 = _mm_loaddup_pd(&B[(l_n*84)+19]);
#endif
    __m128d c19_0 = _mm_load_sd(&C[(l_n*84)+0]);
    __m128d a19_0 = _mm_load_sd(&A[103]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_0 = _mm_add_sd(c19_0, _mm_mul_sd(a19_0, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_0 = _mm_add_sd(c19_0, _mm_mul_sd(a19_0, b19));
#endif
    _mm_store_sd(&C[(l_n*84)+0], c19_0);
    __m128d c19_1 = _mm_load_sd(&C[(l_n*84)+3]);
    __m128d a19_1 = _mm_load_sd(&A[104]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_1 = _mm_add_sd(c19_1, _mm_mul_sd(a19_1, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_1 = _mm_add_sd(c19_1, _mm_mul_sd(a19_1, b19));
#endif
    _mm_store_sd(&C[(l_n*84)+3], c19_1);
    __m128d c19_2 = _mm_load_sd(&C[(l_n*84)+9]);
    __m128d a19_2 = _mm_load_sd(&A[105]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_2 = _mm_add_sd(c19_2, _mm_mul_sd(a19_2, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_2 = _mm_add_sd(c19_2, _mm_mul_sd(a19_2, b19));
#endif
    _mm_store_sd(&C[(l_n*84)+9], c19_2);
    __m128d c19_3 = _mm_load_sd(&C[(l_n*84)+19]);
    __m128d a19_3 = _mm_load_sd(&A[106]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_3 = _mm_add_sd(c19_3, _mm_mul_sd(a19_3, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_3 = _mm_add_sd(c19_3, _mm_mul_sd(a19_3, b19));
#endif
    _mm_store_sd(&C[(l_n*84)+19], c19_3);
    __m128d c19_4 = _mm_load_sd(&C[(l_n*84)+34]);
    __m128d a19_4 = _mm_load_sd(&A[107]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_4 = _mm_add_sd(c19_4, _mm_mul_sd(a19_4, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_4 = _mm_add_sd(c19_4, _mm_mul_sd(a19_4, b19));
#endif
    _mm_store_sd(&C[(l_n*84)+34], c19_4);
    __m128d c19_5 = _mm_load_sd(&C[(l_n*84)+55]);
    __m128d a19_5 = _mm_load_sd(&A[108]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_5 = _mm_add_sd(c19_5, _mm_mul_sd(a19_5, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_5 = _mm_add_sd(c19_5, _mm_mul_sd(a19_5, b19));
#endif
    _mm_store_sd(&C[(l_n*84)+55], c19_5);
    __m128d c19_6 = _mm_load_sd(&C[(l_n*84)+83]);
    __m128d a19_6 = _mm_load_sd(&A[109]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_6 = _mm_add_sd(c19_6, _mm_mul_sd(a19_6, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_6 = _mm_add_sd(c19_6, _mm_mul_sd(a19_6, b19));
#endif
    _mm_store_sd(&C[(l_n*84)+83], c19_6);
#else
    C[(l_n*84)+0] += A[103] * B[(l_n*84)+19];
    C[(l_n*84)+3] += A[104] * B[(l_n*84)+19];
    C[(l_n*84)+9] += A[105] * B[(l_n*84)+19];
    C[(l_n*84)+19] += A[106] * B[(l_n*84)+19];
    C[(l_n*84)+34] += A[107] * B[(l_n*84)+19];
    C[(l_n*84)+55] += A[108] * B[(l_n*84)+19];
    C[(l_n*84)+83] += A[109] * B[(l_n*84)+19];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b20 = _mm256_broadcast_sd(&B[(l_n*84)+20]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b20 = _mm_loaddup_pd(&B[(l_n*84)+20]);
#endif
    __m128d c20_0 = _mm_load_sd(&C[(l_n*84)+20]);
    __m128d a20_0 = _mm_load_sd(&A[110]);
#if defined(__SSE3__) && defined(__AVX__)
    c20_0 = _mm_add_sd(c20_0, _mm_mul_sd(a20_0, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c20_0 = _mm_add_sd(c20_0, _mm_mul_sd(a20_0, b20));
#endif
    _mm_store_sd(&C[(l_n*84)+20], c20_0);
    __m128d c20_1 = _mm_load_sd(&C[(l_n*84)+41]);
    __m128d a20_1 = _mm_load_sd(&A[111]);
#if defined(__SSE3__) && defined(__AVX__)
    c20_1 = _mm_add_sd(c20_1, _mm_mul_sd(a20_1, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c20_1 = _mm_add_sd(c20_1, _mm_mul_sd(a20_1, b20));
#endif
    _mm_store_sd(&C[(l_n*84)+41], c20_1);
    __m128d c20_2 = _mm_load_sd(&C[(l_n*84)+69]);
    __m128d a20_2 = _mm_load_sd(&A[112]);
#if defined(__SSE3__) && defined(__AVX__)
    c20_2 = _mm_add_sd(c20_2, _mm_mul_sd(a20_2, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c20_2 = _mm_add_sd(c20_2, _mm_mul_sd(a20_2, b20));
#endif
    _mm_store_sd(&C[(l_n*84)+69], c20_2);
#else
    C[(l_n*84)+20] += A[110] * B[(l_n*84)+20];
    C[(l_n*84)+41] += A[111] * B[(l_n*84)+20];
    C[(l_n*84)+69] += A[112] * B[(l_n*84)+20];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b21 = _mm256_broadcast_sd(&B[(l_n*84)+21]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b21 = _mm_loaddup_pd(&B[(l_n*84)+21]);
#endif
    __m128d c21_0 = _mm_load_sd(&C[(l_n*84)+21]);
    __m128d a21_0 = _mm_load_sd(&A[113]);
#if defined(__SSE3__) && defined(__AVX__)
    c21_0 = _mm_add_sd(c21_0, _mm_mul_sd(a21_0, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c21_0 = _mm_add_sd(c21_0, _mm_mul_sd(a21_0, b21));
#endif
    _mm_store_sd(&C[(l_n*84)+21], c21_0);
    __m128d c21_1 = _mm_load_sd(&C[(l_n*84)+42]);
    __m128d a21_1 = _mm_load_sd(&A[114]);
#if defined(__SSE3__) && defined(__AVX__)
    c21_1 = _mm_add_sd(c21_1, _mm_mul_sd(a21_1, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c21_1 = _mm_add_sd(c21_1, _mm_mul_sd(a21_1, b21));
#endif
    _mm_store_sd(&C[(l_n*84)+42], c21_1);
    __m128d c21_2 = _mm_load_sd(&C[(l_n*84)+70]);
    __m128d a21_2 = _mm_load_sd(&A[115]);
#if defined(__SSE3__) && defined(__AVX__)
    c21_2 = _mm_add_sd(c21_2, _mm_mul_sd(a21_2, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c21_2 = _mm_add_sd(c21_2, _mm_mul_sd(a21_2, b21));
#endif
    _mm_store_sd(&C[(l_n*84)+70], c21_2);
#else
    C[(l_n*84)+21] += A[113] * B[(l_n*84)+21];
    C[(l_n*84)+42] += A[114] * B[(l_n*84)+21];
    C[(l_n*84)+70] += A[115] * B[(l_n*84)+21];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b22 = _mm256_broadcast_sd(&B[(l_n*84)+22]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b22 = _mm_loaddup_pd(&B[(l_n*84)+22]);
#endif
    __m128d c22_0 = _mm_load_sd(&C[(l_n*84)+22]);
    __m128d a22_0 = _mm_load_sd(&A[116]);
#if defined(__SSE3__) && defined(__AVX__)
    c22_0 = _mm_add_sd(c22_0, _mm_mul_sd(a22_0, _mm256_castpd256_pd128(b22)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c22_0 = _mm_add_sd(c22_0, _mm_mul_sd(a22_0, b22));
#endif
    _mm_store_sd(&C[(l_n*84)+22], c22_0);
    __m128d c22_1 = _mm_load_sd(&C[(l_n*84)+43]);
    __m128d a22_1 = _mm_load_sd(&A[117]);
#if defined(__SSE3__) && defined(__AVX__)
    c22_1 = _mm_add_sd(c22_1, _mm_mul_sd(a22_1, _mm256_castpd256_pd128(b22)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c22_1 = _mm_add_sd(c22_1, _mm_mul_sd(a22_1, b22));
#endif
    _mm_store_sd(&C[(l_n*84)+43], c22_1);
    __m128d c22_2 = _mm_load_sd(&C[(l_n*84)+71]);
    __m128d a22_2 = _mm_load_sd(&A[118]);
#if defined(__SSE3__) && defined(__AVX__)
    c22_2 = _mm_add_sd(c22_2, _mm_mul_sd(a22_2, _mm256_castpd256_pd128(b22)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c22_2 = _mm_add_sd(c22_2, _mm_mul_sd(a22_2, b22));
#endif
    _mm_store_sd(&C[(l_n*84)+71], c22_2);
#else
    C[(l_n*84)+22] += A[116] * B[(l_n*84)+22];
    C[(l_n*84)+43] += A[117] * B[(l_n*84)+22];
    C[(l_n*84)+71] += A[118] * B[(l_n*84)+22];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b23 = _mm256_broadcast_sd(&B[(l_n*84)+23]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b23 = _mm_loaddup_pd(&B[(l_n*84)+23]);
#endif
    __m128d c23_0 = _mm_load_sd(&C[(l_n*84)+23]);
    __m128d a23_0 = _mm_load_sd(&A[119]);
#if defined(__SSE3__) && defined(__AVX__)
    c23_0 = _mm_add_sd(c23_0, _mm_mul_sd(a23_0, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c23_0 = _mm_add_sd(c23_0, _mm_mul_sd(a23_0, b23));
#endif
    _mm_store_sd(&C[(l_n*84)+23], c23_0);
    __m128d c23_1 = _mm_load_sd(&C[(l_n*84)+44]);
    __m128d a23_1 = _mm_load_sd(&A[120]);
#if defined(__SSE3__) && defined(__AVX__)
    c23_1 = _mm_add_sd(c23_1, _mm_mul_sd(a23_1, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c23_1 = _mm_add_sd(c23_1, _mm_mul_sd(a23_1, b23));
#endif
    _mm_store_sd(&C[(l_n*84)+44], c23_1);
    __m128d c23_2 = _mm_load_sd(&C[(l_n*84)+72]);
    __m128d a23_2 = _mm_load_sd(&A[121]);
#if defined(__SSE3__) && defined(__AVX__)
    c23_2 = _mm_add_sd(c23_2, _mm_mul_sd(a23_2, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c23_2 = _mm_add_sd(c23_2, _mm_mul_sd(a23_2, b23));
#endif
    _mm_store_sd(&C[(l_n*84)+72], c23_2);
#else
    C[(l_n*84)+23] += A[119] * B[(l_n*84)+23];
    C[(l_n*84)+44] += A[120] * B[(l_n*84)+23];
    C[(l_n*84)+72] += A[121] * B[(l_n*84)+23];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b24 = _mm256_broadcast_sd(&B[(l_n*84)+24]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b24 = _mm_loaddup_pd(&B[(l_n*84)+24]);
#endif
    __m128d c24_0 = _mm_load_sd(&C[(l_n*84)+24]);
    __m128d a24_0 = _mm_load_sd(&A[122]);
#if defined(__SSE3__) && defined(__AVX__)
    c24_0 = _mm_add_sd(c24_0, _mm_mul_sd(a24_0, _mm256_castpd256_pd128(b24)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c24_0 = _mm_add_sd(c24_0, _mm_mul_sd(a24_0, b24));
#endif
    _mm_store_sd(&C[(l_n*84)+24], c24_0);
    __m128d c24_1 = _mm_load_sd(&C[(l_n*84)+45]);
    __m128d a24_1 = _mm_load_sd(&A[123]);
#if defined(__SSE3__) && defined(__AVX__)
    c24_1 = _mm_add_sd(c24_1, _mm_mul_sd(a24_1, _mm256_castpd256_pd128(b24)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c24_1 = _mm_add_sd(c24_1, _mm_mul_sd(a24_1, b24));
#endif
    _mm_store_sd(&C[(l_n*84)+45], c24_1);
    __m128d c24_2 = _mm_load_sd(&C[(l_n*84)+73]);
    __m128d a24_2 = _mm_load_sd(&A[124]);
#if defined(__SSE3__) && defined(__AVX__)
    c24_2 = _mm_add_sd(c24_2, _mm_mul_sd(a24_2, _mm256_castpd256_pd128(b24)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c24_2 = _mm_add_sd(c24_2, _mm_mul_sd(a24_2, b24));
#endif
    _mm_store_sd(&C[(l_n*84)+73], c24_2);
#else
    C[(l_n*84)+24] += A[122] * B[(l_n*84)+24];
    C[(l_n*84)+45] += A[123] * B[(l_n*84)+24];
    C[(l_n*84)+73] += A[124] * B[(l_n*84)+24];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b25 = _mm256_broadcast_sd(&B[(l_n*84)+25]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b25 = _mm_loaddup_pd(&B[(l_n*84)+25]);
#endif
    __m128d c25_0 = _mm_load_sd(&C[(l_n*84)+10]);
    __m128d a25_0 = _mm_load_sd(&A[125]);
#if defined(__SSE3__) && defined(__AVX__)
    c25_0 = _mm_add_sd(c25_0, _mm_mul_sd(a25_0, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c25_0 = _mm_add_sd(c25_0, _mm_mul_sd(a25_0, b25));
#endif
    _mm_store_sd(&C[(l_n*84)+10], c25_0);
    __m128d c25_1 = _mm_load_sd(&C[(l_n*84)+25]);
    __m128d a25_1 = _mm_load_sd(&A[126]);
#if defined(__SSE3__) && defined(__AVX__)
    c25_1 = _mm_add_sd(c25_1, _mm_mul_sd(a25_1, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c25_1 = _mm_add_sd(c25_1, _mm_mul_sd(a25_1, b25));
#endif
    _mm_store_sd(&C[(l_n*84)+25], c25_1);
    __m128d c25_2 = _mm_load_sd(&C[(l_n*84)+46]);
    __m128d a25_2 = _mm_load_sd(&A[127]);
#if defined(__SSE3__) && defined(__AVX__)
    c25_2 = _mm_add_sd(c25_2, _mm_mul_sd(a25_2, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c25_2 = _mm_add_sd(c25_2, _mm_mul_sd(a25_2, b25));
#endif
    _mm_store_sd(&C[(l_n*84)+46], c25_2);
    __m128d c25_3 = _mm_load_sd(&C[(l_n*84)+74]);
    __m128d a25_3 = _mm_load_sd(&A[128]);
#if defined(__SSE3__) && defined(__AVX__)
    c25_3 = _mm_add_sd(c25_3, _mm_mul_sd(a25_3, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c25_3 = _mm_add_sd(c25_3, _mm_mul_sd(a25_3, b25));
#endif
    _mm_store_sd(&C[(l_n*84)+74], c25_3);
#else
    C[(l_n*84)+10] += A[125] * B[(l_n*84)+25];
    C[(l_n*84)+25] += A[126] * B[(l_n*84)+25];
    C[(l_n*84)+46] += A[127] * B[(l_n*84)+25];
    C[(l_n*84)+74] += A[128] * B[(l_n*84)+25];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b26 = _mm256_broadcast_sd(&B[(l_n*84)+26]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b26 = _mm_loaddup_pd(&B[(l_n*84)+26]);
#endif
    __m128d c26_0 = _mm_load_sd(&C[(l_n*84)+11]);
    __m128d a26_0 = _mm_load_sd(&A[129]);
#if defined(__SSE3__) && defined(__AVX__)
    c26_0 = _mm_add_sd(c26_0, _mm_mul_sd(a26_0, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c26_0 = _mm_add_sd(c26_0, _mm_mul_sd(a26_0, b26));
#endif
    _mm_store_sd(&C[(l_n*84)+11], c26_0);
    __m128d c26_1 = _mm_load_sd(&C[(l_n*84)+26]);
    __m128d a26_1 = _mm_load_sd(&A[130]);
#if defined(__SSE3__) && defined(__AVX__)
    c26_1 = _mm_add_sd(c26_1, _mm_mul_sd(a26_1, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c26_1 = _mm_add_sd(c26_1, _mm_mul_sd(a26_1, b26));
#endif
    _mm_store_sd(&C[(l_n*84)+26], c26_1);
    __m128d c26_2 = _mm_load_sd(&C[(l_n*84)+47]);
    __m128d a26_2 = _mm_load_sd(&A[131]);
#if defined(__SSE3__) && defined(__AVX__)
    c26_2 = _mm_add_sd(c26_2, _mm_mul_sd(a26_2, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c26_2 = _mm_add_sd(c26_2, _mm_mul_sd(a26_2, b26));
#endif
    _mm_store_sd(&C[(l_n*84)+47], c26_2);
    __m128d c26_3 = _mm_load_sd(&C[(l_n*84)+75]);
    __m128d a26_3 = _mm_load_sd(&A[132]);
#if defined(__SSE3__) && defined(__AVX__)
    c26_3 = _mm_add_sd(c26_3, _mm_mul_sd(a26_3, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c26_3 = _mm_add_sd(c26_3, _mm_mul_sd(a26_3, b26));
#endif
    _mm_store_sd(&C[(l_n*84)+75], c26_3);
#else
    C[(l_n*84)+11] += A[129] * B[(l_n*84)+26];
    C[(l_n*84)+26] += A[130] * B[(l_n*84)+26];
    C[(l_n*84)+47] += A[131] * B[(l_n*84)+26];
    C[(l_n*84)+75] += A[132] * B[(l_n*84)+26];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b27 = _mm256_broadcast_sd(&B[(l_n*84)+27]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b27 = _mm_loaddup_pd(&B[(l_n*84)+27]);
#endif
    __m128d c27_0 = _mm_load_sd(&C[(l_n*84)+12]);
    __m128d a27_0 = _mm_load_sd(&A[133]);
#if defined(__SSE3__) && defined(__AVX__)
    c27_0 = _mm_add_sd(c27_0, _mm_mul_sd(a27_0, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c27_0 = _mm_add_sd(c27_0, _mm_mul_sd(a27_0, b27));
#endif
    _mm_store_sd(&C[(l_n*84)+12], c27_0);
    __m128d c27_1 = _mm_load_sd(&C[(l_n*84)+27]);
    __m128d a27_1 = _mm_load_sd(&A[134]);
#if defined(__SSE3__) && defined(__AVX__)
    c27_1 = _mm_add_sd(c27_1, _mm_mul_sd(a27_1, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c27_1 = _mm_add_sd(c27_1, _mm_mul_sd(a27_1, b27));
#endif
    _mm_store_sd(&C[(l_n*84)+27], c27_1);
    __m128d c27_2 = _mm_load_sd(&C[(l_n*84)+48]);
    __m128d a27_2 = _mm_load_sd(&A[135]);
#if defined(__SSE3__) && defined(__AVX__)
    c27_2 = _mm_add_sd(c27_2, _mm_mul_sd(a27_2, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c27_2 = _mm_add_sd(c27_2, _mm_mul_sd(a27_2, b27));
#endif
    _mm_store_sd(&C[(l_n*84)+48], c27_2);
    __m128d c27_3 = _mm_load_sd(&C[(l_n*84)+76]);
    __m128d a27_3 = _mm_load_sd(&A[136]);
#if defined(__SSE3__) && defined(__AVX__)
    c27_3 = _mm_add_sd(c27_3, _mm_mul_sd(a27_3, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c27_3 = _mm_add_sd(c27_3, _mm_mul_sd(a27_3, b27));
#endif
    _mm_store_sd(&C[(l_n*84)+76], c27_3);
#else
    C[(l_n*84)+12] += A[133] * B[(l_n*84)+27];
    C[(l_n*84)+27] += A[134] * B[(l_n*84)+27];
    C[(l_n*84)+48] += A[135] * B[(l_n*84)+27];
    C[(l_n*84)+76] += A[136] * B[(l_n*84)+27];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b28 = _mm256_broadcast_sd(&B[(l_n*84)+28]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b28 = _mm_loaddup_pd(&B[(l_n*84)+28]);
#endif
    __m128d c28_0 = _mm_load_sd(&C[(l_n*84)+13]);
    __m128d a28_0 = _mm_load_sd(&A[137]);
#if defined(__SSE3__) && defined(__AVX__)
    c28_0 = _mm_add_sd(c28_0, _mm_mul_sd(a28_0, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c28_0 = _mm_add_sd(c28_0, _mm_mul_sd(a28_0, b28));
#endif
    _mm_store_sd(&C[(l_n*84)+13], c28_0);
    __m128d c28_1 = _mm_load_sd(&C[(l_n*84)+28]);
    __m128d a28_1 = _mm_load_sd(&A[138]);
#if defined(__SSE3__) && defined(__AVX__)
    c28_1 = _mm_add_sd(c28_1, _mm_mul_sd(a28_1, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c28_1 = _mm_add_sd(c28_1, _mm_mul_sd(a28_1, b28));
#endif
    _mm_store_sd(&C[(l_n*84)+28], c28_1);
    __m128d c28_2 = _mm_load_sd(&C[(l_n*84)+49]);
    __m128d a28_2 = _mm_load_sd(&A[139]);
#if defined(__SSE3__) && defined(__AVX__)
    c28_2 = _mm_add_sd(c28_2, _mm_mul_sd(a28_2, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c28_2 = _mm_add_sd(c28_2, _mm_mul_sd(a28_2, b28));
#endif
    _mm_store_sd(&C[(l_n*84)+49], c28_2);
    __m128d c28_3 = _mm_load_sd(&C[(l_n*84)+77]);
    __m128d a28_3 = _mm_load_sd(&A[140]);
#if defined(__SSE3__) && defined(__AVX__)
    c28_3 = _mm_add_sd(c28_3, _mm_mul_sd(a28_3, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c28_3 = _mm_add_sd(c28_3, _mm_mul_sd(a28_3, b28));
#endif
    _mm_store_sd(&C[(l_n*84)+77], c28_3);
#else
    C[(l_n*84)+13] += A[137] * B[(l_n*84)+28];
    C[(l_n*84)+28] += A[138] * B[(l_n*84)+28];
    C[(l_n*84)+49] += A[139] * B[(l_n*84)+28];
    C[(l_n*84)+77] += A[140] * B[(l_n*84)+28];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b29 = _mm256_broadcast_sd(&B[(l_n*84)+29]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b29 = _mm_loaddup_pd(&B[(l_n*84)+29]);
#endif
    __m128d c29_0 = _mm_load_sd(&C[(l_n*84)+4]);
    __m128d a29_0 = _mm_load_sd(&A[141]);
#if defined(__SSE3__) && defined(__AVX__)
    c29_0 = _mm_add_sd(c29_0, _mm_mul_sd(a29_0, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c29_0 = _mm_add_sd(c29_0, _mm_mul_sd(a29_0, b29));
#endif
    _mm_store_sd(&C[(l_n*84)+4], c29_0);
    __m128d c29_1 = _mm_load_sd(&C[(l_n*84)+14]);
    __m128d a29_1 = _mm_load_sd(&A[142]);
#if defined(__SSE3__) && defined(__AVX__)
    c29_1 = _mm_add_sd(c29_1, _mm_mul_sd(a29_1, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c29_1 = _mm_add_sd(c29_1, _mm_mul_sd(a29_1, b29));
#endif
    _mm_store_sd(&C[(l_n*84)+14], c29_1);
    __m128d c29_2 = _mm_load_sd(&C[(l_n*84)+29]);
    __m128d a29_2 = _mm_load_sd(&A[143]);
#if defined(__SSE3__) && defined(__AVX__)
    c29_2 = _mm_add_sd(c29_2, _mm_mul_sd(a29_2, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c29_2 = _mm_add_sd(c29_2, _mm_mul_sd(a29_2, b29));
#endif
    _mm_store_sd(&C[(l_n*84)+29], c29_2);
    __m128d c29_3 = _mm_load_sd(&C[(l_n*84)+50]);
    __m128d a29_3 = _mm_load_sd(&A[144]);
#if defined(__SSE3__) && defined(__AVX__)
    c29_3 = _mm_add_sd(c29_3, _mm_mul_sd(a29_3, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c29_3 = _mm_add_sd(c29_3, _mm_mul_sd(a29_3, b29));
#endif
    _mm_store_sd(&C[(l_n*84)+50], c29_3);
    __m128d c29_4 = _mm_load_sd(&C[(l_n*84)+78]);
    __m128d a29_4 = _mm_load_sd(&A[145]);
#if defined(__SSE3__) && defined(__AVX__)
    c29_4 = _mm_add_sd(c29_4, _mm_mul_sd(a29_4, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c29_4 = _mm_add_sd(c29_4, _mm_mul_sd(a29_4, b29));
#endif
    _mm_store_sd(&C[(l_n*84)+78], c29_4);
#else
    C[(l_n*84)+4] += A[141] * B[(l_n*84)+29];
    C[(l_n*84)+14] += A[142] * B[(l_n*84)+29];
    C[(l_n*84)+29] += A[143] * B[(l_n*84)+29];
    C[(l_n*84)+50] += A[144] * B[(l_n*84)+29];
    C[(l_n*84)+78] += A[145] * B[(l_n*84)+29];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b30 = _mm256_broadcast_sd(&B[(l_n*84)+30]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b30 = _mm_loaddup_pd(&B[(l_n*84)+30]);
#endif
    __m128d c30_0 = _mm_load_sd(&C[(l_n*84)+5]);
    __m128d a30_0 = _mm_load_sd(&A[146]);
#if defined(__SSE3__) && defined(__AVX__)
    c30_0 = _mm_add_sd(c30_0, _mm_mul_sd(a30_0, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c30_0 = _mm_add_sd(c30_0, _mm_mul_sd(a30_0, b30));
#endif
    _mm_store_sd(&C[(l_n*84)+5], c30_0);
    __m128d c30_1 = _mm_load_sd(&C[(l_n*84)+15]);
    __m128d a30_1 = _mm_load_sd(&A[147]);
#if defined(__SSE3__) && defined(__AVX__)
    c30_1 = _mm_add_sd(c30_1, _mm_mul_sd(a30_1, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c30_1 = _mm_add_sd(c30_1, _mm_mul_sd(a30_1, b30));
#endif
    _mm_store_sd(&C[(l_n*84)+15], c30_1);
    __m128d c30_2 = _mm_load_sd(&C[(l_n*84)+30]);
    __m128d a30_2 = _mm_load_sd(&A[148]);
#if defined(__SSE3__) && defined(__AVX__)
    c30_2 = _mm_add_sd(c30_2, _mm_mul_sd(a30_2, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c30_2 = _mm_add_sd(c30_2, _mm_mul_sd(a30_2, b30));
#endif
    _mm_store_sd(&C[(l_n*84)+30], c30_2);
    __m128d c30_3 = _mm_load_sd(&C[(l_n*84)+51]);
    __m128d a30_3 = _mm_load_sd(&A[149]);
#if defined(__SSE3__) && defined(__AVX__)
    c30_3 = _mm_add_sd(c30_3, _mm_mul_sd(a30_3, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c30_3 = _mm_add_sd(c30_3, _mm_mul_sd(a30_3, b30));
#endif
    _mm_store_sd(&C[(l_n*84)+51], c30_3);
    __m128d c30_4 = _mm_load_sd(&C[(l_n*84)+79]);
    __m128d a30_4 = _mm_load_sd(&A[150]);
#if defined(__SSE3__) && defined(__AVX__)
    c30_4 = _mm_add_sd(c30_4, _mm_mul_sd(a30_4, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c30_4 = _mm_add_sd(c30_4, _mm_mul_sd(a30_4, b30));
#endif
    _mm_store_sd(&C[(l_n*84)+79], c30_4);
#else
    C[(l_n*84)+5] += A[146] * B[(l_n*84)+30];
    C[(l_n*84)+15] += A[147] * B[(l_n*84)+30];
    C[(l_n*84)+30] += A[148] * B[(l_n*84)+30];
    C[(l_n*84)+51] += A[149] * B[(l_n*84)+30];
    C[(l_n*84)+79] += A[150] * B[(l_n*84)+30];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b31 = _mm256_broadcast_sd(&B[(l_n*84)+31]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b31 = _mm_loaddup_pd(&B[(l_n*84)+31]);
#endif
    __m128d c31_0 = _mm_load_sd(&C[(l_n*84)+6]);
    __m128d a31_0 = _mm_load_sd(&A[151]);
#if defined(__SSE3__) && defined(__AVX__)
    c31_0 = _mm_add_sd(c31_0, _mm_mul_sd(a31_0, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c31_0 = _mm_add_sd(c31_0, _mm_mul_sd(a31_0, b31));
#endif
    _mm_store_sd(&C[(l_n*84)+6], c31_0);
    __m128d c31_1 = _mm_load_sd(&C[(l_n*84)+16]);
    __m128d a31_1 = _mm_load_sd(&A[152]);
#if defined(__SSE3__) && defined(__AVX__)
    c31_1 = _mm_add_sd(c31_1, _mm_mul_sd(a31_1, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c31_1 = _mm_add_sd(c31_1, _mm_mul_sd(a31_1, b31));
#endif
    _mm_store_sd(&C[(l_n*84)+16], c31_1);
    __m128d c31_2 = _mm_load_sd(&C[(l_n*84)+31]);
    __m128d a31_2 = _mm_load_sd(&A[153]);
#if defined(__SSE3__) && defined(__AVX__)
    c31_2 = _mm_add_sd(c31_2, _mm_mul_sd(a31_2, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c31_2 = _mm_add_sd(c31_2, _mm_mul_sd(a31_2, b31));
#endif
    _mm_store_sd(&C[(l_n*84)+31], c31_2);
    __m128d c31_3 = _mm_load_sd(&C[(l_n*84)+52]);
    __m128d a31_3 = _mm_load_sd(&A[154]);
#if defined(__SSE3__) && defined(__AVX__)
    c31_3 = _mm_add_sd(c31_3, _mm_mul_sd(a31_3, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c31_3 = _mm_add_sd(c31_3, _mm_mul_sd(a31_3, b31));
#endif
    _mm_store_sd(&C[(l_n*84)+52], c31_3);
    __m128d c31_4 = _mm_load_sd(&C[(l_n*84)+80]);
    __m128d a31_4 = _mm_load_sd(&A[155]);
#if defined(__SSE3__) && defined(__AVX__)
    c31_4 = _mm_add_sd(c31_4, _mm_mul_sd(a31_4, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c31_4 = _mm_add_sd(c31_4, _mm_mul_sd(a31_4, b31));
#endif
    _mm_store_sd(&C[(l_n*84)+80], c31_4);
#else
    C[(l_n*84)+6] += A[151] * B[(l_n*84)+31];
    C[(l_n*84)+16] += A[152] * B[(l_n*84)+31];
    C[(l_n*84)+31] += A[153] * B[(l_n*84)+31];
    C[(l_n*84)+52] += A[154] * B[(l_n*84)+31];
    C[(l_n*84)+80] += A[155] * B[(l_n*84)+31];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b32 = _mm256_broadcast_sd(&B[(l_n*84)+32]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b32 = _mm_loaddup_pd(&B[(l_n*84)+32]);
#endif
    __m128d c32_0 = _mm_load_sd(&C[(l_n*84)+1]);
    __m128d a32_0 = _mm_load_sd(&A[156]);
#if defined(__SSE3__) && defined(__AVX__)
    c32_0 = _mm_add_sd(c32_0, _mm_mul_sd(a32_0, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c32_0 = _mm_add_sd(c32_0, _mm_mul_sd(a32_0, b32));
#endif
    _mm_store_sd(&C[(l_n*84)+1], c32_0);
    __m128d c32_1 = _mm_load_sd(&C[(l_n*84)+7]);
    __m128d a32_1 = _mm_load_sd(&A[157]);
#if defined(__SSE3__) && defined(__AVX__)
    c32_1 = _mm_add_sd(c32_1, _mm_mul_sd(a32_1, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c32_1 = _mm_add_sd(c32_1, _mm_mul_sd(a32_1, b32));
#endif
    _mm_store_sd(&C[(l_n*84)+7], c32_1);
    __m128d c32_2 = _mm_load_sd(&C[(l_n*84)+17]);
    __m128d a32_2 = _mm_load_sd(&A[158]);
#if defined(__SSE3__) && defined(__AVX__)
    c32_2 = _mm_add_sd(c32_2, _mm_mul_sd(a32_2, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c32_2 = _mm_add_sd(c32_2, _mm_mul_sd(a32_2, b32));
#endif
    _mm_store_sd(&C[(l_n*84)+17], c32_2);
    __m128d c32_3 = _mm_load_sd(&C[(l_n*84)+32]);
    __m128d a32_3 = _mm_load_sd(&A[159]);
#if defined(__SSE3__) && defined(__AVX__)
    c32_3 = _mm_add_sd(c32_3, _mm_mul_sd(a32_3, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c32_3 = _mm_add_sd(c32_3, _mm_mul_sd(a32_3, b32));
#endif
    _mm_store_sd(&C[(l_n*84)+32], c32_3);
    __m128d c32_4 = _mm_load_sd(&C[(l_n*84)+53]);
    __m128d a32_4 = _mm_load_sd(&A[160]);
#if defined(__SSE3__) && defined(__AVX__)
    c32_4 = _mm_add_sd(c32_4, _mm_mul_sd(a32_4, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c32_4 = _mm_add_sd(c32_4, _mm_mul_sd(a32_4, b32));
#endif
    _mm_store_sd(&C[(l_n*84)+53], c32_4);
    __m128d c32_5 = _mm_load_sd(&C[(l_n*84)+81]);
    __m128d a32_5 = _mm_load_sd(&A[161]);
#if defined(__SSE3__) && defined(__AVX__)
    c32_5 = _mm_add_sd(c32_5, _mm_mul_sd(a32_5, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c32_5 = _mm_add_sd(c32_5, _mm_mul_sd(a32_5, b32));
#endif
    _mm_store_sd(&C[(l_n*84)+81], c32_5);
#else
    C[(l_n*84)+1] += A[156] * B[(l_n*84)+32];
    C[(l_n*84)+7] += A[157] * B[(l_n*84)+32];
    C[(l_n*84)+17] += A[158] * B[(l_n*84)+32];
    C[(l_n*84)+32] += A[159] * B[(l_n*84)+32];
    C[(l_n*84)+53] += A[160] * B[(l_n*84)+32];
    C[(l_n*84)+81] += A[161] * B[(l_n*84)+32];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b33 = _mm256_broadcast_sd(&B[(l_n*84)+33]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b33 = _mm_loaddup_pd(&B[(l_n*84)+33]);
#endif
    __m128d c33_0 = _mm_load_sd(&C[(l_n*84)+2]);
    __m128d a33_0 = _mm_load_sd(&A[162]);
#if defined(__SSE3__) && defined(__AVX__)
    c33_0 = _mm_add_sd(c33_0, _mm_mul_sd(a33_0, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c33_0 = _mm_add_sd(c33_0, _mm_mul_sd(a33_0, b33));
#endif
    _mm_store_sd(&C[(l_n*84)+2], c33_0);
    __m128d c33_1 = _mm_load_sd(&C[(l_n*84)+8]);
    __m128d a33_1 = _mm_load_sd(&A[163]);
#if defined(__SSE3__) && defined(__AVX__)
    c33_1 = _mm_add_sd(c33_1, _mm_mul_sd(a33_1, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c33_1 = _mm_add_sd(c33_1, _mm_mul_sd(a33_1, b33));
#endif
    _mm_store_sd(&C[(l_n*84)+8], c33_1);
    __m128d c33_2 = _mm_load_sd(&C[(l_n*84)+18]);
    __m128d a33_2 = _mm_load_sd(&A[164]);
#if defined(__SSE3__) && defined(__AVX__)
    c33_2 = _mm_add_sd(c33_2, _mm_mul_sd(a33_2, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c33_2 = _mm_add_sd(c33_2, _mm_mul_sd(a33_2, b33));
#endif
    _mm_store_sd(&C[(l_n*84)+18], c33_2);
    __m128d c33_3 = _mm_load_sd(&C[(l_n*84)+33]);
    __m128d a33_3 = _mm_load_sd(&A[165]);
#if defined(__SSE3__) && defined(__AVX__)
    c33_3 = _mm_add_sd(c33_3, _mm_mul_sd(a33_3, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c33_3 = _mm_add_sd(c33_3, _mm_mul_sd(a33_3, b33));
#endif
    _mm_store_sd(&C[(l_n*84)+33], c33_3);
    __m128d c33_4 = _mm_load_sd(&C[(l_n*84)+54]);
    __m128d a33_4 = _mm_load_sd(&A[166]);
#if defined(__SSE3__) && defined(__AVX__)
    c33_4 = _mm_add_sd(c33_4, _mm_mul_sd(a33_4, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c33_4 = _mm_add_sd(c33_4, _mm_mul_sd(a33_4, b33));
#endif
    _mm_store_sd(&C[(l_n*84)+54], c33_4);
    __m128d c33_5 = _mm_load_sd(&C[(l_n*84)+82]);
    __m128d a33_5 = _mm_load_sd(&A[167]);
#if defined(__SSE3__) && defined(__AVX__)
    c33_5 = _mm_add_sd(c33_5, _mm_mul_sd(a33_5, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c33_5 = _mm_add_sd(c33_5, _mm_mul_sd(a33_5, b33));
#endif
    _mm_store_sd(&C[(l_n*84)+82], c33_5);
#else
    C[(l_n*84)+2] += A[162] * B[(l_n*84)+33];
    C[(l_n*84)+8] += A[163] * B[(l_n*84)+33];
    C[(l_n*84)+18] += A[164] * B[(l_n*84)+33];
    C[(l_n*84)+33] += A[165] * B[(l_n*84)+33];
    C[(l_n*84)+54] += A[166] * B[(l_n*84)+33];
    C[(l_n*84)+82] += A[167] * B[(l_n*84)+33];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b34 = _mm256_broadcast_sd(&B[(l_n*84)+34]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b34 = _mm_loaddup_pd(&B[(l_n*84)+34]);
#endif
    __m128d c34_0 = _mm_load_sd(&C[(l_n*84)+0]);
    __m128d a34_0 = _mm_load_sd(&A[168]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_0 = _mm_add_sd(c34_0, _mm_mul_sd(a34_0, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_0 = _mm_add_sd(c34_0, _mm_mul_sd(a34_0, b34));
#endif
    _mm_store_sd(&C[(l_n*84)+0], c34_0);
    __m128d c34_1 = _mm_load_sd(&C[(l_n*84)+3]);
    __m128d a34_1 = _mm_load_sd(&A[169]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_1 = _mm_add_sd(c34_1, _mm_mul_sd(a34_1, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_1 = _mm_add_sd(c34_1, _mm_mul_sd(a34_1, b34));
#endif
    _mm_store_sd(&C[(l_n*84)+3], c34_1);
    __m128d c34_2 = _mm_load_sd(&C[(l_n*84)+9]);
    __m128d a34_2 = _mm_load_sd(&A[170]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_2 = _mm_add_sd(c34_2, _mm_mul_sd(a34_2, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_2 = _mm_add_sd(c34_2, _mm_mul_sd(a34_2, b34));
#endif
    _mm_store_sd(&C[(l_n*84)+9], c34_2);
    __m128d c34_3 = _mm_load_sd(&C[(l_n*84)+19]);
    __m128d a34_3 = _mm_load_sd(&A[171]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_3 = _mm_add_sd(c34_3, _mm_mul_sd(a34_3, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_3 = _mm_add_sd(c34_3, _mm_mul_sd(a34_3, b34));
#endif
    _mm_store_sd(&C[(l_n*84)+19], c34_3);
    __m128d c34_4 = _mm_load_sd(&C[(l_n*84)+34]);
    __m128d a34_4 = _mm_load_sd(&A[172]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_4 = _mm_add_sd(c34_4, _mm_mul_sd(a34_4, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_4 = _mm_add_sd(c34_4, _mm_mul_sd(a34_4, b34));
#endif
    _mm_store_sd(&C[(l_n*84)+34], c34_4);
    __m128d c34_5 = _mm_load_sd(&C[(l_n*84)+55]);
    __m128d a34_5 = _mm_load_sd(&A[173]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_5 = _mm_add_sd(c34_5, _mm_mul_sd(a34_5, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_5 = _mm_add_sd(c34_5, _mm_mul_sd(a34_5, b34));
#endif
    _mm_store_sd(&C[(l_n*84)+55], c34_5);
    __m128d c34_6 = _mm_load_sd(&C[(l_n*84)+83]);
    __m128d a34_6 = _mm_load_sd(&A[174]);
#if defined(__SSE3__) && defined(__AVX__)
    c34_6 = _mm_add_sd(c34_6, _mm_mul_sd(a34_6, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c34_6 = _mm_add_sd(c34_6, _mm_mul_sd(a34_6, b34));
#endif
    _mm_store_sd(&C[(l_n*84)+83], c34_6);
#else
    C[(l_n*84)+0] += A[168] * B[(l_n*84)+34];
    C[(l_n*84)+3] += A[169] * B[(l_n*84)+34];
    C[(l_n*84)+9] += A[170] * B[(l_n*84)+34];
    C[(l_n*84)+19] += A[171] * B[(l_n*84)+34];
    C[(l_n*84)+34] += A[172] * B[(l_n*84)+34];
    C[(l_n*84)+55] += A[173] * B[(l_n*84)+34];
    C[(l_n*84)+83] += A[174] * B[(l_n*84)+34];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b35 = _mm256_broadcast_sd(&B[(l_n*84)+35]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b35 = _mm_loaddup_pd(&B[(l_n*84)+35]);
#endif
    __m128d c35_0 = _mm_load_sd(&C[(l_n*84)+35]);
    __m128d a35_0 = _mm_load_sd(&A[175]);
#if defined(__SSE3__) && defined(__AVX__)
    c35_0 = _mm_add_sd(c35_0, _mm_mul_sd(a35_0, _mm256_castpd256_pd128(b35)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c35_0 = _mm_add_sd(c35_0, _mm_mul_sd(a35_0, b35));
#endif
    _mm_store_sd(&C[(l_n*84)+35], c35_0);
    __m128d c35_1 = _mm_load_sd(&C[(l_n*84)+63]);
    __m128d a35_1 = _mm_load_sd(&A[176]);
#if defined(__SSE3__) && defined(__AVX__)
    c35_1 = _mm_add_sd(c35_1, _mm_mul_sd(a35_1, _mm256_castpd256_pd128(b35)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c35_1 = _mm_add_sd(c35_1, _mm_mul_sd(a35_1, b35));
#endif
    _mm_store_sd(&C[(l_n*84)+63], c35_1);
#else
    C[(l_n*84)+35] += A[175] * B[(l_n*84)+35];
    C[(l_n*84)+63] += A[176] * B[(l_n*84)+35];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b36 = _mm256_broadcast_sd(&B[(l_n*84)+36]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b36 = _mm_loaddup_pd(&B[(l_n*84)+36]);
#endif
    __m128d c36_0 = _mm_load_sd(&C[(l_n*84)+36]);
    __m128d a36_0 = _mm_load_sd(&A[177]);
#if defined(__SSE3__) && defined(__AVX__)
    c36_0 = _mm_add_sd(c36_0, _mm_mul_sd(a36_0, _mm256_castpd256_pd128(b36)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c36_0 = _mm_add_sd(c36_0, _mm_mul_sd(a36_0, b36));
#endif
    _mm_store_sd(&C[(l_n*84)+36], c36_0);
    __m128d c36_1 = _mm_load_sd(&C[(l_n*84)+64]);
    __m128d a36_1 = _mm_load_sd(&A[178]);
#if defined(__SSE3__) && defined(__AVX__)
    c36_1 = _mm_add_sd(c36_1, _mm_mul_sd(a36_1, _mm256_castpd256_pd128(b36)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c36_1 = _mm_add_sd(c36_1, _mm_mul_sd(a36_1, b36));
#endif
    _mm_store_sd(&C[(l_n*84)+64], c36_1);
#else
    C[(l_n*84)+36] += A[177] * B[(l_n*84)+36];
    C[(l_n*84)+64] += A[178] * B[(l_n*84)+36];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b37 = _mm256_broadcast_sd(&B[(l_n*84)+37]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b37 = _mm_loaddup_pd(&B[(l_n*84)+37]);
#endif
    __m128d c37_0 = _mm_load_sd(&C[(l_n*84)+37]);
    __m128d a37_0 = _mm_load_sd(&A[179]);
#if defined(__SSE3__) && defined(__AVX__)
    c37_0 = _mm_add_sd(c37_0, _mm_mul_sd(a37_0, _mm256_castpd256_pd128(b37)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c37_0 = _mm_add_sd(c37_0, _mm_mul_sd(a37_0, b37));
#endif
    _mm_store_sd(&C[(l_n*84)+37], c37_0);
    __m128d c37_1 = _mm_load_sd(&C[(l_n*84)+65]);
    __m128d a37_1 = _mm_load_sd(&A[180]);
#if defined(__SSE3__) && defined(__AVX__)
    c37_1 = _mm_add_sd(c37_1, _mm_mul_sd(a37_1, _mm256_castpd256_pd128(b37)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c37_1 = _mm_add_sd(c37_1, _mm_mul_sd(a37_1, b37));
#endif
    _mm_store_sd(&C[(l_n*84)+65], c37_1);
#else
    C[(l_n*84)+37] += A[179] * B[(l_n*84)+37];
    C[(l_n*84)+65] += A[180] * B[(l_n*84)+37];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b38 = _mm256_broadcast_sd(&B[(l_n*84)+38]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b38 = _mm_loaddup_pd(&B[(l_n*84)+38]);
#endif
    __m128d c38_0 = _mm_load_sd(&C[(l_n*84)+38]);
    __m128d a38_0 = _mm_load_sd(&A[181]);
#if defined(__SSE3__) && defined(__AVX__)
    c38_0 = _mm_add_sd(c38_0, _mm_mul_sd(a38_0, _mm256_castpd256_pd128(b38)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c38_0 = _mm_add_sd(c38_0, _mm_mul_sd(a38_0, b38));
#endif
    _mm_store_sd(&C[(l_n*84)+38], c38_0);
    __m128d c38_1 = _mm_load_sd(&C[(l_n*84)+66]);
    __m128d a38_1 = _mm_load_sd(&A[182]);
#if defined(__SSE3__) && defined(__AVX__)
    c38_1 = _mm_add_sd(c38_1, _mm_mul_sd(a38_1, _mm256_castpd256_pd128(b38)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c38_1 = _mm_add_sd(c38_1, _mm_mul_sd(a38_1, b38));
#endif
    _mm_store_sd(&C[(l_n*84)+66], c38_1);
#else
    C[(l_n*84)+38] += A[181] * B[(l_n*84)+38];
    C[(l_n*84)+66] += A[182] * B[(l_n*84)+38];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b39 = _mm256_broadcast_sd(&B[(l_n*84)+39]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b39 = _mm_loaddup_pd(&B[(l_n*84)+39]);
#endif
    __m128d c39_0 = _mm_load_sd(&C[(l_n*84)+39]);
    __m128d a39_0 = _mm_load_sd(&A[183]);
#if defined(__SSE3__) && defined(__AVX__)
    c39_0 = _mm_add_sd(c39_0, _mm_mul_sd(a39_0, _mm256_castpd256_pd128(b39)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c39_0 = _mm_add_sd(c39_0, _mm_mul_sd(a39_0, b39));
#endif
    _mm_store_sd(&C[(l_n*84)+39], c39_0);
    __m128d c39_1 = _mm_load_sd(&C[(l_n*84)+67]);
    __m128d a39_1 = _mm_load_sd(&A[184]);
#if defined(__SSE3__) && defined(__AVX__)
    c39_1 = _mm_add_sd(c39_1, _mm_mul_sd(a39_1, _mm256_castpd256_pd128(b39)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c39_1 = _mm_add_sd(c39_1, _mm_mul_sd(a39_1, b39));
#endif
    _mm_store_sd(&C[(l_n*84)+67], c39_1);
#else
    C[(l_n*84)+39] += A[183] * B[(l_n*84)+39];
    C[(l_n*84)+67] += A[184] * B[(l_n*84)+39];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b40 = _mm256_broadcast_sd(&B[(l_n*84)+40]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b40 = _mm_loaddup_pd(&B[(l_n*84)+40]);
#endif
    __m128d c40_0 = _mm_load_sd(&C[(l_n*84)+40]);
    __m128d a40_0 = _mm_load_sd(&A[185]);
#if defined(__SSE3__) && defined(__AVX__)
    c40_0 = _mm_add_sd(c40_0, _mm_mul_sd(a40_0, _mm256_castpd256_pd128(b40)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c40_0 = _mm_add_sd(c40_0, _mm_mul_sd(a40_0, b40));
#endif
    _mm_store_sd(&C[(l_n*84)+40], c40_0);
    __m128d c40_1 = _mm_load_sd(&C[(l_n*84)+68]);
    __m128d a40_1 = _mm_load_sd(&A[186]);
#if defined(__SSE3__) && defined(__AVX__)
    c40_1 = _mm_add_sd(c40_1, _mm_mul_sd(a40_1, _mm256_castpd256_pd128(b40)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c40_1 = _mm_add_sd(c40_1, _mm_mul_sd(a40_1, b40));
#endif
    _mm_store_sd(&C[(l_n*84)+68], c40_1);
#else
    C[(l_n*84)+40] += A[185] * B[(l_n*84)+40];
    C[(l_n*84)+68] += A[186] * B[(l_n*84)+40];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b41 = _mm256_broadcast_sd(&B[(l_n*84)+41]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b41 = _mm_loaddup_pd(&B[(l_n*84)+41]);
#endif
    __m128d c41_0 = _mm_load_sd(&C[(l_n*84)+20]);
    __m128d a41_0 = _mm_load_sd(&A[187]);
#if defined(__SSE3__) && defined(__AVX__)
    c41_0 = _mm_add_sd(c41_0, _mm_mul_sd(a41_0, _mm256_castpd256_pd128(b41)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c41_0 = _mm_add_sd(c41_0, _mm_mul_sd(a41_0, b41));
#endif
    _mm_store_sd(&C[(l_n*84)+20], c41_0);
    __m128d c41_1 = _mm_load_sd(&C[(l_n*84)+41]);
    __m128d a41_1 = _mm_load_sd(&A[188]);
#if defined(__SSE3__) && defined(__AVX__)
    c41_1 = _mm_add_sd(c41_1, _mm_mul_sd(a41_1, _mm256_castpd256_pd128(b41)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c41_1 = _mm_add_sd(c41_1, _mm_mul_sd(a41_1, b41));
#endif
    _mm_store_sd(&C[(l_n*84)+41], c41_1);
    __m128d c41_2 = _mm_load_sd(&C[(l_n*84)+69]);
    __m128d a41_2 = _mm_load_sd(&A[189]);
#if defined(__SSE3__) && defined(__AVX__)
    c41_2 = _mm_add_sd(c41_2, _mm_mul_sd(a41_2, _mm256_castpd256_pd128(b41)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c41_2 = _mm_add_sd(c41_2, _mm_mul_sd(a41_2, b41));
#endif
    _mm_store_sd(&C[(l_n*84)+69], c41_2);
#else
    C[(l_n*84)+20] += A[187] * B[(l_n*84)+41];
    C[(l_n*84)+41] += A[188] * B[(l_n*84)+41];
    C[(l_n*84)+69] += A[189] * B[(l_n*84)+41];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b42 = _mm256_broadcast_sd(&B[(l_n*84)+42]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b42 = _mm_loaddup_pd(&B[(l_n*84)+42]);
#endif
    __m128d c42_0 = _mm_load_sd(&C[(l_n*84)+21]);
    __m128d a42_0 = _mm_load_sd(&A[190]);
#if defined(__SSE3__) && defined(__AVX__)
    c42_0 = _mm_add_sd(c42_0, _mm_mul_sd(a42_0, _mm256_castpd256_pd128(b42)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c42_0 = _mm_add_sd(c42_0, _mm_mul_sd(a42_0, b42));
#endif
    _mm_store_sd(&C[(l_n*84)+21], c42_0);
    __m128d c42_1 = _mm_load_sd(&C[(l_n*84)+42]);
    __m128d a42_1 = _mm_load_sd(&A[191]);
#if defined(__SSE3__) && defined(__AVX__)
    c42_1 = _mm_add_sd(c42_1, _mm_mul_sd(a42_1, _mm256_castpd256_pd128(b42)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c42_1 = _mm_add_sd(c42_1, _mm_mul_sd(a42_1, b42));
#endif
    _mm_store_sd(&C[(l_n*84)+42], c42_1);
    __m128d c42_2 = _mm_load_sd(&C[(l_n*84)+70]);
    __m128d a42_2 = _mm_load_sd(&A[192]);
#if defined(__SSE3__) && defined(__AVX__)
    c42_2 = _mm_add_sd(c42_2, _mm_mul_sd(a42_2, _mm256_castpd256_pd128(b42)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c42_2 = _mm_add_sd(c42_2, _mm_mul_sd(a42_2, b42));
#endif
    _mm_store_sd(&C[(l_n*84)+70], c42_2);
#else
    C[(l_n*84)+21] += A[190] * B[(l_n*84)+42];
    C[(l_n*84)+42] += A[191] * B[(l_n*84)+42];
    C[(l_n*84)+70] += A[192] * B[(l_n*84)+42];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b43 = _mm256_broadcast_sd(&B[(l_n*84)+43]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b43 = _mm_loaddup_pd(&B[(l_n*84)+43]);
#endif
    __m128d c43_0 = _mm_load_sd(&C[(l_n*84)+22]);
    __m128d a43_0 = _mm_load_sd(&A[193]);
#if defined(__SSE3__) && defined(__AVX__)
    c43_0 = _mm_add_sd(c43_0, _mm_mul_sd(a43_0, _mm256_castpd256_pd128(b43)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c43_0 = _mm_add_sd(c43_0, _mm_mul_sd(a43_0, b43));
#endif
    _mm_store_sd(&C[(l_n*84)+22], c43_0);
    __m128d c43_1 = _mm_load_sd(&C[(l_n*84)+43]);
    __m128d a43_1 = _mm_load_sd(&A[194]);
#if defined(__SSE3__) && defined(__AVX__)
    c43_1 = _mm_add_sd(c43_1, _mm_mul_sd(a43_1, _mm256_castpd256_pd128(b43)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c43_1 = _mm_add_sd(c43_1, _mm_mul_sd(a43_1, b43));
#endif
    _mm_store_sd(&C[(l_n*84)+43], c43_1);
    __m128d c43_2 = _mm_load_sd(&C[(l_n*84)+71]);
    __m128d a43_2 = _mm_load_sd(&A[195]);
#if defined(__SSE3__) && defined(__AVX__)
    c43_2 = _mm_add_sd(c43_2, _mm_mul_sd(a43_2, _mm256_castpd256_pd128(b43)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c43_2 = _mm_add_sd(c43_2, _mm_mul_sd(a43_2, b43));
#endif
    _mm_store_sd(&C[(l_n*84)+71], c43_2);
#else
    C[(l_n*84)+22] += A[193] * B[(l_n*84)+43];
    C[(l_n*84)+43] += A[194] * B[(l_n*84)+43];
    C[(l_n*84)+71] += A[195] * B[(l_n*84)+43];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b44 = _mm256_broadcast_sd(&B[(l_n*84)+44]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b44 = _mm_loaddup_pd(&B[(l_n*84)+44]);
#endif
    __m128d c44_0 = _mm_load_sd(&C[(l_n*84)+23]);
    __m128d a44_0 = _mm_load_sd(&A[196]);
#if defined(__SSE3__) && defined(__AVX__)
    c44_0 = _mm_add_sd(c44_0, _mm_mul_sd(a44_0, _mm256_castpd256_pd128(b44)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c44_0 = _mm_add_sd(c44_0, _mm_mul_sd(a44_0, b44));
#endif
    _mm_store_sd(&C[(l_n*84)+23], c44_0);
    __m128d c44_1 = _mm_load_sd(&C[(l_n*84)+44]);
    __m128d a44_1 = _mm_load_sd(&A[197]);
#if defined(__SSE3__) && defined(__AVX__)
    c44_1 = _mm_add_sd(c44_1, _mm_mul_sd(a44_1, _mm256_castpd256_pd128(b44)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c44_1 = _mm_add_sd(c44_1, _mm_mul_sd(a44_1, b44));
#endif
    _mm_store_sd(&C[(l_n*84)+44], c44_1);
    __m128d c44_2 = _mm_load_sd(&C[(l_n*84)+72]);
    __m128d a44_2 = _mm_load_sd(&A[198]);
#if defined(__SSE3__) && defined(__AVX__)
    c44_2 = _mm_add_sd(c44_2, _mm_mul_sd(a44_2, _mm256_castpd256_pd128(b44)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c44_2 = _mm_add_sd(c44_2, _mm_mul_sd(a44_2, b44));
#endif
    _mm_store_sd(&C[(l_n*84)+72], c44_2);
#else
    C[(l_n*84)+23] += A[196] * B[(l_n*84)+44];
    C[(l_n*84)+44] += A[197] * B[(l_n*84)+44];
    C[(l_n*84)+72] += A[198] * B[(l_n*84)+44];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b45 = _mm256_broadcast_sd(&B[(l_n*84)+45]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b45 = _mm_loaddup_pd(&B[(l_n*84)+45]);
#endif
    __m128d c45_0 = _mm_load_sd(&C[(l_n*84)+24]);
    __m128d a45_0 = _mm_load_sd(&A[199]);
#if defined(__SSE3__) && defined(__AVX__)
    c45_0 = _mm_add_sd(c45_0, _mm_mul_sd(a45_0, _mm256_castpd256_pd128(b45)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c45_0 = _mm_add_sd(c45_0, _mm_mul_sd(a45_0, b45));
#endif
    _mm_store_sd(&C[(l_n*84)+24], c45_0);
    __m128d c45_1 = _mm_load_sd(&C[(l_n*84)+45]);
    __m128d a45_1 = _mm_load_sd(&A[200]);
#if defined(__SSE3__) && defined(__AVX__)
    c45_1 = _mm_add_sd(c45_1, _mm_mul_sd(a45_1, _mm256_castpd256_pd128(b45)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c45_1 = _mm_add_sd(c45_1, _mm_mul_sd(a45_1, b45));
#endif
    _mm_store_sd(&C[(l_n*84)+45], c45_1);
    __m128d c45_2 = _mm_load_sd(&C[(l_n*84)+73]);
    __m128d a45_2 = _mm_load_sd(&A[201]);
#if defined(__SSE3__) && defined(__AVX__)
    c45_2 = _mm_add_sd(c45_2, _mm_mul_sd(a45_2, _mm256_castpd256_pd128(b45)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c45_2 = _mm_add_sd(c45_2, _mm_mul_sd(a45_2, b45));
#endif
    _mm_store_sd(&C[(l_n*84)+73], c45_2);
#else
    C[(l_n*84)+24] += A[199] * B[(l_n*84)+45];
    C[(l_n*84)+45] += A[200] * B[(l_n*84)+45];
    C[(l_n*84)+73] += A[201] * B[(l_n*84)+45];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b46 = _mm256_broadcast_sd(&B[(l_n*84)+46]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b46 = _mm_loaddup_pd(&B[(l_n*84)+46]);
#endif
    __m128d c46_0 = _mm_load_sd(&C[(l_n*84)+10]);
    __m128d a46_0 = _mm_load_sd(&A[202]);
#if defined(__SSE3__) && defined(__AVX__)
    c46_0 = _mm_add_sd(c46_0, _mm_mul_sd(a46_0, _mm256_castpd256_pd128(b46)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c46_0 = _mm_add_sd(c46_0, _mm_mul_sd(a46_0, b46));
#endif
    _mm_store_sd(&C[(l_n*84)+10], c46_0);
    __m128d c46_1 = _mm_load_sd(&C[(l_n*84)+25]);
    __m128d a46_1 = _mm_load_sd(&A[203]);
#if defined(__SSE3__) && defined(__AVX__)
    c46_1 = _mm_add_sd(c46_1, _mm_mul_sd(a46_1, _mm256_castpd256_pd128(b46)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c46_1 = _mm_add_sd(c46_1, _mm_mul_sd(a46_1, b46));
#endif
    _mm_store_sd(&C[(l_n*84)+25], c46_1);
    __m128d c46_2 = _mm_load_sd(&C[(l_n*84)+46]);
    __m128d a46_2 = _mm_load_sd(&A[204]);
#if defined(__SSE3__) && defined(__AVX__)
    c46_2 = _mm_add_sd(c46_2, _mm_mul_sd(a46_2, _mm256_castpd256_pd128(b46)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c46_2 = _mm_add_sd(c46_2, _mm_mul_sd(a46_2, b46));
#endif
    _mm_store_sd(&C[(l_n*84)+46], c46_2);
    __m128d c46_3 = _mm_load_sd(&C[(l_n*84)+74]);
    __m128d a46_3 = _mm_load_sd(&A[205]);
#if defined(__SSE3__) && defined(__AVX__)
    c46_3 = _mm_add_sd(c46_3, _mm_mul_sd(a46_3, _mm256_castpd256_pd128(b46)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c46_3 = _mm_add_sd(c46_3, _mm_mul_sd(a46_3, b46));
#endif
    _mm_store_sd(&C[(l_n*84)+74], c46_3);
#else
    C[(l_n*84)+10] += A[202] * B[(l_n*84)+46];
    C[(l_n*84)+25] += A[203] * B[(l_n*84)+46];
    C[(l_n*84)+46] += A[204] * B[(l_n*84)+46];
    C[(l_n*84)+74] += A[205] * B[(l_n*84)+46];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b47 = _mm256_broadcast_sd(&B[(l_n*84)+47]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b47 = _mm_loaddup_pd(&B[(l_n*84)+47]);
#endif
    __m128d c47_0 = _mm_load_sd(&C[(l_n*84)+11]);
    __m128d a47_0 = _mm_load_sd(&A[206]);
#if defined(__SSE3__) && defined(__AVX__)
    c47_0 = _mm_add_sd(c47_0, _mm_mul_sd(a47_0, _mm256_castpd256_pd128(b47)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c47_0 = _mm_add_sd(c47_0, _mm_mul_sd(a47_0, b47));
#endif
    _mm_store_sd(&C[(l_n*84)+11], c47_0);
    __m128d c47_1 = _mm_load_sd(&C[(l_n*84)+26]);
    __m128d a47_1 = _mm_load_sd(&A[207]);
#if defined(__SSE3__) && defined(__AVX__)
    c47_1 = _mm_add_sd(c47_1, _mm_mul_sd(a47_1, _mm256_castpd256_pd128(b47)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c47_1 = _mm_add_sd(c47_1, _mm_mul_sd(a47_1, b47));
#endif
    _mm_store_sd(&C[(l_n*84)+26], c47_1);
    __m128d c47_2 = _mm_load_sd(&C[(l_n*84)+47]);
    __m128d a47_2 = _mm_load_sd(&A[208]);
#if defined(__SSE3__) && defined(__AVX__)
    c47_2 = _mm_add_sd(c47_2, _mm_mul_sd(a47_2, _mm256_castpd256_pd128(b47)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c47_2 = _mm_add_sd(c47_2, _mm_mul_sd(a47_2, b47));
#endif
    _mm_store_sd(&C[(l_n*84)+47], c47_2);
    __m128d c47_3 = _mm_load_sd(&C[(l_n*84)+75]);
    __m128d a47_3 = _mm_load_sd(&A[209]);
#if defined(__SSE3__) && defined(__AVX__)
    c47_3 = _mm_add_sd(c47_3, _mm_mul_sd(a47_3, _mm256_castpd256_pd128(b47)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c47_3 = _mm_add_sd(c47_3, _mm_mul_sd(a47_3, b47));
#endif
    _mm_store_sd(&C[(l_n*84)+75], c47_3);
#else
    C[(l_n*84)+11] += A[206] * B[(l_n*84)+47];
    C[(l_n*84)+26] += A[207] * B[(l_n*84)+47];
    C[(l_n*84)+47] += A[208] * B[(l_n*84)+47];
    C[(l_n*84)+75] += A[209] * B[(l_n*84)+47];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b48 = _mm256_broadcast_sd(&B[(l_n*84)+48]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b48 = _mm_loaddup_pd(&B[(l_n*84)+48]);
#endif
    __m128d c48_0 = _mm_load_sd(&C[(l_n*84)+12]);
    __m128d a48_0 = _mm_load_sd(&A[210]);
#if defined(__SSE3__) && defined(__AVX__)
    c48_0 = _mm_add_sd(c48_0, _mm_mul_sd(a48_0, _mm256_castpd256_pd128(b48)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c48_0 = _mm_add_sd(c48_0, _mm_mul_sd(a48_0, b48));
#endif
    _mm_store_sd(&C[(l_n*84)+12], c48_0);
    __m128d c48_1 = _mm_load_sd(&C[(l_n*84)+27]);
    __m128d a48_1 = _mm_load_sd(&A[211]);
#if defined(__SSE3__) && defined(__AVX__)
    c48_1 = _mm_add_sd(c48_1, _mm_mul_sd(a48_1, _mm256_castpd256_pd128(b48)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c48_1 = _mm_add_sd(c48_1, _mm_mul_sd(a48_1, b48));
#endif
    _mm_store_sd(&C[(l_n*84)+27], c48_1);
    __m128d c48_2 = _mm_load_sd(&C[(l_n*84)+48]);
    __m128d a48_2 = _mm_load_sd(&A[212]);
#if defined(__SSE3__) && defined(__AVX__)
    c48_2 = _mm_add_sd(c48_2, _mm_mul_sd(a48_2, _mm256_castpd256_pd128(b48)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c48_2 = _mm_add_sd(c48_2, _mm_mul_sd(a48_2, b48));
#endif
    _mm_store_sd(&C[(l_n*84)+48], c48_2);
    __m128d c48_3 = _mm_load_sd(&C[(l_n*84)+76]);
    __m128d a48_3 = _mm_load_sd(&A[213]);
#if defined(__SSE3__) && defined(__AVX__)
    c48_3 = _mm_add_sd(c48_3, _mm_mul_sd(a48_3, _mm256_castpd256_pd128(b48)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c48_3 = _mm_add_sd(c48_3, _mm_mul_sd(a48_3, b48));
#endif
    _mm_store_sd(&C[(l_n*84)+76], c48_3);
#else
    C[(l_n*84)+12] += A[210] * B[(l_n*84)+48];
    C[(l_n*84)+27] += A[211] * B[(l_n*84)+48];
    C[(l_n*84)+48] += A[212] * B[(l_n*84)+48];
    C[(l_n*84)+76] += A[213] * B[(l_n*84)+48];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b49 = _mm256_broadcast_sd(&B[(l_n*84)+49]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b49 = _mm_loaddup_pd(&B[(l_n*84)+49]);
#endif
    __m128d c49_0 = _mm_load_sd(&C[(l_n*84)+13]);
    __m128d a49_0 = _mm_load_sd(&A[214]);
#if defined(__SSE3__) && defined(__AVX__)
    c49_0 = _mm_add_sd(c49_0, _mm_mul_sd(a49_0, _mm256_castpd256_pd128(b49)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c49_0 = _mm_add_sd(c49_0, _mm_mul_sd(a49_0, b49));
#endif
    _mm_store_sd(&C[(l_n*84)+13], c49_0);
    __m128d c49_1 = _mm_load_sd(&C[(l_n*84)+28]);
    __m128d a49_1 = _mm_load_sd(&A[215]);
#if defined(__SSE3__) && defined(__AVX__)
    c49_1 = _mm_add_sd(c49_1, _mm_mul_sd(a49_1, _mm256_castpd256_pd128(b49)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c49_1 = _mm_add_sd(c49_1, _mm_mul_sd(a49_1, b49));
#endif
    _mm_store_sd(&C[(l_n*84)+28], c49_1);
    __m128d c49_2 = _mm_load_sd(&C[(l_n*84)+49]);
    __m128d a49_2 = _mm_load_sd(&A[216]);
#if defined(__SSE3__) && defined(__AVX__)
    c49_2 = _mm_add_sd(c49_2, _mm_mul_sd(a49_2, _mm256_castpd256_pd128(b49)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c49_2 = _mm_add_sd(c49_2, _mm_mul_sd(a49_2, b49));
#endif
    _mm_store_sd(&C[(l_n*84)+49], c49_2);
    __m128d c49_3 = _mm_load_sd(&C[(l_n*84)+77]);
    __m128d a49_3 = _mm_load_sd(&A[217]);
#if defined(__SSE3__) && defined(__AVX__)
    c49_3 = _mm_add_sd(c49_3, _mm_mul_sd(a49_3, _mm256_castpd256_pd128(b49)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c49_3 = _mm_add_sd(c49_3, _mm_mul_sd(a49_3, b49));
#endif
    _mm_store_sd(&C[(l_n*84)+77], c49_3);
#else
    C[(l_n*84)+13] += A[214] * B[(l_n*84)+49];
    C[(l_n*84)+28] += A[215] * B[(l_n*84)+49];
    C[(l_n*84)+49] += A[216] * B[(l_n*84)+49];
    C[(l_n*84)+77] += A[217] * B[(l_n*84)+49];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b50 = _mm256_broadcast_sd(&B[(l_n*84)+50]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b50 = _mm_loaddup_pd(&B[(l_n*84)+50]);
#endif
    __m128d c50_0 = _mm_load_sd(&C[(l_n*84)+4]);
    __m128d a50_0 = _mm_load_sd(&A[218]);
#if defined(__SSE3__) && defined(__AVX__)
    c50_0 = _mm_add_sd(c50_0, _mm_mul_sd(a50_0, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c50_0 = _mm_add_sd(c50_0, _mm_mul_sd(a50_0, b50));
#endif
    _mm_store_sd(&C[(l_n*84)+4], c50_0);
    __m128d c50_1 = _mm_load_sd(&C[(l_n*84)+14]);
    __m128d a50_1 = _mm_load_sd(&A[219]);
#if defined(__SSE3__) && defined(__AVX__)
    c50_1 = _mm_add_sd(c50_1, _mm_mul_sd(a50_1, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c50_1 = _mm_add_sd(c50_1, _mm_mul_sd(a50_1, b50));
#endif
    _mm_store_sd(&C[(l_n*84)+14], c50_1);
    __m128d c50_2 = _mm_load_sd(&C[(l_n*84)+29]);
    __m128d a50_2 = _mm_load_sd(&A[220]);
#if defined(__SSE3__) && defined(__AVX__)
    c50_2 = _mm_add_sd(c50_2, _mm_mul_sd(a50_2, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c50_2 = _mm_add_sd(c50_2, _mm_mul_sd(a50_2, b50));
#endif
    _mm_store_sd(&C[(l_n*84)+29], c50_2);
    __m128d c50_3 = _mm_load_sd(&C[(l_n*84)+50]);
    __m128d a50_3 = _mm_load_sd(&A[221]);
#if defined(__SSE3__) && defined(__AVX__)
    c50_3 = _mm_add_sd(c50_3, _mm_mul_sd(a50_3, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c50_3 = _mm_add_sd(c50_3, _mm_mul_sd(a50_3, b50));
#endif
    _mm_store_sd(&C[(l_n*84)+50], c50_3);
    __m128d c50_4 = _mm_load_sd(&C[(l_n*84)+78]);
    __m128d a50_4 = _mm_load_sd(&A[222]);
#if defined(__SSE3__) && defined(__AVX__)
    c50_4 = _mm_add_sd(c50_4, _mm_mul_sd(a50_4, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c50_4 = _mm_add_sd(c50_4, _mm_mul_sd(a50_4, b50));
#endif
    _mm_store_sd(&C[(l_n*84)+78], c50_4);
#else
    C[(l_n*84)+4] += A[218] * B[(l_n*84)+50];
    C[(l_n*84)+14] += A[219] * B[(l_n*84)+50];
    C[(l_n*84)+29] += A[220] * B[(l_n*84)+50];
    C[(l_n*84)+50] += A[221] * B[(l_n*84)+50];
    C[(l_n*84)+78] += A[222] * B[(l_n*84)+50];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b51 = _mm256_broadcast_sd(&B[(l_n*84)+51]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b51 = _mm_loaddup_pd(&B[(l_n*84)+51]);
#endif
    __m128d c51_0 = _mm_load_sd(&C[(l_n*84)+5]);
    __m128d a51_0 = _mm_load_sd(&A[223]);
#if defined(__SSE3__) && defined(__AVX__)
    c51_0 = _mm_add_sd(c51_0, _mm_mul_sd(a51_0, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c51_0 = _mm_add_sd(c51_0, _mm_mul_sd(a51_0, b51));
#endif
    _mm_store_sd(&C[(l_n*84)+5], c51_0);
    __m128d c51_1 = _mm_load_sd(&C[(l_n*84)+15]);
    __m128d a51_1 = _mm_load_sd(&A[224]);
#if defined(__SSE3__) && defined(__AVX__)
    c51_1 = _mm_add_sd(c51_1, _mm_mul_sd(a51_1, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c51_1 = _mm_add_sd(c51_1, _mm_mul_sd(a51_1, b51));
#endif
    _mm_store_sd(&C[(l_n*84)+15], c51_1);
    __m128d c51_2 = _mm_load_sd(&C[(l_n*84)+30]);
    __m128d a51_2 = _mm_load_sd(&A[225]);
#if defined(__SSE3__) && defined(__AVX__)
    c51_2 = _mm_add_sd(c51_2, _mm_mul_sd(a51_2, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c51_2 = _mm_add_sd(c51_2, _mm_mul_sd(a51_2, b51));
#endif
    _mm_store_sd(&C[(l_n*84)+30], c51_2);
    __m128d c51_3 = _mm_load_sd(&C[(l_n*84)+51]);
    __m128d a51_3 = _mm_load_sd(&A[226]);
#if defined(__SSE3__) && defined(__AVX__)
    c51_3 = _mm_add_sd(c51_3, _mm_mul_sd(a51_3, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c51_3 = _mm_add_sd(c51_3, _mm_mul_sd(a51_3, b51));
#endif
    _mm_store_sd(&C[(l_n*84)+51], c51_3);
    __m128d c51_4 = _mm_load_sd(&C[(l_n*84)+79]);
    __m128d a51_4 = _mm_load_sd(&A[227]);
#if defined(__SSE3__) && defined(__AVX__)
    c51_4 = _mm_add_sd(c51_4, _mm_mul_sd(a51_4, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c51_4 = _mm_add_sd(c51_4, _mm_mul_sd(a51_4, b51));
#endif
    _mm_store_sd(&C[(l_n*84)+79], c51_4);
#else
    C[(l_n*84)+5] += A[223] * B[(l_n*84)+51];
    C[(l_n*84)+15] += A[224] * B[(l_n*84)+51];
    C[(l_n*84)+30] += A[225] * B[(l_n*84)+51];
    C[(l_n*84)+51] += A[226] * B[(l_n*84)+51];
    C[(l_n*84)+79] += A[227] * B[(l_n*84)+51];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b52 = _mm256_broadcast_sd(&B[(l_n*84)+52]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b52 = _mm_loaddup_pd(&B[(l_n*84)+52]);
#endif
    __m128d c52_0 = _mm_load_sd(&C[(l_n*84)+6]);
    __m128d a52_0 = _mm_load_sd(&A[228]);
#if defined(__SSE3__) && defined(__AVX__)
    c52_0 = _mm_add_sd(c52_0, _mm_mul_sd(a52_0, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c52_0 = _mm_add_sd(c52_0, _mm_mul_sd(a52_0, b52));
#endif
    _mm_store_sd(&C[(l_n*84)+6], c52_0);
    __m128d c52_1 = _mm_load_sd(&C[(l_n*84)+16]);
    __m128d a52_1 = _mm_load_sd(&A[229]);
#if defined(__SSE3__) && defined(__AVX__)
    c52_1 = _mm_add_sd(c52_1, _mm_mul_sd(a52_1, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c52_1 = _mm_add_sd(c52_1, _mm_mul_sd(a52_1, b52));
#endif
    _mm_store_sd(&C[(l_n*84)+16], c52_1);
    __m128d c52_2 = _mm_load_sd(&C[(l_n*84)+31]);
    __m128d a52_2 = _mm_load_sd(&A[230]);
#if defined(__SSE3__) && defined(__AVX__)
    c52_2 = _mm_add_sd(c52_2, _mm_mul_sd(a52_2, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c52_2 = _mm_add_sd(c52_2, _mm_mul_sd(a52_2, b52));
#endif
    _mm_store_sd(&C[(l_n*84)+31], c52_2);
    __m128d c52_3 = _mm_load_sd(&C[(l_n*84)+52]);
    __m128d a52_3 = _mm_load_sd(&A[231]);
#if defined(__SSE3__) && defined(__AVX__)
    c52_3 = _mm_add_sd(c52_3, _mm_mul_sd(a52_3, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c52_3 = _mm_add_sd(c52_3, _mm_mul_sd(a52_3, b52));
#endif
    _mm_store_sd(&C[(l_n*84)+52], c52_3);
    __m128d c52_4 = _mm_load_sd(&C[(l_n*84)+80]);
    __m128d a52_4 = _mm_load_sd(&A[232]);
#if defined(__SSE3__) && defined(__AVX__)
    c52_4 = _mm_add_sd(c52_4, _mm_mul_sd(a52_4, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c52_4 = _mm_add_sd(c52_4, _mm_mul_sd(a52_4, b52));
#endif
    _mm_store_sd(&C[(l_n*84)+80], c52_4);
#else
    C[(l_n*84)+6] += A[228] * B[(l_n*84)+52];
    C[(l_n*84)+16] += A[229] * B[(l_n*84)+52];
    C[(l_n*84)+31] += A[230] * B[(l_n*84)+52];
    C[(l_n*84)+52] += A[231] * B[(l_n*84)+52];
    C[(l_n*84)+80] += A[232] * B[(l_n*84)+52];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b53 = _mm256_broadcast_sd(&B[(l_n*84)+53]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b53 = _mm_loaddup_pd(&B[(l_n*84)+53]);
#endif
    __m128d c53_0 = _mm_load_sd(&C[(l_n*84)+1]);
    __m128d a53_0 = _mm_load_sd(&A[233]);
#if defined(__SSE3__) && defined(__AVX__)
    c53_0 = _mm_add_sd(c53_0, _mm_mul_sd(a53_0, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c53_0 = _mm_add_sd(c53_0, _mm_mul_sd(a53_0, b53));
#endif
    _mm_store_sd(&C[(l_n*84)+1], c53_0);
    __m128d c53_1 = _mm_load_sd(&C[(l_n*84)+7]);
    __m128d a53_1 = _mm_load_sd(&A[234]);
#if defined(__SSE3__) && defined(__AVX__)
    c53_1 = _mm_add_sd(c53_1, _mm_mul_sd(a53_1, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c53_1 = _mm_add_sd(c53_1, _mm_mul_sd(a53_1, b53));
#endif
    _mm_store_sd(&C[(l_n*84)+7], c53_1);
    __m128d c53_2 = _mm_load_sd(&C[(l_n*84)+17]);
    __m128d a53_2 = _mm_load_sd(&A[235]);
#if defined(__SSE3__) && defined(__AVX__)
    c53_2 = _mm_add_sd(c53_2, _mm_mul_sd(a53_2, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c53_2 = _mm_add_sd(c53_2, _mm_mul_sd(a53_2, b53));
#endif
    _mm_store_sd(&C[(l_n*84)+17], c53_2);
    __m128d c53_3 = _mm_load_sd(&C[(l_n*84)+32]);
    __m128d a53_3 = _mm_load_sd(&A[236]);
#if defined(__SSE3__) && defined(__AVX__)
    c53_3 = _mm_add_sd(c53_3, _mm_mul_sd(a53_3, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c53_3 = _mm_add_sd(c53_3, _mm_mul_sd(a53_3, b53));
#endif
    _mm_store_sd(&C[(l_n*84)+32], c53_3);
    __m128d c53_4 = _mm_load_sd(&C[(l_n*84)+53]);
    __m128d a53_4 = _mm_load_sd(&A[237]);
#if defined(__SSE3__) && defined(__AVX__)
    c53_4 = _mm_add_sd(c53_4, _mm_mul_sd(a53_4, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c53_4 = _mm_add_sd(c53_4, _mm_mul_sd(a53_4, b53));
#endif
    _mm_store_sd(&C[(l_n*84)+53], c53_4);
    __m128d c53_5 = _mm_load_sd(&C[(l_n*84)+81]);
    __m128d a53_5 = _mm_load_sd(&A[238]);
#if defined(__SSE3__) && defined(__AVX__)
    c53_5 = _mm_add_sd(c53_5, _mm_mul_sd(a53_5, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c53_5 = _mm_add_sd(c53_5, _mm_mul_sd(a53_5, b53));
#endif
    _mm_store_sd(&C[(l_n*84)+81], c53_5);
#else
    C[(l_n*84)+1] += A[233] * B[(l_n*84)+53];
    C[(l_n*84)+7] += A[234] * B[(l_n*84)+53];
    C[(l_n*84)+17] += A[235] * B[(l_n*84)+53];
    C[(l_n*84)+32] += A[236] * B[(l_n*84)+53];
    C[(l_n*84)+53] += A[237] * B[(l_n*84)+53];
    C[(l_n*84)+81] += A[238] * B[(l_n*84)+53];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b54 = _mm256_broadcast_sd(&B[(l_n*84)+54]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b54 = _mm_loaddup_pd(&B[(l_n*84)+54]);
#endif
    __m128d c54_0 = _mm_load_sd(&C[(l_n*84)+2]);
    __m128d a54_0 = _mm_load_sd(&A[239]);
#if defined(__SSE3__) && defined(__AVX__)
    c54_0 = _mm_add_sd(c54_0, _mm_mul_sd(a54_0, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c54_0 = _mm_add_sd(c54_0, _mm_mul_sd(a54_0, b54));
#endif
    _mm_store_sd(&C[(l_n*84)+2], c54_0);
    __m128d c54_1 = _mm_load_sd(&C[(l_n*84)+8]);
    __m128d a54_1 = _mm_load_sd(&A[240]);
#if defined(__SSE3__) && defined(__AVX__)
    c54_1 = _mm_add_sd(c54_1, _mm_mul_sd(a54_1, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c54_1 = _mm_add_sd(c54_1, _mm_mul_sd(a54_1, b54));
#endif
    _mm_store_sd(&C[(l_n*84)+8], c54_1);
    __m128d c54_2 = _mm_load_sd(&C[(l_n*84)+18]);
    __m128d a54_2 = _mm_load_sd(&A[241]);
#if defined(__SSE3__) && defined(__AVX__)
    c54_2 = _mm_add_sd(c54_2, _mm_mul_sd(a54_2, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c54_2 = _mm_add_sd(c54_2, _mm_mul_sd(a54_2, b54));
#endif
    _mm_store_sd(&C[(l_n*84)+18], c54_2);
    __m128d c54_3 = _mm_load_sd(&C[(l_n*84)+33]);
    __m128d a54_3 = _mm_load_sd(&A[242]);
#if defined(__SSE3__) && defined(__AVX__)
    c54_3 = _mm_add_sd(c54_3, _mm_mul_sd(a54_3, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c54_3 = _mm_add_sd(c54_3, _mm_mul_sd(a54_3, b54));
#endif
    _mm_store_sd(&C[(l_n*84)+33], c54_3);
    __m128d c54_4 = _mm_load_sd(&C[(l_n*84)+54]);
    __m128d a54_4 = _mm_load_sd(&A[243]);
#if defined(__SSE3__) && defined(__AVX__)
    c54_4 = _mm_add_sd(c54_4, _mm_mul_sd(a54_4, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c54_4 = _mm_add_sd(c54_4, _mm_mul_sd(a54_4, b54));
#endif
    _mm_store_sd(&C[(l_n*84)+54], c54_4);
    __m128d c54_5 = _mm_load_sd(&C[(l_n*84)+82]);
    __m128d a54_5 = _mm_load_sd(&A[244]);
#if defined(__SSE3__) && defined(__AVX__)
    c54_5 = _mm_add_sd(c54_5, _mm_mul_sd(a54_5, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c54_5 = _mm_add_sd(c54_5, _mm_mul_sd(a54_5, b54));
#endif
    _mm_store_sd(&C[(l_n*84)+82], c54_5);
#else
    C[(l_n*84)+2] += A[239] * B[(l_n*84)+54];
    C[(l_n*84)+8] += A[240] * B[(l_n*84)+54];
    C[(l_n*84)+18] += A[241] * B[(l_n*84)+54];
    C[(l_n*84)+33] += A[242] * B[(l_n*84)+54];
    C[(l_n*84)+54] += A[243] * B[(l_n*84)+54];
    C[(l_n*84)+82] += A[244] * B[(l_n*84)+54];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b55 = _mm256_broadcast_sd(&B[(l_n*84)+55]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b55 = _mm_loaddup_pd(&B[(l_n*84)+55]);
#endif
    __m128d c55_0 = _mm_load_sd(&C[(l_n*84)+0]);
    __m128d a55_0 = _mm_load_sd(&A[245]);
#if defined(__SSE3__) && defined(__AVX__)
    c55_0 = _mm_add_sd(c55_0, _mm_mul_sd(a55_0, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c55_0 = _mm_add_sd(c55_0, _mm_mul_sd(a55_0, b55));
#endif
    _mm_store_sd(&C[(l_n*84)+0], c55_0);
    __m128d c55_1 = _mm_load_sd(&C[(l_n*84)+3]);
    __m128d a55_1 = _mm_load_sd(&A[246]);
#if defined(__SSE3__) && defined(__AVX__)
    c55_1 = _mm_add_sd(c55_1, _mm_mul_sd(a55_1, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c55_1 = _mm_add_sd(c55_1, _mm_mul_sd(a55_1, b55));
#endif
    _mm_store_sd(&C[(l_n*84)+3], c55_1);
    __m128d c55_2 = _mm_load_sd(&C[(l_n*84)+9]);
    __m128d a55_2 = _mm_load_sd(&A[247]);
#if defined(__SSE3__) && defined(__AVX__)
    c55_2 = _mm_add_sd(c55_2, _mm_mul_sd(a55_2, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c55_2 = _mm_add_sd(c55_2, _mm_mul_sd(a55_2, b55));
#endif
    _mm_store_sd(&C[(l_n*84)+9], c55_2);
    __m128d c55_3 = _mm_load_sd(&C[(l_n*84)+19]);
    __m128d a55_3 = _mm_load_sd(&A[248]);
#if defined(__SSE3__) && defined(__AVX__)
    c55_3 = _mm_add_sd(c55_3, _mm_mul_sd(a55_3, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c55_3 = _mm_add_sd(c55_3, _mm_mul_sd(a55_3, b55));
#endif
    _mm_store_sd(&C[(l_n*84)+19], c55_3);
    __m128d c55_4 = _mm_load_sd(&C[(l_n*84)+34]);
    __m128d a55_4 = _mm_load_sd(&A[249]);
#if defined(__SSE3__) && defined(__AVX__)
    c55_4 = _mm_add_sd(c55_4, _mm_mul_sd(a55_4, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c55_4 = _mm_add_sd(c55_4, _mm_mul_sd(a55_4, b55));
#endif
    _mm_store_sd(&C[(l_n*84)+34], c55_4);
    __m128d c55_5 = _mm_load_sd(&C[(l_n*84)+55]);
    __m128d a55_5 = _mm_load_sd(&A[250]);
#if defined(__SSE3__) && defined(__AVX__)
    c55_5 = _mm_add_sd(c55_5, _mm_mul_sd(a55_5, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c55_5 = _mm_add_sd(c55_5, _mm_mul_sd(a55_5, b55));
#endif
    _mm_store_sd(&C[(l_n*84)+55], c55_5);
    __m128d c55_6 = _mm_load_sd(&C[(l_n*84)+83]);
    __m128d a55_6 = _mm_load_sd(&A[251]);
#if defined(__SSE3__) && defined(__AVX__)
    c55_6 = _mm_add_sd(c55_6, _mm_mul_sd(a55_6, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c55_6 = _mm_add_sd(c55_6, _mm_mul_sd(a55_6, b55));
#endif
    _mm_store_sd(&C[(l_n*84)+83], c55_6);
#else
    C[(l_n*84)+0] += A[245] * B[(l_n*84)+55];
    C[(l_n*84)+3] += A[246] * B[(l_n*84)+55];
    C[(l_n*84)+9] += A[247] * B[(l_n*84)+55];
    C[(l_n*84)+19] += A[248] * B[(l_n*84)+55];
    C[(l_n*84)+34] += A[249] * B[(l_n*84)+55];
    C[(l_n*84)+55] += A[250] * B[(l_n*84)+55];
    C[(l_n*84)+83] += A[251] * B[(l_n*84)+55];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b56 = _mm256_broadcast_sd(&B[(l_n*84)+56]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b56 = _mm_loaddup_pd(&B[(l_n*84)+56]);
#endif
    __m128d c56_0 = _mm_load_sd(&C[(l_n*84)+56]);
    __m128d a56_0 = _mm_load_sd(&A[252]);
#if defined(__SSE3__) && defined(__AVX__)
    c56_0 = _mm_add_sd(c56_0, _mm_mul_sd(a56_0, _mm256_castpd256_pd128(b56)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c56_0 = _mm_add_sd(c56_0, _mm_mul_sd(a56_0, b56));
#endif
    _mm_store_sd(&C[(l_n*84)+56], c56_0);
#else
    C[(l_n*84)+56] += A[252] * B[(l_n*84)+56];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b57 = _mm256_broadcast_sd(&B[(l_n*84)+57]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b57 = _mm_loaddup_pd(&B[(l_n*84)+57]);
#endif
    __m128d c57_0 = _mm_load_sd(&C[(l_n*84)+57]);
    __m128d a57_0 = _mm_load_sd(&A[253]);
#if defined(__SSE3__) && defined(__AVX__)
    c57_0 = _mm_add_sd(c57_0, _mm_mul_sd(a57_0, _mm256_castpd256_pd128(b57)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c57_0 = _mm_add_sd(c57_0, _mm_mul_sd(a57_0, b57));
#endif
    _mm_store_sd(&C[(l_n*84)+57], c57_0);
#else
    C[(l_n*84)+57] += A[253] * B[(l_n*84)+57];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b58 = _mm256_broadcast_sd(&B[(l_n*84)+58]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b58 = _mm_loaddup_pd(&B[(l_n*84)+58]);
#endif
    __m128d c58_0 = _mm_load_sd(&C[(l_n*84)+58]);
    __m128d a58_0 = _mm_load_sd(&A[254]);
#if defined(__SSE3__) && defined(__AVX__)
    c58_0 = _mm_add_sd(c58_0, _mm_mul_sd(a58_0, _mm256_castpd256_pd128(b58)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c58_0 = _mm_add_sd(c58_0, _mm_mul_sd(a58_0, b58));
#endif
    _mm_store_sd(&C[(l_n*84)+58], c58_0);
#else
    C[(l_n*84)+58] += A[254] * B[(l_n*84)+58];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b59 = _mm256_broadcast_sd(&B[(l_n*84)+59]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b59 = _mm_loaddup_pd(&B[(l_n*84)+59]);
#endif
    __m128d c59_0 = _mm_load_sd(&C[(l_n*84)+59]);
    __m128d a59_0 = _mm_load_sd(&A[255]);
#if defined(__SSE3__) && defined(__AVX__)
    c59_0 = _mm_add_sd(c59_0, _mm_mul_sd(a59_0, _mm256_castpd256_pd128(b59)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c59_0 = _mm_add_sd(c59_0, _mm_mul_sd(a59_0, b59));
#endif
    _mm_store_sd(&C[(l_n*84)+59], c59_0);
#else
    C[(l_n*84)+59] += A[255] * B[(l_n*84)+59];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b60 = _mm256_broadcast_sd(&B[(l_n*84)+60]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b60 = _mm_loaddup_pd(&B[(l_n*84)+60]);
#endif
    __m128d c60_0 = _mm_load_sd(&C[(l_n*84)+60]);
    __m128d a60_0 = _mm_load_sd(&A[256]);
#if defined(__SSE3__) && defined(__AVX__)
    c60_0 = _mm_add_sd(c60_0, _mm_mul_sd(a60_0, _mm256_castpd256_pd128(b60)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c60_0 = _mm_add_sd(c60_0, _mm_mul_sd(a60_0, b60));
#endif
    _mm_store_sd(&C[(l_n*84)+60], c60_0);
#else
    C[(l_n*84)+60] += A[256] * B[(l_n*84)+60];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b61 = _mm256_broadcast_sd(&B[(l_n*84)+61]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b61 = _mm_loaddup_pd(&B[(l_n*84)+61]);
#endif
    __m128d c61_0 = _mm_load_sd(&C[(l_n*84)+61]);
    __m128d a61_0 = _mm_load_sd(&A[257]);
#if defined(__SSE3__) && defined(__AVX__)
    c61_0 = _mm_add_sd(c61_0, _mm_mul_sd(a61_0, _mm256_castpd256_pd128(b61)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c61_0 = _mm_add_sd(c61_0, _mm_mul_sd(a61_0, b61));
#endif
    _mm_store_sd(&C[(l_n*84)+61], c61_0);
#else
    C[(l_n*84)+61] += A[257] * B[(l_n*84)+61];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b62 = _mm256_broadcast_sd(&B[(l_n*84)+62]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b62 = _mm_loaddup_pd(&B[(l_n*84)+62]);
#endif
    __m128d c62_0 = _mm_load_sd(&C[(l_n*84)+62]);
    __m128d a62_0 = _mm_load_sd(&A[258]);
#if defined(__SSE3__) && defined(__AVX__)
    c62_0 = _mm_add_sd(c62_0, _mm_mul_sd(a62_0, _mm256_castpd256_pd128(b62)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c62_0 = _mm_add_sd(c62_0, _mm_mul_sd(a62_0, b62));
#endif
    _mm_store_sd(&C[(l_n*84)+62], c62_0);
#else
    C[(l_n*84)+62] += A[258] * B[(l_n*84)+62];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b63 = _mm256_broadcast_sd(&B[(l_n*84)+63]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b63 = _mm_loaddup_pd(&B[(l_n*84)+63]);
#endif
    __m128d c63_0 = _mm_load_sd(&C[(l_n*84)+35]);
    __m128d a63_0 = _mm_load_sd(&A[259]);
#if defined(__SSE3__) && defined(__AVX__)
    c63_0 = _mm_add_sd(c63_0, _mm_mul_sd(a63_0, _mm256_castpd256_pd128(b63)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c63_0 = _mm_add_sd(c63_0, _mm_mul_sd(a63_0, b63));
#endif
    _mm_store_sd(&C[(l_n*84)+35], c63_0);
    __m128d c63_1 = _mm_load_sd(&C[(l_n*84)+63]);
    __m128d a63_1 = _mm_load_sd(&A[260]);
#if defined(__SSE3__) && defined(__AVX__)
    c63_1 = _mm_add_sd(c63_1, _mm_mul_sd(a63_1, _mm256_castpd256_pd128(b63)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c63_1 = _mm_add_sd(c63_1, _mm_mul_sd(a63_1, b63));
#endif
    _mm_store_sd(&C[(l_n*84)+63], c63_1);
#else
    C[(l_n*84)+35] += A[259] * B[(l_n*84)+63];
    C[(l_n*84)+63] += A[260] * B[(l_n*84)+63];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b64 = _mm256_broadcast_sd(&B[(l_n*84)+64]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b64 = _mm_loaddup_pd(&B[(l_n*84)+64]);
#endif
    __m128d c64_0 = _mm_load_sd(&C[(l_n*84)+36]);
    __m128d a64_0 = _mm_load_sd(&A[261]);
#if defined(__SSE3__) && defined(__AVX__)
    c64_0 = _mm_add_sd(c64_0, _mm_mul_sd(a64_0, _mm256_castpd256_pd128(b64)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c64_0 = _mm_add_sd(c64_0, _mm_mul_sd(a64_0, b64));
#endif
    _mm_store_sd(&C[(l_n*84)+36], c64_0);
    __m128d c64_1 = _mm_load_sd(&C[(l_n*84)+64]);
    __m128d a64_1 = _mm_load_sd(&A[262]);
#if defined(__SSE3__) && defined(__AVX__)
    c64_1 = _mm_add_sd(c64_1, _mm_mul_sd(a64_1, _mm256_castpd256_pd128(b64)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c64_1 = _mm_add_sd(c64_1, _mm_mul_sd(a64_1, b64));
#endif
    _mm_store_sd(&C[(l_n*84)+64], c64_1);
#else
    C[(l_n*84)+36] += A[261] * B[(l_n*84)+64];
    C[(l_n*84)+64] += A[262] * B[(l_n*84)+64];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b65 = _mm256_broadcast_sd(&B[(l_n*84)+65]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b65 = _mm_loaddup_pd(&B[(l_n*84)+65]);
#endif
    __m128d c65_0 = _mm_load_sd(&C[(l_n*84)+37]);
    __m128d a65_0 = _mm_load_sd(&A[263]);
#if defined(__SSE3__) && defined(__AVX__)
    c65_0 = _mm_add_sd(c65_0, _mm_mul_sd(a65_0, _mm256_castpd256_pd128(b65)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c65_0 = _mm_add_sd(c65_0, _mm_mul_sd(a65_0, b65));
#endif
    _mm_store_sd(&C[(l_n*84)+37], c65_0);
    __m128d c65_1 = _mm_load_sd(&C[(l_n*84)+65]);
    __m128d a65_1 = _mm_load_sd(&A[264]);
#if defined(__SSE3__) && defined(__AVX__)
    c65_1 = _mm_add_sd(c65_1, _mm_mul_sd(a65_1, _mm256_castpd256_pd128(b65)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c65_1 = _mm_add_sd(c65_1, _mm_mul_sd(a65_1, b65));
#endif
    _mm_store_sd(&C[(l_n*84)+65], c65_1);
#else
    C[(l_n*84)+37] += A[263] * B[(l_n*84)+65];
    C[(l_n*84)+65] += A[264] * B[(l_n*84)+65];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b66 = _mm256_broadcast_sd(&B[(l_n*84)+66]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b66 = _mm_loaddup_pd(&B[(l_n*84)+66]);
#endif
    __m128d c66_0 = _mm_load_sd(&C[(l_n*84)+38]);
    __m128d a66_0 = _mm_load_sd(&A[265]);
#if defined(__SSE3__) && defined(__AVX__)
    c66_0 = _mm_add_sd(c66_0, _mm_mul_sd(a66_0, _mm256_castpd256_pd128(b66)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c66_0 = _mm_add_sd(c66_0, _mm_mul_sd(a66_0, b66));
#endif
    _mm_store_sd(&C[(l_n*84)+38], c66_0);
    __m128d c66_1 = _mm_load_sd(&C[(l_n*84)+66]);
    __m128d a66_1 = _mm_load_sd(&A[266]);
#if defined(__SSE3__) && defined(__AVX__)
    c66_1 = _mm_add_sd(c66_1, _mm_mul_sd(a66_1, _mm256_castpd256_pd128(b66)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c66_1 = _mm_add_sd(c66_1, _mm_mul_sd(a66_1, b66));
#endif
    _mm_store_sd(&C[(l_n*84)+66], c66_1);
#else
    C[(l_n*84)+38] += A[265] * B[(l_n*84)+66];
    C[(l_n*84)+66] += A[266] * B[(l_n*84)+66];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b67 = _mm256_broadcast_sd(&B[(l_n*84)+67]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b67 = _mm_loaddup_pd(&B[(l_n*84)+67]);
#endif
    __m128d c67_0 = _mm_load_sd(&C[(l_n*84)+39]);
    __m128d a67_0 = _mm_load_sd(&A[267]);
#if defined(__SSE3__) && defined(__AVX__)
    c67_0 = _mm_add_sd(c67_0, _mm_mul_sd(a67_0, _mm256_castpd256_pd128(b67)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c67_0 = _mm_add_sd(c67_0, _mm_mul_sd(a67_0, b67));
#endif
    _mm_store_sd(&C[(l_n*84)+39], c67_0);
    __m128d c67_1 = _mm_load_sd(&C[(l_n*84)+67]);
    __m128d a67_1 = _mm_load_sd(&A[268]);
#if defined(__SSE3__) && defined(__AVX__)
    c67_1 = _mm_add_sd(c67_1, _mm_mul_sd(a67_1, _mm256_castpd256_pd128(b67)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c67_1 = _mm_add_sd(c67_1, _mm_mul_sd(a67_1, b67));
#endif
    _mm_store_sd(&C[(l_n*84)+67], c67_1);
#else
    C[(l_n*84)+39] += A[267] * B[(l_n*84)+67];
    C[(l_n*84)+67] += A[268] * B[(l_n*84)+67];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b68 = _mm256_broadcast_sd(&B[(l_n*84)+68]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b68 = _mm_loaddup_pd(&B[(l_n*84)+68]);
#endif
    __m128d c68_0 = _mm_load_sd(&C[(l_n*84)+40]);
    __m128d a68_0 = _mm_load_sd(&A[269]);
#if defined(__SSE3__) && defined(__AVX__)
    c68_0 = _mm_add_sd(c68_0, _mm_mul_sd(a68_0, _mm256_castpd256_pd128(b68)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c68_0 = _mm_add_sd(c68_0, _mm_mul_sd(a68_0, b68));
#endif
    _mm_store_sd(&C[(l_n*84)+40], c68_0);
    __m128d c68_1 = _mm_load_sd(&C[(l_n*84)+68]);
    __m128d a68_1 = _mm_load_sd(&A[270]);
#if defined(__SSE3__) && defined(__AVX__)
    c68_1 = _mm_add_sd(c68_1, _mm_mul_sd(a68_1, _mm256_castpd256_pd128(b68)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c68_1 = _mm_add_sd(c68_1, _mm_mul_sd(a68_1, b68));
#endif
    _mm_store_sd(&C[(l_n*84)+68], c68_1);
#else
    C[(l_n*84)+40] += A[269] * B[(l_n*84)+68];
    C[(l_n*84)+68] += A[270] * B[(l_n*84)+68];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b69 = _mm256_broadcast_sd(&B[(l_n*84)+69]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b69 = _mm_loaddup_pd(&B[(l_n*84)+69]);
#endif
    __m128d c69_0 = _mm_load_sd(&C[(l_n*84)+20]);
    __m128d a69_0 = _mm_load_sd(&A[271]);
#if defined(__SSE3__) && defined(__AVX__)
    c69_0 = _mm_add_sd(c69_0, _mm_mul_sd(a69_0, _mm256_castpd256_pd128(b69)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c69_0 = _mm_add_sd(c69_0, _mm_mul_sd(a69_0, b69));
#endif
    _mm_store_sd(&C[(l_n*84)+20], c69_0);
    __m128d c69_1 = _mm_load_sd(&C[(l_n*84)+41]);
    __m128d a69_1 = _mm_load_sd(&A[272]);
#if defined(__SSE3__) && defined(__AVX__)
    c69_1 = _mm_add_sd(c69_1, _mm_mul_sd(a69_1, _mm256_castpd256_pd128(b69)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c69_1 = _mm_add_sd(c69_1, _mm_mul_sd(a69_1, b69));
#endif
    _mm_store_sd(&C[(l_n*84)+41], c69_1);
    __m128d c69_2 = _mm_load_sd(&C[(l_n*84)+69]);
    __m128d a69_2 = _mm_load_sd(&A[273]);
#if defined(__SSE3__) && defined(__AVX__)
    c69_2 = _mm_add_sd(c69_2, _mm_mul_sd(a69_2, _mm256_castpd256_pd128(b69)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c69_2 = _mm_add_sd(c69_2, _mm_mul_sd(a69_2, b69));
#endif
    _mm_store_sd(&C[(l_n*84)+69], c69_2);
#else
    C[(l_n*84)+20] += A[271] * B[(l_n*84)+69];
    C[(l_n*84)+41] += A[272] * B[(l_n*84)+69];
    C[(l_n*84)+69] += A[273] * B[(l_n*84)+69];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b70 = _mm256_broadcast_sd(&B[(l_n*84)+70]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b70 = _mm_loaddup_pd(&B[(l_n*84)+70]);
#endif
    __m128d c70_0 = _mm_load_sd(&C[(l_n*84)+21]);
    __m128d a70_0 = _mm_load_sd(&A[274]);
#if defined(__SSE3__) && defined(__AVX__)
    c70_0 = _mm_add_sd(c70_0, _mm_mul_sd(a70_0, _mm256_castpd256_pd128(b70)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c70_0 = _mm_add_sd(c70_0, _mm_mul_sd(a70_0, b70));
#endif
    _mm_store_sd(&C[(l_n*84)+21], c70_0);
    __m128d c70_1 = _mm_load_sd(&C[(l_n*84)+42]);
    __m128d a70_1 = _mm_load_sd(&A[275]);
#if defined(__SSE3__) && defined(__AVX__)
    c70_1 = _mm_add_sd(c70_1, _mm_mul_sd(a70_1, _mm256_castpd256_pd128(b70)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c70_1 = _mm_add_sd(c70_1, _mm_mul_sd(a70_1, b70));
#endif
    _mm_store_sd(&C[(l_n*84)+42], c70_1);
    __m128d c70_2 = _mm_load_sd(&C[(l_n*84)+70]);
    __m128d a70_2 = _mm_load_sd(&A[276]);
#if defined(__SSE3__) && defined(__AVX__)
    c70_2 = _mm_add_sd(c70_2, _mm_mul_sd(a70_2, _mm256_castpd256_pd128(b70)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c70_2 = _mm_add_sd(c70_2, _mm_mul_sd(a70_2, b70));
#endif
    _mm_store_sd(&C[(l_n*84)+70], c70_2);
#else
    C[(l_n*84)+21] += A[274] * B[(l_n*84)+70];
    C[(l_n*84)+42] += A[275] * B[(l_n*84)+70];
    C[(l_n*84)+70] += A[276] * B[(l_n*84)+70];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b71 = _mm256_broadcast_sd(&B[(l_n*84)+71]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b71 = _mm_loaddup_pd(&B[(l_n*84)+71]);
#endif
    __m128d c71_0 = _mm_load_sd(&C[(l_n*84)+22]);
    __m128d a71_0 = _mm_load_sd(&A[277]);
#if defined(__SSE3__) && defined(__AVX__)
    c71_0 = _mm_add_sd(c71_0, _mm_mul_sd(a71_0, _mm256_castpd256_pd128(b71)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c71_0 = _mm_add_sd(c71_0, _mm_mul_sd(a71_0, b71));
#endif
    _mm_store_sd(&C[(l_n*84)+22], c71_0);
    __m128d c71_1 = _mm_load_sd(&C[(l_n*84)+43]);
    __m128d a71_1 = _mm_load_sd(&A[278]);
#if defined(__SSE3__) && defined(__AVX__)
    c71_1 = _mm_add_sd(c71_1, _mm_mul_sd(a71_1, _mm256_castpd256_pd128(b71)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c71_1 = _mm_add_sd(c71_1, _mm_mul_sd(a71_1, b71));
#endif
    _mm_store_sd(&C[(l_n*84)+43], c71_1);
    __m128d c71_2 = _mm_load_sd(&C[(l_n*84)+71]);
    __m128d a71_2 = _mm_load_sd(&A[279]);
#if defined(__SSE3__) && defined(__AVX__)
    c71_2 = _mm_add_sd(c71_2, _mm_mul_sd(a71_2, _mm256_castpd256_pd128(b71)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c71_2 = _mm_add_sd(c71_2, _mm_mul_sd(a71_2, b71));
#endif
    _mm_store_sd(&C[(l_n*84)+71], c71_2);
#else
    C[(l_n*84)+22] += A[277] * B[(l_n*84)+71];
    C[(l_n*84)+43] += A[278] * B[(l_n*84)+71];
    C[(l_n*84)+71] += A[279] * B[(l_n*84)+71];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b72 = _mm256_broadcast_sd(&B[(l_n*84)+72]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b72 = _mm_loaddup_pd(&B[(l_n*84)+72]);
#endif
    __m128d c72_0 = _mm_load_sd(&C[(l_n*84)+23]);
    __m128d a72_0 = _mm_load_sd(&A[280]);
#if defined(__SSE3__) && defined(__AVX__)
    c72_0 = _mm_add_sd(c72_0, _mm_mul_sd(a72_0, _mm256_castpd256_pd128(b72)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c72_0 = _mm_add_sd(c72_0, _mm_mul_sd(a72_0, b72));
#endif
    _mm_store_sd(&C[(l_n*84)+23], c72_0);
    __m128d c72_1 = _mm_load_sd(&C[(l_n*84)+44]);
    __m128d a72_1 = _mm_load_sd(&A[281]);
#if defined(__SSE3__) && defined(__AVX__)
    c72_1 = _mm_add_sd(c72_1, _mm_mul_sd(a72_1, _mm256_castpd256_pd128(b72)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c72_1 = _mm_add_sd(c72_1, _mm_mul_sd(a72_1, b72));
#endif
    _mm_store_sd(&C[(l_n*84)+44], c72_1);
    __m128d c72_2 = _mm_load_sd(&C[(l_n*84)+72]);
    __m128d a72_2 = _mm_load_sd(&A[282]);
#if defined(__SSE3__) && defined(__AVX__)
    c72_2 = _mm_add_sd(c72_2, _mm_mul_sd(a72_2, _mm256_castpd256_pd128(b72)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c72_2 = _mm_add_sd(c72_2, _mm_mul_sd(a72_2, b72));
#endif
    _mm_store_sd(&C[(l_n*84)+72], c72_2);
#else
    C[(l_n*84)+23] += A[280] * B[(l_n*84)+72];
    C[(l_n*84)+44] += A[281] * B[(l_n*84)+72];
    C[(l_n*84)+72] += A[282] * B[(l_n*84)+72];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b73 = _mm256_broadcast_sd(&B[(l_n*84)+73]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b73 = _mm_loaddup_pd(&B[(l_n*84)+73]);
#endif
    __m128d c73_0 = _mm_load_sd(&C[(l_n*84)+24]);
    __m128d a73_0 = _mm_load_sd(&A[283]);
#if defined(__SSE3__) && defined(__AVX__)
    c73_0 = _mm_add_sd(c73_0, _mm_mul_sd(a73_0, _mm256_castpd256_pd128(b73)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c73_0 = _mm_add_sd(c73_0, _mm_mul_sd(a73_0, b73));
#endif
    _mm_store_sd(&C[(l_n*84)+24], c73_0);
    __m128d c73_1 = _mm_load_sd(&C[(l_n*84)+45]);
    __m128d a73_1 = _mm_load_sd(&A[284]);
#if defined(__SSE3__) && defined(__AVX__)
    c73_1 = _mm_add_sd(c73_1, _mm_mul_sd(a73_1, _mm256_castpd256_pd128(b73)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c73_1 = _mm_add_sd(c73_1, _mm_mul_sd(a73_1, b73));
#endif
    _mm_store_sd(&C[(l_n*84)+45], c73_1);
    __m128d c73_2 = _mm_load_sd(&C[(l_n*84)+73]);
    __m128d a73_2 = _mm_load_sd(&A[285]);
#if defined(__SSE3__) && defined(__AVX__)
    c73_2 = _mm_add_sd(c73_2, _mm_mul_sd(a73_2, _mm256_castpd256_pd128(b73)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c73_2 = _mm_add_sd(c73_2, _mm_mul_sd(a73_2, b73));
#endif
    _mm_store_sd(&C[(l_n*84)+73], c73_2);
#else
    C[(l_n*84)+24] += A[283] * B[(l_n*84)+73];
    C[(l_n*84)+45] += A[284] * B[(l_n*84)+73];
    C[(l_n*84)+73] += A[285] * B[(l_n*84)+73];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b74 = _mm256_broadcast_sd(&B[(l_n*84)+74]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b74 = _mm_loaddup_pd(&B[(l_n*84)+74]);
#endif
    __m128d c74_0 = _mm_load_sd(&C[(l_n*84)+10]);
    __m128d a74_0 = _mm_load_sd(&A[286]);
#if defined(__SSE3__) && defined(__AVX__)
    c74_0 = _mm_add_sd(c74_0, _mm_mul_sd(a74_0, _mm256_castpd256_pd128(b74)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c74_0 = _mm_add_sd(c74_0, _mm_mul_sd(a74_0, b74));
#endif
    _mm_store_sd(&C[(l_n*84)+10], c74_0);
    __m128d c74_1 = _mm_load_sd(&C[(l_n*84)+25]);
    __m128d a74_1 = _mm_load_sd(&A[287]);
#if defined(__SSE3__) && defined(__AVX__)
    c74_1 = _mm_add_sd(c74_1, _mm_mul_sd(a74_1, _mm256_castpd256_pd128(b74)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c74_1 = _mm_add_sd(c74_1, _mm_mul_sd(a74_1, b74));
#endif
    _mm_store_sd(&C[(l_n*84)+25], c74_1);
    __m128d c74_2 = _mm_load_sd(&C[(l_n*84)+46]);
    __m128d a74_2 = _mm_load_sd(&A[288]);
#if defined(__SSE3__) && defined(__AVX__)
    c74_2 = _mm_add_sd(c74_2, _mm_mul_sd(a74_2, _mm256_castpd256_pd128(b74)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c74_2 = _mm_add_sd(c74_2, _mm_mul_sd(a74_2, b74));
#endif
    _mm_store_sd(&C[(l_n*84)+46], c74_2);
    __m128d c74_3 = _mm_load_sd(&C[(l_n*84)+74]);
    __m128d a74_3 = _mm_load_sd(&A[289]);
#if defined(__SSE3__) && defined(__AVX__)
    c74_3 = _mm_add_sd(c74_3, _mm_mul_sd(a74_3, _mm256_castpd256_pd128(b74)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c74_3 = _mm_add_sd(c74_3, _mm_mul_sd(a74_3, b74));
#endif
    _mm_store_sd(&C[(l_n*84)+74], c74_3);
#else
    C[(l_n*84)+10] += A[286] * B[(l_n*84)+74];
    C[(l_n*84)+25] += A[287] * B[(l_n*84)+74];
    C[(l_n*84)+46] += A[288] * B[(l_n*84)+74];
    C[(l_n*84)+74] += A[289] * B[(l_n*84)+74];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b75 = _mm256_broadcast_sd(&B[(l_n*84)+75]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b75 = _mm_loaddup_pd(&B[(l_n*84)+75]);
#endif
    __m128d c75_0 = _mm_load_sd(&C[(l_n*84)+11]);
    __m128d a75_0 = _mm_load_sd(&A[290]);
#if defined(__SSE3__) && defined(__AVX__)
    c75_0 = _mm_add_sd(c75_0, _mm_mul_sd(a75_0, _mm256_castpd256_pd128(b75)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c75_0 = _mm_add_sd(c75_0, _mm_mul_sd(a75_0, b75));
#endif
    _mm_store_sd(&C[(l_n*84)+11], c75_0);
    __m128d c75_1 = _mm_load_sd(&C[(l_n*84)+26]);
    __m128d a75_1 = _mm_load_sd(&A[291]);
#if defined(__SSE3__) && defined(__AVX__)
    c75_1 = _mm_add_sd(c75_1, _mm_mul_sd(a75_1, _mm256_castpd256_pd128(b75)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c75_1 = _mm_add_sd(c75_1, _mm_mul_sd(a75_1, b75));
#endif
    _mm_store_sd(&C[(l_n*84)+26], c75_1);
    __m128d c75_2 = _mm_load_sd(&C[(l_n*84)+47]);
    __m128d a75_2 = _mm_load_sd(&A[292]);
#if defined(__SSE3__) && defined(__AVX__)
    c75_2 = _mm_add_sd(c75_2, _mm_mul_sd(a75_2, _mm256_castpd256_pd128(b75)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c75_2 = _mm_add_sd(c75_2, _mm_mul_sd(a75_2, b75));
#endif
    _mm_store_sd(&C[(l_n*84)+47], c75_2);
    __m128d c75_3 = _mm_load_sd(&C[(l_n*84)+75]);
    __m128d a75_3 = _mm_load_sd(&A[293]);
#if defined(__SSE3__) && defined(__AVX__)
    c75_3 = _mm_add_sd(c75_3, _mm_mul_sd(a75_3, _mm256_castpd256_pd128(b75)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c75_3 = _mm_add_sd(c75_3, _mm_mul_sd(a75_3, b75));
#endif
    _mm_store_sd(&C[(l_n*84)+75], c75_3);
#else
    C[(l_n*84)+11] += A[290] * B[(l_n*84)+75];
    C[(l_n*84)+26] += A[291] * B[(l_n*84)+75];
    C[(l_n*84)+47] += A[292] * B[(l_n*84)+75];
    C[(l_n*84)+75] += A[293] * B[(l_n*84)+75];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b76 = _mm256_broadcast_sd(&B[(l_n*84)+76]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b76 = _mm_loaddup_pd(&B[(l_n*84)+76]);
#endif
    __m128d c76_0 = _mm_load_sd(&C[(l_n*84)+12]);
    __m128d a76_0 = _mm_load_sd(&A[294]);
#if defined(__SSE3__) && defined(__AVX__)
    c76_0 = _mm_add_sd(c76_0, _mm_mul_sd(a76_0, _mm256_castpd256_pd128(b76)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c76_0 = _mm_add_sd(c76_0, _mm_mul_sd(a76_0, b76));
#endif
    _mm_store_sd(&C[(l_n*84)+12], c76_0);
    __m128d c76_1 = _mm_load_sd(&C[(l_n*84)+27]);
    __m128d a76_1 = _mm_load_sd(&A[295]);
#if defined(__SSE3__) && defined(__AVX__)
    c76_1 = _mm_add_sd(c76_1, _mm_mul_sd(a76_1, _mm256_castpd256_pd128(b76)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c76_1 = _mm_add_sd(c76_1, _mm_mul_sd(a76_1, b76));
#endif
    _mm_store_sd(&C[(l_n*84)+27], c76_1);
    __m128d c76_2 = _mm_load_sd(&C[(l_n*84)+48]);
    __m128d a76_2 = _mm_load_sd(&A[296]);
#if defined(__SSE3__) && defined(__AVX__)
    c76_2 = _mm_add_sd(c76_2, _mm_mul_sd(a76_2, _mm256_castpd256_pd128(b76)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c76_2 = _mm_add_sd(c76_2, _mm_mul_sd(a76_2, b76));
#endif
    _mm_store_sd(&C[(l_n*84)+48], c76_2);
    __m128d c76_3 = _mm_load_sd(&C[(l_n*84)+76]);
    __m128d a76_3 = _mm_load_sd(&A[297]);
#if defined(__SSE3__) && defined(__AVX__)
    c76_3 = _mm_add_sd(c76_3, _mm_mul_sd(a76_3, _mm256_castpd256_pd128(b76)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c76_3 = _mm_add_sd(c76_3, _mm_mul_sd(a76_3, b76));
#endif
    _mm_store_sd(&C[(l_n*84)+76], c76_3);
#else
    C[(l_n*84)+12] += A[294] * B[(l_n*84)+76];
    C[(l_n*84)+27] += A[295] * B[(l_n*84)+76];
    C[(l_n*84)+48] += A[296] * B[(l_n*84)+76];
    C[(l_n*84)+76] += A[297] * B[(l_n*84)+76];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b77 = _mm256_broadcast_sd(&B[(l_n*84)+77]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b77 = _mm_loaddup_pd(&B[(l_n*84)+77]);
#endif
    __m128d c77_0 = _mm_load_sd(&C[(l_n*84)+13]);
    __m128d a77_0 = _mm_load_sd(&A[298]);
#if defined(__SSE3__) && defined(__AVX__)
    c77_0 = _mm_add_sd(c77_0, _mm_mul_sd(a77_0, _mm256_castpd256_pd128(b77)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c77_0 = _mm_add_sd(c77_0, _mm_mul_sd(a77_0, b77));
#endif
    _mm_store_sd(&C[(l_n*84)+13], c77_0);
    __m128d c77_1 = _mm_load_sd(&C[(l_n*84)+28]);
    __m128d a77_1 = _mm_load_sd(&A[299]);
#if defined(__SSE3__) && defined(__AVX__)
    c77_1 = _mm_add_sd(c77_1, _mm_mul_sd(a77_1, _mm256_castpd256_pd128(b77)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c77_1 = _mm_add_sd(c77_1, _mm_mul_sd(a77_1, b77));
#endif
    _mm_store_sd(&C[(l_n*84)+28], c77_1);
    __m128d c77_2 = _mm_load_sd(&C[(l_n*84)+49]);
    __m128d a77_2 = _mm_load_sd(&A[300]);
#if defined(__SSE3__) && defined(__AVX__)
    c77_2 = _mm_add_sd(c77_2, _mm_mul_sd(a77_2, _mm256_castpd256_pd128(b77)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c77_2 = _mm_add_sd(c77_2, _mm_mul_sd(a77_2, b77));
#endif
    _mm_store_sd(&C[(l_n*84)+49], c77_2);
    __m128d c77_3 = _mm_load_sd(&C[(l_n*84)+77]);
    __m128d a77_3 = _mm_load_sd(&A[301]);
#if defined(__SSE3__) && defined(__AVX__)
    c77_3 = _mm_add_sd(c77_3, _mm_mul_sd(a77_3, _mm256_castpd256_pd128(b77)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c77_3 = _mm_add_sd(c77_3, _mm_mul_sd(a77_3, b77));
#endif
    _mm_store_sd(&C[(l_n*84)+77], c77_3);
#else
    C[(l_n*84)+13] += A[298] * B[(l_n*84)+77];
    C[(l_n*84)+28] += A[299] * B[(l_n*84)+77];
    C[(l_n*84)+49] += A[300] * B[(l_n*84)+77];
    C[(l_n*84)+77] += A[301] * B[(l_n*84)+77];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b78 = _mm256_broadcast_sd(&B[(l_n*84)+78]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b78 = _mm_loaddup_pd(&B[(l_n*84)+78]);
#endif
    __m128d c78_0 = _mm_load_sd(&C[(l_n*84)+4]);
    __m128d a78_0 = _mm_load_sd(&A[302]);
#if defined(__SSE3__) && defined(__AVX__)
    c78_0 = _mm_add_sd(c78_0, _mm_mul_sd(a78_0, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c78_0 = _mm_add_sd(c78_0, _mm_mul_sd(a78_0, b78));
#endif
    _mm_store_sd(&C[(l_n*84)+4], c78_0);
    __m128d c78_1 = _mm_load_sd(&C[(l_n*84)+14]);
    __m128d a78_1 = _mm_load_sd(&A[303]);
#if defined(__SSE3__) && defined(__AVX__)
    c78_1 = _mm_add_sd(c78_1, _mm_mul_sd(a78_1, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c78_1 = _mm_add_sd(c78_1, _mm_mul_sd(a78_1, b78));
#endif
    _mm_store_sd(&C[(l_n*84)+14], c78_1);
    __m128d c78_2 = _mm_load_sd(&C[(l_n*84)+29]);
    __m128d a78_2 = _mm_load_sd(&A[304]);
#if defined(__SSE3__) && defined(__AVX__)
    c78_2 = _mm_add_sd(c78_2, _mm_mul_sd(a78_2, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c78_2 = _mm_add_sd(c78_2, _mm_mul_sd(a78_2, b78));
#endif
    _mm_store_sd(&C[(l_n*84)+29], c78_2);
    __m128d c78_3 = _mm_load_sd(&C[(l_n*84)+50]);
    __m128d a78_3 = _mm_load_sd(&A[305]);
#if defined(__SSE3__) && defined(__AVX__)
    c78_3 = _mm_add_sd(c78_3, _mm_mul_sd(a78_3, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c78_3 = _mm_add_sd(c78_3, _mm_mul_sd(a78_3, b78));
#endif
    _mm_store_sd(&C[(l_n*84)+50], c78_3);
    __m128d c78_4 = _mm_load_sd(&C[(l_n*84)+78]);
    __m128d a78_4 = _mm_load_sd(&A[306]);
#if defined(__SSE3__) && defined(__AVX__)
    c78_4 = _mm_add_sd(c78_4, _mm_mul_sd(a78_4, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c78_4 = _mm_add_sd(c78_4, _mm_mul_sd(a78_4, b78));
#endif
    _mm_store_sd(&C[(l_n*84)+78], c78_4);
#else
    C[(l_n*84)+4] += A[302] * B[(l_n*84)+78];
    C[(l_n*84)+14] += A[303] * B[(l_n*84)+78];
    C[(l_n*84)+29] += A[304] * B[(l_n*84)+78];
    C[(l_n*84)+50] += A[305] * B[(l_n*84)+78];
    C[(l_n*84)+78] += A[306] * B[(l_n*84)+78];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b79 = _mm256_broadcast_sd(&B[(l_n*84)+79]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b79 = _mm_loaddup_pd(&B[(l_n*84)+79]);
#endif
    __m128d c79_0 = _mm_load_sd(&C[(l_n*84)+5]);
    __m128d a79_0 = _mm_load_sd(&A[307]);
#if defined(__SSE3__) && defined(__AVX__)
    c79_0 = _mm_add_sd(c79_0, _mm_mul_sd(a79_0, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c79_0 = _mm_add_sd(c79_0, _mm_mul_sd(a79_0, b79));
#endif
    _mm_store_sd(&C[(l_n*84)+5], c79_0);
    __m128d c79_1 = _mm_load_sd(&C[(l_n*84)+15]);
    __m128d a79_1 = _mm_load_sd(&A[308]);
#if defined(__SSE3__) && defined(__AVX__)
    c79_1 = _mm_add_sd(c79_1, _mm_mul_sd(a79_1, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c79_1 = _mm_add_sd(c79_1, _mm_mul_sd(a79_1, b79));
#endif
    _mm_store_sd(&C[(l_n*84)+15], c79_1);
    __m128d c79_2 = _mm_load_sd(&C[(l_n*84)+30]);
    __m128d a79_2 = _mm_load_sd(&A[309]);
#if defined(__SSE3__) && defined(__AVX__)
    c79_2 = _mm_add_sd(c79_2, _mm_mul_sd(a79_2, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c79_2 = _mm_add_sd(c79_2, _mm_mul_sd(a79_2, b79));
#endif
    _mm_store_sd(&C[(l_n*84)+30], c79_2);
    __m128d c79_3 = _mm_load_sd(&C[(l_n*84)+51]);
    __m128d a79_3 = _mm_load_sd(&A[310]);
#if defined(__SSE3__) && defined(__AVX__)
    c79_3 = _mm_add_sd(c79_3, _mm_mul_sd(a79_3, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c79_3 = _mm_add_sd(c79_3, _mm_mul_sd(a79_3, b79));
#endif
    _mm_store_sd(&C[(l_n*84)+51], c79_3);
    __m128d c79_4 = _mm_load_sd(&C[(l_n*84)+79]);
    __m128d a79_4 = _mm_load_sd(&A[311]);
#if defined(__SSE3__) && defined(__AVX__)
    c79_4 = _mm_add_sd(c79_4, _mm_mul_sd(a79_4, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c79_4 = _mm_add_sd(c79_4, _mm_mul_sd(a79_4, b79));
#endif
    _mm_store_sd(&C[(l_n*84)+79], c79_4);
#else
    C[(l_n*84)+5] += A[307] * B[(l_n*84)+79];
    C[(l_n*84)+15] += A[308] * B[(l_n*84)+79];
    C[(l_n*84)+30] += A[309] * B[(l_n*84)+79];
    C[(l_n*84)+51] += A[310] * B[(l_n*84)+79];
    C[(l_n*84)+79] += A[311] * B[(l_n*84)+79];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b80 = _mm256_broadcast_sd(&B[(l_n*84)+80]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b80 = _mm_loaddup_pd(&B[(l_n*84)+80]);
#endif
    __m128d c80_0 = _mm_load_sd(&C[(l_n*84)+6]);
    __m128d a80_0 = _mm_load_sd(&A[312]);
#if defined(__SSE3__) && defined(__AVX__)
    c80_0 = _mm_add_sd(c80_0, _mm_mul_sd(a80_0, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c80_0 = _mm_add_sd(c80_0, _mm_mul_sd(a80_0, b80));
#endif
    _mm_store_sd(&C[(l_n*84)+6], c80_0);
    __m128d c80_1 = _mm_load_sd(&C[(l_n*84)+16]);
    __m128d a80_1 = _mm_load_sd(&A[313]);
#if defined(__SSE3__) && defined(__AVX__)
    c80_1 = _mm_add_sd(c80_1, _mm_mul_sd(a80_1, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c80_1 = _mm_add_sd(c80_1, _mm_mul_sd(a80_1, b80));
#endif
    _mm_store_sd(&C[(l_n*84)+16], c80_1);
    __m128d c80_2 = _mm_load_sd(&C[(l_n*84)+31]);
    __m128d a80_2 = _mm_load_sd(&A[314]);
#if defined(__SSE3__) && defined(__AVX__)
    c80_2 = _mm_add_sd(c80_2, _mm_mul_sd(a80_2, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c80_2 = _mm_add_sd(c80_2, _mm_mul_sd(a80_2, b80));
#endif
    _mm_store_sd(&C[(l_n*84)+31], c80_2);
    __m128d c80_3 = _mm_load_sd(&C[(l_n*84)+52]);
    __m128d a80_3 = _mm_load_sd(&A[315]);
#if defined(__SSE3__) && defined(__AVX__)
    c80_3 = _mm_add_sd(c80_3, _mm_mul_sd(a80_3, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c80_3 = _mm_add_sd(c80_3, _mm_mul_sd(a80_3, b80));
#endif
    _mm_store_sd(&C[(l_n*84)+52], c80_3);
    __m128d c80_4 = _mm_load_sd(&C[(l_n*84)+80]);
    __m128d a80_4 = _mm_load_sd(&A[316]);
#if defined(__SSE3__) && defined(__AVX__)
    c80_4 = _mm_add_sd(c80_4, _mm_mul_sd(a80_4, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c80_4 = _mm_add_sd(c80_4, _mm_mul_sd(a80_4, b80));
#endif
    _mm_store_sd(&C[(l_n*84)+80], c80_4);
#else
    C[(l_n*84)+6] += A[312] * B[(l_n*84)+80];
    C[(l_n*84)+16] += A[313] * B[(l_n*84)+80];
    C[(l_n*84)+31] += A[314] * B[(l_n*84)+80];
    C[(l_n*84)+52] += A[315] * B[(l_n*84)+80];
    C[(l_n*84)+80] += A[316] * B[(l_n*84)+80];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b81 = _mm256_broadcast_sd(&B[(l_n*84)+81]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b81 = _mm_loaddup_pd(&B[(l_n*84)+81]);
#endif
    __m128d c81_0 = _mm_load_sd(&C[(l_n*84)+1]);
    __m128d a81_0 = _mm_load_sd(&A[317]);
#if defined(__SSE3__) && defined(__AVX__)
    c81_0 = _mm_add_sd(c81_0, _mm_mul_sd(a81_0, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c81_0 = _mm_add_sd(c81_0, _mm_mul_sd(a81_0, b81));
#endif
    _mm_store_sd(&C[(l_n*84)+1], c81_0);
    __m128d c81_1 = _mm_load_sd(&C[(l_n*84)+7]);
    __m128d a81_1 = _mm_load_sd(&A[318]);
#if defined(__SSE3__) && defined(__AVX__)
    c81_1 = _mm_add_sd(c81_1, _mm_mul_sd(a81_1, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c81_1 = _mm_add_sd(c81_1, _mm_mul_sd(a81_1, b81));
#endif
    _mm_store_sd(&C[(l_n*84)+7], c81_1);
    __m128d c81_2 = _mm_load_sd(&C[(l_n*84)+17]);
    __m128d a81_2 = _mm_load_sd(&A[319]);
#if defined(__SSE3__) && defined(__AVX__)
    c81_2 = _mm_add_sd(c81_2, _mm_mul_sd(a81_2, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c81_2 = _mm_add_sd(c81_2, _mm_mul_sd(a81_2, b81));
#endif
    _mm_store_sd(&C[(l_n*84)+17], c81_2);
    __m128d c81_3 = _mm_load_sd(&C[(l_n*84)+32]);
    __m128d a81_3 = _mm_load_sd(&A[320]);
#if defined(__SSE3__) && defined(__AVX__)
    c81_3 = _mm_add_sd(c81_3, _mm_mul_sd(a81_3, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c81_3 = _mm_add_sd(c81_3, _mm_mul_sd(a81_3, b81));
#endif
    _mm_store_sd(&C[(l_n*84)+32], c81_3);
    __m128d c81_4 = _mm_load_sd(&C[(l_n*84)+53]);
    __m128d a81_4 = _mm_load_sd(&A[321]);
#if defined(__SSE3__) && defined(__AVX__)
    c81_4 = _mm_add_sd(c81_4, _mm_mul_sd(a81_4, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c81_4 = _mm_add_sd(c81_4, _mm_mul_sd(a81_4, b81));
#endif
    _mm_store_sd(&C[(l_n*84)+53], c81_4);
    __m128d c81_5 = _mm_load_sd(&C[(l_n*84)+81]);
    __m128d a81_5 = _mm_load_sd(&A[322]);
#if defined(__SSE3__) && defined(__AVX__)
    c81_5 = _mm_add_sd(c81_5, _mm_mul_sd(a81_5, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c81_5 = _mm_add_sd(c81_5, _mm_mul_sd(a81_5, b81));
#endif
    _mm_store_sd(&C[(l_n*84)+81], c81_5);
#else
    C[(l_n*84)+1] += A[317] * B[(l_n*84)+81];
    C[(l_n*84)+7] += A[318] * B[(l_n*84)+81];
    C[(l_n*84)+17] += A[319] * B[(l_n*84)+81];
    C[(l_n*84)+32] += A[320] * B[(l_n*84)+81];
    C[(l_n*84)+53] += A[321] * B[(l_n*84)+81];
    C[(l_n*84)+81] += A[322] * B[(l_n*84)+81];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b82 = _mm256_broadcast_sd(&B[(l_n*84)+82]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b82 = _mm_loaddup_pd(&B[(l_n*84)+82]);
#endif
    __m128d c82_0 = _mm_load_sd(&C[(l_n*84)+2]);
    __m128d a82_0 = _mm_load_sd(&A[323]);
#if defined(__SSE3__) && defined(__AVX__)
    c82_0 = _mm_add_sd(c82_0, _mm_mul_sd(a82_0, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c82_0 = _mm_add_sd(c82_0, _mm_mul_sd(a82_0, b82));
#endif
    _mm_store_sd(&C[(l_n*84)+2], c82_0);
    __m128d c82_1 = _mm_load_sd(&C[(l_n*84)+8]);
    __m128d a82_1 = _mm_load_sd(&A[324]);
#if defined(__SSE3__) && defined(__AVX__)
    c82_1 = _mm_add_sd(c82_1, _mm_mul_sd(a82_1, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c82_1 = _mm_add_sd(c82_1, _mm_mul_sd(a82_1, b82));
#endif
    _mm_store_sd(&C[(l_n*84)+8], c82_1);
    __m128d c82_2 = _mm_load_sd(&C[(l_n*84)+18]);
    __m128d a82_2 = _mm_load_sd(&A[325]);
#if defined(__SSE3__) && defined(__AVX__)
    c82_2 = _mm_add_sd(c82_2, _mm_mul_sd(a82_2, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c82_2 = _mm_add_sd(c82_2, _mm_mul_sd(a82_2, b82));
#endif
    _mm_store_sd(&C[(l_n*84)+18], c82_2);
    __m128d c82_3 = _mm_load_sd(&C[(l_n*84)+33]);
    __m128d a82_3 = _mm_load_sd(&A[326]);
#if defined(__SSE3__) && defined(__AVX__)
    c82_3 = _mm_add_sd(c82_3, _mm_mul_sd(a82_3, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c82_3 = _mm_add_sd(c82_3, _mm_mul_sd(a82_3, b82));
#endif
    _mm_store_sd(&C[(l_n*84)+33], c82_3);
    __m128d c82_4 = _mm_load_sd(&C[(l_n*84)+54]);
    __m128d a82_4 = _mm_load_sd(&A[327]);
#if defined(__SSE3__) && defined(__AVX__)
    c82_4 = _mm_add_sd(c82_4, _mm_mul_sd(a82_4, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c82_4 = _mm_add_sd(c82_4, _mm_mul_sd(a82_4, b82));
#endif
    _mm_store_sd(&C[(l_n*84)+54], c82_4);
    __m128d c82_5 = _mm_load_sd(&C[(l_n*84)+82]);
    __m128d a82_5 = _mm_load_sd(&A[328]);
#if defined(__SSE3__) && defined(__AVX__)
    c82_5 = _mm_add_sd(c82_5, _mm_mul_sd(a82_5, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c82_5 = _mm_add_sd(c82_5, _mm_mul_sd(a82_5, b82));
#endif
    _mm_store_sd(&C[(l_n*84)+82], c82_5);
#else
    C[(l_n*84)+2] += A[323] * B[(l_n*84)+82];
    C[(l_n*84)+8] += A[324] * B[(l_n*84)+82];
    C[(l_n*84)+18] += A[325] * B[(l_n*84)+82];
    C[(l_n*84)+33] += A[326] * B[(l_n*84)+82];
    C[(l_n*84)+54] += A[327] * B[(l_n*84)+82];
    C[(l_n*84)+82] += A[328] * B[(l_n*84)+82];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b83 = _mm256_broadcast_sd(&B[(l_n*84)+83]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b83 = _mm_loaddup_pd(&B[(l_n*84)+83]);
#endif
    __m128d c83_0 = _mm_load_sd(&C[(l_n*84)+0]);
    __m128d a83_0 = _mm_load_sd(&A[329]);
#if defined(__SSE3__) && defined(__AVX__)
    c83_0 = _mm_add_sd(c83_0, _mm_mul_sd(a83_0, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c83_0 = _mm_add_sd(c83_0, _mm_mul_sd(a83_0, b83));
#endif
    _mm_store_sd(&C[(l_n*84)+0], c83_0);
    __m128d c83_1 = _mm_load_sd(&C[(l_n*84)+3]);
    __m128d a83_1 = _mm_load_sd(&A[330]);
#if defined(__SSE3__) && defined(__AVX__)
    c83_1 = _mm_add_sd(c83_1, _mm_mul_sd(a83_1, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c83_1 = _mm_add_sd(c83_1, _mm_mul_sd(a83_1, b83));
#endif
    _mm_store_sd(&C[(l_n*84)+3], c83_1);
    __m128d c83_2 = _mm_load_sd(&C[(l_n*84)+9]);
    __m128d a83_2 = _mm_load_sd(&A[331]);
#if defined(__SSE3__) && defined(__AVX__)
    c83_2 = _mm_add_sd(c83_2, _mm_mul_sd(a83_2, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c83_2 = _mm_add_sd(c83_2, _mm_mul_sd(a83_2, b83));
#endif
    _mm_store_sd(&C[(l_n*84)+9], c83_2);
    __m128d c83_3 = _mm_load_sd(&C[(l_n*84)+19]);
    __m128d a83_3 = _mm_load_sd(&A[332]);
#if defined(__SSE3__) && defined(__AVX__)
    c83_3 = _mm_add_sd(c83_3, _mm_mul_sd(a83_3, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c83_3 = _mm_add_sd(c83_3, _mm_mul_sd(a83_3, b83));
#endif
    _mm_store_sd(&C[(l_n*84)+19], c83_3);
    __m128d c83_4 = _mm_load_sd(&C[(l_n*84)+34]);
    __m128d a83_4 = _mm_load_sd(&A[333]);
#if defined(__SSE3__) && defined(__AVX__)
    c83_4 = _mm_add_sd(c83_4, _mm_mul_sd(a83_4, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c83_4 = _mm_add_sd(c83_4, _mm_mul_sd(a83_4, b83));
#endif
    _mm_store_sd(&C[(l_n*84)+34], c83_4);
    __m128d c83_5 = _mm_load_sd(&C[(l_n*84)+55]);
    __m128d a83_5 = _mm_load_sd(&A[334]);
#if defined(__SSE3__) && defined(__AVX__)
    c83_5 = _mm_add_sd(c83_5, _mm_mul_sd(a83_5, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c83_5 = _mm_add_sd(c83_5, _mm_mul_sd(a83_5, b83));
#endif
    _mm_store_sd(&C[(l_n*84)+55], c83_5);
    __m128d c83_6 = _mm_load_sd(&C[(l_n*84)+83]);
    __m128d a83_6 = _mm_load_sd(&A[335]);
#if defined(__SSE3__) && defined(__AVX__)
    c83_6 = _mm_add_sd(c83_6, _mm_mul_sd(a83_6, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c83_6 = _mm_add_sd(c83_6, _mm_mul_sd(a83_6, b83));
#endif
    _mm_store_sd(&C[(l_n*84)+83], c83_6);
#else
    C[(l_n*84)+0] += A[329] * B[(l_n*84)+83];
    C[(l_n*84)+3] += A[330] * B[(l_n*84)+83];
    C[(l_n*84)+9] += A[331] * B[(l_n*84)+83];
    C[(l_n*84)+19] += A[332] * B[(l_n*84)+83];
    C[(l_n*84)+34] += A[333] * B[(l_n*84)+83];
    C[(l_n*84)+55] += A[334] * B[(l_n*84)+83];
    C[(l_n*84)+83] += A[335] * B[(l_n*84)+83];
#endif

  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 6048;
#endif
}

void dsparse_starMatrix_m120_n9_k9_ldA120_ldBna8_ldC120_beta1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
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
