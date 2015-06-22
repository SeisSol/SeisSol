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
// @date 2015-05-09 22:18:08.079116
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

void dsparse_kXiDivMT_m1_n9_k4_ldAna2_ldB8_ldC8_beta0_pfsigonly(const double* values, const double* B, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma nounroll_and_jam
for (unsigned int i = 0; i < 9; i++)
{
  for (unsigned int m = 0; m < 1; m++) {
    C[(i*8)+m] = 0.0;
  }
#if defined(__SSE3__) || defined(__AVX__)
#else
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b1 = _mm256_broadcast_sd(&B[(i*8)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b1 = _mm_loaddup_pd(&B[(i*8)+1]);
#endif
__m128d c1_0 = _mm_load_sd(&C[(i*8)+0]);
__m128d a1_0 = _mm_load_sd(&values[0]);
#if defined(__SSE3__) && defined(__AVX__)
c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, b1));
#endif
_mm_store_sd(&C[(i*8)+0], c1_0);
#else
C[(i*8)+0] += values[0] * B[(i*8)+1];
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

void dsparse_starMatrix_m1_n9_k9_ldA8_ldBna2_ldC8_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
for (unsigned int i = 0; i < 1; i += 1)
{
C[(i)+(0)] += A[(i)+(48)] * values[0];
C[(i)+(0)] += A[(i)+(56)] * values[1];
C[(i)+(0)] += A[(i)+(64)] * values[2];
C[(i)+(8)] += A[(i)+(48)] * values[3];
C[(i)+(8)] += A[(i)+(56)] * values[4];
C[(i)+(8)] += A[(i)+(64)] * values[5];
C[(i)+(16)] += A[(i)+(48)] * values[6];
C[(i)+(16)] += A[(i)+(56)] * values[7];
C[(i)+(16)] += A[(i)+(64)] * values[8];
C[(i)+(24)] += A[(i)+(48)] * values[9];
C[(i)+(24)] += A[(i)+(56)] * values[10];
C[(i)+(32)] += A[(i)+(56)] * values[11];
C[(i)+(32)] += A[(i)+(64)] * values[12];
C[(i)+(40)] += A[(i)+(48)] * values[13];
C[(i)+(40)] += A[(i)+(64)] * values[14];
C[(i)+(48)] += A[(i)+(0)] * values[15];
C[(i)+(48)] += A[(i)+(24)] * values[16];
C[(i)+(48)] += A[(i)+(40)] * values[17];
C[(i)+(56)] += A[(i)+(8)] * values[18];
C[(i)+(56)] += A[(i)+(24)] * values[19];
C[(i)+(56)] += A[(i)+(32)] * values[20];
C[(i)+(64)] += A[(i)+(16)] * values[21];
C[(i)+(64)] += A[(i)+(32)] * values[22];
C[(i)+(64)] += A[(i)+(40)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 48;
#endif

}

void dsparse_starMatrix_m4_n9_k9_ldA8_ldBna3_ldC8_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(32)
#pragma vector aligned
for (unsigned int i = 0; i < 4; i += 1)
{
C[(i)+(0)] += A[(i)+(48)] * values[0];
C[(i)+(0)] += A[(i)+(56)] * values[1];
C[(i)+(0)] += A[(i)+(64)] * values[2];
C[(i)+(8)] += A[(i)+(48)] * values[3];
C[(i)+(8)] += A[(i)+(56)] * values[4];
C[(i)+(8)] += A[(i)+(64)] * values[5];
C[(i)+(16)] += A[(i)+(48)] * values[6];
C[(i)+(16)] += A[(i)+(56)] * values[7];
C[(i)+(16)] += A[(i)+(64)] * values[8];
C[(i)+(24)] += A[(i)+(48)] * values[9];
C[(i)+(24)] += A[(i)+(56)] * values[10];
C[(i)+(32)] += A[(i)+(56)] * values[11];
C[(i)+(32)] += A[(i)+(64)] * values[12];
C[(i)+(40)] += A[(i)+(48)] * values[13];
C[(i)+(40)] += A[(i)+(64)] * values[14];
C[(i)+(48)] += A[(i)+(0)] * values[15];
C[(i)+(48)] += A[(i)+(24)] * values[16];
C[(i)+(48)] += A[(i)+(40)] * values[17];
C[(i)+(56)] += A[(i)+(8)] * values[18];
C[(i)+(56)] += A[(i)+(24)] * values[19];
C[(i)+(56)] += A[(i)+(32)] * values[20];
C[(i)+(64)] += A[(i)+(16)] * values[21];
C[(i)+(64)] += A[(i)+(32)] * values[22];
C[(i)+(64)] += A[(i)+(40)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 192;
#endif

}

void dsparse_starMatrix_m1_n9_k9_ldA8_ldBna3_ldC8_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
for (unsigned int i = 0; i < 1; i += 1)
{
C[(i)+(0)] += A[(i)+(48)] * values[0];
C[(i)+(0)] += A[(i)+(56)] * values[1];
C[(i)+(0)] += A[(i)+(64)] * values[2];
C[(i)+(8)] += A[(i)+(48)] * values[3];
C[(i)+(8)] += A[(i)+(56)] * values[4];
C[(i)+(8)] += A[(i)+(64)] * values[5];
C[(i)+(16)] += A[(i)+(48)] * values[6];
C[(i)+(16)] += A[(i)+(56)] * values[7];
C[(i)+(16)] += A[(i)+(64)] * values[8];
C[(i)+(24)] += A[(i)+(48)] * values[9];
C[(i)+(24)] += A[(i)+(56)] * values[10];
C[(i)+(32)] += A[(i)+(56)] * values[11];
C[(i)+(32)] += A[(i)+(64)] * values[12];
C[(i)+(40)] += A[(i)+(48)] * values[13];
C[(i)+(40)] += A[(i)+(64)] * values[14];
C[(i)+(48)] += A[(i)+(0)] * values[15];
C[(i)+(48)] += A[(i)+(24)] * values[16];
C[(i)+(48)] += A[(i)+(40)] * values[17];
C[(i)+(56)] += A[(i)+(8)] * values[18];
C[(i)+(56)] += A[(i)+(24)] * values[19];
C[(i)+(56)] += A[(i)+(32)] * values[20];
C[(i)+(64)] += A[(i)+(16)] * values[21];
C[(i)+(64)] += A[(i)+(32)] * values[22];
C[(i)+(64)] += A[(i)+(40)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 48;
#endif

}

void dsparse_starMatrix_m10_n9_k9_ldA16_ldBna4_ldC16_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(32)
#pragma vector aligned
for (unsigned int i = 0; i < 10; i += 1)
{
C[(i)+(0)] += A[(i)+(96)] * values[0];
C[(i)+(0)] += A[(i)+(112)] * values[1];
C[(i)+(0)] += A[(i)+(128)] * values[2];
C[(i)+(16)] += A[(i)+(96)] * values[3];
C[(i)+(16)] += A[(i)+(112)] * values[4];
C[(i)+(16)] += A[(i)+(128)] * values[5];
C[(i)+(32)] += A[(i)+(96)] * values[6];
C[(i)+(32)] += A[(i)+(112)] * values[7];
C[(i)+(32)] += A[(i)+(128)] * values[8];
C[(i)+(48)] += A[(i)+(96)] * values[9];
C[(i)+(48)] += A[(i)+(112)] * values[10];
C[(i)+(64)] += A[(i)+(112)] * values[11];
C[(i)+(64)] += A[(i)+(128)] * values[12];
C[(i)+(80)] += A[(i)+(96)] * values[13];
C[(i)+(80)] += A[(i)+(128)] * values[14];
C[(i)+(96)] += A[(i)+(0)] * values[15];
C[(i)+(96)] += A[(i)+(48)] * values[16];
C[(i)+(96)] += A[(i)+(80)] * values[17];
C[(i)+(112)] += A[(i)+(16)] * values[18];
C[(i)+(112)] += A[(i)+(48)] * values[19];
C[(i)+(112)] += A[(i)+(64)] * values[20];
C[(i)+(128)] += A[(i)+(32)] * values[21];
C[(i)+(128)] += A[(i)+(64)] * values[22];
C[(i)+(128)] += A[(i)+(80)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 480;
#endif

}

void dsparse_starMatrix_m4_n9_k9_ldA8_ldBna4_ldC8_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(32)
#pragma vector aligned
for (unsigned int i = 0; i < 4; i += 1)
{
C[(i)+(0)] += A[(i)+(48)] * values[0];
C[(i)+(0)] += A[(i)+(56)] * values[1];
C[(i)+(0)] += A[(i)+(64)] * values[2];
C[(i)+(8)] += A[(i)+(48)] * values[3];
C[(i)+(8)] += A[(i)+(56)] * values[4];
C[(i)+(8)] += A[(i)+(64)] * values[5];
C[(i)+(16)] += A[(i)+(48)] * values[6];
C[(i)+(16)] += A[(i)+(56)] * values[7];
C[(i)+(16)] += A[(i)+(64)] * values[8];
C[(i)+(24)] += A[(i)+(48)] * values[9];
C[(i)+(24)] += A[(i)+(56)] * values[10];
C[(i)+(32)] += A[(i)+(56)] * values[11];
C[(i)+(32)] += A[(i)+(64)] * values[12];
C[(i)+(40)] += A[(i)+(48)] * values[13];
C[(i)+(40)] += A[(i)+(64)] * values[14];
C[(i)+(48)] += A[(i)+(0)] * values[15];
C[(i)+(48)] += A[(i)+(24)] * values[16];
C[(i)+(48)] += A[(i)+(40)] * values[17];
C[(i)+(56)] += A[(i)+(8)] * values[18];
C[(i)+(56)] += A[(i)+(24)] * values[19];
C[(i)+(56)] += A[(i)+(32)] * values[20];
C[(i)+(64)] += A[(i)+(16)] * values[21];
C[(i)+(64)] += A[(i)+(32)] * values[22];
C[(i)+(64)] += A[(i)+(40)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 192;
#endif

}

void dsparse_starMatrix_m1_n9_k9_ldA8_ldBna4_ldC8_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
for (unsigned int i = 0; i < 1; i += 1)
{
C[(i)+(0)] += A[(i)+(48)] * values[0];
C[(i)+(0)] += A[(i)+(56)] * values[1];
C[(i)+(0)] += A[(i)+(64)] * values[2];
C[(i)+(8)] += A[(i)+(48)] * values[3];
C[(i)+(8)] += A[(i)+(56)] * values[4];
C[(i)+(8)] += A[(i)+(64)] * values[5];
C[(i)+(16)] += A[(i)+(48)] * values[6];
C[(i)+(16)] += A[(i)+(56)] * values[7];
C[(i)+(16)] += A[(i)+(64)] * values[8];
C[(i)+(24)] += A[(i)+(48)] * values[9];
C[(i)+(24)] += A[(i)+(56)] * values[10];
C[(i)+(32)] += A[(i)+(56)] * values[11];
C[(i)+(32)] += A[(i)+(64)] * values[12];
C[(i)+(40)] += A[(i)+(48)] * values[13];
C[(i)+(40)] += A[(i)+(64)] * values[14];
C[(i)+(48)] += A[(i)+(0)] * values[15];
C[(i)+(48)] += A[(i)+(24)] * values[16];
C[(i)+(48)] += A[(i)+(40)] * values[17];
C[(i)+(56)] += A[(i)+(8)] * values[18];
C[(i)+(56)] += A[(i)+(24)] * values[19];
C[(i)+(56)] += A[(i)+(32)] * values[20];
C[(i)+(64)] += A[(i)+(16)] * values[21];
C[(i)+(64)] += A[(i)+(32)] * values[22];
C[(i)+(64)] += A[(i)+(40)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 48;
#endif

}

void dsparse_starMatrix_m20_n9_k9_ldA24_ldBna5_ldC24_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(32)
#pragma vector aligned
for (unsigned int i = 0; i < 20; i += 1)
{
C[(i)+(0)] += A[(i)+(144)] * values[0];
C[(i)+(0)] += A[(i)+(168)] * values[1];
C[(i)+(0)] += A[(i)+(192)] * values[2];
C[(i)+(24)] += A[(i)+(144)] * values[3];
C[(i)+(24)] += A[(i)+(168)] * values[4];
C[(i)+(24)] += A[(i)+(192)] * values[5];
C[(i)+(48)] += A[(i)+(144)] * values[6];
C[(i)+(48)] += A[(i)+(168)] * values[7];
C[(i)+(48)] += A[(i)+(192)] * values[8];
C[(i)+(72)] += A[(i)+(144)] * values[9];
C[(i)+(72)] += A[(i)+(168)] * values[10];
C[(i)+(96)] += A[(i)+(168)] * values[11];
C[(i)+(96)] += A[(i)+(192)] * values[12];
C[(i)+(120)] += A[(i)+(144)] * values[13];
C[(i)+(120)] += A[(i)+(192)] * values[14];
C[(i)+(144)] += A[(i)+(0)] * values[15];
C[(i)+(144)] += A[(i)+(72)] * values[16];
C[(i)+(144)] += A[(i)+(120)] * values[17];
C[(i)+(168)] += A[(i)+(24)] * values[18];
C[(i)+(168)] += A[(i)+(72)] * values[19];
C[(i)+(168)] += A[(i)+(96)] * values[20];
C[(i)+(192)] += A[(i)+(48)] * values[21];
C[(i)+(192)] += A[(i)+(96)] * values[22];
C[(i)+(192)] += A[(i)+(120)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 960;
#endif

}

void dsparse_starMatrix_m10_n9_k9_ldA16_ldBna5_ldC16_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(32)
#pragma vector aligned
for (unsigned int i = 0; i < 10; i += 1)
{
C[(i)+(0)] += A[(i)+(96)] * values[0];
C[(i)+(0)] += A[(i)+(112)] * values[1];
C[(i)+(0)] += A[(i)+(128)] * values[2];
C[(i)+(16)] += A[(i)+(96)] * values[3];
C[(i)+(16)] += A[(i)+(112)] * values[4];
C[(i)+(16)] += A[(i)+(128)] * values[5];
C[(i)+(32)] += A[(i)+(96)] * values[6];
C[(i)+(32)] += A[(i)+(112)] * values[7];
C[(i)+(32)] += A[(i)+(128)] * values[8];
C[(i)+(48)] += A[(i)+(96)] * values[9];
C[(i)+(48)] += A[(i)+(112)] * values[10];
C[(i)+(64)] += A[(i)+(112)] * values[11];
C[(i)+(64)] += A[(i)+(128)] * values[12];
C[(i)+(80)] += A[(i)+(96)] * values[13];
C[(i)+(80)] += A[(i)+(128)] * values[14];
C[(i)+(96)] += A[(i)+(0)] * values[15];
C[(i)+(96)] += A[(i)+(48)] * values[16];
C[(i)+(96)] += A[(i)+(80)] * values[17];
C[(i)+(112)] += A[(i)+(16)] * values[18];
C[(i)+(112)] += A[(i)+(48)] * values[19];
C[(i)+(112)] += A[(i)+(64)] * values[20];
C[(i)+(128)] += A[(i)+(32)] * values[21];
C[(i)+(128)] += A[(i)+(64)] * values[22];
C[(i)+(128)] += A[(i)+(80)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 480;
#endif

}

void dsparse_starMatrix_m4_n9_k9_ldA8_ldBna5_ldC8_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(32)
#pragma vector aligned
for (unsigned int i = 0; i < 4; i += 1)
{
C[(i)+(0)] += A[(i)+(48)] * values[0];
C[(i)+(0)] += A[(i)+(56)] * values[1];
C[(i)+(0)] += A[(i)+(64)] * values[2];
C[(i)+(8)] += A[(i)+(48)] * values[3];
C[(i)+(8)] += A[(i)+(56)] * values[4];
C[(i)+(8)] += A[(i)+(64)] * values[5];
C[(i)+(16)] += A[(i)+(48)] * values[6];
C[(i)+(16)] += A[(i)+(56)] * values[7];
C[(i)+(16)] += A[(i)+(64)] * values[8];
C[(i)+(24)] += A[(i)+(48)] * values[9];
C[(i)+(24)] += A[(i)+(56)] * values[10];
C[(i)+(32)] += A[(i)+(56)] * values[11];
C[(i)+(32)] += A[(i)+(64)] * values[12];
C[(i)+(40)] += A[(i)+(48)] * values[13];
C[(i)+(40)] += A[(i)+(64)] * values[14];
C[(i)+(48)] += A[(i)+(0)] * values[15];
C[(i)+(48)] += A[(i)+(24)] * values[16];
C[(i)+(48)] += A[(i)+(40)] * values[17];
C[(i)+(56)] += A[(i)+(8)] * values[18];
C[(i)+(56)] += A[(i)+(24)] * values[19];
C[(i)+(56)] += A[(i)+(32)] * values[20];
C[(i)+(64)] += A[(i)+(16)] * values[21];
C[(i)+(64)] += A[(i)+(32)] * values[22];
C[(i)+(64)] += A[(i)+(40)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 192;
#endif

}

void dsparse_starMatrix_m1_n9_k9_ldA8_ldBna5_ldC8_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
for (unsigned int i = 0; i < 1; i += 1)
{
C[(i)+(0)] += A[(i)+(48)] * values[0];
C[(i)+(0)] += A[(i)+(56)] * values[1];
C[(i)+(0)] += A[(i)+(64)] * values[2];
C[(i)+(8)] += A[(i)+(48)] * values[3];
C[(i)+(8)] += A[(i)+(56)] * values[4];
C[(i)+(8)] += A[(i)+(64)] * values[5];
C[(i)+(16)] += A[(i)+(48)] * values[6];
C[(i)+(16)] += A[(i)+(56)] * values[7];
C[(i)+(16)] += A[(i)+(64)] * values[8];
C[(i)+(24)] += A[(i)+(48)] * values[9];
C[(i)+(24)] += A[(i)+(56)] * values[10];
C[(i)+(32)] += A[(i)+(56)] * values[11];
C[(i)+(32)] += A[(i)+(64)] * values[12];
C[(i)+(40)] += A[(i)+(48)] * values[13];
C[(i)+(40)] += A[(i)+(64)] * values[14];
C[(i)+(48)] += A[(i)+(0)] * values[15];
C[(i)+(48)] += A[(i)+(24)] * values[16];
C[(i)+(48)] += A[(i)+(40)] * values[17];
C[(i)+(56)] += A[(i)+(8)] * values[18];
C[(i)+(56)] += A[(i)+(24)] * values[19];
C[(i)+(56)] += A[(i)+(32)] * values[20];
C[(i)+(64)] += A[(i)+(16)] * values[21];
C[(i)+(64)] += A[(i)+(32)] * values[22];
C[(i)+(64)] += A[(i)+(40)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 48;
#endif

}

void dsparse_starMatrix_m35_n9_k9_ldA40_ldBna6_ldC40_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(32)
#pragma vector aligned
for (unsigned int i = 0; i < 35; i += 1)
{
C[(i)+(0)] += A[(i)+(240)] * values[0];
C[(i)+(0)] += A[(i)+(280)] * values[1];
C[(i)+(0)] += A[(i)+(320)] * values[2];
C[(i)+(40)] += A[(i)+(240)] * values[3];
C[(i)+(40)] += A[(i)+(280)] * values[4];
C[(i)+(40)] += A[(i)+(320)] * values[5];
C[(i)+(80)] += A[(i)+(240)] * values[6];
C[(i)+(80)] += A[(i)+(280)] * values[7];
C[(i)+(80)] += A[(i)+(320)] * values[8];
C[(i)+(120)] += A[(i)+(240)] * values[9];
C[(i)+(120)] += A[(i)+(280)] * values[10];
C[(i)+(160)] += A[(i)+(280)] * values[11];
C[(i)+(160)] += A[(i)+(320)] * values[12];
C[(i)+(200)] += A[(i)+(240)] * values[13];
C[(i)+(200)] += A[(i)+(320)] * values[14];
C[(i)+(240)] += A[(i)+(0)] * values[15];
C[(i)+(240)] += A[(i)+(120)] * values[16];
C[(i)+(240)] += A[(i)+(200)] * values[17];
C[(i)+(280)] += A[(i)+(40)] * values[18];
C[(i)+(280)] += A[(i)+(120)] * values[19];
C[(i)+(280)] += A[(i)+(160)] * values[20];
C[(i)+(320)] += A[(i)+(80)] * values[21];
C[(i)+(320)] += A[(i)+(160)] * values[22];
C[(i)+(320)] += A[(i)+(200)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 1680;
#endif

}

void dsparse_starMatrix_m20_n9_k9_ldA24_ldBna6_ldC24_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(32)
#pragma vector aligned
for (unsigned int i = 0; i < 20; i += 1)
{
C[(i)+(0)] += A[(i)+(144)] * values[0];
C[(i)+(0)] += A[(i)+(168)] * values[1];
C[(i)+(0)] += A[(i)+(192)] * values[2];
C[(i)+(24)] += A[(i)+(144)] * values[3];
C[(i)+(24)] += A[(i)+(168)] * values[4];
C[(i)+(24)] += A[(i)+(192)] * values[5];
C[(i)+(48)] += A[(i)+(144)] * values[6];
C[(i)+(48)] += A[(i)+(168)] * values[7];
C[(i)+(48)] += A[(i)+(192)] * values[8];
C[(i)+(72)] += A[(i)+(144)] * values[9];
C[(i)+(72)] += A[(i)+(168)] * values[10];
C[(i)+(96)] += A[(i)+(168)] * values[11];
C[(i)+(96)] += A[(i)+(192)] * values[12];
C[(i)+(120)] += A[(i)+(144)] * values[13];
C[(i)+(120)] += A[(i)+(192)] * values[14];
C[(i)+(144)] += A[(i)+(0)] * values[15];
C[(i)+(144)] += A[(i)+(72)] * values[16];
C[(i)+(144)] += A[(i)+(120)] * values[17];
C[(i)+(168)] += A[(i)+(24)] * values[18];
C[(i)+(168)] += A[(i)+(72)] * values[19];
C[(i)+(168)] += A[(i)+(96)] * values[20];
C[(i)+(192)] += A[(i)+(48)] * values[21];
C[(i)+(192)] += A[(i)+(96)] * values[22];
C[(i)+(192)] += A[(i)+(120)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 960;
#endif

}

void dsparse_starMatrix_m10_n9_k9_ldA16_ldBna6_ldC16_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(32)
#pragma vector aligned
for (unsigned int i = 0; i < 10; i += 1)
{
C[(i)+(0)] += A[(i)+(96)] * values[0];
C[(i)+(0)] += A[(i)+(112)] * values[1];
C[(i)+(0)] += A[(i)+(128)] * values[2];
C[(i)+(16)] += A[(i)+(96)] * values[3];
C[(i)+(16)] += A[(i)+(112)] * values[4];
C[(i)+(16)] += A[(i)+(128)] * values[5];
C[(i)+(32)] += A[(i)+(96)] * values[6];
C[(i)+(32)] += A[(i)+(112)] * values[7];
C[(i)+(32)] += A[(i)+(128)] * values[8];
C[(i)+(48)] += A[(i)+(96)] * values[9];
C[(i)+(48)] += A[(i)+(112)] * values[10];
C[(i)+(64)] += A[(i)+(112)] * values[11];
C[(i)+(64)] += A[(i)+(128)] * values[12];
C[(i)+(80)] += A[(i)+(96)] * values[13];
C[(i)+(80)] += A[(i)+(128)] * values[14];
C[(i)+(96)] += A[(i)+(0)] * values[15];
C[(i)+(96)] += A[(i)+(48)] * values[16];
C[(i)+(96)] += A[(i)+(80)] * values[17];
C[(i)+(112)] += A[(i)+(16)] * values[18];
C[(i)+(112)] += A[(i)+(48)] * values[19];
C[(i)+(112)] += A[(i)+(64)] * values[20];
C[(i)+(128)] += A[(i)+(32)] * values[21];
C[(i)+(128)] += A[(i)+(64)] * values[22];
C[(i)+(128)] += A[(i)+(80)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 480;
#endif

}

void dsparse_starMatrix_m4_n9_k9_ldA8_ldBna6_ldC8_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(32)
#pragma vector aligned
for (unsigned int i = 0; i < 4; i += 1)
{
C[(i)+(0)] += A[(i)+(48)] * values[0];
C[(i)+(0)] += A[(i)+(56)] * values[1];
C[(i)+(0)] += A[(i)+(64)] * values[2];
C[(i)+(8)] += A[(i)+(48)] * values[3];
C[(i)+(8)] += A[(i)+(56)] * values[4];
C[(i)+(8)] += A[(i)+(64)] * values[5];
C[(i)+(16)] += A[(i)+(48)] * values[6];
C[(i)+(16)] += A[(i)+(56)] * values[7];
C[(i)+(16)] += A[(i)+(64)] * values[8];
C[(i)+(24)] += A[(i)+(48)] * values[9];
C[(i)+(24)] += A[(i)+(56)] * values[10];
C[(i)+(32)] += A[(i)+(56)] * values[11];
C[(i)+(32)] += A[(i)+(64)] * values[12];
C[(i)+(40)] += A[(i)+(48)] * values[13];
C[(i)+(40)] += A[(i)+(64)] * values[14];
C[(i)+(48)] += A[(i)+(0)] * values[15];
C[(i)+(48)] += A[(i)+(24)] * values[16];
C[(i)+(48)] += A[(i)+(40)] * values[17];
C[(i)+(56)] += A[(i)+(8)] * values[18];
C[(i)+(56)] += A[(i)+(24)] * values[19];
C[(i)+(56)] += A[(i)+(32)] * values[20];
C[(i)+(64)] += A[(i)+(16)] * values[21];
C[(i)+(64)] += A[(i)+(32)] * values[22];
C[(i)+(64)] += A[(i)+(40)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 192;
#endif

}

void dsparse_starMatrix_m1_n9_k9_ldA8_ldBna6_ldC8_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
for (unsigned int i = 0; i < 1; i += 1)
{
C[(i)+(0)] += A[(i)+(48)] * values[0];
C[(i)+(0)] += A[(i)+(56)] * values[1];
C[(i)+(0)] += A[(i)+(64)] * values[2];
C[(i)+(8)] += A[(i)+(48)] * values[3];
C[(i)+(8)] += A[(i)+(56)] * values[4];
C[(i)+(8)] += A[(i)+(64)] * values[5];
C[(i)+(16)] += A[(i)+(48)] * values[6];
C[(i)+(16)] += A[(i)+(56)] * values[7];
C[(i)+(16)] += A[(i)+(64)] * values[8];
C[(i)+(24)] += A[(i)+(48)] * values[9];
C[(i)+(24)] += A[(i)+(56)] * values[10];
C[(i)+(32)] += A[(i)+(56)] * values[11];
C[(i)+(32)] += A[(i)+(64)] * values[12];
C[(i)+(40)] += A[(i)+(48)] * values[13];
C[(i)+(40)] += A[(i)+(64)] * values[14];
C[(i)+(48)] += A[(i)+(0)] * values[15];
C[(i)+(48)] += A[(i)+(24)] * values[16];
C[(i)+(48)] += A[(i)+(40)] * values[17];
C[(i)+(56)] += A[(i)+(8)] * values[18];
C[(i)+(56)] += A[(i)+(24)] * values[19];
C[(i)+(56)] += A[(i)+(32)] * values[20];
C[(i)+(64)] += A[(i)+(16)] * values[21];
C[(i)+(64)] += A[(i)+(32)] * values[22];
C[(i)+(64)] += A[(i)+(40)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 48;
#endif

}

void dsparse_starMatrix_m56_n9_k9_ldA56_ldBna7_ldC56_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(32)
#pragma vector aligned
for (unsigned int i = 0; i < 56; i += 1)
{
C[(i)+(0)] += A[(i)+(336)] * values[0];
C[(i)+(0)] += A[(i)+(392)] * values[1];
C[(i)+(0)] += A[(i)+(448)] * values[2];
C[(i)+(56)] += A[(i)+(336)] * values[3];
C[(i)+(56)] += A[(i)+(392)] * values[4];
C[(i)+(56)] += A[(i)+(448)] * values[5];
C[(i)+(112)] += A[(i)+(336)] * values[6];
C[(i)+(112)] += A[(i)+(392)] * values[7];
C[(i)+(112)] += A[(i)+(448)] * values[8];
C[(i)+(168)] += A[(i)+(336)] * values[9];
C[(i)+(168)] += A[(i)+(392)] * values[10];
C[(i)+(224)] += A[(i)+(392)] * values[11];
C[(i)+(224)] += A[(i)+(448)] * values[12];
C[(i)+(280)] += A[(i)+(336)] * values[13];
C[(i)+(280)] += A[(i)+(448)] * values[14];
C[(i)+(336)] += A[(i)+(0)] * values[15];
C[(i)+(336)] += A[(i)+(168)] * values[16];
C[(i)+(336)] += A[(i)+(280)] * values[17];
C[(i)+(392)] += A[(i)+(56)] * values[18];
C[(i)+(392)] += A[(i)+(168)] * values[19];
C[(i)+(392)] += A[(i)+(224)] * values[20];
C[(i)+(448)] += A[(i)+(112)] * values[21];
C[(i)+(448)] += A[(i)+(224)] * values[22];
C[(i)+(448)] += A[(i)+(280)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 2688;
#endif

}

void dsparse_starMatrix_m35_n9_k9_ldA40_ldBna7_ldC40_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(32)
#pragma vector aligned
for (unsigned int i = 0; i < 35; i += 1)
{
C[(i)+(0)] += A[(i)+(240)] * values[0];
C[(i)+(0)] += A[(i)+(280)] * values[1];
C[(i)+(0)] += A[(i)+(320)] * values[2];
C[(i)+(40)] += A[(i)+(240)] * values[3];
C[(i)+(40)] += A[(i)+(280)] * values[4];
C[(i)+(40)] += A[(i)+(320)] * values[5];
C[(i)+(80)] += A[(i)+(240)] * values[6];
C[(i)+(80)] += A[(i)+(280)] * values[7];
C[(i)+(80)] += A[(i)+(320)] * values[8];
C[(i)+(120)] += A[(i)+(240)] * values[9];
C[(i)+(120)] += A[(i)+(280)] * values[10];
C[(i)+(160)] += A[(i)+(280)] * values[11];
C[(i)+(160)] += A[(i)+(320)] * values[12];
C[(i)+(200)] += A[(i)+(240)] * values[13];
C[(i)+(200)] += A[(i)+(320)] * values[14];
C[(i)+(240)] += A[(i)+(0)] * values[15];
C[(i)+(240)] += A[(i)+(120)] * values[16];
C[(i)+(240)] += A[(i)+(200)] * values[17];
C[(i)+(280)] += A[(i)+(40)] * values[18];
C[(i)+(280)] += A[(i)+(120)] * values[19];
C[(i)+(280)] += A[(i)+(160)] * values[20];
C[(i)+(320)] += A[(i)+(80)] * values[21];
C[(i)+(320)] += A[(i)+(160)] * values[22];
C[(i)+(320)] += A[(i)+(200)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 1680;
#endif

}

void dsparse_starMatrix_m20_n9_k9_ldA24_ldBna7_ldC24_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(32)
#pragma vector aligned
for (unsigned int i = 0; i < 20; i += 1)
{
C[(i)+(0)] += A[(i)+(144)] * values[0];
C[(i)+(0)] += A[(i)+(168)] * values[1];
C[(i)+(0)] += A[(i)+(192)] * values[2];
C[(i)+(24)] += A[(i)+(144)] * values[3];
C[(i)+(24)] += A[(i)+(168)] * values[4];
C[(i)+(24)] += A[(i)+(192)] * values[5];
C[(i)+(48)] += A[(i)+(144)] * values[6];
C[(i)+(48)] += A[(i)+(168)] * values[7];
C[(i)+(48)] += A[(i)+(192)] * values[8];
C[(i)+(72)] += A[(i)+(144)] * values[9];
C[(i)+(72)] += A[(i)+(168)] * values[10];
C[(i)+(96)] += A[(i)+(168)] * values[11];
C[(i)+(96)] += A[(i)+(192)] * values[12];
C[(i)+(120)] += A[(i)+(144)] * values[13];
C[(i)+(120)] += A[(i)+(192)] * values[14];
C[(i)+(144)] += A[(i)+(0)] * values[15];
C[(i)+(144)] += A[(i)+(72)] * values[16];
C[(i)+(144)] += A[(i)+(120)] * values[17];
C[(i)+(168)] += A[(i)+(24)] * values[18];
C[(i)+(168)] += A[(i)+(72)] * values[19];
C[(i)+(168)] += A[(i)+(96)] * values[20];
C[(i)+(192)] += A[(i)+(48)] * values[21];
C[(i)+(192)] += A[(i)+(96)] * values[22];
C[(i)+(192)] += A[(i)+(120)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 960;
#endif

}

void dsparse_starMatrix_m10_n9_k9_ldA16_ldBna7_ldC16_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(32)
#pragma vector aligned
for (unsigned int i = 0; i < 10; i += 1)
{
C[(i)+(0)] += A[(i)+(96)] * values[0];
C[(i)+(0)] += A[(i)+(112)] * values[1];
C[(i)+(0)] += A[(i)+(128)] * values[2];
C[(i)+(16)] += A[(i)+(96)] * values[3];
C[(i)+(16)] += A[(i)+(112)] * values[4];
C[(i)+(16)] += A[(i)+(128)] * values[5];
C[(i)+(32)] += A[(i)+(96)] * values[6];
C[(i)+(32)] += A[(i)+(112)] * values[7];
C[(i)+(32)] += A[(i)+(128)] * values[8];
C[(i)+(48)] += A[(i)+(96)] * values[9];
C[(i)+(48)] += A[(i)+(112)] * values[10];
C[(i)+(64)] += A[(i)+(112)] * values[11];
C[(i)+(64)] += A[(i)+(128)] * values[12];
C[(i)+(80)] += A[(i)+(96)] * values[13];
C[(i)+(80)] += A[(i)+(128)] * values[14];
C[(i)+(96)] += A[(i)+(0)] * values[15];
C[(i)+(96)] += A[(i)+(48)] * values[16];
C[(i)+(96)] += A[(i)+(80)] * values[17];
C[(i)+(112)] += A[(i)+(16)] * values[18];
C[(i)+(112)] += A[(i)+(48)] * values[19];
C[(i)+(112)] += A[(i)+(64)] * values[20];
C[(i)+(128)] += A[(i)+(32)] * values[21];
C[(i)+(128)] += A[(i)+(64)] * values[22];
C[(i)+(128)] += A[(i)+(80)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 480;
#endif

}

void dsparse_starMatrix_m4_n9_k9_ldA8_ldBna7_ldC8_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(32)
#pragma vector aligned
for (unsigned int i = 0; i < 4; i += 1)
{
C[(i)+(0)] += A[(i)+(48)] * values[0];
C[(i)+(0)] += A[(i)+(56)] * values[1];
C[(i)+(0)] += A[(i)+(64)] * values[2];
C[(i)+(8)] += A[(i)+(48)] * values[3];
C[(i)+(8)] += A[(i)+(56)] * values[4];
C[(i)+(8)] += A[(i)+(64)] * values[5];
C[(i)+(16)] += A[(i)+(48)] * values[6];
C[(i)+(16)] += A[(i)+(56)] * values[7];
C[(i)+(16)] += A[(i)+(64)] * values[8];
C[(i)+(24)] += A[(i)+(48)] * values[9];
C[(i)+(24)] += A[(i)+(56)] * values[10];
C[(i)+(32)] += A[(i)+(56)] * values[11];
C[(i)+(32)] += A[(i)+(64)] * values[12];
C[(i)+(40)] += A[(i)+(48)] * values[13];
C[(i)+(40)] += A[(i)+(64)] * values[14];
C[(i)+(48)] += A[(i)+(0)] * values[15];
C[(i)+(48)] += A[(i)+(24)] * values[16];
C[(i)+(48)] += A[(i)+(40)] * values[17];
C[(i)+(56)] += A[(i)+(8)] * values[18];
C[(i)+(56)] += A[(i)+(24)] * values[19];
C[(i)+(56)] += A[(i)+(32)] * values[20];
C[(i)+(64)] += A[(i)+(16)] * values[21];
C[(i)+(64)] += A[(i)+(32)] * values[22];
C[(i)+(64)] += A[(i)+(40)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 192;
#endif

}

void dsparse_starMatrix_m1_n9_k9_ldA8_ldBna7_ldC8_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
for (unsigned int i = 0; i < 1; i += 1)
{
C[(i)+(0)] += A[(i)+(48)] * values[0];
C[(i)+(0)] += A[(i)+(56)] * values[1];
C[(i)+(0)] += A[(i)+(64)] * values[2];
C[(i)+(8)] += A[(i)+(48)] * values[3];
C[(i)+(8)] += A[(i)+(56)] * values[4];
C[(i)+(8)] += A[(i)+(64)] * values[5];
C[(i)+(16)] += A[(i)+(48)] * values[6];
C[(i)+(16)] += A[(i)+(56)] * values[7];
C[(i)+(16)] += A[(i)+(64)] * values[8];
C[(i)+(24)] += A[(i)+(48)] * values[9];
C[(i)+(24)] += A[(i)+(56)] * values[10];
C[(i)+(32)] += A[(i)+(56)] * values[11];
C[(i)+(32)] += A[(i)+(64)] * values[12];
C[(i)+(40)] += A[(i)+(48)] * values[13];
C[(i)+(40)] += A[(i)+(64)] * values[14];
C[(i)+(48)] += A[(i)+(0)] * values[15];
C[(i)+(48)] += A[(i)+(24)] * values[16];
C[(i)+(48)] += A[(i)+(40)] * values[17];
C[(i)+(56)] += A[(i)+(8)] * values[18];
C[(i)+(56)] += A[(i)+(24)] * values[19];
C[(i)+(56)] += A[(i)+(32)] * values[20];
C[(i)+(64)] += A[(i)+(16)] * values[21];
C[(i)+(64)] += A[(i)+(32)] * values[22];
C[(i)+(64)] += A[(i)+(40)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 48;
#endif

}

void dsparse_starMatrix_m4_n9_k9_ldA8_ldBna2_ldC8_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(32)
#pragma vector aligned
for (unsigned int i = 0; i < 4; i += 1)
{
C[(i)+(0)] += A[(i)+(48)] * values[0];
C[(i)+(0)] += A[(i)+(56)] * values[1];
C[(i)+(0)] += A[(i)+(64)] * values[2];
C[(i)+(8)] += A[(i)+(48)] * values[3];
C[(i)+(8)] += A[(i)+(56)] * values[4];
C[(i)+(8)] += A[(i)+(64)] * values[5];
C[(i)+(16)] += A[(i)+(48)] * values[6];
C[(i)+(16)] += A[(i)+(56)] * values[7];
C[(i)+(16)] += A[(i)+(64)] * values[8];
C[(i)+(24)] += A[(i)+(48)] * values[9];
C[(i)+(24)] += A[(i)+(56)] * values[10];
C[(i)+(32)] += A[(i)+(56)] * values[11];
C[(i)+(32)] += A[(i)+(64)] * values[12];
C[(i)+(40)] += A[(i)+(48)] * values[13];
C[(i)+(40)] += A[(i)+(64)] * values[14];
C[(i)+(48)] += A[(i)+(0)] * values[15];
C[(i)+(48)] += A[(i)+(24)] * values[16];
C[(i)+(48)] += A[(i)+(40)] * values[17];
C[(i)+(56)] += A[(i)+(8)] * values[18];
C[(i)+(56)] += A[(i)+(24)] * values[19];
C[(i)+(56)] += A[(i)+(32)] * values[20];
C[(i)+(64)] += A[(i)+(16)] * values[21];
C[(i)+(64)] += A[(i)+(32)] * values[22];
C[(i)+(64)] += A[(i)+(40)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 192;
#endif

}

void dsparse_starMatrix_m10_n9_k9_ldA16_ldBna3_ldC16_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(32)
#pragma vector aligned
for (unsigned int i = 0; i < 10; i += 1)
{
C[(i)+(0)] += A[(i)+(96)] * values[0];
C[(i)+(0)] += A[(i)+(112)] * values[1];
C[(i)+(0)] += A[(i)+(128)] * values[2];
C[(i)+(16)] += A[(i)+(96)] * values[3];
C[(i)+(16)] += A[(i)+(112)] * values[4];
C[(i)+(16)] += A[(i)+(128)] * values[5];
C[(i)+(32)] += A[(i)+(96)] * values[6];
C[(i)+(32)] += A[(i)+(112)] * values[7];
C[(i)+(32)] += A[(i)+(128)] * values[8];
C[(i)+(48)] += A[(i)+(96)] * values[9];
C[(i)+(48)] += A[(i)+(112)] * values[10];
C[(i)+(64)] += A[(i)+(112)] * values[11];
C[(i)+(64)] += A[(i)+(128)] * values[12];
C[(i)+(80)] += A[(i)+(96)] * values[13];
C[(i)+(80)] += A[(i)+(128)] * values[14];
C[(i)+(96)] += A[(i)+(0)] * values[15];
C[(i)+(96)] += A[(i)+(48)] * values[16];
C[(i)+(96)] += A[(i)+(80)] * values[17];
C[(i)+(112)] += A[(i)+(16)] * values[18];
C[(i)+(112)] += A[(i)+(48)] * values[19];
C[(i)+(112)] += A[(i)+(64)] * values[20];
C[(i)+(128)] += A[(i)+(32)] * values[21];
C[(i)+(128)] += A[(i)+(64)] * values[22];
C[(i)+(128)] += A[(i)+(80)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 480;
#endif

}

void dsparse_starMatrix_m20_n9_k9_ldA24_ldBna4_ldC24_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(32)
#pragma vector aligned
for (unsigned int i = 0; i < 20; i += 1)
{
C[(i)+(0)] += A[(i)+(144)] * values[0];
C[(i)+(0)] += A[(i)+(168)] * values[1];
C[(i)+(0)] += A[(i)+(192)] * values[2];
C[(i)+(24)] += A[(i)+(144)] * values[3];
C[(i)+(24)] += A[(i)+(168)] * values[4];
C[(i)+(24)] += A[(i)+(192)] * values[5];
C[(i)+(48)] += A[(i)+(144)] * values[6];
C[(i)+(48)] += A[(i)+(168)] * values[7];
C[(i)+(48)] += A[(i)+(192)] * values[8];
C[(i)+(72)] += A[(i)+(144)] * values[9];
C[(i)+(72)] += A[(i)+(168)] * values[10];
C[(i)+(96)] += A[(i)+(168)] * values[11];
C[(i)+(96)] += A[(i)+(192)] * values[12];
C[(i)+(120)] += A[(i)+(144)] * values[13];
C[(i)+(120)] += A[(i)+(192)] * values[14];
C[(i)+(144)] += A[(i)+(0)] * values[15];
C[(i)+(144)] += A[(i)+(72)] * values[16];
C[(i)+(144)] += A[(i)+(120)] * values[17];
C[(i)+(168)] += A[(i)+(24)] * values[18];
C[(i)+(168)] += A[(i)+(72)] * values[19];
C[(i)+(168)] += A[(i)+(96)] * values[20];
C[(i)+(192)] += A[(i)+(48)] * values[21];
C[(i)+(192)] += A[(i)+(96)] * values[22];
C[(i)+(192)] += A[(i)+(120)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 960;
#endif

}

void dsparse_starMatrix_m35_n9_k9_ldA40_ldBna5_ldC40_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(32)
#pragma vector aligned
for (unsigned int i = 0; i < 35; i += 1)
{
C[(i)+(0)] += A[(i)+(240)] * values[0];
C[(i)+(0)] += A[(i)+(280)] * values[1];
C[(i)+(0)] += A[(i)+(320)] * values[2];
C[(i)+(40)] += A[(i)+(240)] * values[3];
C[(i)+(40)] += A[(i)+(280)] * values[4];
C[(i)+(40)] += A[(i)+(320)] * values[5];
C[(i)+(80)] += A[(i)+(240)] * values[6];
C[(i)+(80)] += A[(i)+(280)] * values[7];
C[(i)+(80)] += A[(i)+(320)] * values[8];
C[(i)+(120)] += A[(i)+(240)] * values[9];
C[(i)+(120)] += A[(i)+(280)] * values[10];
C[(i)+(160)] += A[(i)+(280)] * values[11];
C[(i)+(160)] += A[(i)+(320)] * values[12];
C[(i)+(200)] += A[(i)+(240)] * values[13];
C[(i)+(200)] += A[(i)+(320)] * values[14];
C[(i)+(240)] += A[(i)+(0)] * values[15];
C[(i)+(240)] += A[(i)+(120)] * values[16];
C[(i)+(240)] += A[(i)+(200)] * values[17];
C[(i)+(280)] += A[(i)+(40)] * values[18];
C[(i)+(280)] += A[(i)+(120)] * values[19];
C[(i)+(280)] += A[(i)+(160)] * values[20];
C[(i)+(320)] += A[(i)+(80)] * values[21];
C[(i)+(320)] += A[(i)+(160)] * values[22];
C[(i)+(320)] += A[(i)+(200)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 1680;
#endif

}

void dsparse_starMatrix_m56_n9_k9_ldA56_ldBna6_ldC56_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(32)
#pragma vector aligned
for (unsigned int i = 0; i < 56; i += 1)
{
C[(i)+(0)] += A[(i)+(336)] * values[0];
C[(i)+(0)] += A[(i)+(392)] * values[1];
C[(i)+(0)] += A[(i)+(448)] * values[2];
C[(i)+(56)] += A[(i)+(336)] * values[3];
C[(i)+(56)] += A[(i)+(392)] * values[4];
C[(i)+(56)] += A[(i)+(448)] * values[5];
C[(i)+(112)] += A[(i)+(336)] * values[6];
C[(i)+(112)] += A[(i)+(392)] * values[7];
C[(i)+(112)] += A[(i)+(448)] * values[8];
C[(i)+(168)] += A[(i)+(336)] * values[9];
C[(i)+(168)] += A[(i)+(392)] * values[10];
C[(i)+(224)] += A[(i)+(392)] * values[11];
C[(i)+(224)] += A[(i)+(448)] * values[12];
C[(i)+(280)] += A[(i)+(336)] * values[13];
C[(i)+(280)] += A[(i)+(448)] * values[14];
C[(i)+(336)] += A[(i)+(0)] * values[15];
C[(i)+(336)] += A[(i)+(168)] * values[16];
C[(i)+(336)] += A[(i)+(280)] * values[17];
C[(i)+(392)] += A[(i)+(56)] * values[18];
C[(i)+(392)] += A[(i)+(168)] * values[19];
C[(i)+(392)] += A[(i)+(224)] * values[20];
C[(i)+(448)] += A[(i)+(112)] * values[21];
C[(i)+(448)] += A[(i)+(224)] * values[22];
C[(i)+(448)] += A[(i)+(280)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 2688;
#endif

}

void dsparse_fP113DivM_m56_n9_k56_ldAna6_ldB56_ldC56_beta0_pfsigonly(const double* values, const double* B, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma nounroll_and_jam
for (unsigned int i = 0; i < 9; i++)
{
  #pragma simd
  for (unsigned int m = 0; m < 56; m++) {
    C[(i*56)+m] = 0.0;
  }
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b0 = _mm256_broadcast_sd(&B[(i*56)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b0 = _mm_loaddup_pd(&B[(i*56)+0]);
#endif
__m128d c0_0 = _mm_load_sd(&C[(i*56)+0]);
__m128d a0_0 = _mm_load_sd(&values[0]);
#if defined(__SSE3__) && defined(__AVX__)
c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, b0));
#endif
_mm_store_sd(&C[(i*56)+0], c0_0);
__m128d c0_1 = _mm_load_sd(&C[(i*56)+3]);
__m128d a0_1 = _mm_load_sd(&values[1]);
#if defined(__SSE3__) && defined(__AVX__)
c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, b0));
#endif
_mm_store_sd(&C[(i*56)+3], c0_1);
__m128d c0_2 = _mm_load_sd(&C[(i*56)+9]);
__m128d a0_2 = _mm_load_sd(&values[2]);
#if defined(__SSE3__) && defined(__AVX__)
c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, b0));
#endif
_mm_store_sd(&C[(i*56)+9], c0_2);
__m128d c0_3 = _mm_load_sd(&C[(i*56)+19]);
__m128d a0_3 = _mm_load_sd(&values[3]);
#if defined(__SSE3__) && defined(__AVX__)
c0_3 = _mm_add_sd(c0_3, _mm_mul_sd(a0_3, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_3 = _mm_add_sd(c0_3, _mm_mul_sd(a0_3, b0));
#endif
_mm_store_sd(&C[(i*56)+19], c0_3);
__m128d c0_4 = _mm_load_sd(&C[(i*56)+34]);
__m128d a0_4 = _mm_load_sd(&values[4]);
#if defined(__SSE3__) && defined(__AVX__)
c0_4 = _mm_add_sd(c0_4, _mm_mul_sd(a0_4, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_4 = _mm_add_sd(c0_4, _mm_mul_sd(a0_4, b0));
#endif
_mm_store_sd(&C[(i*56)+34], c0_4);
__m128d c0_5 = _mm_load_sd(&C[(i*56)+55]);
__m128d a0_5 = _mm_load_sd(&values[5]);
#if defined(__SSE3__) && defined(__AVX__)
c0_5 = _mm_add_sd(c0_5, _mm_mul_sd(a0_5, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_5 = _mm_add_sd(c0_5, _mm_mul_sd(a0_5, b0));
#endif
_mm_store_sd(&C[(i*56)+55], c0_5);
#else
C[(i*56)+0] += values[0] * B[(i*56)+0];
C[(i*56)+3] += values[1] * B[(i*56)+0];
C[(i*56)+9] += values[2] * B[(i*56)+0];
C[(i*56)+19] += values[3] * B[(i*56)+0];
C[(i*56)+34] += values[4] * B[(i*56)+0];
C[(i*56)+55] += values[5] * B[(i*56)+0];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b1 = _mm256_broadcast_sd(&B[(i*56)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b1 = _mm_loaddup_pd(&B[(i*56)+1]);
#endif
__m128d c1_0 = _mm_load_sd(&C[(i*56)+1]);
__m128d a1_0 = _mm_load_sd(&values[6]);
#if defined(__SSE3__) && defined(__AVX__)
c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, b1));
#endif
_mm_store_sd(&C[(i*56)+1], c1_0);
__m128d c1_1 = _mm_load_sd(&C[(i*56)+7]);
__m128d a1_1 = _mm_load_sd(&values[7]);
#if defined(__SSE3__) && defined(__AVX__)
c1_1 = _mm_add_sd(c1_1, _mm_mul_sd(a1_1, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_1 = _mm_add_sd(c1_1, _mm_mul_sd(a1_1, b1));
#endif
_mm_store_sd(&C[(i*56)+7], c1_1);
__m128d c1_2 = _mm_load_sd(&C[(i*56)+17]);
__m128d a1_2 = _mm_load_sd(&values[8]);
#if defined(__SSE3__) && defined(__AVX__)
c1_2 = _mm_add_sd(c1_2, _mm_mul_sd(a1_2, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_2 = _mm_add_sd(c1_2, _mm_mul_sd(a1_2, b1));
#endif
_mm_store_sd(&C[(i*56)+17], c1_2);
__m128d c1_3 = _mm_load_sd(&C[(i*56)+32]);
__m128d a1_3 = _mm_load_sd(&values[9]);
#if defined(__SSE3__) && defined(__AVX__)
c1_3 = _mm_add_sd(c1_3, _mm_mul_sd(a1_3, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_3 = _mm_add_sd(c1_3, _mm_mul_sd(a1_3, b1));
#endif
_mm_store_sd(&C[(i*56)+32], c1_3);
__m128d c1_4 = _mm_load_sd(&C[(i*56)+53]);
__m128d a1_4 = _mm_load_sd(&values[10]);
#if defined(__SSE3__) && defined(__AVX__)
c1_4 = _mm_add_sd(c1_4, _mm_mul_sd(a1_4, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_4 = _mm_add_sd(c1_4, _mm_mul_sd(a1_4, b1));
#endif
_mm_store_sd(&C[(i*56)+53], c1_4);
#else
C[(i*56)+1] += values[6] * B[(i*56)+1];
C[(i*56)+7] += values[7] * B[(i*56)+1];
C[(i*56)+17] += values[8] * B[(i*56)+1];
C[(i*56)+32] += values[9] * B[(i*56)+1];
C[(i*56)+53] += values[10] * B[(i*56)+1];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b2 = _mm256_broadcast_sd(&B[(i*56)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b2 = _mm_loaddup_pd(&B[(i*56)+2]);
#endif
__m128d c2_0 = _mm_load_sd(&C[(i*56)+2]);
__m128d a2_0 = _mm_load_sd(&values[11]);
#if defined(__SSE3__) && defined(__AVX__)
c2_0 = _mm_add_sd(c2_0, _mm_mul_sd(a2_0, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_0 = _mm_add_sd(c2_0, _mm_mul_sd(a2_0, b2));
#endif
_mm_store_sd(&C[(i*56)+2], c2_0);
__m128d c2_1 = _mm_load_sd(&C[(i*56)+8]);
__m128d a2_1 = _mm_load_sd(&values[12]);
#if defined(__SSE3__) && defined(__AVX__)
c2_1 = _mm_add_sd(c2_1, _mm_mul_sd(a2_1, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_1 = _mm_add_sd(c2_1, _mm_mul_sd(a2_1, b2));
#endif
_mm_store_sd(&C[(i*56)+8], c2_1);
__m128d c2_2 = _mm_load_sd(&C[(i*56)+18]);
__m128d a2_2 = _mm_load_sd(&values[13]);
#if defined(__SSE3__) && defined(__AVX__)
c2_2 = _mm_add_sd(c2_2, _mm_mul_sd(a2_2, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_2 = _mm_add_sd(c2_2, _mm_mul_sd(a2_2, b2));
#endif
_mm_store_sd(&C[(i*56)+18], c2_2);
__m128d c2_3 = _mm_load_sd(&C[(i*56)+33]);
__m128d a2_3 = _mm_load_sd(&values[14]);
#if defined(__SSE3__) && defined(__AVX__)
c2_3 = _mm_add_sd(c2_3, _mm_mul_sd(a2_3, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_3 = _mm_add_sd(c2_3, _mm_mul_sd(a2_3, b2));
#endif
_mm_store_sd(&C[(i*56)+33], c2_3);
__m128d c2_4 = _mm_load_sd(&C[(i*56)+54]);
__m128d a2_4 = _mm_load_sd(&values[15]);
#if defined(__SSE3__) && defined(__AVX__)
c2_4 = _mm_add_sd(c2_4, _mm_mul_sd(a2_4, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_4 = _mm_add_sd(c2_4, _mm_mul_sd(a2_4, b2));
#endif
_mm_store_sd(&C[(i*56)+54], c2_4);
#else
C[(i*56)+2] += values[11] * B[(i*56)+2];
C[(i*56)+8] += values[12] * B[(i*56)+2];
C[(i*56)+18] += values[13] * B[(i*56)+2];
C[(i*56)+33] += values[14] * B[(i*56)+2];
C[(i*56)+54] += values[15] * B[(i*56)+2];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b3 = _mm256_broadcast_sd(&B[(i*56)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b3 = _mm_loaddup_pd(&B[(i*56)+3]);
#endif
__m128d c3_0 = _mm_load_sd(&C[(i*56)+0]);
__m128d a3_0 = _mm_load_sd(&values[16]);
#if defined(__SSE3__) && defined(__AVX__)
c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, b3));
#endif
_mm_store_sd(&C[(i*56)+0], c3_0);
__m128d c3_1 = _mm_load_sd(&C[(i*56)+3]);
__m128d a3_1 = _mm_load_sd(&values[17]);
#if defined(__SSE3__) && defined(__AVX__)
c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, b3));
#endif
_mm_store_sd(&C[(i*56)+3], c3_1);
__m128d c3_2 = _mm_load_sd(&C[(i*56)+9]);
__m128d a3_2 = _mm_load_sd(&values[18]);
#if defined(__SSE3__) && defined(__AVX__)
c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, b3));
#endif
_mm_store_sd(&C[(i*56)+9], c3_2);
__m128d c3_3 = _mm_load_sd(&C[(i*56)+19]);
__m128d a3_3 = _mm_load_sd(&values[19]);
#if defined(__SSE3__) && defined(__AVX__)
c3_3 = _mm_add_sd(c3_3, _mm_mul_sd(a3_3, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_3 = _mm_add_sd(c3_3, _mm_mul_sd(a3_3, b3));
#endif
_mm_store_sd(&C[(i*56)+19], c3_3);
__m128d c3_4 = _mm_load_sd(&C[(i*56)+34]);
__m128d a3_4 = _mm_load_sd(&values[20]);
#if defined(__SSE3__) && defined(__AVX__)
c3_4 = _mm_add_sd(c3_4, _mm_mul_sd(a3_4, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_4 = _mm_add_sd(c3_4, _mm_mul_sd(a3_4, b3));
#endif
_mm_store_sd(&C[(i*56)+34], c3_4);
__m128d c3_5 = _mm_load_sd(&C[(i*56)+55]);
__m128d a3_5 = _mm_load_sd(&values[21]);
#if defined(__SSE3__) && defined(__AVX__)
c3_5 = _mm_add_sd(c3_5, _mm_mul_sd(a3_5, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_5 = _mm_add_sd(c3_5, _mm_mul_sd(a3_5, b3));
#endif
_mm_store_sd(&C[(i*56)+55], c3_5);
#else
C[(i*56)+0] += values[16] * B[(i*56)+3];
C[(i*56)+3] += values[17] * B[(i*56)+3];
C[(i*56)+9] += values[18] * B[(i*56)+3];
C[(i*56)+19] += values[19] * B[(i*56)+3];
C[(i*56)+34] += values[20] * B[(i*56)+3];
C[(i*56)+55] += values[21] * B[(i*56)+3];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b4 = _mm256_broadcast_sd(&B[(i*56)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b4 = _mm_loaddup_pd(&B[(i*56)+4]);
#endif
__m128d c4_0 = _mm_load_sd(&C[(i*56)+4]);
__m128d a4_0 = _mm_load_sd(&values[22]);
#if defined(__SSE3__) && defined(__AVX__)
c4_0 = _mm_add_sd(c4_0, _mm_mul_sd(a4_0, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_0 = _mm_add_sd(c4_0, _mm_mul_sd(a4_0, b4));
#endif
_mm_store_sd(&C[(i*56)+4], c4_0);
__m128d c4_1 = _mm_load_sd(&C[(i*56)+14]);
__m128d a4_1 = _mm_load_sd(&values[23]);
#if defined(__SSE3__) && defined(__AVX__)
c4_1 = _mm_add_sd(c4_1, _mm_mul_sd(a4_1, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_1 = _mm_add_sd(c4_1, _mm_mul_sd(a4_1, b4));
#endif
_mm_store_sd(&C[(i*56)+14], c4_1);
__m128d c4_2 = _mm_load_sd(&C[(i*56)+29]);
__m128d a4_2 = _mm_load_sd(&values[24]);
#if defined(__SSE3__) && defined(__AVX__)
c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, b4));
#endif
_mm_store_sd(&C[(i*56)+29], c4_2);
__m128d c4_3 = _mm_load_sd(&C[(i*56)+50]);
__m128d a4_3 = _mm_load_sd(&values[25]);
#if defined(__SSE3__) && defined(__AVX__)
c4_3 = _mm_add_sd(c4_3, _mm_mul_sd(a4_3, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_3 = _mm_add_sd(c4_3, _mm_mul_sd(a4_3, b4));
#endif
_mm_store_sd(&C[(i*56)+50], c4_3);
#else
C[(i*56)+4] += values[22] * B[(i*56)+4];
C[(i*56)+14] += values[23] * B[(i*56)+4];
C[(i*56)+29] += values[24] * B[(i*56)+4];
C[(i*56)+50] += values[25] * B[(i*56)+4];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b5 = _mm256_broadcast_sd(&B[(i*56)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b5 = _mm_loaddup_pd(&B[(i*56)+5]);
#endif
__m128d c5_0 = _mm_load_sd(&C[(i*56)+5]);
__m128d a5_0 = _mm_load_sd(&values[26]);
#if defined(__SSE3__) && defined(__AVX__)
c5_0 = _mm_add_sd(c5_0, _mm_mul_sd(a5_0, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_0 = _mm_add_sd(c5_0, _mm_mul_sd(a5_0, b5));
#endif
_mm_store_sd(&C[(i*56)+5], c5_0);
__m128d c5_1 = _mm_load_sd(&C[(i*56)+15]);
__m128d a5_1 = _mm_load_sd(&values[27]);
#if defined(__SSE3__) && defined(__AVX__)
c5_1 = _mm_add_sd(c5_1, _mm_mul_sd(a5_1, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_1 = _mm_add_sd(c5_1, _mm_mul_sd(a5_1, b5));
#endif
_mm_store_sd(&C[(i*56)+15], c5_1);
__m128d c5_2 = _mm_load_sd(&C[(i*56)+30]);
__m128d a5_2 = _mm_load_sd(&values[28]);
#if defined(__SSE3__) && defined(__AVX__)
c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, b5));
#endif
_mm_store_sd(&C[(i*56)+30], c5_2);
__m128d c5_3 = _mm_load_sd(&C[(i*56)+51]);
__m128d a5_3 = _mm_load_sd(&values[29]);
#if defined(__SSE3__) && defined(__AVX__)
c5_3 = _mm_add_sd(c5_3, _mm_mul_sd(a5_3, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_3 = _mm_add_sd(c5_3, _mm_mul_sd(a5_3, b5));
#endif
_mm_store_sd(&C[(i*56)+51], c5_3);
#else
C[(i*56)+5] += values[26] * B[(i*56)+5];
C[(i*56)+15] += values[27] * B[(i*56)+5];
C[(i*56)+30] += values[28] * B[(i*56)+5];
C[(i*56)+51] += values[29] * B[(i*56)+5];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b6 = _mm256_broadcast_sd(&B[(i*56)+6]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b6 = _mm_loaddup_pd(&B[(i*56)+6]);
#endif
__m128d c6_0 = _mm_load_sd(&C[(i*56)+6]);
__m128d a6_0 = _mm_load_sd(&values[30]);
#if defined(__SSE3__) && defined(__AVX__)
c6_0 = _mm_add_sd(c6_0, _mm_mul_sd(a6_0, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_0 = _mm_add_sd(c6_0, _mm_mul_sd(a6_0, b6));
#endif
_mm_store_sd(&C[(i*56)+6], c6_0);
__m128d c6_1 = _mm_load_sd(&C[(i*56)+16]);
__m128d a6_1 = _mm_load_sd(&values[31]);
#if defined(__SSE3__) && defined(__AVX__)
c6_1 = _mm_add_sd(c6_1, _mm_mul_sd(a6_1, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_1 = _mm_add_sd(c6_1, _mm_mul_sd(a6_1, b6));
#endif
_mm_store_sd(&C[(i*56)+16], c6_1);
__m128d c6_2 = _mm_load_sd(&C[(i*56)+31]);
__m128d a6_2 = _mm_load_sd(&values[32]);
#if defined(__SSE3__) && defined(__AVX__)
c6_2 = _mm_add_sd(c6_2, _mm_mul_sd(a6_2, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_2 = _mm_add_sd(c6_2, _mm_mul_sd(a6_2, b6));
#endif
_mm_store_sd(&C[(i*56)+31], c6_2);
__m128d c6_3 = _mm_load_sd(&C[(i*56)+52]);
__m128d a6_3 = _mm_load_sd(&values[33]);
#if defined(__SSE3__) && defined(__AVX__)
c6_3 = _mm_add_sd(c6_3, _mm_mul_sd(a6_3, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_3 = _mm_add_sd(c6_3, _mm_mul_sd(a6_3, b6));
#endif
_mm_store_sd(&C[(i*56)+52], c6_3);
#else
C[(i*56)+6] += values[30] * B[(i*56)+6];
C[(i*56)+16] += values[31] * B[(i*56)+6];
C[(i*56)+31] += values[32] * B[(i*56)+6];
C[(i*56)+52] += values[33] * B[(i*56)+6];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b7 = _mm256_broadcast_sd(&B[(i*56)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b7 = _mm_loaddup_pd(&B[(i*56)+7]);
#endif
__m128d c7_0 = _mm_load_sd(&C[(i*56)+1]);
__m128d a7_0 = _mm_load_sd(&values[34]);
#if defined(__SSE3__) && defined(__AVX__)
c7_0 = _mm_add_sd(c7_0, _mm_mul_sd(a7_0, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_0 = _mm_add_sd(c7_0, _mm_mul_sd(a7_0, b7));
#endif
_mm_store_sd(&C[(i*56)+1], c7_0);
__m128d c7_1 = _mm_load_sd(&C[(i*56)+7]);
__m128d a7_1 = _mm_load_sd(&values[35]);
#if defined(__SSE3__) && defined(__AVX__)
c7_1 = _mm_add_sd(c7_1, _mm_mul_sd(a7_1, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_1 = _mm_add_sd(c7_1, _mm_mul_sd(a7_1, b7));
#endif
_mm_store_sd(&C[(i*56)+7], c7_1);
__m128d c7_2 = _mm_load_sd(&C[(i*56)+17]);
__m128d a7_2 = _mm_load_sd(&values[36]);
#if defined(__SSE3__) && defined(__AVX__)
c7_2 = _mm_add_sd(c7_2, _mm_mul_sd(a7_2, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_2 = _mm_add_sd(c7_2, _mm_mul_sd(a7_2, b7));
#endif
_mm_store_sd(&C[(i*56)+17], c7_2);
__m128d c7_3 = _mm_load_sd(&C[(i*56)+32]);
__m128d a7_3 = _mm_load_sd(&values[37]);
#if defined(__SSE3__) && defined(__AVX__)
c7_3 = _mm_add_sd(c7_3, _mm_mul_sd(a7_3, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_3 = _mm_add_sd(c7_3, _mm_mul_sd(a7_3, b7));
#endif
_mm_store_sd(&C[(i*56)+32], c7_3);
__m128d c7_4 = _mm_load_sd(&C[(i*56)+53]);
__m128d a7_4 = _mm_load_sd(&values[38]);
#if defined(__SSE3__) && defined(__AVX__)
c7_4 = _mm_add_sd(c7_4, _mm_mul_sd(a7_4, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_4 = _mm_add_sd(c7_4, _mm_mul_sd(a7_4, b7));
#endif
_mm_store_sd(&C[(i*56)+53], c7_4);
#else
C[(i*56)+1] += values[34] * B[(i*56)+7];
C[(i*56)+7] += values[35] * B[(i*56)+7];
C[(i*56)+17] += values[36] * B[(i*56)+7];
C[(i*56)+32] += values[37] * B[(i*56)+7];
C[(i*56)+53] += values[38] * B[(i*56)+7];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b8 = _mm256_broadcast_sd(&B[(i*56)+8]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b8 = _mm_loaddup_pd(&B[(i*56)+8]);
#endif
__m128d c8_0 = _mm_load_sd(&C[(i*56)+2]);
__m128d a8_0 = _mm_load_sd(&values[39]);
#if defined(__SSE3__) && defined(__AVX__)
c8_0 = _mm_add_sd(c8_0, _mm_mul_sd(a8_0, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_0 = _mm_add_sd(c8_0, _mm_mul_sd(a8_0, b8));
#endif
_mm_store_sd(&C[(i*56)+2], c8_0);
__m128d c8_1 = _mm_load_sd(&C[(i*56)+8]);
__m128d a8_1 = _mm_load_sd(&values[40]);
#if defined(__SSE3__) && defined(__AVX__)
c8_1 = _mm_add_sd(c8_1, _mm_mul_sd(a8_1, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_1 = _mm_add_sd(c8_1, _mm_mul_sd(a8_1, b8));
#endif
_mm_store_sd(&C[(i*56)+8], c8_1);
__m128d c8_2 = _mm_load_sd(&C[(i*56)+18]);
__m128d a8_2 = _mm_load_sd(&values[41]);
#if defined(__SSE3__) && defined(__AVX__)
c8_2 = _mm_add_sd(c8_2, _mm_mul_sd(a8_2, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_2 = _mm_add_sd(c8_2, _mm_mul_sd(a8_2, b8));
#endif
_mm_store_sd(&C[(i*56)+18], c8_2);
__m128d c8_3 = _mm_load_sd(&C[(i*56)+33]);
__m128d a8_3 = _mm_load_sd(&values[42]);
#if defined(__SSE3__) && defined(__AVX__)
c8_3 = _mm_add_sd(c8_3, _mm_mul_sd(a8_3, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_3 = _mm_add_sd(c8_3, _mm_mul_sd(a8_3, b8));
#endif
_mm_store_sd(&C[(i*56)+33], c8_3);
__m128d c8_4 = _mm_load_sd(&C[(i*56)+54]);
__m128d a8_4 = _mm_load_sd(&values[43]);
#if defined(__SSE3__) && defined(__AVX__)
c8_4 = _mm_add_sd(c8_4, _mm_mul_sd(a8_4, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_4 = _mm_add_sd(c8_4, _mm_mul_sd(a8_4, b8));
#endif
_mm_store_sd(&C[(i*56)+54], c8_4);
#else
C[(i*56)+2] += values[39] * B[(i*56)+8];
C[(i*56)+8] += values[40] * B[(i*56)+8];
C[(i*56)+18] += values[41] * B[(i*56)+8];
C[(i*56)+33] += values[42] * B[(i*56)+8];
C[(i*56)+54] += values[43] * B[(i*56)+8];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b9 = _mm256_broadcast_sd(&B[(i*56)+9]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b9 = _mm_loaddup_pd(&B[(i*56)+9]);
#endif
__m128d c9_0 = _mm_load_sd(&C[(i*56)+0]);
__m128d a9_0 = _mm_load_sd(&values[44]);
#if defined(__SSE3__) && defined(__AVX__)
c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, b9));
#endif
_mm_store_sd(&C[(i*56)+0], c9_0);
__m128d c9_1 = _mm_load_sd(&C[(i*56)+3]);
__m128d a9_1 = _mm_load_sd(&values[45]);
#if defined(__SSE3__) && defined(__AVX__)
c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, b9));
#endif
_mm_store_sd(&C[(i*56)+3], c9_1);
__m128d c9_2 = _mm_load_sd(&C[(i*56)+9]);
__m128d a9_2 = _mm_load_sd(&values[46]);
#if defined(__SSE3__) && defined(__AVX__)
c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, b9));
#endif
_mm_store_sd(&C[(i*56)+9], c9_2);
__m128d c9_3 = _mm_load_sd(&C[(i*56)+19]);
__m128d a9_3 = _mm_load_sd(&values[47]);
#if defined(__SSE3__) && defined(__AVX__)
c9_3 = _mm_add_sd(c9_3, _mm_mul_sd(a9_3, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_3 = _mm_add_sd(c9_3, _mm_mul_sd(a9_3, b9));
#endif
_mm_store_sd(&C[(i*56)+19], c9_3);
__m128d c9_4 = _mm_load_sd(&C[(i*56)+34]);
__m128d a9_4 = _mm_load_sd(&values[48]);
#if defined(__SSE3__) && defined(__AVX__)
c9_4 = _mm_add_sd(c9_4, _mm_mul_sd(a9_4, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_4 = _mm_add_sd(c9_4, _mm_mul_sd(a9_4, b9));
#endif
_mm_store_sd(&C[(i*56)+34], c9_4);
__m128d c9_5 = _mm_load_sd(&C[(i*56)+55]);
__m128d a9_5 = _mm_load_sd(&values[49]);
#if defined(__SSE3__) && defined(__AVX__)
c9_5 = _mm_add_sd(c9_5, _mm_mul_sd(a9_5, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_5 = _mm_add_sd(c9_5, _mm_mul_sd(a9_5, b9));
#endif
_mm_store_sd(&C[(i*56)+55], c9_5);
#else
C[(i*56)+0] += values[44] * B[(i*56)+9];
C[(i*56)+3] += values[45] * B[(i*56)+9];
C[(i*56)+9] += values[46] * B[(i*56)+9];
C[(i*56)+19] += values[47] * B[(i*56)+9];
C[(i*56)+34] += values[48] * B[(i*56)+9];
C[(i*56)+55] += values[49] * B[(i*56)+9];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b10 = _mm256_broadcast_sd(&B[(i*56)+10]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b10 = _mm_loaddup_pd(&B[(i*56)+10]);
#endif
__m128d c10_0 = _mm_load_sd(&C[(i*56)+10]);
__m128d a10_0 = _mm_load_sd(&values[50]);
#if defined(__SSE3__) && defined(__AVX__)
c10_0 = _mm_add_sd(c10_0, _mm_mul_sd(a10_0, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c10_0 = _mm_add_sd(c10_0, _mm_mul_sd(a10_0, b10));
#endif
_mm_store_sd(&C[(i*56)+10], c10_0);
__m128d c10_1 = _mm_load_sd(&C[(i*56)+25]);
__m128d a10_1 = _mm_load_sd(&values[51]);
#if defined(__SSE3__) && defined(__AVX__)
c10_1 = _mm_add_sd(c10_1, _mm_mul_sd(a10_1, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c10_1 = _mm_add_sd(c10_1, _mm_mul_sd(a10_1, b10));
#endif
_mm_store_sd(&C[(i*56)+25], c10_1);
__m128d c10_2 = _mm_load_sd(&C[(i*56)+46]);
__m128d a10_2 = _mm_load_sd(&values[52]);
#if defined(__SSE3__) && defined(__AVX__)
c10_2 = _mm_add_sd(c10_2, _mm_mul_sd(a10_2, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c10_2 = _mm_add_sd(c10_2, _mm_mul_sd(a10_2, b10));
#endif
_mm_store_sd(&C[(i*56)+46], c10_2);
#else
C[(i*56)+10] += values[50] * B[(i*56)+10];
C[(i*56)+25] += values[51] * B[(i*56)+10];
C[(i*56)+46] += values[52] * B[(i*56)+10];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b11 = _mm256_broadcast_sd(&B[(i*56)+11]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b11 = _mm_loaddup_pd(&B[(i*56)+11]);
#endif
__m128d c11_0 = _mm_load_sd(&C[(i*56)+11]);
__m128d a11_0 = _mm_load_sd(&values[53]);
#if defined(__SSE3__) && defined(__AVX__)
c11_0 = _mm_add_sd(c11_0, _mm_mul_sd(a11_0, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c11_0 = _mm_add_sd(c11_0, _mm_mul_sd(a11_0, b11));
#endif
_mm_store_sd(&C[(i*56)+11], c11_0);
__m128d c11_1 = _mm_load_sd(&C[(i*56)+26]);
__m128d a11_1 = _mm_load_sd(&values[54]);
#if defined(__SSE3__) && defined(__AVX__)
c11_1 = _mm_add_sd(c11_1, _mm_mul_sd(a11_1, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c11_1 = _mm_add_sd(c11_1, _mm_mul_sd(a11_1, b11));
#endif
_mm_store_sd(&C[(i*56)+26], c11_1);
__m128d c11_2 = _mm_load_sd(&C[(i*56)+47]);
__m128d a11_2 = _mm_load_sd(&values[55]);
#if defined(__SSE3__) && defined(__AVX__)
c11_2 = _mm_add_sd(c11_2, _mm_mul_sd(a11_2, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c11_2 = _mm_add_sd(c11_2, _mm_mul_sd(a11_2, b11));
#endif
_mm_store_sd(&C[(i*56)+47], c11_2);
#else
C[(i*56)+11] += values[53] * B[(i*56)+11];
C[(i*56)+26] += values[54] * B[(i*56)+11];
C[(i*56)+47] += values[55] * B[(i*56)+11];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b12 = _mm256_broadcast_sd(&B[(i*56)+12]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b12 = _mm_loaddup_pd(&B[(i*56)+12]);
#endif
__m128d c12_0 = _mm_load_sd(&C[(i*56)+12]);
__m128d a12_0 = _mm_load_sd(&values[56]);
#if defined(__SSE3__) && defined(__AVX__)
c12_0 = _mm_add_sd(c12_0, _mm_mul_sd(a12_0, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c12_0 = _mm_add_sd(c12_0, _mm_mul_sd(a12_0, b12));
#endif
_mm_store_sd(&C[(i*56)+12], c12_0);
__m128d c12_1 = _mm_load_sd(&C[(i*56)+27]);
__m128d a12_1 = _mm_load_sd(&values[57]);
#if defined(__SSE3__) && defined(__AVX__)
c12_1 = _mm_add_sd(c12_1, _mm_mul_sd(a12_1, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c12_1 = _mm_add_sd(c12_1, _mm_mul_sd(a12_1, b12));
#endif
_mm_store_sd(&C[(i*56)+27], c12_1);
__m128d c12_2 = _mm_load_sd(&C[(i*56)+48]);
__m128d a12_2 = _mm_load_sd(&values[58]);
#if defined(__SSE3__) && defined(__AVX__)
c12_2 = _mm_add_sd(c12_2, _mm_mul_sd(a12_2, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c12_2 = _mm_add_sd(c12_2, _mm_mul_sd(a12_2, b12));
#endif
_mm_store_sd(&C[(i*56)+48], c12_2);
#else
C[(i*56)+12] += values[56] * B[(i*56)+12];
C[(i*56)+27] += values[57] * B[(i*56)+12];
C[(i*56)+48] += values[58] * B[(i*56)+12];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b13 = _mm256_broadcast_sd(&B[(i*56)+13]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b13 = _mm_loaddup_pd(&B[(i*56)+13]);
#endif
__m128d c13_0 = _mm_load_sd(&C[(i*56)+13]);
__m128d a13_0 = _mm_load_sd(&values[59]);
#if defined(__SSE3__) && defined(__AVX__)
c13_0 = _mm_add_sd(c13_0, _mm_mul_sd(a13_0, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c13_0 = _mm_add_sd(c13_0, _mm_mul_sd(a13_0, b13));
#endif
_mm_store_sd(&C[(i*56)+13], c13_0);
__m128d c13_1 = _mm_load_sd(&C[(i*56)+28]);
__m128d a13_1 = _mm_load_sd(&values[60]);
#if defined(__SSE3__) && defined(__AVX__)
c13_1 = _mm_add_sd(c13_1, _mm_mul_sd(a13_1, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c13_1 = _mm_add_sd(c13_1, _mm_mul_sd(a13_1, b13));
#endif
_mm_store_sd(&C[(i*56)+28], c13_1);
__m128d c13_2 = _mm_load_sd(&C[(i*56)+49]);
__m128d a13_2 = _mm_load_sd(&values[61]);
#if defined(__SSE3__) && defined(__AVX__)
c13_2 = _mm_add_sd(c13_2, _mm_mul_sd(a13_2, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c13_2 = _mm_add_sd(c13_2, _mm_mul_sd(a13_2, b13));
#endif
_mm_store_sd(&C[(i*56)+49], c13_2);
#else
C[(i*56)+13] += values[59] * B[(i*56)+13];
C[(i*56)+28] += values[60] * B[(i*56)+13];
C[(i*56)+49] += values[61] * B[(i*56)+13];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b14 = _mm256_broadcast_sd(&B[(i*56)+14]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b14 = _mm_loaddup_pd(&B[(i*56)+14]);
#endif
__m128d c14_0 = _mm_load_sd(&C[(i*56)+4]);
__m128d a14_0 = _mm_load_sd(&values[62]);
#if defined(__SSE3__) && defined(__AVX__)
c14_0 = _mm_add_sd(c14_0, _mm_mul_sd(a14_0, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_0 = _mm_add_sd(c14_0, _mm_mul_sd(a14_0, b14));
#endif
_mm_store_sd(&C[(i*56)+4], c14_0);
__m128d c14_1 = _mm_load_sd(&C[(i*56)+14]);
__m128d a14_1 = _mm_load_sd(&values[63]);
#if defined(__SSE3__) && defined(__AVX__)
c14_1 = _mm_add_sd(c14_1, _mm_mul_sd(a14_1, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_1 = _mm_add_sd(c14_1, _mm_mul_sd(a14_1, b14));
#endif
_mm_store_sd(&C[(i*56)+14], c14_1);
__m128d c14_2 = _mm_load_sd(&C[(i*56)+29]);
__m128d a14_2 = _mm_load_sd(&values[64]);
#if defined(__SSE3__) && defined(__AVX__)
c14_2 = _mm_add_sd(c14_2, _mm_mul_sd(a14_2, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_2 = _mm_add_sd(c14_2, _mm_mul_sd(a14_2, b14));
#endif
_mm_store_sd(&C[(i*56)+29], c14_2);
__m128d c14_3 = _mm_load_sd(&C[(i*56)+50]);
__m128d a14_3 = _mm_load_sd(&values[65]);
#if defined(__SSE3__) && defined(__AVX__)
c14_3 = _mm_add_sd(c14_3, _mm_mul_sd(a14_3, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_3 = _mm_add_sd(c14_3, _mm_mul_sd(a14_3, b14));
#endif
_mm_store_sd(&C[(i*56)+50], c14_3);
#else
C[(i*56)+4] += values[62] * B[(i*56)+14];
C[(i*56)+14] += values[63] * B[(i*56)+14];
C[(i*56)+29] += values[64] * B[(i*56)+14];
C[(i*56)+50] += values[65] * B[(i*56)+14];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b15 = _mm256_broadcast_sd(&B[(i*56)+15]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b15 = _mm_loaddup_pd(&B[(i*56)+15]);
#endif
__m128d c15_0 = _mm_load_sd(&C[(i*56)+5]);
__m128d a15_0 = _mm_load_sd(&values[66]);
#if defined(__SSE3__) && defined(__AVX__)
c15_0 = _mm_add_sd(c15_0, _mm_mul_sd(a15_0, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_0 = _mm_add_sd(c15_0, _mm_mul_sd(a15_0, b15));
#endif
_mm_store_sd(&C[(i*56)+5], c15_0);
__m128d c15_1 = _mm_load_sd(&C[(i*56)+15]);
__m128d a15_1 = _mm_load_sd(&values[67]);
#if defined(__SSE3__) && defined(__AVX__)
c15_1 = _mm_add_sd(c15_1, _mm_mul_sd(a15_1, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_1 = _mm_add_sd(c15_1, _mm_mul_sd(a15_1, b15));
#endif
_mm_store_sd(&C[(i*56)+15], c15_1);
__m128d c15_2 = _mm_load_sd(&C[(i*56)+30]);
__m128d a15_2 = _mm_load_sd(&values[68]);
#if defined(__SSE3__) && defined(__AVX__)
c15_2 = _mm_add_sd(c15_2, _mm_mul_sd(a15_2, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_2 = _mm_add_sd(c15_2, _mm_mul_sd(a15_2, b15));
#endif
_mm_store_sd(&C[(i*56)+30], c15_2);
__m128d c15_3 = _mm_load_sd(&C[(i*56)+51]);
__m128d a15_3 = _mm_load_sd(&values[69]);
#if defined(__SSE3__) && defined(__AVX__)
c15_3 = _mm_add_sd(c15_3, _mm_mul_sd(a15_3, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_3 = _mm_add_sd(c15_3, _mm_mul_sd(a15_3, b15));
#endif
_mm_store_sd(&C[(i*56)+51], c15_3);
#else
C[(i*56)+5] += values[66] * B[(i*56)+15];
C[(i*56)+15] += values[67] * B[(i*56)+15];
C[(i*56)+30] += values[68] * B[(i*56)+15];
C[(i*56)+51] += values[69] * B[(i*56)+15];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b16 = _mm256_broadcast_sd(&B[(i*56)+16]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b16 = _mm_loaddup_pd(&B[(i*56)+16]);
#endif
__m128d c16_0 = _mm_load_sd(&C[(i*56)+6]);
__m128d a16_0 = _mm_load_sd(&values[70]);
#if defined(__SSE3__) && defined(__AVX__)
c16_0 = _mm_add_sd(c16_0, _mm_mul_sd(a16_0, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_0 = _mm_add_sd(c16_0, _mm_mul_sd(a16_0, b16));
#endif
_mm_store_sd(&C[(i*56)+6], c16_0);
__m128d c16_1 = _mm_load_sd(&C[(i*56)+16]);
__m128d a16_1 = _mm_load_sd(&values[71]);
#if defined(__SSE3__) && defined(__AVX__)
c16_1 = _mm_add_sd(c16_1, _mm_mul_sd(a16_1, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_1 = _mm_add_sd(c16_1, _mm_mul_sd(a16_1, b16));
#endif
_mm_store_sd(&C[(i*56)+16], c16_1);
__m128d c16_2 = _mm_load_sd(&C[(i*56)+31]);
__m128d a16_2 = _mm_load_sd(&values[72]);
#if defined(__SSE3__) && defined(__AVX__)
c16_2 = _mm_add_sd(c16_2, _mm_mul_sd(a16_2, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_2 = _mm_add_sd(c16_2, _mm_mul_sd(a16_2, b16));
#endif
_mm_store_sd(&C[(i*56)+31], c16_2);
__m128d c16_3 = _mm_load_sd(&C[(i*56)+52]);
__m128d a16_3 = _mm_load_sd(&values[73]);
#if defined(__SSE3__) && defined(__AVX__)
c16_3 = _mm_add_sd(c16_3, _mm_mul_sd(a16_3, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_3 = _mm_add_sd(c16_3, _mm_mul_sd(a16_3, b16));
#endif
_mm_store_sd(&C[(i*56)+52], c16_3);
#else
C[(i*56)+6] += values[70] * B[(i*56)+16];
C[(i*56)+16] += values[71] * B[(i*56)+16];
C[(i*56)+31] += values[72] * B[(i*56)+16];
C[(i*56)+52] += values[73] * B[(i*56)+16];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b17 = _mm256_broadcast_sd(&B[(i*56)+17]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b17 = _mm_loaddup_pd(&B[(i*56)+17]);
#endif
__m128d c17_0 = _mm_load_sd(&C[(i*56)+1]);
__m128d a17_0 = _mm_load_sd(&values[74]);
#if defined(__SSE3__) && defined(__AVX__)
c17_0 = _mm_add_sd(c17_0, _mm_mul_sd(a17_0, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_0 = _mm_add_sd(c17_0, _mm_mul_sd(a17_0, b17));
#endif
_mm_store_sd(&C[(i*56)+1], c17_0);
__m128d c17_1 = _mm_load_sd(&C[(i*56)+7]);
__m128d a17_1 = _mm_load_sd(&values[75]);
#if defined(__SSE3__) && defined(__AVX__)
c17_1 = _mm_add_sd(c17_1, _mm_mul_sd(a17_1, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_1 = _mm_add_sd(c17_1, _mm_mul_sd(a17_1, b17));
#endif
_mm_store_sd(&C[(i*56)+7], c17_1);
__m128d c17_2 = _mm_load_sd(&C[(i*56)+17]);
__m128d a17_2 = _mm_load_sd(&values[76]);
#if defined(__SSE3__) && defined(__AVX__)
c17_2 = _mm_add_sd(c17_2, _mm_mul_sd(a17_2, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_2 = _mm_add_sd(c17_2, _mm_mul_sd(a17_2, b17));
#endif
_mm_store_sd(&C[(i*56)+17], c17_2);
__m128d c17_3 = _mm_load_sd(&C[(i*56)+32]);
__m128d a17_3 = _mm_load_sd(&values[77]);
#if defined(__SSE3__) && defined(__AVX__)
c17_3 = _mm_add_sd(c17_3, _mm_mul_sd(a17_3, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_3 = _mm_add_sd(c17_3, _mm_mul_sd(a17_3, b17));
#endif
_mm_store_sd(&C[(i*56)+32], c17_3);
__m128d c17_4 = _mm_load_sd(&C[(i*56)+53]);
__m128d a17_4 = _mm_load_sd(&values[78]);
#if defined(__SSE3__) && defined(__AVX__)
c17_4 = _mm_add_sd(c17_4, _mm_mul_sd(a17_4, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_4 = _mm_add_sd(c17_4, _mm_mul_sd(a17_4, b17));
#endif
_mm_store_sd(&C[(i*56)+53], c17_4);
#else
C[(i*56)+1] += values[74] * B[(i*56)+17];
C[(i*56)+7] += values[75] * B[(i*56)+17];
C[(i*56)+17] += values[76] * B[(i*56)+17];
C[(i*56)+32] += values[77] * B[(i*56)+17];
C[(i*56)+53] += values[78] * B[(i*56)+17];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b18 = _mm256_broadcast_sd(&B[(i*56)+18]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b18 = _mm_loaddup_pd(&B[(i*56)+18]);
#endif
__m128d c18_0 = _mm_load_sd(&C[(i*56)+2]);
__m128d a18_0 = _mm_load_sd(&values[79]);
#if defined(__SSE3__) && defined(__AVX__)
c18_0 = _mm_add_sd(c18_0, _mm_mul_sd(a18_0, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_0 = _mm_add_sd(c18_0, _mm_mul_sd(a18_0, b18));
#endif
_mm_store_sd(&C[(i*56)+2], c18_0);
__m128d c18_1 = _mm_load_sd(&C[(i*56)+8]);
__m128d a18_1 = _mm_load_sd(&values[80]);
#if defined(__SSE3__) && defined(__AVX__)
c18_1 = _mm_add_sd(c18_1, _mm_mul_sd(a18_1, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_1 = _mm_add_sd(c18_1, _mm_mul_sd(a18_1, b18));
#endif
_mm_store_sd(&C[(i*56)+8], c18_1);
__m128d c18_2 = _mm_load_sd(&C[(i*56)+18]);
__m128d a18_2 = _mm_load_sd(&values[81]);
#if defined(__SSE3__) && defined(__AVX__)
c18_2 = _mm_add_sd(c18_2, _mm_mul_sd(a18_2, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_2 = _mm_add_sd(c18_2, _mm_mul_sd(a18_2, b18));
#endif
_mm_store_sd(&C[(i*56)+18], c18_2);
__m128d c18_3 = _mm_load_sd(&C[(i*56)+33]);
__m128d a18_3 = _mm_load_sd(&values[82]);
#if defined(__SSE3__) && defined(__AVX__)
c18_3 = _mm_add_sd(c18_3, _mm_mul_sd(a18_3, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_3 = _mm_add_sd(c18_3, _mm_mul_sd(a18_3, b18));
#endif
_mm_store_sd(&C[(i*56)+33], c18_3);
__m128d c18_4 = _mm_load_sd(&C[(i*56)+54]);
__m128d a18_4 = _mm_load_sd(&values[83]);
#if defined(__SSE3__) && defined(__AVX__)
c18_4 = _mm_add_sd(c18_4, _mm_mul_sd(a18_4, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_4 = _mm_add_sd(c18_4, _mm_mul_sd(a18_4, b18));
#endif
_mm_store_sd(&C[(i*56)+54], c18_4);
#else
C[(i*56)+2] += values[79] * B[(i*56)+18];
C[(i*56)+8] += values[80] * B[(i*56)+18];
C[(i*56)+18] += values[81] * B[(i*56)+18];
C[(i*56)+33] += values[82] * B[(i*56)+18];
C[(i*56)+54] += values[83] * B[(i*56)+18];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b19 = _mm256_broadcast_sd(&B[(i*56)+19]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b19 = _mm_loaddup_pd(&B[(i*56)+19]);
#endif
__m128d c19_0 = _mm_load_sd(&C[(i*56)+0]);
__m128d a19_0 = _mm_load_sd(&values[84]);
#if defined(__SSE3__) && defined(__AVX__)
c19_0 = _mm_add_sd(c19_0, _mm_mul_sd(a19_0, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_0 = _mm_add_sd(c19_0, _mm_mul_sd(a19_0, b19));
#endif
_mm_store_sd(&C[(i*56)+0], c19_0);
__m128d c19_1 = _mm_load_sd(&C[(i*56)+3]);
__m128d a19_1 = _mm_load_sd(&values[85]);
#if defined(__SSE3__) && defined(__AVX__)
c19_1 = _mm_add_sd(c19_1, _mm_mul_sd(a19_1, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_1 = _mm_add_sd(c19_1, _mm_mul_sd(a19_1, b19));
#endif
_mm_store_sd(&C[(i*56)+3], c19_1);
__m128d c19_2 = _mm_load_sd(&C[(i*56)+9]);
__m128d a19_2 = _mm_load_sd(&values[86]);
#if defined(__SSE3__) && defined(__AVX__)
c19_2 = _mm_add_sd(c19_2, _mm_mul_sd(a19_2, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_2 = _mm_add_sd(c19_2, _mm_mul_sd(a19_2, b19));
#endif
_mm_store_sd(&C[(i*56)+9], c19_2);
__m128d c19_3 = _mm_load_sd(&C[(i*56)+19]);
__m128d a19_3 = _mm_load_sd(&values[87]);
#if defined(__SSE3__) && defined(__AVX__)
c19_3 = _mm_add_sd(c19_3, _mm_mul_sd(a19_3, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_3 = _mm_add_sd(c19_3, _mm_mul_sd(a19_3, b19));
#endif
_mm_store_sd(&C[(i*56)+19], c19_3);
__m128d c19_4 = _mm_load_sd(&C[(i*56)+34]);
__m128d a19_4 = _mm_load_sd(&values[88]);
#if defined(__SSE3__) && defined(__AVX__)
c19_4 = _mm_add_sd(c19_4, _mm_mul_sd(a19_4, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_4 = _mm_add_sd(c19_4, _mm_mul_sd(a19_4, b19));
#endif
_mm_store_sd(&C[(i*56)+34], c19_4);
__m128d c19_5 = _mm_load_sd(&C[(i*56)+55]);
__m128d a19_5 = _mm_load_sd(&values[89]);
#if defined(__SSE3__) && defined(__AVX__)
c19_5 = _mm_add_sd(c19_5, _mm_mul_sd(a19_5, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_5 = _mm_add_sd(c19_5, _mm_mul_sd(a19_5, b19));
#endif
_mm_store_sd(&C[(i*56)+55], c19_5);
#else
C[(i*56)+0] += values[84] * B[(i*56)+19];
C[(i*56)+3] += values[85] * B[(i*56)+19];
C[(i*56)+9] += values[86] * B[(i*56)+19];
C[(i*56)+19] += values[87] * B[(i*56)+19];
C[(i*56)+34] += values[88] * B[(i*56)+19];
C[(i*56)+55] += values[89] * B[(i*56)+19];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b20 = _mm256_broadcast_sd(&B[(i*56)+20]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b20 = _mm_loaddup_pd(&B[(i*56)+20]);
#endif
__m128d c20_0 = _mm_load_sd(&C[(i*56)+20]);
__m128d a20_0 = _mm_load_sd(&values[90]);
#if defined(__SSE3__) && defined(__AVX__)
c20_0 = _mm_add_sd(c20_0, _mm_mul_sd(a20_0, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c20_0 = _mm_add_sd(c20_0, _mm_mul_sd(a20_0, b20));
#endif
_mm_store_sd(&C[(i*56)+20], c20_0);
__m128d c20_1 = _mm_load_sd(&C[(i*56)+41]);
__m128d a20_1 = _mm_load_sd(&values[91]);
#if defined(__SSE3__) && defined(__AVX__)
c20_1 = _mm_add_sd(c20_1, _mm_mul_sd(a20_1, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c20_1 = _mm_add_sd(c20_1, _mm_mul_sd(a20_1, b20));
#endif
_mm_store_sd(&C[(i*56)+41], c20_1);
#else
C[(i*56)+20] += values[90] * B[(i*56)+20];
C[(i*56)+41] += values[91] * B[(i*56)+20];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b21 = _mm256_broadcast_sd(&B[(i*56)+21]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b21 = _mm_loaddup_pd(&B[(i*56)+21]);
#endif
__m128d c21_0 = _mm_load_sd(&C[(i*56)+21]);
__m128d a21_0 = _mm_load_sd(&values[92]);
#if defined(__SSE3__) && defined(__AVX__)
c21_0 = _mm_add_sd(c21_0, _mm_mul_sd(a21_0, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c21_0 = _mm_add_sd(c21_0, _mm_mul_sd(a21_0, b21));
#endif
_mm_store_sd(&C[(i*56)+21], c21_0);
__m128d c21_1 = _mm_load_sd(&C[(i*56)+42]);
__m128d a21_1 = _mm_load_sd(&values[93]);
#if defined(__SSE3__) && defined(__AVX__)
c21_1 = _mm_add_sd(c21_1, _mm_mul_sd(a21_1, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c21_1 = _mm_add_sd(c21_1, _mm_mul_sd(a21_1, b21));
#endif
_mm_store_sd(&C[(i*56)+42], c21_1);
#else
C[(i*56)+21] += values[92] * B[(i*56)+21];
C[(i*56)+42] += values[93] * B[(i*56)+21];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b22 = _mm256_broadcast_sd(&B[(i*56)+22]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b22 = _mm_loaddup_pd(&B[(i*56)+22]);
#endif
__m128d c22_0 = _mm_load_sd(&C[(i*56)+22]);
__m128d a22_0 = _mm_load_sd(&values[94]);
#if defined(__SSE3__) && defined(__AVX__)
c22_0 = _mm_add_sd(c22_0, _mm_mul_sd(a22_0, _mm256_castpd256_pd128(b22)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c22_0 = _mm_add_sd(c22_0, _mm_mul_sd(a22_0, b22));
#endif
_mm_store_sd(&C[(i*56)+22], c22_0);
__m128d c22_1 = _mm_load_sd(&C[(i*56)+43]);
__m128d a22_1 = _mm_load_sd(&values[95]);
#if defined(__SSE3__) && defined(__AVX__)
c22_1 = _mm_add_sd(c22_1, _mm_mul_sd(a22_1, _mm256_castpd256_pd128(b22)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c22_1 = _mm_add_sd(c22_1, _mm_mul_sd(a22_1, b22));
#endif
_mm_store_sd(&C[(i*56)+43], c22_1);
#else
C[(i*56)+22] += values[94] * B[(i*56)+22];
C[(i*56)+43] += values[95] * B[(i*56)+22];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b23 = _mm256_broadcast_sd(&B[(i*56)+23]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b23 = _mm_loaddup_pd(&B[(i*56)+23]);
#endif
__m128d c23_0 = _mm_load_sd(&C[(i*56)+23]);
__m128d a23_0 = _mm_load_sd(&values[96]);
#if defined(__SSE3__) && defined(__AVX__)
c23_0 = _mm_add_sd(c23_0, _mm_mul_sd(a23_0, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c23_0 = _mm_add_sd(c23_0, _mm_mul_sd(a23_0, b23));
#endif
_mm_store_sd(&C[(i*56)+23], c23_0);
__m128d c23_1 = _mm_load_sd(&C[(i*56)+44]);
__m128d a23_1 = _mm_load_sd(&values[97]);
#if defined(__SSE3__) && defined(__AVX__)
c23_1 = _mm_add_sd(c23_1, _mm_mul_sd(a23_1, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c23_1 = _mm_add_sd(c23_1, _mm_mul_sd(a23_1, b23));
#endif
_mm_store_sd(&C[(i*56)+44], c23_1);
#else
C[(i*56)+23] += values[96] * B[(i*56)+23];
C[(i*56)+44] += values[97] * B[(i*56)+23];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b24 = _mm256_broadcast_sd(&B[(i*56)+24]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b24 = _mm_loaddup_pd(&B[(i*56)+24]);
#endif
__m128d c24_0 = _mm_load_sd(&C[(i*56)+24]);
__m128d a24_0 = _mm_load_sd(&values[98]);
#if defined(__SSE3__) && defined(__AVX__)
c24_0 = _mm_add_sd(c24_0, _mm_mul_sd(a24_0, _mm256_castpd256_pd128(b24)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c24_0 = _mm_add_sd(c24_0, _mm_mul_sd(a24_0, b24));
#endif
_mm_store_sd(&C[(i*56)+24], c24_0);
__m128d c24_1 = _mm_load_sd(&C[(i*56)+45]);
__m128d a24_1 = _mm_load_sd(&values[99]);
#if defined(__SSE3__) && defined(__AVX__)
c24_1 = _mm_add_sd(c24_1, _mm_mul_sd(a24_1, _mm256_castpd256_pd128(b24)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c24_1 = _mm_add_sd(c24_1, _mm_mul_sd(a24_1, b24));
#endif
_mm_store_sd(&C[(i*56)+45], c24_1);
#else
C[(i*56)+24] += values[98] * B[(i*56)+24];
C[(i*56)+45] += values[99] * B[(i*56)+24];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b25 = _mm256_broadcast_sd(&B[(i*56)+25]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b25 = _mm_loaddup_pd(&B[(i*56)+25]);
#endif
__m128d c25_0 = _mm_load_sd(&C[(i*56)+10]);
__m128d a25_0 = _mm_load_sd(&values[100]);
#if defined(__SSE3__) && defined(__AVX__)
c25_0 = _mm_add_sd(c25_0, _mm_mul_sd(a25_0, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c25_0 = _mm_add_sd(c25_0, _mm_mul_sd(a25_0, b25));
#endif
_mm_store_sd(&C[(i*56)+10], c25_0);
__m128d c25_1 = _mm_load_sd(&C[(i*56)+25]);
__m128d a25_1 = _mm_load_sd(&values[101]);
#if defined(__SSE3__) && defined(__AVX__)
c25_1 = _mm_add_sd(c25_1, _mm_mul_sd(a25_1, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c25_1 = _mm_add_sd(c25_1, _mm_mul_sd(a25_1, b25));
#endif
_mm_store_sd(&C[(i*56)+25], c25_1);
__m128d c25_2 = _mm_load_sd(&C[(i*56)+46]);
__m128d a25_2 = _mm_load_sd(&values[102]);
#if defined(__SSE3__) && defined(__AVX__)
c25_2 = _mm_add_sd(c25_2, _mm_mul_sd(a25_2, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c25_2 = _mm_add_sd(c25_2, _mm_mul_sd(a25_2, b25));
#endif
_mm_store_sd(&C[(i*56)+46], c25_2);
#else
C[(i*56)+10] += values[100] * B[(i*56)+25];
C[(i*56)+25] += values[101] * B[(i*56)+25];
C[(i*56)+46] += values[102] * B[(i*56)+25];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b26 = _mm256_broadcast_sd(&B[(i*56)+26]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b26 = _mm_loaddup_pd(&B[(i*56)+26]);
#endif
__m128d c26_0 = _mm_load_sd(&C[(i*56)+11]);
__m128d a26_0 = _mm_load_sd(&values[103]);
#if defined(__SSE3__) && defined(__AVX__)
c26_0 = _mm_add_sd(c26_0, _mm_mul_sd(a26_0, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c26_0 = _mm_add_sd(c26_0, _mm_mul_sd(a26_0, b26));
#endif
_mm_store_sd(&C[(i*56)+11], c26_0);
__m128d c26_1 = _mm_load_sd(&C[(i*56)+26]);
__m128d a26_1 = _mm_load_sd(&values[104]);
#if defined(__SSE3__) && defined(__AVX__)
c26_1 = _mm_add_sd(c26_1, _mm_mul_sd(a26_1, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c26_1 = _mm_add_sd(c26_1, _mm_mul_sd(a26_1, b26));
#endif
_mm_store_sd(&C[(i*56)+26], c26_1);
__m128d c26_2 = _mm_load_sd(&C[(i*56)+47]);
__m128d a26_2 = _mm_load_sd(&values[105]);
#if defined(__SSE3__) && defined(__AVX__)
c26_2 = _mm_add_sd(c26_2, _mm_mul_sd(a26_2, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c26_2 = _mm_add_sd(c26_2, _mm_mul_sd(a26_2, b26));
#endif
_mm_store_sd(&C[(i*56)+47], c26_2);
#else
C[(i*56)+11] += values[103] * B[(i*56)+26];
C[(i*56)+26] += values[104] * B[(i*56)+26];
C[(i*56)+47] += values[105] * B[(i*56)+26];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b27 = _mm256_broadcast_sd(&B[(i*56)+27]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b27 = _mm_loaddup_pd(&B[(i*56)+27]);
#endif
__m128d c27_0 = _mm_load_sd(&C[(i*56)+12]);
__m128d a27_0 = _mm_load_sd(&values[106]);
#if defined(__SSE3__) && defined(__AVX__)
c27_0 = _mm_add_sd(c27_0, _mm_mul_sd(a27_0, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c27_0 = _mm_add_sd(c27_0, _mm_mul_sd(a27_0, b27));
#endif
_mm_store_sd(&C[(i*56)+12], c27_0);
__m128d c27_1 = _mm_load_sd(&C[(i*56)+27]);
__m128d a27_1 = _mm_load_sd(&values[107]);
#if defined(__SSE3__) && defined(__AVX__)
c27_1 = _mm_add_sd(c27_1, _mm_mul_sd(a27_1, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c27_1 = _mm_add_sd(c27_1, _mm_mul_sd(a27_1, b27));
#endif
_mm_store_sd(&C[(i*56)+27], c27_1);
__m128d c27_2 = _mm_load_sd(&C[(i*56)+48]);
__m128d a27_2 = _mm_load_sd(&values[108]);
#if defined(__SSE3__) && defined(__AVX__)
c27_2 = _mm_add_sd(c27_2, _mm_mul_sd(a27_2, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c27_2 = _mm_add_sd(c27_2, _mm_mul_sd(a27_2, b27));
#endif
_mm_store_sd(&C[(i*56)+48], c27_2);
#else
C[(i*56)+12] += values[106] * B[(i*56)+27];
C[(i*56)+27] += values[107] * B[(i*56)+27];
C[(i*56)+48] += values[108] * B[(i*56)+27];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b28 = _mm256_broadcast_sd(&B[(i*56)+28]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b28 = _mm_loaddup_pd(&B[(i*56)+28]);
#endif
__m128d c28_0 = _mm_load_sd(&C[(i*56)+13]);
__m128d a28_0 = _mm_load_sd(&values[109]);
#if defined(__SSE3__) && defined(__AVX__)
c28_0 = _mm_add_sd(c28_0, _mm_mul_sd(a28_0, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c28_0 = _mm_add_sd(c28_0, _mm_mul_sd(a28_0, b28));
#endif
_mm_store_sd(&C[(i*56)+13], c28_0);
__m128d c28_1 = _mm_load_sd(&C[(i*56)+28]);
__m128d a28_1 = _mm_load_sd(&values[110]);
#if defined(__SSE3__) && defined(__AVX__)
c28_1 = _mm_add_sd(c28_1, _mm_mul_sd(a28_1, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c28_1 = _mm_add_sd(c28_1, _mm_mul_sd(a28_1, b28));
#endif
_mm_store_sd(&C[(i*56)+28], c28_1);
__m128d c28_2 = _mm_load_sd(&C[(i*56)+49]);
__m128d a28_2 = _mm_load_sd(&values[111]);
#if defined(__SSE3__) && defined(__AVX__)
c28_2 = _mm_add_sd(c28_2, _mm_mul_sd(a28_2, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c28_2 = _mm_add_sd(c28_2, _mm_mul_sd(a28_2, b28));
#endif
_mm_store_sd(&C[(i*56)+49], c28_2);
#else
C[(i*56)+13] += values[109] * B[(i*56)+28];
C[(i*56)+28] += values[110] * B[(i*56)+28];
C[(i*56)+49] += values[111] * B[(i*56)+28];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b29 = _mm256_broadcast_sd(&B[(i*56)+29]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b29 = _mm_loaddup_pd(&B[(i*56)+29]);
#endif
__m128d c29_0 = _mm_load_sd(&C[(i*56)+4]);
__m128d a29_0 = _mm_load_sd(&values[112]);
#if defined(__SSE3__) && defined(__AVX__)
c29_0 = _mm_add_sd(c29_0, _mm_mul_sd(a29_0, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_0 = _mm_add_sd(c29_0, _mm_mul_sd(a29_0, b29));
#endif
_mm_store_sd(&C[(i*56)+4], c29_0);
__m128d c29_1 = _mm_load_sd(&C[(i*56)+14]);
__m128d a29_1 = _mm_load_sd(&values[113]);
#if defined(__SSE3__) && defined(__AVX__)
c29_1 = _mm_add_sd(c29_1, _mm_mul_sd(a29_1, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_1 = _mm_add_sd(c29_1, _mm_mul_sd(a29_1, b29));
#endif
_mm_store_sd(&C[(i*56)+14], c29_1);
__m128d c29_2 = _mm_load_sd(&C[(i*56)+29]);
__m128d a29_2 = _mm_load_sd(&values[114]);
#if defined(__SSE3__) && defined(__AVX__)
c29_2 = _mm_add_sd(c29_2, _mm_mul_sd(a29_2, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_2 = _mm_add_sd(c29_2, _mm_mul_sd(a29_2, b29));
#endif
_mm_store_sd(&C[(i*56)+29], c29_2);
__m128d c29_3 = _mm_load_sd(&C[(i*56)+50]);
__m128d a29_3 = _mm_load_sd(&values[115]);
#if defined(__SSE3__) && defined(__AVX__)
c29_3 = _mm_add_sd(c29_3, _mm_mul_sd(a29_3, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_3 = _mm_add_sd(c29_3, _mm_mul_sd(a29_3, b29));
#endif
_mm_store_sd(&C[(i*56)+50], c29_3);
#else
C[(i*56)+4] += values[112] * B[(i*56)+29];
C[(i*56)+14] += values[113] * B[(i*56)+29];
C[(i*56)+29] += values[114] * B[(i*56)+29];
C[(i*56)+50] += values[115] * B[(i*56)+29];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b30 = _mm256_broadcast_sd(&B[(i*56)+30]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b30 = _mm_loaddup_pd(&B[(i*56)+30]);
#endif
__m128d c30_0 = _mm_load_sd(&C[(i*56)+5]);
__m128d a30_0 = _mm_load_sd(&values[116]);
#if defined(__SSE3__) && defined(__AVX__)
c30_0 = _mm_add_sd(c30_0, _mm_mul_sd(a30_0, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_0 = _mm_add_sd(c30_0, _mm_mul_sd(a30_0, b30));
#endif
_mm_store_sd(&C[(i*56)+5], c30_0);
__m128d c30_1 = _mm_load_sd(&C[(i*56)+15]);
__m128d a30_1 = _mm_load_sd(&values[117]);
#if defined(__SSE3__) && defined(__AVX__)
c30_1 = _mm_add_sd(c30_1, _mm_mul_sd(a30_1, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_1 = _mm_add_sd(c30_1, _mm_mul_sd(a30_1, b30));
#endif
_mm_store_sd(&C[(i*56)+15], c30_1);
__m128d c30_2 = _mm_load_sd(&C[(i*56)+30]);
__m128d a30_2 = _mm_load_sd(&values[118]);
#if defined(__SSE3__) && defined(__AVX__)
c30_2 = _mm_add_sd(c30_2, _mm_mul_sd(a30_2, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_2 = _mm_add_sd(c30_2, _mm_mul_sd(a30_2, b30));
#endif
_mm_store_sd(&C[(i*56)+30], c30_2);
__m128d c30_3 = _mm_load_sd(&C[(i*56)+51]);
__m128d a30_3 = _mm_load_sd(&values[119]);
#if defined(__SSE3__) && defined(__AVX__)
c30_3 = _mm_add_sd(c30_3, _mm_mul_sd(a30_3, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_3 = _mm_add_sd(c30_3, _mm_mul_sd(a30_3, b30));
#endif
_mm_store_sd(&C[(i*56)+51], c30_3);
#else
C[(i*56)+5] += values[116] * B[(i*56)+30];
C[(i*56)+15] += values[117] * B[(i*56)+30];
C[(i*56)+30] += values[118] * B[(i*56)+30];
C[(i*56)+51] += values[119] * B[(i*56)+30];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b31 = _mm256_broadcast_sd(&B[(i*56)+31]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b31 = _mm_loaddup_pd(&B[(i*56)+31]);
#endif
__m128d c31_0 = _mm_load_sd(&C[(i*56)+6]);
__m128d a31_0 = _mm_load_sd(&values[120]);
#if defined(__SSE3__) && defined(__AVX__)
c31_0 = _mm_add_sd(c31_0, _mm_mul_sd(a31_0, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_0 = _mm_add_sd(c31_0, _mm_mul_sd(a31_0, b31));
#endif
_mm_store_sd(&C[(i*56)+6], c31_0);
__m128d c31_1 = _mm_load_sd(&C[(i*56)+16]);
__m128d a31_1 = _mm_load_sd(&values[121]);
#if defined(__SSE3__) && defined(__AVX__)
c31_1 = _mm_add_sd(c31_1, _mm_mul_sd(a31_1, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_1 = _mm_add_sd(c31_1, _mm_mul_sd(a31_1, b31));
#endif
_mm_store_sd(&C[(i*56)+16], c31_1);
__m128d c31_2 = _mm_load_sd(&C[(i*56)+31]);
__m128d a31_2 = _mm_load_sd(&values[122]);
#if defined(__SSE3__) && defined(__AVX__)
c31_2 = _mm_add_sd(c31_2, _mm_mul_sd(a31_2, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_2 = _mm_add_sd(c31_2, _mm_mul_sd(a31_2, b31));
#endif
_mm_store_sd(&C[(i*56)+31], c31_2);
__m128d c31_3 = _mm_load_sd(&C[(i*56)+52]);
__m128d a31_3 = _mm_load_sd(&values[123]);
#if defined(__SSE3__) && defined(__AVX__)
c31_3 = _mm_add_sd(c31_3, _mm_mul_sd(a31_3, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_3 = _mm_add_sd(c31_3, _mm_mul_sd(a31_3, b31));
#endif
_mm_store_sd(&C[(i*56)+52], c31_3);
#else
C[(i*56)+6] += values[120] * B[(i*56)+31];
C[(i*56)+16] += values[121] * B[(i*56)+31];
C[(i*56)+31] += values[122] * B[(i*56)+31];
C[(i*56)+52] += values[123] * B[(i*56)+31];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b32 = _mm256_broadcast_sd(&B[(i*56)+32]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b32 = _mm_loaddup_pd(&B[(i*56)+32]);
#endif
__m128d c32_0 = _mm_load_sd(&C[(i*56)+1]);
__m128d a32_0 = _mm_load_sd(&values[124]);
#if defined(__SSE3__) && defined(__AVX__)
c32_0 = _mm_add_sd(c32_0, _mm_mul_sd(a32_0, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_0 = _mm_add_sd(c32_0, _mm_mul_sd(a32_0, b32));
#endif
_mm_store_sd(&C[(i*56)+1], c32_0);
__m128d c32_1 = _mm_load_sd(&C[(i*56)+7]);
__m128d a32_1 = _mm_load_sd(&values[125]);
#if defined(__SSE3__) && defined(__AVX__)
c32_1 = _mm_add_sd(c32_1, _mm_mul_sd(a32_1, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_1 = _mm_add_sd(c32_1, _mm_mul_sd(a32_1, b32));
#endif
_mm_store_sd(&C[(i*56)+7], c32_1);
__m128d c32_2 = _mm_load_sd(&C[(i*56)+17]);
__m128d a32_2 = _mm_load_sd(&values[126]);
#if defined(__SSE3__) && defined(__AVX__)
c32_2 = _mm_add_sd(c32_2, _mm_mul_sd(a32_2, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_2 = _mm_add_sd(c32_2, _mm_mul_sd(a32_2, b32));
#endif
_mm_store_sd(&C[(i*56)+17], c32_2);
__m128d c32_3 = _mm_load_sd(&C[(i*56)+32]);
__m128d a32_3 = _mm_load_sd(&values[127]);
#if defined(__SSE3__) && defined(__AVX__)
c32_3 = _mm_add_sd(c32_3, _mm_mul_sd(a32_3, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_3 = _mm_add_sd(c32_3, _mm_mul_sd(a32_3, b32));
#endif
_mm_store_sd(&C[(i*56)+32], c32_3);
__m128d c32_4 = _mm_load_sd(&C[(i*56)+53]);
__m128d a32_4 = _mm_load_sd(&values[128]);
#if defined(__SSE3__) && defined(__AVX__)
c32_4 = _mm_add_sd(c32_4, _mm_mul_sd(a32_4, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_4 = _mm_add_sd(c32_4, _mm_mul_sd(a32_4, b32));
#endif
_mm_store_sd(&C[(i*56)+53], c32_4);
#else
C[(i*56)+1] += values[124] * B[(i*56)+32];
C[(i*56)+7] += values[125] * B[(i*56)+32];
C[(i*56)+17] += values[126] * B[(i*56)+32];
C[(i*56)+32] += values[127] * B[(i*56)+32];
C[(i*56)+53] += values[128] * B[(i*56)+32];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b33 = _mm256_broadcast_sd(&B[(i*56)+33]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b33 = _mm_loaddup_pd(&B[(i*56)+33]);
#endif
__m128d c33_0 = _mm_load_sd(&C[(i*56)+2]);
__m128d a33_0 = _mm_load_sd(&values[129]);
#if defined(__SSE3__) && defined(__AVX__)
c33_0 = _mm_add_sd(c33_0, _mm_mul_sd(a33_0, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_0 = _mm_add_sd(c33_0, _mm_mul_sd(a33_0, b33));
#endif
_mm_store_sd(&C[(i*56)+2], c33_0);
__m128d c33_1 = _mm_load_sd(&C[(i*56)+8]);
__m128d a33_1 = _mm_load_sd(&values[130]);
#if defined(__SSE3__) && defined(__AVX__)
c33_1 = _mm_add_sd(c33_1, _mm_mul_sd(a33_1, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_1 = _mm_add_sd(c33_1, _mm_mul_sd(a33_1, b33));
#endif
_mm_store_sd(&C[(i*56)+8], c33_1);
__m128d c33_2 = _mm_load_sd(&C[(i*56)+18]);
__m128d a33_2 = _mm_load_sd(&values[131]);
#if defined(__SSE3__) && defined(__AVX__)
c33_2 = _mm_add_sd(c33_2, _mm_mul_sd(a33_2, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_2 = _mm_add_sd(c33_2, _mm_mul_sd(a33_2, b33));
#endif
_mm_store_sd(&C[(i*56)+18], c33_2);
__m128d c33_3 = _mm_load_sd(&C[(i*56)+33]);
__m128d a33_3 = _mm_load_sd(&values[132]);
#if defined(__SSE3__) && defined(__AVX__)
c33_3 = _mm_add_sd(c33_3, _mm_mul_sd(a33_3, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_3 = _mm_add_sd(c33_3, _mm_mul_sd(a33_3, b33));
#endif
_mm_store_sd(&C[(i*56)+33], c33_3);
__m128d c33_4 = _mm_load_sd(&C[(i*56)+54]);
__m128d a33_4 = _mm_load_sd(&values[133]);
#if defined(__SSE3__) && defined(__AVX__)
c33_4 = _mm_add_sd(c33_4, _mm_mul_sd(a33_4, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_4 = _mm_add_sd(c33_4, _mm_mul_sd(a33_4, b33));
#endif
_mm_store_sd(&C[(i*56)+54], c33_4);
#else
C[(i*56)+2] += values[129] * B[(i*56)+33];
C[(i*56)+8] += values[130] * B[(i*56)+33];
C[(i*56)+18] += values[131] * B[(i*56)+33];
C[(i*56)+33] += values[132] * B[(i*56)+33];
C[(i*56)+54] += values[133] * B[(i*56)+33];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b34 = _mm256_broadcast_sd(&B[(i*56)+34]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b34 = _mm_loaddup_pd(&B[(i*56)+34]);
#endif
__m128d c34_0 = _mm_load_sd(&C[(i*56)+0]);
__m128d a34_0 = _mm_load_sd(&values[134]);
#if defined(__SSE3__) && defined(__AVX__)
c34_0 = _mm_add_sd(c34_0, _mm_mul_sd(a34_0, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_0 = _mm_add_sd(c34_0, _mm_mul_sd(a34_0, b34));
#endif
_mm_store_sd(&C[(i*56)+0], c34_0);
__m128d c34_1 = _mm_load_sd(&C[(i*56)+3]);
__m128d a34_1 = _mm_load_sd(&values[135]);
#if defined(__SSE3__) && defined(__AVX__)
c34_1 = _mm_add_sd(c34_1, _mm_mul_sd(a34_1, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_1 = _mm_add_sd(c34_1, _mm_mul_sd(a34_1, b34));
#endif
_mm_store_sd(&C[(i*56)+3], c34_1);
__m128d c34_2 = _mm_load_sd(&C[(i*56)+9]);
__m128d a34_2 = _mm_load_sd(&values[136]);
#if defined(__SSE3__) && defined(__AVX__)
c34_2 = _mm_add_sd(c34_2, _mm_mul_sd(a34_2, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_2 = _mm_add_sd(c34_2, _mm_mul_sd(a34_2, b34));
#endif
_mm_store_sd(&C[(i*56)+9], c34_2);
__m128d c34_3 = _mm_load_sd(&C[(i*56)+19]);
__m128d a34_3 = _mm_load_sd(&values[137]);
#if defined(__SSE3__) && defined(__AVX__)
c34_3 = _mm_add_sd(c34_3, _mm_mul_sd(a34_3, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_3 = _mm_add_sd(c34_3, _mm_mul_sd(a34_3, b34));
#endif
_mm_store_sd(&C[(i*56)+19], c34_3);
__m128d c34_4 = _mm_load_sd(&C[(i*56)+34]);
__m128d a34_4 = _mm_load_sd(&values[138]);
#if defined(__SSE3__) && defined(__AVX__)
c34_4 = _mm_add_sd(c34_4, _mm_mul_sd(a34_4, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_4 = _mm_add_sd(c34_4, _mm_mul_sd(a34_4, b34));
#endif
_mm_store_sd(&C[(i*56)+34], c34_4);
__m128d c34_5 = _mm_load_sd(&C[(i*56)+55]);
__m128d a34_5 = _mm_load_sd(&values[139]);
#if defined(__SSE3__) && defined(__AVX__)
c34_5 = _mm_add_sd(c34_5, _mm_mul_sd(a34_5, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_5 = _mm_add_sd(c34_5, _mm_mul_sd(a34_5, b34));
#endif
_mm_store_sd(&C[(i*56)+55], c34_5);
#else
C[(i*56)+0] += values[134] * B[(i*56)+34];
C[(i*56)+3] += values[135] * B[(i*56)+34];
C[(i*56)+9] += values[136] * B[(i*56)+34];
C[(i*56)+19] += values[137] * B[(i*56)+34];
C[(i*56)+34] += values[138] * B[(i*56)+34];
C[(i*56)+55] += values[139] * B[(i*56)+34];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b35 = _mm256_broadcast_sd(&B[(i*56)+35]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b35 = _mm_loaddup_pd(&B[(i*56)+35]);
#endif
__m128d c35_0 = _mm_load_sd(&C[(i*56)+35]);
__m128d a35_0 = _mm_load_sd(&values[140]);
#if defined(__SSE3__) && defined(__AVX__)
c35_0 = _mm_add_sd(c35_0, _mm_mul_sd(a35_0, _mm256_castpd256_pd128(b35)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c35_0 = _mm_add_sd(c35_0, _mm_mul_sd(a35_0, b35));
#endif
_mm_store_sd(&C[(i*56)+35], c35_0);
#else
C[(i*56)+35] += values[140] * B[(i*56)+35];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b36 = _mm256_broadcast_sd(&B[(i*56)+36]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b36 = _mm_loaddup_pd(&B[(i*56)+36]);
#endif
__m128d c36_0 = _mm_load_sd(&C[(i*56)+36]);
__m128d a36_0 = _mm_load_sd(&values[141]);
#if defined(__SSE3__) && defined(__AVX__)
c36_0 = _mm_add_sd(c36_0, _mm_mul_sd(a36_0, _mm256_castpd256_pd128(b36)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c36_0 = _mm_add_sd(c36_0, _mm_mul_sd(a36_0, b36));
#endif
_mm_store_sd(&C[(i*56)+36], c36_0);
#else
C[(i*56)+36] += values[141] * B[(i*56)+36];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b37 = _mm256_broadcast_sd(&B[(i*56)+37]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b37 = _mm_loaddup_pd(&B[(i*56)+37]);
#endif
__m128d c37_0 = _mm_load_sd(&C[(i*56)+37]);
__m128d a37_0 = _mm_load_sd(&values[142]);
#if defined(__SSE3__) && defined(__AVX__)
c37_0 = _mm_add_sd(c37_0, _mm_mul_sd(a37_0, _mm256_castpd256_pd128(b37)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c37_0 = _mm_add_sd(c37_0, _mm_mul_sd(a37_0, b37));
#endif
_mm_store_sd(&C[(i*56)+37], c37_0);
#else
C[(i*56)+37] += values[142] * B[(i*56)+37];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b38 = _mm256_broadcast_sd(&B[(i*56)+38]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b38 = _mm_loaddup_pd(&B[(i*56)+38]);
#endif
__m128d c38_0 = _mm_load_sd(&C[(i*56)+38]);
__m128d a38_0 = _mm_load_sd(&values[143]);
#if defined(__SSE3__) && defined(__AVX__)
c38_0 = _mm_add_sd(c38_0, _mm_mul_sd(a38_0, _mm256_castpd256_pd128(b38)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c38_0 = _mm_add_sd(c38_0, _mm_mul_sd(a38_0, b38));
#endif
_mm_store_sd(&C[(i*56)+38], c38_0);
#else
C[(i*56)+38] += values[143] * B[(i*56)+38];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b39 = _mm256_broadcast_sd(&B[(i*56)+39]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b39 = _mm_loaddup_pd(&B[(i*56)+39]);
#endif
__m128d c39_0 = _mm_load_sd(&C[(i*56)+39]);
__m128d a39_0 = _mm_load_sd(&values[144]);
#if defined(__SSE3__) && defined(__AVX__)
c39_0 = _mm_add_sd(c39_0, _mm_mul_sd(a39_0, _mm256_castpd256_pd128(b39)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c39_0 = _mm_add_sd(c39_0, _mm_mul_sd(a39_0, b39));
#endif
_mm_store_sd(&C[(i*56)+39], c39_0);
#else
C[(i*56)+39] += values[144] * B[(i*56)+39];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b40 = _mm256_broadcast_sd(&B[(i*56)+40]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b40 = _mm_loaddup_pd(&B[(i*56)+40]);
#endif
__m128d c40_0 = _mm_load_sd(&C[(i*56)+40]);
__m128d a40_0 = _mm_load_sd(&values[145]);
#if defined(__SSE3__) && defined(__AVX__)
c40_0 = _mm_add_sd(c40_0, _mm_mul_sd(a40_0, _mm256_castpd256_pd128(b40)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c40_0 = _mm_add_sd(c40_0, _mm_mul_sd(a40_0, b40));
#endif
_mm_store_sd(&C[(i*56)+40], c40_0);
#else
C[(i*56)+40] += values[145] * B[(i*56)+40];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b41 = _mm256_broadcast_sd(&B[(i*56)+41]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b41 = _mm_loaddup_pd(&B[(i*56)+41]);
#endif
__m128d c41_0 = _mm_load_sd(&C[(i*56)+20]);
__m128d a41_0 = _mm_load_sd(&values[146]);
#if defined(__SSE3__) && defined(__AVX__)
c41_0 = _mm_add_sd(c41_0, _mm_mul_sd(a41_0, _mm256_castpd256_pd128(b41)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c41_0 = _mm_add_sd(c41_0, _mm_mul_sd(a41_0, b41));
#endif
_mm_store_sd(&C[(i*56)+20], c41_0);
__m128d c41_1 = _mm_load_sd(&C[(i*56)+41]);
__m128d a41_1 = _mm_load_sd(&values[147]);
#if defined(__SSE3__) && defined(__AVX__)
c41_1 = _mm_add_sd(c41_1, _mm_mul_sd(a41_1, _mm256_castpd256_pd128(b41)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c41_1 = _mm_add_sd(c41_1, _mm_mul_sd(a41_1, b41));
#endif
_mm_store_sd(&C[(i*56)+41], c41_1);
#else
C[(i*56)+20] += values[146] * B[(i*56)+41];
C[(i*56)+41] += values[147] * B[(i*56)+41];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b42 = _mm256_broadcast_sd(&B[(i*56)+42]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b42 = _mm_loaddup_pd(&B[(i*56)+42]);
#endif
__m128d c42_0 = _mm_load_sd(&C[(i*56)+21]);
__m128d a42_0 = _mm_load_sd(&values[148]);
#if defined(__SSE3__) && defined(__AVX__)
c42_0 = _mm_add_sd(c42_0, _mm_mul_sd(a42_0, _mm256_castpd256_pd128(b42)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c42_0 = _mm_add_sd(c42_0, _mm_mul_sd(a42_0, b42));
#endif
_mm_store_sd(&C[(i*56)+21], c42_0);
__m128d c42_1 = _mm_load_sd(&C[(i*56)+42]);
__m128d a42_1 = _mm_load_sd(&values[149]);
#if defined(__SSE3__) && defined(__AVX__)
c42_1 = _mm_add_sd(c42_1, _mm_mul_sd(a42_1, _mm256_castpd256_pd128(b42)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c42_1 = _mm_add_sd(c42_1, _mm_mul_sd(a42_1, b42));
#endif
_mm_store_sd(&C[(i*56)+42], c42_1);
#else
C[(i*56)+21] += values[148] * B[(i*56)+42];
C[(i*56)+42] += values[149] * B[(i*56)+42];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b43 = _mm256_broadcast_sd(&B[(i*56)+43]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b43 = _mm_loaddup_pd(&B[(i*56)+43]);
#endif
__m128d c43_0 = _mm_load_sd(&C[(i*56)+22]);
__m128d a43_0 = _mm_load_sd(&values[150]);
#if defined(__SSE3__) && defined(__AVX__)
c43_0 = _mm_add_sd(c43_0, _mm_mul_sd(a43_0, _mm256_castpd256_pd128(b43)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c43_0 = _mm_add_sd(c43_0, _mm_mul_sd(a43_0, b43));
#endif
_mm_store_sd(&C[(i*56)+22], c43_0);
__m128d c43_1 = _mm_load_sd(&C[(i*56)+43]);
__m128d a43_1 = _mm_load_sd(&values[151]);
#if defined(__SSE3__) && defined(__AVX__)
c43_1 = _mm_add_sd(c43_1, _mm_mul_sd(a43_1, _mm256_castpd256_pd128(b43)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c43_1 = _mm_add_sd(c43_1, _mm_mul_sd(a43_1, b43));
#endif
_mm_store_sd(&C[(i*56)+43], c43_1);
#else
C[(i*56)+22] += values[150] * B[(i*56)+43];
C[(i*56)+43] += values[151] * B[(i*56)+43];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b44 = _mm256_broadcast_sd(&B[(i*56)+44]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b44 = _mm_loaddup_pd(&B[(i*56)+44]);
#endif
__m128d c44_0 = _mm_load_sd(&C[(i*56)+23]);
__m128d a44_0 = _mm_load_sd(&values[152]);
#if defined(__SSE3__) && defined(__AVX__)
c44_0 = _mm_add_sd(c44_0, _mm_mul_sd(a44_0, _mm256_castpd256_pd128(b44)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c44_0 = _mm_add_sd(c44_0, _mm_mul_sd(a44_0, b44));
#endif
_mm_store_sd(&C[(i*56)+23], c44_0);
__m128d c44_1 = _mm_load_sd(&C[(i*56)+44]);
__m128d a44_1 = _mm_load_sd(&values[153]);
#if defined(__SSE3__) && defined(__AVX__)
c44_1 = _mm_add_sd(c44_1, _mm_mul_sd(a44_1, _mm256_castpd256_pd128(b44)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c44_1 = _mm_add_sd(c44_1, _mm_mul_sd(a44_1, b44));
#endif
_mm_store_sd(&C[(i*56)+44], c44_1);
#else
C[(i*56)+23] += values[152] * B[(i*56)+44];
C[(i*56)+44] += values[153] * B[(i*56)+44];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b45 = _mm256_broadcast_sd(&B[(i*56)+45]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b45 = _mm_loaddup_pd(&B[(i*56)+45]);
#endif
__m128d c45_0 = _mm_load_sd(&C[(i*56)+24]);
__m128d a45_0 = _mm_load_sd(&values[154]);
#if defined(__SSE3__) && defined(__AVX__)
c45_0 = _mm_add_sd(c45_0, _mm_mul_sd(a45_0, _mm256_castpd256_pd128(b45)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c45_0 = _mm_add_sd(c45_0, _mm_mul_sd(a45_0, b45));
#endif
_mm_store_sd(&C[(i*56)+24], c45_0);
__m128d c45_1 = _mm_load_sd(&C[(i*56)+45]);
__m128d a45_1 = _mm_load_sd(&values[155]);
#if defined(__SSE3__) && defined(__AVX__)
c45_1 = _mm_add_sd(c45_1, _mm_mul_sd(a45_1, _mm256_castpd256_pd128(b45)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c45_1 = _mm_add_sd(c45_1, _mm_mul_sd(a45_1, b45));
#endif
_mm_store_sd(&C[(i*56)+45], c45_1);
#else
C[(i*56)+24] += values[154] * B[(i*56)+45];
C[(i*56)+45] += values[155] * B[(i*56)+45];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b46 = _mm256_broadcast_sd(&B[(i*56)+46]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b46 = _mm_loaddup_pd(&B[(i*56)+46]);
#endif
__m128d c46_0 = _mm_load_sd(&C[(i*56)+10]);
__m128d a46_0 = _mm_load_sd(&values[156]);
#if defined(__SSE3__) && defined(__AVX__)
c46_0 = _mm_add_sd(c46_0, _mm_mul_sd(a46_0, _mm256_castpd256_pd128(b46)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c46_0 = _mm_add_sd(c46_0, _mm_mul_sd(a46_0, b46));
#endif
_mm_store_sd(&C[(i*56)+10], c46_0);
__m128d c46_1 = _mm_load_sd(&C[(i*56)+25]);
__m128d a46_1 = _mm_load_sd(&values[157]);
#if defined(__SSE3__) && defined(__AVX__)
c46_1 = _mm_add_sd(c46_1, _mm_mul_sd(a46_1, _mm256_castpd256_pd128(b46)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c46_1 = _mm_add_sd(c46_1, _mm_mul_sd(a46_1, b46));
#endif
_mm_store_sd(&C[(i*56)+25], c46_1);
__m128d c46_2 = _mm_load_sd(&C[(i*56)+46]);
__m128d a46_2 = _mm_load_sd(&values[158]);
#if defined(__SSE3__) && defined(__AVX__)
c46_2 = _mm_add_sd(c46_2, _mm_mul_sd(a46_2, _mm256_castpd256_pd128(b46)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c46_2 = _mm_add_sd(c46_2, _mm_mul_sd(a46_2, b46));
#endif
_mm_store_sd(&C[(i*56)+46], c46_2);
#else
C[(i*56)+10] += values[156] * B[(i*56)+46];
C[(i*56)+25] += values[157] * B[(i*56)+46];
C[(i*56)+46] += values[158] * B[(i*56)+46];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b47 = _mm256_broadcast_sd(&B[(i*56)+47]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b47 = _mm_loaddup_pd(&B[(i*56)+47]);
#endif
__m128d c47_0 = _mm_load_sd(&C[(i*56)+11]);
__m128d a47_0 = _mm_load_sd(&values[159]);
#if defined(__SSE3__) && defined(__AVX__)
c47_0 = _mm_add_sd(c47_0, _mm_mul_sd(a47_0, _mm256_castpd256_pd128(b47)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c47_0 = _mm_add_sd(c47_0, _mm_mul_sd(a47_0, b47));
#endif
_mm_store_sd(&C[(i*56)+11], c47_0);
__m128d c47_1 = _mm_load_sd(&C[(i*56)+26]);
__m128d a47_1 = _mm_load_sd(&values[160]);
#if defined(__SSE3__) && defined(__AVX__)
c47_1 = _mm_add_sd(c47_1, _mm_mul_sd(a47_1, _mm256_castpd256_pd128(b47)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c47_1 = _mm_add_sd(c47_1, _mm_mul_sd(a47_1, b47));
#endif
_mm_store_sd(&C[(i*56)+26], c47_1);
__m128d c47_2 = _mm_load_sd(&C[(i*56)+47]);
__m128d a47_2 = _mm_load_sd(&values[161]);
#if defined(__SSE3__) && defined(__AVX__)
c47_2 = _mm_add_sd(c47_2, _mm_mul_sd(a47_2, _mm256_castpd256_pd128(b47)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c47_2 = _mm_add_sd(c47_2, _mm_mul_sd(a47_2, b47));
#endif
_mm_store_sd(&C[(i*56)+47], c47_2);
#else
C[(i*56)+11] += values[159] * B[(i*56)+47];
C[(i*56)+26] += values[160] * B[(i*56)+47];
C[(i*56)+47] += values[161] * B[(i*56)+47];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b48 = _mm256_broadcast_sd(&B[(i*56)+48]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b48 = _mm_loaddup_pd(&B[(i*56)+48]);
#endif
__m128d c48_0 = _mm_load_sd(&C[(i*56)+12]);
__m128d a48_0 = _mm_load_sd(&values[162]);
#if defined(__SSE3__) && defined(__AVX__)
c48_0 = _mm_add_sd(c48_0, _mm_mul_sd(a48_0, _mm256_castpd256_pd128(b48)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c48_0 = _mm_add_sd(c48_0, _mm_mul_sd(a48_0, b48));
#endif
_mm_store_sd(&C[(i*56)+12], c48_0);
__m128d c48_1 = _mm_load_sd(&C[(i*56)+27]);
__m128d a48_1 = _mm_load_sd(&values[163]);
#if defined(__SSE3__) && defined(__AVX__)
c48_1 = _mm_add_sd(c48_1, _mm_mul_sd(a48_1, _mm256_castpd256_pd128(b48)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c48_1 = _mm_add_sd(c48_1, _mm_mul_sd(a48_1, b48));
#endif
_mm_store_sd(&C[(i*56)+27], c48_1);
__m128d c48_2 = _mm_load_sd(&C[(i*56)+48]);
__m128d a48_2 = _mm_load_sd(&values[164]);
#if defined(__SSE3__) && defined(__AVX__)
c48_2 = _mm_add_sd(c48_2, _mm_mul_sd(a48_2, _mm256_castpd256_pd128(b48)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c48_2 = _mm_add_sd(c48_2, _mm_mul_sd(a48_2, b48));
#endif
_mm_store_sd(&C[(i*56)+48], c48_2);
#else
C[(i*56)+12] += values[162] * B[(i*56)+48];
C[(i*56)+27] += values[163] * B[(i*56)+48];
C[(i*56)+48] += values[164] * B[(i*56)+48];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b49 = _mm256_broadcast_sd(&B[(i*56)+49]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b49 = _mm_loaddup_pd(&B[(i*56)+49]);
#endif
__m128d c49_0 = _mm_load_sd(&C[(i*56)+13]);
__m128d a49_0 = _mm_load_sd(&values[165]);
#if defined(__SSE3__) && defined(__AVX__)
c49_0 = _mm_add_sd(c49_0, _mm_mul_sd(a49_0, _mm256_castpd256_pd128(b49)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c49_0 = _mm_add_sd(c49_0, _mm_mul_sd(a49_0, b49));
#endif
_mm_store_sd(&C[(i*56)+13], c49_0);
__m128d c49_1 = _mm_load_sd(&C[(i*56)+28]);
__m128d a49_1 = _mm_load_sd(&values[166]);
#if defined(__SSE3__) && defined(__AVX__)
c49_1 = _mm_add_sd(c49_1, _mm_mul_sd(a49_1, _mm256_castpd256_pd128(b49)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c49_1 = _mm_add_sd(c49_1, _mm_mul_sd(a49_1, b49));
#endif
_mm_store_sd(&C[(i*56)+28], c49_1);
__m128d c49_2 = _mm_load_sd(&C[(i*56)+49]);
__m128d a49_2 = _mm_load_sd(&values[167]);
#if defined(__SSE3__) && defined(__AVX__)
c49_2 = _mm_add_sd(c49_2, _mm_mul_sd(a49_2, _mm256_castpd256_pd128(b49)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c49_2 = _mm_add_sd(c49_2, _mm_mul_sd(a49_2, b49));
#endif
_mm_store_sd(&C[(i*56)+49], c49_2);
#else
C[(i*56)+13] += values[165] * B[(i*56)+49];
C[(i*56)+28] += values[166] * B[(i*56)+49];
C[(i*56)+49] += values[167] * B[(i*56)+49];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b50 = _mm256_broadcast_sd(&B[(i*56)+50]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b50 = _mm_loaddup_pd(&B[(i*56)+50]);
#endif
__m128d c50_0 = _mm_load_sd(&C[(i*56)+4]);
__m128d a50_0 = _mm_load_sd(&values[168]);
#if defined(__SSE3__) && defined(__AVX__)
c50_0 = _mm_add_sd(c50_0, _mm_mul_sd(a50_0, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c50_0 = _mm_add_sd(c50_0, _mm_mul_sd(a50_0, b50));
#endif
_mm_store_sd(&C[(i*56)+4], c50_0);
__m128d c50_1 = _mm_load_sd(&C[(i*56)+14]);
__m128d a50_1 = _mm_load_sd(&values[169]);
#if defined(__SSE3__) && defined(__AVX__)
c50_1 = _mm_add_sd(c50_1, _mm_mul_sd(a50_1, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c50_1 = _mm_add_sd(c50_1, _mm_mul_sd(a50_1, b50));
#endif
_mm_store_sd(&C[(i*56)+14], c50_1);
__m128d c50_2 = _mm_load_sd(&C[(i*56)+29]);
__m128d a50_2 = _mm_load_sd(&values[170]);
#if defined(__SSE3__) && defined(__AVX__)
c50_2 = _mm_add_sd(c50_2, _mm_mul_sd(a50_2, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c50_2 = _mm_add_sd(c50_2, _mm_mul_sd(a50_2, b50));
#endif
_mm_store_sd(&C[(i*56)+29], c50_2);
__m128d c50_3 = _mm_load_sd(&C[(i*56)+50]);
__m128d a50_3 = _mm_load_sd(&values[171]);
#if defined(__SSE3__) && defined(__AVX__)
c50_3 = _mm_add_sd(c50_3, _mm_mul_sd(a50_3, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c50_3 = _mm_add_sd(c50_3, _mm_mul_sd(a50_3, b50));
#endif
_mm_store_sd(&C[(i*56)+50], c50_3);
#else
C[(i*56)+4] += values[168] * B[(i*56)+50];
C[(i*56)+14] += values[169] * B[(i*56)+50];
C[(i*56)+29] += values[170] * B[(i*56)+50];
C[(i*56)+50] += values[171] * B[(i*56)+50];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b51 = _mm256_broadcast_sd(&B[(i*56)+51]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b51 = _mm_loaddup_pd(&B[(i*56)+51]);
#endif
__m128d c51_0 = _mm_load_sd(&C[(i*56)+5]);
__m128d a51_0 = _mm_load_sd(&values[172]);
#if defined(__SSE3__) && defined(__AVX__)
c51_0 = _mm_add_sd(c51_0, _mm_mul_sd(a51_0, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c51_0 = _mm_add_sd(c51_0, _mm_mul_sd(a51_0, b51));
#endif
_mm_store_sd(&C[(i*56)+5], c51_0);
__m128d c51_1 = _mm_load_sd(&C[(i*56)+15]);
__m128d a51_1 = _mm_load_sd(&values[173]);
#if defined(__SSE3__) && defined(__AVX__)
c51_1 = _mm_add_sd(c51_1, _mm_mul_sd(a51_1, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c51_1 = _mm_add_sd(c51_1, _mm_mul_sd(a51_1, b51));
#endif
_mm_store_sd(&C[(i*56)+15], c51_1);
__m128d c51_2 = _mm_load_sd(&C[(i*56)+30]);
__m128d a51_2 = _mm_load_sd(&values[174]);
#if defined(__SSE3__) && defined(__AVX__)
c51_2 = _mm_add_sd(c51_2, _mm_mul_sd(a51_2, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c51_2 = _mm_add_sd(c51_2, _mm_mul_sd(a51_2, b51));
#endif
_mm_store_sd(&C[(i*56)+30], c51_2);
__m128d c51_3 = _mm_load_sd(&C[(i*56)+51]);
__m128d a51_3 = _mm_load_sd(&values[175]);
#if defined(__SSE3__) && defined(__AVX__)
c51_3 = _mm_add_sd(c51_3, _mm_mul_sd(a51_3, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c51_3 = _mm_add_sd(c51_3, _mm_mul_sd(a51_3, b51));
#endif
_mm_store_sd(&C[(i*56)+51], c51_3);
#else
C[(i*56)+5] += values[172] * B[(i*56)+51];
C[(i*56)+15] += values[173] * B[(i*56)+51];
C[(i*56)+30] += values[174] * B[(i*56)+51];
C[(i*56)+51] += values[175] * B[(i*56)+51];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b52 = _mm256_broadcast_sd(&B[(i*56)+52]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b52 = _mm_loaddup_pd(&B[(i*56)+52]);
#endif
__m128d c52_0 = _mm_load_sd(&C[(i*56)+6]);
__m128d a52_0 = _mm_load_sd(&values[176]);
#if defined(__SSE3__) && defined(__AVX__)
c52_0 = _mm_add_sd(c52_0, _mm_mul_sd(a52_0, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c52_0 = _mm_add_sd(c52_0, _mm_mul_sd(a52_0, b52));
#endif
_mm_store_sd(&C[(i*56)+6], c52_0);
__m128d c52_1 = _mm_load_sd(&C[(i*56)+16]);
__m128d a52_1 = _mm_load_sd(&values[177]);
#if defined(__SSE3__) && defined(__AVX__)
c52_1 = _mm_add_sd(c52_1, _mm_mul_sd(a52_1, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c52_1 = _mm_add_sd(c52_1, _mm_mul_sd(a52_1, b52));
#endif
_mm_store_sd(&C[(i*56)+16], c52_1);
__m128d c52_2 = _mm_load_sd(&C[(i*56)+31]);
__m128d a52_2 = _mm_load_sd(&values[178]);
#if defined(__SSE3__) && defined(__AVX__)
c52_2 = _mm_add_sd(c52_2, _mm_mul_sd(a52_2, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c52_2 = _mm_add_sd(c52_2, _mm_mul_sd(a52_2, b52));
#endif
_mm_store_sd(&C[(i*56)+31], c52_2);
__m128d c52_3 = _mm_load_sd(&C[(i*56)+52]);
__m128d a52_3 = _mm_load_sd(&values[179]);
#if defined(__SSE3__) && defined(__AVX__)
c52_3 = _mm_add_sd(c52_3, _mm_mul_sd(a52_3, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c52_3 = _mm_add_sd(c52_3, _mm_mul_sd(a52_3, b52));
#endif
_mm_store_sd(&C[(i*56)+52], c52_3);
#else
C[(i*56)+6] += values[176] * B[(i*56)+52];
C[(i*56)+16] += values[177] * B[(i*56)+52];
C[(i*56)+31] += values[178] * B[(i*56)+52];
C[(i*56)+52] += values[179] * B[(i*56)+52];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b53 = _mm256_broadcast_sd(&B[(i*56)+53]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b53 = _mm_loaddup_pd(&B[(i*56)+53]);
#endif
__m128d c53_0 = _mm_load_sd(&C[(i*56)+1]);
__m128d a53_0 = _mm_load_sd(&values[180]);
#if defined(__SSE3__) && defined(__AVX__)
c53_0 = _mm_add_sd(c53_0, _mm_mul_sd(a53_0, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c53_0 = _mm_add_sd(c53_0, _mm_mul_sd(a53_0, b53));
#endif
_mm_store_sd(&C[(i*56)+1], c53_0);
__m128d c53_1 = _mm_load_sd(&C[(i*56)+7]);
__m128d a53_1 = _mm_load_sd(&values[181]);
#if defined(__SSE3__) && defined(__AVX__)
c53_1 = _mm_add_sd(c53_1, _mm_mul_sd(a53_1, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c53_1 = _mm_add_sd(c53_1, _mm_mul_sd(a53_1, b53));
#endif
_mm_store_sd(&C[(i*56)+7], c53_1);
__m128d c53_2 = _mm_load_sd(&C[(i*56)+17]);
__m128d a53_2 = _mm_load_sd(&values[182]);
#if defined(__SSE3__) && defined(__AVX__)
c53_2 = _mm_add_sd(c53_2, _mm_mul_sd(a53_2, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c53_2 = _mm_add_sd(c53_2, _mm_mul_sd(a53_2, b53));
#endif
_mm_store_sd(&C[(i*56)+17], c53_2);
__m128d c53_3 = _mm_load_sd(&C[(i*56)+32]);
__m128d a53_3 = _mm_load_sd(&values[183]);
#if defined(__SSE3__) && defined(__AVX__)
c53_3 = _mm_add_sd(c53_3, _mm_mul_sd(a53_3, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c53_3 = _mm_add_sd(c53_3, _mm_mul_sd(a53_3, b53));
#endif
_mm_store_sd(&C[(i*56)+32], c53_3);
__m128d c53_4 = _mm_load_sd(&C[(i*56)+53]);
__m128d a53_4 = _mm_load_sd(&values[184]);
#if defined(__SSE3__) && defined(__AVX__)
c53_4 = _mm_add_sd(c53_4, _mm_mul_sd(a53_4, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c53_4 = _mm_add_sd(c53_4, _mm_mul_sd(a53_4, b53));
#endif
_mm_store_sd(&C[(i*56)+53], c53_4);
#else
C[(i*56)+1] += values[180] * B[(i*56)+53];
C[(i*56)+7] += values[181] * B[(i*56)+53];
C[(i*56)+17] += values[182] * B[(i*56)+53];
C[(i*56)+32] += values[183] * B[(i*56)+53];
C[(i*56)+53] += values[184] * B[(i*56)+53];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b54 = _mm256_broadcast_sd(&B[(i*56)+54]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b54 = _mm_loaddup_pd(&B[(i*56)+54]);
#endif
__m128d c54_0 = _mm_load_sd(&C[(i*56)+2]);
__m128d a54_0 = _mm_load_sd(&values[185]);
#if defined(__SSE3__) && defined(__AVX__)
c54_0 = _mm_add_sd(c54_0, _mm_mul_sd(a54_0, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c54_0 = _mm_add_sd(c54_0, _mm_mul_sd(a54_0, b54));
#endif
_mm_store_sd(&C[(i*56)+2], c54_0);
__m128d c54_1 = _mm_load_sd(&C[(i*56)+8]);
__m128d a54_1 = _mm_load_sd(&values[186]);
#if defined(__SSE3__) && defined(__AVX__)
c54_1 = _mm_add_sd(c54_1, _mm_mul_sd(a54_1, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c54_1 = _mm_add_sd(c54_1, _mm_mul_sd(a54_1, b54));
#endif
_mm_store_sd(&C[(i*56)+8], c54_1);
__m128d c54_2 = _mm_load_sd(&C[(i*56)+18]);
__m128d a54_2 = _mm_load_sd(&values[187]);
#if defined(__SSE3__) && defined(__AVX__)
c54_2 = _mm_add_sd(c54_2, _mm_mul_sd(a54_2, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c54_2 = _mm_add_sd(c54_2, _mm_mul_sd(a54_2, b54));
#endif
_mm_store_sd(&C[(i*56)+18], c54_2);
__m128d c54_3 = _mm_load_sd(&C[(i*56)+33]);
__m128d a54_3 = _mm_load_sd(&values[188]);
#if defined(__SSE3__) && defined(__AVX__)
c54_3 = _mm_add_sd(c54_3, _mm_mul_sd(a54_3, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c54_3 = _mm_add_sd(c54_3, _mm_mul_sd(a54_3, b54));
#endif
_mm_store_sd(&C[(i*56)+33], c54_3);
__m128d c54_4 = _mm_load_sd(&C[(i*56)+54]);
__m128d a54_4 = _mm_load_sd(&values[189]);
#if defined(__SSE3__) && defined(__AVX__)
c54_4 = _mm_add_sd(c54_4, _mm_mul_sd(a54_4, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c54_4 = _mm_add_sd(c54_4, _mm_mul_sd(a54_4, b54));
#endif
_mm_store_sd(&C[(i*56)+54], c54_4);
#else
C[(i*56)+2] += values[185] * B[(i*56)+54];
C[(i*56)+8] += values[186] * B[(i*56)+54];
C[(i*56)+18] += values[187] * B[(i*56)+54];
C[(i*56)+33] += values[188] * B[(i*56)+54];
C[(i*56)+54] += values[189] * B[(i*56)+54];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b55 = _mm256_broadcast_sd(&B[(i*56)+55]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b55 = _mm_loaddup_pd(&B[(i*56)+55]);
#endif
__m128d c55_0 = _mm_load_sd(&C[(i*56)+0]);
__m128d a55_0 = _mm_load_sd(&values[190]);
#if defined(__SSE3__) && defined(__AVX__)
c55_0 = _mm_add_sd(c55_0, _mm_mul_sd(a55_0, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_0 = _mm_add_sd(c55_0, _mm_mul_sd(a55_0, b55));
#endif
_mm_store_sd(&C[(i*56)+0], c55_0);
__m128d c55_1 = _mm_load_sd(&C[(i*56)+3]);
__m128d a55_1 = _mm_load_sd(&values[191]);
#if defined(__SSE3__) && defined(__AVX__)
c55_1 = _mm_add_sd(c55_1, _mm_mul_sd(a55_1, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_1 = _mm_add_sd(c55_1, _mm_mul_sd(a55_1, b55));
#endif
_mm_store_sd(&C[(i*56)+3], c55_1);
__m128d c55_2 = _mm_load_sd(&C[(i*56)+9]);
__m128d a55_2 = _mm_load_sd(&values[192]);
#if defined(__SSE3__) && defined(__AVX__)
c55_2 = _mm_add_sd(c55_2, _mm_mul_sd(a55_2, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_2 = _mm_add_sd(c55_2, _mm_mul_sd(a55_2, b55));
#endif
_mm_store_sd(&C[(i*56)+9], c55_2);
__m128d c55_3 = _mm_load_sd(&C[(i*56)+19]);
__m128d a55_3 = _mm_load_sd(&values[193]);
#if defined(__SSE3__) && defined(__AVX__)
c55_3 = _mm_add_sd(c55_3, _mm_mul_sd(a55_3, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_3 = _mm_add_sd(c55_3, _mm_mul_sd(a55_3, b55));
#endif
_mm_store_sd(&C[(i*56)+19], c55_3);
__m128d c55_4 = _mm_load_sd(&C[(i*56)+34]);
__m128d a55_4 = _mm_load_sd(&values[194]);
#if defined(__SSE3__) && defined(__AVX__)
c55_4 = _mm_add_sd(c55_4, _mm_mul_sd(a55_4, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_4 = _mm_add_sd(c55_4, _mm_mul_sd(a55_4, b55));
#endif
_mm_store_sd(&C[(i*56)+34], c55_4);
__m128d c55_5 = _mm_load_sd(&C[(i*56)+55]);
__m128d a55_5 = _mm_load_sd(&values[195]);
#if defined(__SSE3__) && defined(__AVX__)
c55_5 = _mm_add_sd(c55_5, _mm_mul_sd(a55_5, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_5 = _mm_add_sd(c55_5, _mm_mul_sd(a55_5, b55));
#endif
_mm_store_sd(&C[(i*56)+55], c55_5);
#else
C[(i*56)+0] += values[190] * B[(i*56)+55];
C[(i*56)+3] += values[191] * B[(i*56)+55];
C[(i*56)+9] += values[192] * B[(i*56)+55];
C[(i*56)+19] += values[193] * B[(i*56)+55];
C[(i*56)+34] += values[194] * B[(i*56)+55];
C[(i*56)+55] += values[195] * B[(i*56)+55];
#endif

}

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 3528;
#endif

}

void dsparse_starMatrix_m84_n9_k9_ldA88_ldBna7_ldC88_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(32)
#pragma vector aligned
for (unsigned int i = 0; i < 84; i += 1)
{
C[(i)+(0)] += A[(i)+(528)] * values[0];
C[(i)+(0)] += A[(i)+(616)] * values[1];
C[(i)+(0)] += A[(i)+(704)] * values[2];
C[(i)+(88)] += A[(i)+(528)] * values[3];
C[(i)+(88)] += A[(i)+(616)] * values[4];
C[(i)+(88)] += A[(i)+(704)] * values[5];
C[(i)+(176)] += A[(i)+(528)] * values[6];
C[(i)+(176)] += A[(i)+(616)] * values[7];
C[(i)+(176)] += A[(i)+(704)] * values[8];
C[(i)+(264)] += A[(i)+(528)] * values[9];
C[(i)+(264)] += A[(i)+(616)] * values[10];
C[(i)+(352)] += A[(i)+(616)] * values[11];
C[(i)+(352)] += A[(i)+(704)] * values[12];
C[(i)+(440)] += A[(i)+(528)] * values[13];
C[(i)+(440)] += A[(i)+(704)] * values[14];
C[(i)+(528)] += A[(i)+(0)] * values[15];
C[(i)+(528)] += A[(i)+(264)] * values[16];
C[(i)+(528)] += A[(i)+(440)] * values[17];
C[(i)+(616)] += A[(i)+(88)] * values[18];
C[(i)+(616)] += A[(i)+(264)] * values[19];
C[(i)+(616)] += A[(i)+(352)] * values[20];
C[(i)+(704)] += A[(i)+(176)] * values[21];
C[(i)+(704)] += A[(i)+(352)] * values[22];
C[(i)+(704)] += A[(i)+(440)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 4032;
#endif

}

void dsparse_fP113DivM_m84_n9_k84_ldAna7_ldB88_ldC88_beta0_pfsigonly(const double* values, const double* B, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma nounroll_and_jam
for (unsigned int i = 0; i < 9; i++)
{
  #pragma simd
  for (unsigned int m = 0; m < 84; m++) {
    C[(i*88)+m] = 0.0;
  }
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b0 = _mm256_broadcast_sd(&B[(i*88)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b0 = _mm_loaddup_pd(&B[(i*88)+0]);
#endif
__m128d c0_0 = _mm_load_sd(&C[(i*88)+0]);
__m128d a0_0 = _mm_load_sd(&values[0]);
#if defined(__SSE3__) && defined(__AVX__)
c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, b0));
#endif
_mm_store_sd(&C[(i*88)+0], c0_0);
__m128d c0_1 = _mm_load_sd(&C[(i*88)+3]);
__m128d a0_1 = _mm_load_sd(&values[1]);
#if defined(__SSE3__) && defined(__AVX__)
c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, b0));
#endif
_mm_store_sd(&C[(i*88)+3], c0_1);
__m128d c0_2 = _mm_load_sd(&C[(i*88)+9]);
__m128d a0_2 = _mm_load_sd(&values[2]);
#if defined(__SSE3__) && defined(__AVX__)
c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, b0));
#endif
_mm_store_sd(&C[(i*88)+9], c0_2);
__m128d c0_3 = _mm_load_sd(&C[(i*88)+19]);
__m128d a0_3 = _mm_load_sd(&values[3]);
#if defined(__SSE3__) && defined(__AVX__)
c0_3 = _mm_add_sd(c0_3, _mm_mul_sd(a0_3, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_3 = _mm_add_sd(c0_3, _mm_mul_sd(a0_3, b0));
#endif
_mm_store_sd(&C[(i*88)+19], c0_3);
__m128d c0_4 = _mm_load_sd(&C[(i*88)+34]);
__m128d a0_4 = _mm_load_sd(&values[4]);
#if defined(__SSE3__) && defined(__AVX__)
c0_4 = _mm_add_sd(c0_4, _mm_mul_sd(a0_4, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_4 = _mm_add_sd(c0_4, _mm_mul_sd(a0_4, b0));
#endif
_mm_store_sd(&C[(i*88)+34], c0_4);
__m128d c0_5 = _mm_load_sd(&C[(i*88)+55]);
__m128d a0_5 = _mm_load_sd(&values[5]);
#if defined(__SSE3__) && defined(__AVX__)
c0_5 = _mm_add_sd(c0_5, _mm_mul_sd(a0_5, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_5 = _mm_add_sd(c0_5, _mm_mul_sd(a0_5, b0));
#endif
_mm_store_sd(&C[(i*88)+55], c0_5);
__m128d c0_6 = _mm_load_sd(&C[(i*88)+83]);
__m128d a0_6 = _mm_load_sd(&values[6]);
#if defined(__SSE3__) && defined(__AVX__)
c0_6 = _mm_add_sd(c0_6, _mm_mul_sd(a0_6, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_6 = _mm_add_sd(c0_6, _mm_mul_sd(a0_6, b0));
#endif
_mm_store_sd(&C[(i*88)+83], c0_6);
#else
C[(i*88)+0] += values[0] * B[(i*88)+0];
C[(i*88)+3] += values[1] * B[(i*88)+0];
C[(i*88)+9] += values[2] * B[(i*88)+0];
C[(i*88)+19] += values[3] * B[(i*88)+0];
C[(i*88)+34] += values[4] * B[(i*88)+0];
C[(i*88)+55] += values[5] * B[(i*88)+0];
C[(i*88)+83] += values[6] * B[(i*88)+0];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b1 = _mm256_broadcast_sd(&B[(i*88)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b1 = _mm_loaddup_pd(&B[(i*88)+1]);
#endif
__m128d c1_0 = _mm_load_sd(&C[(i*88)+1]);
__m128d a1_0 = _mm_load_sd(&values[7]);
#if defined(__SSE3__) && defined(__AVX__)
c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, b1));
#endif
_mm_store_sd(&C[(i*88)+1], c1_0);
__m128d c1_1 = _mm_load_sd(&C[(i*88)+7]);
__m128d a1_1 = _mm_load_sd(&values[8]);
#if defined(__SSE3__) && defined(__AVX__)
c1_1 = _mm_add_sd(c1_1, _mm_mul_sd(a1_1, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_1 = _mm_add_sd(c1_1, _mm_mul_sd(a1_1, b1));
#endif
_mm_store_sd(&C[(i*88)+7], c1_1);
__m128d c1_2 = _mm_load_sd(&C[(i*88)+17]);
__m128d a1_2 = _mm_load_sd(&values[9]);
#if defined(__SSE3__) && defined(__AVX__)
c1_2 = _mm_add_sd(c1_2, _mm_mul_sd(a1_2, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_2 = _mm_add_sd(c1_2, _mm_mul_sd(a1_2, b1));
#endif
_mm_store_sd(&C[(i*88)+17], c1_2);
__m128d c1_3 = _mm_load_sd(&C[(i*88)+32]);
__m128d a1_3 = _mm_load_sd(&values[10]);
#if defined(__SSE3__) && defined(__AVX__)
c1_3 = _mm_add_sd(c1_3, _mm_mul_sd(a1_3, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_3 = _mm_add_sd(c1_3, _mm_mul_sd(a1_3, b1));
#endif
_mm_store_sd(&C[(i*88)+32], c1_3);
__m128d c1_4 = _mm_load_sd(&C[(i*88)+53]);
__m128d a1_4 = _mm_load_sd(&values[11]);
#if defined(__SSE3__) && defined(__AVX__)
c1_4 = _mm_add_sd(c1_4, _mm_mul_sd(a1_4, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_4 = _mm_add_sd(c1_4, _mm_mul_sd(a1_4, b1));
#endif
_mm_store_sd(&C[(i*88)+53], c1_4);
__m128d c1_5 = _mm_load_sd(&C[(i*88)+81]);
__m128d a1_5 = _mm_load_sd(&values[12]);
#if defined(__SSE3__) && defined(__AVX__)
c1_5 = _mm_add_sd(c1_5, _mm_mul_sd(a1_5, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_5 = _mm_add_sd(c1_5, _mm_mul_sd(a1_5, b1));
#endif
_mm_store_sd(&C[(i*88)+81], c1_5);
#else
C[(i*88)+1] += values[7] * B[(i*88)+1];
C[(i*88)+7] += values[8] * B[(i*88)+1];
C[(i*88)+17] += values[9] * B[(i*88)+1];
C[(i*88)+32] += values[10] * B[(i*88)+1];
C[(i*88)+53] += values[11] * B[(i*88)+1];
C[(i*88)+81] += values[12] * B[(i*88)+1];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b2 = _mm256_broadcast_sd(&B[(i*88)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b2 = _mm_loaddup_pd(&B[(i*88)+2]);
#endif
__m128d c2_0 = _mm_load_sd(&C[(i*88)+2]);
__m128d a2_0 = _mm_load_sd(&values[13]);
#if defined(__SSE3__) && defined(__AVX__)
c2_0 = _mm_add_sd(c2_0, _mm_mul_sd(a2_0, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_0 = _mm_add_sd(c2_0, _mm_mul_sd(a2_0, b2));
#endif
_mm_store_sd(&C[(i*88)+2], c2_0);
__m128d c2_1 = _mm_load_sd(&C[(i*88)+8]);
__m128d a2_1 = _mm_load_sd(&values[14]);
#if defined(__SSE3__) && defined(__AVX__)
c2_1 = _mm_add_sd(c2_1, _mm_mul_sd(a2_1, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_1 = _mm_add_sd(c2_1, _mm_mul_sd(a2_1, b2));
#endif
_mm_store_sd(&C[(i*88)+8], c2_1);
__m128d c2_2 = _mm_load_sd(&C[(i*88)+18]);
__m128d a2_2 = _mm_load_sd(&values[15]);
#if defined(__SSE3__) && defined(__AVX__)
c2_2 = _mm_add_sd(c2_2, _mm_mul_sd(a2_2, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_2 = _mm_add_sd(c2_2, _mm_mul_sd(a2_2, b2));
#endif
_mm_store_sd(&C[(i*88)+18], c2_2);
__m128d c2_3 = _mm_load_sd(&C[(i*88)+33]);
__m128d a2_3 = _mm_load_sd(&values[16]);
#if defined(__SSE3__) && defined(__AVX__)
c2_3 = _mm_add_sd(c2_3, _mm_mul_sd(a2_3, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_3 = _mm_add_sd(c2_3, _mm_mul_sd(a2_3, b2));
#endif
_mm_store_sd(&C[(i*88)+33], c2_3);
__m128d c2_4 = _mm_load_sd(&C[(i*88)+54]);
__m128d a2_4 = _mm_load_sd(&values[17]);
#if defined(__SSE3__) && defined(__AVX__)
c2_4 = _mm_add_sd(c2_4, _mm_mul_sd(a2_4, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_4 = _mm_add_sd(c2_4, _mm_mul_sd(a2_4, b2));
#endif
_mm_store_sd(&C[(i*88)+54], c2_4);
__m128d c2_5 = _mm_load_sd(&C[(i*88)+82]);
__m128d a2_5 = _mm_load_sd(&values[18]);
#if defined(__SSE3__) && defined(__AVX__)
c2_5 = _mm_add_sd(c2_5, _mm_mul_sd(a2_5, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_5 = _mm_add_sd(c2_5, _mm_mul_sd(a2_5, b2));
#endif
_mm_store_sd(&C[(i*88)+82], c2_5);
#else
C[(i*88)+2] += values[13] * B[(i*88)+2];
C[(i*88)+8] += values[14] * B[(i*88)+2];
C[(i*88)+18] += values[15] * B[(i*88)+2];
C[(i*88)+33] += values[16] * B[(i*88)+2];
C[(i*88)+54] += values[17] * B[(i*88)+2];
C[(i*88)+82] += values[18] * B[(i*88)+2];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b3 = _mm256_broadcast_sd(&B[(i*88)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b3 = _mm_loaddup_pd(&B[(i*88)+3]);
#endif
__m128d c3_0 = _mm_load_sd(&C[(i*88)+0]);
__m128d a3_0 = _mm_load_sd(&values[19]);
#if defined(__SSE3__) && defined(__AVX__)
c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, b3));
#endif
_mm_store_sd(&C[(i*88)+0], c3_0);
__m128d c3_1 = _mm_load_sd(&C[(i*88)+3]);
__m128d a3_1 = _mm_load_sd(&values[20]);
#if defined(__SSE3__) && defined(__AVX__)
c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, b3));
#endif
_mm_store_sd(&C[(i*88)+3], c3_1);
__m128d c3_2 = _mm_load_sd(&C[(i*88)+9]);
__m128d a3_2 = _mm_load_sd(&values[21]);
#if defined(__SSE3__) && defined(__AVX__)
c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, b3));
#endif
_mm_store_sd(&C[(i*88)+9], c3_2);
__m128d c3_3 = _mm_load_sd(&C[(i*88)+19]);
__m128d a3_3 = _mm_load_sd(&values[22]);
#if defined(__SSE3__) && defined(__AVX__)
c3_3 = _mm_add_sd(c3_3, _mm_mul_sd(a3_3, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_3 = _mm_add_sd(c3_3, _mm_mul_sd(a3_3, b3));
#endif
_mm_store_sd(&C[(i*88)+19], c3_3);
__m128d c3_4 = _mm_load_sd(&C[(i*88)+34]);
__m128d a3_4 = _mm_load_sd(&values[23]);
#if defined(__SSE3__) && defined(__AVX__)
c3_4 = _mm_add_sd(c3_4, _mm_mul_sd(a3_4, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_4 = _mm_add_sd(c3_4, _mm_mul_sd(a3_4, b3));
#endif
_mm_store_sd(&C[(i*88)+34], c3_4);
__m128d c3_5 = _mm_load_sd(&C[(i*88)+55]);
__m128d a3_5 = _mm_load_sd(&values[24]);
#if defined(__SSE3__) && defined(__AVX__)
c3_5 = _mm_add_sd(c3_5, _mm_mul_sd(a3_5, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_5 = _mm_add_sd(c3_5, _mm_mul_sd(a3_5, b3));
#endif
_mm_store_sd(&C[(i*88)+55], c3_5);
__m128d c3_6 = _mm_load_sd(&C[(i*88)+83]);
__m128d a3_6 = _mm_load_sd(&values[25]);
#if defined(__SSE3__) && defined(__AVX__)
c3_6 = _mm_add_sd(c3_6, _mm_mul_sd(a3_6, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_6 = _mm_add_sd(c3_6, _mm_mul_sd(a3_6, b3));
#endif
_mm_store_sd(&C[(i*88)+83], c3_6);
#else
C[(i*88)+0] += values[19] * B[(i*88)+3];
C[(i*88)+3] += values[20] * B[(i*88)+3];
C[(i*88)+9] += values[21] * B[(i*88)+3];
C[(i*88)+19] += values[22] * B[(i*88)+3];
C[(i*88)+34] += values[23] * B[(i*88)+3];
C[(i*88)+55] += values[24] * B[(i*88)+3];
C[(i*88)+83] += values[25] * B[(i*88)+3];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b4 = _mm256_broadcast_sd(&B[(i*88)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b4 = _mm_loaddup_pd(&B[(i*88)+4]);
#endif
__m128d c4_0 = _mm_load_sd(&C[(i*88)+4]);
__m128d a4_0 = _mm_load_sd(&values[26]);
#if defined(__SSE3__) && defined(__AVX__)
c4_0 = _mm_add_sd(c4_0, _mm_mul_sd(a4_0, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_0 = _mm_add_sd(c4_0, _mm_mul_sd(a4_0, b4));
#endif
_mm_store_sd(&C[(i*88)+4], c4_0);
__m128d c4_1 = _mm_load_sd(&C[(i*88)+14]);
__m128d a4_1 = _mm_load_sd(&values[27]);
#if defined(__SSE3__) && defined(__AVX__)
c4_1 = _mm_add_sd(c4_1, _mm_mul_sd(a4_1, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_1 = _mm_add_sd(c4_1, _mm_mul_sd(a4_1, b4));
#endif
_mm_store_sd(&C[(i*88)+14], c4_1);
__m128d c4_2 = _mm_load_sd(&C[(i*88)+29]);
__m128d a4_2 = _mm_load_sd(&values[28]);
#if defined(__SSE3__) && defined(__AVX__)
c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, b4));
#endif
_mm_store_sd(&C[(i*88)+29], c4_2);
__m128d c4_3 = _mm_load_sd(&C[(i*88)+50]);
__m128d a4_3 = _mm_load_sd(&values[29]);
#if defined(__SSE3__) && defined(__AVX__)
c4_3 = _mm_add_sd(c4_3, _mm_mul_sd(a4_3, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_3 = _mm_add_sd(c4_3, _mm_mul_sd(a4_3, b4));
#endif
_mm_store_sd(&C[(i*88)+50], c4_3);
__m128d c4_4 = _mm_load_sd(&C[(i*88)+78]);
__m128d a4_4 = _mm_load_sd(&values[30]);
#if defined(__SSE3__) && defined(__AVX__)
c4_4 = _mm_add_sd(c4_4, _mm_mul_sd(a4_4, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_4 = _mm_add_sd(c4_4, _mm_mul_sd(a4_4, b4));
#endif
_mm_store_sd(&C[(i*88)+78], c4_4);
#else
C[(i*88)+4] += values[26] * B[(i*88)+4];
C[(i*88)+14] += values[27] * B[(i*88)+4];
C[(i*88)+29] += values[28] * B[(i*88)+4];
C[(i*88)+50] += values[29] * B[(i*88)+4];
C[(i*88)+78] += values[30] * B[(i*88)+4];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b5 = _mm256_broadcast_sd(&B[(i*88)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b5 = _mm_loaddup_pd(&B[(i*88)+5]);
#endif
__m128d c5_0 = _mm_load_sd(&C[(i*88)+5]);
__m128d a5_0 = _mm_load_sd(&values[31]);
#if defined(__SSE3__) && defined(__AVX__)
c5_0 = _mm_add_sd(c5_0, _mm_mul_sd(a5_0, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_0 = _mm_add_sd(c5_0, _mm_mul_sd(a5_0, b5));
#endif
_mm_store_sd(&C[(i*88)+5], c5_0);
__m128d c5_1 = _mm_load_sd(&C[(i*88)+15]);
__m128d a5_1 = _mm_load_sd(&values[32]);
#if defined(__SSE3__) && defined(__AVX__)
c5_1 = _mm_add_sd(c5_1, _mm_mul_sd(a5_1, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_1 = _mm_add_sd(c5_1, _mm_mul_sd(a5_1, b5));
#endif
_mm_store_sd(&C[(i*88)+15], c5_1);
__m128d c5_2 = _mm_load_sd(&C[(i*88)+30]);
__m128d a5_2 = _mm_load_sd(&values[33]);
#if defined(__SSE3__) && defined(__AVX__)
c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, b5));
#endif
_mm_store_sd(&C[(i*88)+30], c5_2);
__m128d c5_3 = _mm_load_sd(&C[(i*88)+51]);
__m128d a5_3 = _mm_load_sd(&values[34]);
#if defined(__SSE3__) && defined(__AVX__)
c5_3 = _mm_add_sd(c5_3, _mm_mul_sd(a5_3, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_3 = _mm_add_sd(c5_3, _mm_mul_sd(a5_3, b5));
#endif
_mm_store_sd(&C[(i*88)+51], c5_3);
__m128d c5_4 = _mm_load_sd(&C[(i*88)+79]);
__m128d a5_4 = _mm_load_sd(&values[35]);
#if defined(__SSE3__) && defined(__AVX__)
c5_4 = _mm_add_sd(c5_4, _mm_mul_sd(a5_4, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_4 = _mm_add_sd(c5_4, _mm_mul_sd(a5_4, b5));
#endif
_mm_store_sd(&C[(i*88)+79], c5_4);
#else
C[(i*88)+5] += values[31] * B[(i*88)+5];
C[(i*88)+15] += values[32] * B[(i*88)+5];
C[(i*88)+30] += values[33] * B[(i*88)+5];
C[(i*88)+51] += values[34] * B[(i*88)+5];
C[(i*88)+79] += values[35] * B[(i*88)+5];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b6 = _mm256_broadcast_sd(&B[(i*88)+6]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b6 = _mm_loaddup_pd(&B[(i*88)+6]);
#endif
__m128d c6_0 = _mm_load_sd(&C[(i*88)+6]);
__m128d a6_0 = _mm_load_sd(&values[36]);
#if defined(__SSE3__) && defined(__AVX__)
c6_0 = _mm_add_sd(c6_0, _mm_mul_sd(a6_0, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_0 = _mm_add_sd(c6_0, _mm_mul_sd(a6_0, b6));
#endif
_mm_store_sd(&C[(i*88)+6], c6_0);
__m128d c6_1 = _mm_load_sd(&C[(i*88)+16]);
__m128d a6_1 = _mm_load_sd(&values[37]);
#if defined(__SSE3__) && defined(__AVX__)
c6_1 = _mm_add_sd(c6_1, _mm_mul_sd(a6_1, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_1 = _mm_add_sd(c6_1, _mm_mul_sd(a6_1, b6));
#endif
_mm_store_sd(&C[(i*88)+16], c6_1);
__m128d c6_2 = _mm_load_sd(&C[(i*88)+31]);
__m128d a6_2 = _mm_load_sd(&values[38]);
#if defined(__SSE3__) && defined(__AVX__)
c6_2 = _mm_add_sd(c6_2, _mm_mul_sd(a6_2, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_2 = _mm_add_sd(c6_2, _mm_mul_sd(a6_2, b6));
#endif
_mm_store_sd(&C[(i*88)+31], c6_2);
__m128d c6_3 = _mm_load_sd(&C[(i*88)+52]);
__m128d a6_3 = _mm_load_sd(&values[39]);
#if defined(__SSE3__) && defined(__AVX__)
c6_3 = _mm_add_sd(c6_3, _mm_mul_sd(a6_3, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_3 = _mm_add_sd(c6_3, _mm_mul_sd(a6_3, b6));
#endif
_mm_store_sd(&C[(i*88)+52], c6_3);
__m128d c6_4 = _mm_load_sd(&C[(i*88)+80]);
__m128d a6_4 = _mm_load_sd(&values[40]);
#if defined(__SSE3__) && defined(__AVX__)
c6_4 = _mm_add_sd(c6_4, _mm_mul_sd(a6_4, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_4 = _mm_add_sd(c6_4, _mm_mul_sd(a6_4, b6));
#endif
_mm_store_sd(&C[(i*88)+80], c6_4);
#else
C[(i*88)+6] += values[36] * B[(i*88)+6];
C[(i*88)+16] += values[37] * B[(i*88)+6];
C[(i*88)+31] += values[38] * B[(i*88)+6];
C[(i*88)+52] += values[39] * B[(i*88)+6];
C[(i*88)+80] += values[40] * B[(i*88)+6];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b7 = _mm256_broadcast_sd(&B[(i*88)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b7 = _mm_loaddup_pd(&B[(i*88)+7]);
#endif
__m128d c7_0 = _mm_load_sd(&C[(i*88)+1]);
__m128d a7_0 = _mm_load_sd(&values[41]);
#if defined(__SSE3__) && defined(__AVX__)
c7_0 = _mm_add_sd(c7_0, _mm_mul_sd(a7_0, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_0 = _mm_add_sd(c7_0, _mm_mul_sd(a7_0, b7));
#endif
_mm_store_sd(&C[(i*88)+1], c7_0);
__m128d c7_1 = _mm_load_sd(&C[(i*88)+7]);
__m128d a7_1 = _mm_load_sd(&values[42]);
#if defined(__SSE3__) && defined(__AVX__)
c7_1 = _mm_add_sd(c7_1, _mm_mul_sd(a7_1, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_1 = _mm_add_sd(c7_1, _mm_mul_sd(a7_1, b7));
#endif
_mm_store_sd(&C[(i*88)+7], c7_1);
__m128d c7_2 = _mm_load_sd(&C[(i*88)+17]);
__m128d a7_2 = _mm_load_sd(&values[43]);
#if defined(__SSE3__) && defined(__AVX__)
c7_2 = _mm_add_sd(c7_2, _mm_mul_sd(a7_2, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_2 = _mm_add_sd(c7_2, _mm_mul_sd(a7_2, b7));
#endif
_mm_store_sd(&C[(i*88)+17], c7_2);
__m128d c7_3 = _mm_load_sd(&C[(i*88)+32]);
__m128d a7_3 = _mm_load_sd(&values[44]);
#if defined(__SSE3__) && defined(__AVX__)
c7_3 = _mm_add_sd(c7_3, _mm_mul_sd(a7_3, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_3 = _mm_add_sd(c7_3, _mm_mul_sd(a7_3, b7));
#endif
_mm_store_sd(&C[(i*88)+32], c7_3);
__m128d c7_4 = _mm_load_sd(&C[(i*88)+53]);
__m128d a7_4 = _mm_load_sd(&values[45]);
#if defined(__SSE3__) && defined(__AVX__)
c7_4 = _mm_add_sd(c7_4, _mm_mul_sd(a7_4, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_4 = _mm_add_sd(c7_4, _mm_mul_sd(a7_4, b7));
#endif
_mm_store_sd(&C[(i*88)+53], c7_4);
__m128d c7_5 = _mm_load_sd(&C[(i*88)+81]);
__m128d a7_5 = _mm_load_sd(&values[46]);
#if defined(__SSE3__) && defined(__AVX__)
c7_5 = _mm_add_sd(c7_5, _mm_mul_sd(a7_5, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_5 = _mm_add_sd(c7_5, _mm_mul_sd(a7_5, b7));
#endif
_mm_store_sd(&C[(i*88)+81], c7_5);
#else
C[(i*88)+1] += values[41] * B[(i*88)+7];
C[(i*88)+7] += values[42] * B[(i*88)+7];
C[(i*88)+17] += values[43] * B[(i*88)+7];
C[(i*88)+32] += values[44] * B[(i*88)+7];
C[(i*88)+53] += values[45] * B[(i*88)+7];
C[(i*88)+81] += values[46] * B[(i*88)+7];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b8 = _mm256_broadcast_sd(&B[(i*88)+8]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b8 = _mm_loaddup_pd(&B[(i*88)+8]);
#endif
__m128d c8_0 = _mm_load_sd(&C[(i*88)+2]);
__m128d a8_0 = _mm_load_sd(&values[47]);
#if defined(__SSE3__) && defined(__AVX__)
c8_0 = _mm_add_sd(c8_0, _mm_mul_sd(a8_0, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_0 = _mm_add_sd(c8_0, _mm_mul_sd(a8_0, b8));
#endif
_mm_store_sd(&C[(i*88)+2], c8_0);
__m128d c8_1 = _mm_load_sd(&C[(i*88)+8]);
__m128d a8_1 = _mm_load_sd(&values[48]);
#if defined(__SSE3__) && defined(__AVX__)
c8_1 = _mm_add_sd(c8_1, _mm_mul_sd(a8_1, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_1 = _mm_add_sd(c8_1, _mm_mul_sd(a8_1, b8));
#endif
_mm_store_sd(&C[(i*88)+8], c8_1);
__m128d c8_2 = _mm_load_sd(&C[(i*88)+18]);
__m128d a8_2 = _mm_load_sd(&values[49]);
#if defined(__SSE3__) && defined(__AVX__)
c8_2 = _mm_add_sd(c8_2, _mm_mul_sd(a8_2, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_2 = _mm_add_sd(c8_2, _mm_mul_sd(a8_2, b8));
#endif
_mm_store_sd(&C[(i*88)+18], c8_2);
__m128d c8_3 = _mm_load_sd(&C[(i*88)+33]);
__m128d a8_3 = _mm_load_sd(&values[50]);
#if defined(__SSE3__) && defined(__AVX__)
c8_3 = _mm_add_sd(c8_3, _mm_mul_sd(a8_3, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_3 = _mm_add_sd(c8_3, _mm_mul_sd(a8_3, b8));
#endif
_mm_store_sd(&C[(i*88)+33], c8_3);
__m128d c8_4 = _mm_load_sd(&C[(i*88)+54]);
__m128d a8_4 = _mm_load_sd(&values[51]);
#if defined(__SSE3__) && defined(__AVX__)
c8_4 = _mm_add_sd(c8_4, _mm_mul_sd(a8_4, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_4 = _mm_add_sd(c8_4, _mm_mul_sd(a8_4, b8));
#endif
_mm_store_sd(&C[(i*88)+54], c8_4);
__m128d c8_5 = _mm_load_sd(&C[(i*88)+82]);
__m128d a8_5 = _mm_load_sd(&values[52]);
#if defined(__SSE3__) && defined(__AVX__)
c8_5 = _mm_add_sd(c8_5, _mm_mul_sd(a8_5, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_5 = _mm_add_sd(c8_5, _mm_mul_sd(a8_5, b8));
#endif
_mm_store_sd(&C[(i*88)+82], c8_5);
#else
C[(i*88)+2] += values[47] * B[(i*88)+8];
C[(i*88)+8] += values[48] * B[(i*88)+8];
C[(i*88)+18] += values[49] * B[(i*88)+8];
C[(i*88)+33] += values[50] * B[(i*88)+8];
C[(i*88)+54] += values[51] * B[(i*88)+8];
C[(i*88)+82] += values[52] * B[(i*88)+8];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b9 = _mm256_broadcast_sd(&B[(i*88)+9]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b9 = _mm_loaddup_pd(&B[(i*88)+9]);
#endif
__m128d c9_0 = _mm_load_sd(&C[(i*88)+0]);
__m128d a9_0 = _mm_load_sd(&values[53]);
#if defined(__SSE3__) && defined(__AVX__)
c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, b9));
#endif
_mm_store_sd(&C[(i*88)+0], c9_0);
__m128d c9_1 = _mm_load_sd(&C[(i*88)+3]);
__m128d a9_1 = _mm_load_sd(&values[54]);
#if defined(__SSE3__) && defined(__AVX__)
c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, b9));
#endif
_mm_store_sd(&C[(i*88)+3], c9_1);
__m128d c9_2 = _mm_load_sd(&C[(i*88)+9]);
__m128d a9_2 = _mm_load_sd(&values[55]);
#if defined(__SSE3__) && defined(__AVX__)
c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, b9));
#endif
_mm_store_sd(&C[(i*88)+9], c9_2);
__m128d c9_3 = _mm_load_sd(&C[(i*88)+19]);
__m128d a9_3 = _mm_load_sd(&values[56]);
#if defined(__SSE3__) && defined(__AVX__)
c9_3 = _mm_add_sd(c9_3, _mm_mul_sd(a9_3, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_3 = _mm_add_sd(c9_3, _mm_mul_sd(a9_3, b9));
#endif
_mm_store_sd(&C[(i*88)+19], c9_3);
__m128d c9_4 = _mm_load_sd(&C[(i*88)+34]);
__m128d a9_4 = _mm_load_sd(&values[57]);
#if defined(__SSE3__) && defined(__AVX__)
c9_4 = _mm_add_sd(c9_4, _mm_mul_sd(a9_4, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_4 = _mm_add_sd(c9_4, _mm_mul_sd(a9_4, b9));
#endif
_mm_store_sd(&C[(i*88)+34], c9_4);
__m128d c9_5 = _mm_load_sd(&C[(i*88)+55]);
__m128d a9_5 = _mm_load_sd(&values[58]);
#if defined(__SSE3__) && defined(__AVX__)
c9_5 = _mm_add_sd(c9_5, _mm_mul_sd(a9_5, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_5 = _mm_add_sd(c9_5, _mm_mul_sd(a9_5, b9));
#endif
_mm_store_sd(&C[(i*88)+55], c9_5);
__m128d c9_6 = _mm_load_sd(&C[(i*88)+83]);
__m128d a9_6 = _mm_load_sd(&values[59]);
#if defined(__SSE3__) && defined(__AVX__)
c9_6 = _mm_add_sd(c9_6, _mm_mul_sd(a9_6, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_6 = _mm_add_sd(c9_6, _mm_mul_sd(a9_6, b9));
#endif
_mm_store_sd(&C[(i*88)+83], c9_6);
#else
C[(i*88)+0] += values[53] * B[(i*88)+9];
C[(i*88)+3] += values[54] * B[(i*88)+9];
C[(i*88)+9] += values[55] * B[(i*88)+9];
C[(i*88)+19] += values[56] * B[(i*88)+9];
C[(i*88)+34] += values[57] * B[(i*88)+9];
C[(i*88)+55] += values[58] * B[(i*88)+9];
C[(i*88)+83] += values[59] * B[(i*88)+9];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b10 = _mm256_broadcast_sd(&B[(i*88)+10]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b10 = _mm_loaddup_pd(&B[(i*88)+10]);
#endif
__m128d c10_0 = _mm_load_sd(&C[(i*88)+10]);
__m128d a10_0 = _mm_load_sd(&values[60]);
#if defined(__SSE3__) && defined(__AVX__)
c10_0 = _mm_add_sd(c10_0, _mm_mul_sd(a10_0, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c10_0 = _mm_add_sd(c10_0, _mm_mul_sd(a10_0, b10));
#endif
_mm_store_sd(&C[(i*88)+10], c10_0);
__m128d c10_1 = _mm_load_sd(&C[(i*88)+25]);
__m128d a10_1 = _mm_load_sd(&values[61]);
#if defined(__SSE3__) && defined(__AVX__)
c10_1 = _mm_add_sd(c10_1, _mm_mul_sd(a10_1, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c10_1 = _mm_add_sd(c10_1, _mm_mul_sd(a10_1, b10));
#endif
_mm_store_sd(&C[(i*88)+25], c10_1);
__m128d c10_2 = _mm_load_sd(&C[(i*88)+46]);
__m128d a10_2 = _mm_load_sd(&values[62]);
#if defined(__SSE3__) && defined(__AVX__)
c10_2 = _mm_add_sd(c10_2, _mm_mul_sd(a10_2, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c10_2 = _mm_add_sd(c10_2, _mm_mul_sd(a10_2, b10));
#endif
_mm_store_sd(&C[(i*88)+46], c10_2);
__m128d c10_3 = _mm_load_sd(&C[(i*88)+74]);
__m128d a10_3 = _mm_load_sd(&values[63]);
#if defined(__SSE3__) && defined(__AVX__)
c10_3 = _mm_add_sd(c10_3, _mm_mul_sd(a10_3, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c10_3 = _mm_add_sd(c10_3, _mm_mul_sd(a10_3, b10));
#endif
_mm_store_sd(&C[(i*88)+74], c10_3);
#else
C[(i*88)+10] += values[60] * B[(i*88)+10];
C[(i*88)+25] += values[61] * B[(i*88)+10];
C[(i*88)+46] += values[62] * B[(i*88)+10];
C[(i*88)+74] += values[63] * B[(i*88)+10];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b11 = _mm256_broadcast_sd(&B[(i*88)+11]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b11 = _mm_loaddup_pd(&B[(i*88)+11]);
#endif
__m128d c11_0 = _mm_load_sd(&C[(i*88)+11]);
__m128d a11_0 = _mm_load_sd(&values[64]);
#if defined(__SSE3__) && defined(__AVX__)
c11_0 = _mm_add_sd(c11_0, _mm_mul_sd(a11_0, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c11_0 = _mm_add_sd(c11_0, _mm_mul_sd(a11_0, b11));
#endif
_mm_store_sd(&C[(i*88)+11], c11_0);
__m128d c11_1 = _mm_load_sd(&C[(i*88)+26]);
__m128d a11_1 = _mm_load_sd(&values[65]);
#if defined(__SSE3__) && defined(__AVX__)
c11_1 = _mm_add_sd(c11_1, _mm_mul_sd(a11_1, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c11_1 = _mm_add_sd(c11_1, _mm_mul_sd(a11_1, b11));
#endif
_mm_store_sd(&C[(i*88)+26], c11_1);
__m128d c11_2 = _mm_load_sd(&C[(i*88)+47]);
__m128d a11_2 = _mm_load_sd(&values[66]);
#if defined(__SSE3__) && defined(__AVX__)
c11_2 = _mm_add_sd(c11_2, _mm_mul_sd(a11_2, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c11_2 = _mm_add_sd(c11_2, _mm_mul_sd(a11_2, b11));
#endif
_mm_store_sd(&C[(i*88)+47], c11_2);
__m128d c11_3 = _mm_load_sd(&C[(i*88)+75]);
__m128d a11_3 = _mm_load_sd(&values[67]);
#if defined(__SSE3__) && defined(__AVX__)
c11_3 = _mm_add_sd(c11_3, _mm_mul_sd(a11_3, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c11_3 = _mm_add_sd(c11_3, _mm_mul_sd(a11_3, b11));
#endif
_mm_store_sd(&C[(i*88)+75], c11_3);
#else
C[(i*88)+11] += values[64] * B[(i*88)+11];
C[(i*88)+26] += values[65] * B[(i*88)+11];
C[(i*88)+47] += values[66] * B[(i*88)+11];
C[(i*88)+75] += values[67] * B[(i*88)+11];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b12 = _mm256_broadcast_sd(&B[(i*88)+12]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b12 = _mm_loaddup_pd(&B[(i*88)+12]);
#endif
__m128d c12_0 = _mm_load_sd(&C[(i*88)+12]);
__m128d a12_0 = _mm_load_sd(&values[68]);
#if defined(__SSE3__) && defined(__AVX__)
c12_0 = _mm_add_sd(c12_0, _mm_mul_sd(a12_0, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c12_0 = _mm_add_sd(c12_0, _mm_mul_sd(a12_0, b12));
#endif
_mm_store_sd(&C[(i*88)+12], c12_0);
__m128d c12_1 = _mm_load_sd(&C[(i*88)+27]);
__m128d a12_1 = _mm_load_sd(&values[69]);
#if defined(__SSE3__) && defined(__AVX__)
c12_1 = _mm_add_sd(c12_1, _mm_mul_sd(a12_1, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c12_1 = _mm_add_sd(c12_1, _mm_mul_sd(a12_1, b12));
#endif
_mm_store_sd(&C[(i*88)+27], c12_1);
__m128d c12_2 = _mm_load_sd(&C[(i*88)+48]);
__m128d a12_2 = _mm_load_sd(&values[70]);
#if defined(__SSE3__) && defined(__AVX__)
c12_2 = _mm_add_sd(c12_2, _mm_mul_sd(a12_2, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c12_2 = _mm_add_sd(c12_2, _mm_mul_sd(a12_2, b12));
#endif
_mm_store_sd(&C[(i*88)+48], c12_2);
__m128d c12_3 = _mm_load_sd(&C[(i*88)+76]);
__m128d a12_3 = _mm_load_sd(&values[71]);
#if defined(__SSE3__) && defined(__AVX__)
c12_3 = _mm_add_sd(c12_3, _mm_mul_sd(a12_3, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c12_3 = _mm_add_sd(c12_3, _mm_mul_sd(a12_3, b12));
#endif
_mm_store_sd(&C[(i*88)+76], c12_3);
#else
C[(i*88)+12] += values[68] * B[(i*88)+12];
C[(i*88)+27] += values[69] * B[(i*88)+12];
C[(i*88)+48] += values[70] * B[(i*88)+12];
C[(i*88)+76] += values[71] * B[(i*88)+12];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b13 = _mm256_broadcast_sd(&B[(i*88)+13]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b13 = _mm_loaddup_pd(&B[(i*88)+13]);
#endif
__m128d c13_0 = _mm_load_sd(&C[(i*88)+13]);
__m128d a13_0 = _mm_load_sd(&values[72]);
#if defined(__SSE3__) && defined(__AVX__)
c13_0 = _mm_add_sd(c13_0, _mm_mul_sd(a13_0, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c13_0 = _mm_add_sd(c13_0, _mm_mul_sd(a13_0, b13));
#endif
_mm_store_sd(&C[(i*88)+13], c13_0);
__m128d c13_1 = _mm_load_sd(&C[(i*88)+28]);
__m128d a13_1 = _mm_load_sd(&values[73]);
#if defined(__SSE3__) && defined(__AVX__)
c13_1 = _mm_add_sd(c13_1, _mm_mul_sd(a13_1, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c13_1 = _mm_add_sd(c13_1, _mm_mul_sd(a13_1, b13));
#endif
_mm_store_sd(&C[(i*88)+28], c13_1);
__m128d c13_2 = _mm_load_sd(&C[(i*88)+49]);
__m128d a13_2 = _mm_load_sd(&values[74]);
#if defined(__SSE3__) && defined(__AVX__)
c13_2 = _mm_add_sd(c13_2, _mm_mul_sd(a13_2, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c13_2 = _mm_add_sd(c13_2, _mm_mul_sd(a13_2, b13));
#endif
_mm_store_sd(&C[(i*88)+49], c13_2);
__m128d c13_3 = _mm_load_sd(&C[(i*88)+77]);
__m128d a13_3 = _mm_load_sd(&values[75]);
#if defined(__SSE3__) && defined(__AVX__)
c13_3 = _mm_add_sd(c13_3, _mm_mul_sd(a13_3, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c13_3 = _mm_add_sd(c13_3, _mm_mul_sd(a13_3, b13));
#endif
_mm_store_sd(&C[(i*88)+77], c13_3);
#else
C[(i*88)+13] += values[72] * B[(i*88)+13];
C[(i*88)+28] += values[73] * B[(i*88)+13];
C[(i*88)+49] += values[74] * B[(i*88)+13];
C[(i*88)+77] += values[75] * B[(i*88)+13];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b14 = _mm256_broadcast_sd(&B[(i*88)+14]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b14 = _mm_loaddup_pd(&B[(i*88)+14]);
#endif
__m128d c14_0 = _mm_load_sd(&C[(i*88)+4]);
__m128d a14_0 = _mm_load_sd(&values[76]);
#if defined(__SSE3__) && defined(__AVX__)
c14_0 = _mm_add_sd(c14_0, _mm_mul_sd(a14_0, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_0 = _mm_add_sd(c14_0, _mm_mul_sd(a14_0, b14));
#endif
_mm_store_sd(&C[(i*88)+4], c14_0);
__m128d c14_1 = _mm_load_sd(&C[(i*88)+14]);
__m128d a14_1 = _mm_load_sd(&values[77]);
#if defined(__SSE3__) && defined(__AVX__)
c14_1 = _mm_add_sd(c14_1, _mm_mul_sd(a14_1, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_1 = _mm_add_sd(c14_1, _mm_mul_sd(a14_1, b14));
#endif
_mm_store_sd(&C[(i*88)+14], c14_1);
__m128d c14_2 = _mm_load_sd(&C[(i*88)+29]);
__m128d a14_2 = _mm_load_sd(&values[78]);
#if defined(__SSE3__) && defined(__AVX__)
c14_2 = _mm_add_sd(c14_2, _mm_mul_sd(a14_2, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_2 = _mm_add_sd(c14_2, _mm_mul_sd(a14_2, b14));
#endif
_mm_store_sd(&C[(i*88)+29], c14_2);
__m128d c14_3 = _mm_load_sd(&C[(i*88)+50]);
__m128d a14_3 = _mm_load_sd(&values[79]);
#if defined(__SSE3__) && defined(__AVX__)
c14_3 = _mm_add_sd(c14_3, _mm_mul_sd(a14_3, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_3 = _mm_add_sd(c14_3, _mm_mul_sd(a14_3, b14));
#endif
_mm_store_sd(&C[(i*88)+50], c14_3);
__m128d c14_4 = _mm_load_sd(&C[(i*88)+78]);
__m128d a14_4 = _mm_load_sd(&values[80]);
#if defined(__SSE3__) && defined(__AVX__)
c14_4 = _mm_add_sd(c14_4, _mm_mul_sd(a14_4, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_4 = _mm_add_sd(c14_4, _mm_mul_sd(a14_4, b14));
#endif
_mm_store_sd(&C[(i*88)+78], c14_4);
#else
C[(i*88)+4] += values[76] * B[(i*88)+14];
C[(i*88)+14] += values[77] * B[(i*88)+14];
C[(i*88)+29] += values[78] * B[(i*88)+14];
C[(i*88)+50] += values[79] * B[(i*88)+14];
C[(i*88)+78] += values[80] * B[(i*88)+14];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b15 = _mm256_broadcast_sd(&B[(i*88)+15]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b15 = _mm_loaddup_pd(&B[(i*88)+15]);
#endif
__m128d c15_0 = _mm_load_sd(&C[(i*88)+5]);
__m128d a15_0 = _mm_load_sd(&values[81]);
#if defined(__SSE3__) && defined(__AVX__)
c15_0 = _mm_add_sd(c15_0, _mm_mul_sd(a15_0, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_0 = _mm_add_sd(c15_0, _mm_mul_sd(a15_0, b15));
#endif
_mm_store_sd(&C[(i*88)+5], c15_0);
__m128d c15_1 = _mm_load_sd(&C[(i*88)+15]);
__m128d a15_1 = _mm_load_sd(&values[82]);
#if defined(__SSE3__) && defined(__AVX__)
c15_1 = _mm_add_sd(c15_1, _mm_mul_sd(a15_1, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_1 = _mm_add_sd(c15_1, _mm_mul_sd(a15_1, b15));
#endif
_mm_store_sd(&C[(i*88)+15], c15_1);
__m128d c15_2 = _mm_load_sd(&C[(i*88)+30]);
__m128d a15_2 = _mm_load_sd(&values[83]);
#if defined(__SSE3__) && defined(__AVX__)
c15_2 = _mm_add_sd(c15_2, _mm_mul_sd(a15_2, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_2 = _mm_add_sd(c15_2, _mm_mul_sd(a15_2, b15));
#endif
_mm_store_sd(&C[(i*88)+30], c15_2);
__m128d c15_3 = _mm_load_sd(&C[(i*88)+51]);
__m128d a15_3 = _mm_load_sd(&values[84]);
#if defined(__SSE3__) && defined(__AVX__)
c15_3 = _mm_add_sd(c15_3, _mm_mul_sd(a15_3, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_3 = _mm_add_sd(c15_3, _mm_mul_sd(a15_3, b15));
#endif
_mm_store_sd(&C[(i*88)+51], c15_3);
__m128d c15_4 = _mm_load_sd(&C[(i*88)+79]);
__m128d a15_4 = _mm_load_sd(&values[85]);
#if defined(__SSE3__) && defined(__AVX__)
c15_4 = _mm_add_sd(c15_4, _mm_mul_sd(a15_4, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_4 = _mm_add_sd(c15_4, _mm_mul_sd(a15_4, b15));
#endif
_mm_store_sd(&C[(i*88)+79], c15_4);
#else
C[(i*88)+5] += values[81] * B[(i*88)+15];
C[(i*88)+15] += values[82] * B[(i*88)+15];
C[(i*88)+30] += values[83] * B[(i*88)+15];
C[(i*88)+51] += values[84] * B[(i*88)+15];
C[(i*88)+79] += values[85] * B[(i*88)+15];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b16 = _mm256_broadcast_sd(&B[(i*88)+16]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b16 = _mm_loaddup_pd(&B[(i*88)+16]);
#endif
__m128d c16_0 = _mm_load_sd(&C[(i*88)+6]);
__m128d a16_0 = _mm_load_sd(&values[86]);
#if defined(__SSE3__) && defined(__AVX__)
c16_0 = _mm_add_sd(c16_0, _mm_mul_sd(a16_0, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_0 = _mm_add_sd(c16_0, _mm_mul_sd(a16_0, b16));
#endif
_mm_store_sd(&C[(i*88)+6], c16_0);
__m128d c16_1 = _mm_load_sd(&C[(i*88)+16]);
__m128d a16_1 = _mm_load_sd(&values[87]);
#if defined(__SSE3__) && defined(__AVX__)
c16_1 = _mm_add_sd(c16_1, _mm_mul_sd(a16_1, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_1 = _mm_add_sd(c16_1, _mm_mul_sd(a16_1, b16));
#endif
_mm_store_sd(&C[(i*88)+16], c16_1);
__m128d c16_2 = _mm_load_sd(&C[(i*88)+31]);
__m128d a16_2 = _mm_load_sd(&values[88]);
#if defined(__SSE3__) && defined(__AVX__)
c16_2 = _mm_add_sd(c16_2, _mm_mul_sd(a16_2, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_2 = _mm_add_sd(c16_2, _mm_mul_sd(a16_2, b16));
#endif
_mm_store_sd(&C[(i*88)+31], c16_2);
__m128d c16_3 = _mm_load_sd(&C[(i*88)+52]);
__m128d a16_3 = _mm_load_sd(&values[89]);
#if defined(__SSE3__) && defined(__AVX__)
c16_3 = _mm_add_sd(c16_3, _mm_mul_sd(a16_3, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_3 = _mm_add_sd(c16_3, _mm_mul_sd(a16_3, b16));
#endif
_mm_store_sd(&C[(i*88)+52], c16_3);
__m128d c16_4 = _mm_load_sd(&C[(i*88)+80]);
__m128d a16_4 = _mm_load_sd(&values[90]);
#if defined(__SSE3__) && defined(__AVX__)
c16_4 = _mm_add_sd(c16_4, _mm_mul_sd(a16_4, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_4 = _mm_add_sd(c16_4, _mm_mul_sd(a16_4, b16));
#endif
_mm_store_sd(&C[(i*88)+80], c16_4);
#else
C[(i*88)+6] += values[86] * B[(i*88)+16];
C[(i*88)+16] += values[87] * B[(i*88)+16];
C[(i*88)+31] += values[88] * B[(i*88)+16];
C[(i*88)+52] += values[89] * B[(i*88)+16];
C[(i*88)+80] += values[90] * B[(i*88)+16];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b17 = _mm256_broadcast_sd(&B[(i*88)+17]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b17 = _mm_loaddup_pd(&B[(i*88)+17]);
#endif
__m128d c17_0 = _mm_load_sd(&C[(i*88)+1]);
__m128d a17_0 = _mm_load_sd(&values[91]);
#if defined(__SSE3__) && defined(__AVX__)
c17_0 = _mm_add_sd(c17_0, _mm_mul_sd(a17_0, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_0 = _mm_add_sd(c17_0, _mm_mul_sd(a17_0, b17));
#endif
_mm_store_sd(&C[(i*88)+1], c17_0);
__m128d c17_1 = _mm_load_sd(&C[(i*88)+7]);
__m128d a17_1 = _mm_load_sd(&values[92]);
#if defined(__SSE3__) && defined(__AVX__)
c17_1 = _mm_add_sd(c17_1, _mm_mul_sd(a17_1, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_1 = _mm_add_sd(c17_1, _mm_mul_sd(a17_1, b17));
#endif
_mm_store_sd(&C[(i*88)+7], c17_1);
__m128d c17_2 = _mm_load_sd(&C[(i*88)+17]);
__m128d a17_2 = _mm_load_sd(&values[93]);
#if defined(__SSE3__) && defined(__AVX__)
c17_2 = _mm_add_sd(c17_2, _mm_mul_sd(a17_2, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_2 = _mm_add_sd(c17_2, _mm_mul_sd(a17_2, b17));
#endif
_mm_store_sd(&C[(i*88)+17], c17_2);
__m128d c17_3 = _mm_load_sd(&C[(i*88)+32]);
__m128d a17_3 = _mm_load_sd(&values[94]);
#if defined(__SSE3__) && defined(__AVX__)
c17_3 = _mm_add_sd(c17_3, _mm_mul_sd(a17_3, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_3 = _mm_add_sd(c17_3, _mm_mul_sd(a17_3, b17));
#endif
_mm_store_sd(&C[(i*88)+32], c17_3);
__m128d c17_4 = _mm_load_sd(&C[(i*88)+53]);
__m128d a17_4 = _mm_load_sd(&values[95]);
#if defined(__SSE3__) && defined(__AVX__)
c17_4 = _mm_add_sd(c17_4, _mm_mul_sd(a17_4, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_4 = _mm_add_sd(c17_4, _mm_mul_sd(a17_4, b17));
#endif
_mm_store_sd(&C[(i*88)+53], c17_4);
__m128d c17_5 = _mm_load_sd(&C[(i*88)+81]);
__m128d a17_5 = _mm_load_sd(&values[96]);
#if defined(__SSE3__) && defined(__AVX__)
c17_5 = _mm_add_sd(c17_5, _mm_mul_sd(a17_5, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_5 = _mm_add_sd(c17_5, _mm_mul_sd(a17_5, b17));
#endif
_mm_store_sd(&C[(i*88)+81], c17_5);
#else
C[(i*88)+1] += values[91] * B[(i*88)+17];
C[(i*88)+7] += values[92] * B[(i*88)+17];
C[(i*88)+17] += values[93] * B[(i*88)+17];
C[(i*88)+32] += values[94] * B[(i*88)+17];
C[(i*88)+53] += values[95] * B[(i*88)+17];
C[(i*88)+81] += values[96] * B[(i*88)+17];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b18 = _mm256_broadcast_sd(&B[(i*88)+18]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b18 = _mm_loaddup_pd(&B[(i*88)+18]);
#endif
__m128d c18_0 = _mm_load_sd(&C[(i*88)+2]);
__m128d a18_0 = _mm_load_sd(&values[97]);
#if defined(__SSE3__) && defined(__AVX__)
c18_0 = _mm_add_sd(c18_0, _mm_mul_sd(a18_0, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_0 = _mm_add_sd(c18_0, _mm_mul_sd(a18_0, b18));
#endif
_mm_store_sd(&C[(i*88)+2], c18_0);
__m128d c18_1 = _mm_load_sd(&C[(i*88)+8]);
__m128d a18_1 = _mm_load_sd(&values[98]);
#if defined(__SSE3__) && defined(__AVX__)
c18_1 = _mm_add_sd(c18_1, _mm_mul_sd(a18_1, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_1 = _mm_add_sd(c18_1, _mm_mul_sd(a18_1, b18));
#endif
_mm_store_sd(&C[(i*88)+8], c18_1);
__m128d c18_2 = _mm_load_sd(&C[(i*88)+18]);
__m128d a18_2 = _mm_load_sd(&values[99]);
#if defined(__SSE3__) && defined(__AVX__)
c18_2 = _mm_add_sd(c18_2, _mm_mul_sd(a18_2, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_2 = _mm_add_sd(c18_2, _mm_mul_sd(a18_2, b18));
#endif
_mm_store_sd(&C[(i*88)+18], c18_2);
__m128d c18_3 = _mm_load_sd(&C[(i*88)+33]);
__m128d a18_3 = _mm_load_sd(&values[100]);
#if defined(__SSE3__) && defined(__AVX__)
c18_3 = _mm_add_sd(c18_3, _mm_mul_sd(a18_3, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_3 = _mm_add_sd(c18_3, _mm_mul_sd(a18_3, b18));
#endif
_mm_store_sd(&C[(i*88)+33], c18_3);
__m128d c18_4 = _mm_load_sd(&C[(i*88)+54]);
__m128d a18_4 = _mm_load_sd(&values[101]);
#if defined(__SSE3__) && defined(__AVX__)
c18_4 = _mm_add_sd(c18_4, _mm_mul_sd(a18_4, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_4 = _mm_add_sd(c18_4, _mm_mul_sd(a18_4, b18));
#endif
_mm_store_sd(&C[(i*88)+54], c18_4);
__m128d c18_5 = _mm_load_sd(&C[(i*88)+82]);
__m128d a18_5 = _mm_load_sd(&values[102]);
#if defined(__SSE3__) && defined(__AVX__)
c18_5 = _mm_add_sd(c18_5, _mm_mul_sd(a18_5, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_5 = _mm_add_sd(c18_5, _mm_mul_sd(a18_5, b18));
#endif
_mm_store_sd(&C[(i*88)+82], c18_5);
#else
C[(i*88)+2] += values[97] * B[(i*88)+18];
C[(i*88)+8] += values[98] * B[(i*88)+18];
C[(i*88)+18] += values[99] * B[(i*88)+18];
C[(i*88)+33] += values[100] * B[(i*88)+18];
C[(i*88)+54] += values[101] * B[(i*88)+18];
C[(i*88)+82] += values[102] * B[(i*88)+18];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b19 = _mm256_broadcast_sd(&B[(i*88)+19]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b19 = _mm_loaddup_pd(&B[(i*88)+19]);
#endif
__m128d c19_0 = _mm_load_sd(&C[(i*88)+0]);
__m128d a19_0 = _mm_load_sd(&values[103]);
#if defined(__SSE3__) && defined(__AVX__)
c19_0 = _mm_add_sd(c19_0, _mm_mul_sd(a19_0, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_0 = _mm_add_sd(c19_0, _mm_mul_sd(a19_0, b19));
#endif
_mm_store_sd(&C[(i*88)+0], c19_0);
__m128d c19_1 = _mm_load_sd(&C[(i*88)+3]);
__m128d a19_1 = _mm_load_sd(&values[104]);
#if defined(__SSE3__) && defined(__AVX__)
c19_1 = _mm_add_sd(c19_1, _mm_mul_sd(a19_1, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_1 = _mm_add_sd(c19_1, _mm_mul_sd(a19_1, b19));
#endif
_mm_store_sd(&C[(i*88)+3], c19_1);
__m128d c19_2 = _mm_load_sd(&C[(i*88)+9]);
__m128d a19_2 = _mm_load_sd(&values[105]);
#if defined(__SSE3__) && defined(__AVX__)
c19_2 = _mm_add_sd(c19_2, _mm_mul_sd(a19_2, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_2 = _mm_add_sd(c19_2, _mm_mul_sd(a19_2, b19));
#endif
_mm_store_sd(&C[(i*88)+9], c19_2);
__m128d c19_3 = _mm_load_sd(&C[(i*88)+19]);
__m128d a19_3 = _mm_load_sd(&values[106]);
#if defined(__SSE3__) && defined(__AVX__)
c19_3 = _mm_add_sd(c19_3, _mm_mul_sd(a19_3, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_3 = _mm_add_sd(c19_3, _mm_mul_sd(a19_3, b19));
#endif
_mm_store_sd(&C[(i*88)+19], c19_3);
__m128d c19_4 = _mm_load_sd(&C[(i*88)+34]);
__m128d a19_4 = _mm_load_sd(&values[107]);
#if defined(__SSE3__) && defined(__AVX__)
c19_4 = _mm_add_sd(c19_4, _mm_mul_sd(a19_4, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_4 = _mm_add_sd(c19_4, _mm_mul_sd(a19_4, b19));
#endif
_mm_store_sd(&C[(i*88)+34], c19_4);
__m128d c19_5 = _mm_load_sd(&C[(i*88)+55]);
__m128d a19_5 = _mm_load_sd(&values[108]);
#if defined(__SSE3__) && defined(__AVX__)
c19_5 = _mm_add_sd(c19_5, _mm_mul_sd(a19_5, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_5 = _mm_add_sd(c19_5, _mm_mul_sd(a19_5, b19));
#endif
_mm_store_sd(&C[(i*88)+55], c19_5);
__m128d c19_6 = _mm_load_sd(&C[(i*88)+83]);
__m128d a19_6 = _mm_load_sd(&values[109]);
#if defined(__SSE3__) && defined(__AVX__)
c19_6 = _mm_add_sd(c19_6, _mm_mul_sd(a19_6, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_6 = _mm_add_sd(c19_6, _mm_mul_sd(a19_6, b19));
#endif
_mm_store_sd(&C[(i*88)+83], c19_6);
#else
C[(i*88)+0] += values[103] * B[(i*88)+19];
C[(i*88)+3] += values[104] * B[(i*88)+19];
C[(i*88)+9] += values[105] * B[(i*88)+19];
C[(i*88)+19] += values[106] * B[(i*88)+19];
C[(i*88)+34] += values[107] * B[(i*88)+19];
C[(i*88)+55] += values[108] * B[(i*88)+19];
C[(i*88)+83] += values[109] * B[(i*88)+19];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b20 = _mm256_broadcast_sd(&B[(i*88)+20]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b20 = _mm_loaddup_pd(&B[(i*88)+20]);
#endif
__m128d c20_0 = _mm_load_sd(&C[(i*88)+20]);
__m128d a20_0 = _mm_load_sd(&values[110]);
#if defined(__SSE3__) && defined(__AVX__)
c20_0 = _mm_add_sd(c20_0, _mm_mul_sd(a20_0, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c20_0 = _mm_add_sd(c20_0, _mm_mul_sd(a20_0, b20));
#endif
_mm_store_sd(&C[(i*88)+20], c20_0);
__m128d c20_1 = _mm_load_sd(&C[(i*88)+41]);
__m128d a20_1 = _mm_load_sd(&values[111]);
#if defined(__SSE3__) && defined(__AVX__)
c20_1 = _mm_add_sd(c20_1, _mm_mul_sd(a20_1, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c20_1 = _mm_add_sd(c20_1, _mm_mul_sd(a20_1, b20));
#endif
_mm_store_sd(&C[(i*88)+41], c20_1);
__m128d c20_2 = _mm_load_sd(&C[(i*88)+69]);
__m128d a20_2 = _mm_load_sd(&values[112]);
#if defined(__SSE3__) && defined(__AVX__)
c20_2 = _mm_add_sd(c20_2, _mm_mul_sd(a20_2, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c20_2 = _mm_add_sd(c20_2, _mm_mul_sd(a20_2, b20));
#endif
_mm_store_sd(&C[(i*88)+69], c20_2);
#else
C[(i*88)+20] += values[110] * B[(i*88)+20];
C[(i*88)+41] += values[111] * B[(i*88)+20];
C[(i*88)+69] += values[112] * B[(i*88)+20];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b21 = _mm256_broadcast_sd(&B[(i*88)+21]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b21 = _mm_loaddup_pd(&B[(i*88)+21]);
#endif
__m128d c21_0 = _mm_load_sd(&C[(i*88)+21]);
__m128d a21_0 = _mm_load_sd(&values[113]);
#if defined(__SSE3__) && defined(__AVX__)
c21_0 = _mm_add_sd(c21_0, _mm_mul_sd(a21_0, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c21_0 = _mm_add_sd(c21_0, _mm_mul_sd(a21_0, b21));
#endif
_mm_store_sd(&C[(i*88)+21], c21_0);
__m128d c21_1 = _mm_load_sd(&C[(i*88)+42]);
__m128d a21_1 = _mm_load_sd(&values[114]);
#if defined(__SSE3__) && defined(__AVX__)
c21_1 = _mm_add_sd(c21_1, _mm_mul_sd(a21_1, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c21_1 = _mm_add_sd(c21_1, _mm_mul_sd(a21_1, b21));
#endif
_mm_store_sd(&C[(i*88)+42], c21_1);
__m128d c21_2 = _mm_load_sd(&C[(i*88)+70]);
__m128d a21_2 = _mm_load_sd(&values[115]);
#if defined(__SSE3__) && defined(__AVX__)
c21_2 = _mm_add_sd(c21_2, _mm_mul_sd(a21_2, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c21_2 = _mm_add_sd(c21_2, _mm_mul_sd(a21_2, b21));
#endif
_mm_store_sd(&C[(i*88)+70], c21_2);
#else
C[(i*88)+21] += values[113] * B[(i*88)+21];
C[(i*88)+42] += values[114] * B[(i*88)+21];
C[(i*88)+70] += values[115] * B[(i*88)+21];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b22 = _mm256_broadcast_sd(&B[(i*88)+22]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b22 = _mm_loaddup_pd(&B[(i*88)+22]);
#endif
__m128d c22_0 = _mm_load_sd(&C[(i*88)+22]);
__m128d a22_0 = _mm_load_sd(&values[116]);
#if defined(__SSE3__) && defined(__AVX__)
c22_0 = _mm_add_sd(c22_0, _mm_mul_sd(a22_0, _mm256_castpd256_pd128(b22)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c22_0 = _mm_add_sd(c22_0, _mm_mul_sd(a22_0, b22));
#endif
_mm_store_sd(&C[(i*88)+22], c22_0);
__m128d c22_1 = _mm_load_sd(&C[(i*88)+43]);
__m128d a22_1 = _mm_load_sd(&values[117]);
#if defined(__SSE3__) && defined(__AVX__)
c22_1 = _mm_add_sd(c22_1, _mm_mul_sd(a22_1, _mm256_castpd256_pd128(b22)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c22_1 = _mm_add_sd(c22_1, _mm_mul_sd(a22_1, b22));
#endif
_mm_store_sd(&C[(i*88)+43], c22_1);
__m128d c22_2 = _mm_load_sd(&C[(i*88)+71]);
__m128d a22_2 = _mm_load_sd(&values[118]);
#if defined(__SSE3__) && defined(__AVX__)
c22_2 = _mm_add_sd(c22_2, _mm_mul_sd(a22_2, _mm256_castpd256_pd128(b22)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c22_2 = _mm_add_sd(c22_2, _mm_mul_sd(a22_2, b22));
#endif
_mm_store_sd(&C[(i*88)+71], c22_2);
#else
C[(i*88)+22] += values[116] * B[(i*88)+22];
C[(i*88)+43] += values[117] * B[(i*88)+22];
C[(i*88)+71] += values[118] * B[(i*88)+22];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b23 = _mm256_broadcast_sd(&B[(i*88)+23]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b23 = _mm_loaddup_pd(&B[(i*88)+23]);
#endif
__m128d c23_0 = _mm_load_sd(&C[(i*88)+23]);
__m128d a23_0 = _mm_load_sd(&values[119]);
#if defined(__SSE3__) && defined(__AVX__)
c23_0 = _mm_add_sd(c23_0, _mm_mul_sd(a23_0, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c23_0 = _mm_add_sd(c23_0, _mm_mul_sd(a23_0, b23));
#endif
_mm_store_sd(&C[(i*88)+23], c23_0);
__m128d c23_1 = _mm_load_sd(&C[(i*88)+44]);
__m128d a23_1 = _mm_load_sd(&values[120]);
#if defined(__SSE3__) && defined(__AVX__)
c23_1 = _mm_add_sd(c23_1, _mm_mul_sd(a23_1, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c23_1 = _mm_add_sd(c23_1, _mm_mul_sd(a23_1, b23));
#endif
_mm_store_sd(&C[(i*88)+44], c23_1);
__m128d c23_2 = _mm_load_sd(&C[(i*88)+72]);
__m128d a23_2 = _mm_load_sd(&values[121]);
#if defined(__SSE3__) && defined(__AVX__)
c23_2 = _mm_add_sd(c23_2, _mm_mul_sd(a23_2, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c23_2 = _mm_add_sd(c23_2, _mm_mul_sd(a23_2, b23));
#endif
_mm_store_sd(&C[(i*88)+72], c23_2);
#else
C[(i*88)+23] += values[119] * B[(i*88)+23];
C[(i*88)+44] += values[120] * B[(i*88)+23];
C[(i*88)+72] += values[121] * B[(i*88)+23];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b24 = _mm256_broadcast_sd(&B[(i*88)+24]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b24 = _mm_loaddup_pd(&B[(i*88)+24]);
#endif
__m128d c24_0 = _mm_load_sd(&C[(i*88)+24]);
__m128d a24_0 = _mm_load_sd(&values[122]);
#if defined(__SSE3__) && defined(__AVX__)
c24_0 = _mm_add_sd(c24_0, _mm_mul_sd(a24_0, _mm256_castpd256_pd128(b24)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c24_0 = _mm_add_sd(c24_0, _mm_mul_sd(a24_0, b24));
#endif
_mm_store_sd(&C[(i*88)+24], c24_0);
__m128d c24_1 = _mm_load_sd(&C[(i*88)+45]);
__m128d a24_1 = _mm_load_sd(&values[123]);
#if defined(__SSE3__) && defined(__AVX__)
c24_1 = _mm_add_sd(c24_1, _mm_mul_sd(a24_1, _mm256_castpd256_pd128(b24)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c24_1 = _mm_add_sd(c24_1, _mm_mul_sd(a24_1, b24));
#endif
_mm_store_sd(&C[(i*88)+45], c24_1);
__m128d c24_2 = _mm_load_sd(&C[(i*88)+73]);
__m128d a24_2 = _mm_load_sd(&values[124]);
#if defined(__SSE3__) && defined(__AVX__)
c24_2 = _mm_add_sd(c24_2, _mm_mul_sd(a24_2, _mm256_castpd256_pd128(b24)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c24_2 = _mm_add_sd(c24_2, _mm_mul_sd(a24_2, b24));
#endif
_mm_store_sd(&C[(i*88)+73], c24_2);
#else
C[(i*88)+24] += values[122] * B[(i*88)+24];
C[(i*88)+45] += values[123] * B[(i*88)+24];
C[(i*88)+73] += values[124] * B[(i*88)+24];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b25 = _mm256_broadcast_sd(&B[(i*88)+25]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b25 = _mm_loaddup_pd(&B[(i*88)+25]);
#endif
__m128d c25_0 = _mm_load_sd(&C[(i*88)+10]);
__m128d a25_0 = _mm_load_sd(&values[125]);
#if defined(__SSE3__) && defined(__AVX__)
c25_0 = _mm_add_sd(c25_0, _mm_mul_sd(a25_0, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c25_0 = _mm_add_sd(c25_0, _mm_mul_sd(a25_0, b25));
#endif
_mm_store_sd(&C[(i*88)+10], c25_0);
__m128d c25_1 = _mm_load_sd(&C[(i*88)+25]);
__m128d a25_1 = _mm_load_sd(&values[126]);
#if defined(__SSE3__) && defined(__AVX__)
c25_1 = _mm_add_sd(c25_1, _mm_mul_sd(a25_1, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c25_1 = _mm_add_sd(c25_1, _mm_mul_sd(a25_1, b25));
#endif
_mm_store_sd(&C[(i*88)+25], c25_1);
__m128d c25_2 = _mm_load_sd(&C[(i*88)+46]);
__m128d a25_2 = _mm_load_sd(&values[127]);
#if defined(__SSE3__) && defined(__AVX__)
c25_2 = _mm_add_sd(c25_2, _mm_mul_sd(a25_2, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c25_2 = _mm_add_sd(c25_2, _mm_mul_sd(a25_2, b25));
#endif
_mm_store_sd(&C[(i*88)+46], c25_2);
__m128d c25_3 = _mm_load_sd(&C[(i*88)+74]);
__m128d a25_3 = _mm_load_sd(&values[128]);
#if defined(__SSE3__) && defined(__AVX__)
c25_3 = _mm_add_sd(c25_3, _mm_mul_sd(a25_3, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c25_3 = _mm_add_sd(c25_3, _mm_mul_sd(a25_3, b25));
#endif
_mm_store_sd(&C[(i*88)+74], c25_3);
#else
C[(i*88)+10] += values[125] * B[(i*88)+25];
C[(i*88)+25] += values[126] * B[(i*88)+25];
C[(i*88)+46] += values[127] * B[(i*88)+25];
C[(i*88)+74] += values[128] * B[(i*88)+25];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b26 = _mm256_broadcast_sd(&B[(i*88)+26]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b26 = _mm_loaddup_pd(&B[(i*88)+26]);
#endif
__m128d c26_0 = _mm_load_sd(&C[(i*88)+11]);
__m128d a26_0 = _mm_load_sd(&values[129]);
#if defined(__SSE3__) && defined(__AVX__)
c26_0 = _mm_add_sd(c26_0, _mm_mul_sd(a26_0, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c26_0 = _mm_add_sd(c26_0, _mm_mul_sd(a26_0, b26));
#endif
_mm_store_sd(&C[(i*88)+11], c26_0);
__m128d c26_1 = _mm_load_sd(&C[(i*88)+26]);
__m128d a26_1 = _mm_load_sd(&values[130]);
#if defined(__SSE3__) && defined(__AVX__)
c26_1 = _mm_add_sd(c26_1, _mm_mul_sd(a26_1, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c26_1 = _mm_add_sd(c26_1, _mm_mul_sd(a26_1, b26));
#endif
_mm_store_sd(&C[(i*88)+26], c26_1);
__m128d c26_2 = _mm_load_sd(&C[(i*88)+47]);
__m128d a26_2 = _mm_load_sd(&values[131]);
#if defined(__SSE3__) && defined(__AVX__)
c26_2 = _mm_add_sd(c26_2, _mm_mul_sd(a26_2, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c26_2 = _mm_add_sd(c26_2, _mm_mul_sd(a26_2, b26));
#endif
_mm_store_sd(&C[(i*88)+47], c26_2);
__m128d c26_3 = _mm_load_sd(&C[(i*88)+75]);
__m128d a26_3 = _mm_load_sd(&values[132]);
#if defined(__SSE3__) && defined(__AVX__)
c26_3 = _mm_add_sd(c26_3, _mm_mul_sd(a26_3, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c26_3 = _mm_add_sd(c26_3, _mm_mul_sd(a26_3, b26));
#endif
_mm_store_sd(&C[(i*88)+75], c26_3);
#else
C[(i*88)+11] += values[129] * B[(i*88)+26];
C[(i*88)+26] += values[130] * B[(i*88)+26];
C[(i*88)+47] += values[131] * B[(i*88)+26];
C[(i*88)+75] += values[132] * B[(i*88)+26];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b27 = _mm256_broadcast_sd(&B[(i*88)+27]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b27 = _mm_loaddup_pd(&B[(i*88)+27]);
#endif
__m128d c27_0 = _mm_load_sd(&C[(i*88)+12]);
__m128d a27_0 = _mm_load_sd(&values[133]);
#if defined(__SSE3__) && defined(__AVX__)
c27_0 = _mm_add_sd(c27_0, _mm_mul_sd(a27_0, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c27_0 = _mm_add_sd(c27_0, _mm_mul_sd(a27_0, b27));
#endif
_mm_store_sd(&C[(i*88)+12], c27_0);
__m128d c27_1 = _mm_load_sd(&C[(i*88)+27]);
__m128d a27_1 = _mm_load_sd(&values[134]);
#if defined(__SSE3__) && defined(__AVX__)
c27_1 = _mm_add_sd(c27_1, _mm_mul_sd(a27_1, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c27_1 = _mm_add_sd(c27_1, _mm_mul_sd(a27_1, b27));
#endif
_mm_store_sd(&C[(i*88)+27], c27_1);
__m128d c27_2 = _mm_load_sd(&C[(i*88)+48]);
__m128d a27_2 = _mm_load_sd(&values[135]);
#if defined(__SSE3__) && defined(__AVX__)
c27_2 = _mm_add_sd(c27_2, _mm_mul_sd(a27_2, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c27_2 = _mm_add_sd(c27_2, _mm_mul_sd(a27_2, b27));
#endif
_mm_store_sd(&C[(i*88)+48], c27_2);
__m128d c27_3 = _mm_load_sd(&C[(i*88)+76]);
__m128d a27_3 = _mm_load_sd(&values[136]);
#if defined(__SSE3__) && defined(__AVX__)
c27_3 = _mm_add_sd(c27_3, _mm_mul_sd(a27_3, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c27_3 = _mm_add_sd(c27_3, _mm_mul_sd(a27_3, b27));
#endif
_mm_store_sd(&C[(i*88)+76], c27_3);
#else
C[(i*88)+12] += values[133] * B[(i*88)+27];
C[(i*88)+27] += values[134] * B[(i*88)+27];
C[(i*88)+48] += values[135] * B[(i*88)+27];
C[(i*88)+76] += values[136] * B[(i*88)+27];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b28 = _mm256_broadcast_sd(&B[(i*88)+28]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b28 = _mm_loaddup_pd(&B[(i*88)+28]);
#endif
__m128d c28_0 = _mm_load_sd(&C[(i*88)+13]);
__m128d a28_0 = _mm_load_sd(&values[137]);
#if defined(__SSE3__) && defined(__AVX__)
c28_0 = _mm_add_sd(c28_0, _mm_mul_sd(a28_0, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c28_0 = _mm_add_sd(c28_0, _mm_mul_sd(a28_0, b28));
#endif
_mm_store_sd(&C[(i*88)+13], c28_0);
__m128d c28_1 = _mm_load_sd(&C[(i*88)+28]);
__m128d a28_1 = _mm_load_sd(&values[138]);
#if defined(__SSE3__) && defined(__AVX__)
c28_1 = _mm_add_sd(c28_1, _mm_mul_sd(a28_1, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c28_1 = _mm_add_sd(c28_1, _mm_mul_sd(a28_1, b28));
#endif
_mm_store_sd(&C[(i*88)+28], c28_1);
__m128d c28_2 = _mm_load_sd(&C[(i*88)+49]);
__m128d a28_2 = _mm_load_sd(&values[139]);
#if defined(__SSE3__) && defined(__AVX__)
c28_2 = _mm_add_sd(c28_2, _mm_mul_sd(a28_2, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c28_2 = _mm_add_sd(c28_2, _mm_mul_sd(a28_2, b28));
#endif
_mm_store_sd(&C[(i*88)+49], c28_2);
__m128d c28_3 = _mm_load_sd(&C[(i*88)+77]);
__m128d a28_3 = _mm_load_sd(&values[140]);
#if defined(__SSE3__) && defined(__AVX__)
c28_3 = _mm_add_sd(c28_3, _mm_mul_sd(a28_3, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c28_3 = _mm_add_sd(c28_3, _mm_mul_sd(a28_3, b28));
#endif
_mm_store_sd(&C[(i*88)+77], c28_3);
#else
C[(i*88)+13] += values[137] * B[(i*88)+28];
C[(i*88)+28] += values[138] * B[(i*88)+28];
C[(i*88)+49] += values[139] * B[(i*88)+28];
C[(i*88)+77] += values[140] * B[(i*88)+28];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b29 = _mm256_broadcast_sd(&B[(i*88)+29]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b29 = _mm_loaddup_pd(&B[(i*88)+29]);
#endif
__m128d c29_0 = _mm_load_sd(&C[(i*88)+4]);
__m128d a29_0 = _mm_load_sd(&values[141]);
#if defined(__SSE3__) && defined(__AVX__)
c29_0 = _mm_add_sd(c29_0, _mm_mul_sd(a29_0, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_0 = _mm_add_sd(c29_0, _mm_mul_sd(a29_0, b29));
#endif
_mm_store_sd(&C[(i*88)+4], c29_0);
__m128d c29_1 = _mm_load_sd(&C[(i*88)+14]);
__m128d a29_1 = _mm_load_sd(&values[142]);
#if defined(__SSE3__) && defined(__AVX__)
c29_1 = _mm_add_sd(c29_1, _mm_mul_sd(a29_1, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_1 = _mm_add_sd(c29_1, _mm_mul_sd(a29_1, b29));
#endif
_mm_store_sd(&C[(i*88)+14], c29_1);
__m128d c29_2 = _mm_load_sd(&C[(i*88)+29]);
__m128d a29_2 = _mm_load_sd(&values[143]);
#if defined(__SSE3__) && defined(__AVX__)
c29_2 = _mm_add_sd(c29_2, _mm_mul_sd(a29_2, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_2 = _mm_add_sd(c29_2, _mm_mul_sd(a29_2, b29));
#endif
_mm_store_sd(&C[(i*88)+29], c29_2);
__m128d c29_3 = _mm_load_sd(&C[(i*88)+50]);
__m128d a29_3 = _mm_load_sd(&values[144]);
#if defined(__SSE3__) && defined(__AVX__)
c29_3 = _mm_add_sd(c29_3, _mm_mul_sd(a29_3, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_3 = _mm_add_sd(c29_3, _mm_mul_sd(a29_3, b29));
#endif
_mm_store_sd(&C[(i*88)+50], c29_3);
__m128d c29_4 = _mm_load_sd(&C[(i*88)+78]);
__m128d a29_4 = _mm_load_sd(&values[145]);
#if defined(__SSE3__) && defined(__AVX__)
c29_4 = _mm_add_sd(c29_4, _mm_mul_sd(a29_4, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_4 = _mm_add_sd(c29_4, _mm_mul_sd(a29_4, b29));
#endif
_mm_store_sd(&C[(i*88)+78], c29_4);
#else
C[(i*88)+4] += values[141] * B[(i*88)+29];
C[(i*88)+14] += values[142] * B[(i*88)+29];
C[(i*88)+29] += values[143] * B[(i*88)+29];
C[(i*88)+50] += values[144] * B[(i*88)+29];
C[(i*88)+78] += values[145] * B[(i*88)+29];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b30 = _mm256_broadcast_sd(&B[(i*88)+30]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b30 = _mm_loaddup_pd(&B[(i*88)+30]);
#endif
__m128d c30_0 = _mm_load_sd(&C[(i*88)+5]);
__m128d a30_0 = _mm_load_sd(&values[146]);
#if defined(__SSE3__) && defined(__AVX__)
c30_0 = _mm_add_sd(c30_0, _mm_mul_sd(a30_0, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_0 = _mm_add_sd(c30_0, _mm_mul_sd(a30_0, b30));
#endif
_mm_store_sd(&C[(i*88)+5], c30_0);
__m128d c30_1 = _mm_load_sd(&C[(i*88)+15]);
__m128d a30_1 = _mm_load_sd(&values[147]);
#if defined(__SSE3__) && defined(__AVX__)
c30_1 = _mm_add_sd(c30_1, _mm_mul_sd(a30_1, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_1 = _mm_add_sd(c30_1, _mm_mul_sd(a30_1, b30));
#endif
_mm_store_sd(&C[(i*88)+15], c30_1);
__m128d c30_2 = _mm_load_sd(&C[(i*88)+30]);
__m128d a30_2 = _mm_load_sd(&values[148]);
#if defined(__SSE3__) && defined(__AVX__)
c30_2 = _mm_add_sd(c30_2, _mm_mul_sd(a30_2, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_2 = _mm_add_sd(c30_2, _mm_mul_sd(a30_2, b30));
#endif
_mm_store_sd(&C[(i*88)+30], c30_2);
__m128d c30_3 = _mm_load_sd(&C[(i*88)+51]);
__m128d a30_3 = _mm_load_sd(&values[149]);
#if defined(__SSE3__) && defined(__AVX__)
c30_3 = _mm_add_sd(c30_3, _mm_mul_sd(a30_3, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_3 = _mm_add_sd(c30_3, _mm_mul_sd(a30_3, b30));
#endif
_mm_store_sd(&C[(i*88)+51], c30_3);
__m128d c30_4 = _mm_load_sd(&C[(i*88)+79]);
__m128d a30_4 = _mm_load_sd(&values[150]);
#if defined(__SSE3__) && defined(__AVX__)
c30_4 = _mm_add_sd(c30_4, _mm_mul_sd(a30_4, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_4 = _mm_add_sd(c30_4, _mm_mul_sd(a30_4, b30));
#endif
_mm_store_sd(&C[(i*88)+79], c30_4);
#else
C[(i*88)+5] += values[146] * B[(i*88)+30];
C[(i*88)+15] += values[147] * B[(i*88)+30];
C[(i*88)+30] += values[148] * B[(i*88)+30];
C[(i*88)+51] += values[149] * B[(i*88)+30];
C[(i*88)+79] += values[150] * B[(i*88)+30];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b31 = _mm256_broadcast_sd(&B[(i*88)+31]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b31 = _mm_loaddup_pd(&B[(i*88)+31]);
#endif
__m128d c31_0 = _mm_load_sd(&C[(i*88)+6]);
__m128d a31_0 = _mm_load_sd(&values[151]);
#if defined(__SSE3__) && defined(__AVX__)
c31_0 = _mm_add_sd(c31_0, _mm_mul_sd(a31_0, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_0 = _mm_add_sd(c31_0, _mm_mul_sd(a31_0, b31));
#endif
_mm_store_sd(&C[(i*88)+6], c31_0);
__m128d c31_1 = _mm_load_sd(&C[(i*88)+16]);
__m128d a31_1 = _mm_load_sd(&values[152]);
#if defined(__SSE3__) && defined(__AVX__)
c31_1 = _mm_add_sd(c31_1, _mm_mul_sd(a31_1, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_1 = _mm_add_sd(c31_1, _mm_mul_sd(a31_1, b31));
#endif
_mm_store_sd(&C[(i*88)+16], c31_1);
__m128d c31_2 = _mm_load_sd(&C[(i*88)+31]);
__m128d a31_2 = _mm_load_sd(&values[153]);
#if defined(__SSE3__) && defined(__AVX__)
c31_2 = _mm_add_sd(c31_2, _mm_mul_sd(a31_2, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_2 = _mm_add_sd(c31_2, _mm_mul_sd(a31_2, b31));
#endif
_mm_store_sd(&C[(i*88)+31], c31_2);
__m128d c31_3 = _mm_load_sd(&C[(i*88)+52]);
__m128d a31_3 = _mm_load_sd(&values[154]);
#if defined(__SSE3__) && defined(__AVX__)
c31_3 = _mm_add_sd(c31_3, _mm_mul_sd(a31_3, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_3 = _mm_add_sd(c31_3, _mm_mul_sd(a31_3, b31));
#endif
_mm_store_sd(&C[(i*88)+52], c31_3);
__m128d c31_4 = _mm_load_sd(&C[(i*88)+80]);
__m128d a31_4 = _mm_load_sd(&values[155]);
#if defined(__SSE3__) && defined(__AVX__)
c31_4 = _mm_add_sd(c31_4, _mm_mul_sd(a31_4, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_4 = _mm_add_sd(c31_4, _mm_mul_sd(a31_4, b31));
#endif
_mm_store_sd(&C[(i*88)+80], c31_4);
#else
C[(i*88)+6] += values[151] * B[(i*88)+31];
C[(i*88)+16] += values[152] * B[(i*88)+31];
C[(i*88)+31] += values[153] * B[(i*88)+31];
C[(i*88)+52] += values[154] * B[(i*88)+31];
C[(i*88)+80] += values[155] * B[(i*88)+31];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b32 = _mm256_broadcast_sd(&B[(i*88)+32]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b32 = _mm_loaddup_pd(&B[(i*88)+32]);
#endif
__m128d c32_0 = _mm_load_sd(&C[(i*88)+1]);
__m128d a32_0 = _mm_load_sd(&values[156]);
#if defined(__SSE3__) && defined(__AVX__)
c32_0 = _mm_add_sd(c32_0, _mm_mul_sd(a32_0, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_0 = _mm_add_sd(c32_0, _mm_mul_sd(a32_0, b32));
#endif
_mm_store_sd(&C[(i*88)+1], c32_0);
__m128d c32_1 = _mm_load_sd(&C[(i*88)+7]);
__m128d a32_1 = _mm_load_sd(&values[157]);
#if defined(__SSE3__) && defined(__AVX__)
c32_1 = _mm_add_sd(c32_1, _mm_mul_sd(a32_1, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_1 = _mm_add_sd(c32_1, _mm_mul_sd(a32_1, b32));
#endif
_mm_store_sd(&C[(i*88)+7], c32_1);
__m128d c32_2 = _mm_load_sd(&C[(i*88)+17]);
__m128d a32_2 = _mm_load_sd(&values[158]);
#if defined(__SSE3__) && defined(__AVX__)
c32_2 = _mm_add_sd(c32_2, _mm_mul_sd(a32_2, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_2 = _mm_add_sd(c32_2, _mm_mul_sd(a32_2, b32));
#endif
_mm_store_sd(&C[(i*88)+17], c32_2);
__m128d c32_3 = _mm_load_sd(&C[(i*88)+32]);
__m128d a32_3 = _mm_load_sd(&values[159]);
#if defined(__SSE3__) && defined(__AVX__)
c32_3 = _mm_add_sd(c32_3, _mm_mul_sd(a32_3, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_3 = _mm_add_sd(c32_3, _mm_mul_sd(a32_3, b32));
#endif
_mm_store_sd(&C[(i*88)+32], c32_3);
__m128d c32_4 = _mm_load_sd(&C[(i*88)+53]);
__m128d a32_4 = _mm_load_sd(&values[160]);
#if defined(__SSE3__) && defined(__AVX__)
c32_4 = _mm_add_sd(c32_4, _mm_mul_sd(a32_4, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_4 = _mm_add_sd(c32_4, _mm_mul_sd(a32_4, b32));
#endif
_mm_store_sd(&C[(i*88)+53], c32_4);
__m128d c32_5 = _mm_load_sd(&C[(i*88)+81]);
__m128d a32_5 = _mm_load_sd(&values[161]);
#if defined(__SSE3__) && defined(__AVX__)
c32_5 = _mm_add_sd(c32_5, _mm_mul_sd(a32_5, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_5 = _mm_add_sd(c32_5, _mm_mul_sd(a32_5, b32));
#endif
_mm_store_sd(&C[(i*88)+81], c32_5);
#else
C[(i*88)+1] += values[156] * B[(i*88)+32];
C[(i*88)+7] += values[157] * B[(i*88)+32];
C[(i*88)+17] += values[158] * B[(i*88)+32];
C[(i*88)+32] += values[159] * B[(i*88)+32];
C[(i*88)+53] += values[160] * B[(i*88)+32];
C[(i*88)+81] += values[161] * B[(i*88)+32];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b33 = _mm256_broadcast_sd(&B[(i*88)+33]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b33 = _mm_loaddup_pd(&B[(i*88)+33]);
#endif
__m128d c33_0 = _mm_load_sd(&C[(i*88)+2]);
__m128d a33_0 = _mm_load_sd(&values[162]);
#if defined(__SSE3__) && defined(__AVX__)
c33_0 = _mm_add_sd(c33_0, _mm_mul_sd(a33_0, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_0 = _mm_add_sd(c33_0, _mm_mul_sd(a33_0, b33));
#endif
_mm_store_sd(&C[(i*88)+2], c33_0);
__m128d c33_1 = _mm_load_sd(&C[(i*88)+8]);
__m128d a33_1 = _mm_load_sd(&values[163]);
#if defined(__SSE3__) && defined(__AVX__)
c33_1 = _mm_add_sd(c33_1, _mm_mul_sd(a33_1, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_1 = _mm_add_sd(c33_1, _mm_mul_sd(a33_1, b33));
#endif
_mm_store_sd(&C[(i*88)+8], c33_1);
__m128d c33_2 = _mm_load_sd(&C[(i*88)+18]);
__m128d a33_2 = _mm_load_sd(&values[164]);
#if defined(__SSE3__) && defined(__AVX__)
c33_2 = _mm_add_sd(c33_2, _mm_mul_sd(a33_2, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_2 = _mm_add_sd(c33_2, _mm_mul_sd(a33_2, b33));
#endif
_mm_store_sd(&C[(i*88)+18], c33_2);
__m128d c33_3 = _mm_load_sd(&C[(i*88)+33]);
__m128d a33_3 = _mm_load_sd(&values[165]);
#if defined(__SSE3__) && defined(__AVX__)
c33_3 = _mm_add_sd(c33_3, _mm_mul_sd(a33_3, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_3 = _mm_add_sd(c33_3, _mm_mul_sd(a33_3, b33));
#endif
_mm_store_sd(&C[(i*88)+33], c33_3);
__m128d c33_4 = _mm_load_sd(&C[(i*88)+54]);
__m128d a33_4 = _mm_load_sd(&values[166]);
#if defined(__SSE3__) && defined(__AVX__)
c33_4 = _mm_add_sd(c33_4, _mm_mul_sd(a33_4, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_4 = _mm_add_sd(c33_4, _mm_mul_sd(a33_4, b33));
#endif
_mm_store_sd(&C[(i*88)+54], c33_4);
__m128d c33_5 = _mm_load_sd(&C[(i*88)+82]);
__m128d a33_5 = _mm_load_sd(&values[167]);
#if defined(__SSE3__) && defined(__AVX__)
c33_5 = _mm_add_sd(c33_5, _mm_mul_sd(a33_5, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_5 = _mm_add_sd(c33_5, _mm_mul_sd(a33_5, b33));
#endif
_mm_store_sd(&C[(i*88)+82], c33_5);
#else
C[(i*88)+2] += values[162] * B[(i*88)+33];
C[(i*88)+8] += values[163] * B[(i*88)+33];
C[(i*88)+18] += values[164] * B[(i*88)+33];
C[(i*88)+33] += values[165] * B[(i*88)+33];
C[(i*88)+54] += values[166] * B[(i*88)+33];
C[(i*88)+82] += values[167] * B[(i*88)+33];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b34 = _mm256_broadcast_sd(&B[(i*88)+34]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b34 = _mm_loaddup_pd(&B[(i*88)+34]);
#endif
__m128d c34_0 = _mm_load_sd(&C[(i*88)+0]);
__m128d a34_0 = _mm_load_sd(&values[168]);
#if defined(__SSE3__) && defined(__AVX__)
c34_0 = _mm_add_sd(c34_0, _mm_mul_sd(a34_0, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_0 = _mm_add_sd(c34_0, _mm_mul_sd(a34_0, b34));
#endif
_mm_store_sd(&C[(i*88)+0], c34_0);
__m128d c34_1 = _mm_load_sd(&C[(i*88)+3]);
__m128d a34_1 = _mm_load_sd(&values[169]);
#if defined(__SSE3__) && defined(__AVX__)
c34_1 = _mm_add_sd(c34_1, _mm_mul_sd(a34_1, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_1 = _mm_add_sd(c34_1, _mm_mul_sd(a34_1, b34));
#endif
_mm_store_sd(&C[(i*88)+3], c34_1);
__m128d c34_2 = _mm_load_sd(&C[(i*88)+9]);
__m128d a34_2 = _mm_load_sd(&values[170]);
#if defined(__SSE3__) && defined(__AVX__)
c34_2 = _mm_add_sd(c34_2, _mm_mul_sd(a34_2, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_2 = _mm_add_sd(c34_2, _mm_mul_sd(a34_2, b34));
#endif
_mm_store_sd(&C[(i*88)+9], c34_2);
__m128d c34_3 = _mm_load_sd(&C[(i*88)+19]);
__m128d a34_3 = _mm_load_sd(&values[171]);
#if defined(__SSE3__) && defined(__AVX__)
c34_3 = _mm_add_sd(c34_3, _mm_mul_sd(a34_3, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_3 = _mm_add_sd(c34_3, _mm_mul_sd(a34_3, b34));
#endif
_mm_store_sd(&C[(i*88)+19], c34_3);
__m128d c34_4 = _mm_load_sd(&C[(i*88)+34]);
__m128d a34_4 = _mm_load_sd(&values[172]);
#if defined(__SSE3__) && defined(__AVX__)
c34_4 = _mm_add_sd(c34_4, _mm_mul_sd(a34_4, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_4 = _mm_add_sd(c34_4, _mm_mul_sd(a34_4, b34));
#endif
_mm_store_sd(&C[(i*88)+34], c34_4);
__m128d c34_5 = _mm_load_sd(&C[(i*88)+55]);
__m128d a34_5 = _mm_load_sd(&values[173]);
#if defined(__SSE3__) && defined(__AVX__)
c34_5 = _mm_add_sd(c34_5, _mm_mul_sd(a34_5, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_5 = _mm_add_sd(c34_5, _mm_mul_sd(a34_5, b34));
#endif
_mm_store_sd(&C[(i*88)+55], c34_5);
__m128d c34_6 = _mm_load_sd(&C[(i*88)+83]);
__m128d a34_6 = _mm_load_sd(&values[174]);
#if defined(__SSE3__) && defined(__AVX__)
c34_6 = _mm_add_sd(c34_6, _mm_mul_sd(a34_6, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_6 = _mm_add_sd(c34_6, _mm_mul_sd(a34_6, b34));
#endif
_mm_store_sd(&C[(i*88)+83], c34_6);
#else
C[(i*88)+0] += values[168] * B[(i*88)+34];
C[(i*88)+3] += values[169] * B[(i*88)+34];
C[(i*88)+9] += values[170] * B[(i*88)+34];
C[(i*88)+19] += values[171] * B[(i*88)+34];
C[(i*88)+34] += values[172] * B[(i*88)+34];
C[(i*88)+55] += values[173] * B[(i*88)+34];
C[(i*88)+83] += values[174] * B[(i*88)+34];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b35 = _mm256_broadcast_sd(&B[(i*88)+35]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b35 = _mm_loaddup_pd(&B[(i*88)+35]);
#endif
__m128d c35_0 = _mm_load_sd(&C[(i*88)+35]);
__m128d a35_0 = _mm_load_sd(&values[175]);
#if defined(__SSE3__) && defined(__AVX__)
c35_0 = _mm_add_sd(c35_0, _mm_mul_sd(a35_0, _mm256_castpd256_pd128(b35)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c35_0 = _mm_add_sd(c35_0, _mm_mul_sd(a35_0, b35));
#endif
_mm_store_sd(&C[(i*88)+35], c35_0);
__m128d c35_1 = _mm_load_sd(&C[(i*88)+63]);
__m128d a35_1 = _mm_load_sd(&values[176]);
#if defined(__SSE3__) && defined(__AVX__)
c35_1 = _mm_add_sd(c35_1, _mm_mul_sd(a35_1, _mm256_castpd256_pd128(b35)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c35_1 = _mm_add_sd(c35_1, _mm_mul_sd(a35_1, b35));
#endif
_mm_store_sd(&C[(i*88)+63], c35_1);
#else
C[(i*88)+35] += values[175] * B[(i*88)+35];
C[(i*88)+63] += values[176] * B[(i*88)+35];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b36 = _mm256_broadcast_sd(&B[(i*88)+36]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b36 = _mm_loaddup_pd(&B[(i*88)+36]);
#endif
__m128d c36_0 = _mm_load_sd(&C[(i*88)+36]);
__m128d a36_0 = _mm_load_sd(&values[177]);
#if defined(__SSE3__) && defined(__AVX__)
c36_0 = _mm_add_sd(c36_0, _mm_mul_sd(a36_0, _mm256_castpd256_pd128(b36)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c36_0 = _mm_add_sd(c36_0, _mm_mul_sd(a36_0, b36));
#endif
_mm_store_sd(&C[(i*88)+36], c36_0);
__m128d c36_1 = _mm_load_sd(&C[(i*88)+64]);
__m128d a36_1 = _mm_load_sd(&values[178]);
#if defined(__SSE3__) && defined(__AVX__)
c36_1 = _mm_add_sd(c36_1, _mm_mul_sd(a36_1, _mm256_castpd256_pd128(b36)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c36_1 = _mm_add_sd(c36_1, _mm_mul_sd(a36_1, b36));
#endif
_mm_store_sd(&C[(i*88)+64], c36_1);
#else
C[(i*88)+36] += values[177] * B[(i*88)+36];
C[(i*88)+64] += values[178] * B[(i*88)+36];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b37 = _mm256_broadcast_sd(&B[(i*88)+37]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b37 = _mm_loaddup_pd(&B[(i*88)+37]);
#endif
__m128d c37_0 = _mm_load_sd(&C[(i*88)+37]);
__m128d a37_0 = _mm_load_sd(&values[179]);
#if defined(__SSE3__) && defined(__AVX__)
c37_0 = _mm_add_sd(c37_0, _mm_mul_sd(a37_0, _mm256_castpd256_pd128(b37)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c37_0 = _mm_add_sd(c37_0, _mm_mul_sd(a37_0, b37));
#endif
_mm_store_sd(&C[(i*88)+37], c37_0);
__m128d c37_1 = _mm_load_sd(&C[(i*88)+65]);
__m128d a37_1 = _mm_load_sd(&values[180]);
#if defined(__SSE3__) && defined(__AVX__)
c37_1 = _mm_add_sd(c37_1, _mm_mul_sd(a37_1, _mm256_castpd256_pd128(b37)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c37_1 = _mm_add_sd(c37_1, _mm_mul_sd(a37_1, b37));
#endif
_mm_store_sd(&C[(i*88)+65], c37_1);
#else
C[(i*88)+37] += values[179] * B[(i*88)+37];
C[(i*88)+65] += values[180] * B[(i*88)+37];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b38 = _mm256_broadcast_sd(&B[(i*88)+38]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b38 = _mm_loaddup_pd(&B[(i*88)+38]);
#endif
__m128d c38_0 = _mm_load_sd(&C[(i*88)+38]);
__m128d a38_0 = _mm_load_sd(&values[181]);
#if defined(__SSE3__) && defined(__AVX__)
c38_0 = _mm_add_sd(c38_0, _mm_mul_sd(a38_0, _mm256_castpd256_pd128(b38)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c38_0 = _mm_add_sd(c38_0, _mm_mul_sd(a38_0, b38));
#endif
_mm_store_sd(&C[(i*88)+38], c38_0);
__m128d c38_1 = _mm_load_sd(&C[(i*88)+66]);
__m128d a38_1 = _mm_load_sd(&values[182]);
#if defined(__SSE3__) && defined(__AVX__)
c38_1 = _mm_add_sd(c38_1, _mm_mul_sd(a38_1, _mm256_castpd256_pd128(b38)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c38_1 = _mm_add_sd(c38_1, _mm_mul_sd(a38_1, b38));
#endif
_mm_store_sd(&C[(i*88)+66], c38_1);
#else
C[(i*88)+38] += values[181] * B[(i*88)+38];
C[(i*88)+66] += values[182] * B[(i*88)+38];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b39 = _mm256_broadcast_sd(&B[(i*88)+39]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b39 = _mm_loaddup_pd(&B[(i*88)+39]);
#endif
__m128d c39_0 = _mm_load_sd(&C[(i*88)+39]);
__m128d a39_0 = _mm_load_sd(&values[183]);
#if defined(__SSE3__) && defined(__AVX__)
c39_0 = _mm_add_sd(c39_0, _mm_mul_sd(a39_0, _mm256_castpd256_pd128(b39)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c39_0 = _mm_add_sd(c39_0, _mm_mul_sd(a39_0, b39));
#endif
_mm_store_sd(&C[(i*88)+39], c39_0);
__m128d c39_1 = _mm_load_sd(&C[(i*88)+67]);
__m128d a39_1 = _mm_load_sd(&values[184]);
#if defined(__SSE3__) && defined(__AVX__)
c39_1 = _mm_add_sd(c39_1, _mm_mul_sd(a39_1, _mm256_castpd256_pd128(b39)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c39_1 = _mm_add_sd(c39_1, _mm_mul_sd(a39_1, b39));
#endif
_mm_store_sd(&C[(i*88)+67], c39_1);
#else
C[(i*88)+39] += values[183] * B[(i*88)+39];
C[(i*88)+67] += values[184] * B[(i*88)+39];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b40 = _mm256_broadcast_sd(&B[(i*88)+40]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b40 = _mm_loaddup_pd(&B[(i*88)+40]);
#endif
__m128d c40_0 = _mm_load_sd(&C[(i*88)+40]);
__m128d a40_0 = _mm_load_sd(&values[185]);
#if defined(__SSE3__) && defined(__AVX__)
c40_0 = _mm_add_sd(c40_0, _mm_mul_sd(a40_0, _mm256_castpd256_pd128(b40)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c40_0 = _mm_add_sd(c40_0, _mm_mul_sd(a40_0, b40));
#endif
_mm_store_sd(&C[(i*88)+40], c40_0);
__m128d c40_1 = _mm_load_sd(&C[(i*88)+68]);
__m128d a40_1 = _mm_load_sd(&values[186]);
#if defined(__SSE3__) && defined(__AVX__)
c40_1 = _mm_add_sd(c40_1, _mm_mul_sd(a40_1, _mm256_castpd256_pd128(b40)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c40_1 = _mm_add_sd(c40_1, _mm_mul_sd(a40_1, b40));
#endif
_mm_store_sd(&C[(i*88)+68], c40_1);
#else
C[(i*88)+40] += values[185] * B[(i*88)+40];
C[(i*88)+68] += values[186] * B[(i*88)+40];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b41 = _mm256_broadcast_sd(&B[(i*88)+41]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b41 = _mm_loaddup_pd(&B[(i*88)+41]);
#endif
__m128d c41_0 = _mm_load_sd(&C[(i*88)+20]);
__m128d a41_0 = _mm_load_sd(&values[187]);
#if defined(__SSE3__) && defined(__AVX__)
c41_0 = _mm_add_sd(c41_0, _mm_mul_sd(a41_0, _mm256_castpd256_pd128(b41)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c41_0 = _mm_add_sd(c41_0, _mm_mul_sd(a41_0, b41));
#endif
_mm_store_sd(&C[(i*88)+20], c41_0);
__m128d c41_1 = _mm_load_sd(&C[(i*88)+41]);
__m128d a41_1 = _mm_load_sd(&values[188]);
#if defined(__SSE3__) && defined(__AVX__)
c41_1 = _mm_add_sd(c41_1, _mm_mul_sd(a41_1, _mm256_castpd256_pd128(b41)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c41_1 = _mm_add_sd(c41_1, _mm_mul_sd(a41_1, b41));
#endif
_mm_store_sd(&C[(i*88)+41], c41_1);
__m128d c41_2 = _mm_load_sd(&C[(i*88)+69]);
__m128d a41_2 = _mm_load_sd(&values[189]);
#if defined(__SSE3__) && defined(__AVX__)
c41_2 = _mm_add_sd(c41_2, _mm_mul_sd(a41_2, _mm256_castpd256_pd128(b41)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c41_2 = _mm_add_sd(c41_2, _mm_mul_sd(a41_2, b41));
#endif
_mm_store_sd(&C[(i*88)+69], c41_2);
#else
C[(i*88)+20] += values[187] * B[(i*88)+41];
C[(i*88)+41] += values[188] * B[(i*88)+41];
C[(i*88)+69] += values[189] * B[(i*88)+41];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b42 = _mm256_broadcast_sd(&B[(i*88)+42]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b42 = _mm_loaddup_pd(&B[(i*88)+42]);
#endif
__m128d c42_0 = _mm_load_sd(&C[(i*88)+21]);
__m128d a42_0 = _mm_load_sd(&values[190]);
#if defined(__SSE3__) && defined(__AVX__)
c42_0 = _mm_add_sd(c42_0, _mm_mul_sd(a42_0, _mm256_castpd256_pd128(b42)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c42_0 = _mm_add_sd(c42_0, _mm_mul_sd(a42_0, b42));
#endif
_mm_store_sd(&C[(i*88)+21], c42_0);
__m128d c42_1 = _mm_load_sd(&C[(i*88)+42]);
__m128d a42_1 = _mm_load_sd(&values[191]);
#if defined(__SSE3__) && defined(__AVX__)
c42_1 = _mm_add_sd(c42_1, _mm_mul_sd(a42_1, _mm256_castpd256_pd128(b42)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c42_1 = _mm_add_sd(c42_1, _mm_mul_sd(a42_1, b42));
#endif
_mm_store_sd(&C[(i*88)+42], c42_1);
__m128d c42_2 = _mm_load_sd(&C[(i*88)+70]);
__m128d a42_2 = _mm_load_sd(&values[192]);
#if defined(__SSE3__) && defined(__AVX__)
c42_2 = _mm_add_sd(c42_2, _mm_mul_sd(a42_2, _mm256_castpd256_pd128(b42)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c42_2 = _mm_add_sd(c42_2, _mm_mul_sd(a42_2, b42));
#endif
_mm_store_sd(&C[(i*88)+70], c42_2);
#else
C[(i*88)+21] += values[190] * B[(i*88)+42];
C[(i*88)+42] += values[191] * B[(i*88)+42];
C[(i*88)+70] += values[192] * B[(i*88)+42];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b43 = _mm256_broadcast_sd(&B[(i*88)+43]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b43 = _mm_loaddup_pd(&B[(i*88)+43]);
#endif
__m128d c43_0 = _mm_load_sd(&C[(i*88)+22]);
__m128d a43_0 = _mm_load_sd(&values[193]);
#if defined(__SSE3__) && defined(__AVX__)
c43_0 = _mm_add_sd(c43_0, _mm_mul_sd(a43_0, _mm256_castpd256_pd128(b43)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c43_0 = _mm_add_sd(c43_0, _mm_mul_sd(a43_0, b43));
#endif
_mm_store_sd(&C[(i*88)+22], c43_0);
__m128d c43_1 = _mm_load_sd(&C[(i*88)+43]);
__m128d a43_1 = _mm_load_sd(&values[194]);
#if defined(__SSE3__) && defined(__AVX__)
c43_1 = _mm_add_sd(c43_1, _mm_mul_sd(a43_1, _mm256_castpd256_pd128(b43)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c43_1 = _mm_add_sd(c43_1, _mm_mul_sd(a43_1, b43));
#endif
_mm_store_sd(&C[(i*88)+43], c43_1);
__m128d c43_2 = _mm_load_sd(&C[(i*88)+71]);
__m128d a43_2 = _mm_load_sd(&values[195]);
#if defined(__SSE3__) && defined(__AVX__)
c43_2 = _mm_add_sd(c43_2, _mm_mul_sd(a43_2, _mm256_castpd256_pd128(b43)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c43_2 = _mm_add_sd(c43_2, _mm_mul_sd(a43_2, b43));
#endif
_mm_store_sd(&C[(i*88)+71], c43_2);
#else
C[(i*88)+22] += values[193] * B[(i*88)+43];
C[(i*88)+43] += values[194] * B[(i*88)+43];
C[(i*88)+71] += values[195] * B[(i*88)+43];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b44 = _mm256_broadcast_sd(&B[(i*88)+44]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b44 = _mm_loaddup_pd(&B[(i*88)+44]);
#endif
__m128d c44_0 = _mm_load_sd(&C[(i*88)+23]);
__m128d a44_0 = _mm_load_sd(&values[196]);
#if defined(__SSE3__) && defined(__AVX__)
c44_0 = _mm_add_sd(c44_0, _mm_mul_sd(a44_0, _mm256_castpd256_pd128(b44)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c44_0 = _mm_add_sd(c44_0, _mm_mul_sd(a44_0, b44));
#endif
_mm_store_sd(&C[(i*88)+23], c44_0);
__m128d c44_1 = _mm_load_sd(&C[(i*88)+44]);
__m128d a44_1 = _mm_load_sd(&values[197]);
#if defined(__SSE3__) && defined(__AVX__)
c44_1 = _mm_add_sd(c44_1, _mm_mul_sd(a44_1, _mm256_castpd256_pd128(b44)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c44_1 = _mm_add_sd(c44_1, _mm_mul_sd(a44_1, b44));
#endif
_mm_store_sd(&C[(i*88)+44], c44_1);
__m128d c44_2 = _mm_load_sd(&C[(i*88)+72]);
__m128d a44_2 = _mm_load_sd(&values[198]);
#if defined(__SSE3__) && defined(__AVX__)
c44_2 = _mm_add_sd(c44_2, _mm_mul_sd(a44_2, _mm256_castpd256_pd128(b44)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c44_2 = _mm_add_sd(c44_2, _mm_mul_sd(a44_2, b44));
#endif
_mm_store_sd(&C[(i*88)+72], c44_2);
#else
C[(i*88)+23] += values[196] * B[(i*88)+44];
C[(i*88)+44] += values[197] * B[(i*88)+44];
C[(i*88)+72] += values[198] * B[(i*88)+44];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b45 = _mm256_broadcast_sd(&B[(i*88)+45]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b45 = _mm_loaddup_pd(&B[(i*88)+45]);
#endif
__m128d c45_0 = _mm_load_sd(&C[(i*88)+24]);
__m128d a45_0 = _mm_load_sd(&values[199]);
#if defined(__SSE3__) && defined(__AVX__)
c45_0 = _mm_add_sd(c45_0, _mm_mul_sd(a45_0, _mm256_castpd256_pd128(b45)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c45_0 = _mm_add_sd(c45_0, _mm_mul_sd(a45_0, b45));
#endif
_mm_store_sd(&C[(i*88)+24], c45_0);
__m128d c45_1 = _mm_load_sd(&C[(i*88)+45]);
__m128d a45_1 = _mm_load_sd(&values[200]);
#if defined(__SSE3__) && defined(__AVX__)
c45_1 = _mm_add_sd(c45_1, _mm_mul_sd(a45_1, _mm256_castpd256_pd128(b45)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c45_1 = _mm_add_sd(c45_1, _mm_mul_sd(a45_1, b45));
#endif
_mm_store_sd(&C[(i*88)+45], c45_1);
__m128d c45_2 = _mm_load_sd(&C[(i*88)+73]);
__m128d a45_2 = _mm_load_sd(&values[201]);
#if defined(__SSE3__) && defined(__AVX__)
c45_2 = _mm_add_sd(c45_2, _mm_mul_sd(a45_2, _mm256_castpd256_pd128(b45)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c45_2 = _mm_add_sd(c45_2, _mm_mul_sd(a45_2, b45));
#endif
_mm_store_sd(&C[(i*88)+73], c45_2);
#else
C[(i*88)+24] += values[199] * B[(i*88)+45];
C[(i*88)+45] += values[200] * B[(i*88)+45];
C[(i*88)+73] += values[201] * B[(i*88)+45];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b46 = _mm256_broadcast_sd(&B[(i*88)+46]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b46 = _mm_loaddup_pd(&B[(i*88)+46]);
#endif
__m128d c46_0 = _mm_load_sd(&C[(i*88)+10]);
__m128d a46_0 = _mm_load_sd(&values[202]);
#if defined(__SSE3__) && defined(__AVX__)
c46_0 = _mm_add_sd(c46_0, _mm_mul_sd(a46_0, _mm256_castpd256_pd128(b46)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c46_0 = _mm_add_sd(c46_0, _mm_mul_sd(a46_0, b46));
#endif
_mm_store_sd(&C[(i*88)+10], c46_0);
__m128d c46_1 = _mm_load_sd(&C[(i*88)+25]);
__m128d a46_1 = _mm_load_sd(&values[203]);
#if defined(__SSE3__) && defined(__AVX__)
c46_1 = _mm_add_sd(c46_1, _mm_mul_sd(a46_1, _mm256_castpd256_pd128(b46)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c46_1 = _mm_add_sd(c46_1, _mm_mul_sd(a46_1, b46));
#endif
_mm_store_sd(&C[(i*88)+25], c46_1);
__m128d c46_2 = _mm_load_sd(&C[(i*88)+46]);
__m128d a46_2 = _mm_load_sd(&values[204]);
#if defined(__SSE3__) && defined(__AVX__)
c46_2 = _mm_add_sd(c46_2, _mm_mul_sd(a46_2, _mm256_castpd256_pd128(b46)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c46_2 = _mm_add_sd(c46_2, _mm_mul_sd(a46_2, b46));
#endif
_mm_store_sd(&C[(i*88)+46], c46_2);
__m128d c46_3 = _mm_load_sd(&C[(i*88)+74]);
__m128d a46_3 = _mm_load_sd(&values[205]);
#if defined(__SSE3__) && defined(__AVX__)
c46_3 = _mm_add_sd(c46_3, _mm_mul_sd(a46_3, _mm256_castpd256_pd128(b46)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c46_3 = _mm_add_sd(c46_3, _mm_mul_sd(a46_3, b46));
#endif
_mm_store_sd(&C[(i*88)+74], c46_3);
#else
C[(i*88)+10] += values[202] * B[(i*88)+46];
C[(i*88)+25] += values[203] * B[(i*88)+46];
C[(i*88)+46] += values[204] * B[(i*88)+46];
C[(i*88)+74] += values[205] * B[(i*88)+46];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b47 = _mm256_broadcast_sd(&B[(i*88)+47]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b47 = _mm_loaddup_pd(&B[(i*88)+47]);
#endif
__m128d c47_0 = _mm_load_sd(&C[(i*88)+11]);
__m128d a47_0 = _mm_load_sd(&values[206]);
#if defined(__SSE3__) && defined(__AVX__)
c47_0 = _mm_add_sd(c47_0, _mm_mul_sd(a47_0, _mm256_castpd256_pd128(b47)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c47_0 = _mm_add_sd(c47_0, _mm_mul_sd(a47_0, b47));
#endif
_mm_store_sd(&C[(i*88)+11], c47_0);
__m128d c47_1 = _mm_load_sd(&C[(i*88)+26]);
__m128d a47_1 = _mm_load_sd(&values[207]);
#if defined(__SSE3__) && defined(__AVX__)
c47_1 = _mm_add_sd(c47_1, _mm_mul_sd(a47_1, _mm256_castpd256_pd128(b47)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c47_1 = _mm_add_sd(c47_1, _mm_mul_sd(a47_1, b47));
#endif
_mm_store_sd(&C[(i*88)+26], c47_1);
__m128d c47_2 = _mm_load_sd(&C[(i*88)+47]);
__m128d a47_2 = _mm_load_sd(&values[208]);
#if defined(__SSE3__) && defined(__AVX__)
c47_2 = _mm_add_sd(c47_2, _mm_mul_sd(a47_2, _mm256_castpd256_pd128(b47)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c47_2 = _mm_add_sd(c47_2, _mm_mul_sd(a47_2, b47));
#endif
_mm_store_sd(&C[(i*88)+47], c47_2);
__m128d c47_3 = _mm_load_sd(&C[(i*88)+75]);
__m128d a47_3 = _mm_load_sd(&values[209]);
#if defined(__SSE3__) && defined(__AVX__)
c47_3 = _mm_add_sd(c47_3, _mm_mul_sd(a47_3, _mm256_castpd256_pd128(b47)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c47_3 = _mm_add_sd(c47_3, _mm_mul_sd(a47_3, b47));
#endif
_mm_store_sd(&C[(i*88)+75], c47_3);
#else
C[(i*88)+11] += values[206] * B[(i*88)+47];
C[(i*88)+26] += values[207] * B[(i*88)+47];
C[(i*88)+47] += values[208] * B[(i*88)+47];
C[(i*88)+75] += values[209] * B[(i*88)+47];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b48 = _mm256_broadcast_sd(&B[(i*88)+48]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b48 = _mm_loaddup_pd(&B[(i*88)+48]);
#endif
__m128d c48_0 = _mm_load_sd(&C[(i*88)+12]);
__m128d a48_0 = _mm_load_sd(&values[210]);
#if defined(__SSE3__) && defined(__AVX__)
c48_0 = _mm_add_sd(c48_0, _mm_mul_sd(a48_0, _mm256_castpd256_pd128(b48)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c48_0 = _mm_add_sd(c48_0, _mm_mul_sd(a48_0, b48));
#endif
_mm_store_sd(&C[(i*88)+12], c48_0);
__m128d c48_1 = _mm_load_sd(&C[(i*88)+27]);
__m128d a48_1 = _mm_load_sd(&values[211]);
#if defined(__SSE3__) && defined(__AVX__)
c48_1 = _mm_add_sd(c48_1, _mm_mul_sd(a48_1, _mm256_castpd256_pd128(b48)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c48_1 = _mm_add_sd(c48_1, _mm_mul_sd(a48_1, b48));
#endif
_mm_store_sd(&C[(i*88)+27], c48_1);
__m128d c48_2 = _mm_load_sd(&C[(i*88)+48]);
__m128d a48_2 = _mm_load_sd(&values[212]);
#if defined(__SSE3__) && defined(__AVX__)
c48_2 = _mm_add_sd(c48_2, _mm_mul_sd(a48_2, _mm256_castpd256_pd128(b48)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c48_2 = _mm_add_sd(c48_2, _mm_mul_sd(a48_2, b48));
#endif
_mm_store_sd(&C[(i*88)+48], c48_2);
__m128d c48_3 = _mm_load_sd(&C[(i*88)+76]);
__m128d a48_3 = _mm_load_sd(&values[213]);
#if defined(__SSE3__) && defined(__AVX__)
c48_3 = _mm_add_sd(c48_3, _mm_mul_sd(a48_3, _mm256_castpd256_pd128(b48)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c48_3 = _mm_add_sd(c48_3, _mm_mul_sd(a48_3, b48));
#endif
_mm_store_sd(&C[(i*88)+76], c48_3);
#else
C[(i*88)+12] += values[210] * B[(i*88)+48];
C[(i*88)+27] += values[211] * B[(i*88)+48];
C[(i*88)+48] += values[212] * B[(i*88)+48];
C[(i*88)+76] += values[213] * B[(i*88)+48];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b49 = _mm256_broadcast_sd(&B[(i*88)+49]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b49 = _mm_loaddup_pd(&B[(i*88)+49]);
#endif
__m128d c49_0 = _mm_load_sd(&C[(i*88)+13]);
__m128d a49_0 = _mm_load_sd(&values[214]);
#if defined(__SSE3__) && defined(__AVX__)
c49_0 = _mm_add_sd(c49_0, _mm_mul_sd(a49_0, _mm256_castpd256_pd128(b49)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c49_0 = _mm_add_sd(c49_0, _mm_mul_sd(a49_0, b49));
#endif
_mm_store_sd(&C[(i*88)+13], c49_0);
__m128d c49_1 = _mm_load_sd(&C[(i*88)+28]);
__m128d a49_1 = _mm_load_sd(&values[215]);
#if defined(__SSE3__) && defined(__AVX__)
c49_1 = _mm_add_sd(c49_1, _mm_mul_sd(a49_1, _mm256_castpd256_pd128(b49)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c49_1 = _mm_add_sd(c49_1, _mm_mul_sd(a49_1, b49));
#endif
_mm_store_sd(&C[(i*88)+28], c49_1);
__m128d c49_2 = _mm_load_sd(&C[(i*88)+49]);
__m128d a49_2 = _mm_load_sd(&values[216]);
#if defined(__SSE3__) && defined(__AVX__)
c49_2 = _mm_add_sd(c49_2, _mm_mul_sd(a49_2, _mm256_castpd256_pd128(b49)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c49_2 = _mm_add_sd(c49_2, _mm_mul_sd(a49_2, b49));
#endif
_mm_store_sd(&C[(i*88)+49], c49_2);
__m128d c49_3 = _mm_load_sd(&C[(i*88)+77]);
__m128d a49_3 = _mm_load_sd(&values[217]);
#if defined(__SSE3__) && defined(__AVX__)
c49_3 = _mm_add_sd(c49_3, _mm_mul_sd(a49_3, _mm256_castpd256_pd128(b49)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c49_3 = _mm_add_sd(c49_3, _mm_mul_sd(a49_3, b49));
#endif
_mm_store_sd(&C[(i*88)+77], c49_3);
#else
C[(i*88)+13] += values[214] * B[(i*88)+49];
C[(i*88)+28] += values[215] * B[(i*88)+49];
C[(i*88)+49] += values[216] * B[(i*88)+49];
C[(i*88)+77] += values[217] * B[(i*88)+49];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b50 = _mm256_broadcast_sd(&B[(i*88)+50]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b50 = _mm_loaddup_pd(&B[(i*88)+50]);
#endif
__m128d c50_0 = _mm_load_sd(&C[(i*88)+4]);
__m128d a50_0 = _mm_load_sd(&values[218]);
#if defined(__SSE3__) && defined(__AVX__)
c50_0 = _mm_add_sd(c50_0, _mm_mul_sd(a50_0, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c50_0 = _mm_add_sd(c50_0, _mm_mul_sd(a50_0, b50));
#endif
_mm_store_sd(&C[(i*88)+4], c50_0);
__m128d c50_1 = _mm_load_sd(&C[(i*88)+14]);
__m128d a50_1 = _mm_load_sd(&values[219]);
#if defined(__SSE3__) && defined(__AVX__)
c50_1 = _mm_add_sd(c50_1, _mm_mul_sd(a50_1, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c50_1 = _mm_add_sd(c50_1, _mm_mul_sd(a50_1, b50));
#endif
_mm_store_sd(&C[(i*88)+14], c50_1);
__m128d c50_2 = _mm_load_sd(&C[(i*88)+29]);
__m128d a50_2 = _mm_load_sd(&values[220]);
#if defined(__SSE3__) && defined(__AVX__)
c50_2 = _mm_add_sd(c50_2, _mm_mul_sd(a50_2, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c50_2 = _mm_add_sd(c50_2, _mm_mul_sd(a50_2, b50));
#endif
_mm_store_sd(&C[(i*88)+29], c50_2);
__m128d c50_3 = _mm_load_sd(&C[(i*88)+50]);
__m128d a50_3 = _mm_load_sd(&values[221]);
#if defined(__SSE3__) && defined(__AVX__)
c50_3 = _mm_add_sd(c50_3, _mm_mul_sd(a50_3, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c50_3 = _mm_add_sd(c50_3, _mm_mul_sd(a50_3, b50));
#endif
_mm_store_sd(&C[(i*88)+50], c50_3);
__m128d c50_4 = _mm_load_sd(&C[(i*88)+78]);
__m128d a50_4 = _mm_load_sd(&values[222]);
#if defined(__SSE3__) && defined(__AVX__)
c50_4 = _mm_add_sd(c50_4, _mm_mul_sd(a50_4, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c50_4 = _mm_add_sd(c50_4, _mm_mul_sd(a50_4, b50));
#endif
_mm_store_sd(&C[(i*88)+78], c50_4);
#else
C[(i*88)+4] += values[218] * B[(i*88)+50];
C[(i*88)+14] += values[219] * B[(i*88)+50];
C[(i*88)+29] += values[220] * B[(i*88)+50];
C[(i*88)+50] += values[221] * B[(i*88)+50];
C[(i*88)+78] += values[222] * B[(i*88)+50];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b51 = _mm256_broadcast_sd(&B[(i*88)+51]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b51 = _mm_loaddup_pd(&B[(i*88)+51]);
#endif
__m128d c51_0 = _mm_load_sd(&C[(i*88)+5]);
__m128d a51_0 = _mm_load_sd(&values[223]);
#if defined(__SSE3__) && defined(__AVX__)
c51_0 = _mm_add_sd(c51_0, _mm_mul_sd(a51_0, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c51_0 = _mm_add_sd(c51_0, _mm_mul_sd(a51_0, b51));
#endif
_mm_store_sd(&C[(i*88)+5], c51_0);
__m128d c51_1 = _mm_load_sd(&C[(i*88)+15]);
__m128d a51_1 = _mm_load_sd(&values[224]);
#if defined(__SSE3__) && defined(__AVX__)
c51_1 = _mm_add_sd(c51_1, _mm_mul_sd(a51_1, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c51_1 = _mm_add_sd(c51_1, _mm_mul_sd(a51_1, b51));
#endif
_mm_store_sd(&C[(i*88)+15], c51_1);
__m128d c51_2 = _mm_load_sd(&C[(i*88)+30]);
__m128d a51_2 = _mm_load_sd(&values[225]);
#if defined(__SSE3__) && defined(__AVX__)
c51_2 = _mm_add_sd(c51_2, _mm_mul_sd(a51_2, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c51_2 = _mm_add_sd(c51_2, _mm_mul_sd(a51_2, b51));
#endif
_mm_store_sd(&C[(i*88)+30], c51_2);
__m128d c51_3 = _mm_load_sd(&C[(i*88)+51]);
__m128d a51_3 = _mm_load_sd(&values[226]);
#if defined(__SSE3__) && defined(__AVX__)
c51_3 = _mm_add_sd(c51_3, _mm_mul_sd(a51_3, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c51_3 = _mm_add_sd(c51_3, _mm_mul_sd(a51_3, b51));
#endif
_mm_store_sd(&C[(i*88)+51], c51_3);
__m128d c51_4 = _mm_load_sd(&C[(i*88)+79]);
__m128d a51_4 = _mm_load_sd(&values[227]);
#if defined(__SSE3__) && defined(__AVX__)
c51_4 = _mm_add_sd(c51_4, _mm_mul_sd(a51_4, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c51_4 = _mm_add_sd(c51_4, _mm_mul_sd(a51_4, b51));
#endif
_mm_store_sd(&C[(i*88)+79], c51_4);
#else
C[(i*88)+5] += values[223] * B[(i*88)+51];
C[(i*88)+15] += values[224] * B[(i*88)+51];
C[(i*88)+30] += values[225] * B[(i*88)+51];
C[(i*88)+51] += values[226] * B[(i*88)+51];
C[(i*88)+79] += values[227] * B[(i*88)+51];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b52 = _mm256_broadcast_sd(&B[(i*88)+52]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b52 = _mm_loaddup_pd(&B[(i*88)+52]);
#endif
__m128d c52_0 = _mm_load_sd(&C[(i*88)+6]);
__m128d a52_0 = _mm_load_sd(&values[228]);
#if defined(__SSE3__) && defined(__AVX__)
c52_0 = _mm_add_sd(c52_0, _mm_mul_sd(a52_0, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c52_0 = _mm_add_sd(c52_0, _mm_mul_sd(a52_0, b52));
#endif
_mm_store_sd(&C[(i*88)+6], c52_0);
__m128d c52_1 = _mm_load_sd(&C[(i*88)+16]);
__m128d a52_1 = _mm_load_sd(&values[229]);
#if defined(__SSE3__) && defined(__AVX__)
c52_1 = _mm_add_sd(c52_1, _mm_mul_sd(a52_1, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c52_1 = _mm_add_sd(c52_1, _mm_mul_sd(a52_1, b52));
#endif
_mm_store_sd(&C[(i*88)+16], c52_1);
__m128d c52_2 = _mm_load_sd(&C[(i*88)+31]);
__m128d a52_2 = _mm_load_sd(&values[230]);
#if defined(__SSE3__) && defined(__AVX__)
c52_2 = _mm_add_sd(c52_2, _mm_mul_sd(a52_2, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c52_2 = _mm_add_sd(c52_2, _mm_mul_sd(a52_2, b52));
#endif
_mm_store_sd(&C[(i*88)+31], c52_2);
__m128d c52_3 = _mm_load_sd(&C[(i*88)+52]);
__m128d a52_3 = _mm_load_sd(&values[231]);
#if defined(__SSE3__) && defined(__AVX__)
c52_3 = _mm_add_sd(c52_3, _mm_mul_sd(a52_3, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c52_3 = _mm_add_sd(c52_3, _mm_mul_sd(a52_3, b52));
#endif
_mm_store_sd(&C[(i*88)+52], c52_3);
__m128d c52_4 = _mm_load_sd(&C[(i*88)+80]);
__m128d a52_4 = _mm_load_sd(&values[232]);
#if defined(__SSE3__) && defined(__AVX__)
c52_4 = _mm_add_sd(c52_4, _mm_mul_sd(a52_4, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c52_4 = _mm_add_sd(c52_4, _mm_mul_sd(a52_4, b52));
#endif
_mm_store_sd(&C[(i*88)+80], c52_4);
#else
C[(i*88)+6] += values[228] * B[(i*88)+52];
C[(i*88)+16] += values[229] * B[(i*88)+52];
C[(i*88)+31] += values[230] * B[(i*88)+52];
C[(i*88)+52] += values[231] * B[(i*88)+52];
C[(i*88)+80] += values[232] * B[(i*88)+52];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b53 = _mm256_broadcast_sd(&B[(i*88)+53]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b53 = _mm_loaddup_pd(&B[(i*88)+53]);
#endif
__m128d c53_0 = _mm_load_sd(&C[(i*88)+1]);
__m128d a53_0 = _mm_load_sd(&values[233]);
#if defined(__SSE3__) && defined(__AVX__)
c53_0 = _mm_add_sd(c53_0, _mm_mul_sd(a53_0, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c53_0 = _mm_add_sd(c53_0, _mm_mul_sd(a53_0, b53));
#endif
_mm_store_sd(&C[(i*88)+1], c53_0);
__m128d c53_1 = _mm_load_sd(&C[(i*88)+7]);
__m128d a53_1 = _mm_load_sd(&values[234]);
#if defined(__SSE3__) && defined(__AVX__)
c53_1 = _mm_add_sd(c53_1, _mm_mul_sd(a53_1, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c53_1 = _mm_add_sd(c53_1, _mm_mul_sd(a53_1, b53));
#endif
_mm_store_sd(&C[(i*88)+7], c53_1);
__m128d c53_2 = _mm_load_sd(&C[(i*88)+17]);
__m128d a53_2 = _mm_load_sd(&values[235]);
#if defined(__SSE3__) && defined(__AVX__)
c53_2 = _mm_add_sd(c53_2, _mm_mul_sd(a53_2, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c53_2 = _mm_add_sd(c53_2, _mm_mul_sd(a53_2, b53));
#endif
_mm_store_sd(&C[(i*88)+17], c53_2);
__m128d c53_3 = _mm_load_sd(&C[(i*88)+32]);
__m128d a53_3 = _mm_load_sd(&values[236]);
#if defined(__SSE3__) && defined(__AVX__)
c53_3 = _mm_add_sd(c53_3, _mm_mul_sd(a53_3, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c53_3 = _mm_add_sd(c53_3, _mm_mul_sd(a53_3, b53));
#endif
_mm_store_sd(&C[(i*88)+32], c53_3);
__m128d c53_4 = _mm_load_sd(&C[(i*88)+53]);
__m128d a53_4 = _mm_load_sd(&values[237]);
#if defined(__SSE3__) && defined(__AVX__)
c53_4 = _mm_add_sd(c53_4, _mm_mul_sd(a53_4, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c53_4 = _mm_add_sd(c53_4, _mm_mul_sd(a53_4, b53));
#endif
_mm_store_sd(&C[(i*88)+53], c53_4);
__m128d c53_5 = _mm_load_sd(&C[(i*88)+81]);
__m128d a53_5 = _mm_load_sd(&values[238]);
#if defined(__SSE3__) && defined(__AVX__)
c53_5 = _mm_add_sd(c53_5, _mm_mul_sd(a53_5, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c53_5 = _mm_add_sd(c53_5, _mm_mul_sd(a53_5, b53));
#endif
_mm_store_sd(&C[(i*88)+81], c53_5);
#else
C[(i*88)+1] += values[233] * B[(i*88)+53];
C[(i*88)+7] += values[234] * B[(i*88)+53];
C[(i*88)+17] += values[235] * B[(i*88)+53];
C[(i*88)+32] += values[236] * B[(i*88)+53];
C[(i*88)+53] += values[237] * B[(i*88)+53];
C[(i*88)+81] += values[238] * B[(i*88)+53];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b54 = _mm256_broadcast_sd(&B[(i*88)+54]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b54 = _mm_loaddup_pd(&B[(i*88)+54]);
#endif
__m128d c54_0 = _mm_load_sd(&C[(i*88)+2]);
__m128d a54_0 = _mm_load_sd(&values[239]);
#if defined(__SSE3__) && defined(__AVX__)
c54_0 = _mm_add_sd(c54_0, _mm_mul_sd(a54_0, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c54_0 = _mm_add_sd(c54_0, _mm_mul_sd(a54_0, b54));
#endif
_mm_store_sd(&C[(i*88)+2], c54_0);
__m128d c54_1 = _mm_load_sd(&C[(i*88)+8]);
__m128d a54_1 = _mm_load_sd(&values[240]);
#if defined(__SSE3__) && defined(__AVX__)
c54_1 = _mm_add_sd(c54_1, _mm_mul_sd(a54_1, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c54_1 = _mm_add_sd(c54_1, _mm_mul_sd(a54_1, b54));
#endif
_mm_store_sd(&C[(i*88)+8], c54_1);
__m128d c54_2 = _mm_load_sd(&C[(i*88)+18]);
__m128d a54_2 = _mm_load_sd(&values[241]);
#if defined(__SSE3__) && defined(__AVX__)
c54_2 = _mm_add_sd(c54_2, _mm_mul_sd(a54_2, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c54_2 = _mm_add_sd(c54_2, _mm_mul_sd(a54_2, b54));
#endif
_mm_store_sd(&C[(i*88)+18], c54_2);
__m128d c54_3 = _mm_load_sd(&C[(i*88)+33]);
__m128d a54_3 = _mm_load_sd(&values[242]);
#if defined(__SSE3__) && defined(__AVX__)
c54_3 = _mm_add_sd(c54_3, _mm_mul_sd(a54_3, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c54_3 = _mm_add_sd(c54_3, _mm_mul_sd(a54_3, b54));
#endif
_mm_store_sd(&C[(i*88)+33], c54_3);
__m128d c54_4 = _mm_load_sd(&C[(i*88)+54]);
__m128d a54_4 = _mm_load_sd(&values[243]);
#if defined(__SSE3__) && defined(__AVX__)
c54_4 = _mm_add_sd(c54_4, _mm_mul_sd(a54_4, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c54_4 = _mm_add_sd(c54_4, _mm_mul_sd(a54_4, b54));
#endif
_mm_store_sd(&C[(i*88)+54], c54_4);
__m128d c54_5 = _mm_load_sd(&C[(i*88)+82]);
__m128d a54_5 = _mm_load_sd(&values[244]);
#if defined(__SSE3__) && defined(__AVX__)
c54_5 = _mm_add_sd(c54_5, _mm_mul_sd(a54_5, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c54_5 = _mm_add_sd(c54_5, _mm_mul_sd(a54_5, b54));
#endif
_mm_store_sd(&C[(i*88)+82], c54_5);
#else
C[(i*88)+2] += values[239] * B[(i*88)+54];
C[(i*88)+8] += values[240] * B[(i*88)+54];
C[(i*88)+18] += values[241] * B[(i*88)+54];
C[(i*88)+33] += values[242] * B[(i*88)+54];
C[(i*88)+54] += values[243] * B[(i*88)+54];
C[(i*88)+82] += values[244] * B[(i*88)+54];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b55 = _mm256_broadcast_sd(&B[(i*88)+55]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b55 = _mm_loaddup_pd(&B[(i*88)+55]);
#endif
__m128d c55_0 = _mm_load_sd(&C[(i*88)+0]);
__m128d a55_0 = _mm_load_sd(&values[245]);
#if defined(__SSE3__) && defined(__AVX__)
c55_0 = _mm_add_sd(c55_0, _mm_mul_sd(a55_0, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_0 = _mm_add_sd(c55_0, _mm_mul_sd(a55_0, b55));
#endif
_mm_store_sd(&C[(i*88)+0], c55_0);
__m128d c55_1 = _mm_load_sd(&C[(i*88)+3]);
__m128d a55_1 = _mm_load_sd(&values[246]);
#if defined(__SSE3__) && defined(__AVX__)
c55_1 = _mm_add_sd(c55_1, _mm_mul_sd(a55_1, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_1 = _mm_add_sd(c55_1, _mm_mul_sd(a55_1, b55));
#endif
_mm_store_sd(&C[(i*88)+3], c55_1);
__m128d c55_2 = _mm_load_sd(&C[(i*88)+9]);
__m128d a55_2 = _mm_load_sd(&values[247]);
#if defined(__SSE3__) && defined(__AVX__)
c55_2 = _mm_add_sd(c55_2, _mm_mul_sd(a55_2, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_2 = _mm_add_sd(c55_2, _mm_mul_sd(a55_2, b55));
#endif
_mm_store_sd(&C[(i*88)+9], c55_2);
__m128d c55_3 = _mm_load_sd(&C[(i*88)+19]);
__m128d a55_3 = _mm_load_sd(&values[248]);
#if defined(__SSE3__) && defined(__AVX__)
c55_3 = _mm_add_sd(c55_3, _mm_mul_sd(a55_3, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_3 = _mm_add_sd(c55_3, _mm_mul_sd(a55_3, b55));
#endif
_mm_store_sd(&C[(i*88)+19], c55_3);
__m128d c55_4 = _mm_load_sd(&C[(i*88)+34]);
__m128d a55_4 = _mm_load_sd(&values[249]);
#if defined(__SSE3__) && defined(__AVX__)
c55_4 = _mm_add_sd(c55_4, _mm_mul_sd(a55_4, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_4 = _mm_add_sd(c55_4, _mm_mul_sd(a55_4, b55));
#endif
_mm_store_sd(&C[(i*88)+34], c55_4);
__m128d c55_5 = _mm_load_sd(&C[(i*88)+55]);
__m128d a55_5 = _mm_load_sd(&values[250]);
#if defined(__SSE3__) && defined(__AVX__)
c55_5 = _mm_add_sd(c55_5, _mm_mul_sd(a55_5, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_5 = _mm_add_sd(c55_5, _mm_mul_sd(a55_5, b55));
#endif
_mm_store_sd(&C[(i*88)+55], c55_5);
__m128d c55_6 = _mm_load_sd(&C[(i*88)+83]);
__m128d a55_6 = _mm_load_sd(&values[251]);
#if defined(__SSE3__) && defined(__AVX__)
c55_6 = _mm_add_sd(c55_6, _mm_mul_sd(a55_6, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_6 = _mm_add_sd(c55_6, _mm_mul_sd(a55_6, b55));
#endif
_mm_store_sd(&C[(i*88)+83], c55_6);
#else
C[(i*88)+0] += values[245] * B[(i*88)+55];
C[(i*88)+3] += values[246] * B[(i*88)+55];
C[(i*88)+9] += values[247] * B[(i*88)+55];
C[(i*88)+19] += values[248] * B[(i*88)+55];
C[(i*88)+34] += values[249] * B[(i*88)+55];
C[(i*88)+55] += values[250] * B[(i*88)+55];
C[(i*88)+83] += values[251] * B[(i*88)+55];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b56 = _mm256_broadcast_sd(&B[(i*88)+56]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b56 = _mm_loaddup_pd(&B[(i*88)+56]);
#endif
__m128d c56_0 = _mm_load_sd(&C[(i*88)+56]);
__m128d a56_0 = _mm_load_sd(&values[252]);
#if defined(__SSE3__) && defined(__AVX__)
c56_0 = _mm_add_sd(c56_0, _mm_mul_sd(a56_0, _mm256_castpd256_pd128(b56)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c56_0 = _mm_add_sd(c56_0, _mm_mul_sd(a56_0, b56));
#endif
_mm_store_sd(&C[(i*88)+56], c56_0);
#else
C[(i*88)+56] += values[252] * B[(i*88)+56];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b57 = _mm256_broadcast_sd(&B[(i*88)+57]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b57 = _mm_loaddup_pd(&B[(i*88)+57]);
#endif
__m128d c57_0 = _mm_load_sd(&C[(i*88)+57]);
__m128d a57_0 = _mm_load_sd(&values[253]);
#if defined(__SSE3__) && defined(__AVX__)
c57_0 = _mm_add_sd(c57_0, _mm_mul_sd(a57_0, _mm256_castpd256_pd128(b57)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c57_0 = _mm_add_sd(c57_0, _mm_mul_sd(a57_0, b57));
#endif
_mm_store_sd(&C[(i*88)+57], c57_0);
#else
C[(i*88)+57] += values[253] * B[(i*88)+57];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b58 = _mm256_broadcast_sd(&B[(i*88)+58]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b58 = _mm_loaddup_pd(&B[(i*88)+58]);
#endif
__m128d c58_0 = _mm_load_sd(&C[(i*88)+58]);
__m128d a58_0 = _mm_load_sd(&values[254]);
#if defined(__SSE3__) && defined(__AVX__)
c58_0 = _mm_add_sd(c58_0, _mm_mul_sd(a58_0, _mm256_castpd256_pd128(b58)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c58_0 = _mm_add_sd(c58_0, _mm_mul_sd(a58_0, b58));
#endif
_mm_store_sd(&C[(i*88)+58], c58_0);
#else
C[(i*88)+58] += values[254] * B[(i*88)+58];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b59 = _mm256_broadcast_sd(&B[(i*88)+59]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b59 = _mm_loaddup_pd(&B[(i*88)+59]);
#endif
__m128d c59_0 = _mm_load_sd(&C[(i*88)+59]);
__m128d a59_0 = _mm_load_sd(&values[255]);
#if defined(__SSE3__) && defined(__AVX__)
c59_0 = _mm_add_sd(c59_0, _mm_mul_sd(a59_0, _mm256_castpd256_pd128(b59)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c59_0 = _mm_add_sd(c59_0, _mm_mul_sd(a59_0, b59));
#endif
_mm_store_sd(&C[(i*88)+59], c59_0);
#else
C[(i*88)+59] += values[255] * B[(i*88)+59];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b60 = _mm256_broadcast_sd(&B[(i*88)+60]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b60 = _mm_loaddup_pd(&B[(i*88)+60]);
#endif
__m128d c60_0 = _mm_load_sd(&C[(i*88)+60]);
__m128d a60_0 = _mm_load_sd(&values[256]);
#if defined(__SSE3__) && defined(__AVX__)
c60_0 = _mm_add_sd(c60_0, _mm_mul_sd(a60_0, _mm256_castpd256_pd128(b60)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c60_0 = _mm_add_sd(c60_0, _mm_mul_sd(a60_0, b60));
#endif
_mm_store_sd(&C[(i*88)+60], c60_0);
#else
C[(i*88)+60] += values[256] * B[(i*88)+60];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b61 = _mm256_broadcast_sd(&B[(i*88)+61]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b61 = _mm_loaddup_pd(&B[(i*88)+61]);
#endif
__m128d c61_0 = _mm_load_sd(&C[(i*88)+61]);
__m128d a61_0 = _mm_load_sd(&values[257]);
#if defined(__SSE3__) && defined(__AVX__)
c61_0 = _mm_add_sd(c61_0, _mm_mul_sd(a61_0, _mm256_castpd256_pd128(b61)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c61_0 = _mm_add_sd(c61_0, _mm_mul_sd(a61_0, b61));
#endif
_mm_store_sd(&C[(i*88)+61], c61_0);
#else
C[(i*88)+61] += values[257] * B[(i*88)+61];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b62 = _mm256_broadcast_sd(&B[(i*88)+62]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b62 = _mm_loaddup_pd(&B[(i*88)+62]);
#endif
__m128d c62_0 = _mm_load_sd(&C[(i*88)+62]);
__m128d a62_0 = _mm_load_sd(&values[258]);
#if defined(__SSE3__) && defined(__AVX__)
c62_0 = _mm_add_sd(c62_0, _mm_mul_sd(a62_0, _mm256_castpd256_pd128(b62)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c62_0 = _mm_add_sd(c62_0, _mm_mul_sd(a62_0, b62));
#endif
_mm_store_sd(&C[(i*88)+62], c62_0);
#else
C[(i*88)+62] += values[258] * B[(i*88)+62];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b63 = _mm256_broadcast_sd(&B[(i*88)+63]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b63 = _mm_loaddup_pd(&B[(i*88)+63]);
#endif
__m128d c63_0 = _mm_load_sd(&C[(i*88)+35]);
__m128d a63_0 = _mm_load_sd(&values[259]);
#if defined(__SSE3__) && defined(__AVX__)
c63_0 = _mm_add_sd(c63_0, _mm_mul_sd(a63_0, _mm256_castpd256_pd128(b63)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c63_0 = _mm_add_sd(c63_0, _mm_mul_sd(a63_0, b63));
#endif
_mm_store_sd(&C[(i*88)+35], c63_0);
__m128d c63_1 = _mm_load_sd(&C[(i*88)+63]);
__m128d a63_1 = _mm_load_sd(&values[260]);
#if defined(__SSE3__) && defined(__AVX__)
c63_1 = _mm_add_sd(c63_1, _mm_mul_sd(a63_1, _mm256_castpd256_pd128(b63)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c63_1 = _mm_add_sd(c63_1, _mm_mul_sd(a63_1, b63));
#endif
_mm_store_sd(&C[(i*88)+63], c63_1);
#else
C[(i*88)+35] += values[259] * B[(i*88)+63];
C[(i*88)+63] += values[260] * B[(i*88)+63];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b64 = _mm256_broadcast_sd(&B[(i*88)+64]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b64 = _mm_loaddup_pd(&B[(i*88)+64]);
#endif
__m128d c64_0 = _mm_load_sd(&C[(i*88)+36]);
__m128d a64_0 = _mm_load_sd(&values[261]);
#if defined(__SSE3__) && defined(__AVX__)
c64_0 = _mm_add_sd(c64_0, _mm_mul_sd(a64_0, _mm256_castpd256_pd128(b64)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c64_0 = _mm_add_sd(c64_0, _mm_mul_sd(a64_0, b64));
#endif
_mm_store_sd(&C[(i*88)+36], c64_0);
__m128d c64_1 = _mm_load_sd(&C[(i*88)+64]);
__m128d a64_1 = _mm_load_sd(&values[262]);
#if defined(__SSE3__) && defined(__AVX__)
c64_1 = _mm_add_sd(c64_1, _mm_mul_sd(a64_1, _mm256_castpd256_pd128(b64)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c64_1 = _mm_add_sd(c64_1, _mm_mul_sd(a64_1, b64));
#endif
_mm_store_sd(&C[(i*88)+64], c64_1);
#else
C[(i*88)+36] += values[261] * B[(i*88)+64];
C[(i*88)+64] += values[262] * B[(i*88)+64];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b65 = _mm256_broadcast_sd(&B[(i*88)+65]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b65 = _mm_loaddup_pd(&B[(i*88)+65]);
#endif
__m128d c65_0 = _mm_load_sd(&C[(i*88)+37]);
__m128d a65_0 = _mm_load_sd(&values[263]);
#if defined(__SSE3__) && defined(__AVX__)
c65_0 = _mm_add_sd(c65_0, _mm_mul_sd(a65_0, _mm256_castpd256_pd128(b65)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c65_0 = _mm_add_sd(c65_0, _mm_mul_sd(a65_0, b65));
#endif
_mm_store_sd(&C[(i*88)+37], c65_0);
__m128d c65_1 = _mm_load_sd(&C[(i*88)+65]);
__m128d a65_1 = _mm_load_sd(&values[264]);
#if defined(__SSE3__) && defined(__AVX__)
c65_1 = _mm_add_sd(c65_1, _mm_mul_sd(a65_1, _mm256_castpd256_pd128(b65)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c65_1 = _mm_add_sd(c65_1, _mm_mul_sd(a65_1, b65));
#endif
_mm_store_sd(&C[(i*88)+65], c65_1);
#else
C[(i*88)+37] += values[263] * B[(i*88)+65];
C[(i*88)+65] += values[264] * B[(i*88)+65];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b66 = _mm256_broadcast_sd(&B[(i*88)+66]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b66 = _mm_loaddup_pd(&B[(i*88)+66]);
#endif
__m128d c66_0 = _mm_load_sd(&C[(i*88)+38]);
__m128d a66_0 = _mm_load_sd(&values[265]);
#if defined(__SSE3__) && defined(__AVX__)
c66_0 = _mm_add_sd(c66_0, _mm_mul_sd(a66_0, _mm256_castpd256_pd128(b66)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c66_0 = _mm_add_sd(c66_0, _mm_mul_sd(a66_0, b66));
#endif
_mm_store_sd(&C[(i*88)+38], c66_0);
__m128d c66_1 = _mm_load_sd(&C[(i*88)+66]);
__m128d a66_1 = _mm_load_sd(&values[266]);
#if defined(__SSE3__) && defined(__AVX__)
c66_1 = _mm_add_sd(c66_1, _mm_mul_sd(a66_1, _mm256_castpd256_pd128(b66)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c66_1 = _mm_add_sd(c66_1, _mm_mul_sd(a66_1, b66));
#endif
_mm_store_sd(&C[(i*88)+66], c66_1);
#else
C[(i*88)+38] += values[265] * B[(i*88)+66];
C[(i*88)+66] += values[266] * B[(i*88)+66];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b67 = _mm256_broadcast_sd(&B[(i*88)+67]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b67 = _mm_loaddup_pd(&B[(i*88)+67]);
#endif
__m128d c67_0 = _mm_load_sd(&C[(i*88)+39]);
__m128d a67_0 = _mm_load_sd(&values[267]);
#if defined(__SSE3__) && defined(__AVX__)
c67_0 = _mm_add_sd(c67_0, _mm_mul_sd(a67_0, _mm256_castpd256_pd128(b67)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c67_0 = _mm_add_sd(c67_0, _mm_mul_sd(a67_0, b67));
#endif
_mm_store_sd(&C[(i*88)+39], c67_0);
__m128d c67_1 = _mm_load_sd(&C[(i*88)+67]);
__m128d a67_1 = _mm_load_sd(&values[268]);
#if defined(__SSE3__) && defined(__AVX__)
c67_1 = _mm_add_sd(c67_1, _mm_mul_sd(a67_1, _mm256_castpd256_pd128(b67)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c67_1 = _mm_add_sd(c67_1, _mm_mul_sd(a67_1, b67));
#endif
_mm_store_sd(&C[(i*88)+67], c67_1);
#else
C[(i*88)+39] += values[267] * B[(i*88)+67];
C[(i*88)+67] += values[268] * B[(i*88)+67];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b68 = _mm256_broadcast_sd(&B[(i*88)+68]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b68 = _mm_loaddup_pd(&B[(i*88)+68]);
#endif
__m128d c68_0 = _mm_load_sd(&C[(i*88)+40]);
__m128d a68_0 = _mm_load_sd(&values[269]);
#if defined(__SSE3__) && defined(__AVX__)
c68_0 = _mm_add_sd(c68_0, _mm_mul_sd(a68_0, _mm256_castpd256_pd128(b68)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c68_0 = _mm_add_sd(c68_0, _mm_mul_sd(a68_0, b68));
#endif
_mm_store_sd(&C[(i*88)+40], c68_0);
__m128d c68_1 = _mm_load_sd(&C[(i*88)+68]);
__m128d a68_1 = _mm_load_sd(&values[270]);
#if defined(__SSE3__) && defined(__AVX__)
c68_1 = _mm_add_sd(c68_1, _mm_mul_sd(a68_1, _mm256_castpd256_pd128(b68)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c68_1 = _mm_add_sd(c68_1, _mm_mul_sd(a68_1, b68));
#endif
_mm_store_sd(&C[(i*88)+68], c68_1);
#else
C[(i*88)+40] += values[269] * B[(i*88)+68];
C[(i*88)+68] += values[270] * B[(i*88)+68];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b69 = _mm256_broadcast_sd(&B[(i*88)+69]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b69 = _mm_loaddup_pd(&B[(i*88)+69]);
#endif
__m128d c69_0 = _mm_load_sd(&C[(i*88)+20]);
__m128d a69_0 = _mm_load_sd(&values[271]);
#if defined(__SSE3__) && defined(__AVX__)
c69_0 = _mm_add_sd(c69_0, _mm_mul_sd(a69_0, _mm256_castpd256_pd128(b69)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c69_0 = _mm_add_sd(c69_0, _mm_mul_sd(a69_0, b69));
#endif
_mm_store_sd(&C[(i*88)+20], c69_0);
__m128d c69_1 = _mm_load_sd(&C[(i*88)+41]);
__m128d a69_1 = _mm_load_sd(&values[272]);
#if defined(__SSE3__) && defined(__AVX__)
c69_1 = _mm_add_sd(c69_1, _mm_mul_sd(a69_1, _mm256_castpd256_pd128(b69)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c69_1 = _mm_add_sd(c69_1, _mm_mul_sd(a69_1, b69));
#endif
_mm_store_sd(&C[(i*88)+41], c69_1);
__m128d c69_2 = _mm_load_sd(&C[(i*88)+69]);
__m128d a69_2 = _mm_load_sd(&values[273]);
#if defined(__SSE3__) && defined(__AVX__)
c69_2 = _mm_add_sd(c69_2, _mm_mul_sd(a69_2, _mm256_castpd256_pd128(b69)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c69_2 = _mm_add_sd(c69_2, _mm_mul_sd(a69_2, b69));
#endif
_mm_store_sd(&C[(i*88)+69], c69_2);
#else
C[(i*88)+20] += values[271] * B[(i*88)+69];
C[(i*88)+41] += values[272] * B[(i*88)+69];
C[(i*88)+69] += values[273] * B[(i*88)+69];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b70 = _mm256_broadcast_sd(&B[(i*88)+70]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b70 = _mm_loaddup_pd(&B[(i*88)+70]);
#endif
__m128d c70_0 = _mm_load_sd(&C[(i*88)+21]);
__m128d a70_0 = _mm_load_sd(&values[274]);
#if defined(__SSE3__) && defined(__AVX__)
c70_0 = _mm_add_sd(c70_0, _mm_mul_sd(a70_0, _mm256_castpd256_pd128(b70)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c70_0 = _mm_add_sd(c70_0, _mm_mul_sd(a70_0, b70));
#endif
_mm_store_sd(&C[(i*88)+21], c70_0);
__m128d c70_1 = _mm_load_sd(&C[(i*88)+42]);
__m128d a70_1 = _mm_load_sd(&values[275]);
#if defined(__SSE3__) && defined(__AVX__)
c70_1 = _mm_add_sd(c70_1, _mm_mul_sd(a70_1, _mm256_castpd256_pd128(b70)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c70_1 = _mm_add_sd(c70_1, _mm_mul_sd(a70_1, b70));
#endif
_mm_store_sd(&C[(i*88)+42], c70_1);
__m128d c70_2 = _mm_load_sd(&C[(i*88)+70]);
__m128d a70_2 = _mm_load_sd(&values[276]);
#if defined(__SSE3__) && defined(__AVX__)
c70_2 = _mm_add_sd(c70_2, _mm_mul_sd(a70_2, _mm256_castpd256_pd128(b70)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c70_2 = _mm_add_sd(c70_2, _mm_mul_sd(a70_2, b70));
#endif
_mm_store_sd(&C[(i*88)+70], c70_2);
#else
C[(i*88)+21] += values[274] * B[(i*88)+70];
C[(i*88)+42] += values[275] * B[(i*88)+70];
C[(i*88)+70] += values[276] * B[(i*88)+70];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b71 = _mm256_broadcast_sd(&B[(i*88)+71]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b71 = _mm_loaddup_pd(&B[(i*88)+71]);
#endif
__m128d c71_0 = _mm_load_sd(&C[(i*88)+22]);
__m128d a71_0 = _mm_load_sd(&values[277]);
#if defined(__SSE3__) && defined(__AVX__)
c71_0 = _mm_add_sd(c71_0, _mm_mul_sd(a71_0, _mm256_castpd256_pd128(b71)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c71_0 = _mm_add_sd(c71_0, _mm_mul_sd(a71_0, b71));
#endif
_mm_store_sd(&C[(i*88)+22], c71_0);
__m128d c71_1 = _mm_load_sd(&C[(i*88)+43]);
__m128d a71_1 = _mm_load_sd(&values[278]);
#if defined(__SSE3__) && defined(__AVX__)
c71_1 = _mm_add_sd(c71_1, _mm_mul_sd(a71_1, _mm256_castpd256_pd128(b71)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c71_1 = _mm_add_sd(c71_1, _mm_mul_sd(a71_1, b71));
#endif
_mm_store_sd(&C[(i*88)+43], c71_1);
__m128d c71_2 = _mm_load_sd(&C[(i*88)+71]);
__m128d a71_2 = _mm_load_sd(&values[279]);
#if defined(__SSE3__) && defined(__AVX__)
c71_2 = _mm_add_sd(c71_2, _mm_mul_sd(a71_2, _mm256_castpd256_pd128(b71)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c71_2 = _mm_add_sd(c71_2, _mm_mul_sd(a71_2, b71));
#endif
_mm_store_sd(&C[(i*88)+71], c71_2);
#else
C[(i*88)+22] += values[277] * B[(i*88)+71];
C[(i*88)+43] += values[278] * B[(i*88)+71];
C[(i*88)+71] += values[279] * B[(i*88)+71];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b72 = _mm256_broadcast_sd(&B[(i*88)+72]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b72 = _mm_loaddup_pd(&B[(i*88)+72]);
#endif
__m128d c72_0 = _mm_load_sd(&C[(i*88)+23]);
__m128d a72_0 = _mm_load_sd(&values[280]);
#if defined(__SSE3__) && defined(__AVX__)
c72_0 = _mm_add_sd(c72_0, _mm_mul_sd(a72_0, _mm256_castpd256_pd128(b72)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c72_0 = _mm_add_sd(c72_0, _mm_mul_sd(a72_0, b72));
#endif
_mm_store_sd(&C[(i*88)+23], c72_0);
__m128d c72_1 = _mm_load_sd(&C[(i*88)+44]);
__m128d a72_1 = _mm_load_sd(&values[281]);
#if defined(__SSE3__) && defined(__AVX__)
c72_1 = _mm_add_sd(c72_1, _mm_mul_sd(a72_1, _mm256_castpd256_pd128(b72)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c72_1 = _mm_add_sd(c72_1, _mm_mul_sd(a72_1, b72));
#endif
_mm_store_sd(&C[(i*88)+44], c72_1);
__m128d c72_2 = _mm_load_sd(&C[(i*88)+72]);
__m128d a72_2 = _mm_load_sd(&values[282]);
#if defined(__SSE3__) && defined(__AVX__)
c72_2 = _mm_add_sd(c72_2, _mm_mul_sd(a72_2, _mm256_castpd256_pd128(b72)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c72_2 = _mm_add_sd(c72_2, _mm_mul_sd(a72_2, b72));
#endif
_mm_store_sd(&C[(i*88)+72], c72_2);
#else
C[(i*88)+23] += values[280] * B[(i*88)+72];
C[(i*88)+44] += values[281] * B[(i*88)+72];
C[(i*88)+72] += values[282] * B[(i*88)+72];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b73 = _mm256_broadcast_sd(&B[(i*88)+73]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b73 = _mm_loaddup_pd(&B[(i*88)+73]);
#endif
__m128d c73_0 = _mm_load_sd(&C[(i*88)+24]);
__m128d a73_0 = _mm_load_sd(&values[283]);
#if defined(__SSE3__) && defined(__AVX__)
c73_0 = _mm_add_sd(c73_0, _mm_mul_sd(a73_0, _mm256_castpd256_pd128(b73)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c73_0 = _mm_add_sd(c73_0, _mm_mul_sd(a73_0, b73));
#endif
_mm_store_sd(&C[(i*88)+24], c73_0);
__m128d c73_1 = _mm_load_sd(&C[(i*88)+45]);
__m128d a73_1 = _mm_load_sd(&values[284]);
#if defined(__SSE3__) && defined(__AVX__)
c73_1 = _mm_add_sd(c73_1, _mm_mul_sd(a73_1, _mm256_castpd256_pd128(b73)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c73_1 = _mm_add_sd(c73_1, _mm_mul_sd(a73_1, b73));
#endif
_mm_store_sd(&C[(i*88)+45], c73_1);
__m128d c73_2 = _mm_load_sd(&C[(i*88)+73]);
__m128d a73_2 = _mm_load_sd(&values[285]);
#if defined(__SSE3__) && defined(__AVX__)
c73_2 = _mm_add_sd(c73_2, _mm_mul_sd(a73_2, _mm256_castpd256_pd128(b73)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c73_2 = _mm_add_sd(c73_2, _mm_mul_sd(a73_2, b73));
#endif
_mm_store_sd(&C[(i*88)+73], c73_2);
#else
C[(i*88)+24] += values[283] * B[(i*88)+73];
C[(i*88)+45] += values[284] * B[(i*88)+73];
C[(i*88)+73] += values[285] * B[(i*88)+73];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b74 = _mm256_broadcast_sd(&B[(i*88)+74]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b74 = _mm_loaddup_pd(&B[(i*88)+74]);
#endif
__m128d c74_0 = _mm_load_sd(&C[(i*88)+10]);
__m128d a74_0 = _mm_load_sd(&values[286]);
#if defined(__SSE3__) && defined(__AVX__)
c74_0 = _mm_add_sd(c74_0, _mm_mul_sd(a74_0, _mm256_castpd256_pd128(b74)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c74_0 = _mm_add_sd(c74_0, _mm_mul_sd(a74_0, b74));
#endif
_mm_store_sd(&C[(i*88)+10], c74_0);
__m128d c74_1 = _mm_load_sd(&C[(i*88)+25]);
__m128d a74_1 = _mm_load_sd(&values[287]);
#if defined(__SSE3__) && defined(__AVX__)
c74_1 = _mm_add_sd(c74_1, _mm_mul_sd(a74_1, _mm256_castpd256_pd128(b74)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c74_1 = _mm_add_sd(c74_1, _mm_mul_sd(a74_1, b74));
#endif
_mm_store_sd(&C[(i*88)+25], c74_1);
__m128d c74_2 = _mm_load_sd(&C[(i*88)+46]);
__m128d a74_2 = _mm_load_sd(&values[288]);
#if defined(__SSE3__) && defined(__AVX__)
c74_2 = _mm_add_sd(c74_2, _mm_mul_sd(a74_2, _mm256_castpd256_pd128(b74)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c74_2 = _mm_add_sd(c74_2, _mm_mul_sd(a74_2, b74));
#endif
_mm_store_sd(&C[(i*88)+46], c74_2);
__m128d c74_3 = _mm_load_sd(&C[(i*88)+74]);
__m128d a74_3 = _mm_load_sd(&values[289]);
#if defined(__SSE3__) && defined(__AVX__)
c74_3 = _mm_add_sd(c74_3, _mm_mul_sd(a74_3, _mm256_castpd256_pd128(b74)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c74_3 = _mm_add_sd(c74_3, _mm_mul_sd(a74_3, b74));
#endif
_mm_store_sd(&C[(i*88)+74], c74_3);
#else
C[(i*88)+10] += values[286] * B[(i*88)+74];
C[(i*88)+25] += values[287] * B[(i*88)+74];
C[(i*88)+46] += values[288] * B[(i*88)+74];
C[(i*88)+74] += values[289] * B[(i*88)+74];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b75 = _mm256_broadcast_sd(&B[(i*88)+75]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b75 = _mm_loaddup_pd(&B[(i*88)+75]);
#endif
__m128d c75_0 = _mm_load_sd(&C[(i*88)+11]);
__m128d a75_0 = _mm_load_sd(&values[290]);
#if defined(__SSE3__) && defined(__AVX__)
c75_0 = _mm_add_sd(c75_0, _mm_mul_sd(a75_0, _mm256_castpd256_pd128(b75)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c75_0 = _mm_add_sd(c75_0, _mm_mul_sd(a75_0, b75));
#endif
_mm_store_sd(&C[(i*88)+11], c75_0);
__m128d c75_1 = _mm_load_sd(&C[(i*88)+26]);
__m128d a75_1 = _mm_load_sd(&values[291]);
#if defined(__SSE3__) && defined(__AVX__)
c75_1 = _mm_add_sd(c75_1, _mm_mul_sd(a75_1, _mm256_castpd256_pd128(b75)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c75_1 = _mm_add_sd(c75_1, _mm_mul_sd(a75_1, b75));
#endif
_mm_store_sd(&C[(i*88)+26], c75_1);
__m128d c75_2 = _mm_load_sd(&C[(i*88)+47]);
__m128d a75_2 = _mm_load_sd(&values[292]);
#if defined(__SSE3__) && defined(__AVX__)
c75_2 = _mm_add_sd(c75_2, _mm_mul_sd(a75_2, _mm256_castpd256_pd128(b75)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c75_2 = _mm_add_sd(c75_2, _mm_mul_sd(a75_2, b75));
#endif
_mm_store_sd(&C[(i*88)+47], c75_2);
__m128d c75_3 = _mm_load_sd(&C[(i*88)+75]);
__m128d a75_3 = _mm_load_sd(&values[293]);
#if defined(__SSE3__) && defined(__AVX__)
c75_3 = _mm_add_sd(c75_3, _mm_mul_sd(a75_3, _mm256_castpd256_pd128(b75)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c75_3 = _mm_add_sd(c75_3, _mm_mul_sd(a75_3, b75));
#endif
_mm_store_sd(&C[(i*88)+75], c75_3);
#else
C[(i*88)+11] += values[290] * B[(i*88)+75];
C[(i*88)+26] += values[291] * B[(i*88)+75];
C[(i*88)+47] += values[292] * B[(i*88)+75];
C[(i*88)+75] += values[293] * B[(i*88)+75];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b76 = _mm256_broadcast_sd(&B[(i*88)+76]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b76 = _mm_loaddup_pd(&B[(i*88)+76]);
#endif
__m128d c76_0 = _mm_load_sd(&C[(i*88)+12]);
__m128d a76_0 = _mm_load_sd(&values[294]);
#if defined(__SSE3__) && defined(__AVX__)
c76_0 = _mm_add_sd(c76_0, _mm_mul_sd(a76_0, _mm256_castpd256_pd128(b76)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c76_0 = _mm_add_sd(c76_0, _mm_mul_sd(a76_0, b76));
#endif
_mm_store_sd(&C[(i*88)+12], c76_0);
__m128d c76_1 = _mm_load_sd(&C[(i*88)+27]);
__m128d a76_1 = _mm_load_sd(&values[295]);
#if defined(__SSE3__) && defined(__AVX__)
c76_1 = _mm_add_sd(c76_1, _mm_mul_sd(a76_1, _mm256_castpd256_pd128(b76)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c76_1 = _mm_add_sd(c76_1, _mm_mul_sd(a76_1, b76));
#endif
_mm_store_sd(&C[(i*88)+27], c76_1);
__m128d c76_2 = _mm_load_sd(&C[(i*88)+48]);
__m128d a76_2 = _mm_load_sd(&values[296]);
#if defined(__SSE3__) && defined(__AVX__)
c76_2 = _mm_add_sd(c76_2, _mm_mul_sd(a76_2, _mm256_castpd256_pd128(b76)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c76_2 = _mm_add_sd(c76_2, _mm_mul_sd(a76_2, b76));
#endif
_mm_store_sd(&C[(i*88)+48], c76_2);
__m128d c76_3 = _mm_load_sd(&C[(i*88)+76]);
__m128d a76_3 = _mm_load_sd(&values[297]);
#if defined(__SSE3__) && defined(__AVX__)
c76_3 = _mm_add_sd(c76_3, _mm_mul_sd(a76_3, _mm256_castpd256_pd128(b76)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c76_3 = _mm_add_sd(c76_3, _mm_mul_sd(a76_3, b76));
#endif
_mm_store_sd(&C[(i*88)+76], c76_3);
#else
C[(i*88)+12] += values[294] * B[(i*88)+76];
C[(i*88)+27] += values[295] * B[(i*88)+76];
C[(i*88)+48] += values[296] * B[(i*88)+76];
C[(i*88)+76] += values[297] * B[(i*88)+76];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b77 = _mm256_broadcast_sd(&B[(i*88)+77]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b77 = _mm_loaddup_pd(&B[(i*88)+77]);
#endif
__m128d c77_0 = _mm_load_sd(&C[(i*88)+13]);
__m128d a77_0 = _mm_load_sd(&values[298]);
#if defined(__SSE3__) && defined(__AVX__)
c77_0 = _mm_add_sd(c77_0, _mm_mul_sd(a77_0, _mm256_castpd256_pd128(b77)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c77_0 = _mm_add_sd(c77_0, _mm_mul_sd(a77_0, b77));
#endif
_mm_store_sd(&C[(i*88)+13], c77_0);
__m128d c77_1 = _mm_load_sd(&C[(i*88)+28]);
__m128d a77_1 = _mm_load_sd(&values[299]);
#if defined(__SSE3__) && defined(__AVX__)
c77_1 = _mm_add_sd(c77_1, _mm_mul_sd(a77_1, _mm256_castpd256_pd128(b77)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c77_1 = _mm_add_sd(c77_1, _mm_mul_sd(a77_1, b77));
#endif
_mm_store_sd(&C[(i*88)+28], c77_1);
__m128d c77_2 = _mm_load_sd(&C[(i*88)+49]);
__m128d a77_2 = _mm_load_sd(&values[300]);
#if defined(__SSE3__) && defined(__AVX__)
c77_2 = _mm_add_sd(c77_2, _mm_mul_sd(a77_2, _mm256_castpd256_pd128(b77)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c77_2 = _mm_add_sd(c77_2, _mm_mul_sd(a77_2, b77));
#endif
_mm_store_sd(&C[(i*88)+49], c77_2);
__m128d c77_3 = _mm_load_sd(&C[(i*88)+77]);
__m128d a77_3 = _mm_load_sd(&values[301]);
#if defined(__SSE3__) && defined(__AVX__)
c77_3 = _mm_add_sd(c77_3, _mm_mul_sd(a77_3, _mm256_castpd256_pd128(b77)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c77_3 = _mm_add_sd(c77_3, _mm_mul_sd(a77_3, b77));
#endif
_mm_store_sd(&C[(i*88)+77], c77_3);
#else
C[(i*88)+13] += values[298] * B[(i*88)+77];
C[(i*88)+28] += values[299] * B[(i*88)+77];
C[(i*88)+49] += values[300] * B[(i*88)+77];
C[(i*88)+77] += values[301] * B[(i*88)+77];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b78 = _mm256_broadcast_sd(&B[(i*88)+78]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b78 = _mm_loaddup_pd(&B[(i*88)+78]);
#endif
__m128d c78_0 = _mm_load_sd(&C[(i*88)+4]);
__m128d a78_0 = _mm_load_sd(&values[302]);
#if defined(__SSE3__) && defined(__AVX__)
c78_0 = _mm_add_sd(c78_0, _mm_mul_sd(a78_0, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c78_0 = _mm_add_sd(c78_0, _mm_mul_sd(a78_0, b78));
#endif
_mm_store_sd(&C[(i*88)+4], c78_0);
__m128d c78_1 = _mm_load_sd(&C[(i*88)+14]);
__m128d a78_1 = _mm_load_sd(&values[303]);
#if defined(__SSE3__) && defined(__AVX__)
c78_1 = _mm_add_sd(c78_1, _mm_mul_sd(a78_1, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c78_1 = _mm_add_sd(c78_1, _mm_mul_sd(a78_1, b78));
#endif
_mm_store_sd(&C[(i*88)+14], c78_1);
__m128d c78_2 = _mm_load_sd(&C[(i*88)+29]);
__m128d a78_2 = _mm_load_sd(&values[304]);
#if defined(__SSE3__) && defined(__AVX__)
c78_2 = _mm_add_sd(c78_2, _mm_mul_sd(a78_2, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c78_2 = _mm_add_sd(c78_2, _mm_mul_sd(a78_2, b78));
#endif
_mm_store_sd(&C[(i*88)+29], c78_2);
__m128d c78_3 = _mm_load_sd(&C[(i*88)+50]);
__m128d a78_3 = _mm_load_sd(&values[305]);
#if defined(__SSE3__) && defined(__AVX__)
c78_3 = _mm_add_sd(c78_3, _mm_mul_sd(a78_3, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c78_3 = _mm_add_sd(c78_3, _mm_mul_sd(a78_3, b78));
#endif
_mm_store_sd(&C[(i*88)+50], c78_3);
__m128d c78_4 = _mm_load_sd(&C[(i*88)+78]);
__m128d a78_4 = _mm_load_sd(&values[306]);
#if defined(__SSE3__) && defined(__AVX__)
c78_4 = _mm_add_sd(c78_4, _mm_mul_sd(a78_4, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c78_4 = _mm_add_sd(c78_4, _mm_mul_sd(a78_4, b78));
#endif
_mm_store_sd(&C[(i*88)+78], c78_4);
#else
C[(i*88)+4] += values[302] * B[(i*88)+78];
C[(i*88)+14] += values[303] * B[(i*88)+78];
C[(i*88)+29] += values[304] * B[(i*88)+78];
C[(i*88)+50] += values[305] * B[(i*88)+78];
C[(i*88)+78] += values[306] * B[(i*88)+78];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b79 = _mm256_broadcast_sd(&B[(i*88)+79]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b79 = _mm_loaddup_pd(&B[(i*88)+79]);
#endif
__m128d c79_0 = _mm_load_sd(&C[(i*88)+5]);
__m128d a79_0 = _mm_load_sd(&values[307]);
#if defined(__SSE3__) && defined(__AVX__)
c79_0 = _mm_add_sd(c79_0, _mm_mul_sd(a79_0, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c79_0 = _mm_add_sd(c79_0, _mm_mul_sd(a79_0, b79));
#endif
_mm_store_sd(&C[(i*88)+5], c79_0);
__m128d c79_1 = _mm_load_sd(&C[(i*88)+15]);
__m128d a79_1 = _mm_load_sd(&values[308]);
#if defined(__SSE3__) && defined(__AVX__)
c79_1 = _mm_add_sd(c79_1, _mm_mul_sd(a79_1, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c79_1 = _mm_add_sd(c79_1, _mm_mul_sd(a79_1, b79));
#endif
_mm_store_sd(&C[(i*88)+15], c79_1);
__m128d c79_2 = _mm_load_sd(&C[(i*88)+30]);
__m128d a79_2 = _mm_load_sd(&values[309]);
#if defined(__SSE3__) && defined(__AVX__)
c79_2 = _mm_add_sd(c79_2, _mm_mul_sd(a79_2, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c79_2 = _mm_add_sd(c79_2, _mm_mul_sd(a79_2, b79));
#endif
_mm_store_sd(&C[(i*88)+30], c79_2);
__m128d c79_3 = _mm_load_sd(&C[(i*88)+51]);
__m128d a79_3 = _mm_load_sd(&values[310]);
#if defined(__SSE3__) && defined(__AVX__)
c79_3 = _mm_add_sd(c79_3, _mm_mul_sd(a79_3, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c79_3 = _mm_add_sd(c79_3, _mm_mul_sd(a79_3, b79));
#endif
_mm_store_sd(&C[(i*88)+51], c79_3);
__m128d c79_4 = _mm_load_sd(&C[(i*88)+79]);
__m128d a79_4 = _mm_load_sd(&values[311]);
#if defined(__SSE3__) && defined(__AVX__)
c79_4 = _mm_add_sd(c79_4, _mm_mul_sd(a79_4, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c79_4 = _mm_add_sd(c79_4, _mm_mul_sd(a79_4, b79));
#endif
_mm_store_sd(&C[(i*88)+79], c79_4);
#else
C[(i*88)+5] += values[307] * B[(i*88)+79];
C[(i*88)+15] += values[308] * B[(i*88)+79];
C[(i*88)+30] += values[309] * B[(i*88)+79];
C[(i*88)+51] += values[310] * B[(i*88)+79];
C[(i*88)+79] += values[311] * B[(i*88)+79];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b80 = _mm256_broadcast_sd(&B[(i*88)+80]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b80 = _mm_loaddup_pd(&B[(i*88)+80]);
#endif
__m128d c80_0 = _mm_load_sd(&C[(i*88)+6]);
__m128d a80_0 = _mm_load_sd(&values[312]);
#if defined(__SSE3__) && defined(__AVX__)
c80_0 = _mm_add_sd(c80_0, _mm_mul_sd(a80_0, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c80_0 = _mm_add_sd(c80_0, _mm_mul_sd(a80_0, b80));
#endif
_mm_store_sd(&C[(i*88)+6], c80_0);
__m128d c80_1 = _mm_load_sd(&C[(i*88)+16]);
__m128d a80_1 = _mm_load_sd(&values[313]);
#if defined(__SSE3__) && defined(__AVX__)
c80_1 = _mm_add_sd(c80_1, _mm_mul_sd(a80_1, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c80_1 = _mm_add_sd(c80_1, _mm_mul_sd(a80_1, b80));
#endif
_mm_store_sd(&C[(i*88)+16], c80_1);
__m128d c80_2 = _mm_load_sd(&C[(i*88)+31]);
__m128d a80_2 = _mm_load_sd(&values[314]);
#if defined(__SSE3__) && defined(__AVX__)
c80_2 = _mm_add_sd(c80_2, _mm_mul_sd(a80_2, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c80_2 = _mm_add_sd(c80_2, _mm_mul_sd(a80_2, b80));
#endif
_mm_store_sd(&C[(i*88)+31], c80_2);
__m128d c80_3 = _mm_load_sd(&C[(i*88)+52]);
__m128d a80_3 = _mm_load_sd(&values[315]);
#if defined(__SSE3__) && defined(__AVX__)
c80_3 = _mm_add_sd(c80_3, _mm_mul_sd(a80_3, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c80_3 = _mm_add_sd(c80_3, _mm_mul_sd(a80_3, b80));
#endif
_mm_store_sd(&C[(i*88)+52], c80_3);
__m128d c80_4 = _mm_load_sd(&C[(i*88)+80]);
__m128d a80_4 = _mm_load_sd(&values[316]);
#if defined(__SSE3__) && defined(__AVX__)
c80_4 = _mm_add_sd(c80_4, _mm_mul_sd(a80_4, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c80_4 = _mm_add_sd(c80_4, _mm_mul_sd(a80_4, b80));
#endif
_mm_store_sd(&C[(i*88)+80], c80_4);
#else
C[(i*88)+6] += values[312] * B[(i*88)+80];
C[(i*88)+16] += values[313] * B[(i*88)+80];
C[(i*88)+31] += values[314] * B[(i*88)+80];
C[(i*88)+52] += values[315] * B[(i*88)+80];
C[(i*88)+80] += values[316] * B[(i*88)+80];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b81 = _mm256_broadcast_sd(&B[(i*88)+81]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b81 = _mm_loaddup_pd(&B[(i*88)+81]);
#endif
__m128d c81_0 = _mm_load_sd(&C[(i*88)+1]);
__m128d a81_0 = _mm_load_sd(&values[317]);
#if defined(__SSE3__) && defined(__AVX__)
c81_0 = _mm_add_sd(c81_0, _mm_mul_sd(a81_0, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c81_0 = _mm_add_sd(c81_0, _mm_mul_sd(a81_0, b81));
#endif
_mm_store_sd(&C[(i*88)+1], c81_0);
__m128d c81_1 = _mm_load_sd(&C[(i*88)+7]);
__m128d a81_1 = _mm_load_sd(&values[318]);
#if defined(__SSE3__) && defined(__AVX__)
c81_1 = _mm_add_sd(c81_1, _mm_mul_sd(a81_1, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c81_1 = _mm_add_sd(c81_1, _mm_mul_sd(a81_1, b81));
#endif
_mm_store_sd(&C[(i*88)+7], c81_1);
__m128d c81_2 = _mm_load_sd(&C[(i*88)+17]);
__m128d a81_2 = _mm_load_sd(&values[319]);
#if defined(__SSE3__) && defined(__AVX__)
c81_2 = _mm_add_sd(c81_2, _mm_mul_sd(a81_2, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c81_2 = _mm_add_sd(c81_2, _mm_mul_sd(a81_2, b81));
#endif
_mm_store_sd(&C[(i*88)+17], c81_2);
__m128d c81_3 = _mm_load_sd(&C[(i*88)+32]);
__m128d a81_3 = _mm_load_sd(&values[320]);
#if defined(__SSE3__) && defined(__AVX__)
c81_3 = _mm_add_sd(c81_3, _mm_mul_sd(a81_3, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c81_3 = _mm_add_sd(c81_3, _mm_mul_sd(a81_3, b81));
#endif
_mm_store_sd(&C[(i*88)+32], c81_3);
__m128d c81_4 = _mm_load_sd(&C[(i*88)+53]);
__m128d a81_4 = _mm_load_sd(&values[321]);
#if defined(__SSE3__) && defined(__AVX__)
c81_4 = _mm_add_sd(c81_4, _mm_mul_sd(a81_4, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c81_4 = _mm_add_sd(c81_4, _mm_mul_sd(a81_4, b81));
#endif
_mm_store_sd(&C[(i*88)+53], c81_4);
__m128d c81_5 = _mm_load_sd(&C[(i*88)+81]);
__m128d a81_5 = _mm_load_sd(&values[322]);
#if defined(__SSE3__) && defined(__AVX__)
c81_5 = _mm_add_sd(c81_5, _mm_mul_sd(a81_5, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c81_5 = _mm_add_sd(c81_5, _mm_mul_sd(a81_5, b81));
#endif
_mm_store_sd(&C[(i*88)+81], c81_5);
#else
C[(i*88)+1] += values[317] * B[(i*88)+81];
C[(i*88)+7] += values[318] * B[(i*88)+81];
C[(i*88)+17] += values[319] * B[(i*88)+81];
C[(i*88)+32] += values[320] * B[(i*88)+81];
C[(i*88)+53] += values[321] * B[(i*88)+81];
C[(i*88)+81] += values[322] * B[(i*88)+81];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b82 = _mm256_broadcast_sd(&B[(i*88)+82]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b82 = _mm_loaddup_pd(&B[(i*88)+82]);
#endif
__m128d c82_0 = _mm_load_sd(&C[(i*88)+2]);
__m128d a82_0 = _mm_load_sd(&values[323]);
#if defined(__SSE3__) && defined(__AVX__)
c82_0 = _mm_add_sd(c82_0, _mm_mul_sd(a82_0, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c82_0 = _mm_add_sd(c82_0, _mm_mul_sd(a82_0, b82));
#endif
_mm_store_sd(&C[(i*88)+2], c82_0);
__m128d c82_1 = _mm_load_sd(&C[(i*88)+8]);
__m128d a82_1 = _mm_load_sd(&values[324]);
#if defined(__SSE3__) && defined(__AVX__)
c82_1 = _mm_add_sd(c82_1, _mm_mul_sd(a82_1, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c82_1 = _mm_add_sd(c82_1, _mm_mul_sd(a82_1, b82));
#endif
_mm_store_sd(&C[(i*88)+8], c82_1);
__m128d c82_2 = _mm_load_sd(&C[(i*88)+18]);
__m128d a82_2 = _mm_load_sd(&values[325]);
#if defined(__SSE3__) && defined(__AVX__)
c82_2 = _mm_add_sd(c82_2, _mm_mul_sd(a82_2, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c82_2 = _mm_add_sd(c82_2, _mm_mul_sd(a82_2, b82));
#endif
_mm_store_sd(&C[(i*88)+18], c82_2);
__m128d c82_3 = _mm_load_sd(&C[(i*88)+33]);
__m128d a82_3 = _mm_load_sd(&values[326]);
#if defined(__SSE3__) && defined(__AVX__)
c82_3 = _mm_add_sd(c82_3, _mm_mul_sd(a82_3, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c82_3 = _mm_add_sd(c82_3, _mm_mul_sd(a82_3, b82));
#endif
_mm_store_sd(&C[(i*88)+33], c82_3);
__m128d c82_4 = _mm_load_sd(&C[(i*88)+54]);
__m128d a82_4 = _mm_load_sd(&values[327]);
#if defined(__SSE3__) && defined(__AVX__)
c82_4 = _mm_add_sd(c82_4, _mm_mul_sd(a82_4, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c82_4 = _mm_add_sd(c82_4, _mm_mul_sd(a82_4, b82));
#endif
_mm_store_sd(&C[(i*88)+54], c82_4);
__m128d c82_5 = _mm_load_sd(&C[(i*88)+82]);
__m128d a82_5 = _mm_load_sd(&values[328]);
#if defined(__SSE3__) && defined(__AVX__)
c82_5 = _mm_add_sd(c82_5, _mm_mul_sd(a82_5, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c82_5 = _mm_add_sd(c82_5, _mm_mul_sd(a82_5, b82));
#endif
_mm_store_sd(&C[(i*88)+82], c82_5);
#else
C[(i*88)+2] += values[323] * B[(i*88)+82];
C[(i*88)+8] += values[324] * B[(i*88)+82];
C[(i*88)+18] += values[325] * B[(i*88)+82];
C[(i*88)+33] += values[326] * B[(i*88)+82];
C[(i*88)+54] += values[327] * B[(i*88)+82];
C[(i*88)+82] += values[328] * B[(i*88)+82];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b83 = _mm256_broadcast_sd(&B[(i*88)+83]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b83 = _mm_loaddup_pd(&B[(i*88)+83]);
#endif
__m128d c83_0 = _mm_load_sd(&C[(i*88)+0]);
__m128d a83_0 = _mm_load_sd(&values[329]);
#if defined(__SSE3__) && defined(__AVX__)
c83_0 = _mm_add_sd(c83_0, _mm_mul_sd(a83_0, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c83_0 = _mm_add_sd(c83_0, _mm_mul_sd(a83_0, b83));
#endif
_mm_store_sd(&C[(i*88)+0], c83_0);
__m128d c83_1 = _mm_load_sd(&C[(i*88)+3]);
__m128d a83_1 = _mm_load_sd(&values[330]);
#if defined(__SSE3__) && defined(__AVX__)
c83_1 = _mm_add_sd(c83_1, _mm_mul_sd(a83_1, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c83_1 = _mm_add_sd(c83_1, _mm_mul_sd(a83_1, b83));
#endif
_mm_store_sd(&C[(i*88)+3], c83_1);
__m128d c83_2 = _mm_load_sd(&C[(i*88)+9]);
__m128d a83_2 = _mm_load_sd(&values[331]);
#if defined(__SSE3__) && defined(__AVX__)
c83_2 = _mm_add_sd(c83_2, _mm_mul_sd(a83_2, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c83_2 = _mm_add_sd(c83_2, _mm_mul_sd(a83_2, b83));
#endif
_mm_store_sd(&C[(i*88)+9], c83_2);
__m128d c83_3 = _mm_load_sd(&C[(i*88)+19]);
__m128d a83_3 = _mm_load_sd(&values[332]);
#if defined(__SSE3__) && defined(__AVX__)
c83_3 = _mm_add_sd(c83_3, _mm_mul_sd(a83_3, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c83_3 = _mm_add_sd(c83_3, _mm_mul_sd(a83_3, b83));
#endif
_mm_store_sd(&C[(i*88)+19], c83_3);
__m128d c83_4 = _mm_load_sd(&C[(i*88)+34]);
__m128d a83_4 = _mm_load_sd(&values[333]);
#if defined(__SSE3__) && defined(__AVX__)
c83_4 = _mm_add_sd(c83_4, _mm_mul_sd(a83_4, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c83_4 = _mm_add_sd(c83_4, _mm_mul_sd(a83_4, b83));
#endif
_mm_store_sd(&C[(i*88)+34], c83_4);
__m128d c83_5 = _mm_load_sd(&C[(i*88)+55]);
__m128d a83_5 = _mm_load_sd(&values[334]);
#if defined(__SSE3__) && defined(__AVX__)
c83_5 = _mm_add_sd(c83_5, _mm_mul_sd(a83_5, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c83_5 = _mm_add_sd(c83_5, _mm_mul_sd(a83_5, b83));
#endif
_mm_store_sd(&C[(i*88)+55], c83_5);
__m128d c83_6 = _mm_load_sd(&C[(i*88)+83]);
__m128d a83_6 = _mm_load_sd(&values[335]);
#if defined(__SSE3__) && defined(__AVX__)
c83_6 = _mm_add_sd(c83_6, _mm_mul_sd(a83_6, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c83_6 = _mm_add_sd(c83_6, _mm_mul_sd(a83_6, b83));
#endif
_mm_store_sd(&C[(i*88)+83], c83_6);
#else
C[(i*88)+0] += values[329] * B[(i*88)+83];
C[(i*88)+3] += values[330] * B[(i*88)+83];
C[(i*88)+9] += values[331] * B[(i*88)+83];
C[(i*88)+19] += values[332] * B[(i*88)+83];
C[(i*88)+34] += values[333] * B[(i*88)+83];
C[(i*88)+55] += values[334] * B[(i*88)+83];
C[(i*88)+83] += values[335] * B[(i*88)+83];
#endif

}

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 6048;
#endif

}

void dsparse_fP111DivM_m84_n9_k84_ldAna7_ldB88_ldC88_beta0_pfsigonly(const double* values, const double* B, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma nounroll_and_jam
for (unsigned int i = 0; i < 9; i++)
{
  #pragma simd
  for (unsigned int m = 0; m < 84; m++) {
    C[(i*88)+m] = 0.0;
  }
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b0 = _mm256_broadcast_sd(&B[(i*88)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b0 = _mm_loaddup_pd(&B[(i*88)+0]);
#endif
__m128d c0_0 = _mm_load_sd(&C[(i*88)+0]);
__m128d a0_0 = _mm_load_sd(&values[0]);
#if defined(__SSE3__) && defined(__AVX__)
c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, b0));
#endif
_mm_store_sd(&C[(i*88)+0], c0_0);
__m128d c0_1 = _mm_load_sd(&C[(i*88)+3]);
__m128d a0_1 = _mm_load_sd(&values[1]);
#if defined(__SSE3__) && defined(__AVX__)
c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, b0));
#endif
_mm_store_sd(&C[(i*88)+3], c0_1);
__m128d c0_2 = _mm_load_sd(&C[(i*88)+9]);
__m128d a0_2 = _mm_load_sd(&values[2]);
#if defined(__SSE3__) && defined(__AVX__)
c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, b0));
#endif
_mm_store_sd(&C[(i*88)+9], c0_2);
__m128d c0_3 = _mm_load_sd(&C[(i*88)+19]);
__m128d a0_3 = _mm_load_sd(&values[3]);
#if defined(__SSE3__) && defined(__AVX__)
c0_3 = _mm_add_sd(c0_3, _mm_mul_sd(a0_3, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_3 = _mm_add_sd(c0_3, _mm_mul_sd(a0_3, b0));
#endif
_mm_store_sd(&C[(i*88)+19], c0_3);
__m128d c0_4 = _mm_load_sd(&C[(i*88)+34]);
__m128d a0_4 = _mm_load_sd(&values[4]);
#if defined(__SSE3__) && defined(__AVX__)
c0_4 = _mm_add_sd(c0_4, _mm_mul_sd(a0_4, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_4 = _mm_add_sd(c0_4, _mm_mul_sd(a0_4, b0));
#endif
_mm_store_sd(&C[(i*88)+34], c0_4);
__m128d c0_5 = _mm_load_sd(&C[(i*88)+55]);
__m128d a0_5 = _mm_load_sd(&values[5]);
#if defined(__SSE3__) && defined(__AVX__)
c0_5 = _mm_add_sd(c0_5, _mm_mul_sd(a0_5, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_5 = _mm_add_sd(c0_5, _mm_mul_sd(a0_5, b0));
#endif
_mm_store_sd(&C[(i*88)+55], c0_5);
__m128d c0_6 = _mm_load_sd(&C[(i*88)+83]);
__m128d a0_6 = _mm_load_sd(&values[6]);
#if defined(__SSE3__) && defined(__AVX__)
c0_6 = _mm_add_sd(c0_6, _mm_mul_sd(a0_6, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_6 = _mm_add_sd(c0_6, _mm_mul_sd(a0_6, b0));
#endif
_mm_store_sd(&C[(i*88)+83], c0_6);
#else
C[(i*88)+0] += values[0] * B[(i*88)+0];
C[(i*88)+3] += values[1] * B[(i*88)+0];
C[(i*88)+9] += values[2] * B[(i*88)+0];
C[(i*88)+19] += values[3] * B[(i*88)+0];
C[(i*88)+34] += values[4] * B[(i*88)+0];
C[(i*88)+55] += values[5] * B[(i*88)+0];
C[(i*88)+83] += values[6] * B[(i*88)+0];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b1 = _mm256_broadcast_sd(&B[(i*88)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b1 = _mm_loaddup_pd(&B[(i*88)+1]);
#endif
__m128d c1_0 = _mm_loadu_pd(&C[(i*88)+1]);
__m128d a1_0 = _mm_loadu_pd(&values[7]);
#if defined(__SSE3__) && defined(__AVX__)
c1_0 = _mm_add_pd(c1_0, _mm_mul_pd(a1_0, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_0 = _mm_add_pd(c1_0, _mm_mul_pd(a1_0, b1));
#endif
_mm_storeu_pd(&C[(i*88)+1], c1_0);
__m128d c1_2 = _mm_loadu_pd(&C[(i*88)+7]);
__m128d a1_2 = _mm_loadu_pd(&values[9]);
#if defined(__SSE3__) && defined(__AVX__)
c1_2 = _mm_add_pd(c1_2, _mm_mul_pd(a1_2, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_2 = _mm_add_pd(c1_2, _mm_mul_pd(a1_2, b1));
#endif
_mm_storeu_pd(&C[(i*88)+7], c1_2);
__m128d c1_4 = _mm_loadu_pd(&C[(i*88)+17]);
__m128d a1_4 = _mm_loadu_pd(&values[11]);
#if defined(__SSE3__) && defined(__AVX__)
c1_4 = _mm_add_pd(c1_4, _mm_mul_pd(a1_4, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_4 = _mm_add_pd(c1_4, _mm_mul_pd(a1_4, b1));
#endif
_mm_storeu_pd(&C[(i*88)+17], c1_4);
__m128d c1_6 = _mm_loadu_pd(&C[(i*88)+32]);
__m128d a1_6 = _mm_loadu_pd(&values[13]);
#if defined(__SSE3__) && defined(__AVX__)
c1_6 = _mm_add_pd(c1_6, _mm_mul_pd(a1_6, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_6 = _mm_add_pd(c1_6, _mm_mul_pd(a1_6, b1));
#endif
_mm_storeu_pd(&C[(i*88)+32], c1_6);
__m128d c1_8 = _mm_loadu_pd(&C[(i*88)+53]);
__m128d a1_8 = _mm_loadu_pd(&values[15]);
#if defined(__SSE3__) && defined(__AVX__)
c1_8 = _mm_add_pd(c1_8, _mm_mul_pd(a1_8, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_8 = _mm_add_pd(c1_8, _mm_mul_pd(a1_8, b1));
#endif
_mm_storeu_pd(&C[(i*88)+53], c1_8);
__m128d c1_10 = _mm_loadu_pd(&C[(i*88)+81]);
__m128d a1_10 = _mm_loadu_pd(&values[17]);
#if defined(__SSE3__) && defined(__AVX__)
c1_10 = _mm_add_pd(c1_10, _mm_mul_pd(a1_10, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_10 = _mm_add_pd(c1_10, _mm_mul_pd(a1_10, b1));
#endif
_mm_storeu_pd(&C[(i*88)+81], c1_10);
#else
C[(i*88)+1] += values[7] * B[(i*88)+1];
C[(i*88)+2] += values[8] * B[(i*88)+1];
C[(i*88)+7] += values[9] * B[(i*88)+1];
C[(i*88)+8] += values[10] * B[(i*88)+1];
C[(i*88)+17] += values[11] * B[(i*88)+1];
C[(i*88)+18] += values[12] * B[(i*88)+1];
C[(i*88)+32] += values[13] * B[(i*88)+1];
C[(i*88)+33] += values[14] * B[(i*88)+1];
C[(i*88)+53] += values[15] * B[(i*88)+1];
C[(i*88)+54] += values[16] * B[(i*88)+1];
C[(i*88)+81] += values[17] * B[(i*88)+1];
C[(i*88)+82] += values[18] * B[(i*88)+1];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b2 = _mm256_broadcast_sd(&B[(i*88)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b2 = _mm_loaddup_pd(&B[(i*88)+2]);
#endif
__m128d c2_0 = _mm_loadu_pd(&C[(i*88)+1]);
__m128d a2_0 = _mm_loadu_pd(&values[19]);
#if defined(__SSE3__) && defined(__AVX__)
c2_0 = _mm_add_pd(c2_0, _mm_mul_pd(a2_0, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_0 = _mm_add_pd(c2_0, _mm_mul_pd(a2_0, b2));
#endif
_mm_storeu_pd(&C[(i*88)+1], c2_0);
__m128d c2_2 = _mm_loadu_pd(&C[(i*88)+7]);
__m128d a2_2 = _mm_loadu_pd(&values[21]);
#if defined(__SSE3__) && defined(__AVX__)
c2_2 = _mm_add_pd(c2_2, _mm_mul_pd(a2_2, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_2 = _mm_add_pd(c2_2, _mm_mul_pd(a2_2, b2));
#endif
_mm_storeu_pd(&C[(i*88)+7], c2_2);
__m128d c2_4 = _mm_loadu_pd(&C[(i*88)+17]);
__m128d a2_4 = _mm_loadu_pd(&values[23]);
#if defined(__SSE3__) && defined(__AVX__)
c2_4 = _mm_add_pd(c2_4, _mm_mul_pd(a2_4, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_4 = _mm_add_pd(c2_4, _mm_mul_pd(a2_4, b2));
#endif
_mm_storeu_pd(&C[(i*88)+17], c2_4);
__m128d c2_6 = _mm_loadu_pd(&C[(i*88)+32]);
__m128d a2_6 = _mm_loadu_pd(&values[25]);
#if defined(__SSE3__) && defined(__AVX__)
c2_6 = _mm_add_pd(c2_6, _mm_mul_pd(a2_6, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_6 = _mm_add_pd(c2_6, _mm_mul_pd(a2_6, b2));
#endif
_mm_storeu_pd(&C[(i*88)+32], c2_6);
__m128d c2_8 = _mm_loadu_pd(&C[(i*88)+53]);
__m128d a2_8 = _mm_loadu_pd(&values[27]);
#if defined(__SSE3__) && defined(__AVX__)
c2_8 = _mm_add_pd(c2_8, _mm_mul_pd(a2_8, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_8 = _mm_add_pd(c2_8, _mm_mul_pd(a2_8, b2));
#endif
_mm_storeu_pd(&C[(i*88)+53], c2_8);
__m128d c2_10 = _mm_loadu_pd(&C[(i*88)+81]);
__m128d a2_10 = _mm_loadu_pd(&values[29]);
#if defined(__SSE3__) && defined(__AVX__)
c2_10 = _mm_add_pd(c2_10, _mm_mul_pd(a2_10, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_10 = _mm_add_pd(c2_10, _mm_mul_pd(a2_10, b2));
#endif
_mm_storeu_pd(&C[(i*88)+81], c2_10);
#else
C[(i*88)+1] += values[19] * B[(i*88)+2];
C[(i*88)+2] += values[20] * B[(i*88)+2];
C[(i*88)+7] += values[21] * B[(i*88)+2];
C[(i*88)+8] += values[22] * B[(i*88)+2];
C[(i*88)+17] += values[23] * B[(i*88)+2];
C[(i*88)+18] += values[24] * B[(i*88)+2];
C[(i*88)+32] += values[25] * B[(i*88)+2];
C[(i*88)+33] += values[26] * B[(i*88)+2];
C[(i*88)+53] += values[27] * B[(i*88)+2];
C[(i*88)+54] += values[28] * B[(i*88)+2];
C[(i*88)+81] += values[29] * B[(i*88)+2];
C[(i*88)+82] += values[30] * B[(i*88)+2];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b3 = _mm256_broadcast_sd(&B[(i*88)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b3 = _mm_loaddup_pd(&B[(i*88)+3]);
#endif
__m128d c3_0 = _mm_load_sd(&C[(i*88)+0]);
__m128d a3_0 = _mm_load_sd(&values[31]);
#if defined(__SSE3__) && defined(__AVX__)
c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, b3));
#endif
_mm_store_sd(&C[(i*88)+0], c3_0);
__m128d c3_1 = _mm_load_sd(&C[(i*88)+3]);
__m128d a3_1 = _mm_load_sd(&values[32]);
#if defined(__SSE3__) && defined(__AVX__)
c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, b3));
#endif
_mm_store_sd(&C[(i*88)+3], c3_1);
__m128d c3_2 = _mm_load_sd(&C[(i*88)+9]);
__m128d a3_2 = _mm_load_sd(&values[33]);
#if defined(__SSE3__) && defined(__AVX__)
c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, b3));
#endif
_mm_store_sd(&C[(i*88)+9], c3_2);
__m128d c3_3 = _mm_load_sd(&C[(i*88)+19]);
__m128d a3_3 = _mm_load_sd(&values[34]);
#if defined(__SSE3__) && defined(__AVX__)
c3_3 = _mm_add_sd(c3_3, _mm_mul_sd(a3_3, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_3 = _mm_add_sd(c3_3, _mm_mul_sd(a3_3, b3));
#endif
_mm_store_sd(&C[(i*88)+19], c3_3);
__m128d c3_4 = _mm_load_sd(&C[(i*88)+34]);
__m128d a3_4 = _mm_load_sd(&values[35]);
#if defined(__SSE3__) && defined(__AVX__)
c3_4 = _mm_add_sd(c3_4, _mm_mul_sd(a3_4, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_4 = _mm_add_sd(c3_4, _mm_mul_sd(a3_4, b3));
#endif
_mm_store_sd(&C[(i*88)+34], c3_4);
__m128d c3_5 = _mm_load_sd(&C[(i*88)+55]);
__m128d a3_5 = _mm_load_sd(&values[36]);
#if defined(__SSE3__) && defined(__AVX__)
c3_5 = _mm_add_sd(c3_5, _mm_mul_sd(a3_5, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_5 = _mm_add_sd(c3_5, _mm_mul_sd(a3_5, b3));
#endif
_mm_store_sd(&C[(i*88)+55], c3_5);
__m128d c3_6 = _mm_load_sd(&C[(i*88)+83]);
__m128d a3_6 = _mm_load_sd(&values[37]);
#if defined(__SSE3__) && defined(__AVX__)
c3_6 = _mm_add_sd(c3_6, _mm_mul_sd(a3_6, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_6 = _mm_add_sd(c3_6, _mm_mul_sd(a3_6, b3));
#endif
_mm_store_sd(&C[(i*88)+83], c3_6);
#else
C[(i*88)+0] += values[31] * B[(i*88)+3];
C[(i*88)+3] += values[32] * B[(i*88)+3];
C[(i*88)+9] += values[33] * B[(i*88)+3];
C[(i*88)+19] += values[34] * B[(i*88)+3];
C[(i*88)+34] += values[35] * B[(i*88)+3];
C[(i*88)+55] += values[36] * B[(i*88)+3];
C[(i*88)+83] += values[37] * B[(i*88)+3];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b4 = _mm256_broadcast_sd(&B[(i*88)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b4 = _mm_loaddup_pd(&B[(i*88)+4]);
#endif
__m128d c4_0 = _mm_loadu_pd(&C[(i*88)+4]);
__m128d a4_0 = _mm_loadu_pd(&values[38]);
#if defined(__SSE3__) && defined(__AVX__)
c4_0 = _mm_add_pd(c4_0, _mm_mul_pd(a4_0, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_0 = _mm_add_pd(c4_0, _mm_mul_pd(a4_0, b4));
#endif
_mm_storeu_pd(&C[(i*88)+4], c4_0);
__m128d c4_2 = _mm_load_sd(&C[(i*88)+6]);
__m128d a4_2 = _mm_load_sd(&values[40]);
#if defined(__SSE3__) && defined(__AVX__)
c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, b4));
#endif
_mm_store_sd(&C[(i*88)+6], c4_2);
__m128d c4_3 = _mm_loadu_pd(&C[(i*88)+14]);
__m128d a4_3 = _mm_loadu_pd(&values[41]);
#if defined(__SSE3__) && defined(__AVX__)
c4_3 = _mm_add_pd(c4_3, _mm_mul_pd(a4_3, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_3 = _mm_add_pd(c4_3, _mm_mul_pd(a4_3, b4));
#endif
_mm_storeu_pd(&C[(i*88)+14], c4_3);
__m128d c4_5 = _mm_load_sd(&C[(i*88)+16]);
__m128d a4_5 = _mm_load_sd(&values[43]);
#if defined(__SSE3__) && defined(__AVX__)
c4_5 = _mm_add_sd(c4_5, _mm_mul_sd(a4_5, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_5 = _mm_add_sd(c4_5, _mm_mul_sd(a4_5, b4));
#endif
_mm_store_sd(&C[(i*88)+16], c4_5);
__m128d c4_6 = _mm_loadu_pd(&C[(i*88)+29]);
__m128d a4_6 = _mm_loadu_pd(&values[44]);
#if defined(__SSE3__) && defined(__AVX__)
c4_6 = _mm_add_pd(c4_6, _mm_mul_pd(a4_6, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_6 = _mm_add_pd(c4_6, _mm_mul_pd(a4_6, b4));
#endif
_mm_storeu_pd(&C[(i*88)+29], c4_6);
__m128d c4_8 = _mm_load_sd(&C[(i*88)+31]);
__m128d a4_8 = _mm_load_sd(&values[46]);
#if defined(__SSE3__) && defined(__AVX__)
c4_8 = _mm_add_sd(c4_8, _mm_mul_sd(a4_8, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_8 = _mm_add_sd(c4_8, _mm_mul_sd(a4_8, b4));
#endif
_mm_store_sd(&C[(i*88)+31], c4_8);
__m128d c4_9 = _mm_loadu_pd(&C[(i*88)+50]);
__m128d a4_9 = _mm_loadu_pd(&values[47]);
#if defined(__SSE3__) && defined(__AVX__)
c4_9 = _mm_add_pd(c4_9, _mm_mul_pd(a4_9, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_9 = _mm_add_pd(c4_9, _mm_mul_pd(a4_9, b4));
#endif
_mm_storeu_pd(&C[(i*88)+50], c4_9);
__m128d c4_11 = _mm_load_sd(&C[(i*88)+52]);
__m128d a4_11 = _mm_load_sd(&values[49]);
#if defined(__SSE3__) && defined(__AVX__)
c4_11 = _mm_add_sd(c4_11, _mm_mul_sd(a4_11, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_11 = _mm_add_sd(c4_11, _mm_mul_sd(a4_11, b4));
#endif
_mm_store_sd(&C[(i*88)+52], c4_11);
__m128d c4_12 = _mm_loadu_pd(&C[(i*88)+78]);
__m128d a4_12 = _mm_loadu_pd(&values[50]);
#if defined(__SSE3__) && defined(__AVX__)
c4_12 = _mm_add_pd(c4_12, _mm_mul_pd(a4_12, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_12 = _mm_add_pd(c4_12, _mm_mul_pd(a4_12, b4));
#endif
_mm_storeu_pd(&C[(i*88)+78], c4_12);
__m128d c4_14 = _mm_load_sd(&C[(i*88)+80]);
__m128d a4_14 = _mm_load_sd(&values[52]);
#if defined(__SSE3__) && defined(__AVX__)
c4_14 = _mm_add_sd(c4_14, _mm_mul_sd(a4_14, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_14 = _mm_add_sd(c4_14, _mm_mul_sd(a4_14, b4));
#endif
_mm_store_sd(&C[(i*88)+80], c4_14);
#else
C[(i*88)+4] += values[38] * B[(i*88)+4];
C[(i*88)+5] += values[39] * B[(i*88)+4];
C[(i*88)+6] += values[40] * B[(i*88)+4];
C[(i*88)+14] += values[41] * B[(i*88)+4];
C[(i*88)+15] += values[42] * B[(i*88)+4];
C[(i*88)+16] += values[43] * B[(i*88)+4];
C[(i*88)+29] += values[44] * B[(i*88)+4];
C[(i*88)+30] += values[45] * B[(i*88)+4];
C[(i*88)+31] += values[46] * B[(i*88)+4];
C[(i*88)+50] += values[47] * B[(i*88)+4];
C[(i*88)+51] += values[48] * B[(i*88)+4];
C[(i*88)+52] += values[49] * B[(i*88)+4];
C[(i*88)+78] += values[50] * B[(i*88)+4];
C[(i*88)+79] += values[51] * B[(i*88)+4];
C[(i*88)+80] += values[52] * B[(i*88)+4];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b5 = _mm256_broadcast_sd(&B[(i*88)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b5 = _mm_loaddup_pd(&B[(i*88)+5]);
#endif
__m128d c5_0 = _mm_loadu_pd(&C[(i*88)+4]);
__m128d a5_0 = _mm_loadu_pd(&values[53]);
#if defined(__SSE3__) && defined(__AVX__)
c5_0 = _mm_add_pd(c5_0, _mm_mul_pd(a5_0, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_0 = _mm_add_pd(c5_0, _mm_mul_pd(a5_0, b5));
#endif
_mm_storeu_pd(&C[(i*88)+4], c5_0);
__m128d c5_2 = _mm_load_sd(&C[(i*88)+6]);
__m128d a5_2 = _mm_load_sd(&values[55]);
#if defined(__SSE3__) && defined(__AVX__)
c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, b5));
#endif
_mm_store_sd(&C[(i*88)+6], c5_2);
__m128d c5_3 = _mm_loadu_pd(&C[(i*88)+14]);
__m128d a5_3 = _mm_loadu_pd(&values[56]);
#if defined(__SSE3__) && defined(__AVX__)
c5_3 = _mm_add_pd(c5_3, _mm_mul_pd(a5_3, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_3 = _mm_add_pd(c5_3, _mm_mul_pd(a5_3, b5));
#endif
_mm_storeu_pd(&C[(i*88)+14], c5_3);
__m128d c5_5 = _mm_load_sd(&C[(i*88)+16]);
__m128d a5_5 = _mm_load_sd(&values[58]);
#if defined(__SSE3__) && defined(__AVX__)
c5_5 = _mm_add_sd(c5_5, _mm_mul_sd(a5_5, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_5 = _mm_add_sd(c5_5, _mm_mul_sd(a5_5, b5));
#endif
_mm_store_sd(&C[(i*88)+16], c5_5);
__m128d c5_6 = _mm_loadu_pd(&C[(i*88)+29]);
__m128d a5_6 = _mm_loadu_pd(&values[59]);
#if defined(__SSE3__) && defined(__AVX__)
c5_6 = _mm_add_pd(c5_6, _mm_mul_pd(a5_6, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_6 = _mm_add_pd(c5_6, _mm_mul_pd(a5_6, b5));
#endif
_mm_storeu_pd(&C[(i*88)+29], c5_6);
__m128d c5_8 = _mm_load_sd(&C[(i*88)+31]);
__m128d a5_8 = _mm_load_sd(&values[61]);
#if defined(__SSE3__) && defined(__AVX__)
c5_8 = _mm_add_sd(c5_8, _mm_mul_sd(a5_8, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_8 = _mm_add_sd(c5_8, _mm_mul_sd(a5_8, b5));
#endif
_mm_store_sd(&C[(i*88)+31], c5_8);
__m128d c5_9 = _mm_loadu_pd(&C[(i*88)+50]);
__m128d a5_9 = _mm_loadu_pd(&values[62]);
#if defined(__SSE3__) && defined(__AVX__)
c5_9 = _mm_add_pd(c5_9, _mm_mul_pd(a5_9, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_9 = _mm_add_pd(c5_9, _mm_mul_pd(a5_9, b5));
#endif
_mm_storeu_pd(&C[(i*88)+50], c5_9);
__m128d c5_11 = _mm_load_sd(&C[(i*88)+52]);
__m128d a5_11 = _mm_load_sd(&values[64]);
#if defined(__SSE3__) && defined(__AVX__)
c5_11 = _mm_add_sd(c5_11, _mm_mul_sd(a5_11, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_11 = _mm_add_sd(c5_11, _mm_mul_sd(a5_11, b5));
#endif
_mm_store_sd(&C[(i*88)+52], c5_11);
__m128d c5_12 = _mm_loadu_pd(&C[(i*88)+78]);
__m128d a5_12 = _mm_loadu_pd(&values[65]);
#if defined(__SSE3__) && defined(__AVX__)
c5_12 = _mm_add_pd(c5_12, _mm_mul_pd(a5_12, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_12 = _mm_add_pd(c5_12, _mm_mul_pd(a5_12, b5));
#endif
_mm_storeu_pd(&C[(i*88)+78], c5_12);
__m128d c5_14 = _mm_load_sd(&C[(i*88)+80]);
__m128d a5_14 = _mm_load_sd(&values[67]);
#if defined(__SSE3__) && defined(__AVX__)
c5_14 = _mm_add_sd(c5_14, _mm_mul_sd(a5_14, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_14 = _mm_add_sd(c5_14, _mm_mul_sd(a5_14, b5));
#endif
_mm_store_sd(&C[(i*88)+80], c5_14);
#else
C[(i*88)+4] += values[53] * B[(i*88)+5];
C[(i*88)+5] += values[54] * B[(i*88)+5];
C[(i*88)+6] += values[55] * B[(i*88)+5];
C[(i*88)+14] += values[56] * B[(i*88)+5];
C[(i*88)+15] += values[57] * B[(i*88)+5];
C[(i*88)+16] += values[58] * B[(i*88)+5];
C[(i*88)+29] += values[59] * B[(i*88)+5];
C[(i*88)+30] += values[60] * B[(i*88)+5];
C[(i*88)+31] += values[61] * B[(i*88)+5];
C[(i*88)+50] += values[62] * B[(i*88)+5];
C[(i*88)+51] += values[63] * B[(i*88)+5];
C[(i*88)+52] += values[64] * B[(i*88)+5];
C[(i*88)+78] += values[65] * B[(i*88)+5];
C[(i*88)+79] += values[66] * B[(i*88)+5];
C[(i*88)+80] += values[67] * B[(i*88)+5];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b6 = _mm256_broadcast_sd(&B[(i*88)+6]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b6 = _mm_loaddup_pd(&B[(i*88)+6]);
#endif
__m128d c6_0 = _mm_loadu_pd(&C[(i*88)+4]);
__m128d a6_0 = _mm_loadu_pd(&values[68]);
#if defined(__SSE3__) && defined(__AVX__)
c6_0 = _mm_add_pd(c6_0, _mm_mul_pd(a6_0, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_0 = _mm_add_pd(c6_0, _mm_mul_pd(a6_0, b6));
#endif
_mm_storeu_pd(&C[(i*88)+4], c6_0);
__m128d c6_2 = _mm_load_sd(&C[(i*88)+6]);
__m128d a6_2 = _mm_load_sd(&values[70]);
#if defined(__SSE3__) && defined(__AVX__)
c6_2 = _mm_add_sd(c6_2, _mm_mul_sd(a6_2, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_2 = _mm_add_sd(c6_2, _mm_mul_sd(a6_2, b6));
#endif
_mm_store_sd(&C[(i*88)+6], c6_2);
__m128d c6_3 = _mm_loadu_pd(&C[(i*88)+14]);
__m128d a6_3 = _mm_loadu_pd(&values[71]);
#if defined(__SSE3__) && defined(__AVX__)
c6_3 = _mm_add_pd(c6_3, _mm_mul_pd(a6_3, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_3 = _mm_add_pd(c6_3, _mm_mul_pd(a6_3, b6));
#endif
_mm_storeu_pd(&C[(i*88)+14], c6_3);
__m128d c6_5 = _mm_load_sd(&C[(i*88)+16]);
__m128d a6_5 = _mm_load_sd(&values[73]);
#if defined(__SSE3__) && defined(__AVX__)
c6_5 = _mm_add_sd(c6_5, _mm_mul_sd(a6_5, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_5 = _mm_add_sd(c6_5, _mm_mul_sd(a6_5, b6));
#endif
_mm_store_sd(&C[(i*88)+16], c6_5);
__m128d c6_6 = _mm_loadu_pd(&C[(i*88)+29]);
__m128d a6_6 = _mm_loadu_pd(&values[74]);
#if defined(__SSE3__) && defined(__AVX__)
c6_6 = _mm_add_pd(c6_6, _mm_mul_pd(a6_6, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_6 = _mm_add_pd(c6_6, _mm_mul_pd(a6_6, b6));
#endif
_mm_storeu_pd(&C[(i*88)+29], c6_6);
__m128d c6_8 = _mm_load_sd(&C[(i*88)+31]);
__m128d a6_8 = _mm_load_sd(&values[76]);
#if defined(__SSE3__) && defined(__AVX__)
c6_8 = _mm_add_sd(c6_8, _mm_mul_sd(a6_8, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_8 = _mm_add_sd(c6_8, _mm_mul_sd(a6_8, b6));
#endif
_mm_store_sd(&C[(i*88)+31], c6_8);
__m128d c6_9 = _mm_loadu_pd(&C[(i*88)+50]);
__m128d a6_9 = _mm_loadu_pd(&values[77]);
#if defined(__SSE3__) && defined(__AVX__)
c6_9 = _mm_add_pd(c6_9, _mm_mul_pd(a6_9, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_9 = _mm_add_pd(c6_9, _mm_mul_pd(a6_9, b6));
#endif
_mm_storeu_pd(&C[(i*88)+50], c6_9);
__m128d c6_11 = _mm_load_sd(&C[(i*88)+52]);
__m128d a6_11 = _mm_load_sd(&values[79]);
#if defined(__SSE3__) && defined(__AVX__)
c6_11 = _mm_add_sd(c6_11, _mm_mul_sd(a6_11, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_11 = _mm_add_sd(c6_11, _mm_mul_sd(a6_11, b6));
#endif
_mm_store_sd(&C[(i*88)+52], c6_11);
__m128d c6_12 = _mm_loadu_pd(&C[(i*88)+78]);
__m128d a6_12 = _mm_loadu_pd(&values[80]);
#if defined(__SSE3__) && defined(__AVX__)
c6_12 = _mm_add_pd(c6_12, _mm_mul_pd(a6_12, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_12 = _mm_add_pd(c6_12, _mm_mul_pd(a6_12, b6));
#endif
_mm_storeu_pd(&C[(i*88)+78], c6_12);
__m128d c6_14 = _mm_load_sd(&C[(i*88)+80]);
__m128d a6_14 = _mm_load_sd(&values[82]);
#if defined(__SSE3__) && defined(__AVX__)
c6_14 = _mm_add_sd(c6_14, _mm_mul_sd(a6_14, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_14 = _mm_add_sd(c6_14, _mm_mul_sd(a6_14, b6));
#endif
_mm_store_sd(&C[(i*88)+80], c6_14);
#else
C[(i*88)+4] += values[68] * B[(i*88)+6];
C[(i*88)+5] += values[69] * B[(i*88)+6];
C[(i*88)+6] += values[70] * B[(i*88)+6];
C[(i*88)+14] += values[71] * B[(i*88)+6];
C[(i*88)+15] += values[72] * B[(i*88)+6];
C[(i*88)+16] += values[73] * B[(i*88)+6];
C[(i*88)+29] += values[74] * B[(i*88)+6];
C[(i*88)+30] += values[75] * B[(i*88)+6];
C[(i*88)+31] += values[76] * B[(i*88)+6];
C[(i*88)+50] += values[77] * B[(i*88)+6];
C[(i*88)+51] += values[78] * B[(i*88)+6];
C[(i*88)+52] += values[79] * B[(i*88)+6];
C[(i*88)+78] += values[80] * B[(i*88)+6];
C[(i*88)+79] += values[81] * B[(i*88)+6];
C[(i*88)+80] += values[82] * B[(i*88)+6];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b7 = _mm256_broadcast_sd(&B[(i*88)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b7 = _mm_loaddup_pd(&B[(i*88)+7]);
#endif
__m128d c7_0 = _mm_loadu_pd(&C[(i*88)+1]);
__m128d a7_0 = _mm_loadu_pd(&values[83]);
#if defined(__SSE3__) && defined(__AVX__)
c7_0 = _mm_add_pd(c7_0, _mm_mul_pd(a7_0, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_0 = _mm_add_pd(c7_0, _mm_mul_pd(a7_0, b7));
#endif
_mm_storeu_pd(&C[(i*88)+1], c7_0);
__m128d c7_2 = _mm_loadu_pd(&C[(i*88)+7]);
__m128d a7_2 = _mm_loadu_pd(&values[85]);
#if defined(__SSE3__) && defined(__AVX__)
c7_2 = _mm_add_pd(c7_2, _mm_mul_pd(a7_2, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_2 = _mm_add_pd(c7_2, _mm_mul_pd(a7_2, b7));
#endif
_mm_storeu_pd(&C[(i*88)+7], c7_2);
__m128d c7_4 = _mm_loadu_pd(&C[(i*88)+17]);
__m128d a7_4 = _mm_loadu_pd(&values[87]);
#if defined(__SSE3__) && defined(__AVX__)
c7_4 = _mm_add_pd(c7_4, _mm_mul_pd(a7_4, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_4 = _mm_add_pd(c7_4, _mm_mul_pd(a7_4, b7));
#endif
_mm_storeu_pd(&C[(i*88)+17], c7_4);
__m128d c7_6 = _mm_loadu_pd(&C[(i*88)+32]);
__m128d a7_6 = _mm_loadu_pd(&values[89]);
#if defined(__SSE3__) && defined(__AVX__)
c7_6 = _mm_add_pd(c7_6, _mm_mul_pd(a7_6, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_6 = _mm_add_pd(c7_6, _mm_mul_pd(a7_6, b7));
#endif
_mm_storeu_pd(&C[(i*88)+32], c7_6);
__m128d c7_8 = _mm_loadu_pd(&C[(i*88)+53]);
__m128d a7_8 = _mm_loadu_pd(&values[91]);
#if defined(__SSE3__) && defined(__AVX__)
c7_8 = _mm_add_pd(c7_8, _mm_mul_pd(a7_8, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_8 = _mm_add_pd(c7_8, _mm_mul_pd(a7_8, b7));
#endif
_mm_storeu_pd(&C[(i*88)+53], c7_8);
__m128d c7_10 = _mm_loadu_pd(&C[(i*88)+81]);
__m128d a7_10 = _mm_loadu_pd(&values[93]);
#if defined(__SSE3__) && defined(__AVX__)
c7_10 = _mm_add_pd(c7_10, _mm_mul_pd(a7_10, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_10 = _mm_add_pd(c7_10, _mm_mul_pd(a7_10, b7));
#endif
_mm_storeu_pd(&C[(i*88)+81], c7_10);
#else
C[(i*88)+1] += values[83] * B[(i*88)+7];
C[(i*88)+2] += values[84] * B[(i*88)+7];
C[(i*88)+7] += values[85] * B[(i*88)+7];
C[(i*88)+8] += values[86] * B[(i*88)+7];
C[(i*88)+17] += values[87] * B[(i*88)+7];
C[(i*88)+18] += values[88] * B[(i*88)+7];
C[(i*88)+32] += values[89] * B[(i*88)+7];
C[(i*88)+33] += values[90] * B[(i*88)+7];
C[(i*88)+53] += values[91] * B[(i*88)+7];
C[(i*88)+54] += values[92] * B[(i*88)+7];
C[(i*88)+81] += values[93] * B[(i*88)+7];
C[(i*88)+82] += values[94] * B[(i*88)+7];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b8 = _mm256_broadcast_sd(&B[(i*88)+8]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b8 = _mm_loaddup_pd(&B[(i*88)+8]);
#endif
__m128d c8_0 = _mm_loadu_pd(&C[(i*88)+1]);
__m128d a8_0 = _mm_loadu_pd(&values[95]);
#if defined(__SSE3__) && defined(__AVX__)
c8_0 = _mm_add_pd(c8_0, _mm_mul_pd(a8_0, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_0 = _mm_add_pd(c8_0, _mm_mul_pd(a8_0, b8));
#endif
_mm_storeu_pd(&C[(i*88)+1], c8_0);
__m128d c8_2 = _mm_loadu_pd(&C[(i*88)+7]);
__m128d a8_2 = _mm_loadu_pd(&values[97]);
#if defined(__SSE3__) && defined(__AVX__)
c8_2 = _mm_add_pd(c8_2, _mm_mul_pd(a8_2, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_2 = _mm_add_pd(c8_2, _mm_mul_pd(a8_2, b8));
#endif
_mm_storeu_pd(&C[(i*88)+7], c8_2);
__m128d c8_4 = _mm_loadu_pd(&C[(i*88)+17]);
__m128d a8_4 = _mm_loadu_pd(&values[99]);
#if defined(__SSE3__) && defined(__AVX__)
c8_4 = _mm_add_pd(c8_4, _mm_mul_pd(a8_4, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_4 = _mm_add_pd(c8_4, _mm_mul_pd(a8_4, b8));
#endif
_mm_storeu_pd(&C[(i*88)+17], c8_4);
__m128d c8_6 = _mm_loadu_pd(&C[(i*88)+32]);
__m128d a8_6 = _mm_loadu_pd(&values[101]);
#if defined(__SSE3__) && defined(__AVX__)
c8_6 = _mm_add_pd(c8_6, _mm_mul_pd(a8_6, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_6 = _mm_add_pd(c8_6, _mm_mul_pd(a8_6, b8));
#endif
_mm_storeu_pd(&C[(i*88)+32], c8_6);
__m128d c8_8 = _mm_loadu_pd(&C[(i*88)+53]);
__m128d a8_8 = _mm_loadu_pd(&values[103]);
#if defined(__SSE3__) && defined(__AVX__)
c8_8 = _mm_add_pd(c8_8, _mm_mul_pd(a8_8, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_8 = _mm_add_pd(c8_8, _mm_mul_pd(a8_8, b8));
#endif
_mm_storeu_pd(&C[(i*88)+53], c8_8);
__m128d c8_10 = _mm_loadu_pd(&C[(i*88)+81]);
__m128d a8_10 = _mm_loadu_pd(&values[105]);
#if defined(__SSE3__) && defined(__AVX__)
c8_10 = _mm_add_pd(c8_10, _mm_mul_pd(a8_10, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_10 = _mm_add_pd(c8_10, _mm_mul_pd(a8_10, b8));
#endif
_mm_storeu_pd(&C[(i*88)+81], c8_10);
#else
C[(i*88)+1] += values[95] * B[(i*88)+8];
C[(i*88)+2] += values[96] * B[(i*88)+8];
C[(i*88)+7] += values[97] * B[(i*88)+8];
C[(i*88)+8] += values[98] * B[(i*88)+8];
C[(i*88)+17] += values[99] * B[(i*88)+8];
C[(i*88)+18] += values[100] * B[(i*88)+8];
C[(i*88)+32] += values[101] * B[(i*88)+8];
C[(i*88)+33] += values[102] * B[(i*88)+8];
C[(i*88)+53] += values[103] * B[(i*88)+8];
C[(i*88)+54] += values[104] * B[(i*88)+8];
C[(i*88)+81] += values[105] * B[(i*88)+8];
C[(i*88)+82] += values[106] * B[(i*88)+8];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b9 = _mm256_broadcast_sd(&B[(i*88)+9]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b9 = _mm_loaddup_pd(&B[(i*88)+9]);
#endif
__m128d c9_0 = _mm_load_sd(&C[(i*88)+0]);
__m128d a9_0 = _mm_load_sd(&values[107]);
#if defined(__SSE3__) && defined(__AVX__)
c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, b9));
#endif
_mm_store_sd(&C[(i*88)+0], c9_0);
__m128d c9_1 = _mm_load_sd(&C[(i*88)+3]);
__m128d a9_1 = _mm_load_sd(&values[108]);
#if defined(__SSE3__) && defined(__AVX__)
c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, b9));
#endif
_mm_store_sd(&C[(i*88)+3], c9_1);
__m128d c9_2 = _mm_load_sd(&C[(i*88)+9]);
__m128d a9_2 = _mm_load_sd(&values[109]);
#if defined(__SSE3__) && defined(__AVX__)
c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, b9));
#endif
_mm_store_sd(&C[(i*88)+9], c9_2);
__m128d c9_3 = _mm_load_sd(&C[(i*88)+19]);
__m128d a9_3 = _mm_load_sd(&values[110]);
#if defined(__SSE3__) && defined(__AVX__)
c9_3 = _mm_add_sd(c9_3, _mm_mul_sd(a9_3, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_3 = _mm_add_sd(c9_3, _mm_mul_sd(a9_3, b9));
#endif
_mm_store_sd(&C[(i*88)+19], c9_3);
__m128d c9_4 = _mm_load_sd(&C[(i*88)+34]);
__m128d a9_4 = _mm_load_sd(&values[111]);
#if defined(__SSE3__) && defined(__AVX__)
c9_4 = _mm_add_sd(c9_4, _mm_mul_sd(a9_4, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_4 = _mm_add_sd(c9_4, _mm_mul_sd(a9_4, b9));
#endif
_mm_store_sd(&C[(i*88)+34], c9_4);
__m128d c9_5 = _mm_load_sd(&C[(i*88)+55]);
__m128d a9_5 = _mm_load_sd(&values[112]);
#if defined(__SSE3__) && defined(__AVX__)
c9_5 = _mm_add_sd(c9_5, _mm_mul_sd(a9_5, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_5 = _mm_add_sd(c9_5, _mm_mul_sd(a9_5, b9));
#endif
_mm_store_sd(&C[(i*88)+55], c9_5);
__m128d c9_6 = _mm_load_sd(&C[(i*88)+83]);
__m128d a9_6 = _mm_load_sd(&values[113]);
#if defined(__SSE3__) && defined(__AVX__)
c9_6 = _mm_add_sd(c9_6, _mm_mul_sd(a9_6, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_6 = _mm_add_sd(c9_6, _mm_mul_sd(a9_6, b9));
#endif
_mm_store_sd(&C[(i*88)+83], c9_6);
#else
C[(i*88)+0] += values[107] * B[(i*88)+9];
C[(i*88)+3] += values[108] * B[(i*88)+9];
C[(i*88)+9] += values[109] * B[(i*88)+9];
C[(i*88)+19] += values[110] * B[(i*88)+9];
C[(i*88)+34] += values[111] * B[(i*88)+9];
C[(i*88)+55] += values[112] * B[(i*88)+9];
C[(i*88)+83] += values[113] * B[(i*88)+9];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b10 = _mm256_broadcast_sd(&B[(i*88)+10]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b10 = _mm_loaddup_pd(&B[(i*88)+10]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c10_0 = _mm256_loadu_pd(&C[(i*88)+10]);
__m256d a10_0 = _mm256_loadu_pd(&values[114]);
c10_0 = _mm256_add_pd(c10_0, _mm256_mul_pd(a10_0, b10));
_mm256_storeu_pd(&C[(i*88)+10], c10_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c10_0 = _mm_loadu_pd(&C[(i*88)+10]);
__m128d a10_0 = _mm_loadu_pd(&values[114]);
c10_0 = _mm_add_pd(c10_0, _mm_mul_pd(a10_0, b10));
_mm_storeu_pd(&C[(i*88)+10], c10_0);
__m128d c10_2 = _mm_loadu_pd(&C[(i*88)+12]);
__m128d a10_2 = _mm_loadu_pd(&values[116]);
c10_2 = _mm_add_pd(c10_2, _mm_mul_pd(a10_2, b10));
_mm_storeu_pd(&C[(i*88)+12], c10_2);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c10_4 = _mm256_loadu_pd(&C[(i*88)+25]);
__m256d a10_4 = _mm256_loadu_pd(&values[118]);
c10_4 = _mm256_add_pd(c10_4, _mm256_mul_pd(a10_4, b10));
_mm256_storeu_pd(&C[(i*88)+25], c10_4);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c10_4 = _mm_loadu_pd(&C[(i*88)+25]);
__m128d a10_4 = _mm_loadu_pd(&values[118]);
c10_4 = _mm_add_pd(c10_4, _mm_mul_pd(a10_4, b10));
_mm_storeu_pd(&C[(i*88)+25], c10_4);
__m128d c10_6 = _mm_loadu_pd(&C[(i*88)+27]);
__m128d a10_6 = _mm_loadu_pd(&values[120]);
c10_6 = _mm_add_pd(c10_6, _mm_mul_pd(a10_6, b10));
_mm_storeu_pd(&C[(i*88)+27], c10_6);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c10_8 = _mm256_loadu_pd(&C[(i*88)+46]);
__m256d a10_8 = _mm256_loadu_pd(&values[122]);
c10_8 = _mm256_add_pd(c10_8, _mm256_mul_pd(a10_8, b10));
_mm256_storeu_pd(&C[(i*88)+46], c10_8);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c10_8 = _mm_loadu_pd(&C[(i*88)+46]);
__m128d a10_8 = _mm_loadu_pd(&values[122]);
c10_8 = _mm_add_pd(c10_8, _mm_mul_pd(a10_8, b10));
_mm_storeu_pd(&C[(i*88)+46], c10_8);
__m128d c10_10 = _mm_loadu_pd(&C[(i*88)+48]);
__m128d a10_10 = _mm_loadu_pd(&values[124]);
c10_10 = _mm_add_pd(c10_10, _mm_mul_pd(a10_10, b10));
_mm_storeu_pd(&C[(i*88)+48], c10_10);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c10_12 = _mm256_loadu_pd(&C[(i*88)+74]);
__m256d a10_12 = _mm256_loadu_pd(&values[126]);
c10_12 = _mm256_add_pd(c10_12, _mm256_mul_pd(a10_12, b10));
_mm256_storeu_pd(&C[(i*88)+74], c10_12);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c10_12 = _mm_loadu_pd(&C[(i*88)+74]);
__m128d a10_12 = _mm_loadu_pd(&values[126]);
c10_12 = _mm_add_pd(c10_12, _mm_mul_pd(a10_12, b10));
_mm_storeu_pd(&C[(i*88)+74], c10_12);
__m128d c10_14 = _mm_loadu_pd(&C[(i*88)+76]);
__m128d a10_14 = _mm_loadu_pd(&values[128]);
c10_14 = _mm_add_pd(c10_14, _mm_mul_pd(a10_14, b10));
_mm_storeu_pd(&C[(i*88)+76], c10_14);
#endif
#else
C[(i*88)+10] += values[114] * B[(i*88)+10];
C[(i*88)+11] += values[115] * B[(i*88)+10];
C[(i*88)+12] += values[116] * B[(i*88)+10];
C[(i*88)+13] += values[117] * B[(i*88)+10];
C[(i*88)+25] += values[118] * B[(i*88)+10];
C[(i*88)+26] += values[119] * B[(i*88)+10];
C[(i*88)+27] += values[120] * B[(i*88)+10];
C[(i*88)+28] += values[121] * B[(i*88)+10];
C[(i*88)+46] += values[122] * B[(i*88)+10];
C[(i*88)+47] += values[123] * B[(i*88)+10];
C[(i*88)+48] += values[124] * B[(i*88)+10];
C[(i*88)+49] += values[125] * B[(i*88)+10];
C[(i*88)+74] += values[126] * B[(i*88)+10];
C[(i*88)+75] += values[127] * B[(i*88)+10];
C[(i*88)+76] += values[128] * B[(i*88)+10];
C[(i*88)+77] += values[129] * B[(i*88)+10];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b11 = _mm256_broadcast_sd(&B[(i*88)+11]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b11 = _mm_loaddup_pd(&B[(i*88)+11]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c11_0 = _mm256_loadu_pd(&C[(i*88)+10]);
__m256d a11_0 = _mm256_loadu_pd(&values[130]);
c11_0 = _mm256_add_pd(c11_0, _mm256_mul_pd(a11_0, b11));
_mm256_storeu_pd(&C[(i*88)+10], c11_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c11_0 = _mm_loadu_pd(&C[(i*88)+10]);
__m128d a11_0 = _mm_loadu_pd(&values[130]);
c11_0 = _mm_add_pd(c11_0, _mm_mul_pd(a11_0, b11));
_mm_storeu_pd(&C[(i*88)+10], c11_0);
__m128d c11_2 = _mm_loadu_pd(&C[(i*88)+12]);
__m128d a11_2 = _mm_loadu_pd(&values[132]);
c11_2 = _mm_add_pd(c11_2, _mm_mul_pd(a11_2, b11));
_mm_storeu_pd(&C[(i*88)+12], c11_2);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c11_4 = _mm256_loadu_pd(&C[(i*88)+25]);
__m256d a11_4 = _mm256_loadu_pd(&values[134]);
c11_4 = _mm256_add_pd(c11_4, _mm256_mul_pd(a11_4, b11));
_mm256_storeu_pd(&C[(i*88)+25], c11_4);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c11_4 = _mm_loadu_pd(&C[(i*88)+25]);
__m128d a11_4 = _mm_loadu_pd(&values[134]);
c11_4 = _mm_add_pd(c11_4, _mm_mul_pd(a11_4, b11));
_mm_storeu_pd(&C[(i*88)+25], c11_4);
__m128d c11_6 = _mm_loadu_pd(&C[(i*88)+27]);
__m128d a11_6 = _mm_loadu_pd(&values[136]);
c11_6 = _mm_add_pd(c11_6, _mm_mul_pd(a11_6, b11));
_mm_storeu_pd(&C[(i*88)+27], c11_6);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c11_8 = _mm256_loadu_pd(&C[(i*88)+46]);
__m256d a11_8 = _mm256_loadu_pd(&values[138]);
c11_8 = _mm256_add_pd(c11_8, _mm256_mul_pd(a11_8, b11));
_mm256_storeu_pd(&C[(i*88)+46], c11_8);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c11_8 = _mm_loadu_pd(&C[(i*88)+46]);
__m128d a11_8 = _mm_loadu_pd(&values[138]);
c11_8 = _mm_add_pd(c11_8, _mm_mul_pd(a11_8, b11));
_mm_storeu_pd(&C[(i*88)+46], c11_8);
__m128d c11_10 = _mm_loadu_pd(&C[(i*88)+48]);
__m128d a11_10 = _mm_loadu_pd(&values[140]);
c11_10 = _mm_add_pd(c11_10, _mm_mul_pd(a11_10, b11));
_mm_storeu_pd(&C[(i*88)+48], c11_10);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c11_12 = _mm256_loadu_pd(&C[(i*88)+74]);
__m256d a11_12 = _mm256_loadu_pd(&values[142]);
c11_12 = _mm256_add_pd(c11_12, _mm256_mul_pd(a11_12, b11));
_mm256_storeu_pd(&C[(i*88)+74], c11_12);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c11_12 = _mm_loadu_pd(&C[(i*88)+74]);
__m128d a11_12 = _mm_loadu_pd(&values[142]);
c11_12 = _mm_add_pd(c11_12, _mm_mul_pd(a11_12, b11));
_mm_storeu_pd(&C[(i*88)+74], c11_12);
__m128d c11_14 = _mm_loadu_pd(&C[(i*88)+76]);
__m128d a11_14 = _mm_loadu_pd(&values[144]);
c11_14 = _mm_add_pd(c11_14, _mm_mul_pd(a11_14, b11));
_mm_storeu_pd(&C[(i*88)+76], c11_14);
#endif
#else
C[(i*88)+10] += values[130] * B[(i*88)+11];
C[(i*88)+11] += values[131] * B[(i*88)+11];
C[(i*88)+12] += values[132] * B[(i*88)+11];
C[(i*88)+13] += values[133] * B[(i*88)+11];
C[(i*88)+25] += values[134] * B[(i*88)+11];
C[(i*88)+26] += values[135] * B[(i*88)+11];
C[(i*88)+27] += values[136] * B[(i*88)+11];
C[(i*88)+28] += values[137] * B[(i*88)+11];
C[(i*88)+46] += values[138] * B[(i*88)+11];
C[(i*88)+47] += values[139] * B[(i*88)+11];
C[(i*88)+48] += values[140] * B[(i*88)+11];
C[(i*88)+49] += values[141] * B[(i*88)+11];
C[(i*88)+74] += values[142] * B[(i*88)+11];
C[(i*88)+75] += values[143] * B[(i*88)+11];
C[(i*88)+76] += values[144] * B[(i*88)+11];
C[(i*88)+77] += values[145] * B[(i*88)+11];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b12 = _mm256_broadcast_sd(&B[(i*88)+12]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b12 = _mm_loaddup_pd(&B[(i*88)+12]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c12_0 = _mm256_loadu_pd(&C[(i*88)+10]);
__m256d a12_0 = _mm256_loadu_pd(&values[146]);
c12_0 = _mm256_add_pd(c12_0, _mm256_mul_pd(a12_0, b12));
_mm256_storeu_pd(&C[(i*88)+10], c12_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c12_0 = _mm_loadu_pd(&C[(i*88)+10]);
__m128d a12_0 = _mm_loadu_pd(&values[146]);
c12_0 = _mm_add_pd(c12_0, _mm_mul_pd(a12_0, b12));
_mm_storeu_pd(&C[(i*88)+10], c12_0);
__m128d c12_2 = _mm_loadu_pd(&C[(i*88)+12]);
__m128d a12_2 = _mm_loadu_pd(&values[148]);
c12_2 = _mm_add_pd(c12_2, _mm_mul_pd(a12_2, b12));
_mm_storeu_pd(&C[(i*88)+12], c12_2);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c12_4 = _mm256_loadu_pd(&C[(i*88)+25]);
__m256d a12_4 = _mm256_loadu_pd(&values[150]);
c12_4 = _mm256_add_pd(c12_4, _mm256_mul_pd(a12_4, b12));
_mm256_storeu_pd(&C[(i*88)+25], c12_4);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c12_4 = _mm_loadu_pd(&C[(i*88)+25]);
__m128d a12_4 = _mm_loadu_pd(&values[150]);
c12_4 = _mm_add_pd(c12_4, _mm_mul_pd(a12_4, b12));
_mm_storeu_pd(&C[(i*88)+25], c12_4);
__m128d c12_6 = _mm_loadu_pd(&C[(i*88)+27]);
__m128d a12_6 = _mm_loadu_pd(&values[152]);
c12_6 = _mm_add_pd(c12_6, _mm_mul_pd(a12_6, b12));
_mm_storeu_pd(&C[(i*88)+27], c12_6);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c12_8 = _mm256_loadu_pd(&C[(i*88)+46]);
__m256d a12_8 = _mm256_loadu_pd(&values[154]);
c12_8 = _mm256_add_pd(c12_8, _mm256_mul_pd(a12_8, b12));
_mm256_storeu_pd(&C[(i*88)+46], c12_8);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c12_8 = _mm_loadu_pd(&C[(i*88)+46]);
__m128d a12_8 = _mm_loadu_pd(&values[154]);
c12_8 = _mm_add_pd(c12_8, _mm_mul_pd(a12_8, b12));
_mm_storeu_pd(&C[(i*88)+46], c12_8);
__m128d c12_10 = _mm_loadu_pd(&C[(i*88)+48]);
__m128d a12_10 = _mm_loadu_pd(&values[156]);
c12_10 = _mm_add_pd(c12_10, _mm_mul_pd(a12_10, b12));
_mm_storeu_pd(&C[(i*88)+48], c12_10);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c12_12 = _mm256_loadu_pd(&C[(i*88)+74]);
__m256d a12_12 = _mm256_loadu_pd(&values[158]);
c12_12 = _mm256_add_pd(c12_12, _mm256_mul_pd(a12_12, b12));
_mm256_storeu_pd(&C[(i*88)+74], c12_12);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c12_12 = _mm_loadu_pd(&C[(i*88)+74]);
__m128d a12_12 = _mm_loadu_pd(&values[158]);
c12_12 = _mm_add_pd(c12_12, _mm_mul_pd(a12_12, b12));
_mm_storeu_pd(&C[(i*88)+74], c12_12);
__m128d c12_14 = _mm_loadu_pd(&C[(i*88)+76]);
__m128d a12_14 = _mm_loadu_pd(&values[160]);
c12_14 = _mm_add_pd(c12_14, _mm_mul_pd(a12_14, b12));
_mm_storeu_pd(&C[(i*88)+76], c12_14);
#endif
#else
C[(i*88)+10] += values[146] * B[(i*88)+12];
C[(i*88)+11] += values[147] * B[(i*88)+12];
C[(i*88)+12] += values[148] * B[(i*88)+12];
C[(i*88)+13] += values[149] * B[(i*88)+12];
C[(i*88)+25] += values[150] * B[(i*88)+12];
C[(i*88)+26] += values[151] * B[(i*88)+12];
C[(i*88)+27] += values[152] * B[(i*88)+12];
C[(i*88)+28] += values[153] * B[(i*88)+12];
C[(i*88)+46] += values[154] * B[(i*88)+12];
C[(i*88)+47] += values[155] * B[(i*88)+12];
C[(i*88)+48] += values[156] * B[(i*88)+12];
C[(i*88)+49] += values[157] * B[(i*88)+12];
C[(i*88)+74] += values[158] * B[(i*88)+12];
C[(i*88)+75] += values[159] * B[(i*88)+12];
C[(i*88)+76] += values[160] * B[(i*88)+12];
C[(i*88)+77] += values[161] * B[(i*88)+12];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b13 = _mm256_broadcast_sd(&B[(i*88)+13]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b13 = _mm_loaddup_pd(&B[(i*88)+13]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c13_0 = _mm256_loadu_pd(&C[(i*88)+10]);
__m256d a13_0 = _mm256_loadu_pd(&values[162]);
c13_0 = _mm256_add_pd(c13_0, _mm256_mul_pd(a13_0, b13));
_mm256_storeu_pd(&C[(i*88)+10], c13_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c13_0 = _mm_loadu_pd(&C[(i*88)+10]);
__m128d a13_0 = _mm_loadu_pd(&values[162]);
c13_0 = _mm_add_pd(c13_0, _mm_mul_pd(a13_0, b13));
_mm_storeu_pd(&C[(i*88)+10], c13_0);
__m128d c13_2 = _mm_loadu_pd(&C[(i*88)+12]);
__m128d a13_2 = _mm_loadu_pd(&values[164]);
c13_2 = _mm_add_pd(c13_2, _mm_mul_pd(a13_2, b13));
_mm_storeu_pd(&C[(i*88)+12], c13_2);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c13_4 = _mm256_loadu_pd(&C[(i*88)+25]);
__m256d a13_4 = _mm256_loadu_pd(&values[166]);
c13_4 = _mm256_add_pd(c13_4, _mm256_mul_pd(a13_4, b13));
_mm256_storeu_pd(&C[(i*88)+25], c13_4);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c13_4 = _mm_loadu_pd(&C[(i*88)+25]);
__m128d a13_4 = _mm_loadu_pd(&values[166]);
c13_4 = _mm_add_pd(c13_4, _mm_mul_pd(a13_4, b13));
_mm_storeu_pd(&C[(i*88)+25], c13_4);
__m128d c13_6 = _mm_loadu_pd(&C[(i*88)+27]);
__m128d a13_6 = _mm_loadu_pd(&values[168]);
c13_6 = _mm_add_pd(c13_6, _mm_mul_pd(a13_6, b13));
_mm_storeu_pd(&C[(i*88)+27], c13_6);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c13_8 = _mm256_loadu_pd(&C[(i*88)+46]);
__m256d a13_8 = _mm256_loadu_pd(&values[170]);
c13_8 = _mm256_add_pd(c13_8, _mm256_mul_pd(a13_8, b13));
_mm256_storeu_pd(&C[(i*88)+46], c13_8);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c13_8 = _mm_loadu_pd(&C[(i*88)+46]);
__m128d a13_8 = _mm_loadu_pd(&values[170]);
c13_8 = _mm_add_pd(c13_8, _mm_mul_pd(a13_8, b13));
_mm_storeu_pd(&C[(i*88)+46], c13_8);
__m128d c13_10 = _mm_loadu_pd(&C[(i*88)+48]);
__m128d a13_10 = _mm_loadu_pd(&values[172]);
c13_10 = _mm_add_pd(c13_10, _mm_mul_pd(a13_10, b13));
_mm_storeu_pd(&C[(i*88)+48], c13_10);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c13_12 = _mm256_loadu_pd(&C[(i*88)+74]);
__m256d a13_12 = _mm256_loadu_pd(&values[174]);
c13_12 = _mm256_add_pd(c13_12, _mm256_mul_pd(a13_12, b13));
_mm256_storeu_pd(&C[(i*88)+74], c13_12);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c13_12 = _mm_loadu_pd(&C[(i*88)+74]);
__m128d a13_12 = _mm_loadu_pd(&values[174]);
c13_12 = _mm_add_pd(c13_12, _mm_mul_pd(a13_12, b13));
_mm_storeu_pd(&C[(i*88)+74], c13_12);
__m128d c13_14 = _mm_loadu_pd(&C[(i*88)+76]);
__m128d a13_14 = _mm_loadu_pd(&values[176]);
c13_14 = _mm_add_pd(c13_14, _mm_mul_pd(a13_14, b13));
_mm_storeu_pd(&C[(i*88)+76], c13_14);
#endif
#else
C[(i*88)+10] += values[162] * B[(i*88)+13];
C[(i*88)+11] += values[163] * B[(i*88)+13];
C[(i*88)+12] += values[164] * B[(i*88)+13];
C[(i*88)+13] += values[165] * B[(i*88)+13];
C[(i*88)+25] += values[166] * B[(i*88)+13];
C[(i*88)+26] += values[167] * B[(i*88)+13];
C[(i*88)+27] += values[168] * B[(i*88)+13];
C[(i*88)+28] += values[169] * B[(i*88)+13];
C[(i*88)+46] += values[170] * B[(i*88)+13];
C[(i*88)+47] += values[171] * B[(i*88)+13];
C[(i*88)+48] += values[172] * B[(i*88)+13];
C[(i*88)+49] += values[173] * B[(i*88)+13];
C[(i*88)+74] += values[174] * B[(i*88)+13];
C[(i*88)+75] += values[175] * B[(i*88)+13];
C[(i*88)+76] += values[176] * B[(i*88)+13];
C[(i*88)+77] += values[177] * B[(i*88)+13];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b14 = _mm256_broadcast_sd(&B[(i*88)+14]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b14 = _mm_loaddup_pd(&B[(i*88)+14]);
#endif
__m128d c14_0 = _mm_loadu_pd(&C[(i*88)+4]);
__m128d a14_0 = _mm_loadu_pd(&values[178]);
#if defined(__SSE3__) && defined(__AVX__)
c14_0 = _mm_add_pd(c14_0, _mm_mul_pd(a14_0, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_0 = _mm_add_pd(c14_0, _mm_mul_pd(a14_0, b14));
#endif
_mm_storeu_pd(&C[(i*88)+4], c14_0);
__m128d c14_2 = _mm_load_sd(&C[(i*88)+6]);
__m128d a14_2 = _mm_load_sd(&values[180]);
#if defined(__SSE3__) && defined(__AVX__)
c14_2 = _mm_add_sd(c14_2, _mm_mul_sd(a14_2, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_2 = _mm_add_sd(c14_2, _mm_mul_sd(a14_2, b14));
#endif
_mm_store_sd(&C[(i*88)+6], c14_2);
__m128d c14_3 = _mm_loadu_pd(&C[(i*88)+14]);
__m128d a14_3 = _mm_loadu_pd(&values[181]);
#if defined(__SSE3__) && defined(__AVX__)
c14_3 = _mm_add_pd(c14_3, _mm_mul_pd(a14_3, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_3 = _mm_add_pd(c14_3, _mm_mul_pd(a14_3, b14));
#endif
_mm_storeu_pd(&C[(i*88)+14], c14_3);
__m128d c14_5 = _mm_load_sd(&C[(i*88)+16]);
__m128d a14_5 = _mm_load_sd(&values[183]);
#if defined(__SSE3__) && defined(__AVX__)
c14_5 = _mm_add_sd(c14_5, _mm_mul_sd(a14_5, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_5 = _mm_add_sd(c14_5, _mm_mul_sd(a14_5, b14));
#endif
_mm_store_sd(&C[(i*88)+16], c14_5);
__m128d c14_6 = _mm_loadu_pd(&C[(i*88)+29]);
__m128d a14_6 = _mm_loadu_pd(&values[184]);
#if defined(__SSE3__) && defined(__AVX__)
c14_6 = _mm_add_pd(c14_6, _mm_mul_pd(a14_6, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_6 = _mm_add_pd(c14_6, _mm_mul_pd(a14_6, b14));
#endif
_mm_storeu_pd(&C[(i*88)+29], c14_6);
__m128d c14_8 = _mm_load_sd(&C[(i*88)+31]);
__m128d a14_8 = _mm_load_sd(&values[186]);
#if defined(__SSE3__) && defined(__AVX__)
c14_8 = _mm_add_sd(c14_8, _mm_mul_sd(a14_8, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_8 = _mm_add_sd(c14_8, _mm_mul_sd(a14_8, b14));
#endif
_mm_store_sd(&C[(i*88)+31], c14_8);
__m128d c14_9 = _mm_loadu_pd(&C[(i*88)+50]);
__m128d a14_9 = _mm_loadu_pd(&values[187]);
#if defined(__SSE3__) && defined(__AVX__)
c14_9 = _mm_add_pd(c14_9, _mm_mul_pd(a14_9, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_9 = _mm_add_pd(c14_9, _mm_mul_pd(a14_9, b14));
#endif
_mm_storeu_pd(&C[(i*88)+50], c14_9);
__m128d c14_11 = _mm_load_sd(&C[(i*88)+52]);
__m128d a14_11 = _mm_load_sd(&values[189]);
#if defined(__SSE3__) && defined(__AVX__)
c14_11 = _mm_add_sd(c14_11, _mm_mul_sd(a14_11, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_11 = _mm_add_sd(c14_11, _mm_mul_sd(a14_11, b14));
#endif
_mm_store_sd(&C[(i*88)+52], c14_11);
__m128d c14_12 = _mm_loadu_pd(&C[(i*88)+78]);
__m128d a14_12 = _mm_loadu_pd(&values[190]);
#if defined(__SSE3__) && defined(__AVX__)
c14_12 = _mm_add_pd(c14_12, _mm_mul_pd(a14_12, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_12 = _mm_add_pd(c14_12, _mm_mul_pd(a14_12, b14));
#endif
_mm_storeu_pd(&C[(i*88)+78], c14_12);
__m128d c14_14 = _mm_load_sd(&C[(i*88)+80]);
__m128d a14_14 = _mm_load_sd(&values[192]);
#if defined(__SSE3__) && defined(__AVX__)
c14_14 = _mm_add_sd(c14_14, _mm_mul_sd(a14_14, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_14 = _mm_add_sd(c14_14, _mm_mul_sd(a14_14, b14));
#endif
_mm_store_sd(&C[(i*88)+80], c14_14);
#else
C[(i*88)+4] += values[178] * B[(i*88)+14];
C[(i*88)+5] += values[179] * B[(i*88)+14];
C[(i*88)+6] += values[180] * B[(i*88)+14];
C[(i*88)+14] += values[181] * B[(i*88)+14];
C[(i*88)+15] += values[182] * B[(i*88)+14];
C[(i*88)+16] += values[183] * B[(i*88)+14];
C[(i*88)+29] += values[184] * B[(i*88)+14];
C[(i*88)+30] += values[185] * B[(i*88)+14];
C[(i*88)+31] += values[186] * B[(i*88)+14];
C[(i*88)+50] += values[187] * B[(i*88)+14];
C[(i*88)+51] += values[188] * B[(i*88)+14];
C[(i*88)+52] += values[189] * B[(i*88)+14];
C[(i*88)+78] += values[190] * B[(i*88)+14];
C[(i*88)+79] += values[191] * B[(i*88)+14];
C[(i*88)+80] += values[192] * B[(i*88)+14];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b15 = _mm256_broadcast_sd(&B[(i*88)+15]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b15 = _mm_loaddup_pd(&B[(i*88)+15]);
#endif
__m128d c15_0 = _mm_loadu_pd(&C[(i*88)+4]);
__m128d a15_0 = _mm_loadu_pd(&values[193]);
#if defined(__SSE3__) && defined(__AVX__)
c15_0 = _mm_add_pd(c15_0, _mm_mul_pd(a15_0, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_0 = _mm_add_pd(c15_0, _mm_mul_pd(a15_0, b15));
#endif
_mm_storeu_pd(&C[(i*88)+4], c15_0);
__m128d c15_2 = _mm_load_sd(&C[(i*88)+6]);
__m128d a15_2 = _mm_load_sd(&values[195]);
#if defined(__SSE3__) && defined(__AVX__)
c15_2 = _mm_add_sd(c15_2, _mm_mul_sd(a15_2, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_2 = _mm_add_sd(c15_2, _mm_mul_sd(a15_2, b15));
#endif
_mm_store_sd(&C[(i*88)+6], c15_2);
__m128d c15_3 = _mm_loadu_pd(&C[(i*88)+14]);
__m128d a15_3 = _mm_loadu_pd(&values[196]);
#if defined(__SSE3__) && defined(__AVX__)
c15_3 = _mm_add_pd(c15_3, _mm_mul_pd(a15_3, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_3 = _mm_add_pd(c15_3, _mm_mul_pd(a15_3, b15));
#endif
_mm_storeu_pd(&C[(i*88)+14], c15_3);
__m128d c15_5 = _mm_load_sd(&C[(i*88)+16]);
__m128d a15_5 = _mm_load_sd(&values[198]);
#if defined(__SSE3__) && defined(__AVX__)
c15_5 = _mm_add_sd(c15_5, _mm_mul_sd(a15_5, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_5 = _mm_add_sd(c15_5, _mm_mul_sd(a15_5, b15));
#endif
_mm_store_sd(&C[(i*88)+16], c15_5);
__m128d c15_6 = _mm_loadu_pd(&C[(i*88)+29]);
__m128d a15_6 = _mm_loadu_pd(&values[199]);
#if defined(__SSE3__) && defined(__AVX__)
c15_6 = _mm_add_pd(c15_6, _mm_mul_pd(a15_6, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_6 = _mm_add_pd(c15_6, _mm_mul_pd(a15_6, b15));
#endif
_mm_storeu_pd(&C[(i*88)+29], c15_6);
__m128d c15_8 = _mm_load_sd(&C[(i*88)+31]);
__m128d a15_8 = _mm_load_sd(&values[201]);
#if defined(__SSE3__) && defined(__AVX__)
c15_8 = _mm_add_sd(c15_8, _mm_mul_sd(a15_8, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_8 = _mm_add_sd(c15_8, _mm_mul_sd(a15_8, b15));
#endif
_mm_store_sd(&C[(i*88)+31], c15_8);
__m128d c15_9 = _mm_loadu_pd(&C[(i*88)+50]);
__m128d a15_9 = _mm_loadu_pd(&values[202]);
#if defined(__SSE3__) && defined(__AVX__)
c15_9 = _mm_add_pd(c15_9, _mm_mul_pd(a15_9, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_9 = _mm_add_pd(c15_9, _mm_mul_pd(a15_9, b15));
#endif
_mm_storeu_pd(&C[(i*88)+50], c15_9);
__m128d c15_11 = _mm_load_sd(&C[(i*88)+52]);
__m128d a15_11 = _mm_load_sd(&values[204]);
#if defined(__SSE3__) && defined(__AVX__)
c15_11 = _mm_add_sd(c15_11, _mm_mul_sd(a15_11, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_11 = _mm_add_sd(c15_11, _mm_mul_sd(a15_11, b15));
#endif
_mm_store_sd(&C[(i*88)+52], c15_11);
__m128d c15_12 = _mm_loadu_pd(&C[(i*88)+78]);
__m128d a15_12 = _mm_loadu_pd(&values[205]);
#if defined(__SSE3__) && defined(__AVX__)
c15_12 = _mm_add_pd(c15_12, _mm_mul_pd(a15_12, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_12 = _mm_add_pd(c15_12, _mm_mul_pd(a15_12, b15));
#endif
_mm_storeu_pd(&C[(i*88)+78], c15_12);
__m128d c15_14 = _mm_load_sd(&C[(i*88)+80]);
__m128d a15_14 = _mm_load_sd(&values[207]);
#if defined(__SSE3__) && defined(__AVX__)
c15_14 = _mm_add_sd(c15_14, _mm_mul_sd(a15_14, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_14 = _mm_add_sd(c15_14, _mm_mul_sd(a15_14, b15));
#endif
_mm_store_sd(&C[(i*88)+80], c15_14);
#else
C[(i*88)+4] += values[193] * B[(i*88)+15];
C[(i*88)+5] += values[194] * B[(i*88)+15];
C[(i*88)+6] += values[195] * B[(i*88)+15];
C[(i*88)+14] += values[196] * B[(i*88)+15];
C[(i*88)+15] += values[197] * B[(i*88)+15];
C[(i*88)+16] += values[198] * B[(i*88)+15];
C[(i*88)+29] += values[199] * B[(i*88)+15];
C[(i*88)+30] += values[200] * B[(i*88)+15];
C[(i*88)+31] += values[201] * B[(i*88)+15];
C[(i*88)+50] += values[202] * B[(i*88)+15];
C[(i*88)+51] += values[203] * B[(i*88)+15];
C[(i*88)+52] += values[204] * B[(i*88)+15];
C[(i*88)+78] += values[205] * B[(i*88)+15];
C[(i*88)+79] += values[206] * B[(i*88)+15];
C[(i*88)+80] += values[207] * B[(i*88)+15];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b16 = _mm256_broadcast_sd(&B[(i*88)+16]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b16 = _mm_loaddup_pd(&B[(i*88)+16]);
#endif
__m128d c16_0 = _mm_loadu_pd(&C[(i*88)+4]);
__m128d a16_0 = _mm_loadu_pd(&values[208]);
#if defined(__SSE3__) && defined(__AVX__)
c16_0 = _mm_add_pd(c16_0, _mm_mul_pd(a16_0, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_0 = _mm_add_pd(c16_0, _mm_mul_pd(a16_0, b16));
#endif
_mm_storeu_pd(&C[(i*88)+4], c16_0);
__m128d c16_2 = _mm_load_sd(&C[(i*88)+6]);
__m128d a16_2 = _mm_load_sd(&values[210]);
#if defined(__SSE3__) && defined(__AVX__)
c16_2 = _mm_add_sd(c16_2, _mm_mul_sd(a16_2, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_2 = _mm_add_sd(c16_2, _mm_mul_sd(a16_2, b16));
#endif
_mm_store_sd(&C[(i*88)+6], c16_2);
__m128d c16_3 = _mm_loadu_pd(&C[(i*88)+14]);
__m128d a16_3 = _mm_loadu_pd(&values[211]);
#if defined(__SSE3__) && defined(__AVX__)
c16_3 = _mm_add_pd(c16_3, _mm_mul_pd(a16_3, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_3 = _mm_add_pd(c16_3, _mm_mul_pd(a16_3, b16));
#endif
_mm_storeu_pd(&C[(i*88)+14], c16_3);
__m128d c16_5 = _mm_load_sd(&C[(i*88)+16]);
__m128d a16_5 = _mm_load_sd(&values[213]);
#if defined(__SSE3__) && defined(__AVX__)
c16_5 = _mm_add_sd(c16_5, _mm_mul_sd(a16_5, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_5 = _mm_add_sd(c16_5, _mm_mul_sd(a16_5, b16));
#endif
_mm_store_sd(&C[(i*88)+16], c16_5);
__m128d c16_6 = _mm_loadu_pd(&C[(i*88)+29]);
__m128d a16_6 = _mm_loadu_pd(&values[214]);
#if defined(__SSE3__) && defined(__AVX__)
c16_6 = _mm_add_pd(c16_6, _mm_mul_pd(a16_6, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_6 = _mm_add_pd(c16_6, _mm_mul_pd(a16_6, b16));
#endif
_mm_storeu_pd(&C[(i*88)+29], c16_6);
__m128d c16_8 = _mm_load_sd(&C[(i*88)+31]);
__m128d a16_8 = _mm_load_sd(&values[216]);
#if defined(__SSE3__) && defined(__AVX__)
c16_8 = _mm_add_sd(c16_8, _mm_mul_sd(a16_8, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_8 = _mm_add_sd(c16_8, _mm_mul_sd(a16_8, b16));
#endif
_mm_store_sd(&C[(i*88)+31], c16_8);
__m128d c16_9 = _mm_loadu_pd(&C[(i*88)+50]);
__m128d a16_9 = _mm_loadu_pd(&values[217]);
#if defined(__SSE3__) && defined(__AVX__)
c16_9 = _mm_add_pd(c16_9, _mm_mul_pd(a16_9, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_9 = _mm_add_pd(c16_9, _mm_mul_pd(a16_9, b16));
#endif
_mm_storeu_pd(&C[(i*88)+50], c16_9);
__m128d c16_11 = _mm_load_sd(&C[(i*88)+52]);
__m128d a16_11 = _mm_load_sd(&values[219]);
#if defined(__SSE3__) && defined(__AVX__)
c16_11 = _mm_add_sd(c16_11, _mm_mul_sd(a16_11, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_11 = _mm_add_sd(c16_11, _mm_mul_sd(a16_11, b16));
#endif
_mm_store_sd(&C[(i*88)+52], c16_11);
__m128d c16_12 = _mm_loadu_pd(&C[(i*88)+78]);
__m128d a16_12 = _mm_loadu_pd(&values[220]);
#if defined(__SSE3__) && defined(__AVX__)
c16_12 = _mm_add_pd(c16_12, _mm_mul_pd(a16_12, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_12 = _mm_add_pd(c16_12, _mm_mul_pd(a16_12, b16));
#endif
_mm_storeu_pd(&C[(i*88)+78], c16_12);
__m128d c16_14 = _mm_load_sd(&C[(i*88)+80]);
__m128d a16_14 = _mm_load_sd(&values[222]);
#if defined(__SSE3__) && defined(__AVX__)
c16_14 = _mm_add_sd(c16_14, _mm_mul_sd(a16_14, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_14 = _mm_add_sd(c16_14, _mm_mul_sd(a16_14, b16));
#endif
_mm_store_sd(&C[(i*88)+80], c16_14);
#else
C[(i*88)+4] += values[208] * B[(i*88)+16];
C[(i*88)+5] += values[209] * B[(i*88)+16];
C[(i*88)+6] += values[210] * B[(i*88)+16];
C[(i*88)+14] += values[211] * B[(i*88)+16];
C[(i*88)+15] += values[212] * B[(i*88)+16];
C[(i*88)+16] += values[213] * B[(i*88)+16];
C[(i*88)+29] += values[214] * B[(i*88)+16];
C[(i*88)+30] += values[215] * B[(i*88)+16];
C[(i*88)+31] += values[216] * B[(i*88)+16];
C[(i*88)+50] += values[217] * B[(i*88)+16];
C[(i*88)+51] += values[218] * B[(i*88)+16];
C[(i*88)+52] += values[219] * B[(i*88)+16];
C[(i*88)+78] += values[220] * B[(i*88)+16];
C[(i*88)+79] += values[221] * B[(i*88)+16];
C[(i*88)+80] += values[222] * B[(i*88)+16];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b17 = _mm256_broadcast_sd(&B[(i*88)+17]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b17 = _mm_loaddup_pd(&B[(i*88)+17]);
#endif
__m128d c17_0 = _mm_loadu_pd(&C[(i*88)+1]);
__m128d a17_0 = _mm_loadu_pd(&values[223]);
#if defined(__SSE3__) && defined(__AVX__)
c17_0 = _mm_add_pd(c17_0, _mm_mul_pd(a17_0, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_0 = _mm_add_pd(c17_0, _mm_mul_pd(a17_0, b17));
#endif
_mm_storeu_pd(&C[(i*88)+1], c17_0);
__m128d c17_2 = _mm_loadu_pd(&C[(i*88)+7]);
__m128d a17_2 = _mm_loadu_pd(&values[225]);
#if defined(__SSE3__) && defined(__AVX__)
c17_2 = _mm_add_pd(c17_2, _mm_mul_pd(a17_2, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_2 = _mm_add_pd(c17_2, _mm_mul_pd(a17_2, b17));
#endif
_mm_storeu_pd(&C[(i*88)+7], c17_2);
__m128d c17_4 = _mm_loadu_pd(&C[(i*88)+17]);
__m128d a17_4 = _mm_loadu_pd(&values[227]);
#if defined(__SSE3__) && defined(__AVX__)
c17_4 = _mm_add_pd(c17_4, _mm_mul_pd(a17_4, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_4 = _mm_add_pd(c17_4, _mm_mul_pd(a17_4, b17));
#endif
_mm_storeu_pd(&C[(i*88)+17], c17_4);
__m128d c17_6 = _mm_loadu_pd(&C[(i*88)+32]);
__m128d a17_6 = _mm_loadu_pd(&values[229]);
#if defined(__SSE3__) && defined(__AVX__)
c17_6 = _mm_add_pd(c17_6, _mm_mul_pd(a17_6, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_6 = _mm_add_pd(c17_6, _mm_mul_pd(a17_6, b17));
#endif
_mm_storeu_pd(&C[(i*88)+32], c17_6);
__m128d c17_8 = _mm_loadu_pd(&C[(i*88)+53]);
__m128d a17_8 = _mm_loadu_pd(&values[231]);
#if defined(__SSE3__) && defined(__AVX__)
c17_8 = _mm_add_pd(c17_8, _mm_mul_pd(a17_8, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_8 = _mm_add_pd(c17_8, _mm_mul_pd(a17_8, b17));
#endif
_mm_storeu_pd(&C[(i*88)+53], c17_8);
__m128d c17_10 = _mm_loadu_pd(&C[(i*88)+81]);
__m128d a17_10 = _mm_loadu_pd(&values[233]);
#if defined(__SSE3__) && defined(__AVX__)
c17_10 = _mm_add_pd(c17_10, _mm_mul_pd(a17_10, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_10 = _mm_add_pd(c17_10, _mm_mul_pd(a17_10, b17));
#endif
_mm_storeu_pd(&C[(i*88)+81], c17_10);
#else
C[(i*88)+1] += values[223] * B[(i*88)+17];
C[(i*88)+2] += values[224] * B[(i*88)+17];
C[(i*88)+7] += values[225] * B[(i*88)+17];
C[(i*88)+8] += values[226] * B[(i*88)+17];
C[(i*88)+17] += values[227] * B[(i*88)+17];
C[(i*88)+18] += values[228] * B[(i*88)+17];
C[(i*88)+32] += values[229] * B[(i*88)+17];
C[(i*88)+33] += values[230] * B[(i*88)+17];
C[(i*88)+53] += values[231] * B[(i*88)+17];
C[(i*88)+54] += values[232] * B[(i*88)+17];
C[(i*88)+81] += values[233] * B[(i*88)+17];
C[(i*88)+82] += values[234] * B[(i*88)+17];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b18 = _mm256_broadcast_sd(&B[(i*88)+18]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b18 = _mm_loaddup_pd(&B[(i*88)+18]);
#endif
__m128d c18_0 = _mm_loadu_pd(&C[(i*88)+1]);
__m128d a18_0 = _mm_loadu_pd(&values[235]);
#if defined(__SSE3__) && defined(__AVX__)
c18_0 = _mm_add_pd(c18_0, _mm_mul_pd(a18_0, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_0 = _mm_add_pd(c18_0, _mm_mul_pd(a18_0, b18));
#endif
_mm_storeu_pd(&C[(i*88)+1], c18_0);
__m128d c18_2 = _mm_loadu_pd(&C[(i*88)+7]);
__m128d a18_2 = _mm_loadu_pd(&values[237]);
#if defined(__SSE3__) && defined(__AVX__)
c18_2 = _mm_add_pd(c18_2, _mm_mul_pd(a18_2, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_2 = _mm_add_pd(c18_2, _mm_mul_pd(a18_2, b18));
#endif
_mm_storeu_pd(&C[(i*88)+7], c18_2);
__m128d c18_4 = _mm_loadu_pd(&C[(i*88)+17]);
__m128d a18_4 = _mm_loadu_pd(&values[239]);
#if defined(__SSE3__) && defined(__AVX__)
c18_4 = _mm_add_pd(c18_4, _mm_mul_pd(a18_4, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_4 = _mm_add_pd(c18_4, _mm_mul_pd(a18_4, b18));
#endif
_mm_storeu_pd(&C[(i*88)+17], c18_4);
__m128d c18_6 = _mm_loadu_pd(&C[(i*88)+32]);
__m128d a18_6 = _mm_loadu_pd(&values[241]);
#if defined(__SSE3__) && defined(__AVX__)
c18_6 = _mm_add_pd(c18_6, _mm_mul_pd(a18_6, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_6 = _mm_add_pd(c18_6, _mm_mul_pd(a18_6, b18));
#endif
_mm_storeu_pd(&C[(i*88)+32], c18_6);
__m128d c18_8 = _mm_loadu_pd(&C[(i*88)+53]);
__m128d a18_8 = _mm_loadu_pd(&values[243]);
#if defined(__SSE3__) && defined(__AVX__)
c18_8 = _mm_add_pd(c18_8, _mm_mul_pd(a18_8, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_8 = _mm_add_pd(c18_8, _mm_mul_pd(a18_8, b18));
#endif
_mm_storeu_pd(&C[(i*88)+53], c18_8);
__m128d c18_10 = _mm_loadu_pd(&C[(i*88)+81]);
__m128d a18_10 = _mm_loadu_pd(&values[245]);
#if defined(__SSE3__) && defined(__AVX__)
c18_10 = _mm_add_pd(c18_10, _mm_mul_pd(a18_10, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_10 = _mm_add_pd(c18_10, _mm_mul_pd(a18_10, b18));
#endif
_mm_storeu_pd(&C[(i*88)+81], c18_10);
#else
C[(i*88)+1] += values[235] * B[(i*88)+18];
C[(i*88)+2] += values[236] * B[(i*88)+18];
C[(i*88)+7] += values[237] * B[(i*88)+18];
C[(i*88)+8] += values[238] * B[(i*88)+18];
C[(i*88)+17] += values[239] * B[(i*88)+18];
C[(i*88)+18] += values[240] * B[(i*88)+18];
C[(i*88)+32] += values[241] * B[(i*88)+18];
C[(i*88)+33] += values[242] * B[(i*88)+18];
C[(i*88)+53] += values[243] * B[(i*88)+18];
C[(i*88)+54] += values[244] * B[(i*88)+18];
C[(i*88)+81] += values[245] * B[(i*88)+18];
C[(i*88)+82] += values[246] * B[(i*88)+18];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b19 = _mm256_broadcast_sd(&B[(i*88)+19]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b19 = _mm_loaddup_pd(&B[(i*88)+19]);
#endif
__m128d c19_0 = _mm_load_sd(&C[(i*88)+0]);
__m128d a19_0 = _mm_load_sd(&values[247]);
#if defined(__SSE3__) && defined(__AVX__)
c19_0 = _mm_add_sd(c19_0, _mm_mul_sd(a19_0, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_0 = _mm_add_sd(c19_0, _mm_mul_sd(a19_0, b19));
#endif
_mm_store_sd(&C[(i*88)+0], c19_0);
__m128d c19_1 = _mm_load_sd(&C[(i*88)+3]);
__m128d a19_1 = _mm_load_sd(&values[248]);
#if defined(__SSE3__) && defined(__AVX__)
c19_1 = _mm_add_sd(c19_1, _mm_mul_sd(a19_1, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_1 = _mm_add_sd(c19_1, _mm_mul_sd(a19_1, b19));
#endif
_mm_store_sd(&C[(i*88)+3], c19_1);
__m128d c19_2 = _mm_load_sd(&C[(i*88)+9]);
__m128d a19_2 = _mm_load_sd(&values[249]);
#if defined(__SSE3__) && defined(__AVX__)
c19_2 = _mm_add_sd(c19_2, _mm_mul_sd(a19_2, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_2 = _mm_add_sd(c19_2, _mm_mul_sd(a19_2, b19));
#endif
_mm_store_sd(&C[(i*88)+9], c19_2);
__m128d c19_3 = _mm_load_sd(&C[(i*88)+19]);
__m128d a19_3 = _mm_load_sd(&values[250]);
#if defined(__SSE3__) && defined(__AVX__)
c19_3 = _mm_add_sd(c19_3, _mm_mul_sd(a19_3, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_3 = _mm_add_sd(c19_3, _mm_mul_sd(a19_3, b19));
#endif
_mm_store_sd(&C[(i*88)+19], c19_3);
__m128d c19_4 = _mm_load_sd(&C[(i*88)+34]);
__m128d a19_4 = _mm_load_sd(&values[251]);
#if defined(__SSE3__) && defined(__AVX__)
c19_4 = _mm_add_sd(c19_4, _mm_mul_sd(a19_4, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_4 = _mm_add_sd(c19_4, _mm_mul_sd(a19_4, b19));
#endif
_mm_store_sd(&C[(i*88)+34], c19_4);
__m128d c19_5 = _mm_load_sd(&C[(i*88)+55]);
__m128d a19_5 = _mm_load_sd(&values[252]);
#if defined(__SSE3__) && defined(__AVX__)
c19_5 = _mm_add_sd(c19_5, _mm_mul_sd(a19_5, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_5 = _mm_add_sd(c19_5, _mm_mul_sd(a19_5, b19));
#endif
_mm_store_sd(&C[(i*88)+55], c19_5);
__m128d c19_6 = _mm_load_sd(&C[(i*88)+83]);
__m128d a19_6 = _mm_load_sd(&values[253]);
#if defined(__SSE3__) && defined(__AVX__)
c19_6 = _mm_add_sd(c19_6, _mm_mul_sd(a19_6, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_6 = _mm_add_sd(c19_6, _mm_mul_sd(a19_6, b19));
#endif
_mm_store_sd(&C[(i*88)+83], c19_6);
#else
C[(i*88)+0] += values[247] * B[(i*88)+19];
C[(i*88)+3] += values[248] * B[(i*88)+19];
C[(i*88)+9] += values[249] * B[(i*88)+19];
C[(i*88)+19] += values[250] * B[(i*88)+19];
C[(i*88)+34] += values[251] * B[(i*88)+19];
C[(i*88)+55] += values[252] * B[(i*88)+19];
C[(i*88)+83] += values[253] * B[(i*88)+19];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b20 = _mm256_broadcast_sd(&B[(i*88)+20]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b20 = _mm_loaddup_pd(&B[(i*88)+20]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c20_0 = _mm256_loadu_pd(&C[(i*88)+20]);
__m256d a20_0 = _mm256_loadu_pd(&values[254]);
c20_0 = _mm256_add_pd(c20_0, _mm256_mul_pd(a20_0, b20));
_mm256_storeu_pd(&C[(i*88)+20], c20_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c20_0 = _mm_loadu_pd(&C[(i*88)+20]);
__m128d a20_0 = _mm_loadu_pd(&values[254]);
c20_0 = _mm_add_pd(c20_0, _mm_mul_pd(a20_0, b20));
_mm_storeu_pd(&C[(i*88)+20], c20_0);
__m128d c20_2 = _mm_loadu_pd(&C[(i*88)+22]);
__m128d a20_2 = _mm_loadu_pd(&values[256]);
c20_2 = _mm_add_pd(c20_2, _mm_mul_pd(a20_2, b20));
_mm_storeu_pd(&C[(i*88)+22], c20_2);
#endif
__m128d c20_4 = _mm_load_sd(&C[(i*88)+24]);
__m128d a20_4 = _mm_load_sd(&values[258]);
#if defined(__SSE3__) && defined(__AVX__)
c20_4 = _mm_add_sd(c20_4, _mm_mul_sd(a20_4, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c20_4 = _mm_add_sd(c20_4, _mm_mul_sd(a20_4, b20));
#endif
_mm_store_sd(&C[(i*88)+24], c20_4);
#if defined(__SSE3__) && defined(__AVX__)
__m256d c20_5 = _mm256_loadu_pd(&C[(i*88)+41]);
__m256d a20_5 = _mm256_loadu_pd(&values[259]);
c20_5 = _mm256_add_pd(c20_5, _mm256_mul_pd(a20_5, b20));
_mm256_storeu_pd(&C[(i*88)+41], c20_5);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c20_5 = _mm_loadu_pd(&C[(i*88)+41]);
__m128d a20_5 = _mm_loadu_pd(&values[259]);
c20_5 = _mm_add_pd(c20_5, _mm_mul_pd(a20_5, b20));
_mm_storeu_pd(&C[(i*88)+41], c20_5);
__m128d c20_7 = _mm_loadu_pd(&C[(i*88)+43]);
__m128d a20_7 = _mm_loadu_pd(&values[261]);
c20_7 = _mm_add_pd(c20_7, _mm_mul_pd(a20_7, b20));
_mm_storeu_pd(&C[(i*88)+43], c20_7);
#endif
__m128d c20_9 = _mm_load_sd(&C[(i*88)+45]);
__m128d a20_9 = _mm_load_sd(&values[263]);
#if defined(__SSE3__) && defined(__AVX__)
c20_9 = _mm_add_sd(c20_9, _mm_mul_sd(a20_9, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c20_9 = _mm_add_sd(c20_9, _mm_mul_sd(a20_9, b20));
#endif
_mm_store_sd(&C[(i*88)+45], c20_9);
#if defined(__SSE3__) && defined(__AVX__)
__m256d c20_10 = _mm256_loadu_pd(&C[(i*88)+69]);
__m256d a20_10 = _mm256_loadu_pd(&values[264]);
c20_10 = _mm256_add_pd(c20_10, _mm256_mul_pd(a20_10, b20));
_mm256_storeu_pd(&C[(i*88)+69], c20_10);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c20_10 = _mm_loadu_pd(&C[(i*88)+69]);
__m128d a20_10 = _mm_loadu_pd(&values[264]);
c20_10 = _mm_add_pd(c20_10, _mm_mul_pd(a20_10, b20));
_mm_storeu_pd(&C[(i*88)+69], c20_10);
__m128d c20_12 = _mm_loadu_pd(&C[(i*88)+71]);
__m128d a20_12 = _mm_loadu_pd(&values[266]);
c20_12 = _mm_add_pd(c20_12, _mm_mul_pd(a20_12, b20));
_mm_storeu_pd(&C[(i*88)+71], c20_12);
#endif
__m128d c20_14 = _mm_load_sd(&C[(i*88)+73]);
__m128d a20_14 = _mm_load_sd(&values[268]);
#if defined(__SSE3__) && defined(__AVX__)
c20_14 = _mm_add_sd(c20_14, _mm_mul_sd(a20_14, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c20_14 = _mm_add_sd(c20_14, _mm_mul_sd(a20_14, b20));
#endif
_mm_store_sd(&C[(i*88)+73], c20_14);
#else
C[(i*88)+20] += values[254] * B[(i*88)+20];
C[(i*88)+21] += values[255] * B[(i*88)+20];
C[(i*88)+22] += values[256] * B[(i*88)+20];
C[(i*88)+23] += values[257] * B[(i*88)+20];
C[(i*88)+24] += values[258] * B[(i*88)+20];
C[(i*88)+41] += values[259] * B[(i*88)+20];
C[(i*88)+42] += values[260] * B[(i*88)+20];
C[(i*88)+43] += values[261] * B[(i*88)+20];
C[(i*88)+44] += values[262] * B[(i*88)+20];
C[(i*88)+45] += values[263] * B[(i*88)+20];
C[(i*88)+69] += values[264] * B[(i*88)+20];
C[(i*88)+70] += values[265] * B[(i*88)+20];
C[(i*88)+71] += values[266] * B[(i*88)+20];
C[(i*88)+72] += values[267] * B[(i*88)+20];
C[(i*88)+73] += values[268] * B[(i*88)+20];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b21 = _mm256_broadcast_sd(&B[(i*88)+21]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b21 = _mm_loaddup_pd(&B[(i*88)+21]);
#endif
__m128d c21_0 = _mm_loadu_pd(&C[(i*88)+20]);
__m128d a21_0 = _mm_loadu_pd(&values[269]);
#if defined(__SSE3__) && defined(__AVX__)
c21_0 = _mm_add_pd(c21_0, _mm_mul_pd(a21_0, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c21_0 = _mm_add_pd(c21_0, _mm_mul_pd(a21_0, b21));
#endif
_mm_storeu_pd(&C[(i*88)+20], c21_0);
__m128d c21_2 = _mm_load_sd(&C[(i*88)+22]);
__m128d a21_2 = _mm_load_sd(&values[271]);
#if defined(__SSE3__) && defined(__AVX__)
c21_2 = _mm_add_sd(c21_2, _mm_mul_sd(a21_2, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c21_2 = _mm_add_sd(c21_2, _mm_mul_sd(a21_2, b21));
#endif
_mm_store_sd(&C[(i*88)+22], c21_2);
__m128d c21_3 = _mm_load_sd(&C[(i*88)+24]);
__m128d a21_3 = _mm_load_sd(&values[272]);
#if defined(__SSE3__) && defined(__AVX__)
c21_3 = _mm_add_sd(c21_3, _mm_mul_sd(a21_3, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c21_3 = _mm_add_sd(c21_3, _mm_mul_sd(a21_3, b21));
#endif
_mm_store_sd(&C[(i*88)+24], c21_3);
__m128d c21_4 = _mm_loadu_pd(&C[(i*88)+41]);
__m128d a21_4 = _mm_loadu_pd(&values[273]);
#if defined(__SSE3__) && defined(__AVX__)
c21_4 = _mm_add_pd(c21_4, _mm_mul_pd(a21_4, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c21_4 = _mm_add_pd(c21_4, _mm_mul_pd(a21_4, b21));
#endif
_mm_storeu_pd(&C[(i*88)+41], c21_4);
__m128d c21_6 = _mm_load_sd(&C[(i*88)+43]);
__m128d a21_6 = _mm_load_sd(&values[275]);
#if defined(__SSE3__) && defined(__AVX__)
c21_6 = _mm_add_sd(c21_6, _mm_mul_sd(a21_6, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c21_6 = _mm_add_sd(c21_6, _mm_mul_sd(a21_6, b21));
#endif
_mm_store_sd(&C[(i*88)+43], c21_6);
__m128d c21_7 = _mm_load_sd(&C[(i*88)+45]);
__m128d a21_7 = _mm_load_sd(&values[276]);
#if defined(__SSE3__) && defined(__AVX__)
c21_7 = _mm_add_sd(c21_7, _mm_mul_sd(a21_7, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c21_7 = _mm_add_sd(c21_7, _mm_mul_sd(a21_7, b21));
#endif
_mm_store_sd(&C[(i*88)+45], c21_7);
__m128d c21_8 = _mm_loadu_pd(&C[(i*88)+69]);
__m128d a21_8 = _mm_loadu_pd(&values[277]);
#if defined(__SSE3__) && defined(__AVX__)
c21_8 = _mm_add_pd(c21_8, _mm_mul_pd(a21_8, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c21_8 = _mm_add_pd(c21_8, _mm_mul_pd(a21_8, b21));
#endif
_mm_storeu_pd(&C[(i*88)+69], c21_8);
__m128d c21_10 = _mm_load_sd(&C[(i*88)+71]);
__m128d a21_10 = _mm_load_sd(&values[279]);
#if defined(__SSE3__) && defined(__AVX__)
c21_10 = _mm_add_sd(c21_10, _mm_mul_sd(a21_10, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c21_10 = _mm_add_sd(c21_10, _mm_mul_sd(a21_10, b21));
#endif
_mm_store_sd(&C[(i*88)+71], c21_10);
__m128d c21_11 = _mm_load_sd(&C[(i*88)+73]);
__m128d a21_11 = _mm_load_sd(&values[280]);
#if defined(__SSE3__) && defined(__AVX__)
c21_11 = _mm_add_sd(c21_11, _mm_mul_sd(a21_11, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c21_11 = _mm_add_sd(c21_11, _mm_mul_sd(a21_11, b21));
#endif
_mm_store_sd(&C[(i*88)+73], c21_11);
#else
C[(i*88)+20] += values[269] * B[(i*88)+21];
C[(i*88)+21] += values[270] * B[(i*88)+21];
C[(i*88)+22] += values[271] * B[(i*88)+21];
C[(i*88)+24] += values[272] * B[(i*88)+21];
C[(i*88)+41] += values[273] * B[(i*88)+21];
C[(i*88)+42] += values[274] * B[(i*88)+21];
C[(i*88)+43] += values[275] * B[(i*88)+21];
C[(i*88)+45] += values[276] * B[(i*88)+21];
C[(i*88)+69] += values[277] * B[(i*88)+21];
C[(i*88)+70] += values[278] * B[(i*88)+21];
C[(i*88)+71] += values[279] * B[(i*88)+21];
C[(i*88)+73] += values[280] * B[(i*88)+21];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b22 = _mm256_broadcast_sd(&B[(i*88)+22]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b22 = _mm_loaddup_pd(&B[(i*88)+22]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c22_0 = _mm256_loadu_pd(&C[(i*88)+20]);
__m256d a22_0 = _mm256_loadu_pd(&values[281]);
c22_0 = _mm256_add_pd(c22_0, _mm256_mul_pd(a22_0, b22));
_mm256_storeu_pd(&C[(i*88)+20], c22_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c22_0 = _mm_loadu_pd(&C[(i*88)+20]);
__m128d a22_0 = _mm_loadu_pd(&values[281]);
c22_0 = _mm_add_pd(c22_0, _mm_mul_pd(a22_0, b22));
_mm_storeu_pd(&C[(i*88)+20], c22_0);
__m128d c22_2 = _mm_loadu_pd(&C[(i*88)+22]);
__m128d a22_2 = _mm_loadu_pd(&values[283]);
c22_2 = _mm_add_pd(c22_2, _mm_mul_pd(a22_2, b22));
_mm_storeu_pd(&C[(i*88)+22], c22_2);
#endif
__m128d c22_4 = _mm_load_sd(&C[(i*88)+24]);
__m128d a22_4 = _mm_load_sd(&values[285]);
#if defined(__SSE3__) && defined(__AVX__)
c22_4 = _mm_add_sd(c22_4, _mm_mul_sd(a22_4, _mm256_castpd256_pd128(b22)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c22_4 = _mm_add_sd(c22_4, _mm_mul_sd(a22_4, b22));
#endif
_mm_store_sd(&C[(i*88)+24], c22_4);
#if defined(__SSE3__) && defined(__AVX__)
__m256d c22_5 = _mm256_loadu_pd(&C[(i*88)+41]);
__m256d a22_5 = _mm256_loadu_pd(&values[286]);
c22_5 = _mm256_add_pd(c22_5, _mm256_mul_pd(a22_5, b22));
_mm256_storeu_pd(&C[(i*88)+41], c22_5);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c22_5 = _mm_loadu_pd(&C[(i*88)+41]);
__m128d a22_5 = _mm_loadu_pd(&values[286]);
c22_5 = _mm_add_pd(c22_5, _mm_mul_pd(a22_5, b22));
_mm_storeu_pd(&C[(i*88)+41], c22_5);
__m128d c22_7 = _mm_loadu_pd(&C[(i*88)+43]);
__m128d a22_7 = _mm_loadu_pd(&values[288]);
c22_7 = _mm_add_pd(c22_7, _mm_mul_pd(a22_7, b22));
_mm_storeu_pd(&C[(i*88)+43], c22_7);
#endif
__m128d c22_9 = _mm_load_sd(&C[(i*88)+45]);
__m128d a22_9 = _mm_load_sd(&values[290]);
#if defined(__SSE3__) && defined(__AVX__)
c22_9 = _mm_add_sd(c22_9, _mm_mul_sd(a22_9, _mm256_castpd256_pd128(b22)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c22_9 = _mm_add_sd(c22_9, _mm_mul_sd(a22_9, b22));
#endif
_mm_store_sd(&C[(i*88)+45], c22_9);
#if defined(__SSE3__) && defined(__AVX__)
__m256d c22_10 = _mm256_loadu_pd(&C[(i*88)+69]);
__m256d a22_10 = _mm256_loadu_pd(&values[291]);
c22_10 = _mm256_add_pd(c22_10, _mm256_mul_pd(a22_10, b22));
_mm256_storeu_pd(&C[(i*88)+69], c22_10);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c22_10 = _mm_loadu_pd(&C[(i*88)+69]);
__m128d a22_10 = _mm_loadu_pd(&values[291]);
c22_10 = _mm_add_pd(c22_10, _mm_mul_pd(a22_10, b22));
_mm_storeu_pd(&C[(i*88)+69], c22_10);
__m128d c22_12 = _mm_loadu_pd(&C[(i*88)+71]);
__m128d a22_12 = _mm_loadu_pd(&values[293]);
c22_12 = _mm_add_pd(c22_12, _mm_mul_pd(a22_12, b22));
_mm_storeu_pd(&C[(i*88)+71], c22_12);
#endif
__m128d c22_14 = _mm_load_sd(&C[(i*88)+73]);
__m128d a22_14 = _mm_load_sd(&values[295]);
#if defined(__SSE3__) && defined(__AVX__)
c22_14 = _mm_add_sd(c22_14, _mm_mul_sd(a22_14, _mm256_castpd256_pd128(b22)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c22_14 = _mm_add_sd(c22_14, _mm_mul_sd(a22_14, b22));
#endif
_mm_store_sd(&C[(i*88)+73], c22_14);
#else
C[(i*88)+20] += values[281] * B[(i*88)+22];
C[(i*88)+21] += values[282] * B[(i*88)+22];
C[(i*88)+22] += values[283] * B[(i*88)+22];
C[(i*88)+23] += values[284] * B[(i*88)+22];
C[(i*88)+24] += values[285] * B[(i*88)+22];
C[(i*88)+41] += values[286] * B[(i*88)+22];
C[(i*88)+42] += values[287] * B[(i*88)+22];
C[(i*88)+43] += values[288] * B[(i*88)+22];
C[(i*88)+44] += values[289] * B[(i*88)+22];
C[(i*88)+45] += values[290] * B[(i*88)+22];
C[(i*88)+69] += values[291] * B[(i*88)+22];
C[(i*88)+70] += values[292] * B[(i*88)+22];
C[(i*88)+71] += values[293] * B[(i*88)+22];
C[(i*88)+72] += values[294] * B[(i*88)+22];
C[(i*88)+73] += values[295] * B[(i*88)+22];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b23 = _mm256_broadcast_sd(&B[(i*88)+23]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b23 = _mm_loaddup_pd(&B[(i*88)+23]);
#endif
__m128d c23_0 = _mm_load_sd(&C[(i*88)+20]);
__m128d a23_0 = _mm_load_sd(&values[296]);
#if defined(__SSE3__) && defined(__AVX__)
c23_0 = _mm_add_sd(c23_0, _mm_mul_sd(a23_0, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c23_0 = _mm_add_sd(c23_0, _mm_mul_sd(a23_0, b23));
#endif
_mm_store_sd(&C[(i*88)+20], c23_0);
__m128d c23_1 = _mm_loadu_pd(&C[(i*88)+22]);
__m128d a23_1 = _mm_loadu_pd(&values[297]);
#if defined(__SSE3__) && defined(__AVX__)
c23_1 = _mm_add_pd(c23_1, _mm_mul_pd(a23_1, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c23_1 = _mm_add_pd(c23_1, _mm_mul_pd(a23_1, b23));
#endif
_mm_storeu_pd(&C[(i*88)+22], c23_1);
__m128d c23_3 = _mm_load_sd(&C[(i*88)+24]);
__m128d a23_3 = _mm_load_sd(&values[299]);
#if defined(__SSE3__) && defined(__AVX__)
c23_3 = _mm_add_sd(c23_3, _mm_mul_sd(a23_3, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c23_3 = _mm_add_sd(c23_3, _mm_mul_sd(a23_3, b23));
#endif
_mm_store_sd(&C[(i*88)+24], c23_3);
__m128d c23_4 = _mm_load_sd(&C[(i*88)+41]);
__m128d a23_4 = _mm_load_sd(&values[300]);
#if defined(__SSE3__) && defined(__AVX__)
c23_4 = _mm_add_sd(c23_4, _mm_mul_sd(a23_4, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c23_4 = _mm_add_sd(c23_4, _mm_mul_sd(a23_4, b23));
#endif
_mm_store_sd(&C[(i*88)+41], c23_4);
__m128d c23_5 = _mm_loadu_pd(&C[(i*88)+43]);
__m128d a23_5 = _mm_loadu_pd(&values[301]);
#if defined(__SSE3__) && defined(__AVX__)
c23_5 = _mm_add_pd(c23_5, _mm_mul_pd(a23_5, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c23_5 = _mm_add_pd(c23_5, _mm_mul_pd(a23_5, b23));
#endif
_mm_storeu_pd(&C[(i*88)+43], c23_5);
__m128d c23_7 = _mm_load_sd(&C[(i*88)+45]);
__m128d a23_7 = _mm_load_sd(&values[303]);
#if defined(__SSE3__) && defined(__AVX__)
c23_7 = _mm_add_sd(c23_7, _mm_mul_sd(a23_7, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c23_7 = _mm_add_sd(c23_7, _mm_mul_sd(a23_7, b23));
#endif
_mm_store_sd(&C[(i*88)+45], c23_7);
__m128d c23_8 = _mm_load_sd(&C[(i*88)+69]);
__m128d a23_8 = _mm_load_sd(&values[304]);
#if defined(__SSE3__) && defined(__AVX__)
c23_8 = _mm_add_sd(c23_8, _mm_mul_sd(a23_8, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c23_8 = _mm_add_sd(c23_8, _mm_mul_sd(a23_8, b23));
#endif
_mm_store_sd(&C[(i*88)+69], c23_8);
__m128d c23_9 = _mm_loadu_pd(&C[(i*88)+71]);
__m128d a23_9 = _mm_loadu_pd(&values[305]);
#if defined(__SSE3__) && defined(__AVX__)
c23_9 = _mm_add_pd(c23_9, _mm_mul_pd(a23_9, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c23_9 = _mm_add_pd(c23_9, _mm_mul_pd(a23_9, b23));
#endif
_mm_storeu_pd(&C[(i*88)+71], c23_9);
__m128d c23_11 = _mm_load_sd(&C[(i*88)+73]);
__m128d a23_11 = _mm_load_sd(&values[307]);
#if defined(__SSE3__) && defined(__AVX__)
c23_11 = _mm_add_sd(c23_11, _mm_mul_sd(a23_11, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c23_11 = _mm_add_sd(c23_11, _mm_mul_sd(a23_11, b23));
#endif
_mm_store_sd(&C[(i*88)+73], c23_11);
#else
C[(i*88)+20] += values[296] * B[(i*88)+23];
C[(i*88)+22] += values[297] * B[(i*88)+23];
C[(i*88)+23] += values[298] * B[(i*88)+23];
C[(i*88)+24] += values[299] * B[(i*88)+23];
C[(i*88)+41] += values[300] * B[(i*88)+23];
C[(i*88)+43] += values[301] * B[(i*88)+23];
C[(i*88)+44] += values[302] * B[(i*88)+23];
C[(i*88)+45] += values[303] * B[(i*88)+23];
C[(i*88)+69] += values[304] * B[(i*88)+23];
C[(i*88)+71] += values[305] * B[(i*88)+23];
C[(i*88)+72] += values[306] * B[(i*88)+23];
C[(i*88)+73] += values[307] * B[(i*88)+23];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b24 = _mm256_broadcast_sd(&B[(i*88)+24]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b24 = _mm_loaddup_pd(&B[(i*88)+24]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c24_0 = _mm256_loadu_pd(&C[(i*88)+20]);
__m256d a24_0 = _mm256_loadu_pd(&values[308]);
c24_0 = _mm256_add_pd(c24_0, _mm256_mul_pd(a24_0, b24));
_mm256_storeu_pd(&C[(i*88)+20], c24_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c24_0 = _mm_loadu_pd(&C[(i*88)+20]);
__m128d a24_0 = _mm_loadu_pd(&values[308]);
c24_0 = _mm_add_pd(c24_0, _mm_mul_pd(a24_0, b24));
_mm_storeu_pd(&C[(i*88)+20], c24_0);
__m128d c24_2 = _mm_loadu_pd(&C[(i*88)+22]);
__m128d a24_2 = _mm_loadu_pd(&values[310]);
c24_2 = _mm_add_pd(c24_2, _mm_mul_pd(a24_2, b24));
_mm_storeu_pd(&C[(i*88)+22], c24_2);
#endif
__m128d c24_4 = _mm_load_sd(&C[(i*88)+24]);
__m128d a24_4 = _mm_load_sd(&values[312]);
#if defined(__SSE3__) && defined(__AVX__)
c24_4 = _mm_add_sd(c24_4, _mm_mul_sd(a24_4, _mm256_castpd256_pd128(b24)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c24_4 = _mm_add_sd(c24_4, _mm_mul_sd(a24_4, b24));
#endif
_mm_store_sd(&C[(i*88)+24], c24_4);
#if defined(__SSE3__) && defined(__AVX__)
__m256d c24_5 = _mm256_loadu_pd(&C[(i*88)+41]);
__m256d a24_5 = _mm256_loadu_pd(&values[313]);
c24_5 = _mm256_add_pd(c24_5, _mm256_mul_pd(a24_5, b24));
_mm256_storeu_pd(&C[(i*88)+41], c24_5);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c24_5 = _mm_loadu_pd(&C[(i*88)+41]);
__m128d a24_5 = _mm_loadu_pd(&values[313]);
c24_5 = _mm_add_pd(c24_5, _mm_mul_pd(a24_5, b24));
_mm_storeu_pd(&C[(i*88)+41], c24_5);
__m128d c24_7 = _mm_loadu_pd(&C[(i*88)+43]);
__m128d a24_7 = _mm_loadu_pd(&values[315]);
c24_7 = _mm_add_pd(c24_7, _mm_mul_pd(a24_7, b24));
_mm_storeu_pd(&C[(i*88)+43], c24_7);
#endif
__m128d c24_9 = _mm_load_sd(&C[(i*88)+45]);
__m128d a24_9 = _mm_load_sd(&values[317]);
#if defined(__SSE3__) && defined(__AVX__)
c24_9 = _mm_add_sd(c24_9, _mm_mul_sd(a24_9, _mm256_castpd256_pd128(b24)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c24_9 = _mm_add_sd(c24_9, _mm_mul_sd(a24_9, b24));
#endif
_mm_store_sd(&C[(i*88)+45], c24_9);
#if defined(__SSE3__) && defined(__AVX__)
__m256d c24_10 = _mm256_loadu_pd(&C[(i*88)+69]);
__m256d a24_10 = _mm256_loadu_pd(&values[318]);
c24_10 = _mm256_add_pd(c24_10, _mm256_mul_pd(a24_10, b24));
_mm256_storeu_pd(&C[(i*88)+69], c24_10);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c24_10 = _mm_loadu_pd(&C[(i*88)+69]);
__m128d a24_10 = _mm_loadu_pd(&values[318]);
c24_10 = _mm_add_pd(c24_10, _mm_mul_pd(a24_10, b24));
_mm_storeu_pd(&C[(i*88)+69], c24_10);
__m128d c24_12 = _mm_loadu_pd(&C[(i*88)+71]);
__m128d a24_12 = _mm_loadu_pd(&values[320]);
c24_12 = _mm_add_pd(c24_12, _mm_mul_pd(a24_12, b24));
_mm_storeu_pd(&C[(i*88)+71], c24_12);
#endif
__m128d c24_14 = _mm_load_sd(&C[(i*88)+73]);
__m128d a24_14 = _mm_load_sd(&values[322]);
#if defined(__SSE3__) && defined(__AVX__)
c24_14 = _mm_add_sd(c24_14, _mm_mul_sd(a24_14, _mm256_castpd256_pd128(b24)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c24_14 = _mm_add_sd(c24_14, _mm_mul_sd(a24_14, b24));
#endif
_mm_store_sd(&C[(i*88)+73], c24_14);
#else
C[(i*88)+20] += values[308] * B[(i*88)+24];
C[(i*88)+21] += values[309] * B[(i*88)+24];
C[(i*88)+22] += values[310] * B[(i*88)+24];
C[(i*88)+23] += values[311] * B[(i*88)+24];
C[(i*88)+24] += values[312] * B[(i*88)+24];
C[(i*88)+41] += values[313] * B[(i*88)+24];
C[(i*88)+42] += values[314] * B[(i*88)+24];
C[(i*88)+43] += values[315] * B[(i*88)+24];
C[(i*88)+44] += values[316] * B[(i*88)+24];
C[(i*88)+45] += values[317] * B[(i*88)+24];
C[(i*88)+69] += values[318] * B[(i*88)+24];
C[(i*88)+70] += values[319] * B[(i*88)+24];
C[(i*88)+71] += values[320] * B[(i*88)+24];
C[(i*88)+72] += values[321] * B[(i*88)+24];
C[(i*88)+73] += values[322] * B[(i*88)+24];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b25 = _mm256_broadcast_sd(&B[(i*88)+25]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b25 = _mm_loaddup_pd(&B[(i*88)+25]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c25_0 = _mm256_loadu_pd(&C[(i*88)+10]);
__m256d a25_0 = _mm256_loadu_pd(&values[323]);
c25_0 = _mm256_add_pd(c25_0, _mm256_mul_pd(a25_0, b25));
_mm256_storeu_pd(&C[(i*88)+10], c25_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c25_0 = _mm_loadu_pd(&C[(i*88)+10]);
__m128d a25_0 = _mm_loadu_pd(&values[323]);
c25_0 = _mm_add_pd(c25_0, _mm_mul_pd(a25_0, b25));
_mm_storeu_pd(&C[(i*88)+10], c25_0);
__m128d c25_2 = _mm_loadu_pd(&C[(i*88)+12]);
__m128d a25_2 = _mm_loadu_pd(&values[325]);
c25_2 = _mm_add_pd(c25_2, _mm_mul_pd(a25_2, b25));
_mm_storeu_pd(&C[(i*88)+12], c25_2);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c25_4 = _mm256_loadu_pd(&C[(i*88)+25]);
__m256d a25_4 = _mm256_loadu_pd(&values[327]);
c25_4 = _mm256_add_pd(c25_4, _mm256_mul_pd(a25_4, b25));
_mm256_storeu_pd(&C[(i*88)+25], c25_4);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c25_4 = _mm_loadu_pd(&C[(i*88)+25]);
__m128d a25_4 = _mm_loadu_pd(&values[327]);
c25_4 = _mm_add_pd(c25_4, _mm_mul_pd(a25_4, b25));
_mm_storeu_pd(&C[(i*88)+25], c25_4);
__m128d c25_6 = _mm_loadu_pd(&C[(i*88)+27]);
__m128d a25_6 = _mm_loadu_pd(&values[329]);
c25_6 = _mm_add_pd(c25_6, _mm_mul_pd(a25_6, b25));
_mm_storeu_pd(&C[(i*88)+27], c25_6);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c25_8 = _mm256_loadu_pd(&C[(i*88)+46]);
__m256d a25_8 = _mm256_loadu_pd(&values[331]);
c25_8 = _mm256_add_pd(c25_8, _mm256_mul_pd(a25_8, b25));
_mm256_storeu_pd(&C[(i*88)+46], c25_8);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c25_8 = _mm_loadu_pd(&C[(i*88)+46]);
__m128d a25_8 = _mm_loadu_pd(&values[331]);
c25_8 = _mm_add_pd(c25_8, _mm_mul_pd(a25_8, b25));
_mm_storeu_pd(&C[(i*88)+46], c25_8);
__m128d c25_10 = _mm_loadu_pd(&C[(i*88)+48]);
__m128d a25_10 = _mm_loadu_pd(&values[333]);
c25_10 = _mm_add_pd(c25_10, _mm_mul_pd(a25_10, b25));
_mm_storeu_pd(&C[(i*88)+48], c25_10);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c25_12 = _mm256_loadu_pd(&C[(i*88)+74]);
__m256d a25_12 = _mm256_loadu_pd(&values[335]);
c25_12 = _mm256_add_pd(c25_12, _mm256_mul_pd(a25_12, b25));
_mm256_storeu_pd(&C[(i*88)+74], c25_12);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c25_12 = _mm_loadu_pd(&C[(i*88)+74]);
__m128d a25_12 = _mm_loadu_pd(&values[335]);
c25_12 = _mm_add_pd(c25_12, _mm_mul_pd(a25_12, b25));
_mm_storeu_pd(&C[(i*88)+74], c25_12);
__m128d c25_14 = _mm_loadu_pd(&C[(i*88)+76]);
__m128d a25_14 = _mm_loadu_pd(&values[337]);
c25_14 = _mm_add_pd(c25_14, _mm_mul_pd(a25_14, b25));
_mm_storeu_pd(&C[(i*88)+76], c25_14);
#endif
#else
C[(i*88)+10] += values[323] * B[(i*88)+25];
C[(i*88)+11] += values[324] * B[(i*88)+25];
C[(i*88)+12] += values[325] * B[(i*88)+25];
C[(i*88)+13] += values[326] * B[(i*88)+25];
C[(i*88)+25] += values[327] * B[(i*88)+25];
C[(i*88)+26] += values[328] * B[(i*88)+25];
C[(i*88)+27] += values[329] * B[(i*88)+25];
C[(i*88)+28] += values[330] * B[(i*88)+25];
C[(i*88)+46] += values[331] * B[(i*88)+25];
C[(i*88)+47] += values[332] * B[(i*88)+25];
C[(i*88)+48] += values[333] * B[(i*88)+25];
C[(i*88)+49] += values[334] * B[(i*88)+25];
C[(i*88)+74] += values[335] * B[(i*88)+25];
C[(i*88)+75] += values[336] * B[(i*88)+25];
C[(i*88)+76] += values[337] * B[(i*88)+25];
C[(i*88)+77] += values[338] * B[(i*88)+25];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b26 = _mm256_broadcast_sd(&B[(i*88)+26]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b26 = _mm_loaddup_pd(&B[(i*88)+26]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c26_0 = _mm256_loadu_pd(&C[(i*88)+10]);
__m256d a26_0 = _mm256_loadu_pd(&values[339]);
c26_0 = _mm256_add_pd(c26_0, _mm256_mul_pd(a26_0, b26));
_mm256_storeu_pd(&C[(i*88)+10], c26_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c26_0 = _mm_loadu_pd(&C[(i*88)+10]);
__m128d a26_0 = _mm_loadu_pd(&values[339]);
c26_0 = _mm_add_pd(c26_0, _mm_mul_pd(a26_0, b26));
_mm_storeu_pd(&C[(i*88)+10], c26_0);
__m128d c26_2 = _mm_loadu_pd(&C[(i*88)+12]);
__m128d a26_2 = _mm_loadu_pd(&values[341]);
c26_2 = _mm_add_pd(c26_2, _mm_mul_pd(a26_2, b26));
_mm_storeu_pd(&C[(i*88)+12], c26_2);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c26_4 = _mm256_loadu_pd(&C[(i*88)+25]);
__m256d a26_4 = _mm256_loadu_pd(&values[343]);
c26_4 = _mm256_add_pd(c26_4, _mm256_mul_pd(a26_4, b26));
_mm256_storeu_pd(&C[(i*88)+25], c26_4);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c26_4 = _mm_loadu_pd(&C[(i*88)+25]);
__m128d a26_4 = _mm_loadu_pd(&values[343]);
c26_4 = _mm_add_pd(c26_4, _mm_mul_pd(a26_4, b26));
_mm_storeu_pd(&C[(i*88)+25], c26_4);
__m128d c26_6 = _mm_loadu_pd(&C[(i*88)+27]);
__m128d a26_6 = _mm_loadu_pd(&values[345]);
c26_6 = _mm_add_pd(c26_6, _mm_mul_pd(a26_6, b26));
_mm_storeu_pd(&C[(i*88)+27], c26_6);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c26_8 = _mm256_loadu_pd(&C[(i*88)+46]);
__m256d a26_8 = _mm256_loadu_pd(&values[347]);
c26_8 = _mm256_add_pd(c26_8, _mm256_mul_pd(a26_8, b26));
_mm256_storeu_pd(&C[(i*88)+46], c26_8);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c26_8 = _mm_loadu_pd(&C[(i*88)+46]);
__m128d a26_8 = _mm_loadu_pd(&values[347]);
c26_8 = _mm_add_pd(c26_8, _mm_mul_pd(a26_8, b26));
_mm_storeu_pd(&C[(i*88)+46], c26_8);
__m128d c26_10 = _mm_loadu_pd(&C[(i*88)+48]);
__m128d a26_10 = _mm_loadu_pd(&values[349]);
c26_10 = _mm_add_pd(c26_10, _mm_mul_pd(a26_10, b26));
_mm_storeu_pd(&C[(i*88)+48], c26_10);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c26_12 = _mm256_loadu_pd(&C[(i*88)+74]);
__m256d a26_12 = _mm256_loadu_pd(&values[351]);
c26_12 = _mm256_add_pd(c26_12, _mm256_mul_pd(a26_12, b26));
_mm256_storeu_pd(&C[(i*88)+74], c26_12);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c26_12 = _mm_loadu_pd(&C[(i*88)+74]);
__m128d a26_12 = _mm_loadu_pd(&values[351]);
c26_12 = _mm_add_pd(c26_12, _mm_mul_pd(a26_12, b26));
_mm_storeu_pd(&C[(i*88)+74], c26_12);
__m128d c26_14 = _mm_loadu_pd(&C[(i*88)+76]);
__m128d a26_14 = _mm_loadu_pd(&values[353]);
c26_14 = _mm_add_pd(c26_14, _mm_mul_pd(a26_14, b26));
_mm_storeu_pd(&C[(i*88)+76], c26_14);
#endif
#else
C[(i*88)+10] += values[339] * B[(i*88)+26];
C[(i*88)+11] += values[340] * B[(i*88)+26];
C[(i*88)+12] += values[341] * B[(i*88)+26];
C[(i*88)+13] += values[342] * B[(i*88)+26];
C[(i*88)+25] += values[343] * B[(i*88)+26];
C[(i*88)+26] += values[344] * B[(i*88)+26];
C[(i*88)+27] += values[345] * B[(i*88)+26];
C[(i*88)+28] += values[346] * B[(i*88)+26];
C[(i*88)+46] += values[347] * B[(i*88)+26];
C[(i*88)+47] += values[348] * B[(i*88)+26];
C[(i*88)+48] += values[349] * B[(i*88)+26];
C[(i*88)+49] += values[350] * B[(i*88)+26];
C[(i*88)+74] += values[351] * B[(i*88)+26];
C[(i*88)+75] += values[352] * B[(i*88)+26];
C[(i*88)+76] += values[353] * B[(i*88)+26];
C[(i*88)+77] += values[354] * B[(i*88)+26];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b27 = _mm256_broadcast_sd(&B[(i*88)+27]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b27 = _mm_loaddup_pd(&B[(i*88)+27]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c27_0 = _mm256_loadu_pd(&C[(i*88)+10]);
__m256d a27_0 = _mm256_loadu_pd(&values[355]);
c27_0 = _mm256_add_pd(c27_0, _mm256_mul_pd(a27_0, b27));
_mm256_storeu_pd(&C[(i*88)+10], c27_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c27_0 = _mm_loadu_pd(&C[(i*88)+10]);
__m128d a27_0 = _mm_loadu_pd(&values[355]);
c27_0 = _mm_add_pd(c27_0, _mm_mul_pd(a27_0, b27));
_mm_storeu_pd(&C[(i*88)+10], c27_0);
__m128d c27_2 = _mm_loadu_pd(&C[(i*88)+12]);
__m128d a27_2 = _mm_loadu_pd(&values[357]);
c27_2 = _mm_add_pd(c27_2, _mm_mul_pd(a27_2, b27));
_mm_storeu_pd(&C[(i*88)+12], c27_2);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c27_4 = _mm256_loadu_pd(&C[(i*88)+25]);
__m256d a27_4 = _mm256_loadu_pd(&values[359]);
c27_4 = _mm256_add_pd(c27_4, _mm256_mul_pd(a27_4, b27));
_mm256_storeu_pd(&C[(i*88)+25], c27_4);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c27_4 = _mm_loadu_pd(&C[(i*88)+25]);
__m128d a27_4 = _mm_loadu_pd(&values[359]);
c27_4 = _mm_add_pd(c27_4, _mm_mul_pd(a27_4, b27));
_mm_storeu_pd(&C[(i*88)+25], c27_4);
__m128d c27_6 = _mm_loadu_pd(&C[(i*88)+27]);
__m128d a27_6 = _mm_loadu_pd(&values[361]);
c27_6 = _mm_add_pd(c27_6, _mm_mul_pd(a27_6, b27));
_mm_storeu_pd(&C[(i*88)+27], c27_6);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c27_8 = _mm256_loadu_pd(&C[(i*88)+46]);
__m256d a27_8 = _mm256_loadu_pd(&values[363]);
c27_8 = _mm256_add_pd(c27_8, _mm256_mul_pd(a27_8, b27));
_mm256_storeu_pd(&C[(i*88)+46], c27_8);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c27_8 = _mm_loadu_pd(&C[(i*88)+46]);
__m128d a27_8 = _mm_loadu_pd(&values[363]);
c27_8 = _mm_add_pd(c27_8, _mm_mul_pd(a27_8, b27));
_mm_storeu_pd(&C[(i*88)+46], c27_8);
__m128d c27_10 = _mm_loadu_pd(&C[(i*88)+48]);
__m128d a27_10 = _mm_loadu_pd(&values[365]);
c27_10 = _mm_add_pd(c27_10, _mm_mul_pd(a27_10, b27));
_mm_storeu_pd(&C[(i*88)+48], c27_10);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c27_12 = _mm256_loadu_pd(&C[(i*88)+74]);
__m256d a27_12 = _mm256_loadu_pd(&values[367]);
c27_12 = _mm256_add_pd(c27_12, _mm256_mul_pd(a27_12, b27));
_mm256_storeu_pd(&C[(i*88)+74], c27_12);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c27_12 = _mm_loadu_pd(&C[(i*88)+74]);
__m128d a27_12 = _mm_loadu_pd(&values[367]);
c27_12 = _mm_add_pd(c27_12, _mm_mul_pd(a27_12, b27));
_mm_storeu_pd(&C[(i*88)+74], c27_12);
__m128d c27_14 = _mm_loadu_pd(&C[(i*88)+76]);
__m128d a27_14 = _mm_loadu_pd(&values[369]);
c27_14 = _mm_add_pd(c27_14, _mm_mul_pd(a27_14, b27));
_mm_storeu_pd(&C[(i*88)+76], c27_14);
#endif
#else
C[(i*88)+10] += values[355] * B[(i*88)+27];
C[(i*88)+11] += values[356] * B[(i*88)+27];
C[(i*88)+12] += values[357] * B[(i*88)+27];
C[(i*88)+13] += values[358] * B[(i*88)+27];
C[(i*88)+25] += values[359] * B[(i*88)+27];
C[(i*88)+26] += values[360] * B[(i*88)+27];
C[(i*88)+27] += values[361] * B[(i*88)+27];
C[(i*88)+28] += values[362] * B[(i*88)+27];
C[(i*88)+46] += values[363] * B[(i*88)+27];
C[(i*88)+47] += values[364] * B[(i*88)+27];
C[(i*88)+48] += values[365] * B[(i*88)+27];
C[(i*88)+49] += values[366] * B[(i*88)+27];
C[(i*88)+74] += values[367] * B[(i*88)+27];
C[(i*88)+75] += values[368] * B[(i*88)+27];
C[(i*88)+76] += values[369] * B[(i*88)+27];
C[(i*88)+77] += values[370] * B[(i*88)+27];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b28 = _mm256_broadcast_sd(&B[(i*88)+28]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b28 = _mm_loaddup_pd(&B[(i*88)+28]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c28_0 = _mm256_loadu_pd(&C[(i*88)+10]);
__m256d a28_0 = _mm256_loadu_pd(&values[371]);
c28_0 = _mm256_add_pd(c28_0, _mm256_mul_pd(a28_0, b28));
_mm256_storeu_pd(&C[(i*88)+10], c28_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c28_0 = _mm_loadu_pd(&C[(i*88)+10]);
__m128d a28_0 = _mm_loadu_pd(&values[371]);
c28_0 = _mm_add_pd(c28_0, _mm_mul_pd(a28_0, b28));
_mm_storeu_pd(&C[(i*88)+10], c28_0);
__m128d c28_2 = _mm_loadu_pd(&C[(i*88)+12]);
__m128d a28_2 = _mm_loadu_pd(&values[373]);
c28_2 = _mm_add_pd(c28_2, _mm_mul_pd(a28_2, b28));
_mm_storeu_pd(&C[(i*88)+12], c28_2);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c28_4 = _mm256_loadu_pd(&C[(i*88)+25]);
__m256d a28_4 = _mm256_loadu_pd(&values[375]);
c28_4 = _mm256_add_pd(c28_4, _mm256_mul_pd(a28_4, b28));
_mm256_storeu_pd(&C[(i*88)+25], c28_4);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c28_4 = _mm_loadu_pd(&C[(i*88)+25]);
__m128d a28_4 = _mm_loadu_pd(&values[375]);
c28_4 = _mm_add_pd(c28_4, _mm_mul_pd(a28_4, b28));
_mm_storeu_pd(&C[(i*88)+25], c28_4);
__m128d c28_6 = _mm_loadu_pd(&C[(i*88)+27]);
__m128d a28_6 = _mm_loadu_pd(&values[377]);
c28_6 = _mm_add_pd(c28_6, _mm_mul_pd(a28_6, b28));
_mm_storeu_pd(&C[(i*88)+27], c28_6);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c28_8 = _mm256_loadu_pd(&C[(i*88)+46]);
__m256d a28_8 = _mm256_loadu_pd(&values[379]);
c28_8 = _mm256_add_pd(c28_8, _mm256_mul_pd(a28_8, b28));
_mm256_storeu_pd(&C[(i*88)+46], c28_8);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c28_8 = _mm_loadu_pd(&C[(i*88)+46]);
__m128d a28_8 = _mm_loadu_pd(&values[379]);
c28_8 = _mm_add_pd(c28_8, _mm_mul_pd(a28_8, b28));
_mm_storeu_pd(&C[(i*88)+46], c28_8);
__m128d c28_10 = _mm_loadu_pd(&C[(i*88)+48]);
__m128d a28_10 = _mm_loadu_pd(&values[381]);
c28_10 = _mm_add_pd(c28_10, _mm_mul_pd(a28_10, b28));
_mm_storeu_pd(&C[(i*88)+48], c28_10);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c28_12 = _mm256_loadu_pd(&C[(i*88)+74]);
__m256d a28_12 = _mm256_loadu_pd(&values[383]);
c28_12 = _mm256_add_pd(c28_12, _mm256_mul_pd(a28_12, b28));
_mm256_storeu_pd(&C[(i*88)+74], c28_12);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c28_12 = _mm_loadu_pd(&C[(i*88)+74]);
__m128d a28_12 = _mm_loadu_pd(&values[383]);
c28_12 = _mm_add_pd(c28_12, _mm_mul_pd(a28_12, b28));
_mm_storeu_pd(&C[(i*88)+74], c28_12);
__m128d c28_14 = _mm_loadu_pd(&C[(i*88)+76]);
__m128d a28_14 = _mm_loadu_pd(&values[385]);
c28_14 = _mm_add_pd(c28_14, _mm_mul_pd(a28_14, b28));
_mm_storeu_pd(&C[(i*88)+76], c28_14);
#endif
#else
C[(i*88)+10] += values[371] * B[(i*88)+28];
C[(i*88)+11] += values[372] * B[(i*88)+28];
C[(i*88)+12] += values[373] * B[(i*88)+28];
C[(i*88)+13] += values[374] * B[(i*88)+28];
C[(i*88)+25] += values[375] * B[(i*88)+28];
C[(i*88)+26] += values[376] * B[(i*88)+28];
C[(i*88)+27] += values[377] * B[(i*88)+28];
C[(i*88)+28] += values[378] * B[(i*88)+28];
C[(i*88)+46] += values[379] * B[(i*88)+28];
C[(i*88)+47] += values[380] * B[(i*88)+28];
C[(i*88)+48] += values[381] * B[(i*88)+28];
C[(i*88)+49] += values[382] * B[(i*88)+28];
C[(i*88)+74] += values[383] * B[(i*88)+28];
C[(i*88)+75] += values[384] * B[(i*88)+28];
C[(i*88)+76] += values[385] * B[(i*88)+28];
C[(i*88)+77] += values[386] * B[(i*88)+28];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b29 = _mm256_broadcast_sd(&B[(i*88)+29]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b29 = _mm_loaddup_pd(&B[(i*88)+29]);
#endif
__m128d c29_0 = _mm_loadu_pd(&C[(i*88)+4]);
__m128d a29_0 = _mm_loadu_pd(&values[387]);
#if defined(__SSE3__) && defined(__AVX__)
c29_0 = _mm_add_pd(c29_0, _mm_mul_pd(a29_0, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_0 = _mm_add_pd(c29_0, _mm_mul_pd(a29_0, b29));
#endif
_mm_storeu_pd(&C[(i*88)+4], c29_0);
__m128d c29_2 = _mm_load_sd(&C[(i*88)+6]);
__m128d a29_2 = _mm_load_sd(&values[389]);
#if defined(__SSE3__) && defined(__AVX__)
c29_2 = _mm_add_sd(c29_2, _mm_mul_sd(a29_2, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_2 = _mm_add_sd(c29_2, _mm_mul_sd(a29_2, b29));
#endif
_mm_store_sd(&C[(i*88)+6], c29_2);
__m128d c29_3 = _mm_loadu_pd(&C[(i*88)+14]);
__m128d a29_3 = _mm_loadu_pd(&values[390]);
#if defined(__SSE3__) && defined(__AVX__)
c29_3 = _mm_add_pd(c29_3, _mm_mul_pd(a29_3, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_3 = _mm_add_pd(c29_3, _mm_mul_pd(a29_3, b29));
#endif
_mm_storeu_pd(&C[(i*88)+14], c29_3);
__m128d c29_5 = _mm_load_sd(&C[(i*88)+16]);
__m128d a29_5 = _mm_load_sd(&values[392]);
#if defined(__SSE3__) && defined(__AVX__)
c29_5 = _mm_add_sd(c29_5, _mm_mul_sd(a29_5, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_5 = _mm_add_sd(c29_5, _mm_mul_sd(a29_5, b29));
#endif
_mm_store_sd(&C[(i*88)+16], c29_5);
__m128d c29_6 = _mm_loadu_pd(&C[(i*88)+29]);
__m128d a29_6 = _mm_loadu_pd(&values[393]);
#if defined(__SSE3__) && defined(__AVX__)
c29_6 = _mm_add_pd(c29_6, _mm_mul_pd(a29_6, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_6 = _mm_add_pd(c29_6, _mm_mul_pd(a29_6, b29));
#endif
_mm_storeu_pd(&C[(i*88)+29], c29_6);
__m128d c29_8 = _mm_load_sd(&C[(i*88)+31]);
__m128d a29_8 = _mm_load_sd(&values[395]);
#if defined(__SSE3__) && defined(__AVX__)
c29_8 = _mm_add_sd(c29_8, _mm_mul_sd(a29_8, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_8 = _mm_add_sd(c29_8, _mm_mul_sd(a29_8, b29));
#endif
_mm_store_sd(&C[(i*88)+31], c29_8);
__m128d c29_9 = _mm_loadu_pd(&C[(i*88)+50]);
__m128d a29_9 = _mm_loadu_pd(&values[396]);
#if defined(__SSE3__) && defined(__AVX__)
c29_9 = _mm_add_pd(c29_9, _mm_mul_pd(a29_9, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_9 = _mm_add_pd(c29_9, _mm_mul_pd(a29_9, b29));
#endif
_mm_storeu_pd(&C[(i*88)+50], c29_9);
__m128d c29_11 = _mm_load_sd(&C[(i*88)+52]);
__m128d a29_11 = _mm_load_sd(&values[398]);
#if defined(__SSE3__) && defined(__AVX__)
c29_11 = _mm_add_sd(c29_11, _mm_mul_sd(a29_11, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_11 = _mm_add_sd(c29_11, _mm_mul_sd(a29_11, b29));
#endif
_mm_store_sd(&C[(i*88)+52], c29_11);
__m128d c29_12 = _mm_loadu_pd(&C[(i*88)+78]);
__m128d a29_12 = _mm_loadu_pd(&values[399]);
#if defined(__SSE3__) && defined(__AVX__)
c29_12 = _mm_add_pd(c29_12, _mm_mul_pd(a29_12, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_12 = _mm_add_pd(c29_12, _mm_mul_pd(a29_12, b29));
#endif
_mm_storeu_pd(&C[(i*88)+78], c29_12);
__m128d c29_14 = _mm_load_sd(&C[(i*88)+80]);
__m128d a29_14 = _mm_load_sd(&values[401]);
#if defined(__SSE3__) && defined(__AVX__)
c29_14 = _mm_add_sd(c29_14, _mm_mul_sd(a29_14, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_14 = _mm_add_sd(c29_14, _mm_mul_sd(a29_14, b29));
#endif
_mm_store_sd(&C[(i*88)+80], c29_14);
#else
C[(i*88)+4] += values[387] * B[(i*88)+29];
C[(i*88)+5] += values[388] * B[(i*88)+29];
C[(i*88)+6] += values[389] * B[(i*88)+29];
C[(i*88)+14] += values[390] * B[(i*88)+29];
C[(i*88)+15] += values[391] * B[(i*88)+29];
C[(i*88)+16] += values[392] * B[(i*88)+29];
C[(i*88)+29] += values[393] * B[(i*88)+29];
C[(i*88)+30] += values[394] * B[(i*88)+29];
C[(i*88)+31] += values[395] * B[(i*88)+29];
C[(i*88)+50] += values[396] * B[(i*88)+29];
C[(i*88)+51] += values[397] * B[(i*88)+29];
C[(i*88)+52] += values[398] * B[(i*88)+29];
C[(i*88)+78] += values[399] * B[(i*88)+29];
C[(i*88)+79] += values[400] * B[(i*88)+29];
C[(i*88)+80] += values[401] * B[(i*88)+29];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b30 = _mm256_broadcast_sd(&B[(i*88)+30]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b30 = _mm_loaddup_pd(&B[(i*88)+30]);
#endif
__m128d c30_0 = _mm_loadu_pd(&C[(i*88)+4]);
__m128d a30_0 = _mm_loadu_pd(&values[402]);
#if defined(__SSE3__) && defined(__AVX__)
c30_0 = _mm_add_pd(c30_0, _mm_mul_pd(a30_0, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_0 = _mm_add_pd(c30_0, _mm_mul_pd(a30_0, b30));
#endif
_mm_storeu_pd(&C[(i*88)+4], c30_0);
__m128d c30_2 = _mm_load_sd(&C[(i*88)+6]);
__m128d a30_2 = _mm_load_sd(&values[404]);
#if defined(__SSE3__) && defined(__AVX__)
c30_2 = _mm_add_sd(c30_2, _mm_mul_sd(a30_2, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_2 = _mm_add_sd(c30_2, _mm_mul_sd(a30_2, b30));
#endif
_mm_store_sd(&C[(i*88)+6], c30_2);
__m128d c30_3 = _mm_loadu_pd(&C[(i*88)+14]);
__m128d a30_3 = _mm_loadu_pd(&values[405]);
#if defined(__SSE3__) && defined(__AVX__)
c30_3 = _mm_add_pd(c30_3, _mm_mul_pd(a30_3, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_3 = _mm_add_pd(c30_3, _mm_mul_pd(a30_3, b30));
#endif
_mm_storeu_pd(&C[(i*88)+14], c30_3);
__m128d c30_5 = _mm_load_sd(&C[(i*88)+16]);
__m128d a30_5 = _mm_load_sd(&values[407]);
#if defined(__SSE3__) && defined(__AVX__)
c30_5 = _mm_add_sd(c30_5, _mm_mul_sd(a30_5, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_5 = _mm_add_sd(c30_5, _mm_mul_sd(a30_5, b30));
#endif
_mm_store_sd(&C[(i*88)+16], c30_5);
__m128d c30_6 = _mm_loadu_pd(&C[(i*88)+29]);
__m128d a30_6 = _mm_loadu_pd(&values[408]);
#if defined(__SSE3__) && defined(__AVX__)
c30_6 = _mm_add_pd(c30_6, _mm_mul_pd(a30_6, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_6 = _mm_add_pd(c30_6, _mm_mul_pd(a30_6, b30));
#endif
_mm_storeu_pd(&C[(i*88)+29], c30_6);
__m128d c30_8 = _mm_load_sd(&C[(i*88)+31]);
__m128d a30_8 = _mm_load_sd(&values[410]);
#if defined(__SSE3__) && defined(__AVX__)
c30_8 = _mm_add_sd(c30_8, _mm_mul_sd(a30_8, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_8 = _mm_add_sd(c30_8, _mm_mul_sd(a30_8, b30));
#endif
_mm_store_sd(&C[(i*88)+31], c30_8);
__m128d c30_9 = _mm_loadu_pd(&C[(i*88)+50]);
__m128d a30_9 = _mm_loadu_pd(&values[411]);
#if defined(__SSE3__) && defined(__AVX__)
c30_9 = _mm_add_pd(c30_9, _mm_mul_pd(a30_9, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_9 = _mm_add_pd(c30_9, _mm_mul_pd(a30_9, b30));
#endif
_mm_storeu_pd(&C[(i*88)+50], c30_9);
__m128d c30_11 = _mm_load_sd(&C[(i*88)+52]);
__m128d a30_11 = _mm_load_sd(&values[413]);
#if defined(__SSE3__) && defined(__AVX__)
c30_11 = _mm_add_sd(c30_11, _mm_mul_sd(a30_11, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_11 = _mm_add_sd(c30_11, _mm_mul_sd(a30_11, b30));
#endif
_mm_store_sd(&C[(i*88)+52], c30_11);
__m128d c30_12 = _mm_loadu_pd(&C[(i*88)+78]);
__m128d a30_12 = _mm_loadu_pd(&values[414]);
#if defined(__SSE3__) && defined(__AVX__)
c30_12 = _mm_add_pd(c30_12, _mm_mul_pd(a30_12, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_12 = _mm_add_pd(c30_12, _mm_mul_pd(a30_12, b30));
#endif
_mm_storeu_pd(&C[(i*88)+78], c30_12);
__m128d c30_14 = _mm_load_sd(&C[(i*88)+80]);
__m128d a30_14 = _mm_load_sd(&values[416]);
#if defined(__SSE3__) && defined(__AVX__)
c30_14 = _mm_add_sd(c30_14, _mm_mul_sd(a30_14, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_14 = _mm_add_sd(c30_14, _mm_mul_sd(a30_14, b30));
#endif
_mm_store_sd(&C[(i*88)+80], c30_14);
#else
C[(i*88)+4] += values[402] * B[(i*88)+30];
C[(i*88)+5] += values[403] * B[(i*88)+30];
C[(i*88)+6] += values[404] * B[(i*88)+30];
C[(i*88)+14] += values[405] * B[(i*88)+30];
C[(i*88)+15] += values[406] * B[(i*88)+30];
C[(i*88)+16] += values[407] * B[(i*88)+30];
C[(i*88)+29] += values[408] * B[(i*88)+30];
C[(i*88)+30] += values[409] * B[(i*88)+30];
C[(i*88)+31] += values[410] * B[(i*88)+30];
C[(i*88)+50] += values[411] * B[(i*88)+30];
C[(i*88)+51] += values[412] * B[(i*88)+30];
C[(i*88)+52] += values[413] * B[(i*88)+30];
C[(i*88)+78] += values[414] * B[(i*88)+30];
C[(i*88)+79] += values[415] * B[(i*88)+30];
C[(i*88)+80] += values[416] * B[(i*88)+30];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b31 = _mm256_broadcast_sd(&B[(i*88)+31]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b31 = _mm_loaddup_pd(&B[(i*88)+31]);
#endif
__m128d c31_0 = _mm_loadu_pd(&C[(i*88)+4]);
__m128d a31_0 = _mm_loadu_pd(&values[417]);
#if defined(__SSE3__) && defined(__AVX__)
c31_0 = _mm_add_pd(c31_0, _mm_mul_pd(a31_0, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_0 = _mm_add_pd(c31_0, _mm_mul_pd(a31_0, b31));
#endif
_mm_storeu_pd(&C[(i*88)+4], c31_0);
__m128d c31_2 = _mm_load_sd(&C[(i*88)+6]);
__m128d a31_2 = _mm_load_sd(&values[419]);
#if defined(__SSE3__) && defined(__AVX__)
c31_2 = _mm_add_sd(c31_2, _mm_mul_sd(a31_2, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_2 = _mm_add_sd(c31_2, _mm_mul_sd(a31_2, b31));
#endif
_mm_store_sd(&C[(i*88)+6], c31_2);
__m128d c31_3 = _mm_loadu_pd(&C[(i*88)+14]);
__m128d a31_3 = _mm_loadu_pd(&values[420]);
#if defined(__SSE3__) && defined(__AVX__)
c31_3 = _mm_add_pd(c31_3, _mm_mul_pd(a31_3, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_3 = _mm_add_pd(c31_3, _mm_mul_pd(a31_3, b31));
#endif
_mm_storeu_pd(&C[(i*88)+14], c31_3);
__m128d c31_5 = _mm_load_sd(&C[(i*88)+16]);
__m128d a31_5 = _mm_load_sd(&values[422]);
#if defined(__SSE3__) && defined(__AVX__)
c31_5 = _mm_add_sd(c31_5, _mm_mul_sd(a31_5, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_5 = _mm_add_sd(c31_5, _mm_mul_sd(a31_5, b31));
#endif
_mm_store_sd(&C[(i*88)+16], c31_5);
__m128d c31_6 = _mm_loadu_pd(&C[(i*88)+29]);
__m128d a31_6 = _mm_loadu_pd(&values[423]);
#if defined(__SSE3__) && defined(__AVX__)
c31_6 = _mm_add_pd(c31_6, _mm_mul_pd(a31_6, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_6 = _mm_add_pd(c31_6, _mm_mul_pd(a31_6, b31));
#endif
_mm_storeu_pd(&C[(i*88)+29], c31_6);
__m128d c31_8 = _mm_load_sd(&C[(i*88)+31]);
__m128d a31_8 = _mm_load_sd(&values[425]);
#if defined(__SSE3__) && defined(__AVX__)
c31_8 = _mm_add_sd(c31_8, _mm_mul_sd(a31_8, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_8 = _mm_add_sd(c31_8, _mm_mul_sd(a31_8, b31));
#endif
_mm_store_sd(&C[(i*88)+31], c31_8);
__m128d c31_9 = _mm_loadu_pd(&C[(i*88)+50]);
__m128d a31_9 = _mm_loadu_pd(&values[426]);
#if defined(__SSE3__) && defined(__AVX__)
c31_9 = _mm_add_pd(c31_9, _mm_mul_pd(a31_9, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_9 = _mm_add_pd(c31_9, _mm_mul_pd(a31_9, b31));
#endif
_mm_storeu_pd(&C[(i*88)+50], c31_9);
__m128d c31_11 = _mm_load_sd(&C[(i*88)+52]);
__m128d a31_11 = _mm_load_sd(&values[428]);
#if defined(__SSE3__) && defined(__AVX__)
c31_11 = _mm_add_sd(c31_11, _mm_mul_sd(a31_11, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_11 = _mm_add_sd(c31_11, _mm_mul_sd(a31_11, b31));
#endif
_mm_store_sd(&C[(i*88)+52], c31_11);
__m128d c31_12 = _mm_loadu_pd(&C[(i*88)+78]);
__m128d a31_12 = _mm_loadu_pd(&values[429]);
#if defined(__SSE3__) && defined(__AVX__)
c31_12 = _mm_add_pd(c31_12, _mm_mul_pd(a31_12, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_12 = _mm_add_pd(c31_12, _mm_mul_pd(a31_12, b31));
#endif
_mm_storeu_pd(&C[(i*88)+78], c31_12);
__m128d c31_14 = _mm_load_sd(&C[(i*88)+80]);
__m128d a31_14 = _mm_load_sd(&values[431]);
#if defined(__SSE3__) && defined(__AVX__)
c31_14 = _mm_add_sd(c31_14, _mm_mul_sd(a31_14, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_14 = _mm_add_sd(c31_14, _mm_mul_sd(a31_14, b31));
#endif
_mm_store_sd(&C[(i*88)+80], c31_14);
#else
C[(i*88)+4] += values[417] * B[(i*88)+31];
C[(i*88)+5] += values[418] * B[(i*88)+31];
C[(i*88)+6] += values[419] * B[(i*88)+31];
C[(i*88)+14] += values[420] * B[(i*88)+31];
C[(i*88)+15] += values[421] * B[(i*88)+31];
C[(i*88)+16] += values[422] * B[(i*88)+31];
C[(i*88)+29] += values[423] * B[(i*88)+31];
C[(i*88)+30] += values[424] * B[(i*88)+31];
C[(i*88)+31] += values[425] * B[(i*88)+31];
C[(i*88)+50] += values[426] * B[(i*88)+31];
C[(i*88)+51] += values[427] * B[(i*88)+31];
C[(i*88)+52] += values[428] * B[(i*88)+31];
C[(i*88)+78] += values[429] * B[(i*88)+31];
C[(i*88)+79] += values[430] * B[(i*88)+31];
C[(i*88)+80] += values[431] * B[(i*88)+31];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b32 = _mm256_broadcast_sd(&B[(i*88)+32]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b32 = _mm_loaddup_pd(&B[(i*88)+32]);
#endif
__m128d c32_0 = _mm_loadu_pd(&C[(i*88)+1]);
__m128d a32_0 = _mm_loadu_pd(&values[432]);
#if defined(__SSE3__) && defined(__AVX__)
c32_0 = _mm_add_pd(c32_0, _mm_mul_pd(a32_0, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_0 = _mm_add_pd(c32_0, _mm_mul_pd(a32_0, b32));
#endif
_mm_storeu_pd(&C[(i*88)+1], c32_0);
__m128d c32_2 = _mm_loadu_pd(&C[(i*88)+7]);
__m128d a32_2 = _mm_loadu_pd(&values[434]);
#if defined(__SSE3__) && defined(__AVX__)
c32_2 = _mm_add_pd(c32_2, _mm_mul_pd(a32_2, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_2 = _mm_add_pd(c32_2, _mm_mul_pd(a32_2, b32));
#endif
_mm_storeu_pd(&C[(i*88)+7], c32_2);
__m128d c32_4 = _mm_loadu_pd(&C[(i*88)+17]);
__m128d a32_4 = _mm_loadu_pd(&values[436]);
#if defined(__SSE3__) && defined(__AVX__)
c32_4 = _mm_add_pd(c32_4, _mm_mul_pd(a32_4, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_4 = _mm_add_pd(c32_4, _mm_mul_pd(a32_4, b32));
#endif
_mm_storeu_pd(&C[(i*88)+17], c32_4);
__m128d c32_6 = _mm_loadu_pd(&C[(i*88)+32]);
__m128d a32_6 = _mm_loadu_pd(&values[438]);
#if defined(__SSE3__) && defined(__AVX__)
c32_6 = _mm_add_pd(c32_6, _mm_mul_pd(a32_6, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_6 = _mm_add_pd(c32_6, _mm_mul_pd(a32_6, b32));
#endif
_mm_storeu_pd(&C[(i*88)+32], c32_6);
__m128d c32_8 = _mm_loadu_pd(&C[(i*88)+53]);
__m128d a32_8 = _mm_loadu_pd(&values[440]);
#if defined(__SSE3__) && defined(__AVX__)
c32_8 = _mm_add_pd(c32_8, _mm_mul_pd(a32_8, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_8 = _mm_add_pd(c32_8, _mm_mul_pd(a32_8, b32));
#endif
_mm_storeu_pd(&C[(i*88)+53], c32_8);
__m128d c32_10 = _mm_loadu_pd(&C[(i*88)+81]);
__m128d a32_10 = _mm_loadu_pd(&values[442]);
#if defined(__SSE3__) && defined(__AVX__)
c32_10 = _mm_add_pd(c32_10, _mm_mul_pd(a32_10, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_10 = _mm_add_pd(c32_10, _mm_mul_pd(a32_10, b32));
#endif
_mm_storeu_pd(&C[(i*88)+81], c32_10);
#else
C[(i*88)+1] += values[432] * B[(i*88)+32];
C[(i*88)+2] += values[433] * B[(i*88)+32];
C[(i*88)+7] += values[434] * B[(i*88)+32];
C[(i*88)+8] += values[435] * B[(i*88)+32];
C[(i*88)+17] += values[436] * B[(i*88)+32];
C[(i*88)+18] += values[437] * B[(i*88)+32];
C[(i*88)+32] += values[438] * B[(i*88)+32];
C[(i*88)+33] += values[439] * B[(i*88)+32];
C[(i*88)+53] += values[440] * B[(i*88)+32];
C[(i*88)+54] += values[441] * B[(i*88)+32];
C[(i*88)+81] += values[442] * B[(i*88)+32];
C[(i*88)+82] += values[443] * B[(i*88)+32];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b33 = _mm256_broadcast_sd(&B[(i*88)+33]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b33 = _mm_loaddup_pd(&B[(i*88)+33]);
#endif
__m128d c33_0 = _mm_loadu_pd(&C[(i*88)+1]);
__m128d a33_0 = _mm_loadu_pd(&values[444]);
#if defined(__SSE3__) && defined(__AVX__)
c33_0 = _mm_add_pd(c33_0, _mm_mul_pd(a33_0, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_0 = _mm_add_pd(c33_0, _mm_mul_pd(a33_0, b33));
#endif
_mm_storeu_pd(&C[(i*88)+1], c33_0);
__m128d c33_2 = _mm_loadu_pd(&C[(i*88)+7]);
__m128d a33_2 = _mm_loadu_pd(&values[446]);
#if defined(__SSE3__) && defined(__AVX__)
c33_2 = _mm_add_pd(c33_2, _mm_mul_pd(a33_2, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_2 = _mm_add_pd(c33_2, _mm_mul_pd(a33_2, b33));
#endif
_mm_storeu_pd(&C[(i*88)+7], c33_2);
__m128d c33_4 = _mm_loadu_pd(&C[(i*88)+17]);
__m128d a33_4 = _mm_loadu_pd(&values[448]);
#if defined(__SSE3__) && defined(__AVX__)
c33_4 = _mm_add_pd(c33_4, _mm_mul_pd(a33_4, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_4 = _mm_add_pd(c33_4, _mm_mul_pd(a33_4, b33));
#endif
_mm_storeu_pd(&C[(i*88)+17], c33_4);
__m128d c33_6 = _mm_loadu_pd(&C[(i*88)+32]);
__m128d a33_6 = _mm_loadu_pd(&values[450]);
#if defined(__SSE3__) && defined(__AVX__)
c33_6 = _mm_add_pd(c33_6, _mm_mul_pd(a33_6, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_6 = _mm_add_pd(c33_6, _mm_mul_pd(a33_6, b33));
#endif
_mm_storeu_pd(&C[(i*88)+32], c33_6);
__m128d c33_8 = _mm_loadu_pd(&C[(i*88)+53]);
__m128d a33_8 = _mm_loadu_pd(&values[452]);
#if defined(__SSE3__) && defined(__AVX__)
c33_8 = _mm_add_pd(c33_8, _mm_mul_pd(a33_8, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_8 = _mm_add_pd(c33_8, _mm_mul_pd(a33_8, b33));
#endif
_mm_storeu_pd(&C[(i*88)+53], c33_8);
__m128d c33_10 = _mm_loadu_pd(&C[(i*88)+81]);
__m128d a33_10 = _mm_loadu_pd(&values[454]);
#if defined(__SSE3__) && defined(__AVX__)
c33_10 = _mm_add_pd(c33_10, _mm_mul_pd(a33_10, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_10 = _mm_add_pd(c33_10, _mm_mul_pd(a33_10, b33));
#endif
_mm_storeu_pd(&C[(i*88)+81], c33_10);
#else
C[(i*88)+1] += values[444] * B[(i*88)+33];
C[(i*88)+2] += values[445] * B[(i*88)+33];
C[(i*88)+7] += values[446] * B[(i*88)+33];
C[(i*88)+8] += values[447] * B[(i*88)+33];
C[(i*88)+17] += values[448] * B[(i*88)+33];
C[(i*88)+18] += values[449] * B[(i*88)+33];
C[(i*88)+32] += values[450] * B[(i*88)+33];
C[(i*88)+33] += values[451] * B[(i*88)+33];
C[(i*88)+53] += values[452] * B[(i*88)+33];
C[(i*88)+54] += values[453] * B[(i*88)+33];
C[(i*88)+81] += values[454] * B[(i*88)+33];
C[(i*88)+82] += values[455] * B[(i*88)+33];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b34 = _mm256_broadcast_sd(&B[(i*88)+34]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b34 = _mm_loaddup_pd(&B[(i*88)+34]);
#endif
__m128d c34_0 = _mm_load_sd(&C[(i*88)+0]);
__m128d a34_0 = _mm_load_sd(&values[456]);
#if defined(__SSE3__) && defined(__AVX__)
c34_0 = _mm_add_sd(c34_0, _mm_mul_sd(a34_0, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_0 = _mm_add_sd(c34_0, _mm_mul_sd(a34_0, b34));
#endif
_mm_store_sd(&C[(i*88)+0], c34_0);
__m128d c34_1 = _mm_load_sd(&C[(i*88)+3]);
__m128d a34_1 = _mm_load_sd(&values[457]);
#if defined(__SSE3__) && defined(__AVX__)
c34_1 = _mm_add_sd(c34_1, _mm_mul_sd(a34_1, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_1 = _mm_add_sd(c34_1, _mm_mul_sd(a34_1, b34));
#endif
_mm_store_sd(&C[(i*88)+3], c34_1);
__m128d c34_2 = _mm_load_sd(&C[(i*88)+9]);
__m128d a34_2 = _mm_load_sd(&values[458]);
#if defined(__SSE3__) && defined(__AVX__)
c34_2 = _mm_add_sd(c34_2, _mm_mul_sd(a34_2, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_2 = _mm_add_sd(c34_2, _mm_mul_sd(a34_2, b34));
#endif
_mm_store_sd(&C[(i*88)+9], c34_2);
__m128d c34_3 = _mm_load_sd(&C[(i*88)+19]);
__m128d a34_3 = _mm_load_sd(&values[459]);
#if defined(__SSE3__) && defined(__AVX__)
c34_3 = _mm_add_sd(c34_3, _mm_mul_sd(a34_3, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_3 = _mm_add_sd(c34_3, _mm_mul_sd(a34_3, b34));
#endif
_mm_store_sd(&C[(i*88)+19], c34_3);
__m128d c34_4 = _mm_load_sd(&C[(i*88)+34]);
__m128d a34_4 = _mm_load_sd(&values[460]);
#if defined(__SSE3__) && defined(__AVX__)
c34_4 = _mm_add_sd(c34_4, _mm_mul_sd(a34_4, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_4 = _mm_add_sd(c34_4, _mm_mul_sd(a34_4, b34));
#endif
_mm_store_sd(&C[(i*88)+34], c34_4);
__m128d c34_5 = _mm_load_sd(&C[(i*88)+55]);
__m128d a34_5 = _mm_load_sd(&values[461]);
#if defined(__SSE3__) && defined(__AVX__)
c34_5 = _mm_add_sd(c34_5, _mm_mul_sd(a34_5, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_5 = _mm_add_sd(c34_5, _mm_mul_sd(a34_5, b34));
#endif
_mm_store_sd(&C[(i*88)+55], c34_5);
__m128d c34_6 = _mm_load_sd(&C[(i*88)+83]);
__m128d a34_6 = _mm_load_sd(&values[462]);
#if defined(__SSE3__) && defined(__AVX__)
c34_6 = _mm_add_sd(c34_6, _mm_mul_sd(a34_6, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_6 = _mm_add_sd(c34_6, _mm_mul_sd(a34_6, b34));
#endif
_mm_store_sd(&C[(i*88)+83], c34_6);
#else
C[(i*88)+0] += values[456] * B[(i*88)+34];
C[(i*88)+3] += values[457] * B[(i*88)+34];
C[(i*88)+9] += values[458] * B[(i*88)+34];
C[(i*88)+19] += values[459] * B[(i*88)+34];
C[(i*88)+34] += values[460] * B[(i*88)+34];
C[(i*88)+55] += values[461] * B[(i*88)+34];
C[(i*88)+83] += values[462] * B[(i*88)+34];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b35 = _mm256_broadcast_sd(&B[(i*88)+35]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b35 = _mm_loaddup_pd(&B[(i*88)+35]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c35_0 = _mm256_loadu_pd(&C[(i*88)+35]);
__m256d a35_0 = _mm256_loadu_pd(&values[463]);
c35_0 = _mm256_add_pd(c35_0, _mm256_mul_pd(a35_0, b35));
_mm256_storeu_pd(&C[(i*88)+35], c35_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c35_0 = _mm_loadu_pd(&C[(i*88)+35]);
__m128d a35_0 = _mm_loadu_pd(&values[463]);
c35_0 = _mm_add_pd(c35_0, _mm_mul_pd(a35_0, b35));
_mm_storeu_pd(&C[(i*88)+35], c35_0);
__m128d c35_2 = _mm_loadu_pd(&C[(i*88)+37]);
__m128d a35_2 = _mm_loadu_pd(&values[465]);
c35_2 = _mm_add_pd(c35_2, _mm_mul_pd(a35_2, b35));
_mm_storeu_pd(&C[(i*88)+37], c35_2);
#endif
__m128d c35_4 = _mm_loadu_pd(&C[(i*88)+39]);
__m128d a35_4 = _mm_loadu_pd(&values[467]);
#if defined(__SSE3__) && defined(__AVX__)
c35_4 = _mm_add_pd(c35_4, _mm_mul_pd(a35_4, _mm256_castpd256_pd128(b35)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c35_4 = _mm_add_pd(c35_4, _mm_mul_pd(a35_4, b35));
#endif
_mm_storeu_pd(&C[(i*88)+39], c35_4);
#if defined(__SSE3__) && defined(__AVX__)
__m256d c35_6 = _mm256_loadu_pd(&C[(i*88)+63]);
__m256d a35_6 = _mm256_loadu_pd(&values[469]);
c35_6 = _mm256_add_pd(c35_6, _mm256_mul_pd(a35_6, b35));
_mm256_storeu_pd(&C[(i*88)+63], c35_6);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c35_6 = _mm_loadu_pd(&C[(i*88)+63]);
__m128d a35_6 = _mm_loadu_pd(&values[469]);
c35_6 = _mm_add_pd(c35_6, _mm_mul_pd(a35_6, b35));
_mm_storeu_pd(&C[(i*88)+63], c35_6);
__m128d c35_8 = _mm_loadu_pd(&C[(i*88)+65]);
__m128d a35_8 = _mm_loadu_pd(&values[471]);
c35_8 = _mm_add_pd(c35_8, _mm_mul_pd(a35_8, b35));
_mm_storeu_pd(&C[(i*88)+65], c35_8);
#endif
__m128d c35_10 = _mm_loadu_pd(&C[(i*88)+67]);
__m128d a35_10 = _mm_loadu_pd(&values[473]);
#if defined(__SSE3__) && defined(__AVX__)
c35_10 = _mm_add_pd(c35_10, _mm_mul_pd(a35_10, _mm256_castpd256_pd128(b35)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c35_10 = _mm_add_pd(c35_10, _mm_mul_pd(a35_10, b35));
#endif
_mm_storeu_pd(&C[(i*88)+67], c35_10);
#else
C[(i*88)+35] += values[463] * B[(i*88)+35];
C[(i*88)+36] += values[464] * B[(i*88)+35];
C[(i*88)+37] += values[465] * B[(i*88)+35];
C[(i*88)+38] += values[466] * B[(i*88)+35];
C[(i*88)+39] += values[467] * B[(i*88)+35];
C[(i*88)+40] += values[468] * B[(i*88)+35];
C[(i*88)+63] += values[469] * B[(i*88)+35];
C[(i*88)+64] += values[470] * B[(i*88)+35];
C[(i*88)+65] += values[471] * B[(i*88)+35];
C[(i*88)+66] += values[472] * B[(i*88)+35];
C[(i*88)+67] += values[473] * B[(i*88)+35];
C[(i*88)+68] += values[474] * B[(i*88)+35];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b36 = _mm256_broadcast_sd(&B[(i*88)+36]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b36 = _mm_loaddup_pd(&B[(i*88)+36]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c36_0 = _mm256_loadu_pd(&C[(i*88)+35]);
__m256d a36_0 = _mm256_loadu_pd(&values[475]);
c36_0 = _mm256_add_pd(c36_0, _mm256_mul_pd(a36_0, b36));
_mm256_storeu_pd(&C[(i*88)+35], c36_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c36_0 = _mm_loadu_pd(&C[(i*88)+35]);
__m128d a36_0 = _mm_loadu_pd(&values[475]);
c36_0 = _mm_add_pd(c36_0, _mm_mul_pd(a36_0, b36));
_mm_storeu_pd(&C[(i*88)+35], c36_0);
__m128d c36_2 = _mm_loadu_pd(&C[(i*88)+37]);
__m128d a36_2 = _mm_loadu_pd(&values[477]);
c36_2 = _mm_add_pd(c36_2, _mm_mul_pd(a36_2, b36));
_mm_storeu_pd(&C[(i*88)+37], c36_2);
#endif
__m128d c36_4 = _mm_loadu_pd(&C[(i*88)+39]);
__m128d a36_4 = _mm_loadu_pd(&values[479]);
#if defined(__SSE3__) && defined(__AVX__)
c36_4 = _mm_add_pd(c36_4, _mm_mul_pd(a36_4, _mm256_castpd256_pd128(b36)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c36_4 = _mm_add_pd(c36_4, _mm_mul_pd(a36_4, b36));
#endif
_mm_storeu_pd(&C[(i*88)+39], c36_4);
#if defined(__SSE3__) && defined(__AVX__)
__m256d c36_6 = _mm256_loadu_pd(&C[(i*88)+63]);
__m256d a36_6 = _mm256_loadu_pd(&values[481]);
c36_6 = _mm256_add_pd(c36_6, _mm256_mul_pd(a36_6, b36));
_mm256_storeu_pd(&C[(i*88)+63], c36_6);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c36_6 = _mm_loadu_pd(&C[(i*88)+63]);
__m128d a36_6 = _mm_loadu_pd(&values[481]);
c36_6 = _mm_add_pd(c36_6, _mm_mul_pd(a36_6, b36));
_mm_storeu_pd(&C[(i*88)+63], c36_6);
__m128d c36_8 = _mm_loadu_pd(&C[(i*88)+65]);
__m128d a36_8 = _mm_loadu_pd(&values[483]);
c36_8 = _mm_add_pd(c36_8, _mm_mul_pd(a36_8, b36));
_mm_storeu_pd(&C[(i*88)+65], c36_8);
#endif
__m128d c36_10 = _mm_loadu_pd(&C[(i*88)+67]);
__m128d a36_10 = _mm_loadu_pd(&values[485]);
#if defined(__SSE3__) && defined(__AVX__)
c36_10 = _mm_add_pd(c36_10, _mm_mul_pd(a36_10, _mm256_castpd256_pd128(b36)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c36_10 = _mm_add_pd(c36_10, _mm_mul_pd(a36_10, b36));
#endif
_mm_storeu_pd(&C[(i*88)+67], c36_10);
#else
C[(i*88)+35] += values[475] * B[(i*88)+36];
C[(i*88)+36] += values[476] * B[(i*88)+36];
C[(i*88)+37] += values[477] * B[(i*88)+36];
C[(i*88)+38] += values[478] * B[(i*88)+36];
C[(i*88)+39] += values[479] * B[(i*88)+36];
C[(i*88)+40] += values[480] * B[(i*88)+36];
C[(i*88)+63] += values[481] * B[(i*88)+36];
C[(i*88)+64] += values[482] * B[(i*88)+36];
C[(i*88)+65] += values[483] * B[(i*88)+36];
C[(i*88)+66] += values[484] * B[(i*88)+36];
C[(i*88)+67] += values[485] * B[(i*88)+36];
C[(i*88)+68] += values[486] * B[(i*88)+36];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b37 = _mm256_broadcast_sd(&B[(i*88)+37]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b37 = _mm_loaddup_pd(&B[(i*88)+37]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c37_0 = _mm256_loadu_pd(&C[(i*88)+35]);
__m256d a37_0 = _mm256_loadu_pd(&values[487]);
c37_0 = _mm256_add_pd(c37_0, _mm256_mul_pd(a37_0, b37));
_mm256_storeu_pd(&C[(i*88)+35], c37_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c37_0 = _mm_loadu_pd(&C[(i*88)+35]);
__m128d a37_0 = _mm_loadu_pd(&values[487]);
c37_0 = _mm_add_pd(c37_0, _mm_mul_pd(a37_0, b37));
_mm_storeu_pd(&C[(i*88)+35], c37_0);
__m128d c37_2 = _mm_loadu_pd(&C[(i*88)+37]);
__m128d a37_2 = _mm_loadu_pd(&values[489]);
c37_2 = _mm_add_pd(c37_2, _mm_mul_pd(a37_2, b37));
_mm_storeu_pd(&C[(i*88)+37], c37_2);
#endif
__m128d c37_4 = _mm_loadu_pd(&C[(i*88)+39]);
__m128d a37_4 = _mm_loadu_pd(&values[491]);
#if defined(__SSE3__) && defined(__AVX__)
c37_4 = _mm_add_pd(c37_4, _mm_mul_pd(a37_4, _mm256_castpd256_pd128(b37)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c37_4 = _mm_add_pd(c37_4, _mm_mul_pd(a37_4, b37));
#endif
_mm_storeu_pd(&C[(i*88)+39], c37_4);
#if defined(__SSE3__) && defined(__AVX__)
__m256d c37_6 = _mm256_loadu_pd(&C[(i*88)+63]);
__m256d a37_6 = _mm256_loadu_pd(&values[493]);
c37_6 = _mm256_add_pd(c37_6, _mm256_mul_pd(a37_6, b37));
_mm256_storeu_pd(&C[(i*88)+63], c37_6);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c37_6 = _mm_loadu_pd(&C[(i*88)+63]);
__m128d a37_6 = _mm_loadu_pd(&values[493]);
c37_6 = _mm_add_pd(c37_6, _mm_mul_pd(a37_6, b37));
_mm_storeu_pd(&C[(i*88)+63], c37_6);
__m128d c37_8 = _mm_loadu_pd(&C[(i*88)+65]);
__m128d a37_8 = _mm_loadu_pd(&values[495]);
c37_8 = _mm_add_pd(c37_8, _mm_mul_pd(a37_8, b37));
_mm_storeu_pd(&C[(i*88)+65], c37_8);
#endif
__m128d c37_10 = _mm_loadu_pd(&C[(i*88)+67]);
__m128d a37_10 = _mm_loadu_pd(&values[497]);
#if defined(__SSE3__) && defined(__AVX__)
c37_10 = _mm_add_pd(c37_10, _mm_mul_pd(a37_10, _mm256_castpd256_pd128(b37)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c37_10 = _mm_add_pd(c37_10, _mm_mul_pd(a37_10, b37));
#endif
_mm_storeu_pd(&C[(i*88)+67], c37_10);
#else
C[(i*88)+35] += values[487] * B[(i*88)+37];
C[(i*88)+36] += values[488] * B[(i*88)+37];
C[(i*88)+37] += values[489] * B[(i*88)+37];
C[(i*88)+38] += values[490] * B[(i*88)+37];
C[(i*88)+39] += values[491] * B[(i*88)+37];
C[(i*88)+40] += values[492] * B[(i*88)+37];
C[(i*88)+63] += values[493] * B[(i*88)+37];
C[(i*88)+64] += values[494] * B[(i*88)+37];
C[(i*88)+65] += values[495] * B[(i*88)+37];
C[(i*88)+66] += values[496] * B[(i*88)+37];
C[(i*88)+67] += values[497] * B[(i*88)+37];
C[(i*88)+68] += values[498] * B[(i*88)+37];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b38 = _mm256_broadcast_sd(&B[(i*88)+38]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b38 = _mm_loaddup_pd(&B[(i*88)+38]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c38_0 = _mm256_loadu_pd(&C[(i*88)+35]);
__m256d a38_0 = _mm256_loadu_pd(&values[499]);
c38_0 = _mm256_add_pd(c38_0, _mm256_mul_pd(a38_0, b38));
_mm256_storeu_pd(&C[(i*88)+35], c38_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c38_0 = _mm_loadu_pd(&C[(i*88)+35]);
__m128d a38_0 = _mm_loadu_pd(&values[499]);
c38_0 = _mm_add_pd(c38_0, _mm_mul_pd(a38_0, b38));
_mm_storeu_pd(&C[(i*88)+35], c38_0);
__m128d c38_2 = _mm_loadu_pd(&C[(i*88)+37]);
__m128d a38_2 = _mm_loadu_pd(&values[501]);
c38_2 = _mm_add_pd(c38_2, _mm_mul_pd(a38_2, b38));
_mm_storeu_pd(&C[(i*88)+37], c38_2);
#endif
__m128d c38_4 = _mm_loadu_pd(&C[(i*88)+39]);
__m128d a38_4 = _mm_loadu_pd(&values[503]);
#if defined(__SSE3__) && defined(__AVX__)
c38_4 = _mm_add_pd(c38_4, _mm_mul_pd(a38_4, _mm256_castpd256_pd128(b38)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c38_4 = _mm_add_pd(c38_4, _mm_mul_pd(a38_4, b38));
#endif
_mm_storeu_pd(&C[(i*88)+39], c38_4);
#if defined(__SSE3__) && defined(__AVX__)
__m256d c38_6 = _mm256_loadu_pd(&C[(i*88)+63]);
__m256d a38_6 = _mm256_loadu_pd(&values[505]);
c38_6 = _mm256_add_pd(c38_6, _mm256_mul_pd(a38_6, b38));
_mm256_storeu_pd(&C[(i*88)+63], c38_6);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c38_6 = _mm_loadu_pd(&C[(i*88)+63]);
__m128d a38_6 = _mm_loadu_pd(&values[505]);
c38_6 = _mm_add_pd(c38_6, _mm_mul_pd(a38_6, b38));
_mm_storeu_pd(&C[(i*88)+63], c38_6);
__m128d c38_8 = _mm_loadu_pd(&C[(i*88)+65]);
__m128d a38_8 = _mm_loadu_pd(&values[507]);
c38_8 = _mm_add_pd(c38_8, _mm_mul_pd(a38_8, b38));
_mm_storeu_pd(&C[(i*88)+65], c38_8);
#endif
__m128d c38_10 = _mm_loadu_pd(&C[(i*88)+67]);
__m128d a38_10 = _mm_loadu_pd(&values[509]);
#if defined(__SSE3__) && defined(__AVX__)
c38_10 = _mm_add_pd(c38_10, _mm_mul_pd(a38_10, _mm256_castpd256_pd128(b38)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c38_10 = _mm_add_pd(c38_10, _mm_mul_pd(a38_10, b38));
#endif
_mm_storeu_pd(&C[(i*88)+67], c38_10);
#else
C[(i*88)+35] += values[499] * B[(i*88)+38];
C[(i*88)+36] += values[500] * B[(i*88)+38];
C[(i*88)+37] += values[501] * B[(i*88)+38];
C[(i*88)+38] += values[502] * B[(i*88)+38];
C[(i*88)+39] += values[503] * B[(i*88)+38];
C[(i*88)+40] += values[504] * B[(i*88)+38];
C[(i*88)+63] += values[505] * B[(i*88)+38];
C[(i*88)+64] += values[506] * B[(i*88)+38];
C[(i*88)+65] += values[507] * B[(i*88)+38];
C[(i*88)+66] += values[508] * B[(i*88)+38];
C[(i*88)+67] += values[509] * B[(i*88)+38];
C[(i*88)+68] += values[510] * B[(i*88)+38];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b39 = _mm256_broadcast_sd(&B[(i*88)+39]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b39 = _mm_loaddup_pd(&B[(i*88)+39]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c39_0 = _mm256_loadu_pd(&C[(i*88)+35]);
__m256d a39_0 = _mm256_loadu_pd(&values[511]);
c39_0 = _mm256_add_pd(c39_0, _mm256_mul_pd(a39_0, b39));
_mm256_storeu_pd(&C[(i*88)+35], c39_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c39_0 = _mm_loadu_pd(&C[(i*88)+35]);
__m128d a39_0 = _mm_loadu_pd(&values[511]);
c39_0 = _mm_add_pd(c39_0, _mm_mul_pd(a39_0, b39));
_mm_storeu_pd(&C[(i*88)+35], c39_0);
__m128d c39_2 = _mm_loadu_pd(&C[(i*88)+37]);
__m128d a39_2 = _mm_loadu_pd(&values[513]);
c39_2 = _mm_add_pd(c39_2, _mm_mul_pd(a39_2, b39));
_mm_storeu_pd(&C[(i*88)+37], c39_2);
#endif
__m128d c39_4 = _mm_loadu_pd(&C[(i*88)+39]);
__m128d a39_4 = _mm_loadu_pd(&values[515]);
#if defined(__SSE3__) && defined(__AVX__)
c39_4 = _mm_add_pd(c39_4, _mm_mul_pd(a39_4, _mm256_castpd256_pd128(b39)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c39_4 = _mm_add_pd(c39_4, _mm_mul_pd(a39_4, b39));
#endif
_mm_storeu_pd(&C[(i*88)+39], c39_4);
#if defined(__SSE3__) && defined(__AVX__)
__m256d c39_6 = _mm256_loadu_pd(&C[(i*88)+63]);
__m256d a39_6 = _mm256_loadu_pd(&values[517]);
c39_6 = _mm256_add_pd(c39_6, _mm256_mul_pd(a39_6, b39));
_mm256_storeu_pd(&C[(i*88)+63], c39_6);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c39_6 = _mm_loadu_pd(&C[(i*88)+63]);
__m128d a39_6 = _mm_loadu_pd(&values[517]);
c39_6 = _mm_add_pd(c39_6, _mm_mul_pd(a39_6, b39));
_mm_storeu_pd(&C[(i*88)+63], c39_6);
__m128d c39_8 = _mm_loadu_pd(&C[(i*88)+65]);
__m128d a39_8 = _mm_loadu_pd(&values[519]);
c39_8 = _mm_add_pd(c39_8, _mm_mul_pd(a39_8, b39));
_mm_storeu_pd(&C[(i*88)+65], c39_8);
#endif
__m128d c39_10 = _mm_loadu_pd(&C[(i*88)+67]);
__m128d a39_10 = _mm_loadu_pd(&values[521]);
#if defined(__SSE3__) && defined(__AVX__)
c39_10 = _mm_add_pd(c39_10, _mm_mul_pd(a39_10, _mm256_castpd256_pd128(b39)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c39_10 = _mm_add_pd(c39_10, _mm_mul_pd(a39_10, b39));
#endif
_mm_storeu_pd(&C[(i*88)+67], c39_10);
#else
C[(i*88)+35] += values[511] * B[(i*88)+39];
C[(i*88)+36] += values[512] * B[(i*88)+39];
C[(i*88)+37] += values[513] * B[(i*88)+39];
C[(i*88)+38] += values[514] * B[(i*88)+39];
C[(i*88)+39] += values[515] * B[(i*88)+39];
C[(i*88)+40] += values[516] * B[(i*88)+39];
C[(i*88)+63] += values[517] * B[(i*88)+39];
C[(i*88)+64] += values[518] * B[(i*88)+39];
C[(i*88)+65] += values[519] * B[(i*88)+39];
C[(i*88)+66] += values[520] * B[(i*88)+39];
C[(i*88)+67] += values[521] * B[(i*88)+39];
C[(i*88)+68] += values[522] * B[(i*88)+39];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b40 = _mm256_broadcast_sd(&B[(i*88)+40]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b40 = _mm_loaddup_pd(&B[(i*88)+40]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c40_0 = _mm256_loadu_pd(&C[(i*88)+35]);
__m256d a40_0 = _mm256_loadu_pd(&values[523]);
c40_0 = _mm256_add_pd(c40_0, _mm256_mul_pd(a40_0, b40));
_mm256_storeu_pd(&C[(i*88)+35], c40_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c40_0 = _mm_loadu_pd(&C[(i*88)+35]);
__m128d a40_0 = _mm_loadu_pd(&values[523]);
c40_0 = _mm_add_pd(c40_0, _mm_mul_pd(a40_0, b40));
_mm_storeu_pd(&C[(i*88)+35], c40_0);
__m128d c40_2 = _mm_loadu_pd(&C[(i*88)+37]);
__m128d a40_2 = _mm_loadu_pd(&values[525]);
c40_2 = _mm_add_pd(c40_2, _mm_mul_pd(a40_2, b40));
_mm_storeu_pd(&C[(i*88)+37], c40_2);
#endif
__m128d c40_4 = _mm_loadu_pd(&C[(i*88)+39]);
__m128d a40_4 = _mm_loadu_pd(&values[527]);
#if defined(__SSE3__) && defined(__AVX__)
c40_4 = _mm_add_pd(c40_4, _mm_mul_pd(a40_4, _mm256_castpd256_pd128(b40)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c40_4 = _mm_add_pd(c40_4, _mm_mul_pd(a40_4, b40));
#endif
_mm_storeu_pd(&C[(i*88)+39], c40_4);
#if defined(__SSE3__) && defined(__AVX__)
__m256d c40_6 = _mm256_loadu_pd(&C[(i*88)+63]);
__m256d a40_6 = _mm256_loadu_pd(&values[529]);
c40_6 = _mm256_add_pd(c40_6, _mm256_mul_pd(a40_6, b40));
_mm256_storeu_pd(&C[(i*88)+63], c40_6);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c40_6 = _mm_loadu_pd(&C[(i*88)+63]);
__m128d a40_6 = _mm_loadu_pd(&values[529]);
c40_6 = _mm_add_pd(c40_6, _mm_mul_pd(a40_6, b40));
_mm_storeu_pd(&C[(i*88)+63], c40_6);
__m128d c40_8 = _mm_loadu_pd(&C[(i*88)+65]);
__m128d a40_8 = _mm_loadu_pd(&values[531]);
c40_8 = _mm_add_pd(c40_8, _mm_mul_pd(a40_8, b40));
_mm_storeu_pd(&C[(i*88)+65], c40_8);
#endif
__m128d c40_10 = _mm_loadu_pd(&C[(i*88)+67]);
__m128d a40_10 = _mm_loadu_pd(&values[533]);
#if defined(__SSE3__) && defined(__AVX__)
c40_10 = _mm_add_pd(c40_10, _mm_mul_pd(a40_10, _mm256_castpd256_pd128(b40)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c40_10 = _mm_add_pd(c40_10, _mm_mul_pd(a40_10, b40));
#endif
_mm_storeu_pd(&C[(i*88)+67], c40_10);
#else
C[(i*88)+35] += values[523] * B[(i*88)+40];
C[(i*88)+36] += values[524] * B[(i*88)+40];
C[(i*88)+37] += values[525] * B[(i*88)+40];
C[(i*88)+38] += values[526] * B[(i*88)+40];
C[(i*88)+39] += values[527] * B[(i*88)+40];
C[(i*88)+40] += values[528] * B[(i*88)+40];
C[(i*88)+63] += values[529] * B[(i*88)+40];
C[(i*88)+64] += values[530] * B[(i*88)+40];
C[(i*88)+65] += values[531] * B[(i*88)+40];
C[(i*88)+66] += values[532] * B[(i*88)+40];
C[(i*88)+67] += values[533] * B[(i*88)+40];
C[(i*88)+68] += values[534] * B[(i*88)+40];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b41 = _mm256_broadcast_sd(&B[(i*88)+41]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b41 = _mm_loaddup_pd(&B[(i*88)+41]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c41_0 = _mm256_loadu_pd(&C[(i*88)+20]);
__m256d a41_0 = _mm256_loadu_pd(&values[535]);
c41_0 = _mm256_add_pd(c41_0, _mm256_mul_pd(a41_0, b41));
_mm256_storeu_pd(&C[(i*88)+20], c41_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c41_0 = _mm_loadu_pd(&C[(i*88)+20]);
__m128d a41_0 = _mm_loadu_pd(&values[535]);
c41_0 = _mm_add_pd(c41_0, _mm_mul_pd(a41_0, b41));
_mm_storeu_pd(&C[(i*88)+20], c41_0);
__m128d c41_2 = _mm_loadu_pd(&C[(i*88)+22]);
__m128d a41_2 = _mm_loadu_pd(&values[537]);
c41_2 = _mm_add_pd(c41_2, _mm_mul_pd(a41_2, b41));
_mm_storeu_pd(&C[(i*88)+22], c41_2);
#endif
__m128d c41_4 = _mm_load_sd(&C[(i*88)+24]);
__m128d a41_4 = _mm_load_sd(&values[539]);
#if defined(__SSE3__) && defined(__AVX__)
c41_4 = _mm_add_sd(c41_4, _mm_mul_sd(a41_4, _mm256_castpd256_pd128(b41)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c41_4 = _mm_add_sd(c41_4, _mm_mul_sd(a41_4, b41));
#endif
_mm_store_sd(&C[(i*88)+24], c41_4);
#if defined(__SSE3__) && defined(__AVX__)
__m256d c41_5 = _mm256_loadu_pd(&C[(i*88)+41]);
__m256d a41_5 = _mm256_loadu_pd(&values[540]);
c41_5 = _mm256_add_pd(c41_5, _mm256_mul_pd(a41_5, b41));
_mm256_storeu_pd(&C[(i*88)+41], c41_5);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c41_5 = _mm_loadu_pd(&C[(i*88)+41]);
__m128d a41_5 = _mm_loadu_pd(&values[540]);
c41_5 = _mm_add_pd(c41_5, _mm_mul_pd(a41_5, b41));
_mm_storeu_pd(&C[(i*88)+41], c41_5);
__m128d c41_7 = _mm_loadu_pd(&C[(i*88)+43]);
__m128d a41_7 = _mm_loadu_pd(&values[542]);
c41_7 = _mm_add_pd(c41_7, _mm_mul_pd(a41_7, b41));
_mm_storeu_pd(&C[(i*88)+43], c41_7);
#endif
__m128d c41_9 = _mm_load_sd(&C[(i*88)+45]);
__m128d a41_9 = _mm_load_sd(&values[544]);
#if defined(__SSE3__) && defined(__AVX__)
c41_9 = _mm_add_sd(c41_9, _mm_mul_sd(a41_9, _mm256_castpd256_pd128(b41)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c41_9 = _mm_add_sd(c41_9, _mm_mul_sd(a41_9, b41));
#endif
_mm_store_sd(&C[(i*88)+45], c41_9);
#if defined(__SSE3__) && defined(__AVX__)
__m256d c41_10 = _mm256_loadu_pd(&C[(i*88)+69]);
__m256d a41_10 = _mm256_loadu_pd(&values[545]);
c41_10 = _mm256_add_pd(c41_10, _mm256_mul_pd(a41_10, b41));
_mm256_storeu_pd(&C[(i*88)+69], c41_10);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c41_10 = _mm_loadu_pd(&C[(i*88)+69]);
__m128d a41_10 = _mm_loadu_pd(&values[545]);
c41_10 = _mm_add_pd(c41_10, _mm_mul_pd(a41_10, b41));
_mm_storeu_pd(&C[(i*88)+69], c41_10);
__m128d c41_12 = _mm_loadu_pd(&C[(i*88)+71]);
__m128d a41_12 = _mm_loadu_pd(&values[547]);
c41_12 = _mm_add_pd(c41_12, _mm_mul_pd(a41_12, b41));
_mm_storeu_pd(&C[(i*88)+71], c41_12);
#endif
__m128d c41_14 = _mm_load_sd(&C[(i*88)+73]);
__m128d a41_14 = _mm_load_sd(&values[549]);
#if defined(__SSE3__) && defined(__AVX__)
c41_14 = _mm_add_sd(c41_14, _mm_mul_sd(a41_14, _mm256_castpd256_pd128(b41)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c41_14 = _mm_add_sd(c41_14, _mm_mul_sd(a41_14, b41));
#endif
_mm_store_sd(&C[(i*88)+73], c41_14);
#else
C[(i*88)+20] += values[535] * B[(i*88)+41];
C[(i*88)+21] += values[536] * B[(i*88)+41];
C[(i*88)+22] += values[537] * B[(i*88)+41];
C[(i*88)+23] += values[538] * B[(i*88)+41];
C[(i*88)+24] += values[539] * B[(i*88)+41];
C[(i*88)+41] += values[540] * B[(i*88)+41];
C[(i*88)+42] += values[541] * B[(i*88)+41];
C[(i*88)+43] += values[542] * B[(i*88)+41];
C[(i*88)+44] += values[543] * B[(i*88)+41];
C[(i*88)+45] += values[544] * B[(i*88)+41];
C[(i*88)+69] += values[545] * B[(i*88)+41];
C[(i*88)+70] += values[546] * B[(i*88)+41];
C[(i*88)+71] += values[547] * B[(i*88)+41];
C[(i*88)+72] += values[548] * B[(i*88)+41];
C[(i*88)+73] += values[549] * B[(i*88)+41];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b42 = _mm256_broadcast_sd(&B[(i*88)+42]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b42 = _mm_loaddup_pd(&B[(i*88)+42]);
#endif
__m128d c42_0 = _mm_loadu_pd(&C[(i*88)+20]);
__m128d a42_0 = _mm_loadu_pd(&values[550]);
#if defined(__SSE3__) && defined(__AVX__)
c42_0 = _mm_add_pd(c42_0, _mm_mul_pd(a42_0, _mm256_castpd256_pd128(b42)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c42_0 = _mm_add_pd(c42_0, _mm_mul_pd(a42_0, b42));
#endif
_mm_storeu_pd(&C[(i*88)+20], c42_0);
__m128d c42_2 = _mm_load_sd(&C[(i*88)+22]);
__m128d a42_2 = _mm_load_sd(&values[552]);
#if defined(__SSE3__) && defined(__AVX__)
c42_2 = _mm_add_sd(c42_2, _mm_mul_sd(a42_2, _mm256_castpd256_pd128(b42)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c42_2 = _mm_add_sd(c42_2, _mm_mul_sd(a42_2, b42));
#endif
_mm_store_sd(&C[(i*88)+22], c42_2);
__m128d c42_3 = _mm_load_sd(&C[(i*88)+24]);
__m128d a42_3 = _mm_load_sd(&values[553]);
#if defined(__SSE3__) && defined(__AVX__)
c42_3 = _mm_add_sd(c42_3, _mm_mul_sd(a42_3, _mm256_castpd256_pd128(b42)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c42_3 = _mm_add_sd(c42_3, _mm_mul_sd(a42_3, b42));
#endif
_mm_store_sd(&C[(i*88)+24], c42_3);
__m128d c42_4 = _mm_loadu_pd(&C[(i*88)+41]);
__m128d a42_4 = _mm_loadu_pd(&values[554]);
#if defined(__SSE3__) && defined(__AVX__)
c42_4 = _mm_add_pd(c42_4, _mm_mul_pd(a42_4, _mm256_castpd256_pd128(b42)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c42_4 = _mm_add_pd(c42_4, _mm_mul_pd(a42_4, b42));
#endif
_mm_storeu_pd(&C[(i*88)+41], c42_4);
__m128d c42_6 = _mm_load_sd(&C[(i*88)+43]);
__m128d a42_6 = _mm_load_sd(&values[556]);
#if defined(__SSE3__) && defined(__AVX__)
c42_6 = _mm_add_sd(c42_6, _mm_mul_sd(a42_6, _mm256_castpd256_pd128(b42)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c42_6 = _mm_add_sd(c42_6, _mm_mul_sd(a42_6, b42));
#endif
_mm_store_sd(&C[(i*88)+43], c42_6);
__m128d c42_7 = _mm_load_sd(&C[(i*88)+45]);
__m128d a42_7 = _mm_load_sd(&values[557]);
#if defined(__SSE3__) && defined(__AVX__)
c42_7 = _mm_add_sd(c42_7, _mm_mul_sd(a42_7, _mm256_castpd256_pd128(b42)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c42_7 = _mm_add_sd(c42_7, _mm_mul_sd(a42_7, b42));
#endif
_mm_store_sd(&C[(i*88)+45], c42_7);
__m128d c42_8 = _mm_loadu_pd(&C[(i*88)+69]);
__m128d a42_8 = _mm_loadu_pd(&values[558]);
#if defined(__SSE3__) && defined(__AVX__)
c42_8 = _mm_add_pd(c42_8, _mm_mul_pd(a42_8, _mm256_castpd256_pd128(b42)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c42_8 = _mm_add_pd(c42_8, _mm_mul_pd(a42_8, b42));
#endif
_mm_storeu_pd(&C[(i*88)+69], c42_8);
__m128d c42_10 = _mm_load_sd(&C[(i*88)+71]);
__m128d a42_10 = _mm_load_sd(&values[560]);
#if defined(__SSE3__) && defined(__AVX__)
c42_10 = _mm_add_sd(c42_10, _mm_mul_sd(a42_10, _mm256_castpd256_pd128(b42)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c42_10 = _mm_add_sd(c42_10, _mm_mul_sd(a42_10, b42));
#endif
_mm_store_sd(&C[(i*88)+71], c42_10);
__m128d c42_11 = _mm_load_sd(&C[(i*88)+73]);
__m128d a42_11 = _mm_load_sd(&values[561]);
#if defined(__SSE3__) && defined(__AVX__)
c42_11 = _mm_add_sd(c42_11, _mm_mul_sd(a42_11, _mm256_castpd256_pd128(b42)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c42_11 = _mm_add_sd(c42_11, _mm_mul_sd(a42_11, b42));
#endif
_mm_store_sd(&C[(i*88)+73], c42_11);
#else
C[(i*88)+20] += values[550] * B[(i*88)+42];
C[(i*88)+21] += values[551] * B[(i*88)+42];
C[(i*88)+22] += values[552] * B[(i*88)+42];
C[(i*88)+24] += values[553] * B[(i*88)+42];
C[(i*88)+41] += values[554] * B[(i*88)+42];
C[(i*88)+42] += values[555] * B[(i*88)+42];
C[(i*88)+43] += values[556] * B[(i*88)+42];
C[(i*88)+45] += values[557] * B[(i*88)+42];
C[(i*88)+69] += values[558] * B[(i*88)+42];
C[(i*88)+70] += values[559] * B[(i*88)+42];
C[(i*88)+71] += values[560] * B[(i*88)+42];
C[(i*88)+73] += values[561] * B[(i*88)+42];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b43 = _mm256_broadcast_sd(&B[(i*88)+43]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b43 = _mm_loaddup_pd(&B[(i*88)+43]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c43_0 = _mm256_loadu_pd(&C[(i*88)+20]);
__m256d a43_0 = _mm256_loadu_pd(&values[562]);
c43_0 = _mm256_add_pd(c43_0, _mm256_mul_pd(a43_0, b43));
_mm256_storeu_pd(&C[(i*88)+20], c43_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c43_0 = _mm_loadu_pd(&C[(i*88)+20]);
__m128d a43_0 = _mm_loadu_pd(&values[562]);
c43_0 = _mm_add_pd(c43_0, _mm_mul_pd(a43_0, b43));
_mm_storeu_pd(&C[(i*88)+20], c43_0);
__m128d c43_2 = _mm_loadu_pd(&C[(i*88)+22]);
__m128d a43_2 = _mm_loadu_pd(&values[564]);
c43_2 = _mm_add_pd(c43_2, _mm_mul_pd(a43_2, b43));
_mm_storeu_pd(&C[(i*88)+22], c43_2);
#endif
__m128d c43_4 = _mm_load_sd(&C[(i*88)+24]);
__m128d a43_4 = _mm_load_sd(&values[566]);
#if defined(__SSE3__) && defined(__AVX__)
c43_4 = _mm_add_sd(c43_4, _mm_mul_sd(a43_4, _mm256_castpd256_pd128(b43)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c43_4 = _mm_add_sd(c43_4, _mm_mul_sd(a43_4, b43));
#endif
_mm_store_sd(&C[(i*88)+24], c43_4);
#if defined(__SSE3__) && defined(__AVX__)
__m256d c43_5 = _mm256_loadu_pd(&C[(i*88)+41]);
__m256d a43_5 = _mm256_loadu_pd(&values[567]);
c43_5 = _mm256_add_pd(c43_5, _mm256_mul_pd(a43_5, b43));
_mm256_storeu_pd(&C[(i*88)+41], c43_5);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c43_5 = _mm_loadu_pd(&C[(i*88)+41]);
__m128d a43_5 = _mm_loadu_pd(&values[567]);
c43_5 = _mm_add_pd(c43_5, _mm_mul_pd(a43_5, b43));
_mm_storeu_pd(&C[(i*88)+41], c43_5);
__m128d c43_7 = _mm_loadu_pd(&C[(i*88)+43]);
__m128d a43_7 = _mm_loadu_pd(&values[569]);
c43_7 = _mm_add_pd(c43_7, _mm_mul_pd(a43_7, b43));
_mm_storeu_pd(&C[(i*88)+43], c43_7);
#endif
__m128d c43_9 = _mm_load_sd(&C[(i*88)+45]);
__m128d a43_9 = _mm_load_sd(&values[571]);
#if defined(__SSE3__) && defined(__AVX__)
c43_9 = _mm_add_sd(c43_9, _mm_mul_sd(a43_9, _mm256_castpd256_pd128(b43)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c43_9 = _mm_add_sd(c43_9, _mm_mul_sd(a43_9, b43));
#endif
_mm_store_sd(&C[(i*88)+45], c43_9);
#if defined(__SSE3__) && defined(__AVX__)
__m256d c43_10 = _mm256_loadu_pd(&C[(i*88)+69]);
__m256d a43_10 = _mm256_loadu_pd(&values[572]);
c43_10 = _mm256_add_pd(c43_10, _mm256_mul_pd(a43_10, b43));
_mm256_storeu_pd(&C[(i*88)+69], c43_10);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c43_10 = _mm_loadu_pd(&C[(i*88)+69]);
__m128d a43_10 = _mm_loadu_pd(&values[572]);
c43_10 = _mm_add_pd(c43_10, _mm_mul_pd(a43_10, b43));
_mm_storeu_pd(&C[(i*88)+69], c43_10);
__m128d c43_12 = _mm_loadu_pd(&C[(i*88)+71]);
__m128d a43_12 = _mm_loadu_pd(&values[574]);
c43_12 = _mm_add_pd(c43_12, _mm_mul_pd(a43_12, b43));
_mm_storeu_pd(&C[(i*88)+71], c43_12);
#endif
__m128d c43_14 = _mm_load_sd(&C[(i*88)+73]);
__m128d a43_14 = _mm_load_sd(&values[576]);
#if defined(__SSE3__) && defined(__AVX__)
c43_14 = _mm_add_sd(c43_14, _mm_mul_sd(a43_14, _mm256_castpd256_pd128(b43)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c43_14 = _mm_add_sd(c43_14, _mm_mul_sd(a43_14, b43));
#endif
_mm_store_sd(&C[(i*88)+73], c43_14);
#else
C[(i*88)+20] += values[562] * B[(i*88)+43];
C[(i*88)+21] += values[563] * B[(i*88)+43];
C[(i*88)+22] += values[564] * B[(i*88)+43];
C[(i*88)+23] += values[565] * B[(i*88)+43];
C[(i*88)+24] += values[566] * B[(i*88)+43];
C[(i*88)+41] += values[567] * B[(i*88)+43];
C[(i*88)+42] += values[568] * B[(i*88)+43];
C[(i*88)+43] += values[569] * B[(i*88)+43];
C[(i*88)+44] += values[570] * B[(i*88)+43];
C[(i*88)+45] += values[571] * B[(i*88)+43];
C[(i*88)+69] += values[572] * B[(i*88)+43];
C[(i*88)+70] += values[573] * B[(i*88)+43];
C[(i*88)+71] += values[574] * B[(i*88)+43];
C[(i*88)+72] += values[575] * B[(i*88)+43];
C[(i*88)+73] += values[576] * B[(i*88)+43];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b44 = _mm256_broadcast_sd(&B[(i*88)+44]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b44 = _mm_loaddup_pd(&B[(i*88)+44]);
#endif
__m128d c44_0 = _mm_load_sd(&C[(i*88)+20]);
__m128d a44_0 = _mm_load_sd(&values[577]);
#if defined(__SSE3__) && defined(__AVX__)
c44_0 = _mm_add_sd(c44_0, _mm_mul_sd(a44_0, _mm256_castpd256_pd128(b44)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c44_0 = _mm_add_sd(c44_0, _mm_mul_sd(a44_0, b44));
#endif
_mm_store_sd(&C[(i*88)+20], c44_0);
__m128d c44_1 = _mm_loadu_pd(&C[(i*88)+22]);
__m128d a44_1 = _mm_loadu_pd(&values[578]);
#if defined(__SSE3__) && defined(__AVX__)
c44_1 = _mm_add_pd(c44_1, _mm_mul_pd(a44_1, _mm256_castpd256_pd128(b44)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c44_1 = _mm_add_pd(c44_1, _mm_mul_pd(a44_1, b44));
#endif
_mm_storeu_pd(&C[(i*88)+22], c44_1);
__m128d c44_3 = _mm_load_sd(&C[(i*88)+24]);
__m128d a44_3 = _mm_load_sd(&values[580]);
#if defined(__SSE3__) && defined(__AVX__)
c44_3 = _mm_add_sd(c44_3, _mm_mul_sd(a44_3, _mm256_castpd256_pd128(b44)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c44_3 = _mm_add_sd(c44_3, _mm_mul_sd(a44_3, b44));
#endif
_mm_store_sd(&C[(i*88)+24], c44_3);
__m128d c44_4 = _mm_load_sd(&C[(i*88)+41]);
__m128d a44_4 = _mm_load_sd(&values[581]);
#if defined(__SSE3__) && defined(__AVX__)
c44_4 = _mm_add_sd(c44_4, _mm_mul_sd(a44_4, _mm256_castpd256_pd128(b44)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c44_4 = _mm_add_sd(c44_4, _mm_mul_sd(a44_4, b44));
#endif
_mm_store_sd(&C[(i*88)+41], c44_4);
__m128d c44_5 = _mm_loadu_pd(&C[(i*88)+43]);
__m128d a44_5 = _mm_loadu_pd(&values[582]);
#if defined(__SSE3__) && defined(__AVX__)
c44_5 = _mm_add_pd(c44_5, _mm_mul_pd(a44_5, _mm256_castpd256_pd128(b44)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c44_5 = _mm_add_pd(c44_5, _mm_mul_pd(a44_5, b44));
#endif
_mm_storeu_pd(&C[(i*88)+43], c44_5);
__m128d c44_7 = _mm_load_sd(&C[(i*88)+45]);
__m128d a44_7 = _mm_load_sd(&values[584]);
#if defined(__SSE3__) && defined(__AVX__)
c44_7 = _mm_add_sd(c44_7, _mm_mul_sd(a44_7, _mm256_castpd256_pd128(b44)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c44_7 = _mm_add_sd(c44_7, _mm_mul_sd(a44_7, b44));
#endif
_mm_store_sd(&C[(i*88)+45], c44_7);
__m128d c44_8 = _mm_load_sd(&C[(i*88)+69]);
__m128d a44_8 = _mm_load_sd(&values[585]);
#if defined(__SSE3__) && defined(__AVX__)
c44_8 = _mm_add_sd(c44_8, _mm_mul_sd(a44_8, _mm256_castpd256_pd128(b44)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c44_8 = _mm_add_sd(c44_8, _mm_mul_sd(a44_8, b44));
#endif
_mm_store_sd(&C[(i*88)+69], c44_8);
__m128d c44_9 = _mm_loadu_pd(&C[(i*88)+71]);
__m128d a44_9 = _mm_loadu_pd(&values[586]);
#if defined(__SSE3__) && defined(__AVX__)
c44_9 = _mm_add_pd(c44_9, _mm_mul_pd(a44_9, _mm256_castpd256_pd128(b44)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c44_9 = _mm_add_pd(c44_9, _mm_mul_pd(a44_9, b44));
#endif
_mm_storeu_pd(&C[(i*88)+71], c44_9);
__m128d c44_11 = _mm_load_sd(&C[(i*88)+73]);
__m128d a44_11 = _mm_load_sd(&values[588]);
#if defined(__SSE3__) && defined(__AVX__)
c44_11 = _mm_add_sd(c44_11, _mm_mul_sd(a44_11, _mm256_castpd256_pd128(b44)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c44_11 = _mm_add_sd(c44_11, _mm_mul_sd(a44_11, b44));
#endif
_mm_store_sd(&C[(i*88)+73], c44_11);
#else
C[(i*88)+20] += values[577] * B[(i*88)+44];
C[(i*88)+22] += values[578] * B[(i*88)+44];
C[(i*88)+23] += values[579] * B[(i*88)+44];
C[(i*88)+24] += values[580] * B[(i*88)+44];
C[(i*88)+41] += values[581] * B[(i*88)+44];
C[(i*88)+43] += values[582] * B[(i*88)+44];
C[(i*88)+44] += values[583] * B[(i*88)+44];
C[(i*88)+45] += values[584] * B[(i*88)+44];
C[(i*88)+69] += values[585] * B[(i*88)+44];
C[(i*88)+71] += values[586] * B[(i*88)+44];
C[(i*88)+72] += values[587] * B[(i*88)+44];
C[(i*88)+73] += values[588] * B[(i*88)+44];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b45 = _mm256_broadcast_sd(&B[(i*88)+45]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b45 = _mm_loaddup_pd(&B[(i*88)+45]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c45_0 = _mm256_loadu_pd(&C[(i*88)+20]);
__m256d a45_0 = _mm256_loadu_pd(&values[589]);
c45_0 = _mm256_add_pd(c45_0, _mm256_mul_pd(a45_0, b45));
_mm256_storeu_pd(&C[(i*88)+20], c45_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c45_0 = _mm_loadu_pd(&C[(i*88)+20]);
__m128d a45_0 = _mm_loadu_pd(&values[589]);
c45_0 = _mm_add_pd(c45_0, _mm_mul_pd(a45_0, b45));
_mm_storeu_pd(&C[(i*88)+20], c45_0);
__m128d c45_2 = _mm_loadu_pd(&C[(i*88)+22]);
__m128d a45_2 = _mm_loadu_pd(&values[591]);
c45_2 = _mm_add_pd(c45_2, _mm_mul_pd(a45_2, b45));
_mm_storeu_pd(&C[(i*88)+22], c45_2);
#endif
__m128d c45_4 = _mm_load_sd(&C[(i*88)+24]);
__m128d a45_4 = _mm_load_sd(&values[593]);
#if defined(__SSE3__) && defined(__AVX__)
c45_4 = _mm_add_sd(c45_4, _mm_mul_sd(a45_4, _mm256_castpd256_pd128(b45)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c45_4 = _mm_add_sd(c45_4, _mm_mul_sd(a45_4, b45));
#endif
_mm_store_sd(&C[(i*88)+24], c45_4);
#if defined(__SSE3__) && defined(__AVX__)
__m256d c45_5 = _mm256_loadu_pd(&C[(i*88)+41]);
__m256d a45_5 = _mm256_loadu_pd(&values[594]);
c45_5 = _mm256_add_pd(c45_5, _mm256_mul_pd(a45_5, b45));
_mm256_storeu_pd(&C[(i*88)+41], c45_5);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c45_5 = _mm_loadu_pd(&C[(i*88)+41]);
__m128d a45_5 = _mm_loadu_pd(&values[594]);
c45_5 = _mm_add_pd(c45_5, _mm_mul_pd(a45_5, b45));
_mm_storeu_pd(&C[(i*88)+41], c45_5);
__m128d c45_7 = _mm_loadu_pd(&C[(i*88)+43]);
__m128d a45_7 = _mm_loadu_pd(&values[596]);
c45_7 = _mm_add_pd(c45_7, _mm_mul_pd(a45_7, b45));
_mm_storeu_pd(&C[(i*88)+43], c45_7);
#endif
__m128d c45_9 = _mm_load_sd(&C[(i*88)+45]);
__m128d a45_9 = _mm_load_sd(&values[598]);
#if defined(__SSE3__) && defined(__AVX__)
c45_9 = _mm_add_sd(c45_9, _mm_mul_sd(a45_9, _mm256_castpd256_pd128(b45)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c45_9 = _mm_add_sd(c45_9, _mm_mul_sd(a45_9, b45));
#endif
_mm_store_sd(&C[(i*88)+45], c45_9);
#if defined(__SSE3__) && defined(__AVX__)
__m256d c45_10 = _mm256_loadu_pd(&C[(i*88)+69]);
__m256d a45_10 = _mm256_loadu_pd(&values[599]);
c45_10 = _mm256_add_pd(c45_10, _mm256_mul_pd(a45_10, b45));
_mm256_storeu_pd(&C[(i*88)+69], c45_10);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c45_10 = _mm_loadu_pd(&C[(i*88)+69]);
__m128d a45_10 = _mm_loadu_pd(&values[599]);
c45_10 = _mm_add_pd(c45_10, _mm_mul_pd(a45_10, b45));
_mm_storeu_pd(&C[(i*88)+69], c45_10);
__m128d c45_12 = _mm_loadu_pd(&C[(i*88)+71]);
__m128d a45_12 = _mm_loadu_pd(&values[601]);
c45_12 = _mm_add_pd(c45_12, _mm_mul_pd(a45_12, b45));
_mm_storeu_pd(&C[(i*88)+71], c45_12);
#endif
__m128d c45_14 = _mm_load_sd(&C[(i*88)+73]);
__m128d a45_14 = _mm_load_sd(&values[603]);
#if defined(__SSE3__) && defined(__AVX__)
c45_14 = _mm_add_sd(c45_14, _mm_mul_sd(a45_14, _mm256_castpd256_pd128(b45)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c45_14 = _mm_add_sd(c45_14, _mm_mul_sd(a45_14, b45));
#endif
_mm_store_sd(&C[(i*88)+73], c45_14);
#else
C[(i*88)+20] += values[589] * B[(i*88)+45];
C[(i*88)+21] += values[590] * B[(i*88)+45];
C[(i*88)+22] += values[591] * B[(i*88)+45];
C[(i*88)+23] += values[592] * B[(i*88)+45];
C[(i*88)+24] += values[593] * B[(i*88)+45];
C[(i*88)+41] += values[594] * B[(i*88)+45];
C[(i*88)+42] += values[595] * B[(i*88)+45];
C[(i*88)+43] += values[596] * B[(i*88)+45];
C[(i*88)+44] += values[597] * B[(i*88)+45];
C[(i*88)+45] += values[598] * B[(i*88)+45];
C[(i*88)+69] += values[599] * B[(i*88)+45];
C[(i*88)+70] += values[600] * B[(i*88)+45];
C[(i*88)+71] += values[601] * B[(i*88)+45];
C[(i*88)+72] += values[602] * B[(i*88)+45];
C[(i*88)+73] += values[603] * B[(i*88)+45];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b46 = _mm256_broadcast_sd(&B[(i*88)+46]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b46 = _mm_loaddup_pd(&B[(i*88)+46]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c46_0 = _mm256_loadu_pd(&C[(i*88)+10]);
__m256d a46_0 = _mm256_loadu_pd(&values[604]);
c46_0 = _mm256_add_pd(c46_0, _mm256_mul_pd(a46_0, b46));
_mm256_storeu_pd(&C[(i*88)+10], c46_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c46_0 = _mm_loadu_pd(&C[(i*88)+10]);
__m128d a46_0 = _mm_loadu_pd(&values[604]);
c46_0 = _mm_add_pd(c46_0, _mm_mul_pd(a46_0, b46));
_mm_storeu_pd(&C[(i*88)+10], c46_0);
__m128d c46_2 = _mm_loadu_pd(&C[(i*88)+12]);
__m128d a46_2 = _mm_loadu_pd(&values[606]);
c46_2 = _mm_add_pd(c46_2, _mm_mul_pd(a46_2, b46));
_mm_storeu_pd(&C[(i*88)+12], c46_2);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c46_4 = _mm256_loadu_pd(&C[(i*88)+25]);
__m256d a46_4 = _mm256_loadu_pd(&values[608]);
c46_4 = _mm256_add_pd(c46_4, _mm256_mul_pd(a46_4, b46));
_mm256_storeu_pd(&C[(i*88)+25], c46_4);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c46_4 = _mm_loadu_pd(&C[(i*88)+25]);
__m128d a46_4 = _mm_loadu_pd(&values[608]);
c46_4 = _mm_add_pd(c46_4, _mm_mul_pd(a46_4, b46));
_mm_storeu_pd(&C[(i*88)+25], c46_4);
__m128d c46_6 = _mm_loadu_pd(&C[(i*88)+27]);
__m128d a46_6 = _mm_loadu_pd(&values[610]);
c46_6 = _mm_add_pd(c46_6, _mm_mul_pd(a46_6, b46));
_mm_storeu_pd(&C[(i*88)+27], c46_6);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c46_8 = _mm256_loadu_pd(&C[(i*88)+46]);
__m256d a46_8 = _mm256_loadu_pd(&values[612]);
c46_8 = _mm256_add_pd(c46_8, _mm256_mul_pd(a46_8, b46));
_mm256_storeu_pd(&C[(i*88)+46], c46_8);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c46_8 = _mm_loadu_pd(&C[(i*88)+46]);
__m128d a46_8 = _mm_loadu_pd(&values[612]);
c46_8 = _mm_add_pd(c46_8, _mm_mul_pd(a46_8, b46));
_mm_storeu_pd(&C[(i*88)+46], c46_8);
__m128d c46_10 = _mm_loadu_pd(&C[(i*88)+48]);
__m128d a46_10 = _mm_loadu_pd(&values[614]);
c46_10 = _mm_add_pd(c46_10, _mm_mul_pd(a46_10, b46));
_mm_storeu_pd(&C[(i*88)+48], c46_10);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c46_12 = _mm256_loadu_pd(&C[(i*88)+74]);
__m256d a46_12 = _mm256_loadu_pd(&values[616]);
c46_12 = _mm256_add_pd(c46_12, _mm256_mul_pd(a46_12, b46));
_mm256_storeu_pd(&C[(i*88)+74], c46_12);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c46_12 = _mm_loadu_pd(&C[(i*88)+74]);
__m128d a46_12 = _mm_loadu_pd(&values[616]);
c46_12 = _mm_add_pd(c46_12, _mm_mul_pd(a46_12, b46));
_mm_storeu_pd(&C[(i*88)+74], c46_12);
__m128d c46_14 = _mm_loadu_pd(&C[(i*88)+76]);
__m128d a46_14 = _mm_loadu_pd(&values[618]);
c46_14 = _mm_add_pd(c46_14, _mm_mul_pd(a46_14, b46));
_mm_storeu_pd(&C[(i*88)+76], c46_14);
#endif
#else
C[(i*88)+10] += values[604] * B[(i*88)+46];
C[(i*88)+11] += values[605] * B[(i*88)+46];
C[(i*88)+12] += values[606] * B[(i*88)+46];
C[(i*88)+13] += values[607] * B[(i*88)+46];
C[(i*88)+25] += values[608] * B[(i*88)+46];
C[(i*88)+26] += values[609] * B[(i*88)+46];
C[(i*88)+27] += values[610] * B[(i*88)+46];
C[(i*88)+28] += values[611] * B[(i*88)+46];
C[(i*88)+46] += values[612] * B[(i*88)+46];
C[(i*88)+47] += values[613] * B[(i*88)+46];
C[(i*88)+48] += values[614] * B[(i*88)+46];
C[(i*88)+49] += values[615] * B[(i*88)+46];
C[(i*88)+74] += values[616] * B[(i*88)+46];
C[(i*88)+75] += values[617] * B[(i*88)+46];
C[(i*88)+76] += values[618] * B[(i*88)+46];
C[(i*88)+77] += values[619] * B[(i*88)+46];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b47 = _mm256_broadcast_sd(&B[(i*88)+47]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b47 = _mm_loaddup_pd(&B[(i*88)+47]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c47_0 = _mm256_loadu_pd(&C[(i*88)+10]);
__m256d a47_0 = _mm256_loadu_pd(&values[620]);
c47_0 = _mm256_add_pd(c47_0, _mm256_mul_pd(a47_0, b47));
_mm256_storeu_pd(&C[(i*88)+10], c47_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c47_0 = _mm_loadu_pd(&C[(i*88)+10]);
__m128d a47_0 = _mm_loadu_pd(&values[620]);
c47_0 = _mm_add_pd(c47_0, _mm_mul_pd(a47_0, b47));
_mm_storeu_pd(&C[(i*88)+10], c47_0);
__m128d c47_2 = _mm_loadu_pd(&C[(i*88)+12]);
__m128d a47_2 = _mm_loadu_pd(&values[622]);
c47_2 = _mm_add_pd(c47_2, _mm_mul_pd(a47_2, b47));
_mm_storeu_pd(&C[(i*88)+12], c47_2);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c47_4 = _mm256_loadu_pd(&C[(i*88)+25]);
__m256d a47_4 = _mm256_loadu_pd(&values[624]);
c47_4 = _mm256_add_pd(c47_4, _mm256_mul_pd(a47_4, b47));
_mm256_storeu_pd(&C[(i*88)+25], c47_4);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c47_4 = _mm_loadu_pd(&C[(i*88)+25]);
__m128d a47_4 = _mm_loadu_pd(&values[624]);
c47_4 = _mm_add_pd(c47_4, _mm_mul_pd(a47_4, b47));
_mm_storeu_pd(&C[(i*88)+25], c47_4);
__m128d c47_6 = _mm_loadu_pd(&C[(i*88)+27]);
__m128d a47_6 = _mm_loadu_pd(&values[626]);
c47_6 = _mm_add_pd(c47_6, _mm_mul_pd(a47_6, b47));
_mm_storeu_pd(&C[(i*88)+27], c47_6);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c47_8 = _mm256_loadu_pd(&C[(i*88)+46]);
__m256d a47_8 = _mm256_loadu_pd(&values[628]);
c47_8 = _mm256_add_pd(c47_8, _mm256_mul_pd(a47_8, b47));
_mm256_storeu_pd(&C[(i*88)+46], c47_8);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c47_8 = _mm_loadu_pd(&C[(i*88)+46]);
__m128d a47_8 = _mm_loadu_pd(&values[628]);
c47_8 = _mm_add_pd(c47_8, _mm_mul_pd(a47_8, b47));
_mm_storeu_pd(&C[(i*88)+46], c47_8);
__m128d c47_10 = _mm_loadu_pd(&C[(i*88)+48]);
__m128d a47_10 = _mm_loadu_pd(&values[630]);
c47_10 = _mm_add_pd(c47_10, _mm_mul_pd(a47_10, b47));
_mm_storeu_pd(&C[(i*88)+48], c47_10);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c47_12 = _mm256_loadu_pd(&C[(i*88)+74]);
__m256d a47_12 = _mm256_loadu_pd(&values[632]);
c47_12 = _mm256_add_pd(c47_12, _mm256_mul_pd(a47_12, b47));
_mm256_storeu_pd(&C[(i*88)+74], c47_12);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c47_12 = _mm_loadu_pd(&C[(i*88)+74]);
__m128d a47_12 = _mm_loadu_pd(&values[632]);
c47_12 = _mm_add_pd(c47_12, _mm_mul_pd(a47_12, b47));
_mm_storeu_pd(&C[(i*88)+74], c47_12);
__m128d c47_14 = _mm_loadu_pd(&C[(i*88)+76]);
__m128d a47_14 = _mm_loadu_pd(&values[634]);
c47_14 = _mm_add_pd(c47_14, _mm_mul_pd(a47_14, b47));
_mm_storeu_pd(&C[(i*88)+76], c47_14);
#endif
#else
C[(i*88)+10] += values[620] * B[(i*88)+47];
C[(i*88)+11] += values[621] * B[(i*88)+47];
C[(i*88)+12] += values[622] * B[(i*88)+47];
C[(i*88)+13] += values[623] * B[(i*88)+47];
C[(i*88)+25] += values[624] * B[(i*88)+47];
C[(i*88)+26] += values[625] * B[(i*88)+47];
C[(i*88)+27] += values[626] * B[(i*88)+47];
C[(i*88)+28] += values[627] * B[(i*88)+47];
C[(i*88)+46] += values[628] * B[(i*88)+47];
C[(i*88)+47] += values[629] * B[(i*88)+47];
C[(i*88)+48] += values[630] * B[(i*88)+47];
C[(i*88)+49] += values[631] * B[(i*88)+47];
C[(i*88)+74] += values[632] * B[(i*88)+47];
C[(i*88)+75] += values[633] * B[(i*88)+47];
C[(i*88)+76] += values[634] * B[(i*88)+47];
C[(i*88)+77] += values[635] * B[(i*88)+47];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b48 = _mm256_broadcast_sd(&B[(i*88)+48]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b48 = _mm_loaddup_pd(&B[(i*88)+48]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c48_0 = _mm256_loadu_pd(&C[(i*88)+10]);
__m256d a48_0 = _mm256_loadu_pd(&values[636]);
c48_0 = _mm256_add_pd(c48_0, _mm256_mul_pd(a48_0, b48));
_mm256_storeu_pd(&C[(i*88)+10], c48_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c48_0 = _mm_loadu_pd(&C[(i*88)+10]);
__m128d a48_0 = _mm_loadu_pd(&values[636]);
c48_0 = _mm_add_pd(c48_0, _mm_mul_pd(a48_0, b48));
_mm_storeu_pd(&C[(i*88)+10], c48_0);
__m128d c48_2 = _mm_loadu_pd(&C[(i*88)+12]);
__m128d a48_2 = _mm_loadu_pd(&values[638]);
c48_2 = _mm_add_pd(c48_2, _mm_mul_pd(a48_2, b48));
_mm_storeu_pd(&C[(i*88)+12], c48_2);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c48_4 = _mm256_loadu_pd(&C[(i*88)+25]);
__m256d a48_4 = _mm256_loadu_pd(&values[640]);
c48_4 = _mm256_add_pd(c48_4, _mm256_mul_pd(a48_4, b48));
_mm256_storeu_pd(&C[(i*88)+25], c48_4);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c48_4 = _mm_loadu_pd(&C[(i*88)+25]);
__m128d a48_4 = _mm_loadu_pd(&values[640]);
c48_4 = _mm_add_pd(c48_4, _mm_mul_pd(a48_4, b48));
_mm_storeu_pd(&C[(i*88)+25], c48_4);
__m128d c48_6 = _mm_loadu_pd(&C[(i*88)+27]);
__m128d a48_6 = _mm_loadu_pd(&values[642]);
c48_6 = _mm_add_pd(c48_6, _mm_mul_pd(a48_6, b48));
_mm_storeu_pd(&C[(i*88)+27], c48_6);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c48_8 = _mm256_loadu_pd(&C[(i*88)+46]);
__m256d a48_8 = _mm256_loadu_pd(&values[644]);
c48_8 = _mm256_add_pd(c48_8, _mm256_mul_pd(a48_8, b48));
_mm256_storeu_pd(&C[(i*88)+46], c48_8);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c48_8 = _mm_loadu_pd(&C[(i*88)+46]);
__m128d a48_8 = _mm_loadu_pd(&values[644]);
c48_8 = _mm_add_pd(c48_8, _mm_mul_pd(a48_8, b48));
_mm_storeu_pd(&C[(i*88)+46], c48_8);
__m128d c48_10 = _mm_loadu_pd(&C[(i*88)+48]);
__m128d a48_10 = _mm_loadu_pd(&values[646]);
c48_10 = _mm_add_pd(c48_10, _mm_mul_pd(a48_10, b48));
_mm_storeu_pd(&C[(i*88)+48], c48_10);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c48_12 = _mm256_loadu_pd(&C[(i*88)+74]);
__m256d a48_12 = _mm256_loadu_pd(&values[648]);
c48_12 = _mm256_add_pd(c48_12, _mm256_mul_pd(a48_12, b48));
_mm256_storeu_pd(&C[(i*88)+74], c48_12);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c48_12 = _mm_loadu_pd(&C[(i*88)+74]);
__m128d a48_12 = _mm_loadu_pd(&values[648]);
c48_12 = _mm_add_pd(c48_12, _mm_mul_pd(a48_12, b48));
_mm_storeu_pd(&C[(i*88)+74], c48_12);
__m128d c48_14 = _mm_loadu_pd(&C[(i*88)+76]);
__m128d a48_14 = _mm_loadu_pd(&values[650]);
c48_14 = _mm_add_pd(c48_14, _mm_mul_pd(a48_14, b48));
_mm_storeu_pd(&C[(i*88)+76], c48_14);
#endif
#else
C[(i*88)+10] += values[636] * B[(i*88)+48];
C[(i*88)+11] += values[637] * B[(i*88)+48];
C[(i*88)+12] += values[638] * B[(i*88)+48];
C[(i*88)+13] += values[639] * B[(i*88)+48];
C[(i*88)+25] += values[640] * B[(i*88)+48];
C[(i*88)+26] += values[641] * B[(i*88)+48];
C[(i*88)+27] += values[642] * B[(i*88)+48];
C[(i*88)+28] += values[643] * B[(i*88)+48];
C[(i*88)+46] += values[644] * B[(i*88)+48];
C[(i*88)+47] += values[645] * B[(i*88)+48];
C[(i*88)+48] += values[646] * B[(i*88)+48];
C[(i*88)+49] += values[647] * B[(i*88)+48];
C[(i*88)+74] += values[648] * B[(i*88)+48];
C[(i*88)+75] += values[649] * B[(i*88)+48];
C[(i*88)+76] += values[650] * B[(i*88)+48];
C[(i*88)+77] += values[651] * B[(i*88)+48];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b49 = _mm256_broadcast_sd(&B[(i*88)+49]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b49 = _mm_loaddup_pd(&B[(i*88)+49]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c49_0 = _mm256_loadu_pd(&C[(i*88)+10]);
__m256d a49_0 = _mm256_loadu_pd(&values[652]);
c49_0 = _mm256_add_pd(c49_0, _mm256_mul_pd(a49_0, b49));
_mm256_storeu_pd(&C[(i*88)+10], c49_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c49_0 = _mm_loadu_pd(&C[(i*88)+10]);
__m128d a49_0 = _mm_loadu_pd(&values[652]);
c49_0 = _mm_add_pd(c49_0, _mm_mul_pd(a49_0, b49));
_mm_storeu_pd(&C[(i*88)+10], c49_0);
__m128d c49_2 = _mm_loadu_pd(&C[(i*88)+12]);
__m128d a49_2 = _mm_loadu_pd(&values[654]);
c49_2 = _mm_add_pd(c49_2, _mm_mul_pd(a49_2, b49));
_mm_storeu_pd(&C[(i*88)+12], c49_2);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c49_4 = _mm256_loadu_pd(&C[(i*88)+25]);
__m256d a49_4 = _mm256_loadu_pd(&values[656]);
c49_4 = _mm256_add_pd(c49_4, _mm256_mul_pd(a49_4, b49));
_mm256_storeu_pd(&C[(i*88)+25], c49_4);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c49_4 = _mm_loadu_pd(&C[(i*88)+25]);
__m128d a49_4 = _mm_loadu_pd(&values[656]);
c49_4 = _mm_add_pd(c49_4, _mm_mul_pd(a49_4, b49));
_mm_storeu_pd(&C[(i*88)+25], c49_4);
__m128d c49_6 = _mm_loadu_pd(&C[(i*88)+27]);
__m128d a49_6 = _mm_loadu_pd(&values[658]);
c49_6 = _mm_add_pd(c49_6, _mm_mul_pd(a49_6, b49));
_mm_storeu_pd(&C[(i*88)+27], c49_6);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c49_8 = _mm256_loadu_pd(&C[(i*88)+46]);
__m256d a49_8 = _mm256_loadu_pd(&values[660]);
c49_8 = _mm256_add_pd(c49_8, _mm256_mul_pd(a49_8, b49));
_mm256_storeu_pd(&C[(i*88)+46], c49_8);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c49_8 = _mm_loadu_pd(&C[(i*88)+46]);
__m128d a49_8 = _mm_loadu_pd(&values[660]);
c49_8 = _mm_add_pd(c49_8, _mm_mul_pd(a49_8, b49));
_mm_storeu_pd(&C[(i*88)+46], c49_8);
__m128d c49_10 = _mm_loadu_pd(&C[(i*88)+48]);
__m128d a49_10 = _mm_loadu_pd(&values[662]);
c49_10 = _mm_add_pd(c49_10, _mm_mul_pd(a49_10, b49));
_mm_storeu_pd(&C[(i*88)+48], c49_10);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c49_12 = _mm256_loadu_pd(&C[(i*88)+74]);
__m256d a49_12 = _mm256_loadu_pd(&values[664]);
c49_12 = _mm256_add_pd(c49_12, _mm256_mul_pd(a49_12, b49));
_mm256_storeu_pd(&C[(i*88)+74], c49_12);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c49_12 = _mm_loadu_pd(&C[(i*88)+74]);
__m128d a49_12 = _mm_loadu_pd(&values[664]);
c49_12 = _mm_add_pd(c49_12, _mm_mul_pd(a49_12, b49));
_mm_storeu_pd(&C[(i*88)+74], c49_12);
__m128d c49_14 = _mm_loadu_pd(&C[(i*88)+76]);
__m128d a49_14 = _mm_loadu_pd(&values[666]);
c49_14 = _mm_add_pd(c49_14, _mm_mul_pd(a49_14, b49));
_mm_storeu_pd(&C[(i*88)+76], c49_14);
#endif
#else
C[(i*88)+10] += values[652] * B[(i*88)+49];
C[(i*88)+11] += values[653] * B[(i*88)+49];
C[(i*88)+12] += values[654] * B[(i*88)+49];
C[(i*88)+13] += values[655] * B[(i*88)+49];
C[(i*88)+25] += values[656] * B[(i*88)+49];
C[(i*88)+26] += values[657] * B[(i*88)+49];
C[(i*88)+27] += values[658] * B[(i*88)+49];
C[(i*88)+28] += values[659] * B[(i*88)+49];
C[(i*88)+46] += values[660] * B[(i*88)+49];
C[(i*88)+47] += values[661] * B[(i*88)+49];
C[(i*88)+48] += values[662] * B[(i*88)+49];
C[(i*88)+49] += values[663] * B[(i*88)+49];
C[(i*88)+74] += values[664] * B[(i*88)+49];
C[(i*88)+75] += values[665] * B[(i*88)+49];
C[(i*88)+76] += values[666] * B[(i*88)+49];
C[(i*88)+77] += values[667] * B[(i*88)+49];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b50 = _mm256_broadcast_sd(&B[(i*88)+50]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b50 = _mm_loaddup_pd(&B[(i*88)+50]);
#endif
__m128d c50_0 = _mm_loadu_pd(&C[(i*88)+4]);
__m128d a50_0 = _mm_loadu_pd(&values[668]);
#if defined(__SSE3__) && defined(__AVX__)
c50_0 = _mm_add_pd(c50_0, _mm_mul_pd(a50_0, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c50_0 = _mm_add_pd(c50_0, _mm_mul_pd(a50_0, b50));
#endif
_mm_storeu_pd(&C[(i*88)+4], c50_0);
__m128d c50_2 = _mm_load_sd(&C[(i*88)+6]);
__m128d a50_2 = _mm_load_sd(&values[670]);
#if defined(__SSE3__) && defined(__AVX__)
c50_2 = _mm_add_sd(c50_2, _mm_mul_sd(a50_2, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c50_2 = _mm_add_sd(c50_2, _mm_mul_sd(a50_2, b50));
#endif
_mm_store_sd(&C[(i*88)+6], c50_2);
__m128d c50_3 = _mm_loadu_pd(&C[(i*88)+14]);
__m128d a50_3 = _mm_loadu_pd(&values[671]);
#if defined(__SSE3__) && defined(__AVX__)
c50_3 = _mm_add_pd(c50_3, _mm_mul_pd(a50_3, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c50_3 = _mm_add_pd(c50_3, _mm_mul_pd(a50_3, b50));
#endif
_mm_storeu_pd(&C[(i*88)+14], c50_3);
__m128d c50_5 = _mm_load_sd(&C[(i*88)+16]);
__m128d a50_5 = _mm_load_sd(&values[673]);
#if defined(__SSE3__) && defined(__AVX__)
c50_5 = _mm_add_sd(c50_5, _mm_mul_sd(a50_5, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c50_5 = _mm_add_sd(c50_5, _mm_mul_sd(a50_5, b50));
#endif
_mm_store_sd(&C[(i*88)+16], c50_5);
__m128d c50_6 = _mm_loadu_pd(&C[(i*88)+29]);
__m128d a50_6 = _mm_loadu_pd(&values[674]);
#if defined(__SSE3__) && defined(__AVX__)
c50_6 = _mm_add_pd(c50_6, _mm_mul_pd(a50_6, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c50_6 = _mm_add_pd(c50_6, _mm_mul_pd(a50_6, b50));
#endif
_mm_storeu_pd(&C[(i*88)+29], c50_6);
__m128d c50_8 = _mm_load_sd(&C[(i*88)+31]);
__m128d a50_8 = _mm_load_sd(&values[676]);
#if defined(__SSE3__) && defined(__AVX__)
c50_8 = _mm_add_sd(c50_8, _mm_mul_sd(a50_8, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c50_8 = _mm_add_sd(c50_8, _mm_mul_sd(a50_8, b50));
#endif
_mm_store_sd(&C[(i*88)+31], c50_8);
__m128d c50_9 = _mm_loadu_pd(&C[(i*88)+50]);
__m128d a50_9 = _mm_loadu_pd(&values[677]);
#if defined(__SSE3__) && defined(__AVX__)
c50_9 = _mm_add_pd(c50_9, _mm_mul_pd(a50_9, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c50_9 = _mm_add_pd(c50_9, _mm_mul_pd(a50_9, b50));
#endif
_mm_storeu_pd(&C[(i*88)+50], c50_9);
__m128d c50_11 = _mm_load_sd(&C[(i*88)+52]);
__m128d a50_11 = _mm_load_sd(&values[679]);
#if defined(__SSE3__) && defined(__AVX__)
c50_11 = _mm_add_sd(c50_11, _mm_mul_sd(a50_11, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c50_11 = _mm_add_sd(c50_11, _mm_mul_sd(a50_11, b50));
#endif
_mm_store_sd(&C[(i*88)+52], c50_11);
__m128d c50_12 = _mm_loadu_pd(&C[(i*88)+78]);
__m128d a50_12 = _mm_loadu_pd(&values[680]);
#if defined(__SSE3__) && defined(__AVX__)
c50_12 = _mm_add_pd(c50_12, _mm_mul_pd(a50_12, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c50_12 = _mm_add_pd(c50_12, _mm_mul_pd(a50_12, b50));
#endif
_mm_storeu_pd(&C[(i*88)+78], c50_12);
__m128d c50_14 = _mm_load_sd(&C[(i*88)+80]);
__m128d a50_14 = _mm_load_sd(&values[682]);
#if defined(__SSE3__) && defined(__AVX__)
c50_14 = _mm_add_sd(c50_14, _mm_mul_sd(a50_14, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c50_14 = _mm_add_sd(c50_14, _mm_mul_sd(a50_14, b50));
#endif
_mm_store_sd(&C[(i*88)+80], c50_14);
#else
C[(i*88)+4] += values[668] * B[(i*88)+50];
C[(i*88)+5] += values[669] * B[(i*88)+50];
C[(i*88)+6] += values[670] * B[(i*88)+50];
C[(i*88)+14] += values[671] * B[(i*88)+50];
C[(i*88)+15] += values[672] * B[(i*88)+50];
C[(i*88)+16] += values[673] * B[(i*88)+50];
C[(i*88)+29] += values[674] * B[(i*88)+50];
C[(i*88)+30] += values[675] * B[(i*88)+50];
C[(i*88)+31] += values[676] * B[(i*88)+50];
C[(i*88)+50] += values[677] * B[(i*88)+50];
C[(i*88)+51] += values[678] * B[(i*88)+50];
C[(i*88)+52] += values[679] * B[(i*88)+50];
C[(i*88)+78] += values[680] * B[(i*88)+50];
C[(i*88)+79] += values[681] * B[(i*88)+50];
C[(i*88)+80] += values[682] * B[(i*88)+50];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b51 = _mm256_broadcast_sd(&B[(i*88)+51]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b51 = _mm_loaddup_pd(&B[(i*88)+51]);
#endif
__m128d c51_0 = _mm_loadu_pd(&C[(i*88)+4]);
__m128d a51_0 = _mm_loadu_pd(&values[683]);
#if defined(__SSE3__) && defined(__AVX__)
c51_0 = _mm_add_pd(c51_0, _mm_mul_pd(a51_0, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c51_0 = _mm_add_pd(c51_0, _mm_mul_pd(a51_0, b51));
#endif
_mm_storeu_pd(&C[(i*88)+4], c51_0);
__m128d c51_2 = _mm_load_sd(&C[(i*88)+6]);
__m128d a51_2 = _mm_load_sd(&values[685]);
#if defined(__SSE3__) && defined(__AVX__)
c51_2 = _mm_add_sd(c51_2, _mm_mul_sd(a51_2, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c51_2 = _mm_add_sd(c51_2, _mm_mul_sd(a51_2, b51));
#endif
_mm_store_sd(&C[(i*88)+6], c51_2);
__m128d c51_3 = _mm_loadu_pd(&C[(i*88)+14]);
__m128d a51_3 = _mm_loadu_pd(&values[686]);
#if defined(__SSE3__) && defined(__AVX__)
c51_3 = _mm_add_pd(c51_3, _mm_mul_pd(a51_3, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c51_3 = _mm_add_pd(c51_3, _mm_mul_pd(a51_3, b51));
#endif
_mm_storeu_pd(&C[(i*88)+14], c51_3);
__m128d c51_5 = _mm_load_sd(&C[(i*88)+16]);
__m128d a51_5 = _mm_load_sd(&values[688]);
#if defined(__SSE3__) && defined(__AVX__)
c51_5 = _mm_add_sd(c51_5, _mm_mul_sd(a51_5, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c51_5 = _mm_add_sd(c51_5, _mm_mul_sd(a51_5, b51));
#endif
_mm_store_sd(&C[(i*88)+16], c51_5);
__m128d c51_6 = _mm_loadu_pd(&C[(i*88)+29]);
__m128d a51_6 = _mm_loadu_pd(&values[689]);
#if defined(__SSE3__) && defined(__AVX__)
c51_6 = _mm_add_pd(c51_6, _mm_mul_pd(a51_6, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c51_6 = _mm_add_pd(c51_6, _mm_mul_pd(a51_6, b51));
#endif
_mm_storeu_pd(&C[(i*88)+29], c51_6);
__m128d c51_8 = _mm_load_sd(&C[(i*88)+31]);
__m128d a51_8 = _mm_load_sd(&values[691]);
#if defined(__SSE3__) && defined(__AVX__)
c51_8 = _mm_add_sd(c51_8, _mm_mul_sd(a51_8, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c51_8 = _mm_add_sd(c51_8, _mm_mul_sd(a51_8, b51));
#endif
_mm_store_sd(&C[(i*88)+31], c51_8);
__m128d c51_9 = _mm_loadu_pd(&C[(i*88)+50]);
__m128d a51_9 = _mm_loadu_pd(&values[692]);
#if defined(__SSE3__) && defined(__AVX__)
c51_9 = _mm_add_pd(c51_9, _mm_mul_pd(a51_9, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c51_9 = _mm_add_pd(c51_9, _mm_mul_pd(a51_9, b51));
#endif
_mm_storeu_pd(&C[(i*88)+50], c51_9);
__m128d c51_11 = _mm_load_sd(&C[(i*88)+52]);
__m128d a51_11 = _mm_load_sd(&values[694]);
#if defined(__SSE3__) && defined(__AVX__)
c51_11 = _mm_add_sd(c51_11, _mm_mul_sd(a51_11, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c51_11 = _mm_add_sd(c51_11, _mm_mul_sd(a51_11, b51));
#endif
_mm_store_sd(&C[(i*88)+52], c51_11);
__m128d c51_12 = _mm_loadu_pd(&C[(i*88)+78]);
__m128d a51_12 = _mm_loadu_pd(&values[695]);
#if defined(__SSE3__) && defined(__AVX__)
c51_12 = _mm_add_pd(c51_12, _mm_mul_pd(a51_12, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c51_12 = _mm_add_pd(c51_12, _mm_mul_pd(a51_12, b51));
#endif
_mm_storeu_pd(&C[(i*88)+78], c51_12);
__m128d c51_14 = _mm_load_sd(&C[(i*88)+80]);
__m128d a51_14 = _mm_load_sd(&values[697]);
#if defined(__SSE3__) && defined(__AVX__)
c51_14 = _mm_add_sd(c51_14, _mm_mul_sd(a51_14, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c51_14 = _mm_add_sd(c51_14, _mm_mul_sd(a51_14, b51));
#endif
_mm_store_sd(&C[(i*88)+80], c51_14);
#else
C[(i*88)+4] += values[683] * B[(i*88)+51];
C[(i*88)+5] += values[684] * B[(i*88)+51];
C[(i*88)+6] += values[685] * B[(i*88)+51];
C[(i*88)+14] += values[686] * B[(i*88)+51];
C[(i*88)+15] += values[687] * B[(i*88)+51];
C[(i*88)+16] += values[688] * B[(i*88)+51];
C[(i*88)+29] += values[689] * B[(i*88)+51];
C[(i*88)+30] += values[690] * B[(i*88)+51];
C[(i*88)+31] += values[691] * B[(i*88)+51];
C[(i*88)+50] += values[692] * B[(i*88)+51];
C[(i*88)+51] += values[693] * B[(i*88)+51];
C[(i*88)+52] += values[694] * B[(i*88)+51];
C[(i*88)+78] += values[695] * B[(i*88)+51];
C[(i*88)+79] += values[696] * B[(i*88)+51];
C[(i*88)+80] += values[697] * B[(i*88)+51];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b52 = _mm256_broadcast_sd(&B[(i*88)+52]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b52 = _mm_loaddup_pd(&B[(i*88)+52]);
#endif
__m128d c52_0 = _mm_loadu_pd(&C[(i*88)+4]);
__m128d a52_0 = _mm_loadu_pd(&values[698]);
#if defined(__SSE3__) && defined(__AVX__)
c52_0 = _mm_add_pd(c52_0, _mm_mul_pd(a52_0, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c52_0 = _mm_add_pd(c52_0, _mm_mul_pd(a52_0, b52));
#endif
_mm_storeu_pd(&C[(i*88)+4], c52_0);
__m128d c52_2 = _mm_load_sd(&C[(i*88)+6]);
__m128d a52_2 = _mm_load_sd(&values[700]);
#if defined(__SSE3__) && defined(__AVX__)
c52_2 = _mm_add_sd(c52_2, _mm_mul_sd(a52_2, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c52_2 = _mm_add_sd(c52_2, _mm_mul_sd(a52_2, b52));
#endif
_mm_store_sd(&C[(i*88)+6], c52_2);
__m128d c52_3 = _mm_loadu_pd(&C[(i*88)+14]);
__m128d a52_3 = _mm_loadu_pd(&values[701]);
#if defined(__SSE3__) && defined(__AVX__)
c52_3 = _mm_add_pd(c52_3, _mm_mul_pd(a52_3, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c52_3 = _mm_add_pd(c52_3, _mm_mul_pd(a52_3, b52));
#endif
_mm_storeu_pd(&C[(i*88)+14], c52_3);
__m128d c52_5 = _mm_load_sd(&C[(i*88)+16]);
__m128d a52_5 = _mm_load_sd(&values[703]);
#if defined(__SSE3__) && defined(__AVX__)
c52_5 = _mm_add_sd(c52_5, _mm_mul_sd(a52_5, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c52_5 = _mm_add_sd(c52_5, _mm_mul_sd(a52_5, b52));
#endif
_mm_store_sd(&C[(i*88)+16], c52_5);
__m128d c52_6 = _mm_loadu_pd(&C[(i*88)+29]);
__m128d a52_6 = _mm_loadu_pd(&values[704]);
#if defined(__SSE3__) && defined(__AVX__)
c52_6 = _mm_add_pd(c52_6, _mm_mul_pd(a52_6, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c52_6 = _mm_add_pd(c52_6, _mm_mul_pd(a52_6, b52));
#endif
_mm_storeu_pd(&C[(i*88)+29], c52_6);
__m128d c52_8 = _mm_load_sd(&C[(i*88)+31]);
__m128d a52_8 = _mm_load_sd(&values[706]);
#if defined(__SSE3__) && defined(__AVX__)
c52_8 = _mm_add_sd(c52_8, _mm_mul_sd(a52_8, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c52_8 = _mm_add_sd(c52_8, _mm_mul_sd(a52_8, b52));
#endif
_mm_store_sd(&C[(i*88)+31], c52_8);
__m128d c52_9 = _mm_loadu_pd(&C[(i*88)+50]);
__m128d a52_9 = _mm_loadu_pd(&values[707]);
#if defined(__SSE3__) && defined(__AVX__)
c52_9 = _mm_add_pd(c52_9, _mm_mul_pd(a52_9, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c52_9 = _mm_add_pd(c52_9, _mm_mul_pd(a52_9, b52));
#endif
_mm_storeu_pd(&C[(i*88)+50], c52_9);
__m128d c52_11 = _mm_load_sd(&C[(i*88)+52]);
__m128d a52_11 = _mm_load_sd(&values[709]);
#if defined(__SSE3__) && defined(__AVX__)
c52_11 = _mm_add_sd(c52_11, _mm_mul_sd(a52_11, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c52_11 = _mm_add_sd(c52_11, _mm_mul_sd(a52_11, b52));
#endif
_mm_store_sd(&C[(i*88)+52], c52_11);
__m128d c52_12 = _mm_loadu_pd(&C[(i*88)+78]);
__m128d a52_12 = _mm_loadu_pd(&values[710]);
#if defined(__SSE3__) && defined(__AVX__)
c52_12 = _mm_add_pd(c52_12, _mm_mul_pd(a52_12, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c52_12 = _mm_add_pd(c52_12, _mm_mul_pd(a52_12, b52));
#endif
_mm_storeu_pd(&C[(i*88)+78], c52_12);
__m128d c52_14 = _mm_load_sd(&C[(i*88)+80]);
__m128d a52_14 = _mm_load_sd(&values[712]);
#if defined(__SSE3__) && defined(__AVX__)
c52_14 = _mm_add_sd(c52_14, _mm_mul_sd(a52_14, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c52_14 = _mm_add_sd(c52_14, _mm_mul_sd(a52_14, b52));
#endif
_mm_store_sd(&C[(i*88)+80], c52_14);
#else
C[(i*88)+4] += values[698] * B[(i*88)+52];
C[(i*88)+5] += values[699] * B[(i*88)+52];
C[(i*88)+6] += values[700] * B[(i*88)+52];
C[(i*88)+14] += values[701] * B[(i*88)+52];
C[(i*88)+15] += values[702] * B[(i*88)+52];
C[(i*88)+16] += values[703] * B[(i*88)+52];
C[(i*88)+29] += values[704] * B[(i*88)+52];
C[(i*88)+30] += values[705] * B[(i*88)+52];
C[(i*88)+31] += values[706] * B[(i*88)+52];
C[(i*88)+50] += values[707] * B[(i*88)+52];
C[(i*88)+51] += values[708] * B[(i*88)+52];
C[(i*88)+52] += values[709] * B[(i*88)+52];
C[(i*88)+78] += values[710] * B[(i*88)+52];
C[(i*88)+79] += values[711] * B[(i*88)+52];
C[(i*88)+80] += values[712] * B[(i*88)+52];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b53 = _mm256_broadcast_sd(&B[(i*88)+53]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b53 = _mm_loaddup_pd(&B[(i*88)+53]);
#endif
__m128d c53_0 = _mm_loadu_pd(&C[(i*88)+1]);
__m128d a53_0 = _mm_loadu_pd(&values[713]);
#if defined(__SSE3__) && defined(__AVX__)
c53_0 = _mm_add_pd(c53_0, _mm_mul_pd(a53_0, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c53_0 = _mm_add_pd(c53_0, _mm_mul_pd(a53_0, b53));
#endif
_mm_storeu_pd(&C[(i*88)+1], c53_0);
__m128d c53_2 = _mm_loadu_pd(&C[(i*88)+7]);
__m128d a53_2 = _mm_loadu_pd(&values[715]);
#if defined(__SSE3__) && defined(__AVX__)
c53_2 = _mm_add_pd(c53_2, _mm_mul_pd(a53_2, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c53_2 = _mm_add_pd(c53_2, _mm_mul_pd(a53_2, b53));
#endif
_mm_storeu_pd(&C[(i*88)+7], c53_2);
__m128d c53_4 = _mm_loadu_pd(&C[(i*88)+17]);
__m128d a53_4 = _mm_loadu_pd(&values[717]);
#if defined(__SSE3__) && defined(__AVX__)
c53_4 = _mm_add_pd(c53_4, _mm_mul_pd(a53_4, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c53_4 = _mm_add_pd(c53_4, _mm_mul_pd(a53_4, b53));
#endif
_mm_storeu_pd(&C[(i*88)+17], c53_4);
__m128d c53_6 = _mm_loadu_pd(&C[(i*88)+32]);
__m128d a53_6 = _mm_loadu_pd(&values[719]);
#if defined(__SSE3__) && defined(__AVX__)
c53_6 = _mm_add_pd(c53_6, _mm_mul_pd(a53_6, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c53_6 = _mm_add_pd(c53_6, _mm_mul_pd(a53_6, b53));
#endif
_mm_storeu_pd(&C[(i*88)+32], c53_6);
__m128d c53_8 = _mm_loadu_pd(&C[(i*88)+53]);
__m128d a53_8 = _mm_loadu_pd(&values[721]);
#if defined(__SSE3__) && defined(__AVX__)
c53_8 = _mm_add_pd(c53_8, _mm_mul_pd(a53_8, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c53_8 = _mm_add_pd(c53_8, _mm_mul_pd(a53_8, b53));
#endif
_mm_storeu_pd(&C[(i*88)+53], c53_8);
__m128d c53_10 = _mm_loadu_pd(&C[(i*88)+81]);
__m128d a53_10 = _mm_loadu_pd(&values[723]);
#if defined(__SSE3__) && defined(__AVX__)
c53_10 = _mm_add_pd(c53_10, _mm_mul_pd(a53_10, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c53_10 = _mm_add_pd(c53_10, _mm_mul_pd(a53_10, b53));
#endif
_mm_storeu_pd(&C[(i*88)+81], c53_10);
#else
C[(i*88)+1] += values[713] * B[(i*88)+53];
C[(i*88)+2] += values[714] * B[(i*88)+53];
C[(i*88)+7] += values[715] * B[(i*88)+53];
C[(i*88)+8] += values[716] * B[(i*88)+53];
C[(i*88)+17] += values[717] * B[(i*88)+53];
C[(i*88)+18] += values[718] * B[(i*88)+53];
C[(i*88)+32] += values[719] * B[(i*88)+53];
C[(i*88)+33] += values[720] * B[(i*88)+53];
C[(i*88)+53] += values[721] * B[(i*88)+53];
C[(i*88)+54] += values[722] * B[(i*88)+53];
C[(i*88)+81] += values[723] * B[(i*88)+53];
C[(i*88)+82] += values[724] * B[(i*88)+53];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b54 = _mm256_broadcast_sd(&B[(i*88)+54]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b54 = _mm_loaddup_pd(&B[(i*88)+54]);
#endif
__m128d c54_0 = _mm_loadu_pd(&C[(i*88)+1]);
__m128d a54_0 = _mm_loadu_pd(&values[725]);
#if defined(__SSE3__) && defined(__AVX__)
c54_0 = _mm_add_pd(c54_0, _mm_mul_pd(a54_0, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c54_0 = _mm_add_pd(c54_0, _mm_mul_pd(a54_0, b54));
#endif
_mm_storeu_pd(&C[(i*88)+1], c54_0);
__m128d c54_2 = _mm_loadu_pd(&C[(i*88)+7]);
__m128d a54_2 = _mm_loadu_pd(&values[727]);
#if defined(__SSE3__) && defined(__AVX__)
c54_2 = _mm_add_pd(c54_2, _mm_mul_pd(a54_2, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c54_2 = _mm_add_pd(c54_2, _mm_mul_pd(a54_2, b54));
#endif
_mm_storeu_pd(&C[(i*88)+7], c54_2);
__m128d c54_4 = _mm_loadu_pd(&C[(i*88)+17]);
__m128d a54_4 = _mm_loadu_pd(&values[729]);
#if defined(__SSE3__) && defined(__AVX__)
c54_4 = _mm_add_pd(c54_4, _mm_mul_pd(a54_4, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c54_4 = _mm_add_pd(c54_4, _mm_mul_pd(a54_4, b54));
#endif
_mm_storeu_pd(&C[(i*88)+17], c54_4);
__m128d c54_6 = _mm_loadu_pd(&C[(i*88)+32]);
__m128d a54_6 = _mm_loadu_pd(&values[731]);
#if defined(__SSE3__) && defined(__AVX__)
c54_6 = _mm_add_pd(c54_6, _mm_mul_pd(a54_6, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c54_6 = _mm_add_pd(c54_6, _mm_mul_pd(a54_6, b54));
#endif
_mm_storeu_pd(&C[(i*88)+32], c54_6);
__m128d c54_8 = _mm_loadu_pd(&C[(i*88)+53]);
__m128d a54_8 = _mm_loadu_pd(&values[733]);
#if defined(__SSE3__) && defined(__AVX__)
c54_8 = _mm_add_pd(c54_8, _mm_mul_pd(a54_8, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c54_8 = _mm_add_pd(c54_8, _mm_mul_pd(a54_8, b54));
#endif
_mm_storeu_pd(&C[(i*88)+53], c54_8);
__m128d c54_10 = _mm_loadu_pd(&C[(i*88)+81]);
__m128d a54_10 = _mm_loadu_pd(&values[735]);
#if defined(__SSE3__) && defined(__AVX__)
c54_10 = _mm_add_pd(c54_10, _mm_mul_pd(a54_10, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c54_10 = _mm_add_pd(c54_10, _mm_mul_pd(a54_10, b54));
#endif
_mm_storeu_pd(&C[(i*88)+81], c54_10);
#else
C[(i*88)+1] += values[725] * B[(i*88)+54];
C[(i*88)+2] += values[726] * B[(i*88)+54];
C[(i*88)+7] += values[727] * B[(i*88)+54];
C[(i*88)+8] += values[728] * B[(i*88)+54];
C[(i*88)+17] += values[729] * B[(i*88)+54];
C[(i*88)+18] += values[730] * B[(i*88)+54];
C[(i*88)+32] += values[731] * B[(i*88)+54];
C[(i*88)+33] += values[732] * B[(i*88)+54];
C[(i*88)+53] += values[733] * B[(i*88)+54];
C[(i*88)+54] += values[734] * B[(i*88)+54];
C[(i*88)+81] += values[735] * B[(i*88)+54];
C[(i*88)+82] += values[736] * B[(i*88)+54];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b55 = _mm256_broadcast_sd(&B[(i*88)+55]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b55 = _mm_loaddup_pd(&B[(i*88)+55]);
#endif
__m128d c55_0 = _mm_load_sd(&C[(i*88)+0]);
__m128d a55_0 = _mm_load_sd(&values[737]);
#if defined(__SSE3__) && defined(__AVX__)
c55_0 = _mm_add_sd(c55_0, _mm_mul_sd(a55_0, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_0 = _mm_add_sd(c55_0, _mm_mul_sd(a55_0, b55));
#endif
_mm_store_sd(&C[(i*88)+0], c55_0);
__m128d c55_1 = _mm_load_sd(&C[(i*88)+3]);
__m128d a55_1 = _mm_load_sd(&values[738]);
#if defined(__SSE3__) && defined(__AVX__)
c55_1 = _mm_add_sd(c55_1, _mm_mul_sd(a55_1, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_1 = _mm_add_sd(c55_1, _mm_mul_sd(a55_1, b55));
#endif
_mm_store_sd(&C[(i*88)+3], c55_1);
__m128d c55_2 = _mm_load_sd(&C[(i*88)+9]);
__m128d a55_2 = _mm_load_sd(&values[739]);
#if defined(__SSE3__) && defined(__AVX__)
c55_2 = _mm_add_sd(c55_2, _mm_mul_sd(a55_2, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_2 = _mm_add_sd(c55_2, _mm_mul_sd(a55_2, b55));
#endif
_mm_store_sd(&C[(i*88)+9], c55_2);
__m128d c55_3 = _mm_load_sd(&C[(i*88)+19]);
__m128d a55_3 = _mm_load_sd(&values[740]);
#if defined(__SSE3__) && defined(__AVX__)
c55_3 = _mm_add_sd(c55_3, _mm_mul_sd(a55_3, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_3 = _mm_add_sd(c55_3, _mm_mul_sd(a55_3, b55));
#endif
_mm_store_sd(&C[(i*88)+19], c55_3);
__m128d c55_4 = _mm_load_sd(&C[(i*88)+34]);
__m128d a55_4 = _mm_load_sd(&values[741]);
#if defined(__SSE3__) && defined(__AVX__)
c55_4 = _mm_add_sd(c55_4, _mm_mul_sd(a55_4, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_4 = _mm_add_sd(c55_4, _mm_mul_sd(a55_4, b55));
#endif
_mm_store_sd(&C[(i*88)+34], c55_4);
__m128d c55_5 = _mm_load_sd(&C[(i*88)+55]);
__m128d a55_5 = _mm_load_sd(&values[742]);
#if defined(__SSE3__) && defined(__AVX__)
c55_5 = _mm_add_sd(c55_5, _mm_mul_sd(a55_5, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_5 = _mm_add_sd(c55_5, _mm_mul_sd(a55_5, b55));
#endif
_mm_store_sd(&C[(i*88)+55], c55_5);
__m128d c55_6 = _mm_load_sd(&C[(i*88)+83]);
__m128d a55_6 = _mm_load_sd(&values[743]);
#if defined(__SSE3__) && defined(__AVX__)
c55_6 = _mm_add_sd(c55_6, _mm_mul_sd(a55_6, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_6 = _mm_add_sd(c55_6, _mm_mul_sd(a55_6, b55));
#endif
_mm_store_sd(&C[(i*88)+83], c55_6);
#else
C[(i*88)+0] += values[737] * B[(i*88)+55];
C[(i*88)+3] += values[738] * B[(i*88)+55];
C[(i*88)+9] += values[739] * B[(i*88)+55];
C[(i*88)+19] += values[740] * B[(i*88)+55];
C[(i*88)+34] += values[741] * B[(i*88)+55];
C[(i*88)+55] += values[742] * B[(i*88)+55];
C[(i*88)+83] += values[743] * B[(i*88)+55];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b56 = _mm256_broadcast_sd(&B[(i*88)+56]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b56 = _mm_loaddup_pd(&B[(i*88)+56]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c56_0 = _mm256_loadu_pd(&C[(i*88)+56]);
__m256d a56_0 = _mm256_loadu_pd(&values[744]);
c56_0 = _mm256_add_pd(c56_0, _mm256_mul_pd(a56_0, b56));
_mm256_storeu_pd(&C[(i*88)+56], c56_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c56_0 = _mm_loadu_pd(&C[(i*88)+56]);
__m128d a56_0 = _mm_loadu_pd(&values[744]);
c56_0 = _mm_add_pd(c56_0, _mm_mul_pd(a56_0, b56));
_mm_storeu_pd(&C[(i*88)+56], c56_0);
__m128d c56_2 = _mm_loadu_pd(&C[(i*88)+58]);
__m128d a56_2 = _mm_loadu_pd(&values[746]);
c56_2 = _mm_add_pd(c56_2, _mm_mul_pd(a56_2, b56));
_mm_storeu_pd(&C[(i*88)+58], c56_2);
#endif
__m128d c56_4 = _mm_loadu_pd(&C[(i*88)+60]);
__m128d a56_4 = _mm_loadu_pd(&values[748]);
#if defined(__SSE3__) && defined(__AVX__)
c56_4 = _mm_add_pd(c56_4, _mm_mul_pd(a56_4, _mm256_castpd256_pd128(b56)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c56_4 = _mm_add_pd(c56_4, _mm_mul_pd(a56_4, b56));
#endif
_mm_storeu_pd(&C[(i*88)+60], c56_4);
__m128d c56_6 = _mm_load_sd(&C[(i*88)+62]);
__m128d a56_6 = _mm_load_sd(&values[750]);
#if defined(__SSE3__) && defined(__AVX__)
c56_6 = _mm_add_sd(c56_6, _mm_mul_sd(a56_6, _mm256_castpd256_pd128(b56)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c56_6 = _mm_add_sd(c56_6, _mm_mul_sd(a56_6, b56));
#endif
_mm_store_sd(&C[(i*88)+62], c56_6);
#else
C[(i*88)+56] += values[744] * B[(i*88)+56];
C[(i*88)+57] += values[745] * B[(i*88)+56];
C[(i*88)+58] += values[746] * B[(i*88)+56];
C[(i*88)+59] += values[747] * B[(i*88)+56];
C[(i*88)+60] += values[748] * B[(i*88)+56];
C[(i*88)+61] += values[749] * B[(i*88)+56];
C[(i*88)+62] += values[750] * B[(i*88)+56];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b57 = _mm256_broadcast_sd(&B[(i*88)+57]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b57 = _mm_loaddup_pd(&B[(i*88)+57]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c57_0 = _mm256_loadu_pd(&C[(i*88)+56]);
__m256d a57_0 = _mm256_loadu_pd(&values[751]);
c57_0 = _mm256_add_pd(c57_0, _mm256_mul_pd(a57_0, b57));
_mm256_storeu_pd(&C[(i*88)+56], c57_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c57_0 = _mm_loadu_pd(&C[(i*88)+56]);
__m128d a57_0 = _mm_loadu_pd(&values[751]);
c57_0 = _mm_add_pd(c57_0, _mm_mul_pd(a57_0, b57));
_mm_storeu_pd(&C[(i*88)+56], c57_0);
__m128d c57_2 = _mm_loadu_pd(&C[(i*88)+58]);
__m128d a57_2 = _mm_loadu_pd(&values[753]);
c57_2 = _mm_add_pd(c57_2, _mm_mul_pd(a57_2, b57));
_mm_storeu_pd(&C[(i*88)+58], c57_2);
#endif
__m128d c57_4 = _mm_loadu_pd(&C[(i*88)+60]);
__m128d a57_4 = _mm_loadu_pd(&values[755]);
#if defined(__SSE3__) && defined(__AVX__)
c57_4 = _mm_add_pd(c57_4, _mm_mul_pd(a57_4, _mm256_castpd256_pd128(b57)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c57_4 = _mm_add_pd(c57_4, _mm_mul_pd(a57_4, b57));
#endif
_mm_storeu_pd(&C[(i*88)+60], c57_4);
__m128d c57_6 = _mm_load_sd(&C[(i*88)+62]);
__m128d a57_6 = _mm_load_sd(&values[757]);
#if defined(__SSE3__) && defined(__AVX__)
c57_6 = _mm_add_sd(c57_6, _mm_mul_sd(a57_6, _mm256_castpd256_pd128(b57)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c57_6 = _mm_add_sd(c57_6, _mm_mul_sd(a57_6, b57));
#endif
_mm_store_sd(&C[(i*88)+62], c57_6);
#else
C[(i*88)+56] += values[751] * B[(i*88)+57];
C[(i*88)+57] += values[752] * B[(i*88)+57];
C[(i*88)+58] += values[753] * B[(i*88)+57];
C[(i*88)+59] += values[754] * B[(i*88)+57];
C[(i*88)+60] += values[755] * B[(i*88)+57];
C[(i*88)+61] += values[756] * B[(i*88)+57];
C[(i*88)+62] += values[757] * B[(i*88)+57];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b58 = _mm256_broadcast_sd(&B[(i*88)+58]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b58 = _mm_loaddup_pd(&B[(i*88)+58]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c58_0 = _mm256_loadu_pd(&C[(i*88)+56]);
__m256d a58_0 = _mm256_loadu_pd(&values[758]);
c58_0 = _mm256_add_pd(c58_0, _mm256_mul_pd(a58_0, b58));
_mm256_storeu_pd(&C[(i*88)+56], c58_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c58_0 = _mm_loadu_pd(&C[(i*88)+56]);
__m128d a58_0 = _mm_loadu_pd(&values[758]);
c58_0 = _mm_add_pd(c58_0, _mm_mul_pd(a58_0, b58));
_mm_storeu_pd(&C[(i*88)+56], c58_0);
__m128d c58_2 = _mm_loadu_pd(&C[(i*88)+58]);
__m128d a58_2 = _mm_loadu_pd(&values[760]);
c58_2 = _mm_add_pd(c58_2, _mm_mul_pd(a58_2, b58));
_mm_storeu_pd(&C[(i*88)+58], c58_2);
#endif
__m128d c58_4 = _mm_loadu_pd(&C[(i*88)+60]);
__m128d a58_4 = _mm_loadu_pd(&values[762]);
#if defined(__SSE3__) && defined(__AVX__)
c58_4 = _mm_add_pd(c58_4, _mm_mul_pd(a58_4, _mm256_castpd256_pd128(b58)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c58_4 = _mm_add_pd(c58_4, _mm_mul_pd(a58_4, b58));
#endif
_mm_storeu_pd(&C[(i*88)+60], c58_4);
__m128d c58_6 = _mm_load_sd(&C[(i*88)+62]);
__m128d a58_6 = _mm_load_sd(&values[764]);
#if defined(__SSE3__) && defined(__AVX__)
c58_6 = _mm_add_sd(c58_6, _mm_mul_sd(a58_6, _mm256_castpd256_pd128(b58)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c58_6 = _mm_add_sd(c58_6, _mm_mul_sd(a58_6, b58));
#endif
_mm_store_sd(&C[(i*88)+62], c58_6);
#else
C[(i*88)+56] += values[758] * B[(i*88)+58];
C[(i*88)+57] += values[759] * B[(i*88)+58];
C[(i*88)+58] += values[760] * B[(i*88)+58];
C[(i*88)+59] += values[761] * B[(i*88)+58];
C[(i*88)+60] += values[762] * B[(i*88)+58];
C[(i*88)+61] += values[763] * B[(i*88)+58];
C[(i*88)+62] += values[764] * B[(i*88)+58];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b59 = _mm256_broadcast_sd(&B[(i*88)+59]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b59 = _mm_loaddup_pd(&B[(i*88)+59]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c59_0 = _mm256_loadu_pd(&C[(i*88)+56]);
__m256d a59_0 = _mm256_loadu_pd(&values[765]);
c59_0 = _mm256_add_pd(c59_0, _mm256_mul_pd(a59_0, b59));
_mm256_storeu_pd(&C[(i*88)+56], c59_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c59_0 = _mm_loadu_pd(&C[(i*88)+56]);
__m128d a59_0 = _mm_loadu_pd(&values[765]);
c59_0 = _mm_add_pd(c59_0, _mm_mul_pd(a59_0, b59));
_mm_storeu_pd(&C[(i*88)+56], c59_0);
__m128d c59_2 = _mm_loadu_pd(&C[(i*88)+58]);
__m128d a59_2 = _mm_loadu_pd(&values[767]);
c59_2 = _mm_add_pd(c59_2, _mm_mul_pd(a59_2, b59));
_mm_storeu_pd(&C[(i*88)+58], c59_2);
#endif
__m128d c59_4 = _mm_loadu_pd(&C[(i*88)+60]);
__m128d a59_4 = _mm_loadu_pd(&values[769]);
#if defined(__SSE3__) && defined(__AVX__)
c59_4 = _mm_add_pd(c59_4, _mm_mul_pd(a59_4, _mm256_castpd256_pd128(b59)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c59_4 = _mm_add_pd(c59_4, _mm_mul_pd(a59_4, b59));
#endif
_mm_storeu_pd(&C[(i*88)+60], c59_4);
__m128d c59_6 = _mm_load_sd(&C[(i*88)+62]);
__m128d a59_6 = _mm_load_sd(&values[771]);
#if defined(__SSE3__) && defined(__AVX__)
c59_6 = _mm_add_sd(c59_6, _mm_mul_sd(a59_6, _mm256_castpd256_pd128(b59)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c59_6 = _mm_add_sd(c59_6, _mm_mul_sd(a59_6, b59));
#endif
_mm_store_sd(&C[(i*88)+62], c59_6);
#else
C[(i*88)+56] += values[765] * B[(i*88)+59];
C[(i*88)+57] += values[766] * B[(i*88)+59];
C[(i*88)+58] += values[767] * B[(i*88)+59];
C[(i*88)+59] += values[768] * B[(i*88)+59];
C[(i*88)+60] += values[769] * B[(i*88)+59];
C[(i*88)+61] += values[770] * B[(i*88)+59];
C[(i*88)+62] += values[771] * B[(i*88)+59];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b60 = _mm256_broadcast_sd(&B[(i*88)+60]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b60 = _mm_loaddup_pd(&B[(i*88)+60]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c60_0 = _mm256_loadu_pd(&C[(i*88)+56]);
__m256d a60_0 = _mm256_loadu_pd(&values[772]);
c60_0 = _mm256_add_pd(c60_0, _mm256_mul_pd(a60_0, b60));
_mm256_storeu_pd(&C[(i*88)+56], c60_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c60_0 = _mm_loadu_pd(&C[(i*88)+56]);
__m128d a60_0 = _mm_loadu_pd(&values[772]);
c60_0 = _mm_add_pd(c60_0, _mm_mul_pd(a60_0, b60));
_mm_storeu_pd(&C[(i*88)+56], c60_0);
__m128d c60_2 = _mm_loadu_pd(&C[(i*88)+58]);
__m128d a60_2 = _mm_loadu_pd(&values[774]);
c60_2 = _mm_add_pd(c60_2, _mm_mul_pd(a60_2, b60));
_mm_storeu_pd(&C[(i*88)+58], c60_2);
#endif
__m128d c60_4 = _mm_loadu_pd(&C[(i*88)+60]);
__m128d a60_4 = _mm_loadu_pd(&values[776]);
#if defined(__SSE3__) && defined(__AVX__)
c60_4 = _mm_add_pd(c60_4, _mm_mul_pd(a60_4, _mm256_castpd256_pd128(b60)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c60_4 = _mm_add_pd(c60_4, _mm_mul_pd(a60_4, b60));
#endif
_mm_storeu_pd(&C[(i*88)+60], c60_4);
__m128d c60_6 = _mm_load_sd(&C[(i*88)+62]);
__m128d a60_6 = _mm_load_sd(&values[778]);
#if defined(__SSE3__) && defined(__AVX__)
c60_6 = _mm_add_sd(c60_6, _mm_mul_sd(a60_6, _mm256_castpd256_pd128(b60)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c60_6 = _mm_add_sd(c60_6, _mm_mul_sd(a60_6, b60));
#endif
_mm_store_sd(&C[(i*88)+62], c60_6);
#else
C[(i*88)+56] += values[772] * B[(i*88)+60];
C[(i*88)+57] += values[773] * B[(i*88)+60];
C[(i*88)+58] += values[774] * B[(i*88)+60];
C[(i*88)+59] += values[775] * B[(i*88)+60];
C[(i*88)+60] += values[776] * B[(i*88)+60];
C[(i*88)+61] += values[777] * B[(i*88)+60];
C[(i*88)+62] += values[778] * B[(i*88)+60];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b61 = _mm256_broadcast_sd(&B[(i*88)+61]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b61 = _mm_loaddup_pd(&B[(i*88)+61]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c61_0 = _mm256_loadu_pd(&C[(i*88)+56]);
__m256d a61_0 = _mm256_loadu_pd(&values[779]);
c61_0 = _mm256_add_pd(c61_0, _mm256_mul_pd(a61_0, b61));
_mm256_storeu_pd(&C[(i*88)+56], c61_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c61_0 = _mm_loadu_pd(&C[(i*88)+56]);
__m128d a61_0 = _mm_loadu_pd(&values[779]);
c61_0 = _mm_add_pd(c61_0, _mm_mul_pd(a61_0, b61));
_mm_storeu_pd(&C[(i*88)+56], c61_0);
__m128d c61_2 = _mm_loadu_pd(&C[(i*88)+58]);
__m128d a61_2 = _mm_loadu_pd(&values[781]);
c61_2 = _mm_add_pd(c61_2, _mm_mul_pd(a61_2, b61));
_mm_storeu_pd(&C[(i*88)+58], c61_2);
#endif
__m128d c61_4 = _mm_loadu_pd(&C[(i*88)+60]);
__m128d a61_4 = _mm_loadu_pd(&values[783]);
#if defined(__SSE3__) && defined(__AVX__)
c61_4 = _mm_add_pd(c61_4, _mm_mul_pd(a61_4, _mm256_castpd256_pd128(b61)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c61_4 = _mm_add_pd(c61_4, _mm_mul_pd(a61_4, b61));
#endif
_mm_storeu_pd(&C[(i*88)+60], c61_4);
__m128d c61_6 = _mm_load_sd(&C[(i*88)+62]);
__m128d a61_6 = _mm_load_sd(&values[785]);
#if defined(__SSE3__) && defined(__AVX__)
c61_6 = _mm_add_sd(c61_6, _mm_mul_sd(a61_6, _mm256_castpd256_pd128(b61)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c61_6 = _mm_add_sd(c61_6, _mm_mul_sd(a61_6, b61));
#endif
_mm_store_sd(&C[(i*88)+62], c61_6);
#else
C[(i*88)+56] += values[779] * B[(i*88)+61];
C[(i*88)+57] += values[780] * B[(i*88)+61];
C[(i*88)+58] += values[781] * B[(i*88)+61];
C[(i*88)+59] += values[782] * B[(i*88)+61];
C[(i*88)+60] += values[783] * B[(i*88)+61];
C[(i*88)+61] += values[784] * B[(i*88)+61];
C[(i*88)+62] += values[785] * B[(i*88)+61];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b62 = _mm256_broadcast_sd(&B[(i*88)+62]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b62 = _mm_loaddup_pd(&B[(i*88)+62]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c62_0 = _mm256_loadu_pd(&C[(i*88)+56]);
__m256d a62_0 = _mm256_loadu_pd(&values[786]);
c62_0 = _mm256_add_pd(c62_0, _mm256_mul_pd(a62_0, b62));
_mm256_storeu_pd(&C[(i*88)+56], c62_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c62_0 = _mm_loadu_pd(&C[(i*88)+56]);
__m128d a62_0 = _mm_loadu_pd(&values[786]);
c62_0 = _mm_add_pd(c62_0, _mm_mul_pd(a62_0, b62));
_mm_storeu_pd(&C[(i*88)+56], c62_0);
__m128d c62_2 = _mm_loadu_pd(&C[(i*88)+58]);
__m128d a62_2 = _mm_loadu_pd(&values[788]);
c62_2 = _mm_add_pd(c62_2, _mm_mul_pd(a62_2, b62));
_mm_storeu_pd(&C[(i*88)+58], c62_2);
#endif
__m128d c62_4 = _mm_loadu_pd(&C[(i*88)+60]);
__m128d a62_4 = _mm_loadu_pd(&values[790]);
#if defined(__SSE3__) && defined(__AVX__)
c62_4 = _mm_add_pd(c62_4, _mm_mul_pd(a62_4, _mm256_castpd256_pd128(b62)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c62_4 = _mm_add_pd(c62_4, _mm_mul_pd(a62_4, b62));
#endif
_mm_storeu_pd(&C[(i*88)+60], c62_4);
__m128d c62_6 = _mm_load_sd(&C[(i*88)+62]);
__m128d a62_6 = _mm_load_sd(&values[792]);
#if defined(__SSE3__) && defined(__AVX__)
c62_6 = _mm_add_sd(c62_6, _mm_mul_sd(a62_6, _mm256_castpd256_pd128(b62)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c62_6 = _mm_add_sd(c62_6, _mm_mul_sd(a62_6, b62));
#endif
_mm_store_sd(&C[(i*88)+62], c62_6);
#else
C[(i*88)+56] += values[786] * B[(i*88)+62];
C[(i*88)+57] += values[787] * B[(i*88)+62];
C[(i*88)+58] += values[788] * B[(i*88)+62];
C[(i*88)+59] += values[789] * B[(i*88)+62];
C[(i*88)+60] += values[790] * B[(i*88)+62];
C[(i*88)+61] += values[791] * B[(i*88)+62];
C[(i*88)+62] += values[792] * B[(i*88)+62];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b63 = _mm256_broadcast_sd(&B[(i*88)+63]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b63 = _mm_loaddup_pd(&B[(i*88)+63]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c63_0 = _mm256_loadu_pd(&C[(i*88)+35]);
__m256d a63_0 = _mm256_loadu_pd(&values[793]);
c63_0 = _mm256_add_pd(c63_0, _mm256_mul_pd(a63_0, b63));
_mm256_storeu_pd(&C[(i*88)+35], c63_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c63_0 = _mm_loadu_pd(&C[(i*88)+35]);
__m128d a63_0 = _mm_loadu_pd(&values[793]);
c63_0 = _mm_add_pd(c63_0, _mm_mul_pd(a63_0, b63));
_mm_storeu_pd(&C[(i*88)+35], c63_0);
__m128d c63_2 = _mm_loadu_pd(&C[(i*88)+37]);
__m128d a63_2 = _mm_loadu_pd(&values[795]);
c63_2 = _mm_add_pd(c63_2, _mm_mul_pd(a63_2, b63));
_mm_storeu_pd(&C[(i*88)+37], c63_2);
#endif
__m128d c63_4 = _mm_loadu_pd(&C[(i*88)+39]);
__m128d a63_4 = _mm_loadu_pd(&values[797]);
#if defined(__SSE3__) && defined(__AVX__)
c63_4 = _mm_add_pd(c63_4, _mm_mul_pd(a63_4, _mm256_castpd256_pd128(b63)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c63_4 = _mm_add_pd(c63_4, _mm_mul_pd(a63_4, b63));
#endif
_mm_storeu_pd(&C[(i*88)+39], c63_4);
#if defined(__SSE3__) && defined(__AVX__)
__m256d c63_6 = _mm256_loadu_pd(&C[(i*88)+63]);
__m256d a63_6 = _mm256_loadu_pd(&values[799]);
c63_6 = _mm256_add_pd(c63_6, _mm256_mul_pd(a63_6, b63));
_mm256_storeu_pd(&C[(i*88)+63], c63_6);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c63_6 = _mm_loadu_pd(&C[(i*88)+63]);
__m128d a63_6 = _mm_loadu_pd(&values[799]);
c63_6 = _mm_add_pd(c63_6, _mm_mul_pd(a63_6, b63));
_mm_storeu_pd(&C[(i*88)+63], c63_6);
__m128d c63_8 = _mm_loadu_pd(&C[(i*88)+65]);
__m128d a63_8 = _mm_loadu_pd(&values[801]);
c63_8 = _mm_add_pd(c63_8, _mm_mul_pd(a63_8, b63));
_mm_storeu_pd(&C[(i*88)+65], c63_8);
#endif
__m128d c63_10 = _mm_loadu_pd(&C[(i*88)+67]);
__m128d a63_10 = _mm_loadu_pd(&values[803]);
#if defined(__SSE3__) && defined(__AVX__)
c63_10 = _mm_add_pd(c63_10, _mm_mul_pd(a63_10, _mm256_castpd256_pd128(b63)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c63_10 = _mm_add_pd(c63_10, _mm_mul_pd(a63_10, b63));
#endif
_mm_storeu_pd(&C[(i*88)+67], c63_10);
#else
C[(i*88)+35] += values[793] * B[(i*88)+63];
C[(i*88)+36] += values[794] * B[(i*88)+63];
C[(i*88)+37] += values[795] * B[(i*88)+63];
C[(i*88)+38] += values[796] * B[(i*88)+63];
C[(i*88)+39] += values[797] * B[(i*88)+63];
C[(i*88)+40] += values[798] * B[(i*88)+63];
C[(i*88)+63] += values[799] * B[(i*88)+63];
C[(i*88)+64] += values[800] * B[(i*88)+63];
C[(i*88)+65] += values[801] * B[(i*88)+63];
C[(i*88)+66] += values[802] * B[(i*88)+63];
C[(i*88)+67] += values[803] * B[(i*88)+63];
C[(i*88)+68] += values[804] * B[(i*88)+63];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b64 = _mm256_broadcast_sd(&B[(i*88)+64]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b64 = _mm_loaddup_pd(&B[(i*88)+64]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c64_0 = _mm256_loadu_pd(&C[(i*88)+35]);
__m256d a64_0 = _mm256_loadu_pd(&values[805]);
c64_0 = _mm256_add_pd(c64_0, _mm256_mul_pd(a64_0, b64));
_mm256_storeu_pd(&C[(i*88)+35], c64_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c64_0 = _mm_loadu_pd(&C[(i*88)+35]);
__m128d a64_0 = _mm_loadu_pd(&values[805]);
c64_0 = _mm_add_pd(c64_0, _mm_mul_pd(a64_0, b64));
_mm_storeu_pd(&C[(i*88)+35], c64_0);
__m128d c64_2 = _mm_loadu_pd(&C[(i*88)+37]);
__m128d a64_2 = _mm_loadu_pd(&values[807]);
c64_2 = _mm_add_pd(c64_2, _mm_mul_pd(a64_2, b64));
_mm_storeu_pd(&C[(i*88)+37], c64_2);
#endif
__m128d c64_4 = _mm_loadu_pd(&C[(i*88)+39]);
__m128d a64_4 = _mm_loadu_pd(&values[809]);
#if defined(__SSE3__) && defined(__AVX__)
c64_4 = _mm_add_pd(c64_4, _mm_mul_pd(a64_4, _mm256_castpd256_pd128(b64)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c64_4 = _mm_add_pd(c64_4, _mm_mul_pd(a64_4, b64));
#endif
_mm_storeu_pd(&C[(i*88)+39], c64_4);
#if defined(__SSE3__) && defined(__AVX__)
__m256d c64_6 = _mm256_loadu_pd(&C[(i*88)+63]);
__m256d a64_6 = _mm256_loadu_pd(&values[811]);
c64_6 = _mm256_add_pd(c64_6, _mm256_mul_pd(a64_6, b64));
_mm256_storeu_pd(&C[(i*88)+63], c64_6);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c64_6 = _mm_loadu_pd(&C[(i*88)+63]);
__m128d a64_6 = _mm_loadu_pd(&values[811]);
c64_6 = _mm_add_pd(c64_6, _mm_mul_pd(a64_6, b64));
_mm_storeu_pd(&C[(i*88)+63], c64_6);
__m128d c64_8 = _mm_loadu_pd(&C[(i*88)+65]);
__m128d a64_8 = _mm_loadu_pd(&values[813]);
c64_8 = _mm_add_pd(c64_8, _mm_mul_pd(a64_8, b64));
_mm_storeu_pd(&C[(i*88)+65], c64_8);
#endif
__m128d c64_10 = _mm_loadu_pd(&C[(i*88)+67]);
__m128d a64_10 = _mm_loadu_pd(&values[815]);
#if defined(__SSE3__) && defined(__AVX__)
c64_10 = _mm_add_pd(c64_10, _mm_mul_pd(a64_10, _mm256_castpd256_pd128(b64)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c64_10 = _mm_add_pd(c64_10, _mm_mul_pd(a64_10, b64));
#endif
_mm_storeu_pd(&C[(i*88)+67], c64_10);
#else
C[(i*88)+35] += values[805] * B[(i*88)+64];
C[(i*88)+36] += values[806] * B[(i*88)+64];
C[(i*88)+37] += values[807] * B[(i*88)+64];
C[(i*88)+38] += values[808] * B[(i*88)+64];
C[(i*88)+39] += values[809] * B[(i*88)+64];
C[(i*88)+40] += values[810] * B[(i*88)+64];
C[(i*88)+63] += values[811] * B[(i*88)+64];
C[(i*88)+64] += values[812] * B[(i*88)+64];
C[(i*88)+65] += values[813] * B[(i*88)+64];
C[(i*88)+66] += values[814] * B[(i*88)+64];
C[(i*88)+67] += values[815] * B[(i*88)+64];
C[(i*88)+68] += values[816] * B[(i*88)+64];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b65 = _mm256_broadcast_sd(&B[(i*88)+65]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b65 = _mm_loaddup_pd(&B[(i*88)+65]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c65_0 = _mm256_loadu_pd(&C[(i*88)+35]);
__m256d a65_0 = _mm256_loadu_pd(&values[817]);
c65_0 = _mm256_add_pd(c65_0, _mm256_mul_pd(a65_0, b65));
_mm256_storeu_pd(&C[(i*88)+35], c65_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c65_0 = _mm_loadu_pd(&C[(i*88)+35]);
__m128d a65_0 = _mm_loadu_pd(&values[817]);
c65_0 = _mm_add_pd(c65_0, _mm_mul_pd(a65_0, b65));
_mm_storeu_pd(&C[(i*88)+35], c65_0);
__m128d c65_2 = _mm_loadu_pd(&C[(i*88)+37]);
__m128d a65_2 = _mm_loadu_pd(&values[819]);
c65_2 = _mm_add_pd(c65_2, _mm_mul_pd(a65_2, b65));
_mm_storeu_pd(&C[(i*88)+37], c65_2);
#endif
__m128d c65_4 = _mm_loadu_pd(&C[(i*88)+39]);
__m128d a65_4 = _mm_loadu_pd(&values[821]);
#if defined(__SSE3__) && defined(__AVX__)
c65_4 = _mm_add_pd(c65_4, _mm_mul_pd(a65_4, _mm256_castpd256_pd128(b65)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c65_4 = _mm_add_pd(c65_4, _mm_mul_pd(a65_4, b65));
#endif
_mm_storeu_pd(&C[(i*88)+39], c65_4);
#if defined(__SSE3__) && defined(__AVX__)
__m256d c65_6 = _mm256_loadu_pd(&C[(i*88)+63]);
__m256d a65_6 = _mm256_loadu_pd(&values[823]);
c65_6 = _mm256_add_pd(c65_6, _mm256_mul_pd(a65_6, b65));
_mm256_storeu_pd(&C[(i*88)+63], c65_6);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c65_6 = _mm_loadu_pd(&C[(i*88)+63]);
__m128d a65_6 = _mm_loadu_pd(&values[823]);
c65_6 = _mm_add_pd(c65_6, _mm_mul_pd(a65_6, b65));
_mm_storeu_pd(&C[(i*88)+63], c65_6);
__m128d c65_8 = _mm_loadu_pd(&C[(i*88)+65]);
__m128d a65_8 = _mm_loadu_pd(&values[825]);
c65_8 = _mm_add_pd(c65_8, _mm_mul_pd(a65_8, b65));
_mm_storeu_pd(&C[(i*88)+65], c65_8);
#endif
__m128d c65_10 = _mm_loadu_pd(&C[(i*88)+67]);
__m128d a65_10 = _mm_loadu_pd(&values[827]);
#if defined(__SSE3__) && defined(__AVX__)
c65_10 = _mm_add_pd(c65_10, _mm_mul_pd(a65_10, _mm256_castpd256_pd128(b65)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c65_10 = _mm_add_pd(c65_10, _mm_mul_pd(a65_10, b65));
#endif
_mm_storeu_pd(&C[(i*88)+67], c65_10);
#else
C[(i*88)+35] += values[817] * B[(i*88)+65];
C[(i*88)+36] += values[818] * B[(i*88)+65];
C[(i*88)+37] += values[819] * B[(i*88)+65];
C[(i*88)+38] += values[820] * B[(i*88)+65];
C[(i*88)+39] += values[821] * B[(i*88)+65];
C[(i*88)+40] += values[822] * B[(i*88)+65];
C[(i*88)+63] += values[823] * B[(i*88)+65];
C[(i*88)+64] += values[824] * B[(i*88)+65];
C[(i*88)+65] += values[825] * B[(i*88)+65];
C[(i*88)+66] += values[826] * B[(i*88)+65];
C[(i*88)+67] += values[827] * B[(i*88)+65];
C[(i*88)+68] += values[828] * B[(i*88)+65];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b66 = _mm256_broadcast_sd(&B[(i*88)+66]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b66 = _mm_loaddup_pd(&B[(i*88)+66]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c66_0 = _mm256_loadu_pd(&C[(i*88)+35]);
__m256d a66_0 = _mm256_loadu_pd(&values[829]);
c66_0 = _mm256_add_pd(c66_0, _mm256_mul_pd(a66_0, b66));
_mm256_storeu_pd(&C[(i*88)+35], c66_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c66_0 = _mm_loadu_pd(&C[(i*88)+35]);
__m128d a66_0 = _mm_loadu_pd(&values[829]);
c66_0 = _mm_add_pd(c66_0, _mm_mul_pd(a66_0, b66));
_mm_storeu_pd(&C[(i*88)+35], c66_0);
__m128d c66_2 = _mm_loadu_pd(&C[(i*88)+37]);
__m128d a66_2 = _mm_loadu_pd(&values[831]);
c66_2 = _mm_add_pd(c66_2, _mm_mul_pd(a66_2, b66));
_mm_storeu_pd(&C[(i*88)+37], c66_2);
#endif
__m128d c66_4 = _mm_loadu_pd(&C[(i*88)+39]);
__m128d a66_4 = _mm_loadu_pd(&values[833]);
#if defined(__SSE3__) && defined(__AVX__)
c66_4 = _mm_add_pd(c66_4, _mm_mul_pd(a66_4, _mm256_castpd256_pd128(b66)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c66_4 = _mm_add_pd(c66_4, _mm_mul_pd(a66_4, b66));
#endif
_mm_storeu_pd(&C[(i*88)+39], c66_4);
#if defined(__SSE3__) && defined(__AVX__)
__m256d c66_6 = _mm256_loadu_pd(&C[(i*88)+63]);
__m256d a66_6 = _mm256_loadu_pd(&values[835]);
c66_6 = _mm256_add_pd(c66_6, _mm256_mul_pd(a66_6, b66));
_mm256_storeu_pd(&C[(i*88)+63], c66_6);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c66_6 = _mm_loadu_pd(&C[(i*88)+63]);
__m128d a66_6 = _mm_loadu_pd(&values[835]);
c66_6 = _mm_add_pd(c66_6, _mm_mul_pd(a66_6, b66));
_mm_storeu_pd(&C[(i*88)+63], c66_6);
__m128d c66_8 = _mm_loadu_pd(&C[(i*88)+65]);
__m128d a66_8 = _mm_loadu_pd(&values[837]);
c66_8 = _mm_add_pd(c66_8, _mm_mul_pd(a66_8, b66));
_mm_storeu_pd(&C[(i*88)+65], c66_8);
#endif
__m128d c66_10 = _mm_loadu_pd(&C[(i*88)+67]);
__m128d a66_10 = _mm_loadu_pd(&values[839]);
#if defined(__SSE3__) && defined(__AVX__)
c66_10 = _mm_add_pd(c66_10, _mm_mul_pd(a66_10, _mm256_castpd256_pd128(b66)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c66_10 = _mm_add_pd(c66_10, _mm_mul_pd(a66_10, b66));
#endif
_mm_storeu_pd(&C[(i*88)+67], c66_10);
#else
C[(i*88)+35] += values[829] * B[(i*88)+66];
C[(i*88)+36] += values[830] * B[(i*88)+66];
C[(i*88)+37] += values[831] * B[(i*88)+66];
C[(i*88)+38] += values[832] * B[(i*88)+66];
C[(i*88)+39] += values[833] * B[(i*88)+66];
C[(i*88)+40] += values[834] * B[(i*88)+66];
C[(i*88)+63] += values[835] * B[(i*88)+66];
C[(i*88)+64] += values[836] * B[(i*88)+66];
C[(i*88)+65] += values[837] * B[(i*88)+66];
C[(i*88)+66] += values[838] * B[(i*88)+66];
C[(i*88)+67] += values[839] * B[(i*88)+66];
C[(i*88)+68] += values[840] * B[(i*88)+66];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b67 = _mm256_broadcast_sd(&B[(i*88)+67]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b67 = _mm_loaddup_pd(&B[(i*88)+67]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c67_0 = _mm256_loadu_pd(&C[(i*88)+35]);
__m256d a67_0 = _mm256_loadu_pd(&values[841]);
c67_0 = _mm256_add_pd(c67_0, _mm256_mul_pd(a67_0, b67));
_mm256_storeu_pd(&C[(i*88)+35], c67_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c67_0 = _mm_loadu_pd(&C[(i*88)+35]);
__m128d a67_0 = _mm_loadu_pd(&values[841]);
c67_0 = _mm_add_pd(c67_0, _mm_mul_pd(a67_0, b67));
_mm_storeu_pd(&C[(i*88)+35], c67_0);
__m128d c67_2 = _mm_loadu_pd(&C[(i*88)+37]);
__m128d a67_2 = _mm_loadu_pd(&values[843]);
c67_2 = _mm_add_pd(c67_2, _mm_mul_pd(a67_2, b67));
_mm_storeu_pd(&C[(i*88)+37], c67_2);
#endif
__m128d c67_4 = _mm_loadu_pd(&C[(i*88)+39]);
__m128d a67_4 = _mm_loadu_pd(&values[845]);
#if defined(__SSE3__) && defined(__AVX__)
c67_4 = _mm_add_pd(c67_4, _mm_mul_pd(a67_4, _mm256_castpd256_pd128(b67)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c67_4 = _mm_add_pd(c67_4, _mm_mul_pd(a67_4, b67));
#endif
_mm_storeu_pd(&C[(i*88)+39], c67_4);
#if defined(__SSE3__) && defined(__AVX__)
__m256d c67_6 = _mm256_loadu_pd(&C[(i*88)+63]);
__m256d a67_6 = _mm256_loadu_pd(&values[847]);
c67_6 = _mm256_add_pd(c67_6, _mm256_mul_pd(a67_6, b67));
_mm256_storeu_pd(&C[(i*88)+63], c67_6);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c67_6 = _mm_loadu_pd(&C[(i*88)+63]);
__m128d a67_6 = _mm_loadu_pd(&values[847]);
c67_6 = _mm_add_pd(c67_6, _mm_mul_pd(a67_6, b67));
_mm_storeu_pd(&C[(i*88)+63], c67_6);
__m128d c67_8 = _mm_loadu_pd(&C[(i*88)+65]);
__m128d a67_8 = _mm_loadu_pd(&values[849]);
c67_8 = _mm_add_pd(c67_8, _mm_mul_pd(a67_8, b67));
_mm_storeu_pd(&C[(i*88)+65], c67_8);
#endif
__m128d c67_10 = _mm_loadu_pd(&C[(i*88)+67]);
__m128d a67_10 = _mm_loadu_pd(&values[851]);
#if defined(__SSE3__) && defined(__AVX__)
c67_10 = _mm_add_pd(c67_10, _mm_mul_pd(a67_10, _mm256_castpd256_pd128(b67)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c67_10 = _mm_add_pd(c67_10, _mm_mul_pd(a67_10, b67));
#endif
_mm_storeu_pd(&C[(i*88)+67], c67_10);
#else
C[(i*88)+35] += values[841] * B[(i*88)+67];
C[(i*88)+36] += values[842] * B[(i*88)+67];
C[(i*88)+37] += values[843] * B[(i*88)+67];
C[(i*88)+38] += values[844] * B[(i*88)+67];
C[(i*88)+39] += values[845] * B[(i*88)+67];
C[(i*88)+40] += values[846] * B[(i*88)+67];
C[(i*88)+63] += values[847] * B[(i*88)+67];
C[(i*88)+64] += values[848] * B[(i*88)+67];
C[(i*88)+65] += values[849] * B[(i*88)+67];
C[(i*88)+66] += values[850] * B[(i*88)+67];
C[(i*88)+67] += values[851] * B[(i*88)+67];
C[(i*88)+68] += values[852] * B[(i*88)+67];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b68 = _mm256_broadcast_sd(&B[(i*88)+68]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b68 = _mm_loaddup_pd(&B[(i*88)+68]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c68_0 = _mm256_loadu_pd(&C[(i*88)+35]);
__m256d a68_0 = _mm256_loadu_pd(&values[853]);
c68_0 = _mm256_add_pd(c68_0, _mm256_mul_pd(a68_0, b68));
_mm256_storeu_pd(&C[(i*88)+35], c68_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c68_0 = _mm_loadu_pd(&C[(i*88)+35]);
__m128d a68_0 = _mm_loadu_pd(&values[853]);
c68_0 = _mm_add_pd(c68_0, _mm_mul_pd(a68_0, b68));
_mm_storeu_pd(&C[(i*88)+35], c68_0);
__m128d c68_2 = _mm_loadu_pd(&C[(i*88)+37]);
__m128d a68_2 = _mm_loadu_pd(&values[855]);
c68_2 = _mm_add_pd(c68_2, _mm_mul_pd(a68_2, b68));
_mm_storeu_pd(&C[(i*88)+37], c68_2);
#endif
__m128d c68_4 = _mm_loadu_pd(&C[(i*88)+39]);
__m128d a68_4 = _mm_loadu_pd(&values[857]);
#if defined(__SSE3__) && defined(__AVX__)
c68_4 = _mm_add_pd(c68_4, _mm_mul_pd(a68_4, _mm256_castpd256_pd128(b68)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c68_4 = _mm_add_pd(c68_4, _mm_mul_pd(a68_4, b68));
#endif
_mm_storeu_pd(&C[(i*88)+39], c68_4);
#if defined(__SSE3__) && defined(__AVX__)
__m256d c68_6 = _mm256_loadu_pd(&C[(i*88)+63]);
__m256d a68_6 = _mm256_loadu_pd(&values[859]);
c68_6 = _mm256_add_pd(c68_6, _mm256_mul_pd(a68_6, b68));
_mm256_storeu_pd(&C[(i*88)+63], c68_6);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c68_6 = _mm_loadu_pd(&C[(i*88)+63]);
__m128d a68_6 = _mm_loadu_pd(&values[859]);
c68_6 = _mm_add_pd(c68_6, _mm_mul_pd(a68_6, b68));
_mm_storeu_pd(&C[(i*88)+63], c68_6);
__m128d c68_8 = _mm_loadu_pd(&C[(i*88)+65]);
__m128d a68_8 = _mm_loadu_pd(&values[861]);
c68_8 = _mm_add_pd(c68_8, _mm_mul_pd(a68_8, b68));
_mm_storeu_pd(&C[(i*88)+65], c68_8);
#endif
__m128d c68_10 = _mm_loadu_pd(&C[(i*88)+67]);
__m128d a68_10 = _mm_loadu_pd(&values[863]);
#if defined(__SSE3__) && defined(__AVX__)
c68_10 = _mm_add_pd(c68_10, _mm_mul_pd(a68_10, _mm256_castpd256_pd128(b68)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c68_10 = _mm_add_pd(c68_10, _mm_mul_pd(a68_10, b68));
#endif
_mm_storeu_pd(&C[(i*88)+67], c68_10);
#else
C[(i*88)+35] += values[853] * B[(i*88)+68];
C[(i*88)+36] += values[854] * B[(i*88)+68];
C[(i*88)+37] += values[855] * B[(i*88)+68];
C[(i*88)+38] += values[856] * B[(i*88)+68];
C[(i*88)+39] += values[857] * B[(i*88)+68];
C[(i*88)+40] += values[858] * B[(i*88)+68];
C[(i*88)+63] += values[859] * B[(i*88)+68];
C[(i*88)+64] += values[860] * B[(i*88)+68];
C[(i*88)+65] += values[861] * B[(i*88)+68];
C[(i*88)+66] += values[862] * B[(i*88)+68];
C[(i*88)+67] += values[863] * B[(i*88)+68];
C[(i*88)+68] += values[864] * B[(i*88)+68];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b69 = _mm256_broadcast_sd(&B[(i*88)+69]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b69 = _mm_loaddup_pd(&B[(i*88)+69]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c69_0 = _mm256_loadu_pd(&C[(i*88)+20]);
__m256d a69_0 = _mm256_loadu_pd(&values[865]);
c69_0 = _mm256_add_pd(c69_0, _mm256_mul_pd(a69_0, b69));
_mm256_storeu_pd(&C[(i*88)+20], c69_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c69_0 = _mm_loadu_pd(&C[(i*88)+20]);
__m128d a69_0 = _mm_loadu_pd(&values[865]);
c69_0 = _mm_add_pd(c69_0, _mm_mul_pd(a69_0, b69));
_mm_storeu_pd(&C[(i*88)+20], c69_0);
__m128d c69_2 = _mm_loadu_pd(&C[(i*88)+22]);
__m128d a69_2 = _mm_loadu_pd(&values[867]);
c69_2 = _mm_add_pd(c69_2, _mm_mul_pd(a69_2, b69));
_mm_storeu_pd(&C[(i*88)+22], c69_2);
#endif
__m128d c69_4 = _mm_load_sd(&C[(i*88)+24]);
__m128d a69_4 = _mm_load_sd(&values[869]);
#if defined(__SSE3__) && defined(__AVX__)
c69_4 = _mm_add_sd(c69_4, _mm_mul_sd(a69_4, _mm256_castpd256_pd128(b69)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c69_4 = _mm_add_sd(c69_4, _mm_mul_sd(a69_4, b69));
#endif
_mm_store_sd(&C[(i*88)+24], c69_4);
#if defined(__SSE3__) && defined(__AVX__)
__m256d c69_5 = _mm256_loadu_pd(&C[(i*88)+41]);
__m256d a69_5 = _mm256_loadu_pd(&values[870]);
c69_5 = _mm256_add_pd(c69_5, _mm256_mul_pd(a69_5, b69));
_mm256_storeu_pd(&C[(i*88)+41], c69_5);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c69_5 = _mm_loadu_pd(&C[(i*88)+41]);
__m128d a69_5 = _mm_loadu_pd(&values[870]);
c69_5 = _mm_add_pd(c69_5, _mm_mul_pd(a69_5, b69));
_mm_storeu_pd(&C[(i*88)+41], c69_5);
__m128d c69_7 = _mm_loadu_pd(&C[(i*88)+43]);
__m128d a69_7 = _mm_loadu_pd(&values[872]);
c69_7 = _mm_add_pd(c69_7, _mm_mul_pd(a69_7, b69));
_mm_storeu_pd(&C[(i*88)+43], c69_7);
#endif
__m128d c69_9 = _mm_load_sd(&C[(i*88)+45]);
__m128d a69_9 = _mm_load_sd(&values[874]);
#if defined(__SSE3__) && defined(__AVX__)
c69_9 = _mm_add_sd(c69_9, _mm_mul_sd(a69_9, _mm256_castpd256_pd128(b69)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c69_9 = _mm_add_sd(c69_9, _mm_mul_sd(a69_9, b69));
#endif
_mm_store_sd(&C[(i*88)+45], c69_9);
#if defined(__SSE3__) && defined(__AVX__)
__m256d c69_10 = _mm256_loadu_pd(&C[(i*88)+69]);
__m256d a69_10 = _mm256_loadu_pd(&values[875]);
c69_10 = _mm256_add_pd(c69_10, _mm256_mul_pd(a69_10, b69));
_mm256_storeu_pd(&C[(i*88)+69], c69_10);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c69_10 = _mm_loadu_pd(&C[(i*88)+69]);
__m128d a69_10 = _mm_loadu_pd(&values[875]);
c69_10 = _mm_add_pd(c69_10, _mm_mul_pd(a69_10, b69));
_mm_storeu_pd(&C[(i*88)+69], c69_10);
__m128d c69_12 = _mm_loadu_pd(&C[(i*88)+71]);
__m128d a69_12 = _mm_loadu_pd(&values[877]);
c69_12 = _mm_add_pd(c69_12, _mm_mul_pd(a69_12, b69));
_mm_storeu_pd(&C[(i*88)+71], c69_12);
#endif
__m128d c69_14 = _mm_load_sd(&C[(i*88)+73]);
__m128d a69_14 = _mm_load_sd(&values[879]);
#if defined(__SSE3__) && defined(__AVX__)
c69_14 = _mm_add_sd(c69_14, _mm_mul_sd(a69_14, _mm256_castpd256_pd128(b69)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c69_14 = _mm_add_sd(c69_14, _mm_mul_sd(a69_14, b69));
#endif
_mm_store_sd(&C[(i*88)+73], c69_14);
#else
C[(i*88)+20] += values[865] * B[(i*88)+69];
C[(i*88)+21] += values[866] * B[(i*88)+69];
C[(i*88)+22] += values[867] * B[(i*88)+69];
C[(i*88)+23] += values[868] * B[(i*88)+69];
C[(i*88)+24] += values[869] * B[(i*88)+69];
C[(i*88)+41] += values[870] * B[(i*88)+69];
C[(i*88)+42] += values[871] * B[(i*88)+69];
C[(i*88)+43] += values[872] * B[(i*88)+69];
C[(i*88)+44] += values[873] * B[(i*88)+69];
C[(i*88)+45] += values[874] * B[(i*88)+69];
C[(i*88)+69] += values[875] * B[(i*88)+69];
C[(i*88)+70] += values[876] * B[(i*88)+69];
C[(i*88)+71] += values[877] * B[(i*88)+69];
C[(i*88)+72] += values[878] * B[(i*88)+69];
C[(i*88)+73] += values[879] * B[(i*88)+69];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b70 = _mm256_broadcast_sd(&B[(i*88)+70]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b70 = _mm_loaddup_pd(&B[(i*88)+70]);
#endif
__m128d c70_0 = _mm_loadu_pd(&C[(i*88)+20]);
__m128d a70_0 = _mm_loadu_pd(&values[880]);
#if defined(__SSE3__) && defined(__AVX__)
c70_0 = _mm_add_pd(c70_0, _mm_mul_pd(a70_0, _mm256_castpd256_pd128(b70)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c70_0 = _mm_add_pd(c70_0, _mm_mul_pd(a70_0, b70));
#endif
_mm_storeu_pd(&C[(i*88)+20], c70_0);
__m128d c70_2 = _mm_load_sd(&C[(i*88)+22]);
__m128d a70_2 = _mm_load_sd(&values[882]);
#if defined(__SSE3__) && defined(__AVX__)
c70_2 = _mm_add_sd(c70_2, _mm_mul_sd(a70_2, _mm256_castpd256_pd128(b70)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c70_2 = _mm_add_sd(c70_2, _mm_mul_sd(a70_2, b70));
#endif
_mm_store_sd(&C[(i*88)+22], c70_2);
__m128d c70_3 = _mm_load_sd(&C[(i*88)+24]);
__m128d a70_3 = _mm_load_sd(&values[883]);
#if defined(__SSE3__) && defined(__AVX__)
c70_3 = _mm_add_sd(c70_3, _mm_mul_sd(a70_3, _mm256_castpd256_pd128(b70)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c70_3 = _mm_add_sd(c70_3, _mm_mul_sd(a70_3, b70));
#endif
_mm_store_sd(&C[(i*88)+24], c70_3);
__m128d c70_4 = _mm_loadu_pd(&C[(i*88)+41]);
__m128d a70_4 = _mm_loadu_pd(&values[884]);
#if defined(__SSE3__) && defined(__AVX__)
c70_4 = _mm_add_pd(c70_4, _mm_mul_pd(a70_4, _mm256_castpd256_pd128(b70)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c70_4 = _mm_add_pd(c70_4, _mm_mul_pd(a70_4, b70));
#endif
_mm_storeu_pd(&C[(i*88)+41], c70_4);
__m128d c70_6 = _mm_load_sd(&C[(i*88)+43]);
__m128d a70_6 = _mm_load_sd(&values[886]);
#if defined(__SSE3__) && defined(__AVX__)
c70_6 = _mm_add_sd(c70_6, _mm_mul_sd(a70_6, _mm256_castpd256_pd128(b70)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c70_6 = _mm_add_sd(c70_6, _mm_mul_sd(a70_6, b70));
#endif
_mm_store_sd(&C[(i*88)+43], c70_6);
__m128d c70_7 = _mm_load_sd(&C[(i*88)+45]);
__m128d a70_7 = _mm_load_sd(&values[887]);
#if defined(__SSE3__) && defined(__AVX__)
c70_7 = _mm_add_sd(c70_7, _mm_mul_sd(a70_7, _mm256_castpd256_pd128(b70)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c70_7 = _mm_add_sd(c70_7, _mm_mul_sd(a70_7, b70));
#endif
_mm_store_sd(&C[(i*88)+45], c70_7);
__m128d c70_8 = _mm_loadu_pd(&C[(i*88)+69]);
__m128d a70_8 = _mm_loadu_pd(&values[888]);
#if defined(__SSE3__) && defined(__AVX__)
c70_8 = _mm_add_pd(c70_8, _mm_mul_pd(a70_8, _mm256_castpd256_pd128(b70)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c70_8 = _mm_add_pd(c70_8, _mm_mul_pd(a70_8, b70));
#endif
_mm_storeu_pd(&C[(i*88)+69], c70_8);
__m128d c70_10 = _mm_load_sd(&C[(i*88)+71]);
__m128d a70_10 = _mm_load_sd(&values[890]);
#if defined(__SSE3__) && defined(__AVX__)
c70_10 = _mm_add_sd(c70_10, _mm_mul_sd(a70_10, _mm256_castpd256_pd128(b70)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c70_10 = _mm_add_sd(c70_10, _mm_mul_sd(a70_10, b70));
#endif
_mm_store_sd(&C[(i*88)+71], c70_10);
__m128d c70_11 = _mm_load_sd(&C[(i*88)+73]);
__m128d a70_11 = _mm_load_sd(&values[891]);
#if defined(__SSE3__) && defined(__AVX__)
c70_11 = _mm_add_sd(c70_11, _mm_mul_sd(a70_11, _mm256_castpd256_pd128(b70)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c70_11 = _mm_add_sd(c70_11, _mm_mul_sd(a70_11, b70));
#endif
_mm_store_sd(&C[(i*88)+73], c70_11);
#else
C[(i*88)+20] += values[880] * B[(i*88)+70];
C[(i*88)+21] += values[881] * B[(i*88)+70];
C[(i*88)+22] += values[882] * B[(i*88)+70];
C[(i*88)+24] += values[883] * B[(i*88)+70];
C[(i*88)+41] += values[884] * B[(i*88)+70];
C[(i*88)+42] += values[885] * B[(i*88)+70];
C[(i*88)+43] += values[886] * B[(i*88)+70];
C[(i*88)+45] += values[887] * B[(i*88)+70];
C[(i*88)+69] += values[888] * B[(i*88)+70];
C[(i*88)+70] += values[889] * B[(i*88)+70];
C[(i*88)+71] += values[890] * B[(i*88)+70];
C[(i*88)+73] += values[891] * B[(i*88)+70];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b71 = _mm256_broadcast_sd(&B[(i*88)+71]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b71 = _mm_loaddup_pd(&B[(i*88)+71]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c71_0 = _mm256_loadu_pd(&C[(i*88)+20]);
__m256d a71_0 = _mm256_loadu_pd(&values[892]);
c71_0 = _mm256_add_pd(c71_0, _mm256_mul_pd(a71_0, b71));
_mm256_storeu_pd(&C[(i*88)+20], c71_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c71_0 = _mm_loadu_pd(&C[(i*88)+20]);
__m128d a71_0 = _mm_loadu_pd(&values[892]);
c71_0 = _mm_add_pd(c71_0, _mm_mul_pd(a71_0, b71));
_mm_storeu_pd(&C[(i*88)+20], c71_0);
__m128d c71_2 = _mm_loadu_pd(&C[(i*88)+22]);
__m128d a71_2 = _mm_loadu_pd(&values[894]);
c71_2 = _mm_add_pd(c71_2, _mm_mul_pd(a71_2, b71));
_mm_storeu_pd(&C[(i*88)+22], c71_2);
#endif
__m128d c71_4 = _mm_load_sd(&C[(i*88)+24]);
__m128d a71_4 = _mm_load_sd(&values[896]);
#if defined(__SSE3__) && defined(__AVX__)
c71_4 = _mm_add_sd(c71_4, _mm_mul_sd(a71_4, _mm256_castpd256_pd128(b71)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c71_4 = _mm_add_sd(c71_4, _mm_mul_sd(a71_4, b71));
#endif
_mm_store_sd(&C[(i*88)+24], c71_4);
#if defined(__SSE3__) && defined(__AVX__)
__m256d c71_5 = _mm256_loadu_pd(&C[(i*88)+41]);
__m256d a71_5 = _mm256_loadu_pd(&values[897]);
c71_5 = _mm256_add_pd(c71_5, _mm256_mul_pd(a71_5, b71));
_mm256_storeu_pd(&C[(i*88)+41], c71_5);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c71_5 = _mm_loadu_pd(&C[(i*88)+41]);
__m128d a71_5 = _mm_loadu_pd(&values[897]);
c71_5 = _mm_add_pd(c71_5, _mm_mul_pd(a71_5, b71));
_mm_storeu_pd(&C[(i*88)+41], c71_5);
__m128d c71_7 = _mm_loadu_pd(&C[(i*88)+43]);
__m128d a71_7 = _mm_loadu_pd(&values[899]);
c71_7 = _mm_add_pd(c71_7, _mm_mul_pd(a71_7, b71));
_mm_storeu_pd(&C[(i*88)+43], c71_7);
#endif
__m128d c71_9 = _mm_load_sd(&C[(i*88)+45]);
__m128d a71_9 = _mm_load_sd(&values[901]);
#if defined(__SSE3__) && defined(__AVX__)
c71_9 = _mm_add_sd(c71_9, _mm_mul_sd(a71_9, _mm256_castpd256_pd128(b71)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c71_9 = _mm_add_sd(c71_9, _mm_mul_sd(a71_9, b71));
#endif
_mm_store_sd(&C[(i*88)+45], c71_9);
#if defined(__SSE3__) && defined(__AVX__)
__m256d c71_10 = _mm256_loadu_pd(&C[(i*88)+69]);
__m256d a71_10 = _mm256_loadu_pd(&values[902]);
c71_10 = _mm256_add_pd(c71_10, _mm256_mul_pd(a71_10, b71));
_mm256_storeu_pd(&C[(i*88)+69], c71_10);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c71_10 = _mm_loadu_pd(&C[(i*88)+69]);
__m128d a71_10 = _mm_loadu_pd(&values[902]);
c71_10 = _mm_add_pd(c71_10, _mm_mul_pd(a71_10, b71));
_mm_storeu_pd(&C[(i*88)+69], c71_10);
__m128d c71_12 = _mm_loadu_pd(&C[(i*88)+71]);
__m128d a71_12 = _mm_loadu_pd(&values[904]);
c71_12 = _mm_add_pd(c71_12, _mm_mul_pd(a71_12, b71));
_mm_storeu_pd(&C[(i*88)+71], c71_12);
#endif
__m128d c71_14 = _mm_load_sd(&C[(i*88)+73]);
__m128d a71_14 = _mm_load_sd(&values[906]);
#if defined(__SSE3__) && defined(__AVX__)
c71_14 = _mm_add_sd(c71_14, _mm_mul_sd(a71_14, _mm256_castpd256_pd128(b71)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c71_14 = _mm_add_sd(c71_14, _mm_mul_sd(a71_14, b71));
#endif
_mm_store_sd(&C[(i*88)+73], c71_14);
#else
C[(i*88)+20] += values[892] * B[(i*88)+71];
C[(i*88)+21] += values[893] * B[(i*88)+71];
C[(i*88)+22] += values[894] * B[(i*88)+71];
C[(i*88)+23] += values[895] * B[(i*88)+71];
C[(i*88)+24] += values[896] * B[(i*88)+71];
C[(i*88)+41] += values[897] * B[(i*88)+71];
C[(i*88)+42] += values[898] * B[(i*88)+71];
C[(i*88)+43] += values[899] * B[(i*88)+71];
C[(i*88)+44] += values[900] * B[(i*88)+71];
C[(i*88)+45] += values[901] * B[(i*88)+71];
C[(i*88)+69] += values[902] * B[(i*88)+71];
C[(i*88)+70] += values[903] * B[(i*88)+71];
C[(i*88)+71] += values[904] * B[(i*88)+71];
C[(i*88)+72] += values[905] * B[(i*88)+71];
C[(i*88)+73] += values[906] * B[(i*88)+71];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b72 = _mm256_broadcast_sd(&B[(i*88)+72]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b72 = _mm_loaddup_pd(&B[(i*88)+72]);
#endif
__m128d c72_0 = _mm_load_sd(&C[(i*88)+20]);
__m128d a72_0 = _mm_load_sd(&values[907]);
#if defined(__SSE3__) && defined(__AVX__)
c72_0 = _mm_add_sd(c72_0, _mm_mul_sd(a72_0, _mm256_castpd256_pd128(b72)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c72_0 = _mm_add_sd(c72_0, _mm_mul_sd(a72_0, b72));
#endif
_mm_store_sd(&C[(i*88)+20], c72_0);
__m128d c72_1 = _mm_loadu_pd(&C[(i*88)+22]);
__m128d a72_1 = _mm_loadu_pd(&values[908]);
#if defined(__SSE3__) && defined(__AVX__)
c72_1 = _mm_add_pd(c72_1, _mm_mul_pd(a72_1, _mm256_castpd256_pd128(b72)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c72_1 = _mm_add_pd(c72_1, _mm_mul_pd(a72_1, b72));
#endif
_mm_storeu_pd(&C[(i*88)+22], c72_1);
__m128d c72_3 = _mm_load_sd(&C[(i*88)+24]);
__m128d a72_3 = _mm_load_sd(&values[910]);
#if defined(__SSE3__) && defined(__AVX__)
c72_3 = _mm_add_sd(c72_3, _mm_mul_sd(a72_3, _mm256_castpd256_pd128(b72)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c72_3 = _mm_add_sd(c72_3, _mm_mul_sd(a72_3, b72));
#endif
_mm_store_sd(&C[(i*88)+24], c72_3);
__m128d c72_4 = _mm_load_sd(&C[(i*88)+41]);
__m128d a72_4 = _mm_load_sd(&values[911]);
#if defined(__SSE3__) && defined(__AVX__)
c72_4 = _mm_add_sd(c72_4, _mm_mul_sd(a72_4, _mm256_castpd256_pd128(b72)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c72_4 = _mm_add_sd(c72_4, _mm_mul_sd(a72_4, b72));
#endif
_mm_store_sd(&C[(i*88)+41], c72_4);
__m128d c72_5 = _mm_loadu_pd(&C[(i*88)+43]);
__m128d a72_5 = _mm_loadu_pd(&values[912]);
#if defined(__SSE3__) && defined(__AVX__)
c72_5 = _mm_add_pd(c72_5, _mm_mul_pd(a72_5, _mm256_castpd256_pd128(b72)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c72_5 = _mm_add_pd(c72_5, _mm_mul_pd(a72_5, b72));
#endif
_mm_storeu_pd(&C[(i*88)+43], c72_5);
__m128d c72_7 = _mm_load_sd(&C[(i*88)+45]);
__m128d a72_7 = _mm_load_sd(&values[914]);
#if defined(__SSE3__) && defined(__AVX__)
c72_7 = _mm_add_sd(c72_7, _mm_mul_sd(a72_7, _mm256_castpd256_pd128(b72)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c72_7 = _mm_add_sd(c72_7, _mm_mul_sd(a72_7, b72));
#endif
_mm_store_sd(&C[(i*88)+45], c72_7);
__m128d c72_8 = _mm_load_sd(&C[(i*88)+69]);
__m128d a72_8 = _mm_load_sd(&values[915]);
#if defined(__SSE3__) && defined(__AVX__)
c72_8 = _mm_add_sd(c72_8, _mm_mul_sd(a72_8, _mm256_castpd256_pd128(b72)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c72_8 = _mm_add_sd(c72_8, _mm_mul_sd(a72_8, b72));
#endif
_mm_store_sd(&C[(i*88)+69], c72_8);
__m128d c72_9 = _mm_loadu_pd(&C[(i*88)+71]);
__m128d a72_9 = _mm_loadu_pd(&values[916]);
#if defined(__SSE3__) && defined(__AVX__)
c72_9 = _mm_add_pd(c72_9, _mm_mul_pd(a72_9, _mm256_castpd256_pd128(b72)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c72_9 = _mm_add_pd(c72_9, _mm_mul_pd(a72_9, b72));
#endif
_mm_storeu_pd(&C[(i*88)+71], c72_9);
__m128d c72_11 = _mm_load_sd(&C[(i*88)+73]);
__m128d a72_11 = _mm_load_sd(&values[918]);
#if defined(__SSE3__) && defined(__AVX__)
c72_11 = _mm_add_sd(c72_11, _mm_mul_sd(a72_11, _mm256_castpd256_pd128(b72)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c72_11 = _mm_add_sd(c72_11, _mm_mul_sd(a72_11, b72));
#endif
_mm_store_sd(&C[(i*88)+73], c72_11);
#else
C[(i*88)+20] += values[907] * B[(i*88)+72];
C[(i*88)+22] += values[908] * B[(i*88)+72];
C[(i*88)+23] += values[909] * B[(i*88)+72];
C[(i*88)+24] += values[910] * B[(i*88)+72];
C[(i*88)+41] += values[911] * B[(i*88)+72];
C[(i*88)+43] += values[912] * B[(i*88)+72];
C[(i*88)+44] += values[913] * B[(i*88)+72];
C[(i*88)+45] += values[914] * B[(i*88)+72];
C[(i*88)+69] += values[915] * B[(i*88)+72];
C[(i*88)+71] += values[916] * B[(i*88)+72];
C[(i*88)+72] += values[917] * B[(i*88)+72];
C[(i*88)+73] += values[918] * B[(i*88)+72];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b73 = _mm256_broadcast_sd(&B[(i*88)+73]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b73 = _mm_loaddup_pd(&B[(i*88)+73]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c73_0 = _mm256_loadu_pd(&C[(i*88)+20]);
__m256d a73_0 = _mm256_loadu_pd(&values[919]);
c73_0 = _mm256_add_pd(c73_0, _mm256_mul_pd(a73_0, b73));
_mm256_storeu_pd(&C[(i*88)+20], c73_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c73_0 = _mm_loadu_pd(&C[(i*88)+20]);
__m128d a73_0 = _mm_loadu_pd(&values[919]);
c73_0 = _mm_add_pd(c73_0, _mm_mul_pd(a73_0, b73));
_mm_storeu_pd(&C[(i*88)+20], c73_0);
__m128d c73_2 = _mm_loadu_pd(&C[(i*88)+22]);
__m128d a73_2 = _mm_loadu_pd(&values[921]);
c73_2 = _mm_add_pd(c73_2, _mm_mul_pd(a73_2, b73));
_mm_storeu_pd(&C[(i*88)+22], c73_2);
#endif
__m128d c73_4 = _mm_load_sd(&C[(i*88)+24]);
__m128d a73_4 = _mm_load_sd(&values[923]);
#if defined(__SSE3__) && defined(__AVX__)
c73_4 = _mm_add_sd(c73_4, _mm_mul_sd(a73_4, _mm256_castpd256_pd128(b73)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c73_4 = _mm_add_sd(c73_4, _mm_mul_sd(a73_4, b73));
#endif
_mm_store_sd(&C[(i*88)+24], c73_4);
#if defined(__SSE3__) && defined(__AVX__)
__m256d c73_5 = _mm256_loadu_pd(&C[(i*88)+41]);
__m256d a73_5 = _mm256_loadu_pd(&values[924]);
c73_5 = _mm256_add_pd(c73_5, _mm256_mul_pd(a73_5, b73));
_mm256_storeu_pd(&C[(i*88)+41], c73_5);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c73_5 = _mm_loadu_pd(&C[(i*88)+41]);
__m128d a73_5 = _mm_loadu_pd(&values[924]);
c73_5 = _mm_add_pd(c73_5, _mm_mul_pd(a73_5, b73));
_mm_storeu_pd(&C[(i*88)+41], c73_5);
__m128d c73_7 = _mm_loadu_pd(&C[(i*88)+43]);
__m128d a73_7 = _mm_loadu_pd(&values[926]);
c73_7 = _mm_add_pd(c73_7, _mm_mul_pd(a73_7, b73));
_mm_storeu_pd(&C[(i*88)+43], c73_7);
#endif
__m128d c73_9 = _mm_load_sd(&C[(i*88)+45]);
__m128d a73_9 = _mm_load_sd(&values[928]);
#if defined(__SSE3__) && defined(__AVX__)
c73_9 = _mm_add_sd(c73_9, _mm_mul_sd(a73_9, _mm256_castpd256_pd128(b73)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c73_9 = _mm_add_sd(c73_9, _mm_mul_sd(a73_9, b73));
#endif
_mm_store_sd(&C[(i*88)+45], c73_9);
#if defined(__SSE3__) && defined(__AVX__)
__m256d c73_10 = _mm256_loadu_pd(&C[(i*88)+69]);
__m256d a73_10 = _mm256_loadu_pd(&values[929]);
c73_10 = _mm256_add_pd(c73_10, _mm256_mul_pd(a73_10, b73));
_mm256_storeu_pd(&C[(i*88)+69], c73_10);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c73_10 = _mm_loadu_pd(&C[(i*88)+69]);
__m128d a73_10 = _mm_loadu_pd(&values[929]);
c73_10 = _mm_add_pd(c73_10, _mm_mul_pd(a73_10, b73));
_mm_storeu_pd(&C[(i*88)+69], c73_10);
__m128d c73_12 = _mm_loadu_pd(&C[(i*88)+71]);
__m128d a73_12 = _mm_loadu_pd(&values[931]);
c73_12 = _mm_add_pd(c73_12, _mm_mul_pd(a73_12, b73));
_mm_storeu_pd(&C[(i*88)+71], c73_12);
#endif
__m128d c73_14 = _mm_load_sd(&C[(i*88)+73]);
__m128d a73_14 = _mm_load_sd(&values[933]);
#if defined(__SSE3__) && defined(__AVX__)
c73_14 = _mm_add_sd(c73_14, _mm_mul_sd(a73_14, _mm256_castpd256_pd128(b73)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c73_14 = _mm_add_sd(c73_14, _mm_mul_sd(a73_14, b73));
#endif
_mm_store_sd(&C[(i*88)+73], c73_14);
#else
C[(i*88)+20] += values[919] * B[(i*88)+73];
C[(i*88)+21] += values[920] * B[(i*88)+73];
C[(i*88)+22] += values[921] * B[(i*88)+73];
C[(i*88)+23] += values[922] * B[(i*88)+73];
C[(i*88)+24] += values[923] * B[(i*88)+73];
C[(i*88)+41] += values[924] * B[(i*88)+73];
C[(i*88)+42] += values[925] * B[(i*88)+73];
C[(i*88)+43] += values[926] * B[(i*88)+73];
C[(i*88)+44] += values[927] * B[(i*88)+73];
C[(i*88)+45] += values[928] * B[(i*88)+73];
C[(i*88)+69] += values[929] * B[(i*88)+73];
C[(i*88)+70] += values[930] * B[(i*88)+73];
C[(i*88)+71] += values[931] * B[(i*88)+73];
C[(i*88)+72] += values[932] * B[(i*88)+73];
C[(i*88)+73] += values[933] * B[(i*88)+73];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b74 = _mm256_broadcast_sd(&B[(i*88)+74]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b74 = _mm_loaddup_pd(&B[(i*88)+74]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c74_0 = _mm256_loadu_pd(&C[(i*88)+10]);
__m256d a74_0 = _mm256_loadu_pd(&values[934]);
c74_0 = _mm256_add_pd(c74_0, _mm256_mul_pd(a74_0, b74));
_mm256_storeu_pd(&C[(i*88)+10], c74_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c74_0 = _mm_loadu_pd(&C[(i*88)+10]);
__m128d a74_0 = _mm_loadu_pd(&values[934]);
c74_0 = _mm_add_pd(c74_0, _mm_mul_pd(a74_0, b74));
_mm_storeu_pd(&C[(i*88)+10], c74_0);
__m128d c74_2 = _mm_loadu_pd(&C[(i*88)+12]);
__m128d a74_2 = _mm_loadu_pd(&values[936]);
c74_2 = _mm_add_pd(c74_2, _mm_mul_pd(a74_2, b74));
_mm_storeu_pd(&C[(i*88)+12], c74_2);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c74_4 = _mm256_loadu_pd(&C[(i*88)+25]);
__m256d a74_4 = _mm256_loadu_pd(&values[938]);
c74_4 = _mm256_add_pd(c74_4, _mm256_mul_pd(a74_4, b74));
_mm256_storeu_pd(&C[(i*88)+25], c74_4);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c74_4 = _mm_loadu_pd(&C[(i*88)+25]);
__m128d a74_4 = _mm_loadu_pd(&values[938]);
c74_4 = _mm_add_pd(c74_4, _mm_mul_pd(a74_4, b74));
_mm_storeu_pd(&C[(i*88)+25], c74_4);
__m128d c74_6 = _mm_loadu_pd(&C[(i*88)+27]);
__m128d a74_6 = _mm_loadu_pd(&values[940]);
c74_6 = _mm_add_pd(c74_6, _mm_mul_pd(a74_6, b74));
_mm_storeu_pd(&C[(i*88)+27], c74_6);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c74_8 = _mm256_loadu_pd(&C[(i*88)+46]);
__m256d a74_8 = _mm256_loadu_pd(&values[942]);
c74_8 = _mm256_add_pd(c74_8, _mm256_mul_pd(a74_8, b74));
_mm256_storeu_pd(&C[(i*88)+46], c74_8);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c74_8 = _mm_loadu_pd(&C[(i*88)+46]);
__m128d a74_8 = _mm_loadu_pd(&values[942]);
c74_8 = _mm_add_pd(c74_8, _mm_mul_pd(a74_8, b74));
_mm_storeu_pd(&C[(i*88)+46], c74_8);
__m128d c74_10 = _mm_loadu_pd(&C[(i*88)+48]);
__m128d a74_10 = _mm_loadu_pd(&values[944]);
c74_10 = _mm_add_pd(c74_10, _mm_mul_pd(a74_10, b74));
_mm_storeu_pd(&C[(i*88)+48], c74_10);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c74_12 = _mm256_loadu_pd(&C[(i*88)+74]);
__m256d a74_12 = _mm256_loadu_pd(&values[946]);
c74_12 = _mm256_add_pd(c74_12, _mm256_mul_pd(a74_12, b74));
_mm256_storeu_pd(&C[(i*88)+74], c74_12);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c74_12 = _mm_loadu_pd(&C[(i*88)+74]);
__m128d a74_12 = _mm_loadu_pd(&values[946]);
c74_12 = _mm_add_pd(c74_12, _mm_mul_pd(a74_12, b74));
_mm_storeu_pd(&C[(i*88)+74], c74_12);
__m128d c74_14 = _mm_loadu_pd(&C[(i*88)+76]);
__m128d a74_14 = _mm_loadu_pd(&values[948]);
c74_14 = _mm_add_pd(c74_14, _mm_mul_pd(a74_14, b74));
_mm_storeu_pd(&C[(i*88)+76], c74_14);
#endif
#else
C[(i*88)+10] += values[934] * B[(i*88)+74];
C[(i*88)+11] += values[935] * B[(i*88)+74];
C[(i*88)+12] += values[936] * B[(i*88)+74];
C[(i*88)+13] += values[937] * B[(i*88)+74];
C[(i*88)+25] += values[938] * B[(i*88)+74];
C[(i*88)+26] += values[939] * B[(i*88)+74];
C[(i*88)+27] += values[940] * B[(i*88)+74];
C[(i*88)+28] += values[941] * B[(i*88)+74];
C[(i*88)+46] += values[942] * B[(i*88)+74];
C[(i*88)+47] += values[943] * B[(i*88)+74];
C[(i*88)+48] += values[944] * B[(i*88)+74];
C[(i*88)+49] += values[945] * B[(i*88)+74];
C[(i*88)+74] += values[946] * B[(i*88)+74];
C[(i*88)+75] += values[947] * B[(i*88)+74];
C[(i*88)+76] += values[948] * B[(i*88)+74];
C[(i*88)+77] += values[949] * B[(i*88)+74];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b75 = _mm256_broadcast_sd(&B[(i*88)+75]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b75 = _mm_loaddup_pd(&B[(i*88)+75]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c75_0 = _mm256_loadu_pd(&C[(i*88)+10]);
__m256d a75_0 = _mm256_loadu_pd(&values[950]);
c75_0 = _mm256_add_pd(c75_0, _mm256_mul_pd(a75_0, b75));
_mm256_storeu_pd(&C[(i*88)+10], c75_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c75_0 = _mm_loadu_pd(&C[(i*88)+10]);
__m128d a75_0 = _mm_loadu_pd(&values[950]);
c75_0 = _mm_add_pd(c75_0, _mm_mul_pd(a75_0, b75));
_mm_storeu_pd(&C[(i*88)+10], c75_0);
__m128d c75_2 = _mm_loadu_pd(&C[(i*88)+12]);
__m128d a75_2 = _mm_loadu_pd(&values[952]);
c75_2 = _mm_add_pd(c75_2, _mm_mul_pd(a75_2, b75));
_mm_storeu_pd(&C[(i*88)+12], c75_2);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c75_4 = _mm256_loadu_pd(&C[(i*88)+25]);
__m256d a75_4 = _mm256_loadu_pd(&values[954]);
c75_4 = _mm256_add_pd(c75_4, _mm256_mul_pd(a75_4, b75));
_mm256_storeu_pd(&C[(i*88)+25], c75_4);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c75_4 = _mm_loadu_pd(&C[(i*88)+25]);
__m128d a75_4 = _mm_loadu_pd(&values[954]);
c75_4 = _mm_add_pd(c75_4, _mm_mul_pd(a75_4, b75));
_mm_storeu_pd(&C[(i*88)+25], c75_4);
__m128d c75_6 = _mm_loadu_pd(&C[(i*88)+27]);
__m128d a75_6 = _mm_loadu_pd(&values[956]);
c75_6 = _mm_add_pd(c75_6, _mm_mul_pd(a75_6, b75));
_mm_storeu_pd(&C[(i*88)+27], c75_6);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c75_8 = _mm256_loadu_pd(&C[(i*88)+46]);
__m256d a75_8 = _mm256_loadu_pd(&values[958]);
c75_8 = _mm256_add_pd(c75_8, _mm256_mul_pd(a75_8, b75));
_mm256_storeu_pd(&C[(i*88)+46], c75_8);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c75_8 = _mm_loadu_pd(&C[(i*88)+46]);
__m128d a75_8 = _mm_loadu_pd(&values[958]);
c75_8 = _mm_add_pd(c75_8, _mm_mul_pd(a75_8, b75));
_mm_storeu_pd(&C[(i*88)+46], c75_8);
__m128d c75_10 = _mm_loadu_pd(&C[(i*88)+48]);
__m128d a75_10 = _mm_loadu_pd(&values[960]);
c75_10 = _mm_add_pd(c75_10, _mm_mul_pd(a75_10, b75));
_mm_storeu_pd(&C[(i*88)+48], c75_10);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c75_12 = _mm256_loadu_pd(&C[(i*88)+74]);
__m256d a75_12 = _mm256_loadu_pd(&values[962]);
c75_12 = _mm256_add_pd(c75_12, _mm256_mul_pd(a75_12, b75));
_mm256_storeu_pd(&C[(i*88)+74], c75_12);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c75_12 = _mm_loadu_pd(&C[(i*88)+74]);
__m128d a75_12 = _mm_loadu_pd(&values[962]);
c75_12 = _mm_add_pd(c75_12, _mm_mul_pd(a75_12, b75));
_mm_storeu_pd(&C[(i*88)+74], c75_12);
__m128d c75_14 = _mm_loadu_pd(&C[(i*88)+76]);
__m128d a75_14 = _mm_loadu_pd(&values[964]);
c75_14 = _mm_add_pd(c75_14, _mm_mul_pd(a75_14, b75));
_mm_storeu_pd(&C[(i*88)+76], c75_14);
#endif
#else
C[(i*88)+10] += values[950] * B[(i*88)+75];
C[(i*88)+11] += values[951] * B[(i*88)+75];
C[(i*88)+12] += values[952] * B[(i*88)+75];
C[(i*88)+13] += values[953] * B[(i*88)+75];
C[(i*88)+25] += values[954] * B[(i*88)+75];
C[(i*88)+26] += values[955] * B[(i*88)+75];
C[(i*88)+27] += values[956] * B[(i*88)+75];
C[(i*88)+28] += values[957] * B[(i*88)+75];
C[(i*88)+46] += values[958] * B[(i*88)+75];
C[(i*88)+47] += values[959] * B[(i*88)+75];
C[(i*88)+48] += values[960] * B[(i*88)+75];
C[(i*88)+49] += values[961] * B[(i*88)+75];
C[(i*88)+74] += values[962] * B[(i*88)+75];
C[(i*88)+75] += values[963] * B[(i*88)+75];
C[(i*88)+76] += values[964] * B[(i*88)+75];
C[(i*88)+77] += values[965] * B[(i*88)+75];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b76 = _mm256_broadcast_sd(&B[(i*88)+76]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b76 = _mm_loaddup_pd(&B[(i*88)+76]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c76_0 = _mm256_loadu_pd(&C[(i*88)+10]);
__m256d a76_0 = _mm256_loadu_pd(&values[966]);
c76_0 = _mm256_add_pd(c76_0, _mm256_mul_pd(a76_0, b76));
_mm256_storeu_pd(&C[(i*88)+10], c76_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c76_0 = _mm_loadu_pd(&C[(i*88)+10]);
__m128d a76_0 = _mm_loadu_pd(&values[966]);
c76_0 = _mm_add_pd(c76_0, _mm_mul_pd(a76_0, b76));
_mm_storeu_pd(&C[(i*88)+10], c76_0);
__m128d c76_2 = _mm_loadu_pd(&C[(i*88)+12]);
__m128d a76_2 = _mm_loadu_pd(&values[968]);
c76_2 = _mm_add_pd(c76_2, _mm_mul_pd(a76_2, b76));
_mm_storeu_pd(&C[(i*88)+12], c76_2);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c76_4 = _mm256_loadu_pd(&C[(i*88)+25]);
__m256d a76_4 = _mm256_loadu_pd(&values[970]);
c76_4 = _mm256_add_pd(c76_4, _mm256_mul_pd(a76_4, b76));
_mm256_storeu_pd(&C[(i*88)+25], c76_4);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c76_4 = _mm_loadu_pd(&C[(i*88)+25]);
__m128d a76_4 = _mm_loadu_pd(&values[970]);
c76_4 = _mm_add_pd(c76_4, _mm_mul_pd(a76_4, b76));
_mm_storeu_pd(&C[(i*88)+25], c76_4);
__m128d c76_6 = _mm_loadu_pd(&C[(i*88)+27]);
__m128d a76_6 = _mm_loadu_pd(&values[972]);
c76_6 = _mm_add_pd(c76_6, _mm_mul_pd(a76_6, b76));
_mm_storeu_pd(&C[(i*88)+27], c76_6);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c76_8 = _mm256_loadu_pd(&C[(i*88)+46]);
__m256d a76_8 = _mm256_loadu_pd(&values[974]);
c76_8 = _mm256_add_pd(c76_8, _mm256_mul_pd(a76_8, b76));
_mm256_storeu_pd(&C[(i*88)+46], c76_8);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c76_8 = _mm_loadu_pd(&C[(i*88)+46]);
__m128d a76_8 = _mm_loadu_pd(&values[974]);
c76_8 = _mm_add_pd(c76_8, _mm_mul_pd(a76_8, b76));
_mm_storeu_pd(&C[(i*88)+46], c76_8);
__m128d c76_10 = _mm_loadu_pd(&C[(i*88)+48]);
__m128d a76_10 = _mm_loadu_pd(&values[976]);
c76_10 = _mm_add_pd(c76_10, _mm_mul_pd(a76_10, b76));
_mm_storeu_pd(&C[(i*88)+48], c76_10);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c76_12 = _mm256_loadu_pd(&C[(i*88)+74]);
__m256d a76_12 = _mm256_loadu_pd(&values[978]);
c76_12 = _mm256_add_pd(c76_12, _mm256_mul_pd(a76_12, b76));
_mm256_storeu_pd(&C[(i*88)+74], c76_12);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c76_12 = _mm_loadu_pd(&C[(i*88)+74]);
__m128d a76_12 = _mm_loadu_pd(&values[978]);
c76_12 = _mm_add_pd(c76_12, _mm_mul_pd(a76_12, b76));
_mm_storeu_pd(&C[(i*88)+74], c76_12);
__m128d c76_14 = _mm_loadu_pd(&C[(i*88)+76]);
__m128d a76_14 = _mm_loadu_pd(&values[980]);
c76_14 = _mm_add_pd(c76_14, _mm_mul_pd(a76_14, b76));
_mm_storeu_pd(&C[(i*88)+76], c76_14);
#endif
#else
C[(i*88)+10] += values[966] * B[(i*88)+76];
C[(i*88)+11] += values[967] * B[(i*88)+76];
C[(i*88)+12] += values[968] * B[(i*88)+76];
C[(i*88)+13] += values[969] * B[(i*88)+76];
C[(i*88)+25] += values[970] * B[(i*88)+76];
C[(i*88)+26] += values[971] * B[(i*88)+76];
C[(i*88)+27] += values[972] * B[(i*88)+76];
C[(i*88)+28] += values[973] * B[(i*88)+76];
C[(i*88)+46] += values[974] * B[(i*88)+76];
C[(i*88)+47] += values[975] * B[(i*88)+76];
C[(i*88)+48] += values[976] * B[(i*88)+76];
C[(i*88)+49] += values[977] * B[(i*88)+76];
C[(i*88)+74] += values[978] * B[(i*88)+76];
C[(i*88)+75] += values[979] * B[(i*88)+76];
C[(i*88)+76] += values[980] * B[(i*88)+76];
C[(i*88)+77] += values[981] * B[(i*88)+76];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b77 = _mm256_broadcast_sd(&B[(i*88)+77]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b77 = _mm_loaddup_pd(&B[(i*88)+77]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c77_0 = _mm256_loadu_pd(&C[(i*88)+10]);
__m256d a77_0 = _mm256_loadu_pd(&values[982]);
c77_0 = _mm256_add_pd(c77_0, _mm256_mul_pd(a77_0, b77));
_mm256_storeu_pd(&C[(i*88)+10], c77_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c77_0 = _mm_loadu_pd(&C[(i*88)+10]);
__m128d a77_0 = _mm_loadu_pd(&values[982]);
c77_0 = _mm_add_pd(c77_0, _mm_mul_pd(a77_0, b77));
_mm_storeu_pd(&C[(i*88)+10], c77_0);
__m128d c77_2 = _mm_loadu_pd(&C[(i*88)+12]);
__m128d a77_2 = _mm_loadu_pd(&values[984]);
c77_2 = _mm_add_pd(c77_2, _mm_mul_pd(a77_2, b77));
_mm_storeu_pd(&C[(i*88)+12], c77_2);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c77_4 = _mm256_loadu_pd(&C[(i*88)+25]);
__m256d a77_4 = _mm256_loadu_pd(&values[986]);
c77_4 = _mm256_add_pd(c77_4, _mm256_mul_pd(a77_4, b77));
_mm256_storeu_pd(&C[(i*88)+25], c77_4);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c77_4 = _mm_loadu_pd(&C[(i*88)+25]);
__m128d a77_4 = _mm_loadu_pd(&values[986]);
c77_4 = _mm_add_pd(c77_4, _mm_mul_pd(a77_4, b77));
_mm_storeu_pd(&C[(i*88)+25], c77_4);
__m128d c77_6 = _mm_loadu_pd(&C[(i*88)+27]);
__m128d a77_6 = _mm_loadu_pd(&values[988]);
c77_6 = _mm_add_pd(c77_6, _mm_mul_pd(a77_6, b77));
_mm_storeu_pd(&C[(i*88)+27], c77_6);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c77_8 = _mm256_loadu_pd(&C[(i*88)+46]);
__m256d a77_8 = _mm256_loadu_pd(&values[990]);
c77_8 = _mm256_add_pd(c77_8, _mm256_mul_pd(a77_8, b77));
_mm256_storeu_pd(&C[(i*88)+46], c77_8);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c77_8 = _mm_loadu_pd(&C[(i*88)+46]);
__m128d a77_8 = _mm_loadu_pd(&values[990]);
c77_8 = _mm_add_pd(c77_8, _mm_mul_pd(a77_8, b77));
_mm_storeu_pd(&C[(i*88)+46], c77_8);
__m128d c77_10 = _mm_loadu_pd(&C[(i*88)+48]);
__m128d a77_10 = _mm_loadu_pd(&values[992]);
c77_10 = _mm_add_pd(c77_10, _mm_mul_pd(a77_10, b77));
_mm_storeu_pd(&C[(i*88)+48], c77_10);
#endif
#if defined(__SSE3__) && defined(__AVX__)
__m256d c77_12 = _mm256_loadu_pd(&C[(i*88)+74]);
__m256d a77_12 = _mm256_loadu_pd(&values[994]);
c77_12 = _mm256_add_pd(c77_12, _mm256_mul_pd(a77_12, b77));
_mm256_storeu_pd(&C[(i*88)+74], c77_12);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d c77_12 = _mm_loadu_pd(&C[(i*88)+74]);
__m128d a77_12 = _mm_loadu_pd(&values[994]);
c77_12 = _mm_add_pd(c77_12, _mm_mul_pd(a77_12, b77));
_mm_storeu_pd(&C[(i*88)+74], c77_12);
__m128d c77_14 = _mm_loadu_pd(&C[(i*88)+76]);
__m128d a77_14 = _mm_loadu_pd(&values[996]);
c77_14 = _mm_add_pd(c77_14, _mm_mul_pd(a77_14, b77));
_mm_storeu_pd(&C[(i*88)+76], c77_14);
#endif
#else
C[(i*88)+10] += values[982] * B[(i*88)+77];
C[(i*88)+11] += values[983] * B[(i*88)+77];
C[(i*88)+12] += values[984] * B[(i*88)+77];
C[(i*88)+13] += values[985] * B[(i*88)+77];
C[(i*88)+25] += values[986] * B[(i*88)+77];
C[(i*88)+26] += values[987] * B[(i*88)+77];
C[(i*88)+27] += values[988] * B[(i*88)+77];
C[(i*88)+28] += values[989] * B[(i*88)+77];
C[(i*88)+46] += values[990] * B[(i*88)+77];
C[(i*88)+47] += values[991] * B[(i*88)+77];
C[(i*88)+48] += values[992] * B[(i*88)+77];
C[(i*88)+49] += values[993] * B[(i*88)+77];
C[(i*88)+74] += values[994] * B[(i*88)+77];
C[(i*88)+75] += values[995] * B[(i*88)+77];
C[(i*88)+76] += values[996] * B[(i*88)+77];
C[(i*88)+77] += values[997] * B[(i*88)+77];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b78 = _mm256_broadcast_sd(&B[(i*88)+78]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b78 = _mm_loaddup_pd(&B[(i*88)+78]);
#endif
__m128d c78_0 = _mm_loadu_pd(&C[(i*88)+4]);
__m128d a78_0 = _mm_loadu_pd(&values[998]);
#if defined(__SSE3__) && defined(__AVX__)
c78_0 = _mm_add_pd(c78_0, _mm_mul_pd(a78_0, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c78_0 = _mm_add_pd(c78_0, _mm_mul_pd(a78_0, b78));
#endif
_mm_storeu_pd(&C[(i*88)+4], c78_0);
__m128d c78_2 = _mm_load_sd(&C[(i*88)+6]);
__m128d a78_2 = _mm_load_sd(&values[1000]);
#if defined(__SSE3__) && defined(__AVX__)
c78_2 = _mm_add_sd(c78_2, _mm_mul_sd(a78_2, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c78_2 = _mm_add_sd(c78_2, _mm_mul_sd(a78_2, b78));
#endif
_mm_store_sd(&C[(i*88)+6], c78_2);
__m128d c78_3 = _mm_loadu_pd(&C[(i*88)+14]);
__m128d a78_3 = _mm_loadu_pd(&values[1001]);
#if defined(__SSE3__) && defined(__AVX__)
c78_3 = _mm_add_pd(c78_3, _mm_mul_pd(a78_3, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c78_3 = _mm_add_pd(c78_3, _mm_mul_pd(a78_3, b78));
#endif
_mm_storeu_pd(&C[(i*88)+14], c78_3);
__m128d c78_5 = _mm_load_sd(&C[(i*88)+16]);
__m128d a78_5 = _mm_load_sd(&values[1003]);
#if defined(__SSE3__) && defined(__AVX__)
c78_5 = _mm_add_sd(c78_5, _mm_mul_sd(a78_5, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c78_5 = _mm_add_sd(c78_5, _mm_mul_sd(a78_5, b78));
#endif
_mm_store_sd(&C[(i*88)+16], c78_5);
__m128d c78_6 = _mm_loadu_pd(&C[(i*88)+29]);
__m128d a78_6 = _mm_loadu_pd(&values[1004]);
#if defined(__SSE3__) && defined(__AVX__)
c78_6 = _mm_add_pd(c78_6, _mm_mul_pd(a78_6, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c78_6 = _mm_add_pd(c78_6, _mm_mul_pd(a78_6, b78));
#endif
_mm_storeu_pd(&C[(i*88)+29], c78_6);
__m128d c78_8 = _mm_load_sd(&C[(i*88)+31]);
__m128d a78_8 = _mm_load_sd(&values[1006]);
#if defined(__SSE3__) && defined(__AVX__)
c78_8 = _mm_add_sd(c78_8, _mm_mul_sd(a78_8, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c78_8 = _mm_add_sd(c78_8, _mm_mul_sd(a78_8, b78));
#endif
_mm_store_sd(&C[(i*88)+31], c78_8);
__m128d c78_9 = _mm_loadu_pd(&C[(i*88)+50]);
__m128d a78_9 = _mm_loadu_pd(&values[1007]);
#if defined(__SSE3__) && defined(__AVX__)
c78_9 = _mm_add_pd(c78_9, _mm_mul_pd(a78_9, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c78_9 = _mm_add_pd(c78_9, _mm_mul_pd(a78_9, b78));
#endif
_mm_storeu_pd(&C[(i*88)+50], c78_9);
__m128d c78_11 = _mm_load_sd(&C[(i*88)+52]);
__m128d a78_11 = _mm_load_sd(&values[1009]);
#if defined(__SSE3__) && defined(__AVX__)
c78_11 = _mm_add_sd(c78_11, _mm_mul_sd(a78_11, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c78_11 = _mm_add_sd(c78_11, _mm_mul_sd(a78_11, b78));
#endif
_mm_store_sd(&C[(i*88)+52], c78_11);
__m128d c78_12 = _mm_loadu_pd(&C[(i*88)+78]);
__m128d a78_12 = _mm_loadu_pd(&values[1010]);
#if defined(__SSE3__) && defined(__AVX__)
c78_12 = _mm_add_pd(c78_12, _mm_mul_pd(a78_12, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c78_12 = _mm_add_pd(c78_12, _mm_mul_pd(a78_12, b78));
#endif
_mm_storeu_pd(&C[(i*88)+78], c78_12);
__m128d c78_14 = _mm_load_sd(&C[(i*88)+80]);
__m128d a78_14 = _mm_load_sd(&values[1012]);
#if defined(__SSE3__) && defined(__AVX__)
c78_14 = _mm_add_sd(c78_14, _mm_mul_sd(a78_14, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c78_14 = _mm_add_sd(c78_14, _mm_mul_sd(a78_14, b78));
#endif
_mm_store_sd(&C[(i*88)+80], c78_14);
#else
C[(i*88)+4] += values[998] * B[(i*88)+78];
C[(i*88)+5] += values[999] * B[(i*88)+78];
C[(i*88)+6] += values[1000] * B[(i*88)+78];
C[(i*88)+14] += values[1001] * B[(i*88)+78];
C[(i*88)+15] += values[1002] * B[(i*88)+78];
C[(i*88)+16] += values[1003] * B[(i*88)+78];
C[(i*88)+29] += values[1004] * B[(i*88)+78];
C[(i*88)+30] += values[1005] * B[(i*88)+78];
C[(i*88)+31] += values[1006] * B[(i*88)+78];
C[(i*88)+50] += values[1007] * B[(i*88)+78];
C[(i*88)+51] += values[1008] * B[(i*88)+78];
C[(i*88)+52] += values[1009] * B[(i*88)+78];
C[(i*88)+78] += values[1010] * B[(i*88)+78];
C[(i*88)+79] += values[1011] * B[(i*88)+78];
C[(i*88)+80] += values[1012] * B[(i*88)+78];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b79 = _mm256_broadcast_sd(&B[(i*88)+79]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b79 = _mm_loaddup_pd(&B[(i*88)+79]);
#endif
__m128d c79_0 = _mm_loadu_pd(&C[(i*88)+4]);
__m128d a79_0 = _mm_loadu_pd(&values[1013]);
#if defined(__SSE3__) && defined(__AVX__)
c79_0 = _mm_add_pd(c79_0, _mm_mul_pd(a79_0, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c79_0 = _mm_add_pd(c79_0, _mm_mul_pd(a79_0, b79));
#endif
_mm_storeu_pd(&C[(i*88)+4], c79_0);
__m128d c79_2 = _mm_load_sd(&C[(i*88)+6]);
__m128d a79_2 = _mm_load_sd(&values[1015]);
#if defined(__SSE3__) && defined(__AVX__)
c79_2 = _mm_add_sd(c79_2, _mm_mul_sd(a79_2, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c79_2 = _mm_add_sd(c79_2, _mm_mul_sd(a79_2, b79));
#endif
_mm_store_sd(&C[(i*88)+6], c79_2);
__m128d c79_3 = _mm_loadu_pd(&C[(i*88)+14]);
__m128d a79_3 = _mm_loadu_pd(&values[1016]);
#if defined(__SSE3__) && defined(__AVX__)
c79_3 = _mm_add_pd(c79_3, _mm_mul_pd(a79_3, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c79_3 = _mm_add_pd(c79_3, _mm_mul_pd(a79_3, b79));
#endif
_mm_storeu_pd(&C[(i*88)+14], c79_3);
__m128d c79_5 = _mm_load_sd(&C[(i*88)+16]);
__m128d a79_5 = _mm_load_sd(&values[1018]);
#if defined(__SSE3__) && defined(__AVX__)
c79_5 = _mm_add_sd(c79_5, _mm_mul_sd(a79_5, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c79_5 = _mm_add_sd(c79_5, _mm_mul_sd(a79_5, b79));
#endif
_mm_store_sd(&C[(i*88)+16], c79_5);
__m128d c79_6 = _mm_loadu_pd(&C[(i*88)+29]);
__m128d a79_6 = _mm_loadu_pd(&values[1019]);
#if defined(__SSE3__) && defined(__AVX__)
c79_6 = _mm_add_pd(c79_6, _mm_mul_pd(a79_6, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c79_6 = _mm_add_pd(c79_6, _mm_mul_pd(a79_6, b79));
#endif
_mm_storeu_pd(&C[(i*88)+29], c79_6);
__m128d c79_8 = _mm_load_sd(&C[(i*88)+31]);
__m128d a79_8 = _mm_load_sd(&values[1021]);
#if defined(__SSE3__) && defined(__AVX__)
c79_8 = _mm_add_sd(c79_8, _mm_mul_sd(a79_8, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c79_8 = _mm_add_sd(c79_8, _mm_mul_sd(a79_8, b79));
#endif
_mm_store_sd(&C[(i*88)+31], c79_8);
__m128d c79_9 = _mm_loadu_pd(&C[(i*88)+50]);
__m128d a79_9 = _mm_loadu_pd(&values[1022]);
#if defined(__SSE3__) && defined(__AVX__)
c79_9 = _mm_add_pd(c79_9, _mm_mul_pd(a79_9, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c79_9 = _mm_add_pd(c79_9, _mm_mul_pd(a79_9, b79));
#endif
_mm_storeu_pd(&C[(i*88)+50], c79_9);
__m128d c79_11 = _mm_load_sd(&C[(i*88)+52]);
__m128d a79_11 = _mm_load_sd(&values[1024]);
#if defined(__SSE3__) && defined(__AVX__)
c79_11 = _mm_add_sd(c79_11, _mm_mul_sd(a79_11, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c79_11 = _mm_add_sd(c79_11, _mm_mul_sd(a79_11, b79));
#endif
_mm_store_sd(&C[(i*88)+52], c79_11);
__m128d c79_12 = _mm_loadu_pd(&C[(i*88)+78]);
__m128d a79_12 = _mm_loadu_pd(&values[1025]);
#if defined(__SSE3__) && defined(__AVX__)
c79_12 = _mm_add_pd(c79_12, _mm_mul_pd(a79_12, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c79_12 = _mm_add_pd(c79_12, _mm_mul_pd(a79_12, b79));
#endif
_mm_storeu_pd(&C[(i*88)+78], c79_12);
__m128d c79_14 = _mm_load_sd(&C[(i*88)+80]);
__m128d a79_14 = _mm_load_sd(&values[1027]);
#if defined(__SSE3__) && defined(__AVX__)
c79_14 = _mm_add_sd(c79_14, _mm_mul_sd(a79_14, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c79_14 = _mm_add_sd(c79_14, _mm_mul_sd(a79_14, b79));
#endif
_mm_store_sd(&C[(i*88)+80], c79_14);
#else
C[(i*88)+4] += values[1013] * B[(i*88)+79];
C[(i*88)+5] += values[1014] * B[(i*88)+79];
C[(i*88)+6] += values[1015] * B[(i*88)+79];
C[(i*88)+14] += values[1016] * B[(i*88)+79];
C[(i*88)+15] += values[1017] * B[(i*88)+79];
C[(i*88)+16] += values[1018] * B[(i*88)+79];
C[(i*88)+29] += values[1019] * B[(i*88)+79];
C[(i*88)+30] += values[1020] * B[(i*88)+79];
C[(i*88)+31] += values[1021] * B[(i*88)+79];
C[(i*88)+50] += values[1022] * B[(i*88)+79];
C[(i*88)+51] += values[1023] * B[(i*88)+79];
C[(i*88)+52] += values[1024] * B[(i*88)+79];
C[(i*88)+78] += values[1025] * B[(i*88)+79];
C[(i*88)+79] += values[1026] * B[(i*88)+79];
C[(i*88)+80] += values[1027] * B[(i*88)+79];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b80 = _mm256_broadcast_sd(&B[(i*88)+80]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b80 = _mm_loaddup_pd(&B[(i*88)+80]);
#endif
__m128d c80_0 = _mm_loadu_pd(&C[(i*88)+4]);
__m128d a80_0 = _mm_loadu_pd(&values[1028]);
#if defined(__SSE3__) && defined(__AVX__)
c80_0 = _mm_add_pd(c80_0, _mm_mul_pd(a80_0, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c80_0 = _mm_add_pd(c80_0, _mm_mul_pd(a80_0, b80));
#endif
_mm_storeu_pd(&C[(i*88)+4], c80_0);
__m128d c80_2 = _mm_load_sd(&C[(i*88)+6]);
__m128d a80_2 = _mm_load_sd(&values[1030]);
#if defined(__SSE3__) && defined(__AVX__)
c80_2 = _mm_add_sd(c80_2, _mm_mul_sd(a80_2, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c80_2 = _mm_add_sd(c80_2, _mm_mul_sd(a80_2, b80));
#endif
_mm_store_sd(&C[(i*88)+6], c80_2);
__m128d c80_3 = _mm_loadu_pd(&C[(i*88)+14]);
__m128d a80_3 = _mm_loadu_pd(&values[1031]);
#if defined(__SSE3__) && defined(__AVX__)
c80_3 = _mm_add_pd(c80_3, _mm_mul_pd(a80_3, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c80_3 = _mm_add_pd(c80_3, _mm_mul_pd(a80_3, b80));
#endif
_mm_storeu_pd(&C[(i*88)+14], c80_3);
__m128d c80_5 = _mm_load_sd(&C[(i*88)+16]);
__m128d a80_5 = _mm_load_sd(&values[1033]);
#if defined(__SSE3__) && defined(__AVX__)
c80_5 = _mm_add_sd(c80_5, _mm_mul_sd(a80_5, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c80_5 = _mm_add_sd(c80_5, _mm_mul_sd(a80_5, b80));
#endif
_mm_store_sd(&C[(i*88)+16], c80_5);
__m128d c80_6 = _mm_loadu_pd(&C[(i*88)+29]);
__m128d a80_6 = _mm_loadu_pd(&values[1034]);
#if defined(__SSE3__) && defined(__AVX__)
c80_6 = _mm_add_pd(c80_6, _mm_mul_pd(a80_6, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c80_6 = _mm_add_pd(c80_6, _mm_mul_pd(a80_6, b80));
#endif
_mm_storeu_pd(&C[(i*88)+29], c80_6);
__m128d c80_8 = _mm_load_sd(&C[(i*88)+31]);
__m128d a80_8 = _mm_load_sd(&values[1036]);
#if defined(__SSE3__) && defined(__AVX__)
c80_8 = _mm_add_sd(c80_8, _mm_mul_sd(a80_8, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c80_8 = _mm_add_sd(c80_8, _mm_mul_sd(a80_8, b80));
#endif
_mm_store_sd(&C[(i*88)+31], c80_8);
__m128d c80_9 = _mm_loadu_pd(&C[(i*88)+50]);
__m128d a80_9 = _mm_loadu_pd(&values[1037]);
#if defined(__SSE3__) && defined(__AVX__)
c80_9 = _mm_add_pd(c80_9, _mm_mul_pd(a80_9, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c80_9 = _mm_add_pd(c80_9, _mm_mul_pd(a80_9, b80));
#endif
_mm_storeu_pd(&C[(i*88)+50], c80_9);
__m128d c80_11 = _mm_load_sd(&C[(i*88)+52]);
__m128d a80_11 = _mm_load_sd(&values[1039]);
#if defined(__SSE3__) && defined(__AVX__)
c80_11 = _mm_add_sd(c80_11, _mm_mul_sd(a80_11, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c80_11 = _mm_add_sd(c80_11, _mm_mul_sd(a80_11, b80));
#endif
_mm_store_sd(&C[(i*88)+52], c80_11);
__m128d c80_12 = _mm_loadu_pd(&C[(i*88)+78]);
__m128d a80_12 = _mm_loadu_pd(&values[1040]);
#if defined(__SSE3__) && defined(__AVX__)
c80_12 = _mm_add_pd(c80_12, _mm_mul_pd(a80_12, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c80_12 = _mm_add_pd(c80_12, _mm_mul_pd(a80_12, b80));
#endif
_mm_storeu_pd(&C[(i*88)+78], c80_12);
__m128d c80_14 = _mm_load_sd(&C[(i*88)+80]);
__m128d a80_14 = _mm_load_sd(&values[1042]);
#if defined(__SSE3__) && defined(__AVX__)
c80_14 = _mm_add_sd(c80_14, _mm_mul_sd(a80_14, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c80_14 = _mm_add_sd(c80_14, _mm_mul_sd(a80_14, b80));
#endif
_mm_store_sd(&C[(i*88)+80], c80_14);
#else
C[(i*88)+4] += values[1028] * B[(i*88)+80];
C[(i*88)+5] += values[1029] * B[(i*88)+80];
C[(i*88)+6] += values[1030] * B[(i*88)+80];
C[(i*88)+14] += values[1031] * B[(i*88)+80];
C[(i*88)+15] += values[1032] * B[(i*88)+80];
C[(i*88)+16] += values[1033] * B[(i*88)+80];
C[(i*88)+29] += values[1034] * B[(i*88)+80];
C[(i*88)+30] += values[1035] * B[(i*88)+80];
C[(i*88)+31] += values[1036] * B[(i*88)+80];
C[(i*88)+50] += values[1037] * B[(i*88)+80];
C[(i*88)+51] += values[1038] * B[(i*88)+80];
C[(i*88)+52] += values[1039] * B[(i*88)+80];
C[(i*88)+78] += values[1040] * B[(i*88)+80];
C[(i*88)+79] += values[1041] * B[(i*88)+80];
C[(i*88)+80] += values[1042] * B[(i*88)+80];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b81 = _mm256_broadcast_sd(&B[(i*88)+81]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b81 = _mm_loaddup_pd(&B[(i*88)+81]);
#endif
__m128d c81_0 = _mm_loadu_pd(&C[(i*88)+1]);
__m128d a81_0 = _mm_loadu_pd(&values[1043]);
#if defined(__SSE3__) && defined(__AVX__)
c81_0 = _mm_add_pd(c81_0, _mm_mul_pd(a81_0, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c81_0 = _mm_add_pd(c81_0, _mm_mul_pd(a81_0, b81));
#endif
_mm_storeu_pd(&C[(i*88)+1], c81_0);
__m128d c81_2 = _mm_loadu_pd(&C[(i*88)+7]);
__m128d a81_2 = _mm_loadu_pd(&values[1045]);
#if defined(__SSE3__) && defined(__AVX__)
c81_2 = _mm_add_pd(c81_2, _mm_mul_pd(a81_2, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c81_2 = _mm_add_pd(c81_2, _mm_mul_pd(a81_2, b81));
#endif
_mm_storeu_pd(&C[(i*88)+7], c81_2);
__m128d c81_4 = _mm_loadu_pd(&C[(i*88)+17]);
__m128d a81_4 = _mm_loadu_pd(&values[1047]);
#if defined(__SSE3__) && defined(__AVX__)
c81_4 = _mm_add_pd(c81_4, _mm_mul_pd(a81_4, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c81_4 = _mm_add_pd(c81_4, _mm_mul_pd(a81_4, b81));
#endif
_mm_storeu_pd(&C[(i*88)+17], c81_4);
__m128d c81_6 = _mm_loadu_pd(&C[(i*88)+32]);
__m128d a81_6 = _mm_loadu_pd(&values[1049]);
#if defined(__SSE3__) && defined(__AVX__)
c81_6 = _mm_add_pd(c81_6, _mm_mul_pd(a81_6, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c81_6 = _mm_add_pd(c81_6, _mm_mul_pd(a81_6, b81));
#endif
_mm_storeu_pd(&C[(i*88)+32], c81_6);
__m128d c81_8 = _mm_loadu_pd(&C[(i*88)+53]);
__m128d a81_8 = _mm_loadu_pd(&values[1051]);
#if defined(__SSE3__) && defined(__AVX__)
c81_8 = _mm_add_pd(c81_8, _mm_mul_pd(a81_8, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c81_8 = _mm_add_pd(c81_8, _mm_mul_pd(a81_8, b81));
#endif
_mm_storeu_pd(&C[(i*88)+53], c81_8);
__m128d c81_10 = _mm_loadu_pd(&C[(i*88)+81]);
__m128d a81_10 = _mm_loadu_pd(&values[1053]);
#if defined(__SSE3__) && defined(__AVX__)
c81_10 = _mm_add_pd(c81_10, _mm_mul_pd(a81_10, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c81_10 = _mm_add_pd(c81_10, _mm_mul_pd(a81_10, b81));
#endif
_mm_storeu_pd(&C[(i*88)+81], c81_10);
#else
C[(i*88)+1] += values[1043] * B[(i*88)+81];
C[(i*88)+2] += values[1044] * B[(i*88)+81];
C[(i*88)+7] += values[1045] * B[(i*88)+81];
C[(i*88)+8] += values[1046] * B[(i*88)+81];
C[(i*88)+17] += values[1047] * B[(i*88)+81];
C[(i*88)+18] += values[1048] * B[(i*88)+81];
C[(i*88)+32] += values[1049] * B[(i*88)+81];
C[(i*88)+33] += values[1050] * B[(i*88)+81];
C[(i*88)+53] += values[1051] * B[(i*88)+81];
C[(i*88)+54] += values[1052] * B[(i*88)+81];
C[(i*88)+81] += values[1053] * B[(i*88)+81];
C[(i*88)+82] += values[1054] * B[(i*88)+81];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b82 = _mm256_broadcast_sd(&B[(i*88)+82]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b82 = _mm_loaddup_pd(&B[(i*88)+82]);
#endif
__m128d c82_0 = _mm_loadu_pd(&C[(i*88)+1]);
__m128d a82_0 = _mm_loadu_pd(&values[1055]);
#if defined(__SSE3__) && defined(__AVX__)
c82_0 = _mm_add_pd(c82_0, _mm_mul_pd(a82_0, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c82_0 = _mm_add_pd(c82_0, _mm_mul_pd(a82_0, b82));
#endif
_mm_storeu_pd(&C[(i*88)+1], c82_0);
__m128d c82_2 = _mm_loadu_pd(&C[(i*88)+7]);
__m128d a82_2 = _mm_loadu_pd(&values[1057]);
#if defined(__SSE3__) && defined(__AVX__)
c82_2 = _mm_add_pd(c82_2, _mm_mul_pd(a82_2, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c82_2 = _mm_add_pd(c82_2, _mm_mul_pd(a82_2, b82));
#endif
_mm_storeu_pd(&C[(i*88)+7], c82_2);
__m128d c82_4 = _mm_loadu_pd(&C[(i*88)+17]);
__m128d a82_4 = _mm_loadu_pd(&values[1059]);
#if defined(__SSE3__) && defined(__AVX__)
c82_4 = _mm_add_pd(c82_4, _mm_mul_pd(a82_4, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c82_4 = _mm_add_pd(c82_4, _mm_mul_pd(a82_4, b82));
#endif
_mm_storeu_pd(&C[(i*88)+17], c82_4);
__m128d c82_6 = _mm_loadu_pd(&C[(i*88)+32]);
__m128d a82_6 = _mm_loadu_pd(&values[1061]);
#if defined(__SSE3__) && defined(__AVX__)
c82_6 = _mm_add_pd(c82_6, _mm_mul_pd(a82_6, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c82_6 = _mm_add_pd(c82_6, _mm_mul_pd(a82_6, b82));
#endif
_mm_storeu_pd(&C[(i*88)+32], c82_6);
__m128d c82_8 = _mm_loadu_pd(&C[(i*88)+53]);
__m128d a82_8 = _mm_loadu_pd(&values[1063]);
#if defined(__SSE3__) && defined(__AVX__)
c82_8 = _mm_add_pd(c82_8, _mm_mul_pd(a82_8, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c82_8 = _mm_add_pd(c82_8, _mm_mul_pd(a82_8, b82));
#endif
_mm_storeu_pd(&C[(i*88)+53], c82_8);
__m128d c82_10 = _mm_loadu_pd(&C[(i*88)+81]);
__m128d a82_10 = _mm_loadu_pd(&values[1065]);
#if defined(__SSE3__) && defined(__AVX__)
c82_10 = _mm_add_pd(c82_10, _mm_mul_pd(a82_10, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c82_10 = _mm_add_pd(c82_10, _mm_mul_pd(a82_10, b82));
#endif
_mm_storeu_pd(&C[(i*88)+81], c82_10);
#else
C[(i*88)+1] += values[1055] * B[(i*88)+82];
C[(i*88)+2] += values[1056] * B[(i*88)+82];
C[(i*88)+7] += values[1057] * B[(i*88)+82];
C[(i*88)+8] += values[1058] * B[(i*88)+82];
C[(i*88)+17] += values[1059] * B[(i*88)+82];
C[(i*88)+18] += values[1060] * B[(i*88)+82];
C[(i*88)+32] += values[1061] * B[(i*88)+82];
C[(i*88)+33] += values[1062] * B[(i*88)+82];
C[(i*88)+53] += values[1063] * B[(i*88)+82];
C[(i*88)+54] += values[1064] * B[(i*88)+82];
C[(i*88)+81] += values[1065] * B[(i*88)+82];
C[(i*88)+82] += values[1066] * B[(i*88)+82];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b83 = _mm256_broadcast_sd(&B[(i*88)+83]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b83 = _mm_loaddup_pd(&B[(i*88)+83]);
#endif
__m128d c83_0 = _mm_load_sd(&C[(i*88)+0]);
__m128d a83_0 = _mm_load_sd(&values[1067]);
#if defined(__SSE3__) && defined(__AVX__)
c83_0 = _mm_add_sd(c83_0, _mm_mul_sd(a83_0, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c83_0 = _mm_add_sd(c83_0, _mm_mul_sd(a83_0, b83));
#endif
_mm_store_sd(&C[(i*88)+0], c83_0);
__m128d c83_1 = _mm_load_sd(&C[(i*88)+3]);
__m128d a83_1 = _mm_load_sd(&values[1068]);
#if defined(__SSE3__) && defined(__AVX__)
c83_1 = _mm_add_sd(c83_1, _mm_mul_sd(a83_1, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c83_1 = _mm_add_sd(c83_1, _mm_mul_sd(a83_1, b83));
#endif
_mm_store_sd(&C[(i*88)+3], c83_1);
__m128d c83_2 = _mm_load_sd(&C[(i*88)+9]);
__m128d a83_2 = _mm_load_sd(&values[1069]);
#if defined(__SSE3__) && defined(__AVX__)
c83_2 = _mm_add_sd(c83_2, _mm_mul_sd(a83_2, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c83_2 = _mm_add_sd(c83_2, _mm_mul_sd(a83_2, b83));
#endif
_mm_store_sd(&C[(i*88)+9], c83_2);
__m128d c83_3 = _mm_load_sd(&C[(i*88)+19]);
__m128d a83_3 = _mm_load_sd(&values[1070]);
#if defined(__SSE3__) && defined(__AVX__)
c83_3 = _mm_add_sd(c83_3, _mm_mul_sd(a83_3, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c83_3 = _mm_add_sd(c83_3, _mm_mul_sd(a83_3, b83));
#endif
_mm_store_sd(&C[(i*88)+19], c83_3);
__m128d c83_4 = _mm_load_sd(&C[(i*88)+34]);
__m128d a83_4 = _mm_load_sd(&values[1071]);
#if defined(__SSE3__) && defined(__AVX__)
c83_4 = _mm_add_sd(c83_4, _mm_mul_sd(a83_4, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c83_4 = _mm_add_sd(c83_4, _mm_mul_sd(a83_4, b83));
#endif
_mm_store_sd(&C[(i*88)+34], c83_4);
__m128d c83_5 = _mm_load_sd(&C[(i*88)+55]);
__m128d a83_5 = _mm_load_sd(&values[1072]);
#if defined(__SSE3__) && defined(__AVX__)
c83_5 = _mm_add_sd(c83_5, _mm_mul_sd(a83_5, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c83_5 = _mm_add_sd(c83_5, _mm_mul_sd(a83_5, b83));
#endif
_mm_store_sd(&C[(i*88)+55], c83_5);
__m128d c83_6 = _mm_load_sd(&C[(i*88)+83]);
__m128d a83_6 = _mm_load_sd(&values[1073]);
#if defined(__SSE3__) && defined(__AVX__)
c83_6 = _mm_add_sd(c83_6, _mm_mul_sd(a83_6, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c83_6 = _mm_add_sd(c83_6, _mm_mul_sd(a83_6, b83));
#endif
_mm_store_sd(&C[(i*88)+83], c83_6);
#else
C[(i*88)+0] += values[1067] * B[(i*88)+83];
C[(i*88)+3] += values[1068] * B[(i*88)+83];
C[(i*88)+9] += values[1069] * B[(i*88)+83];
C[(i*88)+19] += values[1070] * B[(i*88)+83];
C[(i*88)+34] += values[1071] * B[(i*88)+83];
C[(i*88)+55] += values[1072] * B[(i*88)+83];
C[(i*88)+83] += values[1073] * B[(i*88)+83];
#endif

}

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 19332;
#endif

}

#endif
