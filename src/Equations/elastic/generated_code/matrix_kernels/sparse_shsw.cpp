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
// @date 2015-05-09 22:17:48.633330
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
#ifndef SPARSESHSWCPP
#define SPARSESHSWCPP

#if defined( __SSE3__) || defined(__MIC__)
#include <immintrin.h>
#endif

#include <cstddef>
#ifndef NDEBUG
extern long long libxsmm_num_total_flops;
#endif

void ssparse_starMatrix_m1_n9_k9_ldA8_ldBna2_ldC8_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
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

void ssparse_starMatrix_m4_n9_k9_ldA8_ldBna3_ldC8_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma simd vectorlength(4)
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

void ssparse_starMatrix_m1_n9_k9_ldA8_ldBna3_ldC8_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
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

void ssparse_starMatrix_m10_n9_k9_ldA16_ldBna4_ldC16_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
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

void ssparse_starMatrix_m4_n9_k9_ldA8_ldBna4_ldC8_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma simd vectorlength(4)
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

void ssparse_starMatrix_m1_n9_k9_ldA8_ldBna4_ldC8_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
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

void ssparse_starMatrix_m20_n9_k9_ldA24_ldBna5_ldC24_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
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

void ssparse_starMatrix_m10_n9_k9_ldA16_ldBna5_ldC16_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
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

void ssparse_starMatrix_m4_n9_k9_ldA8_ldBna5_ldC8_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma simd vectorlength(4)
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

void ssparse_starMatrix_m1_n9_k9_ldA8_ldBna5_ldC8_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
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

void ssparse_starMatrix_m35_n9_k9_ldA40_ldBna6_ldC40_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
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

void ssparse_starMatrix_m20_n9_k9_ldA24_ldBna6_ldC24_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
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

void ssparse_starMatrix_m10_n9_k9_ldA16_ldBna6_ldC16_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
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

void ssparse_starMatrix_m4_n9_k9_ldA8_ldBna6_ldC8_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma simd vectorlength(4)
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

void ssparse_starMatrix_m1_n9_k9_ldA8_ldBna6_ldC8_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
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

void ssparse_starMatrix_m56_n9_k9_ldA56_ldBna7_ldC56_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
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

void ssparse_starMatrix_m35_n9_k9_ldA40_ldBna7_ldC40_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
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

void ssparse_starMatrix_m20_n9_k9_ldA24_ldBna7_ldC24_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
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

void ssparse_starMatrix_m10_n9_k9_ldA16_ldBna7_ldC16_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
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

void ssparse_starMatrix_m4_n9_k9_ldA8_ldBna7_ldC8_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma simd vectorlength(4)
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

void ssparse_starMatrix_m1_n9_k9_ldA8_ldBna7_ldC8_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
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

void ssparse_starMatrix_m84_n9_k9_ldA88_ldBna8_ldC88_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
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

void ssparse_starMatrix_m56_n9_k9_ldA56_ldBna8_ldC56_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
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

void ssparse_starMatrix_m35_n9_k9_ldA40_ldBna8_ldC40_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
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

void ssparse_starMatrix_m20_n9_k9_ldA24_ldBna8_ldC24_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
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

void ssparse_starMatrix_m10_n9_k9_ldA16_ldBna8_ldC16_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
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

void ssparse_starMatrix_m4_n9_k9_ldA8_ldBna8_ldC8_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma simd vectorlength(4)
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

void ssparse_starMatrix_m1_n9_k9_ldA8_ldBna8_ldC8_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
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

void ssparse_starMatrix_m4_n9_k9_ldA8_ldBna2_ldC8_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma simd vectorlength(4)
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

void ssparse_starMatrix_m10_n9_k9_ldA16_ldBna3_ldC16_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
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

void ssparse_fP113DivM_m10_n9_k10_ldAna3_ldB16_ldC16_beta0_pfsigonly(const float* values, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma nounroll_and_jam
for (unsigned int i = 0; i < 9; i++)
{
  #pragma simd
  for (unsigned int m = 0; m < 10; m++) {
    C[(i*16)+m] = 0.0;
  }
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b0 = _mm_broadcast_ss(&B[(i*16)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b0 = _mm_load_ss(&B[(i*16)+0]);
b0 = _mm_shuffle_ps(b0, b0, 0x00);
#endif
__m128 c0_0 = _mm_load_ss(&C[(i*16)+0]);
__m128 a0_0 = _mm_load_ss(&values[0]);
c0_0 = _mm_add_ss(c0_0, _mm_mul_ss(a0_0, b0));
_mm_store_ss(&C[(i*16)+0], c0_0);
__m128 c0_1 = _mm_load_ss(&C[(i*16)+3]);
__m128 a0_1 = _mm_load_ss(&values[1]);
c0_1 = _mm_add_ss(c0_1, _mm_mul_ss(a0_1, b0));
_mm_store_ss(&C[(i*16)+3], c0_1);
__m128 c0_2 = _mm_load_ss(&C[(i*16)+9]);
__m128 a0_2 = _mm_load_ss(&values[2]);
c0_2 = _mm_add_ss(c0_2, _mm_mul_ss(a0_2, b0));
_mm_store_ss(&C[(i*16)+9], c0_2);
#else
C[(i*16)+0] += values[0] * B[(i*16)+0];
C[(i*16)+3] += values[1] * B[(i*16)+0];
C[(i*16)+9] += values[2] * B[(i*16)+0];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b1 = _mm_broadcast_ss(&B[(i*16)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b1 = _mm_load_ss(&B[(i*16)+1]);
b1 = _mm_shuffle_ps(b1, b1, 0x00);
#endif
__m128 c1_0 = _mm_load_ss(&C[(i*16)+1]);
__m128 a1_0 = _mm_load_ss(&values[3]);
c1_0 = _mm_add_ss(c1_0, _mm_mul_ss(a1_0, b1));
_mm_store_ss(&C[(i*16)+1], c1_0);
__m128 c1_1 = _mm_load_ss(&C[(i*16)+7]);
__m128 a1_1 = _mm_load_ss(&values[4]);
c1_1 = _mm_add_ss(c1_1, _mm_mul_ss(a1_1, b1));
_mm_store_ss(&C[(i*16)+7], c1_1);
#else
C[(i*16)+1] += values[3] * B[(i*16)+1];
C[(i*16)+7] += values[4] * B[(i*16)+1];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b2 = _mm_broadcast_ss(&B[(i*16)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b2 = _mm_load_ss(&B[(i*16)+2]);
b2 = _mm_shuffle_ps(b2, b2, 0x00);
#endif
__m128 c2_0 = _mm_load_ss(&C[(i*16)+2]);
__m128 a2_0 = _mm_load_ss(&values[5]);
c2_0 = _mm_add_ss(c2_0, _mm_mul_ss(a2_0, b2));
_mm_store_ss(&C[(i*16)+2], c2_0);
__m128 c2_1 = _mm_load_ss(&C[(i*16)+8]);
__m128 a2_1 = _mm_load_ss(&values[6]);
c2_1 = _mm_add_ss(c2_1, _mm_mul_ss(a2_1, b2));
_mm_store_ss(&C[(i*16)+8], c2_1);
#else
C[(i*16)+2] += values[5] * B[(i*16)+2];
C[(i*16)+8] += values[6] * B[(i*16)+2];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b3 = _mm_broadcast_ss(&B[(i*16)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b3 = _mm_load_ss(&B[(i*16)+3]);
b3 = _mm_shuffle_ps(b3, b3, 0x00);
#endif
__m128 c3_0 = _mm_load_ss(&C[(i*16)+0]);
__m128 a3_0 = _mm_load_ss(&values[7]);
c3_0 = _mm_add_ss(c3_0, _mm_mul_ss(a3_0, b3));
_mm_store_ss(&C[(i*16)+0], c3_0);
__m128 c3_1 = _mm_load_ss(&C[(i*16)+3]);
__m128 a3_1 = _mm_load_ss(&values[8]);
c3_1 = _mm_add_ss(c3_1, _mm_mul_ss(a3_1, b3));
_mm_store_ss(&C[(i*16)+3], c3_1);
__m128 c3_2 = _mm_load_ss(&C[(i*16)+9]);
__m128 a3_2 = _mm_load_ss(&values[9]);
c3_2 = _mm_add_ss(c3_2, _mm_mul_ss(a3_2, b3));
_mm_store_ss(&C[(i*16)+9], c3_2);
#else
C[(i*16)+0] += values[7] * B[(i*16)+3];
C[(i*16)+3] += values[8] * B[(i*16)+3];
C[(i*16)+9] += values[9] * B[(i*16)+3];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b4 = _mm_broadcast_ss(&B[(i*16)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b4 = _mm_load_ss(&B[(i*16)+4]);
b4 = _mm_shuffle_ps(b4, b4, 0x00);
#endif
__m128 c4_0 = _mm_load_ss(&C[(i*16)+4]);
__m128 a4_0 = _mm_load_ss(&values[10]);
c4_0 = _mm_add_ss(c4_0, _mm_mul_ss(a4_0, b4));
_mm_store_ss(&C[(i*16)+4], c4_0);
#else
C[(i*16)+4] += values[10] * B[(i*16)+4];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b5 = _mm_broadcast_ss(&B[(i*16)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b5 = _mm_load_ss(&B[(i*16)+5]);
b5 = _mm_shuffle_ps(b5, b5, 0x00);
#endif
__m128 c5_0 = _mm_load_ss(&C[(i*16)+5]);
__m128 a5_0 = _mm_load_ss(&values[11]);
c5_0 = _mm_add_ss(c5_0, _mm_mul_ss(a5_0, b5));
_mm_store_ss(&C[(i*16)+5], c5_0);
#else
C[(i*16)+5] += values[11] * B[(i*16)+5];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b6 = _mm_broadcast_ss(&B[(i*16)+6]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b6 = _mm_load_ss(&B[(i*16)+6]);
b6 = _mm_shuffle_ps(b6, b6, 0x00);
#endif
__m128 c6_0 = _mm_load_ss(&C[(i*16)+6]);
__m128 a6_0 = _mm_load_ss(&values[12]);
c6_0 = _mm_add_ss(c6_0, _mm_mul_ss(a6_0, b6));
_mm_store_ss(&C[(i*16)+6], c6_0);
#else
C[(i*16)+6] += values[12] * B[(i*16)+6];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b7 = _mm_broadcast_ss(&B[(i*16)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b7 = _mm_load_ss(&B[(i*16)+7]);
b7 = _mm_shuffle_ps(b7, b7, 0x00);
#endif
__m128 c7_0 = _mm_load_ss(&C[(i*16)+1]);
__m128 a7_0 = _mm_load_ss(&values[13]);
c7_0 = _mm_add_ss(c7_0, _mm_mul_ss(a7_0, b7));
_mm_store_ss(&C[(i*16)+1], c7_0);
__m128 c7_1 = _mm_load_ss(&C[(i*16)+7]);
__m128 a7_1 = _mm_load_ss(&values[14]);
c7_1 = _mm_add_ss(c7_1, _mm_mul_ss(a7_1, b7));
_mm_store_ss(&C[(i*16)+7], c7_1);
#else
C[(i*16)+1] += values[13] * B[(i*16)+7];
C[(i*16)+7] += values[14] * B[(i*16)+7];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b8 = _mm_broadcast_ss(&B[(i*16)+8]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b8 = _mm_load_ss(&B[(i*16)+8]);
b8 = _mm_shuffle_ps(b8, b8, 0x00);
#endif
__m128 c8_0 = _mm_load_ss(&C[(i*16)+2]);
__m128 a8_0 = _mm_load_ss(&values[15]);
c8_0 = _mm_add_ss(c8_0, _mm_mul_ss(a8_0, b8));
_mm_store_ss(&C[(i*16)+2], c8_0);
__m128 c8_1 = _mm_load_ss(&C[(i*16)+8]);
__m128 a8_1 = _mm_load_ss(&values[16]);
c8_1 = _mm_add_ss(c8_1, _mm_mul_ss(a8_1, b8));
_mm_store_ss(&C[(i*16)+8], c8_1);
#else
C[(i*16)+2] += values[15] * B[(i*16)+8];
C[(i*16)+8] += values[16] * B[(i*16)+8];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b9 = _mm_broadcast_ss(&B[(i*16)+9]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b9 = _mm_load_ss(&B[(i*16)+9]);
b9 = _mm_shuffle_ps(b9, b9, 0x00);
#endif
__m128 c9_0 = _mm_load_ss(&C[(i*16)+0]);
__m128 a9_0 = _mm_load_ss(&values[17]);
c9_0 = _mm_add_ss(c9_0, _mm_mul_ss(a9_0, b9));
_mm_store_ss(&C[(i*16)+0], c9_0);
__m128 c9_1 = _mm_load_ss(&C[(i*16)+3]);
__m128 a9_1 = _mm_load_ss(&values[18]);
c9_1 = _mm_add_ss(c9_1, _mm_mul_ss(a9_1, b9));
_mm_store_ss(&C[(i*16)+3], c9_1);
__m128 c9_2 = _mm_load_ss(&C[(i*16)+9]);
__m128 a9_2 = _mm_load_ss(&values[19]);
c9_2 = _mm_add_ss(c9_2, _mm_mul_ss(a9_2, b9));
_mm_store_ss(&C[(i*16)+9], c9_2);
#else
C[(i*16)+0] += values[17] * B[(i*16)+9];
C[(i*16)+3] += values[18] * B[(i*16)+9];
C[(i*16)+9] += values[19] * B[(i*16)+9];
#endif

}

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 360;
#endif

}

void ssparse_fP111DivM_m10_n9_k10_ldAna3_ldB16_ldC16_beta0_pfsigonly(const float* values, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma nounroll_and_jam
for (unsigned int i = 0; i < 9; i++)
{
  #pragma simd
  for (unsigned int m = 0; m < 10; m++) {
    C[(i*16)+m] = 0.0;
  }
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b0 = _mm_broadcast_ss(&B[(i*16)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b0 = _mm_load_ss(&B[(i*16)+0]);
b0 = _mm_shuffle_ps(b0, b0, 0x00);
#endif
__m128 c0_0 = _mm_load_ss(&C[(i*16)+0]);
__m128 a0_0 = _mm_load_ss(&values[0]);
c0_0 = _mm_add_ss(c0_0, _mm_mul_ss(a0_0, b0));
_mm_store_ss(&C[(i*16)+0], c0_0);
__m128 c0_1 = _mm_load_ss(&C[(i*16)+3]);
__m128 a0_1 = _mm_load_ss(&values[1]);
c0_1 = _mm_add_ss(c0_1, _mm_mul_ss(a0_1, b0));
_mm_store_ss(&C[(i*16)+3], c0_1);
__m128 c0_2 = _mm_load_ss(&C[(i*16)+9]);
__m128 a0_2 = _mm_load_ss(&values[2]);
c0_2 = _mm_add_ss(c0_2, _mm_mul_ss(a0_2, b0));
_mm_store_ss(&C[(i*16)+9], c0_2);
#else
C[(i*16)+0] += values[0] * B[(i*16)+0];
C[(i*16)+3] += values[1] * B[(i*16)+0];
C[(i*16)+9] += values[2] * B[(i*16)+0];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b1 = _mm_broadcast_ss(&B[(i*16)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b1 = _mm_load_ss(&B[(i*16)+1]);
b1 = _mm_shuffle_ps(b1, b1, 0x00);
#endif
__m128 c1_0 = _mm_castpd_ps(_mm_load_sd((double*)(&C[(i*16)+1])));
__m128 a1_0 = _mm_castpd_ps(_mm_load_sd((double*)(&values[3])));
c1_0 = _mm_add_ps(c1_0, _mm_mul_ps(a1_0, b1));
_mm_store_sd((double*)(&C[(i*16)+1]), _mm_castps_pd(c1_0));
__m128 c1_2 = _mm_castpd_ps(_mm_load_sd((double*)(&C[(i*16)+7])));
__m128 a1_2 = _mm_castpd_ps(_mm_load_sd((double*)(&values[5])));
c1_2 = _mm_add_ps(c1_2, _mm_mul_ps(a1_2, b1));
_mm_store_sd((double*)(&C[(i*16)+7]), _mm_castps_pd(c1_2));
#else
C[(i*16)+1] += values[3] * B[(i*16)+1];
C[(i*16)+2] += values[4] * B[(i*16)+1];
C[(i*16)+7] += values[5] * B[(i*16)+1];
C[(i*16)+8] += values[6] * B[(i*16)+1];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b2 = _mm_broadcast_ss(&B[(i*16)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b2 = _mm_load_ss(&B[(i*16)+2]);
b2 = _mm_shuffle_ps(b2, b2, 0x00);
#endif
__m128 c2_0 = _mm_castpd_ps(_mm_load_sd((double*)(&C[(i*16)+1])));
__m128 a2_0 = _mm_castpd_ps(_mm_load_sd((double*)(&values[7])));
c2_0 = _mm_add_ps(c2_0, _mm_mul_ps(a2_0, b2));
_mm_store_sd((double*)(&C[(i*16)+1]), _mm_castps_pd(c2_0));
__m128 c2_2 = _mm_castpd_ps(_mm_load_sd((double*)(&C[(i*16)+7])));
__m128 a2_2 = _mm_castpd_ps(_mm_load_sd((double*)(&values[9])));
c2_2 = _mm_add_ps(c2_2, _mm_mul_ps(a2_2, b2));
_mm_store_sd((double*)(&C[(i*16)+7]), _mm_castps_pd(c2_2));
#else
C[(i*16)+1] += values[7] * B[(i*16)+2];
C[(i*16)+2] += values[8] * B[(i*16)+2];
C[(i*16)+7] += values[9] * B[(i*16)+2];
C[(i*16)+8] += values[10] * B[(i*16)+2];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b3 = _mm_broadcast_ss(&B[(i*16)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b3 = _mm_load_ss(&B[(i*16)+3]);
b3 = _mm_shuffle_ps(b3, b3, 0x00);
#endif
__m128 c3_0 = _mm_load_ss(&C[(i*16)+0]);
__m128 a3_0 = _mm_load_ss(&values[11]);
c3_0 = _mm_add_ss(c3_0, _mm_mul_ss(a3_0, b3));
_mm_store_ss(&C[(i*16)+0], c3_0);
__m128 c3_1 = _mm_load_ss(&C[(i*16)+3]);
__m128 a3_1 = _mm_load_ss(&values[12]);
c3_1 = _mm_add_ss(c3_1, _mm_mul_ss(a3_1, b3));
_mm_store_ss(&C[(i*16)+3], c3_1);
__m128 c3_2 = _mm_load_ss(&C[(i*16)+9]);
__m128 a3_2 = _mm_load_ss(&values[13]);
c3_2 = _mm_add_ss(c3_2, _mm_mul_ss(a3_2, b3));
_mm_store_ss(&C[(i*16)+9], c3_2);
#else
C[(i*16)+0] += values[11] * B[(i*16)+3];
C[(i*16)+3] += values[12] * B[(i*16)+3];
C[(i*16)+9] += values[13] * B[(i*16)+3];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b4 = _mm_broadcast_ss(&B[(i*16)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b4 = _mm_load_ss(&B[(i*16)+4]);
b4 = _mm_shuffle_ps(b4, b4, 0x00);
#endif
__m128 c4_0 = _mm_castpd_ps(_mm_load_sd((double*)(&C[(i*16)+4])));
__m128 a4_0 = _mm_castpd_ps(_mm_load_sd((double*)(&values[14])));
c4_0 = _mm_add_ps(c4_0, _mm_mul_ps(a4_0, b4));
_mm_store_sd((double*)(&C[(i*16)+4]), _mm_castps_pd(c4_0));
__m128 c4_2 = _mm_load_ss(&C[(i*16)+6]);
__m128 a4_2 = _mm_load_ss(&values[16]);
c4_2 = _mm_add_ss(c4_2, _mm_mul_ss(a4_2, b4));
_mm_store_ss(&C[(i*16)+6], c4_2);
#else
C[(i*16)+4] += values[14] * B[(i*16)+4];
C[(i*16)+5] += values[15] * B[(i*16)+4];
C[(i*16)+6] += values[16] * B[(i*16)+4];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b5 = _mm_broadcast_ss(&B[(i*16)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b5 = _mm_load_ss(&B[(i*16)+5]);
b5 = _mm_shuffle_ps(b5, b5, 0x00);
#endif
__m128 c5_0 = _mm_castpd_ps(_mm_load_sd((double*)(&C[(i*16)+4])));
__m128 a5_0 = _mm_castpd_ps(_mm_load_sd((double*)(&values[17])));
c5_0 = _mm_add_ps(c5_0, _mm_mul_ps(a5_0, b5));
_mm_store_sd((double*)(&C[(i*16)+4]), _mm_castps_pd(c5_0));
__m128 c5_2 = _mm_load_ss(&C[(i*16)+6]);
__m128 a5_2 = _mm_load_ss(&values[19]);
c5_2 = _mm_add_ss(c5_2, _mm_mul_ss(a5_2, b5));
_mm_store_ss(&C[(i*16)+6], c5_2);
#else
C[(i*16)+4] += values[17] * B[(i*16)+5];
C[(i*16)+5] += values[18] * B[(i*16)+5];
C[(i*16)+6] += values[19] * B[(i*16)+5];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b6 = _mm_broadcast_ss(&B[(i*16)+6]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b6 = _mm_load_ss(&B[(i*16)+6]);
b6 = _mm_shuffle_ps(b6, b6, 0x00);
#endif
__m128 c6_0 = _mm_castpd_ps(_mm_load_sd((double*)(&C[(i*16)+4])));
__m128 a6_0 = _mm_castpd_ps(_mm_load_sd((double*)(&values[20])));
c6_0 = _mm_add_ps(c6_0, _mm_mul_ps(a6_0, b6));
_mm_store_sd((double*)(&C[(i*16)+4]), _mm_castps_pd(c6_0));
__m128 c6_2 = _mm_load_ss(&C[(i*16)+6]);
__m128 a6_2 = _mm_load_ss(&values[22]);
c6_2 = _mm_add_ss(c6_2, _mm_mul_ss(a6_2, b6));
_mm_store_ss(&C[(i*16)+6], c6_2);
#else
C[(i*16)+4] += values[20] * B[(i*16)+6];
C[(i*16)+5] += values[21] * B[(i*16)+6];
C[(i*16)+6] += values[22] * B[(i*16)+6];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b7 = _mm_broadcast_ss(&B[(i*16)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b7 = _mm_load_ss(&B[(i*16)+7]);
b7 = _mm_shuffle_ps(b7, b7, 0x00);
#endif
__m128 c7_0 = _mm_castpd_ps(_mm_load_sd((double*)(&C[(i*16)+1])));
__m128 a7_0 = _mm_castpd_ps(_mm_load_sd((double*)(&values[23])));
c7_0 = _mm_add_ps(c7_0, _mm_mul_ps(a7_0, b7));
_mm_store_sd((double*)(&C[(i*16)+1]), _mm_castps_pd(c7_0));
__m128 c7_2 = _mm_castpd_ps(_mm_load_sd((double*)(&C[(i*16)+7])));
__m128 a7_2 = _mm_castpd_ps(_mm_load_sd((double*)(&values[25])));
c7_2 = _mm_add_ps(c7_2, _mm_mul_ps(a7_2, b7));
_mm_store_sd((double*)(&C[(i*16)+7]), _mm_castps_pd(c7_2));
#else
C[(i*16)+1] += values[23] * B[(i*16)+7];
C[(i*16)+2] += values[24] * B[(i*16)+7];
C[(i*16)+7] += values[25] * B[(i*16)+7];
C[(i*16)+8] += values[26] * B[(i*16)+7];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b8 = _mm_broadcast_ss(&B[(i*16)+8]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b8 = _mm_load_ss(&B[(i*16)+8]);
b8 = _mm_shuffle_ps(b8, b8, 0x00);
#endif
__m128 c8_0 = _mm_castpd_ps(_mm_load_sd((double*)(&C[(i*16)+1])));
__m128 a8_0 = _mm_castpd_ps(_mm_load_sd((double*)(&values[27])));
c8_0 = _mm_add_ps(c8_0, _mm_mul_ps(a8_0, b8));
_mm_store_sd((double*)(&C[(i*16)+1]), _mm_castps_pd(c8_0));
__m128 c8_2 = _mm_castpd_ps(_mm_load_sd((double*)(&C[(i*16)+7])));
__m128 a8_2 = _mm_castpd_ps(_mm_load_sd((double*)(&values[29])));
c8_2 = _mm_add_ps(c8_2, _mm_mul_ps(a8_2, b8));
_mm_store_sd((double*)(&C[(i*16)+7]), _mm_castps_pd(c8_2));
#else
C[(i*16)+1] += values[27] * B[(i*16)+8];
C[(i*16)+2] += values[28] * B[(i*16)+8];
C[(i*16)+7] += values[29] * B[(i*16)+8];
C[(i*16)+8] += values[30] * B[(i*16)+8];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b9 = _mm_broadcast_ss(&B[(i*16)+9]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b9 = _mm_load_ss(&B[(i*16)+9]);
b9 = _mm_shuffle_ps(b9, b9, 0x00);
#endif
__m128 c9_0 = _mm_load_ss(&C[(i*16)+0]);
__m128 a9_0 = _mm_load_ss(&values[31]);
c9_0 = _mm_add_ss(c9_0, _mm_mul_ss(a9_0, b9));
_mm_store_ss(&C[(i*16)+0], c9_0);
__m128 c9_1 = _mm_load_ss(&C[(i*16)+3]);
__m128 a9_1 = _mm_load_ss(&values[32]);
c9_1 = _mm_add_ss(c9_1, _mm_mul_ss(a9_1, b9));
_mm_store_ss(&C[(i*16)+3], c9_1);
__m128 c9_2 = _mm_load_ss(&C[(i*16)+9]);
__m128 a9_2 = _mm_load_ss(&values[33]);
c9_2 = _mm_add_ss(c9_2, _mm_mul_ss(a9_2, b9));
_mm_store_ss(&C[(i*16)+9], c9_2);
#else
C[(i*16)+0] += values[31] * B[(i*16)+9];
C[(i*16)+3] += values[32] * B[(i*16)+9];
C[(i*16)+9] += values[33] * B[(i*16)+9];
#endif

}

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 612;
#endif

}

void ssparse_starMatrix_m20_n9_k9_ldA24_ldBna4_ldC24_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
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

void ssparse_starMatrix_m35_n9_k9_ldA40_ldBna5_ldC40_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
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

void ssparse_fM1DivM_m35_n9_k35_ldAna5_ldB40_ldC40_beta0_pfsigonly(const float* values, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma nounroll_and_jam
for (unsigned int i = 0; i < 9; i++)
{
  #pragma simd
  for (unsigned int m = 0; m < 35; m++) {
    C[(i*40)+m] = 0.0;
  }
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b0 = _mm_broadcast_ss(&B[(i*40)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b0 = _mm_load_ss(&B[(i*40)+0]);
b0 = _mm_shuffle_ps(b0, b0, 0x00);
#endif
__m128 c0_0 = _mm_load_ss(&C[(i*40)+0]);
__m128 a0_0 = _mm_load_ss(&values[0]);
c0_0 = _mm_add_ss(c0_0, _mm_mul_ss(a0_0, b0));
_mm_store_ss(&C[(i*40)+0], c0_0);
__m128 c0_1 = _mm_load_ss(&C[(i*40)+3]);
__m128 a0_1 = _mm_load_ss(&values[1]);
c0_1 = _mm_add_ss(c0_1, _mm_mul_ss(a0_1, b0));
_mm_store_ss(&C[(i*40)+3], c0_1);
__m128 c0_2 = _mm_load_ss(&C[(i*40)+9]);
__m128 a0_2 = _mm_load_ss(&values[2]);
c0_2 = _mm_add_ss(c0_2, _mm_mul_ss(a0_2, b0));
_mm_store_ss(&C[(i*40)+9], c0_2);
__m128 c0_3 = _mm_load_ss(&C[(i*40)+19]);
__m128 a0_3 = _mm_load_ss(&values[3]);
c0_3 = _mm_add_ss(c0_3, _mm_mul_ss(a0_3, b0));
_mm_store_ss(&C[(i*40)+19], c0_3);
__m128 c0_4 = _mm_load_ss(&C[(i*40)+34]);
__m128 a0_4 = _mm_load_ss(&values[4]);
c0_4 = _mm_add_ss(c0_4, _mm_mul_ss(a0_4, b0));
_mm_store_ss(&C[(i*40)+34], c0_4);
#else
C[(i*40)+0] += values[0] * B[(i*40)+0];
C[(i*40)+3] += values[1] * B[(i*40)+0];
C[(i*40)+9] += values[2] * B[(i*40)+0];
C[(i*40)+19] += values[3] * B[(i*40)+0];
C[(i*40)+34] += values[4] * B[(i*40)+0];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b1 = _mm_broadcast_ss(&B[(i*40)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b1 = _mm_load_ss(&B[(i*40)+1]);
b1 = _mm_shuffle_ps(b1, b1, 0x00);
#endif
__m128 c1_0 = _mm_load_ss(&C[(i*40)+1]);
__m128 a1_0 = _mm_load_ss(&values[5]);
c1_0 = _mm_add_ss(c1_0, _mm_mul_ss(a1_0, b1));
_mm_store_ss(&C[(i*40)+1], c1_0);
__m128 c1_1 = _mm_load_ss(&C[(i*40)+7]);
__m128 a1_1 = _mm_load_ss(&values[6]);
c1_1 = _mm_add_ss(c1_1, _mm_mul_ss(a1_1, b1));
_mm_store_ss(&C[(i*40)+7], c1_1);
__m128 c1_2 = _mm_load_ss(&C[(i*40)+17]);
__m128 a1_2 = _mm_load_ss(&values[7]);
c1_2 = _mm_add_ss(c1_2, _mm_mul_ss(a1_2, b1));
_mm_store_ss(&C[(i*40)+17], c1_2);
__m128 c1_3 = _mm_load_ss(&C[(i*40)+32]);
__m128 a1_3 = _mm_load_ss(&values[8]);
c1_3 = _mm_add_ss(c1_3, _mm_mul_ss(a1_3, b1));
_mm_store_ss(&C[(i*40)+32], c1_3);
#else
C[(i*40)+1] += values[5] * B[(i*40)+1];
C[(i*40)+7] += values[6] * B[(i*40)+1];
C[(i*40)+17] += values[7] * B[(i*40)+1];
C[(i*40)+32] += values[8] * B[(i*40)+1];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b2 = _mm_broadcast_ss(&B[(i*40)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b2 = _mm_load_ss(&B[(i*40)+2]);
b2 = _mm_shuffle_ps(b2, b2, 0x00);
#endif
__m128 c2_0 = _mm_load_ss(&C[(i*40)+2]);
__m128 a2_0 = _mm_load_ss(&values[9]);
c2_0 = _mm_add_ss(c2_0, _mm_mul_ss(a2_0, b2));
_mm_store_ss(&C[(i*40)+2], c2_0);
__m128 c2_1 = _mm_load_ss(&C[(i*40)+8]);
__m128 a2_1 = _mm_load_ss(&values[10]);
c2_1 = _mm_add_ss(c2_1, _mm_mul_ss(a2_1, b2));
_mm_store_ss(&C[(i*40)+8], c2_1);
__m128 c2_2 = _mm_load_ss(&C[(i*40)+18]);
__m128 a2_2 = _mm_load_ss(&values[11]);
c2_2 = _mm_add_ss(c2_2, _mm_mul_ss(a2_2, b2));
_mm_store_ss(&C[(i*40)+18], c2_2);
__m128 c2_3 = _mm_load_ss(&C[(i*40)+33]);
__m128 a2_3 = _mm_load_ss(&values[12]);
c2_3 = _mm_add_ss(c2_3, _mm_mul_ss(a2_3, b2));
_mm_store_ss(&C[(i*40)+33], c2_3);
#else
C[(i*40)+2] += values[9] * B[(i*40)+2];
C[(i*40)+8] += values[10] * B[(i*40)+2];
C[(i*40)+18] += values[11] * B[(i*40)+2];
C[(i*40)+33] += values[12] * B[(i*40)+2];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b3 = _mm_broadcast_ss(&B[(i*40)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b3 = _mm_load_ss(&B[(i*40)+3]);
b3 = _mm_shuffle_ps(b3, b3, 0x00);
#endif
__m128 c3_0 = _mm_load_ss(&C[(i*40)+0]);
__m128 a3_0 = _mm_load_ss(&values[13]);
c3_0 = _mm_add_ss(c3_0, _mm_mul_ss(a3_0, b3));
_mm_store_ss(&C[(i*40)+0], c3_0);
__m128 c3_1 = _mm_load_ss(&C[(i*40)+3]);
__m128 a3_1 = _mm_load_ss(&values[14]);
c3_1 = _mm_add_ss(c3_1, _mm_mul_ss(a3_1, b3));
_mm_store_ss(&C[(i*40)+3], c3_1);
__m128 c3_2 = _mm_load_ss(&C[(i*40)+9]);
__m128 a3_2 = _mm_load_ss(&values[15]);
c3_2 = _mm_add_ss(c3_2, _mm_mul_ss(a3_2, b3));
_mm_store_ss(&C[(i*40)+9], c3_2);
__m128 c3_3 = _mm_load_ss(&C[(i*40)+19]);
__m128 a3_3 = _mm_load_ss(&values[16]);
c3_3 = _mm_add_ss(c3_3, _mm_mul_ss(a3_3, b3));
_mm_store_ss(&C[(i*40)+19], c3_3);
__m128 c3_4 = _mm_load_ss(&C[(i*40)+34]);
__m128 a3_4 = _mm_load_ss(&values[17]);
c3_4 = _mm_add_ss(c3_4, _mm_mul_ss(a3_4, b3));
_mm_store_ss(&C[(i*40)+34], c3_4);
#else
C[(i*40)+0] += values[13] * B[(i*40)+3];
C[(i*40)+3] += values[14] * B[(i*40)+3];
C[(i*40)+9] += values[15] * B[(i*40)+3];
C[(i*40)+19] += values[16] * B[(i*40)+3];
C[(i*40)+34] += values[17] * B[(i*40)+3];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b4 = _mm_broadcast_ss(&B[(i*40)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b4 = _mm_load_ss(&B[(i*40)+4]);
b4 = _mm_shuffle_ps(b4, b4, 0x00);
#endif
__m128 c4_0 = _mm_load_ss(&C[(i*40)+4]);
__m128 a4_0 = _mm_load_ss(&values[18]);
c4_0 = _mm_add_ss(c4_0, _mm_mul_ss(a4_0, b4));
_mm_store_ss(&C[(i*40)+4], c4_0);
__m128 c4_1 = _mm_load_ss(&C[(i*40)+14]);
__m128 a4_1 = _mm_load_ss(&values[19]);
c4_1 = _mm_add_ss(c4_1, _mm_mul_ss(a4_1, b4));
_mm_store_ss(&C[(i*40)+14], c4_1);
__m128 c4_2 = _mm_load_ss(&C[(i*40)+29]);
__m128 a4_2 = _mm_load_ss(&values[20]);
c4_2 = _mm_add_ss(c4_2, _mm_mul_ss(a4_2, b4));
_mm_store_ss(&C[(i*40)+29], c4_2);
#else
C[(i*40)+4] += values[18] * B[(i*40)+4];
C[(i*40)+14] += values[19] * B[(i*40)+4];
C[(i*40)+29] += values[20] * B[(i*40)+4];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b5 = _mm_broadcast_ss(&B[(i*40)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b5 = _mm_load_ss(&B[(i*40)+5]);
b5 = _mm_shuffle_ps(b5, b5, 0x00);
#endif
__m128 c5_0 = _mm_load_ss(&C[(i*40)+5]);
__m128 a5_0 = _mm_load_ss(&values[21]);
c5_0 = _mm_add_ss(c5_0, _mm_mul_ss(a5_0, b5));
_mm_store_ss(&C[(i*40)+5], c5_0);
__m128 c5_1 = _mm_load_ss(&C[(i*40)+15]);
__m128 a5_1 = _mm_load_ss(&values[22]);
c5_1 = _mm_add_ss(c5_1, _mm_mul_ss(a5_1, b5));
_mm_store_ss(&C[(i*40)+15], c5_1);
__m128 c5_2 = _mm_load_ss(&C[(i*40)+30]);
__m128 a5_2 = _mm_load_ss(&values[23]);
c5_2 = _mm_add_ss(c5_2, _mm_mul_ss(a5_2, b5));
_mm_store_ss(&C[(i*40)+30], c5_2);
#else
C[(i*40)+5] += values[21] * B[(i*40)+5];
C[(i*40)+15] += values[22] * B[(i*40)+5];
C[(i*40)+30] += values[23] * B[(i*40)+5];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b6 = _mm_broadcast_ss(&B[(i*40)+6]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b6 = _mm_load_ss(&B[(i*40)+6]);
b6 = _mm_shuffle_ps(b6, b6, 0x00);
#endif
__m128 c6_0 = _mm_load_ss(&C[(i*40)+6]);
__m128 a6_0 = _mm_load_ss(&values[24]);
c6_0 = _mm_add_ss(c6_0, _mm_mul_ss(a6_0, b6));
_mm_store_ss(&C[(i*40)+6], c6_0);
__m128 c6_1 = _mm_load_ss(&C[(i*40)+16]);
__m128 a6_1 = _mm_load_ss(&values[25]);
c6_1 = _mm_add_ss(c6_1, _mm_mul_ss(a6_1, b6));
_mm_store_ss(&C[(i*40)+16], c6_1);
__m128 c6_2 = _mm_load_ss(&C[(i*40)+31]);
__m128 a6_2 = _mm_load_ss(&values[26]);
c6_2 = _mm_add_ss(c6_2, _mm_mul_ss(a6_2, b6));
_mm_store_ss(&C[(i*40)+31], c6_2);
#else
C[(i*40)+6] += values[24] * B[(i*40)+6];
C[(i*40)+16] += values[25] * B[(i*40)+6];
C[(i*40)+31] += values[26] * B[(i*40)+6];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b7 = _mm_broadcast_ss(&B[(i*40)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b7 = _mm_load_ss(&B[(i*40)+7]);
b7 = _mm_shuffle_ps(b7, b7, 0x00);
#endif
__m128 c7_0 = _mm_load_ss(&C[(i*40)+1]);
__m128 a7_0 = _mm_load_ss(&values[27]);
c7_0 = _mm_add_ss(c7_0, _mm_mul_ss(a7_0, b7));
_mm_store_ss(&C[(i*40)+1], c7_0);
__m128 c7_1 = _mm_load_ss(&C[(i*40)+7]);
__m128 a7_1 = _mm_load_ss(&values[28]);
c7_1 = _mm_add_ss(c7_1, _mm_mul_ss(a7_1, b7));
_mm_store_ss(&C[(i*40)+7], c7_1);
__m128 c7_2 = _mm_load_ss(&C[(i*40)+17]);
__m128 a7_2 = _mm_load_ss(&values[29]);
c7_2 = _mm_add_ss(c7_2, _mm_mul_ss(a7_2, b7));
_mm_store_ss(&C[(i*40)+17], c7_2);
__m128 c7_3 = _mm_load_ss(&C[(i*40)+32]);
__m128 a7_3 = _mm_load_ss(&values[30]);
c7_3 = _mm_add_ss(c7_3, _mm_mul_ss(a7_3, b7));
_mm_store_ss(&C[(i*40)+32], c7_3);
#else
C[(i*40)+1] += values[27] * B[(i*40)+7];
C[(i*40)+7] += values[28] * B[(i*40)+7];
C[(i*40)+17] += values[29] * B[(i*40)+7];
C[(i*40)+32] += values[30] * B[(i*40)+7];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b8 = _mm_broadcast_ss(&B[(i*40)+8]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b8 = _mm_load_ss(&B[(i*40)+8]);
b8 = _mm_shuffle_ps(b8, b8, 0x00);
#endif
__m128 c8_0 = _mm_load_ss(&C[(i*40)+2]);
__m128 a8_0 = _mm_load_ss(&values[31]);
c8_0 = _mm_add_ss(c8_0, _mm_mul_ss(a8_0, b8));
_mm_store_ss(&C[(i*40)+2], c8_0);
__m128 c8_1 = _mm_load_ss(&C[(i*40)+8]);
__m128 a8_1 = _mm_load_ss(&values[32]);
c8_1 = _mm_add_ss(c8_1, _mm_mul_ss(a8_1, b8));
_mm_store_ss(&C[(i*40)+8], c8_1);
__m128 c8_2 = _mm_load_ss(&C[(i*40)+18]);
__m128 a8_2 = _mm_load_ss(&values[33]);
c8_2 = _mm_add_ss(c8_2, _mm_mul_ss(a8_2, b8));
_mm_store_ss(&C[(i*40)+18], c8_2);
__m128 c8_3 = _mm_load_ss(&C[(i*40)+33]);
__m128 a8_3 = _mm_load_ss(&values[34]);
c8_3 = _mm_add_ss(c8_3, _mm_mul_ss(a8_3, b8));
_mm_store_ss(&C[(i*40)+33], c8_3);
#else
C[(i*40)+2] += values[31] * B[(i*40)+8];
C[(i*40)+8] += values[32] * B[(i*40)+8];
C[(i*40)+18] += values[33] * B[(i*40)+8];
C[(i*40)+33] += values[34] * B[(i*40)+8];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b9 = _mm_broadcast_ss(&B[(i*40)+9]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b9 = _mm_load_ss(&B[(i*40)+9]);
b9 = _mm_shuffle_ps(b9, b9, 0x00);
#endif
__m128 c9_0 = _mm_load_ss(&C[(i*40)+0]);
__m128 a9_0 = _mm_load_ss(&values[35]);
c9_0 = _mm_add_ss(c9_0, _mm_mul_ss(a9_0, b9));
_mm_store_ss(&C[(i*40)+0], c9_0);
__m128 c9_1 = _mm_load_ss(&C[(i*40)+3]);
__m128 a9_1 = _mm_load_ss(&values[36]);
c9_1 = _mm_add_ss(c9_1, _mm_mul_ss(a9_1, b9));
_mm_store_ss(&C[(i*40)+3], c9_1);
__m128 c9_2 = _mm_load_ss(&C[(i*40)+9]);
__m128 a9_2 = _mm_load_ss(&values[37]);
c9_2 = _mm_add_ss(c9_2, _mm_mul_ss(a9_2, b9));
_mm_store_ss(&C[(i*40)+9], c9_2);
__m128 c9_3 = _mm_load_ss(&C[(i*40)+19]);
__m128 a9_3 = _mm_load_ss(&values[38]);
c9_3 = _mm_add_ss(c9_3, _mm_mul_ss(a9_3, b9));
_mm_store_ss(&C[(i*40)+19], c9_3);
__m128 c9_4 = _mm_load_ss(&C[(i*40)+34]);
__m128 a9_4 = _mm_load_ss(&values[39]);
c9_4 = _mm_add_ss(c9_4, _mm_mul_ss(a9_4, b9));
_mm_store_ss(&C[(i*40)+34], c9_4);
#else
C[(i*40)+0] += values[35] * B[(i*40)+9];
C[(i*40)+3] += values[36] * B[(i*40)+9];
C[(i*40)+9] += values[37] * B[(i*40)+9];
C[(i*40)+19] += values[38] * B[(i*40)+9];
C[(i*40)+34] += values[39] * B[(i*40)+9];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b10 = _mm_broadcast_ss(&B[(i*40)+10]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b10 = _mm_load_ss(&B[(i*40)+10]);
b10 = _mm_shuffle_ps(b10, b10, 0x00);
#endif
__m128 c10_0 = _mm_load_ss(&C[(i*40)+10]);
__m128 a10_0 = _mm_load_ss(&values[40]);
c10_0 = _mm_add_ss(c10_0, _mm_mul_ss(a10_0, b10));
_mm_store_ss(&C[(i*40)+10], c10_0);
__m128 c10_1 = _mm_load_ss(&C[(i*40)+25]);
__m128 a10_1 = _mm_load_ss(&values[41]);
c10_1 = _mm_add_ss(c10_1, _mm_mul_ss(a10_1, b10));
_mm_store_ss(&C[(i*40)+25], c10_1);
#else
C[(i*40)+10] += values[40] * B[(i*40)+10];
C[(i*40)+25] += values[41] * B[(i*40)+10];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b11 = _mm_broadcast_ss(&B[(i*40)+11]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b11 = _mm_load_ss(&B[(i*40)+11]);
b11 = _mm_shuffle_ps(b11, b11, 0x00);
#endif
__m128 c11_0 = _mm_load_ss(&C[(i*40)+11]);
__m128 a11_0 = _mm_load_ss(&values[42]);
c11_0 = _mm_add_ss(c11_0, _mm_mul_ss(a11_0, b11));
_mm_store_ss(&C[(i*40)+11], c11_0);
__m128 c11_1 = _mm_load_ss(&C[(i*40)+26]);
__m128 a11_1 = _mm_load_ss(&values[43]);
c11_1 = _mm_add_ss(c11_1, _mm_mul_ss(a11_1, b11));
_mm_store_ss(&C[(i*40)+26], c11_1);
#else
C[(i*40)+11] += values[42] * B[(i*40)+11];
C[(i*40)+26] += values[43] * B[(i*40)+11];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b12 = _mm_broadcast_ss(&B[(i*40)+12]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b12 = _mm_load_ss(&B[(i*40)+12]);
b12 = _mm_shuffle_ps(b12, b12, 0x00);
#endif
__m128 c12_0 = _mm_load_ss(&C[(i*40)+12]);
__m128 a12_0 = _mm_load_ss(&values[44]);
c12_0 = _mm_add_ss(c12_0, _mm_mul_ss(a12_0, b12));
_mm_store_ss(&C[(i*40)+12], c12_0);
__m128 c12_1 = _mm_load_ss(&C[(i*40)+27]);
__m128 a12_1 = _mm_load_ss(&values[45]);
c12_1 = _mm_add_ss(c12_1, _mm_mul_ss(a12_1, b12));
_mm_store_ss(&C[(i*40)+27], c12_1);
#else
C[(i*40)+12] += values[44] * B[(i*40)+12];
C[(i*40)+27] += values[45] * B[(i*40)+12];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b13 = _mm_broadcast_ss(&B[(i*40)+13]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b13 = _mm_load_ss(&B[(i*40)+13]);
b13 = _mm_shuffle_ps(b13, b13, 0x00);
#endif
__m128 c13_0 = _mm_load_ss(&C[(i*40)+13]);
__m128 a13_0 = _mm_load_ss(&values[46]);
c13_0 = _mm_add_ss(c13_0, _mm_mul_ss(a13_0, b13));
_mm_store_ss(&C[(i*40)+13], c13_0);
__m128 c13_1 = _mm_load_ss(&C[(i*40)+28]);
__m128 a13_1 = _mm_load_ss(&values[47]);
c13_1 = _mm_add_ss(c13_1, _mm_mul_ss(a13_1, b13));
_mm_store_ss(&C[(i*40)+28], c13_1);
#else
C[(i*40)+13] += values[46] * B[(i*40)+13];
C[(i*40)+28] += values[47] * B[(i*40)+13];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b14 = _mm_broadcast_ss(&B[(i*40)+14]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b14 = _mm_load_ss(&B[(i*40)+14]);
b14 = _mm_shuffle_ps(b14, b14, 0x00);
#endif
__m128 c14_0 = _mm_load_ss(&C[(i*40)+4]);
__m128 a14_0 = _mm_load_ss(&values[48]);
c14_0 = _mm_add_ss(c14_0, _mm_mul_ss(a14_0, b14));
_mm_store_ss(&C[(i*40)+4], c14_0);
__m128 c14_1 = _mm_load_ss(&C[(i*40)+14]);
__m128 a14_1 = _mm_load_ss(&values[49]);
c14_1 = _mm_add_ss(c14_1, _mm_mul_ss(a14_1, b14));
_mm_store_ss(&C[(i*40)+14], c14_1);
__m128 c14_2 = _mm_load_ss(&C[(i*40)+29]);
__m128 a14_2 = _mm_load_ss(&values[50]);
c14_2 = _mm_add_ss(c14_2, _mm_mul_ss(a14_2, b14));
_mm_store_ss(&C[(i*40)+29], c14_2);
#else
C[(i*40)+4] += values[48] * B[(i*40)+14];
C[(i*40)+14] += values[49] * B[(i*40)+14];
C[(i*40)+29] += values[50] * B[(i*40)+14];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b15 = _mm_broadcast_ss(&B[(i*40)+15]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b15 = _mm_load_ss(&B[(i*40)+15]);
b15 = _mm_shuffle_ps(b15, b15, 0x00);
#endif
__m128 c15_0 = _mm_load_ss(&C[(i*40)+5]);
__m128 a15_0 = _mm_load_ss(&values[51]);
c15_0 = _mm_add_ss(c15_0, _mm_mul_ss(a15_0, b15));
_mm_store_ss(&C[(i*40)+5], c15_0);
__m128 c15_1 = _mm_load_ss(&C[(i*40)+15]);
__m128 a15_1 = _mm_load_ss(&values[52]);
c15_1 = _mm_add_ss(c15_1, _mm_mul_ss(a15_1, b15));
_mm_store_ss(&C[(i*40)+15], c15_1);
__m128 c15_2 = _mm_load_ss(&C[(i*40)+30]);
__m128 a15_2 = _mm_load_ss(&values[53]);
c15_2 = _mm_add_ss(c15_2, _mm_mul_ss(a15_2, b15));
_mm_store_ss(&C[(i*40)+30], c15_2);
#else
C[(i*40)+5] += values[51] * B[(i*40)+15];
C[(i*40)+15] += values[52] * B[(i*40)+15];
C[(i*40)+30] += values[53] * B[(i*40)+15];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b16 = _mm_broadcast_ss(&B[(i*40)+16]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b16 = _mm_load_ss(&B[(i*40)+16]);
b16 = _mm_shuffle_ps(b16, b16, 0x00);
#endif
__m128 c16_0 = _mm_load_ss(&C[(i*40)+6]);
__m128 a16_0 = _mm_load_ss(&values[54]);
c16_0 = _mm_add_ss(c16_0, _mm_mul_ss(a16_0, b16));
_mm_store_ss(&C[(i*40)+6], c16_0);
__m128 c16_1 = _mm_load_ss(&C[(i*40)+16]);
__m128 a16_1 = _mm_load_ss(&values[55]);
c16_1 = _mm_add_ss(c16_1, _mm_mul_ss(a16_1, b16));
_mm_store_ss(&C[(i*40)+16], c16_1);
__m128 c16_2 = _mm_load_ss(&C[(i*40)+31]);
__m128 a16_2 = _mm_load_ss(&values[56]);
c16_2 = _mm_add_ss(c16_2, _mm_mul_ss(a16_2, b16));
_mm_store_ss(&C[(i*40)+31], c16_2);
#else
C[(i*40)+6] += values[54] * B[(i*40)+16];
C[(i*40)+16] += values[55] * B[(i*40)+16];
C[(i*40)+31] += values[56] * B[(i*40)+16];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b17 = _mm_broadcast_ss(&B[(i*40)+17]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b17 = _mm_load_ss(&B[(i*40)+17]);
b17 = _mm_shuffle_ps(b17, b17, 0x00);
#endif
__m128 c17_0 = _mm_load_ss(&C[(i*40)+1]);
__m128 a17_0 = _mm_load_ss(&values[57]);
c17_0 = _mm_add_ss(c17_0, _mm_mul_ss(a17_0, b17));
_mm_store_ss(&C[(i*40)+1], c17_0);
__m128 c17_1 = _mm_load_ss(&C[(i*40)+7]);
__m128 a17_1 = _mm_load_ss(&values[58]);
c17_1 = _mm_add_ss(c17_1, _mm_mul_ss(a17_1, b17));
_mm_store_ss(&C[(i*40)+7], c17_1);
__m128 c17_2 = _mm_load_ss(&C[(i*40)+17]);
__m128 a17_2 = _mm_load_ss(&values[59]);
c17_2 = _mm_add_ss(c17_2, _mm_mul_ss(a17_2, b17));
_mm_store_ss(&C[(i*40)+17], c17_2);
__m128 c17_3 = _mm_load_ss(&C[(i*40)+32]);
__m128 a17_3 = _mm_load_ss(&values[60]);
c17_3 = _mm_add_ss(c17_3, _mm_mul_ss(a17_3, b17));
_mm_store_ss(&C[(i*40)+32], c17_3);
#else
C[(i*40)+1] += values[57] * B[(i*40)+17];
C[(i*40)+7] += values[58] * B[(i*40)+17];
C[(i*40)+17] += values[59] * B[(i*40)+17];
C[(i*40)+32] += values[60] * B[(i*40)+17];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b18 = _mm_broadcast_ss(&B[(i*40)+18]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b18 = _mm_load_ss(&B[(i*40)+18]);
b18 = _mm_shuffle_ps(b18, b18, 0x00);
#endif
__m128 c18_0 = _mm_load_ss(&C[(i*40)+2]);
__m128 a18_0 = _mm_load_ss(&values[61]);
c18_0 = _mm_add_ss(c18_0, _mm_mul_ss(a18_0, b18));
_mm_store_ss(&C[(i*40)+2], c18_0);
__m128 c18_1 = _mm_load_ss(&C[(i*40)+8]);
__m128 a18_1 = _mm_load_ss(&values[62]);
c18_1 = _mm_add_ss(c18_1, _mm_mul_ss(a18_1, b18));
_mm_store_ss(&C[(i*40)+8], c18_1);
__m128 c18_2 = _mm_load_ss(&C[(i*40)+18]);
__m128 a18_2 = _mm_load_ss(&values[63]);
c18_2 = _mm_add_ss(c18_2, _mm_mul_ss(a18_2, b18));
_mm_store_ss(&C[(i*40)+18], c18_2);
__m128 c18_3 = _mm_load_ss(&C[(i*40)+33]);
__m128 a18_3 = _mm_load_ss(&values[64]);
c18_3 = _mm_add_ss(c18_3, _mm_mul_ss(a18_3, b18));
_mm_store_ss(&C[(i*40)+33], c18_3);
#else
C[(i*40)+2] += values[61] * B[(i*40)+18];
C[(i*40)+8] += values[62] * B[(i*40)+18];
C[(i*40)+18] += values[63] * B[(i*40)+18];
C[(i*40)+33] += values[64] * B[(i*40)+18];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b19 = _mm_broadcast_ss(&B[(i*40)+19]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b19 = _mm_load_ss(&B[(i*40)+19]);
b19 = _mm_shuffle_ps(b19, b19, 0x00);
#endif
__m128 c19_0 = _mm_load_ss(&C[(i*40)+0]);
__m128 a19_0 = _mm_load_ss(&values[65]);
c19_0 = _mm_add_ss(c19_0, _mm_mul_ss(a19_0, b19));
_mm_store_ss(&C[(i*40)+0], c19_0);
__m128 c19_1 = _mm_load_ss(&C[(i*40)+3]);
__m128 a19_1 = _mm_load_ss(&values[66]);
c19_1 = _mm_add_ss(c19_1, _mm_mul_ss(a19_1, b19));
_mm_store_ss(&C[(i*40)+3], c19_1);
__m128 c19_2 = _mm_load_ss(&C[(i*40)+9]);
__m128 a19_2 = _mm_load_ss(&values[67]);
c19_2 = _mm_add_ss(c19_2, _mm_mul_ss(a19_2, b19));
_mm_store_ss(&C[(i*40)+9], c19_2);
__m128 c19_3 = _mm_load_ss(&C[(i*40)+19]);
__m128 a19_3 = _mm_load_ss(&values[68]);
c19_3 = _mm_add_ss(c19_3, _mm_mul_ss(a19_3, b19));
_mm_store_ss(&C[(i*40)+19], c19_3);
__m128 c19_4 = _mm_load_ss(&C[(i*40)+34]);
__m128 a19_4 = _mm_load_ss(&values[69]);
c19_4 = _mm_add_ss(c19_4, _mm_mul_ss(a19_4, b19));
_mm_store_ss(&C[(i*40)+34], c19_4);
#else
C[(i*40)+0] += values[65] * B[(i*40)+19];
C[(i*40)+3] += values[66] * B[(i*40)+19];
C[(i*40)+9] += values[67] * B[(i*40)+19];
C[(i*40)+19] += values[68] * B[(i*40)+19];
C[(i*40)+34] += values[69] * B[(i*40)+19];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b20 = _mm_broadcast_ss(&B[(i*40)+20]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b20 = _mm_load_ss(&B[(i*40)+20]);
b20 = _mm_shuffle_ps(b20, b20, 0x00);
#endif
__m128 c20_0 = _mm_load_ss(&C[(i*40)+20]);
__m128 a20_0 = _mm_load_ss(&values[70]);
c20_0 = _mm_add_ss(c20_0, _mm_mul_ss(a20_0, b20));
_mm_store_ss(&C[(i*40)+20], c20_0);
#else
C[(i*40)+20] += values[70] * B[(i*40)+20];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b21 = _mm_broadcast_ss(&B[(i*40)+21]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b21 = _mm_load_ss(&B[(i*40)+21]);
b21 = _mm_shuffle_ps(b21, b21, 0x00);
#endif
__m128 c21_0 = _mm_load_ss(&C[(i*40)+21]);
__m128 a21_0 = _mm_load_ss(&values[71]);
c21_0 = _mm_add_ss(c21_0, _mm_mul_ss(a21_0, b21));
_mm_store_ss(&C[(i*40)+21], c21_0);
#else
C[(i*40)+21] += values[71] * B[(i*40)+21];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b22 = _mm_broadcast_ss(&B[(i*40)+22]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b22 = _mm_load_ss(&B[(i*40)+22]);
b22 = _mm_shuffle_ps(b22, b22, 0x00);
#endif
__m128 c22_0 = _mm_load_ss(&C[(i*40)+22]);
__m128 a22_0 = _mm_load_ss(&values[72]);
c22_0 = _mm_add_ss(c22_0, _mm_mul_ss(a22_0, b22));
_mm_store_ss(&C[(i*40)+22], c22_0);
#else
C[(i*40)+22] += values[72] * B[(i*40)+22];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b23 = _mm_broadcast_ss(&B[(i*40)+23]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b23 = _mm_load_ss(&B[(i*40)+23]);
b23 = _mm_shuffle_ps(b23, b23, 0x00);
#endif
__m128 c23_0 = _mm_load_ss(&C[(i*40)+23]);
__m128 a23_0 = _mm_load_ss(&values[73]);
c23_0 = _mm_add_ss(c23_0, _mm_mul_ss(a23_0, b23));
_mm_store_ss(&C[(i*40)+23], c23_0);
#else
C[(i*40)+23] += values[73] * B[(i*40)+23];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b24 = _mm_broadcast_ss(&B[(i*40)+24]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b24 = _mm_load_ss(&B[(i*40)+24]);
b24 = _mm_shuffle_ps(b24, b24, 0x00);
#endif
__m128 c24_0 = _mm_load_ss(&C[(i*40)+24]);
__m128 a24_0 = _mm_load_ss(&values[74]);
c24_0 = _mm_add_ss(c24_0, _mm_mul_ss(a24_0, b24));
_mm_store_ss(&C[(i*40)+24], c24_0);
#else
C[(i*40)+24] += values[74] * B[(i*40)+24];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b25 = _mm_broadcast_ss(&B[(i*40)+25]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b25 = _mm_load_ss(&B[(i*40)+25]);
b25 = _mm_shuffle_ps(b25, b25, 0x00);
#endif
__m128 c25_0 = _mm_load_ss(&C[(i*40)+10]);
__m128 a25_0 = _mm_load_ss(&values[75]);
c25_0 = _mm_add_ss(c25_0, _mm_mul_ss(a25_0, b25));
_mm_store_ss(&C[(i*40)+10], c25_0);
__m128 c25_1 = _mm_load_ss(&C[(i*40)+25]);
__m128 a25_1 = _mm_load_ss(&values[76]);
c25_1 = _mm_add_ss(c25_1, _mm_mul_ss(a25_1, b25));
_mm_store_ss(&C[(i*40)+25], c25_1);
#else
C[(i*40)+10] += values[75] * B[(i*40)+25];
C[(i*40)+25] += values[76] * B[(i*40)+25];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b26 = _mm_broadcast_ss(&B[(i*40)+26]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b26 = _mm_load_ss(&B[(i*40)+26]);
b26 = _mm_shuffle_ps(b26, b26, 0x00);
#endif
__m128 c26_0 = _mm_load_ss(&C[(i*40)+11]);
__m128 a26_0 = _mm_load_ss(&values[77]);
c26_0 = _mm_add_ss(c26_0, _mm_mul_ss(a26_0, b26));
_mm_store_ss(&C[(i*40)+11], c26_0);
__m128 c26_1 = _mm_load_ss(&C[(i*40)+26]);
__m128 a26_1 = _mm_load_ss(&values[78]);
c26_1 = _mm_add_ss(c26_1, _mm_mul_ss(a26_1, b26));
_mm_store_ss(&C[(i*40)+26], c26_1);
#else
C[(i*40)+11] += values[77] * B[(i*40)+26];
C[(i*40)+26] += values[78] * B[(i*40)+26];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b27 = _mm_broadcast_ss(&B[(i*40)+27]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b27 = _mm_load_ss(&B[(i*40)+27]);
b27 = _mm_shuffle_ps(b27, b27, 0x00);
#endif
__m128 c27_0 = _mm_load_ss(&C[(i*40)+12]);
__m128 a27_0 = _mm_load_ss(&values[79]);
c27_0 = _mm_add_ss(c27_0, _mm_mul_ss(a27_0, b27));
_mm_store_ss(&C[(i*40)+12], c27_0);
__m128 c27_1 = _mm_load_ss(&C[(i*40)+27]);
__m128 a27_1 = _mm_load_ss(&values[80]);
c27_1 = _mm_add_ss(c27_1, _mm_mul_ss(a27_1, b27));
_mm_store_ss(&C[(i*40)+27], c27_1);
#else
C[(i*40)+12] += values[79] * B[(i*40)+27];
C[(i*40)+27] += values[80] * B[(i*40)+27];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b28 = _mm_broadcast_ss(&B[(i*40)+28]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b28 = _mm_load_ss(&B[(i*40)+28]);
b28 = _mm_shuffle_ps(b28, b28, 0x00);
#endif
__m128 c28_0 = _mm_load_ss(&C[(i*40)+13]);
__m128 a28_0 = _mm_load_ss(&values[81]);
c28_0 = _mm_add_ss(c28_0, _mm_mul_ss(a28_0, b28));
_mm_store_ss(&C[(i*40)+13], c28_0);
__m128 c28_1 = _mm_load_ss(&C[(i*40)+28]);
__m128 a28_1 = _mm_load_ss(&values[82]);
c28_1 = _mm_add_ss(c28_1, _mm_mul_ss(a28_1, b28));
_mm_store_ss(&C[(i*40)+28], c28_1);
#else
C[(i*40)+13] += values[81] * B[(i*40)+28];
C[(i*40)+28] += values[82] * B[(i*40)+28];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b29 = _mm_broadcast_ss(&B[(i*40)+29]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b29 = _mm_load_ss(&B[(i*40)+29]);
b29 = _mm_shuffle_ps(b29, b29, 0x00);
#endif
__m128 c29_0 = _mm_load_ss(&C[(i*40)+4]);
__m128 a29_0 = _mm_load_ss(&values[83]);
c29_0 = _mm_add_ss(c29_0, _mm_mul_ss(a29_0, b29));
_mm_store_ss(&C[(i*40)+4], c29_0);
__m128 c29_1 = _mm_load_ss(&C[(i*40)+14]);
__m128 a29_1 = _mm_load_ss(&values[84]);
c29_1 = _mm_add_ss(c29_1, _mm_mul_ss(a29_1, b29));
_mm_store_ss(&C[(i*40)+14], c29_1);
__m128 c29_2 = _mm_load_ss(&C[(i*40)+29]);
__m128 a29_2 = _mm_load_ss(&values[85]);
c29_2 = _mm_add_ss(c29_2, _mm_mul_ss(a29_2, b29));
_mm_store_ss(&C[(i*40)+29], c29_2);
#else
C[(i*40)+4] += values[83] * B[(i*40)+29];
C[(i*40)+14] += values[84] * B[(i*40)+29];
C[(i*40)+29] += values[85] * B[(i*40)+29];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b30 = _mm_broadcast_ss(&B[(i*40)+30]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b30 = _mm_load_ss(&B[(i*40)+30]);
b30 = _mm_shuffle_ps(b30, b30, 0x00);
#endif
__m128 c30_0 = _mm_load_ss(&C[(i*40)+5]);
__m128 a30_0 = _mm_load_ss(&values[86]);
c30_0 = _mm_add_ss(c30_0, _mm_mul_ss(a30_0, b30));
_mm_store_ss(&C[(i*40)+5], c30_0);
__m128 c30_1 = _mm_load_ss(&C[(i*40)+15]);
__m128 a30_1 = _mm_load_ss(&values[87]);
c30_1 = _mm_add_ss(c30_1, _mm_mul_ss(a30_1, b30));
_mm_store_ss(&C[(i*40)+15], c30_1);
__m128 c30_2 = _mm_load_ss(&C[(i*40)+30]);
__m128 a30_2 = _mm_load_ss(&values[88]);
c30_2 = _mm_add_ss(c30_2, _mm_mul_ss(a30_2, b30));
_mm_store_ss(&C[(i*40)+30], c30_2);
#else
C[(i*40)+5] += values[86] * B[(i*40)+30];
C[(i*40)+15] += values[87] * B[(i*40)+30];
C[(i*40)+30] += values[88] * B[(i*40)+30];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b31 = _mm_broadcast_ss(&B[(i*40)+31]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b31 = _mm_load_ss(&B[(i*40)+31]);
b31 = _mm_shuffle_ps(b31, b31, 0x00);
#endif
__m128 c31_0 = _mm_load_ss(&C[(i*40)+6]);
__m128 a31_0 = _mm_load_ss(&values[89]);
c31_0 = _mm_add_ss(c31_0, _mm_mul_ss(a31_0, b31));
_mm_store_ss(&C[(i*40)+6], c31_0);
__m128 c31_1 = _mm_load_ss(&C[(i*40)+16]);
__m128 a31_1 = _mm_load_ss(&values[90]);
c31_1 = _mm_add_ss(c31_1, _mm_mul_ss(a31_1, b31));
_mm_store_ss(&C[(i*40)+16], c31_1);
__m128 c31_2 = _mm_load_ss(&C[(i*40)+31]);
__m128 a31_2 = _mm_load_ss(&values[91]);
c31_2 = _mm_add_ss(c31_2, _mm_mul_ss(a31_2, b31));
_mm_store_ss(&C[(i*40)+31], c31_2);
#else
C[(i*40)+6] += values[89] * B[(i*40)+31];
C[(i*40)+16] += values[90] * B[(i*40)+31];
C[(i*40)+31] += values[91] * B[(i*40)+31];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b32 = _mm_broadcast_ss(&B[(i*40)+32]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b32 = _mm_load_ss(&B[(i*40)+32]);
b32 = _mm_shuffle_ps(b32, b32, 0x00);
#endif
__m128 c32_0 = _mm_load_ss(&C[(i*40)+1]);
__m128 a32_0 = _mm_load_ss(&values[92]);
c32_0 = _mm_add_ss(c32_0, _mm_mul_ss(a32_0, b32));
_mm_store_ss(&C[(i*40)+1], c32_0);
__m128 c32_1 = _mm_load_ss(&C[(i*40)+7]);
__m128 a32_1 = _mm_load_ss(&values[93]);
c32_1 = _mm_add_ss(c32_1, _mm_mul_ss(a32_1, b32));
_mm_store_ss(&C[(i*40)+7], c32_1);
__m128 c32_2 = _mm_load_ss(&C[(i*40)+17]);
__m128 a32_2 = _mm_load_ss(&values[94]);
c32_2 = _mm_add_ss(c32_2, _mm_mul_ss(a32_2, b32));
_mm_store_ss(&C[(i*40)+17], c32_2);
__m128 c32_3 = _mm_load_ss(&C[(i*40)+32]);
__m128 a32_3 = _mm_load_ss(&values[95]);
c32_3 = _mm_add_ss(c32_3, _mm_mul_ss(a32_3, b32));
_mm_store_ss(&C[(i*40)+32], c32_3);
#else
C[(i*40)+1] += values[92] * B[(i*40)+32];
C[(i*40)+7] += values[93] * B[(i*40)+32];
C[(i*40)+17] += values[94] * B[(i*40)+32];
C[(i*40)+32] += values[95] * B[(i*40)+32];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b33 = _mm_broadcast_ss(&B[(i*40)+33]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b33 = _mm_load_ss(&B[(i*40)+33]);
b33 = _mm_shuffle_ps(b33, b33, 0x00);
#endif
__m128 c33_0 = _mm_load_ss(&C[(i*40)+2]);
__m128 a33_0 = _mm_load_ss(&values[96]);
c33_0 = _mm_add_ss(c33_0, _mm_mul_ss(a33_0, b33));
_mm_store_ss(&C[(i*40)+2], c33_0);
__m128 c33_1 = _mm_load_ss(&C[(i*40)+8]);
__m128 a33_1 = _mm_load_ss(&values[97]);
c33_1 = _mm_add_ss(c33_1, _mm_mul_ss(a33_1, b33));
_mm_store_ss(&C[(i*40)+8], c33_1);
__m128 c33_2 = _mm_load_ss(&C[(i*40)+18]);
__m128 a33_2 = _mm_load_ss(&values[98]);
c33_2 = _mm_add_ss(c33_2, _mm_mul_ss(a33_2, b33));
_mm_store_ss(&C[(i*40)+18], c33_2);
__m128 c33_3 = _mm_load_ss(&C[(i*40)+33]);
__m128 a33_3 = _mm_load_ss(&values[99]);
c33_3 = _mm_add_ss(c33_3, _mm_mul_ss(a33_3, b33));
_mm_store_ss(&C[(i*40)+33], c33_3);
#else
C[(i*40)+2] += values[96] * B[(i*40)+33];
C[(i*40)+8] += values[97] * B[(i*40)+33];
C[(i*40)+18] += values[98] * B[(i*40)+33];
C[(i*40)+33] += values[99] * B[(i*40)+33];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b34 = _mm_broadcast_ss(&B[(i*40)+34]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b34 = _mm_load_ss(&B[(i*40)+34]);
b34 = _mm_shuffle_ps(b34, b34, 0x00);
#endif
__m128 c34_0 = _mm_load_ss(&C[(i*40)+0]);
__m128 a34_0 = _mm_load_ss(&values[100]);
c34_0 = _mm_add_ss(c34_0, _mm_mul_ss(a34_0, b34));
_mm_store_ss(&C[(i*40)+0], c34_0);
__m128 c34_1 = _mm_load_ss(&C[(i*40)+3]);
__m128 a34_1 = _mm_load_ss(&values[101]);
c34_1 = _mm_add_ss(c34_1, _mm_mul_ss(a34_1, b34));
_mm_store_ss(&C[(i*40)+3], c34_1);
__m128 c34_2 = _mm_load_ss(&C[(i*40)+9]);
__m128 a34_2 = _mm_load_ss(&values[102]);
c34_2 = _mm_add_ss(c34_2, _mm_mul_ss(a34_2, b34));
_mm_store_ss(&C[(i*40)+9], c34_2);
__m128 c34_3 = _mm_load_ss(&C[(i*40)+19]);
__m128 a34_3 = _mm_load_ss(&values[103]);
c34_3 = _mm_add_ss(c34_3, _mm_mul_ss(a34_3, b34));
_mm_store_ss(&C[(i*40)+19], c34_3);
__m128 c34_4 = _mm_load_ss(&C[(i*40)+34]);
__m128 a34_4 = _mm_load_ss(&values[104]);
c34_4 = _mm_add_ss(c34_4, _mm_mul_ss(a34_4, b34));
_mm_store_ss(&C[(i*40)+34], c34_4);
#else
C[(i*40)+0] += values[100] * B[(i*40)+34];
C[(i*40)+3] += values[101] * B[(i*40)+34];
C[(i*40)+9] += values[102] * B[(i*40)+34];
C[(i*40)+19] += values[103] * B[(i*40)+34];
C[(i*40)+34] += values[104] * B[(i*40)+34];
#endif

}

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 1890;
#endif

}

void ssparse_starMatrix_m56_n9_k9_ldA56_ldBna6_ldC56_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
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

void ssparse_fM1DivM_m56_n9_k56_ldAna6_ldB56_ldC56_beta0_pfsigonly(const float* values, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
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
__m128 b0 = _mm_broadcast_ss(&B[(i*56)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b0 = _mm_load_ss(&B[(i*56)+0]);
b0 = _mm_shuffle_ps(b0, b0, 0x00);
#endif
__m128 c0_0 = _mm_load_ss(&C[(i*56)+0]);
__m128 a0_0 = _mm_load_ss(&values[0]);
c0_0 = _mm_add_ss(c0_0, _mm_mul_ss(a0_0, b0));
_mm_store_ss(&C[(i*56)+0], c0_0);
__m128 c0_1 = _mm_load_ss(&C[(i*56)+3]);
__m128 a0_1 = _mm_load_ss(&values[1]);
c0_1 = _mm_add_ss(c0_1, _mm_mul_ss(a0_1, b0));
_mm_store_ss(&C[(i*56)+3], c0_1);
__m128 c0_2 = _mm_load_ss(&C[(i*56)+9]);
__m128 a0_2 = _mm_load_ss(&values[2]);
c0_2 = _mm_add_ss(c0_2, _mm_mul_ss(a0_2, b0));
_mm_store_ss(&C[(i*56)+9], c0_2);
__m128 c0_3 = _mm_load_ss(&C[(i*56)+19]);
__m128 a0_3 = _mm_load_ss(&values[3]);
c0_3 = _mm_add_ss(c0_3, _mm_mul_ss(a0_3, b0));
_mm_store_ss(&C[(i*56)+19], c0_3);
__m128 c0_4 = _mm_load_ss(&C[(i*56)+34]);
__m128 a0_4 = _mm_load_ss(&values[4]);
c0_4 = _mm_add_ss(c0_4, _mm_mul_ss(a0_4, b0));
_mm_store_ss(&C[(i*56)+34], c0_4);
__m128 c0_5 = _mm_load_ss(&C[(i*56)+55]);
__m128 a0_5 = _mm_load_ss(&values[5]);
c0_5 = _mm_add_ss(c0_5, _mm_mul_ss(a0_5, b0));
_mm_store_ss(&C[(i*56)+55], c0_5);
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
__m128 b1 = _mm_broadcast_ss(&B[(i*56)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b1 = _mm_load_ss(&B[(i*56)+1]);
b1 = _mm_shuffle_ps(b1, b1, 0x00);
#endif
__m128 c1_0 = _mm_load_ss(&C[(i*56)+1]);
__m128 a1_0 = _mm_load_ss(&values[6]);
c1_0 = _mm_add_ss(c1_0, _mm_mul_ss(a1_0, b1));
_mm_store_ss(&C[(i*56)+1], c1_0);
__m128 c1_1 = _mm_load_ss(&C[(i*56)+7]);
__m128 a1_1 = _mm_load_ss(&values[7]);
c1_1 = _mm_add_ss(c1_1, _mm_mul_ss(a1_1, b1));
_mm_store_ss(&C[(i*56)+7], c1_1);
__m128 c1_2 = _mm_load_ss(&C[(i*56)+17]);
__m128 a1_2 = _mm_load_ss(&values[8]);
c1_2 = _mm_add_ss(c1_2, _mm_mul_ss(a1_2, b1));
_mm_store_ss(&C[(i*56)+17], c1_2);
__m128 c1_3 = _mm_load_ss(&C[(i*56)+32]);
__m128 a1_3 = _mm_load_ss(&values[9]);
c1_3 = _mm_add_ss(c1_3, _mm_mul_ss(a1_3, b1));
_mm_store_ss(&C[(i*56)+32], c1_3);
__m128 c1_4 = _mm_load_ss(&C[(i*56)+53]);
__m128 a1_4 = _mm_load_ss(&values[10]);
c1_4 = _mm_add_ss(c1_4, _mm_mul_ss(a1_4, b1));
_mm_store_ss(&C[(i*56)+53], c1_4);
#else
C[(i*56)+1] += values[6] * B[(i*56)+1];
C[(i*56)+7] += values[7] * B[(i*56)+1];
C[(i*56)+17] += values[8] * B[(i*56)+1];
C[(i*56)+32] += values[9] * B[(i*56)+1];
C[(i*56)+53] += values[10] * B[(i*56)+1];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b2 = _mm_broadcast_ss(&B[(i*56)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b2 = _mm_load_ss(&B[(i*56)+2]);
b2 = _mm_shuffle_ps(b2, b2, 0x00);
#endif
__m128 c2_0 = _mm_load_ss(&C[(i*56)+2]);
__m128 a2_0 = _mm_load_ss(&values[11]);
c2_0 = _mm_add_ss(c2_0, _mm_mul_ss(a2_0, b2));
_mm_store_ss(&C[(i*56)+2], c2_0);
__m128 c2_1 = _mm_load_ss(&C[(i*56)+8]);
__m128 a2_1 = _mm_load_ss(&values[12]);
c2_1 = _mm_add_ss(c2_1, _mm_mul_ss(a2_1, b2));
_mm_store_ss(&C[(i*56)+8], c2_1);
__m128 c2_2 = _mm_load_ss(&C[(i*56)+18]);
__m128 a2_2 = _mm_load_ss(&values[13]);
c2_2 = _mm_add_ss(c2_2, _mm_mul_ss(a2_2, b2));
_mm_store_ss(&C[(i*56)+18], c2_2);
__m128 c2_3 = _mm_load_ss(&C[(i*56)+33]);
__m128 a2_3 = _mm_load_ss(&values[14]);
c2_3 = _mm_add_ss(c2_3, _mm_mul_ss(a2_3, b2));
_mm_store_ss(&C[(i*56)+33], c2_3);
__m128 c2_4 = _mm_load_ss(&C[(i*56)+54]);
__m128 a2_4 = _mm_load_ss(&values[15]);
c2_4 = _mm_add_ss(c2_4, _mm_mul_ss(a2_4, b2));
_mm_store_ss(&C[(i*56)+54], c2_4);
#else
C[(i*56)+2] += values[11] * B[(i*56)+2];
C[(i*56)+8] += values[12] * B[(i*56)+2];
C[(i*56)+18] += values[13] * B[(i*56)+2];
C[(i*56)+33] += values[14] * B[(i*56)+2];
C[(i*56)+54] += values[15] * B[(i*56)+2];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b3 = _mm_broadcast_ss(&B[(i*56)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b3 = _mm_load_ss(&B[(i*56)+3]);
b3 = _mm_shuffle_ps(b3, b3, 0x00);
#endif
__m128 c3_0 = _mm_load_ss(&C[(i*56)+0]);
__m128 a3_0 = _mm_load_ss(&values[16]);
c3_0 = _mm_add_ss(c3_0, _mm_mul_ss(a3_0, b3));
_mm_store_ss(&C[(i*56)+0], c3_0);
__m128 c3_1 = _mm_load_ss(&C[(i*56)+3]);
__m128 a3_1 = _mm_load_ss(&values[17]);
c3_1 = _mm_add_ss(c3_1, _mm_mul_ss(a3_1, b3));
_mm_store_ss(&C[(i*56)+3], c3_1);
__m128 c3_2 = _mm_load_ss(&C[(i*56)+9]);
__m128 a3_2 = _mm_load_ss(&values[18]);
c3_2 = _mm_add_ss(c3_2, _mm_mul_ss(a3_2, b3));
_mm_store_ss(&C[(i*56)+9], c3_2);
__m128 c3_3 = _mm_load_ss(&C[(i*56)+19]);
__m128 a3_3 = _mm_load_ss(&values[19]);
c3_3 = _mm_add_ss(c3_3, _mm_mul_ss(a3_3, b3));
_mm_store_ss(&C[(i*56)+19], c3_3);
__m128 c3_4 = _mm_load_ss(&C[(i*56)+34]);
__m128 a3_4 = _mm_load_ss(&values[20]);
c3_4 = _mm_add_ss(c3_4, _mm_mul_ss(a3_4, b3));
_mm_store_ss(&C[(i*56)+34], c3_4);
__m128 c3_5 = _mm_load_ss(&C[(i*56)+55]);
__m128 a3_5 = _mm_load_ss(&values[21]);
c3_5 = _mm_add_ss(c3_5, _mm_mul_ss(a3_5, b3));
_mm_store_ss(&C[(i*56)+55], c3_5);
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
__m128 b4 = _mm_broadcast_ss(&B[(i*56)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b4 = _mm_load_ss(&B[(i*56)+4]);
b4 = _mm_shuffle_ps(b4, b4, 0x00);
#endif
__m128 c4_0 = _mm_load_ss(&C[(i*56)+4]);
__m128 a4_0 = _mm_load_ss(&values[22]);
c4_0 = _mm_add_ss(c4_0, _mm_mul_ss(a4_0, b4));
_mm_store_ss(&C[(i*56)+4], c4_0);
__m128 c4_1 = _mm_load_ss(&C[(i*56)+14]);
__m128 a4_1 = _mm_load_ss(&values[23]);
c4_1 = _mm_add_ss(c4_1, _mm_mul_ss(a4_1, b4));
_mm_store_ss(&C[(i*56)+14], c4_1);
__m128 c4_2 = _mm_load_ss(&C[(i*56)+29]);
__m128 a4_2 = _mm_load_ss(&values[24]);
c4_2 = _mm_add_ss(c4_2, _mm_mul_ss(a4_2, b4));
_mm_store_ss(&C[(i*56)+29], c4_2);
__m128 c4_3 = _mm_load_ss(&C[(i*56)+50]);
__m128 a4_3 = _mm_load_ss(&values[25]);
c4_3 = _mm_add_ss(c4_3, _mm_mul_ss(a4_3, b4));
_mm_store_ss(&C[(i*56)+50], c4_3);
#else
C[(i*56)+4] += values[22] * B[(i*56)+4];
C[(i*56)+14] += values[23] * B[(i*56)+4];
C[(i*56)+29] += values[24] * B[(i*56)+4];
C[(i*56)+50] += values[25] * B[(i*56)+4];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b5 = _mm_broadcast_ss(&B[(i*56)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b5 = _mm_load_ss(&B[(i*56)+5]);
b5 = _mm_shuffle_ps(b5, b5, 0x00);
#endif
__m128 c5_0 = _mm_load_ss(&C[(i*56)+5]);
__m128 a5_0 = _mm_load_ss(&values[26]);
c5_0 = _mm_add_ss(c5_0, _mm_mul_ss(a5_0, b5));
_mm_store_ss(&C[(i*56)+5], c5_0);
__m128 c5_1 = _mm_load_ss(&C[(i*56)+15]);
__m128 a5_1 = _mm_load_ss(&values[27]);
c5_1 = _mm_add_ss(c5_1, _mm_mul_ss(a5_1, b5));
_mm_store_ss(&C[(i*56)+15], c5_1);
__m128 c5_2 = _mm_load_ss(&C[(i*56)+30]);
__m128 a5_2 = _mm_load_ss(&values[28]);
c5_2 = _mm_add_ss(c5_2, _mm_mul_ss(a5_2, b5));
_mm_store_ss(&C[(i*56)+30], c5_2);
__m128 c5_3 = _mm_load_ss(&C[(i*56)+51]);
__m128 a5_3 = _mm_load_ss(&values[29]);
c5_3 = _mm_add_ss(c5_3, _mm_mul_ss(a5_3, b5));
_mm_store_ss(&C[(i*56)+51], c5_3);
#else
C[(i*56)+5] += values[26] * B[(i*56)+5];
C[(i*56)+15] += values[27] * B[(i*56)+5];
C[(i*56)+30] += values[28] * B[(i*56)+5];
C[(i*56)+51] += values[29] * B[(i*56)+5];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b6 = _mm_broadcast_ss(&B[(i*56)+6]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b6 = _mm_load_ss(&B[(i*56)+6]);
b6 = _mm_shuffle_ps(b6, b6, 0x00);
#endif
__m128 c6_0 = _mm_load_ss(&C[(i*56)+6]);
__m128 a6_0 = _mm_load_ss(&values[30]);
c6_0 = _mm_add_ss(c6_0, _mm_mul_ss(a6_0, b6));
_mm_store_ss(&C[(i*56)+6], c6_0);
__m128 c6_1 = _mm_load_ss(&C[(i*56)+16]);
__m128 a6_1 = _mm_load_ss(&values[31]);
c6_1 = _mm_add_ss(c6_1, _mm_mul_ss(a6_1, b6));
_mm_store_ss(&C[(i*56)+16], c6_1);
__m128 c6_2 = _mm_load_ss(&C[(i*56)+31]);
__m128 a6_2 = _mm_load_ss(&values[32]);
c6_2 = _mm_add_ss(c6_2, _mm_mul_ss(a6_2, b6));
_mm_store_ss(&C[(i*56)+31], c6_2);
__m128 c6_3 = _mm_load_ss(&C[(i*56)+52]);
__m128 a6_3 = _mm_load_ss(&values[33]);
c6_3 = _mm_add_ss(c6_3, _mm_mul_ss(a6_3, b6));
_mm_store_ss(&C[(i*56)+52], c6_3);
#else
C[(i*56)+6] += values[30] * B[(i*56)+6];
C[(i*56)+16] += values[31] * B[(i*56)+6];
C[(i*56)+31] += values[32] * B[(i*56)+6];
C[(i*56)+52] += values[33] * B[(i*56)+6];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b7 = _mm_broadcast_ss(&B[(i*56)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b7 = _mm_load_ss(&B[(i*56)+7]);
b7 = _mm_shuffle_ps(b7, b7, 0x00);
#endif
__m128 c7_0 = _mm_load_ss(&C[(i*56)+1]);
__m128 a7_0 = _mm_load_ss(&values[34]);
c7_0 = _mm_add_ss(c7_0, _mm_mul_ss(a7_0, b7));
_mm_store_ss(&C[(i*56)+1], c7_0);
__m128 c7_1 = _mm_load_ss(&C[(i*56)+7]);
__m128 a7_1 = _mm_load_ss(&values[35]);
c7_1 = _mm_add_ss(c7_1, _mm_mul_ss(a7_1, b7));
_mm_store_ss(&C[(i*56)+7], c7_1);
__m128 c7_2 = _mm_load_ss(&C[(i*56)+17]);
__m128 a7_2 = _mm_load_ss(&values[36]);
c7_2 = _mm_add_ss(c7_2, _mm_mul_ss(a7_2, b7));
_mm_store_ss(&C[(i*56)+17], c7_2);
__m128 c7_3 = _mm_load_ss(&C[(i*56)+32]);
__m128 a7_3 = _mm_load_ss(&values[37]);
c7_3 = _mm_add_ss(c7_3, _mm_mul_ss(a7_3, b7));
_mm_store_ss(&C[(i*56)+32], c7_3);
__m128 c7_4 = _mm_load_ss(&C[(i*56)+53]);
__m128 a7_4 = _mm_load_ss(&values[38]);
c7_4 = _mm_add_ss(c7_4, _mm_mul_ss(a7_4, b7));
_mm_store_ss(&C[(i*56)+53], c7_4);
#else
C[(i*56)+1] += values[34] * B[(i*56)+7];
C[(i*56)+7] += values[35] * B[(i*56)+7];
C[(i*56)+17] += values[36] * B[(i*56)+7];
C[(i*56)+32] += values[37] * B[(i*56)+7];
C[(i*56)+53] += values[38] * B[(i*56)+7];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b8 = _mm_broadcast_ss(&B[(i*56)+8]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b8 = _mm_load_ss(&B[(i*56)+8]);
b8 = _mm_shuffle_ps(b8, b8, 0x00);
#endif
__m128 c8_0 = _mm_load_ss(&C[(i*56)+2]);
__m128 a8_0 = _mm_load_ss(&values[39]);
c8_0 = _mm_add_ss(c8_0, _mm_mul_ss(a8_0, b8));
_mm_store_ss(&C[(i*56)+2], c8_0);
__m128 c8_1 = _mm_load_ss(&C[(i*56)+8]);
__m128 a8_1 = _mm_load_ss(&values[40]);
c8_1 = _mm_add_ss(c8_1, _mm_mul_ss(a8_1, b8));
_mm_store_ss(&C[(i*56)+8], c8_1);
__m128 c8_2 = _mm_load_ss(&C[(i*56)+18]);
__m128 a8_2 = _mm_load_ss(&values[41]);
c8_2 = _mm_add_ss(c8_2, _mm_mul_ss(a8_2, b8));
_mm_store_ss(&C[(i*56)+18], c8_2);
__m128 c8_3 = _mm_load_ss(&C[(i*56)+33]);
__m128 a8_3 = _mm_load_ss(&values[42]);
c8_3 = _mm_add_ss(c8_3, _mm_mul_ss(a8_3, b8));
_mm_store_ss(&C[(i*56)+33], c8_3);
__m128 c8_4 = _mm_load_ss(&C[(i*56)+54]);
__m128 a8_4 = _mm_load_ss(&values[43]);
c8_4 = _mm_add_ss(c8_4, _mm_mul_ss(a8_4, b8));
_mm_store_ss(&C[(i*56)+54], c8_4);
#else
C[(i*56)+2] += values[39] * B[(i*56)+8];
C[(i*56)+8] += values[40] * B[(i*56)+8];
C[(i*56)+18] += values[41] * B[(i*56)+8];
C[(i*56)+33] += values[42] * B[(i*56)+8];
C[(i*56)+54] += values[43] * B[(i*56)+8];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b9 = _mm_broadcast_ss(&B[(i*56)+9]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b9 = _mm_load_ss(&B[(i*56)+9]);
b9 = _mm_shuffle_ps(b9, b9, 0x00);
#endif
__m128 c9_0 = _mm_load_ss(&C[(i*56)+0]);
__m128 a9_0 = _mm_load_ss(&values[44]);
c9_0 = _mm_add_ss(c9_0, _mm_mul_ss(a9_0, b9));
_mm_store_ss(&C[(i*56)+0], c9_0);
__m128 c9_1 = _mm_load_ss(&C[(i*56)+3]);
__m128 a9_1 = _mm_load_ss(&values[45]);
c9_1 = _mm_add_ss(c9_1, _mm_mul_ss(a9_1, b9));
_mm_store_ss(&C[(i*56)+3], c9_1);
__m128 c9_2 = _mm_load_ss(&C[(i*56)+9]);
__m128 a9_2 = _mm_load_ss(&values[46]);
c9_2 = _mm_add_ss(c9_2, _mm_mul_ss(a9_2, b9));
_mm_store_ss(&C[(i*56)+9], c9_2);
__m128 c9_3 = _mm_load_ss(&C[(i*56)+19]);
__m128 a9_3 = _mm_load_ss(&values[47]);
c9_3 = _mm_add_ss(c9_3, _mm_mul_ss(a9_3, b9));
_mm_store_ss(&C[(i*56)+19], c9_3);
__m128 c9_4 = _mm_load_ss(&C[(i*56)+34]);
__m128 a9_4 = _mm_load_ss(&values[48]);
c9_4 = _mm_add_ss(c9_4, _mm_mul_ss(a9_4, b9));
_mm_store_ss(&C[(i*56)+34], c9_4);
__m128 c9_5 = _mm_load_ss(&C[(i*56)+55]);
__m128 a9_5 = _mm_load_ss(&values[49]);
c9_5 = _mm_add_ss(c9_5, _mm_mul_ss(a9_5, b9));
_mm_store_ss(&C[(i*56)+55], c9_5);
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
__m128 b10 = _mm_broadcast_ss(&B[(i*56)+10]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b10 = _mm_load_ss(&B[(i*56)+10]);
b10 = _mm_shuffle_ps(b10, b10, 0x00);
#endif
__m128 c10_0 = _mm_load_ss(&C[(i*56)+10]);
__m128 a10_0 = _mm_load_ss(&values[50]);
c10_0 = _mm_add_ss(c10_0, _mm_mul_ss(a10_0, b10));
_mm_store_ss(&C[(i*56)+10], c10_0);
__m128 c10_1 = _mm_load_ss(&C[(i*56)+25]);
__m128 a10_1 = _mm_load_ss(&values[51]);
c10_1 = _mm_add_ss(c10_1, _mm_mul_ss(a10_1, b10));
_mm_store_ss(&C[(i*56)+25], c10_1);
__m128 c10_2 = _mm_load_ss(&C[(i*56)+46]);
__m128 a10_2 = _mm_load_ss(&values[52]);
c10_2 = _mm_add_ss(c10_2, _mm_mul_ss(a10_2, b10));
_mm_store_ss(&C[(i*56)+46], c10_2);
#else
C[(i*56)+10] += values[50] * B[(i*56)+10];
C[(i*56)+25] += values[51] * B[(i*56)+10];
C[(i*56)+46] += values[52] * B[(i*56)+10];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b11 = _mm_broadcast_ss(&B[(i*56)+11]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b11 = _mm_load_ss(&B[(i*56)+11]);
b11 = _mm_shuffle_ps(b11, b11, 0x00);
#endif
__m128 c11_0 = _mm_load_ss(&C[(i*56)+11]);
__m128 a11_0 = _mm_load_ss(&values[53]);
c11_0 = _mm_add_ss(c11_0, _mm_mul_ss(a11_0, b11));
_mm_store_ss(&C[(i*56)+11], c11_0);
__m128 c11_1 = _mm_load_ss(&C[(i*56)+26]);
__m128 a11_1 = _mm_load_ss(&values[54]);
c11_1 = _mm_add_ss(c11_1, _mm_mul_ss(a11_1, b11));
_mm_store_ss(&C[(i*56)+26], c11_1);
__m128 c11_2 = _mm_load_ss(&C[(i*56)+47]);
__m128 a11_2 = _mm_load_ss(&values[55]);
c11_2 = _mm_add_ss(c11_2, _mm_mul_ss(a11_2, b11));
_mm_store_ss(&C[(i*56)+47], c11_2);
#else
C[(i*56)+11] += values[53] * B[(i*56)+11];
C[(i*56)+26] += values[54] * B[(i*56)+11];
C[(i*56)+47] += values[55] * B[(i*56)+11];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b12 = _mm_broadcast_ss(&B[(i*56)+12]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b12 = _mm_load_ss(&B[(i*56)+12]);
b12 = _mm_shuffle_ps(b12, b12, 0x00);
#endif
__m128 c12_0 = _mm_load_ss(&C[(i*56)+12]);
__m128 a12_0 = _mm_load_ss(&values[56]);
c12_0 = _mm_add_ss(c12_0, _mm_mul_ss(a12_0, b12));
_mm_store_ss(&C[(i*56)+12], c12_0);
__m128 c12_1 = _mm_load_ss(&C[(i*56)+27]);
__m128 a12_1 = _mm_load_ss(&values[57]);
c12_1 = _mm_add_ss(c12_1, _mm_mul_ss(a12_1, b12));
_mm_store_ss(&C[(i*56)+27], c12_1);
__m128 c12_2 = _mm_load_ss(&C[(i*56)+48]);
__m128 a12_2 = _mm_load_ss(&values[58]);
c12_2 = _mm_add_ss(c12_2, _mm_mul_ss(a12_2, b12));
_mm_store_ss(&C[(i*56)+48], c12_2);
#else
C[(i*56)+12] += values[56] * B[(i*56)+12];
C[(i*56)+27] += values[57] * B[(i*56)+12];
C[(i*56)+48] += values[58] * B[(i*56)+12];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b13 = _mm_broadcast_ss(&B[(i*56)+13]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b13 = _mm_load_ss(&B[(i*56)+13]);
b13 = _mm_shuffle_ps(b13, b13, 0x00);
#endif
__m128 c13_0 = _mm_load_ss(&C[(i*56)+13]);
__m128 a13_0 = _mm_load_ss(&values[59]);
c13_0 = _mm_add_ss(c13_0, _mm_mul_ss(a13_0, b13));
_mm_store_ss(&C[(i*56)+13], c13_0);
__m128 c13_1 = _mm_load_ss(&C[(i*56)+28]);
__m128 a13_1 = _mm_load_ss(&values[60]);
c13_1 = _mm_add_ss(c13_1, _mm_mul_ss(a13_1, b13));
_mm_store_ss(&C[(i*56)+28], c13_1);
__m128 c13_2 = _mm_load_ss(&C[(i*56)+49]);
__m128 a13_2 = _mm_load_ss(&values[61]);
c13_2 = _mm_add_ss(c13_2, _mm_mul_ss(a13_2, b13));
_mm_store_ss(&C[(i*56)+49], c13_2);
#else
C[(i*56)+13] += values[59] * B[(i*56)+13];
C[(i*56)+28] += values[60] * B[(i*56)+13];
C[(i*56)+49] += values[61] * B[(i*56)+13];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b14 = _mm_broadcast_ss(&B[(i*56)+14]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b14 = _mm_load_ss(&B[(i*56)+14]);
b14 = _mm_shuffle_ps(b14, b14, 0x00);
#endif
__m128 c14_0 = _mm_load_ss(&C[(i*56)+4]);
__m128 a14_0 = _mm_load_ss(&values[62]);
c14_0 = _mm_add_ss(c14_0, _mm_mul_ss(a14_0, b14));
_mm_store_ss(&C[(i*56)+4], c14_0);
__m128 c14_1 = _mm_load_ss(&C[(i*56)+14]);
__m128 a14_1 = _mm_load_ss(&values[63]);
c14_1 = _mm_add_ss(c14_1, _mm_mul_ss(a14_1, b14));
_mm_store_ss(&C[(i*56)+14], c14_1);
__m128 c14_2 = _mm_load_ss(&C[(i*56)+29]);
__m128 a14_2 = _mm_load_ss(&values[64]);
c14_2 = _mm_add_ss(c14_2, _mm_mul_ss(a14_2, b14));
_mm_store_ss(&C[(i*56)+29], c14_2);
__m128 c14_3 = _mm_load_ss(&C[(i*56)+50]);
__m128 a14_3 = _mm_load_ss(&values[65]);
c14_3 = _mm_add_ss(c14_3, _mm_mul_ss(a14_3, b14));
_mm_store_ss(&C[(i*56)+50], c14_3);
#else
C[(i*56)+4] += values[62] * B[(i*56)+14];
C[(i*56)+14] += values[63] * B[(i*56)+14];
C[(i*56)+29] += values[64] * B[(i*56)+14];
C[(i*56)+50] += values[65] * B[(i*56)+14];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b15 = _mm_broadcast_ss(&B[(i*56)+15]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b15 = _mm_load_ss(&B[(i*56)+15]);
b15 = _mm_shuffle_ps(b15, b15, 0x00);
#endif
__m128 c15_0 = _mm_load_ss(&C[(i*56)+5]);
__m128 a15_0 = _mm_load_ss(&values[66]);
c15_0 = _mm_add_ss(c15_0, _mm_mul_ss(a15_0, b15));
_mm_store_ss(&C[(i*56)+5], c15_0);
__m128 c15_1 = _mm_load_ss(&C[(i*56)+15]);
__m128 a15_1 = _mm_load_ss(&values[67]);
c15_1 = _mm_add_ss(c15_1, _mm_mul_ss(a15_1, b15));
_mm_store_ss(&C[(i*56)+15], c15_1);
__m128 c15_2 = _mm_load_ss(&C[(i*56)+30]);
__m128 a15_2 = _mm_load_ss(&values[68]);
c15_2 = _mm_add_ss(c15_2, _mm_mul_ss(a15_2, b15));
_mm_store_ss(&C[(i*56)+30], c15_2);
__m128 c15_3 = _mm_load_ss(&C[(i*56)+51]);
__m128 a15_3 = _mm_load_ss(&values[69]);
c15_3 = _mm_add_ss(c15_3, _mm_mul_ss(a15_3, b15));
_mm_store_ss(&C[(i*56)+51], c15_3);
#else
C[(i*56)+5] += values[66] * B[(i*56)+15];
C[(i*56)+15] += values[67] * B[(i*56)+15];
C[(i*56)+30] += values[68] * B[(i*56)+15];
C[(i*56)+51] += values[69] * B[(i*56)+15];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b16 = _mm_broadcast_ss(&B[(i*56)+16]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b16 = _mm_load_ss(&B[(i*56)+16]);
b16 = _mm_shuffle_ps(b16, b16, 0x00);
#endif
__m128 c16_0 = _mm_load_ss(&C[(i*56)+6]);
__m128 a16_0 = _mm_load_ss(&values[70]);
c16_0 = _mm_add_ss(c16_0, _mm_mul_ss(a16_0, b16));
_mm_store_ss(&C[(i*56)+6], c16_0);
__m128 c16_1 = _mm_load_ss(&C[(i*56)+16]);
__m128 a16_1 = _mm_load_ss(&values[71]);
c16_1 = _mm_add_ss(c16_1, _mm_mul_ss(a16_1, b16));
_mm_store_ss(&C[(i*56)+16], c16_1);
__m128 c16_2 = _mm_load_ss(&C[(i*56)+31]);
__m128 a16_2 = _mm_load_ss(&values[72]);
c16_2 = _mm_add_ss(c16_2, _mm_mul_ss(a16_2, b16));
_mm_store_ss(&C[(i*56)+31], c16_2);
__m128 c16_3 = _mm_load_ss(&C[(i*56)+52]);
__m128 a16_3 = _mm_load_ss(&values[73]);
c16_3 = _mm_add_ss(c16_3, _mm_mul_ss(a16_3, b16));
_mm_store_ss(&C[(i*56)+52], c16_3);
#else
C[(i*56)+6] += values[70] * B[(i*56)+16];
C[(i*56)+16] += values[71] * B[(i*56)+16];
C[(i*56)+31] += values[72] * B[(i*56)+16];
C[(i*56)+52] += values[73] * B[(i*56)+16];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b17 = _mm_broadcast_ss(&B[(i*56)+17]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b17 = _mm_load_ss(&B[(i*56)+17]);
b17 = _mm_shuffle_ps(b17, b17, 0x00);
#endif
__m128 c17_0 = _mm_load_ss(&C[(i*56)+1]);
__m128 a17_0 = _mm_load_ss(&values[74]);
c17_0 = _mm_add_ss(c17_0, _mm_mul_ss(a17_0, b17));
_mm_store_ss(&C[(i*56)+1], c17_0);
__m128 c17_1 = _mm_load_ss(&C[(i*56)+7]);
__m128 a17_1 = _mm_load_ss(&values[75]);
c17_1 = _mm_add_ss(c17_1, _mm_mul_ss(a17_1, b17));
_mm_store_ss(&C[(i*56)+7], c17_1);
__m128 c17_2 = _mm_load_ss(&C[(i*56)+17]);
__m128 a17_2 = _mm_load_ss(&values[76]);
c17_2 = _mm_add_ss(c17_2, _mm_mul_ss(a17_2, b17));
_mm_store_ss(&C[(i*56)+17], c17_2);
__m128 c17_3 = _mm_load_ss(&C[(i*56)+32]);
__m128 a17_3 = _mm_load_ss(&values[77]);
c17_3 = _mm_add_ss(c17_3, _mm_mul_ss(a17_3, b17));
_mm_store_ss(&C[(i*56)+32], c17_3);
__m128 c17_4 = _mm_load_ss(&C[(i*56)+53]);
__m128 a17_4 = _mm_load_ss(&values[78]);
c17_4 = _mm_add_ss(c17_4, _mm_mul_ss(a17_4, b17));
_mm_store_ss(&C[(i*56)+53], c17_4);
#else
C[(i*56)+1] += values[74] * B[(i*56)+17];
C[(i*56)+7] += values[75] * B[(i*56)+17];
C[(i*56)+17] += values[76] * B[(i*56)+17];
C[(i*56)+32] += values[77] * B[(i*56)+17];
C[(i*56)+53] += values[78] * B[(i*56)+17];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b18 = _mm_broadcast_ss(&B[(i*56)+18]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b18 = _mm_load_ss(&B[(i*56)+18]);
b18 = _mm_shuffle_ps(b18, b18, 0x00);
#endif
__m128 c18_0 = _mm_load_ss(&C[(i*56)+2]);
__m128 a18_0 = _mm_load_ss(&values[79]);
c18_0 = _mm_add_ss(c18_0, _mm_mul_ss(a18_0, b18));
_mm_store_ss(&C[(i*56)+2], c18_0);
__m128 c18_1 = _mm_load_ss(&C[(i*56)+8]);
__m128 a18_1 = _mm_load_ss(&values[80]);
c18_1 = _mm_add_ss(c18_1, _mm_mul_ss(a18_1, b18));
_mm_store_ss(&C[(i*56)+8], c18_1);
__m128 c18_2 = _mm_load_ss(&C[(i*56)+18]);
__m128 a18_2 = _mm_load_ss(&values[81]);
c18_2 = _mm_add_ss(c18_2, _mm_mul_ss(a18_2, b18));
_mm_store_ss(&C[(i*56)+18], c18_2);
__m128 c18_3 = _mm_load_ss(&C[(i*56)+33]);
__m128 a18_3 = _mm_load_ss(&values[82]);
c18_3 = _mm_add_ss(c18_3, _mm_mul_ss(a18_3, b18));
_mm_store_ss(&C[(i*56)+33], c18_3);
__m128 c18_4 = _mm_load_ss(&C[(i*56)+54]);
__m128 a18_4 = _mm_load_ss(&values[83]);
c18_4 = _mm_add_ss(c18_4, _mm_mul_ss(a18_4, b18));
_mm_store_ss(&C[(i*56)+54], c18_4);
#else
C[(i*56)+2] += values[79] * B[(i*56)+18];
C[(i*56)+8] += values[80] * B[(i*56)+18];
C[(i*56)+18] += values[81] * B[(i*56)+18];
C[(i*56)+33] += values[82] * B[(i*56)+18];
C[(i*56)+54] += values[83] * B[(i*56)+18];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b19 = _mm_broadcast_ss(&B[(i*56)+19]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b19 = _mm_load_ss(&B[(i*56)+19]);
b19 = _mm_shuffle_ps(b19, b19, 0x00);
#endif
__m128 c19_0 = _mm_load_ss(&C[(i*56)+0]);
__m128 a19_0 = _mm_load_ss(&values[84]);
c19_0 = _mm_add_ss(c19_0, _mm_mul_ss(a19_0, b19));
_mm_store_ss(&C[(i*56)+0], c19_0);
__m128 c19_1 = _mm_load_ss(&C[(i*56)+3]);
__m128 a19_1 = _mm_load_ss(&values[85]);
c19_1 = _mm_add_ss(c19_1, _mm_mul_ss(a19_1, b19));
_mm_store_ss(&C[(i*56)+3], c19_1);
__m128 c19_2 = _mm_load_ss(&C[(i*56)+9]);
__m128 a19_2 = _mm_load_ss(&values[86]);
c19_2 = _mm_add_ss(c19_2, _mm_mul_ss(a19_2, b19));
_mm_store_ss(&C[(i*56)+9], c19_2);
__m128 c19_3 = _mm_load_ss(&C[(i*56)+19]);
__m128 a19_3 = _mm_load_ss(&values[87]);
c19_3 = _mm_add_ss(c19_3, _mm_mul_ss(a19_3, b19));
_mm_store_ss(&C[(i*56)+19], c19_3);
__m128 c19_4 = _mm_load_ss(&C[(i*56)+34]);
__m128 a19_4 = _mm_load_ss(&values[88]);
c19_4 = _mm_add_ss(c19_4, _mm_mul_ss(a19_4, b19));
_mm_store_ss(&C[(i*56)+34], c19_4);
__m128 c19_5 = _mm_load_ss(&C[(i*56)+55]);
__m128 a19_5 = _mm_load_ss(&values[89]);
c19_5 = _mm_add_ss(c19_5, _mm_mul_ss(a19_5, b19));
_mm_store_ss(&C[(i*56)+55], c19_5);
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
__m128 b20 = _mm_broadcast_ss(&B[(i*56)+20]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b20 = _mm_load_ss(&B[(i*56)+20]);
b20 = _mm_shuffle_ps(b20, b20, 0x00);
#endif
__m128 c20_0 = _mm_load_ss(&C[(i*56)+20]);
__m128 a20_0 = _mm_load_ss(&values[90]);
c20_0 = _mm_add_ss(c20_0, _mm_mul_ss(a20_0, b20));
_mm_store_ss(&C[(i*56)+20], c20_0);
__m128 c20_1 = _mm_load_ss(&C[(i*56)+41]);
__m128 a20_1 = _mm_load_ss(&values[91]);
c20_1 = _mm_add_ss(c20_1, _mm_mul_ss(a20_1, b20));
_mm_store_ss(&C[(i*56)+41], c20_1);
#else
C[(i*56)+20] += values[90] * B[(i*56)+20];
C[(i*56)+41] += values[91] * B[(i*56)+20];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b21 = _mm_broadcast_ss(&B[(i*56)+21]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b21 = _mm_load_ss(&B[(i*56)+21]);
b21 = _mm_shuffle_ps(b21, b21, 0x00);
#endif
__m128 c21_0 = _mm_load_ss(&C[(i*56)+21]);
__m128 a21_0 = _mm_load_ss(&values[92]);
c21_0 = _mm_add_ss(c21_0, _mm_mul_ss(a21_0, b21));
_mm_store_ss(&C[(i*56)+21], c21_0);
__m128 c21_1 = _mm_load_ss(&C[(i*56)+42]);
__m128 a21_1 = _mm_load_ss(&values[93]);
c21_1 = _mm_add_ss(c21_1, _mm_mul_ss(a21_1, b21));
_mm_store_ss(&C[(i*56)+42], c21_1);
#else
C[(i*56)+21] += values[92] * B[(i*56)+21];
C[(i*56)+42] += values[93] * B[(i*56)+21];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b22 = _mm_broadcast_ss(&B[(i*56)+22]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b22 = _mm_load_ss(&B[(i*56)+22]);
b22 = _mm_shuffle_ps(b22, b22, 0x00);
#endif
__m128 c22_0 = _mm_load_ss(&C[(i*56)+22]);
__m128 a22_0 = _mm_load_ss(&values[94]);
c22_0 = _mm_add_ss(c22_0, _mm_mul_ss(a22_0, b22));
_mm_store_ss(&C[(i*56)+22], c22_0);
__m128 c22_1 = _mm_load_ss(&C[(i*56)+43]);
__m128 a22_1 = _mm_load_ss(&values[95]);
c22_1 = _mm_add_ss(c22_1, _mm_mul_ss(a22_1, b22));
_mm_store_ss(&C[(i*56)+43], c22_1);
#else
C[(i*56)+22] += values[94] * B[(i*56)+22];
C[(i*56)+43] += values[95] * B[(i*56)+22];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b23 = _mm_broadcast_ss(&B[(i*56)+23]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b23 = _mm_load_ss(&B[(i*56)+23]);
b23 = _mm_shuffle_ps(b23, b23, 0x00);
#endif
__m128 c23_0 = _mm_load_ss(&C[(i*56)+23]);
__m128 a23_0 = _mm_load_ss(&values[96]);
c23_0 = _mm_add_ss(c23_0, _mm_mul_ss(a23_0, b23));
_mm_store_ss(&C[(i*56)+23], c23_0);
__m128 c23_1 = _mm_load_ss(&C[(i*56)+44]);
__m128 a23_1 = _mm_load_ss(&values[97]);
c23_1 = _mm_add_ss(c23_1, _mm_mul_ss(a23_1, b23));
_mm_store_ss(&C[(i*56)+44], c23_1);
#else
C[(i*56)+23] += values[96] * B[(i*56)+23];
C[(i*56)+44] += values[97] * B[(i*56)+23];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b24 = _mm_broadcast_ss(&B[(i*56)+24]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b24 = _mm_load_ss(&B[(i*56)+24]);
b24 = _mm_shuffle_ps(b24, b24, 0x00);
#endif
__m128 c24_0 = _mm_load_ss(&C[(i*56)+24]);
__m128 a24_0 = _mm_load_ss(&values[98]);
c24_0 = _mm_add_ss(c24_0, _mm_mul_ss(a24_0, b24));
_mm_store_ss(&C[(i*56)+24], c24_0);
__m128 c24_1 = _mm_load_ss(&C[(i*56)+45]);
__m128 a24_1 = _mm_load_ss(&values[99]);
c24_1 = _mm_add_ss(c24_1, _mm_mul_ss(a24_1, b24));
_mm_store_ss(&C[(i*56)+45], c24_1);
#else
C[(i*56)+24] += values[98] * B[(i*56)+24];
C[(i*56)+45] += values[99] * B[(i*56)+24];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b25 = _mm_broadcast_ss(&B[(i*56)+25]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b25 = _mm_load_ss(&B[(i*56)+25]);
b25 = _mm_shuffle_ps(b25, b25, 0x00);
#endif
__m128 c25_0 = _mm_load_ss(&C[(i*56)+10]);
__m128 a25_0 = _mm_load_ss(&values[100]);
c25_0 = _mm_add_ss(c25_0, _mm_mul_ss(a25_0, b25));
_mm_store_ss(&C[(i*56)+10], c25_0);
__m128 c25_1 = _mm_load_ss(&C[(i*56)+25]);
__m128 a25_1 = _mm_load_ss(&values[101]);
c25_1 = _mm_add_ss(c25_1, _mm_mul_ss(a25_1, b25));
_mm_store_ss(&C[(i*56)+25], c25_1);
__m128 c25_2 = _mm_load_ss(&C[(i*56)+46]);
__m128 a25_2 = _mm_load_ss(&values[102]);
c25_2 = _mm_add_ss(c25_2, _mm_mul_ss(a25_2, b25));
_mm_store_ss(&C[(i*56)+46], c25_2);
#else
C[(i*56)+10] += values[100] * B[(i*56)+25];
C[(i*56)+25] += values[101] * B[(i*56)+25];
C[(i*56)+46] += values[102] * B[(i*56)+25];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b26 = _mm_broadcast_ss(&B[(i*56)+26]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b26 = _mm_load_ss(&B[(i*56)+26]);
b26 = _mm_shuffle_ps(b26, b26, 0x00);
#endif
__m128 c26_0 = _mm_load_ss(&C[(i*56)+11]);
__m128 a26_0 = _mm_load_ss(&values[103]);
c26_0 = _mm_add_ss(c26_0, _mm_mul_ss(a26_0, b26));
_mm_store_ss(&C[(i*56)+11], c26_0);
__m128 c26_1 = _mm_load_ss(&C[(i*56)+26]);
__m128 a26_1 = _mm_load_ss(&values[104]);
c26_1 = _mm_add_ss(c26_1, _mm_mul_ss(a26_1, b26));
_mm_store_ss(&C[(i*56)+26], c26_1);
__m128 c26_2 = _mm_load_ss(&C[(i*56)+47]);
__m128 a26_2 = _mm_load_ss(&values[105]);
c26_2 = _mm_add_ss(c26_2, _mm_mul_ss(a26_2, b26));
_mm_store_ss(&C[(i*56)+47], c26_2);
#else
C[(i*56)+11] += values[103] * B[(i*56)+26];
C[(i*56)+26] += values[104] * B[(i*56)+26];
C[(i*56)+47] += values[105] * B[(i*56)+26];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b27 = _mm_broadcast_ss(&B[(i*56)+27]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b27 = _mm_load_ss(&B[(i*56)+27]);
b27 = _mm_shuffle_ps(b27, b27, 0x00);
#endif
__m128 c27_0 = _mm_load_ss(&C[(i*56)+12]);
__m128 a27_0 = _mm_load_ss(&values[106]);
c27_0 = _mm_add_ss(c27_0, _mm_mul_ss(a27_0, b27));
_mm_store_ss(&C[(i*56)+12], c27_0);
__m128 c27_1 = _mm_load_ss(&C[(i*56)+27]);
__m128 a27_1 = _mm_load_ss(&values[107]);
c27_1 = _mm_add_ss(c27_1, _mm_mul_ss(a27_1, b27));
_mm_store_ss(&C[(i*56)+27], c27_1);
__m128 c27_2 = _mm_load_ss(&C[(i*56)+48]);
__m128 a27_2 = _mm_load_ss(&values[108]);
c27_2 = _mm_add_ss(c27_2, _mm_mul_ss(a27_2, b27));
_mm_store_ss(&C[(i*56)+48], c27_2);
#else
C[(i*56)+12] += values[106] * B[(i*56)+27];
C[(i*56)+27] += values[107] * B[(i*56)+27];
C[(i*56)+48] += values[108] * B[(i*56)+27];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b28 = _mm_broadcast_ss(&B[(i*56)+28]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b28 = _mm_load_ss(&B[(i*56)+28]);
b28 = _mm_shuffle_ps(b28, b28, 0x00);
#endif
__m128 c28_0 = _mm_load_ss(&C[(i*56)+13]);
__m128 a28_0 = _mm_load_ss(&values[109]);
c28_0 = _mm_add_ss(c28_0, _mm_mul_ss(a28_0, b28));
_mm_store_ss(&C[(i*56)+13], c28_0);
__m128 c28_1 = _mm_load_ss(&C[(i*56)+28]);
__m128 a28_1 = _mm_load_ss(&values[110]);
c28_1 = _mm_add_ss(c28_1, _mm_mul_ss(a28_1, b28));
_mm_store_ss(&C[(i*56)+28], c28_1);
__m128 c28_2 = _mm_load_ss(&C[(i*56)+49]);
__m128 a28_2 = _mm_load_ss(&values[111]);
c28_2 = _mm_add_ss(c28_2, _mm_mul_ss(a28_2, b28));
_mm_store_ss(&C[(i*56)+49], c28_2);
#else
C[(i*56)+13] += values[109] * B[(i*56)+28];
C[(i*56)+28] += values[110] * B[(i*56)+28];
C[(i*56)+49] += values[111] * B[(i*56)+28];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b29 = _mm_broadcast_ss(&B[(i*56)+29]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b29 = _mm_load_ss(&B[(i*56)+29]);
b29 = _mm_shuffle_ps(b29, b29, 0x00);
#endif
__m128 c29_0 = _mm_load_ss(&C[(i*56)+4]);
__m128 a29_0 = _mm_load_ss(&values[112]);
c29_0 = _mm_add_ss(c29_0, _mm_mul_ss(a29_0, b29));
_mm_store_ss(&C[(i*56)+4], c29_0);
__m128 c29_1 = _mm_load_ss(&C[(i*56)+14]);
__m128 a29_1 = _mm_load_ss(&values[113]);
c29_1 = _mm_add_ss(c29_1, _mm_mul_ss(a29_1, b29));
_mm_store_ss(&C[(i*56)+14], c29_1);
__m128 c29_2 = _mm_load_ss(&C[(i*56)+29]);
__m128 a29_2 = _mm_load_ss(&values[114]);
c29_2 = _mm_add_ss(c29_2, _mm_mul_ss(a29_2, b29));
_mm_store_ss(&C[(i*56)+29], c29_2);
__m128 c29_3 = _mm_load_ss(&C[(i*56)+50]);
__m128 a29_3 = _mm_load_ss(&values[115]);
c29_3 = _mm_add_ss(c29_3, _mm_mul_ss(a29_3, b29));
_mm_store_ss(&C[(i*56)+50], c29_3);
#else
C[(i*56)+4] += values[112] * B[(i*56)+29];
C[(i*56)+14] += values[113] * B[(i*56)+29];
C[(i*56)+29] += values[114] * B[(i*56)+29];
C[(i*56)+50] += values[115] * B[(i*56)+29];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b30 = _mm_broadcast_ss(&B[(i*56)+30]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b30 = _mm_load_ss(&B[(i*56)+30]);
b30 = _mm_shuffle_ps(b30, b30, 0x00);
#endif
__m128 c30_0 = _mm_load_ss(&C[(i*56)+5]);
__m128 a30_0 = _mm_load_ss(&values[116]);
c30_0 = _mm_add_ss(c30_0, _mm_mul_ss(a30_0, b30));
_mm_store_ss(&C[(i*56)+5], c30_0);
__m128 c30_1 = _mm_load_ss(&C[(i*56)+15]);
__m128 a30_1 = _mm_load_ss(&values[117]);
c30_1 = _mm_add_ss(c30_1, _mm_mul_ss(a30_1, b30));
_mm_store_ss(&C[(i*56)+15], c30_1);
__m128 c30_2 = _mm_load_ss(&C[(i*56)+30]);
__m128 a30_2 = _mm_load_ss(&values[118]);
c30_2 = _mm_add_ss(c30_2, _mm_mul_ss(a30_2, b30));
_mm_store_ss(&C[(i*56)+30], c30_2);
__m128 c30_3 = _mm_load_ss(&C[(i*56)+51]);
__m128 a30_3 = _mm_load_ss(&values[119]);
c30_3 = _mm_add_ss(c30_3, _mm_mul_ss(a30_3, b30));
_mm_store_ss(&C[(i*56)+51], c30_3);
#else
C[(i*56)+5] += values[116] * B[(i*56)+30];
C[(i*56)+15] += values[117] * B[(i*56)+30];
C[(i*56)+30] += values[118] * B[(i*56)+30];
C[(i*56)+51] += values[119] * B[(i*56)+30];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b31 = _mm_broadcast_ss(&B[(i*56)+31]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b31 = _mm_load_ss(&B[(i*56)+31]);
b31 = _mm_shuffle_ps(b31, b31, 0x00);
#endif
__m128 c31_0 = _mm_load_ss(&C[(i*56)+6]);
__m128 a31_0 = _mm_load_ss(&values[120]);
c31_0 = _mm_add_ss(c31_0, _mm_mul_ss(a31_0, b31));
_mm_store_ss(&C[(i*56)+6], c31_0);
__m128 c31_1 = _mm_load_ss(&C[(i*56)+16]);
__m128 a31_1 = _mm_load_ss(&values[121]);
c31_1 = _mm_add_ss(c31_1, _mm_mul_ss(a31_1, b31));
_mm_store_ss(&C[(i*56)+16], c31_1);
__m128 c31_2 = _mm_load_ss(&C[(i*56)+31]);
__m128 a31_2 = _mm_load_ss(&values[122]);
c31_2 = _mm_add_ss(c31_2, _mm_mul_ss(a31_2, b31));
_mm_store_ss(&C[(i*56)+31], c31_2);
__m128 c31_3 = _mm_load_ss(&C[(i*56)+52]);
__m128 a31_3 = _mm_load_ss(&values[123]);
c31_3 = _mm_add_ss(c31_3, _mm_mul_ss(a31_3, b31));
_mm_store_ss(&C[(i*56)+52], c31_3);
#else
C[(i*56)+6] += values[120] * B[(i*56)+31];
C[(i*56)+16] += values[121] * B[(i*56)+31];
C[(i*56)+31] += values[122] * B[(i*56)+31];
C[(i*56)+52] += values[123] * B[(i*56)+31];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b32 = _mm_broadcast_ss(&B[(i*56)+32]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b32 = _mm_load_ss(&B[(i*56)+32]);
b32 = _mm_shuffle_ps(b32, b32, 0x00);
#endif
__m128 c32_0 = _mm_load_ss(&C[(i*56)+1]);
__m128 a32_0 = _mm_load_ss(&values[124]);
c32_0 = _mm_add_ss(c32_0, _mm_mul_ss(a32_0, b32));
_mm_store_ss(&C[(i*56)+1], c32_0);
__m128 c32_1 = _mm_load_ss(&C[(i*56)+7]);
__m128 a32_1 = _mm_load_ss(&values[125]);
c32_1 = _mm_add_ss(c32_1, _mm_mul_ss(a32_1, b32));
_mm_store_ss(&C[(i*56)+7], c32_1);
__m128 c32_2 = _mm_load_ss(&C[(i*56)+17]);
__m128 a32_2 = _mm_load_ss(&values[126]);
c32_2 = _mm_add_ss(c32_2, _mm_mul_ss(a32_2, b32));
_mm_store_ss(&C[(i*56)+17], c32_2);
__m128 c32_3 = _mm_load_ss(&C[(i*56)+32]);
__m128 a32_3 = _mm_load_ss(&values[127]);
c32_3 = _mm_add_ss(c32_3, _mm_mul_ss(a32_3, b32));
_mm_store_ss(&C[(i*56)+32], c32_3);
__m128 c32_4 = _mm_load_ss(&C[(i*56)+53]);
__m128 a32_4 = _mm_load_ss(&values[128]);
c32_4 = _mm_add_ss(c32_4, _mm_mul_ss(a32_4, b32));
_mm_store_ss(&C[(i*56)+53], c32_4);
#else
C[(i*56)+1] += values[124] * B[(i*56)+32];
C[(i*56)+7] += values[125] * B[(i*56)+32];
C[(i*56)+17] += values[126] * B[(i*56)+32];
C[(i*56)+32] += values[127] * B[(i*56)+32];
C[(i*56)+53] += values[128] * B[(i*56)+32];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b33 = _mm_broadcast_ss(&B[(i*56)+33]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b33 = _mm_load_ss(&B[(i*56)+33]);
b33 = _mm_shuffle_ps(b33, b33, 0x00);
#endif
__m128 c33_0 = _mm_load_ss(&C[(i*56)+2]);
__m128 a33_0 = _mm_load_ss(&values[129]);
c33_0 = _mm_add_ss(c33_0, _mm_mul_ss(a33_0, b33));
_mm_store_ss(&C[(i*56)+2], c33_0);
__m128 c33_1 = _mm_load_ss(&C[(i*56)+8]);
__m128 a33_1 = _mm_load_ss(&values[130]);
c33_1 = _mm_add_ss(c33_1, _mm_mul_ss(a33_1, b33));
_mm_store_ss(&C[(i*56)+8], c33_1);
__m128 c33_2 = _mm_load_ss(&C[(i*56)+18]);
__m128 a33_2 = _mm_load_ss(&values[131]);
c33_2 = _mm_add_ss(c33_2, _mm_mul_ss(a33_2, b33));
_mm_store_ss(&C[(i*56)+18], c33_2);
__m128 c33_3 = _mm_load_ss(&C[(i*56)+33]);
__m128 a33_3 = _mm_load_ss(&values[132]);
c33_3 = _mm_add_ss(c33_3, _mm_mul_ss(a33_3, b33));
_mm_store_ss(&C[(i*56)+33], c33_3);
__m128 c33_4 = _mm_load_ss(&C[(i*56)+54]);
__m128 a33_4 = _mm_load_ss(&values[133]);
c33_4 = _mm_add_ss(c33_4, _mm_mul_ss(a33_4, b33));
_mm_store_ss(&C[(i*56)+54], c33_4);
#else
C[(i*56)+2] += values[129] * B[(i*56)+33];
C[(i*56)+8] += values[130] * B[(i*56)+33];
C[(i*56)+18] += values[131] * B[(i*56)+33];
C[(i*56)+33] += values[132] * B[(i*56)+33];
C[(i*56)+54] += values[133] * B[(i*56)+33];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b34 = _mm_broadcast_ss(&B[(i*56)+34]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b34 = _mm_load_ss(&B[(i*56)+34]);
b34 = _mm_shuffle_ps(b34, b34, 0x00);
#endif
__m128 c34_0 = _mm_load_ss(&C[(i*56)+0]);
__m128 a34_0 = _mm_load_ss(&values[134]);
c34_0 = _mm_add_ss(c34_0, _mm_mul_ss(a34_0, b34));
_mm_store_ss(&C[(i*56)+0], c34_0);
__m128 c34_1 = _mm_load_ss(&C[(i*56)+3]);
__m128 a34_1 = _mm_load_ss(&values[135]);
c34_1 = _mm_add_ss(c34_1, _mm_mul_ss(a34_1, b34));
_mm_store_ss(&C[(i*56)+3], c34_1);
__m128 c34_2 = _mm_load_ss(&C[(i*56)+9]);
__m128 a34_2 = _mm_load_ss(&values[136]);
c34_2 = _mm_add_ss(c34_2, _mm_mul_ss(a34_2, b34));
_mm_store_ss(&C[(i*56)+9], c34_2);
__m128 c34_3 = _mm_load_ss(&C[(i*56)+19]);
__m128 a34_3 = _mm_load_ss(&values[137]);
c34_3 = _mm_add_ss(c34_3, _mm_mul_ss(a34_3, b34));
_mm_store_ss(&C[(i*56)+19], c34_3);
__m128 c34_4 = _mm_load_ss(&C[(i*56)+34]);
__m128 a34_4 = _mm_load_ss(&values[138]);
c34_4 = _mm_add_ss(c34_4, _mm_mul_ss(a34_4, b34));
_mm_store_ss(&C[(i*56)+34], c34_4);
__m128 c34_5 = _mm_load_ss(&C[(i*56)+55]);
__m128 a34_5 = _mm_load_ss(&values[139]);
c34_5 = _mm_add_ss(c34_5, _mm_mul_ss(a34_5, b34));
_mm_store_ss(&C[(i*56)+55], c34_5);
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
__m128 b35 = _mm_broadcast_ss(&B[(i*56)+35]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b35 = _mm_load_ss(&B[(i*56)+35]);
b35 = _mm_shuffle_ps(b35, b35, 0x00);
#endif
__m128 c35_0 = _mm_load_ss(&C[(i*56)+35]);
__m128 a35_0 = _mm_load_ss(&values[140]);
c35_0 = _mm_add_ss(c35_0, _mm_mul_ss(a35_0, b35));
_mm_store_ss(&C[(i*56)+35], c35_0);
#else
C[(i*56)+35] += values[140] * B[(i*56)+35];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b36 = _mm_broadcast_ss(&B[(i*56)+36]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b36 = _mm_load_ss(&B[(i*56)+36]);
b36 = _mm_shuffle_ps(b36, b36, 0x00);
#endif
__m128 c36_0 = _mm_load_ss(&C[(i*56)+36]);
__m128 a36_0 = _mm_load_ss(&values[141]);
c36_0 = _mm_add_ss(c36_0, _mm_mul_ss(a36_0, b36));
_mm_store_ss(&C[(i*56)+36], c36_0);
#else
C[(i*56)+36] += values[141] * B[(i*56)+36];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b37 = _mm_broadcast_ss(&B[(i*56)+37]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b37 = _mm_load_ss(&B[(i*56)+37]);
b37 = _mm_shuffle_ps(b37, b37, 0x00);
#endif
__m128 c37_0 = _mm_load_ss(&C[(i*56)+37]);
__m128 a37_0 = _mm_load_ss(&values[142]);
c37_0 = _mm_add_ss(c37_0, _mm_mul_ss(a37_0, b37));
_mm_store_ss(&C[(i*56)+37], c37_0);
#else
C[(i*56)+37] += values[142] * B[(i*56)+37];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b38 = _mm_broadcast_ss(&B[(i*56)+38]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b38 = _mm_load_ss(&B[(i*56)+38]);
b38 = _mm_shuffle_ps(b38, b38, 0x00);
#endif
__m128 c38_0 = _mm_load_ss(&C[(i*56)+38]);
__m128 a38_0 = _mm_load_ss(&values[143]);
c38_0 = _mm_add_ss(c38_0, _mm_mul_ss(a38_0, b38));
_mm_store_ss(&C[(i*56)+38], c38_0);
#else
C[(i*56)+38] += values[143] * B[(i*56)+38];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b39 = _mm_broadcast_ss(&B[(i*56)+39]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b39 = _mm_load_ss(&B[(i*56)+39]);
b39 = _mm_shuffle_ps(b39, b39, 0x00);
#endif
__m128 c39_0 = _mm_load_ss(&C[(i*56)+39]);
__m128 a39_0 = _mm_load_ss(&values[144]);
c39_0 = _mm_add_ss(c39_0, _mm_mul_ss(a39_0, b39));
_mm_store_ss(&C[(i*56)+39], c39_0);
#else
C[(i*56)+39] += values[144] * B[(i*56)+39];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b40 = _mm_broadcast_ss(&B[(i*56)+40]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b40 = _mm_load_ss(&B[(i*56)+40]);
b40 = _mm_shuffle_ps(b40, b40, 0x00);
#endif
__m128 c40_0 = _mm_load_ss(&C[(i*56)+40]);
__m128 a40_0 = _mm_load_ss(&values[145]);
c40_0 = _mm_add_ss(c40_0, _mm_mul_ss(a40_0, b40));
_mm_store_ss(&C[(i*56)+40], c40_0);
#else
C[(i*56)+40] += values[145] * B[(i*56)+40];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b41 = _mm_broadcast_ss(&B[(i*56)+41]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b41 = _mm_load_ss(&B[(i*56)+41]);
b41 = _mm_shuffle_ps(b41, b41, 0x00);
#endif
__m128 c41_0 = _mm_load_ss(&C[(i*56)+20]);
__m128 a41_0 = _mm_load_ss(&values[146]);
c41_0 = _mm_add_ss(c41_0, _mm_mul_ss(a41_0, b41));
_mm_store_ss(&C[(i*56)+20], c41_0);
__m128 c41_1 = _mm_load_ss(&C[(i*56)+41]);
__m128 a41_1 = _mm_load_ss(&values[147]);
c41_1 = _mm_add_ss(c41_1, _mm_mul_ss(a41_1, b41));
_mm_store_ss(&C[(i*56)+41], c41_1);
#else
C[(i*56)+20] += values[146] * B[(i*56)+41];
C[(i*56)+41] += values[147] * B[(i*56)+41];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b42 = _mm_broadcast_ss(&B[(i*56)+42]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b42 = _mm_load_ss(&B[(i*56)+42]);
b42 = _mm_shuffle_ps(b42, b42, 0x00);
#endif
__m128 c42_0 = _mm_load_ss(&C[(i*56)+21]);
__m128 a42_0 = _mm_load_ss(&values[148]);
c42_0 = _mm_add_ss(c42_0, _mm_mul_ss(a42_0, b42));
_mm_store_ss(&C[(i*56)+21], c42_0);
__m128 c42_1 = _mm_load_ss(&C[(i*56)+42]);
__m128 a42_1 = _mm_load_ss(&values[149]);
c42_1 = _mm_add_ss(c42_1, _mm_mul_ss(a42_1, b42));
_mm_store_ss(&C[(i*56)+42], c42_1);
#else
C[(i*56)+21] += values[148] * B[(i*56)+42];
C[(i*56)+42] += values[149] * B[(i*56)+42];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b43 = _mm_broadcast_ss(&B[(i*56)+43]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b43 = _mm_load_ss(&B[(i*56)+43]);
b43 = _mm_shuffle_ps(b43, b43, 0x00);
#endif
__m128 c43_0 = _mm_load_ss(&C[(i*56)+22]);
__m128 a43_0 = _mm_load_ss(&values[150]);
c43_0 = _mm_add_ss(c43_0, _mm_mul_ss(a43_0, b43));
_mm_store_ss(&C[(i*56)+22], c43_0);
__m128 c43_1 = _mm_load_ss(&C[(i*56)+43]);
__m128 a43_1 = _mm_load_ss(&values[151]);
c43_1 = _mm_add_ss(c43_1, _mm_mul_ss(a43_1, b43));
_mm_store_ss(&C[(i*56)+43], c43_1);
#else
C[(i*56)+22] += values[150] * B[(i*56)+43];
C[(i*56)+43] += values[151] * B[(i*56)+43];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b44 = _mm_broadcast_ss(&B[(i*56)+44]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b44 = _mm_load_ss(&B[(i*56)+44]);
b44 = _mm_shuffle_ps(b44, b44, 0x00);
#endif
__m128 c44_0 = _mm_load_ss(&C[(i*56)+23]);
__m128 a44_0 = _mm_load_ss(&values[152]);
c44_0 = _mm_add_ss(c44_0, _mm_mul_ss(a44_0, b44));
_mm_store_ss(&C[(i*56)+23], c44_0);
__m128 c44_1 = _mm_load_ss(&C[(i*56)+44]);
__m128 a44_1 = _mm_load_ss(&values[153]);
c44_1 = _mm_add_ss(c44_1, _mm_mul_ss(a44_1, b44));
_mm_store_ss(&C[(i*56)+44], c44_1);
#else
C[(i*56)+23] += values[152] * B[(i*56)+44];
C[(i*56)+44] += values[153] * B[(i*56)+44];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b45 = _mm_broadcast_ss(&B[(i*56)+45]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b45 = _mm_load_ss(&B[(i*56)+45]);
b45 = _mm_shuffle_ps(b45, b45, 0x00);
#endif
__m128 c45_0 = _mm_load_ss(&C[(i*56)+24]);
__m128 a45_0 = _mm_load_ss(&values[154]);
c45_0 = _mm_add_ss(c45_0, _mm_mul_ss(a45_0, b45));
_mm_store_ss(&C[(i*56)+24], c45_0);
__m128 c45_1 = _mm_load_ss(&C[(i*56)+45]);
__m128 a45_1 = _mm_load_ss(&values[155]);
c45_1 = _mm_add_ss(c45_1, _mm_mul_ss(a45_1, b45));
_mm_store_ss(&C[(i*56)+45], c45_1);
#else
C[(i*56)+24] += values[154] * B[(i*56)+45];
C[(i*56)+45] += values[155] * B[(i*56)+45];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b46 = _mm_broadcast_ss(&B[(i*56)+46]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b46 = _mm_load_ss(&B[(i*56)+46]);
b46 = _mm_shuffle_ps(b46, b46, 0x00);
#endif
__m128 c46_0 = _mm_load_ss(&C[(i*56)+10]);
__m128 a46_0 = _mm_load_ss(&values[156]);
c46_0 = _mm_add_ss(c46_0, _mm_mul_ss(a46_0, b46));
_mm_store_ss(&C[(i*56)+10], c46_0);
__m128 c46_1 = _mm_load_ss(&C[(i*56)+25]);
__m128 a46_1 = _mm_load_ss(&values[157]);
c46_1 = _mm_add_ss(c46_1, _mm_mul_ss(a46_1, b46));
_mm_store_ss(&C[(i*56)+25], c46_1);
__m128 c46_2 = _mm_load_ss(&C[(i*56)+46]);
__m128 a46_2 = _mm_load_ss(&values[158]);
c46_2 = _mm_add_ss(c46_2, _mm_mul_ss(a46_2, b46));
_mm_store_ss(&C[(i*56)+46], c46_2);
#else
C[(i*56)+10] += values[156] * B[(i*56)+46];
C[(i*56)+25] += values[157] * B[(i*56)+46];
C[(i*56)+46] += values[158] * B[(i*56)+46];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b47 = _mm_broadcast_ss(&B[(i*56)+47]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b47 = _mm_load_ss(&B[(i*56)+47]);
b47 = _mm_shuffle_ps(b47, b47, 0x00);
#endif
__m128 c47_0 = _mm_load_ss(&C[(i*56)+11]);
__m128 a47_0 = _mm_load_ss(&values[159]);
c47_0 = _mm_add_ss(c47_0, _mm_mul_ss(a47_0, b47));
_mm_store_ss(&C[(i*56)+11], c47_0);
__m128 c47_1 = _mm_load_ss(&C[(i*56)+26]);
__m128 a47_1 = _mm_load_ss(&values[160]);
c47_1 = _mm_add_ss(c47_1, _mm_mul_ss(a47_1, b47));
_mm_store_ss(&C[(i*56)+26], c47_1);
__m128 c47_2 = _mm_load_ss(&C[(i*56)+47]);
__m128 a47_2 = _mm_load_ss(&values[161]);
c47_2 = _mm_add_ss(c47_2, _mm_mul_ss(a47_2, b47));
_mm_store_ss(&C[(i*56)+47], c47_2);
#else
C[(i*56)+11] += values[159] * B[(i*56)+47];
C[(i*56)+26] += values[160] * B[(i*56)+47];
C[(i*56)+47] += values[161] * B[(i*56)+47];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b48 = _mm_broadcast_ss(&B[(i*56)+48]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b48 = _mm_load_ss(&B[(i*56)+48]);
b48 = _mm_shuffle_ps(b48, b48, 0x00);
#endif
__m128 c48_0 = _mm_load_ss(&C[(i*56)+12]);
__m128 a48_0 = _mm_load_ss(&values[162]);
c48_0 = _mm_add_ss(c48_0, _mm_mul_ss(a48_0, b48));
_mm_store_ss(&C[(i*56)+12], c48_0);
__m128 c48_1 = _mm_load_ss(&C[(i*56)+27]);
__m128 a48_1 = _mm_load_ss(&values[163]);
c48_1 = _mm_add_ss(c48_1, _mm_mul_ss(a48_1, b48));
_mm_store_ss(&C[(i*56)+27], c48_1);
__m128 c48_2 = _mm_load_ss(&C[(i*56)+48]);
__m128 a48_2 = _mm_load_ss(&values[164]);
c48_2 = _mm_add_ss(c48_2, _mm_mul_ss(a48_2, b48));
_mm_store_ss(&C[(i*56)+48], c48_2);
#else
C[(i*56)+12] += values[162] * B[(i*56)+48];
C[(i*56)+27] += values[163] * B[(i*56)+48];
C[(i*56)+48] += values[164] * B[(i*56)+48];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b49 = _mm_broadcast_ss(&B[(i*56)+49]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b49 = _mm_load_ss(&B[(i*56)+49]);
b49 = _mm_shuffle_ps(b49, b49, 0x00);
#endif
__m128 c49_0 = _mm_load_ss(&C[(i*56)+13]);
__m128 a49_0 = _mm_load_ss(&values[165]);
c49_0 = _mm_add_ss(c49_0, _mm_mul_ss(a49_0, b49));
_mm_store_ss(&C[(i*56)+13], c49_0);
__m128 c49_1 = _mm_load_ss(&C[(i*56)+28]);
__m128 a49_1 = _mm_load_ss(&values[166]);
c49_1 = _mm_add_ss(c49_1, _mm_mul_ss(a49_1, b49));
_mm_store_ss(&C[(i*56)+28], c49_1);
__m128 c49_2 = _mm_load_ss(&C[(i*56)+49]);
__m128 a49_2 = _mm_load_ss(&values[167]);
c49_2 = _mm_add_ss(c49_2, _mm_mul_ss(a49_2, b49));
_mm_store_ss(&C[(i*56)+49], c49_2);
#else
C[(i*56)+13] += values[165] * B[(i*56)+49];
C[(i*56)+28] += values[166] * B[(i*56)+49];
C[(i*56)+49] += values[167] * B[(i*56)+49];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b50 = _mm_broadcast_ss(&B[(i*56)+50]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b50 = _mm_load_ss(&B[(i*56)+50]);
b50 = _mm_shuffle_ps(b50, b50, 0x00);
#endif
__m128 c50_0 = _mm_load_ss(&C[(i*56)+4]);
__m128 a50_0 = _mm_load_ss(&values[168]);
c50_0 = _mm_add_ss(c50_0, _mm_mul_ss(a50_0, b50));
_mm_store_ss(&C[(i*56)+4], c50_0);
__m128 c50_1 = _mm_load_ss(&C[(i*56)+14]);
__m128 a50_1 = _mm_load_ss(&values[169]);
c50_1 = _mm_add_ss(c50_1, _mm_mul_ss(a50_1, b50));
_mm_store_ss(&C[(i*56)+14], c50_1);
__m128 c50_2 = _mm_load_ss(&C[(i*56)+29]);
__m128 a50_2 = _mm_load_ss(&values[170]);
c50_2 = _mm_add_ss(c50_2, _mm_mul_ss(a50_2, b50));
_mm_store_ss(&C[(i*56)+29], c50_2);
__m128 c50_3 = _mm_load_ss(&C[(i*56)+50]);
__m128 a50_3 = _mm_load_ss(&values[171]);
c50_3 = _mm_add_ss(c50_3, _mm_mul_ss(a50_3, b50));
_mm_store_ss(&C[(i*56)+50], c50_3);
#else
C[(i*56)+4] += values[168] * B[(i*56)+50];
C[(i*56)+14] += values[169] * B[(i*56)+50];
C[(i*56)+29] += values[170] * B[(i*56)+50];
C[(i*56)+50] += values[171] * B[(i*56)+50];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b51 = _mm_broadcast_ss(&B[(i*56)+51]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b51 = _mm_load_ss(&B[(i*56)+51]);
b51 = _mm_shuffle_ps(b51, b51, 0x00);
#endif
__m128 c51_0 = _mm_load_ss(&C[(i*56)+5]);
__m128 a51_0 = _mm_load_ss(&values[172]);
c51_0 = _mm_add_ss(c51_0, _mm_mul_ss(a51_0, b51));
_mm_store_ss(&C[(i*56)+5], c51_0);
__m128 c51_1 = _mm_load_ss(&C[(i*56)+15]);
__m128 a51_1 = _mm_load_ss(&values[173]);
c51_1 = _mm_add_ss(c51_1, _mm_mul_ss(a51_1, b51));
_mm_store_ss(&C[(i*56)+15], c51_1);
__m128 c51_2 = _mm_load_ss(&C[(i*56)+30]);
__m128 a51_2 = _mm_load_ss(&values[174]);
c51_2 = _mm_add_ss(c51_2, _mm_mul_ss(a51_2, b51));
_mm_store_ss(&C[(i*56)+30], c51_2);
__m128 c51_3 = _mm_load_ss(&C[(i*56)+51]);
__m128 a51_3 = _mm_load_ss(&values[175]);
c51_3 = _mm_add_ss(c51_3, _mm_mul_ss(a51_3, b51));
_mm_store_ss(&C[(i*56)+51], c51_3);
#else
C[(i*56)+5] += values[172] * B[(i*56)+51];
C[(i*56)+15] += values[173] * B[(i*56)+51];
C[(i*56)+30] += values[174] * B[(i*56)+51];
C[(i*56)+51] += values[175] * B[(i*56)+51];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b52 = _mm_broadcast_ss(&B[(i*56)+52]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b52 = _mm_load_ss(&B[(i*56)+52]);
b52 = _mm_shuffle_ps(b52, b52, 0x00);
#endif
__m128 c52_0 = _mm_load_ss(&C[(i*56)+6]);
__m128 a52_0 = _mm_load_ss(&values[176]);
c52_0 = _mm_add_ss(c52_0, _mm_mul_ss(a52_0, b52));
_mm_store_ss(&C[(i*56)+6], c52_0);
__m128 c52_1 = _mm_load_ss(&C[(i*56)+16]);
__m128 a52_1 = _mm_load_ss(&values[177]);
c52_1 = _mm_add_ss(c52_1, _mm_mul_ss(a52_1, b52));
_mm_store_ss(&C[(i*56)+16], c52_1);
__m128 c52_2 = _mm_load_ss(&C[(i*56)+31]);
__m128 a52_2 = _mm_load_ss(&values[178]);
c52_2 = _mm_add_ss(c52_2, _mm_mul_ss(a52_2, b52));
_mm_store_ss(&C[(i*56)+31], c52_2);
__m128 c52_3 = _mm_load_ss(&C[(i*56)+52]);
__m128 a52_3 = _mm_load_ss(&values[179]);
c52_3 = _mm_add_ss(c52_3, _mm_mul_ss(a52_3, b52));
_mm_store_ss(&C[(i*56)+52], c52_3);
#else
C[(i*56)+6] += values[176] * B[(i*56)+52];
C[(i*56)+16] += values[177] * B[(i*56)+52];
C[(i*56)+31] += values[178] * B[(i*56)+52];
C[(i*56)+52] += values[179] * B[(i*56)+52];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b53 = _mm_broadcast_ss(&B[(i*56)+53]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b53 = _mm_load_ss(&B[(i*56)+53]);
b53 = _mm_shuffle_ps(b53, b53, 0x00);
#endif
__m128 c53_0 = _mm_load_ss(&C[(i*56)+1]);
__m128 a53_0 = _mm_load_ss(&values[180]);
c53_0 = _mm_add_ss(c53_0, _mm_mul_ss(a53_0, b53));
_mm_store_ss(&C[(i*56)+1], c53_0);
__m128 c53_1 = _mm_load_ss(&C[(i*56)+7]);
__m128 a53_1 = _mm_load_ss(&values[181]);
c53_1 = _mm_add_ss(c53_1, _mm_mul_ss(a53_1, b53));
_mm_store_ss(&C[(i*56)+7], c53_1);
__m128 c53_2 = _mm_load_ss(&C[(i*56)+17]);
__m128 a53_2 = _mm_load_ss(&values[182]);
c53_2 = _mm_add_ss(c53_2, _mm_mul_ss(a53_2, b53));
_mm_store_ss(&C[(i*56)+17], c53_2);
__m128 c53_3 = _mm_load_ss(&C[(i*56)+32]);
__m128 a53_3 = _mm_load_ss(&values[183]);
c53_3 = _mm_add_ss(c53_3, _mm_mul_ss(a53_3, b53));
_mm_store_ss(&C[(i*56)+32], c53_3);
__m128 c53_4 = _mm_load_ss(&C[(i*56)+53]);
__m128 a53_4 = _mm_load_ss(&values[184]);
c53_4 = _mm_add_ss(c53_4, _mm_mul_ss(a53_4, b53));
_mm_store_ss(&C[(i*56)+53], c53_4);
#else
C[(i*56)+1] += values[180] * B[(i*56)+53];
C[(i*56)+7] += values[181] * B[(i*56)+53];
C[(i*56)+17] += values[182] * B[(i*56)+53];
C[(i*56)+32] += values[183] * B[(i*56)+53];
C[(i*56)+53] += values[184] * B[(i*56)+53];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b54 = _mm_broadcast_ss(&B[(i*56)+54]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b54 = _mm_load_ss(&B[(i*56)+54]);
b54 = _mm_shuffle_ps(b54, b54, 0x00);
#endif
__m128 c54_0 = _mm_load_ss(&C[(i*56)+2]);
__m128 a54_0 = _mm_load_ss(&values[185]);
c54_0 = _mm_add_ss(c54_0, _mm_mul_ss(a54_0, b54));
_mm_store_ss(&C[(i*56)+2], c54_0);
__m128 c54_1 = _mm_load_ss(&C[(i*56)+8]);
__m128 a54_1 = _mm_load_ss(&values[186]);
c54_1 = _mm_add_ss(c54_1, _mm_mul_ss(a54_1, b54));
_mm_store_ss(&C[(i*56)+8], c54_1);
__m128 c54_2 = _mm_load_ss(&C[(i*56)+18]);
__m128 a54_2 = _mm_load_ss(&values[187]);
c54_2 = _mm_add_ss(c54_2, _mm_mul_ss(a54_2, b54));
_mm_store_ss(&C[(i*56)+18], c54_2);
__m128 c54_3 = _mm_load_ss(&C[(i*56)+33]);
__m128 a54_3 = _mm_load_ss(&values[188]);
c54_3 = _mm_add_ss(c54_3, _mm_mul_ss(a54_3, b54));
_mm_store_ss(&C[(i*56)+33], c54_3);
__m128 c54_4 = _mm_load_ss(&C[(i*56)+54]);
__m128 a54_4 = _mm_load_ss(&values[189]);
c54_4 = _mm_add_ss(c54_4, _mm_mul_ss(a54_4, b54));
_mm_store_ss(&C[(i*56)+54], c54_4);
#else
C[(i*56)+2] += values[185] * B[(i*56)+54];
C[(i*56)+8] += values[186] * B[(i*56)+54];
C[(i*56)+18] += values[187] * B[(i*56)+54];
C[(i*56)+33] += values[188] * B[(i*56)+54];
C[(i*56)+54] += values[189] * B[(i*56)+54];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b55 = _mm_broadcast_ss(&B[(i*56)+55]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b55 = _mm_load_ss(&B[(i*56)+55]);
b55 = _mm_shuffle_ps(b55, b55, 0x00);
#endif
__m128 c55_0 = _mm_load_ss(&C[(i*56)+0]);
__m128 a55_0 = _mm_load_ss(&values[190]);
c55_0 = _mm_add_ss(c55_0, _mm_mul_ss(a55_0, b55));
_mm_store_ss(&C[(i*56)+0], c55_0);
__m128 c55_1 = _mm_load_ss(&C[(i*56)+3]);
__m128 a55_1 = _mm_load_ss(&values[191]);
c55_1 = _mm_add_ss(c55_1, _mm_mul_ss(a55_1, b55));
_mm_store_ss(&C[(i*56)+3], c55_1);
__m128 c55_2 = _mm_load_ss(&C[(i*56)+9]);
__m128 a55_2 = _mm_load_ss(&values[192]);
c55_2 = _mm_add_ss(c55_2, _mm_mul_ss(a55_2, b55));
_mm_store_ss(&C[(i*56)+9], c55_2);
__m128 c55_3 = _mm_load_ss(&C[(i*56)+19]);
__m128 a55_3 = _mm_load_ss(&values[193]);
c55_3 = _mm_add_ss(c55_3, _mm_mul_ss(a55_3, b55));
_mm_store_ss(&C[(i*56)+19], c55_3);
__m128 c55_4 = _mm_load_ss(&C[(i*56)+34]);
__m128 a55_4 = _mm_load_ss(&values[194]);
c55_4 = _mm_add_ss(c55_4, _mm_mul_ss(a55_4, b55));
_mm_store_ss(&C[(i*56)+34], c55_4);
__m128 c55_5 = _mm_load_ss(&C[(i*56)+55]);
__m128 a55_5 = _mm_load_ss(&values[195]);
c55_5 = _mm_add_ss(c55_5, _mm_mul_ss(a55_5, b55));
_mm_store_ss(&C[(i*56)+55], c55_5);
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

void ssparse_starMatrix_m84_n9_k9_ldA88_ldBna7_ldC88_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
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

void ssparse_fM1DivM_m84_n9_k84_ldAna7_ldB88_ldC88_beta0_pfsigonly(const float* values, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
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
__m128 b0 = _mm_broadcast_ss(&B[(i*88)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b0 = _mm_load_ss(&B[(i*88)+0]);
b0 = _mm_shuffle_ps(b0, b0, 0x00);
#endif
__m128 c0_0 = _mm_load_ss(&C[(i*88)+0]);
__m128 a0_0 = _mm_load_ss(&values[0]);
c0_0 = _mm_add_ss(c0_0, _mm_mul_ss(a0_0, b0));
_mm_store_ss(&C[(i*88)+0], c0_0);
__m128 c0_1 = _mm_load_ss(&C[(i*88)+3]);
__m128 a0_1 = _mm_load_ss(&values[1]);
c0_1 = _mm_add_ss(c0_1, _mm_mul_ss(a0_1, b0));
_mm_store_ss(&C[(i*88)+3], c0_1);
__m128 c0_2 = _mm_load_ss(&C[(i*88)+9]);
__m128 a0_2 = _mm_load_ss(&values[2]);
c0_2 = _mm_add_ss(c0_2, _mm_mul_ss(a0_2, b0));
_mm_store_ss(&C[(i*88)+9], c0_2);
__m128 c0_3 = _mm_load_ss(&C[(i*88)+19]);
__m128 a0_3 = _mm_load_ss(&values[3]);
c0_3 = _mm_add_ss(c0_3, _mm_mul_ss(a0_3, b0));
_mm_store_ss(&C[(i*88)+19], c0_3);
__m128 c0_4 = _mm_load_ss(&C[(i*88)+34]);
__m128 a0_4 = _mm_load_ss(&values[4]);
c0_4 = _mm_add_ss(c0_4, _mm_mul_ss(a0_4, b0));
_mm_store_ss(&C[(i*88)+34], c0_4);
__m128 c0_5 = _mm_load_ss(&C[(i*88)+55]);
__m128 a0_5 = _mm_load_ss(&values[5]);
c0_5 = _mm_add_ss(c0_5, _mm_mul_ss(a0_5, b0));
_mm_store_ss(&C[(i*88)+55], c0_5);
__m128 c0_6 = _mm_load_ss(&C[(i*88)+83]);
__m128 a0_6 = _mm_load_ss(&values[6]);
c0_6 = _mm_add_ss(c0_6, _mm_mul_ss(a0_6, b0));
_mm_store_ss(&C[(i*88)+83], c0_6);
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
__m128 b1 = _mm_broadcast_ss(&B[(i*88)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b1 = _mm_load_ss(&B[(i*88)+1]);
b1 = _mm_shuffle_ps(b1, b1, 0x00);
#endif
__m128 c1_0 = _mm_load_ss(&C[(i*88)+1]);
__m128 a1_0 = _mm_load_ss(&values[7]);
c1_0 = _mm_add_ss(c1_0, _mm_mul_ss(a1_0, b1));
_mm_store_ss(&C[(i*88)+1], c1_0);
__m128 c1_1 = _mm_load_ss(&C[(i*88)+7]);
__m128 a1_1 = _mm_load_ss(&values[8]);
c1_1 = _mm_add_ss(c1_1, _mm_mul_ss(a1_1, b1));
_mm_store_ss(&C[(i*88)+7], c1_1);
__m128 c1_2 = _mm_load_ss(&C[(i*88)+17]);
__m128 a1_2 = _mm_load_ss(&values[9]);
c1_2 = _mm_add_ss(c1_2, _mm_mul_ss(a1_2, b1));
_mm_store_ss(&C[(i*88)+17], c1_2);
__m128 c1_3 = _mm_load_ss(&C[(i*88)+32]);
__m128 a1_3 = _mm_load_ss(&values[10]);
c1_3 = _mm_add_ss(c1_3, _mm_mul_ss(a1_3, b1));
_mm_store_ss(&C[(i*88)+32], c1_3);
__m128 c1_4 = _mm_load_ss(&C[(i*88)+53]);
__m128 a1_4 = _mm_load_ss(&values[11]);
c1_4 = _mm_add_ss(c1_4, _mm_mul_ss(a1_4, b1));
_mm_store_ss(&C[(i*88)+53], c1_4);
__m128 c1_5 = _mm_load_ss(&C[(i*88)+81]);
__m128 a1_5 = _mm_load_ss(&values[12]);
c1_5 = _mm_add_ss(c1_5, _mm_mul_ss(a1_5, b1));
_mm_store_ss(&C[(i*88)+81], c1_5);
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
__m128 b2 = _mm_broadcast_ss(&B[(i*88)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b2 = _mm_load_ss(&B[(i*88)+2]);
b2 = _mm_shuffle_ps(b2, b2, 0x00);
#endif
__m128 c2_0 = _mm_load_ss(&C[(i*88)+2]);
__m128 a2_0 = _mm_load_ss(&values[13]);
c2_0 = _mm_add_ss(c2_0, _mm_mul_ss(a2_0, b2));
_mm_store_ss(&C[(i*88)+2], c2_0);
__m128 c2_1 = _mm_load_ss(&C[(i*88)+8]);
__m128 a2_1 = _mm_load_ss(&values[14]);
c2_1 = _mm_add_ss(c2_1, _mm_mul_ss(a2_1, b2));
_mm_store_ss(&C[(i*88)+8], c2_1);
__m128 c2_2 = _mm_load_ss(&C[(i*88)+18]);
__m128 a2_2 = _mm_load_ss(&values[15]);
c2_2 = _mm_add_ss(c2_2, _mm_mul_ss(a2_2, b2));
_mm_store_ss(&C[(i*88)+18], c2_2);
__m128 c2_3 = _mm_load_ss(&C[(i*88)+33]);
__m128 a2_3 = _mm_load_ss(&values[16]);
c2_3 = _mm_add_ss(c2_3, _mm_mul_ss(a2_3, b2));
_mm_store_ss(&C[(i*88)+33], c2_3);
__m128 c2_4 = _mm_load_ss(&C[(i*88)+54]);
__m128 a2_4 = _mm_load_ss(&values[17]);
c2_4 = _mm_add_ss(c2_4, _mm_mul_ss(a2_4, b2));
_mm_store_ss(&C[(i*88)+54], c2_4);
__m128 c2_5 = _mm_load_ss(&C[(i*88)+82]);
__m128 a2_5 = _mm_load_ss(&values[18]);
c2_5 = _mm_add_ss(c2_5, _mm_mul_ss(a2_5, b2));
_mm_store_ss(&C[(i*88)+82], c2_5);
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
__m128 b3 = _mm_broadcast_ss(&B[(i*88)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b3 = _mm_load_ss(&B[(i*88)+3]);
b3 = _mm_shuffle_ps(b3, b3, 0x00);
#endif
__m128 c3_0 = _mm_load_ss(&C[(i*88)+0]);
__m128 a3_0 = _mm_load_ss(&values[19]);
c3_0 = _mm_add_ss(c3_0, _mm_mul_ss(a3_0, b3));
_mm_store_ss(&C[(i*88)+0], c3_0);
__m128 c3_1 = _mm_load_ss(&C[(i*88)+3]);
__m128 a3_1 = _mm_load_ss(&values[20]);
c3_1 = _mm_add_ss(c3_1, _mm_mul_ss(a3_1, b3));
_mm_store_ss(&C[(i*88)+3], c3_1);
__m128 c3_2 = _mm_load_ss(&C[(i*88)+9]);
__m128 a3_2 = _mm_load_ss(&values[21]);
c3_2 = _mm_add_ss(c3_2, _mm_mul_ss(a3_2, b3));
_mm_store_ss(&C[(i*88)+9], c3_2);
__m128 c3_3 = _mm_load_ss(&C[(i*88)+19]);
__m128 a3_3 = _mm_load_ss(&values[22]);
c3_3 = _mm_add_ss(c3_3, _mm_mul_ss(a3_3, b3));
_mm_store_ss(&C[(i*88)+19], c3_3);
__m128 c3_4 = _mm_load_ss(&C[(i*88)+34]);
__m128 a3_4 = _mm_load_ss(&values[23]);
c3_4 = _mm_add_ss(c3_4, _mm_mul_ss(a3_4, b3));
_mm_store_ss(&C[(i*88)+34], c3_4);
__m128 c3_5 = _mm_load_ss(&C[(i*88)+55]);
__m128 a3_5 = _mm_load_ss(&values[24]);
c3_5 = _mm_add_ss(c3_5, _mm_mul_ss(a3_5, b3));
_mm_store_ss(&C[(i*88)+55], c3_5);
__m128 c3_6 = _mm_load_ss(&C[(i*88)+83]);
__m128 a3_6 = _mm_load_ss(&values[25]);
c3_6 = _mm_add_ss(c3_6, _mm_mul_ss(a3_6, b3));
_mm_store_ss(&C[(i*88)+83], c3_6);
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
__m128 b4 = _mm_broadcast_ss(&B[(i*88)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b4 = _mm_load_ss(&B[(i*88)+4]);
b4 = _mm_shuffle_ps(b4, b4, 0x00);
#endif
__m128 c4_0 = _mm_load_ss(&C[(i*88)+4]);
__m128 a4_0 = _mm_load_ss(&values[26]);
c4_0 = _mm_add_ss(c4_0, _mm_mul_ss(a4_0, b4));
_mm_store_ss(&C[(i*88)+4], c4_0);
__m128 c4_1 = _mm_load_ss(&C[(i*88)+14]);
__m128 a4_1 = _mm_load_ss(&values[27]);
c4_1 = _mm_add_ss(c4_1, _mm_mul_ss(a4_1, b4));
_mm_store_ss(&C[(i*88)+14], c4_1);
__m128 c4_2 = _mm_load_ss(&C[(i*88)+29]);
__m128 a4_2 = _mm_load_ss(&values[28]);
c4_2 = _mm_add_ss(c4_2, _mm_mul_ss(a4_2, b4));
_mm_store_ss(&C[(i*88)+29], c4_2);
__m128 c4_3 = _mm_load_ss(&C[(i*88)+50]);
__m128 a4_3 = _mm_load_ss(&values[29]);
c4_3 = _mm_add_ss(c4_3, _mm_mul_ss(a4_3, b4));
_mm_store_ss(&C[(i*88)+50], c4_3);
__m128 c4_4 = _mm_load_ss(&C[(i*88)+78]);
__m128 a4_4 = _mm_load_ss(&values[30]);
c4_4 = _mm_add_ss(c4_4, _mm_mul_ss(a4_4, b4));
_mm_store_ss(&C[(i*88)+78], c4_4);
#else
C[(i*88)+4] += values[26] * B[(i*88)+4];
C[(i*88)+14] += values[27] * B[(i*88)+4];
C[(i*88)+29] += values[28] * B[(i*88)+4];
C[(i*88)+50] += values[29] * B[(i*88)+4];
C[(i*88)+78] += values[30] * B[(i*88)+4];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b5 = _mm_broadcast_ss(&B[(i*88)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b5 = _mm_load_ss(&B[(i*88)+5]);
b5 = _mm_shuffle_ps(b5, b5, 0x00);
#endif
__m128 c5_0 = _mm_load_ss(&C[(i*88)+5]);
__m128 a5_0 = _mm_load_ss(&values[31]);
c5_0 = _mm_add_ss(c5_0, _mm_mul_ss(a5_0, b5));
_mm_store_ss(&C[(i*88)+5], c5_0);
__m128 c5_1 = _mm_load_ss(&C[(i*88)+15]);
__m128 a5_1 = _mm_load_ss(&values[32]);
c5_1 = _mm_add_ss(c5_1, _mm_mul_ss(a5_1, b5));
_mm_store_ss(&C[(i*88)+15], c5_1);
__m128 c5_2 = _mm_load_ss(&C[(i*88)+30]);
__m128 a5_2 = _mm_load_ss(&values[33]);
c5_2 = _mm_add_ss(c5_2, _mm_mul_ss(a5_2, b5));
_mm_store_ss(&C[(i*88)+30], c5_2);
__m128 c5_3 = _mm_load_ss(&C[(i*88)+51]);
__m128 a5_3 = _mm_load_ss(&values[34]);
c5_3 = _mm_add_ss(c5_3, _mm_mul_ss(a5_3, b5));
_mm_store_ss(&C[(i*88)+51], c5_3);
__m128 c5_4 = _mm_load_ss(&C[(i*88)+79]);
__m128 a5_4 = _mm_load_ss(&values[35]);
c5_4 = _mm_add_ss(c5_4, _mm_mul_ss(a5_4, b5));
_mm_store_ss(&C[(i*88)+79], c5_4);
#else
C[(i*88)+5] += values[31] * B[(i*88)+5];
C[(i*88)+15] += values[32] * B[(i*88)+5];
C[(i*88)+30] += values[33] * B[(i*88)+5];
C[(i*88)+51] += values[34] * B[(i*88)+5];
C[(i*88)+79] += values[35] * B[(i*88)+5];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b6 = _mm_broadcast_ss(&B[(i*88)+6]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b6 = _mm_load_ss(&B[(i*88)+6]);
b6 = _mm_shuffle_ps(b6, b6, 0x00);
#endif
__m128 c6_0 = _mm_load_ss(&C[(i*88)+6]);
__m128 a6_0 = _mm_load_ss(&values[36]);
c6_0 = _mm_add_ss(c6_0, _mm_mul_ss(a6_0, b6));
_mm_store_ss(&C[(i*88)+6], c6_0);
__m128 c6_1 = _mm_load_ss(&C[(i*88)+16]);
__m128 a6_1 = _mm_load_ss(&values[37]);
c6_1 = _mm_add_ss(c6_1, _mm_mul_ss(a6_1, b6));
_mm_store_ss(&C[(i*88)+16], c6_1);
__m128 c6_2 = _mm_load_ss(&C[(i*88)+31]);
__m128 a6_2 = _mm_load_ss(&values[38]);
c6_2 = _mm_add_ss(c6_2, _mm_mul_ss(a6_2, b6));
_mm_store_ss(&C[(i*88)+31], c6_2);
__m128 c6_3 = _mm_load_ss(&C[(i*88)+52]);
__m128 a6_3 = _mm_load_ss(&values[39]);
c6_3 = _mm_add_ss(c6_3, _mm_mul_ss(a6_3, b6));
_mm_store_ss(&C[(i*88)+52], c6_3);
__m128 c6_4 = _mm_load_ss(&C[(i*88)+80]);
__m128 a6_4 = _mm_load_ss(&values[40]);
c6_4 = _mm_add_ss(c6_4, _mm_mul_ss(a6_4, b6));
_mm_store_ss(&C[(i*88)+80], c6_4);
#else
C[(i*88)+6] += values[36] * B[(i*88)+6];
C[(i*88)+16] += values[37] * B[(i*88)+6];
C[(i*88)+31] += values[38] * B[(i*88)+6];
C[(i*88)+52] += values[39] * B[(i*88)+6];
C[(i*88)+80] += values[40] * B[(i*88)+6];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b7 = _mm_broadcast_ss(&B[(i*88)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b7 = _mm_load_ss(&B[(i*88)+7]);
b7 = _mm_shuffle_ps(b7, b7, 0x00);
#endif
__m128 c7_0 = _mm_load_ss(&C[(i*88)+1]);
__m128 a7_0 = _mm_load_ss(&values[41]);
c7_0 = _mm_add_ss(c7_0, _mm_mul_ss(a7_0, b7));
_mm_store_ss(&C[(i*88)+1], c7_0);
__m128 c7_1 = _mm_load_ss(&C[(i*88)+7]);
__m128 a7_1 = _mm_load_ss(&values[42]);
c7_1 = _mm_add_ss(c7_1, _mm_mul_ss(a7_1, b7));
_mm_store_ss(&C[(i*88)+7], c7_1);
__m128 c7_2 = _mm_load_ss(&C[(i*88)+17]);
__m128 a7_2 = _mm_load_ss(&values[43]);
c7_2 = _mm_add_ss(c7_2, _mm_mul_ss(a7_2, b7));
_mm_store_ss(&C[(i*88)+17], c7_2);
__m128 c7_3 = _mm_load_ss(&C[(i*88)+32]);
__m128 a7_3 = _mm_load_ss(&values[44]);
c7_3 = _mm_add_ss(c7_3, _mm_mul_ss(a7_3, b7));
_mm_store_ss(&C[(i*88)+32], c7_3);
__m128 c7_4 = _mm_load_ss(&C[(i*88)+53]);
__m128 a7_4 = _mm_load_ss(&values[45]);
c7_4 = _mm_add_ss(c7_4, _mm_mul_ss(a7_4, b7));
_mm_store_ss(&C[(i*88)+53], c7_4);
__m128 c7_5 = _mm_load_ss(&C[(i*88)+81]);
__m128 a7_5 = _mm_load_ss(&values[46]);
c7_5 = _mm_add_ss(c7_5, _mm_mul_ss(a7_5, b7));
_mm_store_ss(&C[(i*88)+81], c7_5);
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
__m128 b8 = _mm_broadcast_ss(&B[(i*88)+8]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b8 = _mm_load_ss(&B[(i*88)+8]);
b8 = _mm_shuffle_ps(b8, b8, 0x00);
#endif
__m128 c8_0 = _mm_load_ss(&C[(i*88)+2]);
__m128 a8_0 = _mm_load_ss(&values[47]);
c8_0 = _mm_add_ss(c8_0, _mm_mul_ss(a8_0, b8));
_mm_store_ss(&C[(i*88)+2], c8_0);
__m128 c8_1 = _mm_load_ss(&C[(i*88)+8]);
__m128 a8_1 = _mm_load_ss(&values[48]);
c8_1 = _mm_add_ss(c8_1, _mm_mul_ss(a8_1, b8));
_mm_store_ss(&C[(i*88)+8], c8_1);
__m128 c8_2 = _mm_load_ss(&C[(i*88)+18]);
__m128 a8_2 = _mm_load_ss(&values[49]);
c8_2 = _mm_add_ss(c8_2, _mm_mul_ss(a8_2, b8));
_mm_store_ss(&C[(i*88)+18], c8_2);
__m128 c8_3 = _mm_load_ss(&C[(i*88)+33]);
__m128 a8_3 = _mm_load_ss(&values[50]);
c8_3 = _mm_add_ss(c8_3, _mm_mul_ss(a8_3, b8));
_mm_store_ss(&C[(i*88)+33], c8_3);
__m128 c8_4 = _mm_load_ss(&C[(i*88)+54]);
__m128 a8_4 = _mm_load_ss(&values[51]);
c8_4 = _mm_add_ss(c8_4, _mm_mul_ss(a8_4, b8));
_mm_store_ss(&C[(i*88)+54], c8_4);
__m128 c8_5 = _mm_load_ss(&C[(i*88)+82]);
__m128 a8_5 = _mm_load_ss(&values[52]);
c8_5 = _mm_add_ss(c8_5, _mm_mul_ss(a8_5, b8));
_mm_store_ss(&C[(i*88)+82], c8_5);
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
__m128 b9 = _mm_broadcast_ss(&B[(i*88)+9]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b9 = _mm_load_ss(&B[(i*88)+9]);
b9 = _mm_shuffle_ps(b9, b9, 0x00);
#endif
__m128 c9_0 = _mm_load_ss(&C[(i*88)+0]);
__m128 a9_0 = _mm_load_ss(&values[53]);
c9_0 = _mm_add_ss(c9_0, _mm_mul_ss(a9_0, b9));
_mm_store_ss(&C[(i*88)+0], c9_0);
__m128 c9_1 = _mm_load_ss(&C[(i*88)+3]);
__m128 a9_1 = _mm_load_ss(&values[54]);
c9_1 = _mm_add_ss(c9_1, _mm_mul_ss(a9_1, b9));
_mm_store_ss(&C[(i*88)+3], c9_1);
__m128 c9_2 = _mm_load_ss(&C[(i*88)+9]);
__m128 a9_2 = _mm_load_ss(&values[55]);
c9_2 = _mm_add_ss(c9_2, _mm_mul_ss(a9_2, b9));
_mm_store_ss(&C[(i*88)+9], c9_2);
__m128 c9_3 = _mm_load_ss(&C[(i*88)+19]);
__m128 a9_3 = _mm_load_ss(&values[56]);
c9_3 = _mm_add_ss(c9_3, _mm_mul_ss(a9_3, b9));
_mm_store_ss(&C[(i*88)+19], c9_3);
__m128 c9_4 = _mm_load_ss(&C[(i*88)+34]);
__m128 a9_4 = _mm_load_ss(&values[57]);
c9_4 = _mm_add_ss(c9_4, _mm_mul_ss(a9_4, b9));
_mm_store_ss(&C[(i*88)+34], c9_4);
__m128 c9_5 = _mm_load_ss(&C[(i*88)+55]);
__m128 a9_5 = _mm_load_ss(&values[58]);
c9_5 = _mm_add_ss(c9_5, _mm_mul_ss(a9_5, b9));
_mm_store_ss(&C[(i*88)+55], c9_5);
__m128 c9_6 = _mm_load_ss(&C[(i*88)+83]);
__m128 a9_6 = _mm_load_ss(&values[59]);
c9_6 = _mm_add_ss(c9_6, _mm_mul_ss(a9_6, b9));
_mm_store_ss(&C[(i*88)+83], c9_6);
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
__m128 b10 = _mm_broadcast_ss(&B[(i*88)+10]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b10 = _mm_load_ss(&B[(i*88)+10]);
b10 = _mm_shuffle_ps(b10, b10, 0x00);
#endif
__m128 c10_0 = _mm_load_ss(&C[(i*88)+10]);
__m128 a10_0 = _mm_load_ss(&values[60]);
c10_0 = _mm_add_ss(c10_0, _mm_mul_ss(a10_0, b10));
_mm_store_ss(&C[(i*88)+10], c10_0);
__m128 c10_1 = _mm_load_ss(&C[(i*88)+25]);
__m128 a10_1 = _mm_load_ss(&values[61]);
c10_1 = _mm_add_ss(c10_1, _mm_mul_ss(a10_1, b10));
_mm_store_ss(&C[(i*88)+25], c10_1);
__m128 c10_2 = _mm_load_ss(&C[(i*88)+46]);
__m128 a10_2 = _mm_load_ss(&values[62]);
c10_2 = _mm_add_ss(c10_2, _mm_mul_ss(a10_2, b10));
_mm_store_ss(&C[(i*88)+46], c10_2);
__m128 c10_3 = _mm_load_ss(&C[(i*88)+74]);
__m128 a10_3 = _mm_load_ss(&values[63]);
c10_3 = _mm_add_ss(c10_3, _mm_mul_ss(a10_3, b10));
_mm_store_ss(&C[(i*88)+74], c10_3);
#else
C[(i*88)+10] += values[60] * B[(i*88)+10];
C[(i*88)+25] += values[61] * B[(i*88)+10];
C[(i*88)+46] += values[62] * B[(i*88)+10];
C[(i*88)+74] += values[63] * B[(i*88)+10];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b11 = _mm_broadcast_ss(&B[(i*88)+11]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b11 = _mm_load_ss(&B[(i*88)+11]);
b11 = _mm_shuffle_ps(b11, b11, 0x00);
#endif
__m128 c11_0 = _mm_load_ss(&C[(i*88)+11]);
__m128 a11_0 = _mm_load_ss(&values[64]);
c11_0 = _mm_add_ss(c11_0, _mm_mul_ss(a11_0, b11));
_mm_store_ss(&C[(i*88)+11], c11_0);
__m128 c11_1 = _mm_load_ss(&C[(i*88)+26]);
__m128 a11_1 = _mm_load_ss(&values[65]);
c11_1 = _mm_add_ss(c11_1, _mm_mul_ss(a11_1, b11));
_mm_store_ss(&C[(i*88)+26], c11_1);
__m128 c11_2 = _mm_load_ss(&C[(i*88)+47]);
__m128 a11_2 = _mm_load_ss(&values[66]);
c11_2 = _mm_add_ss(c11_2, _mm_mul_ss(a11_2, b11));
_mm_store_ss(&C[(i*88)+47], c11_2);
__m128 c11_3 = _mm_load_ss(&C[(i*88)+75]);
__m128 a11_3 = _mm_load_ss(&values[67]);
c11_3 = _mm_add_ss(c11_3, _mm_mul_ss(a11_3, b11));
_mm_store_ss(&C[(i*88)+75], c11_3);
#else
C[(i*88)+11] += values[64] * B[(i*88)+11];
C[(i*88)+26] += values[65] * B[(i*88)+11];
C[(i*88)+47] += values[66] * B[(i*88)+11];
C[(i*88)+75] += values[67] * B[(i*88)+11];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b12 = _mm_broadcast_ss(&B[(i*88)+12]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b12 = _mm_load_ss(&B[(i*88)+12]);
b12 = _mm_shuffle_ps(b12, b12, 0x00);
#endif
__m128 c12_0 = _mm_load_ss(&C[(i*88)+12]);
__m128 a12_0 = _mm_load_ss(&values[68]);
c12_0 = _mm_add_ss(c12_0, _mm_mul_ss(a12_0, b12));
_mm_store_ss(&C[(i*88)+12], c12_0);
__m128 c12_1 = _mm_load_ss(&C[(i*88)+27]);
__m128 a12_1 = _mm_load_ss(&values[69]);
c12_1 = _mm_add_ss(c12_1, _mm_mul_ss(a12_1, b12));
_mm_store_ss(&C[(i*88)+27], c12_1);
__m128 c12_2 = _mm_load_ss(&C[(i*88)+48]);
__m128 a12_2 = _mm_load_ss(&values[70]);
c12_2 = _mm_add_ss(c12_2, _mm_mul_ss(a12_2, b12));
_mm_store_ss(&C[(i*88)+48], c12_2);
__m128 c12_3 = _mm_load_ss(&C[(i*88)+76]);
__m128 a12_3 = _mm_load_ss(&values[71]);
c12_3 = _mm_add_ss(c12_3, _mm_mul_ss(a12_3, b12));
_mm_store_ss(&C[(i*88)+76], c12_3);
#else
C[(i*88)+12] += values[68] * B[(i*88)+12];
C[(i*88)+27] += values[69] * B[(i*88)+12];
C[(i*88)+48] += values[70] * B[(i*88)+12];
C[(i*88)+76] += values[71] * B[(i*88)+12];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b13 = _mm_broadcast_ss(&B[(i*88)+13]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b13 = _mm_load_ss(&B[(i*88)+13]);
b13 = _mm_shuffle_ps(b13, b13, 0x00);
#endif
__m128 c13_0 = _mm_load_ss(&C[(i*88)+13]);
__m128 a13_0 = _mm_load_ss(&values[72]);
c13_0 = _mm_add_ss(c13_0, _mm_mul_ss(a13_0, b13));
_mm_store_ss(&C[(i*88)+13], c13_0);
__m128 c13_1 = _mm_load_ss(&C[(i*88)+28]);
__m128 a13_1 = _mm_load_ss(&values[73]);
c13_1 = _mm_add_ss(c13_1, _mm_mul_ss(a13_1, b13));
_mm_store_ss(&C[(i*88)+28], c13_1);
__m128 c13_2 = _mm_load_ss(&C[(i*88)+49]);
__m128 a13_2 = _mm_load_ss(&values[74]);
c13_2 = _mm_add_ss(c13_2, _mm_mul_ss(a13_2, b13));
_mm_store_ss(&C[(i*88)+49], c13_2);
__m128 c13_3 = _mm_load_ss(&C[(i*88)+77]);
__m128 a13_3 = _mm_load_ss(&values[75]);
c13_3 = _mm_add_ss(c13_3, _mm_mul_ss(a13_3, b13));
_mm_store_ss(&C[(i*88)+77], c13_3);
#else
C[(i*88)+13] += values[72] * B[(i*88)+13];
C[(i*88)+28] += values[73] * B[(i*88)+13];
C[(i*88)+49] += values[74] * B[(i*88)+13];
C[(i*88)+77] += values[75] * B[(i*88)+13];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b14 = _mm_broadcast_ss(&B[(i*88)+14]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b14 = _mm_load_ss(&B[(i*88)+14]);
b14 = _mm_shuffle_ps(b14, b14, 0x00);
#endif
__m128 c14_0 = _mm_load_ss(&C[(i*88)+4]);
__m128 a14_0 = _mm_load_ss(&values[76]);
c14_0 = _mm_add_ss(c14_0, _mm_mul_ss(a14_0, b14));
_mm_store_ss(&C[(i*88)+4], c14_0);
__m128 c14_1 = _mm_load_ss(&C[(i*88)+14]);
__m128 a14_1 = _mm_load_ss(&values[77]);
c14_1 = _mm_add_ss(c14_1, _mm_mul_ss(a14_1, b14));
_mm_store_ss(&C[(i*88)+14], c14_1);
__m128 c14_2 = _mm_load_ss(&C[(i*88)+29]);
__m128 a14_2 = _mm_load_ss(&values[78]);
c14_2 = _mm_add_ss(c14_2, _mm_mul_ss(a14_2, b14));
_mm_store_ss(&C[(i*88)+29], c14_2);
__m128 c14_3 = _mm_load_ss(&C[(i*88)+50]);
__m128 a14_3 = _mm_load_ss(&values[79]);
c14_3 = _mm_add_ss(c14_3, _mm_mul_ss(a14_3, b14));
_mm_store_ss(&C[(i*88)+50], c14_3);
__m128 c14_4 = _mm_load_ss(&C[(i*88)+78]);
__m128 a14_4 = _mm_load_ss(&values[80]);
c14_4 = _mm_add_ss(c14_4, _mm_mul_ss(a14_4, b14));
_mm_store_ss(&C[(i*88)+78], c14_4);
#else
C[(i*88)+4] += values[76] * B[(i*88)+14];
C[(i*88)+14] += values[77] * B[(i*88)+14];
C[(i*88)+29] += values[78] * B[(i*88)+14];
C[(i*88)+50] += values[79] * B[(i*88)+14];
C[(i*88)+78] += values[80] * B[(i*88)+14];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b15 = _mm_broadcast_ss(&B[(i*88)+15]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b15 = _mm_load_ss(&B[(i*88)+15]);
b15 = _mm_shuffle_ps(b15, b15, 0x00);
#endif
__m128 c15_0 = _mm_load_ss(&C[(i*88)+5]);
__m128 a15_0 = _mm_load_ss(&values[81]);
c15_0 = _mm_add_ss(c15_0, _mm_mul_ss(a15_0, b15));
_mm_store_ss(&C[(i*88)+5], c15_0);
__m128 c15_1 = _mm_load_ss(&C[(i*88)+15]);
__m128 a15_1 = _mm_load_ss(&values[82]);
c15_1 = _mm_add_ss(c15_1, _mm_mul_ss(a15_1, b15));
_mm_store_ss(&C[(i*88)+15], c15_1);
__m128 c15_2 = _mm_load_ss(&C[(i*88)+30]);
__m128 a15_2 = _mm_load_ss(&values[83]);
c15_2 = _mm_add_ss(c15_2, _mm_mul_ss(a15_2, b15));
_mm_store_ss(&C[(i*88)+30], c15_2);
__m128 c15_3 = _mm_load_ss(&C[(i*88)+51]);
__m128 a15_3 = _mm_load_ss(&values[84]);
c15_3 = _mm_add_ss(c15_3, _mm_mul_ss(a15_3, b15));
_mm_store_ss(&C[(i*88)+51], c15_3);
__m128 c15_4 = _mm_load_ss(&C[(i*88)+79]);
__m128 a15_4 = _mm_load_ss(&values[85]);
c15_4 = _mm_add_ss(c15_4, _mm_mul_ss(a15_4, b15));
_mm_store_ss(&C[(i*88)+79], c15_4);
#else
C[(i*88)+5] += values[81] * B[(i*88)+15];
C[(i*88)+15] += values[82] * B[(i*88)+15];
C[(i*88)+30] += values[83] * B[(i*88)+15];
C[(i*88)+51] += values[84] * B[(i*88)+15];
C[(i*88)+79] += values[85] * B[(i*88)+15];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b16 = _mm_broadcast_ss(&B[(i*88)+16]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b16 = _mm_load_ss(&B[(i*88)+16]);
b16 = _mm_shuffle_ps(b16, b16, 0x00);
#endif
__m128 c16_0 = _mm_load_ss(&C[(i*88)+6]);
__m128 a16_0 = _mm_load_ss(&values[86]);
c16_0 = _mm_add_ss(c16_0, _mm_mul_ss(a16_0, b16));
_mm_store_ss(&C[(i*88)+6], c16_0);
__m128 c16_1 = _mm_load_ss(&C[(i*88)+16]);
__m128 a16_1 = _mm_load_ss(&values[87]);
c16_1 = _mm_add_ss(c16_1, _mm_mul_ss(a16_1, b16));
_mm_store_ss(&C[(i*88)+16], c16_1);
__m128 c16_2 = _mm_load_ss(&C[(i*88)+31]);
__m128 a16_2 = _mm_load_ss(&values[88]);
c16_2 = _mm_add_ss(c16_2, _mm_mul_ss(a16_2, b16));
_mm_store_ss(&C[(i*88)+31], c16_2);
__m128 c16_3 = _mm_load_ss(&C[(i*88)+52]);
__m128 a16_3 = _mm_load_ss(&values[89]);
c16_3 = _mm_add_ss(c16_3, _mm_mul_ss(a16_3, b16));
_mm_store_ss(&C[(i*88)+52], c16_3);
__m128 c16_4 = _mm_load_ss(&C[(i*88)+80]);
__m128 a16_4 = _mm_load_ss(&values[90]);
c16_4 = _mm_add_ss(c16_4, _mm_mul_ss(a16_4, b16));
_mm_store_ss(&C[(i*88)+80], c16_4);
#else
C[(i*88)+6] += values[86] * B[(i*88)+16];
C[(i*88)+16] += values[87] * B[(i*88)+16];
C[(i*88)+31] += values[88] * B[(i*88)+16];
C[(i*88)+52] += values[89] * B[(i*88)+16];
C[(i*88)+80] += values[90] * B[(i*88)+16];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b17 = _mm_broadcast_ss(&B[(i*88)+17]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b17 = _mm_load_ss(&B[(i*88)+17]);
b17 = _mm_shuffle_ps(b17, b17, 0x00);
#endif
__m128 c17_0 = _mm_load_ss(&C[(i*88)+1]);
__m128 a17_0 = _mm_load_ss(&values[91]);
c17_0 = _mm_add_ss(c17_0, _mm_mul_ss(a17_0, b17));
_mm_store_ss(&C[(i*88)+1], c17_0);
__m128 c17_1 = _mm_load_ss(&C[(i*88)+7]);
__m128 a17_1 = _mm_load_ss(&values[92]);
c17_1 = _mm_add_ss(c17_1, _mm_mul_ss(a17_1, b17));
_mm_store_ss(&C[(i*88)+7], c17_1);
__m128 c17_2 = _mm_load_ss(&C[(i*88)+17]);
__m128 a17_2 = _mm_load_ss(&values[93]);
c17_2 = _mm_add_ss(c17_2, _mm_mul_ss(a17_2, b17));
_mm_store_ss(&C[(i*88)+17], c17_2);
__m128 c17_3 = _mm_load_ss(&C[(i*88)+32]);
__m128 a17_3 = _mm_load_ss(&values[94]);
c17_3 = _mm_add_ss(c17_3, _mm_mul_ss(a17_3, b17));
_mm_store_ss(&C[(i*88)+32], c17_3);
__m128 c17_4 = _mm_load_ss(&C[(i*88)+53]);
__m128 a17_4 = _mm_load_ss(&values[95]);
c17_4 = _mm_add_ss(c17_4, _mm_mul_ss(a17_4, b17));
_mm_store_ss(&C[(i*88)+53], c17_4);
__m128 c17_5 = _mm_load_ss(&C[(i*88)+81]);
__m128 a17_5 = _mm_load_ss(&values[96]);
c17_5 = _mm_add_ss(c17_5, _mm_mul_ss(a17_5, b17));
_mm_store_ss(&C[(i*88)+81], c17_5);
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
__m128 b18 = _mm_broadcast_ss(&B[(i*88)+18]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b18 = _mm_load_ss(&B[(i*88)+18]);
b18 = _mm_shuffle_ps(b18, b18, 0x00);
#endif
__m128 c18_0 = _mm_load_ss(&C[(i*88)+2]);
__m128 a18_0 = _mm_load_ss(&values[97]);
c18_0 = _mm_add_ss(c18_0, _mm_mul_ss(a18_0, b18));
_mm_store_ss(&C[(i*88)+2], c18_0);
__m128 c18_1 = _mm_load_ss(&C[(i*88)+8]);
__m128 a18_1 = _mm_load_ss(&values[98]);
c18_1 = _mm_add_ss(c18_1, _mm_mul_ss(a18_1, b18));
_mm_store_ss(&C[(i*88)+8], c18_1);
__m128 c18_2 = _mm_load_ss(&C[(i*88)+18]);
__m128 a18_2 = _mm_load_ss(&values[99]);
c18_2 = _mm_add_ss(c18_2, _mm_mul_ss(a18_2, b18));
_mm_store_ss(&C[(i*88)+18], c18_2);
__m128 c18_3 = _mm_load_ss(&C[(i*88)+33]);
__m128 a18_3 = _mm_load_ss(&values[100]);
c18_3 = _mm_add_ss(c18_3, _mm_mul_ss(a18_3, b18));
_mm_store_ss(&C[(i*88)+33], c18_3);
__m128 c18_4 = _mm_load_ss(&C[(i*88)+54]);
__m128 a18_4 = _mm_load_ss(&values[101]);
c18_4 = _mm_add_ss(c18_4, _mm_mul_ss(a18_4, b18));
_mm_store_ss(&C[(i*88)+54], c18_4);
__m128 c18_5 = _mm_load_ss(&C[(i*88)+82]);
__m128 a18_5 = _mm_load_ss(&values[102]);
c18_5 = _mm_add_ss(c18_5, _mm_mul_ss(a18_5, b18));
_mm_store_ss(&C[(i*88)+82], c18_5);
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
__m128 b19 = _mm_broadcast_ss(&B[(i*88)+19]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b19 = _mm_load_ss(&B[(i*88)+19]);
b19 = _mm_shuffle_ps(b19, b19, 0x00);
#endif
__m128 c19_0 = _mm_load_ss(&C[(i*88)+0]);
__m128 a19_0 = _mm_load_ss(&values[103]);
c19_0 = _mm_add_ss(c19_0, _mm_mul_ss(a19_0, b19));
_mm_store_ss(&C[(i*88)+0], c19_0);
__m128 c19_1 = _mm_load_ss(&C[(i*88)+3]);
__m128 a19_1 = _mm_load_ss(&values[104]);
c19_1 = _mm_add_ss(c19_1, _mm_mul_ss(a19_1, b19));
_mm_store_ss(&C[(i*88)+3], c19_1);
__m128 c19_2 = _mm_load_ss(&C[(i*88)+9]);
__m128 a19_2 = _mm_load_ss(&values[105]);
c19_2 = _mm_add_ss(c19_2, _mm_mul_ss(a19_2, b19));
_mm_store_ss(&C[(i*88)+9], c19_2);
__m128 c19_3 = _mm_load_ss(&C[(i*88)+19]);
__m128 a19_3 = _mm_load_ss(&values[106]);
c19_3 = _mm_add_ss(c19_3, _mm_mul_ss(a19_3, b19));
_mm_store_ss(&C[(i*88)+19], c19_3);
__m128 c19_4 = _mm_load_ss(&C[(i*88)+34]);
__m128 a19_4 = _mm_load_ss(&values[107]);
c19_4 = _mm_add_ss(c19_4, _mm_mul_ss(a19_4, b19));
_mm_store_ss(&C[(i*88)+34], c19_4);
__m128 c19_5 = _mm_load_ss(&C[(i*88)+55]);
__m128 a19_5 = _mm_load_ss(&values[108]);
c19_5 = _mm_add_ss(c19_5, _mm_mul_ss(a19_5, b19));
_mm_store_ss(&C[(i*88)+55], c19_5);
__m128 c19_6 = _mm_load_ss(&C[(i*88)+83]);
__m128 a19_6 = _mm_load_ss(&values[109]);
c19_6 = _mm_add_ss(c19_6, _mm_mul_ss(a19_6, b19));
_mm_store_ss(&C[(i*88)+83], c19_6);
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
__m128 b20 = _mm_broadcast_ss(&B[(i*88)+20]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b20 = _mm_load_ss(&B[(i*88)+20]);
b20 = _mm_shuffle_ps(b20, b20, 0x00);
#endif
__m128 c20_0 = _mm_load_ss(&C[(i*88)+20]);
__m128 a20_0 = _mm_load_ss(&values[110]);
c20_0 = _mm_add_ss(c20_0, _mm_mul_ss(a20_0, b20));
_mm_store_ss(&C[(i*88)+20], c20_0);
__m128 c20_1 = _mm_load_ss(&C[(i*88)+41]);
__m128 a20_1 = _mm_load_ss(&values[111]);
c20_1 = _mm_add_ss(c20_1, _mm_mul_ss(a20_1, b20));
_mm_store_ss(&C[(i*88)+41], c20_1);
__m128 c20_2 = _mm_load_ss(&C[(i*88)+69]);
__m128 a20_2 = _mm_load_ss(&values[112]);
c20_2 = _mm_add_ss(c20_2, _mm_mul_ss(a20_2, b20));
_mm_store_ss(&C[(i*88)+69], c20_2);
#else
C[(i*88)+20] += values[110] * B[(i*88)+20];
C[(i*88)+41] += values[111] * B[(i*88)+20];
C[(i*88)+69] += values[112] * B[(i*88)+20];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b21 = _mm_broadcast_ss(&B[(i*88)+21]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b21 = _mm_load_ss(&B[(i*88)+21]);
b21 = _mm_shuffle_ps(b21, b21, 0x00);
#endif
__m128 c21_0 = _mm_load_ss(&C[(i*88)+21]);
__m128 a21_0 = _mm_load_ss(&values[113]);
c21_0 = _mm_add_ss(c21_0, _mm_mul_ss(a21_0, b21));
_mm_store_ss(&C[(i*88)+21], c21_0);
__m128 c21_1 = _mm_load_ss(&C[(i*88)+42]);
__m128 a21_1 = _mm_load_ss(&values[114]);
c21_1 = _mm_add_ss(c21_1, _mm_mul_ss(a21_1, b21));
_mm_store_ss(&C[(i*88)+42], c21_1);
__m128 c21_2 = _mm_load_ss(&C[(i*88)+70]);
__m128 a21_2 = _mm_load_ss(&values[115]);
c21_2 = _mm_add_ss(c21_2, _mm_mul_ss(a21_2, b21));
_mm_store_ss(&C[(i*88)+70], c21_2);
#else
C[(i*88)+21] += values[113] * B[(i*88)+21];
C[(i*88)+42] += values[114] * B[(i*88)+21];
C[(i*88)+70] += values[115] * B[(i*88)+21];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b22 = _mm_broadcast_ss(&B[(i*88)+22]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b22 = _mm_load_ss(&B[(i*88)+22]);
b22 = _mm_shuffle_ps(b22, b22, 0x00);
#endif
__m128 c22_0 = _mm_load_ss(&C[(i*88)+22]);
__m128 a22_0 = _mm_load_ss(&values[116]);
c22_0 = _mm_add_ss(c22_0, _mm_mul_ss(a22_0, b22));
_mm_store_ss(&C[(i*88)+22], c22_0);
__m128 c22_1 = _mm_load_ss(&C[(i*88)+43]);
__m128 a22_1 = _mm_load_ss(&values[117]);
c22_1 = _mm_add_ss(c22_1, _mm_mul_ss(a22_1, b22));
_mm_store_ss(&C[(i*88)+43], c22_1);
__m128 c22_2 = _mm_load_ss(&C[(i*88)+71]);
__m128 a22_2 = _mm_load_ss(&values[118]);
c22_2 = _mm_add_ss(c22_2, _mm_mul_ss(a22_2, b22));
_mm_store_ss(&C[(i*88)+71], c22_2);
#else
C[(i*88)+22] += values[116] * B[(i*88)+22];
C[(i*88)+43] += values[117] * B[(i*88)+22];
C[(i*88)+71] += values[118] * B[(i*88)+22];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b23 = _mm_broadcast_ss(&B[(i*88)+23]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b23 = _mm_load_ss(&B[(i*88)+23]);
b23 = _mm_shuffle_ps(b23, b23, 0x00);
#endif
__m128 c23_0 = _mm_load_ss(&C[(i*88)+23]);
__m128 a23_0 = _mm_load_ss(&values[119]);
c23_0 = _mm_add_ss(c23_0, _mm_mul_ss(a23_0, b23));
_mm_store_ss(&C[(i*88)+23], c23_0);
__m128 c23_1 = _mm_load_ss(&C[(i*88)+44]);
__m128 a23_1 = _mm_load_ss(&values[120]);
c23_1 = _mm_add_ss(c23_1, _mm_mul_ss(a23_1, b23));
_mm_store_ss(&C[(i*88)+44], c23_1);
__m128 c23_2 = _mm_load_ss(&C[(i*88)+72]);
__m128 a23_2 = _mm_load_ss(&values[121]);
c23_2 = _mm_add_ss(c23_2, _mm_mul_ss(a23_2, b23));
_mm_store_ss(&C[(i*88)+72], c23_2);
#else
C[(i*88)+23] += values[119] * B[(i*88)+23];
C[(i*88)+44] += values[120] * B[(i*88)+23];
C[(i*88)+72] += values[121] * B[(i*88)+23];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b24 = _mm_broadcast_ss(&B[(i*88)+24]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b24 = _mm_load_ss(&B[(i*88)+24]);
b24 = _mm_shuffle_ps(b24, b24, 0x00);
#endif
__m128 c24_0 = _mm_load_ss(&C[(i*88)+24]);
__m128 a24_0 = _mm_load_ss(&values[122]);
c24_0 = _mm_add_ss(c24_0, _mm_mul_ss(a24_0, b24));
_mm_store_ss(&C[(i*88)+24], c24_0);
__m128 c24_1 = _mm_load_ss(&C[(i*88)+45]);
__m128 a24_1 = _mm_load_ss(&values[123]);
c24_1 = _mm_add_ss(c24_1, _mm_mul_ss(a24_1, b24));
_mm_store_ss(&C[(i*88)+45], c24_1);
__m128 c24_2 = _mm_load_ss(&C[(i*88)+73]);
__m128 a24_2 = _mm_load_ss(&values[124]);
c24_2 = _mm_add_ss(c24_2, _mm_mul_ss(a24_2, b24));
_mm_store_ss(&C[(i*88)+73], c24_2);
#else
C[(i*88)+24] += values[122] * B[(i*88)+24];
C[(i*88)+45] += values[123] * B[(i*88)+24];
C[(i*88)+73] += values[124] * B[(i*88)+24];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b25 = _mm_broadcast_ss(&B[(i*88)+25]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b25 = _mm_load_ss(&B[(i*88)+25]);
b25 = _mm_shuffle_ps(b25, b25, 0x00);
#endif
__m128 c25_0 = _mm_load_ss(&C[(i*88)+10]);
__m128 a25_0 = _mm_load_ss(&values[125]);
c25_0 = _mm_add_ss(c25_0, _mm_mul_ss(a25_0, b25));
_mm_store_ss(&C[(i*88)+10], c25_0);
__m128 c25_1 = _mm_load_ss(&C[(i*88)+25]);
__m128 a25_1 = _mm_load_ss(&values[126]);
c25_1 = _mm_add_ss(c25_1, _mm_mul_ss(a25_1, b25));
_mm_store_ss(&C[(i*88)+25], c25_1);
__m128 c25_2 = _mm_load_ss(&C[(i*88)+46]);
__m128 a25_2 = _mm_load_ss(&values[127]);
c25_2 = _mm_add_ss(c25_2, _mm_mul_ss(a25_2, b25));
_mm_store_ss(&C[(i*88)+46], c25_2);
__m128 c25_3 = _mm_load_ss(&C[(i*88)+74]);
__m128 a25_3 = _mm_load_ss(&values[128]);
c25_3 = _mm_add_ss(c25_3, _mm_mul_ss(a25_3, b25));
_mm_store_ss(&C[(i*88)+74], c25_3);
#else
C[(i*88)+10] += values[125] * B[(i*88)+25];
C[(i*88)+25] += values[126] * B[(i*88)+25];
C[(i*88)+46] += values[127] * B[(i*88)+25];
C[(i*88)+74] += values[128] * B[(i*88)+25];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b26 = _mm_broadcast_ss(&B[(i*88)+26]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b26 = _mm_load_ss(&B[(i*88)+26]);
b26 = _mm_shuffle_ps(b26, b26, 0x00);
#endif
__m128 c26_0 = _mm_load_ss(&C[(i*88)+11]);
__m128 a26_0 = _mm_load_ss(&values[129]);
c26_0 = _mm_add_ss(c26_0, _mm_mul_ss(a26_0, b26));
_mm_store_ss(&C[(i*88)+11], c26_0);
__m128 c26_1 = _mm_load_ss(&C[(i*88)+26]);
__m128 a26_1 = _mm_load_ss(&values[130]);
c26_1 = _mm_add_ss(c26_1, _mm_mul_ss(a26_1, b26));
_mm_store_ss(&C[(i*88)+26], c26_1);
__m128 c26_2 = _mm_load_ss(&C[(i*88)+47]);
__m128 a26_2 = _mm_load_ss(&values[131]);
c26_2 = _mm_add_ss(c26_2, _mm_mul_ss(a26_2, b26));
_mm_store_ss(&C[(i*88)+47], c26_2);
__m128 c26_3 = _mm_load_ss(&C[(i*88)+75]);
__m128 a26_3 = _mm_load_ss(&values[132]);
c26_3 = _mm_add_ss(c26_3, _mm_mul_ss(a26_3, b26));
_mm_store_ss(&C[(i*88)+75], c26_3);
#else
C[(i*88)+11] += values[129] * B[(i*88)+26];
C[(i*88)+26] += values[130] * B[(i*88)+26];
C[(i*88)+47] += values[131] * B[(i*88)+26];
C[(i*88)+75] += values[132] * B[(i*88)+26];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b27 = _mm_broadcast_ss(&B[(i*88)+27]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b27 = _mm_load_ss(&B[(i*88)+27]);
b27 = _mm_shuffle_ps(b27, b27, 0x00);
#endif
__m128 c27_0 = _mm_load_ss(&C[(i*88)+12]);
__m128 a27_0 = _mm_load_ss(&values[133]);
c27_0 = _mm_add_ss(c27_0, _mm_mul_ss(a27_0, b27));
_mm_store_ss(&C[(i*88)+12], c27_0);
__m128 c27_1 = _mm_load_ss(&C[(i*88)+27]);
__m128 a27_1 = _mm_load_ss(&values[134]);
c27_1 = _mm_add_ss(c27_1, _mm_mul_ss(a27_1, b27));
_mm_store_ss(&C[(i*88)+27], c27_1);
__m128 c27_2 = _mm_load_ss(&C[(i*88)+48]);
__m128 a27_2 = _mm_load_ss(&values[135]);
c27_2 = _mm_add_ss(c27_2, _mm_mul_ss(a27_2, b27));
_mm_store_ss(&C[(i*88)+48], c27_2);
__m128 c27_3 = _mm_load_ss(&C[(i*88)+76]);
__m128 a27_3 = _mm_load_ss(&values[136]);
c27_3 = _mm_add_ss(c27_3, _mm_mul_ss(a27_3, b27));
_mm_store_ss(&C[(i*88)+76], c27_3);
#else
C[(i*88)+12] += values[133] * B[(i*88)+27];
C[(i*88)+27] += values[134] * B[(i*88)+27];
C[(i*88)+48] += values[135] * B[(i*88)+27];
C[(i*88)+76] += values[136] * B[(i*88)+27];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b28 = _mm_broadcast_ss(&B[(i*88)+28]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b28 = _mm_load_ss(&B[(i*88)+28]);
b28 = _mm_shuffle_ps(b28, b28, 0x00);
#endif
__m128 c28_0 = _mm_load_ss(&C[(i*88)+13]);
__m128 a28_0 = _mm_load_ss(&values[137]);
c28_0 = _mm_add_ss(c28_0, _mm_mul_ss(a28_0, b28));
_mm_store_ss(&C[(i*88)+13], c28_0);
__m128 c28_1 = _mm_load_ss(&C[(i*88)+28]);
__m128 a28_1 = _mm_load_ss(&values[138]);
c28_1 = _mm_add_ss(c28_1, _mm_mul_ss(a28_1, b28));
_mm_store_ss(&C[(i*88)+28], c28_1);
__m128 c28_2 = _mm_load_ss(&C[(i*88)+49]);
__m128 a28_2 = _mm_load_ss(&values[139]);
c28_2 = _mm_add_ss(c28_2, _mm_mul_ss(a28_2, b28));
_mm_store_ss(&C[(i*88)+49], c28_2);
__m128 c28_3 = _mm_load_ss(&C[(i*88)+77]);
__m128 a28_3 = _mm_load_ss(&values[140]);
c28_3 = _mm_add_ss(c28_3, _mm_mul_ss(a28_3, b28));
_mm_store_ss(&C[(i*88)+77], c28_3);
#else
C[(i*88)+13] += values[137] * B[(i*88)+28];
C[(i*88)+28] += values[138] * B[(i*88)+28];
C[(i*88)+49] += values[139] * B[(i*88)+28];
C[(i*88)+77] += values[140] * B[(i*88)+28];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b29 = _mm_broadcast_ss(&B[(i*88)+29]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b29 = _mm_load_ss(&B[(i*88)+29]);
b29 = _mm_shuffle_ps(b29, b29, 0x00);
#endif
__m128 c29_0 = _mm_load_ss(&C[(i*88)+4]);
__m128 a29_0 = _mm_load_ss(&values[141]);
c29_0 = _mm_add_ss(c29_0, _mm_mul_ss(a29_0, b29));
_mm_store_ss(&C[(i*88)+4], c29_0);
__m128 c29_1 = _mm_load_ss(&C[(i*88)+14]);
__m128 a29_1 = _mm_load_ss(&values[142]);
c29_1 = _mm_add_ss(c29_1, _mm_mul_ss(a29_1, b29));
_mm_store_ss(&C[(i*88)+14], c29_1);
__m128 c29_2 = _mm_load_ss(&C[(i*88)+29]);
__m128 a29_2 = _mm_load_ss(&values[143]);
c29_2 = _mm_add_ss(c29_2, _mm_mul_ss(a29_2, b29));
_mm_store_ss(&C[(i*88)+29], c29_2);
__m128 c29_3 = _mm_load_ss(&C[(i*88)+50]);
__m128 a29_3 = _mm_load_ss(&values[144]);
c29_3 = _mm_add_ss(c29_3, _mm_mul_ss(a29_3, b29));
_mm_store_ss(&C[(i*88)+50], c29_3);
__m128 c29_4 = _mm_load_ss(&C[(i*88)+78]);
__m128 a29_4 = _mm_load_ss(&values[145]);
c29_4 = _mm_add_ss(c29_4, _mm_mul_ss(a29_4, b29));
_mm_store_ss(&C[(i*88)+78], c29_4);
#else
C[(i*88)+4] += values[141] * B[(i*88)+29];
C[(i*88)+14] += values[142] * B[(i*88)+29];
C[(i*88)+29] += values[143] * B[(i*88)+29];
C[(i*88)+50] += values[144] * B[(i*88)+29];
C[(i*88)+78] += values[145] * B[(i*88)+29];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b30 = _mm_broadcast_ss(&B[(i*88)+30]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b30 = _mm_load_ss(&B[(i*88)+30]);
b30 = _mm_shuffle_ps(b30, b30, 0x00);
#endif
__m128 c30_0 = _mm_load_ss(&C[(i*88)+5]);
__m128 a30_0 = _mm_load_ss(&values[146]);
c30_0 = _mm_add_ss(c30_0, _mm_mul_ss(a30_0, b30));
_mm_store_ss(&C[(i*88)+5], c30_0);
__m128 c30_1 = _mm_load_ss(&C[(i*88)+15]);
__m128 a30_1 = _mm_load_ss(&values[147]);
c30_1 = _mm_add_ss(c30_1, _mm_mul_ss(a30_1, b30));
_mm_store_ss(&C[(i*88)+15], c30_1);
__m128 c30_2 = _mm_load_ss(&C[(i*88)+30]);
__m128 a30_2 = _mm_load_ss(&values[148]);
c30_2 = _mm_add_ss(c30_2, _mm_mul_ss(a30_2, b30));
_mm_store_ss(&C[(i*88)+30], c30_2);
__m128 c30_3 = _mm_load_ss(&C[(i*88)+51]);
__m128 a30_3 = _mm_load_ss(&values[149]);
c30_3 = _mm_add_ss(c30_3, _mm_mul_ss(a30_3, b30));
_mm_store_ss(&C[(i*88)+51], c30_3);
__m128 c30_4 = _mm_load_ss(&C[(i*88)+79]);
__m128 a30_4 = _mm_load_ss(&values[150]);
c30_4 = _mm_add_ss(c30_4, _mm_mul_ss(a30_4, b30));
_mm_store_ss(&C[(i*88)+79], c30_4);
#else
C[(i*88)+5] += values[146] * B[(i*88)+30];
C[(i*88)+15] += values[147] * B[(i*88)+30];
C[(i*88)+30] += values[148] * B[(i*88)+30];
C[(i*88)+51] += values[149] * B[(i*88)+30];
C[(i*88)+79] += values[150] * B[(i*88)+30];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b31 = _mm_broadcast_ss(&B[(i*88)+31]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b31 = _mm_load_ss(&B[(i*88)+31]);
b31 = _mm_shuffle_ps(b31, b31, 0x00);
#endif
__m128 c31_0 = _mm_load_ss(&C[(i*88)+6]);
__m128 a31_0 = _mm_load_ss(&values[151]);
c31_0 = _mm_add_ss(c31_0, _mm_mul_ss(a31_0, b31));
_mm_store_ss(&C[(i*88)+6], c31_0);
__m128 c31_1 = _mm_load_ss(&C[(i*88)+16]);
__m128 a31_1 = _mm_load_ss(&values[152]);
c31_1 = _mm_add_ss(c31_1, _mm_mul_ss(a31_1, b31));
_mm_store_ss(&C[(i*88)+16], c31_1);
__m128 c31_2 = _mm_load_ss(&C[(i*88)+31]);
__m128 a31_2 = _mm_load_ss(&values[153]);
c31_2 = _mm_add_ss(c31_2, _mm_mul_ss(a31_2, b31));
_mm_store_ss(&C[(i*88)+31], c31_2);
__m128 c31_3 = _mm_load_ss(&C[(i*88)+52]);
__m128 a31_3 = _mm_load_ss(&values[154]);
c31_3 = _mm_add_ss(c31_3, _mm_mul_ss(a31_3, b31));
_mm_store_ss(&C[(i*88)+52], c31_3);
__m128 c31_4 = _mm_load_ss(&C[(i*88)+80]);
__m128 a31_4 = _mm_load_ss(&values[155]);
c31_4 = _mm_add_ss(c31_4, _mm_mul_ss(a31_4, b31));
_mm_store_ss(&C[(i*88)+80], c31_4);
#else
C[(i*88)+6] += values[151] * B[(i*88)+31];
C[(i*88)+16] += values[152] * B[(i*88)+31];
C[(i*88)+31] += values[153] * B[(i*88)+31];
C[(i*88)+52] += values[154] * B[(i*88)+31];
C[(i*88)+80] += values[155] * B[(i*88)+31];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b32 = _mm_broadcast_ss(&B[(i*88)+32]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b32 = _mm_load_ss(&B[(i*88)+32]);
b32 = _mm_shuffle_ps(b32, b32, 0x00);
#endif
__m128 c32_0 = _mm_load_ss(&C[(i*88)+1]);
__m128 a32_0 = _mm_load_ss(&values[156]);
c32_0 = _mm_add_ss(c32_0, _mm_mul_ss(a32_0, b32));
_mm_store_ss(&C[(i*88)+1], c32_0);
__m128 c32_1 = _mm_load_ss(&C[(i*88)+7]);
__m128 a32_1 = _mm_load_ss(&values[157]);
c32_1 = _mm_add_ss(c32_1, _mm_mul_ss(a32_1, b32));
_mm_store_ss(&C[(i*88)+7], c32_1);
__m128 c32_2 = _mm_load_ss(&C[(i*88)+17]);
__m128 a32_2 = _mm_load_ss(&values[158]);
c32_2 = _mm_add_ss(c32_2, _mm_mul_ss(a32_2, b32));
_mm_store_ss(&C[(i*88)+17], c32_2);
__m128 c32_3 = _mm_load_ss(&C[(i*88)+32]);
__m128 a32_3 = _mm_load_ss(&values[159]);
c32_3 = _mm_add_ss(c32_3, _mm_mul_ss(a32_3, b32));
_mm_store_ss(&C[(i*88)+32], c32_3);
__m128 c32_4 = _mm_load_ss(&C[(i*88)+53]);
__m128 a32_4 = _mm_load_ss(&values[160]);
c32_4 = _mm_add_ss(c32_4, _mm_mul_ss(a32_4, b32));
_mm_store_ss(&C[(i*88)+53], c32_4);
__m128 c32_5 = _mm_load_ss(&C[(i*88)+81]);
__m128 a32_5 = _mm_load_ss(&values[161]);
c32_5 = _mm_add_ss(c32_5, _mm_mul_ss(a32_5, b32));
_mm_store_ss(&C[(i*88)+81], c32_5);
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
__m128 b33 = _mm_broadcast_ss(&B[(i*88)+33]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b33 = _mm_load_ss(&B[(i*88)+33]);
b33 = _mm_shuffle_ps(b33, b33, 0x00);
#endif
__m128 c33_0 = _mm_load_ss(&C[(i*88)+2]);
__m128 a33_0 = _mm_load_ss(&values[162]);
c33_0 = _mm_add_ss(c33_0, _mm_mul_ss(a33_0, b33));
_mm_store_ss(&C[(i*88)+2], c33_0);
__m128 c33_1 = _mm_load_ss(&C[(i*88)+8]);
__m128 a33_1 = _mm_load_ss(&values[163]);
c33_1 = _mm_add_ss(c33_1, _mm_mul_ss(a33_1, b33));
_mm_store_ss(&C[(i*88)+8], c33_1);
__m128 c33_2 = _mm_load_ss(&C[(i*88)+18]);
__m128 a33_2 = _mm_load_ss(&values[164]);
c33_2 = _mm_add_ss(c33_2, _mm_mul_ss(a33_2, b33));
_mm_store_ss(&C[(i*88)+18], c33_2);
__m128 c33_3 = _mm_load_ss(&C[(i*88)+33]);
__m128 a33_3 = _mm_load_ss(&values[165]);
c33_3 = _mm_add_ss(c33_3, _mm_mul_ss(a33_3, b33));
_mm_store_ss(&C[(i*88)+33], c33_3);
__m128 c33_4 = _mm_load_ss(&C[(i*88)+54]);
__m128 a33_4 = _mm_load_ss(&values[166]);
c33_4 = _mm_add_ss(c33_4, _mm_mul_ss(a33_4, b33));
_mm_store_ss(&C[(i*88)+54], c33_4);
__m128 c33_5 = _mm_load_ss(&C[(i*88)+82]);
__m128 a33_5 = _mm_load_ss(&values[167]);
c33_5 = _mm_add_ss(c33_5, _mm_mul_ss(a33_5, b33));
_mm_store_ss(&C[(i*88)+82], c33_5);
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
__m128 b34 = _mm_broadcast_ss(&B[(i*88)+34]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b34 = _mm_load_ss(&B[(i*88)+34]);
b34 = _mm_shuffle_ps(b34, b34, 0x00);
#endif
__m128 c34_0 = _mm_load_ss(&C[(i*88)+0]);
__m128 a34_0 = _mm_load_ss(&values[168]);
c34_0 = _mm_add_ss(c34_0, _mm_mul_ss(a34_0, b34));
_mm_store_ss(&C[(i*88)+0], c34_0);
__m128 c34_1 = _mm_load_ss(&C[(i*88)+3]);
__m128 a34_1 = _mm_load_ss(&values[169]);
c34_1 = _mm_add_ss(c34_1, _mm_mul_ss(a34_1, b34));
_mm_store_ss(&C[(i*88)+3], c34_1);
__m128 c34_2 = _mm_load_ss(&C[(i*88)+9]);
__m128 a34_2 = _mm_load_ss(&values[170]);
c34_2 = _mm_add_ss(c34_2, _mm_mul_ss(a34_2, b34));
_mm_store_ss(&C[(i*88)+9], c34_2);
__m128 c34_3 = _mm_load_ss(&C[(i*88)+19]);
__m128 a34_3 = _mm_load_ss(&values[171]);
c34_3 = _mm_add_ss(c34_3, _mm_mul_ss(a34_3, b34));
_mm_store_ss(&C[(i*88)+19], c34_3);
__m128 c34_4 = _mm_load_ss(&C[(i*88)+34]);
__m128 a34_4 = _mm_load_ss(&values[172]);
c34_4 = _mm_add_ss(c34_4, _mm_mul_ss(a34_4, b34));
_mm_store_ss(&C[(i*88)+34], c34_4);
__m128 c34_5 = _mm_load_ss(&C[(i*88)+55]);
__m128 a34_5 = _mm_load_ss(&values[173]);
c34_5 = _mm_add_ss(c34_5, _mm_mul_ss(a34_5, b34));
_mm_store_ss(&C[(i*88)+55], c34_5);
__m128 c34_6 = _mm_load_ss(&C[(i*88)+83]);
__m128 a34_6 = _mm_load_ss(&values[174]);
c34_6 = _mm_add_ss(c34_6, _mm_mul_ss(a34_6, b34));
_mm_store_ss(&C[(i*88)+83], c34_6);
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
__m128 b35 = _mm_broadcast_ss(&B[(i*88)+35]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b35 = _mm_load_ss(&B[(i*88)+35]);
b35 = _mm_shuffle_ps(b35, b35, 0x00);
#endif
__m128 c35_0 = _mm_load_ss(&C[(i*88)+35]);
__m128 a35_0 = _mm_load_ss(&values[175]);
c35_0 = _mm_add_ss(c35_0, _mm_mul_ss(a35_0, b35));
_mm_store_ss(&C[(i*88)+35], c35_0);
__m128 c35_1 = _mm_load_ss(&C[(i*88)+63]);
__m128 a35_1 = _mm_load_ss(&values[176]);
c35_1 = _mm_add_ss(c35_1, _mm_mul_ss(a35_1, b35));
_mm_store_ss(&C[(i*88)+63], c35_1);
#else
C[(i*88)+35] += values[175] * B[(i*88)+35];
C[(i*88)+63] += values[176] * B[(i*88)+35];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b36 = _mm_broadcast_ss(&B[(i*88)+36]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b36 = _mm_load_ss(&B[(i*88)+36]);
b36 = _mm_shuffle_ps(b36, b36, 0x00);
#endif
__m128 c36_0 = _mm_load_ss(&C[(i*88)+36]);
__m128 a36_0 = _mm_load_ss(&values[177]);
c36_0 = _mm_add_ss(c36_0, _mm_mul_ss(a36_0, b36));
_mm_store_ss(&C[(i*88)+36], c36_0);
__m128 c36_1 = _mm_load_ss(&C[(i*88)+64]);
__m128 a36_1 = _mm_load_ss(&values[178]);
c36_1 = _mm_add_ss(c36_1, _mm_mul_ss(a36_1, b36));
_mm_store_ss(&C[(i*88)+64], c36_1);
#else
C[(i*88)+36] += values[177] * B[(i*88)+36];
C[(i*88)+64] += values[178] * B[(i*88)+36];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b37 = _mm_broadcast_ss(&B[(i*88)+37]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b37 = _mm_load_ss(&B[(i*88)+37]);
b37 = _mm_shuffle_ps(b37, b37, 0x00);
#endif
__m128 c37_0 = _mm_load_ss(&C[(i*88)+37]);
__m128 a37_0 = _mm_load_ss(&values[179]);
c37_0 = _mm_add_ss(c37_0, _mm_mul_ss(a37_0, b37));
_mm_store_ss(&C[(i*88)+37], c37_0);
__m128 c37_1 = _mm_load_ss(&C[(i*88)+65]);
__m128 a37_1 = _mm_load_ss(&values[180]);
c37_1 = _mm_add_ss(c37_1, _mm_mul_ss(a37_1, b37));
_mm_store_ss(&C[(i*88)+65], c37_1);
#else
C[(i*88)+37] += values[179] * B[(i*88)+37];
C[(i*88)+65] += values[180] * B[(i*88)+37];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b38 = _mm_broadcast_ss(&B[(i*88)+38]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b38 = _mm_load_ss(&B[(i*88)+38]);
b38 = _mm_shuffle_ps(b38, b38, 0x00);
#endif
__m128 c38_0 = _mm_load_ss(&C[(i*88)+38]);
__m128 a38_0 = _mm_load_ss(&values[181]);
c38_0 = _mm_add_ss(c38_0, _mm_mul_ss(a38_0, b38));
_mm_store_ss(&C[(i*88)+38], c38_0);
__m128 c38_1 = _mm_load_ss(&C[(i*88)+66]);
__m128 a38_1 = _mm_load_ss(&values[182]);
c38_1 = _mm_add_ss(c38_1, _mm_mul_ss(a38_1, b38));
_mm_store_ss(&C[(i*88)+66], c38_1);
#else
C[(i*88)+38] += values[181] * B[(i*88)+38];
C[(i*88)+66] += values[182] * B[(i*88)+38];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b39 = _mm_broadcast_ss(&B[(i*88)+39]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b39 = _mm_load_ss(&B[(i*88)+39]);
b39 = _mm_shuffle_ps(b39, b39, 0x00);
#endif
__m128 c39_0 = _mm_load_ss(&C[(i*88)+39]);
__m128 a39_0 = _mm_load_ss(&values[183]);
c39_0 = _mm_add_ss(c39_0, _mm_mul_ss(a39_0, b39));
_mm_store_ss(&C[(i*88)+39], c39_0);
__m128 c39_1 = _mm_load_ss(&C[(i*88)+67]);
__m128 a39_1 = _mm_load_ss(&values[184]);
c39_1 = _mm_add_ss(c39_1, _mm_mul_ss(a39_1, b39));
_mm_store_ss(&C[(i*88)+67], c39_1);
#else
C[(i*88)+39] += values[183] * B[(i*88)+39];
C[(i*88)+67] += values[184] * B[(i*88)+39];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b40 = _mm_broadcast_ss(&B[(i*88)+40]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b40 = _mm_load_ss(&B[(i*88)+40]);
b40 = _mm_shuffle_ps(b40, b40, 0x00);
#endif
__m128 c40_0 = _mm_load_ss(&C[(i*88)+40]);
__m128 a40_0 = _mm_load_ss(&values[185]);
c40_0 = _mm_add_ss(c40_0, _mm_mul_ss(a40_0, b40));
_mm_store_ss(&C[(i*88)+40], c40_0);
__m128 c40_1 = _mm_load_ss(&C[(i*88)+68]);
__m128 a40_1 = _mm_load_ss(&values[186]);
c40_1 = _mm_add_ss(c40_1, _mm_mul_ss(a40_1, b40));
_mm_store_ss(&C[(i*88)+68], c40_1);
#else
C[(i*88)+40] += values[185] * B[(i*88)+40];
C[(i*88)+68] += values[186] * B[(i*88)+40];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b41 = _mm_broadcast_ss(&B[(i*88)+41]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b41 = _mm_load_ss(&B[(i*88)+41]);
b41 = _mm_shuffle_ps(b41, b41, 0x00);
#endif
__m128 c41_0 = _mm_load_ss(&C[(i*88)+20]);
__m128 a41_0 = _mm_load_ss(&values[187]);
c41_0 = _mm_add_ss(c41_0, _mm_mul_ss(a41_0, b41));
_mm_store_ss(&C[(i*88)+20], c41_0);
__m128 c41_1 = _mm_load_ss(&C[(i*88)+41]);
__m128 a41_1 = _mm_load_ss(&values[188]);
c41_1 = _mm_add_ss(c41_1, _mm_mul_ss(a41_1, b41));
_mm_store_ss(&C[(i*88)+41], c41_1);
__m128 c41_2 = _mm_load_ss(&C[(i*88)+69]);
__m128 a41_2 = _mm_load_ss(&values[189]);
c41_2 = _mm_add_ss(c41_2, _mm_mul_ss(a41_2, b41));
_mm_store_ss(&C[(i*88)+69], c41_2);
#else
C[(i*88)+20] += values[187] * B[(i*88)+41];
C[(i*88)+41] += values[188] * B[(i*88)+41];
C[(i*88)+69] += values[189] * B[(i*88)+41];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b42 = _mm_broadcast_ss(&B[(i*88)+42]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b42 = _mm_load_ss(&B[(i*88)+42]);
b42 = _mm_shuffle_ps(b42, b42, 0x00);
#endif
__m128 c42_0 = _mm_load_ss(&C[(i*88)+21]);
__m128 a42_0 = _mm_load_ss(&values[190]);
c42_0 = _mm_add_ss(c42_0, _mm_mul_ss(a42_0, b42));
_mm_store_ss(&C[(i*88)+21], c42_0);
__m128 c42_1 = _mm_load_ss(&C[(i*88)+42]);
__m128 a42_1 = _mm_load_ss(&values[191]);
c42_1 = _mm_add_ss(c42_1, _mm_mul_ss(a42_1, b42));
_mm_store_ss(&C[(i*88)+42], c42_1);
__m128 c42_2 = _mm_load_ss(&C[(i*88)+70]);
__m128 a42_2 = _mm_load_ss(&values[192]);
c42_2 = _mm_add_ss(c42_2, _mm_mul_ss(a42_2, b42));
_mm_store_ss(&C[(i*88)+70], c42_2);
#else
C[(i*88)+21] += values[190] * B[(i*88)+42];
C[(i*88)+42] += values[191] * B[(i*88)+42];
C[(i*88)+70] += values[192] * B[(i*88)+42];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b43 = _mm_broadcast_ss(&B[(i*88)+43]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b43 = _mm_load_ss(&B[(i*88)+43]);
b43 = _mm_shuffle_ps(b43, b43, 0x00);
#endif
__m128 c43_0 = _mm_load_ss(&C[(i*88)+22]);
__m128 a43_0 = _mm_load_ss(&values[193]);
c43_0 = _mm_add_ss(c43_0, _mm_mul_ss(a43_0, b43));
_mm_store_ss(&C[(i*88)+22], c43_0);
__m128 c43_1 = _mm_load_ss(&C[(i*88)+43]);
__m128 a43_1 = _mm_load_ss(&values[194]);
c43_1 = _mm_add_ss(c43_1, _mm_mul_ss(a43_1, b43));
_mm_store_ss(&C[(i*88)+43], c43_1);
__m128 c43_2 = _mm_load_ss(&C[(i*88)+71]);
__m128 a43_2 = _mm_load_ss(&values[195]);
c43_2 = _mm_add_ss(c43_2, _mm_mul_ss(a43_2, b43));
_mm_store_ss(&C[(i*88)+71], c43_2);
#else
C[(i*88)+22] += values[193] * B[(i*88)+43];
C[(i*88)+43] += values[194] * B[(i*88)+43];
C[(i*88)+71] += values[195] * B[(i*88)+43];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b44 = _mm_broadcast_ss(&B[(i*88)+44]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b44 = _mm_load_ss(&B[(i*88)+44]);
b44 = _mm_shuffle_ps(b44, b44, 0x00);
#endif
__m128 c44_0 = _mm_load_ss(&C[(i*88)+23]);
__m128 a44_0 = _mm_load_ss(&values[196]);
c44_0 = _mm_add_ss(c44_0, _mm_mul_ss(a44_0, b44));
_mm_store_ss(&C[(i*88)+23], c44_0);
__m128 c44_1 = _mm_load_ss(&C[(i*88)+44]);
__m128 a44_1 = _mm_load_ss(&values[197]);
c44_1 = _mm_add_ss(c44_1, _mm_mul_ss(a44_1, b44));
_mm_store_ss(&C[(i*88)+44], c44_1);
__m128 c44_2 = _mm_load_ss(&C[(i*88)+72]);
__m128 a44_2 = _mm_load_ss(&values[198]);
c44_2 = _mm_add_ss(c44_2, _mm_mul_ss(a44_2, b44));
_mm_store_ss(&C[(i*88)+72], c44_2);
#else
C[(i*88)+23] += values[196] * B[(i*88)+44];
C[(i*88)+44] += values[197] * B[(i*88)+44];
C[(i*88)+72] += values[198] * B[(i*88)+44];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b45 = _mm_broadcast_ss(&B[(i*88)+45]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b45 = _mm_load_ss(&B[(i*88)+45]);
b45 = _mm_shuffle_ps(b45, b45, 0x00);
#endif
__m128 c45_0 = _mm_load_ss(&C[(i*88)+24]);
__m128 a45_0 = _mm_load_ss(&values[199]);
c45_0 = _mm_add_ss(c45_0, _mm_mul_ss(a45_0, b45));
_mm_store_ss(&C[(i*88)+24], c45_0);
__m128 c45_1 = _mm_load_ss(&C[(i*88)+45]);
__m128 a45_1 = _mm_load_ss(&values[200]);
c45_1 = _mm_add_ss(c45_1, _mm_mul_ss(a45_1, b45));
_mm_store_ss(&C[(i*88)+45], c45_1);
__m128 c45_2 = _mm_load_ss(&C[(i*88)+73]);
__m128 a45_2 = _mm_load_ss(&values[201]);
c45_2 = _mm_add_ss(c45_2, _mm_mul_ss(a45_2, b45));
_mm_store_ss(&C[(i*88)+73], c45_2);
#else
C[(i*88)+24] += values[199] * B[(i*88)+45];
C[(i*88)+45] += values[200] * B[(i*88)+45];
C[(i*88)+73] += values[201] * B[(i*88)+45];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b46 = _mm_broadcast_ss(&B[(i*88)+46]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b46 = _mm_load_ss(&B[(i*88)+46]);
b46 = _mm_shuffle_ps(b46, b46, 0x00);
#endif
__m128 c46_0 = _mm_load_ss(&C[(i*88)+10]);
__m128 a46_0 = _mm_load_ss(&values[202]);
c46_0 = _mm_add_ss(c46_0, _mm_mul_ss(a46_0, b46));
_mm_store_ss(&C[(i*88)+10], c46_0);
__m128 c46_1 = _mm_load_ss(&C[(i*88)+25]);
__m128 a46_1 = _mm_load_ss(&values[203]);
c46_1 = _mm_add_ss(c46_1, _mm_mul_ss(a46_1, b46));
_mm_store_ss(&C[(i*88)+25], c46_1);
__m128 c46_2 = _mm_load_ss(&C[(i*88)+46]);
__m128 a46_2 = _mm_load_ss(&values[204]);
c46_2 = _mm_add_ss(c46_2, _mm_mul_ss(a46_2, b46));
_mm_store_ss(&C[(i*88)+46], c46_2);
__m128 c46_3 = _mm_load_ss(&C[(i*88)+74]);
__m128 a46_3 = _mm_load_ss(&values[205]);
c46_3 = _mm_add_ss(c46_3, _mm_mul_ss(a46_3, b46));
_mm_store_ss(&C[(i*88)+74], c46_3);
#else
C[(i*88)+10] += values[202] * B[(i*88)+46];
C[(i*88)+25] += values[203] * B[(i*88)+46];
C[(i*88)+46] += values[204] * B[(i*88)+46];
C[(i*88)+74] += values[205] * B[(i*88)+46];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b47 = _mm_broadcast_ss(&B[(i*88)+47]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b47 = _mm_load_ss(&B[(i*88)+47]);
b47 = _mm_shuffle_ps(b47, b47, 0x00);
#endif
__m128 c47_0 = _mm_load_ss(&C[(i*88)+11]);
__m128 a47_0 = _mm_load_ss(&values[206]);
c47_0 = _mm_add_ss(c47_0, _mm_mul_ss(a47_0, b47));
_mm_store_ss(&C[(i*88)+11], c47_0);
__m128 c47_1 = _mm_load_ss(&C[(i*88)+26]);
__m128 a47_1 = _mm_load_ss(&values[207]);
c47_1 = _mm_add_ss(c47_1, _mm_mul_ss(a47_1, b47));
_mm_store_ss(&C[(i*88)+26], c47_1);
__m128 c47_2 = _mm_load_ss(&C[(i*88)+47]);
__m128 a47_2 = _mm_load_ss(&values[208]);
c47_2 = _mm_add_ss(c47_2, _mm_mul_ss(a47_2, b47));
_mm_store_ss(&C[(i*88)+47], c47_2);
__m128 c47_3 = _mm_load_ss(&C[(i*88)+75]);
__m128 a47_3 = _mm_load_ss(&values[209]);
c47_3 = _mm_add_ss(c47_3, _mm_mul_ss(a47_3, b47));
_mm_store_ss(&C[(i*88)+75], c47_3);
#else
C[(i*88)+11] += values[206] * B[(i*88)+47];
C[(i*88)+26] += values[207] * B[(i*88)+47];
C[(i*88)+47] += values[208] * B[(i*88)+47];
C[(i*88)+75] += values[209] * B[(i*88)+47];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b48 = _mm_broadcast_ss(&B[(i*88)+48]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b48 = _mm_load_ss(&B[(i*88)+48]);
b48 = _mm_shuffle_ps(b48, b48, 0x00);
#endif
__m128 c48_0 = _mm_load_ss(&C[(i*88)+12]);
__m128 a48_0 = _mm_load_ss(&values[210]);
c48_0 = _mm_add_ss(c48_0, _mm_mul_ss(a48_0, b48));
_mm_store_ss(&C[(i*88)+12], c48_0);
__m128 c48_1 = _mm_load_ss(&C[(i*88)+27]);
__m128 a48_1 = _mm_load_ss(&values[211]);
c48_1 = _mm_add_ss(c48_1, _mm_mul_ss(a48_1, b48));
_mm_store_ss(&C[(i*88)+27], c48_1);
__m128 c48_2 = _mm_load_ss(&C[(i*88)+48]);
__m128 a48_2 = _mm_load_ss(&values[212]);
c48_2 = _mm_add_ss(c48_2, _mm_mul_ss(a48_2, b48));
_mm_store_ss(&C[(i*88)+48], c48_2);
__m128 c48_3 = _mm_load_ss(&C[(i*88)+76]);
__m128 a48_3 = _mm_load_ss(&values[213]);
c48_3 = _mm_add_ss(c48_3, _mm_mul_ss(a48_3, b48));
_mm_store_ss(&C[(i*88)+76], c48_3);
#else
C[(i*88)+12] += values[210] * B[(i*88)+48];
C[(i*88)+27] += values[211] * B[(i*88)+48];
C[(i*88)+48] += values[212] * B[(i*88)+48];
C[(i*88)+76] += values[213] * B[(i*88)+48];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b49 = _mm_broadcast_ss(&B[(i*88)+49]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b49 = _mm_load_ss(&B[(i*88)+49]);
b49 = _mm_shuffle_ps(b49, b49, 0x00);
#endif
__m128 c49_0 = _mm_load_ss(&C[(i*88)+13]);
__m128 a49_0 = _mm_load_ss(&values[214]);
c49_0 = _mm_add_ss(c49_0, _mm_mul_ss(a49_0, b49));
_mm_store_ss(&C[(i*88)+13], c49_0);
__m128 c49_1 = _mm_load_ss(&C[(i*88)+28]);
__m128 a49_1 = _mm_load_ss(&values[215]);
c49_1 = _mm_add_ss(c49_1, _mm_mul_ss(a49_1, b49));
_mm_store_ss(&C[(i*88)+28], c49_1);
__m128 c49_2 = _mm_load_ss(&C[(i*88)+49]);
__m128 a49_2 = _mm_load_ss(&values[216]);
c49_2 = _mm_add_ss(c49_2, _mm_mul_ss(a49_2, b49));
_mm_store_ss(&C[(i*88)+49], c49_2);
__m128 c49_3 = _mm_load_ss(&C[(i*88)+77]);
__m128 a49_3 = _mm_load_ss(&values[217]);
c49_3 = _mm_add_ss(c49_3, _mm_mul_ss(a49_3, b49));
_mm_store_ss(&C[(i*88)+77], c49_3);
#else
C[(i*88)+13] += values[214] * B[(i*88)+49];
C[(i*88)+28] += values[215] * B[(i*88)+49];
C[(i*88)+49] += values[216] * B[(i*88)+49];
C[(i*88)+77] += values[217] * B[(i*88)+49];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b50 = _mm_broadcast_ss(&B[(i*88)+50]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b50 = _mm_load_ss(&B[(i*88)+50]);
b50 = _mm_shuffle_ps(b50, b50, 0x00);
#endif
__m128 c50_0 = _mm_load_ss(&C[(i*88)+4]);
__m128 a50_0 = _mm_load_ss(&values[218]);
c50_0 = _mm_add_ss(c50_0, _mm_mul_ss(a50_0, b50));
_mm_store_ss(&C[(i*88)+4], c50_0);
__m128 c50_1 = _mm_load_ss(&C[(i*88)+14]);
__m128 a50_1 = _mm_load_ss(&values[219]);
c50_1 = _mm_add_ss(c50_1, _mm_mul_ss(a50_1, b50));
_mm_store_ss(&C[(i*88)+14], c50_1);
__m128 c50_2 = _mm_load_ss(&C[(i*88)+29]);
__m128 a50_2 = _mm_load_ss(&values[220]);
c50_2 = _mm_add_ss(c50_2, _mm_mul_ss(a50_2, b50));
_mm_store_ss(&C[(i*88)+29], c50_2);
__m128 c50_3 = _mm_load_ss(&C[(i*88)+50]);
__m128 a50_3 = _mm_load_ss(&values[221]);
c50_3 = _mm_add_ss(c50_3, _mm_mul_ss(a50_3, b50));
_mm_store_ss(&C[(i*88)+50], c50_3);
__m128 c50_4 = _mm_load_ss(&C[(i*88)+78]);
__m128 a50_4 = _mm_load_ss(&values[222]);
c50_4 = _mm_add_ss(c50_4, _mm_mul_ss(a50_4, b50));
_mm_store_ss(&C[(i*88)+78], c50_4);
#else
C[(i*88)+4] += values[218] * B[(i*88)+50];
C[(i*88)+14] += values[219] * B[(i*88)+50];
C[(i*88)+29] += values[220] * B[(i*88)+50];
C[(i*88)+50] += values[221] * B[(i*88)+50];
C[(i*88)+78] += values[222] * B[(i*88)+50];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b51 = _mm_broadcast_ss(&B[(i*88)+51]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b51 = _mm_load_ss(&B[(i*88)+51]);
b51 = _mm_shuffle_ps(b51, b51, 0x00);
#endif
__m128 c51_0 = _mm_load_ss(&C[(i*88)+5]);
__m128 a51_0 = _mm_load_ss(&values[223]);
c51_0 = _mm_add_ss(c51_0, _mm_mul_ss(a51_0, b51));
_mm_store_ss(&C[(i*88)+5], c51_0);
__m128 c51_1 = _mm_load_ss(&C[(i*88)+15]);
__m128 a51_1 = _mm_load_ss(&values[224]);
c51_1 = _mm_add_ss(c51_1, _mm_mul_ss(a51_1, b51));
_mm_store_ss(&C[(i*88)+15], c51_1);
__m128 c51_2 = _mm_load_ss(&C[(i*88)+30]);
__m128 a51_2 = _mm_load_ss(&values[225]);
c51_2 = _mm_add_ss(c51_2, _mm_mul_ss(a51_2, b51));
_mm_store_ss(&C[(i*88)+30], c51_2);
__m128 c51_3 = _mm_load_ss(&C[(i*88)+51]);
__m128 a51_3 = _mm_load_ss(&values[226]);
c51_3 = _mm_add_ss(c51_3, _mm_mul_ss(a51_3, b51));
_mm_store_ss(&C[(i*88)+51], c51_3);
__m128 c51_4 = _mm_load_ss(&C[(i*88)+79]);
__m128 a51_4 = _mm_load_ss(&values[227]);
c51_4 = _mm_add_ss(c51_4, _mm_mul_ss(a51_4, b51));
_mm_store_ss(&C[(i*88)+79], c51_4);
#else
C[(i*88)+5] += values[223] * B[(i*88)+51];
C[(i*88)+15] += values[224] * B[(i*88)+51];
C[(i*88)+30] += values[225] * B[(i*88)+51];
C[(i*88)+51] += values[226] * B[(i*88)+51];
C[(i*88)+79] += values[227] * B[(i*88)+51];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b52 = _mm_broadcast_ss(&B[(i*88)+52]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b52 = _mm_load_ss(&B[(i*88)+52]);
b52 = _mm_shuffle_ps(b52, b52, 0x00);
#endif
__m128 c52_0 = _mm_load_ss(&C[(i*88)+6]);
__m128 a52_0 = _mm_load_ss(&values[228]);
c52_0 = _mm_add_ss(c52_0, _mm_mul_ss(a52_0, b52));
_mm_store_ss(&C[(i*88)+6], c52_0);
__m128 c52_1 = _mm_load_ss(&C[(i*88)+16]);
__m128 a52_1 = _mm_load_ss(&values[229]);
c52_1 = _mm_add_ss(c52_1, _mm_mul_ss(a52_1, b52));
_mm_store_ss(&C[(i*88)+16], c52_1);
__m128 c52_2 = _mm_load_ss(&C[(i*88)+31]);
__m128 a52_2 = _mm_load_ss(&values[230]);
c52_2 = _mm_add_ss(c52_2, _mm_mul_ss(a52_2, b52));
_mm_store_ss(&C[(i*88)+31], c52_2);
__m128 c52_3 = _mm_load_ss(&C[(i*88)+52]);
__m128 a52_3 = _mm_load_ss(&values[231]);
c52_3 = _mm_add_ss(c52_3, _mm_mul_ss(a52_3, b52));
_mm_store_ss(&C[(i*88)+52], c52_3);
__m128 c52_4 = _mm_load_ss(&C[(i*88)+80]);
__m128 a52_4 = _mm_load_ss(&values[232]);
c52_4 = _mm_add_ss(c52_4, _mm_mul_ss(a52_4, b52));
_mm_store_ss(&C[(i*88)+80], c52_4);
#else
C[(i*88)+6] += values[228] * B[(i*88)+52];
C[(i*88)+16] += values[229] * B[(i*88)+52];
C[(i*88)+31] += values[230] * B[(i*88)+52];
C[(i*88)+52] += values[231] * B[(i*88)+52];
C[(i*88)+80] += values[232] * B[(i*88)+52];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b53 = _mm_broadcast_ss(&B[(i*88)+53]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b53 = _mm_load_ss(&B[(i*88)+53]);
b53 = _mm_shuffle_ps(b53, b53, 0x00);
#endif
__m128 c53_0 = _mm_load_ss(&C[(i*88)+1]);
__m128 a53_0 = _mm_load_ss(&values[233]);
c53_0 = _mm_add_ss(c53_0, _mm_mul_ss(a53_0, b53));
_mm_store_ss(&C[(i*88)+1], c53_0);
__m128 c53_1 = _mm_load_ss(&C[(i*88)+7]);
__m128 a53_1 = _mm_load_ss(&values[234]);
c53_1 = _mm_add_ss(c53_1, _mm_mul_ss(a53_1, b53));
_mm_store_ss(&C[(i*88)+7], c53_1);
__m128 c53_2 = _mm_load_ss(&C[(i*88)+17]);
__m128 a53_2 = _mm_load_ss(&values[235]);
c53_2 = _mm_add_ss(c53_2, _mm_mul_ss(a53_2, b53));
_mm_store_ss(&C[(i*88)+17], c53_2);
__m128 c53_3 = _mm_load_ss(&C[(i*88)+32]);
__m128 a53_3 = _mm_load_ss(&values[236]);
c53_3 = _mm_add_ss(c53_3, _mm_mul_ss(a53_3, b53));
_mm_store_ss(&C[(i*88)+32], c53_3);
__m128 c53_4 = _mm_load_ss(&C[(i*88)+53]);
__m128 a53_4 = _mm_load_ss(&values[237]);
c53_4 = _mm_add_ss(c53_4, _mm_mul_ss(a53_4, b53));
_mm_store_ss(&C[(i*88)+53], c53_4);
__m128 c53_5 = _mm_load_ss(&C[(i*88)+81]);
__m128 a53_5 = _mm_load_ss(&values[238]);
c53_5 = _mm_add_ss(c53_5, _mm_mul_ss(a53_5, b53));
_mm_store_ss(&C[(i*88)+81], c53_5);
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
__m128 b54 = _mm_broadcast_ss(&B[(i*88)+54]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b54 = _mm_load_ss(&B[(i*88)+54]);
b54 = _mm_shuffle_ps(b54, b54, 0x00);
#endif
__m128 c54_0 = _mm_load_ss(&C[(i*88)+2]);
__m128 a54_0 = _mm_load_ss(&values[239]);
c54_0 = _mm_add_ss(c54_0, _mm_mul_ss(a54_0, b54));
_mm_store_ss(&C[(i*88)+2], c54_0);
__m128 c54_1 = _mm_load_ss(&C[(i*88)+8]);
__m128 a54_1 = _mm_load_ss(&values[240]);
c54_1 = _mm_add_ss(c54_1, _mm_mul_ss(a54_1, b54));
_mm_store_ss(&C[(i*88)+8], c54_1);
__m128 c54_2 = _mm_load_ss(&C[(i*88)+18]);
__m128 a54_2 = _mm_load_ss(&values[241]);
c54_2 = _mm_add_ss(c54_2, _mm_mul_ss(a54_2, b54));
_mm_store_ss(&C[(i*88)+18], c54_2);
__m128 c54_3 = _mm_load_ss(&C[(i*88)+33]);
__m128 a54_3 = _mm_load_ss(&values[242]);
c54_3 = _mm_add_ss(c54_3, _mm_mul_ss(a54_3, b54));
_mm_store_ss(&C[(i*88)+33], c54_3);
__m128 c54_4 = _mm_load_ss(&C[(i*88)+54]);
__m128 a54_4 = _mm_load_ss(&values[243]);
c54_4 = _mm_add_ss(c54_4, _mm_mul_ss(a54_4, b54));
_mm_store_ss(&C[(i*88)+54], c54_4);
__m128 c54_5 = _mm_load_ss(&C[(i*88)+82]);
__m128 a54_5 = _mm_load_ss(&values[244]);
c54_5 = _mm_add_ss(c54_5, _mm_mul_ss(a54_5, b54));
_mm_store_ss(&C[(i*88)+82], c54_5);
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
__m128 b55 = _mm_broadcast_ss(&B[(i*88)+55]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b55 = _mm_load_ss(&B[(i*88)+55]);
b55 = _mm_shuffle_ps(b55, b55, 0x00);
#endif
__m128 c55_0 = _mm_load_ss(&C[(i*88)+0]);
__m128 a55_0 = _mm_load_ss(&values[245]);
c55_0 = _mm_add_ss(c55_0, _mm_mul_ss(a55_0, b55));
_mm_store_ss(&C[(i*88)+0], c55_0);
__m128 c55_1 = _mm_load_ss(&C[(i*88)+3]);
__m128 a55_1 = _mm_load_ss(&values[246]);
c55_1 = _mm_add_ss(c55_1, _mm_mul_ss(a55_1, b55));
_mm_store_ss(&C[(i*88)+3], c55_1);
__m128 c55_2 = _mm_load_ss(&C[(i*88)+9]);
__m128 a55_2 = _mm_load_ss(&values[247]);
c55_2 = _mm_add_ss(c55_2, _mm_mul_ss(a55_2, b55));
_mm_store_ss(&C[(i*88)+9], c55_2);
__m128 c55_3 = _mm_load_ss(&C[(i*88)+19]);
__m128 a55_3 = _mm_load_ss(&values[248]);
c55_3 = _mm_add_ss(c55_3, _mm_mul_ss(a55_3, b55));
_mm_store_ss(&C[(i*88)+19], c55_3);
__m128 c55_4 = _mm_load_ss(&C[(i*88)+34]);
__m128 a55_4 = _mm_load_ss(&values[249]);
c55_4 = _mm_add_ss(c55_4, _mm_mul_ss(a55_4, b55));
_mm_store_ss(&C[(i*88)+34], c55_4);
__m128 c55_5 = _mm_load_ss(&C[(i*88)+55]);
__m128 a55_5 = _mm_load_ss(&values[250]);
c55_5 = _mm_add_ss(c55_5, _mm_mul_ss(a55_5, b55));
_mm_store_ss(&C[(i*88)+55], c55_5);
__m128 c55_6 = _mm_load_ss(&C[(i*88)+83]);
__m128 a55_6 = _mm_load_ss(&values[251]);
c55_6 = _mm_add_ss(c55_6, _mm_mul_ss(a55_6, b55));
_mm_store_ss(&C[(i*88)+83], c55_6);
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
__m128 b56 = _mm_broadcast_ss(&B[(i*88)+56]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b56 = _mm_load_ss(&B[(i*88)+56]);
b56 = _mm_shuffle_ps(b56, b56, 0x00);
#endif
__m128 c56_0 = _mm_load_ss(&C[(i*88)+56]);
__m128 a56_0 = _mm_load_ss(&values[252]);
c56_0 = _mm_add_ss(c56_0, _mm_mul_ss(a56_0, b56));
_mm_store_ss(&C[(i*88)+56], c56_0);
#else
C[(i*88)+56] += values[252] * B[(i*88)+56];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b57 = _mm_broadcast_ss(&B[(i*88)+57]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b57 = _mm_load_ss(&B[(i*88)+57]);
b57 = _mm_shuffle_ps(b57, b57, 0x00);
#endif
__m128 c57_0 = _mm_load_ss(&C[(i*88)+57]);
__m128 a57_0 = _mm_load_ss(&values[253]);
c57_0 = _mm_add_ss(c57_0, _mm_mul_ss(a57_0, b57));
_mm_store_ss(&C[(i*88)+57], c57_0);
#else
C[(i*88)+57] += values[253] * B[(i*88)+57];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b58 = _mm_broadcast_ss(&B[(i*88)+58]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b58 = _mm_load_ss(&B[(i*88)+58]);
b58 = _mm_shuffle_ps(b58, b58, 0x00);
#endif
__m128 c58_0 = _mm_load_ss(&C[(i*88)+58]);
__m128 a58_0 = _mm_load_ss(&values[254]);
c58_0 = _mm_add_ss(c58_0, _mm_mul_ss(a58_0, b58));
_mm_store_ss(&C[(i*88)+58], c58_0);
#else
C[(i*88)+58] += values[254] * B[(i*88)+58];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b59 = _mm_broadcast_ss(&B[(i*88)+59]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b59 = _mm_load_ss(&B[(i*88)+59]);
b59 = _mm_shuffle_ps(b59, b59, 0x00);
#endif
__m128 c59_0 = _mm_load_ss(&C[(i*88)+59]);
__m128 a59_0 = _mm_load_ss(&values[255]);
c59_0 = _mm_add_ss(c59_0, _mm_mul_ss(a59_0, b59));
_mm_store_ss(&C[(i*88)+59], c59_0);
#else
C[(i*88)+59] += values[255] * B[(i*88)+59];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b60 = _mm_broadcast_ss(&B[(i*88)+60]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b60 = _mm_load_ss(&B[(i*88)+60]);
b60 = _mm_shuffle_ps(b60, b60, 0x00);
#endif
__m128 c60_0 = _mm_load_ss(&C[(i*88)+60]);
__m128 a60_0 = _mm_load_ss(&values[256]);
c60_0 = _mm_add_ss(c60_0, _mm_mul_ss(a60_0, b60));
_mm_store_ss(&C[(i*88)+60], c60_0);
#else
C[(i*88)+60] += values[256] * B[(i*88)+60];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b61 = _mm_broadcast_ss(&B[(i*88)+61]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b61 = _mm_load_ss(&B[(i*88)+61]);
b61 = _mm_shuffle_ps(b61, b61, 0x00);
#endif
__m128 c61_0 = _mm_load_ss(&C[(i*88)+61]);
__m128 a61_0 = _mm_load_ss(&values[257]);
c61_0 = _mm_add_ss(c61_0, _mm_mul_ss(a61_0, b61));
_mm_store_ss(&C[(i*88)+61], c61_0);
#else
C[(i*88)+61] += values[257] * B[(i*88)+61];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b62 = _mm_broadcast_ss(&B[(i*88)+62]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b62 = _mm_load_ss(&B[(i*88)+62]);
b62 = _mm_shuffle_ps(b62, b62, 0x00);
#endif
__m128 c62_0 = _mm_load_ss(&C[(i*88)+62]);
__m128 a62_0 = _mm_load_ss(&values[258]);
c62_0 = _mm_add_ss(c62_0, _mm_mul_ss(a62_0, b62));
_mm_store_ss(&C[(i*88)+62], c62_0);
#else
C[(i*88)+62] += values[258] * B[(i*88)+62];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b63 = _mm_broadcast_ss(&B[(i*88)+63]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b63 = _mm_load_ss(&B[(i*88)+63]);
b63 = _mm_shuffle_ps(b63, b63, 0x00);
#endif
__m128 c63_0 = _mm_load_ss(&C[(i*88)+35]);
__m128 a63_0 = _mm_load_ss(&values[259]);
c63_0 = _mm_add_ss(c63_0, _mm_mul_ss(a63_0, b63));
_mm_store_ss(&C[(i*88)+35], c63_0);
__m128 c63_1 = _mm_load_ss(&C[(i*88)+63]);
__m128 a63_1 = _mm_load_ss(&values[260]);
c63_1 = _mm_add_ss(c63_1, _mm_mul_ss(a63_1, b63));
_mm_store_ss(&C[(i*88)+63], c63_1);
#else
C[(i*88)+35] += values[259] * B[(i*88)+63];
C[(i*88)+63] += values[260] * B[(i*88)+63];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b64 = _mm_broadcast_ss(&B[(i*88)+64]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b64 = _mm_load_ss(&B[(i*88)+64]);
b64 = _mm_shuffle_ps(b64, b64, 0x00);
#endif
__m128 c64_0 = _mm_load_ss(&C[(i*88)+36]);
__m128 a64_0 = _mm_load_ss(&values[261]);
c64_0 = _mm_add_ss(c64_0, _mm_mul_ss(a64_0, b64));
_mm_store_ss(&C[(i*88)+36], c64_0);
__m128 c64_1 = _mm_load_ss(&C[(i*88)+64]);
__m128 a64_1 = _mm_load_ss(&values[262]);
c64_1 = _mm_add_ss(c64_1, _mm_mul_ss(a64_1, b64));
_mm_store_ss(&C[(i*88)+64], c64_1);
#else
C[(i*88)+36] += values[261] * B[(i*88)+64];
C[(i*88)+64] += values[262] * B[(i*88)+64];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b65 = _mm_broadcast_ss(&B[(i*88)+65]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b65 = _mm_load_ss(&B[(i*88)+65]);
b65 = _mm_shuffle_ps(b65, b65, 0x00);
#endif
__m128 c65_0 = _mm_load_ss(&C[(i*88)+37]);
__m128 a65_0 = _mm_load_ss(&values[263]);
c65_0 = _mm_add_ss(c65_0, _mm_mul_ss(a65_0, b65));
_mm_store_ss(&C[(i*88)+37], c65_0);
__m128 c65_1 = _mm_load_ss(&C[(i*88)+65]);
__m128 a65_1 = _mm_load_ss(&values[264]);
c65_1 = _mm_add_ss(c65_1, _mm_mul_ss(a65_1, b65));
_mm_store_ss(&C[(i*88)+65], c65_1);
#else
C[(i*88)+37] += values[263] * B[(i*88)+65];
C[(i*88)+65] += values[264] * B[(i*88)+65];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b66 = _mm_broadcast_ss(&B[(i*88)+66]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b66 = _mm_load_ss(&B[(i*88)+66]);
b66 = _mm_shuffle_ps(b66, b66, 0x00);
#endif
__m128 c66_0 = _mm_load_ss(&C[(i*88)+38]);
__m128 a66_0 = _mm_load_ss(&values[265]);
c66_0 = _mm_add_ss(c66_0, _mm_mul_ss(a66_0, b66));
_mm_store_ss(&C[(i*88)+38], c66_0);
__m128 c66_1 = _mm_load_ss(&C[(i*88)+66]);
__m128 a66_1 = _mm_load_ss(&values[266]);
c66_1 = _mm_add_ss(c66_1, _mm_mul_ss(a66_1, b66));
_mm_store_ss(&C[(i*88)+66], c66_1);
#else
C[(i*88)+38] += values[265] * B[(i*88)+66];
C[(i*88)+66] += values[266] * B[(i*88)+66];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b67 = _mm_broadcast_ss(&B[(i*88)+67]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b67 = _mm_load_ss(&B[(i*88)+67]);
b67 = _mm_shuffle_ps(b67, b67, 0x00);
#endif
__m128 c67_0 = _mm_load_ss(&C[(i*88)+39]);
__m128 a67_0 = _mm_load_ss(&values[267]);
c67_0 = _mm_add_ss(c67_0, _mm_mul_ss(a67_0, b67));
_mm_store_ss(&C[(i*88)+39], c67_0);
__m128 c67_1 = _mm_load_ss(&C[(i*88)+67]);
__m128 a67_1 = _mm_load_ss(&values[268]);
c67_1 = _mm_add_ss(c67_1, _mm_mul_ss(a67_1, b67));
_mm_store_ss(&C[(i*88)+67], c67_1);
#else
C[(i*88)+39] += values[267] * B[(i*88)+67];
C[(i*88)+67] += values[268] * B[(i*88)+67];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b68 = _mm_broadcast_ss(&B[(i*88)+68]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b68 = _mm_load_ss(&B[(i*88)+68]);
b68 = _mm_shuffle_ps(b68, b68, 0x00);
#endif
__m128 c68_0 = _mm_load_ss(&C[(i*88)+40]);
__m128 a68_0 = _mm_load_ss(&values[269]);
c68_0 = _mm_add_ss(c68_0, _mm_mul_ss(a68_0, b68));
_mm_store_ss(&C[(i*88)+40], c68_0);
__m128 c68_1 = _mm_load_ss(&C[(i*88)+68]);
__m128 a68_1 = _mm_load_ss(&values[270]);
c68_1 = _mm_add_ss(c68_1, _mm_mul_ss(a68_1, b68));
_mm_store_ss(&C[(i*88)+68], c68_1);
#else
C[(i*88)+40] += values[269] * B[(i*88)+68];
C[(i*88)+68] += values[270] * B[(i*88)+68];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b69 = _mm_broadcast_ss(&B[(i*88)+69]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b69 = _mm_load_ss(&B[(i*88)+69]);
b69 = _mm_shuffle_ps(b69, b69, 0x00);
#endif
__m128 c69_0 = _mm_load_ss(&C[(i*88)+20]);
__m128 a69_0 = _mm_load_ss(&values[271]);
c69_0 = _mm_add_ss(c69_0, _mm_mul_ss(a69_0, b69));
_mm_store_ss(&C[(i*88)+20], c69_0);
__m128 c69_1 = _mm_load_ss(&C[(i*88)+41]);
__m128 a69_1 = _mm_load_ss(&values[272]);
c69_1 = _mm_add_ss(c69_1, _mm_mul_ss(a69_1, b69));
_mm_store_ss(&C[(i*88)+41], c69_1);
__m128 c69_2 = _mm_load_ss(&C[(i*88)+69]);
__m128 a69_2 = _mm_load_ss(&values[273]);
c69_2 = _mm_add_ss(c69_2, _mm_mul_ss(a69_2, b69));
_mm_store_ss(&C[(i*88)+69], c69_2);
#else
C[(i*88)+20] += values[271] * B[(i*88)+69];
C[(i*88)+41] += values[272] * B[(i*88)+69];
C[(i*88)+69] += values[273] * B[(i*88)+69];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b70 = _mm_broadcast_ss(&B[(i*88)+70]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b70 = _mm_load_ss(&B[(i*88)+70]);
b70 = _mm_shuffle_ps(b70, b70, 0x00);
#endif
__m128 c70_0 = _mm_load_ss(&C[(i*88)+21]);
__m128 a70_0 = _mm_load_ss(&values[274]);
c70_0 = _mm_add_ss(c70_0, _mm_mul_ss(a70_0, b70));
_mm_store_ss(&C[(i*88)+21], c70_0);
__m128 c70_1 = _mm_load_ss(&C[(i*88)+42]);
__m128 a70_1 = _mm_load_ss(&values[275]);
c70_1 = _mm_add_ss(c70_1, _mm_mul_ss(a70_1, b70));
_mm_store_ss(&C[(i*88)+42], c70_1);
__m128 c70_2 = _mm_load_ss(&C[(i*88)+70]);
__m128 a70_2 = _mm_load_ss(&values[276]);
c70_2 = _mm_add_ss(c70_2, _mm_mul_ss(a70_2, b70));
_mm_store_ss(&C[(i*88)+70], c70_2);
#else
C[(i*88)+21] += values[274] * B[(i*88)+70];
C[(i*88)+42] += values[275] * B[(i*88)+70];
C[(i*88)+70] += values[276] * B[(i*88)+70];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b71 = _mm_broadcast_ss(&B[(i*88)+71]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b71 = _mm_load_ss(&B[(i*88)+71]);
b71 = _mm_shuffle_ps(b71, b71, 0x00);
#endif
__m128 c71_0 = _mm_load_ss(&C[(i*88)+22]);
__m128 a71_0 = _mm_load_ss(&values[277]);
c71_0 = _mm_add_ss(c71_0, _mm_mul_ss(a71_0, b71));
_mm_store_ss(&C[(i*88)+22], c71_0);
__m128 c71_1 = _mm_load_ss(&C[(i*88)+43]);
__m128 a71_1 = _mm_load_ss(&values[278]);
c71_1 = _mm_add_ss(c71_1, _mm_mul_ss(a71_1, b71));
_mm_store_ss(&C[(i*88)+43], c71_1);
__m128 c71_2 = _mm_load_ss(&C[(i*88)+71]);
__m128 a71_2 = _mm_load_ss(&values[279]);
c71_2 = _mm_add_ss(c71_2, _mm_mul_ss(a71_2, b71));
_mm_store_ss(&C[(i*88)+71], c71_2);
#else
C[(i*88)+22] += values[277] * B[(i*88)+71];
C[(i*88)+43] += values[278] * B[(i*88)+71];
C[(i*88)+71] += values[279] * B[(i*88)+71];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b72 = _mm_broadcast_ss(&B[(i*88)+72]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b72 = _mm_load_ss(&B[(i*88)+72]);
b72 = _mm_shuffle_ps(b72, b72, 0x00);
#endif
__m128 c72_0 = _mm_load_ss(&C[(i*88)+23]);
__m128 a72_0 = _mm_load_ss(&values[280]);
c72_0 = _mm_add_ss(c72_0, _mm_mul_ss(a72_0, b72));
_mm_store_ss(&C[(i*88)+23], c72_0);
__m128 c72_1 = _mm_load_ss(&C[(i*88)+44]);
__m128 a72_1 = _mm_load_ss(&values[281]);
c72_1 = _mm_add_ss(c72_1, _mm_mul_ss(a72_1, b72));
_mm_store_ss(&C[(i*88)+44], c72_1);
__m128 c72_2 = _mm_load_ss(&C[(i*88)+72]);
__m128 a72_2 = _mm_load_ss(&values[282]);
c72_2 = _mm_add_ss(c72_2, _mm_mul_ss(a72_2, b72));
_mm_store_ss(&C[(i*88)+72], c72_2);
#else
C[(i*88)+23] += values[280] * B[(i*88)+72];
C[(i*88)+44] += values[281] * B[(i*88)+72];
C[(i*88)+72] += values[282] * B[(i*88)+72];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b73 = _mm_broadcast_ss(&B[(i*88)+73]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b73 = _mm_load_ss(&B[(i*88)+73]);
b73 = _mm_shuffle_ps(b73, b73, 0x00);
#endif
__m128 c73_0 = _mm_load_ss(&C[(i*88)+24]);
__m128 a73_0 = _mm_load_ss(&values[283]);
c73_0 = _mm_add_ss(c73_0, _mm_mul_ss(a73_0, b73));
_mm_store_ss(&C[(i*88)+24], c73_0);
__m128 c73_1 = _mm_load_ss(&C[(i*88)+45]);
__m128 a73_1 = _mm_load_ss(&values[284]);
c73_1 = _mm_add_ss(c73_1, _mm_mul_ss(a73_1, b73));
_mm_store_ss(&C[(i*88)+45], c73_1);
__m128 c73_2 = _mm_load_ss(&C[(i*88)+73]);
__m128 a73_2 = _mm_load_ss(&values[285]);
c73_2 = _mm_add_ss(c73_2, _mm_mul_ss(a73_2, b73));
_mm_store_ss(&C[(i*88)+73], c73_2);
#else
C[(i*88)+24] += values[283] * B[(i*88)+73];
C[(i*88)+45] += values[284] * B[(i*88)+73];
C[(i*88)+73] += values[285] * B[(i*88)+73];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b74 = _mm_broadcast_ss(&B[(i*88)+74]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b74 = _mm_load_ss(&B[(i*88)+74]);
b74 = _mm_shuffle_ps(b74, b74, 0x00);
#endif
__m128 c74_0 = _mm_load_ss(&C[(i*88)+10]);
__m128 a74_0 = _mm_load_ss(&values[286]);
c74_0 = _mm_add_ss(c74_0, _mm_mul_ss(a74_0, b74));
_mm_store_ss(&C[(i*88)+10], c74_0);
__m128 c74_1 = _mm_load_ss(&C[(i*88)+25]);
__m128 a74_1 = _mm_load_ss(&values[287]);
c74_1 = _mm_add_ss(c74_1, _mm_mul_ss(a74_1, b74));
_mm_store_ss(&C[(i*88)+25], c74_1);
__m128 c74_2 = _mm_load_ss(&C[(i*88)+46]);
__m128 a74_2 = _mm_load_ss(&values[288]);
c74_2 = _mm_add_ss(c74_2, _mm_mul_ss(a74_2, b74));
_mm_store_ss(&C[(i*88)+46], c74_2);
__m128 c74_3 = _mm_load_ss(&C[(i*88)+74]);
__m128 a74_3 = _mm_load_ss(&values[289]);
c74_3 = _mm_add_ss(c74_3, _mm_mul_ss(a74_3, b74));
_mm_store_ss(&C[(i*88)+74], c74_3);
#else
C[(i*88)+10] += values[286] * B[(i*88)+74];
C[(i*88)+25] += values[287] * B[(i*88)+74];
C[(i*88)+46] += values[288] * B[(i*88)+74];
C[(i*88)+74] += values[289] * B[(i*88)+74];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b75 = _mm_broadcast_ss(&B[(i*88)+75]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b75 = _mm_load_ss(&B[(i*88)+75]);
b75 = _mm_shuffle_ps(b75, b75, 0x00);
#endif
__m128 c75_0 = _mm_load_ss(&C[(i*88)+11]);
__m128 a75_0 = _mm_load_ss(&values[290]);
c75_0 = _mm_add_ss(c75_0, _mm_mul_ss(a75_0, b75));
_mm_store_ss(&C[(i*88)+11], c75_0);
__m128 c75_1 = _mm_load_ss(&C[(i*88)+26]);
__m128 a75_1 = _mm_load_ss(&values[291]);
c75_1 = _mm_add_ss(c75_1, _mm_mul_ss(a75_1, b75));
_mm_store_ss(&C[(i*88)+26], c75_1);
__m128 c75_2 = _mm_load_ss(&C[(i*88)+47]);
__m128 a75_2 = _mm_load_ss(&values[292]);
c75_2 = _mm_add_ss(c75_2, _mm_mul_ss(a75_2, b75));
_mm_store_ss(&C[(i*88)+47], c75_2);
__m128 c75_3 = _mm_load_ss(&C[(i*88)+75]);
__m128 a75_3 = _mm_load_ss(&values[293]);
c75_3 = _mm_add_ss(c75_3, _mm_mul_ss(a75_3, b75));
_mm_store_ss(&C[(i*88)+75], c75_3);
#else
C[(i*88)+11] += values[290] * B[(i*88)+75];
C[(i*88)+26] += values[291] * B[(i*88)+75];
C[(i*88)+47] += values[292] * B[(i*88)+75];
C[(i*88)+75] += values[293] * B[(i*88)+75];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b76 = _mm_broadcast_ss(&B[(i*88)+76]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b76 = _mm_load_ss(&B[(i*88)+76]);
b76 = _mm_shuffle_ps(b76, b76, 0x00);
#endif
__m128 c76_0 = _mm_load_ss(&C[(i*88)+12]);
__m128 a76_0 = _mm_load_ss(&values[294]);
c76_0 = _mm_add_ss(c76_0, _mm_mul_ss(a76_0, b76));
_mm_store_ss(&C[(i*88)+12], c76_0);
__m128 c76_1 = _mm_load_ss(&C[(i*88)+27]);
__m128 a76_1 = _mm_load_ss(&values[295]);
c76_1 = _mm_add_ss(c76_1, _mm_mul_ss(a76_1, b76));
_mm_store_ss(&C[(i*88)+27], c76_1);
__m128 c76_2 = _mm_load_ss(&C[(i*88)+48]);
__m128 a76_2 = _mm_load_ss(&values[296]);
c76_2 = _mm_add_ss(c76_2, _mm_mul_ss(a76_2, b76));
_mm_store_ss(&C[(i*88)+48], c76_2);
__m128 c76_3 = _mm_load_ss(&C[(i*88)+76]);
__m128 a76_3 = _mm_load_ss(&values[297]);
c76_3 = _mm_add_ss(c76_3, _mm_mul_ss(a76_3, b76));
_mm_store_ss(&C[(i*88)+76], c76_3);
#else
C[(i*88)+12] += values[294] * B[(i*88)+76];
C[(i*88)+27] += values[295] * B[(i*88)+76];
C[(i*88)+48] += values[296] * B[(i*88)+76];
C[(i*88)+76] += values[297] * B[(i*88)+76];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b77 = _mm_broadcast_ss(&B[(i*88)+77]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b77 = _mm_load_ss(&B[(i*88)+77]);
b77 = _mm_shuffle_ps(b77, b77, 0x00);
#endif
__m128 c77_0 = _mm_load_ss(&C[(i*88)+13]);
__m128 a77_0 = _mm_load_ss(&values[298]);
c77_0 = _mm_add_ss(c77_0, _mm_mul_ss(a77_0, b77));
_mm_store_ss(&C[(i*88)+13], c77_0);
__m128 c77_1 = _mm_load_ss(&C[(i*88)+28]);
__m128 a77_1 = _mm_load_ss(&values[299]);
c77_1 = _mm_add_ss(c77_1, _mm_mul_ss(a77_1, b77));
_mm_store_ss(&C[(i*88)+28], c77_1);
__m128 c77_2 = _mm_load_ss(&C[(i*88)+49]);
__m128 a77_2 = _mm_load_ss(&values[300]);
c77_2 = _mm_add_ss(c77_2, _mm_mul_ss(a77_2, b77));
_mm_store_ss(&C[(i*88)+49], c77_2);
__m128 c77_3 = _mm_load_ss(&C[(i*88)+77]);
__m128 a77_3 = _mm_load_ss(&values[301]);
c77_3 = _mm_add_ss(c77_3, _mm_mul_ss(a77_3, b77));
_mm_store_ss(&C[(i*88)+77], c77_3);
#else
C[(i*88)+13] += values[298] * B[(i*88)+77];
C[(i*88)+28] += values[299] * B[(i*88)+77];
C[(i*88)+49] += values[300] * B[(i*88)+77];
C[(i*88)+77] += values[301] * B[(i*88)+77];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b78 = _mm_broadcast_ss(&B[(i*88)+78]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b78 = _mm_load_ss(&B[(i*88)+78]);
b78 = _mm_shuffle_ps(b78, b78, 0x00);
#endif
__m128 c78_0 = _mm_load_ss(&C[(i*88)+4]);
__m128 a78_0 = _mm_load_ss(&values[302]);
c78_0 = _mm_add_ss(c78_0, _mm_mul_ss(a78_0, b78));
_mm_store_ss(&C[(i*88)+4], c78_0);
__m128 c78_1 = _mm_load_ss(&C[(i*88)+14]);
__m128 a78_1 = _mm_load_ss(&values[303]);
c78_1 = _mm_add_ss(c78_1, _mm_mul_ss(a78_1, b78));
_mm_store_ss(&C[(i*88)+14], c78_1);
__m128 c78_2 = _mm_load_ss(&C[(i*88)+29]);
__m128 a78_2 = _mm_load_ss(&values[304]);
c78_2 = _mm_add_ss(c78_2, _mm_mul_ss(a78_2, b78));
_mm_store_ss(&C[(i*88)+29], c78_2);
__m128 c78_3 = _mm_load_ss(&C[(i*88)+50]);
__m128 a78_3 = _mm_load_ss(&values[305]);
c78_3 = _mm_add_ss(c78_3, _mm_mul_ss(a78_3, b78));
_mm_store_ss(&C[(i*88)+50], c78_3);
__m128 c78_4 = _mm_load_ss(&C[(i*88)+78]);
__m128 a78_4 = _mm_load_ss(&values[306]);
c78_4 = _mm_add_ss(c78_4, _mm_mul_ss(a78_4, b78));
_mm_store_ss(&C[(i*88)+78], c78_4);
#else
C[(i*88)+4] += values[302] * B[(i*88)+78];
C[(i*88)+14] += values[303] * B[(i*88)+78];
C[(i*88)+29] += values[304] * B[(i*88)+78];
C[(i*88)+50] += values[305] * B[(i*88)+78];
C[(i*88)+78] += values[306] * B[(i*88)+78];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b79 = _mm_broadcast_ss(&B[(i*88)+79]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b79 = _mm_load_ss(&B[(i*88)+79]);
b79 = _mm_shuffle_ps(b79, b79, 0x00);
#endif
__m128 c79_0 = _mm_load_ss(&C[(i*88)+5]);
__m128 a79_0 = _mm_load_ss(&values[307]);
c79_0 = _mm_add_ss(c79_0, _mm_mul_ss(a79_0, b79));
_mm_store_ss(&C[(i*88)+5], c79_0);
__m128 c79_1 = _mm_load_ss(&C[(i*88)+15]);
__m128 a79_1 = _mm_load_ss(&values[308]);
c79_1 = _mm_add_ss(c79_1, _mm_mul_ss(a79_1, b79));
_mm_store_ss(&C[(i*88)+15], c79_1);
__m128 c79_2 = _mm_load_ss(&C[(i*88)+30]);
__m128 a79_2 = _mm_load_ss(&values[309]);
c79_2 = _mm_add_ss(c79_2, _mm_mul_ss(a79_2, b79));
_mm_store_ss(&C[(i*88)+30], c79_2);
__m128 c79_3 = _mm_load_ss(&C[(i*88)+51]);
__m128 a79_3 = _mm_load_ss(&values[310]);
c79_3 = _mm_add_ss(c79_3, _mm_mul_ss(a79_3, b79));
_mm_store_ss(&C[(i*88)+51], c79_3);
__m128 c79_4 = _mm_load_ss(&C[(i*88)+79]);
__m128 a79_4 = _mm_load_ss(&values[311]);
c79_4 = _mm_add_ss(c79_4, _mm_mul_ss(a79_4, b79));
_mm_store_ss(&C[(i*88)+79], c79_4);
#else
C[(i*88)+5] += values[307] * B[(i*88)+79];
C[(i*88)+15] += values[308] * B[(i*88)+79];
C[(i*88)+30] += values[309] * B[(i*88)+79];
C[(i*88)+51] += values[310] * B[(i*88)+79];
C[(i*88)+79] += values[311] * B[(i*88)+79];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b80 = _mm_broadcast_ss(&B[(i*88)+80]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b80 = _mm_load_ss(&B[(i*88)+80]);
b80 = _mm_shuffle_ps(b80, b80, 0x00);
#endif
__m128 c80_0 = _mm_load_ss(&C[(i*88)+6]);
__m128 a80_0 = _mm_load_ss(&values[312]);
c80_0 = _mm_add_ss(c80_0, _mm_mul_ss(a80_0, b80));
_mm_store_ss(&C[(i*88)+6], c80_0);
__m128 c80_1 = _mm_load_ss(&C[(i*88)+16]);
__m128 a80_1 = _mm_load_ss(&values[313]);
c80_1 = _mm_add_ss(c80_1, _mm_mul_ss(a80_1, b80));
_mm_store_ss(&C[(i*88)+16], c80_1);
__m128 c80_2 = _mm_load_ss(&C[(i*88)+31]);
__m128 a80_2 = _mm_load_ss(&values[314]);
c80_2 = _mm_add_ss(c80_2, _mm_mul_ss(a80_2, b80));
_mm_store_ss(&C[(i*88)+31], c80_2);
__m128 c80_3 = _mm_load_ss(&C[(i*88)+52]);
__m128 a80_3 = _mm_load_ss(&values[315]);
c80_3 = _mm_add_ss(c80_3, _mm_mul_ss(a80_3, b80));
_mm_store_ss(&C[(i*88)+52], c80_3);
__m128 c80_4 = _mm_load_ss(&C[(i*88)+80]);
__m128 a80_4 = _mm_load_ss(&values[316]);
c80_4 = _mm_add_ss(c80_4, _mm_mul_ss(a80_4, b80));
_mm_store_ss(&C[(i*88)+80], c80_4);
#else
C[(i*88)+6] += values[312] * B[(i*88)+80];
C[(i*88)+16] += values[313] * B[(i*88)+80];
C[(i*88)+31] += values[314] * B[(i*88)+80];
C[(i*88)+52] += values[315] * B[(i*88)+80];
C[(i*88)+80] += values[316] * B[(i*88)+80];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m128 b81 = _mm_broadcast_ss(&B[(i*88)+81]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b81 = _mm_load_ss(&B[(i*88)+81]);
b81 = _mm_shuffle_ps(b81, b81, 0x00);
#endif
__m128 c81_0 = _mm_load_ss(&C[(i*88)+1]);
__m128 a81_0 = _mm_load_ss(&values[317]);
c81_0 = _mm_add_ss(c81_0, _mm_mul_ss(a81_0, b81));
_mm_store_ss(&C[(i*88)+1], c81_0);
__m128 c81_1 = _mm_load_ss(&C[(i*88)+7]);
__m128 a81_1 = _mm_load_ss(&values[318]);
c81_1 = _mm_add_ss(c81_1, _mm_mul_ss(a81_1, b81));
_mm_store_ss(&C[(i*88)+7], c81_1);
__m128 c81_2 = _mm_load_ss(&C[(i*88)+17]);
__m128 a81_2 = _mm_load_ss(&values[319]);
c81_2 = _mm_add_ss(c81_2, _mm_mul_ss(a81_2, b81));
_mm_store_ss(&C[(i*88)+17], c81_2);
__m128 c81_3 = _mm_load_ss(&C[(i*88)+32]);
__m128 a81_3 = _mm_load_ss(&values[320]);
c81_3 = _mm_add_ss(c81_3, _mm_mul_ss(a81_3, b81));
_mm_store_ss(&C[(i*88)+32], c81_3);
__m128 c81_4 = _mm_load_ss(&C[(i*88)+53]);
__m128 a81_4 = _mm_load_ss(&values[321]);
c81_4 = _mm_add_ss(c81_4, _mm_mul_ss(a81_4, b81));
_mm_store_ss(&C[(i*88)+53], c81_4);
__m128 c81_5 = _mm_load_ss(&C[(i*88)+81]);
__m128 a81_5 = _mm_load_ss(&values[322]);
c81_5 = _mm_add_ss(c81_5, _mm_mul_ss(a81_5, b81));
_mm_store_ss(&C[(i*88)+81], c81_5);
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
__m128 b82 = _mm_broadcast_ss(&B[(i*88)+82]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b82 = _mm_load_ss(&B[(i*88)+82]);
b82 = _mm_shuffle_ps(b82, b82, 0x00);
#endif
__m128 c82_0 = _mm_load_ss(&C[(i*88)+2]);
__m128 a82_0 = _mm_load_ss(&values[323]);
c82_0 = _mm_add_ss(c82_0, _mm_mul_ss(a82_0, b82));
_mm_store_ss(&C[(i*88)+2], c82_0);
__m128 c82_1 = _mm_load_ss(&C[(i*88)+8]);
__m128 a82_1 = _mm_load_ss(&values[324]);
c82_1 = _mm_add_ss(c82_1, _mm_mul_ss(a82_1, b82));
_mm_store_ss(&C[(i*88)+8], c82_1);
__m128 c82_2 = _mm_load_ss(&C[(i*88)+18]);
__m128 a82_2 = _mm_load_ss(&values[325]);
c82_2 = _mm_add_ss(c82_2, _mm_mul_ss(a82_2, b82));
_mm_store_ss(&C[(i*88)+18], c82_2);
__m128 c82_3 = _mm_load_ss(&C[(i*88)+33]);
__m128 a82_3 = _mm_load_ss(&values[326]);
c82_3 = _mm_add_ss(c82_3, _mm_mul_ss(a82_3, b82));
_mm_store_ss(&C[(i*88)+33], c82_3);
__m128 c82_4 = _mm_load_ss(&C[(i*88)+54]);
__m128 a82_4 = _mm_load_ss(&values[327]);
c82_4 = _mm_add_ss(c82_4, _mm_mul_ss(a82_4, b82));
_mm_store_ss(&C[(i*88)+54], c82_4);
__m128 c82_5 = _mm_load_ss(&C[(i*88)+82]);
__m128 a82_5 = _mm_load_ss(&values[328]);
c82_5 = _mm_add_ss(c82_5, _mm_mul_ss(a82_5, b82));
_mm_store_ss(&C[(i*88)+82], c82_5);
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
__m128 b83 = _mm_broadcast_ss(&B[(i*88)+83]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128 b83 = _mm_load_ss(&B[(i*88)+83]);
b83 = _mm_shuffle_ps(b83, b83, 0x00);
#endif
__m128 c83_0 = _mm_load_ss(&C[(i*88)+0]);
__m128 a83_0 = _mm_load_ss(&values[329]);
c83_0 = _mm_add_ss(c83_0, _mm_mul_ss(a83_0, b83));
_mm_store_ss(&C[(i*88)+0], c83_0);
__m128 c83_1 = _mm_load_ss(&C[(i*88)+3]);
__m128 a83_1 = _mm_load_ss(&values[330]);
c83_1 = _mm_add_ss(c83_1, _mm_mul_ss(a83_1, b83));
_mm_store_ss(&C[(i*88)+3], c83_1);
__m128 c83_2 = _mm_load_ss(&C[(i*88)+9]);
__m128 a83_2 = _mm_load_ss(&values[331]);
c83_2 = _mm_add_ss(c83_2, _mm_mul_ss(a83_2, b83));
_mm_store_ss(&C[(i*88)+9], c83_2);
__m128 c83_3 = _mm_load_ss(&C[(i*88)+19]);
__m128 a83_3 = _mm_load_ss(&values[332]);
c83_3 = _mm_add_ss(c83_3, _mm_mul_ss(a83_3, b83));
_mm_store_ss(&C[(i*88)+19], c83_3);
__m128 c83_4 = _mm_load_ss(&C[(i*88)+34]);
__m128 a83_4 = _mm_load_ss(&values[333]);
c83_4 = _mm_add_ss(c83_4, _mm_mul_ss(a83_4, b83));
_mm_store_ss(&C[(i*88)+34], c83_4);
__m128 c83_5 = _mm_load_ss(&C[(i*88)+55]);
__m128 a83_5 = _mm_load_ss(&values[334]);
c83_5 = _mm_add_ss(c83_5, _mm_mul_ss(a83_5, b83));
_mm_store_ss(&C[(i*88)+55], c83_5);
__m128 c83_6 = _mm_load_ss(&C[(i*88)+83]);
__m128 a83_6 = _mm_load_ss(&values[335]);
c83_6 = _mm_add_ss(c83_6, _mm_mul_ss(a83_6, b83));
_mm_store_ss(&C[(i*88)+83], c83_6);
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

void ssparse_starMatrix_m120_n9_k9_ldA120_ldBna8_ldC120_beta1_pfsigonly(const float* A, const float* values, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
#pragma vector aligned
for (unsigned int i = 0; i < 120; i += 1)
{
C[(i)+(0)] += A[(i)+(720)] * values[0];
C[(i)+(0)] += A[(i)+(840)] * values[1];
C[(i)+(0)] += A[(i)+(960)] * values[2];
C[(i)+(120)] += A[(i)+(720)] * values[3];
C[(i)+(120)] += A[(i)+(840)] * values[4];
C[(i)+(120)] += A[(i)+(960)] * values[5];
C[(i)+(240)] += A[(i)+(720)] * values[6];
C[(i)+(240)] += A[(i)+(840)] * values[7];
C[(i)+(240)] += A[(i)+(960)] * values[8];
C[(i)+(360)] += A[(i)+(720)] * values[9];
C[(i)+(360)] += A[(i)+(840)] * values[10];
C[(i)+(480)] += A[(i)+(840)] * values[11];
C[(i)+(480)] += A[(i)+(960)] * values[12];
C[(i)+(600)] += A[(i)+(720)] * values[13];
C[(i)+(600)] += A[(i)+(960)] * values[14];
C[(i)+(720)] += A[(i)+(0)] * values[15];
C[(i)+(720)] += A[(i)+(360)] * values[16];
C[(i)+(720)] += A[(i)+(600)] * values[17];
C[(i)+(840)] += A[(i)+(120)] * values[18];
C[(i)+(840)] += A[(i)+(360)] * values[19];
C[(i)+(840)] += A[(i)+(480)] * values[20];
C[(i)+(960)] += A[(i)+(240)] * values[21];
C[(i)+(960)] += A[(i)+(480)] * values[22];
C[(i)+(960)] += A[(i)+(600)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 5760;
#endif

}

#endif
