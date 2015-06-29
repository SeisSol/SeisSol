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
// @date 2015-05-09 22:18:08.041117
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

void dsparse_starMatrix_m1_n9_k9_ldA4_ldBna2_ldC4_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
for (unsigned int i = 0; i < 1; i += 1)
{
C[(i)+(0)] += A[(i)+(24)] * values[0];
C[(i)+(0)] += A[(i)+(28)] * values[1];
C[(i)+(0)] += A[(i)+(32)] * values[2];
C[(i)+(4)] += A[(i)+(24)] * values[3];
C[(i)+(4)] += A[(i)+(28)] * values[4];
C[(i)+(4)] += A[(i)+(32)] * values[5];
C[(i)+(8)] += A[(i)+(24)] * values[6];
C[(i)+(8)] += A[(i)+(28)] * values[7];
C[(i)+(8)] += A[(i)+(32)] * values[8];
C[(i)+(12)] += A[(i)+(24)] * values[9];
C[(i)+(12)] += A[(i)+(28)] * values[10];
C[(i)+(16)] += A[(i)+(28)] * values[11];
C[(i)+(16)] += A[(i)+(32)] * values[12];
C[(i)+(20)] += A[(i)+(24)] * values[13];
C[(i)+(20)] += A[(i)+(32)] * values[14];
C[(i)+(24)] += A[(i)+(0)] * values[15];
C[(i)+(24)] += A[(i)+(12)] * values[16];
C[(i)+(24)] += A[(i)+(20)] * values[17];
C[(i)+(28)] += A[(i)+(4)] * values[18];
C[(i)+(28)] += A[(i)+(12)] * values[19];
C[(i)+(28)] += A[(i)+(16)] * values[20];
C[(i)+(32)] += A[(i)+(8)] * values[21];
C[(i)+(32)] += A[(i)+(16)] * values[22];
C[(i)+(32)] += A[(i)+(20)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 48;
#endif

}

void dsparse_starMatrix_m4_n9_k9_ldA4_ldBna3_ldC4_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(4)
#pragma vector aligned
for (unsigned int i = 0; i < 4; i += 1)
{
C[(i)+(0)] += A[(i)+(24)] * values[0];
C[(i)+(0)] += A[(i)+(28)] * values[1];
C[(i)+(0)] += A[(i)+(32)] * values[2];
C[(i)+(4)] += A[(i)+(24)] * values[3];
C[(i)+(4)] += A[(i)+(28)] * values[4];
C[(i)+(4)] += A[(i)+(32)] * values[5];
C[(i)+(8)] += A[(i)+(24)] * values[6];
C[(i)+(8)] += A[(i)+(28)] * values[7];
C[(i)+(8)] += A[(i)+(32)] * values[8];
C[(i)+(12)] += A[(i)+(24)] * values[9];
C[(i)+(12)] += A[(i)+(28)] * values[10];
C[(i)+(16)] += A[(i)+(28)] * values[11];
C[(i)+(16)] += A[(i)+(32)] * values[12];
C[(i)+(20)] += A[(i)+(24)] * values[13];
C[(i)+(20)] += A[(i)+(32)] * values[14];
C[(i)+(24)] += A[(i)+(0)] * values[15];
C[(i)+(24)] += A[(i)+(12)] * values[16];
C[(i)+(24)] += A[(i)+(20)] * values[17];
C[(i)+(28)] += A[(i)+(4)] * values[18];
C[(i)+(28)] += A[(i)+(12)] * values[19];
C[(i)+(28)] += A[(i)+(16)] * values[20];
C[(i)+(32)] += A[(i)+(8)] * values[21];
C[(i)+(32)] += A[(i)+(16)] * values[22];
C[(i)+(32)] += A[(i)+(20)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 192;
#endif

}

void dsparse_starMatrix_m1_n9_k9_ldA4_ldBna3_ldC4_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
for (unsigned int i = 0; i < 1; i += 1)
{
C[(i)+(0)] += A[(i)+(24)] * values[0];
C[(i)+(0)] += A[(i)+(28)] * values[1];
C[(i)+(0)] += A[(i)+(32)] * values[2];
C[(i)+(4)] += A[(i)+(24)] * values[3];
C[(i)+(4)] += A[(i)+(28)] * values[4];
C[(i)+(4)] += A[(i)+(32)] * values[5];
C[(i)+(8)] += A[(i)+(24)] * values[6];
C[(i)+(8)] += A[(i)+(28)] * values[7];
C[(i)+(8)] += A[(i)+(32)] * values[8];
C[(i)+(12)] += A[(i)+(24)] * values[9];
C[(i)+(12)] += A[(i)+(28)] * values[10];
C[(i)+(16)] += A[(i)+(28)] * values[11];
C[(i)+(16)] += A[(i)+(32)] * values[12];
C[(i)+(20)] += A[(i)+(24)] * values[13];
C[(i)+(20)] += A[(i)+(32)] * values[14];
C[(i)+(24)] += A[(i)+(0)] * values[15];
C[(i)+(24)] += A[(i)+(12)] * values[16];
C[(i)+(24)] += A[(i)+(20)] * values[17];
C[(i)+(28)] += A[(i)+(4)] * values[18];
C[(i)+(28)] += A[(i)+(12)] * values[19];
C[(i)+(28)] += A[(i)+(16)] * values[20];
C[(i)+(32)] += A[(i)+(8)] * values[21];
C[(i)+(32)] += A[(i)+(16)] * values[22];
C[(i)+(32)] += A[(i)+(20)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 48;
#endif

}

void dsparse_starMatrix_m10_n9_k9_ldA12_ldBna4_ldC12_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
#pragma vector aligned
for (unsigned int i = 0; i < 10; i += 1)
{
C[(i)+(0)] += A[(i)+(72)] * values[0];
C[(i)+(0)] += A[(i)+(84)] * values[1];
C[(i)+(0)] += A[(i)+(96)] * values[2];
C[(i)+(12)] += A[(i)+(72)] * values[3];
C[(i)+(12)] += A[(i)+(84)] * values[4];
C[(i)+(12)] += A[(i)+(96)] * values[5];
C[(i)+(24)] += A[(i)+(72)] * values[6];
C[(i)+(24)] += A[(i)+(84)] * values[7];
C[(i)+(24)] += A[(i)+(96)] * values[8];
C[(i)+(36)] += A[(i)+(72)] * values[9];
C[(i)+(36)] += A[(i)+(84)] * values[10];
C[(i)+(48)] += A[(i)+(84)] * values[11];
C[(i)+(48)] += A[(i)+(96)] * values[12];
C[(i)+(60)] += A[(i)+(72)] * values[13];
C[(i)+(60)] += A[(i)+(96)] * values[14];
C[(i)+(72)] += A[(i)+(0)] * values[15];
C[(i)+(72)] += A[(i)+(36)] * values[16];
C[(i)+(72)] += A[(i)+(60)] * values[17];
C[(i)+(84)] += A[(i)+(12)] * values[18];
C[(i)+(84)] += A[(i)+(36)] * values[19];
C[(i)+(84)] += A[(i)+(48)] * values[20];
C[(i)+(96)] += A[(i)+(24)] * values[21];
C[(i)+(96)] += A[(i)+(48)] * values[22];
C[(i)+(96)] += A[(i)+(60)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 480;
#endif

}

void dsparse_starMatrix_m4_n9_k9_ldA4_ldBna4_ldC4_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(4)
#pragma vector aligned
for (unsigned int i = 0; i < 4; i += 1)
{
C[(i)+(0)] += A[(i)+(24)] * values[0];
C[(i)+(0)] += A[(i)+(28)] * values[1];
C[(i)+(0)] += A[(i)+(32)] * values[2];
C[(i)+(4)] += A[(i)+(24)] * values[3];
C[(i)+(4)] += A[(i)+(28)] * values[4];
C[(i)+(4)] += A[(i)+(32)] * values[5];
C[(i)+(8)] += A[(i)+(24)] * values[6];
C[(i)+(8)] += A[(i)+(28)] * values[7];
C[(i)+(8)] += A[(i)+(32)] * values[8];
C[(i)+(12)] += A[(i)+(24)] * values[9];
C[(i)+(12)] += A[(i)+(28)] * values[10];
C[(i)+(16)] += A[(i)+(28)] * values[11];
C[(i)+(16)] += A[(i)+(32)] * values[12];
C[(i)+(20)] += A[(i)+(24)] * values[13];
C[(i)+(20)] += A[(i)+(32)] * values[14];
C[(i)+(24)] += A[(i)+(0)] * values[15];
C[(i)+(24)] += A[(i)+(12)] * values[16];
C[(i)+(24)] += A[(i)+(20)] * values[17];
C[(i)+(28)] += A[(i)+(4)] * values[18];
C[(i)+(28)] += A[(i)+(12)] * values[19];
C[(i)+(28)] += A[(i)+(16)] * values[20];
C[(i)+(32)] += A[(i)+(8)] * values[21];
C[(i)+(32)] += A[(i)+(16)] * values[22];
C[(i)+(32)] += A[(i)+(20)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 192;
#endif

}

void dsparse_starMatrix_m1_n9_k9_ldA4_ldBna4_ldC4_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
for (unsigned int i = 0; i < 1; i += 1)
{
C[(i)+(0)] += A[(i)+(24)] * values[0];
C[(i)+(0)] += A[(i)+(28)] * values[1];
C[(i)+(0)] += A[(i)+(32)] * values[2];
C[(i)+(4)] += A[(i)+(24)] * values[3];
C[(i)+(4)] += A[(i)+(28)] * values[4];
C[(i)+(4)] += A[(i)+(32)] * values[5];
C[(i)+(8)] += A[(i)+(24)] * values[6];
C[(i)+(8)] += A[(i)+(28)] * values[7];
C[(i)+(8)] += A[(i)+(32)] * values[8];
C[(i)+(12)] += A[(i)+(24)] * values[9];
C[(i)+(12)] += A[(i)+(28)] * values[10];
C[(i)+(16)] += A[(i)+(28)] * values[11];
C[(i)+(16)] += A[(i)+(32)] * values[12];
C[(i)+(20)] += A[(i)+(24)] * values[13];
C[(i)+(20)] += A[(i)+(32)] * values[14];
C[(i)+(24)] += A[(i)+(0)] * values[15];
C[(i)+(24)] += A[(i)+(12)] * values[16];
C[(i)+(24)] += A[(i)+(20)] * values[17];
C[(i)+(28)] += A[(i)+(4)] * values[18];
C[(i)+(28)] += A[(i)+(12)] * values[19];
C[(i)+(28)] += A[(i)+(16)] * values[20];
C[(i)+(32)] += A[(i)+(8)] * values[21];
C[(i)+(32)] += A[(i)+(16)] * values[22];
C[(i)+(32)] += A[(i)+(20)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 48;
#endif

}

void dsparse_starMatrix_m20_n9_k9_ldA20_ldBna5_ldC20_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
#pragma vector aligned
for (unsigned int i = 0; i < 20; i += 1)
{
C[(i)+(0)] += A[(i)+(120)] * values[0];
C[(i)+(0)] += A[(i)+(140)] * values[1];
C[(i)+(0)] += A[(i)+(160)] * values[2];
C[(i)+(20)] += A[(i)+(120)] * values[3];
C[(i)+(20)] += A[(i)+(140)] * values[4];
C[(i)+(20)] += A[(i)+(160)] * values[5];
C[(i)+(40)] += A[(i)+(120)] * values[6];
C[(i)+(40)] += A[(i)+(140)] * values[7];
C[(i)+(40)] += A[(i)+(160)] * values[8];
C[(i)+(60)] += A[(i)+(120)] * values[9];
C[(i)+(60)] += A[(i)+(140)] * values[10];
C[(i)+(80)] += A[(i)+(140)] * values[11];
C[(i)+(80)] += A[(i)+(160)] * values[12];
C[(i)+(100)] += A[(i)+(120)] * values[13];
C[(i)+(100)] += A[(i)+(160)] * values[14];
C[(i)+(120)] += A[(i)+(0)] * values[15];
C[(i)+(120)] += A[(i)+(60)] * values[16];
C[(i)+(120)] += A[(i)+(100)] * values[17];
C[(i)+(140)] += A[(i)+(20)] * values[18];
C[(i)+(140)] += A[(i)+(60)] * values[19];
C[(i)+(140)] += A[(i)+(80)] * values[20];
C[(i)+(160)] += A[(i)+(40)] * values[21];
C[(i)+(160)] += A[(i)+(80)] * values[22];
C[(i)+(160)] += A[(i)+(100)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 960;
#endif

}

void dsparse_starMatrix_m10_n9_k9_ldA12_ldBna5_ldC12_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
#pragma vector aligned
for (unsigned int i = 0; i < 10; i += 1)
{
C[(i)+(0)] += A[(i)+(72)] * values[0];
C[(i)+(0)] += A[(i)+(84)] * values[1];
C[(i)+(0)] += A[(i)+(96)] * values[2];
C[(i)+(12)] += A[(i)+(72)] * values[3];
C[(i)+(12)] += A[(i)+(84)] * values[4];
C[(i)+(12)] += A[(i)+(96)] * values[5];
C[(i)+(24)] += A[(i)+(72)] * values[6];
C[(i)+(24)] += A[(i)+(84)] * values[7];
C[(i)+(24)] += A[(i)+(96)] * values[8];
C[(i)+(36)] += A[(i)+(72)] * values[9];
C[(i)+(36)] += A[(i)+(84)] * values[10];
C[(i)+(48)] += A[(i)+(84)] * values[11];
C[(i)+(48)] += A[(i)+(96)] * values[12];
C[(i)+(60)] += A[(i)+(72)] * values[13];
C[(i)+(60)] += A[(i)+(96)] * values[14];
C[(i)+(72)] += A[(i)+(0)] * values[15];
C[(i)+(72)] += A[(i)+(36)] * values[16];
C[(i)+(72)] += A[(i)+(60)] * values[17];
C[(i)+(84)] += A[(i)+(12)] * values[18];
C[(i)+(84)] += A[(i)+(36)] * values[19];
C[(i)+(84)] += A[(i)+(48)] * values[20];
C[(i)+(96)] += A[(i)+(24)] * values[21];
C[(i)+(96)] += A[(i)+(48)] * values[22];
C[(i)+(96)] += A[(i)+(60)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 480;
#endif

}

void dsparse_starMatrix_m4_n9_k9_ldA4_ldBna5_ldC4_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(4)
#pragma vector aligned
for (unsigned int i = 0; i < 4; i += 1)
{
C[(i)+(0)] += A[(i)+(24)] * values[0];
C[(i)+(0)] += A[(i)+(28)] * values[1];
C[(i)+(0)] += A[(i)+(32)] * values[2];
C[(i)+(4)] += A[(i)+(24)] * values[3];
C[(i)+(4)] += A[(i)+(28)] * values[4];
C[(i)+(4)] += A[(i)+(32)] * values[5];
C[(i)+(8)] += A[(i)+(24)] * values[6];
C[(i)+(8)] += A[(i)+(28)] * values[7];
C[(i)+(8)] += A[(i)+(32)] * values[8];
C[(i)+(12)] += A[(i)+(24)] * values[9];
C[(i)+(12)] += A[(i)+(28)] * values[10];
C[(i)+(16)] += A[(i)+(28)] * values[11];
C[(i)+(16)] += A[(i)+(32)] * values[12];
C[(i)+(20)] += A[(i)+(24)] * values[13];
C[(i)+(20)] += A[(i)+(32)] * values[14];
C[(i)+(24)] += A[(i)+(0)] * values[15];
C[(i)+(24)] += A[(i)+(12)] * values[16];
C[(i)+(24)] += A[(i)+(20)] * values[17];
C[(i)+(28)] += A[(i)+(4)] * values[18];
C[(i)+(28)] += A[(i)+(12)] * values[19];
C[(i)+(28)] += A[(i)+(16)] * values[20];
C[(i)+(32)] += A[(i)+(8)] * values[21];
C[(i)+(32)] += A[(i)+(16)] * values[22];
C[(i)+(32)] += A[(i)+(20)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 192;
#endif

}

void dsparse_starMatrix_m1_n9_k9_ldA4_ldBna5_ldC4_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
for (unsigned int i = 0; i < 1; i += 1)
{
C[(i)+(0)] += A[(i)+(24)] * values[0];
C[(i)+(0)] += A[(i)+(28)] * values[1];
C[(i)+(0)] += A[(i)+(32)] * values[2];
C[(i)+(4)] += A[(i)+(24)] * values[3];
C[(i)+(4)] += A[(i)+(28)] * values[4];
C[(i)+(4)] += A[(i)+(32)] * values[5];
C[(i)+(8)] += A[(i)+(24)] * values[6];
C[(i)+(8)] += A[(i)+(28)] * values[7];
C[(i)+(8)] += A[(i)+(32)] * values[8];
C[(i)+(12)] += A[(i)+(24)] * values[9];
C[(i)+(12)] += A[(i)+(28)] * values[10];
C[(i)+(16)] += A[(i)+(28)] * values[11];
C[(i)+(16)] += A[(i)+(32)] * values[12];
C[(i)+(20)] += A[(i)+(24)] * values[13];
C[(i)+(20)] += A[(i)+(32)] * values[14];
C[(i)+(24)] += A[(i)+(0)] * values[15];
C[(i)+(24)] += A[(i)+(12)] * values[16];
C[(i)+(24)] += A[(i)+(20)] * values[17];
C[(i)+(28)] += A[(i)+(4)] * values[18];
C[(i)+(28)] += A[(i)+(12)] * values[19];
C[(i)+(28)] += A[(i)+(16)] * values[20];
C[(i)+(32)] += A[(i)+(8)] * values[21];
C[(i)+(32)] += A[(i)+(16)] * values[22];
C[(i)+(32)] += A[(i)+(20)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 48;
#endif

}

void dsparse_starMatrix_m35_n9_k9_ldA36_ldBna6_ldC36_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
#pragma vector aligned
for (unsigned int i = 0; i < 35; i += 1)
{
C[(i)+(0)] += A[(i)+(216)] * values[0];
C[(i)+(0)] += A[(i)+(252)] * values[1];
C[(i)+(0)] += A[(i)+(288)] * values[2];
C[(i)+(36)] += A[(i)+(216)] * values[3];
C[(i)+(36)] += A[(i)+(252)] * values[4];
C[(i)+(36)] += A[(i)+(288)] * values[5];
C[(i)+(72)] += A[(i)+(216)] * values[6];
C[(i)+(72)] += A[(i)+(252)] * values[7];
C[(i)+(72)] += A[(i)+(288)] * values[8];
C[(i)+(108)] += A[(i)+(216)] * values[9];
C[(i)+(108)] += A[(i)+(252)] * values[10];
C[(i)+(144)] += A[(i)+(252)] * values[11];
C[(i)+(144)] += A[(i)+(288)] * values[12];
C[(i)+(180)] += A[(i)+(216)] * values[13];
C[(i)+(180)] += A[(i)+(288)] * values[14];
C[(i)+(216)] += A[(i)+(0)] * values[15];
C[(i)+(216)] += A[(i)+(108)] * values[16];
C[(i)+(216)] += A[(i)+(180)] * values[17];
C[(i)+(252)] += A[(i)+(36)] * values[18];
C[(i)+(252)] += A[(i)+(108)] * values[19];
C[(i)+(252)] += A[(i)+(144)] * values[20];
C[(i)+(288)] += A[(i)+(72)] * values[21];
C[(i)+(288)] += A[(i)+(144)] * values[22];
C[(i)+(288)] += A[(i)+(180)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 1680;
#endif

}

void dsparse_starMatrix_m20_n9_k9_ldA20_ldBna6_ldC20_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
#pragma vector aligned
for (unsigned int i = 0; i < 20; i += 1)
{
C[(i)+(0)] += A[(i)+(120)] * values[0];
C[(i)+(0)] += A[(i)+(140)] * values[1];
C[(i)+(0)] += A[(i)+(160)] * values[2];
C[(i)+(20)] += A[(i)+(120)] * values[3];
C[(i)+(20)] += A[(i)+(140)] * values[4];
C[(i)+(20)] += A[(i)+(160)] * values[5];
C[(i)+(40)] += A[(i)+(120)] * values[6];
C[(i)+(40)] += A[(i)+(140)] * values[7];
C[(i)+(40)] += A[(i)+(160)] * values[8];
C[(i)+(60)] += A[(i)+(120)] * values[9];
C[(i)+(60)] += A[(i)+(140)] * values[10];
C[(i)+(80)] += A[(i)+(140)] * values[11];
C[(i)+(80)] += A[(i)+(160)] * values[12];
C[(i)+(100)] += A[(i)+(120)] * values[13];
C[(i)+(100)] += A[(i)+(160)] * values[14];
C[(i)+(120)] += A[(i)+(0)] * values[15];
C[(i)+(120)] += A[(i)+(60)] * values[16];
C[(i)+(120)] += A[(i)+(100)] * values[17];
C[(i)+(140)] += A[(i)+(20)] * values[18];
C[(i)+(140)] += A[(i)+(60)] * values[19];
C[(i)+(140)] += A[(i)+(80)] * values[20];
C[(i)+(160)] += A[(i)+(40)] * values[21];
C[(i)+(160)] += A[(i)+(80)] * values[22];
C[(i)+(160)] += A[(i)+(100)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 960;
#endif

}

void dsparse_starMatrix_m10_n9_k9_ldA12_ldBna6_ldC12_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
#pragma vector aligned
for (unsigned int i = 0; i < 10; i += 1)
{
C[(i)+(0)] += A[(i)+(72)] * values[0];
C[(i)+(0)] += A[(i)+(84)] * values[1];
C[(i)+(0)] += A[(i)+(96)] * values[2];
C[(i)+(12)] += A[(i)+(72)] * values[3];
C[(i)+(12)] += A[(i)+(84)] * values[4];
C[(i)+(12)] += A[(i)+(96)] * values[5];
C[(i)+(24)] += A[(i)+(72)] * values[6];
C[(i)+(24)] += A[(i)+(84)] * values[7];
C[(i)+(24)] += A[(i)+(96)] * values[8];
C[(i)+(36)] += A[(i)+(72)] * values[9];
C[(i)+(36)] += A[(i)+(84)] * values[10];
C[(i)+(48)] += A[(i)+(84)] * values[11];
C[(i)+(48)] += A[(i)+(96)] * values[12];
C[(i)+(60)] += A[(i)+(72)] * values[13];
C[(i)+(60)] += A[(i)+(96)] * values[14];
C[(i)+(72)] += A[(i)+(0)] * values[15];
C[(i)+(72)] += A[(i)+(36)] * values[16];
C[(i)+(72)] += A[(i)+(60)] * values[17];
C[(i)+(84)] += A[(i)+(12)] * values[18];
C[(i)+(84)] += A[(i)+(36)] * values[19];
C[(i)+(84)] += A[(i)+(48)] * values[20];
C[(i)+(96)] += A[(i)+(24)] * values[21];
C[(i)+(96)] += A[(i)+(48)] * values[22];
C[(i)+(96)] += A[(i)+(60)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 480;
#endif

}

void dsparse_starMatrix_m4_n9_k9_ldA4_ldBna6_ldC4_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(4)
#pragma vector aligned
for (unsigned int i = 0; i < 4; i += 1)
{
C[(i)+(0)] += A[(i)+(24)] * values[0];
C[(i)+(0)] += A[(i)+(28)] * values[1];
C[(i)+(0)] += A[(i)+(32)] * values[2];
C[(i)+(4)] += A[(i)+(24)] * values[3];
C[(i)+(4)] += A[(i)+(28)] * values[4];
C[(i)+(4)] += A[(i)+(32)] * values[5];
C[(i)+(8)] += A[(i)+(24)] * values[6];
C[(i)+(8)] += A[(i)+(28)] * values[7];
C[(i)+(8)] += A[(i)+(32)] * values[8];
C[(i)+(12)] += A[(i)+(24)] * values[9];
C[(i)+(12)] += A[(i)+(28)] * values[10];
C[(i)+(16)] += A[(i)+(28)] * values[11];
C[(i)+(16)] += A[(i)+(32)] * values[12];
C[(i)+(20)] += A[(i)+(24)] * values[13];
C[(i)+(20)] += A[(i)+(32)] * values[14];
C[(i)+(24)] += A[(i)+(0)] * values[15];
C[(i)+(24)] += A[(i)+(12)] * values[16];
C[(i)+(24)] += A[(i)+(20)] * values[17];
C[(i)+(28)] += A[(i)+(4)] * values[18];
C[(i)+(28)] += A[(i)+(12)] * values[19];
C[(i)+(28)] += A[(i)+(16)] * values[20];
C[(i)+(32)] += A[(i)+(8)] * values[21];
C[(i)+(32)] += A[(i)+(16)] * values[22];
C[(i)+(32)] += A[(i)+(20)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 192;
#endif

}

void dsparse_starMatrix_m1_n9_k9_ldA4_ldBna6_ldC4_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
for (unsigned int i = 0; i < 1; i += 1)
{
C[(i)+(0)] += A[(i)+(24)] * values[0];
C[(i)+(0)] += A[(i)+(28)] * values[1];
C[(i)+(0)] += A[(i)+(32)] * values[2];
C[(i)+(4)] += A[(i)+(24)] * values[3];
C[(i)+(4)] += A[(i)+(28)] * values[4];
C[(i)+(4)] += A[(i)+(32)] * values[5];
C[(i)+(8)] += A[(i)+(24)] * values[6];
C[(i)+(8)] += A[(i)+(28)] * values[7];
C[(i)+(8)] += A[(i)+(32)] * values[8];
C[(i)+(12)] += A[(i)+(24)] * values[9];
C[(i)+(12)] += A[(i)+(28)] * values[10];
C[(i)+(16)] += A[(i)+(28)] * values[11];
C[(i)+(16)] += A[(i)+(32)] * values[12];
C[(i)+(20)] += A[(i)+(24)] * values[13];
C[(i)+(20)] += A[(i)+(32)] * values[14];
C[(i)+(24)] += A[(i)+(0)] * values[15];
C[(i)+(24)] += A[(i)+(12)] * values[16];
C[(i)+(24)] += A[(i)+(20)] * values[17];
C[(i)+(28)] += A[(i)+(4)] * values[18];
C[(i)+(28)] += A[(i)+(12)] * values[19];
C[(i)+(28)] += A[(i)+(16)] * values[20];
C[(i)+(32)] += A[(i)+(8)] * values[21];
C[(i)+(32)] += A[(i)+(16)] * values[22];
C[(i)+(32)] += A[(i)+(20)] * values[23];
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

void dsparse_starMatrix_m35_n9_k9_ldA36_ldBna7_ldC36_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
#pragma vector aligned
for (unsigned int i = 0; i < 35; i += 1)
{
C[(i)+(0)] += A[(i)+(216)] * values[0];
C[(i)+(0)] += A[(i)+(252)] * values[1];
C[(i)+(0)] += A[(i)+(288)] * values[2];
C[(i)+(36)] += A[(i)+(216)] * values[3];
C[(i)+(36)] += A[(i)+(252)] * values[4];
C[(i)+(36)] += A[(i)+(288)] * values[5];
C[(i)+(72)] += A[(i)+(216)] * values[6];
C[(i)+(72)] += A[(i)+(252)] * values[7];
C[(i)+(72)] += A[(i)+(288)] * values[8];
C[(i)+(108)] += A[(i)+(216)] * values[9];
C[(i)+(108)] += A[(i)+(252)] * values[10];
C[(i)+(144)] += A[(i)+(252)] * values[11];
C[(i)+(144)] += A[(i)+(288)] * values[12];
C[(i)+(180)] += A[(i)+(216)] * values[13];
C[(i)+(180)] += A[(i)+(288)] * values[14];
C[(i)+(216)] += A[(i)+(0)] * values[15];
C[(i)+(216)] += A[(i)+(108)] * values[16];
C[(i)+(216)] += A[(i)+(180)] * values[17];
C[(i)+(252)] += A[(i)+(36)] * values[18];
C[(i)+(252)] += A[(i)+(108)] * values[19];
C[(i)+(252)] += A[(i)+(144)] * values[20];
C[(i)+(288)] += A[(i)+(72)] * values[21];
C[(i)+(288)] += A[(i)+(144)] * values[22];
C[(i)+(288)] += A[(i)+(180)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 1680;
#endif

}

void dsparse_starMatrix_m20_n9_k9_ldA20_ldBna7_ldC20_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
#pragma vector aligned
for (unsigned int i = 0; i < 20; i += 1)
{
C[(i)+(0)] += A[(i)+(120)] * values[0];
C[(i)+(0)] += A[(i)+(140)] * values[1];
C[(i)+(0)] += A[(i)+(160)] * values[2];
C[(i)+(20)] += A[(i)+(120)] * values[3];
C[(i)+(20)] += A[(i)+(140)] * values[4];
C[(i)+(20)] += A[(i)+(160)] * values[5];
C[(i)+(40)] += A[(i)+(120)] * values[6];
C[(i)+(40)] += A[(i)+(140)] * values[7];
C[(i)+(40)] += A[(i)+(160)] * values[8];
C[(i)+(60)] += A[(i)+(120)] * values[9];
C[(i)+(60)] += A[(i)+(140)] * values[10];
C[(i)+(80)] += A[(i)+(140)] * values[11];
C[(i)+(80)] += A[(i)+(160)] * values[12];
C[(i)+(100)] += A[(i)+(120)] * values[13];
C[(i)+(100)] += A[(i)+(160)] * values[14];
C[(i)+(120)] += A[(i)+(0)] * values[15];
C[(i)+(120)] += A[(i)+(60)] * values[16];
C[(i)+(120)] += A[(i)+(100)] * values[17];
C[(i)+(140)] += A[(i)+(20)] * values[18];
C[(i)+(140)] += A[(i)+(60)] * values[19];
C[(i)+(140)] += A[(i)+(80)] * values[20];
C[(i)+(160)] += A[(i)+(40)] * values[21];
C[(i)+(160)] += A[(i)+(80)] * values[22];
C[(i)+(160)] += A[(i)+(100)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 960;
#endif

}

void dsparse_starMatrix_m10_n9_k9_ldA12_ldBna7_ldC12_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
#pragma vector aligned
for (unsigned int i = 0; i < 10; i += 1)
{
C[(i)+(0)] += A[(i)+(72)] * values[0];
C[(i)+(0)] += A[(i)+(84)] * values[1];
C[(i)+(0)] += A[(i)+(96)] * values[2];
C[(i)+(12)] += A[(i)+(72)] * values[3];
C[(i)+(12)] += A[(i)+(84)] * values[4];
C[(i)+(12)] += A[(i)+(96)] * values[5];
C[(i)+(24)] += A[(i)+(72)] * values[6];
C[(i)+(24)] += A[(i)+(84)] * values[7];
C[(i)+(24)] += A[(i)+(96)] * values[8];
C[(i)+(36)] += A[(i)+(72)] * values[9];
C[(i)+(36)] += A[(i)+(84)] * values[10];
C[(i)+(48)] += A[(i)+(84)] * values[11];
C[(i)+(48)] += A[(i)+(96)] * values[12];
C[(i)+(60)] += A[(i)+(72)] * values[13];
C[(i)+(60)] += A[(i)+(96)] * values[14];
C[(i)+(72)] += A[(i)+(0)] * values[15];
C[(i)+(72)] += A[(i)+(36)] * values[16];
C[(i)+(72)] += A[(i)+(60)] * values[17];
C[(i)+(84)] += A[(i)+(12)] * values[18];
C[(i)+(84)] += A[(i)+(36)] * values[19];
C[(i)+(84)] += A[(i)+(48)] * values[20];
C[(i)+(96)] += A[(i)+(24)] * values[21];
C[(i)+(96)] += A[(i)+(48)] * values[22];
C[(i)+(96)] += A[(i)+(60)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 480;
#endif

}

void dsparse_starMatrix_m4_n9_k9_ldA4_ldBna7_ldC4_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(4)
#pragma vector aligned
for (unsigned int i = 0; i < 4; i += 1)
{
C[(i)+(0)] += A[(i)+(24)] * values[0];
C[(i)+(0)] += A[(i)+(28)] * values[1];
C[(i)+(0)] += A[(i)+(32)] * values[2];
C[(i)+(4)] += A[(i)+(24)] * values[3];
C[(i)+(4)] += A[(i)+(28)] * values[4];
C[(i)+(4)] += A[(i)+(32)] * values[5];
C[(i)+(8)] += A[(i)+(24)] * values[6];
C[(i)+(8)] += A[(i)+(28)] * values[7];
C[(i)+(8)] += A[(i)+(32)] * values[8];
C[(i)+(12)] += A[(i)+(24)] * values[9];
C[(i)+(12)] += A[(i)+(28)] * values[10];
C[(i)+(16)] += A[(i)+(28)] * values[11];
C[(i)+(16)] += A[(i)+(32)] * values[12];
C[(i)+(20)] += A[(i)+(24)] * values[13];
C[(i)+(20)] += A[(i)+(32)] * values[14];
C[(i)+(24)] += A[(i)+(0)] * values[15];
C[(i)+(24)] += A[(i)+(12)] * values[16];
C[(i)+(24)] += A[(i)+(20)] * values[17];
C[(i)+(28)] += A[(i)+(4)] * values[18];
C[(i)+(28)] += A[(i)+(12)] * values[19];
C[(i)+(28)] += A[(i)+(16)] * values[20];
C[(i)+(32)] += A[(i)+(8)] * values[21];
C[(i)+(32)] += A[(i)+(16)] * values[22];
C[(i)+(32)] += A[(i)+(20)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 192;
#endif

}

void dsparse_starMatrix_m1_n9_k9_ldA4_ldBna7_ldC4_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
for (unsigned int i = 0; i < 1; i += 1)
{
C[(i)+(0)] += A[(i)+(24)] * values[0];
C[(i)+(0)] += A[(i)+(28)] * values[1];
C[(i)+(0)] += A[(i)+(32)] * values[2];
C[(i)+(4)] += A[(i)+(24)] * values[3];
C[(i)+(4)] += A[(i)+(28)] * values[4];
C[(i)+(4)] += A[(i)+(32)] * values[5];
C[(i)+(8)] += A[(i)+(24)] * values[6];
C[(i)+(8)] += A[(i)+(28)] * values[7];
C[(i)+(8)] += A[(i)+(32)] * values[8];
C[(i)+(12)] += A[(i)+(24)] * values[9];
C[(i)+(12)] += A[(i)+(28)] * values[10];
C[(i)+(16)] += A[(i)+(28)] * values[11];
C[(i)+(16)] += A[(i)+(32)] * values[12];
C[(i)+(20)] += A[(i)+(24)] * values[13];
C[(i)+(20)] += A[(i)+(32)] * values[14];
C[(i)+(24)] += A[(i)+(0)] * values[15];
C[(i)+(24)] += A[(i)+(12)] * values[16];
C[(i)+(24)] += A[(i)+(20)] * values[17];
C[(i)+(28)] += A[(i)+(4)] * values[18];
C[(i)+(28)] += A[(i)+(12)] * values[19];
C[(i)+(28)] += A[(i)+(16)] * values[20];
C[(i)+(32)] += A[(i)+(8)] * values[21];
C[(i)+(32)] += A[(i)+(16)] * values[22];
C[(i)+(32)] += A[(i)+(20)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 48;
#endif

}

void dsparse_starMatrix_m84_n9_k9_ldA84_ldBna8_ldC84_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
#pragma vector aligned
for (unsigned int i = 0; i < 84; i += 1)
{
C[(i)+(0)] += A[(i)+(504)] * values[0];
C[(i)+(0)] += A[(i)+(588)] * values[1];
C[(i)+(0)] += A[(i)+(672)] * values[2];
C[(i)+(84)] += A[(i)+(504)] * values[3];
C[(i)+(84)] += A[(i)+(588)] * values[4];
C[(i)+(84)] += A[(i)+(672)] * values[5];
C[(i)+(168)] += A[(i)+(504)] * values[6];
C[(i)+(168)] += A[(i)+(588)] * values[7];
C[(i)+(168)] += A[(i)+(672)] * values[8];
C[(i)+(252)] += A[(i)+(504)] * values[9];
C[(i)+(252)] += A[(i)+(588)] * values[10];
C[(i)+(336)] += A[(i)+(588)] * values[11];
C[(i)+(336)] += A[(i)+(672)] * values[12];
C[(i)+(420)] += A[(i)+(504)] * values[13];
C[(i)+(420)] += A[(i)+(672)] * values[14];
C[(i)+(504)] += A[(i)+(0)] * values[15];
C[(i)+(504)] += A[(i)+(252)] * values[16];
C[(i)+(504)] += A[(i)+(420)] * values[17];
C[(i)+(588)] += A[(i)+(84)] * values[18];
C[(i)+(588)] += A[(i)+(252)] * values[19];
C[(i)+(588)] += A[(i)+(336)] * values[20];
C[(i)+(672)] += A[(i)+(168)] * values[21];
C[(i)+(672)] += A[(i)+(336)] * values[22];
C[(i)+(672)] += A[(i)+(420)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 4032;
#endif

}

void dsparse_starMatrix_m56_n9_k9_ldA56_ldBna8_ldC56_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
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

void dsparse_starMatrix_m35_n9_k9_ldA36_ldBna8_ldC36_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
#pragma vector aligned
for (unsigned int i = 0; i < 35; i += 1)
{
C[(i)+(0)] += A[(i)+(216)] * values[0];
C[(i)+(0)] += A[(i)+(252)] * values[1];
C[(i)+(0)] += A[(i)+(288)] * values[2];
C[(i)+(36)] += A[(i)+(216)] * values[3];
C[(i)+(36)] += A[(i)+(252)] * values[4];
C[(i)+(36)] += A[(i)+(288)] * values[5];
C[(i)+(72)] += A[(i)+(216)] * values[6];
C[(i)+(72)] += A[(i)+(252)] * values[7];
C[(i)+(72)] += A[(i)+(288)] * values[8];
C[(i)+(108)] += A[(i)+(216)] * values[9];
C[(i)+(108)] += A[(i)+(252)] * values[10];
C[(i)+(144)] += A[(i)+(252)] * values[11];
C[(i)+(144)] += A[(i)+(288)] * values[12];
C[(i)+(180)] += A[(i)+(216)] * values[13];
C[(i)+(180)] += A[(i)+(288)] * values[14];
C[(i)+(216)] += A[(i)+(0)] * values[15];
C[(i)+(216)] += A[(i)+(108)] * values[16];
C[(i)+(216)] += A[(i)+(180)] * values[17];
C[(i)+(252)] += A[(i)+(36)] * values[18];
C[(i)+(252)] += A[(i)+(108)] * values[19];
C[(i)+(252)] += A[(i)+(144)] * values[20];
C[(i)+(288)] += A[(i)+(72)] * values[21];
C[(i)+(288)] += A[(i)+(144)] * values[22];
C[(i)+(288)] += A[(i)+(180)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 1680;
#endif

}

void dsparse_starMatrix_m20_n9_k9_ldA20_ldBna8_ldC20_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
#pragma vector aligned
for (unsigned int i = 0; i < 20; i += 1)
{
C[(i)+(0)] += A[(i)+(120)] * values[0];
C[(i)+(0)] += A[(i)+(140)] * values[1];
C[(i)+(0)] += A[(i)+(160)] * values[2];
C[(i)+(20)] += A[(i)+(120)] * values[3];
C[(i)+(20)] += A[(i)+(140)] * values[4];
C[(i)+(20)] += A[(i)+(160)] * values[5];
C[(i)+(40)] += A[(i)+(120)] * values[6];
C[(i)+(40)] += A[(i)+(140)] * values[7];
C[(i)+(40)] += A[(i)+(160)] * values[8];
C[(i)+(60)] += A[(i)+(120)] * values[9];
C[(i)+(60)] += A[(i)+(140)] * values[10];
C[(i)+(80)] += A[(i)+(140)] * values[11];
C[(i)+(80)] += A[(i)+(160)] * values[12];
C[(i)+(100)] += A[(i)+(120)] * values[13];
C[(i)+(100)] += A[(i)+(160)] * values[14];
C[(i)+(120)] += A[(i)+(0)] * values[15];
C[(i)+(120)] += A[(i)+(60)] * values[16];
C[(i)+(120)] += A[(i)+(100)] * values[17];
C[(i)+(140)] += A[(i)+(20)] * values[18];
C[(i)+(140)] += A[(i)+(60)] * values[19];
C[(i)+(140)] += A[(i)+(80)] * values[20];
C[(i)+(160)] += A[(i)+(40)] * values[21];
C[(i)+(160)] += A[(i)+(80)] * values[22];
C[(i)+(160)] += A[(i)+(100)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 960;
#endif

}

void dsparse_starMatrix_m10_n9_k9_ldA12_ldBna8_ldC12_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
#pragma vector aligned
for (unsigned int i = 0; i < 10; i += 1)
{
C[(i)+(0)] += A[(i)+(72)] * values[0];
C[(i)+(0)] += A[(i)+(84)] * values[1];
C[(i)+(0)] += A[(i)+(96)] * values[2];
C[(i)+(12)] += A[(i)+(72)] * values[3];
C[(i)+(12)] += A[(i)+(84)] * values[4];
C[(i)+(12)] += A[(i)+(96)] * values[5];
C[(i)+(24)] += A[(i)+(72)] * values[6];
C[(i)+(24)] += A[(i)+(84)] * values[7];
C[(i)+(24)] += A[(i)+(96)] * values[8];
C[(i)+(36)] += A[(i)+(72)] * values[9];
C[(i)+(36)] += A[(i)+(84)] * values[10];
C[(i)+(48)] += A[(i)+(84)] * values[11];
C[(i)+(48)] += A[(i)+(96)] * values[12];
C[(i)+(60)] += A[(i)+(72)] * values[13];
C[(i)+(60)] += A[(i)+(96)] * values[14];
C[(i)+(72)] += A[(i)+(0)] * values[15];
C[(i)+(72)] += A[(i)+(36)] * values[16];
C[(i)+(72)] += A[(i)+(60)] * values[17];
C[(i)+(84)] += A[(i)+(12)] * values[18];
C[(i)+(84)] += A[(i)+(36)] * values[19];
C[(i)+(84)] += A[(i)+(48)] * values[20];
C[(i)+(96)] += A[(i)+(24)] * values[21];
C[(i)+(96)] += A[(i)+(48)] * values[22];
C[(i)+(96)] += A[(i)+(60)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 480;
#endif

}

void dsparse_starMatrix_m4_n9_k9_ldA4_ldBna8_ldC4_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(4)
#pragma vector aligned
for (unsigned int i = 0; i < 4; i += 1)
{
C[(i)+(0)] += A[(i)+(24)] * values[0];
C[(i)+(0)] += A[(i)+(28)] * values[1];
C[(i)+(0)] += A[(i)+(32)] * values[2];
C[(i)+(4)] += A[(i)+(24)] * values[3];
C[(i)+(4)] += A[(i)+(28)] * values[4];
C[(i)+(4)] += A[(i)+(32)] * values[5];
C[(i)+(8)] += A[(i)+(24)] * values[6];
C[(i)+(8)] += A[(i)+(28)] * values[7];
C[(i)+(8)] += A[(i)+(32)] * values[8];
C[(i)+(12)] += A[(i)+(24)] * values[9];
C[(i)+(12)] += A[(i)+(28)] * values[10];
C[(i)+(16)] += A[(i)+(28)] * values[11];
C[(i)+(16)] += A[(i)+(32)] * values[12];
C[(i)+(20)] += A[(i)+(24)] * values[13];
C[(i)+(20)] += A[(i)+(32)] * values[14];
C[(i)+(24)] += A[(i)+(0)] * values[15];
C[(i)+(24)] += A[(i)+(12)] * values[16];
C[(i)+(24)] += A[(i)+(20)] * values[17];
C[(i)+(28)] += A[(i)+(4)] * values[18];
C[(i)+(28)] += A[(i)+(12)] * values[19];
C[(i)+(28)] += A[(i)+(16)] * values[20];
C[(i)+(32)] += A[(i)+(8)] * values[21];
C[(i)+(32)] += A[(i)+(16)] * values[22];
C[(i)+(32)] += A[(i)+(20)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 192;
#endif

}

void dsparse_starMatrix_m1_n9_k9_ldA4_ldBna8_ldC4_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
for (unsigned int i = 0; i < 1; i += 1)
{
C[(i)+(0)] += A[(i)+(24)] * values[0];
C[(i)+(0)] += A[(i)+(28)] * values[1];
C[(i)+(0)] += A[(i)+(32)] * values[2];
C[(i)+(4)] += A[(i)+(24)] * values[3];
C[(i)+(4)] += A[(i)+(28)] * values[4];
C[(i)+(4)] += A[(i)+(32)] * values[5];
C[(i)+(8)] += A[(i)+(24)] * values[6];
C[(i)+(8)] += A[(i)+(28)] * values[7];
C[(i)+(8)] += A[(i)+(32)] * values[8];
C[(i)+(12)] += A[(i)+(24)] * values[9];
C[(i)+(12)] += A[(i)+(28)] * values[10];
C[(i)+(16)] += A[(i)+(28)] * values[11];
C[(i)+(16)] += A[(i)+(32)] * values[12];
C[(i)+(20)] += A[(i)+(24)] * values[13];
C[(i)+(20)] += A[(i)+(32)] * values[14];
C[(i)+(24)] += A[(i)+(0)] * values[15];
C[(i)+(24)] += A[(i)+(12)] * values[16];
C[(i)+(24)] += A[(i)+(20)] * values[17];
C[(i)+(28)] += A[(i)+(4)] * values[18];
C[(i)+(28)] += A[(i)+(12)] * values[19];
C[(i)+(28)] += A[(i)+(16)] * values[20];
C[(i)+(32)] += A[(i)+(8)] * values[21];
C[(i)+(32)] += A[(i)+(16)] * values[22];
C[(i)+(32)] += A[(i)+(20)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 48;
#endif

}

void dsparse_starMatrix_m4_n9_k9_ldA4_ldBna2_ldC4_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(4)
#pragma vector aligned
for (unsigned int i = 0; i < 4; i += 1)
{
C[(i)+(0)] += A[(i)+(24)] * values[0];
C[(i)+(0)] += A[(i)+(28)] * values[1];
C[(i)+(0)] += A[(i)+(32)] * values[2];
C[(i)+(4)] += A[(i)+(24)] * values[3];
C[(i)+(4)] += A[(i)+(28)] * values[4];
C[(i)+(4)] += A[(i)+(32)] * values[5];
C[(i)+(8)] += A[(i)+(24)] * values[6];
C[(i)+(8)] += A[(i)+(28)] * values[7];
C[(i)+(8)] += A[(i)+(32)] * values[8];
C[(i)+(12)] += A[(i)+(24)] * values[9];
C[(i)+(12)] += A[(i)+(28)] * values[10];
C[(i)+(16)] += A[(i)+(28)] * values[11];
C[(i)+(16)] += A[(i)+(32)] * values[12];
C[(i)+(20)] += A[(i)+(24)] * values[13];
C[(i)+(20)] += A[(i)+(32)] * values[14];
C[(i)+(24)] += A[(i)+(0)] * values[15];
C[(i)+(24)] += A[(i)+(12)] * values[16];
C[(i)+(24)] += A[(i)+(20)] * values[17];
C[(i)+(28)] += A[(i)+(4)] * values[18];
C[(i)+(28)] += A[(i)+(12)] * values[19];
C[(i)+(28)] += A[(i)+(16)] * values[20];
C[(i)+(32)] += A[(i)+(8)] * values[21];
C[(i)+(32)] += A[(i)+(16)] * values[22];
C[(i)+(32)] += A[(i)+(20)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 192;
#endif

}

void dsparse_starMatrix_m10_n9_k9_ldA12_ldBna3_ldC12_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
#pragma vector aligned
for (unsigned int i = 0; i < 10; i += 1)
{
C[(i)+(0)] += A[(i)+(72)] * values[0];
C[(i)+(0)] += A[(i)+(84)] * values[1];
C[(i)+(0)] += A[(i)+(96)] * values[2];
C[(i)+(12)] += A[(i)+(72)] * values[3];
C[(i)+(12)] += A[(i)+(84)] * values[4];
C[(i)+(12)] += A[(i)+(96)] * values[5];
C[(i)+(24)] += A[(i)+(72)] * values[6];
C[(i)+(24)] += A[(i)+(84)] * values[7];
C[(i)+(24)] += A[(i)+(96)] * values[8];
C[(i)+(36)] += A[(i)+(72)] * values[9];
C[(i)+(36)] += A[(i)+(84)] * values[10];
C[(i)+(48)] += A[(i)+(84)] * values[11];
C[(i)+(48)] += A[(i)+(96)] * values[12];
C[(i)+(60)] += A[(i)+(72)] * values[13];
C[(i)+(60)] += A[(i)+(96)] * values[14];
C[(i)+(72)] += A[(i)+(0)] * values[15];
C[(i)+(72)] += A[(i)+(36)] * values[16];
C[(i)+(72)] += A[(i)+(60)] * values[17];
C[(i)+(84)] += A[(i)+(12)] * values[18];
C[(i)+(84)] += A[(i)+(36)] * values[19];
C[(i)+(84)] += A[(i)+(48)] * values[20];
C[(i)+(96)] += A[(i)+(24)] * values[21];
C[(i)+(96)] += A[(i)+(48)] * values[22];
C[(i)+(96)] += A[(i)+(60)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 480;
#endif

}

void dsparse_fP113DivM_m10_n9_k10_ldAna3_ldB12_ldC12_beta0_pfsigonly(const double* values, const double* B, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma nounroll_and_jam
for (unsigned int i = 0; i < 9; i++)
{
  #pragma simd
  for (unsigned int m = 0; m < 10; m++) {
    C[(i*12)+m] = 0.0;
  }
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b0 = _mm256_broadcast_sd(&B[(i*12)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b0 = _mm_loaddup_pd(&B[(i*12)+0]);
#endif
__m128d c0_0 = _mm_load_sd(&C[(i*12)+0]);
__m128d a0_0 = _mm_load_sd(&values[0]);
#if defined(__SSE3__) && defined(__AVX__)
c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, b0));
#endif
_mm_store_sd(&C[(i*12)+0], c0_0);
__m128d c0_1 = _mm_load_sd(&C[(i*12)+3]);
__m128d a0_1 = _mm_load_sd(&values[1]);
#if defined(__SSE3__) && defined(__AVX__)
c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, b0));
#endif
_mm_store_sd(&C[(i*12)+3], c0_1);
__m128d c0_2 = _mm_load_sd(&C[(i*12)+9]);
__m128d a0_2 = _mm_load_sd(&values[2]);
#if defined(__SSE3__) && defined(__AVX__)
c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, b0));
#endif
_mm_store_sd(&C[(i*12)+9], c0_2);
#else
C[(i*12)+0] += values[0] * B[(i*12)+0];
C[(i*12)+3] += values[1] * B[(i*12)+0];
C[(i*12)+9] += values[2] * B[(i*12)+0];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b1 = _mm256_broadcast_sd(&B[(i*12)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b1 = _mm_loaddup_pd(&B[(i*12)+1]);
#endif
__m128d c1_0 = _mm_load_sd(&C[(i*12)+1]);
__m128d a1_0 = _mm_load_sd(&values[3]);
#if defined(__SSE3__) && defined(__AVX__)
c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, b1));
#endif
_mm_store_sd(&C[(i*12)+1], c1_0);
__m128d c1_1 = _mm_load_sd(&C[(i*12)+7]);
__m128d a1_1 = _mm_load_sd(&values[4]);
#if defined(__SSE3__) && defined(__AVX__)
c1_1 = _mm_add_sd(c1_1, _mm_mul_sd(a1_1, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_1 = _mm_add_sd(c1_1, _mm_mul_sd(a1_1, b1));
#endif
_mm_store_sd(&C[(i*12)+7], c1_1);
#else
C[(i*12)+1] += values[3] * B[(i*12)+1];
C[(i*12)+7] += values[4] * B[(i*12)+1];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b2 = _mm256_broadcast_sd(&B[(i*12)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b2 = _mm_loaddup_pd(&B[(i*12)+2]);
#endif
__m128d c2_0 = _mm_load_sd(&C[(i*12)+2]);
__m128d a2_0 = _mm_load_sd(&values[5]);
#if defined(__SSE3__) && defined(__AVX__)
c2_0 = _mm_add_sd(c2_0, _mm_mul_sd(a2_0, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_0 = _mm_add_sd(c2_0, _mm_mul_sd(a2_0, b2));
#endif
_mm_store_sd(&C[(i*12)+2], c2_0);
__m128d c2_1 = _mm_load_sd(&C[(i*12)+8]);
__m128d a2_1 = _mm_load_sd(&values[6]);
#if defined(__SSE3__) && defined(__AVX__)
c2_1 = _mm_add_sd(c2_1, _mm_mul_sd(a2_1, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_1 = _mm_add_sd(c2_1, _mm_mul_sd(a2_1, b2));
#endif
_mm_store_sd(&C[(i*12)+8], c2_1);
#else
C[(i*12)+2] += values[5] * B[(i*12)+2];
C[(i*12)+8] += values[6] * B[(i*12)+2];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b3 = _mm256_broadcast_sd(&B[(i*12)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b3 = _mm_loaddup_pd(&B[(i*12)+3]);
#endif
__m128d c3_0 = _mm_load_sd(&C[(i*12)+0]);
__m128d a3_0 = _mm_load_sd(&values[7]);
#if defined(__SSE3__) && defined(__AVX__)
c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, b3));
#endif
_mm_store_sd(&C[(i*12)+0], c3_0);
__m128d c3_1 = _mm_load_sd(&C[(i*12)+3]);
__m128d a3_1 = _mm_load_sd(&values[8]);
#if defined(__SSE3__) && defined(__AVX__)
c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, b3));
#endif
_mm_store_sd(&C[(i*12)+3], c3_1);
__m128d c3_2 = _mm_load_sd(&C[(i*12)+9]);
__m128d a3_2 = _mm_load_sd(&values[9]);
#if defined(__SSE3__) && defined(__AVX__)
c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, b3));
#endif
_mm_store_sd(&C[(i*12)+9], c3_2);
#else
C[(i*12)+0] += values[7] * B[(i*12)+3];
C[(i*12)+3] += values[8] * B[(i*12)+3];
C[(i*12)+9] += values[9] * B[(i*12)+3];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b4 = _mm256_broadcast_sd(&B[(i*12)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b4 = _mm_loaddup_pd(&B[(i*12)+4]);
#endif
__m128d c4_0 = _mm_load_sd(&C[(i*12)+4]);
__m128d a4_0 = _mm_load_sd(&values[10]);
#if defined(__SSE3__) && defined(__AVX__)
c4_0 = _mm_add_sd(c4_0, _mm_mul_sd(a4_0, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_0 = _mm_add_sd(c4_0, _mm_mul_sd(a4_0, b4));
#endif
_mm_store_sd(&C[(i*12)+4], c4_0);
#else
C[(i*12)+4] += values[10] * B[(i*12)+4];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b5 = _mm256_broadcast_sd(&B[(i*12)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b5 = _mm_loaddup_pd(&B[(i*12)+5]);
#endif
__m128d c5_0 = _mm_load_sd(&C[(i*12)+5]);
__m128d a5_0 = _mm_load_sd(&values[11]);
#if defined(__SSE3__) && defined(__AVX__)
c5_0 = _mm_add_sd(c5_0, _mm_mul_sd(a5_0, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_0 = _mm_add_sd(c5_0, _mm_mul_sd(a5_0, b5));
#endif
_mm_store_sd(&C[(i*12)+5], c5_0);
#else
C[(i*12)+5] += values[11] * B[(i*12)+5];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b6 = _mm256_broadcast_sd(&B[(i*12)+6]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b6 = _mm_loaddup_pd(&B[(i*12)+6]);
#endif
__m128d c6_0 = _mm_load_sd(&C[(i*12)+6]);
__m128d a6_0 = _mm_load_sd(&values[12]);
#if defined(__SSE3__) && defined(__AVX__)
c6_0 = _mm_add_sd(c6_0, _mm_mul_sd(a6_0, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_0 = _mm_add_sd(c6_0, _mm_mul_sd(a6_0, b6));
#endif
_mm_store_sd(&C[(i*12)+6], c6_0);
#else
C[(i*12)+6] += values[12] * B[(i*12)+6];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b7 = _mm256_broadcast_sd(&B[(i*12)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b7 = _mm_loaddup_pd(&B[(i*12)+7]);
#endif
__m128d c7_0 = _mm_load_sd(&C[(i*12)+1]);
__m128d a7_0 = _mm_load_sd(&values[13]);
#if defined(__SSE3__) && defined(__AVX__)
c7_0 = _mm_add_sd(c7_0, _mm_mul_sd(a7_0, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_0 = _mm_add_sd(c7_0, _mm_mul_sd(a7_0, b7));
#endif
_mm_store_sd(&C[(i*12)+1], c7_0);
__m128d c7_1 = _mm_load_sd(&C[(i*12)+7]);
__m128d a7_1 = _mm_load_sd(&values[14]);
#if defined(__SSE3__) && defined(__AVX__)
c7_1 = _mm_add_sd(c7_1, _mm_mul_sd(a7_1, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_1 = _mm_add_sd(c7_1, _mm_mul_sd(a7_1, b7));
#endif
_mm_store_sd(&C[(i*12)+7], c7_1);
#else
C[(i*12)+1] += values[13] * B[(i*12)+7];
C[(i*12)+7] += values[14] * B[(i*12)+7];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b8 = _mm256_broadcast_sd(&B[(i*12)+8]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b8 = _mm_loaddup_pd(&B[(i*12)+8]);
#endif
__m128d c8_0 = _mm_load_sd(&C[(i*12)+2]);
__m128d a8_0 = _mm_load_sd(&values[15]);
#if defined(__SSE3__) && defined(__AVX__)
c8_0 = _mm_add_sd(c8_0, _mm_mul_sd(a8_0, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_0 = _mm_add_sd(c8_0, _mm_mul_sd(a8_0, b8));
#endif
_mm_store_sd(&C[(i*12)+2], c8_0);
__m128d c8_1 = _mm_load_sd(&C[(i*12)+8]);
__m128d a8_1 = _mm_load_sd(&values[16]);
#if defined(__SSE3__) && defined(__AVX__)
c8_1 = _mm_add_sd(c8_1, _mm_mul_sd(a8_1, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_1 = _mm_add_sd(c8_1, _mm_mul_sd(a8_1, b8));
#endif
_mm_store_sd(&C[(i*12)+8], c8_1);
#else
C[(i*12)+2] += values[15] * B[(i*12)+8];
C[(i*12)+8] += values[16] * B[(i*12)+8];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b9 = _mm256_broadcast_sd(&B[(i*12)+9]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b9 = _mm_loaddup_pd(&B[(i*12)+9]);
#endif
__m128d c9_0 = _mm_load_sd(&C[(i*12)+0]);
__m128d a9_0 = _mm_load_sd(&values[17]);
#if defined(__SSE3__) && defined(__AVX__)
c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, b9));
#endif
_mm_store_sd(&C[(i*12)+0], c9_0);
__m128d c9_1 = _mm_load_sd(&C[(i*12)+3]);
__m128d a9_1 = _mm_load_sd(&values[18]);
#if defined(__SSE3__) && defined(__AVX__)
c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, b9));
#endif
_mm_store_sd(&C[(i*12)+3], c9_1);
__m128d c9_2 = _mm_load_sd(&C[(i*12)+9]);
__m128d a9_2 = _mm_load_sd(&values[19]);
#if defined(__SSE3__) && defined(__AVX__)
c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, b9));
#endif
_mm_store_sd(&C[(i*12)+9], c9_2);
#else
C[(i*12)+0] += values[17] * B[(i*12)+9];
C[(i*12)+3] += values[18] * B[(i*12)+9];
C[(i*12)+9] += values[19] * B[(i*12)+9];
#endif

}

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 360;
#endif

}

void dsparse_fP111DivM_m10_n9_k10_ldAna3_ldB12_ldC12_beta0_pfsigonly(const double* values, const double* B, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma nounroll_and_jam
for (unsigned int i = 0; i < 9; i++)
{
  #pragma simd
  for (unsigned int m = 0; m < 10; m++) {
    C[(i*12)+m] = 0.0;
  }
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b0 = _mm256_broadcast_sd(&B[(i*12)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b0 = _mm_loaddup_pd(&B[(i*12)+0]);
#endif
__m128d c0_0 = _mm_load_sd(&C[(i*12)+0]);
__m128d a0_0 = _mm_load_sd(&values[0]);
#if defined(__SSE3__) && defined(__AVX__)
c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, b0));
#endif
_mm_store_sd(&C[(i*12)+0], c0_0);
__m128d c0_1 = _mm_load_sd(&C[(i*12)+3]);
__m128d a0_1 = _mm_load_sd(&values[1]);
#if defined(__SSE3__) && defined(__AVX__)
c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, b0));
#endif
_mm_store_sd(&C[(i*12)+3], c0_1);
__m128d c0_2 = _mm_load_sd(&C[(i*12)+9]);
__m128d a0_2 = _mm_load_sd(&values[2]);
#if defined(__SSE3__) && defined(__AVX__)
c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, b0));
#endif
_mm_store_sd(&C[(i*12)+9], c0_2);
#else
C[(i*12)+0] += values[0] * B[(i*12)+0];
C[(i*12)+3] += values[1] * B[(i*12)+0];
C[(i*12)+9] += values[2] * B[(i*12)+0];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b1 = _mm256_broadcast_sd(&B[(i*12)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b1 = _mm_loaddup_pd(&B[(i*12)+1]);
#endif
__m128d c1_0 = _mm_loadu_pd(&C[(i*12)+1]);
__m128d a1_0 = _mm_loadu_pd(&values[3]);
#if defined(__SSE3__) && defined(__AVX__)
c1_0 = _mm_add_pd(c1_0, _mm_mul_pd(a1_0, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_0 = _mm_add_pd(c1_0, _mm_mul_pd(a1_0, b1));
#endif
_mm_storeu_pd(&C[(i*12)+1], c1_0);
__m128d c1_2 = _mm_loadu_pd(&C[(i*12)+7]);
__m128d a1_2 = _mm_loadu_pd(&values[5]);
#if defined(__SSE3__) && defined(__AVX__)
c1_2 = _mm_add_pd(c1_2, _mm_mul_pd(a1_2, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_2 = _mm_add_pd(c1_2, _mm_mul_pd(a1_2, b1));
#endif
_mm_storeu_pd(&C[(i*12)+7], c1_2);
#else
C[(i*12)+1] += values[3] * B[(i*12)+1];
C[(i*12)+2] += values[4] * B[(i*12)+1];
C[(i*12)+7] += values[5] * B[(i*12)+1];
C[(i*12)+8] += values[6] * B[(i*12)+1];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b2 = _mm256_broadcast_sd(&B[(i*12)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b2 = _mm_loaddup_pd(&B[(i*12)+2]);
#endif
__m128d c2_0 = _mm_loadu_pd(&C[(i*12)+1]);
__m128d a2_0 = _mm_loadu_pd(&values[7]);
#if defined(__SSE3__) && defined(__AVX__)
c2_0 = _mm_add_pd(c2_0, _mm_mul_pd(a2_0, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_0 = _mm_add_pd(c2_0, _mm_mul_pd(a2_0, b2));
#endif
_mm_storeu_pd(&C[(i*12)+1], c2_0);
__m128d c2_2 = _mm_loadu_pd(&C[(i*12)+7]);
__m128d a2_2 = _mm_loadu_pd(&values[9]);
#if defined(__SSE3__) && defined(__AVX__)
c2_2 = _mm_add_pd(c2_2, _mm_mul_pd(a2_2, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_2 = _mm_add_pd(c2_2, _mm_mul_pd(a2_2, b2));
#endif
_mm_storeu_pd(&C[(i*12)+7], c2_2);
#else
C[(i*12)+1] += values[7] * B[(i*12)+2];
C[(i*12)+2] += values[8] * B[(i*12)+2];
C[(i*12)+7] += values[9] * B[(i*12)+2];
C[(i*12)+8] += values[10] * B[(i*12)+2];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b3 = _mm256_broadcast_sd(&B[(i*12)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b3 = _mm_loaddup_pd(&B[(i*12)+3]);
#endif
__m128d c3_0 = _mm_load_sd(&C[(i*12)+0]);
__m128d a3_0 = _mm_load_sd(&values[11]);
#if defined(__SSE3__) && defined(__AVX__)
c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, b3));
#endif
_mm_store_sd(&C[(i*12)+0], c3_0);
__m128d c3_1 = _mm_load_sd(&C[(i*12)+3]);
__m128d a3_1 = _mm_load_sd(&values[12]);
#if defined(__SSE3__) && defined(__AVX__)
c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, b3));
#endif
_mm_store_sd(&C[(i*12)+3], c3_1);
__m128d c3_2 = _mm_load_sd(&C[(i*12)+9]);
__m128d a3_2 = _mm_load_sd(&values[13]);
#if defined(__SSE3__) && defined(__AVX__)
c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, b3));
#endif
_mm_store_sd(&C[(i*12)+9], c3_2);
#else
C[(i*12)+0] += values[11] * B[(i*12)+3];
C[(i*12)+3] += values[12] * B[(i*12)+3];
C[(i*12)+9] += values[13] * B[(i*12)+3];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b4 = _mm256_broadcast_sd(&B[(i*12)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b4 = _mm_loaddup_pd(&B[(i*12)+4]);
#endif
__m128d c4_0 = _mm_loadu_pd(&C[(i*12)+4]);
__m128d a4_0 = _mm_loadu_pd(&values[14]);
#if defined(__SSE3__) && defined(__AVX__)
c4_0 = _mm_add_pd(c4_0, _mm_mul_pd(a4_0, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_0 = _mm_add_pd(c4_0, _mm_mul_pd(a4_0, b4));
#endif
_mm_storeu_pd(&C[(i*12)+4], c4_0);
__m128d c4_2 = _mm_load_sd(&C[(i*12)+6]);
__m128d a4_2 = _mm_load_sd(&values[16]);
#if defined(__SSE3__) && defined(__AVX__)
c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, b4));
#endif
_mm_store_sd(&C[(i*12)+6], c4_2);
#else
C[(i*12)+4] += values[14] * B[(i*12)+4];
C[(i*12)+5] += values[15] * B[(i*12)+4];
C[(i*12)+6] += values[16] * B[(i*12)+4];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b5 = _mm256_broadcast_sd(&B[(i*12)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b5 = _mm_loaddup_pd(&B[(i*12)+5]);
#endif
__m128d c5_0 = _mm_loadu_pd(&C[(i*12)+4]);
__m128d a5_0 = _mm_loadu_pd(&values[17]);
#if defined(__SSE3__) && defined(__AVX__)
c5_0 = _mm_add_pd(c5_0, _mm_mul_pd(a5_0, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_0 = _mm_add_pd(c5_0, _mm_mul_pd(a5_0, b5));
#endif
_mm_storeu_pd(&C[(i*12)+4], c5_0);
__m128d c5_2 = _mm_load_sd(&C[(i*12)+6]);
__m128d a5_2 = _mm_load_sd(&values[19]);
#if defined(__SSE3__) && defined(__AVX__)
c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, b5));
#endif
_mm_store_sd(&C[(i*12)+6], c5_2);
#else
C[(i*12)+4] += values[17] * B[(i*12)+5];
C[(i*12)+5] += values[18] * B[(i*12)+5];
C[(i*12)+6] += values[19] * B[(i*12)+5];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b6 = _mm256_broadcast_sd(&B[(i*12)+6]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b6 = _mm_loaddup_pd(&B[(i*12)+6]);
#endif
__m128d c6_0 = _mm_loadu_pd(&C[(i*12)+4]);
__m128d a6_0 = _mm_loadu_pd(&values[20]);
#if defined(__SSE3__) && defined(__AVX__)
c6_0 = _mm_add_pd(c6_0, _mm_mul_pd(a6_0, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_0 = _mm_add_pd(c6_0, _mm_mul_pd(a6_0, b6));
#endif
_mm_storeu_pd(&C[(i*12)+4], c6_0);
__m128d c6_2 = _mm_load_sd(&C[(i*12)+6]);
__m128d a6_2 = _mm_load_sd(&values[22]);
#if defined(__SSE3__) && defined(__AVX__)
c6_2 = _mm_add_sd(c6_2, _mm_mul_sd(a6_2, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_2 = _mm_add_sd(c6_2, _mm_mul_sd(a6_2, b6));
#endif
_mm_store_sd(&C[(i*12)+6], c6_2);
#else
C[(i*12)+4] += values[20] * B[(i*12)+6];
C[(i*12)+5] += values[21] * B[(i*12)+6];
C[(i*12)+6] += values[22] * B[(i*12)+6];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b7 = _mm256_broadcast_sd(&B[(i*12)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b7 = _mm_loaddup_pd(&B[(i*12)+7]);
#endif
__m128d c7_0 = _mm_loadu_pd(&C[(i*12)+1]);
__m128d a7_0 = _mm_loadu_pd(&values[23]);
#if defined(__SSE3__) && defined(__AVX__)
c7_0 = _mm_add_pd(c7_0, _mm_mul_pd(a7_0, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_0 = _mm_add_pd(c7_0, _mm_mul_pd(a7_0, b7));
#endif
_mm_storeu_pd(&C[(i*12)+1], c7_0);
__m128d c7_2 = _mm_loadu_pd(&C[(i*12)+7]);
__m128d a7_2 = _mm_loadu_pd(&values[25]);
#if defined(__SSE3__) && defined(__AVX__)
c7_2 = _mm_add_pd(c7_2, _mm_mul_pd(a7_2, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_2 = _mm_add_pd(c7_2, _mm_mul_pd(a7_2, b7));
#endif
_mm_storeu_pd(&C[(i*12)+7], c7_2);
#else
C[(i*12)+1] += values[23] * B[(i*12)+7];
C[(i*12)+2] += values[24] * B[(i*12)+7];
C[(i*12)+7] += values[25] * B[(i*12)+7];
C[(i*12)+8] += values[26] * B[(i*12)+7];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b8 = _mm256_broadcast_sd(&B[(i*12)+8]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b8 = _mm_loaddup_pd(&B[(i*12)+8]);
#endif
__m128d c8_0 = _mm_loadu_pd(&C[(i*12)+1]);
__m128d a8_0 = _mm_loadu_pd(&values[27]);
#if defined(__SSE3__) && defined(__AVX__)
c8_0 = _mm_add_pd(c8_0, _mm_mul_pd(a8_0, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_0 = _mm_add_pd(c8_0, _mm_mul_pd(a8_0, b8));
#endif
_mm_storeu_pd(&C[(i*12)+1], c8_0);
__m128d c8_2 = _mm_loadu_pd(&C[(i*12)+7]);
__m128d a8_2 = _mm_loadu_pd(&values[29]);
#if defined(__SSE3__) && defined(__AVX__)
c8_2 = _mm_add_pd(c8_2, _mm_mul_pd(a8_2, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_2 = _mm_add_pd(c8_2, _mm_mul_pd(a8_2, b8));
#endif
_mm_storeu_pd(&C[(i*12)+7], c8_2);
#else
C[(i*12)+1] += values[27] * B[(i*12)+8];
C[(i*12)+2] += values[28] * B[(i*12)+8];
C[(i*12)+7] += values[29] * B[(i*12)+8];
C[(i*12)+8] += values[30] * B[(i*12)+8];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b9 = _mm256_broadcast_sd(&B[(i*12)+9]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b9 = _mm_loaddup_pd(&B[(i*12)+9]);
#endif
__m128d c9_0 = _mm_load_sd(&C[(i*12)+0]);
__m128d a9_0 = _mm_load_sd(&values[31]);
#if defined(__SSE3__) && defined(__AVX__)
c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, b9));
#endif
_mm_store_sd(&C[(i*12)+0], c9_0);
__m128d c9_1 = _mm_load_sd(&C[(i*12)+3]);
__m128d a9_1 = _mm_load_sd(&values[32]);
#if defined(__SSE3__) && defined(__AVX__)
c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, b9));
#endif
_mm_store_sd(&C[(i*12)+3], c9_1);
__m128d c9_2 = _mm_load_sd(&C[(i*12)+9]);
__m128d a9_2 = _mm_load_sd(&values[33]);
#if defined(__SSE3__) && defined(__AVX__)
c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, b9));
#endif
_mm_store_sd(&C[(i*12)+9], c9_2);
#else
C[(i*12)+0] += values[31] * B[(i*12)+9];
C[(i*12)+3] += values[32] * B[(i*12)+9];
C[(i*12)+9] += values[33] * B[(i*12)+9];
#endif

}

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 612;
#endif

}

void dsparse_starMatrix_m20_n9_k9_ldA20_ldBna4_ldC20_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
#pragma vector aligned
for (unsigned int i = 0; i < 20; i += 1)
{
C[(i)+(0)] += A[(i)+(120)] * values[0];
C[(i)+(0)] += A[(i)+(140)] * values[1];
C[(i)+(0)] += A[(i)+(160)] * values[2];
C[(i)+(20)] += A[(i)+(120)] * values[3];
C[(i)+(20)] += A[(i)+(140)] * values[4];
C[(i)+(20)] += A[(i)+(160)] * values[5];
C[(i)+(40)] += A[(i)+(120)] * values[6];
C[(i)+(40)] += A[(i)+(140)] * values[7];
C[(i)+(40)] += A[(i)+(160)] * values[8];
C[(i)+(60)] += A[(i)+(120)] * values[9];
C[(i)+(60)] += A[(i)+(140)] * values[10];
C[(i)+(80)] += A[(i)+(140)] * values[11];
C[(i)+(80)] += A[(i)+(160)] * values[12];
C[(i)+(100)] += A[(i)+(120)] * values[13];
C[(i)+(100)] += A[(i)+(160)] * values[14];
C[(i)+(120)] += A[(i)+(0)] * values[15];
C[(i)+(120)] += A[(i)+(60)] * values[16];
C[(i)+(120)] += A[(i)+(100)] * values[17];
C[(i)+(140)] += A[(i)+(20)] * values[18];
C[(i)+(140)] += A[(i)+(60)] * values[19];
C[(i)+(140)] += A[(i)+(80)] * values[20];
C[(i)+(160)] += A[(i)+(40)] * values[21];
C[(i)+(160)] += A[(i)+(80)] * values[22];
C[(i)+(160)] += A[(i)+(100)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 960;
#endif

}

void dsparse_starMatrix_m35_n9_k9_ldA36_ldBna5_ldC36_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
#pragma vector aligned
for (unsigned int i = 0; i < 35; i += 1)
{
C[(i)+(0)] += A[(i)+(216)] * values[0];
C[(i)+(0)] += A[(i)+(252)] * values[1];
C[(i)+(0)] += A[(i)+(288)] * values[2];
C[(i)+(36)] += A[(i)+(216)] * values[3];
C[(i)+(36)] += A[(i)+(252)] * values[4];
C[(i)+(36)] += A[(i)+(288)] * values[5];
C[(i)+(72)] += A[(i)+(216)] * values[6];
C[(i)+(72)] += A[(i)+(252)] * values[7];
C[(i)+(72)] += A[(i)+(288)] * values[8];
C[(i)+(108)] += A[(i)+(216)] * values[9];
C[(i)+(108)] += A[(i)+(252)] * values[10];
C[(i)+(144)] += A[(i)+(252)] * values[11];
C[(i)+(144)] += A[(i)+(288)] * values[12];
C[(i)+(180)] += A[(i)+(216)] * values[13];
C[(i)+(180)] += A[(i)+(288)] * values[14];
C[(i)+(216)] += A[(i)+(0)] * values[15];
C[(i)+(216)] += A[(i)+(108)] * values[16];
C[(i)+(216)] += A[(i)+(180)] * values[17];
C[(i)+(252)] += A[(i)+(36)] * values[18];
C[(i)+(252)] += A[(i)+(108)] * values[19];
C[(i)+(252)] += A[(i)+(144)] * values[20];
C[(i)+(288)] += A[(i)+(72)] * values[21];
C[(i)+(288)] += A[(i)+(144)] * values[22];
C[(i)+(288)] += A[(i)+(180)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 1680;
#endif

}

void dsparse_fM1DivM_m35_n9_k35_ldAna5_ldB36_ldC36_beta0_pfsigonly(const double* values, const double* B, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma nounroll_and_jam
for (unsigned int i = 0; i < 9; i++)
{
  #pragma simd
  for (unsigned int m = 0; m < 35; m++) {
    C[(i*36)+m] = 0.0;
  }
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b0 = _mm256_broadcast_sd(&B[(i*36)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b0 = _mm_loaddup_pd(&B[(i*36)+0]);
#endif
__m128d c0_0 = _mm_load_sd(&C[(i*36)+0]);
__m128d a0_0 = _mm_load_sd(&values[0]);
#if defined(__SSE3__) && defined(__AVX__)
c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, b0));
#endif
_mm_store_sd(&C[(i*36)+0], c0_0);
__m128d c0_1 = _mm_load_sd(&C[(i*36)+3]);
__m128d a0_1 = _mm_load_sd(&values[1]);
#if defined(__SSE3__) && defined(__AVX__)
c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, b0));
#endif
_mm_store_sd(&C[(i*36)+3], c0_1);
__m128d c0_2 = _mm_load_sd(&C[(i*36)+9]);
__m128d a0_2 = _mm_load_sd(&values[2]);
#if defined(__SSE3__) && defined(__AVX__)
c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, b0));
#endif
_mm_store_sd(&C[(i*36)+9], c0_2);
__m128d c0_3 = _mm_load_sd(&C[(i*36)+19]);
__m128d a0_3 = _mm_load_sd(&values[3]);
#if defined(__SSE3__) && defined(__AVX__)
c0_3 = _mm_add_sd(c0_3, _mm_mul_sd(a0_3, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_3 = _mm_add_sd(c0_3, _mm_mul_sd(a0_3, b0));
#endif
_mm_store_sd(&C[(i*36)+19], c0_3);
__m128d c0_4 = _mm_load_sd(&C[(i*36)+34]);
__m128d a0_4 = _mm_load_sd(&values[4]);
#if defined(__SSE3__) && defined(__AVX__)
c0_4 = _mm_add_sd(c0_4, _mm_mul_sd(a0_4, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_4 = _mm_add_sd(c0_4, _mm_mul_sd(a0_4, b0));
#endif
_mm_store_sd(&C[(i*36)+34], c0_4);
#else
C[(i*36)+0] += values[0] * B[(i*36)+0];
C[(i*36)+3] += values[1] * B[(i*36)+0];
C[(i*36)+9] += values[2] * B[(i*36)+0];
C[(i*36)+19] += values[3] * B[(i*36)+0];
C[(i*36)+34] += values[4] * B[(i*36)+0];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b1 = _mm256_broadcast_sd(&B[(i*36)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b1 = _mm_loaddup_pd(&B[(i*36)+1]);
#endif
__m128d c1_0 = _mm_load_sd(&C[(i*36)+1]);
__m128d a1_0 = _mm_load_sd(&values[5]);
#if defined(__SSE3__) && defined(__AVX__)
c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, b1));
#endif
_mm_store_sd(&C[(i*36)+1], c1_0);
__m128d c1_1 = _mm_load_sd(&C[(i*36)+7]);
__m128d a1_1 = _mm_load_sd(&values[6]);
#if defined(__SSE3__) && defined(__AVX__)
c1_1 = _mm_add_sd(c1_1, _mm_mul_sd(a1_1, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_1 = _mm_add_sd(c1_1, _mm_mul_sd(a1_1, b1));
#endif
_mm_store_sd(&C[(i*36)+7], c1_1);
__m128d c1_2 = _mm_load_sd(&C[(i*36)+17]);
__m128d a1_2 = _mm_load_sd(&values[7]);
#if defined(__SSE3__) && defined(__AVX__)
c1_2 = _mm_add_sd(c1_2, _mm_mul_sd(a1_2, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_2 = _mm_add_sd(c1_2, _mm_mul_sd(a1_2, b1));
#endif
_mm_store_sd(&C[(i*36)+17], c1_2);
__m128d c1_3 = _mm_load_sd(&C[(i*36)+32]);
__m128d a1_3 = _mm_load_sd(&values[8]);
#if defined(__SSE3__) && defined(__AVX__)
c1_3 = _mm_add_sd(c1_3, _mm_mul_sd(a1_3, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_3 = _mm_add_sd(c1_3, _mm_mul_sd(a1_3, b1));
#endif
_mm_store_sd(&C[(i*36)+32], c1_3);
#else
C[(i*36)+1] += values[5] * B[(i*36)+1];
C[(i*36)+7] += values[6] * B[(i*36)+1];
C[(i*36)+17] += values[7] * B[(i*36)+1];
C[(i*36)+32] += values[8] * B[(i*36)+1];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b2 = _mm256_broadcast_sd(&B[(i*36)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b2 = _mm_loaddup_pd(&B[(i*36)+2]);
#endif
__m128d c2_0 = _mm_load_sd(&C[(i*36)+2]);
__m128d a2_0 = _mm_load_sd(&values[9]);
#if defined(__SSE3__) && defined(__AVX__)
c2_0 = _mm_add_sd(c2_0, _mm_mul_sd(a2_0, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_0 = _mm_add_sd(c2_0, _mm_mul_sd(a2_0, b2));
#endif
_mm_store_sd(&C[(i*36)+2], c2_0);
__m128d c2_1 = _mm_load_sd(&C[(i*36)+8]);
__m128d a2_1 = _mm_load_sd(&values[10]);
#if defined(__SSE3__) && defined(__AVX__)
c2_1 = _mm_add_sd(c2_1, _mm_mul_sd(a2_1, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_1 = _mm_add_sd(c2_1, _mm_mul_sd(a2_1, b2));
#endif
_mm_store_sd(&C[(i*36)+8], c2_1);
__m128d c2_2 = _mm_load_sd(&C[(i*36)+18]);
__m128d a2_2 = _mm_load_sd(&values[11]);
#if defined(__SSE3__) && defined(__AVX__)
c2_2 = _mm_add_sd(c2_2, _mm_mul_sd(a2_2, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_2 = _mm_add_sd(c2_2, _mm_mul_sd(a2_2, b2));
#endif
_mm_store_sd(&C[(i*36)+18], c2_2);
__m128d c2_3 = _mm_load_sd(&C[(i*36)+33]);
__m128d a2_3 = _mm_load_sd(&values[12]);
#if defined(__SSE3__) && defined(__AVX__)
c2_3 = _mm_add_sd(c2_3, _mm_mul_sd(a2_3, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_3 = _mm_add_sd(c2_3, _mm_mul_sd(a2_3, b2));
#endif
_mm_store_sd(&C[(i*36)+33], c2_3);
#else
C[(i*36)+2] += values[9] * B[(i*36)+2];
C[(i*36)+8] += values[10] * B[(i*36)+2];
C[(i*36)+18] += values[11] * B[(i*36)+2];
C[(i*36)+33] += values[12] * B[(i*36)+2];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b3 = _mm256_broadcast_sd(&B[(i*36)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b3 = _mm_loaddup_pd(&B[(i*36)+3]);
#endif
__m128d c3_0 = _mm_load_sd(&C[(i*36)+0]);
__m128d a3_0 = _mm_load_sd(&values[13]);
#if defined(__SSE3__) && defined(__AVX__)
c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, b3));
#endif
_mm_store_sd(&C[(i*36)+0], c3_0);
__m128d c3_1 = _mm_load_sd(&C[(i*36)+3]);
__m128d a3_1 = _mm_load_sd(&values[14]);
#if defined(__SSE3__) && defined(__AVX__)
c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, b3));
#endif
_mm_store_sd(&C[(i*36)+3], c3_1);
__m128d c3_2 = _mm_load_sd(&C[(i*36)+9]);
__m128d a3_2 = _mm_load_sd(&values[15]);
#if defined(__SSE3__) && defined(__AVX__)
c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, b3));
#endif
_mm_store_sd(&C[(i*36)+9], c3_2);
__m128d c3_3 = _mm_load_sd(&C[(i*36)+19]);
__m128d a3_3 = _mm_load_sd(&values[16]);
#if defined(__SSE3__) && defined(__AVX__)
c3_3 = _mm_add_sd(c3_3, _mm_mul_sd(a3_3, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_3 = _mm_add_sd(c3_3, _mm_mul_sd(a3_3, b3));
#endif
_mm_store_sd(&C[(i*36)+19], c3_3);
__m128d c3_4 = _mm_load_sd(&C[(i*36)+34]);
__m128d a3_4 = _mm_load_sd(&values[17]);
#if defined(__SSE3__) && defined(__AVX__)
c3_4 = _mm_add_sd(c3_4, _mm_mul_sd(a3_4, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_4 = _mm_add_sd(c3_4, _mm_mul_sd(a3_4, b3));
#endif
_mm_store_sd(&C[(i*36)+34], c3_4);
#else
C[(i*36)+0] += values[13] * B[(i*36)+3];
C[(i*36)+3] += values[14] * B[(i*36)+3];
C[(i*36)+9] += values[15] * B[(i*36)+3];
C[(i*36)+19] += values[16] * B[(i*36)+3];
C[(i*36)+34] += values[17] * B[(i*36)+3];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b4 = _mm256_broadcast_sd(&B[(i*36)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b4 = _mm_loaddup_pd(&B[(i*36)+4]);
#endif
__m128d c4_0 = _mm_load_sd(&C[(i*36)+4]);
__m128d a4_0 = _mm_load_sd(&values[18]);
#if defined(__SSE3__) && defined(__AVX__)
c4_0 = _mm_add_sd(c4_0, _mm_mul_sd(a4_0, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_0 = _mm_add_sd(c4_0, _mm_mul_sd(a4_0, b4));
#endif
_mm_store_sd(&C[(i*36)+4], c4_0);
__m128d c4_1 = _mm_load_sd(&C[(i*36)+14]);
__m128d a4_1 = _mm_load_sd(&values[19]);
#if defined(__SSE3__) && defined(__AVX__)
c4_1 = _mm_add_sd(c4_1, _mm_mul_sd(a4_1, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_1 = _mm_add_sd(c4_1, _mm_mul_sd(a4_1, b4));
#endif
_mm_store_sd(&C[(i*36)+14], c4_1);
__m128d c4_2 = _mm_load_sd(&C[(i*36)+29]);
__m128d a4_2 = _mm_load_sd(&values[20]);
#if defined(__SSE3__) && defined(__AVX__)
c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, b4));
#endif
_mm_store_sd(&C[(i*36)+29], c4_2);
#else
C[(i*36)+4] += values[18] * B[(i*36)+4];
C[(i*36)+14] += values[19] * B[(i*36)+4];
C[(i*36)+29] += values[20] * B[(i*36)+4];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b5 = _mm256_broadcast_sd(&B[(i*36)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b5 = _mm_loaddup_pd(&B[(i*36)+5]);
#endif
__m128d c5_0 = _mm_load_sd(&C[(i*36)+5]);
__m128d a5_0 = _mm_load_sd(&values[21]);
#if defined(__SSE3__) && defined(__AVX__)
c5_0 = _mm_add_sd(c5_0, _mm_mul_sd(a5_0, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_0 = _mm_add_sd(c5_0, _mm_mul_sd(a5_0, b5));
#endif
_mm_store_sd(&C[(i*36)+5], c5_0);
__m128d c5_1 = _mm_load_sd(&C[(i*36)+15]);
__m128d a5_1 = _mm_load_sd(&values[22]);
#if defined(__SSE3__) && defined(__AVX__)
c5_1 = _mm_add_sd(c5_1, _mm_mul_sd(a5_1, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_1 = _mm_add_sd(c5_1, _mm_mul_sd(a5_1, b5));
#endif
_mm_store_sd(&C[(i*36)+15], c5_1);
__m128d c5_2 = _mm_load_sd(&C[(i*36)+30]);
__m128d a5_2 = _mm_load_sd(&values[23]);
#if defined(__SSE3__) && defined(__AVX__)
c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, b5));
#endif
_mm_store_sd(&C[(i*36)+30], c5_2);
#else
C[(i*36)+5] += values[21] * B[(i*36)+5];
C[(i*36)+15] += values[22] * B[(i*36)+5];
C[(i*36)+30] += values[23] * B[(i*36)+5];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b6 = _mm256_broadcast_sd(&B[(i*36)+6]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b6 = _mm_loaddup_pd(&B[(i*36)+6]);
#endif
__m128d c6_0 = _mm_load_sd(&C[(i*36)+6]);
__m128d a6_0 = _mm_load_sd(&values[24]);
#if defined(__SSE3__) && defined(__AVX__)
c6_0 = _mm_add_sd(c6_0, _mm_mul_sd(a6_0, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_0 = _mm_add_sd(c6_0, _mm_mul_sd(a6_0, b6));
#endif
_mm_store_sd(&C[(i*36)+6], c6_0);
__m128d c6_1 = _mm_load_sd(&C[(i*36)+16]);
__m128d a6_1 = _mm_load_sd(&values[25]);
#if defined(__SSE3__) && defined(__AVX__)
c6_1 = _mm_add_sd(c6_1, _mm_mul_sd(a6_1, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_1 = _mm_add_sd(c6_1, _mm_mul_sd(a6_1, b6));
#endif
_mm_store_sd(&C[(i*36)+16], c6_1);
__m128d c6_2 = _mm_load_sd(&C[(i*36)+31]);
__m128d a6_2 = _mm_load_sd(&values[26]);
#if defined(__SSE3__) && defined(__AVX__)
c6_2 = _mm_add_sd(c6_2, _mm_mul_sd(a6_2, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_2 = _mm_add_sd(c6_2, _mm_mul_sd(a6_2, b6));
#endif
_mm_store_sd(&C[(i*36)+31], c6_2);
#else
C[(i*36)+6] += values[24] * B[(i*36)+6];
C[(i*36)+16] += values[25] * B[(i*36)+6];
C[(i*36)+31] += values[26] * B[(i*36)+6];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b7 = _mm256_broadcast_sd(&B[(i*36)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b7 = _mm_loaddup_pd(&B[(i*36)+7]);
#endif
__m128d c7_0 = _mm_load_sd(&C[(i*36)+1]);
__m128d a7_0 = _mm_load_sd(&values[27]);
#if defined(__SSE3__) && defined(__AVX__)
c7_0 = _mm_add_sd(c7_0, _mm_mul_sd(a7_0, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_0 = _mm_add_sd(c7_0, _mm_mul_sd(a7_0, b7));
#endif
_mm_store_sd(&C[(i*36)+1], c7_0);
__m128d c7_1 = _mm_load_sd(&C[(i*36)+7]);
__m128d a7_1 = _mm_load_sd(&values[28]);
#if defined(__SSE3__) && defined(__AVX__)
c7_1 = _mm_add_sd(c7_1, _mm_mul_sd(a7_1, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_1 = _mm_add_sd(c7_1, _mm_mul_sd(a7_1, b7));
#endif
_mm_store_sd(&C[(i*36)+7], c7_1);
__m128d c7_2 = _mm_load_sd(&C[(i*36)+17]);
__m128d a7_2 = _mm_load_sd(&values[29]);
#if defined(__SSE3__) && defined(__AVX__)
c7_2 = _mm_add_sd(c7_2, _mm_mul_sd(a7_2, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_2 = _mm_add_sd(c7_2, _mm_mul_sd(a7_2, b7));
#endif
_mm_store_sd(&C[(i*36)+17], c7_2);
__m128d c7_3 = _mm_load_sd(&C[(i*36)+32]);
__m128d a7_3 = _mm_load_sd(&values[30]);
#if defined(__SSE3__) && defined(__AVX__)
c7_3 = _mm_add_sd(c7_3, _mm_mul_sd(a7_3, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_3 = _mm_add_sd(c7_3, _mm_mul_sd(a7_3, b7));
#endif
_mm_store_sd(&C[(i*36)+32], c7_3);
#else
C[(i*36)+1] += values[27] * B[(i*36)+7];
C[(i*36)+7] += values[28] * B[(i*36)+7];
C[(i*36)+17] += values[29] * B[(i*36)+7];
C[(i*36)+32] += values[30] * B[(i*36)+7];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b8 = _mm256_broadcast_sd(&B[(i*36)+8]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b8 = _mm_loaddup_pd(&B[(i*36)+8]);
#endif
__m128d c8_0 = _mm_load_sd(&C[(i*36)+2]);
__m128d a8_0 = _mm_load_sd(&values[31]);
#if defined(__SSE3__) && defined(__AVX__)
c8_0 = _mm_add_sd(c8_0, _mm_mul_sd(a8_0, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_0 = _mm_add_sd(c8_0, _mm_mul_sd(a8_0, b8));
#endif
_mm_store_sd(&C[(i*36)+2], c8_0);
__m128d c8_1 = _mm_load_sd(&C[(i*36)+8]);
__m128d a8_1 = _mm_load_sd(&values[32]);
#if defined(__SSE3__) && defined(__AVX__)
c8_1 = _mm_add_sd(c8_1, _mm_mul_sd(a8_1, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_1 = _mm_add_sd(c8_1, _mm_mul_sd(a8_1, b8));
#endif
_mm_store_sd(&C[(i*36)+8], c8_1);
__m128d c8_2 = _mm_load_sd(&C[(i*36)+18]);
__m128d a8_2 = _mm_load_sd(&values[33]);
#if defined(__SSE3__) && defined(__AVX__)
c8_2 = _mm_add_sd(c8_2, _mm_mul_sd(a8_2, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_2 = _mm_add_sd(c8_2, _mm_mul_sd(a8_2, b8));
#endif
_mm_store_sd(&C[(i*36)+18], c8_2);
__m128d c8_3 = _mm_load_sd(&C[(i*36)+33]);
__m128d a8_3 = _mm_load_sd(&values[34]);
#if defined(__SSE3__) && defined(__AVX__)
c8_3 = _mm_add_sd(c8_3, _mm_mul_sd(a8_3, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_3 = _mm_add_sd(c8_3, _mm_mul_sd(a8_3, b8));
#endif
_mm_store_sd(&C[(i*36)+33], c8_3);
#else
C[(i*36)+2] += values[31] * B[(i*36)+8];
C[(i*36)+8] += values[32] * B[(i*36)+8];
C[(i*36)+18] += values[33] * B[(i*36)+8];
C[(i*36)+33] += values[34] * B[(i*36)+8];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b9 = _mm256_broadcast_sd(&B[(i*36)+9]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b9 = _mm_loaddup_pd(&B[(i*36)+9]);
#endif
__m128d c9_0 = _mm_load_sd(&C[(i*36)+0]);
__m128d a9_0 = _mm_load_sd(&values[35]);
#if defined(__SSE3__) && defined(__AVX__)
c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, b9));
#endif
_mm_store_sd(&C[(i*36)+0], c9_0);
__m128d c9_1 = _mm_load_sd(&C[(i*36)+3]);
__m128d a9_1 = _mm_load_sd(&values[36]);
#if defined(__SSE3__) && defined(__AVX__)
c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, b9));
#endif
_mm_store_sd(&C[(i*36)+3], c9_1);
__m128d c9_2 = _mm_load_sd(&C[(i*36)+9]);
__m128d a9_2 = _mm_load_sd(&values[37]);
#if defined(__SSE3__) && defined(__AVX__)
c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, b9));
#endif
_mm_store_sd(&C[(i*36)+9], c9_2);
__m128d c9_3 = _mm_load_sd(&C[(i*36)+19]);
__m128d a9_3 = _mm_load_sd(&values[38]);
#if defined(__SSE3__) && defined(__AVX__)
c9_3 = _mm_add_sd(c9_3, _mm_mul_sd(a9_3, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_3 = _mm_add_sd(c9_3, _mm_mul_sd(a9_3, b9));
#endif
_mm_store_sd(&C[(i*36)+19], c9_3);
__m128d c9_4 = _mm_load_sd(&C[(i*36)+34]);
__m128d a9_4 = _mm_load_sd(&values[39]);
#if defined(__SSE3__) && defined(__AVX__)
c9_4 = _mm_add_sd(c9_4, _mm_mul_sd(a9_4, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_4 = _mm_add_sd(c9_4, _mm_mul_sd(a9_4, b9));
#endif
_mm_store_sd(&C[(i*36)+34], c9_4);
#else
C[(i*36)+0] += values[35] * B[(i*36)+9];
C[(i*36)+3] += values[36] * B[(i*36)+9];
C[(i*36)+9] += values[37] * B[(i*36)+9];
C[(i*36)+19] += values[38] * B[(i*36)+9];
C[(i*36)+34] += values[39] * B[(i*36)+9];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b10 = _mm256_broadcast_sd(&B[(i*36)+10]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b10 = _mm_loaddup_pd(&B[(i*36)+10]);
#endif
__m128d c10_0 = _mm_load_sd(&C[(i*36)+10]);
__m128d a10_0 = _mm_load_sd(&values[40]);
#if defined(__SSE3__) && defined(__AVX__)
c10_0 = _mm_add_sd(c10_0, _mm_mul_sd(a10_0, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c10_0 = _mm_add_sd(c10_0, _mm_mul_sd(a10_0, b10));
#endif
_mm_store_sd(&C[(i*36)+10], c10_0);
__m128d c10_1 = _mm_load_sd(&C[(i*36)+25]);
__m128d a10_1 = _mm_load_sd(&values[41]);
#if defined(__SSE3__) && defined(__AVX__)
c10_1 = _mm_add_sd(c10_1, _mm_mul_sd(a10_1, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c10_1 = _mm_add_sd(c10_1, _mm_mul_sd(a10_1, b10));
#endif
_mm_store_sd(&C[(i*36)+25], c10_1);
#else
C[(i*36)+10] += values[40] * B[(i*36)+10];
C[(i*36)+25] += values[41] * B[(i*36)+10];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b11 = _mm256_broadcast_sd(&B[(i*36)+11]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b11 = _mm_loaddup_pd(&B[(i*36)+11]);
#endif
__m128d c11_0 = _mm_load_sd(&C[(i*36)+11]);
__m128d a11_0 = _mm_load_sd(&values[42]);
#if defined(__SSE3__) && defined(__AVX__)
c11_0 = _mm_add_sd(c11_0, _mm_mul_sd(a11_0, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c11_0 = _mm_add_sd(c11_0, _mm_mul_sd(a11_0, b11));
#endif
_mm_store_sd(&C[(i*36)+11], c11_0);
__m128d c11_1 = _mm_load_sd(&C[(i*36)+26]);
__m128d a11_1 = _mm_load_sd(&values[43]);
#if defined(__SSE3__) && defined(__AVX__)
c11_1 = _mm_add_sd(c11_1, _mm_mul_sd(a11_1, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c11_1 = _mm_add_sd(c11_1, _mm_mul_sd(a11_1, b11));
#endif
_mm_store_sd(&C[(i*36)+26], c11_1);
#else
C[(i*36)+11] += values[42] * B[(i*36)+11];
C[(i*36)+26] += values[43] * B[(i*36)+11];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b12 = _mm256_broadcast_sd(&B[(i*36)+12]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b12 = _mm_loaddup_pd(&B[(i*36)+12]);
#endif
__m128d c12_0 = _mm_load_sd(&C[(i*36)+12]);
__m128d a12_0 = _mm_load_sd(&values[44]);
#if defined(__SSE3__) && defined(__AVX__)
c12_0 = _mm_add_sd(c12_0, _mm_mul_sd(a12_0, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c12_0 = _mm_add_sd(c12_0, _mm_mul_sd(a12_0, b12));
#endif
_mm_store_sd(&C[(i*36)+12], c12_0);
__m128d c12_1 = _mm_load_sd(&C[(i*36)+27]);
__m128d a12_1 = _mm_load_sd(&values[45]);
#if defined(__SSE3__) && defined(__AVX__)
c12_1 = _mm_add_sd(c12_1, _mm_mul_sd(a12_1, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c12_1 = _mm_add_sd(c12_1, _mm_mul_sd(a12_1, b12));
#endif
_mm_store_sd(&C[(i*36)+27], c12_1);
#else
C[(i*36)+12] += values[44] * B[(i*36)+12];
C[(i*36)+27] += values[45] * B[(i*36)+12];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b13 = _mm256_broadcast_sd(&B[(i*36)+13]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b13 = _mm_loaddup_pd(&B[(i*36)+13]);
#endif
__m128d c13_0 = _mm_load_sd(&C[(i*36)+13]);
__m128d a13_0 = _mm_load_sd(&values[46]);
#if defined(__SSE3__) && defined(__AVX__)
c13_0 = _mm_add_sd(c13_0, _mm_mul_sd(a13_0, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c13_0 = _mm_add_sd(c13_0, _mm_mul_sd(a13_0, b13));
#endif
_mm_store_sd(&C[(i*36)+13], c13_0);
__m128d c13_1 = _mm_load_sd(&C[(i*36)+28]);
__m128d a13_1 = _mm_load_sd(&values[47]);
#if defined(__SSE3__) && defined(__AVX__)
c13_1 = _mm_add_sd(c13_1, _mm_mul_sd(a13_1, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c13_1 = _mm_add_sd(c13_1, _mm_mul_sd(a13_1, b13));
#endif
_mm_store_sd(&C[(i*36)+28], c13_1);
#else
C[(i*36)+13] += values[46] * B[(i*36)+13];
C[(i*36)+28] += values[47] * B[(i*36)+13];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b14 = _mm256_broadcast_sd(&B[(i*36)+14]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b14 = _mm_loaddup_pd(&B[(i*36)+14]);
#endif
__m128d c14_0 = _mm_load_sd(&C[(i*36)+4]);
__m128d a14_0 = _mm_load_sd(&values[48]);
#if defined(__SSE3__) && defined(__AVX__)
c14_0 = _mm_add_sd(c14_0, _mm_mul_sd(a14_0, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_0 = _mm_add_sd(c14_0, _mm_mul_sd(a14_0, b14));
#endif
_mm_store_sd(&C[(i*36)+4], c14_0);
__m128d c14_1 = _mm_load_sd(&C[(i*36)+14]);
__m128d a14_1 = _mm_load_sd(&values[49]);
#if defined(__SSE3__) && defined(__AVX__)
c14_1 = _mm_add_sd(c14_1, _mm_mul_sd(a14_1, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_1 = _mm_add_sd(c14_1, _mm_mul_sd(a14_1, b14));
#endif
_mm_store_sd(&C[(i*36)+14], c14_1);
__m128d c14_2 = _mm_load_sd(&C[(i*36)+29]);
__m128d a14_2 = _mm_load_sd(&values[50]);
#if defined(__SSE3__) && defined(__AVX__)
c14_2 = _mm_add_sd(c14_2, _mm_mul_sd(a14_2, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_2 = _mm_add_sd(c14_2, _mm_mul_sd(a14_2, b14));
#endif
_mm_store_sd(&C[(i*36)+29], c14_2);
#else
C[(i*36)+4] += values[48] * B[(i*36)+14];
C[(i*36)+14] += values[49] * B[(i*36)+14];
C[(i*36)+29] += values[50] * B[(i*36)+14];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b15 = _mm256_broadcast_sd(&B[(i*36)+15]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b15 = _mm_loaddup_pd(&B[(i*36)+15]);
#endif
__m128d c15_0 = _mm_load_sd(&C[(i*36)+5]);
__m128d a15_0 = _mm_load_sd(&values[51]);
#if defined(__SSE3__) && defined(__AVX__)
c15_0 = _mm_add_sd(c15_0, _mm_mul_sd(a15_0, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_0 = _mm_add_sd(c15_0, _mm_mul_sd(a15_0, b15));
#endif
_mm_store_sd(&C[(i*36)+5], c15_0);
__m128d c15_1 = _mm_load_sd(&C[(i*36)+15]);
__m128d a15_1 = _mm_load_sd(&values[52]);
#if defined(__SSE3__) && defined(__AVX__)
c15_1 = _mm_add_sd(c15_1, _mm_mul_sd(a15_1, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_1 = _mm_add_sd(c15_1, _mm_mul_sd(a15_1, b15));
#endif
_mm_store_sd(&C[(i*36)+15], c15_1);
__m128d c15_2 = _mm_load_sd(&C[(i*36)+30]);
__m128d a15_2 = _mm_load_sd(&values[53]);
#if defined(__SSE3__) && defined(__AVX__)
c15_2 = _mm_add_sd(c15_2, _mm_mul_sd(a15_2, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_2 = _mm_add_sd(c15_2, _mm_mul_sd(a15_2, b15));
#endif
_mm_store_sd(&C[(i*36)+30], c15_2);
#else
C[(i*36)+5] += values[51] * B[(i*36)+15];
C[(i*36)+15] += values[52] * B[(i*36)+15];
C[(i*36)+30] += values[53] * B[(i*36)+15];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b16 = _mm256_broadcast_sd(&B[(i*36)+16]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b16 = _mm_loaddup_pd(&B[(i*36)+16]);
#endif
__m128d c16_0 = _mm_load_sd(&C[(i*36)+6]);
__m128d a16_0 = _mm_load_sd(&values[54]);
#if defined(__SSE3__) && defined(__AVX__)
c16_0 = _mm_add_sd(c16_0, _mm_mul_sd(a16_0, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_0 = _mm_add_sd(c16_0, _mm_mul_sd(a16_0, b16));
#endif
_mm_store_sd(&C[(i*36)+6], c16_0);
__m128d c16_1 = _mm_load_sd(&C[(i*36)+16]);
__m128d a16_1 = _mm_load_sd(&values[55]);
#if defined(__SSE3__) && defined(__AVX__)
c16_1 = _mm_add_sd(c16_1, _mm_mul_sd(a16_1, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_1 = _mm_add_sd(c16_1, _mm_mul_sd(a16_1, b16));
#endif
_mm_store_sd(&C[(i*36)+16], c16_1);
__m128d c16_2 = _mm_load_sd(&C[(i*36)+31]);
__m128d a16_2 = _mm_load_sd(&values[56]);
#if defined(__SSE3__) && defined(__AVX__)
c16_2 = _mm_add_sd(c16_2, _mm_mul_sd(a16_2, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_2 = _mm_add_sd(c16_2, _mm_mul_sd(a16_2, b16));
#endif
_mm_store_sd(&C[(i*36)+31], c16_2);
#else
C[(i*36)+6] += values[54] * B[(i*36)+16];
C[(i*36)+16] += values[55] * B[(i*36)+16];
C[(i*36)+31] += values[56] * B[(i*36)+16];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b17 = _mm256_broadcast_sd(&B[(i*36)+17]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b17 = _mm_loaddup_pd(&B[(i*36)+17]);
#endif
__m128d c17_0 = _mm_load_sd(&C[(i*36)+1]);
__m128d a17_0 = _mm_load_sd(&values[57]);
#if defined(__SSE3__) && defined(__AVX__)
c17_0 = _mm_add_sd(c17_0, _mm_mul_sd(a17_0, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_0 = _mm_add_sd(c17_0, _mm_mul_sd(a17_0, b17));
#endif
_mm_store_sd(&C[(i*36)+1], c17_0);
__m128d c17_1 = _mm_load_sd(&C[(i*36)+7]);
__m128d a17_1 = _mm_load_sd(&values[58]);
#if defined(__SSE3__) && defined(__AVX__)
c17_1 = _mm_add_sd(c17_1, _mm_mul_sd(a17_1, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_1 = _mm_add_sd(c17_1, _mm_mul_sd(a17_1, b17));
#endif
_mm_store_sd(&C[(i*36)+7], c17_1);
__m128d c17_2 = _mm_load_sd(&C[(i*36)+17]);
__m128d a17_2 = _mm_load_sd(&values[59]);
#if defined(__SSE3__) && defined(__AVX__)
c17_2 = _mm_add_sd(c17_2, _mm_mul_sd(a17_2, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_2 = _mm_add_sd(c17_2, _mm_mul_sd(a17_2, b17));
#endif
_mm_store_sd(&C[(i*36)+17], c17_2);
__m128d c17_3 = _mm_load_sd(&C[(i*36)+32]);
__m128d a17_3 = _mm_load_sd(&values[60]);
#if defined(__SSE3__) && defined(__AVX__)
c17_3 = _mm_add_sd(c17_3, _mm_mul_sd(a17_3, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_3 = _mm_add_sd(c17_3, _mm_mul_sd(a17_3, b17));
#endif
_mm_store_sd(&C[(i*36)+32], c17_3);
#else
C[(i*36)+1] += values[57] * B[(i*36)+17];
C[(i*36)+7] += values[58] * B[(i*36)+17];
C[(i*36)+17] += values[59] * B[(i*36)+17];
C[(i*36)+32] += values[60] * B[(i*36)+17];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b18 = _mm256_broadcast_sd(&B[(i*36)+18]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b18 = _mm_loaddup_pd(&B[(i*36)+18]);
#endif
__m128d c18_0 = _mm_load_sd(&C[(i*36)+2]);
__m128d a18_0 = _mm_load_sd(&values[61]);
#if defined(__SSE3__) && defined(__AVX__)
c18_0 = _mm_add_sd(c18_0, _mm_mul_sd(a18_0, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_0 = _mm_add_sd(c18_0, _mm_mul_sd(a18_0, b18));
#endif
_mm_store_sd(&C[(i*36)+2], c18_0);
__m128d c18_1 = _mm_load_sd(&C[(i*36)+8]);
__m128d a18_1 = _mm_load_sd(&values[62]);
#if defined(__SSE3__) && defined(__AVX__)
c18_1 = _mm_add_sd(c18_1, _mm_mul_sd(a18_1, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_1 = _mm_add_sd(c18_1, _mm_mul_sd(a18_1, b18));
#endif
_mm_store_sd(&C[(i*36)+8], c18_1);
__m128d c18_2 = _mm_load_sd(&C[(i*36)+18]);
__m128d a18_2 = _mm_load_sd(&values[63]);
#if defined(__SSE3__) && defined(__AVX__)
c18_2 = _mm_add_sd(c18_2, _mm_mul_sd(a18_2, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_2 = _mm_add_sd(c18_2, _mm_mul_sd(a18_2, b18));
#endif
_mm_store_sd(&C[(i*36)+18], c18_2);
__m128d c18_3 = _mm_load_sd(&C[(i*36)+33]);
__m128d a18_3 = _mm_load_sd(&values[64]);
#if defined(__SSE3__) && defined(__AVX__)
c18_3 = _mm_add_sd(c18_3, _mm_mul_sd(a18_3, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_3 = _mm_add_sd(c18_3, _mm_mul_sd(a18_3, b18));
#endif
_mm_store_sd(&C[(i*36)+33], c18_3);
#else
C[(i*36)+2] += values[61] * B[(i*36)+18];
C[(i*36)+8] += values[62] * B[(i*36)+18];
C[(i*36)+18] += values[63] * B[(i*36)+18];
C[(i*36)+33] += values[64] * B[(i*36)+18];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b19 = _mm256_broadcast_sd(&B[(i*36)+19]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b19 = _mm_loaddup_pd(&B[(i*36)+19]);
#endif
__m128d c19_0 = _mm_load_sd(&C[(i*36)+0]);
__m128d a19_0 = _mm_load_sd(&values[65]);
#if defined(__SSE3__) && defined(__AVX__)
c19_0 = _mm_add_sd(c19_0, _mm_mul_sd(a19_0, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_0 = _mm_add_sd(c19_0, _mm_mul_sd(a19_0, b19));
#endif
_mm_store_sd(&C[(i*36)+0], c19_0);
__m128d c19_1 = _mm_load_sd(&C[(i*36)+3]);
__m128d a19_1 = _mm_load_sd(&values[66]);
#if defined(__SSE3__) && defined(__AVX__)
c19_1 = _mm_add_sd(c19_1, _mm_mul_sd(a19_1, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_1 = _mm_add_sd(c19_1, _mm_mul_sd(a19_1, b19));
#endif
_mm_store_sd(&C[(i*36)+3], c19_1);
__m128d c19_2 = _mm_load_sd(&C[(i*36)+9]);
__m128d a19_2 = _mm_load_sd(&values[67]);
#if defined(__SSE3__) && defined(__AVX__)
c19_2 = _mm_add_sd(c19_2, _mm_mul_sd(a19_2, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_2 = _mm_add_sd(c19_2, _mm_mul_sd(a19_2, b19));
#endif
_mm_store_sd(&C[(i*36)+9], c19_2);
__m128d c19_3 = _mm_load_sd(&C[(i*36)+19]);
__m128d a19_3 = _mm_load_sd(&values[68]);
#if defined(__SSE3__) && defined(__AVX__)
c19_3 = _mm_add_sd(c19_3, _mm_mul_sd(a19_3, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_3 = _mm_add_sd(c19_3, _mm_mul_sd(a19_3, b19));
#endif
_mm_store_sd(&C[(i*36)+19], c19_3);
__m128d c19_4 = _mm_load_sd(&C[(i*36)+34]);
__m128d a19_4 = _mm_load_sd(&values[69]);
#if defined(__SSE3__) && defined(__AVX__)
c19_4 = _mm_add_sd(c19_4, _mm_mul_sd(a19_4, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_4 = _mm_add_sd(c19_4, _mm_mul_sd(a19_4, b19));
#endif
_mm_store_sd(&C[(i*36)+34], c19_4);
#else
C[(i*36)+0] += values[65] * B[(i*36)+19];
C[(i*36)+3] += values[66] * B[(i*36)+19];
C[(i*36)+9] += values[67] * B[(i*36)+19];
C[(i*36)+19] += values[68] * B[(i*36)+19];
C[(i*36)+34] += values[69] * B[(i*36)+19];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b20 = _mm256_broadcast_sd(&B[(i*36)+20]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b20 = _mm_loaddup_pd(&B[(i*36)+20]);
#endif
__m128d c20_0 = _mm_load_sd(&C[(i*36)+20]);
__m128d a20_0 = _mm_load_sd(&values[70]);
#if defined(__SSE3__) && defined(__AVX__)
c20_0 = _mm_add_sd(c20_0, _mm_mul_sd(a20_0, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c20_0 = _mm_add_sd(c20_0, _mm_mul_sd(a20_0, b20));
#endif
_mm_store_sd(&C[(i*36)+20], c20_0);
#else
C[(i*36)+20] += values[70] * B[(i*36)+20];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b21 = _mm256_broadcast_sd(&B[(i*36)+21]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b21 = _mm_loaddup_pd(&B[(i*36)+21]);
#endif
__m128d c21_0 = _mm_load_sd(&C[(i*36)+21]);
__m128d a21_0 = _mm_load_sd(&values[71]);
#if defined(__SSE3__) && defined(__AVX__)
c21_0 = _mm_add_sd(c21_0, _mm_mul_sd(a21_0, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c21_0 = _mm_add_sd(c21_0, _mm_mul_sd(a21_0, b21));
#endif
_mm_store_sd(&C[(i*36)+21], c21_0);
#else
C[(i*36)+21] += values[71] * B[(i*36)+21];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b22 = _mm256_broadcast_sd(&B[(i*36)+22]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b22 = _mm_loaddup_pd(&B[(i*36)+22]);
#endif
__m128d c22_0 = _mm_load_sd(&C[(i*36)+22]);
__m128d a22_0 = _mm_load_sd(&values[72]);
#if defined(__SSE3__) && defined(__AVX__)
c22_0 = _mm_add_sd(c22_0, _mm_mul_sd(a22_0, _mm256_castpd256_pd128(b22)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c22_0 = _mm_add_sd(c22_0, _mm_mul_sd(a22_0, b22));
#endif
_mm_store_sd(&C[(i*36)+22], c22_0);
#else
C[(i*36)+22] += values[72] * B[(i*36)+22];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b23 = _mm256_broadcast_sd(&B[(i*36)+23]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b23 = _mm_loaddup_pd(&B[(i*36)+23]);
#endif
__m128d c23_0 = _mm_load_sd(&C[(i*36)+23]);
__m128d a23_0 = _mm_load_sd(&values[73]);
#if defined(__SSE3__) && defined(__AVX__)
c23_0 = _mm_add_sd(c23_0, _mm_mul_sd(a23_0, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c23_0 = _mm_add_sd(c23_0, _mm_mul_sd(a23_0, b23));
#endif
_mm_store_sd(&C[(i*36)+23], c23_0);
#else
C[(i*36)+23] += values[73] * B[(i*36)+23];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b24 = _mm256_broadcast_sd(&B[(i*36)+24]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b24 = _mm_loaddup_pd(&B[(i*36)+24]);
#endif
__m128d c24_0 = _mm_load_sd(&C[(i*36)+24]);
__m128d a24_0 = _mm_load_sd(&values[74]);
#if defined(__SSE3__) && defined(__AVX__)
c24_0 = _mm_add_sd(c24_0, _mm_mul_sd(a24_0, _mm256_castpd256_pd128(b24)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c24_0 = _mm_add_sd(c24_0, _mm_mul_sd(a24_0, b24));
#endif
_mm_store_sd(&C[(i*36)+24], c24_0);
#else
C[(i*36)+24] += values[74] * B[(i*36)+24];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b25 = _mm256_broadcast_sd(&B[(i*36)+25]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b25 = _mm_loaddup_pd(&B[(i*36)+25]);
#endif
__m128d c25_0 = _mm_load_sd(&C[(i*36)+10]);
__m128d a25_0 = _mm_load_sd(&values[75]);
#if defined(__SSE3__) && defined(__AVX__)
c25_0 = _mm_add_sd(c25_0, _mm_mul_sd(a25_0, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c25_0 = _mm_add_sd(c25_0, _mm_mul_sd(a25_0, b25));
#endif
_mm_store_sd(&C[(i*36)+10], c25_0);
__m128d c25_1 = _mm_load_sd(&C[(i*36)+25]);
__m128d a25_1 = _mm_load_sd(&values[76]);
#if defined(__SSE3__) && defined(__AVX__)
c25_1 = _mm_add_sd(c25_1, _mm_mul_sd(a25_1, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c25_1 = _mm_add_sd(c25_1, _mm_mul_sd(a25_1, b25));
#endif
_mm_store_sd(&C[(i*36)+25], c25_1);
#else
C[(i*36)+10] += values[75] * B[(i*36)+25];
C[(i*36)+25] += values[76] * B[(i*36)+25];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b26 = _mm256_broadcast_sd(&B[(i*36)+26]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b26 = _mm_loaddup_pd(&B[(i*36)+26]);
#endif
__m128d c26_0 = _mm_load_sd(&C[(i*36)+11]);
__m128d a26_0 = _mm_load_sd(&values[77]);
#if defined(__SSE3__) && defined(__AVX__)
c26_0 = _mm_add_sd(c26_0, _mm_mul_sd(a26_0, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c26_0 = _mm_add_sd(c26_0, _mm_mul_sd(a26_0, b26));
#endif
_mm_store_sd(&C[(i*36)+11], c26_0);
__m128d c26_1 = _mm_load_sd(&C[(i*36)+26]);
__m128d a26_1 = _mm_load_sd(&values[78]);
#if defined(__SSE3__) && defined(__AVX__)
c26_1 = _mm_add_sd(c26_1, _mm_mul_sd(a26_1, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c26_1 = _mm_add_sd(c26_1, _mm_mul_sd(a26_1, b26));
#endif
_mm_store_sd(&C[(i*36)+26], c26_1);
#else
C[(i*36)+11] += values[77] * B[(i*36)+26];
C[(i*36)+26] += values[78] * B[(i*36)+26];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b27 = _mm256_broadcast_sd(&B[(i*36)+27]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b27 = _mm_loaddup_pd(&B[(i*36)+27]);
#endif
__m128d c27_0 = _mm_load_sd(&C[(i*36)+12]);
__m128d a27_0 = _mm_load_sd(&values[79]);
#if defined(__SSE3__) && defined(__AVX__)
c27_0 = _mm_add_sd(c27_0, _mm_mul_sd(a27_0, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c27_0 = _mm_add_sd(c27_0, _mm_mul_sd(a27_0, b27));
#endif
_mm_store_sd(&C[(i*36)+12], c27_0);
__m128d c27_1 = _mm_load_sd(&C[(i*36)+27]);
__m128d a27_1 = _mm_load_sd(&values[80]);
#if defined(__SSE3__) && defined(__AVX__)
c27_1 = _mm_add_sd(c27_1, _mm_mul_sd(a27_1, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c27_1 = _mm_add_sd(c27_1, _mm_mul_sd(a27_1, b27));
#endif
_mm_store_sd(&C[(i*36)+27], c27_1);
#else
C[(i*36)+12] += values[79] * B[(i*36)+27];
C[(i*36)+27] += values[80] * B[(i*36)+27];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b28 = _mm256_broadcast_sd(&B[(i*36)+28]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b28 = _mm_loaddup_pd(&B[(i*36)+28]);
#endif
__m128d c28_0 = _mm_load_sd(&C[(i*36)+13]);
__m128d a28_0 = _mm_load_sd(&values[81]);
#if defined(__SSE3__) && defined(__AVX__)
c28_0 = _mm_add_sd(c28_0, _mm_mul_sd(a28_0, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c28_0 = _mm_add_sd(c28_0, _mm_mul_sd(a28_0, b28));
#endif
_mm_store_sd(&C[(i*36)+13], c28_0);
__m128d c28_1 = _mm_load_sd(&C[(i*36)+28]);
__m128d a28_1 = _mm_load_sd(&values[82]);
#if defined(__SSE3__) && defined(__AVX__)
c28_1 = _mm_add_sd(c28_1, _mm_mul_sd(a28_1, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c28_1 = _mm_add_sd(c28_1, _mm_mul_sd(a28_1, b28));
#endif
_mm_store_sd(&C[(i*36)+28], c28_1);
#else
C[(i*36)+13] += values[81] * B[(i*36)+28];
C[(i*36)+28] += values[82] * B[(i*36)+28];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b29 = _mm256_broadcast_sd(&B[(i*36)+29]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b29 = _mm_loaddup_pd(&B[(i*36)+29]);
#endif
__m128d c29_0 = _mm_load_sd(&C[(i*36)+4]);
__m128d a29_0 = _mm_load_sd(&values[83]);
#if defined(__SSE3__) && defined(__AVX__)
c29_0 = _mm_add_sd(c29_0, _mm_mul_sd(a29_0, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_0 = _mm_add_sd(c29_0, _mm_mul_sd(a29_0, b29));
#endif
_mm_store_sd(&C[(i*36)+4], c29_0);
__m128d c29_1 = _mm_load_sd(&C[(i*36)+14]);
__m128d a29_1 = _mm_load_sd(&values[84]);
#if defined(__SSE3__) && defined(__AVX__)
c29_1 = _mm_add_sd(c29_1, _mm_mul_sd(a29_1, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_1 = _mm_add_sd(c29_1, _mm_mul_sd(a29_1, b29));
#endif
_mm_store_sd(&C[(i*36)+14], c29_1);
__m128d c29_2 = _mm_load_sd(&C[(i*36)+29]);
__m128d a29_2 = _mm_load_sd(&values[85]);
#if defined(__SSE3__) && defined(__AVX__)
c29_2 = _mm_add_sd(c29_2, _mm_mul_sd(a29_2, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_2 = _mm_add_sd(c29_2, _mm_mul_sd(a29_2, b29));
#endif
_mm_store_sd(&C[(i*36)+29], c29_2);
#else
C[(i*36)+4] += values[83] * B[(i*36)+29];
C[(i*36)+14] += values[84] * B[(i*36)+29];
C[(i*36)+29] += values[85] * B[(i*36)+29];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b30 = _mm256_broadcast_sd(&B[(i*36)+30]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b30 = _mm_loaddup_pd(&B[(i*36)+30]);
#endif
__m128d c30_0 = _mm_load_sd(&C[(i*36)+5]);
__m128d a30_0 = _mm_load_sd(&values[86]);
#if defined(__SSE3__) && defined(__AVX__)
c30_0 = _mm_add_sd(c30_0, _mm_mul_sd(a30_0, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_0 = _mm_add_sd(c30_0, _mm_mul_sd(a30_0, b30));
#endif
_mm_store_sd(&C[(i*36)+5], c30_0);
__m128d c30_1 = _mm_load_sd(&C[(i*36)+15]);
__m128d a30_1 = _mm_load_sd(&values[87]);
#if defined(__SSE3__) && defined(__AVX__)
c30_1 = _mm_add_sd(c30_1, _mm_mul_sd(a30_1, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_1 = _mm_add_sd(c30_1, _mm_mul_sd(a30_1, b30));
#endif
_mm_store_sd(&C[(i*36)+15], c30_1);
__m128d c30_2 = _mm_load_sd(&C[(i*36)+30]);
__m128d a30_2 = _mm_load_sd(&values[88]);
#if defined(__SSE3__) && defined(__AVX__)
c30_2 = _mm_add_sd(c30_2, _mm_mul_sd(a30_2, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_2 = _mm_add_sd(c30_2, _mm_mul_sd(a30_2, b30));
#endif
_mm_store_sd(&C[(i*36)+30], c30_2);
#else
C[(i*36)+5] += values[86] * B[(i*36)+30];
C[(i*36)+15] += values[87] * B[(i*36)+30];
C[(i*36)+30] += values[88] * B[(i*36)+30];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b31 = _mm256_broadcast_sd(&B[(i*36)+31]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b31 = _mm_loaddup_pd(&B[(i*36)+31]);
#endif
__m128d c31_0 = _mm_load_sd(&C[(i*36)+6]);
__m128d a31_0 = _mm_load_sd(&values[89]);
#if defined(__SSE3__) && defined(__AVX__)
c31_0 = _mm_add_sd(c31_0, _mm_mul_sd(a31_0, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_0 = _mm_add_sd(c31_0, _mm_mul_sd(a31_0, b31));
#endif
_mm_store_sd(&C[(i*36)+6], c31_0);
__m128d c31_1 = _mm_load_sd(&C[(i*36)+16]);
__m128d a31_1 = _mm_load_sd(&values[90]);
#if defined(__SSE3__) && defined(__AVX__)
c31_1 = _mm_add_sd(c31_1, _mm_mul_sd(a31_1, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_1 = _mm_add_sd(c31_1, _mm_mul_sd(a31_1, b31));
#endif
_mm_store_sd(&C[(i*36)+16], c31_1);
__m128d c31_2 = _mm_load_sd(&C[(i*36)+31]);
__m128d a31_2 = _mm_load_sd(&values[91]);
#if defined(__SSE3__) && defined(__AVX__)
c31_2 = _mm_add_sd(c31_2, _mm_mul_sd(a31_2, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_2 = _mm_add_sd(c31_2, _mm_mul_sd(a31_2, b31));
#endif
_mm_store_sd(&C[(i*36)+31], c31_2);
#else
C[(i*36)+6] += values[89] * B[(i*36)+31];
C[(i*36)+16] += values[90] * B[(i*36)+31];
C[(i*36)+31] += values[91] * B[(i*36)+31];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b32 = _mm256_broadcast_sd(&B[(i*36)+32]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b32 = _mm_loaddup_pd(&B[(i*36)+32]);
#endif
__m128d c32_0 = _mm_load_sd(&C[(i*36)+1]);
__m128d a32_0 = _mm_load_sd(&values[92]);
#if defined(__SSE3__) && defined(__AVX__)
c32_0 = _mm_add_sd(c32_0, _mm_mul_sd(a32_0, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_0 = _mm_add_sd(c32_0, _mm_mul_sd(a32_0, b32));
#endif
_mm_store_sd(&C[(i*36)+1], c32_0);
__m128d c32_1 = _mm_load_sd(&C[(i*36)+7]);
__m128d a32_1 = _mm_load_sd(&values[93]);
#if defined(__SSE3__) && defined(__AVX__)
c32_1 = _mm_add_sd(c32_1, _mm_mul_sd(a32_1, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_1 = _mm_add_sd(c32_1, _mm_mul_sd(a32_1, b32));
#endif
_mm_store_sd(&C[(i*36)+7], c32_1);
__m128d c32_2 = _mm_load_sd(&C[(i*36)+17]);
__m128d a32_2 = _mm_load_sd(&values[94]);
#if defined(__SSE3__) && defined(__AVX__)
c32_2 = _mm_add_sd(c32_2, _mm_mul_sd(a32_2, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_2 = _mm_add_sd(c32_2, _mm_mul_sd(a32_2, b32));
#endif
_mm_store_sd(&C[(i*36)+17], c32_2);
__m128d c32_3 = _mm_load_sd(&C[(i*36)+32]);
__m128d a32_3 = _mm_load_sd(&values[95]);
#if defined(__SSE3__) && defined(__AVX__)
c32_3 = _mm_add_sd(c32_3, _mm_mul_sd(a32_3, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_3 = _mm_add_sd(c32_3, _mm_mul_sd(a32_3, b32));
#endif
_mm_store_sd(&C[(i*36)+32], c32_3);
#else
C[(i*36)+1] += values[92] * B[(i*36)+32];
C[(i*36)+7] += values[93] * B[(i*36)+32];
C[(i*36)+17] += values[94] * B[(i*36)+32];
C[(i*36)+32] += values[95] * B[(i*36)+32];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b33 = _mm256_broadcast_sd(&B[(i*36)+33]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b33 = _mm_loaddup_pd(&B[(i*36)+33]);
#endif
__m128d c33_0 = _mm_load_sd(&C[(i*36)+2]);
__m128d a33_0 = _mm_load_sd(&values[96]);
#if defined(__SSE3__) && defined(__AVX__)
c33_0 = _mm_add_sd(c33_0, _mm_mul_sd(a33_0, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_0 = _mm_add_sd(c33_0, _mm_mul_sd(a33_0, b33));
#endif
_mm_store_sd(&C[(i*36)+2], c33_0);
__m128d c33_1 = _mm_load_sd(&C[(i*36)+8]);
__m128d a33_1 = _mm_load_sd(&values[97]);
#if defined(__SSE3__) && defined(__AVX__)
c33_1 = _mm_add_sd(c33_1, _mm_mul_sd(a33_1, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_1 = _mm_add_sd(c33_1, _mm_mul_sd(a33_1, b33));
#endif
_mm_store_sd(&C[(i*36)+8], c33_1);
__m128d c33_2 = _mm_load_sd(&C[(i*36)+18]);
__m128d a33_2 = _mm_load_sd(&values[98]);
#if defined(__SSE3__) && defined(__AVX__)
c33_2 = _mm_add_sd(c33_2, _mm_mul_sd(a33_2, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_2 = _mm_add_sd(c33_2, _mm_mul_sd(a33_2, b33));
#endif
_mm_store_sd(&C[(i*36)+18], c33_2);
__m128d c33_3 = _mm_load_sd(&C[(i*36)+33]);
__m128d a33_3 = _mm_load_sd(&values[99]);
#if defined(__SSE3__) && defined(__AVX__)
c33_3 = _mm_add_sd(c33_3, _mm_mul_sd(a33_3, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_3 = _mm_add_sd(c33_3, _mm_mul_sd(a33_3, b33));
#endif
_mm_store_sd(&C[(i*36)+33], c33_3);
#else
C[(i*36)+2] += values[96] * B[(i*36)+33];
C[(i*36)+8] += values[97] * B[(i*36)+33];
C[(i*36)+18] += values[98] * B[(i*36)+33];
C[(i*36)+33] += values[99] * B[(i*36)+33];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b34 = _mm256_broadcast_sd(&B[(i*36)+34]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b34 = _mm_loaddup_pd(&B[(i*36)+34]);
#endif
__m128d c34_0 = _mm_load_sd(&C[(i*36)+0]);
__m128d a34_0 = _mm_load_sd(&values[100]);
#if defined(__SSE3__) && defined(__AVX__)
c34_0 = _mm_add_sd(c34_0, _mm_mul_sd(a34_0, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_0 = _mm_add_sd(c34_0, _mm_mul_sd(a34_0, b34));
#endif
_mm_store_sd(&C[(i*36)+0], c34_0);
__m128d c34_1 = _mm_load_sd(&C[(i*36)+3]);
__m128d a34_1 = _mm_load_sd(&values[101]);
#if defined(__SSE3__) && defined(__AVX__)
c34_1 = _mm_add_sd(c34_1, _mm_mul_sd(a34_1, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_1 = _mm_add_sd(c34_1, _mm_mul_sd(a34_1, b34));
#endif
_mm_store_sd(&C[(i*36)+3], c34_1);
__m128d c34_2 = _mm_load_sd(&C[(i*36)+9]);
__m128d a34_2 = _mm_load_sd(&values[102]);
#if defined(__SSE3__) && defined(__AVX__)
c34_2 = _mm_add_sd(c34_2, _mm_mul_sd(a34_2, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_2 = _mm_add_sd(c34_2, _mm_mul_sd(a34_2, b34));
#endif
_mm_store_sd(&C[(i*36)+9], c34_2);
__m128d c34_3 = _mm_load_sd(&C[(i*36)+19]);
__m128d a34_3 = _mm_load_sd(&values[103]);
#if defined(__SSE3__) && defined(__AVX__)
c34_3 = _mm_add_sd(c34_3, _mm_mul_sd(a34_3, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_3 = _mm_add_sd(c34_3, _mm_mul_sd(a34_3, b34));
#endif
_mm_store_sd(&C[(i*36)+19], c34_3);
__m128d c34_4 = _mm_load_sd(&C[(i*36)+34]);
__m128d a34_4 = _mm_load_sd(&values[104]);
#if defined(__SSE3__) && defined(__AVX__)
c34_4 = _mm_add_sd(c34_4, _mm_mul_sd(a34_4, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_4 = _mm_add_sd(c34_4, _mm_mul_sd(a34_4, b34));
#endif
_mm_store_sd(&C[(i*36)+34], c34_4);
#else
C[(i*36)+0] += values[100] * B[(i*36)+34];
C[(i*36)+3] += values[101] * B[(i*36)+34];
C[(i*36)+9] += values[102] * B[(i*36)+34];
C[(i*36)+19] += values[103] * B[(i*36)+34];
C[(i*36)+34] += values[104] * B[(i*36)+34];
#endif

}

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 1890;
#endif

}

void dsparse_fP113DivM_m35_n9_k35_ldAna5_ldB36_ldC36_beta0_pfsigonly(const double* values, const double* B, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma nounroll_and_jam
for (unsigned int i = 0; i < 9; i++)
{
  #pragma simd
  for (unsigned int m = 0; m < 35; m++) {
    C[(i*36)+m] = 0.0;
  }
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b0 = _mm256_broadcast_sd(&B[(i*36)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b0 = _mm_loaddup_pd(&B[(i*36)+0]);
#endif
__m128d c0_0 = _mm_load_sd(&C[(i*36)+0]);
__m128d a0_0 = _mm_load_sd(&values[0]);
#if defined(__SSE3__) && defined(__AVX__)
c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, b0));
#endif
_mm_store_sd(&C[(i*36)+0], c0_0);
__m128d c0_1 = _mm_load_sd(&C[(i*36)+3]);
__m128d a0_1 = _mm_load_sd(&values[1]);
#if defined(__SSE3__) && defined(__AVX__)
c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, b0));
#endif
_mm_store_sd(&C[(i*36)+3], c0_1);
__m128d c0_2 = _mm_load_sd(&C[(i*36)+9]);
__m128d a0_2 = _mm_load_sd(&values[2]);
#if defined(__SSE3__) && defined(__AVX__)
c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, b0));
#endif
_mm_store_sd(&C[(i*36)+9], c0_2);
__m128d c0_3 = _mm_load_sd(&C[(i*36)+19]);
__m128d a0_3 = _mm_load_sd(&values[3]);
#if defined(__SSE3__) && defined(__AVX__)
c0_3 = _mm_add_sd(c0_3, _mm_mul_sd(a0_3, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_3 = _mm_add_sd(c0_3, _mm_mul_sd(a0_3, b0));
#endif
_mm_store_sd(&C[(i*36)+19], c0_3);
__m128d c0_4 = _mm_load_sd(&C[(i*36)+34]);
__m128d a0_4 = _mm_load_sd(&values[4]);
#if defined(__SSE3__) && defined(__AVX__)
c0_4 = _mm_add_sd(c0_4, _mm_mul_sd(a0_4, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_4 = _mm_add_sd(c0_4, _mm_mul_sd(a0_4, b0));
#endif
_mm_store_sd(&C[(i*36)+34], c0_4);
#else
C[(i*36)+0] += values[0] * B[(i*36)+0];
C[(i*36)+3] += values[1] * B[(i*36)+0];
C[(i*36)+9] += values[2] * B[(i*36)+0];
C[(i*36)+19] += values[3] * B[(i*36)+0];
C[(i*36)+34] += values[4] * B[(i*36)+0];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b1 = _mm256_broadcast_sd(&B[(i*36)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b1 = _mm_loaddup_pd(&B[(i*36)+1]);
#endif
__m128d c1_0 = _mm_load_sd(&C[(i*36)+1]);
__m128d a1_0 = _mm_load_sd(&values[5]);
#if defined(__SSE3__) && defined(__AVX__)
c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, b1));
#endif
_mm_store_sd(&C[(i*36)+1], c1_0);
__m128d c1_1 = _mm_load_sd(&C[(i*36)+7]);
__m128d a1_1 = _mm_load_sd(&values[6]);
#if defined(__SSE3__) && defined(__AVX__)
c1_1 = _mm_add_sd(c1_1, _mm_mul_sd(a1_1, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_1 = _mm_add_sd(c1_1, _mm_mul_sd(a1_1, b1));
#endif
_mm_store_sd(&C[(i*36)+7], c1_1);
__m128d c1_2 = _mm_load_sd(&C[(i*36)+17]);
__m128d a1_2 = _mm_load_sd(&values[7]);
#if defined(__SSE3__) && defined(__AVX__)
c1_2 = _mm_add_sd(c1_2, _mm_mul_sd(a1_2, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_2 = _mm_add_sd(c1_2, _mm_mul_sd(a1_2, b1));
#endif
_mm_store_sd(&C[(i*36)+17], c1_2);
__m128d c1_3 = _mm_load_sd(&C[(i*36)+32]);
__m128d a1_3 = _mm_load_sd(&values[8]);
#if defined(__SSE3__) && defined(__AVX__)
c1_3 = _mm_add_sd(c1_3, _mm_mul_sd(a1_3, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_3 = _mm_add_sd(c1_3, _mm_mul_sd(a1_3, b1));
#endif
_mm_store_sd(&C[(i*36)+32], c1_3);
#else
C[(i*36)+1] += values[5] * B[(i*36)+1];
C[(i*36)+7] += values[6] * B[(i*36)+1];
C[(i*36)+17] += values[7] * B[(i*36)+1];
C[(i*36)+32] += values[8] * B[(i*36)+1];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b2 = _mm256_broadcast_sd(&B[(i*36)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b2 = _mm_loaddup_pd(&B[(i*36)+2]);
#endif
__m128d c2_0 = _mm_load_sd(&C[(i*36)+2]);
__m128d a2_0 = _mm_load_sd(&values[9]);
#if defined(__SSE3__) && defined(__AVX__)
c2_0 = _mm_add_sd(c2_0, _mm_mul_sd(a2_0, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_0 = _mm_add_sd(c2_0, _mm_mul_sd(a2_0, b2));
#endif
_mm_store_sd(&C[(i*36)+2], c2_0);
__m128d c2_1 = _mm_load_sd(&C[(i*36)+8]);
__m128d a2_1 = _mm_load_sd(&values[10]);
#if defined(__SSE3__) && defined(__AVX__)
c2_1 = _mm_add_sd(c2_1, _mm_mul_sd(a2_1, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_1 = _mm_add_sd(c2_1, _mm_mul_sd(a2_1, b2));
#endif
_mm_store_sd(&C[(i*36)+8], c2_1);
__m128d c2_2 = _mm_load_sd(&C[(i*36)+18]);
__m128d a2_2 = _mm_load_sd(&values[11]);
#if defined(__SSE3__) && defined(__AVX__)
c2_2 = _mm_add_sd(c2_2, _mm_mul_sd(a2_2, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_2 = _mm_add_sd(c2_2, _mm_mul_sd(a2_2, b2));
#endif
_mm_store_sd(&C[(i*36)+18], c2_2);
__m128d c2_3 = _mm_load_sd(&C[(i*36)+33]);
__m128d a2_3 = _mm_load_sd(&values[12]);
#if defined(__SSE3__) && defined(__AVX__)
c2_3 = _mm_add_sd(c2_3, _mm_mul_sd(a2_3, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_3 = _mm_add_sd(c2_3, _mm_mul_sd(a2_3, b2));
#endif
_mm_store_sd(&C[(i*36)+33], c2_3);
#else
C[(i*36)+2] += values[9] * B[(i*36)+2];
C[(i*36)+8] += values[10] * B[(i*36)+2];
C[(i*36)+18] += values[11] * B[(i*36)+2];
C[(i*36)+33] += values[12] * B[(i*36)+2];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b3 = _mm256_broadcast_sd(&B[(i*36)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b3 = _mm_loaddup_pd(&B[(i*36)+3]);
#endif
__m128d c3_0 = _mm_load_sd(&C[(i*36)+0]);
__m128d a3_0 = _mm_load_sd(&values[13]);
#if defined(__SSE3__) && defined(__AVX__)
c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, b3));
#endif
_mm_store_sd(&C[(i*36)+0], c3_0);
__m128d c3_1 = _mm_load_sd(&C[(i*36)+3]);
__m128d a3_1 = _mm_load_sd(&values[14]);
#if defined(__SSE3__) && defined(__AVX__)
c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, b3));
#endif
_mm_store_sd(&C[(i*36)+3], c3_1);
__m128d c3_2 = _mm_load_sd(&C[(i*36)+9]);
__m128d a3_2 = _mm_load_sd(&values[15]);
#if defined(__SSE3__) && defined(__AVX__)
c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, b3));
#endif
_mm_store_sd(&C[(i*36)+9], c3_2);
__m128d c3_3 = _mm_load_sd(&C[(i*36)+19]);
__m128d a3_3 = _mm_load_sd(&values[16]);
#if defined(__SSE3__) && defined(__AVX__)
c3_3 = _mm_add_sd(c3_3, _mm_mul_sd(a3_3, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_3 = _mm_add_sd(c3_3, _mm_mul_sd(a3_3, b3));
#endif
_mm_store_sd(&C[(i*36)+19], c3_3);
__m128d c3_4 = _mm_load_sd(&C[(i*36)+34]);
__m128d a3_4 = _mm_load_sd(&values[17]);
#if defined(__SSE3__) && defined(__AVX__)
c3_4 = _mm_add_sd(c3_4, _mm_mul_sd(a3_4, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_4 = _mm_add_sd(c3_4, _mm_mul_sd(a3_4, b3));
#endif
_mm_store_sd(&C[(i*36)+34], c3_4);
#else
C[(i*36)+0] += values[13] * B[(i*36)+3];
C[(i*36)+3] += values[14] * B[(i*36)+3];
C[(i*36)+9] += values[15] * B[(i*36)+3];
C[(i*36)+19] += values[16] * B[(i*36)+3];
C[(i*36)+34] += values[17] * B[(i*36)+3];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b4 = _mm256_broadcast_sd(&B[(i*36)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b4 = _mm_loaddup_pd(&B[(i*36)+4]);
#endif
__m128d c4_0 = _mm_load_sd(&C[(i*36)+4]);
__m128d a4_0 = _mm_load_sd(&values[18]);
#if defined(__SSE3__) && defined(__AVX__)
c4_0 = _mm_add_sd(c4_0, _mm_mul_sd(a4_0, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_0 = _mm_add_sd(c4_0, _mm_mul_sd(a4_0, b4));
#endif
_mm_store_sd(&C[(i*36)+4], c4_0);
__m128d c4_1 = _mm_load_sd(&C[(i*36)+14]);
__m128d a4_1 = _mm_load_sd(&values[19]);
#if defined(__SSE3__) && defined(__AVX__)
c4_1 = _mm_add_sd(c4_1, _mm_mul_sd(a4_1, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_1 = _mm_add_sd(c4_1, _mm_mul_sd(a4_1, b4));
#endif
_mm_store_sd(&C[(i*36)+14], c4_1);
__m128d c4_2 = _mm_load_sd(&C[(i*36)+29]);
__m128d a4_2 = _mm_load_sd(&values[20]);
#if defined(__SSE3__) && defined(__AVX__)
c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, b4));
#endif
_mm_store_sd(&C[(i*36)+29], c4_2);
#else
C[(i*36)+4] += values[18] * B[(i*36)+4];
C[(i*36)+14] += values[19] * B[(i*36)+4];
C[(i*36)+29] += values[20] * B[(i*36)+4];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b5 = _mm256_broadcast_sd(&B[(i*36)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b5 = _mm_loaddup_pd(&B[(i*36)+5]);
#endif
__m128d c5_0 = _mm_load_sd(&C[(i*36)+5]);
__m128d a5_0 = _mm_load_sd(&values[21]);
#if defined(__SSE3__) && defined(__AVX__)
c5_0 = _mm_add_sd(c5_0, _mm_mul_sd(a5_0, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_0 = _mm_add_sd(c5_0, _mm_mul_sd(a5_0, b5));
#endif
_mm_store_sd(&C[(i*36)+5], c5_0);
__m128d c5_1 = _mm_load_sd(&C[(i*36)+15]);
__m128d a5_1 = _mm_load_sd(&values[22]);
#if defined(__SSE3__) && defined(__AVX__)
c5_1 = _mm_add_sd(c5_1, _mm_mul_sd(a5_1, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_1 = _mm_add_sd(c5_1, _mm_mul_sd(a5_1, b5));
#endif
_mm_store_sd(&C[(i*36)+15], c5_1);
__m128d c5_2 = _mm_load_sd(&C[(i*36)+30]);
__m128d a5_2 = _mm_load_sd(&values[23]);
#if defined(__SSE3__) && defined(__AVX__)
c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, b5));
#endif
_mm_store_sd(&C[(i*36)+30], c5_2);
#else
C[(i*36)+5] += values[21] * B[(i*36)+5];
C[(i*36)+15] += values[22] * B[(i*36)+5];
C[(i*36)+30] += values[23] * B[(i*36)+5];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b6 = _mm256_broadcast_sd(&B[(i*36)+6]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b6 = _mm_loaddup_pd(&B[(i*36)+6]);
#endif
__m128d c6_0 = _mm_load_sd(&C[(i*36)+6]);
__m128d a6_0 = _mm_load_sd(&values[24]);
#if defined(__SSE3__) && defined(__AVX__)
c6_0 = _mm_add_sd(c6_0, _mm_mul_sd(a6_0, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_0 = _mm_add_sd(c6_0, _mm_mul_sd(a6_0, b6));
#endif
_mm_store_sd(&C[(i*36)+6], c6_0);
__m128d c6_1 = _mm_load_sd(&C[(i*36)+16]);
__m128d a6_1 = _mm_load_sd(&values[25]);
#if defined(__SSE3__) && defined(__AVX__)
c6_1 = _mm_add_sd(c6_1, _mm_mul_sd(a6_1, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_1 = _mm_add_sd(c6_1, _mm_mul_sd(a6_1, b6));
#endif
_mm_store_sd(&C[(i*36)+16], c6_1);
__m128d c6_2 = _mm_load_sd(&C[(i*36)+31]);
__m128d a6_2 = _mm_load_sd(&values[26]);
#if defined(__SSE3__) && defined(__AVX__)
c6_2 = _mm_add_sd(c6_2, _mm_mul_sd(a6_2, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_2 = _mm_add_sd(c6_2, _mm_mul_sd(a6_2, b6));
#endif
_mm_store_sd(&C[(i*36)+31], c6_2);
#else
C[(i*36)+6] += values[24] * B[(i*36)+6];
C[(i*36)+16] += values[25] * B[(i*36)+6];
C[(i*36)+31] += values[26] * B[(i*36)+6];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b7 = _mm256_broadcast_sd(&B[(i*36)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b7 = _mm_loaddup_pd(&B[(i*36)+7]);
#endif
__m128d c7_0 = _mm_load_sd(&C[(i*36)+1]);
__m128d a7_0 = _mm_load_sd(&values[27]);
#if defined(__SSE3__) && defined(__AVX__)
c7_0 = _mm_add_sd(c7_0, _mm_mul_sd(a7_0, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_0 = _mm_add_sd(c7_0, _mm_mul_sd(a7_0, b7));
#endif
_mm_store_sd(&C[(i*36)+1], c7_0);
__m128d c7_1 = _mm_load_sd(&C[(i*36)+7]);
__m128d a7_1 = _mm_load_sd(&values[28]);
#if defined(__SSE3__) && defined(__AVX__)
c7_1 = _mm_add_sd(c7_1, _mm_mul_sd(a7_1, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_1 = _mm_add_sd(c7_1, _mm_mul_sd(a7_1, b7));
#endif
_mm_store_sd(&C[(i*36)+7], c7_1);
__m128d c7_2 = _mm_load_sd(&C[(i*36)+17]);
__m128d a7_2 = _mm_load_sd(&values[29]);
#if defined(__SSE3__) && defined(__AVX__)
c7_2 = _mm_add_sd(c7_2, _mm_mul_sd(a7_2, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_2 = _mm_add_sd(c7_2, _mm_mul_sd(a7_2, b7));
#endif
_mm_store_sd(&C[(i*36)+17], c7_2);
__m128d c7_3 = _mm_load_sd(&C[(i*36)+32]);
__m128d a7_3 = _mm_load_sd(&values[30]);
#if defined(__SSE3__) && defined(__AVX__)
c7_3 = _mm_add_sd(c7_3, _mm_mul_sd(a7_3, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_3 = _mm_add_sd(c7_3, _mm_mul_sd(a7_3, b7));
#endif
_mm_store_sd(&C[(i*36)+32], c7_3);
#else
C[(i*36)+1] += values[27] * B[(i*36)+7];
C[(i*36)+7] += values[28] * B[(i*36)+7];
C[(i*36)+17] += values[29] * B[(i*36)+7];
C[(i*36)+32] += values[30] * B[(i*36)+7];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b8 = _mm256_broadcast_sd(&B[(i*36)+8]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b8 = _mm_loaddup_pd(&B[(i*36)+8]);
#endif
__m128d c8_0 = _mm_load_sd(&C[(i*36)+2]);
__m128d a8_0 = _mm_load_sd(&values[31]);
#if defined(__SSE3__) && defined(__AVX__)
c8_0 = _mm_add_sd(c8_0, _mm_mul_sd(a8_0, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_0 = _mm_add_sd(c8_0, _mm_mul_sd(a8_0, b8));
#endif
_mm_store_sd(&C[(i*36)+2], c8_0);
__m128d c8_1 = _mm_load_sd(&C[(i*36)+8]);
__m128d a8_1 = _mm_load_sd(&values[32]);
#if defined(__SSE3__) && defined(__AVX__)
c8_1 = _mm_add_sd(c8_1, _mm_mul_sd(a8_1, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_1 = _mm_add_sd(c8_1, _mm_mul_sd(a8_1, b8));
#endif
_mm_store_sd(&C[(i*36)+8], c8_1);
__m128d c8_2 = _mm_load_sd(&C[(i*36)+18]);
__m128d a8_2 = _mm_load_sd(&values[33]);
#if defined(__SSE3__) && defined(__AVX__)
c8_2 = _mm_add_sd(c8_2, _mm_mul_sd(a8_2, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_2 = _mm_add_sd(c8_2, _mm_mul_sd(a8_2, b8));
#endif
_mm_store_sd(&C[(i*36)+18], c8_2);
__m128d c8_3 = _mm_load_sd(&C[(i*36)+33]);
__m128d a8_3 = _mm_load_sd(&values[34]);
#if defined(__SSE3__) && defined(__AVX__)
c8_3 = _mm_add_sd(c8_3, _mm_mul_sd(a8_3, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_3 = _mm_add_sd(c8_3, _mm_mul_sd(a8_3, b8));
#endif
_mm_store_sd(&C[(i*36)+33], c8_3);
#else
C[(i*36)+2] += values[31] * B[(i*36)+8];
C[(i*36)+8] += values[32] * B[(i*36)+8];
C[(i*36)+18] += values[33] * B[(i*36)+8];
C[(i*36)+33] += values[34] * B[(i*36)+8];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b9 = _mm256_broadcast_sd(&B[(i*36)+9]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b9 = _mm_loaddup_pd(&B[(i*36)+9]);
#endif
__m128d c9_0 = _mm_load_sd(&C[(i*36)+0]);
__m128d a9_0 = _mm_load_sd(&values[35]);
#if defined(__SSE3__) && defined(__AVX__)
c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, b9));
#endif
_mm_store_sd(&C[(i*36)+0], c9_0);
__m128d c9_1 = _mm_load_sd(&C[(i*36)+3]);
__m128d a9_1 = _mm_load_sd(&values[36]);
#if defined(__SSE3__) && defined(__AVX__)
c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, b9));
#endif
_mm_store_sd(&C[(i*36)+3], c9_1);
__m128d c9_2 = _mm_load_sd(&C[(i*36)+9]);
__m128d a9_2 = _mm_load_sd(&values[37]);
#if defined(__SSE3__) && defined(__AVX__)
c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, b9));
#endif
_mm_store_sd(&C[(i*36)+9], c9_2);
__m128d c9_3 = _mm_load_sd(&C[(i*36)+19]);
__m128d a9_3 = _mm_load_sd(&values[38]);
#if defined(__SSE3__) && defined(__AVX__)
c9_3 = _mm_add_sd(c9_3, _mm_mul_sd(a9_3, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_3 = _mm_add_sd(c9_3, _mm_mul_sd(a9_3, b9));
#endif
_mm_store_sd(&C[(i*36)+19], c9_3);
__m128d c9_4 = _mm_load_sd(&C[(i*36)+34]);
__m128d a9_4 = _mm_load_sd(&values[39]);
#if defined(__SSE3__) && defined(__AVX__)
c9_4 = _mm_add_sd(c9_4, _mm_mul_sd(a9_4, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_4 = _mm_add_sd(c9_4, _mm_mul_sd(a9_4, b9));
#endif
_mm_store_sd(&C[(i*36)+34], c9_4);
#else
C[(i*36)+0] += values[35] * B[(i*36)+9];
C[(i*36)+3] += values[36] * B[(i*36)+9];
C[(i*36)+9] += values[37] * B[(i*36)+9];
C[(i*36)+19] += values[38] * B[(i*36)+9];
C[(i*36)+34] += values[39] * B[(i*36)+9];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b10 = _mm256_broadcast_sd(&B[(i*36)+10]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b10 = _mm_loaddup_pd(&B[(i*36)+10]);
#endif
__m128d c10_0 = _mm_load_sd(&C[(i*36)+10]);
__m128d a10_0 = _mm_load_sd(&values[40]);
#if defined(__SSE3__) && defined(__AVX__)
c10_0 = _mm_add_sd(c10_0, _mm_mul_sd(a10_0, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c10_0 = _mm_add_sd(c10_0, _mm_mul_sd(a10_0, b10));
#endif
_mm_store_sd(&C[(i*36)+10], c10_0);
__m128d c10_1 = _mm_load_sd(&C[(i*36)+25]);
__m128d a10_1 = _mm_load_sd(&values[41]);
#if defined(__SSE3__) && defined(__AVX__)
c10_1 = _mm_add_sd(c10_1, _mm_mul_sd(a10_1, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c10_1 = _mm_add_sd(c10_1, _mm_mul_sd(a10_1, b10));
#endif
_mm_store_sd(&C[(i*36)+25], c10_1);
#else
C[(i*36)+10] += values[40] * B[(i*36)+10];
C[(i*36)+25] += values[41] * B[(i*36)+10];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b11 = _mm256_broadcast_sd(&B[(i*36)+11]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b11 = _mm_loaddup_pd(&B[(i*36)+11]);
#endif
__m128d c11_0 = _mm_load_sd(&C[(i*36)+11]);
__m128d a11_0 = _mm_load_sd(&values[42]);
#if defined(__SSE3__) && defined(__AVX__)
c11_0 = _mm_add_sd(c11_0, _mm_mul_sd(a11_0, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c11_0 = _mm_add_sd(c11_0, _mm_mul_sd(a11_0, b11));
#endif
_mm_store_sd(&C[(i*36)+11], c11_0);
__m128d c11_1 = _mm_load_sd(&C[(i*36)+26]);
__m128d a11_1 = _mm_load_sd(&values[43]);
#if defined(__SSE3__) && defined(__AVX__)
c11_1 = _mm_add_sd(c11_1, _mm_mul_sd(a11_1, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c11_1 = _mm_add_sd(c11_1, _mm_mul_sd(a11_1, b11));
#endif
_mm_store_sd(&C[(i*36)+26], c11_1);
#else
C[(i*36)+11] += values[42] * B[(i*36)+11];
C[(i*36)+26] += values[43] * B[(i*36)+11];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b12 = _mm256_broadcast_sd(&B[(i*36)+12]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b12 = _mm_loaddup_pd(&B[(i*36)+12]);
#endif
__m128d c12_0 = _mm_load_sd(&C[(i*36)+12]);
__m128d a12_0 = _mm_load_sd(&values[44]);
#if defined(__SSE3__) && defined(__AVX__)
c12_0 = _mm_add_sd(c12_0, _mm_mul_sd(a12_0, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c12_0 = _mm_add_sd(c12_0, _mm_mul_sd(a12_0, b12));
#endif
_mm_store_sd(&C[(i*36)+12], c12_0);
__m128d c12_1 = _mm_load_sd(&C[(i*36)+27]);
__m128d a12_1 = _mm_load_sd(&values[45]);
#if defined(__SSE3__) && defined(__AVX__)
c12_1 = _mm_add_sd(c12_1, _mm_mul_sd(a12_1, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c12_1 = _mm_add_sd(c12_1, _mm_mul_sd(a12_1, b12));
#endif
_mm_store_sd(&C[(i*36)+27], c12_1);
#else
C[(i*36)+12] += values[44] * B[(i*36)+12];
C[(i*36)+27] += values[45] * B[(i*36)+12];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b13 = _mm256_broadcast_sd(&B[(i*36)+13]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b13 = _mm_loaddup_pd(&B[(i*36)+13]);
#endif
__m128d c13_0 = _mm_load_sd(&C[(i*36)+13]);
__m128d a13_0 = _mm_load_sd(&values[46]);
#if defined(__SSE3__) && defined(__AVX__)
c13_0 = _mm_add_sd(c13_0, _mm_mul_sd(a13_0, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c13_0 = _mm_add_sd(c13_0, _mm_mul_sd(a13_0, b13));
#endif
_mm_store_sd(&C[(i*36)+13], c13_0);
__m128d c13_1 = _mm_load_sd(&C[(i*36)+28]);
__m128d a13_1 = _mm_load_sd(&values[47]);
#if defined(__SSE3__) && defined(__AVX__)
c13_1 = _mm_add_sd(c13_1, _mm_mul_sd(a13_1, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c13_1 = _mm_add_sd(c13_1, _mm_mul_sd(a13_1, b13));
#endif
_mm_store_sd(&C[(i*36)+28], c13_1);
#else
C[(i*36)+13] += values[46] * B[(i*36)+13];
C[(i*36)+28] += values[47] * B[(i*36)+13];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b14 = _mm256_broadcast_sd(&B[(i*36)+14]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b14 = _mm_loaddup_pd(&B[(i*36)+14]);
#endif
__m128d c14_0 = _mm_load_sd(&C[(i*36)+4]);
__m128d a14_0 = _mm_load_sd(&values[48]);
#if defined(__SSE3__) && defined(__AVX__)
c14_0 = _mm_add_sd(c14_0, _mm_mul_sd(a14_0, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_0 = _mm_add_sd(c14_0, _mm_mul_sd(a14_0, b14));
#endif
_mm_store_sd(&C[(i*36)+4], c14_0);
__m128d c14_1 = _mm_load_sd(&C[(i*36)+14]);
__m128d a14_1 = _mm_load_sd(&values[49]);
#if defined(__SSE3__) && defined(__AVX__)
c14_1 = _mm_add_sd(c14_1, _mm_mul_sd(a14_1, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_1 = _mm_add_sd(c14_1, _mm_mul_sd(a14_1, b14));
#endif
_mm_store_sd(&C[(i*36)+14], c14_1);
__m128d c14_2 = _mm_load_sd(&C[(i*36)+29]);
__m128d a14_2 = _mm_load_sd(&values[50]);
#if defined(__SSE3__) && defined(__AVX__)
c14_2 = _mm_add_sd(c14_2, _mm_mul_sd(a14_2, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_2 = _mm_add_sd(c14_2, _mm_mul_sd(a14_2, b14));
#endif
_mm_store_sd(&C[(i*36)+29], c14_2);
#else
C[(i*36)+4] += values[48] * B[(i*36)+14];
C[(i*36)+14] += values[49] * B[(i*36)+14];
C[(i*36)+29] += values[50] * B[(i*36)+14];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b15 = _mm256_broadcast_sd(&B[(i*36)+15]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b15 = _mm_loaddup_pd(&B[(i*36)+15]);
#endif
__m128d c15_0 = _mm_load_sd(&C[(i*36)+5]);
__m128d a15_0 = _mm_load_sd(&values[51]);
#if defined(__SSE3__) && defined(__AVX__)
c15_0 = _mm_add_sd(c15_0, _mm_mul_sd(a15_0, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_0 = _mm_add_sd(c15_0, _mm_mul_sd(a15_0, b15));
#endif
_mm_store_sd(&C[(i*36)+5], c15_0);
__m128d c15_1 = _mm_load_sd(&C[(i*36)+15]);
__m128d a15_1 = _mm_load_sd(&values[52]);
#if defined(__SSE3__) && defined(__AVX__)
c15_1 = _mm_add_sd(c15_1, _mm_mul_sd(a15_1, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_1 = _mm_add_sd(c15_1, _mm_mul_sd(a15_1, b15));
#endif
_mm_store_sd(&C[(i*36)+15], c15_1);
__m128d c15_2 = _mm_load_sd(&C[(i*36)+30]);
__m128d a15_2 = _mm_load_sd(&values[53]);
#if defined(__SSE3__) && defined(__AVX__)
c15_2 = _mm_add_sd(c15_2, _mm_mul_sd(a15_2, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_2 = _mm_add_sd(c15_2, _mm_mul_sd(a15_2, b15));
#endif
_mm_store_sd(&C[(i*36)+30], c15_2);
#else
C[(i*36)+5] += values[51] * B[(i*36)+15];
C[(i*36)+15] += values[52] * B[(i*36)+15];
C[(i*36)+30] += values[53] * B[(i*36)+15];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b16 = _mm256_broadcast_sd(&B[(i*36)+16]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b16 = _mm_loaddup_pd(&B[(i*36)+16]);
#endif
__m128d c16_0 = _mm_load_sd(&C[(i*36)+6]);
__m128d a16_0 = _mm_load_sd(&values[54]);
#if defined(__SSE3__) && defined(__AVX__)
c16_0 = _mm_add_sd(c16_0, _mm_mul_sd(a16_0, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_0 = _mm_add_sd(c16_0, _mm_mul_sd(a16_0, b16));
#endif
_mm_store_sd(&C[(i*36)+6], c16_0);
__m128d c16_1 = _mm_load_sd(&C[(i*36)+16]);
__m128d a16_1 = _mm_load_sd(&values[55]);
#if defined(__SSE3__) && defined(__AVX__)
c16_1 = _mm_add_sd(c16_1, _mm_mul_sd(a16_1, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_1 = _mm_add_sd(c16_1, _mm_mul_sd(a16_1, b16));
#endif
_mm_store_sd(&C[(i*36)+16], c16_1);
__m128d c16_2 = _mm_load_sd(&C[(i*36)+31]);
__m128d a16_2 = _mm_load_sd(&values[56]);
#if defined(__SSE3__) && defined(__AVX__)
c16_2 = _mm_add_sd(c16_2, _mm_mul_sd(a16_2, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_2 = _mm_add_sd(c16_2, _mm_mul_sd(a16_2, b16));
#endif
_mm_store_sd(&C[(i*36)+31], c16_2);
#else
C[(i*36)+6] += values[54] * B[(i*36)+16];
C[(i*36)+16] += values[55] * B[(i*36)+16];
C[(i*36)+31] += values[56] * B[(i*36)+16];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b17 = _mm256_broadcast_sd(&B[(i*36)+17]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b17 = _mm_loaddup_pd(&B[(i*36)+17]);
#endif
__m128d c17_0 = _mm_load_sd(&C[(i*36)+1]);
__m128d a17_0 = _mm_load_sd(&values[57]);
#if defined(__SSE3__) && defined(__AVX__)
c17_0 = _mm_add_sd(c17_0, _mm_mul_sd(a17_0, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_0 = _mm_add_sd(c17_0, _mm_mul_sd(a17_0, b17));
#endif
_mm_store_sd(&C[(i*36)+1], c17_0);
__m128d c17_1 = _mm_load_sd(&C[(i*36)+7]);
__m128d a17_1 = _mm_load_sd(&values[58]);
#if defined(__SSE3__) && defined(__AVX__)
c17_1 = _mm_add_sd(c17_1, _mm_mul_sd(a17_1, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_1 = _mm_add_sd(c17_1, _mm_mul_sd(a17_1, b17));
#endif
_mm_store_sd(&C[(i*36)+7], c17_1);
__m128d c17_2 = _mm_load_sd(&C[(i*36)+17]);
__m128d a17_2 = _mm_load_sd(&values[59]);
#if defined(__SSE3__) && defined(__AVX__)
c17_2 = _mm_add_sd(c17_2, _mm_mul_sd(a17_2, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_2 = _mm_add_sd(c17_2, _mm_mul_sd(a17_2, b17));
#endif
_mm_store_sd(&C[(i*36)+17], c17_2);
__m128d c17_3 = _mm_load_sd(&C[(i*36)+32]);
__m128d a17_3 = _mm_load_sd(&values[60]);
#if defined(__SSE3__) && defined(__AVX__)
c17_3 = _mm_add_sd(c17_3, _mm_mul_sd(a17_3, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_3 = _mm_add_sd(c17_3, _mm_mul_sd(a17_3, b17));
#endif
_mm_store_sd(&C[(i*36)+32], c17_3);
#else
C[(i*36)+1] += values[57] * B[(i*36)+17];
C[(i*36)+7] += values[58] * B[(i*36)+17];
C[(i*36)+17] += values[59] * B[(i*36)+17];
C[(i*36)+32] += values[60] * B[(i*36)+17];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b18 = _mm256_broadcast_sd(&B[(i*36)+18]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b18 = _mm_loaddup_pd(&B[(i*36)+18]);
#endif
__m128d c18_0 = _mm_load_sd(&C[(i*36)+2]);
__m128d a18_0 = _mm_load_sd(&values[61]);
#if defined(__SSE3__) && defined(__AVX__)
c18_0 = _mm_add_sd(c18_0, _mm_mul_sd(a18_0, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_0 = _mm_add_sd(c18_0, _mm_mul_sd(a18_0, b18));
#endif
_mm_store_sd(&C[(i*36)+2], c18_0);
__m128d c18_1 = _mm_load_sd(&C[(i*36)+8]);
__m128d a18_1 = _mm_load_sd(&values[62]);
#if defined(__SSE3__) && defined(__AVX__)
c18_1 = _mm_add_sd(c18_1, _mm_mul_sd(a18_1, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_1 = _mm_add_sd(c18_1, _mm_mul_sd(a18_1, b18));
#endif
_mm_store_sd(&C[(i*36)+8], c18_1);
__m128d c18_2 = _mm_load_sd(&C[(i*36)+18]);
__m128d a18_2 = _mm_load_sd(&values[63]);
#if defined(__SSE3__) && defined(__AVX__)
c18_2 = _mm_add_sd(c18_2, _mm_mul_sd(a18_2, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_2 = _mm_add_sd(c18_2, _mm_mul_sd(a18_2, b18));
#endif
_mm_store_sd(&C[(i*36)+18], c18_2);
__m128d c18_3 = _mm_load_sd(&C[(i*36)+33]);
__m128d a18_3 = _mm_load_sd(&values[64]);
#if defined(__SSE3__) && defined(__AVX__)
c18_3 = _mm_add_sd(c18_3, _mm_mul_sd(a18_3, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_3 = _mm_add_sd(c18_3, _mm_mul_sd(a18_3, b18));
#endif
_mm_store_sd(&C[(i*36)+33], c18_3);
#else
C[(i*36)+2] += values[61] * B[(i*36)+18];
C[(i*36)+8] += values[62] * B[(i*36)+18];
C[(i*36)+18] += values[63] * B[(i*36)+18];
C[(i*36)+33] += values[64] * B[(i*36)+18];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b19 = _mm256_broadcast_sd(&B[(i*36)+19]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b19 = _mm_loaddup_pd(&B[(i*36)+19]);
#endif
__m128d c19_0 = _mm_load_sd(&C[(i*36)+0]);
__m128d a19_0 = _mm_load_sd(&values[65]);
#if defined(__SSE3__) && defined(__AVX__)
c19_0 = _mm_add_sd(c19_0, _mm_mul_sd(a19_0, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_0 = _mm_add_sd(c19_0, _mm_mul_sd(a19_0, b19));
#endif
_mm_store_sd(&C[(i*36)+0], c19_0);
__m128d c19_1 = _mm_load_sd(&C[(i*36)+3]);
__m128d a19_1 = _mm_load_sd(&values[66]);
#if defined(__SSE3__) && defined(__AVX__)
c19_1 = _mm_add_sd(c19_1, _mm_mul_sd(a19_1, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_1 = _mm_add_sd(c19_1, _mm_mul_sd(a19_1, b19));
#endif
_mm_store_sd(&C[(i*36)+3], c19_1);
__m128d c19_2 = _mm_load_sd(&C[(i*36)+9]);
__m128d a19_2 = _mm_load_sd(&values[67]);
#if defined(__SSE3__) && defined(__AVX__)
c19_2 = _mm_add_sd(c19_2, _mm_mul_sd(a19_2, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_2 = _mm_add_sd(c19_2, _mm_mul_sd(a19_2, b19));
#endif
_mm_store_sd(&C[(i*36)+9], c19_2);
__m128d c19_3 = _mm_load_sd(&C[(i*36)+19]);
__m128d a19_3 = _mm_load_sd(&values[68]);
#if defined(__SSE3__) && defined(__AVX__)
c19_3 = _mm_add_sd(c19_3, _mm_mul_sd(a19_3, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_3 = _mm_add_sd(c19_3, _mm_mul_sd(a19_3, b19));
#endif
_mm_store_sd(&C[(i*36)+19], c19_3);
__m128d c19_4 = _mm_load_sd(&C[(i*36)+34]);
__m128d a19_4 = _mm_load_sd(&values[69]);
#if defined(__SSE3__) && defined(__AVX__)
c19_4 = _mm_add_sd(c19_4, _mm_mul_sd(a19_4, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_4 = _mm_add_sd(c19_4, _mm_mul_sd(a19_4, b19));
#endif
_mm_store_sd(&C[(i*36)+34], c19_4);
#else
C[(i*36)+0] += values[65] * B[(i*36)+19];
C[(i*36)+3] += values[66] * B[(i*36)+19];
C[(i*36)+9] += values[67] * B[(i*36)+19];
C[(i*36)+19] += values[68] * B[(i*36)+19];
C[(i*36)+34] += values[69] * B[(i*36)+19];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b20 = _mm256_broadcast_sd(&B[(i*36)+20]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b20 = _mm_loaddup_pd(&B[(i*36)+20]);
#endif
__m128d c20_0 = _mm_load_sd(&C[(i*36)+20]);
__m128d a20_0 = _mm_load_sd(&values[70]);
#if defined(__SSE3__) && defined(__AVX__)
c20_0 = _mm_add_sd(c20_0, _mm_mul_sd(a20_0, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c20_0 = _mm_add_sd(c20_0, _mm_mul_sd(a20_0, b20));
#endif
_mm_store_sd(&C[(i*36)+20], c20_0);
#else
C[(i*36)+20] += values[70] * B[(i*36)+20];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b21 = _mm256_broadcast_sd(&B[(i*36)+21]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b21 = _mm_loaddup_pd(&B[(i*36)+21]);
#endif
__m128d c21_0 = _mm_load_sd(&C[(i*36)+21]);
__m128d a21_0 = _mm_load_sd(&values[71]);
#if defined(__SSE3__) && defined(__AVX__)
c21_0 = _mm_add_sd(c21_0, _mm_mul_sd(a21_0, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c21_0 = _mm_add_sd(c21_0, _mm_mul_sd(a21_0, b21));
#endif
_mm_store_sd(&C[(i*36)+21], c21_0);
#else
C[(i*36)+21] += values[71] * B[(i*36)+21];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b22 = _mm256_broadcast_sd(&B[(i*36)+22]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b22 = _mm_loaddup_pd(&B[(i*36)+22]);
#endif
__m128d c22_0 = _mm_load_sd(&C[(i*36)+22]);
__m128d a22_0 = _mm_load_sd(&values[72]);
#if defined(__SSE3__) && defined(__AVX__)
c22_0 = _mm_add_sd(c22_0, _mm_mul_sd(a22_0, _mm256_castpd256_pd128(b22)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c22_0 = _mm_add_sd(c22_0, _mm_mul_sd(a22_0, b22));
#endif
_mm_store_sd(&C[(i*36)+22], c22_0);
#else
C[(i*36)+22] += values[72] * B[(i*36)+22];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b23 = _mm256_broadcast_sd(&B[(i*36)+23]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b23 = _mm_loaddup_pd(&B[(i*36)+23]);
#endif
__m128d c23_0 = _mm_load_sd(&C[(i*36)+23]);
__m128d a23_0 = _mm_load_sd(&values[73]);
#if defined(__SSE3__) && defined(__AVX__)
c23_0 = _mm_add_sd(c23_0, _mm_mul_sd(a23_0, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c23_0 = _mm_add_sd(c23_0, _mm_mul_sd(a23_0, b23));
#endif
_mm_store_sd(&C[(i*36)+23], c23_0);
#else
C[(i*36)+23] += values[73] * B[(i*36)+23];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b24 = _mm256_broadcast_sd(&B[(i*36)+24]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b24 = _mm_loaddup_pd(&B[(i*36)+24]);
#endif
__m128d c24_0 = _mm_load_sd(&C[(i*36)+24]);
__m128d a24_0 = _mm_load_sd(&values[74]);
#if defined(__SSE3__) && defined(__AVX__)
c24_0 = _mm_add_sd(c24_0, _mm_mul_sd(a24_0, _mm256_castpd256_pd128(b24)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c24_0 = _mm_add_sd(c24_0, _mm_mul_sd(a24_0, b24));
#endif
_mm_store_sd(&C[(i*36)+24], c24_0);
#else
C[(i*36)+24] += values[74] * B[(i*36)+24];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b25 = _mm256_broadcast_sd(&B[(i*36)+25]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b25 = _mm_loaddup_pd(&B[(i*36)+25]);
#endif
__m128d c25_0 = _mm_load_sd(&C[(i*36)+10]);
__m128d a25_0 = _mm_load_sd(&values[75]);
#if defined(__SSE3__) && defined(__AVX__)
c25_0 = _mm_add_sd(c25_0, _mm_mul_sd(a25_0, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c25_0 = _mm_add_sd(c25_0, _mm_mul_sd(a25_0, b25));
#endif
_mm_store_sd(&C[(i*36)+10], c25_0);
__m128d c25_1 = _mm_load_sd(&C[(i*36)+25]);
__m128d a25_1 = _mm_load_sd(&values[76]);
#if defined(__SSE3__) && defined(__AVX__)
c25_1 = _mm_add_sd(c25_1, _mm_mul_sd(a25_1, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c25_1 = _mm_add_sd(c25_1, _mm_mul_sd(a25_1, b25));
#endif
_mm_store_sd(&C[(i*36)+25], c25_1);
#else
C[(i*36)+10] += values[75] * B[(i*36)+25];
C[(i*36)+25] += values[76] * B[(i*36)+25];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b26 = _mm256_broadcast_sd(&B[(i*36)+26]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b26 = _mm_loaddup_pd(&B[(i*36)+26]);
#endif
__m128d c26_0 = _mm_load_sd(&C[(i*36)+11]);
__m128d a26_0 = _mm_load_sd(&values[77]);
#if defined(__SSE3__) && defined(__AVX__)
c26_0 = _mm_add_sd(c26_0, _mm_mul_sd(a26_0, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c26_0 = _mm_add_sd(c26_0, _mm_mul_sd(a26_0, b26));
#endif
_mm_store_sd(&C[(i*36)+11], c26_0);
__m128d c26_1 = _mm_load_sd(&C[(i*36)+26]);
__m128d a26_1 = _mm_load_sd(&values[78]);
#if defined(__SSE3__) && defined(__AVX__)
c26_1 = _mm_add_sd(c26_1, _mm_mul_sd(a26_1, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c26_1 = _mm_add_sd(c26_1, _mm_mul_sd(a26_1, b26));
#endif
_mm_store_sd(&C[(i*36)+26], c26_1);
#else
C[(i*36)+11] += values[77] * B[(i*36)+26];
C[(i*36)+26] += values[78] * B[(i*36)+26];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b27 = _mm256_broadcast_sd(&B[(i*36)+27]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b27 = _mm_loaddup_pd(&B[(i*36)+27]);
#endif
__m128d c27_0 = _mm_load_sd(&C[(i*36)+12]);
__m128d a27_0 = _mm_load_sd(&values[79]);
#if defined(__SSE3__) && defined(__AVX__)
c27_0 = _mm_add_sd(c27_0, _mm_mul_sd(a27_0, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c27_0 = _mm_add_sd(c27_0, _mm_mul_sd(a27_0, b27));
#endif
_mm_store_sd(&C[(i*36)+12], c27_0);
__m128d c27_1 = _mm_load_sd(&C[(i*36)+27]);
__m128d a27_1 = _mm_load_sd(&values[80]);
#if defined(__SSE3__) && defined(__AVX__)
c27_1 = _mm_add_sd(c27_1, _mm_mul_sd(a27_1, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c27_1 = _mm_add_sd(c27_1, _mm_mul_sd(a27_1, b27));
#endif
_mm_store_sd(&C[(i*36)+27], c27_1);
#else
C[(i*36)+12] += values[79] * B[(i*36)+27];
C[(i*36)+27] += values[80] * B[(i*36)+27];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b28 = _mm256_broadcast_sd(&B[(i*36)+28]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b28 = _mm_loaddup_pd(&B[(i*36)+28]);
#endif
__m128d c28_0 = _mm_load_sd(&C[(i*36)+13]);
__m128d a28_0 = _mm_load_sd(&values[81]);
#if defined(__SSE3__) && defined(__AVX__)
c28_0 = _mm_add_sd(c28_0, _mm_mul_sd(a28_0, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c28_0 = _mm_add_sd(c28_0, _mm_mul_sd(a28_0, b28));
#endif
_mm_store_sd(&C[(i*36)+13], c28_0);
__m128d c28_1 = _mm_load_sd(&C[(i*36)+28]);
__m128d a28_1 = _mm_load_sd(&values[82]);
#if defined(__SSE3__) && defined(__AVX__)
c28_1 = _mm_add_sd(c28_1, _mm_mul_sd(a28_1, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c28_1 = _mm_add_sd(c28_1, _mm_mul_sd(a28_1, b28));
#endif
_mm_store_sd(&C[(i*36)+28], c28_1);
#else
C[(i*36)+13] += values[81] * B[(i*36)+28];
C[(i*36)+28] += values[82] * B[(i*36)+28];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b29 = _mm256_broadcast_sd(&B[(i*36)+29]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b29 = _mm_loaddup_pd(&B[(i*36)+29]);
#endif
__m128d c29_0 = _mm_load_sd(&C[(i*36)+4]);
__m128d a29_0 = _mm_load_sd(&values[83]);
#if defined(__SSE3__) && defined(__AVX__)
c29_0 = _mm_add_sd(c29_0, _mm_mul_sd(a29_0, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_0 = _mm_add_sd(c29_0, _mm_mul_sd(a29_0, b29));
#endif
_mm_store_sd(&C[(i*36)+4], c29_0);
__m128d c29_1 = _mm_load_sd(&C[(i*36)+14]);
__m128d a29_1 = _mm_load_sd(&values[84]);
#if defined(__SSE3__) && defined(__AVX__)
c29_1 = _mm_add_sd(c29_1, _mm_mul_sd(a29_1, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_1 = _mm_add_sd(c29_1, _mm_mul_sd(a29_1, b29));
#endif
_mm_store_sd(&C[(i*36)+14], c29_1);
__m128d c29_2 = _mm_load_sd(&C[(i*36)+29]);
__m128d a29_2 = _mm_load_sd(&values[85]);
#if defined(__SSE3__) && defined(__AVX__)
c29_2 = _mm_add_sd(c29_2, _mm_mul_sd(a29_2, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_2 = _mm_add_sd(c29_2, _mm_mul_sd(a29_2, b29));
#endif
_mm_store_sd(&C[(i*36)+29], c29_2);
#else
C[(i*36)+4] += values[83] * B[(i*36)+29];
C[(i*36)+14] += values[84] * B[(i*36)+29];
C[(i*36)+29] += values[85] * B[(i*36)+29];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b30 = _mm256_broadcast_sd(&B[(i*36)+30]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b30 = _mm_loaddup_pd(&B[(i*36)+30]);
#endif
__m128d c30_0 = _mm_load_sd(&C[(i*36)+5]);
__m128d a30_0 = _mm_load_sd(&values[86]);
#if defined(__SSE3__) && defined(__AVX__)
c30_0 = _mm_add_sd(c30_0, _mm_mul_sd(a30_0, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_0 = _mm_add_sd(c30_0, _mm_mul_sd(a30_0, b30));
#endif
_mm_store_sd(&C[(i*36)+5], c30_0);
__m128d c30_1 = _mm_load_sd(&C[(i*36)+15]);
__m128d a30_1 = _mm_load_sd(&values[87]);
#if defined(__SSE3__) && defined(__AVX__)
c30_1 = _mm_add_sd(c30_1, _mm_mul_sd(a30_1, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_1 = _mm_add_sd(c30_1, _mm_mul_sd(a30_1, b30));
#endif
_mm_store_sd(&C[(i*36)+15], c30_1);
__m128d c30_2 = _mm_load_sd(&C[(i*36)+30]);
__m128d a30_2 = _mm_load_sd(&values[88]);
#if defined(__SSE3__) && defined(__AVX__)
c30_2 = _mm_add_sd(c30_2, _mm_mul_sd(a30_2, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_2 = _mm_add_sd(c30_2, _mm_mul_sd(a30_2, b30));
#endif
_mm_store_sd(&C[(i*36)+30], c30_2);
#else
C[(i*36)+5] += values[86] * B[(i*36)+30];
C[(i*36)+15] += values[87] * B[(i*36)+30];
C[(i*36)+30] += values[88] * B[(i*36)+30];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b31 = _mm256_broadcast_sd(&B[(i*36)+31]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b31 = _mm_loaddup_pd(&B[(i*36)+31]);
#endif
__m128d c31_0 = _mm_load_sd(&C[(i*36)+6]);
__m128d a31_0 = _mm_load_sd(&values[89]);
#if defined(__SSE3__) && defined(__AVX__)
c31_0 = _mm_add_sd(c31_0, _mm_mul_sd(a31_0, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_0 = _mm_add_sd(c31_0, _mm_mul_sd(a31_0, b31));
#endif
_mm_store_sd(&C[(i*36)+6], c31_0);
__m128d c31_1 = _mm_load_sd(&C[(i*36)+16]);
__m128d a31_1 = _mm_load_sd(&values[90]);
#if defined(__SSE3__) && defined(__AVX__)
c31_1 = _mm_add_sd(c31_1, _mm_mul_sd(a31_1, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_1 = _mm_add_sd(c31_1, _mm_mul_sd(a31_1, b31));
#endif
_mm_store_sd(&C[(i*36)+16], c31_1);
__m128d c31_2 = _mm_load_sd(&C[(i*36)+31]);
__m128d a31_2 = _mm_load_sd(&values[91]);
#if defined(__SSE3__) && defined(__AVX__)
c31_2 = _mm_add_sd(c31_2, _mm_mul_sd(a31_2, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_2 = _mm_add_sd(c31_2, _mm_mul_sd(a31_2, b31));
#endif
_mm_store_sd(&C[(i*36)+31], c31_2);
#else
C[(i*36)+6] += values[89] * B[(i*36)+31];
C[(i*36)+16] += values[90] * B[(i*36)+31];
C[(i*36)+31] += values[91] * B[(i*36)+31];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b32 = _mm256_broadcast_sd(&B[(i*36)+32]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b32 = _mm_loaddup_pd(&B[(i*36)+32]);
#endif
__m128d c32_0 = _mm_load_sd(&C[(i*36)+1]);
__m128d a32_0 = _mm_load_sd(&values[92]);
#if defined(__SSE3__) && defined(__AVX__)
c32_0 = _mm_add_sd(c32_0, _mm_mul_sd(a32_0, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_0 = _mm_add_sd(c32_0, _mm_mul_sd(a32_0, b32));
#endif
_mm_store_sd(&C[(i*36)+1], c32_0);
__m128d c32_1 = _mm_load_sd(&C[(i*36)+7]);
__m128d a32_1 = _mm_load_sd(&values[93]);
#if defined(__SSE3__) && defined(__AVX__)
c32_1 = _mm_add_sd(c32_1, _mm_mul_sd(a32_1, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_1 = _mm_add_sd(c32_1, _mm_mul_sd(a32_1, b32));
#endif
_mm_store_sd(&C[(i*36)+7], c32_1);
__m128d c32_2 = _mm_load_sd(&C[(i*36)+17]);
__m128d a32_2 = _mm_load_sd(&values[94]);
#if defined(__SSE3__) && defined(__AVX__)
c32_2 = _mm_add_sd(c32_2, _mm_mul_sd(a32_2, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_2 = _mm_add_sd(c32_2, _mm_mul_sd(a32_2, b32));
#endif
_mm_store_sd(&C[(i*36)+17], c32_2);
__m128d c32_3 = _mm_load_sd(&C[(i*36)+32]);
__m128d a32_3 = _mm_load_sd(&values[95]);
#if defined(__SSE3__) && defined(__AVX__)
c32_3 = _mm_add_sd(c32_3, _mm_mul_sd(a32_3, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_3 = _mm_add_sd(c32_3, _mm_mul_sd(a32_3, b32));
#endif
_mm_store_sd(&C[(i*36)+32], c32_3);
#else
C[(i*36)+1] += values[92] * B[(i*36)+32];
C[(i*36)+7] += values[93] * B[(i*36)+32];
C[(i*36)+17] += values[94] * B[(i*36)+32];
C[(i*36)+32] += values[95] * B[(i*36)+32];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b33 = _mm256_broadcast_sd(&B[(i*36)+33]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b33 = _mm_loaddup_pd(&B[(i*36)+33]);
#endif
__m128d c33_0 = _mm_load_sd(&C[(i*36)+2]);
__m128d a33_0 = _mm_load_sd(&values[96]);
#if defined(__SSE3__) && defined(__AVX__)
c33_0 = _mm_add_sd(c33_0, _mm_mul_sd(a33_0, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_0 = _mm_add_sd(c33_0, _mm_mul_sd(a33_0, b33));
#endif
_mm_store_sd(&C[(i*36)+2], c33_0);
__m128d c33_1 = _mm_load_sd(&C[(i*36)+8]);
__m128d a33_1 = _mm_load_sd(&values[97]);
#if defined(__SSE3__) && defined(__AVX__)
c33_1 = _mm_add_sd(c33_1, _mm_mul_sd(a33_1, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_1 = _mm_add_sd(c33_1, _mm_mul_sd(a33_1, b33));
#endif
_mm_store_sd(&C[(i*36)+8], c33_1);
__m128d c33_2 = _mm_load_sd(&C[(i*36)+18]);
__m128d a33_2 = _mm_load_sd(&values[98]);
#if defined(__SSE3__) && defined(__AVX__)
c33_2 = _mm_add_sd(c33_2, _mm_mul_sd(a33_2, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_2 = _mm_add_sd(c33_2, _mm_mul_sd(a33_2, b33));
#endif
_mm_store_sd(&C[(i*36)+18], c33_2);
__m128d c33_3 = _mm_load_sd(&C[(i*36)+33]);
__m128d a33_3 = _mm_load_sd(&values[99]);
#if defined(__SSE3__) && defined(__AVX__)
c33_3 = _mm_add_sd(c33_3, _mm_mul_sd(a33_3, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_3 = _mm_add_sd(c33_3, _mm_mul_sd(a33_3, b33));
#endif
_mm_store_sd(&C[(i*36)+33], c33_3);
#else
C[(i*36)+2] += values[96] * B[(i*36)+33];
C[(i*36)+8] += values[97] * B[(i*36)+33];
C[(i*36)+18] += values[98] * B[(i*36)+33];
C[(i*36)+33] += values[99] * B[(i*36)+33];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b34 = _mm256_broadcast_sd(&B[(i*36)+34]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b34 = _mm_loaddup_pd(&B[(i*36)+34]);
#endif
__m128d c34_0 = _mm_load_sd(&C[(i*36)+0]);
__m128d a34_0 = _mm_load_sd(&values[100]);
#if defined(__SSE3__) && defined(__AVX__)
c34_0 = _mm_add_sd(c34_0, _mm_mul_sd(a34_0, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_0 = _mm_add_sd(c34_0, _mm_mul_sd(a34_0, b34));
#endif
_mm_store_sd(&C[(i*36)+0], c34_0);
__m128d c34_1 = _mm_load_sd(&C[(i*36)+3]);
__m128d a34_1 = _mm_load_sd(&values[101]);
#if defined(__SSE3__) && defined(__AVX__)
c34_1 = _mm_add_sd(c34_1, _mm_mul_sd(a34_1, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_1 = _mm_add_sd(c34_1, _mm_mul_sd(a34_1, b34));
#endif
_mm_store_sd(&C[(i*36)+3], c34_1);
__m128d c34_2 = _mm_load_sd(&C[(i*36)+9]);
__m128d a34_2 = _mm_load_sd(&values[102]);
#if defined(__SSE3__) && defined(__AVX__)
c34_2 = _mm_add_sd(c34_2, _mm_mul_sd(a34_2, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_2 = _mm_add_sd(c34_2, _mm_mul_sd(a34_2, b34));
#endif
_mm_store_sd(&C[(i*36)+9], c34_2);
__m128d c34_3 = _mm_load_sd(&C[(i*36)+19]);
__m128d a34_3 = _mm_load_sd(&values[103]);
#if defined(__SSE3__) && defined(__AVX__)
c34_3 = _mm_add_sd(c34_3, _mm_mul_sd(a34_3, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_3 = _mm_add_sd(c34_3, _mm_mul_sd(a34_3, b34));
#endif
_mm_store_sd(&C[(i*36)+19], c34_3);
__m128d c34_4 = _mm_load_sd(&C[(i*36)+34]);
__m128d a34_4 = _mm_load_sd(&values[104]);
#if defined(__SSE3__) && defined(__AVX__)
c34_4 = _mm_add_sd(c34_4, _mm_mul_sd(a34_4, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_4 = _mm_add_sd(c34_4, _mm_mul_sd(a34_4, b34));
#endif
_mm_store_sd(&C[(i*36)+34], c34_4);
#else
C[(i*36)+0] += values[100] * B[(i*36)+34];
C[(i*36)+3] += values[101] * B[(i*36)+34];
C[(i*36)+9] += values[102] * B[(i*36)+34];
C[(i*36)+19] += values[103] * B[(i*36)+34];
C[(i*36)+34] += values[104] * B[(i*36)+34];
#endif

}

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 1890;
#endif

}

void dsparse_starMatrix_m56_n9_k9_ldA56_ldBna6_ldC56_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
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

void dsparse_fM1DivM_m56_n9_k56_ldAna6_ldB56_ldC56_beta0_pfsigonly(const double* values, const double* B, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
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

void dsparse_starMatrix_m84_n9_k9_ldA84_ldBna7_ldC84_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
#pragma vector aligned
for (unsigned int i = 0; i < 84; i += 1)
{
C[(i)+(0)] += A[(i)+(504)] * values[0];
C[(i)+(0)] += A[(i)+(588)] * values[1];
C[(i)+(0)] += A[(i)+(672)] * values[2];
C[(i)+(84)] += A[(i)+(504)] * values[3];
C[(i)+(84)] += A[(i)+(588)] * values[4];
C[(i)+(84)] += A[(i)+(672)] * values[5];
C[(i)+(168)] += A[(i)+(504)] * values[6];
C[(i)+(168)] += A[(i)+(588)] * values[7];
C[(i)+(168)] += A[(i)+(672)] * values[8];
C[(i)+(252)] += A[(i)+(504)] * values[9];
C[(i)+(252)] += A[(i)+(588)] * values[10];
C[(i)+(336)] += A[(i)+(588)] * values[11];
C[(i)+(336)] += A[(i)+(672)] * values[12];
C[(i)+(420)] += A[(i)+(504)] * values[13];
C[(i)+(420)] += A[(i)+(672)] * values[14];
C[(i)+(504)] += A[(i)+(0)] * values[15];
C[(i)+(504)] += A[(i)+(252)] * values[16];
C[(i)+(504)] += A[(i)+(420)] * values[17];
C[(i)+(588)] += A[(i)+(84)] * values[18];
C[(i)+(588)] += A[(i)+(252)] * values[19];
C[(i)+(588)] += A[(i)+(336)] * values[20];
C[(i)+(672)] += A[(i)+(168)] * values[21];
C[(i)+(672)] += A[(i)+(336)] * values[22];
C[(i)+(672)] += A[(i)+(420)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 4032;
#endif

}

void dsparse_fM1DivM_m84_n9_k84_ldAna7_ldB84_ldC84_beta0_pfsigonly(const double* values, const double* B, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma nounroll_and_jam
for (unsigned int i = 0; i < 9; i++)
{
  #pragma simd
  for (unsigned int m = 0; m < 84; m++) {
    C[(i*84)+m] = 0.0;
  }
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b0 = _mm256_broadcast_sd(&B[(i*84)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b0 = _mm_loaddup_pd(&B[(i*84)+0]);
#endif
__m128d c0_0 = _mm_load_sd(&C[(i*84)+0]);
__m128d a0_0 = _mm_load_sd(&values[0]);
#if defined(__SSE3__) && defined(__AVX__)
c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, b0));
#endif
_mm_store_sd(&C[(i*84)+0], c0_0);
__m128d c0_1 = _mm_load_sd(&C[(i*84)+3]);
__m128d a0_1 = _mm_load_sd(&values[1]);
#if defined(__SSE3__) && defined(__AVX__)
c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, b0));
#endif
_mm_store_sd(&C[(i*84)+3], c0_1);
__m128d c0_2 = _mm_load_sd(&C[(i*84)+9]);
__m128d a0_2 = _mm_load_sd(&values[2]);
#if defined(__SSE3__) && defined(__AVX__)
c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, b0));
#endif
_mm_store_sd(&C[(i*84)+9], c0_2);
__m128d c0_3 = _mm_load_sd(&C[(i*84)+19]);
__m128d a0_3 = _mm_load_sd(&values[3]);
#if defined(__SSE3__) && defined(__AVX__)
c0_3 = _mm_add_sd(c0_3, _mm_mul_sd(a0_3, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_3 = _mm_add_sd(c0_3, _mm_mul_sd(a0_3, b0));
#endif
_mm_store_sd(&C[(i*84)+19], c0_3);
__m128d c0_4 = _mm_load_sd(&C[(i*84)+34]);
__m128d a0_4 = _mm_load_sd(&values[4]);
#if defined(__SSE3__) && defined(__AVX__)
c0_4 = _mm_add_sd(c0_4, _mm_mul_sd(a0_4, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_4 = _mm_add_sd(c0_4, _mm_mul_sd(a0_4, b0));
#endif
_mm_store_sd(&C[(i*84)+34], c0_4);
__m128d c0_5 = _mm_load_sd(&C[(i*84)+55]);
__m128d a0_5 = _mm_load_sd(&values[5]);
#if defined(__SSE3__) && defined(__AVX__)
c0_5 = _mm_add_sd(c0_5, _mm_mul_sd(a0_5, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_5 = _mm_add_sd(c0_5, _mm_mul_sd(a0_5, b0));
#endif
_mm_store_sd(&C[(i*84)+55], c0_5);
__m128d c0_6 = _mm_load_sd(&C[(i*84)+83]);
__m128d a0_6 = _mm_load_sd(&values[6]);
#if defined(__SSE3__) && defined(__AVX__)
c0_6 = _mm_add_sd(c0_6, _mm_mul_sd(a0_6, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_6 = _mm_add_sd(c0_6, _mm_mul_sd(a0_6, b0));
#endif
_mm_store_sd(&C[(i*84)+83], c0_6);
#else
C[(i*84)+0] += values[0] * B[(i*84)+0];
C[(i*84)+3] += values[1] * B[(i*84)+0];
C[(i*84)+9] += values[2] * B[(i*84)+0];
C[(i*84)+19] += values[3] * B[(i*84)+0];
C[(i*84)+34] += values[4] * B[(i*84)+0];
C[(i*84)+55] += values[5] * B[(i*84)+0];
C[(i*84)+83] += values[6] * B[(i*84)+0];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b1 = _mm256_broadcast_sd(&B[(i*84)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b1 = _mm_loaddup_pd(&B[(i*84)+1]);
#endif
__m128d c1_0 = _mm_load_sd(&C[(i*84)+1]);
__m128d a1_0 = _mm_load_sd(&values[7]);
#if defined(__SSE3__) && defined(__AVX__)
c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, b1));
#endif
_mm_store_sd(&C[(i*84)+1], c1_0);
__m128d c1_1 = _mm_load_sd(&C[(i*84)+7]);
__m128d a1_1 = _mm_load_sd(&values[8]);
#if defined(__SSE3__) && defined(__AVX__)
c1_1 = _mm_add_sd(c1_1, _mm_mul_sd(a1_1, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_1 = _mm_add_sd(c1_1, _mm_mul_sd(a1_1, b1));
#endif
_mm_store_sd(&C[(i*84)+7], c1_1);
__m128d c1_2 = _mm_load_sd(&C[(i*84)+17]);
__m128d a1_2 = _mm_load_sd(&values[9]);
#if defined(__SSE3__) && defined(__AVX__)
c1_2 = _mm_add_sd(c1_2, _mm_mul_sd(a1_2, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_2 = _mm_add_sd(c1_2, _mm_mul_sd(a1_2, b1));
#endif
_mm_store_sd(&C[(i*84)+17], c1_2);
__m128d c1_3 = _mm_load_sd(&C[(i*84)+32]);
__m128d a1_3 = _mm_load_sd(&values[10]);
#if defined(__SSE3__) && defined(__AVX__)
c1_3 = _mm_add_sd(c1_3, _mm_mul_sd(a1_3, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_3 = _mm_add_sd(c1_3, _mm_mul_sd(a1_3, b1));
#endif
_mm_store_sd(&C[(i*84)+32], c1_3);
__m128d c1_4 = _mm_load_sd(&C[(i*84)+53]);
__m128d a1_4 = _mm_load_sd(&values[11]);
#if defined(__SSE3__) && defined(__AVX__)
c1_4 = _mm_add_sd(c1_4, _mm_mul_sd(a1_4, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_4 = _mm_add_sd(c1_4, _mm_mul_sd(a1_4, b1));
#endif
_mm_store_sd(&C[(i*84)+53], c1_4);
__m128d c1_5 = _mm_load_sd(&C[(i*84)+81]);
__m128d a1_5 = _mm_load_sd(&values[12]);
#if defined(__SSE3__) && defined(__AVX__)
c1_5 = _mm_add_sd(c1_5, _mm_mul_sd(a1_5, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_5 = _mm_add_sd(c1_5, _mm_mul_sd(a1_5, b1));
#endif
_mm_store_sd(&C[(i*84)+81], c1_5);
#else
C[(i*84)+1] += values[7] * B[(i*84)+1];
C[(i*84)+7] += values[8] * B[(i*84)+1];
C[(i*84)+17] += values[9] * B[(i*84)+1];
C[(i*84)+32] += values[10] * B[(i*84)+1];
C[(i*84)+53] += values[11] * B[(i*84)+1];
C[(i*84)+81] += values[12] * B[(i*84)+1];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b2 = _mm256_broadcast_sd(&B[(i*84)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b2 = _mm_loaddup_pd(&B[(i*84)+2]);
#endif
__m128d c2_0 = _mm_load_sd(&C[(i*84)+2]);
__m128d a2_0 = _mm_load_sd(&values[13]);
#if defined(__SSE3__) && defined(__AVX__)
c2_0 = _mm_add_sd(c2_0, _mm_mul_sd(a2_0, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_0 = _mm_add_sd(c2_0, _mm_mul_sd(a2_0, b2));
#endif
_mm_store_sd(&C[(i*84)+2], c2_0);
__m128d c2_1 = _mm_load_sd(&C[(i*84)+8]);
__m128d a2_1 = _mm_load_sd(&values[14]);
#if defined(__SSE3__) && defined(__AVX__)
c2_1 = _mm_add_sd(c2_1, _mm_mul_sd(a2_1, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_1 = _mm_add_sd(c2_1, _mm_mul_sd(a2_1, b2));
#endif
_mm_store_sd(&C[(i*84)+8], c2_1);
__m128d c2_2 = _mm_load_sd(&C[(i*84)+18]);
__m128d a2_2 = _mm_load_sd(&values[15]);
#if defined(__SSE3__) && defined(__AVX__)
c2_2 = _mm_add_sd(c2_2, _mm_mul_sd(a2_2, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_2 = _mm_add_sd(c2_2, _mm_mul_sd(a2_2, b2));
#endif
_mm_store_sd(&C[(i*84)+18], c2_2);
__m128d c2_3 = _mm_load_sd(&C[(i*84)+33]);
__m128d a2_3 = _mm_load_sd(&values[16]);
#if defined(__SSE3__) && defined(__AVX__)
c2_3 = _mm_add_sd(c2_3, _mm_mul_sd(a2_3, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_3 = _mm_add_sd(c2_3, _mm_mul_sd(a2_3, b2));
#endif
_mm_store_sd(&C[(i*84)+33], c2_3);
__m128d c2_4 = _mm_load_sd(&C[(i*84)+54]);
__m128d a2_4 = _mm_load_sd(&values[17]);
#if defined(__SSE3__) && defined(__AVX__)
c2_4 = _mm_add_sd(c2_4, _mm_mul_sd(a2_4, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_4 = _mm_add_sd(c2_4, _mm_mul_sd(a2_4, b2));
#endif
_mm_store_sd(&C[(i*84)+54], c2_4);
__m128d c2_5 = _mm_load_sd(&C[(i*84)+82]);
__m128d a2_5 = _mm_load_sd(&values[18]);
#if defined(__SSE3__) && defined(__AVX__)
c2_5 = _mm_add_sd(c2_5, _mm_mul_sd(a2_5, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_5 = _mm_add_sd(c2_5, _mm_mul_sd(a2_5, b2));
#endif
_mm_store_sd(&C[(i*84)+82], c2_5);
#else
C[(i*84)+2] += values[13] * B[(i*84)+2];
C[(i*84)+8] += values[14] * B[(i*84)+2];
C[(i*84)+18] += values[15] * B[(i*84)+2];
C[(i*84)+33] += values[16] * B[(i*84)+2];
C[(i*84)+54] += values[17] * B[(i*84)+2];
C[(i*84)+82] += values[18] * B[(i*84)+2];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b3 = _mm256_broadcast_sd(&B[(i*84)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b3 = _mm_loaddup_pd(&B[(i*84)+3]);
#endif
__m128d c3_0 = _mm_load_sd(&C[(i*84)+0]);
__m128d a3_0 = _mm_load_sd(&values[19]);
#if defined(__SSE3__) && defined(__AVX__)
c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, b3));
#endif
_mm_store_sd(&C[(i*84)+0], c3_0);
__m128d c3_1 = _mm_load_sd(&C[(i*84)+3]);
__m128d a3_1 = _mm_load_sd(&values[20]);
#if defined(__SSE3__) && defined(__AVX__)
c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, b3));
#endif
_mm_store_sd(&C[(i*84)+3], c3_1);
__m128d c3_2 = _mm_load_sd(&C[(i*84)+9]);
__m128d a3_2 = _mm_load_sd(&values[21]);
#if defined(__SSE3__) && defined(__AVX__)
c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, b3));
#endif
_mm_store_sd(&C[(i*84)+9], c3_2);
__m128d c3_3 = _mm_load_sd(&C[(i*84)+19]);
__m128d a3_3 = _mm_load_sd(&values[22]);
#if defined(__SSE3__) && defined(__AVX__)
c3_3 = _mm_add_sd(c3_3, _mm_mul_sd(a3_3, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_3 = _mm_add_sd(c3_3, _mm_mul_sd(a3_3, b3));
#endif
_mm_store_sd(&C[(i*84)+19], c3_3);
__m128d c3_4 = _mm_load_sd(&C[(i*84)+34]);
__m128d a3_4 = _mm_load_sd(&values[23]);
#if defined(__SSE3__) && defined(__AVX__)
c3_4 = _mm_add_sd(c3_4, _mm_mul_sd(a3_4, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_4 = _mm_add_sd(c3_4, _mm_mul_sd(a3_4, b3));
#endif
_mm_store_sd(&C[(i*84)+34], c3_4);
__m128d c3_5 = _mm_load_sd(&C[(i*84)+55]);
__m128d a3_5 = _mm_load_sd(&values[24]);
#if defined(__SSE3__) && defined(__AVX__)
c3_5 = _mm_add_sd(c3_5, _mm_mul_sd(a3_5, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_5 = _mm_add_sd(c3_5, _mm_mul_sd(a3_5, b3));
#endif
_mm_store_sd(&C[(i*84)+55], c3_5);
__m128d c3_6 = _mm_load_sd(&C[(i*84)+83]);
__m128d a3_6 = _mm_load_sd(&values[25]);
#if defined(__SSE3__) && defined(__AVX__)
c3_6 = _mm_add_sd(c3_6, _mm_mul_sd(a3_6, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_6 = _mm_add_sd(c3_6, _mm_mul_sd(a3_6, b3));
#endif
_mm_store_sd(&C[(i*84)+83], c3_6);
#else
C[(i*84)+0] += values[19] * B[(i*84)+3];
C[(i*84)+3] += values[20] * B[(i*84)+3];
C[(i*84)+9] += values[21] * B[(i*84)+3];
C[(i*84)+19] += values[22] * B[(i*84)+3];
C[(i*84)+34] += values[23] * B[(i*84)+3];
C[(i*84)+55] += values[24] * B[(i*84)+3];
C[(i*84)+83] += values[25] * B[(i*84)+3];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b4 = _mm256_broadcast_sd(&B[(i*84)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b4 = _mm_loaddup_pd(&B[(i*84)+4]);
#endif
__m128d c4_0 = _mm_load_sd(&C[(i*84)+4]);
__m128d a4_0 = _mm_load_sd(&values[26]);
#if defined(__SSE3__) && defined(__AVX__)
c4_0 = _mm_add_sd(c4_0, _mm_mul_sd(a4_0, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_0 = _mm_add_sd(c4_0, _mm_mul_sd(a4_0, b4));
#endif
_mm_store_sd(&C[(i*84)+4], c4_0);
__m128d c4_1 = _mm_load_sd(&C[(i*84)+14]);
__m128d a4_1 = _mm_load_sd(&values[27]);
#if defined(__SSE3__) && defined(__AVX__)
c4_1 = _mm_add_sd(c4_1, _mm_mul_sd(a4_1, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_1 = _mm_add_sd(c4_1, _mm_mul_sd(a4_1, b4));
#endif
_mm_store_sd(&C[(i*84)+14], c4_1);
__m128d c4_2 = _mm_load_sd(&C[(i*84)+29]);
__m128d a4_2 = _mm_load_sd(&values[28]);
#if defined(__SSE3__) && defined(__AVX__)
c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, b4));
#endif
_mm_store_sd(&C[(i*84)+29], c4_2);
__m128d c4_3 = _mm_load_sd(&C[(i*84)+50]);
__m128d a4_3 = _mm_load_sd(&values[29]);
#if defined(__SSE3__) && defined(__AVX__)
c4_3 = _mm_add_sd(c4_3, _mm_mul_sd(a4_3, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_3 = _mm_add_sd(c4_3, _mm_mul_sd(a4_3, b4));
#endif
_mm_store_sd(&C[(i*84)+50], c4_3);
__m128d c4_4 = _mm_load_sd(&C[(i*84)+78]);
__m128d a4_4 = _mm_load_sd(&values[30]);
#if defined(__SSE3__) && defined(__AVX__)
c4_4 = _mm_add_sd(c4_4, _mm_mul_sd(a4_4, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_4 = _mm_add_sd(c4_4, _mm_mul_sd(a4_4, b4));
#endif
_mm_store_sd(&C[(i*84)+78], c4_4);
#else
C[(i*84)+4] += values[26] * B[(i*84)+4];
C[(i*84)+14] += values[27] * B[(i*84)+4];
C[(i*84)+29] += values[28] * B[(i*84)+4];
C[(i*84)+50] += values[29] * B[(i*84)+4];
C[(i*84)+78] += values[30] * B[(i*84)+4];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b5 = _mm256_broadcast_sd(&B[(i*84)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b5 = _mm_loaddup_pd(&B[(i*84)+5]);
#endif
__m128d c5_0 = _mm_load_sd(&C[(i*84)+5]);
__m128d a5_0 = _mm_load_sd(&values[31]);
#if defined(__SSE3__) && defined(__AVX__)
c5_0 = _mm_add_sd(c5_0, _mm_mul_sd(a5_0, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_0 = _mm_add_sd(c5_0, _mm_mul_sd(a5_0, b5));
#endif
_mm_store_sd(&C[(i*84)+5], c5_0);
__m128d c5_1 = _mm_load_sd(&C[(i*84)+15]);
__m128d a5_1 = _mm_load_sd(&values[32]);
#if defined(__SSE3__) && defined(__AVX__)
c5_1 = _mm_add_sd(c5_1, _mm_mul_sd(a5_1, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_1 = _mm_add_sd(c5_1, _mm_mul_sd(a5_1, b5));
#endif
_mm_store_sd(&C[(i*84)+15], c5_1);
__m128d c5_2 = _mm_load_sd(&C[(i*84)+30]);
__m128d a5_2 = _mm_load_sd(&values[33]);
#if defined(__SSE3__) && defined(__AVX__)
c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, b5));
#endif
_mm_store_sd(&C[(i*84)+30], c5_2);
__m128d c5_3 = _mm_load_sd(&C[(i*84)+51]);
__m128d a5_3 = _mm_load_sd(&values[34]);
#if defined(__SSE3__) && defined(__AVX__)
c5_3 = _mm_add_sd(c5_3, _mm_mul_sd(a5_3, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_3 = _mm_add_sd(c5_3, _mm_mul_sd(a5_3, b5));
#endif
_mm_store_sd(&C[(i*84)+51], c5_3);
__m128d c5_4 = _mm_load_sd(&C[(i*84)+79]);
__m128d a5_4 = _mm_load_sd(&values[35]);
#if defined(__SSE3__) && defined(__AVX__)
c5_4 = _mm_add_sd(c5_4, _mm_mul_sd(a5_4, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_4 = _mm_add_sd(c5_4, _mm_mul_sd(a5_4, b5));
#endif
_mm_store_sd(&C[(i*84)+79], c5_4);
#else
C[(i*84)+5] += values[31] * B[(i*84)+5];
C[(i*84)+15] += values[32] * B[(i*84)+5];
C[(i*84)+30] += values[33] * B[(i*84)+5];
C[(i*84)+51] += values[34] * B[(i*84)+5];
C[(i*84)+79] += values[35] * B[(i*84)+5];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b6 = _mm256_broadcast_sd(&B[(i*84)+6]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b6 = _mm_loaddup_pd(&B[(i*84)+6]);
#endif
__m128d c6_0 = _mm_load_sd(&C[(i*84)+6]);
__m128d a6_0 = _mm_load_sd(&values[36]);
#if defined(__SSE3__) && defined(__AVX__)
c6_0 = _mm_add_sd(c6_0, _mm_mul_sd(a6_0, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_0 = _mm_add_sd(c6_0, _mm_mul_sd(a6_0, b6));
#endif
_mm_store_sd(&C[(i*84)+6], c6_0);
__m128d c6_1 = _mm_load_sd(&C[(i*84)+16]);
__m128d a6_1 = _mm_load_sd(&values[37]);
#if defined(__SSE3__) && defined(__AVX__)
c6_1 = _mm_add_sd(c6_1, _mm_mul_sd(a6_1, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_1 = _mm_add_sd(c6_1, _mm_mul_sd(a6_1, b6));
#endif
_mm_store_sd(&C[(i*84)+16], c6_1);
__m128d c6_2 = _mm_load_sd(&C[(i*84)+31]);
__m128d a6_2 = _mm_load_sd(&values[38]);
#if defined(__SSE3__) && defined(__AVX__)
c6_2 = _mm_add_sd(c6_2, _mm_mul_sd(a6_2, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_2 = _mm_add_sd(c6_2, _mm_mul_sd(a6_2, b6));
#endif
_mm_store_sd(&C[(i*84)+31], c6_2);
__m128d c6_3 = _mm_load_sd(&C[(i*84)+52]);
__m128d a6_3 = _mm_load_sd(&values[39]);
#if defined(__SSE3__) && defined(__AVX__)
c6_3 = _mm_add_sd(c6_3, _mm_mul_sd(a6_3, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_3 = _mm_add_sd(c6_3, _mm_mul_sd(a6_3, b6));
#endif
_mm_store_sd(&C[(i*84)+52], c6_3);
__m128d c6_4 = _mm_load_sd(&C[(i*84)+80]);
__m128d a6_4 = _mm_load_sd(&values[40]);
#if defined(__SSE3__) && defined(__AVX__)
c6_4 = _mm_add_sd(c6_4, _mm_mul_sd(a6_4, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_4 = _mm_add_sd(c6_4, _mm_mul_sd(a6_4, b6));
#endif
_mm_store_sd(&C[(i*84)+80], c6_4);
#else
C[(i*84)+6] += values[36] * B[(i*84)+6];
C[(i*84)+16] += values[37] * B[(i*84)+6];
C[(i*84)+31] += values[38] * B[(i*84)+6];
C[(i*84)+52] += values[39] * B[(i*84)+6];
C[(i*84)+80] += values[40] * B[(i*84)+6];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b7 = _mm256_broadcast_sd(&B[(i*84)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b7 = _mm_loaddup_pd(&B[(i*84)+7]);
#endif
__m128d c7_0 = _mm_load_sd(&C[(i*84)+1]);
__m128d a7_0 = _mm_load_sd(&values[41]);
#if defined(__SSE3__) && defined(__AVX__)
c7_0 = _mm_add_sd(c7_0, _mm_mul_sd(a7_0, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_0 = _mm_add_sd(c7_0, _mm_mul_sd(a7_0, b7));
#endif
_mm_store_sd(&C[(i*84)+1], c7_0);
__m128d c7_1 = _mm_load_sd(&C[(i*84)+7]);
__m128d a7_1 = _mm_load_sd(&values[42]);
#if defined(__SSE3__) && defined(__AVX__)
c7_1 = _mm_add_sd(c7_1, _mm_mul_sd(a7_1, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_1 = _mm_add_sd(c7_1, _mm_mul_sd(a7_1, b7));
#endif
_mm_store_sd(&C[(i*84)+7], c7_1);
__m128d c7_2 = _mm_load_sd(&C[(i*84)+17]);
__m128d a7_2 = _mm_load_sd(&values[43]);
#if defined(__SSE3__) && defined(__AVX__)
c7_2 = _mm_add_sd(c7_2, _mm_mul_sd(a7_2, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_2 = _mm_add_sd(c7_2, _mm_mul_sd(a7_2, b7));
#endif
_mm_store_sd(&C[(i*84)+17], c7_2);
__m128d c7_3 = _mm_load_sd(&C[(i*84)+32]);
__m128d a7_3 = _mm_load_sd(&values[44]);
#if defined(__SSE3__) && defined(__AVX__)
c7_3 = _mm_add_sd(c7_3, _mm_mul_sd(a7_3, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_3 = _mm_add_sd(c7_3, _mm_mul_sd(a7_3, b7));
#endif
_mm_store_sd(&C[(i*84)+32], c7_3);
__m128d c7_4 = _mm_load_sd(&C[(i*84)+53]);
__m128d a7_4 = _mm_load_sd(&values[45]);
#if defined(__SSE3__) && defined(__AVX__)
c7_4 = _mm_add_sd(c7_4, _mm_mul_sd(a7_4, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_4 = _mm_add_sd(c7_4, _mm_mul_sd(a7_4, b7));
#endif
_mm_store_sd(&C[(i*84)+53], c7_4);
__m128d c7_5 = _mm_load_sd(&C[(i*84)+81]);
__m128d a7_5 = _mm_load_sd(&values[46]);
#if defined(__SSE3__) && defined(__AVX__)
c7_5 = _mm_add_sd(c7_5, _mm_mul_sd(a7_5, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_5 = _mm_add_sd(c7_5, _mm_mul_sd(a7_5, b7));
#endif
_mm_store_sd(&C[(i*84)+81], c7_5);
#else
C[(i*84)+1] += values[41] * B[(i*84)+7];
C[(i*84)+7] += values[42] * B[(i*84)+7];
C[(i*84)+17] += values[43] * B[(i*84)+7];
C[(i*84)+32] += values[44] * B[(i*84)+7];
C[(i*84)+53] += values[45] * B[(i*84)+7];
C[(i*84)+81] += values[46] * B[(i*84)+7];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b8 = _mm256_broadcast_sd(&B[(i*84)+8]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b8 = _mm_loaddup_pd(&B[(i*84)+8]);
#endif
__m128d c8_0 = _mm_load_sd(&C[(i*84)+2]);
__m128d a8_0 = _mm_load_sd(&values[47]);
#if defined(__SSE3__) && defined(__AVX__)
c8_0 = _mm_add_sd(c8_0, _mm_mul_sd(a8_0, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_0 = _mm_add_sd(c8_0, _mm_mul_sd(a8_0, b8));
#endif
_mm_store_sd(&C[(i*84)+2], c8_0);
__m128d c8_1 = _mm_load_sd(&C[(i*84)+8]);
__m128d a8_1 = _mm_load_sd(&values[48]);
#if defined(__SSE3__) && defined(__AVX__)
c8_1 = _mm_add_sd(c8_1, _mm_mul_sd(a8_1, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_1 = _mm_add_sd(c8_1, _mm_mul_sd(a8_1, b8));
#endif
_mm_store_sd(&C[(i*84)+8], c8_1);
__m128d c8_2 = _mm_load_sd(&C[(i*84)+18]);
__m128d a8_2 = _mm_load_sd(&values[49]);
#if defined(__SSE3__) && defined(__AVX__)
c8_2 = _mm_add_sd(c8_2, _mm_mul_sd(a8_2, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_2 = _mm_add_sd(c8_2, _mm_mul_sd(a8_2, b8));
#endif
_mm_store_sd(&C[(i*84)+18], c8_2);
__m128d c8_3 = _mm_load_sd(&C[(i*84)+33]);
__m128d a8_3 = _mm_load_sd(&values[50]);
#if defined(__SSE3__) && defined(__AVX__)
c8_3 = _mm_add_sd(c8_3, _mm_mul_sd(a8_3, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_3 = _mm_add_sd(c8_3, _mm_mul_sd(a8_3, b8));
#endif
_mm_store_sd(&C[(i*84)+33], c8_3);
__m128d c8_4 = _mm_load_sd(&C[(i*84)+54]);
__m128d a8_4 = _mm_load_sd(&values[51]);
#if defined(__SSE3__) && defined(__AVX__)
c8_4 = _mm_add_sd(c8_4, _mm_mul_sd(a8_4, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_4 = _mm_add_sd(c8_4, _mm_mul_sd(a8_4, b8));
#endif
_mm_store_sd(&C[(i*84)+54], c8_4);
__m128d c8_5 = _mm_load_sd(&C[(i*84)+82]);
__m128d a8_5 = _mm_load_sd(&values[52]);
#if defined(__SSE3__) && defined(__AVX__)
c8_5 = _mm_add_sd(c8_5, _mm_mul_sd(a8_5, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_5 = _mm_add_sd(c8_5, _mm_mul_sd(a8_5, b8));
#endif
_mm_store_sd(&C[(i*84)+82], c8_5);
#else
C[(i*84)+2] += values[47] * B[(i*84)+8];
C[(i*84)+8] += values[48] * B[(i*84)+8];
C[(i*84)+18] += values[49] * B[(i*84)+8];
C[(i*84)+33] += values[50] * B[(i*84)+8];
C[(i*84)+54] += values[51] * B[(i*84)+8];
C[(i*84)+82] += values[52] * B[(i*84)+8];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b9 = _mm256_broadcast_sd(&B[(i*84)+9]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b9 = _mm_loaddup_pd(&B[(i*84)+9]);
#endif
__m128d c9_0 = _mm_load_sd(&C[(i*84)+0]);
__m128d a9_0 = _mm_load_sd(&values[53]);
#if defined(__SSE3__) && defined(__AVX__)
c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, b9));
#endif
_mm_store_sd(&C[(i*84)+0], c9_0);
__m128d c9_1 = _mm_load_sd(&C[(i*84)+3]);
__m128d a9_1 = _mm_load_sd(&values[54]);
#if defined(__SSE3__) && defined(__AVX__)
c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, b9));
#endif
_mm_store_sd(&C[(i*84)+3], c9_1);
__m128d c9_2 = _mm_load_sd(&C[(i*84)+9]);
__m128d a9_2 = _mm_load_sd(&values[55]);
#if defined(__SSE3__) && defined(__AVX__)
c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, b9));
#endif
_mm_store_sd(&C[(i*84)+9], c9_2);
__m128d c9_3 = _mm_load_sd(&C[(i*84)+19]);
__m128d a9_3 = _mm_load_sd(&values[56]);
#if defined(__SSE3__) && defined(__AVX__)
c9_3 = _mm_add_sd(c9_3, _mm_mul_sd(a9_3, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_3 = _mm_add_sd(c9_3, _mm_mul_sd(a9_3, b9));
#endif
_mm_store_sd(&C[(i*84)+19], c9_3);
__m128d c9_4 = _mm_load_sd(&C[(i*84)+34]);
__m128d a9_4 = _mm_load_sd(&values[57]);
#if defined(__SSE3__) && defined(__AVX__)
c9_4 = _mm_add_sd(c9_4, _mm_mul_sd(a9_4, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_4 = _mm_add_sd(c9_4, _mm_mul_sd(a9_4, b9));
#endif
_mm_store_sd(&C[(i*84)+34], c9_4);
__m128d c9_5 = _mm_load_sd(&C[(i*84)+55]);
__m128d a9_5 = _mm_load_sd(&values[58]);
#if defined(__SSE3__) && defined(__AVX__)
c9_5 = _mm_add_sd(c9_5, _mm_mul_sd(a9_5, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_5 = _mm_add_sd(c9_5, _mm_mul_sd(a9_5, b9));
#endif
_mm_store_sd(&C[(i*84)+55], c9_5);
__m128d c9_6 = _mm_load_sd(&C[(i*84)+83]);
__m128d a9_6 = _mm_load_sd(&values[59]);
#if defined(__SSE3__) && defined(__AVX__)
c9_6 = _mm_add_sd(c9_6, _mm_mul_sd(a9_6, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_6 = _mm_add_sd(c9_6, _mm_mul_sd(a9_6, b9));
#endif
_mm_store_sd(&C[(i*84)+83], c9_6);
#else
C[(i*84)+0] += values[53] * B[(i*84)+9];
C[(i*84)+3] += values[54] * B[(i*84)+9];
C[(i*84)+9] += values[55] * B[(i*84)+9];
C[(i*84)+19] += values[56] * B[(i*84)+9];
C[(i*84)+34] += values[57] * B[(i*84)+9];
C[(i*84)+55] += values[58] * B[(i*84)+9];
C[(i*84)+83] += values[59] * B[(i*84)+9];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b10 = _mm256_broadcast_sd(&B[(i*84)+10]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b10 = _mm_loaddup_pd(&B[(i*84)+10]);
#endif
__m128d c10_0 = _mm_load_sd(&C[(i*84)+10]);
__m128d a10_0 = _mm_load_sd(&values[60]);
#if defined(__SSE3__) && defined(__AVX__)
c10_0 = _mm_add_sd(c10_0, _mm_mul_sd(a10_0, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c10_0 = _mm_add_sd(c10_0, _mm_mul_sd(a10_0, b10));
#endif
_mm_store_sd(&C[(i*84)+10], c10_0);
__m128d c10_1 = _mm_load_sd(&C[(i*84)+25]);
__m128d a10_1 = _mm_load_sd(&values[61]);
#if defined(__SSE3__) && defined(__AVX__)
c10_1 = _mm_add_sd(c10_1, _mm_mul_sd(a10_1, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c10_1 = _mm_add_sd(c10_1, _mm_mul_sd(a10_1, b10));
#endif
_mm_store_sd(&C[(i*84)+25], c10_1);
__m128d c10_2 = _mm_load_sd(&C[(i*84)+46]);
__m128d a10_2 = _mm_load_sd(&values[62]);
#if defined(__SSE3__) && defined(__AVX__)
c10_2 = _mm_add_sd(c10_2, _mm_mul_sd(a10_2, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c10_2 = _mm_add_sd(c10_2, _mm_mul_sd(a10_2, b10));
#endif
_mm_store_sd(&C[(i*84)+46], c10_2);
__m128d c10_3 = _mm_load_sd(&C[(i*84)+74]);
__m128d a10_3 = _mm_load_sd(&values[63]);
#if defined(__SSE3__) && defined(__AVX__)
c10_3 = _mm_add_sd(c10_3, _mm_mul_sd(a10_3, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c10_3 = _mm_add_sd(c10_3, _mm_mul_sd(a10_3, b10));
#endif
_mm_store_sd(&C[(i*84)+74], c10_3);
#else
C[(i*84)+10] += values[60] * B[(i*84)+10];
C[(i*84)+25] += values[61] * B[(i*84)+10];
C[(i*84)+46] += values[62] * B[(i*84)+10];
C[(i*84)+74] += values[63] * B[(i*84)+10];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b11 = _mm256_broadcast_sd(&B[(i*84)+11]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b11 = _mm_loaddup_pd(&B[(i*84)+11]);
#endif
__m128d c11_0 = _mm_load_sd(&C[(i*84)+11]);
__m128d a11_0 = _mm_load_sd(&values[64]);
#if defined(__SSE3__) && defined(__AVX__)
c11_0 = _mm_add_sd(c11_0, _mm_mul_sd(a11_0, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c11_0 = _mm_add_sd(c11_0, _mm_mul_sd(a11_0, b11));
#endif
_mm_store_sd(&C[(i*84)+11], c11_0);
__m128d c11_1 = _mm_load_sd(&C[(i*84)+26]);
__m128d a11_1 = _mm_load_sd(&values[65]);
#if defined(__SSE3__) && defined(__AVX__)
c11_1 = _mm_add_sd(c11_1, _mm_mul_sd(a11_1, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c11_1 = _mm_add_sd(c11_1, _mm_mul_sd(a11_1, b11));
#endif
_mm_store_sd(&C[(i*84)+26], c11_1);
__m128d c11_2 = _mm_load_sd(&C[(i*84)+47]);
__m128d a11_2 = _mm_load_sd(&values[66]);
#if defined(__SSE3__) && defined(__AVX__)
c11_2 = _mm_add_sd(c11_2, _mm_mul_sd(a11_2, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c11_2 = _mm_add_sd(c11_2, _mm_mul_sd(a11_2, b11));
#endif
_mm_store_sd(&C[(i*84)+47], c11_2);
__m128d c11_3 = _mm_load_sd(&C[(i*84)+75]);
__m128d a11_3 = _mm_load_sd(&values[67]);
#if defined(__SSE3__) && defined(__AVX__)
c11_3 = _mm_add_sd(c11_3, _mm_mul_sd(a11_3, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c11_3 = _mm_add_sd(c11_3, _mm_mul_sd(a11_3, b11));
#endif
_mm_store_sd(&C[(i*84)+75], c11_3);
#else
C[(i*84)+11] += values[64] * B[(i*84)+11];
C[(i*84)+26] += values[65] * B[(i*84)+11];
C[(i*84)+47] += values[66] * B[(i*84)+11];
C[(i*84)+75] += values[67] * B[(i*84)+11];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b12 = _mm256_broadcast_sd(&B[(i*84)+12]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b12 = _mm_loaddup_pd(&B[(i*84)+12]);
#endif
__m128d c12_0 = _mm_load_sd(&C[(i*84)+12]);
__m128d a12_0 = _mm_load_sd(&values[68]);
#if defined(__SSE3__) && defined(__AVX__)
c12_0 = _mm_add_sd(c12_0, _mm_mul_sd(a12_0, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c12_0 = _mm_add_sd(c12_0, _mm_mul_sd(a12_0, b12));
#endif
_mm_store_sd(&C[(i*84)+12], c12_0);
__m128d c12_1 = _mm_load_sd(&C[(i*84)+27]);
__m128d a12_1 = _mm_load_sd(&values[69]);
#if defined(__SSE3__) && defined(__AVX__)
c12_1 = _mm_add_sd(c12_1, _mm_mul_sd(a12_1, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c12_1 = _mm_add_sd(c12_1, _mm_mul_sd(a12_1, b12));
#endif
_mm_store_sd(&C[(i*84)+27], c12_1);
__m128d c12_2 = _mm_load_sd(&C[(i*84)+48]);
__m128d a12_2 = _mm_load_sd(&values[70]);
#if defined(__SSE3__) && defined(__AVX__)
c12_2 = _mm_add_sd(c12_2, _mm_mul_sd(a12_2, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c12_2 = _mm_add_sd(c12_2, _mm_mul_sd(a12_2, b12));
#endif
_mm_store_sd(&C[(i*84)+48], c12_2);
__m128d c12_3 = _mm_load_sd(&C[(i*84)+76]);
__m128d a12_3 = _mm_load_sd(&values[71]);
#if defined(__SSE3__) && defined(__AVX__)
c12_3 = _mm_add_sd(c12_3, _mm_mul_sd(a12_3, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c12_3 = _mm_add_sd(c12_3, _mm_mul_sd(a12_3, b12));
#endif
_mm_store_sd(&C[(i*84)+76], c12_3);
#else
C[(i*84)+12] += values[68] * B[(i*84)+12];
C[(i*84)+27] += values[69] * B[(i*84)+12];
C[(i*84)+48] += values[70] * B[(i*84)+12];
C[(i*84)+76] += values[71] * B[(i*84)+12];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b13 = _mm256_broadcast_sd(&B[(i*84)+13]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b13 = _mm_loaddup_pd(&B[(i*84)+13]);
#endif
__m128d c13_0 = _mm_load_sd(&C[(i*84)+13]);
__m128d a13_0 = _mm_load_sd(&values[72]);
#if defined(__SSE3__) && defined(__AVX__)
c13_0 = _mm_add_sd(c13_0, _mm_mul_sd(a13_0, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c13_0 = _mm_add_sd(c13_0, _mm_mul_sd(a13_0, b13));
#endif
_mm_store_sd(&C[(i*84)+13], c13_0);
__m128d c13_1 = _mm_load_sd(&C[(i*84)+28]);
__m128d a13_1 = _mm_load_sd(&values[73]);
#if defined(__SSE3__) && defined(__AVX__)
c13_1 = _mm_add_sd(c13_1, _mm_mul_sd(a13_1, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c13_1 = _mm_add_sd(c13_1, _mm_mul_sd(a13_1, b13));
#endif
_mm_store_sd(&C[(i*84)+28], c13_1);
__m128d c13_2 = _mm_load_sd(&C[(i*84)+49]);
__m128d a13_2 = _mm_load_sd(&values[74]);
#if defined(__SSE3__) && defined(__AVX__)
c13_2 = _mm_add_sd(c13_2, _mm_mul_sd(a13_2, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c13_2 = _mm_add_sd(c13_2, _mm_mul_sd(a13_2, b13));
#endif
_mm_store_sd(&C[(i*84)+49], c13_2);
__m128d c13_3 = _mm_load_sd(&C[(i*84)+77]);
__m128d a13_3 = _mm_load_sd(&values[75]);
#if defined(__SSE3__) && defined(__AVX__)
c13_3 = _mm_add_sd(c13_3, _mm_mul_sd(a13_3, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c13_3 = _mm_add_sd(c13_3, _mm_mul_sd(a13_3, b13));
#endif
_mm_store_sd(&C[(i*84)+77], c13_3);
#else
C[(i*84)+13] += values[72] * B[(i*84)+13];
C[(i*84)+28] += values[73] * B[(i*84)+13];
C[(i*84)+49] += values[74] * B[(i*84)+13];
C[(i*84)+77] += values[75] * B[(i*84)+13];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b14 = _mm256_broadcast_sd(&B[(i*84)+14]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b14 = _mm_loaddup_pd(&B[(i*84)+14]);
#endif
__m128d c14_0 = _mm_load_sd(&C[(i*84)+4]);
__m128d a14_0 = _mm_load_sd(&values[76]);
#if defined(__SSE3__) && defined(__AVX__)
c14_0 = _mm_add_sd(c14_0, _mm_mul_sd(a14_0, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_0 = _mm_add_sd(c14_0, _mm_mul_sd(a14_0, b14));
#endif
_mm_store_sd(&C[(i*84)+4], c14_0);
__m128d c14_1 = _mm_load_sd(&C[(i*84)+14]);
__m128d a14_1 = _mm_load_sd(&values[77]);
#if defined(__SSE3__) && defined(__AVX__)
c14_1 = _mm_add_sd(c14_1, _mm_mul_sd(a14_1, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_1 = _mm_add_sd(c14_1, _mm_mul_sd(a14_1, b14));
#endif
_mm_store_sd(&C[(i*84)+14], c14_1);
__m128d c14_2 = _mm_load_sd(&C[(i*84)+29]);
__m128d a14_2 = _mm_load_sd(&values[78]);
#if defined(__SSE3__) && defined(__AVX__)
c14_2 = _mm_add_sd(c14_2, _mm_mul_sd(a14_2, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_2 = _mm_add_sd(c14_2, _mm_mul_sd(a14_2, b14));
#endif
_mm_store_sd(&C[(i*84)+29], c14_2);
__m128d c14_3 = _mm_load_sd(&C[(i*84)+50]);
__m128d a14_3 = _mm_load_sd(&values[79]);
#if defined(__SSE3__) && defined(__AVX__)
c14_3 = _mm_add_sd(c14_3, _mm_mul_sd(a14_3, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_3 = _mm_add_sd(c14_3, _mm_mul_sd(a14_3, b14));
#endif
_mm_store_sd(&C[(i*84)+50], c14_3);
__m128d c14_4 = _mm_load_sd(&C[(i*84)+78]);
__m128d a14_4 = _mm_load_sd(&values[80]);
#if defined(__SSE3__) && defined(__AVX__)
c14_4 = _mm_add_sd(c14_4, _mm_mul_sd(a14_4, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_4 = _mm_add_sd(c14_4, _mm_mul_sd(a14_4, b14));
#endif
_mm_store_sd(&C[(i*84)+78], c14_4);
#else
C[(i*84)+4] += values[76] * B[(i*84)+14];
C[(i*84)+14] += values[77] * B[(i*84)+14];
C[(i*84)+29] += values[78] * B[(i*84)+14];
C[(i*84)+50] += values[79] * B[(i*84)+14];
C[(i*84)+78] += values[80] * B[(i*84)+14];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b15 = _mm256_broadcast_sd(&B[(i*84)+15]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b15 = _mm_loaddup_pd(&B[(i*84)+15]);
#endif
__m128d c15_0 = _mm_load_sd(&C[(i*84)+5]);
__m128d a15_0 = _mm_load_sd(&values[81]);
#if defined(__SSE3__) && defined(__AVX__)
c15_0 = _mm_add_sd(c15_0, _mm_mul_sd(a15_0, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_0 = _mm_add_sd(c15_0, _mm_mul_sd(a15_0, b15));
#endif
_mm_store_sd(&C[(i*84)+5], c15_0);
__m128d c15_1 = _mm_load_sd(&C[(i*84)+15]);
__m128d a15_1 = _mm_load_sd(&values[82]);
#if defined(__SSE3__) && defined(__AVX__)
c15_1 = _mm_add_sd(c15_1, _mm_mul_sd(a15_1, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_1 = _mm_add_sd(c15_1, _mm_mul_sd(a15_1, b15));
#endif
_mm_store_sd(&C[(i*84)+15], c15_1);
__m128d c15_2 = _mm_load_sd(&C[(i*84)+30]);
__m128d a15_2 = _mm_load_sd(&values[83]);
#if defined(__SSE3__) && defined(__AVX__)
c15_2 = _mm_add_sd(c15_2, _mm_mul_sd(a15_2, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_2 = _mm_add_sd(c15_2, _mm_mul_sd(a15_2, b15));
#endif
_mm_store_sd(&C[(i*84)+30], c15_2);
__m128d c15_3 = _mm_load_sd(&C[(i*84)+51]);
__m128d a15_3 = _mm_load_sd(&values[84]);
#if defined(__SSE3__) && defined(__AVX__)
c15_3 = _mm_add_sd(c15_3, _mm_mul_sd(a15_3, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_3 = _mm_add_sd(c15_3, _mm_mul_sd(a15_3, b15));
#endif
_mm_store_sd(&C[(i*84)+51], c15_3);
__m128d c15_4 = _mm_load_sd(&C[(i*84)+79]);
__m128d a15_4 = _mm_load_sd(&values[85]);
#if defined(__SSE3__) && defined(__AVX__)
c15_4 = _mm_add_sd(c15_4, _mm_mul_sd(a15_4, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_4 = _mm_add_sd(c15_4, _mm_mul_sd(a15_4, b15));
#endif
_mm_store_sd(&C[(i*84)+79], c15_4);
#else
C[(i*84)+5] += values[81] * B[(i*84)+15];
C[(i*84)+15] += values[82] * B[(i*84)+15];
C[(i*84)+30] += values[83] * B[(i*84)+15];
C[(i*84)+51] += values[84] * B[(i*84)+15];
C[(i*84)+79] += values[85] * B[(i*84)+15];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b16 = _mm256_broadcast_sd(&B[(i*84)+16]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b16 = _mm_loaddup_pd(&B[(i*84)+16]);
#endif
__m128d c16_0 = _mm_load_sd(&C[(i*84)+6]);
__m128d a16_0 = _mm_load_sd(&values[86]);
#if defined(__SSE3__) && defined(__AVX__)
c16_0 = _mm_add_sd(c16_0, _mm_mul_sd(a16_0, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_0 = _mm_add_sd(c16_0, _mm_mul_sd(a16_0, b16));
#endif
_mm_store_sd(&C[(i*84)+6], c16_0);
__m128d c16_1 = _mm_load_sd(&C[(i*84)+16]);
__m128d a16_1 = _mm_load_sd(&values[87]);
#if defined(__SSE3__) && defined(__AVX__)
c16_1 = _mm_add_sd(c16_1, _mm_mul_sd(a16_1, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_1 = _mm_add_sd(c16_1, _mm_mul_sd(a16_1, b16));
#endif
_mm_store_sd(&C[(i*84)+16], c16_1);
__m128d c16_2 = _mm_load_sd(&C[(i*84)+31]);
__m128d a16_2 = _mm_load_sd(&values[88]);
#if defined(__SSE3__) && defined(__AVX__)
c16_2 = _mm_add_sd(c16_2, _mm_mul_sd(a16_2, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_2 = _mm_add_sd(c16_2, _mm_mul_sd(a16_2, b16));
#endif
_mm_store_sd(&C[(i*84)+31], c16_2);
__m128d c16_3 = _mm_load_sd(&C[(i*84)+52]);
__m128d a16_3 = _mm_load_sd(&values[89]);
#if defined(__SSE3__) && defined(__AVX__)
c16_3 = _mm_add_sd(c16_3, _mm_mul_sd(a16_3, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_3 = _mm_add_sd(c16_3, _mm_mul_sd(a16_3, b16));
#endif
_mm_store_sd(&C[(i*84)+52], c16_3);
__m128d c16_4 = _mm_load_sd(&C[(i*84)+80]);
__m128d a16_4 = _mm_load_sd(&values[90]);
#if defined(__SSE3__) && defined(__AVX__)
c16_4 = _mm_add_sd(c16_4, _mm_mul_sd(a16_4, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_4 = _mm_add_sd(c16_4, _mm_mul_sd(a16_4, b16));
#endif
_mm_store_sd(&C[(i*84)+80], c16_4);
#else
C[(i*84)+6] += values[86] * B[(i*84)+16];
C[(i*84)+16] += values[87] * B[(i*84)+16];
C[(i*84)+31] += values[88] * B[(i*84)+16];
C[(i*84)+52] += values[89] * B[(i*84)+16];
C[(i*84)+80] += values[90] * B[(i*84)+16];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b17 = _mm256_broadcast_sd(&B[(i*84)+17]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b17 = _mm_loaddup_pd(&B[(i*84)+17]);
#endif
__m128d c17_0 = _mm_load_sd(&C[(i*84)+1]);
__m128d a17_0 = _mm_load_sd(&values[91]);
#if defined(__SSE3__) && defined(__AVX__)
c17_0 = _mm_add_sd(c17_0, _mm_mul_sd(a17_0, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_0 = _mm_add_sd(c17_0, _mm_mul_sd(a17_0, b17));
#endif
_mm_store_sd(&C[(i*84)+1], c17_0);
__m128d c17_1 = _mm_load_sd(&C[(i*84)+7]);
__m128d a17_1 = _mm_load_sd(&values[92]);
#if defined(__SSE3__) && defined(__AVX__)
c17_1 = _mm_add_sd(c17_1, _mm_mul_sd(a17_1, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_1 = _mm_add_sd(c17_1, _mm_mul_sd(a17_1, b17));
#endif
_mm_store_sd(&C[(i*84)+7], c17_1);
__m128d c17_2 = _mm_load_sd(&C[(i*84)+17]);
__m128d a17_2 = _mm_load_sd(&values[93]);
#if defined(__SSE3__) && defined(__AVX__)
c17_2 = _mm_add_sd(c17_2, _mm_mul_sd(a17_2, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_2 = _mm_add_sd(c17_2, _mm_mul_sd(a17_2, b17));
#endif
_mm_store_sd(&C[(i*84)+17], c17_2);
__m128d c17_3 = _mm_load_sd(&C[(i*84)+32]);
__m128d a17_3 = _mm_load_sd(&values[94]);
#if defined(__SSE3__) && defined(__AVX__)
c17_3 = _mm_add_sd(c17_3, _mm_mul_sd(a17_3, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_3 = _mm_add_sd(c17_3, _mm_mul_sd(a17_3, b17));
#endif
_mm_store_sd(&C[(i*84)+32], c17_3);
__m128d c17_4 = _mm_load_sd(&C[(i*84)+53]);
__m128d a17_4 = _mm_load_sd(&values[95]);
#if defined(__SSE3__) && defined(__AVX__)
c17_4 = _mm_add_sd(c17_4, _mm_mul_sd(a17_4, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_4 = _mm_add_sd(c17_4, _mm_mul_sd(a17_4, b17));
#endif
_mm_store_sd(&C[(i*84)+53], c17_4);
__m128d c17_5 = _mm_load_sd(&C[(i*84)+81]);
__m128d a17_5 = _mm_load_sd(&values[96]);
#if defined(__SSE3__) && defined(__AVX__)
c17_5 = _mm_add_sd(c17_5, _mm_mul_sd(a17_5, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_5 = _mm_add_sd(c17_5, _mm_mul_sd(a17_5, b17));
#endif
_mm_store_sd(&C[(i*84)+81], c17_5);
#else
C[(i*84)+1] += values[91] * B[(i*84)+17];
C[(i*84)+7] += values[92] * B[(i*84)+17];
C[(i*84)+17] += values[93] * B[(i*84)+17];
C[(i*84)+32] += values[94] * B[(i*84)+17];
C[(i*84)+53] += values[95] * B[(i*84)+17];
C[(i*84)+81] += values[96] * B[(i*84)+17];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b18 = _mm256_broadcast_sd(&B[(i*84)+18]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b18 = _mm_loaddup_pd(&B[(i*84)+18]);
#endif
__m128d c18_0 = _mm_load_sd(&C[(i*84)+2]);
__m128d a18_0 = _mm_load_sd(&values[97]);
#if defined(__SSE3__) && defined(__AVX__)
c18_0 = _mm_add_sd(c18_0, _mm_mul_sd(a18_0, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_0 = _mm_add_sd(c18_0, _mm_mul_sd(a18_0, b18));
#endif
_mm_store_sd(&C[(i*84)+2], c18_0);
__m128d c18_1 = _mm_load_sd(&C[(i*84)+8]);
__m128d a18_1 = _mm_load_sd(&values[98]);
#if defined(__SSE3__) && defined(__AVX__)
c18_1 = _mm_add_sd(c18_1, _mm_mul_sd(a18_1, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_1 = _mm_add_sd(c18_1, _mm_mul_sd(a18_1, b18));
#endif
_mm_store_sd(&C[(i*84)+8], c18_1);
__m128d c18_2 = _mm_load_sd(&C[(i*84)+18]);
__m128d a18_2 = _mm_load_sd(&values[99]);
#if defined(__SSE3__) && defined(__AVX__)
c18_2 = _mm_add_sd(c18_2, _mm_mul_sd(a18_2, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_2 = _mm_add_sd(c18_2, _mm_mul_sd(a18_2, b18));
#endif
_mm_store_sd(&C[(i*84)+18], c18_2);
__m128d c18_3 = _mm_load_sd(&C[(i*84)+33]);
__m128d a18_3 = _mm_load_sd(&values[100]);
#if defined(__SSE3__) && defined(__AVX__)
c18_3 = _mm_add_sd(c18_3, _mm_mul_sd(a18_3, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_3 = _mm_add_sd(c18_3, _mm_mul_sd(a18_3, b18));
#endif
_mm_store_sd(&C[(i*84)+33], c18_3);
__m128d c18_4 = _mm_load_sd(&C[(i*84)+54]);
__m128d a18_4 = _mm_load_sd(&values[101]);
#if defined(__SSE3__) && defined(__AVX__)
c18_4 = _mm_add_sd(c18_4, _mm_mul_sd(a18_4, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_4 = _mm_add_sd(c18_4, _mm_mul_sd(a18_4, b18));
#endif
_mm_store_sd(&C[(i*84)+54], c18_4);
__m128d c18_5 = _mm_load_sd(&C[(i*84)+82]);
__m128d a18_5 = _mm_load_sd(&values[102]);
#if defined(__SSE3__) && defined(__AVX__)
c18_5 = _mm_add_sd(c18_5, _mm_mul_sd(a18_5, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_5 = _mm_add_sd(c18_5, _mm_mul_sd(a18_5, b18));
#endif
_mm_store_sd(&C[(i*84)+82], c18_5);
#else
C[(i*84)+2] += values[97] * B[(i*84)+18];
C[(i*84)+8] += values[98] * B[(i*84)+18];
C[(i*84)+18] += values[99] * B[(i*84)+18];
C[(i*84)+33] += values[100] * B[(i*84)+18];
C[(i*84)+54] += values[101] * B[(i*84)+18];
C[(i*84)+82] += values[102] * B[(i*84)+18];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b19 = _mm256_broadcast_sd(&B[(i*84)+19]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b19 = _mm_loaddup_pd(&B[(i*84)+19]);
#endif
__m128d c19_0 = _mm_load_sd(&C[(i*84)+0]);
__m128d a19_0 = _mm_load_sd(&values[103]);
#if defined(__SSE3__) && defined(__AVX__)
c19_0 = _mm_add_sd(c19_0, _mm_mul_sd(a19_0, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_0 = _mm_add_sd(c19_0, _mm_mul_sd(a19_0, b19));
#endif
_mm_store_sd(&C[(i*84)+0], c19_0);
__m128d c19_1 = _mm_load_sd(&C[(i*84)+3]);
__m128d a19_1 = _mm_load_sd(&values[104]);
#if defined(__SSE3__) && defined(__AVX__)
c19_1 = _mm_add_sd(c19_1, _mm_mul_sd(a19_1, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_1 = _mm_add_sd(c19_1, _mm_mul_sd(a19_1, b19));
#endif
_mm_store_sd(&C[(i*84)+3], c19_1);
__m128d c19_2 = _mm_load_sd(&C[(i*84)+9]);
__m128d a19_2 = _mm_load_sd(&values[105]);
#if defined(__SSE3__) && defined(__AVX__)
c19_2 = _mm_add_sd(c19_2, _mm_mul_sd(a19_2, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_2 = _mm_add_sd(c19_2, _mm_mul_sd(a19_2, b19));
#endif
_mm_store_sd(&C[(i*84)+9], c19_2);
__m128d c19_3 = _mm_load_sd(&C[(i*84)+19]);
__m128d a19_3 = _mm_load_sd(&values[106]);
#if defined(__SSE3__) && defined(__AVX__)
c19_3 = _mm_add_sd(c19_3, _mm_mul_sd(a19_3, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_3 = _mm_add_sd(c19_3, _mm_mul_sd(a19_3, b19));
#endif
_mm_store_sd(&C[(i*84)+19], c19_3);
__m128d c19_4 = _mm_load_sd(&C[(i*84)+34]);
__m128d a19_4 = _mm_load_sd(&values[107]);
#if defined(__SSE3__) && defined(__AVX__)
c19_4 = _mm_add_sd(c19_4, _mm_mul_sd(a19_4, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_4 = _mm_add_sd(c19_4, _mm_mul_sd(a19_4, b19));
#endif
_mm_store_sd(&C[(i*84)+34], c19_4);
__m128d c19_5 = _mm_load_sd(&C[(i*84)+55]);
__m128d a19_5 = _mm_load_sd(&values[108]);
#if defined(__SSE3__) && defined(__AVX__)
c19_5 = _mm_add_sd(c19_5, _mm_mul_sd(a19_5, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_5 = _mm_add_sd(c19_5, _mm_mul_sd(a19_5, b19));
#endif
_mm_store_sd(&C[(i*84)+55], c19_5);
__m128d c19_6 = _mm_load_sd(&C[(i*84)+83]);
__m128d a19_6 = _mm_load_sd(&values[109]);
#if defined(__SSE3__) && defined(__AVX__)
c19_6 = _mm_add_sd(c19_6, _mm_mul_sd(a19_6, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_6 = _mm_add_sd(c19_6, _mm_mul_sd(a19_6, b19));
#endif
_mm_store_sd(&C[(i*84)+83], c19_6);
#else
C[(i*84)+0] += values[103] * B[(i*84)+19];
C[(i*84)+3] += values[104] * B[(i*84)+19];
C[(i*84)+9] += values[105] * B[(i*84)+19];
C[(i*84)+19] += values[106] * B[(i*84)+19];
C[(i*84)+34] += values[107] * B[(i*84)+19];
C[(i*84)+55] += values[108] * B[(i*84)+19];
C[(i*84)+83] += values[109] * B[(i*84)+19];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b20 = _mm256_broadcast_sd(&B[(i*84)+20]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b20 = _mm_loaddup_pd(&B[(i*84)+20]);
#endif
__m128d c20_0 = _mm_load_sd(&C[(i*84)+20]);
__m128d a20_0 = _mm_load_sd(&values[110]);
#if defined(__SSE3__) && defined(__AVX__)
c20_0 = _mm_add_sd(c20_0, _mm_mul_sd(a20_0, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c20_0 = _mm_add_sd(c20_0, _mm_mul_sd(a20_0, b20));
#endif
_mm_store_sd(&C[(i*84)+20], c20_0);
__m128d c20_1 = _mm_load_sd(&C[(i*84)+41]);
__m128d a20_1 = _mm_load_sd(&values[111]);
#if defined(__SSE3__) && defined(__AVX__)
c20_1 = _mm_add_sd(c20_1, _mm_mul_sd(a20_1, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c20_1 = _mm_add_sd(c20_1, _mm_mul_sd(a20_1, b20));
#endif
_mm_store_sd(&C[(i*84)+41], c20_1);
__m128d c20_2 = _mm_load_sd(&C[(i*84)+69]);
__m128d a20_2 = _mm_load_sd(&values[112]);
#if defined(__SSE3__) && defined(__AVX__)
c20_2 = _mm_add_sd(c20_2, _mm_mul_sd(a20_2, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c20_2 = _mm_add_sd(c20_2, _mm_mul_sd(a20_2, b20));
#endif
_mm_store_sd(&C[(i*84)+69], c20_2);
#else
C[(i*84)+20] += values[110] * B[(i*84)+20];
C[(i*84)+41] += values[111] * B[(i*84)+20];
C[(i*84)+69] += values[112] * B[(i*84)+20];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b21 = _mm256_broadcast_sd(&B[(i*84)+21]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b21 = _mm_loaddup_pd(&B[(i*84)+21]);
#endif
__m128d c21_0 = _mm_load_sd(&C[(i*84)+21]);
__m128d a21_0 = _mm_load_sd(&values[113]);
#if defined(__SSE3__) && defined(__AVX__)
c21_0 = _mm_add_sd(c21_0, _mm_mul_sd(a21_0, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c21_0 = _mm_add_sd(c21_0, _mm_mul_sd(a21_0, b21));
#endif
_mm_store_sd(&C[(i*84)+21], c21_0);
__m128d c21_1 = _mm_load_sd(&C[(i*84)+42]);
__m128d a21_1 = _mm_load_sd(&values[114]);
#if defined(__SSE3__) && defined(__AVX__)
c21_1 = _mm_add_sd(c21_1, _mm_mul_sd(a21_1, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c21_1 = _mm_add_sd(c21_1, _mm_mul_sd(a21_1, b21));
#endif
_mm_store_sd(&C[(i*84)+42], c21_1);
__m128d c21_2 = _mm_load_sd(&C[(i*84)+70]);
__m128d a21_2 = _mm_load_sd(&values[115]);
#if defined(__SSE3__) && defined(__AVX__)
c21_2 = _mm_add_sd(c21_2, _mm_mul_sd(a21_2, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c21_2 = _mm_add_sd(c21_2, _mm_mul_sd(a21_2, b21));
#endif
_mm_store_sd(&C[(i*84)+70], c21_2);
#else
C[(i*84)+21] += values[113] * B[(i*84)+21];
C[(i*84)+42] += values[114] * B[(i*84)+21];
C[(i*84)+70] += values[115] * B[(i*84)+21];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b22 = _mm256_broadcast_sd(&B[(i*84)+22]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b22 = _mm_loaddup_pd(&B[(i*84)+22]);
#endif
__m128d c22_0 = _mm_load_sd(&C[(i*84)+22]);
__m128d a22_0 = _mm_load_sd(&values[116]);
#if defined(__SSE3__) && defined(__AVX__)
c22_0 = _mm_add_sd(c22_0, _mm_mul_sd(a22_0, _mm256_castpd256_pd128(b22)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c22_0 = _mm_add_sd(c22_0, _mm_mul_sd(a22_0, b22));
#endif
_mm_store_sd(&C[(i*84)+22], c22_0);
__m128d c22_1 = _mm_load_sd(&C[(i*84)+43]);
__m128d a22_1 = _mm_load_sd(&values[117]);
#if defined(__SSE3__) && defined(__AVX__)
c22_1 = _mm_add_sd(c22_1, _mm_mul_sd(a22_1, _mm256_castpd256_pd128(b22)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c22_1 = _mm_add_sd(c22_1, _mm_mul_sd(a22_1, b22));
#endif
_mm_store_sd(&C[(i*84)+43], c22_1);
__m128d c22_2 = _mm_load_sd(&C[(i*84)+71]);
__m128d a22_2 = _mm_load_sd(&values[118]);
#if defined(__SSE3__) && defined(__AVX__)
c22_2 = _mm_add_sd(c22_2, _mm_mul_sd(a22_2, _mm256_castpd256_pd128(b22)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c22_2 = _mm_add_sd(c22_2, _mm_mul_sd(a22_2, b22));
#endif
_mm_store_sd(&C[(i*84)+71], c22_2);
#else
C[(i*84)+22] += values[116] * B[(i*84)+22];
C[(i*84)+43] += values[117] * B[(i*84)+22];
C[(i*84)+71] += values[118] * B[(i*84)+22];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b23 = _mm256_broadcast_sd(&B[(i*84)+23]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b23 = _mm_loaddup_pd(&B[(i*84)+23]);
#endif
__m128d c23_0 = _mm_load_sd(&C[(i*84)+23]);
__m128d a23_0 = _mm_load_sd(&values[119]);
#if defined(__SSE3__) && defined(__AVX__)
c23_0 = _mm_add_sd(c23_0, _mm_mul_sd(a23_0, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c23_0 = _mm_add_sd(c23_0, _mm_mul_sd(a23_0, b23));
#endif
_mm_store_sd(&C[(i*84)+23], c23_0);
__m128d c23_1 = _mm_load_sd(&C[(i*84)+44]);
__m128d a23_1 = _mm_load_sd(&values[120]);
#if defined(__SSE3__) && defined(__AVX__)
c23_1 = _mm_add_sd(c23_1, _mm_mul_sd(a23_1, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c23_1 = _mm_add_sd(c23_1, _mm_mul_sd(a23_1, b23));
#endif
_mm_store_sd(&C[(i*84)+44], c23_1);
__m128d c23_2 = _mm_load_sd(&C[(i*84)+72]);
__m128d a23_2 = _mm_load_sd(&values[121]);
#if defined(__SSE3__) && defined(__AVX__)
c23_2 = _mm_add_sd(c23_2, _mm_mul_sd(a23_2, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c23_2 = _mm_add_sd(c23_2, _mm_mul_sd(a23_2, b23));
#endif
_mm_store_sd(&C[(i*84)+72], c23_2);
#else
C[(i*84)+23] += values[119] * B[(i*84)+23];
C[(i*84)+44] += values[120] * B[(i*84)+23];
C[(i*84)+72] += values[121] * B[(i*84)+23];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b24 = _mm256_broadcast_sd(&B[(i*84)+24]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b24 = _mm_loaddup_pd(&B[(i*84)+24]);
#endif
__m128d c24_0 = _mm_load_sd(&C[(i*84)+24]);
__m128d a24_0 = _mm_load_sd(&values[122]);
#if defined(__SSE3__) && defined(__AVX__)
c24_0 = _mm_add_sd(c24_0, _mm_mul_sd(a24_0, _mm256_castpd256_pd128(b24)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c24_0 = _mm_add_sd(c24_0, _mm_mul_sd(a24_0, b24));
#endif
_mm_store_sd(&C[(i*84)+24], c24_0);
__m128d c24_1 = _mm_load_sd(&C[(i*84)+45]);
__m128d a24_1 = _mm_load_sd(&values[123]);
#if defined(__SSE3__) && defined(__AVX__)
c24_1 = _mm_add_sd(c24_1, _mm_mul_sd(a24_1, _mm256_castpd256_pd128(b24)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c24_1 = _mm_add_sd(c24_1, _mm_mul_sd(a24_1, b24));
#endif
_mm_store_sd(&C[(i*84)+45], c24_1);
__m128d c24_2 = _mm_load_sd(&C[(i*84)+73]);
__m128d a24_2 = _mm_load_sd(&values[124]);
#if defined(__SSE3__) && defined(__AVX__)
c24_2 = _mm_add_sd(c24_2, _mm_mul_sd(a24_2, _mm256_castpd256_pd128(b24)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c24_2 = _mm_add_sd(c24_2, _mm_mul_sd(a24_2, b24));
#endif
_mm_store_sd(&C[(i*84)+73], c24_2);
#else
C[(i*84)+24] += values[122] * B[(i*84)+24];
C[(i*84)+45] += values[123] * B[(i*84)+24];
C[(i*84)+73] += values[124] * B[(i*84)+24];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b25 = _mm256_broadcast_sd(&B[(i*84)+25]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b25 = _mm_loaddup_pd(&B[(i*84)+25]);
#endif
__m128d c25_0 = _mm_load_sd(&C[(i*84)+10]);
__m128d a25_0 = _mm_load_sd(&values[125]);
#if defined(__SSE3__) && defined(__AVX__)
c25_0 = _mm_add_sd(c25_0, _mm_mul_sd(a25_0, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c25_0 = _mm_add_sd(c25_0, _mm_mul_sd(a25_0, b25));
#endif
_mm_store_sd(&C[(i*84)+10], c25_0);
__m128d c25_1 = _mm_load_sd(&C[(i*84)+25]);
__m128d a25_1 = _mm_load_sd(&values[126]);
#if defined(__SSE3__) && defined(__AVX__)
c25_1 = _mm_add_sd(c25_1, _mm_mul_sd(a25_1, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c25_1 = _mm_add_sd(c25_1, _mm_mul_sd(a25_1, b25));
#endif
_mm_store_sd(&C[(i*84)+25], c25_1);
__m128d c25_2 = _mm_load_sd(&C[(i*84)+46]);
__m128d a25_2 = _mm_load_sd(&values[127]);
#if defined(__SSE3__) && defined(__AVX__)
c25_2 = _mm_add_sd(c25_2, _mm_mul_sd(a25_2, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c25_2 = _mm_add_sd(c25_2, _mm_mul_sd(a25_2, b25));
#endif
_mm_store_sd(&C[(i*84)+46], c25_2);
__m128d c25_3 = _mm_load_sd(&C[(i*84)+74]);
__m128d a25_3 = _mm_load_sd(&values[128]);
#if defined(__SSE3__) && defined(__AVX__)
c25_3 = _mm_add_sd(c25_3, _mm_mul_sd(a25_3, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c25_3 = _mm_add_sd(c25_3, _mm_mul_sd(a25_3, b25));
#endif
_mm_store_sd(&C[(i*84)+74], c25_3);
#else
C[(i*84)+10] += values[125] * B[(i*84)+25];
C[(i*84)+25] += values[126] * B[(i*84)+25];
C[(i*84)+46] += values[127] * B[(i*84)+25];
C[(i*84)+74] += values[128] * B[(i*84)+25];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b26 = _mm256_broadcast_sd(&B[(i*84)+26]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b26 = _mm_loaddup_pd(&B[(i*84)+26]);
#endif
__m128d c26_0 = _mm_load_sd(&C[(i*84)+11]);
__m128d a26_0 = _mm_load_sd(&values[129]);
#if defined(__SSE3__) && defined(__AVX__)
c26_0 = _mm_add_sd(c26_0, _mm_mul_sd(a26_0, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c26_0 = _mm_add_sd(c26_0, _mm_mul_sd(a26_0, b26));
#endif
_mm_store_sd(&C[(i*84)+11], c26_0);
__m128d c26_1 = _mm_load_sd(&C[(i*84)+26]);
__m128d a26_1 = _mm_load_sd(&values[130]);
#if defined(__SSE3__) && defined(__AVX__)
c26_1 = _mm_add_sd(c26_1, _mm_mul_sd(a26_1, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c26_1 = _mm_add_sd(c26_1, _mm_mul_sd(a26_1, b26));
#endif
_mm_store_sd(&C[(i*84)+26], c26_1);
__m128d c26_2 = _mm_load_sd(&C[(i*84)+47]);
__m128d a26_2 = _mm_load_sd(&values[131]);
#if defined(__SSE3__) && defined(__AVX__)
c26_2 = _mm_add_sd(c26_2, _mm_mul_sd(a26_2, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c26_2 = _mm_add_sd(c26_2, _mm_mul_sd(a26_2, b26));
#endif
_mm_store_sd(&C[(i*84)+47], c26_2);
__m128d c26_3 = _mm_load_sd(&C[(i*84)+75]);
__m128d a26_3 = _mm_load_sd(&values[132]);
#if defined(__SSE3__) && defined(__AVX__)
c26_3 = _mm_add_sd(c26_3, _mm_mul_sd(a26_3, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c26_3 = _mm_add_sd(c26_3, _mm_mul_sd(a26_3, b26));
#endif
_mm_store_sd(&C[(i*84)+75], c26_3);
#else
C[(i*84)+11] += values[129] * B[(i*84)+26];
C[(i*84)+26] += values[130] * B[(i*84)+26];
C[(i*84)+47] += values[131] * B[(i*84)+26];
C[(i*84)+75] += values[132] * B[(i*84)+26];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b27 = _mm256_broadcast_sd(&B[(i*84)+27]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b27 = _mm_loaddup_pd(&B[(i*84)+27]);
#endif
__m128d c27_0 = _mm_load_sd(&C[(i*84)+12]);
__m128d a27_0 = _mm_load_sd(&values[133]);
#if defined(__SSE3__) && defined(__AVX__)
c27_0 = _mm_add_sd(c27_0, _mm_mul_sd(a27_0, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c27_0 = _mm_add_sd(c27_0, _mm_mul_sd(a27_0, b27));
#endif
_mm_store_sd(&C[(i*84)+12], c27_0);
__m128d c27_1 = _mm_load_sd(&C[(i*84)+27]);
__m128d a27_1 = _mm_load_sd(&values[134]);
#if defined(__SSE3__) && defined(__AVX__)
c27_1 = _mm_add_sd(c27_1, _mm_mul_sd(a27_1, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c27_1 = _mm_add_sd(c27_1, _mm_mul_sd(a27_1, b27));
#endif
_mm_store_sd(&C[(i*84)+27], c27_1);
__m128d c27_2 = _mm_load_sd(&C[(i*84)+48]);
__m128d a27_2 = _mm_load_sd(&values[135]);
#if defined(__SSE3__) && defined(__AVX__)
c27_2 = _mm_add_sd(c27_2, _mm_mul_sd(a27_2, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c27_2 = _mm_add_sd(c27_2, _mm_mul_sd(a27_2, b27));
#endif
_mm_store_sd(&C[(i*84)+48], c27_2);
__m128d c27_3 = _mm_load_sd(&C[(i*84)+76]);
__m128d a27_3 = _mm_load_sd(&values[136]);
#if defined(__SSE3__) && defined(__AVX__)
c27_3 = _mm_add_sd(c27_3, _mm_mul_sd(a27_3, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c27_3 = _mm_add_sd(c27_3, _mm_mul_sd(a27_3, b27));
#endif
_mm_store_sd(&C[(i*84)+76], c27_3);
#else
C[(i*84)+12] += values[133] * B[(i*84)+27];
C[(i*84)+27] += values[134] * B[(i*84)+27];
C[(i*84)+48] += values[135] * B[(i*84)+27];
C[(i*84)+76] += values[136] * B[(i*84)+27];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b28 = _mm256_broadcast_sd(&B[(i*84)+28]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b28 = _mm_loaddup_pd(&B[(i*84)+28]);
#endif
__m128d c28_0 = _mm_load_sd(&C[(i*84)+13]);
__m128d a28_0 = _mm_load_sd(&values[137]);
#if defined(__SSE3__) && defined(__AVX__)
c28_0 = _mm_add_sd(c28_0, _mm_mul_sd(a28_0, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c28_0 = _mm_add_sd(c28_0, _mm_mul_sd(a28_0, b28));
#endif
_mm_store_sd(&C[(i*84)+13], c28_0);
__m128d c28_1 = _mm_load_sd(&C[(i*84)+28]);
__m128d a28_1 = _mm_load_sd(&values[138]);
#if defined(__SSE3__) && defined(__AVX__)
c28_1 = _mm_add_sd(c28_1, _mm_mul_sd(a28_1, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c28_1 = _mm_add_sd(c28_1, _mm_mul_sd(a28_1, b28));
#endif
_mm_store_sd(&C[(i*84)+28], c28_1);
__m128d c28_2 = _mm_load_sd(&C[(i*84)+49]);
__m128d a28_2 = _mm_load_sd(&values[139]);
#if defined(__SSE3__) && defined(__AVX__)
c28_2 = _mm_add_sd(c28_2, _mm_mul_sd(a28_2, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c28_2 = _mm_add_sd(c28_2, _mm_mul_sd(a28_2, b28));
#endif
_mm_store_sd(&C[(i*84)+49], c28_2);
__m128d c28_3 = _mm_load_sd(&C[(i*84)+77]);
__m128d a28_3 = _mm_load_sd(&values[140]);
#if defined(__SSE3__) && defined(__AVX__)
c28_3 = _mm_add_sd(c28_3, _mm_mul_sd(a28_3, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c28_3 = _mm_add_sd(c28_3, _mm_mul_sd(a28_3, b28));
#endif
_mm_store_sd(&C[(i*84)+77], c28_3);
#else
C[(i*84)+13] += values[137] * B[(i*84)+28];
C[(i*84)+28] += values[138] * B[(i*84)+28];
C[(i*84)+49] += values[139] * B[(i*84)+28];
C[(i*84)+77] += values[140] * B[(i*84)+28];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b29 = _mm256_broadcast_sd(&B[(i*84)+29]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b29 = _mm_loaddup_pd(&B[(i*84)+29]);
#endif
__m128d c29_0 = _mm_load_sd(&C[(i*84)+4]);
__m128d a29_0 = _mm_load_sd(&values[141]);
#if defined(__SSE3__) && defined(__AVX__)
c29_0 = _mm_add_sd(c29_0, _mm_mul_sd(a29_0, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_0 = _mm_add_sd(c29_0, _mm_mul_sd(a29_0, b29));
#endif
_mm_store_sd(&C[(i*84)+4], c29_0);
__m128d c29_1 = _mm_load_sd(&C[(i*84)+14]);
__m128d a29_1 = _mm_load_sd(&values[142]);
#if defined(__SSE3__) && defined(__AVX__)
c29_1 = _mm_add_sd(c29_1, _mm_mul_sd(a29_1, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_1 = _mm_add_sd(c29_1, _mm_mul_sd(a29_1, b29));
#endif
_mm_store_sd(&C[(i*84)+14], c29_1);
__m128d c29_2 = _mm_load_sd(&C[(i*84)+29]);
__m128d a29_2 = _mm_load_sd(&values[143]);
#if defined(__SSE3__) && defined(__AVX__)
c29_2 = _mm_add_sd(c29_2, _mm_mul_sd(a29_2, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_2 = _mm_add_sd(c29_2, _mm_mul_sd(a29_2, b29));
#endif
_mm_store_sd(&C[(i*84)+29], c29_2);
__m128d c29_3 = _mm_load_sd(&C[(i*84)+50]);
__m128d a29_3 = _mm_load_sd(&values[144]);
#if defined(__SSE3__) && defined(__AVX__)
c29_3 = _mm_add_sd(c29_3, _mm_mul_sd(a29_3, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_3 = _mm_add_sd(c29_3, _mm_mul_sd(a29_3, b29));
#endif
_mm_store_sd(&C[(i*84)+50], c29_3);
__m128d c29_4 = _mm_load_sd(&C[(i*84)+78]);
__m128d a29_4 = _mm_load_sd(&values[145]);
#if defined(__SSE3__) && defined(__AVX__)
c29_4 = _mm_add_sd(c29_4, _mm_mul_sd(a29_4, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_4 = _mm_add_sd(c29_4, _mm_mul_sd(a29_4, b29));
#endif
_mm_store_sd(&C[(i*84)+78], c29_4);
#else
C[(i*84)+4] += values[141] * B[(i*84)+29];
C[(i*84)+14] += values[142] * B[(i*84)+29];
C[(i*84)+29] += values[143] * B[(i*84)+29];
C[(i*84)+50] += values[144] * B[(i*84)+29];
C[(i*84)+78] += values[145] * B[(i*84)+29];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b30 = _mm256_broadcast_sd(&B[(i*84)+30]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b30 = _mm_loaddup_pd(&B[(i*84)+30]);
#endif
__m128d c30_0 = _mm_load_sd(&C[(i*84)+5]);
__m128d a30_0 = _mm_load_sd(&values[146]);
#if defined(__SSE3__) && defined(__AVX__)
c30_0 = _mm_add_sd(c30_0, _mm_mul_sd(a30_0, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_0 = _mm_add_sd(c30_0, _mm_mul_sd(a30_0, b30));
#endif
_mm_store_sd(&C[(i*84)+5], c30_0);
__m128d c30_1 = _mm_load_sd(&C[(i*84)+15]);
__m128d a30_1 = _mm_load_sd(&values[147]);
#if defined(__SSE3__) && defined(__AVX__)
c30_1 = _mm_add_sd(c30_1, _mm_mul_sd(a30_1, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_1 = _mm_add_sd(c30_1, _mm_mul_sd(a30_1, b30));
#endif
_mm_store_sd(&C[(i*84)+15], c30_1);
__m128d c30_2 = _mm_load_sd(&C[(i*84)+30]);
__m128d a30_2 = _mm_load_sd(&values[148]);
#if defined(__SSE3__) && defined(__AVX__)
c30_2 = _mm_add_sd(c30_2, _mm_mul_sd(a30_2, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_2 = _mm_add_sd(c30_2, _mm_mul_sd(a30_2, b30));
#endif
_mm_store_sd(&C[(i*84)+30], c30_2);
__m128d c30_3 = _mm_load_sd(&C[(i*84)+51]);
__m128d a30_3 = _mm_load_sd(&values[149]);
#if defined(__SSE3__) && defined(__AVX__)
c30_3 = _mm_add_sd(c30_3, _mm_mul_sd(a30_3, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_3 = _mm_add_sd(c30_3, _mm_mul_sd(a30_3, b30));
#endif
_mm_store_sd(&C[(i*84)+51], c30_3);
__m128d c30_4 = _mm_load_sd(&C[(i*84)+79]);
__m128d a30_4 = _mm_load_sd(&values[150]);
#if defined(__SSE3__) && defined(__AVX__)
c30_4 = _mm_add_sd(c30_4, _mm_mul_sd(a30_4, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_4 = _mm_add_sd(c30_4, _mm_mul_sd(a30_4, b30));
#endif
_mm_store_sd(&C[(i*84)+79], c30_4);
#else
C[(i*84)+5] += values[146] * B[(i*84)+30];
C[(i*84)+15] += values[147] * B[(i*84)+30];
C[(i*84)+30] += values[148] * B[(i*84)+30];
C[(i*84)+51] += values[149] * B[(i*84)+30];
C[(i*84)+79] += values[150] * B[(i*84)+30];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b31 = _mm256_broadcast_sd(&B[(i*84)+31]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b31 = _mm_loaddup_pd(&B[(i*84)+31]);
#endif
__m128d c31_0 = _mm_load_sd(&C[(i*84)+6]);
__m128d a31_0 = _mm_load_sd(&values[151]);
#if defined(__SSE3__) && defined(__AVX__)
c31_0 = _mm_add_sd(c31_0, _mm_mul_sd(a31_0, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_0 = _mm_add_sd(c31_0, _mm_mul_sd(a31_0, b31));
#endif
_mm_store_sd(&C[(i*84)+6], c31_0);
__m128d c31_1 = _mm_load_sd(&C[(i*84)+16]);
__m128d a31_1 = _mm_load_sd(&values[152]);
#if defined(__SSE3__) && defined(__AVX__)
c31_1 = _mm_add_sd(c31_1, _mm_mul_sd(a31_1, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_1 = _mm_add_sd(c31_1, _mm_mul_sd(a31_1, b31));
#endif
_mm_store_sd(&C[(i*84)+16], c31_1);
__m128d c31_2 = _mm_load_sd(&C[(i*84)+31]);
__m128d a31_2 = _mm_load_sd(&values[153]);
#if defined(__SSE3__) && defined(__AVX__)
c31_2 = _mm_add_sd(c31_2, _mm_mul_sd(a31_2, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_2 = _mm_add_sd(c31_2, _mm_mul_sd(a31_2, b31));
#endif
_mm_store_sd(&C[(i*84)+31], c31_2);
__m128d c31_3 = _mm_load_sd(&C[(i*84)+52]);
__m128d a31_3 = _mm_load_sd(&values[154]);
#if defined(__SSE3__) && defined(__AVX__)
c31_3 = _mm_add_sd(c31_3, _mm_mul_sd(a31_3, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_3 = _mm_add_sd(c31_3, _mm_mul_sd(a31_3, b31));
#endif
_mm_store_sd(&C[(i*84)+52], c31_3);
__m128d c31_4 = _mm_load_sd(&C[(i*84)+80]);
__m128d a31_4 = _mm_load_sd(&values[155]);
#if defined(__SSE3__) && defined(__AVX__)
c31_4 = _mm_add_sd(c31_4, _mm_mul_sd(a31_4, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_4 = _mm_add_sd(c31_4, _mm_mul_sd(a31_4, b31));
#endif
_mm_store_sd(&C[(i*84)+80], c31_4);
#else
C[(i*84)+6] += values[151] * B[(i*84)+31];
C[(i*84)+16] += values[152] * B[(i*84)+31];
C[(i*84)+31] += values[153] * B[(i*84)+31];
C[(i*84)+52] += values[154] * B[(i*84)+31];
C[(i*84)+80] += values[155] * B[(i*84)+31];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b32 = _mm256_broadcast_sd(&B[(i*84)+32]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b32 = _mm_loaddup_pd(&B[(i*84)+32]);
#endif
__m128d c32_0 = _mm_load_sd(&C[(i*84)+1]);
__m128d a32_0 = _mm_load_sd(&values[156]);
#if defined(__SSE3__) && defined(__AVX__)
c32_0 = _mm_add_sd(c32_0, _mm_mul_sd(a32_0, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_0 = _mm_add_sd(c32_0, _mm_mul_sd(a32_0, b32));
#endif
_mm_store_sd(&C[(i*84)+1], c32_0);
__m128d c32_1 = _mm_load_sd(&C[(i*84)+7]);
__m128d a32_1 = _mm_load_sd(&values[157]);
#if defined(__SSE3__) && defined(__AVX__)
c32_1 = _mm_add_sd(c32_1, _mm_mul_sd(a32_1, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_1 = _mm_add_sd(c32_1, _mm_mul_sd(a32_1, b32));
#endif
_mm_store_sd(&C[(i*84)+7], c32_1);
__m128d c32_2 = _mm_load_sd(&C[(i*84)+17]);
__m128d a32_2 = _mm_load_sd(&values[158]);
#if defined(__SSE3__) && defined(__AVX__)
c32_2 = _mm_add_sd(c32_2, _mm_mul_sd(a32_2, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_2 = _mm_add_sd(c32_2, _mm_mul_sd(a32_2, b32));
#endif
_mm_store_sd(&C[(i*84)+17], c32_2);
__m128d c32_3 = _mm_load_sd(&C[(i*84)+32]);
__m128d a32_3 = _mm_load_sd(&values[159]);
#if defined(__SSE3__) && defined(__AVX__)
c32_3 = _mm_add_sd(c32_3, _mm_mul_sd(a32_3, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_3 = _mm_add_sd(c32_3, _mm_mul_sd(a32_3, b32));
#endif
_mm_store_sd(&C[(i*84)+32], c32_3);
__m128d c32_4 = _mm_load_sd(&C[(i*84)+53]);
__m128d a32_4 = _mm_load_sd(&values[160]);
#if defined(__SSE3__) && defined(__AVX__)
c32_4 = _mm_add_sd(c32_4, _mm_mul_sd(a32_4, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_4 = _mm_add_sd(c32_4, _mm_mul_sd(a32_4, b32));
#endif
_mm_store_sd(&C[(i*84)+53], c32_4);
__m128d c32_5 = _mm_load_sd(&C[(i*84)+81]);
__m128d a32_5 = _mm_load_sd(&values[161]);
#if defined(__SSE3__) && defined(__AVX__)
c32_5 = _mm_add_sd(c32_5, _mm_mul_sd(a32_5, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_5 = _mm_add_sd(c32_5, _mm_mul_sd(a32_5, b32));
#endif
_mm_store_sd(&C[(i*84)+81], c32_5);
#else
C[(i*84)+1] += values[156] * B[(i*84)+32];
C[(i*84)+7] += values[157] * B[(i*84)+32];
C[(i*84)+17] += values[158] * B[(i*84)+32];
C[(i*84)+32] += values[159] * B[(i*84)+32];
C[(i*84)+53] += values[160] * B[(i*84)+32];
C[(i*84)+81] += values[161] * B[(i*84)+32];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b33 = _mm256_broadcast_sd(&B[(i*84)+33]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b33 = _mm_loaddup_pd(&B[(i*84)+33]);
#endif
__m128d c33_0 = _mm_load_sd(&C[(i*84)+2]);
__m128d a33_0 = _mm_load_sd(&values[162]);
#if defined(__SSE3__) && defined(__AVX__)
c33_0 = _mm_add_sd(c33_0, _mm_mul_sd(a33_0, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_0 = _mm_add_sd(c33_0, _mm_mul_sd(a33_0, b33));
#endif
_mm_store_sd(&C[(i*84)+2], c33_0);
__m128d c33_1 = _mm_load_sd(&C[(i*84)+8]);
__m128d a33_1 = _mm_load_sd(&values[163]);
#if defined(__SSE3__) && defined(__AVX__)
c33_1 = _mm_add_sd(c33_1, _mm_mul_sd(a33_1, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_1 = _mm_add_sd(c33_1, _mm_mul_sd(a33_1, b33));
#endif
_mm_store_sd(&C[(i*84)+8], c33_1);
__m128d c33_2 = _mm_load_sd(&C[(i*84)+18]);
__m128d a33_2 = _mm_load_sd(&values[164]);
#if defined(__SSE3__) && defined(__AVX__)
c33_2 = _mm_add_sd(c33_2, _mm_mul_sd(a33_2, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_2 = _mm_add_sd(c33_2, _mm_mul_sd(a33_2, b33));
#endif
_mm_store_sd(&C[(i*84)+18], c33_2);
__m128d c33_3 = _mm_load_sd(&C[(i*84)+33]);
__m128d a33_3 = _mm_load_sd(&values[165]);
#if defined(__SSE3__) && defined(__AVX__)
c33_3 = _mm_add_sd(c33_3, _mm_mul_sd(a33_3, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_3 = _mm_add_sd(c33_3, _mm_mul_sd(a33_3, b33));
#endif
_mm_store_sd(&C[(i*84)+33], c33_3);
__m128d c33_4 = _mm_load_sd(&C[(i*84)+54]);
__m128d a33_4 = _mm_load_sd(&values[166]);
#if defined(__SSE3__) && defined(__AVX__)
c33_4 = _mm_add_sd(c33_4, _mm_mul_sd(a33_4, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_4 = _mm_add_sd(c33_4, _mm_mul_sd(a33_4, b33));
#endif
_mm_store_sd(&C[(i*84)+54], c33_4);
__m128d c33_5 = _mm_load_sd(&C[(i*84)+82]);
__m128d a33_5 = _mm_load_sd(&values[167]);
#if defined(__SSE3__) && defined(__AVX__)
c33_5 = _mm_add_sd(c33_5, _mm_mul_sd(a33_5, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_5 = _mm_add_sd(c33_5, _mm_mul_sd(a33_5, b33));
#endif
_mm_store_sd(&C[(i*84)+82], c33_5);
#else
C[(i*84)+2] += values[162] * B[(i*84)+33];
C[(i*84)+8] += values[163] * B[(i*84)+33];
C[(i*84)+18] += values[164] * B[(i*84)+33];
C[(i*84)+33] += values[165] * B[(i*84)+33];
C[(i*84)+54] += values[166] * B[(i*84)+33];
C[(i*84)+82] += values[167] * B[(i*84)+33];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b34 = _mm256_broadcast_sd(&B[(i*84)+34]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b34 = _mm_loaddup_pd(&B[(i*84)+34]);
#endif
__m128d c34_0 = _mm_load_sd(&C[(i*84)+0]);
__m128d a34_0 = _mm_load_sd(&values[168]);
#if defined(__SSE3__) && defined(__AVX__)
c34_0 = _mm_add_sd(c34_0, _mm_mul_sd(a34_0, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_0 = _mm_add_sd(c34_0, _mm_mul_sd(a34_0, b34));
#endif
_mm_store_sd(&C[(i*84)+0], c34_0);
__m128d c34_1 = _mm_load_sd(&C[(i*84)+3]);
__m128d a34_1 = _mm_load_sd(&values[169]);
#if defined(__SSE3__) && defined(__AVX__)
c34_1 = _mm_add_sd(c34_1, _mm_mul_sd(a34_1, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_1 = _mm_add_sd(c34_1, _mm_mul_sd(a34_1, b34));
#endif
_mm_store_sd(&C[(i*84)+3], c34_1);
__m128d c34_2 = _mm_load_sd(&C[(i*84)+9]);
__m128d a34_2 = _mm_load_sd(&values[170]);
#if defined(__SSE3__) && defined(__AVX__)
c34_2 = _mm_add_sd(c34_2, _mm_mul_sd(a34_2, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_2 = _mm_add_sd(c34_2, _mm_mul_sd(a34_2, b34));
#endif
_mm_store_sd(&C[(i*84)+9], c34_2);
__m128d c34_3 = _mm_load_sd(&C[(i*84)+19]);
__m128d a34_3 = _mm_load_sd(&values[171]);
#if defined(__SSE3__) && defined(__AVX__)
c34_3 = _mm_add_sd(c34_3, _mm_mul_sd(a34_3, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_3 = _mm_add_sd(c34_3, _mm_mul_sd(a34_3, b34));
#endif
_mm_store_sd(&C[(i*84)+19], c34_3);
__m128d c34_4 = _mm_load_sd(&C[(i*84)+34]);
__m128d a34_4 = _mm_load_sd(&values[172]);
#if defined(__SSE3__) && defined(__AVX__)
c34_4 = _mm_add_sd(c34_4, _mm_mul_sd(a34_4, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_4 = _mm_add_sd(c34_4, _mm_mul_sd(a34_4, b34));
#endif
_mm_store_sd(&C[(i*84)+34], c34_4);
__m128d c34_5 = _mm_load_sd(&C[(i*84)+55]);
__m128d a34_5 = _mm_load_sd(&values[173]);
#if defined(__SSE3__) && defined(__AVX__)
c34_5 = _mm_add_sd(c34_5, _mm_mul_sd(a34_5, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_5 = _mm_add_sd(c34_5, _mm_mul_sd(a34_5, b34));
#endif
_mm_store_sd(&C[(i*84)+55], c34_5);
__m128d c34_6 = _mm_load_sd(&C[(i*84)+83]);
__m128d a34_6 = _mm_load_sd(&values[174]);
#if defined(__SSE3__) && defined(__AVX__)
c34_6 = _mm_add_sd(c34_6, _mm_mul_sd(a34_6, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_6 = _mm_add_sd(c34_6, _mm_mul_sd(a34_6, b34));
#endif
_mm_store_sd(&C[(i*84)+83], c34_6);
#else
C[(i*84)+0] += values[168] * B[(i*84)+34];
C[(i*84)+3] += values[169] * B[(i*84)+34];
C[(i*84)+9] += values[170] * B[(i*84)+34];
C[(i*84)+19] += values[171] * B[(i*84)+34];
C[(i*84)+34] += values[172] * B[(i*84)+34];
C[(i*84)+55] += values[173] * B[(i*84)+34];
C[(i*84)+83] += values[174] * B[(i*84)+34];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b35 = _mm256_broadcast_sd(&B[(i*84)+35]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b35 = _mm_loaddup_pd(&B[(i*84)+35]);
#endif
__m128d c35_0 = _mm_load_sd(&C[(i*84)+35]);
__m128d a35_0 = _mm_load_sd(&values[175]);
#if defined(__SSE3__) && defined(__AVX__)
c35_0 = _mm_add_sd(c35_0, _mm_mul_sd(a35_0, _mm256_castpd256_pd128(b35)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c35_0 = _mm_add_sd(c35_0, _mm_mul_sd(a35_0, b35));
#endif
_mm_store_sd(&C[(i*84)+35], c35_0);
__m128d c35_1 = _mm_load_sd(&C[(i*84)+63]);
__m128d a35_1 = _mm_load_sd(&values[176]);
#if defined(__SSE3__) && defined(__AVX__)
c35_1 = _mm_add_sd(c35_1, _mm_mul_sd(a35_1, _mm256_castpd256_pd128(b35)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c35_1 = _mm_add_sd(c35_1, _mm_mul_sd(a35_1, b35));
#endif
_mm_store_sd(&C[(i*84)+63], c35_1);
#else
C[(i*84)+35] += values[175] * B[(i*84)+35];
C[(i*84)+63] += values[176] * B[(i*84)+35];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b36 = _mm256_broadcast_sd(&B[(i*84)+36]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b36 = _mm_loaddup_pd(&B[(i*84)+36]);
#endif
__m128d c36_0 = _mm_load_sd(&C[(i*84)+36]);
__m128d a36_0 = _mm_load_sd(&values[177]);
#if defined(__SSE3__) && defined(__AVX__)
c36_0 = _mm_add_sd(c36_0, _mm_mul_sd(a36_0, _mm256_castpd256_pd128(b36)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c36_0 = _mm_add_sd(c36_0, _mm_mul_sd(a36_0, b36));
#endif
_mm_store_sd(&C[(i*84)+36], c36_0);
__m128d c36_1 = _mm_load_sd(&C[(i*84)+64]);
__m128d a36_1 = _mm_load_sd(&values[178]);
#if defined(__SSE3__) && defined(__AVX__)
c36_1 = _mm_add_sd(c36_1, _mm_mul_sd(a36_1, _mm256_castpd256_pd128(b36)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c36_1 = _mm_add_sd(c36_1, _mm_mul_sd(a36_1, b36));
#endif
_mm_store_sd(&C[(i*84)+64], c36_1);
#else
C[(i*84)+36] += values[177] * B[(i*84)+36];
C[(i*84)+64] += values[178] * B[(i*84)+36];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b37 = _mm256_broadcast_sd(&B[(i*84)+37]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b37 = _mm_loaddup_pd(&B[(i*84)+37]);
#endif
__m128d c37_0 = _mm_load_sd(&C[(i*84)+37]);
__m128d a37_0 = _mm_load_sd(&values[179]);
#if defined(__SSE3__) && defined(__AVX__)
c37_0 = _mm_add_sd(c37_0, _mm_mul_sd(a37_0, _mm256_castpd256_pd128(b37)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c37_0 = _mm_add_sd(c37_0, _mm_mul_sd(a37_0, b37));
#endif
_mm_store_sd(&C[(i*84)+37], c37_0);
__m128d c37_1 = _mm_load_sd(&C[(i*84)+65]);
__m128d a37_1 = _mm_load_sd(&values[180]);
#if defined(__SSE3__) && defined(__AVX__)
c37_1 = _mm_add_sd(c37_1, _mm_mul_sd(a37_1, _mm256_castpd256_pd128(b37)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c37_1 = _mm_add_sd(c37_1, _mm_mul_sd(a37_1, b37));
#endif
_mm_store_sd(&C[(i*84)+65], c37_1);
#else
C[(i*84)+37] += values[179] * B[(i*84)+37];
C[(i*84)+65] += values[180] * B[(i*84)+37];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b38 = _mm256_broadcast_sd(&B[(i*84)+38]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b38 = _mm_loaddup_pd(&B[(i*84)+38]);
#endif
__m128d c38_0 = _mm_load_sd(&C[(i*84)+38]);
__m128d a38_0 = _mm_load_sd(&values[181]);
#if defined(__SSE3__) && defined(__AVX__)
c38_0 = _mm_add_sd(c38_0, _mm_mul_sd(a38_0, _mm256_castpd256_pd128(b38)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c38_0 = _mm_add_sd(c38_0, _mm_mul_sd(a38_0, b38));
#endif
_mm_store_sd(&C[(i*84)+38], c38_0);
__m128d c38_1 = _mm_load_sd(&C[(i*84)+66]);
__m128d a38_1 = _mm_load_sd(&values[182]);
#if defined(__SSE3__) && defined(__AVX__)
c38_1 = _mm_add_sd(c38_1, _mm_mul_sd(a38_1, _mm256_castpd256_pd128(b38)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c38_1 = _mm_add_sd(c38_1, _mm_mul_sd(a38_1, b38));
#endif
_mm_store_sd(&C[(i*84)+66], c38_1);
#else
C[(i*84)+38] += values[181] * B[(i*84)+38];
C[(i*84)+66] += values[182] * B[(i*84)+38];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b39 = _mm256_broadcast_sd(&B[(i*84)+39]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b39 = _mm_loaddup_pd(&B[(i*84)+39]);
#endif
__m128d c39_0 = _mm_load_sd(&C[(i*84)+39]);
__m128d a39_0 = _mm_load_sd(&values[183]);
#if defined(__SSE3__) && defined(__AVX__)
c39_0 = _mm_add_sd(c39_0, _mm_mul_sd(a39_0, _mm256_castpd256_pd128(b39)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c39_0 = _mm_add_sd(c39_0, _mm_mul_sd(a39_0, b39));
#endif
_mm_store_sd(&C[(i*84)+39], c39_0);
__m128d c39_1 = _mm_load_sd(&C[(i*84)+67]);
__m128d a39_1 = _mm_load_sd(&values[184]);
#if defined(__SSE3__) && defined(__AVX__)
c39_1 = _mm_add_sd(c39_1, _mm_mul_sd(a39_1, _mm256_castpd256_pd128(b39)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c39_1 = _mm_add_sd(c39_1, _mm_mul_sd(a39_1, b39));
#endif
_mm_store_sd(&C[(i*84)+67], c39_1);
#else
C[(i*84)+39] += values[183] * B[(i*84)+39];
C[(i*84)+67] += values[184] * B[(i*84)+39];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b40 = _mm256_broadcast_sd(&B[(i*84)+40]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b40 = _mm_loaddup_pd(&B[(i*84)+40]);
#endif
__m128d c40_0 = _mm_load_sd(&C[(i*84)+40]);
__m128d a40_0 = _mm_load_sd(&values[185]);
#if defined(__SSE3__) && defined(__AVX__)
c40_0 = _mm_add_sd(c40_0, _mm_mul_sd(a40_0, _mm256_castpd256_pd128(b40)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c40_0 = _mm_add_sd(c40_0, _mm_mul_sd(a40_0, b40));
#endif
_mm_store_sd(&C[(i*84)+40], c40_0);
__m128d c40_1 = _mm_load_sd(&C[(i*84)+68]);
__m128d a40_1 = _mm_load_sd(&values[186]);
#if defined(__SSE3__) && defined(__AVX__)
c40_1 = _mm_add_sd(c40_1, _mm_mul_sd(a40_1, _mm256_castpd256_pd128(b40)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c40_1 = _mm_add_sd(c40_1, _mm_mul_sd(a40_1, b40));
#endif
_mm_store_sd(&C[(i*84)+68], c40_1);
#else
C[(i*84)+40] += values[185] * B[(i*84)+40];
C[(i*84)+68] += values[186] * B[(i*84)+40];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b41 = _mm256_broadcast_sd(&B[(i*84)+41]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b41 = _mm_loaddup_pd(&B[(i*84)+41]);
#endif
__m128d c41_0 = _mm_load_sd(&C[(i*84)+20]);
__m128d a41_0 = _mm_load_sd(&values[187]);
#if defined(__SSE3__) && defined(__AVX__)
c41_0 = _mm_add_sd(c41_0, _mm_mul_sd(a41_0, _mm256_castpd256_pd128(b41)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c41_0 = _mm_add_sd(c41_0, _mm_mul_sd(a41_0, b41));
#endif
_mm_store_sd(&C[(i*84)+20], c41_0);
__m128d c41_1 = _mm_load_sd(&C[(i*84)+41]);
__m128d a41_1 = _mm_load_sd(&values[188]);
#if defined(__SSE3__) && defined(__AVX__)
c41_1 = _mm_add_sd(c41_1, _mm_mul_sd(a41_1, _mm256_castpd256_pd128(b41)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c41_1 = _mm_add_sd(c41_1, _mm_mul_sd(a41_1, b41));
#endif
_mm_store_sd(&C[(i*84)+41], c41_1);
__m128d c41_2 = _mm_load_sd(&C[(i*84)+69]);
__m128d a41_2 = _mm_load_sd(&values[189]);
#if defined(__SSE3__) && defined(__AVX__)
c41_2 = _mm_add_sd(c41_2, _mm_mul_sd(a41_2, _mm256_castpd256_pd128(b41)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c41_2 = _mm_add_sd(c41_2, _mm_mul_sd(a41_2, b41));
#endif
_mm_store_sd(&C[(i*84)+69], c41_2);
#else
C[(i*84)+20] += values[187] * B[(i*84)+41];
C[(i*84)+41] += values[188] * B[(i*84)+41];
C[(i*84)+69] += values[189] * B[(i*84)+41];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b42 = _mm256_broadcast_sd(&B[(i*84)+42]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b42 = _mm_loaddup_pd(&B[(i*84)+42]);
#endif
__m128d c42_0 = _mm_load_sd(&C[(i*84)+21]);
__m128d a42_0 = _mm_load_sd(&values[190]);
#if defined(__SSE3__) && defined(__AVX__)
c42_0 = _mm_add_sd(c42_0, _mm_mul_sd(a42_0, _mm256_castpd256_pd128(b42)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c42_0 = _mm_add_sd(c42_0, _mm_mul_sd(a42_0, b42));
#endif
_mm_store_sd(&C[(i*84)+21], c42_0);
__m128d c42_1 = _mm_load_sd(&C[(i*84)+42]);
__m128d a42_1 = _mm_load_sd(&values[191]);
#if defined(__SSE3__) && defined(__AVX__)
c42_1 = _mm_add_sd(c42_1, _mm_mul_sd(a42_1, _mm256_castpd256_pd128(b42)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c42_1 = _mm_add_sd(c42_1, _mm_mul_sd(a42_1, b42));
#endif
_mm_store_sd(&C[(i*84)+42], c42_1);
__m128d c42_2 = _mm_load_sd(&C[(i*84)+70]);
__m128d a42_2 = _mm_load_sd(&values[192]);
#if defined(__SSE3__) && defined(__AVX__)
c42_2 = _mm_add_sd(c42_2, _mm_mul_sd(a42_2, _mm256_castpd256_pd128(b42)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c42_2 = _mm_add_sd(c42_2, _mm_mul_sd(a42_2, b42));
#endif
_mm_store_sd(&C[(i*84)+70], c42_2);
#else
C[(i*84)+21] += values[190] * B[(i*84)+42];
C[(i*84)+42] += values[191] * B[(i*84)+42];
C[(i*84)+70] += values[192] * B[(i*84)+42];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b43 = _mm256_broadcast_sd(&B[(i*84)+43]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b43 = _mm_loaddup_pd(&B[(i*84)+43]);
#endif
__m128d c43_0 = _mm_load_sd(&C[(i*84)+22]);
__m128d a43_0 = _mm_load_sd(&values[193]);
#if defined(__SSE3__) && defined(__AVX__)
c43_0 = _mm_add_sd(c43_0, _mm_mul_sd(a43_0, _mm256_castpd256_pd128(b43)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c43_0 = _mm_add_sd(c43_0, _mm_mul_sd(a43_0, b43));
#endif
_mm_store_sd(&C[(i*84)+22], c43_0);
__m128d c43_1 = _mm_load_sd(&C[(i*84)+43]);
__m128d a43_1 = _mm_load_sd(&values[194]);
#if defined(__SSE3__) && defined(__AVX__)
c43_1 = _mm_add_sd(c43_1, _mm_mul_sd(a43_1, _mm256_castpd256_pd128(b43)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c43_1 = _mm_add_sd(c43_1, _mm_mul_sd(a43_1, b43));
#endif
_mm_store_sd(&C[(i*84)+43], c43_1);
__m128d c43_2 = _mm_load_sd(&C[(i*84)+71]);
__m128d a43_2 = _mm_load_sd(&values[195]);
#if defined(__SSE3__) && defined(__AVX__)
c43_2 = _mm_add_sd(c43_2, _mm_mul_sd(a43_2, _mm256_castpd256_pd128(b43)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c43_2 = _mm_add_sd(c43_2, _mm_mul_sd(a43_2, b43));
#endif
_mm_store_sd(&C[(i*84)+71], c43_2);
#else
C[(i*84)+22] += values[193] * B[(i*84)+43];
C[(i*84)+43] += values[194] * B[(i*84)+43];
C[(i*84)+71] += values[195] * B[(i*84)+43];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b44 = _mm256_broadcast_sd(&B[(i*84)+44]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b44 = _mm_loaddup_pd(&B[(i*84)+44]);
#endif
__m128d c44_0 = _mm_load_sd(&C[(i*84)+23]);
__m128d a44_0 = _mm_load_sd(&values[196]);
#if defined(__SSE3__) && defined(__AVX__)
c44_0 = _mm_add_sd(c44_0, _mm_mul_sd(a44_0, _mm256_castpd256_pd128(b44)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c44_0 = _mm_add_sd(c44_0, _mm_mul_sd(a44_0, b44));
#endif
_mm_store_sd(&C[(i*84)+23], c44_0);
__m128d c44_1 = _mm_load_sd(&C[(i*84)+44]);
__m128d a44_1 = _mm_load_sd(&values[197]);
#if defined(__SSE3__) && defined(__AVX__)
c44_1 = _mm_add_sd(c44_1, _mm_mul_sd(a44_1, _mm256_castpd256_pd128(b44)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c44_1 = _mm_add_sd(c44_1, _mm_mul_sd(a44_1, b44));
#endif
_mm_store_sd(&C[(i*84)+44], c44_1);
__m128d c44_2 = _mm_load_sd(&C[(i*84)+72]);
__m128d a44_2 = _mm_load_sd(&values[198]);
#if defined(__SSE3__) && defined(__AVX__)
c44_2 = _mm_add_sd(c44_2, _mm_mul_sd(a44_2, _mm256_castpd256_pd128(b44)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c44_2 = _mm_add_sd(c44_2, _mm_mul_sd(a44_2, b44));
#endif
_mm_store_sd(&C[(i*84)+72], c44_2);
#else
C[(i*84)+23] += values[196] * B[(i*84)+44];
C[(i*84)+44] += values[197] * B[(i*84)+44];
C[(i*84)+72] += values[198] * B[(i*84)+44];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b45 = _mm256_broadcast_sd(&B[(i*84)+45]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b45 = _mm_loaddup_pd(&B[(i*84)+45]);
#endif
__m128d c45_0 = _mm_load_sd(&C[(i*84)+24]);
__m128d a45_0 = _mm_load_sd(&values[199]);
#if defined(__SSE3__) && defined(__AVX__)
c45_0 = _mm_add_sd(c45_0, _mm_mul_sd(a45_0, _mm256_castpd256_pd128(b45)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c45_0 = _mm_add_sd(c45_0, _mm_mul_sd(a45_0, b45));
#endif
_mm_store_sd(&C[(i*84)+24], c45_0);
__m128d c45_1 = _mm_load_sd(&C[(i*84)+45]);
__m128d a45_1 = _mm_load_sd(&values[200]);
#if defined(__SSE3__) && defined(__AVX__)
c45_1 = _mm_add_sd(c45_1, _mm_mul_sd(a45_1, _mm256_castpd256_pd128(b45)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c45_1 = _mm_add_sd(c45_1, _mm_mul_sd(a45_1, b45));
#endif
_mm_store_sd(&C[(i*84)+45], c45_1);
__m128d c45_2 = _mm_load_sd(&C[(i*84)+73]);
__m128d a45_2 = _mm_load_sd(&values[201]);
#if defined(__SSE3__) && defined(__AVX__)
c45_2 = _mm_add_sd(c45_2, _mm_mul_sd(a45_2, _mm256_castpd256_pd128(b45)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c45_2 = _mm_add_sd(c45_2, _mm_mul_sd(a45_2, b45));
#endif
_mm_store_sd(&C[(i*84)+73], c45_2);
#else
C[(i*84)+24] += values[199] * B[(i*84)+45];
C[(i*84)+45] += values[200] * B[(i*84)+45];
C[(i*84)+73] += values[201] * B[(i*84)+45];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b46 = _mm256_broadcast_sd(&B[(i*84)+46]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b46 = _mm_loaddup_pd(&B[(i*84)+46]);
#endif
__m128d c46_0 = _mm_load_sd(&C[(i*84)+10]);
__m128d a46_0 = _mm_load_sd(&values[202]);
#if defined(__SSE3__) && defined(__AVX__)
c46_0 = _mm_add_sd(c46_0, _mm_mul_sd(a46_0, _mm256_castpd256_pd128(b46)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c46_0 = _mm_add_sd(c46_0, _mm_mul_sd(a46_0, b46));
#endif
_mm_store_sd(&C[(i*84)+10], c46_0);
__m128d c46_1 = _mm_load_sd(&C[(i*84)+25]);
__m128d a46_1 = _mm_load_sd(&values[203]);
#if defined(__SSE3__) && defined(__AVX__)
c46_1 = _mm_add_sd(c46_1, _mm_mul_sd(a46_1, _mm256_castpd256_pd128(b46)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c46_1 = _mm_add_sd(c46_1, _mm_mul_sd(a46_1, b46));
#endif
_mm_store_sd(&C[(i*84)+25], c46_1);
__m128d c46_2 = _mm_load_sd(&C[(i*84)+46]);
__m128d a46_2 = _mm_load_sd(&values[204]);
#if defined(__SSE3__) && defined(__AVX__)
c46_2 = _mm_add_sd(c46_2, _mm_mul_sd(a46_2, _mm256_castpd256_pd128(b46)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c46_2 = _mm_add_sd(c46_2, _mm_mul_sd(a46_2, b46));
#endif
_mm_store_sd(&C[(i*84)+46], c46_2);
__m128d c46_3 = _mm_load_sd(&C[(i*84)+74]);
__m128d a46_3 = _mm_load_sd(&values[205]);
#if defined(__SSE3__) && defined(__AVX__)
c46_3 = _mm_add_sd(c46_3, _mm_mul_sd(a46_3, _mm256_castpd256_pd128(b46)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c46_3 = _mm_add_sd(c46_3, _mm_mul_sd(a46_3, b46));
#endif
_mm_store_sd(&C[(i*84)+74], c46_3);
#else
C[(i*84)+10] += values[202] * B[(i*84)+46];
C[(i*84)+25] += values[203] * B[(i*84)+46];
C[(i*84)+46] += values[204] * B[(i*84)+46];
C[(i*84)+74] += values[205] * B[(i*84)+46];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b47 = _mm256_broadcast_sd(&B[(i*84)+47]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b47 = _mm_loaddup_pd(&B[(i*84)+47]);
#endif
__m128d c47_0 = _mm_load_sd(&C[(i*84)+11]);
__m128d a47_0 = _mm_load_sd(&values[206]);
#if defined(__SSE3__) && defined(__AVX__)
c47_0 = _mm_add_sd(c47_0, _mm_mul_sd(a47_0, _mm256_castpd256_pd128(b47)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c47_0 = _mm_add_sd(c47_0, _mm_mul_sd(a47_0, b47));
#endif
_mm_store_sd(&C[(i*84)+11], c47_0);
__m128d c47_1 = _mm_load_sd(&C[(i*84)+26]);
__m128d a47_1 = _mm_load_sd(&values[207]);
#if defined(__SSE3__) && defined(__AVX__)
c47_1 = _mm_add_sd(c47_1, _mm_mul_sd(a47_1, _mm256_castpd256_pd128(b47)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c47_1 = _mm_add_sd(c47_1, _mm_mul_sd(a47_1, b47));
#endif
_mm_store_sd(&C[(i*84)+26], c47_1);
__m128d c47_2 = _mm_load_sd(&C[(i*84)+47]);
__m128d a47_2 = _mm_load_sd(&values[208]);
#if defined(__SSE3__) && defined(__AVX__)
c47_2 = _mm_add_sd(c47_2, _mm_mul_sd(a47_2, _mm256_castpd256_pd128(b47)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c47_2 = _mm_add_sd(c47_2, _mm_mul_sd(a47_2, b47));
#endif
_mm_store_sd(&C[(i*84)+47], c47_2);
__m128d c47_3 = _mm_load_sd(&C[(i*84)+75]);
__m128d a47_3 = _mm_load_sd(&values[209]);
#if defined(__SSE3__) && defined(__AVX__)
c47_3 = _mm_add_sd(c47_3, _mm_mul_sd(a47_3, _mm256_castpd256_pd128(b47)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c47_3 = _mm_add_sd(c47_3, _mm_mul_sd(a47_3, b47));
#endif
_mm_store_sd(&C[(i*84)+75], c47_3);
#else
C[(i*84)+11] += values[206] * B[(i*84)+47];
C[(i*84)+26] += values[207] * B[(i*84)+47];
C[(i*84)+47] += values[208] * B[(i*84)+47];
C[(i*84)+75] += values[209] * B[(i*84)+47];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b48 = _mm256_broadcast_sd(&B[(i*84)+48]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b48 = _mm_loaddup_pd(&B[(i*84)+48]);
#endif
__m128d c48_0 = _mm_load_sd(&C[(i*84)+12]);
__m128d a48_0 = _mm_load_sd(&values[210]);
#if defined(__SSE3__) && defined(__AVX__)
c48_0 = _mm_add_sd(c48_0, _mm_mul_sd(a48_0, _mm256_castpd256_pd128(b48)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c48_0 = _mm_add_sd(c48_0, _mm_mul_sd(a48_0, b48));
#endif
_mm_store_sd(&C[(i*84)+12], c48_0);
__m128d c48_1 = _mm_load_sd(&C[(i*84)+27]);
__m128d a48_1 = _mm_load_sd(&values[211]);
#if defined(__SSE3__) && defined(__AVX__)
c48_1 = _mm_add_sd(c48_1, _mm_mul_sd(a48_1, _mm256_castpd256_pd128(b48)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c48_1 = _mm_add_sd(c48_1, _mm_mul_sd(a48_1, b48));
#endif
_mm_store_sd(&C[(i*84)+27], c48_1);
__m128d c48_2 = _mm_load_sd(&C[(i*84)+48]);
__m128d a48_2 = _mm_load_sd(&values[212]);
#if defined(__SSE3__) && defined(__AVX__)
c48_2 = _mm_add_sd(c48_2, _mm_mul_sd(a48_2, _mm256_castpd256_pd128(b48)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c48_2 = _mm_add_sd(c48_2, _mm_mul_sd(a48_2, b48));
#endif
_mm_store_sd(&C[(i*84)+48], c48_2);
__m128d c48_3 = _mm_load_sd(&C[(i*84)+76]);
__m128d a48_3 = _mm_load_sd(&values[213]);
#if defined(__SSE3__) && defined(__AVX__)
c48_3 = _mm_add_sd(c48_3, _mm_mul_sd(a48_3, _mm256_castpd256_pd128(b48)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c48_3 = _mm_add_sd(c48_3, _mm_mul_sd(a48_3, b48));
#endif
_mm_store_sd(&C[(i*84)+76], c48_3);
#else
C[(i*84)+12] += values[210] * B[(i*84)+48];
C[(i*84)+27] += values[211] * B[(i*84)+48];
C[(i*84)+48] += values[212] * B[(i*84)+48];
C[(i*84)+76] += values[213] * B[(i*84)+48];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b49 = _mm256_broadcast_sd(&B[(i*84)+49]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b49 = _mm_loaddup_pd(&B[(i*84)+49]);
#endif
__m128d c49_0 = _mm_load_sd(&C[(i*84)+13]);
__m128d a49_0 = _mm_load_sd(&values[214]);
#if defined(__SSE3__) && defined(__AVX__)
c49_0 = _mm_add_sd(c49_0, _mm_mul_sd(a49_0, _mm256_castpd256_pd128(b49)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c49_0 = _mm_add_sd(c49_0, _mm_mul_sd(a49_0, b49));
#endif
_mm_store_sd(&C[(i*84)+13], c49_0);
__m128d c49_1 = _mm_load_sd(&C[(i*84)+28]);
__m128d a49_1 = _mm_load_sd(&values[215]);
#if defined(__SSE3__) && defined(__AVX__)
c49_1 = _mm_add_sd(c49_1, _mm_mul_sd(a49_1, _mm256_castpd256_pd128(b49)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c49_1 = _mm_add_sd(c49_1, _mm_mul_sd(a49_1, b49));
#endif
_mm_store_sd(&C[(i*84)+28], c49_1);
__m128d c49_2 = _mm_load_sd(&C[(i*84)+49]);
__m128d a49_2 = _mm_load_sd(&values[216]);
#if defined(__SSE3__) && defined(__AVX__)
c49_2 = _mm_add_sd(c49_2, _mm_mul_sd(a49_2, _mm256_castpd256_pd128(b49)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c49_2 = _mm_add_sd(c49_2, _mm_mul_sd(a49_2, b49));
#endif
_mm_store_sd(&C[(i*84)+49], c49_2);
__m128d c49_3 = _mm_load_sd(&C[(i*84)+77]);
__m128d a49_3 = _mm_load_sd(&values[217]);
#if defined(__SSE3__) && defined(__AVX__)
c49_3 = _mm_add_sd(c49_3, _mm_mul_sd(a49_3, _mm256_castpd256_pd128(b49)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c49_3 = _mm_add_sd(c49_3, _mm_mul_sd(a49_3, b49));
#endif
_mm_store_sd(&C[(i*84)+77], c49_3);
#else
C[(i*84)+13] += values[214] * B[(i*84)+49];
C[(i*84)+28] += values[215] * B[(i*84)+49];
C[(i*84)+49] += values[216] * B[(i*84)+49];
C[(i*84)+77] += values[217] * B[(i*84)+49];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b50 = _mm256_broadcast_sd(&B[(i*84)+50]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b50 = _mm_loaddup_pd(&B[(i*84)+50]);
#endif
__m128d c50_0 = _mm_load_sd(&C[(i*84)+4]);
__m128d a50_0 = _mm_load_sd(&values[218]);
#if defined(__SSE3__) && defined(__AVX__)
c50_0 = _mm_add_sd(c50_0, _mm_mul_sd(a50_0, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c50_0 = _mm_add_sd(c50_0, _mm_mul_sd(a50_0, b50));
#endif
_mm_store_sd(&C[(i*84)+4], c50_0);
__m128d c50_1 = _mm_load_sd(&C[(i*84)+14]);
__m128d a50_1 = _mm_load_sd(&values[219]);
#if defined(__SSE3__) && defined(__AVX__)
c50_1 = _mm_add_sd(c50_1, _mm_mul_sd(a50_1, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c50_1 = _mm_add_sd(c50_1, _mm_mul_sd(a50_1, b50));
#endif
_mm_store_sd(&C[(i*84)+14], c50_1);
__m128d c50_2 = _mm_load_sd(&C[(i*84)+29]);
__m128d a50_2 = _mm_load_sd(&values[220]);
#if defined(__SSE3__) && defined(__AVX__)
c50_2 = _mm_add_sd(c50_2, _mm_mul_sd(a50_2, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c50_2 = _mm_add_sd(c50_2, _mm_mul_sd(a50_2, b50));
#endif
_mm_store_sd(&C[(i*84)+29], c50_2);
__m128d c50_3 = _mm_load_sd(&C[(i*84)+50]);
__m128d a50_3 = _mm_load_sd(&values[221]);
#if defined(__SSE3__) && defined(__AVX__)
c50_3 = _mm_add_sd(c50_3, _mm_mul_sd(a50_3, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c50_3 = _mm_add_sd(c50_3, _mm_mul_sd(a50_3, b50));
#endif
_mm_store_sd(&C[(i*84)+50], c50_3);
__m128d c50_4 = _mm_load_sd(&C[(i*84)+78]);
__m128d a50_4 = _mm_load_sd(&values[222]);
#if defined(__SSE3__) && defined(__AVX__)
c50_4 = _mm_add_sd(c50_4, _mm_mul_sd(a50_4, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c50_4 = _mm_add_sd(c50_4, _mm_mul_sd(a50_4, b50));
#endif
_mm_store_sd(&C[(i*84)+78], c50_4);
#else
C[(i*84)+4] += values[218] * B[(i*84)+50];
C[(i*84)+14] += values[219] * B[(i*84)+50];
C[(i*84)+29] += values[220] * B[(i*84)+50];
C[(i*84)+50] += values[221] * B[(i*84)+50];
C[(i*84)+78] += values[222] * B[(i*84)+50];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b51 = _mm256_broadcast_sd(&B[(i*84)+51]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b51 = _mm_loaddup_pd(&B[(i*84)+51]);
#endif
__m128d c51_0 = _mm_load_sd(&C[(i*84)+5]);
__m128d a51_0 = _mm_load_sd(&values[223]);
#if defined(__SSE3__) && defined(__AVX__)
c51_0 = _mm_add_sd(c51_0, _mm_mul_sd(a51_0, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c51_0 = _mm_add_sd(c51_0, _mm_mul_sd(a51_0, b51));
#endif
_mm_store_sd(&C[(i*84)+5], c51_0);
__m128d c51_1 = _mm_load_sd(&C[(i*84)+15]);
__m128d a51_1 = _mm_load_sd(&values[224]);
#if defined(__SSE3__) && defined(__AVX__)
c51_1 = _mm_add_sd(c51_1, _mm_mul_sd(a51_1, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c51_1 = _mm_add_sd(c51_1, _mm_mul_sd(a51_1, b51));
#endif
_mm_store_sd(&C[(i*84)+15], c51_1);
__m128d c51_2 = _mm_load_sd(&C[(i*84)+30]);
__m128d a51_2 = _mm_load_sd(&values[225]);
#if defined(__SSE3__) && defined(__AVX__)
c51_2 = _mm_add_sd(c51_2, _mm_mul_sd(a51_2, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c51_2 = _mm_add_sd(c51_2, _mm_mul_sd(a51_2, b51));
#endif
_mm_store_sd(&C[(i*84)+30], c51_2);
__m128d c51_3 = _mm_load_sd(&C[(i*84)+51]);
__m128d a51_3 = _mm_load_sd(&values[226]);
#if defined(__SSE3__) && defined(__AVX__)
c51_3 = _mm_add_sd(c51_3, _mm_mul_sd(a51_3, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c51_3 = _mm_add_sd(c51_3, _mm_mul_sd(a51_3, b51));
#endif
_mm_store_sd(&C[(i*84)+51], c51_3);
__m128d c51_4 = _mm_load_sd(&C[(i*84)+79]);
__m128d a51_4 = _mm_load_sd(&values[227]);
#if defined(__SSE3__) && defined(__AVX__)
c51_4 = _mm_add_sd(c51_4, _mm_mul_sd(a51_4, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c51_4 = _mm_add_sd(c51_4, _mm_mul_sd(a51_4, b51));
#endif
_mm_store_sd(&C[(i*84)+79], c51_4);
#else
C[(i*84)+5] += values[223] * B[(i*84)+51];
C[(i*84)+15] += values[224] * B[(i*84)+51];
C[(i*84)+30] += values[225] * B[(i*84)+51];
C[(i*84)+51] += values[226] * B[(i*84)+51];
C[(i*84)+79] += values[227] * B[(i*84)+51];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b52 = _mm256_broadcast_sd(&B[(i*84)+52]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b52 = _mm_loaddup_pd(&B[(i*84)+52]);
#endif
__m128d c52_0 = _mm_load_sd(&C[(i*84)+6]);
__m128d a52_0 = _mm_load_sd(&values[228]);
#if defined(__SSE3__) && defined(__AVX__)
c52_0 = _mm_add_sd(c52_0, _mm_mul_sd(a52_0, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c52_0 = _mm_add_sd(c52_0, _mm_mul_sd(a52_0, b52));
#endif
_mm_store_sd(&C[(i*84)+6], c52_0);
__m128d c52_1 = _mm_load_sd(&C[(i*84)+16]);
__m128d a52_1 = _mm_load_sd(&values[229]);
#if defined(__SSE3__) && defined(__AVX__)
c52_1 = _mm_add_sd(c52_1, _mm_mul_sd(a52_1, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c52_1 = _mm_add_sd(c52_1, _mm_mul_sd(a52_1, b52));
#endif
_mm_store_sd(&C[(i*84)+16], c52_1);
__m128d c52_2 = _mm_load_sd(&C[(i*84)+31]);
__m128d a52_2 = _mm_load_sd(&values[230]);
#if defined(__SSE3__) && defined(__AVX__)
c52_2 = _mm_add_sd(c52_2, _mm_mul_sd(a52_2, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c52_2 = _mm_add_sd(c52_2, _mm_mul_sd(a52_2, b52));
#endif
_mm_store_sd(&C[(i*84)+31], c52_2);
__m128d c52_3 = _mm_load_sd(&C[(i*84)+52]);
__m128d a52_3 = _mm_load_sd(&values[231]);
#if defined(__SSE3__) && defined(__AVX__)
c52_3 = _mm_add_sd(c52_3, _mm_mul_sd(a52_3, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c52_3 = _mm_add_sd(c52_3, _mm_mul_sd(a52_3, b52));
#endif
_mm_store_sd(&C[(i*84)+52], c52_3);
__m128d c52_4 = _mm_load_sd(&C[(i*84)+80]);
__m128d a52_4 = _mm_load_sd(&values[232]);
#if defined(__SSE3__) && defined(__AVX__)
c52_4 = _mm_add_sd(c52_4, _mm_mul_sd(a52_4, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c52_4 = _mm_add_sd(c52_4, _mm_mul_sd(a52_4, b52));
#endif
_mm_store_sd(&C[(i*84)+80], c52_4);
#else
C[(i*84)+6] += values[228] * B[(i*84)+52];
C[(i*84)+16] += values[229] * B[(i*84)+52];
C[(i*84)+31] += values[230] * B[(i*84)+52];
C[(i*84)+52] += values[231] * B[(i*84)+52];
C[(i*84)+80] += values[232] * B[(i*84)+52];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b53 = _mm256_broadcast_sd(&B[(i*84)+53]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b53 = _mm_loaddup_pd(&B[(i*84)+53]);
#endif
__m128d c53_0 = _mm_load_sd(&C[(i*84)+1]);
__m128d a53_0 = _mm_load_sd(&values[233]);
#if defined(__SSE3__) && defined(__AVX__)
c53_0 = _mm_add_sd(c53_0, _mm_mul_sd(a53_0, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c53_0 = _mm_add_sd(c53_0, _mm_mul_sd(a53_0, b53));
#endif
_mm_store_sd(&C[(i*84)+1], c53_0);
__m128d c53_1 = _mm_load_sd(&C[(i*84)+7]);
__m128d a53_1 = _mm_load_sd(&values[234]);
#if defined(__SSE3__) && defined(__AVX__)
c53_1 = _mm_add_sd(c53_1, _mm_mul_sd(a53_1, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c53_1 = _mm_add_sd(c53_1, _mm_mul_sd(a53_1, b53));
#endif
_mm_store_sd(&C[(i*84)+7], c53_1);
__m128d c53_2 = _mm_load_sd(&C[(i*84)+17]);
__m128d a53_2 = _mm_load_sd(&values[235]);
#if defined(__SSE3__) && defined(__AVX__)
c53_2 = _mm_add_sd(c53_2, _mm_mul_sd(a53_2, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c53_2 = _mm_add_sd(c53_2, _mm_mul_sd(a53_2, b53));
#endif
_mm_store_sd(&C[(i*84)+17], c53_2);
__m128d c53_3 = _mm_load_sd(&C[(i*84)+32]);
__m128d a53_3 = _mm_load_sd(&values[236]);
#if defined(__SSE3__) && defined(__AVX__)
c53_3 = _mm_add_sd(c53_3, _mm_mul_sd(a53_3, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c53_3 = _mm_add_sd(c53_3, _mm_mul_sd(a53_3, b53));
#endif
_mm_store_sd(&C[(i*84)+32], c53_3);
__m128d c53_4 = _mm_load_sd(&C[(i*84)+53]);
__m128d a53_4 = _mm_load_sd(&values[237]);
#if defined(__SSE3__) && defined(__AVX__)
c53_4 = _mm_add_sd(c53_4, _mm_mul_sd(a53_4, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c53_4 = _mm_add_sd(c53_4, _mm_mul_sd(a53_4, b53));
#endif
_mm_store_sd(&C[(i*84)+53], c53_4);
__m128d c53_5 = _mm_load_sd(&C[(i*84)+81]);
__m128d a53_5 = _mm_load_sd(&values[238]);
#if defined(__SSE3__) && defined(__AVX__)
c53_5 = _mm_add_sd(c53_5, _mm_mul_sd(a53_5, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c53_5 = _mm_add_sd(c53_5, _mm_mul_sd(a53_5, b53));
#endif
_mm_store_sd(&C[(i*84)+81], c53_5);
#else
C[(i*84)+1] += values[233] * B[(i*84)+53];
C[(i*84)+7] += values[234] * B[(i*84)+53];
C[(i*84)+17] += values[235] * B[(i*84)+53];
C[(i*84)+32] += values[236] * B[(i*84)+53];
C[(i*84)+53] += values[237] * B[(i*84)+53];
C[(i*84)+81] += values[238] * B[(i*84)+53];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b54 = _mm256_broadcast_sd(&B[(i*84)+54]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b54 = _mm_loaddup_pd(&B[(i*84)+54]);
#endif
__m128d c54_0 = _mm_load_sd(&C[(i*84)+2]);
__m128d a54_0 = _mm_load_sd(&values[239]);
#if defined(__SSE3__) && defined(__AVX__)
c54_0 = _mm_add_sd(c54_0, _mm_mul_sd(a54_0, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c54_0 = _mm_add_sd(c54_0, _mm_mul_sd(a54_0, b54));
#endif
_mm_store_sd(&C[(i*84)+2], c54_0);
__m128d c54_1 = _mm_load_sd(&C[(i*84)+8]);
__m128d a54_1 = _mm_load_sd(&values[240]);
#if defined(__SSE3__) && defined(__AVX__)
c54_1 = _mm_add_sd(c54_1, _mm_mul_sd(a54_1, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c54_1 = _mm_add_sd(c54_1, _mm_mul_sd(a54_1, b54));
#endif
_mm_store_sd(&C[(i*84)+8], c54_1);
__m128d c54_2 = _mm_load_sd(&C[(i*84)+18]);
__m128d a54_2 = _mm_load_sd(&values[241]);
#if defined(__SSE3__) && defined(__AVX__)
c54_2 = _mm_add_sd(c54_2, _mm_mul_sd(a54_2, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c54_2 = _mm_add_sd(c54_2, _mm_mul_sd(a54_2, b54));
#endif
_mm_store_sd(&C[(i*84)+18], c54_2);
__m128d c54_3 = _mm_load_sd(&C[(i*84)+33]);
__m128d a54_3 = _mm_load_sd(&values[242]);
#if defined(__SSE3__) && defined(__AVX__)
c54_3 = _mm_add_sd(c54_3, _mm_mul_sd(a54_3, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c54_3 = _mm_add_sd(c54_3, _mm_mul_sd(a54_3, b54));
#endif
_mm_store_sd(&C[(i*84)+33], c54_3);
__m128d c54_4 = _mm_load_sd(&C[(i*84)+54]);
__m128d a54_4 = _mm_load_sd(&values[243]);
#if defined(__SSE3__) && defined(__AVX__)
c54_4 = _mm_add_sd(c54_4, _mm_mul_sd(a54_4, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c54_4 = _mm_add_sd(c54_4, _mm_mul_sd(a54_4, b54));
#endif
_mm_store_sd(&C[(i*84)+54], c54_4);
__m128d c54_5 = _mm_load_sd(&C[(i*84)+82]);
__m128d a54_5 = _mm_load_sd(&values[244]);
#if defined(__SSE3__) && defined(__AVX__)
c54_5 = _mm_add_sd(c54_5, _mm_mul_sd(a54_5, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c54_5 = _mm_add_sd(c54_5, _mm_mul_sd(a54_5, b54));
#endif
_mm_store_sd(&C[(i*84)+82], c54_5);
#else
C[(i*84)+2] += values[239] * B[(i*84)+54];
C[(i*84)+8] += values[240] * B[(i*84)+54];
C[(i*84)+18] += values[241] * B[(i*84)+54];
C[(i*84)+33] += values[242] * B[(i*84)+54];
C[(i*84)+54] += values[243] * B[(i*84)+54];
C[(i*84)+82] += values[244] * B[(i*84)+54];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b55 = _mm256_broadcast_sd(&B[(i*84)+55]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b55 = _mm_loaddup_pd(&B[(i*84)+55]);
#endif
__m128d c55_0 = _mm_load_sd(&C[(i*84)+0]);
__m128d a55_0 = _mm_load_sd(&values[245]);
#if defined(__SSE3__) && defined(__AVX__)
c55_0 = _mm_add_sd(c55_0, _mm_mul_sd(a55_0, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_0 = _mm_add_sd(c55_0, _mm_mul_sd(a55_0, b55));
#endif
_mm_store_sd(&C[(i*84)+0], c55_0);
__m128d c55_1 = _mm_load_sd(&C[(i*84)+3]);
__m128d a55_1 = _mm_load_sd(&values[246]);
#if defined(__SSE3__) && defined(__AVX__)
c55_1 = _mm_add_sd(c55_1, _mm_mul_sd(a55_1, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_1 = _mm_add_sd(c55_1, _mm_mul_sd(a55_1, b55));
#endif
_mm_store_sd(&C[(i*84)+3], c55_1);
__m128d c55_2 = _mm_load_sd(&C[(i*84)+9]);
__m128d a55_2 = _mm_load_sd(&values[247]);
#if defined(__SSE3__) && defined(__AVX__)
c55_2 = _mm_add_sd(c55_2, _mm_mul_sd(a55_2, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_2 = _mm_add_sd(c55_2, _mm_mul_sd(a55_2, b55));
#endif
_mm_store_sd(&C[(i*84)+9], c55_2);
__m128d c55_3 = _mm_load_sd(&C[(i*84)+19]);
__m128d a55_3 = _mm_load_sd(&values[248]);
#if defined(__SSE3__) && defined(__AVX__)
c55_3 = _mm_add_sd(c55_3, _mm_mul_sd(a55_3, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_3 = _mm_add_sd(c55_3, _mm_mul_sd(a55_3, b55));
#endif
_mm_store_sd(&C[(i*84)+19], c55_3);
__m128d c55_4 = _mm_load_sd(&C[(i*84)+34]);
__m128d a55_4 = _mm_load_sd(&values[249]);
#if defined(__SSE3__) && defined(__AVX__)
c55_4 = _mm_add_sd(c55_4, _mm_mul_sd(a55_4, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_4 = _mm_add_sd(c55_4, _mm_mul_sd(a55_4, b55));
#endif
_mm_store_sd(&C[(i*84)+34], c55_4);
__m128d c55_5 = _mm_load_sd(&C[(i*84)+55]);
__m128d a55_5 = _mm_load_sd(&values[250]);
#if defined(__SSE3__) && defined(__AVX__)
c55_5 = _mm_add_sd(c55_5, _mm_mul_sd(a55_5, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_5 = _mm_add_sd(c55_5, _mm_mul_sd(a55_5, b55));
#endif
_mm_store_sd(&C[(i*84)+55], c55_5);
__m128d c55_6 = _mm_load_sd(&C[(i*84)+83]);
__m128d a55_6 = _mm_load_sd(&values[251]);
#if defined(__SSE3__) && defined(__AVX__)
c55_6 = _mm_add_sd(c55_6, _mm_mul_sd(a55_6, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_6 = _mm_add_sd(c55_6, _mm_mul_sd(a55_6, b55));
#endif
_mm_store_sd(&C[(i*84)+83], c55_6);
#else
C[(i*84)+0] += values[245] * B[(i*84)+55];
C[(i*84)+3] += values[246] * B[(i*84)+55];
C[(i*84)+9] += values[247] * B[(i*84)+55];
C[(i*84)+19] += values[248] * B[(i*84)+55];
C[(i*84)+34] += values[249] * B[(i*84)+55];
C[(i*84)+55] += values[250] * B[(i*84)+55];
C[(i*84)+83] += values[251] * B[(i*84)+55];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b56 = _mm256_broadcast_sd(&B[(i*84)+56]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b56 = _mm_loaddup_pd(&B[(i*84)+56]);
#endif
__m128d c56_0 = _mm_load_sd(&C[(i*84)+56]);
__m128d a56_0 = _mm_load_sd(&values[252]);
#if defined(__SSE3__) && defined(__AVX__)
c56_0 = _mm_add_sd(c56_0, _mm_mul_sd(a56_0, _mm256_castpd256_pd128(b56)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c56_0 = _mm_add_sd(c56_0, _mm_mul_sd(a56_0, b56));
#endif
_mm_store_sd(&C[(i*84)+56], c56_0);
#else
C[(i*84)+56] += values[252] * B[(i*84)+56];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b57 = _mm256_broadcast_sd(&B[(i*84)+57]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b57 = _mm_loaddup_pd(&B[(i*84)+57]);
#endif
__m128d c57_0 = _mm_load_sd(&C[(i*84)+57]);
__m128d a57_0 = _mm_load_sd(&values[253]);
#if defined(__SSE3__) && defined(__AVX__)
c57_0 = _mm_add_sd(c57_0, _mm_mul_sd(a57_0, _mm256_castpd256_pd128(b57)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c57_0 = _mm_add_sd(c57_0, _mm_mul_sd(a57_0, b57));
#endif
_mm_store_sd(&C[(i*84)+57], c57_0);
#else
C[(i*84)+57] += values[253] * B[(i*84)+57];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b58 = _mm256_broadcast_sd(&B[(i*84)+58]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b58 = _mm_loaddup_pd(&B[(i*84)+58]);
#endif
__m128d c58_0 = _mm_load_sd(&C[(i*84)+58]);
__m128d a58_0 = _mm_load_sd(&values[254]);
#if defined(__SSE3__) && defined(__AVX__)
c58_0 = _mm_add_sd(c58_0, _mm_mul_sd(a58_0, _mm256_castpd256_pd128(b58)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c58_0 = _mm_add_sd(c58_0, _mm_mul_sd(a58_0, b58));
#endif
_mm_store_sd(&C[(i*84)+58], c58_0);
#else
C[(i*84)+58] += values[254] * B[(i*84)+58];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b59 = _mm256_broadcast_sd(&B[(i*84)+59]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b59 = _mm_loaddup_pd(&B[(i*84)+59]);
#endif
__m128d c59_0 = _mm_load_sd(&C[(i*84)+59]);
__m128d a59_0 = _mm_load_sd(&values[255]);
#if defined(__SSE3__) && defined(__AVX__)
c59_0 = _mm_add_sd(c59_0, _mm_mul_sd(a59_0, _mm256_castpd256_pd128(b59)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c59_0 = _mm_add_sd(c59_0, _mm_mul_sd(a59_0, b59));
#endif
_mm_store_sd(&C[(i*84)+59], c59_0);
#else
C[(i*84)+59] += values[255] * B[(i*84)+59];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b60 = _mm256_broadcast_sd(&B[(i*84)+60]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b60 = _mm_loaddup_pd(&B[(i*84)+60]);
#endif
__m128d c60_0 = _mm_load_sd(&C[(i*84)+60]);
__m128d a60_0 = _mm_load_sd(&values[256]);
#if defined(__SSE3__) && defined(__AVX__)
c60_0 = _mm_add_sd(c60_0, _mm_mul_sd(a60_0, _mm256_castpd256_pd128(b60)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c60_0 = _mm_add_sd(c60_0, _mm_mul_sd(a60_0, b60));
#endif
_mm_store_sd(&C[(i*84)+60], c60_0);
#else
C[(i*84)+60] += values[256] * B[(i*84)+60];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b61 = _mm256_broadcast_sd(&B[(i*84)+61]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b61 = _mm_loaddup_pd(&B[(i*84)+61]);
#endif
__m128d c61_0 = _mm_load_sd(&C[(i*84)+61]);
__m128d a61_0 = _mm_load_sd(&values[257]);
#if defined(__SSE3__) && defined(__AVX__)
c61_0 = _mm_add_sd(c61_0, _mm_mul_sd(a61_0, _mm256_castpd256_pd128(b61)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c61_0 = _mm_add_sd(c61_0, _mm_mul_sd(a61_0, b61));
#endif
_mm_store_sd(&C[(i*84)+61], c61_0);
#else
C[(i*84)+61] += values[257] * B[(i*84)+61];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b62 = _mm256_broadcast_sd(&B[(i*84)+62]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b62 = _mm_loaddup_pd(&B[(i*84)+62]);
#endif
__m128d c62_0 = _mm_load_sd(&C[(i*84)+62]);
__m128d a62_0 = _mm_load_sd(&values[258]);
#if defined(__SSE3__) && defined(__AVX__)
c62_0 = _mm_add_sd(c62_0, _mm_mul_sd(a62_0, _mm256_castpd256_pd128(b62)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c62_0 = _mm_add_sd(c62_0, _mm_mul_sd(a62_0, b62));
#endif
_mm_store_sd(&C[(i*84)+62], c62_0);
#else
C[(i*84)+62] += values[258] * B[(i*84)+62];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b63 = _mm256_broadcast_sd(&B[(i*84)+63]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b63 = _mm_loaddup_pd(&B[(i*84)+63]);
#endif
__m128d c63_0 = _mm_load_sd(&C[(i*84)+35]);
__m128d a63_0 = _mm_load_sd(&values[259]);
#if defined(__SSE3__) && defined(__AVX__)
c63_0 = _mm_add_sd(c63_0, _mm_mul_sd(a63_0, _mm256_castpd256_pd128(b63)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c63_0 = _mm_add_sd(c63_0, _mm_mul_sd(a63_0, b63));
#endif
_mm_store_sd(&C[(i*84)+35], c63_0);
__m128d c63_1 = _mm_load_sd(&C[(i*84)+63]);
__m128d a63_1 = _mm_load_sd(&values[260]);
#if defined(__SSE3__) && defined(__AVX__)
c63_1 = _mm_add_sd(c63_1, _mm_mul_sd(a63_1, _mm256_castpd256_pd128(b63)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c63_1 = _mm_add_sd(c63_1, _mm_mul_sd(a63_1, b63));
#endif
_mm_store_sd(&C[(i*84)+63], c63_1);
#else
C[(i*84)+35] += values[259] * B[(i*84)+63];
C[(i*84)+63] += values[260] * B[(i*84)+63];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b64 = _mm256_broadcast_sd(&B[(i*84)+64]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b64 = _mm_loaddup_pd(&B[(i*84)+64]);
#endif
__m128d c64_0 = _mm_load_sd(&C[(i*84)+36]);
__m128d a64_0 = _mm_load_sd(&values[261]);
#if defined(__SSE3__) && defined(__AVX__)
c64_0 = _mm_add_sd(c64_0, _mm_mul_sd(a64_0, _mm256_castpd256_pd128(b64)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c64_0 = _mm_add_sd(c64_0, _mm_mul_sd(a64_0, b64));
#endif
_mm_store_sd(&C[(i*84)+36], c64_0);
__m128d c64_1 = _mm_load_sd(&C[(i*84)+64]);
__m128d a64_1 = _mm_load_sd(&values[262]);
#if defined(__SSE3__) && defined(__AVX__)
c64_1 = _mm_add_sd(c64_1, _mm_mul_sd(a64_1, _mm256_castpd256_pd128(b64)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c64_1 = _mm_add_sd(c64_1, _mm_mul_sd(a64_1, b64));
#endif
_mm_store_sd(&C[(i*84)+64], c64_1);
#else
C[(i*84)+36] += values[261] * B[(i*84)+64];
C[(i*84)+64] += values[262] * B[(i*84)+64];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b65 = _mm256_broadcast_sd(&B[(i*84)+65]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b65 = _mm_loaddup_pd(&B[(i*84)+65]);
#endif
__m128d c65_0 = _mm_load_sd(&C[(i*84)+37]);
__m128d a65_0 = _mm_load_sd(&values[263]);
#if defined(__SSE3__) && defined(__AVX__)
c65_0 = _mm_add_sd(c65_0, _mm_mul_sd(a65_0, _mm256_castpd256_pd128(b65)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c65_0 = _mm_add_sd(c65_0, _mm_mul_sd(a65_0, b65));
#endif
_mm_store_sd(&C[(i*84)+37], c65_0);
__m128d c65_1 = _mm_load_sd(&C[(i*84)+65]);
__m128d a65_1 = _mm_load_sd(&values[264]);
#if defined(__SSE3__) && defined(__AVX__)
c65_1 = _mm_add_sd(c65_1, _mm_mul_sd(a65_1, _mm256_castpd256_pd128(b65)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c65_1 = _mm_add_sd(c65_1, _mm_mul_sd(a65_1, b65));
#endif
_mm_store_sd(&C[(i*84)+65], c65_1);
#else
C[(i*84)+37] += values[263] * B[(i*84)+65];
C[(i*84)+65] += values[264] * B[(i*84)+65];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b66 = _mm256_broadcast_sd(&B[(i*84)+66]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b66 = _mm_loaddup_pd(&B[(i*84)+66]);
#endif
__m128d c66_0 = _mm_load_sd(&C[(i*84)+38]);
__m128d a66_0 = _mm_load_sd(&values[265]);
#if defined(__SSE3__) && defined(__AVX__)
c66_0 = _mm_add_sd(c66_0, _mm_mul_sd(a66_0, _mm256_castpd256_pd128(b66)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c66_0 = _mm_add_sd(c66_0, _mm_mul_sd(a66_0, b66));
#endif
_mm_store_sd(&C[(i*84)+38], c66_0);
__m128d c66_1 = _mm_load_sd(&C[(i*84)+66]);
__m128d a66_1 = _mm_load_sd(&values[266]);
#if defined(__SSE3__) && defined(__AVX__)
c66_1 = _mm_add_sd(c66_1, _mm_mul_sd(a66_1, _mm256_castpd256_pd128(b66)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c66_1 = _mm_add_sd(c66_1, _mm_mul_sd(a66_1, b66));
#endif
_mm_store_sd(&C[(i*84)+66], c66_1);
#else
C[(i*84)+38] += values[265] * B[(i*84)+66];
C[(i*84)+66] += values[266] * B[(i*84)+66];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b67 = _mm256_broadcast_sd(&B[(i*84)+67]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b67 = _mm_loaddup_pd(&B[(i*84)+67]);
#endif
__m128d c67_0 = _mm_load_sd(&C[(i*84)+39]);
__m128d a67_0 = _mm_load_sd(&values[267]);
#if defined(__SSE3__) && defined(__AVX__)
c67_0 = _mm_add_sd(c67_0, _mm_mul_sd(a67_0, _mm256_castpd256_pd128(b67)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c67_0 = _mm_add_sd(c67_0, _mm_mul_sd(a67_0, b67));
#endif
_mm_store_sd(&C[(i*84)+39], c67_0);
__m128d c67_1 = _mm_load_sd(&C[(i*84)+67]);
__m128d a67_1 = _mm_load_sd(&values[268]);
#if defined(__SSE3__) && defined(__AVX__)
c67_1 = _mm_add_sd(c67_1, _mm_mul_sd(a67_1, _mm256_castpd256_pd128(b67)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c67_1 = _mm_add_sd(c67_1, _mm_mul_sd(a67_1, b67));
#endif
_mm_store_sd(&C[(i*84)+67], c67_1);
#else
C[(i*84)+39] += values[267] * B[(i*84)+67];
C[(i*84)+67] += values[268] * B[(i*84)+67];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b68 = _mm256_broadcast_sd(&B[(i*84)+68]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b68 = _mm_loaddup_pd(&B[(i*84)+68]);
#endif
__m128d c68_0 = _mm_load_sd(&C[(i*84)+40]);
__m128d a68_0 = _mm_load_sd(&values[269]);
#if defined(__SSE3__) && defined(__AVX__)
c68_0 = _mm_add_sd(c68_0, _mm_mul_sd(a68_0, _mm256_castpd256_pd128(b68)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c68_0 = _mm_add_sd(c68_0, _mm_mul_sd(a68_0, b68));
#endif
_mm_store_sd(&C[(i*84)+40], c68_0);
__m128d c68_1 = _mm_load_sd(&C[(i*84)+68]);
__m128d a68_1 = _mm_load_sd(&values[270]);
#if defined(__SSE3__) && defined(__AVX__)
c68_1 = _mm_add_sd(c68_1, _mm_mul_sd(a68_1, _mm256_castpd256_pd128(b68)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c68_1 = _mm_add_sd(c68_1, _mm_mul_sd(a68_1, b68));
#endif
_mm_store_sd(&C[(i*84)+68], c68_1);
#else
C[(i*84)+40] += values[269] * B[(i*84)+68];
C[(i*84)+68] += values[270] * B[(i*84)+68];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b69 = _mm256_broadcast_sd(&B[(i*84)+69]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b69 = _mm_loaddup_pd(&B[(i*84)+69]);
#endif
__m128d c69_0 = _mm_load_sd(&C[(i*84)+20]);
__m128d a69_0 = _mm_load_sd(&values[271]);
#if defined(__SSE3__) && defined(__AVX__)
c69_0 = _mm_add_sd(c69_0, _mm_mul_sd(a69_0, _mm256_castpd256_pd128(b69)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c69_0 = _mm_add_sd(c69_0, _mm_mul_sd(a69_0, b69));
#endif
_mm_store_sd(&C[(i*84)+20], c69_0);
__m128d c69_1 = _mm_load_sd(&C[(i*84)+41]);
__m128d a69_1 = _mm_load_sd(&values[272]);
#if defined(__SSE3__) && defined(__AVX__)
c69_1 = _mm_add_sd(c69_1, _mm_mul_sd(a69_1, _mm256_castpd256_pd128(b69)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c69_1 = _mm_add_sd(c69_1, _mm_mul_sd(a69_1, b69));
#endif
_mm_store_sd(&C[(i*84)+41], c69_1);
__m128d c69_2 = _mm_load_sd(&C[(i*84)+69]);
__m128d a69_2 = _mm_load_sd(&values[273]);
#if defined(__SSE3__) && defined(__AVX__)
c69_2 = _mm_add_sd(c69_2, _mm_mul_sd(a69_2, _mm256_castpd256_pd128(b69)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c69_2 = _mm_add_sd(c69_2, _mm_mul_sd(a69_2, b69));
#endif
_mm_store_sd(&C[(i*84)+69], c69_2);
#else
C[(i*84)+20] += values[271] * B[(i*84)+69];
C[(i*84)+41] += values[272] * B[(i*84)+69];
C[(i*84)+69] += values[273] * B[(i*84)+69];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b70 = _mm256_broadcast_sd(&B[(i*84)+70]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b70 = _mm_loaddup_pd(&B[(i*84)+70]);
#endif
__m128d c70_0 = _mm_load_sd(&C[(i*84)+21]);
__m128d a70_0 = _mm_load_sd(&values[274]);
#if defined(__SSE3__) && defined(__AVX__)
c70_0 = _mm_add_sd(c70_0, _mm_mul_sd(a70_0, _mm256_castpd256_pd128(b70)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c70_0 = _mm_add_sd(c70_0, _mm_mul_sd(a70_0, b70));
#endif
_mm_store_sd(&C[(i*84)+21], c70_0);
__m128d c70_1 = _mm_load_sd(&C[(i*84)+42]);
__m128d a70_1 = _mm_load_sd(&values[275]);
#if defined(__SSE3__) && defined(__AVX__)
c70_1 = _mm_add_sd(c70_1, _mm_mul_sd(a70_1, _mm256_castpd256_pd128(b70)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c70_1 = _mm_add_sd(c70_1, _mm_mul_sd(a70_1, b70));
#endif
_mm_store_sd(&C[(i*84)+42], c70_1);
__m128d c70_2 = _mm_load_sd(&C[(i*84)+70]);
__m128d a70_2 = _mm_load_sd(&values[276]);
#if defined(__SSE3__) && defined(__AVX__)
c70_2 = _mm_add_sd(c70_2, _mm_mul_sd(a70_2, _mm256_castpd256_pd128(b70)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c70_2 = _mm_add_sd(c70_2, _mm_mul_sd(a70_2, b70));
#endif
_mm_store_sd(&C[(i*84)+70], c70_2);
#else
C[(i*84)+21] += values[274] * B[(i*84)+70];
C[(i*84)+42] += values[275] * B[(i*84)+70];
C[(i*84)+70] += values[276] * B[(i*84)+70];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b71 = _mm256_broadcast_sd(&B[(i*84)+71]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b71 = _mm_loaddup_pd(&B[(i*84)+71]);
#endif
__m128d c71_0 = _mm_load_sd(&C[(i*84)+22]);
__m128d a71_0 = _mm_load_sd(&values[277]);
#if defined(__SSE3__) && defined(__AVX__)
c71_0 = _mm_add_sd(c71_0, _mm_mul_sd(a71_0, _mm256_castpd256_pd128(b71)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c71_0 = _mm_add_sd(c71_0, _mm_mul_sd(a71_0, b71));
#endif
_mm_store_sd(&C[(i*84)+22], c71_0);
__m128d c71_1 = _mm_load_sd(&C[(i*84)+43]);
__m128d a71_1 = _mm_load_sd(&values[278]);
#if defined(__SSE3__) && defined(__AVX__)
c71_1 = _mm_add_sd(c71_1, _mm_mul_sd(a71_1, _mm256_castpd256_pd128(b71)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c71_1 = _mm_add_sd(c71_1, _mm_mul_sd(a71_1, b71));
#endif
_mm_store_sd(&C[(i*84)+43], c71_1);
__m128d c71_2 = _mm_load_sd(&C[(i*84)+71]);
__m128d a71_2 = _mm_load_sd(&values[279]);
#if defined(__SSE3__) && defined(__AVX__)
c71_2 = _mm_add_sd(c71_2, _mm_mul_sd(a71_2, _mm256_castpd256_pd128(b71)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c71_2 = _mm_add_sd(c71_2, _mm_mul_sd(a71_2, b71));
#endif
_mm_store_sd(&C[(i*84)+71], c71_2);
#else
C[(i*84)+22] += values[277] * B[(i*84)+71];
C[(i*84)+43] += values[278] * B[(i*84)+71];
C[(i*84)+71] += values[279] * B[(i*84)+71];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b72 = _mm256_broadcast_sd(&B[(i*84)+72]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b72 = _mm_loaddup_pd(&B[(i*84)+72]);
#endif
__m128d c72_0 = _mm_load_sd(&C[(i*84)+23]);
__m128d a72_0 = _mm_load_sd(&values[280]);
#if defined(__SSE3__) && defined(__AVX__)
c72_0 = _mm_add_sd(c72_0, _mm_mul_sd(a72_0, _mm256_castpd256_pd128(b72)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c72_0 = _mm_add_sd(c72_0, _mm_mul_sd(a72_0, b72));
#endif
_mm_store_sd(&C[(i*84)+23], c72_0);
__m128d c72_1 = _mm_load_sd(&C[(i*84)+44]);
__m128d a72_1 = _mm_load_sd(&values[281]);
#if defined(__SSE3__) && defined(__AVX__)
c72_1 = _mm_add_sd(c72_1, _mm_mul_sd(a72_1, _mm256_castpd256_pd128(b72)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c72_1 = _mm_add_sd(c72_1, _mm_mul_sd(a72_1, b72));
#endif
_mm_store_sd(&C[(i*84)+44], c72_1);
__m128d c72_2 = _mm_load_sd(&C[(i*84)+72]);
__m128d a72_2 = _mm_load_sd(&values[282]);
#if defined(__SSE3__) && defined(__AVX__)
c72_2 = _mm_add_sd(c72_2, _mm_mul_sd(a72_2, _mm256_castpd256_pd128(b72)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c72_2 = _mm_add_sd(c72_2, _mm_mul_sd(a72_2, b72));
#endif
_mm_store_sd(&C[(i*84)+72], c72_2);
#else
C[(i*84)+23] += values[280] * B[(i*84)+72];
C[(i*84)+44] += values[281] * B[(i*84)+72];
C[(i*84)+72] += values[282] * B[(i*84)+72];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b73 = _mm256_broadcast_sd(&B[(i*84)+73]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b73 = _mm_loaddup_pd(&B[(i*84)+73]);
#endif
__m128d c73_0 = _mm_load_sd(&C[(i*84)+24]);
__m128d a73_0 = _mm_load_sd(&values[283]);
#if defined(__SSE3__) && defined(__AVX__)
c73_0 = _mm_add_sd(c73_0, _mm_mul_sd(a73_0, _mm256_castpd256_pd128(b73)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c73_0 = _mm_add_sd(c73_0, _mm_mul_sd(a73_0, b73));
#endif
_mm_store_sd(&C[(i*84)+24], c73_0);
__m128d c73_1 = _mm_load_sd(&C[(i*84)+45]);
__m128d a73_1 = _mm_load_sd(&values[284]);
#if defined(__SSE3__) && defined(__AVX__)
c73_1 = _mm_add_sd(c73_1, _mm_mul_sd(a73_1, _mm256_castpd256_pd128(b73)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c73_1 = _mm_add_sd(c73_1, _mm_mul_sd(a73_1, b73));
#endif
_mm_store_sd(&C[(i*84)+45], c73_1);
__m128d c73_2 = _mm_load_sd(&C[(i*84)+73]);
__m128d a73_2 = _mm_load_sd(&values[285]);
#if defined(__SSE3__) && defined(__AVX__)
c73_2 = _mm_add_sd(c73_2, _mm_mul_sd(a73_2, _mm256_castpd256_pd128(b73)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c73_2 = _mm_add_sd(c73_2, _mm_mul_sd(a73_2, b73));
#endif
_mm_store_sd(&C[(i*84)+73], c73_2);
#else
C[(i*84)+24] += values[283] * B[(i*84)+73];
C[(i*84)+45] += values[284] * B[(i*84)+73];
C[(i*84)+73] += values[285] * B[(i*84)+73];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b74 = _mm256_broadcast_sd(&B[(i*84)+74]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b74 = _mm_loaddup_pd(&B[(i*84)+74]);
#endif
__m128d c74_0 = _mm_load_sd(&C[(i*84)+10]);
__m128d a74_0 = _mm_load_sd(&values[286]);
#if defined(__SSE3__) && defined(__AVX__)
c74_0 = _mm_add_sd(c74_0, _mm_mul_sd(a74_0, _mm256_castpd256_pd128(b74)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c74_0 = _mm_add_sd(c74_0, _mm_mul_sd(a74_0, b74));
#endif
_mm_store_sd(&C[(i*84)+10], c74_0);
__m128d c74_1 = _mm_load_sd(&C[(i*84)+25]);
__m128d a74_1 = _mm_load_sd(&values[287]);
#if defined(__SSE3__) && defined(__AVX__)
c74_1 = _mm_add_sd(c74_1, _mm_mul_sd(a74_1, _mm256_castpd256_pd128(b74)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c74_1 = _mm_add_sd(c74_1, _mm_mul_sd(a74_1, b74));
#endif
_mm_store_sd(&C[(i*84)+25], c74_1);
__m128d c74_2 = _mm_load_sd(&C[(i*84)+46]);
__m128d a74_2 = _mm_load_sd(&values[288]);
#if defined(__SSE3__) && defined(__AVX__)
c74_2 = _mm_add_sd(c74_2, _mm_mul_sd(a74_2, _mm256_castpd256_pd128(b74)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c74_2 = _mm_add_sd(c74_2, _mm_mul_sd(a74_2, b74));
#endif
_mm_store_sd(&C[(i*84)+46], c74_2);
__m128d c74_3 = _mm_load_sd(&C[(i*84)+74]);
__m128d a74_3 = _mm_load_sd(&values[289]);
#if defined(__SSE3__) && defined(__AVX__)
c74_3 = _mm_add_sd(c74_3, _mm_mul_sd(a74_3, _mm256_castpd256_pd128(b74)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c74_3 = _mm_add_sd(c74_3, _mm_mul_sd(a74_3, b74));
#endif
_mm_store_sd(&C[(i*84)+74], c74_3);
#else
C[(i*84)+10] += values[286] * B[(i*84)+74];
C[(i*84)+25] += values[287] * B[(i*84)+74];
C[(i*84)+46] += values[288] * B[(i*84)+74];
C[(i*84)+74] += values[289] * B[(i*84)+74];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b75 = _mm256_broadcast_sd(&B[(i*84)+75]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b75 = _mm_loaddup_pd(&B[(i*84)+75]);
#endif
__m128d c75_0 = _mm_load_sd(&C[(i*84)+11]);
__m128d a75_0 = _mm_load_sd(&values[290]);
#if defined(__SSE3__) && defined(__AVX__)
c75_0 = _mm_add_sd(c75_0, _mm_mul_sd(a75_0, _mm256_castpd256_pd128(b75)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c75_0 = _mm_add_sd(c75_0, _mm_mul_sd(a75_0, b75));
#endif
_mm_store_sd(&C[(i*84)+11], c75_0);
__m128d c75_1 = _mm_load_sd(&C[(i*84)+26]);
__m128d a75_1 = _mm_load_sd(&values[291]);
#if defined(__SSE3__) && defined(__AVX__)
c75_1 = _mm_add_sd(c75_1, _mm_mul_sd(a75_1, _mm256_castpd256_pd128(b75)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c75_1 = _mm_add_sd(c75_1, _mm_mul_sd(a75_1, b75));
#endif
_mm_store_sd(&C[(i*84)+26], c75_1);
__m128d c75_2 = _mm_load_sd(&C[(i*84)+47]);
__m128d a75_2 = _mm_load_sd(&values[292]);
#if defined(__SSE3__) && defined(__AVX__)
c75_2 = _mm_add_sd(c75_2, _mm_mul_sd(a75_2, _mm256_castpd256_pd128(b75)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c75_2 = _mm_add_sd(c75_2, _mm_mul_sd(a75_2, b75));
#endif
_mm_store_sd(&C[(i*84)+47], c75_2);
__m128d c75_3 = _mm_load_sd(&C[(i*84)+75]);
__m128d a75_3 = _mm_load_sd(&values[293]);
#if defined(__SSE3__) && defined(__AVX__)
c75_3 = _mm_add_sd(c75_3, _mm_mul_sd(a75_3, _mm256_castpd256_pd128(b75)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c75_3 = _mm_add_sd(c75_3, _mm_mul_sd(a75_3, b75));
#endif
_mm_store_sd(&C[(i*84)+75], c75_3);
#else
C[(i*84)+11] += values[290] * B[(i*84)+75];
C[(i*84)+26] += values[291] * B[(i*84)+75];
C[(i*84)+47] += values[292] * B[(i*84)+75];
C[(i*84)+75] += values[293] * B[(i*84)+75];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b76 = _mm256_broadcast_sd(&B[(i*84)+76]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b76 = _mm_loaddup_pd(&B[(i*84)+76]);
#endif
__m128d c76_0 = _mm_load_sd(&C[(i*84)+12]);
__m128d a76_0 = _mm_load_sd(&values[294]);
#if defined(__SSE3__) && defined(__AVX__)
c76_0 = _mm_add_sd(c76_0, _mm_mul_sd(a76_0, _mm256_castpd256_pd128(b76)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c76_0 = _mm_add_sd(c76_0, _mm_mul_sd(a76_0, b76));
#endif
_mm_store_sd(&C[(i*84)+12], c76_0);
__m128d c76_1 = _mm_load_sd(&C[(i*84)+27]);
__m128d a76_1 = _mm_load_sd(&values[295]);
#if defined(__SSE3__) && defined(__AVX__)
c76_1 = _mm_add_sd(c76_1, _mm_mul_sd(a76_1, _mm256_castpd256_pd128(b76)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c76_1 = _mm_add_sd(c76_1, _mm_mul_sd(a76_1, b76));
#endif
_mm_store_sd(&C[(i*84)+27], c76_1);
__m128d c76_2 = _mm_load_sd(&C[(i*84)+48]);
__m128d a76_2 = _mm_load_sd(&values[296]);
#if defined(__SSE3__) && defined(__AVX__)
c76_2 = _mm_add_sd(c76_2, _mm_mul_sd(a76_2, _mm256_castpd256_pd128(b76)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c76_2 = _mm_add_sd(c76_2, _mm_mul_sd(a76_2, b76));
#endif
_mm_store_sd(&C[(i*84)+48], c76_2);
__m128d c76_3 = _mm_load_sd(&C[(i*84)+76]);
__m128d a76_3 = _mm_load_sd(&values[297]);
#if defined(__SSE3__) && defined(__AVX__)
c76_3 = _mm_add_sd(c76_3, _mm_mul_sd(a76_3, _mm256_castpd256_pd128(b76)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c76_3 = _mm_add_sd(c76_3, _mm_mul_sd(a76_3, b76));
#endif
_mm_store_sd(&C[(i*84)+76], c76_3);
#else
C[(i*84)+12] += values[294] * B[(i*84)+76];
C[(i*84)+27] += values[295] * B[(i*84)+76];
C[(i*84)+48] += values[296] * B[(i*84)+76];
C[(i*84)+76] += values[297] * B[(i*84)+76];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b77 = _mm256_broadcast_sd(&B[(i*84)+77]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b77 = _mm_loaddup_pd(&B[(i*84)+77]);
#endif
__m128d c77_0 = _mm_load_sd(&C[(i*84)+13]);
__m128d a77_0 = _mm_load_sd(&values[298]);
#if defined(__SSE3__) && defined(__AVX__)
c77_0 = _mm_add_sd(c77_0, _mm_mul_sd(a77_0, _mm256_castpd256_pd128(b77)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c77_0 = _mm_add_sd(c77_0, _mm_mul_sd(a77_0, b77));
#endif
_mm_store_sd(&C[(i*84)+13], c77_0);
__m128d c77_1 = _mm_load_sd(&C[(i*84)+28]);
__m128d a77_1 = _mm_load_sd(&values[299]);
#if defined(__SSE3__) && defined(__AVX__)
c77_1 = _mm_add_sd(c77_1, _mm_mul_sd(a77_1, _mm256_castpd256_pd128(b77)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c77_1 = _mm_add_sd(c77_1, _mm_mul_sd(a77_1, b77));
#endif
_mm_store_sd(&C[(i*84)+28], c77_1);
__m128d c77_2 = _mm_load_sd(&C[(i*84)+49]);
__m128d a77_2 = _mm_load_sd(&values[300]);
#if defined(__SSE3__) && defined(__AVX__)
c77_2 = _mm_add_sd(c77_2, _mm_mul_sd(a77_2, _mm256_castpd256_pd128(b77)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c77_2 = _mm_add_sd(c77_2, _mm_mul_sd(a77_2, b77));
#endif
_mm_store_sd(&C[(i*84)+49], c77_2);
__m128d c77_3 = _mm_load_sd(&C[(i*84)+77]);
__m128d a77_3 = _mm_load_sd(&values[301]);
#if defined(__SSE3__) && defined(__AVX__)
c77_3 = _mm_add_sd(c77_3, _mm_mul_sd(a77_3, _mm256_castpd256_pd128(b77)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c77_3 = _mm_add_sd(c77_3, _mm_mul_sd(a77_3, b77));
#endif
_mm_store_sd(&C[(i*84)+77], c77_3);
#else
C[(i*84)+13] += values[298] * B[(i*84)+77];
C[(i*84)+28] += values[299] * B[(i*84)+77];
C[(i*84)+49] += values[300] * B[(i*84)+77];
C[(i*84)+77] += values[301] * B[(i*84)+77];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b78 = _mm256_broadcast_sd(&B[(i*84)+78]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b78 = _mm_loaddup_pd(&B[(i*84)+78]);
#endif
__m128d c78_0 = _mm_load_sd(&C[(i*84)+4]);
__m128d a78_0 = _mm_load_sd(&values[302]);
#if defined(__SSE3__) && defined(__AVX__)
c78_0 = _mm_add_sd(c78_0, _mm_mul_sd(a78_0, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c78_0 = _mm_add_sd(c78_0, _mm_mul_sd(a78_0, b78));
#endif
_mm_store_sd(&C[(i*84)+4], c78_0);
__m128d c78_1 = _mm_load_sd(&C[(i*84)+14]);
__m128d a78_1 = _mm_load_sd(&values[303]);
#if defined(__SSE3__) && defined(__AVX__)
c78_1 = _mm_add_sd(c78_1, _mm_mul_sd(a78_1, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c78_1 = _mm_add_sd(c78_1, _mm_mul_sd(a78_1, b78));
#endif
_mm_store_sd(&C[(i*84)+14], c78_1);
__m128d c78_2 = _mm_load_sd(&C[(i*84)+29]);
__m128d a78_2 = _mm_load_sd(&values[304]);
#if defined(__SSE3__) && defined(__AVX__)
c78_2 = _mm_add_sd(c78_2, _mm_mul_sd(a78_2, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c78_2 = _mm_add_sd(c78_2, _mm_mul_sd(a78_2, b78));
#endif
_mm_store_sd(&C[(i*84)+29], c78_2);
__m128d c78_3 = _mm_load_sd(&C[(i*84)+50]);
__m128d a78_3 = _mm_load_sd(&values[305]);
#if defined(__SSE3__) && defined(__AVX__)
c78_3 = _mm_add_sd(c78_3, _mm_mul_sd(a78_3, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c78_3 = _mm_add_sd(c78_3, _mm_mul_sd(a78_3, b78));
#endif
_mm_store_sd(&C[(i*84)+50], c78_3);
__m128d c78_4 = _mm_load_sd(&C[(i*84)+78]);
__m128d a78_4 = _mm_load_sd(&values[306]);
#if defined(__SSE3__) && defined(__AVX__)
c78_4 = _mm_add_sd(c78_4, _mm_mul_sd(a78_4, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c78_4 = _mm_add_sd(c78_4, _mm_mul_sd(a78_4, b78));
#endif
_mm_store_sd(&C[(i*84)+78], c78_4);
#else
C[(i*84)+4] += values[302] * B[(i*84)+78];
C[(i*84)+14] += values[303] * B[(i*84)+78];
C[(i*84)+29] += values[304] * B[(i*84)+78];
C[(i*84)+50] += values[305] * B[(i*84)+78];
C[(i*84)+78] += values[306] * B[(i*84)+78];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b79 = _mm256_broadcast_sd(&B[(i*84)+79]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b79 = _mm_loaddup_pd(&B[(i*84)+79]);
#endif
__m128d c79_0 = _mm_load_sd(&C[(i*84)+5]);
__m128d a79_0 = _mm_load_sd(&values[307]);
#if defined(__SSE3__) && defined(__AVX__)
c79_0 = _mm_add_sd(c79_0, _mm_mul_sd(a79_0, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c79_0 = _mm_add_sd(c79_0, _mm_mul_sd(a79_0, b79));
#endif
_mm_store_sd(&C[(i*84)+5], c79_0);
__m128d c79_1 = _mm_load_sd(&C[(i*84)+15]);
__m128d a79_1 = _mm_load_sd(&values[308]);
#if defined(__SSE3__) && defined(__AVX__)
c79_1 = _mm_add_sd(c79_1, _mm_mul_sd(a79_1, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c79_1 = _mm_add_sd(c79_1, _mm_mul_sd(a79_1, b79));
#endif
_mm_store_sd(&C[(i*84)+15], c79_1);
__m128d c79_2 = _mm_load_sd(&C[(i*84)+30]);
__m128d a79_2 = _mm_load_sd(&values[309]);
#if defined(__SSE3__) && defined(__AVX__)
c79_2 = _mm_add_sd(c79_2, _mm_mul_sd(a79_2, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c79_2 = _mm_add_sd(c79_2, _mm_mul_sd(a79_2, b79));
#endif
_mm_store_sd(&C[(i*84)+30], c79_2);
__m128d c79_3 = _mm_load_sd(&C[(i*84)+51]);
__m128d a79_3 = _mm_load_sd(&values[310]);
#if defined(__SSE3__) && defined(__AVX__)
c79_3 = _mm_add_sd(c79_3, _mm_mul_sd(a79_3, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c79_3 = _mm_add_sd(c79_3, _mm_mul_sd(a79_3, b79));
#endif
_mm_store_sd(&C[(i*84)+51], c79_3);
__m128d c79_4 = _mm_load_sd(&C[(i*84)+79]);
__m128d a79_4 = _mm_load_sd(&values[311]);
#if defined(__SSE3__) && defined(__AVX__)
c79_4 = _mm_add_sd(c79_4, _mm_mul_sd(a79_4, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c79_4 = _mm_add_sd(c79_4, _mm_mul_sd(a79_4, b79));
#endif
_mm_store_sd(&C[(i*84)+79], c79_4);
#else
C[(i*84)+5] += values[307] * B[(i*84)+79];
C[(i*84)+15] += values[308] * B[(i*84)+79];
C[(i*84)+30] += values[309] * B[(i*84)+79];
C[(i*84)+51] += values[310] * B[(i*84)+79];
C[(i*84)+79] += values[311] * B[(i*84)+79];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b80 = _mm256_broadcast_sd(&B[(i*84)+80]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b80 = _mm_loaddup_pd(&B[(i*84)+80]);
#endif
__m128d c80_0 = _mm_load_sd(&C[(i*84)+6]);
__m128d a80_0 = _mm_load_sd(&values[312]);
#if defined(__SSE3__) && defined(__AVX__)
c80_0 = _mm_add_sd(c80_0, _mm_mul_sd(a80_0, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c80_0 = _mm_add_sd(c80_0, _mm_mul_sd(a80_0, b80));
#endif
_mm_store_sd(&C[(i*84)+6], c80_0);
__m128d c80_1 = _mm_load_sd(&C[(i*84)+16]);
__m128d a80_1 = _mm_load_sd(&values[313]);
#if defined(__SSE3__) && defined(__AVX__)
c80_1 = _mm_add_sd(c80_1, _mm_mul_sd(a80_1, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c80_1 = _mm_add_sd(c80_1, _mm_mul_sd(a80_1, b80));
#endif
_mm_store_sd(&C[(i*84)+16], c80_1);
__m128d c80_2 = _mm_load_sd(&C[(i*84)+31]);
__m128d a80_2 = _mm_load_sd(&values[314]);
#if defined(__SSE3__) && defined(__AVX__)
c80_2 = _mm_add_sd(c80_2, _mm_mul_sd(a80_2, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c80_2 = _mm_add_sd(c80_2, _mm_mul_sd(a80_2, b80));
#endif
_mm_store_sd(&C[(i*84)+31], c80_2);
__m128d c80_3 = _mm_load_sd(&C[(i*84)+52]);
__m128d a80_3 = _mm_load_sd(&values[315]);
#if defined(__SSE3__) && defined(__AVX__)
c80_3 = _mm_add_sd(c80_3, _mm_mul_sd(a80_3, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c80_3 = _mm_add_sd(c80_3, _mm_mul_sd(a80_3, b80));
#endif
_mm_store_sd(&C[(i*84)+52], c80_3);
__m128d c80_4 = _mm_load_sd(&C[(i*84)+80]);
__m128d a80_4 = _mm_load_sd(&values[316]);
#if defined(__SSE3__) && defined(__AVX__)
c80_4 = _mm_add_sd(c80_4, _mm_mul_sd(a80_4, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c80_4 = _mm_add_sd(c80_4, _mm_mul_sd(a80_4, b80));
#endif
_mm_store_sd(&C[(i*84)+80], c80_4);
#else
C[(i*84)+6] += values[312] * B[(i*84)+80];
C[(i*84)+16] += values[313] * B[(i*84)+80];
C[(i*84)+31] += values[314] * B[(i*84)+80];
C[(i*84)+52] += values[315] * B[(i*84)+80];
C[(i*84)+80] += values[316] * B[(i*84)+80];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b81 = _mm256_broadcast_sd(&B[(i*84)+81]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b81 = _mm_loaddup_pd(&B[(i*84)+81]);
#endif
__m128d c81_0 = _mm_load_sd(&C[(i*84)+1]);
__m128d a81_0 = _mm_load_sd(&values[317]);
#if defined(__SSE3__) && defined(__AVX__)
c81_0 = _mm_add_sd(c81_0, _mm_mul_sd(a81_0, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c81_0 = _mm_add_sd(c81_0, _mm_mul_sd(a81_0, b81));
#endif
_mm_store_sd(&C[(i*84)+1], c81_0);
__m128d c81_1 = _mm_load_sd(&C[(i*84)+7]);
__m128d a81_1 = _mm_load_sd(&values[318]);
#if defined(__SSE3__) && defined(__AVX__)
c81_1 = _mm_add_sd(c81_1, _mm_mul_sd(a81_1, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c81_1 = _mm_add_sd(c81_1, _mm_mul_sd(a81_1, b81));
#endif
_mm_store_sd(&C[(i*84)+7], c81_1);
__m128d c81_2 = _mm_load_sd(&C[(i*84)+17]);
__m128d a81_2 = _mm_load_sd(&values[319]);
#if defined(__SSE3__) && defined(__AVX__)
c81_2 = _mm_add_sd(c81_2, _mm_mul_sd(a81_2, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c81_2 = _mm_add_sd(c81_2, _mm_mul_sd(a81_2, b81));
#endif
_mm_store_sd(&C[(i*84)+17], c81_2);
__m128d c81_3 = _mm_load_sd(&C[(i*84)+32]);
__m128d a81_3 = _mm_load_sd(&values[320]);
#if defined(__SSE3__) && defined(__AVX__)
c81_3 = _mm_add_sd(c81_3, _mm_mul_sd(a81_3, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c81_3 = _mm_add_sd(c81_3, _mm_mul_sd(a81_3, b81));
#endif
_mm_store_sd(&C[(i*84)+32], c81_3);
__m128d c81_4 = _mm_load_sd(&C[(i*84)+53]);
__m128d a81_4 = _mm_load_sd(&values[321]);
#if defined(__SSE3__) && defined(__AVX__)
c81_4 = _mm_add_sd(c81_4, _mm_mul_sd(a81_4, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c81_4 = _mm_add_sd(c81_4, _mm_mul_sd(a81_4, b81));
#endif
_mm_store_sd(&C[(i*84)+53], c81_4);
__m128d c81_5 = _mm_load_sd(&C[(i*84)+81]);
__m128d a81_5 = _mm_load_sd(&values[322]);
#if defined(__SSE3__) && defined(__AVX__)
c81_5 = _mm_add_sd(c81_5, _mm_mul_sd(a81_5, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c81_5 = _mm_add_sd(c81_5, _mm_mul_sd(a81_5, b81));
#endif
_mm_store_sd(&C[(i*84)+81], c81_5);
#else
C[(i*84)+1] += values[317] * B[(i*84)+81];
C[(i*84)+7] += values[318] * B[(i*84)+81];
C[(i*84)+17] += values[319] * B[(i*84)+81];
C[(i*84)+32] += values[320] * B[(i*84)+81];
C[(i*84)+53] += values[321] * B[(i*84)+81];
C[(i*84)+81] += values[322] * B[(i*84)+81];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b82 = _mm256_broadcast_sd(&B[(i*84)+82]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b82 = _mm_loaddup_pd(&B[(i*84)+82]);
#endif
__m128d c82_0 = _mm_load_sd(&C[(i*84)+2]);
__m128d a82_0 = _mm_load_sd(&values[323]);
#if defined(__SSE3__) && defined(__AVX__)
c82_0 = _mm_add_sd(c82_0, _mm_mul_sd(a82_0, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c82_0 = _mm_add_sd(c82_0, _mm_mul_sd(a82_0, b82));
#endif
_mm_store_sd(&C[(i*84)+2], c82_0);
__m128d c82_1 = _mm_load_sd(&C[(i*84)+8]);
__m128d a82_1 = _mm_load_sd(&values[324]);
#if defined(__SSE3__) && defined(__AVX__)
c82_1 = _mm_add_sd(c82_1, _mm_mul_sd(a82_1, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c82_1 = _mm_add_sd(c82_1, _mm_mul_sd(a82_1, b82));
#endif
_mm_store_sd(&C[(i*84)+8], c82_1);
__m128d c82_2 = _mm_load_sd(&C[(i*84)+18]);
__m128d a82_2 = _mm_load_sd(&values[325]);
#if defined(__SSE3__) && defined(__AVX__)
c82_2 = _mm_add_sd(c82_2, _mm_mul_sd(a82_2, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c82_2 = _mm_add_sd(c82_2, _mm_mul_sd(a82_2, b82));
#endif
_mm_store_sd(&C[(i*84)+18], c82_2);
__m128d c82_3 = _mm_load_sd(&C[(i*84)+33]);
__m128d a82_3 = _mm_load_sd(&values[326]);
#if defined(__SSE3__) && defined(__AVX__)
c82_3 = _mm_add_sd(c82_3, _mm_mul_sd(a82_3, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c82_3 = _mm_add_sd(c82_3, _mm_mul_sd(a82_3, b82));
#endif
_mm_store_sd(&C[(i*84)+33], c82_3);
__m128d c82_4 = _mm_load_sd(&C[(i*84)+54]);
__m128d a82_4 = _mm_load_sd(&values[327]);
#if defined(__SSE3__) && defined(__AVX__)
c82_4 = _mm_add_sd(c82_4, _mm_mul_sd(a82_4, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c82_4 = _mm_add_sd(c82_4, _mm_mul_sd(a82_4, b82));
#endif
_mm_store_sd(&C[(i*84)+54], c82_4);
__m128d c82_5 = _mm_load_sd(&C[(i*84)+82]);
__m128d a82_5 = _mm_load_sd(&values[328]);
#if defined(__SSE3__) && defined(__AVX__)
c82_5 = _mm_add_sd(c82_5, _mm_mul_sd(a82_5, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c82_5 = _mm_add_sd(c82_5, _mm_mul_sd(a82_5, b82));
#endif
_mm_store_sd(&C[(i*84)+82], c82_5);
#else
C[(i*84)+2] += values[323] * B[(i*84)+82];
C[(i*84)+8] += values[324] * B[(i*84)+82];
C[(i*84)+18] += values[325] * B[(i*84)+82];
C[(i*84)+33] += values[326] * B[(i*84)+82];
C[(i*84)+54] += values[327] * B[(i*84)+82];
C[(i*84)+82] += values[328] * B[(i*84)+82];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b83 = _mm256_broadcast_sd(&B[(i*84)+83]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b83 = _mm_loaddup_pd(&B[(i*84)+83]);
#endif
__m128d c83_0 = _mm_load_sd(&C[(i*84)+0]);
__m128d a83_0 = _mm_load_sd(&values[329]);
#if defined(__SSE3__) && defined(__AVX__)
c83_0 = _mm_add_sd(c83_0, _mm_mul_sd(a83_0, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c83_0 = _mm_add_sd(c83_0, _mm_mul_sd(a83_0, b83));
#endif
_mm_store_sd(&C[(i*84)+0], c83_0);
__m128d c83_1 = _mm_load_sd(&C[(i*84)+3]);
__m128d a83_1 = _mm_load_sd(&values[330]);
#if defined(__SSE3__) && defined(__AVX__)
c83_1 = _mm_add_sd(c83_1, _mm_mul_sd(a83_1, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c83_1 = _mm_add_sd(c83_1, _mm_mul_sd(a83_1, b83));
#endif
_mm_store_sd(&C[(i*84)+3], c83_1);
__m128d c83_2 = _mm_load_sd(&C[(i*84)+9]);
__m128d a83_2 = _mm_load_sd(&values[331]);
#if defined(__SSE3__) && defined(__AVX__)
c83_2 = _mm_add_sd(c83_2, _mm_mul_sd(a83_2, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c83_2 = _mm_add_sd(c83_2, _mm_mul_sd(a83_2, b83));
#endif
_mm_store_sd(&C[(i*84)+9], c83_2);
__m128d c83_3 = _mm_load_sd(&C[(i*84)+19]);
__m128d a83_3 = _mm_load_sd(&values[332]);
#if defined(__SSE3__) && defined(__AVX__)
c83_3 = _mm_add_sd(c83_3, _mm_mul_sd(a83_3, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c83_3 = _mm_add_sd(c83_3, _mm_mul_sd(a83_3, b83));
#endif
_mm_store_sd(&C[(i*84)+19], c83_3);
__m128d c83_4 = _mm_load_sd(&C[(i*84)+34]);
__m128d a83_4 = _mm_load_sd(&values[333]);
#if defined(__SSE3__) && defined(__AVX__)
c83_4 = _mm_add_sd(c83_4, _mm_mul_sd(a83_4, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c83_4 = _mm_add_sd(c83_4, _mm_mul_sd(a83_4, b83));
#endif
_mm_store_sd(&C[(i*84)+34], c83_4);
__m128d c83_5 = _mm_load_sd(&C[(i*84)+55]);
__m128d a83_5 = _mm_load_sd(&values[334]);
#if defined(__SSE3__) && defined(__AVX__)
c83_5 = _mm_add_sd(c83_5, _mm_mul_sd(a83_5, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c83_5 = _mm_add_sd(c83_5, _mm_mul_sd(a83_5, b83));
#endif
_mm_store_sd(&C[(i*84)+55], c83_5);
__m128d c83_6 = _mm_load_sd(&C[(i*84)+83]);
__m128d a83_6 = _mm_load_sd(&values[335]);
#if defined(__SSE3__) && defined(__AVX__)
c83_6 = _mm_add_sd(c83_6, _mm_mul_sd(a83_6, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c83_6 = _mm_add_sd(c83_6, _mm_mul_sd(a83_6, b83));
#endif
_mm_store_sd(&C[(i*84)+83], c83_6);
#else
C[(i*84)+0] += values[329] * B[(i*84)+83];
C[(i*84)+3] += values[330] * B[(i*84)+83];
C[(i*84)+9] += values[331] * B[(i*84)+83];
C[(i*84)+19] += values[332] * B[(i*84)+83];
C[(i*84)+34] += values[333] * B[(i*84)+83];
C[(i*84)+55] += values[334] * B[(i*84)+83];
C[(i*84)+83] += values[335] * B[(i*84)+83];
#endif

}

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 6048;
#endif

}

void dsparse_fP113DivM_m84_n9_k84_ldAna7_ldB84_ldC84_beta0_pfsigonly(const double* values, const double* B, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma nounroll_and_jam
for (unsigned int i = 0; i < 9; i++)
{
  #pragma simd
  for (unsigned int m = 0; m < 84; m++) {
    C[(i*84)+m] = 0.0;
  }
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b0 = _mm256_broadcast_sd(&B[(i*84)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b0 = _mm_loaddup_pd(&B[(i*84)+0]);
#endif
__m128d c0_0 = _mm_load_sd(&C[(i*84)+0]);
__m128d a0_0 = _mm_load_sd(&values[0]);
#if defined(__SSE3__) && defined(__AVX__)
c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, b0));
#endif
_mm_store_sd(&C[(i*84)+0], c0_0);
__m128d c0_1 = _mm_load_sd(&C[(i*84)+3]);
__m128d a0_1 = _mm_load_sd(&values[1]);
#if defined(__SSE3__) && defined(__AVX__)
c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, b0));
#endif
_mm_store_sd(&C[(i*84)+3], c0_1);
__m128d c0_2 = _mm_load_sd(&C[(i*84)+9]);
__m128d a0_2 = _mm_load_sd(&values[2]);
#if defined(__SSE3__) && defined(__AVX__)
c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, b0));
#endif
_mm_store_sd(&C[(i*84)+9], c0_2);
__m128d c0_3 = _mm_load_sd(&C[(i*84)+19]);
__m128d a0_3 = _mm_load_sd(&values[3]);
#if defined(__SSE3__) && defined(__AVX__)
c0_3 = _mm_add_sd(c0_3, _mm_mul_sd(a0_3, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_3 = _mm_add_sd(c0_3, _mm_mul_sd(a0_3, b0));
#endif
_mm_store_sd(&C[(i*84)+19], c0_3);
__m128d c0_4 = _mm_load_sd(&C[(i*84)+34]);
__m128d a0_4 = _mm_load_sd(&values[4]);
#if defined(__SSE3__) && defined(__AVX__)
c0_4 = _mm_add_sd(c0_4, _mm_mul_sd(a0_4, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_4 = _mm_add_sd(c0_4, _mm_mul_sd(a0_4, b0));
#endif
_mm_store_sd(&C[(i*84)+34], c0_4);
__m128d c0_5 = _mm_load_sd(&C[(i*84)+55]);
__m128d a0_5 = _mm_load_sd(&values[5]);
#if defined(__SSE3__) && defined(__AVX__)
c0_5 = _mm_add_sd(c0_5, _mm_mul_sd(a0_5, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_5 = _mm_add_sd(c0_5, _mm_mul_sd(a0_5, b0));
#endif
_mm_store_sd(&C[(i*84)+55], c0_5);
__m128d c0_6 = _mm_load_sd(&C[(i*84)+83]);
__m128d a0_6 = _mm_load_sd(&values[6]);
#if defined(__SSE3__) && defined(__AVX__)
c0_6 = _mm_add_sd(c0_6, _mm_mul_sd(a0_6, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c0_6 = _mm_add_sd(c0_6, _mm_mul_sd(a0_6, b0));
#endif
_mm_store_sd(&C[(i*84)+83], c0_6);
#else
C[(i*84)+0] += values[0] * B[(i*84)+0];
C[(i*84)+3] += values[1] * B[(i*84)+0];
C[(i*84)+9] += values[2] * B[(i*84)+0];
C[(i*84)+19] += values[3] * B[(i*84)+0];
C[(i*84)+34] += values[4] * B[(i*84)+0];
C[(i*84)+55] += values[5] * B[(i*84)+0];
C[(i*84)+83] += values[6] * B[(i*84)+0];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b1 = _mm256_broadcast_sd(&B[(i*84)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b1 = _mm_loaddup_pd(&B[(i*84)+1]);
#endif
__m128d c1_0 = _mm_load_sd(&C[(i*84)+1]);
__m128d a1_0 = _mm_load_sd(&values[7]);
#if defined(__SSE3__) && defined(__AVX__)
c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, b1));
#endif
_mm_store_sd(&C[(i*84)+1], c1_0);
__m128d c1_1 = _mm_load_sd(&C[(i*84)+7]);
__m128d a1_1 = _mm_load_sd(&values[8]);
#if defined(__SSE3__) && defined(__AVX__)
c1_1 = _mm_add_sd(c1_1, _mm_mul_sd(a1_1, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_1 = _mm_add_sd(c1_1, _mm_mul_sd(a1_1, b1));
#endif
_mm_store_sd(&C[(i*84)+7], c1_1);
__m128d c1_2 = _mm_load_sd(&C[(i*84)+17]);
__m128d a1_2 = _mm_load_sd(&values[9]);
#if defined(__SSE3__) && defined(__AVX__)
c1_2 = _mm_add_sd(c1_2, _mm_mul_sd(a1_2, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_2 = _mm_add_sd(c1_2, _mm_mul_sd(a1_2, b1));
#endif
_mm_store_sd(&C[(i*84)+17], c1_2);
__m128d c1_3 = _mm_load_sd(&C[(i*84)+32]);
__m128d a1_3 = _mm_load_sd(&values[10]);
#if defined(__SSE3__) && defined(__AVX__)
c1_3 = _mm_add_sd(c1_3, _mm_mul_sd(a1_3, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_3 = _mm_add_sd(c1_3, _mm_mul_sd(a1_3, b1));
#endif
_mm_store_sd(&C[(i*84)+32], c1_3);
__m128d c1_4 = _mm_load_sd(&C[(i*84)+53]);
__m128d a1_4 = _mm_load_sd(&values[11]);
#if defined(__SSE3__) && defined(__AVX__)
c1_4 = _mm_add_sd(c1_4, _mm_mul_sd(a1_4, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_4 = _mm_add_sd(c1_4, _mm_mul_sd(a1_4, b1));
#endif
_mm_store_sd(&C[(i*84)+53], c1_4);
__m128d c1_5 = _mm_load_sd(&C[(i*84)+81]);
__m128d a1_5 = _mm_load_sd(&values[12]);
#if defined(__SSE3__) && defined(__AVX__)
c1_5 = _mm_add_sd(c1_5, _mm_mul_sd(a1_5, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c1_5 = _mm_add_sd(c1_5, _mm_mul_sd(a1_5, b1));
#endif
_mm_store_sd(&C[(i*84)+81], c1_5);
#else
C[(i*84)+1] += values[7] * B[(i*84)+1];
C[(i*84)+7] += values[8] * B[(i*84)+1];
C[(i*84)+17] += values[9] * B[(i*84)+1];
C[(i*84)+32] += values[10] * B[(i*84)+1];
C[(i*84)+53] += values[11] * B[(i*84)+1];
C[(i*84)+81] += values[12] * B[(i*84)+1];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b2 = _mm256_broadcast_sd(&B[(i*84)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b2 = _mm_loaddup_pd(&B[(i*84)+2]);
#endif
__m128d c2_0 = _mm_load_sd(&C[(i*84)+2]);
__m128d a2_0 = _mm_load_sd(&values[13]);
#if defined(__SSE3__) && defined(__AVX__)
c2_0 = _mm_add_sd(c2_0, _mm_mul_sd(a2_0, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_0 = _mm_add_sd(c2_0, _mm_mul_sd(a2_0, b2));
#endif
_mm_store_sd(&C[(i*84)+2], c2_0);
__m128d c2_1 = _mm_load_sd(&C[(i*84)+8]);
__m128d a2_1 = _mm_load_sd(&values[14]);
#if defined(__SSE3__) && defined(__AVX__)
c2_1 = _mm_add_sd(c2_1, _mm_mul_sd(a2_1, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_1 = _mm_add_sd(c2_1, _mm_mul_sd(a2_1, b2));
#endif
_mm_store_sd(&C[(i*84)+8], c2_1);
__m128d c2_2 = _mm_load_sd(&C[(i*84)+18]);
__m128d a2_2 = _mm_load_sd(&values[15]);
#if defined(__SSE3__) && defined(__AVX__)
c2_2 = _mm_add_sd(c2_2, _mm_mul_sd(a2_2, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_2 = _mm_add_sd(c2_2, _mm_mul_sd(a2_2, b2));
#endif
_mm_store_sd(&C[(i*84)+18], c2_2);
__m128d c2_3 = _mm_load_sd(&C[(i*84)+33]);
__m128d a2_3 = _mm_load_sd(&values[16]);
#if defined(__SSE3__) && defined(__AVX__)
c2_3 = _mm_add_sd(c2_3, _mm_mul_sd(a2_3, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_3 = _mm_add_sd(c2_3, _mm_mul_sd(a2_3, b2));
#endif
_mm_store_sd(&C[(i*84)+33], c2_3);
__m128d c2_4 = _mm_load_sd(&C[(i*84)+54]);
__m128d a2_4 = _mm_load_sd(&values[17]);
#if defined(__SSE3__) && defined(__AVX__)
c2_4 = _mm_add_sd(c2_4, _mm_mul_sd(a2_4, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_4 = _mm_add_sd(c2_4, _mm_mul_sd(a2_4, b2));
#endif
_mm_store_sd(&C[(i*84)+54], c2_4);
__m128d c2_5 = _mm_load_sd(&C[(i*84)+82]);
__m128d a2_5 = _mm_load_sd(&values[18]);
#if defined(__SSE3__) && defined(__AVX__)
c2_5 = _mm_add_sd(c2_5, _mm_mul_sd(a2_5, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c2_5 = _mm_add_sd(c2_5, _mm_mul_sd(a2_5, b2));
#endif
_mm_store_sd(&C[(i*84)+82], c2_5);
#else
C[(i*84)+2] += values[13] * B[(i*84)+2];
C[(i*84)+8] += values[14] * B[(i*84)+2];
C[(i*84)+18] += values[15] * B[(i*84)+2];
C[(i*84)+33] += values[16] * B[(i*84)+2];
C[(i*84)+54] += values[17] * B[(i*84)+2];
C[(i*84)+82] += values[18] * B[(i*84)+2];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b3 = _mm256_broadcast_sd(&B[(i*84)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b3 = _mm_loaddup_pd(&B[(i*84)+3]);
#endif
__m128d c3_0 = _mm_load_sd(&C[(i*84)+0]);
__m128d a3_0 = _mm_load_sd(&values[19]);
#if defined(__SSE3__) && defined(__AVX__)
c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, b3));
#endif
_mm_store_sd(&C[(i*84)+0], c3_0);
__m128d c3_1 = _mm_load_sd(&C[(i*84)+3]);
__m128d a3_1 = _mm_load_sd(&values[20]);
#if defined(__SSE3__) && defined(__AVX__)
c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_1 = _mm_add_sd(c3_1, _mm_mul_sd(a3_1, b3));
#endif
_mm_store_sd(&C[(i*84)+3], c3_1);
__m128d c3_2 = _mm_load_sd(&C[(i*84)+9]);
__m128d a3_2 = _mm_load_sd(&values[21]);
#if defined(__SSE3__) && defined(__AVX__)
c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, b3));
#endif
_mm_store_sd(&C[(i*84)+9], c3_2);
__m128d c3_3 = _mm_load_sd(&C[(i*84)+19]);
__m128d a3_3 = _mm_load_sd(&values[22]);
#if defined(__SSE3__) && defined(__AVX__)
c3_3 = _mm_add_sd(c3_3, _mm_mul_sd(a3_3, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_3 = _mm_add_sd(c3_3, _mm_mul_sd(a3_3, b3));
#endif
_mm_store_sd(&C[(i*84)+19], c3_3);
__m128d c3_4 = _mm_load_sd(&C[(i*84)+34]);
__m128d a3_4 = _mm_load_sd(&values[23]);
#if defined(__SSE3__) && defined(__AVX__)
c3_4 = _mm_add_sd(c3_4, _mm_mul_sd(a3_4, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_4 = _mm_add_sd(c3_4, _mm_mul_sd(a3_4, b3));
#endif
_mm_store_sd(&C[(i*84)+34], c3_4);
__m128d c3_5 = _mm_load_sd(&C[(i*84)+55]);
__m128d a3_5 = _mm_load_sd(&values[24]);
#if defined(__SSE3__) && defined(__AVX__)
c3_5 = _mm_add_sd(c3_5, _mm_mul_sd(a3_5, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_5 = _mm_add_sd(c3_5, _mm_mul_sd(a3_5, b3));
#endif
_mm_store_sd(&C[(i*84)+55], c3_5);
__m128d c3_6 = _mm_load_sd(&C[(i*84)+83]);
__m128d a3_6 = _mm_load_sd(&values[25]);
#if defined(__SSE3__) && defined(__AVX__)
c3_6 = _mm_add_sd(c3_6, _mm_mul_sd(a3_6, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c3_6 = _mm_add_sd(c3_6, _mm_mul_sd(a3_6, b3));
#endif
_mm_store_sd(&C[(i*84)+83], c3_6);
#else
C[(i*84)+0] += values[19] * B[(i*84)+3];
C[(i*84)+3] += values[20] * B[(i*84)+3];
C[(i*84)+9] += values[21] * B[(i*84)+3];
C[(i*84)+19] += values[22] * B[(i*84)+3];
C[(i*84)+34] += values[23] * B[(i*84)+3];
C[(i*84)+55] += values[24] * B[(i*84)+3];
C[(i*84)+83] += values[25] * B[(i*84)+3];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b4 = _mm256_broadcast_sd(&B[(i*84)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b4 = _mm_loaddup_pd(&B[(i*84)+4]);
#endif
__m128d c4_0 = _mm_load_sd(&C[(i*84)+4]);
__m128d a4_0 = _mm_load_sd(&values[26]);
#if defined(__SSE3__) && defined(__AVX__)
c4_0 = _mm_add_sd(c4_0, _mm_mul_sd(a4_0, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_0 = _mm_add_sd(c4_0, _mm_mul_sd(a4_0, b4));
#endif
_mm_store_sd(&C[(i*84)+4], c4_0);
__m128d c4_1 = _mm_load_sd(&C[(i*84)+14]);
__m128d a4_1 = _mm_load_sd(&values[27]);
#if defined(__SSE3__) && defined(__AVX__)
c4_1 = _mm_add_sd(c4_1, _mm_mul_sd(a4_1, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_1 = _mm_add_sd(c4_1, _mm_mul_sd(a4_1, b4));
#endif
_mm_store_sd(&C[(i*84)+14], c4_1);
__m128d c4_2 = _mm_load_sd(&C[(i*84)+29]);
__m128d a4_2 = _mm_load_sd(&values[28]);
#if defined(__SSE3__) && defined(__AVX__)
c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, b4));
#endif
_mm_store_sd(&C[(i*84)+29], c4_2);
__m128d c4_3 = _mm_load_sd(&C[(i*84)+50]);
__m128d a4_3 = _mm_load_sd(&values[29]);
#if defined(__SSE3__) && defined(__AVX__)
c4_3 = _mm_add_sd(c4_3, _mm_mul_sd(a4_3, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_3 = _mm_add_sd(c4_3, _mm_mul_sd(a4_3, b4));
#endif
_mm_store_sd(&C[(i*84)+50], c4_3);
__m128d c4_4 = _mm_load_sd(&C[(i*84)+78]);
__m128d a4_4 = _mm_load_sd(&values[30]);
#if defined(__SSE3__) && defined(__AVX__)
c4_4 = _mm_add_sd(c4_4, _mm_mul_sd(a4_4, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c4_4 = _mm_add_sd(c4_4, _mm_mul_sd(a4_4, b4));
#endif
_mm_store_sd(&C[(i*84)+78], c4_4);
#else
C[(i*84)+4] += values[26] * B[(i*84)+4];
C[(i*84)+14] += values[27] * B[(i*84)+4];
C[(i*84)+29] += values[28] * B[(i*84)+4];
C[(i*84)+50] += values[29] * B[(i*84)+4];
C[(i*84)+78] += values[30] * B[(i*84)+4];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b5 = _mm256_broadcast_sd(&B[(i*84)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b5 = _mm_loaddup_pd(&B[(i*84)+5]);
#endif
__m128d c5_0 = _mm_load_sd(&C[(i*84)+5]);
__m128d a5_0 = _mm_load_sd(&values[31]);
#if defined(__SSE3__) && defined(__AVX__)
c5_0 = _mm_add_sd(c5_0, _mm_mul_sd(a5_0, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_0 = _mm_add_sd(c5_0, _mm_mul_sd(a5_0, b5));
#endif
_mm_store_sd(&C[(i*84)+5], c5_0);
__m128d c5_1 = _mm_load_sd(&C[(i*84)+15]);
__m128d a5_1 = _mm_load_sd(&values[32]);
#if defined(__SSE3__) && defined(__AVX__)
c5_1 = _mm_add_sd(c5_1, _mm_mul_sd(a5_1, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_1 = _mm_add_sd(c5_1, _mm_mul_sd(a5_1, b5));
#endif
_mm_store_sd(&C[(i*84)+15], c5_1);
__m128d c5_2 = _mm_load_sd(&C[(i*84)+30]);
__m128d a5_2 = _mm_load_sd(&values[33]);
#if defined(__SSE3__) && defined(__AVX__)
c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, b5));
#endif
_mm_store_sd(&C[(i*84)+30], c5_2);
__m128d c5_3 = _mm_load_sd(&C[(i*84)+51]);
__m128d a5_3 = _mm_load_sd(&values[34]);
#if defined(__SSE3__) && defined(__AVX__)
c5_3 = _mm_add_sd(c5_3, _mm_mul_sd(a5_3, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_3 = _mm_add_sd(c5_3, _mm_mul_sd(a5_3, b5));
#endif
_mm_store_sd(&C[(i*84)+51], c5_3);
__m128d c5_4 = _mm_load_sd(&C[(i*84)+79]);
__m128d a5_4 = _mm_load_sd(&values[35]);
#if defined(__SSE3__) && defined(__AVX__)
c5_4 = _mm_add_sd(c5_4, _mm_mul_sd(a5_4, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c5_4 = _mm_add_sd(c5_4, _mm_mul_sd(a5_4, b5));
#endif
_mm_store_sd(&C[(i*84)+79], c5_4);
#else
C[(i*84)+5] += values[31] * B[(i*84)+5];
C[(i*84)+15] += values[32] * B[(i*84)+5];
C[(i*84)+30] += values[33] * B[(i*84)+5];
C[(i*84)+51] += values[34] * B[(i*84)+5];
C[(i*84)+79] += values[35] * B[(i*84)+5];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b6 = _mm256_broadcast_sd(&B[(i*84)+6]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b6 = _mm_loaddup_pd(&B[(i*84)+6]);
#endif
__m128d c6_0 = _mm_load_sd(&C[(i*84)+6]);
__m128d a6_0 = _mm_load_sd(&values[36]);
#if defined(__SSE3__) && defined(__AVX__)
c6_0 = _mm_add_sd(c6_0, _mm_mul_sd(a6_0, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_0 = _mm_add_sd(c6_0, _mm_mul_sd(a6_0, b6));
#endif
_mm_store_sd(&C[(i*84)+6], c6_0);
__m128d c6_1 = _mm_load_sd(&C[(i*84)+16]);
__m128d a6_1 = _mm_load_sd(&values[37]);
#if defined(__SSE3__) && defined(__AVX__)
c6_1 = _mm_add_sd(c6_1, _mm_mul_sd(a6_1, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_1 = _mm_add_sd(c6_1, _mm_mul_sd(a6_1, b6));
#endif
_mm_store_sd(&C[(i*84)+16], c6_1);
__m128d c6_2 = _mm_load_sd(&C[(i*84)+31]);
__m128d a6_2 = _mm_load_sd(&values[38]);
#if defined(__SSE3__) && defined(__AVX__)
c6_2 = _mm_add_sd(c6_2, _mm_mul_sd(a6_2, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_2 = _mm_add_sd(c6_2, _mm_mul_sd(a6_2, b6));
#endif
_mm_store_sd(&C[(i*84)+31], c6_2);
__m128d c6_3 = _mm_load_sd(&C[(i*84)+52]);
__m128d a6_3 = _mm_load_sd(&values[39]);
#if defined(__SSE3__) && defined(__AVX__)
c6_3 = _mm_add_sd(c6_3, _mm_mul_sd(a6_3, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_3 = _mm_add_sd(c6_3, _mm_mul_sd(a6_3, b6));
#endif
_mm_store_sd(&C[(i*84)+52], c6_3);
__m128d c6_4 = _mm_load_sd(&C[(i*84)+80]);
__m128d a6_4 = _mm_load_sd(&values[40]);
#if defined(__SSE3__) && defined(__AVX__)
c6_4 = _mm_add_sd(c6_4, _mm_mul_sd(a6_4, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c6_4 = _mm_add_sd(c6_4, _mm_mul_sd(a6_4, b6));
#endif
_mm_store_sd(&C[(i*84)+80], c6_4);
#else
C[(i*84)+6] += values[36] * B[(i*84)+6];
C[(i*84)+16] += values[37] * B[(i*84)+6];
C[(i*84)+31] += values[38] * B[(i*84)+6];
C[(i*84)+52] += values[39] * B[(i*84)+6];
C[(i*84)+80] += values[40] * B[(i*84)+6];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b7 = _mm256_broadcast_sd(&B[(i*84)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b7 = _mm_loaddup_pd(&B[(i*84)+7]);
#endif
__m128d c7_0 = _mm_load_sd(&C[(i*84)+1]);
__m128d a7_0 = _mm_load_sd(&values[41]);
#if defined(__SSE3__) && defined(__AVX__)
c7_0 = _mm_add_sd(c7_0, _mm_mul_sd(a7_0, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_0 = _mm_add_sd(c7_0, _mm_mul_sd(a7_0, b7));
#endif
_mm_store_sd(&C[(i*84)+1], c7_0);
__m128d c7_1 = _mm_load_sd(&C[(i*84)+7]);
__m128d a7_1 = _mm_load_sd(&values[42]);
#if defined(__SSE3__) && defined(__AVX__)
c7_1 = _mm_add_sd(c7_1, _mm_mul_sd(a7_1, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_1 = _mm_add_sd(c7_1, _mm_mul_sd(a7_1, b7));
#endif
_mm_store_sd(&C[(i*84)+7], c7_1);
__m128d c7_2 = _mm_load_sd(&C[(i*84)+17]);
__m128d a7_2 = _mm_load_sd(&values[43]);
#if defined(__SSE3__) && defined(__AVX__)
c7_2 = _mm_add_sd(c7_2, _mm_mul_sd(a7_2, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_2 = _mm_add_sd(c7_2, _mm_mul_sd(a7_2, b7));
#endif
_mm_store_sd(&C[(i*84)+17], c7_2);
__m128d c7_3 = _mm_load_sd(&C[(i*84)+32]);
__m128d a7_3 = _mm_load_sd(&values[44]);
#if defined(__SSE3__) && defined(__AVX__)
c7_3 = _mm_add_sd(c7_3, _mm_mul_sd(a7_3, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_3 = _mm_add_sd(c7_3, _mm_mul_sd(a7_3, b7));
#endif
_mm_store_sd(&C[(i*84)+32], c7_3);
__m128d c7_4 = _mm_load_sd(&C[(i*84)+53]);
__m128d a7_4 = _mm_load_sd(&values[45]);
#if defined(__SSE3__) && defined(__AVX__)
c7_4 = _mm_add_sd(c7_4, _mm_mul_sd(a7_4, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_4 = _mm_add_sd(c7_4, _mm_mul_sd(a7_4, b7));
#endif
_mm_store_sd(&C[(i*84)+53], c7_4);
__m128d c7_5 = _mm_load_sd(&C[(i*84)+81]);
__m128d a7_5 = _mm_load_sd(&values[46]);
#if defined(__SSE3__) && defined(__AVX__)
c7_5 = _mm_add_sd(c7_5, _mm_mul_sd(a7_5, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c7_5 = _mm_add_sd(c7_5, _mm_mul_sd(a7_5, b7));
#endif
_mm_store_sd(&C[(i*84)+81], c7_5);
#else
C[(i*84)+1] += values[41] * B[(i*84)+7];
C[(i*84)+7] += values[42] * B[(i*84)+7];
C[(i*84)+17] += values[43] * B[(i*84)+7];
C[(i*84)+32] += values[44] * B[(i*84)+7];
C[(i*84)+53] += values[45] * B[(i*84)+7];
C[(i*84)+81] += values[46] * B[(i*84)+7];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b8 = _mm256_broadcast_sd(&B[(i*84)+8]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b8 = _mm_loaddup_pd(&B[(i*84)+8]);
#endif
__m128d c8_0 = _mm_load_sd(&C[(i*84)+2]);
__m128d a8_0 = _mm_load_sd(&values[47]);
#if defined(__SSE3__) && defined(__AVX__)
c8_0 = _mm_add_sd(c8_0, _mm_mul_sd(a8_0, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_0 = _mm_add_sd(c8_0, _mm_mul_sd(a8_0, b8));
#endif
_mm_store_sd(&C[(i*84)+2], c8_0);
__m128d c8_1 = _mm_load_sd(&C[(i*84)+8]);
__m128d a8_1 = _mm_load_sd(&values[48]);
#if defined(__SSE3__) && defined(__AVX__)
c8_1 = _mm_add_sd(c8_1, _mm_mul_sd(a8_1, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_1 = _mm_add_sd(c8_1, _mm_mul_sd(a8_1, b8));
#endif
_mm_store_sd(&C[(i*84)+8], c8_1);
__m128d c8_2 = _mm_load_sd(&C[(i*84)+18]);
__m128d a8_2 = _mm_load_sd(&values[49]);
#if defined(__SSE3__) && defined(__AVX__)
c8_2 = _mm_add_sd(c8_2, _mm_mul_sd(a8_2, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_2 = _mm_add_sd(c8_2, _mm_mul_sd(a8_2, b8));
#endif
_mm_store_sd(&C[(i*84)+18], c8_2);
__m128d c8_3 = _mm_load_sd(&C[(i*84)+33]);
__m128d a8_3 = _mm_load_sd(&values[50]);
#if defined(__SSE3__) && defined(__AVX__)
c8_3 = _mm_add_sd(c8_3, _mm_mul_sd(a8_3, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_3 = _mm_add_sd(c8_3, _mm_mul_sd(a8_3, b8));
#endif
_mm_store_sd(&C[(i*84)+33], c8_3);
__m128d c8_4 = _mm_load_sd(&C[(i*84)+54]);
__m128d a8_4 = _mm_load_sd(&values[51]);
#if defined(__SSE3__) && defined(__AVX__)
c8_4 = _mm_add_sd(c8_4, _mm_mul_sd(a8_4, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_4 = _mm_add_sd(c8_4, _mm_mul_sd(a8_4, b8));
#endif
_mm_store_sd(&C[(i*84)+54], c8_4);
__m128d c8_5 = _mm_load_sd(&C[(i*84)+82]);
__m128d a8_5 = _mm_load_sd(&values[52]);
#if defined(__SSE3__) && defined(__AVX__)
c8_5 = _mm_add_sd(c8_5, _mm_mul_sd(a8_5, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c8_5 = _mm_add_sd(c8_5, _mm_mul_sd(a8_5, b8));
#endif
_mm_store_sd(&C[(i*84)+82], c8_5);
#else
C[(i*84)+2] += values[47] * B[(i*84)+8];
C[(i*84)+8] += values[48] * B[(i*84)+8];
C[(i*84)+18] += values[49] * B[(i*84)+8];
C[(i*84)+33] += values[50] * B[(i*84)+8];
C[(i*84)+54] += values[51] * B[(i*84)+8];
C[(i*84)+82] += values[52] * B[(i*84)+8];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b9 = _mm256_broadcast_sd(&B[(i*84)+9]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b9 = _mm_loaddup_pd(&B[(i*84)+9]);
#endif
__m128d c9_0 = _mm_load_sd(&C[(i*84)+0]);
__m128d a9_0 = _mm_load_sd(&values[53]);
#if defined(__SSE3__) && defined(__AVX__)
c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, b9));
#endif
_mm_store_sd(&C[(i*84)+0], c9_0);
__m128d c9_1 = _mm_load_sd(&C[(i*84)+3]);
__m128d a9_1 = _mm_load_sd(&values[54]);
#if defined(__SSE3__) && defined(__AVX__)
c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_1 = _mm_add_sd(c9_1, _mm_mul_sd(a9_1, b9));
#endif
_mm_store_sd(&C[(i*84)+3], c9_1);
__m128d c9_2 = _mm_load_sd(&C[(i*84)+9]);
__m128d a9_2 = _mm_load_sd(&values[55]);
#if defined(__SSE3__) && defined(__AVX__)
c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_2 = _mm_add_sd(c9_2, _mm_mul_sd(a9_2, b9));
#endif
_mm_store_sd(&C[(i*84)+9], c9_2);
__m128d c9_3 = _mm_load_sd(&C[(i*84)+19]);
__m128d a9_3 = _mm_load_sd(&values[56]);
#if defined(__SSE3__) && defined(__AVX__)
c9_3 = _mm_add_sd(c9_3, _mm_mul_sd(a9_3, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_3 = _mm_add_sd(c9_3, _mm_mul_sd(a9_3, b9));
#endif
_mm_store_sd(&C[(i*84)+19], c9_3);
__m128d c9_4 = _mm_load_sd(&C[(i*84)+34]);
__m128d a9_4 = _mm_load_sd(&values[57]);
#if defined(__SSE3__) && defined(__AVX__)
c9_4 = _mm_add_sd(c9_4, _mm_mul_sd(a9_4, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_4 = _mm_add_sd(c9_4, _mm_mul_sd(a9_4, b9));
#endif
_mm_store_sd(&C[(i*84)+34], c9_4);
__m128d c9_5 = _mm_load_sd(&C[(i*84)+55]);
__m128d a9_5 = _mm_load_sd(&values[58]);
#if defined(__SSE3__) && defined(__AVX__)
c9_5 = _mm_add_sd(c9_5, _mm_mul_sd(a9_5, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_5 = _mm_add_sd(c9_5, _mm_mul_sd(a9_5, b9));
#endif
_mm_store_sd(&C[(i*84)+55], c9_5);
__m128d c9_6 = _mm_load_sd(&C[(i*84)+83]);
__m128d a9_6 = _mm_load_sd(&values[59]);
#if defined(__SSE3__) && defined(__AVX__)
c9_6 = _mm_add_sd(c9_6, _mm_mul_sd(a9_6, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c9_6 = _mm_add_sd(c9_6, _mm_mul_sd(a9_6, b9));
#endif
_mm_store_sd(&C[(i*84)+83], c9_6);
#else
C[(i*84)+0] += values[53] * B[(i*84)+9];
C[(i*84)+3] += values[54] * B[(i*84)+9];
C[(i*84)+9] += values[55] * B[(i*84)+9];
C[(i*84)+19] += values[56] * B[(i*84)+9];
C[(i*84)+34] += values[57] * B[(i*84)+9];
C[(i*84)+55] += values[58] * B[(i*84)+9];
C[(i*84)+83] += values[59] * B[(i*84)+9];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b10 = _mm256_broadcast_sd(&B[(i*84)+10]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b10 = _mm_loaddup_pd(&B[(i*84)+10]);
#endif
__m128d c10_0 = _mm_load_sd(&C[(i*84)+10]);
__m128d a10_0 = _mm_load_sd(&values[60]);
#if defined(__SSE3__) && defined(__AVX__)
c10_0 = _mm_add_sd(c10_0, _mm_mul_sd(a10_0, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c10_0 = _mm_add_sd(c10_0, _mm_mul_sd(a10_0, b10));
#endif
_mm_store_sd(&C[(i*84)+10], c10_0);
__m128d c10_1 = _mm_load_sd(&C[(i*84)+25]);
__m128d a10_1 = _mm_load_sd(&values[61]);
#if defined(__SSE3__) && defined(__AVX__)
c10_1 = _mm_add_sd(c10_1, _mm_mul_sd(a10_1, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c10_1 = _mm_add_sd(c10_1, _mm_mul_sd(a10_1, b10));
#endif
_mm_store_sd(&C[(i*84)+25], c10_1);
__m128d c10_2 = _mm_load_sd(&C[(i*84)+46]);
__m128d a10_2 = _mm_load_sd(&values[62]);
#if defined(__SSE3__) && defined(__AVX__)
c10_2 = _mm_add_sd(c10_2, _mm_mul_sd(a10_2, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c10_2 = _mm_add_sd(c10_2, _mm_mul_sd(a10_2, b10));
#endif
_mm_store_sd(&C[(i*84)+46], c10_2);
__m128d c10_3 = _mm_load_sd(&C[(i*84)+74]);
__m128d a10_3 = _mm_load_sd(&values[63]);
#if defined(__SSE3__) && defined(__AVX__)
c10_3 = _mm_add_sd(c10_3, _mm_mul_sd(a10_3, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c10_3 = _mm_add_sd(c10_3, _mm_mul_sd(a10_3, b10));
#endif
_mm_store_sd(&C[(i*84)+74], c10_3);
#else
C[(i*84)+10] += values[60] * B[(i*84)+10];
C[(i*84)+25] += values[61] * B[(i*84)+10];
C[(i*84)+46] += values[62] * B[(i*84)+10];
C[(i*84)+74] += values[63] * B[(i*84)+10];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b11 = _mm256_broadcast_sd(&B[(i*84)+11]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b11 = _mm_loaddup_pd(&B[(i*84)+11]);
#endif
__m128d c11_0 = _mm_load_sd(&C[(i*84)+11]);
__m128d a11_0 = _mm_load_sd(&values[64]);
#if defined(__SSE3__) && defined(__AVX__)
c11_0 = _mm_add_sd(c11_0, _mm_mul_sd(a11_0, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c11_0 = _mm_add_sd(c11_0, _mm_mul_sd(a11_0, b11));
#endif
_mm_store_sd(&C[(i*84)+11], c11_0);
__m128d c11_1 = _mm_load_sd(&C[(i*84)+26]);
__m128d a11_1 = _mm_load_sd(&values[65]);
#if defined(__SSE3__) && defined(__AVX__)
c11_1 = _mm_add_sd(c11_1, _mm_mul_sd(a11_1, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c11_1 = _mm_add_sd(c11_1, _mm_mul_sd(a11_1, b11));
#endif
_mm_store_sd(&C[(i*84)+26], c11_1);
__m128d c11_2 = _mm_load_sd(&C[(i*84)+47]);
__m128d a11_2 = _mm_load_sd(&values[66]);
#if defined(__SSE3__) && defined(__AVX__)
c11_2 = _mm_add_sd(c11_2, _mm_mul_sd(a11_2, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c11_2 = _mm_add_sd(c11_2, _mm_mul_sd(a11_2, b11));
#endif
_mm_store_sd(&C[(i*84)+47], c11_2);
__m128d c11_3 = _mm_load_sd(&C[(i*84)+75]);
__m128d a11_3 = _mm_load_sd(&values[67]);
#if defined(__SSE3__) && defined(__AVX__)
c11_3 = _mm_add_sd(c11_3, _mm_mul_sd(a11_3, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c11_3 = _mm_add_sd(c11_3, _mm_mul_sd(a11_3, b11));
#endif
_mm_store_sd(&C[(i*84)+75], c11_3);
#else
C[(i*84)+11] += values[64] * B[(i*84)+11];
C[(i*84)+26] += values[65] * B[(i*84)+11];
C[(i*84)+47] += values[66] * B[(i*84)+11];
C[(i*84)+75] += values[67] * B[(i*84)+11];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b12 = _mm256_broadcast_sd(&B[(i*84)+12]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b12 = _mm_loaddup_pd(&B[(i*84)+12]);
#endif
__m128d c12_0 = _mm_load_sd(&C[(i*84)+12]);
__m128d a12_0 = _mm_load_sd(&values[68]);
#if defined(__SSE3__) && defined(__AVX__)
c12_0 = _mm_add_sd(c12_0, _mm_mul_sd(a12_0, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c12_0 = _mm_add_sd(c12_0, _mm_mul_sd(a12_0, b12));
#endif
_mm_store_sd(&C[(i*84)+12], c12_0);
__m128d c12_1 = _mm_load_sd(&C[(i*84)+27]);
__m128d a12_1 = _mm_load_sd(&values[69]);
#if defined(__SSE3__) && defined(__AVX__)
c12_1 = _mm_add_sd(c12_1, _mm_mul_sd(a12_1, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c12_1 = _mm_add_sd(c12_1, _mm_mul_sd(a12_1, b12));
#endif
_mm_store_sd(&C[(i*84)+27], c12_1);
__m128d c12_2 = _mm_load_sd(&C[(i*84)+48]);
__m128d a12_2 = _mm_load_sd(&values[70]);
#if defined(__SSE3__) && defined(__AVX__)
c12_2 = _mm_add_sd(c12_2, _mm_mul_sd(a12_2, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c12_2 = _mm_add_sd(c12_2, _mm_mul_sd(a12_2, b12));
#endif
_mm_store_sd(&C[(i*84)+48], c12_2);
__m128d c12_3 = _mm_load_sd(&C[(i*84)+76]);
__m128d a12_3 = _mm_load_sd(&values[71]);
#if defined(__SSE3__) && defined(__AVX__)
c12_3 = _mm_add_sd(c12_3, _mm_mul_sd(a12_3, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c12_3 = _mm_add_sd(c12_3, _mm_mul_sd(a12_3, b12));
#endif
_mm_store_sd(&C[(i*84)+76], c12_3);
#else
C[(i*84)+12] += values[68] * B[(i*84)+12];
C[(i*84)+27] += values[69] * B[(i*84)+12];
C[(i*84)+48] += values[70] * B[(i*84)+12];
C[(i*84)+76] += values[71] * B[(i*84)+12];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b13 = _mm256_broadcast_sd(&B[(i*84)+13]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b13 = _mm_loaddup_pd(&B[(i*84)+13]);
#endif
__m128d c13_0 = _mm_load_sd(&C[(i*84)+13]);
__m128d a13_0 = _mm_load_sd(&values[72]);
#if defined(__SSE3__) && defined(__AVX__)
c13_0 = _mm_add_sd(c13_0, _mm_mul_sd(a13_0, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c13_0 = _mm_add_sd(c13_0, _mm_mul_sd(a13_0, b13));
#endif
_mm_store_sd(&C[(i*84)+13], c13_0);
__m128d c13_1 = _mm_load_sd(&C[(i*84)+28]);
__m128d a13_1 = _mm_load_sd(&values[73]);
#if defined(__SSE3__) && defined(__AVX__)
c13_1 = _mm_add_sd(c13_1, _mm_mul_sd(a13_1, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c13_1 = _mm_add_sd(c13_1, _mm_mul_sd(a13_1, b13));
#endif
_mm_store_sd(&C[(i*84)+28], c13_1);
__m128d c13_2 = _mm_load_sd(&C[(i*84)+49]);
__m128d a13_2 = _mm_load_sd(&values[74]);
#if defined(__SSE3__) && defined(__AVX__)
c13_2 = _mm_add_sd(c13_2, _mm_mul_sd(a13_2, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c13_2 = _mm_add_sd(c13_2, _mm_mul_sd(a13_2, b13));
#endif
_mm_store_sd(&C[(i*84)+49], c13_2);
__m128d c13_3 = _mm_load_sd(&C[(i*84)+77]);
__m128d a13_3 = _mm_load_sd(&values[75]);
#if defined(__SSE3__) && defined(__AVX__)
c13_3 = _mm_add_sd(c13_3, _mm_mul_sd(a13_3, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c13_3 = _mm_add_sd(c13_3, _mm_mul_sd(a13_3, b13));
#endif
_mm_store_sd(&C[(i*84)+77], c13_3);
#else
C[(i*84)+13] += values[72] * B[(i*84)+13];
C[(i*84)+28] += values[73] * B[(i*84)+13];
C[(i*84)+49] += values[74] * B[(i*84)+13];
C[(i*84)+77] += values[75] * B[(i*84)+13];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b14 = _mm256_broadcast_sd(&B[(i*84)+14]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b14 = _mm_loaddup_pd(&B[(i*84)+14]);
#endif
__m128d c14_0 = _mm_load_sd(&C[(i*84)+4]);
__m128d a14_0 = _mm_load_sd(&values[76]);
#if defined(__SSE3__) && defined(__AVX__)
c14_0 = _mm_add_sd(c14_0, _mm_mul_sd(a14_0, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_0 = _mm_add_sd(c14_0, _mm_mul_sd(a14_0, b14));
#endif
_mm_store_sd(&C[(i*84)+4], c14_0);
__m128d c14_1 = _mm_load_sd(&C[(i*84)+14]);
__m128d a14_1 = _mm_load_sd(&values[77]);
#if defined(__SSE3__) && defined(__AVX__)
c14_1 = _mm_add_sd(c14_1, _mm_mul_sd(a14_1, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_1 = _mm_add_sd(c14_1, _mm_mul_sd(a14_1, b14));
#endif
_mm_store_sd(&C[(i*84)+14], c14_1);
__m128d c14_2 = _mm_load_sd(&C[(i*84)+29]);
__m128d a14_2 = _mm_load_sd(&values[78]);
#if defined(__SSE3__) && defined(__AVX__)
c14_2 = _mm_add_sd(c14_2, _mm_mul_sd(a14_2, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_2 = _mm_add_sd(c14_2, _mm_mul_sd(a14_2, b14));
#endif
_mm_store_sd(&C[(i*84)+29], c14_2);
__m128d c14_3 = _mm_load_sd(&C[(i*84)+50]);
__m128d a14_3 = _mm_load_sd(&values[79]);
#if defined(__SSE3__) && defined(__AVX__)
c14_3 = _mm_add_sd(c14_3, _mm_mul_sd(a14_3, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_3 = _mm_add_sd(c14_3, _mm_mul_sd(a14_3, b14));
#endif
_mm_store_sd(&C[(i*84)+50], c14_3);
__m128d c14_4 = _mm_load_sd(&C[(i*84)+78]);
__m128d a14_4 = _mm_load_sd(&values[80]);
#if defined(__SSE3__) && defined(__AVX__)
c14_4 = _mm_add_sd(c14_4, _mm_mul_sd(a14_4, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c14_4 = _mm_add_sd(c14_4, _mm_mul_sd(a14_4, b14));
#endif
_mm_store_sd(&C[(i*84)+78], c14_4);
#else
C[(i*84)+4] += values[76] * B[(i*84)+14];
C[(i*84)+14] += values[77] * B[(i*84)+14];
C[(i*84)+29] += values[78] * B[(i*84)+14];
C[(i*84)+50] += values[79] * B[(i*84)+14];
C[(i*84)+78] += values[80] * B[(i*84)+14];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b15 = _mm256_broadcast_sd(&B[(i*84)+15]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b15 = _mm_loaddup_pd(&B[(i*84)+15]);
#endif
__m128d c15_0 = _mm_load_sd(&C[(i*84)+5]);
__m128d a15_0 = _mm_load_sd(&values[81]);
#if defined(__SSE3__) && defined(__AVX__)
c15_0 = _mm_add_sd(c15_0, _mm_mul_sd(a15_0, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_0 = _mm_add_sd(c15_0, _mm_mul_sd(a15_0, b15));
#endif
_mm_store_sd(&C[(i*84)+5], c15_0);
__m128d c15_1 = _mm_load_sd(&C[(i*84)+15]);
__m128d a15_1 = _mm_load_sd(&values[82]);
#if defined(__SSE3__) && defined(__AVX__)
c15_1 = _mm_add_sd(c15_1, _mm_mul_sd(a15_1, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_1 = _mm_add_sd(c15_1, _mm_mul_sd(a15_1, b15));
#endif
_mm_store_sd(&C[(i*84)+15], c15_1);
__m128d c15_2 = _mm_load_sd(&C[(i*84)+30]);
__m128d a15_2 = _mm_load_sd(&values[83]);
#if defined(__SSE3__) && defined(__AVX__)
c15_2 = _mm_add_sd(c15_2, _mm_mul_sd(a15_2, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_2 = _mm_add_sd(c15_2, _mm_mul_sd(a15_2, b15));
#endif
_mm_store_sd(&C[(i*84)+30], c15_2);
__m128d c15_3 = _mm_load_sd(&C[(i*84)+51]);
__m128d a15_3 = _mm_load_sd(&values[84]);
#if defined(__SSE3__) && defined(__AVX__)
c15_3 = _mm_add_sd(c15_3, _mm_mul_sd(a15_3, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_3 = _mm_add_sd(c15_3, _mm_mul_sd(a15_3, b15));
#endif
_mm_store_sd(&C[(i*84)+51], c15_3);
__m128d c15_4 = _mm_load_sd(&C[(i*84)+79]);
__m128d a15_4 = _mm_load_sd(&values[85]);
#if defined(__SSE3__) && defined(__AVX__)
c15_4 = _mm_add_sd(c15_4, _mm_mul_sd(a15_4, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c15_4 = _mm_add_sd(c15_4, _mm_mul_sd(a15_4, b15));
#endif
_mm_store_sd(&C[(i*84)+79], c15_4);
#else
C[(i*84)+5] += values[81] * B[(i*84)+15];
C[(i*84)+15] += values[82] * B[(i*84)+15];
C[(i*84)+30] += values[83] * B[(i*84)+15];
C[(i*84)+51] += values[84] * B[(i*84)+15];
C[(i*84)+79] += values[85] * B[(i*84)+15];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b16 = _mm256_broadcast_sd(&B[(i*84)+16]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b16 = _mm_loaddup_pd(&B[(i*84)+16]);
#endif
__m128d c16_0 = _mm_load_sd(&C[(i*84)+6]);
__m128d a16_0 = _mm_load_sd(&values[86]);
#if defined(__SSE3__) && defined(__AVX__)
c16_0 = _mm_add_sd(c16_0, _mm_mul_sd(a16_0, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_0 = _mm_add_sd(c16_0, _mm_mul_sd(a16_0, b16));
#endif
_mm_store_sd(&C[(i*84)+6], c16_0);
__m128d c16_1 = _mm_load_sd(&C[(i*84)+16]);
__m128d a16_1 = _mm_load_sd(&values[87]);
#if defined(__SSE3__) && defined(__AVX__)
c16_1 = _mm_add_sd(c16_1, _mm_mul_sd(a16_1, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_1 = _mm_add_sd(c16_1, _mm_mul_sd(a16_1, b16));
#endif
_mm_store_sd(&C[(i*84)+16], c16_1);
__m128d c16_2 = _mm_load_sd(&C[(i*84)+31]);
__m128d a16_2 = _mm_load_sd(&values[88]);
#if defined(__SSE3__) && defined(__AVX__)
c16_2 = _mm_add_sd(c16_2, _mm_mul_sd(a16_2, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_2 = _mm_add_sd(c16_2, _mm_mul_sd(a16_2, b16));
#endif
_mm_store_sd(&C[(i*84)+31], c16_2);
__m128d c16_3 = _mm_load_sd(&C[(i*84)+52]);
__m128d a16_3 = _mm_load_sd(&values[89]);
#if defined(__SSE3__) && defined(__AVX__)
c16_3 = _mm_add_sd(c16_3, _mm_mul_sd(a16_3, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_3 = _mm_add_sd(c16_3, _mm_mul_sd(a16_3, b16));
#endif
_mm_store_sd(&C[(i*84)+52], c16_3);
__m128d c16_4 = _mm_load_sd(&C[(i*84)+80]);
__m128d a16_4 = _mm_load_sd(&values[90]);
#if defined(__SSE3__) && defined(__AVX__)
c16_4 = _mm_add_sd(c16_4, _mm_mul_sd(a16_4, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c16_4 = _mm_add_sd(c16_4, _mm_mul_sd(a16_4, b16));
#endif
_mm_store_sd(&C[(i*84)+80], c16_4);
#else
C[(i*84)+6] += values[86] * B[(i*84)+16];
C[(i*84)+16] += values[87] * B[(i*84)+16];
C[(i*84)+31] += values[88] * B[(i*84)+16];
C[(i*84)+52] += values[89] * B[(i*84)+16];
C[(i*84)+80] += values[90] * B[(i*84)+16];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b17 = _mm256_broadcast_sd(&B[(i*84)+17]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b17 = _mm_loaddup_pd(&B[(i*84)+17]);
#endif
__m128d c17_0 = _mm_load_sd(&C[(i*84)+1]);
__m128d a17_0 = _mm_load_sd(&values[91]);
#if defined(__SSE3__) && defined(__AVX__)
c17_0 = _mm_add_sd(c17_0, _mm_mul_sd(a17_0, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_0 = _mm_add_sd(c17_0, _mm_mul_sd(a17_0, b17));
#endif
_mm_store_sd(&C[(i*84)+1], c17_0);
__m128d c17_1 = _mm_load_sd(&C[(i*84)+7]);
__m128d a17_1 = _mm_load_sd(&values[92]);
#if defined(__SSE3__) && defined(__AVX__)
c17_1 = _mm_add_sd(c17_1, _mm_mul_sd(a17_1, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_1 = _mm_add_sd(c17_1, _mm_mul_sd(a17_1, b17));
#endif
_mm_store_sd(&C[(i*84)+7], c17_1);
__m128d c17_2 = _mm_load_sd(&C[(i*84)+17]);
__m128d a17_2 = _mm_load_sd(&values[93]);
#if defined(__SSE3__) && defined(__AVX__)
c17_2 = _mm_add_sd(c17_2, _mm_mul_sd(a17_2, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_2 = _mm_add_sd(c17_2, _mm_mul_sd(a17_2, b17));
#endif
_mm_store_sd(&C[(i*84)+17], c17_2);
__m128d c17_3 = _mm_load_sd(&C[(i*84)+32]);
__m128d a17_3 = _mm_load_sd(&values[94]);
#if defined(__SSE3__) && defined(__AVX__)
c17_3 = _mm_add_sd(c17_3, _mm_mul_sd(a17_3, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_3 = _mm_add_sd(c17_3, _mm_mul_sd(a17_3, b17));
#endif
_mm_store_sd(&C[(i*84)+32], c17_3);
__m128d c17_4 = _mm_load_sd(&C[(i*84)+53]);
__m128d a17_4 = _mm_load_sd(&values[95]);
#if defined(__SSE3__) && defined(__AVX__)
c17_4 = _mm_add_sd(c17_4, _mm_mul_sd(a17_4, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_4 = _mm_add_sd(c17_4, _mm_mul_sd(a17_4, b17));
#endif
_mm_store_sd(&C[(i*84)+53], c17_4);
__m128d c17_5 = _mm_load_sd(&C[(i*84)+81]);
__m128d a17_5 = _mm_load_sd(&values[96]);
#if defined(__SSE3__) && defined(__AVX__)
c17_5 = _mm_add_sd(c17_5, _mm_mul_sd(a17_5, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c17_5 = _mm_add_sd(c17_5, _mm_mul_sd(a17_5, b17));
#endif
_mm_store_sd(&C[(i*84)+81], c17_5);
#else
C[(i*84)+1] += values[91] * B[(i*84)+17];
C[(i*84)+7] += values[92] * B[(i*84)+17];
C[(i*84)+17] += values[93] * B[(i*84)+17];
C[(i*84)+32] += values[94] * B[(i*84)+17];
C[(i*84)+53] += values[95] * B[(i*84)+17];
C[(i*84)+81] += values[96] * B[(i*84)+17];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b18 = _mm256_broadcast_sd(&B[(i*84)+18]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b18 = _mm_loaddup_pd(&B[(i*84)+18]);
#endif
__m128d c18_0 = _mm_load_sd(&C[(i*84)+2]);
__m128d a18_0 = _mm_load_sd(&values[97]);
#if defined(__SSE3__) && defined(__AVX__)
c18_0 = _mm_add_sd(c18_0, _mm_mul_sd(a18_0, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_0 = _mm_add_sd(c18_0, _mm_mul_sd(a18_0, b18));
#endif
_mm_store_sd(&C[(i*84)+2], c18_0);
__m128d c18_1 = _mm_load_sd(&C[(i*84)+8]);
__m128d a18_1 = _mm_load_sd(&values[98]);
#if defined(__SSE3__) && defined(__AVX__)
c18_1 = _mm_add_sd(c18_1, _mm_mul_sd(a18_1, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_1 = _mm_add_sd(c18_1, _mm_mul_sd(a18_1, b18));
#endif
_mm_store_sd(&C[(i*84)+8], c18_1);
__m128d c18_2 = _mm_load_sd(&C[(i*84)+18]);
__m128d a18_2 = _mm_load_sd(&values[99]);
#if defined(__SSE3__) && defined(__AVX__)
c18_2 = _mm_add_sd(c18_2, _mm_mul_sd(a18_2, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_2 = _mm_add_sd(c18_2, _mm_mul_sd(a18_2, b18));
#endif
_mm_store_sd(&C[(i*84)+18], c18_2);
__m128d c18_3 = _mm_load_sd(&C[(i*84)+33]);
__m128d a18_3 = _mm_load_sd(&values[100]);
#if defined(__SSE3__) && defined(__AVX__)
c18_3 = _mm_add_sd(c18_3, _mm_mul_sd(a18_3, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_3 = _mm_add_sd(c18_3, _mm_mul_sd(a18_3, b18));
#endif
_mm_store_sd(&C[(i*84)+33], c18_3);
__m128d c18_4 = _mm_load_sd(&C[(i*84)+54]);
__m128d a18_4 = _mm_load_sd(&values[101]);
#if defined(__SSE3__) && defined(__AVX__)
c18_4 = _mm_add_sd(c18_4, _mm_mul_sd(a18_4, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_4 = _mm_add_sd(c18_4, _mm_mul_sd(a18_4, b18));
#endif
_mm_store_sd(&C[(i*84)+54], c18_4);
__m128d c18_5 = _mm_load_sd(&C[(i*84)+82]);
__m128d a18_5 = _mm_load_sd(&values[102]);
#if defined(__SSE3__) && defined(__AVX__)
c18_5 = _mm_add_sd(c18_5, _mm_mul_sd(a18_5, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c18_5 = _mm_add_sd(c18_5, _mm_mul_sd(a18_5, b18));
#endif
_mm_store_sd(&C[(i*84)+82], c18_5);
#else
C[(i*84)+2] += values[97] * B[(i*84)+18];
C[(i*84)+8] += values[98] * B[(i*84)+18];
C[(i*84)+18] += values[99] * B[(i*84)+18];
C[(i*84)+33] += values[100] * B[(i*84)+18];
C[(i*84)+54] += values[101] * B[(i*84)+18];
C[(i*84)+82] += values[102] * B[(i*84)+18];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b19 = _mm256_broadcast_sd(&B[(i*84)+19]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b19 = _mm_loaddup_pd(&B[(i*84)+19]);
#endif
__m128d c19_0 = _mm_load_sd(&C[(i*84)+0]);
__m128d a19_0 = _mm_load_sd(&values[103]);
#if defined(__SSE3__) && defined(__AVX__)
c19_0 = _mm_add_sd(c19_0, _mm_mul_sd(a19_0, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_0 = _mm_add_sd(c19_0, _mm_mul_sd(a19_0, b19));
#endif
_mm_store_sd(&C[(i*84)+0], c19_0);
__m128d c19_1 = _mm_load_sd(&C[(i*84)+3]);
__m128d a19_1 = _mm_load_sd(&values[104]);
#if defined(__SSE3__) && defined(__AVX__)
c19_1 = _mm_add_sd(c19_1, _mm_mul_sd(a19_1, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_1 = _mm_add_sd(c19_1, _mm_mul_sd(a19_1, b19));
#endif
_mm_store_sd(&C[(i*84)+3], c19_1);
__m128d c19_2 = _mm_load_sd(&C[(i*84)+9]);
__m128d a19_2 = _mm_load_sd(&values[105]);
#if defined(__SSE3__) && defined(__AVX__)
c19_2 = _mm_add_sd(c19_2, _mm_mul_sd(a19_2, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_2 = _mm_add_sd(c19_2, _mm_mul_sd(a19_2, b19));
#endif
_mm_store_sd(&C[(i*84)+9], c19_2);
__m128d c19_3 = _mm_load_sd(&C[(i*84)+19]);
__m128d a19_3 = _mm_load_sd(&values[106]);
#if defined(__SSE3__) && defined(__AVX__)
c19_3 = _mm_add_sd(c19_3, _mm_mul_sd(a19_3, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_3 = _mm_add_sd(c19_3, _mm_mul_sd(a19_3, b19));
#endif
_mm_store_sd(&C[(i*84)+19], c19_3);
__m128d c19_4 = _mm_load_sd(&C[(i*84)+34]);
__m128d a19_4 = _mm_load_sd(&values[107]);
#if defined(__SSE3__) && defined(__AVX__)
c19_4 = _mm_add_sd(c19_4, _mm_mul_sd(a19_4, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_4 = _mm_add_sd(c19_4, _mm_mul_sd(a19_4, b19));
#endif
_mm_store_sd(&C[(i*84)+34], c19_4);
__m128d c19_5 = _mm_load_sd(&C[(i*84)+55]);
__m128d a19_5 = _mm_load_sd(&values[108]);
#if defined(__SSE3__) && defined(__AVX__)
c19_5 = _mm_add_sd(c19_5, _mm_mul_sd(a19_5, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_5 = _mm_add_sd(c19_5, _mm_mul_sd(a19_5, b19));
#endif
_mm_store_sd(&C[(i*84)+55], c19_5);
__m128d c19_6 = _mm_load_sd(&C[(i*84)+83]);
__m128d a19_6 = _mm_load_sd(&values[109]);
#if defined(__SSE3__) && defined(__AVX__)
c19_6 = _mm_add_sd(c19_6, _mm_mul_sd(a19_6, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c19_6 = _mm_add_sd(c19_6, _mm_mul_sd(a19_6, b19));
#endif
_mm_store_sd(&C[(i*84)+83], c19_6);
#else
C[(i*84)+0] += values[103] * B[(i*84)+19];
C[(i*84)+3] += values[104] * B[(i*84)+19];
C[(i*84)+9] += values[105] * B[(i*84)+19];
C[(i*84)+19] += values[106] * B[(i*84)+19];
C[(i*84)+34] += values[107] * B[(i*84)+19];
C[(i*84)+55] += values[108] * B[(i*84)+19];
C[(i*84)+83] += values[109] * B[(i*84)+19];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b20 = _mm256_broadcast_sd(&B[(i*84)+20]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b20 = _mm_loaddup_pd(&B[(i*84)+20]);
#endif
__m128d c20_0 = _mm_load_sd(&C[(i*84)+20]);
__m128d a20_0 = _mm_load_sd(&values[110]);
#if defined(__SSE3__) && defined(__AVX__)
c20_0 = _mm_add_sd(c20_0, _mm_mul_sd(a20_0, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c20_0 = _mm_add_sd(c20_0, _mm_mul_sd(a20_0, b20));
#endif
_mm_store_sd(&C[(i*84)+20], c20_0);
__m128d c20_1 = _mm_load_sd(&C[(i*84)+41]);
__m128d a20_1 = _mm_load_sd(&values[111]);
#if defined(__SSE3__) && defined(__AVX__)
c20_1 = _mm_add_sd(c20_1, _mm_mul_sd(a20_1, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c20_1 = _mm_add_sd(c20_1, _mm_mul_sd(a20_1, b20));
#endif
_mm_store_sd(&C[(i*84)+41], c20_1);
__m128d c20_2 = _mm_load_sd(&C[(i*84)+69]);
__m128d a20_2 = _mm_load_sd(&values[112]);
#if defined(__SSE3__) && defined(__AVX__)
c20_2 = _mm_add_sd(c20_2, _mm_mul_sd(a20_2, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c20_2 = _mm_add_sd(c20_2, _mm_mul_sd(a20_2, b20));
#endif
_mm_store_sd(&C[(i*84)+69], c20_2);
#else
C[(i*84)+20] += values[110] * B[(i*84)+20];
C[(i*84)+41] += values[111] * B[(i*84)+20];
C[(i*84)+69] += values[112] * B[(i*84)+20];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b21 = _mm256_broadcast_sd(&B[(i*84)+21]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b21 = _mm_loaddup_pd(&B[(i*84)+21]);
#endif
__m128d c21_0 = _mm_load_sd(&C[(i*84)+21]);
__m128d a21_0 = _mm_load_sd(&values[113]);
#if defined(__SSE3__) && defined(__AVX__)
c21_0 = _mm_add_sd(c21_0, _mm_mul_sd(a21_0, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c21_0 = _mm_add_sd(c21_0, _mm_mul_sd(a21_0, b21));
#endif
_mm_store_sd(&C[(i*84)+21], c21_0);
__m128d c21_1 = _mm_load_sd(&C[(i*84)+42]);
__m128d a21_1 = _mm_load_sd(&values[114]);
#if defined(__SSE3__) && defined(__AVX__)
c21_1 = _mm_add_sd(c21_1, _mm_mul_sd(a21_1, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c21_1 = _mm_add_sd(c21_1, _mm_mul_sd(a21_1, b21));
#endif
_mm_store_sd(&C[(i*84)+42], c21_1);
__m128d c21_2 = _mm_load_sd(&C[(i*84)+70]);
__m128d a21_2 = _mm_load_sd(&values[115]);
#if defined(__SSE3__) && defined(__AVX__)
c21_2 = _mm_add_sd(c21_2, _mm_mul_sd(a21_2, _mm256_castpd256_pd128(b21)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c21_2 = _mm_add_sd(c21_2, _mm_mul_sd(a21_2, b21));
#endif
_mm_store_sd(&C[(i*84)+70], c21_2);
#else
C[(i*84)+21] += values[113] * B[(i*84)+21];
C[(i*84)+42] += values[114] * B[(i*84)+21];
C[(i*84)+70] += values[115] * B[(i*84)+21];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b22 = _mm256_broadcast_sd(&B[(i*84)+22]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b22 = _mm_loaddup_pd(&B[(i*84)+22]);
#endif
__m128d c22_0 = _mm_load_sd(&C[(i*84)+22]);
__m128d a22_0 = _mm_load_sd(&values[116]);
#if defined(__SSE3__) && defined(__AVX__)
c22_0 = _mm_add_sd(c22_0, _mm_mul_sd(a22_0, _mm256_castpd256_pd128(b22)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c22_0 = _mm_add_sd(c22_0, _mm_mul_sd(a22_0, b22));
#endif
_mm_store_sd(&C[(i*84)+22], c22_0);
__m128d c22_1 = _mm_load_sd(&C[(i*84)+43]);
__m128d a22_1 = _mm_load_sd(&values[117]);
#if defined(__SSE3__) && defined(__AVX__)
c22_1 = _mm_add_sd(c22_1, _mm_mul_sd(a22_1, _mm256_castpd256_pd128(b22)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c22_1 = _mm_add_sd(c22_1, _mm_mul_sd(a22_1, b22));
#endif
_mm_store_sd(&C[(i*84)+43], c22_1);
__m128d c22_2 = _mm_load_sd(&C[(i*84)+71]);
__m128d a22_2 = _mm_load_sd(&values[118]);
#if defined(__SSE3__) && defined(__AVX__)
c22_2 = _mm_add_sd(c22_2, _mm_mul_sd(a22_2, _mm256_castpd256_pd128(b22)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c22_2 = _mm_add_sd(c22_2, _mm_mul_sd(a22_2, b22));
#endif
_mm_store_sd(&C[(i*84)+71], c22_2);
#else
C[(i*84)+22] += values[116] * B[(i*84)+22];
C[(i*84)+43] += values[117] * B[(i*84)+22];
C[(i*84)+71] += values[118] * B[(i*84)+22];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b23 = _mm256_broadcast_sd(&B[(i*84)+23]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b23 = _mm_loaddup_pd(&B[(i*84)+23]);
#endif
__m128d c23_0 = _mm_load_sd(&C[(i*84)+23]);
__m128d a23_0 = _mm_load_sd(&values[119]);
#if defined(__SSE3__) && defined(__AVX__)
c23_0 = _mm_add_sd(c23_0, _mm_mul_sd(a23_0, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c23_0 = _mm_add_sd(c23_0, _mm_mul_sd(a23_0, b23));
#endif
_mm_store_sd(&C[(i*84)+23], c23_0);
__m128d c23_1 = _mm_load_sd(&C[(i*84)+44]);
__m128d a23_1 = _mm_load_sd(&values[120]);
#if defined(__SSE3__) && defined(__AVX__)
c23_1 = _mm_add_sd(c23_1, _mm_mul_sd(a23_1, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c23_1 = _mm_add_sd(c23_1, _mm_mul_sd(a23_1, b23));
#endif
_mm_store_sd(&C[(i*84)+44], c23_1);
__m128d c23_2 = _mm_load_sd(&C[(i*84)+72]);
__m128d a23_2 = _mm_load_sd(&values[121]);
#if defined(__SSE3__) && defined(__AVX__)
c23_2 = _mm_add_sd(c23_2, _mm_mul_sd(a23_2, _mm256_castpd256_pd128(b23)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c23_2 = _mm_add_sd(c23_2, _mm_mul_sd(a23_2, b23));
#endif
_mm_store_sd(&C[(i*84)+72], c23_2);
#else
C[(i*84)+23] += values[119] * B[(i*84)+23];
C[(i*84)+44] += values[120] * B[(i*84)+23];
C[(i*84)+72] += values[121] * B[(i*84)+23];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b24 = _mm256_broadcast_sd(&B[(i*84)+24]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b24 = _mm_loaddup_pd(&B[(i*84)+24]);
#endif
__m128d c24_0 = _mm_load_sd(&C[(i*84)+24]);
__m128d a24_0 = _mm_load_sd(&values[122]);
#if defined(__SSE3__) && defined(__AVX__)
c24_0 = _mm_add_sd(c24_0, _mm_mul_sd(a24_0, _mm256_castpd256_pd128(b24)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c24_0 = _mm_add_sd(c24_0, _mm_mul_sd(a24_0, b24));
#endif
_mm_store_sd(&C[(i*84)+24], c24_0);
__m128d c24_1 = _mm_load_sd(&C[(i*84)+45]);
__m128d a24_1 = _mm_load_sd(&values[123]);
#if defined(__SSE3__) && defined(__AVX__)
c24_1 = _mm_add_sd(c24_1, _mm_mul_sd(a24_1, _mm256_castpd256_pd128(b24)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c24_1 = _mm_add_sd(c24_1, _mm_mul_sd(a24_1, b24));
#endif
_mm_store_sd(&C[(i*84)+45], c24_1);
__m128d c24_2 = _mm_load_sd(&C[(i*84)+73]);
__m128d a24_2 = _mm_load_sd(&values[124]);
#if defined(__SSE3__) && defined(__AVX__)
c24_2 = _mm_add_sd(c24_2, _mm_mul_sd(a24_2, _mm256_castpd256_pd128(b24)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c24_2 = _mm_add_sd(c24_2, _mm_mul_sd(a24_2, b24));
#endif
_mm_store_sd(&C[(i*84)+73], c24_2);
#else
C[(i*84)+24] += values[122] * B[(i*84)+24];
C[(i*84)+45] += values[123] * B[(i*84)+24];
C[(i*84)+73] += values[124] * B[(i*84)+24];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b25 = _mm256_broadcast_sd(&B[(i*84)+25]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b25 = _mm_loaddup_pd(&B[(i*84)+25]);
#endif
__m128d c25_0 = _mm_load_sd(&C[(i*84)+10]);
__m128d a25_0 = _mm_load_sd(&values[125]);
#if defined(__SSE3__) && defined(__AVX__)
c25_0 = _mm_add_sd(c25_0, _mm_mul_sd(a25_0, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c25_0 = _mm_add_sd(c25_0, _mm_mul_sd(a25_0, b25));
#endif
_mm_store_sd(&C[(i*84)+10], c25_0);
__m128d c25_1 = _mm_load_sd(&C[(i*84)+25]);
__m128d a25_1 = _mm_load_sd(&values[126]);
#if defined(__SSE3__) && defined(__AVX__)
c25_1 = _mm_add_sd(c25_1, _mm_mul_sd(a25_1, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c25_1 = _mm_add_sd(c25_1, _mm_mul_sd(a25_1, b25));
#endif
_mm_store_sd(&C[(i*84)+25], c25_1);
__m128d c25_2 = _mm_load_sd(&C[(i*84)+46]);
__m128d a25_2 = _mm_load_sd(&values[127]);
#if defined(__SSE3__) && defined(__AVX__)
c25_2 = _mm_add_sd(c25_2, _mm_mul_sd(a25_2, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c25_2 = _mm_add_sd(c25_2, _mm_mul_sd(a25_2, b25));
#endif
_mm_store_sd(&C[(i*84)+46], c25_2);
__m128d c25_3 = _mm_load_sd(&C[(i*84)+74]);
__m128d a25_3 = _mm_load_sd(&values[128]);
#if defined(__SSE3__) && defined(__AVX__)
c25_3 = _mm_add_sd(c25_3, _mm_mul_sd(a25_3, _mm256_castpd256_pd128(b25)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c25_3 = _mm_add_sd(c25_3, _mm_mul_sd(a25_3, b25));
#endif
_mm_store_sd(&C[(i*84)+74], c25_3);
#else
C[(i*84)+10] += values[125] * B[(i*84)+25];
C[(i*84)+25] += values[126] * B[(i*84)+25];
C[(i*84)+46] += values[127] * B[(i*84)+25];
C[(i*84)+74] += values[128] * B[(i*84)+25];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b26 = _mm256_broadcast_sd(&B[(i*84)+26]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b26 = _mm_loaddup_pd(&B[(i*84)+26]);
#endif
__m128d c26_0 = _mm_load_sd(&C[(i*84)+11]);
__m128d a26_0 = _mm_load_sd(&values[129]);
#if defined(__SSE3__) && defined(__AVX__)
c26_0 = _mm_add_sd(c26_0, _mm_mul_sd(a26_0, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c26_0 = _mm_add_sd(c26_0, _mm_mul_sd(a26_0, b26));
#endif
_mm_store_sd(&C[(i*84)+11], c26_0);
__m128d c26_1 = _mm_load_sd(&C[(i*84)+26]);
__m128d a26_1 = _mm_load_sd(&values[130]);
#if defined(__SSE3__) && defined(__AVX__)
c26_1 = _mm_add_sd(c26_1, _mm_mul_sd(a26_1, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c26_1 = _mm_add_sd(c26_1, _mm_mul_sd(a26_1, b26));
#endif
_mm_store_sd(&C[(i*84)+26], c26_1);
__m128d c26_2 = _mm_load_sd(&C[(i*84)+47]);
__m128d a26_2 = _mm_load_sd(&values[131]);
#if defined(__SSE3__) && defined(__AVX__)
c26_2 = _mm_add_sd(c26_2, _mm_mul_sd(a26_2, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c26_2 = _mm_add_sd(c26_2, _mm_mul_sd(a26_2, b26));
#endif
_mm_store_sd(&C[(i*84)+47], c26_2);
__m128d c26_3 = _mm_load_sd(&C[(i*84)+75]);
__m128d a26_3 = _mm_load_sd(&values[132]);
#if defined(__SSE3__) && defined(__AVX__)
c26_3 = _mm_add_sd(c26_3, _mm_mul_sd(a26_3, _mm256_castpd256_pd128(b26)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c26_3 = _mm_add_sd(c26_3, _mm_mul_sd(a26_3, b26));
#endif
_mm_store_sd(&C[(i*84)+75], c26_3);
#else
C[(i*84)+11] += values[129] * B[(i*84)+26];
C[(i*84)+26] += values[130] * B[(i*84)+26];
C[(i*84)+47] += values[131] * B[(i*84)+26];
C[(i*84)+75] += values[132] * B[(i*84)+26];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b27 = _mm256_broadcast_sd(&B[(i*84)+27]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b27 = _mm_loaddup_pd(&B[(i*84)+27]);
#endif
__m128d c27_0 = _mm_load_sd(&C[(i*84)+12]);
__m128d a27_0 = _mm_load_sd(&values[133]);
#if defined(__SSE3__) && defined(__AVX__)
c27_0 = _mm_add_sd(c27_0, _mm_mul_sd(a27_0, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c27_0 = _mm_add_sd(c27_0, _mm_mul_sd(a27_0, b27));
#endif
_mm_store_sd(&C[(i*84)+12], c27_0);
__m128d c27_1 = _mm_load_sd(&C[(i*84)+27]);
__m128d a27_1 = _mm_load_sd(&values[134]);
#if defined(__SSE3__) && defined(__AVX__)
c27_1 = _mm_add_sd(c27_1, _mm_mul_sd(a27_1, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c27_1 = _mm_add_sd(c27_1, _mm_mul_sd(a27_1, b27));
#endif
_mm_store_sd(&C[(i*84)+27], c27_1);
__m128d c27_2 = _mm_load_sd(&C[(i*84)+48]);
__m128d a27_2 = _mm_load_sd(&values[135]);
#if defined(__SSE3__) && defined(__AVX__)
c27_2 = _mm_add_sd(c27_2, _mm_mul_sd(a27_2, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c27_2 = _mm_add_sd(c27_2, _mm_mul_sd(a27_2, b27));
#endif
_mm_store_sd(&C[(i*84)+48], c27_2);
__m128d c27_3 = _mm_load_sd(&C[(i*84)+76]);
__m128d a27_3 = _mm_load_sd(&values[136]);
#if defined(__SSE3__) && defined(__AVX__)
c27_3 = _mm_add_sd(c27_3, _mm_mul_sd(a27_3, _mm256_castpd256_pd128(b27)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c27_3 = _mm_add_sd(c27_3, _mm_mul_sd(a27_3, b27));
#endif
_mm_store_sd(&C[(i*84)+76], c27_3);
#else
C[(i*84)+12] += values[133] * B[(i*84)+27];
C[(i*84)+27] += values[134] * B[(i*84)+27];
C[(i*84)+48] += values[135] * B[(i*84)+27];
C[(i*84)+76] += values[136] * B[(i*84)+27];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b28 = _mm256_broadcast_sd(&B[(i*84)+28]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b28 = _mm_loaddup_pd(&B[(i*84)+28]);
#endif
__m128d c28_0 = _mm_load_sd(&C[(i*84)+13]);
__m128d a28_0 = _mm_load_sd(&values[137]);
#if defined(__SSE3__) && defined(__AVX__)
c28_0 = _mm_add_sd(c28_0, _mm_mul_sd(a28_0, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c28_0 = _mm_add_sd(c28_0, _mm_mul_sd(a28_0, b28));
#endif
_mm_store_sd(&C[(i*84)+13], c28_0);
__m128d c28_1 = _mm_load_sd(&C[(i*84)+28]);
__m128d a28_1 = _mm_load_sd(&values[138]);
#if defined(__SSE3__) && defined(__AVX__)
c28_1 = _mm_add_sd(c28_1, _mm_mul_sd(a28_1, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c28_1 = _mm_add_sd(c28_1, _mm_mul_sd(a28_1, b28));
#endif
_mm_store_sd(&C[(i*84)+28], c28_1);
__m128d c28_2 = _mm_load_sd(&C[(i*84)+49]);
__m128d a28_2 = _mm_load_sd(&values[139]);
#if defined(__SSE3__) && defined(__AVX__)
c28_2 = _mm_add_sd(c28_2, _mm_mul_sd(a28_2, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c28_2 = _mm_add_sd(c28_2, _mm_mul_sd(a28_2, b28));
#endif
_mm_store_sd(&C[(i*84)+49], c28_2);
__m128d c28_3 = _mm_load_sd(&C[(i*84)+77]);
__m128d a28_3 = _mm_load_sd(&values[140]);
#if defined(__SSE3__) && defined(__AVX__)
c28_3 = _mm_add_sd(c28_3, _mm_mul_sd(a28_3, _mm256_castpd256_pd128(b28)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c28_3 = _mm_add_sd(c28_3, _mm_mul_sd(a28_3, b28));
#endif
_mm_store_sd(&C[(i*84)+77], c28_3);
#else
C[(i*84)+13] += values[137] * B[(i*84)+28];
C[(i*84)+28] += values[138] * B[(i*84)+28];
C[(i*84)+49] += values[139] * B[(i*84)+28];
C[(i*84)+77] += values[140] * B[(i*84)+28];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b29 = _mm256_broadcast_sd(&B[(i*84)+29]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b29 = _mm_loaddup_pd(&B[(i*84)+29]);
#endif
__m128d c29_0 = _mm_load_sd(&C[(i*84)+4]);
__m128d a29_0 = _mm_load_sd(&values[141]);
#if defined(__SSE3__) && defined(__AVX__)
c29_0 = _mm_add_sd(c29_0, _mm_mul_sd(a29_0, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_0 = _mm_add_sd(c29_0, _mm_mul_sd(a29_0, b29));
#endif
_mm_store_sd(&C[(i*84)+4], c29_0);
__m128d c29_1 = _mm_load_sd(&C[(i*84)+14]);
__m128d a29_1 = _mm_load_sd(&values[142]);
#if defined(__SSE3__) && defined(__AVX__)
c29_1 = _mm_add_sd(c29_1, _mm_mul_sd(a29_1, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_1 = _mm_add_sd(c29_1, _mm_mul_sd(a29_1, b29));
#endif
_mm_store_sd(&C[(i*84)+14], c29_1);
__m128d c29_2 = _mm_load_sd(&C[(i*84)+29]);
__m128d a29_2 = _mm_load_sd(&values[143]);
#if defined(__SSE3__) && defined(__AVX__)
c29_2 = _mm_add_sd(c29_2, _mm_mul_sd(a29_2, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_2 = _mm_add_sd(c29_2, _mm_mul_sd(a29_2, b29));
#endif
_mm_store_sd(&C[(i*84)+29], c29_2);
__m128d c29_3 = _mm_load_sd(&C[(i*84)+50]);
__m128d a29_3 = _mm_load_sd(&values[144]);
#if defined(__SSE3__) && defined(__AVX__)
c29_3 = _mm_add_sd(c29_3, _mm_mul_sd(a29_3, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_3 = _mm_add_sd(c29_3, _mm_mul_sd(a29_3, b29));
#endif
_mm_store_sd(&C[(i*84)+50], c29_3);
__m128d c29_4 = _mm_load_sd(&C[(i*84)+78]);
__m128d a29_4 = _mm_load_sd(&values[145]);
#if defined(__SSE3__) && defined(__AVX__)
c29_4 = _mm_add_sd(c29_4, _mm_mul_sd(a29_4, _mm256_castpd256_pd128(b29)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c29_4 = _mm_add_sd(c29_4, _mm_mul_sd(a29_4, b29));
#endif
_mm_store_sd(&C[(i*84)+78], c29_4);
#else
C[(i*84)+4] += values[141] * B[(i*84)+29];
C[(i*84)+14] += values[142] * B[(i*84)+29];
C[(i*84)+29] += values[143] * B[(i*84)+29];
C[(i*84)+50] += values[144] * B[(i*84)+29];
C[(i*84)+78] += values[145] * B[(i*84)+29];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b30 = _mm256_broadcast_sd(&B[(i*84)+30]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b30 = _mm_loaddup_pd(&B[(i*84)+30]);
#endif
__m128d c30_0 = _mm_load_sd(&C[(i*84)+5]);
__m128d a30_0 = _mm_load_sd(&values[146]);
#if defined(__SSE3__) && defined(__AVX__)
c30_0 = _mm_add_sd(c30_0, _mm_mul_sd(a30_0, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_0 = _mm_add_sd(c30_0, _mm_mul_sd(a30_0, b30));
#endif
_mm_store_sd(&C[(i*84)+5], c30_0);
__m128d c30_1 = _mm_load_sd(&C[(i*84)+15]);
__m128d a30_1 = _mm_load_sd(&values[147]);
#if defined(__SSE3__) && defined(__AVX__)
c30_1 = _mm_add_sd(c30_1, _mm_mul_sd(a30_1, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_1 = _mm_add_sd(c30_1, _mm_mul_sd(a30_1, b30));
#endif
_mm_store_sd(&C[(i*84)+15], c30_1);
__m128d c30_2 = _mm_load_sd(&C[(i*84)+30]);
__m128d a30_2 = _mm_load_sd(&values[148]);
#if defined(__SSE3__) && defined(__AVX__)
c30_2 = _mm_add_sd(c30_2, _mm_mul_sd(a30_2, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_2 = _mm_add_sd(c30_2, _mm_mul_sd(a30_2, b30));
#endif
_mm_store_sd(&C[(i*84)+30], c30_2);
__m128d c30_3 = _mm_load_sd(&C[(i*84)+51]);
__m128d a30_3 = _mm_load_sd(&values[149]);
#if defined(__SSE3__) && defined(__AVX__)
c30_3 = _mm_add_sd(c30_3, _mm_mul_sd(a30_3, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_3 = _mm_add_sd(c30_3, _mm_mul_sd(a30_3, b30));
#endif
_mm_store_sd(&C[(i*84)+51], c30_3);
__m128d c30_4 = _mm_load_sd(&C[(i*84)+79]);
__m128d a30_4 = _mm_load_sd(&values[150]);
#if defined(__SSE3__) && defined(__AVX__)
c30_4 = _mm_add_sd(c30_4, _mm_mul_sd(a30_4, _mm256_castpd256_pd128(b30)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c30_4 = _mm_add_sd(c30_4, _mm_mul_sd(a30_4, b30));
#endif
_mm_store_sd(&C[(i*84)+79], c30_4);
#else
C[(i*84)+5] += values[146] * B[(i*84)+30];
C[(i*84)+15] += values[147] * B[(i*84)+30];
C[(i*84)+30] += values[148] * B[(i*84)+30];
C[(i*84)+51] += values[149] * B[(i*84)+30];
C[(i*84)+79] += values[150] * B[(i*84)+30];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b31 = _mm256_broadcast_sd(&B[(i*84)+31]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b31 = _mm_loaddup_pd(&B[(i*84)+31]);
#endif
__m128d c31_0 = _mm_load_sd(&C[(i*84)+6]);
__m128d a31_0 = _mm_load_sd(&values[151]);
#if defined(__SSE3__) && defined(__AVX__)
c31_0 = _mm_add_sd(c31_0, _mm_mul_sd(a31_0, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_0 = _mm_add_sd(c31_0, _mm_mul_sd(a31_0, b31));
#endif
_mm_store_sd(&C[(i*84)+6], c31_0);
__m128d c31_1 = _mm_load_sd(&C[(i*84)+16]);
__m128d a31_1 = _mm_load_sd(&values[152]);
#if defined(__SSE3__) && defined(__AVX__)
c31_1 = _mm_add_sd(c31_1, _mm_mul_sd(a31_1, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_1 = _mm_add_sd(c31_1, _mm_mul_sd(a31_1, b31));
#endif
_mm_store_sd(&C[(i*84)+16], c31_1);
__m128d c31_2 = _mm_load_sd(&C[(i*84)+31]);
__m128d a31_2 = _mm_load_sd(&values[153]);
#if defined(__SSE3__) && defined(__AVX__)
c31_2 = _mm_add_sd(c31_2, _mm_mul_sd(a31_2, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_2 = _mm_add_sd(c31_2, _mm_mul_sd(a31_2, b31));
#endif
_mm_store_sd(&C[(i*84)+31], c31_2);
__m128d c31_3 = _mm_load_sd(&C[(i*84)+52]);
__m128d a31_3 = _mm_load_sd(&values[154]);
#if defined(__SSE3__) && defined(__AVX__)
c31_3 = _mm_add_sd(c31_3, _mm_mul_sd(a31_3, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_3 = _mm_add_sd(c31_3, _mm_mul_sd(a31_3, b31));
#endif
_mm_store_sd(&C[(i*84)+52], c31_3);
__m128d c31_4 = _mm_load_sd(&C[(i*84)+80]);
__m128d a31_4 = _mm_load_sd(&values[155]);
#if defined(__SSE3__) && defined(__AVX__)
c31_4 = _mm_add_sd(c31_4, _mm_mul_sd(a31_4, _mm256_castpd256_pd128(b31)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c31_4 = _mm_add_sd(c31_4, _mm_mul_sd(a31_4, b31));
#endif
_mm_store_sd(&C[(i*84)+80], c31_4);
#else
C[(i*84)+6] += values[151] * B[(i*84)+31];
C[(i*84)+16] += values[152] * B[(i*84)+31];
C[(i*84)+31] += values[153] * B[(i*84)+31];
C[(i*84)+52] += values[154] * B[(i*84)+31];
C[(i*84)+80] += values[155] * B[(i*84)+31];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b32 = _mm256_broadcast_sd(&B[(i*84)+32]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b32 = _mm_loaddup_pd(&B[(i*84)+32]);
#endif
__m128d c32_0 = _mm_load_sd(&C[(i*84)+1]);
__m128d a32_0 = _mm_load_sd(&values[156]);
#if defined(__SSE3__) && defined(__AVX__)
c32_0 = _mm_add_sd(c32_0, _mm_mul_sd(a32_0, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_0 = _mm_add_sd(c32_0, _mm_mul_sd(a32_0, b32));
#endif
_mm_store_sd(&C[(i*84)+1], c32_0);
__m128d c32_1 = _mm_load_sd(&C[(i*84)+7]);
__m128d a32_1 = _mm_load_sd(&values[157]);
#if defined(__SSE3__) && defined(__AVX__)
c32_1 = _mm_add_sd(c32_1, _mm_mul_sd(a32_1, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_1 = _mm_add_sd(c32_1, _mm_mul_sd(a32_1, b32));
#endif
_mm_store_sd(&C[(i*84)+7], c32_1);
__m128d c32_2 = _mm_load_sd(&C[(i*84)+17]);
__m128d a32_2 = _mm_load_sd(&values[158]);
#if defined(__SSE3__) && defined(__AVX__)
c32_2 = _mm_add_sd(c32_2, _mm_mul_sd(a32_2, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_2 = _mm_add_sd(c32_2, _mm_mul_sd(a32_2, b32));
#endif
_mm_store_sd(&C[(i*84)+17], c32_2);
__m128d c32_3 = _mm_load_sd(&C[(i*84)+32]);
__m128d a32_3 = _mm_load_sd(&values[159]);
#if defined(__SSE3__) && defined(__AVX__)
c32_3 = _mm_add_sd(c32_3, _mm_mul_sd(a32_3, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_3 = _mm_add_sd(c32_3, _mm_mul_sd(a32_3, b32));
#endif
_mm_store_sd(&C[(i*84)+32], c32_3);
__m128d c32_4 = _mm_load_sd(&C[(i*84)+53]);
__m128d a32_4 = _mm_load_sd(&values[160]);
#if defined(__SSE3__) && defined(__AVX__)
c32_4 = _mm_add_sd(c32_4, _mm_mul_sd(a32_4, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_4 = _mm_add_sd(c32_4, _mm_mul_sd(a32_4, b32));
#endif
_mm_store_sd(&C[(i*84)+53], c32_4);
__m128d c32_5 = _mm_load_sd(&C[(i*84)+81]);
__m128d a32_5 = _mm_load_sd(&values[161]);
#if defined(__SSE3__) && defined(__AVX__)
c32_5 = _mm_add_sd(c32_5, _mm_mul_sd(a32_5, _mm256_castpd256_pd128(b32)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c32_5 = _mm_add_sd(c32_5, _mm_mul_sd(a32_5, b32));
#endif
_mm_store_sd(&C[(i*84)+81], c32_5);
#else
C[(i*84)+1] += values[156] * B[(i*84)+32];
C[(i*84)+7] += values[157] * B[(i*84)+32];
C[(i*84)+17] += values[158] * B[(i*84)+32];
C[(i*84)+32] += values[159] * B[(i*84)+32];
C[(i*84)+53] += values[160] * B[(i*84)+32];
C[(i*84)+81] += values[161] * B[(i*84)+32];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b33 = _mm256_broadcast_sd(&B[(i*84)+33]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b33 = _mm_loaddup_pd(&B[(i*84)+33]);
#endif
__m128d c33_0 = _mm_load_sd(&C[(i*84)+2]);
__m128d a33_0 = _mm_load_sd(&values[162]);
#if defined(__SSE3__) && defined(__AVX__)
c33_0 = _mm_add_sd(c33_0, _mm_mul_sd(a33_0, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_0 = _mm_add_sd(c33_0, _mm_mul_sd(a33_0, b33));
#endif
_mm_store_sd(&C[(i*84)+2], c33_0);
__m128d c33_1 = _mm_load_sd(&C[(i*84)+8]);
__m128d a33_1 = _mm_load_sd(&values[163]);
#if defined(__SSE3__) && defined(__AVX__)
c33_1 = _mm_add_sd(c33_1, _mm_mul_sd(a33_1, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_1 = _mm_add_sd(c33_1, _mm_mul_sd(a33_1, b33));
#endif
_mm_store_sd(&C[(i*84)+8], c33_1);
__m128d c33_2 = _mm_load_sd(&C[(i*84)+18]);
__m128d a33_2 = _mm_load_sd(&values[164]);
#if defined(__SSE3__) && defined(__AVX__)
c33_2 = _mm_add_sd(c33_2, _mm_mul_sd(a33_2, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_2 = _mm_add_sd(c33_2, _mm_mul_sd(a33_2, b33));
#endif
_mm_store_sd(&C[(i*84)+18], c33_2);
__m128d c33_3 = _mm_load_sd(&C[(i*84)+33]);
__m128d a33_3 = _mm_load_sd(&values[165]);
#if defined(__SSE3__) && defined(__AVX__)
c33_3 = _mm_add_sd(c33_3, _mm_mul_sd(a33_3, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_3 = _mm_add_sd(c33_3, _mm_mul_sd(a33_3, b33));
#endif
_mm_store_sd(&C[(i*84)+33], c33_3);
__m128d c33_4 = _mm_load_sd(&C[(i*84)+54]);
__m128d a33_4 = _mm_load_sd(&values[166]);
#if defined(__SSE3__) && defined(__AVX__)
c33_4 = _mm_add_sd(c33_4, _mm_mul_sd(a33_4, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_4 = _mm_add_sd(c33_4, _mm_mul_sd(a33_4, b33));
#endif
_mm_store_sd(&C[(i*84)+54], c33_4);
__m128d c33_5 = _mm_load_sd(&C[(i*84)+82]);
__m128d a33_5 = _mm_load_sd(&values[167]);
#if defined(__SSE3__) && defined(__AVX__)
c33_5 = _mm_add_sd(c33_5, _mm_mul_sd(a33_5, _mm256_castpd256_pd128(b33)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c33_5 = _mm_add_sd(c33_5, _mm_mul_sd(a33_5, b33));
#endif
_mm_store_sd(&C[(i*84)+82], c33_5);
#else
C[(i*84)+2] += values[162] * B[(i*84)+33];
C[(i*84)+8] += values[163] * B[(i*84)+33];
C[(i*84)+18] += values[164] * B[(i*84)+33];
C[(i*84)+33] += values[165] * B[(i*84)+33];
C[(i*84)+54] += values[166] * B[(i*84)+33];
C[(i*84)+82] += values[167] * B[(i*84)+33];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b34 = _mm256_broadcast_sd(&B[(i*84)+34]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b34 = _mm_loaddup_pd(&B[(i*84)+34]);
#endif
__m128d c34_0 = _mm_load_sd(&C[(i*84)+0]);
__m128d a34_0 = _mm_load_sd(&values[168]);
#if defined(__SSE3__) && defined(__AVX__)
c34_0 = _mm_add_sd(c34_0, _mm_mul_sd(a34_0, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_0 = _mm_add_sd(c34_0, _mm_mul_sd(a34_0, b34));
#endif
_mm_store_sd(&C[(i*84)+0], c34_0);
__m128d c34_1 = _mm_load_sd(&C[(i*84)+3]);
__m128d a34_1 = _mm_load_sd(&values[169]);
#if defined(__SSE3__) && defined(__AVX__)
c34_1 = _mm_add_sd(c34_1, _mm_mul_sd(a34_1, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_1 = _mm_add_sd(c34_1, _mm_mul_sd(a34_1, b34));
#endif
_mm_store_sd(&C[(i*84)+3], c34_1);
__m128d c34_2 = _mm_load_sd(&C[(i*84)+9]);
__m128d a34_2 = _mm_load_sd(&values[170]);
#if defined(__SSE3__) && defined(__AVX__)
c34_2 = _mm_add_sd(c34_2, _mm_mul_sd(a34_2, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_2 = _mm_add_sd(c34_2, _mm_mul_sd(a34_2, b34));
#endif
_mm_store_sd(&C[(i*84)+9], c34_2);
__m128d c34_3 = _mm_load_sd(&C[(i*84)+19]);
__m128d a34_3 = _mm_load_sd(&values[171]);
#if defined(__SSE3__) && defined(__AVX__)
c34_3 = _mm_add_sd(c34_3, _mm_mul_sd(a34_3, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_3 = _mm_add_sd(c34_3, _mm_mul_sd(a34_3, b34));
#endif
_mm_store_sd(&C[(i*84)+19], c34_3);
__m128d c34_4 = _mm_load_sd(&C[(i*84)+34]);
__m128d a34_4 = _mm_load_sd(&values[172]);
#if defined(__SSE3__) && defined(__AVX__)
c34_4 = _mm_add_sd(c34_4, _mm_mul_sd(a34_4, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_4 = _mm_add_sd(c34_4, _mm_mul_sd(a34_4, b34));
#endif
_mm_store_sd(&C[(i*84)+34], c34_4);
__m128d c34_5 = _mm_load_sd(&C[(i*84)+55]);
__m128d a34_5 = _mm_load_sd(&values[173]);
#if defined(__SSE3__) && defined(__AVX__)
c34_5 = _mm_add_sd(c34_5, _mm_mul_sd(a34_5, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_5 = _mm_add_sd(c34_5, _mm_mul_sd(a34_5, b34));
#endif
_mm_store_sd(&C[(i*84)+55], c34_5);
__m128d c34_6 = _mm_load_sd(&C[(i*84)+83]);
__m128d a34_6 = _mm_load_sd(&values[174]);
#if defined(__SSE3__) && defined(__AVX__)
c34_6 = _mm_add_sd(c34_6, _mm_mul_sd(a34_6, _mm256_castpd256_pd128(b34)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c34_6 = _mm_add_sd(c34_6, _mm_mul_sd(a34_6, b34));
#endif
_mm_store_sd(&C[(i*84)+83], c34_6);
#else
C[(i*84)+0] += values[168] * B[(i*84)+34];
C[(i*84)+3] += values[169] * B[(i*84)+34];
C[(i*84)+9] += values[170] * B[(i*84)+34];
C[(i*84)+19] += values[171] * B[(i*84)+34];
C[(i*84)+34] += values[172] * B[(i*84)+34];
C[(i*84)+55] += values[173] * B[(i*84)+34];
C[(i*84)+83] += values[174] * B[(i*84)+34];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b35 = _mm256_broadcast_sd(&B[(i*84)+35]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b35 = _mm_loaddup_pd(&B[(i*84)+35]);
#endif
__m128d c35_0 = _mm_load_sd(&C[(i*84)+35]);
__m128d a35_0 = _mm_load_sd(&values[175]);
#if defined(__SSE3__) && defined(__AVX__)
c35_0 = _mm_add_sd(c35_0, _mm_mul_sd(a35_0, _mm256_castpd256_pd128(b35)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c35_0 = _mm_add_sd(c35_0, _mm_mul_sd(a35_0, b35));
#endif
_mm_store_sd(&C[(i*84)+35], c35_0);
__m128d c35_1 = _mm_load_sd(&C[(i*84)+63]);
__m128d a35_1 = _mm_load_sd(&values[176]);
#if defined(__SSE3__) && defined(__AVX__)
c35_1 = _mm_add_sd(c35_1, _mm_mul_sd(a35_1, _mm256_castpd256_pd128(b35)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c35_1 = _mm_add_sd(c35_1, _mm_mul_sd(a35_1, b35));
#endif
_mm_store_sd(&C[(i*84)+63], c35_1);
#else
C[(i*84)+35] += values[175] * B[(i*84)+35];
C[(i*84)+63] += values[176] * B[(i*84)+35];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b36 = _mm256_broadcast_sd(&B[(i*84)+36]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b36 = _mm_loaddup_pd(&B[(i*84)+36]);
#endif
__m128d c36_0 = _mm_load_sd(&C[(i*84)+36]);
__m128d a36_0 = _mm_load_sd(&values[177]);
#if defined(__SSE3__) && defined(__AVX__)
c36_0 = _mm_add_sd(c36_0, _mm_mul_sd(a36_0, _mm256_castpd256_pd128(b36)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c36_0 = _mm_add_sd(c36_0, _mm_mul_sd(a36_0, b36));
#endif
_mm_store_sd(&C[(i*84)+36], c36_0);
__m128d c36_1 = _mm_load_sd(&C[(i*84)+64]);
__m128d a36_1 = _mm_load_sd(&values[178]);
#if defined(__SSE3__) && defined(__AVX__)
c36_1 = _mm_add_sd(c36_1, _mm_mul_sd(a36_1, _mm256_castpd256_pd128(b36)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c36_1 = _mm_add_sd(c36_1, _mm_mul_sd(a36_1, b36));
#endif
_mm_store_sd(&C[(i*84)+64], c36_1);
#else
C[(i*84)+36] += values[177] * B[(i*84)+36];
C[(i*84)+64] += values[178] * B[(i*84)+36];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b37 = _mm256_broadcast_sd(&B[(i*84)+37]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b37 = _mm_loaddup_pd(&B[(i*84)+37]);
#endif
__m128d c37_0 = _mm_load_sd(&C[(i*84)+37]);
__m128d a37_0 = _mm_load_sd(&values[179]);
#if defined(__SSE3__) && defined(__AVX__)
c37_0 = _mm_add_sd(c37_0, _mm_mul_sd(a37_0, _mm256_castpd256_pd128(b37)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c37_0 = _mm_add_sd(c37_0, _mm_mul_sd(a37_0, b37));
#endif
_mm_store_sd(&C[(i*84)+37], c37_0);
__m128d c37_1 = _mm_load_sd(&C[(i*84)+65]);
__m128d a37_1 = _mm_load_sd(&values[180]);
#if defined(__SSE3__) && defined(__AVX__)
c37_1 = _mm_add_sd(c37_1, _mm_mul_sd(a37_1, _mm256_castpd256_pd128(b37)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c37_1 = _mm_add_sd(c37_1, _mm_mul_sd(a37_1, b37));
#endif
_mm_store_sd(&C[(i*84)+65], c37_1);
#else
C[(i*84)+37] += values[179] * B[(i*84)+37];
C[(i*84)+65] += values[180] * B[(i*84)+37];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b38 = _mm256_broadcast_sd(&B[(i*84)+38]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b38 = _mm_loaddup_pd(&B[(i*84)+38]);
#endif
__m128d c38_0 = _mm_load_sd(&C[(i*84)+38]);
__m128d a38_0 = _mm_load_sd(&values[181]);
#if defined(__SSE3__) && defined(__AVX__)
c38_0 = _mm_add_sd(c38_0, _mm_mul_sd(a38_0, _mm256_castpd256_pd128(b38)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c38_0 = _mm_add_sd(c38_0, _mm_mul_sd(a38_0, b38));
#endif
_mm_store_sd(&C[(i*84)+38], c38_0);
__m128d c38_1 = _mm_load_sd(&C[(i*84)+66]);
__m128d a38_1 = _mm_load_sd(&values[182]);
#if defined(__SSE3__) && defined(__AVX__)
c38_1 = _mm_add_sd(c38_1, _mm_mul_sd(a38_1, _mm256_castpd256_pd128(b38)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c38_1 = _mm_add_sd(c38_1, _mm_mul_sd(a38_1, b38));
#endif
_mm_store_sd(&C[(i*84)+66], c38_1);
#else
C[(i*84)+38] += values[181] * B[(i*84)+38];
C[(i*84)+66] += values[182] * B[(i*84)+38];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b39 = _mm256_broadcast_sd(&B[(i*84)+39]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b39 = _mm_loaddup_pd(&B[(i*84)+39]);
#endif
__m128d c39_0 = _mm_load_sd(&C[(i*84)+39]);
__m128d a39_0 = _mm_load_sd(&values[183]);
#if defined(__SSE3__) && defined(__AVX__)
c39_0 = _mm_add_sd(c39_0, _mm_mul_sd(a39_0, _mm256_castpd256_pd128(b39)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c39_0 = _mm_add_sd(c39_0, _mm_mul_sd(a39_0, b39));
#endif
_mm_store_sd(&C[(i*84)+39], c39_0);
__m128d c39_1 = _mm_load_sd(&C[(i*84)+67]);
__m128d a39_1 = _mm_load_sd(&values[184]);
#if defined(__SSE3__) && defined(__AVX__)
c39_1 = _mm_add_sd(c39_1, _mm_mul_sd(a39_1, _mm256_castpd256_pd128(b39)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c39_1 = _mm_add_sd(c39_1, _mm_mul_sd(a39_1, b39));
#endif
_mm_store_sd(&C[(i*84)+67], c39_1);
#else
C[(i*84)+39] += values[183] * B[(i*84)+39];
C[(i*84)+67] += values[184] * B[(i*84)+39];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b40 = _mm256_broadcast_sd(&B[(i*84)+40]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b40 = _mm_loaddup_pd(&B[(i*84)+40]);
#endif
__m128d c40_0 = _mm_load_sd(&C[(i*84)+40]);
__m128d a40_0 = _mm_load_sd(&values[185]);
#if defined(__SSE3__) && defined(__AVX__)
c40_0 = _mm_add_sd(c40_0, _mm_mul_sd(a40_0, _mm256_castpd256_pd128(b40)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c40_0 = _mm_add_sd(c40_0, _mm_mul_sd(a40_0, b40));
#endif
_mm_store_sd(&C[(i*84)+40], c40_0);
__m128d c40_1 = _mm_load_sd(&C[(i*84)+68]);
__m128d a40_1 = _mm_load_sd(&values[186]);
#if defined(__SSE3__) && defined(__AVX__)
c40_1 = _mm_add_sd(c40_1, _mm_mul_sd(a40_1, _mm256_castpd256_pd128(b40)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c40_1 = _mm_add_sd(c40_1, _mm_mul_sd(a40_1, b40));
#endif
_mm_store_sd(&C[(i*84)+68], c40_1);
#else
C[(i*84)+40] += values[185] * B[(i*84)+40];
C[(i*84)+68] += values[186] * B[(i*84)+40];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b41 = _mm256_broadcast_sd(&B[(i*84)+41]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b41 = _mm_loaddup_pd(&B[(i*84)+41]);
#endif
__m128d c41_0 = _mm_load_sd(&C[(i*84)+20]);
__m128d a41_0 = _mm_load_sd(&values[187]);
#if defined(__SSE3__) && defined(__AVX__)
c41_0 = _mm_add_sd(c41_0, _mm_mul_sd(a41_0, _mm256_castpd256_pd128(b41)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c41_0 = _mm_add_sd(c41_0, _mm_mul_sd(a41_0, b41));
#endif
_mm_store_sd(&C[(i*84)+20], c41_0);
__m128d c41_1 = _mm_load_sd(&C[(i*84)+41]);
__m128d a41_1 = _mm_load_sd(&values[188]);
#if defined(__SSE3__) && defined(__AVX__)
c41_1 = _mm_add_sd(c41_1, _mm_mul_sd(a41_1, _mm256_castpd256_pd128(b41)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c41_1 = _mm_add_sd(c41_1, _mm_mul_sd(a41_1, b41));
#endif
_mm_store_sd(&C[(i*84)+41], c41_1);
__m128d c41_2 = _mm_load_sd(&C[(i*84)+69]);
__m128d a41_2 = _mm_load_sd(&values[189]);
#if defined(__SSE3__) && defined(__AVX__)
c41_2 = _mm_add_sd(c41_2, _mm_mul_sd(a41_2, _mm256_castpd256_pd128(b41)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c41_2 = _mm_add_sd(c41_2, _mm_mul_sd(a41_2, b41));
#endif
_mm_store_sd(&C[(i*84)+69], c41_2);
#else
C[(i*84)+20] += values[187] * B[(i*84)+41];
C[(i*84)+41] += values[188] * B[(i*84)+41];
C[(i*84)+69] += values[189] * B[(i*84)+41];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b42 = _mm256_broadcast_sd(&B[(i*84)+42]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b42 = _mm_loaddup_pd(&B[(i*84)+42]);
#endif
__m128d c42_0 = _mm_load_sd(&C[(i*84)+21]);
__m128d a42_0 = _mm_load_sd(&values[190]);
#if defined(__SSE3__) && defined(__AVX__)
c42_0 = _mm_add_sd(c42_0, _mm_mul_sd(a42_0, _mm256_castpd256_pd128(b42)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c42_0 = _mm_add_sd(c42_0, _mm_mul_sd(a42_0, b42));
#endif
_mm_store_sd(&C[(i*84)+21], c42_0);
__m128d c42_1 = _mm_load_sd(&C[(i*84)+42]);
__m128d a42_1 = _mm_load_sd(&values[191]);
#if defined(__SSE3__) && defined(__AVX__)
c42_1 = _mm_add_sd(c42_1, _mm_mul_sd(a42_1, _mm256_castpd256_pd128(b42)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c42_1 = _mm_add_sd(c42_1, _mm_mul_sd(a42_1, b42));
#endif
_mm_store_sd(&C[(i*84)+42], c42_1);
__m128d c42_2 = _mm_load_sd(&C[(i*84)+70]);
__m128d a42_2 = _mm_load_sd(&values[192]);
#if defined(__SSE3__) && defined(__AVX__)
c42_2 = _mm_add_sd(c42_2, _mm_mul_sd(a42_2, _mm256_castpd256_pd128(b42)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c42_2 = _mm_add_sd(c42_2, _mm_mul_sd(a42_2, b42));
#endif
_mm_store_sd(&C[(i*84)+70], c42_2);
#else
C[(i*84)+21] += values[190] * B[(i*84)+42];
C[(i*84)+42] += values[191] * B[(i*84)+42];
C[(i*84)+70] += values[192] * B[(i*84)+42];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b43 = _mm256_broadcast_sd(&B[(i*84)+43]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b43 = _mm_loaddup_pd(&B[(i*84)+43]);
#endif
__m128d c43_0 = _mm_load_sd(&C[(i*84)+22]);
__m128d a43_0 = _mm_load_sd(&values[193]);
#if defined(__SSE3__) && defined(__AVX__)
c43_0 = _mm_add_sd(c43_0, _mm_mul_sd(a43_0, _mm256_castpd256_pd128(b43)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c43_0 = _mm_add_sd(c43_0, _mm_mul_sd(a43_0, b43));
#endif
_mm_store_sd(&C[(i*84)+22], c43_0);
__m128d c43_1 = _mm_load_sd(&C[(i*84)+43]);
__m128d a43_1 = _mm_load_sd(&values[194]);
#if defined(__SSE3__) && defined(__AVX__)
c43_1 = _mm_add_sd(c43_1, _mm_mul_sd(a43_1, _mm256_castpd256_pd128(b43)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c43_1 = _mm_add_sd(c43_1, _mm_mul_sd(a43_1, b43));
#endif
_mm_store_sd(&C[(i*84)+43], c43_1);
__m128d c43_2 = _mm_load_sd(&C[(i*84)+71]);
__m128d a43_2 = _mm_load_sd(&values[195]);
#if defined(__SSE3__) && defined(__AVX__)
c43_2 = _mm_add_sd(c43_2, _mm_mul_sd(a43_2, _mm256_castpd256_pd128(b43)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c43_2 = _mm_add_sd(c43_2, _mm_mul_sd(a43_2, b43));
#endif
_mm_store_sd(&C[(i*84)+71], c43_2);
#else
C[(i*84)+22] += values[193] * B[(i*84)+43];
C[(i*84)+43] += values[194] * B[(i*84)+43];
C[(i*84)+71] += values[195] * B[(i*84)+43];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b44 = _mm256_broadcast_sd(&B[(i*84)+44]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b44 = _mm_loaddup_pd(&B[(i*84)+44]);
#endif
__m128d c44_0 = _mm_load_sd(&C[(i*84)+23]);
__m128d a44_0 = _mm_load_sd(&values[196]);
#if defined(__SSE3__) && defined(__AVX__)
c44_0 = _mm_add_sd(c44_0, _mm_mul_sd(a44_0, _mm256_castpd256_pd128(b44)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c44_0 = _mm_add_sd(c44_0, _mm_mul_sd(a44_0, b44));
#endif
_mm_store_sd(&C[(i*84)+23], c44_0);
__m128d c44_1 = _mm_load_sd(&C[(i*84)+44]);
__m128d a44_1 = _mm_load_sd(&values[197]);
#if defined(__SSE3__) && defined(__AVX__)
c44_1 = _mm_add_sd(c44_1, _mm_mul_sd(a44_1, _mm256_castpd256_pd128(b44)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c44_1 = _mm_add_sd(c44_1, _mm_mul_sd(a44_1, b44));
#endif
_mm_store_sd(&C[(i*84)+44], c44_1);
__m128d c44_2 = _mm_load_sd(&C[(i*84)+72]);
__m128d a44_2 = _mm_load_sd(&values[198]);
#if defined(__SSE3__) && defined(__AVX__)
c44_2 = _mm_add_sd(c44_2, _mm_mul_sd(a44_2, _mm256_castpd256_pd128(b44)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c44_2 = _mm_add_sd(c44_2, _mm_mul_sd(a44_2, b44));
#endif
_mm_store_sd(&C[(i*84)+72], c44_2);
#else
C[(i*84)+23] += values[196] * B[(i*84)+44];
C[(i*84)+44] += values[197] * B[(i*84)+44];
C[(i*84)+72] += values[198] * B[(i*84)+44];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b45 = _mm256_broadcast_sd(&B[(i*84)+45]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b45 = _mm_loaddup_pd(&B[(i*84)+45]);
#endif
__m128d c45_0 = _mm_load_sd(&C[(i*84)+24]);
__m128d a45_0 = _mm_load_sd(&values[199]);
#if defined(__SSE3__) && defined(__AVX__)
c45_0 = _mm_add_sd(c45_0, _mm_mul_sd(a45_0, _mm256_castpd256_pd128(b45)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c45_0 = _mm_add_sd(c45_0, _mm_mul_sd(a45_0, b45));
#endif
_mm_store_sd(&C[(i*84)+24], c45_0);
__m128d c45_1 = _mm_load_sd(&C[(i*84)+45]);
__m128d a45_1 = _mm_load_sd(&values[200]);
#if defined(__SSE3__) && defined(__AVX__)
c45_1 = _mm_add_sd(c45_1, _mm_mul_sd(a45_1, _mm256_castpd256_pd128(b45)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c45_1 = _mm_add_sd(c45_1, _mm_mul_sd(a45_1, b45));
#endif
_mm_store_sd(&C[(i*84)+45], c45_1);
__m128d c45_2 = _mm_load_sd(&C[(i*84)+73]);
__m128d a45_2 = _mm_load_sd(&values[201]);
#if defined(__SSE3__) && defined(__AVX__)
c45_2 = _mm_add_sd(c45_2, _mm_mul_sd(a45_2, _mm256_castpd256_pd128(b45)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c45_2 = _mm_add_sd(c45_2, _mm_mul_sd(a45_2, b45));
#endif
_mm_store_sd(&C[(i*84)+73], c45_2);
#else
C[(i*84)+24] += values[199] * B[(i*84)+45];
C[(i*84)+45] += values[200] * B[(i*84)+45];
C[(i*84)+73] += values[201] * B[(i*84)+45];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b46 = _mm256_broadcast_sd(&B[(i*84)+46]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b46 = _mm_loaddup_pd(&B[(i*84)+46]);
#endif
__m128d c46_0 = _mm_load_sd(&C[(i*84)+10]);
__m128d a46_0 = _mm_load_sd(&values[202]);
#if defined(__SSE3__) && defined(__AVX__)
c46_0 = _mm_add_sd(c46_0, _mm_mul_sd(a46_0, _mm256_castpd256_pd128(b46)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c46_0 = _mm_add_sd(c46_0, _mm_mul_sd(a46_0, b46));
#endif
_mm_store_sd(&C[(i*84)+10], c46_0);
__m128d c46_1 = _mm_load_sd(&C[(i*84)+25]);
__m128d a46_1 = _mm_load_sd(&values[203]);
#if defined(__SSE3__) && defined(__AVX__)
c46_1 = _mm_add_sd(c46_1, _mm_mul_sd(a46_1, _mm256_castpd256_pd128(b46)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c46_1 = _mm_add_sd(c46_1, _mm_mul_sd(a46_1, b46));
#endif
_mm_store_sd(&C[(i*84)+25], c46_1);
__m128d c46_2 = _mm_load_sd(&C[(i*84)+46]);
__m128d a46_2 = _mm_load_sd(&values[204]);
#if defined(__SSE3__) && defined(__AVX__)
c46_2 = _mm_add_sd(c46_2, _mm_mul_sd(a46_2, _mm256_castpd256_pd128(b46)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c46_2 = _mm_add_sd(c46_2, _mm_mul_sd(a46_2, b46));
#endif
_mm_store_sd(&C[(i*84)+46], c46_2);
__m128d c46_3 = _mm_load_sd(&C[(i*84)+74]);
__m128d a46_3 = _mm_load_sd(&values[205]);
#if defined(__SSE3__) && defined(__AVX__)
c46_3 = _mm_add_sd(c46_3, _mm_mul_sd(a46_3, _mm256_castpd256_pd128(b46)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c46_3 = _mm_add_sd(c46_3, _mm_mul_sd(a46_3, b46));
#endif
_mm_store_sd(&C[(i*84)+74], c46_3);
#else
C[(i*84)+10] += values[202] * B[(i*84)+46];
C[(i*84)+25] += values[203] * B[(i*84)+46];
C[(i*84)+46] += values[204] * B[(i*84)+46];
C[(i*84)+74] += values[205] * B[(i*84)+46];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b47 = _mm256_broadcast_sd(&B[(i*84)+47]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b47 = _mm_loaddup_pd(&B[(i*84)+47]);
#endif
__m128d c47_0 = _mm_load_sd(&C[(i*84)+11]);
__m128d a47_0 = _mm_load_sd(&values[206]);
#if defined(__SSE3__) && defined(__AVX__)
c47_0 = _mm_add_sd(c47_0, _mm_mul_sd(a47_0, _mm256_castpd256_pd128(b47)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c47_0 = _mm_add_sd(c47_0, _mm_mul_sd(a47_0, b47));
#endif
_mm_store_sd(&C[(i*84)+11], c47_0);
__m128d c47_1 = _mm_load_sd(&C[(i*84)+26]);
__m128d a47_1 = _mm_load_sd(&values[207]);
#if defined(__SSE3__) && defined(__AVX__)
c47_1 = _mm_add_sd(c47_1, _mm_mul_sd(a47_1, _mm256_castpd256_pd128(b47)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c47_1 = _mm_add_sd(c47_1, _mm_mul_sd(a47_1, b47));
#endif
_mm_store_sd(&C[(i*84)+26], c47_1);
__m128d c47_2 = _mm_load_sd(&C[(i*84)+47]);
__m128d a47_2 = _mm_load_sd(&values[208]);
#if defined(__SSE3__) && defined(__AVX__)
c47_2 = _mm_add_sd(c47_2, _mm_mul_sd(a47_2, _mm256_castpd256_pd128(b47)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c47_2 = _mm_add_sd(c47_2, _mm_mul_sd(a47_2, b47));
#endif
_mm_store_sd(&C[(i*84)+47], c47_2);
__m128d c47_3 = _mm_load_sd(&C[(i*84)+75]);
__m128d a47_3 = _mm_load_sd(&values[209]);
#if defined(__SSE3__) && defined(__AVX__)
c47_3 = _mm_add_sd(c47_3, _mm_mul_sd(a47_3, _mm256_castpd256_pd128(b47)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c47_3 = _mm_add_sd(c47_3, _mm_mul_sd(a47_3, b47));
#endif
_mm_store_sd(&C[(i*84)+75], c47_3);
#else
C[(i*84)+11] += values[206] * B[(i*84)+47];
C[(i*84)+26] += values[207] * B[(i*84)+47];
C[(i*84)+47] += values[208] * B[(i*84)+47];
C[(i*84)+75] += values[209] * B[(i*84)+47];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b48 = _mm256_broadcast_sd(&B[(i*84)+48]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b48 = _mm_loaddup_pd(&B[(i*84)+48]);
#endif
__m128d c48_0 = _mm_load_sd(&C[(i*84)+12]);
__m128d a48_0 = _mm_load_sd(&values[210]);
#if defined(__SSE3__) && defined(__AVX__)
c48_0 = _mm_add_sd(c48_0, _mm_mul_sd(a48_0, _mm256_castpd256_pd128(b48)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c48_0 = _mm_add_sd(c48_0, _mm_mul_sd(a48_0, b48));
#endif
_mm_store_sd(&C[(i*84)+12], c48_0);
__m128d c48_1 = _mm_load_sd(&C[(i*84)+27]);
__m128d a48_1 = _mm_load_sd(&values[211]);
#if defined(__SSE3__) && defined(__AVX__)
c48_1 = _mm_add_sd(c48_1, _mm_mul_sd(a48_1, _mm256_castpd256_pd128(b48)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c48_1 = _mm_add_sd(c48_1, _mm_mul_sd(a48_1, b48));
#endif
_mm_store_sd(&C[(i*84)+27], c48_1);
__m128d c48_2 = _mm_load_sd(&C[(i*84)+48]);
__m128d a48_2 = _mm_load_sd(&values[212]);
#if defined(__SSE3__) && defined(__AVX__)
c48_2 = _mm_add_sd(c48_2, _mm_mul_sd(a48_2, _mm256_castpd256_pd128(b48)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c48_2 = _mm_add_sd(c48_2, _mm_mul_sd(a48_2, b48));
#endif
_mm_store_sd(&C[(i*84)+48], c48_2);
__m128d c48_3 = _mm_load_sd(&C[(i*84)+76]);
__m128d a48_3 = _mm_load_sd(&values[213]);
#if defined(__SSE3__) && defined(__AVX__)
c48_3 = _mm_add_sd(c48_3, _mm_mul_sd(a48_3, _mm256_castpd256_pd128(b48)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c48_3 = _mm_add_sd(c48_3, _mm_mul_sd(a48_3, b48));
#endif
_mm_store_sd(&C[(i*84)+76], c48_3);
#else
C[(i*84)+12] += values[210] * B[(i*84)+48];
C[(i*84)+27] += values[211] * B[(i*84)+48];
C[(i*84)+48] += values[212] * B[(i*84)+48];
C[(i*84)+76] += values[213] * B[(i*84)+48];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b49 = _mm256_broadcast_sd(&B[(i*84)+49]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b49 = _mm_loaddup_pd(&B[(i*84)+49]);
#endif
__m128d c49_0 = _mm_load_sd(&C[(i*84)+13]);
__m128d a49_0 = _mm_load_sd(&values[214]);
#if defined(__SSE3__) && defined(__AVX__)
c49_0 = _mm_add_sd(c49_0, _mm_mul_sd(a49_0, _mm256_castpd256_pd128(b49)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c49_0 = _mm_add_sd(c49_0, _mm_mul_sd(a49_0, b49));
#endif
_mm_store_sd(&C[(i*84)+13], c49_0);
__m128d c49_1 = _mm_load_sd(&C[(i*84)+28]);
__m128d a49_1 = _mm_load_sd(&values[215]);
#if defined(__SSE3__) && defined(__AVX__)
c49_1 = _mm_add_sd(c49_1, _mm_mul_sd(a49_1, _mm256_castpd256_pd128(b49)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c49_1 = _mm_add_sd(c49_1, _mm_mul_sd(a49_1, b49));
#endif
_mm_store_sd(&C[(i*84)+28], c49_1);
__m128d c49_2 = _mm_load_sd(&C[(i*84)+49]);
__m128d a49_2 = _mm_load_sd(&values[216]);
#if defined(__SSE3__) && defined(__AVX__)
c49_2 = _mm_add_sd(c49_2, _mm_mul_sd(a49_2, _mm256_castpd256_pd128(b49)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c49_2 = _mm_add_sd(c49_2, _mm_mul_sd(a49_2, b49));
#endif
_mm_store_sd(&C[(i*84)+49], c49_2);
__m128d c49_3 = _mm_load_sd(&C[(i*84)+77]);
__m128d a49_3 = _mm_load_sd(&values[217]);
#if defined(__SSE3__) && defined(__AVX__)
c49_3 = _mm_add_sd(c49_3, _mm_mul_sd(a49_3, _mm256_castpd256_pd128(b49)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c49_3 = _mm_add_sd(c49_3, _mm_mul_sd(a49_3, b49));
#endif
_mm_store_sd(&C[(i*84)+77], c49_3);
#else
C[(i*84)+13] += values[214] * B[(i*84)+49];
C[(i*84)+28] += values[215] * B[(i*84)+49];
C[(i*84)+49] += values[216] * B[(i*84)+49];
C[(i*84)+77] += values[217] * B[(i*84)+49];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b50 = _mm256_broadcast_sd(&B[(i*84)+50]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b50 = _mm_loaddup_pd(&B[(i*84)+50]);
#endif
__m128d c50_0 = _mm_load_sd(&C[(i*84)+4]);
__m128d a50_0 = _mm_load_sd(&values[218]);
#if defined(__SSE3__) && defined(__AVX__)
c50_0 = _mm_add_sd(c50_0, _mm_mul_sd(a50_0, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c50_0 = _mm_add_sd(c50_0, _mm_mul_sd(a50_0, b50));
#endif
_mm_store_sd(&C[(i*84)+4], c50_0);
__m128d c50_1 = _mm_load_sd(&C[(i*84)+14]);
__m128d a50_1 = _mm_load_sd(&values[219]);
#if defined(__SSE3__) && defined(__AVX__)
c50_1 = _mm_add_sd(c50_1, _mm_mul_sd(a50_1, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c50_1 = _mm_add_sd(c50_1, _mm_mul_sd(a50_1, b50));
#endif
_mm_store_sd(&C[(i*84)+14], c50_1);
__m128d c50_2 = _mm_load_sd(&C[(i*84)+29]);
__m128d a50_2 = _mm_load_sd(&values[220]);
#if defined(__SSE3__) && defined(__AVX__)
c50_2 = _mm_add_sd(c50_2, _mm_mul_sd(a50_2, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c50_2 = _mm_add_sd(c50_2, _mm_mul_sd(a50_2, b50));
#endif
_mm_store_sd(&C[(i*84)+29], c50_2);
__m128d c50_3 = _mm_load_sd(&C[(i*84)+50]);
__m128d a50_3 = _mm_load_sd(&values[221]);
#if defined(__SSE3__) && defined(__AVX__)
c50_3 = _mm_add_sd(c50_3, _mm_mul_sd(a50_3, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c50_3 = _mm_add_sd(c50_3, _mm_mul_sd(a50_3, b50));
#endif
_mm_store_sd(&C[(i*84)+50], c50_3);
__m128d c50_4 = _mm_load_sd(&C[(i*84)+78]);
__m128d a50_4 = _mm_load_sd(&values[222]);
#if defined(__SSE3__) && defined(__AVX__)
c50_4 = _mm_add_sd(c50_4, _mm_mul_sd(a50_4, _mm256_castpd256_pd128(b50)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c50_4 = _mm_add_sd(c50_4, _mm_mul_sd(a50_4, b50));
#endif
_mm_store_sd(&C[(i*84)+78], c50_4);
#else
C[(i*84)+4] += values[218] * B[(i*84)+50];
C[(i*84)+14] += values[219] * B[(i*84)+50];
C[(i*84)+29] += values[220] * B[(i*84)+50];
C[(i*84)+50] += values[221] * B[(i*84)+50];
C[(i*84)+78] += values[222] * B[(i*84)+50];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b51 = _mm256_broadcast_sd(&B[(i*84)+51]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b51 = _mm_loaddup_pd(&B[(i*84)+51]);
#endif
__m128d c51_0 = _mm_load_sd(&C[(i*84)+5]);
__m128d a51_0 = _mm_load_sd(&values[223]);
#if defined(__SSE3__) && defined(__AVX__)
c51_0 = _mm_add_sd(c51_0, _mm_mul_sd(a51_0, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c51_0 = _mm_add_sd(c51_0, _mm_mul_sd(a51_0, b51));
#endif
_mm_store_sd(&C[(i*84)+5], c51_0);
__m128d c51_1 = _mm_load_sd(&C[(i*84)+15]);
__m128d a51_1 = _mm_load_sd(&values[224]);
#if defined(__SSE3__) && defined(__AVX__)
c51_1 = _mm_add_sd(c51_1, _mm_mul_sd(a51_1, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c51_1 = _mm_add_sd(c51_1, _mm_mul_sd(a51_1, b51));
#endif
_mm_store_sd(&C[(i*84)+15], c51_1);
__m128d c51_2 = _mm_load_sd(&C[(i*84)+30]);
__m128d a51_2 = _mm_load_sd(&values[225]);
#if defined(__SSE3__) && defined(__AVX__)
c51_2 = _mm_add_sd(c51_2, _mm_mul_sd(a51_2, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c51_2 = _mm_add_sd(c51_2, _mm_mul_sd(a51_2, b51));
#endif
_mm_store_sd(&C[(i*84)+30], c51_2);
__m128d c51_3 = _mm_load_sd(&C[(i*84)+51]);
__m128d a51_3 = _mm_load_sd(&values[226]);
#if defined(__SSE3__) && defined(__AVX__)
c51_3 = _mm_add_sd(c51_3, _mm_mul_sd(a51_3, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c51_3 = _mm_add_sd(c51_3, _mm_mul_sd(a51_3, b51));
#endif
_mm_store_sd(&C[(i*84)+51], c51_3);
__m128d c51_4 = _mm_load_sd(&C[(i*84)+79]);
__m128d a51_4 = _mm_load_sd(&values[227]);
#if defined(__SSE3__) && defined(__AVX__)
c51_4 = _mm_add_sd(c51_4, _mm_mul_sd(a51_4, _mm256_castpd256_pd128(b51)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c51_4 = _mm_add_sd(c51_4, _mm_mul_sd(a51_4, b51));
#endif
_mm_store_sd(&C[(i*84)+79], c51_4);
#else
C[(i*84)+5] += values[223] * B[(i*84)+51];
C[(i*84)+15] += values[224] * B[(i*84)+51];
C[(i*84)+30] += values[225] * B[(i*84)+51];
C[(i*84)+51] += values[226] * B[(i*84)+51];
C[(i*84)+79] += values[227] * B[(i*84)+51];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b52 = _mm256_broadcast_sd(&B[(i*84)+52]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b52 = _mm_loaddup_pd(&B[(i*84)+52]);
#endif
__m128d c52_0 = _mm_load_sd(&C[(i*84)+6]);
__m128d a52_0 = _mm_load_sd(&values[228]);
#if defined(__SSE3__) && defined(__AVX__)
c52_0 = _mm_add_sd(c52_0, _mm_mul_sd(a52_0, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c52_0 = _mm_add_sd(c52_0, _mm_mul_sd(a52_0, b52));
#endif
_mm_store_sd(&C[(i*84)+6], c52_0);
__m128d c52_1 = _mm_load_sd(&C[(i*84)+16]);
__m128d a52_1 = _mm_load_sd(&values[229]);
#if defined(__SSE3__) && defined(__AVX__)
c52_1 = _mm_add_sd(c52_1, _mm_mul_sd(a52_1, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c52_1 = _mm_add_sd(c52_1, _mm_mul_sd(a52_1, b52));
#endif
_mm_store_sd(&C[(i*84)+16], c52_1);
__m128d c52_2 = _mm_load_sd(&C[(i*84)+31]);
__m128d a52_2 = _mm_load_sd(&values[230]);
#if defined(__SSE3__) && defined(__AVX__)
c52_2 = _mm_add_sd(c52_2, _mm_mul_sd(a52_2, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c52_2 = _mm_add_sd(c52_2, _mm_mul_sd(a52_2, b52));
#endif
_mm_store_sd(&C[(i*84)+31], c52_2);
__m128d c52_3 = _mm_load_sd(&C[(i*84)+52]);
__m128d a52_3 = _mm_load_sd(&values[231]);
#if defined(__SSE3__) && defined(__AVX__)
c52_3 = _mm_add_sd(c52_3, _mm_mul_sd(a52_3, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c52_3 = _mm_add_sd(c52_3, _mm_mul_sd(a52_3, b52));
#endif
_mm_store_sd(&C[(i*84)+52], c52_3);
__m128d c52_4 = _mm_load_sd(&C[(i*84)+80]);
__m128d a52_4 = _mm_load_sd(&values[232]);
#if defined(__SSE3__) && defined(__AVX__)
c52_4 = _mm_add_sd(c52_4, _mm_mul_sd(a52_4, _mm256_castpd256_pd128(b52)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c52_4 = _mm_add_sd(c52_4, _mm_mul_sd(a52_4, b52));
#endif
_mm_store_sd(&C[(i*84)+80], c52_4);
#else
C[(i*84)+6] += values[228] * B[(i*84)+52];
C[(i*84)+16] += values[229] * B[(i*84)+52];
C[(i*84)+31] += values[230] * B[(i*84)+52];
C[(i*84)+52] += values[231] * B[(i*84)+52];
C[(i*84)+80] += values[232] * B[(i*84)+52];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b53 = _mm256_broadcast_sd(&B[(i*84)+53]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b53 = _mm_loaddup_pd(&B[(i*84)+53]);
#endif
__m128d c53_0 = _mm_load_sd(&C[(i*84)+1]);
__m128d a53_0 = _mm_load_sd(&values[233]);
#if defined(__SSE3__) && defined(__AVX__)
c53_0 = _mm_add_sd(c53_0, _mm_mul_sd(a53_0, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c53_0 = _mm_add_sd(c53_0, _mm_mul_sd(a53_0, b53));
#endif
_mm_store_sd(&C[(i*84)+1], c53_0);
__m128d c53_1 = _mm_load_sd(&C[(i*84)+7]);
__m128d a53_1 = _mm_load_sd(&values[234]);
#if defined(__SSE3__) && defined(__AVX__)
c53_1 = _mm_add_sd(c53_1, _mm_mul_sd(a53_1, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c53_1 = _mm_add_sd(c53_1, _mm_mul_sd(a53_1, b53));
#endif
_mm_store_sd(&C[(i*84)+7], c53_1);
__m128d c53_2 = _mm_load_sd(&C[(i*84)+17]);
__m128d a53_2 = _mm_load_sd(&values[235]);
#if defined(__SSE3__) && defined(__AVX__)
c53_2 = _mm_add_sd(c53_2, _mm_mul_sd(a53_2, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c53_2 = _mm_add_sd(c53_2, _mm_mul_sd(a53_2, b53));
#endif
_mm_store_sd(&C[(i*84)+17], c53_2);
__m128d c53_3 = _mm_load_sd(&C[(i*84)+32]);
__m128d a53_3 = _mm_load_sd(&values[236]);
#if defined(__SSE3__) && defined(__AVX__)
c53_3 = _mm_add_sd(c53_3, _mm_mul_sd(a53_3, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c53_3 = _mm_add_sd(c53_3, _mm_mul_sd(a53_3, b53));
#endif
_mm_store_sd(&C[(i*84)+32], c53_3);
__m128d c53_4 = _mm_load_sd(&C[(i*84)+53]);
__m128d a53_4 = _mm_load_sd(&values[237]);
#if defined(__SSE3__) && defined(__AVX__)
c53_4 = _mm_add_sd(c53_4, _mm_mul_sd(a53_4, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c53_4 = _mm_add_sd(c53_4, _mm_mul_sd(a53_4, b53));
#endif
_mm_store_sd(&C[(i*84)+53], c53_4);
__m128d c53_5 = _mm_load_sd(&C[(i*84)+81]);
__m128d a53_5 = _mm_load_sd(&values[238]);
#if defined(__SSE3__) && defined(__AVX__)
c53_5 = _mm_add_sd(c53_5, _mm_mul_sd(a53_5, _mm256_castpd256_pd128(b53)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c53_5 = _mm_add_sd(c53_5, _mm_mul_sd(a53_5, b53));
#endif
_mm_store_sd(&C[(i*84)+81], c53_5);
#else
C[(i*84)+1] += values[233] * B[(i*84)+53];
C[(i*84)+7] += values[234] * B[(i*84)+53];
C[(i*84)+17] += values[235] * B[(i*84)+53];
C[(i*84)+32] += values[236] * B[(i*84)+53];
C[(i*84)+53] += values[237] * B[(i*84)+53];
C[(i*84)+81] += values[238] * B[(i*84)+53];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b54 = _mm256_broadcast_sd(&B[(i*84)+54]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b54 = _mm_loaddup_pd(&B[(i*84)+54]);
#endif
__m128d c54_0 = _mm_load_sd(&C[(i*84)+2]);
__m128d a54_0 = _mm_load_sd(&values[239]);
#if defined(__SSE3__) && defined(__AVX__)
c54_0 = _mm_add_sd(c54_0, _mm_mul_sd(a54_0, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c54_0 = _mm_add_sd(c54_0, _mm_mul_sd(a54_0, b54));
#endif
_mm_store_sd(&C[(i*84)+2], c54_0);
__m128d c54_1 = _mm_load_sd(&C[(i*84)+8]);
__m128d a54_1 = _mm_load_sd(&values[240]);
#if defined(__SSE3__) && defined(__AVX__)
c54_1 = _mm_add_sd(c54_1, _mm_mul_sd(a54_1, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c54_1 = _mm_add_sd(c54_1, _mm_mul_sd(a54_1, b54));
#endif
_mm_store_sd(&C[(i*84)+8], c54_1);
__m128d c54_2 = _mm_load_sd(&C[(i*84)+18]);
__m128d a54_2 = _mm_load_sd(&values[241]);
#if defined(__SSE3__) && defined(__AVX__)
c54_2 = _mm_add_sd(c54_2, _mm_mul_sd(a54_2, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c54_2 = _mm_add_sd(c54_2, _mm_mul_sd(a54_2, b54));
#endif
_mm_store_sd(&C[(i*84)+18], c54_2);
__m128d c54_3 = _mm_load_sd(&C[(i*84)+33]);
__m128d a54_3 = _mm_load_sd(&values[242]);
#if defined(__SSE3__) && defined(__AVX__)
c54_3 = _mm_add_sd(c54_3, _mm_mul_sd(a54_3, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c54_3 = _mm_add_sd(c54_3, _mm_mul_sd(a54_3, b54));
#endif
_mm_store_sd(&C[(i*84)+33], c54_3);
__m128d c54_4 = _mm_load_sd(&C[(i*84)+54]);
__m128d a54_4 = _mm_load_sd(&values[243]);
#if defined(__SSE3__) && defined(__AVX__)
c54_4 = _mm_add_sd(c54_4, _mm_mul_sd(a54_4, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c54_4 = _mm_add_sd(c54_4, _mm_mul_sd(a54_4, b54));
#endif
_mm_store_sd(&C[(i*84)+54], c54_4);
__m128d c54_5 = _mm_load_sd(&C[(i*84)+82]);
__m128d a54_5 = _mm_load_sd(&values[244]);
#if defined(__SSE3__) && defined(__AVX__)
c54_5 = _mm_add_sd(c54_5, _mm_mul_sd(a54_5, _mm256_castpd256_pd128(b54)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c54_5 = _mm_add_sd(c54_5, _mm_mul_sd(a54_5, b54));
#endif
_mm_store_sd(&C[(i*84)+82], c54_5);
#else
C[(i*84)+2] += values[239] * B[(i*84)+54];
C[(i*84)+8] += values[240] * B[(i*84)+54];
C[(i*84)+18] += values[241] * B[(i*84)+54];
C[(i*84)+33] += values[242] * B[(i*84)+54];
C[(i*84)+54] += values[243] * B[(i*84)+54];
C[(i*84)+82] += values[244] * B[(i*84)+54];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b55 = _mm256_broadcast_sd(&B[(i*84)+55]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b55 = _mm_loaddup_pd(&B[(i*84)+55]);
#endif
__m128d c55_0 = _mm_load_sd(&C[(i*84)+0]);
__m128d a55_0 = _mm_load_sd(&values[245]);
#if defined(__SSE3__) && defined(__AVX__)
c55_0 = _mm_add_sd(c55_0, _mm_mul_sd(a55_0, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_0 = _mm_add_sd(c55_0, _mm_mul_sd(a55_0, b55));
#endif
_mm_store_sd(&C[(i*84)+0], c55_0);
__m128d c55_1 = _mm_load_sd(&C[(i*84)+3]);
__m128d a55_1 = _mm_load_sd(&values[246]);
#if defined(__SSE3__) && defined(__AVX__)
c55_1 = _mm_add_sd(c55_1, _mm_mul_sd(a55_1, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_1 = _mm_add_sd(c55_1, _mm_mul_sd(a55_1, b55));
#endif
_mm_store_sd(&C[(i*84)+3], c55_1);
__m128d c55_2 = _mm_load_sd(&C[(i*84)+9]);
__m128d a55_2 = _mm_load_sd(&values[247]);
#if defined(__SSE3__) && defined(__AVX__)
c55_2 = _mm_add_sd(c55_2, _mm_mul_sd(a55_2, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_2 = _mm_add_sd(c55_2, _mm_mul_sd(a55_2, b55));
#endif
_mm_store_sd(&C[(i*84)+9], c55_2);
__m128d c55_3 = _mm_load_sd(&C[(i*84)+19]);
__m128d a55_3 = _mm_load_sd(&values[248]);
#if defined(__SSE3__) && defined(__AVX__)
c55_3 = _mm_add_sd(c55_3, _mm_mul_sd(a55_3, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_3 = _mm_add_sd(c55_3, _mm_mul_sd(a55_3, b55));
#endif
_mm_store_sd(&C[(i*84)+19], c55_3);
__m128d c55_4 = _mm_load_sd(&C[(i*84)+34]);
__m128d a55_4 = _mm_load_sd(&values[249]);
#if defined(__SSE3__) && defined(__AVX__)
c55_4 = _mm_add_sd(c55_4, _mm_mul_sd(a55_4, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_4 = _mm_add_sd(c55_4, _mm_mul_sd(a55_4, b55));
#endif
_mm_store_sd(&C[(i*84)+34], c55_4);
__m128d c55_5 = _mm_load_sd(&C[(i*84)+55]);
__m128d a55_5 = _mm_load_sd(&values[250]);
#if defined(__SSE3__) && defined(__AVX__)
c55_5 = _mm_add_sd(c55_5, _mm_mul_sd(a55_5, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_5 = _mm_add_sd(c55_5, _mm_mul_sd(a55_5, b55));
#endif
_mm_store_sd(&C[(i*84)+55], c55_5);
__m128d c55_6 = _mm_load_sd(&C[(i*84)+83]);
__m128d a55_6 = _mm_load_sd(&values[251]);
#if defined(__SSE3__) && defined(__AVX__)
c55_6 = _mm_add_sd(c55_6, _mm_mul_sd(a55_6, _mm256_castpd256_pd128(b55)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c55_6 = _mm_add_sd(c55_6, _mm_mul_sd(a55_6, b55));
#endif
_mm_store_sd(&C[(i*84)+83], c55_6);
#else
C[(i*84)+0] += values[245] * B[(i*84)+55];
C[(i*84)+3] += values[246] * B[(i*84)+55];
C[(i*84)+9] += values[247] * B[(i*84)+55];
C[(i*84)+19] += values[248] * B[(i*84)+55];
C[(i*84)+34] += values[249] * B[(i*84)+55];
C[(i*84)+55] += values[250] * B[(i*84)+55];
C[(i*84)+83] += values[251] * B[(i*84)+55];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b56 = _mm256_broadcast_sd(&B[(i*84)+56]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b56 = _mm_loaddup_pd(&B[(i*84)+56]);
#endif
__m128d c56_0 = _mm_load_sd(&C[(i*84)+56]);
__m128d a56_0 = _mm_load_sd(&values[252]);
#if defined(__SSE3__) && defined(__AVX__)
c56_0 = _mm_add_sd(c56_0, _mm_mul_sd(a56_0, _mm256_castpd256_pd128(b56)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c56_0 = _mm_add_sd(c56_0, _mm_mul_sd(a56_0, b56));
#endif
_mm_store_sd(&C[(i*84)+56], c56_0);
#else
C[(i*84)+56] += values[252] * B[(i*84)+56];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b57 = _mm256_broadcast_sd(&B[(i*84)+57]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b57 = _mm_loaddup_pd(&B[(i*84)+57]);
#endif
__m128d c57_0 = _mm_load_sd(&C[(i*84)+57]);
__m128d a57_0 = _mm_load_sd(&values[253]);
#if defined(__SSE3__) && defined(__AVX__)
c57_0 = _mm_add_sd(c57_0, _mm_mul_sd(a57_0, _mm256_castpd256_pd128(b57)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c57_0 = _mm_add_sd(c57_0, _mm_mul_sd(a57_0, b57));
#endif
_mm_store_sd(&C[(i*84)+57], c57_0);
#else
C[(i*84)+57] += values[253] * B[(i*84)+57];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b58 = _mm256_broadcast_sd(&B[(i*84)+58]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b58 = _mm_loaddup_pd(&B[(i*84)+58]);
#endif
__m128d c58_0 = _mm_load_sd(&C[(i*84)+58]);
__m128d a58_0 = _mm_load_sd(&values[254]);
#if defined(__SSE3__) && defined(__AVX__)
c58_0 = _mm_add_sd(c58_0, _mm_mul_sd(a58_0, _mm256_castpd256_pd128(b58)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c58_0 = _mm_add_sd(c58_0, _mm_mul_sd(a58_0, b58));
#endif
_mm_store_sd(&C[(i*84)+58], c58_0);
#else
C[(i*84)+58] += values[254] * B[(i*84)+58];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b59 = _mm256_broadcast_sd(&B[(i*84)+59]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b59 = _mm_loaddup_pd(&B[(i*84)+59]);
#endif
__m128d c59_0 = _mm_load_sd(&C[(i*84)+59]);
__m128d a59_0 = _mm_load_sd(&values[255]);
#if defined(__SSE3__) && defined(__AVX__)
c59_0 = _mm_add_sd(c59_0, _mm_mul_sd(a59_0, _mm256_castpd256_pd128(b59)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c59_0 = _mm_add_sd(c59_0, _mm_mul_sd(a59_0, b59));
#endif
_mm_store_sd(&C[(i*84)+59], c59_0);
#else
C[(i*84)+59] += values[255] * B[(i*84)+59];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b60 = _mm256_broadcast_sd(&B[(i*84)+60]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b60 = _mm_loaddup_pd(&B[(i*84)+60]);
#endif
__m128d c60_0 = _mm_load_sd(&C[(i*84)+60]);
__m128d a60_0 = _mm_load_sd(&values[256]);
#if defined(__SSE3__) && defined(__AVX__)
c60_0 = _mm_add_sd(c60_0, _mm_mul_sd(a60_0, _mm256_castpd256_pd128(b60)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c60_0 = _mm_add_sd(c60_0, _mm_mul_sd(a60_0, b60));
#endif
_mm_store_sd(&C[(i*84)+60], c60_0);
#else
C[(i*84)+60] += values[256] * B[(i*84)+60];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b61 = _mm256_broadcast_sd(&B[(i*84)+61]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b61 = _mm_loaddup_pd(&B[(i*84)+61]);
#endif
__m128d c61_0 = _mm_load_sd(&C[(i*84)+61]);
__m128d a61_0 = _mm_load_sd(&values[257]);
#if defined(__SSE3__) && defined(__AVX__)
c61_0 = _mm_add_sd(c61_0, _mm_mul_sd(a61_0, _mm256_castpd256_pd128(b61)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c61_0 = _mm_add_sd(c61_0, _mm_mul_sd(a61_0, b61));
#endif
_mm_store_sd(&C[(i*84)+61], c61_0);
#else
C[(i*84)+61] += values[257] * B[(i*84)+61];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b62 = _mm256_broadcast_sd(&B[(i*84)+62]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b62 = _mm_loaddup_pd(&B[(i*84)+62]);
#endif
__m128d c62_0 = _mm_load_sd(&C[(i*84)+62]);
__m128d a62_0 = _mm_load_sd(&values[258]);
#if defined(__SSE3__) && defined(__AVX__)
c62_0 = _mm_add_sd(c62_0, _mm_mul_sd(a62_0, _mm256_castpd256_pd128(b62)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c62_0 = _mm_add_sd(c62_0, _mm_mul_sd(a62_0, b62));
#endif
_mm_store_sd(&C[(i*84)+62], c62_0);
#else
C[(i*84)+62] += values[258] * B[(i*84)+62];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b63 = _mm256_broadcast_sd(&B[(i*84)+63]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b63 = _mm_loaddup_pd(&B[(i*84)+63]);
#endif
__m128d c63_0 = _mm_load_sd(&C[(i*84)+35]);
__m128d a63_0 = _mm_load_sd(&values[259]);
#if defined(__SSE3__) && defined(__AVX__)
c63_0 = _mm_add_sd(c63_0, _mm_mul_sd(a63_0, _mm256_castpd256_pd128(b63)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c63_0 = _mm_add_sd(c63_0, _mm_mul_sd(a63_0, b63));
#endif
_mm_store_sd(&C[(i*84)+35], c63_0);
__m128d c63_1 = _mm_load_sd(&C[(i*84)+63]);
__m128d a63_1 = _mm_load_sd(&values[260]);
#if defined(__SSE3__) && defined(__AVX__)
c63_1 = _mm_add_sd(c63_1, _mm_mul_sd(a63_1, _mm256_castpd256_pd128(b63)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c63_1 = _mm_add_sd(c63_1, _mm_mul_sd(a63_1, b63));
#endif
_mm_store_sd(&C[(i*84)+63], c63_1);
#else
C[(i*84)+35] += values[259] * B[(i*84)+63];
C[(i*84)+63] += values[260] * B[(i*84)+63];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b64 = _mm256_broadcast_sd(&B[(i*84)+64]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b64 = _mm_loaddup_pd(&B[(i*84)+64]);
#endif
__m128d c64_0 = _mm_load_sd(&C[(i*84)+36]);
__m128d a64_0 = _mm_load_sd(&values[261]);
#if defined(__SSE3__) && defined(__AVX__)
c64_0 = _mm_add_sd(c64_0, _mm_mul_sd(a64_0, _mm256_castpd256_pd128(b64)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c64_0 = _mm_add_sd(c64_0, _mm_mul_sd(a64_0, b64));
#endif
_mm_store_sd(&C[(i*84)+36], c64_0);
__m128d c64_1 = _mm_load_sd(&C[(i*84)+64]);
__m128d a64_1 = _mm_load_sd(&values[262]);
#if defined(__SSE3__) && defined(__AVX__)
c64_1 = _mm_add_sd(c64_1, _mm_mul_sd(a64_1, _mm256_castpd256_pd128(b64)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c64_1 = _mm_add_sd(c64_1, _mm_mul_sd(a64_1, b64));
#endif
_mm_store_sd(&C[(i*84)+64], c64_1);
#else
C[(i*84)+36] += values[261] * B[(i*84)+64];
C[(i*84)+64] += values[262] * B[(i*84)+64];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b65 = _mm256_broadcast_sd(&B[(i*84)+65]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b65 = _mm_loaddup_pd(&B[(i*84)+65]);
#endif
__m128d c65_0 = _mm_load_sd(&C[(i*84)+37]);
__m128d a65_0 = _mm_load_sd(&values[263]);
#if defined(__SSE3__) && defined(__AVX__)
c65_0 = _mm_add_sd(c65_0, _mm_mul_sd(a65_0, _mm256_castpd256_pd128(b65)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c65_0 = _mm_add_sd(c65_0, _mm_mul_sd(a65_0, b65));
#endif
_mm_store_sd(&C[(i*84)+37], c65_0);
__m128d c65_1 = _mm_load_sd(&C[(i*84)+65]);
__m128d a65_1 = _mm_load_sd(&values[264]);
#if defined(__SSE3__) && defined(__AVX__)
c65_1 = _mm_add_sd(c65_1, _mm_mul_sd(a65_1, _mm256_castpd256_pd128(b65)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c65_1 = _mm_add_sd(c65_1, _mm_mul_sd(a65_1, b65));
#endif
_mm_store_sd(&C[(i*84)+65], c65_1);
#else
C[(i*84)+37] += values[263] * B[(i*84)+65];
C[(i*84)+65] += values[264] * B[(i*84)+65];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b66 = _mm256_broadcast_sd(&B[(i*84)+66]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b66 = _mm_loaddup_pd(&B[(i*84)+66]);
#endif
__m128d c66_0 = _mm_load_sd(&C[(i*84)+38]);
__m128d a66_0 = _mm_load_sd(&values[265]);
#if defined(__SSE3__) && defined(__AVX__)
c66_0 = _mm_add_sd(c66_0, _mm_mul_sd(a66_0, _mm256_castpd256_pd128(b66)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c66_0 = _mm_add_sd(c66_0, _mm_mul_sd(a66_0, b66));
#endif
_mm_store_sd(&C[(i*84)+38], c66_0);
__m128d c66_1 = _mm_load_sd(&C[(i*84)+66]);
__m128d a66_1 = _mm_load_sd(&values[266]);
#if defined(__SSE3__) && defined(__AVX__)
c66_1 = _mm_add_sd(c66_1, _mm_mul_sd(a66_1, _mm256_castpd256_pd128(b66)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c66_1 = _mm_add_sd(c66_1, _mm_mul_sd(a66_1, b66));
#endif
_mm_store_sd(&C[(i*84)+66], c66_1);
#else
C[(i*84)+38] += values[265] * B[(i*84)+66];
C[(i*84)+66] += values[266] * B[(i*84)+66];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b67 = _mm256_broadcast_sd(&B[(i*84)+67]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b67 = _mm_loaddup_pd(&B[(i*84)+67]);
#endif
__m128d c67_0 = _mm_load_sd(&C[(i*84)+39]);
__m128d a67_0 = _mm_load_sd(&values[267]);
#if defined(__SSE3__) && defined(__AVX__)
c67_0 = _mm_add_sd(c67_0, _mm_mul_sd(a67_0, _mm256_castpd256_pd128(b67)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c67_0 = _mm_add_sd(c67_0, _mm_mul_sd(a67_0, b67));
#endif
_mm_store_sd(&C[(i*84)+39], c67_0);
__m128d c67_1 = _mm_load_sd(&C[(i*84)+67]);
__m128d a67_1 = _mm_load_sd(&values[268]);
#if defined(__SSE3__) && defined(__AVX__)
c67_1 = _mm_add_sd(c67_1, _mm_mul_sd(a67_1, _mm256_castpd256_pd128(b67)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c67_1 = _mm_add_sd(c67_1, _mm_mul_sd(a67_1, b67));
#endif
_mm_store_sd(&C[(i*84)+67], c67_1);
#else
C[(i*84)+39] += values[267] * B[(i*84)+67];
C[(i*84)+67] += values[268] * B[(i*84)+67];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b68 = _mm256_broadcast_sd(&B[(i*84)+68]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b68 = _mm_loaddup_pd(&B[(i*84)+68]);
#endif
__m128d c68_0 = _mm_load_sd(&C[(i*84)+40]);
__m128d a68_0 = _mm_load_sd(&values[269]);
#if defined(__SSE3__) && defined(__AVX__)
c68_0 = _mm_add_sd(c68_0, _mm_mul_sd(a68_0, _mm256_castpd256_pd128(b68)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c68_0 = _mm_add_sd(c68_0, _mm_mul_sd(a68_0, b68));
#endif
_mm_store_sd(&C[(i*84)+40], c68_0);
__m128d c68_1 = _mm_load_sd(&C[(i*84)+68]);
__m128d a68_1 = _mm_load_sd(&values[270]);
#if defined(__SSE3__) && defined(__AVX__)
c68_1 = _mm_add_sd(c68_1, _mm_mul_sd(a68_1, _mm256_castpd256_pd128(b68)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c68_1 = _mm_add_sd(c68_1, _mm_mul_sd(a68_1, b68));
#endif
_mm_store_sd(&C[(i*84)+68], c68_1);
#else
C[(i*84)+40] += values[269] * B[(i*84)+68];
C[(i*84)+68] += values[270] * B[(i*84)+68];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b69 = _mm256_broadcast_sd(&B[(i*84)+69]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b69 = _mm_loaddup_pd(&B[(i*84)+69]);
#endif
__m128d c69_0 = _mm_load_sd(&C[(i*84)+20]);
__m128d a69_0 = _mm_load_sd(&values[271]);
#if defined(__SSE3__) && defined(__AVX__)
c69_0 = _mm_add_sd(c69_0, _mm_mul_sd(a69_0, _mm256_castpd256_pd128(b69)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c69_0 = _mm_add_sd(c69_0, _mm_mul_sd(a69_0, b69));
#endif
_mm_store_sd(&C[(i*84)+20], c69_0);
__m128d c69_1 = _mm_load_sd(&C[(i*84)+41]);
__m128d a69_1 = _mm_load_sd(&values[272]);
#if defined(__SSE3__) && defined(__AVX__)
c69_1 = _mm_add_sd(c69_1, _mm_mul_sd(a69_1, _mm256_castpd256_pd128(b69)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c69_1 = _mm_add_sd(c69_1, _mm_mul_sd(a69_1, b69));
#endif
_mm_store_sd(&C[(i*84)+41], c69_1);
__m128d c69_2 = _mm_load_sd(&C[(i*84)+69]);
__m128d a69_2 = _mm_load_sd(&values[273]);
#if defined(__SSE3__) && defined(__AVX__)
c69_2 = _mm_add_sd(c69_2, _mm_mul_sd(a69_2, _mm256_castpd256_pd128(b69)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c69_2 = _mm_add_sd(c69_2, _mm_mul_sd(a69_2, b69));
#endif
_mm_store_sd(&C[(i*84)+69], c69_2);
#else
C[(i*84)+20] += values[271] * B[(i*84)+69];
C[(i*84)+41] += values[272] * B[(i*84)+69];
C[(i*84)+69] += values[273] * B[(i*84)+69];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b70 = _mm256_broadcast_sd(&B[(i*84)+70]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b70 = _mm_loaddup_pd(&B[(i*84)+70]);
#endif
__m128d c70_0 = _mm_load_sd(&C[(i*84)+21]);
__m128d a70_0 = _mm_load_sd(&values[274]);
#if defined(__SSE3__) && defined(__AVX__)
c70_0 = _mm_add_sd(c70_0, _mm_mul_sd(a70_0, _mm256_castpd256_pd128(b70)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c70_0 = _mm_add_sd(c70_0, _mm_mul_sd(a70_0, b70));
#endif
_mm_store_sd(&C[(i*84)+21], c70_0);
__m128d c70_1 = _mm_load_sd(&C[(i*84)+42]);
__m128d a70_1 = _mm_load_sd(&values[275]);
#if defined(__SSE3__) && defined(__AVX__)
c70_1 = _mm_add_sd(c70_1, _mm_mul_sd(a70_1, _mm256_castpd256_pd128(b70)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c70_1 = _mm_add_sd(c70_1, _mm_mul_sd(a70_1, b70));
#endif
_mm_store_sd(&C[(i*84)+42], c70_1);
__m128d c70_2 = _mm_load_sd(&C[(i*84)+70]);
__m128d a70_2 = _mm_load_sd(&values[276]);
#if defined(__SSE3__) && defined(__AVX__)
c70_2 = _mm_add_sd(c70_2, _mm_mul_sd(a70_2, _mm256_castpd256_pd128(b70)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c70_2 = _mm_add_sd(c70_2, _mm_mul_sd(a70_2, b70));
#endif
_mm_store_sd(&C[(i*84)+70], c70_2);
#else
C[(i*84)+21] += values[274] * B[(i*84)+70];
C[(i*84)+42] += values[275] * B[(i*84)+70];
C[(i*84)+70] += values[276] * B[(i*84)+70];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b71 = _mm256_broadcast_sd(&B[(i*84)+71]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b71 = _mm_loaddup_pd(&B[(i*84)+71]);
#endif
__m128d c71_0 = _mm_load_sd(&C[(i*84)+22]);
__m128d a71_0 = _mm_load_sd(&values[277]);
#if defined(__SSE3__) && defined(__AVX__)
c71_0 = _mm_add_sd(c71_0, _mm_mul_sd(a71_0, _mm256_castpd256_pd128(b71)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c71_0 = _mm_add_sd(c71_0, _mm_mul_sd(a71_0, b71));
#endif
_mm_store_sd(&C[(i*84)+22], c71_0);
__m128d c71_1 = _mm_load_sd(&C[(i*84)+43]);
__m128d a71_1 = _mm_load_sd(&values[278]);
#if defined(__SSE3__) && defined(__AVX__)
c71_1 = _mm_add_sd(c71_1, _mm_mul_sd(a71_1, _mm256_castpd256_pd128(b71)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c71_1 = _mm_add_sd(c71_1, _mm_mul_sd(a71_1, b71));
#endif
_mm_store_sd(&C[(i*84)+43], c71_1);
__m128d c71_2 = _mm_load_sd(&C[(i*84)+71]);
__m128d a71_2 = _mm_load_sd(&values[279]);
#if defined(__SSE3__) && defined(__AVX__)
c71_2 = _mm_add_sd(c71_2, _mm_mul_sd(a71_2, _mm256_castpd256_pd128(b71)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c71_2 = _mm_add_sd(c71_2, _mm_mul_sd(a71_2, b71));
#endif
_mm_store_sd(&C[(i*84)+71], c71_2);
#else
C[(i*84)+22] += values[277] * B[(i*84)+71];
C[(i*84)+43] += values[278] * B[(i*84)+71];
C[(i*84)+71] += values[279] * B[(i*84)+71];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b72 = _mm256_broadcast_sd(&B[(i*84)+72]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b72 = _mm_loaddup_pd(&B[(i*84)+72]);
#endif
__m128d c72_0 = _mm_load_sd(&C[(i*84)+23]);
__m128d a72_0 = _mm_load_sd(&values[280]);
#if defined(__SSE3__) && defined(__AVX__)
c72_0 = _mm_add_sd(c72_0, _mm_mul_sd(a72_0, _mm256_castpd256_pd128(b72)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c72_0 = _mm_add_sd(c72_0, _mm_mul_sd(a72_0, b72));
#endif
_mm_store_sd(&C[(i*84)+23], c72_0);
__m128d c72_1 = _mm_load_sd(&C[(i*84)+44]);
__m128d a72_1 = _mm_load_sd(&values[281]);
#if defined(__SSE3__) && defined(__AVX__)
c72_1 = _mm_add_sd(c72_1, _mm_mul_sd(a72_1, _mm256_castpd256_pd128(b72)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c72_1 = _mm_add_sd(c72_1, _mm_mul_sd(a72_1, b72));
#endif
_mm_store_sd(&C[(i*84)+44], c72_1);
__m128d c72_2 = _mm_load_sd(&C[(i*84)+72]);
__m128d a72_2 = _mm_load_sd(&values[282]);
#if defined(__SSE3__) && defined(__AVX__)
c72_2 = _mm_add_sd(c72_2, _mm_mul_sd(a72_2, _mm256_castpd256_pd128(b72)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c72_2 = _mm_add_sd(c72_2, _mm_mul_sd(a72_2, b72));
#endif
_mm_store_sd(&C[(i*84)+72], c72_2);
#else
C[(i*84)+23] += values[280] * B[(i*84)+72];
C[(i*84)+44] += values[281] * B[(i*84)+72];
C[(i*84)+72] += values[282] * B[(i*84)+72];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b73 = _mm256_broadcast_sd(&B[(i*84)+73]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b73 = _mm_loaddup_pd(&B[(i*84)+73]);
#endif
__m128d c73_0 = _mm_load_sd(&C[(i*84)+24]);
__m128d a73_0 = _mm_load_sd(&values[283]);
#if defined(__SSE3__) && defined(__AVX__)
c73_0 = _mm_add_sd(c73_0, _mm_mul_sd(a73_0, _mm256_castpd256_pd128(b73)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c73_0 = _mm_add_sd(c73_0, _mm_mul_sd(a73_0, b73));
#endif
_mm_store_sd(&C[(i*84)+24], c73_0);
__m128d c73_1 = _mm_load_sd(&C[(i*84)+45]);
__m128d a73_1 = _mm_load_sd(&values[284]);
#if defined(__SSE3__) && defined(__AVX__)
c73_1 = _mm_add_sd(c73_1, _mm_mul_sd(a73_1, _mm256_castpd256_pd128(b73)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c73_1 = _mm_add_sd(c73_1, _mm_mul_sd(a73_1, b73));
#endif
_mm_store_sd(&C[(i*84)+45], c73_1);
__m128d c73_2 = _mm_load_sd(&C[(i*84)+73]);
__m128d a73_2 = _mm_load_sd(&values[285]);
#if defined(__SSE3__) && defined(__AVX__)
c73_2 = _mm_add_sd(c73_2, _mm_mul_sd(a73_2, _mm256_castpd256_pd128(b73)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c73_2 = _mm_add_sd(c73_2, _mm_mul_sd(a73_2, b73));
#endif
_mm_store_sd(&C[(i*84)+73], c73_2);
#else
C[(i*84)+24] += values[283] * B[(i*84)+73];
C[(i*84)+45] += values[284] * B[(i*84)+73];
C[(i*84)+73] += values[285] * B[(i*84)+73];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b74 = _mm256_broadcast_sd(&B[(i*84)+74]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b74 = _mm_loaddup_pd(&B[(i*84)+74]);
#endif
__m128d c74_0 = _mm_load_sd(&C[(i*84)+10]);
__m128d a74_0 = _mm_load_sd(&values[286]);
#if defined(__SSE3__) && defined(__AVX__)
c74_0 = _mm_add_sd(c74_0, _mm_mul_sd(a74_0, _mm256_castpd256_pd128(b74)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c74_0 = _mm_add_sd(c74_0, _mm_mul_sd(a74_0, b74));
#endif
_mm_store_sd(&C[(i*84)+10], c74_0);
__m128d c74_1 = _mm_load_sd(&C[(i*84)+25]);
__m128d a74_1 = _mm_load_sd(&values[287]);
#if defined(__SSE3__) && defined(__AVX__)
c74_1 = _mm_add_sd(c74_1, _mm_mul_sd(a74_1, _mm256_castpd256_pd128(b74)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c74_1 = _mm_add_sd(c74_1, _mm_mul_sd(a74_1, b74));
#endif
_mm_store_sd(&C[(i*84)+25], c74_1);
__m128d c74_2 = _mm_load_sd(&C[(i*84)+46]);
__m128d a74_2 = _mm_load_sd(&values[288]);
#if defined(__SSE3__) && defined(__AVX__)
c74_2 = _mm_add_sd(c74_2, _mm_mul_sd(a74_2, _mm256_castpd256_pd128(b74)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c74_2 = _mm_add_sd(c74_2, _mm_mul_sd(a74_2, b74));
#endif
_mm_store_sd(&C[(i*84)+46], c74_2);
__m128d c74_3 = _mm_load_sd(&C[(i*84)+74]);
__m128d a74_3 = _mm_load_sd(&values[289]);
#if defined(__SSE3__) && defined(__AVX__)
c74_3 = _mm_add_sd(c74_3, _mm_mul_sd(a74_3, _mm256_castpd256_pd128(b74)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c74_3 = _mm_add_sd(c74_3, _mm_mul_sd(a74_3, b74));
#endif
_mm_store_sd(&C[(i*84)+74], c74_3);
#else
C[(i*84)+10] += values[286] * B[(i*84)+74];
C[(i*84)+25] += values[287] * B[(i*84)+74];
C[(i*84)+46] += values[288] * B[(i*84)+74];
C[(i*84)+74] += values[289] * B[(i*84)+74];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b75 = _mm256_broadcast_sd(&B[(i*84)+75]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b75 = _mm_loaddup_pd(&B[(i*84)+75]);
#endif
__m128d c75_0 = _mm_load_sd(&C[(i*84)+11]);
__m128d a75_0 = _mm_load_sd(&values[290]);
#if defined(__SSE3__) && defined(__AVX__)
c75_0 = _mm_add_sd(c75_0, _mm_mul_sd(a75_0, _mm256_castpd256_pd128(b75)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c75_0 = _mm_add_sd(c75_0, _mm_mul_sd(a75_0, b75));
#endif
_mm_store_sd(&C[(i*84)+11], c75_0);
__m128d c75_1 = _mm_load_sd(&C[(i*84)+26]);
__m128d a75_1 = _mm_load_sd(&values[291]);
#if defined(__SSE3__) && defined(__AVX__)
c75_1 = _mm_add_sd(c75_1, _mm_mul_sd(a75_1, _mm256_castpd256_pd128(b75)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c75_1 = _mm_add_sd(c75_1, _mm_mul_sd(a75_1, b75));
#endif
_mm_store_sd(&C[(i*84)+26], c75_1);
__m128d c75_2 = _mm_load_sd(&C[(i*84)+47]);
__m128d a75_2 = _mm_load_sd(&values[292]);
#if defined(__SSE3__) && defined(__AVX__)
c75_2 = _mm_add_sd(c75_2, _mm_mul_sd(a75_2, _mm256_castpd256_pd128(b75)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c75_2 = _mm_add_sd(c75_2, _mm_mul_sd(a75_2, b75));
#endif
_mm_store_sd(&C[(i*84)+47], c75_2);
__m128d c75_3 = _mm_load_sd(&C[(i*84)+75]);
__m128d a75_3 = _mm_load_sd(&values[293]);
#if defined(__SSE3__) && defined(__AVX__)
c75_3 = _mm_add_sd(c75_3, _mm_mul_sd(a75_3, _mm256_castpd256_pd128(b75)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c75_3 = _mm_add_sd(c75_3, _mm_mul_sd(a75_3, b75));
#endif
_mm_store_sd(&C[(i*84)+75], c75_3);
#else
C[(i*84)+11] += values[290] * B[(i*84)+75];
C[(i*84)+26] += values[291] * B[(i*84)+75];
C[(i*84)+47] += values[292] * B[(i*84)+75];
C[(i*84)+75] += values[293] * B[(i*84)+75];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b76 = _mm256_broadcast_sd(&B[(i*84)+76]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b76 = _mm_loaddup_pd(&B[(i*84)+76]);
#endif
__m128d c76_0 = _mm_load_sd(&C[(i*84)+12]);
__m128d a76_0 = _mm_load_sd(&values[294]);
#if defined(__SSE3__) && defined(__AVX__)
c76_0 = _mm_add_sd(c76_0, _mm_mul_sd(a76_0, _mm256_castpd256_pd128(b76)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c76_0 = _mm_add_sd(c76_0, _mm_mul_sd(a76_0, b76));
#endif
_mm_store_sd(&C[(i*84)+12], c76_0);
__m128d c76_1 = _mm_load_sd(&C[(i*84)+27]);
__m128d a76_1 = _mm_load_sd(&values[295]);
#if defined(__SSE3__) && defined(__AVX__)
c76_1 = _mm_add_sd(c76_1, _mm_mul_sd(a76_1, _mm256_castpd256_pd128(b76)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c76_1 = _mm_add_sd(c76_1, _mm_mul_sd(a76_1, b76));
#endif
_mm_store_sd(&C[(i*84)+27], c76_1);
__m128d c76_2 = _mm_load_sd(&C[(i*84)+48]);
__m128d a76_2 = _mm_load_sd(&values[296]);
#if defined(__SSE3__) && defined(__AVX__)
c76_2 = _mm_add_sd(c76_2, _mm_mul_sd(a76_2, _mm256_castpd256_pd128(b76)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c76_2 = _mm_add_sd(c76_2, _mm_mul_sd(a76_2, b76));
#endif
_mm_store_sd(&C[(i*84)+48], c76_2);
__m128d c76_3 = _mm_load_sd(&C[(i*84)+76]);
__m128d a76_3 = _mm_load_sd(&values[297]);
#if defined(__SSE3__) && defined(__AVX__)
c76_3 = _mm_add_sd(c76_3, _mm_mul_sd(a76_3, _mm256_castpd256_pd128(b76)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c76_3 = _mm_add_sd(c76_3, _mm_mul_sd(a76_3, b76));
#endif
_mm_store_sd(&C[(i*84)+76], c76_3);
#else
C[(i*84)+12] += values[294] * B[(i*84)+76];
C[(i*84)+27] += values[295] * B[(i*84)+76];
C[(i*84)+48] += values[296] * B[(i*84)+76];
C[(i*84)+76] += values[297] * B[(i*84)+76];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b77 = _mm256_broadcast_sd(&B[(i*84)+77]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b77 = _mm_loaddup_pd(&B[(i*84)+77]);
#endif
__m128d c77_0 = _mm_load_sd(&C[(i*84)+13]);
__m128d a77_0 = _mm_load_sd(&values[298]);
#if defined(__SSE3__) && defined(__AVX__)
c77_0 = _mm_add_sd(c77_0, _mm_mul_sd(a77_0, _mm256_castpd256_pd128(b77)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c77_0 = _mm_add_sd(c77_0, _mm_mul_sd(a77_0, b77));
#endif
_mm_store_sd(&C[(i*84)+13], c77_0);
__m128d c77_1 = _mm_load_sd(&C[(i*84)+28]);
__m128d a77_1 = _mm_load_sd(&values[299]);
#if defined(__SSE3__) && defined(__AVX__)
c77_1 = _mm_add_sd(c77_1, _mm_mul_sd(a77_1, _mm256_castpd256_pd128(b77)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c77_1 = _mm_add_sd(c77_1, _mm_mul_sd(a77_1, b77));
#endif
_mm_store_sd(&C[(i*84)+28], c77_1);
__m128d c77_2 = _mm_load_sd(&C[(i*84)+49]);
__m128d a77_2 = _mm_load_sd(&values[300]);
#if defined(__SSE3__) && defined(__AVX__)
c77_2 = _mm_add_sd(c77_2, _mm_mul_sd(a77_2, _mm256_castpd256_pd128(b77)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c77_2 = _mm_add_sd(c77_2, _mm_mul_sd(a77_2, b77));
#endif
_mm_store_sd(&C[(i*84)+49], c77_2);
__m128d c77_3 = _mm_load_sd(&C[(i*84)+77]);
__m128d a77_3 = _mm_load_sd(&values[301]);
#if defined(__SSE3__) && defined(__AVX__)
c77_3 = _mm_add_sd(c77_3, _mm_mul_sd(a77_3, _mm256_castpd256_pd128(b77)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c77_3 = _mm_add_sd(c77_3, _mm_mul_sd(a77_3, b77));
#endif
_mm_store_sd(&C[(i*84)+77], c77_3);
#else
C[(i*84)+13] += values[298] * B[(i*84)+77];
C[(i*84)+28] += values[299] * B[(i*84)+77];
C[(i*84)+49] += values[300] * B[(i*84)+77];
C[(i*84)+77] += values[301] * B[(i*84)+77];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b78 = _mm256_broadcast_sd(&B[(i*84)+78]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b78 = _mm_loaddup_pd(&B[(i*84)+78]);
#endif
__m128d c78_0 = _mm_load_sd(&C[(i*84)+4]);
__m128d a78_0 = _mm_load_sd(&values[302]);
#if defined(__SSE3__) && defined(__AVX__)
c78_0 = _mm_add_sd(c78_0, _mm_mul_sd(a78_0, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c78_0 = _mm_add_sd(c78_0, _mm_mul_sd(a78_0, b78));
#endif
_mm_store_sd(&C[(i*84)+4], c78_0);
__m128d c78_1 = _mm_load_sd(&C[(i*84)+14]);
__m128d a78_1 = _mm_load_sd(&values[303]);
#if defined(__SSE3__) && defined(__AVX__)
c78_1 = _mm_add_sd(c78_1, _mm_mul_sd(a78_1, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c78_1 = _mm_add_sd(c78_1, _mm_mul_sd(a78_1, b78));
#endif
_mm_store_sd(&C[(i*84)+14], c78_1);
__m128d c78_2 = _mm_load_sd(&C[(i*84)+29]);
__m128d a78_2 = _mm_load_sd(&values[304]);
#if defined(__SSE3__) && defined(__AVX__)
c78_2 = _mm_add_sd(c78_2, _mm_mul_sd(a78_2, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c78_2 = _mm_add_sd(c78_2, _mm_mul_sd(a78_2, b78));
#endif
_mm_store_sd(&C[(i*84)+29], c78_2);
__m128d c78_3 = _mm_load_sd(&C[(i*84)+50]);
__m128d a78_3 = _mm_load_sd(&values[305]);
#if defined(__SSE3__) && defined(__AVX__)
c78_3 = _mm_add_sd(c78_3, _mm_mul_sd(a78_3, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c78_3 = _mm_add_sd(c78_3, _mm_mul_sd(a78_3, b78));
#endif
_mm_store_sd(&C[(i*84)+50], c78_3);
__m128d c78_4 = _mm_load_sd(&C[(i*84)+78]);
__m128d a78_4 = _mm_load_sd(&values[306]);
#if defined(__SSE3__) && defined(__AVX__)
c78_4 = _mm_add_sd(c78_4, _mm_mul_sd(a78_4, _mm256_castpd256_pd128(b78)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c78_4 = _mm_add_sd(c78_4, _mm_mul_sd(a78_4, b78));
#endif
_mm_store_sd(&C[(i*84)+78], c78_4);
#else
C[(i*84)+4] += values[302] * B[(i*84)+78];
C[(i*84)+14] += values[303] * B[(i*84)+78];
C[(i*84)+29] += values[304] * B[(i*84)+78];
C[(i*84)+50] += values[305] * B[(i*84)+78];
C[(i*84)+78] += values[306] * B[(i*84)+78];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b79 = _mm256_broadcast_sd(&B[(i*84)+79]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b79 = _mm_loaddup_pd(&B[(i*84)+79]);
#endif
__m128d c79_0 = _mm_load_sd(&C[(i*84)+5]);
__m128d a79_0 = _mm_load_sd(&values[307]);
#if defined(__SSE3__) && defined(__AVX__)
c79_0 = _mm_add_sd(c79_0, _mm_mul_sd(a79_0, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c79_0 = _mm_add_sd(c79_0, _mm_mul_sd(a79_0, b79));
#endif
_mm_store_sd(&C[(i*84)+5], c79_0);
__m128d c79_1 = _mm_load_sd(&C[(i*84)+15]);
__m128d a79_1 = _mm_load_sd(&values[308]);
#if defined(__SSE3__) && defined(__AVX__)
c79_1 = _mm_add_sd(c79_1, _mm_mul_sd(a79_1, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c79_1 = _mm_add_sd(c79_1, _mm_mul_sd(a79_1, b79));
#endif
_mm_store_sd(&C[(i*84)+15], c79_1);
__m128d c79_2 = _mm_load_sd(&C[(i*84)+30]);
__m128d a79_2 = _mm_load_sd(&values[309]);
#if defined(__SSE3__) && defined(__AVX__)
c79_2 = _mm_add_sd(c79_2, _mm_mul_sd(a79_2, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c79_2 = _mm_add_sd(c79_2, _mm_mul_sd(a79_2, b79));
#endif
_mm_store_sd(&C[(i*84)+30], c79_2);
__m128d c79_3 = _mm_load_sd(&C[(i*84)+51]);
__m128d a79_3 = _mm_load_sd(&values[310]);
#if defined(__SSE3__) && defined(__AVX__)
c79_3 = _mm_add_sd(c79_3, _mm_mul_sd(a79_3, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c79_3 = _mm_add_sd(c79_3, _mm_mul_sd(a79_3, b79));
#endif
_mm_store_sd(&C[(i*84)+51], c79_3);
__m128d c79_4 = _mm_load_sd(&C[(i*84)+79]);
__m128d a79_4 = _mm_load_sd(&values[311]);
#if defined(__SSE3__) && defined(__AVX__)
c79_4 = _mm_add_sd(c79_4, _mm_mul_sd(a79_4, _mm256_castpd256_pd128(b79)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c79_4 = _mm_add_sd(c79_4, _mm_mul_sd(a79_4, b79));
#endif
_mm_store_sd(&C[(i*84)+79], c79_4);
#else
C[(i*84)+5] += values[307] * B[(i*84)+79];
C[(i*84)+15] += values[308] * B[(i*84)+79];
C[(i*84)+30] += values[309] * B[(i*84)+79];
C[(i*84)+51] += values[310] * B[(i*84)+79];
C[(i*84)+79] += values[311] * B[(i*84)+79];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b80 = _mm256_broadcast_sd(&B[(i*84)+80]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b80 = _mm_loaddup_pd(&B[(i*84)+80]);
#endif
__m128d c80_0 = _mm_load_sd(&C[(i*84)+6]);
__m128d a80_0 = _mm_load_sd(&values[312]);
#if defined(__SSE3__) && defined(__AVX__)
c80_0 = _mm_add_sd(c80_0, _mm_mul_sd(a80_0, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c80_0 = _mm_add_sd(c80_0, _mm_mul_sd(a80_0, b80));
#endif
_mm_store_sd(&C[(i*84)+6], c80_0);
__m128d c80_1 = _mm_load_sd(&C[(i*84)+16]);
__m128d a80_1 = _mm_load_sd(&values[313]);
#if defined(__SSE3__) && defined(__AVX__)
c80_1 = _mm_add_sd(c80_1, _mm_mul_sd(a80_1, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c80_1 = _mm_add_sd(c80_1, _mm_mul_sd(a80_1, b80));
#endif
_mm_store_sd(&C[(i*84)+16], c80_1);
__m128d c80_2 = _mm_load_sd(&C[(i*84)+31]);
__m128d a80_2 = _mm_load_sd(&values[314]);
#if defined(__SSE3__) && defined(__AVX__)
c80_2 = _mm_add_sd(c80_2, _mm_mul_sd(a80_2, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c80_2 = _mm_add_sd(c80_2, _mm_mul_sd(a80_2, b80));
#endif
_mm_store_sd(&C[(i*84)+31], c80_2);
__m128d c80_3 = _mm_load_sd(&C[(i*84)+52]);
__m128d a80_3 = _mm_load_sd(&values[315]);
#if defined(__SSE3__) && defined(__AVX__)
c80_3 = _mm_add_sd(c80_3, _mm_mul_sd(a80_3, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c80_3 = _mm_add_sd(c80_3, _mm_mul_sd(a80_3, b80));
#endif
_mm_store_sd(&C[(i*84)+52], c80_3);
__m128d c80_4 = _mm_load_sd(&C[(i*84)+80]);
__m128d a80_4 = _mm_load_sd(&values[316]);
#if defined(__SSE3__) && defined(__AVX__)
c80_4 = _mm_add_sd(c80_4, _mm_mul_sd(a80_4, _mm256_castpd256_pd128(b80)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c80_4 = _mm_add_sd(c80_4, _mm_mul_sd(a80_4, b80));
#endif
_mm_store_sd(&C[(i*84)+80], c80_4);
#else
C[(i*84)+6] += values[312] * B[(i*84)+80];
C[(i*84)+16] += values[313] * B[(i*84)+80];
C[(i*84)+31] += values[314] * B[(i*84)+80];
C[(i*84)+52] += values[315] * B[(i*84)+80];
C[(i*84)+80] += values[316] * B[(i*84)+80];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b81 = _mm256_broadcast_sd(&B[(i*84)+81]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b81 = _mm_loaddup_pd(&B[(i*84)+81]);
#endif
__m128d c81_0 = _mm_load_sd(&C[(i*84)+1]);
__m128d a81_0 = _mm_load_sd(&values[317]);
#if defined(__SSE3__) && defined(__AVX__)
c81_0 = _mm_add_sd(c81_0, _mm_mul_sd(a81_0, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c81_0 = _mm_add_sd(c81_0, _mm_mul_sd(a81_0, b81));
#endif
_mm_store_sd(&C[(i*84)+1], c81_0);
__m128d c81_1 = _mm_load_sd(&C[(i*84)+7]);
__m128d a81_1 = _mm_load_sd(&values[318]);
#if defined(__SSE3__) && defined(__AVX__)
c81_1 = _mm_add_sd(c81_1, _mm_mul_sd(a81_1, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c81_1 = _mm_add_sd(c81_1, _mm_mul_sd(a81_1, b81));
#endif
_mm_store_sd(&C[(i*84)+7], c81_1);
__m128d c81_2 = _mm_load_sd(&C[(i*84)+17]);
__m128d a81_2 = _mm_load_sd(&values[319]);
#if defined(__SSE3__) && defined(__AVX__)
c81_2 = _mm_add_sd(c81_2, _mm_mul_sd(a81_2, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c81_2 = _mm_add_sd(c81_2, _mm_mul_sd(a81_2, b81));
#endif
_mm_store_sd(&C[(i*84)+17], c81_2);
__m128d c81_3 = _mm_load_sd(&C[(i*84)+32]);
__m128d a81_3 = _mm_load_sd(&values[320]);
#if defined(__SSE3__) && defined(__AVX__)
c81_3 = _mm_add_sd(c81_3, _mm_mul_sd(a81_3, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c81_3 = _mm_add_sd(c81_3, _mm_mul_sd(a81_3, b81));
#endif
_mm_store_sd(&C[(i*84)+32], c81_3);
__m128d c81_4 = _mm_load_sd(&C[(i*84)+53]);
__m128d a81_4 = _mm_load_sd(&values[321]);
#if defined(__SSE3__) && defined(__AVX__)
c81_4 = _mm_add_sd(c81_4, _mm_mul_sd(a81_4, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c81_4 = _mm_add_sd(c81_4, _mm_mul_sd(a81_4, b81));
#endif
_mm_store_sd(&C[(i*84)+53], c81_4);
__m128d c81_5 = _mm_load_sd(&C[(i*84)+81]);
__m128d a81_5 = _mm_load_sd(&values[322]);
#if defined(__SSE3__) && defined(__AVX__)
c81_5 = _mm_add_sd(c81_5, _mm_mul_sd(a81_5, _mm256_castpd256_pd128(b81)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c81_5 = _mm_add_sd(c81_5, _mm_mul_sd(a81_5, b81));
#endif
_mm_store_sd(&C[(i*84)+81], c81_5);
#else
C[(i*84)+1] += values[317] * B[(i*84)+81];
C[(i*84)+7] += values[318] * B[(i*84)+81];
C[(i*84)+17] += values[319] * B[(i*84)+81];
C[(i*84)+32] += values[320] * B[(i*84)+81];
C[(i*84)+53] += values[321] * B[(i*84)+81];
C[(i*84)+81] += values[322] * B[(i*84)+81];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b82 = _mm256_broadcast_sd(&B[(i*84)+82]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b82 = _mm_loaddup_pd(&B[(i*84)+82]);
#endif
__m128d c82_0 = _mm_load_sd(&C[(i*84)+2]);
__m128d a82_0 = _mm_load_sd(&values[323]);
#if defined(__SSE3__) && defined(__AVX__)
c82_0 = _mm_add_sd(c82_0, _mm_mul_sd(a82_0, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c82_0 = _mm_add_sd(c82_0, _mm_mul_sd(a82_0, b82));
#endif
_mm_store_sd(&C[(i*84)+2], c82_0);
__m128d c82_1 = _mm_load_sd(&C[(i*84)+8]);
__m128d a82_1 = _mm_load_sd(&values[324]);
#if defined(__SSE3__) && defined(__AVX__)
c82_1 = _mm_add_sd(c82_1, _mm_mul_sd(a82_1, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c82_1 = _mm_add_sd(c82_1, _mm_mul_sd(a82_1, b82));
#endif
_mm_store_sd(&C[(i*84)+8], c82_1);
__m128d c82_2 = _mm_load_sd(&C[(i*84)+18]);
__m128d a82_2 = _mm_load_sd(&values[325]);
#if defined(__SSE3__) && defined(__AVX__)
c82_2 = _mm_add_sd(c82_2, _mm_mul_sd(a82_2, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c82_2 = _mm_add_sd(c82_2, _mm_mul_sd(a82_2, b82));
#endif
_mm_store_sd(&C[(i*84)+18], c82_2);
__m128d c82_3 = _mm_load_sd(&C[(i*84)+33]);
__m128d a82_3 = _mm_load_sd(&values[326]);
#if defined(__SSE3__) && defined(__AVX__)
c82_3 = _mm_add_sd(c82_3, _mm_mul_sd(a82_3, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c82_3 = _mm_add_sd(c82_3, _mm_mul_sd(a82_3, b82));
#endif
_mm_store_sd(&C[(i*84)+33], c82_3);
__m128d c82_4 = _mm_load_sd(&C[(i*84)+54]);
__m128d a82_4 = _mm_load_sd(&values[327]);
#if defined(__SSE3__) && defined(__AVX__)
c82_4 = _mm_add_sd(c82_4, _mm_mul_sd(a82_4, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c82_4 = _mm_add_sd(c82_4, _mm_mul_sd(a82_4, b82));
#endif
_mm_store_sd(&C[(i*84)+54], c82_4);
__m128d c82_5 = _mm_load_sd(&C[(i*84)+82]);
__m128d a82_5 = _mm_load_sd(&values[328]);
#if defined(__SSE3__) && defined(__AVX__)
c82_5 = _mm_add_sd(c82_5, _mm_mul_sd(a82_5, _mm256_castpd256_pd128(b82)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c82_5 = _mm_add_sd(c82_5, _mm_mul_sd(a82_5, b82));
#endif
_mm_store_sd(&C[(i*84)+82], c82_5);
#else
C[(i*84)+2] += values[323] * B[(i*84)+82];
C[(i*84)+8] += values[324] * B[(i*84)+82];
C[(i*84)+18] += values[325] * B[(i*84)+82];
C[(i*84)+33] += values[326] * B[(i*84)+82];
C[(i*84)+54] += values[327] * B[(i*84)+82];
C[(i*84)+82] += values[328] * B[(i*84)+82];
#endif
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
__m256d b83 = _mm256_broadcast_sd(&B[(i*84)+83]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
__m128d b83 = _mm_loaddup_pd(&B[(i*84)+83]);
#endif
__m128d c83_0 = _mm_load_sd(&C[(i*84)+0]);
__m128d a83_0 = _mm_load_sd(&values[329]);
#if defined(__SSE3__) && defined(__AVX__)
c83_0 = _mm_add_sd(c83_0, _mm_mul_sd(a83_0, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c83_0 = _mm_add_sd(c83_0, _mm_mul_sd(a83_0, b83));
#endif
_mm_store_sd(&C[(i*84)+0], c83_0);
__m128d c83_1 = _mm_load_sd(&C[(i*84)+3]);
__m128d a83_1 = _mm_load_sd(&values[330]);
#if defined(__SSE3__) && defined(__AVX__)
c83_1 = _mm_add_sd(c83_1, _mm_mul_sd(a83_1, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c83_1 = _mm_add_sd(c83_1, _mm_mul_sd(a83_1, b83));
#endif
_mm_store_sd(&C[(i*84)+3], c83_1);
__m128d c83_2 = _mm_load_sd(&C[(i*84)+9]);
__m128d a83_2 = _mm_load_sd(&values[331]);
#if defined(__SSE3__) && defined(__AVX__)
c83_2 = _mm_add_sd(c83_2, _mm_mul_sd(a83_2, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c83_2 = _mm_add_sd(c83_2, _mm_mul_sd(a83_2, b83));
#endif
_mm_store_sd(&C[(i*84)+9], c83_2);
__m128d c83_3 = _mm_load_sd(&C[(i*84)+19]);
__m128d a83_3 = _mm_load_sd(&values[332]);
#if defined(__SSE3__) && defined(__AVX__)
c83_3 = _mm_add_sd(c83_3, _mm_mul_sd(a83_3, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c83_3 = _mm_add_sd(c83_3, _mm_mul_sd(a83_3, b83));
#endif
_mm_store_sd(&C[(i*84)+19], c83_3);
__m128d c83_4 = _mm_load_sd(&C[(i*84)+34]);
__m128d a83_4 = _mm_load_sd(&values[333]);
#if defined(__SSE3__) && defined(__AVX__)
c83_4 = _mm_add_sd(c83_4, _mm_mul_sd(a83_4, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c83_4 = _mm_add_sd(c83_4, _mm_mul_sd(a83_4, b83));
#endif
_mm_store_sd(&C[(i*84)+34], c83_4);
__m128d c83_5 = _mm_load_sd(&C[(i*84)+55]);
__m128d a83_5 = _mm_load_sd(&values[334]);
#if defined(__SSE3__) && defined(__AVX__)
c83_5 = _mm_add_sd(c83_5, _mm_mul_sd(a83_5, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c83_5 = _mm_add_sd(c83_5, _mm_mul_sd(a83_5, b83));
#endif
_mm_store_sd(&C[(i*84)+55], c83_5);
__m128d c83_6 = _mm_load_sd(&C[(i*84)+83]);
__m128d a83_6 = _mm_load_sd(&values[335]);
#if defined(__SSE3__) && defined(__AVX__)
c83_6 = _mm_add_sd(c83_6, _mm_mul_sd(a83_6, _mm256_castpd256_pd128(b83)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
c83_6 = _mm_add_sd(c83_6, _mm_mul_sd(a83_6, b83));
#endif
_mm_store_sd(&C[(i*84)+83], c83_6);
#else
C[(i*84)+0] += values[329] * B[(i*84)+83];
C[(i*84)+3] += values[330] * B[(i*84)+83];
C[(i*84)+9] += values[331] * B[(i*84)+83];
C[(i*84)+19] += values[332] * B[(i*84)+83];
C[(i*84)+34] += values[333] * B[(i*84)+83];
C[(i*84)+55] += values[334] * B[(i*84)+83];
C[(i*84)+83] += values[335] * B[(i*84)+83];
#endif

}

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 6048;
#endif

}

void dsparse_starMatrix_m120_n9_k9_ldA120_ldBna8_ldC120_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
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
