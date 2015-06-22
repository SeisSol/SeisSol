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
// @date 2015-05-09 22:18:08.091636
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
#ifndef SPARSEDNOARCHCPP
#define SPARSEDNOARCHCPP

#if defined( __SSE3__) || defined(__MIC__)
#include <immintrin.h>
#endif

#include <cstddef>
#ifndef NDEBUG
extern long long libxsmm_num_total_flops;
#endif

void dsparse_starMatrix_m1_n9_k9_ldA2_ldBna2_ldC2_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
for (unsigned int i = 0; i < 1; i += 1)
{
C[(i)+(0)] += A[(i)+(12)] * values[0];
C[(i)+(0)] += A[(i)+(14)] * values[1];
C[(i)+(0)] += A[(i)+(16)] * values[2];
C[(i)+(2)] += A[(i)+(12)] * values[3];
C[(i)+(2)] += A[(i)+(14)] * values[4];
C[(i)+(2)] += A[(i)+(16)] * values[5];
C[(i)+(4)] += A[(i)+(12)] * values[6];
C[(i)+(4)] += A[(i)+(14)] * values[7];
C[(i)+(4)] += A[(i)+(16)] * values[8];
C[(i)+(6)] += A[(i)+(12)] * values[9];
C[(i)+(6)] += A[(i)+(14)] * values[10];
C[(i)+(8)] += A[(i)+(14)] * values[11];
C[(i)+(8)] += A[(i)+(16)] * values[12];
C[(i)+(10)] += A[(i)+(12)] * values[13];
C[(i)+(10)] += A[(i)+(16)] * values[14];
C[(i)+(12)] += A[(i)+(0)] * values[15];
C[(i)+(12)] += A[(i)+(6)] * values[16];
C[(i)+(12)] += A[(i)+(10)] * values[17];
C[(i)+(14)] += A[(i)+(2)] * values[18];
C[(i)+(14)] += A[(i)+(6)] * values[19];
C[(i)+(14)] += A[(i)+(8)] * values[20];
C[(i)+(16)] += A[(i)+(4)] * values[21];
C[(i)+(16)] += A[(i)+(8)] * values[22];
C[(i)+(16)] += A[(i)+(10)] * values[23];
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

void dsparse_starMatrix_m1_n9_k9_ldA2_ldBna3_ldC2_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
for (unsigned int i = 0; i < 1; i += 1)
{
C[(i)+(0)] += A[(i)+(12)] * values[0];
C[(i)+(0)] += A[(i)+(14)] * values[1];
C[(i)+(0)] += A[(i)+(16)] * values[2];
C[(i)+(2)] += A[(i)+(12)] * values[3];
C[(i)+(2)] += A[(i)+(14)] * values[4];
C[(i)+(2)] += A[(i)+(16)] * values[5];
C[(i)+(4)] += A[(i)+(12)] * values[6];
C[(i)+(4)] += A[(i)+(14)] * values[7];
C[(i)+(4)] += A[(i)+(16)] * values[8];
C[(i)+(6)] += A[(i)+(12)] * values[9];
C[(i)+(6)] += A[(i)+(14)] * values[10];
C[(i)+(8)] += A[(i)+(14)] * values[11];
C[(i)+(8)] += A[(i)+(16)] * values[12];
C[(i)+(10)] += A[(i)+(12)] * values[13];
C[(i)+(10)] += A[(i)+(16)] * values[14];
C[(i)+(12)] += A[(i)+(0)] * values[15];
C[(i)+(12)] += A[(i)+(6)] * values[16];
C[(i)+(12)] += A[(i)+(10)] * values[17];
C[(i)+(14)] += A[(i)+(2)] * values[18];
C[(i)+(14)] += A[(i)+(6)] * values[19];
C[(i)+(14)] += A[(i)+(8)] * values[20];
C[(i)+(16)] += A[(i)+(4)] * values[21];
C[(i)+(16)] += A[(i)+(8)] * values[22];
C[(i)+(16)] += A[(i)+(10)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 48;
#endif

}

void dsparse_starMatrix_m10_n9_k9_ldA10_ldBna4_ldC10_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
#pragma vector aligned
for (unsigned int i = 0; i < 10; i += 1)
{
C[(i)+(0)] += A[(i)+(60)] * values[0];
C[(i)+(0)] += A[(i)+(70)] * values[1];
C[(i)+(0)] += A[(i)+(80)] * values[2];
C[(i)+(10)] += A[(i)+(60)] * values[3];
C[(i)+(10)] += A[(i)+(70)] * values[4];
C[(i)+(10)] += A[(i)+(80)] * values[5];
C[(i)+(20)] += A[(i)+(60)] * values[6];
C[(i)+(20)] += A[(i)+(70)] * values[7];
C[(i)+(20)] += A[(i)+(80)] * values[8];
C[(i)+(30)] += A[(i)+(60)] * values[9];
C[(i)+(30)] += A[(i)+(70)] * values[10];
C[(i)+(40)] += A[(i)+(70)] * values[11];
C[(i)+(40)] += A[(i)+(80)] * values[12];
C[(i)+(50)] += A[(i)+(60)] * values[13];
C[(i)+(50)] += A[(i)+(80)] * values[14];
C[(i)+(60)] += A[(i)+(0)] * values[15];
C[(i)+(60)] += A[(i)+(30)] * values[16];
C[(i)+(60)] += A[(i)+(50)] * values[17];
C[(i)+(70)] += A[(i)+(10)] * values[18];
C[(i)+(70)] += A[(i)+(30)] * values[19];
C[(i)+(70)] += A[(i)+(40)] * values[20];
C[(i)+(80)] += A[(i)+(20)] * values[21];
C[(i)+(80)] += A[(i)+(40)] * values[22];
C[(i)+(80)] += A[(i)+(50)] * values[23];
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

void dsparse_starMatrix_m1_n9_k9_ldA2_ldBna4_ldC2_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
for (unsigned int i = 0; i < 1; i += 1)
{
C[(i)+(0)] += A[(i)+(12)] * values[0];
C[(i)+(0)] += A[(i)+(14)] * values[1];
C[(i)+(0)] += A[(i)+(16)] * values[2];
C[(i)+(2)] += A[(i)+(12)] * values[3];
C[(i)+(2)] += A[(i)+(14)] * values[4];
C[(i)+(2)] += A[(i)+(16)] * values[5];
C[(i)+(4)] += A[(i)+(12)] * values[6];
C[(i)+(4)] += A[(i)+(14)] * values[7];
C[(i)+(4)] += A[(i)+(16)] * values[8];
C[(i)+(6)] += A[(i)+(12)] * values[9];
C[(i)+(6)] += A[(i)+(14)] * values[10];
C[(i)+(8)] += A[(i)+(14)] * values[11];
C[(i)+(8)] += A[(i)+(16)] * values[12];
C[(i)+(10)] += A[(i)+(12)] * values[13];
C[(i)+(10)] += A[(i)+(16)] * values[14];
C[(i)+(12)] += A[(i)+(0)] * values[15];
C[(i)+(12)] += A[(i)+(6)] * values[16];
C[(i)+(12)] += A[(i)+(10)] * values[17];
C[(i)+(14)] += A[(i)+(2)] * values[18];
C[(i)+(14)] += A[(i)+(6)] * values[19];
C[(i)+(14)] += A[(i)+(8)] * values[20];
C[(i)+(16)] += A[(i)+(4)] * values[21];
C[(i)+(16)] += A[(i)+(8)] * values[22];
C[(i)+(16)] += A[(i)+(10)] * values[23];
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

void dsparse_starMatrix_m10_n9_k9_ldA10_ldBna5_ldC10_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
#pragma vector aligned
for (unsigned int i = 0; i < 10; i += 1)
{
C[(i)+(0)] += A[(i)+(60)] * values[0];
C[(i)+(0)] += A[(i)+(70)] * values[1];
C[(i)+(0)] += A[(i)+(80)] * values[2];
C[(i)+(10)] += A[(i)+(60)] * values[3];
C[(i)+(10)] += A[(i)+(70)] * values[4];
C[(i)+(10)] += A[(i)+(80)] * values[5];
C[(i)+(20)] += A[(i)+(60)] * values[6];
C[(i)+(20)] += A[(i)+(70)] * values[7];
C[(i)+(20)] += A[(i)+(80)] * values[8];
C[(i)+(30)] += A[(i)+(60)] * values[9];
C[(i)+(30)] += A[(i)+(70)] * values[10];
C[(i)+(40)] += A[(i)+(70)] * values[11];
C[(i)+(40)] += A[(i)+(80)] * values[12];
C[(i)+(50)] += A[(i)+(60)] * values[13];
C[(i)+(50)] += A[(i)+(80)] * values[14];
C[(i)+(60)] += A[(i)+(0)] * values[15];
C[(i)+(60)] += A[(i)+(30)] * values[16];
C[(i)+(60)] += A[(i)+(50)] * values[17];
C[(i)+(70)] += A[(i)+(10)] * values[18];
C[(i)+(70)] += A[(i)+(30)] * values[19];
C[(i)+(70)] += A[(i)+(40)] * values[20];
C[(i)+(80)] += A[(i)+(20)] * values[21];
C[(i)+(80)] += A[(i)+(40)] * values[22];
C[(i)+(80)] += A[(i)+(50)] * values[23];
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

void dsparse_starMatrix_m1_n9_k9_ldA2_ldBna5_ldC2_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
for (unsigned int i = 0; i < 1; i += 1)
{
C[(i)+(0)] += A[(i)+(12)] * values[0];
C[(i)+(0)] += A[(i)+(14)] * values[1];
C[(i)+(0)] += A[(i)+(16)] * values[2];
C[(i)+(2)] += A[(i)+(12)] * values[3];
C[(i)+(2)] += A[(i)+(14)] * values[4];
C[(i)+(2)] += A[(i)+(16)] * values[5];
C[(i)+(4)] += A[(i)+(12)] * values[6];
C[(i)+(4)] += A[(i)+(14)] * values[7];
C[(i)+(4)] += A[(i)+(16)] * values[8];
C[(i)+(6)] += A[(i)+(12)] * values[9];
C[(i)+(6)] += A[(i)+(14)] * values[10];
C[(i)+(8)] += A[(i)+(14)] * values[11];
C[(i)+(8)] += A[(i)+(16)] * values[12];
C[(i)+(10)] += A[(i)+(12)] * values[13];
C[(i)+(10)] += A[(i)+(16)] * values[14];
C[(i)+(12)] += A[(i)+(0)] * values[15];
C[(i)+(12)] += A[(i)+(6)] * values[16];
C[(i)+(12)] += A[(i)+(10)] * values[17];
C[(i)+(14)] += A[(i)+(2)] * values[18];
C[(i)+(14)] += A[(i)+(6)] * values[19];
C[(i)+(14)] += A[(i)+(8)] * values[20];
C[(i)+(16)] += A[(i)+(4)] * values[21];
C[(i)+(16)] += A[(i)+(8)] * values[22];
C[(i)+(16)] += A[(i)+(10)] * values[23];
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

void dsparse_starMatrix_m10_n9_k9_ldA10_ldBna6_ldC10_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
#pragma vector aligned
for (unsigned int i = 0; i < 10; i += 1)
{
C[(i)+(0)] += A[(i)+(60)] * values[0];
C[(i)+(0)] += A[(i)+(70)] * values[1];
C[(i)+(0)] += A[(i)+(80)] * values[2];
C[(i)+(10)] += A[(i)+(60)] * values[3];
C[(i)+(10)] += A[(i)+(70)] * values[4];
C[(i)+(10)] += A[(i)+(80)] * values[5];
C[(i)+(20)] += A[(i)+(60)] * values[6];
C[(i)+(20)] += A[(i)+(70)] * values[7];
C[(i)+(20)] += A[(i)+(80)] * values[8];
C[(i)+(30)] += A[(i)+(60)] * values[9];
C[(i)+(30)] += A[(i)+(70)] * values[10];
C[(i)+(40)] += A[(i)+(70)] * values[11];
C[(i)+(40)] += A[(i)+(80)] * values[12];
C[(i)+(50)] += A[(i)+(60)] * values[13];
C[(i)+(50)] += A[(i)+(80)] * values[14];
C[(i)+(60)] += A[(i)+(0)] * values[15];
C[(i)+(60)] += A[(i)+(30)] * values[16];
C[(i)+(60)] += A[(i)+(50)] * values[17];
C[(i)+(70)] += A[(i)+(10)] * values[18];
C[(i)+(70)] += A[(i)+(30)] * values[19];
C[(i)+(70)] += A[(i)+(40)] * values[20];
C[(i)+(80)] += A[(i)+(20)] * values[21];
C[(i)+(80)] += A[(i)+(40)] * values[22];
C[(i)+(80)] += A[(i)+(50)] * values[23];
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

void dsparse_starMatrix_m1_n9_k9_ldA2_ldBna6_ldC2_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
for (unsigned int i = 0; i < 1; i += 1)
{
C[(i)+(0)] += A[(i)+(12)] * values[0];
C[(i)+(0)] += A[(i)+(14)] * values[1];
C[(i)+(0)] += A[(i)+(16)] * values[2];
C[(i)+(2)] += A[(i)+(12)] * values[3];
C[(i)+(2)] += A[(i)+(14)] * values[4];
C[(i)+(2)] += A[(i)+(16)] * values[5];
C[(i)+(4)] += A[(i)+(12)] * values[6];
C[(i)+(4)] += A[(i)+(14)] * values[7];
C[(i)+(4)] += A[(i)+(16)] * values[8];
C[(i)+(6)] += A[(i)+(12)] * values[9];
C[(i)+(6)] += A[(i)+(14)] * values[10];
C[(i)+(8)] += A[(i)+(14)] * values[11];
C[(i)+(8)] += A[(i)+(16)] * values[12];
C[(i)+(10)] += A[(i)+(12)] * values[13];
C[(i)+(10)] += A[(i)+(16)] * values[14];
C[(i)+(12)] += A[(i)+(0)] * values[15];
C[(i)+(12)] += A[(i)+(6)] * values[16];
C[(i)+(12)] += A[(i)+(10)] * values[17];
C[(i)+(14)] += A[(i)+(2)] * values[18];
C[(i)+(14)] += A[(i)+(6)] * values[19];
C[(i)+(14)] += A[(i)+(8)] * values[20];
C[(i)+(16)] += A[(i)+(4)] * values[21];
C[(i)+(16)] += A[(i)+(8)] * values[22];
C[(i)+(16)] += A[(i)+(10)] * values[23];
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

void dsparse_starMatrix_m10_n9_k9_ldA10_ldBna7_ldC10_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
#pragma vector aligned
for (unsigned int i = 0; i < 10; i += 1)
{
C[(i)+(0)] += A[(i)+(60)] * values[0];
C[(i)+(0)] += A[(i)+(70)] * values[1];
C[(i)+(0)] += A[(i)+(80)] * values[2];
C[(i)+(10)] += A[(i)+(60)] * values[3];
C[(i)+(10)] += A[(i)+(70)] * values[4];
C[(i)+(10)] += A[(i)+(80)] * values[5];
C[(i)+(20)] += A[(i)+(60)] * values[6];
C[(i)+(20)] += A[(i)+(70)] * values[7];
C[(i)+(20)] += A[(i)+(80)] * values[8];
C[(i)+(30)] += A[(i)+(60)] * values[9];
C[(i)+(30)] += A[(i)+(70)] * values[10];
C[(i)+(40)] += A[(i)+(70)] * values[11];
C[(i)+(40)] += A[(i)+(80)] * values[12];
C[(i)+(50)] += A[(i)+(60)] * values[13];
C[(i)+(50)] += A[(i)+(80)] * values[14];
C[(i)+(60)] += A[(i)+(0)] * values[15];
C[(i)+(60)] += A[(i)+(30)] * values[16];
C[(i)+(60)] += A[(i)+(50)] * values[17];
C[(i)+(70)] += A[(i)+(10)] * values[18];
C[(i)+(70)] += A[(i)+(30)] * values[19];
C[(i)+(70)] += A[(i)+(40)] * values[20];
C[(i)+(80)] += A[(i)+(20)] * values[21];
C[(i)+(80)] += A[(i)+(40)] * values[22];
C[(i)+(80)] += A[(i)+(50)] * values[23];
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

void dsparse_starMatrix_m1_n9_k9_ldA2_ldBna7_ldC2_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
for (unsigned int i = 0; i < 1; i += 1)
{
C[(i)+(0)] += A[(i)+(12)] * values[0];
C[(i)+(0)] += A[(i)+(14)] * values[1];
C[(i)+(0)] += A[(i)+(16)] * values[2];
C[(i)+(2)] += A[(i)+(12)] * values[3];
C[(i)+(2)] += A[(i)+(14)] * values[4];
C[(i)+(2)] += A[(i)+(16)] * values[5];
C[(i)+(4)] += A[(i)+(12)] * values[6];
C[(i)+(4)] += A[(i)+(14)] * values[7];
C[(i)+(4)] += A[(i)+(16)] * values[8];
C[(i)+(6)] += A[(i)+(12)] * values[9];
C[(i)+(6)] += A[(i)+(14)] * values[10];
C[(i)+(8)] += A[(i)+(14)] * values[11];
C[(i)+(8)] += A[(i)+(16)] * values[12];
C[(i)+(10)] += A[(i)+(12)] * values[13];
C[(i)+(10)] += A[(i)+(16)] * values[14];
C[(i)+(12)] += A[(i)+(0)] * values[15];
C[(i)+(12)] += A[(i)+(6)] * values[16];
C[(i)+(12)] += A[(i)+(10)] * values[17];
C[(i)+(14)] += A[(i)+(2)] * values[18];
C[(i)+(14)] += A[(i)+(6)] * values[19];
C[(i)+(14)] += A[(i)+(8)] * values[20];
C[(i)+(16)] += A[(i)+(4)] * values[21];
C[(i)+(16)] += A[(i)+(8)] * values[22];
C[(i)+(16)] += A[(i)+(10)] * values[23];
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

void dsparse_starMatrix_m10_n9_k9_ldA10_ldBna8_ldC10_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
#pragma vector aligned
for (unsigned int i = 0; i < 10; i += 1)
{
C[(i)+(0)] += A[(i)+(60)] * values[0];
C[(i)+(0)] += A[(i)+(70)] * values[1];
C[(i)+(0)] += A[(i)+(80)] * values[2];
C[(i)+(10)] += A[(i)+(60)] * values[3];
C[(i)+(10)] += A[(i)+(70)] * values[4];
C[(i)+(10)] += A[(i)+(80)] * values[5];
C[(i)+(20)] += A[(i)+(60)] * values[6];
C[(i)+(20)] += A[(i)+(70)] * values[7];
C[(i)+(20)] += A[(i)+(80)] * values[8];
C[(i)+(30)] += A[(i)+(60)] * values[9];
C[(i)+(30)] += A[(i)+(70)] * values[10];
C[(i)+(40)] += A[(i)+(70)] * values[11];
C[(i)+(40)] += A[(i)+(80)] * values[12];
C[(i)+(50)] += A[(i)+(60)] * values[13];
C[(i)+(50)] += A[(i)+(80)] * values[14];
C[(i)+(60)] += A[(i)+(0)] * values[15];
C[(i)+(60)] += A[(i)+(30)] * values[16];
C[(i)+(60)] += A[(i)+(50)] * values[17];
C[(i)+(70)] += A[(i)+(10)] * values[18];
C[(i)+(70)] += A[(i)+(30)] * values[19];
C[(i)+(70)] += A[(i)+(40)] * values[20];
C[(i)+(80)] += A[(i)+(20)] * values[21];
C[(i)+(80)] += A[(i)+(40)] * values[22];
C[(i)+(80)] += A[(i)+(50)] * values[23];
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

void dsparse_starMatrix_m1_n9_k9_ldA2_ldBna8_ldC2_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
for (unsigned int i = 0; i < 1; i += 1)
{
C[(i)+(0)] += A[(i)+(12)] * values[0];
C[(i)+(0)] += A[(i)+(14)] * values[1];
C[(i)+(0)] += A[(i)+(16)] * values[2];
C[(i)+(2)] += A[(i)+(12)] * values[3];
C[(i)+(2)] += A[(i)+(14)] * values[4];
C[(i)+(2)] += A[(i)+(16)] * values[5];
C[(i)+(4)] += A[(i)+(12)] * values[6];
C[(i)+(4)] += A[(i)+(14)] * values[7];
C[(i)+(4)] += A[(i)+(16)] * values[8];
C[(i)+(6)] += A[(i)+(12)] * values[9];
C[(i)+(6)] += A[(i)+(14)] * values[10];
C[(i)+(8)] += A[(i)+(14)] * values[11];
C[(i)+(8)] += A[(i)+(16)] * values[12];
C[(i)+(10)] += A[(i)+(12)] * values[13];
C[(i)+(10)] += A[(i)+(16)] * values[14];
C[(i)+(12)] += A[(i)+(0)] * values[15];
C[(i)+(12)] += A[(i)+(6)] * values[16];
C[(i)+(12)] += A[(i)+(10)] * values[17];
C[(i)+(14)] += A[(i)+(2)] * values[18];
C[(i)+(14)] += A[(i)+(6)] * values[19];
C[(i)+(14)] += A[(i)+(8)] * values[20];
C[(i)+(16)] += A[(i)+(4)] * values[21];
C[(i)+(16)] += A[(i)+(8)] * values[22];
C[(i)+(16)] += A[(i)+(10)] * values[23];
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

void dsparse_starMatrix_m10_n9_k9_ldA10_ldBna3_ldC10_beta1_pfsigonly(const double* A, const double* values, double* C, const double* A_prefetch = NULL, const double* B_prefetch = NULL, const double* C_prefetch = NULL)
{
#pragma simd vectorlength(8)
#pragma vector aligned
for (unsigned int i = 0; i < 10; i += 1)
{
C[(i)+(0)] += A[(i)+(60)] * values[0];
C[(i)+(0)] += A[(i)+(70)] * values[1];
C[(i)+(0)] += A[(i)+(80)] * values[2];
C[(i)+(10)] += A[(i)+(60)] * values[3];
C[(i)+(10)] += A[(i)+(70)] * values[4];
C[(i)+(10)] += A[(i)+(80)] * values[5];
C[(i)+(20)] += A[(i)+(60)] * values[6];
C[(i)+(20)] += A[(i)+(70)] * values[7];
C[(i)+(20)] += A[(i)+(80)] * values[8];
C[(i)+(30)] += A[(i)+(60)] * values[9];
C[(i)+(30)] += A[(i)+(70)] * values[10];
C[(i)+(40)] += A[(i)+(70)] * values[11];
C[(i)+(40)] += A[(i)+(80)] * values[12];
C[(i)+(50)] += A[(i)+(60)] * values[13];
C[(i)+(50)] += A[(i)+(80)] * values[14];
C[(i)+(60)] += A[(i)+(0)] * values[15];
C[(i)+(60)] += A[(i)+(30)] * values[16];
C[(i)+(60)] += A[(i)+(50)] * values[17];
C[(i)+(70)] += A[(i)+(10)] * values[18];
C[(i)+(70)] += A[(i)+(30)] * values[19];
C[(i)+(70)] += A[(i)+(40)] * values[20];
C[(i)+(80)] += A[(i)+(20)] * values[21];
C[(i)+(80)] += A[(i)+(40)] * values[22];
C[(i)+(80)] += A[(i)+(50)] * values[23];
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 480;
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
