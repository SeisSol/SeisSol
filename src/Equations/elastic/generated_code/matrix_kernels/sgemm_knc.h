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
// @date 2015-10-20 16:04:07.186005
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
#ifndef SGEMMKNCH
#define SGEMMKNCH

#if defined( __SSE3__) || defined(__MIC__)
#include <immintrin.h>
#endif

#include <cstddef>
#ifndef NDEBUG
extern long long libxsmm_num_total_flops;
#endif

void sgemm_m16_n9_k4_ldA64_ldB16_ldC16_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m48_n9_k9_ldA48_ldB9_ldC48_beta1_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m128_n9_k84_ldA128_ldB128_ldC128_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m16_n9_k10_ldA48_ldB16_ldC16_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m128_n9_k9_ldA128_ldB9_ldC128_beta1_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m96_n9_k9_ldA96_ldB9_ldC96_beta1_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m32_n9_k10_ldA32_ldB32_ldC32_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m96_n9_k120_ldA96_ldB128_ldC96_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m16_n9_k4_ldA48_ldB16_ldC16_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m64_n9_k35_ldA64_ldB64_ldC64_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m16_n9_k20_ldA48_ldB32_ldC16_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m96_n9_k56_ldA96_ldB96_ldC96_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m16_n9_k20_ldA96_ldB32_ldC16_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m32_n9_k35_ldA64_ldB48_ldC32_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m32_n9_k35_ldA48_ldB48_ldC32_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m48_n9_k56_ldA64_ldB64_ldC48_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m64_n9_k84_ldA64_ldB96_ldC64_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m16_n9_k1_ldA16_ldB16_ldC16_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m48_n9_k56_ldA96_ldB64_ldC48_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m16_n9_k20_ldA64_ldB32_ldC16_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m16_n9_k10_ldA96_ldB16_ldC16_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m16_n9_k4_ldA96_ldB16_ldC16_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m48_n9_k20_ldA48_ldB48_ldC48_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m16_n9_k4_ldA32_ldB16_ldC16_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m32_n9_k9_ldA32_ldB9_ldC32_beta1_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m32_n9_k35_ldA96_ldB48_ldC32_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m16_n9_k10_ldA64_ldB16_ldC16_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m16_n9_k20_ldA16_ldB32_ldC16_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m16_n9_k10_ldA32_ldB16_ldC16_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m64_n9_k9_ldA64_ldB9_ldC64_beta1_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m16_n9_k20_ldA32_ldB32_ldC16_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m48_n9_k56_ldA48_ldB64_ldC48_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m16_n9_k0_ldA16_ldB16_ldC16_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m32_n9_k35_ldA32_ldB48_ldC32_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void sgemm_m64_n9_k84_ldA96_ldB96_ldC64_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
#endif
