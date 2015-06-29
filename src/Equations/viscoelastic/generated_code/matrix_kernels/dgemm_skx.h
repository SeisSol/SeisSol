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
// @date 2015-06-29 16:03:31.307646
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
#ifndef DGEMMSKXH
#define DGEMMSKXH

#if defined( __SSE3__) || defined(__MIC__)
#include <immintrin.h>
#endif

#include <cstddef>
#ifndef NDEBUG
extern long long libxsmm_num_total_flops;
#endif

void dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m40_n27_k20_ldA40_ldB40_ldC40_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m8_n27_k1_ldA8_ldB8_ldC8_beta0_BL2viaC(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_BL2viaC(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m24_n27_k35_ldA40_ldB40_ldC24_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m8_n27_k10_ldA8_ldB16_ldC8_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m24_n27_k35_ldA88_ldB40_ldC24_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m56_n27_k27_ldA56_ldB27_ldC56_beta1_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m16_n27_k20_ldA24_ldB24_ldC16_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m24_n27_k10_ldA24_ldB24_ldC24_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m8_n27_k4_ldA16_ldB8_ldC8_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m40_n27_k56_ldA88_ldB56_ldC40_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m24_n27_k35_ldA24_ldB40_ldC24_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m24_n27_k35_ldA56_ldB40_ldC24_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m88_n27_k56_ldA88_ldB88_ldC88_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m8_n27_k27_ldA8_ldB27_ldC8_beta1_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m8_n27_k4_ldA56_ldB8_ldC8_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m16_n27_k20_ldA40_ldB24_ldC16_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_BL2viaC(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m8_n27_k10_ldA88_ldB16_ldC8_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m16_n27_k20_ldA56_ldB24_ldC16_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m16_n27_k20_ldA16_ldB24_ldC16_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m8_n27_k0_ldA8_ldB8_ldC8_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m40_n27_k56_ldA56_ldB56_ldC40_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m8_n27_k1_ldA8_ldB8_ldC8_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m8_n27_k10_ldA24_ldB16_ldC8_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_BL2viaC(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m56_n27_k84_ldA88_ldB88_ldC56_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m120_n27_k84_ldA120_ldB120_ldC120_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m8_n27_k10_ldA40_ldB16_ldC8_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m56_n27_k84_ldA56_ldB88_ldC56_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m8_n27_k4_ldA88_ldB8_ldC8_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m16_n27_k20_ldA88_ldB24_ldC16_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m8_n27_k4_ldA40_ldB8_ldC8_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m8_n27_k10_ldA16_ldB16_ldC8_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m56_n27_k35_ldA56_ldB56_ldC56_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m120_n27_k27_ldA120_ldB27_ldC120_beta1_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m88_n27_k120_ldA88_ldB120_ldC88_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_BL2viaC(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m88_n27_k27_ldA88_ldB27_ldC88_beta1_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m8_n27_k10_ldA56_ldB16_ldC8_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m40_n27_k56_ldA40_ldB56_ldC40_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m40_n27_k27_ldA40_ldB27_ldC40_beta1_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m24_n27_k27_ldA24_ldB27_ldC24_beta1_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
void dgemm_m8_n27_k4_ldA24_ldB8_ldC8_beta0_pfsigonly(const real *i_A, const real *i_B, real *io_C,const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );
#endif
