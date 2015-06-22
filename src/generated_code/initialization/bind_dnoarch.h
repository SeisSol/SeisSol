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
// @date 2015-05-09 22:18:40.240054
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
#if ALIGNMENT!=16
#error alignment-architecture mismatch
#endif

#if CONVERGENCE_ORDER==2

#ifdef TIME_KERNEL
m_nonZeroFlops[0] = 18;
m_matrixKernels[0] = dgemm_m2_n9_k4_ldA2_ldB4_ldC2_beta0_pfsigonly;
m_hardwareFlops[0] = 144;
m_nonZeroFlops[1] = 36;
m_matrixKernels[1] = dgemm_m2_n9_k4_ldA2_ldB4_ldC2_beta0_pfsigonly;
m_hardwareFlops[1] = 144;
m_nonZeroFlops[2] = 54;
m_matrixKernels[2] = dgemm_m2_n9_k4_ldA2_ldB4_ldC2_beta0_pfsigonly;
m_hardwareFlops[2] = 144;
m_nonZeroFlops[3] = 48;
m_hardwareFlops[3] = 48;
m_matrixKernels[3] = dsparse_starMatrix_m1_n9_k9_ldA2_ldBna2_ldC2_beta1_pfsigonly;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 18;
m_matrixKernels[0] = dgemm_m4_n9_k1_ldA4_ldB4_ldC4_beta0_pfsigonly;
m_hardwareFlops[0] = 72;
m_nonZeroFlops[1] = 36;
m_matrixKernels[1] = dgemm_m4_n9_k1_ldA4_ldB4_ldC4_beta0_pfsigonly;
m_hardwareFlops[1] = 72;
m_nonZeroFlops[2] = 54;
m_matrixKernels[2] = dgemm_m4_n9_k1_ldA4_ldB4_ldC4_beta0_pfsigonly;
m_hardwareFlops[2] = 72;
m_nonZeroFlops[3] = 192;
m_hardwareFlops[3] = 192;
m_matrixKernels[3] = dsparse_starMatrix_m4_n9_k9_ldA4_ldBna2_ldC4_beta1_pfsigonly;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 108;
m_hardwareFlops[0] = 288;
m_matrixKernels[0] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
m_nonZeroFlops[1] = 144;
m_hardwareFlops[1] = 288;
m_matrixKernels[1] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
m_nonZeroFlops[2] = 180;
m_hardwareFlops[2] = 288;
m_matrixKernels[2] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
m_nonZeroFlops[3] = 180;
m_hardwareFlops[3] = 288;
m_matrixKernels[3] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
m_nonZeroFlops[4] = 144;
m_hardwareFlops[4] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[4] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 144;
m_hardwareFlops[5] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[5] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 108;
m_hardwareFlops[6] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[6] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 162;
m_hardwareFlops[7] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[7] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 216;
m_hardwareFlops[8] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[8] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 216;
m_hardwareFlops[9] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[9] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 252;
m_hardwareFlops[10] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[10] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 234;
m_hardwareFlops[11] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[11] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 234;
m_hardwareFlops[12] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[12] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 234;
m_hardwareFlops[13] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[13] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 234;
m_hardwareFlops[14] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[14] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 252;
m_hardwareFlops[15] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[15] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 162;
m_hardwareFlops[16] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[16] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 216;
m_hardwareFlops[17] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[17] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 216;
m_hardwareFlops[18] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[18] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 252;
m_hardwareFlops[19] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[19] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 144;
m_hardwareFlops[20] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[20] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 252;
m_hardwareFlops[21] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[21] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 198;
m_hardwareFlops[22] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[22] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 252;
m_hardwareFlops[23] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[23] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 252;
m_hardwareFlops[24] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[24] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 252;
m_hardwareFlops[25] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[25] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 198;
m_hardwareFlops[26] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[26] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 252;
m_hardwareFlops[27] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[27] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 252;
m_hardwareFlops[28] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[28] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 234;
m_hardwareFlops[29] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[29] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 234;
m_hardwareFlops[30] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[30] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 198;
m_hardwareFlops[31] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[31] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 252;
m_hardwareFlops[32] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[32] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 252;
m_hardwareFlops[33] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[33] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 216;
m_hardwareFlops[34] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[34] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 288;
m_hardwareFlops[35] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[35] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 216;
m_hardwareFlops[36] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[36] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 180;
m_hardwareFlops[37] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[37] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 270;
m_hardwareFlops[38] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[38] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 270;
m_hardwareFlops[39] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[39] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 234;
m_hardwareFlops[40] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[40] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 234;
m_hardwareFlops[41] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[41] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 252;
m_hardwareFlops[42] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[42] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 252;
m_hardwareFlops[43] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[43] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 198;
m_hardwareFlops[44] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[44] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 252;
m_hardwareFlops[45] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[45] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 180;
m_hardwareFlops[46] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[46] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 270;
m_hardwareFlops[47] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[47] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 270;
m_hardwareFlops[48] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[48] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 216;
m_hardwareFlops[49] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[49] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 216;
m_hardwareFlops[50] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[50] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 288;
m_hardwareFlops[51] = 288;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#else
m_matrixKernels[51] = dgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 648;
m_hardwareFlops[52] = 648;
m_matrixKernels[52] = dgemm_m4_n9_k9_ldA4_ldB9_ldC4_beta1_pfsigonly;
m_nonZeroFlops[53] = 648;
m_hardwareFlops[53] = 648;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = dgemm_m4_n9_k9_ldA4_ldB9_ldC4_beta1_pfsigonly;
#else
m_matrixKernels[53] = dgemm_m4_n9_k9_ldA4_ldB9_ldC4_beta1_pfsigonly;
#endif
#endif

#ifdef SPARSE_SWITCH
m_sparseSwitch[0] = -1; 
m_sparseSwitch[1] = -1; 
m_sparseSwitch[2] = -1; 
m_sparseSwitch[3] = -1; 
m_sparseSwitch[4] = -1; 
m_sparseSwitch[5] = -1; 
m_sparseSwitch[6] = -1; 
m_sparseSwitch[7] = -1; 
m_sparseSwitch[8] = -1; 
m_sparseSwitch[9] = -1; 
m_sparseSwitch[10] = -1; 
m_sparseSwitch[11] = -1; 
m_sparseSwitch[12] = -1; 
m_sparseSwitch[13] = -1; 
m_sparseSwitch[14] = -1; 
m_sparseSwitch[15] = -1; 
m_sparseSwitch[16] = -1; 
m_sparseSwitch[17] = -1; 
m_sparseSwitch[18] = -1; 
m_sparseSwitch[19] = -1; 
m_sparseSwitch[20] = -1; 
m_sparseSwitch[21] = -1; 
m_sparseSwitch[22] = -1; 
m_sparseSwitch[23] = -1; 
m_sparseSwitch[24] = -1; 
m_sparseSwitch[25] = -1; 
m_sparseSwitch[26] = -1; 
m_sparseSwitch[27] = -1; 
m_sparseSwitch[28] = -1; 
m_sparseSwitch[29] = -1; 
m_sparseSwitch[30] = -1; 
m_sparseSwitch[31] = -1; 
m_sparseSwitch[32] = -1; 
m_sparseSwitch[33] = -1; 
m_sparseSwitch[34] = -1; 
m_sparseSwitch[35] = -1; 
m_sparseSwitch[36] = -1; 
m_sparseSwitch[37] = -1; 
m_sparseSwitch[38] = -1; 
m_sparseSwitch[39] = -1; 
m_sparseSwitch[40] = -1; 
m_sparseSwitch[41] = -1; 
m_sparseSwitch[42] = -1; 
m_sparseSwitch[43] = -1; 
m_sparseSwitch[44] = -1; 
m_sparseSwitch[45] = -1; 
m_sparseSwitch[46] = -1; 
m_sparseSwitch[47] = -1; 
m_sparseSwitch[48] = -1; 
m_sparseSwitch[49] = -1; 
m_sparseSwitch[50] = -1; 
m_sparseSwitch[51] = -1; 
m_sparseSwitch[52] = -1; 
m_sparseSwitch[53] = -1; 
m_sparseSwitch[54] = -1; 
m_sparseSwitch[55] = -1; 
m_sparseSwitch[56] = -1; 
m_sparseSwitch[57] = -1; 
m_sparseSwitch[58] = -1; 
m_sparseSwitch[59] = 24; 
#endif

#define STAR_NNZ 24

#endif

#if CONVERGENCE_ORDER==3

#ifdef TIME_KERNEL
m_nonZeroFlops[0] = 126;
m_matrixKernels[0] = dgemm_m4_n9_k10_ldA4_ldB10_ldC4_beta0_pfsigonly;
m_hardwareFlops[0] = 720;
m_nonZeroFlops[1] = 306;
m_matrixKernels[1] = dgemm_m4_n9_k10_ldA4_ldB10_ldC4_beta0_pfsigonly;
m_hardwareFlops[1] = 720;
m_nonZeroFlops[2] = 396;
m_matrixKernels[2] = dgemm_m4_n9_k10_ldA4_ldB10_ldC4_beta0_pfsigonly;
m_hardwareFlops[2] = 720;
m_nonZeroFlops[3] = 192;
m_hardwareFlops[3] = 192;
m_matrixKernels[3] = dsparse_starMatrix_m4_n9_k9_ldA4_ldBna3_ldC4_beta1_pfsigonly;
m_nonZeroFlops[4] = 18;
m_matrixKernels[4] = dgemm_m2_n9_k4_ldA4_ldB4_ldC2_beta0_pfsigonly;
m_hardwareFlops[4] = 144;
m_nonZeroFlops[5] = 36;
m_matrixKernels[5] = dgemm_m2_n9_k4_ldA4_ldB4_ldC2_beta0_pfsigonly;
m_hardwareFlops[5] = 144;
m_nonZeroFlops[6] = 54;
m_matrixKernels[6] = dgemm_m2_n9_k4_ldA4_ldB4_ldC2_beta0_pfsigonly;
m_hardwareFlops[6] = 144;
m_nonZeroFlops[7] = 48;
m_hardwareFlops[7] = 48;
m_matrixKernels[7] = dsparse_starMatrix_m1_n9_k9_ldA2_ldBna3_ldC2_beta1_pfsigonly;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 126;
m_matrixKernels[0] = dgemm_m10_n9_k4_ldA10_ldB10_ldC10_beta0_pfsigonly;
m_hardwareFlops[0] = 720;
m_nonZeroFlops[1] = 306;
m_matrixKernels[1] = dgemm_m10_n9_k4_ldA10_ldB10_ldC10_beta0_pfsigonly;
m_hardwareFlops[1] = 720;
m_nonZeroFlops[2] = 396;
m_matrixKernels[2] = dgemm_m10_n9_k4_ldA10_ldB10_ldC10_beta0_pfsigonly;
m_hardwareFlops[2] = 720;
m_nonZeroFlops[3] = 480;
m_hardwareFlops[3] = 480;
m_matrixKernels[3] = dsparse_starMatrix_m10_n9_k9_ldA10_ldBna3_ldC10_beta1_pfsigonly;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 360;
m_hardwareFlops[0] = 1800;
m_matrixKernels[0] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
m_nonZeroFlops[1] = 612;
m_hardwareFlops[1] = 1800;
m_matrixKernels[1] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
m_nonZeroFlops[2] = 972;
m_hardwareFlops[2] = 1800;
m_matrixKernels[2] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
m_nonZeroFlops[3] = 972;
m_hardwareFlops[3] = 1800;
m_matrixKernels[3] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
m_nonZeroFlops[4] = 612;
m_hardwareFlops[4] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[4] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 612;
m_hardwareFlops[5] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[5] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 360;
m_hardwareFlops[6] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[6] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 720;
m_hardwareFlops[7] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[7] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 1224;
m_hardwareFlops[8] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[8] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 1224;
m_hardwareFlops[9] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[9] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 1458;
m_hardwareFlops[10] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[10] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 1368;
m_hardwareFlops[11] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[11] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 1368;
m_hardwareFlops[12] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[12] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 1368;
m_hardwareFlops[13] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[13] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 1368;
m_hardwareFlops[14] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[14] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 1458;
m_hardwareFlops[15] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[15] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 720;
m_hardwareFlops[16] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[16] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 1224;
m_hardwareFlops[17] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[17] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 1224;
m_hardwareFlops[18] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[18] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 1512;
m_hardwareFlops[19] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[19] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 612;
m_hardwareFlops[20] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[20] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 1512;
m_hardwareFlops[21] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[21] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 1098;
m_hardwareFlops[22] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[22] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 1512;
m_hardwareFlops[23] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[23] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 1512;
m_hardwareFlops[24] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[24] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 1512;
m_hardwareFlops[25] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[25] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 1098;
m_hardwareFlops[26] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[26] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 1512;
m_hardwareFlops[27] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[27] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 1458;
m_hardwareFlops[28] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[28] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 1368;
m_hardwareFlops[29] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[29] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 1368;
m_hardwareFlops[30] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[30] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 1098;
m_hardwareFlops[31] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[31] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 1512;
m_hardwareFlops[32] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[32] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 1512;
m_hardwareFlops[33] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[33] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 1224;
m_hardwareFlops[34] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[34] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 1764;
m_hardwareFlops[35] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[35] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 1260;
m_hardwareFlops[36] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[36] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 972;
m_hardwareFlops[37] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[37] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 1674;
m_hardwareFlops[38] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[38] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 1674;
m_hardwareFlops[39] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[39] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 1368;
m_hardwareFlops[40] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[40] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 1368;
m_hardwareFlops[41] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[41] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 1458;
m_hardwareFlops[42] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[42] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 1512;
m_hardwareFlops[43] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[43] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 1098;
m_hardwareFlops[44] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[44] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 1512;
m_hardwareFlops[45] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[45] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 972;
m_hardwareFlops[46] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[46] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 1674;
m_hardwareFlops[47] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[47] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 1674;
m_hardwareFlops[48] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[48] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 1224;
m_hardwareFlops[49] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[49] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 1260;
m_hardwareFlops[50] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[50] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 1764;
m_hardwareFlops[51] = 1800;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#else
m_matrixKernels[51] = dgemm_m10_n9_k10_ldA10_ldB10_ldC10_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 1620;
m_hardwareFlops[52] = 1620;
m_matrixKernels[52] = dgemm_m10_n9_k9_ldA10_ldB9_ldC10_beta1_pfsigonly;
m_nonZeroFlops[53] = 1620;
m_hardwareFlops[53] = 1620;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = dgemm_m10_n9_k9_ldA10_ldB9_ldC10_beta1_pfsigonly;
#else
m_matrixKernels[53] = dgemm_m10_n9_k9_ldA10_ldB9_ldC10_beta1_pfsigonly;
#endif
#endif

#ifdef SPARSE_SWITCH
m_sparseSwitch[0] = -1; 
m_sparseSwitch[1] = -1; 
m_sparseSwitch[2] = -1; 
m_sparseSwitch[3] = -1; 
m_sparseSwitch[4] = -1; 
m_sparseSwitch[5] = -1; 
m_sparseSwitch[6] = -1; 
m_sparseSwitch[7] = -1; 
m_sparseSwitch[8] = -1; 
m_sparseSwitch[9] = -1; 
m_sparseSwitch[10] = -1; 
m_sparseSwitch[11] = -1; 
m_sparseSwitch[12] = -1; 
m_sparseSwitch[13] = -1; 
m_sparseSwitch[14] = -1; 
m_sparseSwitch[15] = -1; 
m_sparseSwitch[16] = -1; 
m_sparseSwitch[17] = -1; 
m_sparseSwitch[18] = -1; 
m_sparseSwitch[19] = -1; 
m_sparseSwitch[20] = -1; 
m_sparseSwitch[21] = -1; 
m_sparseSwitch[22] = -1; 
m_sparseSwitch[23] = -1; 
m_sparseSwitch[24] = -1; 
m_sparseSwitch[25] = -1; 
m_sparseSwitch[26] = -1; 
m_sparseSwitch[27] = -1; 
m_sparseSwitch[28] = -1; 
m_sparseSwitch[29] = -1; 
m_sparseSwitch[30] = -1; 
m_sparseSwitch[31] = -1; 
m_sparseSwitch[32] = -1; 
m_sparseSwitch[33] = -1; 
m_sparseSwitch[34] = -1; 
m_sparseSwitch[35] = -1; 
m_sparseSwitch[36] = -1; 
m_sparseSwitch[37] = -1; 
m_sparseSwitch[38] = -1; 
m_sparseSwitch[39] = -1; 
m_sparseSwitch[40] = -1; 
m_sparseSwitch[41] = -1; 
m_sparseSwitch[42] = -1; 
m_sparseSwitch[43] = -1; 
m_sparseSwitch[44] = -1; 
m_sparseSwitch[45] = -1; 
m_sparseSwitch[46] = -1; 
m_sparseSwitch[47] = -1; 
m_sparseSwitch[48] = -1; 
m_sparseSwitch[49] = -1; 
m_sparseSwitch[50] = -1; 
m_sparseSwitch[51] = -1; 
m_sparseSwitch[52] = -1; 
m_sparseSwitch[53] = -1; 
m_sparseSwitch[54] = -1; 
m_sparseSwitch[55] = -1; 
m_sparseSwitch[56] = -1; 
m_sparseSwitch[57] = -1; 
m_sparseSwitch[58] = -1; 
m_sparseSwitch[59] = 24; 
#endif

#define STAR_NNZ 24

#endif

#if CONVERGENCE_ORDER==4

#ifdef TIME_KERNEL
m_nonZeroFlops[0] = 594;
m_matrixKernels[0] = dgemm_m10_n9_k20_ldA10_ldB20_ldC10_beta0_pfsigonly;
m_hardwareFlops[0] = 3600;
m_nonZeroFlops[1] = 1386;
m_matrixKernels[1] = dgemm_m10_n9_k20_ldA10_ldB20_ldC10_beta0_pfsigonly;
m_hardwareFlops[1] = 3600;
m_nonZeroFlops[2] = 1656;
m_matrixKernels[2] = dgemm_m10_n9_k20_ldA10_ldB20_ldC10_beta0_pfsigonly;
m_hardwareFlops[2] = 3600;
m_nonZeroFlops[3] = 480;
m_hardwareFlops[3] = 480;
m_matrixKernels[3] = dsparse_starMatrix_m10_n9_k9_ldA10_ldBna4_ldC10_beta1_pfsigonly;
m_nonZeroFlops[4] = 126;
m_matrixKernels[4] = dgemm_m4_n9_k10_ldA10_ldB10_ldC4_beta0_pfsigonly;
m_hardwareFlops[4] = 720;
m_nonZeroFlops[5] = 306;
m_matrixKernels[5] = dgemm_m4_n9_k10_ldA10_ldB10_ldC4_beta0_pfsigonly;
m_hardwareFlops[5] = 720;
m_nonZeroFlops[6] = 396;
m_matrixKernels[6] = dgemm_m4_n9_k10_ldA10_ldB10_ldC4_beta0_pfsigonly;
m_hardwareFlops[6] = 720;
m_nonZeroFlops[7] = 192;
m_hardwareFlops[7] = 192;
m_matrixKernels[7] = dsparse_starMatrix_m4_n9_k9_ldA4_ldBna4_ldC4_beta1_pfsigonly;
m_nonZeroFlops[8] = 18;
m_matrixKernels[8] = dgemm_m2_n9_k4_ldA10_ldB4_ldC2_beta0_pfsigonly;
m_hardwareFlops[8] = 144;
m_nonZeroFlops[9] = 36;
m_matrixKernels[9] = dgemm_m2_n9_k4_ldA10_ldB4_ldC2_beta0_pfsigonly;
m_hardwareFlops[9] = 144;
m_nonZeroFlops[10] = 54;
m_matrixKernels[10] = dgemm_m2_n9_k4_ldA10_ldB4_ldC2_beta0_pfsigonly;
m_hardwareFlops[10] = 144;
m_nonZeroFlops[11] = 48;
m_hardwareFlops[11] = 48;
m_matrixKernels[11] = dsparse_starMatrix_m1_n9_k9_ldA2_ldBna4_ldC2_beta1_pfsigonly;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 594;
m_matrixKernels[0] = dgemm_m20_n9_k10_ldA20_ldB20_ldC20_beta0_pfsigonly;
m_hardwareFlops[0] = 3600;
m_nonZeroFlops[1] = 1386;
m_matrixKernels[1] = dgemm_m20_n9_k10_ldA20_ldB20_ldC20_beta0_pfsigonly;
m_hardwareFlops[1] = 3600;
m_nonZeroFlops[2] = 1656;
m_matrixKernels[2] = dgemm_m20_n9_k10_ldA20_ldB20_ldC20_beta0_pfsigonly;
m_hardwareFlops[2] = 3600;
m_nonZeroFlops[3] = 960;
m_hardwareFlops[3] = 960;
m_matrixKernels[3] = dsparse_starMatrix_m20_n9_k9_ldA20_ldBna4_ldC20_beta1_pfsigonly;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 900;
m_hardwareFlops[0] = 7200;
m_matrixKernels[0] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
m_nonZeroFlops[1] = 1872;
m_hardwareFlops[1] = 7200;
m_matrixKernels[1] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
m_nonZeroFlops[2] = 3672;
m_hardwareFlops[2] = 7200;
m_matrixKernels[2] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
m_nonZeroFlops[3] = 3672;
m_hardwareFlops[3] = 7200;
m_matrixKernels[3] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
m_nonZeroFlops[4] = 1872;
m_hardwareFlops[4] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[4] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 1872;
m_hardwareFlops[5] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[5] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 900;
m_hardwareFlops[6] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[6] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 2250;
m_hardwareFlops[7] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[7] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 4680;
m_hardwareFlops[8] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[8] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 4680;
m_hardwareFlops[9] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[9] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 5760;
m_hardwareFlops[10] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[10] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 5148;
m_hardwareFlops[11] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[11] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 5310;
m_hardwareFlops[12] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[12] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 5310;
m_hardwareFlops[13] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[13] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 5148;
m_hardwareFlops[14] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[14] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 5760;
m_hardwareFlops[15] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[15] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 2250;
m_hardwareFlops[16] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[16] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 4680;
m_hardwareFlops[17] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[17] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 4680;
m_hardwareFlops[18] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[18] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 6084;
m_hardwareFlops[19] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[19] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 1872;
m_hardwareFlops[20] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[20] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 6084;
m_hardwareFlops[21] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[21] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 4176;
m_hardwareFlops[22] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[22] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 6120;
m_hardwareFlops[23] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[23] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 5850;
m_hardwareFlops[24] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[24] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 6120;
m_hardwareFlops[25] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[25] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 4176;
m_hardwareFlops[26] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[26] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 5850;
m_hardwareFlops[27] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[27] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 5760;
m_hardwareFlops[28] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[28] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 5148;
m_hardwareFlops[29] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[29] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 5310;
m_hardwareFlops[30] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[30] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 4176;
m_hardwareFlops[31] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[31] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 6120;
m_hardwareFlops[32] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[32] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 5850;
m_hardwareFlops[33] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[33] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 4644;
m_hardwareFlops[34] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[34] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 7092;
m_hardwareFlops[35] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[35] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 4932;
m_hardwareFlops[36] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[36] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 3672;
m_hardwareFlops[37] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[37] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 6678;
m_hardwareFlops[38] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[38] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 6678;
m_hardwareFlops[39] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[39] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 5310;
m_hardwareFlops[40] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[40] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 5148;
m_hardwareFlops[41] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[41] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 5760;
m_hardwareFlops[42] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[42] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 6120;
m_hardwareFlops[43] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[43] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 4176;
m_hardwareFlops[44] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[44] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 5850;
m_hardwareFlops[45] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[45] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 3672;
m_hardwareFlops[46] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[46] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 6678;
m_hardwareFlops[47] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[47] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 6678;
m_hardwareFlops[48] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[48] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 4644;
m_hardwareFlops[49] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[49] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 4932;
m_hardwareFlops[50] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[50] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 7092;
m_hardwareFlops[51] = 7200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#else
m_matrixKernels[51] = dgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 3240;
m_hardwareFlops[52] = 3240;
m_matrixKernels[52] = dgemm_m20_n9_k9_ldA20_ldB9_ldC20_beta1_pfsigonly;
m_nonZeroFlops[53] = 3240;
m_hardwareFlops[53] = 3240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = dgemm_m20_n9_k9_ldA20_ldB9_ldC20_beta1_pfsigonly;
#else
m_matrixKernels[53] = dgemm_m20_n9_k9_ldA20_ldB9_ldC20_beta1_pfsigonly;
#endif
#endif

#ifdef SPARSE_SWITCH
m_sparseSwitch[0] = -1; 
m_sparseSwitch[1] = -1; 
m_sparseSwitch[2] = -1; 
m_sparseSwitch[3] = -1; 
m_sparseSwitch[4] = -1; 
m_sparseSwitch[5] = -1; 
m_sparseSwitch[6] = -1; 
m_sparseSwitch[7] = -1; 
m_sparseSwitch[8] = -1; 
m_sparseSwitch[9] = -1; 
m_sparseSwitch[10] = -1; 
m_sparseSwitch[11] = -1; 
m_sparseSwitch[12] = -1; 
m_sparseSwitch[13] = -1; 
m_sparseSwitch[14] = -1; 
m_sparseSwitch[15] = -1; 
m_sparseSwitch[16] = -1; 
m_sparseSwitch[17] = -1; 
m_sparseSwitch[18] = -1; 
m_sparseSwitch[19] = -1; 
m_sparseSwitch[20] = -1; 
m_sparseSwitch[21] = -1; 
m_sparseSwitch[22] = -1; 
m_sparseSwitch[23] = -1; 
m_sparseSwitch[24] = -1; 
m_sparseSwitch[25] = -1; 
m_sparseSwitch[26] = -1; 
m_sparseSwitch[27] = -1; 
m_sparseSwitch[28] = -1; 
m_sparseSwitch[29] = -1; 
m_sparseSwitch[30] = -1; 
m_sparseSwitch[31] = -1; 
m_sparseSwitch[32] = -1; 
m_sparseSwitch[33] = -1; 
m_sparseSwitch[34] = -1; 
m_sparseSwitch[35] = -1; 
m_sparseSwitch[36] = -1; 
m_sparseSwitch[37] = -1; 
m_sparseSwitch[38] = -1; 
m_sparseSwitch[39] = -1; 
m_sparseSwitch[40] = -1; 
m_sparseSwitch[41] = -1; 
m_sparseSwitch[42] = -1; 
m_sparseSwitch[43] = -1; 
m_sparseSwitch[44] = -1; 
m_sparseSwitch[45] = -1; 
m_sparseSwitch[46] = -1; 
m_sparseSwitch[47] = -1; 
m_sparseSwitch[48] = -1; 
m_sparseSwitch[49] = -1; 
m_sparseSwitch[50] = -1; 
m_sparseSwitch[51] = -1; 
m_sparseSwitch[52] = -1; 
m_sparseSwitch[53] = -1; 
m_sparseSwitch[54] = -1; 
m_sparseSwitch[55] = -1; 
m_sparseSwitch[56] = -1; 
m_sparseSwitch[57] = -1; 
m_sparseSwitch[58] = -1; 
m_sparseSwitch[59] = 24; 
#endif

#define STAR_NNZ 24

#endif

#if CONVERGENCE_ORDER==5

#ifdef TIME_KERNEL
m_nonZeroFlops[0] = 1944;
m_matrixKernels[0] = dgemm_m20_n9_k35_ldA20_ldB36_ldC20_beta0_pfsigonly;
m_hardwareFlops[0] = 12600;
m_nonZeroFlops[1] = 4536;
m_matrixKernels[1] = dgemm_m20_n9_k35_ldA20_ldB36_ldC20_beta0_pfsigonly;
m_hardwareFlops[1] = 12600;
m_nonZeroFlops[2] = 5166;
m_matrixKernels[2] = dgemm_m20_n9_k35_ldA20_ldB36_ldC20_beta0_pfsigonly;
m_hardwareFlops[2] = 12600;
m_nonZeroFlops[3] = 960;
m_hardwareFlops[3] = 960;
m_matrixKernels[3] = dsparse_starMatrix_m20_n9_k9_ldA20_ldBna5_ldC20_beta1_pfsigonly;
m_nonZeroFlops[4] = 594;
m_matrixKernels[4] = dgemm_m10_n9_k20_ldA20_ldB20_ldC10_beta0_pfsigonly;
m_hardwareFlops[4] = 3600;
m_nonZeroFlops[5] = 1386;
m_matrixKernels[5] = dgemm_m10_n9_k20_ldA20_ldB20_ldC10_beta0_pfsigonly;
m_hardwareFlops[5] = 3600;
m_nonZeroFlops[6] = 1656;
m_matrixKernels[6] = dgemm_m10_n9_k20_ldA20_ldB20_ldC10_beta0_pfsigonly;
m_hardwareFlops[6] = 3600;
m_nonZeroFlops[7] = 480;
m_hardwareFlops[7] = 480;
m_matrixKernels[7] = dsparse_starMatrix_m10_n9_k9_ldA10_ldBna5_ldC10_beta1_pfsigonly;
m_nonZeroFlops[8] = 126;
m_matrixKernels[8] = dgemm_m4_n9_k10_ldA20_ldB10_ldC4_beta0_pfsigonly;
m_hardwareFlops[8] = 720;
m_nonZeroFlops[9] = 306;
m_matrixKernels[9] = dgemm_m4_n9_k10_ldA20_ldB10_ldC4_beta0_pfsigonly;
m_hardwareFlops[9] = 720;
m_nonZeroFlops[10] = 396;
m_matrixKernels[10] = dgemm_m4_n9_k10_ldA20_ldB10_ldC4_beta0_pfsigonly;
m_hardwareFlops[10] = 720;
m_nonZeroFlops[11] = 192;
m_hardwareFlops[11] = 192;
m_matrixKernels[11] = dsparse_starMatrix_m4_n9_k9_ldA4_ldBna5_ldC4_beta1_pfsigonly;
m_nonZeroFlops[12] = 18;
m_matrixKernels[12] = dgemm_m2_n9_k4_ldA20_ldB4_ldC2_beta0_pfsigonly;
m_hardwareFlops[12] = 144;
m_nonZeroFlops[13] = 36;
m_matrixKernels[13] = dgemm_m2_n9_k4_ldA20_ldB4_ldC2_beta0_pfsigonly;
m_hardwareFlops[13] = 144;
m_nonZeroFlops[14] = 54;
m_matrixKernels[14] = dgemm_m2_n9_k4_ldA20_ldB4_ldC2_beta0_pfsigonly;
m_hardwareFlops[14] = 144;
m_nonZeroFlops[15] = 48;
m_hardwareFlops[15] = 48;
m_matrixKernels[15] = dsparse_starMatrix_m1_n9_k9_ldA2_ldBna5_ldC2_beta1_pfsigonly;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 1944;
m_matrixKernels[0] = dgemm_m36_n9_k20_ldA36_ldB36_ldC36_beta0_pfsigonly;
m_hardwareFlops[0] = 12960;
m_nonZeroFlops[1] = 4536;
m_matrixKernels[1] = dgemm_m36_n9_k20_ldA36_ldB36_ldC36_beta0_pfsigonly;
m_hardwareFlops[1] = 12960;
m_nonZeroFlops[2] = 5166;
m_matrixKernels[2] = dgemm_m36_n9_k20_ldA36_ldB36_ldC36_beta0_pfsigonly;
m_hardwareFlops[2] = 12960;
m_nonZeroFlops[3] = 1680;
m_hardwareFlops[3] = 1680;
m_matrixKernels[3] = dsparse_starMatrix_m35_n9_k9_ldA36_ldBna5_ldC36_beta1_pfsigonly;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 1890;
m_hardwareFlops[0] = 22680;
m_matrixKernels[0] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
m_nonZeroFlops[1] = 4662;
m_hardwareFlops[1] = 22680;
m_matrixKernels[1] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
m_nonZeroFlops[2] = 10962;
m_hardwareFlops[2] = 22680;
m_matrixKernels[2] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
m_nonZeroFlops[3] = 10962;
m_hardwareFlops[3] = 22680;
m_matrixKernels[3] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
m_nonZeroFlops[4] = 4626;
m_hardwareFlops[4] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[4] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 4626;
m_hardwareFlops[5] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[5] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 1890;
m_hardwareFlops[6] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[6] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 5670;
m_hardwareFlops[7] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[7] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 13878;
m_hardwareFlops[8] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[8] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 13878;
m_hardwareFlops[9] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[9] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 17622;
m_hardwareFlops[10] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[10] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 15480;
m_hardwareFlops[11] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[11] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 15768;
m_hardwareFlops[12] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[12] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 15768;
m_hardwareFlops[13] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[13] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 15480;
m_hardwareFlops[14] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[14] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 17622;
m_hardwareFlops[15] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[15] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 5670;
m_hardwareFlops[16] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[16] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 13878;
m_hardwareFlops[17] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[17] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 13878;
m_hardwareFlops[18] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[18] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 18810;
m_hardwareFlops[19] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[19] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 4662;
m_hardwareFlops[20] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[20] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 18810;
m_hardwareFlops[21] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[21] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 12312;
m_hardwareFlops[22] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[22] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 18702;
m_hardwareFlops[23] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[23] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 17640;
m_hardwareFlops[24] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[24] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 18702;
m_hardwareFlops[25] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[25] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 12312;
m_hardwareFlops[26] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[26] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 17640;
m_hardwareFlops[27] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[27] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 17622;
m_hardwareFlops[28] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[28] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 15480;
m_hardwareFlops[29] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[29] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 15768;
m_hardwareFlops[30] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[30] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 12312;
m_hardwareFlops[31] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[31] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 18702;
m_hardwareFlops[32] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[32] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 17640;
m_hardwareFlops[33] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[33] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 13914;
m_hardwareFlops[34] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[34] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 21834;
m_hardwareFlops[35] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[35] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 14958;
m_hardwareFlops[36] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[36] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 10962;
m_hardwareFlops[37] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[37] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 20358;
m_hardwareFlops[38] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[38] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 20358;
m_hardwareFlops[39] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[39] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 15768;
m_hardwareFlops[40] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[40] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 15480;
m_hardwareFlops[41] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[41] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 17622;
m_hardwareFlops[42] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[42] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 18702;
m_hardwareFlops[43] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[43] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 12312;
m_hardwareFlops[44] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[44] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 17640;
m_hardwareFlops[45] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[45] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 10962;
m_hardwareFlops[46] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[46] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 20358;
m_hardwareFlops[47] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[47] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 20358;
m_hardwareFlops[48] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[48] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 13914;
m_hardwareFlops[49] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[49] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 14958;
m_hardwareFlops[50] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[50] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 21834;
m_hardwareFlops[51] = 22680;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#else
m_matrixKernels[51] = dgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 5670;
m_hardwareFlops[52] = 5832;
m_matrixKernels[52] = dgemm_m36_n9_k9_ldA36_ldB9_ldC36_beta1_pfsigonly;
m_nonZeroFlops[53] = 5670;
m_hardwareFlops[53] = 5832;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = dgemm_m36_n9_k9_ldA36_ldB9_ldC36_beta1_pfsigonly;
#else
m_matrixKernels[53] = dgemm_m36_n9_k9_ldA36_ldB9_ldC36_beta1_pfsigonly;
#endif
#endif

#ifdef SPARSE_SWITCH
m_sparseSwitch[0] = -1; 
m_sparseSwitch[1] = -1; 
m_sparseSwitch[2] = -1; 
m_sparseSwitch[3] = -1; 
m_sparseSwitch[4] = -1; 
m_sparseSwitch[5] = -1; 
m_sparseSwitch[6] = -1; 
m_sparseSwitch[7] = -1; 
m_sparseSwitch[8] = -1; 
m_sparseSwitch[9] = -1; 
m_sparseSwitch[10] = -1; 
m_sparseSwitch[11] = -1; 
m_sparseSwitch[12] = -1; 
m_sparseSwitch[13] = -1; 
m_sparseSwitch[14] = -1; 
m_sparseSwitch[15] = -1; 
m_sparseSwitch[16] = -1; 
m_sparseSwitch[17] = -1; 
m_sparseSwitch[18] = -1; 
m_sparseSwitch[19] = -1; 
m_sparseSwitch[20] = -1; 
m_sparseSwitch[21] = -1; 
m_sparseSwitch[22] = -1; 
m_sparseSwitch[23] = -1; 
m_sparseSwitch[24] = -1; 
m_sparseSwitch[25] = -1; 
m_sparseSwitch[26] = -1; 
m_sparseSwitch[27] = -1; 
m_sparseSwitch[28] = -1; 
m_sparseSwitch[29] = -1; 
m_sparseSwitch[30] = -1; 
m_sparseSwitch[31] = -1; 
m_sparseSwitch[32] = -1; 
m_sparseSwitch[33] = -1; 
m_sparseSwitch[34] = -1; 
m_sparseSwitch[35] = -1; 
m_sparseSwitch[36] = -1; 
m_sparseSwitch[37] = -1; 
m_sparseSwitch[38] = -1; 
m_sparseSwitch[39] = -1; 
m_sparseSwitch[40] = -1; 
m_sparseSwitch[41] = -1; 
m_sparseSwitch[42] = -1; 
m_sparseSwitch[43] = -1; 
m_sparseSwitch[44] = -1; 
m_sparseSwitch[45] = -1; 
m_sparseSwitch[46] = -1; 
m_sparseSwitch[47] = -1; 
m_sparseSwitch[48] = -1; 
m_sparseSwitch[49] = -1; 
m_sparseSwitch[50] = -1; 
m_sparseSwitch[51] = -1; 
m_sparseSwitch[52] = -1; 
m_sparseSwitch[53] = -1; 
m_sparseSwitch[54] = -1; 
m_sparseSwitch[55] = -1; 
m_sparseSwitch[56] = -1; 
m_sparseSwitch[57] = -1; 
m_sparseSwitch[58] = -1; 
m_sparseSwitch[59] = 24; 
#endif

#define STAR_NNZ 24

#endif

#if CONVERGENCE_ORDER==6

#ifdef TIME_KERNEL
m_nonZeroFlops[0] = 5292;
m_matrixKernels[0] = dgemm_m36_n9_k56_ldA36_ldB56_ldC36_beta0_pfsigonly;
m_hardwareFlops[0] = 36288;
m_nonZeroFlops[1] = 12096;
m_matrixKernels[1] = dgemm_m36_n9_k56_ldA36_ldB56_ldC36_beta0_pfsigonly;
m_hardwareFlops[1] = 36288;
m_nonZeroFlops[2] = 13356;
m_matrixKernels[2] = dgemm_m36_n9_k56_ldA36_ldB56_ldC36_beta0_pfsigonly;
m_hardwareFlops[2] = 36288;
m_nonZeroFlops[3] = 1680;
m_hardwareFlops[3] = 1680;
m_matrixKernels[3] = dsparse_starMatrix_m35_n9_k9_ldA36_ldBna6_ldC36_beta1_pfsigonly;
m_nonZeroFlops[4] = 1944;
m_matrixKernels[4] = dgemm_m20_n9_k35_ldA36_ldB36_ldC20_beta0_pfsigonly;
m_hardwareFlops[4] = 12600;
m_nonZeroFlops[5] = 4536;
m_matrixKernels[5] = dgemm_m20_n9_k35_ldA36_ldB36_ldC20_beta0_pfsigonly;
m_hardwareFlops[5] = 12600;
m_nonZeroFlops[6] = 5166;
m_matrixKernels[6] = dgemm_m20_n9_k35_ldA36_ldB36_ldC20_beta0_pfsigonly;
m_hardwareFlops[6] = 12600;
m_nonZeroFlops[7] = 960;
m_hardwareFlops[7] = 960;
m_matrixKernels[7] = dsparse_starMatrix_m20_n9_k9_ldA20_ldBna6_ldC20_beta1_pfsigonly;
m_nonZeroFlops[8] = 594;
m_matrixKernels[8] = dgemm_m10_n9_k20_ldA36_ldB20_ldC10_beta0_pfsigonly;
m_hardwareFlops[8] = 3600;
m_nonZeroFlops[9] = 1386;
m_matrixKernels[9] = dgemm_m10_n9_k20_ldA36_ldB20_ldC10_beta0_pfsigonly;
m_hardwareFlops[9] = 3600;
m_nonZeroFlops[10] = 1656;
m_matrixKernels[10] = dgemm_m10_n9_k20_ldA36_ldB20_ldC10_beta0_pfsigonly;
m_hardwareFlops[10] = 3600;
m_nonZeroFlops[11] = 480;
m_hardwareFlops[11] = 480;
m_matrixKernels[11] = dsparse_starMatrix_m10_n9_k9_ldA10_ldBna6_ldC10_beta1_pfsigonly;
m_nonZeroFlops[12] = 126;
m_matrixKernels[12] = dgemm_m4_n9_k10_ldA36_ldB10_ldC4_beta0_pfsigonly;
m_hardwareFlops[12] = 720;
m_nonZeroFlops[13] = 306;
m_matrixKernels[13] = dgemm_m4_n9_k10_ldA36_ldB10_ldC4_beta0_pfsigonly;
m_hardwareFlops[13] = 720;
m_nonZeroFlops[14] = 396;
m_matrixKernels[14] = dgemm_m4_n9_k10_ldA36_ldB10_ldC4_beta0_pfsigonly;
m_hardwareFlops[14] = 720;
m_nonZeroFlops[15] = 192;
m_hardwareFlops[15] = 192;
m_matrixKernels[15] = dsparse_starMatrix_m4_n9_k9_ldA4_ldBna6_ldC4_beta1_pfsigonly;
m_nonZeroFlops[16] = 18;
m_matrixKernels[16] = dgemm_m2_n9_k4_ldA36_ldB4_ldC2_beta0_pfsigonly;
m_hardwareFlops[16] = 144;
m_nonZeroFlops[17] = 36;
m_matrixKernels[17] = dgemm_m2_n9_k4_ldA36_ldB4_ldC2_beta0_pfsigonly;
m_hardwareFlops[17] = 144;
m_nonZeroFlops[18] = 54;
m_matrixKernels[18] = dgemm_m2_n9_k4_ldA36_ldB4_ldC2_beta0_pfsigonly;
m_hardwareFlops[18] = 144;
m_nonZeroFlops[19] = 48;
m_hardwareFlops[19] = 48;
m_matrixKernels[19] = dsparse_starMatrix_m1_n9_k9_ldA2_ldBna6_ldC2_beta1_pfsigonly;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 5292;
m_matrixKernels[0] = dgemm_m56_n9_k35_ldA56_ldB56_ldC56_beta0_pfsigonly;
m_hardwareFlops[0] = 35280;
m_nonZeroFlops[1] = 12096;
m_matrixKernels[1] = dgemm_m56_n9_k35_ldA56_ldB56_ldC56_beta0_pfsigonly;
m_hardwareFlops[1] = 35280;
m_nonZeroFlops[2] = 13356;
m_matrixKernels[2] = dgemm_m56_n9_k35_ldA56_ldB56_ldC56_beta0_pfsigonly;
m_hardwareFlops[2] = 35280;
m_nonZeroFlops[3] = 2688;
m_hardwareFlops[3] = 2688;
m_matrixKernels[3] = dsparse_starMatrix_m56_n9_k9_ldA56_ldBna6_ldC56_beta1_pfsigonly;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 3528;
m_hardwareFlops[0] = 56448;
m_matrixKernels[0] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
m_nonZeroFlops[1] = 10080;
m_hardwareFlops[1] = 56448;
m_matrixKernels[1] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
m_nonZeroFlops[2] = 27720;
m_hardwareFlops[2] = 56448;
m_matrixKernels[2] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
m_nonZeroFlops[3] = 27720;
m_hardwareFlops[3] = 56448;
m_matrixKernels[3] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
m_nonZeroFlops[4] = 9936;
m_hardwareFlops[4] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[4] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 9936;
m_hardwareFlops[5] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[5] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 3528;
m_hardwareFlops[6] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[6] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 12348;
m_hardwareFlops[7] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[7] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 34776;
m_hardwareFlops[8] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[8] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 34776;
m_hardwareFlops[9] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[9] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 45324;
m_hardwareFlops[10] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[10] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 39852;
m_hardwareFlops[11] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[11] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 39852;
m_hardwareFlops[12] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[12] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 39852;
m_hardwareFlops[13] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[13] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 39852;
m_hardwareFlops[14] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[14] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 45324;
m_hardwareFlops[15] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[15] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 12348;
m_hardwareFlops[16] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[16] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 34776;
m_hardwareFlops[17] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[17] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 34776;
m_hardwareFlops[18] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[18] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 48564;
m_hardwareFlops[19] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[19] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 10080;
m_hardwareFlops[20] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[20] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 48564;
m_hardwareFlops[21] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[21] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 30996;
m_hardwareFlops[22] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[22] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 48348;
m_hardwareFlops[23] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[23] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 45162;
m_hardwareFlops[24] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[24] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 48348;
m_hardwareFlops[25] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[25] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 30996;
m_hardwareFlops[26] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[26] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 45162;
m_hardwareFlops[27] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[27] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 45324;
m_hardwareFlops[28] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[28] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 39852;
m_hardwareFlops[29] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[29] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 39852;
m_hardwareFlops[30] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[30] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 30996;
m_hardwareFlops[31] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[31] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 48348;
m_hardwareFlops[32] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[32] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 45162;
m_hardwareFlops[33] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[33] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 35028;
m_hardwareFlops[34] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[34] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 56232;
m_hardwareFlops[35] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[35] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 38160;
m_hardwareFlops[36] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[36] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 27720;
m_hardwareFlops[37] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[37] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 52290;
m_hardwareFlops[38] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[38] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 52290;
m_hardwareFlops[39] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[39] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 39852;
m_hardwareFlops[40] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[40] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 39852;
m_hardwareFlops[41] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[41] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 45324;
m_hardwareFlops[42] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[42] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 48348;
m_hardwareFlops[43] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[43] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 30996;
m_hardwareFlops[44] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[44] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 45162;
m_hardwareFlops[45] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[45] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 27720;
m_hardwareFlops[46] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[46] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 52290;
m_hardwareFlops[47] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[47] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 52290;
m_hardwareFlops[48] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[48] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 35028;
m_hardwareFlops[49] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[49] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 38160;
m_hardwareFlops[50] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[50] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 56232;
m_hardwareFlops[51] = 56448;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#else
m_matrixKernels[51] = dgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 9072;
m_hardwareFlops[52] = 9072;
m_matrixKernels[52] = dgemm_m56_n9_k9_ldA56_ldB9_ldC56_beta1_pfsigonly;
m_nonZeroFlops[53] = 9072;
m_hardwareFlops[53] = 9072;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = dgemm_m56_n9_k9_ldA56_ldB9_ldC56_beta1_pfsigonly;
#else
m_matrixKernels[53] = dgemm_m56_n9_k9_ldA56_ldB9_ldC56_beta1_pfsigonly;
#endif
#endif

#ifdef SPARSE_SWITCH
m_sparseSwitch[0] = -1; 
m_sparseSwitch[1] = -1; 
m_sparseSwitch[2] = -1; 
m_sparseSwitch[3] = -1; 
m_sparseSwitch[4] = -1; 
m_sparseSwitch[5] = -1; 
m_sparseSwitch[6] = -1; 
m_sparseSwitch[7] = -1; 
m_sparseSwitch[8] = -1; 
m_sparseSwitch[9] = -1; 
m_sparseSwitch[10] = -1; 
m_sparseSwitch[11] = -1; 
m_sparseSwitch[12] = -1; 
m_sparseSwitch[13] = -1; 
m_sparseSwitch[14] = -1; 
m_sparseSwitch[15] = -1; 
m_sparseSwitch[16] = -1; 
m_sparseSwitch[17] = -1; 
m_sparseSwitch[18] = -1; 
m_sparseSwitch[19] = -1; 
m_sparseSwitch[20] = -1; 
m_sparseSwitch[21] = -1; 
m_sparseSwitch[22] = -1; 
m_sparseSwitch[23] = -1; 
m_sparseSwitch[24] = -1; 
m_sparseSwitch[25] = -1; 
m_sparseSwitch[26] = -1; 
m_sparseSwitch[27] = -1; 
m_sparseSwitch[28] = -1; 
m_sparseSwitch[29] = -1; 
m_sparseSwitch[30] = -1; 
m_sparseSwitch[31] = -1; 
m_sparseSwitch[32] = -1; 
m_sparseSwitch[33] = -1; 
m_sparseSwitch[34] = -1; 
m_sparseSwitch[35] = -1; 
m_sparseSwitch[36] = -1; 
m_sparseSwitch[37] = -1; 
m_sparseSwitch[38] = -1; 
m_sparseSwitch[39] = -1; 
m_sparseSwitch[40] = -1; 
m_sparseSwitch[41] = -1; 
m_sparseSwitch[42] = -1; 
m_sparseSwitch[43] = -1; 
m_sparseSwitch[44] = -1; 
m_sparseSwitch[45] = -1; 
m_sparseSwitch[46] = -1; 
m_sparseSwitch[47] = -1; 
m_sparseSwitch[48] = -1; 
m_sparseSwitch[49] = -1; 
m_sparseSwitch[50] = -1; 
m_sparseSwitch[51] = -1; 
m_sparseSwitch[52] = -1; 
m_sparseSwitch[53] = -1; 
m_sparseSwitch[54] = -1; 
m_sparseSwitch[55] = -1; 
m_sparseSwitch[56] = -1; 
m_sparseSwitch[57] = -1; 
m_sparseSwitch[58] = -1; 
m_sparseSwitch[59] = 24; 
#endif

#define STAR_NNZ 24

#endif

#if CONVERGENCE_ORDER==7

#ifdef TIME_KERNEL
m_nonZeroFlops[0] = 12348;
m_matrixKernels[0] = dgemm_m56_n9_k84_ldA56_ldB84_ldC56_beta0_pfsigonly;
m_hardwareFlops[0] = 84672;
m_nonZeroFlops[1] = 27972;
m_matrixKernels[1] = dgemm_m56_n9_k84_ldA56_ldB84_ldC56_beta0_pfsigonly;
m_hardwareFlops[1] = 84672;
m_nonZeroFlops[2] = 30240;
m_matrixKernels[2] = dgemm_m56_n9_k84_ldA56_ldB84_ldC56_beta0_pfsigonly;
m_hardwareFlops[2] = 84672;
m_nonZeroFlops[3] = 2688;
m_hardwareFlops[3] = 2688;
m_matrixKernels[3] = dsparse_starMatrix_m56_n9_k9_ldA56_ldBna7_ldC56_beta1_pfsigonly;
m_nonZeroFlops[4] = 5292;
m_matrixKernels[4] = dgemm_m36_n9_k56_ldA56_ldB56_ldC36_beta0_pfsigonly;
m_hardwareFlops[4] = 36288;
m_nonZeroFlops[5] = 12096;
m_matrixKernels[5] = dgemm_m36_n9_k56_ldA56_ldB56_ldC36_beta0_pfsigonly;
m_hardwareFlops[5] = 36288;
m_nonZeroFlops[6] = 13356;
m_matrixKernels[6] = dgemm_m36_n9_k56_ldA56_ldB56_ldC36_beta0_pfsigonly;
m_hardwareFlops[6] = 36288;
m_nonZeroFlops[7] = 1680;
m_hardwareFlops[7] = 1680;
m_matrixKernels[7] = dsparse_starMatrix_m35_n9_k9_ldA36_ldBna7_ldC36_beta1_pfsigonly;
m_nonZeroFlops[8] = 1944;
m_matrixKernels[8] = dgemm_m20_n9_k35_ldA56_ldB36_ldC20_beta0_pfsigonly;
m_hardwareFlops[8] = 12600;
m_nonZeroFlops[9] = 4536;
m_matrixKernels[9] = dgemm_m20_n9_k35_ldA56_ldB36_ldC20_beta0_pfsigonly;
m_hardwareFlops[9] = 12600;
m_nonZeroFlops[10] = 5166;
m_matrixKernels[10] = dgemm_m20_n9_k35_ldA56_ldB36_ldC20_beta0_pfsigonly;
m_hardwareFlops[10] = 12600;
m_nonZeroFlops[11] = 960;
m_hardwareFlops[11] = 960;
m_matrixKernels[11] = dsparse_starMatrix_m20_n9_k9_ldA20_ldBna7_ldC20_beta1_pfsigonly;
m_nonZeroFlops[12] = 594;
m_matrixKernels[12] = dgemm_m10_n9_k20_ldA56_ldB20_ldC10_beta0_pfsigonly;
m_hardwareFlops[12] = 3600;
m_nonZeroFlops[13] = 1386;
m_matrixKernels[13] = dgemm_m10_n9_k20_ldA56_ldB20_ldC10_beta0_pfsigonly;
m_hardwareFlops[13] = 3600;
m_nonZeroFlops[14] = 1656;
m_matrixKernels[14] = dgemm_m10_n9_k20_ldA56_ldB20_ldC10_beta0_pfsigonly;
m_hardwareFlops[14] = 3600;
m_nonZeroFlops[15] = 480;
m_hardwareFlops[15] = 480;
m_matrixKernels[15] = dsparse_starMatrix_m10_n9_k9_ldA10_ldBna7_ldC10_beta1_pfsigonly;
m_nonZeroFlops[16] = 126;
m_matrixKernels[16] = dgemm_m4_n9_k10_ldA56_ldB10_ldC4_beta0_pfsigonly;
m_hardwareFlops[16] = 720;
m_nonZeroFlops[17] = 306;
m_matrixKernels[17] = dgemm_m4_n9_k10_ldA56_ldB10_ldC4_beta0_pfsigonly;
m_hardwareFlops[17] = 720;
m_nonZeroFlops[18] = 396;
m_matrixKernels[18] = dgemm_m4_n9_k10_ldA56_ldB10_ldC4_beta0_pfsigonly;
m_hardwareFlops[18] = 720;
m_nonZeroFlops[19] = 192;
m_hardwareFlops[19] = 192;
m_matrixKernels[19] = dsparse_starMatrix_m4_n9_k9_ldA4_ldBna7_ldC4_beta1_pfsigonly;
m_nonZeroFlops[20] = 18;
m_matrixKernels[20] = dgemm_m2_n9_k4_ldA56_ldB4_ldC2_beta0_pfsigonly;
m_hardwareFlops[20] = 144;
m_nonZeroFlops[21] = 36;
m_matrixKernels[21] = dgemm_m2_n9_k4_ldA56_ldB4_ldC2_beta0_pfsigonly;
m_hardwareFlops[21] = 144;
m_nonZeroFlops[22] = 54;
m_matrixKernels[22] = dgemm_m2_n9_k4_ldA56_ldB4_ldC2_beta0_pfsigonly;
m_hardwareFlops[22] = 144;
m_nonZeroFlops[23] = 48;
m_hardwareFlops[23] = 48;
m_matrixKernels[23] = dsparse_starMatrix_m1_n9_k9_ldA2_ldBna7_ldC2_beta1_pfsigonly;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 12348;
m_matrixKernels[0] = dgemm_m84_n9_k56_ldA84_ldB84_ldC84_beta0_pfsigonly;
m_hardwareFlops[0] = 84672;
m_nonZeroFlops[1] = 27972;
m_matrixKernels[1] = dgemm_m84_n9_k56_ldA84_ldB84_ldC84_beta0_pfsigonly;
m_hardwareFlops[1] = 84672;
m_nonZeroFlops[2] = 30240;
m_matrixKernels[2] = dgemm_m84_n9_k56_ldA84_ldB84_ldC84_beta0_pfsigonly;
m_hardwareFlops[2] = 84672;
m_nonZeroFlops[3] = 4032;
m_hardwareFlops[3] = 4032;
m_matrixKernels[3] = dsparse_starMatrix_m84_n9_k9_ldA84_ldBna7_ldC84_beta1_pfsigonly;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 6048;
m_hardwareFlops[0] = 127008;
m_matrixKernels[0] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
m_nonZeroFlops[1] = 19656;
m_hardwareFlops[1] = 127008;
m_matrixKernels[1] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
m_nonZeroFlops[2] = 61992;
m_hardwareFlops[2] = 127008;
m_matrixKernels[2] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
m_nonZeroFlops[3] = 61992;
m_hardwareFlops[3] = 127008;
m_matrixKernels[3] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
m_nonZeroFlops[4] = 19332;
m_hardwareFlops[4] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[4] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 19332;
m_hardwareFlops[5] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[5] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 6048;
m_hardwareFlops[6] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[6] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 24192;
m_hardwareFlops[7] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[7] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 77328;
m_hardwareFlops[8] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[8] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 77328;
m_hardwareFlops[9] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[9] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 101268;
m_hardwareFlops[10] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[10] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 89802;
m_hardwareFlops[11] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[11] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 89190;
m_hardwareFlops[12] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[12] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 89190;
m_hardwareFlops[13] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[13] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 89802;
m_hardwareFlops[14] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[14] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 101268;
m_hardwareFlops[15] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[15] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 24192;
m_hardwareFlops[16] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[16] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 77328;
m_hardwareFlops[17] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[17] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 77328;
m_hardwareFlops[18] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[18] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 109620;
m_hardwareFlops[19] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[19] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 19656;
m_hardwareFlops[20] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[20] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 109620;
m_hardwareFlops[21] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[21] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 69210;
m_hardwareFlops[22] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[22] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 109242;
m_hardwareFlops[23] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[23] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 101844;
m_hardwareFlops[24] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[24] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 109242;
m_hardwareFlops[25] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[25] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 69210;
m_hardwareFlops[26] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[26] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 101844;
m_hardwareFlops[27] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[27] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 101268;
m_hardwareFlops[28] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[28] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 89802;
m_hardwareFlops[29] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[29] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 89190;
m_hardwareFlops[30] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[30] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 69210;
m_hardwareFlops[31] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[31] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 109242;
m_hardwareFlops[32] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[32] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 101844;
m_hardwareFlops[33] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[33] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 78156;
m_hardwareFlops[34] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[34] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 126396;
m_hardwareFlops[35] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[35] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 85788;
m_hardwareFlops[36] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[36] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 61992;
m_hardwareFlops[37] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[37] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 118008;
m_hardwareFlops[38] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[38] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 118008;
m_hardwareFlops[39] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[39] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 89190;
m_hardwareFlops[40] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[40] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 89802;
m_hardwareFlops[41] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[41] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 101268;
m_hardwareFlops[42] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[42] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 109242;
m_hardwareFlops[43] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[43] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 69210;
m_hardwareFlops[44] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[44] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 101844;
m_hardwareFlops[45] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[45] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 61992;
m_hardwareFlops[46] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[46] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 118008;
m_hardwareFlops[47] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[47] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 118008;
m_hardwareFlops[48] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[48] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 78156;
m_hardwareFlops[49] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[49] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 85788;
m_hardwareFlops[50] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[50] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 126396;
m_hardwareFlops[51] = 127008;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#else
m_matrixKernels[51] = dgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 13608;
m_hardwareFlops[52] = 13608;
m_matrixKernels[52] = dgemm_m84_n9_k9_ldA84_ldB9_ldC84_beta1_pfsigonly;
m_nonZeroFlops[53] = 13608;
m_hardwareFlops[53] = 13608;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = dgemm_m84_n9_k9_ldA84_ldB9_ldC84_beta1_pfsigonly;
#else
m_matrixKernels[53] = dgemm_m84_n9_k9_ldA84_ldB9_ldC84_beta1_pfsigonly;
#endif
#endif

#ifdef SPARSE_SWITCH
m_sparseSwitch[0] = -1; 
m_sparseSwitch[1] = -1; 
m_sparseSwitch[2] = -1; 
m_sparseSwitch[3] = -1; 
m_sparseSwitch[4] = -1; 
m_sparseSwitch[5] = -1; 
m_sparseSwitch[6] = -1; 
m_sparseSwitch[7] = -1; 
m_sparseSwitch[8] = -1; 
m_sparseSwitch[9] = -1; 
m_sparseSwitch[10] = -1; 
m_sparseSwitch[11] = -1; 
m_sparseSwitch[12] = -1; 
m_sparseSwitch[13] = -1; 
m_sparseSwitch[14] = -1; 
m_sparseSwitch[15] = -1; 
m_sparseSwitch[16] = -1; 
m_sparseSwitch[17] = -1; 
m_sparseSwitch[18] = -1; 
m_sparseSwitch[19] = -1; 
m_sparseSwitch[20] = -1; 
m_sparseSwitch[21] = -1; 
m_sparseSwitch[22] = -1; 
m_sparseSwitch[23] = -1; 
m_sparseSwitch[24] = -1; 
m_sparseSwitch[25] = -1; 
m_sparseSwitch[26] = -1; 
m_sparseSwitch[27] = -1; 
m_sparseSwitch[28] = -1; 
m_sparseSwitch[29] = -1; 
m_sparseSwitch[30] = -1; 
m_sparseSwitch[31] = -1; 
m_sparseSwitch[32] = -1; 
m_sparseSwitch[33] = -1; 
m_sparseSwitch[34] = -1; 
m_sparseSwitch[35] = -1; 
m_sparseSwitch[36] = -1; 
m_sparseSwitch[37] = -1; 
m_sparseSwitch[38] = -1; 
m_sparseSwitch[39] = -1; 
m_sparseSwitch[40] = -1; 
m_sparseSwitch[41] = -1; 
m_sparseSwitch[42] = -1; 
m_sparseSwitch[43] = -1; 
m_sparseSwitch[44] = -1; 
m_sparseSwitch[45] = -1; 
m_sparseSwitch[46] = -1; 
m_sparseSwitch[47] = -1; 
m_sparseSwitch[48] = -1; 
m_sparseSwitch[49] = -1; 
m_sparseSwitch[50] = -1; 
m_sparseSwitch[51] = -1; 
m_sparseSwitch[52] = -1; 
m_sparseSwitch[53] = -1; 
m_sparseSwitch[54] = -1; 
m_sparseSwitch[55] = -1; 
m_sparseSwitch[56] = -1; 
m_sparseSwitch[57] = -1; 
m_sparseSwitch[58] = -1; 
m_sparseSwitch[59] = 24; 
#endif

#define STAR_NNZ 24

#endif

#if CONVERGENCE_ORDER==8

#ifdef TIME_KERNEL
m_nonZeroFlops[0] = 26028;
m_matrixKernels[0] = dgemm_m84_n9_k120_ldA84_ldB120_ldC84_beta0_pfsigonly;
m_hardwareFlops[0] = 181440;
m_nonZeroFlops[1] = 58104;
m_matrixKernels[1] = dgemm_m84_n9_k120_ldA84_ldB120_ldC84_beta0_pfsigonly;
m_hardwareFlops[1] = 181440;
m_nonZeroFlops[2] = 61884;
m_matrixKernels[2] = dgemm_m84_n9_k120_ldA84_ldB120_ldC84_beta0_pfsigonly;
m_hardwareFlops[2] = 181440;
m_nonZeroFlops[3] = 4032;
m_hardwareFlops[3] = 4032;
m_matrixKernels[3] = dsparse_starMatrix_m84_n9_k9_ldA84_ldBna8_ldC84_beta1_pfsigonly;
m_nonZeroFlops[4] = 12348;
m_matrixKernels[4] = dgemm_m56_n9_k84_ldA84_ldB84_ldC56_beta0_pfsigonly;
m_hardwareFlops[4] = 84672;
m_nonZeroFlops[5] = 27972;
m_matrixKernels[5] = dgemm_m56_n9_k84_ldA84_ldB84_ldC56_beta0_pfsigonly;
m_hardwareFlops[5] = 84672;
m_nonZeroFlops[6] = 30240;
m_matrixKernels[6] = dgemm_m56_n9_k84_ldA84_ldB84_ldC56_beta0_pfsigonly;
m_hardwareFlops[6] = 84672;
m_nonZeroFlops[7] = 2688;
m_hardwareFlops[7] = 2688;
m_matrixKernels[7] = dsparse_starMatrix_m56_n9_k9_ldA56_ldBna8_ldC56_beta1_pfsigonly;
m_nonZeroFlops[8] = 5292;
m_matrixKernels[8] = dgemm_m36_n9_k56_ldA84_ldB56_ldC36_beta0_pfsigonly;
m_hardwareFlops[8] = 36288;
m_nonZeroFlops[9] = 12096;
m_matrixKernels[9] = dgemm_m36_n9_k56_ldA84_ldB56_ldC36_beta0_pfsigonly;
m_hardwareFlops[9] = 36288;
m_nonZeroFlops[10] = 13356;
m_matrixKernels[10] = dgemm_m36_n9_k56_ldA84_ldB56_ldC36_beta0_pfsigonly;
m_hardwareFlops[10] = 36288;
m_nonZeroFlops[11] = 1680;
m_hardwareFlops[11] = 1680;
m_matrixKernels[11] = dsparse_starMatrix_m35_n9_k9_ldA36_ldBna8_ldC36_beta1_pfsigonly;
m_nonZeroFlops[12] = 1944;
m_matrixKernels[12] = dgemm_m20_n9_k35_ldA84_ldB36_ldC20_beta0_pfsigonly;
m_hardwareFlops[12] = 12600;
m_nonZeroFlops[13] = 4536;
m_matrixKernels[13] = dgemm_m20_n9_k35_ldA84_ldB36_ldC20_beta0_pfsigonly;
m_hardwareFlops[13] = 12600;
m_nonZeroFlops[14] = 5166;
m_matrixKernels[14] = dgemm_m20_n9_k35_ldA84_ldB36_ldC20_beta0_pfsigonly;
m_hardwareFlops[14] = 12600;
m_nonZeroFlops[15] = 960;
m_hardwareFlops[15] = 960;
m_matrixKernels[15] = dsparse_starMatrix_m20_n9_k9_ldA20_ldBna8_ldC20_beta1_pfsigonly;
m_nonZeroFlops[16] = 594;
m_matrixKernels[16] = dgemm_m10_n9_k20_ldA84_ldB20_ldC10_beta0_pfsigonly;
m_hardwareFlops[16] = 3600;
m_nonZeroFlops[17] = 1386;
m_matrixKernels[17] = dgemm_m10_n9_k20_ldA84_ldB20_ldC10_beta0_pfsigonly;
m_hardwareFlops[17] = 3600;
m_nonZeroFlops[18] = 1656;
m_matrixKernels[18] = dgemm_m10_n9_k20_ldA84_ldB20_ldC10_beta0_pfsigonly;
m_hardwareFlops[18] = 3600;
m_nonZeroFlops[19] = 480;
m_hardwareFlops[19] = 480;
m_matrixKernels[19] = dsparse_starMatrix_m10_n9_k9_ldA10_ldBna8_ldC10_beta1_pfsigonly;
m_nonZeroFlops[20] = 126;
m_matrixKernels[20] = dgemm_m4_n9_k10_ldA84_ldB10_ldC4_beta0_pfsigonly;
m_hardwareFlops[20] = 720;
m_nonZeroFlops[21] = 306;
m_matrixKernels[21] = dgemm_m4_n9_k10_ldA84_ldB10_ldC4_beta0_pfsigonly;
m_hardwareFlops[21] = 720;
m_nonZeroFlops[22] = 396;
m_matrixKernels[22] = dgemm_m4_n9_k10_ldA84_ldB10_ldC4_beta0_pfsigonly;
m_hardwareFlops[22] = 720;
m_nonZeroFlops[23] = 192;
m_hardwareFlops[23] = 192;
m_matrixKernels[23] = dsparse_starMatrix_m4_n9_k9_ldA4_ldBna8_ldC4_beta1_pfsigonly;
m_nonZeroFlops[24] = 18;
m_matrixKernels[24] = dgemm_m2_n9_k4_ldA84_ldB4_ldC2_beta0_pfsigonly;
m_hardwareFlops[24] = 144;
m_nonZeroFlops[25] = 36;
m_matrixKernels[25] = dgemm_m2_n9_k4_ldA84_ldB4_ldC2_beta0_pfsigonly;
m_hardwareFlops[25] = 144;
m_nonZeroFlops[26] = 54;
m_matrixKernels[26] = dgemm_m2_n9_k4_ldA84_ldB4_ldC2_beta0_pfsigonly;
m_hardwareFlops[26] = 144;
m_nonZeroFlops[27] = 48;
m_hardwareFlops[27] = 48;
m_matrixKernels[27] = dsparse_starMatrix_m1_n9_k9_ldA2_ldBna8_ldC2_beta1_pfsigonly;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 26028;
m_matrixKernels[0] = dgemm_m120_n9_k84_ldA120_ldB120_ldC120_beta0_pfsigonly;
m_hardwareFlops[0] = 181440;
m_nonZeroFlops[1] = 58104;
m_matrixKernels[1] = dgemm_m120_n9_k84_ldA120_ldB120_ldC120_beta0_pfsigonly;
m_hardwareFlops[1] = 181440;
m_nonZeroFlops[2] = 61884;
m_matrixKernels[2] = dgemm_m120_n9_k84_ldA120_ldB120_ldC120_beta0_pfsigonly;
m_hardwareFlops[2] = 181440;
m_nonZeroFlops[3] = 5760;
m_hardwareFlops[3] = 5760;
m_matrixKernels[3] = dsparse_starMatrix_m120_n9_k9_ldA120_ldBna8_ldC120_beta1_pfsigonly;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 9720;
m_hardwareFlops[0] = 259200;
m_matrixKernels[0] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
m_nonZeroFlops[1] = 35424;
m_hardwareFlops[1] = 259200;
m_matrixKernels[1] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
m_nonZeroFlops[2] = 125928;
m_hardwareFlops[2] = 259200;
m_matrixKernels[2] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
m_nonZeroFlops[3] = 125928;
m_hardwareFlops[3] = 259200;
m_matrixKernels[3] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
m_nonZeroFlops[4] = 34848;
m_hardwareFlops[4] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[4] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 34848;
m_hardwareFlops[5] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[5] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 9720;
m_hardwareFlops[6] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[6] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 43740;
m_hardwareFlops[7] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[7] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 156816;
m_hardwareFlops[8] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[8] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 156816;
m_hardwareFlops[9] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[9] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 207162;
m_hardwareFlops[10] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[10] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 181980;
m_hardwareFlops[11] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[11] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 181620;
m_hardwareFlops[12] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[12] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 181620;
m_hardwareFlops[13] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[13] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 181980;
m_hardwareFlops[14] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[14] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 207162;
m_hardwareFlops[15] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[15] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 43740;
m_hardwareFlops[16] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[16] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 156816;
m_hardwareFlops[17] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[17] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 156816;
m_hardwareFlops[18] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[18] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 224676;
m_hardwareFlops[19] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[19] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 35424;
m_hardwareFlops[20] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[20] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 224676;
m_hardwareFlops[21] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[21] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 140760;
m_hardwareFlops[22] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[22] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 224388;
m_hardwareFlops[23] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[23] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 207144;
m_hardwareFlops[24] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[24] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 224388;
m_hardwareFlops[25] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[25] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 140760;
m_hardwareFlops[26] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[26] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 207144;
m_hardwareFlops[27] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[27] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 207162;
m_hardwareFlops[28] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[28] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 181980;
m_hardwareFlops[29] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[29] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 181620;
m_hardwareFlops[30] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[30] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 140760;
m_hardwareFlops[31] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[31] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 224388;
m_hardwareFlops[32] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[32] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 207144;
m_hardwareFlops[33] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[33] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 158688;
m_hardwareFlops[34] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[34] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 258444;
m_hardwareFlops[35] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[35] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 175284;
m_hardwareFlops[36] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[36] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 125928;
m_hardwareFlops[37] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[37] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 241578;
m_hardwareFlops[38] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[38] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 241578;
m_hardwareFlops[39] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[39] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 181620;
m_hardwareFlops[40] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[40] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 181980;
m_hardwareFlops[41] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[41] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 207162;
m_hardwareFlops[42] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[42] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 224388;
m_hardwareFlops[43] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[43] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 140760;
m_hardwareFlops[44] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[44] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 207144;
m_hardwareFlops[45] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[45] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 125928;
m_hardwareFlops[46] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[46] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 241578;
m_hardwareFlops[47] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[47] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 241578;
m_hardwareFlops[48] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[48] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 158688;
m_hardwareFlops[49] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[49] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 175284;
m_hardwareFlops[50] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[50] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 258444;
m_hardwareFlops[51] = 259200;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#else
m_matrixKernels[51] = dgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 19440;
m_hardwareFlops[52] = 19440;
m_matrixKernels[52] = dgemm_m120_n9_k9_ldA120_ldB9_ldC120_beta1_pfsigonly;
m_nonZeroFlops[53] = 19440;
m_hardwareFlops[53] = 19440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = dgemm_m120_n9_k9_ldA120_ldB9_ldC120_beta1_pfsigonly;
#else
m_matrixKernels[53] = dgemm_m120_n9_k9_ldA120_ldB9_ldC120_beta1_pfsigonly;
#endif
#endif

#ifdef SPARSE_SWITCH
m_sparseSwitch[0] = -1; 
m_sparseSwitch[1] = -1; 
m_sparseSwitch[2] = -1; 
m_sparseSwitch[3] = -1; 
m_sparseSwitch[4] = -1; 
m_sparseSwitch[5] = -1; 
m_sparseSwitch[6] = -1; 
m_sparseSwitch[7] = -1; 
m_sparseSwitch[8] = -1; 
m_sparseSwitch[9] = -1; 
m_sparseSwitch[10] = -1; 
m_sparseSwitch[11] = -1; 
m_sparseSwitch[12] = -1; 
m_sparseSwitch[13] = -1; 
m_sparseSwitch[14] = -1; 
m_sparseSwitch[15] = -1; 
m_sparseSwitch[16] = -1; 
m_sparseSwitch[17] = -1; 
m_sparseSwitch[18] = -1; 
m_sparseSwitch[19] = -1; 
m_sparseSwitch[20] = -1; 
m_sparseSwitch[21] = -1; 
m_sparseSwitch[22] = -1; 
m_sparseSwitch[23] = -1; 
m_sparseSwitch[24] = -1; 
m_sparseSwitch[25] = -1; 
m_sparseSwitch[26] = -1; 
m_sparseSwitch[27] = -1; 
m_sparseSwitch[28] = -1; 
m_sparseSwitch[29] = -1; 
m_sparseSwitch[30] = -1; 
m_sparseSwitch[31] = -1; 
m_sparseSwitch[32] = -1; 
m_sparseSwitch[33] = -1; 
m_sparseSwitch[34] = -1; 
m_sparseSwitch[35] = -1; 
m_sparseSwitch[36] = -1; 
m_sparseSwitch[37] = -1; 
m_sparseSwitch[38] = -1; 
m_sparseSwitch[39] = -1; 
m_sparseSwitch[40] = -1; 
m_sparseSwitch[41] = -1; 
m_sparseSwitch[42] = -1; 
m_sparseSwitch[43] = -1; 
m_sparseSwitch[44] = -1; 
m_sparseSwitch[45] = -1; 
m_sparseSwitch[46] = -1; 
m_sparseSwitch[47] = -1; 
m_sparseSwitch[48] = -1; 
m_sparseSwitch[49] = -1; 
m_sparseSwitch[50] = -1; 
m_sparseSwitch[51] = -1; 
m_sparseSwitch[52] = -1; 
m_sparseSwitch[53] = -1; 
m_sparseSwitch[54] = -1; 
m_sparseSwitch[55] = -1; 
m_sparseSwitch[56] = -1; 
m_sparseSwitch[57] = -1; 
m_sparseSwitch[58] = -1; 
m_sparseSwitch[59] = 24; 
#endif

#define STAR_NNZ 24

#endif

