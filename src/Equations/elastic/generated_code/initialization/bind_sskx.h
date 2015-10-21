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
// @date 2015-10-20 16:05:15.386273
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
#if ALIGNMENT!=64
#error alignment-architecture mismatch
#endif

#if CONVERGENCE_ORDER==2

#ifdef TIME_KERNEL
m_nonZeroFlops[0] = 18;
m_matrixKernels[0] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[0] = 1152;
m_nonZeroFlops[1] = 36;
m_matrixKernels[1] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[1] = 1152;
m_nonZeroFlops[2] = 54;
m_matrixKernels[2] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[2] = 1152;
m_nonZeroFlops[3] = 48;
m_matrixKernels[3] = sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly;
m_hardwareFlops[3] = 2592;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 18;
m_matrixKernels[0] = sgemm_m16_n9_k1_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[0] = 288;
m_nonZeroFlops[1] = 36;
m_matrixKernels[1] = sgemm_m16_n9_k1_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[1] = 288;
m_nonZeroFlops[2] = 54;
m_matrixKernels[2] = sgemm_m16_n9_k1_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[2] = 288;
m_nonZeroFlops[3] = 192;
m_matrixKernels[3] = sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly;
m_hardwareFlops[3] = 2592;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 108;
m_hardwareFlops[0] = 1152;
m_matrixKernels[0] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_nonZeroFlops[1] = 144;
m_hardwareFlops[1] = 1152;
m_matrixKernels[1] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_nonZeroFlops[2] = 180;
m_hardwareFlops[2] = 1152;
m_matrixKernels[2] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_nonZeroFlops[3] = 180;
m_hardwareFlops[3] = 1152;
m_matrixKernels[3] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_nonZeroFlops[4] = 144;
m_hardwareFlops[4] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[4] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 144;
m_hardwareFlops[5] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[5] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 108;
m_hardwareFlops[6] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[6] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 162;
m_hardwareFlops[7] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[7] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 216;
m_hardwareFlops[8] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[8] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 216;
m_hardwareFlops[9] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[9] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 252;
m_hardwareFlops[10] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[10] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 234;
m_hardwareFlops[11] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[11] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 234;
m_hardwareFlops[12] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[12] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 234;
m_hardwareFlops[13] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[13] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 234;
m_hardwareFlops[14] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[14] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 252;
m_hardwareFlops[15] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[15] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 162;
m_hardwareFlops[16] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[16] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 216;
m_hardwareFlops[17] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[17] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 216;
m_hardwareFlops[18] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[18] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 252;
m_hardwareFlops[19] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[19] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 144;
m_hardwareFlops[20] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[20] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 252;
m_hardwareFlops[21] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[21] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 198;
m_hardwareFlops[22] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[22] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 252;
m_hardwareFlops[23] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[23] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 252;
m_hardwareFlops[24] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[24] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 252;
m_hardwareFlops[25] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[25] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 198;
m_hardwareFlops[26] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[26] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 252;
m_hardwareFlops[27] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[27] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 252;
m_hardwareFlops[28] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[28] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 234;
m_hardwareFlops[29] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[29] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 234;
m_hardwareFlops[30] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[30] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 198;
m_hardwareFlops[31] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[31] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 252;
m_hardwareFlops[32] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[32] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 252;
m_hardwareFlops[33] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[33] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 216;
m_hardwareFlops[34] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[34] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 288;
m_hardwareFlops[35] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[35] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 216;
m_hardwareFlops[36] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[36] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 180;
m_hardwareFlops[37] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[37] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 270;
m_hardwareFlops[38] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[38] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 270;
m_hardwareFlops[39] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[39] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 234;
m_hardwareFlops[40] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[40] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 234;
m_hardwareFlops[41] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[41] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 252;
m_hardwareFlops[42] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[42] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 252;
m_hardwareFlops[43] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[43] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 198;
m_hardwareFlops[44] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[44] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 252;
m_hardwareFlops[45] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[45] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 180;
m_hardwareFlops[46] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[46] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 270;
m_hardwareFlops[47] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[47] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 270;
m_hardwareFlops[48] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[48] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 216;
m_hardwareFlops[49] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[49] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 216;
m_hardwareFlops[50] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[50] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 288;
m_hardwareFlops[51] = 1152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[51] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 648;
m_hardwareFlops[52] = 2592;
m_matrixKernels[52] = sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly;
m_nonZeroFlops[53] = 648;
m_hardwareFlops[53] = 2592;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly;
#else
m_matrixKernels[53] = sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly;
#endif
m_nonZeroFlops[54] = 648;
m_hardwareFlops[54] = 2592;
#ifdef ENABLE_STREAM_MATRIX_PREFETCH
m_matrixKernels[54] = sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly;
#else
m_matrixKernels[54] = sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly;
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
m_sparseSwitch[59] = -1; 
#endif

#define STAR_NNZ 81

#endif

#if CONVERGENCE_ORDER==3

#ifdef TIME_KERNEL
m_nonZeroFlops[0] = 126;
m_matrixKernels[0] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[0] = 2880;
m_nonZeroFlops[1] = 306;
m_matrixKernels[1] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[1] = 2880;
m_nonZeroFlops[2] = 396;
m_matrixKernels[2] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[2] = 2880;
m_nonZeroFlops[3] = 192;
m_matrixKernels[3] = sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly;
m_hardwareFlops[3] = 2592;
m_nonZeroFlops[4] = 18;
m_matrixKernels[4] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[4] = 1152;
m_nonZeroFlops[5] = 36;
m_matrixKernels[5] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[5] = 1152;
m_nonZeroFlops[6] = 54;
m_matrixKernels[6] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[6] = 1152;
m_nonZeroFlops[7] = 48;
m_matrixKernels[7] = sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly;
m_hardwareFlops[7] = 2592;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 126;
m_matrixKernels[0] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[0] = 1152;
m_nonZeroFlops[1] = 306;
m_matrixKernels[1] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[1] = 1152;
m_nonZeroFlops[2] = 396;
m_matrixKernels[2] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[2] = 1152;
m_nonZeroFlops[3] = 480;
m_matrixKernels[3] = sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly;
m_hardwareFlops[3] = 2592;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 360;
m_hardwareFlops[0] = 2880;
m_matrixKernels[0] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_nonZeroFlops[1] = 612;
m_hardwareFlops[1] = 2880;
m_matrixKernels[1] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_nonZeroFlops[2] = 972;
m_hardwareFlops[2] = 2880;
m_matrixKernels[2] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_nonZeroFlops[3] = 972;
m_hardwareFlops[3] = 2880;
m_matrixKernels[3] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_nonZeroFlops[4] = 612;
m_hardwareFlops[4] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[4] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 612;
m_hardwareFlops[5] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[5] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 360;
m_hardwareFlops[6] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[6] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 720;
m_hardwareFlops[7] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[7] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 1224;
m_hardwareFlops[8] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[8] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 1224;
m_hardwareFlops[9] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[9] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 1458;
m_hardwareFlops[10] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[10] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 1368;
m_hardwareFlops[11] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[11] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 1368;
m_hardwareFlops[12] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[12] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 1368;
m_hardwareFlops[13] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[13] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 1368;
m_hardwareFlops[14] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[14] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 1458;
m_hardwareFlops[15] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[15] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 720;
m_hardwareFlops[16] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[16] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 1224;
m_hardwareFlops[17] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[17] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 1224;
m_hardwareFlops[18] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[18] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 1512;
m_hardwareFlops[19] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[19] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 612;
m_hardwareFlops[20] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[20] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 1512;
m_hardwareFlops[21] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[21] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 1098;
m_hardwareFlops[22] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[22] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 1512;
m_hardwareFlops[23] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[23] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 1512;
m_hardwareFlops[24] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[24] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 1512;
m_hardwareFlops[25] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[25] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 1098;
m_hardwareFlops[26] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[26] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 1512;
m_hardwareFlops[27] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[27] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 1458;
m_hardwareFlops[28] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[28] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 1368;
m_hardwareFlops[29] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[29] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 1368;
m_hardwareFlops[30] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[30] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 1098;
m_hardwareFlops[31] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[31] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 1512;
m_hardwareFlops[32] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[32] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 1512;
m_hardwareFlops[33] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[33] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 1224;
m_hardwareFlops[34] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[34] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 1764;
m_hardwareFlops[35] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[35] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 1260;
m_hardwareFlops[36] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[36] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 972;
m_hardwareFlops[37] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[37] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 1674;
m_hardwareFlops[38] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[38] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 1674;
m_hardwareFlops[39] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[39] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 1368;
m_hardwareFlops[40] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[40] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 1368;
m_hardwareFlops[41] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[41] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 1458;
m_hardwareFlops[42] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[42] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 1512;
m_hardwareFlops[43] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[43] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 1098;
m_hardwareFlops[44] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[44] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 1512;
m_hardwareFlops[45] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[45] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 972;
m_hardwareFlops[46] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[46] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 1674;
m_hardwareFlops[47] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[47] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 1674;
m_hardwareFlops[48] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[48] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 1224;
m_hardwareFlops[49] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[49] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 1260;
m_hardwareFlops[50] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[50] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 1764;
m_hardwareFlops[51] = 2880;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[51] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 1620;
m_hardwareFlops[52] = 2592;
m_matrixKernels[52] = sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly;
m_nonZeroFlops[53] = 1620;
m_hardwareFlops[53] = 2592;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly;
#else
m_matrixKernels[53] = sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly;
#endif
m_nonZeroFlops[54] = 1620;
m_hardwareFlops[54] = 2592;
#ifdef ENABLE_STREAM_MATRIX_PREFETCH
m_matrixKernels[54] = sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly;
#else
m_matrixKernels[54] = sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly;
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
m_sparseSwitch[59] = -1; 
#endif

#define STAR_NNZ 81

#endif

#if CONVERGENCE_ORDER==4

#ifdef TIME_KERNEL
m_nonZeroFlops[0] = 594;
m_matrixKernels[0] = sgemm_m16_n9_k20_ldA16_ldB32_ldC16_beta0_pfsigonly;
m_hardwareFlops[0] = 5760;
m_nonZeroFlops[1] = 1386;
m_matrixKernels[1] = sgemm_m16_n9_k20_ldA16_ldB32_ldC16_beta0_pfsigonly;
m_hardwareFlops[1] = 5760;
m_nonZeroFlops[2] = 1656;
m_matrixKernels[2] = sgemm_m16_n9_k20_ldA16_ldB32_ldC16_beta0_pfsigonly;
m_hardwareFlops[2] = 5760;
m_nonZeroFlops[3] = 480;
m_matrixKernels[3] = sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly;
m_hardwareFlops[3] = 2592;
m_nonZeroFlops[4] = 126;
m_matrixKernels[4] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[4] = 2880;
m_nonZeroFlops[5] = 306;
m_matrixKernels[5] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[5] = 2880;
m_nonZeroFlops[6] = 396;
m_matrixKernels[6] = sgemm_m16_n9_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[6] = 2880;
m_nonZeroFlops[7] = 192;
m_matrixKernels[7] = sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly;
m_hardwareFlops[7] = 2592;
m_nonZeroFlops[8] = 18;
m_matrixKernels[8] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[8] = 1152;
m_nonZeroFlops[9] = 36;
m_matrixKernels[9] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[9] = 1152;
m_nonZeroFlops[10] = 54;
m_matrixKernels[10] = sgemm_m16_n9_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[10] = 1152;
m_nonZeroFlops[11] = 48;
m_matrixKernels[11] = sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly;
m_hardwareFlops[11] = 2592;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 594;
m_matrixKernels[0] = sgemm_m32_n9_k10_ldA32_ldB32_ldC32_beta0_pfsigonly;
m_hardwareFlops[0] = 5760;
m_nonZeroFlops[1] = 1386;
m_matrixKernels[1] = sgemm_m32_n9_k10_ldA32_ldB32_ldC32_beta0_pfsigonly;
m_hardwareFlops[1] = 5760;
m_nonZeroFlops[2] = 1656;
m_matrixKernels[2] = sgemm_m32_n9_k10_ldA32_ldB32_ldC32_beta0_pfsigonly;
m_hardwareFlops[2] = 5760;
m_nonZeroFlops[3] = 960;
m_matrixKernels[3] = sgemm_m32_n9_k9_ldA32_ldB9_ldC32_beta1_pfsigonly;
m_hardwareFlops[3] = 5184;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 900;
m_hardwareFlops[0] = 11520;
m_matrixKernels[0] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
m_nonZeroFlops[1] = 1872;
m_hardwareFlops[1] = 11520;
m_matrixKernels[1] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
m_nonZeroFlops[2] = 3672;
m_hardwareFlops[2] = 11520;
m_matrixKernels[2] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
m_nonZeroFlops[3] = 3672;
m_hardwareFlops[3] = 11520;
m_matrixKernels[3] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
m_nonZeroFlops[4] = 1872;
m_hardwareFlops[4] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[4] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 1872;
m_hardwareFlops[5] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[5] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 900;
m_hardwareFlops[6] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[6] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 2250;
m_hardwareFlops[7] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[7] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 4680;
m_hardwareFlops[8] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[8] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 4680;
m_hardwareFlops[9] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[9] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 5760;
m_hardwareFlops[10] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[10] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 5148;
m_hardwareFlops[11] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[11] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 5310;
m_hardwareFlops[12] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[12] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 5310;
m_hardwareFlops[13] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[13] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 5148;
m_hardwareFlops[14] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[14] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 5760;
m_hardwareFlops[15] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[15] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 2250;
m_hardwareFlops[16] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[16] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 4680;
m_hardwareFlops[17] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[17] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 4680;
m_hardwareFlops[18] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[18] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 6084;
m_hardwareFlops[19] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[19] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 1872;
m_hardwareFlops[20] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[20] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 6084;
m_hardwareFlops[21] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[21] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 4176;
m_hardwareFlops[22] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[22] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 6120;
m_hardwareFlops[23] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[23] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 5850;
m_hardwareFlops[24] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[24] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 6120;
m_hardwareFlops[25] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[25] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 4176;
m_hardwareFlops[26] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[26] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 5850;
m_hardwareFlops[27] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[27] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 5760;
m_hardwareFlops[28] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[28] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 5148;
m_hardwareFlops[29] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[29] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 5310;
m_hardwareFlops[30] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[30] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 4176;
m_hardwareFlops[31] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[31] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 6120;
m_hardwareFlops[32] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[32] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 5850;
m_hardwareFlops[33] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[33] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 4644;
m_hardwareFlops[34] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[34] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 7092;
m_hardwareFlops[35] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[35] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 4932;
m_hardwareFlops[36] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[36] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 3672;
m_hardwareFlops[37] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[37] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 6678;
m_hardwareFlops[38] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[38] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 6678;
m_hardwareFlops[39] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[39] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 5310;
m_hardwareFlops[40] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[40] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 5148;
m_hardwareFlops[41] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[41] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 5760;
m_hardwareFlops[42] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[42] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 6120;
m_hardwareFlops[43] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[43] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 4176;
m_hardwareFlops[44] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[44] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 5850;
m_hardwareFlops[45] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[45] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 3672;
m_hardwareFlops[46] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[46] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 6678;
m_hardwareFlops[47] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[47] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 6678;
m_hardwareFlops[48] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[48] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 4644;
m_hardwareFlops[49] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[49] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 4932;
m_hardwareFlops[50] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[50] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 7092;
m_hardwareFlops[51] = 11520;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[51] = sgemm_m32_n9_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 3240;
m_hardwareFlops[52] = 5184;
m_matrixKernels[52] = sgemm_m32_n9_k9_ldA32_ldB9_ldC32_beta1_pfsigonly;
m_nonZeroFlops[53] = 3240;
m_hardwareFlops[53] = 5184;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = sgemm_m32_n9_k9_ldA32_ldB9_ldC32_beta1_pfsigonly;
#else
m_matrixKernels[53] = sgemm_m32_n9_k9_ldA32_ldB9_ldC32_beta1_pfsigonly;
#endif
m_nonZeroFlops[54] = 3240;
m_hardwareFlops[54] = 5184;
#ifdef ENABLE_STREAM_MATRIX_PREFETCH
m_matrixKernels[54] = sgemm_m32_n9_k9_ldA32_ldB9_ldC32_beta1_pfsigonly;
#else
m_matrixKernels[54] = sgemm_m32_n9_k9_ldA32_ldB9_ldC32_beta1_pfsigonly;
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
m_sparseSwitch[59] = -1; 
#endif

#define STAR_NNZ 81

#endif

#if CONVERGENCE_ORDER==5

#ifdef TIME_KERNEL
m_nonZeroFlops[0] = 1944;
m_matrixKernels[0] = sgemm_m32_n9_k35_ldA32_ldB48_ldC32_beta0_pfsigonly;
m_hardwareFlops[0] = 20160;
m_nonZeroFlops[1] = 4536;
m_matrixKernels[1] = sgemm_m32_n9_k35_ldA32_ldB48_ldC32_beta0_pfsigonly;
m_hardwareFlops[1] = 20160;
m_nonZeroFlops[2] = 5166;
m_matrixKernels[2] = sgemm_m32_n9_k35_ldA32_ldB48_ldC32_beta0_pfsigonly;
m_hardwareFlops[2] = 20160;
m_nonZeroFlops[3] = 960;
m_matrixKernels[3] = sgemm_m32_n9_k9_ldA32_ldB9_ldC32_beta1_pfsigonly;
m_hardwareFlops[3] = 5184;
m_nonZeroFlops[4] = 594;
m_matrixKernels[4] = sgemm_m16_n9_k20_ldA32_ldB32_ldC16_beta0_pfsigonly;
m_hardwareFlops[4] = 5760;
m_nonZeroFlops[5] = 1386;
m_matrixKernels[5] = sgemm_m16_n9_k20_ldA32_ldB32_ldC16_beta0_pfsigonly;
m_hardwareFlops[5] = 5760;
m_nonZeroFlops[6] = 1656;
m_matrixKernels[6] = sgemm_m16_n9_k20_ldA32_ldB32_ldC16_beta0_pfsigonly;
m_hardwareFlops[6] = 5760;
m_nonZeroFlops[7] = 480;
m_matrixKernels[7] = sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly;
m_hardwareFlops[7] = 2592;
m_nonZeroFlops[8] = 126;
m_matrixKernels[8] = sgemm_m16_n9_k10_ldA32_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[8] = 2880;
m_nonZeroFlops[9] = 306;
m_matrixKernels[9] = sgemm_m16_n9_k10_ldA32_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[9] = 2880;
m_nonZeroFlops[10] = 396;
m_matrixKernels[10] = sgemm_m16_n9_k10_ldA32_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[10] = 2880;
m_nonZeroFlops[11] = 192;
m_matrixKernels[11] = sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly;
m_hardwareFlops[11] = 2592;
m_nonZeroFlops[12] = 18;
m_matrixKernels[12] = sgemm_m16_n9_k4_ldA32_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[12] = 1152;
m_nonZeroFlops[13] = 36;
m_matrixKernels[13] = sgemm_m16_n9_k4_ldA32_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[13] = 1152;
m_nonZeroFlops[14] = 54;
m_matrixKernels[14] = sgemm_m16_n9_k4_ldA32_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[14] = 1152;
m_nonZeroFlops[15] = 48;
m_matrixKernels[15] = sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly;
m_hardwareFlops[15] = 2592;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 1944;
m_matrixKernels[0] = sgemm_m48_n9_k20_ldA48_ldB48_ldC48_beta0_pfsigonly;
m_hardwareFlops[0] = 17280;
m_nonZeroFlops[1] = 4536;
m_matrixKernels[1] = sgemm_m48_n9_k20_ldA48_ldB48_ldC48_beta0_pfsigonly;
m_hardwareFlops[1] = 17280;
m_nonZeroFlops[2] = 5166;
m_matrixKernels[2] = sgemm_m48_n9_k20_ldA48_ldB48_ldC48_beta0_pfsigonly;
m_hardwareFlops[2] = 17280;
m_nonZeroFlops[3] = 1680;
m_matrixKernels[3] = sgemm_m48_n9_k9_ldA48_ldB9_ldC48_beta1_pfsigonly;
m_hardwareFlops[3] = 7776;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 1890;
m_hardwareFlops[0] = 30240;
m_matrixKernels[0] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
m_nonZeroFlops[1] = 4662;
m_hardwareFlops[1] = 30240;
m_matrixKernels[1] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
m_nonZeroFlops[2] = 10962;
m_hardwareFlops[2] = 30240;
m_matrixKernels[2] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
m_nonZeroFlops[3] = 10962;
m_hardwareFlops[3] = 30240;
m_matrixKernels[3] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
m_nonZeroFlops[4] = 4626;
m_hardwareFlops[4] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[4] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 4626;
m_hardwareFlops[5] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[5] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 1890;
m_hardwareFlops[6] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[6] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 5670;
m_hardwareFlops[7] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[7] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 13878;
m_hardwareFlops[8] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[8] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 13878;
m_hardwareFlops[9] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[9] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 17622;
m_hardwareFlops[10] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[10] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 15480;
m_hardwareFlops[11] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[11] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 15768;
m_hardwareFlops[12] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[12] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 15768;
m_hardwareFlops[13] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[13] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 15480;
m_hardwareFlops[14] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[14] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 17622;
m_hardwareFlops[15] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[15] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 5670;
m_hardwareFlops[16] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[16] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 13878;
m_hardwareFlops[17] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[17] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 13878;
m_hardwareFlops[18] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[18] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 18810;
m_hardwareFlops[19] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[19] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 4662;
m_hardwareFlops[20] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[20] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 18810;
m_hardwareFlops[21] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[21] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 12312;
m_hardwareFlops[22] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[22] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 18702;
m_hardwareFlops[23] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[23] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 17640;
m_hardwareFlops[24] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[24] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 18702;
m_hardwareFlops[25] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[25] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 12312;
m_hardwareFlops[26] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[26] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 17640;
m_hardwareFlops[27] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[27] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 17622;
m_hardwareFlops[28] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[28] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 15480;
m_hardwareFlops[29] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[29] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 15768;
m_hardwareFlops[30] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[30] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 12312;
m_hardwareFlops[31] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[31] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 18702;
m_hardwareFlops[32] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[32] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 17640;
m_hardwareFlops[33] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[33] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 13914;
m_hardwareFlops[34] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[34] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 21834;
m_hardwareFlops[35] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[35] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 14958;
m_hardwareFlops[36] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[36] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 10962;
m_hardwareFlops[37] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[37] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 20358;
m_hardwareFlops[38] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[38] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 20358;
m_hardwareFlops[39] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[39] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 15768;
m_hardwareFlops[40] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[40] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 15480;
m_hardwareFlops[41] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[41] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 17622;
m_hardwareFlops[42] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[42] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 18702;
m_hardwareFlops[43] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[43] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 12312;
m_hardwareFlops[44] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[44] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 17640;
m_hardwareFlops[45] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[45] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 10962;
m_hardwareFlops[46] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[46] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 20358;
m_hardwareFlops[47] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[47] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 20358;
m_hardwareFlops[48] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[48] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 13914;
m_hardwareFlops[49] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[49] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 14958;
m_hardwareFlops[50] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[50] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 21834;
m_hardwareFlops[51] = 30240;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[51] = sgemm_m48_n9_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 5670;
m_hardwareFlops[52] = 7776;
m_matrixKernels[52] = sgemm_m48_n9_k9_ldA48_ldB9_ldC48_beta1_pfsigonly;
m_nonZeroFlops[53] = 5670;
m_hardwareFlops[53] = 7776;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = sgemm_m48_n9_k9_ldA48_ldB9_ldC48_beta1_pfsigonly;
#else
m_matrixKernels[53] = sgemm_m48_n9_k9_ldA48_ldB9_ldC48_beta1_pfsigonly;
#endif
m_nonZeroFlops[54] = 5670;
m_hardwareFlops[54] = 7776;
#ifdef ENABLE_STREAM_MATRIX_PREFETCH
m_matrixKernels[54] = sgemm_m48_n9_k9_ldA48_ldB9_ldC48_beta1_pfsigonly;
#else
m_matrixKernels[54] = sgemm_m48_n9_k9_ldA48_ldB9_ldC48_beta1_pfsigonly;
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
m_sparseSwitch[59] = -1; 
#endif

#define STAR_NNZ 81

#endif

#if CONVERGENCE_ORDER==6

#ifdef TIME_KERNEL
m_nonZeroFlops[0] = 5292;
m_matrixKernels[0] = sgemm_m48_n9_k56_ldA48_ldB64_ldC48_beta0_pfsigonly;
m_hardwareFlops[0] = 48384;
m_nonZeroFlops[1] = 12096;
m_matrixKernels[1] = sgemm_m48_n9_k56_ldA48_ldB64_ldC48_beta0_pfsigonly;
m_hardwareFlops[1] = 48384;
m_nonZeroFlops[2] = 13356;
m_matrixKernels[2] = sgemm_m48_n9_k56_ldA48_ldB64_ldC48_beta0_pfsigonly;
m_hardwareFlops[2] = 48384;
m_nonZeroFlops[3] = 1680;
m_matrixKernels[3] = sgemm_m48_n9_k9_ldA48_ldB9_ldC48_beta1_pfsigonly;
m_hardwareFlops[3] = 7776;
m_nonZeroFlops[4] = 1944;
m_matrixKernels[4] = sgemm_m32_n9_k35_ldA48_ldB48_ldC32_beta0_pfsigonly;
m_hardwareFlops[4] = 20160;
m_nonZeroFlops[5] = 4536;
m_matrixKernels[5] = sgemm_m32_n9_k35_ldA48_ldB48_ldC32_beta0_pfsigonly;
m_hardwareFlops[5] = 20160;
m_nonZeroFlops[6] = 5166;
m_matrixKernels[6] = sgemm_m32_n9_k35_ldA48_ldB48_ldC32_beta0_pfsigonly;
m_hardwareFlops[6] = 20160;
m_nonZeroFlops[7] = 960;
m_matrixKernels[7] = sgemm_m32_n9_k9_ldA32_ldB9_ldC32_beta1_pfsigonly;
m_hardwareFlops[7] = 5184;
m_nonZeroFlops[8] = 594;
m_matrixKernels[8] = sgemm_m16_n9_k20_ldA48_ldB32_ldC16_beta0_pfsigonly;
m_hardwareFlops[8] = 5760;
m_nonZeroFlops[9] = 1386;
m_matrixKernels[9] = sgemm_m16_n9_k20_ldA48_ldB32_ldC16_beta0_pfsigonly;
m_hardwareFlops[9] = 5760;
m_nonZeroFlops[10] = 1656;
m_matrixKernels[10] = sgemm_m16_n9_k20_ldA48_ldB32_ldC16_beta0_pfsigonly;
m_hardwareFlops[10] = 5760;
m_nonZeroFlops[11] = 480;
m_matrixKernels[11] = sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly;
m_hardwareFlops[11] = 2592;
m_nonZeroFlops[12] = 126;
m_matrixKernels[12] = sgemm_m16_n9_k10_ldA48_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[12] = 2880;
m_nonZeroFlops[13] = 306;
m_matrixKernels[13] = sgemm_m16_n9_k10_ldA48_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[13] = 2880;
m_nonZeroFlops[14] = 396;
m_matrixKernels[14] = sgemm_m16_n9_k10_ldA48_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[14] = 2880;
m_nonZeroFlops[15] = 192;
m_matrixKernels[15] = sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly;
m_hardwareFlops[15] = 2592;
m_nonZeroFlops[16] = 18;
m_matrixKernels[16] = sgemm_m16_n9_k4_ldA48_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[16] = 1152;
m_nonZeroFlops[17] = 36;
m_matrixKernels[17] = sgemm_m16_n9_k4_ldA48_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[17] = 1152;
m_nonZeroFlops[18] = 54;
m_matrixKernels[18] = sgemm_m16_n9_k4_ldA48_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[18] = 1152;
m_nonZeroFlops[19] = 48;
m_matrixKernels[19] = sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly;
m_hardwareFlops[19] = 2592;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 5292;
m_matrixKernels[0] = sgemm_m64_n9_k35_ldA64_ldB64_ldC64_beta0_pfsigonly;
m_hardwareFlops[0] = 40320;
m_nonZeroFlops[1] = 12096;
m_matrixKernels[1] = sgemm_m64_n9_k35_ldA64_ldB64_ldC64_beta0_pfsigonly;
m_hardwareFlops[1] = 40320;
m_nonZeroFlops[2] = 13356;
m_matrixKernels[2] = sgemm_m64_n9_k35_ldA64_ldB64_ldC64_beta0_pfsigonly;
m_hardwareFlops[2] = 40320;
m_nonZeroFlops[3] = 2688;
m_matrixKernels[3] = sgemm_m64_n9_k9_ldA64_ldB9_ldC64_beta1_pfsigonly;
m_hardwareFlops[3] = 10368;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 3528;
m_hardwareFlops[0] = 64512;
m_matrixKernels[0] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
m_nonZeroFlops[1] = 10080;
m_hardwareFlops[1] = 64512;
m_matrixKernels[1] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
m_nonZeroFlops[2] = 27720;
m_hardwareFlops[2] = 64512;
m_matrixKernels[2] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
m_nonZeroFlops[3] = 27720;
m_hardwareFlops[3] = 64512;
m_matrixKernels[3] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
m_nonZeroFlops[4] = 9936;
m_hardwareFlops[4] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[4] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 9936;
m_hardwareFlops[5] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[5] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 3528;
m_hardwareFlops[6] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[6] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 12348;
m_hardwareFlops[7] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[7] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 34776;
m_hardwareFlops[8] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[8] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 34776;
m_hardwareFlops[9] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[9] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 45324;
m_hardwareFlops[10] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[10] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 39852;
m_hardwareFlops[11] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[11] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 39852;
m_hardwareFlops[12] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[12] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 39852;
m_hardwareFlops[13] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[13] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 39852;
m_hardwareFlops[14] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[14] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 45324;
m_hardwareFlops[15] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[15] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 12348;
m_hardwareFlops[16] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[16] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 34776;
m_hardwareFlops[17] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[17] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 34776;
m_hardwareFlops[18] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[18] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 48564;
m_hardwareFlops[19] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[19] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 10080;
m_hardwareFlops[20] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[20] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 48564;
m_hardwareFlops[21] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[21] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 30996;
m_hardwareFlops[22] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[22] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 48348;
m_hardwareFlops[23] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[23] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 45162;
m_hardwareFlops[24] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[24] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 48348;
m_hardwareFlops[25] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[25] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 30996;
m_hardwareFlops[26] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[26] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 45162;
m_hardwareFlops[27] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[27] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 45324;
m_hardwareFlops[28] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[28] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 39852;
m_hardwareFlops[29] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[29] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 39852;
m_hardwareFlops[30] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[30] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 30996;
m_hardwareFlops[31] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[31] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 48348;
m_hardwareFlops[32] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[32] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 45162;
m_hardwareFlops[33] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[33] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 35028;
m_hardwareFlops[34] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[34] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 56232;
m_hardwareFlops[35] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[35] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 38160;
m_hardwareFlops[36] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[36] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 27720;
m_hardwareFlops[37] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[37] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 52290;
m_hardwareFlops[38] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[38] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 52290;
m_hardwareFlops[39] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[39] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 39852;
m_hardwareFlops[40] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[40] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 39852;
m_hardwareFlops[41] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[41] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 45324;
m_hardwareFlops[42] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[42] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 48348;
m_hardwareFlops[43] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[43] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 30996;
m_hardwareFlops[44] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[44] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 45162;
m_hardwareFlops[45] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[45] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 27720;
m_hardwareFlops[46] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[46] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 52290;
m_hardwareFlops[47] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[47] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 52290;
m_hardwareFlops[48] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[48] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 35028;
m_hardwareFlops[49] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[49] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 38160;
m_hardwareFlops[50] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[50] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 56232;
m_hardwareFlops[51] = 64512;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[51] = sgemm_m64_n9_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 9072;
m_hardwareFlops[52] = 10368;
m_matrixKernels[52] = sgemm_m64_n9_k9_ldA64_ldB9_ldC64_beta1_pfsigonly;
m_nonZeroFlops[53] = 9072;
m_hardwareFlops[53] = 10368;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = sgemm_m64_n9_k9_ldA64_ldB9_ldC64_beta1_pfsigonly;
#else
m_matrixKernels[53] = sgemm_m64_n9_k9_ldA64_ldB9_ldC64_beta1_pfsigonly;
#endif
m_nonZeroFlops[54] = 9072;
m_hardwareFlops[54] = 10368;
#ifdef ENABLE_STREAM_MATRIX_PREFETCH
m_matrixKernels[54] = sgemm_m64_n9_k9_ldA64_ldB9_ldC64_beta1_pfsigonly;
#else
m_matrixKernels[54] = sgemm_m64_n9_k9_ldA64_ldB9_ldC64_beta1_pfsigonly;
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
m_sparseSwitch[59] = -1; 
#endif

#define STAR_NNZ 81

#endif

#if CONVERGENCE_ORDER==7

#ifdef TIME_KERNEL
m_nonZeroFlops[0] = 12348;
m_matrixKernels[0] = sgemm_m64_n9_k84_ldA64_ldB96_ldC64_beta0_pfsigonly;
m_hardwareFlops[0] = 96768;
m_nonZeroFlops[1] = 27972;
m_matrixKernels[1] = sgemm_m64_n9_k84_ldA64_ldB96_ldC64_beta0_pfsigonly;
m_hardwareFlops[1] = 96768;
m_nonZeroFlops[2] = 30240;
m_matrixKernels[2] = sgemm_m64_n9_k84_ldA64_ldB96_ldC64_beta0_pfsigonly;
m_hardwareFlops[2] = 96768;
m_nonZeroFlops[3] = 2688;
m_matrixKernels[3] = sgemm_m64_n9_k9_ldA64_ldB9_ldC64_beta1_pfsigonly;
m_hardwareFlops[3] = 10368;
m_nonZeroFlops[4] = 5292;
m_matrixKernels[4] = sgemm_m48_n9_k56_ldA64_ldB64_ldC48_beta0_pfsigonly;
m_hardwareFlops[4] = 48384;
m_nonZeroFlops[5] = 12096;
m_matrixKernels[5] = sgemm_m48_n9_k56_ldA64_ldB64_ldC48_beta0_pfsigonly;
m_hardwareFlops[5] = 48384;
m_nonZeroFlops[6] = 13356;
m_matrixKernels[6] = sgemm_m48_n9_k56_ldA64_ldB64_ldC48_beta0_pfsigonly;
m_hardwareFlops[6] = 48384;
m_nonZeroFlops[7] = 1680;
m_matrixKernels[7] = sgemm_m48_n9_k9_ldA48_ldB9_ldC48_beta1_pfsigonly;
m_hardwareFlops[7] = 7776;
m_nonZeroFlops[8] = 1944;
m_matrixKernels[8] = sgemm_m32_n9_k35_ldA64_ldB48_ldC32_beta0_pfsigonly;
m_hardwareFlops[8] = 20160;
m_nonZeroFlops[9] = 4536;
m_matrixKernels[9] = sgemm_m32_n9_k35_ldA64_ldB48_ldC32_beta0_pfsigonly;
m_hardwareFlops[9] = 20160;
m_nonZeroFlops[10] = 5166;
m_matrixKernels[10] = sgemm_m32_n9_k35_ldA64_ldB48_ldC32_beta0_pfsigonly;
m_hardwareFlops[10] = 20160;
m_nonZeroFlops[11] = 960;
m_matrixKernels[11] = sgemm_m32_n9_k9_ldA32_ldB9_ldC32_beta1_pfsigonly;
m_hardwareFlops[11] = 5184;
m_nonZeroFlops[12] = 594;
m_matrixKernels[12] = sgemm_m16_n9_k20_ldA64_ldB32_ldC16_beta0_pfsigonly;
m_hardwareFlops[12] = 5760;
m_nonZeroFlops[13] = 1386;
m_matrixKernels[13] = sgemm_m16_n9_k20_ldA64_ldB32_ldC16_beta0_pfsigonly;
m_hardwareFlops[13] = 5760;
m_nonZeroFlops[14] = 1656;
m_matrixKernels[14] = sgemm_m16_n9_k20_ldA64_ldB32_ldC16_beta0_pfsigonly;
m_hardwareFlops[14] = 5760;
m_nonZeroFlops[15] = 480;
m_matrixKernels[15] = sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly;
m_hardwareFlops[15] = 2592;
m_nonZeroFlops[16] = 126;
m_matrixKernels[16] = sgemm_m16_n9_k10_ldA64_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[16] = 2880;
m_nonZeroFlops[17] = 306;
m_matrixKernels[17] = sgemm_m16_n9_k10_ldA64_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[17] = 2880;
m_nonZeroFlops[18] = 396;
m_matrixKernels[18] = sgemm_m16_n9_k10_ldA64_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[18] = 2880;
m_nonZeroFlops[19] = 192;
m_matrixKernels[19] = sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly;
m_hardwareFlops[19] = 2592;
m_nonZeroFlops[20] = 18;
m_matrixKernels[20] = sgemm_m16_n9_k4_ldA64_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[20] = 1152;
m_nonZeroFlops[21] = 36;
m_matrixKernels[21] = sgemm_m16_n9_k4_ldA64_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[21] = 1152;
m_nonZeroFlops[22] = 54;
m_matrixKernels[22] = sgemm_m16_n9_k4_ldA64_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[22] = 1152;
m_nonZeroFlops[23] = 48;
m_matrixKernels[23] = sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly;
m_hardwareFlops[23] = 2592;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 12348;
m_matrixKernels[0] = sgemm_m96_n9_k56_ldA96_ldB96_ldC96_beta0_pfsigonly;
m_hardwareFlops[0] = 96768;
m_nonZeroFlops[1] = 27972;
m_matrixKernels[1] = sgemm_m96_n9_k56_ldA96_ldB96_ldC96_beta0_pfsigonly;
m_hardwareFlops[1] = 96768;
m_nonZeroFlops[2] = 30240;
m_matrixKernels[2] = sgemm_m96_n9_k56_ldA96_ldB96_ldC96_beta0_pfsigonly;
m_hardwareFlops[2] = 96768;
m_nonZeroFlops[3] = 4032;
m_matrixKernels[3] = sgemm_m96_n9_k9_ldA96_ldB9_ldC96_beta1_pfsigonly;
m_hardwareFlops[3] = 15552;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 6048;
m_hardwareFlops[0] = 145152;
m_matrixKernels[0] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
m_nonZeroFlops[1] = 19656;
m_hardwareFlops[1] = 145152;
m_matrixKernels[1] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
m_nonZeroFlops[2] = 61992;
m_hardwareFlops[2] = 145152;
m_matrixKernels[2] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
m_nonZeroFlops[3] = 61992;
m_hardwareFlops[3] = 145152;
m_matrixKernels[3] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
m_nonZeroFlops[4] = 19332;
m_hardwareFlops[4] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[4] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 19332;
m_hardwareFlops[5] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[5] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 6048;
m_hardwareFlops[6] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[6] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 24192;
m_hardwareFlops[7] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[7] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 77328;
m_hardwareFlops[8] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[8] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 77328;
m_hardwareFlops[9] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[9] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 101268;
m_hardwareFlops[10] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[10] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 89802;
m_hardwareFlops[11] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[11] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 89190;
m_hardwareFlops[12] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[12] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 89190;
m_hardwareFlops[13] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[13] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 89802;
m_hardwareFlops[14] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[14] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 101268;
m_hardwareFlops[15] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[15] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 24192;
m_hardwareFlops[16] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[16] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 77328;
m_hardwareFlops[17] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[17] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 77328;
m_hardwareFlops[18] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[18] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 109620;
m_hardwareFlops[19] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[19] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 19656;
m_hardwareFlops[20] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[20] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 109620;
m_hardwareFlops[21] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[21] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 69210;
m_hardwareFlops[22] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[22] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 109242;
m_hardwareFlops[23] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[23] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 101844;
m_hardwareFlops[24] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[24] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 109242;
m_hardwareFlops[25] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[25] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 69210;
m_hardwareFlops[26] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[26] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 101844;
m_hardwareFlops[27] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[27] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 101268;
m_hardwareFlops[28] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[28] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 89802;
m_hardwareFlops[29] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[29] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 89190;
m_hardwareFlops[30] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[30] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 69210;
m_hardwareFlops[31] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[31] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 109242;
m_hardwareFlops[32] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[32] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 101844;
m_hardwareFlops[33] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[33] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 78156;
m_hardwareFlops[34] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[34] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 126396;
m_hardwareFlops[35] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[35] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 85788;
m_hardwareFlops[36] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[36] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 61992;
m_hardwareFlops[37] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[37] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 118008;
m_hardwareFlops[38] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[38] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 118008;
m_hardwareFlops[39] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[39] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 89190;
m_hardwareFlops[40] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[40] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 89802;
m_hardwareFlops[41] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[41] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 101268;
m_hardwareFlops[42] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[42] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 109242;
m_hardwareFlops[43] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[43] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 69210;
m_hardwareFlops[44] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[44] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 101844;
m_hardwareFlops[45] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[45] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 61992;
m_hardwareFlops[46] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[46] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 118008;
m_hardwareFlops[47] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[47] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 118008;
m_hardwareFlops[48] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[48] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 78156;
m_hardwareFlops[49] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[49] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 85788;
m_hardwareFlops[50] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[50] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 126396;
m_hardwareFlops[51] = 145152;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[51] = sgemm_m96_n9_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 13608;
m_hardwareFlops[52] = 15552;
m_matrixKernels[52] = sgemm_m96_n9_k9_ldA96_ldB9_ldC96_beta1_pfsigonly;
m_nonZeroFlops[53] = 13608;
m_hardwareFlops[53] = 15552;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = sgemm_m96_n9_k9_ldA96_ldB9_ldC96_beta1_pfsigonly;
#else
m_matrixKernels[53] = sgemm_m96_n9_k9_ldA96_ldB9_ldC96_beta1_pfsigonly;
#endif
m_nonZeroFlops[54] = 13608;
m_hardwareFlops[54] = 15552;
#ifdef ENABLE_STREAM_MATRIX_PREFETCH
m_matrixKernels[54] = sgemm_m96_n9_k9_ldA96_ldB9_ldC96_beta1_pfsigonly;
#else
m_matrixKernels[54] = sgemm_m96_n9_k9_ldA96_ldB9_ldC96_beta1_pfsigonly;
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
m_sparseSwitch[59] = -1; 
#endif

#define STAR_NNZ 81

#endif

#if CONVERGENCE_ORDER==8

#ifdef TIME_KERNEL
m_nonZeroFlops[0] = 26028;
m_matrixKernels[0] = sgemm_m96_n9_k120_ldA96_ldB128_ldC96_beta0_pfsigonly;
m_hardwareFlops[0] = 207360;
m_nonZeroFlops[1] = 58104;
m_matrixKernels[1] = sgemm_m96_n9_k120_ldA96_ldB128_ldC96_beta0_pfsigonly;
m_hardwareFlops[1] = 207360;
m_nonZeroFlops[2] = 61884;
m_matrixKernels[2] = sgemm_m96_n9_k120_ldA96_ldB128_ldC96_beta0_pfsigonly;
m_hardwareFlops[2] = 207360;
m_nonZeroFlops[3] = 4032;
m_matrixKernels[3] = sgemm_m96_n9_k9_ldA96_ldB9_ldC96_beta1_pfsigonly;
m_hardwareFlops[3] = 15552;
m_nonZeroFlops[4] = 12348;
m_matrixKernels[4] = sgemm_m64_n9_k84_ldA96_ldB96_ldC64_beta0_pfsigonly;
m_hardwareFlops[4] = 96768;
m_nonZeroFlops[5] = 27972;
m_matrixKernels[5] = sgemm_m64_n9_k84_ldA96_ldB96_ldC64_beta0_pfsigonly;
m_hardwareFlops[5] = 96768;
m_nonZeroFlops[6] = 30240;
m_matrixKernels[6] = sgemm_m64_n9_k84_ldA96_ldB96_ldC64_beta0_pfsigonly;
m_hardwareFlops[6] = 96768;
m_nonZeroFlops[7] = 2688;
m_matrixKernels[7] = sgemm_m64_n9_k9_ldA64_ldB9_ldC64_beta1_pfsigonly;
m_hardwareFlops[7] = 10368;
m_nonZeroFlops[8] = 5292;
m_matrixKernels[8] = sgemm_m48_n9_k56_ldA96_ldB64_ldC48_beta0_pfsigonly;
m_hardwareFlops[8] = 48384;
m_nonZeroFlops[9] = 12096;
m_matrixKernels[9] = sgemm_m48_n9_k56_ldA96_ldB64_ldC48_beta0_pfsigonly;
m_hardwareFlops[9] = 48384;
m_nonZeroFlops[10] = 13356;
m_matrixKernels[10] = sgemm_m48_n9_k56_ldA96_ldB64_ldC48_beta0_pfsigonly;
m_hardwareFlops[10] = 48384;
m_nonZeroFlops[11] = 1680;
m_matrixKernels[11] = sgemm_m48_n9_k9_ldA48_ldB9_ldC48_beta1_pfsigonly;
m_hardwareFlops[11] = 7776;
m_nonZeroFlops[12] = 1944;
m_matrixKernels[12] = sgemm_m32_n9_k35_ldA96_ldB48_ldC32_beta0_pfsigonly;
m_hardwareFlops[12] = 20160;
m_nonZeroFlops[13] = 4536;
m_matrixKernels[13] = sgemm_m32_n9_k35_ldA96_ldB48_ldC32_beta0_pfsigonly;
m_hardwareFlops[13] = 20160;
m_nonZeroFlops[14] = 5166;
m_matrixKernels[14] = sgemm_m32_n9_k35_ldA96_ldB48_ldC32_beta0_pfsigonly;
m_hardwareFlops[14] = 20160;
m_nonZeroFlops[15] = 960;
m_matrixKernels[15] = sgemm_m32_n9_k9_ldA32_ldB9_ldC32_beta1_pfsigonly;
m_hardwareFlops[15] = 5184;
m_nonZeroFlops[16] = 594;
m_matrixKernels[16] = sgemm_m16_n9_k20_ldA96_ldB32_ldC16_beta0_pfsigonly;
m_hardwareFlops[16] = 5760;
m_nonZeroFlops[17] = 1386;
m_matrixKernels[17] = sgemm_m16_n9_k20_ldA96_ldB32_ldC16_beta0_pfsigonly;
m_hardwareFlops[17] = 5760;
m_nonZeroFlops[18] = 1656;
m_matrixKernels[18] = sgemm_m16_n9_k20_ldA96_ldB32_ldC16_beta0_pfsigonly;
m_hardwareFlops[18] = 5760;
m_nonZeroFlops[19] = 480;
m_matrixKernels[19] = sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly;
m_hardwareFlops[19] = 2592;
m_nonZeroFlops[20] = 126;
m_matrixKernels[20] = sgemm_m16_n9_k10_ldA96_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[20] = 2880;
m_nonZeroFlops[21] = 306;
m_matrixKernels[21] = sgemm_m16_n9_k10_ldA96_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[21] = 2880;
m_nonZeroFlops[22] = 396;
m_matrixKernels[22] = sgemm_m16_n9_k10_ldA96_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[22] = 2880;
m_nonZeroFlops[23] = 192;
m_matrixKernels[23] = sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly;
m_hardwareFlops[23] = 2592;
m_nonZeroFlops[24] = 18;
m_matrixKernels[24] = sgemm_m16_n9_k4_ldA96_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[24] = 1152;
m_nonZeroFlops[25] = 36;
m_matrixKernels[25] = sgemm_m16_n9_k4_ldA96_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[25] = 1152;
m_nonZeroFlops[26] = 54;
m_matrixKernels[26] = sgemm_m16_n9_k4_ldA96_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[26] = 1152;
m_nonZeroFlops[27] = 48;
m_matrixKernels[27] = sgemm_m16_n9_k9_ldA16_ldB9_ldC16_beta1_pfsigonly;
m_hardwareFlops[27] = 2592;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 26028;
m_matrixKernels[0] = sgemm_m128_n9_k84_ldA128_ldB128_ldC128_beta0_pfsigonly;
m_hardwareFlops[0] = 193536;
m_nonZeroFlops[1] = 58104;
m_matrixKernels[1] = sgemm_m128_n9_k84_ldA128_ldB128_ldC128_beta0_pfsigonly;
m_hardwareFlops[1] = 193536;
m_nonZeroFlops[2] = 61884;
m_matrixKernels[2] = sgemm_m128_n9_k84_ldA128_ldB128_ldC128_beta0_pfsigonly;
m_hardwareFlops[2] = 193536;
m_nonZeroFlops[3] = 5760;
m_matrixKernels[3] = sgemm_m128_n9_k9_ldA128_ldB9_ldC128_beta1_pfsigonly;
m_hardwareFlops[3] = 20736;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 9720;
m_hardwareFlops[0] = 276480;
m_matrixKernels[0] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
m_nonZeroFlops[1] = 35424;
m_hardwareFlops[1] = 276480;
m_matrixKernels[1] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
m_nonZeroFlops[2] = 125928;
m_hardwareFlops[2] = 276480;
m_matrixKernels[2] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
m_nonZeroFlops[3] = 125928;
m_hardwareFlops[3] = 276480;
m_matrixKernels[3] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
m_nonZeroFlops[4] = 34848;
m_hardwareFlops[4] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[4] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 34848;
m_hardwareFlops[5] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[5] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 9720;
m_hardwareFlops[6] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[6] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 43740;
m_hardwareFlops[7] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[7] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 156816;
m_hardwareFlops[8] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[8] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 156816;
m_hardwareFlops[9] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[9] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 207162;
m_hardwareFlops[10] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[10] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 181980;
m_hardwareFlops[11] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[11] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 181620;
m_hardwareFlops[12] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[12] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 181620;
m_hardwareFlops[13] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[13] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 181980;
m_hardwareFlops[14] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[14] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 207162;
m_hardwareFlops[15] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[15] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 43740;
m_hardwareFlops[16] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[16] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 156816;
m_hardwareFlops[17] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[17] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 156816;
m_hardwareFlops[18] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[18] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 224676;
m_hardwareFlops[19] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[19] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 35424;
m_hardwareFlops[20] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[20] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 224676;
m_hardwareFlops[21] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[21] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 140760;
m_hardwareFlops[22] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[22] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 224388;
m_hardwareFlops[23] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[23] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 207144;
m_hardwareFlops[24] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[24] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 224388;
m_hardwareFlops[25] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[25] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 140760;
m_hardwareFlops[26] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[26] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 207144;
m_hardwareFlops[27] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[27] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 207162;
m_hardwareFlops[28] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[28] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 181980;
m_hardwareFlops[29] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[29] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 181620;
m_hardwareFlops[30] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[30] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 140760;
m_hardwareFlops[31] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[31] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 224388;
m_hardwareFlops[32] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[32] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 207144;
m_hardwareFlops[33] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[33] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 158688;
m_hardwareFlops[34] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[34] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 258444;
m_hardwareFlops[35] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[35] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 175284;
m_hardwareFlops[36] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[36] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 125928;
m_hardwareFlops[37] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[37] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 241578;
m_hardwareFlops[38] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[38] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 241578;
m_hardwareFlops[39] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[39] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 181620;
m_hardwareFlops[40] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[40] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 181980;
m_hardwareFlops[41] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[41] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 207162;
m_hardwareFlops[42] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[42] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 224388;
m_hardwareFlops[43] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[43] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 140760;
m_hardwareFlops[44] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[44] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 207144;
m_hardwareFlops[45] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[45] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 125928;
m_hardwareFlops[46] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[46] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 241578;
m_hardwareFlops[47] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[47] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 241578;
m_hardwareFlops[48] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[48] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 158688;
m_hardwareFlops[49] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[49] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 175284;
m_hardwareFlops[50] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[50] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 258444;
m_hardwareFlops[51] = 276480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[51] = sgemm_m128_n9_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 19440;
m_hardwareFlops[52] = 20736;
m_matrixKernels[52] = sgemm_m128_n9_k9_ldA128_ldB9_ldC128_beta1_pfsigonly;
m_nonZeroFlops[53] = 19440;
m_hardwareFlops[53] = 20736;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = sgemm_m128_n9_k9_ldA128_ldB9_ldC128_beta1_pfsigonly;
#else
m_matrixKernels[53] = sgemm_m128_n9_k9_ldA128_ldB9_ldC128_beta1_pfsigonly;
#endif
m_nonZeroFlops[54] = 19440;
m_hardwareFlops[54] = 20736;
#ifdef ENABLE_STREAM_MATRIX_PREFETCH
m_matrixKernels[54] = sgemm_m128_n9_k9_ldA128_ldB9_ldC128_beta1_pfsigonly;
#else
m_matrixKernels[54] = sgemm_m128_n9_k9_ldA128_ldB9_ldC128_beta1_pfsigonly;
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
m_sparseSwitch[59] = -1; 
#endif

#define STAR_NNZ 81

#endif

