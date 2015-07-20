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
// @date 2015-07-13 15:29:43.565018
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
m_nonZeroFlops[0] = 54;
m_matrixKernels[0] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[0] = 3456;
m_nonZeroFlops[1] = 108;
m_matrixKernels[1] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[1] = 3456;
m_nonZeroFlops[2] = 162;
m_matrixKernels[2] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[2] = 3456;
m_nonZeroFlops[3] = 48;
m_matrixKernels[3] = sgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
m_hardwareFlops[3] = 23328;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 54;
m_matrixKernels[0] = sgemm_m16_n27_k1_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[0] = 864;
m_nonZeroFlops[1] = 108;
m_matrixKernels[1] = sgemm_m16_n27_k1_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[1] = 864;
m_nonZeroFlops[2] = 162;
m_matrixKernels[2] = sgemm_m16_n27_k1_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[2] = 864;
m_nonZeroFlops[3] = 192;
m_matrixKernels[3] = sgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
m_hardwareFlops[3] = 23328;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 324;
m_hardwareFlops[0] = 3456;
m_matrixKernels[0] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_nonZeroFlops[1] = 432;
m_hardwareFlops[1] = 3456;
m_matrixKernels[1] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_nonZeroFlops[2] = 540;
m_hardwareFlops[2] = 3456;
m_matrixKernels[2] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_nonZeroFlops[3] = 540;
m_hardwareFlops[3] = 3456;
m_matrixKernels[3] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_nonZeroFlops[4] = 432;
m_hardwareFlops[4] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[4] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 432;
m_hardwareFlops[5] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[5] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 324;
m_hardwareFlops[6] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[6] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 486;
m_hardwareFlops[7] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[7] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 648;
m_hardwareFlops[8] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[8] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 648;
m_hardwareFlops[9] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[9] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 756;
m_hardwareFlops[10] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[10] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 702;
m_hardwareFlops[11] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[11] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 702;
m_hardwareFlops[12] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[12] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 702;
m_hardwareFlops[13] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[13] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 702;
m_hardwareFlops[14] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[14] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 756;
m_hardwareFlops[15] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[15] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 486;
m_hardwareFlops[16] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[16] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 648;
m_hardwareFlops[17] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[17] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 648;
m_hardwareFlops[18] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[18] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 756;
m_hardwareFlops[19] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[19] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 432;
m_hardwareFlops[20] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[20] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 756;
m_hardwareFlops[21] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[21] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 594;
m_hardwareFlops[22] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[22] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 756;
m_hardwareFlops[23] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[23] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 756;
m_hardwareFlops[24] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[24] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 756;
m_hardwareFlops[25] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[25] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 594;
m_hardwareFlops[26] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[26] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 756;
m_hardwareFlops[27] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[27] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 756;
m_hardwareFlops[28] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[28] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 702;
m_hardwareFlops[29] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[29] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 702;
m_hardwareFlops[30] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[30] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 594;
m_hardwareFlops[31] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[31] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 756;
m_hardwareFlops[32] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[32] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 756;
m_hardwareFlops[33] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[33] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 648;
m_hardwareFlops[34] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[34] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 864;
m_hardwareFlops[35] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[35] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 648;
m_hardwareFlops[36] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[36] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 540;
m_hardwareFlops[37] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[37] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 810;
m_hardwareFlops[38] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[38] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 810;
m_hardwareFlops[39] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[39] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 702;
m_hardwareFlops[40] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[40] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 702;
m_hardwareFlops[41] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[41] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 756;
m_hardwareFlops[42] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[42] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 756;
m_hardwareFlops[43] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[43] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 594;
m_hardwareFlops[44] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[44] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 756;
m_hardwareFlops[45] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[45] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 540;
m_hardwareFlops[46] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[46] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 810;
m_hardwareFlops[47] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[47] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 810;
m_hardwareFlops[48] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[48] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 648;
m_hardwareFlops[49] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[49] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 648;
m_hardwareFlops[50] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[50] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 864;
m_hardwareFlops[51] = 3456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[51] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 648;
m_hardwareFlops[52] = 23328;
m_matrixKernels[52] = sgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
m_nonZeroFlops[53] = 648;
m_hardwareFlops[53] = 23328;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = sgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
#else
m_matrixKernels[53] = sgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
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

#define STAR_NNZ 729

#endif

#if CONVERGENCE_ORDER==3

#ifdef TIME_KERNEL
m_nonZeroFlops[0] = 378;
m_matrixKernels[0] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[0] = 8640;
m_nonZeroFlops[1] = 918;
m_matrixKernels[1] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[1] = 8640;
m_nonZeroFlops[2] = 1188;
m_matrixKernels[2] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[2] = 8640;
m_nonZeroFlops[3] = 192;
m_matrixKernels[3] = sgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
m_hardwareFlops[3] = 23328;
m_nonZeroFlops[4] = 54;
m_matrixKernels[4] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[4] = 3456;
m_nonZeroFlops[5] = 108;
m_matrixKernels[5] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[5] = 3456;
m_nonZeroFlops[6] = 162;
m_matrixKernels[6] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[6] = 3456;
m_nonZeroFlops[7] = 48;
m_matrixKernels[7] = sgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
m_hardwareFlops[7] = 23328;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 378;
m_matrixKernels[0] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[0] = 3456;
m_nonZeroFlops[1] = 918;
m_matrixKernels[1] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[1] = 3456;
m_nonZeroFlops[2] = 1188;
m_matrixKernels[2] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[2] = 3456;
m_nonZeroFlops[3] = 480;
m_matrixKernels[3] = sgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
m_hardwareFlops[3] = 23328;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 1080;
m_hardwareFlops[0] = 8640;
m_matrixKernels[0] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_nonZeroFlops[1] = 1836;
m_hardwareFlops[1] = 8640;
m_matrixKernels[1] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_nonZeroFlops[2] = 2916;
m_hardwareFlops[2] = 8640;
m_matrixKernels[2] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_nonZeroFlops[3] = 2916;
m_hardwareFlops[3] = 8640;
m_matrixKernels[3] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_nonZeroFlops[4] = 1836;
m_hardwareFlops[4] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[4] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 1836;
m_hardwareFlops[5] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[5] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 1080;
m_hardwareFlops[6] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[6] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 2160;
m_hardwareFlops[7] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[7] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 3672;
m_hardwareFlops[8] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[8] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 3672;
m_hardwareFlops[9] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[9] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 4374;
m_hardwareFlops[10] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[10] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 4104;
m_hardwareFlops[11] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[11] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 4104;
m_hardwareFlops[12] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[12] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 4104;
m_hardwareFlops[13] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[13] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 4104;
m_hardwareFlops[14] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[14] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 4374;
m_hardwareFlops[15] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[15] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 2160;
m_hardwareFlops[16] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[16] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 3672;
m_hardwareFlops[17] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[17] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 3672;
m_hardwareFlops[18] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[18] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 4536;
m_hardwareFlops[19] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[19] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 1836;
m_hardwareFlops[20] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[20] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 4536;
m_hardwareFlops[21] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[21] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 3294;
m_hardwareFlops[22] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[22] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 4536;
m_hardwareFlops[23] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[23] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 4536;
m_hardwareFlops[24] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[24] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 4536;
m_hardwareFlops[25] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[25] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 3294;
m_hardwareFlops[26] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[26] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 4536;
m_hardwareFlops[27] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[27] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 4374;
m_hardwareFlops[28] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[28] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 4104;
m_hardwareFlops[29] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[29] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 4104;
m_hardwareFlops[30] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[30] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 3294;
m_hardwareFlops[31] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[31] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 4536;
m_hardwareFlops[32] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[32] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 4536;
m_hardwareFlops[33] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[33] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 3672;
m_hardwareFlops[34] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[34] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 5292;
m_hardwareFlops[35] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[35] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 3780;
m_hardwareFlops[36] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[36] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 2916;
m_hardwareFlops[37] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[37] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 5022;
m_hardwareFlops[38] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[38] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 5022;
m_hardwareFlops[39] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[39] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 4104;
m_hardwareFlops[40] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[40] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 4104;
m_hardwareFlops[41] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[41] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 4374;
m_hardwareFlops[42] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[42] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 4536;
m_hardwareFlops[43] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[43] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 3294;
m_hardwareFlops[44] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[44] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 4536;
m_hardwareFlops[45] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[45] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 2916;
m_hardwareFlops[46] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[46] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 5022;
m_hardwareFlops[47] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[47] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 5022;
m_hardwareFlops[48] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[48] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 3672;
m_hardwareFlops[49] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[49] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 3780;
m_hardwareFlops[50] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[50] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 5292;
m_hardwareFlops[51] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_BL2viaC;
#else
m_matrixKernels[51] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 1620;
m_hardwareFlops[52] = 23328;
m_matrixKernels[52] = sgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
m_nonZeroFlops[53] = 1620;
m_hardwareFlops[53] = 23328;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = sgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
#else
m_matrixKernels[53] = sgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
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

#define STAR_NNZ 729

#endif

#if CONVERGENCE_ORDER==4

#ifdef TIME_KERNEL
m_nonZeroFlops[0] = 1782;
m_matrixKernels[0] = sgemm_m16_n27_k20_ldA16_ldB32_ldC16_beta0_pfsigonly;
m_hardwareFlops[0] = 17280;
m_nonZeroFlops[1] = 4158;
m_matrixKernels[1] = sgemm_m16_n27_k20_ldA16_ldB32_ldC16_beta0_pfsigonly;
m_hardwareFlops[1] = 17280;
m_nonZeroFlops[2] = 4968;
m_matrixKernels[2] = sgemm_m16_n27_k20_ldA16_ldB32_ldC16_beta0_pfsigonly;
m_hardwareFlops[2] = 17280;
m_nonZeroFlops[3] = 480;
m_matrixKernels[3] = sgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
m_hardwareFlops[3] = 23328;
m_nonZeroFlops[4] = 378;
m_matrixKernels[4] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[4] = 8640;
m_nonZeroFlops[5] = 918;
m_matrixKernels[5] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[5] = 8640;
m_nonZeroFlops[6] = 1188;
m_matrixKernels[6] = sgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[6] = 8640;
m_nonZeroFlops[7] = 192;
m_matrixKernels[7] = sgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
m_hardwareFlops[7] = 23328;
m_nonZeroFlops[8] = 54;
m_matrixKernels[8] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[8] = 3456;
m_nonZeroFlops[9] = 108;
m_matrixKernels[9] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[9] = 3456;
m_nonZeroFlops[10] = 162;
m_matrixKernels[10] = sgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[10] = 3456;
m_nonZeroFlops[11] = 48;
m_matrixKernels[11] = sgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
m_hardwareFlops[11] = 23328;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 1782;
m_matrixKernels[0] = sgemm_m32_n27_k10_ldA32_ldB32_ldC32_beta0_pfsigonly;
m_hardwareFlops[0] = 17280;
m_nonZeroFlops[1] = 4158;
m_matrixKernels[1] = sgemm_m32_n27_k10_ldA32_ldB32_ldC32_beta0_pfsigonly;
m_hardwareFlops[1] = 17280;
m_nonZeroFlops[2] = 4968;
m_matrixKernels[2] = sgemm_m32_n27_k10_ldA32_ldB32_ldC32_beta0_pfsigonly;
m_hardwareFlops[2] = 17280;
m_nonZeroFlops[3] = 960;
m_matrixKernels[3] = sgemm_m32_n27_k27_ldA32_ldB27_ldC32_beta1_pfsigonly;
m_hardwareFlops[3] = 46656;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 2700;
m_hardwareFlops[0] = 34560;
m_matrixKernels[0] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
m_nonZeroFlops[1] = 5616;
m_hardwareFlops[1] = 34560;
m_matrixKernels[1] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
m_nonZeroFlops[2] = 11016;
m_hardwareFlops[2] = 34560;
m_matrixKernels[2] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
m_nonZeroFlops[3] = 11016;
m_hardwareFlops[3] = 34560;
m_matrixKernels[3] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
m_nonZeroFlops[4] = 5616;
m_hardwareFlops[4] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[4] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 5616;
m_hardwareFlops[5] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[5] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 2700;
m_hardwareFlops[6] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[6] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 6750;
m_hardwareFlops[7] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[7] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 14040;
m_hardwareFlops[8] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[8] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 14040;
m_hardwareFlops[9] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[9] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 17280;
m_hardwareFlops[10] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[10] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 15444;
m_hardwareFlops[11] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[11] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 15930;
m_hardwareFlops[12] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[12] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 15930;
m_hardwareFlops[13] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[13] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 15444;
m_hardwareFlops[14] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[14] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 17280;
m_hardwareFlops[15] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[15] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 6750;
m_hardwareFlops[16] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[16] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 14040;
m_hardwareFlops[17] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[17] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 14040;
m_hardwareFlops[18] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[18] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 18252;
m_hardwareFlops[19] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[19] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 5616;
m_hardwareFlops[20] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[20] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 18252;
m_hardwareFlops[21] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[21] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 12528;
m_hardwareFlops[22] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[22] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 18360;
m_hardwareFlops[23] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[23] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 17550;
m_hardwareFlops[24] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[24] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 18360;
m_hardwareFlops[25] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[25] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 12528;
m_hardwareFlops[26] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[26] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 17550;
m_hardwareFlops[27] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[27] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 17280;
m_hardwareFlops[28] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[28] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 15444;
m_hardwareFlops[29] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[29] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 15930;
m_hardwareFlops[30] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[30] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 12528;
m_hardwareFlops[31] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[31] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 18360;
m_hardwareFlops[32] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[32] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 17550;
m_hardwareFlops[33] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[33] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 13932;
m_hardwareFlops[34] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[34] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 21276;
m_hardwareFlops[35] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[35] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 14796;
m_hardwareFlops[36] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[36] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 11016;
m_hardwareFlops[37] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[37] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 20034;
m_hardwareFlops[38] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[38] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 20034;
m_hardwareFlops[39] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[39] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 15930;
m_hardwareFlops[40] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[40] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 15444;
m_hardwareFlops[41] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[41] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 17280;
m_hardwareFlops[42] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[42] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 18360;
m_hardwareFlops[43] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[43] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 12528;
m_hardwareFlops[44] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[44] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 17550;
m_hardwareFlops[45] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[45] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 11016;
m_hardwareFlops[46] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[46] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 20034;
m_hardwareFlops[47] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[47] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 20034;
m_hardwareFlops[48] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[48] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 13932;
m_hardwareFlops[49] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[49] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 14796;
m_hardwareFlops[50] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[50] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 21276;
m_hardwareFlops[51] = 34560;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_BL2viaC;
#else
m_matrixKernels[51] = sgemm_m32_n27_k20_ldA32_ldB32_ldC32_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 3240;
m_hardwareFlops[52] = 46656;
m_matrixKernels[52] = sgemm_m32_n27_k27_ldA32_ldB27_ldC32_beta1_pfsigonly;
m_nonZeroFlops[53] = 3240;
m_hardwareFlops[53] = 46656;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = sgemm_m32_n27_k27_ldA32_ldB27_ldC32_beta1_pfsigonly;
#else
m_matrixKernels[53] = sgemm_m32_n27_k27_ldA32_ldB27_ldC32_beta1_pfsigonly;
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

#define STAR_NNZ 729

#endif

#if CONVERGENCE_ORDER==5

#ifdef TIME_KERNEL
m_nonZeroFlops[0] = 5832;
m_matrixKernels[0] = sgemm_m32_n27_k35_ldA32_ldB48_ldC32_beta0_pfsigonly;
m_hardwareFlops[0] = 60480;
m_nonZeroFlops[1] = 13608;
m_matrixKernels[1] = sgemm_m32_n27_k35_ldA32_ldB48_ldC32_beta0_pfsigonly;
m_hardwareFlops[1] = 60480;
m_nonZeroFlops[2] = 15498;
m_matrixKernels[2] = sgemm_m32_n27_k35_ldA32_ldB48_ldC32_beta0_pfsigonly;
m_hardwareFlops[2] = 60480;
m_nonZeroFlops[3] = 960;
m_matrixKernels[3] = sgemm_m32_n27_k27_ldA32_ldB27_ldC32_beta1_pfsigonly;
m_hardwareFlops[3] = 46656;
m_nonZeroFlops[4] = 1782;
m_matrixKernels[4] = sgemm_m16_n27_k20_ldA32_ldB32_ldC16_beta0_pfsigonly;
m_hardwareFlops[4] = 17280;
m_nonZeroFlops[5] = 4158;
m_matrixKernels[5] = sgemm_m16_n27_k20_ldA32_ldB32_ldC16_beta0_pfsigonly;
m_hardwareFlops[5] = 17280;
m_nonZeroFlops[6] = 4968;
m_matrixKernels[6] = sgemm_m16_n27_k20_ldA32_ldB32_ldC16_beta0_pfsigonly;
m_hardwareFlops[6] = 17280;
m_nonZeroFlops[7] = 480;
m_matrixKernels[7] = sgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
m_hardwareFlops[7] = 23328;
m_nonZeroFlops[8] = 378;
m_matrixKernels[8] = sgemm_m16_n27_k10_ldA32_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[8] = 8640;
m_nonZeroFlops[9] = 918;
m_matrixKernels[9] = sgemm_m16_n27_k10_ldA32_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[9] = 8640;
m_nonZeroFlops[10] = 1188;
m_matrixKernels[10] = sgemm_m16_n27_k10_ldA32_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[10] = 8640;
m_nonZeroFlops[11] = 192;
m_matrixKernels[11] = sgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
m_hardwareFlops[11] = 23328;
m_nonZeroFlops[12] = 54;
m_matrixKernels[12] = sgemm_m16_n27_k4_ldA32_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[12] = 3456;
m_nonZeroFlops[13] = 108;
m_matrixKernels[13] = sgemm_m16_n27_k4_ldA32_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[13] = 3456;
m_nonZeroFlops[14] = 162;
m_matrixKernels[14] = sgemm_m16_n27_k4_ldA32_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[14] = 3456;
m_nonZeroFlops[15] = 48;
m_matrixKernels[15] = sgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
m_hardwareFlops[15] = 23328;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 5832;
m_matrixKernels[0] = sgemm_m48_n27_k20_ldA48_ldB48_ldC48_beta0_pfsigonly;
m_hardwareFlops[0] = 51840;
m_nonZeroFlops[1] = 13608;
m_matrixKernels[1] = sgemm_m48_n27_k20_ldA48_ldB48_ldC48_beta0_pfsigonly;
m_hardwareFlops[1] = 51840;
m_nonZeroFlops[2] = 15498;
m_matrixKernels[2] = sgemm_m48_n27_k20_ldA48_ldB48_ldC48_beta0_pfsigonly;
m_hardwareFlops[2] = 51840;
m_nonZeroFlops[3] = 1680;
m_matrixKernels[3] = sgemm_m48_n27_k27_ldA48_ldB27_ldC48_beta1_pfsigonly;
m_hardwareFlops[3] = 69984;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 5670;
m_hardwareFlops[0] = 90720;
m_matrixKernels[0] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
m_nonZeroFlops[1] = 13986;
m_hardwareFlops[1] = 90720;
m_matrixKernels[1] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
m_nonZeroFlops[2] = 32886;
m_hardwareFlops[2] = 90720;
m_matrixKernels[2] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
m_nonZeroFlops[3] = 32886;
m_hardwareFlops[3] = 90720;
m_matrixKernels[3] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
m_nonZeroFlops[4] = 13878;
m_hardwareFlops[4] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[4] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 13878;
m_hardwareFlops[5] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[5] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 5670;
m_hardwareFlops[6] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[6] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 17010;
m_hardwareFlops[7] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[7] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 41634;
m_hardwareFlops[8] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[8] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 41634;
m_hardwareFlops[9] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[9] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 52866;
m_hardwareFlops[10] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[10] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 46440;
m_hardwareFlops[11] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[11] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 47304;
m_hardwareFlops[12] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[12] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 47304;
m_hardwareFlops[13] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[13] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 46440;
m_hardwareFlops[14] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[14] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 52866;
m_hardwareFlops[15] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[15] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 17010;
m_hardwareFlops[16] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[16] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 41634;
m_hardwareFlops[17] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[17] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 41634;
m_hardwareFlops[18] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[18] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 56430;
m_hardwareFlops[19] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[19] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 13986;
m_hardwareFlops[20] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[20] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 56430;
m_hardwareFlops[21] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[21] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 36936;
m_hardwareFlops[22] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[22] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 56106;
m_hardwareFlops[23] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[23] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 52920;
m_hardwareFlops[24] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[24] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 56106;
m_hardwareFlops[25] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[25] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 36936;
m_hardwareFlops[26] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[26] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 52920;
m_hardwareFlops[27] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[27] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 52866;
m_hardwareFlops[28] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[28] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 46440;
m_hardwareFlops[29] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[29] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 47304;
m_hardwareFlops[30] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[30] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 36936;
m_hardwareFlops[31] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[31] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 56106;
m_hardwareFlops[32] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[32] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 52920;
m_hardwareFlops[33] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[33] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 41742;
m_hardwareFlops[34] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[34] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 65502;
m_hardwareFlops[35] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[35] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 44874;
m_hardwareFlops[36] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[36] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 32886;
m_hardwareFlops[37] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[37] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 61074;
m_hardwareFlops[38] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[38] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 61074;
m_hardwareFlops[39] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[39] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 47304;
m_hardwareFlops[40] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[40] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 46440;
m_hardwareFlops[41] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[41] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 52866;
m_hardwareFlops[42] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[42] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 56106;
m_hardwareFlops[43] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[43] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 36936;
m_hardwareFlops[44] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[44] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 52920;
m_hardwareFlops[45] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[45] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 32886;
m_hardwareFlops[46] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[46] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 61074;
m_hardwareFlops[47] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[47] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 61074;
m_hardwareFlops[48] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[48] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 41742;
m_hardwareFlops[49] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[49] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 44874;
m_hardwareFlops[50] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[50] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 65502;
m_hardwareFlops[51] = 90720;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_BL2viaC;
#else
m_matrixKernels[51] = sgemm_m48_n27_k35_ldA48_ldB48_ldC48_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 5670;
m_hardwareFlops[52] = 69984;
m_matrixKernels[52] = sgemm_m48_n27_k27_ldA48_ldB27_ldC48_beta1_pfsigonly;
m_nonZeroFlops[53] = 5670;
m_hardwareFlops[53] = 69984;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = sgemm_m48_n27_k27_ldA48_ldB27_ldC48_beta1_pfsigonly;
#else
m_matrixKernels[53] = sgemm_m48_n27_k27_ldA48_ldB27_ldC48_beta1_pfsigonly;
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

#define STAR_NNZ 729

#endif

#if CONVERGENCE_ORDER==6

#ifdef TIME_KERNEL
m_nonZeroFlops[0] = 15876;
m_matrixKernels[0] = sgemm_m48_n27_k56_ldA48_ldB64_ldC48_beta0_pfsigonly;
m_hardwareFlops[0] = 145152;
m_nonZeroFlops[1] = 36288;
m_matrixKernels[1] = sgemm_m48_n27_k56_ldA48_ldB64_ldC48_beta0_pfsigonly;
m_hardwareFlops[1] = 145152;
m_nonZeroFlops[2] = 40068;
m_matrixKernels[2] = sgemm_m48_n27_k56_ldA48_ldB64_ldC48_beta0_pfsigonly;
m_hardwareFlops[2] = 145152;
m_nonZeroFlops[3] = 1680;
m_matrixKernels[3] = sgemm_m48_n27_k27_ldA48_ldB27_ldC48_beta1_pfsigonly;
m_hardwareFlops[3] = 69984;
m_nonZeroFlops[4] = 5832;
m_matrixKernels[4] = sgemm_m32_n27_k35_ldA48_ldB48_ldC32_beta0_pfsigonly;
m_hardwareFlops[4] = 60480;
m_nonZeroFlops[5] = 13608;
m_matrixKernels[5] = sgemm_m32_n27_k35_ldA48_ldB48_ldC32_beta0_pfsigonly;
m_hardwareFlops[5] = 60480;
m_nonZeroFlops[6] = 15498;
m_matrixKernels[6] = sgemm_m32_n27_k35_ldA48_ldB48_ldC32_beta0_pfsigonly;
m_hardwareFlops[6] = 60480;
m_nonZeroFlops[7] = 960;
m_matrixKernels[7] = sgemm_m32_n27_k27_ldA32_ldB27_ldC32_beta1_pfsigonly;
m_hardwareFlops[7] = 46656;
m_nonZeroFlops[8] = 1782;
m_matrixKernels[8] = sgemm_m16_n27_k20_ldA48_ldB32_ldC16_beta0_pfsigonly;
m_hardwareFlops[8] = 17280;
m_nonZeroFlops[9] = 4158;
m_matrixKernels[9] = sgemm_m16_n27_k20_ldA48_ldB32_ldC16_beta0_pfsigonly;
m_hardwareFlops[9] = 17280;
m_nonZeroFlops[10] = 4968;
m_matrixKernels[10] = sgemm_m16_n27_k20_ldA48_ldB32_ldC16_beta0_pfsigonly;
m_hardwareFlops[10] = 17280;
m_nonZeroFlops[11] = 480;
m_matrixKernels[11] = sgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
m_hardwareFlops[11] = 23328;
m_nonZeroFlops[12] = 378;
m_matrixKernels[12] = sgemm_m16_n27_k10_ldA48_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[12] = 8640;
m_nonZeroFlops[13] = 918;
m_matrixKernels[13] = sgemm_m16_n27_k10_ldA48_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[13] = 8640;
m_nonZeroFlops[14] = 1188;
m_matrixKernels[14] = sgemm_m16_n27_k10_ldA48_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[14] = 8640;
m_nonZeroFlops[15] = 192;
m_matrixKernels[15] = sgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
m_hardwareFlops[15] = 23328;
m_nonZeroFlops[16] = 54;
m_matrixKernels[16] = sgemm_m16_n27_k4_ldA48_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[16] = 3456;
m_nonZeroFlops[17] = 108;
m_matrixKernels[17] = sgemm_m16_n27_k4_ldA48_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[17] = 3456;
m_nonZeroFlops[18] = 162;
m_matrixKernels[18] = sgemm_m16_n27_k4_ldA48_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[18] = 3456;
m_nonZeroFlops[19] = 48;
m_matrixKernels[19] = sgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
m_hardwareFlops[19] = 23328;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 15876;
m_matrixKernels[0] = sgemm_m64_n27_k35_ldA64_ldB64_ldC64_beta0_pfsigonly;
m_hardwareFlops[0] = 120960;
m_nonZeroFlops[1] = 36288;
m_matrixKernels[1] = sgemm_m64_n27_k35_ldA64_ldB64_ldC64_beta0_pfsigonly;
m_hardwareFlops[1] = 120960;
m_nonZeroFlops[2] = 40068;
m_matrixKernels[2] = sgemm_m64_n27_k35_ldA64_ldB64_ldC64_beta0_pfsigonly;
m_hardwareFlops[2] = 120960;
m_nonZeroFlops[3] = 2688;
m_matrixKernels[3] = sgemm_m64_n27_k27_ldA64_ldB27_ldC64_beta1_pfsigonly;
m_hardwareFlops[3] = 93312;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 10584;
m_hardwareFlops[0] = 193536;
m_matrixKernels[0] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
m_nonZeroFlops[1] = 30240;
m_hardwareFlops[1] = 193536;
m_matrixKernels[1] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
m_nonZeroFlops[2] = 83160;
m_hardwareFlops[2] = 193536;
m_matrixKernels[2] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
m_nonZeroFlops[3] = 83160;
m_hardwareFlops[3] = 193536;
m_matrixKernels[3] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
m_nonZeroFlops[4] = 29808;
m_hardwareFlops[4] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[4] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 29808;
m_hardwareFlops[5] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[5] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 10584;
m_hardwareFlops[6] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[6] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 37044;
m_hardwareFlops[7] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[7] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 104328;
m_hardwareFlops[8] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[8] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 104328;
m_hardwareFlops[9] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[9] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 135972;
m_hardwareFlops[10] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[10] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 119556;
m_hardwareFlops[11] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[11] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 119556;
m_hardwareFlops[12] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[12] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 119556;
m_hardwareFlops[13] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[13] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 119556;
m_hardwareFlops[14] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[14] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 135972;
m_hardwareFlops[15] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[15] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 37044;
m_hardwareFlops[16] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[16] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 104328;
m_hardwareFlops[17] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[17] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 104328;
m_hardwareFlops[18] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[18] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 145692;
m_hardwareFlops[19] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[19] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 30240;
m_hardwareFlops[20] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[20] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 145692;
m_hardwareFlops[21] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[21] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 92988;
m_hardwareFlops[22] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[22] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 145044;
m_hardwareFlops[23] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[23] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 135486;
m_hardwareFlops[24] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[24] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 145044;
m_hardwareFlops[25] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[25] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 92988;
m_hardwareFlops[26] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[26] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 135486;
m_hardwareFlops[27] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[27] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 135972;
m_hardwareFlops[28] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[28] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 119556;
m_hardwareFlops[29] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[29] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 119556;
m_hardwareFlops[30] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[30] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 92988;
m_hardwareFlops[31] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[31] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 145044;
m_hardwareFlops[32] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[32] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 135486;
m_hardwareFlops[33] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[33] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 105084;
m_hardwareFlops[34] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[34] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 168696;
m_hardwareFlops[35] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[35] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 114480;
m_hardwareFlops[36] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[36] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 83160;
m_hardwareFlops[37] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[37] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 156870;
m_hardwareFlops[38] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[38] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 156870;
m_hardwareFlops[39] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[39] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 119556;
m_hardwareFlops[40] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[40] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 119556;
m_hardwareFlops[41] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[41] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 135972;
m_hardwareFlops[42] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[42] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 145044;
m_hardwareFlops[43] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[43] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 92988;
m_hardwareFlops[44] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[44] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 135486;
m_hardwareFlops[45] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[45] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 83160;
m_hardwareFlops[46] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[46] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 156870;
m_hardwareFlops[47] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[47] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 156870;
m_hardwareFlops[48] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[48] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 105084;
m_hardwareFlops[49] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[49] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 114480;
m_hardwareFlops[50] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[50] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 168696;
m_hardwareFlops[51] = 193536;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_BL2viaC;
#else
m_matrixKernels[51] = sgemm_m64_n27_k56_ldA64_ldB64_ldC64_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 9072;
m_hardwareFlops[52] = 93312;
m_matrixKernels[52] = sgemm_m64_n27_k27_ldA64_ldB27_ldC64_beta1_pfsigonly;
m_nonZeroFlops[53] = 9072;
m_hardwareFlops[53] = 93312;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = sgemm_m64_n27_k27_ldA64_ldB27_ldC64_beta1_pfsigonly;
#else
m_matrixKernels[53] = sgemm_m64_n27_k27_ldA64_ldB27_ldC64_beta1_pfsigonly;
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

#define STAR_NNZ 729

#endif

#if CONVERGENCE_ORDER==7

#ifdef TIME_KERNEL
m_nonZeroFlops[0] = 37044;
m_matrixKernels[0] = sgemm_m64_n27_k84_ldA64_ldB96_ldC64_beta0_pfsigonly;
m_hardwareFlops[0] = 290304;
m_nonZeroFlops[1] = 83916;
m_matrixKernels[1] = sgemm_m64_n27_k84_ldA64_ldB96_ldC64_beta0_pfsigonly;
m_hardwareFlops[1] = 290304;
m_nonZeroFlops[2] = 90720;
m_matrixKernels[2] = sgemm_m64_n27_k84_ldA64_ldB96_ldC64_beta0_pfsigonly;
m_hardwareFlops[2] = 290304;
m_nonZeroFlops[3] = 2688;
m_matrixKernels[3] = sgemm_m64_n27_k27_ldA64_ldB27_ldC64_beta1_pfsigonly;
m_hardwareFlops[3] = 93312;
m_nonZeroFlops[4] = 15876;
m_matrixKernels[4] = sgemm_m48_n27_k56_ldA64_ldB64_ldC48_beta0_pfsigonly;
m_hardwareFlops[4] = 145152;
m_nonZeroFlops[5] = 36288;
m_matrixKernels[5] = sgemm_m48_n27_k56_ldA64_ldB64_ldC48_beta0_pfsigonly;
m_hardwareFlops[5] = 145152;
m_nonZeroFlops[6] = 40068;
m_matrixKernels[6] = sgemm_m48_n27_k56_ldA64_ldB64_ldC48_beta0_pfsigonly;
m_hardwareFlops[6] = 145152;
m_nonZeroFlops[7] = 1680;
m_matrixKernels[7] = sgemm_m48_n27_k27_ldA48_ldB27_ldC48_beta1_pfsigonly;
m_hardwareFlops[7] = 69984;
m_nonZeroFlops[8] = 5832;
m_matrixKernels[8] = sgemm_m32_n27_k35_ldA64_ldB48_ldC32_beta0_pfsigonly;
m_hardwareFlops[8] = 60480;
m_nonZeroFlops[9] = 13608;
m_matrixKernels[9] = sgemm_m32_n27_k35_ldA64_ldB48_ldC32_beta0_pfsigonly;
m_hardwareFlops[9] = 60480;
m_nonZeroFlops[10] = 15498;
m_matrixKernels[10] = sgemm_m32_n27_k35_ldA64_ldB48_ldC32_beta0_pfsigonly;
m_hardwareFlops[10] = 60480;
m_nonZeroFlops[11] = 960;
m_matrixKernels[11] = sgemm_m32_n27_k27_ldA32_ldB27_ldC32_beta1_pfsigonly;
m_hardwareFlops[11] = 46656;
m_nonZeroFlops[12] = 1782;
m_matrixKernels[12] = sgemm_m16_n27_k20_ldA64_ldB32_ldC16_beta0_pfsigonly;
m_hardwareFlops[12] = 17280;
m_nonZeroFlops[13] = 4158;
m_matrixKernels[13] = sgemm_m16_n27_k20_ldA64_ldB32_ldC16_beta0_pfsigonly;
m_hardwareFlops[13] = 17280;
m_nonZeroFlops[14] = 4968;
m_matrixKernels[14] = sgemm_m16_n27_k20_ldA64_ldB32_ldC16_beta0_pfsigonly;
m_hardwareFlops[14] = 17280;
m_nonZeroFlops[15] = 480;
m_matrixKernels[15] = sgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
m_hardwareFlops[15] = 23328;
m_nonZeroFlops[16] = 378;
m_matrixKernels[16] = sgemm_m16_n27_k10_ldA64_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[16] = 8640;
m_nonZeroFlops[17] = 918;
m_matrixKernels[17] = sgemm_m16_n27_k10_ldA64_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[17] = 8640;
m_nonZeroFlops[18] = 1188;
m_matrixKernels[18] = sgemm_m16_n27_k10_ldA64_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[18] = 8640;
m_nonZeroFlops[19] = 192;
m_matrixKernels[19] = sgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
m_hardwareFlops[19] = 23328;
m_nonZeroFlops[20] = 54;
m_matrixKernels[20] = sgemm_m16_n27_k4_ldA64_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[20] = 3456;
m_nonZeroFlops[21] = 108;
m_matrixKernels[21] = sgemm_m16_n27_k4_ldA64_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[21] = 3456;
m_nonZeroFlops[22] = 162;
m_matrixKernels[22] = sgemm_m16_n27_k4_ldA64_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[22] = 3456;
m_nonZeroFlops[23] = 48;
m_matrixKernels[23] = sgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
m_hardwareFlops[23] = 23328;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 37044;
m_matrixKernels[0] = sgemm_m96_n27_k56_ldA96_ldB96_ldC96_beta0_pfsigonly;
m_hardwareFlops[0] = 290304;
m_nonZeroFlops[1] = 83916;
m_matrixKernels[1] = sgemm_m96_n27_k56_ldA96_ldB96_ldC96_beta0_pfsigonly;
m_hardwareFlops[1] = 290304;
m_nonZeroFlops[2] = 90720;
m_matrixKernels[2] = sgemm_m96_n27_k56_ldA96_ldB96_ldC96_beta0_pfsigonly;
m_hardwareFlops[2] = 290304;
m_nonZeroFlops[3] = 4032;
m_matrixKernels[3] = sgemm_m96_n27_k27_ldA96_ldB27_ldC96_beta1_pfsigonly;
m_hardwareFlops[3] = 139968;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 18144;
m_hardwareFlops[0] = 435456;
m_matrixKernels[0] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
m_nonZeroFlops[1] = 58968;
m_hardwareFlops[1] = 435456;
m_matrixKernels[1] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
m_nonZeroFlops[2] = 185976;
m_hardwareFlops[2] = 435456;
m_matrixKernels[2] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
m_nonZeroFlops[3] = 185976;
m_hardwareFlops[3] = 435456;
m_matrixKernels[3] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
m_nonZeroFlops[4] = 57996;
m_hardwareFlops[4] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[4] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 57996;
m_hardwareFlops[5] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[5] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 18144;
m_hardwareFlops[6] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[6] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 72576;
m_hardwareFlops[7] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[7] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 231984;
m_hardwareFlops[8] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[8] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 231984;
m_hardwareFlops[9] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[9] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 303804;
m_hardwareFlops[10] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[10] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 269406;
m_hardwareFlops[11] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[11] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 267570;
m_hardwareFlops[12] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[12] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 267570;
m_hardwareFlops[13] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[13] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 269406;
m_hardwareFlops[14] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[14] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 303804;
m_hardwareFlops[15] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[15] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 72576;
m_hardwareFlops[16] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[16] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 231984;
m_hardwareFlops[17] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[17] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 231984;
m_hardwareFlops[18] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[18] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 328860;
m_hardwareFlops[19] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[19] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 58968;
m_hardwareFlops[20] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[20] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 328860;
m_hardwareFlops[21] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[21] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 207630;
m_hardwareFlops[22] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[22] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 327726;
m_hardwareFlops[23] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[23] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 305532;
m_hardwareFlops[24] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[24] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 327726;
m_hardwareFlops[25] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[25] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 207630;
m_hardwareFlops[26] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[26] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 305532;
m_hardwareFlops[27] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[27] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 303804;
m_hardwareFlops[28] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[28] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 269406;
m_hardwareFlops[29] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[29] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 267570;
m_hardwareFlops[30] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[30] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 207630;
m_hardwareFlops[31] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[31] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 327726;
m_hardwareFlops[32] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[32] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 305532;
m_hardwareFlops[33] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[33] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 234468;
m_hardwareFlops[34] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[34] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 379188;
m_hardwareFlops[35] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[35] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 257364;
m_hardwareFlops[36] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[36] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 185976;
m_hardwareFlops[37] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[37] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 354024;
m_hardwareFlops[38] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[38] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 354024;
m_hardwareFlops[39] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[39] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 267570;
m_hardwareFlops[40] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[40] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 269406;
m_hardwareFlops[41] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[41] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 303804;
m_hardwareFlops[42] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[42] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 327726;
m_hardwareFlops[43] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[43] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 207630;
m_hardwareFlops[44] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[44] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 305532;
m_hardwareFlops[45] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[45] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 185976;
m_hardwareFlops[46] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[46] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 354024;
m_hardwareFlops[47] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[47] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 354024;
m_hardwareFlops[48] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[48] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 234468;
m_hardwareFlops[49] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[49] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 257364;
m_hardwareFlops[50] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[50] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 379188;
m_hardwareFlops[51] = 435456;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_BL2viaC;
#else
m_matrixKernels[51] = sgemm_m96_n27_k84_ldA96_ldB96_ldC96_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 13608;
m_hardwareFlops[52] = 139968;
m_matrixKernels[52] = sgemm_m96_n27_k27_ldA96_ldB27_ldC96_beta1_pfsigonly;
m_nonZeroFlops[53] = 13608;
m_hardwareFlops[53] = 139968;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = sgemm_m96_n27_k27_ldA96_ldB27_ldC96_beta1_pfsigonly;
#else
m_matrixKernels[53] = sgemm_m96_n27_k27_ldA96_ldB27_ldC96_beta1_pfsigonly;
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

#define STAR_NNZ 729

#endif

#if CONVERGENCE_ORDER==8

#ifdef TIME_KERNEL
m_nonZeroFlops[0] = 78084;
m_matrixKernels[0] = sgemm_m96_n27_k120_ldA96_ldB128_ldC96_beta0_pfsigonly;
m_hardwareFlops[0] = 622080;
m_nonZeroFlops[1] = 174312;
m_matrixKernels[1] = sgemm_m96_n27_k120_ldA96_ldB128_ldC96_beta0_pfsigonly;
m_hardwareFlops[1] = 622080;
m_nonZeroFlops[2] = 185652;
m_matrixKernels[2] = sgemm_m96_n27_k120_ldA96_ldB128_ldC96_beta0_pfsigonly;
m_hardwareFlops[2] = 622080;
m_nonZeroFlops[3] = 4032;
m_matrixKernels[3] = sgemm_m96_n27_k27_ldA96_ldB27_ldC96_beta1_pfsigonly;
m_hardwareFlops[3] = 139968;
m_nonZeroFlops[4] = 37044;
m_matrixKernels[4] = sgemm_m64_n27_k84_ldA96_ldB96_ldC64_beta0_pfsigonly;
m_hardwareFlops[4] = 290304;
m_nonZeroFlops[5] = 83916;
m_matrixKernels[5] = sgemm_m64_n27_k84_ldA96_ldB96_ldC64_beta0_pfsigonly;
m_hardwareFlops[5] = 290304;
m_nonZeroFlops[6] = 90720;
m_matrixKernels[6] = sgemm_m64_n27_k84_ldA96_ldB96_ldC64_beta0_pfsigonly;
m_hardwareFlops[6] = 290304;
m_nonZeroFlops[7] = 2688;
m_matrixKernels[7] = sgemm_m64_n27_k27_ldA64_ldB27_ldC64_beta1_pfsigonly;
m_hardwareFlops[7] = 93312;
m_nonZeroFlops[8] = 15876;
m_matrixKernels[8] = sgemm_m48_n27_k56_ldA96_ldB64_ldC48_beta0_pfsigonly;
m_hardwareFlops[8] = 145152;
m_nonZeroFlops[9] = 36288;
m_matrixKernels[9] = sgemm_m48_n27_k56_ldA96_ldB64_ldC48_beta0_pfsigonly;
m_hardwareFlops[9] = 145152;
m_nonZeroFlops[10] = 40068;
m_matrixKernels[10] = sgemm_m48_n27_k56_ldA96_ldB64_ldC48_beta0_pfsigonly;
m_hardwareFlops[10] = 145152;
m_nonZeroFlops[11] = 1680;
m_matrixKernels[11] = sgemm_m48_n27_k27_ldA48_ldB27_ldC48_beta1_pfsigonly;
m_hardwareFlops[11] = 69984;
m_nonZeroFlops[12] = 5832;
m_matrixKernels[12] = sgemm_m32_n27_k35_ldA96_ldB48_ldC32_beta0_pfsigonly;
m_hardwareFlops[12] = 60480;
m_nonZeroFlops[13] = 13608;
m_matrixKernels[13] = sgemm_m32_n27_k35_ldA96_ldB48_ldC32_beta0_pfsigonly;
m_hardwareFlops[13] = 60480;
m_nonZeroFlops[14] = 15498;
m_matrixKernels[14] = sgemm_m32_n27_k35_ldA96_ldB48_ldC32_beta0_pfsigonly;
m_hardwareFlops[14] = 60480;
m_nonZeroFlops[15] = 960;
m_matrixKernels[15] = sgemm_m32_n27_k27_ldA32_ldB27_ldC32_beta1_pfsigonly;
m_hardwareFlops[15] = 46656;
m_nonZeroFlops[16] = 1782;
m_matrixKernels[16] = sgemm_m16_n27_k20_ldA96_ldB32_ldC16_beta0_pfsigonly;
m_hardwareFlops[16] = 17280;
m_nonZeroFlops[17] = 4158;
m_matrixKernels[17] = sgemm_m16_n27_k20_ldA96_ldB32_ldC16_beta0_pfsigonly;
m_hardwareFlops[17] = 17280;
m_nonZeroFlops[18] = 4968;
m_matrixKernels[18] = sgemm_m16_n27_k20_ldA96_ldB32_ldC16_beta0_pfsigonly;
m_hardwareFlops[18] = 17280;
m_nonZeroFlops[19] = 480;
m_matrixKernels[19] = sgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
m_hardwareFlops[19] = 23328;
m_nonZeroFlops[20] = 378;
m_matrixKernels[20] = sgemm_m16_n27_k10_ldA96_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[20] = 8640;
m_nonZeroFlops[21] = 918;
m_matrixKernels[21] = sgemm_m16_n27_k10_ldA96_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[21] = 8640;
m_nonZeroFlops[22] = 1188;
m_matrixKernels[22] = sgemm_m16_n27_k10_ldA96_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[22] = 8640;
m_nonZeroFlops[23] = 192;
m_matrixKernels[23] = sgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
m_hardwareFlops[23] = 23328;
m_nonZeroFlops[24] = 54;
m_matrixKernels[24] = sgemm_m16_n27_k4_ldA96_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[24] = 3456;
m_nonZeroFlops[25] = 108;
m_matrixKernels[25] = sgemm_m16_n27_k4_ldA96_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[25] = 3456;
m_nonZeroFlops[26] = 162;
m_matrixKernels[26] = sgemm_m16_n27_k4_ldA96_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[26] = 3456;
m_nonZeroFlops[27] = 48;
m_matrixKernels[27] = sgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
m_hardwareFlops[27] = 23328;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 78084;
m_matrixKernels[0] = sgemm_m128_n27_k84_ldA128_ldB128_ldC128_beta0_pfsigonly;
m_hardwareFlops[0] = 580608;
m_nonZeroFlops[1] = 174312;
m_matrixKernels[1] = sgemm_m128_n27_k84_ldA128_ldB128_ldC128_beta0_pfsigonly;
m_hardwareFlops[1] = 580608;
m_nonZeroFlops[2] = 185652;
m_matrixKernels[2] = sgemm_m128_n27_k84_ldA128_ldB128_ldC128_beta0_pfsigonly;
m_hardwareFlops[2] = 580608;
m_nonZeroFlops[3] = 5760;
m_matrixKernels[3] = sgemm_m128_n27_k27_ldA128_ldB27_ldC128_beta1_pfsigonly;
m_hardwareFlops[3] = 186624;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 29160;
m_hardwareFlops[0] = 829440;
m_matrixKernels[0] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
m_nonZeroFlops[1] = 106272;
m_hardwareFlops[1] = 829440;
m_matrixKernels[1] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
m_nonZeroFlops[2] = 377784;
m_hardwareFlops[2] = 829440;
m_matrixKernels[2] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
m_nonZeroFlops[3] = 377784;
m_hardwareFlops[3] = 829440;
m_matrixKernels[3] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
m_nonZeroFlops[4] = 104544;
m_hardwareFlops[4] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[4] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 104544;
m_hardwareFlops[5] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[5] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 29160;
m_hardwareFlops[6] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[6] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 131220;
m_hardwareFlops[7] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[7] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 470448;
m_hardwareFlops[8] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[8] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 470448;
m_hardwareFlops[9] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[9] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 621486;
m_hardwareFlops[10] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[10] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 545940;
m_hardwareFlops[11] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[11] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 544860;
m_hardwareFlops[12] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[12] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 544860;
m_hardwareFlops[13] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[13] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 545940;
m_hardwareFlops[14] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[14] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 621486;
m_hardwareFlops[15] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[15] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 131220;
m_hardwareFlops[16] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[16] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 470448;
m_hardwareFlops[17] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[17] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 470448;
m_hardwareFlops[18] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[18] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 674028;
m_hardwareFlops[19] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[19] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 106272;
m_hardwareFlops[20] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[20] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 674028;
m_hardwareFlops[21] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[21] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 422280;
m_hardwareFlops[22] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[22] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 673164;
m_hardwareFlops[23] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[23] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 621432;
m_hardwareFlops[24] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[24] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 673164;
m_hardwareFlops[25] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[25] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 422280;
m_hardwareFlops[26] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[26] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 621432;
m_hardwareFlops[27] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[27] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 621486;
m_hardwareFlops[28] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[28] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 545940;
m_hardwareFlops[29] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[29] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 544860;
m_hardwareFlops[30] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[30] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 422280;
m_hardwareFlops[31] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[31] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 673164;
m_hardwareFlops[32] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[32] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 621432;
m_hardwareFlops[33] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[33] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 476064;
m_hardwareFlops[34] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[34] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 775332;
m_hardwareFlops[35] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[35] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 525852;
m_hardwareFlops[36] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[36] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 377784;
m_hardwareFlops[37] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[37] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 724734;
m_hardwareFlops[38] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[38] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 724734;
m_hardwareFlops[39] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[39] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 544860;
m_hardwareFlops[40] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[40] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 545940;
m_hardwareFlops[41] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[41] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 621486;
m_hardwareFlops[42] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[42] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 673164;
m_hardwareFlops[43] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[43] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 422280;
m_hardwareFlops[44] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[44] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 621432;
m_hardwareFlops[45] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[45] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 377784;
m_hardwareFlops[46] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[46] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 724734;
m_hardwareFlops[47] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[47] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 724734;
m_hardwareFlops[48] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[48] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 476064;
m_hardwareFlops[49] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[49] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 525852;
m_hardwareFlops[50] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[50] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 775332;
m_hardwareFlops[51] = 829440;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_BL2viaC;
#else
m_matrixKernels[51] = sgemm_m128_n27_k120_ldA128_ldB128_ldC128_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 19440;
m_hardwareFlops[52] = 186624;
m_matrixKernels[52] = sgemm_m128_n27_k27_ldA128_ldB27_ldC128_beta1_pfsigonly;
m_nonZeroFlops[53] = 19440;
m_hardwareFlops[53] = 186624;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = sgemm_m128_n27_k27_ldA128_ldB27_ldC128_beta1_pfsigonly;
#else
m_matrixKernels[53] = sgemm_m128_n27_k27_ldA128_ldB27_ldC128_beta1_pfsigonly;
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

#define STAR_NNZ 729

#endif

