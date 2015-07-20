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
// @date 2015-07-13 15:29:44.247587
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
m_matrixKernels[0] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
m_hardwareFlops[0] = 1728;
m_nonZeroFlops[1] = 108;
m_matrixKernels[1] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
m_hardwareFlops[1] = 1728;
m_nonZeroFlops[2] = 162;
m_matrixKernels[2] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
m_hardwareFlops[2] = 1728;
m_nonZeroFlops[3] = 48;
m_matrixKernels[3] = dgemm_m8_n27_k27_ldA8_ldB27_ldC8_beta1_pfsigonly;
m_hardwareFlops[3] = 11664;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 54;
m_matrixKernels[0] = dgemm_m8_n27_k1_ldA8_ldB8_ldC8_beta0_pfsigonly;
m_hardwareFlops[0] = 432;
m_nonZeroFlops[1] = 108;
m_matrixKernels[1] = dgemm_m8_n27_k1_ldA8_ldB8_ldC8_beta0_pfsigonly;
m_hardwareFlops[1] = 432;
m_nonZeroFlops[2] = 162;
m_matrixKernels[2] = dgemm_m8_n27_k1_ldA8_ldB8_ldC8_beta0_pfsigonly;
m_hardwareFlops[2] = 432;
m_nonZeroFlops[3] = 192;
m_matrixKernels[3] = dgemm_m8_n27_k27_ldA8_ldB27_ldC8_beta1_pfsigonly;
m_hardwareFlops[3] = 11664;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 324;
m_hardwareFlops[0] = 1728;
m_matrixKernels[0] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
m_nonZeroFlops[1] = 432;
m_hardwareFlops[1] = 1728;
m_matrixKernels[1] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
m_nonZeroFlops[2] = 540;
m_hardwareFlops[2] = 1728;
m_matrixKernels[2] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
m_nonZeroFlops[3] = 540;
m_hardwareFlops[3] = 1728;
m_matrixKernels[3] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
m_nonZeroFlops[4] = 432;
m_hardwareFlops[4] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[4] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 432;
m_hardwareFlops[5] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[5] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 324;
m_hardwareFlops[6] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[6] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 486;
m_hardwareFlops[7] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[7] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 648;
m_hardwareFlops[8] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[8] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 648;
m_hardwareFlops[9] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[9] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 756;
m_hardwareFlops[10] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[10] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 702;
m_hardwareFlops[11] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[11] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 702;
m_hardwareFlops[12] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[12] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 702;
m_hardwareFlops[13] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[13] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 702;
m_hardwareFlops[14] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[14] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 756;
m_hardwareFlops[15] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[15] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 486;
m_hardwareFlops[16] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[16] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 648;
m_hardwareFlops[17] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[17] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 648;
m_hardwareFlops[18] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[18] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 756;
m_hardwareFlops[19] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[19] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 432;
m_hardwareFlops[20] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[20] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 756;
m_hardwareFlops[21] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[21] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 594;
m_hardwareFlops[22] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[22] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 756;
m_hardwareFlops[23] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[23] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 756;
m_hardwareFlops[24] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[24] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 756;
m_hardwareFlops[25] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[25] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 594;
m_hardwareFlops[26] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[26] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 756;
m_hardwareFlops[27] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[27] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 756;
m_hardwareFlops[28] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[28] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 702;
m_hardwareFlops[29] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[29] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 702;
m_hardwareFlops[30] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[30] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 594;
m_hardwareFlops[31] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[31] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 756;
m_hardwareFlops[32] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[32] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 756;
m_hardwareFlops[33] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[33] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 648;
m_hardwareFlops[34] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[34] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 864;
m_hardwareFlops[35] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[35] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 648;
m_hardwareFlops[36] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[36] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 540;
m_hardwareFlops[37] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[37] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 810;
m_hardwareFlops[38] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[38] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 810;
m_hardwareFlops[39] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[39] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 702;
m_hardwareFlops[40] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[40] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 702;
m_hardwareFlops[41] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[41] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 756;
m_hardwareFlops[42] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[42] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 756;
m_hardwareFlops[43] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[43] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 594;
m_hardwareFlops[44] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[44] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 756;
m_hardwareFlops[45] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[45] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 540;
m_hardwareFlops[46] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[46] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 810;
m_hardwareFlops[47] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[47] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 810;
m_hardwareFlops[48] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[48] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 648;
m_hardwareFlops[49] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[49] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 648;
m_hardwareFlops[50] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[50] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 864;
m_hardwareFlops[51] = 1728;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[51] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 648;
m_hardwareFlops[52] = 11664;
m_matrixKernels[52] = dgemm_m8_n27_k27_ldA8_ldB27_ldC8_beta1_pfsigonly;
m_nonZeroFlops[53] = 648;
m_hardwareFlops[53] = 11664;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = dgemm_m8_n27_k27_ldA8_ldB27_ldC8_beta1_AL2jpst_BL2viaC;
#else
m_matrixKernels[53] = dgemm_m8_n27_k27_ldA8_ldB27_ldC8_beta1_pfsigonly;
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
m_matrixKernels[0] = dgemm_m8_n27_k10_ldA8_ldB16_ldC8_beta0_pfsigonly;
m_hardwareFlops[0] = 4320;
m_nonZeroFlops[1] = 918;
m_matrixKernels[1] = dgemm_m8_n27_k10_ldA8_ldB16_ldC8_beta0_pfsigonly;
m_hardwareFlops[1] = 4320;
m_nonZeroFlops[2] = 1188;
m_matrixKernels[2] = dgemm_m8_n27_k10_ldA8_ldB16_ldC8_beta0_pfsigonly;
m_hardwareFlops[2] = 4320;
m_nonZeroFlops[3] = 192;
m_matrixKernels[3] = dgemm_m8_n27_k27_ldA8_ldB27_ldC8_beta1_pfsigonly;
m_hardwareFlops[3] = 11664;
m_nonZeroFlops[4] = 54;
m_matrixKernels[4] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
m_hardwareFlops[4] = 1728;
m_nonZeroFlops[5] = 108;
m_matrixKernels[5] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
m_hardwareFlops[5] = 1728;
m_nonZeroFlops[6] = 162;
m_matrixKernels[6] = dgemm_m8_n27_k4_ldA8_ldB8_ldC8_beta0_pfsigonly;
m_hardwareFlops[6] = 1728;
m_nonZeroFlops[7] = 48;
m_matrixKernels[7] = dgemm_m8_n27_k27_ldA8_ldB27_ldC8_beta1_pfsigonly;
m_hardwareFlops[7] = 11664;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 378;
m_matrixKernels[0] = dgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[0] = 3456;
m_nonZeroFlops[1] = 918;
m_matrixKernels[1] = dgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[1] = 3456;
m_nonZeroFlops[2] = 1188;
m_matrixKernels[2] = dgemm_m16_n27_k4_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_hardwareFlops[2] = 3456;
m_nonZeroFlops[3] = 480;
m_matrixKernels[3] = dgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
m_hardwareFlops[3] = 23328;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 1080;
m_hardwareFlops[0] = 8640;
m_matrixKernels[0] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_nonZeroFlops[1] = 1836;
m_hardwareFlops[1] = 8640;
m_matrixKernels[1] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_nonZeroFlops[2] = 2916;
m_hardwareFlops[2] = 8640;
m_matrixKernels[2] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_nonZeroFlops[3] = 2916;
m_hardwareFlops[3] = 8640;
m_matrixKernels[3] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
m_nonZeroFlops[4] = 1836;
m_hardwareFlops[4] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[4] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 1836;
m_hardwareFlops[5] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[5] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 1080;
m_hardwareFlops[6] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[6] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 2160;
m_hardwareFlops[7] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[7] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 3672;
m_hardwareFlops[8] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[8] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 3672;
m_hardwareFlops[9] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[9] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 4374;
m_hardwareFlops[10] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[10] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 4104;
m_hardwareFlops[11] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[11] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 4104;
m_hardwareFlops[12] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[12] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 4104;
m_hardwareFlops[13] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[13] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 4104;
m_hardwareFlops[14] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[14] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 4374;
m_hardwareFlops[15] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[15] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 2160;
m_hardwareFlops[16] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[16] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 3672;
m_hardwareFlops[17] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[17] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 3672;
m_hardwareFlops[18] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[18] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 4536;
m_hardwareFlops[19] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[19] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 1836;
m_hardwareFlops[20] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[20] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 4536;
m_hardwareFlops[21] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[21] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 3294;
m_hardwareFlops[22] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[22] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 4536;
m_hardwareFlops[23] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[23] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 4536;
m_hardwareFlops[24] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[24] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 4536;
m_hardwareFlops[25] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[25] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 3294;
m_hardwareFlops[26] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[26] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 4536;
m_hardwareFlops[27] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[27] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 4374;
m_hardwareFlops[28] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[28] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 4104;
m_hardwareFlops[29] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[29] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 4104;
m_hardwareFlops[30] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[30] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 3294;
m_hardwareFlops[31] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[31] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 4536;
m_hardwareFlops[32] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[32] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 4536;
m_hardwareFlops[33] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[33] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 3672;
m_hardwareFlops[34] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[34] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 5292;
m_hardwareFlops[35] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[35] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 3780;
m_hardwareFlops[36] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[36] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 2916;
m_hardwareFlops[37] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[37] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 5022;
m_hardwareFlops[38] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[38] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 5022;
m_hardwareFlops[39] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[39] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 4104;
m_hardwareFlops[40] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[40] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 4104;
m_hardwareFlops[41] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[41] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 4374;
m_hardwareFlops[42] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[42] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 4536;
m_hardwareFlops[43] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[43] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 3294;
m_hardwareFlops[44] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[44] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 4536;
m_hardwareFlops[45] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[45] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 2916;
m_hardwareFlops[46] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[46] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 5022;
m_hardwareFlops[47] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[47] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 5022;
m_hardwareFlops[48] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[48] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 3672;
m_hardwareFlops[49] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[49] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 3780;
m_hardwareFlops[50] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[50] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 5292;
m_hardwareFlops[51] = 8640;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[51] = dgemm_m16_n27_k10_ldA16_ldB16_ldC16_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 1620;
m_hardwareFlops[52] = 23328;
m_matrixKernels[52] = dgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
m_nonZeroFlops[53] = 1620;
m_hardwareFlops[53] = 23328;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = dgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_AL2jpst_BL2viaC;
#else
m_matrixKernels[53] = dgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
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
m_matrixKernels[0] = dgemm_m16_n27_k20_ldA16_ldB24_ldC16_beta0_pfsigonly;
m_hardwareFlops[0] = 17280;
m_nonZeroFlops[1] = 4158;
m_matrixKernels[1] = dgemm_m16_n27_k20_ldA16_ldB24_ldC16_beta0_pfsigonly;
m_hardwareFlops[1] = 17280;
m_nonZeroFlops[2] = 4968;
m_matrixKernels[2] = dgemm_m16_n27_k20_ldA16_ldB24_ldC16_beta0_pfsigonly;
m_hardwareFlops[2] = 17280;
m_nonZeroFlops[3] = 480;
m_matrixKernels[3] = dgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
m_hardwareFlops[3] = 23328;
m_nonZeroFlops[4] = 378;
m_matrixKernels[4] = dgemm_m8_n27_k10_ldA16_ldB16_ldC8_beta0_pfsigonly;
m_hardwareFlops[4] = 4320;
m_nonZeroFlops[5] = 918;
m_matrixKernels[5] = dgemm_m8_n27_k10_ldA16_ldB16_ldC8_beta0_pfsigonly;
m_hardwareFlops[5] = 4320;
m_nonZeroFlops[6] = 1188;
m_matrixKernels[6] = dgemm_m8_n27_k10_ldA16_ldB16_ldC8_beta0_pfsigonly;
m_hardwareFlops[6] = 4320;
m_nonZeroFlops[7] = 192;
m_matrixKernels[7] = dgemm_m8_n27_k27_ldA8_ldB27_ldC8_beta1_pfsigonly;
m_hardwareFlops[7] = 11664;
m_nonZeroFlops[8] = 54;
m_matrixKernels[8] = dgemm_m8_n27_k4_ldA16_ldB8_ldC8_beta0_pfsigonly;
m_hardwareFlops[8] = 1728;
m_nonZeroFlops[9] = 108;
m_matrixKernels[9] = dgemm_m8_n27_k4_ldA16_ldB8_ldC8_beta0_pfsigonly;
m_hardwareFlops[9] = 1728;
m_nonZeroFlops[10] = 162;
m_matrixKernels[10] = dgemm_m8_n27_k4_ldA16_ldB8_ldC8_beta0_pfsigonly;
m_hardwareFlops[10] = 1728;
m_nonZeroFlops[11] = 48;
m_matrixKernels[11] = dgemm_m8_n27_k27_ldA8_ldB27_ldC8_beta1_pfsigonly;
m_hardwareFlops[11] = 11664;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 1782;
m_matrixKernels[0] = dgemm_m24_n27_k10_ldA24_ldB24_ldC24_beta0_pfsigonly;
m_hardwareFlops[0] = 12960;
m_nonZeroFlops[1] = 4158;
m_matrixKernels[1] = dgemm_m24_n27_k10_ldA24_ldB24_ldC24_beta0_pfsigonly;
m_hardwareFlops[1] = 12960;
m_nonZeroFlops[2] = 4968;
m_matrixKernels[2] = dgemm_m24_n27_k10_ldA24_ldB24_ldC24_beta0_pfsigonly;
m_hardwareFlops[2] = 12960;
m_nonZeroFlops[3] = 960;
m_matrixKernels[3] = dgemm_m24_n27_k27_ldA24_ldB27_ldC24_beta1_pfsigonly;
m_hardwareFlops[3] = 34992;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 2700;
m_hardwareFlops[0] = 25920;
m_matrixKernels[0] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
m_nonZeroFlops[1] = 5616;
m_hardwareFlops[1] = 25920;
m_matrixKernels[1] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
m_nonZeroFlops[2] = 11016;
m_hardwareFlops[2] = 25920;
m_matrixKernels[2] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
m_nonZeroFlops[3] = 11016;
m_hardwareFlops[3] = 25920;
m_matrixKernels[3] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
m_nonZeroFlops[4] = 5616;
m_hardwareFlops[4] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[4] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 5616;
m_hardwareFlops[5] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[5] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 2700;
m_hardwareFlops[6] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[6] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 6750;
m_hardwareFlops[7] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[7] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 14040;
m_hardwareFlops[8] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[8] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 14040;
m_hardwareFlops[9] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[9] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 17280;
m_hardwareFlops[10] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[10] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 15444;
m_hardwareFlops[11] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[11] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 15930;
m_hardwareFlops[12] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[12] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 15930;
m_hardwareFlops[13] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[13] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 15444;
m_hardwareFlops[14] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[14] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 17280;
m_hardwareFlops[15] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[15] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 6750;
m_hardwareFlops[16] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[16] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 14040;
m_hardwareFlops[17] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[17] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 14040;
m_hardwareFlops[18] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[18] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 18252;
m_hardwareFlops[19] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[19] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 5616;
m_hardwareFlops[20] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[20] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 18252;
m_hardwareFlops[21] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[21] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 12528;
m_hardwareFlops[22] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[22] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 18360;
m_hardwareFlops[23] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[23] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 17550;
m_hardwareFlops[24] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[24] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 18360;
m_hardwareFlops[25] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[25] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 12528;
m_hardwareFlops[26] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[26] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 17550;
m_hardwareFlops[27] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[27] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 17280;
m_hardwareFlops[28] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[28] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 15444;
m_hardwareFlops[29] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[29] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 15930;
m_hardwareFlops[30] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[30] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 12528;
m_hardwareFlops[31] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[31] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 18360;
m_hardwareFlops[32] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[32] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 17550;
m_hardwareFlops[33] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[33] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 13932;
m_hardwareFlops[34] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[34] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 21276;
m_hardwareFlops[35] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[35] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 14796;
m_hardwareFlops[36] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[36] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 11016;
m_hardwareFlops[37] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[37] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 20034;
m_hardwareFlops[38] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[38] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 20034;
m_hardwareFlops[39] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[39] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 15930;
m_hardwareFlops[40] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[40] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 15444;
m_hardwareFlops[41] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[41] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 17280;
m_hardwareFlops[42] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[42] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 18360;
m_hardwareFlops[43] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[43] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 12528;
m_hardwareFlops[44] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[44] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 17550;
m_hardwareFlops[45] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[45] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 11016;
m_hardwareFlops[46] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[46] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 20034;
m_hardwareFlops[47] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[47] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 20034;
m_hardwareFlops[48] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[48] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 13932;
m_hardwareFlops[49] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[49] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 14796;
m_hardwareFlops[50] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[50] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 21276;
m_hardwareFlops[51] = 25920;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[51] = dgemm_m24_n27_k20_ldA24_ldB24_ldC24_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 3240;
m_hardwareFlops[52] = 34992;
m_matrixKernels[52] = dgemm_m24_n27_k27_ldA24_ldB27_ldC24_beta1_pfsigonly;
m_nonZeroFlops[53] = 3240;
m_hardwareFlops[53] = 34992;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = dgemm_m24_n27_k27_ldA24_ldB27_ldC24_beta1_AL2jpst_BL2viaC;
#else
m_matrixKernels[53] = dgemm_m24_n27_k27_ldA24_ldB27_ldC24_beta1_pfsigonly;
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
m_matrixKernels[0] = dgemm_m24_n27_k35_ldA24_ldB40_ldC24_beta0_pfsigonly;
m_hardwareFlops[0] = 45360;
m_nonZeroFlops[1] = 13608;
m_matrixKernels[1] = dgemm_m24_n27_k35_ldA24_ldB40_ldC24_beta0_pfsigonly;
m_hardwareFlops[1] = 45360;
m_nonZeroFlops[2] = 15498;
m_matrixKernels[2] = dgemm_m24_n27_k35_ldA24_ldB40_ldC24_beta0_pfsigonly;
m_hardwareFlops[2] = 45360;
m_nonZeroFlops[3] = 960;
m_matrixKernels[3] = dgemm_m24_n27_k27_ldA24_ldB27_ldC24_beta1_pfsigonly;
m_hardwareFlops[3] = 34992;
m_nonZeroFlops[4] = 1782;
m_matrixKernels[4] = dgemm_m16_n27_k20_ldA24_ldB24_ldC16_beta0_pfsigonly;
m_hardwareFlops[4] = 17280;
m_nonZeroFlops[5] = 4158;
m_matrixKernels[5] = dgemm_m16_n27_k20_ldA24_ldB24_ldC16_beta0_pfsigonly;
m_hardwareFlops[5] = 17280;
m_nonZeroFlops[6] = 4968;
m_matrixKernels[6] = dgemm_m16_n27_k20_ldA24_ldB24_ldC16_beta0_pfsigonly;
m_hardwareFlops[6] = 17280;
m_nonZeroFlops[7] = 480;
m_matrixKernels[7] = dgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
m_hardwareFlops[7] = 23328;
m_nonZeroFlops[8] = 378;
m_matrixKernels[8] = dgemm_m8_n27_k10_ldA24_ldB16_ldC8_beta0_pfsigonly;
m_hardwareFlops[8] = 4320;
m_nonZeroFlops[9] = 918;
m_matrixKernels[9] = dgemm_m8_n27_k10_ldA24_ldB16_ldC8_beta0_pfsigonly;
m_hardwareFlops[9] = 4320;
m_nonZeroFlops[10] = 1188;
m_matrixKernels[10] = dgemm_m8_n27_k10_ldA24_ldB16_ldC8_beta0_pfsigonly;
m_hardwareFlops[10] = 4320;
m_nonZeroFlops[11] = 192;
m_matrixKernels[11] = dgemm_m8_n27_k27_ldA8_ldB27_ldC8_beta1_pfsigonly;
m_hardwareFlops[11] = 11664;
m_nonZeroFlops[12] = 54;
m_matrixKernels[12] = dgemm_m8_n27_k4_ldA24_ldB8_ldC8_beta0_pfsigonly;
m_hardwareFlops[12] = 1728;
m_nonZeroFlops[13] = 108;
m_matrixKernels[13] = dgemm_m8_n27_k4_ldA24_ldB8_ldC8_beta0_pfsigonly;
m_hardwareFlops[13] = 1728;
m_nonZeroFlops[14] = 162;
m_matrixKernels[14] = dgemm_m8_n27_k4_ldA24_ldB8_ldC8_beta0_pfsigonly;
m_hardwareFlops[14] = 1728;
m_nonZeroFlops[15] = 48;
m_matrixKernels[15] = dgemm_m8_n27_k27_ldA8_ldB27_ldC8_beta1_pfsigonly;
m_hardwareFlops[15] = 11664;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 5832;
m_matrixKernels[0] = dgemm_m40_n27_k20_ldA40_ldB40_ldC40_beta0_pfsigonly;
m_hardwareFlops[0] = 43200;
m_nonZeroFlops[1] = 13608;
m_matrixKernels[1] = dgemm_m40_n27_k20_ldA40_ldB40_ldC40_beta0_pfsigonly;
m_hardwareFlops[1] = 43200;
m_nonZeroFlops[2] = 15498;
m_matrixKernels[2] = dgemm_m40_n27_k20_ldA40_ldB40_ldC40_beta0_pfsigonly;
m_hardwareFlops[2] = 43200;
m_nonZeroFlops[3] = 1680;
m_matrixKernels[3] = dgemm_m40_n27_k27_ldA40_ldB27_ldC40_beta1_pfsigonly;
m_hardwareFlops[3] = 58320;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 5670;
m_hardwareFlops[0] = 75600;
m_matrixKernels[0] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
m_nonZeroFlops[1] = 13986;
m_hardwareFlops[1] = 75600;
m_matrixKernels[1] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
m_nonZeroFlops[2] = 32886;
m_hardwareFlops[2] = 75600;
m_matrixKernels[2] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
m_nonZeroFlops[3] = 32886;
m_hardwareFlops[3] = 75600;
m_matrixKernels[3] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
m_nonZeroFlops[4] = 13878;
m_hardwareFlops[4] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[4] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 13878;
m_hardwareFlops[5] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[5] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 5670;
m_hardwareFlops[6] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[6] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 17010;
m_hardwareFlops[7] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[7] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 41634;
m_hardwareFlops[8] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[8] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 41634;
m_hardwareFlops[9] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[9] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 52866;
m_hardwareFlops[10] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[10] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 46440;
m_hardwareFlops[11] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[11] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 47304;
m_hardwareFlops[12] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[12] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 47304;
m_hardwareFlops[13] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[13] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 46440;
m_hardwareFlops[14] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[14] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 52866;
m_hardwareFlops[15] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[15] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 17010;
m_hardwareFlops[16] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[16] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 41634;
m_hardwareFlops[17] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[17] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 41634;
m_hardwareFlops[18] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[18] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 56430;
m_hardwareFlops[19] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[19] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 13986;
m_hardwareFlops[20] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[20] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 56430;
m_hardwareFlops[21] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[21] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 36936;
m_hardwareFlops[22] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[22] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 56106;
m_hardwareFlops[23] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[23] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 52920;
m_hardwareFlops[24] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[24] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 56106;
m_hardwareFlops[25] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[25] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 36936;
m_hardwareFlops[26] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[26] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 52920;
m_hardwareFlops[27] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[27] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 52866;
m_hardwareFlops[28] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[28] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 46440;
m_hardwareFlops[29] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[29] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 47304;
m_hardwareFlops[30] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[30] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 36936;
m_hardwareFlops[31] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[31] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 56106;
m_hardwareFlops[32] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[32] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 52920;
m_hardwareFlops[33] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[33] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 41742;
m_hardwareFlops[34] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[34] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 65502;
m_hardwareFlops[35] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[35] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 44874;
m_hardwareFlops[36] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[36] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 32886;
m_hardwareFlops[37] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[37] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 61074;
m_hardwareFlops[38] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[38] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 61074;
m_hardwareFlops[39] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[39] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 47304;
m_hardwareFlops[40] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[40] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 46440;
m_hardwareFlops[41] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[41] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 52866;
m_hardwareFlops[42] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[42] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 56106;
m_hardwareFlops[43] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[43] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 36936;
m_hardwareFlops[44] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[44] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 52920;
m_hardwareFlops[45] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[45] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 32886;
m_hardwareFlops[46] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[46] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 61074;
m_hardwareFlops[47] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[47] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 61074;
m_hardwareFlops[48] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[48] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 41742;
m_hardwareFlops[49] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[49] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 44874;
m_hardwareFlops[50] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[50] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 65502;
m_hardwareFlops[51] = 75600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[51] = dgemm_m40_n27_k35_ldA40_ldB40_ldC40_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 5670;
m_hardwareFlops[52] = 58320;
m_matrixKernels[52] = dgemm_m40_n27_k27_ldA40_ldB27_ldC40_beta1_pfsigonly;
m_nonZeroFlops[53] = 5670;
m_hardwareFlops[53] = 58320;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = dgemm_m40_n27_k27_ldA40_ldB27_ldC40_beta1_AL2jpst_BL2viaC;
#else
m_matrixKernels[53] = dgemm_m40_n27_k27_ldA40_ldB27_ldC40_beta1_pfsigonly;
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
m_matrixKernels[0] = dgemm_m40_n27_k56_ldA40_ldB56_ldC40_beta0_pfsigonly;
m_hardwareFlops[0] = 120960;
m_nonZeroFlops[1] = 36288;
m_matrixKernels[1] = dgemm_m40_n27_k56_ldA40_ldB56_ldC40_beta0_pfsigonly;
m_hardwareFlops[1] = 120960;
m_nonZeroFlops[2] = 40068;
m_matrixKernels[2] = dgemm_m40_n27_k56_ldA40_ldB56_ldC40_beta0_pfsigonly;
m_hardwareFlops[2] = 120960;
m_nonZeroFlops[3] = 1680;
m_matrixKernels[3] = dgemm_m40_n27_k27_ldA40_ldB27_ldC40_beta1_pfsigonly;
m_hardwareFlops[3] = 58320;
m_nonZeroFlops[4] = 5832;
m_matrixKernels[4] = dgemm_m24_n27_k35_ldA40_ldB40_ldC24_beta0_pfsigonly;
m_hardwareFlops[4] = 45360;
m_nonZeroFlops[5] = 13608;
m_matrixKernels[5] = dgemm_m24_n27_k35_ldA40_ldB40_ldC24_beta0_pfsigonly;
m_hardwareFlops[5] = 45360;
m_nonZeroFlops[6] = 15498;
m_matrixKernels[6] = dgemm_m24_n27_k35_ldA40_ldB40_ldC24_beta0_pfsigonly;
m_hardwareFlops[6] = 45360;
m_nonZeroFlops[7] = 960;
m_matrixKernels[7] = dgemm_m24_n27_k27_ldA24_ldB27_ldC24_beta1_pfsigonly;
m_hardwareFlops[7] = 34992;
m_nonZeroFlops[8] = 1782;
m_matrixKernels[8] = dgemm_m16_n27_k20_ldA40_ldB24_ldC16_beta0_pfsigonly;
m_hardwareFlops[8] = 17280;
m_nonZeroFlops[9] = 4158;
m_matrixKernels[9] = dgemm_m16_n27_k20_ldA40_ldB24_ldC16_beta0_pfsigonly;
m_hardwareFlops[9] = 17280;
m_nonZeroFlops[10] = 4968;
m_matrixKernels[10] = dgemm_m16_n27_k20_ldA40_ldB24_ldC16_beta0_pfsigonly;
m_hardwareFlops[10] = 17280;
m_nonZeroFlops[11] = 480;
m_matrixKernels[11] = dgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
m_hardwareFlops[11] = 23328;
m_nonZeroFlops[12] = 378;
m_matrixKernels[12] = dgemm_m8_n27_k10_ldA40_ldB16_ldC8_beta0_pfsigonly;
m_hardwareFlops[12] = 4320;
m_nonZeroFlops[13] = 918;
m_matrixKernels[13] = dgemm_m8_n27_k10_ldA40_ldB16_ldC8_beta0_pfsigonly;
m_hardwareFlops[13] = 4320;
m_nonZeroFlops[14] = 1188;
m_matrixKernels[14] = dgemm_m8_n27_k10_ldA40_ldB16_ldC8_beta0_pfsigonly;
m_hardwareFlops[14] = 4320;
m_nonZeroFlops[15] = 192;
m_matrixKernels[15] = dgemm_m8_n27_k27_ldA8_ldB27_ldC8_beta1_pfsigonly;
m_hardwareFlops[15] = 11664;
m_nonZeroFlops[16] = 54;
m_matrixKernels[16] = dgemm_m8_n27_k4_ldA40_ldB8_ldC8_beta0_pfsigonly;
m_hardwareFlops[16] = 1728;
m_nonZeroFlops[17] = 108;
m_matrixKernels[17] = dgemm_m8_n27_k4_ldA40_ldB8_ldC8_beta0_pfsigonly;
m_hardwareFlops[17] = 1728;
m_nonZeroFlops[18] = 162;
m_matrixKernels[18] = dgemm_m8_n27_k4_ldA40_ldB8_ldC8_beta0_pfsigonly;
m_hardwareFlops[18] = 1728;
m_nonZeroFlops[19] = 48;
m_matrixKernels[19] = dgemm_m8_n27_k27_ldA8_ldB27_ldC8_beta1_pfsigonly;
m_hardwareFlops[19] = 11664;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 15876;
m_matrixKernels[0] = dgemm_m56_n27_k35_ldA56_ldB56_ldC56_beta0_pfsigonly;
m_hardwareFlops[0] = 105840;
m_nonZeroFlops[1] = 36288;
m_matrixKernels[1] = dgemm_m56_n27_k35_ldA56_ldB56_ldC56_beta0_pfsigonly;
m_hardwareFlops[1] = 105840;
m_nonZeroFlops[2] = 40068;
m_matrixKernels[2] = dgemm_m56_n27_k35_ldA56_ldB56_ldC56_beta0_pfsigonly;
m_hardwareFlops[2] = 105840;
m_nonZeroFlops[3] = 2688;
m_matrixKernels[3] = dgemm_m56_n27_k27_ldA56_ldB27_ldC56_beta1_pfsigonly;
m_hardwareFlops[3] = 81648;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 10584;
m_hardwareFlops[0] = 169344;
m_matrixKernels[0] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
m_nonZeroFlops[1] = 30240;
m_hardwareFlops[1] = 169344;
m_matrixKernels[1] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
m_nonZeroFlops[2] = 83160;
m_hardwareFlops[2] = 169344;
m_matrixKernels[2] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
m_nonZeroFlops[3] = 83160;
m_hardwareFlops[3] = 169344;
m_matrixKernels[3] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
m_nonZeroFlops[4] = 29808;
m_hardwareFlops[4] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[4] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 29808;
m_hardwareFlops[5] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[5] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 10584;
m_hardwareFlops[6] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[6] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 37044;
m_hardwareFlops[7] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[7] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 104328;
m_hardwareFlops[8] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[8] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 104328;
m_hardwareFlops[9] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[9] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 135972;
m_hardwareFlops[10] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[10] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 119556;
m_hardwareFlops[11] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[11] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 119556;
m_hardwareFlops[12] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[12] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 119556;
m_hardwareFlops[13] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[13] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 119556;
m_hardwareFlops[14] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[14] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 135972;
m_hardwareFlops[15] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[15] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 37044;
m_hardwareFlops[16] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[16] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 104328;
m_hardwareFlops[17] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[17] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 104328;
m_hardwareFlops[18] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[18] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 145692;
m_hardwareFlops[19] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[19] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 30240;
m_hardwareFlops[20] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[20] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 145692;
m_hardwareFlops[21] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[21] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 92988;
m_hardwareFlops[22] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[22] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 145044;
m_hardwareFlops[23] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[23] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 135486;
m_hardwareFlops[24] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[24] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 145044;
m_hardwareFlops[25] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[25] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 92988;
m_hardwareFlops[26] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[26] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 135486;
m_hardwareFlops[27] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[27] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 135972;
m_hardwareFlops[28] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[28] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 119556;
m_hardwareFlops[29] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[29] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 119556;
m_hardwareFlops[30] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[30] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 92988;
m_hardwareFlops[31] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[31] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 145044;
m_hardwareFlops[32] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[32] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 135486;
m_hardwareFlops[33] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[33] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 105084;
m_hardwareFlops[34] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[34] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 168696;
m_hardwareFlops[35] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[35] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 114480;
m_hardwareFlops[36] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[36] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 83160;
m_hardwareFlops[37] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[37] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 156870;
m_hardwareFlops[38] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[38] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 156870;
m_hardwareFlops[39] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[39] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 119556;
m_hardwareFlops[40] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[40] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 119556;
m_hardwareFlops[41] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[41] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 135972;
m_hardwareFlops[42] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[42] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 145044;
m_hardwareFlops[43] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[43] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 92988;
m_hardwareFlops[44] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[44] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 135486;
m_hardwareFlops[45] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[45] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 83160;
m_hardwareFlops[46] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[46] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 156870;
m_hardwareFlops[47] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[47] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 156870;
m_hardwareFlops[48] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[48] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 105084;
m_hardwareFlops[49] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[49] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 114480;
m_hardwareFlops[50] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[50] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 168696;
m_hardwareFlops[51] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[51] = dgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 9072;
m_hardwareFlops[52] = 81648;
m_matrixKernels[52] = dgemm_m56_n27_k27_ldA56_ldB27_ldC56_beta1_pfsigonly;
m_nonZeroFlops[53] = 9072;
m_hardwareFlops[53] = 81648;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = dgemm_m56_n27_k27_ldA56_ldB27_ldC56_beta1_AL2jpst_BL2viaC;
#else
m_matrixKernels[53] = dgemm_m56_n27_k27_ldA56_ldB27_ldC56_beta1_pfsigonly;
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
m_matrixKernels[0] = dgemm_m56_n27_k84_ldA56_ldB88_ldC56_beta0_pfsigonly;
m_hardwareFlops[0] = 254016;
m_nonZeroFlops[1] = 83916;
m_matrixKernels[1] = dgemm_m56_n27_k84_ldA56_ldB88_ldC56_beta0_pfsigonly;
m_hardwareFlops[1] = 254016;
m_nonZeroFlops[2] = 90720;
m_matrixKernels[2] = dgemm_m56_n27_k84_ldA56_ldB88_ldC56_beta0_pfsigonly;
m_hardwareFlops[2] = 254016;
m_nonZeroFlops[3] = 2688;
m_matrixKernels[3] = dgemm_m56_n27_k27_ldA56_ldB27_ldC56_beta1_pfsigonly;
m_hardwareFlops[3] = 81648;
m_nonZeroFlops[4] = 15876;
m_matrixKernels[4] = dgemm_m40_n27_k56_ldA56_ldB56_ldC40_beta0_pfsigonly;
m_hardwareFlops[4] = 120960;
m_nonZeroFlops[5] = 36288;
m_matrixKernels[5] = dgemm_m40_n27_k56_ldA56_ldB56_ldC40_beta0_pfsigonly;
m_hardwareFlops[5] = 120960;
m_nonZeroFlops[6] = 40068;
m_matrixKernels[6] = dgemm_m40_n27_k56_ldA56_ldB56_ldC40_beta0_pfsigonly;
m_hardwareFlops[6] = 120960;
m_nonZeroFlops[7] = 1680;
m_matrixKernels[7] = dgemm_m40_n27_k27_ldA40_ldB27_ldC40_beta1_pfsigonly;
m_hardwareFlops[7] = 58320;
m_nonZeroFlops[8] = 5832;
m_matrixKernels[8] = dgemm_m24_n27_k35_ldA56_ldB40_ldC24_beta0_pfsigonly;
m_hardwareFlops[8] = 45360;
m_nonZeroFlops[9] = 13608;
m_matrixKernels[9] = dgemm_m24_n27_k35_ldA56_ldB40_ldC24_beta0_pfsigonly;
m_hardwareFlops[9] = 45360;
m_nonZeroFlops[10] = 15498;
m_matrixKernels[10] = dgemm_m24_n27_k35_ldA56_ldB40_ldC24_beta0_pfsigonly;
m_hardwareFlops[10] = 45360;
m_nonZeroFlops[11] = 960;
m_matrixKernels[11] = dgemm_m24_n27_k27_ldA24_ldB27_ldC24_beta1_pfsigonly;
m_hardwareFlops[11] = 34992;
m_nonZeroFlops[12] = 1782;
m_matrixKernels[12] = dgemm_m16_n27_k20_ldA56_ldB24_ldC16_beta0_pfsigonly;
m_hardwareFlops[12] = 17280;
m_nonZeroFlops[13] = 4158;
m_matrixKernels[13] = dgemm_m16_n27_k20_ldA56_ldB24_ldC16_beta0_pfsigonly;
m_hardwareFlops[13] = 17280;
m_nonZeroFlops[14] = 4968;
m_matrixKernels[14] = dgemm_m16_n27_k20_ldA56_ldB24_ldC16_beta0_pfsigonly;
m_hardwareFlops[14] = 17280;
m_nonZeroFlops[15] = 480;
m_matrixKernels[15] = dgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
m_hardwareFlops[15] = 23328;
m_nonZeroFlops[16] = 378;
m_matrixKernels[16] = dgemm_m8_n27_k10_ldA56_ldB16_ldC8_beta0_pfsigonly;
m_hardwareFlops[16] = 4320;
m_nonZeroFlops[17] = 918;
m_matrixKernels[17] = dgemm_m8_n27_k10_ldA56_ldB16_ldC8_beta0_pfsigonly;
m_hardwareFlops[17] = 4320;
m_nonZeroFlops[18] = 1188;
m_matrixKernels[18] = dgemm_m8_n27_k10_ldA56_ldB16_ldC8_beta0_pfsigonly;
m_hardwareFlops[18] = 4320;
m_nonZeroFlops[19] = 192;
m_matrixKernels[19] = dgemm_m8_n27_k27_ldA8_ldB27_ldC8_beta1_pfsigonly;
m_hardwareFlops[19] = 11664;
m_nonZeroFlops[20] = 54;
m_matrixKernels[20] = dgemm_m8_n27_k4_ldA56_ldB8_ldC8_beta0_pfsigonly;
m_hardwareFlops[20] = 1728;
m_nonZeroFlops[21] = 108;
m_matrixKernels[21] = dgemm_m8_n27_k4_ldA56_ldB8_ldC8_beta0_pfsigonly;
m_hardwareFlops[21] = 1728;
m_nonZeroFlops[22] = 162;
m_matrixKernels[22] = dgemm_m8_n27_k4_ldA56_ldB8_ldC8_beta0_pfsigonly;
m_hardwareFlops[22] = 1728;
m_nonZeroFlops[23] = 48;
m_matrixKernels[23] = dgemm_m8_n27_k27_ldA8_ldB27_ldC8_beta1_pfsigonly;
m_hardwareFlops[23] = 11664;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 37044;
m_matrixKernels[0] = dgemm_m88_n27_k56_ldA88_ldB88_ldC88_beta0_pfsigonly;
m_hardwareFlops[0] = 266112;
m_nonZeroFlops[1] = 83916;
m_matrixKernels[1] = dgemm_m88_n27_k56_ldA88_ldB88_ldC88_beta0_pfsigonly;
m_hardwareFlops[1] = 266112;
m_nonZeroFlops[2] = 90720;
m_matrixKernels[2] = dgemm_m88_n27_k56_ldA88_ldB88_ldC88_beta0_pfsigonly;
m_hardwareFlops[2] = 266112;
m_nonZeroFlops[3] = 4032;
m_matrixKernels[3] = dgemm_m88_n27_k27_ldA88_ldB27_ldC88_beta1_pfsigonly;
m_hardwareFlops[3] = 128304;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 18144;
m_hardwareFlops[0] = 399168;
m_matrixKernels[0] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
m_nonZeroFlops[1] = 58968;
m_hardwareFlops[1] = 399168;
m_matrixKernels[1] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
m_nonZeroFlops[2] = 185976;
m_hardwareFlops[2] = 399168;
m_matrixKernels[2] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
m_nonZeroFlops[3] = 185976;
m_hardwareFlops[3] = 399168;
m_matrixKernels[3] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
m_nonZeroFlops[4] = 57996;
m_hardwareFlops[4] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[4] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 57996;
m_hardwareFlops[5] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[5] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 18144;
m_hardwareFlops[6] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[6] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 72576;
m_hardwareFlops[7] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[7] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 231984;
m_hardwareFlops[8] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[8] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 231984;
m_hardwareFlops[9] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[9] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 303804;
m_hardwareFlops[10] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[10] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 269406;
m_hardwareFlops[11] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[11] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 267570;
m_hardwareFlops[12] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[12] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 267570;
m_hardwareFlops[13] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[13] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 269406;
m_hardwareFlops[14] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[14] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 303804;
m_hardwareFlops[15] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[15] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 72576;
m_hardwareFlops[16] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[16] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 231984;
m_hardwareFlops[17] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[17] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 231984;
m_hardwareFlops[18] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[18] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 328860;
m_hardwareFlops[19] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[19] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 58968;
m_hardwareFlops[20] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[20] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 328860;
m_hardwareFlops[21] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[21] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 207630;
m_hardwareFlops[22] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[22] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 327726;
m_hardwareFlops[23] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[23] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 305532;
m_hardwareFlops[24] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[24] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 327726;
m_hardwareFlops[25] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[25] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 207630;
m_hardwareFlops[26] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[26] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 305532;
m_hardwareFlops[27] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[27] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 303804;
m_hardwareFlops[28] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[28] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 269406;
m_hardwareFlops[29] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[29] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 267570;
m_hardwareFlops[30] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[30] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 207630;
m_hardwareFlops[31] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[31] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 327726;
m_hardwareFlops[32] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[32] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 305532;
m_hardwareFlops[33] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[33] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 234468;
m_hardwareFlops[34] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[34] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 379188;
m_hardwareFlops[35] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[35] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 257364;
m_hardwareFlops[36] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[36] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 185976;
m_hardwareFlops[37] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[37] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 354024;
m_hardwareFlops[38] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[38] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 354024;
m_hardwareFlops[39] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[39] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 267570;
m_hardwareFlops[40] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[40] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 269406;
m_hardwareFlops[41] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[41] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 303804;
m_hardwareFlops[42] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[42] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 327726;
m_hardwareFlops[43] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[43] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 207630;
m_hardwareFlops[44] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[44] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 305532;
m_hardwareFlops[45] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[45] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 185976;
m_hardwareFlops[46] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[46] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 354024;
m_hardwareFlops[47] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[47] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 354024;
m_hardwareFlops[48] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[48] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 234468;
m_hardwareFlops[49] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[49] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 257364;
m_hardwareFlops[50] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[50] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 379188;
m_hardwareFlops[51] = 399168;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[51] = dgemm_m88_n27_k84_ldA88_ldB88_ldC88_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 13608;
m_hardwareFlops[52] = 128304;
m_matrixKernels[52] = dgemm_m88_n27_k27_ldA88_ldB27_ldC88_beta1_pfsigonly;
m_nonZeroFlops[53] = 13608;
m_hardwareFlops[53] = 128304;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = dgemm_m88_n27_k27_ldA88_ldB27_ldC88_beta1_AL2jpst_BL2viaC;
#else
m_matrixKernels[53] = dgemm_m88_n27_k27_ldA88_ldB27_ldC88_beta1_pfsigonly;
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
m_matrixKernels[0] = dgemm_m88_n27_k120_ldA88_ldB120_ldC88_beta0_pfsigonly;
m_hardwareFlops[0] = 570240;
m_nonZeroFlops[1] = 174312;
m_matrixKernels[1] = dgemm_m88_n27_k120_ldA88_ldB120_ldC88_beta0_pfsigonly;
m_hardwareFlops[1] = 570240;
m_nonZeroFlops[2] = 185652;
m_matrixKernels[2] = dgemm_m88_n27_k120_ldA88_ldB120_ldC88_beta0_pfsigonly;
m_hardwareFlops[2] = 570240;
m_nonZeroFlops[3] = 4032;
m_matrixKernels[3] = dgemm_m88_n27_k27_ldA88_ldB27_ldC88_beta1_pfsigonly;
m_hardwareFlops[3] = 128304;
m_nonZeroFlops[4] = 37044;
m_matrixKernels[4] = dgemm_m56_n27_k84_ldA88_ldB88_ldC56_beta0_pfsigonly;
m_hardwareFlops[4] = 254016;
m_nonZeroFlops[5] = 83916;
m_matrixKernels[5] = dgemm_m56_n27_k84_ldA88_ldB88_ldC56_beta0_pfsigonly;
m_hardwareFlops[5] = 254016;
m_nonZeroFlops[6] = 90720;
m_matrixKernels[6] = dgemm_m56_n27_k84_ldA88_ldB88_ldC56_beta0_pfsigonly;
m_hardwareFlops[6] = 254016;
m_nonZeroFlops[7] = 2688;
m_matrixKernels[7] = dgemm_m56_n27_k27_ldA56_ldB27_ldC56_beta1_pfsigonly;
m_hardwareFlops[7] = 81648;
m_nonZeroFlops[8] = 15876;
m_matrixKernels[8] = dgemm_m40_n27_k56_ldA88_ldB56_ldC40_beta0_pfsigonly;
m_hardwareFlops[8] = 120960;
m_nonZeroFlops[9] = 36288;
m_matrixKernels[9] = dgemm_m40_n27_k56_ldA88_ldB56_ldC40_beta0_pfsigonly;
m_hardwareFlops[9] = 120960;
m_nonZeroFlops[10] = 40068;
m_matrixKernels[10] = dgemm_m40_n27_k56_ldA88_ldB56_ldC40_beta0_pfsigonly;
m_hardwareFlops[10] = 120960;
m_nonZeroFlops[11] = 1680;
m_matrixKernels[11] = dgemm_m40_n27_k27_ldA40_ldB27_ldC40_beta1_pfsigonly;
m_hardwareFlops[11] = 58320;
m_nonZeroFlops[12] = 5832;
m_matrixKernels[12] = dgemm_m24_n27_k35_ldA88_ldB40_ldC24_beta0_pfsigonly;
m_hardwareFlops[12] = 45360;
m_nonZeroFlops[13] = 13608;
m_matrixKernels[13] = dgemm_m24_n27_k35_ldA88_ldB40_ldC24_beta0_pfsigonly;
m_hardwareFlops[13] = 45360;
m_nonZeroFlops[14] = 15498;
m_matrixKernels[14] = dgemm_m24_n27_k35_ldA88_ldB40_ldC24_beta0_pfsigonly;
m_hardwareFlops[14] = 45360;
m_nonZeroFlops[15] = 960;
m_matrixKernels[15] = dgemm_m24_n27_k27_ldA24_ldB27_ldC24_beta1_pfsigonly;
m_hardwareFlops[15] = 34992;
m_nonZeroFlops[16] = 1782;
m_matrixKernels[16] = dgemm_m16_n27_k20_ldA88_ldB24_ldC16_beta0_pfsigonly;
m_hardwareFlops[16] = 17280;
m_nonZeroFlops[17] = 4158;
m_matrixKernels[17] = dgemm_m16_n27_k20_ldA88_ldB24_ldC16_beta0_pfsigonly;
m_hardwareFlops[17] = 17280;
m_nonZeroFlops[18] = 4968;
m_matrixKernels[18] = dgemm_m16_n27_k20_ldA88_ldB24_ldC16_beta0_pfsigonly;
m_hardwareFlops[18] = 17280;
m_nonZeroFlops[19] = 480;
m_matrixKernels[19] = dgemm_m16_n27_k27_ldA16_ldB27_ldC16_beta1_pfsigonly;
m_hardwareFlops[19] = 23328;
m_nonZeroFlops[20] = 378;
m_matrixKernels[20] = dgemm_m8_n27_k10_ldA88_ldB16_ldC8_beta0_pfsigonly;
m_hardwareFlops[20] = 4320;
m_nonZeroFlops[21] = 918;
m_matrixKernels[21] = dgemm_m8_n27_k10_ldA88_ldB16_ldC8_beta0_pfsigonly;
m_hardwareFlops[21] = 4320;
m_nonZeroFlops[22] = 1188;
m_matrixKernels[22] = dgemm_m8_n27_k10_ldA88_ldB16_ldC8_beta0_pfsigonly;
m_hardwareFlops[22] = 4320;
m_nonZeroFlops[23] = 192;
m_matrixKernels[23] = dgemm_m8_n27_k27_ldA8_ldB27_ldC8_beta1_pfsigonly;
m_hardwareFlops[23] = 11664;
m_nonZeroFlops[24] = 54;
m_matrixKernels[24] = dgemm_m8_n27_k4_ldA88_ldB8_ldC8_beta0_pfsigonly;
m_hardwareFlops[24] = 1728;
m_nonZeroFlops[25] = 108;
m_matrixKernels[25] = dgemm_m8_n27_k4_ldA88_ldB8_ldC8_beta0_pfsigonly;
m_hardwareFlops[25] = 1728;
m_nonZeroFlops[26] = 162;
m_matrixKernels[26] = dgemm_m8_n27_k4_ldA88_ldB8_ldC8_beta0_pfsigonly;
m_hardwareFlops[26] = 1728;
m_nonZeroFlops[27] = 48;
m_matrixKernels[27] = dgemm_m8_n27_k27_ldA8_ldB27_ldC8_beta1_pfsigonly;
m_hardwareFlops[27] = 11664;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 78084;
m_matrixKernels[0] = dgemm_m120_n27_k84_ldA120_ldB120_ldC120_beta0_pfsigonly;
m_hardwareFlops[0] = 544320;
m_nonZeroFlops[1] = 174312;
m_matrixKernels[1] = dgemm_m120_n27_k84_ldA120_ldB120_ldC120_beta0_pfsigonly;
m_hardwareFlops[1] = 544320;
m_nonZeroFlops[2] = 185652;
m_matrixKernels[2] = dgemm_m120_n27_k84_ldA120_ldB120_ldC120_beta0_pfsigonly;
m_hardwareFlops[2] = 544320;
m_nonZeroFlops[3] = 5760;
m_matrixKernels[3] = dgemm_m120_n27_k27_ldA120_ldB27_ldC120_beta1_pfsigonly;
m_hardwareFlops[3] = 174960;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 29160;
m_hardwareFlops[0] = 777600;
m_matrixKernels[0] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
m_nonZeroFlops[1] = 106272;
m_hardwareFlops[1] = 777600;
m_matrixKernels[1] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
m_nonZeroFlops[2] = 377784;
m_hardwareFlops[2] = 777600;
m_matrixKernels[2] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
m_nonZeroFlops[3] = 377784;
m_hardwareFlops[3] = 777600;
m_matrixKernels[3] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
m_nonZeroFlops[4] = 104544;
m_hardwareFlops[4] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[4] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 104544;
m_hardwareFlops[5] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[5] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 29160;
m_hardwareFlops[6] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[6] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 131220;
m_hardwareFlops[7] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[7] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 470448;
m_hardwareFlops[8] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[8] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 470448;
m_hardwareFlops[9] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[9] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 621486;
m_hardwareFlops[10] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[10] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 545940;
m_hardwareFlops[11] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[11] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 544860;
m_hardwareFlops[12] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[12] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 544860;
m_hardwareFlops[13] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[13] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 545940;
m_hardwareFlops[14] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[14] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 621486;
m_hardwareFlops[15] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[15] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 131220;
m_hardwareFlops[16] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[16] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 470448;
m_hardwareFlops[17] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[17] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 470448;
m_hardwareFlops[18] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[18] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 674028;
m_hardwareFlops[19] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[19] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 106272;
m_hardwareFlops[20] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[20] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 674028;
m_hardwareFlops[21] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[21] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 422280;
m_hardwareFlops[22] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[22] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 673164;
m_hardwareFlops[23] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[23] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 621432;
m_hardwareFlops[24] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[24] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 673164;
m_hardwareFlops[25] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[25] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 422280;
m_hardwareFlops[26] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[26] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 621432;
m_hardwareFlops[27] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[27] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 621486;
m_hardwareFlops[28] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[28] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 545940;
m_hardwareFlops[29] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[29] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 544860;
m_hardwareFlops[30] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[30] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 422280;
m_hardwareFlops[31] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[31] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 673164;
m_hardwareFlops[32] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[32] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 621432;
m_hardwareFlops[33] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[33] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 476064;
m_hardwareFlops[34] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[34] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 775332;
m_hardwareFlops[35] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[35] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 525852;
m_hardwareFlops[36] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[36] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 377784;
m_hardwareFlops[37] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[37] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 724734;
m_hardwareFlops[38] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[38] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 724734;
m_hardwareFlops[39] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[39] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 544860;
m_hardwareFlops[40] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[40] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 545940;
m_hardwareFlops[41] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[41] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 621486;
m_hardwareFlops[42] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[42] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 673164;
m_hardwareFlops[43] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[43] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 422280;
m_hardwareFlops[44] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[44] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 621432;
m_hardwareFlops[45] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[45] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 377784;
m_hardwareFlops[46] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[46] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 724734;
m_hardwareFlops[47] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[47] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 724734;
m_hardwareFlops[48] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[48] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 476064;
m_hardwareFlops[49] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[49] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 525852;
m_hardwareFlops[50] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[50] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 775332;
m_hardwareFlops[51] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_curAL2_BL2viaC;
#else
m_matrixKernels[51] = dgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 19440;
m_hardwareFlops[52] = 174960;
m_matrixKernels[52] = dgemm_m120_n27_k27_ldA120_ldB27_ldC120_beta1_pfsigonly;
m_nonZeroFlops[53] = 19440;
m_hardwareFlops[53] = 174960;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = dgemm_m120_n27_k27_ldA120_ldB27_ldC120_beta1_AL2jpst_BL2viaC;
#else
m_matrixKernels[53] = dgemm_m120_n27_k27_ldA120_ldB27_ldC120_beta1_pfsigonly;
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

