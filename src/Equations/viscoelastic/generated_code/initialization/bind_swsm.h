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
// @date 2015-07-13 15:29:43.279785
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
m_nonZeroFlops[0] = 54;
m_matrixKernels[0] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
m_hardwareFlops[0] = 864;
m_nonZeroFlops[1] = 108;
m_matrixKernels[1] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
m_hardwareFlops[1] = 864;
m_nonZeroFlops[2] = 162;
m_matrixKernels[2] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
m_hardwareFlops[2] = 864;
m_nonZeroFlops[3] = 48;
m_matrixKernels[3] = sgemm_m4_n27_k27_ldA4_ldB27_ldC4_beta1_pfsigonly;
m_hardwareFlops[3] = 5832;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 54;
m_matrixKernels[0] = sgemm_m4_n27_k1_ldA4_ldB4_ldC4_beta0_pfsigonly;
m_hardwareFlops[0] = 216;
m_nonZeroFlops[1] = 108;
m_matrixKernels[1] = sgemm_m4_n27_k1_ldA4_ldB4_ldC4_beta0_pfsigonly;
m_hardwareFlops[1] = 216;
m_nonZeroFlops[2] = 162;
m_matrixKernels[2] = sgemm_m4_n27_k1_ldA4_ldB4_ldC4_beta0_pfsigonly;
m_hardwareFlops[2] = 216;
m_nonZeroFlops[3] = 192;
m_matrixKernels[3] = sgemm_m4_n27_k27_ldA4_ldB27_ldC4_beta1_pfsigonly;
m_hardwareFlops[3] = 5832;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 324;
m_hardwareFlops[0] = 864;
m_matrixKernels[0] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
m_nonZeroFlops[1] = 432;
m_hardwareFlops[1] = 864;
m_matrixKernels[1] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
m_nonZeroFlops[2] = 540;
m_hardwareFlops[2] = 864;
m_matrixKernels[2] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
m_nonZeroFlops[3] = 540;
m_hardwareFlops[3] = 864;
m_matrixKernels[3] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
m_nonZeroFlops[4] = 432;
m_hardwareFlops[4] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[4] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 432;
m_hardwareFlops[5] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[5] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 324;
m_hardwareFlops[6] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[6] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 486;
m_hardwareFlops[7] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[7] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 648;
m_hardwareFlops[8] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[8] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 648;
m_hardwareFlops[9] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[9] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 756;
m_hardwareFlops[10] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[10] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 702;
m_hardwareFlops[11] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[11] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 702;
m_hardwareFlops[12] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[12] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 702;
m_hardwareFlops[13] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[13] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 702;
m_hardwareFlops[14] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[14] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 756;
m_hardwareFlops[15] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[15] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 486;
m_hardwareFlops[16] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[16] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 648;
m_hardwareFlops[17] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[17] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 648;
m_hardwareFlops[18] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[18] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 756;
m_hardwareFlops[19] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[19] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 432;
m_hardwareFlops[20] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[20] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 756;
m_hardwareFlops[21] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[21] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 594;
m_hardwareFlops[22] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[22] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 756;
m_hardwareFlops[23] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[23] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 756;
m_hardwareFlops[24] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[24] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 756;
m_hardwareFlops[25] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[25] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 594;
m_hardwareFlops[26] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[26] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 756;
m_hardwareFlops[27] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[27] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 756;
m_hardwareFlops[28] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[28] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 702;
m_hardwareFlops[29] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[29] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 702;
m_hardwareFlops[30] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[30] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 594;
m_hardwareFlops[31] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[31] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 756;
m_hardwareFlops[32] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[32] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 756;
m_hardwareFlops[33] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[33] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 648;
m_hardwareFlops[34] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[34] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 864;
m_hardwareFlops[35] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[35] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 648;
m_hardwareFlops[36] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[36] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 540;
m_hardwareFlops[37] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[37] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 810;
m_hardwareFlops[38] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[38] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 810;
m_hardwareFlops[39] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[39] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 702;
m_hardwareFlops[40] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[40] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 702;
m_hardwareFlops[41] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[41] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 756;
m_hardwareFlops[42] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[42] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 756;
m_hardwareFlops[43] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[43] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 594;
m_hardwareFlops[44] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[44] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 756;
m_hardwareFlops[45] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[45] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 540;
m_hardwareFlops[46] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[46] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 810;
m_hardwareFlops[47] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[47] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 810;
m_hardwareFlops[48] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[48] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 648;
m_hardwareFlops[49] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[49] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 648;
m_hardwareFlops[50] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[50] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 864;
m_hardwareFlops[51] = 864;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_BL2viaC;
#else
m_matrixKernels[51] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 648;
m_hardwareFlops[52] = 5832;
m_matrixKernels[52] = sgemm_m4_n27_k27_ldA4_ldB27_ldC4_beta1_pfsigonly;
m_nonZeroFlops[53] = 648;
m_hardwareFlops[53] = 5832;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = sgemm_m4_n27_k27_ldA4_ldB27_ldC4_beta1_pfsigonly;
#else
m_matrixKernels[53] = sgemm_m4_n27_k27_ldA4_ldB27_ldC4_beta1_pfsigonly;
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
m_matrixKernels[0] = sgemm_m4_n27_k10_ldA4_ldB12_ldC4_beta0_pfsigonly;
m_hardwareFlops[0] = 2160;
m_nonZeroFlops[1] = 918;
m_matrixKernels[1] = sgemm_m4_n27_k10_ldA4_ldB12_ldC4_beta0_pfsigonly;
m_hardwareFlops[1] = 2160;
m_nonZeroFlops[2] = 1188;
m_matrixKernels[2] = sgemm_m4_n27_k10_ldA4_ldB12_ldC4_beta0_pfsigonly;
m_hardwareFlops[2] = 2160;
m_nonZeroFlops[3] = 192;
m_matrixKernels[3] = sgemm_m4_n27_k27_ldA4_ldB27_ldC4_beta1_pfsigonly;
m_hardwareFlops[3] = 5832;
m_nonZeroFlops[4] = 54;
m_matrixKernels[4] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
m_hardwareFlops[4] = 864;
m_nonZeroFlops[5] = 108;
m_matrixKernels[5] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
m_hardwareFlops[5] = 864;
m_nonZeroFlops[6] = 162;
m_matrixKernels[6] = sgemm_m4_n27_k4_ldA4_ldB4_ldC4_beta0_pfsigonly;
m_hardwareFlops[6] = 864;
m_nonZeroFlops[7] = 48;
m_matrixKernels[7] = sgemm_m4_n27_k27_ldA4_ldB27_ldC4_beta1_pfsigonly;
m_hardwareFlops[7] = 5832;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 378;
m_matrixKernels[0] = sgemm_m12_n27_k4_ldA12_ldB12_ldC12_beta0_pfsigonly;
m_hardwareFlops[0] = 2592;
m_nonZeroFlops[1] = 918;
m_matrixKernels[1] = sgemm_m12_n27_k4_ldA12_ldB12_ldC12_beta0_pfsigonly;
m_hardwareFlops[1] = 2592;
m_nonZeroFlops[2] = 1188;
m_matrixKernels[2] = sgemm_m12_n27_k4_ldA12_ldB12_ldC12_beta0_pfsigonly;
m_hardwareFlops[2] = 2592;
m_nonZeroFlops[3] = 480;
m_matrixKernels[3] = sgemm_m12_n27_k27_ldA12_ldB27_ldC12_beta1_pfsigonly;
m_hardwareFlops[3] = 17496;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 1080;
m_hardwareFlops[0] = 6480;
m_matrixKernels[0] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
m_nonZeroFlops[1] = 1836;
m_hardwareFlops[1] = 6480;
m_matrixKernels[1] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
m_nonZeroFlops[2] = 2916;
m_hardwareFlops[2] = 6480;
m_matrixKernels[2] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
m_nonZeroFlops[3] = 2916;
m_hardwareFlops[3] = 6480;
m_matrixKernels[3] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
m_nonZeroFlops[4] = 1836;
m_hardwareFlops[4] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[4] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 1836;
m_hardwareFlops[5] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[5] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 1080;
m_hardwareFlops[6] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[6] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 2160;
m_hardwareFlops[7] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[7] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 3672;
m_hardwareFlops[8] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[8] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 3672;
m_hardwareFlops[9] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[9] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 4374;
m_hardwareFlops[10] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[10] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 4104;
m_hardwareFlops[11] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[11] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 4104;
m_hardwareFlops[12] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[12] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 4104;
m_hardwareFlops[13] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[13] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 4104;
m_hardwareFlops[14] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[14] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 4374;
m_hardwareFlops[15] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[15] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 2160;
m_hardwareFlops[16] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[16] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 3672;
m_hardwareFlops[17] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[17] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 3672;
m_hardwareFlops[18] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[18] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 4536;
m_hardwareFlops[19] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[19] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 1836;
m_hardwareFlops[20] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[20] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 4536;
m_hardwareFlops[21] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[21] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 3294;
m_hardwareFlops[22] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[22] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 4536;
m_hardwareFlops[23] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[23] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 4536;
m_hardwareFlops[24] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[24] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 4536;
m_hardwareFlops[25] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[25] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 3294;
m_hardwareFlops[26] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[26] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 4536;
m_hardwareFlops[27] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[27] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 4374;
m_hardwareFlops[28] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[28] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 4104;
m_hardwareFlops[29] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[29] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 4104;
m_hardwareFlops[30] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[30] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 3294;
m_hardwareFlops[31] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[31] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 4536;
m_hardwareFlops[32] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[32] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 4536;
m_hardwareFlops[33] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[33] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 3672;
m_hardwareFlops[34] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[34] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 5292;
m_hardwareFlops[35] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[35] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 3780;
m_hardwareFlops[36] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[36] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 2916;
m_hardwareFlops[37] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[37] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 5022;
m_hardwareFlops[38] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[38] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 5022;
m_hardwareFlops[39] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[39] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 4104;
m_hardwareFlops[40] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[40] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 4104;
m_hardwareFlops[41] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[41] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 4374;
m_hardwareFlops[42] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[42] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 4536;
m_hardwareFlops[43] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[43] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 3294;
m_hardwareFlops[44] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[44] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 4536;
m_hardwareFlops[45] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[45] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 2916;
m_hardwareFlops[46] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[46] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 5022;
m_hardwareFlops[47] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[47] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 5022;
m_hardwareFlops[48] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[48] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 3672;
m_hardwareFlops[49] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[49] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 3780;
m_hardwareFlops[50] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[50] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 5292;
m_hardwareFlops[51] = 6480;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_BL2viaC;
#else
m_matrixKernels[51] = sgemm_m12_n27_k10_ldA12_ldB12_ldC12_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 1620;
m_hardwareFlops[52] = 17496;
m_matrixKernels[52] = sgemm_m12_n27_k27_ldA12_ldB27_ldC12_beta1_pfsigonly;
m_nonZeroFlops[53] = 1620;
m_hardwareFlops[53] = 17496;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = sgemm_m12_n27_k27_ldA12_ldB27_ldC12_beta1_pfsigonly;
#else
m_matrixKernels[53] = sgemm_m12_n27_k27_ldA12_ldB27_ldC12_beta1_pfsigonly;
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
m_matrixKernels[0] = sgemm_m12_n27_k20_ldA12_ldB20_ldC12_beta0_pfsigonly;
m_hardwareFlops[0] = 12960;
m_nonZeroFlops[1] = 4158;
m_matrixKernels[1] = sgemm_m12_n27_k20_ldA12_ldB20_ldC12_beta0_pfsigonly;
m_hardwareFlops[1] = 12960;
m_nonZeroFlops[2] = 4968;
m_matrixKernels[2] = sgemm_m12_n27_k20_ldA12_ldB20_ldC12_beta0_pfsigonly;
m_hardwareFlops[2] = 12960;
m_nonZeroFlops[3] = 480;
m_matrixKernels[3] = sgemm_m12_n27_k27_ldA12_ldB27_ldC12_beta1_pfsigonly;
m_hardwareFlops[3] = 17496;
m_nonZeroFlops[4] = 378;
m_matrixKernels[4] = sgemm_m4_n27_k10_ldA12_ldB12_ldC4_beta0_pfsigonly;
m_hardwareFlops[4] = 2160;
m_nonZeroFlops[5] = 918;
m_matrixKernels[5] = sgemm_m4_n27_k10_ldA12_ldB12_ldC4_beta0_pfsigonly;
m_hardwareFlops[5] = 2160;
m_nonZeroFlops[6] = 1188;
m_matrixKernels[6] = sgemm_m4_n27_k10_ldA12_ldB12_ldC4_beta0_pfsigonly;
m_hardwareFlops[6] = 2160;
m_nonZeroFlops[7] = 192;
m_matrixKernels[7] = sgemm_m4_n27_k27_ldA4_ldB27_ldC4_beta1_pfsigonly;
m_hardwareFlops[7] = 5832;
m_nonZeroFlops[8] = 54;
m_matrixKernels[8] = sgemm_m4_n27_k4_ldA12_ldB4_ldC4_beta0_pfsigonly;
m_hardwareFlops[8] = 864;
m_nonZeroFlops[9] = 108;
m_matrixKernels[9] = sgemm_m4_n27_k4_ldA12_ldB4_ldC4_beta0_pfsigonly;
m_hardwareFlops[9] = 864;
m_nonZeroFlops[10] = 162;
m_matrixKernels[10] = sgemm_m4_n27_k4_ldA12_ldB4_ldC4_beta0_pfsigonly;
m_hardwareFlops[10] = 864;
m_nonZeroFlops[11] = 48;
m_matrixKernels[11] = sgemm_m4_n27_k27_ldA4_ldB27_ldC4_beta1_pfsigonly;
m_hardwareFlops[11] = 5832;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 1782;
m_matrixKernels[0] = sgemm_m20_n27_k10_ldA20_ldB20_ldC20_beta0_pfsigonly;
m_hardwareFlops[0] = 10800;
m_nonZeroFlops[1] = 4158;
m_matrixKernels[1] = sgemm_m20_n27_k10_ldA20_ldB20_ldC20_beta0_pfsigonly;
m_hardwareFlops[1] = 10800;
m_nonZeroFlops[2] = 4968;
m_matrixKernels[2] = sgemm_m20_n27_k10_ldA20_ldB20_ldC20_beta0_pfsigonly;
m_hardwareFlops[2] = 10800;
m_nonZeroFlops[3] = 960;
m_matrixKernels[3] = sgemm_m20_n27_k27_ldA20_ldB27_ldC20_beta1_pfsigonly;
m_hardwareFlops[3] = 29160;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 2700;
m_hardwareFlops[0] = 21600;
m_matrixKernels[0] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
m_nonZeroFlops[1] = 5616;
m_hardwareFlops[1] = 21600;
m_matrixKernels[1] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
m_nonZeroFlops[2] = 11016;
m_hardwareFlops[2] = 21600;
m_matrixKernels[2] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
m_nonZeroFlops[3] = 11016;
m_hardwareFlops[3] = 21600;
m_matrixKernels[3] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
m_nonZeroFlops[4] = 5616;
m_hardwareFlops[4] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[4] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 5616;
m_hardwareFlops[5] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[5] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 2700;
m_hardwareFlops[6] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[6] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 6750;
m_hardwareFlops[7] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[7] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 14040;
m_hardwareFlops[8] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[8] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 14040;
m_hardwareFlops[9] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[9] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 17280;
m_hardwareFlops[10] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[10] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 15444;
m_hardwareFlops[11] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[11] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 15930;
m_hardwareFlops[12] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[12] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 15930;
m_hardwareFlops[13] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[13] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 15444;
m_hardwareFlops[14] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[14] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 17280;
m_hardwareFlops[15] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[15] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 6750;
m_hardwareFlops[16] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[16] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 14040;
m_hardwareFlops[17] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[17] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 14040;
m_hardwareFlops[18] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[18] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 18252;
m_hardwareFlops[19] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[19] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 5616;
m_hardwareFlops[20] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[20] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 18252;
m_hardwareFlops[21] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[21] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 12528;
m_hardwareFlops[22] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[22] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 18360;
m_hardwareFlops[23] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[23] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 17550;
m_hardwareFlops[24] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[24] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 18360;
m_hardwareFlops[25] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[25] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 12528;
m_hardwareFlops[26] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[26] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 17550;
m_hardwareFlops[27] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[27] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 17280;
m_hardwareFlops[28] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[28] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 15444;
m_hardwareFlops[29] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[29] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 15930;
m_hardwareFlops[30] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[30] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 12528;
m_hardwareFlops[31] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[31] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 18360;
m_hardwareFlops[32] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[32] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 17550;
m_hardwareFlops[33] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[33] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 13932;
m_hardwareFlops[34] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[34] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 21276;
m_hardwareFlops[35] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[35] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 14796;
m_hardwareFlops[36] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[36] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 11016;
m_hardwareFlops[37] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[37] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 20034;
m_hardwareFlops[38] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[38] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 20034;
m_hardwareFlops[39] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[39] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 15930;
m_hardwareFlops[40] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[40] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 15444;
m_hardwareFlops[41] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[41] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 17280;
m_hardwareFlops[42] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[42] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 18360;
m_hardwareFlops[43] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[43] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 12528;
m_hardwareFlops[44] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[44] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 17550;
m_hardwareFlops[45] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[45] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 11016;
m_hardwareFlops[46] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[46] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 20034;
m_hardwareFlops[47] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[47] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 20034;
m_hardwareFlops[48] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[48] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 13932;
m_hardwareFlops[49] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[49] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 14796;
m_hardwareFlops[50] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[50] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 21276;
m_hardwareFlops[51] = 21600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_BL2viaC;
#else
m_matrixKernels[51] = sgemm_m20_n27_k20_ldA20_ldB20_ldC20_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 3240;
m_hardwareFlops[52] = 29160;
m_matrixKernels[52] = sgemm_m20_n27_k27_ldA20_ldB27_ldC20_beta1_pfsigonly;
m_nonZeroFlops[53] = 3240;
m_hardwareFlops[53] = 29160;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = sgemm_m20_n27_k27_ldA20_ldB27_ldC20_beta1_pfsigonly;
#else
m_matrixKernels[53] = sgemm_m20_n27_k27_ldA20_ldB27_ldC20_beta1_pfsigonly;
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
m_matrixKernels[0] = sgemm_m20_n27_k35_ldA20_ldB36_ldC20_beta0_pfsigonly;
m_hardwareFlops[0] = 37800;
m_nonZeroFlops[1] = 13608;
m_matrixKernels[1] = sgemm_m20_n27_k35_ldA20_ldB36_ldC20_beta0_pfsigonly;
m_hardwareFlops[1] = 37800;
m_nonZeroFlops[2] = 15498;
m_matrixKernels[2] = sgemm_m20_n27_k35_ldA20_ldB36_ldC20_beta0_pfsigonly;
m_hardwareFlops[2] = 37800;
m_nonZeroFlops[3] = 960;
m_matrixKernels[3] = sgemm_m20_n27_k27_ldA20_ldB27_ldC20_beta1_pfsigonly;
m_hardwareFlops[3] = 29160;
m_nonZeroFlops[4] = 1782;
m_matrixKernels[4] = sgemm_m12_n27_k20_ldA20_ldB20_ldC12_beta0_pfsigonly;
m_hardwareFlops[4] = 12960;
m_nonZeroFlops[5] = 4158;
m_matrixKernels[5] = sgemm_m12_n27_k20_ldA20_ldB20_ldC12_beta0_pfsigonly;
m_hardwareFlops[5] = 12960;
m_nonZeroFlops[6] = 4968;
m_matrixKernels[6] = sgemm_m12_n27_k20_ldA20_ldB20_ldC12_beta0_pfsigonly;
m_hardwareFlops[6] = 12960;
m_nonZeroFlops[7] = 480;
m_matrixKernels[7] = sgemm_m12_n27_k27_ldA12_ldB27_ldC12_beta1_pfsigonly;
m_hardwareFlops[7] = 17496;
m_nonZeroFlops[8] = 378;
m_matrixKernels[8] = sgemm_m4_n27_k10_ldA20_ldB12_ldC4_beta0_pfsigonly;
m_hardwareFlops[8] = 2160;
m_nonZeroFlops[9] = 918;
m_matrixKernels[9] = sgemm_m4_n27_k10_ldA20_ldB12_ldC4_beta0_pfsigonly;
m_hardwareFlops[9] = 2160;
m_nonZeroFlops[10] = 1188;
m_matrixKernels[10] = sgemm_m4_n27_k10_ldA20_ldB12_ldC4_beta0_pfsigonly;
m_hardwareFlops[10] = 2160;
m_nonZeroFlops[11] = 192;
m_matrixKernels[11] = sgemm_m4_n27_k27_ldA4_ldB27_ldC4_beta1_pfsigonly;
m_hardwareFlops[11] = 5832;
m_nonZeroFlops[12] = 54;
m_matrixKernels[12] = sgemm_m4_n27_k4_ldA20_ldB4_ldC4_beta0_pfsigonly;
m_hardwareFlops[12] = 864;
m_nonZeroFlops[13] = 108;
m_matrixKernels[13] = sgemm_m4_n27_k4_ldA20_ldB4_ldC4_beta0_pfsigonly;
m_hardwareFlops[13] = 864;
m_nonZeroFlops[14] = 162;
m_matrixKernels[14] = sgemm_m4_n27_k4_ldA20_ldB4_ldC4_beta0_pfsigonly;
m_hardwareFlops[14] = 864;
m_nonZeroFlops[15] = 48;
m_matrixKernels[15] = sgemm_m4_n27_k27_ldA4_ldB27_ldC4_beta1_pfsigonly;
m_hardwareFlops[15] = 5832;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 5832;
m_matrixKernels[0] = sgemm_m36_n27_k20_ldA36_ldB36_ldC36_beta0_pfsigonly;
m_hardwareFlops[0] = 38880;
m_nonZeroFlops[1] = 13608;
m_matrixKernels[1] = sgemm_m36_n27_k20_ldA36_ldB36_ldC36_beta0_pfsigonly;
m_hardwareFlops[1] = 38880;
m_nonZeroFlops[2] = 15498;
m_matrixKernels[2] = sgemm_m36_n27_k20_ldA36_ldB36_ldC36_beta0_pfsigonly;
m_hardwareFlops[2] = 38880;
m_nonZeroFlops[3] = 1680;
m_matrixKernels[3] = sgemm_m36_n27_k27_ldA36_ldB27_ldC36_beta1_pfsigonly;
m_hardwareFlops[3] = 52488;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 5670;
m_hardwareFlops[0] = 68040;
m_matrixKernels[0] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
m_nonZeroFlops[1] = 13986;
m_hardwareFlops[1] = 68040;
m_matrixKernels[1] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
m_nonZeroFlops[2] = 32886;
m_hardwareFlops[2] = 68040;
m_matrixKernels[2] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
m_nonZeroFlops[3] = 32886;
m_hardwareFlops[3] = 68040;
m_matrixKernels[3] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
m_nonZeroFlops[4] = 13878;
m_hardwareFlops[4] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[4] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 13878;
m_hardwareFlops[5] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[5] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 5670;
m_hardwareFlops[6] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[6] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 17010;
m_hardwareFlops[7] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[7] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 41634;
m_hardwareFlops[8] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[8] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 41634;
m_hardwareFlops[9] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[9] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 52866;
m_hardwareFlops[10] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[10] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 46440;
m_hardwareFlops[11] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[11] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 47304;
m_hardwareFlops[12] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[12] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 47304;
m_hardwareFlops[13] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[13] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 46440;
m_hardwareFlops[14] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[14] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 52866;
m_hardwareFlops[15] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[15] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 17010;
m_hardwareFlops[16] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[16] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 41634;
m_hardwareFlops[17] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[17] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 41634;
m_hardwareFlops[18] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[18] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 56430;
m_hardwareFlops[19] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[19] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 13986;
m_hardwareFlops[20] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[20] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 56430;
m_hardwareFlops[21] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[21] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 36936;
m_hardwareFlops[22] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[22] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 56106;
m_hardwareFlops[23] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[23] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 52920;
m_hardwareFlops[24] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[24] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 56106;
m_hardwareFlops[25] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[25] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 36936;
m_hardwareFlops[26] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[26] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 52920;
m_hardwareFlops[27] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[27] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 52866;
m_hardwareFlops[28] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[28] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 46440;
m_hardwareFlops[29] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[29] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 47304;
m_hardwareFlops[30] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[30] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 36936;
m_hardwareFlops[31] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[31] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 56106;
m_hardwareFlops[32] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[32] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 52920;
m_hardwareFlops[33] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[33] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 41742;
m_hardwareFlops[34] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[34] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 65502;
m_hardwareFlops[35] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[35] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 44874;
m_hardwareFlops[36] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[36] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 32886;
m_hardwareFlops[37] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[37] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 61074;
m_hardwareFlops[38] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[38] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 61074;
m_hardwareFlops[39] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[39] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 47304;
m_hardwareFlops[40] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[40] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 46440;
m_hardwareFlops[41] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[41] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 52866;
m_hardwareFlops[42] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[42] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 56106;
m_hardwareFlops[43] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[43] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 36936;
m_hardwareFlops[44] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[44] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 52920;
m_hardwareFlops[45] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[45] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 32886;
m_hardwareFlops[46] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[46] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 61074;
m_hardwareFlops[47] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[47] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 61074;
m_hardwareFlops[48] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[48] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 41742;
m_hardwareFlops[49] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[49] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 44874;
m_hardwareFlops[50] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[50] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 65502;
m_hardwareFlops[51] = 68040;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_BL2viaC;
#else
m_matrixKernels[51] = sgemm_m36_n27_k35_ldA36_ldB36_ldC36_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 5670;
m_hardwareFlops[52] = 52488;
m_matrixKernels[52] = sgemm_m36_n27_k27_ldA36_ldB27_ldC36_beta1_pfsigonly;
m_nonZeroFlops[53] = 5670;
m_hardwareFlops[53] = 52488;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = sgemm_m36_n27_k27_ldA36_ldB27_ldC36_beta1_pfsigonly;
#else
m_matrixKernels[53] = sgemm_m36_n27_k27_ldA36_ldB27_ldC36_beta1_pfsigonly;
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
m_matrixKernels[0] = sgemm_m36_n27_k56_ldA36_ldB56_ldC36_beta0_pfsigonly;
m_hardwareFlops[0] = 108864;
m_nonZeroFlops[1] = 36288;
m_matrixKernels[1] = sgemm_m36_n27_k56_ldA36_ldB56_ldC36_beta0_pfsigonly;
m_hardwareFlops[1] = 108864;
m_nonZeroFlops[2] = 40068;
m_matrixKernels[2] = sgemm_m36_n27_k56_ldA36_ldB56_ldC36_beta0_pfsigonly;
m_hardwareFlops[2] = 108864;
m_nonZeroFlops[3] = 1680;
m_matrixKernels[3] = sgemm_m36_n27_k27_ldA36_ldB27_ldC36_beta1_pfsigonly;
m_hardwareFlops[3] = 52488;
m_nonZeroFlops[4] = 5832;
m_matrixKernels[4] = sgemm_m20_n27_k35_ldA36_ldB36_ldC20_beta0_pfsigonly;
m_hardwareFlops[4] = 37800;
m_nonZeroFlops[5] = 13608;
m_matrixKernels[5] = sgemm_m20_n27_k35_ldA36_ldB36_ldC20_beta0_pfsigonly;
m_hardwareFlops[5] = 37800;
m_nonZeroFlops[6] = 15498;
m_matrixKernels[6] = sgemm_m20_n27_k35_ldA36_ldB36_ldC20_beta0_pfsigonly;
m_hardwareFlops[6] = 37800;
m_nonZeroFlops[7] = 960;
m_matrixKernels[7] = sgemm_m20_n27_k27_ldA20_ldB27_ldC20_beta1_pfsigonly;
m_hardwareFlops[7] = 29160;
m_nonZeroFlops[8] = 1782;
m_matrixKernels[8] = sgemm_m12_n27_k20_ldA36_ldB20_ldC12_beta0_pfsigonly;
m_hardwareFlops[8] = 12960;
m_nonZeroFlops[9] = 4158;
m_matrixKernels[9] = sgemm_m12_n27_k20_ldA36_ldB20_ldC12_beta0_pfsigonly;
m_hardwareFlops[9] = 12960;
m_nonZeroFlops[10] = 4968;
m_matrixKernels[10] = sgemm_m12_n27_k20_ldA36_ldB20_ldC12_beta0_pfsigonly;
m_hardwareFlops[10] = 12960;
m_nonZeroFlops[11] = 480;
m_matrixKernels[11] = sgemm_m12_n27_k27_ldA12_ldB27_ldC12_beta1_pfsigonly;
m_hardwareFlops[11] = 17496;
m_nonZeroFlops[12] = 378;
m_matrixKernels[12] = sgemm_m4_n27_k10_ldA36_ldB12_ldC4_beta0_pfsigonly;
m_hardwareFlops[12] = 2160;
m_nonZeroFlops[13] = 918;
m_matrixKernels[13] = sgemm_m4_n27_k10_ldA36_ldB12_ldC4_beta0_pfsigonly;
m_hardwareFlops[13] = 2160;
m_nonZeroFlops[14] = 1188;
m_matrixKernels[14] = sgemm_m4_n27_k10_ldA36_ldB12_ldC4_beta0_pfsigonly;
m_hardwareFlops[14] = 2160;
m_nonZeroFlops[15] = 192;
m_matrixKernels[15] = sgemm_m4_n27_k27_ldA4_ldB27_ldC4_beta1_pfsigonly;
m_hardwareFlops[15] = 5832;
m_nonZeroFlops[16] = 54;
m_matrixKernels[16] = sgemm_m4_n27_k4_ldA36_ldB4_ldC4_beta0_pfsigonly;
m_hardwareFlops[16] = 864;
m_nonZeroFlops[17] = 108;
m_matrixKernels[17] = sgemm_m4_n27_k4_ldA36_ldB4_ldC4_beta0_pfsigonly;
m_hardwareFlops[17] = 864;
m_nonZeroFlops[18] = 162;
m_matrixKernels[18] = sgemm_m4_n27_k4_ldA36_ldB4_ldC4_beta0_pfsigonly;
m_hardwareFlops[18] = 864;
m_nonZeroFlops[19] = 48;
m_matrixKernels[19] = sgemm_m4_n27_k27_ldA4_ldB27_ldC4_beta1_pfsigonly;
m_hardwareFlops[19] = 5832;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 15876;
m_matrixKernels[0] = sgemm_m56_n27_k35_ldA56_ldB56_ldC56_beta0_pfsigonly;
m_hardwareFlops[0] = 105840;
m_nonZeroFlops[1] = 36288;
m_matrixKernels[1] = sgemm_m56_n27_k35_ldA56_ldB56_ldC56_beta0_pfsigonly;
m_hardwareFlops[1] = 105840;
m_nonZeroFlops[2] = 40068;
m_matrixKernels[2] = sgemm_m56_n27_k35_ldA56_ldB56_ldC56_beta0_pfsigonly;
m_hardwareFlops[2] = 105840;
m_nonZeroFlops[3] = 2688;
m_matrixKernels[3] = sgemm_m56_n27_k27_ldA56_ldB27_ldC56_beta1_pfsigonly;
m_hardwareFlops[3] = 81648;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 10584;
m_hardwareFlops[0] = 169344;
m_matrixKernels[0] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
m_nonZeroFlops[1] = 30240;
m_hardwareFlops[1] = 169344;
m_matrixKernels[1] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
m_nonZeroFlops[2] = 83160;
m_hardwareFlops[2] = 169344;
m_matrixKernels[2] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
m_nonZeroFlops[3] = 83160;
m_hardwareFlops[3] = 169344;
m_matrixKernels[3] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
m_nonZeroFlops[4] = 29808;
m_hardwareFlops[4] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[4] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 29808;
m_hardwareFlops[5] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[5] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 10584;
m_hardwareFlops[6] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[6] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 37044;
m_hardwareFlops[7] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[7] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 104328;
m_hardwareFlops[8] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[8] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 104328;
m_hardwareFlops[9] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[9] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 135972;
m_hardwareFlops[10] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[10] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 119556;
m_hardwareFlops[11] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[11] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 119556;
m_hardwareFlops[12] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[12] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 119556;
m_hardwareFlops[13] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[13] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 119556;
m_hardwareFlops[14] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[14] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 135972;
m_hardwareFlops[15] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[15] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 37044;
m_hardwareFlops[16] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[16] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 104328;
m_hardwareFlops[17] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[17] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 104328;
m_hardwareFlops[18] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[18] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 145692;
m_hardwareFlops[19] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[19] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 30240;
m_hardwareFlops[20] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[20] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 145692;
m_hardwareFlops[21] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[21] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 92988;
m_hardwareFlops[22] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[22] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 145044;
m_hardwareFlops[23] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[23] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 135486;
m_hardwareFlops[24] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[24] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 145044;
m_hardwareFlops[25] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[25] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 92988;
m_hardwareFlops[26] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[26] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 135486;
m_hardwareFlops[27] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[27] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 135972;
m_hardwareFlops[28] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[28] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 119556;
m_hardwareFlops[29] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[29] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 119556;
m_hardwareFlops[30] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[30] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 92988;
m_hardwareFlops[31] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[31] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 145044;
m_hardwareFlops[32] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[32] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 135486;
m_hardwareFlops[33] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[33] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 105084;
m_hardwareFlops[34] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[34] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 168696;
m_hardwareFlops[35] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[35] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 114480;
m_hardwareFlops[36] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[36] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 83160;
m_hardwareFlops[37] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[37] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 156870;
m_hardwareFlops[38] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[38] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 156870;
m_hardwareFlops[39] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[39] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 119556;
m_hardwareFlops[40] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[40] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 119556;
m_hardwareFlops[41] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[41] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 135972;
m_hardwareFlops[42] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[42] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 145044;
m_hardwareFlops[43] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[43] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 92988;
m_hardwareFlops[44] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[44] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 135486;
m_hardwareFlops[45] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[45] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 83160;
m_hardwareFlops[46] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[46] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 156870;
m_hardwareFlops[47] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[47] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 156870;
m_hardwareFlops[48] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[48] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 105084;
m_hardwareFlops[49] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[49] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 114480;
m_hardwareFlops[50] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[50] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 168696;
m_hardwareFlops[51] = 169344;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_BL2viaC;
#else
m_matrixKernels[51] = sgemm_m56_n27_k56_ldA56_ldB56_ldC56_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 9072;
m_hardwareFlops[52] = 81648;
m_matrixKernels[52] = sgemm_m56_n27_k27_ldA56_ldB27_ldC56_beta1_pfsigonly;
m_nonZeroFlops[53] = 9072;
m_hardwareFlops[53] = 81648;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = sgemm_m56_n27_k27_ldA56_ldB27_ldC56_beta1_pfsigonly;
#else
m_matrixKernels[53] = sgemm_m56_n27_k27_ldA56_ldB27_ldC56_beta1_pfsigonly;
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
m_matrixKernels[0] = sgemm_m56_n27_k84_ldA56_ldB84_ldC56_beta0_pfsigonly;
m_hardwareFlops[0] = 254016;
m_nonZeroFlops[1] = 83916;
m_matrixKernels[1] = sgemm_m56_n27_k84_ldA56_ldB84_ldC56_beta0_pfsigonly;
m_hardwareFlops[1] = 254016;
m_nonZeroFlops[2] = 90720;
m_matrixKernels[2] = sgemm_m56_n27_k84_ldA56_ldB84_ldC56_beta0_pfsigonly;
m_hardwareFlops[2] = 254016;
m_nonZeroFlops[3] = 2688;
m_matrixKernels[3] = sgemm_m56_n27_k27_ldA56_ldB27_ldC56_beta1_pfsigonly;
m_hardwareFlops[3] = 81648;
m_nonZeroFlops[4] = 15876;
m_matrixKernels[4] = sgemm_m36_n27_k56_ldA56_ldB56_ldC36_beta0_pfsigonly;
m_hardwareFlops[4] = 108864;
m_nonZeroFlops[5] = 36288;
m_matrixKernels[5] = sgemm_m36_n27_k56_ldA56_ldB56_ldC36_beta0_pfsigonly;
m_hardwareFlops[5] = 108864;
m_nonZeroFlops[6] = 40068;
m_matrixKernels[6] = sgemm_m36_n27_k56_ldA56_ldB56_ldC36_beta0_pfsigonly;
m_hardwareFlops[6] = 108864;
m_nonZeroFlops[7] = 1680;
m_matrixKernels[7] = sgemm_m36_n27_k27_ldA36_ldB27_ldC36_beta1_pfsigonly;
m_hardwareFlops[7] = 52488;
m_nonZeroFlops[8] = 5832;
m_matrixKernels[8] = sgemm_m20_n27_k35_ldA56_ldB36_ldC20_beta0_pfsigonly;
m_hardwareFlops[8] = 37800;
m_nonZeroFlops[9] = 13608;
m_matrixKernels[9] = sgemm_m20_n27_k35_ldA56_ldB36_ldC20_beta0_pfsigonly;
m_hardwareFlops[9] = 37800;
m_nonZeroFlops[10] = 15498;
m_matrixKernels[10] = sgemm_m20_n27_k35_ldA56_ldB36_ldC20_beta0_pfsigonly;
m_hardwareFlops[10] = 37800;
m_nonZeroFlops[11] = 960;
m_matrixKernels[11] = sgemm_m20_n27_k27_ldA20_ldB27_ldC20_beta1_pfsigonly;
m_hardwareFlops[11] = 29160;
m_nonZeroFlops[12] = 1782;
m_matrixKernels[12] = sgemm_m12_n27_k20_ldA56_ldB20_ldC12_beta0_pfsigonly;
m_hardwareFlops[12] = 12960;
m_nonZeroFlops[13] = 4158;
m_matrixKernels[13] = sgemm_m12_n27_k20_ldA56_ldB20_ldC12_beta0_pfsigonly;
m_hardwareFlops[13] = 12960;
m_nonZeroFlops[14] = 4968;
m_matrixKernels[14] = sgemm_m12_n27_k20_ldA56_ldB20_ldC12_beta0_pfsigonly;
m_hardwareFlops[14] = 12960;
m_nonZeroFlops[15] = 480;
m_matrixKernels[15] = sgemm_m12_n27_k27_ldA12_ldB27_ldC12_beta1_pfsigonly;
m_hardwareFlops[15] = 17496;
m_nonZeroFlops[16] = 378;
m_matrixKernels[16] = sgemm_m4_n27_k10_ldA56_ldB12_ldC4_beta0_pfsigonly;
m_hardwareFlops[16] = 2160;
m_nonZeroFlops[17] = 918;
m_matrixKernels[17] = sgemm_m4_n27_k10_ldA56_ldB12_ldC4_beta0_pfsigonly;
m_hardwareFlops[17] = 2160;
m_nonZeroFlops[18] = 1188;
m_matrixKernels[18] = sgemm_m4_n27_k10_ldA56_ldB12_ldC4_beta0_pfsigonly;
m_hardwareFlops[18] = 2160;
m_nonZeroFlops[19] = 192;
m_matrixKernels[19] = sgemm_m4_n27_k27_ldA4_ldB27_ldC4_beta1_pfsigonly;
m_hardwareFlops[19] = 5832;
m_nonZeroFlops[20] = 54;
m_matrixKernels[20] = sgemm_m4_n27_k4_ldA56_ldB4_ldC4_beta0_pfsigonly;
m_hardwareFlops[20] = 864;
m_nonZeroFlops[21] = 108;
m_matrixKernels[21] = sgemm_m4_n27_k4_ldA56_ldB4_ldC4_beta0_pfsigonly;
m_hardwareFlops[21] = 864;
m_nonZeroFlops[22] = 162;
m_matrixKernels[22] = sgemm_m4_n27_k4_ldA56_ldB4_ldC4_beta0_pfsigonly;
m_hardwareFlops[22] = 864;
m_nonZeroFlops[23] = 48;
m_matrixKernels[23] = sgemm_m4_n27_k27_ldA4_ldB27_ldC4_beta1_pfsigonly;
m_hardwareFlops[23] = 5832;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 37044;
m_matrixKernels[0] = sgemm_m84_n27_k56_ldA84_ldB84_ldC84_beta0_pfsigonly;
m_hardwareFlops[0] = 254016;
m_nonZeroFlops[1] = 83916;
m_matrixKernels[1] = sgemm_m84_n27_k56_ldA84_ldB84_ldC84_beta0_pfsigonly;
m_hardwareFlops[1] = 254016;
m_nonZeroFlops[2] = 90720;
m_matrixKernels[2] = sgemm_m84_n27_k56_ldA84_ldB84_ldC84_beta0_pfsigonly;
m_hardwareFlops[2] = 254016;
m_nonZeroFlops[3] = 4032;
m_matrixKernels[3] = sgemm_m84_n27_k27_ldA84_ldB27_ldC84_beta1_pfsigonly;
m_hardwareFlops[3] = 122472;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 18144;
m_hardwareFlops[0] = 381024;
m_matrixKernels[0] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
m_nonZeroFlops[1] = 58968;
m_hardwareFlops[1] = 381024;
m_matrixKernels[1] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
m_nonZeroFlops[2] = 185976;
m_hardwareFlops[2] = 381024;
m_matrixKernels[2] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
m_nonZeroFlops[3] = 185976;
m_hardwareFlops[3] = 381024;
m_matrixKernels[3] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
m_nonZeroFlops[4] = 57996;
m_hardwareFlops[4] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[4] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 57996;
m_hardwareFlops[5] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[5] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 18144;
m_hardwareFlops[6] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[6] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 72576;
m_hardwareFlops[7] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[7] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 231984;
m_hardwareFlops[8] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[8] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 231984;
m_hardwareFlops[9] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[9] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 303804;
m_hardwareFlops[10] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[10] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 269406;
m_hardwareFlops[11] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[11] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 267570;
m_hardwareFlops[12] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[12] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 267570;
m_hardwareFlops[13] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[13] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 269406;
m_hardwareFlops[14] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[14] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 303804;
m_hardwareFlops[15] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[15] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 72576;
m_hardwareFlops[16] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[16] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 231984;
m_hardwareFlops[17] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[17] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 231984;
m_hardwareFlops[18] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[18] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 328860;
m_hardwareFlops[19] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[19] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 58968;
m_hardwareFlops[20] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[20] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 328860;
m_hardwareFlops[21] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[21] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 207630;
m_hardwareFlops[22] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[22] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 327726;
m_hardwareFlops[23] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[23] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 305532;
m_hardwareFlops[24] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[24] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 327726;
m_hardwareFlops[25] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[25] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 207630;
m_hardwareFlops[26] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[26] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 305532;
m_hardwareFlops[27] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[27] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 303804;
m_hardwareFlops[28] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[28] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 269406;
m_hardwareFlops[29] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[29] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 267570;
m_hardwareFlops[30] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[30] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 207630;
m_hardwareFlops[31] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[31] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 327726;
m_hardwareFlops[32] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[32] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 305532;
m_hardwareFlops[33] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[33] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 234468;
m_hardwareFlops[34] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[34] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 379188;
m_hardwareFlops[35] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[35] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 257364;
m_hardwareFlops[36] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[36] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 185976;
m_hardwareFlops[37] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[37] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 354024;
m_hardwareFlops[38] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[38] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 354024;
m_hardwareFlops[39] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[39] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 267570;
m_hardwareFlops[40] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[40] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 269406;
m_hardwareFlops[41] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[41] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 303804;
m_hardwareFlops[42] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[42] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 327726;
m_hardwareFlops[43] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[43] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 207630;
m_hardwareFlops[44] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[44] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 305532;
m_hardwareFlops[45] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[45] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 185976;
m_hardwareFlops[46] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[46] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 354024;
m_hardwareFlops[47] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[47] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 354024;
m_hardwareFlops[48] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[48] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 234468;
m_hardwareFlops[49] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[49] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 257364;
m_hardwareFlops[50] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[50] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 379188;
m_hardwareFlops[51] = 381024;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_BL2viaC;
#else
m_matrixKernels[51] = sgemm_m84_n27_k84_ldA84_ldB84_ldC84_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 13608;
m_hardwareFlops[52] = 122472;
m_matrixKernels[52] = sgemm_m84_n27_k27_ldA84_ldB27_ldC84_beta1_pfsigonly;
m_nonZeroFlops[53] = 13608;
m_hardwareFlops[53] = 122472;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = sgemm_m84_n27_k27_ldA84_ldB27_ldC84_beta1_pfsigonly;
#else
m_matrixKernels[53] = sgemm_m84_n27_k27_ldA84_ldB27_ldC84_beta1_pfsigonly;
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
m_matrixKernels[0] = sgemm_m84_n27_k120_ldA84_ldB120_ldC84_beta0_pfsigonly;
m_hardwareFlops[0] = 544320;
m_nonZeroFlops[1] = 174312;
m_matrixKernels[1] = sgemm_m84_n27_k120_ldA84_ldB120_ldC84_beta0_pfsigonly;
m_hardwareFlops[1] = 544320;
m_nonZeroFlops[2] = 185652;
m_matrixKernels[2] = sgemm_m84_n27_k120_ldA84_ldB120_ldC84_beta0_pfsigonly;
m_hardwareFlops[2] = 544320;
m_nonZeroFlops[3] = 4032;
m_matrixKernels[3] = sgemm_m84_n27_k27_ldA84_ldB27_ldC84_beta1_pfsigonly;
m_hardwareFlops[3] = 122472;
m_nonZeroFlops[4] = 37044;
m_matrixKernels[4] = sgemm_m56_n27_k84_ldA84_ldB84_ldC56_beta0_pfsigonly;
m_hardwareFlops[4] = 254016;
m_nonZeroFlops[5] = 83916;
m_matrixKernels[5] = sgemm_m56_n27_k84_ldA84_ldB84_ldC56_beta0_pfsigonly;
m_hardwareFlops[5] = 254016;
m_nonZeroFlops[6] = 90720;
m_matrixKernels[6] = sgemm_m56_n27_k84_ldA84_ldB84_ldC56_beta0_pfsigonly;
m_hardwareFlops[6] = 254016;
m_nonZeroFlops[7] = 2688;
m_matrixKernels[7] = sgemm_m56_n27_k27_ldA56_ldB27_ldC56_beta1_pfsigonly;
m_hardwareFlops[7] = 81648;
m_nonZeroFlops[8] = 15876;
m_matrixKernels[8] = sgemm_m36_n27_k56_ldA84_ldB56_ldC36_beta0_pfsigonly;
m_hardwareFlops[8] = 108864;
m_nonZeroFlops[9] = 36288;
m_matrixKernels[9] = sgemm_m36_n27_k56_ldA84_ldB56_ldC36_beta0_pfsigonly;
m_hardwareFlops[9] = 108864;
m_nonZeroFlops[10] = 40068;
m_matrixKernels[10] = sgemm_m36_n27_k56_ldA84_ldB56_ldC36_beta0_pfsigonly;
m_hardwareFlops[10] = 108864;
m_nonZeroFlops[11] = 1680;
m_matrixKernels[11] = sgemm_m36_n27_k27_ldA36_ldB27_ldC36_beta1_pfsigonly;
m_hardwareFlops[11] = 52488;
m_nonZeroFlops[12] = 5832;
m_matrixKernels[12] = sgemm_m20_n27_k35_ldA84_ldB36_ldC20_beta0_pfsigonly;
m_hardwareFlops[12] = 37800;
m_nonZeroFlops[13] = 13608;
m_matrixKernels[13] = sgemm_m20_n27_k35_ldA84_ldB36_ldC20_beta0_pfsigonly;
m_hardwareFlops[13] = 37800;
m_nonZeroFlops[14] = 15498;
m_matrixKernels[14] = sgemm_m20_n27_k35_ldA84_ldB36_ldC20_beta0_pfsigonly;
m_hardwareFlops[14] = 37800;
m_nonZeroFlops[15] = 960;
m_matrixKernels[15] = sgemm_m20_n27_k27_ldA20_ldB27_ldC20_beta1_pfsigonly;
m_hardwareFlops[15] = 29160;
m_nonZeroFlops[16] = 1782;
m_matrixKernels[16] = sgemm_m12_n27_k20_ldA84_ldB20_ldC12_beta0_pfsigonly;
m_hardwareFlops[16] = 12960;
m_nonZeroFlops[17] = 4158;
m_matrixKernels[17] = sgemm_m12_n27_k20_ldA84_ldB20_ldC12_beta0_pfsigonly;
m_hardwareFlops[17] = 12960;
m_nonZeroFlops[18] = 4968;
m_matrixKernels[18] = sgemm_m12_n27_k20_ldA84_ldB20_ldC12_beta0_pfsigonly;
m_hardwareFlops[18] = 12960;
m_nonZeroFlops[19] = 480;
m_matrixKernels[19] = sgemm_m12_n27_k27_ldA12_ldB27_ldC12_beta1_pfsigonly;
m_hardwareFlops[19] = 17496;
m_nonZeroFlops[20] = 378;
m_matrixKernels[20] = sgemm_m4_n27_k10_ldA84_ldB12_ldC4_beta0_pfsigonly;
m_hardwareFlops[20] = 2160;
m_nonZeroFlops[21] = 918;
m_matrixKernels[21] = sgemm_m4_n27_k10_ldA84_ldB12_ldC4_beta0_pfsigonly;
m_hardwareFlops[21] = 2160;
m_nonZeroFlops[22] = 1188;
m_matrixKernels[22] = sgemm_m4_n27_k10_ldA84_ldB12_ldC4_beta0_pfsigonly;
m_hardwareFlops[22] = 2160;
m_nonZeroFlops[23] = 192;
m_matrixKernels[23] = sgemm_m4_n27_k27_ldA4_ldB27_ldC4_beta1_pfsigonly;
m_hardwareFlops[23] = 5832;
m_nonZeroFlops[24] = 54;
m_matrixKernels[24] = sgemm_m4_n27_k4_ldA84_ldB4_ldC4_beta0_pfsigonly;
m_hardwareFlops[24] = 864;
m_nonZeroFlops[25] = 108;
m_matrixKernels[25] = sgemm_m4_n27_k4_ldA84_ldB4_ldC4_beta0_pfsigonly;
m_hardwareFlops[25] = 864;
m_nonZeroFlops[26] = 162;
m_matrixKernels[26] = sgemm_m4_n27_k4_ldA84_ldB4_ldC4_beta0_pfsigonly;
m_hardwareFlops[26] = 864;
m_nonZeroFlops[27] = 48;
m_matrixKernels[27] = sgemm_m4_n27_k27_ldA4_ldB27_ldC4_beta1_pfsigonly;
m_hardwareFlops[27] = 5832;
#endif


#ifdef VOLUME_KERNEL
m_nonZeroFlops[0] = 78084;
m_matrixKernels[0] = sgemm_m120_n27_k84_ldA120_ldB120_ldC120_beta0_pfsigonly;
m_hardwareFlops[0] = 544320;
m_nonZeroFlops[1] = 174312;
m_matrixKernels[1] = sgemm_m120_n27_k84_ldA120_ldB120_ldC120_beta0_pfsigonly;
m_hardwareFlops[1] = 544320;
m_nonZeroFlops[2] = 185652;
m_matrixKernels[2] = sgemm_m120_n27_k84_ldA120_ldB120_ldC120_beta0_pfsigonly;
m_hardwareFlops[2] = 544320;
m_nonZeroFlops[3] = 5760;
m_matrixKernels[3] = sgemm_m120_n27_k27_ldA120_ldB27_ldC120_beta1_pfsigonly;
m_hardwareFlops[3] = 174960;
#endif


#ifdef BOUNDARY_KERNEL
m_nonZeroFlops[0] = 29160;
m_hardwareFlops[0] = 777600;
m_matrixKernels[0] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
m_nonZeroFlops[1] = 106272;
m_hardwareFlops[1] = 777600;
m_matrixKernels[1] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
m_nonZeroFlops[2] = 377784;
m_hardwareFlops[2] = 777600;
m_matrixKernels[2] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
m_nonZeroFlops[3] = 377784;
m_hardwareFlops[3] = 777600;
m_matrixKernels[3] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
m_nonZeroFlops[4] = 104544;
m_hardwareFlops[4] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[4] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[4] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[5] = 104544;
m_hardwareFlops[5] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[5] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[5] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[6] = 29160;
m_hardwareFlops[6] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[6] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[6] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[7] = 131220;
m_hardwareFlops[7] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[7] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[7] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[8] = 470448;
m_hardwareFlops[8] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[8] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[8] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[9] = 470448;
m_hardwareFlops[9] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[9] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[9] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[10] = 621486;
m_hardwareFlops[10] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[10] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[10] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[11] = 545940;
m_hardwareFlops[11] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[11] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[11] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[12] = 544860;
m_hardwareFlops[12] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[12] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[12] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[13] = 544860;
m_hardwareFlops[13] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[13] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[13] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[14] = 545940;
m_hardwareFlops[14] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[14] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[14] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[15] = 621486;
m_hardwareFlops[15] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[15] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[15] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[16] = 131220;
m_hardwareFlops[16] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[16] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[16] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[17] = 470448;
m_hardwareFlops[17] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[17] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[17] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[18] = 470448;
m_hardwareFlops[18] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[18] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[18] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[19] = 674028;
m_hardwareFlops[19] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[19] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[19] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[20] = 106272;
m_hardwareFlops[20] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[20] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[20] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[21] = 674028;
m_hardwareFlops[21] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[21] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[21] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[22] = 422280;
m_hardwareFlops[22] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[22] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[22] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[23] = 673164;
m_hardwareFlops[23] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[23] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[23] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[24] = 621432;
m_hardwareFlops[24] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[24] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[24] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[25] = 673164;
m_hardwareFlops[25] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[25] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[25] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[26] = 422280;
m_hardwareFlops[26] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[26] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[26] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[27] = 621432;
m_hardwareFlops[27] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[27] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[27] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[28] = 621486;
m_hardwareFlops[28] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[28] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[28] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[29] = 545940;
m_hardwareFlops[29] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[29] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[29] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[30] = 544860;
m_hardwareFlops[30] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[30] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[30] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[31] = 422280;
m_hardwareFlops[31] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[31] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[31] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[32] = 673164;
m_hardwareFlops[32] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[32] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[32] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[33] = 621432;
m_hardwareFlops[33] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[33] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[33] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[34] = 476064;
m_hardwareFlops[34] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[34] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[34] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[35] = 775332;
m_hardwareFlops[35] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[35] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[35] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[36] = 525852;
m_hardwareFlops[36] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[36] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[36] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[37] = 377784;
m_hardwareFlops[37] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[37] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[37] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[38] = 724734;
m_hardwareFlops[38] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[38] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[38] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[39] = 724734;
m_hardwareFlops[39] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[39] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[39] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[40] = 544860;
m_hardwareFlops[40] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[40] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[40] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[41] = 545940;
m_hardwareFlops[41] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[41] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[41] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[42] = 621486;
m_hardwareFlops[42] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[42] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[42] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[43] = 673164;
m_hardwareFlops[43] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[43] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[43] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[44] = 422280;
m_hardwareFlops[44] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[44] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[44] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[45] = 621432;
m_hardwareFlops[45] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[45] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[45] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[46] = 377784;
m_hardwareFlops[46] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[46] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[46] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[47] = 724734;
m_hardwareFlops[47] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[47] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[47] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[48] = 724734;
m_hardwareFlops[48] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[48] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[48] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[49] = 476064;
m_hardwareFlops[49] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[49] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[49] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[50] = 525852;
m_hardwareFlops[50] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[50] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[50] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[51] = 775332;
m_hardwareFlops[51] = 777600;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[51] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_BL2viaC;
#else
m_matrixKernels[51] = sgemm_m120_n27_k120_ldA120_ldB120_ldC120_beta0_pfsigonly;
#endif
m_nonZeroFlops[52] = 19440;
m_hardwareFlops[52] = 174960;
m_matrixKernels[52] = sgemm_m120_n27_k27_ldA120_ldB27_ldC120_beta1_pfsigonly;
m_nonZeroFlops[53] = 19440;
m_hardwareFlops[53] = 174960;
#ifdef ENABLE_MATRIX_PREFETCH
m_matrixKernels[53] = sgemm_m120_n27_k27_ldA120_ldB27_ldC120_beta1_pfsigonly;
#else
m_matrixKernels[53] = sgemm_m120_n27_k27_ldA120_ldB27_ldC120_beta1_pfsigonly;
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

