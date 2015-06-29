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
// @date 2015-05-09 22:17:48.595721
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
#ifndef SGEMMNOARCHCPP
#define SGEMMNOARCHCPP

#if defined( __SSE3__) || defined(__MIC__)
#include <immintrin.h>
#endif

#include <cstddef>
#ifndef NDEBUG
extern long long libxsmm_num_total_flops;
#endif

void sgemm_m20_n9_k10_ldA20_ldB20_ldC20_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 20; m++) { C[(n*20)+m] = 0.0; }

  for (unsigned int k = 0; k < 10; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 20; m++) {
      C[(n*20)+m] += A[(k*20)+m] * B[(n*20)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 3600;
#endif

}

void sgemm_m4_n9_k1_ldA4_ldB4_ldC4_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 4; m++) { C[(n*4)+m] = 0.0; }

  for (unsigned int k = 0; k < 1; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 4; m++) {
      C[(n*4)+m] += A[(k*4)+m] * B[(n*4)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 72;
#endif

}

void sgemm_m20_n9_k35_ldA20_ldB36_ldC20_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 20; m++) { C[(n*20)+m] = 0.0; }

  for (unsigned int k = 0; k < 35; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 20; m++) {
      C[(n*20)+m] += A[(k*20)+m] * B[(n*36)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 12600;
#endif

}

void sgemm_m4_n9_k10_ldA84_ldB12_ldC4_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 4; m++) { C[(n*4)+m] = 0.0; }

  for (unsigned int k = 0; k < 10; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 4; m++) {
      C[(n*4)+m] += A[(k*84)+m] * B[(n*12)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 720;
#endif

}

void sgemm_m4_n9_k10_ldA56_ldB12_ldC4_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 4; m++) { C[(n*4)+m] = 0.0; }

  for (unsigned int k = 0; k < 10; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 4; m++) {
      C[(n*4)+m] += A[(k*56)+m] * B[(n*12)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 720;
#endif

}

void sgemm_m4_n9_k4_ldA84_ldB4_ldC4_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 4; m++) { C[(n*4)+m] = 0.0; }

  for (unsigned int k = 0; k < 4; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 4; m++) {
      C[(n*4)+m] += A[(k*84)+m] * B[(n*4)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 288;
#endif

}

void sgemm_m56_n9_k35_ldA56_ldB56_ldC56_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 56; m++) { C[(n*56)+m] = 0.0; }

  for (unsigned int k = 0; k < 35; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 56; m++) {
      C[(n*56)+m] += A[(k*56)+m] * B[(n*56)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 35280;
#endif

}

void sgemm_m12_n9_k20_ldA20_ldB20_ldC12_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 12; m++) { C[(n*12)+m] = 0.0; }

  for (unsigned int k = 0; k < 20; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 12; m++) {
      C[(n*12)+m] += A[(k*20)+m] * B[(n*20)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 4320;
#endif

}

void sgemm_m4_n9_k10_ldA12_ldB12_ldC4_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 4; m++) { C[(n*4)+m] = 0.0; }

  for (unsigned int k = 0; k < 10; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 4; m++) {
      C[(n*4)+m] += A[(k*12)+m] * B[(n*12)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 720;
#endif

}

void sgemm_m4_n9_k4_ldA20_ldB4_ldC4_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 4; m++) { C[(n*4)+m] = 0.0; }

  for (unsigned int k = 0; k < 4; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 4; m++) {
      C[(n*4)+m] += A[(k*20)+m] * B[(n*4)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 288;
#endif

}

void sgemm_m12_n9_k10_ldA12_ldB12_ldC12_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 12; m++) { C[(n*12)+m] = 0.0; }

  for (unsigned int k = 0; k < 10; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 12; m++) {
      C[(n*12)+m] += A[(k*12)+m] * B[(n*12)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 2160;
#endif

}

void sgemm_m56_n9_k84_ldA56_ldB84_ldC56_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 56; m++) { C[(n*56)+m] = 0.0; }

  for (unsigned int k = 0; k < 84; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 56; m++) {
      C[(n*56)+m] += A[(k*56)+m] * B[(n*84)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 84672;
#endif

}

void sgemm_m36_n9_k56_ldA84_ldB56_ldC36_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 36; m++) { C[(n*36)+m] = 0.0; }

  for (unsigned int k = 0; k < 56; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 36; m++) {
      C[(n*36)+m] += A[(k*84)+m] * B[(n*56)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 36288;
#endif

}

void sgemm_m4_n9_k4_ldA4_ldB4_ldC4_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 4; m++) { C[(n*4)+m] = 0.0; }

  for (unsigned int k = 0; k < 4; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 4; m++) {
      C[(n*4)+m] += A[(k*4)+m] * B[(n*4)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 288;
#endif

}

void sgemm_m84_n9_k9_ldA84_ldB9_ldC84_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for (unsigned int k = 0; k < 9; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 84; m++) {
      C[(n*84)+m] += A[(k*84)+m] * B[(n*9)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 13608;
#endif

}

void sgemm_m12_n9_k20_ldA12_ldB20_ldC12_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 12; m++) { C[(n*12)+m] = 0.0; }

  for (unsigned int k = 0; k < 20; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 12; m++) {
      C[(n*12)+m] += A[(k*12)+m] * B[(n*20)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 4320;
#endif

}

void sgemm_m4_n9_k10_ldA4_ldB12_ldC4_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 4; m++) { C[(n*4)+m] = 0.0; }

  for (unsigned int k = 0; k < 10; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 4; m++) {
      C[(n*4)+m] += A[(k*4)+m] * B[(n*12)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 720;
#endif

}

void sgemm_m120_n9_k9_ldA120_ldB9_ldC120_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for (unsigned int k = 0; k < 9; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 120; m++) {
      C[(n*120)+m] += A[(k*120)+m] * B[(n*9)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 19440;
#endif

}

void sgemm_m84_n9_k56_ldA84_ldB84_ldC84_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 84; m++) { C[(n*84)+m] = 0.0; }

  for (unsigned int k = 0; k < 56; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 84; m++) {
      C[(n*84)+m] += A[(k*84)+m] * B[(n*84)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 84672;
#endif

}

void sgemm_m20_n9_k9_ldA20_ldB9_ldC20_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for (unsigned int k = 0; k < 9; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 20; m++) {
      C[(n*20)+m] += A[(k*20)+m] * B[(n*9)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 3240;
#endif

}

void sgemm_m36_n9_k20_ldA36_ldB36_ldC36_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 36; m++) { C[(n*36)+m] = 0.0; }

  for (unsigned int k = 0; k < 20; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 36; m++) {
      C[(n*36)+m] += A[(k*36)+m] * B[(n*36)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 12960;
#endif

}

void sgemm_m12_n9_k20_ldA36_ldB20_ldC12_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 12; m++) { C[(n*12)+m] = 0.0; }

  for (unsigned int k = 0; k < 20; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 12; m++) {
      C[(n*12)+m] += A[(k*36)+m] * B[(n*20)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 4320;
#endif

}

void sgemm_m120_n9_k84_ldA120_ldB120_ldC120_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 120; m++) { C[(n*120)+m] = 0.0; }

  for (unsigned int k = 0; k < 84; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 120; m++) {
      C[(n*120)+m] += A[(k*120)+m] * B[(n*120)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 181440;
#endif

}

void sgemm_m20_n9_k35_ldA56_ldB36_ldC20_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 20; m++) { C[(n*20)+m] = 0.0; }

  for (unsigned int k = 0; k < 35; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 20; m++) {
      C[(n*20)+m] += A[(k*56)+m] * B[(n*36)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 12600;
#endif

}

void sgemm_m4_n9_k10_ldA20_ldB12_ldC4_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 4; m++) { C[(n*4)+m] = 0.0; }

  for (unsigned int k = 0; k < 10; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 4; m++) {
      C[(n*4)+m] += A[(k*20)+m] * B[(n*12)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 720;
#endif

}

void sgemm_m4_n9_k4_ldA56_ldB4_ldC4_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 4; m++) { C[(n*4)+m] = 0.0; }

  for (unsigned int k = 0; k < 4; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 4; m++) {
      C[(n*4)+m] += A[(k*56)+m] * B[(n*4)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 288;
#endif

}

void sgemm_m20_n9_k35_ldA84_ldB36_ldC20_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 20; m++) { C[(n*20)+m] = 0.0; }

  for (unsigned int k = 0; k < 35; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 20; m++) {
      C[(n*20)+m] += A[(k*84)+m] * B[(n*36)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 12600;
#endif

}

void sgemm_m12_n9_k4_ldA12_ldB12_ldC12_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 12; m++) { C[(n*12)+m] = 0.0; }

  for (unsigned int k = 0; k < 4; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 12; m++) {
      C[(n*12)+m] += A[(k*12)+m] * B[(n*12)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 864;
#endif

}

void sgemm_m84_n9_k84_ldA84_ldB84_ldC84_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 84; m++) { C[(n*84)+m] = 0.0; }

  for (unsigned int k = 0; k < 84; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 84; m++) {
      C[(n*84)+m] += A[(k*84)+m] * B[(n*84)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 127008;
#endif

}

void sgemm_m56_n9_k9_ldA56_ldB9_ldC56_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for (unsigned int k = 0; k < 9; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 56; m++) {
      C[(n*56)+m] += A[(k*56)+m] * B[(n*9)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 9072;
#endif

}

void sgemm_m36_n9_k56_ldA56_ldB56_ldC36_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 36; m++) { C[(n*36)+m] = 0.0; }

  for (unsigned int k = 0; k < 56; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 36; m++) {
      C[(n*36)+m] += A[(k*56)+m] * B[(n*56)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 36288;
#endif

}

void sgemm_m4_n9_k4_ldA36_ldB4_ldC4_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 4; m++) { C[(n*4)+m] = 0.0; }

  for (unsigned int k = 0; k < 4; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 4; m++) {
      C[(n*4)+m] += A[(k*36)+m] * B[(n*4)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 288;
#endif

}

void sgemm_m20_n9_k35_ldA36_ldB36_ldC20_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 20; m++) { C[(n*20)+m] = 0.0; }

  for (unsigned int k = 0; k < 35; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 20; m++) {
      C[(n*20)+m] += A[(k*36)+m] * B[(n*36)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 12600;
#endif

}

void sgemm_m20_n9_k20_ldA20_ldB20_ldC20_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 20; m++) { C[(n*20)+m] = 0.0; }

  for (unsigned int k = 0; k < 20; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 20; m++) {
      C[(n*20)+m] += A[(k*20)+m] * B[(n*20)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 7200;
#endif

}

void sgemm_m36_n9_k35_ldA36_ldB36_ldC36_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 36; m++) { C[(n*36)+m] = 0.0; }

  for (unsigned int k = 0; k < 35; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 36; m++) {
      C[(n*36)+m] += A[(k*36)+m] * B[(n*36)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 22680;
#endif

}

void sgemm_m56_n9_k56_ldA56_ldB56_ldC56_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 56; m++) { C[(n*56)+m] = 0.0; }

  for (unsigned int k = 0; k < 56; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 56; m++) {
      C[(n*56)+m] += A[(k*56)+m] * B[(n*56)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 56448;
#endif

}

void sgemm_m120_n9_k120_ldA120_ldB120_ldC120_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 120; m++) { C[(n*120)+m] = 0.0; }

  for (unsigned int k = 0; k < 120; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 120; m++) {
      C[(n*120)+m] += A[(k*120)+m] * B[(n*120)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 259200;
#endif

}

void sgemm_m84_n9_k120_ldA84_ldB120_ldC84_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 84; m++) { C[(n*84)+m] = 0.0; }

  for (unsigned int k = 0; k < 120; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 84; m++) {
      C[(n*84)+m] += A[(k*84)+m] * B[(n*120)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 181440;
#endif

}

void sgemm_m12_n9_k9_ldA12_ldB9_ldC12_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for (unsigned int k = 0; k < 9; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 12; m++) {
      C[(n*12)+m] += A[(k*12)+m] * B[(n*9)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 1944;
#endif

}

void sgemm_m12_n9_k20_ldA56_ldB20_ldC12_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 12; m++) { C[(n*12)+m] = 0.0; }

  for (unsigned int k = 0; k < 20; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 12; m++) {
      C[(n*12)+m] += A[(k*56)+m] * B[(n*20)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 4320;
#endif

}

void sgemm_m56_n9_k84_ldA84_ldB84_ldC56_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 56; m++) { C[(n*56)+m] = 0.0; }

  for (unsigned int k = 0; k < 84; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 56; m++) {
      C[(n*56)+m] += A[(k*84)+m] * B[(n*84)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 84672;
#endif

}

void sgemm_m4_n9_k9_ldA4_ldB9_ldC4_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for (unsigned int k = 0; k < 9; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 4; m++) {
      C[(n*4)+m] += A[(k*4)+m] * B[(n*9)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 648;
#endif

}

void sgemm_m4_n9_k4_ldA12_ldB4_ldC4_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 4; m++) { C[(n*4)+m] = 0.0; }

  for (unsigned int k = 0; k < 4; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 4; m++) {
      C[(n*4)+m] += A[(k*12)+m] * B[(n*4)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 288;
#endif

}

void sgemm_m36_n9_k9_ldA36_ldB9_ldC36_beta1_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for (unsigned int k = 0; k < 9; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 36; m++) {
      C[(n*36)+m] += A[(k*36)+m] * B[(n*9)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 5832;
#endif

}

void sgemm_m4_n9_k0_ldA4_ldB4_ldC4_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 4; m++) { C[(n*4)+m] = 0.0; }

  for (unsigned int k = 0; k < 0; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 4; m++) {
      C[(n*4)+m] += A[(k*4)+m] * B[(n*4)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 0;
#endif

}

void sgemm_m4_n9_k10_ldA36_ldB12_ldC4_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 4; m++) { C[(n*4)+m] = 0.0; }

  for (unsigned int k = 0; k < 10; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 4; m++) {
      C[(n*4)+m] += A[(k*36)+m] * B[(n*12)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 720;
#endif

}

void sgemm_m12_n9_k20_ldA84_ldB20_ldC12_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 12; m++) { C[(n*12)+m] = 0.0; }

  for (unsigned int k = 0; k < 20; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 12; m++) {
      C[(n*12)+m] += A[(k*84)+m] * B[(n*20)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 4320;
#endif

}

void sgemm_m36_n9_k56_ldA36_ldB56_ldC36_beta0_pfsigonly(const float* A, const float* B, float* C, const float* A_prefetch = NULL, const float* B_prefetch = NULL, const float* C_prefetch = NULL)
{
#pragma message ("KERNEL COMPILATION WARNING: compiling arch-independent gemm kernel in: " __FILE__)

for (unsigned int n = 0; n < 9; n++) {
  for(unsigned int m = 0; m < 36; m++) { C[(n*36)+m] = 0.0; }

  for (unsigned int k = 0; k < 56; k++) {
    #pragma simd
    for(unsigned int m = 0; m < 36; m++) {
      C[(n*36)+m] += A[(k*36)+m] * B[(n*56)+k];
    }
  }
}
#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 36288;
#endif

}

#endif
