/*
Copyright (c) 2015, Intel Corporation

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
 * Copyright (c) 2013-2014, SeisSol Group
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 **/
 
long long libxsmm_num_total_flops = 0;

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <immintrin.h>
#include <sys/time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef USE_MEMKIND
#include <hbwmalloc.h>
//#define USE_HBM_DOFS
//#define USE_HBM_TDOFS
//#define USE_HBM_DERS
//#define USE_HBM_CELLLOCAL_LOCAL
//#define USE_HBM_CELLLOCAL_NEIGH
//#define USE_HBM_GLOBALDATA
#endif

#ifdef __MIC__
#define __USE_RDTSC
#endif

double derive_cycles_from_time(double time) {
  // first try to read proxy env variable with freq
  char* p_freq;
  double d_freq;
  double cycles = 1.0;
  p_freq = getenv ("SEISSOL_PROXY_FREQUENCY");
  if (p_freq !=NULL ) {
    d_freq = atof(p_freq);
    printf("detected frequency (SEISSOL_PROXY_FREQUENCY): %f\n", d_freq);
    cycles = time * d_freq * 1.0e6;
  } else {
    FILE* fp;
    fp = popen("lscpu | grep MHz | awk '{print $3}'", "r");
    if(fp > 0) {
      char tmp_buffer[20];
      fread(tmp_buffer, 20, 1, fp);
      d_freq = atof(tmp_buffer);
      printf("detected frequency (lscpu): %f\n", d_freq);
      cycles = time * d_freq * 1.0e6;
      pclose(fp);
    } else {
      cycles = 1.0;
      printf("detected frequency (lscpu) FAILED!\n");
    }
  }
  return cycles;
}

#include <generated_code/init.h>
#include <generated_code/flops.h>
#include <Initializer/typedefs.hpp>
#include <Initializer/MemoryAllocator.h>

#include <Kernels/TimeCommon.h>
#include <Kernels/Time.h>
#ifdef REQUIRE_SOURCE_MATRIX
#include <Kernels/Local.h>
#include <Kernels/Neighbor.h>
#else
#include <Kernels/Volume.h>
#include <Kernels/Boundary.h>
#endif

#include <omp.h>

// seissol_kernel includes
#include "proxy_seissol_allocator.hpp"
#include "proxy_seissol_flops.hpp"
#include "proxy_seissol_bytes.hpp"
#include "proxy_seissol_integrators.hpp"

inline double sec(struct timeval start, struct timeval end) {
  return ((double)(((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)))) / 1.0e6;
}

int main(int argc, char* argv[]) {
  if (argc != 4) {
    printf("Wrong parameters!\n");
    printf(" #cells #timesteps kernel\n");
    printf("   kernel-values: all, local, neigh, ader, localwoader\n");
    return -1;
  }

  unsigned int i_cells = atoi(argv[1]);
  unsigned int i_timesteps = atoi(argv[2]);
  std::string s_part;
  s_part.assign(argv[3]);

  // double-check if the selected kernel exists
  if ( (s_part.compare("all") != 0) &&
       (s_part.compare("local") != 0) &&
       (s_part.compare("neigh") != 0) &&
       (s_part.compare("ader") != 0) &&
       (s_part.compare("localwoader")) )
  {
    printf("Wrong parameters!\n");
    printf(" #cells #timesteps kernel\n");
    printf("   kernel-values: all, local, neigh, ader, localwoader\n");
    return -1;
  }

  printf("Allocating fake data...\n");
  i_cells = init_data_structures(i_cells);
  printf("...done\n\n");

  struct timeval start_time, end_time;
  size_t cycles_start, cycles_end;
  double total = 0.0;
  double total_cycles = 0.0;

  // init OpenMP and LLC
  if (s_part.compare("all") == 0) {
    computeLocalIntegration();
    computeNeighboringIntegration();
  } else if (s_part.compare("local") == 0) {
    computeLocalIntegration();
  } else if (s_part.compare("neigh") == 0) {
    computeNeighboringIntegration();
  } else if (s_part.compare("ader") == 0) {
    computeAderIntegration();
  } else {
    computeLocalWithoutAderIntegration();
  }

  gettimeofday(&start_time, NULL);
#ifdef __USE_RDTSC
  cycles_start = __rdtsc();
#endif

  if (s_part.compare("all") == 0) {
    for (unsigned int t = 0; t < i_timesteps; t++) {
      computeLocalIntegration();
      computeNeighboringIntegration();
    }
  } else if (s_part.compare("local") == 0) {
    for (unsigned int t = 0; t < i_timesteps; t++) {
      computeLocalIntegration();
    }
  } else if (s_part.compare("neigh") == 0) {
    for (unsigned int t = 0; t < i_timesteps; t++) {
      computeNeighboringIntegration();
    }
  } else if (s_part.compare("ader") == 0) {
    for (unsigned int t = 0; t < i_timesteps; t++) {
      computeAderIntegration();
    }
  } else {
    for (unsigned int t = 0; t < i_timesteps; t++) {
      computeLocalWithoutAderIntegration();
    }
  }
#ifdef __USE_RDTSC  
  cycles_end = __rdtsc();
#endif
  gettimeofday(&end_time, NULL);
  total = sec(start_time, end_time);
#ifdef __USE_RDTSC
  printf("Cycles via __rdtsc()!\n");
  total_cycles = (double)(cycles_end-cycles_start);
#else
  total_cycles = derive_cycles_from_time(total);
#endif

  printf("=================================================\n");
  printf("===            PERFORMANCE SUMMARY            ===\n");
  printf("=================================================\n");
  printf("seissol proxy mode                  : %s\n", s_part.c_str());
  printf("time for seissol proxy              : %f\n", total);
  printf("cycles                              : %f\n\n", total_cycles);

  seissol_flops (*flop_fun)(unsigned);
  double (*bytes_fun)(unsigned);
  if (s_part.compare("all") == 0) {
    flop_fun = &flops_all_actual;
    bytes_fun = &bytes_all;
  } else if (s_part.compare("local") == 0) {
    flop_fun = &flops_local_actual;
    bytes_fun = &bytes_local;
  } else if (s_part.compare("neigh") == 0) {
    flop_fun = &flops_neigh_actual;
    bytes_fun = &bytes_neigh;
  } else if (s_part.compare("ader") == 0) {
    flop_fun = &flops_ader_actual;
    bytes_fun = &noestimate;
  } else {
    flop_fun = &flops_localWithoutAder_actual;
    bytes_fun = &noestimate;
  }
  
  seissol_flops actual_flops = (*flop_fun)(i_timesteps);
  double bytes_estimate = (*bytes_fun)(i_timesteps);
  printf("GFLOP (non-zero) for seissol proxy  : %f\n", actual_flops.d_nonZeroFlops/(1e9));
  printf("GFLOP (hardware) for seissol proxy  : %f\n", actual_flops.d_hardwareFlops/(1e9));
  printf("GiB (estimate) for seissol proxy    : %f\n\n", bytes_estimate/(1024.0*1024.0*1024.0));
  printf("FLOPS/cycle (non-zero)              : %f\n", actual_flops.d_nonZeroFlops/total_cycles);
  printf("FLOPS/cycle (hardware)              : %f\n", actual_flops.d_hardwareFlops/total_cycles);
  printf("Bytes/cycle (estimate)              : %f\n\n", bytes_estimate/total_cycles);
  printf("GFLOPS (non-zero) for seissol proxy : %f\n", (actual_flops.d_nonZeroFlops/(1e9))/total);
  printf("GFLOPS (hardware) for seissol proxy : %f\n", (actual_flops.d_hardwareFlops/(1e9))/total);
  printf("GiB/s (estimate) for seissol proxy  : %f\n", (bytes_estimate/(1024.0*1024.0*1024.0))/total);
  printf("=================================================\n");
  printf("\n");
  
  return 0;
}

