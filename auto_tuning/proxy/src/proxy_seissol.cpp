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
 
extern long long libxsmm_num_total_flops;
extern long long pspamm_num_total_flops;

#include <sys/time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef USE_MEMKIND
#include <hbwmalloc.h>
#endif

#include <utils/args.h>

#ifdef __MIC__
#define __USE_RDTSC
#endif

#include <Kernels/TimeCommon.h>
#include <Kernels/Time.h>
#include <Kernels/Local.h>
#include <Kernels/Neighbor.h>
#include <Kernels/DynamicRupture.h>
#include "utils/logger.h"

// seissol_kernel includes
#include "proxy_seissol_tools.hpp"
#include "proxy_seissol_allocator.hpp"
#include "proxy_seissol_flops.hpp"
#include "proxy_seissol_bytes.hpp"
#include "proxy_seissol_integrators.hpp"
#ifdef ACL_DEVICE
#include "proxy_seissol_device_integrators.hpp"
#endif

#ifdef ACL_DEVICE
using namespace proxy::device;
#else
using namespace proxy::cpu;
#endif

enum Kernel { all = 0, local, neigh, ader, localwoader, neigh_dr, godunov_dr };
char const* Kernels[] = {"all", "local", "neigh", "ader", "localwoader", "neigh_dr", "godunov_dr"};

void testKernel(unsigned kernel, unsigned timesteps) {
  unsigned t = 0;
  switch (kernel) {
    case all:
      for (; t < timesteps; ++t) {
        computeLocalIntegration();
        computeNeighboringIntegration();
      }
      break;
    case local:
      for (; t < timesteps; ++t) {
        computeLocalIntegration();
      }
      break;
    case neigh:
    case neigh_dr:
      for (; t < timesteps; ++t) {
        computeNeighboringIntegration();
      }
      break;
    case ader:
      for (; t < timesteps; ++t) {
        computeAderIntegration();
      }
      break;
    case localwoader:
      for (; t < timesteps; ++t) {
        computeLocalWithoutAderIntegration();
      }
      break;    
    case godunov_dr:
      for (; t < timesteps; ++t) {
        computeDynRupGodunovState();
      }
      break;
    default:
      break;
  }
}

int main(int argc, char* argv[]) {
  std::stringstream kernelHelp;
  kernelHelp << "Kernel: " << Kernels[0];
  for (size_t k = 1; k < sizeof(Kernels)/sizeof(char*); ++k) {
    kernelHelp << ", " << Kernels[k];
  }

  utils::Args args;
  args.addAdditionalOption("cells", "Number of cells");
  args.addAdditionalOption("timesteps", "Number of timesteps");
  args.addAdditionalOption("kernel", kernelHelp.str());
  
  if (args.parse(argc, argv) != utils::Args::Success) {
    return -1;
  }
  
  unsigned cells = args.getAdditionalArgument<unsigned>("cells");
  unsigned timesteps = args.getAdditionalArgument<unsigned>("timesteps");
  std::string kernelStr = args.getAdditionalArgument<std::string>("kernel");
  unsigned kernel = 0;
  for (; kernel < sizeof(Kernels)/sizeof(char*); ++kernel) {
    if (kernelStr.compare(Kernels[kernel]) == 0) {
      break;
    }
  }
  if (kernel >= sizeof(Kernels)/sizeof(char*)) {
    std::cerr << "Unknown kernel " << kernelStr << std::endl;
    return -1;
  }
  
  bool enableDynamicRupture = false;
  if (kernel == neigh_dr || kernel == godunov_dr) {
    enableDynamicRupture = true;
  }

#ifdef ACL_DEVICE
  deviceT &device = deviceT::getInstance();
  device.api->initialize();
  device.api->setDevice(0);
  device.api->allocateStackMem();
#endif

  m_ltsTree = new seissol::initializers::LTSTree;
  m_dynRupTree = new seissol::initializers::LTSTree;
  m_allocator = new seissol::memory::ManagedAllocator;

  print_hostname();

  printf("Allocating fake data...\n");
  initGlobalData();
  cells = initDataStructures(cells, enableDynamicRupture);
#ifdef ACL_DEVICE
  initDataStructuresOnDevice(enableDynamicRupture);
#endif // ACL_DEVICE
  printf("...done\n\n");

  struct timeval start_time, end_time;
#ifdef __USE_RDTSC
  size_t cycles_start, cycles_end;
#endif
  double total = 0.0;
  double total_cycles = 0.0;

  // init OpenMP and LLC
  testKernel(kernel, 1);
  
  libxsmm_num_total_flops = 0;
  pspamm_num_total_flops = 0;

  gettimeofday(&start_time, NULL);
#ifdef __USE_RDTSC
  cycles_start = __rdtsc();
#endif

  testKernel(kernel, timesteps);

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

  seissol_flops (*flop_fun)(unsigned);
  double (*bytes_fun)(unsigned);
  switch (kernel) {
    case all:
      flop_fun = &flops_all_actual;
      bytes_fun = &bytes_all;
      break;
    case local:
      flop_fun = &flops_local_actual;
      bytes_fun = &bytes_local;
      break;
    case neigh:
    case neigh_dr:
      flop_fun = &flops_neigh_actual;
      bytes_fun = &bytes_neigh;
      break;
    case ader:
      flop_fun = &flops_ader_actual;
      bytes_fun = &noestimate;
      break;
    case localwoader:
      flop_fun = &flops_localWithoutAder_actual;
      bytes_fun = &noestimate;
      break;
    case godunov_dr:
      flop_fun = &flops_drgod_actual;
      bytes_fun = &noestimate;
      break;
  }
  
  seissol_flops actual_flops = (*flop_fun)(timesteps);
  double bytes_estimate = (*bytes_fun)(timesteps);

  printf("=================================================\n");
  printf("===            PERFORMANCE SUMMARY            ===\n");
  printf("=================================================\n");
  printf("seissol proxy mode                  : %s\n", kernelStr.c_str());
  printf("time for seissol proxy              : %f\n", total);
  printf("cycles                              : %f\n\n", total_cycles);
  printf("GFLOP (libxsmm)                     : %f\n", libxsmm_num_total_flops      * 1.e-9);
  printf("GFLOP (pspamm)                      : %f\n", pspamm_num_total_flops           * 1.e-9);
  printf("GFLOP (libxsmm + pspamm)            : %f\n", (libxsmm_num_total_flops + pspamm_num_total_flops) * 1.e-9);
  printf("GFLOP (non-zero) for seissol proxy  : %f\n", actual_flops.d_nonZeroFlops  * 1.e-9);
  printf("GFLOP (hardware) for seissol proxy  : %f\n", actual_flops.d_hardwareFlops * 1.e-9);
  printf("GiB (estimate) for seissol proxy    : %f\n\n", bytes_estimate/(1024.0*1024.0*1024.0));
  printf("FLOPS/cycle (non-zero)              : %f\n", actual_flops.d_nonZeroFlops/total_cycles);
  printf("FLOPS/cycle (hardware)              : %f\n", actual_flops.d_hardwareFlops/total_cycles);
  printf("Bytes/cycle (estimate)              : %f\n\n", bytes_estimate/total_cycles);
  printf("GFLOPS (non-zero) for seissol proxy : %f\n", (actual_flops.d_nonZeroFlops  * 1.e-9)/total);
  printf("GFLOPS (hardware) for seissol proxy : %f\n", (actual_flops.d_hardwareFlops * 1.e-9)/total);
  printf("GiB/s (estimate) for seissol proxy  : %f\n", (bytes_estimate/(1024.0*1024.0*1024.0))/total);
  printf("=================================================\n");
  printf("\n");

  delete m_ltsTree;
  delete m_dynRupTree;
  delete m_allocator;

#ifdef ACL_DEVICE
  device.finalize();
#endif
  return 0;
}

