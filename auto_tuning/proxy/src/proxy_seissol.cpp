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
 
#include <Parallel/Runtime/Stream.h>
#include <sys/time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef USE_MEMKIND
#include <hbwmalloc.h>
#endif

#include "likwid_wrapper.h"
#include <utils/args.h>
#include "proxy_common.hpp"

#ifdef __MIC__
#define __USE_RDTSC
#endif

#include "Kernels/TimeCommon.h"
#include "Kernels/Time.h"
#include "Kernels/Local.h"
#include "Kernels/Neighbor.h"
#include "Kernels/DynamicRupture.h"
#include "Monitoring/FlopCounter.h"
#include "utils/logger.h"
#include <cassert>

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

namespace {
static inline void advanceTd(unsigned t, unsigned phase, bool skip) {
#ifdef USE_TARGETDART
  if (!skip) {
    td_advance(1);
    if (t % phase == 0) {
      td_phase_progress(1);
    }
  }
#endif
}
} // namespace

void testKernel(Kernel kernel, unsigned timesteps, unsigned phase, bool skip) {
  unsigned t = 0;
  switch (kernel) {
    case Kernel::All:
      for (; t < timesteps; ++t) {
        computeLocalIntegration();
        computeNeighboringIntegration();
        advanceTd(t, phase, skip);
      }
      break;
    case Kernel::Local:
      for (; t < timesteps; ++t) {
        computeLocalIntegration();
        advanceTd(t, phase, skip);
      }
      break;
    case Kernel::Neighbor:
    case Kernel::NeighborDR:
      for (; t < timesteps; ++t) {
        computeNeighboringIntegration();
        advanceTd(t, phase, skip);
      }
      break;
    case Kernel::Ader:
      for (; t < timesteps; ++t) {
        computeAderIntegration();
        advanceTd(t, phase, skip);
      }
      break;
    case Kernel::LocalWOAder:
      for (; t < timesteps; ++t) {
        computeLocalWithoutAderIntegration();
        advanceTd(t, phase, skip);
      }
      break;    
    case Kernel::GodunovDR:
      for (; t < timesteps; ++t) {
        computeDynRupGodunovState();
        advanceTd(t, phase, skip);
      }
      break;
    case dynrup:
      for (; t < timesteps; ++t) {
        computeDynamicRupture();
        advanceTd(t, phase, skip);
      }
      break;
    default:
      break;
  }
}


ProxyOutput runProxy(ProxyConfig config) {
  LIKWID_MARKER_INIT;

  registerMarkers();

  bool enableDynamicRupture = false;
  if (config.kernel == Kernel::NeighborDR || config.kernel == Kernel::GodunovDR || config.kernel == Kernel::Dynrup) {
    enableDynamicRupture = true;
  }

#ifdef ACL_DEVICE
#ifdef USE_MPI
  seissol::MPI::mpi.bindAcceleratorDevice();
#endif // USE_MPI
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
  device.api->initialize();
  device.api->allocateStackMem();
#endif // ACL_DEVICE

  m_ltsTree = new seissol::initializer::LTSTree;
  m_dynRupTree = new seissol::initializer::LTSTree;
  m_allocator = new seissol::memory::ManagedAllocator;

  print_hostname();

  if (config.verbose)
    printf("Allocating fake data...\n");

  initGlobalData();
  config.cells = initDataStructures(config.cells, enableDynamicRupture, config.fault);
#ifdef ACL_DEVICE
  initDataStructuresOnDevice(enableDynamicRupture);
#endif // ACL_DEVICE

  runtime = new seissol::parallel::runtime::StreamRuntime();

  if (config.verbose) {
    printf("...done\n\n");
  }

  struct timeval start_time, end_time;
#ifdef __USE_RDTSC
  size_t cycles_start, cycles_end;
#endif
  double total = 0.0;
  double total_cycles = 0.0;

  // run Godunov DR kernel first
  if (config.kernel == dynrup) {
    testKernel(godunov_dr, 1, 10, true);
  }

  // init OpenMP and LLC
  testKernel(config.kernel, 1, 10, true);

  runtime->wait();

  seissol::monitoring::FlopCounter flopCounter;

  gettimeofday(&start_time, NULL);
#ifdef __USE_RDTSC
  cycles_start = __rdtsc();
#endif

  testKernel(config.kernel, config.timesteps, config.phase, false);

  runtime->wait();

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

  seissol_flops (*flop_fun)(unsigned) = nullptr;
  double (*bytes_fun)(unsigned) = nullptr;
  switch (config.kernel) {
    case Kernel::All:
      flop_fun = &flops_all_actual;
      bytes_fun = &bytes_all;
      break;
    case Kernel::Local:
      flop_fun = &flops_local_actual;
      bytes_fun = &bytes_local;
      break;
    case Kernel::Neighbor:
    case Kernel::NeighborDR:
      flop_fun = &flops_neigh_actual;
      bytes_fun = &bytes_neigh;
      break;
    case Kernel::Ader:
      flop_fun = &flops_ader_actual;
      bytes_fun = &bytes_noestimate;
      break;
    case Kernel::LocalWOAder:
      flop_fun = &flops_localWithoutAder_actual;
      bytes_fun = &bytes_noestimate;
      break;
    case Kernel::GodunovDR:
      flop_fun = &flops_drgod_actual;
      bytes_fun = &bytes_noestimate;
      break;
    case dynrup:
      flop_fun = &flops_noestimate;
      bytes_fun = &bytes_noestimate;
      break;
  } 

  assert(flop_fun != nullptr);
  assert(bytes_fun != nullptr);

  seissol_flops actual_flops = (*flop_fun)(config.timesteps);
  double bytes_estimate = (*bytes_fun)(config.timesteps);

  ProxyOutput output{};
  output.time = total;
  output.cycles = total_cycles;
  output.libxsmmNumTotalGFlop = static_cast<double>(libxsmm_num_total_flops) * 1.e-9;
  output.pspammNumTotalGFlop = static_cast<double>(pspamm_num_total_flops) * 1.e-9;
  output.libxsmmAndpspammNumTotalGFlop = static_cast<double>(libxsmm_num_total_flops + pspamm_num_total_flops) * 1.e-9;
  output.actualNonZeroGFlop = static_cast<double>(actual_flops.d_nonZeroFlops)  * 1.e-9;
  output.actualHardwareGFlop = static_cast<double>(actual_flops.d_hardwareFlops) * 1.e-9;
  output.gib = bytes_estimate/(1024.0*1024.0*1024.0);
  output.nonZeroFlopPerCycle = static_cast<double>(actual_flops.d_nonZeroFlops)/total_cycles;
  output.hardwareFlopPerCycle = static_cast<double>(actual_flops.d_hardwareFlops)/total_cycles;
  output.bytesPerCycle = bytes_estimate/total_cycles;
  output.nonZeroGFlops = (static_cast<double>(actual_flops.d_nonZeroFlops)  * 1.e-9)/total;
  output.hardwareGFlops = (static_cast<double>(actual_flops.d_hardwareFlops) * 1.e-9)/total;
  output.gibPerSecond = (bytes_estimate/(1024.0*1024.0*1024.0))/total;

  m_frictionSolver.reset();
  delete m_ltsTree;
  delete m_dynRupTree;
  delete m_allocator;

  delete runtime;

#ifdef ACL_DEVICE
  device.finalize();
#endif

  LIKWID_MARKER_CLOSE;
  return output;
}

