// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Runner.h"

#include "GeneratedCode/kernel.h"
#include "KernelDevice.h"
#include "KernelHost.h"
#include "Kernels/Common.h"
#include "Parallel/Runtime/Stream.h"
#include "Proxy/Kernel.h"

#include <cstddef>
#include <iostream>
#include <memory>
#include <sys/time.h>
#include <vector>

#ifdef USE_MEMKIND
#include <hbwmalloc.h>
#endif

#include "Common.h"

#ifdef __MIC__
#define __USE_RDTSC
#endif

#include "Monitoring/FlopCounter.h"

#include <cassert>

// seissol_kernel includes
#include "Allocator.h"
#include "Tools.h"

namespace seissol::proxy {

namespace {

void testKernel(std::shared_ptr<ProxyData>& data,
                std::shared_ptr<parallel::runtime::StreamRuntime>& runtime,
                std::shared_ptr<ProxyKernel>& kernel,
                std::uint64_t timesteps) {
  for (std::uint64_t i = 0; i < timesteps; ++i) {
    kernel->run(*data, *runtime);
  }
}

} // namespace

auto runProxy(const ProxyConfig& config) -> ProxyOutput {
  auto kernel = [&]() {
    std::vector<std::shared_ptr<ProxyKernel>> subkernels;
    for (const auto& kernelName : config.kernels) {
      if constexpr (isDeviceOn()) {
        subkernels.emplace_back(getProxyKernelDevice(kernelName));
      } else {
        subkernels.emplace_back(getProxyKernelHost(kernelName));
      }
    }
    return std::dynamic_pointer_cast<ProxyKernel>(std::make_shared<ChainKernel>(subkernels));
  }();

  const bool enableDynamicRupture = kernel->needsDR();

  if (config.verbose) {
    std::cerr << "Allocating fake data... ";
  }

  auto data = std::make_shared<ProxyData>(config.cells, enableDynamicRupture);

  auto runtime = std::make_shared<seissol::parallel::runtime::StreamRuntime>();

  if (config.verbose) {
    std::cerr << "...done" << std::endl;
  }

  struct timeval startTime{};
  struct timeval endTime{};
  double total = 0.0;
  double totalCycles = 0.0;

  // init OpenMP and LLC
  testKernel(data, runtime, kernel, 1);

  runtime->wait();

  gettimeofday(&startTime, nullptr);
  const auto cyclesStart = getCycles();

  testKernel(data, runtime, kernel, config.timesteps);

  runtime->wait();

  const auto cyclesEnd = getCycles();
  gettimeofday(&endTime, nullptr);
  total = sec(startTime, endTime);

  if (cyclesEnd - cyclesStart > 0) {
    totalCycles = static_cast<double>(cyclesEnd - cyclesStart);
  } else {
    totalCycles = deriveCyclesFromTime(total);
  }

  const auto performanceEstimate = kernel->performanceEstimate(*data);

  const double hardwareFlops = config.timesteps * performanceEstimate.hardwareFlop;
  const double nonzeroFlops = config.timesteps * performanceEstimate.nonzeroFlop;
  const double bytesEstimate = config.timesteps * performanceEstimate.bytes;

  ProxyOutput output{};
  output.time = total;
  output.cycles = totalCycles;
  output.libxsmmNumTotalGFlop = static_cast<double>(libxsmm_num_total_flops) * 1.e-9;
  output.pspammNumTotalGFlop = static_cast<double>(pspamm_num_total_flops) * 1.e-9;
  output.libxsmmAndpspammNumTotalGFlop =
      static_cast<double>(libxsmm_num_total_flops + pspamm_num_total_flops) * 1.e-9;
  output.actualNonZeroGFlop = static_cast<double>(nonzeroFlops) * 1.e-9;
  output.actualHardwareGFlop = static_cast<double>(hardwareFlops) * 1.e-9;
  output.gib = bytesEstimate / (1024.0 * 1024.0 * 1024.0);
  output.nonZeroFlopPerCycle = static_cast<double>(nonzeroFlops) / totalCycles;
  output.hardwareFlopPerCycle = static_cast<double>(hardwareFlops) / totalCycles;
  output.bytesPerCycle = bytesEstimate / totalCycles;
  output.nonZeroGFlops = (static_cast<double>(nonzeroFlops) * 1.e-9) / total;
  output.hardwareGFlops = (static_cast<double>(hardwareFlops) * 1.e-9) / total;
  output.gibPerSecond = (bytesEstimate / (1024.0 * 1024.0 * 1024.0)) / total;

  runtime.reset();

  data.reset();

  return output;
}

} // namespace seissol::proxy
