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
#include "Proxy/Cycles.h"
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
                std::size_t timesteps) {
  for (std::size_t i = 0; i < timesteps; ++i) {
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
  std::uint64_t cyclesStart = 0;
  std::uint64_t cyclesEnd = 0;
  double total = 0.0;
  double totalCycles = 0.0;

  // init OpenMP and LLC
  testKernel(data, runtime, kernel, 1);

  runtime->wait();

  const seissol::monitoring::FlopCounter flopCounter{};

  gettimeofday(&startTime, nullptr);
  cyclesStart = readCycles();

  testKernel(data, runtime, kernel, config.timesteps);

  runtime->wait();

  cyclesEnd = readCycles();
  gettimeofday(&endTime, nullptr);
  total = sec(startTime, endTime);
  totalCycles = static_cast<double>(cyclesEnd - cyclesStart);

  const auto performanceEstimate = kernel->performanceEstimate(*data);

  const double hardwareFlops = config.timesteps * performanceEstimate.hardwareFlop;
  const double nonzeroFlops = config.timesteps * performanceEstimate.nonzeroFlop;
  const double bytesEstimate = config.timesteps * performanceEstimate.bytes;
  const double bytesKernel = config.timesteps * performanceEstimate.kernelBytes;

  ProxyOutput output{};
  output.time = total;
  output.cycles = totalCycles;
  output.cycleSource = cycleSourceName();
  output.libxsmmNumTotalGFlop = static_cast<double>(libxsmm_num_total_flops) * 1.e-9;
  output.pspammNumTotalGFlop = static_cast<double>(pspamm_num_total_flops) * 1.e-9;
  output.libxsmmAndpspammNumTotalGFlop =
      static_cast<double>(libxsmm_num_total_flops + pspamm_num_total_flops) * 1.e-9;
  output.actualNonZeroGFlop = static_cast<double>(nonzeroFlops) * 1.e-9;
  output.actualHardwareGFlop = static_cast<double>(hardwareFlops) * 1.e-9;
  output.gib = bytesEstimate / (1024.0 * 1024.0 * 1024.0);
  output.kernelGib = bytesKernel / (1024.0 * 1024.0 * 1024.0);

  // Without a tick counter there is nothing to divide by, and on a run short
  // enough to fit between two ticks the elapsed count is legitimately zero.
  // Reporting zero keeps every field finite, which is what the JSON output
  // needs to stay parseable; cycleSource tells a reader which case this is.
  const auto perCycle = [totalCycles](double value) {
    return totalCycles > 0.0 ? value / totalCycles : 0.0;
  };
  output.nonZeroFlopPerCycle = perCycle(static_cast<double>(nonzeroFlops));
  output.hardwareFlopPerCycle = perCycle(static_cast<double>(hardwareFlops));
  output.bytesPerCycle = perCycle(bytesEstimate);
  output.kernelBytesPerCycle = perCycle(bytesKernel);

  const auto perSecond = [total](double value) { return total > 0.0 ? value / total : 0.0; };
  output.nonZeroGFlops = perSecond(static_cast<double>(nonzeroFlops) * 1.e-9);
  output.hardwareGFlops = perSecond(static_cast<double>(hardwareFlops) * 1.e-9);
  output.gibPerSecond = perSecond(bytesEstimate / (1024.0 * 1024.0 * 1024.0));
  output.kernelGibPerSecond = perSecond(bytesKernel / (1024.0 * 1024.0 * 1024.0));

  runtime.reset();

  data.reset();

  return output;
}

} // namespace seissol::proxy
