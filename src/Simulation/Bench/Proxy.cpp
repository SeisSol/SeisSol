// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Proxy.h"

#include "Core/Common/Json.h"
#include "Core/Monitoring/FlopCounter.h"
#include "Core/Monitoring/Metric.h"
#include "Simulation/Bench/Cycles.h"
#include "Simulation/Bench/LoopBench.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <utils/stringutils.h>
#include <vector>

namespace seissol::proxy {

namespace {

using KernelName = std::pair<Kernel, std::string_view>;

/// the kernels by their names, in the order the help lists them
constexpr std::array<KernelName, 6> KernelNames{{
    {Kernel::All, "all"},
    {Kernel::AllDR, "all_dr"},
    {Kernel::Local, "local"},
    {Kernel::Neighbor, "neigh"},
    {Kernel::NeighborDR, "neigh_dr"},
    {Kernel::GodunovDR, "godunov_dr"},
}};

/// The entry of KernelNames that `matches`; nullptr if there is none.
template <typename F>
auto findKernel(const F& matches) -> const KernelName* {
  for (const auto& entry : KernelNames) {
    if (matches(entry)) {
      return &entry;
    }
  }
  return nullptr;
}

auto makeBench(const ProxyConfig& config) -> std::unique_ptr<LoopBench> {
  return makeLoopBench(
      config.config, config.cells, needsDynamicRupture(config.kernels), config.executor);
}

/// The work of one step of the kernels of `config`.
auto stepEstimate(const LoopBench& bench, const ProxyConfig& config) -> PerformanceEstimate {
  PerformanceEstimate estimate{};
  for (const auto kernel : config.kernels) {
    estimate += bench.estimate(kernel);
  }
  return estimate;
}

} // namespace

auto Aux::kernel2str(const std::vector<Kernel>& kernels) -> std::string {
  std::string output;
  for (std::size_t i = 0; i < kernels.size(); ++i) {
    if (i > 0) {
      output += ',';
    }
    const auto* entry =
        findKernel([&](const KernelName& candidate) { return candidate.first == kernels[i]; });
    if (entry == nullptr) {
      throw std::runtime_error("Unknown kernel type");
    }
    output += entry->second;
  }
  return output;
}

auto Aux::str2kernel(const std::string& kernelStr) -> std::vector<Kernel> {
  auto parts = utils::StringUtils::split(kernelStr, ',');
  for (auto& part : parts) {
    utils::StringUtils::trim(part);
  }

  std::vector<Kernel> kernels(parts.size());
  for (std::size_t i = 0; i < parts.size(); ++i) {
    const auto* entry =
        findKernel([&](const KernelName& candidate) { return candidate.second == parts[i]; });
    if (entry == nullptr) {
      throw std::runtime_error("Unknown kernel string: " + parts[i]);
    }
    kernels[i] = entry->first;
  }

  return kernels;
}

auto Aux::getAllowedKernels() -> std::vector<std::string> {
  std::vector<std::string> kernels{};
  kernels.reserve(KernelNames.size());
  for (const auto& [kernel, name] : KernelNames) {
    kernels.emplace_back(name);
  }
  return kernels;
}

void Aux::writeOutput(std::ostream& stream,
                      const ProxyOutput& output,
                      const std::string& kernelStr,
                      OutputFormat format) {
  if (format == OutputFormat::Plain) {
    stream << "=================================================\n";
    stream << "===            PERFORMANCE SUMMARY            ===\n";
    stream << "=================================================\n";
    stream << "seissol proxy mode                  : " << kernelStr << '\n';
    stream << "time for seissol proxy              : " << output.time << '\n';
    stream << "cycles                              : " << output.cycles << '\n';
    stream << '\n';
    stream << "GFLOP (libxsmm)                     : " << output.libxsmmNumTotalGFlop << '\n';
    stream << "GFLOP (pspamm)                      : " << output.pspammNumTotalGFlop << '\n';
    stream << "GFLOP (libxsmm + pspamm)            : " << output.libxsmmAndpspammNumTotalGFlop
           << '\n';
    stream << "GFLOP (non-zero) for seissol proxy  : " << output.actualNonZeroGFlop << '\n';
    stream << "GFLOP (hardware) for seissol proxy  : " << output.actualHardwareGFlop << '\n';
    stream << "GiB (estimate) for seissol proxy    : " << output.gib << '\n';
    stream << "GiB (kernels) for seissol proxy     : " << output.kernelGib << '\n';
    stream << '\n';
    stream << "FLOPS/cycle (non-zero)              : " << output.nonZeroFlopPerCycle << '\n';
    stream << "FLOPS/cycle (hardware)              : " << output.hardwareFlopPerCycle << '\n';
    stream << "Bytes/cycle (estimate)              : " << output.bytesPerCycle << '\n';
    stream << "Bytes/cycle (kernel)                : " << output.kernelBytesPerCycle << '\n';
    stream << '\n';
    stream << "GFLOPS (non-zero) for seissol proxy : " << output.nonZeroGFlops << '\n';
    stream << "GFLOPS (hardware) for seissol proxy : " << output.hardwareGFlops << '\n';
    stream << "GiB/s (estimate) for seissol proxy  : " << output.gibPerSecond << '\n';
    stream << "GiB/s (kernel) for seissol proxy    : " << output.kernelGibPerSecond << '\n';
    stream << "=================================================\n";
    stream << '\n';
  } else {
    // Written for machine consumption, on one line. No non-finite value may reach it: JSON has no
    // spelling for those, and the readers expect a number for every field. runProxy keeps every
    // field finite; this is the second line of defence.
    const auto finite = [](double value) { return std::isfinite(value) ? value : 0.0; };
    json::Value document;
    document["name"] = kernelStr;
    document["cycle-source"] = output.cycleSource;
    document["time"] = finite(output.time);
    document["cycles"] = finite(static_cast<double>(output.cycles));
    document["gflop-libxsmm"] = finite(output.libxsmmNumTotalGFlop);
    document["gflop-pspamm"] = finite(output.pspammNumTotalGFlop);
    document["gflop-libxsmm-pspamm"] = finite(output.libxsmmAndpspammNumTotalGFlop);
    document["gflop-nz"] = finite(output.actualNonZeroGFlop);
    document["gflop-hw"] = finite(output.actualHardwareGFlop);
    document["gib"] = finite(output.gib);
    document["gib-kernel"] = finite(output.kernelGib);
    document["gflopcycle-nz"] = finite(output.nonZeroFlopPerCycle);
    document["gflopcycle-hw"] = finite(output.hardwareFlopPerCycle);
    document["gibcycle"] = finite(output.bytesPerCycle);
    document["gibcycle-kernel"] = finite(output.kernelBytesPerCycle);
    document["gflops-nz"] = finite(output.nonZeroGFlops);
    document["gflops-hw"] = finite(output.hardwareGFlops);
    document["gibs"] = finite(output.gibPerSecond);
    document["gibs-kernel"] = finite(output.kernelGibPerSecond);
    stream << document.dump(-1) << '\n';
  }
}

void Aux::displayOutput(const ProxyOutput& output, const std::string& kernelStr) {
  writeOutput(std::cout, output, kernelStr, OutputFormat::Plain);
}

auto needsDynamicRupture(const std::vector<Kernel>& kernels) -> bool {
  return std::ranges::any_of(kernels, [](Kernel kernel) {
    return kernel == Kernel::AllDR || kernel == Kernel::NeighborDR || kernel == Kernel::GodunovDR;
  });
}

auto estimateProxy(const ProxyConfig& config) -> PerformanceEstimate {
  return stepEstimate(*makeBench(config), config);
}

auto runProxy(const ProxyConfig& config) -> ProxyOutput {
  if (config.verbose) {
    std::cerr << "Allocating fake data... ";
  }

  const auto bench = makeBench(config);

  if (config.verbose) {
    std::cerr << "...done" << '\n';
  }

  const auto runSteps = [&](std::size_t steps) {
    for (std::size_t step = 0; step < steps; ++step) {
      for (const auto kernel : config.kernels) {
        bench->run(kernel);
      }
    }
    bench->wait();
  };

  // init OpenMP and LLC
  runSteps(1);

  if (config.warmup > 0) {
    const auto warmupStart = std::chrono::steady_clock::now();
    while (std::chrono::duration<double>(std::chrono::steady_clock::now() - warmupStart).count() <
           config.warmup) {
      runSteps(10);
    }
  }

  const auto startTime = std::chrono::steady_clock::now();
  const auto cyclesStart = readCycles();

  runSteps(config.timesteps);

  const auto cyclesEnd = readCycles();
  const auto endTime = std::chrono::steady_clock::now();
  const double total = std::chrono::duration<double>(endTime - startTime).count();
  const auto totalCycles = static_cast<double>(cyclesEnd - cyclesStart);

  const auto performanceEstimate = stepEstimate(*bench, config);

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

  return output;
}

} // namespace seissol::proxy
