// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Runner.h"

#include "KernelDevice.h"
#include "KernelHost.h"
#include "Kernels/Common.h"
#include "Monitoring/Energy/Discovery.h"
#include "Monitoring/Energy/Sensor.h"
#include "Monitoring/Metric.h"
#include "Parallel/Runtime/Stream.h"
#include "Proxy/Cycles.h"
#include "Proxy/Kernel.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <sys/time.h>
#include <utils/env.h>
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

//! @brief The sensors of the node, and what they read at the start of the run.
struct EnergyMeasurement {
  std::vector<std::unique_ptr<monitoring::energy::Sensor>> sensors{};
  std::vector<std::optional<monitoring::energy::EnergyReading>> start{};
};

EnergyMeasurement startEnergy() {
  EnergyMeasurement measurement;
  if (utils::Env("SEISSOL_").get<bool>("ENERGY", true)) {
    measurement.sensors = monitoring::energy::nodeSensors("/");
  }
  for (auto& sensor : measurement.sensors) {
    measurement.start.push_back(sensor->read());
  }
  return measurement;
}

/**
 * @brief Puts the energy since startEnergy() into @p output , over @p time seconds of the run.
 *
 * A part whose sensors share a clock that did not advance during the run is not counted: the run
 * was too short for them to update.
 */
void endEnergy(EnergyMeasurement& measurement,
               double time,
               double nonzeroFlops,
               double hardwareFlops,
               ProxyOutput& output) {
  using monitoring::energy::Part;
  constexpr std::array<Part, 5> Parts{
      Part::Node, Part::Cpu, Part::Memory, Part::Accelerators, Part::Devices};
  // two spans of the clocks of sensors closer than this, in seconds, are the same
  constexpr double SameTime = 1e-6;
  std::array<std::optional<double>, Parts.size()> joules{};
  // the span of the clock of the sensors of a part, if they all have the same, as the files of the
  // HPE Cray counters have
  std::array<std::optional<double>, Parts.size()> seconds{};
  std::array<bool, Parts.size()> shared{};
  shared.fill(true);
  for (std::size_t i = 0; i < measurement.sensors.size(); ++i) {
    const auto end = measurement.sensors[i]->read();
    const auto& start = measurement.start[i];
    if (!start.has_value() || !end.has_value()) {
      continue;
    }
    for (std::size_t part = 0; part < Parts.size(); ++part) {
      if (measurement.sensors[i]->part() != Parts[part]) {
        continue;
      }
      joules[part] = joules[part].value_or(0.0) + end.value().joules - start.value().joules;
      if (!start.value().seconds.has_value() || !end.value().seconds.has_value()) {
        shared[part] = false;
        continue;
      }
      const double span = end.value().seconds.value() - start.value().seconds.value();
      if (seconds[part].has_value() && std::abs(seconds[part].value() - span) > SameTime) {
        shared[part] = false;
      }
      seconds[part] = span;
    }
  }
  // a part whose sensors share a clock takes its power over it, and is not counted if they did
  // not update during the run, which was too short then
  std::array<bool, Parts.size()> counted{};
  for (std::size_t part = 0; part < Parts.size(); ++part) {
    const bool clocked = shared[part] && seconds[part].has_value();
    counted[part] = joules[part].has_value() && !(clocked && seconds[part].value() <= 0.0);
    if (counted[part]) {
      output.partJoules.emplace_back(monitoring::energy::partName(Parts[part]),
                                     joules[part].value());
    }
  }
  // the basis: the whole node, else its GPUs, else its CPUs, as in the run summary
  for (const auto basis : {Part::Node, Part::Devices, Part::Cpu}) {
    for (std::size_t part = 0; part < Parts.size(); ++part) {
      if (Parts[part] != basis || !counted[part] || joules[part].value() <= 0.0) {
        continue;
      }
      const double value = joules[part].value();
      const bool clocked = shared[part] && seconds[part].has_value();
      output.energyBasis = monitoring::energy::partName(basis);
      output.joules = value;
      const double span = clocked ? seconds[part].value() : time;
      output.watts = span > 0.0 ? value / span : 0.0;
      output.nonZeroGFlopPerJoule = nonzeroFlops * 1.e-9 / value;
      output.hardwareGFlopPerJoule = hardwareFlops * 1.e-9 / value;
      return;
    }
  }
}

void testKernel(const std::shared_ptr<ProxyData>& data,
                const std::shared_ptr<parallel::runtime::StreamRuntime>& runtime,
                const std::shared_ptr<ProxyKernel>& kernel,
                std::size_t timesteps) {
  for (std::size_t i = 0; i < timesteps; ++i) {
    kernel->run(*data, *runtime);
  }
}

auto makeKernel(const ProxyConfig& config) -> std::shared_ptr<ProxyKernel> {
  std::vector<std::shared_ptr<ProxyKernel>> subkernels;
  for (const auto& kernelName : config.kernels) {
    if constexpr (isDeviceOn()) {
      subkernels.emplace_back(getProxyKernelDevice(kernelName));
    } else {
      subkernels.emplace_back(getProxyKernelHost(kernelName));
    }
  }
  return std::dynamic_pointer_cast<ProxyKernel>(std::make_shared<ChainKernel>(subkernels));
}

} // namespace

auto estimateProxy(const ProxyConfig& config) -> PerformanceEstimate {
  const auto kernel = makeKernel(config);
  const auto data = makeProxyData(config.config, config.cells, kernel->needsDR());
  return kernel->performanceEstimate(*data);
}

auto runProxy(const ProxyConfig& config) -> ProxyOutput {
  const auto kernel = makeKernel(config);

  const bool enableDynamicRupture = kernel->needsDR();

  if (config.verbose) {
    std::cerr << "Allocating fake data... ";
  }

  auto data = makeProxyData(config.config, config.cells, enableDynamicRupture);

  auto runtime = std::make_shared<seissol::parallel::runtime::StreamRuntime>();

  if (config.verbose) {
    std::cerr << "...done" << '\n';
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

  auto energy = startEnergy();
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

  endEnergy(energy, total, nonzeroFlops, hardwareFlops, output);

  runtime.reset();

  data.reset();

  return output;
}

} // namespace seissol::proxy
