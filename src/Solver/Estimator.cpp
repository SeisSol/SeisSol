// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: David Schneller

#include "Estimator.h"

#include "Common/ConfigRegistry.h"
#include "Common/ConfigValue.h"
#include "Common/Executor.h"
#include "Common/Real.h"
#include "Kernels/Common.h"
#include "Numerical/Statistics.h"
#include "Parallel/MPI.h"
#include "Proxy/Common.h"
#include "Proxy/Runner.h"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <utils/logger.h>
#include <vector>

namespace seissol::solver {

namespace {

/// The proxy run that compares the configurations of a run: the time, local and neighbor kernels
/// of `cells` cells over a few steps.
auto costProxyConfig(ConfigId config, unsigned cells) -> proxy::ProxyConfig {
  constexpr unsigned Timesteps = 5;
  return proxy::ProxyConfig{.cells = cells,
                            .timesteps = Timesteps,
                            .kernels = {proxy::Kernel::All},
                            .verbose = false,
                            .executor = isDeviceOn() ? Executor::Device : Executor::Host,
                            .config = config};
}

/// The time the cost run of `config` takes, the best of a few runs against noise on the machine.
auto measuredCost(ConfigId config) -> double {
  // enough cells to leave the caches of a host, and to keep a device busy
  constexpr unsigned Cells = isDeviceOn() ? 50000 : 10000;
  constexpr int Repetitions = 3;
  auto best = std::numeric_limits<double>::infinity();
  for (int i = 0; i < Repetitions; ++i) {
    best = std::min(best, proxy::runProxy(costProxyConfig(config, Cells)).time);
  }
  return best;
}

/// The hardware FLOPs of the cost run of `config`, weighted with the size of its reals: a vector
/// register holds half as many reals in double precision as in single precision.
auto estimatedCost(ConfigId config) -> double {
  // the count per cell barely depends on the number of cells
  constexpr unsigned Cells = 1000;
  const auto flop = proxy::estimateProxy(costProxyConfig(config, Cells)).hardwareFlop;
  const auto realSize = sizeOfRealType(configValue(config).precision);
  return static_cast<double>(flop) * static_cast<double>(realSize);
}

/// One step of the cost run of `config` with `cells` cells: the best of `repetitions` runs of
/// `steps` steps each.
auto layerStepTime(ConfigId config, unsigned cells, unsigned steps, int repetitions) -> double {
  auto proxyConfig = costProxyConfig(config, cells);
  proxyConfig.timesteps = steps;
  auto best = std::numeric_limits<double>::infinity();
  for (int i = 0; i < repetitions; ++i) {
    best = std::min(best, proxy::runProxy(proxyConfig).time / steps);
  }
  return best;
}

} // namespace

auto layerCombineThreshold(ConfigId config) -> std::size_t {
  if constexpr (!isDeviceOn()) {
    return 0;
  }
  // a device that has been idle raises its clocks only under load (see miniSeisSol)
  // the devices measured so far (November 2026) fill up below 2000 cells, so the time per cell
  // comes from two sizes well above that, which stay small next to the memory of a device
  constexpr unsigned SmallCells = 64;
  constexpr unsigned MediumCells = 8192;
  constexpr unsigned LargeCells = 32768;
  auto warmup = costProxyConfig(config, LargeCells);
  warmup.timesteps = 1;
  warmup.warmup = 0.3;
  proxy::runProxy(warmup);

  const auto fixedTime = layerStepTime(config, SmallCells, 50, 2);
  const auto cellTime =
      (layerStepTime(config, LargeCells, 10, 2) - layerStepTime(config, MediumCells, 10, 2)) /
      static_cast<double>(LargeCells - MediumCells);
  double cells = cellTime > 0 ? fixedTime / cellTime : 0.0;

  // the same threshold on all ranks
  const auto summary = statistics::parallelSummary(cells);
  cells = summary.median;
  Mpi::mpi.broadcast(&cells, 0);
  logInfo() << "A step of the cell kernels of" << configName(configValue(config)) << "takes"
            << fixedTime * 1e6 << "us with" << SmallCells << "cells, and" << cellTime * 1e9
            << "ns per further cell: below" << cells << "cells (median over the ranks; min"
            << summary.min << ", max" << summary.max
            << "), its fixed part weighs more than its cells.";
  return static_cast<std::size_t>(cells);
}

auto miniSeisSol(ConfigId config) -> double {
  // A device that has been idle raises its clocks within some hundred milliseconds under load; ten
  // steps (~20 ms on a B200) measured how fast it did that, not how fast it computes: on Horizon
  // (November 2026), some nodes came out 1.5 times as fast as the others, but were 3 % faster in
  // 1000 steps, and their ranks got 1.5 times the cells. So the kernels run for a while first.
  constexpr bool OnDevice = isDeviceOn();
  auto proxyConfig = proxy::ProxyConfig{.cells = 50000,
                                        .timesteps = OnDevice ? 100U : 10U,
                                        .kernels = {seissol::proxy::Kernel::Local},
                                        .verbose = false,
                                        .executor = OnDevice ? Executor::Device : Executor::Host,
                                        .config = config};
  proxyConfig.warmup = OnDevice ? 0.5 : 0.0;

  if (OnDevice) {
    logInfo() << "Running MiniSeisSol with" << proxyConfig.cells << "cells and"
              << proxyConfig.timesteps << "repetitions, after" << proxyConfig.warmup
              << "s of warm-up.";
  } else {
    logInfo() << "Running MiniSeisSol with" << proxyConfig.cells << "cells and"
              << proxyConfig.timesteps << "repetitions.";
  }
  const auto proxyResult = seissol::proxy::runProxy(proxyConfig);

  const auto summary = statistics::parallelSummary(proxyResult.time);
  logInfo() << "Runtime results:" << "min:" << summary.min << "max" << summary.max
            << "mean:" << summary.mean << "median:" << summary.median << "stddev:" << summary.std;

  return proxyResult.time;
}

auto configCostFactors(const std::vector<ConfigId>& configs, ConfigId reference, bool measure)
    -> std::vector<double> {
  const auto cost = [&](ConfigId config) {
    return measure ? measuredCost(config) : estimatedCost(config);
  };

  logInfo() << "Determining the cost of a cell per configuration relative to"
            << configName(configValue(reference))
            << (measure ? "by running MiniSeisSol." : "from the FLOPs of its kernels.");

  const auto referenceCost = cost(reference);
  std::vector<double> factors(builtConfigCount(), 1.0);
  for (const auto config : configs) {
    if (config == reference) {
      continue;
    }
    auto factor = cost(config) / referenceCost;
    if (measure) {
      // the same factors on all ranks, so that they all weight the cells alike
      const auto summary = statistics::parallelSummary(factor);
      factor = summary.median;
      Mpi::mpi.broadcast(&factor, 0);
      logInfo() << "Relative cost of a cell of" << configName(configValue(config))
                << ": median =" << factor << " min =" << summary.min << " max =" << summary.max;
    } else {
      logInfo() << "Relative cost of a cell of" << configName(configValue(config)) << ":" << factor;
    }
    factors[config] = factor;
  }
  return factors;
}

auto hostDeviceSwitch() -> int {
  if constexpr (!isDeviceOn()) {
    return 0;
  }

  unsigned clusterSize = 1;
  bool found = false;
  logInfo() << "Running host-device switchpoint detection test";
  for (int i = 0; i < 20; ++i) {
    auto config =
        proxy::ProxyConfig{.cells = clusterSize,
                           .timesteps = static_cast<unsigned int>(clusterSize < 100 ? 100 : 10),
                           .kernels = {seissol::proxy::Kernel::Local},
                           .verbose = false,
                           .executor = Executor::Host};
    const auto resultHost = proxy::runProxy(config);
    config.executor = Executor::Device;
    const auto resultDevice = proxy::runProxy(config);

    if (resultHost.time > resultDevice.time) {
      clusterSize -= 1;
      found = true;
      break;
    }

    clusterSize *= 2;
  }

  if (!found) {
    clusterSize = 0;
  }

  const auto summary = statistics::parallelSummary(clusterSize);

  logInfo() << "Switchpoint summary:" << "min:" << summary.min << "max" << summary.max
            << "mean:" << summary.mean << "median:" << summary.median << "stddev:" << summary.std;

  return 0;
}

} // namespace seissol::solver
