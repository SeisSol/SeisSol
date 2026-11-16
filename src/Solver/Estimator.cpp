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
#include <mpi.h>
#include <utility>
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

LayerCost::LayerCost(std::vector<std::pair<std::size_t, double>> points)
    : points_(std::move(points)) {
  // noise must not make a larger layer cheaper
  for (std::size_t i = 1; i < points_.size(); ++i) {
    points_[i].second = std::max(points_[i].second, points_[i - 1].second);
  }
}

double LayerCost::operator()(std::size_t cells) const {
  if (cells == 0 || points_.empty()) {
    return 0;
  }
  if (points_.size() == 1) {
    return points_.front().second;
  }
  // the stretch the size falls into; below the smallest size, the first one, and beyond the largest
  // size, the last one
  std::size_t upper = 1;
  while (upper + 1 < points_.size() && cells > points_[upper].first) {
    ++upper;
  }
  const auto [cells0, time0] = points_[upper - 1];
  const auto [cells1, time1] = points_[upper];
  const auto slope = (time1 - time0) / static_cast<double>(cells1 - cells0);
  return std::max(0.0, time0 + slope * (static_cast<double>(cells) - static_cast<double>(cells0)));
}

double LayerCost::combinedSaving(std::size_t first, std::size_t second) const {
  if (first == 0 || second == 0) {
    return 0;
  }
  return (*this)(first) + (*this)(second) - (*this)(first + second);
}

auto measureLayerCost(ConfigId config) -> LayerCost {
  // a device that has been idle raises its clocks only under load (see miniSeisSol)
  constexpr unsigned SmallestCells = 8;
  constexpr unsigned LargestCells = 16384;
  auto warmup = costProxyConfig(config, LargestCells);
  warmup.timesteps = 1;
  warmup.warmup = 0.3;
  proxy::runProxy(warmup);

  // from where the fixed part of a layer dominates to where its cells do: up to where a step takes
  // longer than RunTime on all ranks, beyond which the fixed part weighs nothing and the time grows
  // along a line. Each size runs for about RunTime (by the time of the size before: twice the cells
  // take at most twice as long), with a few steps at least, the best of two runs against noise on
  // the machine.
  constexpr double RunTime = 0.02;
  constexpr unsigned ProbeSteps = 10;
  constexpr double MinSteps = 3;
  constexpr double MaxSteps = 10000;
  constexpr int Repetitions = 2;
  std::vector<std::size_t> sizes;
  std::vector<double> times;
  auto expected = layerStepTime(config, SmallestCells, ProbeSteps, 1);
  for (unsigned cells = SmallestCells; cells <= LargestCells; cells *= 2) {
    const auto steps = std::clamp(RunTime / expected, MinSteps, MaxSteps);
    sizes.push_back(cells);
    times.push_back(layerStepTime(config, cells, static_cast<unsigned>(steps), Repetitions));
    if (Mpi::mpi.allreduce(times.back(), MPI_MIN) > RunTime) {
      break;
    }
    expected = 2 * times.back();
  }

  // the same on all ranks: the median over them, size by size
  const auto ofRanks = Mpi::mpi.collectContainer(times);
  if (Mpi::mpi.rank() == 0) {
    for (std::size_t i = 0; i < times.size(); ++i) {
      std::vector<double> ofSize;
      ofSize.reserve(ofRanks.size());
      for (const auto& rankTimes : ofRanks) {
        ofSize.push_back(rankTimes[i]);
      }
      std::ranges::sort(ofSize);
      const auto middle = ofSize.size() / 2;
      times[i] =
          ofSize.size() % 2 == 1 ? ofSize[middle] : (ofSize[middle - 1] + ofSize[middle]) / 2;
    }
  }
  Mpi::mpi.broadcastContainer(times, 0);

  std::vector<std::pair<std::size_t, double>> points;
  points.reserve(sizes.size());
  for (std::size_t i = 0; i < sizes.size(); ++i) {
    points.emplace_back(sizes[i], times[i]);
  }
  return LayerCost(std::move(points));
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
