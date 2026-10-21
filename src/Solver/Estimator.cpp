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
#include <limits>
#include <utils/logger.h>
#include <vector>

namespace seissol::solver {

namespace {

/// The proxy run that compares the configurations of a run: the time, local and neighbor kernels
/// of `cells` cells over a few steps.
auto costProxyConfig(ConfigId config, unsigned cells) -> proxy::ProxyConfig {
  constexpr unsigned Timesteps = 5;
  return proxy::ProxyConfig{cells,
                            Timesteps,
                            {proxy::Kernel::All},
                            false,
                            isDeviceOn() ? Executor::Device : Executor::Host,
                            config};
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

} // namespace

auto miniSeisSol(ConfigId config) -> double {
  const auto proxyConfig = proxy::ProxyConfig{50000,
                                              10,
                                              {seissol::proxy::Kernel::Local},
                                              false,
                                              isDeviceOn() ? Executor::Device : Executor::Host,
                                              config};

  logInfo() << "Running MiniSeisSol with" << proxyConfig.cells << "cells and"
            << proxyConfig.timesteps << "repetitions.";
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
    auto config = proxy::ProxyConfig{clusterSize,
                                     static_cast<unsigned int>(clusterSize < 100 ? 100 : 10),
                                     {seissol::proxy::Kernel::Local},
                                     false,
                                     Executor::Host};
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
