// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: David Schneller

#include "Estimator.h"
#include <Common/Executor.h>
#include <Kernels/Common.h>
#include <Numerical/Statistics.h>
#include <Proxy/Common.h>
#include <Proxy/Runner.h>

#include <utils/logger.h>

namespace seissol::solver {

auto miniSeisSol() -> double {
  const auto config = proxy::ProxyConfig{50000,
                                         10,
                                         {seissol::proxy::Kernel::Local},
                                         false,
                                         isDeviceOn() ? Executor::Device : Executor::Host};

  logInfo() << "Running MiniSeisSol with" << config.cells << "cells and" << config.timesteps
            << "repetitions.";
  const auto proxyResult = seissol::proxy::runProxy(config);

  const auto summary = statistics::parallelSummary(proxyResult.time);
  logInfo() << "Runtime results:" << "min:" << summary.min << "max" << summary.max
            << "mean:" << summary.mean << "median:" << summary.median << "stddev:" << summary.std;

  return proxyResult.time;
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
