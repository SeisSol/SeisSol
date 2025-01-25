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
#include <Proxy/Common.h>
#include <Proxy/Runner.h>

namespace seissol::solver {

auto miniSeisSol() -> double {
  const auto config = proxy::ProxyConfig{50000,
                                         10,
                                         {seissol::proxy::Kernel::Local},
                                         false,
                                         isDeviceOn() ? Executor::Device : Executor::Host};
  const auto proxyResult = seissol::proxy::runProxy(config);

  return proxyResult.time;
}

auto hostDeviceSwitch() -> int {
  if constexpr (!isDeviceOn()) {
    return 0;
  }

  unsigned clusterSize = 1;
  for (int i = 0; i < 20; ++i) {
    auto config = proxy::ProxyConfig{clusterSize,
                                     clusterSize < 100 ? 100 : 10,
                                     {seissol::proxy::Kernel::Local},
                                     false,
                                     Executor::Host};
    const auto resultHost = proxy::runProxy(config);
    config.executor = Executor::Device;
    const auto resultDevice = proxy::runProxy(config);

    if (resultHost.time > resultDevice.time) {
      return clusterSize;
    }

    clusterSize *= 2;
  }

  return 0;
}

} // namespace seissol::solver
