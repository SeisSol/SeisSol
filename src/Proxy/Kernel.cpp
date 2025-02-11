// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Kernel.h"
#include "Allocator.h"
#include <Parallel/Runtime/Stream.h>
#include <algorithm>
#include <memory>
#include <vector>

namespace seissol::proxy {

auto PerformanceEstimate::operator+(const PerformanceEstimate& other) const -> PerformanceEstimate {
  return PerformanceEstimate{
      hardwareFlop + other.hardwareFlop, nonzeroFlop + other.nonzeroFlop, bytes + other.bytes};
}

ChainKernel::ChainKernel(const std::vector<std::shared_ptr<ProxyKernel>>& kernels)
    : kernels(kernels) {}

void ChainKernel::run(ProxyData& data, seissol::parallel::runtime::StreamRuntime& runtime) const {
  for (const auto& kernel : kernels) {
    kernel->run(data, runtime);
  }
}

auto ChainKernel::performanceEstimate(ProxyData& data) const -> PerformanceEstimate {
  PerformanceEstimate estimate{};

  for (const auto& kernel : kernels) {
    estimate = estimate + kernel->performanceEstimate(data);
  }

  return estimate;
}

auto ChainKernel::needsDR() const -> bool {
  return std::any_of(
      kernels.begin(), kernels.end(), [](const auto& kernel) { return kernel->needsDR(); });
}

} // namespace seissol::proxy
