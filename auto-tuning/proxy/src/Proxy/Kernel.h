// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_AUTO_TUNING_PROXY_SRC_PROXY_KERNEL_H_
#define SEISSOL_AUTO_TUNING_PROXY_SRC_PROXY_KERNEL_H_

#include "Allocator.h"
#include <Parallel/Runtime/Stream.h>
#include <type_traits>
namespace seissol::proxy {

struct PerformanceEstimate {
  std::size_t hardwareFlop{0};
  std::size_t nonzeroFlop{0};
  std::size_t bytes{0};

  auto operator+(const PerformanceEstimate& other) const -> PerformanceEstimate;
};

class ProxyKernel {
  public:
  virtual void run(ProxyData& data, seissol::parallel::runtime::StreamRuntime& runtime) const = 0;
  virtual auto performanceEstimate(ProxyData& data) const -> PerformanceEstimate = 0;
  [[nodiscard]] virtual auto needsDR() const -> bool = 0;
  virtual ~ProxyKernel() = default;
};

template <typename... Kernels>
class CompoundKernel : public ProxyKernel {
  public:
  static_assert((std::is_base_of_v<ProxyKernel, Kernels> && ...),
                "All classes in a CompoundKernel need to inherit from ProxyKernel");

  void run(ProxyData& data, seissol::parallel::runtime::StreamRuntime& runtime) const override {
    (Kernels().run(data, runtime), ...);
  }

  auto performanceEstimate(ProxyData& data) const -> PerformanceEstimate override {
    return (Kernels().performanceEstimate(data) + ...);
  }

  [[nodiscard]] auto needsDR() const -> bool override { return (Kernels().needsDR() || ...); }
};

class ChainKernel : public ProxyKernel {
  public:
  ChainKernel(const std::vector<std::shared_ptr<ProxyKernel>>& kernels);

  void run(ProxyData& data, seissol::parallel::runtime::StreamRuntime& runtime) const override;

  auto performanceEstimate(ProxyData& data) const -> PerformanceEstimate override;

  [[nodiscard]] auto needsDR() const -> bool override;

  private:
  std::vector<std::shared_ptr<ProxyKernel>> kernels;
};

} // namespace seissol::proxy

#endif // SEISSOL_AUTO_TUNING_PROXY_SRC_PROXY_KERNEL_H_
