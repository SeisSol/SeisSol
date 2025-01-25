// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_AUTO_TUNING_PROXY_SRC_PROXY_KERNELHOST_H_
#define SEISSOL_AUTO_TUNING_PROXY_SRC_PROXY_KERNELHOST_H_

#include "Common.h"
#include "Kernel.h"
namespace seissol::proxy {

class ProxyKernelHostAder : public ProxyKernel {
  public:
  void run(ProxyData& data, seissol::parallel::runtime::StreamRuntime& runtime) const override;
  auto performanceEstimate(ProxyData& data) const -> PerformanceEstimate override;
  [[nodiscard]] auto needsDR() const -> bool override;
};

class ProxyKernelHostLocalWOAder : public ProxyKernel {
  public:
  void run(ProxyData& data, seissol::parallel::runtime::StreamRuntime& runtime) const override;
  auto performanceEstimate(ProxyData& data) const -> PerformanceEstimate override;
  [[nodiscard]] auto needsDR() const -> bool override;
};

class ProxyKernelHostLocal
    : public CompoundKernel<ProxyKernelHostAder, ProxyKernelHostLocalWOAder> {
  public:
  void run(ProxyData& data, seissol::parallel::runtime::StreamRuntime& runtime) const override;
};

class ProxyKernelHostNeighbor : public ProxyKernel {
  public:
  void run(ProxyData& data, seissol::parallel::runtime::StreamRuntime& runtime) const override;
  auto performanceEstimate(ProxyData& data) const -> PerformanceEstimate override;
  [[nodiscard]] auto needsDR() const -> bool override;
};

class ProxyKernelHostNeighborDR : public ProxyKernelHostNeighbor {
  public:
  [[nodiscard]] auto needsDR() const -> bool override;
};

class ProxyKernelHostGodunovDR : public ProxyKernel {
  public:
  void run(ProxyData& data, seissol::parallel::runtime::StreamRuntime& runtime) const override;
  auto performanceEstimate(ProxyData& data) const -> PerformanceEstimate override;
  [[nodiscard]] auto needsDR() const -> bool override;
};

using ProxyKernelHostAll = CompoundKernel<ProxyKernelHostLocal, ProxyKernelHostNeighbor>;
using ProxyKernelHostAllDR =
    CompoundKernel<ProxyKernelHostLocal, ProxyKernelHostGodunovDR, ProxyKernelHostNeighborDR>;

std::shared_ptr<ProxyKernel> getProxyKernelHost(Kernel kernel);

} // namespace seissol::proxy

#endif // SEISSOL_AUTO_TUNING_PROXY_SRC_PROXY_KERNELHOST_H_
