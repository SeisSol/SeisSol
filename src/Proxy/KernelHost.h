// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PROXY_KERNELHOST_H_
#define SEISSOL_SRC_PROXY_KERNELHOST_H_

#include "Common.h"
#include "Kernel.h"
namespace seissol::proxy {

template <typename Cfg>
class ProxyKernelHostAder : public ProxyKernel {
  public:
  void run(ProxyData& predata, seissol::parallel::runtime::StreamRuntime& runtime) const override;
  auto performanceEstimate(ProxyData& predata) const -> PerformanceEstimate override;
  [[nodiscard]] auto needsDR() const -> bool override;
};

template <typename Cfg>
class ProxyKernelHostLocalWOAder : public ProxyKernel {
  public:
  void run(ProxyData& predata, seissol::parallel::runtime::StreamRuntime& runtime) const override;
  auto performanceEstimate(ProxyData& predata) const -> PerformanceEstimate override;
  [[nodiscard]] auto needsDR() const -> bool override;
};

template <typename Cfg>
class ProxyKernelHostLocal
    : public CompoundKernel<ProxyKernelHostAder<Cfg>, ProxyKernelHostLocalWOAder<Cfg>> {
  public:
  void run(ProxyData& predata, seissol::parallel::runtime::StreamRuntime& runtime) const override;
};

template <typename Cfg>
class ProxyKernelHostNeighbor : public ProxyKernel {
  public:
  void run(ProxyData& predata, seissol::parallel::runtime::StreamRuntime& runtime) const override;
  auto performanceEstimate(ProxyData& predata) const -> PerformanceEstimate override;
  [[nodiscard]] auto needsDR() const -> bool override;
};

template <typename Cfg>
class ProxyKernelHostNeighborDR : public ProxyKernelHostNeighbor<Cfg> {
  public:
  [[nodiscard]] auto needsDR() const -> bool override;
};

template <typename Cfg>
class ProxyKernelHostGodunovDR : public ProxyKernel {
  public:
  void run(ProxyData& predata, seissol::parallel::runtime::StreamRuntime& runtime) const override;
  auto performanceEstimate(ProxyData& predata) const -> PerformanceEstimate override;
  [[nodiscard]] auto needsDR() const -> bool override;
};

template <typename Cfg>
using ProxyKernelHostAll = CompoundKernel<ProxyKernelHostLocal<Cfg>, ProxyKernelHostNeighbor<Cfg>>;

template <typename Cfg>
using ProxyKernelHostAllDR = CompoundKernel<ProxyKernelHostLocal<Cfg>,
                                            ProxyKernelHostGodunovDR<Cfg>,
                                            ProxyKernelHostNeighborDR<Cfg>>;

std::shared_ptr<ProxyKernel> getProxyKernelHost(Kernel kernel, ConfigVariant variant);

} // namespace seissol::proxy

#endif // SEISSOL_SRC_PROXY_KERNELHOST_H_
