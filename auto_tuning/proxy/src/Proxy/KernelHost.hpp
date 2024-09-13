#pragma once

#include "Common.hpp"
#include "Kernel.hpp"
#include <cstddef>
namespace seissol::proxy {

class ProxyKernelHostAder : public ProxyKernel {
  public:
  void run(ProxyData& data, seissol::parallel::runtime::StreamRuntime& runtime) const override;
  auto performanceEstimate(ProxyData& data) const -> PerformanceEstimate override;
  auto needsDR() const -> bool override;
};

class ProxyKernelHostLocalWOAder : public ProxyKernel {
  public:
  void run(ProxyData& data, seissol::parallel::runtime::StreamRuntime& runtime) const override;
  auto performanceEstimate(ProxyData& data) const -> PerformanceEstimate override;
  auto needsDR() const -> bool override;
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
  auto needsDR() const -> bool override;
};

class ProxyKernelHostNeighborDR : public ProxyKernelHostNeighbor {
  public:
  auto needsDR() const -> bool override;
};

class ProxyKernelHostGodunovDR : public ProxyKernel {
  public:
  void run(ProxyData& data, seissol::parallel::runtime::StreamRuntime& runtime) const override;
  auto performanceEstimate(ProxyData& data) const -> PerformanceEstimate override;
  auto needsDR() const -> bool override;
};

using ProxyKernelHostAll = CompoundKernel<ProxyKernelHostLocal, ProxyKernelHostNeighbor>;
using ProxyKernelHostAllDR =
    CompoundKernel<ProxyKernelHostLocal, ProxyKernelHostGodunovDR, ProxyKernelHostNeighborDR>;

std::shared_ptr<ProxyKernel> getProxyKernelHost(Kernel kernel);

} // namespace seissol::proxy
