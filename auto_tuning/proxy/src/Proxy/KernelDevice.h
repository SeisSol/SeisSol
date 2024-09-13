// SPDX-FileCopyrightInfo: 2020-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "KernelHost.h"

namespace seissol::proxy {

class ProxyKernelDeviceAder : public ProxyKernelHostAder {
  public:
  void run(ProxyData& data, seissol::parallel::runtime::StreamRuntime& runtime) const override;
};

class ProxyKernelDeviceLocalWOAder : public ProxyKernelHostLocalWOAder {
  public:
  void run(ProxyData& data, seissol::parallel::runtime::StreamRuntime& runtime) const override;
};

class ProxyKernelDeviceLocal : public ProxyKernelHostLocal {
  public:
  void run(ProxyData& data, seissol::parallel::runtime::StreamRuntime& runtime) const override;
};

class ProxyKernelDeviceNeighbor : public ProxyKernelHostNeighbor {
  public:
  void run(ProxyData& data, seissol::parallel::runtime::StreamRuntime& runtime) const override;
};

class ProxyKernelDeviceNeighborDR : public ProxyKernelDeviceNeighbor {
  public:
  auto needsDR() const -> bool override;
};

class ProxyKernelDeviceGodunovDR : public ProxyKernelHostGodunovDR {
  public:
  void run(ProxyData& data, seissol::parallel::runtime::StreamRuntime& runtime) const override;
};

using ProxyKernelDeviceAll = CompoundKernel<ProxyKernelDeviceLocal, ProxyKernelDeviceNeighbor>;
using ProxyKernelDeviceAllDR =
    CompoundKernel<ProxyKernelDeviceLocal, ProxyKernelDeviceGodunovDR, ProxyKernelDeviceNeighborDR>;

std::shared_ptr<ProxyKernel> getProxyKernelDevice(Kernel kernel);

} // namespace seissol::proxy
