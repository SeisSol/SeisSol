// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PROXY_KERNELDEVICE_H_
#define SEISSOL_SRC_PROXY_KERNELDEVICE_H_

#include "KernelHost.h"

namespace seissol::proxy {

template <typename Cfg>
class ProxyKernelDeviceAder : public ProxyKernelHostAder<Cfg> {
  public:
  void run(ProxyData& predata, seissol::parallel::runtime::StreamRuntime& runtime) const override;
};

template <typename Cfg>
class ProxyKernelDeviceLocalWOAder : public ProxyKernelHostLocalWOAder<Cfg> {
  public:
  void run(ProxyData& predata, seissol::parallel::runtime::StreamRuntime& runtime) const override;
};

template <typename Cfg>
class ProxyKernelDeviceLocal : public ProxyKernelHostLocal<Cfg> {
  public:
  void run(ProxyData& predata, seissol::parallel::runtime::StreamRuntime& runtime) const override;
};

template <typename Cfg>
class ProxyKernelDeviceNeighbor : public ProxyKernelHostNeighbor<Cfg> {
  public:
  void run(ProxyData& predata, seissol::parallel::runtime::StreamRuntime& runtime) const override;
};

template <typename Cfg>
class ProxyKernelDeviceNeighborDR : public ProxyKernelDeviceNeighbor<Cfg> {
  public:
  [[nodiscard]] auto needsDR() const -> bool override;
};

template <typename Cfg>
class ProxyKernelDeviceGodunovDR : public ProxyKernelHostGodunovDR<Cfg> {
  public:
  void run(ProxyData& predata, seissol::parallel::runtime::StreamRuntime& runtime) const override;
};

template <typename Cfg>
using ProxyKernelDeviceAll =
    CompoundKernel<ProxyKernelDeviceLocal<Cfg>, ProxyKernelDeviceNeighbor<Cfg>>;

template <typename Cfg>
using ProxyKernelDeviceAllDR = CompoundKernel<ProxyKernelDeviceLocal<Cfg>,
                                              ProxyKernelDeviceGodunovDR<Cfg>,
                                              ProxyKernelDeviceNeighborDR<Cfg>>;

std::shared_ptr<ProxyKernel> getProxyKernelDevice(Kernel kernel, ConfigVariant variant);

} // namespace seissol::proxy

#endif // SEISSOL_SRC_PROXY_KERNELDEVICE_H_
