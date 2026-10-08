// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PROXY_ALLOCATOR_H_
#define SEISSOL_SRC_PROXY_ALLOCATOR_H_

#include "Common/ConfigDispatch.h"
#include "Common/ConfigRegistry.h"
#include "Common/Real.h"
#include "Kernels/DynamicRupture.h"
#include "Kernels/Local.h"
#include "Kernels/Neighbor.h"
#include "Kernels/Solver.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/GlobalData.h"
#include "Memory/Tree/LTSTree.h"
#include "Memory/Tree/Layer.h"
#include "Parallel/Runtime/Stream.h"

#include <cstddef>
#include <functional>
#include <memory>
#include <unordered_set>
#include <yateto.h>

#ifdef ACL_DEVICE
#include "Initializer/BatchRecorders/Recorders.h"

#include <Device/device.h>
#endif

namespace seissol::proxy {

/// The data the kernels of the proxy run on: one layer of cells and, if needed, of fault faces, in
/// the configuration `config`.
struct ProxyData {
  ProxyData(std::size_t cellCount, ConfigId config);
  virtual ~ProxyData() = default;

  ProxyData(const ProxyData&) = delete;
  ProxyData(ProxyData&&) = delete;
  auto operator=(const ProxyData&) = delete;
  auto operator=(ProxyData&&) = delete;

  std::size_t cellCount;
  ConfigId config;

  LTS::Storage ltsStorage;
  DynamicRupture::Storage drStorage;

  seissol::memory::ManagedAllocator allocator;

  initializer::LayerIdentifier layerId;
};

/// The data of the proxy in the configuration `Cfg`, together with the kernels and the global data
/// of that configuration.
template <typename Cfg>
struct ProxyDataImpl : public ProxyData {
  using real = Real<Cfg>; // NOLINT(readability-identifier-naming)

  ProxyDataImpl(std::size_t cellCount, bool enableDR);

  GlobalData<Cfg> globalDataOnHost;
  GlobalData<Cfg> globalDataOnDevice;

  real* fakeDerivatives = nullptr;
  real* fakeDerivativesHost = nullptr;

  kernels::TimeBasis<Cfg> timeBasis{Cfg::ConvergenceOrder};

  kernels::Spacetime<Cfg> spacetimeKernel;
  kernels::Time<Cfg> timeKernel;
  kernels::Local<Cfg> localKernel;
  kernels::Neighbor<Cfg> neighborKernel;
  kernels::DynamicRupture<Cfg> dynRupKernel;

  private:
  void initGlobalData();
  void initDataStructures(bool enableDR);
  void initDataStructuresOnDevice(bool enableDR);
};

/// Allocates the data of the proxy in the configuration `config`.
std::shared_ptr<ProxyData> makeProxyData(ConfigId config, std::size_t cellCount, bool enableDR);

/// Calls `function` with `data` as the data of its configuration, i.e. as a `ProxyDataImpl<Cfg>`.
template <typename F>
decltype(auto) dispatchProxyData(ProxyData& data, F&& function) {
  return dispatchConfig(data.config, [&](auto cfg) -> decltype(auto) {
    using Cfg = decltype(cfg);
    return std::invoke(std::forward<F>(function), static_cast<ProxyDataImpl<Cfg>&>(data));
  });
}

} // namespace seissol::proxy

#endif // SEISSOL_SRC_PROXY_ALLOCATOR_H_
