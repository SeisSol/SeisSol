// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PROXY_ALLOCATOR_H_
#define SEISSOL_SRC_PROXY_ALLOCATOR_H_

#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/GlobalData.h"
#include "Memory/Tree/LTSTree.h"
#include <Config.h>
#include <Kernels/DynamicRupture.h>
#include <Kernels/Local.h>
#include <Kernels/Neighbor.h>
#include <Kernels/Solver.h>
#include <Memory/Tree/Layer.h>
#include <Parallel/Runtime/Stream.h>
#include <unordered_set>
#include <yateto.h>

#ifdef ACL_DEVICE
#include "Initializer/BatchRecorders/Recorders.h"
#include <device.h>
#endif

namespace seissol::proxy {

struct ProxyData {
  static std::shared_ptr<ProxyData>
      get(ConfigVariant variant, std::size_t cellCount, bool enableDR);
};

template <typename Cfg>
struct ProxyDataImpl : public ProxyData {
  std::size_t cellCount;

  LTS::Storage ltsStorage;
  DynamicRupture::Storage drStorage;

  GlobalData globalData;

  Real<Cfg>* fakeDerivatives = nullptr;
  Real<Cfg>* fakeDerivativesHost = nullptr;

  typename kernels::Solver<Cfg>::template TimeBasis<Real<Cfg>> timeBasis{Cfg::ConvergenceOrder};

  kernels::Spacetime<Cfg> spacetimeKernel;
  kernels::Time<Cfg> timeKernel;
  kernels::Local<Cfg> localKernel;
  kernels::Neighbor<Cfg> neighborKernel;
  kernels::DynamicRupture<Cfg> dynRupKernel;

  seissol::memory::ManagedAllocator allocator;

  ProxyDataImpl(std::size_t cellCount, bool enableDR);

  initializer::LayerIdentifier layerId{HaloType::Interior, Cfg(), 0};

  // TODO: check copyability (probably not)

  private:
  void initGlobalData();
  void initDataStructures(bool enableDR);
  void initDataStructuresOnDevice(bool enableDR);
};

} // namespace seissol::proxy

#endif // SEISSOL_SRC_PROXY_ALLOCATOR_H_
