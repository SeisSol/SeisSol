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
  std::size_t cellCount;

  LTS::Storage ltsStorage;
  DynamicRupture::Storage drStorage;

  GlobalData globalDataOnHost;
  GlobalData globalDataOnDevice;

  real* fakeDerivatives = nullptr;
  real* fakeDerivativesHost = nullptr;

  typename kernels::Solver<Config>::TimeBasis<real> timeBasis{Config::ConvergenceOrder};

  kernels::Spacetime<Cfg> spacetimeKernel;
  kernels::Time<Cfg> timeKernel;
  kernels::Local<Cfg> localKernel;
  kernels::Neighbor<Cfg> neighborKernel;
  kernels::DynamicRupture<Cfg> dynRupKernel;

  seissol::memory::ManagedAllocator allocator;

  ProxyData(std::size_t cellCount, bool enableDR);

  initializer::LayerIdentifier layerId{HaloType::Interior, Config(), 0};

  // TODO: check copyability (probably not)

  private:
  void initGlobalData();
  void initDataStructures(bool enableDR);
  void initDataStructuresOnDevice(bool enableDR);
};

} // namespace seissol::proxy

#endif // SEISSOL_SRC_PROXY_ALLOCATOR_H_
