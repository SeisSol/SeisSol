// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SEISSOL_AUTO_TUNING_PROXY_SRC_PROXY_ALLOCATOR_H_
#define SEISSOL_AUTO_TUNING_PROXY_SRC_PROXY_ALLOCATOR_H_

#include "Initializer/DynamicRupture.h"
#include "Initializer/GlobalData.h"
#include "Initializer/LTS.h"
#include "Initializer/Tree/LTSTree.h"
#include <Kernels/DynamicRupture.h>
#include <Kernels/Local.h>
#include <Kernels/Neighbor.h>
#include <Parallel/Runtime/Stream.h>
#include <unordered_set>
#include <yateto.h>

#ifdef ACL_DEVICE
#include "Initializer/BatchRecorders/Recorders.h"
#include <device.h>
#include <unordered_set>
#endif

namespace seissol::proxy {

struct ProxyData {
  std::size_t cellCount;

  seissol::initializer::LTSTree ltsTree;
  seissol::initializer::LTS lts;
  seissol::initializer::LTSTree dynRupTree;
  seissol::initializer::DynamicRupture dynRup;

  GlobalData globalDataOnHost;
  GlobalData globalDataOnDevice;

  real* fakeDerivatives = nullptr;
  real* fakeDerivativesHost = nullptr;

  kernels::Time timeKernel;
  kernels::Local localKernel;
  kernels::Neighbor neighborKernel;
  kernels::DynamicRupture dynRupKernel;

  seissol::memory::ManagedAllocator allocator;

  ProxyData(std::size_t cellCount, bool enableDR);

  // TODO: check copyability (probably not)

  private:
  void initGlobalData();
  void initDataStructures(bool enableDR);
  void initDataStructuresOnDevice(bool enableDR);
};

} // namespace seissol::proxy

#endif // SEISSOL_AUTO_TUNING_PROXY_SRC_PROXY_ALLOCATOR_H_