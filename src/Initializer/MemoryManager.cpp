// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#include "MemoryManager.h"

#include "DynamicRupture/Factory.h"
#include "DynamicRupture/Misc.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Kernels/Common.h"
#include "Memory/Descriptor/Boundary.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Descriptor/Surface.h"
#include "Memory/GlobalData.h"
#include "Memory/MemoryAllocator.h"
#include "Memory/Tree/Layer.h"
#include "SeisSol.h"

#include <memory>
#include <utility>
#include <utils/logger.h>
#include <yateto.h>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif // ACL_DEVICE

namespace seissol::initializer {

MemoryManager::MemoryManager(seissol::SeisSol& instance) : seissolInstance_(instance) {}

void MemoryManager::initialize() {
  // initialize global matrices
  GlobalDataInitializerOnHost::init(globalDataOnHost_, memoryAllocator_, memory::Memkind::Standard);
  if constexpr (seissol::isDeviceOn()) {
    GlobalDataInitializerOnDevice::init(
        globalDataOnDevice_, memoryAllocator_, memory::Memkind::DeviceGlobalMemory);
  }
}

void MemoryManager::initializeFrictionLaw() {
  const auto& params = seissolInstance_.parameters().drParameters;
  const auto drParameters = std::make_shared<parameters::DRParameters>(params);
  logInfo() << "Initialize Friction Model";

  logInfo() << "Friction law:" << dr::misc::frictionLawName(drParameters->frictionLawType).c_str()
            << "(" << static_cast<int>(drParameters->frictionLawType) << ")";
  logInfo() << "Thermal pressurization:" << (drParameters->isThermalPressureOn ? "on" : "off");

  const auto factory = seissol::dr::factory::getFactory(drParameters, seissolInstance_);
  auto product = factory->produce();
  dynRup_ = std::move(product.storage);
  drInitializer_ = std::move(product.initializer);
  frictionLaw_ = std::move(product.frictionLaw);
  frictionLawDevice_ = std::move(product.frictionLawDevice);
  faultOutputManager_ = std::move(product.output);
}

void MemoryManager::initFrictionData() {
  const auto& params = seissolInstance_.parameters().drParameters;
  if (params.isDynamicRuptureEnabled) {

    drInitializer_->initializeFault(drStorage_);
  }
}

void MemoryManager::synchronizeTo(AllocationPlace place) {
  if (place == AllocationPlace::Device) {
    logInfo() << "Synchronizing data... (host->device)";
  } else {
    logInfo() << "Synchronizing data... (device->host)";
  }

  void* defaultStream = nullptr;

#ifdef ACL_DEVICE
  defaultStream = device::DeviceInstance::getInstance().api->getDefaultStream();
#endif
  ltsStorage_.synchronizeTo(place, defaultStream);
  drStorage_.synchronizeTo(place, defaultStream);
  boundaryStorage_.synchronizeTo(place, defaultStream);
  surfaceStorage_.synchronizeTo(place, defaultStream);

#ifdef ACL_DEVICE
  device::DeviceInstance::getInstance().api->syncDefaultStreamWithHost();
#endif
}

} // namespace seissol::initializer
