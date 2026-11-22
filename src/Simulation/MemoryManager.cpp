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

#include "Core/Common/ConfigDispatch.h"
#include "Core/Common/ConfigRegistry.h"
#include "Core/Kernels/Common.h"
#include "Core/Memory/Descriptor/Boundary.h"
#include "Core/Memory/Descriptor/DynamicRupture.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/Descriptor/Surface.h"
#include "Core/Memory/GlobalData.h"
#include "Core/Memory/HostData.h"
#include "Core/Memory/MemoryAllocator.h"
#include "Core/Memory/Tree/Layer.h"
#include "Features/DynamicRupture/DRParameters.h"
#include "Features/DynamicRupture/Misc.h"
#include "Features/DynamicRupture/Output/FaultOutputParameters.h"
#include "Features/DynamicRupture/Output/OutputManager.h"
#include "Features/DynamicRupture/Setup/Factory.h"
#include "Features/DynamicRupture/Setup/Initializer/BaseDRInitializer.h"
#include "Simulation/Parameters/SeisSolParameters.h"
#include "Simulation/Reader/ReaderBuilder.h"
#include "Simulation/SeisSol.h"

#include <memory>
#include <string>
#include <utility>
#include <utils/logger.h>
#include <variant>
#include <vector>
#include <yateto.h>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif // ACL_DEVICE

namespace seissol::initializer {

MemoryManager::MemoryManager(seissol::SeisSol& instance)
    : seissolInstance_(instance), iniConds_(builtConfigCount()) {}

void MemoryManager::initialize() {
  // initialize global matrices
  forEachConfig([&](auto cfg) {
    using Cfg = decltype(cfg);
    GlobalDataInitializerOnHost::init<Cfg>(
        globalData_.onHost<Cfg>(), memoryAllocator_, memory::Memkind::Standard);
    if constexpr (seissol::isDeviceOn()) {
      GlobalDataInitializerOnDevice::init<Cfg>(
          globalData_.onDevice<Cfg>(), memoryAllocator_, memory::Memkind::DeviceGlobalMemory);
    }
  });
}

void MemoryManager::initializeFrictionLaw() {
  const auto& params = seissolInstance_.parameters().drParameters;
  const auto drParameters = std::make_shared<parameters::DRParameters>(params);
  logInfo() << "Initialize Friction Model";

  logInfo() << "Friction law:" << dr::misc::frictionLawName(drParameters->frictionLawType) << "("
            << static_cast<int>(drParameters->frictionLawType) << ")";
  logInfo() << "Thermal pressurization:" << (drParameters->isThermalPressureOn ? "on" : "off");

  const auto factory = seissol::dr::factory::getFactory(
      drParameters,
      dr::initializer::FaultMesh{.mesh = seissolInstance_.meshReader(), .backmap = drBackmap_},
      seissol::data::modelLoader());
  auto product = factory->produce();
  dynRup_ = std::move(product.storage);
  drInitializer_ = std::move(product.initializer);
  frictionLaw_ = std::move(product.frictionLaw);
  frictionLawDevice_ = std::move(product.frictionLawDevice);
  faultOutputManager_ =
      std::make_unique<dr::output::OutputManager>(std::move(product.output), seissolInstance_);
  // the outputs of the fault: its whole surface, which a run has at most once, and its receivers
  const auto& output = seissolInstance_.parameters().output;
  std::vector<std::pair<std::string, parameters::PickpointParameters>> pickpoints;
  for (const auto& entry : output.outputs) {
    if (const auto* settings = std::get_if<parameters::PickpointParameters>(&entry.settings)) {
      pickpoints.emplace_back(entry.name, *settings);
    }
  }
  faultOutputManager_->setOutputs(
      output.singleOutputOf<parameters::ElementwiseFaultParameters>("fault"), pickpoints);
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
  defaultStream = device::DeviceInstance::instance().api().getDefaultStream();
#endif
  ltsStorage_.synchronizeTo(place, defaultStream);
  drStorage_.synchronizeTo(place, defaultStream);
  boundaryStorage_.synchronizeTo(place, defaultStream);
  surfaceStorage_.synchronizeTo(place, defaultStream);

#ifdef ACL_DEVICE
  device::DeviceInstance::instance().api().syncDefaultStreamWithHost();
#endif
}

void MemoryManager::synchronizeToHost(const HostData& data) {
  if (data.isEverything()) {
    synchronizeTo(AllocationPlace::Host);
    return;
  }
  if (data.empty()) {
    return;
  }
  logInfo() << "Synchronizing data... (device->host;" << data.size() << "variables)";
#ifdef ACL_DEVICE
  auto& api = device::DeviceInstance::instance().api();
  data.copyToHost(api.getDefaultStream());
  api.syncDefaultStreamWithHost();
#endif
}

} // namespace seissol::initializer
