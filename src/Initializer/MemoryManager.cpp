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

#include "BatchRecorders/Recorders.h"
#include "Common/Constants.h"
#include "Common/Iterator.h"
#include "DynamicRupture/Factory.h"
#include "DynamicRupture/Misc.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BatchRecorders/Recorders.h"
#include "Initializer/BoundaryHelper.h"
#include "Initializer/CellLocalInformation.h"
#include "Initializer/InitProcedure/Internal/Scratchpads.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Memory/Descriptor/Boundary.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Descriptor/Surface.h"
#include "Memory/GlobalData.h"
#include "Memory/MemoryAllocator.h"
#include "Memory/Tree/Layer.h"
#include "SeisSol.h"

#include <array>
#include <cstddef>
#include <memory>
#include <utility>
#include <utils/logger.h>
#include <yateto.h>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif // ACL_DEVICE

namespace seissol::initializer {

void MemoryManager::initialize() {
  // initialize global matrices
  GlobalDataInitializerOnHost::init(globalDataOnHost_, memoryAllocator_, memory::Memkind::Standard);
  if constexpr (seissol::isDeviceOn()) {
    GlobalDataInitializerOnDevice::init(
        globalDataOnDevice_, memoryAllocator_, memory::Memkind::DeviceGlobalMemory);
  }
}

void MemoryManager::fixateLtsStorage() {
  if constexpr (isDeviceOn()) {
    seissol::initializer::internal::deriveRequiredScratchpadMemoryForDr(drStorage_);
    drStorage_.allocateScratchPads();
  }
}

void MemoryManager::fixateBoundaryStorage() {
  const LayerMask ghostMask(Ghost);

  boundaryTree_.setName("boundary");

  // Boundary face storage
  Boundary::addTo(boundaryTree_);
  boundaryTree_.setLayerCount(ltsStorage_.getColorMap());
  boundaryTree_.fixate();

  // Iterate over layers of standard lts storage and face lts storage together.
  for (auto [layer, boundaryLayer] :
       seissol::common::zip(ltsStorage_.leaves(ghostMask), boundaryTree_.leaves(ghostMask))) {
    CellLocalInformation* cellInformation = layer.var<LTS::CellInformation>();

    std::size_t numberOfBoundaryFaces = 0;
    const auto layerSize = layer.size();

#pragma omp parallel for schedule(static) reduction(+ : numberOfBoundaryFaces)
    for (std::size_t cell = 0; cell < layerSize; ++cell) {
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        if (requiresNodalFlux(cellInformation[cell].faceTypes[face])) {
          ++numberOfBoundaryFaces;
        }
      }
    }
    boundaryLayer.setNumberOfCells(numberOfBoundaryFaces);
  }
  boundaryTree_.allocateVariables();
  boundaryTree_.touchVariables();

  // The boundary storage is now allocated, now we only need to map from cell lts
  // to face lts.
  // We do this by, once again, iterating over both storages at the same time.
  for (auto [layer, boundaryLayer] :
       seissol::common::zip(ltsStorage_.leaves(ghostMask), boundaryTree_.leaves(ghostMask))) {
    auto* cellInformation = layer.var<LTS::CellInformation>();
    auto* boundaryMapping = layer.var<LTS::BoundaryMapping>();
    auto* boundaryMappingDevice = layer.var<LTS::BoundaryMappingDevice>();
    auto* faceInformation = boundaryLayer.var<Boundary::FaceInformation>(AllocationPlace::Host);
    auto* faceInformationDevice =
        boundaryLayer.var<Boundary::FaceInformation>(AllocationPlace::Device);

    std::size_t boundaryFace = 0;
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        if (requiresNodalFlux(cellInformation[cell].faceTypes[face])) {
          boundaryMapping[cell][face] = CellBoundaryMapping(faceInformation[boundaryFace]);
          boundaryMappingDevice[cell][face] =
              CellBoundaryMapping(faceInformationDevice[boundaryFace]);
          ++boundaryFace;
        } else {
          boundaryMapping[cell][face] = CellBoundaryMapping();
          boundaryMappingDevice[cell][face] = CellBoundaryMapping();
        }
      }
    }
  }

  surfaceStorage_.setName("surface");
  SurfaceLTS::addTo(surfaceStorage_);

  int refinement = 0;
  const auto& outputParams = seissolInstance_.getSeisSolParameters().output;
  if (outputParams.freeSurfaceParameters.enabled &&
      outputParams.freeSurfaceParameters.vtkorder < 0) {
    refinement = outputParams.freeSurfaceParameters.refinement;
  }
  seissolInstance_.freeSurfaceIntegrator().initialize(refinement, ltsStorage_, surfaceStorage_);
}

void MemoryManager::initializeMemoryLayout() {
  if constexpr (isDeviceOn()) {
    seissol::initializer::internal::deriveRequiredScratchpadMemoryForWp(
        seissolInstance_.getSeisSolParameters().model.plasticity, ltsStorage_);
    ltsStorage_.allocateScratchPads();
  }
}

void MemoryManager::initializeEasiBoundaryReader(const char* fileName) {
  const auto fileNameStr = std::string{fileName};
  if (!fileNameStr.empty()) {
    easiBoundary_ = EasiBoundary(fileNameStr);
  }
}

void MemoryManager::recordExecutionPaths(bool usePlasticity) {
  recording::CompositeRecorder<LTS::LTSVarmap> recorder;
  recorder.addRecorder(new recording::LocalIntegrationRecorder);
  recorder.addRecorder(new recording::NeighIntegrationRecorder);

  if (usePlasticity) {
    recorder.addRecorder(new recording::PlasticityRecorder);
  }

  for (auto& layer : ltsStorage_.leaves(Ghost)) {
    recorder.record(layer);
  }

  recording::CompositeRecorder<DynamicRupture::DynrupVarmap> drRecorder;
  drRecorder.addRecorder(new recording::DynamicRuptureRecorder);
  for (auto& layer : drStorage_.leaves(Ghost)) {
    drRecorder.record(layer);
  }
}

void MemoryManager::initializeFrictionLaw() {
  const auto& params = seissolInstance_.getSeisSolParameters().drParameters;
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

void MemoryManager::initFaultOutputManager(const std::string& backupTimeStamp) {
  const auto& params = seissolInstance_.getSeisSolParameters().drParameters;
  // TODO: switch dynRup_ to shared or weak pointer
  if (params.isDynamicRuptureEnabled) {
    faultOutputManager_->setInputParam(seissolInstance_.meshReader());
    faultOutputManager_->setLtsData(ltsStorage_, backmap_, drStorage_);
    faultOutputManager_->setBackupTimeStamp(backupTimeStamp);
    faultOutputManager_->init();
  }
}

void MemoryManager::initFrictionData() {
  const auto& params = seissolInstance_.getSeisSolParameters().drParameters;
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
  boundaryTree_.synchronizeTo(place, defaultStream);
  surfaceStorage_.synchronizeTo(place, defaultStream);

#ifdef ACL_DEVICE
  device::DeviceInstance::getInstance().api->syncDefaultStreamWithHost();
#endif
}

} // namespace seissol::initializer
