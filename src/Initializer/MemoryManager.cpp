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

#include "Common/Constants.h"
#include "Common/Iterator.h"
#include "DynamicRupture/Factory.h"
#include "DynamicRupture/Misc.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/CellLocalInformation.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/Boundary.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Descriptor/Surface.h"
#include "Memory/GlobalData.h"
#include "Memory/MemoryAllocator.h"
#include "Memory/Tree/Layer.h"
#include "SeisSol.h"

#include <array>
#include <cstddef>
#include <limits>
#include <memory>
#include <utility>
#include <utils/logger.h>
#include <yateto.h>

#ifdef ACL_DEVICE
#include "BatchRecorders/Recorders.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"
#include "Solver/MultipleSimulations.h"

#include <Device/device.h>
#endif // ACL_DEVICE

namespace seissol::initializer {

void MemoryManager::initialize() {
  // initialize global matrices
  GlobalDataInitializerOnHost::init(
      m_globalDataOnHost, m_memoryAllocator, memory::Memkind::Standard);
  if constexpr (seissol::isDeviceOn()) {
    GlobalDataInitializerOnDevice::init(
        m_globalDataOnDevice, m_memoryAllocator, memory::Memkind::DeviceGlobalMemory);
  }
}

void MemoryManager::fixateLtsStorage() {
#ifdef ACL_DEVICE
  MemoryManager::deriveRequiredScratchpadMemoryForDr(drStorage);
  drStorage.allocateScratchPads();
#endif
}

void MemoryManager::fixateBoundaryStorage() {
  const LayerMask ghostMask(Ghost);

  m_boundaryTree.setName("boundary");

  // Boundary face storage
  Boundary::addTo(m_boundaryTree);
  m_boundaryTree.setLayerCount(ltsStorage.getColorMap());
  m_boundaryTree.fixate();

  // Iterate over layers of standard lts storage and face lts storage together.
  for (auto [layer, boundaryLayer] :
       seissol::common::zip(ltsStorage.leaves(ghostMask), m_boundaryTree.leaves(ghostMask))) {
    CellLocalInformation* cellInformation = layer.var<LTS::CellInformation>();

    std::size_t numberOfBoundaryFaces = 0;
    const auto layerSize = layer.size();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+ : numberOfBoundaryFaces)
#endif // _OPENMP
    for (std::size_t cell = 0; cell < layerSize; ++cell) {
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        if (requiresNodalFlux(cellInformation[cell].faceTypes[face])) {
          ++numberOfBoundaryFaces;
        }
      }
    }
    boundaryLayer.setNumberOfCells(numberOfBoundaryFaces);
  }
  m_boundaryTree.allocateVariables();
  m_boundaryTree.touchVariables();

  // The boundary storage is now allocated, now we only need to map from cell lts
  // to face lts.
  // We do this by, once again, iterating over both storages at the same time.
  for (auto [layer, boundaryLayer] :
       seissol::common::zip(ltsStorage.leaves(ghostMask), m_boundaryTree.leaves(ghostMask))) {
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

  surfaceStorage.setName("surface");
  SurfaceLTS::addTo(surfaceStorage);

  int refinement = 0;
  const auto& outputParams = seissolInstance.getSeisSolParameters().output;
  if (outputParams.freeSurfaceParameters.enabled &&
      outputParams.freeSurfaceParameters.vtkorder < 0) {
    refinement = outputParams.freeSurfaceParameters.refinement;
  }
  seissolInstance.freeSurfaceIntegrator().initialize(refinement, ltsStorage, surfaceStorage);
}

#ifdef ACL_DEVICE
void MemoryManager::deriveRequiredScratchpadMemoryForWp(bool plasticity, LTS::Storage& ltsStorage) {
  constexpr size_t totalDerivativesSize = kernels::Solver::DerivativesSize;
  constexpr size_t nodalDisplacementsSize = tensor::averageNormalDisplacement::size();

  for (auto& layer : ltsStorage.leaves(Ghost)) {

    CellLocalInformation* cellInformation = layer.var<LTS::CellInformation>();
    std::unordered_set<real*> registry{};
    real*(*faceNeighbors)[Cell::NumFaces] = layer.var<LTS::FaceNeighborsDevice>();

    std::size_t derivativesCounter{0};
    std::size_t integratedDofsCounter{0};
    std::size_t nodalDisplacementsCounter{0};
    std::size_t analyticCounter = 0;

    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      bool needsScratchMemForDerivatives = !cellInformation[cell].ltsSetup.hasDerivatives();
      if (needsScratchMemForDerivatives) {
        ++derivativesCounter;
      }
      ++integratedDofsCounter;

      // include data provided by ghost layers
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        real* neighborBuffer = faceNeighbors[cell][face];

        // check whether a neighbor element idofs has not been counted twice
        if ((registry.find(neighborBuffer) == registry.end())) {

          // maybe, because of BCs, a pointer can be a nullptr, i.e. skip it
          if (neighborBuffer != nullptr) {
            if (cellInformation[cell].faceTypes[face] != FaceType::Outflow &&
                cellInformation[cell].faceTypes[face] != FaceType::DynamicRupture) {

              const bool isNeighbProvidesDerivatives =
                  cellInformation[cell].ltsSetup.neighborHasDerivatives(face);
              if (isNeighbProvidesDerivatives) {
                ++integratedDofsCounter;
              }
              registry.insert(neighborBuffer);
            }
          }
        }

        if (cellInformation[cell].faceTypes[face] == FaceType::FreeSurfaceGravity) {
          ++nodalDisplacementsCounter;
        }

        if (cellInformation[cell].faceTypes[face] == FaceType::Analytical) {
          ++analyticCounter;
        }
      }
    }

    layer.setEntrySize<LTS::IntegratedDofsScratch>(integratedDofsCounter * tensor::I::size() *
                                                   sizeof(real));
    layer.setEntrySize<LTS::DerivativesScratch>(derivativesCounter * totalDerivativesSize *
                                                sizeof(real));
    layer.setEntrySize<LTS::NodalAvgDisplacements>(nodalDisplacementsCounter *
                                                   nodalDisplacementsSize * sizeof(real));
    layer.setEntrySize<LTS::FSGData>(nodalDisplacementsCounter * 2 * sizeof(real));
#ifdef USE_VISCOELASTIC2
    layer.setEntrySize<LTS::IDofsAneScratch>(layer.size() * tensor::Iane::size() * sizeof(real));
    layer.setEntrySize<LTS::DerivativesExtScratch>(
        layer.size() * (tensor::dQext::size(1) + tensor::dQext::size(2)) * sizeof(real));
    layer.setEntrySize<LTS::DerivativesAneScratch>(
        layer.size() * (tensor::dQane::size(1) + tensor::dQane::size(2)) * sizeof(real));
    layer.setEntrySize<LTS::DofsExtScratch>(layer.size() * tensor::Qext::size() * sizeof(real));
#endif
    layer.setEntrySize<LTS::AnalyticScratch>(analyticCounter * tensor::INodal::size() *
                                             sizeof(real));
    if (plasticity) {
      layer.setEntrySize<LTS::FlagScratch>(layer.size() * sizeof(unsigned));
      layer.setEntrySize<LTS::PrevDofsScratch>(layer.size() * tensor::Q::Size * sizeof(real));
      layer.setEntrySize<LTS::QEtaNodalScratch>(layer.size() * tensor::QEtaNodal::Size *
                                                sizeof(real));
      layer.setEntrySize<LTS::QStressNodalScratch>(layer.size() * tensor::QStressNodal::Size *
                                                   sizeof(real));
    }

#ifdef USE_POROELASTIC
    layer.setEntrySize<LTS::ZinvExtra>(layer.size() * yateto::computeFamilySize<tensor::Zinv>() *
                                       sizeof(real));
#endif
  }
}

void MemoryManager::deriveRequiredScratchpadMemoryForDr(DynamicRupture::Storage& drStorage) {
  constexpr size_t idofsSize = tensor::Q::size() * sizeof(real);
  for (auto& layer : drStorage.leaves()) {
    const auto layerSize = layer.size();
    layer.setEntrySize<DynamicRupture::IdofsPlusOnDevice>(idofsSize * layerSize);
    layer.setEntrySize<DynamicRupture::IdofsMinusOnDevice>(idofsSize * layerSize);
  }
}
#endif

void MemoryManager::initializeMemoryLayout() {
#ifdef ACL_DEVICE
  MemoryManager::deriveRequiredScratchpadMemoryForWp(
      seissolInstance.getSeisSolParameters().model.plasticity, ltsStorage);
  ltsStorage.allocateScratchPads();
#endif
}

void MemoryManager::initializeEasiBoundaryReader(const char* fileName) {
  const auto fileNameStr = std::string{fileName};
  if (!fileNameStr.empty()) {
    m_easiBoundary = EasiBoundary(fileNameStr);
  }
}

#ifdef ACL_DEVICE
void MemoryManager::recordExecutionPaths(bool usePlasticity) {
  recording::CompositeRecorder<LTS::LTSVarmap> recorder;
  recorder.addRecorder(new recording::LocalIntegrationRecorder);
  recorder.addRecorder(new recording::NeighIntegrationRecorder);

  if (usePlasticity) {
    recorder.addRecorder(new recording::PlasticityRecorder);
  }

  for (auto& layer : ltsStorage.leaves(Ghost)) {
    recorder.record(layer);
  }

  recording::CompositeRecorder<DynamicRupture::DynrupVarmap> drRecorder;
  drRecorder.addRecorder(new recording::DynamicRuptureRecorder);
  for (auto& layer : drStorage.leaves(Ghost)) {
    drRecorder.record(layer);
  }
}
#endif // ACL_DEVICE

bool isAcousticSideOfElasticAcousticInterface(CellMaterialData& material, std::size_t face) {
  constexpr auto Eps = std::numeric_limits<real>::epsilon();
  return material.neighbor[face]->getMuBar() > Eps && material.local->getMuBar() < Eps;
}
bool isElasticSideOfElasticAcousticInterface(CellMaterialData& material, std::size_t face) {
  constexpr auto Eps = std::numeric_limits<real>::epsilon();
  return material.local->getMuBar() > Eps && material.neighbor[face]->getMuBar() < Eps;
}

bool isAtElasticAcousticInterface(CellMaterialData& material, std::size_t face) {
  // We define the interface cells as all cells that are in the elastic domain but have a
  // neighbor with acoustic material.
  return material.local != nullptr && material.neighbor[face] != nullptr &&
         (isAcousticSideOfElasticAcousticInterface(material, face) ||
          isElasticSideOfElasticAcousticInterface(material, face));
}

bool requiresDisplacement(CellLocalInformation cellLocalInformation,
                          CellMaterialData& material,
                          std::size_t face) {
  const auto faceType = cellLocalInformation.faceTypes[face];
  return faceType == FaceType::FreeSurface || faceType == FaceType::FreeSurfaceGravity ||
         isAtElasticAcousticInterface(material, face);
}

bool requiresNodalFlux(FaceType f) {
  return (f == FaceType::FreeSurfaceGravity || f == FaceType::Dirichlet ||
          f == FaceType::Analytical);
}

void MemoryManager::initializeFrictionLaw() {
  const auto& params = seissolInstance.getSeisSolParameters().drParameters;
  const auto drParameters = std::make_shared<parameters::DRParameters>(params);
  logInfo() << "Initialize Friction Model";

  logInfo() << "Friction law:" << dr::misc::frictionLawName(drParameters->frictionLawType).c_str()
            << "(" << static_cast<int>(drParameters->frictionLawType) << ")";
  logInfo() << "Thermal pressurization:" << (drParameters->isThermalPressureOn ? "on" : "off");

  const auto factory = seissol::dr::factory::getFactory(drParameters, seissolInstance);
  auto product = factory->produce();
  m_dynRup = std::move(product.storage);
  m_DRInitializer = std::move(product.initializer);
  m_FrictionLaw = std::move(product.frictionLaw);
  m_FrictionLawDevice = std::move(product.frictionLawDevice);
  m_faultOutputManager = std::move(product.output);
}

void MemoryManager::initFaultOutputManager(const std::string& backupTimeStamp) {
  const auto& params = seissolInstance.getSeisSolParameters().drParameters;
  // TODO: switch m_dynRup to shared or weak pointer
  if (params.isDynamicRuptureEnabled) {
    m_faultOutputManager->setInputParam(seissolInstance.meshReader());
    m_faultOutputManager->setLtsData(ltsStorage, backmap, drStorage);
    m_faultOutputManager->setBackupTimeStamp(backupTimeStamp);
    m_faultOutputManager->init();
  }
}

void MemoryManager::initFrictionData() {
  const auto& params = seissolInstance.getSeisSolParameters().drParameters;
  if (params.isDynamicRuptureEnabled) {

    m_DRInitializer->initializeFault(drStorage);
  }
}

void MemoryManager::synchronizeTo(AllocationPlace place) {
#ifdef ACL_DEVICE
  if (place == AllocationPlace::Device) {
    logInfo() << "Synchronizing data... (host->device)";
  } else {
    logInfo() << "Synchronizing data... (device->host)";
  }
  const auto& defaultStream = device::DeviceInstance::getInstance().api->getDefaultStream();
  ltsStorage.synchronizeTo(place, defaultStream);
  drStorage.synchronizeTo(place, defaultStream);
  m_boundaryTree.synchronizeTo(place, defaultStream);
  surfaceStorage.synchronizeTo(place, defaultStream);
  device::DeviceInstance::getInstance().api->syncDefaultStreamWithHost();
#endif
}

} // namespace seissol::initializer
