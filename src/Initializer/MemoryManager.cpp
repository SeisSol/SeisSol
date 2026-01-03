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
  GlobalDataInitializerOnHost::init(globalDataOnHost_, memoryAllocator_, memory::Memkind::Standard);
  if constexpr (seissol::isDeviceOn()) {
    GlobalDataInitializerOnDevice::init(
        globalDataOnDevice_, memoryAllocator_, memory::Memkind::DeviceGlobalMemory);
  }
}

void MemoryManager::fixateLtsStorage() {
#ifdef ACL_DEVICE
  MemoryManager::deriveRequiredScratchpadMemoryForDr(drStorage_);
  drStorage_.allocateScratchPads();
#endif
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
    std::size_t numPlasticCells = 0;

    std::array<std::size_t, 4> freeSurfacePerFace{};
    std::array<std::size_t, 4> dirichletPerFace{};

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

        if (cellInformation[cell].faceTypes[face] == FaceType::FreeSurfaceGravity) {
          ++freeSurfacePerFace[face];
        }

        if (cellInformation[cell].faceTypes[face] == FaceType::Dirichlet) {
          ++dirichletPerFace[face];
        }

        if (cellInformation[cell].plasticityEnabled) {
          ++numPlasticCells;
        }
      }
    }
    const auto freeSurfaceCount =
        *std::max_element(freeSurfacePerFace.begin(), freeSurfacePerFace.end());
    const auto dirichletCountPre =
        *std::max_element(dirichletPerFace.begin(), dirichletPerFace.end());

    // FSG also counts as Dirichlet
    const auto dirichletCount = std::max(dirichletCountPre, freeSurfaceCount);

    layer.setEntrySize<LTS::IntegratedDofsScratch>(integratedDofsCounter * tensor::I::size() *
                                                   sizeof(real));
    layer.setEntrySize<LTS::DerivativesScratch>(derivativesCounter * totalDerivativesSize *
                                                sizeof(real));
    layer.setEntrySize<LTS::NodalAvgDisplacements>(nodalDisplacementsCounter *
                                                   nodalDisplacementsSize * sizeof(real));
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
      layer.setEntrySize<LTS::FlagScratch>(numPlasticCells * sizeof(unsigned));
      layer.setEntrySize<LTS::PrevDofsScratch>(numPlasticCells * tensor::Q::Size * sizeof(real));
      layer.setEntrySize<LTS::QEtaNodalScratch>(numPlasticCells * tensor::QEtaNodal::Size *
                                                sizeof(real));
      layer.setEntrySize<LTS::QStressNodalScratch>(numPlasticCells * tensor::QStressNodal::Size *
                                                   sizeof(real));
    }

    layer.setEntrySize<LTS::DofsFaceBoundaryNodalScratch>(sizeof(real) * dirichletCount *
                                                          tensor::INodal::size());

    layer.setEntrySize<LTS::RotateDisplacementToFaceNormalScratch>(
        sizeof(real) * freeSurfaceCount * init::displacementRotationMatrix::Size);
    layer.setEntrySize<LTS::RotateDisplacementToGlobalScratch>(
        sizeof(real) * freeSurfaceCount * init::displacementRotationMatrix::Size);
    layer.setEntrySize<LTS::RotatedFaceDisplacementScratch>(sizeof(real) * freeSurfaceCount *
                                                            init::rotatedFaceDisplacement::Size);
    layer.setEntrySize<LTS::DofsFaceNodalScratch>(sizeof(real) * freeSurfaceCount *
                                                  tensor::INodal::size());
    layer.setEntrySize<LTS::PrevCoefficientsScratch>(
        sizeof(real) * freeSurfaceCount *
        nodal::tensor::nodes2D::Shape[multisim::BasisFunctionDimension]);
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
      seissolInstance_.getSeisSolParameters().model.plasticity, ltsStorage_);
  ltsStorage_.allocateScratchPads();
#endif
}

void MemoryManager::initializeEasiBoundaryReader(const char* fileName) {
  const auto fileNameStr = std::string{fileName};
  if (!fileNameStr.empty()) {
    easiBoundary_ = EasiBoundary(fileNameStr);
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

  for (auto& layer : ltsStorage_.leaves(Ghost)) {
    recorder.record(layer);
  }

  recording::CompositeRecorder<DynamicRupture::DynrupVarmap> drRecorder;
  drRecorder.addRecorder(new recording::DynamicRuptureRecorder);
  for (auto& layer : drStorage_.leaves(Ghost)) {
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
  const auto& params = seissolInstance_.getSeisSolParameters().drParameters;
  const auto drParameters = std::make_shared<parameters::DRParameters>(params);
  logInfo() << "Initialize Friction Model";

  logInfo() << "Friction law:" << dr::misc::frictionLawName(drParameters->frictionLawType).c_str()
            << "(" << static_cast<int>(drParameters->frictionLawType) << ")";
  logInfo() << "Thermal pressurization:" << (drParameters->isThermalPressureOn ? "on" : "off");

  const auto factory = seissol::dr::factory::getFactory(drParameters, seissolInstance_);
  auto product = factory->produce();
  dynRup_ = std::move(product.storage);
  DRInitializer_ = std::move(product.initializer);
  FrictionLaw_ = std::move(product.frictionLaw);
  FrictionLawDevice_ = std::move(product.frictionLawDevice);
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

    DRInitializer_->initializeFault(drStorage_);
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
  ltsStorage_.synchronizeTo(place, defaultStream);
  drStorage_.synchronizeTo(place, defaultStream);
  boundaryTree_.synchronizeTo(place, defaultStream);
  surfaceStorage_.synchronizeTo(place, defaultStream);
  device::DeviceInstance::getInstance().api->syncDefaultStreamWithHost();
#endif
}

} // namespace seissol::initializer
