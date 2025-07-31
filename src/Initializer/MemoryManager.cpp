// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#include <Initializer/BasicTypedefs.h>
#include <Initializer/Typedefs.h>
#include <Memory/Tree/Colormap.h>
#include <Solver/MultipleSimulations.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "Memory/MemoryAllocator.h"
#include "SeisSol.h"
#include "MemoryManager.h"
#include "InternalState.h"
#include "Memory/Tree/Layer.h"
#include <algorithm>
#include <array>
#include <cstddef>
#include <yateto.h>
#include <unordered_set>
#include <cmath>
#include <type_traits>
#include "Memory/GlobalData.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Kernels/Common.h"
#include "Kernels/Touch.h"

#include <DynamicRupture/Misc.h>

#include "Common/Iterator.h"

#include "generated_code/tensor.h"

#ifdef ACL_DEVICE
#include "BatchRecorders/Recorders.h"
#include "device.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"
#endif // ACL_DEVICE

void seissol::initializer::MemoryManager::initialize()
{
  // initialize global matrices
  GlobalDataInitializerOnHost::init(m_globalDataOnHost, m_memoryAllocator, memory::Standard);
  if constexpr (seissol::isDeviceOn()) {
    GlobalDataInitializerOnDevice::init(m_globalDataOnDevice, m_memoryAllocator, memory::DeviceGlobalMemory);
  }
}

void seissol::initializer::MemoryManager::fixateLtsStorage(struct ClusterLayout& clusterLayout,
                                                         const std::vector<std::size_t>& volumeSizes,
                                                         const std::vector<std::size_t>& drSizes,
                                                         bool usePlasticity) {
  // store mesh structure and the number of time clusters

  ltsStorage.setName("cluster");

  // Setup storage variables
  LTS::addTo(ltsStorage, usePlasticity);
  seissolInstance.postProcessor().allocateMemory(ltsStorage);

  this->clusterLayout = clusterLayout;

  std::vector<std::size_t> clusterMap(clusterLayout.globalClusterCount);
    std::iota(clusterMap.begin(), clusterMap.end(), 0);

  LTSColorMap map(
                 initializer::EnumLayer<HaloType>(
                     {HaloType::Ghost, HaloType::Copy, HaloType::Interior}),
                 initializer::EnumLayer<std::size_t>(clusterMap),
                     initializer::TraitLayer<initializer::ConfigVariant>({Config()}));

  ltsStorage.setLayerCount(map);

  /// From this point, the storage layout, variables, and buckets cannot be changed anymore
  ltsStorage.fixate();

  // Set number of cells and bucket sizes in ltstre
  for (auto [i, layer] : common::enumerate(ltsStorage.leaves())) {
    layer.setNumberOfCells(volumeSizes[i]);
  }

  ltsStorage.allocateVariables();
  ltsStorage.touchVariables();

  drStorage.setName("dr");

  /// Dynamic rupture storage
  m_dynRup->addTo(drStorage);

  drStorage.setLayerCount(ltsStorage.getColorMap());
  drStorage.fixate();

  for (auto [i, layer] : common::enumerate(drStorage.leaves())) {
    layer.setNumberOfCells(drSizes[i]);
  }

  drStorage.allocateVariables();
  drStorage.touchVariables();

#ifdef ACL_DEVICE
  MemoryManager::deriveRequiredScratchpadMemoryForDr(drStorage);
  drStorage.allocateScratchPads();
#endif
}

void seissol::initializer::MemoryManager::fixateBoundaryStorage() {
  seissol::initializer::LayerMask ghostMask(Ghost);

  m_boundaryTree.setName("boundary");

  // Boundary face storage
  Boundary::addTo(m_boundaryTree);
  m_boundaryTree.setLayerCount(ltsStorage.getColorMap());
  m_boundaryTree.fixate();

  // Iterate over layers of standard lts storage and face lts storage together.
  for (auto [layer, boundaryLayer] : seissol::common::zip(ltsStorage.leaves(ghostMask), m_boundaryTree.leaves(ghostMask))) {
    CellLocalInformation* cellInformation = layer.var<LTS::CellInformation>();

    unsigned numberOfBoundaryFaces = 0;
    const auto layerSize = layer.size();
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) reduction(+ : numberOfBoundaryFaces)
#endif // _OPENMP
    for (unsigned cell = 0; cell < layerSize; ++cell) {
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
  for (auto [layer, boundaryLayer] : seissol::common::zip(ltsStorage.leaves(ghostMask), m_boundaryTree.leaves(ghostMask))) {
    auto* cellInformation = layer.var<LTS::CellInformation>();
    auto* boundaryMapping = layer.var<LTS::BoundaryMapping>();
    auto* boundaryMappingDevice = layer.var<LTS::BoundaryMappingDevice>();
    auto* faceInformation = boundaryLayer.var<Boundary::FaceInformation>(AllocationPlace::Host);
    auto* faceInformationDevice = boundaryLayer.var<Boundary::FaceInformation>(AllocationPlace::Device);

    std::size_t boundaryFace = 0;
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        if (requiresNodalFlux(cellInformation[cell].faceTypes[face])) {
          boundaryMapping[cell][face] = CellBoundaryMapping(faceInformation[boundaryFace]);
          boundaryMappingDevice[cell][face] = CellBoundaryMapping(faceInformationDevice[boundaryFace]);
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
  if (outputParams.freeSurfaceParameters.enabled && outputParams.freeSurfaceParameters.vtkorder < 0) {
    refinement = outputParams.freeSurfaceParameters.refinement;
  }
  seissolInstance.freeSurfaceIntegrator().initialize(refinement, &m_globalDataOnHost, ltsStorage, surfaceStorage);
}

#ifdef ACL_DEVICE
void seissol::initializer::MemoryManager::deriveRequiredScratchpadMemoryForWp(bool plasticity, LTS::Storage& ltsStorage) {
  constexpr size_t totalDerivativesSize = yateto::computeFamilySize<tensor::dQ>();
  constexpr size_t nodalDisplacementsSize = tensor::averageNormalDisplacement::size();

  for (auto& layer : ltsStorage.leaves(Ghost)) {

    CellLocalInformation *cellInformation = layer.var<LTS::CellInformation>();
    std::unordered_set<real *> registry{};
    real *(*faceNeighbors)[Cell::NumFaces] = layer.var<LTS::FaceNeighborsDevice>();

    std::size_t derivativesCounter{0};
    std::size_t integratedDofsCounter{0};
    std::size_t nodalDisplacementsCounter{0};
    std::size_t analyticCounter = 0;

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
        real *neighborBuffer = faceNeighbors[cell][face];

        // check whether a neighbor element idofs has not been counted twice
        if ((registry.find(neighborBuffer) == registry.end())) {

          // maybe, because of BCs, a pointer can be a nullptr, i.e. skip it
          if (neighborBuffer != nullptr) {
            if (cellInformation[cell].faceTypes[face] != FaceType::Outflow &&
                cellInformation[cell].faceTypes[face] != FaceType::DynamicRupture) {

              const bool isNeighbProvidesDerivatives = cellInformation[cell].ltsSetup.neighborHasDerivatives(face);
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
      }
    }
    const auto freeSurfaceCount = *std::max_element(freeSurfacePerFace.begin(), freeSurfacePerFace.end());
    const auto dirichletCountPre = *std::max_element(dirichletPerFace.begin(), dirichletPerFace.end());

    // FSG also counts as Dirichlet
    const auto dirichletCount = std::max(dirichletCountPre, freeSurfaceCount);

    layer.setEntrySize<LTS::IntegratedDofsScratch>(
                             integratedDofsCounter * tensor::I::size() * sizeof(real));
    layer.setEntrySize<LTS::DerivativesScratch>(
                             derivativesCounter * totalDerivativesSize * sizeof(real));
    layer.setEntrySize<LTS::NodalAvgDisplacements>(
                             nodalDisplacementsCounter * nodalDisplacementsSize * sizeof(real));
#ifdef USE_VISCOELASTIC2
    layer.setEntrySize<LTS::IDofsAneScratch>(
                             layer.size() * tensor::Iane::size() * sizeof(real));
    layer.setEntrySize<LTS::DerivativesExtScratch>(
                              layer.size() * (tensor::dQext::size(1) + tensor::dQext::size(2)) * sizeof(real));
    layer.setEntrySize<LTS::DerivativesAneScratch>(
                             layer.size() * (tensor::dQane::size(1) + tensor::dQane::size(2)) * sizeof(real));
    layer.setEntrySize<LTS::DofsExtScratch>(
                             layer.size() * tensor::Qext::size() * sizeof(real));
#endif
    layer.setEntrySize<LTS::AnalyticScratch>(
                             analyticCounter * tensor::INodal::size() * sizeof(real));
    if (plasticity) {
      layer.setEntrySize<LTS::FlagScratch>(
                                layer.size() * sizeof(unsigned));
      layer.setEntrySize<LTS::PrevDofsScratch>(
                                layer.size() * tensor::Q::Size * sizeof(real));
      layer.setEntrySize<LTS::QEtaNodalScratch>(
                                layer.size() * tensor::QEtaNodal::Size * sizeof(real));
      layer.setEntrySize<LTS::QStressNodalScratch>(
                                layer.size() * tensor::QStressNodal::Size * sizeof(real));
    }

    layer.setEntrySize<LTS::DofsFaceBoundaryNodalScratch>(sizeof(real) * dirichletCount * tensor::INodal::size());

    layer.setEntrySize<LTS::RotateDisplacementToFaceNormalScratch>(
      sizeof(real) * freeSurfaceCount * init::displacementRotationMatrix::Size);
    layer.setEntrySize<LTS::RotateDisplacementToGlobalScratch>( 
      sizeof(real) * freeSurfaceCount * init::displacementRotationMatrix::Size);
    layer.setEntrySize<LTS::RotatedFaceDisplacementScratch>( 
      sizeof(real) * freeSurfaceCount * init::rotatedFaceDisplacement::Size);
    layer.setEntrySize<LTS::DofsFaceNodalScratch>(
      sizeof(real) * freeSurfaceCount * tensor::INodal::size());
    layer.setEntrySize<LTS::PrevCoefficientsScratch>(
      sizeof(real) * freeSurfaceCount * nodal::tensor::nodes2D::Shape[0]);
  }
}

void seissol::initializer::MemoryManager::deriveRequiredScratchpadMemoryForDr(
    DynamicRupture::Storage& drStorage) {
  constexpr size_t idofsSize = tensor::Q::size() * sizeof(real);
  for (auto& layer : drStorage.leaves()) {
    const auto layerSize = layer.size();
    layer.setEntrySize<DynamicRupture::IdofsPlusOnDevice>(idofsSize * layerSize);
    layer.setEntrySize<DynamicRupture::IdofsMinusOnDevice>(idofsSize * layerSize);
  }
}
#endif

void seissol::initializer::MemoryManager::initializeMemoryLayout()
{
#ifdef ACL_DEVICE
  seissol::initializer::MemoryManager::deriveRequiredScratchpadMemoryForWp(seissolInstance.getSeisSolParameters().model.plasticity, ltsStorage);
  ltsStorage.allocateScratchPads();
#endif
}

void seissol::initializer::MemoryManager::initializeEasiBoundaryReader(const char* fileName) {
  const auto fileNameStr = std::string{fileName};
  if (fileNameStr != "") {
    m_easiBoundary = EasiBoundary(fileNameStr);
  }
}


#ifdef ACL_DEVICE
void seissol::initializer::MemoryManager::recordExecutionPaths(bool usePlasticity) {
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

bool seissol::initializer::isAcousticSideOfElasticAcousticInterface(CellMaterialData &material,
                                              unsigned int face) {
  constexpr auto eps = std::numeric_limits<real>::epsilon();
  return material.neighbor[face].getMuBar() > eps && material.local.getMuBar() < eps;
}
bool seissol::initializer::isElasticSideOfElasticAcousticInterface(CellMaterialData &material,
                                             unsigned int face) {
  constexpr auto eps = std::numeric_limits<real>::epsilon();
  return material.local.getMuBar() > eps && material.neighbor[face].getMuBar() < eps;
}

bool seissol::initializer::isAtElasticAcousticInterface(CellMaterialData &material, unsigned int face) {
  // We define the interface cells as all cells that are in the elastic domain but have a
  // neighbor with acoustic material.
  return isAcousticSideOfElasticAcousticInterface(material, face) || isElasticSideOfElasticAcousticInterface(material, face);
}


bool seissol::initializer::requiresDisplacement(CellLocalInformation cellLocalInformation,
                                                 CellMaterialData &material,
                                                 unsigned int face) {
  const auto faceType = cellLocalInformation.faceTypes[face];
  return faceType == FaceType::FreeSurface
  || faceType == FaceType::FreeSurfaceGravity
  || isAtElasticAcousticInterface(material, face);
}

bool seissol::initializer::requiresNodalFlux(FaceType f) {
  return (f == FaceType::FreeSurfaceGravity
          || f == FaceType::Dirichlet
          || f == FaceType::Analytical);
}

void seissol::initializer::MemoryManager::initializeFrictionLaw() {
  const auto& params = seissolInstance.getSeisSolParameters().drParameters;
  const auto drParameters = std::make_shared<seissol::initializer::parameters::DRParameters>(params);
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

void seissol::initializer::MemoryManager::initFaultOutputManager(const std::string& backupTimeStamp) {
  const auto& params = seissolInstance.getSeisSolParameters().drParameters;
  // TODO: switch m_dynRup to shared or weak pointer
  if (params.isDynamicRuptureEnabled) {
    m_faultOutputManager->setInputParam(seissolInstance.meshReader());
    m_faultOutputManager->setLtsData(ltsStorage,
                                     backmap,
                                     drStorage);
    m_faultOutputManager->setBackupTimeStamp(backupTimeStamp);
    m_faultOutputManager->init();

  }
}


void seissol::initializer::MemoryManager::initFrictionData() {
  const auto& params = seissolInstance.getSeisSolParameters().drParameters;
  if (params.isDynamicRuptureEnabled) {

    m_DRInitializer->initializeFault(drStorage);

#ifdef ACL_DEVICE
    if (auto* impl = dynamic_cast<dr::friction_law::gpu::FrictionSolverInterface*>(m_FrictionLawDevice.get())) {
      impl->allocateAuxiliaryMemory();
    }
#endif // ACL_DEVICE
  }
}

void seissol::initializer::MemoryManager::synchronizeTo(seissol::initializer::AllocationPlace place) {
#ifdef ACL_DEVICE
  if (place == seissol::initializer::AllocationPlace::Device) {
    logInfo() << "Synchronizing data... (host->device)";
  }
  else {
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

