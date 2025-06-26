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
#include "InternalState.h"
#include "Kernels/Common.h"
#include "Kernels/Touch.h"
#include "Memory/Tree/Layer.h"
#include "SeisSol.h"
#include <DynamicRupture/Factory.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/CellLocalInformation.h>
#include <Initializer/Parameters/DRParameters.h>
#include <Initializer/TimeStepping/ClusterLayout.h>
#include <Initializer/Typedefs.h>
#include <Kernels/Precision.h>
#include <algorithm>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <utils/logger.h>
#include <yateto.h>

#include <DynamicRupture/Misc.h>

#include "generated_code/tensor.h"

#ifdef ACL_DEVICE
#include "BatchRecorders/Recorders.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"
#include "device.h"
#endif // ACL_DEVICE

void seissol::initializer::MemoryManager::initialize() {
  // initialize global matrices
  GlobalDataInitializerOnHost::init(m_globalDataOnHost, m_memoryAllocator, memory::Standard);
  if constexpr (seissol::isDeviceOn()) {
    GlobalDataInitializerOnDevice::init(
        m_globalDataOnDevice, m_memoryAllocator, memory::DeviceGlobalMemory);
  }
}

void seissol::initializer::MemoryManager::fixateLtsTree(struct ClusterLayout& clusterLayout,
                                                        struct MeshStructure* meshStructure,
                                                        unsigned* numberOfDRCopyFaces,
                                                        unsigned* numberOfDRInteriorFaces,
                                                        bool usePlasticity) {
  // store mesh structure and the number of time clusters
  m_meshStructure = meshStructure;

  m_dynRupTree.setName("dr");

  /// Dynamic rupture tree
  m_dynRup->addTo(m_dynRupTree);

  // FIXME: m_dynRupTree.setLayerCount(m_ltsTree.numTimeClusters(), m_ltsTree.getConfigs());
  // FIXME: m_dynRupTree.fixate();
  /* FIXME:
    for (unsigned tc = 0; tc < m_dynRupTree.numTimeClusters(); ++tc) {
      m_dynRupTree.layer(initializer::LayerIdentifier(HaloType::Ghost, Config(), tc))
          .setNumberOfCells(0);
      if (tc < timeStepping.numberOfLocalClusters) {
        m_dynRupTree.layer(initializer::LayerIdentifier(HaloType::Copy, Config(), tc))
            .setNumberOfCells(numberOfDRCopyFaces[tc]);
        m_dynRupTree.layer(initializer::LayerIdentifier(HaloType::Interior, Config(), tc))
            .setNumberOfCells(numberOfDRInteriorFaces[tc]);
      }
    }
  */
  m_dynRupTree.allocateVariables();
  m_dynRupTree.touchVariables();

#ifdef ACL_DEVICE
  MemoryManager::deriveRequiredScratchpadMemoryForDr(m_dynRupTree, *m_dynRup.get());
  m_dynRupTree.allocateScratchPads();
#endif
}

void seissol::initializer::MemoryManager::fixateBoundaryLtsTree() {
  const seissol::initializer::LayerMask ghostMask(Ghost);

  m_boundaryTree.setName("boundary");

  // Boundary face tree
  m_boundary.addTo(m_boundaryTree);
  // FIXME: m_boundaryTree.setLayerCount(m_ltsTree.numTimeClusters(), m_ltsTree.getConfigs());
  // FIXME: m_boundaryTree.fixate();

  // Iterate over layers of standard lts tree and face lts tree together.
  for (auto [layer, boundaryLayer] :
       seissol::common::zip(m_ltsTree.leaves(ghostMask), m_boundaryTree.leaves(ghostMask))) {
    CellLocalInformation* cellInformation = layer.var(m_lts.cellInformation);

    unsigned numberOfBoundaryFaces = 0;
    const auto layerSize = layer.size();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+ : numberOfBoundaryFaces)
#endif // _OPENMP
    for (unsigned cell = 0; cell < layerSize; ++cell) {
      for (unsigned face = 0; face < 4; ++face) {
        if (requiresNodalFlux(cellInformation[cell].faceTypes[face])) {
          ++numberOfBoundaryFaces;
        }
      }
    }
    // FIXME: boundaryLayer.setNumberOfCells(numberOfBoundaryFaces);
  }
  m_boundaryTree.allocateVariables();
  m_boundaryTree.touchVariables();

  // The boundary tree is now allocated, now we only need to map from cell lts
  // to face lts.
  // We do this by, once again, iterating over both trees at the same time.
  for (auto [layer, boundaryLayer] :
       seissol::common::zip(m_ltsTree.leaves(ghostMask), m_boundaryTree.leaves(ghostMask))) {
    auto* cellInformation = layer.var(m_lts.cellInformation);
    auto* boundaryMapping = layer.var(m_lts.boundaryMapping);
    auto* boundaryMappingDevice = layer.var(m_lts.boundaryMappingDevice);
    auto* faceInformation = boundaryLayer.var(m_boundary.faceInformation, AllocationPlace::Host);
    auto* faceInformationDevice =
        boundaryLayer.var(m_boundary.faceInformation, AllocationPlace::Device);

    auto boundaryFace = 0;
    for (unsigned cell = 0; cell < layer.size(); ++cell) {
      for (unsigned face = 0; face < 4; ++face) {
        if (requiresNodalFlux(cellInformation[cell].faceTypes[face])) {
          boundaryMapping[cell][face].nodes = faceInformation[boundaryFace].nodes;
          boundaryMapping[cell][face].dataT = faceInformation[boundaryFace].dataT;
          boundaryMapping[cell][face].dataTinv = faceInformation[boundaryFace].dataTinv;
          boundaryMapping[cell][face].easiBoundaryMap =
              faceInformation[boundaryFace].easiBoundaryMap;
          boundaryMapping[cell][face].easiBoundaryConstant =
              faceInformation[boundaryFace].easiBoundaryConstant;
          boundaryMappingDevice[cell][face].nodes = faceInformationDevice[boundaryFace].nodes;
          boundaryMappingDevice[cell][face].dataT = faceInformationDevice[boundaryFace].dataT;
          boundaryMappingDevice[cell][face].dataTinv = faceInformationDevice[boundaryFace].dataTinv;
          boundaryMappingDevice[cell][face].easiBoundaryMap =
              faceInformationDevice[boundaryFace].easiBoundaryMap;
          boundaryMappingDevice[cell][face].easiBoundaryConstant =
              faceInformationDevice[boundaryFace].easiBoundaryConstant;
          ++boundaryFace;
        } else {
          boundaryMapping[cell][face].nodes = nullptr;
          boundaryMapping[cell][face].dataT = nullptr;
          boundaryMapping[cell][face].dataTinv = nullptr;
          boundaryMapping[cell][face].easiBoundaryMap = nullptr;
          boundaryMapping[cell][face].easiBoundaryConstant = nullptr;
          boundaryMappingDevice[cell][face].nodes = nullptr;
          boundaryMappingDevice[cell][face].dataT = nullptr;
          boundaryMappingDevice[cell][face].dataTinv = nullptr;
          boundaryMappingDevice[cell][face].easiBoundaryMap = nullptr;
          boundaryMappingDevice[cell][face].easiBoundaryConstant = nullptr;
        }
      }
    }
  }
}

#ifdef ACL_DEVICE
void seissol::initializer::MemoryManager::deriveRequiredScratchpadMemoryForWp(bool plasticity,
                                                                              LTSTree& ltsTree,
                                                                              LTS& lts) {
  constexpr size_t totalDerivativesSize = yateto::computeFamilySize<tensor::dQ>();
  constexpr size_t nodalDisplacementsSize = tensor::averageNormalDisplacement::size();

  for (auto& layer : ltsTree.leaves(Ghost)) {

    CellLocalInformation* cellInformation = layer.var(lts.cellInformation);
    std::unordered_set<real*> registry{};
    real*(*faceNeighbors)[4] = layer.var(lts.faceNeighborsDevice);

    std::size_t derivativesCounter{0};
    std::size_t integratedDofsCounter{0};
    std::size_t nodalDisplacementsCounter{0};
    std::size_t analyticCounter = 0;

    std::array<std::size_t, 4> freeSurfacePerFace{};
    std::array<std::size_t, 4> dirichletPerFace{};

    for (unsigned cell = 0; cell < layer.size(); ++cell) {
      bool needsScratchMemForDerivatives = (cellInformation[cell].ltsSetup >> 9) % 2 == 0;
      if (needsScratchMemForDerivatives) {
        ++derivativesCounter;
      }
      ++integratedDofsCounter;

      // include data provided by ghost layers
      for (int face = 0; face < 4; ++face) {
        real* neighborBuffer = faceNeighbors[cell][face];

        // check whether a neighbor element idofs has not been counted twice
        if ((registry.find(neighborBuffer) == registry.end())) {

          // maybe, because of BCs, a pointer can be a nullptr, i.e. skip it
          if (neighborBuffer != nullptr) {
            if (cellInformation[cell].faceTypes[face] != FaceType::Outflow &&
                cellInformation[cell].faceTypes[face] != FaceType::DynamicRupture) {

              bool isNeighbProvidesDerivatives =
                  ((cellInformation[cell].ltsSetup >> face) % 2) == 1;
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
    const auto freeSurfaceCount =
        *std::max_element(freeSurfacePerFace.begin(), freeSurfacePerFace.end());
    const auto dirichletCount = *std::max_element(dirichletPerFace.begin(), dirichletPerFace.end());

    layer.setEntrySize(lts.integratedDofsScratch,
                       integratedDofsCounter * tensor::I::size() * sizeof(real));
    layer.setEntrySize(lts.derivativesScratch,
                       derivativesCounter * totalDerivativesSize * sizeof(real));
    layer.setEntrySize(lts.nodalAvgDisplacements,
                       nodalDisplacementsCounter * nodalDisplacementsSize * sizeof(real));
#ifdef USE_VISCOELASTIC2
    layer.setEntrySize(lts.idofsAneScratch, layer.size() * tensor::Iane::size() * sizeof(real));
    layer.setEntrySize(lts.derivativesExtScratch,
                       layer.size() * (tensor::dQext::size(1) + tensor::dQext::size(2)) *
                           sizeof(real));
    layer.setEntrySize(lts.derivativesAneScratch,
                       layer.size() * (tensor::dQane::size(1) + tensor::dQane::size(2)) *
                           sizeof(real));
    layer.setEntrySize(lts.dofsExtScratch, layer.size() * tensor::Qext::size() * sizeof(real));
#endif
    layer.setEntrySize(lts.analyticScratch,
                       analyticCounter * tensor::INodal::size() * sizeof(real));
    if (plasticity) {
      layer.setEntrySize(lts.flagScratch, layer.size() * sizeof(unsigned));
      layer.setEntrySize(lts.prevDofsScratch, layer.size() * tensor::Q::Size * sizeof(real));
      layer.setEntrySize(lts.qEtaNodalScratch,
                         layer.size() * tensor::QEtaNodal::Size * sizeof(real));
      layer.setEntrySize(lts.qStressNodalScratch,
                         layer.size() * tensor::QStressNodal::Size * sizeof(real));
    }

    layer.setEntrySize(lts.dofsFaceBoundaryNodalScratch,
                       sizeof(real) * dirichletCount * tensor::INodal::size());

    layer.setEntrySize(lts.rotateDisplacementToFaceNormalScratch,
                       sizeof(real) * freeSurfaceCount * init::displacementRotationMatrix::Size);
    layer.setEntrySize(lts.rotateDisplacementToGlobalScratch,
                       sizeof(real) * freeSurfaceCount * init::displacementRotationMatrix::Size);
    layer.setEntrySize(lts.rotatedFaceDisplacementScratch,
                       sizeof(real) * freeSurfaceCount * init::rotatedFaceDisplacement::Size);
    layer.setEntrySize(lts.dofsFaceNodalScratch,
                       sizeof(real) * freeSurfaceCount * tensor::INodal::size());
    layer.setEntrySize(lts.prevCoefficientsScratch,
                       sizeof(real) * freeSurfaceCount * nodal::tensor::nodes2D::Shape[0]);
  }
}

void seissol::initializer::MemoryManager::deriveRequiredScratchpadMemoryForDr(
    LTSTree& ltsTree, DynamicRupture& dynRup) {
  constexpr size_t idofsSize = tensor::Q::size() * sizeof(real);
  for (auto& layer : ltsTree.leaves()) {
    const auto layerSize = layer.size();
    layer.setEntrySize(dynRup.idofsPlusOnDevice, idofsSize * layerSize);
    layer.setEntrySize(dynRup.idofsMinusOnDevice, idofsSize * layerSize);
  }
}
#endif

void seissol::initializer::MemoryManager::initializeMemoryLayout() {
#ifdef ACL_DEVICE
  void* stream = device::DeviceInstance::getInstance().api->getDefaultStream();
  for (auto& layer : m_ltsTree.leaves()) {
    if (layer.getEntrySize(m_lts.buffersDerivatives) > 0) {
      void* data =
          layer.var(m_lts.buffersDerivatives, seissol::initializer::AllocationPlace::Device);
      device::DeviceInstance::getInstance().algorithms.touchMemory(
          reinterpret_cast<real*>(data),
          layer.getEntrySize(m_lts.buffersDerivatives) / sizeof(real),
          true,
          stream);
    }
  }
  device::DeviceInstance::getInstance().api->syncDefaultStreamWithHost();
#endif
  for (auto& layer : m_ltsTree.leaves()) {
    auto* buffers = layer.var(m_lts.buffers);
    auto* derivatives = layer.var(m_lts.derivatives);
    kernels::touchBuffersDerivatives(buffers, derivatives, layer.size());
  }

#ifdef ACL_DEVICE
  seissol::initializer::MemoryManager::deriveRequiredScratchpadMemoryForWp(
      seissolInstance.getSeisSolParameters().model.plasticity, m_ltsTree, m_lts);
  m_ltsTree.allocateScratchPads();
#endif
}

std::pair<MeshStructure*, CompoundGlobalData>
    seissol::initializer::MemoryManager::getMemoryLayout(unsigned int cluster) {
  MeshStructure* meshStructure = m_meshStructure + cluster;

  return std::make_pair(meshStructure, memoryContainer().globalDataStorage);
}

void seissol::initializer::MemoryManager::initializeEasiBoundaryReader(const char* fileName) {
  const auto fileNameStr = std::string{fileName};
  if (!fileNameStr.empty()) {
    m_easiBoundary = EasiBoundary(fileNameStr);
  }
}

#ifdef ACL_DEVICE
void seissol::initializer::MemoryManager::recordExecutionPaths(bool usePlasticity) {
  recording::CompositeRecorder<seissol::initializer::LTS> recorder;
  recorder.addRecorder(new recording::LocalIntegrationRecorder);
  recorder.addRecorder(new recording::NeighIntegrationRecorder);

  if (usePlasticity) {
    recorder.addRecorder(new recording::PlasticityRecorder);
  }

  for (auto& layer : m_ltsTree.leaves(Ghost)) {
    recorder.record(m_lts, layer);
  }

  recording::CompositeRecorder<seissol::initializer::DynamicRupture> drRecorder;
  drRecorder.addRecorder(new recording::DynamicRuptureRecorder);
  for (auto& layer : m_dynRupTree.leaves(Ghost)) {
    drRecorder.record(*m_dynRup, layer);
  }
}
#endif // ACL_DEVICE

bool seissol::initializer::isAcousticSideOfElasticAcousticInterface(
    const CellMaterialData& material, unsigned int face) {
  constexpr auto Eps = std::numeric_limits<real>::epsilon();
  return material.neighbor[face]->getMuBar() > Eps && material.local->getMuBar() < Eps;
}
bool seissol::initializer::isElasticSideOfElasticAcousticInterface(const CellMaterialData& material,
                                                                   unsigned int face) {
  constexpr auto Eps = std::numeric_limits<real>::epsilon();
  return material.local->getMuBar() > Eps && material.neighbor[face]->getMuBar() < Eps;
}

bool seissol::initializer::isAtElasticAcousticInterface(const CellMaterialData& material,
                                                        unsigned int face) {
  // We define the interface cells as all cells that are in the elastic domain but have a
  // neighbor with acoustic material.
  return isAcousticSideOfElasticAcousticInterface(material, face) ||
         isElasticSideOfElasticAcousticInterface(material, face);
}

bool seissol::initializer::requiresDisplacement(CellLocalInformation cellLocalInformation,
                                                const CellMaterialData& material,
                                                unsigned int face) {
  const auto faceType = cellLocalInformation.faceTypes[face];
  return faceType == FaceType::FreeSurface || faceType == FaceType::FreeSurfaceGravity ||
         isAtElasticAcousticInterface(material, face);
}

bool seissol::initializer::requiresNodalFlux(FaceType f) {
  return (f == FaceType::FreeSurfaceGravity || f == FaceType::Dirichlet ||
          f == FaceType::Analytical);
}

void seissol::initializer::MemoryManager::initializeFrictionLaw() {
  const auto drParameters = std::make_shared<seissol::initializer::parameters::DRParameters>(
      m_seissolParams->drParameters);
  logInfo() << "Initialize Friction Model";

  logInfo() << "Friction law:" << dr::misc::frictionLawName(drParameters->frictionLawType).c_str()
            << "(" << static_cast<int>(drParameters->frictionLawType) << ")";
  logInfo() << "Thermal pressurization:" << (drParameters->isThermalPressureOn ? "on" : "off");

  const auto factory = seissol::dr::factory::getFactory(drParameters, seissolInstance);
  auto product = factory->produce();
  m_dynRup = std::move(product.ltsTree);
  m_DRInitializer = std::move(product.initializer);
  m_FrictionLaw = std::move(product.frictionLaw);
  m_FrictionLawDevice = std::move(product.frictionLawDevice);
  m_faultOutputManager = std::move(product.output);
}

void seissol::initializer::MemoryManager::initFaultOutputManager(
    const std::string& backupTimeStamp) {
  // TODO: switch m_dynRup to shared or weak pointer
  if (m_seissolParams->drParameters.isDynamicRuptureEnabled) {
    m_faultOutputManager->setInputParam(seissolInstance.meshReader());
    m_faultOutputManager->setLtsData(&memoryContainer());
    m_faultOutputManager->setBackupTimeStamp(backupTimeStamp);
    m_faultOutputManager->init();
  }
}

void seissol::initializer::MemoryManager::initFrictionData() {
  if (m_seissolParams->drParameters.isDynamicRuptureEnabled) {

    m_DRInitializer->initializeFault(m_dynRup.get(), &m_dynRupTree);

#ifdef ACL_DEVICE
    if (auto* impl = dynamic_cast<dr::friction_law::gpu::FrictionSolverInterface*>(
            m_FrictionLawDevice.get())) {
      impl->allocateAuxiliaryMemory();
    }
#endif // ACL_DEVICE
  }
}

void seissol::initializer::MemoryManager::synchronizeTo(
    seissol::initializer::AllocationPlace place) {
#ifdef ACL_DEVICE
  if (place == seissol::initializer::AllocationPlace::Device) {
    logInfo() << "Synchronizing data... (host->device)";
  } else {
    logInfo() << "Synchronizing data... (device->host)";
  }
  const auto& defaultStream = device::DeviceInstance::getInstance().api->getDefaultStream();
  m_ltsTree.synchronizeTo(place, defaultStream);
  m_dynRupTree.synchronizeTo(place, defaultStream);
  m_boundaryTree.synchronizeTo(place, defaultStream);
  device::DeviceInstance::getInstance().api->syncDefaultStreamWithHost();
#endif
}
