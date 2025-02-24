// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Equations/Datastructures.h"
#include "Initializer/CellLocalMatrices.h"
#include "Initializer/ParameterDB.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Initializer/TimeStepping/Common.h"
#include "Initializer/Typedefs.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/LTSSync.h"
#include "Memory/Tree/LTSTree.h"
#include "Memory/Tree/Lut.h"
#include <Common/Constants.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/MemoryManager.h>
#include <Initializer/Parameters/ModelParameters.h>
#include <Kernels/Common.h>
#include <Kernels/Precision.h>
#include <Memory/Tree/Layer.h>
#include <Model/CommonDatastructures.h>
#include <Model/Plasticity.h>
#include <Modules/Modules.h>
#include <Monitoring/Stopwatch.h>
#include <Physics/InstantaneousTimeMirrorManager.h>
#include <Solver/Estimator.h>
#include <array>
#include <cassert>
#include <cstddef>
#include <memory>
#include <string>
#include <unordered_map>
#include <utils/env.h>
#include <utils/logger.h>
#include <utils/stringutils.h>
#include <vector>

#include "InitModel.h"
#include "SeisSol.h"

#include "Parallel/MPI.h"

#include <cmath>

#if defined(USE_VISCOELASTIC) || defined(USE_VISCOELASTIC2)
#include "Physics/Attenuation.h"
#endif

using namespace seissol::initializer;

namespace {

using MaterialT = seissol::model::MaterialT;
using Plasticity = seissol::model::Plasticity;

template <typename T>
std::vector<T> queryDB(const std::shared_ptr<seissol::initializer::QueryGenerator>& queryGen,
                       const std::string& fileName,
                       size_t size) {
  std::vector<T> vectorDB(size);
  seissol::initializer::MaterialParameterDB<T> parameterDB;
  parameterDB.setMaterialVector(&vectorDB);
  parameterDB.evaluateModel(fileName, *queryGen);
  return vectorDB;
}

void initializeCellMaterial(seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.getSeisSolParameters();
  const auto& meshReader = seissolInstance.meshReader();
  initializer::MemoryManager& memoryManager = seissolInstance.getMemoryManager();

  // unpack ghost layer (merely a re-ordering operation, since the CellToVertexArray right now
  // requires an vector there)
  std::vector<std::array<std::array<double, 3>, 4>> ghostVertices;
  std::vector<int> ghostGroups;
  std::unordered_map<int, std::vector<unsigned>> ghostIdxMap;
  for (const auto& neighbor : meshReader.getGhostlayerMetadata()) {
    ghostIdxMap[neighbor.first].reserve(neighbor.second.size());
    for (const auto& metadata : neighbor.second) {
      ghostIdxMap[neighbor.first].push_back(ghostVertices.size());
      std::array<std::array<double, 3>, 4> vertices{};
      for (size_t i = 0; i < 4; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          vertices[i][j] = metadata.vertices[i][j];
        }
      }
      ghostVertices.emplace_back(vertices);
      ghostGroups.push_back(metadata.group);
    }
  }

  // just a helper function for better readability
  const auto getBestQueryGenerator = [&](const seissol::initializer::CellToVertexArray& ctvArray) {
    return seissol::initializer::getBestQueryGenerator(
        seissolParams.model.plasticity, seissolParams.model.useCellHomogenizedMaterial, ctvArray);
  };

  // material retrieval for copy+interior layers
  const auto queryGen =
      getBestQueryGenerator(seissol::initializer::CellToVertexArray::fromMeshReader(meshReader));
  auto materialsDB = queryDB<MaterialT>(
      queryGen, seissolParams.model.materialFileName, meshReader.getElements().size());

  // plasticity (if needed)
  std::vector<Plasticity> plasticityDB;
  if (seissolParams.model.plasticity) {
    // plasticity information is only needed on all interior+copy cells.
    plasticityDB = queryDB<Plasticity>(
        queryGen, seissolParams.model.materialFileName, meshReader.getElements().size());
  }

  // material retrieval for ghost layers
  auto queryGenGhost = getBestQueryGenerator(
      seissol::initializer::CellToVertexArray::fromVectors(ghostVertices, ghostGroups));
  auto materialsDBGhost =
      queryDB<MaterialT>(queryGenGhost, seissolParams.model.materialFileName, ghostVertices.size());

#if defined(USE_VISCOELASTIC) || defined(USE_VISCOELASTIC2)
  // we need to compute all model parameters before we can use them...
  // TODO(David): integrate this with the Viscoelastic material class or the ParameterDB directly?
  logDebug() << "Initializing attenuation.";
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (size_t i = 0; i < materialsDB.size(); ++i) {
    auto& cellMat = materialsDB[i];
    seissol::physics::fitAttenuation(
        cellMat, seissolParams.model.freqCentral, seissolParams.model.freqRatio);
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (size_t i = 0; i < materialsDBGhost.size(); ++i) {
    auto& cellMat = materialsDBGhost[i];
    seissol::physics::fitAttenuation(
        cellMat, seissolParams.model.freqCentral, seissolParams.model.freqRatio);
  }
#endif

  logDebug() << "Setting cell materials in the LTS tree (for interior and copy layers).";
  const auto& elements = meshReader.getElements();

  for (auto& layer : memoryManager.getLtsTree()->leaves(Ghost)) {
    auto* cellInformation = layer.var(memoryManager.getLts()->cellInformation);
    auto* secondaryInformation = layer.var(memoryManager.getLts()->secondaryInformation);
    auto* materialArray = layer.var(memoryManager.getLts()->material);
    auto* plasticityArray =
        seissolParams.model.plasticity ? layer.var(memoryManager.getLts()->plasticity) : nullptr;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t cell = 0; cell < layer.getNumberOfCells(); ++cell) {
      // set the materials for the cell volume and its faces
      auto meshId = secondaryInformation[cell].meshId;
      auto& material = materialArray[cell];
      const auto& localMaterial = materialsDB[meshId];
      const auto& element = elements[meshId];
      const auto& localCellInformation = cellInformation[cell];

      initAssign(material.local, localMaterial);
      for (std::size_t side = 0; side < 4; ++side) {
        if (isInternalFaceType(localCellInformation.faceTypes[side])) {
          // use the neighbor face material info in case that we are not at a boundary
          if (element.neighborRanks[side] == seissol::MPI::mpi.rank()) {
            // material from interior or copy
            auto neighbor = element.neighbors[side];
            initAssign(material.neighbor[side], materialsDB[neighbor]);
          } else {
            // material from ghost layer (computed locally)
            auto neighborRank = element.neighborRanks[side];
            auto neighborRankIdx = element.mpiIndices[side];
            auto materialGhostIdx = ghostIdxMap.at(neighborRank)[neighborRankIdx];
            initAssign(material.neighbor[side], materialsDBGhost[materialGhostIdx]);
          }
        } else {
          // otherwise, use the material from the own cell
          initAssign(material.neighbor[side], localMaterial);
        }
      }

      // if enabled, set up the plasticity as well
      if (seissolParams.model.plasticity) {
        auto& plasticity = plasticityArray[cell];
        const auto& localPlasticity = plasticityDB[meshId];

        initAssign(plasticity, seissol::model::PlasticityData(localPlasticity, &material.local));
      }
    }
  }
}

struct LtsInfo {
  unsigned* ltsMeshToFace = nullptr;
  MeshStructure* meshStructure = nullptr;
  TimeStepping timeStepping{};

  // IMPORTANT: DO NOT DEALLOCATE THE ABOVE POINTERS... THEY ARE PASSED ON AND REQUIRED DURING
  // RUNTIME
};

void initializeCellMatrices(LtsInfo& ltsInfo, seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.getSeisSolParameters();

  // \todo Move this to some common initialization place
  auto& meshReader = seissolInstance.meshReader();
  auto& memoryManager = seissolInstance.getMemoryManager();

  seissol::initializer::initializeCellLocalMatrices(meshReader,
                                                    memoryManager.getLtsTree(),
                                                    memoryManager.getLts(),
                                                    memoryManager.getLtsLut(),
                                                    ltsInfo.timeStepping,
                                                    seissolParams.model);

  if (seissolParams.drParameters.etaHack != 1.0) {
    logWarning()
        << "The \"eta hack\" has been enabled to mitigate quasi-divergent solutions in the "
           "friction law. The results may not conform to the existing benchmarks.";
  }

  seissol::initializer::initializeDynamicRuptureMatrices(meshReader,
                                                         memoryManager.getLtsTree(),
                                                         memoryManager.getLts(),
                                                         memoryManager.getLtsLut(),
                                                         memoryManager.getDynamicRuptureTree(),
                                                         memoryManager.getDynamicRupture(),
                                                         ltsInfo.ltsMeshToFace,
                                                         *memoryManager.getGlobalDataOnHost(),
                                                         seissolParams.drParameters.etaHack);

  memoryManager.initFrictionData();

  seissol::initializer::initializeBoundaryMappings(meshReader,
                                                   memoryManager.getEasiBoundaryReader(),
                                                   memoryManager.getLtsTree(),
                                                   memoryManager.getLts(),
                                                   memoryManager.getLtsLut());

#ifdef ACL_DEVICE
  memoryManager.recordExecutionPaths(seissolParams.model.plasticity);
#endif

  auto itmParameters = seissolInstance.getSeisSolParameters().model.itmParameters;

  if (itmParameters.itmEnabled) {
    auto& timeMirrorManagers = seissolInstance.getTimeMirrorManagers();
    const double scalingFactor = itmParameters.itmVelocityScalingFactor;
    const double startingTime = itmParameters.itmStartingTime;

    auto* mLtsTree = memoryManager.getLtsTree();
    auto* mLts = memoryManager.getLts();
    auto* mLtsLut = memoryManager.getLtsLut();
    const auto* mTimeStepping = seissolInstance.timeManager().getTimeStepping();

    initializeTimeMirrorManagers(scalingFactor,
                                 startingTime,
                                 &meshReader,
                                 mLtsTree,
                                 mLts,
                                 mLtsLut,
                                 timeMirrorManagers.first,
                                 timeMirrorManagers.second,
                                 seissolInstance,
                                 mTimeStepping);
  }
}

void hostDeviceCoexecution(seissol::SeisSol& seissolInstance) {
  if constexpr (isDeviceOn()) {
    const auto hdswitch = utils::Env::get<std::string>("SEISSOL_DEVICE_HOST_SWITCH", "none");
    bool hdenabled = false;
    if (hdswitch == "none") {
      hdenabled = false;
      logInfo() << "No host-device switching. Everything runs on the GPU.";
    } else if (hdswitch == "auto") {
      hdenabled = true;
      logInfo() << "Automatic host-device switchpoint detection.";
      const auto hdswitchInt = solver::hostDeviceSwitch();
      seissolInstance.setExecutionPlaceCutoff(hdswitchInt);
    } else {
      hdenabled = true;
      const auto hdswitchInt = utils::StringUtils::parse<int>(hdswitch);
      logInfo() << "Manual host-device cutoff set to" << hdswitchInt << ".";
      seissolInstance.setExecutionPlaceCutoff(hdswitchInt);
    }

#ifdef ACL_DEVICE
    const bool usmDefault = useUSM();
#else
    const bool usmDefault = false;
#endif

    if (!usmDefault && hdenabled) {
      logWarning() << "Using the host-device execution on non-USM systems is not fully supported "
                      "yet. Expect incorrect results.";
    }
  }
}

void initializeClusteredLts(LtsInfo& ltsInfo, seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.getSeisSolParameters();

  assert(seissolParams.timeStepping.lts.getRate() > 0);

  seissolInstance.getLtsLayout().deriveLayout(TimeClustering::MultiRate,
                                              seissolParams.timeStepping.lts.getRate());

  seissolInstance.getLtsLayout().getMeshStructure(ltsInfo.meshStructure);
  seissolInstance.getLtsLayout().getCrossClusterTimeStepping(ltsInfo.timeStepping);

  seissolInstance.getMemoryManager().initializeFrictionLaw();

  unsigned* numberOfDRCopyFaces = nullptr;
  unsigned* numberOfDRInteriorFaces = nullptr;

  seissolInstance.getLtsLayout().getDynamicRuptureInformation(
      ltsInfo.ltsMeshToFace, numberOfDRCopyFaces, numberOfDRInteriorFaces);

  seissolInstance.getMemoryManager().fixateLtsTree(ltsInfo.timeStepping,
                                                   ltsInfo.meshStructure,
                                                   numberOfDRCopyFaces,
                                                   numberOfDRInteriorFaces,
                                                   seissolParams.model.plasticity);

  seissolInstance.getMemoryManager().setLtsToFace(ltsInfo.ltsMeshToFace);

  delete[] numberOfDRCopyFaces;
  delete[] numberOfDRInteriorFaces;

  const auto& ltsTree = seissolInstance.getMemoryManager().getLtsTree();
  const auto& lts = seissolInstance.getMemoryManager().getLts();

  unsigned* ltsToMesh = nullptr;
  unsigned numberOfMeshCells = 0;

  seissolInstance.getLtsLayout().getCellInformation(ltsTree->var(lts->cellInformation),
                                                    ltsTree->var(lts->secondaryInformation),
                                                    ltsToMesh,
                                                    numberOfMeshCells);

  // TODO(David): move all of this method to the MemoryManager
  seissolInstance.getMemoryManager().getLtsLutUnsafe().createLuts(
      ltsTree, ltsToMesh, numberOfMeshCells);

  delete[] ltsToMesh;

  seissol::initializer::time_stepping::deriveLtsSetups(ltsInfo.timeStepping.numberOfLocalClusters,
                                                       ltsInfo.meshStructure,
                                                       ltsTree->var(lts->cellInformation),
                                                       ltsTree->var(lts->secondaryInformation));
}

void initializeMemoryLayout(LtsInfo& ltsInfo, seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.getSeisSolParameters();

  seissolInstance.getMemoryManager().initializeMemoryLayout();

  seissolInstance.timeManager().addClusters(ltsInfo.timeStepping,
                                            ltsInfo.meshStructure,
                                            seissolInstance.getMemoryManager(),
                                            seissolParams.model.plasticity);

  // set tv for all time clusters (this needs to be done, after the time clusters start existing)
  if (seissolParams.model.plasticity) {
    seissolInstance.timeManager().setTv(seissolParams.model.tv);
  }

  seissolInstance.getMemoryManager().fixateBoundaryLtsTree();
}

} // namespace

void seissol::initializer::initprocedure::initModel(seissol::SeisSol& seissolInstance) {
  SCOREP_USER_REGION("init_model", SCOREP_USER_REGION_TYPE_FUNCTION);

  logInfo() << "Begin init model.";

  // Call the pre mesh initialization hook
  seissol::Modules::callHook<ModuleHook::PreModel>();

  seissol::Stopwatch watch;
  watch.start();

  LtsInfo ltsInfo;

  // these four methods need to be called in this order.
  logInfo() << "Model info:";
  logInfo() << "Material:" << MaterialT::Text.c_str();
  logInfo() << "Order:" << ConvergenceOrder;
  logInfo() << "Precision:" << (sizeof(real) == 4 ? "single (f32)" : "double (f64)");
  logInfo() << "Plasticity:"
            << (seissolInstance.getSeisSolParameters().model.plasticity ? "on" : "off");
  logInfo() << "Flux:"
            << parameters::fluxToString(seissolInstance.getSeisSolParameters().model.flux).c_str();
  logInfo() << "Flux near fault:"
            << parameters::fluxToString(seissolInstance.getSeisSolParameters().model.fluxNearFault)
                   .c_str();

  // init LTS
  logInfo() << "Initialize LTS.";
  initializeClusteredLts(ltsInfo, seissolInstance);

  // init cell materials (needs LTS, to place the material in; this part was translated from
  // FORTRAN)
  logInfo() << "Initialize cell material parameters.";
  initializeCellMaterial(seissolInstance);

  if constexpr (isDeviceOn()) {
    logInfo() << "Determine Host-Device switchpoint";
    hostDeviceCoexecution(seissolInstance);
  }

  // init memory layout (needs cell material values to initialize e.g. displacements correctly)
  logInfo() << "Initialize Memory layout.";
  initializeMemoryLayout(ltsInfo, seissolInstance);

  // init cell matrices
  logInfo() << "Initialize cell-local matrices.";
  initializeCellMatrices(ltsInfo, seissolInstance);

  watch.pause();
  watch.printTime("Model initialized in:");

  // Call the post mesh initialization hook
  seissol::Modules::callHook<ModuleHook::PostModel>();

  logInfo() << "End init model.";
}
