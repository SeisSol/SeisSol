// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "InitModel.h"

#include "Common/Constants.h"
#include "Common/Real.h"
#include "Config.h"
#include "Equations/Datastructures.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/CellLocalMatrices.h"
#include "Initializer/MemoryManager.h"
#include "Initializer/ParameterDB.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Initializer/TimeStepping/ClusterLayout.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/LTSTree.h"
#include "Memory/Tree/Layer.h"
#include "Model/CommonDatastructures.h"
#include "Model/Plasticity.h"
#include "Modules/Modules.h"
#include "Monitoring/Stopwatch.h"
#include "Parallel/Helper.h"
#include "Physics/InstantaneousTimeMirrorManager.h"
#include "SeisSol.h"
#include "Solver/Estimator.h"
#include "Solver/MultipleSimulations.h"

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

using namespace seissol::initializer;

namespace {

using namespace seissol;

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
  std::vector<std::array<std::array<double, Cell::Dim>, Cell::NumVertices>> ghostVertices;
  std::vector<int> ghostGroups;
  std::unordered_map<int, std::vector<unsigned>> ghostIdxMap;
  for (const auto& neighbor : meshReader.getGhostlayerMetadata()) {
    ghostIdxMap[neighbor.first].reserve(neighbor.second.size());
    for (const auto& metadata : neighbor.second) {
      ghostIdxMap[neighbor.first].push_back(ghostVertices.size());
      std::array<std::array<double, Cell::Dim>, Cell::NumVertices> vertices{};
      for (size_t i = 0; i < Cell::NumVertices; ++i) {
        for (size_t j = 0; j < Cell::Dim; ++j) {
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
  std::array<std::vector<Plasticity>, seissol::multisim::NumSimulations> plasticityDB;
  if (seissolParams.model.plasticity) {
    // plasticity information is only needed on all interior+copy cells.
    for (size_t i = 0; i < seissol::multisim::NumSimulations; i++) {
      plasticityDB[i] = queryDB<Plasticity>(
          queryGen, seissolParams.model.plasticityFileNames[i], meshReader.getElements().size());
    }
  }

  // material retrieval for ghost layers
  auto queryGenGhost = getBestQueryGenerator(
      seissol::initializer::CellToVertexArray::fromVectors(ghostVertices, ghostGroups));
  auto materialsDBGhost =
      queryDB<MaterialT>(queryGenGhost, seissolParams.model.materialFileName, ghostVertices.size());

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (size_t i = 0; i < materialsDB.size(); ++i) {
    auto& cellMat = materialsDB[i];
    cellMat.initialize(seissolParams.model);
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (size_t i = 0; i < materialsDBGhost.size(); ++i) {
    auto& cellMat = materialsDBGhost[i];
    cellMat.initialize(seissolParams.model);
  }

  logDebug() << "Setting cell materials in the storage (for interior and copy layers).";

  for (auto& layer : memoryManager.getLtsStorage().leaves()) {
    auto* cellInformation = layer.var<LTS::CellInformation>();
    auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();
    auto* materialDataArray = layer.var<LTS::MaterialData>();

    if (layer.getIdentifier().halo == HaloType::Ghost) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (std::size_t cell = 0; cell < layer.size(); ++cell) {
        const auto& localSecondaryInformation = secondaryInformation[cell];
        const auto meshId = localSecondaryInformation.meshId;
        const auto& linear = meshReader.linearGhostlayer()[meshId];

        // explicitly use polymorphic pointer arithmetic here
        // NOLINTNEXTLINE
        auto& materialData = materialDataArray[cell];

        const auto neighborRank = linear.rank;
        const auto neighborRankIdx = linear.inRankIndices[0];
        const auto materialGhostIdx = ghostIdxMap.at(neighborRank)[neighborRankIdx];
        const auto& localMaterial = materialsDBGhost[materialGhostIdx];
        initAssign(materialData, localMaterial);
      }
    } else {
      auto* materialArray = layer.var<LTS::Material>();
      auto* plasticityArray =
          seissolParams.model.plasticity ? layer.var<LTS::Plasticity>() : nullptr;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (std::size_t cell = 0; cell < layer.size(); ++cell) {
        // set the materials for the cell volume and its faces
        const auto& localSecondaryInformation = secondaryInformation[cell];
        const auto meshId = localSecondaryInformation.meshId;
        auto& material = materialArray[cell];
        const auto& localMaterial = materialsDB[meshId];
        const auto& localCellInformation = cellInformation[cell];

        // explicitly use polymorphic pointer arithmetic here
        // NOLINTNEXTLINE
        auto& materialData = materialDataArray[cell];
        initAssign(materialData, localMaterial);
        material.local = &materialData;

        for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
          if (isInternalFaceType(localCellInformation.faceTypes[side])) {
            // use the neighbor face material info in case that we are not at a boundary
            const auto& globalNeighborIndex = localSecondaryInformation.faceNeighbors[side];

            auto* materialNeighbor =
                &memoryManager.getLtsStorage().lookup<LTS::MaterialData>(globalNeighborIndex);
            material.neighbor[side] = materialNeighbor;
          } else {
            // otherwise, use the material from the own cell
            material.neighbor[side] = material.local;
          }
        }

        // if enabled, set up the plasticity as well
        if (seissolParams.model.plasticity) {
          auto& plasticity = plasticityArray[cell];
          assert(plasticityDB.size() == seissol::multisim::NumSimulations &&
                 "Plasticity database size mismatch with number of simulations");
          std::array<Plasticity, seissol::multisim::NumSimulations> localPlasticity;
          for (size_t i = 0; i < seissol::multisim::NumSimulations; ++i) {
            localPlasticity[i] = plasticityDB[i][meshId];
          }
          initAssign(plasticity, seissol::model::PlasticityData(localPlasticity, material.local));
        }
      }
    }
  }
}

void initializeCellMatrices(seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.getSeisSolParameters();

  // \todo Move this to some common initialization place
  auto& meshReader = seissolInstance.meshReader();
  auto& memoryManager = seissolInstance.getMemoryManager();

  seissol::initializer::initializeCellLocalMatrices(meshReader,
                                                    memoryManager.getLtsStorage(),
                                                    memoryManager.clusterLayout(),
                                                    seissolParams.model);

  if (seissolParams.drParameters.etaDamp != 1.0) {
    logWarning() << "The \"eta damp\" (=" << seissolParams.drParameters.etaDamp
                 << ") has been enabled in the timeframe [0,"
                 << seissolParams.drParameters.etaDampEnd
                 << ") to mitigate quasi-divergent solutions in the "
                    "friction law. The results may not conform to the existing benchmarks (which "
                    "are (mostly) computed with \"eta damp\" = 1).";
  }

  seissol::initializer::initializeDynamicRuptureMatrices(meshReader,
                                                         memoryManager.getLtsStorage(),
                                                         memoryManager.getBackmap(),
                                                         memoryManager.getDRStorage(),
                                                         *memoryManager.getGlobalData().onHost,
                                                         seissolParams.drParameters.etaDamp);

  memoryManager.initFrictionData();

  seissol::initializer::initializeBoundaryMappings(
      meshReader, memoryManager.getEasiBoundaryReader(), memoryManager.getLtsStorage());

#ifdef ACL_DEVICE
  memoryManager.recordExecutionPaths(seissolParams.model.plasticity);
#endif

  auto itmParameters = seissolInstance.getSeisSolParameters().model.itmParameters;

  if (itmParameters.itmEnabled) {
    auto& timeMirrorManagers = seissolInstance.getTimeMirrorManagers();
    const double scalingFactor = itmParameters.itmVelocityScalingFactor;
    const double startingTime = itmParameters.itmStartingTime;

    auto& ltsStorage = memoryManager.getLtsStorage();
    const auto* timeStepping = &seissolInstance.timeManager().getClusterLayout();

    initializeTimeMirrorManagers(scalingFactor,
                                 startingTime,
                                 &meshReader,
                                 ltsStorage,
                                 timeMirrorManagers.first,
                                 timeMirrorManagers.second,
                                 seissolInstance,
                                 timeStepping);
  }
}

void hostDeviceCoexecution(seissol::SeisSol& seissolInstance) {
  if constexpr (isDeviceOn()) {
    logInfo() << "Determine Host-Device switchpoint";

    const auto hdswitch = seissolInstance.env().get<std::string>("DEVICE_HOST_SWITCH", "none");
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

    const bool usmDefault = useUSM();

    if (!usmDefault && hdenabled) {
      logWarning() << "Using the host-device execution on non-USM systems is not fully supported "
                      "yet. Expect incorrect results.";
    }
  }
}

void initializeMemoryLayout(seissol::SeisSol& seissolInstance) {
  seissolInstance.getMemoryManager().initializeMemoryLayout();

  seissolInstance.getMemoryManager().fixateBoundaryStorage();
}

} // namespace

void seissol::initializer::initprocedure::initModel(seissol::SeisSol& seissolInstance) {
  SCOREP_USER_REGION("init_model", SCOREP_USER_REGION_TYPE_FUNCTION);

  logInfo() << "Begin init model.";

  // Call the pre mesh initialization hook
  seissol::Modules::callHook<ModuleHook::PreModel>();

  seissol::Stopwatch watch;
  watch.start();

  // these four methods need to be called in this order.
  logInfo() << "Model info:";
  logInfo() << "Material:" << MaterialT::Text.c_str();
  logInfo() << "Order:" << ConvergenceOrder;
  logInfo() << "Precision:"
            << (Config::Precision == RealType::F32 ? "single (f32)" : "double (f64)");
  logInfo() << "Number of simulations: " << Config::NumSimulations;
  logInfo() << "Plasticity:"
            << (seissolInstance.getSeisSolParameters().model.plasticity ? "on" : "off");
  logInfo() << "Flux:"
            << parameters::fluxToString(seissolInstance.getSeisSolParameters().model.flux).c_str();
  logInfo() << "Flux near fault:"
            << parameters::fluxToString(seissolInstance.getSeisSolParameters().model.fluxNearFault)
                   .c_str();

  // init cell materials (needs LTS, to place the material in; this part was translated from
  // FORTRAN)
  logInfo() << "Initialize cell material parameters.";
  initializeCellMaterial(seissolInstance);

  hostDeviceCoexecution(seissolInstance);

  // init memory layout (needs cell material values to initialize e.g. displacements correctly)
  logInfo() << "Initialize Memory layout.";
  initializeMemoryLayout(seissolInstance);

  // init cell matrices
  logInfo() << "Initialize cell-local matrices.";
  initializeCellMatrices(seissolInstance);

  watch.pause();
  watch.printTime("Model initialized in:");

  // Call the post mesh initialization hook
  seissol::Modules::callHook<ModuleHook::PostModel>();

  logInfo() << "End init model.";
}
