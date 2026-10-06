// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "InitModel.h"

#include "Common/ConfigDispatch.h"
#include "Common/ConfigRegistry.h"
#include "Common/ConfigValue.h"
#include "Common/Constants.h"
#include "Equations/Datastructures.h"
#include "Equations/Energy.h" // IWYU pragma: keep
#include "Equations/EnergyBase.h"
#include "GeneratedCode/init.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/InitProcedure/DerivedOutput.h"
#include "Initializer/InitProcedure/Internal/Boundary.h"
#include "Initializer/InitProcedure/Internal/ConfigBoundaryCheck.h"
#include "Initializer/InitProcedure/Internal/FaceTypeCheck.h"
#include "Initializer/InitProcedure/Internal/Recording.h"
#include "Initializer/InitProcedure/Internal/Scratchpads.h"
#include "Initializer/MemoryManager.h"
#include "Initializer/Model/BoundaryMappings.h"
#include "Initializer/Model/CellLocalMatrices.h"
#include "Initializer/Model/DynamicRuptureMatrices.h"
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
#include "Monitoring/Instrumentation.h"
#include "Monitoring/Stopwatch.h"
#include "Parallel/Helper.h"
#include "Physics/InstantaneousTimeMirrorManager.h"
#include "SeisSol.h"
#include "Solver/Estimator.h"

#include <array>
#include <cassert>
#include <cstddef>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <utils/env.h>
#include <utils/logger.h>
#include <utils/stringutils.h>
#include <vector>

namespace seissol::initializer::initprocedure {

namespace {

using Plasticity = seissol::model::Plasticity;

/// The nodes of the plasticity of the configuration `Cfg`, in reference coordinates.
template <typename Cfg>
std::vector<std::array<double, Cell::Dim>> plasticityNodes() {
  const auto nodes = init::vNodes<Cfg>::view::create(init::vNodes<Cfg>::Values);
  std::vector<std::array<double, Cell::Dim>> points(model::PlasticityData<Cfg>::PointCount);
  for (std::size_t i = 0; i < points.size(); ++i) {
    for (std::size_t j = 0; j < Cell::Dim; ++j) {
      if (nodes.isInRange(i, j)) {
        points[i][j] = nodes(i, j);
      }
    }
  }
  return points;
}

template <typename T>
std::vector<T> queryDB(const std::shared_ptr<seissol::initializer::QueryGenerator>& queryGen,
                       const std::string& fileName) {
  std::vector<T> vectorDB;
  seissol::initializer::MaterialParameterDB<T> parameterDB;
  parameterDB.setMaterialVector(&vectorDB);
  parameterDB.evaluateModel(fileName, *queryGen);
  return vectorDB;
}

/// Sets the materials of the cells of the configuration `Cfg`. The cells of the mesh are given by
/// `ctv`, the ghost cells following the inner ones from `ghostOffset` on; `cells` are the ones of
/// the configuration, in increasing order. The material file is queried for them only, since
/// another material need not be defined in the groups of the other cells.
template <typename Cfg>
void initializeCellMaterialOfConfig(
    seissol::SeisSol& seissolInstance,
    const seissol::initializer::CellToVertexArray& ctv,
    const std::vector<std::size_t>& cells,
    std::size_t ghostOffset,
    const std::unordered_map<int, std::vector<unsigned>>& ghostIdxMap) {
  using MaterialT = seissol::model::MaterialOf<Cfg>;
  const auto& seissolParams = seissolInstance.parameters();
  const auto& meshReader = seissolInstance.meshReader();
  initializer::MemoryManager& memoryManager = seissolInstance.memoryManager();

  // the position of a cell of the mesh among the queried ones
  const bool allCells = cells.size() == ctv.size;
  const auto ctvOfConfig =
      allCells ? ctv : seissol::initializer::CellToVertexArray::subset(ctv, cells);
  std::vector<std::size_t> positions;
  if (!allCells) {
    positions.resize(ctv.size);
    for (std::size_t i = 0; i < cells.size(); ++i) {
      positions[cells[i]] = i;
    }
  }
  const auto position = [&](std::size_t cell) { return allCells ? cell : positions[cell]; };

  const auto queryGen = seissol::initializer::getBestQueryGenerator<MaterialT>(
      seissolParams.model.useCellHomogenizedMaterial, ctvOfConfig, Cfg::ConvergenceOrder);
  auto materialsDB = queryDB<MaterialT>(queryGen, seissolParams.model.materialFileName);

  // plasticity (if needed)

  const auto plasticityPointwise = seissolParams.model.plasticityPointwise;

  std::array<std::vector<Plasticity>, Cfg::NumSimulations> plasticityDB;

  if (seissolParams.model.plasticity) {

    // plasticity information is only needed on all interior+copy cells.
    const auto plasticityGen = std::make_shared<PlasticityPointGenerator>(
        ctvOfConfig, plasticityNodes<Cfg>(), plasticityPointwise);
    for (size_t i = 0; i < Cfg::NumSimulations; i++) {
      plasticityDB[i] =
          queryDB<Plasticity>(plasticityGen, seissolParams.model.plasticityFileNames[i]);
    }
  }

#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < materialsDB.size(); ++i) {
    auto& cellMat = materialsDB[i];
    cellMat.initialize(seissolParams.model);
  }

  logDebug() << "Setting cell materials in the storage (for interior and copy layers).";

  for (auto& layer : memoryManager.ltsStorage().leaves()) {
    if (layer.getIdentifier().config != configIdOf<Cfg>()) {
      continue;
    }

    auto* cellInformation = layer.var<LTS::CellInformation>();
    auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();
    auto* materialDataArray = layer.var<LTS::MaterialData>(Cfg());

    if (layer.getIdentifier().halo == HaloType::Ghost) {

#pragma omp parallel for schedule(static)
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
        const auto& localMaterial = materialsDB[position(materialGhostIdx + ghostOffset)];
        initAssign(materialData, localMaterial);
      }
    } else {
      auto* materialArray = layer.var<LTS::Material>();
      auto* plasticityArray =
          seissolParams.model.plasticity ? layer.var<LTS::Plasticity>(Cfg()) : nullptr;
      auto* energyDataArray = layer.var<LTS::EnergyData>(Cfg());

#pragma omp parallel for schedule(static)
      for (std::size_t cell = 0; cell < layer.size(); ++cell) {
        // set the materials for the cell volume and its faces
        const auto& localSecondaryInformation = secondaryInformation[cell];
        const auto meshId = localSecondaryInformation.meshId;
        auto& material = materialArray[cell];
        const auto& localMaterial = materialsDB[position(meshId)];
        const auto& localCellInformation = cellInformation[cell];

        // explicitly use polymorphic pointer arithmetic here
        // NOLINTNEXTLINE
        auto& materialData = materialDataArray[cell];
        initAssign(materialData, localMaterial);
        material.local = &materialData;

        energyDataArray[cell] =
            model::EnergyCompute<MaterialT>::template initEnergyData<Cfg>(materialData);

        for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
          if (isInternalFaceType(localCellInformation.faceTypes[side])) {
            // use the neighbor face material info in case that we are not at a boundary
            const auto& globalNeighborIndex = localSecondaryInformation.faceNeighbors[side];

            // the neighbor holds the material of its own configuration
            auto& storage = memoryManager.ltsStorage();
            const auto neighborConfig =
                storage.layer(globalNeighborIndex.color).getIdentifier().config;
            dispatchConfig(neighborConfig, [&](auto config) {
              material.neighbor[side] =
                  &storage.lookup<LTS::MaterialData>(config, globalNeighborIndex);
            });
          } else {
            // otherwise, use the material from the own cell
            material.neighbor[side] = material.local;
          }
        }

        // if enabled, set up the plasticity as well
        if (seissolParams.model.plasticity) {
          auto& plasticity = plasticityArray[cell];
          assert(plasticityDB.size() == Cfg::NumSimulations &&
                 "Plasticity database size mismatch with number of simulations");
          std::array<const Plasticity*, Cfg::NumSimulations> localPlasticity{};
          for (size_t i = 0; i < Cfg::NumSimulations; ++i) {
            const auto pointsPerCell =
                plasticityPointwise ? model::PlasticityData<Cfg>::PointCount : 1;
            localPlasticity[i] = &plasticityDB[i][position(meshId) * pointsPerCell];
          }
          initAssign(plasticity,
                     seissol::model::PlasticityData<Cfg>(
                         localPlasticity, material.local, plasticityPointwise));
        }
      }
    }
  }
}

void initializeCellMaterial(seissol::SeisSol& seissolInstance) {
  const auto& meshReader = seissolInstance.meshReader();
  initializer::MemoryManager& memoryManager = seissolInstance.memoryManager();

  // unpack ghost layer (merely a re-ordering operation, since the CellToVertexArray right now
  // requires an vector there)
  std::vector<std::array<std::array<double, Cell::Dim>, Cell::NumVertices>> ghostVertices;
  std::vector<int> ghostGroups;
  std::unordered_map<int, std::vector<unsigned>> ghostIdxMap;
  for (const auto& neighbor : meshReader.getGhostlayerMetadata()) {
    ghostIdxMap[neighbor.first].reserve(neighbor.second.size());
    for (const auto& metadata : neighbor.second) {
      ghostIdxMap[neighbor.first].push_back(ghostVertices.size());
      auto& vertices = ghostVertices.emplace_back();
      for (size_t i = 0; i < Cell::NumVertices; ++i) {
        for (size_t j = 0; j < Cell::Dim; ++j) {
          vertices[i][j] = metadata.vertices[i][j];
        }
      }
      ghostGroups.push_back(metadata.group);
    }
  }

  // material retrieval for copy+interior layers
  const auto ctvInner = seissol::initializer::CellToVertexArray::fromMeshReader(meshReader);
  const auto ctvGhost =
      seissol::initializer::CellToVertexArray::fromVectors(ghostVertices, ghostGroups);
  const auto ctv = seissol::initializer::CellToVertexArray::join({ctvInner, ctvGhost});
  const auto ghostOffset = ctvInner.size;

  // every configuration that some cells compute in sets the materials of its cells, which are
  // the ones of its mesh groups
  const auto& model = seissolInstance.parameters().model;
  const bool singleConfig = model.configs().size() == 1;
  forEachConfig([&](auto cfg) {
    using Cfg = decltype(cfg);
    bool hasCells = false;
    for (const auto& layer : memoryManager.ltsStorage().leaves()) {
      hasCells =
          hasCells || (layer.getIdentifier().config == configIdOf<Cfg>() && layer.size() > 0);
    }
    if (hasCells) {
      std::vector<std::size_t> cells;
      for (std::size_t cell = 0; cell < ctv.size; ++cell) {
        if (singleConfig || model.configOfGroup(ctv.elementGroups(cell)) == configIdOf<Cfg>()) {
          cells.push_back(cell);
        }
      }
      initializeCellMaterialOfConfig<Cfg>(seissolInstance, ctv, cells, ghostOffset, ghostIdxMap);
    }
  });
}

void initializeCellMatrices(seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.parameters();

  // \todo Move this to some common initialization place
  auto& meshReader = seissolInstance.meshReader();
  auto& memoryManager = seissolInstance.memoryManager();

  std::optional<DirichletCondition> dirichletCondition;
  if (seissolParams.model.hasBoundaryFile) {
    dirichletCondition = DirichletCondition(seissolParams.model.boundaryFileName);
  }

  // the boundary mappings carry the Dirichlet map, which the flux solvers absorb
  seissol::initializer::initializeBoundaryMappings(
      meshReader, dirichletCondition, memoryManager.ltsStorage());

  seissol::initializer::initializeCellLocalMatrices(
      meshReader, memoryManager.ltsStorage(), memoryManager.clusterLayout(), seissolParams.model);

  if (seissolParams.drParameters.etaDamp != 1.0) {
    logWarning() << "The \"eta damp\" (=" << seissolParams.drParameters.etaDamp
                 << ") has been enabled in the timeframe [0,"
                 << seissolParams.drParameters.etaDampEnd
                 << ") to mitigate quasi-divergent solutions in the "
                    "friction law. The results may not conform to the existing benchmarks (which "
                    "are (mostly) computed with \"eta damp\" = 1).";
  }

  seissol::initializer::initializeDynamicRuptureMatrices(
      meshReader, memoryManager.ltsStorage(), memoryManager.backmap(), memoryManager.drStorage());

  memoryManager.initFrictionData();

  internal::setupRecorders(memoryManager.ltsStorage(),
                           memoryManager.drStorage(),
                           seissolParams.model.plasticity,
                           seissolInstance.gravitationSetup().acceleration);

  auto itmParameters = seissolInstance.parameters().model.itmParameters;

  if (itmParameters.itmEnabled) {
    auto& timeMirrorManagers = seissolInstance.getTimeMirrorManagers();
    const double scalingFactor = itmParameters.itmVelocityScalingFactor;
    const double startingTime = itmParameters.itmStartingTime;

    auto& ltsStorage = memoryManager.ltsStorage();
    const auto* clusterLayout = &seissolInstance.timeManager().getClusterLayout();

    // fixed to one inversion (if you want more, add a loop here)
    auto& increaseManager = timeMirrorManagers.emplace_back(seissolInstance);
    auto& decreaseManager = timeMirrorManagers.emplace_back(seissolInstance);

    initializeTimeMirrorManagers(scalingFactor,
                                 startingTime,
                                 &meshReader,
                                 ltsStorage,
                                 increaseManager,
                                 decreaseManager,
                                 seissolInstance,
                                 clusterLayout);
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

  // set up scratchpads for WP (i.e. mostly for boundary conditions).
  // has to happen after the buckets/buffers are initialized.

  if constexpr (isDeviceOn()) {

    auto& ltsStorage = seissolInstance.memoryManager().ltsStorage();

    seissol::initializer::internal::deriveRequiredScratchpadMemoryForWp(
        seissolInstance.parameters().model.plasticity, ltsStorage);
    ltsStorage.allocateScratchPads();
  }

  auto& mm = seissolInstance.memoryManager();

  internal::initBoundaryStorage(mm.boundaryStorage(), mm.ltsStorage());
  internal::initSurfaceStorage(
      mm.surfaceStorage(),
      mm.ltsStorage(),
      seissolInstance.freeSurfaceIntegrator(),
      derivedStateLayout(seissolInstance, DerivedOutputKind::Surface).size());
}

} // namespace

void initModel(seissol::SeisSol& seissolInstance) {
  SCOREP_USER_REGION("init_model", SCOREP_USER_REGION_TYPE_FUNCTION);

  logInfo() << "Begin init model.";

  // Call the pre mesh initialization hook
  seissol::Modules::callHook<ModuleHook::PreModel>();

  seissol::Stopwatch watch;
  watch.start();

  // these four methods need to be called in this order.
  logInfo() << "Model info:";
  const auto& model = seissolInstance.parameters().model;
  const auto configs = model.configs();
  if (configs.size() == 1) {
    logInfo() << "Configuration:" << configName(configValue(model.config)).c_str();
  } else {
    // numbered as in the output
    const std::map<int, ConfigId> groupConfigs(model.groupConfigs.begin(),
                                               model.groupConfigs.end());
    for (std::size_t i = 0; i < configs.size(); ++i) {
      std::string groups = i == 0 ? "the mesh groups without one of their own" : "the mesh groups";
      for (const auto& [group, config] : groupConfigs) {
        if (config == configs[i]) {
          groups += " " + std::to_string(group);
        }
      }
      logInfo() << "Configuration" << i << ":" << configName(configValue(configs[i])).c_str()
                << "for" << groups.c_str();
    }
  }
  logInfo() << "Plasticity:" << (seissolInstance.parameters().model.plasticity ? "on" : "off");
  logInfo() << "Flux:" << parameters::fluxToString(seissolInstance.parameters().model.flux).c_str();
  logInfo() << "Flux near fault:"
            << parameters::fluxToString(seissolInstance.parameters().model.fluxNearFault).c_str();

  internal::checkFaceTypeSupport(seissolInstance.memoryManager().ltsStorage(),
                                 seissolInstance.parameters().initialization.type);
  internal::checkConfigBoundaries(seissolInstance.memoryManager().ltsStorage());

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

} // namespace seissol::initializer::initprocedure
