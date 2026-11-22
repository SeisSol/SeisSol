// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "InitModel.h"

#include "Core/Common/BasicTypedefs.h"
#include "Core/Common/ConfigDispatch.h"
#include "Core/Common/ConfigRegistry.h"
#include "Core/Common/ConfigValue.h"
#include "Core/Common/Constants.h"
#include "Core/Data/SimulationFiles.h"
#include "Core/Geometry/CellToVertexArray.h"
#include "Core/Geometry/QueryGenerator.h"
#include "Core/Kernels/Common.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/Tree/Backmap.h"
#include "Core/Memory/Tree/LTSTree.h"
#include "Core/Memory/Tree/Layer.h"
#include "Core/Memory/Typedefs.h"
#include "Core/Modules/Modules.h"
#include "Core/Monitoring/Instrumentation.h"
#include "Core/Monitoring/RunSummary.h"
#include "Core/Monitoring/Stopwatch.h"
#include "Core/Parallel/Helper.h"
#include "Core/TimeStepping/Layout/ClusterLayout.h"
#include "Features/DynamicRupture/Setup/DynamicRuptureMatrices.h"
#include "Features/Plasticity/Physics/Plasticity.h"
#include "Features/Plasticity/Setup/PlasticityPointGenerator.h"
#include "Features/Plasticity/Solver/Storage.h"
#include "Features/TimeMirror/Solver/InstantaneousTimeMirrorManager.h"
#include "GeneratedCode/init.h"
#include "Simulation/InitProcedure/Internal/Boundary.h"
#include "Simulation/InitProcedure/Internal/ConfigBoundaryCheck.h"
#include "Simulation/InitProcedure/Internal/FaceTypeCheck.h"
#include "Simulation/InitProcedure/Internal/Recording.h"
#include "Simulation/InitProcedure/Internal/Scratchpads.h"
#include "Simulation/MemoryManager.h"
#include "Simulation/Parameters/ModelParameters.h"
#include "Simulation/Parameters/SeisSolParameters.h"
#include "Simulation/Reader/ReaderBuilder.h"
#include "Simulation/SeisSol.h"
#include "WavePropagation/Physics/CommonDatastructures.h"
#include "WavePropagation/Physics/Datastructures.h"
#include "WavePropagation/Physics/Energy.h" // IWYU pragma: keep
#include "WavePropagation/Physics/EnergyBase.h"
#include "WavePropagation/Setup/BoundaryMappings.h"
#include "WavePropagation/Setup/CellLocalMatrices.h"
#include "WavePropagation/Setup/DirichletCondition.h"
#include "WavePropagation/Setup/MaterialParameterDB.h"
#include "WavePropagation/WavePropagationParameters.h"

#include <array>
#include <cstddef>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <utils/logger.h>
#include <vector>

namespace seissol::initializer::initprocedure {

namespace {

using seissol::model::Plasticity;

// The material of the cell at `position`, which holds it in its own configuration. A function of
// its own rather than part of the body of the loop over the cells: with the lambda in an OpenMP
// loop, NVHPC 24.9 runs it wrongly with several threads.
model::Material* materialAt(LTS::Storage& storage, const StoragePosition& position) {
  model::Material* material = nullptr;
  const auto config = storage.layer(position.color).getIdentifier().config;
  dispatchConfig(config,
                 [&](auto cfg) { material = &storage.lookup<LTS::MaterialData>(cfg, position); });
  return material;
}

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
  parameterDB.evaluate(*seissol::data::buildReader(fileName, {"x", "y", "z"}), *queryGen);
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

  // The material of every simulation, from each of their files once; a copy of the cells takes the
  // one of the simulations it holds, which all name the same file.
  const parameters::SimulationFiles materialFiles(seissolParams.model.materialFileNames);
  const auto& simulationSettings = memoryManager.simulationSettings();
  const auto materialOfCopy = [&](std::size_t copy) {
    return materialFiles.fileOfSimulation[copy * simulationSettings.fusedCount];
  };
  const auto queryGen = seissol::initializer::getBestQueryGenerator<MaterialT>(
      seissolParams.wavePropagation.useCellHomogenizedMaterial, ctvOfConfig, Cfg::ConvergenceOrder);
  std::vector<std::vector<MaterialT>> materialsDB;
  materialsDB.reserve(materialFiles.files.size());
  for (const auto& materialFile : materialFiles.files) {
    materialsDB.push_back(queryDB<MaterialT>(queryGen, materialFile));
  }

  // plasticity (if needed); a material without the stresses of a solid has none, so its file need
  // not define the plasticity parameters
  const bool withPlasticity = seissolParams.plasticity.enabled && MaterialT::SupportsPlasticity;
  const auto plasticityPointwise = seissolParams.plasticity.pointwise;

  // the plasticity parameters of every simulation of the run, from each of their files once; the
  // fused simulations which only fill up the last copy of the cells repeat the first simulations
  const parameters::SimulationFiles plasticityFiles(seissolParams.model.plasticityFileNames);
  std::vector<std::vector<Plasticity>> plasticityDB(plasticityFiles.files.size());

  if (withPlasticity) {

    // plasticity information is only needed on all interior+copy cells.
    const auto plasticityGen = std::make_shared<PlasticityPointGenerator>(
        ctvOfConfig, plasticityNodes<Cfg>(), plasticityPointwise);
    for (size_t i = 0; i < plasticityDB.size(); i++) {
      plasticityDB[i] = queryDB<Plasticity>(plasticityGen, plasticityFiles.files[i]);
    }
  }

  for (auto& materials : materialsDB) {
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < materials.size(); ++i) {
      materials[i].initialize(seissolParams.wavePropagation);
    }
  }

  logDebug() << "Setting cell materials in the storage (for interior and copy layers).";

  for (auto& layer : memoryManager.ltsStorage().leaves()) {
    if (layer.getIdentifier().config != configIdOf<Cfg>()) {
      continue;
    }

    const auto* cellInformation = layer.var<LTS::CellInformation>();
    const auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();
    auto* materialDataArray = layer.var<LTS::MaterialData>(Cfg());

    if (layer.getIdentifier().halo == HaloType::Ghost) {

#pragma omp parallel for schedule(static)
      for (std::size_t cell = 0; cell < layer.size(); ++cell) {
        const auto& localSecondaryInformation = secondaryInformation[cell];
        const auto meshId = localSecondaryInformation.meshId;
        const auto& linear = meshReader.linearGhostlayer()[meshId];

        // the layer holds materials of exactly this type
        // NOLINTNEXTLINE(bugprone-pointer-arithmetic-on-polymorphic-object)
        auto& materialData = materialDataArray[cell];

        const auto neighborRank = linear.rank;
        const auto neighborRankIdx = linear.inRankIndices[0];
        const auto materialGhostIdx = ghostIdxMap.at(neighborRank)[neighborRankIdx];
        const auto& localMaterial = materialsDB[materialOfCopy(localSecondaryInformation.copyId)]
                                               [position(materialGhostIdx + ghostOffset)];
        initAssign(materialData, localMaterial);
      }
    } else {
      auto* materialArray = layer.var<LTS::Material>();
      auto* plasticityArray =
          withPlasticity ? layer.var<LTSPlasticity::Plasticity>(Cfg()) : nullptr;
      auto* energyDataArray = layer.var<LTS::EnergyData>(Cfg());

#pragma omp parallel for schedule(static)
      for (std::size_t cell = 0; cell < layer.size(); ++cell) {
        // set the materials for the cell volume and its faces
        const auto& localSecondaryInformation = secondaryInformation[cell];
        const auto meshId = localSecondaryInformation.meshId;
        auto& material = materialArray[cell];
        const auto& localMaterial =
            materialsDB[materialOfCopy(localSecondaryInformation.copyId)][position(meshId)];
        const auto& localCellInformation = cellInformation[cell];

        // the layer holds materials of exactly this type
        // NOLINTNEXTLINE(bugprone-pointer-arithmetic-on-polymorphic-object)
        auto& materialData = materialDataArray[cell];
        initAssign(materialData, localMaterial);
        material.local = &materialData;

        energyDataArray[cell] =
            model::EnergyCompute<MaterialT>::template initEnergyData<Cfg>(materialData);

        for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
          if (isInternalFaceType(localCellInformation.faceTypes[side])) {
            // use the neighbor face material info in case that we are not at a boundary
            material.neighbor[side] = materialAt(memoryManager.ltsStorage(),
                                                 localSecondaryInformation.faceNeighbors[side]);
          } else {
            // otherwise, use the material from the own cell
            material.neighbor[side] = material.local;
          }
        }

        // if enabled, set up the plasticity as well
        if (withPlasticity) {
          auto& plasticity = plasticityArray[cell];
          // the copy of the cell holds the simulations from copyId * NumSimulations on
          const auto firstSimulation = localSecondaryInformation.copyId * Cfg::NumSimulations;
          std::array<const Plasticity*, Cfg::NumSimulations> localPlasticity{};
          for (size_t i = 0; i < Cfg::NumSimulations; ++i) {
            const auto pointsPerCell =
                plasticityPointwise ? model::PlasticityData<Cfg>::PointCount : 1;
            const auto simulation = (firstSimulation + i) % plasticityFiles.fileOfSimulation.size();
            localPlasticity[i] = &plasticityDB[plasticityFiles.fileOfSimulation[simulation]]
                                              [position(meshId) * pointsPerCell];
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
  const auto& runConfigs = seissolInstance.parameters().runConfigs;
  const bool singleConfig = runConfigs.configs().size() == 1;
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
        if (singleConfig ||
            runConfigs.configOfGroup(ctv.elementGroups(cell)) == configIdOf<Cfg>()) {
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

  // the Dirichlet conditions of the simulations of the run, from each of their files once
  const parameters::SimulationFiles boundaryFiles(seissolParams.model.boundaryFileNames);
  std::vector<std::optional<DirichletCondition>> dirichletConditions(boundaryFiles.files.size());
  for (std::size_t i = 0; i < boundaryFiles.files.size(); ++i) {
    if (!boundaryFiles.files[i].empty()) {
      dirichletConditions[i] =
          DirichletCondition(seissol::data::buildReader(boundaryFiles.files[i], {"x", "y", "z"}));
    }
  }

  // the boundary mappings carry the Dirichlet map, which the flux solvers absorb
  seissol::initializer::initializeBoundaryMappings(
      meshReader, dirichletConditions, boundaryFiles.fileOfSimulation, memoryManager.ltsStorage());

  seissol::initializer::initializeCellLocalMatrices(meshReader,
                                                    memoryManager.ltsStorage(),
                                                    memoryManager.clusterLayout(),
                                                    seissolParams.wavePropagation);

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
                           seissolParams.plasticity.enabled,
                           seissolInstance.gravitationSetup().acceleration);

  const auto itmParameters = seissolInstance.parameters().timeMirror;

  if (itmParameters.itmEnabled) {
    auto& timeMirrorManagers = seissolInstance.getTimeMirrorManagers();
    const double scalingFactor = itmParameters.itmVelocityScalingFactor;
    const double startingTime = itmParameters.itmStartingTime;

    auto& ltsStorage = memoryManager.ltsStorage();
    const auto* clusterLayout = &seissolInstance.timeManager().getClusterLayout();

    // fixed to one inversion (if you want more, add a loop here)
    auto& increaseManager = timeMirrorManagers.emplace_back(
        seissolInstance, itmParameters, seissolParams.wavePropagation);
    auto& decreaseManager = timeMirrorManagers.emplace_back(
        seissolInstance, itmParameters, seissolParams.wavePropagation);

    initializeTimeMirrorManagers(scalingFactor,
                                 startingTime,
                                 &meshReader,
                                 ltsStorage,
                                 increaseManager,
                                 decreaseManager,
                                 itmParameters,
                                 clusterLayout);
  }
}

void initializeMemoryLayout(seissol::SeisSol& seissolInstance) {

  // set up scratchpads for WP (i.e. mostly for boundary conditions).
  // has to happen after the buckets/buffers are initialized.

  if constexpr (isDeviceOn()) {

    auto& ltsStorage = seissolInstance.memoryManager().ltsStorage();

    seissol::initializer::internal::deriveRequiredScratchpadMemoryForWp(
        seissolInstance.parameters().plasticity.enabled, ltsStorage);
    ltsStorage.allocateScratchPads(useScratchpadPerLayer(seissolInstance.env()) ||
                                           useConcurrentClusters(seissolInstance.env())
                                       ? initializer::ScratchpadSharing::PerLayer
                                       : initializer::ScratchpadSharing::Shared);
  }

  auto& mm = seissolInstance.memoryManager();

  internal::initBoundaryStorage(mm.boundaryStorage(), mm.ltsStorage());
  internal::initSurfaceStorage(
      mm.surfaceStorage(), mm.ltsStorage(), seissolInstance.freeSurfaceIntegrator());
}

} // namespace

void initModel(seissol::SeisSol& seissolInstance) {
  SCOREP_USER_REGION("init_model", SCOREP_USER_REGION_TYPE_FUNCTION);

  logInfo() << "Begin init model.";

  // Call the pre mesh initialization hook
  seissolInstance.modules().callHook<ModuleHook::PreModel>();

  seissol::Stopwatch watch;
  watch.start();

  // these four methods need to be called in this order.
  logInfo() << "Model info:";
  const auto& runConfigs = seissolInstance.parameters().runConfigs;
  const auto configs = runConfigs.configs();
  if (configs.size() == 1) {
    logInfo() << "Configuration:" << configName(configValue(runConfigs.config)).c_str();
  } else {
    // numbered as in the output
    const std::map<int, ConfigId> groupConfigs(runConfigs.groupConfigs.begin(),
                                               runConfigs.groupConfigs.end());
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
  logInfo() << "Plasticity:" << (seissolInstance.parameters().plasticity.enabled ? "on" : "off");
  logInfo() << "Flux:"
            << parameters::fluxToString(seissolInstance.parameters().wavePropagation.flux).c_str();
  logInfo() << "Flux near fault:"
            << parameters::fluxToString(seissolInstance.parameters().wavePropagation.fluxNearFault)
                   .c_str();

  internal::checkFaceTypeSupport(seissolInstance.memoryManager().ltsStorage(),
                                 seissolInstance.parameters().initialization.type);
  internal::checkConfigBoundaries(seissolInstance.memoryManager().ltsStorage());

  // init cell materials (needs LTS, to place the material in; this part was translated from
  // FORTRAN)
  logInfo() << "Initialize cell material parameters.";
  initializeCellMaterial(seissolInstance);

  // init memory layout (needs cell material values to initialize e.g. displacements correctly)
  logInfo() << "Initialize Memory layout.";
  initializeMemoryLayout(seissolInstance);

  // init cell matrices
  logInfo() << "Initialize cell-local matrices.";
  initializeCellMatrices(seissolInstance);

  auto& summary = seissolInstance.runSummary();
  summary.setRankValue("initModel", watch.pause());
  summary["time"]["init"]["model"] = monitoring::toJson(watch.printTime("Model initialized in:"));

  // Call the post mesh initialization hook
  seissolInstance.modules().callHook<ModuleHook::PostModel>();

  logInfo() << "End init model.";
}

} // namespace seissol::initializer::initprocedure
