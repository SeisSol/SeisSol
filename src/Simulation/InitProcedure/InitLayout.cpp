// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "InitLayout.h"

#include "Core/Common/BasicTypedefs.h"
#include "Core/Common/ConfigDispatch.h"
#include "Core/Common/ConfigRegistry.h"
#include "Core/Common/ConfigValue.h"
#include "Core/Common/Constants.h"
#include "Core/Common/Iterator.h"
#include "Core/Common/Json.h"
#include "Core/Geometry/MeshReader.h"
#include "Core/Kernels/Common.h"
#include "Core/Memory/Descriptor/DynamicRupture.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/Tree/Backmap.h"
#include "Core/Memory/Tree/Colormap.h"
#include "Core/Memory/Tree/LTSTree.h"
#include "Core/Memory/Tree/Layer.h"
#include "Core/Parallel/Helper.h"
#include "Core/Parallel/MPI.h"
#include "Core/TimeStepping/Compute/FaceClusterFactory.h"
#include "Core/TimeStepping/Layout/ClusterLayout.h"
#include "Core/TimeStepping/Layout/Halo.h"
#include "Core/TimeStepping/Settings.h"
#include "Features/DynamicRupture/Solver/DynamicRuptureFeature.h"
#include "Features/Plasticity/Solver/PlasticityFeature.h"
#include "Internal/MeshLayout.h"
#include "Simulation/Estimator.h"
#include "Simulation/InitProcedure/Internal/Buckets.h"
#include "Simulation/InitProcedure/Internal/LtsSetup.h"
#include "Simulation/InitProcedure/Internal/Scratchpads.h"
#include "Simulation/Parameters/OutputParameters.h"
#include "Simulation/SeisSol.h"
#include "WavePropagation/Output/IntegrationFeature.h"
#include "WavePropagation/Physics/Datastructures.h"
#include "WavePropagation/Solver/WaveClusterFactory.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <mpi.h>
#include <numeric>
#include <string>
#include <unordered_map>
#include <utils/logger.h>
#include <utils/stringutils.h>
#include <vector>

namespace seissol::initializer::initprocedure {

namespace {

void verifyHaloSetup(const LTS::Storage& ltsStorage, const std::vector<ClusterMap>& meshLayout) {
  // just verify everything. I.e. exchange the global IDs of copy and ghost layers, and see if they
  // match (or not).

  std::vector<std::size_t> copyGlobalIds(ltsStorage.size(Ghost | Interior));
  std::vector<std::size_t> ghostGlobalIds(ltsStorage.size(Copy | Interior));
  std::vector<MPI_Request> requests;

  std::size_t copyI = 0;
  std::size_t ghostI = 0;
  for (const auto& layer : ltsStorage.leaves(Interior)) {
    const auto& layout = meshLayout[layer.id()];
    std::size_t copyIL = 0;
    for (const auto& region : layout.regions) {
      auto& request = requests.emplace_back(MPI_REQUEST_NULL);
      if (layer.getIdentifier().halo == HaloType::Copy) {
        for (std::size_t i = 0; i < region.count; ++i) {
          copyGlobalIds[copyI + i] = layer.var<LTS::SecondaryInformation>()[copyIL + i].globalId;
        }
        MPI_Isend(&copyGlobalIds[copyI],
                  region.count,
                  seissol::Mpi::castToMpiType<std::size_t>(),
                  region.rank,
                  region.tag,
                  seissol::Mpi::mpi.comm(),
                  &request);
        copyI += region.count;
        copyIL += region.count;
      }
      if (layer.getIdentifier().halo == HaloType::Ghost) {
        MPI_Irecv(&ghostGlobalIds[ghostI],
                  region.count,
                  seissol::Mpi::castToMpiType<std::size_t>(),
                  region.rank,
                  region.tag,
                  seissol::Mpi::mpi.comm(),
                  &request);
        ghostI += region.count;
      }
    }
  }

  MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

  std::size_t ghostC = 0;
  for (const auto& layer : ltsStorage.leaves(Copy | Interior)) {
    for (std::size_t i = 0; i < layer.size(); ++i) {
      if (layer.var<LTS::SecondaryInformation>()[i].globalId != ghostGlobalIds[ghostC]) {
        logError() << "Internal error: halo setup. Global IDs"
                   << layer.var<LTS::SecondaryInformation>()[i].globalId << "vs"
                   << ghostGlobalIds[ghostC];
      }
      ++ghostC;
    }
  }
}

//! @brief The clusters of the run, for its summary.
json::Value ltsSection(const ClusterLayout& layout) {
  json::Value section;
  section["minimumTimestep"] = layout.minimumTimestep;
  section["rates"] = json::Value::arrayOf(layout.rates);
  if (layout.report.has_value()) {
    section["minimumCellTimestep"] = layout.report->minimumCellTimestep;
    section["wiggleFactor"] = layout.report->wiggleFactor;
    section["theoreticalSpeedup"]["elementwise"] = layout.report->elementwiseSpeedup;
    section["theoreticalSpeedup"]["clustered"] = layout.report->clusteredSpeedup;
  }
  auto& clusters = section["clusters"] = json::Value::array();
  for (std::size_t id = 0; id < layout.globalClusterCount; ++id) {
    auto& cluster = clusters.push(json::Value::object());
    cluster["id"] = id;
    cluster["timestep"] = layout.timestepRate(id);
    cluster["updateFactor"] = layout.clusterRate(id);
    if (layout.report.has_value()) {
      cluster["cells"] = layout.report->cells.at(id);
      cluster["dynamicRuptureFaces"] = layout.report->dynamicRuptureFaces.at(id);
    }
  }
  return section;
}

void setupMemory(seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.parameters();
  const auto& meshReader = seissolInstance.meshReader();

  const auto rank = seissol::Mpi::mpi.rank();

  logInfo() << "Determining cell colors...";

  // Mesh readers that never construct LtsWeights (e.g. the cube generator) leave the channel
  // unset; those produce a single cluster, so falling back to the parameter vector is safe.
  const auto effectiveRates = seissolInstance.getEffectiveClusterRates().empty()
                                  ? seissolParams.timeStepping.lts.getRate()
                                  : seissolInstance.getEffectiveClusterRates();

  const auto clusterLayout = ClusterLayout::fromMesh(
      effectiveRates, seissolInstance.meshReader(), seissolInstance.timestepScale(), true);

  seissolInstance.memoryManager().setClusterLayout(clusterLayout);
  seissolInstance.runSummary()["lts"] = ltsSection(clusterLayout);

  std::vector<std::size_t> clusterMap(clusterLayout.globalClusterCount);
  std::iota(clusterMap.begin(), clusterMap.end(), 0);

  // every cell computes in the configuration of its mesh group, ghost cells included
  const auto& runConfigs = seissolParams.runConfigs;

  const LTSColorMap colorMap(
      initializer::EnumLayer<HaloType>({HaloType::Ghost, HaloType::Copy, HaloType::Interior}),
      initializer::EnumLayer<std::size_t>(clusterMap),
      initializer::EnumLayer<ConfigId>(runConfigs.configs()));

  std::vector<std::size_t> colors(meshReader.getElements().size());
  for (std::size_t i = 0; i < colors.size(); ++i) {
    const auto& element = meshReader.getElements()[i];
    const auto halo = geometry::isCopy(element, rank) ? HaloType::Copy : HaloType::Interior;
    colors[i] = colorMap.color(halo, element.clusterId, runConfigs.configOfGroup(element.group));
  }

  const auto ghostSize = meshReader.linearGhostlayer().size();

  std::vector<std::size_t> colorsGhost(ghostSize);
  for (const auto [i, linearGhost] : common::enumerate(meshReader.linearGhostlayer())) {
    const auto& element =
        meshReader.getGhostlayerMetadata().at(linearGhost.rank)[linearGhost.inRankIndices[0]];
    const auto halo = HaloType::Ghost;
    colorsGhost[i] =
        colorMap.color(halo, element.clusterId, runConfigs.configOfGroup(element.group));
  }

  // the time integrals of the cells, for the wave field outputs that write them or the strain
  const auto waveFields =
      seissolParams.output.outputsOf<seissol::initializer::parameters::WaveFieldOutputParameters>();
  const auto needsIntegration = std::ranges::any_of(waveFields, [](const auto* waveField) {
    return waveField->computeStrain ||
           std::ranges::any_of(waveField->integrationMask, [](const auto& mask) {
             return std::ranges::find(mask, true) != mask.end();
           });
  });
  const auto settings = SimulationSettings(seissolParams.plasticity.enabled,
                                           needsIntegration,
                                           runConfigs.simulationCount,
                                           configValue(runConfigs.config).numSimulations);

  seissolInstance.memoryManager().setSimulationSettings(settings);

  // if we need more simulations than the cells fuse, use copies of the cells for the others
  const auto copyCount = settings.copyCount();
  auto& simulations = seissolInstance.runSummary()["simulations"];
  simulations["count"] = settings.simulationCount;
  simulations["fused"] = settings.fusedCount;
  simulations["copies"] = copyCount;
  logInfo() << "Simulations:" << settings.simulationCount << "in" << copyCount
            << "copies of the cells, which fuse" << settings.fusedCount << "simulations each.";
  if (copyCount * settings.fusedCount > settings.simulationCount) {
    logInfo() << "The last" << copyCount * settings.fusedCount - settings.simulationCount
              << "fused simulations only fill up the last copy; they are not written.";
  }

  logInfo() << "Creating mesh layout...";

  const auto meshLayout = internal::layoutCells(colors, colorsGhost, colorMap, meshReader);
  auto meshLayoutCopy = meshLayout;
  for (auto& cluster : meshLayoutCopy) {
    for (auto& region : cluster.regions) {
      region.count *= copyCount;
    }
  }

  // the features of the run, in the order their work runs in: e.g. the plasticity before the time
  // integral of the DOFs for the outputs, which integrates the DOFs after the plasticity
  auto& memoryManager = seissolInstance.memoryManager();
  auto& features = seissolInstance.features();
  features.clear();
  features.push_back(std::make_unique<seissol::kernels::PlasticityFeature>(
      seissolInstance, seissolParams.plasticity, memoryManager.globalDataStore()));
  features.push_back(std::make_unique<seissol::kernels::IntegrationFeature>(settings.integrate));

  auto& ltsStorage = memoryManager.ltsStorage();
  auto& backmap = memoryManager.backmap();
  LTS::addTo(ltsStorage, settings);
  for (auto& feature : features) {
    feature->addStorage(ltsStorage);
  }
  ltsStorage.setName("cluster");
  ltsStorage.setLayerCount(colorMap);
  // for the outputs on a thread of their own, which read the data while the device goes on
  ltsStorage.setSnapshots(keepOutputSnapshots(seissolInstance.env()));
  ltsStorage.fixate();
  for (auto& layer : ltsStorage.leaves()) {
    layer.setNumberOfCells(meshLayout[layer.id()].cellMap.size() * copyCount);
  }
  ltsStorage.allocateVariables();
  ltsStorage.touchVariables();
  backmap.setSize(meshReader.getElements().size(), copyCount);

  StorageBackmap<Cell::NumFaces> backmapGhost;
  backmapGhost.setSize(ghostSize, copyCount);
  std::vector<std::vector<std::unordered_map<std::size_t, StoragePosition>>> backmapGhostMap(
      copyCount, std::vector<std::unordered_map<std::size_t, StoragePosition>>(ghostSize));

  logInfo() << "Setting up cell storage...";
  // just need pick a test variable that is available everywhere
  const auto* zero = ltsStorage.var<LTS::SecondaryInformation>();
  for (auto& layer : ltsStorage.leaves()) {
    auto* zeroLayer = layer.var<LTS::SecondaryInformation>();
    const auto& layout = meshLayout[layer.id()];
    auto regionIt = layout.regions.begin();
    std::size_t regionOffset = 0;

    const auto addToBackmapCopyInterior = [&](auto cell, auto index, auto copy) {
      return backmap.addElement(layer.id(), zero, zeroLayer, cell, copy, index);
    };
    const auto addToBackmapGhost = [&](auto cell, auto index, auto copy) {
      assert(regionIt != layout.regions.end());
      const auto dup = backmapGhost.addElement(layer.id(), zero, zeroLayer, cell, copy, index);
      backmapGhostMap.at(copy).at(cell)[regionIt->remoteId] =
          backmapGhost.getDup(cell, copy, dup).value();
      return dup;
    };
    const auto& addToBackmap = [&](auto cell, auto index, auto copy) {
      if (layer.getIdentifier().halo == HaloType::Ghost) {
        return addToBackmapGhost(cell, index, copy);
      } else {
        return addToBackmapCopyInterior(cell, index, copy);
      }
    };

    for (const auto [i, cell] : common::enumerate(layout.cellMap)) {
      while (regionIt != layout.regions.end() && i >= regionOffset + regionIt->count) {
        regionOffset += regionIt->count;
        ++regionIt;
      }

      for (std::size_t copy = 0; copy < copyCount; ++copy) {
        const auto idx = i * copyCount + copy;

        const auto dup = addToBackmap(cell, idx, copy);

        zeroLayer[idx].meshId = cell;
        zeroLayer[idx].copyId = copy;
        zeroLayer[idx].configId = layer.getIdentifier().config;
        zeroLayer[idx].clusterId = layer.getIdentifier().lts;
        zeroLayer[idx].duplicate = dup;
        zeroLayer[idx].halo = layer.getIdentifier().halo;
        zeroLayer[idx].rank = -1;
        zeroLayer[idx].color = layer.id();
        for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
          zeroLayer[idx].neighborRanks[face] = -1;
          zeroLayer[idx].faceNeighbors[face] = StoragePosition::NullPosition;
        }
      }
    }
  }

  logInfo() << "Setting up cell information data...";
  for (auto& layer : ltsStorage.leaves()) {
    const auto& layout = meshLayout[layer.id()];
    auto* secondaryCellInformation = layer.var<LTS::SecondaryInformation>();
    auto* cellInformation = layer.var<LTS::CellInformation>();
    const auto& cells = layout.cellMap;
    if (layer.getIdentifier().halo == HaloType::Interior ||
        layer.getIdentifier().halo == HaloType::Copy) {

#pragma omp parallel for schedule(static)
      for (std::size_t i = 0; i < layer.size(); ++i) {
        const auto cell = cells[i / copyCount];
        const auto index = i;
        const auto& element = meshReader.getElements()[cell];

        const auto copy = secondaryCellInformation[i].copyId;

        secondaryCellInformation[index].rank = rank;
        secondaryCellInformation[index].globalId = element.globalId;
        secondaryCellInformation[index].group = element.group;
        for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
          secondaryCellInformation[index].neighborRanks[face] = element.neighborRanks[face];
          cellInformation[index].faceTypes[face] = element.boundaries[face];
          cellInformation[index].faceRelations[face][0] = element.neighborSides[face];
          cellInformation[index].faceRelations[face][1] = element.sideOrientations[face];

          if (element.neighbors[face].hasValue() || element.neighborRanks[face] != rank) {
            const auto& neighbor = [&] {
              const bool ghostNeighbor = element.neighborRanks[face] != rank;
              if (ghostNeighbor) {
                const auto rank = element.neighborRanks[face];
                const auto mpiIndex = element.mpiIndices[face];
                const auto linear = meshReader.toLinearGhostlayer().at({rank, mpiIndex});
                // ghost layer
                return backmapGhostMap.at(copy).at(linear).at(layer.id());
              } else {
                // copy/interior layer
                assert(element.neighbors[face].hasValue());
                return backmap.get(element.neighbors[face].value(), copy);
              }
            }();

            secondaryCellInformation[index].faceNeighbors[face] = neighbor;
            cellInformation[index].neighborConfigIds[face] =
                ltsStorage.layer(neighbor.color).getIdentifier().config;
          } else {
            secondaryCellInformation[index].faceNeighbors[face] = StoragePosition::NullPosition;
            cellInformation[index].neighborConfigIds[face] =
                std::numeric_limits<std::uint32_t>::max();
          }
        }
      }
    } else if (layer.getIdentifier().halo == HaloType::Ghost) {

#pragma omp parallel for schedule(static)
      for (std::size_t i = 0; i < layer.size(); ++i) {
        const auto cell = cells[i / copyCount];

        const auto& linear = meshReader.linearGhostlayer()[cell];
        secondaryCellInformation[i].rank = linear.rank;

        const auto copy = secondaryCellInformation[i].copyId;

        for (const auto& index : linear.inRankIndices) {
          const auto& boundaryElement =
              meshReader.getMPINeighbors().at(linear.rank).elements[index];
          const auto neighbor = backmap.get(boundaryElement.localElement, copy);
          const auto& elementNeighbor = meshReader.getElements()[boundaryElement.localElement];

          const auto& metadata = meshReader.getGhostlayerMetadata().at(linear.rank)[index];
          secondaryCellInformation[i].globalId = metadata.globalId;
          secondaryCellInformation[i].group = metadata.group;

          // equivalent to elementNeighbor.neighborSides[boundaryElement.localSide]
          const auto face = boundaryElement.neighborSide;

          secondaryCellInformation[i].faceNeighbors[face] = neighbor;

          secondaryCellInformation[i].neighborRanks[face] = rank;
          cellInformation[i].neighborConfigIds[face] =
              ltsStorage.layer(neighbor.color).getIdentifier().config;
          cellInformation[i].faceTypes[face] =
              elementNeighbor.boundaries[boundaryElement.localSide];
          cellInformation[i].faceRelations[face][0] = boundaryElement.localSide;
          cellInformation[i].faceRelations[face][1] =
              elementNeighbor.sideOrientations[boundaryElement.localSide];
        }
      }
    }
  }

  logInfo() << "Verify the halo setup...";
  verifyHaloSetup(ltsStorage, meshLayoutCopy);

  // pass 3: LTS setup
  logInfo() << "Setting up LTS configuration...";
  internal::deriveLtsSetups(meshLayoutCopy, ltsStorage);

  if (seissolParams.plasticity.enabled) {
    // remove disabled plasticity groups from the list, and the cells of a material without the
    // stresses of a solid, where nothing yields
    const auto& pdis = seissolParams.plasticity.disabledGroups;
    for (auto& layer : ltsStorage.leaves(Ghost)) {
      const std::size_t size = layer.size();
      auto* cellInfo = layer.var<LTS::CellInformation>();
      const auto* cellInfo2 = layer.var<LTS::SecondaryInformation>();
      const bool supported = dispatchConfig(layer.getIdentifier().config, [](auto config) {
        return model::MaterialOf<decltype(config)>::SupportsPlasticity;
      });

#pragma omp parallel for schedule(static)
      for (std::size_t i = 0; i < size; ++i) {
        cellInfo[i].plasticityEnabled = supported && !pdis.contains(cellInfo2[i].group);
      }
    }
  }

  seissolInstance.memoryManager().initializeFrictionLaw();

  logInfo() << "Setting up DR...";

  const auto drLayout = internal::layoutDR(colors, colorsGhost, colorMap, meshReader);

  auto& drStorage = seissolInstance.memoryManager().drStorage();

  drStorage.setName("dr");

  /// Dynamic rupture storage
  seissolInstance.memoryManager().drDescriptor().addTo(drStorage);

  drStorage.setLayerCount(ltsStorage.getColorMap());
  drStorage.setSnapshots(ltsStorage.snapshots());
  drStorage.fixate();

  for (auto& layer : drStorage.leaves()) {
    layer.setNumberOfCells(drLayout[layer.id()].size() * copyCount);
  }

  drStorage.allocateVariables();
  drStorage.touchVariables();

  auto& drBackmap = seissolInstance.memoryManager().drBackmap();
  drBackmap.setSize(meshReader.getFault().size(), copyCount);

  const auto* zeroDR = drStorage.var<DynamicRupture::FaceInformation>();
  for (auto& layer : drStorage.leaves()) {
    auto* zeroDRLayer = layer.var<DynamicRupture::FaceInformation>();
    assert(layer.getIdentifier().halo != HaloType::Ghost || layer.size() == 0);
    for (std::size_t i = 0; i < layer.size(); ++i) {
      const auto cell = i / copyCount;
      const auto copy = i % copyCount;
      zeroDRLayer[i].meshFace = drLayout[layer.id()][cell];
      zeroDRLayer[i].copy = copy;
      drBackmap.addElement(layer.id(), zeroDR, zeroDRLayer, drLayout[layer.id()][cell], copy, i);
    }
  }

  if constexpr (isDeviceOn()) {
    seissol::initializer::internal::deriveRequiredScratchpadMemoryForDr(drStorage);
    drStorage.allocateScratchPads(useScratchpadPerLayer(seissolInstance.env()) ||
                                          useConcurrentClusters(seissolInstance.env())
                                      ? initializer::ScratchpadSharing::PerLayer
                                      : initializer::ScratchpadSharing::Shared);
  }

  // pass 4: correct LTS setup, again. Do bucket setup, determine communication datastructures
  logInfo() << "Setting up data exchange and face displacements (buckets)...";
  const auto haloCommunication = internal::bucketsAndCommunication(ltsStorage, meshLayoutCopy);

  logInfo() << "Setting up kernel clusters...";
  // the dynamic rupture, once its friction law and its outputs exist
  features.push_back(
      std::make_unique<seissol::dr::DynamicRuptureFeature>(seissolInstance,
                                                           memoryManager.drStorage(),
                                                           memoryManager.drDescriptor(),
                                                           meshReader,
                                                           settings.copyCount(),
                                                           memoryManager.globalDataStore(),
                                                           memoryManager.frictionLaw(),
                                                           memoryManager.frictionLawDevice(),
                                                           memoryManager.faultOutputManager(),
                                                           seissolParams.timeStepping.endTime));

  // the clusters of the faces of the features, e.g. of the dynamic rupture
  std::vector<seissol::solver::FaceClusterFactory*> faceClusters;
  for (auto& feature : features) {
    if (auto* factory = feature->faceClusters(); factory != nullptr) {
      faceClusters.push_back(factory);
    }
  }
  // the clusters of the cells, with the numerics of the wave propagation
  seissol::solver::WaveClusterFactory cellClusters(
      seissolInstance,
      memoryManager.globalDataStore(),
      seissolParams.runConfigs.configs(),
      memoryManager.initialConditions(),
      seissolInstance.gravitationSetup().acceleration,
      [&seissolInstance](std::size_t cells) { return seissolInstance.executionPlace(cells); });
  seissolInstance.timeManager().addClusters(
      clusterLayout, haloCommunication, ltsStorage, cellClusters, faceClusters);

  // the work of the features on the clusters of the cells
  for (auto& feature : features) {
    feature->addStages(seissolInstance.timeManager(), ltsStorage);
  }

  seissolInstance.dofSync().setup(meshLayoutCopy, &ltsStorage);
}

/**
 * Decides which clusters run on the host; before the clusters get set up, which fix their executor.
 */
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
      logWarning() << "Using the host-device execution without unified memory is not fully "
                      "supported yet; expect incorrect results. SEISSOL_USM=1 enables unified "
                      "memory.";
    }
  }
}

} // namespace

void initLayout(seissol::SeisSol& seissolInstance) {
  logInfo() << "Begin init layout.";

  hostDeviceCoexecution(seissolInstance);

  setupMemory(seissolInstance);

  logInfo() << "End init layout.";
}
} // namespace seissol::initializer::initprocedure
