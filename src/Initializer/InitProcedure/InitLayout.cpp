// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "InitLayout.h"

#include "Common/Constants.h"
#include "Common/Iterator.h"
#include "Geometry/MeshReader.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/InitProcedure/Internal/Buckets.h"
#include "Initializer/InitProcedure/Internal/LtsSetup.h"
#include "Initializer/TimeStepping/ClusterLayout.h"
#include "Initializer/TimeStepping/Halo.h"
#include "Internal/MeshLayout.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Backmap.h"
#include "Memory/Tree/Colormap.h"
#include "Memory/Tree/LTSTree.h"
#include "Memory/Tree/Layer.h"
#include "Parallel/MPI.h"
#include "SeisSol.h"

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <mpi.h>
#include <numeric>
#include <unordered_map>
#include <utils/logger.h>
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

void setupMemory(seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.getSeisSolParameters();
  const auto& meshReader = seissolInstance.meshReader();

  const auto rank = seissol::Mpi::mpi.rank();

  logInfo() << "Determining cell colors...";

  const auto clusterLayout = ClusterLayout::fromMesh(seissolParams.timeStepping.lts.getRate(),
                                                     seissolInstance.meshReader(),
                                                     seissolInstance.getTimestepScale(),
                                                     true);

  seissolInstance.getMemoryManager().setClusterLayout(clusterLayout);

  std::vector<std::size_t> clusterMap(clusterLayout.globalClusterCount);
  std::iota(clusterMap.begin(), clusterMap.end(), 0);

  const LTSColorMap colorMap(
      initializer::EnumLayer<HaloType>({HaloType::Ghost, HaloType::Copy, HaloType::Interior}),
      initializer::EnumLayer<std::size_t>(clusterMap),
      initializer::TraitLayer<initializer::ConfigVariant>({initializer::ConfigVariant(Config())}));

  std::vector<std::size_t> colors(meshReader.getElements().size());
  for (std::size_t i = 0; i < colors.size(); ++i) {
    const auto& element = meshReader.getElements()[i];
    const auto halo = geometry::isCopy(element, rank) ? HaloType::Copy : HaloType::Interior;
    colors[i] = colorMap.color(halo, element.clusterId, Config());
  }

  const auto ghostSize = meshReader.linearGhostlayer().size();

  std::vector<std::size_t> colorsGhost(ghostSize);
  for (const auto [i, linearGhost] : common::enumerate(meshReader.linearGhostlayer())) {
    const auto& element =
        meshReader.getGhostlayerMetadata().at(linearGhost.rank)[linearGhost.inRankIndices[0]];
    const auto halo = HaloType::Ghost;
    colorsGhost[i] = colorMap.color(halo, element.clusterId, Config());
  }

  logInfo() << "Creating mesh layout...";

  const auto meshLayout = internal::layoutCells(colors, colorsGhost, colorMap, meshReader);

  auto& ltsStorage = seissolInstance.getMemoryManager().getLtsStorage();
  auto& backmap = seissolInstance.getMemoryManager().getBackmap();
  LTS::addTo(ltsStorage, seissolInstance.getSeisSolParameters().model.plasticity);
  seissolInstance.postProcessor().allocateMemory(ltsStorage);
  ltsStorage.setName("cluster");
  ltsStorage.setLayerCount(colorMap);
  ltsStorage.fixate();
  for (auto& layer : ltsStorage.leaves()) {
    layer.setNumberOfCells(meshLayout[layer.id()].cellMap.size());
  }
  ltsStorage.allocateVariables();
  ltsStorage.touchVariables();
  backmap.setSize(meshReader.getElements().size());

  StorageBackmap<Cell::NumFaces> backmapGhost;
  backmapGhost.setSize(ghostSize);
  std::vector<std::unordered_map<std::size_t, StoragePosition>> backmapGhostMap(ghostSize);

  logInfo() << "Setting up cell storage...";
  // just need pick a test variable that is available everywhere
  const auto* zero = ltsStorage.var<LTS::SecondaryInformation>();
  for (auto& layer : ltsStorage.leaves()) {
    auto* zeroLayer = layer.var<LTS::SecondaryInformation>();
    const auto& layout = meshLayout[layer.id()];
    auto regionIt = layout.regions.begin();
    std::size_t regionOffset = 0;

    const auto addToBackmapCopyInterior = [&](auto cell, auto index) {
      return backmap.addElement(layer.id(), zero, zeroLayer, cell, index);
    };
    const auto addToBackmapGhost = [&](auto cell, auto index) {
      assert(regionIt != layout.regions.end());
      const auto dup = backmapGhost.addElement(layer.id(), zero, zeroLayer, cell, index);
      backmapGhostMap.at(cell)[regionIt->remoteId] = backmapGhost.getDup(cell, dup).value();
      return dup;
    };
    const auto& addToBackmap = [&](auto cell, auto index) {
      if (layer.getIdentifier().halo == HaloType::Ghost) {
        return addToBackmapGhost(cell, index);
      } else {
        return addToBackmapCopyInterior(cell, index);
      }
    };

    for (const auto [i, cell] : common::enumerate(layout.cellMap)) {
      while (regionIt != layout.regions.end() && i >= regionOffset + regionIt->count) {
        regionOffset += regionIt->count;
        ++regionIt;
      }

      const auto dup = addToBackmap(cell, i);

      zeroLayer[i].meshId = cell;
      zeroLayer[i].configId = layer.getIdentifier().config.index();
      zeroLayer[i].clusterId = layer.getIdentifier().lts;
      zeroLayer[i].duplicate = dup;
      zeroLayer[i].halo = layer.getIdentifier().halo;
      zeroLayer[i].rank = -1;
      zeroLayer[i].color = layer.id();
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        zeroLayer[i].neighborRanks[face] = -1;
        zeroLayer[i].faceNeighbors[face] = StoragePosition::NullPosition;
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
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (std::size_t i = 0; i < cells.size(); ++i) {
        const auto cell = cells[i];
        const auto index = i;
        const auto& element = meshReader.getElements()[cell];

        secondaryCellInformation[index].rank = rank;
        secondaryCellInformation[index].globalId = element.globalId;
        secondaryCellInformation[index].group = element.group;
        for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
          secondaryCellInformation[index].neighborRanks[face] = element.neighborRanks[face];
          cellInformation[index].faceTypes[face] = static_cast<FaceType>(element.boundaries[face]);
          cellInformation[index].faceRelations[face][0] = element.neighborSides[face];
          cellInformation[index].faceRelations[face][1] = element.sideOrientations[face];

          if (element.neighbors[face].hasValue() || element.neighborRanks[face] != rank) {
            const auto& neighbor = [&]() {
              const bool ghostNeighbor = element.neighborRanks[face] != rank;
              if (ghostNeighbor) {
                const auto rank = element.neighborRanks[face];
                const auto mpiIndex = element.mpiIndices[face];
                const auto linear = meshReader.toLinearGhostlayer().at({rank, mpiIndex});
                // ghost layer
                return backmapGhostMap.at(linear).at(layer.id());
              } else {
                // copy/interior layer
                assert(element.neighbors[face].hasValue());
                return backmap.get(element.neighbors[face].value());
              }
            }();

            secondaryCellInformation[index].faceNeighbors[face] = neighbor;
            cellInformation[index].neighborConfigIds[face] = 0;
          } else {
            secondaryCellInformation[index].faceNeighbors[face] = StoragePosition::NullPosition;
            cellInformation[index].neighborConfigIds[face] =
                std::numeric_limits<std::uint32_t>::max();
          }
        }
      }
    } else if (layer.getIdentifier().halo == HaloType::Ghost) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (std::size_t i = 0; i < cells.size(); ++i) {
        const auto cell = cells[i];

        const auto& linear = meshReader.linearGhostlayer()[cell];
        secondaryCellInformation[i].rank = linear.rank;

        for (const auto& index : linear.inRankIndices) {
          const auto& boundaryElement =
              meshReader.getMPINeighbors().at(linear.rank).elements[index];
          const auto neighbor = backmap.get(boundaryElement.localElement);
          const auto& elementNeighbor = meshReader.getElements()[boundaryElement.localElement];

          secondaryCellInformation[i].globalId =
              meshReader.getGhostlayerMetadata().at(linear.rank)[index].globalId;

          // equivalent to elementNeighbor.neighborSides[boundaryElement.localSide]
          const auto face = boundaryElement.neighborSide;

          secondaryCellInformation[i].faceNeighbors[face] = neighbor;

          secondaryCellInformation[i].neighborRanks[face] = rank;
          cellInformation[i].neighborConfigIds[face] = neighbor.color;
          cellInformation[i].faceTypes[face] =
              static_cast<FaceType>(elementNeighbor.boundaries[boundaryElement.localSide]);
          cellInformation[i].faceRelations[face][0] = boundaryElement.localSide;
          cellInformation[i].faceRelations[face][1] =
              elementNeighbor.sideOrientations[boundaryElement.localSide];
        }
      }
    }
  }

  logInfo() << "Verify the halo setup...";
  verifyHaloSetup(ltsStorage, meshLayout);

  // pass 3: LTS setup
  logInfo() << "Setting up LTS configuration...";
  internal::deriveLtsSetups(meshLayout, ltsStorage);

  if (seissolParams.model.plasticity) {
    // remove disabled plasticity groups from the list
    const auto& pdis = seissolParams.model.plasticityDisabledGroups;
    for (auto& layer : ltsStorage.leaves(Ghost)) {
      const std::size_t size = layer.size();
      auto* cellInfo = layer.var<LTS::CellInformation>();
      const auto* cellInfo2 = layer.var<LTS::SecondaryInformation>();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (std::size_t i = 0; i < size; ++i) {
        cellInfo[i].plasticityEnabled = pdis.find(cellInfo2[i].group) == pdis.end();
      }
    }
  }

  seissolInstance.getMemoryManager().initializeFrictionLaw();

  logInfo() << "Setting up DR...";

  const auto drLayout = internal::layoutDR(colors, colorsGhost, colorMap, meshReader);

  auto& drStorage = seissolInstance.getMemoryManager().getDRStorage();

  drStorage.setName("dr");

  /// Dynamic rupture storage
  seissolInstance.getMemoryManager().getDynamicRupture().addTo(drStorage);

  drStorage.setLayerCount(ltsStorage.getColorMap());
  drStorage.fixate();

  for (auto& layer : drStorage.leaves()) {
    layer.setNumberOfCells(drLayout[layer.id()].size());
  }

  drStorage.allocateVariables();
  drStorage.touchVariables();

  auto& drBackmap = seissolInstance.getMemoryManager().getDRBackmap();
  drBackmap.setSize(meshReader.getFault().size());

  const auto* zeroDR = drStorage.var<DynamicRupture::FaceInformation>();
  for (auto& layer : drStorage.leaves()) {
    auto* zeroDRLayer = layer.var<DynamicRupture::FaceInformation>();
    assert(layer.getIdentifier().halo != HaloType::Ghost || layer.size() == 0);
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      zeroDRLayer[cell].meshFace = drLayout[layer.id()][cell];
      drBackmap.addElement(layer.id(), zeroDR, zeroDRLayer, drLayout[layer.id()][cell], cell);
    }
  }

  seissolInstance.getMemoryManager().fixateLtsStorage();

  // pass 4: correct LTS setup, again. Do bucket setup, determine communication datastructures
  logInfo() << "Setting up data exchange and face displacements (buckets)...";
  const auto haloCommunication = internal::bucketsAndCommunication(ltsStorage, meshLayout);

  logInfo() << "Setting up kernel clusters...";
  seissolInstance.timeManager().addClusters(clusterLayout,
                                            haloCommunication,
                                            seissolInstance.getMemoryManager(),
                                            seissolParams.model.plasticity);
}

} // namespace

void initLayout(seissol::SeisSol& seissolInstance) {
  logInfo() << "Begin init layout.";

  setupMemory(seissolInstance);

  logInfo() << "End init layout.";
}
} // namespace seissol::initializer::initprocedure
