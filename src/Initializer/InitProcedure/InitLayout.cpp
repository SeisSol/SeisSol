// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "InitLayout.h"
#include "Internal/MeshLayout.h"
#include "SeisSol.h"
#include <Common/Iterator.h>
#include <Equations/Datastructures.h>
#include <Geometry/MeshReader.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/InitProcedure/Internal/Buckets.h>
#include <Initializer/InitProcedure/Internal/LtsSetup.h>
#include <Initializer/ParameterDB.h>
#include <Initializer/TimeStepping/ClusterLayout.h>
#include <Initializer/TimeStepping/Halo.h>
#include <Memory/Tree/Backmap.h>
#include <Memory/Tree/Colormap.h>
#include <Memory/Tree/LTSTree.h>
#include <Model/CommonDatastructures.h>
#include <Model/Plasticity.h>
#include <Parallel/MPI.h>
#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>
#include <mpi.h>
#include <unordered_map>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace {
using namespace seissol::initializer;

void setupMemory(seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.getSeisSolParameters();
  const auto& meshReader = seissolInstance.meshReader();

  seissolInstance.getLtsLayout().setMesh(meshReader);
  seissolInstance.getLtsLayout().deriveLayout(seissolParams.timeStepping.lts.getRate());
  const auto clusterLayout = seissolInstance.getLtsLayout().clusterLayout();

  auto& ltsStorage = seissolInstance.getMemoryManager().getLtsStorage();
  auto& backmap = seissolInstance.getMemoryManager().getBackmap();

  std::vector<std::size_t> clusterMap(clusterLayout.globalClusterCount);
  std::iota(clusterMap.begin(), clusterMap.end(), 0);

  LTSColorMap colorMap(
      initializer::EnumLayer<HaloType>({HaloType::Ghost, HaloType::Copy, HaloType::Interior}),
      initializer::EnumLayer<std::size_t>(clusterMap),
      initializer::TraitLayer<initializer::ConfigVariant>({Config()}));

  std::vector<std::size_t> colors(meshReader.getElements().size());
  for (std::size_t i = 0; i < colors.size(); ++i) {
    const auto& element = meshReader.getElements()[i];
    const auto halo =
        geometry::isCopy(element, MPI::mpi.rank()) ? HaloType::Copy : HaloType::Interior;
    colors[i] = colorMap.color(halo, element.clusterId, Config());
  }

  std::size_t ghostSize = 0;
  for (const auto& [_, neighbor] : meshReader.getGhostlayerMetadata()) {
    ghostSize += neighbor.size();
  }

  std::vector<std::size_t> colorsGhost(ghostSize);
  std::size_t linearId = 0;
  std::map<std::pair<int, std::size_t>, std::size_t> toLinear;
  std::vector<std::pair<int, std::size_t>> fromLinear(ghostSize);
  for (const auto& [rank, neighbor] : meshReader.getGhostlayerMetadata()) {
    for (const auto [i, element] : common::enumerate(neighbor)) {
      const auto halo = HaloType::Ghost;
      colorsGhost[linearId] = colorMap.color(halo, element.clusterId, Config());
      toLinear[std::pair<int, std::size_t>(rank, i)] = linearId;
      fromLinear[linearId] = std::pair<int, std::size_t>(rank, i);
      ++linearId;
    }
  }

  const auto meshLayout = internal::layoutCells(colors, colorsGhost, colorMap, meshReader);

  std::vector<std::size_t> wpStructure(colorMap.size());
  for (std::size_t i = 0; i < wpStructure.size(); ++i) {
    wpStructure[i] = meshLayout[i].cellMap.size();
  }

  LTS::addTo(ltsStorage, seissolInstance.getSeisSolParameters().model.plasticity);
  seissolInstance.postProcessor().allocateMemory(ltsStorage);
  ltsStorage.setName("cluster");
  ltsStorage.setLayerCount(colorMap);
  ltsStorage.fixate();
  for (auto [i, layer] : common::enumerate(ltsStorage.leaves())) {
    layer.setNumberOfCells(meshLayout[i].cellMap.size());
  }
  ltsStorage.allocateVariables();
  ltsStorage.touchVariables();
  backmap.setSize(meshReader.getElements().size());

  StorageBackmap<1> backmapGhost;
  backmapGhost.setSize(ghostSize);

  logInfo() << "Setting up cell storage...";
  // just need pick a test variable that is available everywhere

  backmap.setSize(seissolInstance.meshReader().getElements().size());
  const auto* zero = ltsStorage.var<LTS::SecondaryInformation>();
  for (auto& layer : ltsStorage.leaves()) {
    auto* zeroLayer = layer.var<LTS::SecondaryInformation>();

    const auto addToBackmapCopyInterior = [&](auto cell, auto index) {
      return backmap.addElement(layer.id(), zero, zeroLayer, cell, index);
    };
    const auto addToBackmapGhost = [&](auto cell, auto index) {
      return backmapGhost.addElement(layer.id(), zero, zeroLayer, cell, index);
    };
    const auto& addToBackmap = [&](auto cell, auto index) {
      if (layer.getIdentifier().halo == HaloType::Ghost) {
        return addToBackmapGhost(cell, index);
      } else {
        return addToBackmapCopyInterior(cell, index);
      }
    };

    const auto& layout = meshLayout[layer.id()];
    for (const auto [i, cell] : common::enumerate(layout.cellMap)) {
      // TODO: two boundary elements with the same source (relevant for Ghost)
      const auto dup = addToBackmap(cell, i);

      zeroLayer[i].meshId = cell;
      zeroLayer[i].configId = layer.getIdentifier().config.index();
      zeroLayer[i].clusterId = layer.getIdentifier().lts;
      zeroLayer[i].duplicate = dup;
      zeroLayer[i].rank = -1;
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        zeroLayer[i].neighborRanks[face] = -1;
        zeroLayer[i].faceNeighbors[face] = StoragePosition::NullPosition;
      }
    }
  }

  logInfo() << "Setting up cell information data...";
  for (auto& layer : ltsStorage.leaves()) {
    const auto& layout = meshLayout[colorMap.colorId(layer.getIdentifier())];
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

        secondaryCellInformation[index].rank = MPI::mpi.rank();
        for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
          const auto& element = meshReader.getElements()[cell];
          secondaryCellInformation[index].neighborRanks[face] = element.neighborRanks[face];
          cellInformation[index].faceTypes[face] = static_cast<FaceType>(element.boundaries[face]);
          cellInformation[index].faceRelations[face][0] = element.neighborSides[face];
          cellInformation[index].faceRelations[face][1] = element.sideOrientations[face];

          if (element.neighbors[face] != meshReader.getElements().size() ||
              element.neighborRanks[face] != MPI::mpi.rank()) {
            const auto& neighbor = [&]() {
              const bool ghostNeighbor = element.neighborRanks[face] != MPI::mpi.rank();
              if (ghostNeighbor) {
                const auto rank = element.neighborRanks[face];
                const auto mpiIndex = element.mpiIndices[face];
                // ghost layer
                return backmapGhost.get(toLinear.at({rank, mpiIndex}));
              } else {
                // copy/interior layer
                return backmap.get(element.neighbors[face]);
              }
            }();

            secondaryCellInformation[index].faceNeighbors[face] = neighbor;
            cellInformation[index].neighborConfigIds[face] = 0;
          } else {
            secondaryCellInformation[index].faceNeighbors[face] = StoragePosition::NullPosition;
            cellInformation[index].neighborConfigIds[face] = -1;
          }
        }
      }
    } else if (layer.getIdentifier().halo == HaloType::Ghost) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (std::size_t i = 0; i < cells.size(); ++i) {
        const auto cell = cells[i];
        const auto index = i;

        const auto delinear = fromLinear[cell];
        const auto& boundaryElement =
            meshReader.getMPINeighbors().at(delinear.first).elements[delinear.second];
        secondaryCellInformation[index].rank = delinear.first;
        const auto neighbor = backmap.get(boundaryElement.localElement);
        const auto& elementNeighbor = meshReader.getElements()[boundaryElement.localElement];

        secondaryCellInformation[index].faceNeighbors[boundaryElement.localSide] = neighbor;

        const auto face = boundaryElement.neighborSide;

        secondaryCellInformation[index].neighborRanks[face] = MPI::mpi.rank();
        cellInformation[index].neighborConfigIds[face] = neighbor.color;
        cellInformation[index].faceTypes[face] =
            static_cast<FaceType>(elementNeighbor.boundaries[boundaryElement.localSide]);
        cellInformation[index].faceRelations[face][0] = boundaryElement.localSide;
        cellInformation[index].faceRelations[face][1] =
            elementNeighbor.sideOrientations[boundaryElement.localSide];
      }
    }
  }

  // pass 3: LTS setup
  logInfo() << "Setting up LTS configuration...";
  internal::deriveLtsSetups(meshLayout, ltsStorage);

  // intermediate
  const auto faceSizes = seissolInstance.getLtsLayout().drSizes();

  seissolInstance.getMemoryManager().initializeFrictionLaw();

  seissolInstance.getMemoryManager().fixateLtsStorage(faceSizes, seissolParams.model.plasticity);

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

namespace seissol::initializer::initprocedure {
void initLayout(seissol::SeisSol& seissolInstance) {
  logInfo() << "Begin init layout.";

  setupMemory(seissolInstance);

  logInfo() << "End init layout.";
}
} // namespace seissol::initializer::initprocedure
