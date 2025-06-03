// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "InitLayout.h"
#include "Internal/MeshLayout.h"
#include "SeisSol.h"
#include <Initializer/InitProcedure/Internal/Buckets.h>
#include <Initializer/InitProcedure/Internal/LtsSetup.h>
#include <Memory/Tree/LTSTree.h>

namespace {
using namespace seissol::initializer;

std::vector<initializer::LayerDefinition> determineWPLayerStructure() {}

std::vector<initializer::LayerDefinition> determineDRLayerStructure() {}

void setupMemory(seissol::SeisSol& seissolInstance) {
  const auto& meshReader = seissolInstance.meshReader();
  auto& container = seissolInstance.getMemoryManager().memoryContainer();

  std::vector<std::size_t> colors;
  std::vector<std::size_t> colorsGhost;

  const auto meshLayout =
      internal::layoutCells(colors, colorsGhost, container.colorMap.size(), meshReader);

  const auto wpStructure = determineWPLayerStructure();
  const auto faceStructure = determineDRLayerStructure();
  container.wpdesc.addTo(container.volume, seissolInstance.getSeisSolParameters().model.plasticity);
  container.volume.initialize(container.colorMap, wpStructure);

  logInfo() << "Setting up cell storage...";
  // just need pick a test variable that is available everywhere
  const auto& testVariable = container.wpdesc.cellInformation;
  const auto& testTree = container.volume.var(testVariable);
  for (auto& layer : container.volume.leaves()) {
    const auto& layout = meshLayout[container.colorMap.colorId(layer.getIdentifier())];
    const auto& testLayer = layer.var(testVariable);
    auto* secondaryCellInformation = layer.var(container.wpdesc.secondaryInformation);

    auto addToBackmapCopyInterior = [&](auto cell, auto index) {
      return container.clusterBackmap.addElement(layer.id(), testTree, testLayer, cell, index);
    };
    auto addToBackmapGhost = [&](auto cell, auto index) {
      return container.ghostClusterBackmap.addElement(layer.id(), testTree, testLayer, cell, index);
    };
    const auto& addToBackmap = [&](auto cell, auto index) {
      if (layer.getIdentifier().halo == HaloType::Interior ||
          layer.getIdentifier().halo == HaloType::Copy) {
        return addToBackmapCopyInterior(cell, index);
      } else {
        return addToBackmapGhost(cell, index);
      }
    };

    std::size_t index = 0;
    for (const auto& cell : layout.cells(layer.getIdentifier().halo)) {
      // TODO: two boundary elements with the same source (relevant for Ghost)
      auto dup = addToBackmap(cell, index);

      secondaryCellInformation[index].meshId = cell;
      secondaryCellInformation[index].configId = layer.getIdentifier().config.index();
      secondaryCellInformation[index].clusterId = layer.getIdentifier().lts;
      secondaryCellInformation[index].duplicate = dup;
      secondaryCellInformation[index].rank = -1;
      for (int face = 0; face < 4; ++face) {
        secondaryCellInformation[index].neighborRanks[face] = -1;
      }
      ++index;
    }
  }

  logInfo() << "Retrieving material data for the last time...";
  auto [materialsDB, plasticityDB] = seissol::model::queryMaterial(
      configs, seissol::initializer::CellToVertexArray::fromMeshReader(meshReader), true);
  auto [materialsDBGhost, plasticityDBGhost] = seissol::model::queryMaterial(
      configs,
      seissol::initializer::CellToVertexArray::fromVectors(ghostVertices, ghostGroups),
      true);
  logInfo() << "Setting up cell data...";
  for (auto& layer : container.volume.leaves()) {
    const auto& layout = meshLayout[container.colorMap.colorId(layer.getIdentifier())];
    auto* materialData = layer.var(container.wpdesc.materialData);
    auto* secondaryCellInformation = layer.var(container.wpdesc.secondaryInformation);
    auto* cellInformation = layer.var(container.wpdesc.cellInformation);
    const auto& cells = layout.cells(layer.getIdentifier().halo);
    if (layer.getIdentifier().halo == HaloType::Interior ||
        layer.getIdentifier().halo == HaloType::Copy) {
      auto* material = layer.var(container.wpdesc.material);
      auto* plasticity = layer.var(container.wpdesc.plasticity);
      std::size_t index = 0;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (std::size_t i = 0; i < cells.size(); ++i) {
        const auto cell = cells[i];

        // TODO: parametrize
        new (materialData[index]) model::MaterialT(materialsDB[cell]);
        new (plasticity[index]) model::PlasticityData(plasticityDB[cell], &materialData[index]);

        material[index].local = &materialData[index];
        secondaryCellInformation[index].meshId = cell;
        secondaryCellInformation[index].configId = layer.getIdentifier().config.index();
        secondaryCellInformation[index].clusterId = layer.getIdentifier().lts;
        secondaryCellInformation[index].rank = MPI::mpi.rank();
        for (int face = 0; face < 4; ++face) {
          bool ghostNeighbor =
              meshReader.getElements()[cell].neighbors[face] == meshReader.getElements().size();
          const auto& neighbor = [&]() {
            if (ghostNeighbor) {
              // ghost layer
              return container.ghostClusterBackmap.storagePositionLookup(
                  meshReader.getBoundaryElements()[cell].boundaryElements[face]);
            } else {
              // copy/interior layer
              return container.clusterBackmap.storagePositionLookup(
                  meshReader.getElements()[cell].neighbors[face]);
            }
          }();
          secondaryCellInformation[index].neighborRanks[face] =
              meshReader.getElements()[cell].neighborRanks[face];
          secondaryCellInformation[index].faceNeighborIds[face] = neighbor.second;
          cellInformation[index].neighborConfigIds[face] = neighbor.first;
          cellInformation[index].faceTypes[face] =
              static_cast<FaceType>(meshReader.getElements()[cell].boundaries[face]);
          cellInformation[index].faceRelations[face][0] =
              meshReader.getElements()[cell].neighborSides[face];
          cellInformation[index].faceRelations[face][1] =
              meshReader.getElements()[cell].sideOrientations[face];

          // TODO: parametrize
          const auto& materialData2 = ltsview2.layer(neighbor.first).var(lts.materialData);
          material[index].neighbor[face] = &materialData2[neighbor.second];
        }
        ++index;
      }
    } else if (layer.getIdentifier().halo == HaloType::Ghost) {
      std::size_t index = 0;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (std::size_t i = 0; i < cells.size(); ++i) {
        const auto cell = cells[i];

        // TODO: parametrize
        new (materialData[index]) seissol::model::MaterialT(materialsDBGhost[cell]);
        secondaryCellInformation[index].meshId = cell;
        secondaryCellInformation[index].configId = layer.getIdentifier().config.index();
        secondaryCellInformation[index].clusterId = layer.getIdentifier().lts;

        const auto& boundaryElement = meshReader.getBoundaryElements()[cell].get();
        secondaryCellInformation[index].rank = boundaryElement.neighborRank;
        // auto neighbor = container.clusterBackmap.layerPosition(boundaryElement.localElement);
        secondaryCellInformation[index].faceNeighborIds[boundaryElement.localSide] =
            boundaryElement.localElement;

        const int face = boundaryElement.neighborSide;

        secondaryCellInformation[index].neighborRanks[face] = MPI::mpi.rank();
        cellInformation[index].neighborConfigIds[face] = configs[boundaryElement.group].configId;
        cellInformation[index].faceTypes[face] =
            static_cast<FaceType>(meshReader.getElements()[boundaryElement.localElement]
                                      .boundaries[boundaryElement.localSide]);
        cellInformation[index].faceRelations[face][0] = boundaryElement.localSide;
        cellInformation[index].faceRelations[face][1] =
            meshReader.getElements()[boundaryElement.localElement]
                .sideOrientations[boundaryElement.localSide];
        ++index;
      }
    }
  }

  // pass 3: LTS setup
  logInfo() << "Setting up LTS configuration...";
  internal::handleLtsSetup(container);

  // pass 4: correct LTS setup, again. Do bucket setup, determine communication datastructures
  logInfo() << "Setting up data exchange and face displacements (buckets)...";
  internal::bucketsAndCommunication(container, meshLayout);

  // TODO init DR and boundary
  logInfo() << "Initializing Dynamic Rupture tree and data structures...";

  // allocate boundary tree (effectively behaves like a bucket, but is actually better like that
  // (checkpointing))
  logInfo() << "Initializing Boundary tree...";

  logInfo() << "Allocate scratch pads on GPUs...";
}

} // namespace

namespace seissol::initializer::initprocedure {
void initLayout(seissol::SeisSol& seissolInstance) {
  logInfo() << "Begin init layout.";

  setupMemory(seissolInstance);

  logInfo() << "End init layout.";
}
} // namespace seissol::initializer::initprocedure
