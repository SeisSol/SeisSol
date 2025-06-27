// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "InitLayout.h"
#include "Internal/MeshLayout.h"
#include "SeisSol.h"
#include <Equations/Datastructures.h>
#include <Geometry/MeshReader.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/InitProcedure/Internal/Buckets.h>
#include <Initializer/InitProcedure/Internal/LtsSetup.h>
#include <Initializer/ParameterDB.h>
#include <Initializer/TimeStepping/ClusterLayout.h>
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
  const auto& meshReader = seissolInstance.meshReader();

  std::size_t maxLtsId = 0;
  double minimumTimestep = std::numeric_limits<double>::max();
  for (const auto& element : meshReader.getElements()) {
    maxLtsId = std::max(maxLtsId, static_cast<std::size_t>(element.clusterId));
    minimumTimestep = std::min(minimumTimestep, element.timestep);
  }
  MPI_Allreduce(
      MPI_IN_PLACE, &maxLtsId, 1, MPI::castToMpiType<std::size_t>(), MPI_MAX, MPI::mpi.comm());
  MPI_Allreduce(
      MPI_IN_PLACE, &minimumTimestep, 1, MPI_DOUBLE, MPI_MIN, MPI::mpi.comm());
  ClusterLayout layout({seissolInstance.getSeisSolParameters().timeStepping.lts.getRate()}, minimumTimestep, maxLtsId + 1);

  seissolInstance.getMemoryManager().setupMemoryContainer(maxLtsId, {Config()});
  auto& container = seissolInstance.getMemoryManager().memoryContainer();

  std::vector<std::size_t> colors(meshReader.getElements().size());
  for (std::size_t i = 0; i < colors.size(); ++i) {
    const auto& element = meshReader.getElements()[i];
    const auto halo =
        geometry::isCopy(element, MPI::mpi.rank()) ? HaloType::Copy : HaloType::Interior;
    colors[i] = container.colorMap.color(Config(), element.clusterId, halo);
  }

  std::vector<std::size_t> colorsGhost(meshReader.getBoundaryElements().size());
  for (std::size_t i = 0; i < colorsGhost.size(); ++i) {
    const auto& element = meshReader.getBoundaryElements()[i];
    const auto halo = HaloType::Ghost;
    colorsGhost[i] = container.colorMap.color(Config(), element.clusterId, halo);
  }

  const auto meshLayout =
      internal::layoutCells(colors, colorsGhost, container.colorMap.size(), meshReader);

  std::vector<std::size_t> wpStructure(container.colorMap.size());
  for (std::size_t i = 0; i < wpStructure.size(); ++i) {
    wpStructure[i] = meshLayout[i].cellMap.size();
  }

  container.wpdesc.addTo(container.volume, seissolInstance.getSeisSolParameters().model.plasticity);
  container.volume.initialize(container.colorMap, wpStructure);
  container.volume.allocateTouchVariables();
  container.clusterBackmap.setSize(meshReader.getElements().size());

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
    for (const auto& cell : layout.cellMap) {
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

  const auto& seissolParams = seissolInstance.getSeisSolParameters();
  const auto getBestQueryGenerator = [&](const seissol::initializer::CellToVertexArray& ctvArray) {
    return seissol::initializer::getBestQueryGenerator(
        seissolParams.model.plasticity, seissolParams.model.useCellHomogenizedMaterial, ctvArray);
  };

  const auto queryMaterial = [&](const auto& data, bool plasicity)
      -> std::pair<std::vector<model::MaterialT>, std::vector<model::Plasticity>> {
    const auto queryGen = getBestQueryGenerator(data);
    std::vector<model::MaterialT> vectorDB(data.size);
    seissol::initializer::MaterialParameterDB<model::MaterialT> parameterDB;
    parameterDB.setMaterialVector(&vectorDB);
    parameterDB.evaluateModel(seissolParams.model.materialFileName, *queryGen);

    std::vector<model::Plasticity> vectorDBP(plasicity ? data.size : 0);
    if (plasicity) {
      seissol::initializer::MaterialParameterDB<model::Plasticity> parameterDBP;
      parameterDBP.setMaterialVector(&vectorDBP);
      parameterDBP.evaluateModel(seissolParams.model.materialFileName, *queryGen);
    }
    return {vectorDB, vectorDBP};
  };

  auto innerDB = queryMaterial(seissol::initializer::CellToVertexArray::fromMeshReader(meshReader),
                               seissolParams.model.plasticity);
  auto ghostDB = queryMaterial(
      seissol::initializer::CellToVertexArray::fromVectors(ghostVertices, ghostGroups), false);
  logInfo() << "Setting up cell data...";
  for (auto& layer : container.volume.leaves()) {
    const auto& layout = meshLayout[container.colorMap.colorId(layer.getIdentifier())];
    auto* materialData = layer.var(container.wpdesc.materialData);
    auto* secondaryCellInformation = layer.var(container.wpdesc.secondaryInformation);
    auto* cellInformation = layer.var(container.wpdesc.cellInformation);
    const auto& cells = layout.cellMap;
    if (layer.getIdentifier().halo == HaloType::Interior ||
        layer.getIdentifier().halo == HaloType::Copy) {
      auto* material = layer.var(container.wpdesc.material);
      auto* plasticity = layer.var(container.wpdesc.plasticity);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (std::size_t i = 0; i < cells.size(); ++i) {
        const auto cell = cells[i];
        const auto index = i;

        // TODO: parametrize
        new (&materialData[index]) model::MaterialT(innerDB.first[cell]);
        if (seissolParams.model.plasticity) {
          new (&plasticity[index])
              model::PlasticityData(innerDB.second[cell], &materialData[index]);
        }

        material[index].local = &materialData[index];
        secondaryCellInformation[index].rank = MPI::mpi.rank();
        for (int face = 0; face < 4; ++face) {
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
                return container.ghostClusterBackmap.storagePositionLookup(
                    meshReader.getGhostlayerMetadata().at(rank)[mpiIndex].linearId);
              } else {
                // copy/interior layer
                return container.clusterBackmap.storagePositionLookup(element.neighbors[face]);
              }
            }();

            secondaryCellInformation[index].faceNeighborIds[face] = neighbor.global;
            cellInformation[index].neighborConfigIds[face] = neighbor.color;
            // TODO: parametrize
            auto& neighborLayer = container.volume.layer(neighbor.color);
            const auto& materialData2 = neighborLayer.var(container.wpdesc.materialData);
            material[index].neighbor[face] = &materialData2[neighbor.cell];
          } else {
            secondaryCellInformation[index].faceNeighborIds[face] =
                std::numeric_limits<std::size_t>::max();
            cellInformation[index].neighborConfigIds[face] = -1;
            material[index].neighbor[face] = nullptr;
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

        // TODO: parametrize
        new (&materialData[index]) seissol::model::MaterialT(ghostDB.first[cell]);

        const auto& boundaryElement = meshReader.getBoundaryElements()[cell];
        secondaryCellInformation[index].rank = boundaryElement.rank;
        const auto neighbor =
            container.clusterBackmap.storagePositionLookup(boundaryElement.localElement);
        secondaryCellInformation[index].faceNeighborIds[boundaryElement.localSide] =
            neighbor.global;

        const int face = boundaryElement.neighborSide;

        secondaryCellInformation[index].neighborRanks[face] = MPI::mpi.rank();
        cellInformation[index].neighborConfigIds[face] = neighbor.color;
        cellInformation[index].faceTypes[face] =
            static_cast<FaceType>(meshReader.getElements()[boundaryElement.localElement]
                                      .boundaries[boundaryElement.localSide]);
        cellInformation[index].faceRelations[face][0] = boundaryElement.localSide;
        cellInformation[index].faceRelations[face][1] =
            meshReader.getElements()[boundaryElement.localElement]
                .sideOrientations[boundaryElement.localSide];
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
