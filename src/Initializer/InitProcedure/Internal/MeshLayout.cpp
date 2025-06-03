// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "MeshLayout.h"

#include <Geometry/MeshReader.h>
#include <Parallel/MPI.h>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace seissol::initializer::internal {

std::vector<ClusterLayout> layoutCells(const std::vector<std::size_t>& color,
                                       const std::vector<std::size_t>& ghostColor,
                                       std::size_t maxColors,
                                       const geometry::MeshReader& meshReader) {
  // add cells to color clusters

  std::vector<ClusterLayout> clusters(maxColors);
  std::vector<std::unordered_map<std::size_t, std::vector<std::size_t>>> clusterCopy(maxColors);
  std::vector<std::unordered_map<std::size_t, std::vector<std::size_t>>> clusterGhost(maxColors);

  std::size_t linearId = 0;
  for (const auto& [rank, cells] : meshReader.getMPINeighbors()) {
    for (const auto& cell : cells.elements) {
      clusterGhost[ghostColor[linearId]][rank].push_back(cell.localElement);
      ++linearId;
    }
  }
  for (std::size_t i = 0; i < meshReader.getElements().size(); ++i) {
    clusters[color[i]].cellMap.push_back(i);
  }

  // TODO: regions (take neighbor maximum)

  // sort copy/ghost regions by MPI index
  for (auto& cluster : clusters) {
    /*std::sort(cluster.copy.begin(), cluster.copy.end(), [&](const auto& first, const auto& second)
    { return meshReader.getElements()[first].mpiIndices <
             meshReader.getElements()[second].mpiIndices;
    });
    std::sort(
        cluster.ghost.begin(), cluster.ghost.end(), [&](const auto& first, const auto& second) {
          return meshReader.getElements()[first].mpiIndices <
                 meshReader.getElements()[second].mpiIndices;
        });*/
  }
  return clusters;
}

} // namespace seissol::initializer::internal
