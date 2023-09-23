#include "MeshLayout.hpp"

#include <Geometry/MeshReader.h>
#include <Initializer/tree/Layer.hpp>
#include <Initializer/typedefs.hpp>
#include <mpi.h>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace seissol::initializer::internal {

std::vector<ClusterLayout> layoutCells(const std::vector<int>& color, const std::vector<int>& ghostColor, std::size_t maxColors, const geometry::MeshReader& meshReader) {
    // add cells to color clusters

    std::vector<ClusterLayout> clusters(maxColors);
    std::vector<std::unordered_map<int, std::vector<int>>> clusterCopy(maxColors);
    std::vector<std::unordered_map<int, std::vector<int>>> clusterGhost(maxColors);
    for (const auto& [rank, cells] : meshReader.getMPINeighbors()) {
        for (const auto& cell : cells) {
            // TODO(David): transfer color to neighbors
            clusterGhost[ghostColor[cell.linearId]][rank].push_back(cell.localElement);
        }
    }
    for (int i = 0; i < meshReader.getElements().size(); ++i) {
        const auto& element = meshReader.getElements()[i];
        bool interiorCell = true;
        for (int j = 0; j < 4; ++j) {
            if (element.neighborRanks[j] != MPI::mpi.rank()) {
                clusterCopy[color[i]][element.neighborRanks[j]].push_back(i);
                interiorCell = false;
            }
        }
        if (interiorCell) {
            clusters[color[i]].interior.push_back(i);
        }
    }

    // sort copy/ghost regions by MPI index
    for (auto& cluster : clusters) {
        std::sort(cluster.copy.begin(), cluster.copy.end(), [&](const auto& first, const auto& second) {
            return meshReader.getElements()[first].mpiIndices < meshReader.getElements()[second].mpiIndices;
        });
        std::sort(cluster.ghost.begin(), cluster.ghost.end(), [&](const auto& first, const auto& second) {
            return meshReader.getElements()[first].mpiIndices < meshReader.getElements()[second].mpiIndices;
        });
    }
    return clusters;
}

} // namespace seissol::initializer::internal
