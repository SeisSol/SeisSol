#pragma once

#include <Geometry/MeshReader.h>
#include <Initializer/tree/Layer.hpp>
#include <Initializer/typedefs.hpp>
#include <mpi.h>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace seissol::initializer::internal {

struct TransferRegion {
    int rank;
    std::size_t start;
    std::size_t size;
};

struct ClusterLayout {
    std::vector<int> interior;
    std::vector<int> copy;
    std::vector<int> ghost;
    std::vector<TransferRegion> copyRegions;
    std::vector<TransferRegion> ghostRegions;

    const std::vector<int>& cells(LayerType icg) const {
        switch (icg) {
            case Interior:
                return interior;
            case Copy:
                return copy;
            case Ghost:
                return ghost;
            default:
                throw std::runtime_error("Invalid LayerType at setup");
        }
    }

    bool empty() const {
        return interior.empty() && copy.empty() && ghost.empty();
    }

    bool interiorOnly() const {
        return !interior.empty() && copy.empty() && ghost.empty();
    }
};

std::vector<ClusterLayout> layoutCells(const std::vector<int>& color, const std::vector<int>& ghostColor, std::size_t maxColors, const geometry::MeshReader& meshReader);

} // namespace seissol::initializer::internal
