#pragma once

#include <Initializer/Layout/Internal/MeshLayout.hpp>
#include <Initializer/Layout/Memory.hpp>
namespace seissol::initializer::internal {
    void setupCellInfo(MemoryContainer& container, const std::vector<ClusterLayout>& meshLayout);
} // namespace seissol::initializer::internal
