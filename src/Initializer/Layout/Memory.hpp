#pragma once

#include <Common/configs.hpp>
#include <Geometry/MeshReader.h>
#include <Initializer/GlobalData.h>
#include <Initializer/tree/Backmap.hpp>
#include <Initializer/tree/LTSForest.hpp>
#include <Initializer/tree/LayerMap.hpp>
#include <Initializer/typedefs.hpp>
#include <Solver/time_stepping/TimeManager.h>

namespace seissol::initializer {
    struct MemoryContainer {
        initializers::ClusterColorMap colorMap;

        initializers::ClusterLTSForest cluster;
        initializers::DynRupLTSForest dynrup;
        initializers::BoundaryLTSForest boundary;

        ClusterBackmap clusterBackmap;
        initializers::StorageBackmap<1, initializers::ColorMap<initializers::EnumLayer<SupportedConfigs>>> ghostClusterBackmap;
        initializers::StorageBackmap<1, initializers::ColorMap<initializers::EnumLayer<SupportedConfigs>>> dynrupBackmap;

        initializers::GlobalDataStorage globalDataStorage;
    };

    MemoryContainer setupMemory(const geometry::MeshReader& meshReader, seissol::time_stepping::TimeManager& timeManager);
} // namespace seissol::initializer
