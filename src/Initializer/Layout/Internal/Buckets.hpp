#pragma once

#include <Initializer/Layout/Internal/MeshLayout.hpp>
#include <Initializer/Layout/Memory.hpp>
#include <Initializer/tree/LTSForest.hpp>
#include <Initializer/tree/Layer.hpp>
#include <vector>

namespace seissol::initializer::internal {
    struct CommunicationInfo {
        void* buffer;
        MPI_Datatype datatype;
        int count;
        int tag;
        int rank;
        MPI_Request request;
    };

    template<typename Config>
    std::vector<CommunicationInfo> allocateTransferInfo(seissol::initializers::LTSLayerRef<seissol::initializers::LTS<Config>>& layerView);

    template<typename Config>
    void allocateFaceDisplacements(seissol::initializers::LTSLayerRef<seissol::initializers::LTS<Config>>& layerView);

    template<typename Config>
    void touchBuckets(seissol::initializers::LTSLayerRef<seissol::initializers::LTS<Config>>& layerView);

    void bucketsAndCommunication(MemoryContainer& container, const std::vector<ClusterLayout>& meshLayout);
} // namespace seissol::initializer::internal
