#include "Buckets.hpp"
#include <Common/precision.hpp>
#include <Initializer/Layout/Internal/MeshLayout.hpp>
#include <Initializer/tree/Layer.hpp>
#include <type_traits>
#include <yateto/InitTools.h>
#include <cstring>

namespace {
    class BucketManager {
    private:
        std::size_t dataSize = 0;
    public:
        template<typename T>
        T* markAllocate(bool allocate) {
            if (allocate) {
                std::size_t offset = dataSize;
                dataSize += sizeof(T);
                return reinterpret_cast<T*>(offset + 1);
            }
            else {
                return nullptr;
            }
        }

        std::size_t position() const {
            return dataSize + 1;
        }

        std::size_t size() const {
            return dataSize;
        }
    };

    template<typename T>
    static void initBucketItem(T*& data, void* bucket) {
        if (data != nullptr) {
            std::size_t offset = reinterpret_cast<std::size_t>(data) - 1;
            std::size_t bucketPos = reinterpret_cast<std::size_t>(bucket);
            data = reinterpret_cast<T*>(bucketPos + offset);
            std::memset(data, 0, sizeof(T));
        }
    }
} // namespace

namespace seissol::initializer::internal {
    template<typename Config>
    void allocateTransferInfo(seissol::initializers::LTSLayerRef<seissol::initializers::LTS<Config>>& layerView, const std::vector<TransferRegion>& regions) {
        using RealT = typename Config::RealT;

        auto* buffers = layerView.var(layerView.lts.buffers);
        auto* derivatives = layerView.var(layerView.lts.derivatives);
        const auto* cellInformation = layerView.var(layerView.lts.cellInformation);
        const auto* secondaryCellInformation = layerView.var(layerView.lts.secondaryCellInformation);
        BucketManager manager;

        auto allocate = [&](auto index, bool derivatives) {
            bool hasBuffers = (cellInformation[index].ltsSetup >> 8) != 0;
            bool hasDerivatives = (cellInformation[index].ltsSetup >> 9) != 0;
            if (derivatives) {
                derivatives[index] = manager.markAllocate<typename seissol::initializers::LTS<Config>::DerivType>(hasDerivatives);
            }
            else {
                buffers[index] = manager.markAllocate<typename seissol::initializers::LTS<Config>::DofsType>(hasBuffers);
            }
        };

        auto faceDerivatives = [&](auto index, auto face) {
            return (cellInformation[index].ltsSetup >> (4 + face)) != 0;
        };

        auto useBuffersDerivatives = [&](auto index, auto rank) {
            bool buffers = false;
                    bool derivatives = false;
                    for (int j = 0; j < 4; ++j) {
                        if ((secondaryCellInformation[index].rank == rank && secondaryCellInformation[index].neighborRanks[j] >= 0)
                            || secondaryCellInformation[index].neighborRanks[j] == rank
                            ) {
                                if (faceDerivatives(index, j)) {
                                    derivatives = true;
                                }
                                else {
                                    buffers = true;
                                }
                            }
                    }
            return std::pair<bool, bool>{buffers, derivatives};
        };

        std::vector<CommunicationInfo> communicationInfo;
        communicationInfo.resize(regions.size());

        if (regions.empty()) {
            for (std::size_t index = 0; index < layerView.layer.getNumberOfCells(); ++index) {
                allocate(index, false);
                allocate(index, true);
            }
        }
        else {
            for (const auto& region : regions) {
                auto startPosition = manager.position();
                for (std::size_t i = 0; i < region.size; ++i) {
                    auto index = i + region.start;
                    auto [buffers, derivatives] = useBuffersDerivatives(index, region.rank);
                    if (buffers) {
                        allocate(index, false);
                    }
                    if (derivatives) {
                        allocate(index, true);
                    }
                }
                auto endPosition = manager.position();
                auto size = endPosition - startPosition;
                communicationInfo.emplace_back(CommunicationInfo{
                    reinterpret_cast<void*>(startPosition),
                    PrecisionFromType<RealT>::MPIType,
                    size / sizeof(RealT),
                    TODO,
                    region.rank,
                    MPI_REQUEST_NULL
                });

                for (std::size_t i = 0; i < region.size; ++i) {
                    auto index = i + region.start;
                    auto [buffers, derivatives] = useBuffersDerivatives(index, region.rank);
                    if (!buffers) {
                        allocate(index, false);
                    }
                    if (!derivatives) {
                        allocate(index, true);
                    }
                }
            }
        }

        layerView.layer.setBucketSize(layerView.lts.buffersDerivatives, manager.size() + layerView.layer.getbucketSize(layerView.lts.buffersDerivatives));
    }

    template<typename Config>
    void allocateFaceDisplacements(seissol::initializers::LTSLayerRef<seissol::initializers::LTS<Config>>& layerView) {
        using RealT = typename Config::RealT;

        const auto* cellInformation = layerView.var(layerView.lts.cellInformation);
        const auto* cellMaterialData = layerView.var(layerView.lts.cellMaterialData);
        auto* faceDisplacements = layerView.var(layerView.lts.faceDisplacements);
        BucketManager manager;

        for (unsigned cell = 0; cell < layerView.layer.getNumberOfCells(); ++cell) {
            for (int face = 0; face < 4; ++face) {
                faceDisplacements[cell][face] = manager.markAllocate<RealT[Yateto<Config>::Tensor::faceDisplacement::size()]>(requiresDisplacement(
                    cellInformation[cell],
                    cellMaterialData[cell],
                    face
                ));
            }
        }

        layerView.layer.setBucketSize(layerView.lts.faceDisplacementsBuffer, manager.size());
    }

    template<typename Config>
    void touchBuckets(seissol::initializers::LTSLayerRef<seissol::initializers::LTS<Config>>& layerView) {
        auto* buffers = layerView.var(layerView.lts.buffers);
        auto* derivatives = layerView.var(layerView.lts.derivatives);
        auto* faceDisplacements = layerView.var(layerView.lts.faceDisplacements);

        auto* buffersDerivatives = layerView.var(layerView.lts.buffersDerivatives);
        auto* faceDisplacementsBuffer = layerView.var(layerView.lts.faceDisplacementsBuffer);

#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
        for (unsigned cell = 0; cell < layerView.layer.getNumberOfCells(); ++cell) {
            initBucketItem(buffers[cell], buffersDerivatives);
            initBucketItem(derivatives[cell], buffersDerivatives);
            for (int face = 0; face < 4; ++face) {
                initBucketItem(faceDisplacements[cell][face], faceDisplacementsBuffer);
            }
        }
    }

    template<typename Config>
    void setupFaceNeighbors(MemoryContainer &container, seissol::initializers::LTSLayerRef<seissol::initializers::LTS<Config>>& layerView) {
        const auto* cellInformation = layerView.var(layerView.lts.cellInformation);
        const auto* secondaryCellInformation = layerView.var(layerView.lts.secondaryCellInformation);

        auto* faceNeighbors = layerView.var(layerView.lts.faceNeighbors);

#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
        for (unsigned cell = 0; cell < layerView.layer.getNumberOfCells(); ++cell) {
            for (int face = 0; face < 4; ++face) {
                container.cluster.visitIdx(cellInformation[cell].neighborConfigIds[face], [&](auto&& ltsview){
                    if (((cellInformation[cell].ltsSetup >> (4 + face)) & 1) == 0) {
                        faceNeighbors[cell][face] = ltsview.var(ltsview.lts.buffers)[secondaryCellInformation[cell].faceNeighborIds[face]];
                    }
                    else {
                        faceNeighbors[cell][face] = ltsview.var(ltsview.lts.derivatives)[secondaryCellInformation[cell].faceNeighborIds[face]];
                    }
                });
            }
        }
    }

    template<typename Config>
    std::vector<CommunicationInfo> setupCommunication(seissol::initializers::LTSLayerRef<seissol::initializers::LTS<Config>>& layerView) {
        
    }

    void bucketsAndCommunication(MemoryContainer &container, const std::vector<ClusterLayout>& meshLayout) {

        container.cluster.visitLayers([&](auto&& layerview) {
            const auto& clusterLayout = meshLayout[container.colorMap.id(layerview.config, layerview.cluster)];
            // allocate buffers/derivatives
            if (layerview.icg == Ghost) {
                allocateTransferInfo(layerview, clusterLayout.ghostRegions);
            }
            if (layerview.icg == Copy) {
                // TODO: for all regions, setup
                allocateTransferInfo(layerview, clusterLayout.copyRegions);
            }
            if (layerview.icg == Interior) {
                // no transfer regions
                allocateTransferInfo(layerview, {});
            }

            seissol::initializer::internal::allocateFaceDisplacements(layerview);
        });

        container.cluster.allocateBuckets();

        container.cluster.visitLayers([&](auto&& layerview) {
            seissol::initializer::internal::touchBuckets(layerview);
            setupFaceNeighbors(container, layerview);
        });
    }
} // namespace seissol::initializer::internal

