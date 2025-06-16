// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Buckets.h"
#include "Parallel/MPI.h"
#include <Config.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/InitProcedure/Internal/MeshLayout.h>
#include <Initializer/MemoryManager.h>
#include <Kernels/Precision.h>
#include <Memory/Descriptor/LTS.h>
#include <Memory/MemoryContainer.h>
#include <Memory/Tree/Layer.h>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <type_traits>
#include <utility>
#include <vector>
#include <yateto/InitTools.h>

namespace {
using namespace seissol::initializer;
using namespace seissol::initializer::internal;

class BucketManager {
  private:
  std::size_t dataSize = 0;

  public:
  template <typename T>
  std::remove_extent_t<T>* markAllocate(bool allocate) {
    if (allocate) {
      const uintptr_t offset = dataSize;
      dataSize += sizeof(T);
      return reinterpret_cast<std::remove_extent_t<T>*>(offset + 1);
    } else {
      return nullptr;
    }
  }

  [[nodiscard]] std::size_t position() const { return dataSize + 1; }

  [[nodiscard]] std::size_t size() const { return dataSize; }
};

template <typename T>
void initBucketItem(T*& data, void* bucket) {
  if (data != nullptr) {
    const auto ddata = reinterpret_cast<uintptr_t>(data);
    const auto offset = ddata - 1;
    const auto bucketPos = reinterpret_cast<uintptr_t>(bucket);
    // this rather strange offset behavior is required by clang-tidy (and the reason makes sense)
    data += bucketPos + offset - ddata;
    std::memset(data, 0, sizeof(T));
  }
}

std::vector<CommunicationInfo>
    allocateTransferInfo(Layer& layer, LTS& lts, const std::vector<TransferRegion>& regions) {
  auto* buffers = layer.var(lts.buffers);
  auto* derivatives = layer.var(lts.derivatives);
  const auto* cellInformation = layer.var(lts.cellInformation);
  const auto* secondaryCellInformation = layer.var(lts.secondaryInformation);
  BucketManager manager;

  auto allocate = [&](std::size_t index, bool useDerivatives) {
    const bool hasBuffers = (cellInformation[index].ltsSetup >> 8) != 0;
    const bool hasDerivatives = (cellInformation[index].ltsSetup >> 9) != 0;
    if (useDerivatives) {
      derivatives[index] = manager.markAllocate<LTS::DerivativeT>(hasDerivatives);
    } else {
      buffers[index] = manager.markAllocate<LTS::BufferT>(hasBuffers);
    }
  };

  auto faceDerivatives = [&](std::size_t index, int face) {
    return (cellInformation[index].ltsSetup >> (4 + face)) != 0;
  };

  auto useBuffersDerivatives = [&](std::size_t index, int rank) {
    bool buffers = false;
    bool derivatives = false;
    for (int j = 0; j < 4; ++j) {
      if ((secondaryCellInformation[index].rank == rank &&
           secondaryCellInformation[index].neighborRanks[j] >= 0) ||
          secondaryCellInformation[index].neighborRanks[j] == rank) {
        if (faceDerivatives(index, j)) {
          derivatives = true;
        } else {
          buffers = true;
        }
      }
    }
    return std::pair<bool, bool>{buffers, derivatives};
  };

  std::vector<CommunicationInfo> communicationInfo;
  communicationInfo.resize(regions.size());

  if (regions.empty()) {
    for (std::size_t index = 0; index < layer.size(); ++index) {
      allocate(index, false);
      allocate(index, true);
    }
  } else {
    for (const auto& region : regions) {
      auto startPosition = manager.position();

      // allocate all to-transfer buffers/derivatives first (note: region.rank)
      for (std::size_t i = 0; i < region.size; ++i) {
        const auto index = i + region.start;
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
      assert(size % sizeof(real) == 0);
      communicationInfo.emplace_back(CommunicationInfo{
          nullptr, startPosition, Config::Precision, size / sizeof(real), region.tag, region.rank});

      // allocate local (non-transfer) buffers/derivatives
      for (std::size_t i = 0; i < region.size; ++i) {
        auto index = i + region.start;
        auto [buffers, derivatives] = useBuffersDerivatives(index, seissol::MPI::mpi.rank());
        if (!buffers) {
          allocate(index, false);
        }
        if (!derivatives) {
          allocate(index, true);
        }
      }
    }
  }

  layer.setEntrySize(lts.buffersDerivatives,
                     manager.size() + layer.getEntrySize(lts.buffersDerivatives));
  return communicationInfo;
}

void allocateFaceDisplacements(Layer& layer, LTS& lts) {
  const auto* cellInformation = layer.var(lts.cellInformation);
  const auto* cellMaterialData = layer.var(lts.material);
  auto* faceDisplacements = layer.var(lts.faceDisplacements);
  BucketManager manager;

  for (unsigned cell = 0; cell < layer.size(); ++cell) {
    for (int face = 0; face < 4; ++face) {
      faceDisplacements[cell][face] =
          manager.markAllocate<LTS::FaceDisplacementT>(seissol::initializer::requiresDisplacement(
              cellInformation[cell], cellMaterialData[cell], face));
    }
  }

  layer.setEntrySize(lts.faceDisplacementsBuffer, manager.size());
}

void setupBuckets(Layer& layer, LTS& lts, std::vector<CommunicationInfo>& comm) {
  auto* buffers = layer.var(lts.buffers);
  auto* derivatives = layer.var(lts.derivatives);
  auto* faceDisplacements = layer.var(lts.faceDisplacements);

  auto* buffersDerivatives = layer.var(lts.buffersDerivatives);
  auto* faceDisplacementsBuffer = layer.var(lts.faceDisplacementsBuffer);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (unsigned cell = 0; cell < layer.size(); ++cell) {
    initBucketItem(buffers[cell], buffersDerivatives);
    initBucketItem(derivatives[cell], buffersDerivatives);
    for (int face = 0; face < 4; ++face) {
      initBucketItem(faceDisplacements[cell][face], faceDisplacementsBuffer);
    }
  }

  for (auto& info : comm) {
    info.buffer = buffersDerivatives + info.offset;
  }
}

void setupFaceNeighbors(memory::MemoryContainer& container, Layer& layer, LTS& lts) {
  const auto* cellInformation = layer.var(lts.cellInformation);
  const auto* secondaryCellInformation = layer.var(lts.secondaryInformation);

  auto* faceNeighbors = layer.var(lts.faceNeighbors);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (unsigned cell = 0; cell < layer.size(); ++cell) {
    for (int face = 0; face < 4; ++face) {
      if (((cellInformation[cell].ltsSetup >> (4 + face)) & 1) == 0) {
        faceNeighbors[cell][face] = reinterpret_cast<real*>(
            layer.var(lts.buffers)[secondaryCellInformation[cell].faceNeighborIds[face]]);
      } else {
        faceNeighbors[cell][face] = reinterpret_cast<real*>(
            layer.var(lts.derivatives)[secondaryCellInformation[cell].faceNeighborIds[face]]);
      }
    }
  }
}
} // namespace

namespace seissol::initializer::internal {
void bucketsAndCommunication(memory::MemoryContainer& container,
                             const std::vector<ClusterLayout>& meshLayout) {

  std::vector<std::vector<CommunicationInfo>> commInfo(container.volume.numChildren());
  for (auto& layer : container.volume.leaves()) {
    const auto& clusterLayout = meshLayout[container.colorMap.colorId(layer.getIdentifier())];
    // allocate buffers/derivatives
    if (layer.getIdentifier().halo == HaloType::Ghost) {
      commInfo[layer.id()] = allocateTransferInfo(layer, container.wpdesc, clusterLayout.regions);
    }
    if (layer.getIdentifier().halo == HaloType::Copy) {
      // TODO: for all regions, setup
      commInfo[layer.id()] = allocateTransferInfo(layer, container.wpdesc, clusterLayout.regions);
    }
    if (layer.getIdentifier().halo == HaloType::Interior) {
      // no transfer regions
      commInfo[layer.id()] = allocateTransferInfo(layer, container.wpdesc, {});
    }

    allocateFaceDisplacements(layer, container.wpdesc);
  }

  container.volume.allocateBuckets();

  for (auto& layer : container.volume.leaves()) {
    setupBuckets(layer, container.wpdesc, commInfo[layer.id()]);
    setupFaceNeighbors(container, layer, container.wpdesc);
  }
}
} // namespace seissol::initializer::internal
