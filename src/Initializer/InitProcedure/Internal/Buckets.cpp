// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Buckets.h"
#include "Parallel/MPI.h"
#include <Common/Real.h>
#include <Config.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/MemoryManager.h>
#include <Initializer/TimeStepping/Halo.h>
#include <Kernels/Common.h>
#include <Kernels/Precision.h>
#include <Memory/Descriptor/LTS.h>
#include <Memory/Tree/LTSTree.h>
#include <Memory/Tree/Layer.h>
#include <Solver/TimeStepping/HaloCommunication.h>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <limits>
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
  real* markAllocate(std::size_t size, bool allocate) {
    if (allocate) {
      const uintptr_t offset = dataSize;
      dataSize += size;

      // the following "hack" was copied from the MemoryManager. Add +1 to pointers to differentiate
      // from nullptr NOLINTNEXTLINE
      return reinterpret_cast<real*>(offset + 1);
    } else {
      return nullptr;
    }
  }

  [[nodiscard]] std::size_t position() const { return dataSize + 1; }

  [[nodiscard]] std::size_t size() const { return dataSize; }
};

template <typename T>
void initBucketItem(T*& data, void* bucket, bool memset) {
  if (data != nullptr) {
    const auto ddata = reinterpret_cast<uintptr_t>(data);
    const auto offset = ddata - 1;
    auto* bucketPtr = reinterpret_cast<char*>(bucket);
    // this rather strange offset behavior is required by clang-tidy (and the reason makes sense)
    data = reinterpret_cast<T*>(bucketPtr + offset);
    if (memset) {
      std::memset(data, 0, sizeof(T));
    }
  }
}

std::vector<solver::RemoteCluster>
    allocateTransferInfo(LTS::Layer& layer, const std::vector<RemoteCellRegion>& regions) {
  auto* buffers = layer.var<LTS::Buffers>();
  auto* derivatives = layer.var<LTS::Derivatives>();
  auto* buffersDevice = layer.var<LTS::BuffersDevice>();
  auto* derivativesDevice = layer.var<LTS::DerivativesDevice>();
  const auto* cellInformation = layer.var<LTS::CellInformation>();
  const auto* secondaryCellInformation = layer.var<LTS::SecondaryInformation>();
  BucketManager manager;

  const auto datatype = Config::Precision;
  const auto typeSize = sizeOfRealType(datatype);

  const auto bufferSize = typeSize * yateto::computeFamilySize<tensor::dQ>();
  const auto derivativeSize = typeSize * tensor::I::size();

  const auto allocate = [&](std::size_t index, bool useDerivatives) {
    const bool hasBuffers = cellInformation[index].ltsSetup.hasBuffers();
    const bool hasDerivatives = cellInformation[index].ltsSetup.hasDerivatives();

    if (useDerivatives) {
      derivatives[index] = manager.markAllocate(bufferSize, hasDerivatives);
      derivativesDevice[index] = derivatives[index];
    } else {
      buffers[index] = manager.markAllocate(derivativeSize, hasBuffers);
      buffersDevice[index] = buffers[index];
    }
  };

  const auto useBuffersDerivatives = [&](std::size_t index, int rank) {
    bool buffers = false;
    bool derivatives = false;
    for (std::size_t j = 0; j < Cell::NumFaces; ++j) {
      if ((secondaryCellInformation[index].rank == rank &&
           secondaryCellInformation[index].neighborRanks[j] >= 0) ||
          secondaryCellInformation[index].neighborRanks[j] == rank) {
        if (cellInformation[index].ltsSetup.neighborHasDerivatives(j)) {
          derivatives = true;
        } else {
          buffers = true;
        }
      }
    }
    return std::pair<bool, bool>{buffers, derivatives};
  };

  std::vector<solver::RemoteCluster> remoteClusters;
  remoteClusters.reserve(regions.size());

  if (regions.empty()) {
    for (std::size_t index = 0; index < layer.size(); ++index) {
      allocate(index, false);
      allocate(index, true);
    }
  } else {
    std::size_t counter = 0;

    // allocate all to-transfer buffers/derivatives first (note: region.rank)
    for (const auto& region : regions) {
      auto startPosition = manager.position();

      for (std::size_t i = 0; i < region.count; ++i) {
        const auto index = i + counter;
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
      assert(size % typeSize == 0);

      // NOLINTNEXTLINE
      auto* startPtr = reinterpret_cast<void*>(startPosition);

      remoteClusters.emplace_back(startPtr, size / typeSize, datatype, region.tag, region.rank);

      counter += region.count;
    }

    // allocate local (non-transfer) buffers/derivatives
    counter = 0;
    for (const auto& region : regions) {
      for (std::size_t i = 0; i < region.count; ++i) {
        auto index = i + counter;
        auto [buffers, derivatives] = useBuffersDerivatives(index, seissol::MPI::mpi.rank());
        if (!buffers) {
          allocate(index, false);
        }
        if (!derivatives) {
          allocate(index, true);
        }
      }

      counter += region.count;
    }
  }

  layer.setEntrySize<LTS::BuffersDerivatives>(manager.size() +
                                              layer.getEntrySize<LTS::BuffersDerivatives>());

  return remoteClusters;
}

void setupBuckets(LTS::Layer& layer, std::vector<solver::RemoteCluster>& comm) {
  auto* buffers = layer.var<LTS::Buffers>();
  auto* derivatives = layer.var<LTS::Derivatives>();

  auto* buffersDerivatives = layer.var<LTS::BuffersDerivatives>();

  auto* buffersDevice = layer.var<LTS::BuffersDevice>();
  auto* derivativesDevice = layer.var<LTS::DerivativesDevice>();

  auto* buffersDerivativesDevice = layer.var<LTS::BuffersDerivatives>(AllocationPlace::Device);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (std::size_t cell = 0; cell < layer.size(); ++cell) {
    initBucketItem(buffers[cell], buffersDerivatives, true);
    initBucketItem(derivatives[cell], buffersDerivatives, true);

    if constexpr (isDeviceOn()) {
      initBucketItem(buffersDevice[cell], buffersDerivativesDevice, false);
      initBucketItem(derivativesDevice[cell], buffersDerivativesDevice, false);
    }
  }

  for (auto& info : comm) {
    if constexpr (isDeviceOn()) {
      info.data = reinterpret_cast<intptr_t>(info.data) + buffersDerivativesDevice;
    } else {
      info.data = reinterpret_cast<intptr_t>(info.data) + buffersDerivatives;
    }
  }
}

void setupFaceNeighbors(LTS::Storage& storage, LTS::Layer& layer) {
  const auto* cellInformation = layer.var<LTS::CellInformation>();
  const auto* secondaryCellInformation = layer.var<LTS::SecondaryInformation>();

  auto* faceNeighbors = layer.var<LTS::FaceNeighbors>();
  auto* faceNeighborsDevice = layer.var<LTS::FaceNeighborsDevice>();

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (std::size_t cell = 0; cell < layer.size(); ++cell) {
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      const auto& faceNeighbor = secondaryCellInformation[cell].faceNeighbors[face];
      if (faceNeighbor != StoragePosition::NullPosition) {
        if (cellInformation[cell].ltsSetup.neighborHasDerivatives(face)) {
          faceNeighbors[cell][face] = storage.lookup<LTS::Derivatives>(faceNeighbor);
          if constexpr (isDeviceOn()) {
            faceNeighborsDevice[cell][face] = storage.lookup<LTS::DerivativesDevice>(faceNeighbor);
          }
        } else {
          faceNeighbors[cell][face] = storage.lookup<LTS::Buffers>(faceNeighbor);
          if constexpr (isDeviceOn()) {
            faceNeighborsDevice[cell][face] = storage.lookup<LTS::BuffersDevice>(faceNeighbor);
          }
        }

        assert(faceNeighbors[cell][face] != nullptr);
        if constexpr (isDeviceOn()) {
          assert(faceNeighborsDevice[cell][face] != nullptr);
        }
      }
    }
  }
}
} // namespace

namespace seissol::initializer::internal {
solver::HaloCommunication bucketsAndCommunication(LTS::Storage& storage, const MeshLayout& layout) {
  std::vector<std::vector<solver::RemoteCluster>> commInfo(storage.numChildren());

  for (auto& layer : storage.leaves()) {
    commInfo[layer.id()] = allocateTransferInfo(layer, layout[layer.id()].regions);
  }

  storage.allocateBuckets();

  for (auto& layer : storage.leaves()) {
    setupBuckets(layer, commInfo[layer.id()]);
  }
  for (auto& layer : storage.leaves()) {
    setupFaceNeighbors(storage, layer);
  }

  // TODO: refactor. the following code is far too complicated for what it wants to do.
  std::vector<std::vector<std::vector<solver::RemoteCluster>>> copy;
  std::vector<std::vector<std::vector<solver::RemoteCluster>>> ghost;

  solver::HaloCommunication communication;

  communication.resize(storage.numChildren());
  copy.resize(storage.numChildren());
  ghost.resize(storage.numChildren());
  for (auto& commGhost : communication) {
    commGhost.resize(storage.numChildren());
  }
  for (auto& commGhost : copy) {
    commGhost.resize(storage.numChildren());
  }
  for (auto& commGhost : ghost) {
    commGhost.resize(storage.numChildren());
  }

  for (const auto& layer : storage.leaves()) {
    if (layer.getIdentifier().halo == HaloType::Copy) {
      const auto& idInfo = layout[layer.id()].regions;
      for (std::size_t i = 0; i < idInfo.size(); ++i) {
        copy[idInfo[i].localId][idInfo[i].remoteId].emplace_back(commInfo[layer.id()][i]);
      }
    }
    if (layer.getIdentifier().halo == HaloType::Ghost) {
      const auto& idInfo = layout[layer.id()].regions;
      for (std::size_t i = 0; i < idInfo.size(); ++i) {
        ghost[idInfo[i].localId][idInfo[i].remoteId].emplace_back(commInfo[layer.id()][i]);
      }
    }
  }

  for (std::size_t i = 0; i < storage.numChildren(); ++i) {
    for (std::size_t j = 0; j < storage.numChildren(); ++j) {
      assert(copy[i][j].size() == ghost[i][j].size());
      communication[i][j].reserve(copy[i][j].size());
      for (std::size_t k = 0; k < copy[i][j].size(); ++k) {
        communication[i][j].emplace_back(solver::RemoteClusterPair{copy[i][j][k], ghost[i][j][k]});
      }
    }
  }

  return communication;
}
} // namespace seissol::initializer::internal
