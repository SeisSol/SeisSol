// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Buckets.h"
#include "Parallel/MPI.h"
#include <Common/Constants.h>
#include <Common/Real.h>
#include <Config.h>
#include <GeneratedCode/tensor.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/TimeStepping/Halo.h>
#include <Kernels/Common.h>
#include <Kernels/Precision.h>
#include <Kernels/Solver.h>
#include <Memory/Descriptor/LTS.h>
#include <Memory/Tree/Layer.h>
#include <Solver/TimeStepping/HaloCommunication.h>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <utility>
#include <vector>
#include <yateto/InitTools.h>

namespace {
using namespace seissol::initializer;
using namespace seissol::initializer::internal;

class BucketManager {
  private:
  std::size_t dataSize{0};

  public:
  real* markAllocate(std::size_t size) {
    const uintptr_t offset = this->dataSize;
    this->dataSize += size;

    // the following "hack" was copied from the MemoryManager. Add +1 to pointers to differentiate
    // from nullptr NOLINTNEXTLINE
    return reinterpret_cast<real*>(offset + 1);
  }

  [[nodiscard]] std::size_t position() const { return dataSize; }

  [[nodiscard]] std::size_t size() const { return dataSize; }
};

template <typename T>
void initBucketItem(T*& data, void* bucket, bool memsetCpu) {
  if (data != nullptr) {
    const auto ddata = reinterpret_cast<uintptr_t>(data);
    const auto offset = ddata - 1;
    auto* bucketPtr = reinterpret_cast<char*>(bucket);
    // this rather strange offset behavior is required by clang-tidy (and the reason makes sense)
    data = reinterpret_cast<T*>(bucketPtr + offset);
    if (memsetCpu) {
      std::memset(data, 0, sizeof(T));
    } else {
#ifdef ACL_DEVICE

#endif
    }
  }
}

template<typename Cfg>
std::vector<solver::RemoteCluster>
    allocateTransferInfo(Cfg cfg, LTS::Layer& layer, const std::vector<RemoteCellRegion>& regions) {
  auto* buffers = layer.var<LTS::Buffers>(cfg);
  auto* derivatives = layer.var<LTS::Derivatives>(cfg);
  auto* buffersDevice = layer.var<LTS::BuffersDevice>(cfg);
  auto* derivativesDevice = layer.var<LTS::DerivativesDevice>(cfg);
  const auto* cellInformation = layer.var<LTS::CellInformation>();
  const auto* secondaryCellInformation = layer.var<LTS::SecondaryInformation>();
  BucketManager manager;

  const auto datatype = Cfg::Precision;
  const auto typeSize = sizeOfRealType(datatype);

  const auto bufferSize = typeSize * tensor::I<Cfg>::size();
  const auto derivativeSize = typeSize * kernels::Solver<Cfg>::DerivativesSize;

  const auto allocate = [&](std::size_t index, bool useDerivatives) {
    if (useDerivatives) {
      const bool hasDerivatives = cellInformation[index].ltsSetup.hasDerivatives();
      if (hasDerivatives) {
        derivatives[index] = manager.markAllocate(derivativeSize);
        derivativesDevice[index] = derivatives[index];
      }
    } else {
      const bool hasBuffers = cellInformation[index].ltsSetup.hasBuffers();
      if (hasBuffers) {
        buffers[index] = manager.markAllocate(bufferSize);
        buffersDevice[index] = buffers[index];
      }
    }
  };

  const auto useBuffersDerivatives = [&](std::size_t index, int rank) {
    bool buffers = false;
    bool derivatives = false;

    // TODO: optimize away again
    for (std::size_t j = 0; j < Cell::NumFaces; ++j) {
      if ((secondaryCellInformation[index].rank == rank &&
           secondaryCellInformation[index].neighborRanks[j] >= 0) ||
          (secondaryCellInformation[index].neighborRanks[j] == rank &&
           secondaryCellInformation[index].rank == seissol::MPI::mpi.rank())) {
        if (cellInformation[index].ltsSetup.neighborHasDerivatives(j)) {
          derivatives = true;
        } else {
          buffers = true;
        }
      }
    }
    buffers = true;
    derivatives = true;
    return std::pair<bool, bool>{buffers, derivatives};
  };

  std::vector<solver::RemoteCluster> remoteClusters;
  remoteClusters.reserve(regions.size());

  if (regions.empty()) {
    for (std::size_t index = 0; index < layer.size(); ++index) {
      allocate(index, false);
    }
    for (std::size_t index = 0; index < layer.size(); ++index) {
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
      }
      for (std::size_t i = 0; i < region.count; ++i) {
        const auto index = i + counter;
        auto [buffers, derivatives] = useBuffersDerivatives(index, region.rank);
        if (derivatives) {
          allocate(index, true);
        }
      }
      auto endPosition = manager.position();
      auto size = endPosition - startPosition;
      assert(size % typeSize == 0);

      // NOLINTNEXTLINE
      auto* startPtr = reinterpret_cast<void*>(startPosition);

      remoteClusters.emplace_back(startPtr, size / typeSize, datatype, region.rank, region.tag);

      for (std::size_t i = 0; i < region.count; ++i) {
        const auto index = i + counter;
        auto [buffers, derivatives] = useBuffersDerivatives(index, region.rank);
        if (!buffers) {
          allocate(index, false);
        }
      }
      for (std::size_t i = 0; i < region.count; ++i) {
        const auto index = i + counter;
        auto [buffers, derivatives] = useBuffersDerivatives(index, region.rank);
        if (!derivatives) {
          allocate(index, true);
        }
      }

      counter += region.count;
    }

    assert(counter == layer.size());
  }

  layer.setEntrySize<LTS::BuffersDerivatives>(manager.size());

  return remoteClusters;
}

template<typename Cfg>
void setupBuckets(Cfg cfg, LTS::Layer& layer, std::vector<solver::RemoteCluster>& comm) {
  auto* buffers = layer.var<LTS::Buffers>(cfg);
  auto* derivatives = layer.var<LTS::Derivatives>(cfg);

  auto* buffersDerivatives = layer.var<LTS::BuffersDerivatives>(cfg);

  auto* buffersDevice = layer.var<LTS::BuffersDevice>(cfg);
  auto* derivativesDevice = layer.var<LTS::DerivativesDevice>(cfg);

  auto* buffersDerivativesDevice = layer.var<LTS::BuffersDerivatives>(cfg, AllocationPlace::Device);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (std::size_t cell = 0; cell < layer.size(); ++cell) {
    initBucketItem(buffers[cell], buffersDerivatives, true);
    initBucketItem(derivatives[cell], buffersDerivatives, true);

    assert(!layer.var<LTS::CellInformation>()[cell].ltsSetup.hasBuffers() ||
           buffers[cell] != nullptr);
    assert(!layer.var<LTS::CellInformation>()[cell].ltsSetup.hasDerivatives() ||
           derivatives[cell] != nullptr);

    if constexpr (isDeviceOn()) {
      initBucketItem(buffersDevice[cell], buffersDerivativesDevice, false);
      initBucketItem(derivativesDevice[cell], buffersDerivativesDevice, false);

      assert(!layer.var<LTS::CellInformation>()[cell].ltsSetup.hasBuffers() ||
             buffersDevice[cell] != nullptr);
      assert(!layer.var<LTS::CellInformation>()[cell].ltsSetup.hasDerivatives() ||
             derivativesDevice[cell] != nullptr);
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

template<typename Cfg>
void setupFaceNeighbors(Cfg cfg, LTS::Storage& storage, LTS::Layer& layer) {
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
      if (cellInformation[cell].faceTypes[face] != FaceType::Outflow) {
        if (faceNeighbor == StoragePosition::NullPosition) {
          if (cellInformation[cell].ltsSetup.neighborHasDerivatives(face)) {
            faceNeighbors[cell][face] = layer.var<LTS::Derivatives>(cfg)[cell];
            if constexpr (isDeviceOn()) {
              faceNeighborsDevice[cell][face] = layer.var<LTS::DerivativesDevice>(cfg)[cell];
            }
          } else {
            faceNeighbors[cell][face] = layer.var<LTS::Buffers>(cfg)[cell];
            if constexpr (isDeviceOn()) {
              faceNeighborsDevice[cell][face] = layer.var<LTS::BuffersDevice>(cfg)[cell];
            }
          }
        } else {
          if (cellInformation[cell].ltsSetup.neighborHasDerivatives(face)) {
            faceNeighbors[cell][face] = storage.lookup<LTS::Derivatives>(faceNeighbor);
            if constexpr (isDeviceOn()) {
              faceNeighborsDevice[cell][face] =
                  storage.lookup<LTS::DerivativesDevice>(faceNeighbor);
            }
          } else {
            faceNeighbors[cell][face] = storage.lookup<LTS::Buffers>(faceNeighbor);
            if constexpr (isDeviceOn()) {
              faceNeighborsDevice[cell][face] = storage.lookup<LTS::BuffersDevice>(faceNeighbor);
            }
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
    layer.wrap([&](auto cfg) {
      commInfo[layer.id()] = allocateTransferInfo(cfg, layer, layout[layer.id()].regions);
    });
  }

  storage.allocateBuckets();

  for (auto& layer : storage.leaves()) {
    layer.wrap([&](auto cfg) {
      setupBuckets(cfg, layer, commInfo[layer.id()]);
    });
  }
  for (auto& layer : storage.leaves(Ghost)) {
    layer.wrap([&](auto cfg) {
      setupFaceNeighbors(cfg, storage, layer);
    });
  }

#ifdef ACL_DEVICE

#endif

  solver::HaloCommunication communication;

  communication.resize(storage.numChildren());
  for (auto& comm : communication) {
    comm.resize(storage.numChildren());
  }

  const auto colorAdjust = [&](std::size_t color, HaloType halo) {
    auto id = storage.getColorMap().argument(color);
    id.halo = halo;
    return storage.getColorMap().colorId(id);
  };

  for (const auto& layer : storage.leaves()) {
    const auto& idInfo = layout[layer.id()].regions;
    if (layer.getIdentifier().halo == HaloType::Copy) {
      for (std::size_t i = 0; i < idInfo.size(); ++i) {
        const auto localId = colorAdjust(idInfo[i].localId, HaloType::Copy);
        const auto remoteId = colorAdjust(idInfo[i].remoteId, HaloType::Ghost);
        communication[localId][remoteId].copy.emplace_back(commInfo[layer.id()][i]);
      }
    }
    if (layer.getIdentifier().halo == HaloType::Ghost) {
      for (std::size_t i = 0; i < idInfo.size(); ++i) {
        const auto localId = colorAdjust(idInfo[i].remoteId, HaloType::Copy);
        const auto remoteId = colorAdjust(idInfo[i].localId, HaloType::Ghost);
        communication[localId][remoteId].ghost.emplace_back(commInfo[layer.id()][i]);
      }
    }
  }

  return communication;
}
} // namespace seissol::initializer::internal
