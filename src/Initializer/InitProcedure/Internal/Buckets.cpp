// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Buckets.h"

#include "Alignment.h"
#include "Common/Constants.h"
#include "Common/Real.h"
#include "Config.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/CellLocalInformation.h"
#include "Initializer/LtsSetup.h"
#include "Initializer/TimeStepping/Halo.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"
#include "Kernels/Solver.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Backmap.h"
#include "Memory/Tree/Layer.h"
#include "Solver/TimeStepping/HaloCommunication.h"

#include <array>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <utils/logger.h>
#include <vector>
#include <yateto/InitTools.h>

namespace seissol::initializer::internal {

namespace {

class BucketManager {
  private:
  std::size_t dataSize_{0};

  public:
  void align() {
    // round up by Alignment
    this->dataSize_ = ((this->dataSize_ + Alignment - 1) / Alignment) * Alignment;
  }

  real* markAllocate(std::size_t size) {
    const uintptr_t offset = this->dataSize_;
    this->dataSize_ += size;

    // the following "hack" was copied from the MemoryManager. Add +1 to pointers to differentiate
    // from nullptr NOLINTNEXTLINE
    return reinterpret_cast<real*>(offset + 1);
  }

  [[nodiscard]] std::size_t position() const { return dataSize_; }

  [[nodiscard]] std::size_t size() const { return dataSize_; }
};

template <typename T>
void initBucketItem(T*& data, void* bucket, std::size_t count, bool init) {
  if (data != nullptr) {
    const auto ddata = reinterpret_cast<uintptr_t>(data);
    const auto offset = ddata - 1;
    auto* bucketPtr = reinterpret_cast<uint8_t*>(bucket);
    // this rather strange offset behavior is required by clang-tidy (and the reason makes sense)
    data = reinterpret_cast<T*>(bucketPtr + offset);

    if (init) {
      std::memset(data, 0, count * sizeof(T));
    }
  }
}

auto useBuffersDerivatives(const LTS::Storage& storage,
                           const LTS::Layer& layer,
                           std::size_t index,
                           const RemoteCellRegion& region) {
  std::array<bool, BufferCount> bufferPresence{};

  const auto& myPrimary = layer.var<LTS::CellInformation>()[index];
  const auto& mySecondary = layer.var<LTS::SecondaryInformation>()[index];

  // what we do here: we check whether one of our neighbor cells demands derivatives from us which
  // is in the same region. i.e. we jump onto the neighboring cell (if existent), look for the same
  // face we jumped over from the other side, and check the `neighborHasDerivatives` for the
  // original cell.
  for (std::size_t j = 0; j < Cell::NumFaces; ++j) {
    if (mySecondary.faceNeighbors[j] != StoragePosition::NullPosition) {
      const auto& primary = storage.lookup<LTS::CellInformation>(mySecondary.faceNeighbors[j]);
      const auto& secondary =
          storage.lookup<LTS::SecondaryInformation>(mySecondary.faceNeighbors[j]);

      const auto isCopyGhost = secondary.halo != mySecondary.halo &&
                               mySecondary.halo != HaloType::Interior &&
                               secondary.halo != HaloType::Interior;

      const auto colorCorrect = secondary.color == region.remoteId;

      if (colorCorrect && isCopyGhost &&
          (secondary.rank == region.rank || mySecondary.rank == region.rank)) {
        for (std::size_t k = 0; k < Cell::NumFaces; ++k) {
          if (secondary.faceNeighbors[k] != StoragePosition::NullPosition) {
            const auto& reverseSecondary =
                storage.lookup<LTS::SecondaryInformation>(secondary.faceNeighbors[k]);
            if (reverseSecondary.globalId == mySecondary.globalId) {
              const auto bufferType = primary.ltsSetup.neighborBuffer(k);
              if (!myPrimary.ltsSetup.hasBuffer(bufferType)) {
                logError()
                    << "Setup error: the given buffer type was requested, but are not present.";
              }
              bufferPresence[static_cast<std::size_t>(bufferType)] = true;
            }
          }
        }
      }
    }
  }

  // if there are suspected correctness issues, enable
  // bufferPresence = true;
  return bufferPresence;
}

template <typename TypeP, typename TypeDeviceP>
struct LTSTypes {
  using Type = TypeP;
  using TypeDevice = TypeDeviceP;
};

std::vector<solver::RemoteCluster> allocateTransferInfo(
    const LTS::Storage& storage, LTS::Layer& layer, const std::vector<RemoteCellRegion>& regions) {
  const auto* cellInformation = layer.var<LTS::CellInformation>();
  BucketManager manager;

  const auto datatype = Config::Precision;
  const auto typeSize = sizeOfRealType(datatype);

  const auto bufferSize = typeSize * kernels::Solver::BuffersSize;
  const auto derivativeSize = typeSize * kernels::Solver::DerivativesSize;

  std::array<real**, BufferCount> pointers{};
  std::array<real**, BufferCount> pointersDevice{};
  std::array<std::size_t, BufferCount> sizes{};

  pointers[static_cast<std::size_t>(BufferType::Derivatives)] = layer.var<LTS::Derivatives>();
  pointers[static_cast<std::size_t>(BufferType::StepIntegrals)] = layer.var<LTS::StepIntegrals>();
  pointers[static_cast<std::size_t>(BufferType::AccumulatedIntegrals)] =
      layer.var<LTS::AccumulatedIntegrals>();

  pointersDevice[static_cast<std::size_t>(BufferType::Derivatives)] =
      layer.var<LTS::DerivativesDevice>();
  pointersDevice[static_cast<std::size_t>(BufferType::StepIntegrals)] =
      layer.var<LTS::StepIntegralsDevice>();
  pointersDevice[static_cast<std::size_t>(BufferType::AccumulatedIntegrals)] =
      layer.var<LTS::AccumulatedIntegralsDevice>();

  sizes[static_cast<std::size_t>(BufferType::Derivatives)] = derivativeSize;
  sizes[static_cast<std::size_t>(BufferType::StepIntegrals)] = bufferSize;
  sizes[static_cast<std::size_t>(BufferType::AccumulatedIntegrals)] = bufferSize;

  const auto allocate = [&](std::size_t index, BufferType type) {
    const bool hasBuffer = cellInformation[index].ltsSetup.hasBuffer(type);
    if (hasBuffer) {
      const auto bufferTypeIndex = static_cast<std::size_t>(type);
      pointers[bufferTypeIndex][index] = manager.markAllocate(sizes[bufferTypeIndex]);
      pointersDevice[bufferTypeIndex][index] = pointers[bufferTypeIndex][index];
    }
  };

  const auto allocationPass =
      [&](std::size_t counter, const RemoteCellRegion& region, bool needed, BufferType type) {
        for (std::size_t i = 0; i < region.count; ++i) {
          const auto index = i + counter;
          const auto buffers = useBuffersDerivatives(storage, layer, index, region);
          if ((needed || layer.getIdentifier().halo != HaloType::Ghost) &&
              buffers[static_cast<std::size_t>(type)] == needed) {
            allocate(index, type);
          }
        }
      };

  std::vector<solver::RemoteCluster> remoteClusters;
  remoteClusters.reserve(regions.size());

  if (regions.empty()) {
    for (std::size_t index = 0; index < layer.size(); ++index) {
      allocate(index, BufferType::AccumulatedIntegrals);
    }
    for (std::size_t index = 0; index < layer.size(); ++index) {
      allocate(index, BufferType::StepIntegrals);
    }
    for (std::size_t index = 0; index < layer.size(); ++index) {
      allocate(index, BufferType::Derivatives);
    }
  } else {
    std::size_t counter = 0;

    // allocate all to-transfer buffers/derivatives contiguously (note: region.rank)
    // (thus do non-relevant buffers before and non-relevant derivatives afterwards)
    for (const auto& region : regions) {
      allocationPass(counter, region, false, BufferType::AccumulatedIntegrals);
      allocationPass(counter, region, false, BufferType::StepIntegrals);

      // transfer allocation
      manager.align();
      auto startPosition = manager.position();
      allocationPass(counter, region, true, BufferType::StepIntegrals);
      allocationPass(counter, region, true, BufferType::AccumulatedIntegrals);
      allocationPass(counter, region, true, BufferType::Derivatives);
      manager.align();
      auto endPosition = manager.position();
      auto size = endPosition - startPosition;
      assert(size % typeSize == 0);

      // NOLINTNEXTLINE
      auto* startPtr = reinterpret_cast<void*>(startPosition);

      remoteClusters.emplace_back(startPtr, size / typeSize, datatype, region.rank, region.tag);

      allocationPass(counter, region, false, BufferType::Derivatives);

      counter += region.count;
    }

    assert(counter == layer.size());
  }

  layer.setEntrySize<LTS::Buffers>(manager.size());

  return remoteClusters;
}

void setupBuckets(LTS::Layer& layer, std::vector<solver::RemoteCluster>& comm) {
  auto* buffers = layer.var<LTS::Buffers>();
  auto* buffersDevice = layer.var<LTS::Buffers>(AllocationPlace::Device);

  const auto bufferSize = kernels::Solver::BuffersSize;
  const auto derivativeSize = kernels::Solver::DerivativesSize;

  std::array<real**, BufferCount> pointers{};
  std::array<real**, BufferCount> pointersDevice{};
  std::array<std::size_t, BufferCount> sizes{};

  pointers[static_cast<std::size_t>(BufferType::Derivatives)] = layer.var<LTS::Derivatives>();
  pointers[static_cast<std::size_t>(BufferType::StepIntegrals)] = layer.var<LTS::StepIntegrals>();
  pointers[static_cast<std::size_t>(BufferType::AccumulatedIntegrals)] =
      layer.var<LTS::AccumulatedIntegrals>();

  pointersDevice[static_cast<std::size_t>(BufferType::Derivatives)] =
      layer.var<LTS::DerivativesDevice>();
  pointersDevice[static_cast<std::size_t>(BufferType::StepIntegrals)] =
      layer.var<LTS::StepIntegralsDevice>();
  pointersDevice[static_cast<std::size_t>(BufferType::AccumulatedIntegrals)] =
      layer.var<LTS::AccumulatedIntegralsDevice>();

  sizes[static_cast<std::size_t>(BufferType::Derivatives)] = derivativeSize;
  sizes[static_cast<std::size_t>(BufferType::StepIntegrals)] = bufferSize;
  sizes[static_cast<std::size_t>(BufferType::AccumulatedIntegrals)] = bufferSize;

#pragma omp parallel for schedule(static)
  for (std::size_t cell = 0; cell < layer.size(); ++cell) {
    for (std::size_t type = 0; type < pointers.size(); ++type) {
      initBucketItem(pointers[type][cell], buffers, sizes[type], true);
      assert(!layer.var<LTS::CellInformation>()[cell].ltsSetup.hasBuffer(
                 static_cast<BufferType>(type)) ||
             pointers[type][cell] != nullptr || layer.getIdentifier().halo == HaloType::Ghost);
    }

    if constexpr (isDeviceOn()) {
      for (std::size_t type = 0; type < pointersDevice.size(); ++type) {
        initBucketItem(pointersDevice[type][cell], buffersDevice, bufferSize, false);
        assert(!layer.var<LTS::CellInformation>()[cell].ltsSetup.hasBuffer(
                   static_cast<BufferType>(type)) ||
               pointersDevice[type][cell] != nullptr ||
               layer.getIdentifier().halo == HaloType::Ghost);
      }
    }
  }

  for (auto& info : comm) {
    const auto offset = reinterpret_cast<intptr_t>(info.data);
    uint8_t* base = nullptr;
    if constexpr (isDeviceOn()) {
      base = reinterpret_cast<uint8_t*>(buffersDevice);
    } else {
      base = reinterpret_cast<uint8_t*>(buffers);
    }
    info.data = reinterpret_cast<void*>(base + offset);
  }
}

void setupFaceNeighbors(LTS::Storage& storage, LTS::Layer& layer) {
  const auto* cellInformation = layer.var<LTS::CellInformation>();
  const auto* secondaryCellInformation = layer.var<LTS::SecondaryInformation>();

  auto* faceNeighbors = layer.var<LTS::FaceNeighbors>();
  auto* faceNeighborsDevice = layer.var<LTS::FaceNeighborsDevice>();

#pragma omp parallel for schedule(static)
  for (std::size_t cell = 0; cell < layer.size(); ++cell) {
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      const auto& faceNeighbor = secondaryCellInformation[cell].faceNeighbors[face];
      if (getBCType(cellInformation[cell].faceTypes[face]) != BCType::External) {
        const auto type = cellInformation[cell].ltsSetup.neighborBuffer(face);

        const auto setFaceNeighbors = [&](auto types) {
          using Types = decltype(types);

          if (faceNeighbor == StoragePosition::NullPosition) {
            faceNeighbors[cell][face] = layer.var<typename Types::Type>()[cell];
            if constexpr (isDeviceOn()) {
              faceNeighborsDevice[cell][face] = layer.var<typename Types::TypeDevice>()[cell];
            }
          } else {
            faceNeighbors[cell][face] = storage.lookup<typename Types::Type>(faceNeighbor);
            if constexpr (isDeviceOn()) {
              faceNeighborsDevice[cell][face] =
                  storage.lookup<typename Types::TypeDevice>(faceNeighbor);
            }
          }
        };

        if (type == BufferType::Derivatives) {
          setFaceNeighbors(LTSTypes<LTS::Derivatives, LTS::DerivativesDevice>{});
        }
        if (type == BufferType::StepIntegrals) {
          setFaceNeighbors(LTSTypes<LTS::StepIntegrals, LTS::StepIntegralsDevice>{});
        }
        if (type == BufferType::AccumulatedIntegrals) {
          setFaceNeighbors(LTSTypes<LTS::AccumulatedIntegrals, LTS::AccumulatedIntegralsDevice>{});
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

solver::HaloCommunication bucketsAndCommunication(LTS::Storage& storage, const MeshLayout& layout) {
  std::vector<std::vector<solver::RemoteCluster>> commInfo(storage.getColorMap().size());

  for (auto& layer : storage.leaves()) {
    commInfo[layer.id()] = allocateTransferInfo(storage, layer, layout[layer.id()].regions);
  }

  storage.allocateBuckets();

  for (auto& layer : storage.leaves()) {
    setupBuckets(layer, commInfo[layer.id()]);
  }
  for (auto& layer : storage.leaves(Ghost)) {
    setupFaceNeighbors(storage, layer);
  }

  solver::HaloCommunication communication;

  communication.resize(commInfo.size());
  for (auto& comm : communication) {
    comm.resize(commInfo.size());
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
