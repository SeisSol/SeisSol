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
#include "Initializer/TimeStepping/Halo.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Backmap.h"
#include "Memory/Tree/Layer.h"
#include "Solver/TimeStepping/HaloCommunication.h"

#include <cassert>
#include <cstdint>
#include <cstring>
#include <utility>
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
  bool buffers = false;
  bool derivatives = false;

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
              if (primary.ltsSetup.neighborHasDerivatives(k)) {
                if (!myPrimary.ltsSetup.hasDerivatives()) {
                  logError() << "Setup error: derivatives were requested, but are not present.";
                }
                derivatives = true;
              } else {
                if (!myPrimary.ltsSetup.hasBuffers()) {
                  logError() << "Setup error: buffers were requested, but are not present.";
                }
                buffers = true;
              }
            }
          }
        }
      }
    }
  }

  // if there are suspected correctness issues, enable
  // buffers = true;
  // derivatives = true;
  return std::pair<bool, bool>{buffers, derivatives};
}

std::vector<solver::RemoteCluster> allocateTransferInfo(
    const LTS::Storage& storage, LTS::Layer& layer, const std::vector<RemoteCellRegion>& regions) {
  auto* buffers = layer.var<LTS::Buffers>();
  auto* derivatives = layer.var<LTS::Derivatives>();
  auto* buffersDevice = layer.var<LTS::BuffersDevice>();
  auto* derivativesDevice = layer.var<LTS::DerivativesDevice>();
  const auto* cellInformation = layer.var<LTS::CellInformation>();
  BucketManager manager;

  const auto datatype = Config::Precision;
  const auto typeSize = sizeOfRealType(datatype);

  const auto bufferSize = typeSize * tensor::I::size();
  const auto derivativeSize = typeSize * yateto::computeFamilySize<tensor::dQ>();

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

  const auto allocationPass =
      [&](std::size_t counter, const RemoteCellRegion& region, bool needed, bool useDerivatives) {
        for (std::size_t i = 0; i < region.count; ++i) {
          const auto index = i + counter;
          const auto [buffers, derivatives] = useBuffersDerivatives(storage, layer, index, region);
          if (needed || layer.getIdentifier().halo != HaloType::Ghost) {
            if (!useDerivatives && buffers == needed) {
              allocate(index, false);
            }
            if (useDerivatives && derivatives == needed) {
              allocate(index, true);
            }
          }
        }
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

    // allocate all to-transfer buffers/derivatives contiguously (note: region.rank)
    // (thus do non-relevant buffers before and non-relevant derivatives afterwards)
    for (const auto& region : regions) {
      allocationPass(counter, region, false, false);

      // transfer allocation
      manager.align();
      auto startPosition = manager.position();
      allocationPass(counter, region, true, false);
      allocationPass(counter, region, true, true);
      manager.align();
      auto endPosition = manager.position();
      auto size = endPosition - startPosition;
      assert(size % typeSize == 0);

      // NOLINTNEXTLINE
      auto* startPtr = reinterpret_cast<void*>(startPosition);

      remoteClusters.emplace_back(startPtr, size / typeSize, datatype, region.rank, region.tag);

      allocationPass(counter, region, false, true);

      counter += region.count;
    }

    assert(counter == layer.size());
  }

  layer.setEntrySize<LTS::BuffersDerivatives>(manager.size());

  return remoteClusters;
}

void setupBuckets(LTS::Layer& layer, std::vector<solver::RemoteCluster>& comm) {
  auto* buffers = layer.var<LTS::Buffers>();
  auto* derivatives = layer.var<LTS::Derivatives>();

  auto* buffersDerivatives = layer.var<LTS::BuffersDerivatives>();

  auto* buffersDevice = layer.var<LTS::BuffersDevice>();
  auto* derivativesDevice = layer.var<LTS::DerivativesDevice>();

  auto* buffersDerivativesDevice = layer.var<LTS::BuffersDerivatives>(AllocationPlace::Device);

  const auto bufferSize = tensor::I::size();
  const auto derivativeSize = yateto::computeFamilySize<tensor::dQ>();

#pragma omp parallel for schedule(static)
  for (std::size_t cell = 0; cell < layer.size(); ++cell) {
    initBucketItem(buffers[cell], buffersDerivatives, bufferSize, true);
    initBucketItem(derivatives[cell], buffersDerivatives, derivativeSize, true);

    assert(!layer.var<LTS::CellInformation>()[cell].ltsSetup.hasBuffers() ||
           buffers[cell] != nullptr || layer.getIdentifier().halo == HaloType::Ghost);
    assert(!layer.var<LTS::CellInformation>()[cell].ltsSetup.hasDerivatives() ||
           derivatives[cell] != nullptr || layer.getIdentifier().halo == HaloType::Ghost);

    if constexpr (isDeviceOn()) {
      initBucketItem(buffersDevice[cell], buffersDerivativesDevice, bufferSize, false);
      initBucketItem(derivativesDevice[cell], buffersDerivativesDevice, derivativeSize, false);

      assert(!layer.var<LTS::CellInformation>()[cell].ltsSetup.hasBuffers() ||
             buffersDevice[cell] != nullptr || layer.getIdentifier().halo == HaloType::Ghost);
      assert(!layer.var<LTS::CellInformation>()[cell].ltsSetup.hasDerivatives() ||
             derivativesDevice[cell] != nullptr || layer.getIdentifier().halo == HaloType::Ghost);
    }
  }

  for (auto& info : comm) {
    const auto offset = reinterpret_cast<intptr_t>(info.data);
    uint8_t* base = nullptr;
    if constexpr (isDeviceOn()) {
      base = reinterpret_cast<uint8_t*>(buffersDerivativesDevice);
    } else {
      base = reinterpret_cast<uint8_t*>(buffersDerivatives);
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
      if (cellInformation[cell].faceTypes[face] != FaceType::Outflow) {
        if (faceNeighbor == StoragePosition::NullPosition) {
          if (cellInformation[cell].ltsSetup.neighborHasDerivatives(face)) {
            faceNeighbors[cell][face] = layer.var<LTS::Derivatives>()[cell];
            if constexpr (isDeviceOn()) {
              faceNeighborsDevice[cell][face] = layer.var<LTS::DerivativesDevice>()[cell];
            }
          } else {
            faceNeighbors[cell][face] = layer.var<LTS::Buffers>()[cell];
            if constexpr (isDeviceOn()) {
              faceNeighborsDevice[cell][face] = layer.var<LTS::BuffersDevice>()[cell];
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
