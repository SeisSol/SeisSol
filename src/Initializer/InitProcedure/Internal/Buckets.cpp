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
#include "Kernels/FaceTraces.h"
#include "Kernels/Precision.h"
#include "Kernels/Solver.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Backmap.h"
#include "Memory/Tree/Layer.h"
#include "Solver/TimeStepping/HaloCommunication.h"

#include <algorithm>
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

/// How many reals a buffer of the given type takes on a face.
constexpr std::size_t bufferSize(BufferType type) {
  return type == BufferType::Derivatives ? kernels::Solver::FaceDerivativesSize
                                         : kernels::Solver::FaceTraceSize;
}

/// Whether the neighbor across `face` of the cell `index` of `layer` is one of the cells `region`
/// pairs this layer with on another rank, so that what the cell writes on that face crosses to that
/// rank (a copy cell) or what the neighbor writes there comes from it (a ghost cell).
bool crossesTo(const LTS::Storage& storage,
               const LTS::Layer& layer,
               std::size_t index,
               std::size_t face,
               const RemoteCellRegion& region) {
  const auto& mySecondary = layer.var<LTS::SecondaryInformation>()[index];
  const auto& neighbor = mySecondary.faceNeighbors[face];
  if (neighbor == StoragePosition::NullPosition) {
    return false;
  }
  const auto& secondary = storage.lookup<LTS::SecondaryInformation>(neighbor);
  const auto isCopyGhost = secondary.halo != mySecondary.halo &&
                           mySecondary.halo != HaloType::Interior &&
                           secondary.halo != HaloType::Interior;
  const auto colorCorrect = secondary.color == region.remoteId;
  return colorCorrect && isCopyGhost &&
         (secondary.rank == region.rank || mySecondary.rank == region.rank);
}

/// Checks that what a cell writes on a face crossing to another rank is what the neighbor there
/// reads, both as this rank knows them.
void checkCrossing(const LTS::Storage& storage,
                   const LTS::Layer& layer,
                   std::size_t index,
                   std::size_t face) {
  const auto& primary = layer.var<LTS::CellInformation>()[index];
  const auto& secondary = layer.var<LTS::SecondaryInformation>()[index];
  const auto& neighbor = storage.lookup<LTS::CellInformation>(secondary.faceNeighbors[face]);
  const auto neighborFace = primary.faceRelations[face][0];
  if (primary.ltsSetup.faceBuffer(face) != neighbor.ltsSetup.neighborBuffer(neighborFace)) {
    logError() << "Setup error: a cell writes another buffer on a face than its neighbor reads.";
  }
}

std::vector<solver::RemoteCluster> allocateTransferInfo(
    const LTS::Storage& storage, LTS::Layer& layer, const std::vector<RemoteCellRegion>& regions) {
  const auto* cellInformation = layer.var<LTS::CellInformation>();
  auto* faceTraces = layer.var<LTS::FaceTraces>();
  auto* faceBuffers = layer.var<LTS::FaceBuffers>();
  auto* faceTracesDevice = layer.var<LTS::FaceTracesDevice>();
  auto* faceBuffersDevice = layer.var<LTS::FaceBuffersDevice>();
  // a ghost cell computes nothing here: it has what its neighbors read of it, and no trace of its
  // own
  const bool ghost = layer.getIdentifier().halo == HaloType::Ghost;
  BucketManager manager;

  const auto datatype = Config::Precision;
  const auto typeSize = sizeOfRealType(datatype);

  for (std::size_t index = 0; index < layer.size(); ++index) {
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      faceTraces[index][face] = nullptr;
      faceBuffers[index][face] = nullptr;
      faceTracesDevice[index][face] = nullptr;
      faceBuffersDevice[index][face] = nullptr;
    }
  }

  // what a cell writes on a face for the neighbor across it; the trace of its time integral is
  // that buffer where the neighbor runs at the same time step
  const auto allocateBuffer = [&](std::size_t index, std::size_t face) {
    const auto type = cellInformation[index].ltsSetup.faceBuffer(face).value();
    auto* offset = manager.markAllocate(bufferSize(type) * typeSize);
    faceBuffers[index][face] = offset;
    faceBuffersDevice[index][face] = offset;
    if (type == BufferType::StepIntegrals && !ghost) {
      faceTraces[index][face] = offset;
      faceTracesDevice[index][face] = offset;
    }
  };
  const auto allocateTrace = [&](std::size_t index, std::size_t face) {
    auto* offset = manager.markAllocate(kernels::Solver::FaceTraceSize * typeSize);
    faceTraces[index][face] = offset;
    faceTracesDevice[index][face] = offset;
  };

  // Everything of the cells [begin, end) that stays on this rank: the traces of a cell, and what it
  // writes on a face whose neighbor is not one of the cells `crossing` says cross to another rank.
  const auto localPass = [&](std::size_t begin, std::size_t end, const auto& crossing) {
    if (ghost) {
      return;
    }
    for (std::size_t index = begin; index < end; ++index) {
      const auto& ltsSetup = cellInformation[index].ltsSetup;
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        const auto type = ltsSetup.faceBuffer(face);
        const bool crosses = crossing(index, face);
        if (type.has_value() && !crosses) {
          allocateBuffer(index, face);
        }
        if (type != BufferType::StepIntegrals) {
          allocateTrace(index, face);
        }
      }
    }
  };
  // what crosses: what a copy cell writes for a neighbor on the other rank, or what a neighbor
  // there writes for a cell here, which this rank keeps as the ghost cell's
  const auto transferPass = [&](std::size_t begin, std::size_t end, const auto& crossing) {
    for (std::size_t index = begin; index < end; ++index) {
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        if (crossing(index, face)) {
          checkCrossing(storage, layer, index, face);
          if (cellInformation[index].ltsSetup.faceBuffer(face).has_value()) {
            allocateBuffer(index, face);
          }
        }
      }
    }
  };

  std::vector<solver::RemoteCluster> remoteClusters;
  remoteClusters.reserve(regions.size());

  if (regions.empty()) {
    localPass(0, layer.size(), [](std::size_t, std::size_t) { return false; });
  } else {
    std::size_t counter = 0;

    // allocate all to-transfer buffers contiguously per region (note: region.rank)
    for (const auto& region : regions) {
      const auto crossing = [&](std::size_t index, std::size_t face) {
        return crossesTo(storage, layer, index, face, region);
      };
      localPass(counter, counter + region.count, crossing);

      // transfer allocation
      manager.align();
      auto startPosition = manager.position();
      transferPass(counter, counter + region.count, crossing);
      manager.align();
      auto endPosition = manager.position();
      auto size = endPosition - startPosition;
      assert(size % typeSize == 0);

      // NOLINTNEXTLINE
      auto* startPtr = reinterpret_cast<void*>(startPosition);

      remoteClusters.emplace_back(startPtr, size / typeSize, datatype, region.rank, region.tag);

      counter += region.count;
    }

    assert(counter == layer.size());
  }

  layer.setEntrySize<LTS::Buffers>(manager.size());
  if (!ghost) {
    layer.setEntrySize<LTS::ZeroFaceTrace>(kernels::Solver::FaceTraceSize * sizeof(real));
    layer.setEntrySize<LTS::TraceSideOperators>(2 * kernels::TraceSideOperatorSize * sizeof(real));
  }

  return remoteClusters;
}

void setupBuckets(LTS::Layer& layer, std::vector<solver::RemoteCluster>& comm) {
  if (layer.getIdentifier().halo != HaloType::Ghost) {
    std::fill_n(layer.var<LTS::ZeroFaceTrace>(), kernels::Solver::FaceTraceSize, real{0});
    kernels::writeTraceSideOperators(layer.var<LTS::TraceSideOperators>());
  }

  auto* buffers = layer.var<LTS::Buffers>();
  auto* buffersDevice = layer.var<LTS::Buffers>(AllocationPlace::Device);
  auto* faceTraces = layer.var<LTS::FaceTraces>();
  auto* faceBuffers = layer.var<LTS::FaceBuffers>();
  auto* faceTracesDevice = layer.var<LTS::FaceTracesDevice>();
  auto* faceBuffersDevice = layer.var<LTS::FaceBuffersDevice>();
  const auto* cellInformation = layer.var<LTS::CellInformation>();

#pragma omp parallel for schedule(static)
  for (std::size_t cell = 0; cell < layer.size(); ++cell) {
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      const auto type = cellInformation[cell].ltsSetup.faceBuffer(face);
      initBucketItem(faceTraces[cell][face], buffers, kernels::Solver::FaceTraceSize, true);
      if (type.has_value()) {
        initBucketItem(faceBuffers[cell][face], buffers, bufferSize(type.value()), true);
      }
      if constexpr (isDeviceOn()) {
        initBucketItem(
            faceTracesDevice[cell][face], buffersDevice, kernels::Solver::FaceTraceSize, false);
        if (type.has_value()) {
          initBucketItem(
              faceBuffersDevice[cell][face], buffersDevice, bufferSize(type.value()), false);
        }
      }
      assert(faceTraces[cell][face] != nullptr || layer.getIdentifier().halo == HaloType::Ghost);
    }
  }

#ifndef NDEBUG
  // every pointer we handed out has to have been resolved against the bucket base by now
  const auto bucketBase = reinterpret_cast<std::uintptr_t>(buffers);
  const auto bucketSize = layer.getEntrySize<LTS::Buffers>();
  const auto inside = [&](const real* pointer, std::size_t size) {
    const auto address = reinterpret_cast<std::uintptr_t>(pointer);
    return address == 0 ||
           (address >= bucketBase && address + size * sizeof(real) <= bucketBase + bucketSize);
  };
  for (std::size_t cell = 0; cell < layer.size(); ++cell) {
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      const auto type = cellInformation[cell].ltsSetup.faceBuffer(face);
      assert(inside(faceTraces[cell][face], kernels::Solver::FaceTraceSize));
      assert(!type.has_value() || inside(faceBuffers[cell][face], bufferSize(type.value())));
    }
  }
#endif

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
      faceNeighbors[cell][face] = nullptr;
      faceNeighborsDevice[cell][face] = nullptr;
      if (getBCType(cellInformation[cell].faceTypes[face]) != BCType::External) {
        if (faceNeighbor == StoragePosition::NullPosition) {
          logError() << "A face that needs the Neighbor kernel has no face neighbor.";
        } else {
          // what the neighbor writes on the same face, which it knows by its own number
          const auto neighborFace = cellInformation[cell].faceRelations[face][0];
          faceNeighbors[cell][face] = storage.lookup<LTS::FaceBuffers>(faceNeighbor)[neighborFace];
          if constexpr (isDeviceOn()) {
            faceNeighborsDevice[cell][face] =
                storage.lookup<LTS::FaceBuffersDevice>(faceNeighbor)[neighborFace];
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
