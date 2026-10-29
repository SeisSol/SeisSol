// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "DataCollector.h"

#include "Kernels/Common.h"
#include "Memory/MemoryAllocator.h"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>

#ifdef ACL_DEVICE
#include "Kernels/DeviceAux/GatherAux.h"

#include <Device/device.h>
#endif

namespace seissol::parallel {

namespace {
// where the gathered data goes: memory the host reads directly
memory::Memkind gatherMemkind(bool hostAccessible) {
  return hostAccessible ? memory::Memkind::DeviceUnifiedMemory : memory::Memkind::PinnedMemory;
}
} // namespace

DataCollectorUntyped::DataCollectorUntyped(const std::vector<void*>& indexDataHost,
                                           size_t elemSize,
                                           bool hostAccessible)
    : DataCollectorUntyped(std::vector<std::vector<void*>>{indexDataHost},
                           std::vector<size_t>{elemSize},
                           hostAccessible) {}

DataCollectorUntyped::DataCollectorUntyped(const std::vector<std::vector<void*>>& sets,
                                           const std::vector<size_t>& sizes,
                                           bool hostAccessible)
    : hostAccessible_(hostAccessible) {
  assert(sets.size() == sizes.size());
  for (size_t set = 0; set < sets.size(); ++set) {
    sets_.push_back({indexDataHost_.size(), sets[set].size(), sizes[set], gatheredBytes_});
    indexDataHost_.insert(indexDataHost_.end(), sets[set].begin(), sets[set].end());
    gatheredBytes_ += sets[set].size() * sizes[set];
  }
  indexCount_ = indexDataHost_.size();

  // in case we want to use this class in a host-only scenario
  if constexpr (!isDeviceOn()) {
    direct_ = true;
  }
  if (!direct_ && indexCount_ > 0) {
    indexDataDevice_ =
        memory::allocTyped<void*>(indexCount_, 1, memory::Memkind::DeviceGlobalMemory);
    copiedData_ = memory::allocate(gatheredBytes_, 1, gatherMemkind(hostAccessible_));

    copiedDataDevice_ = memory::hostToDevicePointer(copiedData_, gatherMemkind(hostAccessible_));
    memory::memcopyTyped(indexDataDevice_,
                         indexDataHost_.data(),
                         indexCount_,
                         memory::Memkind::DeviceGlobalMemory,
                         memory::Memkind::Standard);
  }
}

DataCollectorUntyped::~DataCollectorUntyped() {
  if (!direct_ && indexCount_ > 0) {
    memory::freeTyped(indexDataDevice_, memory::Memkind::DeviceGlobalMemory);
    memory::free(copiedData_, gatherMemkind(hostAccessible_));
  }
}

void DataCollectorUntyped::gatherToHost(void* stream) { gatherTo(copiedDataDevice_, stream); }

// NOLINTNEXTLINE (silence errors due to ifdef ACL_DEVICE)
void DataCollectorUntyped::gatherTo(void* destination, void* stream) {
  if (!direct_ && indexCount_ > 0) {
#ifdef ACL_DEVICE
    if (sets_.size() == 1) {
      // a block of threads per block, which suits large blocks (e.g. the DOFs of cells)
      device::DeviceInstance::instance().algorithms().copyScatterToUniformI(
          const_cast<const void**>(indexDataDevice_),
          destination,
          sets_[0].size,
          sets_[0].size,
          indexCount_,
          stream);
    } else {
      kernels::device::aux::gather::copyBlocks(
          sets_.data(), sets_.size(), indexDataDevice_, destination, false, stream);
    }
#endif
  }
}

const void* DataCollectorUntyped::get(size_t set, size_t index, const void* gathered) const {
  const auto& blocks = sets_[set];
  if (direct_) {
    return indexDataHost_[blocks.begin + index];
  }
  return static_cast<const uint8_t*>(gathered) + blocks.offset + index * blocks.size;
}

// NOLINTNEXTLINE (silence errors due to ifdef ACL_DEVICE)
void DataCollectorUntyped::scatterFromHost(void* stream) {
  if (!direct_ && indexCount_ > 0) {
#ifdef ACL_DEVICE
    if (sets_.size() == 1) {
      device::DeviceInstance::instance().algorithms().copyUniformToScatterI(
          const_cast<const void*>(copiedDataDevice_),
          indexDataDevice_,
          sets_[0].size,
          sets_[0].size,
          indexCount_,
          stream);
    } else {
      kernels::device::aux::gather::copyBlocks(
          sets_.data(), sets_.size(), indexDataDevice_, copiedDataDevice_, true, stream);
    }
#endif
  }
}

} // namespace seissol::parallel
