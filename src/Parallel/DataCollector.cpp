// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "DataCollector.h"

#include "Kernels/Common.h"
#include "Memory/MemoryAllocator.h"

#include <cstddef>
#include <vector>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::parallel {

DataCollectorUntyped::DataCollectorUntyped(const std::vector<void*>& indexDataHost,
                                           size_t elemSize,
                                           bool hostAccessible)
    : hostAccessible_(hostAccessible), indexDataHost_(indexDataHost),
      indexCount_(indexDataHost.size()), elemSize_(elemSize) {
  // in case we want to use this class in a host-only scenario
  if constexpr (!isDeviceOn()) {
    this->hostAccessible_ = true;
  }
  if (!this->hostAccessible_ && indexCount_ > 0) {
    indexDataDevice_ =
        memory::allocTyped<void*>(indexCount_, 1, memory::Memkind::DeviceGlobalMemory);
    copiedData_ = memory::allocate(elemSize * indexCount_, 1, memory::Memkind::PinnedMemory);

    copiedDataDevice_ =
        memory::hostToDevicePointerTyped(copiedData_, memory::Memkind::PinnedMemory);
    memory::memcopyTyped(indexDataDevice_,
                         indexDataHost.data(),
                         indexCount_,
                         memory::Memkind::DeviceGlobalMemory,
                         memory::Memkind::Standard);
  }
}

DataCollectorUntyped::~DataCollectorUntyped() {
  if (!hostAccessible_ && indexCount_ > 0) {
    memory::free(static_cast<void*>(indexDataDevice_), memory::Memkind::DeviceGlobalMemory);
    memory::free(copiedData_, memory::Memkind::PinnedMemory);
  }
}

// NOLINTNEXTLINE (silence errors due to ifdef ACL_DEVICE)
void DataCollectorUntyped::gatherToHost(void* stream) {
  if (!hostAccessible_ && indexCount_ > 0) {
#ifdef ACL_DEVICE
    device::DeviceInstance::getInstance().algorithms.copyScatterToUniformI(
        const_cast<const void**>(indexDataDevice_),
        copiedDataDevice_,
        elemSize_,
        elemSize_,
        indexCount_,
        stream);
#endif
  }
}

// NOLINTNEXTLINE (silence errors due to ifdef ACL_DEVICE)
void DataCollectorUntyped::scatterFromHost(void* stream) {
  if (!hostAccessible_ && indexCount_ > 0) {
#ifdef ACL_DEVICE
    device::DeviceInstance::getInstance().algorithms.copyUniformToScatterI(
        const_cast<const void*>(copiedDataDevice_),
        indexDataDevice_,
        elemSize_,
        elemSize_,
        indexCount_,
        stream);
#endif
  }
}

} // namespace seissol::parallel
