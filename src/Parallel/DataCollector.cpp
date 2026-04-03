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
    : hostAccessible(hostAccessible), indexDataHost(indexDataHost),
      indexCount(indexDataHost.size()), elemSize(elemSize) {
  // in case we want to use this class in a host-only scenario
  if constexpr (!isDeviceOn()) {
    this->hostAccessible = true;
  }
  if (!this->hostAccessible && indexCount > 0) {
    indexDataDevice = memory::allocTyped<void*>(indexCount, 1, memory::Memkind::DeviceGlobalMemory);
    copiedData = memory::allocate(elemSize * indexCount, 1, memory::Memkind::PinnedMemory);

    copiedDataDevice = memory::hostToDevicePointerTyped(copiedData, memory::Memkind::PinnedMemory);
    memory::memcopyTyped(indexDataDevice,
                         indexDataHost.data(),
                         indexCount,
                         memory::Memkind::DeviceGlobalMemory,
                         memory::Memkind::Standard);
  }
}

DataCollectorUntyped::~DataCollectorUntyped() {
  if (!hostAccessible && indexCount > 0) {
    memory::free(static_cast<void*>(indexDataDevice), memory::Memkind::DeviceGlobalMemory);
    memory::free(copiedData, memory::Memkind::PinnedMemory);
  }
}

// NOLINTNEXTLINE (silence errors due to ifdef ACL_DEVICE)
void DataCollectorUntyped::gatherToHost(void* stream) {
  if (!hostAccessible && indexCount > 0) {
#ifdef ACL_DEVICE
    device::DeviceInstance::getInstance().algorithms.copyScatterToUniformI(
        const_cast<const void**>(indexDataDevice),
        copiedDataDevice,
        elemSize,
        elemSize,
        indexCount,
        stream);
#endif
  }
}

// NOLINTNEXTLINE (silence errors due to ifdef ACL_DEVICE)
void DataCollectorUntyped::scatterFromHost(void* stream) {
  if (!hostAccessible && indexCount > 0) {
#ifdef ACL_DEVICE
    device::DeviceInstance::getInstance().algorithms.copyUniformToScatterI(
        const_cast<const void*>(copiedDataDevice),
        indexDataDevice,
        elemSize,
        elemSize,
        indexCount,
        stream);
#endif
  }
}

} // namespace seissol::parallel
