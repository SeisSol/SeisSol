#include "DataCollector.h"

#ifdef ACL_DEVICE
#include <device.h>
#endif

namespace seissol::parallel {

DataCollectorUntyped::DataCollectorUntyped(const std::vector<void*>& indexDataHost,
                                           size_t elemSize,
                                           bool hostAccessible)
    : indexCount(indexDataHost.size()), indexDataHost(indexDataHost), elemSize(elemSize),
      hostAccessible(hostAccessible) {
  // in case we want to use this class in a host-only scenario
  if constexpr (!isDeviceOn()) {
    this->hostAccessible = true;
  }
  if (!this->hostAccessible && indexCount > 0) {
    indexDataDevice = memory::allocTyped<void*>(indexCount, 1, memory::DeviceGlobalMemory);
    copiedData = memory::allocate(elemSize * indexCount, 1, memory::PinnedMemory);

    copiedDataDevice = memory::hostToDevicePointerTyped(copiedData, memory::PinnedMemory);
    memory::memcopyTyped(indexDataDevice,
                         indexDataHost.data(),
                         indexCount,
                         memory::DeviceGlobalMemory,
                         memory::Standard);
  }
}

DataCollectorUntyped::~DataCollectorUntyped() {
  if (!hostAccessible && indexCount > 0) {
    memory::free(static_cast<void*>(indexDataDevice), memory::DeviceGlobalMemory);
    memory::free(copiedData, memory::PinnedMemory);
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
