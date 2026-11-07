// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "HostSnapshot.h"

#include "Memory/Tree/Layer.h"

#include <vector>

#ifdef ACL_DEVICE
#include <Device/device.h>
#include <cassert>
#endif

namespace seissol::initializer {

#ifdef ACL_DEVICE

HostSnapshot::~HostSnapshot() {
  if (stream_ != nullptr) {
    device::DeviceInstance::instance().api().destroyGenericStream(stream_);
  }
}

void HostSnapshot::take(const std::vector<const DualMemoryContainer*>& containers, void* stream) {
  assert(pending_.empty());
  auto& api = device::DeviceInstance::instance().api();
  for (const auto* container : containers) {
    if (!container->copiesToHost()) {
      continue;
    }
    if (container->snapshot != nullptr) {
      container->takeSnapshot(stream);
      if (container->snapshotApartFromHost()) {
        pending_.push_back(container);
      }
    } else {
      // without a snapshot, the data gets to the host before the device goes on
      container->synchronizeTo(AllocationPlace::Host, stream);
    }
  }
  if (!pending_.empty() && stream_ == nullptr) {
    stream_ = api.createStream();
  }
  api.syncStreamWithHost(stream);
}

void HostSnapshot::copyToHost() {
  if (pending_.empty()) {
    return;
  }
  auto& api = device::DeviceInstance::instance().api();
  // a thread that has not used the device yet selects the one of the process here
  api.getDeviceId();
  for (const auto* container : pending_) {
    container->snapshotToHost(stream_);
  }
  api.syncStreamWithHost(stream_);
  pending_.clear();
}

#else

void HostSnapshot::take(const std::vector<const DualMemoryContainer*>& containers,
                        void* /*stream*/) {
  for (const auto* container : containers) {
    container->synchronizeTo(AllocationPlace::Host, nullptr);
  }
}

void HostSnapshot::copyToHost() {}

#endif

} // namespace seissol::initializer
