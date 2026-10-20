// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "HostSnapshot.h"

#include "Memory/HostData.h"

#include <vector>

#ifdef ACL_DEVICE
#include "Monitoring/Unit.h"

#include <Device/device.h>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <utils/logger.h>
#endif

namespace seissol::initializer {

#ifdef ACL_DEVICE

namespace {
// the copies start at offsets the device aligns its allocations to
constexpr std::size_t Alignment = 256;

std::size_t aligned(std::size_t size) { return (size + Alignment - 1) / Alignment * Alignment; }
} // namespace

HostSnapshot::~HostSnapshot() {
  auto& api = device::DeviceInstance::instance().api();
  if (buffer_ != nullptr) {
    api.freeGlobMem(buffer_);
  }
  if (stream_ != nullptr) {
    api.destroyGenericStream(stream_);
  }
}

bool HostSnapshot::take(const std::vector<HostCopy>& copies, void* stream) {
  assert(pending_.empty());
  auto& api = device::DeviceInstance::instance().api();
  std::size_t size = 0;
  for (const auto& copy : copies) {
    size += aligned(copy.size);
    if (!copy.pinned) {
      if (!warnedPinned_) {
        logWarning() << "The host memory of data for the outputs is not pinned; the data gets"
                     << "copied to the host before the next epoch.";
        warnedPinned_ = true;
      }
      return false;
    }
  }

  if (size > capacity_) {
    // what SeisSol has allocated besides
    const auto others = api.getCurrentlyOccupiedMem() - capacity_;
    const auto total = api.getMaxAvailableMem();
    constexpr std::size_t MinimumLeft = 2UL << 30U;
    if (others + size + std::max(total / 10, MinimumLeft) > total) {
      if (!warnedSize_) {
        logWarning() << "A snapshot of the data for the outputs would take"
                     << UnitByte.formatPrefix(size).c_str()
                     << "on the device, more than fits; the data gets copied to the host before"
                     << "the next epoch.";
        warnedSize_ = true;
      }
      return false;
    }
    if (buffer_ != nullptr) {
      api.freeGlobMem(buffer_);
    }
    buffer_ = api.allocGlobMem(size);
    capacity_ = size;
    logInfo() << "The snapshots of the data for the outputs take"
              << UnitByte.formatPrefix(size).c_str() << "on the device.";
  }
  if (stream_ == nullptr) {
    stream_ = api.createStream();
  }

  std::size_t offset = 0;
  for (const auto& copy : copies) {
    auto* target = static_cast<char*>(buffer_) + offset;
    api.copyBetweenAsync(target, copy.device, copy.size, stream);
    pending_.push_back(HostCopy{copy.host, target, copy.size, copy.pinned});
    offset += aligned(copy.size);
  }
  api.syncStreamWithHost(stream);
  return true;
}

void HostSnapshot::copyToHost() {
  if (pending_.empty()) {
    return;
  }
  auto& api = device::DeviceInstance::instance().api();
  // a thread that has not used the device yet selects the one of the process here
  api.getDeviceId();
  for (const auto& copy : pending_) {
    api.copyFromAsync(copy.host, copy.device, copy.size, stream_);
  }
  api.syncStreamWithHost(stream_);
  pending_.clear();
}

#else

bool HostSnapshot::take(const std::vector<HostCopy>& /*copies*/, void* /*stream*/) { return false; }

void HostSnapshot::copyToHost() {}

#endif

} // namespace seissol::initializer
