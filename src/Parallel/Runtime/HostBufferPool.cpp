// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "HostBufferPool.h"

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <mutex>
#include <utils/logger.h>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::parallel::runtime {

namespace {
#ifndef ACL_DEVICE
constexpr std::size_t PageSize = 4096;
#endif
} // namespace

HostBufferPool::HostBufferPool(std::size_t bytes) : bytes_(std::max(bytes, std::size_t{1})) {}

HostBufferPool::~HostBufferPool() {
  if (free_.size() != all_.size()) {
    logWarning() << "A pool of host buffers goes while" << all_.size() - free_.size()
                 << "of its buffers are still in use.";
  }
  for (const auto& buffer : all_) {
#ifdef ACL_DEVICE
    device::DeviceInstance::instance().api().freePinnedMem(buffer.host);
#else
    std::free(buffer.host);
#endif
  }
}

HostBufferPool::Buffer HostBufferPool::take() {
  {
    const std::scoped_lock lock(mutex_);
    if (!free_.empty()) {
      const auto buffer = free_.back();
      free_.pop_back();
      return buffer;
    }
  }
  // (only the thread that takes the buffers changes all of them)
  Buffer buffer;
#ifdef ACL_DEVICE
  auto& api = device::DeviceInstance::instance().api();
  buffer.host = api.allocPinnedMem(bytes_);
  buffer.device = api.devicePointer(buffer.host);
#else
  // aligned like pinned memory, for what the buffers hold (e.g. DOFs read by vector loads)
  buffer.host = std::aligned_alloc(PageSize, (bytes_ + PageSize - 1) / PageSize * PageSize);
  buffer.device = buffer.host;
#endif
  all_.push_back(buffer);
  return buffer;
}

void HostBufferPool::give(const Buffer& buffer) {
  const std::scoped_lock lock(mutex_);
  free_.push_back(buffer);
}

} // namespace seissol::parallel::runtime
