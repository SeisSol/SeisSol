// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Stream.h"

#include <utility>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace {
#ifdef ACL_DEVICE
device::DeviceInstance& dev() { return device::DeviceInstance::getInstance(); }
#endif
} // namespace

namespace seissol::parallel::runtime {

ManagedStream::ManagedStream() {
#ifdef ACL_DEVICE
  streamPtr_ = dev().api->createStream();
#endif
}

auto ManagedStream::operator=(ManagedStream&& old) noexcept -> ManagedStream& {
  if (this != &old) {
#ifdef ACL_DEVICE
    if (streamPtr_ != nullptr) {
      dev().api->destroyGenericStream(streamPtr_);
    }
#endif
    streamPtr_ = std::exchange(old.streamPtr_, nullptr);
  }
  return *this;
}

ManagedStream::~ManagedStream() {
#ifdef ACL_DEVICE
  if (streamPtr_ != nullptr) {
    dev().api->destroyGenericStream(streamPtr_);
  }
#endif
}

ManagedEvent::ManagedEvent() {
#ifdef ACL_DEVICE
  eventPtr_ = dev().api->createEvent();
#endif
}

auto ManagedEvent::operator=(ManagedEvent&& old) noexcept -> ManagedEvent& {
  if (this != &old) {
#ifdef ACL_DEVICE
    if (eventPtr_ != nullptr) {
      dev().api->destroyEvent(eventPtr_);
    }
#endif
    eventPtr_ = std::exchange(old.eventPtr_, nullptr);
  }
  return *this;
}

ManagedEvent::~ManagedEvent() {
#ifdef ACL_DEVICE
  if (eventPtr_ != nullptr) {
    dev().api->destroyEvent(eventPtr_);
  }
#endif
}

} // namespace seissol::parallel::runtime
