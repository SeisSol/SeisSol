// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Stream.h"

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

ManagedStream::~ManagedStream() {
#ifdef ACL_DEVICE
  dev().api->destroyGenericStream(streamPtr_);
#endif
}

ManagedEvent::ManagedEvent() {
#ifdef ACL_DEVICE
  eventPtr_ = dev().api->createEvent();
#endif
}

ManagedEvent::~ManagedEvent() {
#ifdef ACL_DEVICE
  dev().api->destroyEvent(eventPtr_);
#endif
}

} // namespace seissol::parallel::runtime
