// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Stream.h"

namespace {
#ifdef ACL_DEVICE
device::DeviceInstance& dev() { return device::DeviceInstance::getInstance(); }
#endif
} // namespace

namespace seissol::parallel::runtime {

ManagedStream::ManagedStream() {
#ifdef ACL_DEVICE
  streamPtr = dev().api->createStream();
#endif
}

ManagedStream::~ManagedStream() {
#ifdef ACL_DEVICE
  dev().api->destroyGenericStream(streamPtr);
#endif
}

ManagedEvent::ManagedEvent() {
#ifdef ACL_DEVICE
  eventPtr = dev().api->createEvent();
#endif
}

ManagedEvent::~ManagedEvent() {
#ifdef ACL_DEVICE
  dev().api->destroyEvent(eventPtr);
#endif
}

} // namespace seissol::parallel::runtime
