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
device::DeviceInstance& dev() { return device::DeviceInstance::instance(); }
#endif
} // namespace

namespace seissol::parallel::runtime {

#ifdef ACL_DEVICE
ManagedStream::ManagedStream() : streamPtr_(dev().api().createStream()) {}

ManagedStream::ManagedStream(double priority) : streamPtr_(dev().api().createStream(priority)) {}
#else
ManagedStream::ManagedStream() = default;

ManagedStream::ManagedStream(double /*priority*/) {}
#endif

auto ManagedStream::operator=(ManagedStream&& old) noexcept -> ManagedStream& {
  if (this != &old) {
#ifdef ACL_DEVICE
    if (streamPtr_ != nullptr) {
      dev().api().destroyGenericStream(streamPtr_);
    }
#endif
    streamPtr_ = std::exchange(old.streamPtr_, nullptr);
  }
  return *this;
}

ManagedStream::~ManagedStream() {
#ifdef ACL_DEVICE
  if (streamPtr_ != nullptr) {
    dev().api().destroyGenericStream(streamPtr_);
  }
#endif
}

} // namespace seissol::parallel::runtime
