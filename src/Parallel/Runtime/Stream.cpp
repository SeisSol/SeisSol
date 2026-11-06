// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Stream.h"

#include <optional>
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

ManagedStream::ManagedStream(OnHsaQueue queue)
    : streamPtr_(dev().api().createStreamOnHsaQueue(queue.queue, queue.priority)) {}
#else
ManagedStream::ManagedStream() = default;

ManagedStream::ManagedStream(double /*priority*/) {}

ManagedStream::ManagedStream(OnHsaQueue /*queue*/) {}
#endif

int hsaQueueOf(const initializer::LayerIdentifier& identifier, bool faces) {
  if (hsaQueueCount() <= 0) {
    return SharedHsaQueue;
  }
  const bool copy = identifier.halo == HaloType::Copy;
  if (faces && hsaFaceQueueCount() > 0) {
    const auto layer = identifier.lts * 2 + (copy ? 1 : 0);
    return FirstFaceHsaQueue -
           static_cast<int>(layer % static_cast<std::size_t>(hsaFaceQueueCount()));
  }
  if (hsaQueueLayout() == HsaQueueLayout::Halo) {
    return hsaQueueCount() > 1 && !copy ? 2 : 1;
  }
  const auto layer = identifier.lts * 2 + (copy ? 1 : 0);
  return 1 + static_cast<int>(layer % static_cast<std::size_t>(hsaQueueCount()));
}

std::optional<double> hsaQueuePriorityOf(const initializer::LayerIdentifier& identifier,
                                         bool faces) {
  if (hsaQueueCount() <= 0) {
    if (streamPriorities().copy) {
      return identifier.halo == HaloType::Copy ? 1.0 : 0.0;
    }
    return std::nullopt;
  }
  if (faces && hsaFaceQueueCount() > 0) {
    return hsaQueuePriorities().faces ? std::optional<double>(1.0) : std::nullopt;
  }
  // (the queues of the layers hold either copy or interior layers: with an even number of them in
  // the layers layout, or with the halo layout)
  const bool byLayer = hsaQueueLayout() == HsaQueueLayout::Halo
                           ? hsaQueueCount() > 1
                           : hsaQueuePriorities().copy && hsaQueueCount() % 2 == 0;
  if (!byLayer) {
    return std::nullopt;
  }
  return identifier.halo == HaloType::Copy ? 1.0 : 0.0;
}

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
