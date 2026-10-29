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

ManagedStream::ManagedStream(OnOwnQueue queue)
    : streamPtr_(dev().api().createStreamOnOwnQueue(queue.queue, queue.priority)) {}
#else
ManagedStream::ManagedStream() = default;

ManagedStream::ManagedStream(double /*priority*/) {}

ManagedStream::ManagedStream(OnOwnQueue /*queue*/) {}
#endif

int ownQueueOf(const initializer::LayerIdentifier& identifier, bool faces) {
  if (ownQueueCount() <= 0) {
    return SharedOwnQueue;
  }
  const bool copy = identifier.halo == HaloType::Copy;
  if (faces && ownFaceQueueCount() > 0) {
    const auto layer = identifier.lts * 2 + (copy ? 1 : 0);
    return FirstFaceOwnQueue -
           static_cast<int>(layer % static_cast<std::size_t>(ownFaceQueueCount()));
  }
  if (ownQueueLayout() == OwnQueueLayout::Halo) {
    return ownQueueCount() > 1 && !copy ? 2 : 1;
  }
  const auto layer = identifier.lts * 2 + (copy ? 1 : 0);
  return 1 + static_cast<int>(layer % static_cast<std::size_t>(ownQueueCount()));
}

std::optional<double> ownQueuePriorityOf(const initializer::LayerIdentifier& identifier,
                                         bool faces) {
  if (ownQueueCount() <= 0) {
    return std::nullopt;
  }
  if (faces && ownFaceQueueCount() > 0) {
    return ownQueuePriorities().faces ? std::optional<double>(1.0) : std::nullopt;
  }
  // (the queues of the layers hold either copy or interior layers: with an even number of them in
  // the layers layout, or with the halo layout)
  const bool byLayer = ownQueueLayout() == OwnQueueLayout::Halo
                           ? ownQueueCount() > 1
                           : ownQueuePriorities().copy && ownQueueCount() % 2 == 0;
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
