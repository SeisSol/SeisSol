// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Stream.h"

#include "Initializer/BasicTypedefs.h"
#include "Memory/Tree/Colormap.h"

#include <cstddef>
#include <optional>
#include <utility>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace {
#ifdef ACL_DEVICE
device::DeviceInstance& dev() { return device::DeviceInstance::instance(); }
#endif

/// the HSA queues of useHsaQueues(), the only one to change them
seissol::parallel::runtime::HsaQueues& hsaQueueSetup() {
  static seissol::parallel::runtime::HsaQueues queues;
  return queues;
}
} // namespace

namespace seissol::parallel::runtime {

#ifdef ACL_DEVICE
ManagedStream::ManagedStream() : streamPtr_(dev().api().createStream()) {}

ManagedStream::ManagedStream(double priority) : streamPtr_(dev().api().createStream(priority)) {}

ManagedStream::ManagedStream(OnHsaQueue queue)
    : streamPtr_(dev().api().createStreamOnHardwareQueue(queue.queue, queue.priority)) {}
#else
ManagedStream::ManagedStream() = default;

ManagedStream::ManagedStream(double /*priority*/) {}

ManagedStream::ManagedStream(OnHsaQueue /*queue*/) {}
#endif

const HsaQueues& hsaQueues() { return hsaQueueSetup(); }

void useHsaQueues(const HsaQueues& queues) { hsaQueueSetup() = queues; }

int hsaQueueOf(const initializer::LayerIdentifier& identifier, bool faces) {
  const auto& queues = hsaQueues();
  if (queues.count <= 0) {
    return SharedHsaQueue;
  }
  const bool copy = identifier.halo == HaloType::Copy;
  if (faces && queues.faceCount > 0) {
    const auto layer = identifier.lts * 2 + (copy ? 1 : 0);
    return FirstFaceHsaQueue - static_cast<int>(layer % static_cast<std::size_t>(queues.faceCount));
  }
  if (queues.layout == HsaQueueLayout::Halo) {
    return queues.count > 1 && !copy ? 2 : 1;
  }
  const auto layer = identifier.lts * 2 + (copy ? 1 : 0);
  return 1 + static_cast<int>(layer % static_cast<std::size_t>(queues.count));
}

std::optional<double> hsaQueuePriorityOf(const initializer::LayerIdentifier& identifier,
                                         bool faces) {
  const auto& queues = hsaQueues();
  if (queues.count <= 0) {
    return std::nullopt;
  }
  if (faces && queues.faceCount > 0) {
    return queues.priorities.faces ? std::optional<double>(1.0) : std::nullopt;
  }
  // (the queues of the layers hold either copy or interior layers: with an even number of them in
  // the layers layout, or with the halo layout)
  const bool byLayer = queues.layout == HsaQueueLayout::Halo
                           ? queues.count > 1
                           : queues.priorities.copy && queues.count % 2 == 0;
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
