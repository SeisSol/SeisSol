// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "EventPool.h"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <memory>
#include <utility>
#include <utils/logger.h>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::parallel::runtime {

namespace {
void* createDeviceEvent() {
#ifdef ACL_DEVICE
  return device::DeviceInstance::instance().api().createEvent();
#else
  return nullptr;
#endif
}

void destroyDeviceEvent([[maybe_unused]] void* event) {
#ifdef ACL_DEVICE
  device::DeviceInstance::instance().api().destroyEvent(event);
#endif
}
} // namespace

EventPool::EventPool() : EventPool(createDeviceEvent, destroyDeviceEvent) {}

EventPool::EventPool(Create create, Destroy destroy, std::size_t initialSize)
    : create_(std::move(create)), destroy_(std::move(destroy)),
      referenced_(std::make_unique<std::atomic<std::size_t>>(0)) {
  grow(initialSize);
}

EventPool::~EventPool() { dispose(); }

void EventPool::grow(std::size_t count) {
  slots_.reserve(slots_.size() + count);
  for (std::size_t i = 0; i < count; ++i) {
    auto slot = std::make_unique<internal::EventSlot>();
    slot->event = create_();
    slot->referenced = referenced_.get();
    slots_.push_back(std::move(slot));
  }
}

EventRef EventPool::next() {
  if (disposed_) {
    logError() << "An event pool was asked for an event after it had been disposed.";
  }

  for (std::size_t probe = 0; probe < slots_.size(); ++probe) {
    auto* slot = slots_[position_].get();
    position_ = (position_ + 1) % slots_.size();

    if (slot->refs.load(std::memory_order_acquire) == 0) {
      EventRef event(slot);
      peakReferenced_ = std::max(peakReferenced_, referenced());
      return event;
    }
  }

  // every event is still spoken for
  const auto oldSize = slots_.size();
  if (oldSize >= MaxSize) {
    logError() << "Ran out of device events (" << oldSize
               << " referenced at once). This points at event references being kept alive "
                  "past the wait they belong to.";
  }
  grow(std::min(Growth, MaxSize - oldSize));
  position_ = oldSize;
  return next();
}

std::size_t EventPool::dispose() {
  if (disposed_) {
    return 0;
  }
  disposed_ = true;

  std::size_t leaked = 0;
  for (auto& slot : slots_) {
    if (slot->refs.load(std::memory_order_acquire) == 0) {
      if (slot->event != nullptr) {
        destroy_(slot->event);
      }
    } else {
      // Somebody may still wait for the event, and will drop the reference to it later. So the
      // event stays as it is, and so does its slot.
      ++leaked;
      [[maybe_unused]] auto* kept = slot.release();
    }
  }
  slots_.clear();

  if (leaked > 0) {
    // the kept slots still count their references here
    [[maybe_unused]] auto* kept = referenced_.release();
    logWarning(true) << "An event pool was disposed while" << leaked
                     << "of its events were still referenced; they are left alone.";
  }
  return leaked;
}

} // namespace seissol::parallel::runtime
