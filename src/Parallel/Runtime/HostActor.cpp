// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "HostActor.h"

#include "Parallel/Runtime/Event.h"
#include "Parallel/Runtime/HostRuntime.h"
#include "Parallel/Runtime/Stream.h"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <functional>
#include <thread>
#include <utility>
#include <utils/logger.h>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::parallel::runtime {

HostActor::Queue::Queue() : Queue(HostActor::instance()) {}

HostActor::Queue::Queue(HostActor& actor) : actor_(actor) { actor_.queues_.push_back(this); }

HostActor::Queue::~Queue() {
  // whatever the work reads or writes may go right after this; it has to be done before
  drain();
  auto& queues = actor_.queues_;
  queues.erase(std::remove(queues.begin(), queues.end(), this), queues.end());
}

void HostActor::Queue::post(Event after, std::function<void()> work) {
  if (after.graph() != nullptr) {
    logError() << "Work on the host cannot wait for a node of a graph that is being built.";
  }
  items_.push_back(Item{std::move(after), std::move(work)});
  ++posted_;
}

bool HostActor::Queue::progress() {
  bool handedOver = false;
  while (!items_.empty() && actor_.ready(items_.front().after)) {
    // out of the queue first: the work may post more
    auto task = [this, work = std::move(items_.front().work)]() {
      std::invoke(work);
      done_.fetch_add(1, std::memory_order_release);
    };
    items_.pop_front();
    if (actor_.executor_) {
      actor_.executor_(std::move(task));
    } else {
      // the work changes only what belongs to the queue
      runtime_.action(HostAccess{}, std::move(task));
    }
    handedOver = true;
  }
  return handedOver;
}

void HostActor::Queue::drain() {
  if (!items_.empty() && recordingOuterGraph()) {
    logError() << "Work on the host cannot be waited for while a graph is being recorded.";
  }
  while (!items_.empty()) {
    if (!progress()) {
      std::this_thread::yield();
    }
  }
  // the work handed over to tasks or to the executor
  while (!empty()) {
    HostRuntime::yield();
    std::this_thread::yield();
  }
}

HostActor::HostActor()
    : HostActor([]([[maybe_unused]] const Event& event) {
#ifdef ACL_DEVICE
        return device::DeviceInstance::instance().api().isEventCompleted(event.get());
#else
        return true;
#endif
      }) {
}

HostActor::HostActor(Completed completed) : completed_(std::move(completed)) {}

HostActor& HostActor::instance() {
  static HostActor actor;
  return actor;
}

void HostActor::setExecutor(Executor executor) {
  if (pending() > 0) {
    logError() << "The host actor cannot change where its work runs while some is pending.";
  }
  executor_ = std::move(executor);
}

bool HostActor::progress() {
  // (no event may be queried while a graph is being recorded)
  if (recordingOuterGraph()) {
    return false;
  }
  bool handedOver = false;
  // by index: the work may set up queues
  for (std::size_t i = 0; i < queues_.size(); ++i) {
    handedOver = queues_[i]->progress() || handedOver;
  }
  return handedOver;
}

void HostActor::drain() {
  for (std::size_t i = 0; i < queues_.size(); ++i) {
    queues_[i]->drain();
  }
}

std::size_t HostActor::pending() const {
  std::size_t pending = 0;
  for (const auto* queue : queues_) {
    pending += queue->pending();
  }
  return pending;
}

bool HostActor::ready(const Event& event) const { return !event || completed_(event); }

} // namespace seissol::parallel::runtime
