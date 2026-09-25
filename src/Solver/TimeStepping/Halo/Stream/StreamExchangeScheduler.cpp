// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "StreamExchangeScheduler.h"

#include "Solver/TimeStepping/Halo/Stream/ExchangeScheduler.h"

#include <algorithm>
#include <cstddef>
#include <utils/logger.h>
#include <vector>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::solver {

std::size_t StreamExchangeScheduler::slot(std::size_t from, std::size_t to) const {
  return launchOrder() == LaunchOrder::Global ? 0 : from * clusterCount_ + to;
}

std::vector<std::size_t> StreamExchangeScheduler::usedSlots() const {
  std::vector<std::size_t> slots;
  for (std::size_t from = 0; from < clusterCount_; ++from) {
    for (std::size_t to = from == 0 ? 0 : from - 1; to < std::min(from + 2, clusterCount_); ++to) {
      const auto current = slot(from, to);
      if (std::find(slots.begin(), slots.end(), current) == slots.end()) {
        slots.push_back(current);
      }
    }
  }
  return slots;
}

#ifdef ACL_DEVICE

namespace {
device::DeviceInstance& deviceInstance() { return device::DeviceInstance::getInstance(); }
} // namespace

StreamExchangeScheduler::StreamExchangeScheduler(std::size_t clusterCount, LaunchOrder order)
    : ExchangeScheduler(clusterCount, order), clusterCount_(clusterCount),
      streams_(order == LaunchOrder::Global ? 1 : clusterCount * clusterCount, nullptr) {
  for (const auto current : usedSlots()) {
    streams_[current] = deviceInstance().api->createStream();
  }
}

StreamExchangeScheduler::~StreamExchangeScheduler() {
  synchronize();
  for (const auto& [ticket, event] : pendingEvents_) {
    deviceInstance().api->destroyEvent(event);
  }
  for (auto* event : launchedEvents_) {
    deviceInstance().api->destroyEvent(event);
  }
  for (auto* stream : streams_) {
    if (stream != nullptr) {
      deviceInstance().api->destroyGenericStream(stream);
    }
  }
}

void StreamExchangeScheduler::synchronize() {
  for (auto* stream : streams_) {
    if (stream != nullptr) {
      deviceInstance().api->syncStreamWithHost(stream);
    }
  }
}

ExchangeScheduler::Ticket StreamExchangeScheduler::launch(std::size_t from,
                                                          std::size_t to,
                                                          std::size_t exchange,
                                                          const ScheduledTransport* sender,
                                                          const ScheduledTransport* receiver,
                                                          const std::vector<void*>& after) {
  const auto current = slot(from, to);
  auto* stream = streams_[current];
  if (stream == nullptr) {
    logError() << "There is no stream for the halo exchange from cluster" << from << "to cluster"
               << to;
  }
  for (auto* event : after) {
    deviceInstance().api->syncStreamWithEvent(stream, event);
  }

  enqueueGroup(current, from, to, exchange, sender, receiver);

  auto* event = deviceInstance().api->createEvent();
  deviceInstance().api->recordEventOnStream(event, stream);
  const auto ticket = nextTicket_++;
  if (streamOrdered()) {
    // the dependent work waits for the event on the device; it stays until the device has completed
    launchedEvents_.push_back(event);
    latestEvent_ = event;
  } else {
    pendingEvents_[ticket] = event;
  }
  return ticket;
}

bool StreamExchangeScheduler::completed(Ticket ticket) {
  const auto pending = pendingEvents_.find(ticket);
  if (pending == pendingEvents_.end()) {
    // tickets are handed out in increasing order; only completed ones are forgotten
    return ticket < nextTicket_;
  }
  if (deviceInstance().api->isEventCompleted(pending->second)) {
    deviceInstance().api->destroyEvent(pending->second);
    pendingEvents_.erase(pending);
    return true;
  }
  return false;
}

std::vector<void*> StreamExchangeScheduler::streams() const {
  std::vector<void*> result;
  for (auto* stream : streams_) {
    if (stream != nullptr) {
      result.push_back(stream);
    }
  }
  return result;
}

void StreamExchangeScheduler::releaseEvents() {
  bool ownsLatest = false;
  for (auto* event : launchedEvents_) {
    if (event != latestEvent_) {
      deviceInstance().api->destroyEvent(event);
    } else {
      ownsLatest = true;
    }
  }
  launchedEvents_.clear();
  if (ownsLatest) {
    launchedEvents_.push_back(latestEvent_);
  }
}

#else

StreamExchangeScheduler::StreamExchangeScheduler(std::size_t clusterCount, LaunchOrder order)
    : ExchangeScheduler(clusterCount, order), clusterCount_(clusterCount) {
  logError() << "Exchanging the halo data on device streams needs a device build.";
}

StreamExchangeScheduler::~StreamExchangeScheduler() = default;

void StreamExchangeScheduler::synchronize() {}

ExchangeScheduler::Ticket StreamExchangeScheduler::launch(std::size_t /*from*/,
                                                          std::size_t /*to*/,
                                                          std::size_t /*exchange*/,
                                                          const ScheduledTransport* /*sender*/,
                                                          const ScheduledTransport* /*receiver*/,
                                                          const std::vector<void*>& /*after*/) {
  return 0;
}

bool StreamExchangeScheduler::completed(Ticket /*ticket*/) { return true; }

std::vector<void*> StreamExchangeScheduler::streams() const { return {}; }

void StreamExchangeScheduler::releaseEvents() {}

#endif

} // namespace seissol::solver
