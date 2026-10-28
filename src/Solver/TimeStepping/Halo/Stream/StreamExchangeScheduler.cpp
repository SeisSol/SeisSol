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
#include <numeric>
#include <utils/logger.h>
#include <vector>

#ifdef ACL_DEVICE
#include "Solver/TimeStepping/Actor/ActorState.h"

#include <Device/device.h>
#include <memory>
#include <utility>
#endif

namespace seissol::solver {

std::size_t StreamExchangeScheduler::slot(std::size_t from, std::size_t to) const {
  return launchOrder() == LaunchOrder::Global ? 0 : from * clusterCount_ + to;
}

std::size_t
    StreamExchangeScheduler::slot(std::size_t from, std::size_t to, std::size_t exchange) const {
  if (launchOrder() == LaunchOrder::Global) {
    // the exchanges of a direction take turns on the lanes, and the directions start apart
    return (from * clusterCount_ + to + exchange) % lanes_;
  }
  return slot(from, to);
}

std::size_t StreamExchangeScheduler::slotCount() const {
  return launchOrder() == LaunchOrder::Global ? lanes_ : clusterCount_ * clusterCount_;
}

std::vector<std::size_t> StreamExchangeScheduler::usedSlots() const {
  if (launchOrder() == LaunchOrder::Global) {
    std::vector<std::size_t> lanes(lanes_);
    std::iota(lanes.begin(), lanes.end(), 0);
    return lanes;
  }
  std::vector<std::size_t> slots;
  for (std::size_t from = 0; from < clusterCount_; ++from) {
    for (std::size_t to = 0; to < clusterCount_; ++to) {
      // the halo exchanges go between neighboring time clusters only
      if (timeCluster(to) + 1 < timeCluster(from) || timeCluster(from) + 1 < timeCluster(to)) {
        continue;
      }
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
device::DeviceInstance& deviceInstance() { return device::DeviceInstance::instance(); }
} // namespace

StreamExchangeScheduler::StreamExchangeScheduler(std::size_t clusterCount,
                                                 LaunchOrder order,
                                                 std::size_t configCount,
                                                 std::size_t lanes)
    : ExchangeScheduler(clusterCount, order, configCount), clusterCount_(clusterCount),
      lanes_(std::max<std::size_t>(lanes, 1)) {
  if (lanes_ > 1 && order != LaunchOrder::Global) {
    logError() << "Lanes of the halo exchange need the global order.";
  }
  runtimes_.resize(slotCount());
  for (const auto current : usedSlots()) {
    // lanes on hardware queues of their own (HIP has a pool of them for each priority)
    runtimes_[current] = lanes_ > 1 ? std::make_unique<parallel::runtime::StreamRuntime>(0, 1.0)
                                    : std::make_unique<parallel::runtime::StreamRuntime>();
  }
  if (lanes_ > 1) {
    join_ = std::make_unique<parallel::runtime::StreamRuntime>(0, 1.0);
  }
}

StreamExchangeScheduler::~StreamExchangeScheduler() {
  synchronize();
  // the pools go before the members of the base class, and their events with them
  forgetEvents();
  pendingEvents_.clear();
  for (auto& runtime : runtimes_) {
    if (runtime != nullptr) {
      runtime->dispose();
    }
  }
  if (join_ != nullptr) {
    join_->dispose();
  }
}

void StreamExchangeScheduler::synchronize() {
  for (auto& runtime : runtimes_) {
    if (runtime != nullptr) {
      runtime->wait();
    }
  }
  if (join_ != nullptr) {
    join_->wait();
  }
}

void StreamExchangeScheduler::isolateGroups() {
  for (auto& runtime : runtimes_) {
    if (runtime != nullptr) {
      runtime->isolateNodes();
    }
  }
}

ExchangeScheduler::Ticket StreamExchangeScheduler::launch(std::size_t from,
                                                          std::size_t to,
                                                          std::size_t exchange,
                                                          const ScheduledTransport* sender,
                                                          const ScheduledTransport* receiver,
                                                          const std::vector<ActorEvent>& after) {
  const auto current = slot(from, to, exchange);
  if (runtimes_[current] == nullptr) {
    logError() << "There is no stream for the halo exchange from cluster" << from << "to cluster"
               << to;
  }
  auto& runtime = *runtimes_[current];
  for (const auto& event : after) {
    runtime.eventSync(event);
  }

  // the time from where the group can start, the waits for the remote side included
  const bool timed = timer_ != nullptr;
  const auto interval = timed ? timer_->begin(runtime.stream()) : DeviceTimer::Interval{};
  enqueueGroup(current, from, to, exchange, sender, receiver);
  if (timed) {
    timer_->end(interval, runtime.stream(), deviceTime_);
  }

  auto event = runtime.eventRecord();
  const auto ticket = nextTicket_++;
  if (streamOrdered() && transportEvents() && lanes_ == 1) {
    // on the one stream, the event of a group stands for all groups before it
    for (const auto* transport : {sender, receiver}) {
      if (transport != nullptr) {
        setTransportEvent(*transport, event);
      }
    }
  }
  if (streamOrdered() && join_ != nullptr) {
    // the latest event stands for the groups of all lanes
    join_->eventSync(event);
    event = join_->eventRecord();
  }
  if (streamOrdered()) {
    // the dependent work waits for the event on the device; whoever gets it to wait for keeps it
    setLatestEvent(std::move(event));
  } else {
    pendingEvents_.emplace(ticket, std::move(event));
  }
  return ticket;
}

bool StreamExchangeScheduler::completed(Ticket ticket) {
  const auto pending = pendingEvents_.find(ticket);
  if (pending == pendingEvents_.end()) {
    // tickets are handed out in increasing order; only completed ones are forgotten
    return ticket < nextTicket_;
  }
  if (deviceInstance().api().isEventCompleted(pending->second.get())) {
    pendingEvents_.erase(pending);
    return true;
  }
  return false;
}

std::vector<parallel::runtime::StreamRuntime*> StreamExchangeScheduler::runtimes() const {
  std::vector<parallel::runtime::StreamRuntime*> result;
  for (const auto& runtime : runtimes_) {
    if (runtime != nullptr) {
      result.push_back(runtime.get());
    }
  }
  if (join_ != nullptr) {
    result.push_back(join_.get());
  }
  return result;
}

ExchangeScheduler::EventUsage StreamExchangeScheduler::eventUsage() const {
  // summed over the runtimes; the peaks need not have been at the same time
  EventUsage usage;
  for (const auto* runtime : runtimes()) {
    usage.pool += runtime->eventPool().size();
    usage.peak += runtime->eventPool().peakReferenced();
  }
  return usage;
}

#else

StreamExchangeScheduler::StreamExchangeScheduler(std::size_t clusterCount,
                                                 LaunchOrder order,
                                                 std::size_t configCount,
                                                 std::size_t lanes)
    : ExchangeScheduler(clusterCount, order, configCount), clusterCount_(clusterCount),
      lanes_(std::max<std::size_t>(lanes, 1)) {
  logError() << "Exchanging the halo data on device streams needs a device build.";
}

StreamExchangeScheduler::~StreamExchangeScheduler() = default;

void StreamExchangeScheduler::synchronize() {}

ExchangeScheduler::Ticket
    StreamExchangeScheduler::launch(std::size_t /*from*/,
                                    std::size_t /*to*/,
                                    std::size_t /*exchange*/,
                                    const ScheduledTransport* /*sender*/,
                                    const ScheduledTransport* /*receiver*/,
                                    const std::vector<ActorEvent>& /*after*/) {
  return 0;
}

bool StreamExchangeScheduler::completed(Ticket /*ticket*/) { return true; }

void StreamExchangeScheduler::isolateGroups() {}

std::vector<parallel::runtime::StreamRuntime*> StreamExchangeScheduler::runtimes() const {
  return {};
}

ExchangeScheduler::EventUsage StreamExchangeScheduler::eventUsage() const { return {}; }

#endif

} // namespace seissol::solver
