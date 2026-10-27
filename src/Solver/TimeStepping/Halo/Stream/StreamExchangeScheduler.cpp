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
#include "Solver/TimeStepping/Actor/ActorState.h"

#include <Device/device.h>
#include <memory>
#include <utility>
#endif

namespace seissol::solver {

std::size_t StreamExchangeScheduler::slot(std::size_t from, std::size_t to) const {
  return launchOrder() == LaunchOrder::Global ? 0 : from * clusterCount_ + to;
}

std::vector<std::size_t> StreamExchangeScheduler::usedSlots() const {
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
                                                 std::size_t configCount)
    : ExchangeScheduler(clusterCount, order, configCount), clusterCount_(clusterCount),
      runtimes_(order == LaunchOrder::Global ? 1 : clusterCount * clusterCount) {
  for (const auto current : usedSlots()) {
    runtimes_[current] = std::make_unique<parallel::runtime::StreamRuntime>();
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
}

void StreamExchangeScheduler::synchronize() {
  for (auto& runtime : runtimes_) {
    if (runtime != nullptr) {
      runtime->wait();
    }
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
  const auto current = slot(from, to);
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
  return result;
}

ExchangeScheduler::EventUsage StreamExchangeScheduler::eventUsage() const {
  // summed over the runtimes; the peaks need not have been at the same time
  EventUsage usage;
  for (const auto& runtime : runtimes_) {
    if (runtime != nullptr) {
      usage.pool += runtime->eventPool().size();
      usage.peak += runtime->eventPool().peakReferenced();
    }
  }
  return usage;
}

#else

StreamExchangeScheduler::StreamExchangeScheduler(std::size_t clusterCount,
                                                 LaunchOrder order,
                                                 std::size_t configCount)
    : ExchangeScheduler(clusterCount, order, configCount), clusterCount_(clusterCount) {
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
