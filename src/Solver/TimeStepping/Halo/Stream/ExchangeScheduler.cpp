// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ExchangeScheduler.h"

#include "Solver/TimeStepping/Halo/HaloCommunication.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iterator>
#include <tuple>
#include <utility>
#include <utils/logger.h>

namespace seissol::time_stepping {

namespace {

/**
 * Orders the regions by tag, and by peer within a tag; both sides of an exchange thus issue the
 * operations towards the same peer in the same order.
 */
void sortRegions(std::vector<solver::RemoteCluster>& regions) {
  std::sort(regions.begin(), regions.end(), [](const auto& a, const auto& b) {
    return a.tag < b.tag || (a.tag == b.tag && a.rank < b.rank);
  });
}

} // namespace

ExchangeScheduler::ExchangeScheduler(std::size_t clusterCount, LaunchOrder order)
    : clusterCount_(clusterCount), order_(order), directions_(clusterCount * clusterCount) {}

ExchangeScheduler::Direction& ExchangeScheduler::direction(std::size_t from, std::size_t to) {
  assert(from < clusterCount_ && to < clusterCount_);
  return directions_[from * clusterCount_ + to];
}

void ExchangeScheduler::add(ScheduledTransport& transport) {
  auto& outgoing = direction(transport.cluster(), transport.otherCluster());
  auto& incoming = direction(transport.otherCluster(), transport.cluster());
  if (outgoing.sender != nullptr || incoming.receiver != nullptr) {
    logError() << "Two halo transports of one process exchange between the clusters"
               << transport.cluster() << "and" << transport.otherCluster();
  }
  outgoing.sender = &transport;
  incoming.receiver = &transport;
  ++transports_;
  added(transport);
}

void ExchangeScheduler::startInterval(const ScheduledTransport& transport,
                                      const ExchangeInterval& interval) {
  auto& outgoing = direction(transport.cluster(), transport.otherCluster());
  outgoing.sendRate = interval.sendRate;
  outgoing.sendSteps = interval.sendSteps;
  outgoing.exchangePeriod = interval.exchangePeriod;
  auto& incoming = direction(transport.otherCluster(), transport.cluster());
  incoming.sendRate = interval.receiveRate;
  incoming.sendSteps = interval.receiveSteps;
  incoming.exchangePeriod = interval.exchangePeriod;

  if (order_ == LaunchOrder::Global) {
    if (announced_ == 0 && launched_ != sequence_.size()) {
      logError() << "Not all halo exchanges of the last interval went out:" << launched_ << "of"
                 << sequence_.size();
    }
    ++announced_;
    if (announced_ == transports_) {
      announced_ = 0;
      orderInterval();
      launchInOrder();
    }
  }
}

void ExchangeScheduler::orderInterval() {
  struct Group {
    long time;
    std::size_t from;
    std::size_t to;
  };
  std::vector<Group> groups;
  // (the ticks are counted in steps of the smallest cluster since the synchronization point)
  for (std::size_t from = 0; from < clusterCount_; ++from) {
    for (std::size_t to = 0; to < clusterCount_; ++to) {
      const auto& current = direction(from, to);
      if (current.sender == nullptr && current.receiver == nullptr) {
        continue;
      }
      // the sending cluster completes the data of exchange k with its prediction that ends at
      // min(k * period, final), where its final prediction ends with the last step it takes
      const auto rate = current.sendRate;
      const auto period = current.exchangePeriod;
      const auto final = (current.sendSteps + rate - 1) / rate * rate;
      const auto count = (current.sendSteps + period - 1) / period;
      for (long exchange = 1; exchange <= count; ++exchange) {
        groups.push_back({std::min(exchange * period, final) - rate, from, to});
      }
    }
  }
  std::sort(groups.begin(), groups.end(), [](const auto& a, const auto& b) {
    return std::tie(a.time, a.from, a.to) < std::tie(b.time, b.from, b.to);
  });

  sequence_.clear();
  for (const auto& group : groups) {
    sequence_.emplace_back(group.time, group.from, group.to);
  }
  launched_ = 0;
}

void ExchangeScheduler::setStreamOrdered(bool streamOrdered) {
  assert(!streamOrdered || order_ == LaunchOrder::Global);
  streamOrdered_ = streamOrdered;
}

void ExchangeScheduler::setHorizon(long time) {
  horizon_ = time;
  if (order_ == LaunchOrder::Global) {
    launchInOrder();
  }
}

bool ExchangeScheduler::launchedBefore(long time) const {
  if (order_ != LaunchOrder::Global || announced_ != 0) {
    return true;
  }
  // the groups are ordered by time
  const auto first = std::find_if(sequence_.begin(), sequence_.end(), [&](const auto& group) {
    return std::get<0>(group) >= time;
  });
  return launched_ >= static_cast<std::size_t>(std::distance(sequence_.begin(), first));
}

std::size_t ExchangeScheduler::readySend(const ScheduledTransport& transport, void* after) {
  auto& outgoing = direction(transport.cluster(), transport.otherCluster());
  assert(outgoing.sender == &transport);
  outgoing.sendsAfter.push_back(after);
  const auto exchange = outgoing.readySends++;
  launchReady(transport.cluster(), transport.otherCluster());
  return exchange;
}

std::size_t ExchangeScheduler::readyReceive(const ScheduledTransport& transport, void* after) {
  auto& incoming = direction(transport.otherCluster(), transport.cluster());
  assert(incoming.receiver == &transport);
  incoming.receivesAfter.push_back(after);
  const auto exchange = incoming.readyReceives++;
  launchReady(transport.otherCluster(), transport.cluster());
  return exchange;
}

bool ExchangeScheduler::ready(const Direction& direction) {
  return (direction.sender == nullptr || direction.readySends > direction.groups.size()) &&
         (direction.receiver == nullptr || direction.readyReceives > direction.groups.size());
}

void ExchangeScheduler::launchReady(std::size_t from, std::size_t to) {
  if (order_ == LaunchOrder::Global) {
    launchInOrder();
    return;
  }
  while (ready(direction(from, to))) {
    launchNext(from, to);
  }
}

void ExchangeScheduler::launchNext(std::size_t from, std::size_t to) {
  auto& current = direction(from, to);
  std::vector<void*> after;
  for (auto* queue : {&current.sendsAfter, &current.receivesAfter}) {
    if (!queue->empty()) {
      if (queue->front() != nullptr) {
        after.push_back(queue->front());
      }
      queue->pop_front();
    }
  }
  if (launching_) {
    const auto exchange = current.groups.size();
    current.groups.push_back(launch(from, to, exchange, current.sender, current.receiver, after));
  } else {
    // the operations come from elsewhere; the group only counts as launched
    current.groups.push_back(0);
  }
}

void ExchangeScheduler::launchInOrder() {
  // nothing goes out while the interval is still being announced
  while (announced_ == 0 && launched_ < sequence_.size()) {
    const auto& group = sequence_[launched_];
    const auto from = std::get<1>(group);
    const auto to = std::get<2>(group);
    if (std::get<0>(group) >= horizon_ || !ready(direction(from, to))) {
      break;
    }
    launchNext(from, to);
    ++launched_;
  }
}

bool ExchangeScheduler::groupCompleted(Direction& direction, std::size_t exchange) {
  // ordered on the device, the work that depends on a group waits for it there
  return exchange < direction.groups.size() &&
         (streamOrdered_ || !launching_ || completed(direction.groups[exchange]));
}

bool ExchangeScheduler::sendCompleted(const ScheduledTransport& transport, std::size_t exchange) {
  return groupCompleted(direction(transport.cluster(), transport.otherCluster()), exchange);
}

bool ExchangeScheduler::receiveCompleted(const ScheduledTransport& transport,
                                         std::size_t exchange) {
  return groupCompleted(direction(transport.otherCluster(), transport.cluster()), exchange);
}

ScheduledTransport::ScheduledTransport(ExchangeScheduler& scheduler,
                                       const solver::RemoteClusterPair& regions,
                                       std::size_t cluster,
                                       std::size_t otherCluster)
    : scheduler_(scheduler), regions_(regions), cluster_(cluster), otherCluster_(otherCluster) {
  sortRegions(regions_.copy);
  sortRegions(regions_.ghost);
  scheduler_.add(*this);
}

void ScheduledTransport::startInterval(const ExchangeInterval& interval) {
  scheduler_.startInterval(*this, interval);
}

void ScheduledTransport::startSendAfter(void* event) {
  assert(!sending_);
  sending_ = true;
  sendExchange_ = scheduler_.readySend(*this, event);
}

void ScheduledTransport::startReceiveAfter(void* event) {
  assert(!receiving_);
  receiving_ = true;
  receiveExchange_ = scheduler_.readyReceive(*this, event);
}

void ScheduledTransport::startSend() {
  assert(!sending_);
  sending_ = true;
  sendExchange_ = scheduler_.readySend(*this);
}

bool ScheduledTransport::testSend() {
  if (sending_ && scheduler_.sendCompleted(*this, sendExchange_)) {
    sending_ = false;
  }
  return !sending_;
}

void ScheduledTransport::startReceive() {
  assert(!receiving_);
  receiving_ = true;
  receiveExchange_ = scheduler_.readyReceive(*this);
}

bool ScheduledTransport::testReceive() {
  if (receiving_ && scheduler_.receiveCompleted(*this, receiveExchange_)) {
    receiving_ = false;
  }
  return !receiving_;
}

} // namespace seissol::time_stepping
