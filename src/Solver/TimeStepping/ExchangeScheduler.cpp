// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ExchangeScheduler.h"

#include "Solver/TimeStepping/HaloCommunication.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
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

ExchangeScheduler::ExchangeScheduler(std::size_t clusterCount)
    : clusterCount_(clusterCount), directions_(clusterCount * clusterCount) {}

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
  added(transport);
}

std::size_t ExchangeScheduler::readySend(const ScheduledTransport& transport) {
  auto& outgoing = direction(transport.cluster(), transport.otherCluster());
  assert(outgoing.sender == &transport);
  const auto exchange = outgoing.readySends++;
  launchReady(transport.cluster(), transport.otherCluster());
  return exchange;
}

std::size_t ExchangeScheduler::readyReceive(const ScheduledTransport& transport) {
  auto& incoming = direction(transport.otherCluster(), transport.cluster());
  assert(incoming.receiver == &transport);
  const auto exchange = incoming.readyReceives++;
  launchReady(transport.otherCluster(), transport.cluster());
  return exchange;
}

void ExchangeScheduler::launchReady(std::size_t from, std::size_t to) {
  auto& current = direction(from, to);
  while ((current.sender == nullptr || current.readySends > current.groups.size()) &&
         (current.receiver == nullptr || current.readyReceives > current.groups.size())) {
    current.groups.push_back(launch(from, to, current.sender, current.receiver));
  }
}

bool ExchangeScheduler::groupCompleted(Direction& direction, std::size_t exchange) {
  return exchange < direction.groups.size() && completed(direction.groups[exchange]);
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
