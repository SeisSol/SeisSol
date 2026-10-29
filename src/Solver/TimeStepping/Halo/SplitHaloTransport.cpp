// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "SplitHaloTransport.h"

#include "Solver/TimeStepping/Actor/ActorState.h"
#include "Solver/TimeStepping/Halo/HaloCommunication.h"
#include "Solver/TimeStepping/Halo/HaloTransport.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::solver {

namespace {
bool completed([[maybe_unused]] const ActorEvent& event) {
#ifdef ACL_DEVICE
  return !event || device::DeviceInstance::instance().api().isEventCompleted(event.get());
#else
  return true;
#endif
}
} // namespace

SplitHaloTransport::SplitHaloTransport(std::unique_ptr<HaloTransport> local,
                                       std::unique_ptr<HaloTransport> remote)
    : local_(std::move(local)), remote_(std::move(remote)) {}

std::pair<RemoteClusterPair, RemoteClusterPair>
    SplitHaloTransport::split(const RemoteClusterPair& regions, const std::vector<bool>& local) {
  RemoteClusterPair onNode;
  RemoteClusterPair offNode;
  for (const auto& region : regions.copy) {
    (local.at(region.rank) ? onNode : offNode).copy.push_back(region);
  }
  for (const auto& region : regions.ghost) {
    (local.at(region.rank) ? onNode : offNode).ghost.push_back(region);
  }
  return {onNode, offNode};
}

void SplitHaloTransport::startSend() {
  local_->startSend();
  remote_->startSend();
}

void SplitHaloTransport::startReceive() {
  local_->startReceive();
  remote_->startReceive();
}

void SplitHaloTransport::startSendAfter(const ActorEvent& event) {
  local_->startSendAfter(event);
  // the host sends once the device has written the copy layer
  sendDeferred_ = true;
  sendAfter_ = event;
  startDeferred();
}

void SplitHaloTransport::startReceiveAfter(const ActorEvent& event) {
  local_->startReceiveAfter(event);
  // the host receives once the device has read the ghost layer
  receiveDeferred_ = true;
  receiveAfter_ = event;
  startDeferred();
}

void SplitHaloTransport::startDeferred() {
  if (sendDeferred_ && completed(sendAfter_)) {
    sendDeferred_ = false;
    sendAfter_ = ActorEvent();
    remote_->startSend();
  }
  if (receiveDeferred_ && completed(receiveAfter_)) {
    receiveDeferred_ = false;
    receiveAfter_ = ActorEvent();
    remote_->startReceive();
  }
}

bool SplitHaloTransport::testSend() {
  startDeferred();
  // both get tested, so that both progress
  const bool local = local_->testSend();
  const bool remote = !sendDeferred_ && remote_->testSend();
  return local && remote;
}

bool SplitHaloTransport::testReceive() {
  startDeferred();
  const bool local = local_->testReceive();
  const bool remote = !receiveDeferred_ && remote_->testReceive();
  return local && remote;
}

std::string SplitHaloTransport::describe() const {
  std::string description = "on the node: " + local_->describe();
  description += "; to other nodes: " + remote_->describe();
  if (sendDeferred_ || receiveDeferred_) {
    description += " (waiting for the device to start";
    description += sendDeferred_ ? " the send" : "";
    description += sendDeferred_ && receiveDeferred_ ? " and" : "";
    description += receiveDeferred_ ? " the receive" : "";
    description += ')';
  }
  return description;
}

void SplitHaloTransport::startInterval(const ExchangeInterval& interval) {
  local_->startInterval(interval);
  remote_->startInterval(interval);
}

void SplitHaloTransport::finalize() {
  local_->finalize();
  remote_->finalize();
}

void SplitHaloTransport::setTraceName(const std::string& name) {
  local_->setTraceName(name);
  remote_->setTraceName(name);
}

} // namespace seissol::solver
