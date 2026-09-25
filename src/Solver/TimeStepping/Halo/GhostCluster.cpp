// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "GhostCluster.h"

#include "Common/Executor.h"
#include "Kernels/Common.h"
#include "Monitoring/Instrumentation.h"
#include "Solver/TimeStepping/Actor/AbstractTimeCluster.h"
#include "Solver/TimeStepping/Actor/ActorState.h"
#include "Solver/TimeStepping/Halo/HaloCommunication.h"
#include "Solver/TimeStepping/Halo/HaloTransport.h"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <utils/logger.h>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::time_stepping {
void GhostCluster::sendCopyLayer(long target, void* after) {
  SCOREP_USER_REGION("sendCopyLayer", SCOREP_USER_REGION_TYPE_FUNCTION)
  assert(!sending_);
  sending_ = true;
  sendTarget_ = target;
  transport_->startSendAfter(after);
  sentMessages_ += copyRegionCount_;
}

void GhostCluster::receiveGhostLayer(long target, void* after) {
  SCOREP_USER_REGION("receiveGhostLayer", SCOREP_USER_REGION_TYPE_FUNCTION)
  assert(!receiving_);
  receiving_ = true;
  receiveTarget_ = target;
  transport_->startReceiveAfter(after);
  receivedMessages_ += ghostRegionCount_;
}

bool GhostCluster::testForCopyLayerSends() {
  SCOREP_USER_REGION("testForCopyLayerSends", SCOREP_USER_REGION_TYPE_FUNCTION)
  return transport_->testSend();
}

bool GhostCluster::testForGhostLayerReceives() {
  SCOREP_USER_REGION("testForGhostLayerReceives", SCOREP_USER_REGION_TYPE_FUNCTION)
  return transport_->testReceive();
}

long GhostCluster::finalSteps() const {
  return (ct_.stepsUntilSync + ct_.timeStepRate - 1) / ct_.timeStepRate * ct_.timeStepRate;
}

long GhostCluster::exchangePeriod() const {
  // the only cluster this one follows is the copy layer
  assert(neighbors_.size() == 1);
  return std::max(ct_.timeStepRate, neighbors_.front().ct.timeStepRate);
}

void GhostCluster::advanceReceived(long target) {
  while (ct_.predictionsSinceLastSync < target) {
    ct_.predictionTime += std::min(syncTime_ - ct_.predictionTime, ct_.maxTimeStepSize);
    ct_.predictionsSinceLastSync += ct_.timeStepRate;
    ct_.predictionsSinceStart += ct_.timeStepRate;
  }
  publishTransportEvent();
  publishProgress();
}

void GhostCluster::advanceSent(long target) {
  while (ct_.stepsSinceLastSync < target) {
    ct_.correctionTime += timeStepSize();
    ct_.stepsSinceLastSync += ct_.timeStepRate;
    ct_.stepsSinceStart += ct_.timeStepRate;
    ++numberOfTimeSteps_;
  }
  publishTransportEvent();
  publishProgress();
}

void GhostCluster::publishTransportEvent() {
  if (transport_->streamOrdered()) {
    // the copy layer waits for the exchanges on the device
    publishEvent(transport_->latestEvent());
  }
}

ActResult GhostCluster::act() {
  const auto stateBefore = state_;
  const auto receivedBefore = ct_.predictionsSinceLastSync;
  const auto sentBefore = ct_.stepsSinceLastSync;

  if (state_ == ActorState::Synced) {
    // restart after the synchronization point, once reset
    if (ct_.stepsSinceLastSync == 0) {
      start();
      state_ = ActorState::Corrected;
    }
  } else {
    // the progress of the copy layer starts the sends and receives
    refreshNeighbors();
    startDeferred();

    if (receiving_ && testForGhostLayerReceives()) {
      receiving_ = false;
      advanceReceived(receiveTarget_);
    }
    if (sending_ && testForCopyLayerSends()) {
      sending_ = false;
      advanceSent(sendTarget_);
    }

    const auto finalSteps = this->finalSteps();
    if (!receiving_ && !sending_ && !receiveDeferred_ && !sendDeferred_ &&
        ct_.predictionsSinceLastSync >= finalSteps && ct_.stepsSinceLastSync >= finalSteps) {
      state_ = ActorState::Synced;
    }
  }

  ActResult result;
  result.isStateChanged = state_ != stateBefore || ct_.predictionsSinceLastSync != receivedBefore ||
                          ct_.stepsSinceLastSync != sentBefore;
  trackProgress(result.isStateChanged);
  return result;
}

void GhostCluster::startDeferred() {
  const auto completed = [](void* event) {
#ifdef ACL_DEVICE
    return event == nullptr || device::DeviceInstance::getInstance().api->isEventCompleted(event);
#else
    return true;
#endif
  };
  if (sendDeferred_ && completed(deferredSendEvent_)) {
    sendDeferred_ = false;
    sendCopyLayer(deferredSendTarget_);
  }
  if (receiveDeferred_ && completed(deferredReceiveEvent_)) {
    receiveDeferred_ = false;
    receiveGhostLayer(deferredReceiveTarget_);
  }
}

void GhostCluster::start() {
  receiveGhostLayer(std::min(exchangePeriod(), finalSteps()),
                    neighbors_.front().progress->event.load(std::memory_order_relaxed));
}

void GhostCluster::handleNeighborPrediction(const NeighborCluster& neighbor) {
  // The copy layer has new data for the remote cluster once it has predicted at least up to the end
  // of the current step of the remote cluster, and at the synchronization point.
  const auto predictions = neighbor.ct.predictionsSinceLastSync;
  const bool copyAtSync =
      neighbor.progress->stepsUntilSync.load(std::memory_order_relaxed) <= predictions;
  if (copyAtSync || predictions >= ct_.nextCorrectionSteps()) {
    const auto rate = ct_.timeStepRate;
    const auto target = std::min((predictions + rate - 1) / rate * rate, finalSteps());
    if (transport_->streamOrdered()) {
      // the send waits on the device until the copy layer has written the data
      sendCopyLayer(target, neighbor.progress->event.load(std::memory_order_relaxed));
    } else if (concurrent()) {
      assert(!sendDeferred_);
      sendDeferred_ = true;
      deferredSendTarget_ = target;
      deferredSendEvent_ = neighbor.progress->event.load(std::memory_order_relaxed);
    } else {
      sendCopyLayer(target);
    }
  }
}

void GhostCluster::handleNeighborCorrection(const NeighborCluster& neighbor) {
  // The ghost layer may be overwritten once the copy layer has corrected up to the data received
  // last, and at the synchronization point. Before a correction, the copy layer has predicted
  // exactly as far as it has corrected afterwards.
  const auto corrections = neighbor.ct.stepsSinceLastSync;
  const bool copyAtSync =
      neighbor.progress->stepsUntilSync.load(std::memory_order_relaxed) <= corrections;
  if (!copyAtSync && corrections < ct_.predictionsSinceLastSync) {
    return;
  }

  // Once everything up to the synchronization point has arrived, no further receive is posted;
  // start() posts the first one of the next interval.
  const auto finalSteps = this->finalSteps();
  if (ct_.predictionsSinceLastSync < finalSteps) {
    const auto target = std::min(ct_.predictionsSinceLastSync + exchangePeriod(), finalSteps);
    if (transport_->streamOrdered()) {
      // the receive waits on the device until the copy layer has read the last ghost data
      receiveGhostLayer(target, neighbor.progress->event.load(std::memory_order_relaxed));
    } else if (concurrent()) {
      assert(!receiveDeferred_);
      receiveDeferred_ = true;
      deferredReceiveTarget_ = target;
      deferredReceiveEvent_ = neighbor.progress->event.load(std::memory_order_relaxed);
    } else {
      receiveGhostLayer(target);
    }
  }
}

GhostCluster::GhostCluster(double maxTimeStepSize,
                           std::uint64_t timeStepRate,
                           const std::string& displayName,
                           const std::string& otherDisplayName,
                           const solver::RemoteClusterPair& regions,
                           std::unique_ptr<HaloTransport> transport)
    : AbstractTimeCluster(
          maxTimeStepSize, timeStepRate, isDeviceOn() ? Executor::Device : Executor::Host),
      transport_(std::move(transport)), copyRegionCount_(regions.copy.size()),
      ghostRegionCount_(regions.ghost.size()), displayName_(displayName),
      otherDisplayName_(otherDisplayName) {}

void GhostCluster::reset() {
  // all transfers of the previous interval have completed when this cluster synchronized
  assert(!receiving_ && !sending_ && !receiveDeferred_ && !sendDeferred_);
  AbstractTimeCluster::reset();

  // the copy layer has been reset before and published its steps until the synchronization point
  const auto& copy = neighbors_.front();
  transport_->startInterval({copy.ct.timeStepRate,
                             copy.progress->stepsUntilSync.load(std::memory_order_relaxed),
                             ct_.timeStepRate,
                             ct_.stepsUntilSync,
                             exchangePeriod()});
}

std::string GhostCluster::description() const {
  return "comm-" + displayName_ + "-" + otherDisplayName_;
}

bool GhostCluster::timeoutFail() const { return true; }

void GhostCluster::printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) {
  logWarning(true) << "Cluster" << description() << "received up to" << ct_.predictionsSinceLastSync
                   << (receiving_ ? "(receiving up to " + std::to_string(receiveTarget_) + ")"
                                  : std::string("(not receiving)"))
                   << "and sent up to" << ct_.stepsSinceLastSync
                   << (sending_ ? "(sending up to " + std::to_string(sendTarget_) + ")"
                                : std::string("(not sending)"))
                   << "of" << finalSteps() << "steps.";
  AbstractTimeCluster::printTimeoutMessage(timeSinceLastUpdate);
}

void GhostCluster::finalize() { transport_->finalize(); }

std::size_t GhostCluster::sentMessages() const { return sentMessages_; }

std::size_t GhostCluster::receivedMessages() const { return receivedMessages_; }

} // namespace seissol::time_stepping
