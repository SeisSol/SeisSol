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
#include "Solver/TimeStepping/AbstractTimeCluster.h"
#include "Solver/TimeStepping/ActorState.h"
#include "Solver/TimeStepping/HaloCommunication.h"
#include "Solver/TimeStepping/HaloTransport.h"

#include <atomic>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>

namespace seissol::time_stepping {
void GhostCluster::sendCopyLayer() {
  SCOREP_USER_REGION("sendCopyLayer", SCOREP_USER_REGION_TYPE_FUNCTION)
  assert(ct_.correctionTime > lastSendTime_);
  lastSendTime_ = ct_.correctionTime;
  transport_->startSend();
  sentMessages_ += copyRegionCount_;
}

void GhostCluster::receiveGhostLayer() {
  SCOREP_USER_REGION("receiveGhostLayer", SCOREP_USER_REGION_TYPE_FUNCTION)
  assert(ct_.predictionTime >= lastSendTime_);
  transport_->startReceive();
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

ActResult GhostCluster::act() {
  // Always check for receives/send for quicker MPI progression.
  testForGhostLayerReceives();
  testForCopyLayerSends();
  return AbstractTimeCluster::act();
}

void GhostCluster::start() {
  assert(testForGhostLayerReceives());
  receiveGhostLayer();
}

void GhostCluster::predict() {
  // Doesn't do anything
}

void GhostCluster::correct() {
  // Doesn't do anything
}

bool GhostCluster::mayCorrect() {
  return testForCopyLayerSends() && AbstractTimeCluster::mayCorrect();
}

bool GhostCluster::mayPredict() {
  return testForGhostLayerReceives() && AbstractTimeCluster::mayPredict();
}

bool GhostCluster::maySync() {
  return testForGhostLayerReceives() && testForCopyLayerSends() && AbstractTimeCluster::maySync();
}

void GhostCluster::handleNeighborPrediction(const NeighborCluster& neighbor) {
  // The copy layer has new data for the remote cluster once it has predicted at least up to the end
  // of the current step of the remote cluster, and at the synchronization point.
  const bool copyAtSync = neighbor.progress->stepsUntilSync.load(std::memory_order_relaxed) <=
                          neighbor.ct.predictionsSinceLastSync;
  if (copyAtSync || neighbor.ct.predictionsSinceLastSync >= ct_.nextCorrectionSteps()) {
    assert(testForCopyLayerSends());
    sendCopyLayer();
  }
}

void GhostCluster::handleNeighborCorrection(const NeighborCluster& neighbor) {
  // The ghost layer may be overwritten once the copy layer has corrected up to the data received
  // last, and at the synchronization point. Before a correction, the copy layer has predicted
  // exactly as far as it has corrected afterwards.
  const bool copyAtSync = neighbor.progress->stepsUntilSync.load(std::memory_order_relaxed) <=
                          neighbor.ct.stepsSinceLastSync;
  if (!copyAtSync && neighbor.ct.stepsSinceLastSync < ct_.predictionsSinceLastSync) {
    return;
  }

  assert(testForGhostLayerReceives());
  auto upcomingCorrectionSteps = ct_.stepsSinceLastSync;
  if (state_ == ActorState::Predicted) {
    upcomingCorrectionSteps = ct_.nextCorrectionSteps();
  }
  const bool atSync = upcomingCorrectionSteps >= ct_.stepsUntilSync;
  // If we are already at a sync point, we must not post an additional receive, as otherwise start()
  // posts an additional request! This is also true for the last sync point (i.e. end of
  // simulation), as in this case we do not want to have any hanging request.
  if (!atSync) {
    receiveGhostLayer();
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
  AbstractTimeCluster::reset();
  assert(testForGhostLayerReceives());
  lastSendTime_ = -1;
}

std::string GhostCluster::description() const {
  return "comm-" + displayName_ + "-" + otherDisplayName_;
}

bool GhostCluster::timeoutFail() const { return true; }

void GhostCluster::finalize() { transport_->finalize(); }

std::size_t GhostCluster::sentMessages() const { return sentMessages_; }

std::size_t GhostCluster::receivedMessages() const { return receivedMessages_; }

} // namespace seissol::time_stepping
