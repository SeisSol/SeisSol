// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Solver/TimeStepping/AbstractGhostTimeCluster.h"

#include "Common/Executor.h"
#include "Kernels/Common.h"
#include "Monitoring/Instrumentation.h"
#include "Solver/TimeStepping/AbstractTimeCluster.h"
#include "Solver/TimeStepping/ActorState.h"
#include "Solver/TimeStepping/HaloCommunication.h"

#include <atomic>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <list>
#include <mpi.h>
#include <string>

namespace seissol::time_stepping {
bool AbstractGhostTimeCluster::testQueue(MPI_Request* requests, std::list<std::size_t>& regions) {
  for (auto region = regions.begin(); region != regions.end();) {
    MPI_Request* request = &requests[*region];
    int testSuccess = 0;
    MPI_Test(request, &testSuccess, MPI_STATUS_IGNORE);
    if (testSuccess != 0) {
      region = regions.erase(region);
    } else {
      ++region;
    }
  }
  return regions.empty();
}

bool AbstractGhostTimeCluster::testForCopyLayerSends() {
  SCOREP_USER_REGION("testForCopyLayerSends", SCOREP_USER_REGION_TYPE_FUNCTION)
  return testQueue(sendRequests_.data(), sendQueue_);
}

ActResult AbstractGhostTimeCluster::act() {
  // Always check for receives/send for quicker MPI progression.
  testForGhostLayerReceives();
  testForCopyLayerSends();
  return AbstractTimeCluster::act();
}

void AbstractGhostTimeCluster::start() {
  assert(testForGhostLayerReceives());
  receiveGhostLayer();
}

void AbstractGhostTimeCluster::predict() {
  // Doesn't do anything
}

void AbstractGhostTimeCluster::correct() {
  // Doesn't do anything
}

bool AbstractGhostTimeCluster::mayCorrect() {
  return testForCopyLayerSends() && AbstractTimeCluster::mayCorrect();
}

bool AbstractGhostTimeCluster::mayPredict() {
  return testForGhostLayerReceives() && AbstractTimeCluster::mayPredict();
}

bool AbstractGhostTimeCluster::maySync() {
  return testForGhostLayerReceives() && testForCopyLayerSends() && AbstractTimeCluster::maySync();
}

void AbstractGhostTimeCluster::handleNeighborPrediction(const NeighborCluster& neighbor) {
  // The copy layer has new data for the remote cluster once it has predicted at least up to the end
  // of the current step of the remote cluster, and at the synchronization point.
  const bool copyAtSync = neighbor.progress->stepsUntilSync.load(std::memory_order_relaxed) <=
                          neighbor.ct.predictionsSinceLastSync;
  if (copyAtSync || neighbor.ct.predictionsSinceLastSync >= ct_.nextCorrectionSteps()) {
    assert(testForCopyLayerSends());
    sendCopyLayer();
  }
}

void AbstractGhostTimeCluster::handleNeighborCorrection(const NeighborCluster& neighbor) {
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

AbstractGhostTimeCluster::AbstractGhostTimeCluster(
    double maxTimeStepSize,
    std::uint64_t timeStepRate,
    std::size_t globalTimeClusterId,
    std::size_t otherGlobalTimeClusterId,
    const std::string& displayName,
    const std::string& otherDisplayName,
    const seissol::solver::HaloCommunication& meshStructure)
    : AbstractTimeCluster(
          maxTimeStepSize, timeStepRate, isDeviceOn() ? Executor::Device : Executor::Host),
      globalClusterId_(globalTimeClusterId), otherGlobalClusterId_(otherGlobalTimeClusterId),
      meshStructure_(meshStructure.at(globalTimeClusterId).at(otherGlobalTimeClusterId)),
      sendRequests_(meshStructure.at(globalTimeClusterId).at(otherGlobalTimeClusterId).copy.size()),
      recvRequests_(
          meshStructure.at(globalTimeClusterId).at(otherGlobalTimeClusterId).ghost.size()),
      displayName_(displayName), otherDisplayName_(otherDisplayName) {}

void AbstractGhostTimeCluster::reset() {
  AbstractTimeCluster::reset();
  assert(testForGhostLayerReceives());
  lastSendTime_ = -1;
}

std::string AbstractGhostTimeCluster::description() const {
  return "comm-" + displayName_ + "-" + otherDisplayName_;
}

bool AbstractGhostTimeCluster::timeoutFail() const { return true; }

} // namespace seissol::time_stepping
