// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Solver/TimeStepping/AbstractGhostTimeCluster.h"
#include <Common/Executor.h>
#include <Kernels/Common.h>
#include <Monitoring/Instrumentation.h>
#include <Solver/TimeStepping/AbstractTimeCluster.h>
#include <Solver/TimeStepping/ActorState.h>
#include <Solver/TimeStepping/HaloCommunication.h>
#include <cassert>
#include <chrono>
#include <list>
#include <mpi.h>
#include <string>
#include <utils/logger.h>

namespace seissol::time_stepping {
bool AbstractGhostTimeCluster::testQueue(MPI_Request* requests, std::list<unsigned int>& regions) {
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
  return testQueue(sendRequests.data(), sendQueue);
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

void AbstractGhostTimeCluster::handleAdvancedPredictionTimeMessage(
    const NeighborCluster& /*neighborCluster*/) {
  assert(testForCopyLayerSends());
  sendCopyLayer();
}

void AbstractGhostTimeCluster::handleAdvancedCorrectionTimeMessage(
    const NeighborCluster& neighborCluster) {
  assert(testForGhostLayerReceives());

  auto upcomingCorrectionSteps = ct.stepsSinceLastSync;
  if (state == ActorState::Predicted) {
    upcomingCorrectionSteps = ct.nextCorrectionSteps();
  }

  const bool ignoreMessage = upcomingCorrectionSteps >= ct.stepsUntilSync;

  // If we are already at a sync point, we must not post an additional receive, as otherwise start()
  // posts an additional request! This is also true for the last sync point (i.e. end of
  // simulation), as in this case we do not want to have any hanging request.
  if (!ignoreMessage) {
    receiveGhostLayer();
  }
}

AbstractGhostTimeCluster::AbstractGhostTimeCluster(
    double maxTimeStepSize,
    int timeStepRate,
    int globalTimeClusterId,
    int otherGlobalTimeClusterId,
    const seissol::solver::HaloCommunication& meshStructure)
    : AbstractTimeCluster(
          maxTimeStepSize, timeStepRate, isDeviceOn() ? Executor::Device : Executor::Host),
      globalClusterId(globalTimeClusterId), otherGlobalClusterId(otherGlobalTimeClusterId),
      meshStructure(meshStructure.at(globalTimeClusterId).at(otherGlobalTimeClusterId)),
      sendRequests(meshStructure.at(globalTimeClusterId).at(otherGlobalTimeClusterId).copy.size()),
      recvRequests(
          meshStructure.at(globalTimeClusterId).at(otherGlobalTimeClusterId).ghost.size()) {}

void AbstractGhostTimeCluster::reset() {
  AbstractTimeCluster::reset();
  assert(testForGhostLayerReceives());
  lastSendTime = -1;
}

void AbstractGhostTimeCluster::printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) {
  logError() << "Ghost: No update since " << timeSinceLastUpdate.count()
             << "[s] for global cluster " << globalClusterId << " with other cluster id "
             << otherGlobalClusterId << " at state " << actorStateToString(state)
             << " mayPredict = " << mayPredict()
             << " mayPredict (steps) = " << AbstractTimeCluster::mayPredict()
             << " mayCorrect = " << mayCorrect()
             << " mayCorrect (steps) = " << AbstractTimeCluster::mayCorrect()
             << " maySync = " << maySync();
  for (auto& neighbor : neighbors) {
    logError() << "Neighbor with rate = " << neighbor.ct.timeStepRate
               << "PredTime = " << neighbor.ct.predictionTime
               << "CorrTime = " << neighbor.ct.correctionTime
               << "predictionsSinceSync = " << neighbor.ct.predictionsSinceLastSync
               << "correctionsSinceSync = " << neighbor.ct.stepsSinceLastSync;
  }
}

std::string AbstractGhostTimeCluster::description() const {
  return "comm-" + std::to_string(globalClusterId) + "-" + std::to_string(otherGlobalClusterId);
}

} // namespace seissol::time_stepping
