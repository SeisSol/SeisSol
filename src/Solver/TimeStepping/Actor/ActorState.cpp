// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ActorState.h"

#include "Common/Executor.h"

#include <algorithm>
#include <cmath>
#include <string>

namespace seissol::solver {

std::string actorStateToString(ActorState state) {
  switch (state) {
  case ActorState::Corrected:
    return "Corrected";
  case ActorState::Predicted:
    return "Predicted";
  case ActorState::Synced:
    return "Synced";
  }
  throw;
}

double ClusterTimes::nextCorrectionTime(double syncTime) const {
  return std::min(syncTime, correctionTime + maxTimeStepSize);
}

long ClusterTimes::nextCorrectionSteps() const {
  return std::min(stepsSinceLastSync + timeStepRate, stepsUntilSync);
}

double ClusterTimes::timeStepSize(double time, long steps, double syncTime) const {
  if (isLastStep(steps)) {
    // may exceed the maximum by the tolerance of the step count
    return syncTime - time;
  }
  return std::min(syncTime - time, maxTimeStepSize);
}

double ClusterTimes::timeStepSize(double syncTime) const {
  return timeStepSize(correctionTime, stepsSinceLastSync, syncTime);
}

long ClusterTimes::computeStepsUntilSyncTime(double oldSyncTime, double newSyncTime) const {
  const double timeDiff = newSyncTime - oldSyncTime;
  if (timeDiff <= 0) {
    return 0;
  }
  const double ticks = timeStepRate * timeDiff / maxTimeStepSize;
  return std::max(1L, static_cast<long>(std::ceil(ticks - TickTolerance)));
}

NeighborCluster::NeighborCluster(double maxTimeStepSize,
                                 int timeStepRate,
                                 Executor neighborExecutor)
    : executor(neighborExecutor) {
  ct.maxTimeStepSize = maxTimeStepSize;
  ct.timeStepRate = timeStepRate;
}

} // namespace seissol::solver
