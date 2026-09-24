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
#include <cstddef>
#include <string>

namespace seissol::time_stepping {

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

double ClusterTimes::timeStepSize(double syncTime) const {
  return std::min(syncTime - correctionTime, maxTimeStepSize);
}

long ClusterTimes::computeStepsUntilSyncTime(double oldSyncTime, double newSyncTime) const {
  const double timeDiff = newSyncTime - oldSyncTime;
  return static_cast<long>(std::ceil(timeStepRate * timeDiff / maxTimeStepSize));
}

NeighborCluster::NeighborCluster(double maxTimeStepSize,
                                 int timeStepRate,
                                 Executor neighborExecutor)
    : executor(neighborExecutor) {
  ct.maxTimeStepSize = maxTimeStepSize;
  ct.timeStepRate = timeStepRate;
}

} // namespace seissol::time_stepping
