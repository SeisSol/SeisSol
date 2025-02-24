// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <iostream>
#include <cmath>
#include <memory>
#include <algorithm>
#include <type_traits>

#include "ActorState.h"

namespace seissol::time_stepping {

inline std::ostream& operator<<(std::ostream& stream, const Message& message) {
  std::visit([&stream](auto&& msg) {
    using T = std::decay_t<decltype(msg)>;
    if constexpr (std::is_same_v<T, AdvancedPredictionTimeMessage>) {
      stream << "AdvancedPredictionTimeMessage, t = " << msg.time;
    } else if constexpr (std::is_same_v<T, AdvancedCorrectionTimeMessage>) {
      stream << "AdvancedCorrectionTimeMessage, t = " << msg.time;
    } else {
      static_assert(always_false<T>::value, "non-exhaustive visitor");
    }
  }, message);
  return stream;
}

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


void MessageQueue::push(const Message& message) {
  std::lock_guard lock{mutex};
  queue.push(message);
}

Message MessageQueue::pop() {
  std::lock_guard lock{mutex};
  const Message message = queue.front();
  queue.pop();
  return message;
}

bool MessageQueue::hasMessages() const {
  return !queue.empty();
}

size_t MessageQueue::size() const {
  return queue.size();
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
  const double timeDiff = newSyncTime-oldSyncTime;
  return static_cast<long>(std::ceil(timeStepRate*timeDiff/maxTimeStepSize));
}

NeighborCluster::NeighborCluster(double maxTimeStepSize, int timeStepRate, Executor neighborExecutor) {
  ct.maxTimeStepSize = maxTimeStepSize;
  ct.timeStepRate = timeStepRate;
  executor = neighborExecutor;
}

DynamicRuptureScheduler::DynamicRuptureScheduler(long numberOfDynamicRuptureFaces, bool isFirstDynamicRuptureCluster) :
    numberOfDynamicRuptureFaces(numberOfDynamicRuptureFaces),
      firstClusterWithDynamicRuptureFaces(isFirstDynamicRuptureCluster) {}

bool DynamicRuptureScheduler::mayComputeInterior(long curCorrectionSteps) const {
  return curCorrectionSteps > lastCorrectionStepsInterior;
}

bool DynamicRuptureScheduler::mayComputeFaultOutput(long curCorrectionSteps) const {
  return curCorrectionSteps == lastCorrectionStepsInterior
         && curCorrectionSteps == lastCorrectionStepsCopy
         && curCorrectionSteps > lastFaultOutput;
}

void DynamicRuptureScheduler::setLastCorrectionStepsInterior(long steps) {
  lastCorrectionStepsInterior = steps;
}

void DynamicRuptureScheduler::setLastCorrectionStepsCopy(long steps) {
  lastCorrectionStepsCopy = steps;
}

void DynamicRuptureScheduler::setLastFaultOutput(long steps) {
  lastFaultOutput = steps;
}

bool DynamicRuptureScheduler::hasDynamicRuptureFaces() const {
  return numberOfDynamicRuptureFaces > 0;
}


bool DynamicRuptureScheduler::isFirstClusterWithDynamicRuptureFaces() const {
  return firstClusterWithDynamicRuptureFaces;
}
} // namespace seissol::time_stepping

