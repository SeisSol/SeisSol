// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ActorState.h"

#include "Common/Executor.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <mutex>
#include <string>
#include <type_traits>
#include <variant>

namespace seissol::time_stepping {

inline std::ostream& operator<<(std::ostream& stream, const Message& message) {
  std::visit(
      [&stream](auto&& msg) {
        using T = std::decay_t<decltype(msg)>;
        if constexpr (std::is_same_v<T, AdvancedPredictionTimeMessage>) {
          stream << "AdvancedPredictionTimeMessage, t = " << msg.time;
        } else if constexpr (std::is_same_v<T, AdvancedCorrectionTimeMessage>) {
          stream << "AdvancedCorrectionTimeMessage, t = " << msg.time;
        } else {
          static_assert(sizeof(T) == 0, "non-exhaustive visitor");
        }
      },
      message);
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
  std::abort();
}

void MessageQueue::push(const Message& message) {
  const std::scoped_lock lock{mutex_};
  queue_.push(message);
  size_.store(queue_.size(), std::memory_order_release);
}

Message MessageQueue::pop() {
  const std::scoped_lock lock{mutex_};
  const Message message = queue_.front();
  queue_.pop();
  size_.store(queue_.size(), std::memory_order_release);
  return message;
}

bool MessageQueue::hasMessages() const { return size_.load(std::memory_order_acquire) > 0; }

size_t MessageQueue::size() const { return size_.load(std::memory_order_acquire); }

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

DynamicRuptureScheduler::DynamicRuptureScheduler(long numberOfDynamicRuptureFaces,
                                                 double outputTimestep)
    : numberOfDynamicRuptureFaces_(numberOfDynamicRuptureFaces), outputTimestep_(outputTimestep) {}

bool DynamicRuptureScheduler::mayComputeInterior(long curCorrectionSteps) const {
  return curCorrectionSteps > lastCorrectionStepsInterior_;
}

void DynamicRuptureScheduler::setLastCorrectionStepsInterior(long steps) {
  lastCorrectionStepsInterior_ = steps;
}

void DynamicRuptureScheduler::setLastCorrectionStepsCopy(long steps) {
  lastCorrectionStepsCopy_ = steps;
}

bool DynamicRuptureScheduler::hasDynamicRuptureFaces() const {
  return numberOfDynamicRuptureFaces_ > 0;
}

double DynamicRuptureScheduler::getOutputTimestep() const { return outputTimestep_; }
} // namespace seissol::time_stepping
