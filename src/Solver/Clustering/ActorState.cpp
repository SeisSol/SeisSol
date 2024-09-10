#include "ActorState.h"

#include <Common/Executor.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <mutex>
#include <string>

namespace seissol::solver::clustering {

inline std::ostream& operator<<(std::ostream& stream, const Message& message) {
  stream << "Message, t = " << message.time;
  return stream;
}

std::string actorStateToString(ActorState state) {
  const std::string stepstring = [&]() {
    if (state.step == ComputeStep::Interact) {
      return "Interact";
    }
    if (state.step == ComputeStep::Predict) {
      return "Predict";
    }
    if (state.step == ComputeStep::Correct) {
      return "Correct";
    }
    throw;
  }();
  if (state.type == StateType::Synchronized) {
    return "Synchronized";
  }
  if (state.type == StateType::ComputeStart) {
    return "ComputeStart: " + stepstring;
  }
  if (state.type == StateType::ComputeDone) {
    return "ComputeDone: " + stepstring;
  }
  throw;
}

void MessageQueue::push(const Message& message) {
  const std::lock_guard lock{mutex};
  queue.push(message);
}

Message MessageQueue::pop() {
  const std::lock_guard lock{mutex};
  const Message message = queue.front();
  queue.pop();
  return message;
}

bool MessageQueue::hasMessages() const { return !queue.empty(); }

size_t MessageQueue::size() const { return queue.size(); }

double ClusterTimes::nextComputeTime(ComputeStep step, double syncTime) const {
  return std::min(syncTime, time.at(step) + maxTimeStepSize);
}

long ClusterTimes::nextSteps() const {
  return std::min(computeSinceLastSync.at(ComputeStep::Correct) + timeStepRate, stepsUntilSync);
}

double ClusterTimes::timeStepSize(double syncTime) const {
  return std::min(syncTime - time.at(ComputeStep::Correct), maxTimeStepSize);
}

long ClusterTimes::computeStepsUntilSyncTime(double oldSyncTime, double newSyncTime) const {
  const double timeDiff = newSyncTime - oldSyncTime;
  return static_cast<long>(std::ceil(timeStepRate * timeDiff / maxTimeStepSize));
}

NeighborCluster::NeighborCluster(double maxTimeStepSize,
                                 int timeStepRate,
                                 Executor neighborExecutor) {
  ct.maxTimeStepSize = maxTimeStepSize;
  ct.timeStepRate = timeStepRate;
  executor = neighborExecutor;
}

} // namespace seissol::time_stepping
