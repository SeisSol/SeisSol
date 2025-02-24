// SPDX-FileCopyrightText: 2020-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <algorithm>
#include <iostream>
#include <cassert>
#include "utils/logger.h"

#include "Parallel/MPI.h"
#include "AbstractTimeCluster.h"

namespace seissol::time_stepping {
double AbstractTimeCluster::timeStepSize() const {
  return ct.timeStepSize(syncTime);
}

AbstractTimeCluster::AbstractTimeCluster(double maxTimeStepSize, long timeStepRate, Executor executor)
    : timeOfLastStageChange(std::chrono::steady_clock::now()),
      timeStepRate(timeStepRate), numberOfTimeSteps(0), executor(executor) {
  ct.maxTimeStepSize = maxTimeStepSize;
  ct.timeStepRate = timeStepRate;
}

ActorAction AbstractTimeCluster::getNextLegalAction() {
  processMessages();
  switch (state) {
    case ActorState::Corrected: {
      if (maySync()) {
        return ActorAction::Sync;
      } else if (mayPredict()) {
        return ActorAction::Predict;
      }
      break;
    }
    case ActorState::Predicted: {
      if (mayCorrect()) {
        return ActorAction::Correct;
      }
      break;
    }
    case ActorState::Synced: {
      if (ct.stepsSinceLastSync == 0) {
        return ActorAction::RestartAfterSync;
      }
      break;
    }
    default:
      logError() << "Invalid actor state in getNextLegalAction()" << static_cast<int>(state);
  }
  return ActorAction::Nothing;
}

void AbstractTimeCluster::unsafePerformAction(ActorAction action) {
  switch (action) {
    case ActorAction::Nothing:
      break;
    case ActorAction::Correct:
      assert(state == ActorState::Predicted);
      correct();
      ct.correctionTime += timeStepSize();
      ++numberOfTimeSteps;
      ct.stepsSinceLastSync += ct.timeStepRate;
      ct.stepsSinceStart += ct.timeStepRate;
      for (auto &neighbor : neighbors) {
        const bool justBeforeSync = ct.stepsUntilSync <= ct.predictionsSinceLastSync;
        const bool sendMessage = justBeforeSync
                                 || ct.stepsSinceLastSync >= neighbor.ct.predictionsSinceLastSync;
        if (sendMessage) {
          AdvancedCorrectionTimeMessage message{};
          message.time = ct.correctionTime;
          message.stepsSinceSync = ct.stepsSinceLastSync;
          neighbor.outbox->push(message);
        }
      }
      state = ActorState::Corrected;
      break;
    case ActorAction::Predict:
      assert(state == ActorState::Corrected);
      predict();
      ct.predictionsSinceLastSync += ct.timeStepRate;
      ct.predictionsSinceStart += ct.timeStepRate;
      ct.predictionTime += timeStepSize();

      for (auto &neighbor : neighbors) {
        // Maybe check also how many steps neighbor has to sync!
        const bool justBeforeSync = ct.stepsUntilSync <= ct.predictionsSinceLastSync;
        const bool sendMessage = justBeforeSync
                                 || ct.predictionsSinceLastSync >= neighbor.ct.nextCorrectionSteps();
        if (sendMessage) {
          AdvancedPredictionTimeMessage message{};
          message.time = ct.predictionTime;
          message.stepsSinceSync = ct.predictionsSinceLastSync;
          neighbor.outbox->push(message);
        }
      }
      state = ActorState::Predicted;
      break;
    case ActorAction::Sync:
      assert(state == ActorState::Corrected);
      logDebug() << "synced at" << syncTime
                                << ", corrTime =" << ct.correctionTime
                                << "stepsSinceLastSync" << ct.stepsSinceLastSync
                                << "stepsUntilLastSync" << ct.stepsUntilSync
                                << std::endl;
      state = ActorState::Synced;
      break;
    case ActorAction::RestartAfterSync:
      start();
      state = ActorState::Corrected;
      break;
    default:
      logError() << "Invalid actor action in getNextLegalAction()" << static_cast<int>(state);
      break;
  }
}

ActResult AbstractTimeCluster::act() {
  ActResult result;
  auto stateBefore = state;
  auto nextAction = getNextLegalAction();
  unsafePerformAction(nextAction);

  const auto currentTime = std::chrono::steady_clock::now();
  result.isStateChanged = stateBefore != state;
  if (!result.isStateChanged) {
    const auto timeSinceLastUpdate = currentTime - timeOfLastStageChange;
    if (timeSinceLastUpdate > timeout && !alreadyPrintedTimeOut) {
        alreadyPrintedTimeOut = true;
        printTimeoutMessage(std::chrono::duration_cast<std::chrono::seconds>(timeSinceLastUpdate));
    }
  } else {
    timeOfLastStageChange = currentTime;
    alreadyPrintedTimeOut = false;
  }
  return result;
}


bool AbstractTimeCluster::processMessages() {
  bool processed = false;
  for (auto& neighbor : neighbors) {
    if (neighbor.inbox->hasMessages()) {
      processed = true;
      Message message = neighbor.inbox->pop();
      std::visit([&neighbor, this](auto&& msg) {
        using T = std::decay_t<decltype(msg)>;
        if constexpr (std::is_same_v<T, AdvancedPredictionTimeMessage>) {
          assert(msg.time > neighbor.ct.predictionTime);
          neighbor.ct.predictionTime = msg.time;
          neighbor.ct.predictionsSinceLastSync = msg.stepsSinceSync;
          handleAdvancedPredictionTimeMessage(neighbor);
        } else if constexpr (std::is_same_v<T, AdvancedCorrectionTimeMessage>) {
          assert(msg.time > neighbor.ct.correctionTime);
          neighbor.ct.correctionTime = msg.time;
          neighbor.ct.stepsSinceLastSync = msg.stepsSinceSync;
          handleAdvancedCorrectionTimeMessage(neighbor);
        } else {
          static_assert(always_false<T>::value, "non-exhaustive visitor!");
        }
      }, message);
    }
  }
  return processed;
}

bool AbstractTimeCluster::mayPredict() {
  // We can predict, if our prediction time is smaller/equals than the next correction time of all neighbors.
    const auto minNeighborSteps = std::min_element(
            neighbors.begin(), neighbors.end(),
            [](NeighborCluster const &a, NeighborCluster const &b) {
                return a.ct.nextCorrectionSteps() < b.ct.nextCorrectionSteps();
            });
    bool stepBasedPredict = minNeighborSteps == neighbors.end()
            || ct.predictionsSinceLastSync < minNeighborSteps->ct.nextCorrectionSteps();
    return stepBasedPredict;
}

bool AbstractTimeCluster::mayCorrect() {
  // We can correct, if our prediction time is smaller than the one of all neighbors.
  bool stepBasedCorrect = true;
  for (auto& neighbor : neighbors) {
      const bool isSynced = neighbor.ct.stepsUntilSync <= neighbor.ct.predictionsSinceLastSync;
      stepBasedCorrect = stepBasedCorrect
              && (isSynced || (ct.predictionsSinceLastSync <= neighbor.ct.predictionsSinceLastSync));
  }
  return stepBasedCorrect;
}

Executor AbstractTimeCluster::getExecutor() const {
  return executor;
}

bool AbstractTimeCluster::maySync() {
    return ct.stepsSinceLastSync >= ct.stepsUntilSync;
}

void AbstractTimeCluster::connect(AbstractTimeCluster &other) {
  neighbors.emplace_back(other.ct.maxTimeStepSize, other.ct.timeStepRate, other.executor);
  other.neighbors.emplace_back(ct.maxTimeStepSize, ct.timeStepRate, executor);
  neighbors.back().inbox = std::make_shared<MessageQueue>();
  other.neighbors.back().inbox = std::make_shared<MessageQueue>();
  neighbors.back().outbox = other.neighbors.back().inbox;
  other.neighbors.back().outbox = neighbors.back().inbox;
}

void AbstractTimeCluster::setSyncTime(double newSyncTime) {
  assert(newSyncTime > syncTime);
  assert(state == ActorState::Synced);
  syncTime = newSyncTime;
}

bool AbstractTimeCluster::synced() const {
  return state == ActorState::Synced;
}
void AbstractTimeCluster::reset() {
  assert(state == ActorState::Synced);

  // There can be pending messages from before the sync point
  processMessages();
  for ([[maybe_unused]] const auto& neighbor : neighbors) {
    assert(!neighbor.inbox->hasMessages());
  }
  ct.stepsSinceLastSync = 0;
  ct.predictionsSinceLastSync = 0;
  ct.stepsUntilSync = ct.computeStepsUntilSyncTime(ct.correctionTime, syncTime);

  for (auto& neighbor : neighbors) {
    neighbor.ct.stepsUntilSync = neighbor.ct.computeStepsUntilSyncTime(ct.correctionTime, syncTime);
    neighbor.ct.stepsSinceLastSync = 0;
    neighbor.ct.predictionsSinceLastSync = 0;
  }

}

ActorPriority AbstractTimeCluster::getPriority() const {
  return priority;
}

void AbstractTimeCluster::setPriority(ActorPriority newPriority) {
  this->priority = newPriority;
}

ActorState AbstractTimeCluster::getState() const {
  return state;
}

void AbstractTimeCluster::setPredictionTime(double time) {
  ct.predictionTime = time;
  for (auto& neighbor : neighbors) {
    neighbor.ct.predictionTime = time;
  }
}

void AbstractTimeCluster::setCorrectionTime(double time) {
  ct.correctionTime = time;
  for (auto& neighbor : neighbors) {
    neighbor.ct.correctionTime = time;
  }
}

long AbstractTimeCluster::getTimeStepRate() {
  return timeStepRate;
}

void AbstractTimeCluster::finalize() {}

double AbstractTimeCluster::getClusterTimes(){
  return ct.getTimeStepSize();
}

void AbstractTimeCluster::setClusterTimes(double newTimeStepSize) {
  ct.setTimeStepSize(newTimeStepSize);
}

std::vector<NeighborCluster>* AbstractTimeCluster::getNeighborClusters(){
  return &neighbors;
}

bool AbstractTimeCluster::hasDifferentExecutorNeighbor() {
  return std::any_of(neighbors.begin(), neighbors.end(), [&](auto& neighbor) {
    return neighbor.executor != executor;
  });
}

} // namespace seissol::time_stepping

