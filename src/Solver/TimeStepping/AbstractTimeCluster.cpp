// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "AbstractTimeCluster.h"

#include "Common/Executor.h"
#include "Solver/TimeStepping/ActorState.h"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <iostream>
#include <memory>
#include <type_traits>
#include <utils/logger.h>
#include <variant>
#include <vector>

namespace seissol::time_stepping {
double AbstractTimeCluster::timeStepSize() const { return ct_.timeStepSize(syncTime_); }

AbstractTimeCluster::AbstractTimeCluster(double maxTimeStepSize,
                                         long timeStepRate,
                                         Executor executor)
    : timeOfLastStageChange_(std::chrono::steady_clock::now()), timeStepRate_(timeStepRate),
      executor_(executor) {
  ct_.maxTimeStepSize = maxTimeStepSize;
  ct_.timeStepRate = timeStepRate;
}

ActorAction AbstractTimeCluster::getNextLegalAction() {
  processMessages();
  switch (state_) {
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
    if (ct_.stepsSinceLastSync == 0) {
      return ActorAction::RestartAfterSync;
    }
    break;
  }
  default:
    logError() << "Invalid actor state in getNextLegalAction()" << static_cast<int>(state_);
  }
  return ActorAction::Nothing;
}

void AbstractTimeCluster::unsafePerformAction(ActorAction action) {
  switch (action) {
  case ActorAction::Nothing:
    break;
  case ActorAction::Correct:
    assert(state_ == ActorState::Predicted);
    correct();
    ct_.correctionTime += timeStepSize();
    ++numberOfTimeSteps_;
    ct_.stepsSinceLastSync += ct_.timeStepRate;
    ct_.stepsSinceStart += ct_.timeStepRate;
    for (auto& neighbor : neighbors_) {
      const bool justBeforeSync = ct_.stepsUntilSync <= ct_.predictionsSinceLastSync;
      const bool sendMessage =
          justBeforeSync || ct_.stepsSinceLastSync >= neighbor.ct.predictionsSinceLastSync;
      if (sendMessage) {
        AdvancedCorrectionTimeMessage message{};
        message.time = ct_.correctionTime;
        message.stepsSinceSync = ct_.stepsSinceLastSync;
        neighbor.outbox->push(message);
      }
    }
    state_ = ActorState::Corrected;
    break;
  case ActorAction::Predict:
    assert(state_ == ActorState::Corrected);
    predict();
    ct_.predictionsSinceLastSync += ct_.timeStepRate;
    ct_.predictionsSinceStart += ct_.timeStepRate;
    ct_.predictionTime += timeStepSize();

    for (auto& neighbor : neighbors_) {
      // Maybe check also how many steps neighbor has to sync!
      const bool justBeforeSync = ct_.stepsUntilSync <= ct_.predictionsSinceLastSync;
      const bool sendMessage =
          justBeforeSync || ct_.predictionsSinceLastSync >= neighbor.ct.nextCorrectionSteps();
      if (sendMessage) {
        AdvancedPredictionTimeMessage message{};
        message.time = ct_.predictionTime;
        message.stepsSinceSync = ct_.predictionsSinceLastSync;
        neighbor.outbox->push(message);
      }
    }
    state_ = ActorState::Predicted;
    break;
  case ActorAction::Sync:
    assert(state_ == ActorState::Corrected);
    logDebug() << "synced at" << syncTime_ << ", corrTime =" << ct_.correctionTime
               << "stepsSinceLastSync" << ct_.stepsSinceLastSync << "stepsUntilLastSync"
               << ct_.stepsUntilSync << std::endl;
    state_ = ActorState::Synced;
    break;
  case ActorAction::RestartAfterSync:
    start();
    state_ = ActorState::Corrected;
    break;
  default:
    logError() << "Invalid actor action in getNextLegalAction()" << static_cast<int>(state_);
    break;
  }
}

ActResult AbstractTimeCluster::act() {
  ActResult result;
  auto stateBefore = state_;
  auto nextAction = getNextLegalAction();
  unsafePerformAction(nextAction);

  const auto currentTime = std::chrono::steady_clock::now();
  result.isStateChanged = stateBefore != state_;
  if (!result.isStateChanged) {
    const auto timeSinceLastUpdate = currentTime - timeOfLastStageChange_;
    if (timeSinceLastUpdate > timeout && !alreadyPrintedTimeOut_) {
      alreadyPrintedTimeOut_ = true;
      printTimeoutMessage(std::chrono::duration_cast<std::chrono::seconds>(timeSinceLastUpdate));
    }
  } else {
    timeOfLastStageChange_ = currentTime;
    alreadyPrintedTimeOut_ = false;
  }
  return result;
}

bool AbstractTimeCluster::processMessages() {
  bool processed = false;
  for (auto& neighbor : neighbors_) {
    if (neighbor.inbox->hasMessages()) {
      processed = true;
      Message message = neighbor.inbox->pop();
      std::visit(
          [&neighbor, this](auto&& msg) {
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
              static_assert(sizeof(T) == 0, "non-exhaustive visitor!");
            }
          },
          message);
    }
  }
  return processed;
}

bool AbstractTimeCluster::mayPredict() {
  // We can predict, if our prediction time is smaller/equals than the next correction time of all
  // neighbors.
  const auto minNeighborSteps = std::min_element(
      neighbors_.begin(), neighbors_.end(), [](const NeighborCluster& a, const NeighborCluster& b) {
        return a.ct.nextCorrectionSteps() < b.ct.nextCorrectionSteps();
      });
  const bool stepBasedPredict =
      minNeighborSteps == neighbors_.end() ||
      ct_.predictionsSinceLastSync < minNeighborSteps->ct.nextCorrectionSteps();
  return stepBasedPredict;
}

bool AbstractTimeCluster::mayCorrect() {
  // We can correct, if our prediction time is smaller than the one of all neighbors.
  bool stepBasedCorrect = true;
  for (auto& neighbor : neighbors_) {
    const bool isSynced = neighbor.ct.stepsUntilSync <= neighbor.ct.predictionsSinceLastSync;
    stepBasedCorrect =
        stepBasedCorrect &&
        (isSynced || (ct_.predictionsSinceLastSync <= neighbor.ct.predictionsSinceLastSync));
  }
  return stepBasedCorrect;
}

Executor AbstractTimeCluster::getExecutor() const { return executor_; }

bool AbstractTimeCluster::maySync() { return ct_.stepsSinceLastSync >= ct_.stepsUntilSync; }

void AbstractTimeCluster::connect(AbstractTimeCluster& other) {
  neighbors_.emplace_back(other.ct_.maxTimeStepSize, other.ct_.timeStepRate, other.executor_);
  other.neighbors_.emplace_back(ct_.maxTimeStepSize, ct_.timeStepRate, executor_);
  neighbors_.back().inbox = std::make_shared<MessageQueue>();
  other.neighbors_.back().inbox = std::make_shared<MessageQueue>();
  neighbors_.back().outbox = other.neighbors_.back().inbox;
  other.neighbors_.back().outbox = neighbors_.back().inbox;
}

void AbstractTimeCluster::setSyncTime(double newSyncTime) {
  assert(newSyncTime > syncTime_);
  assert(state_ == ActorState::Synced);
  syncTime_ = newSyncTime;
}

bool AbstractTimeCluster::synced() const { return state_ == ActorState::Synced; }
void AbstractTimeCluster::reset() {
  assert(state_ == ActorState::Synced);

  // There can be pending messages from before the sync point
  processMessages();
  for ([[maybe_unused]] const auto& neighbor : neighbors_) {
    assert(!neighbor.inbox->hasMessages());
  }
  ct_.stepsSinceLastSync = 0;
  ct_.predictionsSinceLastSync = 0;
  ct_.stepsUntilSync = ct_.computeStepsUntilSyncTime(ct_.correctionTime, syncTime_);

  for (auto& neighbor : neighbors_) {
    neighbor.ct.stepsUntilSync =
        neighbor.ct.computeStepsUntilSyncTime(ct_.correctionTime, syncTime_);
    neighbor.ct.stepsSinceLastSync = 0;
    neighbor.ct.predictionsSinceLastSync = 0;
  }
}

ActorPriority AbstractTimeCluster::getPriority() const { return priority_; }

void AbstractTimeCluster::setPriority(ActorPriority newPriority) { this->priority_ = newPriority; }

ActorState AbstractTimeCluster::getState() const { return state_; }

void AbstractTimeCluster::setTime(double time) {
  ct_.predictionTime = time;
  ct_.correctionTime = time;
  for (auto& neighbor : neighbors_) {
    neighbor.ct.predictionTime = time;
    neighbor.ct.correctionTime = time;
  }
}

long AbstractTimeCluster::getTimeStepRate() const { return timeStepRate_; }

void AbstractTimeCluster::finalize() {}

double AbstractTimeCluster::getClusterTimes() { return ct_.getTimeStepSize(); }

void AbstractTimeCluster::setClusterTimes(double newTimeStepSize) {
  ct_.setTimeStepSize(newTimeStepSize);
}

std::vector<NeighborCluster>* AbstractTimeCluster::getNeighborClusters() { return &neighbors_; }

bool AbstractTimeCluster::hasDifferentExecutorNeighbor() {
  return std::any_of(neighbors_.begin(), neighbors_.end(), [&](auto& neighbor) {
    return neighbor.executor != executor_;
  });
}

void AbstractTimeCluster::finishPhase() {}

} // namespace seissol::time_stepping
