// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "AbstractTimeCluster.h"

#include "Common/Executor.h"
#include "Solver/TimeStepping/ActorState.h"
#include "Solver/TimeStepping/StepParams.h"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <iostream>
#include <utils/logger.h>
#include <vector>

namespace seissol::time_stepping {
double AbstractTimeCluster::timeStepSize() const { return ct_.timeStepSize(syncTime_); }

StepParams AbstractTimeCluster::stepParams() const {
  return computeStepParams(ct_, syncTime_, stepContext_);
}

AbstractTimeCluster::AbstractTimeCluster(double maxTimeStepSize,
                                         long timeStepRate,
                                         Executor executor)
    : timeOfLastStageChange_(std::chrono::steady_clock::now()), timeStepRate_(timeStepRate),
      executor_(executor) {
  ct_.maxTimeStepSize = maxTimeStepSize;
  ct_.timeStepRate = timeStepRate;
}

ActorAction AbstractTimeCluster::getNextLegalAction() {
  refreshNeighbors();
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
    publishProgress();
    state_ = ActorState::Corrected;
    break;
  case ActorAction::Predict:
    assert(state_ == ActorState::Corrected);
    predict();
    ct_.predictionsSinceLastSync += ct_.timeStepRate;
    ct_.predictionsSinceStart += ct_.timeStepRate;
    ct_.predictionTime += timeStepSize();
    publishProgress();
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

  result.isStateChanged = stateBefore != state_;
  trackProgress(result.isStateChanged);
  return result;
}

void AbstractTimeCluster::trackProgress(bool progressed) {
  const auto currentTime = std::chrono::steady_clock::now();
  if (!progressed) {
    const auto timeSinceLastUpdate = currentTime - timeOfLastStageChange_;
    if (timeSinceLastUpdate > timeout && !alreadyPrintedTimeOut_) {
      alreadyPrintedTimeOut_ = true;
      printTimeoutMessage(std::chrono::duration_cast<std::chrono::seconds>(timeSinceLastUpdate));
    }
  } else {
    timeOfLastStageChange_ = currentTime;
    alreadyPrintedTimeOut_ = false;
  }
}

void AbstractTimeCluster::publishProgress() {
  progress_.predictionTime.store(ct_.predictionTime, std::memory_order_relaxed);
  progress_.correctionTime.store(ct_.correctionTime, std::memory_order_relaxed);
  progress_.stepsUntilSync.store(ct_.stepsUntilSync, std::memory_order_relaxed);
  progress_.stepsSinceLastSync.store(ct_.stepsSinceLastSync, std::memory_order_release);
  progress_.predictionsSinceLastSync.store(ct_.predictionsSinceLastSync, std::memory_order_release);
}

void AbstractTimeCluster::refreshNeighbors() {
  for (auto& neighbor : neighbors_) {
    const auto predictions =
        neighbor.progress->predictionsSinceLastSync.load(std::memory_order_acquire);
    const auto corrections = neighbor.progress->stepsSinceLastSync.load(std::memory_order_acquire);

    // a cluster corrects a step only after it has predicted it; the correction thus concerns an
    // older step and is handled first
    if (corrections > neighbor.ct.stepsSinceLastSync) {
      neighbor.ct.stepsSinceLastSync = corrections;
      neighbor.ct.correctionTime =
          neighbor.progress->correctionTime.load(std::memory_order_relaxed);
      handleNeighborCorrection(neighbor);
    }
    if (predictions > neighbor.ct.predictionsSinceLastSync) {
      neighbor.ct.predictionsSinceLastSync = predictions;
      neighbor.ct.predictionTime =
          neighbor.progress->predictionTime.load(std::memory_order_relaxed);
      handleNeighborPrediction(neighbor);
    }
  }
}

bool AbstractTimeCluster::mayPredict() {
  // We can predict, if our prediction time is smaller/equals than the next correction time of all
  // neighbors.
  bool stepBasedPredict = true;
  for (const auto& neighbor : neighbors_) {
    stepBasedPredict =
        stepBasedPredict && ct_.predictionsSinceLastSync < neighbor.ct.nextCorrectionSteps();
  }
  return stepBasedPredict;
}

bool AbstractTimeCluster::mayCorrect() {
  // We can correct, if our prediction time is smaller than the one of all neighbors.
  bool stepBasedCorrect = true;
  for (auto& neighbor : neighbors_) {
    // the progress up to which the neighbor has made its data available
    const auto provided = neighbor.dataReadiness == DataReadiness::AfterPrediction
                              ? neighbor.ct.predictionsSinceLastSync
                              : neighbor.ct.stepsSinceLastSync;
    const bool isSynced = neighbor.ct.stepsUntilSync <= provided;
    stepBasedCorrect = stepBasedCorrect && (isSynced || (ct_.predictionsSinceLastSync <= provided));
  }
  return stepBasedCorrect;
}

Executor AbstractTimeCluster::getExecutor() const { return executor_; }

bool AbstractTimeCluster::maySync() { return ct_.stepsSinceLastSync >= ct_.stepsUntilSync; }

void AbstractTimeCluster::connect(AbstractTimeCluster& other) {
  observe(other);
  other.observe(*this);
}

void AbstractTimeCluster::observe(AbstractTimeCluster& other) {
  auto& neighbor =
      neighbors_.emplace_back(other.ct_.maxTimeStepSize, other.ct_.timeStepRate, other.executor_);
  neighbor.progress = &other.progress_;
  neighbor.dataReadiness = other.dataReadiness();
}

DataReadiness AbstractTimeCluster::dataReadiness() const { return DataReadiness::AfterPrediction; }

void AbstractTimeCluster::setSyncTime(double newSyncTime) {
  assert(newSyncTime > syncTime_);
  assert(state_ == ActorState::Synced);
  syncTime_ = newSyncTime;
}

bool AbstractTimeCluster::synced() const { return state_ == ActorState::Synced; }
void AbstractTimeCluster::reset() {
  assert(state_ == ActorState::Synced);

  ct_.stepsSinceLastSync = 0;
  ct_.predictionsSinceLastSync = 0;
  ct_.stepsUntilSync = ct_.computeStepsUntilSyncTime(ct_.correctionTime, syncTime_);

  for (auto& neighbor : neighbors_) {
    neighbor.ct.stepsUntilSync =
        neighbor.ct.computeStepsUntilSyncTime(ct_.correctionTime, syncTime_);
    neighbor.ct.stepsSinceLastSync = 0;
    neighbor.ct.predictionsSinceLastSync = 0;
  }

  stepContext_ = computeStepContext(ct_, neighbors_);

  publishProgress();
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
  publishProgress();
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

void AbstractTimeCluster::printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) {
  logWarning(true) << "No update since " << timeSinceLastUpdate.count() << "[s] for cluster "
                   << description() << " at state " << actorStateToString(state_)
                   << " mayPredict = " << mayPredict()
                   << " mayPredict (steps) = " << AbstractTimeCluster::mayPredict()
                   << " mayCorrect = " << mayCorrect()
                   << " mayCorrect (steps) = " << AbstractTimeCluster::mayCorrect()
                   << " maySync = " << maySync();
  for (auto& neighbor : neighbors_) {
    logWarning(true) << "Neighbor with rate = " << neighbor.ct.timeStepRate
                     << "PredTime = " << neighbor.ct.predictionTime
                     << "CorrTime = " << neighbor.ct.correctionTime
                     << "predictionsSinceSync = " << neighbor.ct.predictionsSinceLastSync
                     << "correctionsSinceSync = " << neighbor.ct.stepsSinceLastSync;
  }
  if (timeoutFail()) {
    logError() << "Cluster" << description() << "timed out. Aborting simulation.";
  }
}

bool AbstractTimeCluster::timeoutFail() const { return false; }

} // namespace seissol::time_stepping
