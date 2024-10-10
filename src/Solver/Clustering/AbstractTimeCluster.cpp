// SPDX-FileCopyrightText: 2020-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#include "Solver/Clustering/AbstractTimeCluster.h"

#include "ActorState.h"
#include "utils/logger.h"
#include <Common/Executor.h>
#include <Parallel/Helper.h>
#include <Parallel/Host/CpuExecutor.h>
#include <Parallel/Runtime/Stream.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <iostream>
#include <memory>
#include <vector>

#include "Parallel/MPI.h"

namespace seissol::solver::clustering {
double AbstractTimeCluster::timeStepSize() const { return ct.timeStepSize(syncTime); }

AbstractTimeCluster::AbstractTimeCluster(double maxTimeStepSize,
                                         long timeStepRate,
                                         Executor executor,
                                         const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor)
    : timeOfLastStageChange(std::chrono::steady_clock::now()), timeStepRate(timeStepRate),
      numberOfTimeSteps(0), executor(executor), streamRuntime(cpuExecutor) {
  ct.maxTimeStepSize = maxTimeStepSize;
  ct.timeStepRate = timeStepRate;

  concurrent = concurrentClusters();
}

bool AbstractTimeCluster::advanceState() {
  processMessages();
  if (state.type == StateType::Synchronized) {
    if (ct.computeSinceLastSync[lastStep()] == 0) {
      restart();
      return true;
    }
  } else if (state.type == StateType::ComputeDone) {
    auto next = nextStep(state.step);
    if (state.step == lastStep() && maySynchronize()) {
      synchronize();
      return true;
    } else if (allComputed(state.step)) {
      preCompute(state.step);
      runCompute(next);
      if (!concurrent) {
        streamRuntime.wait();
      }
      state.step = next;
      state.type = StateType::ComputeStart;
      return true;
    } else if (emptyStep(next)) {
      state.step = next;
      state.type = StateType::ComputeStart;
      return true;
    }
  } else if (state.type == StateType::ComputeStart) {
    if (emptyStep(state.step) || pollCompute(state.step)) {
      postCompute(state.step);
      state.type = StateType::ComputeDone;
      return true;
    }
  }
  return false;
}

void AbstractTimeCluster::synchronize() {
  assert(state.type == StateType::ComputeDone);
  logDebug(MPI::mpi.rank()) << "synced at" << syncTime
                            << ", corrTime =" << ct.time[ComputeStep::Correct]
                            << "computeSinceLastSync.at(ComputeStep::Correct)"
                            << ct.computeSinceLastSync.at(ComputeStep::Correct)
                            << "stepsUntilLastSync" << ct.stepsUntilSync << std::endl;
  state.type = StateType::Synchronized;
}

void AbstractTimeCluster::restart() {
  start();
  state.type = StateType::ComputeDone;
  state.step = lastStep();
}

AbstractTimeCluster::~AbstractTimeCluster() = default;

void AbstractTimeCluster::postCompute(ComputeStep step) {
  ct.computeSinceLastSync[step] += ct.timeStepRate;
  ct.computeSinceStart[step] += ct.timeStepRate;
  ct.time[step] += timeStepSize();

  if (!emptyStep(step) && concurrent) {
    events.push(streamRuntime.recordEvent());
  }

  for (auto& neighbor : neighbors) {
    // TODO: maybe check also how many steps the neighbor has to sync?
    const bool justBeforeSync = ct.stepsUntilSync <= ct.computeSinceLastSync[step];

    // for now, keep manual
    const bool usefulUpdate = (ct.computeSinceLastSync.at(step) >= neighbor.ct.nextSteps() &&
                               step == ComputeStep::Predict) ||
                              (ct.computeSinceLastSync.at(step) >=
                                   neighbor.ct.computeSinceLastSync.at(ComputeStep::Predict) &&
                               step == ComputeStep::Interact) ||
                              (ct.computeSinceLastSync.at(step) >=
                                   neighbor.ct.computeSinceLastSync.at(ComputeStep::Interact) &&
                               step == ComputeStep::Correct);

    const bool sendMessage = justBeforeSync || usefulUpdate;
    if (sendMessage) {
      Message message{};
      message.step = step;
      message.time = ct.time[step];
      message.stepsSinceSync = ct.computeSinceLastSync[step];
      if (events.empty()) {
        message.completionEvent = nullptr;
      } else {
        message.completionEvent = events.back();
      }
      neighbor.outbox->push(message);
    }
  }
}

ActResult AbstractTimeCluster::act() {
  ActResult result;
  auto stateBefore = state;
  auto changed = advanceState();

  const auto currentTime = std::chrono::steady_clock::now();
  result.isStateChanged = changed;
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
      const Message message = neighbor.inbox->pop();
      assert(message.time > neighbor.ct.time[message.step]);
      neighbor.ct.time[message.step] = message.time;
      neighbor.ct.computeSinceLastSync[message.step] = message.stepsSinceSync;
      neighbor.events[message.step] = message.completionEvent;
      handleAdvancedComputeTimeMessage(message.step, neighbor);
    }
  }
  return processed;
}

bool AbstractTimeCluster::allComputed(ComputeStep step) {
  // special case for the last step: we need to look into the future
  if (step == lastStep()) {
    const auto minNeighborSteps = std::min_element(
        neighbors.begin(), neighbors.end(), [](const NeighborCluster& a, const NeighborCluster& b) {
          return a.ct.nextSteps() < b.ct.nextSteps();
        });
    const bool stepBasedPredict =
        minNeighborSteps == neighbors.end() ||
        ct.computeSinceLastSync.at(step) < minNeighborSteps->ct.nextSteps();
    return stepBasedPredict;
  } else {
    bool stepBasedCorrect = true;
    for (auto& neighbor : neighbors) {
      const bool isSynced = neighbor.ct.stepsUntilSync <= neighbor.ct.computeSinceLastSync.at(step);
      const bool isAdvanced =
          ct.computeSinceLastSync.at(step) <= neighbor.ct.computeSinceLastSync.at(step);
      stepBasedCorrect = stepBasedCorrect && (isSynced || isAdvanced);
    }
    return stepBasedCorrect;
  }
}

void AbstractTimeCluster::preCompute(ComputeStep step) {
  if (concurrent) {
    for (auto& neighbor : neighbors) {
      // (only) wait upon all events that we haven't waited upon yet
      auto eventfind = neighbor.events.find(step);
      if (eventfind != neighbor.events.end() && neighbor.events.at(step) != nullptr) {
        streamRuntime.waitEvent(eventfind->second);

        // forget about events that have been waited upon already
        eventfind->second = nullptr;
      }
    }
  }
}

Executor AbstractTimeCluster::getExecutor() const { return executor; }

bool AbstractTimeCluster::maySynchronize() {
  return ct.computeSinceLastSync.at(lastStep()) >= ct.stepsUntilSync;
}

void AbstractTimeCluster::connect(AbstractTimeCluster& other) {
  neighbors.emplace_back(other.ct.maxTimeStepSize, other.ct.timeStepRate, other.executor);
  other.neighbors.emplace_back(ct.maxTimeStepSize, ct.timeStepRate, executor);
  neighbors.back().inbox = std::make_shared<MessageQueue>();
  other.neighbors.back().inbox = std::make_shared<MessageQueue>();
  neighbors.back().outbox = other.neighbors.back().inbox;
  other.neighbors.back().outbox = neighbors.back().inbox;
}

void AbstractTimeCluster::setSyncTime(double newSyncTime) {
  assert(newSyncTime > syncTime);
  assert(state.type == StateType::Synchronized);
  syncTime = newSyncTime;
}

bool AbstractTimeCluster::synchronized() const { return state.type == StateType::Synchronized; }
void AbstractTimeCluster::reset() {
  assert(state.type == StateType::Synchronized);

  // There can be pending messages from before the sync point
  processMessages();
  for (auto& neighbor : neighbors) {
    assert(!neighbor.inbox->hasMessages());
  }
  ct.computeSinceLastSync.clear();
  ct.computeSinceLastSync[ComputeStep::Predict] = 0;
  ct.computeSinceLastSync[ComputeStep::Interact] = 0;
  ct.computeSinceLastSync[ComputeStep::Correct] = 0;
  ct.stepsUntilSync = ct.computeStepsUntilSyncTime(ct.time[lastStep()], syncTime);

  for (auto& neighbor : neighbors) {
    neighbor.ct.stepsUntilSync =
        neighbor.ct.computeStepsUntilSyncTime(ct.time[lastStep()], syncTime);
    neighbor.ct.computeSinceLastSync.clear();
    neighbor.ct.computeSinceLastSync[ComputeStep::Predict] = 0;
    neighbor.ct.computeSinceLastSync[ComputeStep::Interact] = 0;
    neighbor.ct.computeSinceLastSync[ComputeStep::Correct] = 0;
    neighbor.events.clear();
  }

  while (!events.empty()) {
    streamRuntime.recycleEvent(events.front());
    events.pop();
  }
}

ActorPriority AbstractTimeCluster::getPriority() const { return priority; }

void AbstractTimeCluster::setPriority(ActorPriority newPriority) { this->priority = newPriority; }

ActorState AbstractTimeCluster::getState() const { return state; }

void AbstractTimeCluster::setTime(double time) {
  ct.time[ComputeStep::Predict] = time;
  ct.time[ComputeStep::Interact] = time;
  ct.time[ComputeStep::Correct] = time;
}

long AbstractTimeCluster::getTimeStepRate() { return timeStepRate; }

void AbstractTimeCluster::finalize() { streamRuntime.dispose(); }

double AbstractTimeCluster::getClusterTimes() { return ct.getTimeStepSize(); }

void AbstractTimeCluster::setClusterTimes(double newTimeStepSize) {
  ct.setTimeStepSize(newTimeStepSize);
}

std::vector<NeighborCluster>* AbstractTimeCluster::getNeighborClusters() { return &neighbors; }

bool AbstractTimeCluster::hasDifferentExecutorNeighbor() {
  return std::any_of(neighbors.begin(), neighbors.end(), [&](auto& neighbor) {
    return neighbor.executor != executor;
  });
}

CellCluster::CellCluster(double maxTimeStepSize, long timeStepRate, Executor executor, const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor)
    : AbstractTimeCluster(maxTimeStepSize, timeStepRate, executor, cpuExecutor) {}

FaceCluster::FaceCluster(double maxTimeStepSize, long timeStepRate, Executor executor, const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor)
    : AbstractTimeCluster(maxTimeStepSize, timeStepRate, executor, cpuExecutor) {}

CellCluster::~CellCluster() = default;

FaceCluster::~FaceCluster() = default;

} // namespace seissol::solver::clustering
