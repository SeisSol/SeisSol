#include <iostream>
#include <Parallel/MPI.h>
#include "Solver/time_stepping/AbstractTimeCluster.h"

namespace seissol::time_stepping {
double AbstractTimeCluster::timeStepSize() const {
  return ct.timeStepSize(syncTime);
}

AbstractTimeCluster::AbstractTimeCluster(double maxTimeStepSize, double timeTolerance, long timeStepRate)
    : timeTolerance(timeTolerance), timeStepRate(timeStepRate), numberOfTimeSteps(0),
    lastStateChange(std::chrono::steady_clock::now()) {
  ct.maxTimeStepSize = maxTimeStepSize;
  ct.timeStepRate = timeStepRate;
}

ActResult AbstractTimeCluster::act() {
  ActResult result;
  auto stateBefore = state;
  switch (state) {
  case ActorState::Corrected: {
    if (maySync()) {
        logDebug(MPI::mpi.rank()) << "synced at " << syncTime
        << ", corrTIme =" << ct.correctionTime
        << ", time tolerence " << timeTolerance
        << " stepsSinceLastSync " << ct.stepsSinceLastSync
        << " stepsUntilLastSync " << ct.stepsSinceLastSync
        << std::endl;
      state = ActorState::Synced;
    } else if (mayPredict()) {
      predict();
      ct.predictionsSinceLastSync += ct.timeStepRate;
      ct.predictionsSinceStart += ct.timeStepRate;
      ct.predictionTime += timeStepSize();

      for (auto &neighbor : neighbors) {
          // TODO(Lukas) Does this handle sync points correctly?
          // Maybe check also how many steps neighbor has to sync!
          const bool justBeforeSync = ct.stepsUntilSync <= ct.predictionsSinceLastSync;
          const bool sendMessageSteps = justBeforeSync
                  || ct.predictionsSinceLastSync >= neighbor.ct.nextCorrectionSteps();
        if (sendMessageSteps) {
          AdvancedPredictionTimeMessage message{};
          message.time = ct.predictionTime;
          message.stepsSinceSync = ct.predictionsSinceLastSync;
          neighbor.outbox->push(message);
        }
      }
      state = ActorState::Predicted;
    } else {
      result.yield = !processMessages();
    }
    break;
  }
  case ActorState::Predicted: {
    if (mayCorrect()) {
      correct();
      ct.correctionTime += timeStepSize();
      ++numberOfTimeSteps;
      ct.stepsSinceLastSync += ct.timeStepRate;
      ct.stepsSinceStart += ct.timeStepRate;
      for (auto &neighbor : neighbors) {
          const bool sendMessageTime = ct.correctionTime >= neighbor.ct.predictionTime - timeTolerance;
          const bool justBeforeSync = ct.stepsUntilSync <= ct.predictionsSinceLastSync;
          const bool sendMessageSteps = justBeforeSync
                  || ct.stepsSinceLastSync >= neighbor.ct.predictionsSinceLastSync;
        if (sendMessageSteps) {
          AdvancedCorrectionTimeMessage message{};
          message.time = ct.correctionTime;
          message.stepsSinceSync = ct.stepsSinceLastSync;
          neighbor.outbox->push(message);
        }
      }
      state = ActorState::Corrected;
    } else {
      result.yield = !processMessages();
    }
    break;
  }
  case ActorState::Synced:
    if (ct.stepsSinceLastSync == 0) {
      start();
      state = ActorState::Corrected;
    } else {
      result.yield = true;
    }
    break;
  default: throw;
  }
  const auto currentTime = std::chrono::steady_clock::now();
  result.isStateChanged = stateBefore != state;
  if (!result.isStateChanged) {
    const auto timeSinceLastUpdate = currentTime - lastStateChange;
    if (timeSinceLastUpdate > timeout && !alreadyPrintedTimeOut) {
        alreadyPrintedTimeOut = true;
        printTimeoutMessage(std::chrono::duration_cast<std::chrono::seconds>(timeSinceLastUpdate));
    }
  } else {
    lastStateChange = currentTime;
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
        // TODO(Lukas) Add asserts to check if we're in the correct state.
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


bool AbstractTimeCluster::maySync() {
    processMessages(); // TODO(Lukas) Do we actually need to do this here?
    return ct.stepsSinceLastSync >= ct.stepsUntilSync;
}

void AbstractTimeCluster::connect(AbstractTimeCluster &other) {
  neighbors.emplace_back(other.ct.maxTimeStepSize, other.ct.timeStepRate);
  other.neighbors.emplace_back(ct.maxTimeStepSize, ct.timeStepRate);
  neighbors.back().inbox = std::make_shared<MessageQueue>();
  other.neighbors.back().inbox = std::make_shared<MessageQueue>();
  neighbors.back().outbox = other.neighbors.back().inbox;
  other.neighbors.back().outbox = neighbors.back().inbox;
}

void AbstractTimeCluster::updateSyncTime(double newSyncTime) {
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
  for (auto& neighbor : neighbors) {
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

int AbstractTimeCluster::getPriority() const {
  return priority;
}

void AbstractTimeCluster::setPriority(int newPriority) {
  this->priority = newPriority;
}

ActorState AbstractTimeCluster::getState() const {
  return state;
}

}
