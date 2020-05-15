#include <iostream>
#include "Solver/time_stepping/AbstractTimeCluster.h"

namespace seissol::time_stepping {
double AbstractTimeCluster::timeStepSize() const {
  return ct.timeStepSize(syncTime);
}

AbstractTimeCluster::AbstractTimeCluster(double maxTimeStepSize, double timeTolerance, int timeStepRate)
    : timeTolerance(timeTolerance), timeStepRate(timeStepRate), numberOfTimeSteps(0) {
  ct.maxTimeStepSize = maxTimeStepSize;
  ct.timeStepRate = timeStepRate;
}

bool AbstractTimeCluster::act() {
  bool yield = false;
  switch (state) {
  case ActorState::Corrected: {
    if (maySync()) {
        std::cout << "synced at " << syncTime
        << ", corrTIme =" << ct.correctionTime
        << ", time tolerence " << timeTolerance
        << " stepsSinceLastSync " << ct.stepsSinceLastSync
        << " stepsUntilLastSync " << ct.stepsSinceLastSync
        << std::endl;
      state = ActorState::Synced;
    } else if (mayPredict()) {
      predict();
      ct.predictionsSinceLastSync += ct.timeStepRate;
      ct.predictionTime += timeStepSize();

      for (auto &neighbor : neighbors) {
          const bool sendMessageTime = ct.predictionTime >= (neighbor.ct.nextCorrectionTime(syncTime) - timeTolerance);
          // TODO(Lukas) Does this handle sync points correctly?
          // Maybe check also how many steps neighbor has to sync!
          const bool justBeforeSync = (ct.stepsUntilSync - ct.predictionsSinceLastSync) == 0;
          const bool sendMessageSteps = justBeforeSync
                  || ct.predictionsSinceLastSync >= (neighbor.ct.stepsSinceLastSync + neighbor.ct.timeStepRate);
          std::cerr << "AdvancedPredictionTimeMessage: justBeforeSync = " << justBeforeSync
          << " our preds = " << ct.predictionsSinceLastSync
          << " their preds = " << neighbor.ct.stepsSinceLastSync
          << " their new preds = " << neighbor.ct.stepsSinceLastSync + neighbor.ct.timeStepRate
          << " our pred time = " << ct.predictionTime
          << " their next corr time = " << neighbor.ct.nextCorrectionTime(syncTime)
          << std::endl;
          assert(sendMessageSteps == sendMessageTime);
        if (sendMessageTime) {
          AdvancedPredictionTimeMessage message{};
          message.time = ct.predictionTime;
          message.stepsSinceSync = ct.predictionsSinceLastSync;
          neighbor.outbox->push(message);
        }
      }
      state = ActorState::Predicted;
    } else {
      yield = !processMessages();
    }
    break;
  }
  case ActorState::Predicted: {
    if (mayCorrect()) {
      correct();
      ct.correctionTime += timeStepSize();
      ++numberOfTimeSteps;
      ct.stepsSinceLastSync += ct.timeStepRate;
      for (auto &neighbor : neighbors) {
          const bool sendMessageTime = ct.correctionTime >= neighbor.ct.predictionTime - timeTolerance;
          const bool justBeforeSync = (ct.stepsUntilSync - ct.predictionsSinceLastSync) == 0;
          const bool sendMessageSteps = justBeforeSync
                  || ct.stepsSinceLastSync >= neighbor.ct.predictionsSinceLastSync;
          std::cerr << "AdvancedCorrectionTimeMessage: justBeforeSync = " << justBeforeSync
                    << " our corrs = " << ct.stepsUntilSync
                    << " their preds = " << neighbor.ct.predictionsSinceLastSync
                    << " our time = " << ct.correctionTime
                    << " their time = " << neighbor.ct.correctionTime
                    << " their new preds = " << neighbor.ct.stepsSinceLastSync + neighbor.ct.timeStepRate
                    << std::endl;
          assert(sendMessageTime == sendMessageSteps);
        if (sendMessageTime) {
          AdvancedCorrectionTimeMessage message{};
          message.time = ct.correctionTime;
          message.stepsSinceSync = ct.stepsSinceLastSync;
          neighbor.outbox->push(message);
        }
      }
      state = ActorState::Corrected;
    } else {
      yield = !processMessages();
    }
    break;
  }
  case ActorState::Synced:
    if (ct.correctionTime + timeTolerance < syncTime) {
      start();
      state = ActorState::Corrected;
    } else {
      yield = true;
    }
    break;
  default: throw;
  }
  return yield;
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
          //neighbor.ct.predictionsSinceLastSync += neighbor.ct.timeStepRate;
          neighbor.ct.predictionsSinceLastSync = msg.stepsSinceSync;
          handleAdvancedPredictionTimeMessage(neighbor);
        } else if constexpr (std::is_same_v<T, AdvancedCorrectionTimeMessage>) {
          assert(msg.time > neighbor.ct.correctionTime);
          neighbor.ct.correctionTime = msg.time;
          //neighbor.ct.stepsSinceLastSync += neighbor.ct.timeStepRate;
          neighbor.ct.stepsSinceLastSync = msg.stepsSinceSync;
          handleAdvancedCorrectionTimeMessage(neighbor);
          std::cout << "Neighbor corrected, rate = " << neighbor.ct.timeStepRate
          <<  " our rate " << ct.timeStepRate
          << " neighbor steps since sync " << neighbor.ct.stepsSinceLastSync
          << " out steps since sync " << ct.stepsSinceLastSync
          << std::endl;
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
  const auto minNeighbor = std::min_element(
      neighbors.begin(), neighbors.end(),
      [this](NeighborCluster const &a, NeighborCluster const &b) {
        return a.ct.nextCorrectionTime(syncTime) < b.ct.nextCorrectionTime(syncTime);
      });
  const bool timeBasedPredict =
          minNeighbor == neighbors.end()
          || ct.predictionTime < minNeighbor->ct.nextCorrectionTime(syncTime) - timeTolerance;
    bool stepBasedPredict = true;
    for (const auto& neighbor : neighbors) {
        const auto neighStepsAfterUpdate = std::min(
                neighbor.ct.stepsSinceLastSync + neighbor.ct.timeStepRate,
                ct.stepsUntilSync
                );
        const bool curStepBasedPredict = (ct.predictionsSinceLastSync < neighStepsAfterUpdate);
        stepBasedPredict = stepBasedPredict && curStepBasedPredict;
    }
    assert(timeBasedPredict == stepBasedPredict);

    return stepBasedPredict;
}

bool AbstractTimeCluster::mayCorrect() {
  // We can correct, if our prediction time is smaller than the one of all neighbors.
  const auto minNeighbor = std::min_element(
      neighbors.begin(), neighbors.end(),
      [](NeighborCluster const &a, NeighborCluster const &b) {
        return a.ct.predictionTime < b.ct.predictionTime;
      });
  const bool timeBasedCorrect =
      minNeighbor == neighbors.end()
      || ct.predictionTime <= minNeighbor->ct.predictionTime + timeTolerance;

  bool stepBasedCorrect = true;
  for (auto& neighbor : neighbors) {
      //const double rateDiff = (1.0*neighbor.ct.timeStepRate) / ct.timeStepRate;
      //const int ourPredictions = static_cast<int>(std::round(
              //rateDiff * ct.predictionsSinceLastSync));
      const bool isSynced = (neighbor.ct.stepsUntilSync - neighbor.ct.predictionsSinceLastSync) == 0;
      stepBasedCorrect = stepBasedCorrect
              && (isSynced || (ct.predictionsSinceLastSync <= neighbor.ct.predictionsSinceLastSync));
      //std::cout << stepBasedCorrect << std::endl;
      /*
      std::cout
              << "rateDiff = " << rateDiff
              << " ourRate = " << ct.timeStepRate
              << " theirRate = " << neighbor.ct.timeStepRate
              << " ourDt = " << ct.maxTimeStepSize
              << " theirDt = " << neighbor.ct.maxTimeStepSize
              << " ourPreds = " << ourPredictions
              << " ourPredsWOFact = " << ct.predictionsSinceLastSync
              << " neighborPreds = " << neighbor.ct.predictionsSinceLastSync
              << " our pred.time = " <<  ct.predictionTime
              << " their pred.time = " << neighbor.ct.predictionTime
              << " mayCorrect (step based) = " << stepBasedCorrect
              << std::endl;
      const auto ourPredComp = ct.predictionsSinceLastSync * ct.maxTimeStepSize / ct.timeStepRate;
      const auto neighPredComp = neighbor.ct.predictionsSinceLastSync * neighbor.ct.maxTimeStepSize / ct.timeStepRate;

      std::cout
              << "our predict * rate * timesteps = " << ourPredComp
              << " neigh predict * rate * timesteps = " << neighPredComp
              << std::endl;
              */
  }
  /*std::cout
  << "stepBased = " << stepBasedCorrect
  << " timeBased = " << timeBasedCorrect
  << std::endl;
   */
  assert(timeBasedCorrect == stepBasedCorrect);
  return stepBasedCorrect;
}


bool AbstractTimeCluster::maySync() {
    const bool timeBasedSync = ct.correctionTime + timeTolerance >= syncTime;
    const bool stepBasedSync = (ct.stepsUntilSync - ct.stepsSinceLastSync) == 0;
    assert(timeBasedSync == stepBasedSync);
    return stepBasedSync && processMessages();
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
  processMessages();
    for (auto& neighbor : neighbors) {
        assert(!neighbor.inbox->hasMessages());
    }
  /*
  for (auto& neighbor : neighbors) {
    // TODO(Lukas) Think this through!
    neighbor.inbox->clear();
    neighbor.outbox->clear();
    neighbor.ct.predictionTime = ct.predictionTime;
    neighbor.ct.correctionTime = ct.correctionTime;
  }
   */

  // TODO(Lukas): Don't -> just to try out scheduling
  // If this works, create numberOfTimeStepsSinceSyn

  // Always use our time, as messages from neighbors may have to be processed.
  numberOfTimeSteps = 0;
  ct.stepsSinceLastSync = 0;
  ct.predictionsSinceLastSync = 0;
  ct.stepsUntilSync = ct.computeStepsUntilSyncTime(ct.correctionTime, syncTime);
  for (auto& neighbor : neighbors) {
    neighbor.ct.stepsUntilSync = neighbor.ct.computeStepsUntilSyncTime(ct.correctionTime, syncTime);
    neighbor.ct.stepsSinceLastSync = 0;
    neighbor.ct.predictionsSinceLastSync = 0;
  }

}

}
