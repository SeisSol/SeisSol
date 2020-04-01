#include "Solver/time_stepping/AbstractTimeCluster.h"

namespace seissol::time_stepping {
double AbstractTimeCluster::timeStepSize() const {
  return ct.timeStepSize(syncTime);
}

AbstractTimeCluster::AbstractTimeCluster(double maxTimeStepSize, double timeTolerance)
    : timeTolerance(timeTolerance) {
  ct.maxTimeStepSize = maxTimeStepSize;
}

bool AbstractTimeCluster::act() {
  bool yield = false;
  switch (state) {
  case ActorState::Corrected: {
    if (ct.correctionTime + timeTolerance >= syncTime) {
      state = ActorState::Synced;
    } else if (mayPredict()) {
      predict();
      //std::cout << "pred dt_max=" << ct.maxTimeStepSize  << " dt=" << timeStepSize() <<
      //" t_p=" << ct.predictionTime << " t_c=" << ct.correctionTime << " reset=" << resetBuffers <<
      //" t_minnext=" << minNeighbor->ct.nextCorrectionTime(syncTime) <<
      //" t_delta=" << ct.predictionTime - minNeighbor->ct.nextCorrectionTime(syncTime) << std::endl;
      ct.predictionTime += timeStepSize();
      for (auto &neighbor : neighbors) {
        /*std::cout << ct.maxTimeStepSize << " sends?? " << ct.predictionTime << " " << neighbor.ct.nextCorrectionTime(syncTime) <<
        " " << (ct.predictionTime >= neighbor.ct.nextCorrectionTime(syncTime)) << std::endl;*/
        if (ct.predictionTime >= neighbor.ct.nextCorrectionTime(syncTime) - timeTolerance) {
          AdvancedPredictionTimeMessage message{};
          message.time = ct.predictionTime;
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
      //std::cout << "corr dt_max=" << ct.maxTimeStepSize << " dt=" << timeStepSize()
      //<< " t_p=" << ct.predictionTime << " t_c=" << ct.correctionTime
      //<< " t_sub=" << subTimeStart << std::endl;
      ct.correctionTime += timeStepSize();
      for (auto &neighbor : neighbors) {
        if (ct.correctionTime >= neighbor.ct.predictionTime - timeTolerance) {
          AdvancedCorrectionTimeMessage message{};
          message.time = ct.correctionTime;
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
          handleAdvancedPredictionTimeMessage(neighbor);
        } else if constexpr (std::is_same_v<T, AdvancedCorrectionTimeMessage>) {
          assert(msg.time > neighbor.ct.correctionTime);
          neighbor.ct.correctionTime = msg.time;
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
  const auto minNeighbor = std::min_element(
      neighbors.begin(), neighbors.end(),
      [this](NeighborCluster const &a, NeighborCluster const &b) {
        return a.ct.nextCorrectionTime(syncTime) < b.ct.nextCorrectionTime(syncTime);
      });
  return minNeighbor == neighbors.end()
      || ct.predictionTime < minNeighbor->ct.nextCorrectionTime(syncTime) - timeTolerance;
}

bool AbstractTimeCluster::mayCorrect() {
  const auto minNeighbor = std::min_element(
      neighbors.begin(), neighbors.end(),
      [](NeighborCluster const &a, NeighborCluster const &b) {
        return a.ct.predictionTime < b.ct.predictionTime;
      });
  return minNeighbor == neighbors.end()
      || ct.predictionTime <= minNeighbor->ct.predictionTime + timeTolerance;
}

void AbstractTimeCluster::connect(AbstractTimeCluster &other) {
  neighbors.emplace_back(other.ct.maxTimeStepSize);
  other.neighbors.emplace_back(ct.maxTimeStepSize);
  neighbors.back().inbox = std::make_shared<MessageQueue>();
  other.neighbors.back().inbox = std::make_shared<MessageQueue>();
  neighbors.back().outbox = other.neighbors.back().inbox;
  other.neighbors.back().outbox = neighbors.back().inbox;
}

void AbstractTimeCluster::updateSyncTime(double newSyncTime) {
  assert(newSyncTime > syncTime);
  syncTime = newSyncTime;
}

bool AbstractTimeCluster::synced() const {
  return state == ActorState::Synced;
}

}
