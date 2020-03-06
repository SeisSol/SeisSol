#include "GhostTimeCluster.h"

namespace seissol::time_stepping {
  bool GhostTimeCluster::act() {
    bool yield = false;
    /*switch (state) {
      case State::Progress: {
        processMessages()
        break;
      }
      case ActorState::Synced: {
        for (auto& neighbor : neighbors) {
          if (neighbor.inbox->hasMessages()) {
            Message message = neighbor.inbox->pop();
            std::visit([&neighbor,this](auto&& msg) {
              using T = std::decay_t<decltype(msg)>;
              if constexpr (std::is_same_v<T, AdvancedPredictionTimeMessage>) {
                assert(msg.time > neighbor.ct.predictionTime);
                neighbor.ct.predictionTime = msg.time;

              } else if constexpr (std::is_same_v<T, AdvancedCorrectionTimeMessage>) {
                assert(msg.time > neighbor.ct.correctionTime);
                neighbor.ct.correctionTime = msg.time;
              } else {
                static_assert(always_false<T>::value, "non-exhaustive visitor!");
              }
            }, message);
          }
        }
        break;
      }
      default:
        throw;
    }*/
    return yield;
  }

  /*bool GhostTimeCluster::processMessages() {
    bool processed = false;
    for (auto& neighbor : neighbors) {
      if (neighbor.inbox->hasMessages()) {
        processed = true;
        Message message = neighbor.inbox->pop();
        std::visit([&neighbor,this](auto&& msg) {
          using T = std::decay_t<decltype(msg)>;
          if constexpr (std::is_same_v<T, AdvancedPredictionTimeMessage>) {
            assert(msg.time > neighbor.ct.predictionTime);
            neighbor.ct.predictionTime = msg.time;
            if (neighbor.ct.maxTimeStepSize > ct.maxTimeStepSize) {
              lastSubTime = neighbor.ct.correctionTime;
            }
            //std::cout << "AdvancedPred " << ct.maxTimeStepSize << " " << msg.time << std::endl;
          } else if constexpr (std::is_same_v<T, AdvancedCorrectionTimeMessage>) {
            assert(msg.time > neighbor.ct.correctionTime);
            neighbor.ct.correctionTime = msg.time;
            //std::cout << "AdvancedCorr " << ct.maxTimeStepSize << " " << msg.time << std::endl;
          } else {
            static_assert(always_false<T>::value, "non-exhaustive visitor!");
          }
        }, message);
      }
    }
    return processed;
}

  bool TimeCluster::synced() const {
    return state == ActorState::Synced;
  }*/
}
