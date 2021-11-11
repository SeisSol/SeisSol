#ifndef SEISSOL_ACTORSTATESTATISTICS_H
#define SEISSOL_ACTORSTATESTATISTICS_H

#include <unordered_map>
#include <vector>
#include <optional>
#include "Solver/time_stepping/ActorState.h"

namespace seissol {
class ActorStateStatistics {
public:
  ActorStateStatistics() : currentSample(time_stepping::ActorState::Synced) {
  }

  void enter(time_stepping::ActorState actorState) {
    if (actorState == currentSample.state) {
      ++currentSample.numEnteredRegion;
    } else {
      currentSample.finish();
      samples.push_back(currentSample);
      currentSample = Sample(actorState);
    }
  }

  void addToLoopStatistics(unsigned globalClusterId, LoopStatistics& loopStatistics) {
    currentSample.finish();
    samples.push_back(currentSample);
    for (const auto& sample : samples) {
      const auto state = sample.state;
      const auto region = loopStatistics.getRegion(
          time_stepping::actorStateToString(state)
          );
      loopStatistics.addSample(region,
                                sample.numEnteredRegion,
                                globalClusterId,
                                sample.begin,
                                sample.end.value());
    }
  }
private:

  struct Sample {
    explicit Sample(time_stepping::ActorState state) : state(state), end(std::nullopt), numEnteredRegion(0) {
      clock_gettime(CLOCK_MONOTONIC, &begin);
    }
    void finish() {
      timespec endTime;
      clock_gettime(CLOCK_MONOTONIC, &endTime);
      end = endTime;
    }
    time_stepping::ActorState state;
    timespec begin;
    std::optional<timespec> end;
    int numEnteredRegion;
    Sample() = delete;
  };

  Sample currentSample;
  std::vector<Sample> samples;


};

class ActorStateStatisticsManager {
public:
  ActorStateStatisticsManager() = default;
  ActorStateStatistics& addCluster(unsigned globalClusterId) {
    return stateStatisticsMap[globalClusterId];
  }

  void addToLoopStatistics(LoopStatistics& loopStatistics) {
    loopStatistics.addRegion(time_stepping::actorStateToString(time_stepping::ActorState::Synced), false);
    loopStatistics.addRegion(time_stepping::actorStateToString(time_stepping::ActorState::Corrected), false);
    loopStatistics.addRegion(time_stepping::actorStateToString(time_stepping::ActorState::Predicted), false);

    for (auto& [globalClusterId, stateStatistics] : stateStatisticsMap) {
      stateStatistics.addToLoopStatistics(globalClusterId, loopStatistics);
    }
  }
private:
  std::unordered_map<unsigned, ActorStateStatistics> stateStatisticsMap{};
};
}

#endif //SEISSOL_ACTORSTATESTATISTICS_H
