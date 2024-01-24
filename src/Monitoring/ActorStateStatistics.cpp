#include "ActorStateStatistics.h"
#include "LoopStatistics.h"

namespace seissol {

ActorStateStatistics::ActorStateStatistics(unsigned globalClusterId, LoopStatistics& loopStatistics)
    : currentSample(time_stepping::ActorState::Synced), loopStatistics(loopStatistics),
      globalClusterId(globalClusterId) {}

void ActorStateStatistics::enter(time_stepping::ActorState actorState) {
  if (actorState == currentSample.state) {
    ++currentSample.numEnteredRegion;
  } else {
    exit();
    currentSample = Sample(actorState);
  }
}

void ActorStateStatistics::exit() {
  currentSample.finish();
  const auto state = currentSample.state;
  const auto region = loopStatistics.getRegion(seissol::time_stepping::actorStateToString(state));
  loopStatistics.addSample(
      region, 1, globalClusterId, currentSample.begin, currentSample.end.value());
}

ActorStateStatistics::Sample::Sample(seissol::time_stepping::ActorState state)
    : state(state), end(std::nullopt), numEnteredRegion(0) {
  clock_gettime(CLOCK_MONOTONIC, &begin);
}
void ActorStateStatistics::Sample::finish() {
  timespec endTime;
  clock_gettime(CLOCK_MONOTONIC, &endTime);
  end = endTime;
}

ActorStateStatisticsManager::ActorStateStatisticsManager(LoopStatistics& loopStatistics)
    : loopStatistics(loopStatistics) {
  loopStatistics.addRegion(time_stepping::actorStateToString(time_stepping::ActorState::Synced),
                           false);
  loopStatistics.addRegion(time_stepping::actorStateToString(time_stepping::ActorState::Corrected),
                           false);
  loopStatistics.addRegion(time_stepping::actorStateToString(time_stepping::ActorState::Predicted),
                           false);
}

ActorStateStatistics& ActorStateStatisticsManager::addCluster(unsigned globalClusterId) {
  return stateStatistics.emplace_back(globalClusterId, loopStatistics);
}

void ActorStateStatisticsManager::finish() {
  for (auto& stateStatistics : stateStatistics) {
    stateStatistics.exit();
  }
}

} // namespace seissol
