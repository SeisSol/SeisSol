#include "ActorStateStatistics.h"
#include "LoopStatistics.h"
#include <Solver/Clustering/ActorState.h>
#include <optional>
#include <time.h>
namespace seissol {

ActorStateStatistics::ActorStateStatistics(unsigned globalClusterId, LoopStatistics& loopStatistics)
    : currentSample(time_stepping::ComputeStep::Correct), globalClusterId(globalClusterId),
      loopStatistics(loopStatistics) {}

void ActorStateStatistics::enter(time_stepping::ComputeStep actorState) {
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
  const auto region = loopStatistics.getRegion(seissol::time_stepping::actorStateToString(
      time_stepping::ActorState{time_stepping::StateType::ComputeDone, state}));
  loopStatistics.addSample(
      region, 1, globalClusterId, currentSample.begin, currentSample.end.value());
}

ActorStateStatistics::Sample::Sample(seissol::time_stepping::ComputeStep state)
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
  loopStatistics.addRegion(
      time_stepping::actorStateToString(time_stepping::ActorState{
          time_stepping::StateType::ComputeDone, time_stepping::ComputeStep::Predict}),
      false);
  loopStatistics.addRegion(
      time_stepping::actorStateToString(time_stepping::ActorState{
          time_stepping::StateType::ComputeDone, time_stepping::ComputeStep::Interact}),
      false);
  loopStatistics.addRegion(
      time_stepping::actorStateToString(time_stepping::ActorState{
          time_stepping::StateType::ComputeDone, time_stepping::ComputeStep::Correct}),
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
