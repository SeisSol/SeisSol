// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ActorStateStatistics.h"
#include "LoopStatistics.h"
#include <Solver/Clustering/ActorState.h>
#include <optional>
#include <time.h>
namespace seissol {

ActorStateStatistics::ActorStateStatistics(unsigned globalClusterId, LoopStatistics& loopStatistics)
    : currentSample(seissol::solver::clustering::ComputeStep::Correct),
      globalClusterId(globalClusterId), loopStatistics(loopStatistics) {}

void ActorStateStatistics::enter(seissol::solver::clustering::ComputeStep actorState) {
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
  const auto region = loopStatistics.getRegion(
      seissol::solver::clustering::actorStateToString(seissol::solver::clustering::ActorState{
          seissol::solver::clustering::StateType::ComputeDone, state}));
  loopStatistics.addSample(
      region, 1, globalClusterId, currentSample.begin, currentSample.end.value());
}

ActorStateStatistics::Sample::Sample(seissol::solver::clustering::ComputeStep state)
    : state(state), end(std::nullopt), numEnteredRegion(0) {
  clock_gettime(CLOCK_MONOTONIC, &begin);
}
void ActorStateStatistics::Sample::finish() {
  timespec endTime{};
  clock_gettime(CLOCK_MONOTONIC, &endTime);
  end = endTime;
}

ActorStateStatisticsManager::ActorStateStatisticsManager(LoopStatistics& loopStatistics)
    : loopStatistics(loopStatistics) {
  loopStatistics.addRegion(
      seissol::solver::clustering::actorStateToString(seissol::solver::clustering::ActorState{
          seissol::solver::clustering::StateType::ComputeDone,
          seissol::solver::clustering::ComputeStep::Predict}),
      false);
  loopStatistics.addRegion(
      seissol::solver::clustering::actorStateToString(seissol::solver::clustering::ActorState{
          seissol::solver::clustering::StateType::ComputeDone,
          seissol::solver::clustering::ComputeStep::Interact}),
      false);
  loopStatistics.addRegion(
      seissol::solver::clustering::actorStateToString(seissol::solver::clustering::ActorState{
          seissol::solver::clustering::StateType::ComputeDone,
          seissol::solver::clustering::ComputeStep::Correct}),
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
