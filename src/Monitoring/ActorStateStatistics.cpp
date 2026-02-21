// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ActorStateStatistics.h"

#include "LoopStatistics.h"
#include "Solver/TimeStepping/ActorState.h"

#include <optional>
#include <time.h>

namespace seissol {

ActorStateStatistics::ActorStateStatistics(unsigned globalClusterId, LoopStatistics& loopStatistics)
    : currentSample_(time_stepping::ActorState::Synced), globalClusterId_(globalClusterId),
      loopStatistics_(loopStatistics) {}

void ActorStateStatistics::enter(time_stepping::ActorState actorState) {
  if (actorState == currentSample_.state) {
    ++currentSample_.numEnteredRegion;
  } else {
    exit();
    currentSample_ = Sample(actorState);
  }
}

void ActorStateStatistics::exit() {
  currentSample_.finish();
  const auto state = currentSample_.state;
  const auto region = loopStatistics_.getRegion(seissol::time_stepping::actorStateToString(state));
  loopStatistics_.addSample(
      region, 1, globalClusterId_, currentSample_.begin, currentSample_.end.value());
}

ActorStateStatistics::Sample::Sample(seissol::time_stepping::ActorState state)
    : state(state), end(std::nullopt), numEnteredRegion(0) {
  (void)clock_gettime(CLOCK_MONOTONIC, &begin);
}
void ActorStateStatistics::Sample::finish() {
  timespec endTime{};
  (void)clock_gettime(CLOCK_MONOTONIC, &endTime);
  end = endTime;
}

ActorStateStatisticsManager::ActorStateStatisticsManager(LoopStatistics& loopStatistics)
    : loopStatistics_(loopStatistics) {
  loopStatistics.addRegion(time_stepping::actorStateToString(time_stepping::ActorState::Synced),
                           false);
  loopStatistics.addRegion(time_stepping::actorStateToString(time_stepping::ActorState::Corrected),
                           false);
  loopStatistics.addRegion(time_stepping::actorStateToString(time_stepping::ActorState::Predicted),
                           false);
}

ActorStateStatistics& ActorStateStatisticsManager::addCluster(unsigned globalClusterId) {
  return stateStatistics_.emplace_back(globalClusterId, loopStatistics_);
}

void ActorStateStatisticsManager::finish() {
  for (auto& stateStatistics : stateStatistics_) {
    stateStatistics.exit();
  }
}

} // namespace seissol
