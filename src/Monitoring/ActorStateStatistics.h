// SPDX-FileCopyrightText: 2021-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MONITORING_ACTORSTATESTATISTICS_H_
#define SEISSOL_SRC_MONITORING_ACTORSTATESTATISTICS_H_

#include "LoopStatistics.h"
#include "Solver/time_stepping/ActorState.h"
#include <list>
#include <optional>
#include <unordered_map>

namespace seissol {
class ActorStateStatistics {
  public:
  ActorStateStatistics(unsigned globalClusterId, LoopStatistics& loopStatistics);

  void enter(time_stepping::ActorState actorState);
  void exit();

  private:
  struct Sample {
    explicit Sample(seissol::time_stepping::ActorState state);
    void finish();
    seissol::time_stepping::ActorState state;
    timespec begin{};
    std::optional<timespec> end;
    int numEnteredRegion;
    Sample() = delete;
  };

  Sample currentSample;

  unsigned globalClusterId;
  LoopStatistics& loopStatistics;
};

class ActorStateStatisticsManager {
  public:
  ActorStateStatisticsManager(LoopStatistics& loopStatistics);
  ActorStateStatistics& addCluster(unsigned globalClusterId);

  void finish();

  private:
  // for now, use an std::list, as this is merely a storage for all ActorStateStatistics objects
  std::list<ActorStateStatistics> stateStatistics;
  LoopStatistics& loopStatistics;
};
} // namespace seissol

#endif // SEISSOL_SRC_MONITORING_ACTORSTATESTATISTICS_H_
