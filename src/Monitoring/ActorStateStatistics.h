#ifndef SEISSOL_ACTORSTATESTATISTICS_H
#define SEISSOL_ACTORSTATESTATISTICS_H

#include "LoopStatistics.h"
#include "Solver/Clustering/ActorState.h"
#include <list>
#include <optional>
#include <unordered_map>

namespace seissol {
class ActorStateStatistics {
  public:
  ActorStateStatistics(unsigned globalClusterId, LoopStatistics& loopStatistics);

  void enter(seissol::solver::clustering::ComputeStep actorState);
  void exit();

  private:
  struct Sample {
    explicit Sample(seissol::solver::clustering::ComputeStep state);
    void finish();
    seissol::solver::clustering::ComputeStep state;
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

#endif // SEISSOL_ACTORSTATESTATISTICS_H
