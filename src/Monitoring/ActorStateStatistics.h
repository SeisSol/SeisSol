#ifndef SEISSOL_ACTORSTATESTATISTICS_H
#define SEISSOL_ACTORSTATESTATISTICS_H

#include <unordered_map>
#include <list>
#include <optional>
#include "LoopStatistics.h"
#include "Solver/time_stepping/ActorState.h"

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
    timespec begin;
    std::optional<timespec> end;
    int numEnteredRegion;
    Sample() = delete;
  };

  Sample currentSample;
  bool started = false;

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
