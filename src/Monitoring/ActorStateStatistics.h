#ifndef SEISSOL_ACTORSTATESTATISTICS_H
#define SEISSOL_ACTORSTATESTATISTICS_H

#include "Initializer/tree/Layer.hpp"
#include "Solver/time_stepping/ActorState.h"
#include <chrono>
#include <functional>
#include <unordered_map>
#include <vector>

namespace seissol {
using TimeClusterKey_t = std::pair<unsigned, LayerType>;
}; // namespace seissol

namespace std {

template <>
struct hash<seissol::TimeClusterKey_t> {
  std::size_t operator()(const seissol::TimeClusterKey_t& k) const;
};

} // namespace std

namespace seissol {
class ActorStateStatistics {
  public:
  enum class EventType { Start, Stop };
  struct Event {
    Event();
    void updateTime();
    EventType type;
    time_stepping::ActorState state;
    // TODO(Lukas) Timestamp
    unsigned threadId;
    std::chrono::time_point<std::chrono::steady_clock> time;
  };

  class Guard {
public:
    Guard(ActorStateStatistics& statistics, Event eventTemplate);
    ~Guard();

private:
    ActorStateStatistics& statistics;
    const Event eventTemplate;
  };

  ActorStateStatistics();
  ActorStateStatistics(unsigned globalClusterId, LayerType type);

  ~ActorStateStatistics();

  void emitEvent(Event event);
  [[nodiscard]] Guard enterWithGuard(Event eventTemplate);

  std::string formatEvents();

  private:
  unsigned globalClusterId;
  LayerType type;

  omp_lock_t eventsLock;
  std::vector<Event> events;
};

class ActorStateStatisticsManager {
  public:
  ActorStateStatisticsManager() = default;
  ActorStateStatistics& addCluster(unsigned globalClusterId, LayerType type);
  void printSummary();

  private:
  std::unordered_map<TimeClusterKey_t, ActorStateStatistics> stateStatisticsMap{};
};
} // namespace seissol

#endif // SEISSOL_ACTORSTATESTATISTICS_H
