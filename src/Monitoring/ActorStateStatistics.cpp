#include "ActorStateStatistics.h"
#include <fstream>
#include <string>
#include <iostream>
#include "Parallel/MPI.h"

namespace std {
std::size_t hash<seissol::TimeClusterKey_t>::operator()(const TimeClusterKey_t& k) const {
  using std::hash;

  return ((hash<unsigned>()(k.first) ^ (hash<LayerType>()(k.second) << 1)) >> 1);
}
} // namespace std

namespace seissol {
void ActorStateStatistics::Event::updateTime() { time = std::chrono::steady_clock::now(); }

ActorStateStatistics::Event::Event() { updateTime(); }

ActorStateStatistics::Guard::Guard(ActorStateStatistics& statistics,
                                   ActorStateStatistics::Event eventTemplate)
    : statistics(statistics), eventTemplate(eventTemplate) {
  auto event = eventTemplate;
  event.type = EventType::Start;
  event.threadId = omp_get_thread_num();
  statistics.emitEvent(event);
}

ActorStateStatistics::Guard::~Guard() {
  auto event = eventTemplate;
  event.updateTime();
  event.type = EventType::Stop;
  event.threadId = omp_get_thread_num(); // Note: May be different to start thread!
  statistics.emitEvent(event);
}

ActorStateStatistics::ActorStateStatistics() {}

ActorStateStatistics::ActorStateStatistics(unsigned int globalClusterId, LayerType type)
    : globalClusterId(globalClusterId), type(type) {}
ActorStateStatistics::~ActorStateStatistics() {}
void ActorStateStatistics::emitEvent(ActorStateStatistics::Event event) {
#pragma omp critical
  events.emplace_back(event);
}
ActorStateStatistics::Guard
    ActorStateStatistics::enterWithGuard(ActorStateStatistics::Event eventTemplate) {
  return Guard(*this, eventTemplate);
}

std::string ActorStateStatistics::formatEvents() {
  using namespace std::chrono;

  std::stringstream formattedEvents;
  assert(type != LayerType::Ghost);
  const std::string layerTypeStr = type == LayerType::Interior ? "interior" : "copy";
  const auto rank = seissol::MPI::mpi.rank();

  for (const auto& event : events) {
    const std::string typeStr = event.type == EventType::Start ? "start" : "stop";
    const auto actorStateStr = time_stepping::actorStateToString(event.state);
    auto timestamp = duration_cast<microseconds>(event.time.time_since_epoch());

    formattedEvents << rank << ",";
    formattedEvents << globalClusterId << ",";
    formattedEvents << layerTypeStr << ",";
    formattedEvents << event.threadId << ",";
    formattedEvents << typeStr << ",";
    formattedEvents << actorStateStr << ",";
    formattedEvents << timestamp.count() << "\n";
  }
  std::cout << formattedEvents.str();
  return formattedEvents.str();
}

ActorStateStatistics& ActorStateStatisticsManager::addCluster(unsigned int globalClusterId,
                                                              LayerType type) {
  auto key = TimeClusterKey_t{globalClusterId, type};
  if (auto it = stateStatisticsMap.find(key); it != stateStatisticsMap.end())
    return it->second;
  else {
    auto [insertIt, _] =
        stateStatisticsMap.insert({key, ActorStateStatistics(globalClusterId, type)});
    return insertIt->second;
  }
}
void ActorStateStatisticsManager::printSummary() {
  const auto rank = seissol::MPI::mpi.size();
  const auto filename = std::string("actor-") + std::to_string(rank) + ".csv";
  std::ofstream out(filename);
  // Print header
  out << "rank,globalClusterId,layerType,threadId,eventType,actorState,timestamp\n";

  for (auto& [key, value] : stateStatisticsMap) {
    out << value.formatEvents();
  }
}
} // namespace seissol
