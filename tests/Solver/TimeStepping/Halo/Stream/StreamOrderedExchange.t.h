// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/Executor.h"
#include "Kernels/Precision.h"
#include "Parallel/Runtime/EventPool.h"
#include "Solver/TimeStepping/Actor/AbstractTimeCluster.h"
#include "Solver/TimeStepping/Actor/ActorState.h"
#include "Solver/TimeStepping/Halo/GhostCluster.h"
#include "Solver/TimeStepping/Halo/HaloCommunication.h"
#include "Solver/TimeStepping/Halo/Stream/ExchangeScheduler.h"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

namespace seissol::unit_test {
using namespace seissol::solver;

namespace {

constexpr auto NoActor = std::numeric_limits<std::size_t>::max();

/**
 * The devices of several processes, and the network between them. Each device runs the items of
 * its streams one after another: work of a cluster, waits for events, records of events, and groups
 * of point-to-point operations, which complete once all their counterparts have been posted (no
 * buffering). Which stream runs next is random. Whenever work or a group runs, the rules of the
 * actor model are checked against what the devices have run so far.
 *
 * Events behave as with CUDA and HIP: a wait refers to the latest record of the event at the time
 * it is enqueued, and is satisfied right away if there is none; and a new event may get the handle
 * of one destroyed before, the one destroyed last first. The devices count the uses of destroyed
 * events, and the records of events that somebody still means to wait for: either would be a wait
 * for other work than meant, or a crash.
 *
 * While a process records its streams into a graph, they may only wait for events recorded inside
 * of the same recording; all other streams may only wait for events recorded outside of any (as
 * with CUDA graphs, where anything else fails). The model runs the recorded items right away; the
 * recorder has to order them as their replay would.
 */
class StreamDevices {
  public:
  using Event = void*;

  struct Actor {
    long rate{1};
    long stepsUntilSync{0};
    // for ghost clusters: the copy layer they stand next to
    std::size_t copy{NoActor};
    std::vector<std::size_t> neighbors;
    // what the device has run: predictions and corrections, for ghost clusters the progress of the
    // completed receives and sends
    long predictions{0};
    long corrections{0};
  };

  struct Group {
    int process;
    std::vector<std::pair<bool, int>> operations;
    std::vector<std::size_t> sequence;
    bool posted{false};
    std::size_t sender{NoActor};
    long sendTarget{0};
    std::size_t receiver{NoActor};
    long receiveTarget{0};
  };

  std::vector<Actor> actors;
  std::size_t violations{0};
  // waits across the border of a recording
  std::size_t isolationViolations{0};
  // waits for, and records of, destroyed events
  std::size_t destroyedUses{0};
  // records of events that somebody still means to wait for
  std::size_t prematureRecords{0};

  /// whether somebody still means to wait for the event (e.g. a neighbor that it gets published to)
  std::function<bool(Event)> awaited;

  /// a new event; may be one destroyed before, as with CUDA
  Event create() {
    if (!freed_.empty()) {
      auto* event = freed_.back();
      freed_.pop_back();
      destroyed_.erase(event);
      latestRecord_.erase(event);
      return event;
    }
    return &handles_.emplace_back();
  }

  void destroy(Event event) {
    destroyed_.insert(event);
    freed_.push_back(event);
  }

  /// a pool of events of the devices, as a stream runtime has one
  parallel::runtime::EventPool pool() {
    return {[this]() { return create(); }, [this](void* event) { destroy(event); }};
  }

  /// makes the stream one of the process, for its recordings
  void assign(std::size_t stream, int process) { processes_[stream] = process; }

  void beginRecording(int process) { recordings_[process] = ++lastRecording_; }
  void endRecording(int process) { recordings_.erase(process); }

  void work(std::size_t stream, std::size_t actor, ActorAction action) {
    streams_[stream].push_back({Kind::Work, actor, action, 0, 0});
  }
  void wait(std::size_t stream, Event event) {
    if (destroyed_.count(event) > 0) {
      ++destroyedUses;
      return;
    }
    // the latest record so far; none (0) means that there is nothing to wait for
    const auto found = latestRecord_.find(event);
    const auto record = found == latestRecord_.end() ? 0 : found->second;
    if ((record == 0 ? 0 : recordedIn_.at(record)) != recording(stream)) {
      ++isolationViolations;
    }
    streams_[stream].push_back({Kind::Wait, 0, ActorAction::Nothing, record, 0});
  }
  void record(std::size_t stream, Event event) {
    if (destroyed_.count(event) > 0) {
      ++destroyedUses;
    }
    if (awaited && awaited(event)) {
      ++prematureRecords;
    }
    const auto record = ++lastRecord_;
    latestRecord_[event] = record;
    recordedIn_[record] = recording(stream);
    streams_[stream].push_back({Kind::Record, 0, ActorAction::Nothing, record, 0});
  }
  void group(std::size_t stream, Group group) {
    groups_.push_back(std::move(group));
    streams_[stream].push_back({Kind::Group, 0, ActorAction::Nothing, 0, groups_.size() - 1});
  }

  /**
   * Runs one item that can run, chosen at random; returns false if there is none.
   */
  bool runOne(std::mt19937& random) {
    std::vector<std::size_t> runnable;
    for (auto& [stream, items] : streams_) {
      if (items.empty()) {
        continue;
      }
      auto& item = items.front();
      if (item.kind == Kind::Group && !groups_[item.group].posted) {
        post(groups_[item.group]);
      }
      if (canRun(item)) {
        runnable.push_back(stream);
      }
    }
    if (runnable.empty()) {
      return false;
    }
    const auto stream =
        runnable[std::uniform_int_distribution<std::size_t>(0, runnable.size() - 1)(random)];
    auto& items = streams_[stream];
    const auto item = items.front();
    items.pop_front();
    run(item);
    return true;
  }

  /**
   * Runs everything; returns false if the devices get stuck.
   */
  bool drain(std::mt19937& random) {
    while (runOne(random)) {
    }
    return std::all_of(
        streams_.begin(), streams_.end(), [](const auto& stream) { return stream.second.empty(); });
  }

  void resetInterval() {
    for (auto& actor : actors) {
      actor.predictions = 0;
      actor.corrections = 0;
    }
  }

  private:
  enum class Kind { Work, Wait, Record, Group };

  /// the recording the stream is part of; 0 if none
  [[nodiscard]] std::size_t recording(std::size_t stream) const {
    const auto process = processes_.find(stream);
    if (process == processes_.end()) {
      return 0;
    }
    const auto found = recordings_.find(process->second);
    return found == recordings_.end() ? 0 : found->second;
  }
  struct Item {
    Kind kind;
    std::size_t actor;
    ActorAction action;
    // for waits and records: the record
    std::size_t record;
    std::size_t group;
  };

  void post(Group& group) {
    for (const auto& [send, peer] : group.operations) {
      auto& counter =
          send ? postedSends_[{group.process, peer}] : postedReceives_[{peer, group.process}];
      group.sequence.push_back(counter++);
    }
    group.posted = true;
  }

  [[nodiscard]] bool canRun(const Item& item) const {
    switch (item.kind) {
    case Kind::Wait:
      return item.record == 0 || completed_.count(item.record) > 0;
    case Kind::Group: {
      const auto& group = groups_[item.group];
      for (std::size_t i = 0; i < group.operations.size(); ++i) {
        const auto& [send, peer] = group.operations[i];
        const auto& counterparts = send ? postedReceives_ : postedSends_;
        const auto key = send ? std::pair{group.process, peer} : std::pair{peer, group.process};
        const auto found = counterparts.find(key);
        if (found == counterparts.end() || found->second <= group.sequence[i]) {
          return false;
        }
      }
      return true;
    }
    default:
      return true;
    }
  }

  [[nodiscard]] long nextCorrection(const Actor& actor) const {
    return std::min(actor.corrections + actor.rate, actor.stepsUntilSync);
  }

  void run(const Item& item) {
    switch (item.kind) {
    case Kind::Record:
      completed_.insert(item.record);
      break;
    case Kind::Wait:
      break;
    case Kind::Work: {
      auto& actor = actors[item.actor];
      if (item.action == ActorAction::Predict) {
        // the neighbors (and the sends of the ghost clusters) are done with what gets overwritten
        for (const auto neighbor : actor.neighbors) {
          if (actor.predictions >= nextCorrection(actors[neighbor])) {
            ++violations;
          }
        }
        actor.predictions += actor.rate;
      } else {
        // the neighbors (and the receives of the ghost clusters) have provided what gets read
        for (const auto neighbor : actor.neighbors) {
          const auto& other = actors[neighbor];
          if (other.stepsUntilSync > other.predictions && actor.predictions > other.predictions) {
            ++violations;
          }
        }
        actor.corrections += actor.rate;
      }
      break;
    }
    case Kind::Group: {
      const auto& group = groups_[item.group];
      if (group.receiver != NoActor) {
        // the copy layer has read the ghost data that gets overwritten
        auto& ghost = actors[group.receiver];
        if (actors[ghost.copy].corrections < ghost.predictions) {
          ++violations;
        }
        ghost.predictions = group.receiveTarget;
      }
      if (group.sender != NoActor) {
        // the copy layer has written the data that goes out
        auto& ghost = actors[group.sender];
        const auto& copy = actors[ghost.copy];
        const auto lastPrediction = (copy.stepsUntilSync + copy.rate - 1) / copy.rate * copy.rate;
        if (copy.predictions < std::min(group.sendTarget, lastPrediction)) {
          ++violations;
        }
        ghost.corrections = group.sendTarget;
      }
      break;
    }
    }
  }

  std::map<std::size_t, std::deque<Item>> streams_;
  std::vector<Group> groups_;
  std::map<std::pair<int, int>, std::size_t> postedSends_;
  std::map<std::pair<int, int>, std::size_t> postedReceives_;
  std::map<std::size_t, int> processes_;
  std::map<int, std::size_t> recordings_;
  std::size_t lastRecording_{0};

  // the handles of the events (a deque keeps their addresses), the destroyed ones, and the ones
  // free to be handed out again, the one destroyed last at the back
  std::deque<char> handles_;
  std::set<Event> destroyed_;
  std::vector<Event> freed_;
  // the records: the latest one of each event, the recording each one belongs to, and the ones the
  // devices have reached
  std::size_t lastRecord_{0};
  std::map<Event, std::size_t> latestRecord_;
  std::map<std::size_t, std::size_t> recordedIn_;
  std::set<std::size_t> completed_;
};

/**
 * A copy layer that only enqueues its work on a stream of its own, with events from a pool as a
 * stream runtime has one.
 */
class StreamCopyCluster : public AbstractTimeCluster {
  public:
  StreamCopyCluster(long timeStepRate, std::size_t actor, StreamDevices& devices)
      : AbstractTimeCluster(static_cast<double>(timeStepRate), timeStepRate, Executor::Device),
        actor_(actor), devices_(devices), events_(devices.pool()) {
    setConcurrent(true);
  }

  /// as the cell clusters do in the end; returns the number of events still referenced
  std::size_t dispose() {
    finalize();
    return events_.dispose();
  }

  protected:
  void start() override {}
  void predict() override { action_ = ActorAction::Predict; }
  void correct() override { action_ = ActorAction::Correct; }
  void handleNeighborPrediction(const NeighborCluster& /*neighbor*/) override {}
  void handleNeighborCorrection(const NeighborCluster& /*neighbor*/) override {}
  void printTimeoutMessage(std::chrono::seconds /*timeSinceLastUpdate*/) override {}

  ActorEvent recordActionEvent() override {
    devices_.work(actor_, actor_, action_);
    auto event = events_.next();
    devices_.record(actor_, event.get());
    return ActorEvent(std::move(event));
  }

  void waitForEvent(void* event) override { devices_.wait(actor_, event); }

  private:
  std::size_t actor_;
  StreamDevices& devices_;
  parallel::runtime::EventPool events_;
  ActorAction action_{ActorAction::Nothing};
};

/**
 * How the exchange scheduler of a process deals with the events of its groups.
 */
enum class EventHandling {
  /// from a pool, and referenced by whoever gets them to wait for (as the stream exchange
  /// schedulers do)
  Pool,
  /// as the stream exchange scheduler did before: a new event for each group, handed on without
  /// any reference; at each synchronization point, all but the latest one get destroyed
  CreateAndDestroy
};

/**
 * The exchange scheduler of one process, in the global order and ordered on the device: one
 * communicator and one stream.
 */
class StreamScheduler : public ExchangeScheduler {
  public:
  StreamScheduler(StreamDevices& devices,
                  int process,
                  std::size_t clusterCount,
                  EventHandling handling)
      : ExchangeScheduler(clusterCount, LaunchOrder::Global), devices_(devices), process_(process),
        handling_(handling), events_(devices.pool()) {
    setStreamOrdered(true);
  }

  std::map<const ScheduledTransport*, std::size_t> ghosts;

  /**
   * At a synchronization point, once the devices have completed all work: with
   * EventHandling::CreateAndDestroy, destroys the events of all groups but the latest one.
   */
  void releaseEvents() {
    auto* latest = latestEvent().get();
    bool ownsLatest = false;
    for (auto* event : launched_) {
      if (event != latest) {
        devices_.destroy(event);
      } else {
        ownsLatest = true;
      }
    }
    launched_.clear();
    if (ownsLatest) {
      launched_.push_back(latest);
    }
  }

  /// as in the end; returns the number of events still referenced
  std::size_t dispose() {
    forgetEvents();
    for (auto* event : launched_) {
      devices_.destroy(event);
    }
    launched_.clear();
    return events_.dispose();
  }

  [[nodiscard]] std::size_t poolSize() const { return events_.size(); }

  protected:
  Ticket launch(std::size_t /*from*/,
                std::size_t /*to*/,
                std::size_t /*exchange*/,
                const ScheduledTransport* sender,
                const ScheduledTransport* receiver,
                const std::vector<void*>& after) override {
    const auto stream = StreamOffset + static_cast<std::size_t>(process_);
    for (auto* event : after) {
      devices_.wait(stream, event);
    }
    StreamDevices::Group group{process_, {}, {}, false, NoActor, 0, NoActor, 0};
    if (sender != nullptr) {
      for (const auto& region : sender->regions().copy) {
        group.operations.emplace_back(true, region.rank);
      }
      group.sender = ghosts.at(sender);
      group.sendTarget = ghostClusters->at(group.sender)->sendTarget();
    }
    if (receiver != nullptr) {
      for (const auto& region : receiver->regions().ghost) {
        group.operations.emplace_back(false, region.rank);
      }
      group.receiver = ghosts.at(receiver);
      group.receiveTarget = ghostClusters->at(group.receiver)->receiveTarget();
    }
    devices_.group(stream, std::move(group));
    if (handling_ == EventHandling::Pool) {
      auto event = events_.next();
      devices_.record(stream, event.get());
      setLatestEvent(ActorEvent(std::move(event)));
    } else {
      auto* event = devices_.create();
      devices_.record(stream, event);
      launched_.push_back(event);
      setLatestEvent(ActorEvent(event));
    }
    return ++tickets_;
  }

  bool completed(Ticket /*ticket*/) override { return true; }

  public:
  static constexpr std::size_t StreamOffset = 1000000;
  static constexpr std::size_t RecorderOffset = 2000000;
  const std::map<std::size_t, GhostCluster*>* ghostClusters{nullptr};

  private:
  StreamDevices& devices_;
  int process_;
  EventHandling handling_;
  parallel::runtime::EventPool events_;
  // with EventHandling::CreateAndDestroy: the events of the groups since the last release
  std::vector<StreamDevices::Event> launched_;
  Ticket tickets_{0};
};

struct Neighborhood {
  int process;
  int otherProcess;
  std::size_t cluster;
  std::size_t otherCluster;
};

/**
 * Records the work of one process into a graph, from any point on, as the `SuperStepRecorder` and
 * the `TimeManager` do on the device and with the scheduler; the replay follows right after the
 * recording. Its events come from a pool, as those of the `SuperStepRecorder` do.
 */
class ProcessRecorder {
  public:
  ProcessRecorder(StreamDevices& devices, int process, StreamScheduler& scheduler)
      : devices_(devices), process_(process), scheduler_(scheduler), events_(devices.pool()) {
    devices_.assign(recorderStream(), process);
    devices_.assign(schedulerStream(), process);
  }

  /// a copy layer of the process, on its stream
  void addCopy(std::size_t stream, AbstractTimeCluster& cluster) {
    devices_.assign(stream, process_);
    copyStreams_.push_back(stream);
    clusters_.push_back(&cluster);
  }

  void addGhost(AbstractTimeCluster& cluster) { clusters_.push_back(&cluster); }

  [[nodiscard]] bool recording() const { return recording_; }

  /// whether the replay still has to wait for the event
  [[nodiscard]] bool awaits(StreamDevices::Event event) const {
    return std::any_of(
        waitFor_.begin(), waitFor_.end(), [&](const auto& held) { return held.get() == event; });
  }

  void begin() {
    // The replay waits for the latest work of the clusters and the streams before the recording,
    // as noted now; the recorder keeps the events until then (SuperStepRecorder::beginReplay()).
    // The model runs the recorded work right away, so it enqueues the waits right away as well.
    waitFor_.clear();
    for (auto* cluster : clusters_) {
      waitFor_.push_back(cluster->latestEvent());
    }
    auto before = events_.next();
    devices_.record(schedulerStream(), before.get());
    waitFor_.emplace_back(std::move(before));
    for (const auto& event : waitFor_) {
      if (event) {
        devices_.wait(recorderStream(), event.get());
      }
    }

    devices_.beginRecording(process_);
    recording_ = true;
    fork_ = ActorEvent(events_.next());
    devices_.record(recorderStream(), fork_.get());
    for (auto* cluster : clusters_) {
      cluster->joinEvent(fork_.get());
      cluster->publishEvent(fork_);
    }
    devices_.wait(schedulerStream(), fork_.get());
    scheduler_.restartAfter(fork_);
  }

  void end() {
    for (const auto stream : copyStreams_) {
      join(stream);
    }
    join(schedulerStream());
    devices_.endRecording(process_);
    recording_ = false;

    // the replay has enqueued its waits
    waitFor_.clear();

    // everything enqueued from now on comes after the replay
    last_ = ActorEvent(events_.next());
    devices_.record(recorderStream(), last_.get());
    for (auto* cluster : clusters_) {
      cluster->joinEvent(last_.get());
      cluster->publishEvent(last_);
    }
    devices_.wait(schedulerStream(), last_.get());
    scheduler_.restartAfter(last_);
  }

  /// as in the end; returns the number of events still referenced
  std::size_t dispose() {
    waitFor_.clear();
    fork_ = ActorEvent();
    last_ = ActorEvent();
    return events_.dispose();
  }

  private:
  [[nodiscard]] std::size_t recorderStream() const {
    return StreamScheduler::RecorderOffset + static_cast<std::size_t>(process_);
  }
  [[nodiscard]] std::size_t schedulerStream() const {
    return StreamScheduler::StreamOffset + static_cast<std::size_t>(process_);
  }

  void join(std::size_t stream) {
    const auto event = events_.next();
    devices_.record(stream, event.get());
    devices_.wait(recorderStream(), event.get());
  }

  StreamDevices& devices_;
  int process_;
  StreamScheduler& scheduler_;
  parallel::runtime::EventPool events_;
  std::vector<std::size_t> copyStreams_;
  std::vector<AbstractTimeCluster*> clusters_;
  bool recording_{false};
  std::vector<ActorEvent> waitFor_;
  ActorEvent fork_;
  ActorEvent last_;
};

/**
 * What happened on the devices, summed over several runs.
 */
struct StreamOutcome {
  // runs in which the devices got stuck
  std::size_t stuck{0};
  std::size_t violations{0};
  std::size_t isolationViolations{0};
  std::size_t destroyedUses{0};
  std::size_t prematureRecords{0};
  // events still referenced when their pools were disposed of
  std::size_t leakedEvents{0};
  // the largest pool of events of a scheduler
  std::size_t schedulerPool{0};
};

/**
 * Enqueues all steps of all processes through the synchronization points, while the devices run
 * behind at random; with recordings, the processes also start and end recording their work at
 * random. At each synchronization point, the devices complete everything, and the schedulers
 * release their events as TimeManager::advanceInTime() did. In the end, everybody lets go of the
 * events as in TimeManager::freeDynamicResources().
 */
void runStreamOrdered(int processes,
                      std::size_t clusterCount,
                      const std::vector<Neighborhood>& neighborhoods,
                      const std::vector<long>& syncTimes,
                      std::mt19937& random,
                      StreamOutcome& outcome,
                      bool recordings = false,
                      EventHandling handling = EventHandling::Pool) {
  StreamDevices devices;
  const auto rate = [](std::size_t cluster) { return 1L << cluster; };

  std::vector<std::unique_ptr<StreamScheduler>> schedulers;
  for (int process = 0; process < processes; ++process) {
    schedulers.push_back(
        std::make_unique<StreamScheduler>(devices, process, clusterCount, handling));
  }

  std::map<std::pair<int, std::size_t>, std::size_t> copyIndex;
  std::map<std::tuple<int, std::size_t, std::size_t>, std::vector<int>> ghostPeers;
  for (const auto& neighborhood : neighborhoods) {
    for (const auto& [process, cluster, peer, otherCluster] :
         {std::tuple{neighborhood.process,
                     neighborhood.cluster,
                     neighborhood.otherProcess,
                     neighborhood.otherCluster},
          std::tuple{neighborhood.otherProcess,
                     neighborhood.otherCluster,
                     neighborhood.process,
                     neighborhood.cluster}}) {
      if (copyIndex.count({process, cluster}) == 0) {
        copyIndex[{process, cluster}] = devices.actors.size();
        devices.actors.push_back({rate(cluster), 0, NoActor, {}, 0, 0});
      }
      ghostPeers[{process, cluster, otherCluster}].push_back(peer);
    }
  }

  std::map<std::size_t, std::unique_ptr<StreamCopyCluster>> copies;
  for (const auto& [key, index] : copyIndex) {
    copies[index] = std::make_unique<StreamCopyCluster>(rate(key.second), index, devices);
  }
  for (const auto& [key, index] : copyIndex) {
    for (const auto& [otherKey, otherIndex] : copyIndex) {
      if (key.first == otherKey.first && key.second < otherKey.second &&
          otherKey.second - key.second <= 1) {
        copies[index]->connect(*copies[otherIndex]);
        devices.actors[index].neighbors.push_back(otherIndex);
        devices.actors[otherIndex].neighbors.push_back(index);
      }
    }
  }

  std::map<std::size_t, std::unique_ptr<GhostCluster>> ghosts;
  std::map<std::size_t, GhostCluster*> ghostPointers;
  std::map<std::size_t, int> ghostProcesses;
  for (auto& [key, peers] : ghostPeers) {
    const auto& [process, cluster, otherCluster] = key;
    std::sort(peers.begin(), peers.end());
    solver::RemoteClusterPair regions;
    for (const auto peer : peers) {
      regions.copy.emplace_back(nullptr, 1, RealType::F64, peer, 0);
      regions.ghost.emplace_back(nullptr, 1, RealType::F64, peer, 0);
    }
    const auto index = devices.actors.size();
    const auto copy = copyIndex.at({process, cluster});
    devices.actors.push_back({rate(otherCluster), 0, copy, {}, 0, 0});
    devices.actors[copy].neighbors.push_back(index);

    auto transport =
        std::make_unique<ScheduledTransport>(*schedulers[process], regions, cluster, otherCluster);
    schedulers[process]->ghosts[transport.get()] = index;
    auto ghost = std::make_unique<GhostCluster>(static_cast<double>(rate(otherCluster)),
                                                rate(otherCluster),
                                                "copy",
                                                "ghost",
                                                regions,
                                                std::move(transport));
    ghost->setConcurrent(true);
    copies.at(copy)->connect(*ghost);
    ghostPointers[index] = ghost.get();
    ghostProcesses[index] = process;
    ghosts[index] = std::move(ghost);
  }
  for (auto& scheduler : schedulers) {
    scheduler->ghostClusters = &ghostPointers;
  }

  std::vector<std::unique_ptr<ProcessRecorder>> recorders;
  recorders.reserve(static_cast<std::size_t>(processes));
  for (int process = 0; process < processes; ++process) {
    recorders.push_back(std::make_unique<ProcessRecorder>(devices, process, *schedulers[process]));
  }
  for (const auto& [key, index] : copyIndex) {
    recorders[static_cast<std::size_t>(key.first)]->addCopy(index, *copies.at(index));
  }
  for (auto& [index, ghost] : ghosts) {
    recorders[static_cast<std::size_t>(ghostProcesses.at(index))]->addGhost(*ghost);
  }

  std::vector<AbstractTimeCluster*> actors;
  for (auto& [index, copy] : copies) {
    actors.push_back(copy.get());
  }
  for (auto& [index, ghost] : ghosts) {
    actors.push_back(ghost.get());
  }

  // an event is still meant to be waited for while a cluster publishes it, or while a recorder
  // has yet to wait for it
  devices.awaited = [&](StreamDevices::Event event) {
    return std::any_of(actors.begin(),
                       actors.end(),
                       [&](const auto* actor) { return actor->latestEvent().get() == event; }) ||
           std::any_of(recorders.begin(), recorders.end(), [&](const auto& recorder) {
             return recorder->awaits(event);
           });
  };

  std::bernoulli_distribution deviceRuns(0.3);
  std::bernoulli_distribution toggleRecording(0.05);
  std::uniform_int_distribution<int> anyProcess(0, processes - 1);
  bool stuck = false;
  for (const auto syncTime : syncTimes) {
    if (!devices.drain(random)) {
      stuck = true;
      break;
    }
    if (handling == EventHandling::CreateAndDestroy) {
      for (auto& scheduler : schedulers) {
        scheduler->releaseEvents();
      }
    }
    devices.resetInterval();
    // the copy layers first: the ghost clusters announce their exchanges from them
    for (auto& [index, copy] : copies) {
      copy->setSyncTime(static_cast<double>(syncTime));
      copy->reset();
      devices.actors[index].stepsUntilSync = copy->getStepsUntilSync();
    }
    for (auto& [index, ghost] : ghosts) {
      ghost->setSyncTime(static_cast<double>(syncTime));
      ghost->reset();
      devices.actors[index].stepsUntilSync = ghost->getStepsUntilSync();
    }

    std::size_t rounds = 0;
    bool first = true;
    while (first || !std::all_of(actors.begin(), actors.end(), [](auto* actor) {
             return actor->synced();
           })) {
      first = false;
      // the host never waits for the devices
      REQUIRE(++rounds < 100000);
      std::shuffle(actors.begin(), actors.end(), random);
      for (auto* actor : actors) {
        actor->act();
        while (deviceRuns(random) && devices.runOne(random)) {
        }
        if (recordings && toggleRecording(random)) {
          auto& recorder = recorders[static_cast<std::size_t>(anyProcess(random))];
          if (recorder->recording()) {
            recorder->end();
          } else {
            recorder->begin();
          }
        }
      }
    }
    for (auto& recorder : recorders) {
      if (recorder->recording()) {
        recorder->end();
      }
    }
  }
  if (stuck || !devices.drain(random)) {
    ++outcome.stuck;
  }

  // the end
  for (auto& scheduler : schedulers) {
    scheduler->forgetEvents();
  }
  for (auto& [index, ghost] : ghosts) {
    ghost->finalize();
  }
  for (auto& [index, copy] : copies) {
    outcome.leakedEvents += copy->dispose();
  }
  for (auto& recorder : recorders) {
    outcome.leakedEvents += recorder->dispose();
  }
  for (auto& scheduler : schedulers) {
    outcome.schedulerPool = std::max(outcome.schedulerPool, scheduler->poolSize());
    outcome.leakedEvents += scheduler->dispose();
  }

  outcome.violations += devices.violations;
  outcome.isolationViolations += devices.isolationViolations;
  outcome.destroyedUses += devices.destroyedUses;
  outcome.prematureRecords += devices.prematureRecords;
}

std::vector<Neighborhood>
    randomNeighborhoods(int processes, std::size_t clusters, std::mt19937& random) {
  std::bernoulli_distribution neighboring(0.4);
  std::vector<Neighborhood> neighborhoods;
  for (int process = 0; process < processes; ++process) {
    for (int otherProcess = process + 1; otherProcess < processes; ++otherProcess) {
      for (std::size_t cluster = 0; cluster < clusters; ++cluster) {
        for (std::size_t otherCluster = cluster == 0 ? 0 : cluster - 1;
             otherCluster < std::min(cluster + 2, clusters);
             ++otherCluster) {
          if (neighboring(random)) {
            neighborhoods.push_back({process, otherProcess, cluster, otherCluster});
          }
        }
      }
    }
  }
  return neighborhoods;
}

} // namespace

TEST_CASE("Exchanges ordered on the device neither get stuck nor break the dependencies" *
          doctest::test_suite("solver")) {
  const std::vector<long> syncTimes{6, 13, 16, 29};
  std::mt19937 random(8642);
  StreamOutcome outcome;
  for (int trial = 0; trial < 150; ++trial) {
    const auto processes = std::uniform_int_distribution<int>(2, 4)(random);
    const auto clusters = std::uniform_int_distribution<std::size_t>(1, 3)(random);
    const auto neighborhoods = randomNeighborhoods(processes, clusters, random);
    runStreamOrdered(processes, clusters, neighborhoods, syncTimes, random, outcome);
  }
  CHECK(outcome.stuck == 0);
  CHECK(outcome.violations == 0);

  // the events stay reserved while anybody means to wait for them, also across the
  // synchronization points, and the scheduler needs no more of them than a pool has at first
  CHECK(outcome.destroyedUses == 0);
  CHECK(outcome.prematureRecords == 0);
  CHECK(outcome.leakedEvents == 0);
  CHECK(outcome.schedulerPool == parallel::runtime::EventPool::InitialSize);
}

TEST_CASE("Recordings of exchanges ordered on the device only wait for their own work" *
          doctest::test_suite("solver")) {
  // e.g. a ghost cluster may acknowledge a group launched before a recording only inside of it
  const std::vector<long> syncTimes{6, 13, 16, 29};
  std::mt19937 random(9753);
  StreamOutcome outcome;
  for (int trial = 0; trial < 150; ++trial) {
    const auto processes = std::uniform_int_distribution<int>(2, 4)(random);
    const auto clusters = std::uniform_int_distribution<std::size_t>(1, 3)(random);
    const auto neighborhoods = randomNeighborhoods(processes, clusters, random);
    runStreamOrdered(processes, clusters, neighborhoods, syncTimes, random, outcome, true);
  }
  CHECK(outcome.stuck == 0);
  CHECK(outcome.violations == 0);
  CHECK(outcome.isolationViolations == 0);

  // also a replay waits for the events it has noted before the recording
  CHECK(outcome.destroyedUses == 0);
  CHECK(outcome.prematureRecords == 0);
  CHECK(outcome.leakedEvents == 0);
  CHECK(outcome.schedulerPool == parallel::runtime::EventPool::InitialSize);
}

TEST_CASE("Exchange events that get destroyed at synchronization points are waited for too late" *
          doctest::test_suite("solver")) {
  // The scheduler as it was: its events were owned by the scheduler alone, and destroyed at the
  // synchronization points, while the ghost clusters still published them. After a
  // synchronization point, a copy layer (or a replay) then waited for a destroyed event, or for
  // one that had been created and recorded anew (CUDA hands out the handles of destroyed events
  // again). This makes sure that the model above sees such events.
  const std::vector<long> syncTimes{6, 13, 16, 29};
  for (const bool recordings : {false, true}) {
    CAPTURE(recordings);
    std::mt19937 random(2468);
    StreamOutcome outcome;
    for (int trial = 0; trial < 30; ++trial) {
      const auto processes = std::uniform_int_distribution<int>(2, 4)(random);
      const auto clusters = std::uniform_int_distribution<std::size_t>(1, 3)(random);
      const auto neighborhoods = randomNeighborhoods(processes, clusters, random);
      runStreamOrdered(processes,
                       clusters,
                       neighborhoods,
                       syncTimes,
                       random,
                       outcome,
                       recordings,
                       EventHandling::CreateAndDestroy);
    }
    MESSAGE("destroyed events used " << outcome.destroyedUses << " times, recorded too early "
                                     << outcome.prematureRecords << " times");
    CHECK(outcome.destroyedUses > 0);
    CHECK(outcome.prematureRecords > 0);
  }
}

} // namespace seissol::unit_test
