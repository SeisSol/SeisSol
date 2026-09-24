// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/Executor.h"
#include "Kernels/Precision.h"
#include "Solver/TimeStepping/AbstractTimeCluster.h"
#include "Solver/TimeStepping/ActorState.h"
#include "Solver/TimeStepping/ExchangeScheduler.h"
#include "Solver/TimeStepping/GhostCluster.h"
#include "Solver/TimeStepping/HaloCommunication.h"
#include "Solver/TimeStepping/HaloTransport.h"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <deque>
#include <functional>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

namespace seissol::unit_test {
using namespace seissol::time_stepping;

namespace {

/**
 * The point-to-point semantics of NCCL and its relatives, for a number of processes, without any
 * buffering: each stream of a process runs its groups of operations one after another, a send
 * completes once the peer has posted the matching receive and vice versa, and operations match in
 * the order they are posted, per communicator, sender and receiver.
 */
class RendezvousNetwork {
  public:
  struct Operation {
    bool send;
    int peer;
    std::size_t communicator;
  };

  std::size_t enqueue(int process, std::size_t stream, std::vector<Operation> operations) {
    const auto ticket = nextTicket_++;
    queues_[{process, stream}].push_back(Group{std::move(operations), ticket, {}, false});
    return ticket;
  }

  [[nodiscard]] bool completed(std::size_t ticket) const {
    return completedTickets_.count(ticket) > 0;
  }

  [[nodiscard]] bool idle() const {
    return std::all_of(
        queues_.begin(), queues_.end(), [](const auto& queue) { return queue.second.empty(); });
  }

  /**
   * Posts the groups that are first in line and removes the completed ones. Returns whether
   * anything happened.
   */
  bool progress() {
    bool changed = false;
    bool again = true;
    while (again) {
      again = false;
      for (auto& [key, queue] : queues_) {
        const auto process = key.first;
        if (queue.empty()) {
          continue;
        }
        auto& group = queue.front();
        if (!group.posted) {
          post(process, group);
          changed = again = true;
        }
        if (done(process, group)) {
          completedTickets_.insert(group.ticket);
          queue.pop_front();
          changed = again = true;
        }
      }
    }
    return changed;
  }

  private:
  struct Group {
    std::vector<Operation> operations;
    std::size_t ticket;
    std::vector<std::size_t> sequence;
    bool posted;
  };

  using Channel = std::tuple<std::size_t, int, int>;

  static Channel channel(int process, const Operation& operation) {
    return operation.send ? Channel{operation.communicator, process, operation.peer}
                          : Channel{operation.communicator, operation.peer, process};
  }

  void post(int process, Group& group) {
    for (const auto& operation : group.operations) {
      auto& counter =
          (operation.send ? postedSends_ : postedReceives_)[channel(process, operation)];
      group.sequence.push_back(counter++);
    }
    group.posted = true;
  }

  [[nodiscard]] bool done(int process, const Group& group) const {
    for (std::size_t i = 0; i < group.operations.size(); ++i) {
      const auto& operation = group.operations[i];
      const auto& counterparts = operation.send ? postedReceives_ : postedSends_;
      const auto counterpart = counterparts.find(channel(process, operation));
      if (counterpart == counterparts.end() || counterpart->second <= group.sequence[i]) {
        return false;
      }
    }
    return true;
  }

  std::map<std::pair<int, std::size_t>, std::deque<Group>> queues_;
  std::map<Channel, std::size_t> postedSends_;
  std::map<Channel, std::size_t> postedReceives_;
  std::set<std::size_t> completedTickets_;
  std::size_t nextTicket_{1};
};

std::vector<RendezvousNetwork::Operation> operations(bool sends,
                                                     bool receives,
                                                     const std::vector<int>& peers,
                                                     std::size_t sendCommunicator,
                                                     std::size_t receiveCommunicator) {
  std::vector<RendezvousNetwork::Operation> result;
  for (const auto peer : peers) {
    if (sends) {
      result.push_back({true, peer, sendCommunicator});
    }
  }
  for (const auto peer : peers) {
    if (receives) {
      result.push_back({false, peer, receiveCommunicator});
    }
  }
  return result;
}

/**
 * How the simulated scheduler of a process launches its groups.
 */
enum class Launching {
  /// one communicator and one stream per direction
  PerDirection,
  /// one communicator and one stream, in the global order
  Global,
  /// one communicator and one stream, but each direction whenever it is ready
  SingleStreamPerDirection
};

/**
 * The exchange scheduler of one process on the simulated network.
 */
class SimulatedScheduler : public ExchangeScheduler {
  public:
  SimulatedScheduler(RendezvousNetwork& network,
                     int process,
                     std::size_t clusterCount,
                     Launching launching)
      : ExchangeScheduler(clusterCount,
                          launching == Launching::Global ? LaunchOrder::Global
                                                         : LaunchOrder::PerDirection),
        network_(network), process_(process), clusterCount_(clusterCount),
        single_(launching != Launching::PerDirection) {}

  protected:
  Ticket launch(std::size_t from,
                std::size_t to,
                const ScheduledTransport* sender,
                const ScheduledTransport* receiver) override {
    const auto communicator = single_ ? 0 : from * clusterCount_ + to;
    std::vector<RendezvousNetwork::Operation> group;
    if (sender != nullptr) {
      for (const auto& region : sender->regions().copy) {
        group.push_back({true, region.rank, communicator});
      }
    }
    if (receiver != nullptr) {
      for (const auto& region : receiver->regions().ghost) {
        group.push_back({false, region.rank, communicator});
      }
    }
    return network_.enqueue(process_, communicator, std::move(group));
  }

  bool completed(Ticket ticket) override { return network_.completed(ticket); }

  private:
  RendezvousNetwork& network_;
  int process_;
  std::size_t clusterCount_;
  bool single_;
};

/**
 * A transport on the simulated network that launches everything right away.
 */
class ImmediateTransport : public HaloTransport {
  public:
  ImmediateTransport(RendezvousNetwork& network,
                     int process,
                     std::size_t cluster,
                     std::size_t otherCluster,
                     std::size_t clusterCount,
                     std::vector<int> peers)
      : network_(network), process_(process), peers_(std::move(peers)),
        sendCommunicator_(cluster * clusterCount + otherCluster),
        receiveCommunicator_(otherCluster * clusterCount + cluster) {}

  void startSend() override {
    sendTicket_ = network_.enqueue(
        process_, 0, operations(true, false, peers_, sendCommunicator_, receiveCommunicator_));
  }
  bool testSend() override { return sendTicket_ == 0 || network_.completed(sendTicket_); }
  void startReceive() override {
    receiveTicket_ = network_.enqueue(
        process_, 0, operations(false, true, peers_, sendCommunicator_, receiveCommunicator_));
  }
  bool testReceive() override { return receiveTicket_ == 0 || network_.completed(receiveTicket_); }

  private:
  RendezvousNetwork& network_;
  int process_;
  std::vector<int> peers_;
  std::size_t sendCommunicator_;
  std::size_t receiveCommunicator_;
  std::size_t sendTicket_{0};
  std::size_t receiveTicket_{0};
};

class SimulatedCopyCluster : public AbstractTimeCluster {
  public:
  explicit SimulatedCopyCluster(long timeStepRate)
      : AbstractTimeCluster(static_cast<double>(timeStepRate), timeStepRate, Executor::Host) {}

  protected:
  void start() override {}
  void predict() override {}
  void correct() override {}
  void handleNeighborPrediction(const NeighborCluster& /*neighbor*/) override {}
  void handleNeighborCorrection(const NeighborCluster& /*neighbor*/) override {}
  void printTimeoutMessage(std::chrono::seconds /*timeSinceLastUpdate*/) override {}
};

/// Cells of `cluster` on `process` lie next to cells of `otherCluster` on `otherProcess`.
struct Adjacency {
  int process;
  int otherProcess;
  std::size_t cluster;
  std::size_t otherCluster;
};

/**
 * Runs the copy and ghost clusters of all processes through the synchronization points. Returns
 * false if they get stuck.
 */
enum class Transport { Scheduled, Immediate };

bool simulate(Transport transport,
              Launching launching,
              int processes,
              std::size_t clusterCount,
              const std::vector<Adjacency>& adjacencies,
              const std::vector<long>& syncTimes,
              std::mt19937& random) {
  RendezvousNetwork network;
  std::vector<std::unique_ptr<SimulatedScheduler>> schedulers;
  for (int process = 0; process < processes; ++process) {
    schedulers.push_back(
        std::make_unique<SimulatedScheduler>(network, process, clusterCount, launching));
  }
  const auto makeTransport = [&](int process,
                                 std::size_t cluster,
                                 std::size_t otherCluster,
                                 const solver::RemoteClusterPair& regions,
                                 const std::vector<int>& peers) -> std::unique_ptr<HaloTransport> {
    if (transport == Transport::Scheduled) {
      return std::make_unique<ScheduledTransport>(
          *schedulers[process], regions, cluster, otherCluster);
    }
    return std::make_unique<ImmediateTransport>(
        network, process, cluster, otherCluster, clusterCount, peers);
  };
  const auto rate = [](std::size_t cluster) { return 1L << cluster; };

  std::map<std::pair<int, std::size_t>, std::unique_ptr<SimulatedCopyCluster>> copies;
  std::map<std::tuple<int, std::size_t, std::size_t>, std::vector<int>> ghostPeers;
  for (const auto& adjacency : adjacencies) {
    for (const auto& [process, cluster, peer, otherCluster] :
         {std::tuple{
              adjacency.process, adjacency.cluster, adjacency.otherProcess, adjacency.otherCluster},
          std::tuple{adjacency.otherProcess,
                     adjacency.otherCluster,
                     adjacency.process,
                     adjacency.cluster}}) {
      auto& copy = copies[{process, cluster}];
      if (copy == nullptr) {
        copy = std::make_unique<SimulatedCopyCluster>(rate(cluster));
      }
      ghostPeers[{process, cluster, otherCluster}].push_back(peer);
    }
  }

  for (auto& [key, copy] : copies) {
    for (auto& [otherKey, otherCopy] : copies) {
      if (key.first == otherKey.first && key.second < otherKey.second &&
          otherKey.second - key.second <= 1) {
        copy->connect(*otherCopy);
      }
    }
  }

  std::vector<std::unique_ptr<GhostCluster>> ghosts;
  for (auto& [key, peers] : ghostPeers) {
    const auto& [process, cluster, otherCluster] = key;
    std::sort(peers.begin(), peers.end());
    solver::RemoteClusterPair regions;
    for (const auto peer : peers) {
      regions.copy.emplace_back(nullptr, 1, RealType::F64, peer, 0);
      regions.ghost.emplace_back(nullptr, 1, RealType::F64, peer, 0);
    }
    ghosts.push_back(std::make_unique<GhostCluster>(
        static_cast<double>(rate(otherCluster)),
        rate(otherCluster),
        "copy",
        "ghost",
        regions,
        makeTransport(process, cluster, otherCluster, regions, peers)));
    copies.at({process, cluster})->connect(*ghosts.back());
  }

  std::vector<AbstractTimeCluster*> actors;
  for (auto& [key, copy] : copies) {
    actors.push_back(copy.get());
  }
  for (auto& ghost : ghosts) {
    actors.push_back(ghost.get());
  }

  for (const auto syncTime : syncTimes) {
    // the copy layers first: the ghost clusters announce their exchanges from them
    for (auto& [key, copy] : copies) {
      copy->setSyncTime(static_cast<double>(syncTime));
      copy->reset();
    }
    for (auto& ghost : ghosts) {
      ghost->setSyncTime(static_cast<double>(syncTime));
      ghost->reset();
    }
    std::size_t idleRounds = 0;
    bool first = true;
    while (first || !std::all_of(actors.begin(), actors.end(), [](auto* actor) {
             return actor->synced();
           })) {
      first = false;
      std::shuffle(actors.begin(), actors.end(), random);
      bool changed = false;
      for (auto* actor : actors) {
        changed = actor->act().isStateChanged || changed;
        changed = network.progress() || changed;
      }
      idleRounds = changed ? 0 : idleRounds + 1;
      if (idleRounds > 10) {
        return false;
      }
    }
    if (!network.idle()) {
      return false;
    }
  }
  return true;
}

std::vector<Adjacency>
    randomAdjacencies(int processes, std::size_t clusterCount, std::mt19937& random) {
  std::bernoulli_distribution adjacent(0.4);
  std::vector<Adjacency> adjacencies;
  for (int process = 0; process < processes; ++process) {
    for (int otherProcess = process + 1; otherProcess < processes; ++otherProcess) {
      for (std::size_t cluster = 0; cluster < clusterCount; ++cluster) {
        for (std::size_t otherCluster = cluster == 0 ? 0 : cluster - 1;
             otherCluster < std::min(cluster + 2, clusterCount);
             ++otherCluster) {
          if (adjacent(random)) {
            adjacencies.push_back({process, otherProcess, cluster, otherCluster});
          }
        }
      }
    }
  }
  return adjacencies;
}

} // namespace

TEST_CASE("Blocking point-to-point transports do not get stuck" * doctest::test_suite("solver")) {
  // synchronization points on and between the steps of the larger clusters
  const std::vector<long> syncTimes{6, 13, 16, 29};
  for (const auto launching : {Launching::PerDirection, Launching::Global}) {
    std::mt19937 random(1234);
    for (int trial = 0; trial < 200; ++trial) {
      std::uniform_int_distribution<int> processCount(2, 5);
      std::uniform_int_distribution<std::size_t> clusterCount(1, 3);
      const auto processes = processCount(random);
      const auto clusters = clusterCount(random);
      const auto adjacencies = randomAdjacencies(processes, clusters, random);
      CAPTURE(launching == Launching::Global);
      CAPTURE(trial);
      CHECK(simulate(
          Transport::Scheduled, launching, processes, clusters, adjacencies, syncTimes, random));
    }
  }
}

TEST_CASE("One stream without the global order gets stuck" * doctest::test_suite("solver")) {
  const std::vector<long> syncTimes{6, 13, 16, 29};
  std::mt19937 random(1234);
  int stuck = 0;
  for (int trial = 0; trial < 200; ++trial) {
    std::uniform_int_distribution<int> processCount(2, 5);
    std::uniform_int_distribution<std::size_t> clusterCount(1, 3);
    const auto processes = processCount(random);
    const auto clusters = clusterCount(random);
    const auto adjacencies = randomAdjacencies(processes, clusters, random);
    if (!simulate(Transport::Scheduled,
                  Launching::SingleStreamPerDirection,
                  processes,
                  clusters,
                  adjacencies,
                  syncTimes,
                  random)) {
      ++stuck;
    }
  }
  MESSAGE("stuck in " << stuck << " of 200 layouts");
  CHECK(stuck > 0);
}

TEST_CASE("Launching receives right away gets stuck" * doctest::test_suite("solver")) {
  // two processes exchanging within one time cluster: both queue their receive first
  std::mt19937 random(1);
  CHECK_FALSE(
      simulate(Transport::Immediate, Launching::PerDirection, 2, 1, {{0, 1, 0, 0}}, {8}, random));
  CHECK(simulate(Transport::Scheduled, Launching::PerDirection, 2, 1, {{0, 1, 0, 0}}, {8}, random));
  CHECK(simulate(Transport::Scheduled, Launching::Global, 2, 1, {{0, 1, 0, 0}}, {8}, random));
}

} // namespace seissol::unit_test
