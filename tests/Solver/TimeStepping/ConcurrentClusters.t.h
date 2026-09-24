// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/Executor.h"
#include "Solver/TimeStepping/AbstractTimeCluster.h"
#include "Solver/TimeStepping/ActorState.h"
#include "Solver/TimeStepping/TimeSteppingPlan.h"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <memory>
#include <random>
#include <set>
#include <vector>

namespace seissol::unit_test {
using namespace seissol::time_stepping;

namespace {

/**
 * A device with one stream per cluster. The clusters only enqueue their actions; the device runs
 * them later, in any order that keeps the order of each stream and the event waits. Whenever an
 * action runs, the rules of the actor model are checked against what the device has run so far.
 */
class SimulatedDevice {
  public:
  struct Cluster {
    long rate;
    long stepsUntilSync{0};
    DataReadiness readiness;
    std::vector<std::size_t> neighbors;
    long predictions{0};
    long corrections{0};
  };

  std::vector<Cluster> clusters;
  std::size_t violations{0};

  void enqueue(std::size_t cluster,
               ActorAction action,
               std::vector<std::uintptr_t> waits,
               std::uintptr_t event) {
    queues_.resize(clusters.size());
    queues_[cluster].push_back({action, std::move(waits), event});
  }

  /**
   * Runs one action whose waits have completed, chosen at random. Returns false if there is none.
   */
  bool runOne(std::mt19937& random) {
    std::vector<std::size_t> runnable;
    for (std::size_t cluster = 0; cluster < queues_.size(); ++cluster) {
      if (!queues_[cluster].empty() &&
          std::all_of(queues_[cluster].front().waits.begin(),
                      queues_[cluster].front().waits.end(),
                      [&](auto event) { return completed_.count(event) > 0; })) {
        runnable.push_back(cluster);
      }
    }
    if (runnable.empty()) {
      return false;
    }
    const auto cluster =
        runnable[std::uniform_int_distribution<std::size_t>(0, runnable.size() - 1)(random)];
    const auto item = queues_[cluster].front();
    queues_[cluster].pop_front();
    run(cluster, item.action);
    completed_.insert(item.event);
    return true;
  }

  void drain(std::mt19937& random) {
    while (runOne(random)) {
    }
  }

  void resetInterval() {
    for (auto& cluster : clusters) {
      cluster.predictions = 0;
      cluster.corrections = 0;
    }
  }

  private:
  struct Item {
    ActorAction action;
    std::vector<std::uintptr_t> waits;
    std::uintptr_t event;
  };

  void run(std::size_t index, ActorAction action) {
    auto& cluster = clusters[index];
    if (action == ActorAction::Predict) {
      // the neighbors must have read what this prediction overwrites
      for (const auto neighbor : cluster.neighbors) {
        const auto& other = clusters[neighbor];
        const auto next = std::min(other.corrections + other.rate, other.stepsUntilSync);
        if (cluster.predictions >= next) {
          ++violations;
        }
      }
      cluster.predictions += cluster.rate;
    } else {
      // the neighbors must have provided what this correction reads
      for (const auto neighbor : cluster.neighbors) {
        const auto& other = clusters[neighbor];
        const auto provided = other.readiness == DataReadiness::AfterPrediction ? other.predictions
                                                                                : other.corrections;
        if (other.stepsUntilSync > provided && cluster.predictions > provided) {
          ++violations;
        }
      }
      cluster.corrections += cluster.rate;
    }
  }

  std::vector<std::deque<Item>> queues_;
  std::set<std::uintptr_t> completed_;
};

class EnqueueingCluster : public AbstractTimeCluster {
  public:
  EnqueueingCluster(std::size_t index,
                    long timeStepRate,
                    DataReadiness readiness,
                    SimulatedDevice& device,
                    std::uintptr_t& nextEvent)
      : AbstractTimeCluster(static_cast<double>(timeStepRate), timeStepRate, Executor::Device),
        index_(index), readiness_(readiness), device_(device), nextEvent_(nextEvent) {
    setConcurrent(true);
  }

  [[nodiscard]] DataReadiness dataReadiness() const override { return readiness_; }

  protected:
  void start() override {}
  void predict() override { action_ = ActorAction::Predict; }
  void correct() override { action_ = ActorAction::Correct; }
  void handleNeighborPrediction(const NeighborCluster& /*neighbor*/) override {}
  void handleNeighborCorrection(const NeighborCluster& /*neighbor*/) override {}
  void printTimeoutMessage(std::chrono::seconds /*timeSinceLastUpdate*/) override {}

  void* recordActionEvent() override {
    const auto event = nextEvent_++;
    device_.enqueue(index_, action_, waits_, event);
    waits_.clear();
    return reinterpret_cast<void*>(event);
  }

  void waitForEvent(void* event) override {
    waits_.push_back(reinterpret_cast<std::uintptr_t>(event));
  }

  private:
  std::size_t index_;
  DataReadiness readiness_;
  SimulatedDevice& device_;
  std::uintptr_t& nextEvent_;
  ActorAction action_{ActorAction::Nothing};
  std::vector<std::uintptr_t> waits_;
};

/**
 * Cell and face clusters of one process, connected as the time manager connects them.
 */
struct ConcurrentLayout {
  SimulatedDevice device;
  std::uintptr_t nextEvent{1};
  std::vector<std::unique_ptr<EnqueueingCluster>> clusters;

  ConcurrentLayout(std::size_t clusterCount, long ratio, std::mt19937& random) {
    std::bernoulli_distribution withFaces(0.6);
    std::vector<std::size_t> cells;
    std::vector<long> cellCluster;
    long rate = 1;
    for (std::size_t cluster = 0; cluster < clusterCount; ++cluster) {
      const auto interior = add(rate, DataReadiness::AfterPrediction, ActorPriority::Low);
      const auto copy = add(rate, DataReadiness::AfterPrediction, ActorPriority::High);
      for (const auto cell : {interior, copy}) {
        for (std::size_t i = 0; i < cells.size(); ++i) {
          if (static_cast<long>(cluster) - cellCluster[i] <= 1) {
            connect(cell, cells[i]);
          }
        }
        cells.push_back(cell);
        cellCluster.push_back(static_cast<long>(cluster));
      }
      if (withFaces(random)) {
        const auto face = add(rate, DataReadiness::AfterCorrection, ActorPriority::Low);
        connect(face, interior);
        connect(face, copy);
      }
      if (withFaces(random)) {
        const auto face = add(rate, DataReadiness::AfterCorrection, ActorPriority::High);
        connect(face, copy);
      }
      rate *= ratio;
    }
  }

  std::size_t add(long rate, DataReadiness readiness, ActorPriority priority) {
    const auto index = clusters.size();
    clusters.push_back(
        std::make_unique<EnqueueingCluster>(index, rate, readiness, device, nextEvent));
    clusters.back()->setPriority(priority);
    device.clusters.push_back({rate, 0, readiness, {}});
    return index;
  }

  void connect(std::size_t a, std::size_t b) {
    clusters[a]->connect(*clusters[b]);
    device.clusters[a].neighbors.push_back(b);
    device.clusters[b].neighbors.push_back(a);
  }
};

/**
 * Enqueues all steps up to each synchronization point, either along the plan or whenever a cluster
 * is ready (in random order), while the simulated device runs behind at random. Returns the number
 * of actions the device ran before their data was ready or while it was still being read.
 */
std::size_t runConcurrently(ConcurrentLayout& layout,
                            const std::vector<long>& syncTimes,
                            bool alongPlan,
                            std::mt19937& random) {
  std::bernoulli_distribution deviceRuns(0.5);
  auto& device = layout.device;
  for (const auto syncTime : syncTimes) {
    for (auto& cluster : layout.clusters) {
      cluster->setSyncTime(static_cast<double>(syncTime));
      cluster->reset();
    }
    device.resetInterval();
    std::vector<PlannedCluster> planned;
    for (std::size_t i = 0; i < layout.clusters.size(); ++i) {
      auto& cluster = *layout.clusters[i];
      REQUIRE(cluster.getNextLegalAction() == ActorAction::RestartAfterSync);
      cluster.act();
      device.clusters[i].stepsUntilSync = cluster.getStepsUntilSync();
      planned.push_back({cluster.getTimeStepRate(),
                         cluster.getStepsUntilSync(),
                         cluster.dataReadiness(),
                         cluster.getPriority()});
    }

    if (alongPlan) {
      for (const auto& step : planTimeSteps(planned)) {
        auto& cluster = *layout.clusters[step.cluster];
        REQUIRE(cluster.getNextLegalAction() == step.action);
        cluster.act();
        while (deviceRuns(random) && device.runOne(random)) {
        }
      }
    } else {
      std::vector<std::size_t> order(layout.clusters.size());
      for (std::size_t i = 0; i < order.size(); ++i) {
        order[i] = i;
      }
      bool acted = true;
      while (acted) {
        acted = false;
        std::shuffle(order.begin(), order.end(), random);
        for (const auto i : order) {
          auto& cluster = *layout.clusters[i];
          const auto action = cluster.getNextLegalAction();
          if (action == ActorAction::Predict || action == ActorAction::Correct) {
            cluster.act();
            acted = true;
          }
          while (deviceRuns(random) && device.runOne(random)) {
          }
        }
      }
    }

    for (auto& cluster : layout.clusters) {
      REQUIRE(cluster->getNextLegalAction() == ActorAction::Sync);
      cluster->act();
    }
    // the synchronization point waits for the device
    device.drain(random);
  }
  return device.violations;
}

} // namespace

TEST_CASE("Concurrent clusters wait for their neighbors on the device" *
          doctest::test_suite("solver")) {
  std::mt19937 random(2468);
  for (const auto alongPlan : {true, false}) {
    std::size_t violations = 0;
    for (int trial = 0; trial < 200; ++trial) {
      std::uniform_int_distribution<std::size_t> clusterCount(1, 4);
      std::uniform_int_distribution<long> ratio(2, 3);
      std::uniform_int_distribution<long> interval(1, 30);
      ConcurrentLayout layout(clusterCount(random), ratio(random), random);
      std::vector<long> syncTimes;
      long time = 0;
      for (int i = 0; i < 3; ++i) {
        time += interval(random);
        syncTimes.push_back(time);
      }
      CAPTURE(alongPlan);
      CAPTURE(trial);
      const auto trialViolations = runConcurrently(layout, syncTimes, alongPlan, random);
      CHECK(trialViolations == 0);
      violations += trialViolations;
    }
    CHECK(violations == 0);
  }
}

} // namespace seissol::unit_test
