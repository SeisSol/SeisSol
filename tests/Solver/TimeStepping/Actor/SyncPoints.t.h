// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/Executor.h"
#include "Solver/TimeStepping/Actor/AbstractTimeCluster.h"
#include "Solver/TimeStepping/Actor/ActorState.h"
#include "Solver/TimeStepping/Plan/TimeSteppingPlan.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

namespace seissol::unit_test {
using namespace seissol::solver;

namespace {

/**
 * A cluster that checks the sizes of its steps.
 */
class SyncTestCluster : public AbstractTimeCluster {
  public:
  SyncTestCluster(double maxTimeStepSize, long timeStepRate, DataReadiness readiness)
      : AbstractTimeCluster(maxTimeStepSize, timeStepRate, Executor::Host), readiness_(readiness) {}

  [[nodiscard]] DataReadiness dataReadiness() const override { return readiness_; }

  [[nodiscard]] double correctionTime() const { return ct_.correctionTime; }

  //! steps that do not advance the time
  std::size_t emptySteps{0};

  //! steps before the last one of an interval that are shorter than the maximum
  std::size_t shortSteps{0};

  protected:
  void start() override {}
  void predict() override {
    const auto size = timeStepSize();
    emptySteps += size > 0 ? 0 : 1;
    shortSteps += ct_.isLastStep(ct_.stepsSinceLastSync) || size == ct_.maxTimeStepSize ? 0 : 1;
  }
  void correct() override {}
  void handleNeighborPrediction(const NeighborCluster& /*neighbor*/) override {}
  void handleNeighborCorrection(const NeighborCluster& /*neighbor*/) override {}
  void printTimeoutMessage(std::chrono::seconds /*timeSinceLastUpdate*/) override {}

  private:
  DataReadiness readiness_;
};

/**
 * The cell clusters (interior and copy) of `levels` time clusters, and face clusters on the two
 * smallest ones; connected as the time manager connects them. The time step sizes are those of
 * ClusterLadder: the update rate times the base time step.
 */
struct Ladder {
  std::vector<std::unique_ptr<SyncTestCluster>> clusters;

  Ladder(double baseTimeStepSize, long ratio, std::size_t levels) {
    std::vector<std::size_t> cells;
    std::vector<std::size_t> cellLevel;
    long rate = 1;
    for (std::size_t level = 0; level < levels; ++level) {
      const auto interior = add(baseTimeStepSize, rate, DataReadiness::AfterPrediction);
      const auto copy = add(baseTimeStepSize, rate, DataReadiness::AfterPrediction);
      clusters[copy]->setPriority(ActorPriority::High);
      for (const auto cell : {interior, copy}) {
        for (std::size_t i = 0; i < cells.size(); ++i) {
          if (level - cellLevel[i] <= 1) {
            clusters[cell]->connect(*clusters[cells[i]]);
          }
        }
        cells.push_back(cell);
        cellLevel.push_back(level);
      }
      if (level < 2) {
        const auto face = add(baseTimeStepSize, rate, DataReadiness::AfterCorrection);
        clusters[face]->connect(*clusters[interior]);
        clusters[face]->connect(*clusters[copy]);
      }
      rate *= ratio;
    }
  }

  std::size_t add(double baseTimeStepSize, long rate, DataReadiness readiness) {
    clusters.push_back(std::make_unique<SyncTestCluster>(
        static_cast<double>(rate) * baseTimeStepSize, rate, readiness));
    return clusters.size() - 1;
  }

  [[nodiscard]] bool synced() const {
    return std::all_of(clusters.begin(), clusters.end(), [](const auto& c) { return c->synced(); });
  }

  /**
   * Acts on any cluster that has an action to take, like the time manager does without a plan.
   * Returns false if the clusters get stuck before the synchronization point.
   */
  bool runFree() {
    for (auto& cluster : clusters) {
      REQUIRE(cluster->getNextLegalAction() == ActorAction::RestartAfterSync);
      cluster->act();
    }
    while (!synced()) {
      bool acted = false;
      for (auto& cluster : clusters) {
        if (!cluster->synced() && cluster->getNextLegalAction() != ActorAction::Nothing) {
          cluster->act();
          acted = true;
        }
      }
      if (!acted) {
        return false;
      }
    }
    return true;
  }

  /**
   * Takes the steps strictly along the time stepping plan. Returns false if a cluster is not
   * ready for its planned action.
   */
  bool runPlan() {
    std::vector<PlannedCluster> planned;
    for (auto& cluster : clusters) {
      REQUIRE(cluster->getNextLegalAction() == ActorAction::RestartAfterSync);
      cluster->act();
      planned.push_back({cluster->getTimeStepRate(),
                         cluster->getStepsUntilSync(),
                         cluster->dataReadiness(),
                         cluster->getPriority()});
    }
    for (const auto& action : planTimeSteps(planned)) {
      auto& cluster = *clusters[action.cluster];
      if (cluster.getNextLegalAction() != action.action) {
        return false;
      }
      cluster.act();
    }
    for (auto& cluster : clusters) {
      if (cluster->getNextLegalAction() != ActorAction::Sync) {
        return false;
      }
      cluster->act();
    }
    return true;
  }
};

struct Outcome {
  std::size_t intervals{0};
  std::size_t stuck{0};
  std::size_t disagreements{0};
  std::size_t emptySteps{0};
  std::size_t shortSteps{0};
  std::size_t offSyncPoint{0};
};

/**
 * Runs the clusters through the synchronization points of an output interval up to the end time;
 * the synchronization times accumulate as the modules accumulate them.
 */
Outcome runIntervals(double baseTimeStepSize,
                     long ratio,
                     std::size_t levels,
                     double interval,
                     double endTime,
                     bool plan) {
  Ladder ladder(baseTimeStepSize, ratio, levels);
  Outcome outcome;
  const double tolerance = ClusterTimes::TickTolerance * baseTimeStepSize;
  double current = 0;
  double next = interval;
  while (endTime > current + tolerance) {
    const double syncTime = std::min(endTime, next);
    for (auto& cluster : ladder.clusters) {
      cluster->setSyncTime(syncTime);
      // dereference first due to a clang-tidy recommendation
      (*cluster).reset();
    }
    const auto ticks = ladder.clusters.front()->getStepsUntilSync();
    for (const auto& cluster : ladder.clusters) {
      outcome.disagreements += cluster->getStepsUntilSync() == ticks ? 0 : 1;
    }

    ++outcome.intervals;
    if (!(plan ? ladder.runPlan() : ladder.runFree())) {
      ++outcome.stuck;
      break;
    }
    for (const auto& cluster : ladder.clusters) {
      outcome.offSyncPoint += cluster->correctionTime() == syncTime ? 0 : 1;
    }

    current = syncTime;
    if (std::abs(current - next) < tolerance) {
      next += interval;
    }
  }
  for (const auto& cluster : ladder.clusters) {
    outcome.emptySteps += cluster->emptySteps;
    outcome.shortSteps += cluster->shortSteps;
  }
  return outcome;
}

void checkIntervals(
    double baseTimeStepSize, long ratio, std::size_t levels, double interval, double endTime) {
  for (const bool plan : {false, true}) {
    CAPTURE(plan);
    const auto outcome = runIntervals(baseTimeStepSize, ratio, levels, interval, endTime, plan);
    CAPTURE(outcome.intervals);
    CHECK(outcome.stuck == 0);
    CHECK(outcome.disagreements == 0);
    CHECK(outcome.emptySteps == 0);
    CHECK(outcome.shortSteps == 0);
    CHECK(outcome.offSyncPoint == 0);
  }
}

// the base time steps of tpv5 in a Release and a Debug build
constexpr double BaseRelease = 1.4412823797093819e-03;
constexpr double BaseDebug = 1.4412823797093841e-03;

} // namespace

TEST_CASE("The clusters agree on the steps until a synchronization point at a multiple of the "
          "time step" *
          doctest::test_suite("solver")) {
  // Each of these intervals got some clusters to count one step more than their neighbors, which
  // then waited for each other forever.
  SUBCASE("10 base steps") { checkIntervals(BaseRelease, 2, 6, 0.01441282379709382, 0.1); }
  SUBCASE("32 base steps") { checkIntervals(BaseDebug, 2, 6, 0.04612103615070029, 0.35); }
  SUBCASE("27 base steps, rate 3") {
    checkIntervals(BaseRelease, 3, 4, 0.038914624252153314, 0.1);
    checkIntervals(BaseDebug, 3, 4, 0.03891462425215337, 0.25);
  }
}

TEST_CASE("The clusters agree on the steps until synchronization points" *
          doctest::test_suite("solver")) {
  for (const auto base : {BaseRelease, BaseDebug}) {
    for (long multiple = 1; multiple <= 24; ++multiple) {
      for (const long unit : {1L, 32L}) {
        CAPTURE(base);
        CAPTURE(multiple);
        CAPTURE(unit);
        const double interval = static_cast<double>(multiple * unit) * base;
        checkIntervals(base, 2, 6, interval, 6 * interval);
      }
    }
  }

  SUBCASE("A remainder within the tolerance does not add a step") {
    checkIntervals(BaseRelease, 2, 4, (40 + 1e-7) * BaseRelease, 0.5);
    checkIntervals(BaseRelease, 2, 4, (40 - 1e-7) * BaseRelease, 0.5);
  }

  SUBCASE("A remainder beyond the tolerance adds a step") {
    checkIntervals(BaseRelease, 2, 4, (40 + 1e-3) * BaseRelease, 0.5);
    checkIntervals(BaseRelease, 2, 4, 0.1, 1.0);
  }
}

} // namespace seissol::unit_test
