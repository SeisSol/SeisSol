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

#include <chrono>
#include <cstddef>
#include <memory>
#include <random>
#include <tuple>
#include <vector>

namespace seissol::unit_test {
using namespace seissol::time_stepping;

namespace {

class PlanTestCluster : public AbstractTimeCluster {
  public:
  PlanTestCluster(long timeStepRate, DataReadiness readiness)
      : AbstractTimeCluster(static_cast<double>(timeStepRate), timeStepRate, Executor::Host),
        readiness_(readiness) {}

  [[nodiscard]] DataReadiness dataReadiness() const override { return readiness_; }

  protected:
  void start() override {}
  void predict() override {}
  void correct() override {}
  void handleNeighborPrediction(const NeighborCluster& /*neighbor*/) override {}
  void handleNeighborCorrection(const NeighborCluster& /*neighbor*/) override {}
  void printTimeoutMessage(std::chrono::seconds /*timeSinceLastUpdate*/) override {}

  private:
  DataReadiness readiness_;
};

/**
 * The cell and face clusters of one process, connected as the time manager connects them.
 */
struct Layout {
  std::vector<std::unique_ptr<PlanTestCluster>> clusters;
  std::vector<ActorPriority> priorities;

  Layout(std::size_t clusterCount, long ratio, std::mt19937& random) {
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
            clusters[cell]->connect(*clusters[cells[i]]);
          }
        }
        cells.push_back(cell);
        cellCluster.push_back(static_cast<long>(cluster));
      }
      if (withFaces(random)) {
        const auto face = add(rate, DataReadiness::AfterCorrection, ActorPriority::Low);
        clusters[face]->connect(*clusters[interior]);
        clusters[face]->connect(*clusters[copy]);
      }
      if (withFaces(random)) {
        const auto face = add(rate, DataReadiness::AfterCorrection, ActorPriority::High);
        clusters[face]->connect(*clusters[copy]);
      }
      rate *= ratio;
    }
  }

  std::size_t add(long rate, DataReadiness readiness, ActorPriority priority) {
    clusters.push_back(std::make_unique<PlanTestCluster>(rate, readiness));
    clusters.back()->setPriority(priority);
    priorities.push_back(priority);
    return clusters.size() - 1;
  }
};

/**
 * Takes the clusters through the synchronization points strictly along the plan. Returns the
 * number of planned actions the clusters were not ready for.
 */
std::size_t followPlan(Layout& layout, const std::vector<long>& syncTimes) {
  std::size_t notReady = 0;
  for (const auto syncTime : syncTimes) {
    std::vector<PlannedCluster> planned;
    for (auto& cluster : layout.clusters) {
      cluster->setSyncTime(static_cast<double>(syncTime));
      cluster->reset();
    }
    for (auto& cluster : layout.clusters) {
      REQUIRE(cluster->getNextLegalAction() == ActorAction::RestartAfterSync);
      cluster->act();
      planned.push_back({cluster->getTimeStepRate(),
                         cluster->getStepsUntilSync(),
                         cluster->dataReadiness(),
                         cluster->getPriority()});
    }

    for (const auto& action : planTimeSteps(planned)) {
      auto& cluster = *layout.clusters[action.cluster];
      if (cluster.getNextLegalAction() != action.action) {
        ++notReady;
        return notReady;
      }
      cluster.act();
    }

    for (auto& cluster : layout.clusters) {
      // after the plan, nothing is left but the synchronization
      CHECK(cluster->getNextLegalAction() == ActorAction::Sync);
      cluster->act();
      CHECK(cluster->synced());
    }
  }
  return notReady;
}

} // namespace

TEST_CASE("The time stepping plan is a legal order of all steps" * doctest::test_suite("solver")) {
  std::mt19937 random(4321);
  std::size_t notReady = 0;
  for (int trial = 0; trial < 300; ++trial) {
    std::uniform_int_distribution<std::size_t> clusterCount(1, 4);
    std::uniform_int_distribution<long> ratio(2, 3);
    std::uniform_int_distribution<long> interval(1, 30);
    Layout layout(clusterCount(random), ratio(random), random);
    std::vector<long> syncTimes;
    long time = 0;
    for (int i = 0; i < 4; ++i) {
      time += interval(random);
      syncTimes.push_back(time);
    }
    CAPTURE(trial);
    const auto trialNotReady = followPlan(layout, syncTimes);
    CHECK(trialNotReady == 0);
    notReady += trialNotReady;
  }
  CHECK(notReady == 0);
}

TEST_CASE("The time stepping plan orders a pair of clusters by logical time" *
          doctest::test_suite("solver")) {
  // a cluster with twice the time step next to the smallest one, up to 3 steps of the smallest
  const auto plan = planTimeSteps({{1, 3, DataReadiness::AfterPrediction, ActorPriority::Low},
                                   {2, 3, DataReadiness::AfterPrediction, ActorPriority::Low}});
  const std::vector<std::tuple<std::size_t, ActorAction, long>> expected{
      {0, ActorAction::Predict, 0},
      {1, ActorAction::Predict, 0},
      {0, ActorAction::Correct, 0},
      {0, ActorAction::Predict, 1},
      {0, ActorAction::Correct, 1},
      {1, ActorAction::Correct, 0},
      {0, ActorAction::Predict, 2},
      {1, ActorAction::Predict, 1},
      {0, ActorAction::Correct, 2},
      {1, ActorAction::Correct, 1}};
  REQUIRE(plan.size() == expected.size());
  for (std::size_t i = 0; i < plan.size(); ++i) {
    CAPTURE(i);
    CHECK(plan[i].cluster == std::get<0>(expected[i]));
    CHECK(plan[i].action == std::get<1>(expected[i]));
    CHECK(plan[i].step == std::get<2>(expected[i]));
  }
}

} // namespace seissol::unit_test
