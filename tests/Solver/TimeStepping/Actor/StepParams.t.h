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
#include "Solver/TimeStepping/Actor/StepParams.h"

#include <chrono>
#include <cstddef>
#include <vector>

namespace seissol::unit_test {
using namespace seissol::solver;

namespace {

ClusterTimes clusterTimes(double maxTimeStepSize, long timeStepRate) {
  ClusterTimes times;
  times.maxTimeStepSize = maxTimeStepSize;
  times.timeStepRate = timeStepRate;
  return times;
}

/**
 * A cluster that records the step parameters of each prediction, together with what the
 * progress of its neighbors says about the same step at that moment.
 */
class RecordingCluster : public AbstractTimeCluster {
  public:
  struct Record {
    StepParams params;
    bool resetFromProgress{true};
    double subTimeStartFromProgress{0};
    bool hasLargerNeighbor{false};
  };

  RecordingCluster(double maxTimeStepSize, long timeStepRate)
      : AbstractTimeCluster(maxTimeStepSize, timeStepRate, Executor::Host) {}

  [[nodiscard]] const std::vector<Record>& records() const { return records_; }

  protected:
  void start() override {}

  void predict() override {
    Record record{};
    record.params = stepParams();
    for (const auto& neighbor : neighbors_) {
      if (neighbor.ct.timeStepRate > ct_.timeStepRate) {
        // the larger neighbor cannot correct its current step before this cluster has predicted
        // all of its steps inside it; hence its correction time is the start of that step
        record.hasLargerNeighbor = true;
        record.subTimeStartFromProgress = ct_.correctionTime - neighbor.ct.correctionTime;
        if (ct_.stepsSinceLastSync > neighbor.ct.stepsSinceLastSync) {
          record.resetFromProgress = false;
        }
      }
    }
    records_.push_back(record);
  }

  void correct() override {}
  void handleNeighborPrediction(const NeighborCluster& /*neighborCluster*/) override {}
  void handleNeighborCorrection(const NeighborCluster& /*neighborCluster*/) override {}
  void printTimeoutMessage(std::chrono::seconds /*timeSinceLastUpdate*/) override {}

  private:
  std::vector<Record> records_;
};

void advanceTo(const std::vector<RecordingCluster*>& clusters, double syncTime) {
  for (auto* cluster : clusters) {
    cluster->setSyncTime(syncTime);
    // dereference first due to a clang-tidy recommendation
    (*cluster).reset();
  }
  bool finished = false;
  std::size_t iterations = 0;
  while (!finished) {
    REQUIRE(iterations < 10000);
    ++iterations;
    finished = true;
    for (auto* cluster : clusters) {
      cluster->act();
      finished = finished && cluster->synced();
    }
  }
}

} // namespace

TEST_CASE("computeStepContext" * doctest::test_suite("solver")) {
  const auto times = clusterTimes(2.0, 2);

  SUBCASE("No neighbors") {
    const auto context = computeStepContext(times, {});
    CHECK(context.largerNeighborRate == 0);
    CHECK(context.largestTimeStepSize == 2.0);
  }

  SUBCASE("Smaller, equal and larger neighbors") {
    std::vector<NeighborCluster> neighbors;
    neighbors.emplace_back(1.0, 1, Executor::Host);
    neighbors.emplace_back(2.0, 2, Executor::Host);
    neighbors.emplace_back(4.0, 4, Executor::Host);
    const auto context = computeStepContext(times, neighbors);
    CHECK(context.largerNeighborRate == 4);
    CHECK(context.largestTimeStepSize == 4.0);
  }

  SUBCASE("Only smaller neighbors") {
    std::vector<NeighborCluster> neighbors;
    neighbors.emplace_back(1.0, 1, Executor::Host);
    const auto context = computeStepContext(times, neighbors);
    CHECK(context.largerNeighborRate == 0);
    CHECK(context.largestTimeStepSize == 2.0);
  }
}

TEST_CASE("computeStepParams" * doctest::test_suite("solver")) {
  SUBCASE("Without a larger neighbor, every step starts anew") {
    auto times = clusterTimes(0.5, 1);
    const StepContext context{0, 0.5};
    for (long step = 0; step < 4; ++step) {
      times.stepsSinceLastSync = step;
      times.correctionTime = 0.5 * static_cast<double>(step);
      const auto params = computeStepParams(times, 10.0, context);
      CHECK(params.resetBuffers);
      CHECK(params.subTimeStart == 0.0);
      CHECK(params.time == times.correctionTime);
      CHECK(params.timeStepSize == 0.5);
      CHECK(params.neighborTimeStepSize == 0.5);
    }
  }

  SUBCASE("Position inside the step of a neighbor with rate ratio 3") {
    auto times = clusterTimes(0.25, 2);
    const StepContext context{6, 0.75};
    for (long step = 0; step < 7; ++step) {
      times.stepsSinceLastSync = 2 * step;
      const auto params = computeStepParams(times, 10.0, context);
      const auto substep = step % 3;
      CHECK(params.resetBuffers == (substep == 0));
      CHECK(params.subTimeStart == 0.25 * static_cast<double>(substep));
      CHECK(params.neighborTimeStepSize == 0.75);
    }
  }

  SUBCASE("The step before the synchronization time is cut off") {
    auto times = clusterTimes(1.0, 1);
    times.correctionTime = 6.5;
    times.stepsSinceLastSync = 1;
    const auto params = computeStepParams(times, 7.0, StepContext{2, 2.0});
    CHECK(params.timeStepSize == doctest::Approx(0.5));
    CHECK(params.subTimeStart == 1.0);
    CHECK_FALSE(params.resetBuffers);
  }
}

TEST_CASE("Step parameters follow the progress of the neighbors" * doctest::test_suite("solver")) {
  // For each prediction, the step parameters must say the same as the progress of the larger
  // neighbor at that moment. Several synchronization intervals are covered, including ones that
  // do not end on a step boundary of the larger cluster.
  const auto check = [](double dt, long ratio, const std::vector<double>& syncTimes) {
    RecordingCluster small(dt, 1);
    RecordingCluster large(dt * static_cast<double>(ratio), ratio);
    small.connect(large);
    const std::vector<RecordingCluster*> clusters{&small, &large};
    for (const auto syncTime : syncTimes) {
      advanceTo(clusters, syncTime);
    }

    REQUIRE_FALSE(small.records().empty());
    for (const auto& record : small.records()) {
      REQUIRE(record.hasLargerNeighbor);
      CHECK(record.params.resetBuffers == record.resetFromProgress);
      CHECK(record.params.subTimeStart ==
            doctest::Approx(record.subTimeStartFromProgress).epsilon(1e-12));
    }
    for (const auto& record : large.records()) {
      CHECK_FALSE(record.hasLargerNeighbor);
      CHECK(record.params.resetBuffers);
      CHECK(record.params.subTimeStart == 0.0);
    }
  };

  SUBCASE("Ratio 2") { check(1.0, 2, {3.0, 7.0, 7.5}); }
  SUBCASE("Ratio 3") { check(1.0, 3, {4.0, 10.0}); }
  SUBCASE("Ratio 2, inexact time step") { check(0.1, 2, {0.3, 0.7, 1.35}); }
}

} // namespace seissol::unit_test
