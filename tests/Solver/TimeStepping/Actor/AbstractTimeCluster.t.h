// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Solver/TimeStepping/Actor/AbstractTimeCluster.h"
#include "TestHelper.h"

#include <iostream>
#include <string>
#include <utility>
#include <vector>
namespace seissol::unit_test {
using namespace time_stepping;

class MockTimeCluster : public time_stepping::AbstractTimeCluster {
  public:
  MockTimeCluster(double maxTimeStepSize, long timeStepRate)
      : AbstractTimeCluster(maxTimeStepSize, timeStepRate, Executor::Host) {}

  // NOLINTNEXTLINE
  MAKE_MOCK0(start, void(void), override);
  // NOLINTNEXTLINE
  MAKE_MOCK0(predict, void(void), override);
  // NOLINTNEXTLINE
  MAKE_MOCK0(correct, void(void), override);
  // NOLINTNEXTLINE
  MAKE_MOCK1(handleNeighborPrediction, void(const NeighborCluster&), override);
  // NOLINTNEXTLINE
  MAKE_MOCK1(handleNeighborCorrection, void(const NeighborCluster&), override);
  // NOLINTNEXTLINE
  MAKE_MOCK1(printTimeoutMessage, void(std::chrono::seconds), override);
};

TEST_CASE("TimeCluster" * doctest::test_suite("solver")) {
  auto cluster = MockTimeCluster(1.0, 1);
  cluster.setSyncTime(10);
  cluster.reset();

  SUBCASE("Cluster start synced") {
    CHECK(cluster.synced());
    CHECK(cluster.getState() == ActorState::Synced);
  }

  SUBCASE("Cluster predicts after sync") {
    REQUIRE(cluster.getNextLegalAction() == ActorAction::RestartAfterSync);
    REQUIRE_CALL(cluster, start());
    auto result = cluster.act();
    CHECK(result.isStateChanged);
    CHECK(cluster.getState() == ActorState::Corrected);

    REQUIRE(cluster.getNextLegalAction() == ActorAction::Predict);
    REQUIRE_CALL(cluster, predict());
    result = cluster.act();
    CHECK(result.isStateChanged);
    CHECK(cluster.getState() == ActorState::Predicted);
  }
}

TEST_CASE("GTS Timesteping works" * doctest::test_suite("solver")) {
  const double dt = 1.0;
  const auto numberOfIterations = 10;
  const double endTime = dt * numberOfIterations;
  auto cluster1 = MockTimeCluster(dt, 1);
  auto cluster2 = MockTimeCluster(dt, 1);
  auto clusters = std::vector<MockTimeCluster*>{
      &cluster1,
      &cluster2,
  };

  cluster1.connect(cluster2);

  for (auto* cluster : clusters) {
    cluster->setSyncTime(endTime);

    // derefernce first due to a clang-tidy recommendation
    (*cluster).reset();
  }

  // First, move from synced -> corrected
  for (auto& cluster : clusters) {
    REQUIRE_CALL(*cluster, start());
    cluster->act();
    CHECK(cluster->getState() == ActorState::Corrected);
  }

  bool isFinished = false;
  auto iteration = 0;
  while (!isFinished) {
    isFinished = true;

    ALLOW_CALL(cluster1, handleNeighborCorrection(ANY(NeighborCluster)));
    ALLOW_CALL(cluster2, handleNeighborCorrection(ANY(NeighborCluster)));
    ALLOW_CALL(cluster1, handleNeighborPrediction(ANY(NeighborCluster)));
    ALLOW_CALL(cluster2, handleNeighborPrediction(ANY(NeighborCluster)));

    for (auto& cluster : clusters) {
      REQUIRE(cluster->getState() == ActorState::Corrected);
      if (iteration < numberOfIterations) {
        REQUIRE(cluster->getNextLegalAction() == ActorAction::Predict);
        REQUIRE_CALL(*cluster, predict());
        cluster->act();
        CHECK(cluster->getState() == ActorState::Predicted);
      } else {
        CHECK(cluster->getNextLegalAction() == ActorAction::Sync);
        cluster->act();
        CHECK(cluster->getState() == ActorState::Synced);
      }
    }
    for (auto& cluster : clusters) {
      if (!cluster->synced()) {
        REQUIRE_CALL(*cluster, correct());
        CHECK(cluster->getNextLegalAction() == ActorAction::Correct);
        cluster->act();
        CHECK(cluster->getState() == ActorState::Corrected);
        isFinished = false;
      }
    }
    ++iteration;
  }

  for (auto& cluster : clusters) {
    CHECK(cluster->synced());
  }
}

TEST_CASE("LTS Timesteping works" * doctest::test_suite("solver")) {
  const double dt = 1.0;
  const auto numberOfIterations = 2;
  const double endTime = dt * numberOfIterations;
  auto cluster1 = MockTimeCluster(dt, 1);
  auto cluster2 = MockTimeCluster(2 * dt, 2);
  auto clusters = std::vector<MockTimeCluster*>{
      &cluster1,
      &cluster2,
  };

  cluster1.connect(cluster2);

  for (auto* cluster : clusters) {
    cluster->setSyncTime(endTime);

    // derefernce first due to a clang-tidy recommendation
    (*cluster).reset();
  }

  ALLOW_CALL(cluster1, handleNeighborCorrection(ANY(NeighborCluster)));
  ALLOW_CALL(cluster2, handleNeighborCorrection(ANY(NeighborCluster)));
  ALLOW_CALL(cluster1, handleNeighborPrediction(ANY(NeighborCluster)));
  ALLOW_CALL(cluster2, handleNeighborPrediction(ANY(NeighborCluster)));

  // First, move from synced -> corrected -> predicted
  for (auto& cluster : clusters) {
    REQUIRE_CALL(*cluster, start());
    cluster->act();
    REQUIRE(cluster->getState() == ActorState::Corrected);

    REQUIRE_CALL(*cluster, predict());
    cluster->act();
    REQUIRE(cluster->getState() == ActorState::Predicted);

    // Cluster should now be blocked by progress of other cluster
    CHECK(cluster->getNextLegalAction() == ActorAction::Nothing);
  }

  // Now, the second cluster should not be able to predict.
  CHECK(cluster2.getNextLegalAction() == ActorAction::Nothing);

  // The first should do: correction -> prediction -> correction
  REQUIRE_CALL(cluster1, correct());
  cluster1.act();
  REQUIRE(cluster1.getState() == ActorState::Corrected);

  REQUIRE_CALL(cluster1, predict());
  cluster1.act();
  REQUIRE(cluster1.getState() == ActorState::Predicted);

  REQUIRE_CALL(cluster1, correct());
  cluster1.act();
  REQUIRE(cluster1.getState() == ActorState::Corrected);

  cluster1.act();
  REQUIRE(cluster1.getState() == ActorState::Synced);

  REQUIRE_CALL(cluster2, correct());
  cluster2.act();
  REQUIRE(cluster2.getState() == ActorState::Corrected);

  cluster2.act();
  CHECK(cluster2.getState() == ActorState::Synced);
}

namespace {

struct ActionRecord {
  std::string cluster;
  char action;
  long step;
};

/**
 * A cluster that logs its predictions and corrections.
 */
class LoggingCluster : public time_stepping::AbstractTimeCluster {
  public:
  LoggingCluster(std::string name,
                 double maxTimeStepSize,
                 long timeStepRate,
                 std::vector<ActionRecord>& log,
                 DataReadiness readiness = DataReadiness::AfterPrediction)
      : AbstractTimeCluster(maxTimeStepSize, timeStepRate, Executor::Host), name_(std::move(name)),
        log_(log), readiness_(readiness) {}

  protected:
  [[nodiscard]] DataReadiness dataReadiness() const override { return readiness_; }
  void start() override {}
  void predict() override {
    log_.push_back({name_, 'P', ct_.predictionsSinceLastSync / ct_.timeStepRate});
  }
  void correct() override {
    log_.push_back({name_, 'C', ct_.stepsSinceLastSync / ct_.timeStepRate});
  }
  void handleNeighborPrediction(const NeighborCluster& /*neighbor*/) override {}
  void handleNeighborCorrection(const NeighborCluster& /*neighbor*/) override {}
  void printTimeoutMessage(std::chrono::seconds /*timeSinceLastUpdate*/) override {}

  private:
  std::string name_;
  std::vector<ActionRecord>& log_;
  DataReadiness readiness_;
};

void runUntilSync(const std::vector<LoggingCluster*>& clusters, double syncTime) {
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

std::size_t position(const std::vector<ActionRecord>& log,
                     const std::string& cluster,
                     char action,
                     long step) {
  for (std::size_t i = 0; i < log.size(); ++i) {
    if (log[i].cluster == cluster && log[i].action == action && log[i].step == step) {
      return i;
    }
  }
  FAIL("action not found: " << cluster << " " << action << " " << step);
  return log.size();
}

} // namespace

TEST_CASE("A face cluster runs between the predictions and corrections of its cells" *
          doctest::test_suite("solver")) {
  constexpr long Steps = 5;
  // different scheduling orders must lead to the same dependencies
  const std::vector<std::vector<std::size_t>> orders{{0, 1, 2}, {2, 0, 1}, {1, 2, 0}};
  for (const auto& order : orders) {
    std::vector<ActionRecord> log;
    LoggingCluster interior("interior", 1.0, 1, log);
    LoggingCluster copy("copy", 1.0, 1, log);
    LoggingCluster face("face", 1.0, 1, log, DataReadiness::AfterCorrection);
    interior.connect(copy);
    face.connect(interior);
    face.connect(copy);

    const std::vector<LoggingCluster*> all{&interior, &copy, &face};
    std::vector<LoggingCluster*> clusters;
    for (const auto index : order) {
      clusters.push_back(all[index]);
    }
    runUntilSync(clusters, static_cast<double>(Steps));

    for (long step = 0; step < Steps; ++step) {
      const auto faceStep = position(log, "face", 'C', step);
      for (const std::string cell : {"interior", "copy"}) {
        // the face needs the predictions of both sides
        CHECK(position(log, cell, 'P', step) < faceStep);
        // the cells need the result of the face
        CHECK(faceStep < position(log, cell, 'C', step));
        if (step + 1 < Steps) {
          // the next prediction overwrites what the face reads
          CHECK(faceStep < position(log, cell, 'P', step + 1));
          // the next face step overwrites what the cells read
          CHECK(position(log, cell, 'C', step) < position(log, "face", 'C', step + 1));
        }
      }
    }
  }
}

TEST_CASE("An observing cluster waits without being waited for" * doctest::test_suite("solver")) {
  constexpr long Steps = 4;
  std::vector<ActionRecord> log;
  LoggingCluster copy("copy", 1.0, 1, log);
  LoggingCluster ghost("ghost", 1.0, 1, log);
  LoggingCluster face("face", 1.0, 1, log, DataReadiness::AfterCorrection);
  copy.connect(ghost);
  face.connect(copy);
  face.observe(ghost);

  // the face waits for the copy and the ghost cluster, while the ghost cluster does not even know
  // the face cluster
  CHECK(face.getNeighborClusters()->size() == 2);
  CHECK(ghost.getNeighborClusters()->size() == 1);

  // the copy cluster sees the face as a cluster that provides its data with the correction
  const auto& copyNeighbors = *copy.getNeighborClusters();
  REQUIRE(copyNeighbors.size() == 2);
  CHECK(copyNeighbors[0].dataReadiness == DataReadiness::AfterPrediction);
  CHECK(copyNeighbors[1].dataReadiness == DataReadiness::AfterCorrection);

  runUntilSync({&ghost, &face, &copy}, static_cast<double>(Steps));

  for (long step = 0; step < Steps; ++step) {
    // the face needs the ghost data of the step
    CHECK(position(log, "ghost", 'P', step) < position(log, "face", 'C', step));
    CHECK(position(log, "face", 'C', step) < position(log, "copy", 'C', step));
  }
}

} // namespace seissol::unit_test
