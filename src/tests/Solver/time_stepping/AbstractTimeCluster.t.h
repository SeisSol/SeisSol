// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Solver/time_stepping/AbstractTimeCluster.h"
#include "tests/TestHelper.h"
#include <iostream>
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
  MAKE_MOCK1(handleAdvancedPredictionTimeMessage, void(const NeighborCluster&), override);
  // NOLINTNEXTLINE
  MAKE_MOCK1(handleAdvancedCorrectionTimeMessage, void(const NeighborCluster&), override);
  // NOLINTNEXTLINE
  MAKE_MOCK1(printTimeoutMessage, void(std::chrono::seconds), override);
};

TEST_CASE("TimeCluster") {
  auto cluster = MockTimeCluster(1.0, 1);
  cluster.setSyncTime(10);
  cluster.reset();

  SUBCASE("Cluster start synced") {
    REQUIRE(cluster.synced());
    REQUIRE(cluster.getState() == ActorState::Synced);
  }

  SUBCASE("Cluster predicts after sync") {
    REQUIRE(cluster.getNextLegalAction() == ActorAction::RestartAfterSync);
    REQUIRE_CALL(cluster, start());
    auto result = cluster.act();
    REQUIRE(result.isStateChanged);
    REQUIRE(cluster.getState() == ActorState::Corrected);

    REQUIRE(cluster.getNextLegalAction() == ActorAction::Predict);
    REQUIRE_CALL(cluster, predict());
    result = cluster.act();
    REQUIRE(result.isStateChanged);
    REQUIRE(cluster.getState() == ActorState::Predicted);
  }
}

TEST_CASE("GTS Timesteping works") {
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
    cluster->reset();
  }

  // First, move from synced -> corrected
  for (auto& cluster : clusters) {
    REQUIRE_CALL(*cluster, start());
    cluster->act();
    REQUIRE(cluster->getState() == ActorState::Corrected);
  }

  bool isFinished = false;
  auto iteration = 0;
  while (!isFinished) {
    isFinished = true;

    ALLOW_CALL(cluster1, handleAdvancedCorrectionTimeMessage(ANY(NeighborCluster)));
    ALLOW_CALL(cluster2, handleAdvancedCorrectionTimeMessage(ANY(NeighborCluster)));
    ALLOW_CALL(cluster1, handleAdvancedPredictionTimeMessage(ANY(NeighborCluster)));
    ALLOW_CALL(cluster2, handleAdvancedPredictionTimeMessage(ANY(NeighborCluster)));

    for (auto& cluster : clusters) {
      REQUIRE(cluster->getState() == ActorState::Corrected);
      if (iteration < numberOfIterations) {
        REQUIRE(cluster->getNextLegalAction() == ActorAction::Predict);
        REQUIRE_CALL(*cluster, predict());
        cluster->act();
        REQUIRE(cluster->getState() == ActorState::Predicted);
      } else {
        REQUIRE(cluster->getNextLegalAction() == ActorAction::Sync);
        cluster->act();
        REQUIRE(cluster->getState() == ActorState::Synced);
      }
    }
    for (auto& cluster : clusters) {
      if (!cluster->synced()) {
        REQUIRE_CALL(*cluster, correct());
        REQUIRE(cluster->getNextLegalAction() == ActorAction::Correct);
        cluster->act();
        REQUIRE(cluster->getState() == ActorState::Corrected);
        isFinished = false;
      }
    }
    ++iteration;
  }

  for (auto& cluster : clusters) {
    REQUIRE(cluster->synced());
  }
}

TEST_CASE("LTS Timesteping works") {
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
    cluster->reset();
  }

  ALLOW_CALL(cluster1, handleAdvancedCorrectionTimeMessage(ANY(NeighborCluster)));
  ALLOW_CALL(cluster2, handleAdvancedCorrectionTimeMessage(ANY(NeighborCluster)));
  ALLOW_CALL(cluster1, handleAdvancedPredictionTimeMessage(ANY(NeighborCluster)));
  ALLOW_CALL(cluster2, handleAdvancedPredictionTimeMessage(ANY(NeighborCluster)));

  // First, move from synced -> corrected -> predicted
  for (auto& cluster : clusters) {
    REQUIRE_CALL(*cluster, start());
    cluster->act();
    REQUIRE(cluster->getState() == ActorState::Corrected);

    REQUIRE_CALL(*cluster, predict());
    cluster->act();
    REQUIRE(cluster->getState() == ActorState::Predicted);

    // Cluster should now be blocked by progress of other cluster
    REQUIRE(cluster->getNextLegalAction() == ActorAction::Nothing);
  }

  // Now, the second cluster should not be able to predict.
  REQUIRE(cluster2.getNextLegalAction() == ActorAction::Nothing);

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
  REQUIRE(cluster2.getState() == ActorState::Synced);
}

} // namespace seissol::unit_test
