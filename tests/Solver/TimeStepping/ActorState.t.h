// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Solver/TimeStepping/ActorState.h"

#include <cmath>
#include <limits>
#include <variant>

namespace seissol::unit_test {
using namespace seissol::time_stepping;

// ---------------------------------------------------------------------------
// actorStateToString
// ---------------------------------------------------------------------------

TEST_CASE("actorStateToString" * doctest::test_suite("solver")) {
  CHECK(actorStateToString(ActorState::Corrected) == "Corrected");
  CHECK(actorStateToString(ActorState::Predicted) == "Predicted");
  CHECK(actorStateToString(ActorState::Synced) == "Synced");
}

// ---------------------------------------------------------------------------
// ClusterTimes
// ---------------------------------------------------------------------------

TEST_CASE("ClusterTimes::nextCorrectionTime" * doctest::test_suite("solver")) {
  ClusterTimes ct;
  ct.correctionTime = 5.0;
  ct.maxTimeStepSize = 1.0;

  SUBCASE("syncTime far away") {
    // correctionTime + maxTimeStepSize = 6.0 < 100.0
    CHECK(ct.nextCorrectionTime(100.0) == doctest::Approx(6.0));
  }

  SUBCASE("syncTime is the limit") {
    // correctionTime + maxTimeStepSize = 6.0 > 5.5
    CHECK(ct.nextCorrectionTime(5.5) == doctest::Approx(5.5));
  }

  SUBCASE("syncTime equals next step") {
    CHECK(ct.nextCorrectionTime(6.0) == doctest::Approx(6.0));
  }

  SUBCASE("Already at syncTime") { CHECK(ct.nextCorrectionTime(5.0) == doctest::Approx(5.0)); }
}

TEST_CASE("ClusterTimes::nextCorrectionSteps" * doctest::test_suite("solver")) {
  ClusterTimes ct;
  ct.stepsSinceLastSync = 3;
  ct.timeStepRate = 2;
  ct.stepsUntilSync = 10;

  SUBCASE("Normal case") {
    // stepsSinceLastSync + timeStepRate = 5 < 10
    CHECK(ct.nextCorrectionSteps() == 5);
  }

  SUBCASE("Clamped by stepsUntilSync") {
    ct.stepsUntilSync = 4;
    // stepsSinceLastSync + timeStepRate = 5 > 4
    CHECK(ct.nextCorrectionSteps() == 4);
  }

  SUBCASE("Exactly at sync") {
    ct.stepsUntilSync = 5;
    CHECK(ct.nextCorrectionSteps() == 5);
  }
}

TEST_CASE("ClusterTimes::timeStepSize" * doctest::test_suite("solver")) {
  ClusterTimes ct;
  ct.correctionTime = 5.0;
  ct.maxTimeStepSize = 1.0;

  SUBCASE("Regular step") {
    // syncTime - correctionTime = 5.0 > maxTimeStepSize = 1.0
    CHECK(ct.timeStepSize(10.0) == doctest::Approx(1.0));
  }

  SUBCASE("Reduced step near sync") {
    // syncTime - correctionTime = 0.3 < maxTimeStepSize = 1.0
    CHECK(ct.timeStepSize(5.3) == doctest::Approx(0.3));
  }

  SUBCASE("Zero step at sync") { CHECK(ct.timeStepSize(5.0) == doctest::Approx(0.0)); }
}

TEST_CASE("ClusterTimes::computeStepsUntilSyncTime" * doctest::test_suite("solver")) {
  ClusterTimes ct;
  ct.maxTimeStepSize = 0.5;
  ct.timeStepRate = 1;

  SUBCASE("Exact division") {
    // timeDiff = 2.0, timeStepRate * timeDiff / maxTimeStepSize = 1 * 2.0 / 0.5 = 4
    CHECK(ct.computeStepsUntilSyncTime(0.0, 2.0) == 4);
  }

  SUBCASE("Needs ceiling") {
    // timeDiff = 1.3, 1 * 1.3 / 0.5 = 2.6 → ceil → 3
    CHECK(ct.computeStepsUntilSyncTime(0.0, 1.3) == 3);
  }

  SUBCASE("Zero time difference") { CHECK(ct.computeStepsUntilSyncTime(5.0, 5.0) == 0); }

  SUBCASE("With timeStepRate > 1") {
    ct.timeStepRate = 4;
    // timeDiff = 1.0, 4 * 1.0 / 0.5 = 8
    CHECK(ct.computeStepsUntilSyncTime(0.0, 1.0) == 8);
  }
}

TEST_CASE("ClusterTimes get/setTimeStepSize" * doctest::test_suite("solver")) {
  ClusterTimes ct;
  ct.setTimeStepSize(0.42);
  CHECK(ct.getTimeStepSize() == doctest::Approx(0.42));
}

// ---------------------------------------------------------------------------
// NeighborCluster
// ---------------------------------------------------------------------------

TEST_CASE("NeighborCluster construction" * doctest::test_suite("solver")) {
  NeighborCluster nc(0.1, 2, Executor::Host);
  CHECK(nc.ct.maxTimeStepSize == doctest::Approx(0.1));
  CHECK(nc.ct.timeStepRate == 2);
  CHECK(nc.executor == Executor::Host);
  CHECK(nc.progress == nullptr);
  CHECK(nc.dataReadiness == DataReadiness::AfterPrediction);
}

} // namespace seissol::unit_test
