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
// MessageQueue
// ---------------------------------------------------------------------------

TEST_CASE("MessageQueue basic operations" * doctest::test_suite("solver")) {
  MessageQueue q;

  SUBCASE("Empty queue") {
    CHECK_FALSE(q.hasMessages());
    CHECK(q.size() == 0);
  }

  SUBCASE("Push and pop prediction message") {
    AdvancedPredictionTimeMessage msg{1.5, 3};
    q.push(msg);
    CHECK(q.hasMessages());
    CHECK(q.size() == 1);

    Message popped = q.pop();
    CHECK(std::holds_alternative<AdvancedPredictionTimeMessage>(popped));
    auto& result = std::get<AdvancedPredictionTimeMessage>(popped);
    CHECK(result.time == doctest::Approx(1.5));
    CHECK(result.stepsSinceSync == 3);
    CHECK_FALSE(q.hasMessages());
  }

  SUBCASE("Push and pop correction message") {
    AdvancedCorrectionTimeMessage msg{2.0, 5};
    q.push(msg);

    Message popped = q.pop();
    CHECK(std::holds_alternative<AdvancedCorrectionTimeMessage>(popped));
    auto& result = std::get<AdvancedCorrectionTimeMessage>(popped);
    CHECK(result.time == doctest::Approx(2.0));
    CHECK(result.stepsSinceSync == 5);
  }

  SUBCASE("FIFO ordering") {
    q.push(AdvancedPredictionTimeMessage{1.0, 1});
    q.push(AdvancedCorrectionTimeMessage{2.0, 2});
    q.push(AdvancedPredictionTimeMessage{3.0, 3});
    CHECK(q.size() == 3);

    auto m1 = q.pop();
    CHECK(std::holds_alternative<AdvancedPredictionTimeMessage>(m1));
    CHECK(std::get<AdvancedPredictionTimeMessage>(m1).time == doctest::Approx(1.0));

    auto m2 = q.pop();
    CHECK(std::holds_alternative<AdvancedCorrectionTimeMessage>(m2));
    CHECK(std::get<AdvancedCorrectionTimeMessage>(m2).time == doctest::Approx(2.0));

    auto m3 = q.pop();
    CHECK(std::holds_alternative<AdvancedPredictionTimeMessage>(m3));
    CHECK(std::get<AdvancedPredictionTimeMessage>(m3).time == doctest::Approx(3.0));

    CHECK(q.size() == 0);
    CHECK_FALSE(q.hasMessages());
  }
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
  CHECK(nc.inbox == nullptr);
  CHECK(nc.outbox == nullptr);
}

// ---------------------------------------------------------------------------
// DynamicRuptureScheduler
// ---------------------------------------------------------------------------

TEST_CASE("DynamicRuptureScheduler" * doctest::test_suite("solver")) {
  SUBCASE("Has DR faces") {
    const DynamicRuptureScheduler sched(100, 0.01);
    CHECK(sched.hasDynamicRuptureFaces());
    CHECK(sched.getOutputTimestep() == doctest::Approx(0.01));
  }

  SUBCASE("No DR faces") {
    const DynamicRuptureScheduler sched(0, 0.0);
    CHECK_FALSE(sched.hasDynamicRuptureFaces());
  }

  SUBCASE("mayComputeInterior starts at -1") {
    const DynamicRuptureScheduler sched(10, 0.01);
    // lastCorrectionStepsInterior is initialized to -1
    CHECK(sched.mayComputeInterior(0));
    CHECK(sched.mayComputeInterior(1));
  }

  SUBCASE("mayComputeInterior updates") {
    DynamicRuptureScheduler sched(10, 0.01);
    sched.setLastCorrectionStepsInterior(5);
    CHECK_FALSE(sched.mayComputeInterior(5));
    CHECK_FALSE(sched.mayComputeInterior(4));
    CHECK(sched.mayComputeInterior(6));
  }

  SUBCASE("setLastCorrectionStepsCopy does not affect interior") {
    DynamicRuptureScheduler sched(10, 0.01);
    sched.setLastCorrectionStepsCopy(5);
    // Interior is still at -1
    CHECK(sched.mayComputeInterior(0));
  }
}

} // namespace seissol::unit_test
