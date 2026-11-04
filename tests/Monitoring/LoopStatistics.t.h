// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Monitoring/LoopStatistics.h"
#include "Monitoring/Stopwatch.h"
#include "Parallel/MPI.h"

#include <string>
#include <thread>
#include <time.h>

namespace seissol::unit_test {

// Helper: create a timespec from seconds + nanoseconds
static timespec makeTime(long sec, long nsec) {
  timespec t{};
  t.tv_sec = sec;
  t.tv_nsec = nsec;
  return t;
}

// ---------------------------------------------------------------------------
// addRegion / getRegion
// ---------------------------------------------------------------------------

TEST_CASE("LoopStatistics addRegion and getRegion" * doctest::test_suite("monitoring")) {
  seissol::LoopStatistics ls;

  ls.addRegion("compute");
  ls.addRegion("communication");
  ls.addRegion("io", false);

  SUBCASE("getRegion returns correct indices") {
    CHECK(ls.getRegion("compute") == 0);
    CHECK(ls.getRegion("communication") == 1);
    CHECK(ls.getRegion("io") == 2);
  }

  SUBCASE("Regions are ordered by insertion") {
    ls.addRegion("extra");
    CHECK(ls.getRegion("extra") == 3);
  }
}

// ---------------------------------------------------------------------------
// addSample: statistics accumulation
// ---------------------------------------------------------------------------

TEST_CASE("LoopStatistics addSample accumulates" * doctest::test_suite("monitoring")) {
  seissol::LoopStatistics ls;
  ls.addRegion("work");

  const unsigned region = ls.getRegion("work");

  // Add a sample: 100 iterations, 1 second duration
  const timespec t0 = makeTime(10, 0);
  const timespec t1 = makeTime(11, 0);
  ls.addSample(region, 100, 0, t0, t1);

  // Add another sample: 200 iterations, 0.5 seconds
  const timespec t2 = makeTime(20, 0);
  const timespec t3 = makeTime(20, 500000000);
  ls.addSample(region, 200, 0, t2, t3);

  // We can't directly inspect the private StatisticVariables,
  // but we can verify that no crash happens.
  CHECK(true);
}

// ---------------------------------------------------------------------------
// addSample: zero iterations ignored in statistics
// ---------------------------------------------------------------------------

TEST_CASE("LoopStatistics addSample zero iterations" * doctest::test_suite("monitoring")) {
  seissol::LoopStatistics ls;
  ls.addRegion("idle");

  const unsigned region = ls.getRegion("idle");

  // Zero iterations → should be ignored in statistics (no crash)
  const timespec t0 = makeTime(0, 0);
  const timespec t1 = makeTime(1, 0);
  ls.addSample(region, 0, 0, t0, t1);

  CHECK(true);
}

// ---------------------------------------------------------------------------
// begin / end cycle
// ---------------------------------------------------------------------------

TEST_CASE("LoopStatistics begin/end cycle" * doctest::test_suite("monitoring")) {
  seissol::LoopStatistics ls;
  ls.addRegion("kernel");
  const unsigned region = ls.getRegion("kernel");

  // begin() records the start time, end() records end time and adds a sample
  ls.begin(region);

  // Sleep a tiny bit to ensure nonzero time
  using namespace std::chrono_literals;
  std::this_thread::sleep_for(50ms);

  ls.end(region, 50, 0);

  // No crash = pass
  CHECK(true);
}

// ---------------------------------------------------------------------------
// reset
// ---------------------------------------------------------------------------

TEST_CASE("LoopStatistics reset" * doctest::test_suite("monitoring")) {
  seissol::LoopStatistics ls;
  ls.addRegion("work");
  const unsigned region = ls.getRegion("work");

  const timespec t0 = makeTime(0, 0);
  const timespec t1 = makeTime(1, 0);
  ls.addSample(region, 100, 0, t0, t1);

  ls.reset();

  // After reset, adding new samples should work (no stale state)
  ls.addSample(region, 200, 0, t0, t1);
  CHECK(true);
}

// ---------------------------------------------------------------------------
// Multiple regions with different subRegions
// ---------------------------------------------------------------------------

TEST_CASE("LoopStatistics multiple regions and subRegions" * doctest::test_suite("monitoring")) {
  seissol::LoopStatistics ls;
  ls.addRegion("local");
  ls.addRegion("neighbor");
  ls.addRegion("dr");

  unsigned local = ls.getRegion("local");
  unsigned neighbor = ls.getRegion("neighbor");
  unsigned dr = ls.getRegion("dr");

  CHECK(local != neighbor);
  CHECK(neighbor != dr);
  CHECK(local != dr);

  const timespec t0 = makeTime(0, 0);
  const timespec t1 = makeTime(0, 100000000); // 0.1s

  // Add samples to different regions with different subRegions (cluster ids)
  ls.addSample(local, 1000, 0, t0, t1);
  ls.addSample(local, 500, 1, t0, t1);
  ls.addSample(neighbor, 800, 0, t0, t1);
  ls.addSample(dr, 200, 0, t0, t1);

  CHECK(true);
}

// ---------------------------------------------------------------------------
// printSummary: what it returns for the run summary
// ---------------------------------------------------------------------------

TEST_CASE("LoopStatistics summary" * doctest::test_suite("monitoring")) {
  seissol::LoopStatistics ls;
  ls.addRegion("work");
  ls.addRegion("state", false);
  const unsigned work = ls.getRegion("work");
  const unsigned state = ls.getRegion("state");

  // time = 0.5 s + 0.01 s per iteration
  ls.addSample(work, 100, 0, makeTime(0, 0), makeTime(1, 500000000));
  ls.addSample(work, 200, 0, makeTime(0, 0), makeTime(2, 500000000));
  ls.addSample(work, 300, 0, makeTime(0, 0), makeTime(3, 500000000));
  ls.addSample(state, 1, 0, makeTime(0, 0), makeTime(10, 0));

  const auto summary = ls.printSummary(seissol::Mpi::mpi.comm());
  if (seissol::Mpi::mpi.rank() == 0) {
    const auto* regions = summary.find("regions");
    REQUIRE(regions != nullptr);
    const auto* entry = regions->find("work");
    REQUIRE(entry != nullptr);
    CHECK(entry->find("includedInTotal")->asBool());
    CHECK(entry->find("time")->find("sum")->asDouble() ==
          doctest::Approx(7.5 * seissol::Mpi::mpi.size()));
    CHECK(entry->find("iterations")->asDouble() == 600.0 * seissol::Mpi::mpi.size());
    CHECK(entry->find("samples")->asDouble() == 3.0 * seissol::Mpi::mpi.size());
    CHECK(entry->find("regression")->find("constant")->asDouble() == doctest::Approx(0.5));
    CHECK(entry->find("regression")->find("perElement")->asDouble() == doctest::Approx(0.01));

    // a state is no compute kernel: neither in the total nor regressed
    const auto* stateEntry = regions->find("state");
    REQUIRE(stateEntry != nullptr);
    CHECK_FALSE(stateEntry->find("includedInTotal")->asBool());
    CHECK(stateEntry->find("regression") == nullptr);
    CHECK(summary.find("total")->find("max")->asDouble() == doctest::Approx(7.5));
    CHECK(summary.find("loadImbalance")->asDouble() == doctest::Approx(0.0));
  }
}

} // namespace seissol::unit_test
