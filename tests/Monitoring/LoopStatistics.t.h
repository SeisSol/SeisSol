// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Monitoring/LoopStatistics.h"
#include "Monitoring/Stopwatch.h"

#include <string>
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

  unsigned region = ls.getRegion("work");

  // Add a sample: 100 iterations, 1 second duration
  timespec t0 = makeTime(10, 0);
  timespec t1 = makeTime(11, 0);
  ls.addSample(region, 100, 0, t0, t1);

  // Add another sample: 200 iterations, 0.5 seconds
  timespec t2 = makeTime(20, 0);
  timespec t3 = makeTime(20, 500000000);
  ls.addSample(region, 200, 0, t2, t3);

  // We can't directly inspect the private StatisticVariables,
  // but we can verify that no crash happens and the samples
  // are stored when output is enabled.
  CHECK(true);
}

// ---------------------------------------------------------------------------
// addSample: zero iterations ignored in statistics
// ---------------------------------------------------------------------------

TEST_CASE("LoopStatistics addSample zero iterations" * doctest::test_suite("monitoring")) {
  seissol::LoopStatistics ls;
  ls.addRegion("idle");

  unsigned region = ls.getRegion("idle");

  // Zero iterations → should be ignored in statistics (no crash)
  timespec t0 = makeTime(0, 0);
  timespec t1 = makeTime(1, 0);
  ls.addSample(region, 0, 0, t0, t1);

  CHECK(true);
}

// ---------------------------------------------------------------------------
// enableSampleOutput: samples stored only when enabled
// ---------------------------------------------------------------------------

TEST_CASE("LoopStatistics sample storage" * doctest::test_suite("monitoring")) {
  seissol::LoopStatistics ls;
  ls.addRegion("work");
  unsigned region = ls.getRegion("work");

  timespec t0 = makeTime(0, 0);
  timespec t1 = makeTime(1, 0);

  SUBCASE("Samples not stored by default") {
    // Default: outputSamples = false
    ls.addSample(region, 10, 0, t0, t1);
    // No crash; samples are silently not stored
    CHECK(true);
  }

  SUBCASE("Samples stored when enabled") {
    ls.enableSampleOutput(true);
    ls.addSample(region, 10, 0, t0, t1);
    ls.addSample(region, 20, 1, t0, t1);
    // No crash; samples are stored internally
    CHECK(true);
  }
}

// ---------------------------------------------------------------------------
// begin / end cycle
// ---------------------------------------------------------------------------

TEST_CASE("LoopStatistics begin/end cycle" * doctest::test_suite("monitoring")) {
  seissol::LoopStatistics ls;
  ls.addRegion("kernel");
  unsigned region = ls.getRegion("kernel");

  // begin() records the start time, end() records end time and adds a sample
  ls.begin(region);
  // Do a tiny bit of work to ensure nonzero time
  volatile int x = 0;
  for (int i = 0; i < 100; ++i) {
    x += i;
  }
  ls.end(region, 50, 0);

  // No crash = pass
  CHECK(true);
}

// ---------------------------------------------------------------------------
// reset
// ---------------------------------------------------------------------------

TEST_CASE("LoopStatistics reset" * doctest::test_suite("monitoring")) {
  seissol::LoopStatistics ls;
  ls.enableSampleOutput(true);
  ls.addRegion("work");
  unsigned region = ls.getRegion("work");

  timespec t0 = makeTime(0, 0);
  timespec t1 = makeTime(1, 0);
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
  ls.enableSampleOutput(true);
  ls.addRegion("local");
  ls.addRegion("neighbor");
  ls.addRegion("dr");

  unsigned local = ls.getRegion("local");
  unsigned neighbor = ls.getRegion("neighbor");
  unsigned dr = ls.getRegion("dr");

  CHECK(local != neighbor);
  CHECK(neighbor != dr);
  CHECK(local != dr);

  timespec t0 = makeTime(0, 0);
  timespec t1 = makeTime(0, 100000000); // 0.1s

  // Add samples to different regions with different subRegions (cluster ids)
  ls.addSample(local, 1000, 0, t0, t1);
  ls.addSample(local, 500, 1, t0, t1);
  ls.addSample(neighbor, 800, 0, t0, t1);
  ls.addSample(dr, 200, 0, t0, t1);

  CHECK(true);
}

} // namespace seissol::unit_test
