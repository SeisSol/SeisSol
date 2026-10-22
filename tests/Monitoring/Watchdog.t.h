// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Monitoring/Watchdog.h"

#include <chrono>

namespace seissol::unit_test {

namespace watchdogtest {

inline monitoring::Watchdog::Settings watchdogSettings() {
  monitoring::Watchdog::Settings settings;
  settings.timeout = std::chrono::seconds(10);
  settings.period = std::chrono::seconds(1);
  return settings;
}

} // namespace watchdogtest

using namespace watchdogtest;

// the thread of the watchdog does not run here: the tests take its looks themselves, at times of
// their choice, and the watchdog never aborts by itself

TEST_CASE("The watchdog aborts a stall in the time stepping a period after reporting it" *
          doctest::test_suite("monitoring")) {
  using Finding = monitoring::Watchdog::Finding;
  using Phase = monitoring::Watchdog::Phase;
  monitoring::Watchdog watchdog(watchdogSettings());
  double progress = 0;
  watchdog.addProbe("progress", [&progress] { return progress; });
  const auto start = monitoring::Watchdog::Clock::now();
  watchdog.setPhase(Phase::Stepping, "the time stepping");

  CHECK(watchdog.watch(start + std::chrono::seconds(5)) == Finding::Progress);
  // progress puts the stall off
  progress = 1;
  CHECK(watchdog.watch(start + std::chrono::seconds(8)) == Finding::Progress);
  CHECK(watchdog.watch(start + std::chrono::seconds(17)) == Finding::Progress);
  CHECK(watchdog.watch(start + std::chrono::seconds(19)) == Finding::Stall);
  CHECK(watchdog.watch(start + std::chrono::seconds(20)) == Finding::Abort);
}

TEST_CASE("The watchdog only reports a phase outside of the time stepping that takes long" *
          doctest::test_suite("monitoring")) {
  using Finding = monitoring::Watchdog::Finding;
  using Phase = monitoring::Watchdog::Phase;
  monitoring::Watchdog watchdog(watchdogSettings());
  const auto start = monitoring::Watchdog::Clock::now();

  // e.g. a checkpoint that takes longer than the timeout
  watchdog.setPhase(Phase::Other, "a synchronization point");
  CHECK(watchdog.watch(start + std::chrono::seconds(5)) == Finding::Progress);
  CHECK(watchdog.watch(start + std::chrono::seconds(11)) == Finding::Stall);
  CHECK(watchdog.watch(start + std::chrono::seconds(12)) == Finding::Stall);
  CHECK(watchdog.watch(start + std::chrono::seconds(3600)) == Finding::Stall);

  // the time stepping after it gets its own timeout, and its stalls abort again
  const auto stepping = monitoring::Watchdog::Clock::now();
  watchdog.setPhase(Phase::Stepping, "the time stepping");
  CHECK(watchdog.watch(stepping + std::chrono::seconds(5)) == Finding::Progress);
  CHECK(watchdog.watch(stepping + std::chrono::seconds(11)) == Finding::Stall);
  CHECK(watchdog.watch(stepping + std::chrono::seconds(12)) == Finding::Abort);
}

TEST_CASE("The watchdog does not watch an idle simulation" * doctest::test_suite("monitoring")) {
  using Finding = monitoring::Watchdog::Finding;
  using Phase = monitoring::Watchdog::Phase;
  monitoring::Watchdog watchdog(watchdogSettings());
  const auto start = monitoring::Watchdog::Clock::now();
  CHECK(watchdog.watch(start + std::chrono::seconds(3600)) == Finding::Progress);
  watchdog.setPhase(Phase::Idle, "the end");
  CHECK(watchdog.watch(start + std::chrono::seconds(3600)) == Finding::Progress);
}

} // namespace seissol::unit_test
