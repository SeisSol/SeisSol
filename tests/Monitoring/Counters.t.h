// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Monitoring/Counters.h"
#include "Monitoring/FlopCounter.h"
#include "Monitoring/Metric.h"
#include "Parallel/MPI.h"

#include <mpi.h>
#include <optional>
#include <vector>

namespace seissol::unit_test {

namespace counterstest {
inline int nodeCount() {
  int leader = seissol::Mpi::mpi.sharedMemMpiRank() == 0 ? 1 : 0;
  MPI_Allreduce(MPI_IN_PLACE, &leader, 1, MPI_INT, MPI_SUM, seissol::Mpi::mpi.comm());
  return leader;
}
} // namespace counterstest

using namespace counterstest;

TEST_CASE("Counters over a stretch of the run" * doctest::test_suite("monitoring")) {
  const auto ranks = seissol::Mpi::mpi.size();
  monitoring::Counters counters;

  double work = 1.0;
  const auto workCounter =
      counters.add(monitoring::CounterInfo{.name = "work", .rateName = "workPerSecond"},
                   [&] { return monitoring::CounterValue{.amount = work, .time = {}}; });

  // a counter of the nodes, with a clock of its own, which the first rank of each node reads
  double energy = 100.0;
  double clock = 10.0;
  const auto nodeCounter = counters.add(
      monitoring::CounterInfo{
          .name = "energy", .rateName = "power", .scope = monitoring::CounterScope::Node},
      [&]() -> std::optional<monitoring::CounterValue> {
        if (seissol::Mpi::mpi.sharedMemMpiRank() != 0) {
          return {};
        }
        return monitoring::CounterValue{.amount = energy, .time = clock};
      });

  CHECK(counters.size() == 2);
  CHECK(counters.find("energy") == nodeCounter);
  CHECK_FALSE(counters.find("none").has_value());
  CHECK(counters.info(workCounter).rateName == "workPerSecond");

  const auto from = counters.snapshot(1.0);
  REQUIRE(from.values.size() == 2);

  // registered after the first snapshot: counts from zero
  double late = 5.0;
  const auto lateCounter =
      counters.add(monitoring::CounterInfo{.name = "late", .rateName = "latePerSecond"},
                   [&] { return monitoring::CounterValue{.amount = late, .time = {}}; });

  work = 7.0;
  energy = 160.0;
  clock = 14.0;
  const auto to = counters.snapshot(3.0);
  const auto totals = counters.reduce(from, to);

  SUBCASE("A counter of the ranks sums over the ranks") {
    CHECK(totals.amounts[workCounter] == doctest::Approx(6.0 * ranks));
    // over the time of the ranks, two seconds
    CHECK(totals.rates[workCounter] == doctest::Approx(3.0 * ranks));
    CHECK(totals.contributors[workCounter] == ranks);
    CHECK(totals.counted(workCounter));
  }

  SUBCASE("A counter of the nodes sums over the nodes, on its own clock") {
    const auto nodes = nodeCount();
    CHECK(totals.amounts[nodeCounter] == doctest::Approx(60.0 * nodes));
    // over the clock of the counter, four seconds
    CHECK(totals.rates[nodeCounter] == doctest::Approx(15.0 * nodes));
    CHECK(totals.contributors[nodeCounter] == nodes);
  }

  SUBCASE("A counter registered later counts from zero") {
    CHECK(totals.amounts[lateCounter] == doctest::Approx(5.0 * ranks));
  }

  SUBCASE("Several stretches reduce at once") {
    const auto again = counters.snapshot(5.0);
    const auto both = counters.reduce({{&from, &to}, {&to, &again}});
    REQUIRE(both.size() == 2);
    CHECK(both[0].amounts[workCounter] == doctest::Approx(totals.amounts[workCounter]));
    CHECK(both[1].amounts[workCounter] == doctest::Approx(0.0));
    CHECK(both[1].contributors[workCounter] == ranks);
  }

  SUBCASE("A counter that no rank counts is not counted") {
    monitoring::Counters none;
    none.add(monitoring::CounterInfo{.name = "nothing", .rateName = "nothingPerSecond"},
             []() -> std::optional<monitoring::CounterValue> { return {}; });
    const auto result = none.reduce(none.snapshot(0.0), none.snapshot(1.0));
    CHECK_FALSE(result.counted(0));
    CHECK(result.amounts[0] == 0.0);
  }
}

TEST_CASE("The FLOP counter registers its totals as counters" * doctest::test_suite("monitoring")) {
  monitoring::FlopCounter flops;
  monitoring::Counters counters;
  const auto one = flops.addMetric("one", "WP");
  const auto two = flops.addMetric("two", "DR");
  flops.registerCounters(counters);

  flops.incrementMetric(
      one, PerformanceEstimate{.hardwareFlop = 10, .nonzeroFlop = 4, .bytes = 1, .kernelBytes = 8});
  flops.incrementMetric(
      two, PerformanceEstimate{.hardwareFlop = 5, .nonzeroFlop = 2, .bytes = 1, .kernelBytes = 3});

  const auto snapshot = counters.snapshot(0.0);
  const auto value = [&](const char* name) {
    return snapshot.values.at(counters.find(name).value()).value().amount;
  };
  CHECK(value("hardwareFlop") == 15.0);
  CHECK(value("nonzeroFlop") == 6.0);
  CHECK(value("kernelBytes") == 11.0);
  CHECK(counters.info(counters.find("hardwareFlop").value()).rateName == "hardwareFlopPerSecond");
}

} // namespace seissol::unit_test
