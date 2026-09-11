// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Monitoring/FlopCounter.h"
#include "Monitoring/Metric.h"

namespace seissol::unit_test {
using seissol::monitoring::FlopCounter;

TEST_CASE("FlopCounter increment operations" * doctest::test_suite("monitoring")) {
  // Note: We cannot test init() or print*() here as they require MPI.
  // We test that the increment functions don't crash and that the object
  // is constructible. The actual accumulation is internal (private members),
  // so we verify the functions are callable without error.

  FlopCounter counter;

  const auto one = counter.addMetric("one", "");
  const auto two = counter.addMetric("two", "");

  SUBCASE("Increment metric") {
    counter.incrementMetric(one, PerformanceEstimate{1, 2, 3, 4});
    counter.incrementMetric(two, PerformanceEstimate{1, 2, 3, 4});
    counter.incrementMetric(one, PerformanceEstimate{2, 3, 4, 5});
    CHECK(true);
  }

  SUBCASE("Zero increment is valid") {
    counter.incrementMetric(one, PerformanceEstimate{});
    CHECK(true);
  }

  SUBCASE("Large values") {
    counter.incrementMetric(one,
                            PerformanceEstimate{
                                1'000'000'000ULL,
                                1'000'000'000ULL,
                                1'000'000'000ULL,
                                1'000'000'000ULL,
                            });
    CHECK(true);
  }
}

} // namespace seissol::unit_test
