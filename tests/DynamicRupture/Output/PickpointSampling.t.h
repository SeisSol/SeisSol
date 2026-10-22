// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "DynamicRupture/Output/PickpointSampling.h"

#include <cstddef>
#include <limits>
#include <vector>

namespace seissol::unit_test {
using namespace seissol::dr::output;

namespace {
/// The samples of consecutive steps of size `stepSize` from `start` up to `end`; the last step ends
/// at `end` exactly, as at a synchronization point.
std::vector<double>
    sampleSteps(const PickpointSampling& sampling, double start, double end, double stepSize) {
  std::vector<double> samples;
  double time = start;
  while (time < end) {
    const double size = (time + stepSize >= end * (1 - 1e-12)) ? end - time : stepSize;
    const auto times = sampling.times(time, size);
    CHECK(sampling.any(time, size) == !times.empty());
    samples.insert(samples.end(), times.begin(), times.end());
    time += size;
  }
  return samples;
}
} // namespace

TEST_CASE("PickpointSampling takes each sample in exactly one step" *
          doctest::test_suite("dynamicrupture")) {
  const double interval = 0.005;
  const double end = 0.1;
  const PickpointSampling sampling(interval, end);

  // step sizes that divide the interval, that do not, and that are larger than it
  for (const double stepSize : {0.001, 0.0007, 0.0014, 0.0123, 0.0333}) {
    CAPTURE(stepSize);
    const auto samples = sampleSteps(sampling, 0.0, end, stepSize);
    REQUIRE(samples.size() == 20);
    for (std::size_t i = 0; i < samples.size(); ++i) {
      CHECK(samples[i] == doctest::Approx((i + 1) * interval));
    }
  }
}

TEST_CASE("PickpointSampling is the same for all clusters" *
          doctest::test_suite("dynamicrupture")) {
  const PickpointSampling sampling(0.0025, 1.0);
  const double base = 0.00031;
  // a finer and a coarser cluster over the same interval between two synchronization points
  const auto fine = sampleSteps(sampling, 0.25, 0.5, base);
  const auto coarse = sampleSteps(sampling, 0.25, 0.5, 4 * base);
  REQUIRE(fine.size() == coarse.size());
  for (std::size_t i = 0; i < fine.size(); ++i) {
    CHECK(fine[i] == coarse[i]);
  }
  // 0.2525 to 0.5: the sample at the synchronization point belongs to the step that ends there
  CHECK(fine.size() == 100);
  CHECK(fine.back() == doctest::Approx(0.5));
}

TEST_CASE("PickpointSampling at the boundaries of steps" * doctest::test_suite("dynamicrupture")) {
  const PickpointSampling sampling(0.1, 10.0);

  SUBCASE("a sample at the end of a step belongs to that step") {
    // 0.1 + 0.2 is 0.30000000000000004, three times 0.1 is 0.30000000000000004 as well, and
    // 0.3 itself lies below both
    const auto first = sampling.times(0.1, 0.2);
    REQUIRE(first.size() == 2);
    CHECK(first[0] == doctest::Approx(0.2));
    CHECK(first[1] == doctest::Approx(0.3));
    CHECK(sampling.times(0.3, 0.05).empty());
    CHECK(sampling.times(0.1 + 0.2, 0.05).empty());
  }

  SUBCASE("no sample at the start of a step") {
    CHECK(sampling.times(0.0, 0.05).empty());
    CHECK(sampling.times(0.2, 0.05).empty());
  }

  SUBCASE("a step without a sample") {
    CHECK(sampling.times(0.21, 0.05).empty());
    CHECK_FALSE(sampling.any(0.21, 0.05));
  }
}

TEST_CASE("PickpointSampling at the end of the simulation" *
          doctest::test_suite("dynamicrupture")) {
  SUBCASE("an end between two samples is a sample of its own") {
    const PickpointSampling sampling(0.1, 0.95);
    const auto times = sampling.times(0.85, 0.1);
    REQUIRE(times.size() == 2);
    CHECK(times[0] == doctest::Approx(0.9));
    CHECK(times[1] == 0.95);
  }

  SUBCASE("an end on the grid is sampled once") {
    const PickpointSampling sampling(0.1, 1.0);
    const auto samples = sampleSteps(sampling, 0.0, 1.0, 0.03);
    REQUIRE(samples.size() == 10);
    CHECK(samples.back() == doctest::Approx(1.0));
  }

  SUBCASE("without a sampling interval, only the end") {
    for (const double interval : {0.0, std::numeric_limits<double>::infinity()}) {
      const PickpointSampling sampling(interval, 0.5);
      const auto samples = sampleSteps(sampling, 0.0, 0.5, 0.07);
      REQUIRE(samples.size() == 1);
      CHECK(samples[0] == 0.5);
    }
  }
}

} // namespace seissol::unit_test
