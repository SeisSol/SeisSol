// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Numerical/Statistics.h"

#include <cmath>
#include <vector>

namespace seissol::unit_test {
using seissol::statistics::Summary;

// ---------------------------------------------------------------------------
// Summary from single scalar
// ---------------------------------------------------------------------------

TEST_CASE("Summary from scalar" * doctest::test_suite("monitoring")) {
  Summary s(42.0);
  CHECK(s.mean == doctest::Approx(42.0));
  CHECK(s.min == doctest::Approx(42.0));
  CHECK(s.max == doctest::Approx(42.0));
  CHECK(s.median == doctest::Approx(42.0));
  // std and sum are default-initialized (0) for the scalar constructor
}

// ---------------------------------------------------------------------------
// Summary from vector: basic statistics
// ---------------------------------------------------------------------------

TEST_CASE("Summary from vector" * doctest::test_suite("monitoring")) {
  const std::vector<double> values = {1.0, 2.0, 3.0, 4.0, 5.0};
  Summary s(values);

  CHECK(s.min == doctest::Approx(1.0));
  CHECK(s.max == doctest::Approx(5.0));
  CHECK(s.mean == doctest::Approx(3.0));
  CHECK(s.median == doctest::Approx(3.0)); // odd count → middle element
  CHECK(s.sum == doctest::Approx(15.0));

  // std = sqrt(mean of squared deviations)
  // variance = ((1-3)^2 + (2-3)^2 + (3-3)^2 + (4-3)^2 + (5-3)^2) / 5
  //          = (4 + 1 + 0 + 1 + 4) / 5 = 2.0
  CHECK(s.std == doctest::Approx(std::sqrt(2.0)));
}

// ---------------------------------------------------------------------------
// Summary from even-length vector: median is average of two middle
// ---------------------------------------------------------------------------

TEST_CASE("Summary even count median" * doctest::test_suite("monitoring")) {
  const std::vector<double> values = {10.0, 20.0, 30.0, 40.0};
  Summary s(values);

  CHECK(s.median == doctest::Approx(25.0)); // (20+30)/2
  CHECK(s.min == doctest::Approx(10.0));
  CHECK(s.max == doctest::Approx(40.0));
  CHECK(s.mean == doctest::Approx(25.0));
  CHECK(s.sum == doctest::Approx(100.0));
}

// ---------------------------------------------------------------------------
// Summary: all identical values
// ---------------------------------------------------------------------------

TEST_CASE("Summary all identical" * doctest::test_suite("monitoring")) {
  const std::vector<double> values = {7.0, 7.0, 7.0};
  Summary s(values);

  CHECK(s.min == doctest::Approx(7.0));
  CHECK(s.max == doctest::Approx(7.0));
  CHECK(s.mean == doctest::Approx(7.0));
  CHECK(s.median == doctest::Approx(7.0));
  CHECK(s.std == doctest::Approx(0.0));
  CHECK(s.sum == doctest::Approx(21.0));
}

// ---------------------------------------------------------------------------
// Summary: two elements
// ---------------------------------------------------------------------------

TEST_CASE("Summary two elements" * doctest::test_suite("monitoring")) {
  const std::vector<double> values = {3.0, 7.0};
  Summary s(values);

  CHECK(s.min == doctest::Approx(3.0));
  CHECK(s.max == doctest::Approx(7.0));
  CHECK(s.mean == doctest::Approx(5.0));
  CHECK(s.median == doctest::Approx(5.0)); // (3+7)/2
  CHECK(s.sum == doctest::Approx(10.0));
  // variance = ((3-5)^2 + (7-5)^2) / 2 = (4+4)/2 = 4 → std = 2
  CHECK(s.std == doctest::Approx(2.0));
}

// ---------------------------------------------------------------------------
// Summary: unsorted input gets sorted internally
// ---------------------------------------------------------------------------

TEST_CASE("Summary handles unsorted input" * doctest::test_suite("monitoring")) {
  const std::vector<double> values = {5.0, 1.0, 3.0, 4.0, 2.0};
  Summary s(values);

  CHECK(s.min == doctest::Approx(1.0));
  CHECK(s.max == doctest::Approx(5.0));
  CHECK(s.median == doctest::Approx(3.0));
}

// ---------------------------------------------------------------------------
// Summary: single element vector
// ---------------------------------------------------------------------------

TEST_CASE("Summary single element vector" * doctest::test_suite("monitoring")) {
  const std::vector<double> values = {99.0};
  Summary s(values);

  CHECK(s.min == doctest::Approx(99.0));
  CHECK(s.max == doctest::Approx(99.0));
  CHECK(s.mean == doctest::Approx(99.0));
  CHECK(s.median == doctest::Approx(99.0));
  CHECK(s.sum == doctest::Approx(99.0));
  CHECK(s.std == doctest::Approx(0.0));
}

} // namespace seissol::unit_test
