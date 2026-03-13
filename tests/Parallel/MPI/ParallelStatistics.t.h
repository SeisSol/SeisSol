// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Numerical/Statistics.h"
#include "Parallel/MPI.h"

#include <cmath>

namespace seissol::unit_test {
using seissol::statistics::Summary;

// ---------------------------------------------------------------------------
// parallelSummary: each rank contributes its rank as a value
// ---------------------------------------------------------------------------

TEST_CASE("parallelSummary with distinct rank values" * doctest::test_suite("mpi")) {
  const auto& mpi = seissol::Mpi::mpi;
  const int rank = mpi.rank();
  const int size = mpi.size();

  // Each rank contributes its rank index as a double
  double value = static_cast<double>(rank);
  auto summary = seissol::statistics::parallelSummary(value);

  if (rank == 0) {
    CHECK(summary.min == doctest::Approx(0.0));
    CHECK(summary.max == doctest::Approx(size - 1));

    // mean = (0 + 1 + ... + (size-1)) / size = (size-1)/2
    double expectedMean = (size - 1) / 2.0;
    CHECK(summary.mean == doctest::Approx(expectedMean));

    // sum = 0 + 1 + ... + (size-1) = size*(size-1)/2
    double expectedSum = size * (size - 1) / 2.0;
    CHECK(summary.sum == doctest::Approx(expectedSum));

    // median for size ranks: depends on parity
    if (size % 2 == 1) {
      // Odd: median is the middle value
      CHECK(summary.median == doctest::Approx((size - 1) / 2.0));
    } else {
      // Even: average of two middle values
      double mid = (size / 2 - 1 + size / 2) / 2.0;
      CHECK(summary.median == doctest::Approx(mid));
    }
  }
}

// ---------------------------------------------------------------------------
// parallelSummary: all ranks send the same value
// ---------------------------------------------------------------------------

TEST_CASE("parallelSummary identical values" * doctest::test_suite("mpi")) {
  const auto& mpi = seissol::Mpi::mpi;
  const int rank = mpi.rank();
  const int size = mpi.size();

  auto summary = seissol::statistics::parallelSummary(42.0);

  if (rank == 0) {
    CHECK(summary.min == doctest::Approx(42.0));
    CHECK(summary.max == doctest::Approx(42.0));
    CHECK(summary.mean == doctest::Approx(42.0));
    CHECK(summary.median == doctest::Approx(42.0));
    CHECK(summary.std == doctest::Approx(0.0));
    CHECK(summary.sum == doctest::Approx(42.0 * size));
  }
}

// ---------------------------------------------------------------------------
// parallelSummary: two distinct groups (rank-dependent)
// ---------------------------------------------------------------------------

TEST_CASE("parallelSummary binary values" * doctest::test_suite("mpi")) {
  const auto& mpi = seissol::Mpi::mpi;
  const int rank = mpi.rank();
  const int size = mpi.size();

  // Even ranks send 0.0, odd ranks send 100.0
  double value = (rank % 2 == 0) ? 0.0 : 100.0;
  auto summary = seissol::statistics::parallelSummary(value);

  if (rank == 0) {
    CHECK(summary.min == doctest::Approx(0.0));

    if (size > 1) {
      CHECK(summary.max == doctest::Approx(100.0));
    } else {
      // Single rank: only even → value is 0.0
      CHECK(summary.max == doctest::Approx(0.0));
    }
  }
}

} // namespace seissol::unit_test
