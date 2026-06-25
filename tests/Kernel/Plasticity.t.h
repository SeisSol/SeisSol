// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Kernels/Plasticity.h"

#include <cmath>
#include <cstdint>

namespace seissol::unit_test {
using seissol::kernels::Plasticity;

// ---------------------------------------------------------------------------
// computeRelaxTime: pure constexpr math
// ---------------------------------------------------------------------------

TEST_CASE("Plasticity computeRelaxTime" * doctest::test_suite("kernel")) {
  SUBCASE("tV = 0 → factor = 1 (instantaneous relaxation)") {
    CHECK(Plasticity::computeRelaxTime(0.0, 1.0) == doctest::Approx(1.0));
    CHECK(Plasticity::computeRelaxTime(0.0, 0.5) == doctest::Approx(1.0));
    CHECK(Plasticity::computeRelaxTime(0.0, 0.001) == doctest::Approx(1.0));
  }

  SUBCASE("tV > 0 → exponential decay factor") {
    // computeRelaxTime(tV, dt) = -expm1(-dt/tV) = 1 - exp(-dt/tV)
    const double tV = 0.1;
    const double dt = 0.05;
    const double expected = 1.0 - std::exp(-dt / tV);
    CHECK(Plasticity::computeRelaxTime(tV, dt) == doctest::Approx(expected));
  }

  SUBCASE("Small dt/tV ratio (nearly no relaxation)") {
    const double tV = 100.0;
    const double dt = 0.001;
    double result = Plasticity::computeRelaxTime(tV, dt);
    // For small dt/tV: result ≈ dt/tV
    CHECK(result == doctest::Approx(dt / tV).epsilon(1e-4));
    CHECK(result > 0.0);
    CHECK(result < 1.0);
  }

  SUBCASE("Large dt/tV ratio (full relaxation)") {
    const double tV = 0.001;
    const double dt = 10.0;
    double result = Plasticity::computeRelaxTime(tV, dt);
    // For large dt/tV: result ≈ 1
    CHECK(result == doctest::Approx(1.0).epsilon(1e-6));
  }

  SUBCASE("dt = tV → known value") {
    const double tV = 1.0;
    const double dt = 1.0;
    // 1 - exp(-1) ≈ 0.6321
    CHECK(Plasticity::computeRelaxTime(tV, dt) == doctest::Approx(1.0 - std::exp(-1.0)));
  }

  SUBCASE("Result is always in (0, 1] for tV > 0") {
    for (const double tV : {0.01, 0.1, 1.0, 10.0, 100.0}) {
      for (const double dt : {0.001, 0.01, 0.1, 1.0, 10.0}) {
        double result = Plasticity::computeRelaxTime(tV, dt);
        CHECK(result > 0.0);
        CHECK(result <= 1.0);
      }
    }
  }

  SUBCASE("Monotone in dt for fixed tV") {
    const double tV = 1.0;
    double prev = 0.0;
    for (const double dt : {0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0}) {
      double result = Plasticity::computeRelaxTime(tV, dt);
      CHECK(result > prev);
      prev = result;
    }
  }
}

// ---------------------------------------------------------------------------
// flopsPlasticity: generated code constants
// ---------------------------------------------------------------------------

TEST_CASE("Plasticity flopsPlasticity" * doctest::test_suite("kernel")) {
  std::uint64_t nzCheck = 0;
  std::uint64_t hwCheck = 0;
  std::uint64_t nzYield = 0;
  std::uint64_t hwYield = 0;

  Plasticity::flopsPlasticity(nzCheck, hwCheck, nzYield, hwYield);

  SUBCASE("Check flops are positive") {
    CHECK(nzCheck > 0);
    CHECK(hwCheck > 0);
  }

  SUBCASE("Yield flops are positive") {
    CHECK(nzYield > 0);
    CHECK(hwYield > 0);
  }

  SUBCASE("Hardware flops >= nonzero flops") {
    CHECK(hwCheck >= nzCheck);
    CHECK(hwYield >= nzYield);
  }
}

} // namespace seissol::unit_test
