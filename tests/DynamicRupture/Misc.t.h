// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "DynamicRupture/Misc.h"
#include "TestHelper.h"

#include <cmath>
#include <limits>

namespace seissol::unit_test {

using namespace seissol::dr::misc;

// ---------------------------------------------------------------------------
// power<Exp>
// ---------------------------------------------------------------------------

TEST_CASE("DR Misc power" * doctest::test_suite("dynamicrupture")) {

  SUBCASE("power<0> is always 1") {
    CHECK(power<0>(0.0) == 1.0);
    CHECK(power<0>(5.0) == 1.0);
    CHECK(power<0>(-3.0) == 1.0);
  }

  SUBCASE("power<1> is identity") {
    CHECK(power<1>(7.0) == 7.0);
    CHECK(power<1>(-2.5) == -2.5);
  }

  SUBCASE("power<2> is square") {
    CHECK(power<2>(3.0) == doctest::Approx(9.0));
    CHECK(power<2>(-4.0) == doctest::Approx(16.0));
  }

  SUBCASE("power<5>") { CHECK(power<5>(2.0) == doctest::Approx(32.0)); }

  SUBCASE("power<8>") { CHECK(power<8>(2.0) == doctest::Approx(256.0)); }

  SUBCASE("power with float type") { CHECK(power<3>(2.0f) == doctest::Approx(8.0f)); }
}

// ---------------------------------------------------------------------------
// square
// ---------------------------------------------------------------------------

TEST_CASE("DR Misc square" * doctest::test_suite("dynamicrupture")) {

  SUBCASE("Single argument") {
    CHECK(square(3.0) == doctest::Approx(9.0));
    CHECK(square(-5.0) == doctest::Approx(25.0));
    CHECK(square(0.0) == doctest::Approx(0.0));
  }

  SUBCASE("Two arguments — sum of squares") {
    // 3^2 + 4^2 = 25
    CHECK(square(3.0, 4.0) == doctest::Approx(25.0));
  }

  SUBCASE("Three arguments") {
    // 1^2 + 2^2 + 3^2 = 14
    CHECK(square(1.0, 2.0, 3.0) == doctest::Approx(14.0));
  }
}

// ---------------------------------------------------------------------------
// magnitude
// ---------------------------------------------------------------------------

TEST_CASE("DR Misc magnitude" * doctest::test_suite("dynamicrupture")) {

  SUBCASE("2D Pythagorean triple") { CHECK(magnitude(3.0, 4.0) == doctest::Approx(5.0)); }

  SUBCASE("3D known result") {
    // sqrt(1 + 4 + 4) = 3
    CHECK(magnitude(1.0, 2.0, 2.0) == doctest::Approx(3.0));
  }

  SUBCASE("Zero vector") { CHECK(magnitude(0.0, 0.0) == doctest::Approx(0.0)); }

  SUBCASE("Unit vector") { CHECK(magnitude(1.0, 0.0, 0.0) == doctest::Approx(1.0)); }
}

// ---------------------------------------------------------------------------
// clamp
// ---------------------------------------------------------------------------

TEST_CASE("DR Misc clamp" * doctest::test_suite("dynamicrupture")) {

  SUBCASE("Value within range") { CHECK(clamp(5.0, 0.0, 10.0) == doctest::Approx(5.0)); }

  SUBCASE("Value below minimum") { CHECK(clamp(-3.0, 0.0, 10.0) == doctest::Approx(0.0)); }

  SUBCASE("Value above maximum") { CHECK(clamp(15.0, 0.0, 10.0) == doctest::Approx(10.0)); }

  SUBCASE("Value at boundary") {
    CHECK(clamp(0.0, 0.0, 10.0) == doctest::Approx(0.0));
    CHECK(clamp(10.0, 0.0, 10.0) == doctest::Approx(10.0));
  }

  SUBCASE("Integer variant") {
    CHECK(clamp(5, 0, 10) == 5);
    CHECK(clamp(-1, 0, 10) == 0);
    CHECK(clamp(20, 0, 10) == 10);
  }
}

// ---------------------------------------------------------------------------
// computeStrikeAndDipVectors
// ---------------------------------------------------------------------------

TEST_CASE("DR Misc computeStrikeAndDipVectors" * doctest::test_suite("dynamicrupture")) {
  constexpr double Eps = 1e-12;

  auto dotProduct = [](const double* a, const double* b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  };

  auto vecLength = [](const double* v) {
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  };

  SUBCASE("Normal along y-axis [0,1,0]") {
    double normal[3] = {0.0, 1.0, 0.0};
    double strike[3] = {};
    double dip[3] = {};

    computeStrikeAndDipVectors(normal, strike, dip);

    // strike and dip should be unit vectors
    CHECK(vecLength(strike) == doctest::Approx(1.0).epsilon(Eps));
    CHECK(vecLength(dip) == doctest::Approx(1.0).epsilon(Eps));

    // strike and dip should be orthogonal to normal
    CHECK(dotProduct(strike, normal) == doctest::Approx(0.0).epsilon(Eps));
    CHECK(dotProduct(dip, normal) == doctest::Approx(0.0).epsilon(Eps));

    // strike and dip should be orthogonal to each other
    CHECK(dotProduct(strike, dip) == doctest::Approx(0.0).epsilon(Eps));
  }

  SUBCASE("Normal along x-axis [1,0,0]") {
    double normal[3] = {1.0, 0.0, 0.0};
    double strike[3] = {};
    double dip[3] = {};

    computeStrikeAndDipVectors(normal, strike, dip);

    CHECK(vecLength(strike) == doctest::Approx(1.0).epsilon(Eps));
    CHECK(vecLength(dip) == doctest::Approx(1.0).epsilon(Eps));
    CHECK(dotProduct(strike, normal) == doctest::Approx(0.0).epsilon(Eps));
    CHECK(dotProduct(dip, normal) == doctest::Approx(0.0).epsilon(Eps));
    CHECK(dotProduct(strike, dip) == doctest::Approx(0.0).epsilon(Eps));
  }

  SUBCASE("Arbitrary oblique normal") {
    double normal[3] = {1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)};
    double strike[3] = {};
    double dip[3] = {};

    computeStrikeAndDipVectors(normal, strike, dip);

    CHECK(vecLength(strike) == doctest::Approx(1.0).epsilon(Eps));
    CHECK(vecLength(dip) == doctest::Approx(1.0).epsilon(Eps));
    CHECK(dotProduct(strike, normal) == doctest::Approx(0.0).epsilon(Eps));
    CHECK(dotProduct(dip, normal) == doctest::Approx(0.0).epsilon(Eps));
    CHECK(dotProduct(strike, dip) == doctest::Approx(0.0).epsilon(Eps));
  }

  SUBCASE("Vertical normal [0,0,1]") {
    double normal[3] = {0.0, 0.0, 1.0};
    double strike[3] = {};
    double dip[3] = {};

    computeStrikeAndDipVectors(normal, strike, dip);

    CHECK(std::isnan(vecLength(strike)));
    CHECK(std::isnan(vecLength(dip)));
  }
}

// ---------------------------------------------------------------------------
// quantity_indices sanity checks
// ---------------------------------------------------------------------------

TEST_CASE("DR Misc quantity_indices" * doctest::test_suite("dynamicrupture")) {
  using namespace quantity_indices;

  // Velocities are at indices 6,7,8
  CHECK(U == 6);
  CHECK(V == 7);
  CHECK(W == 8);

  // Normal stress overlaps with XX
  CHECK(N == XX);
  CHECK(N == 0);

  // Traction1 overlaps with XY
  CHECK(T1 == XY);
  CHECK(T1 == 3);

  // Traction2 overlaps with XZ
  CHECK(T2 == XZ);
  CHECK(T2 == 5);
}

} // namespace seissol::unit_test
