// SPDX-FileCopyrightText: 2017-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Numerical/Functions.h"
#include "doctest.h"
#include "tests/TestHelper.h"
#include <Kernels/Precision.h>

namespace seissol::unit_test {

TEST_CASE("Test Jacobi polynomials") {
  constexpr double Epsilon = 10 * std::numeric_limits<real>::epsilon();
  // Compare to Maple reference solution
  REQUIRE(seissol::functions::JacobiP(0, 1, 0, 0.5) == AbsApprox(1.0).epsilon(Epsilon));
  REQUIRE(seissol::functions::JacobiP(1, 0, 0, -0.3) ==
          AbsApprox(-.3000000000000000).epsilon(Epsilon));
  REQUIRE(seissol::functions::JacobiP(1, 13, 4, -1.0) ==
          AbsApprox(-5.000000000000000).epsilon(Epsilon));
  REQUIRE(seissol::functions::JacobiP(10, 0, 4, -1.0) == AbsApprox(1001.).epsilon(Epsilon));
  REQUIRE(seissol::functions::JacobiP(10, 0, 4, -0.3) ==
          AbsApprox(-1.870307686685156).epsilon(Epsilon));
  REQUIRE(seissol::functions::JacobiP(10, 0, 4, 0.0) ==
          AbsApprox(0.8671875000000000).epsilon(Epsilon));
  REQUIRE(seissol::functions::JacobiP(10, 0, 4, 0.3) ==
          AbsApprox(-.3404032083257812).epsilon(Epsilon));
  REQUIRE(seissol::functions::JacobiP(10, 0, 4, 1.0) == AbsApprox(1.).epsilon(Epsilon));
  REQUIRE(seissol::functions::JacobiP(2, 1, 4, -1.0) == AbsApprox(15.).epsilon(Epsilon));
  REQUIRE(seissol::functions::JacobiP(3, 2, 3, -0.3) ==
          AbsApprox(1.249375000000000).epsilon(Epsilon));
  REQUIRE(seissol::functions::JacobiP(4, 3, 2, 0.0) ==
          AbsApprox(0.9375000000000002).epsilon(Epsilon));
  REQUIRE(seissol::functions::JacobiP(5, 4, 1, 0.3) ==
          AbsApprox(-1.982012812500000).epsilon(Epsilon));
  REQUIRE(seissol::functions::JacobiP(6, 5, 0, 1.0) ==
          AbsApprox(462.0000000000000).epsilon(Epsilon));
}

TEST_CASE("Test first derivative of Jacobi polynomials") {
  constexpr double Epsilon = 100 * std::numeric_limits<real>::epsilon();
  // Compare to Maple reference solution
  REQUIRE(seissol::functions::JacobiPDerivative(0, 1, 0, 0.5) == AbsApprox(0.0).epsilon(Epsilon));
  REQUIRE(seissol::functions::JacobiPDerivative(1, 0, 0, -0.3) ==
          AbsApprox(1.000000000000000).epsilon(Epsilon));
  REQUIRE(seissol::functions::JacobiPDerivative(1, 13, 4, -0.9) ==
          AbsApprox(9.500000000000000).epsilon(Epsilon));
  REQUIRE(seissol::functions::JacobiPDerivative(7, 0, 4, -0.3) ==
          AbsApprox(-9.122014125000001).epsilon(Epsilon));
  REQUIRE(seissol::functions::JacobiPDerivative(7, 0, 4, 0.0) ==
          AbsApprox(7.8750000000000).epsilon(Epsilon));
  REQUIRE(seissol::functions::JacobiPDerivative(7, 0, 4, 0.3) ==
          AbsApprox(-6.067144125000002).epsilon(Epsilon));
  REQUIRE(seissol::functions::JacobiPDerivative(7, 0, 4, 0.9) ==
          AbsApprox(-2.71497712500000).epsilon(Epsilon));
  REQUIRE(seissol::functions::JacobiPDerivative(2, 1, 4, -0.9) ==
          AbsApprox(-22.20000000000000).epsilon(Epsilon));
  REQUIRE(seissol::functions::JacobiPDerivative(3, 2, 3, -0.3) ==
          AbsApprox(3.318750000000000).epsilon(Epsilon));
  REQUIRE(seissol::functions::JacobiPDerivative(4, 3, 2, 0.0) ==
          AbsApprox(-3.749999999999998).epsilon(Epsilon));
  REQUIRE(seissol::functions::JacobiPDerivative(5, 4, 1, 0.3) ==
          AbsApprox(-15.37232812500000).epsilon(Epsilon));
}

TEST_CASE("Test Dubiner polynomials") {
  constexpr double Epsilon = 10 * std::numeric_limits<real>::epsilon();
  const std::vector<std::array<double, 7>> tests = {
      {0, 0, 0, 0.25, 0.25, 0.25, 1},
      {1, 0, 0, 0.25, 0.25, 0.25, 0.},
      {0, 1, 0, 0.25, 0.25, 0.25, 0.},
      {0, 0, 1, 0.25, 0.25, 0.25, 0.},
      {2, 0, 0, 0.25, 0.25, 0.25, -0.1250},
      {1, 1, 0, 0.25, 0.25, 0.25, 0.},
      {0, 2, 0, 0.25, 0.25, 0.25, -0.3125},
      {1, 0, 1, 0.25, 0.25, 0.25, 0.},
      {0, 1, 1, 0.25, 0.25, 0.25, 0.},
      {0, 0, 2, 0.25, 0.25, 0.25, -0.5625},
      {3, 0, 0, 0.25, 0.25, 0.25, -0.},
      {2, 1, 0, 0.25, 0.25, 0.25, -0.125000},
      {1, 2, 0, 0.25, 0.25, 0.25, -0.},
      {0, 3, 0, 0.25, 0.25, 0.25, 0.125000},
      {2, 0, 1, 0.25, 0.25, 0.25, -0.125000},
      {1, 1, 1, 0.25, 0.25, 0.25, 0.},
      {0, 2, 1, 0.25, 0.25, 0.25, -0.31250000000000000000},
      {1, 0, 2, 0.25, 0.25, 0.25, -0.},
      {0, 1, 2, 0.25, 0.25, 0.25, -0.},
      {0, 0, 3, 0.25, 0.25, 0.25, 0.437500},
      {0, 0, 0, 0., 0.5, 0.5, 1.},
      {1, 0, 0, 0., 0.5, 0.5, 0.},
      {0, 1, 0, 0., 0.5, 0.5, 1.0},
      {0, 0, 1, 0., 0.5, 0.5, 1.0},
      {2, 0, 0, 0., 0.5, 0.5, 0.},
      {1, 1, 0, 0., 0.5, 0.5, 0.},
      {0, 2, 0, 0., 0.5, 0.5, 0.75},
      {1, 0, 1, 0., 0.5, 0.5, 0.},
      {0, 1, 1, 0., 0.5, 0.5, 2.00},
      {0, 0, 2, 0., 0.5, 0.5, -0.25},
      {3, 0, 0, 0., 0.5, 0.5, 0.},
      {2, 1, 0, 0., 0.5, 0.5, 0.},
      {1, 2, 0, 0., 0.5, 0.5, 0.},
      {0, 3, 0, 0., 0.5, 0.5, 0.500},
      {2, 0, 1, 0., 0.5, 0.5, 0.},
      {1, 1, 1, 0., 0.5, 0.5, 0.},
      {0, 2, 1, 0., 0.5, 0.5, 2.2500000000000000000},
      {1, 0, 2, 0., 0.5, 0.5, 0.},
      {0, 1, 2, 0., 0.5, 0.5, 1.0},
      {0, 0, 3, 0., 0.5, 0.5, -0.750},
      {0, 0, 0, 0., 1.0, 0., 1.},
      {1, 0, 0, 0., 1.0, 0., 0.},
      {0, 1, 0, 0., 1.0, 0., 2.0},
      {0, 0, 1, 0., 1.0, 0., -1.},
      {2, 0, 0, 0., 1.0, 0., 0.},
      {1, 1, 0, 0., 1.0, 0., 0.},
      {0, 2, 0, 0., 1.0, 0., 3.00},
      {1, 0, 1, 0., 1.0, 0., -0.},
      {0, 1, 1, 0., 1.0, 0., -2.0},
      {0, 0, 2, 0., 1.0, 0., 1.},
      {3, 0, 0, 0., 1.0, 0., 0.},
      {2, 1, 0, 0., 1.0, 0., 0.},
      {1, 2, 0, 0., 1.0, 0., 0.},
      {0, 3, 0, 0., 1.0, 0., 4.000},
      {2, 0, 1, 0., 1.0, 0., -0.},
      {1, 1, 1, 0., 1.0, 0., -0.},
      {0, 2, 1, 0., 1.0, 0., -3.0000000000000000000},
      {1, 0, 2, 0., 1.0, 0., 0.},
      {0, 1, 2, 0., 1.0, 0., 2.0},
      {0, 0, 3, 0., 1.0, 0., -1.}};
  for (const auto& t : tests) {
    REQUIRE(
        seissol::functions::TetraDubinerP(
            {static_cast<unsigned>(t[0]), static_cast<unsigned>(t[1]), static_cast<unsigned>(t[2])},
            {t[3], t[4], t[5]}) == AbsApprox(t[6]).epsilon(Epsilon));
  }
}

} // namespace seissol::unit_test
