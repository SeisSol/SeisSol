// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Equations/Datastructures.h"
#include "Equations/viscoelastic2/Model/Datastructures.h"
#include "Physics/Attenuation.h"
#include "TestHelper.h"

#include <cstdlib>

namespace seissol::unit_test {
TEST_CASE("Attenuation") {
  seissol::model::ViscoElasticMaterialParametrized<3> vm;

  // The test data used here was picked from one cell in a coarse discretization of the Northridge
  // scenario.

  double rho = 2787.216864690955845;
  vm.rho = rho;
  vm.mu = 30969716301.94932938;
  vm.lambda = 30969716301.94934082;

  vm.qp = 58.40512125971850566;
  vm.qs = 29.20256062985925283;

  constexpr double FreqCentral = 0.3;
  constexpr double FreqRatio = 100;

  seissol::physics::fitAttenuation(vm, FreqCentral, FreqRatio);

  REQUIRE(vm.rho == rho);

  REQUIRE(vm.mu == doctest::Approx(33628014790.452877));
  REQUIRE(vm.lambda == doctest::Approx(29578827039.733589));

  REQUIRE(vm.omega[0] == doctest::Approx(0.18849555921538758));
  REQUIRE(vm.omega[1] == doctest::Approx(1.8849555921538763));
  REQUIRE(vm.omega[2] == doctest::Approx(18.849555921538762));

  REQUIRE(vm.theta[0][0] == doctest::Approx(-2662632131.4096899));
  REQUIRE(vm.theta[1][0] == doctest::Approx(-2233089891.3828750));
  REQUIRE(vm.theta[2][0] == doctest::Approx(-2833117392.9633269));
  REQUIRE(vm.theta[0][1] == doctest::Approx(777200013.84026122));
  REQUIRE(vm.theta[1][1] == doctest::Approx(739358274.28576028));
  REQUIRE(vm.theta[2][1] == doctest::Approx(1061064406.7129726));
  REQUIRE(vm.theta[0][2] == doctest::Approx(-3439832145.2499514));
  REQUIRE(vm.theta[1][2] == doctest::Approx(-2972448165.6686354));
  REQUIRE(vm.theta[2][2] == doctest::Approx(-3894181799.6762996));
}
} // namespace seissol::unit_test
