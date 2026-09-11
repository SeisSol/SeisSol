// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Equations/acoustic/Model/Datastructures.h"

#include <cmath>
#include <vector>

namespace seissol::unit_test {
using seissol::model::AcousticMaterial;
using seissol::model::MaterialType;

TEST_CASE("AcousticMaterial default construction" * doctest::test_suite("equations")) {
  AcousticMaterial m;
  CHECK(m.rho == 0.0);
  CHECK(m.lambda == 0.0);
}

TEST_CASE("AcousticMaterial vector construction" * doctest::test_suite("equations")) {
  const std::vector<double> vals = {1000.0, 2.25e9};
  AcousticMaterial m(vals);
  CHECK(m.rho == doctest::Approx(1000.0));
  CHECK(m.lambda == doctest::Approx(2.25e9));
}

TEST_CASE("AcousticMaterial static properties" * doctest::test_suite("equations")) {
  CHECK(AcousticMaterial::NumQuantities == 4);
  CHECK(AcousticMaterial::Type == MaterialType::Acoustic);
  CHECK(AcousticMaterial::Text == "acoustic");
  CHECK(AcousticMaterial::SupportsDR == false);
}

TEST_CASE("AcousticMaterial wave speeds" * doctest::test_suite("equations")) {
  const std::vector<double> vals = {1000.0, 2.25e9};
  AcousticMaterial m(vals);

  SUBCASE("P-wave = sqrt(lambda/rho) = 1500") {
    CHECK(m.getPWaveSpeed() == doctest::Approx(1500.0));
  }
  SUBCASE("S-wave is zero") { CHECK(m.getSWaveSpeed() == doctest::Approx(0.0)); }
  SUBCASE("MaxWaveSpeed = P") { CHECK(m.getMaxWaveSpeed() == doctest::Approx(m.getPWaveSpeed())); }
  SUBCASE("MuBar = 0") { CHECK(m.getMuBar() == doctest::Approx(0.0)); }
  SUBCASE("LambdaBar = lambda") { CHECK(m.getLambdaBar() == doctest::Approx(2.25e9)); }
}

TEST_CASE("AcousticMaterial getMaterialType" * doctest::test_suite("equations")) {
  AcousticMaterial m;
  CHECK(m.getMaterialType() == MaterialType::Acoustic);
}

TEST_CASE("AcousticMaterial setLameParameters" * doctest::test_suite("equations")) {
  AcousticMaterial m;
  m.rho = 1000.0;
  m.setLameParameters(0.0, 4.0e9);
  CHECK(m.lambda == doctest::Approx(4.0e9));
  CHECK(m.getPWaveSpeed() == doctest::Approx(std::sqrt(4.0e9 / 1000.0)));
}

} // namespace seissol::unit_test
