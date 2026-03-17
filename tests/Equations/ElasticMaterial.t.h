// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Equations/elastic/Model/Datastructures.h"

#include <array>
#include <cmath>
#include <vector>

namespace seissol::unit_test {
using seissol::model::ElasticMaterial;
using seissol::model::MaterialType;

TEST_CASE("ElasticMaterial default construction" * doctest::test_suite("equations")) {
  ElasticMaterial m;
  CHECK(m.rho == 0.0);
  CHECK(m.lambda == 0.0);
  CHECK(m.mu == 0.0);
}

TEST_CASE("ElasticMaterial vector construction" * doctest::test_suite("equations")) {
  std::vector<double> vals = {2700.0, 3.24e10, 3.24e10};
  ElasticMaterial m(vals);
  CHECK(m.rho == doctest::Approx(2700.0));
  CHECK(m.mu == doctest::Approx(3.24e10));
  CHECK(m.lambda == doctest::Approx(3.24e10));
}

TEST_CASE("ElasticMaterial static properties" * doctest::test_suite("equations")) {
  CHECK(ElasticMaterial::NumQuantities == 9);
  CHECK(ElasticMaterial::Mechanisms == 0);
  CHECK(ElasticMaterial::Type == MaterialType::Elastic);
  CHECK(ElasticMaterial::Text == "elastic");
  CHECK(ElasticMaterial::SupportsDR == true);
  CHECK(ElasticMaterial::SupportsLTS == true);
}

TEST_CASE("ElasticMaterial wave speeds" * doctest::test_suite("equations")) {
  std::vector<double> vals = {2700.0, 3.24e10, 3.24e10};
  ElasticMaterial m(vals);
  double vp = m.getPWaveSpeed();
  double vs = m.getSWaveSpeed();

  SUBCASE("P-wave = sqrt((lambda+2*mu)/rho)") {
    CHECK(vp == doctest::Approx(std::sqrt((3.24e10 + 2.0 * 3.24e10) / 2700.0)));
  }
  SUBCASE("S-wave = sqrt(mu/rho)") { CHECK(vs == doctest::Approx(std::sqrt(3.24e10 / 2700.0))); }
  SUBCASE("P > S") { CHECK(vp > vs); }
  SUBCASE("MaxWaveSpeed = P") { CHECK(m.getMaxWaveSpeed() == doctest::Approx(vp)); }
  SUBCASE("Both positive") {
    CHECK(vp > 0.0);
    CHECK(vs > 0.0);
  }
}

TEST_CASE("ElasticMaterial getMaterialType" * doctest::test_suite("equations")) {
  ElasticMaterial m;
  CHECK(m.getMaterialType() == MaterialType::Elastic);
}

TEST_CASE("ElasticMaterial getLambdaBar getMuBar" * doctest::test_suite("equations")) {
  std::vector<double> vals = {2000.0, 1e10, 2e10};
  ElasticMaterial m(vals);
  CHECK(m.getLambdaBar() == doctest::Approx(2e10));
  CHECK(m.getMuBar() == doctest::Approx(1e10));
}

TEST_CASE("ElasticMaterial setLameParameters" * doctest::test_suite("equations")) {
  ElasticMaterial m;
  m.rho = 1000.0;
  m.setLameParameters(5e9, 8e9);
  CHECK(m.mu == doctest::Approx(5e9));
  CHECK(m.lambda == doctest::Approx(8e9));
  CHECK(m.getPWaveSpeed() == doctest::Approx(std::sqrt((8e9 + 2.0 * 5e9) / 1000.0)));
}

TEST_CASE("ElasticMaterial stiffness tensor" * doctest::test_suite("equations")) {
  double lambda = 3.0e10, mu = 2.0e10;
  std::vector<double> vals = {2700.0, mu, lambda};
  ElasticMaterial m(vals);
  std::array<double, 81> tensor{};
  m.getFullStiffnessTensor(tensor);
  auto idx = [](int i, int j, int k, int l) { return i * 27 + j * 9 + k * 3 + l; };

  SUBCASE("C_iiii = lambda+2*mu") {
    CHECK(tensor[idx(0, 0, 0, 0)] == doctest::Approx(lambda + 2 * mu));
    CHECK(tensor[idx(1, 1, 1, 1)] == doctest::Approx(lambda + 2 * mu));
    CHECK(tensor[idx(2, 2, 2, 2)] == doctest::Approx(lambda + 2 * mu));
  }
  SUBCASE("C_iijj = lambda") {
    CHECK(tensor[idx(0, 0, 1, 1)] == doctest::Approx(lambda));
    CHECK(tensor[idx(1, 1, 0, 0)] == doctest::Approx(lambda));
  }
  SUBCASE("C_ijij = mu") {
    CHECK(tensor[idx(0, 1, 0, 1)] == doctest::Approx(mu));
    CHECK(tensor[idx(0, 2, 0, 2)] == doctest::Approx(mu));
    CHECK(tensor[idx(1, 2, 1, 2)] == doctest::Approx(mu));
  }
  SUBCASE("Zero entries") { CHECK(tensor[idx(0, 0, 1, 2)] == doctest::Approx(0.0)); }
}

} // namespace seissol::unit_test
