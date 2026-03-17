// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Equations/anisotropic/Model/Datastructures.h"
#include "Equations/elastic/Model/Datastructures.h"

#include <cmath>
#include <vector>

namespace seissol::unit_test {
using seissol::model::AnisotropicMaterial;
using seissol::model::ElasticMaterial;
using seissol::model::MaterialType;

TEST_CASE("AnisotropicMaterial default construction" * doctest::test_suite("equations")) {
  AnisotropicMaterial m;
  CHECK(m.rho == 0.0);
  CHECK(m.c11 == 0.0);
}

TEST_CASE("AnisotropicMaterial static properties" * doctest::test_suite("equations")) {
  CHECK(AnisotropicMaterial::NumQuantities == 9);
  CHECK(AnisotropicMaterial::Type == MaterialType::Anisotropic);
  CHECK(AnisotropicMaterial::Text == "anisotropic");
  CHECK(AnisotropicMaterial::Parameters == 22);
}

TEST_CASE("AnisotropicMaterial from ElasticMaterial" * doctest::test_suite("equations")) {
  double rho = 2700.0, mu = 3.24e10, lambda = 3.24e10;
  ElasticMaterial em({rho, mu, lambda});
  AnisotropicMaterial am(em);

  SUBCASE("Density preserved") { CHECK(am.rho == doctest::Approx(rho)); }

  SUBCASE("Isotropic diagonal: c_ii = lambda+2*mu") {
    CHECK(am.c11 == doctest::Approx(lambda + 2 * mu));
    CHECK(am.c22 == doctest::Approx(lambda + 2 * mu));
    CHECK(am.c33 == doctest::Approx(lambda + 2 * mu));
  }

  SUBCASE("Off-diagonal normal: c_ij = lambda") {
    CHECK(am.c12 == doctest::Approx(lambda));
    CHECK(am.c13 == doctest::Approx(lambda));
    CHECK(am.c23 == doctest::Approx(lambda));
  }

  SUBCASE("Shear: c_44=c_55=c_66 = mu") {
    CHECK(am.c44 == doctest::Approx(mu));
    CHECK(am.c55 == doctest::Approx(mu));
    CHECK(am.c66 == doctest::Approx(mu));
  }

  SUBCASE("Cross terms zero") {
    CHECK(am.c14 == doctest::Approx(0.0));
    CHECK(am.c15 == doctest::Approx(0.0));
    CHECK(am.c16 == doctest::Approx(0.0));
    CHECK(am.c45 == doctest::Approx(0.0));
    CHECK(am.c56 == doctest::Approx(0.0));
  }
}

TEST_CASE("AnisotropicMaterial averaged properties match elastic" *
          doctest::test_suite("equations")) {
  double rho = 2700.0, mu = 3.24e10, lambda = 3.24e10;
  ElasticMaterial em({rho, mu, lambda});
  AnisotropicMaterial am(em);

  CHECK(am.getMuBar() == doctest::Approx(mu));
  CHECK(am.getLambdaBar() == doctest::Approx(lambda));
  CHECK(am.getPWaveSpeed() == doctest::Approx(em.getPWaveSpeed()));
  CHECK(am.getSWaveSpeed() == doctest::Approx(em.getSWaveSpeed()));
}

TEST_CASE("AnisotropicMaterial getMaterialType" * doctest::test_suite("equations")) {
  AnisotropicMaterial m;
  CHECK(m.getMaterialType() == MaterialType::Anisotropic);
}

TEST_CASE("AnisotropicMaterial getMaxWaveSpeed isotropic case" * doctest::test_suite("equations")) {
  ElasticMaterial em({2700.0, 3.24e10, 3.24e10});
  AnisotropicMaterial am(em);
  double vmax = am.getMaxWaveSpeed();
  CHECK(vmax == doctest::Approx(em.getPWaveSpeed()).epsilon(0.01));
  CHECK(vmax > 0.0);
}

} // namespace seissol::unit_test
