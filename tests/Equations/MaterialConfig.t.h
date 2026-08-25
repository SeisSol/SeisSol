// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Config.h"
#include "Equations/Datastructures.h"
#include "Model/CommonDatastructures.h"

#include <limits>

namespace seissol::unit_test {
using seissol::model::MaterialT;
using seissol::model::MaterialType;

TEST_CASE("MaterialTypeSelector consistency" * doctest::test_suite("equations")) {
  CHECK(MaterialT::Type == seissol::Config::MaterialType);
}

TEST_CASE("MaterialT static properties valid" * doctest::test_suite("equations")) {
  CHECK(MaterialT::NumQuantities > 0);
  CHECK(MaterialT::Quantities.size() > 0);
  CHECK_FALSE(MaterialT::Text.empty());
  CHECK(MaterialT::Parameters >= 1);
}

TEST_CASE("MaterialT default density zero" * doctest::test_suite("equations")) {
  MaterialT m;
  CHECK(m.getDensity() == 0.0);
  m.setDensity(2700.0);
  CHECK(m.getDensity() == doctest::Approx(2700.0));
}

TEST_CASE("MaterialT getMaterialType matches static" * doctest::test_suite("equations")) {
  MaterialT m;
  CHECK(m.getMaterialType() == MaterialT::Type);
}

TEST_CASE("MaterialT maximumTimestep default infinity" * doctest::test_suite("equations")) {
  MaterialT m;
  CHECK(m.maximumTimestep() == std::numeric_limits<double>::infinity());
}

TEST_CASE("Config check" * doctest::test_suite("equations")) {
  if (seissol::Config::MaterialType == MaterialType::Elastic) {
    CHECK(MaterialT::NumQuantities == 9);
    CHECK(MaterialT::SupportsDR == true);
    CHECK(MaterialT::Mechanisms == 0);
  }
  if (seissol::Config::MaterialType == MaterialType::Acoustic) {
    CHECK(MaterialT::NumQuantities == 4);
    CHECK(MaterialT::SupportsDR == false);
    CHECK(MaterialT::Mechanisms == 0);
  }
  if (seissol::Config::MaterialType == MaterialType::Anisotropic) {
    CHECK(MaterialT::NumQuantities == 9);
    CHECK(MaterialT::SupportsDR == false);
    CHECK(MaterialT::Parameters == 22);
    CHECK(MaterialT::Mechanisms == 0);
  }
  if (seissol::Config::MaterialType == MaterialType::Viscoelastic) {
    CHECK(MaterialT::NumQuantities == 9 + 6 * MaterialT::Mechanisms);
    CHECK(MaterialT::Mechanisms > 0);
    CHECK(MaterialT::SupportsDR == true);
    CHECK(MaterialT::Mechanisms == seissol::Config::RelaxationMechanisms);
  }
  if (seissol::Config::MaterialType == MaterialType::Poroelastic) {
    CHECK(MaterialT::NumQuantities == 13);
    CHECK(MaterialT::SupportsDR == true);
    CHECK(MaterialT::Mechanisms == 0);
  }
}

} // namespace seissol::unit_test
