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
#include "Model/MaterialType.h"
#include "TestConfigs.h"

#include <limits>

namespace seissol::unit_test {
using seissol::model::MaterialType;

TEST_CASE_TEMPLATE("MaterialTypeSelector consistency" * doctest::test_suite("equations"),
                   Cfg,
                   SEISSOL_CONFIG_TYPES) {
  using MaterialT = model::MaterialOf<Cfg>;
  CHECK(MaterialT::Type == Cfg::MaterialType);
}

TEST_CASE_TEMPLATE_DEFINE("MaterialT static properties valid" * doctest::test_suite("equations"),
                          Cfg,
                          MaterialTStaticPropertiesValid) {
  using MaterialT = model::MaterialOf<Cfg>;
  CHECK(MaterialT::NumQuantities > 0);
  CHECK(MaterialT::Quantities.size() > 0);
  CHECK_FALSE(MaterialT::Text.empty());
  CHECK(MaterialT::Parameters >= 1);
}

TEST_CASE_TEMPLATE_APPLY(MaterialTStaticPropertiesValid, MaterialConfigs);

TEST_CASE_TEMPLATE_DEFINE("MaterialT default density zero" * doctest::test_suite("equations"),
                          Cfg,
                          MaterialTDefaultDensityZero) {
  using MaterialT = model::MaterialOf<Cfg>;
  MaterialT m;
  CHECK(m.getDensity() == 0.0);
  m.setDensity(2700.0);
  CHECK(m.getDensity() == doctest::Approx(2700.0));
}

TEST_CASE_TEMPLATE_APPLY(MaterialTDefaultDensityZero, MaterialConfigs);

TEST_CASE_TEMPLATE_DEFINE("MaterialT getMaterialType matches static" *
                              doctest::test_suite("equations"),
                          Cfg,
                          MaterialTGetMaterialTypeMatchesStatic) {
  using MaterialT = model::MaterialOf<Cfg>;
  const MaterialT m{};
  CHECK(m.getMaterialType() == MaterialT::Type);
}

TEST_CASE_TEMPLATE_APPLY(MaterialTGetMaterialTypeMatchesStatic, MaterialConfigs);

TEST_CASE_TEMPLATE_DEFINE("MaterialT maximumTimestep default infinity" *
                              doctest::test_suite("equations"),
                          Cfg,
                          MaterialTMaximumTimestepDefaultInfinity) {
  using MaterialT = model::MaterialOf<Cfg>;
  const MaterialT m{};
  CHECK(m.maximumTimestep() == std::numeric_limits<double>::infinity());
}

TEST_CASE_TEMPLATE_APPLY(MaterialTMaximumTimestepDefaultInfinity, MaterialConfigs);

TEST_CASE_TEMPLATE("Config check" * doctest::test_suite("equations"), Cfg, SEISSOL_CONFIG_TYPES) {
  using MaterialT = model::MaterialOf<Cfg>;
  if (Cfg::MaterialType == MaterialType::Elastic) {
    CHECK(MaterialT::NumQuantities == 9);
    CHECK(MaterialT::SupportsDR == true);
    CHECK(MaterialT::Mechanisms == 0);
  }
  if (Cfg::MaterialType == MaterialType::Acoustic) {
    CHECK(MaterialT::NumQuantities == 4);
    CHECK(MaterialT::SupportsDR == false);
    CHECK(MaterialT::Mechanisms == 0);
  }
  if (Cfg::MaterialType == MaterialType::Anisotropic) {
    CHECK(MaterialT::NumQuantities == 9);
    CHECK(MaterialT::SupportsDR == true);
    CHECK(MaterialT::Parameters == 22);
    CHECK(MaterialT::Mechanisms == 0);
  }
  if (Cfg::MaterialType == MaterialType::Viscoelastic) {
    CHECK(MaterialT::NumQuantities == 9 + 6 * MaterialT::Mechanisms);
    CHECK(MaterialT::Mechanisms > 0);
    CHECK(MaterialT::SupportsDR == true);
    CHECK(MaterialT::Mechanisms == Cfg::RelaxationMechanisms);
  }
  if (Cfg::MaterialType == MaterialType::Poroelastic) {
    CHECK(MaterialT::NumQuantities == 13);
    CHECK(MaterialT::SupportsDR == true);
    CHECK(MaterialT::Mechanisms == 0);
  }
  if (Cfg::MaterialType == MaterialType::Viscoacoustic) {
    CHECK(MaterialT::NumQuantities == 4 + 1 * MaterialT::Mechanisms);
    CHECK(MaterialT::Mechanisms > 0);
    CHECK(MaterialT::SupportsDR == false);
    CHECK(MaterialT::Mechanisms == Cfg::RelaxationMechanisms);
  }
}

} // namespace seissol::unit_test
