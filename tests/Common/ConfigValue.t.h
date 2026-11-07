// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Common/ConfigValue.h"
#include "Common/Real.h"
#include "Model/MaterialType.h"

#include <array>
#include <cstddef>
#include <set>
#include <string>

namespace seissol::unit_test {

using namespace seissol;

TEST_CASE("Configuration names" * doctest::test_suite("common")) {
  SUBCASE("Every field appears in the name") {
    const ConfigValue elastic{.convergenceOrder = 6,
                              .relaxationMechanisms = 0,
                              .materialType = model::MaterialType::Elastic,
                              .precision = RealType::F64,
                              .numSimulations = 1};
    CHECK(configName(elastic) == "elastic-o6-f64");

    const ConfigValue viscoelastic{.convergenceOrder = 4,
                                   .relaxationMechanisms = 3,
                                   .materialType = model::MaterialType::Viscoelastic,
                                   .precision = RealType::F32,
                                   .numSimulations = 8};
    CHECK(configName(viscoelastic) == "viscoelastic-m3-o4-f32-s8");
  }

  SUBCASE("Distinct values have distinct names, which read back") {
    constexpr std::array Materials{model::MaterialType::Solid,
                                   model::MaterialType::Acoustic,
                                   model::MaterialType::Elastic,
                                   model::MaterialType::Viscoelastic,
                                   model::MaterialType::Viscoacoustic,
                                   model::MaterialType::Anisotropic,
                                   model::MaterialType::Poroelastic};
    constexpr std::array Precisions{RealType::F32, RealType::F64};
    constexpr std::array<std::size_t, 4> Mechanisms{0, 1, 3, 10};
    constexpr std::array<std::size_t, 4> Orders{1, 2, 6, 11};
    constexpr std::array<std::size_t, 4> Simulations{0, 1, 2, 16};

    std::set<std::string> names;
    std::size_t count = 0;
    for (const auto material : Materials) {
      for (const auto mechanisms : Mechanisms) {
        for (const auto order : Orders) {
          for (const auto precision : Precisions) {
            for (const auto simulations : Simulations) {
              const ConfigValue value{.convergenceOrder = order,
                                      .relaxationMechanisms = mechanisms,
                                      .materialType = material,
                                      .precision = precision,
                                      .numSimulations = simulations};
              const auto name = configName(value);
              const auto parsed = parseConfigName(name);
              REQUIRE(parsed.has_value());
              CHECK(parsed.value() == value);
              names.insert(name);
              ++count;
            }
          }
        }
      }
    }
    CHECK(names.size() == count);
  }

  SUBCASE("Only the canonical spelling names a configuration") {
    CHECK(parseConfigName("elastic-o6-f64").has_value());
    CHECK(parseConfigName("viscoelastic-m3-o6-f64-s8").has_value());

    for (const auto* other : {"",
                              "elastic",
                              "elastic-o6",
                              "elastic-o6-f64-",
                              "-elastic-o6-f64",
                              "Elastic-o6-f64",
                              "elastic-o06-f64",
                              "elastic-o+6-f64",
                              "elastic-o-6-f64",
                              "elastic-m0-o6-f64",
                              "elastic-o6-f64-s1",
                              "elastic-o6-f64-s2-s2",
                              "elastic-f64-o6",
                              "elastic-o6-f16",
                              "elastic-o99999999999999999999999-f64",
                              "viscoelastic-3-o6-f64",
                              "viscoelastic-o6-m3-f64",
                              "elastic-linearck-o6-f64-stroud"}) {
      CAPTURE(other);
      CHECK_FALSE(parseConfigName(other).has_value());
    }
  }

  SUBCASE("The description lists every field") {
    const ConfigValue poroelastic{.convergenceOrder = 3,
                                  .relaxationMechanisms = 0,
                                  .materialType = model::MaterialType::Poroelastic,
                                  .precision = RealType::F32,
                                  .numSimulations = 1};
    CHECK(describeConfig(poroelastic) == "Material: poroelastic\n"
                                         "Relaxation mechanisms: 0\n"
                                         "Convergence order: 3\n"
                                         "Precision: f32\n"
                                         "Fused simulations: 1\n");
  }
}

} // namespace seissol::unit_test
