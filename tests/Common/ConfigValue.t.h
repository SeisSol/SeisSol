// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Common/ConfigValue.h"
#include "Common/Real.h"
#include "Common/Typedefs.h"
#include "Model/MaterialType.h"

#include <array>
#include <cstddef>
#include <set>
#include <string>

namespace seissol::unit_test {

using namespace seissol;

TEST_CASE("Configuration names" * doctest::test_suite("common")) {
  SUBCASE("Every field appears in the name") {
    const ConfigValue elastic{6,
                              0,
                              model::MaterialType::Elastic,
                              RealType::F64,
                              SolverType::LinearCK,
                              DRQuadRuleType::Stroud,
                              1};
    CHECK(configName(elastic) == "elastic-linearck-o6-f64-stroud");

    const ConfigValue viscoelastic{4,
                                   3,
                                   model::MaterialType::Viscoelastic,
                                   RealType::F32,
                                   SolverType::LinearCKAnelastic,
                                   DRQuadRuleType::Dunavant,
                                   8};
    CHECK(configName(viscoelastic) == "viscoelastic-linearckanelastic-m3-o4-f32-dunavant-s8");
  }

  SUBCASE("Distinct values have distinct names, which read back") {
    constexpr std::array Materials{model::MaterialType::Solid,
                                   model::MaterialType::Acoustic,
                                   model::MaterialType::Elastic,
                                   model::MaterialType::Viscoelastic,
                                   model::MaterialType::Viscoacoustic,
                                   model::MaterialType::Anisotropic,
                                   model::MaterialType::Poroelastic};
    constexpr std::array Solvers{
        SolverType::LinearCK, SolverType::LinearCKAnelastic, SolverType::STP};
    constexpr std::array Rules{
        DRQuadRuleType::Stroud, DRQuadRuleType::Dunavant, DRQuadRuleType::WitherdenVincent};
    constexpr std::array Precisions{RealType::F32, RealType::F64};
    constexpr std::array<std::size_t, 4> Mechanisms{0, 1, 3, 10};
    constexpr std::array<std::size_t, 4> Orders{1, 2, 6, 11};
    constexpr std::array<std::size_t, 4> Simulations{0, 1, 2, 16};

    std::set<std::string> names;
    std::size_t count = 0;
    for (const auto material : Materials) {
      for (const auto solver : Solvers) {
        for (const auto mechanisms : Mechanisms) {
          for (const auto order : Orders) {
            for (const auto precision : Precisions) {
              for (const auto rule : Rules) {
                for (const auto simulations : Simulations) {
                  const ConfigValue value{
                      order, mechanisms, material, precision, solver, rule, simulations};
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
      }
    }
    CHECK(names.size() == count);
  }

  SUBCASE("Only the canonical spelling names a configuration") {
    CHECK(parseConfigName("elastic-linearck-o6-f64-stroud").has_value());

    for (const auto* other : {"",
                              "elastic",
                              "elastic-linearck",
                              "elastic-linearck-o6-f64",
                              "elastic-linearck-o6-f64-stroud-",
                              "-elastic-linearck-o6-f64-stroud",
                              "Elastic-linearck-o6-f64-stroud",
                              "elastic-linearck-o06-f64-stroud",
                              "elastic-linearck-o+6-f64-stroud",
                              "elastic-linearck-o-6-f64-stroud",
                              "elastic-linearck-m0-o6-f64-stroud",
                              "elastic-linearck-o6-f64-stroud-s1",
                              "elastic-linearck-o6-f64-stroud-s2-s2",
                              "elastic-linearck-f64-o6-stroud",
                              "elastic-linearck-o6-f16-stroud",
                              "elastic-linearck-o6-f64-gauss",
                              "elastic-linearck-o99999999999999999999999-f64-stroud",
                              "elastic-stroud-o6-f64-linearck"}) {
      CAPTURE(other);
      CHECK_FALSE(parseConfigName(other).has_value());
    }
  }

  SUBCASE("The description lists every field") {
    const ConfigValue poroelastic{3,
                                  0,
                                  model::MaterialType::Poroelastic,
                                  RealType::F32,
                                  SolverType::STP,
                                  DRQuadRuleType::WitherdenVincent,
                                  1};
    CHECK(describeConfig(poroelastic) == "Material: poroelastic\n"
                                         "Solver: stp\n"
                                         "Relaxation mechanisms: 0\n"
                                         "Convergence order: 3\n"
                                         "Precision: f32\n"
                                         "Dynamic rupture quadrature rule: witherdenvincent\n"
                                         "Fused simulations: 1\n");
  }
}

} // namespace seissol::unit_test
