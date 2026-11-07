// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/ConfigRegistry.h"
#include "Expr/Program.h"
#include "Expr/SderivFrontend.h"
#include "Physics/ScriptScope.h"

#include <map>
#include <set>
#include <string>
#include <vector>

namespace seissol::unit_test::script_scope {

/// The names of two configurations of a run, of a solid (0) and of a fluid (1), by rule; the ids
/// are not dispatched, so they need not be built.
inline physics::KnowsName namesOf(const std::map<ConfigId, std::set<std::string>>& names) {
  return [names](ConfigId config, const std::string& name) {
    return names.at(config).count(name) > 0;
  };
}

inline std::vector<std::string> outputsOf(const expr::Program& program) {
  std::vector<std::string> names;
  for (const auto& output : program.outputs()) {
    names.push_back(output.name);
  }
  return names;
}

TEST_CASE("ScriptScope: a script leaves out in a configuration what is of another material" *
          doctest::test_suite("physics")) {
  const auto provides = namesOf({{0, {"rho", "mu", "lambda"}}, {1, {"rho", "lambda"}}});
  const auto takes = namesOf({{0, {"s_xx", "v1"}}, {1, {"pprime", "v1"}}});
  const physics::ScriptScope scope({0, 1}, provides, takes);

  // a name of the other configuration only is left out; a name no configuration knows stays, for
  // the consumer to report
  CHECK(scope.keepsInput(0, "mu"));
  CHECK_FALSE(scope.keepsInput(1, "mu"));
  CHECK(scope.keepsInput(1, "rho"));
  CHECK(scope.keepsInput(1, "x"));
  CHECK(scope.keepsOutput(0, "s_xx"));
  CHECK_FALSE(scope.keepsOutput(1, "s_xx"));
  CHECK_FALSE(scope.keepsOutput(0, "pprime"));
  CHECK(scope.keepsOutput(1, "v1"));
  CHECK(scope.keepsOutput(1, "frame"));

  const auto program = expr::compileSderivModule(
      "out def s_xx = mu * 2.0\nout def pprime = -lambda * x\nout def v1 = rho * x\n"
      "out def frame = 1.0\n");
  CHECK(outputsOf(scope.restrict(program, 0, "test", false)) ==
        std::vector<std::string>{"s_xx", "v1", "frame"});
  CHECK(outputsOf(scope.restrict(program, 1, "test", false)) ==
        std::vector<std::string>{"pprime", "v1", "frame"});

  // the velocity of the fluid would read the shear modulus, which it has not
  const auto shear = expr::compileSderivModule("out def v1 = mu + rho\n");
  CHECK(outputsOf(scope.restrict(shear, 0, "test", false)) == std::vector<std::string>{"v1"});
  CHECK(scope.restrict(shear, 1, "test", false).outputs().empty());
  // kept in the solid, so no error
  scope.checkKept(shear, "test");

  // the filter is the restriction
  const auto filter = scope.filter(1, "test", false);
  CHECK(outputsOf(filter(program)) == outputsOf(scope.restrict(program, 1, "test", false)));
}

TEST_CASE("ScriptScope: a run of one configuration runs a script as it is" *
          doctest::test_suite("physics")) {
  const physics::ScriptScope scope(
      {0}, namesOf({{0, {"rho", "lambda"}}}), namesOf({{0, {"pprime", "v1"}}}));
  CHECK(scope.keepsInput(0, "mu"));
  CHECK(scope.keepsOutput(0, "s_xx"));
  const auto program = expr::compileSderivModule("out def s_xx = mu * 2.0\nout def v1 = rho\n");
  CHECK(outputsOf(scope.restrict(program, 0, "test", false)) ==
        std::vector<std::string>{"s_xx", "v1"});
}

} // namespace seissol::unit_test::script_scope
