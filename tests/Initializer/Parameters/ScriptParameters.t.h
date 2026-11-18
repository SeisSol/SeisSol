// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Expr/Program.h"
#include "Initializer/Parameters/ParameterReader.h"
#include "Initializer/Parameters/ScriptParameters.h"

#include <string>
#include <utility>
#include <yaml-cpp/yaml.h>

namespace seissol::unit_test {
using namespace seissol::initializer::parameters;

TEST_CASE("readScriptParameters: the precision and the strictness of the scripts of a run" *
          doctest::test_suite("initializer")) {
  // without the section, as most parameter files are: each configuration's precision, contracted
  {
    const YAML::Node node = YAML::Load("equations: {}\n");
    ParameterReader reader(node, "", false);
    const auto parameters = readScriptParameters(&reader);
    CHECK(!parameters.precision.has_value());
    CHECK(!parameters.strict);
  }
  {
    const YAML::Node node = YAML::Load(R"(
      scripts:
        precision: 'auto'
    )");
    ParameterReader reader(node, "", false);
    CHECK(!readScriptParameters(&reader).precision.has_value());
  }
  for (const auto& [spelling, type] :
       {std::pair{"'f32'", expr::ComputeType::F32}, std::pair{"'f64'", expr::ComputeType::F64}}) {
    CAPTURE(spelling);
    const YAML::Node node =
        YAML::Load(std::string("scripts:\n  precision: ") + spelling + "\n  strict: 1\n");
    ParameterReader reader(node, "", false);
    const auto parameters = readScriptParameters(&reader);
    REQUIRE(parameters.precision.has_value());
    CHECK(*parameters.precision == type);
    CHECK(parameters.strict);
  }
}

} // namespace seissol::unit_test
