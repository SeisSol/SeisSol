// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Initializer/Parameters/InitializationParameters.h"
#include "Initializer/Parameters/ParameterReader.h"

#include <yaml-cpp/yaml.h>

namespace seissol::unit_test {
using namespace seissol::initializer::parameters;

TEST_CASE("readInitializationParameters accepts the acoustic travelling wave with ITM" *
          doctest::test_suite("initializer")) {
  // spelled as in the documentation; the parser compares in lower case
  const YAML::Node node = YAML::Load(R"(
    inicondition:
      cictype: AcousticTravellingwithITM
      k: 6.283
  )");
  ParameterReader reader(node, "", false);
  const auto params = readInitializationParameters(&reader);

  CHECK(params.type == InitializationType::AcousticTravellingWithITM);
  CHECK(params.k == doctest::Approx(6.283));
}

} // namespace seissol::unit_test
