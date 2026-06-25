// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Initializer/Parameters/ModelParameters.h"
#include "Initializer/Parameters/ParameterReader.h"

#include <string>
#include <yaml-cpp/yaml.h>

namespace seissol::unit_test {
using namespace seissol::initializer::parameters;

// ---------------------------------------------------------------------------
// fluxToString
// ---------------------------------------------------------------------------

TEST_CASE("fluxToString" * doctest::test_suite("initializer")) {
  CHECK(fluxToString(NumericalFlux::Godunov) == "Godunov flux");
  CHECK(fluxToString(NumericalFlux::Rusanov) == "Rusanov flux");
}

// ---------------------------------------------------------------------------
// ITMParameters defaults
// ---------------------------------------------------------------------------

TEST_CASE("readITMParameters defaults" * doctest::test_suite("initializer")) {
  // Provide an empty equations section → all defaults should apply
  const YAML::Node node = YAML::Load(R"(
    equations: {}
  )");
  ParameterReader reader(node, "", false);
  auto itm = readITMParameters(&reader);

  CHECK(itm.itmEnabled == false);
  CHECK(itm.itmStartingTime == doctest::Approx(0.0));
  CHECK(itm.itmDuration == doctest::Approx(0.0));
  CHECK(itm.itmVelocityScalingFactor == doctest::Approx(1.0));
  CHECK(itm.itmReflectionType == ReflectionType::BothWaves);
}

TEST_CASE("readITMParameters custom values" * doctest::test_suite("initializer")) {
  const YAML::Node node = YAML::Load(R"(
    equations:
      itmenable: 0
      itmstartingtime: 2.5
      itmtime: 5.0
      itmvelocityscalingfactor: 0.8
      itmreflectiontype: 3
  )");
  ParameterReader reader(node, "", false);
  auto itm = readITMParameters(&reader);

  // Even though values are provided, itmEnabled=false means they're marked unused
  // but still parsed. We can verify the parsing happened.
  CHECK(itm.itmEnabled == false);
  CHECK(itm.itmStartingTime == doctest::Approx(2.5));
  CHECK(itm.itmDuration == doctest::Approx(5.0));
  CHECK(itm.itmVelocityScalingFactor == doctest::Approx(0.8));
  CHECK(itm.itmReflectionType == ReflectionType::Pwave);
}

// ---------------------------------------------------------------------------
// readModelParameters
// ---------------------------------------------------------------------------

TEST_CASE("readModelParameters parses YAML" * doctest::test_suite("initializer")) {
  // Provide a minimal-but-complete equations section.
  // materialfilename is required, so it must be present.
  const YAML::Node node = YAML::Load(R"(
    equations:
      materialfilename: material.yaml
      plasticity: 1
      plasticitypointwise: 0
      usecellhomogenizedmaterial: 0
      gravitationalacceleration: 10.0
      tv: 0.2
      numflux: Rusanov
      numfluxnearfault: Godunov
      qp: 1
      qs: 1
  )");
  ParameterReader reader(node, "", false);
  auto params = readModelParameters(&reader);

  CHECK(params.materialFileName == "material.yaml");
  CHECK(params.plasticity == true);
  CHECK(params.plasticityPointwise == false);
  CHECK(params.useCellHomogenizedMaterial == false);
  CHECK(params.gravitationalAcceleration == doctest::Approx(10.0));
  CHECK(params.tv == doctest::Approx(0.2));
  CHECK(params.flux == NumericalFlux::Rusanov);
  CHECK(params.fluxNearFault == NumericalFlux::Godunov);
  CHECK(params.hasBoundaryFile == false);
}

TEST_CASE("readModelParameters defaults" * doctest::test_suite("initializer")) {
  // (qp and qs are necessary for visco to not fail)

  const YAML::Node node = YAML::Load(R"(
    equations:
      materialfilename: mat.yaml
      qp: 1
      qs: 1
  )");
  ParameterReader reader(node, "", false);
  auto params = readModelParameters(&reader);

  // Check all defaults
  CHECK(params.plasticity == false);
  CHECK(params.plasticityPointwise == true);
  CHECK(params.useCellHomogenizedMaterial == true);
  CHECK(params.gravitationalAcceleration == doctest::Approx(9.81));
  CHECK(params.tv == doctest::Approx(0.1));
  CHECK(params.flux == NumericalFlux::Godunov);
  CHECK(params.fluxNearFault == NumericalFlux::Godunov);
}

} // namespace seissol::unit_test
