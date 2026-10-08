// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/ConfigRegistry.h"
#include "Common/ConfigValue.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Initializer/Parameters/ParameterReader.h"
#include "TestHelper.h"

#include <cstddef>
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
// readConfig
// ---------------------------------------------------------------------------

TEST_CASE("readConfig takes the first configuration if none is named" *
          doctest::test_suite("initializer")) {
  const YAML::Node node = YAML::Load(R"(
    equations:
      materialfilename: mat.yaml
  )");
  ParameterReader reader(node, "", false);
  CHECK(readConfig(&reader) == defaultConfig());
}

TEST_CASE("readConfig finds every configuration built by its name" *
          doctest::test_suite("initializer")) {
  for (std::size_t id = 0; id < builtConfigCount(); ++id) {
    const auto name = configName(configValue(static_cast<ConfigId>(id)));
    CAPTURE(name);
    const YAML::Node node = YAML::Load("equations:\n  configuration: ' " + name + " '\n");
    ParameterReader reader(node, "", false);
    CHECK(readConfig(&reader) == id);
  }
}

TEST_CASE("readConfig takes the configuration of SEISSOL_CONFIGURATION if the file names none" *
          doctest::test_suite("initializer")) {
  const auto last = static_cast<ConfigId>(builtConfigCount() - 1);
  const ScopedEnvironment environment("SEISSOL_CONFIGURATION", configName(configValue(last)));

  SUBCASE("the parameter file names no configuration") {
    const YAML::Node node = YAML::Load("equations:\n  materialfilename: mat.yaml\n");
    ParameterReader reader(node, "", false);
    CHECK(readConfig(&reader) == last);
  }

  SUBCASE("the configuration of the parameter file applies") {
    const YAML::Node node =
        YAML::Load("equations:\n  configuration: " + configName(configValue(0)) + "\n");
    ParameterReader reader(node, "", false);
    CHECK(readConfig(&reader) == 0);
  }
}

// ---------------------------------------------------------------------------
// readGroupConfigs
// ---------------------------------------------------------------------------

TEST_CASE("readGroupConfigs gives no group a configuration of its own by default" *
          doctest::test_suite("initializer")) {
  const YAML::Node node = YAML::Load(R"(
    equations:
      materialfilename: mat.yaml
  )");
  ParameterReader reader(node, "", false);
  CHECK(readGroupConfigs(&reader, defaultConfig()).empty());
}

TEST_CASE("readGroupConfigs reads the groups of each configuration" *
          doctest::test_suite("initializer")) {
  // the last configuration that may share a run with the default one: all configurations of a run
  // fuse the same number of simulations
  auto last = defaultConfig();
  for (std::size_t id = 0; id < builtConfigCount(); ++id) {
    if (configValue(static_cast<ConfigId>(id)).numSimulations ==
        configValue(defaultConfig()).numSimulations) {
      last = static_cast<ConfigId>(id);
    }
  }
  const auto name = configName(configValue(last));
  const auto defaultName = configName(configValue(defaultConfig()));
  const YAML::Node node =
      YAML::Load("equations:\n  configmap: ' 1, 2 : " + name + " ; 7:" + defaultName + ";'\n");
  ParameterReader reader(node, "", false);
  const auto groupConfigs = readGroupConfigs(&reader, defaultConfig());
  CHECK(groupConfigs.size() == 3);
  CHECK(groupConfigs.at(1) == last);
  CHECK(groupConfigs.at(2) == last);
  CHECK(groupConfigs.at(7) == defaultConfig());
}

// ---------------------------------------------------------------------------
// ModelParameters::configOfGroup, ModelParameters::configs
// ---------------------------------------------------------------------------

TEST_CASE("The groups without a configuration of their own take the one of the run" *
          doctest::test_suite("initializer")) {
  ModelParameters params{};
  params.config = defaultConfig();
  params.groupConfigs = {{3, 5}, {1, 4}, {2, 5}};
  CHECK(params.configOfGroup(0) == defaultConfig());
  CHECK(params.configOfGroup(1) == 4);
  CHECK(params.configOfGroup(3) == 5);
  const auto configs = params.configs();
  REQUIRE(configs.size() == 3);
  CHECK(configs[0] == defaultConfig());
  CHECK(configs[1] == 4);
  CHECK(configs[2] == 5);
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
      freqcentral: 1
      freqratio: 1
  )");
  ParameterReader reader(node, "", false);
  auto params = readModelParameters(&reader, defaultConfig(), {});

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
      freqcentral: 1
      freqratio: 1
  )");
  ParameterReader reader(node, "", false);
  auto params = readModelParameters(&reader, defaultConfig(), {});

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
