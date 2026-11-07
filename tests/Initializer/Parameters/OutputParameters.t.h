// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/ConfigLayout.h"
#include "Common/ConfigRegistry.h"
#include "Initializer/Parameters/OutputParameters.h"
#include "Initializer/Parameters/ParameterReader.h"

#include <cstddef>
#include <numeric>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace seissol::unit_test {
using namespace seissol::initializer::parameters;

TEST_CASE("readOutputParameters: without a time series setting, the format of an output decides" *
          doctest::test_suite("initializer")) {
  const YAML::Node node = YAML::Load(R"(
    output:
      outputfile: 'out'
      receiveroutput: 0
      ioutputmask: '1'
  )");
  ParameterReader reader(node, "", false);
  const auto params = readOutputParameters(&reader, {defaultConfig()});

  CHECK_FALSE(params.waveFieldParameters.timeSeries.has_value());
  CHECK_FALSE(params.freeSurfaceParameters.timeSeries.has_value());
  CHECK_FALSE(params.elementwiseParameters.timeSeries.has_value());
}

TEST_CASE("readOutputParameters: epochs end at the outputs only by default" *
          doctest::test_suite("initializer")) {
  const YAML::Node node = YAML::Load(R"(
    output:
      outputfile: 'out'
      receiveroutput: 0
      ioutputmask: '1'
  )");
  ParameterReader reader(node, "", false);
  const auto params = readOutputParameters(&reader, {defaultConfig()});

  CHECK(params.epochInterval == 0.0);
}

TEST_CASE("readOutputParameters: an epoch interval of its own" *
          doctest::test_suite("initializer")) {
  const YAML::Node node = YAML::Load(R"(
    output:
      outputfile: 'out'
      receiveroutput: 0
      ioutputmask: '1'
      epochinterval: 0.25
  )");
  ParameterReader reader(node, "", false);
  const auto params = readOutputParameters(&reader, {defaultConfig()});

  CHECK(params.epochInterval == 0.25);
}

TEST_CASE("readOutputParameters: one setting reaches every mesh output" *
          doctest::test_suite("initializer")) {
  const YAML::Node node = YAML::Load(R"(
    output:
      outputfile: 'out'
      receiveroutput: 0
      ioutputmask: '1'
      outputtimeseries: 'monolith'
  )");
  ParameterReader reader(node, "", false);
  const auto params = readOutputParameters(&reader, {defaultConfig()});

  CHECK(params.waveFieldParameters.timeSeries == TimeSeriesMode::Monolith);
  CHECK(params.freeSurfaceParameters.timeSeries == TimeSeriesMode::Monolith);
  CHECK(params.elementwiseParameters.timeSeries == TimeSeriesMode::Monolith);
}

TEST_CASE("readOutputParameters: an output keeps its own setting" *
          doctest::test_suite("initializer")) {
  const YAML::Node node = YAML::Load(R"(
    output:
      outputfile: 'out'
      receiveroutput: 0
      ioutputmask: '1'
      outputtimeseries: 'incremental'
      wavefieldtimeseries: 'monolith'
      surfacetimeseries: 'snapshot'
    elementwise:
      timeseries: 'monolith'
  )");
  ParameterReader reader(node, "", false);
  const auto params = readOutputParameters(&reader, {defaultConfig()});

  CHECK(params.waveFieldParameters.timeSeries == TimeSeriesMode::Monolith);
  CHECK(params.freeSurfaceParameters.timeSeries == TimeSeriesMode::Snapshot);
  CHECK(params.elementwiseParameters.timeSeries == TimeSeriesMode::Monolith);
}

TEST_CASE("readOutputParameters: the fault output follows the shared setting" *
          doctest::test_suite("initializer")) {
  // The fault output reads its override from a section of its own, so the fallback has to travel
  // there as well as to the two that sit next to the shared field.
  const YAML::Node node = YAML::Load(R"(
    output:
      outputfile: 'out'
      receiveroutput: 0
      ioutputmask: '1'
      outputtimeseries: 'incremental'
    elementwise:
      vtkorder: 2
  )");
  ParameterReader reader(node, "", false);
  const auto params = readOutputParameters(&reader, {defaultConfig()});

  CHECK(params.elementwiseParameters.timeSeries == TimeSeriesMode::Incremental);
}

TEST_CASE("readOutputParameters: the output regions come from a file of their own" *
          doctest::test_suite("initializer")) {
  const YAML::Node none = YAML::Load(R"(
    output:
      outputfile: 'out'
      receiveroutput: 0
      ioutputmask: '1'
  )");
  ParameterReader noneReader(none, "", false);
  CHECK(readOutputParameters(&noneReader, {defaultConfig()}).regionFileName.empty());

  const YAML::Node node = YAML::Load(R"(
    output:
      outputfile: 'out'
      receiveroutput: 0
      ioutputmask: '1'
      outputregionfilename: 'sderiv:region.sderiv'
  )");
  ParameterReader reader(node, "", false);
  // the kind of the script stays in front of the path
  CHECK(readOutputParameters(&reader, {defaultConfig()}).regionFileName == "sderiv:region.sderiv");
}

TEST_CASE("readOutputParameters: every configuration of the run masks its own quantities" *
          doctest::test_suite("initializer")) {
  const YAML::Node node = YAML::Load(R"(
    output:
      outputfile: 'out'
      receiveroutput: 0
      ioutputmask: '1 1 1 1 1 1 1 1 1'
      integrationmask: '0 1'
  )");
  ParameterReader reader(node, "", false);
  std::vector<ConfigId> configs(builtConfigCount());
  std::iota(configs.begin(), configs.end(), 0);
  const auto params = readOutputParameters(&reader, configs);

  for (const auto config : configs) {
    const auto& layout = configLayout(config);
    const auto& outputMask = params.waveFieldParameters.outputMask.at(config);
    const auto& integrationMask = params.waveFieldParameters.integrationMask.at(config);
    const auto& surfaceMask = params.freeSurfaceParameters.outputMask.at(config);
    REQUIRE(outputMask.size() == layout.numQuantities);
    REQUIRE(integrationMask.size() == layout.numQuantities);
    REQUIRE(surfaceMask.size() == layout.numQuantities);
    for (std::size_t quantity = 0; quantity < layout.numQuantities; ++quantity) {
      CHECK(outputMask[quantity] == (quantity < 9));
      CHECK(integrationMask[quantity] == (quantity == 1));
      // the free surface writes the velocities by default
      CHECK(surfaceMask[quantity] ==
            (quantity >= layout.velocityOffset && quantity < layout.velocityOffset + 3));
    }
  }
}

TEST_CASE("readOutputParameters: the receivers and the checkpoints are on where their parameters "
          "say what to write" *
          doctest::test_suite("initializer")) {
  SUBCASE("neither of them without their parameters") {
    const YAML::Node node = YAML::Load(R"(
      output:
        outputfile: 'out'
        ioutputmask: '1'
    )");
    ParameterReader reader(node, "", false);
    const auto params = readOutputParameters(&reader, {defaultConfig()});

    CHECK_FALSE(params.receiverParameters.enabled);
    CHECK_FALSE(params.checkpointParameters.enabled);
  }

  SUBCASE("the receivers with their points, the checkpoints with their interval") {
    const YAML::Node node = YAML::Load(R"(
      output:
        outputfile: 'out'
        ioutputmask: '1'
        rfilename: 'receivers.dat'
        checkpointinterval: 0.5
    )");
    ParameterReader reader(node, "", false);
    const auto params = readOutputParameters(&reader, {defaultConfig()});

    CHECK(params.receiverParameters.enabled);
    CHECK(params.checkpointParameters.enabled);
    CHECK(params.checkpointParameters.interval == 0.5);
  }

  SUBCASE("both of them switched off in spite of their parameters") {
    const YAML::Node node = YAML::Load(R"(
      output:
        outputfile: 'out'
        ioutputmask: '1'
        receiveroutput: 0
        rfilename: 'receivers.dat'
        checkpoint: 0
        checkpointinterval: 0.5
    )");
    ParameterReader reader(node, "", false);
    const auto params = readOutputParameters(&reader, {defaultConfig()});

    CHECK_FALSE(params.receiverParameters.enabled);
    CHECK_FALSE(params.checkpointParameters.enabled);
  }
}

} // namespace seissol::unit_test
