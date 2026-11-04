// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/ConfigRegistry.h"
#include "Initializer/Parameters/OutputParameters.h"
#include "Initializer/Parameters/ParameterReader.h"

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
  const auto params = readOutputParameters(&reader, defaultConfig());

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
  const auto params = readOutputParameters(&reader, defaultConfig());

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
  const auto params = readOutputParameters(&reader, defaultConfig());

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
  const auto params = readOutputParameters(&reader, defaultConfig());

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
  const auto params = readOutputParameters(&reader, defaultConfig());

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
  const auto params = readOutputParameters(&reader, defaultConfig());

  CHECK(params.elementwiseParameters.timeSeries == TimeSeriesMode::Incremental);
}

} // namespace seissol::unit_test
