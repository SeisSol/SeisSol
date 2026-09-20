// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Initializer/Parameters/OutputParameters.h"
#include "Initializer/Parameters/ParameterReader.h"

#include <yaml-cpp/yaml.h>

namespace seissol::unit_test {
using namespace seissol::initializer::parameters;

TEST_CASE("readOutputParameters: time series defaults to snapshots" *
          doctest::test_suite("initializer")) {
  const YAML::Node node = YAML::Load(R"(
    output:
      outputfile: 'out'
  )");
  ParameterReader reader(node, "", false);
  const auto params = readOutputParameters(&reader);

  CHECK(params.waveFieldParameters.timeSeries == TimeSeriesMode::Snapshot);
  CHECK(params.freeSurfaceParameters.timeSeries == TimeSeriesMode::Snapshot);
  CHECK(params.elementwiseParameters.timeSeries == TimeSeriesMode::Snapshot);
}

TEST_CASE("readOutputParameters: one setting reaches every mesh output" *
          doctest::test_suite("initializer")) {
  const YAML::Node node = YAML::Load(R"(
    output:
      outputfile: 'out'
      outputtimeseries: 'monolith'
  )");
  ParameterReader reader(node, "", false);
  const auto params = readOutputParameters(&reader);

  CHECK(params.waveFieldParameters.timeSeries == TimeSeriesMode::Monolith);
  CHECK(params.freeSurfaceParameters.timeSeries == TimeSeriesMode::Monolith);
  CHECK(params.elementwiseParameters.timeSeries == TimeSeriesMode::Monolith);
}

TEST_CASE("readOutputParameters: an output keeps its own setting" *
          doctest::test_suite("initializer")) {
  const YAML::Node node = YAML::Load(R"(
    output:
      outputfile: 'out'
      outputtimeseries: 'incremental'
      wavefieldtimeseries: 'monolith'
      surfacetimeseries: 'snapshot'
    elementwise:
      timeseries: 'monolith'
  )");
  ParameterReader reader(node, "", false);
  const auto params = readOutputParameters(&reader);

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
      outputtimeseries: 'incremental'
    elementwise:
      vtkorder: 2
  )");
  ParameterReader reader(node, "", false);
  const auto params = readOutputParameters(&reader);

  CHECK(params.elementwiseParameters.timeSeries == TimeSeriesMode::Incremental);
}

} // namespace seissol::unit_test
