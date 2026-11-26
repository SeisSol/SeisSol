// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Core/Common/ConfigRegistry.h"
#include "Core/IO/WriterParameters.h"
#include "Features/DynamicRupture/Output/FaultOutputParameters.h"
#include "Simulation/ConfigLayout.h"
#include "Simulation/Parameters/OutputParameters.h"
#include "Simulation/Parameters/ParameterReader.h"
#include "WavePropagation/Output/ReceiverParameters.h"

#include <cstddef>
#include <numeric>
#include <optional>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace seissol::unit_test {
using namespace seissol::initializer::parameters;

// the settings of the only output of the kind T
template <typename T>
const T& onlyOutputOf(const OutputParameters& params) {
  const auto found = params.outputsOf<T>();
  REQUIRE(found.size() == 1);
  return *found.front();
}

// the names of the outputs, in the order of the list
inline std::vector<std::string> outputNames(const OutputParameters& params) {
  std::vector<std::string> names;
  names.reserve(params.outputs.size());
  for (const auto& entry : params.outputs) {
    names.push_back(entry.name);
  }
  return names;
}

TEST_CASE("readOutputParameters: without a time series setting, the format of an output decides" *
          doctest::test_suite("initializer")) {
  const YAML::Node node = YAML::Load(R"(
    output:
      outputfile: 'out'
      receiveroutput: 0
      wavefieldoutput: 1
      ioutputmask: '1'
      surfaceoutput: 1
    dynamicrupture:
      outputpointtype: 4
  )");
  ParameterReader reader(node, "", false);
  const auto params = readOutputParameters(&reader, {defaultConfig()}, true);

  CHECK_FALSE(onlyOutputOf<WaveFieldOutputParameters>(params).timeSeries.has_value());
  CHECK_FALSE(onlyOutputOf<FreeSurfaceOutputParameters>(params).timeSeries.has_value());
  CHECK_FALSE(onlyOutputOf<ElementwiseFaultParameters>(params).timeSeries.has_value());
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
  const auto params = readOutputParameters(&reader, {defaultConfig()}, false);

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
  const auto params = readOutputParameters(&reader, {defaultConfig()}, false);

  CHECK(params.epochInterval == 0.25);
}

TEST_CASE("readOutputParameters: one setting reaches every mesh output" *
          doctest::test_suite("initializer")) {
  const YAML::Node node = YAML::Load(R"(
    output:
      outputfile: 'out'
      receiveroutput: 0
      wavefieldoutput: 1
      ioutputmask: '1'
      surfaceoutput: 1
      outputtimeseries: 'monolith'
    dynamicrupture:
      outputpointtype: 4
  )");
  ParameterReader reader(node, "", false);
  const auto params = readOutputParameters(&reader, {defaultConfig()}, true);

  CHECK(onlyOutputOf<WaveFieldOutputParameters>(params).timeSeries == TimeSeriesMode::Monolith);
  CHECK(onlyOutputOf<FreeSurfaceOutputParameters>(params).timeSeries == TimeSeriesMode::Monolith);
  CHECK(onlyOutputOf<ElementwiseFaultParameters>(params).timeSeries == TimeSeriesMode::Monolith);
}

TEST_CASE("readOutputParameters: an output keeps its own setting" *
          doctest::test_suite("initializer")) {
  const YAML::Node node = YAML::Load(R"(
    output:
      outputfile: 'out'
      receiveroutput: 0
      wavefieldoutput: 1
      ioutputmask: '1'
      surfaceoutput: 1
      outputtimeseries: 'incremental'
      wavefieldtimeseries: 'monolith'
      surfacetimeseries: 'snapshot'
    elementwise:
      timeseries: 'monolith'
    dynamicrupture:
      outputpointtype: 4
  )");
  ParameterReader reader(node, "", false);
  const auto params = readOutputParameters(&reader, {defaultConfig()}, true);

  CHECK(onlyOutputOf<WaveFieldOutputParameters>(params).timeSeries == TimeSeriesMode::Monolith);
  CHECK(onlyOutputOf<FreeSurfaceOutputParameters>(params).timeSeries == TimeSeriesMode::Snapshot);
  CHECK(onlyOutputOf<ElementwiseFaultParameters>(params).timeSeries == TimeSeriesMode::Monolith);
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
    dynamicrupture:
      outputpointtype: 4
  )");
  ParameterReader reader(node, "", false);
  const auto params = readOutputParameters(&reader, {defaultConfig()}, true);

  CHECK(onlyOutputOf<ElementwiseFaultParameters>(params).timeSeries == TimeSeriesMode::Incremental);
}

TEST_CASE("readOutputParameters: every configuration of the run masks its own quantities" *
          doctest::test_suite("initializer")) {
  const YAML::Node node = YAML::Load(R"(
    output:
      outputfile: 'out'
      receiveroutput: 0
      wavefieldoutput: 1
      ioutputmask: '1 1 1 1 1 1 1 1 1'
      integrationmask: '0 1'
      surfaceoutput: 1
  )");
  ParameterReader reader(node, "", false);
  std::vector<ConfigId> configs(builtConfigCount());
  std::iota(configs.begin(), configs.end(), 0);
  const auto params = readOutputParameters(&reader, configs, false);

  const auto& waveField = onlyOutputOf<WaveFieldOutputParameters>(params);
  const auto& surface = onlyOutputOf<FreeSurfaceOutputParameters>(params);
  for (const auto config : configs) {
    const auto& layout = configLayout(config);
    const auto& outputMask = waveField.outputMask.at(config);
    const auto& integrationMask = waveField.integrationMask.at(config);
    const auto& surfaceMask = surface.outputMask.at(config);
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
    const auto params = readOutputParameters(&reader, {defaultConfig()}, false);

    CHECK(params.outputsOf<ReceiverOutputParameters>().empty());
    CHECK_FALSE(params.checkpoint.has_value());
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
    const auto params = readOutputParameters(&reader, {defaultConfig()}, false);

    CHECK(onlyOutputOf<ReceiverOutputParameters>(params).fileName == "receivers.dat");
    REQUIRE(params.checkpoint.has_value());
    CHECK(params.checkpoint.value().interval == 0.5);
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
    const auto params = readOutputParameters(&reader, {defaultConfig()}, false);

    CHECK(params.outputsOf<ReceiverOutputParameters>().empty());
    CHECK_FALSE(params.checkpoint.has_value());
  }
}

TEST_CASE("readOutputParameters: the outputs of a run, each named as the files it writes" *
          doctest::test_suite("initializer")) {
  SUBCASE("none") {
    const YAML::Node node = YAML::Load(R"(
      output:
        outputfile: 'out'
        wavefieldoutput: 0
        ioutputmask: '1'
    )");
    ParameterReader reader(node, "", false);
    const auto params = readOutputParameters(&reader, {defaultConfig()}, false);

    CHECK(params.outputs.empty());
    CHECK_FALSE(params.energy.has_value());
  }

  SUBCASE("all of them") {
    const YAML::Node node = YAML::Load(R"(
      output:
        outputfile: 'out'
        wavefieldoutput: 1
        ioutputmask: '1'
        surfaceoutput: 1
        rfilename: 'receivers.dat'
        energyoutput: 1
      dynamicrupture:
        outputpointtype: 5
    )");
    ParameterReader reader(node, "", false);
    const auto params = readOutputParameters(&reader, {defaultConfig()}, true);

    CHECK(outputNames(params) ==
          std::vector<std::string>{"wavefield", "surface", "receivers", "fault", "faultreceivers"});
    CHECK(params.energy.has_value());
  }

  SUBCASE("those of the fault only with a fault") {
    const YAML::Node node = YAML::Load(R"(
      output:
        outputfile: 'out'
        ioutputmask: '1'
      dynamicrupture:
        outputpointtype: 5
    )");
    ParameterReader reader(node, "", false);
    const auto params = readOutputParameters(&reader, {defaultConfig()}, false);

    CHECK(params.outputs.empty());
  }
}

TEST_CASE("readOutputParameters: outputs in blocks of their own" *
          doctest::test_suite("initializer")) {
  SUBCASE("named after their kind and number, ordered by both") {
    const YAML::Node node = YAML::Load(R"(
      output:
        outputfile: 'out'
        ioutputmask: '1'
      outputsurface3:
        interval: 1.0
      outputwavefield2:
        mask: '1 1 1'
        interval: 0.5
        vtkorder: 2
      outputwavefield:
        mask: '1'
      outputwavefieldx:
        mask: '1'
    )");
    ParameterReader reader(node, "", false);
    const auto params = readOutputParameters(&reader, {defaultConfig()}, false);

    CHECK(outputNames(params) == std::vector<std::string>{"wavefield", "wavefield2", "surface3"});
    const auto waveFields = params.outputsOf<WaveFieldOutputParameters>();
    REQUIRE(waveFields.size() == 2);
    CHECK(waveFields[1]->interval == 0.5);
    CHECK(waveFields[1]->vtkorder == 2);
    CHECK(waveFields[1]->outputMask.at(defaultConfig()).at(2));
    CHECK_FALSE(waveFields[0]->outputMask.at(defaultConfig()).at(2));
    CHECK(onlyOutputOf<FreeSurfaceOutputParameters>(params).interval == 1.0);
  }

  SUBCASE("next to the one given the way it has always been, under a name of its own") {
    const YAML::Node node = YAML::Load(R"(
      output:
        outputfile: 'out'
        wavefieldoutput: 1
        ioutputmask: '1'
      outputwavefield:
        name: 'near-fault'
        mask: '1 1 1'
    )");
    ParameterReader reader(node, "", false);
    const auto params = readOutputParameters(&reader, {defaultConfig()}, false);

    CHECK(outputNames(params) == std::vector<std::string>{"wavefield", "near-fault"});
  }

  SUBCASE("receivers with their file and their own sampling") {
    const YAML::Node node = YAML::Load(R"(
      output:
        outputfile: 'out'
        ioutputmask: '1'
      outputreceivers:
        file: 'receivers.dat'
        samplinginterval: 0.01
        computerotation: 1
    )");
    ParameterReader reader(node, "", false);
    const auto params = readOutputParameters(&reader, {defaultConfig()}, false);

    const auto& receivers = onlyOutputOf<ReceiverOutputParameters>(params);
    CHECK(receivers.fileName == "receivers.dat");
    CHECK(receivers.samplingInterval == 0.01);
    CHECK(receivers.computeRotation);
    CHECK(outputNames(params) == std::vector<std::string>{"receivers"});
  }

  SUBCASE("several outputs of receivers") {
    const YAML::Node node = YAML::Load(R"(
      output:
        outputfile: 'out'
        ioutputmask: '1'
        rfilename: 'receivers.dat'
      outputreceivers2:
        file: 'more.dat'
        samplinginterval: 0.1
    )");
    ParameterReader reader(node, "", false);
    const auto params = readOutputParameters(&reader, {defaultConfig()}, false);

    CHECK(outputNames(params) == std::vector<std::string>{"receivers", "receivers2"});
    const auto receivers = params.outputsOf<ReceiverOutputParameters>();
    REQUIRE(receivers.size() == 2);
    CHECK(receivers[1]->fileName == "more.dat");
    CHECK(receivers[1]->samplingInterval == 0.1);
  }

  SUBCASE("those of the fault only with a fault") {
    const YAML::Node node = YAML::Load(R"(
      output:
        outputfile: 'out'
        ioutputmask: '1'
      outputfault:
        interval: 0.1
        mask: '1 1'
      outputfaultreceivers:
        file: 'points.dat'
        samplingsteps: 2
    )");
    SUBCASE("without one") {
      ParameterReader reader(node, "", false);
      const auto params = readOutputParameters(&reader, {defaultConfig()}, false);
      CHECK(params.outputs.empty());
    }
    SUBCASE("with one") {
      ParameterReader reader(node, "", false);
      const auto params = readOutputParameters(&reader, {defaultConfig()}, true);
      CHECK(outputNames(params) == std::vector<std::string>{"faultreceivers", "fault"});
      const auto& fault = onlyOutputOf<ElementwiseFaultParameters>(params);
      CHECK(fault.printTimeIntervalSec == 0.1);
      CHECK(fault.outputMask[1]);
      CHECK_FALSE(fault.outputMask[2]);
      const auto& points = onlyOutputOf<PickpointParameters>(params);
      CHECK(points.pickpointFileName == std::optional<std::string>("points.dat"));
      CHECK(points.printTimeInterval == 2);
    }
  }
}

TEST_CASE("readOutputParameters: the directory and the name of the outputs" *
          doctest::test_suite("initializer")) {
  SUBCASE("from outdir and outname") {
    const YAML::Node node = YAML::Load(R"(
      output:
        outdir: 'runs/strong-16'
        outname: 'data'
        receiveroutput: 0
        ioutputmask: '1'
    )");
    ParameterReader reader(node, "", false);
    const auto params = readOutputParameters(&reader, {defaultConfig()}, false);

    CHECK(params.outputPath == "runs/strong-16");
    CHECK(params.outputName == "data");
    CHECK_FALSE(params.outputDirPerSim);
  }

  SUBCASE("from the prefix in outputfile") {
    const YAML::Node node = YAML::Load(R"(
      output:
        outputfile: 'runs/strong-16/data'
        receiveroutput: 0
        ioutputmask: '1'
    )");
    ParameterReader reader(node, "", false);
    const auto params = readOutputParameters(&reader, {defaultConfig()}, false);

    CHECK(params.outputPath == "runs/strong-16");
    CHECK(params.outputName == "data");
  }

  SUBCASE("in the working directory, from a prefix without a directory") {
    const YAML::Node node = YAML::Load(R"(
      output:
        outputfile: 'data'
        receiveroutput: 0
        ioutputmask: '1'
    )");
    ParameterReader reader(node, "", false);
    const auto params = readOutputParameters(&reader, {defaultConfig()}, false);

    CHECK(params.outputPath.empty());
    CHECK(params.outputName == "data");
  }

  SUBCASE("outdir in place of the directory of outputfile") {
    const YAML::Node node = YAML::Load(R"(
      output:
        outputfile: 'runs/strong-16/data'
        outdir: 'runs/strong-32'
        receiveroutput: 0
        ioutputmask: '1'
    )");
    ParameterReader reader(node, "", false);
    const auto params = readOutputParameters(&reader, {defaultConfig()}, false);

    CHECK(params.outputPath == "runs/strong-32");
    CHECK(params.outputName == "data");
  }

  SUBCASE("outname in place of the name of outputfile") {
    const YAML::Node node = YAML::Load(R"(
      output:
        outputfile: 'runs/strong-16/data'
        outname: 'other'
        receiveroutput: 0
        ioutputmask: '1'
    )");
    ParameterReader reader(node, "", false);
    const auto params = readOutputParameters(&reader, {defaultConfig()}, false);

    CHECK(params.outputPath == "runs/strong-16");
    CHECK(params.outputName == "other");
  }

  SUBCASE("outname alone, in the working directory") {
    const YAML::Node node = YAML::Load(R"(
      output:
        outname: 'data'
        receiveroutput: 0
        ioutputmask: '1'
    )");
    ParameterReader reader(node, "", false);
    const auto params = readOutputParameters(&reader, {defaultConfig()}, false);

    CHECK(params.outputPath.empty());
    CHECK(params.outputName == "data");
  }

  SUBCASE("outdir and outname before outputfile") {
    const YAML::Node node = YAML::Load(R"(
      output:
        outputfile: 'elsewhere/other'
        outdir: 'runs'
        outname: 'data'
        outsimdir: 1
        receiveroutput: 0
        ioutputmask: '1'
    )");
    ParameterReader reader(node, "", false);
    const auto params = readOutputParameters(&reader, {defaultConfig()}, false);

    CHECK(params.outputPath == "runs");
    CHECK(params.outputName == "data");
    CHECK(params.outputDirPerSim);
  }
}

} // namespace seissol::unit_test
