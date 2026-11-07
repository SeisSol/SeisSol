// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Initializer/Parameters/ParameterOverrides.h"
#include "Initializer/Parameters/ParameterReader.h"

#include <fty/fty.hpp>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace seissol::unit_test {
using namespace seissol::initializer::parameters;

//! @brief The parameter file @p lines, read as SeisSol reads one.
inline YAML::Node readParameters(const fty::StringsT& lines) {
  fty::Converter<fty::AsLowercase> converter;
  return converter.convertToYAML(lines);
}

inline ParameterOverride makeOverride(const std::string& text) {
  const auto parsed = parseOverride(text);
  REQUIRE(parsed.has_value());
  return parsed.value();
}

TEST_CASE("Parsing parameter overrides" * doctest::test_suite("initializer")) {
  SUBCASE("section and field are taken in lower case, the value as it is") {
    const auto parsed = makeOverride("Output.OutputFile=Runs/A=1/data");
    CHECK(parsed.path == std::vector<std::string>{"output", "outputfile"});
    CHECK(parsed.value == "Runs/A=1/data");
    CHECK(parsed.name() == "output.outputfile");
  }

  SUBCASE("deeper paths") {
    CHECK(makeOverride("a.b.c=1").path == std::vector<std::string>{"a", "b", "c"});
  }

  SUBCASE("malformed overrides") {
    CHECK_FALSE(parseOverride("output.outputfile").has_value());
    CHECK_FALSE(parseOverride("outputfile=data").has_value());
    CHECK_FALSE(parseOverride(".outputfile=data").has_value());
    CHECK_FALSE(parseOverride("output.=data").has_value());
    CHECK_FALSE(parseOverride("output..file=data").has_value());
    CHECK_FALSE(parseOverride("out-put.file=data").has_value());
    CHECK_FALSE(parseOverride("=data").has_value());
  }
}

TEST_CASE("Applying parameter overrides" * doctest::test_suite("initializer")) {
  auto parameters = readParameters({"&output",
                                    "OutputFile = 'output/data'",
                                    "WavefieldOutput = 1",
                                    "/",
                                    "&AbortCriteria",
                                    "EndTime = 15.0",
                                    "/"});

  SUBCASE("values are read as in the parameter file") {
    applyOverrides(parameters,
                   {makeOverride("output.outputfile='runs/a/data'"),
                    makeOverride("abortcriteria.endtime=2.5d-1")});
    CHECK(parameters["output"]["outputfile"].as<std::string>() == "runs/a/data");
    CHECK(parameters["abortcriteria"]["endtime"].as<double>() == doctest::Approx(0.25));
    // the other fields stay
    CHECK(parameters["output"]["wavefieldoutput"].as<int>() == 1);
  }

  SUBCASE("a later override wins") {
    applyOverrides(
        parameters,
        {makeOverride("output.wavefieldoutput=0"), makeOverride("output.wavefieldoutput=1")});
    CHECK(parameters["output"]["wavefieldoutput"].as<int>() == 1);
  }

  SUBCASE("a section the file does not have is added") {
    applyOverrides(parameters, {makeOverride("discretization.clusteredlts=2")});
    CHECK(parameters["discretization"]["clusteredlts"].as<int>() == 2);
    CHECK(parameters["output"]["outputfile"].as<std::string>() == "output/data");
  }
}

TEST_CASE("Overrides of parameters that are not read" * doctest::test_suite("initializer")) {
  const auto parameters = readParameters(
      {"&output", "OutputFile = 'output/data'", "/", "&AbortCriteria", "EndTime = 15.0", "/"});
  const std::vector<ParameterOverride> overrides{makeOverride("output.outputfile=runs/a"),
                                                 makeOverride("output.outputfle=runs/a"),
                                                 makeOverride("abortcriteria.endtime=1"),
                                                 makeOverride("abortcriteria.maxtime=1"),
                                                 makeOverride("nosuchsection.field=1")};
  applyOverrides(parameters, overrides);

  ParameterReader reader(parameters, "parameters.par", false);
  auto* output = reader.readSubNode("output");
  CHECK(output->readWithDefault<std::string>("outputfile", "") == "runs/a");
  auto* abortCriteria = reader.readSubNode("abortcriteria");
  CHECK(abortCriteria->readWithDefault("endtime", 0.0) == doctest::Approx(1.0));
  // a field marked unused has been looked at
  abortCriteria->markUnused({"maxtime"});

  const auto unknown = unknownOverrides(reader, overrides);
  REQUIRE(unknown.size() == 2);
  CHECK(unknown[0].name() == "output.outputfle");
  CHECK(unknown[1].name() == "nosuchsection.field");
}

} // namespace seissol::unit_test
