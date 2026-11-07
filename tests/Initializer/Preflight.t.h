// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/Filesystem.h"
#include "Initializer/Preflight.h"

#include <array>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace seissol::unit_test {

TEST_CASE("The check option is taken out of the command line" *
          doctest::test_suite("initializer")) {
  SUBCASE("Not given") {
    std::array<std::string, 2> words{"seissol", "parameters.par"};
    std::array<char*, 3> argv{words[0].data(), words[1].data(), nullptr};
    int argc = 2;
    CHECK_FALSE(initializer::extractCheckOption(argc, argv.data()).has_value());
    CHECK(argc == 2);
  }
  SUBCASE("Without a report") {
    std::array<std::string, 3> words{"seissol", "--check", "parameters.par"};
    std::array<char*, 4> argv{words[0].data(), words[1].data(), words[2].data(), nullptr};
    int argc = 3;
    const auto report = initializer::extractCheckOption(argc, argv.data());
    REQUIRE(report.has_value());
    CHECK(report.value().empty());
    REQUIRE(argc == 2);
    CHECK(std::string(argv[1]) == "parameters.par");
  }
  SUBCASE("With a report") {
    std::array<std::string, 3> words{"seissol", "parameters.par", "--check=report.json"};
    std::array<char*, 4> argv{words[0].data(), words[1].data(), words[2].data(), nullptr};
    int argc = 3;
    const auto report = initializer::extractCheckOption(argc, argv.data());
    REQUIRE(report.has_value());
    CHECK(report.value() == "report.json");
    CHECK(argc == 2);
  }
}

TEST_CASE("The check of a run finds the files it misses" * doctest::test_suite("initializer")) {
  const std::string directory = "preflight-test";
  std::filesystem::create_directories(directory);
  {
    std::ofstream file(directory + "/parameters.par");
    file << "&equations\nMaterialFileName = 'material.yaml'\n/\n"
         << "&MeshNml\nMeshFile = 'mesh.puml.h5'\nmeshgenerator = 'PUML'\n/\n"
         << "&Output\nOutputFile = 'output/data'\nReceiverOutput = 0\n"
         << "iOutputMask = 1 1 1 1 1 1 1 1 1\n/\n"
         << "&AbortCriteria\nEndTime = 1.0\n/\n";
  }
  {
    std::ofstream file(directory + "/material.yaml");
    file << "!ConstantMap\nmap:\n  rho: [1\n";
  }
  const auto reportPath = directory + "/report.json";
  const int code = initializer::preflight(
      directory + "/parameters.par", {"abortcriteria.endtime=2.0"}, std::nullopt, reportPath);
  CHECK(code == 1);

  const std::ifstream stream(reportPath);
  std::stringstream content;
  content << stream.rdbuf();
  const auto report = content.str();
  CHECK(report.find("\"schema\": \"seissol-check/1\"") != std::string::npos);
  CHECK(report.find("\"status\": \"failed\"") != std::string::npos);
  CHECK(report.find("\"endTime\": 2") != std::string::npos);
  CHECK(report.find("\"abortcriteria.endtime\": \"2.0\"") != std::string::npos);
  // the material file is there, but does not parse; the mesh is not there
  CHECK(report.find("names does not parse") != std::string::npos);
  CHECK(report.find("mesh.puml.h5 that meshnml.meshfile names does not exist") !=
        std::string::npos);
  // the output directory does not exist yet, but can be created
  CHECK(report.find("\"writable\": true") != std::string::npos);

  std::filesystem::remove_all(directory);
}

} // namespace seissol::unit_test
