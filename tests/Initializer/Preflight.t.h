// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/Filesystem.h"
#include "Initializer/Preflight.h"
#include "TestHelper.h"

#include <array>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// an error ends the process that checks, so that the check runs in a child; Windows has no fork()
#ifndef _WIN32
#include <sys/wait.h>
#include <unistd.h>
#endif

namespace seissol::unit_test {

TEST_CASE("The check options are taken out of the command line" *
          doctest::test_suite("initializer")) {
  SUBCASE("Not given") {
    std::array<std::string, 2> words{"seissol", "parameters.par"};
    std::array<char*, 3> argv{words[0].data(), words[1].data(), nullptr};
    int argc = 2;
    CHECK_FALSE(initializer::extractCheckOptions(argc, argv.data()).has_value());
    CHECK(argc == 2);
  }
  SUBCASE("Without a report") {
    std::array<std::string, 3> words{"seissol", "--check", "parameters.par"};
    std::array<char*, 4> argv{words[0].data(), words[1].data(), words[2].data(), nullptr};
    int argc = 3;
    const auto options = initializer::extractCheckOptions(argc, argv.data());
    REQUIRE(options.has_value());
    CHECK(options->reportPath.empty());
    CHECK(options->depth == initializer::CheckDepth::Mesh);
    REQUIRE(argc == 2);
    CHECK(std::string(argv[1]) == "parameters.par");
  }
  SUBCASE("With a report") {
    std::array<std::string, 3> words{"seissol", "parameters.par", "--check=report.json"};
    std::array<char*, 4> argv{words[0].data(), words[1].data(), words[2].data(), nullptr};
    int argc = 3;
    const auto options = initializer::extractCheckOptions(argc, argv.data());
    REQUIRE(options.has_value());
    CHECK(options->reportPath == "report.json");
    CHECK(argc == 2);
  }
  SUBCASE("With a depth") {
    std::array<std::string, 4> words{
        "seissol", "--check-depth=files", "parameters.par", "--check=report.json"};
    std::array<char*, 5> argv{
        words[0].data(), words[1].data(), words[2].data(), words[3].data(), nullptr};
    int argc = 4;
    const auto options = initializer::extractCheckOptions(argc, argv.data());
    REQUIRE(options.has_value());
    CHECK(options->reportPath == "report.json");
    CHECK(options->depth == initializer::CheckDepth::Files);
    REQUIRE(argc == 2);
    CHECK(std::string(argv[1]) == "parameters.par");
  }
  SUBCASE("A depth alone asks for the check") {
    std::array<std::string, 3> words{"seissol", "--check-depth=mesh", "parameters.par"};
    std::array<char*, 4> argv{words[0].data(), words[1].data(), words[2].data(), nullptr};
    int argc = 3;
    const auto options = initializer::extractCheckOptions(argc, argv.data());
    REQUIRE(options.has_value());
    CHECK(options->depth == initializer::CheckDepth::Mesh);
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
  const int code = initializer::preflight(directory + "/parameters.par",
                                          {"abortcriteria.endtime=2.0"},
                                          std::nullopt,
                                          {.reportPath = reportPath});
  CHECK(code == 1);

  std::stringstream content;
  {
    // (closed again before the file is removed; Windows cannot remove an open file)
    const std::ifstream stream(reportPath);
    content << stream.rdbuf();
  }
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

#ifndef _WIN32
TEST_CASE("An error stops the check, and its report says which and where" *
          doctest::test_suite("initializer")) {
  const std::string directory = "preflight-error-test";
  std::filesystem::create_directories(directory);
  {
    std::ofstream file(directory + "/parameters.par");
    file << "&equations\nMaterialFileName = 'material.yaml'\n/\n"
         << "&MeshNml\nMeshFile = 'mesh.puml.h5'\nmeshgenerator = 'PUML'\n/\n"
         << "&Output\nOutputFile = 'output/data'\nReceiverOutput = 0\n"
         << "iOutputMask = 1 1 1 1 1 1 1 1 1\n/\n"
         << "&AbortCriteria\nEndTime = 1.0\n/\n";
  }
  const auto reportPath = directory + "/report.json";

  // the error ends the process that checks, so the check runs in a child; what is still buffered
  // would be written by both processes
  static_cast<void>(std::fflush(nullptr));
  const pid_t child = fork();
  if (child == 0) {
    static_cast<void>(initializer::preflight(directory + "/parameters.par",
                                             {"nosuchsection.field=1"},
                                             std::nullopt,
                                             {.reportPath = reportPath}));
    std::_Exit(0);
  }
  int status = 0;
  REQUIRE(waitpid(child, &status, 0) == child);
  CHECK(WIFEXITED(status));
  CHECK(WEXITSTATUS(status) == 1);

  const std::ifstream stream(reportPath);
  std::stringstream content;
  content << stream.rdbuf();
  const auto report = content.str();
  CHECK(report.find("\"status\": \"failed\"") != std::string::npos);
  CHECK(report.find("\"stoppedAt\": \"parameters\"") != std::string::npos);
  CHECK(report.find("names a parameter that SeisSol does not read") != std::string::npos);
  CHECK(report.find("Preflight.cpp") != std::string::npos);

  std::filesystem::remove_all(directory);
}
#endif

TEST_CASE("The check of a run finds the receivers that lie in no cell" *
          doctest::test_suite("initializer")) {
  const std::string directory = "preflight-mesh-test";
  std::filesystem::create_directories(directory);
  // the mesh is the unit cube, with absorbing faces all around
  {
    std::ofstream file(directory + "/parameters.par");
    file << "&equations\nMaterialFileName = 'material.yaml'\n/\n"
         << "&MeshNml\nMeshFile = '" << std::filesystem::absolute(tpath("Testing/mesh.h5")).string()
         << "'\nmeshgenerator = 'PUML'\n/\n"
         << "&Output\nOutputFile = 'output/data'\nReceiverOutput = 1\n"
         << "RFileName = 'receivers.dat'\npickdt = 0.01\n"
         << "iOutputMask = 1 1 1 1 1 1 1 1 1\n/\n"
         << "&AbortCriteria\nEndTime = 1.0\n/\n";
  }
  {
    std::ofstream file(directory + "/material.yaml");
    file << "!ConstantMap\nmap:\n  rho: 1\n  mu: 1\n  lambda: 1\n";
  }
  {
    // inside, on the boundary (which the tolerance of the receivers takes), and beyond it
    std::ofstream file(directory + "/receivers.dat");
    file << "0.5 0.5 0.5\n\n1.0 0.25 0.75\n0.5 0.5 1.25\n";
  }
  const auto reportPath = directory + "/report.json";

  SUBCASE("With the whole mesh") {
    const int code = initializer::preflight(
        directory + "/parameters.par", {}, std::nullopt, {.reportPath = reportPath});
    CHECK(code == 1);
    std::stringstream content;
    {
      const std::ifstream stream(reportPath);
      content << stream.rdbuf();
    }
    const auto report = content.str();
    CHECK(report.find("\"depth\": \"mesh\"") != std::string::npos);
    CHECK(report.find("The receiver 2 at (0.5, 0.5, 1.25)") != std::string::npos);
    CHECK(report.find("lies in no cell of the mesh: it lies outside the mesh, z = 1.25 not in "
                      "[0, 1]") != std::string::npos);
    CHECK(report.find("\"found\": 2") != std::string::npos);
    CHECK(report.find("\"5\": 24") != std::string::npos);
  }
  SUBCASE("With a face map that does not know a tag of the mesh") {
    {
      std::ofstream file(directory + "/facemap.yaml");
      file << "regular: \"0\"\nfreeSurface: \"1\"\n";
    }
    const int code = initializer::preflight(directory + "/parameters.par",
                                            {"meshnml.facemapfile=facemap.yaml"},
                                            std::nullopt,
                                            {.reportPath = reportPath});
    CHECK(code == 1);
    std::stringstream content;
    {
      const std::ifstream stream(reportPath);
      content << stream.rdbuf();
    }
    CHECK(content.str().find("The mesh has 24 faces with the boundary condition tag 5, which the "
                             "face map does not know") != std::string::npos);
  }
  SUBCASE("Without reading the whole mesh") {
    const int code =
        initializer::preflight(directory + "/parameters.par",
                               {},
                               std::nullopt,
                               {.reportPath = reportPath, .depth = initializer::CheckDepth::Files});
    CHECK(code == 0);
  }

  std::filesystem::remove_all(directory);
}

} // namespace seissol::unit_test
