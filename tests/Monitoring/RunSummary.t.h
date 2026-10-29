// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/Filesystem.h"
#include "IO/Instance/Point/Csv.h"
#include "Monitoring/RunSummary.h"
#include "Numerical/Statistics.h"
#include "Parallel/MPI.h"

#include <cstdlib>
#include <fstream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

namespace seissol::unit_test {

inline std::string readFile(const std::string& path) {
  const std::ifstream stream(path);
  std::stringstream content;
  content << stream.rdbuf();
  return content.str();
}

TEST_CASE("Run summary statistics" * doctest::test_suite("monitoring")) {
  const auto value = monitoring::toJson(statistics::Summary(std::vector<double>{1.0, 2.0, 6.0}));
  CHECK(value.find("min")->asDouble() == 1.0);
  CHECK(value.find("max")->asDouble() == 6.0);
  CHECK(value.find("mean")->asDouble() == 3.0);
  CHECK(value.find("median")->asDouble() == 2.0);
  CHECK(value.find("sum")->asDouble() == 9.0);
  CHECK(value.find("std") != nullptr);
}

TEST_CASE("Run summary timestamps" * doctest::test_suite("monitoring")) {
  const std::regex iso8601(R"(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z)");
  CHECK(std::regex_match(monitoring::utcTimestamp(), iso8601));
}

// NOLINTBEGIN(concurrency-mt-unsafe)
TEST_CASE("Run summary environment" * doctest::test_suite("monitoring")) {
  setenv("SEISSOL_SUMMARY_TEST_VARIABLE", "value", 1);
  setenv("OMPI_MCA_orte_summary_test", "internal", 1);
  setenv("UNRELATED_SUMMARY_TEST_VARIABLE", "value", 1);
  const auto environment = monitoring::relevantEnvironment();
  REQUIRE(environment.find("SEISSOL_SUMMARY_TEST_VARIABLE") != nullptr);
  CHECK(environment.find("SEISSOL_SUMMARY_TEST_VARIABLE")->asString() == "value");
  CHECK(environment.find("OMPI_MCA_orte_summary_test") == nullptr);
  CHECK(environment.find("UNRELATED_SUMMARY_TEST_VARIABLE") == nullptr);
  unsetenv("SEISSOL_SUMMARY_TEST_VARIABLE");
  unsetenv("OMPI_MCA_orte_summary_test");
  unsetenv("UNRELATED_SUMMARY_TEST_VARIABLE");
}
// NOLINTEND(concurrency-mt-unsafe)

TEST_CASE("Run summary files" * doctest::test_suite("monitoring")) {
  const std::string prefix = "runsummary-test";
  const std::string runId = "123e4567-e89b-42d3-a456-426614174000";
  const auto summaryPath = prefix + "-summary.json";
  const auto epochPath = prefix + "-epochs.csv";
  std::filesystem::remove(summaryPath);
  std::filesystem::remove(epochPath);

  monitoring::RunSummary summary;
  summary["invocation"]["parameterFile"] = "parameters.par";
  CHECK_FALSE(summary.started());
  // nothing is written before the start
  summary.setStage("init.mesh");
  CHECK_FALSE(std::filesystem::exists(summaryPath));

  summary.start(prefix, runId);
  REQUIRE(std::filesystem::exists(summaryPath));
  {
    const auto text = readFile(summaryPath);
    // the header comes first, what was there before after it
    CHECK(text.find("\"schema\": \"seissol-run-summary/1\"") < text.find("\"invocation\""));
    CHECK(text.find("\"runId\": \"" + runId + "\"") != std::string::npos);
    CHECK(text.find("\"status\": \"running\"") != std::string::npos);
    CHECK(text.find("\"parameterFile\": \"parameters.par\"") != std::string::npos);
  }
  CHECK(summary.content().find("parallel")->find("ranks")->asDouble() >= 1);
  CHECK(summary.content().find("build")->find("configurations")->size() >= 1);

  summary.setStage("simulation");
  CHECK(readFile(summaryPath).find("\"stage\": \"simulation\"") != std::string::npos);

  monitoring::EpochRecord epoch;
  epoch.timeFrom = 0.0;
  epoch.timeTo = 0.5;
  epoch.wall = statistics::Summary(std::vector<double>{1.0, 3.0});
  summary.addEpoch(epoch, 2.0);
  epoch.timeFrom = 0.5;
  epoch.timeTo = 1.0;
  summary.addEpoch(epoch, 2.0);
  CHECK(summary.content().find("simulation")->find("progress")->asDouble() == 0.5);
  CHECK(summary.content().find("simulation")->find("epochs")->asDouble() == 2);

  summary.finish();
  const auto text = readFile(summaryPath);
  CHECK(text.find("\"status\": \"completed\"") != std::string::npos);
  CHECK(text.find("\"stage\": \"finished\"") != std::string::npos);

  const auto table = io::instance::point::readCsv(epochPath);
  REQUIRE(table.rows.size() == 2);
  CHECK(table.rows[1][table.column("epoch")] == "2");
  CHECK(std::stod(table.rows[1][table.column("timeTo")]) == 1.0);
  CHECK(std::stod(table.rows[0][table.column("wallMax")]) == 3.0);
  CHECK(table.rows[0][table.column("runId")] == runId);

  std::filesystem::remove(summaryPath);
  std::filesystem::remove(epochPath);
}

TEST_CASE("Run summary rank table" * doctest::test_suite("monitoring")) {
  const std::string prefix = "runsummary-ranks-test";
  const std::string runId = "123e4567-e89b-42d3-a456-426614174000";
  const auto rankPath = prefix + "-ranks.csv";
  std::filesystem::remove(prefix + "-summary.json");
  std::filesystem::remove(rankPath);

  monitoring::RunSummary summary;
  summary.start(prefix, runId);
  summary.setRankValue("cells", 1.0);
  const auto statistics = summary.overRanks("cellUpdates", 2.0);
  CHECK(statistics.find("min")->asDouble() == 2.0);
  // recording a column again replaces its value, in its place
  summary.setRankValue("cells", 3.0);
  summary.writeRanks();

  const auto table = io::instance::point::readCsv(rankPath);
  REQUIRE(table.header.size() == 6);
  CHECK(table.header[0] == "rank");
  CHECK(table.header[1] == "host");
  CHECK(table.header[2] == "localRank");
  CHECK(table.header[3] == "cells");
  CHECK(table.header[4] == "cellUpdates");
  CHECK(table.header[5] == "runId");
  REQUIRE(table.rows.size() == static_cast<std::size_t>(seissol::Mpi::mpi.size()));
  CHECK(table.rows[0][table.column("rank")] == "0");
  CHECK(std::stod(table.rows[0][table.column("cells")]) == 3.0);
  CHECK(std::stod(table.rows[0][table.column("cellUpdates")]) == 2.0);
  CHECK(table.rows[0][table.column("runId")] == runId);

  std::filesystem::remove(prefix + "-summary.json");
  std::filesystem::remove(rankPath);
}

} // namespace seissol::unit_test
