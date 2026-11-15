// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_MONITORING_FAILUREREPORTMPI_T_H_
#define SEISSOL_TESTS_MONITORING_FAILUREREPORTMPI_T_H_

#include <doctest.h>

#include "Monitoring/FailureReport.h"
#include "Monitoring/RunSummary.h"
#include "Parallel/MPI.h"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <utils/logger.h>

namespace seissol::unit_test {

TEST_CASE("A failure on a rank that does not write the summary leaves its own record" *
          doctest::test_suite("monitoring")) {
  const int rank = seissol::Mpi::mpi.rank();
  const std::string prefix = "monitoringmpi-failure-test";
  const auto summaryPath = prefix + "-summary.json";

  // start the summary on every rank (collective), then fail on the ranks that do not write it
  monitoring::RunSummary summary;
  summary.start(prefix, "123e4567-e89b-42d3-a456-426614174000");
  summary.setStage("simulation");

  // a request to stop leaves no record there, as it reaches every rank
  if (rank != 0) {
    monitoring::writeTerminationReport(R"({"signal":"SIGTERM"})");
  }
  seissol::Mpi::barrier(seissol::Mpi::mpi.comm());
  if (rank == 0) {
    for (int other = 1; other < seissol::Mpi::mpi.size(); ++other) {
      CHECK_FALSE(
          std::filesystem::exists(prefix + "-rank" + std::to_string(other) + "-failure.json"));
    }
  }
  seissol::Mpi::barrier(seissol::Mpi::mpi.comm());

  // whereas a failure does
  if (rank != 0) {
    utils::LogRecord record;
    record.level = utils::LogLevel::Error;
    record.text = "A failure away from rank 0. ";
    record.location = {.file = "Cluster.cpp", .function = "compute", .line = 7};
    record.time = {.seconds = 1700000000, .nanoseconds = 0};
    monitoring::writeFailureReport(monitoring::failureOf(record));
  }

  seissol::Mpi::barrier(seissol::Mpi::mpi.comm());

  // rank 0 checks, so that the result is reported (only rank 0 lists the tests)
  if (rank == 0) {
    CHECK(std::filesystem::exists(summaryPath));
    // the display rank writes the summary, not a record of its own
    CHECK_FALSE(std::filesystem::exists(prefix + "-rank0-failure.json"));
    for (int other = 1; other < seissol::Mpi::mpi.size(); ++other) {
      const auto otherPath = prefix + "-rank" + std::to_string(other) + "-failure.json";
      REQUIRE(std::filesystem::exists(otherPath));
      const std::ifstream stream(otherPath);
      std::stringstream content;
      content << stream.rdbuf();
      const auto text = content.str();
      CHECK(text.find("\"schema\": \"seissol-rank-failure/1\"") != std::string::npos);
      CHECK(text.find("\"rank\": " + std::to_string(other)) != std::string::npos);
      CHECK(text.find("\"status\": \"failed\"") != std::string::npos);
      CHECK(text.find("\"stage\": \"simulation\"") != std::string::npos);
      CHECK(text.find("A failure away from rank 0.") != std::string::npos);
    }
  }

  seissol::Mpi::barrier(seissol::Mpi::mpi.comm());
  if (rank != 0) {
    std::filesystem::remove(prefix + "-rank" + std::to_string(rank) + "-failure.json");
  } else {
    std::filesystem::remove(summaryPath);
    std::filesystem::remove(prefix + "-epochs.csv");
  }
}

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_MONITORING_FAILUREREPORTMPI_T_H_
