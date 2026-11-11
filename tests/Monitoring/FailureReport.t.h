// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Monitoring/FailureReport.h"
#include "Monitoring/RunSummary.h"

#include <csignal>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <sstream>
#include <string>
#include <utils/crash.h>
#include <utils/logger.h>

// a report is written once per process, so that the tests that write one run it in a child; Windows
// has no fork()
#ifndef _WIN32
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

namespace seissol::unit_test {

inline std::string fileContent(const std::string& path) {
  const std::ifstream stream(path);
  std::stringstream content;
  content << stream.rdbuf();
  return content.str();
}

inline utils::LogRecord fatalError() {
  utils::LogRecord record;
  record.level = utils::LogLevel::Error;
  record.text = "The error. ";
  record.location = {.file = "Init.cpp", .function = "initMesh", .line = 12};
  record.time = {.seconds = 1700000000, .nanoseconds = 123456000};
  return record;
}

TEST_CASE("A fatal error as a failure" * doctest::test_suite("monitoring")) {
  CHECK(monitoring::failureOf(fatalError()) ==
        R"({"message":"The error.","source":{"file":"Init.cpp","line":12,"function":"initMesh"},)"
        R"("time":"2023-11-14T22:13:20.123456Z"})");
}

#ifndef _WIN32
/**
 * @brief Runs @p function in a child process, which ends after it, and waits for it: a failure is
 * written once per process.
 */
template <typename F>
void inChild(F function) {
  // what is still buffered would be written by both processes
  static_cast<void>(std::fflush(nullptr));
  const pid_t child = fork();
  if (child == 0) {
    function();
    _exit(0);
  }
  int status = 0;
  REQUIRE(waitpid(child, &status, 0) == child);
  CHECK(WIFEXITED(status));
}

inline json::Value reportTemplate(const char* status, const char* key) {
  auto report = json::Value::object();
  report["status"] = std::string(status);
  report[key] = monitoring::FailurePlaceholder;
  return report;
}

TEST_CASE("The report is written by the first failure or termination" *
          doctest::test_suite("monitoring")) {
  const std::string path = "failure-report-test.json";
  std::filesystem::remove(path);

  SUBCASE("The first failure writes it, and no more") {
    inChild([&path] {
      monitoring::setFailureReport(path, reportTemplate("failed", "error"));
      monitoring::writeFailureReport(R"({"message":"first"})");
      monitoring::writeFailureReport(R"({"message":"second"})");
    });
    CHECK(fileContent(path) ==
          "{\n  \"status\": \"failed\",\n  \"error\": {\"message\":\"first\"}\n}\n");
    CHECK_FALSE(std::filesystem::exists(path + ".failure.partial"));
  }
  SUBCASE("A failure wins over a later termination") {
    inChild([&path] {
      monitoring::setFailureReport(path, reportTemplate("failed", "error"));
      monitoring::setTerminationReport(path, reportTemplate("terminated", "terminatedBy"));
      monitoring::writeFailureReport(R"({"message":"the error"})");
      monitoring::writeTerminationReport(R"({"message":"the signal"})");
    });
    const auto text = fileContent(path);
    CHECK(text.find("\"status\": \"failed\"") != std::string::npos);
    CHECK(text.find("\"status\": \"terminated\"") == std::string::npos);
  }
  SUBCASE("A termination writes the terminated report") {
    inChild([&path] {
      monitoring::setTerminationReport(path, reportTemplate("terminated", "terminatedBy"));
      monitoring::writeTerminationReport(R"({"signal":"SIGTERM"})");
    });
    const auto text = fileContent(path);
    CHECK(text.find("\"status\": \"terminated\"") != std::string::npos);
    CHECK(text.find("\"terminatedBy\": {\"signal\":\"SIGTERM\"}") != std::string::npos);
  }
  SUBCASE("A termination with no report set writes nothing") {
    inChild([&path] {
      monitoring::setFailureReport(path, reportTemplate("failed", "error"));
      monitoring::writeTerminationReport(R"({"signal":"SIGTERM"})");
    });
    CHECK_FALSE(std::filesystem::exists(path));
  }
  SUBCASE("Cleared") {
    inChild([&path] {
      monitoring::setFailureReport(path, reportTemplate("failed", "error"));
      monitoring::clearFailureReports();
      monitoring::writeFailureReport(R"({"message":"first"})");
    });
    CHECK_FALSE(std::filesystem::exists(path));
  }

  std::filesystem::remove(path);
}

TEST_CASE("A failure marks the run summary as failed" * doctest::test_suite("monitoring")) {
  const std::string prefix = "failure-summary-test";
  const auto path = prefix + "-summary.json";

  monitoring::RunSummary summary;
  summary.start(prefix, "123e4567-e89b-42d3-a456-426614174000");
  summary.setStage("simulation");
  inChild([] { monitoring::writeFailureReport(monitoring::failureOf(fatalError())); });
  {
    const auto text = fileContent(path);
    CHECK(text.find("\"status\": \"failed\"") != std::string::npos);
    // after the stage, which says how far the run got
    CHECK(text.find("\"stage\": \"simulation\",\n  \"error\": {\"message\":\"The error.\",") !=
          std::string::npos);
    CHECK(text.find("\"updated\": ") != std::string::npos);
  }

  // a failure after the end leaves the run completed
  summary.finish();
  inChild([] { monitoring::writeFailureReport(monitoring::failureOf(fatalError())); });
  {
    const auto text = fileContent(path);
    CHECK(text.find("\"status\": \"completed\"") != std::string::npos);
    CHECK(text.find("\"error\"") == std::string::npos);
  }

  std::filesystem::remove(path);
  std::filesystem::remove(prefix + "-epochs.csv");
}

TEST_CASE("A termination marks the run summary as terminated" * doctest::test_suite("monitoring")) {
  const std::string prefix = "termination-summary-test";
  const auto path = prefix + "-summary.json";

  monitoring::RunSummary summary;
  summary.start(prefix, "123e4567-e89b-42d3-a456-426614174000");
  summary.setStage("simulation");
  inChild([] {
    monitoring::writeTerminationReport(R"({"message":"The process received the signal SIGTERM )"
                                       R"((termination request).","signal":"SIGTERM"})");
  });
  const auto text = fileContent(path);
  CHECK(text.find("\"status\": \"terminated\"") != std::string::npos);
  CHECK(text.find("\"stage\": \"simulation\",\n  \"terminatedBy\": {\"message\":") !=
        std::string::npos);
  CHECK(text.find("\"signal\":\"SIGTERM\"") != std::string::npos);
  CHECK(text.find("\"status\": \"failed\"") == std::string::npos);

  summary.finish();
  std::filesystem::remove(path);
  std::filesystem::remove(prefix + "-epochs.csv");
}

TEST_CASE("A crash marks the run summary as failed" * doctest::test_suite("monitoring")) {
  const std::string prefix = "crash-summary-test";
  constexpr std::uintptr_t FaultAddress = 16;
  const auto path = prefix + "-summary.json";

  monitoring::RunSummary summary;
  summary.start(prefix, "123e4567-e89b-42d3-a456-426614174000");
  summary.setStage("simulation");
  static_cast<void>(std::fflush(nullptr));
  const pid_t child = fork();
  if (child == 0) {
    const rlimit noCore{.rlim_cur = 0, .rlim_max = 0};
    setrlimit(RLIMIT_CORE, &noCore);
    // doctest has handlers for these, which catchSignals() would keep
    for (const int number : {SIGSEGV, SIGBUS, SIGFPE, SIGILL, SIGABRT}) {
      static_cast<void>(std::signal(number, SIG_DFL));
    }
    utils::catchSignals(
        [](const utils::LogSignal& signal) { monitoring::writeFailureReport(signal.json); });
    // a store into the lowest page, which is never mapped; a store to the null pointer would be
    // stopped by the sanitizer of undefined behavior before it faults
    // NOLINTNEXTLINE(performance-no-int-to-ptr)
    volatile int* volatile pointer = reinterpret_cast<volatile int*>(FaultAddress);
    *pointer = 1;
    _exit(0);
  }
  int status = 0;
  REQUIRE(waitpid(child, &status, 0) == child);
  CHECK(WIFSIGNALED(status));
  CHECK(WTERMSIG(status) == SIGSEGV);
  const auto text = fileContent(path);
  CHECK(text.find("\"status\": \"failed\"") != std::string::npos);
  // (the cause in the message, such as "address not mapped", is the system's to tell)
  CHECK(text.find(R"("stage": "simulation",)"
                  "\n"
                  R"(  "error": {"message":"The process received the signal SIGSEGV )") !=
        std::string::npos);
  CHECK(text.find(R"("signal":"SIGSEGV","address":"0x10","time":")") != std::string::npos);

  summary.finish();
  std::filesystem::remove(path);
  std::filesystem::remove(prefix + "-epochs.csv");
}
#endif

} // namespace seissol::unit_test
