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
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <sstream>
#include <string>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>
#include <utils/crash.h>
#include <utils/logger.h>

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

TEST_CASE("A fatal error as a failure" * doctest::test_suite("monitoring")) {
  CHECK(monitoring::failureOf(fatalError()) ==
        R"({"message":"The error.","source":{"file":"Init.cpp","line":12,"function":"initMesh"},)"
        R"("time":"2023-11-14T22:13:20.123456Z"})");
}

TEST_CASE("The failure report is written by the first failure" *
          doctest::test_suite("monitoring")) {
  const std::string path = "failure-report-test.json";
  std::filesystem::remove(path);

  SUBCASE("Set") {
    inChild([&path] {
      monitoring::setFailureReport(path, R"({"status": "failed", "error": )", "}\n");
      monitoring::writeFailureReport(R"({"message":"first"})");
      monitoring::writeFailureReport(R"({"message":"second"})");
    });
    CHECK(fileContent(path) == "{\"status\": \"failed\", \"error\": {\"message\":\"first\"}}\n");
    CHECK_FALSE(std::filesystem::exists(path + ".failure.partial"));
  }
  SUBCASE("Cleared") {
    inChild([&path] {
      monitoring::setFailureReport(path, "{", "}");
      monitoring::clearFailureReport();
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

TEST_CASE("A crash marks the run summary as failed" * doctest::test_suite("monitoring")) {
  const std::string prefix = "crash-summary-test";
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
    volatile int* volatile pointer = nullptr;
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
  CHECK(text.find(R"("signal":"SIGSEGV","address":"0x0","time":")") != std::string::npos);

  summary.finish();
  std::filesystem::remove(path);
  std::filesystem::remove(prefix + "-epochs.csv");
}

} // namespace seissol::unit_test
