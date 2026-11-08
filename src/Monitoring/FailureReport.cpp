// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "FailureReport.h"

#include "Common/Json.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cstddef>
#include <ctime>
#include <string>
#include <string_view>
#include <utils/logger.h>

// A report is written with what a signal handler may call, which are functions of POSIX. Windows
// has no signals that utils catches (see utils/crash.h); there, the report is written with the
// files of the C++ library.
#ifdef _WIN32
#include <filesystem>
#include <fstream>
#include <ios>
#include <system_error>
#else
#include <cerrno>
#include <cstdio>
#include <fcntl.h>
#include <unistd.h>
#endif

namespace seissol::monitoring {

namespace {

struct Report {
  std::string path;
  std::string temporary;
  std::string before;
  std::string after;
};

// (constant-initialized, so that a signal handler finds them in any case)
// the report a failure writes, and the one a termination writes; each replaced rather than changed
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::atomic<const Report*> failureReport{nullptr};
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::atomic<const Report*> terminationReport{nullptr};
// whether a failure or a termination has written a report, or is writing it; it stays set
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::atomic<bool> claimed{false};

void replaceReport(std::atomic<const Report*>& slot, const Report* report) {
  const Report* previous = slot.exchange(report);
  // a signal handler that writes the report now may be reading the one before; then it is kept
  if (!claimed.load()) {
    delete previous;
  }
}

#ifndef _WIN32
void writeAll(int file, std::string_view text) {
  while (!text.empty()) {
    const auto written = ::write(file, text.data(), text.size());
    if (written < 0 && errno == EINTR) {
      continue;
    }
    if (written <= 0) {
      return;
    }
    text.remove_prefix(static_cast<std::size_t>(written));
  }
}
#endif

//! @brief Splits @p report at the placeholder and puts it in @p slot .
void setReport(std::atomic<const Report*>& slot,
               const std::string& path,
               const json::Value& report) {
  const auto text = report.dump() + "\n";
  const auto placeholder = json::quote(FailurePlaceholder);
  const auto at = text.find(placeholder);
  if (at == std::string::npos) {
    replaceReport(slot, nullptr);
    return;
  }
  // a name of its own, so that it does not meet the file of writeFileAtomically() for the path
  replaceReport(slot,
                new Report{.path = path,
                           .temporary = path + ".failure.partial",
                           .before = text.substr(0, at),
                           .after = text.substr(at + placeholder.size())});
}

//! @brief Writes the report in @p slot with @p detail , if no report has been written yet.
void writeReport(const std::atomic<const Report*>& slot, std::string_view detail) {
  const Report* report = slot.load();
  if (report == nullptr) {
    return;
  }
  // the first failure or termination writes it, and no more
  if (claimed.exchange(true)) {
    return;
  }
  // next to the report and renamed then, as writeFileAtomically() does
#ifdef _WIN32
  {
    std::ofstream stream(report->temporary, std::ios::out | std::ios::trunc | std::ios::binary);
    stream << report->before << detail << report->after;
    if (!stream.good()) {
      return;
    }
  }
  // (std::rename() does not replace a file on Windows)
  std::error_code error;
  std::filesystem::rename(report->temporary, report->path, error);
#else
  // with what a signal handler may call
  constexpr auto Flags = O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC;
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-vararg)
  const int file = ::open(report->temporary.c_str(), Flags, 0644);
  if (file < 0) {
    return;
  }
  writeAll(file, report->before);
  writeAll(file, detail);
  writeAll(file, report->after);
  ::close(file);
  static_cast<void>(std::rename(report->temporary.c_str(), report->path.c_str()));
#endif
}

} // namespace

void setFailureReport(const std::string& path, const json::Value& report) {
  setReport(failureReport, path, report);
}

void setTerminationReport(const std::string& path, const json::Value& report) {
  setReport(terminationReport, path, report);
}

void clearFailureReports() {
  replaceReport(failureReport, nullptr);
  replaceReport(terminationReport, nullptr);
}

void writeFailureReport(std::string_view failure) { writeReport(failureReport, failure); }

void writeTerminationReport(std::string_view termination) {
  writeReport(terminationReport, termination);
}

std::string failureOf(const utils::LogRecord& record) {
  auto failure = json::Value::object();
  failure["message"] = std::string(record.message());
  auto& source = failure["source"];
  source["file"] = std::string(record.location.file);
  source["line"] = record.location.line;
  source["function"] = std::string(record.location.function);

  const auto seconds = static_cast<std::time_t>(record.time.seconds);
  std::tm utc{};
#ifdef _WIN32
  gmtime_s(&utc, &seconds);
#else
  gmtime_r(&seconds, &utc);
#endif
  std::array<char, 32> buffer{};
  const auto length = std::strftime(buffer.data(), buffer.size(), "%Y-%m-%dT%H:%M:%S", &utc);
  // to the microsecond, as the JSON lines of the log
  auto microseconds = std::to_string(record.time.nanoseconds / 1000);
  microseconds.insert(0, 6 - std::min<std::size_t>(6, microseconds.size()), '0');
  failure["time"] = std::string(buffer.data(), length) + "." + microseconds + "Z";
  return failure.dump(-1);
}

} // namespace seissol::monitoring
