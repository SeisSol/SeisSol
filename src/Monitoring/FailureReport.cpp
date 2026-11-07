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
#include <cerrno>
#include <cstddef>
#include <cstdio>
#include <ctime>
#include <fcntl.h>
#include <string>
#include <string_view>
#include <unistd.h>
#include <utility>
#include <utils/logger.h>

namespace seissol::monitoring {

namespace {

struct Report {
  std::string path;
  std::string temporary;
  std::string before;
  std::string after;
};

// (constant-initialized, so that a signal handler finds them in any case)
// the report a failure writes; replaced rather than changed
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::atomic<const Report*> current{nullptr};
// whether a failure has written the report, or is writing it; it stays set
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::atomic<bool> claimed{false};

void replaceReport(const Report* report) {
  const Report* previous = current.exchange(report);
  // a failure that writes the report now may be reading the one before; then it is kept
  if (!claimed.load()) {
    delete previous;
  }
}

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

} // namespace

void setFailureReport(const std::string& path, std::string before, std::string after) {
  // a name of its own, so that it does not meet the file of writeFileAtomically() for the path
  replaceReport(new Report{.path = path,
                           .temporary = path + ".failure.partial",
                           .before = std::move(before),
                           .after = std::move(after)});
}

void setFailureReport(const std::string& path, const json::Value& report) {
  const auto text = report.dump() + "\n";
  const auto placeholder = json::quote(FailurePlaceholder);
  const auto at = text.find(placeholder);
  if (at == std::string::npos) {
    clearFailureReport();
    return;
  }
  setFailureReport(path, text.substr(0, at), text.substr(at + placeholder.size()));
}

void clearFailureReport() { replaceReport(nullptr); }

void writeFailureReport(std::string_view failure) {
  if (claimed.exchange(true)) {
    return;
  }
  const Report* report = current.load();
  if (report == nullptr) {
    return;
  }
  // next to the report and renamed then, as writeFileAtomically() does, with what a signal
  // handler may call
  constexpr auto Flags = O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC;
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-vararg)
  const int file = ::open(report->temporary.c_str(), Flags, 0644);
  if (file < 0) {
    return;
  }
  writeAll(file, report->before);
  writeAll(file, failure);
  writeAll(file, report->after);
  ::close(file);
  static_cast<void>(std::rename(report->temporary.c_str(), report->path.c_str()));
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
