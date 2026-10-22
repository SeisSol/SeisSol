// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

// Every rank records a few spans, inside and outside of the window of the trace and beyond its
// limit, and all ranks write them into one file, which rank 0 then reads back.

#include <doctest.h>

#include "Monitoring/Trace.h"
#include "Parallel/MPI.h"

#include <cstddef>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <string>

namespace seissol::unit_test {

namespace tracetest {

inline std::size_t occurrences(const std::string& text, const std::string& part) {
  std::size_t count = 0;
  for (auto at = text.find(part); at != std::string::npos; at = text.find(part, at + 1)) {
    ++count;
  }
  return count;
}

} // namespace tracetest

TEST_CASE("The trace of all ranks goes into one file, within its window and its limit" *
          doctest::test_suite("mpi")) {
  using monitoring::Trace;
  const auto& mpi = seissol::Mpi::mpi;
  const int rank = mpi.rank();
  const int size = mpi.size();
  const std::string path = "seissol-test-trace.json";

  // the spans from 1 s to 2 s after the start, three per rank at most
  Trace::Settings settings;
  settings.enabled = true;
  settings.start = 1;
  settings.duration = 1;
  settings.maxSpans = 3;
  Trace::start(settings, mpi.comm());
  REQUIRE(Trace::recording());

  const auto name = Trace::name("test span");
  CHECK(Trace::name("test span") == name);
  const int key = 0;
  const auto track = Trace::track(&key, "test track");
  CHECK(Trace::track(&key, "another name") == track);

  // before the window
  Trace::span(name, Trace::threadTrack(), 0.5e6, 0.6e6);
  // into the window: it shows from its start
  Trace::span(name, track, 0.9e6, 1.1e6);
  Trace::span(name, track, 1.2e6, 1.3e6);
  Trace::span(name, Trace::threadTrack(), 1.5e6, 1.6e6);
  // beyond the limit
  Trace::span(name, track, 1.7e6, 1.8e6);
  // after the window
  Trace::span(name, track, 2.5e6, 2.6e6);

  Trace::write(path, mpi.comm());
  CHECK_FALSE(Trace::recording());

  if (rank == 0) {
    std::ifstream file(path);
    const std::string text{std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
    file.close();
    CHECK(text.rfind("{\"traceEvents\":[\n", 0) == 0);
    const std::string end = "\n],\"displayTimeUnit\":\"ms\"}\n";
    REQUIRE(text.size() > end.size());
    CHECK(text.compare(text.size() - end.size(), end.size(), end) == 0);
    for (int process = 0; process < size; ++process) {
      CAPTURE(process);
      CHECK(tracetest::occurrences(
                text, "\"args\":{\"name\":\"rank " + std::to_string(process) + "\"}") == 1);
    }
    CHECK(tracetest::occurrences(text, "\"ph\":\"X\"") == static_cast<std::size_t>(3 * size));
    CHECK(tracetest::occurrences(text, "{\"name\":\"test track\"}") ==
          static_cast<std::size_t>(size));
    CHECK(tracetest::occurrences(text, "\"ts\":1000000.000,\"dur\":100000.000") ==
          static_cast<std::size_t>(size));
    CHECK(tracetest::occurrences(text, "\"ts\":1200000.000,\"dur\":100000.000") ==
          static_cast<std::size_t>(size));
    CHECK(tracetest::occurrences(text, "\"ts\":500000.000") == 0);
    CHECK(tracetest::occurrences(text, "\"ts\":1700000.000") == 0);
    CHECK(tracetest::occurrences(text, "\"ts\":2500000.000") == 0);
    CHECK(std::remove(path.c_str()) == 0);
  }
}

} // namespace seissol::unit_test
