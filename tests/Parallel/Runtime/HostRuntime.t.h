// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Parallel/Runtime/HostRuntime.h"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <mutex>
#include <utility>
#include <vector>

namespace seissol::unit_test {

using parallel::runtime::HostRuntime;

TEST_CASE("A static split covers the range once, evenly" * doctest::test_suite("parallel")) {
  for (const std::size_t count : {0, 1, 5, 12, 13, 100}) {
    for (const std::size_t threads : {1, 3, 12, 16}) {
      std::size_t next = 0;
      for (std::size_t thread = 0; thread < threads; ++thread) {
        const auto [begin, end] = HostRuntime::staticRange(count, thread, threads);
        // the ranges follow each other, and differ in size by one at most
        REQUIRE(begin == next);
        REQUIRE(end >= begin);
        CHECK(end - begin <= count / threads + 1);
        CHECK(end - begin >= count / threads);
        next = end;
      }
      CHECK(next == count);
    }
  }
}

TEST_CASE("The loops of the host runtime visit each element once" *
          doctest::test_suite("parallel")) {
  for (const std::size_t count : {0, 1, 7, 1000}) {
    std::vector<std::atomic<int>> visits(count);
    HostRuntime::forEach(count, [&](std::size_t i) { visits[i].fetch_add(1); });
    CHECK(std::all_of(visits.begin(), visits.end(), [](const auto& v) { return v.load() == 1; }));

    std::mutex mutex;
    std::vector<std::pair<std::size_t, std::size_t>> ranges;
    HostRuntime::forRange(count, [&](std::size_t begin, std::size_t end) {
      const std::scoped_lock lock(mutex);
      ranges.emplace_back(begin, end);
    });
    std::sort(ranges.begin(), ranges.end());
    std::size_t next = 0;
    for (const auto& [begin, end] : ranges) {
      // no empty ranges, none overlapping, none missing
      REQUIRE(begin == next);
      REQUIRE(end > begin);
      next = end;
    }
    CHECK(next == count);
  }
}

TEST_CASE("The work of an action runs right away" * doctest::test_suite("parallel")) {
  HostRuntime runtime;
  int ran = 0;
  runtime.action({/*readsNeighbors=*/true, /*writesShared=*/false}, [&ran]() { ++ran; });
  CHECK(ran == 1);
}

} // namespace seissol::unit_test
