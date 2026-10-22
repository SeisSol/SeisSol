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
#include <chrono>
#include <cstddef>
#include <mutex>
#include <thread>
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

TEST_CASE("The host runtime counts the time of the threads in loops" *
          doctest::test_suite("parallel")) {
  HostRuntime::resetLoopTime();
  HostRuntime::forRange(4, [](std::size_t /*begin*/, std::size_t /*end*/) {
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
  });
  // at least one range ran, for at least its sleep
  CHECK(HostRuntime::loopTime() >= 0.01);
  HostRuntime::resetLoopTime();
  CHECK(HostRuntime::loopTime() == doctest::Approx(0.0));
}

TEST_CASE("The work of an action runs right away" * doctest::test_suite("parallel")) {
  HostRuntime runtime;
  int ran = 0;
  runtime.action({/*readsNeighbors=*/true, /*writesShared=*/false}, [&ran]() { ++ran; });
  CHECK(ran == 1);
}

namespace hostruntime {
/**
 * Lets all host runtimes schedule in a task graph for the scope of a test.
 */
struct TaskGraphScope {
  explicit TaskGraphScope(std::size_t taskSize) {
    HostRuntime::setScheduling(parallel::runtime::HostScheduling::TaskGraph, taskSize);
  }
  ~TaskGraphScope() { HostRuntime::setScheduling(parallel::runtime::HostScheduling::Static, 32); }
  TaskGraphScope(const TaskGraphScope&) = delete;
  TaskGraphScope(TaskGraphScope&&) = delete;
  TaskGraphScope& operator=(const TaskGraphScope&) = delete;
  TaskGraphScope& operator=(TaskGraphScope&&) = delete;
};
} // namespace hostruntime

TEST_CASE("In a task graph, the loops of the work run in tasks of the task size" *
          doctest::test_suite("parallel")) {
  if (!HostRuntime::taskGraphAvailable()) {
    return;
  }
  const hostruntime::TaskGraphScope scope(7);
  HostRuntime runtime;
  std::mutex mutex;
  std::vector<std::pair<std::size_t, std::size_t>> ranges;
  HostRuntime::epoch([&]() {
    runtime.action({}, [&]() {
      HostRuntime::forRange(100, [&](std::size_t begin, std::size_t end) {
        const std::scoped_lock lock(mutex);
        ranges.emplace_back(begin, end);
      });
    });
  });
  std::sort(ranges.begin(), ranges.end());
  REQUIRE(ranges.size() == 15);
  for (std::size_t i = 0; i < ranges.size(); ++i) {
    CHECK(ranges[i].first == 7 * i);
    CHECK(ranges[i].second == std::min<std::size_t>(100, 7 * (i + 1)));
  }
}

TEST_CASE("In a task graph, the work of the actions follows what it reads and writes" *
          doctest::test_suite("parallel")) {
  if (!HostRuntime::taskGraphAvailable()) {
    return;
  }
  const hostruntime::TaskGraphScope scope(32);
  HostRuntime first;
  HostRuntime second;
  first.observe(second);
  second.observe(first);

  std::atomic<int> step{0};
  std::vector<int> order(5, -1);
  const auto record = [&](std::size_t index) { order[index] = step.fetch_add(1); };
  constexpr parallel::runtime::HostAccess Predict{/*readsNeighbors=*/false, /*writesShared=*/true};
  constexpr parallel::runtime::HostAccess Correct{/*readsNeighbors=*/true, /*writesShared=*/false};

  HostRuntime::epoch([&]() {
    // the first cluster predicts slowly; the second one corrects with what it shares
    first.action(Predict, [&]() {
      std::this_thread::sleep_for(std::chrono::milliseconds(20));
      record(0);
    });
    second.action(Correct, [&]() { record(1); });
    // overwriting what the second cluster has read waits for it
    first.action(Predict, [&]() { record(2); });
    // the second cluster's own work comes in order
    second.action(Predict, [&]() { record(3); });
    second.action(Correct, [&]() { record(4); });
  });

  CHECK(order[0] < order[1]);
  CHECK(order[1] < order[2]);
  CHECK(order[1] < order[3]);
  CHECK(order[3] < order[4]);
  CHECK(order[2] < order[4]);
}

TEST_CASE("In a task graph, the main thread stays within the limit of waiting actions" *
          doctest::test_suite("parallel")) {
  if (!HostRuntime::taskGraphAvailable()) {
    return;
  }
  const hostruntime::TaskGraphScope scope(32);
  HostRuntime::setActionLimit(2);
  HostRuntime runtime;
  std::size_t mostWaiting = 0;
  std::atomic<int> done{0};
  HostRuntime::epoch([&]() {
    for (int i = 0; i < 20; ++i) {
      HostRuntime::throttle([]() {});
      mostWaiting = std::max(mostWaiting, HostRuntime::outstandingActions());
      runtime.action({}, [&done]() {
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
        ++done;
      });
    }
  });
  HostRuntime::setActionLimit(0);
  CHECK(mostWaiting < 2);
  CHECK(done == 20);
  CHECK(runtime.createdActions() == 20);
  CHECK(runtime.completedActions() == 20);
  CHECK(HostRuntime::outstandingActions() == 0);
}

} // namespace seissol::unit_test
