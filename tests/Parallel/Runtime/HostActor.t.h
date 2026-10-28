// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Parallel/Runtime/Event.h"
#include "Parallel/Runtime/HostActor.h"
#include "Parallel/Runtime/HostBufferPool.h"
#include "Parallel/Runtime/Stream.h"

#include <array>
#include <atomic>
#include <cstddef>
#include <functional>
#include <future>
#include <mutex>
#include <set>
#include <string>
#include <utility>
#include <vector>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::unit_test {

using parallel::runtime::Event;
using parallel::runtime::HostActor;
using parallel::runtime::HostBufferPool;

namespace hostactortest {
/// events that complete when the test says so, standing for the work of a device before them
struct Device {
  std::array<int, 8> events{};
  std::set<const void*> completed;

  Event event(std::size_t index) { return Event(&events.at(index)); }

  void complete(std::size_t index) { completed.insert(&events.at(index)); }
};

inline void setUpDevice() {
#ifdef ACL_DEVICE
  auto& device = ::device::DeviceInstance::instance();
  device.api().setDevice(0);
  device.api().initialize();
#endif
}
} // namespace hostactortest

using hostactortest::Device;

TEST_CASE("The host actor runs the work of a queue in order, once its events have completed" *
          doctest::test_suite("parallel")) {
  Device device;
  HostActor actor([&](const Event& event) { return device.completed.count(event.get()) > 0; });
  HostActor::Queue queue(actor);
  std::vector<int> ran;
  queue.post(device.event(0), [&]() { ran.push_back(0); });
  queue.post(device.event(1), [&]() { ran.push_back(1); });
  queue.post(Event(), [&]() { ran.push_back(2); });

  CHECK_FALSE(actor.progress());
  CHECK(ran.empty());
  // the second one waits for the first one
  device.complete(1);
  CHECK_FALSE(actor.progress());
  CHECK(ran.empty());
  device.complete(0);
  CHECK(actor.progress());
  // the one without an event runs right after the ones before it
  CHECK(ran == std::vector<int>{0, 1, 2});
  CHECK(queue.empty());
  CHECK(actor.pending() == 0);
}

TEST_CASE("The queues of the host actor do not wait for each other" *
          doctest::test_suite("parallel")) {
  Device device;
  HostActor actor([&](const Event& event) { return device.completed.count(event.get()) > 0; });
  HostActor::Queue first(actor);
  HostActor::Queue second(actor);
  std::vector<std::string> ran;
  first.post(device.event(0), [&]() { ran.emplace_back("first"); });
  second.post(device.event(1), [&]() { ran.emplace_back("second"); });
  CHECK(actor.pending() == 2);
  device.complete(1);
  CHECK(actor.progress());
  CHECK(ran == std::vector<std::string>{"second"});
  CHECK(first.pending() == 1);
  device.complete(0);
  first.drain();
  CHECK(ran == std::vector<std::string>{"second", "first"});
}

TEST_CASE("Draining the host actor waits for the device" * doctest::test_suite("parallel")) {
  Device device;
  std::size_t queries = 0;
  // the device completes the event after a few queries, as if it got there meanwhile
  HostActor actor([&](const Event& /*event*/) { return ++queries > 3; });
  std::vector<int> ran;
  {
    HostActor::Queue queue(actor);
    queue.post(device.event(0), [&]() { ran.push_back(0); });
    queue.post(device.event(1), [&]() { ran.push_back(1); });
    actor.drain();
    CHECK(ran == std::vector<int>{0, 1});
    CHECK(queries > 3);
    // and a queue that goes runs its work before
    queue.post(device.event(2), [&]() { ran.push_back(2); });
  }
  CHECK(ran == std::vector<int>{0, 1, 2});
  CHECK(actor.pending() == 0);
}

TEST_CASE("Work of the host actor may post more work" * doctest::test_suite("parallel")) {
  HostActor actor([](const Event& /*event*/) { return true; });
  HostActor::Queue queue(actor);
  std::vector<int> ran;
  queue.post(Event(), [&]() {
    ran.push_back(0);
    queue.post(Event(), [&]() { ran.push_back(1); });
  });
  actor.drain();
  CHECK(ran == std::vector<int>{0, 1});
}

TEST_CASE("The host actor hands ready work to its executor, and waits for it to be done" *
          doctest::test_suite("parallel")) {
  HostActor actor([](const Event& /*event*/) { return true; });
  std::vector<std::future<void>> running;
  // runs each piece of work on a thread of its own, one after the other
  std::mutex order;
  actor.setExecutor([&](std::function<void()> task) {
    running.push_back(std::async(std::launch::async, [&order, task = std::move(task)]() {
      const std::lock_guard lock(order);
      task();
    }));
  });
  CHECK(actor.hasExecutor());
  HostActor::Queue queue(actor);
  std::atomic<int> ran{0};
  for (int i = 0; i < 4; ++i) {
    queue.post(Event(), [&ran]() { ran.fetch_add(1); });
  }
  queue.drain();
  CHECK(ran.load() == 4);
  CHECK(queue.empty());
  for (auto& task : running) {
    task.get();
  }
  actor.setExecutor({});
  CHECK_FALSE(actor.hasExecutor());
}

TEST_CASE("The host actor hands each piece of work the data staged for it" *
          doctest::test_suite("parallel")) {
  hostactortest::setUpDevice();
  parallel::runtime::StreamRuntime runtime;
  {
    HostActor::Queue queue(sizeof(int));
    std::vector<int> seen;
    for (int value = 0; value < 3; ++value) {
      // as if the device wrote it
      const auto staging = queue.take(runtime);
      *static_cast<int*>(staging.host) = value;
      queue.post(runtime, staging, [&seen](const void* data) {
        seen.push_back(*static_cast<const int*>(data));
      });
    }
    queue.drain();
    CHECK(seen == std::vector<int>{0, 1, 2});
  }
  runtime.dispose();
}

#ifdef ACL_DEVICE
TEST_CASE("The work for the host in a recording runs after each launch, with the data of that "
          "launch" *
          doctest::test_suite("parallel")) {
  hostactortest::setUpDevice();
  auto& api = ::device::DeviceInstance::instance().api();
  parallel::runtime::StreamRuntime runtime;
  {
    HostActor::Queue first(sizeof(int));
    HostActor::Queue second(2 * sizeof(int));
    HostActor::Recording recording;
    // by queue: the queues do not wait for each other
    std::vector<int> seenFirst;
    std::vector<int> seenSecond;

    // what the recorded work stages, as a graph would keep its addresses
    HostActor::instance().beginRecording(recording);
    const auto a = first.take(runtime);
    const auto b = second.take(runtime);
    first.post(runtime, a, [&seenFirst](const void* data) {
      seenFirst.push_back(*static_cast<const int*>(data));
    });
    second.post(runtime, b, [&seenSecond](const void* data) {
      seenSecond.push_back(static_cast<const int*>(data)[0] + static_cast<const int*>(data)[1]);
    });
    HostActor::instance().endRecording();
    CHECK(a.host == nullptr);
    CHECK_FALSE(recording.empty());
    CHECK(seenFirst.empty());
    CHECK(seenSecond.empty());

    // two launches, each of which writes other data, as a replay would
    for (int launch = 0; launch < 2; ++launch) {
      const int value = 10 * (launch + 1);
      const std::array<int, 2> pair{value, 1};
      api.copyTo(a.device, &value, sizeof(int));
      api.copyTo(b.device, pair.data(), sizeof(pair));
      recording.launched(runtime);
    }
    HostActor::instance().drain();
    CHECK(seenFirst == std::vector<int>{10, 20});
    CHECK(seenSecond == std::vector<int>{11, 21});
  }
  runtime.dispose();
}
#endif

TEST_CASE("A host buffer pool grows only when all of its buffers are in use" *
          doctest::test_suite("parallel")) {
  hostactortest::setUpDevice();
  HostBufferPool pool(100);
  const auto first = pool.take();
  const auto second = pool.take();
  CHECK(first.host != second.host);
  CHECK(pool.size() == 2);
  pool.give(first);
  const auto third = pool.take();
  CHECK(third.host == first.host);
  CHECK(pool.size() == 2);
  pool.give(second);
  pool.give(third);
}

} // namespace seissol::unit_test
