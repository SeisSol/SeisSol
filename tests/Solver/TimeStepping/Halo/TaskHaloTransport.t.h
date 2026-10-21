// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Parallel/Runtime/HostRuntime.h"
#include "Solver/TimeStepping/Halo/HaloTransport.h"
#include "Solver/TimeStepping/Halo/TaskHaloTransport.h"

#include <atomic>
#include <chrono>
#include <cstddef>
#include <memory>
#include <optional>
#include <thread>

#ifdef SEISSOL_HOST_TASK_GRAPH

namespace seissol::unit_test {
using namespace seissol::solver;
using parallel::runtime::HostAccess;
using parallel::runtime::HostRuntime;
using parallel::runtime::HostScheduling;

namespace taskhalotransport {

/**
 * A transport whose sends and receives complete after a given number of tests each, and which
 * records what it has seen of the data when it started.
 */
class SlowTransport : public HaloTransport {
  public:
  explicit SlowTransport(const std::atomic<int>& data) : data_(data) {}

  void startSend() override {
    REQUIRE_FALSE(sending);
    sending = true;
    sendTests_ = 0;
    dataAtSend = data_.load();
  }
  bool testSend() override {
    if (sending && ++sendTests_ > Latency) {
      sending = false;
      ++completedSends;
    }
    return !sending;
  }
  void startReceive() override {
    REQUIRE_FALSE(receiving);
    receiving = true;
    receiveTests_ = 0;
  }
  bool testReceive() override {
    if (receiving && ++receiveTests_ > Latency) {
      receiving = false;
      ++completedReceives;
    }
    return !receiving;
  }

  static constexpr int Latency = 50;
  bool sending{false};
  bool receiving{false};
  std::atomic<int> completedSends{0};
  std::atomic<int> completedReceives{0};
  int dataAtSend{-1};

  private:
  const std::atomic<int>& data_;
  int sendTests_{0};
  int receiveTests_{0};
};

/**
 * Lets all host runtimes schedule in a task graph for the scope of a test.
 */
struct TaskGraphScope {
  TaskGraphScope() { HostRuntime::setScheduling(HostScheduling::TaskGraph, 32); }
  ~TaskGraphScope() { HostRuntime::setScheduling(HostScheduling::Static, 32); }
  TaskGraphScope(const TaskGraphScope&) = delete;
  TaskGraphScope(TaskGraphScope&&) = delete;
  TaskGraphScope& operator=(const TaskGraphScope&) = delete;
  TaskGraphScope& operator=(TaskGraphScope&&) = delete;
};

/// the main thread keeps the exchanges going between the tasks it takes
inline void awaitExchanges(TaskHaloEngine& engine) {
  while (engine.outstanding() > 0) {
    engine.poll();
#pragma omp taskyield
  }
}

/**
 * Where the OpenMP runtime may wait for a task of an exchange where it gets created, a thread of
 * its own polls the exchanges, as the communication thread does.
 */
class PollingThread {
  public:
  PollingThread(TaskHaloEngine& engine, const std::atomic<bool>& created)
      : thread_([&engine, &created]() {
          while (!created.load() || engine.outstanding() > 0) {
            engine.poll();
          }
        }) {}
  ~PollingThread() { thread_.join(); }
  PollingThread(const PollingThread&) = delete;
  PollingThread(PollingThread&&) = delete;
  PollingThread& operator=(const PollingThread&) = delete;
  PollingThread& operator=(PollingThread&&) = delete;

  private:
  std::thread thread_;
};

} // namespace taskhalotransport

TEST_CASE("In a task graph, the halo exchange orders itself with the work on the layers" *
          doctest::test_suite("solver")) {
  using namespace taskhalotransport;
  const TaskGraphScope scope;

  // the copy layer reads the ghost layer in its corrections; the ghost cluster sends the copy layer
  HostRuntime copy;
  HostRuntime ghost;
  copy.observe(ghost);
  ghost.observe(copy);

  std::atomic<int> copyData{0};
  TaskHaloEngine engine;
  auto slow = std::make_unique<SlowTransport>(copyData);
  auto* transport = slow.get();
  TaskHaloTransport tasks(std::move(slow), engine);
  tasks.setHostTokens(ghost.neighborToken(0), ghost.sharedToken());
  CHECK(tasks.streamOrdered());

  std::atomic<int> receivedWhenCorrected{-1};
  std::atomic<int> sentWhenOverwritten{-1};
  std::atomic<bool> created{false};
  std::optional<PollingThread> polling;
  if (!HostRuntime::detachedTasksNeverBlock()) {
    polling.emplace(engine, created);
  }
  HostRuntime::epoch([&]() {
    // the copy layer gets written slowly; the send has to wait for it
    copy.action({/*readsNeighbors=*/false, /*writesShared=*/true}, [&]() {
      std::this_thread::sleep_for(std::chrono::milliseconds(20));
      copyData = 1;
    });
    tasks.startSend();
    CHECK(tasks.testSend());

    // the correction reads the ghost layer, which has to be there by then
    tasks.startReceive();
    CHECK(tasks.testReceive());
    copy.action({/*readsNeighbors=*/true, /*writesShared=*/false},
                [&]() { receivedWhenCorrected = transport->completedReceives.load(); });

    // overwriting the copy layer waits until it has been sent
    copy.action({/*readsNeighbors=*/false, /*writesShared=*/true}, [&]() {
      sentWhenOverwritten = transport->completedSends.load();
      copyData = 2;
    });

    created = true;
    if (!polling.has_value()) {
      awaitExchanges(engine);
    }
  });
  polling.reset();

  CHECK(transport->dataAtSend == 1);
  CHECK(receivedWhenCorrected == 1);
  CHECK(sentWhenOverwritten == 1);
  CHECK(engine.outstanding() == 0);
}

} // namespace seissol::unit_test

#endif // SEISSOL_HOST_TASK_GRAPH
