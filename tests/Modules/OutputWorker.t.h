// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Modules/OutputWorker.h"

#include <async/as/Pin.h>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <pthread.h>
#include <thread>
#include <vector>

namespace seissol::unit_test {

TEST_CASE("OutputWorker" * doctest::test_suite("modules")) {
  async::as::CpuMask cpus{};
  cpus.getaffinity_np(pthread_self());
  CHECK(OutputWorker::cpuCount(cpus) >= 1);

  std::vector<std::size_t> order;
  std::atomic<bool> otherThread{true};
  const auto mainThread = std::this_thread::get_id();
  {
    OutputWorker worker(cpus, 1);

    // the tasks run one after the other, in the order they were enqueued in, on the thread of the
    // worker
    for (std::size_t i = 0; i < 16; ++i) {
      worker.enqueue([&order, &otherThread, mainThread, i] {
        if (i % 4 == 0) {
          std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
        otherThread = otherThread && std::this_thread::get_id() != mainThread;
        order.push_back(i);
      });
    }
    worker.wait();
    REQUIRE(order.size() == 16);
    for (std::size_t i = 0; i < order.size(); ++i) {
      CHECK(order[i] == i);
    }
    CHECK(otherThread);

    // waiting without tasks returns, and the worker takes further tasks after a wait; the ones left
    // when it ends run before
    worker.wait();
    worker.enqueue([&order] { order.push_back(16); });
    worker.enqueue([&order] {
      std::this_thread::sleep_for(std::chrono::milliseconds(2));
      order.push_back(17);
    });
  }
  REQUIRE(order.size() == 18);
  CHECK(order[16] == 16);
  CHECK(order[17] == 17);
}

} // namespace seissol::unit_test
