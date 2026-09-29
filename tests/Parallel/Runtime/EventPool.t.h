// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Parallel/Runtime/EventPool.h"

#include <cstddef>
#include <deque>
#include <set>
#include <vector>

namespace seissol::unit_test {

/**
 * Events that are only addresses; keeps track of which ones exist.
 */
struct ModelEventFactory {
  std::deque<char> storage;
  std::set<void*> alive;
  std::set<void*> destroyed;

  parallel::runtime::EventPool pool(std::size_t initialSize) {
    return {[this]() {
              void* event = &storage.emplace_back();
              alive.insert(event);
              return event;
            },
            [this](void* event) {
              alive.erase(event);
              destroyed.insert(event);
            },
            initialSize};
  }
};

TEST_CASE("An event pool hands out an event only while nobody refers to it" *
          doctest::test_suite("parallel")) {
  ModelEventFactory events;
  auto pool = events.pool(4);
  CHECK(pool.size() == 4);
  CHECK(events.alive.size() == 4);

  auto held = pool.next();
  auto copy = held;
  CHECK(pool.referenced() == 1);
  for (int i = 0; i < 20; ++i) {
    const auto other = pool.next();
    REQUIRE(other.get() != held.get());
  }

  // the event comes back once the last reference is gone, and nothing has grown meanwhile
  held = {};
  CHECK(pool.referenced() == 1);
  const void* event = copy.get();
  copy = {};
  CHECK(pool.referenced() == 0);
  bool reused = false;
  for (int i = 0; i < 4; ++i) {
    reused = reused || pool.next().get() == event;
  }
  CHECK(reused);
  CHECK(pool.size() == 4);
  CHECK(pool.peakReferenced() == 2);
}

TEST_CASE("An event pool grows when all its events are referenced" *
          doctest::test_suite("parallel")) {
  ModelEventFactory events;
  auto pool = events.pool(2);
  std::vector<parallel::runtime::EventRef> held;
  held.push_back(pool.next());
  held.push_back(pool.next());
  held.push_back(pool.next());
  CHECK(pool.size() == 2 + parallel::runtime::EventPool::Growth);
  CHECK(pool.referenced() == 3);
  CHECK(pool.peakReferenced() == 3);
  CHECK(std::set<void*>{held[0].get(), held[1].get(), held[2].get()}.size() == 3);
}

TEST_CASE("Disposing an event pool leaves referenced events alone" *
          doctest::test_suite("parallel")) {
  ModelEventFactory events;
  auto pool = events.pool(3);
  auto held = pool.next();
  auto* event = held.get();

  // reported, neither destroyed nor freed
  CHECK(pool.dispose() == 1);
  CHECK(events.destroyed.size() == 2);
  CHECK(events.alive == std::set<void*>{event});
  CHECK(held.get() == event);

  // the reference can still be passed on and dropped
  auto copy = held;
  held = {};
  CHECK(copy.get() == event);
  copy = {};

  // once is enough
  CHECK(pool.dispose() == 0);
  CHECK(events.destroyed.size() == 2);
}

TEST_CASE("Disposing an event pool destroys its events" * doctest::test_suite("parallel")) {
  ModelEventFactory events;
  {
    auto pool = events.pool(3);
    {
      [[maybe_unused]] const auto event = pool.next();
    }
    CHECK(pool.dispose() == 0);
    CHECK(events.alive.empty());
  }
  CHECK(events.destroyed.size() == 3);

  // also when it goes away
  {
    auto pool = events.pool(2);
    [[maybe_unused]] const auto event = pool.next();
  }
  CHECK(events.alive.empty());
  CHECK(events.destroyed.size() == 5);
}

} // namespace seissol::unit_test
