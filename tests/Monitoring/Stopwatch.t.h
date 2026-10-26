// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Monitoring/Stopwatch.h"
#include "TestHelper.h"

#include <thread>

namespace seissol::unit_test {

TEST_CASE("difftime helper" * doctest::test_suite("monitoring")) {
  timespec start{};
  timespec end{};

  SUBCASE("Same time") {
    start.tv_sec = 10;
    start.tv_nsec = 0;
    end.tv_sec = 10;
    end.tv_nsec = 0;
    CHECK(seissol::difftime(start, end) == 0);
  }

  SUBCASE("One second") {
    start.tv_sec = 10;
    start.tv_nsec = 0;
    end.tv_sec = 11;
    end.tv_nsec = 0;
    CHECK(seissol::difftime(start, end) == 1000000000LL);
  }

  SUBCASE("Subsecond") {
    start.tv_sec = 10;
    start.tv_nsec = 500000000;
    end.tv_sec = 10;
    end.tv_nsec = 750000000;
    CHECK(seissol::difftime(start, end) == 250000000LL);
  }

  SUBCASE("Cross-second boundary") {
    start.tv_sec = 10;
    start.tv_nsec = 900000000;
    end.tv_sec = 11;
    end.tv_nsec = 100000000;
    CHECK(seissol::difftime(start, end) == 200000000LL);
  }
}

TEST_CASE("seconds conversion" * doctest::test_suite("monitoring")) {
  CHECK(seissol::seconds(1000000000LL) == doctest::Approx(1.0));
  CHECK(seissol::seconds(500000000LL) == doctest::Approx(0.5));
  CHECK(seissol::seconds(0LL) == doctest::Approx(0.0));
}

TEST_CASE("Stopwatch basic operations" * doctest::test_suite("monitoring")) {

  SUBCASE("Split returns positive after start") {
    Stopwatch sw;
    sw.start();
    // Tiny sleep to ensure measurable time passes
    std::this_thread::sleep_for(std::chrono::milliseconds(5));
    const double elapsed = sw.split();
    CHECK(elapsed > 0.0);
    CHECK(elapsed < 1.0); // sanity: should not be more than 1s
  }

  SUBCASE("Pause accumulates") {
    Stopwatch sw;
    sw.start();
    std::this_thread::sleep_for(std::chrono::milliseconds(5));
    const double first = sw.pause();
    CHECK(first > 0.0);

    sw.start();
    std::this_thread::sleep_for(std::chrono::milliseconds(5));
    const double second = sw.pause();
    CHECK(second > first);
  }

  SUBCASE("Reset clears accumulated time") {
    Stopwatch sw;
    sw.start();
    std::this_thread::sleep_for(std::chrono::milliseconds(5));
    sw.pause();

    sw.reset();
    sw.start();
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
    const double afterReset = sw.pause();
    // After reset, time should be very small (just the 1ms sleep)
    CHECK(afterReset < 0.05);
  }
}

} // namespace seissol::unit_test
