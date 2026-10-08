// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Reader/Datafield/Grid.h"
#include "Reader/Datafield/Interpolation.h"

#include <cstdint>

namespace seissol::unit_test {

using namespace seissol::reader::datafield;

TEST_CASE("Datafield time window placement") {
  const auto cubic = kernelOf(Interpolation::Cubic);   // width 4, offset -1
  const auto linear = kernelOf(Interpolation::Linear); // width 2, offset 0

  SUBCASE("a window that covers the whole axis never moves") {
    const auto w = placeTimeWindow(7, cubic, 100, 20);
    REQUIRE(w.start == 0);
    REQUIRE(w.count == 20);
    for (std::int64_t base = 0; base < 20; ++base) {
      REQUIRE(windowServes(w, base, cubic, 20));
    }
  }

  SUBCASE("the query lands on the first usable slot, headroom ahead") {
    // Cubic needs slot -offset == 1, so the window starts one slice early.
    const auto w = placeTimeWindow(10, cubic, 5, 100);
    REQUIRE(w.start == 9);
    REQUIRE(w.count == 5);
    REQUIRE(windowServes(w, 10, cubic, 100));
    // Usable span is count - width == 1 cell forward, which is exactly what
    // suggestedSyncInterval() promises.
    REQUIRE(windowServes(w, 11, cubic, 100));
    REQUIRE_FALSE(windowServes(w, 12, cubic, 100));
    REQUIRE_FALSE(windowServes(w, 9, cubic, 100));
  }

  SUBCASE("a wider window buys proportionally more headroom") {
    const auto w = placeTimeWindow(10, cubic, 12, 100);
    REQUIRE(w.start == 9);
    for (std::int64_t base = 10; base <= 10 + (12 - 4); ++base) {
      REQUIRE(windowServes(w, base, cubic, 100));
    }
    REQUIRE_FALSE(windowServes(w, 10 + (12 - 4) + 1, cubic, 100));
  }

  SUBCASE("linear needs no lead-in") {
    const auto w = placeTimeWindow(10, linear, 5, 100);
    REQUIRE(w.start == 10);
    REQUIRE(windowServes(w, 10, linear, 100));
    REQUIRE(windowServes(w, 13, linear, 100));
    REQUIRE_FALSE(windowServes(w, 14, linear, 100));
  }

  SUBCASE("placement is clamped to the file, and the ends do not thrash") {
    const auto atStart = placeTimeWindow(0, cubic, 5, 100);
    REQUIRE(atStart.start == 0);
    // Base 0 cannot sit on slot 1 at the very beginning; the edge fallback is
    // the right answer, so the window must report itself as serving rather than
    // chase a placement that does not exist.
    REQUIRE(windowServes(atStart, 0, cubic, 100));

    const auto atEnd = placeTimeWindow(99, cubic, 5, 100);
    REQUIRE(atEnd.start == 95);
    REQUIRE(windowServes(atEnd, 99, cubic, 100));
    REQUIRE(windowServes(atEnd, 98, cubic, 100));
  }

  SUBCASE("placement is idempotent under repeated advance") {
    auto w = placeTimeWindow(40, cubic, 8, 100);
    const auto again = placeTimeWindow(40, cubic, 8, 100);
    REQUIRE(w == again);
  }

  SUBCASE("monotone time reloads once per (count - width) cells") {
    const std::size_t count = 8;
    const std::size_t slices = 200;
    auto w = placeTimeWindow(0, cubic, count, slices);
    std::size_t reloads = 0;
    for (std::int64_t base = 0; base < 150; ++base) {
      if (!windowServes(w, base, cubic, slices)) {
        w = placeTimeWindow(base, cubic, count, slices);
        ++reloads;
        REQUIRE(windowServes(w, base, cubic, slices));
      }
    }
    // A freshly placed window serves the query cell plus (count - width) more,
    // so it covers (count - width + 1) cells before it has to move again. The
    // sync interval is quoted as (count - width) * dt, i.e. one cell of slack —
    // deliberately, so the reload happens before expiry rather than exactly at
    // it.
    const std::size_t span = count - cubic.width + 1;
    REQUIRE(reloads <= 150 / span + 1);
    REQUIRE(reloads >= 150 / span - 1);
  }
}

} // namespace seissol::unit_test
