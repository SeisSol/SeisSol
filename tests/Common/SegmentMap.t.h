// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Common/SegmentMap.h"

namespace seissol::unit_test {

using namespace seissol;

TEST_CASE("SegmentMap" * doctest::test_suite("common")) {
  SegmentMap<int32_t, int32_t> map;

  map.addRange(0, 10, 123);
  map.addRange(14, 14, 23);
  map.addRange(23, 100, 543);
  map.addRange(1000, {}, -9);
  map.addRange({}, -10, -100);

  CHECK(map.at(0) == 123);
  CHECK(map.at(1) == 123);
  CHECK(map.at(9) == 123);
  CHECK(map.at(10) == 123);
  CHECK_FALSE(map.at(11).has_value());
  CHECK_FALSE(map.at(12).has_value());
  CHECK_FALSE(map.at(13).has_value());
  CHECK(map.at(14) == 23);
  CHECK_FALSE(map.at(15).has_value());
  CHECK_FALSE(map.at(16).has_value());
  CHECK_FALSE(map.at(22).has_value());
  CHECK(map.at(23) == 543);
  CHECK(map.at(24) == 543);
  CHECK(map.at(99) == 543);
  CHECK(map.at(100) == 543);
  CHECK_FALSE(map.at(101).has_value());
  CHECK_FALSE(map.at(999).has_value());
  CHECK(map.at(1000) == -9);
  CHECK(map.at(1000000) == -9);

  CHECK_FALSE(map.at(-1).has_value());
  CHECK_FALSE(map.at(-9).has_value());
  CHECK(map.at(-10) == -100);
  CHECK(map.at(-11) == -100);
  CHECK(map.at(-1000000) == -100);
}

} // namespace seissol::unit_test
