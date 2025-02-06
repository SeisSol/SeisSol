// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Parallel/Pin.h"

namespace seissol::unit_test {

TEST_CASE("Online Mask Parsing") {
  using namespace seissol::parallel;
  SUBCASE("Single") {
    const std::string mask = "1";
    auto is = Pinning::parseOnlineCpuMask(mask, 3);
    auto should = std::deque<bool>{false, true, false};
    REQUIRE(is == should);
  }

  SUBCASE("Range") {
    const std::string mask = "1-2";
    auto is = Pinning::parseOnlineCpuMask(mask, 3);
    auto should = std::deque<bool>{false, true, true};
    REQUIRE(is == should);
  }

  SUBCASE("Two single") {
    const std::string mask = "0,3";
    auto is = Pinning::parseOnlineCpuMask(mask, 4);
    auto should = std::deque<bool>{true, false, false, true};
    REQUIRE(is == should);
  }

  SUBCASE("Two ranges") {
    const std::string mask = "0-1,3-4";
    auto is = Pinning::parseOnlineCpuMask(mask, 5);
    auto should = std::deque<bool>{true, true, false, true, true};
    REQUIRE(is == should);
  }

  SUBCASE("Single and range") {
    const std::string mask = "0,3-4";
    auto is = Pinning::parseOnlineCpuMask(mask, 5);
    auto should = std::deque<bool>{true, false, false, true, true};
    REQUIRE(is == should);
  }
}

} // namespace seissol::unit_test
