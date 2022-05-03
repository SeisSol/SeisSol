#ifndef SEISSOL_PINTEST_T_H
#define SEISSOL_PINTEST_T_H

#include "Parallel/Pin.h"

namespace seissol::unit_test {

TEST_CASE("Online Mask Parsing") {
  using namespace seissol::parallel;
  SUBCASE("Single") {
    auto mask = "1";
    auto is = Pinning::parseOnlineCpuMask(mask, 3);
    auto should = std::deque<bool>{false, true, false};
    REQUIRE(is == should);
  }

  SUBCASE("Range") {
    auto mask = "1-2";
    auto is = Pinning::parseOnlineCpuMask(mask, 3);
    auto should = std::deque<bool>{false, true, true};
    REQUIRE(is == should);
  }

  SUBCASE("Two single") {
    auto mask = "0,3";
    auto is = Pinning::parseOnlineCpuMask(mask, 4);
    auto should = std::deque<bool>{true, false, false, true};
    REQUIRE(is == should);
  }

  SUBCASE("Two ranges") {
    auto mask = "0-1,3-4";
    auto is = Pinning::parseOnlineCpuMask(mask, 5);
    auto should = std::deque<bool>{true, true, false, true, true};
    REQUIRE(is == should);
  }

  SUBCASE("Single and range") {
    auto mask = "0,3-4";
    auto is = Pinning::parseOnlineCpuMask(mask, 5);
    auto should = std::deque<bool>{true, false, false, true, true};
    REQUIRE(is == should);
  }
}

} // namespace seissol::unit_test
#endif // SEISSOL_PINTEST_T_H
