// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Parallel/Pin.h"

#include <async/as/Pin.h>
#include <initializer_list>
#include <optional>
#include <sched.h>

namespace seissol::unit_test {

TEST_CASE("Online Mask Parsing" * doctest::test_suite("parallel")) {
  using namespace seissol::parallel;
  SUBCASE("Single") {
    const std::string mask = "1";
    auto is = Pinning::parseOnlineCpuMask(mask, 3);
    auto should = std::deque<bool>{false, true, false};
    CHECK(is == should);
  }

  SUBCASE("Range") {
    const std::string mask = "1-2";
    auto is = Pinning::parseOnlineCpuMask(mask, 3);
    auto should = std::deque<bool>{false, true, true};
    CHECK(is == should);
  }

  SUBCASE("Two single") {
    const std::string mask = "0,3";
    auto is = Pinning::parseOnlineCpuMask(mask, 4);
    auto should = std::deque<bool>{true, false, false, true};
    CHECK(is == should);
  }

  SUBCASE("Two ranges") {
    const std::string mask = "0-1,3-4";
    auto is = Pinning::parseOnlineCpuMask(mask, 5);
    auto should = std::deque<bool>{true, true, false, true, true};
    CHECK(is == should);
  }

  SUBCASE("Single and range") {
    const std::string mask = "0,3-4";
    auto is = Pinning::parseOnlineCpuMask(mask, 5);
    auto should = std::deque<bool>{true, false, false, true, true};
    CHECK(is == should);
  }
}

#ifndef __APPLE__
TEST_CASE("First CPU of a mask" * doctest::test_suite("parallel")) {
  using namespace seissol::parallel;
  const auto maskOf = [](std::initializer_list<int> cpus) {
    async::as::CpuMask cpuMask{};
    CPU_ZERO(&cpuMask.set);
    for (const int cpu : cpus) {
      CPU_SET(cpu, &cpuMask.set);
    }
    return cpuMask;
  };

  CHECK(Pinning::firstCpu(maskOf({})) == std::nullopt);
  CHECK(Pinning::firstCpu(maskOf({0})) == std::optional<int>(0));
  CHECK(Pinning::firstCpu(maskOf({1})) == std::optional<int>(1));
  // ids above the number of CPUs in the mask, such as the last cores of a node
  CHECK(Pinning::firstCpu(maskOf({47})) == std::optional<int>(47));
  CHECK(Pinning::firstCpu(maskOf({46, 47})) == std::optional<int>(46));
  CHECK(Pinning::firstCpu(maskOf({3, 2})) == std::optional<int>(2));
  CHECK(Pinning::firstCpu(maskOf({CPU_SETSIZE - 1})) == std::optional<int>(CPU_SETSIZE - 1));
}
#endif // __APPLE__

} // namespace seissol::unit_test
