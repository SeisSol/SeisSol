// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/ConfigRegistry.h"
#include "Kernels/Common.h"
#include "Solver/Estimator.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <variant>
#include <vector>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::unit_test {

TEST_CASE("Cost factors of the configurations from their FLOPs" * doctest::test_suite("solver")) {
#ifdef ACL_DEVICE
  // the proxy runs on the device, as in a run
  auto& device = ::device::DeviceInstance::instance();
  device.api().setDevice(0);
  device.api().initialize();
#endif

  const auto count = builtConfigCount();
  std::vector<ConfigId> configs(count);
  std::iota(configs.begin(), configs.end(), ConfigId{0});
  const auto last = static_cast<ConfigId>(count - 1);

  // every configuration of the run gets its factor; the reference costs 1
  const auto factors = solver::configCostFactors(configs, 0, false);
  REQUIRE(factors.size() == count);
  CHECK(factors[0] == 1.0);
  for (std::size_t config = 0; config < count; ++config) {
    CAPTURE(config);
    CHECK(std::isfinite(factors[config]));
    CHECK(factors[config] > 0.0);
  }

  // the factors only depend on the reference through its cost
  const auto factorsOfLast = solver::configCostFactors(configs, last, false);
  CHECK(factorsOfLast[last] == 1.0);
  for (std::size_t config = 0; config < count; ++config) {
    CAPTURE(config);
    CHECK(factorsOfLast[config] * factors[last] == doctest::Approx(factors[config]));
  }

  // the configurations that are not in the run cost as much as the reference
  CHECK(solver::configCostFactors({last}, last, false) == std::vector<double>(count, 1.0));

#ifdef ACL_DEVICE
  // as in main(); otherwise, the device is torn down only by static destructors at exit
  device.api().finalize();
#endif
}

/*

// disabled for now, due to time constraints

// ---------------------------------------------------------------------------
// Mini SeisSol
// ---------------------------------------------------------------------------

TEST_CASE("Run mini SeisSol" * doctest::test_suite("solver")) {
  // only check if it runs in a reasonable time (cf. SeisSol proxy)
  const auto time = seissol::solver::miniSeisSol(defaultConfig());

  // let it take less than 10000 s
  CHECK(time < 10000);
}

// ---------------------------------------------------------------------------
// Host-device switch (dummy)
// ---------------------------------------------------------------------------

TEST_CASE("Host-device switch" * doctest::test_suite("solver")) {
  const auto switchpoint = seissol::solver::hostDeviceSwitch();
  if constexpr (isDeviceOn()) {
    // alas, we can only "run" it here and check that the result is reasonable
    // check that it is smaller than 2**21 (the range that is check right now)
    CHECK(switchpoint < 2097152);
  } else {
    // disabled; should always return zero
    CHECK(switchpoint == 0);
  }
}

*/

} // namespace seissol::unit_test
