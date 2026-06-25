// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Kernels/Common.h"
#include "Solver/Estimator.h"

#include <cmath>
#include <limits>
#include <variant>

namespace seissol::unit_test {

// ---------------------------------------------------------------------------
// Mini SeisSol
// ---------------------------------------------------------------------------

TEST_CASE("Run mini SeisSol" * doctest::test_suite("solver")) {
  // only check if it runs in a reasonable time (cf. SeisSol proxy)
  const auto time = seissol::solver::miniSeisSol();

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

} // namespace seissol::unit_test
