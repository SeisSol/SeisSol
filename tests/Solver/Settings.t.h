// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Solver/Settings.h"

namespace seissol::unit_test {

TEST_CASE("Copies of the cells for the simulations" * doctest::test_suite("solver")) {
  // the fused simulations hold all simulations
  CHECK(SimulationSettings(false, false, 1, 1).copyCount() == 1);
  CHECK(SimulationSettings(false, false, 8, 8).copyCount() == 1);
  CHECK(SimulationSettings(false, false, 5, 8).copyCount() == 1);

  // one copy per simulation without fused simulations
  CHECK(SimulationSettings(false, false, 8, 1).copyCount() == 8);

  // the last copy may be filled up only partially
  CHECK(SimulationSettings(false, false, 8, 4).copyCount() == 2);
  CHECK(SimulationSettings(false, false, 9, 4).copyCount() == 3);
}

} // namespace seissol::unit_test
