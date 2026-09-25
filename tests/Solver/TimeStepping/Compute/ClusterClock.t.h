// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Parallel/Runtime/Stream.h"
#include "Solver/TimeStepping/Compute/ClusterClock.h"

namespace seissol::unit_test {

TEST_CASE("The cluster clock adds up the steps like the cluster time" *
          doctest::test_suite("solver")) {
  parallel::runtime::StreamRuntime runtime;
  time_stepping::ClusterClock clock;
  CHECK(*clock.host() == 0);

  clock.set(0.25, runtime);
  CHECK(*clock.host() == 0.25);

  // the same additions as the time of a cluster, including a shortened last step
  double time = 0.25;
  for (const double timeStepSize : {0.1, 0.1, 0.1, 0.037, 1.0e-3, 0.1}) {
    clock.advance(timeStepSize, runtime);
    time += timeStepSize;
    CHECK(*clock.host() == time);
  }

  clock.set(1.0 / 3.0, runtime);
  CHECK(*clock.host() == 1.0 / 3.0);

  // without a device, there is no device copy
  CHECK(clock.device() == nullptr);
  clock.dispose();
}

} // namespace seissol::unit_test
