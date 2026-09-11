// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "SourceTerm/FSRMReader.h"
#include "TestHelper.h"

#include <cstdlib>

namespace seissol::unit_test {
TEST_CASE("FSRM Reader 1" * doctest::test_suite("reader")) {
  seissol::sourceterm::FSRMSource fsrm;
  fsrm.read(tpath("Testing/fsrm_source1.dat"));

  CHECK(fsrm.momentTensor[0][0] == doctest::Approx(0.0));
  CHECK(fsrm.momentTensor[0][1] == doctest::Approx(-0.1e+19));
  CHECK(fsrm.momentTensor[0][2] == doctest::Approx(0.0));
  CHECK(fsrm.momentTensor[1][0] == doctest::Approx(-0.1e+19));
  CHECK(fsrm.momentTensor[1][1] == doctest::Approx(0.0));
  CHECK(fsrm.momentTensor[1][2] == doctest::Approx(0.0));
  CHECK(fsrm.momentTensor[2][0] == doctest::Approx(0.0));
  CHECK(fsrm.momentTensor[2][1] == doctest::Approx(0.0));
  CHECK(fsrm.momentTensor[2][2] == doctest::Approx(0.0));

  CHECK(fsrm.numberOfSources == 1);

  CHECK(fsrm.solidVelocityComponent[0] == doctest::Approx(0.0));
  CHECK(fsrm.solidVelocityComponent[1] == doctest::Approx(0.0));
  CHECK(fsrm.solidVelocityComponent[2] == doctest::Approx(0.0));
  CHECK(fsrm.pressureComponent == doctest::Approx(0.0));
  CHECK(fsrm.fluidVelocityComponent[0] == doctest::Approx(0.0));
  CHECK(fsrm.fluidVelocityComponent[1] == doctest::Approx(0.0));
  CHECK(fsrm.fluidVelocityComponent[2] == doctest::Approx(0.0));

  CHECK(fsrm.centers[0](0) == doctest::Approx(0.0));
  CHECK(fsrm.centers[0](1) == doctest::Approx(0.0));
  CHECK(fsrm.centers[0](2) == doctest::Approx(2000.0));
  CHECK(fsrm.strikes[0] == doctest::Approx(0.0));
  CHECK(fsrm.dips[0] == doctest::Approx(0.0));
  CHECK(fsrm.rakes[0] == doctest::Approx(0.0));
  CHECK(fsrm.areas[0] == doctest::Approx(1.0));
  CHECK(fsrm.onsets[0] == doctest::Approx(0.0));

  CHECK(fsrm.timestep == doctest::Approx(0.02));
  CHECK(fsrm.numberOfSamples == 200);

  CHECK(fsrm.timeHistories.size() == fsrm.numberOfSources);
  CHECK(fsrm.timeHistories[0].size() == fsrm.numberOfSamples);

  // (we won't check the individual samples here)
}

TEST_CASE("FSRM Reader 2" * doctest::test_suite("reader")) {
  seissol::sourceterm::FSRMSource fsrm;
  fsrm.read(tpath("Testing/fsrm_source2.dat"));

  CHECK(fsrm.momentTensor[0][0] == doctest::Approx(0.01));
  CHECK(fsrm.momentTensor[0][1] == doctest::Approx(0.02));
  CHECK(fsrm.momentTensor[0][2] == doctest::Approx(0.03));
  CHECK(fsrm.momentTensor[1][0] == doctest::Approx(0.04));
  CHECK(fsrm.momentTensor[1][1] == doctest::Approx(0.05));
  CHECK(fsrm.momentTensor[1][2] == doctest::Approx(0.06));
  CHECK(fsrm.momentTensor[2][0] == doctest::Approx(0.07));
  CHECK(fsrm.momentTensor[2][1] == doctest::Approx(0.08));
  CHECK(fsrm.momentTensor[2][2] == doctest::Approx(0.09));

  CHECK(fsrm.numberOfSources == 3);

  CHECK(fsrm.solidVelocityComponent[0] == doctest::Approx(1.01));
  CHECK(fsrm.solidVelocityComponent[1] == doctest::Approx(1.02));
  CHECK(fsrm.solidVelocityComponent[2] == doctest::Approx(1.03));
  CHECK(fsrm.pressureComponent == doctest::Approx(2.01));
  CHECK(fsrm.fluidVelocityComponent[0] == doctest::Approx(3.01));
  CHECK(fsrm.fluidVelocityComponent[1] == doctest::Approx(3.02));
  CHECK(fsrm.fluidVelocityComponent[2] == doctest::Approx(3.03));

  CHECK(fsrm.centers[0](0) == doctest::Approx(0.1));
  CHECK(fsrm.centers[0](1) == doctest::Approx(0.2));
  CHECK(fsrm.centers[0](2) == doctest::Approx(0.3));
  CHECK(fsrm.strikes[0] == doctest::Approx(0.4));
  CHECK(fsrm.dips[0] == doctest::Approx(0.5));
  CHECK(fsrm.rakes[0] == doctest::Approx(0.6));
  CHECK(fsrm.areas[0] == doctest::Approx(0.7));
  CHECK(fsrm.onsets[0] == doctest::Approx(0.8));
  CHECK(fsrm.centers[1](0) == doctest::Approx(1.1));
  CHECK(fsrm.centers[1](1) == doctest::Approx(1.2));
  CHECK(fsrm.centers[1](2) == doctest::Approx(1.3));
  CHECK(fsrm.strikes[1] == doctest::Approx(1.4));
  CHECK(fsrm.dips[1] == doctest::Approx(1.5));
  CHECK(fsrm.rakes[1] == doctest::Approx(1.6));
  CHECK(fsrm.areas[1] == doctest::Approx(1.7));
  CHECK(fsrm.onsets[1] == doctest::Approx(1.8));
  CHECK(fsrm.centers[2](0) == doctest::Approx(2.1));
  CHECK(fsrm.centers[2](1) == doctest::Approx(2.2));
  CHECK(fsrm.centers[2](2) == doctest::Approx(2.3));
  CHECK(fsrm.strikes[2] == doctest::Approx(2.4));
  CHECK(fsrm.dips[2] == doctest::Approx(2.5));
  CHECK(fsrm.rakes[2] == doctest::Approx(2.6));
  CHECK(fsrm.areas[2] == doctest::Approx(2.7));
  CHECK(fsrm.onsets[2] == doctest::Approx(2.8));

  CHECK(fsrm.timestep == doctest::Approx(0.9));
  CHECK(fsrm.numberOfSamples == 4);

  CHECK(fsrm.timeHistories[0][0] == doctest::Approx(0.0));
  CHECK(fsrm.timeHistories[0][1] == doctest::Approx(1.0));
  CHECK(fsrm.timeHistories[0][2] == doctest::Approx(2.0));
  CHECK(fsrm.timeHistories[0][3] == doctest::Approx(3.0));
  CHECK(fsrm.timeHistories[1][0] == doctest::Approx(4.0));
  CHECK(fsrm.timeHistories[1][1] == doctest::Approx(5.0));
  CHECK(fsrm.timeHistories[1][2] == doctest::Approx(6.0));
  CHECK(fsrm.timeHistories[1][3] == doctest::Approx(7.0));
  CHECK(fsrm.timeHistories[2][0] == doctest::Approx(8.0));
  CHECK(fsrm.timeHistories[2][1] == doctest::Approx(9.0));
  CHECK(fsrm.timeHistories[2][2] == doctest::Approx(10.0));
  CHECK(fsrm.timeHistories[2][3] == doctest::Approx(11.0));

  CHECK(fsrm.timeHistories.size() == fsrm.numberOfSources);
  CHECK(fsrm.timeHistories[0].size() == fsrm.numberOfSamples);
  CHECK(fsrm.timeHistories[1].size() == fsrm.numberOfSamples);
  CHECK(fsrm.timeHistories[2].size() == fsrm.numberOfSamples);
}

TEST_CASE("FSRM Reader File Not Found" * doctest::test_suite("reader")) {
  seissol::sourceterm::FSRMSource fsrm;
  REQUIRE_THROWS(fsrm.read(tpath("Testing/unknown-file.dat")));
}
} // namespace seissol::unit_test
