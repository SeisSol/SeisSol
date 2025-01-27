// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "tests/TestHelper.h"
#include <cstdlib>

#include "SourceTerm/FSRMReader.h"

namespace seissol::unit_test {
TEST_CASE("FSRM Reader 1") {
  seissol::sourceterm::FSRMSource fsrm;
  fsrm.read(std::string("Testing/fsrm_source1.dat"));

  REQUIRE(fsrm.momentTensor[0][0] == doctest::Approx(0.0));
  REQUIRE(fsrm.momentTensor[0][1] == doctest::Approx(-0.1e+19));
  REQUIRE(fsrm.momentTensor[0][2] == doctest::Approx(0.0));
  REQUIRE(fsrm.momentTensor[1][0] == doctest::Approx(-0.1e+19));
  REQUIRE(fsrm.momentTensor[1][1] == doctest::Approx(0.0));
  REQUIRE(fsrm.momentTensor[1][2] == doctest::Approx(0.0));
  REQUIRE(fsrm.momentTensor[2][0] == doctest::Approx(0.0));
  REQUIRE(fsrm.momentTensor[2][1] == doctest::Approx(0.0));
  REQUIRE(fsrm.momentTensor[2][2] == doctest::Approx(0.0));

  REQUIRE(fsrm.numberOfSources == 1);

  REQUIRE(fsrm.solidVelocityComponent[0] == doctest::Approx(0.0));
  REQUIRE(fsrm.solidVelocityComponent[1] == doctest::Approx(0.0));
  REQUIRE(fsrm.solidVelocityComponent[2] == doctest::Approx(0.0));
  REQUIRE(fsrm.pressureComponent == doctest::Approx(0.0));
  REQUIRE(fsrm.fluidVelocityComponent[0] == doctest::Approx(0.0));
  REQUIRE(fsrm.fluidVelocityComponent[1] == doctest::Approx(0.0));
  REQUIRE(fsrm.fluidVelocityComponent[2] == doctest::Approx(0.0));

  REQUIRE(fsrm.centers[0](0) == doctest::Approx(0.0));
  REQUIRE(fsrm.centers[0](1) == doctest::Approx(0.0));
  REQUIRE(fsrm.centers[0](2) == doctest::Approx(2000.0));
  REQUIRE(fsrm.strikes[0] == doctest::Approx(0.0));
  REQUIRE(fsrm.dips[0] == doctest::Approx(0.0));
  REQUIRE(fsrm.rakes[0] == doctest::Approx(0.0));
  REQUIRE(fsrm.areas[0] == doctest::Approx(1.0));
  REQUIRE(fsrm.onsets[0] == doctest::Approx(0.0));

  REQUIRE(fsrm.timestep == doctest::Approx(0.02));
  REQUIRE(fsrm.numberOfSamples == 200);

  REQUIRE(fsrm.timeHistories.size() == fsrm.numberOfSources);
  REQUIRE(fsrm.timeHistories[0].size() == fsrm.numberOfSamples);

  // (we won't check the individual samples here)
}

TEST_CASE("FSRM Reader 2") {
  seissol::sourceterm::FSRMSource fsrm;
  fsrm.read(std::string("Testing/fsrm_source2.dat"));

  REQUIRE(fsrm.momentTensor[0][0] == doctest::Approx(0.01));
  REQUIRE(fsrm.momentTensor[0][1] == doctest::Approx(0.02));
  REQUIRE(fsrm.momentTensor[0][2] == doctest::Approx(0.03));
  REQUIRE(fsrm.momentTensor[1][0] == doctest::Approx(0.04));
  REQUIRE(fsrm.momentTensor[1][1] == doctest::Approx(0.05));
  REQUIRE(fsrm.momentTensor[1][2] == doctest::Approx(0.06));
  REQUIRE(fsrm.momentTensor[2][0] == doctest::Approx(0.07));
  REQUIRE(fsrm.momentTensor[2][1] == doctest::Approx(0.08));
  REQUIRE(fsrm.momentTensor[2][2] == doctest::Approx(0.09));

  REQUIRE(fsrm.numberOfSources == 3);

  REQUIRE(fsrm.solidVelocityComponent[0] == doctest::Approx(1.01));
  REQUIRE(fsrm.solidVelocityComponent[1] == doctest::Approx(1.02));
  REQUIRE(fsrm.solidVelocityComponent[2] == doctest::Approx(1.03));
  REQUIRE(fsrm.pressureComponent == doctest::Approx(2.01));
  REQUIRE(fsrm.fluidVelocityComponent[0] == doctest::Approx(3.01));
  REQUIRE(fsrm.fluidVelocityComponent[1] == doctest::Approx(3.02));
  REQUIRE(fsrm.fluidVelocityComponent[2] == doctest::Approx(3.03));

  REQUIRE(fsrm.centers[0](0) == doctest::Approx(0.1));
  REQUIRE(fsrm.centers[0](1) == doctest::Approx(0.2));
  REQUIRE(fsrm.centers[0](2) == doctest::Approx(0.3));
  REQUIRE(fsrm.strikes[0] == doctest::Approx(0.4));
  REQUIRE(fsrm.dips[0] == doctest::Approx(0.5));
  REQUIRE(fsrm.rakes[0] == doctest::Approx(0.6));
  REQUIRE(fsrm.areas[0] == doctest::Approx(0.7));
  REQUIRE(fsrm.onsets[0] == doctest::Approx(0.8));
  REQUIRE(fsrm.centers[1](0) == doctest::Approx(1.1));
  REQUIRE(fsrm.centers[1](1) == doctest::Approx(1.2));
  REQUIRE(fsrm.centers[1](2) == doctest::Approx(1.3));
  REQUIRE(fsrm.strikes[1] == doctest::Approx(1.4));
  REQUIRE(fsrm.dips[1] == doctest::Approx(1.5));
  REQUIRE(fsrm.rakes[1] == doctest::Approx(1.6));
  REQUIRE(fsrm.areas[1] == doctest::Approx(1.7));
  REQUIRE(fsrm.onsets[1] == doctest::Approx(1.8));
  REQUIRE(fsrm.centers[2](0) == doctest::Approx(2.1));
  REQUIRE(fsrm.centers[2](1) == doctest::Approx(2.2));
  REQUIRE(fsrm.centers[2](2) == doctest::Approx(2.3));
  REQUIRE(fsrm.strikes[2] == doctest::Approx(2.4));
  REQUIRE(fsrm.dips[2] == doctest::Approx(2.5));
  REQUIRE(fsrm.rakes[2] == doctest::Approx(2.6));
  REQUIRE(fsrm.areas[2] == doctest::Approx(2.7));
  REQUIRE(fsrm.onsets[2] == doctest::Approx(2.8));

  REQUIRE(fsrm.timestep == doctest::Approx(0.9));
  REQUIRE(fsrm.numberOfSamples == 4);

  REQUIRE(fsrm.timeHistories[0][0] == doctest::Approx(0.0));
  REQUIRE(fsrm.timeHistories[0][1] == doctest::Approx(1.0));
  REQUIRE(fsrm.timeHistories[0][2] == doctest::Approx(2.0));
  REQUIRE(fsrm.timeHistories[0][3] == doctest::Approx(3.0));
  REQUIRE(fsrm.timeHistories[1][0] == doctest::Approx(4.0));
  REQUIRE(fsrm.timeHistories[1][1] == doctest::Approx(5.0));
  REQUIRE(fsrm.timeHistories[1][2] == doctest::Approx(6.0));
  REQUIRE(fsrm.timeHistories[1][3] == doctest::Approx(7.0));
  REQUIRE(fsrm.timeHistories[2][0] == doctest::Approx(8.0));
  REQUIRE(fsrm.timeHistories[2][1] == doctest::Approx(9.0));
  REQUIRE(fsrm.timeHistories[2][2] == doctest::Approx(10.0));
  REQUIRE(fsrm.timeHistories[2][3] == doctest::Approx(11.0));

  REQUIRE(fsrm.timeHistories.size() == fsrm.numberOfSources);
  REQUIRE(fsrm.timeHistories[0].size() == fsrm.numberOfSamples);
  REQUIRE(fsrm.timeHistories[1].size() == fsrm.numberOfSamples);
  REQUIRE(fsrm.timeHistories[2].size() == fsrm.numberOfSamples);
}

TEST_CASE("FSRM Reader File Not Found") {
  seissol::sourceterm::FSRMSource fsrm;
  REQUIRE_THROWS(fsrm.read(std::string("Testing/unknown-file.dat")));
}
} // namespace seissol::unit_test
