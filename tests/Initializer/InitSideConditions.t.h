// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Initializer/InitProcedure/InitSideConditions.h"
#include "Initializer/Parameters/InitializationParameters.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Initializer/Parameters/ParameterReader.h"
#include "Initializer/Typedefs.h"

#include <yaml-cpp/yaml.h>

namespace seissol::unit_test {

TEST_CASE("The acoustic travelling wave passes the time mirror only if it is enabled" *
          doctest::test_suite("initializer")) {
  using namespace seissol::initializer::parameters;
  using seissol::initializer::initprocedure::getAcousticTravellingWaveITMInformation;

  // mirror settings as in docs/parameters.par
  YAML::Node node = YAML::Load(R"(
    equations:
      itmenable: 1
      itmstartingtime: 2.0
      itmtime: 0.01
      itmvelocityscalingfactor: 2
    inicondition:
      cictype: AcousticTravellingwithITM
      k: 6.283
  )");
  const auto waveFor = [](const YAML::Node& parameterFile) {
    ParameterReader reader(parameterFile, "", false);
    const ITMParameters itm = readITMParameters(&reader);
    const InitializationParameters initialization = readInitializationParameters(&reader);
    return getAcousticTravellingWaveITMInformation(initialization, itm);
  };

  SUBCASE("enabled") {
    const AcousticTravellingWaveParametersITM wave = waveFor(node);
    CHECK(wave.k == doctest::Approx(6.283));
    CHECK(wave.itmStartingTime == doctest::Approx(2.0));
    CHECK(wave.itmDuration == doctest::Approx(0.01));
    CHECK(wave.itmVelocityScalingFactor == doctest::Approx(2.0));
  }

  SUBCASE("disabled") {
    // The file keeps its mirror settings, but they are ignored now. The reference solution must
    // not model a mirror the run does not have.
    node["equations"]["itmenable"] = 0;
    const AcousticTravellingWaveParametersITM wave = waveFor(node);
    CHECK(wave.k == doctest::Approx(6.283));
    CHECK(wave.itmVelocityScalingFactor == doctest::Approx(1.0));
  }
}

} // namespace seissol::unit_test
