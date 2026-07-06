// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Equations/Datastructures.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Physics/InstantaneousTimeMirrorManager.h"
#include "TestHelper.h"

#include <type_traits>

namespace seissol::unit_test {

TEST_CASE("Anisotropic Instantaneous Time Mirror supports only BothWaves reflection type" *
          doctest::skip(!std::is_same_v<model::MaterialT, model::AnisotropicMaterial>)) {
  using seissol::initializer::parameters::ReflectionType;

  CHECK(seissol::ITM::isAnisotropicReflectionTypeSupported(ReflectionType::BothWaves));
  CHECK_FALSE(seissol::ITM::isAnisotropicReflectionTypeSupported(ReflectionType::Pwave));
  CHECK_FALSE(seissol::ITM::isAnisotropicReflectionTypeSupported(ReflectionType::Swave));
  CHECK_FALSE(
      seissol::ITM::isAnisotropicReflectionTypeSupported(ReflectionType::BothWavesVelocity));
}

TEST_CASE("Instantaneous Time Mirror Swave lambda scaling matches implementation formula") {
  const double lambda = 8.0;
  const double mu = 3.0;
  const double velocityScalingFactor = 2.0;
  const double expected =
      (lambda + 2.0 * mu) / velocityScalingFactor - 2.0 * velocityScalingFactor * mu;

  CHECK(seissol::ITM::getSwaveScaledLambda(lambda, mu, velocityScalingFactor) ==
        AbsApprox(expected).epsilon(1.0e-15));
}

TEST_CASE("Elastic Instantaneous Time Mirror time-step scaling depends on reflection type") {
  using seissol::initializer::parameters::ReflectionType;

  const double velocityScalingFactor = 2.5;
  CHECK(seissol::ITM::getElasticTimeStepScalingFactor(ReflectionType::BothWaves,
                                                      velocityScalingFactor) ==
        AbsApprox(1.0 / velocityScalingFactor).epsilon(1.0e-15));
  CHECK(
      seissol::ITM::getElasticTimeStepScalingFactor(ReflectionType::Pwave, velocityScalingFactor) ==
      AbsApprox(1.0 / velocityScalingFactor).epsilon(1.0e-15));
  CHECK(
      seissol::ITM::getElasticTimeStepScalingFactor(ReflectionType::Swave, velocityScalingFactor) ==
      AbsApprox(velocityScalingFactor).epsilon(1.0e-15));
  CHECK(seissol::ITM::getElasticTimeStepScalingFactor(ReflectionType::BothWavesVelocity,
                                                      velocityScalingFactor) ==
        AbsApprox(1.0).epsilon(1.0e-15));
}
} // namespace seissol::unit_test
