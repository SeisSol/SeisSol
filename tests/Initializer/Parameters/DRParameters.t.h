// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Initializer/Parameters/DRParameters.h"
#include "Initializer/Parameters/ParameterReader.h"

#include <set>
#include <yaml-cpp/yaml.h>

namespace seissol::unit_test {
using namespace seissol::initializer::parameters;

// ---------------------------------------------------------------------------
// FrictionLawType: all enum values have unique integer representations
// ---------------------------------------------------------------------------

TEST_CASE("FrictionLawType enum values are unique" * doctest::test_suite("initializer")) {
  std::set<unsigned int> values;
  auto insert = [&](FrictionLawType t) {
    auto [it, inserted] = values.insert(static_cast<unsigned int>(t));
    CHECK_MESSAGE(inserted, "Duplicate FrictionLawType value: ", static_cast<unsigned int>(t));
  };

  insert(FrictionLawType::NoFault);
  insert(FrictionLawType::LinearSlipWeakeningLegacy);
  insert(FrictionLawType::LinearSlipWeakening);
  insert(FrictionLawType::LinearSlipWeakeningBimaterial);
  insert(FrictionLawType::LinearSlipWeakeningTPApprox);
  insert(FrictionLawType::RateAndStateAgingLaw);
  insert(FrictionLawType::RateAndStateSlipLaw);
  insert(FrictionLawType::RateAndStateFastVelocityWeakening);
  insert(FrictionLawType::ImposedSlipRatesYoffe);
  insert(FrictionLawType::ImposedSlipRatesGaussian);
  insert(FrictionLawType::ImposedSlipRatesDelta);
  insert(FrictionLawType::RateAndStateSevereVelocityWeakening);
  insert(FrictionLawType::RateAndStateAgingNucleation);

  // 13 distinct values
  CHECK(values.size() == 13);
}

// ---------------------------------------------------------------------------
// OutputType enum values
// ---------------------------------------------------------------------------

TEST_CASE("OutputType enum values" * doctest::test_suite("initializer")) {
  CHECK(static_cast<int>(OutputType::None) == 0);
  CHECK(static_cast<int>(OutputType::AtPickpoint) == 3);
  CHECK(static_cast<int>(OutputType::Elementwise) == 4);
  CHECK(static_cast<int>(OutputType::AtPickpointAndElementwise) == 5);
}

// ---------------------------------------------------------------------------
// SlipRateOutputType enum values
// ---------------------------------------------------------------------------

TEST_CASE("SlipRateOutputType enum values" * doctest::test_suite("initializer")) {
  CHECK(static_cast<int>(SlipRateOutputType::VelocityDifference) == 0);
  CHECK(static_cast<int>(SlipRateOutputType::TractionsAndFailure) == 1);
}

// ---------------------------------------------------------------------------
// readDRParameters: defaults with empty section
// ---------------------------------------------------------------------------

TEST_CASE("readDRParameters defaults" * doctest::test_suite("initializer")) {
  const YAML::Node node = YAML::Load(R"(
    dynamicrupture: {}
  )");
  ParameterReader reader(node, "", false);
  auto params = readDRParameters(&reader);

  CHECK(params.frictionLawType == FrictionLawType::NoFault);
  CHECK(params.outputPointType == OutputType::None);
  CHECK(params.refPointMethod == RefPointMethod::Point);
  CHECK(params.slipRateOutputType == SlipRateOutputType::TractionsAndFailure);
  CHECK(params.isThermalPressureOn == false);
  CHECK(params.referencePoint[0] == doctest::Approx(0.0));
  CHECK(params.referencePoint[1] == doctest::Approx(0.0));
  CHECK(params.referencePoint[2] == doctest::Approx(0.0));
  CHECK(params.healingThreshold == doctest::Approx(-1.0));
  CHECK(params.nucleationCount == 1);
}

// ---------------------------------------------------------------------------
// readDRParameters: specific friction law
// ---------------------------------------------------------------------------

TEST_CASE("readDRParameters with linear slip weakening" * doctest::test_suite("initializer")) {
  const YAML::Node node = YAML::Load(R"(
    dynamicrupture:
      fl: 16
      xref: 1.0
      yref: 2.0
      zref: 3.0
  )");
  ParameterReader reader(node, "", false);
  auto params = readDRParameters(&reader);

  CHECK(params.frictionLawType == FrictionLawType::LinearSlipWeakening);
  CHECK(params.referencePoint[0] == doctest::Approx(1.0));
  CHECK(params.referencePoint[1] == doctest::Approx(2.0));
  CHECK(params.referencePoint[2] == doctest::Approx(3.0));
}

// ---------------------------------------------------------------------------
// readDRParameters: imposed slip rate → forces SlipRateOutputType correction
// ---------------------------------------------------------------------------

TEST_CASE("readDRParameters imposed slip rates forces output type" *
          doctest::test_suite("initializer")) {
  // ImposedSlipRatesYoffe = 33; default slipRateOutputType = TractionsAndFailure(1)
  // The code should auto-correct this to VelocityDifference(0)
  const YAML::Node node = YAML::Load(R"(
    dynamicrupture:
      fl: 33
  )");
  ParameterReader reader(node, "", false);
  auto params = readDRParameters(&reader);

  CHECK(params.frictionLawType == FrictionLawType::ImposedSlipRatesYoffe);
  CHECK(params.slipRateOutputType == SlipRateOutputType::VelocityDifference);
}

TEST_CASE("readDRParameters Gaussian imposed also corrects" * doctest::test_suite("initializer")) {
  const YAML::Node node = YAML::Load(R"(
    dynamicrupture:
      fl: 34
  )");
  ParameterReader reader(node, "", false);
  auto params = readDRParameters(&reader);

  CHECK(params.frictionLawType == FrictionLawType::ImposedSlipRatesGaussian);
  CHECK(params.slipRateOutputType == SlipRateOutputType::VelocityDifference);
}

// ---------------------------------------------------------------------------
// readDRParameters: rate and state requires extra parameters
// ---------------------------------------------------------------------------

TEST_CASE("readDRParameters rate and state" * doctest::test_suite("initializer")) {
  const YAML::Node node = YAML::Load(R"(
    dynamicrupture:
      fl: 3
      rs_f0: 0.6
      rs_b: 0.012
      rs_sr0: 1.0e-6
      rs_inisliprate1: 1.0e-16
      rs_inisliprate2: 0.0
  )");
  ParameterReader reader(node, "", false);
  auto params = readDRParameters(&reader);

  CHECK(params.frictionLawType == FrictionLawType::RateAndStateAgingLaw);
  CHECK(params.rsF0 == doctest::Approx(0.6));
  CHECK(params.rsB == doctest::Approx(0.012));
  CHECK(params.rsSr0 == doctest::Approx(1.0e-6));
}

} // namespace seissol::unit_test
