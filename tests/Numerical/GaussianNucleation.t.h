// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Numerical/GaussianNucleationFunction.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>

namespace seissol::unit_test {
namespace gaussiannucleation {

using seissol::gaussianNucleationFunction::smoothStep;
using seissol::gaussianNucleationFunction::smoothStepIncrement;

struct Sample {
  double t0;
  double dt;
  double currentTime;
  double increment;
  double singleIncrement;
};

/**
 * f(t) - f(t - dt) at 60 decimal digits, rounded once to the nearest double. The samples walk
 * the ramp for three nucleation times and three step sizes and include the step that crosses
 * the end of the ramp and the first step after it.
 *
 * t - dt is not representable, and an error there reappears divided by dt, so the reference
 * uses the rounded previous time. singleIncrement is the same quantity for arguments and
 * previous time rounded to single precision.
 */
constexpr std::size_t SampleCount = 54;
constexpr Sample Samples[SampleCount] = {
    {0x1.0000000000000p+0,
     0x1.999999999999ap-4,
     0x1.3333333333333p-2,
     0x1.b5691b271b877p-3,
     0x1.b56919709bd30p-3},
    {0x1.0000000000000p+0,
     0x1.999999999999ap-4,
     0x1.3333333333333p-1,
     0x1.c2b3251202670p-4,
     0x1.c2b32b116f2d4p-4},
    {0x1.0000000000000p+0,
     0x1.999999999999ap-4,
     0x1.ccccccccccccdp-1,
     0x1.f7fa5ef0427a9p-6,
     0x1.f7fa6f1517080p-6},
    {0x1.0000000000000p+0,
     0x1.999999999999ap-4,
     0x1.ff7ced916872bp-1,
     0x1.4ff1be23a076cp-7,
     0x1.4ff1c30030d5bp-7},
    {0x1.0000000000000p+0,
     0x1.999999999999ap-4,
     0x1.0000000000000p+0,
     0x1.4952e7a5a8517p-7,
     0x1.4952f1fd715d6p-7},
    {0x1.0000000000000p+0,
     0x1.999999999999ap-4,
     0x1.0cccccccccccdp+0,
     0x1.48170661903a8p-9,
     0x1.481743f995e03p-9},
    {0x1.0000000000000p+0,
     0x1.47ae147ae147bp-7,
     0x1.3333333333333p-2,
     0x1.53e875e8eb0cfp-6,
     0x1.53e860210bb69p-6},
    {0x1.0000000000000p+0,
     0x1.47ae147ae147bp-7,
     0x1.3333333333333p-1,
     0x1.38236d6636817p-7,
     0x1.382357fd1f878p-7},
    {0x1.0000000000000p+0,
     0x1.47ae147ae147bp-7,
     0x1.ccccccccccccdp-1,
     0x1.164f3c5f1629ep-9,
     0x1.164f2e5cca814p-9},
    {0x1.0000000000000p+0,
     0x1.47ae147ae147bp-7,
     0x1.ff7ced916872bp-1,
     0x1.f758e1054db3ap-14,
     0x1.f758608dc98a7p-14},
    {0x1.0000000000000p+0,
     0x1.47ae147ae147bp-7,
     0x1.0000000000000p+0,
     0x1.a3738d213712ep-14,
     0x1.a37358b21b421p-14},
    {0x1.0000000000000p+0,
     0x1.47ae147ae147bp-7,
     0x1.0147ae147ae14p+0,
     0x1.a36f864b6e114p-16,
     0x1.a36f51dd53e2ap-16},
    {0x1.0000000000000p+0,
     0x1.0624dd2f1a9fcp-10,
     0x1.3333333333333p-2,
     0x1.0e204f4009f34p-9,
     0x1.0e1f6ad5b85afp-9},
    {0x1.0000000000000p+0,
     0x1.0624dd2f1a9fcp-10,
     0x1.3333333333333p-1,
     0x1.ec23da9f3019dp-11,
     0x1.ec22382fb576ap-11},
    {0x1.0000000000000p+0,
     0x1.0624dd2f1a9fcp-10,
     0x1.ccccccccccccdp-1,
     0x1.a9ce7e83b7ab2p-13,
     0x1.a9cd1c2b9571ap-13},
    {0x1.0000000000000p+0,
     0x1.0624dd2f1a9fcp-10,
     0x1.ff7ced916872bp-1,
     0x1.92a7790993f80p-19,
     0x1.92a4d18fa70c4p-19},
    {0x1.0000000000000p+0,
     0x1.0624dd2f1a9fcp-10,
     0x1.0000000000000p+0,
     0x1.0c6f82d72bcb2p-20,
     0x1.0c6dbddbaf422p-20},
    {0x1.0000000000000p+0,
     0x1.0624dd2f1a9fcp-10,
     0x1.0020c49ba5e35p+0,
     0x1.0c6f7c3e524e2p-22,
     0x1.0c75e8730e23dp-22},
    {0x1.0000000000000p-1,
     0x1.999999999999ap-4,
     0x1.3333333333333p-3,
     0x1.795bf99e5f61ep-2,
     0x1.795bfad84dd1fp-2},
    {0x1.0000000000000p-1,
     0x1.999999999999ap-4,
     0x1.3333333333333p-2,
     0x1.06f205723b9d1p-2,
     0x1.06f202f07ae49p-2},
    {0x1.0000000000000p-1,
     0x1.999999999999ap-4,
     0x1.ccccccccccccdp-2,
     0x1.588ba2d789564p-4,
     0x1.588ba3a638affp-4},
    {0x1.0000000000000p-1,
     0x1.999999999999ap-4,
     0x1.ff7ced916872bp-2,
     0x1.51bb3f795c051p-5,
     0x1.51bb39ea7c6f4p-5},
    {0x1.0000000000000p-1,
     0x1.999999999999ap-4,
     0x1.0000000000000p-1,
     0x1.4e51e9618b51ap-5,
     0x1.4e51e6b77436ap-5},
    {0x1.0000000000000p-1,
     0x1.999999999999ap-4,
     0x1.199999999999ap-1,
     0x1.4952e7a5a850ap-7,
     0x1.4952d821fb0f1p-7},
    {0x1.0000000000000p-1,
     0x1.47ae147ae147bp-7,
     0x1.3333333333333p-3,
     0x1.563bab238ff13p-5,
     0x1.563bb6f6e01cep-5},
    {0x1.0000000000000p-1,
     0x1.47ae147ae147bp-7,
     0x1.3333333333333p-2,
     0x1.3d4130efcc52dp-6,
     0x1.3d411addab9efp-6},
    {0x1.0000000000000p-1,
     0x1.47ae147ae147bp-7,
     0x1.ccccccccccccdp-2,
     0x1.23e5e14dd439cp-8,
     0x1.23e5d19acc7a3p-8},
    {0x1.0000000000000p-1,
     0x1.47ae147ae147bp-7,
     0x1.ff7ced916872bp-2,
     0x1.cd79b50727b8ep-12,
     0x1.cd795a8c33438p-12},
    {0x1.0000000000000p-1,
     0x1.47ae147ae147bp-7,
     0x1.0000000000000p-1,
     0x1.a383a8fc47f9ep-12,
     0x1.a3837489251cbp-12},
    {0x1.0000000000000p-1,
     0x1.47ae147ae147bp-7,
     0x1.028f5c28f5c29p-1,
     0x1.a3738d213712ep-14,
     0x1.a37358b21b421p-14},
    {0x1.0000000000000p-1,
     0x1.0624dd2f1a9fcp-10,
     0x1.3333333333333p-3,
     0x1.0e54e552928d9p-8,
     0x1.0e5508f59b6e1p-8},
    {0x1.0000000000000p-1,
     0x1.0624dd2f1a9fcp-10,
     0x1.3333333333333p-2,
     0x1.ecf1f2e248f1ap-10,
     0x1.ecf04f161babfp-10},
    {0x1.0000000000000p-1,
     0x1.0624dd2f1a9fcp-10,
     0x1.ccccccccccccdp-2,
     0x1.abf7ef2887930p-12,
     0x1.abf6892ab58d5p-12},
    {0x1.0000000000000p-1,
     0x1.0624dd2f1a9fcp-10,
     0x1.ff7ced916872bp-2,
     0x1.0c6fd2016fdbap-17,
     0x1.0c6e0d04e8253p-17},
    {0x1.0000000000000p-1,
     0x1.0624dd2f1a9fcp-10,
     0x1.0000000000000p-1,
     0x1.0c6f9d3a94eecp-18,
     0x1.0c6dd83ebf56ap-18},
    {0x1.0000000000000p-1,
     0x1.0624dd2f1a9fcp-10,
     0x1.004189374bc6ap-1,
     0x1.0c6f82d72c0cbp-20,
     0x1.0c658ccb25d63p-20},
    {0x1.0000000000000p+2,
     0x1.999999999999ap-4,
     0x1.3333333333333p+0,
     0x1.ad25b218e19d8p-5,
     0x1.ad25b877dd50bp-5},
    {0x1.0000000000000p+2,
     0x1.999999999999ap-4,
     0x1.3333333333333p+1,
     0x1.8fcbbe8c6a2f2p-6,
     0x1.8fcba28aa1572p-6},
    {0x1.0000000000000p+2,
     0x1.999999999999ap-4,
     0x1.ccccccccccccdp+1,
     0x1.7564d00f24285p-8,
     0x1.7564bb5f92107p-8},
    {0x1.0000000000000p+2,
     0x1.999999999999ap-4,
     0x1.ff7ced916872bp+1,
     0x1.6203a3f017e96p-11,
     0x1.6203632d3f25cp-11},
    {0x1.0000000000000p+2,
     0x1.999999999999ap-4,
     0x1.0000000000000p+2,
     0x1.47c84cc3a802ep-11,
     0x1.47c823c7587c4p-11},
    {0x1.0000000000000p+2,
     0x1.999999999999ap-4,
     0x1.0333333333333p+2,
     0x1.47b4a249fa78dp-13,
     0x1.47b3ac7dc46b6p-13},
    {0x1.0000000000000p+2,
     0x1.47ae147ae147bp-7,
     0x1.3333333333333p+0,
     0x1.520ad49faef4bp-8,
     0x1.520abf0e56b41p-8},
    {0x1.0000000000000p+2,
     0x1.47ae147ae147bp-7,
     0x1.3333333333333p+1,
     0x1.3457ae9bd9441p-9,
     0x1.345799af80c15p-9},
    {0x1.0000000000000p+2,
     0x1.47ae147ae147bp-7,
     0x1.ccccccccccccdp+1,
     0x1.0c27f5f58dffcp-11,
     0x1.0c27e93616938p-11},
    {0x1.0000000000000p+2,
     0x1.47ae147ae147bp-7,
     0x1.ff7ced916872bp+1,
     0x1.797d67858690fp-17,
     0x1.797cb542d062bp-17},
    {0x1.0000000000000p+2,
     0x1.47ae147ae147bp-7,
     0x1.0000000000000p+2,
     0x1.a36e84980b625p-18,
     0x1.a36e502a31c99p-18},
    {0x1.0000000000000p+2,
     0x1.47ae147ae147bp-7,
     0x1.0051eb851eb85p+2,
     0x1.a36e442b53bd0p-20,
     0x1.a368f10904090p-20},
    {0x1.0000000000000p+2,
     0x1.0624dd2f1a9fcp-10,
     0x1.3333333333333p+0,
     0x1.0df8a7948df18p-11,
     0x1.0dfbe233cbabfp-11},
    {0x1.0000000000000p+2,
     0x1.0624dd2f1a9fcp-10,
     0x1.3333333333333p+1,
     0x1.eb8972fe78ee6p-13,
     0x1.eb8050ba48fd2p-13},
    {0x1.0000000000000p+2,
     0x1.0624dd2f1a9fcp-10,
     0x1.ccccccccccccdp+1,
     0x1.a82f8f038caa0p-15,
     0x1.a827b45148af7p-15},
    {0x1.0000000000000p+2,
     0x1.0624dd2f1a9fcp-10,
     0x1.ff7ced916872bp+1,
     0x1.2dfd82a84fe3bp-21,
     0x1.2df6665aa98ddp-21},
    {0x1.0000000000000p+2,
     0x1.0624dd2f1a9fcp-10,
     0x1.0000000000000p+2,
     0x1.0c6f7a981b641p-24,
     0x1.0c65848cb25a8p-24},
    {0x1.0000000000000p+2,
     0x1.0624dd2f1a9fcp-10,
     0x1.00083126e978dp+2,
     0x1.0c6f7a2e8ed10p-26,
     0x1.0c2400231b6cap-26},
};

/// Deviation from the reference, relative, in units of the epsilon of the tested type.
template <typename T>
double ulpDistance(T computed, double reference) {
  if (reference == 0.0) {
    return computed == T{0} ? 0.0 : std::numeric_limits<double>::infinity();
  }
  return std::abs((static_cast<double>(computed) - reference) / reference) /
         static_cast<double>(std::numeric_limits<T>::epsilon());
}

template <typename T>
double reference(const Sample& sample) {
  return std::is_same_v<T, double> ? sample.increment : sample.singleIncrement;
}

} // namespace gaussiannucleation

TEST_CASE_TEMPLATE("Gaussian nucleation increment" * doctest::test_suite("numerical"),
                   RealT,
                   float,
                   double) {
  using namespace gaussiannucleation;
  constexpr double MaxUlp = 32.0;
  for (const auto& sample : Samples) {
    const auto increment = smoothStepIncrement<RealT>(static_cast<RealT>(sample.currentTime),
                                                      static_cast<RealT>(sample.dt),
                                                      static_cast<RealT>(sample.t0));
    REQUIRE(ulpDistance(increment, reference<RealT>(sample)) < MaxUlp);
  }
}

/**
 * The caller sums the increments into the initial stress, so what the sum telescopes to is what
 * the nucleation applies, and a step that goes backwards would unload the fault.
 */
TEST_CASE_TEMPLATE("Gaussian nucleation increments accumulate to the ramp" *
                       doctest::test_suite("numerical"),
                   RealT,
                   float,
                   double) {
  using namespace gaussiannucleation;
  for (const double t0 : {0.5, 1.0, 4.0}) {
    for (const double dt : {1e-1, 1e-2, 1e-3}) {
      const auto steps = static_cast<int>(1.2 * t0 / dt);
      auto sum = static_cast<RealT>(0);
      for (int i = 1; i <= steps; ++i) {
        const auto increment = smoothStepIncrement<RealT>(
            static_cast<RealT>(i * dt), static_cast<RealT>(dt), static_cast<RealT>(t0));
        REQUIRE(increment >= static_cast<RealT>(0));
        sum += increment;
      }
      const auto end = smoothStep<RealT>(static_cast<RealT>(steps * dt), static_cast<RealT>(t0));
      REQUIRE(std::abs(sum - end) < 64 * steps * std::numeric_limits<RealT>::epsilon());
    }
  }
}

TEST_CASE_TEMPLATE("Gaussian nucleation outside the ramp" * doctest::test_suite("numerical"),
                   RealT,
                   float,
                   double) {
  using namespace gaussiannucleation;
  constexpr auto T0 = static_cast<RealT>(2);
  // exactly representable, so that t - dt lands on the ends of the ramp rather than beside them
  constexpr auto Dt = static_cast<RealT>(0.25);
  constexpr auto Zero = static_cast<RealT>(0);
  constexpr auto One = static_cast<RealT>(1);
  REQUIRE(smoothStepIncrement<RealT>(-One, Dt, T0) == Zero);
  REQUIRE(smoothStepIncrement<RealT>(Zero, Dt, T0) == Zero);
  // the ramp is flat at both ends, so the steps just inside it still apply nothing measurable
  REQUIRE(smoothStepIncrement<RealT>(Dt, Dt, T0) >= Zero);
  REQUIRE(smoothStepIncrement<RealT>(T0 + Dt, Dt, T0) == Zero);
  REQUIRE(smoothStepIncrement<RealT>(10 * T0, Dt, T0) == Zero);
  // a step spanning the whole ramp applies all of it, in one go
  REQUIRE(smoothStepIncrement<RealT>(T0, 10 * T0, T0) == One);
  REQUIRE(smoothStepIncrement<RealT>(T0 + Dt, 10 * T0, T0) == One);
}

} // namespace seissol::unit_test
