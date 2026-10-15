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
 * t - dt is not representable, but the caller divides the increment by dt, so the reference is
 * the increment over the exact dt, not over the rounded previous time. singleIncrement is the
 * same quantity for t, dt and t0 rounded to single precision.
 */
constexpr std::size_t SampleCount = 54;
constexpr Sample Samples[SampleCount] = {
    {0x1.0000000000000p+0,
     0x1.999999999999ap-4,
     0x1.3333333333333p-2,
     0x1.b5691b271b877p-3,
     0x1.b5691b86c634fp-3},
    {0x1.0000000000000p+0,
     0x1.999999999999ap-4,
     0x1.3333333333333p-1,
     0x1.c2b3251202673p-4,
     0x1.c2b3236cd380ap-4},
    {0x1.0000000000000p+0,
     0x1.999999999999ap-4,
     0x1.ccccccccccccdp-1,
     0x1.f7fa5ef0427acp-6,
     0x1.f7fa6517403b9p-6},
    {0x1.0000000000000p+0,
     0x1.999999999999ap-4,
     0x1.ff7ced916872bp-1,
     0x1.4ff1be23a076fp-7,
     0x1.4ff1b93490be2p-7},
    {0x1.0000000000000p+0,
     0x1.999999999999ap-4,
     0x1.0000000000000p+0,
     0x1.4952e7a5a851ap-7,
     0x1.4952e84b24e11p-7},
    {0x1.0000000000000p+0,
     0x1.999999999999ap-4,
     0x1.0cccccccccccdp+0,
     0x1.48170661903aep-9,
     0x1.481730ba137d2p-9},
    {0x1.0000000000000p+0,
     0x1.47ae147ae147bp-7,
     0x1.3333333333333p-2,
     0x1.53e875e8eb0c9p-6,
     0x1.53e87505f180ap-6},
    {0x1.0000000000000p+0,
     0x1.47ae147ae147bp-7,
     0x1.3333333333333p-1,
     0x1.38236d6636813p-7,
     0x1.38236b59bb8bfp-7},
    {0x1.0000000000000p+0,
     0x1.47ae147ae147bp-7,
     0x1.ccccccccccccdp-1,
     0x1.164f3c5f16299p-9,
     0x1.164f402d3cef2p-9},
    {0x1.0000000000000p+0,
     0x1.47ae147ae147bp-7,
     0x1.ff7ced916872bp-1,
     0x1.f758e1054db2cp-14,
     0x1.f75898e16dbcdp-14},
    {0x1.0000000000000p+0,
     0x1.47ae147ae147bp-7,
     0x1.0000000000000p+0,
     0x1.a3738d2137121p-14,
     0x1.a3738be69c619p-14},
    {0x1.0000000000000p+0,
     0x1.47ae147ae147bp-7,
     0x1.0147ae147ae14p+0,
     0x1.a36f864b6e0fap-16,
     0x1.a36fb84461e8ep-16},
    {0x1.0000000000000p+0,
     0x1.0624dd2f1a9fcp-10,
     0x1.3333333333333p-2,
     0x1.0e204f4009f30p-9,
     0x1.0e204fc4e8c75p-9},
    {0x1.0000000000000p+0,
     0x1.0624dd2f1a9fcp-10,
     0x1.3333333333333p-1,
     0x1.ec23da9f30195p-11,
     0x1.ec23d9a48d044p-11},
    {0x1.0000000000000p+0,
     0x1.0624dd2f1a9fcp-10,
     0x1.ccccccccccccdp-1,
     0x1.a9ce7e83b7aacp-13,
     0x1.a9ce8699bd14ap-13},
    {0x1.0000000000000p+0,
     0x1.0624dd2f1a9fcp-10,
     0x1.ff7ced916872bp-1,
     0x1.92a7790993f77p-19,
     0x1.92a69836ee7cep-19},
    {0x1.0000000000000p+0,
     0x1.0624dd2f1a9fcp-10,
     0x1.0000000000000p+0,
     0x1.0c6f82d72bcaap-20,
     0x1.0c6f8482fd91dp-20},
    {0x1.0000000000000p+0,
     0x1.0624dd2f1a9fcp-10,
     0x1.0020c49ba5e35p+0,
     0x1.0c6f7c3e524d1p-22,
     0x1.0c7975d0df1a7p-22},
    {0x1.0000000000000p-1,
     0x1.999999999999ap-4,
     0x1.3333333333333p-3,
     0x1.795bf99e5f61ep-2,
     0x1.795bfad84dd1fp-2},
    {0x1.0000000000000p-1,
     0x1.999999999999ap-4,
     0x1.3333333333333p-2,
     0x1.06f205723b9d1p-2,
     0x1.06f2049bd1300p-2},
    {0x1.0000000000000p-1,
     0x1.999999999999ap-4,
     0x1.ccccccccccccdp-2,
     0x1.588ba2d789562p-4,
     0x1.588ba6464b4c5p-4},
    {0x1.0000000000000p-1,
     0x1.999999999999ap-4,
     0x1.ff7ced916872bp-2,
     0x1.51bb3f795c04ep-5,
     0x1.51bb3d43b479ap-5},
    {0x1.0000000000000p-1,
     0x1.999999999999ap-4,
     0x1.0000000000000p-1,
     0x1.4e51e9618b51dp-5,
     0x1.4e51ea0c11191p-5},
    {0x1.0000000000000p-1,
     0x1.999999999999ap-4,
     0x1.199999999999ap-1,
     0x1.4952e7a5a8510p-7,
     0x1.4952de98d889ep-7},
    {0x1.0000000000000p-1,
     0x1.47ae147ae147bp-7,
     0x1.3333333333333p-3,
     0x1.563bab238ff0ep-5,
     0x1.563baa44be5f0p-5},
    {0x1.0000000000000p-1,
     0x1.47ae147ae147bp-7,
     0x1.3333333333333p-2,
     0x1.3d4130efcc528p-6,
     0x1.3d412edba69b7p-6},
    {0x1.0000000000000p-1,
     0x1.47ae147ae147bp-7,
     0x1.ccccccccccccdp-2,
     0x1.23e5e14dd4397p-8,
     0x1.23e5e5154903bp-8},
    {0x1.0000000000000p-1,
     0x1.47ae147ae147bp-7,
     0x1.ff7ced916872bp-2,
     0x1.cd79b50727b81p-12,
     0x1.cd799054d0e6dp-12},
    {0x1.0000000000000p-1,
     0x1.47ae147ae147bp-7,
     0x1.0000000000000p-1,
     0x1.a383a8fc47f91p-12,
     0x1.a383a7c1951e5p-12},
    {0x1.0000000000000p-1,
     0x1.47ae147ae147bp-7,
     0x1.028f5c28f5c29p-1,
     0x1.a3738d2137114p-14,
     0x1.a373bf1b20a15p-14},
    {0x1.0000000000000p-1,
     0x1.0624dd2f1a9fcp-10,
     0x1.3333333333333p-3,
     0x1.0e54e552928d5p-8,
     0x1.0e54e5d827988p-8},
    {0x1.0000000000000p-1,
     0x1.0624dd2f1a9fcp-10,
     0x1.3333333333333p-2,
     0x1.ecf1f2e248f12p-10,
     0x1.ecf1f1e856862p-10},
    {0x1.0000000000000p-1,
     0x1.0624dd2f1a9fcp-10,
     0x1.ccccccccccccdp-2,
     0x1.abf7ef288792ap-12,
     0x1.abf7f74283a0fp-12},
    {0x1.0000000000000p-1,
     0x1.0624dd2f1a9fcp-10,
     0x1.ff7ced916872bp-2,
     0x1.0c6fd2016fdb4p-17,
     0x1.0c6f6202e5862p-17},
    {0x1.0000000000000p-1,
     0x1.0624dd2f1a9fcp-10,
     0x1.0000000000000p-1,
     0x1.0c6f9d3a94ee4p-18,
     0x1.0c6f9ee667099p-18},
    {0x1.0000000000000p-1,
     0x1.0624dd2f1a9fcp-10,
     0x1.004189374bc6ap-1,
     0x1.0c6f82d72c0bap-20,
     0x1.0c691a0d637b4p-20},
    {0x1.0000000000000p+2,
     0x1.999999999999ap-4,
     0x1.3333333333333p+0,
     0x1.ad25b218e19d1p-5,
     0x1.ad25b215cac09p-5},
    {0x1.0000000000000p+2,
     0x1.999999999999ap-4,
     0x1.3333333333333p+1,
     0x1.8fcbbe8c6a2ecp-6,
     0x1.8fcbbcf25e369p-6},
    {0x1.0000000000000p+2,
     0x1.999999999999ap-4,
     0x1.ccccccccccccdp+1,
     0x1.7564d00f2427fp-8,
     0x1.7564d5c861124p-8},
    {0x1.0000000000000p+2,
     0x1.999999999999ap-4,
     0x1.ff7ced916872bp+1,
     0x1.6203a3f017e8cp-11,
     0x1.62038e784f4e6p-11},
    {0x1.0000000000000p+2,
     0x1.999999999999ap-4,
     0x1.0000000000000p+2,
     0x1.47c84cc3a8025p-11,
     0x1.47c84d6799459p-11},
    {0x1.0000000000000p+2,
     0x1.999999999999ap-4,
     0x1.0333333333333p+2,
     0x1.47b4a249fa779p-13,
     0x1.47b3ffb431b92p-13},
    {0x1.0000000000000p+2,
     0x1.47ae147ae147bp-7,
     0x1.3333333333333p+0,
     0x1.520ad49faef46p-8,
     0x1.520ad3ba375aep-8},
    {0x1.0000000000000p+2,
     0x1.47ae147ae147bp-7,
     0x1.3333333333333p+1,
     0x1.3457ae9bd94b5p-9,
     0x1.3457ac95087e4p-9},
    {0x1.0000000000000p+2,
     0x1.47ae147ae147bp-7,
     0x1.ccccccccccccdp+1,
     0x1.0c27f5f58e062p-11,
     0x1.0c27f9c8d4567p-11},
    {0x1.0000000000000p+2,
     0x1.47ae147ae147bp-7,
     0x1.ff7ced916872bp+1,
     0x1.797d6785869ecp-17,
     0x1.797cd919edf08p-17},
    {0x1.0000000000000p+2,
     0x1.47ae147ae147bp-7,
     0x1.0000000000000p+2,
     0x1.a36e84980b75fp-18,
     0x1.a36e835d78525p-18},
    {0x1.0000000000000p+2,
     0x1.47ae147ae147bp-7,
     0x1.0051eb851eb85p+2,
     0x1.a36e442b53e45p-20,
     0x1.a369576ed4c53p-20},
    {0x1.0000000000000p+2,
     0x1.0624dd2f1a9fcp-10,
     0x1.3333333333333p+0,
     0x1.0df8a7948e123p-11,
     0x1.0df8a818e49aep-11},
    {0x1.0000000000000p+2,
     0x1.0624dd2f1a9fcp-10,
     0x1.3333333333333p+1,
     0x1.eb8972fe7929fp-13,
     0x1.eb89720351d0bp-13},
    {0x1.0000000000000p+2,
     0x1.0624dd2f1a9fcp-10,
     0x1.ccccccccccccdp+1,
     0x1.a82f8f038cdd7p-15,
     0x1.a82f97169a3c4p-15},
    {0x1.0000000000000p+2,
     0x1.0624dd2f1a9fcp-10,
     0x1.ff7ced916872bp+1,
     0x1.2dfd82a8500c5p-21,
     0x1.2dfca1356aed0p-21},
    {0x1.0000000000000p+2,
     0x1.0624dd2f1a9fcp-10,
     0x1.0000000000000p+2,
     0x1.0c6f7a981ba51p-24,
     0x1.0c6f7c43ed520p-24},
    {0x1.0000000000000p+2,
     0x1.0624dd2f1a9fcp-10,
     0x1.00083126e978dp+2,
     0x1.0c6f7a2e8f530p-26,
     0x1.0c37ed51e2b4ap-26},
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

  constexpr double MaxUlp = 16.0;
  for (const auto& sample : Samples) {
    CAPTURE(sample.currentTime);
    CAPTURE(sample.dt);
    CAPTURE(sample.t0);
    const auto increment = smoothStepIncrement<RealT>(static_cast<RealT>(sample.currentTime),
                                                      static_cast<RealT>(sample.dt),
                                                      static_cast<RealT>(sample.t0));
    CHECK(ulpDistance(increment, reference<RealT>(sample)) < MaxUlp);
  }
}

/**
 * The caller divides the increments by dt and integrates the rate over time, so the slip it
 * imposes is their sum, and a step that goes backwards would reverse it.
 */
TEST_CASE_TEMPLATE("Gaussian nucleation increments accumulate to the ramp" *
                       doctest::test_suite("numerical"),
                   RealT,
                   float,
                   double) {
  using namespace gaussiannucleation;
  for (const double t0 : {0.5, 1.0, 4.0}) {
    CAPTURE(t0);
    for (const double dt : {1e-1, 1e-2, 1e-3}) {
      CAPTURE(dt);
      const auto steps = static_cast<int>(1.2 * t0 / dt);
      auto sum = static_cast<RealT>(0);
      for (int i = 1; i <= steps; ++i) {
        const auto increment = smoothStepIncrement<RealT>(
            static_cast<RealT>(i * dt), static_cast<RealT>(dt), static_cast<RealT>(t0));
        CHECK(increment >= static_cast<RealT>(0));
        sum += increment;
      }
      const auto end = smoothStep<RealT>(static_cast<RealT>(steps * dt), static_cast<RealT>(t0));
      CHECK(std::abs(sum - end) < 64 * steps * std::numeric_limits<RealT>::epsilon());
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
  CHECK(smoothStepIncrement<RealT>(-One, Dt, T0) == Zero);
  CHECK(smoothStepIncrement<RealT>(Zero, Dt, T0) == Zero);
  // the ramp is flat at both ends, so the steps just inside it still apply nothing measurable
  CHECK(smoothStepIncrement<RealT>(Dt, Dt, T0) >= Zero);
  CHECK(smoothStepIncrement<RealT>(T0 + Dt, Dt, T0) == Zero);
  CHECK(smoothStepIncrement<RealT>(10 * T0, Dt, T0) == Zero);
  // a step spanning the whole ramp applies all of it, in one go
  CHECK(smoothStepIncrement<RealT>(T0, 10 * T0, T0) == One);
  CHECK(smoothStepIncrement<RealT>(T0 + Dt, 10 * T0, T0) == One);
}

} // namespace seissol::unit_test
