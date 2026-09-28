// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Numerical/RegularizedYoffe.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <utility>

namespace seissol::unit_test {
namespace regularizedyoffe {

using seissol::regularizedYoffe::regularizedYoffe;

struct Sample {
  double tauS;
  double tauR;
  double time;
  double value;
  double singleValue;
};

/**
 * The convolution of the Yoffe function with the triangle, Tinti et al. 2005 eq 3, evaluated at
 * 40 decimal digits with the kinks and the endpoint singularities as split points. That is an
 * independent path to the value: the closed form under test is the analytically integrated
 * version of the same convolution, and the two agree to 1e-22 where both are exact.
 *
 * The sample times are the six points at which the closed form switches branches, plus the
 * midpoint of every interval between them. Both regimes are covered, tauR > 2 tauS and
 * tauR <= 2 tauS, the latter being the branch that is easy to get wrong because the intervals
 * change order. singleValue is the same quantity for arguments rounded to single precision.
 */
constexpr std::size_t SampleCount = 48;
constexpr Sample Samples[SampleCount] = {
    {0x1.4189374bc6a7fp-3,
     0x1.2c083126e978dp-1,
     0x1.4189374bc6a7fp-4,
     0x1.f3b07c542300bp-1,
     0x1.f3b07b1476228p-1},
    {0x1.4189374bc6a7fp-3,
     0x1.2c083126e978dp-1,
     0x1.4189374bc6a7fp-3,
     0x1.5c4ef17c7cb45p+1,
     0x1.5c4ef09e8774bp+1},
    {0x1.4189374bc6a7fp-3,
     0x1.2c083126e978dp-1,
     0x1.e24dd2f1a9fbep-3,
     0x1.7c7afc6215285p+1,
     0x1.7c7afc43b2422p+1},
    {0x1.4189374bc6a7fp-3,
     0x1.2c083126e978dp-1,
     0x1.4189374bc6a7fp-2,
     0x1.02881ae0c2e55p+1,
     0x1.02881a41a2be3p+1},
    {0x1.4189374bc6a7fp-3,
     0x1.2c083126e978dp-1,
     0x1.cccccccccccccp-2,
     0x1.1d6dc23383044p+0,
     0x1.1d6dc365032edp+0},
    {0x1.4189374bc6a7fp-3,
     0x1.2c083126e978dp-1,
     0x1.2c083126e978dp-1,
     0x1.5092b08309109p-1,
     0x1.5092af6764bc3p-1},
    {0x1.4189374bc6a7fp-3,
     0x1.2c083126e978dp-1,
     0x1.54395810624ddp-1,
     0x1.9b78a1175690ap-2,
     0x1.9b78a66218cffp-2},
    {0x1.4189374bc6a7fp-3,
     0x1.2c083126e978dp-1,
     0x1.7c6a7ef9db22dp-1,
     0x1.47092c72b7008p-3,
     0x1.4709404f8473ep-3},
    {0x1.4189374bc6a7fp-3,
     0x1.2c083126e978dp-1,
     0x1.a49ba5e353f7cp-1,
     0x1.bf862adbd4ef4p-6,
     0x1.bf8645a0c52d1p-6},
    {0x1.4189374bc6a7fp-3,
     0x1.2c083126e978dp-1,
     0x1.cccccccccccccp-1,
     0x1.eb4ea2e22f5d5p-132,
     0x1.eb4e9f59be610p-57},
    {0x1.999999999999ap-5,
     0x1.0000000000000p+1,
     0x1.999999999999ap-6,
     0x1.e549d3023ba44p-1,
     0x1.e549d2c56b700p-1},
    {0x1.999999999999ap-5,
     0x1.0000000000000p+1,
     0x1.999999999999ap-5,
     0x1.56b830396b92fp+1,
     0x1.56b8300e5d47bp+1},
    {0x1.999999999999ap-5,
     0x1.0000000000000p+1,
     0x1.3333333333334p-4,
     0x1.822d6763bc95fp+1,
     0x1.822d66f462cf9p+1},
    {0x1.999999999999ap-5,
     0x1.0000000000000p+1,
     0x1.999999999999ap-4,
     0x1.19780297a060dp+1,
     0x1.19780273a58a8p+1},
    {0x1.999999999999ap-5,
     0x1.0000000000000p+1,
     0x1.0cccccccccccdp+0,
     0x1.46046c7f7b767p-2,
     0x1.46046d88840e1p-2},
    {0x1.999999999999ap-5,
     0x1.0000000000000p+1,
     0x1.0000000000000p+1,
     0x1.97abfda046568p-5,
     0x1.97abfdd4abd3dp-5},
    {0x1.999999999999ap-5,
     0x1.0000000000000p+1,
     0x1.0333333333333p+1,
     0x1.0a78df4b91a93p-5,
     0x1.0a78b81fb1035p-5},
    {0x1.999999999999ap-5,
     0x1.0000000000000p+1,
     0x1.0666666666666p+1,
     0x1.ba2a0a541c277p-7,
     0x1.ba2a50de03cfep-7},
    {0x1.999999999999ap-5,
     0x1.0000000000000p+1,
     0x1.099999999999ap+1,
     0x1.37cfa66c8c354p-9,
     0x1.37cf475138594p-9},
    {0x1.999999999999ap-5,
     0x1.0000000000000p+1,
     0x1.0cccccccccccdp+1,
     0x0.0p+0,
     0x1.434390f49b134p-54},
    {0x1.47ae147ae147bp-6,
     0x1.0000000000000p+2,
     0x1.47ae147ae147bp-7,
     0x1.0f8e7790577adp+0,
     0x1.0f8e77c348b70p+0},
    {0x1.47ae147ae147bp-6,
     0x1.0000000000000p+2,
     0x1.47ae147ae147bp-6,
     0x1.7ff15ac0f9131p+1,
     0x1.7ff15b0908c97p+1},
    {0x1.47ae147ae147bp-6,
     0x1.0000000000000p+2,
     0x1.eb851eb851eb8p-6,
     0x1.b19d61b7b0e3dp+1,
     0x1.b19d62092ac4fp+1},
    {0x1.47ae147ae147bp-6,
     0x1.0000000000000p+2,
     0x1.47ae147ae147bp-5,
     0x1.3d86457a8179bp+1,
     0x1.3d8645b64e436p+1},
    {0x1.47ae147ae147bp-6,
     0x1.0000000000000p+2,
     0x1.028f5c28f5c29p+1,
     0x1.45f3b8e0ff3c7p-3,
     0x1.45f3b913ed037p-3},
    {0x1.47ae147ae147bp-6,
     0x1.0000000000000p+2,
     0x1.0000000000000p+2,
     0x1.68987a49fef37p-7,
     0x1.68987a0603473p-7},
    {0x1.47ae147ae147bp-6,
     0x1.0000000000000p+2,
     0x1.00a3d70a3d70ap+2,
     0x1.d9476b4c54fc3p-8,
     0x1.d945c87e4af27p-8},
    {0x1.47ae147ae147bp-6,
     0x1.0000000000000p+2,
     0x1.0147ae147ae14p+2,
     0x1.89c6cacbd647dp-9,
     0x1.89c7069e5d424p-9},
    {0x1.47ae147ae147bp-6,
     0x1.0000000000000p+2,
     0x1.01eb851eb851ep+2,
     0x1.164ae6c0a661fp-11,
     0x1.164726ab64948p-11},
    {0x1.47ae147ae147bp-6,
     0x1.0000000000000p+2,
     0x1.028f5c28f5c29p+2,
     0x0.0p+0,
     0x1.06219e88ae71fp-56},
    {0x1.999999999999ap-2,
     0x1.0000000000000p-1,
     0x1.999999999999ap-3,
     0x1.492b897df4391p-1,
     0x1.492b895103b53p-1},
    {0x1.999999999999ap-2,
     0x1.0000000000000p-1,
     0x1.999999999999ap-2,
     0x1.ba8b179e01b07p+0,
     0x1.ba8b175a4078cp+0},
    {0x1.999999999999ap-2,
     0x1.0000000000000p-1,
     0x1.0000000000000p-1,
     0x1.e102b36480f43p+0,
     0x1.e102b323f075bp+0},
    {0x1.999999999999ap-2,
     0x1.0000000000000p-1,
     0x1.3333333333334p-1,
     0x1.aed476820bc6ep+0,
     0x1.aed4755cc2064p+0},
    {0x1.999999999999ap-2,
     0x1.0000000000000p-1,
     0x1.6666666666666p-1,
     0x1.49a76a382e781p+0,
     0x1.49a76b9c49dd7p+0},
    {0x1.999999999999ap-2,
     0x1.0000000000000p-1,
     0x1.999999999999ap-1,
     0x1.85d3a187f93e1p-1,
     0x1.85d3a0defe1d4p-1},
    {0x1.999999999999ap-2,
     0x1.0000000000000p-1,
     0x1.ccccccccccccdp-1,
     0x1.8dfa9936fe179p-2,
     0x1.8dfa9f30b68bdp-2},
    {0x1.999999999999ap-2,
     0x1.0000000000000p-1,
     0x1.199999999999ap+0,
     0x1.e5895c7d18803p-5,
     0x1.e58956944a588p-5},
    {0x1.999999999999ap-2,
     0x1.0000000000000p-1,
     0x1.4cccccccccccdp+0,
     0x0.0p+0,
     0x1.80228cbbcfce4p-59},
    {0x1.3333333333333p-2,
     0x1.0000000000000p-1,
     0x1.3333333333333p-3,
     0x1.806977b97fcc2p-1,
     0x1.80697730d7474p-1},
    {0x1.3333333333333p-2,
     0x1.0000000000000p-1,
     0x1.3333333333333p-2,
     0x1.06604291b2756p+1,
     0x1.0660422d072e8p+1},
    {0x1.3333333333333p-2,
     0x1.0000000000000p-1,
     0x1.cccccccccccccp-2,
     0x1.0e68123262807p+1,
     0x1.0e6812a11e745p+1},
    {0x1.3333333333333p-2,
     0x1.0000000000000p-1,
     0x1.0000000000000p-1,
     0x1.e179b63c86b70p+0,
     0x1.e179b70d5f1cbp+0},
    {0x1.3333333333333p-2,
     0x1.0000000000000p-1,
     0x1.3333333333333p-1,
     0x1.2d9b6780529c6p+0,
     0x1.2d9b667b4b92dp+0},
    {0x1.3333333333333p-2,
     0x1.0000000000000p-1,
     0x1.4cccccccccccdp-1,
     0x1.c31ecb51e7605p-1,
     0x1.c31ed1a694a83p-1},
    {0x1.3333333333333p-2,
     0x1.0000000000000p-1,
     0x1.999999999999ap-1,
     0x1.3d0c9386f291dp-2,
     0x1.3d0c941e254fbp-2},
    {0x1.3333333333333p-2,
     0x1.0000000000000p-1,
     0x1.e666666666667p-1,
     0x1.99017c641f13ap-5,
     0x1.99018abefcc63p-5},
    {0x1.3333333333333p-2, 0x1.0000000000000p-1, 0x1.199999999999ap+0, 0x0.0p+0, 0x0.0p+0},
};

/// Deviation from the reference, absolute, measured against the peak of the function.
template <typename T>
double peakDistance(T computed, double reference, double peak) {
  return std::abs((static_cast<double>(computed) - reference) / peak) /
         static_cast<double>(std::numeric_limits<T>::epsilon());
}

template <typename T>
double reference(const Sample& sample) {
  return std::is_same_v<T, double> ? sample.value : sample.singleValue;
}

/// 2 / (pi tauR) sqrt((tauR - t) / t) convolved with a unit triangle peaks near this
template <typename T>
double peakOf(double tauS, double tauR) {
  double peak = 0;
  for (int i = 1; i < 200; ++i) {
    const double time = (tauR + 2 * tauS) * i / 200.0;
    peak = std::max(peak,
                    static_cast<double>(regularizedYoffe<T>(
                        static_cast<T>(time), static_cast<T>(tauS), static_cast<T>(tauR))));
  }
  return peak;
}

} // namespace regularizedyoffe

TEST_CASE_TEMPLATE("Regularized Yoffe function" * doctest::test_suite("numerical"),
                   RealT,
                   float,
                   double) {
  using namespace regularizedyoffe;
  // the closed form is a second difference of step tauS over a function of scale tauR, so it
  // gives up (tauR / tauS)^2 in relative accuracy; the bound admits that factor and no more,
  // which is what keeps the test honest about where the accuracy actually goes
  for (const auto& sample : Samples) {
    CAPTURE(sample);
    const double ratio = sample.tauR / sample.tauS;
    const double maxUlp = 8.0 * (1.0 + ratio * ratio);
    const double peak = peakOf<RealT>(sample.tauS, sample.tauR);
    const auto value = regularizedYoffe<RealT>(static_cast<RealT>(sample.time),
                                               static_cast<RealT>(sample.tauS),
                                               static_cast<RealT>(sample.tauR));
    CHECK(peakDistance(value, reference<RealT>(sample), peak) < maxUlp);
  }
}

/**
 * Every branch boundary is a point at which the closed form divides by a difference that is
 * zero there, and the value has to come out the same from either side.
 */
TEST_CASE_TEMPLATE("Regularized Yoffe is continuous across its branches" *
                       doctest::test_suite("numerical"),
                   RealT,
                   float,
                   double) {
  using namespace regularizedyoffe;
  for (const auto& pair : {std::pair{0.157, 0.586}, {0.05, 2.0}, {0.4, 0.5}, {0.3, 0.5}}) {
    const auto [tauS, tauR] = pair;
    const double peak = peakOf<RealT>(tauS, tauR);
    const double ratio = tauR / tauS;
    const double tolerance = 8.0 * (1.0 + ratio * ratio) * peak *
                             static_cast<double>(std::numeric_limits<RealT>::epsilon());
    for (const double breakpoint : {0.0, tauS, 2 * tauS, tauR, tauR + tauS, tauR + 2 * tauS}) {
      const auto at = static_cast<RealT>(breakpoint);
      const auto below = std::nextafter(at, static_cast<RealT>(-1));
      const auto above = std::nextafter(at, static_cast<RealT>(1e30));
      for (const auto time : {below, at, above}) {
        const auto value =
            regularizedYoffe<RealT>(time, static_cast<RealT>(tauS), static_cast<RealT>(tauR));
        CHECK(std::isfinite(static_cast<double>(value)));
        // the ends of the support are genuine zeros, where the formula rounds to either sign
        CHECK(static_cast<double>(value) > -tolerance);
      }
      const auto jump =
          regularizedYoffe<RealT>(above, static_cast<RealT>(tauS), static_cast<RealT>(tauR)) -
          regularizedYoffe<RealT>(below, static_cast<RealT>(tauS), static_cast<RealT>(tauR));
      CHECK(std::abs(static_cast<double>(jump)) < tolerance);
    }
  }
}

/**
 * The source time function carries unit slip and vanishes outside [0, tauR + 2 tauS]; both are
 * properties of the convolution rather than of the closed form, so they catch a wrong branch.
 */
TEST_CASE_TEMPLATE("Regularized Yoffe has unit integral and compact support" *
                       doctest::test_suite("numerical"),
                   RealT,
                   float,
                   double) {
  using namespace regularizedyoffe;
  for (const auto& pair : {std::pair{0.157, 0.586}, {0.05, 2.0}, {0.4, 0.5}, {0.3, 0.5}}) {
    const auto [tauS, tauR] = pair;
    const double support = tauR + 2 * tauS;
    constexpr std::size_t Intervals = 20000;
    const double h = support / Intervals;
    double integral = 0;
    for (std::size_t i = 1; i < Intervals; ++i) {
      integral += static_cast<double>(regularizedYoffe<RealT>(
          static_cast<RealT>(i * h), static_cast<RealT>(tauS), static_cast<RealT>(tauR)));
    }
    CHECK(integral * h == doctest::Approx(1.0).epsilon(1e-3));
    for (const double outside : {-1.0, -1e-3, support + 1e-3, 10 * support}) {
      CHECK(regularizedYoffe<RealT>(static_cast<RealT>(outside),
                                    static_cast<RealT>(tauS),
                                    static_cast<RealT>(tauR)) == static_cast<RealT>(0));
    }
  }
}

/**
 * The wider accumulator is what a single precision run reaches for when tauR / tauS is large
 * enough that the second difference eats its mantissa. It has to stay available and has to
 * pin the error at the final rounding, independently of the ratio.
 */
TEST_CASE("Regularized Yoffe in a wider accumulator" * doctest::test_suite("numerical")) {
  using namespace regularizedyoffe;
  // the final rounding of the result, and nothing of the second difference
  constexpr double MaxUlp = 16.0;
  for (const auto& sample : Samples) {
    CAPTURE(sample);
    const double peak = peakOf<float>(sample.tauS, sample.tauR);
    const auto value = regularizedYoffe<float, double>(static_cast<float>(sample.time),
                                                       static_cast<float>(sample.tauS),
                                                       static_cast<float>(sample.tauR));
    CHECK(peakDistance(value, sample.singleValue, peak) < MaxUlp);
  }
}

} // namespace seissol::unit_test
