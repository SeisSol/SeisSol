// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "DynamicRupture/Misc.h"
#include "Initializer/Parameters/DRParameters.h"
#include "TestHelper.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>

namespace seissol::unit_test {

namespace stresssources {

using seissol::dr::FrictionLawParameters;
using seissol::dr::MaxStressSources;
using seissol::dr::stressAtTime;
using seissol::dr::stressSourceCount;
using seissol::dr::stressSourceFraction;

struct RampSample {
  double time;
  double riseTime;
  double onset;
  double fraction;
};

/**
 * The ramp at 50 decimal digits, rounded once to the nearest double. Rise time and onset differ
 * within every sample, so that passing them the wrong way round is caught.
 */
constexpr std::array<RampSample, 8> RampSamples{{
    {0x1.6666666666666p-1, 0x1.0000000000000p+1, 0x1.0000000000000p-1, 0x1.cd4cf1c6c4d43p-7},
    {0x1.8000000000000p+0, 0x1.0000000000000p+1, 0x1.0000000000000p-1, 0x1.6edd3122f2ea5p-1},
    {0x1.0000000000000p+1, 0x1.0000000000000p+1, 0x1.0000000000000p-1, 0x1.defac583c7dd4p-1},
    {0x1.0000000000000p-2, 0x1.0000000000000p+0, 0x0.0p+0, 0x1.1b1681e89d3a2p-2},
    {0x1.3333333333333p-1, 0x1.0000000000000p+0, 0x0.0p+0, 0x1.a73395c533373p-1},
    {0x1.ccccccccccccdp-1, 0x1.0000000000000p+0, 0x0.0p+0, 0x1.fadab461695ecp-1},
    {0x1.8000000000000p+1, 0x1.0000000000000p+2, 0x1.0000000000000p+0, 0x1.6edd3122f2ea5p-1},
    {0x1.2000000000000p+2, 0x1.0000000000000p+2, 0x1.0000000000000p+0, 0x1.f7efeac14b167p-1},
}};

/// a parameter set with the given rise times, the initial state appended as usual
FrictionLawParameters withSources(const std::vector<std::pair<double, double>>& nucleations) {
  seissol::initializer::parameters::DRParameters parameters{};
  parameters.nucleationCount = static_cast<std::uint32_t>(nucleations.size());
  for (std::size_t i = 0; i < nucleations.size(); ++i) {
    parameters.t0[i] = nucleations[i].first;
    parameters.s0[i] = nucleations[i].second;
  }
  return FrictionLawParameters(parameters);
}

/// a field of a face, as the initializer lays it out: the nucleations, then the initial state
void withField(real (*field)[seissol::dr::misc::NumPaddedPoints],
               std::uint32_t sourceCount,
               const std::vector<double>& nucleations) {
  for (std::uint32_t source = 0; source < sourceCount; ++source) {
    const auto value = source + 1 < sourceCount ? nucleations[source] : 0.0;
    for (std::uint32_t point = 0; point < seissol::dr::misc::NumPaddedPoints; ++point) {
      field[source][point] = static_cast<real>(value);
    }
  }
}

} // namespace stresssources

TEST_CASE("Stress source ramp" * doctest::test_suite("dynamicrupture")) {
  using namespace stresssources;

  for (const auto& sample : RampSamples) {
    const auto fraction = stressSourceFraction(static_cast<real>(sample.time),
                                               static_cast<real>(sample.riseTime),
                                               static_cast<real>(sample.onset));
    REQUIRE(static_cast<double>(fraction) ==
            AbsApprox(sample.fraction).epsilon(32 * std::numeric_limits<real>::epsilon()));
  }

  // the ends of the ramp are exact, and it stays within them
  constexpr auto RiseTime = static_cast<real>(2.0);
  constexpr auto Onset = static_cast<real>(0.5);
  REQUIRE(stressSourceFraction(Onset, RiseTime, Onset) == static_cast<real>(0.0));
  REQUIRE(stressSourceFraction(Onset - RiseTime, RiseTime, Onset) == static_cast<real>(0.0));
  REQUIRE(stressSourceFraction(Onset + RiseTime, RiseTime, Onset) == static_cast<real>(1.0));
  REQUIRE(stressSourceFraction(Onset + 10 * RiseTime, RiseTime, Onset) == static_cast<real>(1.0));

  auto previous = static_cast<real>(0.0);
  for (int i = 0; i <= 200; ++i) {
    const auto time = Onset + RiseTime * static_cast<real>(i) / 200;
    const auto fraction = stressSourceFraction(time, RiseTime, Onset);
    REQUIRE(fraction >= previous);
    REQUIRE(fraction <= static_cast<real>(1.0));
    previous = fraction;
  }
}

/**
 * A source without a rise time is what the initial state is, so the ramp has to be in full effect
 * at the onset itself and not one step after it.
 */
TEST_CASE("Stress source without a rise time" * doctest::test_suite("dynamicrupture")) {
  using namespace stresssources;
  for (const double onset : {-1.0, 0.0, 2.5}) {
    const auto s0 = static_cast<real>(onset);
    for (const auto riseTime : {static_cast<real>(0.0), static_cast<real>(-1.0)}) {
      REQUIRE(stressSourceFraction(s0, riseTime, s0) == static_cast<real>(1.0));
      REQUIRE(stressSourceFraction(s0 + static_cast<real>(1.0), riseTime, s0) ==
              static_cast<real>(1.0));
      REQUIRE(stressSourceFraction(std::nextafter(s0, static_cast<real>(-1e30)), riseTime, s0) ==
              static_cast<real>(0.0));
    }
  }
}

TEST_CASE("Stress sources of a parameter set" * doctest::test_suite("dynamicrupture")) {
  using namespace stresssources;

  const auto parameters = withSources({{1.0, 0.5}, {2.0, 3.0}});
  REQUIRE(parameters.sourceCount == 3);
  // the forced rupture ramp is none of the sources and keeps the rise time of the first nucleation
  REQUIRE(parameters.forcedRuptureRiseTime == static_cast<real>(1.0));

  REQUIRE(withSources({}).sourceCount == 1);
  // the initial state has neither a rise time nor an onset, so it is in effect from the start
  REQUIRE(stressSourceFraction(static_cast<real>(0.0), static_cast<real>(0.0),
                               static_cast<real>(0.0)) == static_cast<real>(1.0));
}

TEST_CASE("Stress of a point over its sources" * doctest::test_suite("dynamicrupture")) {
  using namespace stresssources;
  using seissol::dr::misc::NumPaddedPoints;

  const auto parameters = withSources({{1.0, 0.5}, {2.0, 3.0}});
  constexpr std::uint32_t Point = 3;

  alignas(Alignment) static real sources[MaxStressSources][6][NumPaddedPoints]{};
  alignas(Alignment) static real onsets[MaxStressSources][NumPaddedPoints]{};
  alignas(Alignment) static real riseTimes[MaxStressSources][NumPaddedPoints]{};
  withField(onsets, parameters.sourceCount, {0.5, 3.0});
  withField(riseTimes, parameters.sourceCount, {1.0, 2.0});
  for (std::uint32_t source = 0; source < parameters.sourceCount; ++source) {
    for (std::size_t component = 0; component < 6; ++component) {
      sources[source][component][Point] =
          static_cast<real>(1 + source) * static_cast<real>(1 + component);
    }
  }

  // before every nucleation has started, a point carries the initial state and nothing else
  {
    const auto stress = stressAtTime(
        sources, riseTimes, onsets, parameters.sourceCount, Point, static_cast<real>(0.0));
    for (std::size_t component = 0; component < 6; ++component) {
      REQUIRE(stress[component] == sources[parameters.sourceCount - 1][component][Point]);
    }
  }

  // once they have, the fractions they have reached are what enters
  {
    const auto time = static_cast<real>(4.0);
    const auto stress =
        stressAtTime(sources, riseTimes, onsets, parameters.sourceCount, Point, time);
    for (std::size_t component = 0; component < 6; ++component) {
      real expected = 0;
      for (std::uint32_t source = 0; source < parameters.sourceCount; ++source) {
        expected += sources[source][component][Point] *
                    stressSourceFraction(
                        time, riseTimes[source][Point], onsets[source][Point]);
      }
      REQUIRE(stress[component] == expected);
    }
  }
}

/**
 * The stress of a point is a function of the time alone. Walking the times in one order and then
 * in another has to give the very same values, down to the last bit: anything that accumulated
 * across calls, or that remembered where it was last asked, would show up here.
 */
TEST_CASE("Stress of a point does not depend on the order it is asked in" *
          doctest::test_suite("dynamicrupture")) {
  using namespace stresssources;
  using seissol::dr::misc::NumPaddedPoints;

  const auto parameters = withSources({{1.0, 0.5}, {2.0, 3.0}, {0.0, 2.0}});
  constexpr std::uint32_t Point = 5;

  alignas(Alignment) static real sources[MaxStressSources][6][NumPaddedPoints]{};
  alignas(Alignment) static real onsets[MaxStressSources][NumPaddedPoints]{};
  alignas(Alignment) static real riseTimes[MaxStressSources][NumPaddedPoints]{};
  withField(onsets, parameters.sourceCount, {0.5, 3.0, 2.0});
  withField(riseTimes, parameters.sourceCount, {1.0, 2.0, 0.0});
  for (std::uint32_t source = 0; source < parameters.sourceCount; ++source) {
    for (std::size_t component = 0; component < 6; ++component) {
      sources[source][component][Point] = static_cast<real>(7 * source + component) - 5;
    }
  }

  constexpr std::size_t Steps = 64;
  std::array<std::array<real, 6>, Steps> forward{};
  for (std::size_t step = 0; step < Steps; ++step) {
    const auto time = static_cast<real>(6.0 * static_cast<double>(step) / Steps);
    forward[step] =
        stressAtTime(sources, riseTimes, onsets, parameters.sourceCount, Point, time);
  }
  for (std::size_t step = Steps; step-- > 0;) {
    const auto time = static_cast<real>(6.0 * static_cast<double>(step) / Steps);
    const auto backward =
        stressAtTime(sources, riseTimes, onsets, parameters.sourceCount, Point, time);
    for (std::size_t component = 0; component < 6; ++component) {
      REQUIRE(backward[component] == forward[step][component]);
    }
  }
}

/// the onset is a field, so two points of a face may reach the same source at different times
TEST_CASE("Stress sources with an onset per point" * doctest::test_suite("dynamicrupture")) {
  using namespace stresssources;
  using seissol::dr::misc::NumPaddedPoints;

  const auto parameters = withSources({{1.0, 0.0}});
  alignas(Alignment) static real sources[MaxStressSources][6][NumPaddedPoints]{};
  alignas(Alignment) static real onsets[MaxStressSources][NumPaddedPoints]{};
  alignas(Alignment) static real riseTimes[MaxStressSources][NumPaddedPoints]{};

  for (std::uint32_t point = 0; point < NumPaddedPoints; ++point) {
    sources[0][0][point] = static_cast<real>(4.0);
    sources[1][0][point] = static_cast<real>(1.0);
    // the nucleation sweeps across the face, and takes longer the further it gets
    onsets[0][point] = static_cast<real>(point);
    riseTimes[0][point] = static_cast<real>(1.0) + static_cast<real>(point) / 8;
    onsets[1][point] = static_cast<real>(0.0);
    riseTimes[1][point] = static_cast<real>(0.0);
  }

  const auto time = static_cast<real>(3.5);
  const auto count = parameters.sourceCount;
  for (std::uint32_t point = 0; point < NumPaddedPoints; ++point) {
    const auto stress = stressAtTime(sources, riseTimes, onsets, count, point, time);
    const auto expected =
        static_cast<real>(1.0) +
        static_cast<real>(4.0) *
            stressSourceFraction(time, riseTimes[0][point], onsets[0][point]);
    REQUIRE(stress[0] == expected);
  }
  // it has passed the first points and not yet reached the last
  REQUIRE(stressAtTime(sources, riseTimes, onsets, count, 0, time)[0] == static_cast<real>(5.0));
  REQUIRE(stressAtTime(sources, riseTimes, onsets, count, NumPaddedPoints - 1, time)[0] ==
          static_cast<real>(1.0));
}

} // namespace seissol::unit_test
