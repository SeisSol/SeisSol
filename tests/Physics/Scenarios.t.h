// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/ConfigRegistry.h"
#include "Common/Real.h"
#include "Equations/Datastructures.h"
#include "Initializer/Parameters/InitializationParameters.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Initializer/Parameters/ParameterReader.h"
#include "Initializer/Typedefs.h"
#include "Model/MaterialType.h"
#include "Physics/Scenario/Registry.h"
#include "Physics/Scenario/Scenarios.h"
#include "TestConfigs.h"
#include "TestHelper.h"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <math.h>
#include <type_traits>
#include <vector>
#include <yaml-cpp/yaml.h>
#include <yateto.h>

namespace seissol::unit_test {

namespace scenarios {

/// The configurations that hold a fluid and whose material is set up from rho, lambda and mu
/// alone, i.e. those the fluid scenarios are written for.
template <typename Cfg>
struct HoldsFluid : std::bool_constant<Cfg::MaterialType == model::MaterialType::Elastic ||
                                       Cfg::MaterialType == model::MaterialType::Acoustic> {};

/// The configurations whose material is set up from rho, lambda and mu alone. The attenuating ones
/// then carry no relaxation (omega = theta = 0), which leaves their memory variables decoupled.
template <typename Cfg>
struct SetUpByLame : std::bool_constant<HoldsFluid<Cfg>::value ||
                                        Cfg::MaterialType == model::MaterialType::Viscoelastic ||
                                        Cfg::MaterialType == model::MaterialType::Viscoacoustic> {};

/// Fills every entry before a scenario is evaluated.
constexpr double Sentinel = -1234.5;

/// The values of `pointCount` points of the configuration `Cfg`, followed by a guard as large as
/// the elastic quantities (nine) would need in addition. Nothing may be written into the guard.
template <typename Cfg>
std::vector<Real<Cfg>> guardedBuffer(std::size_t pointCount) {
  constexpr auto NumQuantities = model::MaterialOf<Cfg>::NumQuantities;
  return std::vector<Real<Cfg>>(pointCount * (NumQuantities + 9), Sentinel);
}

template <typename Cfg>
yateto::DenseTensorView<2, Real<Cfg>, unsigned> dofsView(std::vector<Real<Cfg>>& buffer,
                                                         std::size_t pointCount) {
  constexpr auto NumQuantities = model::MaterialOf<Cfg>::NumQuantities;
  return yateto::DenseTensorView<2, Real<Cfg>, unsigned>(
      buffer.data(), {static_cast<unsigned>(pointCount), static_cast<unsigned>(NumQuantities)});
}

template <typename Cfg>
void checkGuard(const std::vector<Real<Cfg>>& buffer, std::size_t pointCount) {
  constexpr auto NumQuantities = model::MaterialOf<Cfg>::NumQuantities;
  for (std::size_t i = pointCount * NumQuantities; i < buffer.size(); ++i) {
    CAPTURE(i);
    CHECK(buffer[i] == Sentinel);
  }
}

/// Checks the state of a fluid at `point`: `stress` (i.e. -p) on the normal stresses, no shear
/// stress, and `velocity` in the velocity of the configuration.
template <typename Cfg>
void checkFluidState(const yateto::DenseTensorView<2, Real<Cfg>, unsigned>& dofs,
                     std::size_t point,
                     double stress,
                     const std::array<double, 3>& velocity,
                     double tolerance) {
  using MaterialT = model::MaterialOf<Cfg>;
  CAPTURE(point);
  // a stress tensor has three normal components; the acoustic materials carry a single one
  constexpr std::size_t NormalComponents = std::min<std::size_t>(3, MaterialT::TractionComponents);
  for (std::size_t j = 0; j < MaterialT::TractionComponents; ++j) {
    CAPTURE(j);
    CHECK(dofs(point, j) == AbsApprox(j < NormalComponents ? stress : 0.0).epsilon(tolerance));
  }
  for (std::size_t d = 0; d < velocity.size(); ++d) {
    CAPTURE(d);
    CHECK(dofs(point, MaterialT::VelocityOffset + d) == AbsApprox(velocity[d]).epsilon(tolerance));
  }
}

template <typename Cfg>
model::MaterialOf<Cfg> fluid(double rho, double lambda) {
  model::MaterialOf<Cfg> material;
  material.rho = rho;
  material.setLameParameters(0.0, lambda);
  return material;
}

} // namespace scenarios

TEST_CASE_TEMPLATE_DEFINE(
    "AcousticTravellingWaveITM writes a plane wave into the pressure and velocity of the "
    "configuration" *
        doctest::test_suite("physics"),
    Cfg,
    TravellingWaveItmPlaneWave) {
  using namespace scenarios;
  constexpr double Tolerance = 100 * std::numeric_limits<Real<Cfg>>::epsilon();

  constexpr double Rho = 2.0;
  constexpr double Lambda = 8.0;
  constexpr double K = 3.0;
  const double c = std::sqrt(Lambda / Rho);

  auto material = fluid<Cfg>(Rho, Lambda);
  const CellMaterialData materialData{&material, {}};

  // With a factor of one, the mirror leaves the wave as it is. So before, during and after it, the
  // plain travelling wave p = rho c cos(k x - c k t), u = cos(k x - c k t) has to come out.
  AcousticTravellingWaveParametersITM parameters{};
  parameters.k = K;
  parameters.itmStartingTime = 1.0;
  parameters.itmDuration = 0.5;
  parameters.itmVelocityScalingFactor = 1.0;
  const physics::AcousticTravellingWaveITM wave(materialData, parameters, configIdOf<Cfg>());

  const std::array<std::array<double, 3>, 3> points{
      {{0.1, 0.2, 0.3}, {0.45, -0.7, 1.1}, {0.9, 0.5, -0.4}}};
  for (const double time : {0.5, 1.25, 2.0}) {
    CAPTURE(time);
    auto buffer = guardedBuffer<Cfg>(points.size());
    auto dofs = dofsView<Cfg>(buffer, points.size());
    wave.evaluate(time, points.data(), points.size(), materialData, dofs);
    for (std::size_t i = 0; i < points.size(); ++i) {
      const double phase = K * points[i][0] - c * K * time;
      checkFluidState<Cfg>(
          dofs, i, -Rho * c * std::cos(phase), {std::cos(phase), 0.0, 0.0}, Tolerance);
    }
    checkGuard<Cfg>(buffer, points.size());
  }
}

TEST_CASE_TEMPLATE_DEFINE(
    "AcousticTravellingWaveITM keeps the pressure and scales the velocity at the mirror" *
        doctest::test_suite("physics"),
    Cfg,
    TravellingWaveItmMirror) {
  using namespace scenarios;
  constexpr auto VelocityOffset = model::MaterialOf<Cfg>::VelocityOffset;

  constexpr double Start = 1.0;
  constexpr double Duration = 0.5;
  constexpr double N = 2.0;
  // the values right after a switch are taken one double later than the switch itself
  constexpr double SwitchTolerance =
      std::max(1.0e-10, 100 * static_cast<double>(std::numeric_limits<Real<Cfg>>::epsilon()));

  auto material = fluid<Cfg>(2.0, 8.0);
  const CellMaterialData materialData{&material, {}};

  AcousticTravellingWaveParametersITM parameters{};
  parameters.k = 3.0;
  parameters.itmStartingTime = Start;
  parameters.itmDuration = Duration;
  parameters.itmVelocityScalingFactor = N;
  const physics::AcousticTravellingWaveITM wave(materialData, parameters, configIdOf<Cfg>());

  const std::array<std::array<double, 3>, 2> points{{{0.1, 0.2, 0.3}, {0.9, 0.5, -0.4}}};
  const auto evaluateAt = [&](double time) {
    auto buffer = guardedBuffer<Cfg>(points.size());
    auto dofs = dofsView<Cfg>(buffer, points.size());
    wave.evaluate(time, points.data(), points.size(), materialData, dofs);
    checkGuard<Cfg>(buffer, points.size());
    return buffer;
  };
  // The mirror scales the velocity by the factor n when it starts, and by 1 / n when it ends; the
  // pressure passes both switches continuously.
  const auto checkSwitch = [&](double switchTime, double velocityFactor) {
    auto before = evaluateAt(switchTime);
    auto after = evaluateAt(std::nextafter(switchTime, std::numeric_limits<double>::infinity()));
    const auto dofsBefore = dofsView<Cfg>(before, points.size());
    const auto dofsAfter = dofsView<Cfg>(after, points.size());
    for (std::size_t i = 0; i < points.size(); ++i) {
      CAPTURE(i);
      CHECK(dofsAfter(i, 0) == AbsApprox(dofsBefore(i, 0)).epsilon(SwitchTolerance));
      CHECK(dofsAfter(i, VelocityOffset) ==
            AbsApprox(velocityFactor * dofsBefore(i, VelocityOffset)).epsilon(SwitchTolerance));
    }
  };

  SUBCASE("start") { checkSwitch(Start, N); }
  SUBCASE("end") { checkSwitch(Start + Duration, 1.0 / N); }
}

TEST_CASE_TEMPLATE_DEFINE("Ocean writes a fluid stress and the velocity of the configuration" *
                              doctest::test_suite("physics"),
                          Cfg,
                          OceanLayout) {
  using namespace scenarios;
  constexpr double Tolerance = 100 * std::numeric_limits<Real<Cfg>>::epsilon();

  // the medium the ocean modes are computed for: c = 1.5, g = 9.81e-3
  constexpr double Rho = 1.0;
  constexpr double G = 9.81e-3;
  constexpr double Time = 0.3;

  auto material = fluid<Cfg>(Rho, 2.25);
  const CellMaterialData materialData{&material, {}};

  const physics::Ocean ocean(1, G, configIdOf<Cfg>());

  // no component vanishes here: x and y are no multiples of 5, kStar z is not near a multiple of
  // pi / 2
  const std::array<std::array<double, 3>, 3> points{
      {{1.3, 2.7, -0.4}, {3.1, 7.9, -1.2}, {8.2, 4.4, -2.5}}};
  auto buffer = guardedBuffer<Cfg>(points.size());
  auto dofs = dofsView<Cfg>(buffer, points.size());
  ocean.evaluate(Time, points.data(), points.size(), materialData, dofs);

  // mode 1 of the scenario, cf. Ocean::evaluateIn
  const double kX = std::acos(-1.0) / 10.0;
  const double kY = kX;
  constexpr double KStar = 1.5733628061766445;
  constexpr double Omega = 2.4523337594491745;
  constexpr double B = G * KStar / (Omega * Omega);
  for (std::size_t i = 0; i < points.size(); ++i) {
    const auto [x, y, z] = points[i];
    const double depth = std::sin(KStar * z) + B * std::cos(KStar * z);
    const double stress = -std::sin(kX * x) * std::sin(kY * y) * std::sin(Omega * Time) * depth;
    const std::array<double, 3> velocity{
        kX / (Omega * Rho) * std::cos(kX * x) * std::sin(kY * y) * std::cos(Omega * Time) * depth,
        kY / (Omega * Rho) * std::sin(kX * x) * std::cos(kY * y) * std::cos(Omega * Time) * depth,
        KStar / (Omega * Rho) * std::sin(kX * x) * std::sin(kY * y) * std::cos(Omega * Time) *
            (std::cos(KStar * z) - B * std::sin(KStar * z))};
    checkFluidState<Cfg>(dofs, i, stress, velocity, Tolerance);
  }
  checkGuard<Cfg>(buffer, points.size());
}

TEST_CASE_TEMPLATE_DEFINE("SuperimposedPlanarwave is the sum of its three planar waves" *
                              doctest::test_suite("physics"),
                          Cfg,
                          SuperimposedPlanarwaveSum) {
  using namespace scenarios;
  using real = Real<Cfg>;
  constexpr auto NumQuantities = model::MaterialOf<Cfg>::NumQuantities;
  constexpr double Tolerance = 100 * std::numeric_limits<real>::epsilon();

  model::MaterialOf<Cfg> material;
  material.rho = 1.0;
  material.setLameParameters(1.0, 2.0);
  const CellMaterialData materialData{&material, {}};

  // as many points as the initial field projection evaluates, i.e. more than basis functions
  constexpr std::size_t PointCount =
      (Cfg::ConvergenceOrder + 1) * (Cfg::ConvergenceOrder + 1) * (Cfg::ConvergenceOrder + 1);
  std::vector<std::array<double, 3>> points(PointCount);
  for (std::size_t i = 0; i < PointCount; ++i) {
    const auto s = static_cast<double>(i);
    points[i] = {std::fmod(0.37 * s, 1.0), std::fmod(0.61 * s, 1.0), std::fmod(0.83 * s, 1.0)};
  }

  constexpr double Phase = 0.5;
  constexpr double Time = 0.25;

  std::vector<real> actual(PointCount * NumQuantities);
  auto actualDofs = dofsView<Cfg>(actual, PointCount);
  const physics::SuperimposedPlanarwave superimposed(materialData, configIdOf<Cfg>(), Phase);
  superimposed.evaluate(Time, points.data(), PointCount, materialData, actualDofs);

  std::vector<real> expected(PointCount * NumQuantities, 0);
  std::vector<real> single(PointCount * NumQuantities);
  auto expectedDofs = dofsView<Cfg>(expected, PointCount);
  auto singleDofs = dofsView<Cfg>(single, PointCount);
  const std::array<Eigen::Vector3d, 3> kVecs{Eigen::Vector3d(M_PI, 0.0, 0.0),
                                             Eigen::Vector3d(0.0, M_PI, 0.0),
                                             Eigen::Vector3d(0.0, 0.0, M_PI)};
  for (const auto& kVec : kVecs) {
    const physics::Planarwave planarwave(materialData, configIdOf<Cfg>(), Phase, kVec);
    planarwave.evaluate(Time, points.data(), PointCount, materialData, singleDofs);
    for (std::size_t j = 0; j < NumQuantities; ++j) {
      for (std::size_t i = 0; i < PointCount; ++i) {
        expectedDofs(i, j) += singleDofs(i, j);
      }
    }
  }

  std::size_t mismatches = 0;
  for (std::size_t i = 0; i < actual.size(); ++i) {
    if (!std::isfinite(actual[i]) || std::abs(actual[i] - expected[i]) > Tolerance) {
      ++mismatches;
    }
  }
  CHECK(mismatches == 0);
}

TEST_CASE_TEMPLATE("Scenarios are offered only for materials that can represent them" *
                       doctest::test_suite("physics"),
                   Cfg,
                   SEISSOL_CONFIG_TYPES) {
  using initializer::parameters::InitializationType;
  const auto available = [](InitializationType type) {
    return physics::scenario::availability(type, configIdOf<Cfg>()).available;
  };

  constexpr auto Type = Cfg::MaterialType;
  constexpr bool Elastic = Type == model::MaterialType::Elastic;
  constexpr bool Fluid = Elastic || Type == model::MaterialType::Acoustic;

  // built from the material of the configuration itself
  CHECK(available(InitializationType::Zero));
  CHECK(available(InitializationType::Planarwave));
  CHECK(available(InitializationType::SuperimposedPlanarwave));
  CHECK(available(InitializationType::Easi));
  CHECK(available(InitializationType::Travelling) == (Cfg::RelaxationMechanisms == 0));

  // a solid next to a fluid
  CHECK(available(InitializationType::Scholte) == Elastic);
  CHECK(available(InitializationType::Snell) == Elastic);

  // a fluid
  CHECK(available(InitializationType::AcousticTravellingWithITM) == Fluid);
  CHECK(available(InitializationType::Ocean0) == Fluid);
  CHECK(available(InitializationType::Ocean1) == Fluid);
  CHECK(available(InitializationType::Ocean2) == Fluid);

  CHECK(available(InitializationType::PressureInjection) ==
        (Type == model::MaterialType::Poroelastic));
}

TEST_CASE("The acoustic travelling wave passes the time mirror only if it is enabled" *
          doctest::test_suite("physics")) {
  using namespace seissol::initializer::parameters;

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
    // the configuration only gives the default of the amplitudes, which the wave does not read
    const InitializationParameters initialization = readInitializationParameters(&reader, 0);
    return physics::scenario::acousticTravellingWaveParameters(initialization, itm);
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

TEST_CASE_TEMPLATE_APPLY(TravellingWaveItmPlaneWave, ConfigsWhere<scenarios::HoldsFluid>);
TEST_CASE_TEMPLATE_APPLY(TravellingWaveItmMirror, ConfigsWhere<scenarios::HoldsFluid>);
TEST_CASE_TEMPLATE_APPLY(OceanLayout, ConfigsWhere<scenarios::HoldsFluid>);
TEST_CASE_TEMPLATE_APPLY(SuperimposedPlanarwaveSum, ConfigsWhere<scenarios::SetUpByLame>);

} // namespace seissol::unit_test
