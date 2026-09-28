// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/Constants.h"
#include "Equations/Datastructures.h"
#include "Initializer/Parameters/InitializationParameters.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"
#include "Model/CommonDatastructures.h"
#include "Physics/InitialField.h"
#include "TestHelper.h"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <math.h>
#include <vector>
#include <yateto.h>

namespace seissol::unit_test {

namespace initialfield {

constexpr auto Type = model::MaterialT::Type;
constexpr std::size_t NumQuantities = model::MaterialT::NumQuantities;

/// The builds that hold a fluid and whose material is set up from rho, lambda and mu alone, i.e.
/// the builds the fluid scenarios are written for. Elsewhere, their value tests are skipped (the
/// anisotropic setLameParameters, for example, calls logError).
constexpr bool FluidBuild =
    Type == model::MaterialType::Elastic || Type == model::MaterialType::Acoustic;

/// The builds whose material is set up from rho, lambda and mu alone. The attenuating ones then
/// carry no relaxation (omega = theta = 0), which leaves their memory variables decoupled.
constexpr bool LameBuild = FluidBuild || Type == model::MaterialType::Viscoelastic ||
                           Type == model::MaterialType::Viscoacoustic;

/// Absolute tolerance for values of order one, computed in double and stored in real.
constexpr double Tolerance = 100 * std::numeric_limits<real>::epsilon();

/// Fills every entry before an initial condition is evaluated.
constexpr real Sentinel = -1234.5;

/// Values for `pointCount` points of the build, followed by a guard as large as the elastic layout
/// (nine quantities) would need in addition. Nothing may be written into the guard.
inline std::vector<real> guardedBuffer(std::size_t pointCount) {
  return std::vector<real>(pointCount * (NumQuantities + 9), Sentinel);
}

inline yateto::DenseTensorView<2, real, unsigned> dofsView(std::vector<real>& buffer,
                                                           std::size_t pointCount) {
  return yateto::DenseTensorView<2, real, unsigned>(
      buffer.data(), {static_cast<unsigned>(pointCount), static_cast<unsigned>(NumQuantities)});
}

inline void checkGuard(const std::vector<real>& buffer, std::size_t pointCount) {
  for (std::size_t i = pointCount * NumQuantities; i < buffer.size(); ++i) {
    CAPTURE(i);
    CHECK(buffer[i] == Sentinel);
  }
}

/// Checks the state of a fluid at `point`: `stress` (i.e. -p) on the normal stresses, no shear
/// stress, and `velocity` in the velocity of the build.
inline void checkFluidState(const yateto::DenseTensorView<2, real, unsigned>& dofs,
                            std::size_t point,
                            double stress,
                            const std::array<double, 3>& velocity,
                            double tolerance) {
  CAPTURE(point);
  // a stress tensor has three normal components; the acoustic equations carry a single one
  constexpr std::size_t NormalComponents =
      std::min<std::size_t>(3, model::MaterialT::TractionComponents);
  for (std::size_t j = 0; j < model::MaterialT::TractionComponents; ++j) {
    CAPTURE(j);
    CHECK(dofs(point, j) == AbsApprox(j < NormalComponents ? stress : 0.0).epsilon(tolerance));
  }
  for (std::size_t d = 0; d < velocity.size(); ++d) {
    CAPTURE(d);
    CHECK(dofs(point, model::MaterialT::VelocityOffset + d) ==
          AbsApprox(velocity[d]).epsilon(tolerance));
  }
}

} // namespace initialfield

TEST_CASE("AcousticTravellingWaveITM writes a plane wave into the pressure and velocity of the "
          "build" *
          doctest::skip(!initialfield::FluidBuild) * doctest::test_suite("physics")) {
  using namespace initialfield;

  constexpr double Rho = 2.0;
  constexpr double Lambda = 8.0;
  constexpr double K = 3.0;
  const double c = std::sqrt(Lambda / Rho);

  model::MaterialT material;
  material.rho = Rho;
  material.setLameParameters(0.0, Lambda);
  const CellMaterialData materialData{&material, {}};

  // With a factor of one, the mirror leaves the wave as it is. So before, during and after it, the
  // plain travelling wave p = rho c cos(k x - c k t), u = cos(k x - c k t) has to come out.
  AcousticTravellingWaveParametersITM parameters{};
  parameters.k = K;
  parameters.itmStartingTime = 1.0;
  parameters.itmDuration = 0.5;
  parameters.itmVelocityScalingFactor = 1.0;
  const physics::AcousticTravellingWaveITM wave(materialData, parameters);

  const std::array<std::array<double, 3>, 3> points{
      {{0.1, 0.2, 0.3}, {0.45, -0.7, 1.1}, {0.9, 0.5, -0.4}}};
  for (const double time : {0.5, 1.25, 2.0}) {
    CAPTURE(time);
    auto buffer = guardedBuffer(points.size());
    auto dofs = dofsView(buffer, points.size());
    wave.evaluate(time, points.data(), points.size(), materialData, dofs);
    for (std::size_t i = 0; i < points.size(); ++i) {
      const double phase = K * points[i][0] - c * K * time;
      checkFluidState(dofs, i, -Rho * c * std::cos(phase), {std::cos(phase), 0.0, 0.0}, Tolerance);
    }
    checkGuard(buffer, points.size());
  }
}

TEST_CASE("AcousticTravellingWaveITM keeps the pressure and scales the velocity at the mirror" *
          doctest::skip(!initialfield::FluidBuild) * doctest::test_suite("physics")) {
  using namespace initialfield;

  constexpr double Start = 1.0;
  constexpr double Duration = 0.5;
  constexpr double N = 2.0;
  // the values right after a switch are taken one double later than the switch itself
  constexpr double SwitchTolerance = std::max(1.0e-10, Tolerance);

  model::MaterialT material;
  material.rho = 2.0;
  material.setLameParameters(0.0, 8.0);
  const CellMaterialData materialData{&material, {}};

  AcousticTravellingWaveParametersITM parameters{};
  parameters.k = 3.0;
  parameters.itmStartingTime = Start;
  parameters.itmDuration = Duration;
  parameters.itmVelocityScalingFactor = N;
  const physics::AcousticTravellingWaveITM wave(materialData, parameters);

  const std::array<std::array<double, 3>, 2> points{{{0.1, 0.2, 0.3}, {0.9, 0.5, -0.4}}};
  const auto evaluateAt = [&](double time) {
    auto buffer = guardedBuffer(points.size());
    auto dofs = dofsView(buffer, points.size());
    wave.evaluate(time, points.data(), points.size(), materialData, dofs);
    checkGuard(buffer, points.size());
    return buffer;
  };
  // The mirror scales the velocity by the factor n when it starts, and by 1 / n when it ends; the
  // pressure passes both switches continuously.
  const auto checkSwitch = [&](double switchTime, double velocityFactor) {
    auto before = evaluateAt(switchTime);
    auto after = evaluateAt(std::nextafter(switchTime, std::numeric_limits<double>::infinity()));
    const auto dofsBefore = dofsView(before, points.size());
    const auto dofsAfter = dofsView(after, points.size());
    for (std::size_t i = 0; i < points.size(); ++i) {
      CAPTURE(i);
      CHECK(dofsAfter(i, 0) == AbsApprox(dofsBefore(i, 0)).epsilon(SwitchTolerance));
      CHECK(dofsAfter(i, model::MaterialT::VelocityOffset) ==
            AbsApprox(velocityFactor * dofsBefore(i, model::MaterialT::VelocityOffset))
                .epsilon(SwitchTolerance));
    }
  };

  SUBCASE("start") { checkSwitch(Start, N); }
  SUBCASE("end") { checkSwitch(Start + Duration, 1.0 / N); }
}

TEST_CASE("Ocean writes a fluid stress and the velocity of the build" *
          doctest::skip(!initialfield::FluidBuild) * doctest::test_suite("physics")) {
  using namespace initialfield;

  // the medium the ocean modes are computed for: c = 1.5, g = 9.81e-3
  constexpr double Rho = 1.0;
  constexpr double G = 9.81e-3;
  constexpr double Time = 0.3;

  model::MaterialT material;
  material.rho = Rho;
  material.setLameParameters(0.0, 2.25);
  const CellMaterialData materialData{&material, {}};

  const physics::Ocean ocean(1, G);

  // no component vanishes here: x and y are no multiples of 5, kStar z is not near a multiple of
  // pi / 2
  const std::array<std::array<double, 3>, 3> points{
      {{1.3, 2.7, -0.4}, {3.1, 7.9, -1.2}, {8.2, 4.4, -2.5}}};
  auto buffer = guardedBuffer(points.size());
  auto dofs = dofsView(buffer, points.size());
  ocean.evaluate(Time, points.data(), points.size(), materialData, dofs);

  // mode 1 of the scenario, cf. Ocean::evaluate
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
    checkFluidState(dofs, i, stress, velocity, Tolerance);
  }
  checkGuard(buffer, points.size());
}

TEST_CASE("SuperimposedPlanarwave is the sum of its three planar waves" *
          doctest::skip(!initialfield::LameBuild) * doctest::test_suite("physics")) {
  using namespace initialfield;

  model::MaterialT material;
  material.rho = 1.0;
  material.setLameParameters(1.0, 2.0);
  const CellMaterialData materialData{&material, {}};

  // as many points as the initial field projection evaluates, i.e. more than basis functions
  constexpr std::size_t PointCount =
      (ConvergenceOrder + 1) * (ConvergenceOrder + 1) * (ConvergenceOrder + 1);
  std::vector<std::array<double, 3>> points(PointCount);
  for (std::size_t i = 0; i < PointCount; ++i) {
    const auto s = static_cast<double>(i);
    points[i] = {std::fmod(0.37 * s, 1.0), std::fmod(0.61 * s, 1.0), std::fmod(0.83 * s, 1.0)};
  }

  constexpr real Phase = 0.5;
  constexpr double Time = 0.25;

  std::vector<real> actual(PointCount * NumQuantities);
  auto actualDofs = dofsView(actual, PointCount);
  const physics::SuperimposedPlanarwave superimposed(materialData, Phase);
  superimposed.evaluate(Time, points.data(), PointCount, materialData, actualDofs);

  std::vector<real> expected(PointCount * NumQuantities, 0);
  std::vector<real> single(PointCount * NumQuantities);
  auto expectedDofs = dofsView(expected, PointCount);
  auto singleDofs = dofsView(single, PointCount);
  const std::array<Eigen::Vector3d, 3> kVecs{Eigen::Vector3d(M_PI, 0.0, 0.0),
                                             Eigen::Vector3d(0.0, M_PI, 0.0),
                                             Eigen::Vector3d(0.0, 0.0, M_PI)};
  for (const auto& kVec : kVecs) {
    const physics::Planarwave planarwave(materialData, Phase, kVec);
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

TEST_CASE("Hard-coded initial conditions are offered only for equations that can represent them" *
          doctest::test_suite("physics")) {
  using initializer::parameters::InitializationType;
  using physics::isInitialConditionSupported;

  constexpr auto Type = model::MaterialT::Type;
  constexpr bool Elastic = Type == model::MaterialType::Elastic;
  constexpr bool Fluid = Elastic || Type == model::MaterialType::Acoustic;

  // built from the build's own equations
  CHECK(isInitialConditionSupported(InitializationType::Zero));
  CHECK(isInitialConditionSupported(InitializationType::Planarwave));
  CHECK(isInitialConditionSupported(InitializationType::SuperimposedPlanarwave));
  CHECK(isInitialConditionSupported(InitializationType::Easi));
  CHECK(isInitialConditionSupported(InitializationType::Travelling) ==
        (model::MaterialT::Mechanisms == 0));

  // a solid next to a fluid
  CHECK(isInitialConditionSupported(InitializationType::Scholte) == Elastic);
  CHECK(isInitialConditionSupported(InitializationType::Snell) == Elastic);

  // a fluid
  CHECK(isInitialConditionSupported(InitializationType::AcousticTravellingWithITM) == Fluid);
  CHECK(isInitialConditionSupported(InitializationType::Ocean0) == Fluid);
  CHECK(isInitialConditionSupported(InitializationType::Ocean1) == Fluid);
  CHECK(isInitialConditionSupported(InitializationType::Ocean2) == Fluid);

  CHECK(isInitialConditionSupported(InitializationType::PressureInjection) ==
        (Type == model::MaterialType::Poroelastic));
}

} // namespace seissol::unit_test
