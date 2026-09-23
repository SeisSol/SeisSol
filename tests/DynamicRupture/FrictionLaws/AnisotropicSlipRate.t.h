// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_DYNAMICRUPTURE_FRICTIONLAWS_ANISOTROPICSLIPRATE_T_H_
#define SEISSOL_TESTS_DYNAMICRUPTURE_FRICTIONLAWS_ANISOTROPICSLIPRATE_T_H_

// common::solveSlipRate projects the impedance only for an anisotropic MaterialT; everywhere else
// every projection collapses to the scalar impedance and there is nothing to check here.

#include <doctest.h>

#include "DynamicRupture/FrictionLaws/FrictionSolverCommon.h"
#include "DynamicRupture/Typedefs.h"
#include "Equations/Datastructures.h"
#include "GeneratedCode/init.h"
#include "Kernels/Precision.h"

#include <array>
#include <cmath>
#include <type_traits>

namespace seissol::unit_test {

namespace {

using seissol::dr::ImpedanceMatrices;
using seissol::dr::ImpedancesAndEta;

/// eta = (Y+ + Y-)^-1 of a homogeneous fault in a VTI tilted out of the fault plane, rounded.
/// Symmetric positive definite, with both a shear/shear and a normal/shear coupling.
ImpedanceMatrices testImpedance() {
  ImpedanceMatrices impedanceMatrices;
  auto eta = init::eta::view::create(impedanceMatrices.eta);
  eta(0, 0) = 3.818e6;
  eta(1, 1) = 2.106e6;
  eta(2, 2) = 2.337e6;
  eta(0, 1) = eta(1, 0) = 2.11e5;
  eta(0, 2) = eta(2, 0) = 1.9e4;
  eta(1, 2) = eta(2, 1) = -2.1e4;
  return impedanceMatrices;
}

/// eta with a deliberately asymmetric shear block. Physical impedances are self-adjoint, which
/// makes eta and its transpose interchangeable -- this one tells them apart.
ImpedanceMatrices asymmetricImpedance() {
  ImpedanceMatrices impedanceMatrices;
  auto eta = init::eta::view::create(impedanceMatrices.eta);
  eta(0, 0) = 3.8e6;
  eta(0, 1) = 1.1e5;
  eta(0, 2) = 2.2e5;
  eta(1, 0) = 3.3e5;
  eta(1, 1) = 2.1e6;
  eta(1, 2) = 4.4e5;
  eta(2, 0) = 5.5e5;
  eta(2, 1) = 6.6e5;
  eta(2, 2) = 2.3e6;
  return impedanceMatrices;
}

ImpedanceMatrices isotropicImpedance(real etaS) {
  ImpedanceMatrices impedanceMatrices;
  auto eta = init::eta::view::create(impedanceMatrices.eta);
  eta(0, 0) = 3.818e6;
  eta(1, 1) = etaS;
  eta(2, 2) = etaS;
  return impedanceMatrices;
}

/// Residual of tau0 = (S I + V eta_ss) n with S = strength + slope * V * (eta n)_n, relative to
/// the trial traction. Zero for the exact solution, whatever route produced it.
real slipRateResidual(ImpedanceMatrices impedanceMatrices,
                      const seissol::dr::friction_law::common::SlipRateSolution& solution,
                      real traction1,
                      real traction2,
                      real strength,
                      real strengthSlope) {
  const auto eta = init::eta::view::create(impedanceMatrices.eta);
  const real slip1 = solution.slipRate * solution.direction1;
  const real slip2 = solution.slipRate * solution.direction2;

  const real normalTraction = eta(0, 1) * slip1 + eta(0, 2) * slip2;
  const real localStrength = strength + strengthSlope * normalTraction;

  const real residual1 =
      traction1 - (eta(1, 1) * slip1 + eta(1, 2) * slip2) - localStrength * solution.direction1;
  const real residual2 =
      traction2 - (eta(2, 1) * slip1 + eta(2, 2) * slip2) - localStrength * solution.direction2;

  return std::sqrt(residual1 * residual1 + residual2 * residual2) /
         std::sqrt(traction1 * traction1 + traction2 * traction2);
}

} // namespace

// ---------------------------------------------------------------------------
// The solve the linear slip weakening laws and the fault receiver output share. It is checked
// against the equations it is supposed to satisfy, not against a second implementation of the
// same sweep, so a regression in either direction shows up here.
// ---------------------------------------------------------------------------
TEST_CASE("Anisotropic slip rate solve" *
          doctest::skip(!std::is_same_v<model::MaterialT, model::AnisotropicMaterial>) *
          doctest::test_suite("dynamicrupture")) {
  using seissol::dr::friction_law::common::solveSlipRate;

  const ImpedancesAndEta impAndEta{};
  constexpr real Strength = 30.0e6;
  constexpr real Slope = 0.6;
  // two sweeps leave a truncation error below 1e-4 for an overshoot up to 100 percent, single
  // precision adds a few 1e-6 on top; dropping the normal coupling misses by 2e-2
  constexpr real ResidualBar = 5e-4;

  SUBCASE("solves its defining equations") {
    const auto impedanceMatrices = testImpedance();

    for (const real overshoot : {0.02, 0.1, 0.5, 1.0}) {
      for (const real angle : {0.0, 0.7, 2.4, 4.9}) {
        const real magnitude = Strength * (1 + overshoot);
        const real traction1 = magnitude * std::cos(angle);
        const real traction2 = magnitude * std::sin(angle);

        const auto solution = solveSlipRate(
            impAndEta, impedanceMatrices, traction1, traction2, magnitude, Strength, Slope);

        REQUIRE(solution.slipRate > 0);
        CHECK(std::abs(std::sqrt(solution.direction1 * solution.direction1 +
                                 solution.direction2 * solution.direction2) -
                       1) < 1e-6);
        CHECK(slipRateResidual(impedanceMatrices, solution, traction1, traction2, Strength, Slope) <
              ResidualBar);
      }
    }
  }

  SUBCASE("the normal coupling is what closes the equations") {
    // dropping the strength slope is the state this solve replaces; it has to miss the residual
    // bar by a wide margin, otherwise the check above would pass for the wrong reason
    const auto impedanceMatrices = testImpedance();
    const real traction1 = Strength * 1.5 * std::cos(0.7);
    const real traction2 = Strength * 1.5 * std::sin(0.7);

    const auto uncoupled = solveSlipRate(
        impAndEta, impedanceMatrices, traction1, traction2, Strength * 1.5, Strength, 0.0);

    CHECK(slipRateResidual(impedanceMatrices, uncoupled, traction1, traction2, Strength, Slope) >
          10 * ResidualBar);
  }

  SUBCASE("isotropic impedance") {
    // no shear/shear and no normal/shear coupling: the slip is parallel to the trial traction and
    // the slip rate is the classical one, for any strength slope
    constexpr real EtaS = 2.2e6;
    const auto impedanceMatrices = isotropicImpedance(EtaS);
    const real magnitude = Strength * 1.5;
    const real angle = 0.7;
    const real traction1 = magnitude * std::cos(angle);
    const real traction2 = magnitude * std::sin(angle);

    const auto solution = solveSlipRate(
        impAndEta, impedanceMatrices, traction1, traction2, magnitude, Strength, Slope);

    CHECK(solution.slipRate == doctest::Approx((magnitude - Strength) / EtaS).epsilon(1e-5));
    CHECK(solution.direction1 == doctest::Approx(std::cos(angle)).epsilon(1e-6));
    CHECK(solution.direction2 == doctest::Approx(std::sin(angle)).epsilon(1e-6));
    CHECK(solution.etaEff == doctest::Approx(EtaS).epsilon(1e-6));
  }

  SUBCASE("locked fault") {
    const auto impedanceMatrices = testImpedance();
    const real magnitude = Strength * 0.5;

    const auto solution = solveSlipRate(impAndEta,
                                        impedanceMatrices,
                                        magnitude * std::cos(0.7),
                                        magnitude * std::sin(0.7),
                                        magnitude,
                                        Strength,
                                        Slope);

    CHECK(solution.slipRate == static_cast<real>(0.0));
  }

  SUBCASE("vanishing trial traction") {
    const auto impedanceMatrices = testImpedance();

    const auto solution =
        solveSlipRate(impAndEta, impedanceMatrices, 0.0, 0.0, 0.0, Strength, Slope);

    CHECK(solution.slipRate == static_cast<real>(0.0));
    CHECK(std::isfinite(solution.direction1));
    CHECK(std::isfinite(solution.direction2));
  }
}

// ---------------------------------------------------------------------------
// The projections of eta the friction laws use. They index a flat, column-major array, so the
// checks below state which index is the row and which the column.
// ---------------------------------------------------------------------------
TEST_CASE("Anisotropic impedance projections" *
          doctest::skip(!std::is_same_v<model::MaterialT, model::AnisotropicMaterial>) *
          doctest::test_suite("dynamicrupture")) {
  namespace common = seissol::dr::friction_law::common;

  const ImpedancesAndEta impAndEta{};
  auto impedanceMatrices = asymmetricImpedance();
  const auto eta = init::eta::view::create(impedanceMatrices.eta);

  constexpr real V1 = 0.37;
  constexpr real V2 = -0.91;
  const real magnitude = std::sqrt(V1 * V1 + V2 * V2);
  const real n1 = V1 / magnitude;
  const real n2 = V2 / magnitude;

  SUBCASE("matmulEta applies eta, not its transpose") {
    const auto [w1, w2] = common::matmulEta(impAndEta, impedanceMatrices, V1, V2);

    CHECK(w1 == doctest::Approx(eta(1, 1) * V1 + eta(1, 2) * V2).epsilon(1e-5));
    CHECK(w2 == doctest::Approx(eta(2, 1) * V1 + eta(2, 2) * V2).epsilon(1e-5));
  }

  SUBCASE("the normal coupling reads the fault-normal row") {
    const auto wn = common::matmulEtaNormal(impAndEta, impedanceMatrices, V1, V2);

    CHECK(wn == doctest::Approx(eta(0, 1) * V1 + eta(0, 2) * V2).epsilon(1e-5));
  }

  SUBCASE("projectEta is the quadratic form of the shear block") {
    const auto [etaProj, invEtaProj] =
        common::projectEta(impAndEta, impedanceMatrices, V1, V2, magnitude);

    const real expected =
        eta(1, 1) * n1 * n1 + (eta(1, 2) + eta(2, 1)) * n1 * n2 + eta(2, 2) * n2 * n2;
    CHECK(etaProj == doctest::Approx(expected).epsilon(1e-5));
    CHECK(invEtaProj == doctest::Approx(1.0 / expected).epsilon(1e-5));
  }
}

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_DYNAMICRUPTURE_FRICTIONLAWS_ANISOTROPICSLIPRATE_T_H_
