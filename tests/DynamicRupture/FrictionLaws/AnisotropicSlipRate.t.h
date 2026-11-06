// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_DYNAMICRUPTURE_FRICTIONLAWS_ANISOTROPICSLIPRATE_T_H_
#define SEISSOL_TESTS_DYNAMICRUPTURE_FRICTIONLAWS_ANISOTROPICSLIPRATE_T_H_

// common::solveSlipRate projects the impedance only where it couples the modes of a face -- for an
// anisotropic MaterialT, and wherever the impedance belongs to the node; everywhere else every
// projection collapses to the scalar impedance and there is nothing to check here.

#include <doctest.h>

#include "DynamicRupture/FrictionLaws/FrictionSolverCommon.h"
#include "DynamicRupture/Typedefs.h"
#include "GeneratedCode/init.h"
#include "Kernels/Precision.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>

namespace seissol::unit_test {

namespace anisotropicsliprate {

using seissol::dr::FaultStresses;
using seissol::dr::ImpedanceMatrices;
using seissol::dr::ImpedancesAndEta;

/// eta = (Y+ + Y-)^-1 of a homogeneous fault in a VTI tilted out of the fault plane, rounded.
/// Symmetric positive definite, with both a shear/shear and a normal/shear coupling.
inline ImpedanceMatrices testImpedance() {
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
inline ImpedanceMatrices asymmetricImpedance() {
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

inline ImpedanceMatrices isotropicImpedance(real etaS) {
  ImpedanceMatrices impedanceMatrices;
  auto eta = init::eta::view::create(impedanceMatrices.eta);
  eta(0, 0) = 3.818e6;
  eta(1, 1) = etaS;
  eta(2, 2) = etaS;
  return impedanceMatrices;
}

/// The node the impedance is put at where it belongs to the node; any but the first, so that a
/// solve which ignores the index cannot pass.
constexpr std::uint32_t Node = 3;
static_assert(Node < seissol::dr::misc::NumPaddedPoints);

/// Puts eta where the material keeps it: with the face, or with the node of the stresses where the
/// impedance belongs to the node. Every place it does not belong holds a far softer impedance, so
/// that a solve reading the wrong one misses every check below. A template, so that the impedance
/// of the stresses is looked up once they are put in: they only have one where it belongs to the
/// node.
template <typename StressesT>
void placeEta(const ImpedanceMatrices& source,
              ImpedanceMatrices& face,
              [[maybe_unused]] StressesT& stresses) {
  if constexpr (seissol::dr::NodalImpedance) {
    const auto decoy = isotropicImpedance(1.0e3);
    face = decoy;
    for (std::size_t k = 0; k < 9; ++k) {
      for (std::size_t i = 0; i < seissol::dr::misc::NumPaddedPoints; ++i) {
        stresses.eta[k][i] = decoy.eta[k];
      }
      stresses.eta[k][Node] = source.eta[k];
    }
  } else {
    face = source;
  }
}

/// Residual of tau0 = (S I + V eta_ss) n with S = strength + slope * V * (eta n)_n, relative to
/// the trial traction. Zero for the exact solution, whatever route produced it.
inline real slipRateResidual(ImpedanceMatrices impedanceMatrices,
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

} // namespace anisotropicsliprate

// ---------------------------------------------------------------------------
// The solve the linear slip weakening laws and the fault receiver output share. It is checked
// against the equations it is supposed to satisfy, not against a second implementation of the
// same sweep, so a regression in either direction shows up here.
// ---------------------------------------------------------------------------
TEST_CASE("Slip rate solve of a coupled impedance" * doctest::skip(!seissol::dr::CoupledImpedance) *
          doctest::test_suite("dynamicrupture")) {
  using namespace anisotropicsliprate;

  using seissol::dr::friction_law::common::solveSlipRate;

  const ImpedancesAndEta impAndEta{};
  ImpedanceMatrices impedanceMatrices;
  FaultStresses<Executor::Host> faultStresses{};
  constexpr real Strength = 30.0e6;
  constexpr real Slope = 0.6;
  // two sweeps leave a truncation error below 1e-4 for an overshoot up to 100 percent, single
  // precision adds a few 1e-6 on top; dropping the normal coupling misses by 2e-2
  constexpr real ResidualBar = 5e-4;

  SUBCASE("solves its defining equations") {
    const auto source = testImpedance();
    placeEta(source, impedanceMatrices, faultStresses);

    for (const real overshoot : {0.02, 0.1, 0.5, 1.0}) {
      for (const real angle : {0.0, 0.7, 2.4, 4.9}) {
        const real magnitude = Strength * (1 + overshoot);
        const real traction1 = magnitude * std::cos(angle);
        const real traction2 = magnitude * std::sin(angle);

        const auto solution = solveSlipRate(impAndEta,
                                            impedanceMatrices,
                                            faultStresses,
                                            Node,
                                            traction1,
                                            traction2,
                                            magnitude,
                                            Strength,
                                            Slope);

        REQUIRE(solution.slipRate > 0);
        CHECK(std::abs(std::sqrt(solution.direction1 * solution.direction1 +
                                 solution.direction2 * solution.direction2) -
                       1) < 1e-6);
        CHECK(slipRateResidual(source, solution, traction1, traction2, Strength, Slope) <
              ResidualBar);
      }
    }
  }

  SUBCASE("the normal coupling is what closes the equations") {
    // dropping the strength slope is the state this solve replaces; it has to miss the residual
    // bar by a wide margin, otherwise the check above would pass for the wrong reason
    const auto source = testImpedance();
    placeEta(source, impedanceMatrices, faultStresses);
    const real traction1 = Strength * 1.5 * std::cos(0.7);
    const real traction2 = Strength * 1.5 * std::sin(0.7);

    const auto uncoupled = solveSlipRate(impAndEta,
                                         impedanceMatrices,
                                         faultStresses,
                                         Node,
                                         traction1,
                                         traction2,
                                         Strength * 1.5,
                                         Strength,
                                         0.0);

    CHECK(slipRateResidual(source, uncoupled, traction1, traction2, Strength, Slope) >
          10 * ResidualBar);
  }

  SUBCASE("isotropic impedance") {
    // no shear/shear and no normal/shear coupling: the slip is parallel to the trial traction and
    // the slip rate is the classical one, for any strength slope
    constexpr real EtaS = 2.2e6;
    placeEta(isotropicImpedance(EtaS), impedanceMatrices, faultStresses);
    const real magnitude = Strength * 1.5;
    const real angle = 0.7;
    const real traction1 = magnitude * std::cos(angle);
    const real traction2 = magnitude * std::sin(angle);

    const auto solution = solveSlipRate(impAndEta,
                                        impedanceMatrices,
                                        faultStresses,
                                        Node,
                                        traction1,
                                        traction2,
                                        magnitude,
                                        Strength,
                                        Slope);

    CHECK(solution.slipRate == doctest::Approx((magnitude - Strength) / EtaS).epsilon(1e-5));
    CHECK(solution.direction1 == doctest::Approx(std::cos(angle)).epsilon(1e-6));
    CHECK(solution.direction2 == doctest::Approx(std::sin(angle)).epsilon(1e-6));
    CHECK(solution.etaEff == doctest::Approx(EtaS).epsilon(1e-6));
  }

  SUBCASE("locked fault") {
    placeEta(testImpedance(), impedanceMatrices, faultStresses);
    const real magnitude = Strength * 0.5;

    const auto solution = solveSlipRate(impAndEta,
                                        impedanceMatrices,
                                        faultStresses,
                                        Node,
                                        magnitude * std::cos(0.7),
                                        magnitude * std::sin(0.7),
                                        magnitude,
                                        Strength,
                                        Slope);

    CHECK(solution.slipRate == static_cast<real>(0.0));
  }

  SUBCASE("vanishing trial traction") {
    placeEta(testImpedance(), impedanceMatrices, faultStresses);

    const auto solution = solveSlipRate(
        impAndEta, impedanceMatrices, faultStresses, Node, 0.0, 0.0, 0.0, Strength, Slope);

    CHECK(solution.slipRate == static_cast<real>(0.0));
    CHECK(std::isfinite(solution.direction1));
    CHECK(std::isfinite(solution.direction2));
  }
}

// ---------------------------------------------------------------------------
// The projections of eta the friction laws use. They index a flat, column-major array, so the
// checks below state which index is the row and which the column.
// ---------------------------------------------------------------------------
TEST_CASE("Projections of a coupled impedance" * doctest::skip(!seissol::dr::CoupledImpedance) *
          doctest::test_suite("dynamicrupture")) {
  using namespace anisotropicsliprate;

  namespace common = seissol::dr::friction_law::common;

  const ImpedancesAndEta impAndEta{};
  auto source = asymmetricImpedance();
  const auto eta = init::eta::view::create(source.eta);
  ImpedanceMatrices impedanceMatrices;
  FaultStresses<Executor::Host> faultStresses{};
  placeEta(source, impedanceMatrices, faultStresses);

  constexpr real V1 = 0.37;
  constexpr real V2 = -0.91;
  const real magnitude = std::sqrt(V1 * V1 + V2 * V2);
  const real n1 = V1 / magnitude;
  const real n2 = V2 / magnitude;

  SUBCASE("matmulEta applies eta, not its transpose") {
    const auto [w1, w2] =
        common::matmulEta(impAndEta, impedanceMatrices, faultStresses, Node, V1, V2);

    CHECK(w1 == doctest::Approx(eta(1, 1) * V1 + eta(1, 2) * V2).epsilon(1e-5));
    CHECK(w2 == doctest::Approx(eta(2, 1) * V1 + eta(2, 2) * V2).epsilon(1e-5));
  }

  SUBCASE("the normal coupling reads the fault-normal row") {
    const auto wn =
        common::matmulEtaNormal(impAndEta, impedanceMatrices, faultStresses, Node, V1, V2);
    const auto cn = common::projectEtaNormal(
        impAndEta, impedanceMatrices, faultStresses, Node, V1, V2, magnitude);

    CHECK(wn == doctest::Approx(eta(0, 1) * V1 + eta(0, 2) * V2).epsilon(1e-5));
    CHECK(cn == doctest::Approx(eta(0, 1) * n1 + eta(0, 2) * n2).epsilon(1e-5));
  }

  SUBCASE("the slip direction solves with the shear block") {
    // the adjugate route of updateSlipDirection, checked against the equation it solves:
    // (S I + V eta_ss) n has to point along the trial traction
    constexpr real Strength = 2.0e6;
    constexpr real SlipRate = 1.3;
    const auto [d1, d2] = common::updateSlipDirection(
        impAndEta, impedanceMatrices, faultStresses, Node, Strength, SlipRate, V1, V2, magnitude);
    const real image1 = Strength * d1 + SlipRate * (eta(1, 1) * d1 + eta(1, 2) * d2);
    const real image2 = Strength * d2 + SlipRate * (eta(2, 1) * d1 + eta(2, 2) * d2);

    CHECK(std::sqrt(d1 * d1 + d2 * d2) == doctest::Approx(1.0).epsilon(1e-6));
    CHECK(image1 * V2 - image2 * V1 == doctest::Approx(0.0).scale(magnitude * 1.0e6));
    CHECK(image1 * V1 + image2 * V2 > 0);
  }

  SUBCASE("projectEta is the quadratic form of the shear block") {
    const auto [etaProj, invEtaProj] =
        common::projectEta(impAndEta, impedanceMatrices, faultStresses, Node, V1, V2, magnitude);

    const real expected =
        eta(1, 1) * n1 * n1 + (eta(1, 2) + eta(2, 1)) * n1 * n2 + eta(2, 2) * n2 * n2;
    CHECK(etaProj == doctest::Approx(expected).epsilon(1e-5));
    CHECK(invEtaProj == doctest::Approx(1.0 / expected).epsilon(1e-5));
  }
}

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_DYNAMICRUPTURE_FRICTIONLAWS_ANISOTROPICSLIPRATE_T_H_
