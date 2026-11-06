// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "DynamicRupture/Output/Geometry.h"
#include "DynamicRupture/Typedefs.h"

#include <array>
#include <cmath>
#include <cstddef>

namespace seissol::unit_test {
using namespace seissol::dr;

// ---------------------------------------------------------------------------
// ImpedancesAndEta: physical consistency
// ---------------------------------------------------------------------------

TEST_CASE("ImpedancesAndEta default zero" * doctest::test_suite("dynamicrupture")) {
  ImpedancesAndEta imp{};
  CHECK(imp.zp == doctest::Approx(0.0));
  CHECK(imp.zs == doctest::Approx(0.0));
  CHECK(imp.zpNeig == doctest::Approx(0.0));
  CHECK(imp.zsNeig == doctest::Approx(0.0));
  CHECK(imp.etaP == doctest::Approx(0.0));
  CHECK(imp.etaS == doctest::Approx(0.0));
}

TEST_CASE("ImpedancesAndEta physical setup" * doctest::test_suite("dynamicrupture")) {
  // Typical setup from the existing FrictionSolverCommon test
  ImpedancesAndEta imp;
  imp.zp = 10.0;
  imp.zs = 20.0;
  imp.zpNeig = 15.0;
  imp.zsNeig = 25.0;
  imp.etaP = imp.zp * imp.zpNeig / (imp.zp + imp.zpNeig);
  imp.etaS = imp.zs * imp.zsNeig / (imp.zs + imp.zsNeig);
  imp.invZp = 1.0 / imp.zp;
  imp.invZs = 1.0 / imp.zs;
  imp.invZpNeig = 1.0 / imp.zpNeig;
  imp.invZsNeig = 1.0 / imp.zsNeig;

  SUBCASE("Eta is harmonic mean") {
    // etaP = zp * zpNeig / (zp + zpNeig) = 10*15/25 = 6
    CHECK(imp.etaP == doctest::Approx(6.0));
    // etaS = zs * zsNeig / (zs + zsNeig) = 20*25/45 = 500/45 ≈ 11.111
    CHECK(imp.etaS == doctest::Approx(500.0 / 45.0));
  }

  SUBCASE("Inverses are correct") {
    CHECK(imp.invZp == doctest::Approx(0.1));
    CHECK(imp.invZs == doctest::Approx(0.05));
    CHECK(imp.invZpNeig == doctest::Approx(1.0 / 15.0));
    CHECK(imp.invZsNeig == doctest::Approx(0.04));
  }

  SUBCASE("Eta between the two impedances") {
    // Harmonic mean is always <= arithmetic mean
    CHECK(imp.etaP <= (imp.zp + imp.zpNeig) / 2.0);
    CHECK(imp.etaS <= (imp.zs + imp.zsNeig) / 2.0);
    CHECK(imp.etaP > 0.0);
    CHECK(imp.etaS > 0.0);
  }
}

TEST_CASE("ImpedancesAndEta equal impedances" * doctest::test_suite("dynamicrupture")) {
  ImpedancesAndEta imp;
  imp.zp = 10.0;
  imp.zpNeig = 10.0;
  imp.etaP = imp.zp * imp.zpNeig / (imp.zp + imp.zpNeig);

  // Harmonic mean of equal values = value / 2
  CHECK(imp.etaP == doctest::Approx(5.0));
}

// ---------------------------------------------------------------------------
// Receiver: default initialization
// ---------------------------------------------------------------------------

TEST_CASE("Receiver defaults" * doctest::test_suite("dynamicrupture")) {
  Receiver rp;
  CHECK_FALSE(rp.faultFaceIndex.hasValue());
  CHECK_FALSE(rp.localFaceSideId.hasValue());
  CHECK_FALSE(rp.localNeighborFaceSideId.hasValue());
  CHECK_FALSE(rp.elementIndex.hasValue());
  CHECK_FALSE(rp.elementGlobalIndex.hasValue());
  CHECK_FALSE(rp.elementNeighborGlobalIndex.hasValue());
  CHECK_FALSE(rp.globalReceiverIndex.hasValue());
  CHECK(rp.isInside == false);
  CHECK(rp.nearestGpIndex == -1);
  CHECK(rp.faultTag == -1);
  CHECK(rp.simIndex == 0);
}

// ---------------------------------------------------------------------------
// FaultDirections: default zero
// ---------------------------------------------------------------------------

TEST_CASE("FaultDirections defaults" * doctest::test_suite("dynamicrupture")) {
  FaultDirections fd;
  for (int i = 0; i < 3; ++i) {
    CHECK(fd.faceNormal[i] == doctest::Approx(0.0));
    CHECK(fd.tangent1[i] == doctest::Approx(0.0));
    CHECK(fd.tangent2[i] == doctest::Approx(0.0));
    CHECK(fd.strike[i] == doctest::Approx(0.0));
    CHECK(fd.dip[i] == doctest::Approx(0.0));
  }
}

// ---------------------------------------------------------------------------
// The impedance of a node of a damaged material
// ---------------------------------------------------------------------------

namespace {
/// Both sides of a face alike, without a background strain.
NodalImpedanceParameters damagedNodeParameters(double rho,
                                               double lambda0,
                                               double mu0,
                                               double gammaR,
                                               double xi0,
                                               const std::array<double, 4>& granular) {
  NodalImpedanceParameters params{};
  params.rhoPlus = params.rhoMinus = static_cast<real>(rho);
  params.lambda0Plus = params.lambda0Minus = static_cast<real>(lambda0);
  params.mu0Plus = params.mu0Minus = static_cast<real>(mu0);
  params.gammaRPlus = params.gammaRMinus = static_cast<real>(gammaR);
  params.xi0Plus = params.xi0Minus = static_cast<real>(xi0);
  for (std::size_t k = 0; k < 4; ++k) {
    params.granularPlus[k] = params.granularMinus[k] = static_cast<real>(granular[k]);
  }
  return params;
}

/// Checks the admittance, eta and lateral map of a node whose tangent is the
/// isotropic pair (lambda, mu), with the same material on both sides.
void checkIsotropicNode(
    const NodalImpedanceParameters& params, real alpha, real breakage, double lambda, double mu) {
  const real strain[6] = {static_cast<real>(-2.0e-4),
                          static_cast<real>(1.0e-4),
                          static_cast<real>(-0.5e-4),
                          static_cast<real>(3.0e-4),
                          static_cast<real>(-1.0e-4),
                          static_cast<real>(2.0e-4)};
  real admittancePlus[9];
  real admittanceMinus[9];
  real eta[9];
  real lateralStress[9];
  nodalAdmittance(params,
                  alpha,
                  breakage,
                  strain,
                  alpha,
                  breakage,
                  strain,
                  admittancePlus,
                  admittanceMinus,
                  eta,
                  lateralStress);

  const double rho = params.rhoPlus;
  const double zp = std::sqrt(rho * (lambda + 2.0 * mu));
  const double zs = std::sqrt(rho * mu);
  const double z[3] = {zp, zs, zs};
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t k = 0; k < 3; ++k) {
      const double diagonal = i == k ? 1.0 : 0.0;
      CHECK(admittancePlus[3 * i + k] * z[i] == doctest::Approx(diagonal));
      CHECK(admittanceMinus[3 * i + k] * z[i] == doctest::Approx(diagonal));
      CHECK(eta[3 * i + k] / z[i] == doctest::Approx(0.5 * diagonal));
    }
  }
  // Lateral rows yy, zz, yz against the traction columns nn, nt1, nt2: only the
  // two normal rows read the normal traction, with lambda / (lambda + 2 mu).
  const double ratio = lambda / (lambda + 2.0 * mu);
  const double lateral[9] = {ratio, 0.0, 0.0, ratio, 0.0, 0.0, 0.0, 0.0, 0.0};
  for (std::size_t k = 0; k < 9; ++k) {
    CHECK(lateralStress[k] == doctest::Approx(lateral[k]));
  }
}
} // namespace

TEST_CASE("Nodal impedance of a damaged material" * doctest::test_suite("dynamicrupture")) {
  constexpr double Rho = 2700.0;
  constexpr double Lambda0 = 3.2e10;
  constexpr double Mu0 = 3.0e10;
  constexpr double GammaR = 3.5e10;
  constexpr double Xi0 = 0.8;
  constexpr double A0 = 1.2e10;
  constexpr double A2 = 0.9e10;

  SUBCASE("Undamaged and unbroken: the elastic impedances") {
    const auto params = damagedNodeParameters(Rho, Lambda0, Mu0, GammaR, Xi0, {A0, 0.0, A2, 0.0});
    checkIsotropicNode(params, 0.0, 0.0, Lambda0, Mu0);
  }

  SUBCASE("Fully broken: the isotropic pair of the granular branch, whatever the damage") {
    // With aB1 = aB3 = 0 the granular branch is Hooke's law with mu = aB0 and
    // lambda = 2 aB2.
    const auto params = damagedNodeParameters(Rho, Lambda0, Mu0, GammaR, Xi0, {A0, 0.0, A2, 0.0});
    checkIsotropicNode(params, 0.3, 1.0, 2.0 * A2, A0);
  }

  SUBCASE("Half broken and undamaged: the mean of the two pairs") {
    const auto params = damagedNodeParameters(Rho, Lambda0, Mu0, GammaR, Xi0, {A0, 0.0, A2, 0.0});
    checkIsotropicNode(params, 0.0, 0.5, 0.5 * (Lambda0 + 2.0 * A2), 0.5 * (Mu0 + A0));
  }

  SUBCASE("Damaged and partly broken: the admittance is the inverse root") {
    // An anisotropic tangent, with every group of both branches in it.
    const auto params =
        damagedNodeParameters(Rho, Lambda0, Mu0, GammaR, Xi0, {A0, 0.2e10, A2, -0.1e10});
    const real alpha = 0.2;
    const real breakage = 0.3;
    const real strain[6] = {static_cast<real>(-2.0e-4),
                            static_cast<real>(1.0e-4),
                            static_cast<real>(-0.5e-4),
                            static_cast<real>(3.0e-4),
                            static_cast<real>(-1.0e-4),
                            static_cast<real>(2.0e-4)};

    real tangent[36];
    damageTangent(params.lambda0Plus,
                  params.mu0Plus,
                  params.gammaRPlus,
                  params.xi0Plus,
                  params.granularPlus,
                  alpha,
                  breakage,
                  strain,
                  tangent);
    for (std::size_t a = 0; a < 6; ++a) {
      for (std::size_t b = 0; b < 6; ++b) {
        CHECK(tangent[6 * a + b] == doctest::Approx(tangent[6 * b + a]));
      }
    }

    real admittancePlus[9];
    real admittanceMinus[9];
    real eta[9];
    real lateralStress[9];
    nodalAdmittance(params,
                    alpha,
                    breakage,
                    strain,
                    alpha,
                    breakage,
                    strain,
                    admittancePlus,
                    admittanceMinus,
                    eta,
                    lateralStress);

    // Y rho Gamma Y is the identity, with Gamma_ik = C_(in)(kn) and the
    // traction indices nn, nt1, nt2 of the face frame.
    constexpr std::size_t Traction[3] = {0, 3, 5};
    double scaled[9];
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t k = 0; k < 3; ++k) {
        double sum = 0.0;
        for (std::size_t j = 0; j < 3; ++j) {
          for (std::size_t l = 0; l < 3; ++l) {
            sum += admittancePlus[3 * i + j] * Rho * tangent[6 * Traction[j] + Traction[l]] *
                   admittancePlus[3 * l + k];
          }
        }
        scaled[3 * i + k] = sum;
      }
    }
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t k = 0; k < 3; ++k) {
        CHECK(scaled[3 * i + k] == doctest::Approx(i == k ? 1.0 : 0.0).epsilon(1e-4));
        CHECK(admittancePlus[3 * i + k] == doctest::Approx(admittancePlus[3 * k + i]));
      }
    }
    // the shear modes couple: the tangent is not isotropic here
    CHECK(std::abs(admittancePlus[3 * 1 + 2]) > 0.0);
  }
}

} // namespace seissol::unit_test
