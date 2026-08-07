// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_MODEL_POROELASTICIMPEDANCE_T_H_
#define SEISSOL_TESTS_MODEL_POROELASTICIMPEDANCE_T_H_

#ifdef USE_POROELASTIC

#include <doctest.h>

#include "Equations/Datastructures.h"
#include "Equations/Setup.h"
#include "Initializer/Model/DynamicRuptureImpedance.h"

#include <Eigen/Dense>
#include <cmath>
#include <vector>

namespace seissol::unit_test {

namespace {

/// bulkSolid, rho, lambda, mu, porosity, permeability, tortuosity, bulkFluid, rhoFluid, viscosity
model::PoroElasticMaterial testPoroMaterial(double porosity, double tortuosity) {
  return model::PoroElasticMaterial(std::vector<double>{
      3.60e10, 2650.0, 4.0e9, 6.0e9, porosity, 1.0e-13, tortuosity, 2.2e9, 1000.0, 1.0e-3});
}

} // namespace

// ---------------------------------------------------------------------------
// The poroelastic interface has four variables, (sigma_nn, sigma_ns, sigma_nd, p) against
// (v_n, v_s, v_d, q_n). SeisSol's poroelastic frame is isotropic, so the generalized wave
// impedance Z = Mass # Gamma has a closed form -- a 2x2 fast/slow P block plus two shear scalars.
// Checking it against the general eigendecomposition validates both routes against each other.
// ---------------------------------------------------------------------------
TEST_CASE("Poroelastic DR impedance closed form" * doctest::test_suite("dynamicrupture")) {
  using seissol::initializer::model::checkFaultImpedance;
  using seissol::initializer::model::computeAdmittance;
  using seissol::initializer::model::computeFaultImpedance;
  using seissol::initializer::model::DrLateralMatrix;
  using seissol::initializer::model::DrMatrix;

  // the eigendecomposition of the 13x13 Jacobian is comparatively ill conditioned here
  // (cond(R_t) reaches ~1e6 at high porosity), so this is the accuracy we can demand of it
  constexpr double Epsilon = 1e-7;

  for (const double porosity : {0.05, 0.2, 0.4}) {
    for (const double tortuosity : {1.0, 1.5, 3.0}) {
      const auto sweepMaterial = testPoroMaterial(porosity, tortuosity);

      DrLateralMatrix lateralClosed = DrLateralMatrix::Zero();
      const auto admittanceClosed = computeAdmittance(sweepMaterial, &lateralClosed);

      DrLateralMatrix lateralEigen = DrLateralMatrix::Zero();
      const auto admittanceEigen =
          seissol::initializer::model::impedance_detail::admittanceFromEigendecomposition(
              sweepMaterial, &lateralEigen);

      CHECK((admittanceClosed - admittanceEigen).cwiseAbs().maxCoeff() <
            Epsilon * admittanceEigen.cwiseAbs().maxCoeff());
      CHECK((lateralClosed - lateralEigen).cwiseAbs().maxCoeff() <
            Epsilon * lateralEigen.cwiseAbs().maxCoeff());
    }
  }

  const auto material = testPoroMaterial(0.2, 1.5);
  const auto impedance = computeFaultImpedance(material, material);

  SUBCASE("impedance invariants") {
    // self-adjointness now holds too: it is checked against the signature matrix, which accounts
    // for the fourth traction component being stored as +p while its conjugate partner is -p
    const auto violation = checkFaultImpedance(impedance);
    REQUIRE_MESSAGE(!violation.has_value(), violation.value_or(""));

    // ... and the raw admittance is genuinely *not* self-adjoint without that correction
    const DrMatrix raw = impedance.admittancePlus;
    CHECK((raw - raw.transpose()).cwiseAbs().maxCoeff() > 0.01 * raw.cwiseAbs().maxCoeff());

    // homogeneous fault: the traction is split evenly
    const DrMatrix half = DrMatrix::Identity() * 0.5;
    CHECK((impedance.bPlus - half).cwiseAbs().maxCoeff() < 1e-9);
  }

  SUBCASE("shear impedance equals the Biot value") {
    // Z_s = sqrt(mu * rho1) with the condensed density rho1 = rhoBar - rhoFluid^2 / m, i.e. the
    // static condensation of the tangential filtration velocity
    const double m = material.rhoFluid * material.tortuosity / material.porosity;
    const double rhoBar =
        (1 - material.porosity) * material.rho + material.porosity * material.rhoFluid;
    const double rho1 = rhoBar - material.rhoFluid * material.rhoFluid / m;
    const double shearImpedance = std::sqrt(material.mu * rho1);

    CHECK(impedance.admittancePlus(1, 1) == doctest::Approx(1.0 / shearImpedance).epsilon(1e-12));
    CHECK(impedance.admittancePlus(2, 2) == doctest::Approx(1.0 / shearImpedance).epsilon(1e-12));
  }

  SUBCASE("lateral stress structure of an isotropic solid frame") {
    const auto& lateral = impedance.lateralStressPlus;
    const double scale = lateral.cwiseAbs().maxCoeff();
    REQUIRE(scale > 0.0);

    // sigma_sd is not driven by anything travelling along the fault normal
    CHECK(lateral.row(2).cwiseAbs().maxCoeff() < 1e-12 * scale);
    // sigma_ss and sigma_dd react identically
    CHECK((lateral.row(0) - lateral.row(1)).cwiseAbs().maxCoeff() < 1e-12 * scale);
    // the two shear tractions do not produce any lateral normal stress
    CHECK(lateral.col(1).cwiseAbs().maxCoeff() < 1e-12 * scale);
    CHECK(lateral.col(2).cwiseAbs().maxCoeff() < 1e-12 * scale);
    // ... while sigma_nn and the pore pressure do
    CHECK(std::abs(lateral(0, 0)) > 1e-3 * scale);
    CHECK(std::abs(lateral(0, 3)) > 1e-3 * scale);
  }
}

} // namespace seissol::unit_test

#endif // USE_POROELASTIC

#endif // SEISSOL_TESTS_MODEL_POROELASTICIMPEDANCE_T_H_
