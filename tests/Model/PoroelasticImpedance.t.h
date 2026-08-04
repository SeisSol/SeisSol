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

model::PoroElasticMaterial testPoroMaterial() {
  // bulkSolid, rho, lambda, mu, porosity, permeability, tortuosity, bulkFluid, rhoFluid, viscosity
  return model::PoroElasticMaterial(
      std::vector<double>{3.60e10, 2650.0, 4.0e9, 6.0e9, 0.2, 1.0e-13, 1.5, 2.2e9, 1000.0, 1.0e-3});
}

} // namespace

// ---------------------------------------------------------------------------
// The poroelastic interface has four variables, (sigma_nn, sigma_ns, sigma_nd, p) against
// (v_n, v_s, v_d, q_n). The reconstruction of the remaining stress components therefore is a
// 3 x 4 matrix. SeisSol's poroelastic material has an isotropic solid frame, which pins down the
// structure of that matrix completely -- and thereby the row/column bookkeeping in
// admittanceFromEigendecomposition, which cannot be checked against a closed form here.
// ---------------------------------------------------------------------------
TEST_CASE("Poroelastic DR impedance and lateral stress structure" *
          doctest::test_suite("dynamicrupture")) {
  using seissol::initializer::model::checkFaultImpedance;
  using seissol::initializer::model::computeFaultImpedance;

  const auto material = testPoroMaterial();
  const auto impedance = computeFaultImpedance(material, material);

  SUBCASE("impedance invariants") {
    const auto violation = checkFaultImpedance(impedance);
    REQUIRE_MESSAGE(!violation.has_value(), violation.value_or(""));
    // homogeneous fault: the traction is split evenly
    const auto half = seissol::initializer::model::DrMatrix::Identity() * 0.5;
    CHECK((impedance.bPlus - half).cwiseAbs().maxCoeff() < 1e-9);
  }

  SUBCASE("lateral stress structure of an isotropic solid frame") {
    const auto& lateral = impedance.lateralStressPlus;
    const double scale = lateral.cwiseAbs().maxCoeff();
    REQUIRE(scale > 0.0);

    // sigma_sd is not driven by anything travelling along the fault normal
    CHECK(lateral.row(2).cwiseAbs().maxCoeff() < 1e-9 * scale);
    // sigma_ss and sigma_dd react identically
    CHECK((lateral.row(0) - lateral.row(1)).cwiseAbs().maxCoeff() < 1e-9 * scale);
    // the two shear tractions do not produce any lateral normal stress
    CHECK(lateral.col(1).cwiseAbs().maxCoeff() < 1e-9 * scale);
    CHECK(lateral.col(2).cwiseAbs().maxCoeff() < 1e-9 * scale);
    // ... while sigma_nn and the fluid pressure do
    CHECK(std::abs(lateral(0, 0)) > 1e-3 * scale);
    CHECK(std::abs(lateral(0, 3)) > 1e-3 * scale);
  }
}

} // namespace seissol::unit_test

#endif // USE_POROELASTIC

#endif // SEISSOL_TESTS_MODEL_POROELASTICIMPEDANCE_T_H_
