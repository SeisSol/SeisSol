// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_MODEL_POROELASTICIMPEDANCE_T_H_
#define SEISSOL_TESTS_MODEL_POROELASTICIMPEDANCE_T_H_

// The closed form only depends on the material parameters and runs in every build. Comparing it
// with the eigendecomposition needs MaterialSetup<PoroElasticMaterial>, and the traction matrix
// pattern is generated code -- those two only run in a poroelastic build.

#include <doctest.h>

#include "Alignment.h"
#include "Equations/Datastructures.h"
#include "Equations/Impedance.h"
#include "Equations/ImpedanceBase.h"
#include "Equations/Setup.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "ImpedanceReference.h"
#include "Initializer/Model/DynamicRuptureImpedance.h"
#include "Kernels/Precision.h"

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace seissol::unit_test {

using PoroelasticImpedance = seissol::model::ImpedanceCompute<seissol::model::PoroElasticMaterial>;

/// bulkSolid, rho, lambda, mu, porosity, permeability, tortuosity, bulkFluid, rhoFluid, viscosity
inline model::PoroElasticMaterial testPoroMaterial(double porosity, double tortuosity) {
  return model::PoroElasticMaterial(std::vector<double>{
      3.60e10, 2650.0, 4.0e9, 6.0e9, porosity, 1.0e-13, tortuosity, 2.2e9, 1000.0, 1.0e-3});
}

#ifdef USE_POROELASTIC
// ---------------------------------------------------------------------------
// The poroelastic interface has four variables, (sigma_nn, sigma_ns, sigma_nd, p) against
// (v_n, v_s, v_d, q_n). SeisSol's poroelastic frame is isotropic, so the generalized wave
// impedance Z = Mass # Gamma has a closed form -- a 2x2 fast/slow P block plus two shear scalars.
// Checking it against the general eigendecomposition validates both routes against each other.
// ---------------------------------------------------------------------------
TEST_CASE("Poroelastic DR impedance closed form agrees with the eigendecomposition" *
          doctest::test_suite("dynamicrupture")) {
  using LateralMatrix = PoroelasticImpedance::LateralMatrix;

  // the eigendecomposition of the 13x13 Jacobian is comparatively ill conditioned here
  // (cond(R_t) reaches ~1e6 at high porosity), so this is the accuracy we can demand of it
  constexpr double Epsilon = 1e-7;

  for (const double porosity : {0.05, 0.2, 0.4}) {
    for (const double tortuosity : {1.0, 1.5, 3.0}) {
      const auto sweepMaterial = testPoroMaterial(porosity, tortuosity);

      LateralMatrix lateralClosed = LateralMatrix::Zero();
      const auto admittanceClosed =
          seissol::model::computeAdmittance(sweepMaterial, &lateralClosed);

      LateralMatrix lateralEigen = LateralMatrix::Zero();
      const auto admittanceEigen = admittanceFromEigendecomposition(sweepMaterial, &lateralEigen);

      CHECK((admittanceClosed - admittanceEigen).cwiseAbs().maxCoeff() <
            Epsilon * admittanceEigen.cwiseAbs().maxCoeff());
      CHECK((lateralClosed - lateralEigen).cwiseAbs().maxCoeff() <
            Epsilon * lateralEigen.cwiseAbs().maxCoeff());
    }
  }
}
#endif // USE_POROELASTIC

TEST_CASE("Poroelastic DR impedance closed form" * doctest::test_suite("dynamicrupture")) {
  using seissol::initializer::model::checkFaultImpedance;
  using seissol::initializer::model::computeFaultImpedance;
  using DrMatrix = PoroelasticImpedance::Matrix;

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

  SUBCASE("scalar shear impedances of a bimaterial fault") {
    // the values initializeDynamicRuptureMatrices copies into ImpedancesAndEta, which the friction
    // update, the slip accumulation and the receiver output read
    const auto plus = testPoroMaterial(0.1, 1.0);
    const auto minus = testPoroMaterial(0.3, 2.0);
    const auto bimaterial = computeFaultImpedance(plus, minus);

    const auto shearImpedance = [](const model::PoroElasticMaterial& material) {
      const double m = material.rhoFluid * material.tortuosity / material.porosity;
      const double rhoBar =
          (1 - material.porosity) * material.rho + material.porosity * material.rhoFluid;
      return std::sqrt(material.mu * (rhoBar - material.rhoFluid * material.rhoFluid / m));
    };
    const double zs = shearImpedance(plus);
    const double zsNeig = shearImpedance(minus);

    CHECK(bimaterial.admittancePlus(1, 1) == doctest::Approx(1.0 / zs).epsilon(1e-12));
    CHECK(bimaterial.admittanceMinus(1, 1) == doctest::Approx(1.0 / zsNeig).epsilon(1e-12));
    CHECK(bimaterial.eta(1, 1) == doctest::Approx(zs * zsNeig / (zs + zsNeig)).epsilon(1e-12));

    // a scalar is exact only because the shear rows carry the same entry twice and do not couple
    // to the fault normal or to the fluid
    CHECK(bimaterial.eta(2, 2) == doctest::Approx(bimaterial.eta(1, 1)).epsilon(1e-12));
    CHECK(std::abs(bimaterial.eta(1, 0)) < 1e-10 * bimaterial.eta(1, 1));
    CHECK(std::abs(bimaterial.eta(1, 3)) < 1e-10 * bimaterial.eta(1, 1));
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

#ifdef USE_POROELASTIC
// ---------------------------------------------------------------------------
// The traction averaging matrices map all four interface tractions to the three components the
// frictional work is computed with, so their sparsity pattern has to hold the fluid pressure row.
// Writing a four row matrix into a pattern that only has three leaves the pattern and silently
// overwrites a neighboring entry, which is what this pins down.
// ---------------------------------------------------------------------------
TEST_CASE("Poroelastic traction matrix pattern" * doctest::test_suite("dynamicrupture")) {
  constexpr std::array<std::size_t, 4> StoredRows{0, 3, 5, 9};
  constexpr std::size_t Rows = StoredRows.size();
  constexpr std::size_t Columns = 3;

  REQUIRE(tensor::tractionPlusMatrix::size() == Rows * Columns);
  REQUIRE(tensor::tractionMinusMatrix::size() == Rows * Columns);

  alignas(Alignment) real data[tensor::tractionPlusMatrix::size()]{};
  auto view = init::tractionPlusMatrix::view::create(data);
  view.setZero();
  for (std::size_t col = 0; col < Columns; ++col) {
    for (std::size_t row = 0; row < Rows; ++row) {
      view(StoredRows[row], col) = static_cast<real>(10 * col + row);
    }
  }

  // column major within the stored rows, the layout the friction energy indexing relies on
  for (std::size_t col = 0; col < Columns; ++col) {
    for (std::size_t row = 0; row < Rows; ++row) {
      CHECK(data[Rows * col + row] == doctest::Approx(10.0 * col + row));
    }
  }
}

#endif // USE_POROELASTIC

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_MODEL_POROELASTICIMPEDANCE_T_H_
