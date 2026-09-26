// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EQUATIONS_POROELASTIC_MODEL_IMPEDANCE_H_
#define SEISSOL_SRC_EQUATIONS_POROELASTIC_MODEL_IMPEDANCE_H_

#include "Equations/ImpedanceBase.h"
#include "Equations/poroelastic/Model/Datastructures.h"
#include "Equations/poroelastic/Model/Helper.h"

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <cstddef>

namespace seissol::model {

template <>
struct ImpedanceCompute<PoroElasticMaterial> {
  /// the fourth pair is the fluid pressure against the normal filtration velocity
  static constexpr std::array<std::size_t, 4> TractionIndices{0, 3, 5, 9};
  static constexpr std::array<std::size_t, 4> VelocityIndices{6, 7, 8, 10};
  static constexpr std::size_t Dim = TractionIndices.size();
  using Matrix = AdmittanceMatrix<Dim>;
  using LateralMatrix = LateralStressMatrix<Dim>;

  /**
   * Sign convention of the stored interface traction relative to the energy conjugate one.
   *
   * The energy conjugate partner of the normal filtration velocity is *minus* the pore pressure,
   * while SeisSol carries +p as the fourth traction component. The stored admittance is therefore
   * Y_stored = Y_physical * signature(), and only Y_stored * signature() is self-adjoint.
   */
  static Matrix signature() {
    Matrix signature = Matrix::Identity();
    signature(Dim - 1, Dim - 1) = -1.0;
    return signature;
  }

  /// Closed form square root of a 2x2 matrix with positive eigenvalues (Cayley-Hamilton).
  static Eigen::Matrix2d matrixSqrt2x2(const Eigen::Matrix2d& matrix) {
    const double determinant = std::sqrt(matrix.determinant());
    return (matrix + determinant * Eigen::Matrix2d::Identity()) /
           std::sqrt(matrix.trace() + 2 * determinant);
  }

  /**
   * Generalized wave impedance of a poroelastic half space, in closed form.
   *
   * The interface variables are T = (sigma_nn, sigma_ns, sigma_nd, -p) against
   * U = (v_n, v_s, v_d, q_n), and the Biot system is symmetric hyperbolic in them:
   *
   *   Mass * dU/dt = dT/dx ,   dT/dt = Gamma * dU/dx
   *
   * with a symmetric positive definite mass matrix and a symmetric Gamma. Both are read straight
   * off MaterialSetup<PoroElasticMaterial>::getTransposedCoefficientMatrix(., 0, .); note that the
   * static condensation of the two tangential filtration velocities (which carry no flux, hence
   * zero wave speed) is already baked into rho1 = rhoBar - rhoFluid^2 / m.
   *
   * The impedance is then the matrix geometric mean, characterised by Z * Mass^-1 * Z = Gamma:
   *
   *   Z = Mass # Gamma = Mass^1/2 (Mass^-1/2 Gamma Mass^-1/2)^1/2 Mass^1/2
   *
   * which for Mass = rho * I collapses to the elastic (rho * Gamma)^1/2. Because SeisSol's
   * poroelastic frame is isotropic, Mass and Gamma are block diagonal -- a 2x2 fast/slow P block on
   * (v_n, q_n) plus two identical shear scalars -- so no eigensolver is needed at all.
   *
   * Wave speeds recovered from this: c_S = sqrt(mu / rho1), the two Biot P speeds from the 2x2
   * block.
   */
  static Matrix admittance(const PoroElasticMaterial& materialLocal,
                           LateralMatrix* lateralStress = nullptr) {
    const auto params = getAdditionalParameters(materialLocal);

    // Gamma, rows/columns ordered as (n, s, d, fluid)
    Matrix gamma;
    gamma << params.cBar(0, 0), params.cBar(0, 5), params.cBar(0, 4), params.M * params.alpha(0),
        params.cBar(5, 0), params.cBar(5, 5), params.cBar(5, 4), params.M * params.alpha(5),
        params.cBar(4, 0), params.cBar(4, 5), params.cBar(4, 4), params.M * params.alpha(4),
        params.M * params.alpha(0), params.M * params.alpha(5), params.M * params.alpha(4),
        params.M;

    // fast/slow P block on (v_n, q_n); the shear rows decouple for an isotropic frame
    Eigen::Matrix2d massP;
    massP << params.rhoBar, materialLocal.rhoFluid, materialLocal.rhoFluid, params.m;
    Eigen::Matrix2d gammaP;
    gammaP << gamma(0, 0), gamma(0, 3), gamma(3, 0), gamma(3, 3);
    const Eigen::Matrix2d impedanceP = massP * matrixSqrt2x2(massP.inverse() * gammaP);

    Matrix impedance = Matrix::Zero();
    impedance(0, 0) = impedanceP(0, 0);
    impedance(0, 3) = impedanceP(0, 1);
    impedance(3, 0) = impedanceP(1, 0);
    impedance(3, 3) = impedanceP(1, 1);
    impedance(1, 1) = std::sqrt(params.rho1 * gamma(1, 1));
    impedance(2, 2) = std::sqrt(params.rho1 * gamma(2, 2));

    if (lateralStress != nullptr) {
      // same construction as the anisotropic lateralStressFromChristoffel, with the pore pressure
      // as a fourth column
      LateralMatrix rest;
      rest << params.cBar(1, 0), params.cBar(1, 5), params.cBar(1, 4), params.M * params.alpha(1),
          params.cBar(2, 0), params.cBar(2, 5), params.cBar(2, 4), params.M * params.alpha(2),
          params.cBar(3, 0), params.cBar(3, 5), params.cBar(3, 4), params.M * params.alpha(3);
      *lateralStress = rest * gamma.inverse() * signature();
    }

    return impedance.inverse() * signature();
  }
};

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_POROELASTIC_MODEL_IMPEDANCE_H_
