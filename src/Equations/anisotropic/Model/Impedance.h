// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EQUATIONS_ANISOTROPIC_MODEL_IMPEDANCE_H_
#define SEISSOL_SRC_EQUATIONS_ANISOTROPIC_MODEL_IMPEDANCE_H_

#include "Equations/ImpedanceBase.h"
#include "Equations/anisotropic/Model/Datastructures.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <array>
#include <cstddef>

namespace seissol::model {

template <>
struct ImpedanceCompute<AnisotropicMaterial> {
  static constexpr std::array<std::size_t, 3> TractionIndices{0, 3, 5};
  static constexpr std::array<std::size_t, 3> VelocityIndices{6, 7, 8};
  static constexpr std::size_t Dim = TractionIndices.size();
  using Matrix = AdmittanceMatrix<Dim>;
  using LateralMatrix = LateralStressMatrix<Dim>;

  static Matrix signature() { return Matrix::Identity(); }

  /**
   * Christoffel matrix for propagation along the local x axis, i.e. along the
   * fault normal once the material has been rotated with the Bond matrix:
   *
   *   Gamma_ik = C_{i1k1}
   *
   * In Voigt notation this is the (11, 16, 15 / 16, 66, 56 / 15, 56, 55) block --
   * the same coefficients that appear in rows 6..8 of
   * getTransposedCoefficientMatrix(material, 0, .).
   */
  static Eigen::Matrix3d christoffelMatrix(const AnisotropicMaterial& materialLocal) {
    Eigen::Matrix3d gamma;
    gamma(0, 0) = materialLocal.c11;
    gamma(0, 1) = materialLocal.c16;
    gamma(0, 2) = materialLocal.c15;
    gamma(1, 1) = materialLocal.c66;
    gamma(1, 2) = materialLocal.c56;
    gamma(2, 2) = materialLocal.c55;
    gamma(1, 0) = gamma(0, 1);
    gamma(2, 0) = gamma(0, 2);
    gamma(2, 1) = gamma(1, 2);
    return gamma;
  }

  /**
   * Reconstructs the stress components which do not take part in the fault-normal Riemann problem
   * from the traction difference across the fault.
   *
   * A plane wave travelling along the local x axis only has the strain rates
   * (eps_1, eps_6, eps_5) = (du/dx, dv/dx, dw/dx), so in Voigt notation
   *
   *   [d sigma_1; d sigma_6; d sigma_5] = Gamma * [d eps_1; d eps_6; d eps_5]
   *   [d sigma_2; d sigma_3; d sigma_4] = C[{2,3,4},{1,6,5}] * [d eps_1; d eps_6; d eps_5]
   *
   * with the very same Gamma as above. Eliminating the strain rates gives
   *
   *   [d sigma_ss; d sigma_dd; d sigma_sd] = C[{2,3,4},{1,6,5}] * Gamma^-1 * [traction difference].
   *
   * For an isotropic material this reduces to d sigma_ss = d sigma_dd = d sigma_nn * (1 - 2
   * (cs/cp)^2) and d sigma_sd = 0.
   */
  static Eigen::Matrix3d lateralStressFromChristoffel(const AnisotropicMaterial& materialLocal) {
    Eigen::Matrix3d rest;
    // rows: Voigt 2, 3, 4 (sigma_ss, sigma_dd, sigma_sd); columns: Voigt 1, 6, 5
    rest(0, 0) = materialLocal.c12;
    rest(0, 1) = materialLocal.c26;
    rest(0, 2) = materialLocal.c25;
    rest(1, 0) = materialLocal.c13;
    rest(1, 1) = materialLocal.c36;
    rest(1, 2) = materialLocal.c35;
    rest(2, 0) = materialLocal.c14;
    rest(2, 1) = materialLocal.c46;
    rest(2, 2) = materialLocal.c45;
    return rest * christoffelMatrix(materialLocal).inverse();
  }

  /**
   * Admittance of an elastic/anisotropic half space in closed form.
   *
   *   Z = sum_alpha rho * v_alpha * p_alpha (x) p_alpha = (rho * Gamma)^(1/2)
   *   Y = Z^-1
   *
   * because Gamma p_alpha = rho v_alpha^2 p_alpha with an orthonormal p_alpha.
   * Gamma is symmetric positive definite, so a self-adjoint eigensolver is used
   * and the result is symmetric positive definite by construction -- also when
   * qS1 and qS2 are degenerate, which is exactly where the general 9x9 complex
   * eigendecomposition loses accuracy.
   */
  static Eigen::Matrix3d admittanceFromChristoffel(const Eigen::Matrix3d& gamma, double rho) {
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(rho * gamma);
    const Eigen::Vector3d impedanceEigenvalues = solver.eigenvalues().cwiseSqrt();
    return solver.eigenvectors() * impedanceEigenvalues.cwiseInverse().asDiagonal() *
           solver.eigenvectors().transpose();
  }

  /**
   * Christoffel matrix of the fault normal direction, recovered from the admittance stored for one
   * side of the face.
   *
   * `admittanceFromChristoffel` yields Y = (rho * Gamma)^-1/2, hence Gamma = (Y * Y)^-1 / rho. The
   * point of going back is the shear block: for a unit slip direction d in the fault plane,
   *
   *   d^T Gamma d = C_ijkl n_i d_j n_k d_l
   *
   * is the stiffness a shear dislocation works against, i.e. the modulus turning potency into
   * seismic moment. It equals mu for an isotropic material and any fault orientation.
   */
  static Eigen::Matrix3d christoffelFromAdmittance(const Eigen::Matrix3d& admittance, double rho) {
    return (admittance * admittance).inverse() / rho;
  }

  static Matrix admittance(const AnisotropicMaterial& materialLocal,
                           LateralMatrix* lateralStress = nullptr) {
    if (lateralStress != nullptr) {
      *lateralStress = lateralStressFromChristoffel(materialLocal);
    }
    return admittanceFromChristoffel(christoffelMatrix(materialLocal), materialLocal.rho);
  }
};

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_ANISOTROPIC_MODEL_IMPEDANCE_H_
