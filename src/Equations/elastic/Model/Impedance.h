// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EQUATIONS_ELASTIC_MODEL_IMPEDANCE_H_
#define SEISSOL_SRC_EQUATIONS_ELASTIC_MODEL_IMPEDANCE_H_

#include "Equations/ImpedanceBase.h"
#include "Equations/elastic/Model/Datastructures.h"

#include <Eigen/Dense>
#include <array>
#include <cstddef>

namespace seissol::model {

template <>
struct ImpedanceCompute<ElasticMaterial> {
  static constexpr std::array<std::size_t, 3> TractionIndices{0, 3, 5};
  static constexpr std::array<std::size_t, 3> VelocityIndices{6, 7, 8};
  static constexpr std::size_t Dim = TractionIndices.size();
  using Matrix = AdmittanceMatrix<Dim>;
  using LateralMatrix = LateralStressMatrix<Dim>;

  static Matrix signature() { return Matrix::Identity(); }

  /**
   * Isotropic closed form. The Christoffel matrix is diag(lambda + 2 mu, mu, mu) for every
   * direction, hence
   *
   *   Y = diag(1 / Zp, 1 / Zs, 1 / Zs),   Zp = rho * cp,   Zs = rho * cs.
   *
   * A plane wave along the normal changes sigma_ss and sigma_dd by
   * lambda / (lambda + 2 mu) = 1 - 2 (cs/cp)^2 times sigma_nn, and leaves sigma_sd alone.
   *
   * The isotropic branch of initializeDynamicRuptureMatrices and the fault receiver output do not
   * need the matrices: they work with the scalar impedances Zp and Zs (computed the same way as
   * here) and with the isotropic lateral stress formula directly.
   *
   * Not defined for a fluid (mu = 0): the shear admittance is infinite there.
   */
  static Matrix admittance(const ElasticMaterial& materialLocal,
                           LateralMatrix* lateralStress = nullptr) {
    const double zp = materialLocal.getDensity() * materialLocal.getPWaveSpeed();
    const double zs = materialLocal.getDensity() * materialLocal.getSWaveSpeed();

    if (lateralStress != nullptr) {
      const double ratio = materialLocal.lambda / (materialLocal.lambda + 2 * materialLocal.mu);
      *lateralStress = LateralMatrix::Zero();
      (*lateralStress)(0, 0) = ratio;
      (*lateralStress)(1, 0) = ratio;
    }

    return Eigen::Vector3d(1 / zp, 1 / zs, 1 / zs).asDiagonal();
  }
};

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_ELASTIC_MODEL_IMPEDANCE_H_
