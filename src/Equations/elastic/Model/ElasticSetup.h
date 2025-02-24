// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_EQUATIONS_ELASTIC_MODEL_ELASTICSETUP_H_
#define SEISSOL_SRC_EQUATIONS_ELASTIC_MODEL_ELASTICSETUP_H_

#include "Kernels/Common.h"
#include "Model/Common.h"
#include "Numerical/Eigenvalues.h"
#include "Numerical/Transformation.h"
#include "generated_code/init.h"

namespace seissol::model {
using Matrix99 = Eigen::Matrix<double, 9, 9>;

template <typename T>
inline void
    getTransposedCoefficientMatrix(const ElasticMaterial& i_material, unsigned i_dim, T& o_M) {
  o_M.setZero();

  real lambda2mu = i_material.lambda + 2.0 * i_material.mu;
  real rhoInv = 1.0 / i_material.rho;

  switch (i_dim) {
  case 0:
    o_M(6, 0) = -lambda2mu;
    o_M(6, 1) = -i_material.lambda;
    o_M(6, 2) = -i_material.lambda;
    o_M(7, 3) = -i_material.mu;
    o_M(8, 5) = -i_material.mu;
    o_M(0, 6) = -rhoInv;
    if (!testIfAcoustic(i_material.mu)) {
      o_M(3, 7) = -rhoInv;
      o_M(5, 8) = -rhoInv;
    }
    break;

  case 1:
    o_M(7, 0) = -i_material.lambda;
    o_M(7, 1) = -lambda2mu;
    o_M(7, 2) = -i_material.lambda;
    o_M(6, 3) = -i_material.mu;
    o_M(8, 4) = -i_material.mu;
    o_M(1, 7) = -rhoInv;
    if (!testIfAcoustic(i_material.mu)) {
      o_M(3, 6) = -rhoInv;
      o_M(4, 8) = -rhoInv;
    }
    break;

  case 2:
    o_M(8, 0) = -i_material.lambda;
    o_M(8, 1) = -i_material.lambda;
    o_M(8, 2) = -lambda2mu;
    o_M(7, 4) = -i_material.mu;
    o_M(6, 5) = -i_material.mu;
    o_M(2, 8) = -rhoInv;
    if (!testIfAcoustic(i_material.mu)) {
      o_M(5, 6) = -rhoInv;
      o_M(4, 7) = -rhoInv;
    }
    break;

  default:
    break;
  }
}

template <typename Tloc, typename Tneigh>
inline void getTransposedGodunovState(const ElasticMaterial& local,
                                      const ElasticMaterial& neighbor,
                                      FaceType faceType,
                                      Tloc& QgodLocal,
                                      Tneigh& QgodNeighbor) {
  QgodNeighbor.setZero();

  // Eigenvectors are precomputed
  Matrix99 R = Matrix99::Zero();

  if (testIfAcoustic(local.mu)) {
    R(0, 0) = local.lambda;
    R(1, 0) = local.lambda;
    R(2, 0) = local.lambda;
    R(6, 0) = std::sqrt((local.lambda) / local.rho);

    // scale for better condition number of R
    R(3, 1) = local.lambda;
    R(5, 2) = local.lambda;
  } else {
    R(0, 0) = local.lambda + 2 * local.mu;
    R(1, 0) = local.lambda;
    R(2, 0) = local.lambda;
    R(6, 0) = std::sqrt((local.lambda + 2 * local.mu) / local.rho);

    R(3, 1) = local.mu;
    R(7, 1) = std::sqrt(local.mu / local.rho);

    R(5, 2) = local.mu;
    R(8, 2) = std::sqrt(local.mu / local.rho);
  }

  // scale for better condition number of R
  R(4, 3) = local.lambda + 2 * local.mu;
  R(1, 4) = local.lambda + 2 * local.mu;
  R(2, 5) = local.lambda + 2 * local.mu;

  if (testIfAcoustic(neighbor.mu)) {
    // scale for better condition number of R
    R(7, 6) = neighbor.lambda;
    R(8, 7) = neighbor.lambda;

    R(0, 8) = neighbor.lambda;
    R(1, 8) = neighbor.lambda;
    R(2, 8) = neighbor.lambda;
    R(6, 8) = -std::sqrt((neighbor.lambda + 2 * neighbor.mu) / neighbor.rho);
  } else {
    R(5, 6) = neighbor.mu;
    R(8, 6) = -std::sqrt(neighbor.mu / neighbor.rho);

    R(3, 7) = neighbor.mu;
    R(7, 7) = -std::sqrt(neighbor.mu / neighbor.rho);

    R(0, 8) = neighbor.lambda + 2 * neighbor.mu;
    R(1, 8) = neighbor.lambda;
    R(2, 8) = neighbor.lambda;
    R(6, 8) = -std::sqrt((neighbor.lambda + 2 * neighbor.mu) / neighbor.rho);
  }

  if (faceType == FaceType::FreeSurface) {
    MaterialType materialtype =
        testIfAcoustic(local.mu) ? MaterialType::Acoustic : MaterialType::Elastic;
    getTransposedFreeSurfaceGodunovState(materialtype, QgodLocal, QgodNeighbor, R);
  } else {
    Matrix99 chi = Matrix99::Zero();
    if (!testIfAcoustic(local.mu)) {
      chi(2, 2) = 1.0;
      chi(1, 1) = 1.0;
    }
    chi(0, 0) = 1.0;

    const auto godunov = ((R * chi) * R.inverse()).eval();

    // QgodLocal = I - QgodNeighbor
    for (unsigned i = 0; i < godunov.cols(); ++i) {
      for (unsigned j = 0; j < godunov.rows(); ++j) {
        QgodLocal(i, j) = -godunov(j, i);
        QgodNeighbor(i, j) = godunov(j, i);
      }
    }
    for (unsigned idx = 0; idx < 9; ++idx) {
      QgodLocal(idx, idx) += 1.0;
    }
  }
}
} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_ELASTIC_MODEL_ELASTICSETUP_H_
