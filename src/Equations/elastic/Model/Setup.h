// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_EQUATIONS_ELASTIC_MODEL_SETUP_H_
#define SEISSOL_SRC_EQUATIONS_ELASTIC_MODEL_SETUP_H_

#include "Equations/elastic/Model/Datastructures.h"
#include "GeneratedCode/init.h"
#include "Kernels/Common.h"
#include "Model/Common.h"
#include "Numerical/Eigenvalues.h"
#include "Numerical/Transformation.h"

namespace seissol::model {
using Matrix99 = Eigen::Matrix<double, 9, 9>;

template <>
struct MaterialSetup<ElasticMaterial> {
  template <typename T>
  static void
      getTransposedCoefficientMatrix(const ElasticMaterial& material, unsigned dim, T& matM) {
    matM.setZero();

    const auto lambda2mu = material.lambda + 2.0 * material.mu;
    const auto rhoInv = 1.0 / material.rho;

    switch (dim) {
    case 0:
      matM(6, 0) = -lambda2mu;
      matM(6, 1) = -material.lambda;
      matM(6, 2) = -material.lambda;
      matM(7, 3) = -material.mu;
      matM(8, 5) = -material.mu;
      matM(0, 6) = -rhoInv;
      if (!testIfAcoustic(material.mu)) {
        matM(3, 7) = -rhoInv;
        matM(5, 8) = -rhoInv;
      }
      break;

    case 1:
      matM(7, 0) = -material.lambda;
      matM(7, 1) = -lambda2mu;
      matM(7, 2) = -material.lambda;
      matM(6, 3) = -material.mu;
      matM(8, 4) = -material.mu;
      matM(1, 7) = -rhoInv;
      if (!testIfAcoustic(material.mu)) {
        matM(3, 6) = -rhoInv;
        matM(4, 8) = -rhoInv;
      }
      break;

    case 2:
      matM(8, 0) = -material.lambda;
      matM(8, 1) = -material.lambda;
      matM(8, 2) = -lambda2mu;
      matM(7, 4) = -material.mu;
      matM(6, 5) = -material.mu;
      matM(2, 8) = -rhoInv;
      if (!testIfAcoustic(material.mu)) {
        matM(5, 6) = -rhoInv;
        matM(4, 7) = -rhoInv;
      }
      break;

    default:
      break;
    }
  }

  template <typename Tloc, typename Tneigh>
  static void getTransposedGodunovState(const ElasticMaterial& local,
                                        const ElasticMaterial& neighbor,
                                        FaceType faceType,
                                        Tloc& qGodLocal,
                                        Tneigh& qGodNeighbor) {
    qGodNeighbor.setZero();

    // matR == eigenvector matrix

    // Eigenvectors are precomputed
    Matrix99 matR = Matrix99::Zero();

    if (testIfAcoustic(local.mu)) {
      matR(0, 0) = local.lambda;
      matR(1, 0) = local.lambda;
      matR(2, 0) = local.lambda;
      matR(6, 0) = std::sqrt((local.lambda) / local.rho);

      // scale for better condition number of matR
      matR(3, 1) = local.lambda;
      matR(5, 2) = local.lambda;
    } else {
      matR(0, 0) = local.lambda + 2 * local.mu;
      matR(1, 0) = local.lambda;
      matR(2, 0) = local.lambda;
      matR(6, 0) = std::sqrt((local.lambda + 2 * local.mu) / local.rho);

      matR(3, 1) = local.mu;
      matR(7, 1) = std::sqrt(local.mu / local.rho);

      matR(5, 2) = local.mu;
      matR(8, 2) = std::sqrt(local.mu / local.rho);
    }

    // scale for better condition number of matR
    matR(4, 3) = local.lambda + 2 * local.mu;
    matR(1, 4) = local.lambda + 2 * local.mu;
    matR(2, 5) = local.lambda + 2 * local.mu;

    if (testIfAcoustic(neighbor.mu)) {
      // scale for better condition number of matR
      matR(7, 6) = neighbor.lambda;
      matR(8, 7) = neighbor.lambda;

      matR(0, 8) = neighbor.lambda;
      matR(1, 8) = neighbor.lambda;
      matR(2, 8) = neighbor.lambda;
      matR(6, 8) = -std::sqrt((neighbor.lambda + 2 * neighbor.mu) / neighbor.rho);
    } else {
      matR(5, 6) = neighbor.mu;
      matR(8, 6) = -std::sqrt(neighbor.mu / neighbor.rho);

      matR(3, 7) = neighbor.mu;
      matR(7, 7) = -std::sqrt(neighbor.mu / neighbor.rho);

      matR(0, 8) = neighbor.lambda + 2 * neighbor.mu;
      matR(1, 8) = neighbor.lambda;
      matR(2, 8) = neighbor.lambda;
      matR(6, 8) = -std::sqrt((neighbor.lambda + 2 * neighbor.mu) / neighbor.rho);
    }

    if (faceType == FaceType::FreeSurface) {
      const MaterialType materialtype =
          testIfAcoustic(local.mu) ? MaterialType::Acoustic : MaterialType::Elastic;
      getTransposedFreeSurfaceGodunovState(materialtype, qGodLocal, qGodNeighbor, matR);
    } else {
      // chi selects the waves travelling away from the local cell; everything else (the incoming
      // waves and the non-propagating null-space modes) belongs to the local subsystem. Since the
      // two selectors are complementary by construction, only one projector has to be computed --
      // the other one follows exactly as I - godunov. Deriving it this way is both cheaper and
      // strictly more accurate than a second solve: it makes qGodLocal + qGodNeighbor == I hold to
      // the last bit for every material, whereas two independent solves leave a residual of order
      // cond(matR) * eps (~1e-9 for realistic material contrasts).
      Matrix99 chi = Matrix99::Zero();
      chi(0, 0) = 1.0;
      if (!testIfAcoustic(local.mu)) {
        // shear waves only propagate in solids; in a fluid they degenerate into static modes and
        // stay in the local subsystem
        chi(1, 1) = 1.0;
        chi(2, 2) = 1.0;
      }

      const auto matRT = matR.transpose();
      // Deliberately partialPivLu and not a rank-revealing factorization: matR is a regular
      // eigenvector basis, so a solver that declares it rank deficient and zeroes part of the
      // solution can only ever be wrong -- and silently so. Partial pivoting also preserves the
      // sparsity of matR (measured growth factor 1.0 across all equation sets and materials),
      // which makes it markedly more accurate here than Householder QR.
      const auto matRlu = matRT.partialPivLu();
      const auto godunov = matRlu.solve(chi * matRT).eval();

      // qGodLocal = I - qGodNeighbor
      for (unsigned i = 0; i < godunov.cols(); ++i) {
        for (unsigned j = 0; j < godunov.rows(); ++j) {
          const double identity = (i == j) ? 1.0 : 0.0;
          qGodLocal(i, j) = identity - godunov(i, j);
          qGodNeighbor(i, j) = godunov(i, j);
        }
      }
    }
  }


  static ElasticMaterial
      getRotatedMaterialCoefficients(const std::array<double, 36>& /*rotationParameters*/,
                                     ElasticMaterial& material) {
    return material;
  }

  static void initializeSpecificLocalData(const ElasticMaterial& material,
                                          double timeStepWidth,
                                          typename ElasticMaterial::Solver::LocalData* localData) {}

  static void
      initializeSpecificNeighborData(const ElasticMaterial& material,
                                     typename ElasticMaterial::Solver::NeighborData* neighborData) {
  }
  static void getPlaneWaveOperator(
      const ElasticMaterial& material,
      const double n[3],
      std::complex<double> mdata[ElasticMaterial::NumQuantities * ElasticMaterial::NumQuantities]) {
    getElasticPlaneWaveOperator(material, n, mdata);
  }

  template <typename T>
  static void getTransposedSourceCoefficientTensor(const ElasticMaterial& material,
                                                   T& sourceMatrix) {}
};

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_ELASTIC_MODEL_SETUP_H_
