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
#include <Model/Datastructures.h>
#include <Model/IntegrationData.h>

namespace seissol::model {
using Matrix99 = Eigen::Matrix<double, 9, 9>;

template <>
struct MaterialSetup<ElasticMaterial> {
  template <typename T>
  static void
      getTransposedCoefficientMatrix(const ElasticMaterial& material, unsigned dim, T& matM) {
    matM.setZero();

    real lambda2mu = material.lambda + 2.0 * material.mu;
    real rhoInv = 1.0 / material.rho;

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

  static void getFaceRotationMatrix(const VrtxCoords normal,
                                    const VrtxCoords tangent1,
                                    const VrtxCoords tangent2,
                                    init::T::view::type& matT,
                                    init::Tinv::view::type& matTinv) {
    matT.setZero();
    matTinv.setZero();

    seissol::transformations::symmetricTensor2RotationMatrix(
        normal, tangent1, tangent2, matT, 0, 0);
    seissol::transformations::tensor1RotationMatrix(normal, tangent1, tangent2, matT, 6, 6);

    seissol::transformations::inverseSymmetricTensor2RotationMatrix(
        normal, tangent1, tangent2, matTinv, 0, 0);
    seissol::transformations::inverseTensor1RotationMatrix(
        normal, tangent1, tangent2, matTinv, 6, 6);
  }

  static ElasticMaterial getRotatedMaterialCoefficients(real rotationParameters[36],
                                                        ElasticMaterial& material) {
    return material;
  }

  static void initializeSpecificLocalData(const ElasticMaterial& material,
                                          real timeStepWidth,
                                          ElasticLocalData* localData) {}

  static void initializeSpecificNeighborData(const ElasticMaterial& material,
                                             ElasticNeighborData* localData) {}
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

#endif // SEISSOL_SRC_EQUATIONS_ELASTIC_MODEL_ELASTICSETUP_H_
