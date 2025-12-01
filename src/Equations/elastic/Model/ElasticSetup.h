// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_EQUATIONS_ELASTIC_MODEL_ELASTICSETUP_H_
#define SEISSOL_SRC_EQUATIONS_ELASTIC_MODEL_ELASTICSETUP_H_

#include "Equations/elastic/Model/Datastructures.h"
#include "Equations/elastic/Model/IntegrationData.h"
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

    // Eigenvectors are precomputed
    Matrix99 mR = Matrix99::Zero();

    if (testIfAcoustic(local.mu)) {
      mR(0, 0) = local.lambda;
      mR(1, 0) = local.lambda;
      mR(2, 0) = local.lambda;
      mR(6, 0) = std::sqrt((local.lambda) / local.rho);

      // scale for better condition number of mR
      mR(3, 1) = local.lambda;
      mR(5, 2) = local.lambda;
    } else {
      mR(0, 0) = local.lambda + 2 * local.mu;
      mR(1, 0) = local.lambda;
      mR(2, 0) = local.lambda;
      mR(6, 0) = std::sqrt((local.lambda + 2 * local.mu) / local.rho);

      mR(3, 1) = local.mu;
      mR(7, 1) = std::sqrt(local.mu / local.rho);

      mR(5, 2) = local.mu;
      mR(8, 2) = std::sqrt(local.mu / local.rho);
    }

    // scale for better condition number of mR
    mR(4, 3) = local.lambda + 2 * local.mu;
    mR(1, 4) = local.lambda + 2 * local.mu;
    mR(2, 5) = local.lambda + 2 * local.mu;

    if (testIfAcoustic(neighbor.mu)) {
      // scale for better condition number of mR
      mR(7, 6) = neighbor.lambda;
      mR(8, 7) = neighbor.lambda;

      mR(0, 8) = neighbor.lambda;
      mR(1, 8) = neighbor.lambda;
      mR(2, 8) = neighbor.lambda;
      mR(6, 8) = -std::sqrt((neighbor.lambda + 2 * neighbor.mu) / neighbor.rho);
    } else {
      mR(5, 6) = neighbor.mu;
      mR(8, 6) = -std::sqrt(neighbor.mu / neighbor.rho);

      mR(3, 7) = neighbor.mu;
      mR(7, 7) = -std::sqrt(neighbor.mu / neighbor.rho);

      mR(0, 8) = neighbor.lambda + 2 * neighbor.mu;
      mR(1, 8) = neighbor.lambda;
      mR(2, 8) = neighbor.lambda;
      mR(6, 8) = -std::sqrt((neighbor.lambda + 2 * neighbor.mu) / neighbor.rho);
    }

    if (faceType == FaceType::FreeSurface) {
      const MaterialType materialtype =
          testIfAcoustic(local.mu) ? MaterialType::Acoustic : MaterialType::Elastic;
      getTransposedFreeSurfaceGodunovState(materialtype, qGodLocal, qGodNeighbor, mR);
    } else {
      Matrix99 chi = Matrix99::Zero();
      if (!testIfAcoustic(local.mu)) {
        chi(2, 2) = 1.0;
        chi(1, 1) = 1.0;
      }
      chi(0, 0) = 1.0;

      const auto godunov = ((mR * chi) * mR.inverse()).eval();

      // qGodLocal = I - qGodNeighbor
      for (unsigned i = 0; i < godunov.cols(); ++i) {
        for (unsigned j = 0; j < godunov.rows(); ++j) {
          qGodLocal(i, j) = -godunov(j, i);
          qGodNeighbor(i, j) = godunov(j, i);
        }
      }
      for (unsigned idx = 0; idx < 9; ++idx) {
        qGodLocal(idx, idx) += 1.0;
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

  static ElasticMaterial
      getRotatedMaterialCoefficients(const std::array<double, 36>& /*rotationParameters*/,
                                     ElasticMaterial& material) {
    return material;
  }

  static void initializeSpecificLocalData(const ElasticMaterial& material,
                                          double timeStepWidth,
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
