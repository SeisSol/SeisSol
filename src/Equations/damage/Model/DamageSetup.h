// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf
// SPDX-FileContributor: Zihua Niu

#ifndef SEISSOL_SRC_EQUATIONS_DAMAGE_MODEL_DAMAGESETUP_H_
#define SEISSOL_SRC_EQUATIONS_DAMAGE_MODEL_DAMAGESETUP_H_

#include "Kernels/Common.h"
#include "Model/Common.h"
#include "Numerical/Eigenvalues.h"
#include "Numerical/Transformation.h"
#include "generated_code/init.h"
#include <Equations/damage/Model/Datastructures.h>
#include <Equations/damage/Model/IntegrationData.h>
#include <Geometry/MeshDefinition.h>

namespace seissol::model {
using Matrix1010 = Eigen::Matrix<double, 10, 10>;

template <>
struct MaterialSetup<DamageMaterial> {
  template <typename T>
  static void
      getTransposedCoefficientMatrix(const DamageMaterial& material, unsigned dim, T& matM) {
    matM.setZero();

    real lambda2muInvRho = (material.lambda + 2.0 * material.mu) / material.rho;
    real lambdaInvRho = (material.lambda) / material.rho;
    real muInvRho = (material.mu) / material.rho;

    /* For strain-vel formula:
    {eps_xx, eps_yy, eps_zz, eps_xy, eps_yz, eps_zz, vx, vy, vz, alpha} */
    switch (dim) {
    case 0:
      matM(6, 0) = -1.0;
      matM(7, 3) = -0.5;
      matM(8, 5) = -0.5;
      matM(0, 6) = -lambda2muInvRho;
      matM(1, 6) = -lambdaInvRho;
      matM(2, 6) = -lambdaInvRho;
      if (!testIfAcoustic(material.mu)) {
        matM(3, 7) = -2.0 * muInvRho;
        matM(5, 8) = -2.0 * muInvRho;
      }
      matM(9, 6) = 0.0;
      matM(9, 7) = 0.0;
      matM(9, 8) = 0.0;
      break;

    case 1:
      matM(7, 1) = -1.0;
      matM(6, 3) = -0.5;
      matM(8, 4) = -0.5;
      matM(1, 7) = -lambda2muInvRho;
      matM(0, 7) = -lambdaInvRho;
      matM(2, 7) = -lambdaInvRho;
      if (!testIfAcoustic(material.mu)) {
        matM(3, 6) = -2.0 * muInvRho;
        matM(4, 8) = -2.0 * muInvRho;
      }
      matM(9, 6) = 0.0;
      matM(9, 7) = 0.0;
      matM(9, 8) = 0.0;
      break;

    case 2:
      matM(8, 2) = -1.0;
      matM(7, 4) = -0.5;
      matM(6, 5) = -0.5;
      matM(2, 8) = -lambda2muInvRho;
      matM(0, 8) = -lambdaInvRho;
      matM(1, 8) = -lambdaInvRho;
      if (!testIfAcoustic(material.mu)) {
        matM(5, 6) = -2.0 * muInvRho;
        matM(4, 7) = -2.0 * muInvRho;
      }
      matM(9, 6) = 0.0;
      matM(9, 7) = 0.0;
      matM(9, 8) = 0.0;
      break;

    default:
      break;
    }
  }

  template <typename Tloc, typename Tneigh>
  static void getTransposedGodunovState(const DamageMaterial& local,
                                        const DamageMaterial& neighbor,
                                        FaceType faceType,
                                        Tloc& QgodLocal,
                                        Tneigh& QgodNeighbor) {
    QgodNeighbor.setZero();

    // Eigenvectors are precomputed
    Matrix1010 R = Matrix1010::Zero();

    /* For strain-vel formula:
    {-cp, -cs, -cs, 0, 0, 0, cs, cs, cp, 0} */
    if (testIfAcoustic(local.mu)) {
      R(0, 0) = local.lambda;
      R(1, 0) = local.lambda;
      R(2, 0) = local.lambda;
      R(6, 0) = std::sqrt((local.lambda) / local.rho);

      // scale for better condition number of R
      R(3, 1) = local.lambda;
      R(5, 2) = local.lambda;
    } else {
      R(0, 0) = 1.0;
      R(6, 0) = std::sqrt((local.lambda + 2 * local.mu) / local.rho);

      R(5, 1) = 0.5;
      R(8, 1) = std::sqrt(local.mu / local.rho);

      R(3, 2) = 0.5;
      R(7, 2) = std::sqrt(local.mu / local.rho);
    }

    // scale for better condition number of R
    R(4, 3) = 1.0;

    R(0, 4) = -1.0;
    R(2, 4) = (local.lambda + 2 * local.mu) / local.lambda;

    R(0, 5) = -1.0;
    R(1, 5) = (local.lambda + 2 * local.mu) / local.lambda;

    if (testIfAcoustic(neighbor.mu)) {
      // scale for better condition number of R
      R(7, 6) = neighbor.lambda;
      R(8, 7) = neighbor.lambda;

      R(0, 8) = neighbor.lambda;
      R(1, 8) = neighbor.lambda;
      R(2, 8) = neighbor.lambda;
      R(6, 8) = -std::sqrt((neighbor.lambda + 2 * neighbor.mu) / neighbor.rho);
    } else {
      R(5, 6) = 0.5;
      R(8, 6) = -std::sqrt(neighbor.mu / neighbor.rho);

      R(3, 7) = 0.5;
      R(7, 7) = -std::sqrt(neighbor.mu / neighbor.rho);

      R(0, 8) = 1.0;
      R(6, 8) = -std::sqrt((neighbor.lambda + 2 * neighbor.mu) / neighbor.rho);
    }
    R(9, 9) = 1.0;

    Matrix1010 C = Matrix1010::Zero();
    C(0, 0) = local.lambda + 2.0 * local.mu;
    C(0, 1) = local.lambda;
    C(0, 2) = local.lambda;
    C(1, 0) = local.lambda;
    C(1, 1) = local.lambda + 2.0 * local.mu;
    C(1, 2) = local.lambda;
    C(2, 0) = local.lambda;
    C(2, 1) = local.lambda;
    C(2, 2) = local.lambda + 2.0 * local.mu;
    C(3, 3) = 2.0 * local.mu;
    C(4, 4) = 2.0 * local.mu;
    C(5, 5) = 2.0 * local.mu;
    C(6, 6) = 1;
    C(7, 7) = 1;
    C(8, 8) = 1;
    C(9, 9) = 1;
    // Convert to stress for free-surface condition
    Matrix1010 R_sig = (C * R).eval();

    if (faceType == FaceType::FreeSurface) {
      MaterialType materialtype =
          testIfAcoustic(local.mu) ? MaterialType::Acoustic : MaterialType::Elastic;
      getTransposedFreeSurfaceGodunovState(materialtype, QgodLocal, QgodNeighbor, R_sig);

      Matrix1010 Qgod = Matrix1010::Zero();
      for (unsigned i = 0; i < Qgod.cols(); ++i) {
        for (unsigned j = 0; j < Qgod.rows(); ++j) {
          Qgod(i, j) = -QgodLocal(j, i);
          Qgod(i, j) = QgodLocal(j, i);
        }
      }
      // Convert to back to strain with the free-surface constraints
      Matrix1010 Qgod_temp = Matrix1010::Zero();
      Qgod_temp = ((C.inverse() * Qgod) * C).eval();

      for (unsigned i = 0; i < Qgod.cols(); ++i) {
        for (unsigned j = 0; j < Qgod.rows(); ++j) {
          QgodLocal(i, j) = -Qgod_temp(j, i);
          QgodLocal(i, j) = Qgod_temp(j, i);
        }
      }

    } else {
      // Currently, QgodLocal and QgodNeighbor will not be used in the corrector step.
      // No need to account for heterogeneous effect
      // (in the same way as free-surf) here for Rusanov flux.
      Matrix1010 chi = Matrix1010::Zero();
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
      for (unsigned idx = 0; idx < godunov.cols(); ++idx) {
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

    // damage
    matT(9, 9) = 1;
    matTinv(9, 9) = 1;
  }

  static DamageMaterial getRotatedMaterialCoefficients(double rotationParameters[36],
                                                       DamageMaterial& material) {
    return material;
  }

  static void initializeSpecificLocalData(const DamageMaterial& material,
                                          double timeStepWidth,
                                          Vertex localVertices[4],
                                          real& localVolume,
                                          real localSurfaces[4],
                                          std::array<std::array<double, 3>, 4>& localNormal,
                                          std::array<std::array<double, 3>, 4>& localTangent1,
                                          std::array<std::array<double, 3>, 4>& localTangent2,
                                          DamageLocalData* localData) {
    // currently, it takes new memory. Later on can switch to pointers
    for (int i_v = 0; i_v < 4; i_v++) {
      localData->localVertices[i_v] = localVertices[i_v];
      localData->localSurfaces[i_v] = localSurfaces[i_v];
      for (int i_c = 0; i_c < 3; i_c++) {
        localData->localNormal[i_v][i_c] = localNormal[i_v][i_c];
        localData->localTangent1[i_v][i_c] = localTangent1[i_v][i_c];
        localData->localTangent2[i_v][i_c] = localTangent2[i_v][i_c];
      }
    }
    localData->localVolume = localVolume;

    double x[4];
    double y[4];
    double z[4];
    double gradXi[3];
    double gradEta[3];
    double gradZeta[3];

    for (unsigned vertex = 0; vertex < 4; ++vertex) {
      const VrtxCoords& coords = localData->localVertices[vertex].coords;
      x[vertex] = coords[0];
      y[vertex] = coords[1];
      z[vertex] = coords[2];
    }

    seissol::transformations::tetrahedronGlobalToReferenceJacobian(
        x, y, z, gradXi, gradEta, gradZeta);
    // localData->globalMeshId = meshId;
    for (unsigned int i_x = 0; i_x < 3; i_x++) {
      localData->gradXiEtaZeta[i_x][0] = gradXi[i_x];
      localData->gradXiEtaZeta[i_x][1] = gradEta[i_x];
      localData->gradXiEtaZeta[i_x][2] = gradZeta[i_x];
    }
  }

  // TODO (Zihua): currently, it uses the same input as the LocalData
  //  there might be some redundant variables that may
  //  do not need to be passed. Better remove for cleaning the code.
  static void initializeSpecificNeighborData(const DamageMaterial& material,
                                             DamageMaterial* matNeighbor,
                                             double timeStepWidth,
                                             Vertex localVertices[4],
                                             real& localVolume,
                                             real localSurfaces[4],
                                             std::array<std::array<double, 3>, 4>& localNormal,
                                             std::array<std::array<double, 3>, 4>& localTangent1,
                                             std::array<std::array<double, 3>, 4>& localTangent2,
                                             DamageNeighborData* localData) {
    // currently, it takes new memory. Later on can switch to pointers
    double cpMinus = std::sqrt((material.lambda0 + 2.0 * material.mu0) / material.rho);
    for (int i_v = 0; i_v < 4; i_v++) {
      double cpPlus =
          std::sqrt((matNeighbor[i_v].lambda0 + 2.0 * matNeighbor[i_v].mu0) / matNeighbor[i_v].rho);
      localData->maxWavespeeds[i_v] = std::max(cpMinus, cpPlus);
      localData->localVertices[i_v] = localVertices[i_v];
      localData->localSurfaces[i_v] = localSurfaces[i_v];
      for (int i_c = 0; i_c < 3; i_c++) {
        localData->localNormal[i_v][i_c] = localNormal[i_v][i_c];
        localData->localTangent1[i_v][i_c] = localTangent1[i_v][i_c];
        localData->localTangent2[i_v][i_c] = localTangent2[i_v][i_c];
      }
    }
    localData->localVolume = localVolume;

    double x[4];
    double y[4];
    double z[4];
    double gradXi[3];
    double gradEta[3];
    double gradZeta[3];

    for (unsigned vertex = 0; vertex < 4; ++vertex) {
      const VrtxCoords& coords = localData->localVertices[vertex].coords;
      x[vertex] = coords[0];
      y[vertex] = coords[1];
      z[vertex] = coords[2];
    }

    seissol::transformations::tetrahedronGlobalToReferenceJacobian(
        x, y, z, gradXi, gradEta, gradZeta);
    // localData->globalMeshId = meshId;
    for (unsigned int i_x = 0; i_x < 3; i_x++) {
      localData->gradXiEtaZeta[i_x][0] = gradXi[i_x];
      localData->gradXiEtaZeta[i_x][1] = gradEta[i_x];
      localData->gradXiEtaZeta[i_x][2] = gradZeta[i_x];
    }
  }

  static void getPlaneWaveOperator(
      const DamageMaterial& material,
      const double n[3],
      std::complex<double> mdata[DamageMaterial::NumQuantities * DamageMaterial::NumQuantities]) {
    getElasticPlaneWaveOperator(material, n, mdata);
  }

  template <typename T>
  static void getTransposedSourceCoefficientTensor(const DamageMaterial& material,
                                                   T& sourceMatrix) {}
};

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_DAMAGE_MODEL_ELASTICSETUP_H_
