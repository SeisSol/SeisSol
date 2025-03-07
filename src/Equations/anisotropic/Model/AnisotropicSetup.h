// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_EQUATIONS_ANISOTROPIC_MODEL_ANISOTROPICSETUP_H_
#define SEISSOL_SRC_EQUATIONS_ANISOTROPIC_MODEL_ANISOTROPICSETUP_H_

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Equations/anisotropic/Model/IntegrationData.h>

#include "Datastructures.h"
#include "Kernels/Common.h"
#include "Model/Common.h"
#include "Numerical/Transformation.h"
#include "generated_code/init.h"

namespace seissol::model {
using Matrix99 = Eigen::Matrix<double, 9, 9>;

template <>
struct MaterialSetup<AnisotropicMaterial> {
  template <typename T>
  static void
      getTransposedCoefficientMatrix(const AnisotropicMaterial& material, unsigned dim, T& M) {
    M.setZero();

    real rhoInv = 1.0 / material.rho;

    switch (dim) {
    case 0:
      M(6, 0) = -material.c11;
      M(7, 0) = -material.c16;
      M(8, 0) = -material.c15;
      M(6, 1) = -material.c12;
      M(7, 1) = -material.c26;
      M(8, 1) = -material.c25;
      M(6, 2) = -material.c13;
      M(7, 2) = -material.c36;
      M(8, 2) = -material.c35;
      M(6, 3) = -material.c16;
      M(7, 3) = -material.c66;
      M(8, 3) = -material.c56;
      M(6, 4) = -material.c14;
      M(7, 4) = -material.c46;
      M(8, 4) = -material.c45;
      M(6, 5) = -material.c15;
      M(7, 5) = -material.c56;
      M(8, 5) = -material.c55;
      M(0, 6) = -rhoInv;
      M(3, 7) = -rhoInv;
      M(5, 8) = -rhoInv;
      break;

    case 1:
      M(6, 0) = -material.c16;
      M(7, 0) = -material.c12;
      M(8, 0) = -material.c14;
      M(6, 1) = -material.c26;
      M(7, 1) = -material.c22;
      M(8, 1) = -material.c24;
      M(6, 2) = -material.c36;
      M(7, 2) = -material.c23;
      M(8, 2) = -material.c34;
      M(6, 3) = -material.c66;
      M(7, 3) = -material.c26;
      M(8, 3) = -material.c46;
      M(6, 4) = -material.c46;
      M(7, 4) = -material.c24;
      M(8, 4) = -material.c44;
      M(6, 5) = -material.c56;
      M(7, 5) = -material.c25;
      M(8, 5) = -material.c45;
      M(3, 6) = -rhoInv;
      M(1, 7) = -rhoInv;
      M(4, 8) = -rhoInv;
      break;

    case 2:
      M(6, 0) = -material.c15;
      M(7, 0) = -material.c14;
      M(8, 0) = -material.c13;
      M(6, 1) = -material.c25;
      M(7, 1) = -material.c24;
      M(8, 1) = -material.c23;
      M(6, 2) = -material.c35;
      M(7, 2) = -material.c34;
      M(8, 2) = -material.c33;
      M(6, 3) = -material.c56;
      M(7, 3) = -material.c46;
      M(8, 3) = -material.c36;
      M(6, 4) = -material.c45;
      M(7, 4) = -material.c44;
      M(8, 4) = -material.c34;
      M(6, 5) = -material.c55;
      M(7, 5) = -material.c45;
      M(8, 5) = -material.c35;
      M(5, 6) = -rhoInv;
      M(4, 7) = -rhoInv;
      M(2, 8) = -rhoInv;
      break;

    default:
      break;
    }
  }

  static void getEigenBasisForAnisotropicMaterial(const AnisotropicMaterial& local,
                                                  const AnisotropicMaterial& neighbor,
                                                  Matrix99& R) {
    using Matrix33 = Eigen::Matrix<double, 3, 3, Eigen::ColMajor>;
    using Matrix63 = Eigen::Matrix<double, 6, 3, Eigen::ColMajor>;

    /* Calculate Eigenvectors and Eigenvalues
     * We want to solve
     * /0  A\  /s\ = l /s\
     * \R  0/  \u/     \u/
     * which is equivalent to
     * R * A * u = l*l * u && s = 1/l A * u
     * Here A has shape 6x3 and R has shape 3x6
     */
    Eigen::SelfAdjointEigenSolver<Matrix33> saes;

    double raLocal[9];
    raLocal[0] = local.c11 / local.rho;
    raLocal[1] = local.c16 / local.rho;
    raLocal[2] = local.c15 / local.rho;
    raLocal[3] = local.c16 / local.rho;
    raLocal[4] = local.c66 / local.rho;
    raLocal[5] = local.c56 / local.rho;
    raLocal[6] = local.c15 / local.rho;
    raLocal[7] = local.c56 / local.rho;
    raLocal[8] = local.c55 / local.rho;
    Matrix33 RALocal(raLocal);
    saes.compute(RALocal);
    auto eigenvaluesLocal = saes.eigenvalues();
    auto eigenvectorsLocal = saes.eigenvectors();

    double raNeighbor[9];
    raNeighbor[0] = neighbor.c11 / neighbor.rho;
    raNeighbor[1] = neighbor.c16 / neighbor.rho;
    raNeighbor[2] = neighbor.c15 / neighbor.rho;
    raNeighbor[3] = neighbor.c16 / neighbor.rho;
    raNeighbor[4] = neighbor.c66 / neighbor.rho;
    raNeighbor[5] = neighbor.c56 / neighbor.rho;
    raNeighbor[6] = neighbor.c15 / neighbor.rho;
    raNeighbor[7] = neighbor.c56 / neighbor.rho;
    raNeighbor[8] = neighbor.c55 / neighbor.rho;
    Matrix33 RANeighbor(raNeighbor);
    saes.compute(RANeighbor);
    auto eigenvaluesNeighbor = saes.eigenvalues();
    auto eigenvectorsNeighbor = saes.eigenvectors();

    double aLocal[18];
    aLocal[0] = -local.c11;
    aLocal[1] = -local.c12;
    aLocal[2] = -local.c13;
    aLocal[3] = -local.c16;
    aLocal[4] = -local.c14;
    aLocal[5] = -local.c15;
    aLocal[6] = -local.c16;
    aLocal[7] = -local.c26;
    aLocal[8] = -local.c36;
    aLocal[9] = -local.c66;
    aLocal[10] = -local.c46;
    aLocal[11] = -local.c56;
    aLocal[12] = -local.c15;
    aLocal[13] = -local.c25;
    aLocal[14] = -local.c35;
    aLocal[15] = -local.c46;
    aLocal[16] = -local.c45;
    aLocal[17] = -local.c55;
    Matrix63 ALocal(aLocal);

    double aNeighbor[18];
    aNeighbor[0] = -neighbor.c11;
    aNeighbor[1] = -neighbor.c12;
    aNeighbor[2] = -neighbor.c13;
    aNeighbor[3] = -neighbor.c16;
    aNeighbor[4] = -neighbor.c14;
    aNeighbor[5] = -neighbor.c15;
    aNeighbor[6] = -neighbor.c16;
    aNeighbor[7] = -neighbor.c26;
    aNeighbor[8] = -neighbor.c36;
    aNeighbor[9] = -neighbor.c66;
    aNeighbor[10] = -neighbor.c46;
    aNeighbor[11] = -neighbor.c56;
    aNeighbor[12] = -neighbor.c15;
    aNeighbor[13] = -neighbor.c25;
    aNeighbor[14] = -neighbor.c35;
    aNeighbor[15] = -neighbor.c46;
    aNeighbor[16] = -neighbor.c45;
    aNeighbor[17] = -neighbor.c55;
    Matrix63 ANeighbor(aNeighbor);

    // remember that the eigenvalues of the complete system are the square roots
    // of the eigenvalues of the reduced system
    for (unsigned i = 0; i < 3; i++) {
      eigenvaluesLocal(i) = sqrt(eigenvaluesLocal(i));
    }
    auto lambdaLocal = Eigen::Matrix<double, 3, 3, Eigen::ColMajor>(eigenvaluesLocal.asDiagonal());
    for (unsigned i = 0; i < 3; i++) {
      eigenvaluesNeighbor(i) = sqrt(eigenvaluesNeighbor(i));
    }
    auto lambdaNeighbor = Matrix33(eigenvaluesNeighbor.asDiagonal());

    Matrix63 nullSpaceEigenvectors = Matrix63::Zero();
    nullSpaceEigenvectors(1, 0) = 1;
    nullSpaceEigenvectors(2, 1) = 1;
    nullSpaceEigenvectors(4, 2) = 1;

    R << ALocal * eigenvectorsLocal, nullSpaceEigenvectors, ANeighbor * eigenvectorsNeighbor,
        -eigenvectorsLocal * lambdaLocal, Matrix33::Zero(), eigenvectorsNeighbor * lambdaNeighbor;
  }

  static void getTransposedGodunovState(const AnisotropicMaterial& local,
                                        const AnisotropicMaterial& neighbor,
                                        FaceType faceType,
                                        init::QgodLocal::view::type& QgodLocal,
                                        init::QgodNeighbor::view::type& QgodNeighbor) {

    Matrix99 R = Matrix99::Zero();
    getEigenBasisForAnisotropicMaterial(local, neighbor, R);

    if (faceType == FaceType::FreeSurface) {
      getTransposedFreeSurfaceGodunovState(MaterialType::Anisotropic, QgodLocal, QgodNeighbor, R);

    } else {
      Matrix99 chi = Matrix99::Zero();
      chi(0, 0) = 1.0;
      chi(1, 1) = 1.0;
      chi(2, 2) = 1.0;

      const auto godunov = ((R * chi) * R.inverse()).eval();

      // QgodLocal = I - QgodNeighbor
      for (unsigned i = 0; i < QgodLocal.shape(1); ++i) {
        for (unsigned j = 0; j < QgodLocal.shape(0); ++j) {
          QgodLocal(i, j) = -godunov(j, i);
          QgodNeighbor(i, j) = godunov(j, i);
        }
      }
      for (unsigned idx = 0; idx < QgodLocal.shape(0) && idx < QgodLocal.shape(1); ++idx) {
        QgodLocal(idx, idx) += 1.0;
      }
    }
  }

  static AnisotropicMaterial getRotatedMaterialCoefficients(real rotationParameters[36],
                                                            AnisotropicMaterial& material) {
    AnisotropicMaterial rotatedMaterial;
    rotatedMaterial.rho = material.rho;
    using Matrix66 = Eigen::Matrix<real, 6, 6>;
    Matrix66 N = Matrix66(rotationParameters);
    Matrix66 C = Matrix66();
    C(0, 0) = material.c11;
    C(0, 1) = material.c12;
    C(0, 2) = material.c13;
    C(0, 3) = material.c14;
    C(0, 4) = material.c15;
    C(0, 5) = material.c16;
    C(1, 0) = material.c12;
    C(1, 1) = material.c22;
    C(1, 2) = material.c23;
    C(1, 3) = material.c24;
    C(1, 4) = material.c25;
    C(1, 5) = material.c26;
    C(2, 0) = material.c13;
    C(2, 1) = material.c23;
    C(2, 2) = material.c33;
    C(2, 3) = material.c34;
    C(2, 4) = material.c35;
    C(2, 5) = material.c36;
    C(3, 0) = material.c14;
    C(3, 1) = material.c24;
    C(3, 2) = material.c34;
    C(3, 3) = material.c44;
    C(3, 4) = material.c45;
    C(3, 5) = material.c46;
    C(4, 0) = material.c15;
    C(4, 1) = material.c25;
    C(4, 2) = material.c35;
    C(4, 3) = material.c45;
    C(4, 4) = material.c55;
    C(4, 5) = material.c56;
    C(5, 0) = material.c16;
    C(5, 1) = material.c26;
    C(5, 2) = material.c36;
    C(5, 3) = material.c46;
    C(5, 4) = material.c56;
    C(5, 5) = material.c66;

    Matrix66 rotatedC = N.transpose() * C * N;

    rotatedMaterial.c11 = rotatedC(0, 0);
    rotatedMaterial.c12 = rotatedC(0, 1);
    rotatedMaterial.c13 = rotatedC(0, 2);
    rotatedMaterial.c14 = rotatedC(0, 3);
    rotatedMaterial.c15 = rotatedC(0, 4);
    rotatedMaterial.c16 = rotatedC(0, 5);
    rotatedMaterial.c22 = rotatedC(1, 1);
    rotatedMaterial.c23 = rotatedC(1, 2);
    rotatedMaterial.c24 = rotatedC(1, 3);
    rotatedMaterial.c25 = rotatedC(1, 4);
    rotatedMaterial.c26 = rotatedC(1, 5);
    rotatedMaterial.c33 = rotatedC(2, 2);
    rotatedMaterial.c34 = rotatedC(2, 3);
    rotatedMaterial.c35 = rotatedC(2, 4);
    rotatedMaterial.c36 = rotatedC(2, 5);
    rotatedMaterial.c44 = rotatedC(3, 3);
    rotatedMaterial.c45 = rotatedC(3, 4);
    rotatedMaterial.c46 = rotatedC(3, 5);
    rotatedMaterial.c55 = rotatedC(4, 4);
    rotatedMaterial.c56 = rotatedC(4, 5);
    rotatedMaterial.c66 = rotatedC(5, 5);
    return rotatedMaterial;
  }

  static void initializeSpecificLocalData(const AnisotropicMaterial& material,
                                          real timeStepWidth,
                                          AnisotropicLocalData* localData) {}

  static void initializeSpecificNeighborData(const AnisotropicMaterial& material,
                                             AnisotropicNeighborData* localData) {}

  static void getPlaneWaveOperator(const AnisotropicMaterial& material,
                                   const double n[3],
                                   std::complex<double> mdata[AnisotropicMaterial::NumQuantities *
                                                              AnisotropicMaterial::NumQuantities]) {
    getElasticPlaneWaveOperator(material, n, mdata);
  }
  template <typename T>
  static void getTransposedSourceCoefficientTensor(const AnisotropicMaterial& material,
                                                   T& sourceMatrix) {}

  static void getFaceRotationMatrix(const VrtxCoords normal,
                                    const VrtxCoords tangent1,
                                    const VrtxCoords tangent2,
                                    init::T::view::type& matT,
                                    init::Tinv::view::type& matTinv) {
    ::seissol::model::getFaceRotationMatrix<ElasticMaterial>(
        normal, tangent1, tangent2, matT, matTinv);
  }
};

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_ANISOTROPIC_MODEL_ANISOTROPICSETUP_H_
