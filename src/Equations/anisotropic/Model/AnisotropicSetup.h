// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_EQUATIONS_ANISOTROPIC_MODEL_ANISOTROPICSETUP_H_
#define SEISSOL_SRC_EQUATIONS_ANISOTROPIC_MODEL_ANISOTROPICSETUP_H_

#include "Datastructures.h"
#include "Equations/anisotropic/Model/IntegrationData.h"
#include "GeneratedCode/init.h"
#include "Kernels/Common.h"
#include "Model/Common.h"
#include "Numerical/Transformation.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace seissol::model {
using Matrix99 = Eigen::Matrix<double, 9, 9>;

template <>
struct MaterialSetup<AnisotropicMaterial> {
  template <typename T>
  static void
      getTransposedCoefficientMatrix(const AnisotropicMaterial& material, unsigned dim, T& matM) {
    matM.setZero();

    const auto rhoInv = 1.0 / material.rho;

    switch (dim) {
    case 0:
      matM(6, 0) = -material.c11;
      matM(7, 0) = -material.c16;
      matM(8, 0) = -material.c15;
      matM(6, 1) = -material.c12;
      matM(7, 1) = -material.c26;
      matM(8, 1) = -material.c25;
      matM(6, 2) = -material.c13;
      matM(7, 2) = -material.c36;
      matM(8, 2) = -material.c35;
      matM(6, 3) = -material.c16;
      matM(7, 3) = -material.c66;
      matM(8, 3) = -material.c56;
      matM(6, 4) = -material.c14;
      matM(7, 4) = -material.c46;
      matM(8, 4) = -material.c45;
      matM(6, 5) = -material.c15;
      matM(7, 5) = -material.c56;
      matM(8, 5) = -material.c55;
      matM(0, 6) = -rhoInv;
      matM(3, 7) = -rhoInv;
      matM(5, 8) = -rhoInv;
      break;

    case 1:
      matM(6, 0) = -material.c16;
      matM(7, 0) = -material.c12;
      matM(8, 0) = -material.c14;
      matM(6, 1) = -material.c26;
      matM(7, 1) = -material.c22;
      matM(8, 1) = -material.c24;
      matM(6, 2) = -material.c36;
      matM(7, 2) = -material.c23;
      matM(8, 2) = -material.c34;
      matM(6, 3) = -material.c66;
      matM(7, 3) = -material.c26;
      matM(8, 3) = -material.c46;
      matM(6, 4) = -material.c46;
      matM(7, 4) = -material.c24;
      matM(8, 4) = -material.c44;
      matM(6, 5) = -material.c56;
      matM(7, 5) = -material.c25;
      matM(8, 5) = -material.c45;
      matM(3, 6) = -rhoInv;
      matM(1, 7) = -rhoInv;
      matM(4, 8) = -rhoInv;
      break;

    case 2:
      matM(6, 0) = -material.c15;
      matM(7, 0) = -material.c14;
      matM(8, 0) = -material.c13;
      matM(6, 1) = -material.c25;
      matM(7, 1) = -material.c24;
      matM(8, 1) = -material.c23;
      matM(6, 2) = -material.c35;
      matM(7, 2) = -material.c34;
      matM(8, 2) = -material.c33;
      matM(6, 3) = -material.c56;
      matM(7, 3) = -material.c46;
      matM(8, 3) = -material.c36;
      matM(6, 4) = -material.c45;
      matM(7, 4) = -material.c44;
      matM(8, 4) = -material.c34;
      matM(6, 5) = -material.c55;
      matM(7, 5) = -material.c45;
      matM(8, 5) = -material.c35;
      matM(5, 6) = -rhoInv;
      matM(4, 7) = -rhoInv;
      matM(2, 8) = -rhoInv;
      break;

    default:
      break;
    }
  }

  static void getEigenBasisForAnisotropicMaterial(const AnisotropicMaterial& local,
                                                  const AnisotropicMaterial& neighbor,
                                                  Matrix99& matR) {
    using Matrix33 = Eigen::Matrix<double, 3, 3, Eigen::ColMajor>;
    using Matrix63 = Eigen::Matrix<double, 6, 3, Eigen::ColMajor>;

    // matR == eigenvector matrix

    /* Calculate Eigenvectors and Eigenvalues
     * We want to solve
     * /0  A\  /s\ = l /s\
     * \matR  0/  \u/     \u/
     * which is equivalent to
     * matR * A * u = l*l * u && s = 1/l A * u
     * Here A has shape 6x3 and matR has shape 3x6
     */
    Eigen::SelfAdjointEigenSolver<Matrix33> saes;

    std::array<double, 9> localRAData{};
    localRAData[0] = local.c11 / local.rho;
    localRAData[1] = local.c16 / local.rho;
    localRAData[2] = local.c15 / local.rho;
    localRAData[3] = local.c16 / local.rho;
    localRAData[4] = local.c66 / local.rho;
    localRAData[5] = local.c56 / local.rho;
    localRAData[6] = local.c15 / local.rho;
    localRAData[7] = local.c56 / local.rho;
    localRAData[8] = local.c55 / local.rho;
    const Matrix33 localRA(localRAData.data());
    saes.compute(localRA);
    auto eigenvaluesLocal = saes.eigenvalues();
    auto eigenvectorsLocal = saes.eigenvectors();

    std::array<double, 9> neighborRAData{};
    neighborRAData[0] = neighbor.c11 / neighbor.rho;
    neighborRAData[1] = neighbor.c16 / neighbor.rho;
    neighborRAData[2] = neighbor.c15 / neighbor.rho;
    neighborRAData[3] = neighbor.c16 / neighbor.rho;
    neighborRAData[4] = neighbor.c66 / neighbor.rho;
    neighborRAData[5] = neighbor.c56 / neighbor.rho;
    neighborRAData[6] = neighbor.c15 / neighbor.rho;
    neighborRAData[7] = neighbor.c56 / neighbor.rho;
    neighborRAData[8] = neighbor.c55 / neighbor.rho;
    const Matrix33 neighborRA(neighborRAData.data());
    saes.compute(neighborRA);
    auto eigenvaluesNeighbor = saes.eigenvalues();
    const auto eigenvectorsNeighbor = saes.eigenvectors();

    std::array<double, 18> localAData{};
    localAData[0] = -local.c11;
    localAData[1] = -local.c12;
    localAData[2] = -local.c13;
    localAData[3] = -local.c16;
    localAData[4] = -local.c14;
    localAData[5] = -local.c15;
    localAData[6] = -local.c16;
    localAData[7] = -local.c26;
    localAData[8] = -local.c36;
    localAData[9] = -local.c66;
    localAData[10] = -local.c46;
    localAData[11] = -local.c56;
    localAData[12] = -local.c15;
    localAData[13] = -local.c25;
    localAData[14] = -local.c35;
    localAData[15] = -local.c46;
    localAData[16] = -local.c45;
    localAData[17] = -local.c55;
    const Matrix63 localA(localAData.data());

    std::array<double, 18> neighborAData{};
    neighborAData[0] = -neighbor.c11;
    neighborAData[1] = -neighbor.c12;
    neighborAData[2] = -neighbor.c13;
    neighborAData[3] = -neighbor.c16;
    neighborAData[4] = -neighbor.c14;
    neighborAData[5] = -neighbor.c15;
    neighborAData[6] = -neighbor.c16;
    neighborAData[7] = -neighbor.c26;
    neighborAData[8] = -neighbor.c36;
    neighborAData[9] = -neighbor.c66;
    neighborAData[10] = -neighbor.c46;
    neighborAData[11] = -neighbor.c56;
    neighborAData[12] = -neighbor.c15;
    neighborAData[13] = -neighbor.c25;
    neighborAData[14] = -neighbor.c35;
    neighborAData[15] = -neighbor.c46;
    neighborAData[16] = -neighbor.c45;
    neighborAData[17] = -neighbor.c55;
    const Matrix63 neighborA(neighborAData.data());

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

    matR << localA * eigenvectorsLocal, nullSpaceEigenvectors, neighborA * eigenvectorsNeighbor,
        -eigenvectorsLocal * lambdaLocal, Matrix33::Zero(), eigenvectorsNeighbor * lambdaNeighbor;
  }

  static void getTransposedGodunovState(const AnisotropicMaterial& local,
                                        const AnisotropicMaterial& neighbor,
                                        FaceType faceType,
                                        init::QgodLocal::view::type& qGodLocal,
                                        init::QgodNeighbor::view::type& qGodNeighbor) {

    Matrix99 matR = Matrix99::Zero();
    getEigenBasisForAnisotropicMaterial(local, neighbor, matR);

    if (faceType == FaceType::FreeSurface) {
      getTransposedFreeSurfaceGodunovState(
          MaterialType::Anisotropic, qGodLocal, qGodNeighbor, matR);

    } else {
      Matrix99 chi = Matrix99::Zero();
      Matrix99 chiI = Matrix99::Zero();
      chi(0, 0) = 1.0;
      chi(1, 1) = 1.0;
      chi(2, 2) = 1.0;
      for (std::size_t i = 3; i < 9; ++i) {
        chiI(i, i) = 1.0;
      }

      auto matRT = matR.transpose();
      auto matRlu = matRT.lu();
      const auto godunov = matRlu.solve(chi * matRT).eval();
      const auto godunovI = matRlu.solve(chiI * matRT).eval();

      // QgodLocal = I - QgodNeighbor
      for (unsigned i = 0; i < godunov.cols(); ++i) {
        for (unsigned j = 0; j < godunov.rows(); ++j) {
          qGodLocal(i, j) = godunovI(i, j);
          qGodNeighbor(i, j) = godunov(i, j);
        }
      }
    }
  }

  static AnisotropicMaterial
      getRotatedMaterialCoefficients(const std::array<double, 36>& rotationParameters,
                                     AnisotropicMaterial& material) {
    AnisotropicMaterial rotatedMaterial;
    rotatedMaterial.rho = material.rho;
    using Matrix66 = Eigen::Matrix<double, 6, 6>;
    const Matrix66 mN(rotationParameters.data());
    Matrix66 mC = Matrix66();
    mC(0, 0) = material.c11;
    mC(0, 1) = material.c12;
    mC(0, 2) = material.c13;
    mC(0, 3) = material.c14;
    mC(0, 4) = material.c15;
    mC(0, 5) = material.c16;
    mC(1, 0) = material.c12;
    mC(1, 1) = material.c22;
    mC(1, 2) = material.c23;
    mC(1, 3) = material.c24;
    mC(1, 4) = material.c25;
    mC(1, 5) = material.c26;
    mC(2, 0) = material.c13;
    mC(2, 1) = material.c23;
    mC(2, 2) = material.c33;
    mC(2, 3) = material.c34;
    mC(2, 4) = material.c35;
    mC(2, 5) = material.c36;
    mC(3, 0) = material.c14;
    mC(3, 1) = material.c24;
    mC(3, 2) = material.c34;
    mC(3, 3) = material.c44;
    mC(3, 4) = material.c45;
    mC(3, 5) = material.c46;
    mC(4, 0) = material.c15;
    mC(4, 1) = material.c25;
    mC(4, 2) = material.c35;
    mC(4, 3) = material.c45;
    mC(4, 4) = material.c55;
    mC(4, 5) = material.c56;
    mC(5, 0) = material.c16;
    mC(5, 1) = material.c26;
    mC(5, 2) = material.c36;
    mC(5, 3) = material.c46;
    mC(5, 4) = material.c56;
    mC(5, 5) = material.c66;

    Matrix66 rotatedC = mN.transpose() * mC * mN;

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
                                          double timeStepWidth,
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
