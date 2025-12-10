// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_MODEL_COMMON_H_
#define SEISSOL_SRC_MODEL_COMMON_H_

#include "Equations/Datastructures.h"
#include "GeneratedCode/init.h"
#include "Geometry/MeshTools.h"
#include "Initializer/Typedefs.h"
#include "Model/CommonDatastructures.h"
#include "Numerical/Eigenvalues.h"
#include "Numerical/Transformation.h"

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <limits>
#include <utils/logger.h>

namespace seissol::model {

template <typename Cfg, typename = void>
struct MaterialSetup;

template <typename T>
constexpr bool testIfAcoustic(T mu) {
  return std::abs(mu) <= std::numeric_limits<T>::epsilon();
}

template <typename Cfg, typename Tmatrix>
void getTransposedCoefficientMatrix(const model::MaterialTT<Cfg>& material,
                                    unsigned dim,
                                    Tmatrix& matM) {
  MaterialSetup<Cfg>::getTransposedCoefficientMatrix(material, dim, matM);
}

template <typename Cfg, typename T>
void getTransposedSourceCoefficientTensor(const model::MaterialTT<Cfg>& material, T& matE) {
  MaterialSetup<Cfg>::getTransposedSourceCoefficientTensor(material, matE);
}

template <typename Cfg>
seissol::eigenvalues::Eigenpair<std::complex<double>, model::MaterialTT<Cfg>::NumQuantities>
    getEigenDecomposition(const model::MaterialTT<Cfg>& material, double zeroThreshold = 1e-7);

template <typename Cfg, typename Tloc, typename Tneigh>
void getTransposedGodunovState(const model::MaterialTT<Cfg>& local,
                               const model::MaterialTT<Cfg>& neighbor,
                               FaceType faceType,
                               Tloc& qGodLocal,
                               Tneigh& qGodNeighbor) {
  MaterialSetup<Cfg>::getTransposedGodunovState(local, neighbor, faceType, qGodLocal, qGodNeighbor);
}

// TODO: move to materials (currently not possible due to the acoustic-in-elastic "hack")
template <typename MaterialT, typename T, typename Tmatrix>
void getTransposedFreeSurfaceGodunovState(MaterialType materialtype,
                                          T& qGodLocal,
                                          T& qGodNeighbor,
                                          Tmatrix& matR);

template <typename Cfg>
void getPlaneWaveOperator(const model::MaterialTT<Cfg>& material,
                          const double n[3],
                          std::complex<double> mdata[model::MaterialTT<Cfg>::NumQuantities *
                                                     model::MaterialTT<Cfg>::NumQuantities]) {
  MaterialSetup<Cfg>::getPlaneWaveOperator(material, n, mdata);
}

template <typename Cfg>
void initializeSpecificLocalData(
    const model::MaterialTT<Cfg>& material,
    double timeStepWidth,
    typename model::MaterialTT<Cfg>::template LocalSpecificData<Cfg>* localData) {
  MaterialSetup<Cfg>::initializeSpecificLocalData(material, timeStepWidth, localData);
}

template <typename Cfg>
void initializeSpecificNeighborData(
    const model::MaterialTT<Cfg>& material,
    typename model::MaterialTT<Cfg>::template NeighborSpecificData<Cfg>* neighborData) {
  MaterialSetup<Cfg>::initializeSpecificNeighborData(material, neighborData);
}

/*
 * Calculates the so called Bond matrix. Anisotropic materials are characterized by
 * 21 different material parameters. Due to the directional dependence of anisotropic
 * materials the parameters are not independet of the choice of the coordinate system.
 * The Bond matrix transforms materials from one orthogonal coordinate system to
 * another one.
 * c.f. 10.1111/j.1365-246X.2007.03381.x
 */
void getBondMatrix(const VrtxCoords normal,
                   const VrtxCoords tangent1,
                   const VrtxCoords tangent2,
                   std::array<double, 36>& matN);

template <typename Cfg>
void getFaceRotationMatrix(const Eigen::Vector3d& normal,
                           const Eigen::Vector3d& tangent1,
                           const Eigen::Vector3d& tangent2,
                           typename init::T<Cfg>::view::type& matT,
                           typename init::Tinv<Cfg>::view::type& matTinv) {
  const VrtxCoords n = {normal(0), normal(1), normal(2)};
  const VrtxCoords s = {tangent1(0), tangent1(1), tangent1(2)};
  const VrtxCoords t = {tangent2(0), tangent2(1), tangent2(2)};
  getFaceRotationMatrix<Cfg>(n, s, t, matT, matTinv);
}

template <typename Cfg>
void getFaceRotationMatrix(const VrtxCoords normal,
                           const VrtxCoords tangent1,
                           const VrtxCoords tangent2,
                           typename init::T<Cfg>::view::type& matT,
                           typename init::Tinv<Cfg>::view::type& matTinv) {
  MaterialSetup<Cfg>::getFaceRotationMatrix(normal, tangent1, tangent2, matT, matTinv);
}

template <typename Cfg>
model::MaterialTT<Cfg>
    getRotatedMaterialCoefficients(const std::array<double, 36>& rotationParameters,
                                   model::MaterialTT<Cfg>& material) {
  return MaterialSetup<Cfg>::getRotatedMaterialCoefficients(rotationParameters, material);
}

template <typename Cfg>
void getElasticPlaneWaveOperator(
    const model::MaterialTT<Cfg>& material,
    const double n[3],
    std::complex<double>
        mdata[model::MaterialTT<Cfg>::NumQuantities * model::MaterialTT<Cfg>::NumQuantities]) {
  using MaterialT = model::MaterialTT<Cfg>;
  yateto::DenseTensorView<2, std::complex<double>> matM(
      mdata, {MaterialT::NumQuantities, MaterialT::NumQuantities});
  matM.setZero();

  double data[MaterialT::NumQuantities * MaterialT::NumQuantities];
  yateto::DenseTensorView<2, double> coeff(data,
                                           {MaterialT::NumQuantities, MaterialT::NumQuantities});

  for (unsigned d = 0; d < 3; ++d) {
    coeff.setZero();
    getTransposedCoefficientMatrix<Cfg>(material, d, coeff);

    for (unsigned i = 0; i < MaterialT::NumQuantities; ++i) {
      for (unsigned j = 0; j < MaterialT::NumQuantities; ++j) {
        matM(i, j) += n[d] * coeff(j, i);
      }
    }
  }
  coeff.setZero();
  getTransposedSourceCoefficientTensor<Cfg>(material, coeff);

  for (unsigned i = 0; i < MaterialT::NumQuantities; ++i) {
    for (unsigned j = 0; j < MaterialT::NumQuantities; ++j) {
      matM(i, j) -= std::complex<double>(0.0, coeff(j, i));
    }
  }
}

} // namespace seissol::model

template <typename T, typename Tmatrix, typename Tarray1, typename Tarray2>
void setBlocks(T qGodLocal, Tmatrix mS, Tarray1 tractionIndices, Tarray2 velocityIndices) {
  // set lower left block
  int col = 0;
  for (auto& t : tractionIndices) {
    int row = 0;
    for (auto& v : velocityIndices) {
      qGodLocal(t, v) = mS(row, col);
      row++;
    }
    col++;
  }

  // set lower right block
  for (auto& v : velocityIndices) {
    qGodLocal(v, v) = 1.0;
  }
}

template <typename Cfg>
seissol::eigenvalues::Eigenpair<std::complex<double>,
                                seissol::model::MaterialTT<Cfg>::NumQuantities>
    seissol::model::getEigenDecomposition(const model::MaterialTT<Cfg>& material,
                                          double zeroThreshold) {
  using Tmaterial = model::MaterialTT<Cfg>;
  std::array<std::complex<double>, Tmaterial::NumQuantities * Tmaterial::NumQuantities> dataAT;
  auto viewAT = yateto::DenseTensorView<2, std::complex<double>>(
      dataAT.data(), {Tmaterial::NumQuantities, Tmaterial::NumQuantities});
  getTransposedCoefficientMatrix<Cfg>(material, 0, viewAT);
  std::array<std::complex<double>, Tmaterial::NumQuantities * Tmaterial::NumQuantities> dataA;
  // transpose AT to get A
  for (std::size_t i = 0; i < Tmaterial::NumQuantities; i++) {
    for (std::size_t j = 0; j < Tmaterial::NumQuantities; j++) {
      dataA[i + Tmaterial::NumQuantities * j] = dataAT[Tmaterial::NumQuantities * i + j];
    }
  }
  seissol::eigenvalues::Eigenpair<std::complex<double>, Tmaterial::NumQuantities> eigenpair;

  seissol::eigenvalues::computeEigenvalues(dataA, eigenpair);
#ifndef NDEBUG
  using CMatrix =
      Eigen::Matrix<std::complex<double>, Tmaterial::NumQuantities, Tmaterial::NumQuantities>;
  using CVector = Eigen::Matrix<std::complex<double>, Tmaterial::NumQuantities, 1>;
  const CMatrix eigenvectors = CMatrix(eigenpair.vectors.data());
  const CVector eigenvalues = CVector(eigenpair.values.data());
  // check number of eigenvalues
  // also check that the imaginary parts are zero
  int evNeg = 0;
  int evPos = 0;
  for (std::size_t i = 0; i < Tmaterial::NumQuantities; ++i) {
    assert(std::abs(eigenvalues(i).imag()) < zeroThreshold);
    if (eigenvalues(i).real() < -zeroThreshold) {
      ++evNeg;
    } else if (eigenvalues(i).real() > zeroThreshold) {
      ++evPos;
    }
  }
  constexpr int ExpectedEv = Tmaterial::Type == MaterialType::Poroelastic ? 4 : 3;
  assert(evNeg == ExpectedEv);
  assert(evPos == ExpectedEv);

  // check whether eigensolver is good enough
  const CMatrix coeff(dataA.data());
  const CMatrix matrixMult = coeff * eigenvectors;
  CMatrix eigenvalueMatrix = CMatrix::Zero();
  for (std::size_t i = 0; i < Tmaterial::NumQuantities; i++) {
    eigenvalueMatrix(i, i) = eigenvalues(i);
  }
  const CMatrix vectorMult = eigenvectors * eigenvalueMatrix;
  const CMatrix diff = matrixMult - vectorMult;
  const double norm = diff.norm();

  std::stringstream messageStream;
  messageStream << "Residual " << norm << " is larger than " << zeroThreshold
                << ": Eigensolver is not accurate enough";
  assert(norm < zeroThreshold && messageStream.str().c_str());
#endif
  return eigenpair;
};

template <typename MaterialT, typename T, typename Tmatrix>
void seissol::model::getTransposedFreeSurfaceGodunovState(MaterialType materialtype,
                                                          T& qGodLocal,
                                                          T& qGodNeighbor,
                                                          Tmatrix& matR) {
  for (size_t i = 0; i < MaterialT::NumElasticQuantities; i++) {
    for (size_t j = 0; j < MaterialT::NumElasticQuantities; j++) {
      qGodNeighbor(i, j) = std::numeric_limits<double>::signaling_NaN();
    }
  }

  // matR == eigenvector matrix

  qGodLocal.setZero();
  switch (materialtype) {
  case MaterialType::Acoustic: {
    // Acoustic material only has one traction (=pressure) and one velocity comp.
    // relevant to the Riemann problem
    qGodLocal(0, 6) = -1 * matR(6, 0) * 1 / matR(0, 0); // mS
    qGodLocal(6, 6) = 1.0;
    break;
  }
  case MaterialType::Poroelastic: {
    using Matrix44 = Eigen::Matrix<double, 4, 4>;
    using Matrix64 = Eigen::Matrix<double, 6, 4>;

    const std::array<int, 4> tractionIndices = {0, 3, 5, 9};
    const std::array<int, 6> velocityIndices = {6, 7, 8, 10, 11, 12};
    const std::array<int, 4> columnIndices = {0, 1, 2, 3};
    const Matrix44 matR11 = matR(tractionIndices, columnIndices);
    const Matrix64 matR21 = matR(velocityIndices, columnIndices);
    const Matrix64 matS = (-(matR21 * matR11.inverse())).eval();
    setBlocks(qGodLocal, matS, tractionIndices, velocityIndices);
    break;
  }
  default: {
    const std::array<int, 3> tractionIndices = {0, 3, 5};
    const std::array<int, 3> velocityIndices = {6, 7, 8};
    using Matrix33 = Eigen::Matrix<double, 3, 3>;
    const Matrix33 matR11 = matR(tractionIndices, {0, 1, 2});
    const Matrix33 matR21 = matR(velocityIndices, {0, 1, 2});
    const Matrix33 matS = (-(matR21 * matR11.inverse())).eval();
    setBlocks(qGodLocal, matS, tractionIndices, velocityIndices);
    break;
  }
  }
}

#endif // SEISSOL_SRC_MODEL_COMMON_H_
