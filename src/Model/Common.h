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

template <typename MaterialT>
struct MaterialSetup;

template <typename T>
constexpr bool testIfAcoustic(T mu) {
  return std::abs(mu) <= std::numeric_limits<T>::epsilon();
}

template <typename Tmaterial, typename Tmatrix>
void getTransposedCoefficientMatrix(const Tmaterial& material, unsigned dim, Tmatrix& matM) {
  MaterialSetup<Tmaterial>::getTransposedCoefficientMatrix(material, dim, matM);
}

template <typename Tmaterial, typename T>
void getTransposedSourceCoefficientTensor(const Tmaterial& material, T& mE) {
  MaterialSetup<Tmaterial>::getTransposedSourceCoefficientTensor(material, mE);
}

template <typename Tmaterial>
seissol::eigenvalues::Eigenpair<std::complex<double>, seissol::model::MaterialT::NumQuantities>
    getEigenDecomposition(const Tmaterial& material, double zeroThreshold = 1e-7);

template <typename Tmaterial, typename Tloc, typename Tneigh>
void getTransposedGodunovState(const Tmaterial& local,
                               const Tmaterial& neighbor,
                               FaceType faceType,
                               Tloc& qGodLocal,
                               Tneigh& qGodNeighbor) {
  MaterialSetup<Tmaterial>::getTransposedGodunovState(
      local, neighbor, faceType, qGodLocal, qGodNeighbor);
}

// TODO: move to materials (currently not possible due to the acoustic-in-elastic "hack")
template <typename T, typename Tmatrix>
void getTransposedFreeSurfaceGodunovState(MaterialType materialtype,
                                          T& qGodLocal,
                                          T& qGodNeighbor,
                                          Tmatrix& matR);

template <typename T>
void getPlaneWaveOperator(const T& material,
                          const double n[3],
                          std::complex<double> mdata[T::NumQuantities * T::NumQuantities]) {
  MaterialSetup<T>::getPlaneWaveOperator(material, n, mdata);
}

template <typename T>
void initializeSpecificLocalData(const T& material,
                                 double timeStepWidth,
                                 typename T::LocalSpecificData* localData) {
  MaterialSetup<T>::initializeSpecificLocalData(material, timeStepWidth, localData);
}

template <typename T>
void initializeSpecificNeighborData(const T& material,
                                    typename T::NeighborSpecificData* neighborData) {
  MaterialSetup<T>::initializeSpecificNeighborData(material, neighborData);
}

/*
 * Calculates the so called Bond matrix. Anisotropic materials are characterized by
 * 21 different material parameters. Due to the directional dependence of anisotropic
 * materials the parameters are not independet of the choice of the coordinate system.
 * The Bond matrix transforms materials from one orthogonal coordinate system to
 * another one.
 * c.f. 10.1111/j.1365-246X.2007.03381.x
 * This method is not needed for isotropic materials.
 */
void getBondMatrix(const VrtxCoords normal,
                   const VrtxCoords tangent1,
                   const VrtxCoords tangent2,
                   std::array<double, 36>& matN);

template <typename MaterialT = seissol::model::MaterialT>
void getFaceRotationMatrix(const Eigen::Vector3d& normal,
                           const Eigen::Vector3d& tangent1,
                           const Eigen::Vector3d& tangent2,
                           init::T::view::type& matT,
                           init::Tinv::view::type& matTinv) {
  const VrtxCoords n = {normal(0), normal(1), normal(2)};
  const VrtxCoords s = {tangent1(0), tangent1(1), tangent1(2)};
  const VrtxCoords t = {tangent2(0), tangent2(1), tangent2(2)};
  getFaceRotationMatrix<MaterialT>(n, s, t, matT, matTinv);
}

template <typename MaterialT = seissol::model::MaterialT>
void getFaceRotationMatrix(const VrtxCoords normal,
                           const VrtxCoords tangent1,
                           const VrtxCoords tangent2,
                           init::T::view::type& matT,
                           init::Tinv::view::type& matTinv) {
  MaterialSetup<MaterialT>::getFaceRotationMatrix(normal, tangent1, tangent2, matT, matTinv);
}

template <typename MaterialT>
MaterialT getRotatedMaterialCoefficients(const std::array<double, 36>& rotationParameters,
                                         MaterialT& material) {
  return MaterialSetup<MaterialT>::getRotatedMaterialCoefficients(rotationParameters, material);
}

template <typename MaterialT>
void getElasticPlaneWaveOperator(
    const MaterialT& material,
    const double n[3],
    std::complex<double> mdata[MaterialT::NumQuantities * MaterialT::NumQuantities]) {
  yateto::DenseTensorView<2, std::complex<double>> matM(
      mdata, {MaterialT::NumQuantities, MaterialT::NumQuantities});
  matM.setZero();

  double data[MaterialT::NumQuantities * MaterialT::NumQuantities];
  yateto::DenseTensorView<2, double> coeff(data,
                                           {MaterialT::NumQuantities, MaterialT::NumQuantities});

  for (unsigned d = 0; d < 3; ++d) {
    coeff.setZero();
    getTransposedCoefficientMatrix(material, d, coeff);

    for (unsigned i = 0; i < MaterialT::NumQuantities; ++i) {
      for (unsigned j = 0; j < MaterialT::NumQuantities; ++j) {
        matM(i, j) += n[d] * coeff(j, i);
      }
    }
  }
  coeff.setZero();
  getTransposedSourceCoefficientTensor(material, coeff);

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

template <typename Tmaterial>
seissol::eigenvalues::Eigenpair<std::complex<double>, seissol::model::MaterialT::NumQuantities>
    seissol::model::getEigenDecomposition(const Tmaterial& material,
                                          [[maybe_unused]] double zeroThreshold) {
  std::array<std::complex<double>,
             seissol::model::MaterialT::NumQuantities * seissol::model::MaterialT::NumQuantities>
      dataAT;
  auto viewAT = yateto::DenseTensorView<2, std::complex<double>>(
      dataAT.data(),
      {seissol::model::MaterialT::NumQuantities, seissol::model::MaterialT::NumQuantities});
  getTransposedCoefficientMatrix(material, 0, viewAT);
  std::array<std::complex<double>,
             seissol::model::MaterialT::NumQuantities * seissol::model::MaterialT::NumQuantities>
      dataA;
  // transpose dataAT to get dataA
  for (std::size_t i = 0; i < seissol::model::MaterialT::NumQuantities; i++) {
    for (std::size_t j = 0; j < seissol::model::MaterialT::NumQuantities; j++) {
      dataA[i + seissol::model::MaterialT::NumQuantities * j] =
          dataAT[seissol::model::MaterialT::NumQuantities * i + j];
    }
  }
  seissol::eigenvalues::Eigenpair<std::complex<double>, seissol::model::MaterialT::NumQuantities>
      eigenpair;

  seissol::eigenvalues::computeEigenvalues(dataA, eigenpair);
#ifndef NDEBUG
  using CMatrix = Eigen::Matrix<std::complex<double>,
                                seissol::model::MaterialT::NumQuantities,
                                seissol::model::MaterialT::NumQuantities>;
  using CVector = Eigen::Matrix<std::complex<double>, seissol::model::MaterialT::NumQuantities, 1>;
  const CMatrix eigenvectors = CMatrix(eigenpair.vectors.data());
  const CVector eigenvalues = CVector(eigenpair.values.data());
  // check number of eigenvalues
  // also check that the imaginary parts are zero
  int evNeg = 0;
  int evPos = 0;
  for (std::size_t i = 0; i < seissol::model::MaterialT::NumQuantities; ++i) {
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
  for (std::size_t i = 0; i < seissol::model::MaterialT::NumQuantities; i++) {
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

template <typename T, typename Tmatrix>
void seissol::model::getTransposedFreeSurfaceGodunovState(MaterialType materialtype,
                                                          T& qGodLocal,
                                                          T& qGodNeighbor,
                                                          Tmatrix& matR) {
  for (size_t i = 0; i < seissol::model::MaterialT::NumElasticQuantities; i++) {
    for (size_t j = 0; j < seissol::model::MaterialT::NumElasticQuantities; j++) {
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
