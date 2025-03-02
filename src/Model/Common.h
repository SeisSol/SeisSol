// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_MODEL_COMMON_H_
#define SEISSOL_SRC_MODEL_COMMON_H_

#include <Eigen/Dense>
#include <Equations/Datastructures.h>

#include "Geometry/MeshTools.h"
#include "Initializer/Typedefs.h"
#include "Numerical/Eigenvalues.h"
#include "Numerical/Transformation.h"
#include "generated_code/init.h"
#include "utils/logger.h"

#include "Model/CommonDatastructures.h"

namespace seissol::model {

template <typename MaterialT>
struct MaterialSetup;

bool testIfAcoustic(real mu);

template <typename Tmaterial, typename Tmatrix>
void getTransposedCoefficientMatrix(const Tmaterial& material, unsigned dim, Tmatrix& M) {
  MaterialSetup<Tmaterial>::getTransposedCoefficientMatrix(material, dim, M);
}

template <typename Tmaterial, typename T>
void getTransposedSourceCoefficientTensor(const Tmaterial& material, T& E) {
  MaterialSetup<Tmaterial>::getTransposedSourceCoefficientTensor(material, E);
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
                                          Tmatrix& R);

template <typename T>
void getPlaneWaveOperator(const T& material,
                          const double n[3],
                          std::complex<double> mdata[T::NumQuantities * T::NumQuantities]) {
  MaterialSetup<T>::getPlaneWaveOperator(material, n, mdata);
}

template <typename T>
void initializeSpecificLocalData(const T& material,
                                 real timeStepWidth,
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
 */
void getBondMatrix(const VrtxCoords normal,
                   const VrtxCoords tangent1,
                   const VrtxCoords tangent2,
                   real* matN);

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
MaterialT getRotatedMaterialCoefficients(real rotationParameters[36], MaterialT& material) {
  return MaterialSetup<MaterialT>::getRotatedMaterialCoefficients(rotationParameters, material);
}

template <typename MaterialT>
void getElasticPlaneWaveOperator(
    const MaterialT& material,
    const double n[3],
    std::complex<double> mdata[MaterialT::NumQuantities * MaterialT::NumQuantities]) {
  yateto::DenseTensorView<2, std::complex<double>> M(
      mdata, {MaterialT::NumQuantities, MaterialT::NumQuantities});
  M.setZero();

  double data[MaterialT::NumQuantities * MaterialT::NumQuantities];
  yateto::DenseTensorView<2, double> coeff(data,
                                           {MaterialT::NumQuantities, MaterialT::NumQuantities});

  for (unsigned d = 0; d < 3; ++d) {
    coeff.setZero();
    getTransposedCoefficientMatrix(material, d, coeff);

    for (unsigned i = 0; i < MaterialT::NumQuantities; ++i) {
      for (unsigned j = 0; j < MaterialT::NumQuantities; ++j) {
        M(i, j) += n[d] * coeff(j, i);
      }
    }
  }
  coeff.setZero();
  getTransposedSourceCoefficientTensor(material, coeff);

  for (unsigned i = 0; i < MaterialT::NumQuantities; ++i) {
    for (unsigned j = 0; j < MaterialT::NumQuantities; ++j) {
      M(i, j) -= std::complex<real>(0.0, coeff(j, i));
    }
  }
}

} // namespace seissol::model

template <typename T, typename Tmatrix, typename Tarray1, typename Tarray2>
void setBlocks(T qGodLocal, Tmatrix S, Tarray1 tractionIndices, Tarray2 velocityIndices) {
  // set lower left block
  int col = 0;
  for (auto& t : tractionIndices) {
    int row = 0;
    for (auto& v : velocityIndices) {
      qGodLocal(t, v) = S(row, col);
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
    seissol::model::getEigenDecomposition(const Tmaterial& material, double zeroThreshold) {
  std::array<std::complex<double>,
             seissol::model::MaterialT::NumQuantities * seissol::model::MaterialT::NumQuantities>
      AT;
  auto ATView = yateto::DenseTensorView<2, std::complex<double>>(
      AT.data(),
      {seissol::model::MaterialT::NumQuantities, seissol::model::MaterialT::NumQuantities});
  getTransposedCoefficientMatrix(material, 0, ATView);
  std::array<std::complex<double>,
             seissol::model::MaterialT::NumQuantities * seissol::model::MaterialT::NumQuantities>
      A;
  // transpose AT to get A
  for (int i = 0; i < seissol::model::MaterialT::NumQuantities; i++) {
    for (int j = 0; j < seissol::model::MaterialT::NumQuantities; j++) {
      A[i + seissol::model::MaterialT::NumQuantities * j] =
          AT[seissol::model::MaterialT::NumQuantities * i + j];
    }
  }
  seissol::eigenvalues::Eigenpair<std::complex<double>, seissol::model::MaterialT::NumQuantities>
      eigenpair;

  seissol::eigenvalues::computeEigenvalues(A, eigenpair);
#ifndef NDEBUG
  using CMatrix = Eigen::Matrix<std::complex<double>,
                                seissol::model::MaterialT::NumQuantities,
                                seissol::model::MaterialT::NumQuantities>;
  using CVector = Eigen::Matrix<std::complex<double>, seissol::model::MaterialT::NumQuantities, 1>;
  CMatrix eigenvectors = CMatrix(eigenpair.vectors.data());
  CVector eigenvalues = CVector(eigenpair.values.data());
  // check number of eigenvalues
  // also check that the imaginary parts are zero
  int evNeg = 0;
  int evPos = 0;
  for (int i = 0; i < seissol::model::MaterialT::NumQuantities; ++i) {
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
  CMatrix coeff(A.data());
  const CMatrix matrixMult = coeff * eigenvectors;
  CMatrix eigenvalueMatrix = CMatrix::Zero();
  for (size_t i = 0; i < seissol::model::MaterialT::NumQuantities; i++) {
    eigenvalueMatrix(i, i) = eigenvalues(i);
  }
  const CMatrix vectorMult = eigenvectors * eigenvalueMatrix;
  const CMatrix diff = matrixMult - vectorMult;
  const double norm = diff.norm();

  std::stringstream messageStream;
  messageStream << "Residual " << norm << " is larger than " << zeroThreshold
                << ": Eigensolver is not accurate enough";
  assert((messageStream.str().c_str(), norm < zeroThreshold));
#endif
  return eigenpair;
};

template <typename T, typename Tmatrix>
void seissol::model::getTransposedFreeSurfaceGodunovState(MaterialType materialtype,
                                                          T& qGodLocal,
                                                          T& qGodNeighbor,
                                                          Tmatrix& R) {
  for (size_t i = 0; i < seissol::model::MaterialT::NumElasticQuantities; i++) {
    for (size_t j = 0; j < seissol::model::MaterialT::NumElasticQuantities; j++) {
      qGodNeighbor(i, j) = std::numeric_limits<double>::signaling_NaN();
    }
  }

  qGodLocal.setZero();
  switch (materialtype) {
  case MaterialType::Acoustic: {
    // Acoustic material only has one traction (=pressure) and one velocity comp.
    // relevant to the Riemann problem
    qGodLocal(0, 6) = -1 * R(6, 0) * 1 / R(0, 0); // S
    qGodLocal(6, 6) = 1.0;
    break;
  }
  case MaterialType::Poroelastic: {
    using Matrix44 = Eigen::Matrix<double, 4, 4>;
    using Matrix64 = Eigen::Matrix<double, 6, 4>;

    std::array<int, 4> tractionIndices = {0, 3, 5, 9};
    std::array<int, 6> velocityIndices = {6, 7, 8, 10, 11, 12};
    std::array<int, 4> columnIndices = {0, 1, 2, 3};
    Matrix44 R11 = R(tractionIndices, columnIndices);
    Matrix64 R21 = R(velocityIndices, columnIndices);
    Matrix64 S = (-(R21 * R11.inverse())).eval();
    setBlocks(qGodLocal, S, tractionIndices, velocityIndices);
    break;
  }
  default: {
    std::array<int, 3> tractionIndices = {0, 3, 5};
    std::array<int, 3> velocityIndices = {6, 7, 8};
    using Matrix33 = Eigen::Matrix<double, 3, 3>;
    Matrix33 R11 = R(tractionIndices, {0, 1, 2});
    Matrix33 R21 = R(velocityIndices, {0, 1, 2});
    Matrix33 S = (-(R21 * R11.inverse())).eval();
    setBlocks(qGodLocal, S, tractionIndices, velocityIndices);
    break;
  }
  }
}

#endif // SEISSOL_SRC_MODEL_COMMON_H_
