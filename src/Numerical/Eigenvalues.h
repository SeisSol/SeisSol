// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_NUMERICAL_EIGENVALUES_H_
#define SEISSOL_SRC_NUMERICAL_EIGENVALUES_H_

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <complex>
#include <numeric>
#include <utils/logger.h>
#include <vector>
#include <yateto.h>

namespace seissol::eigenvalues {

/**
 * Stores an eigenpair of a dim x dim matrix of type T, i.e. a vector of eigenvalues and a matrix of
 * eigenvectors
 */
template <typename T, size_t Dim>
struct Eigenpair {
  /**
   * vectors: Matrix of Eigenvectors in column-major format
   */
  std::array<T, Dim * Dim> vectors;
  /**
   * values: Vector of eigenvalues, in the same ordering as the eigenvectors
   */
  std::array<T, Dim> values;
  /**
   * returns: eigenvectors as Eigen3 matrix
   */
  Eigen::Matrix<T, Dim, Dim> getVectorsAsMatrix() {
    return Eigen::Matrix<T, Dim, Dim>(vectors.data());
  }
  /**
   * returns: eigenvalues as Eigen3 vector
   */
  Eigen::Matrix<T, Dim, 1> getValuesAsVector() { return Eigen::Matrix<T, Dim, 1>(values.data()); }
};

/// The diagonal of a diagonal matrix that scales the rows or columns of a MatrixT.
template <typename MatrixT>
using Scaling = Eigen::Matrix<typename Eigen::NumTraits<typename MatrixT::Scalar>::Real,
                              MatrixT::RowsAtCompileTime,
                              1>;

/**
 * Balances the square matrix m in place: with the diagonal matrix D of powers of two that this
 * returns, m becomes D^-1 m D, whose rows and columns have off-diagonal parts of about equal norm
 * (Parlett and Reinsch; LAPACK's xGEBAL without the permutations). Its eigenvalues are those of m,
 * and an eigenvector y of it is the eigenvector D y of m. The scaling is exact in floating point.
 *
 * An eigensolver computes eigenvalues and eigenvectors to a precision relative to the norm of the
 * matrix. The Jacobians of a solid hold stiffnesses of ~1e10 in the rows of the stresses and the
 * inverse of the density, ~1e-4, in the rows of the velocities. Unless such a matrix is balanced,
 * the speeds of its waves come out with relative errors of the order of 1e-3 wherever the waves
 * couple all components, as they do in an anisotropic medium.
 */
template <typename MatrixT>
Scaling<MatrixT> balance(MatrixT& m) {
  using RealT = Eigen::NumTraits<typename MatrixT::Scalar>::Real;
  constexpr RealT Radix = 2;
  constexpr RealT Sufficient = 0.95;

  Scaling<MatrixT> scale = Scaling<MatrixT>::Ones(m.rows());
  bool balanced = false;
  while (!balanced) {
    balanced = true;
    for (Eigen::Index i = 0; i < m.rows(); ++i) {
      RealT column = 0;
      RealT row = 0;
      for (Eigen::Index j = 0; j < m.rows(); ++j) {
        if (j != i) {
          column += std::abs(m(j, i));
          row += std::abs(m(i, j));
        }
      }
      // a row or column without off-diagonal entries has nothing to balance, and one with
      // entries that are not finite cannot be balanced
      if (column == 0 || row == 0 || !std::isfinite(column + row)) {
        continue;
      }
      const RealT sum = column + row;
      // the power of two f that brings the norms column * f and row / f closest together
      RealT factor = 1;
      while (column < row / Radix) {
        factor *= Radix;
        column *= Radix * Radix;
      }
      while (column > row * Radix) {
        factor /= Radix;
        column /= Radix * Radix;
      }
      if ((column + row) / factor < Sufficient * sum) {
        balanced = false;
        scale(i) *= factor;
        m.row(i) /= factor;
        m.col(i) *= factor;
      }
    }
  }
  return scale;
}

/**
 * Computes the eigenvalue decomposition of a dim x dim matrix of type T
 * The eigenvalues are sorted by their real parts, and those whose real parts agree to 1e-8 of the
 * largest modulus by their imaginary parts: the eigenvalues of modes that decay or stay without
 * travelling have no real part, and their order would otherwise be down to its rounding. The
 * eigenvectors have unit length, and their largest component is real and positive, as LAPACK's
 * xGEEV normalizes them: an eigensolver determines an eigenvector only up to a complex factor, and
 * this one fixes it, so that the eigenvector of a simple real eigenvalue of a real matrix is real.
 * The decomposition is the one of the eigen3 library, of the balanced matrix (see balance)
 * @param M: Dense matrix of size dim x dim, stored in column-major format
 * @param output: Reference to an Eigenpair to store the eigenvalue decomposition
 */
template <typename T, size_t Dim>
void computeEigenvalues(std::array<std::complex<T>, Dim * Dim>& m,
                        Eigenpair<std::complex<T>, Dim>& output) {
  using Matrix = Eigen::Matrix<std::complex<T>, Dim, Dim, Eigen::ColMajor>;
  Matrix op(m.data());
  const auto scale = balance(op);
  Eigen::ComplexEigenSolver<Matrix> ces;
  ces.compute(op);

  // sort eigenvalues so that we know which eigenvalue corresponds to which mode
  auto eigenvalues = ces.eigenvalues();
  std::vector<size_t> sortedIndices(Dim);
  std::iota(sortedIndices.begin(), sortedIndices.end(), 0);
  std::sort(sortedIndices.begin(), sortedIndices.end(), [&eigenvalues](size_t a, size_t b) {
    return eigenvalues[a].real() < eigenvalues[b].real();
  });
  T largestModulus = 0;
  for (size_t i = 0; i < Dim; ++i) {
    largestModulus = std::max(largestModulus, std::abs(eigenvalues(i)));
  }
  const T tolerance = static_cast<T>(1e-8) * largestModulus;
  for (size_t begin = 0; begin < Dim;) {
    size_t end = begin + 1;
    while (end < Dim &&
           eigenvalues[sortedIndices[end]].real() - eigenvalues[sortedIndices[end - 1]].real() <=
               tolerance) {
      ++end;
    }
    std::stable_sort(sortedIndices.begin() + begin,
                     sortedIndices.begin() + end,
                     [&eigenvalues](size_t a, size_t b) {
                       return eigenvalues[a].imag() < eigenvalues[b].imag();
                     });
    begin = end;
  }

  for (size_t i = 0; i < Dim; ++i) {
    output.values[i] = eigenvalues(sortedIndices[i], 0);
  }

  Matrix eigenvectors = scale.template cast<std::complex<T>>().asDiagonal() * ces.eigenvectors();
  for (size_t j = 0; j < Dim; ++j) {
    auto column = eigenvectors.col(j);
    column.normalize();
    Eigen::Index largest = 0;
    column.cwiseAbs().maxCoeff(&largest);
    column *= std::conj(column(largest)) / std::abs(column(largest));
  }

  auto r = yateto::DenseTensorView<2, std::complex<T>>(output.vectors.data(), {Dim, Dim});
  for (size_t j = 0; j < Dim; ++j) {
    for (size_t i = 0; i < Dim; ++i) {
      r(i, j) = eigenvectors(i, sortedIndices[j]);
    }
  }
}
} // namespace seissol::eigenvalues

#endif // SEISSOL_SRC_NUMERICAL_EIGENVALUES_H_
