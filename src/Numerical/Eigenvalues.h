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
#include <cmath>
#include <numeric>
#include <utils/logger.h>
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
  using RealT = typename Eigen::NumTraits<typename MatrixT::Scalar>::Real;
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
 * The eigenvalues are sorted by their real parts, and the eigenvectors have unit length
 * As a backend this function uses the eigen3 library, on the balanced matrix (see balance), as
 * LAPACK's xGEEV does
 * @param M: Dense matrix of size dim x dim, stored in column-major format
 * @param output: Reference to an Eigenpair to store the eigenvalue decomposition
 */
template <typename T, size_t Dim>
void computeEigenvaluesWithEigen3(std::array<std::complex<T>, Dim * Dim>& m,
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

  for (size_t i = 0; i < Dim; ++i) {
    output.values[i] = eigenvalues(sortedIndices[i], 0);
  }

  Matrix eigenvectors = scale.template cast<std::complex<T>>().asDiagonal() * ces.eigenvectors();
  for (size_t j = 0; j < Dim; ++j) {
    eigenvectors.col(j).normalize();
  }

  auto r = yateto::DenseTensorView<2, std::complex<T>>(output.vectors.data(), {Dim, Dim});
  for (size_t j = 0; j < Dim; ++j) {
    for (size_t i = 0; i < Dim; ++i) {
      r(i, j) = eigenvectors(i, sortedIndices[j]);
    }
  }
}
} // namespace seissol::eigenvalues

#ifdef USE_LAPACK
#include "FC.h"

#include <complex>

namespace seissol::eigenvalues {
/*
 * LAPACK distinguishes between double and single precision (i.e. zgeev for double,
 * and cgeev for single).
 */
extern "C" {
#define FC_zgeev FC_GLOBAL(zgeev, ZGEEV)
extern void FC_zgeev(char* jobVl,
                     char* jobVr,
                     int* n,
                     std::complex<double>* matA,
                     int* lda,
                     std::complex<double>* ev,
                     std::complex<double>* vl,
                     int* ldvl,
                     std::complex<double>* vr,
                     int* ldvr,
                     std::complex<double>* work,
                     int* lwork,
                     double* rwork,
                     int* info);
#define FC_cgeev FC_GLOBAL(cgeev, CGEEV)
extern void FC_cgeev(char* jobVl,
                     char* jobVr,
                     int* n,
                     std::complex<float>* matA,
                     int* lda,
                     std::complex<float>* ev,
                     std::complex<float>* vl,
                     int* ldvl,
                     std::complex<float>* vr,
                     int* ldvr,
                     std::complex<float>* work,
                     int* lwork,
                     float* rwork,
                     int* info);
}

template <typename T>
void callLapackEigenvalueRoutine(char* jobVl,
                                 char* jobVr,
                                 int* n,
                                 std::complex<T>* matA,
                                 int* lda,
                                 std::complex<T>* ev,
                                 std::complex<T>* vl,
                                 int* ldvl,
                                 std::complex<T>* vr,
                                 int* ldvr,
                                 std::complex<T>* work,
                                 int* lwork,
                                 T* rwork,
                                 int* info) {}

template <>
inline void callLapackEigenvalueRoutine(char* jobVl,
                                        char* jobVr,
                                        int* n,
                                        std::complex<double>* matA,
                                        int* lda,
                                        std::complex<double>* ev,
                                        std::complex<double>* vl,
                                        int* ldvl,
                                        std::complex<double>* vr,
                                        int* ldvr,
                                        std::complex<double>* work,
                                        int* lwork,
                                        double* rwork,
                                        int* info) {
  FC_zgeev(jobVl, jobVr, n, matA, lda, ev, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
}

template <>
inline void callLapackEigenvalueRoutine(char* jobVl,
                                        char* jobVr,
                                        int* n,
                                        std::complex<float>* matA,
                                        int* lda,
                                        std::complex<float>* ev,
                                        std::complex<float>* vl,
                                        int* ldvl,
                                        std::complex<float>* vr,
                                        int* ldvr,
                                        std::complex<float>* work,
                                        int* lwork,
                                        float* rwork,
                                        int* info) {
  FC_cgeev(jobVl, jobVr, n, matA, lda, ev, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
};

/**
 * Computes the eigenvalue decomposition of a dim x dim matrix of type T
 * The eigenvalues are sorted by their real parts
 * As a backend this function uses the LAPACK.
 * @param M: Dense matrix of size dim x dim, stored in column-major format
 * @param output: Reference to an Eigenpair to store the eigenvalue decomposition
 */
template <typename T, size_t Dim>
void computeEigenvaluesWithLapack(std::array<std::complex<T>, Dim * Dim>& m,
                                  Eigenpair<std::complex<T>, Dim>& output) {
  // set up lapack variables
  int n = Dim;
  int lda = Dim;
  int ldvl = Dim;
  int ldvr = Dim;
  int info = 0;
  int lwork = 2 * Dim;
  T rwork[2 * Dim];
  std::complex<T> w[Dim];
  std::complex<T> vl[Dim * Dim];
  std::complex<T> vr[Dim * Dim];
  std::complex<T> work[2 * Dim];
  char computeVectors = 'V';
  char dontComputeVectors = 'N';

  // lapack overrides matrix, so copy to auxiliary array
  std::complex<T> a[Dim * Dim];
  std::copy(m.begin(), m.end(), a);

  // Call LAPACK
  callLapackEigenvalueRoutine<T>(&dontComputeVectors,
                                 &computeVectors,
                                 &n,
                                 a,
                                 &lda,
                                 w,
                                 vl,
                                 &ldvl,
                                 vr,
                                 &ldvr,
                                 work,
                                 &lwork,
                                 rwork,
                                 &info);

  std::vector<size_t> sortedIndices(Dim);
  std::iota(sortedIndices.begin(), sortedIndices.end(), 0);
  std::sort(sortedIndices.begin(), sortedIndices.end(), [&w](size_t a, size_t b) {
    return w[a].real() < w[b].real();
  });

  for (size_t i = 0; i < Dim; ++i) {
    output.values[i] = w[sortedIndices[i]];
  }

  auto r = yateto::DenseTensorView<2, std::complex<T>>(output.vectors.data(), {Dim, Dim});
  for (size_t j = 0; j < Dim; ++j) {
    for (size_t i = 0; i < Dim; ++i) {
      const size_t sortedJ = sortedIndices[j];
      r(i, j) = vr[sortedJ * Dim + i];
    }
  }
}
} // namespace seissol::eigenvalues
#endif // USE_LAPACK

namespace seissol::eigenvalues {

/**
 * The libraries computeEigenvalues can use. LAPACK is linked into a build with a poroelastic
 * configuration only.
 */
enum class Backend { Eigen3, Lapack };

#ifdef USE_LAPACK
constexpr Backend LapackIfLinked = Backend::Lapack;
#else
constexpr Backend LapackIfLinked = Backend::Eigen3;
#endif

/**
 * Computes the eigenvalue decomposition of a dim x dim matrix of type T with the backend B
 * The eigenvalues are sorted by their real parts
 * @param M: Dense matrix of size dim x dim, stored in column-major format
 * @param output: Reference to an Eigenpair to store the eigenvalue decomposition
 */
template <Backend B, typename T, size_t Dim>
void computeEigenvalues(std::array<std::complex<T>, Dim * Dim>& m,
                        Eigenpair<std::complex<T>, Dim>& output) {
  if constexpr (B == Backend::Lapack) {
#ifdef USE_LAPACK
    computeEigenvaluesWithLapack(m, output);
#else
    static_assert(B != Backend::Lapack, "LAPACK is not linked into this build.");
#endif
  } else {
    computeEigenvaluesWithEigen3(m, output);
  }
}
} // namespace seissol::eigenvalues

#endif // SEISSOL_SRC_NUMERICAL_EIGENVALUES_H_
