#ifndef EIGENVALUES_H
#define EIGENVALUES_H

#include <numeric>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>

#include <yateto/TensorView.h>
#include <utils/logger.h>

#include "Kernels/precision.hpp"

namespace seissol::eigenvalues{

  /**
   * Stores an eigenpair of a dim x dim matrix of type T, i.e. a vector of eigenvalues and a matrix of eigenvectors
   */
  template<typename T, size_t dim>
  struct Eigenpair{
    /** 
     * vectors: Matrix of Eigenvectors in column-major format
     */
    std::array<T, dim*dim> vectors;
    /**
     * values: Vector of eigenvalues, in the same ordering as the eigenvectors
     */
    std::array<T, dim> values;
    /**
     * returns: eigenvectors as Eigen3 matrix
     */
     Eigen::Matrix<T, dim, dim> getVectorsAsMatrix() {
       return Eigen::Matrix<T, dim, dim>(vectors.data());
     }
    /**
     * returns: eigenvalues as Eigen3 vector
     */
    Eigen::Matrix<T, dim, 1> getValuesAsVector() {
      return Eigen::Matrix<T, dim, 1>(values.data());
    }
  };

  /**
   * Computes the eigenvalue decomposition of a dim x dim matrix of type T
   * The eigenvalues are sorted by their real parts
   * As a backend this function uses the eigen3 library
   * @param M: Dense matrix of size dim x dim, stored in column-major format
   * @param output: Reference to an Eigenpair to store the eigenvalue decomposition
   */
  template<typename T, size_t dim>
  void computeEigenvaluesWithEigen3(std::array<std::complex<T>, dim*dim>& M, Eigenpair<std::complex<T>, dim>& output) {
    using Matrix = Eigen::Matrix<std::complex<T>, dim, dim, Eigen::ColMajor>;
    Matrix op(M.data());
    Eigen::ComplexEigenSolver<Matrix> ces;
    ces.compute(op);

    //sort eigenvalues so that we know which eigenvalue corresponds to which mode
    auto eigenvalues = ces.eigenvalues();
    std::vector<size_t> sortedIndices(dim);
    std::iota(sortedIndices.begin(), sortedIndices.end(), 0);
    std::sort(sortedIndices.begin(), sortedIndices.end(), [&eigenvalues](size_t a, size_t b) {
      return eigenvalues[a].real() < eigenvalues[b].real();
    });

    for (size_t i = 0; i < dim; ++i) {
      output.values[i] = eigenvalues(sortedIndices[i],0);
    }

    auto eigenvectors = ces.eigenvectors();

    auto R = yateto::DenseTensorView<2,std::complex<T>>(output.vectors.data(), {dim, dim});
    for (size_t j = 0; j < dim; ++j) {
      for (size_t i = 0; i < dim; ++i) {
        R(i,j) = eigenvectors(i,sortedIndices[j]);
      }
    }
  }
} //namespace seissol::eigenvalues

#ifdef USE_POROELASTIC
#include <complex>
#include "FC.h"

namespace seissol::eigenvalues{
  /*
   * LAPACK distinguishes between double and single precision (i.e. zgeev for double,
   * and cgeev for single).
   */
  extern "C" {
#define FC_zgeev FC_GLOBAL(zgeev, ZGEEV)
    extern void FC_zgeev(char* jobVl, char* jobVr, int* n, std::complex<double>* A,
        int* lda, std::complex<double>* ev, std::complex<double>* vl, int* ldvl, 
        std::complex<double>* vr, int* ldvr, std::complex<double>* work, 
        int* lwork, double* rwork, int* info );
#define FC_cgeev FC_GLOBAL(cgeev, CGEEV)
    extern void FC_cgeev(char* jobVl, char* jobVr, int* n, std::complex<float>* A,
        int* lda, std::complex<float>* ev, std::complex<float>* vl, int* ldvl, 
        std::complex<float>* vr, int* ldvr, std::complex<float>* work, 
        int* lwork, float* rwork, int* info );
  }

template<typename T>
void callLapackEigenvalueRoutine(char* jobVl, char* jobVr, int* n, std::complex<T>* A,
        int* lda, std::complex<T>* ev, std::complex<T>* vl, int* ldvl, 
        std::complex<T>* vr, int* ldvr, std::complex<T>* work, 
        int* lwork, T* rwork, int* info) {}

template<>
inline void callLapackEigenvalueRoutine(char* jobVl, char* jobVr, int* n, std::complex<double>* A,
        int* lda, std::complex<double>* ev, std::complex<double>* vl, int* ldvl,
        std::complex<double>* vr, int* ldvr, std::complex<double>* work, 
        int* lwork, double* rwork, int* info) {
  FC_zgeev(jobVl, jobVr, n, A, lda, ev, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
}

template<>
inline void callLapackEigenvalueRoutine(char* jobVl, char* jobVr, int* n, std::complex<float>* A,
        int* lda, std::complex<float>* ev, std::complex<float>* vl, int* ldvl, 
        std::complex<float>* vr, int* ldvr, std::complex<float>* work, 
        int* lwork, float* rwork, int* info) {
  FC_cgeev(jobVl, jobVr, n, A, lda, ev, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
};

  /**
   * Computes the eigenvalue decomposition of a dim x dim matrix of type T
   * The eigenvalues are sorted by their real parts
   * As a backend this function uses the LAPACK.
   * @param M: Dense matrix of size dim x dim, stored in column-major format
   * @param output: Reference to an Eigenpair to store the eigenvalue decomposition
   */
  template<typename T, size_t dim>
  void computeEigenvaluesWithLapack(std::array<std::complex<T>, dim*dim>& M, Eigenpair<std::complex<T>, dim>& output) {
    //set up lapack variables
    int n = dim, lda = dim, ldvl = dim, ldvr = dim; 
    int info;
    int lwork = 2*dim;
    T rwork[2*dim];
    std::complex<T> w[dim], vl[dim*dim], vr[dim*dim];
    std::complex<T> work[2*dim];
    char computeVectors = 'V';
    char dontComputeVectors = 'N';

    //lapack overrides matrix, so copy to auxiliary array
    std::complex<T> a[dim*dim];
    std::copy(M.begin(), M.end(), a);

    //Call LAPACK
    callLapackEigenvalueRoutine<T>(&dontComputeVectors, &computeVectors, &n, a, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info );

    std::vector<size_t> sortedIndices(dim);
    std::iota(sortedIndices.begin(), sortedIndices.end(), 0);
    std::sort(sortedIndices.begin(), sortedIndices.end(), [&w](size_t a, size_t b) {
        return w[a].real() < w[b].real();
        });

    for (size_t i = 0; i < dim; ++i) {
      output.values[i] = w[sortedIndices[i]]; 
    }

    auto R = yateto::DenseTensorView<2,std::complex<T>>(output.vectors.data(), {dim, dim});
    for (size_t j = 0; j < dim; ++j) {
      for (size_t i = 0; i < dim; ++i) {
        size_t sorted_j = sortedIndices[j];
        R(i,j) = vr[sorted_j*dim+i];
      }
    }
  }
} //namespace seissol::eigenvalues
#endif //USE_POROELASTIC
#endif //EIGENVALUES_H
