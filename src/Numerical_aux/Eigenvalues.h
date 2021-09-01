#ifndef EIGENVALUES_H
#define EIGENVALUES_H

#include <numeric>
#include <Eigen/Eigenvalues>

#include <yateto/TensorView.h>

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
  };

  /**
   * Computes the eigenvalue decomposition of a dim x dim matrix of type T
   * The eigenvalues are sorted by their real parts
   * As a backend this function uses the eigen3 library
   * @param M: Dense matrix of size dim x dim, stored in column-major format
   * @param output: Reference to an Eigenpair to store the eigenvalue decomposition
   */
  template<typename T, size_t dim>
  void computeEigenvaluesWithEigen3(std::array<T, dim*dim>& M, Eigenpair<T, dim>& output) {
    using Matrix = Eigen::Matrix<T, dim, dim, Eigen::ColMajor>;
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

    auto R = yateto::DenseTensorView<2,T>(output.vectors.data(), {dim, dim});
    for (size_t j = 0; j < dim; ++j) {
      for (size_t i = 0; i < dim; ++i) {
        R(i,j) = eigenvectors(i,sortedIndices[j]);
      }
    }
  }
} //namespace seissol::eigenvalues

#ifdef USE_POROELASTIC
#include <armadillo>
namespace seissol::eigenvalues{
  /**
   * Computes the eigenvalue decomposition of a dim x dim matrix of type T
   * The eigenvalues are sorted by their real parts
   * As a backend this function uses the armadillo library
   * @param M: Dense matrix of size dim x dim, stored in column-major format
   * @param output: Reference to an Eigenpair to store the eigenvalue decomposition
   */
  template<typename T, size_t dim>
  void computeEigenvaluesWithArmadillo(std::array<T, dim*dim>& M, Eigenpair<T, dim>& output) {
    using Matrix = typename arma::Mat<T>::template fixed<dim, dim>;
    using Vector = typename arma::Col<T>::template fixed<dim>;

    Matrix op(M.data());
    Matrix armaEigenvectors(arma::fill::zeros);
    Vector armaEigenvalues(arma::fill::zeros);
    arma::eig_gen(armaEigenvalues, armaEigenvectors, op, "balance");
    
    //sort eigenvalues
    std::array<size_t, dim> sortedIndices;
    std::iota(sortedIndices.begin(), sortedIndices.end(), 0);
    std::sort(sortedIndices.begin(), sortedIndices.end(), [&armaEigenvalues](size_t a, size_t b) {
       return armaEigenvalues[a].real() < armaEigenvalues[b].real();
    });

    for (size_t i = 0; i < dim; ++i) {
      output.values[i] = armaEigenvalues(sortedIndices[i]);
    }

    auto R = yateto::DenseTensorView<2,T>(output.vectors.data(), {dim, dim});
    for (size_t j = 0; j < dim; ++j) {
      for (size_t i = 0; i < dim; ++i) {
        R(i,j) = armaEigenvectors(i,sortedIndices[j]);
      }
    }
  }
} //namespace seissol::eigenvalues
#endif //USE_POROELASTIC
#endif //EIGENVALUES_H
