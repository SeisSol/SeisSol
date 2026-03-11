// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_NUMERICAL_BASISFUNCTION_H_
#define SEISSOL_SRC_NUMERICAL_BASISFUNCTION_H_

#include "Common/Constants.h"
#include "Functions.h"
#include "Transformation.h"

#include <cmath>
#include <numeric>
#include <type_traits>
#include <vector>

namespace seissol::basisFunction {

//------------------------------------------------------------------------------

/**
 * Functor to sample Basis functions for a given point in the reference
 * tetrahedron.
 * @param T the type to calculate internally.
 */
template <class T>
class BasisFunctionGenerator {
  private:
  T xi_, eta_, zeta_;

  /**
   * Function to sample a Jacobi Polynomial.
   * Source : https://github.com/pjabardo/Jacobi.jl/blob/master/src/jac_poly.jl
   * @param x Sampling point
   */
  static T sampleJacobiPolynomial(T x, std::uint32_t n, std::uint32_t a, std::uint32_t b) {
    return seissol::functions::JacobiP(n, a, b, x);
  }

  public:
  /**
   * Constructs a BasisFunctionGenerator which fixes it to a specific point
   * in the reference tetrahedron.
   * @param eta The eta coordinate in the reference terahedron.
   * @param zeta The zeta coordinate in the reference terahedron.
   * @param xi The xi coordinate in the reference terahedron.
   */
  BasisFunctionGenerator(T xi, T eta, T zeta) : xi_(xi), eta_(eta), zeta_(zeta) {}

  /**
   * Allows functor style call. Generates the sampled Basis function given
   * polynomial information for the generation of the Basis function.
   * This algorithm uses the calculation method from the paper :
   * "Seismic Wave Simulation for Complex Rheologies on Unstructured Meshes"
   * by Josep de la Puente
   * @param i Polynomial index information
   * @param j Polynomial index information
   * @param k Polynomial index information
   */
  T operator()(std::uint32_t i, std::uint32_t j, std::uint32_t k) const {
    return static_cast<T>(functions::TetraDubinerP({i, j, k}, {xi_, eta_, zeta_}));
  }
};
//------------------------------------------------------------------------------

/**
 * Functor to sample Basis function derivatives for a given point in the reference
 * tetrahedron.
 * @param T the type to calculate internally.
 */
template <class T>
class BasisFunctionDerivativeGenerator {
  private:
  T xi_, eta_, zeta_;

  public:
  /**
   * Constructs a BasisFunctionDerivativeGenerator which fixes it to a specific point
   * in the reference tetrahedron.
   * @param eta The eta coordinate in the reference terahedron.
   * @param zeta The zeta coordinate in the reference terahedron.
   * @param xi The xi coordinate in the reference terahedron.
   */
  BasisFunctionDerivativeGenerator(T xi, T eta, T zeta) : xi_(xi), eta_(eta), zeta_(zeta) {}

  /**
   * Allows functor style call.
   * @param i Polynomial index information
   * @param j Polynomial index information
   * @param k Polynomial index information
   */
  std::array<T, 3> operator()(std::uint32_t i, std::uint32_t j, std::uint32_t k) const {
    std::array<double, 3> gradEvaluated =
        functions::gradTetraDubinerP({i, j, k}, {xi_, eta_, zeta_});
    return {static_cast<T>(gradEvaluated[0]),
            static_cast<T>(gradEvaluated[1]),
            static_cast<T>(gradEvaluated[2])};
  }
};

inline std::uint32_t basisFunctionsForOrder(std::uint32_t order) {
  return (order * (order + 1) * (order + 2)) / 6;
}

/**
 * This class represents a vector Basis functions sampled at a specific point.
 * @param T denotes the type to calculate internally.
 */
template <class T>
class SampledBasisFunctions {
  private:
  static_assert(std::is_arithmetic_v<T>, "Type T for SampledBasisFunctions must be arithmetic.");

  /** The basis function samples */
  std::vector<T> data_;

  public:
  SampledBasisFunctions() = default;
  /**
   * Constructor to generate the sampled basis functions of given order
   * and at a given point in the reference tetrahedron.
   * @param order The order of the computation. It determines how many Basis
   * functions are generated. @see getBasisFunctionsPerOrder
   * @param eta The eta coordinate in the reference tetrahedron.
   * @param zeta The zeta coordinate in the reference tetrahedron.
   * @param xi The xi coordinate in the reference tetrahedron.
   */
  SampledBasisFunctions(std::uint32_t order, T xi, T eta, T zeta)
      : data_(basisFunctionsForOrder(order)) {
    const BasisFunctionGenerator<T> gen(xi, eta, zeta);

    std::uint32_t i = 0;
    for (std::uint32_t ord = 0; ord < order; ord++) {
      for (std::uint32_t k = 0; k <= ord; k++) {
        for (std::uint32_t j = 0; j <= ord - k; j++) {
          data_[i++] = gen(ord - j - k, j, k);
        }
      }
    }
  }

  /**
   * Function to evaluate the samples by multiplying the sampled Basis
   * function with its coefficient and summing up the products.
   * res = c0 * bf0 + c1 * bf1 + ... + cn * bfn
   * @param coeffIter the const iterator to read the coefficients
   */
  template <class ConstIterator>
  T evalWithCoeffs(ConstIterator coeffIter) const {
    return std::inner_product(data_.begin(), data_.end(), coeffIter, static_cast<T>(0));
  }

  /**
   * Returns the amount of Basis functions this class represents.
   */
  [[nodiscard]] auto size() const { return data_.size(); }

  /**
   * Returns the basis function vector directly.
   */
  [[nodiscard]] const std::vector<T>& data() const { return data_; }
};

//------------------------------------------------------------------------------

/**
 * This class represents the derivatives of vector Basis functions sampled at a specific point.
 * @param T denotes the type to calculate internally.
 */
template <class T>
struct SampledBasisFunctionDerivatives {
  private:
  static_assert(std::is_arithmetic_v<T>, "Type T for SampledBasisFunctions must be arithmetic.");

  /**
   * The basis function derivative samples w.r.t. the three spatial dimension
   * Use DenseTensorView to access data
   */
  std::vector<T> data_;

  public:
  SampledBasisFunctionDerivatives() = default;
  /**
   * Constructor to generate the sampled basis functions of given order
   * and at a given point in the reference tetrahedron.
   * @param order The order of the computation. It determines how many Basis
   * functions are generated. @see getBasisFunctionsPerOrder
   * @param eta The eta coordinate in the reference tetrahedron.
   * @param zeta The zeta coordinate in the reference tetrahedron.
   * @param xi The xi coordinate in the reference tetrahedron.
   */
  SampledBasisFunctionDerivatives(std::uint32_t order, T xi, T eta, T zeta)
      : data_(3 * basisFunctionsForOrder(order)) {
    const BasisFunctionDerivativeGenerator<T> gen(xi, eta, zeta);
    const auto funs = basisFunctionsForOrder(order);

    std::uint32_t i = 0;
    for (std::uint32_t ord = 0; ord < order; ord++) {
      for (std::uint32_t k = 0; k <= ord; k++) {
        for (std::uint32_t j = 0; j <= ord - k; j++) {
          const auto derivatives = gen(ord - j - k, j, k);
          for (std::uint32_t direction = 0; direction < 3; direction++) {
            data_[i + funs * direction] = derivatives[direction];
          }
          i++;
        }
      }
    }
  }

  /**
   * After a call to the constructor data contains the sampled derivatives
   * w.r.t. xi, eta, zeta on the reference triangle. Use this function to
   * transform the derivatives to derivatives w.r.t. to x, y, z in a physical
   * tetrahedron.
   * @param coords coords[i] contains the 3 coordinates of the ith vertex of the
   * physical tetrahedron.
   */
  void transformToGlobalCoordinates(const double* coords[Cell::NumVertices]) {
    std::array<double, Cell::NumVertices> xCoords{};
    std::array<double, Cell::NumVertices> yCoords{};
    std::array<double, Cell::NumVertices> zCoords{};
    for (size_t i = 0; i < Cell::NumVertices; ++i) {
      xCoords[i] = coords[i][0];
      yCoords[i] = coords[i][1];
      zCoords[i] = coords[i][2];
    }

    std::array<double, 3> gradXi{};
    std::array<double, 3> gradEta{};
    std::array<double, 3> gradZeta{};

    seissol::transformations::tetrahedronGlobalToReferenceJacobian(xCoords.data(),
                                                                   yCoords.data(),
                                                                   zCoords.data(),
                                                                   gradXi.data(),
                                                                   gradEta.data(),
                                                                   gradZeta.data());
    std::vector<T> oldData = data_;
    const auto funs = data_.size() / 3;

    for (size_t i = 0; i < funs; ++i) {
      for (size_t direction = 0; direction < 3; ++direction) {
        // dpsi / di = dphi / dxi * dxi / di + dphi / deta * deta / di + dphi / dzeta * dzeta / di
        data_[i + funs * direction] = oldData[i + 0 * funs] * gradXi[direction] +
                                      oldData[i + 1 * funs] * gradEta[direction] +
                                      oldData[i + 2 * funs] * gradZeta[direction];
      }
    }
  }

  /**
   * Returns the amount of Basis functions this class represents.
   */
  [[nodiscard]] auto size() const { return data_.size(); }

  /**
   * Returns the basis function vectors directly.
   * The layout is: [valuesDX..., valuesDY..., valuesDZ...]; i.e. all values for one derivative
   * direction before another starts.
   */
  [[nodiscard]] const std::vector<T>& data() const { return data_; }
};

//==============================================================================

template <class T>
class TimeBasisFunctionGenerator {
  private:
  T tau_;

  static T sampleJacobiPolynomial(T x, std::uint32_t n) {
    return seissol::functions::JacobiP(n, 0, 0, x);
  }

  public:
  explicit TimeBasisFunctionGenerator(T tau) : tau_(tau) {}

  T operator()(std::uint32_t i) const { return functions::DubinerP<1>({i}, {tau_}); }
};

template <class T>
struct SampledTimeBasisFunctions {
  private:
  static_assert(std::is_arithmetic_v<T>,
                "Type T for SampledTimeBasisFunctions must be arithmetic.");

  std::vector<T> data_;

  public:
  SampledTimeBasisFunctions(std::uint32_t order, T tau) : data_(order) {
    TimeBasisFunctionGenerator<T> gen(tau);

    for (std::uint32_t ord = 0; ord < order; ord++) {
      data_[ord] = gen(ord);
    }
  }

  template <class ConstIterator>
  T evalWithCoeffs(ConstIterator coeffIter) const {
    return std::inner_product(data_.begin(), data_.end(), coeffIter, static_cast<T>(0));
  }

  [[nodiscard]] auto size() const { return data_.size(); }

  [[nodiscard]] const std::vector<T>& data() const { return data_; }
};

namespace tri_dubiner {
inline void evaluatePolynomials(double* phis, double xi, double eta, std::int32_t numPoly) {
  assert(numPoly > 0);
  std::uint32_t idx = 0;
  for (std::uint32_t d = 0; d <= static_cast<std::uint32_t>(numPoly); ++d) {
    for (std::uint32_t j = 0; j <= d; ++j) {
      phis[idx++] = seissol::functions::TriDubinerP({d - j, j}, {xi, eta});
    }
  }
}

inline void evaluateGradPolynomials(double* phis, double xi, double eta, std::int32_t numPoly) {
  assert(numPoly > 0);
  std::uint32_t idx = 0;
  for (std::uint32_t d = 0; d <= static_cast<std::uint32_t>(numPoly); ++d) {
    for (std::uint32_t j = 0; j <= d; ++j) {
      const auto grad = seissol::functions::gradTriDubinerP({d - j, j}, {xi, eta});
      for (const auto& g : grad) {
        phis[idx++] = g;
      }
    }
  }
}
} // namespace tri_dubiner
} // namespace seissol::basisFunction

#endif // SEISSOL_SRC_NUMERICAL_BASISFUNCTION_H_
