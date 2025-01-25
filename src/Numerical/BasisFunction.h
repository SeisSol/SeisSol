// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_NUMERICAL_BASISFUNCTION_H_
#define SEISSOL_SRC_NUMERICAL_BASISFUNCTION_H_

#include "generated_code/init.h"
#include <cmath>
#include <numeric>
#include <type_traits>
#include <vector>

#include "Functions.h"
#include "Transformation.h"

namespace seissol {
namespace basisFunction {

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
  static T sampleJacobiPolynomial(T x, unsigned int n, unsigned int a, unsigned int b) {
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
  T operator()(unsigned int i, unsigned int j, unsigned int k) const {
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
  std::array<T, 3> operator()(unsigned int i, unsigned int j, unsigned int k) const {
    std::array<double, 3> gradEvaluated =
        functions::gradTetraDubinerP({i, j, k}, {xi_, eta_, zeta_});
    return {static_cast<T>(gradEvaluated[0]),
            static_cast<T>(gradEvaluated[1]),
            static_cast<T>(gradEvaluated[2])};
  }
};

inline unsigned int basisFunctionsForOrder(unsigned int order) {
  return (order) * (order + 1) * (order + 2) / 6;
}

/**
 * This class represents a vector Basis functions sampled at a specific point.
 * @param T denotes the type to calculate internally.
 */
template <class T>
class SampledBasisFunctions {
  static_assert(std::is_arithmetic<T>::value,
                "Type T for SampledBasisFunctions must be arithmetic.");

  public:
  /** The basis function samples */
  std::vector<T> m_data{};

  public:
  SampledBasisFunctions() {};
  /**
   * Constructor to generate the sampled basis functions of given order
   * and at a given point in the reference tetrahedron.
   * @param order The order of the computation. It determines how many Basis
   * functions are generated. @see getBasisFunctionsPerOrder
   * @param eta The eta coordinate in the reference tetrahedron.
   * @param zeta The zeta coordinate in the reference tetrahedron.
   * @param xi The xi coordinate in the reference tetrahedron.
   */
  SampledBasisFunctions(unsigned int order, T xi, T eta, T zeta)
      : m_data(basisFunctionsForOrder(order)) {
    BasisFunctionGenerator<T> gen(xi, eta, zeta);

    unsigned int i = 0;
    for (unsigned int ord = 0; ord < order; ord++)
      for (unsigned int k = 0; k <= ord; k++)
        for (unsigned int j = 0; j <= ord - k; j++)
          m_data[i++] = gen(ord - j - k, j, k);
  }

  public:
  /**
   * Function to evaluate the samples by multiplying the sampled Basis
   * function with its coefficient and summing up the products.
   * res = c0 * bf0 + c1 * bf1 + ... + cn * bfn
   * @param coeffIter the const iterator to read the coefficients
   */
  template <class ConstIterator>
  T evalWithCoeffs(ConstIterator coeffIter) const {
    return std::inner_product(m_data.begin(), m_data.end(), coeffIter, static_cast<T>(0));
  }

  /**
   * Returns the amount of Basis functions this class represents.
   */
  unsigned int getSize() const { return m_data.size(); }
};

//------------------------------------------------------------------------------

/**
 * This class represents the derivatives of vector Basis functions sampled at a specific point.
 * @param T denotes the type to calculate internally.
 */
template <class T>
class SampledBasisFunctionDerivatives {
  static_assert(std::is_arithmetic<T>::value,
                "Type T for SampledBasisFunctions must be arithmetic.");

  public:
  /**
   * The basis function derivative samples w.r.t. the three spatial dimension
   * Use DenseTensorView to access data
   */
  std::vector<T> m_data{};

  public:
  SampledBasisFunctionDerivatives() {};
  /**
   * Constructor to generate the sampled basis functions of given order
   * and at a given point in the reference tetrahedron.
   * @param order The order of the computation. It determines how many Basis
   * functions are generated. @see getBasisFunctionsPerOrder
   * @param eta The eta coordinate in the reference tetrahedron.
   * @param zeta The zeta coordinate in the reference tetrahedron.
   * @param xi The xi coordinate in the reference tetrahedron.
   */
  SampledBasisFunctionDerivatives(unsigned int order, T xi, T eta, T zeta)
      : m_data(3 * basisFunctionsForOrder(order)) {
    BasisFunctionDerivativeGenerator<T> gen(xi, eta, zeta);
    auto dataView = init::basisFunctionDerivativesAtPoint::view::create(m_data.data());

    unsigned int i = 0;
    for (unsigned int ord = 0; ord < order; ord++) {
      for (unsigned int k = 0; k <= ord; k++) {
        for (unsigned int j = 0; j <= ord - k; j++) {
          const auto derivatives = gen(ord - j - k, j, k);
          for (unsigned int direction = 0; direction < 3; direction++) {
            dataView(i, direction) = derivatives[direction];
          }
          i++;
        }
      }
    }
  }

  /**
   * After a call to the constructor m_data contains the sampled derivatives
   * w.r.t. xi, eta, zeta on the reference triangle. Use this function to
   * transform the derivatives to derivatives w.r.t. to x, y, z in a physical
   * tetrahedron.
   * @param coords coords[i] contains the 3 coordinates of the ith vertex of the
   * physical tetrahedron.
   */
  void transformToGlobalCoordinates(const double* coords[4]) {
    real xCoords[4];
    real yCoords[4];
    real zCoords[4];
    for (size_t i = 0; i < 4; ++i) {
      xCoords[i] = coords[i][0];
      yCoords[i] = coords[i][1];
      zCoords[i] = coords[i][2];
    }

    real gradXi[3];
    real gradEta[3];
    real gradZeta[3];

    seissol::transformations::tetrahedronGlobalToReferenceJacobian(
        xCoords, yCoords, zCoords, gradXi, gradEta, gradZeta);
    std::vector<T> oldData = m_data;

    auto oldView = init::basisFunctionDerivativesAtPoint::view::create(oldData.data());
    auto newView = init::basisFunctionDerivativesAtPoint::view::create(m_data.data());
    for (size_t i = 0; i < init::basisFunctionDerivativesAtPoint::Shape[0]; ++i) {
      for (size_t direction = 0; direction < init::basisFunctionDerivativesAtPoint::Shape[1];
           ++direction) {
        // dpsi / di = dphi / dxi * dxi / di + dphi / deta * deta / di + dphi / dzeta * dzeta / di
        newView(i, direction) = oldView(i, 0) * gradXi[direction] +
                                oldView(i, 1) * gradEta[direction] +
                                oldView(i, 2) * gradZeta[direction];
      }
    }
  }

  /**
   * Returns the amount of Basis functions this class represents.
   */
  unsigned int getSize() const { return m_data.size(); }
};

//==============================================================================

template <class T>
class TimeBasisFunctionGenerator {
  private:
  T tau_;

  static T sampleJacobiPolynomial(T x, unsigned int n) {
    return seissol::functions::JacobiP(n, 0, 0, x);
  }

  public:
  TimeBasisFunctionGenerator(T tau) : tau_(tau) {}

  T operator()(unsigned int i) const { return functions::DubinerP<1>({i}, {tau_}); }
};

template <class T>
class SampledTimeBasisFunctions {
  static_assert(std::is_arithmetic<T>::value,
                "Type T for SampledTimeBasisFunctions must be arithmetic.");

  public:
  std::vector<T> m_data;

  public:
  SampledTimeBasisFunctions(unsigned int order, T tau) : m_data(order) {
    TimeBasisFunctionGenerator<T> gen(tau);

    for (unsigned int ord = 0; ord < order; ord++) {
      m_data[ord] = gen(ord);
    }
  }

  template <class ConstIterator>
  T evalWithCoeffs(ConstIterator coeffIter) const {
    return std::inner_product(m_data.begin(), m_data.end(), coeffIter, static_cast<T>(0));
  }

  unsigned int getSize() const { return m_data.size(); }
};

namespace tri_dubiner {
inline void evaluatePolynomials(double* phis, double xi, double eta, int numPoly) {
  assert(numPoly > 0);
  unsigned idx = 0;
  for (unsigned int d = 0; d <= static_cast<unsigned>(numPoly); ++d) {
    for (unsigned int j = 0; j <= d; ++j) {
      phis[idx++] = seissol::functions::TriDubinerP({d - j, j}, {xi, eta});
    }
  }
}

inline void evaluateGradPolynomials(double* phis, double xi, double eta, int numPoly) {
  assert(numPoly > 0);
  unsigned idx = 0;
  for (unsigned int d = 0; d <= static_cast<unsigned>(numPoly); ++d) {
    for (unsigned int j = 0; j <= d; ++j) {
      const auto grad = seissol::functions::gradTriDubinerP({d - j, j}, {xi, eta});
      for (const auto& g : grad) {
        phis[idx++] = g;
      }
    }
  }
}
} // namespace tri_dubiner
} // namespace basisFunction
} // namespace seissol

#endif // SEISSOL_SRC_NUMERICAL_BASISFUNCTION_H_
