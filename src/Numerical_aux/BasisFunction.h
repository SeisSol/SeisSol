/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2015, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 */

#ifndef BASIS_FUNCTION_H
#define BASIS_FUNCTION_H

#include <cmath>
#include <numeric>
#include <vector>
#include <type_traits>

#include "Functions.h"

namespace seissol {
namespace basisFunction {

//------------------------------------------------------------------------------

/**
 * Functor to sample Basis functions for a given point in the reference
 * tetrahedron.
 * @param T the type to calculate internally.
 */
template<class T>
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
    BasisFunctionGenerator(T xi, T eta, T zeta) :
      xi_(xi), eta_(eta), zeta_(zeta)
    {}

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
      return functions::TetraDubinerP({i, j, k}, {xi_, eta_, zeta_});
    }

};
//------------------------------------------------------------------------------

/**
 * Functor to sample Basis function derivatives for a given point in the reference
 * tetrahedron.
 * @param T the type to calculate internally.
 */
template<class T>
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
  BasisFunctionDerivativeGenerator(T xi, T eta, T zeta) :
      xi_(xi), eta_(eta), zeta_(zeta)
  {}

  /**
   * Allows functor style call.
   * @param i Polynomial index information
   * @param j Polynomial index information
   * @param k Polynomial index information
   */
  std::array<T, 3> operator()(unsigned int i, unsigned int j, unsigned int k) const {
    return functions::gradTetraDubinerP({i, j, k}, {xi_, eta_, zeta_});
  }

};

inline unsigned int basisFunctionsForOrder(unsigned int order)
{
  return (order)*(order+1)*(order+2)/6;
}

/**
 * This class represents a vector Basis functions sampled at a specific point.
 * @param T denotes the type to calculate internally.
 */
template<class T>
class SampledBasisFunctions {
  static_assert(std::is_arithmetic<T>::value, "Type T for SampledBasisFunctions must be arithmetic.");

public:
  /** The basis function samples */
  std::vector<T> m_data;

public:
  /**
   * Constructor to generate the sampled basis functions of given order
   * and at a given point in the reference tetrahedron.
   * @param order The order of the computation. It determines how many Basis
   * functions are generated. @see getBasisFunctionsPerOrder
   * @param eta The eta coordinate in the reference tetrahedron.
   * @param zeta The zeta coordinate in the reference tetrahedron.
   * @param xi The xi coordinate in the reference tetrahedron.
   */
  SampledBasisFunctions(unsigned int order,
                        T xi, T eta, T zeta)
      : m_data(basisFunctionsForOrder(order))
  {
    BasisFunctionGenerator<T> gen(xi, eta, zeta);

    unsigned int i = 0;
    for (unsigned int ord = 0; ord < order; ord++)
      for (unsigned int k = 0; k <= ord; k++)
        for (unsigned int j = 0; j <= ord-k; j++)
          m_data[i++] = gen(ord-j-k, j, k);
  }

public:
  /**
   * Function to evaluate the samples by multiplying the sampled Basis
   * function with its coefficient and summing up the products.
   * res = c0 * bf0 + c1 * bf1 + ... + cn * bfn
   * @param coeffIter the const iterator to read the coefficients
   */
  template<class ConstIterator>
  T evalWithCoeffs(ConstIterator coeffIter) const
  {
    return std::inner_product(m_data.begin(), m_data.end(), coeffIter, static_cast<T>(0));
  }

  /**
   * Returns the amount of Basis functions this class represents.
   */
  unsigned int getSize() const
  {
    return m_data.size();
  }
};

//------------------------------------------------------------------------------

/**
 * This class represents the derivatives of vector Basis functions sampled at a specific point.
 * @param T denotes the type to calculate internally.
 */
template<class T>
class SampledBasisFunctionDerivatives {
  static_assert(std::is_arithmetic<T>::value, "Type T for SampledBasisFunctions must be arithmetic.");
  
public:
  /** The basis function derivative samples w.r.t. xi, eta, zeta */
  std::array<std::vector<T>, 3> m_data;

public:
    /**
     * Constructor to generate the sampled basis functions of given order
     * and at a given point in the reference tetrahedron.
     * @param order The order of the computation. It determines how many Basis
     * functions are generated. @see getBasisFunctionsPerOrder
     * @param eta The eta coordinate in the reference tetrahedron.
     * @param zeta The zeta coordinate in the reference tetrahedron.
     * @param xi The xi coordinate in the reference tetrahedron.
     */
    SampledBasisFunctionDerivatives(unsigned int order,
                                    T xi, T eta, T zeta)
    	: m_data({std::vector<T>(basisFunctionsForOrder(order)), std::vector<T>(basisFunctionsForOrder(order)), std::vector<T>(basisFunctionsForOrder(order))})
    {
        BasisFunctionDerivativeGenerator<T> gen(xi, eta, zeta);

        unsigned int i = 0;
        for (unsigned int ord = 0; ord < order; ord++) {
          for (unsigned int k = 0; k <= ord; k++) {
            for (unsigned int j = 0; j <= ord - k; j++) {
              const auto derivatives = gen(ord - j - k, j, k);
              for (unsigned int direction = 0; direction < 3; direction++) {
                m_data[direction][i] = derivatives[direction];
              }
              i++;
            }
          }
        }
    }

public:
    /**
     * Function to evaluate the samples by multiplying the sampled Basis
     * function with its coefficient and summing up the products.
     * res = c0 * bf0 + c1 * bf1 + ... + cn * bfn
     * @param coeffIter the const iterator to read the coefficients
     */
    template<class ConstIterator>
    std::array<T, 3> evalWithCoeffs(ConstIterator coeffIter) const
    {
      return {std::inner_product(m_data[0].begin(), m_data[0].end(), coeffIter, static_cast<T>(0)),
              std::inner_product(m_data[1].begin(), m_data[1].end(), coeffIter, static_cast<T>(0)),
              std::inner_product(m_data[2].begin(), m_data[2].end(), coeffIter, static_cast<T>(0))};
    }

    /**
     * Returns the amount of Basis functions this class represents.
     */
    unsigned int getSize() const
    {
        return m_data[0].size();
    }
};

//==============================================================================

template<class T>
class TimeBasisFunctionGenerator {
private:
    T tau_;

    static T sampleJacobiPolynomial(T x, unsigned int n) {
      return seissol::functions::JacobiP(n, 0, 0, x);
    }

public:
    TimeBasisFunctionGenerator(T tau) :
      tau_(tau)
    {}

    T operator()(unsigned int i) const {
      return functions::DubinerP<1>({i}, {tau_});
    }

};

template<class T>
class SampledTimeBasisFunctions {
  static_assert(std::is_arithmetic<T>::value, "Type T for SampledTimeBasisFunctions must be arithmetic.");
  
public:
    std::vector<T> m_data;

public:
    SampledTimeBasisFunctions(unsigned int order, T tau)
    	: m_data(order)
    {
        TimeBasisFunctionGenerator<T> gen(tau);

        for (unsigned int ord = 0; ord < order; ord++) {
          m_data[ord] = gen(ord);
        }
    }

    template<class ConstIterator>
    T evalWithCoeffs(ConstIterator coeffIter) const
    {
        return std::inner_product(m_data.begin(), m_data.end(), coeffIter, static_cast<T>(0));
    }

    unsigned int getSize() const
    {
        return m_data.size();
    }
};

} // namespace basisFunction
} //namespace seissol

#endif // BASIS_FUNCTION_H
