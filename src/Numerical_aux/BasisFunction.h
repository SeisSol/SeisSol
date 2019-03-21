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
    static T sampleJacobiPolynomial(T x, unsigned int n, unsigned int a, unsigned int b);

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
    T operator()(unsigned int i, unsigned int j, unsigned int k) const;

};

//------------------------------------------------------------------------------

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
     * @param eta The eta coordinate in the reference terahedron.
     * @param zeta The zeta coordinate in the reference terahedron.
     * @param xi The xi coordinate in the reference terahedron.
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
    
private:
    static inline unsigned int basisFunctionsForOrder(unsigned int order)
    {
        return (order)*(order+1)*(order+2)/6;
    }
};

//==============================================================================

template<class T>
T BasisFunctionGenerator<T>::sampleJacobiPolynomial(T x, unsigned int n,
						    unsigned int a,
						    unsigned b) {
  return seissol::functions::JacobiP(n, a, b, x);
  /*
    if (n == 0) {
        return T(1);
      }
    T Pm_2;
    T Pm_1 = T(1);
    T Pm = 0.5*a - 0.5*b + (T(1) + 0.5*(a+b)) * x;
    T a2_b2 = T(a*a)-T(b*b);
    for (unsigned m = 2; m <= n; ++m) {
      Pm_2 = Pm_1;
      Pm_1 = Pm;
      Pm = ( (T(2)*m+a+b-T(1))*(a2_b2 + (T(2)*m+a+b)*(T(2)*m+a+b-T(2))*x)*Pm_1 - T(2)*(m+a-T(1))*(m+b-T(1))*(T(2)*m+a+b)*Pm_2 ) / (T(2)*m*(m+a+b)*(T(2)*m+a+b-T(2)));
    }      
    return Pm;
  */
}


//------------------------------------------------------------------------------

template<class T>
T BasisFunctionGenerator<T>::operator()(unsigned int i,
					unsigned int j,
					unsigned int k) const {
  /*
    const T r = (eta_-1.0+zeta_+2.0*xi_)/(1.0-eta_-zeta_);
    const T s = (2.0*eta_-1.0+zeta_)/(1.0-zeta_);
    const T t = 2.0*zeta_-1.0;

    const T theta_a = sampleJacobiPolynomial(r, i, 0, 0);
    const T theta_b = std::pow((1-s)/2, i) * sampleJacobiPolynomial(s, j, 2*i+1, 0);
    const T theta_c = std::pow((1-t)/2, i+j) * sampleJacobiPolynomial(t, k, 2*i+2*j+2, 0);

    return  theta_a * theta_b * theta_c;
  */
  return functions::TetraDubinerP(i, j, k, xi_, eta_, zeta_);
}

//------------------------------------------------------------------------------

} // namespace basisFunction
} //namespace seissol

#endif // BASIS_FUNCTION_H
