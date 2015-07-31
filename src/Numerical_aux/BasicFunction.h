#pragma once

#include <cmath>
#include <vector>
#include <numeric>

namespace BasicFunction {

//------------------------------------------------------------------------------

template<class T>
inline T getBasicFunctionsPerOrder(T order) 
{
    return (order+1)*(order+2)*(order+3)/6;
}

//------------------------------------------------------------------------------

/**
 * This class represents a vector basic functions sampled at a specific point.
 * @param T denotes the type to calculate internally.
 */
template<class T>
class SampledBasicFunctions {
private:
    /** Holds the samples **/
    std::vector<T> data_;

    /**
     * Private contstructor @see sampleAt
     * @param numBasicFunctions sets the number of basic functions to sample.
     */
    SampledBasicFunctions(unsigned int numBasicFunctions) :
        data_(numBasicFunctions)
    {}

public:
    /**
     * Function to evaluate the samples by multiplying the sampled basic
     * function with its coefficient and summing up the products.
     * res = c0 * bf0 + c1 * bf1 + ... + cn * bfn
     * @param iter the const iterator to read the coefficients
     */
    template<class ConstIterator>
    T evalWithCoefs(ConstIterator iter) const
    {
        return std::inner_product(data_.begin(), data_.end(), iter, 0);
    }

    /**
     * Returns the amount of basic functions this class represents.
     */
    std::size_t getSize() const
    {
        return data_.size();
    }
    
    /**
     * Factoryfunction to generate the sampled basic functions of given order
     * and at a given point in the reference tetrahedron.
     * @param order The order of the computation. It determines how many basic
     * functions are generated. @see getBasicFunctionsPerOrder
     * @param eta The eta coordinate in the reference terahedron.
     * @param zeta The zeta coordinate in the reference terahedron.
     * @param xi The xi coordinate in the reference terahedron.
     */
    static SampledBasicFunctions sampleAt( unsigned int order,
            T eta, T zeta, T xi); 
};

//------------------------------------------------------------------------------

/**
 * Functor to sample basic functions for a given point in the reference
 * tetrahedron.
 * @param T the type to calculate internally.
 */
template<class T>
class BasicFunctionGenerator {
private:
    T eta_, zeta_, xi_;

    /**
     * Function to sample a Jacobi Polynomial.
     * Source : https://github.com/pjabardo/Jacobi.jl/blob/master/src/jac_poly.jl
     * @param x Sampling point
     */
    static T sampleJacobiPolynomial(T x, unsigned int n, T a, T b);

public:
    /**
     * Constructs a BasicFunctionGenerator which fixes it to a specific point
     * in the reference tetrahedron.
     * @param eta The eta coordinate in the reference terahedron.
     * @param zeta The zeta coordinate in the reference terahedron.
     * @param xi The xi coordinate in the reference terahedron.
     */
    BasicFunctionGenerator(T eta, T zeta, T xi) :
        eta_(eta), zeta_(zeta), xi_(xi)
    {}

    /**
     * Allows functor style call. Generates the sampled basic function given 
     * polynomial information for the generation of the basic function.
     * This algorithm uses the calculation method from the paper :
     * "Seismic Wave Simulation for Complex Rheologies on Unstructured Meshes"
     * by Josep de la Puente
     * @param i Polynomial index information
     * @param j Polynomial index information
     * @param k Polynomial index information
     */
    T operator()(unsigned int i, unsigned int j, unsigned int k) const;

};

//==============================================================================

template<class T>
T BasicFunctionGenerator<T>::sampleJacobiPolynomial(T x, unsigned int n, T a, T b)
{
    if (n==0)
        return 1.0;
    if (n==1)
        return 0.5 * (a - b + (a + b + 2.0)*x);

    T p0 = 1.0;
    T p1 = 0.5 * (a - b + (a + b + 2.0)*x);
    T p2 = 0.0;

    for (T i = 1; i < n; i++) {
        T a1 = 2.0*(i+1.0)*(i+a+b+1.0)*(2.0*i+a+b);
        T a2 = (2.0*i+a+b+1.0)*(a*a-b*b);
        T a3 = (2.0*i+a+b)*(2.0*i+a+b+1.0)*(2.0*i+a+b+2.0);
        T a4 = 2.0*(i+a)*(i+b)*(2.0*i+a+b+2.0);
        p2 = 1.0/a1*((a2 + a3*x)*p1 - a4*p0);
        p0 = p1;
        p1 = p2;
    }

    return p2;
}

//------------------------------------------------------------------------------

template<class T>
T BasicFunctionGenerator<T>::operator()(unsigned int i, unsigned int j,
        unsigned int k) const
{
    const T r = (eta_-1.0+zeta_+2.0*xi_)/(1.0-eta_-zeta_);
    const T s = (2.0*eta_-1.0+zeta_)/(1.0-zeta_);
    const T t = 2.0*zeta_-1.0;

    const T theta_a = sampleJacobiPolynomial(r, i, 0, 0);
    const T theta_b = std::pow((1-s)/2, i) * sampleJacobiPolynomial(s, j, 2*i+1, 0);
    const T theta_c = std::pow((1-t)/2, i+j) * sampleJacobiPolynomial(t, k, 2*i+2*j+2, 0);

    return  theta_a * theta_b * theta_c;
}

//==============================================================================

template <class T>
SampledBasicFunctions<T> SampledBasicFunctions<T>::sampleAt(
        unsigned int order,
        T eta, T zeta, T xi)
{
    SampledBasicFunctions sbf(getBasicFunctionsPerOrder(order));

    BasicFunctionGenerator<T> gen(eta, zeta, xi);

    typename std::vector<T>::iterator outiter = sbf.data_.begin();

    for (unsigned int ord = 0; ord <= order; ord++)
        for (unsigned int k = 0; k <= ord; k++)
            for (unsigned int j = 0; j <= ord; j++)
                for (unsigned int i = 0; i <= ord; i++)
                    if (i+j+k == ord)
                        *(outiter++) = gen(i, j, k);
    return sbf;
}

//------------------------------------------------------------------------------

} // namespace BasicFunction
