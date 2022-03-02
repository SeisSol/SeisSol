/**
 * Copyright (c) 2017, SeisSol Group
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
 **/
#ifndef FUNCTIONS_20201026_H
#define FUNCTIONS_20201026_H

#include <array>
#include <cstddef>
#include <cstdint>

namespace seissol {
namespace functions {

/**
 * @brief Computes \prod_{i=from}^{to} i. Returns 1 if from > to.
 */
uint64_t rangeProduct(uint64_t from, uint64_t to);

/**
 * @brief Factorial operation
 *
 * @param n
 *
 * @return n!
 */
inline uint64_t factorial(uint64_t n) { return rangeProduct(1, n); }

/**
 * @brief Evaluates Jacobi polynomial P_{n}^{(a,b)}(x).
 */
double JacobiP(unsigned n, unsigned a, unsigned b, double x);

/**
 * @brief Evaluates derivative of Jacobi polynomial.
 */
double JacobiPDerivative(unsigned n, unsigned a, unsigned b, double x);

/**
 * @brief Factors in recursion formulas used by SingularityFreeJacobiP functions.
 */
std::array<double, 5> SingularityFreeJacobiPFactors(unsigned m, unsigned a, unsigned b);

/**
 * @brief Computes JacobiP(n, a, b, x/y) * y^n.
 *
 * @param Output of SingularityFreeJacobiPFactors
 * @param Pm_1 JacobiP(n-1, a, b, x/y) * y^{n-1}
 * @param Pm_2 JacobiP(n-2, a, b, x/y) * y^{n-2}
 *
 */
double SingularityFreeJacobiPRecursion(double x, double y, std::array<double, 5> const& cm,
                                       double Pm_1, double Pm_2);

/**
 * @brief Computes JacobiP(n, a, b, x/y) * y^n.
 *
 * Works for y = 0.
 */
double SingularityFreeJacobiP(unsigned n, unsigned a, unsigned b, double x, double y);

/**
 * @brief Singularity free Jacobi polynomial evaluation with derivatives.
 *
 * Computes K_{a,b}^n(x,y) = JacobiP(n, a, b, x/y) * y^n, dK_{a,b}^n/dx, and dK_{a,b}^n/dy.
 *
 * return {K, dKdx, dKdy}
 */
std::array<double, 3> SingularityFreeJacobiPAndDerivatives(unsigned n, unsigned a, unsigned b,
                                                           double x, double y);

/**
 * @brief Evaluate Dubiner basis on reference triangle
 *
 * Reference triangle is (0,0), (1,0), (0,1).
 *
 * @param i multi-index specifying the polynomial degree
 * @param point in reference triangle
 *
 * @return Function value at xi
 */
double TriDubinerP(std::array<unsigned, 2> const& i, std::array<double, 2> const& xi);

/**
 * @brief Gradient of Dubiner basis on triangle
 *
 * See TriDubinerP.
 *
 * @return Gradient at xi
 */
std::array<double, 2> gradTriDubinerP(std::array<unsigned, 2> const& i,
                                      std::array<double, 2> const& xi);

/**
 * @brief Evaluate Dubiner basis on reference tetrahedron
 *
 * Reference tetrahedron is (0,0,0), (1,0,0), (0,1,0), (0,0,1).
 *
 * Singularity-free variant inspired by
 * R. C. Kirby, "Singularity-free evaluation of collapsed-coordinate orthogonal polynomials",
 * ACM TOMS 37.1, Article 5, doi: 10.1145/1644001.1644006
 *
 * @param i multi-index specifying the polynomial degree
 * @param point in reference tetrahedron
 *
 * @return Function value at xi
 */
double TetraDubinerP(std::array<unsigned, 3> const& i, std::array<double, 3> const& xi);

/**
 * @brief Gradient of Dubiner basis on tetrahedron
 *
 * See TetraDubinerP.
 *
 * @return Gradient at xi
 */
std::array<double, 3> gradTetraDubinerP(std::array<unsigned, 3> const& i,
                                        std::array<double, 3> const& xi);

/**
 * @brief Templated Dubiner basis for D=1,2,3.
 *
 * Reference element is given by vertices
 * D = 1: (0), (1)
 * D = 2: (0,0), (1,0), (0,1)
 * D = 3: (0,0,0), (1,0,0), (0,1,0), (0,0,1)
 */
template <std::size_t D>
double DubinerP(std::array<unsigned, D> const& i, std::array<double, D> const& xi);

/**
 * @brief Templated gradient for D=1,2,3.
 */
template <std::size_t D>
std::array<double, D> gradDubinerP(std::array<unsigned, D> const& i,
                                   std::array<double, D> const& xi);

} // namespace functions
} // namespace seissol

#endif // FUNCTIONS_20201026_H
