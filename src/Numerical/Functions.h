// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_NUMERICAL_FUNCTIONS_H_
#define SEISSOL_SRC_NUMERICAL_FUNCTIONS_H_

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>

#include "Common/Marker.h"

namespace seissol::functions {

/**
 * @brief host standard math functions used in template metaprogramming
 */
struct HostStdFunctions {
#pragma omp declare simd
  template <typename T>
  SEISSOL_HOSTDEVICE static T exp(T value) {
    return std::exp(value);
  }
#pragma omp declare simd
  template <typename T1, typename... T>
  SEISSOL_HOSTDEVICE static T1 max(T1 value1, T... value) {
    return std::max(value1, value...);
  }
#pragma omp declare simd
  template <typename T1, typename... T>
  SEISSOL_HOSTDEVICE static T1 min(T1 value1, T... value) {
    return std::min(value1, value...);
  }
#pragma omp declare simd
  template <typename T>
  SEISSOL_HOSTDEVICE static T ceil(T value) {
    return std::ceil(value);
  }
#pragma omp declare simd
  template <typename T>
  SEISSOL_HOSTDEVICE static T floor(T value) {
    return std::floor(value);
  }
};

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
double SingularityFreeJacobiPRecursion(
    double x, double y, const std::array<double, 5>& cm, double pm1, double pm2);

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
std::array<double, 3>
    SingularityFreeJacobiPAndDerivatives(unsigned n, unsigned a, unsigned b, double x, double y);

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
double TriDubinerP(const std::array<unsigned, 2>& i, const std::array<double, 2>& xi);

/**
 * @brief Gradient of Dubiner basis on triangle
 *
 * See TriDubinerP.
 *
 * @return Gradient at xi
 */
std::array<double, 2> gradTriDubinerP(const std::array<unsigned, 2>& i,
                                      const std::array<double, 2>& xi);

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
double TetraDubinerP(const std::array<unsigned, 3>& i, const std::array<double, 3>& xi);

/**
 * @brief Gradient of Dubiner basis on tetrahedron
 *
 * See TetraDubinerP.
 *
 * @return Gradient at xi
 */
std::array<double, 3> gradTetraDubinerP(const std::array<unsigned, 3>& i,
                                        const std::array<double, 3>& xi);

/**
 * @brief Templated Dubiner basis for D=1,2,3.
 *
 * Reference element is given by vertices
 * D = 1: (0), (1)
 * D = 2: (0,0), (1,0), (0,1)
 * D = 3: (0,0,0), (1,0,0), (0,1,0), (0,0,1)
 */
template <std::size_t D>
double DubinerP(const std::array<unsigned, D>& i, const std::array<double, D>& xi);

/**
 * @brief Templated gradient for D=1,2,3.
 */
template <std::size_t D>
std::array<double, D> gradDubinerP(const std::array<unsigned, D>& i,
                                   const std::array<double, D>& xi);

} // namespace seissol::functions

#endif // SEISSOL_SRC_NUMERICAL_FUNCTIONS_H_
