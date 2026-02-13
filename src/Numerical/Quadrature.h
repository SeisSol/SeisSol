// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_NUMERICAL_QUADRATURE_H_
#define SEISSOL_SRC_NUMERICAL_QUADRATURE_H_

// NOLINTNEXTLINE
#define _USE_MATH_DEFINES

#include "Geometry/MeshDefinition.h"
#include "Numerical/Functions.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <utils/logger.h>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace seissol::quadrature {
static const unsigned MaxIterations = 100;
static const double Tolerance = 10. * std::numeric_limits<double>::epsilon();

/** Returns quadrature points for the interval [-1,1], i.e.
 *  int_{-1}^{1} f(y)dy = sum_{i=0}^{n-1} f(points[i]) * weights[i]
 *  For other intervals use
 *  int_{a}^{b} f(y)dy = (b-a)/2. * sum_{i=0}^{n-1} f( ((b-a) * points[i] + a + b) / 2.) *
 * weights[i]
 */
inline void GaussLegendre(double* points, double* weights, unsigned n) {
  // The polynomials are symmetric, thus we only need to find the first half.
  const unsigned nh = (n + 1) / 2;
  for (unsigned i = 1; i <= nh; ++i) {
    // x = Initial guess for polynomial root
    double x = cos(M_PI * (4. * i - 1.) / (4. * n + 2.));
    double w = 0;
    double pn = 0.0;
    double dpn = 1.0;
    double pn2 = 0.0;
    double pn1 = 0.0;
    unsigned it = 0;
    // Refine polynomial roots with Newton iteration
    do {
      x -= pn / dpn;
      pn1 = 0.0;
      pn = 1.0;
      // Recursive procedure to calculate the n-th Legendre polynomial at x
      for (unsigned j = 1; j <= n; ++j) {
        pn2 = pn1;
        pn1 = pn;
        pn = ((2. * j - 1.) * x * pn1 - (j - 1.) * pn2) / j;
      }
      // Derivative at x
      dpn = (n * pn1 - n * x * pn) / (1. - x * x);
    } while (std::abs(pn) > Tolerance && ++it < MaxIterations);
    // Weight = 2 / [(1-x^2) * P_n'(x)^2]
    w = 2. / ((1 - x * x) * dpn * dpn);
    points[i - 1] = -x;
    points[n - i] = x;
    weights[i - 1] = w;
    weights[n - i] = w;
  }
}

/**
  Return the quadrature points and weights with Gauss Legendre with n points; shifted to the
  interval [a, b].

  @return [quadrature points, quadrature weights]
*/
inline std::pair<std::vector<double>, std::vector<double>>
    ShiftedGaussLegendre(unsigned n, double a, double b) {
  std::pair<std::vector<double>, std::vector<double>> output{std::vector<double>(n),
                                                             std::vector<double>(n)};
  GaussLegendre(output.first.data(), output.second.data(), n);
  const auto ba = b - a;
  const auto ab = b + a;
  for (auto& point : output.first) {
    point = (point * ba + ab) / 2;
  }
  for (auto& weight : output.second) {
    weight = (weight * ba) / 2;
  }
  return output;
}

/** Returns quadrature points for the interval [-1,1] with weight function (1-x)^a * (1+x)^b, i.e.
 *  int_{-1}^{1} f(y)dy = sum_{i=0}^{n-1} f(points[i]) * weights[i]
 *  For other intervals use
 *  int_{a}^{b} f(y)dy = (b-a)/2. * sum_{i=0}^{n-1} f( ((b-a) * points[i] + a + b) / 2.) *
 * weights[i]
 *
 *  Note: Initial guess ported from Fortran gauss_jacobi routine.
 */
inline void GaussJacobi(double* points, double* weights, unsigned n, unsigned a, unsigned b) {
  const double weightFactor =
      -(2.0 * n + a + b + 2) * functions::factorial(n + a) * functions::factorial(n + b) *
      (1 << (a + b)) /
      ((n + a + b + 1.0) * functions::factorial(n + a + b) * functions::factorial(n + 1));
  for (unsigned i = 1; i <= n; ++i) {
    // x = Initial guess for polynomial root
    double x = cos(M_PI * (0.5 * a + i - 0.25) / (0.5 * (1.0 + a + b) + n));
    double pn = 0.0;
    double dpn = 1.0;
    unsigned it = 0;
    // Refine polynomial roots with Newton iteration
    do {
      x -= pn / dpn;
      pn = functions::JacobiP(n, a, b, x);
      dpn = functions::JacobiPDerivative(n, a, b, x);
    } while (std::abs(pn) > Tolerance && ++it < MaxIterations);
    points[i - 1] = x;
    weights[i - 1] = weightFactor / (functions::JacobiP(n + 1, a, b, x) * dpn);
  }
}

/** Returns quadrature points for the reference triangle with vertices (0,0),(1,0),(0,1), i.e.
 *  int_{0}^{1} int_{0}^{1-y} f(x,y)dxdy = sum_{i=0}^{n^2} f(points[i][0], points[i][1]) *
 * weights[i]
 *
 *  n is the polynomial degree. Make sure that points and weights have space for n^2 entries.
 */
inline void TriangleQuadrature(double (*points)[2], double* weights, std::size_t n) {
  auto points0 = std::vector<double>(n);
  auto weights0 = std::vector<double>(n);
  auto points1 = std::vector<double>(n);
  auto weights1 = std::vector<double>(n);

  GaussJacobi(points0.data(), weights0.data(), n, 0, 0);
  GaussJacobi(points1.data(), weights1.data(), n, 1, 0);

  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      const std::size_t idx = i * n + j;
      points[idx][0] = 0.5 * (1.0 + points1[i]);
      points[idx][1] = 0.25 * (1.0 + points0[j]) * (1.0 - points1[i]);
      weights[idx] = weights1[i] * weights0[j] * 0.125;
    }
  }
}

/** Quadrature formula of arbitrary accuracy on the reference tetrahedron
 *  consisting of the nodes (0,0,0), (1,0,0), (0,1,0), (0,0,1)
 */
inline void TetrahedronQuadrature(double (*points)[3], double* weights, std::size_t n) {
  // This is a port of similarly named fortran method
  // (TetrahedronQuadraturePoints) in quadpoints.f90.
  // Note:
  // Our quadrature points are defined by the conical product of 1D Gauss-Jacobi
  // formulas with M quadrature points. See Stround, p. 28ff for details.

  auto points0 = std::vector<double>(n);
  auto points1 = std::vector<double>(n);
  auto points2 = std::vector<double>(n);

  auto weights0 = std::vector<double>(n);
  auto weights1 = std::vector<double>(n);
  auto weights2 = std::vector<double>(n);

  // Get the Gauss-Jacobi positions and weights.
  GaussJacobi(points0.data(), weights0.data(), n, 2, 0);
  GaussJacobi(points1.data(), weights1.data(), n, 1, 0);
  GaussJacobi(points2.data(), weights2.data(), n, 0, 0);

  // Shift and rescale positions because Stroud
  // integrates over [0,1] and gaujac of num. recipes
  // considers [-1,1].
  for (size_t i = 0; i < n; ++i) {
    points0[i] = 0.5 * points0[i] + 0.5;
    points1[i] = 0.5 * points1[i] + 0.5;
    points2[i] = 0.5 * points2[i] + 0.5;

    weights0[i] = 0.5 * 0.5 * 0.5 * weights0[i];
    weights1[i] = 0.5 * 0.5 * weights1[i];
    weights2[i] = 0.5 * weights2[i];
  }

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      for (size_t k = 0; k < n; ++k) {
        const auto curIndex = i * n * n + j * n + k;
        points[curIndex][0] = points0[i];
        points[curIndex][1] = points1[j] * (1 - points0[i]);
        points[curIndex][2] = points2[k] * (1 - points1[j]) * (1 - points0[i]);
        weights[curIndex] = weights0[i] * weights1[j] * weights2[k];
      }
    }
  }

  // TODO(Lukas) Only when debugging?
  constexpr double Tol = 1e-6;
  double sumWeights = 0.0;
  for (size_t i = 0; i < n * n * n; ++i) {
    sumWeights += weights[i];
  }
  if (std::abs(sumWeights - 1. / 6.) > Tol) {
    logError() << "Sum of tetrahedron quadrature weights are " << sumWeights << " /= " << 1. / 6.;
  }
}

template <std::size_t Dim>
std::pair<std::vector<std::array<double, Dim>>, std::vector<double>>
    quadrature(std::size_t degree) {
  static_assert(Dim >= 1, "Dim needs to be non-negative");
  static_assert(Dim <= 3, "Higher dimensions than 3 are not implemented for the quadrature");
  std::size_t arraysize = 1;
  for (std::size_t _ = 0; _ < Dim; ++_) {
    arraysize *= degree;
  }
  std::vector<std::array<double, Dim>> points(arraysize);
  std::vector<double> weights(arraysize);
  if constexpr (Dim == 3) {
    TetrahedronQuadrature(reinterpret_cast<double (*)[3]>(points.data()), weights.data(), degree);
  }
  if constexpr (Dim == 2) {
    TriangleQuadrature(reinterpret_cast<double (*)[2]>(points.data()), weights.data(), degree);
  }
  if constexpr (Dim == 1) {
    GaussLegendre(points.data(), weights.data(), degree);

    // convert to shifted GL
    for (auto& point : points) {
      point = (point + 1) / 2;
    }
    for (auto& weight : weights) {
      weight /= 2;
    }
  }
  return {points, weights};
}

} // namespace seissol::quadrature

#endif // SEISSOL_SRC_NUMERICAL_QUADRATURE_H_
