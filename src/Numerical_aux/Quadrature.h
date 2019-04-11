/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2015-2017, SeisSol Group
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
 **/

#ifndef QUADRATURE_H_
#define QUADRATURE_H_

#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>
#include <vector>

#include "utils/logger.h"
#include "Numerical_aux/Functions.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace seissol {
  namespace quadrature {
    static unsigned const MaxIterations = 100;
    static double const Tolerance = 10.*std::numeric_limits<double>::epsilon();

    /** Returns quadrature points for the interval [-1,1], i.e.
     *  int_{-1}^{1} f(y)dy = sum_{i=0}^{n-1} f(points[i]) * weights[i]
     *  For other intervals use 
     *  int_{a}^{b} f(y)dy = (b-a)/2. * sum_{i=0}^{n-1} f( ((b-a) * points[i] + a + b) / 2.) * weights[i]
     */
    inline void GaussLegendre(double* points, double* weights, unsigned n)
    {
      // The polynomials are symmetric, thus we only need to find the first half.
      unsigned nh = (n+1)/2;
      for (unsigned i = 1; i <= nh; ++i) {
        // x = Initial guess for polynomial root
        double x = cos(M_PI * (4.*i - 1.) / (4.*n+2.)), w;
        double Pn = 0.0, dPn = 1.0, Pn_2, Pn_1;
        unsigned it = 0;
        // Refine polynomial roots with Newton iteration
        do {
          x -= Pn/dPn;
          Pn_1 = 0.0;
          Pn = 1.0;
          // Recursive procedure to calculate the n-th Legendre polynomial at x
          for (unsigned j = 1; j <= n; ++j) {
            Pn_2 = Pn_1;
            Pn_1 = Pn;
            Pn = ((2.*j-1.) * x * Pn_1 - (j-1.) * Pn_2) / j;
          }
          // Derivative at x
          dPn = (n*Pn_1 - n*x*Pn) / (1.-x*x);
        } while (fabs(Pn) > Tolerance && ++it < MaxIterations);
        // Weight = 2 / [(1-x^2) * P_n'(x)^2]
        w = 2./((1-x*x) * dPn * dPn);
        points[i-1] = -x;
        points[n-i] = x;
        weights[i-1] = w;
        weights[n-i] = w;
      }
    }
    
    /** Returns quadrature points for the interval [-1,1] with weight function (1-x)^a * (1+x)^b, i.e.
     *  int_{-1}^{1} f(y)dy = sum_{i=0}^{n-1} f(points[i]) * weights[i]
     *  For other intervals use 
     *  int_{a}^{b} f(y)dy = (b-a)/2. * sum_{i=0}^{n-1} f( ((b-a) * points[i] + a + b) / 2.) * weights[i]
     * 
     *  Note: Initial guess ported from Fortran gauss_jacobi routine.
     */     
    inline void GaussJacobi(double* points, double* weights, unsigned n, unsigned a, unsigned b)
    {
      double weightFactor = -(2.0*n+a+b+2) * functions::factorial(n+a) * functions::factorial(n+b) * (1 << (a+b)) /
                             ((n+a+b+1.0) * functions::factorial(n+a+b) * functions::factorial(n+1));
      for (unsigned i = 1; i <= n; ++i) {
        // x = Initial guess for polynomial root
        double x = cos( M_PI * (0.5*a + i - 0.25) / (0.5*(1.0 + a + b) + n) );
        double Pn = 0.0, dPn = 1.0;
        unsigned it = 0;
        // Refine polynomial roots with Newton iteration
        do {
          x -= Pn/dPn;
          Pn = functions::JacobiP(n, a, b, x);
          dPn = functions::JacobiPFirstDerivative(n, a, b, x);
        } while (fabs(Pn) > Tolerance && ++it < MaxIterations);
        points[i-1] = x;
        weights[i-1] = weightFactor / (functions::JacobiP(n+1, a, b, x) * dPn);
      }
    }
    
    /** Returns quadrature points for the reference triangle with vertices (0,0),(1,0),(0,1), i.e.
     *  int_{0}^{1} int_{0}^{1-y} f(x,y)dxdy = sum_{i=0}^{n^2} f(points[i][0], points[i][1]) * weights[i]
     *  
     *  n is the polynomial degree. Make sure that points and weights have space for n^2 entries.
     */     
    inline void TriangleQuadrature(double (*points)[2], double* weights, unsigned n)
    {
      double* points0  = new double[n];
      double* weights0 = new double[n];
      double* points1  = new double[n];
      double* weights1 = new double[n];
      
      GaussJacobi(points0, weights0, n, 0, 0);
      GaussJacobi(points1, weights1, n, 1, 0);
      
      for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
          unsigned idx = i*n + j;
          points[idx][0] = 0.5 * (1.0 + points1[i]);
          points[idx][1] = 0.25 * (1.0 + points0[j]) * (1.0 - points1[i]);
          weights[idx] = weights1[i] * weights0[j] * 0.125;
        }        
      }
            
      delete[] points0; delete[] weights0;
      delete[] points1; delete[] weights1;
    }

    /** Quadrature formula of arbitrary accuracy on the reference tetrahedron
     *  consisting of the nodes (0,0,0), (1,0,0), (0,1,0), (0,0,1)
     */
    inline void TetrahedronQuadrature(double (*points)[3], double* weights, unsigned int n) {
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
	    points[curIndex][2] = points2[k] * (1 - points1[j]) *
	      (1 - points0[i]);
	    weights[curIndex] = weights0[i] * weights1[j] * weights2[k];
          }
        }
      }

      // TODO(Lukas) Only when debugging?
      const double tol = 1e-6;
      double sumWeights = 0.0;
      for (size_t i = 0; i < n*n*n; ++i) {
	sumWeights += weights[i];
      }
      if (std::abs(sumWeights - 1./6.) > tol) {
	logError() << "Sum of tetrahedron quadrature weights are " << sumWeights << " /= " << 1./6.;
      }
    }
  }
}

#endif
