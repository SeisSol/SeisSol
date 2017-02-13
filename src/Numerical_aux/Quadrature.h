/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
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
 **/

#ifndef QUADRATURE_H_
#define QUADRATURE_H_

#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace seissol {
  namespace quadrature {
    static unsigned const MaxIterations = 100;

    /** Returns quadrature points for the interval [-1,1], i.e.
     *  int_{-1}^{1} f(y)dy = sum_{i=0}^{n-1} f(points[i]) * weights[i]
     *  For other intervals use 
     *  int_{a}^{b} f(y)dy = (b-a)/2. * sum_{i=0}^{n-1} f( ((b-a) * points[i] + a + b) / 2.) * weights[i]
     */
    void GaussLegendre(double* points, double* weights, unsigned n)
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
        } while (fabs(Pn) > 10.*std::numeric_limits<double>::epsilon() && ++it < MaxIterations);
        // Weight = 2 / [(1-x^2) * P_n'(x)^2]
        w = 2./((1-x*x) * dPn * dPn);
        points[i-1] = -x;
        points[n-i] = x;
        weights[i-1] = w;
        weights[n-i] = w;
      }
    }
  }
}

#endif
