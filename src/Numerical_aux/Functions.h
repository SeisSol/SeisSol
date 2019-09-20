/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
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
 *
 * @section DESCRIPTION
 **/

//  TODO: Merge this file with BasisFunction.h
#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include <cmath>

namespace seissol {
  namespace functions {
    /** Computes n!
     */
    inline unsigned factorial(unsigned n)
    {
      unsigned f = 1;
      while (n > 0) {
        f *= n;
        --n;
      }
      return f;
    }
    
    /** Calculate Jacobi polynomial recursively at x.
     *  See recurrence relation on https://en.wikipedia.org/wiki/Jacobi_polynomials
     */
    inline double JacobiP(unsigned n, unsigned a, unsigned b, double x) {
      if (n == 0) {
        return 1.0;
      }
      double Pm_2;
      double Pm_1 = 1.0;
      double Pm = 0.5*a - 0.5*b + (1.0 + 0.5*(a+b)) * x;
      double a2_b2 = static_cast<double>(a*a)-static_cast<double>(b*b);
      for (unsigned m = 2; m <= n; ++m) {
        Pm_2 = Pm_1;
        Pm_1 = Pm;
        Pm = ( (2.0*m+a+b-1.0)*(a2_b2 + (2.0*m+a+b)*(2.0*m+a+b-2.0)*x)*Pm_1 - 2.0*(m+a-1.0)*(m+b-1.0)*(2.0*m+a+b)*Pm_2 ) / (2.0*m*(m+a+b)*(2.0*m+a+b-2.0));
      }      
      return Pm;
    }
    
    /** Computes JacobiP(n, a, b, x/y) * y^n.
     *  Works for y = 0.
     */
    inline double SingularityFreeJacobiP(unsigned n, unsigned a, unsigned b, double x, double y) {
      if (n == 0) {
        return 1.0;
      }
      double Pm_2;
      double Pm_1 = 1.0;
      double Pm = (0.5*a - 0.5*b) * y + (1.0 + 0.5*(a+b)) * x;
      double a2_b2 = static_cast<double>(a*a)-static_cast<double>(b*b);
      for (unsigned m = 2; m <= n; ++m) {
        Pm_2 = Pm_1;
        Pm_1 = Pm;
        Pm = ( (2.0*m+a+b-1.0)*(a2_b2*y + (2.0*m+a+b)*(2.0*m+a+b-2.0)*x)*Pm_1 - 2.0*y*y*(m+a-1.0)*(m+b-1.0)*(2.0*m+a+b)*Pm_2 ) / (2.0*m*(m+a+b)*(2.0*m+a+b-2.0));
      }
      return Pm;
    }

     /** Calculate first derivative of Jacobi polynomial at x.
     */
    inline double JacobiPFirstDerivative(unsigned n, unsigned a, unsigned b, double x) {
      if (n == 0) {
        return 0.0;
      }
      return 0.5*(n+a+b+1.0)*JacobiP(n-1, a+1, b+1, x);
    }
    
    /** Evaluate Dubiner basis on triangle
     *
     * Singularity-free variant inspired by
     * R. C. Kirby, "Singularity-free evaluation of collapsed-coordinate orthogonal polynomials",
     * ACM TOMS 37.1, Article 5, doi: 10.1145/1644001.1644006
     */
    inline double TetraDubinerP(unsigned i, unsigned j, unsigned k, double xi, double eta, double zeta) {
      double r_num = 2.0 * xi - 1.0 + eta + zeta;
      double s_num = 2.0 * eta - 1.0 + zeta;
      double t = 2.0 * zeta - 1.0;
      double sigmatheta = 1.0 - eta - zeta;
      double theta = 1.0 - zeta;

      double ti = SingularityFreeJacobiP(i, 0, 0, r_num, sigmatheta);
      double tij = SingularityFreeJacobiP(j, 2*i+1, 0, s_num, theta);
      double tijk = SingularityFreeJacobiP(k, 2*i+2*j+2, 0, t, 1.0);

      return ti * tij * tijk;
      // The following variant does not work for some points, that is, when the denominator is zero.
      /*double r = 2.0 * xi / (1.0 - eta - zeta) - 1.0;
      double s = 2.0 * eta / (1.0 - zeta) - 1.0;
      double t = 2.0 * zeta - 1.0;
      
      double ti = JacobiP(i, 0, 0, r);
      double tij = JacobiP(j, 2*i+1, 0, s) * std::pow(0.5*(1.0-s), i);
      double tijk = JacobiP(k, 2*i+2*j+2, 0, t) * std::pow(0.5*(1.0-t), i+j);
      
      return ti * tij * tijk;*/
    }
  }
}

#endif
