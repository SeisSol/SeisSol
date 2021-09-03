/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Wolf (wolf.sebastian AT tum.de, https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2015 - 2020, SeisSol Group
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
 * Point source computation.
 **/

#ifndef SOURCETERM_POINTSOURCE_H_
#define SOURCETERM_POINTSOURCE_H_

#include <Initializer/typedefs.hpp>
#include "SourceTerm/typedefs.hpp"

namespace seissol {
  namespace sourceterm {
    /** The local moment tensor shall be transformed into the global coordinate system.
     * 
     * The second order tensor (matrix) can be understood as a transform
     * on a vector, e.g. p_L = T_L * q_L. (Let L = Local moment tensor, M = Moment tensor.)
     * We are looking for the transformed tensor p_M = T_M * q_M, i.e.
     * the local moment tensor rotated by strike, dip, and rake.
     * Assume x_L = R * x_M, where R is an orthogonal matrix. Then
     * p_M = R^T * p_L = R^T * T_L * R * q_M and hence
     *   T_M = R^T * T_L * R.
     * Thus, the rotation matrix R is the transformation from the global (x,y,z)
     * coordinate system to local (fault plane) coordinate system and is obtained
     * by the successive rotations  strike (s) -> dip (d) -> rake (l).
     * 
     *                   |  cos l  sin l    | | 1               |  |  cos s -sin s    |
     * R_l * R_d * R_s = | -sin l  cos l    | |    cos d -sin d |  |  sin s  cos s    |
     *                   |                1 | |    sin d  cos d |  |                1 |
     *
     **/    
    void transformMomentTensor(real const i_localMomentTensor[3][3],
                               real const i_localVelocityComponent[3],
                               real strike,
                               real dip,
                               real rake,
                               real o_forceComponents[seissol::sourceterm::PointSources::TensorSize]);

    /** Converts equally spaced time samples to a one-dimensional
     *  piecewise linear function.
     */
    template<typename real_from>
    void samplesToPiecewiseLinearFunction1D(real_from const* i_samples,
                                            unsigned i_numberOfSamples,
                                            real i_onsetTime,
                                            real i_samplingInterval,
                                            PiecewiseLinearFunction1D* o_pwLF)
    {
      if (i_numberOfSamples == 0) {
        o_pwLF->numberOfPieces = 0;
        o_pwLF->slopes = NULL;
        o_pwLF->intercepts = NULL;
        return;        
      }

      unsigned l_np = i_numberOfSamples - 1;
      
      o_pwLF->slopes = new real[l_np];
      o_pwLF->intercepts = new real[l_np];  
      o_pwLF->onsetTime = i_onsetTime;
      o_pwLF->numberOfPieces = l_np;
      o_pwLF->samplingInterval = i_samplingInterval;
      
      
      /* The piecewise linear function shall be f(t) = m_j * t + n_j,
       * where I_j is the half-open interval [t_o + j*dt, t_o + j*dt).
       * For the j-th sample (say S[j]) we have f(t_o + (j+1)*dt) := S[j].
       * 
       * Hence, S[j] = m_j * (t_o + j*dt) + n_j.
       * Due to the continuity requirement of the PwLF we have
       *        S[j+1] = m_j * (t_o + (j+1)*dt) + n_j
       * Hence, S[j+1] - S[j] = (j+1)*dt*m_j - j*dt*m_j = m_j*dt and thus
       *   
       *   m_j = (S[j+1] - S[j]) / dt;
       * 
       * Further, S[j] = m_j * (t_o + j*dt) + n_j
       * 
       *   n_j = S[j] - m_j * (t_o + j*dt)
       * 
       */
      for (unsigned j = 0; j < l_np; ++j) {
        real m = (i_samples[j+1] - i_samples[j]) / i_samplingInterval;
        o_pwLF->slopes[j] = m;
        o_pwLF->intercepts[j] = i_samples[j] - m * (i_onsetTime + j * i_samplingInterval);
      }
    }

    /** Returns integral_fromTime^toTime i_pwLF dt. */
    real computePwLFTimeIntegral(PiecewiseLinearFunction1D const& i_pwLF,
                                 double i_fromTime,
                                 double i_toTime);

    void addTimeIntegratedPointSourceNRF( real const i_mInvJInvPhisAtSources[tensor::mInvJInvPhisAtSources::size()],
                                          real const faultBasis[9],
                                          real A,
                                          std::array<real, 81> const &stiffnessTensor,
                                          std::array<PiecewiseLinearFunction1D, 3> const &slipRates,
                                          double i_fromTime,
                                          double i_toTime,
                                          real o_dofUpdate[tensor::Q::size()],
                                          unsigned int sourceNumber);
    /**
     * Point sources in SeisSol (\delta(x-x_s) * S(t)).
     * 
     * Computes Q_kl += int_a^b S(t) dt * phiAtSource[k] * momentTensor[l], where
     * Q_kl is a DOF. phiAtSource times momentTensor is to be understood
     * as outer product of two vectors (i.e. yields a rank-1 dof-update-matrix that shall
     * be scaled with the time integral of the source term).
     **/                                      
    void addTimeIntegratedPointSourceFSRM( real const i_mInvJInvPhisAtSources[tensor::mInvJInvPhisAtSources::size()],
                                           real const i_forceComponents[tensor::momentFSRM::size()],
                                           PiecewiseLinearFunction1D const& i_pwLF,
                                           double i_fromTime,
                                           double i_toTime,
                                           real o_dofUpdate[tensor::Q::size()],
                                           unsigned int sourceNumber);
  }
}

#endif
