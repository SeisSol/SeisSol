// Copyright (c) 2015-2020 SeisSol Group
// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SOURCETERM_PIECEWISELINEARFUNCTION1D_H
#define SOURCETERM_PIECEWISELINEARFUNCTION1D_H

#include <Kernels/precision.hpp>

#include <memory>
#include <vector>
#include <algorithm>

namespace seissol::sourceterm {

/** A piecewise linear function.
 *
 *  Say t \in I_j, then
 *    f(t) = m_j * t + n_j,
 *  where I_j is the half-open interval [t_o + j*dt, t_o + (j+1)*dt).
 *  j runs through 0,...,n-1.
 **/
template <class Allocator>
struct PiecewiseLinearFunction1D {
  using RealAllocator = typename std::allocator_traits<Allocator>::template rebind_alloc<real>;

  /** slopes[i] = m_i */
  std::vector<real, RealAllocator> slopes;

  /** intercepts[i] = n_i */
  std::vector<real, RealAllocator> intercepts;

  /** onsetTime = t_o */
  real onsetTime;

  /** samplingInterval = dt */
  real samplingInterval;

  PiecewiseLinearFunction1D(Allocator const& alloc)
      : slopes(RealAllocator(alloc)), intercepts(RealAllocator(alloc)) {}
  template <typename real_from>
  PiecewiseLinearFunction1D(real_from const* i_samples,
                            unsigned i_numberOfSamples,
                            real i_onsetTime,
                            real i_samplingInterval,
                            Allocator const& alloc)
      : slopes(RealAllocator(alloc)), intercepts(RealAllocator(alloc)) {
    if (i_numberOfSamples > 0) {
      unsigned l_np = i_numberOfSamples - 1;

      slopes.resize(l_np);
      intercepts.resize(l_np);
      onsetTime = i_onsetTime;
      samplingInterval = i_samplingInterval;

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
        real m = (i_samples[j + 1] - i_samples[j]) / i_samplingInterval;
        slopes[j] = m;
        intercepts[j] = i_samples[j] - m * (i_onsetTime + j * i_samplingInterval);
      }
    }
  }

  /** Returns integral_fromTime^toTime i_pwLF dt. */
  #pragma omp declare target
  real timeIntegral(double i_fromTime, double i_toTime) const {
    real l_integral;
    // j_{from} := \argmax_j s.t. t_{from} >= t_{onset} + j*dt   =   floor[(t_{from} - t_{onset}) /
    // dt]
    int l_fromIndex = (i_fromTime - onsetTime) / samplingInterval;
    // j_{to}   := \argmin_j s.t. t_{to}   >= t_{onset} + j*dt   =   floor[(t_{to} - t_{onset}) /
    // dt]
    int l_toIndex = (i_toTime - onsetTime) / samplingInterval;

    l_fromIndex = std::max(0, l_fromIndex);
    l_toIndex = std::min(static_cast<int>(slopes.size()) - 1, l_toIndex);

    /* The indefinite integral of the j-th linear function is
     * int m_j * t + n_j dt = 1 / 2 * m_j * t^2 + n_j * t
     */
    real l_time = onsetTime + l_fromIndex * samplingInterval;
    l_integral = 0.0;
    for (int j = l_fromIndex; j <= l_toIndex; ++j) {
      real tFrom = std::max((real)i_fromTime, l_time);
      l_time += samplingInterval;
      real tTo = std::min((real)i_toTime, l_time);
      l_integral += 0.5 * slopes[j] * (tTo * tTo - tFrom * tFrom) + intercepts[j] * (tTo - tFrom);
    }

    return l_integral;
  }
  #pragma omp end declare target
};

} // namespace seissol::sourceterm

#endif // SOURCETERM_PIECEWISELINEARFUNCTION1D_H
