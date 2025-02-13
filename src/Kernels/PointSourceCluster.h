// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
// SPDX-FileCopyrightText: 2023 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_POINTSOURCECLUSTER_H_
#define SEISSOL_SRC_KERNELS_POINTSOURCECLUSTER_H_

#include "Kernels/Precision.h"
#include "Numerical/Functions.h"
#include "Parallel/Runtime/Stream.h"
#include "SourceTerm/Typedefs.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>

namespace seissol::kernels {
class PointSourceCluster {
  public:
  virtual ~PointSourceCluster() = default;
  virtual void addTimeIntegratedPointSources(
      double from, double to, seissol::parallel::runtime::StreamRuntime& runtime) = 0;
  [[nodiscard]] virtual unsigned size() const = 0;
};

struct PointSourceClusterPair {
  std::shared_ptr<kernels::PointSourceCluster> host{nullptr};
  std::shared_ptr<kernels::PointSourceCluster> device{nullptr};
};

/**
 * @brief Integrate sample in time
 *
 * @param from Integration start time
 * @param to Integration end time
 * @param onsetTime Onset time of sample
 * @param samplingInterval Interval length (inverse of sampling rate)
 * @param sample Pointer to sample
 * @param sampleSize Size of the sample
 */
template <typename MathFunctions = seissol::functions::HostStdFunctions>
inline real computeSampleTimeIntegral(double from,
                                      double to,
                                      const double onsetTime,
                                      const double samplingInterval,
                                      const real* sample,
                                      std::size_t sampleSize) {
  const auto integrate = [&samplingInterval, &sample](std::size_t index, double tFrom, double tTo) {
    /* We have f(t) = S0 (t1 - t) / dt + S1 (t - t0) / dt, hence
     * int f(t) dt =  S0 (t1 t - 0.5 t^2) / dt + S1 (0.5 t^2 - t0 t) / dt + const, thus
     * int_tFrom^tTo f(t) dt = S0 (t1 (tTo - tFrom) - 0.5 (tTo^2 - tFrom^2)) / dt
     *                       + S1 (0.5 (tTo^2 - tFrom^2) - t0 (tTo - tFrom)) / dt
     */
    const auto t0 = index * samplingInterval;
    const auto t1 = t0 + samplingInterval;
    const auto s0 = sample[index];
    const auto s1 = sample[index + 1];
    const auto tdiff = tTo - tFrom;
    const auto tdiff2 = 0.5 * (tTo * tTo - tFrom * tFrom);
    return (s0 * (t1 * tdiff - tdiff2) + s1 * (tdiff2 - t0 * tdiff)) / samplingInterval;
  };

  if (sampleSize == 0) {
    return 0.0;
  }

  // Shift time such that t = 0 corresponds to onsetTime
  from -= onsetTime;
  to -= onsetTime;
  // Adjust integration interval to sample time interval
  // Sample is implicitly zero outside of sample time interval
  from = MathFunctions::max(from, 0.0);
  to = MathFunctions::min(to, (sampleSize - 1) * samplingInterval);

  // j_{from} := \argmax_j s.t. t_{from} >= j*dt = floor[t_{from} / dt]
  long fromIndex = MathFunctions::floor(from / samplingInterval);
  // j_{to}   := \argmin_j s.t. t_{to}   <= j*dt =  ceil[t_{to}   / dt]
  long toIndex = MathFunctions::ceil(to / samplingInterval);

  fromIndex = MathFunctions::max(0L, fromIndex);
  toIndex = MathFunctions::min(static_cast<long>(sampleSize) - 1, toIndex);
  // Return zero if there is no overlap between integration interval and sample time interval
  if (fromIndex >= toIndex) {
    return 0.0;
  }

  if (toIndex - fromIndex == 1L) {
    return integrate(fromIndex, from, to);
  }

  real integral = 0.0;
  integral += integrate(fromIndex, from, (fromIndex + 1) * samplingInterval);
  for (auto j = fromIndex + 1; j < toIndex - 1; ++j) {
    integral += 0.5 * samplingInterval * (sample[j] + sample[j + 1]);
  }
  integral += integrate(toIndex - 1, (toIndex - 1) * samplingInterval, to);
  return integral;
}

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_POINTSOURCECLUSTER_H_
