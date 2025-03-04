// SPDX-FileCopyrightText: 2015-2025 SeisSol Group
// SPDX-FileCopyrightText: 2023 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_POINTSOURCECLUSTER_H_
#define SEISSOL_SRC_KERNELS_POINTSOURCECLUSTER_H_

#include "Common/Marker.h"
#include "Kernels/Precision.h"
#include "Numerical/Functions.h"
#include "Parallel/Runtime/Stream.h"
#include "SourceTerm/Typedefs.h"

#include <Memory/MemoryAllocator.h>
#include <init.h>

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
  std::unique_ptr<kernels::PointSourceCluster> host{nullptr};
  std::unique_ptr<kernels::PointSourceCluster> device{nullptr};
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
SEISSOL_HOSTDEVICE inline real computeSampleTimeIntegral(double from,
                                                         double to,
                                                         const double onsetTime,
                                                         const double samplingInterval,
                                                         const real* __restrict sample,
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

// workaround for NVHPC (using constexpr arrays directly caused errors in 24.01)
constexpr std::size_t QSpan = init::Q::Stop[0] - init::Q::Start[0];
constexpr std::size_t MomentFsrmSpan = tensor::momentFSRM::Shape[0];
constexpr std::size_t MInvJInvPhisAtSourcesSpan = tensor::mInvJInvPhisAtSources::Shape[0];

SEISSOL_HOSTDEVICE inline void
    addTimeIntegratedPointSourceNRF(const memory::AlignedArray<real, 3>& __restrict slip,
                                    const real* __restrict mInvJInvPhisAtSources,
                                    const real* __restrict tensor,
                                    real a,
                                    const real* __restrict stiffnessTensor,
                                    double from,
                                    double to,
                                    real dofs[tensor::Q::size()]) {
  real rotatedSlip[3] = {real(0.0)};
  for (unsigned i = 0; i < 3; ++i) {
    for (unsigned j = 0; j < 3; ++j) {
      rotatedSlip[j] += tensor[j + i * 3] * slip[i];
    }
  }

  const auto mom = [&](unsigned p, unsigned q) {
    real m = 0.0;
    for (unsigned i = 0; i < 3; ++i) {
      for (unsigned j = 0; j < 3; ++j) {
        m += -a * stiffnessTensor[p + 3 * q + 9 * i + 27 * j] * rotatedSlip[i] * tensor[6 + j];
      }
    }
    return m;
  };

  const real moment[6] = {mom(0, 0), mom(1, 1), mom(2, 2), mom(0, 1), mom(1, 2), mom(0, 2)};
  for (unsigned t = 0; t < 6; ++t) {
    for (unsigned k = 0; k < MInvJInvPhisAtSourcesSpan; ++k) {
      dofs[k + t * QSpan] += mInvJInvPhisAtSources[k] * moment[t];
    }
  }
}

SEISSOL_HOSTDEVICE inline void pointSourceKernelNRF(
    int index,
    double from,
    double to,
    sourceterm::CellToPointSourcesMapping* __restrict mappingPtr,
    const seissol::memory::
        AlignedArray<real, tensor::mInvJInvPhisAtSources::size()>* __restrict mInvJInvPhisAtSources,
    const seissol::memory::AlignedArray<real,
                                        sourceterm::PointSources::TensorSize>* __restrict tensor,
    const real* __restrict a,
    const seissol::memory::AlignedArray<real, 81>* __restrict stiffnessTensor,
    const double* __restrict onsetTime,
    const double* __restrict samplingInterval,
    const memory::AlignedArray<const std::size_t* __restrict, 3> sampleOffsets,
    const memory::AlignedArray<const real* __restrict, 3> sample) {
  const unsigned startSource = mappingPtr[index].pointSourcesOffset;
  const unsigned endSource =
      mappingPtr[index].pointSourcesOffset + mappingPtr[index].numberOfPointSources;
  for (unsigned source = startSource; source < endSource; ++source) {
    memory::AlignedArray<real, 3> slip;
    for (int i = 0; i < 3; ++i) {
      auto o0 = sampleOffsets[i][source];
      auto o1 = sampleOffsets[i][source + 1];
      slip[i] = computeSampleTimeIntegral(
          from, to, onsetTime[source], samplingInterval[source], sample[i] + o0, o1 - o0);
    }

    addTimeIntegratedPointSourceNRF(slip,
                                    mInvJInvPhisAtSources[source].data(),
                                    tensor[source].data(),
                                    a[source],
                                    stiffnessTensor[source].data(),
                                    from,
                                    to,
                                    *mappingPtr[index].dofs);
  }
}

SEISSOL_HOSTDEVICE inline void
    addTimeIntegratedPointSourceFSRM(real slip,
                                     const real* __restrict mInvJInvPhisAtSources,
                                     const real* __restrict tensor,
                                     double from,
                                     double to,
                                     real* __restrict dofs) {
  for (unsigned p = 0; p < MomentFsrmSpan; ++p) {
    for (unsigned k = 0; k < MInvJInvPhisAtSourcesSpan; ++k) {
      dofs[k + p * QSpan] += slip * mInvJInvPhisAtSources[k] * tensor[p];
    }
  }
}

SEISSOL_HOSTDEVICE inline void pointSourceKernelFSRM(
    int index,
    double from,
    double to,
    sourceterm::CellToPointSourcesMapping* __restrict mappingPtr,
    const seissol::memory::
        AlignedArray<real, tensor::mInvJInvPhisAtSources::size()>* __restrict mInvJInvPhisAtSources,
    const seissol::memory::AlignedArray<real,
                                        sourceterm::PointSources::TensorSize>* __restrict tensor,
    const real* __restrict a,
    const seissol::memory::AlignedArray<real, 81>* __restrict stiffnessTensor,
    const double* __restrict onsetTime,
    const double* __restrict samplingInterval,
    const memory::AlignedArray<const std::size_t* __restrict, 3> sampleOffsets,
    const memory::AlignedArray<const real* __restrict, 3> sample) {
  const unsigned startSource = mappingPtr[index].pointSourcesOffset;
  const unsigned endSource =
      mappingPtr[index].pointSourcesOffset + mappingPtr[index].numberOfPointSources;
  for (unsigned source = startSource; source < endSource; ++source) {
    auto o0 = sampleOffsets[0][source];
    auto o1 = sampleOffsets[0][source + 1];
    const real slip = computeSampleTimeIntegral(
        from, to, onsetTime[source], samplingInterval[source], sample[0] + o0, o1 - o0);
    addTimeIntegratedPointSourceFSRM(slip,
                                     mInvJInvPhisAtSources[source].data(),
                                     tensor[source].data(),
                                     from,
                                     to,
                                     *mappingPtr[index].dofs);
  }
}

void pointSourceKernel(sourceterm::ClusterMapping& clusterMapping,
                       sourceterm::PointSources& sources,
                       double from,
                       double to,
                       seissol::parallel::runtime::StreamRuntime& runtime);

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_POINTSOURCECLUSTER_H_
