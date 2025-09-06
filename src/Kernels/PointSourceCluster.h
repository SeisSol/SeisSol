// SPDX-FileCopyrightText: 2015 SeisSol Group
// SPDX-FileCopyrightText: 2023 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_POINTSOURCECLUSTER_H_
#define SEISSOL_SRC_KERNELS_POINTSOURCECLUSTER_H_

#include "Common/Marker.h"
#include "GeneratedCode/init.h"
#include "Kernels/Precision.h"
#include "Parallel/Runtime/Stream.h"
#include "SourceTerm/Typedefs.h"
#include <Equations/Datastructures.h>
#include <Memory/MemoryAllocator.h>
#include <Solver/MultipleSimulations.h>

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
  [[nodiscard]] virtual std::size_t size() const = 0;
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
template <typename RealT>
SEISSOL_HOSTDEVICE inline RealT computeSampleTimeIntegral(double from,
                                                          double to,
                                                          const double onsetTime,
                                                          const double samplingInterval,
                                                          const RealT* __restrict sample,
                                                          ssize_t sampleSize) {
  const auto integrate = [&samplingInterval, &sample](ssize_t index, double tFrom, double tTo) {
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
  from = std::max(from, 0.0);
  to = std::min(to, (sampleSize - 1) * samplingInterval);

  // j_{from} := \argmax_j s.t. t_{from} >= j*dt = floor[t_{from} / dt]
  ssize_t fromIndex = std::floor(from / samplingInterval);
  // j_{to}   := \argmin_j s.t. t_{to}   <= j*dt =  ceil[t_{to}   / dt]
  ssize_t toIndex = std::ceil(to / samplingInterval);

  fromIndex = std::max(static_cast<ssize_t>(0), fromIndex);
  toIndex = std::min(sampleSize - 1, toIndex);
  // Return zero if there is no overlap between integration interval and sample time interval
  if (fromIndex >= toIndex) {
    return 0.0;
  }

  if (toIndex - fromIndex == static_cast<ssize_t>(1)) {
    return integrate(fromIndex, from, to);
  }

  RealT integral = 0.0;
  integral += integrate(fromIndex, from, (fromIndex + 1) * samplingInterval);
  for (auto j = fromIndex + 1; j < toIndex - 1; ++j) {
    integral += 0.5 * samplingInterval * (sample[j] + sample[j + 1]);
  }
  integral += integrate(toIndex - 1, (toIndex - 1) * samplingInterval, to);
  return integral;
}

// workaround for NVHPC (using constexpr arrays directly caused errors in 24.01)
template <typename Cfg>
constexpr std::size_t QSpan =
    init::Q<Cfg>::Stop[multisim::BasisDim<Cfg>] - init::Q<Cfg>::Start[multisim::BasisDim<Cfg>];
template <typename Cfg>
constexpr std::size_t QMultiSpan = init::Q<Cfg>::Stop[0] - init::Q<Cfg>::Start[0];
template <typename Cfg>
constexpr std::size_t MomentFsrmSpan = tensor::update<Cfg>::Shape[0];
template <typename Cfg>
constexpr std::size_t MInvJInvPhisAtSourcesSpan = tensor::mInvJInvPhisAtSources<Cfg>::Shape[0];

template <typename Cfg>
constexpr std::size_t Quantities = MomentFsrmSpan<Cfg>;

template <typename Cfg>
SEISSOL_HOSTDEVICE constexpr auto&
    dofsAccessor(Real<Cfg>* __restrict dofs, std::uint32_t k, std::uint32_t t, std::uint32_t f) {
  if constexpr (seissol::multisim::MultisimEnabled<Cfg>) {
    return dofs[(k + t * QSpan<Cfg>)*QMultiSpan<Cfg> + f];
  } else {
    return dofs[k + t * QSpan<Cfg>];
  }
}

template <typename Cfg, std::uint32_t Block>
SEISSOL_HOSTDEVICE inline void pointSourceKernelDevice(
    std::uint32_t thread,
    std::size_t index,
    double from,
    double to,
    sourceterm::CellToPointSourcesMapping* __restrict mappingPtr,
    const seissol::memory::AlignedArray<
        Real<Cfg>,
        tensor::mInvJInvPhisAtSources<Cfg>::size()>* __restrict mInvJInvPhisAtSources,
    const std::uint32_t* __restrict simulationIndex,
    const Real<Cfg>* __restrict tensor,
    const double* __restrict onsetTime,
    const double* __restrict samplingInterval,
    const std::size_t* __restrict sampleRange,
    const std::size_t* __restrict sampleOffsets,
    const Real<Cfg>* __restrict sample) {
  const auto startSource = mappingPtr[index].pointSourcesOffset;
  const auto endSource =
      mappingPtr[index].pointSourcesOffset + mappingPtr[index].numberOfPointSources;

  auto* __restrict dofs = reinterpret_cast<Real<Cfg>*>(mappingPtr[index].dofs);
  for (std::size_t source = startSource; source < endSource; ++source) {
    const auto base = sampleRange[source];
    const std::uint32_t localSamples = sampleRange[source + 1] - base;

    const auto* __restrict tensorLocal = tensor + base * Quantities<Cfg>;

    std::array<Real<Cfg>, Quantities<Cfg>> update{};

#pragma unroll 3
    for (std::uint32_t i = 0; i < localSamples; ++i) {
      const auto o0 = sampleOffsets[i + base];
      const auto o1 = sampleOffsets[i + base + 1];
      const auto slip = computeSampleTimeIntegral(
          from, to, onsetTime[source], samplingInterval[source], sample + o0, o1 - o0);

#pragma unroll
      for (std::uint32_t t = 0; t < Quantities<Cfg>; ++t) {
        update[t] += slip * tensorLocal[t + i * Quantities<Cfg>];
      }
    }

#pragma unroll
    for (std::uint32_t t = 0; t < Quantities<Cfg>; ++t) {
      for (std::uint32_t k = thread; k < MInvJInvPhisAtSourcesSpan<Cfg>; k += Block) {
        dofsAccessor<Cfg>(dofs, k, t, simulationIndex[source]) +=
            mInvJInvPhisAtSources[source][k] * update[t];
      }
    }
  }
}

template <typename Cfg>
void pointSourceKernel(sourceterm::ClusterMapping& clusterMapping,
                       sourceterm::PointSources<Cfg>& sources,
                       double from,
                       double to,
                       seissol::parallel::runtime::StreamRuntime& runtime);

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_POINTSOURCECLUSTER_H_
