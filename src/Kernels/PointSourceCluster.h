// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#ifndef KERNELS_POINTSOURCECLUSTER_H_
#define KERNELS_POINTSOURCECLUSTER_H_

#include "SourceTerm/typedefs.hpp"
#include "Kernels/precision.hpp"

#include <cstdint>

namespace seissol::kernels {
class PointSourceCluster {
  public:
  virtual ~PointSourceCluster() = default;
  virtual void addTimeIntegratedPointSources(double from, double to) = 0;
  virtual unsigned size() const = 0;
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
real computeSampleTimeIntegral(double from,
                               double to,
                               double const onsetTime,
                               double const samplingInterval,
                               real* sample,
                               std::size_t sampleSize);

} // namespace seissol::kernels

#endif // KERNELS_POINTSOURCECLUSTER_H_
