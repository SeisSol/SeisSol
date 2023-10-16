// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#ifndef KERNELS_POINTSOURCECLUSTER_H_
#define KERNELS_POINTSOURCECLUSTER_H_

#include "SourceTerm/typedefs.hpp"

namespace seissol::kernels {
class PointSourceCluster {
  public:
  virtual ~PointSourceCluster() = default;
  virtual void addTimeIntegratedPointSources(double from, double to) = 0;
  virtual unsigned size() const = 0;
};
} // namespace seissol::kernels

#endif // KERNELS_POINTSOURCECLUSTER_H_
