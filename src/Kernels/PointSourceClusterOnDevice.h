// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#ifndef KERNELS_POINTSOURCECLUSTERONDEVICE_H_
#define KERNELS_POINTSOURCECLUSTERONDEVICE_H_

#include "PointSourceCluster.h"

#include <SourceTerm/typedefs.hpp>

namespace seissol::kernels {
class PointSourceClusterOnDevice : public PointSourceCluster {
  public:
  PointSourceClusterOnDevice(sourceterm::ClusterMapping mapping, sourceterm::PointSources sources);
  void addTimeIntegratedPointSources(double from, double to) override;

  private:
  static void addTimeIntegratedPointSourceNRF(
      std::array<sourceterm::PiecewiseLinearFunction1D<sourceterm::AllocatorT> const*, 3> slipRates,
      real* mInvJInvPhisAtSources,
      real* tensor,
      real A,
      real* stiffnessTensor,
      double from,
      double to,
      real* dofs);
  static void addTimeIntegratedPointSourceFSRM(
      sourceterm::PiecewiseLinearFunction1D<sourceterm::AllocatorT> const* slipRate0,
      real* mInvJInvPhisAtSources,
      real* tensor,
      double from,
      double to,
      real* dofs);

  sourceterm::ClusterMapping clusterMapping_;
  sourceterm::PointSources sources_;
};
} // namespace seissol::kernels

#endif // KERNELS_POINTSOURCECLUSTERONDEVICE_H_
