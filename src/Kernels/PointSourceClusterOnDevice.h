// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#ifndef KERNELS_POINTSOURCECLUSTERONDEVICE_H_
#define KERNELS_POINTSOURCECLUSTERONDEVICE_H_

#include "PointSourceCluster.h"
#include <SourceTerm/typedefs.hpp>

#include <array>

namespace seissol::kernels {
class PointSourceClusterOnDevice : public PointSourceCluster {
  public:
  PointSourceClusterOnDevice(sourceterm::ClusterMapping mapping, sourceterm::PointSources sources);
  void addTimeIntegratedPointSources(double from, double to) override;
  unsigned size() const override;

  private:
  static void addTimeIntegratedPointSourceNRF(std::array<real, 3> const& slip,
                                              real* mInvJInvPhisAtSources,
                                              real* tensor,
                                              real A,
                                              real* stiffnessTensor,
                                              double from,
                                              double to,
                                              real* dofs);
  static void addTimeIntegratedPointSourceFSRM(
      real slip, real* mInvJInvPhisAtSources, real* tensor, double from, double to, real* dofs);

  sourceterm::ClusterMapping clusterMapping_;
  sourceterm::PointSources sources_;
};
} // namespace seissol::kernels

#endif // KERNELS_POINTSOURCECLUSTERONDEVICE_H_
