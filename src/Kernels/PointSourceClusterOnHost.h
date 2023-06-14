// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#ifndef KERNELS_POINTSOURCECLUSTERONHOST_H_
#define KERNELS_POINTSOURCECLUSTERONHOST_H_

#include "PointSourceCluster.h"

#include <SourceTerm/typedefs.hpp>

namespace seissol::kernels {
class PointSourceClusterOnHost : public PointSourceCluster {
  public:
  PointSourceClusterOnHost(sourceterm::ClusterMapping mapping, sourceterm::PointSources sources);
  void addTimeIntegratedPointSources(double from, double to) override;

  private:
  void addTimeIntegratedPointSourceNRF(unsigned source,
                                       double from,
                                       double to,
                                       real dofs[tensor::Q::size()]);
  void addTimeIntegratedPointSourceFSRM(unsigned source,
                                        double from,
                                        double to,
                                        real dofs[tensor::Q::size()]);

  sourceterm::ClusterMapping clusterMapping_;
  sourceterm::PointSources sources_;
};
} // namespace seissol::kernels

#endif // KERNELS_POINTSOURCECLUSTERONHOST_H_
