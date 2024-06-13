// Copyright (c) 2024 SeisSol Group
// Copyright (c) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#ifndef KERNELS_POINTSOURCECLUSTERONHOST_H_
#define KERNELS_POINTSOURCECLUSTERONHOST_H_

#include "PointSourceCluster.h"

#include "SourceTerm/typedefs.hpp"

namespace seissol::kernels {
class PointSourceClusterOnHost : public PointSourceCluster {
  public:
  PointSourceClusterOnHost(std::shared_ptr<sourceterm::ClusterMapping> mapping,
                           std::shared_ptr<sourceterm::PointSources> sources);
  void addTimeIntegratedPointSources(double from,
                                     double to,
                                     seissol::parallel::runtime::StreamRuntime& runtime) override;
  unsigned size() const override;

  private:
  void addTimeIntegratedPointSourceNRF(unsigned source,
                                       double from,
                                       double to,
                                       real dofs[tensor::Q::size()]);
  void addTimeIntegratedPointSourceFSRM(unsigned source,
                                        double from,
                                        double to,
                                        real dofs[tensor::Q::size()]);

  std::shared_ptr<sourceterm::ClusterMapping> clusterMapping_;
  std::shared_ptr<sourceterm::PointSources> sources_;
};
} // namespace seissol::kernels

#endif // KERNELS_POINTSOURCECLUSTERONHOST_H_
