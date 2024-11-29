// Copyright (c) 2024 SeisSol Group
// Copyright (c) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#ifndef KERNELS_POINTSOURCECLUSTERONDEVICE_H_
#define KERNELS_POINTSOURCECLUSTERONDEVICE_H_

#include "PointSourceCluster.h"
#include "SourceTerm/Typedefs.h"

#include <array>

namespace seissol::kernels {
class PointSourceClusterOnDevice : public PointSourceCluster {
  public:
  PointSourceClusterOnDevice(std::shared_ptr<sourceterm::ClusterMapping> mapping,
                             std::shared_ptr<sourceterm::PointSources> sources);
  void addTimeIntegratedPointSources(double from,
                                     double to,
                                     seissol::parallel::runtime::StreamRuntime& runtime) override;
  unsigned size() const override;

  private:
  #pragma omp declare target
  static void addTimeIntegratedPointSourceNRF(const std::array<real, 3>& slip,
                                              const real* mInvJInvPhisAtSources,
                                              const real* tensor,
                                              real a,
                                              const real* stiffnessTensor,
                                              double from,
                                              double to,
                                              real* dofs);
  static void addTimeIntegratedPointSourceFSRM(real slip,
                                               const real* mInvJInvPhisAtSources,
                                               const real* tensor,
                                               double from,
                                               double to,
                                               real* dofs);
  #pragma omp end declare target

  std::shared_ptr<sourceterm::ClusterMapping> clusterMapping_;
  std::shared_ptr<sourceterm::PointSources> sources_;
};
} // namespace seissol::kernels

#endif // KERNELS_POINTSOURCECLUSTERONDEVICE_H_
