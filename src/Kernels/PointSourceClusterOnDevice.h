// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
// SPDX-FileCopyrightText: 2023 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_POINTSOURCECLUSTERONDEVICE_H_
#define SEISSOL_SRC_KERNELS_POINTSOURCECLUSTERONDEVICE_H_

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
  std::shared_ptr<sourceterm::ClusterMapping> clusterMapping_;
  std::shared_ptr<sourceterm::PointSources> sources_;
};
} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_POINTSOURCECLUSTERONDEVICE_H_
