// SPDX-FileCopyrightText: 2024 SeisSol Group
// SPDX-FileCopyrightText: 2023 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_POINTSOURCECLUSTERONHOST_H_
#define SEISSOL_SRC_KERNELS_POINTSOURCECLUSTERONHOST_H_

#include "Common/Real.h"
#include "GeneratedCode/tensor.h"
#include "PointSourceCluster.h"
#include "SourceTerm/Typedefs.h"

namespace seissol::kernels {
/// The point sources of a time cluster whose cells compute in the configuration `Cfg`, on the host.
template <typename Cfg>
class PointSourceClusterOnHost : public PointSourceCluster {
  public:
  using real = Real<Cfg>; // NOLINT(readability-identifier-naming)

  PointSourceClusterOnHost(std::shared_ptr<sourceterm::ClusterMapping> mapping,
                           std::shared_ptr<sourceterm::PointSources<Cfg>> sources);
  void addTimeIntegratedPointSources(double from,
                                     double to,
                                     seissol::parallel::runtime::StreamRuntime& runtime) override;
  [[nodiscard]] std::size_t size() const override;

  private:
  void addTimeIntegratedPointSource(std::size_t source,
                                    double from,
                                    double to,
                                    real dofs[tensor::Q<Cfg>::size()]);

  std::shared_ptr<sourceterm::ClusterMapping> clusterMapping_;
  std::shared_ptr<sourceterm::PointSources<Cfg>> sources_;
};
} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_POINTSOURCECLUSTERONHOST_H_
