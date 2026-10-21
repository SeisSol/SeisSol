// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_NEIGHBOR_H_
#define SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_NEIGHBOR_H_

#include "GeneratedCode/kernel.h"
#include "Kernels/Neighbor.h"
#include "Kernels/WeightAdjustment.h"
#include "Model/OperatorLayout.h"
#include "Monitoring/Metric.h"

namespace seissol::kernels::solver::linearckanelastic {
class Neighbor : public NeighborKernel {
  public:
  void setGlobalData(const CompoundGlobalData& global) override;

  void computeNeighborsIntegral(
      LTS::Ref& data, const std::array<const real*, Cell::NumFaces>& neighborTraces) override;

  void computeBatchedNeighborsIntegral(recording::ConditionalPointersToRealsTable& table,
                                       seissol::parallel::runtime::StreamRuntime& runtime) override;

  [[nodiscard]] std::pair<PerformanceEstimate, PerformanceEstimate>
      metrics(const std::array<FaceType, Cell::NumFaces>& faceTypes,
              const std::array<CellDRMapping, Cell::NumFaces>& cellDrMapping) const override;

  [[nodiscard]] PerformanceEstimate weightAdjustMetrics() const override {
    return decltype(weightAdjustment_)::metrics();
  }

  protected:
  kernel::faceCorrector correctorKrnlPrototype_;
  kernel::neighbor nKrnlPrototype_;
  dynamicRupture::kernel::nodalFlux drKrnlPrototype_;
  WeightAdjustmentT weightAdjustment_;

#ifdef ACL_DEVICE
  kernel::gpu_faceCorrector deviceCorrectorKrnlPrototype_;
  kernel::gpu_neighbor deviceNKrnlPrototype_;
  dynamicRupture::kernel::gpu_nodalFlux deviceDrKrnlPrototype_;
#endif
};
} // namespace seissol::kernels::solver::linearckanelastic

#endif // SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_NEIGHBOR_H_
