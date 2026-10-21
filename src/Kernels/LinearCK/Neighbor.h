// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2014-2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#ifndef SEISSOL_SRC_KERNELS_LINEARCK_NEIGHBOR_H_
#define SEISSOL_SRC_KERNELS_LINEARCK_NEIGHBOR_H_

#include "Common/Constants.h"
#include "GeneratedCode/kernel.h"
#include "Kernels/Neighbor.h"
#include "Kernels/StarOperands.h"
#include "Kernels/WeightAdjustment.h"
#include "Memory/Descriptor/LTS.h"
#include "Model/OperatorLayout.h"
#include "Monitoring/Metric.h"
#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::kernels::solver::linearck {

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
  dynamicRupture::kernel::nodalFlux drKrnlPrototype_;
  WeightAdjustmentT weightAdjustment_;

#ifdef ACL_DEVICE
  kernel::gpu_faceCorrector deviceCorrectorKrnlPrototype_;
  dynamicRupture::kernel::gpu_nodalFlux deviceDrKrnlPrototype_;
  device::DeviceInstance& device_ = device::DeviceInstance::getInstance();
#endif
};

} // namespace seissol::kernels::solver::linearck

#endif // SEISSOL_SRC_KERNELS_LINEARCK_NEIGHBOR_H_
