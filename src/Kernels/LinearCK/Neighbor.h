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
#include "Memory/Descriptor/LTS.h"
#include "Model/OperatorLayout.h"
#include "Monitoring/Metric.h"
#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::kernels::solver::linearck {

/// Adds the inverse weight of the medium applied to the flux buffer of a weight-adjusted cell to
/// its field (see WeightAdjustable). A class of its own, so that a build without the kernel never
/// instantiates what names its operands.
template <bool Enabled, typename KernelT>
class WeightAdjustment {
  public:
  void setGlobalData(const CompoundGlobalData& global) { prototype_.bindGlobals(*global.onHost); }

  void apply(LTS::Ref& data) const {
    KernelT krnl = prototype_;
    krnl.Q = data.get<LTS::Dofs>();
    krnl.fluxBuffer = data.get<LTS::FluxBuffer>();
    kernels::bindWeightOperands(krnl, data.get<LTS::LocalIntegration>());
    krnl.execute();
  }

  [[nodiscard]] static PerformanceEstimate metrics() {
    return PerformanceEstimate::fromKernel<KernelT>();
  }

  private:
  KernelT prototype_;
};

template <typename KernelT>
class WeightAdjustment<false, KernelT> {
  public:
  void setGlobalData(const CompoundGlobalData& /*global*/) {}
  void apply(LTS::Ref& /*data*/) const {}
  [[nodiscard]] static PerformanceEstimate metrics() { return {}; }
};

class Neighbor : public NeighborKernel {
  public:
  void setGlobalData(const CompoundGlobalData& global) override;

  void computeNeighborsIntegral(
      LTS::Ref& data,
      const std::array<real*, Cell::NumFaces>& timeIntegrated,
      const std::array<real*, Cell::NumFaces>& faceNeighborsPrefetch) override;

  void computeBatchedNeighborsIntegral(recording::ConditionalPointersToRealsTable& table,
                                       seissol::parallel::runtime::StreamRuntime& runtime) override;

  [[nodiscard]] std::pair<PerformanceEstimate, PerformanceEstimate>
      metrics(const std::array<FaceType, Cell::NumFaces>& faceTypes,
              const std::array<std::array<uint8_t, 2>, Cell::NumFaces>& neighboringIndices,
              const std::array<CellDRMapping, Cell::NumFaces>& cellDrMapping) const override;

  [[nodiscard]] PerformanceEstimate weightAdjustMetrics() const override {
    return decltype(weightAdjustment_)::metrics();
  }

  protected:
  kernel::neighboringFlux nfKrnlPrototype_;
  dynamicRupture::kernel::nodalFlux drKrnlPrototype_;
  WeightAdjustment<WeightAdjustable, kernel::weightAdjust> weightAdjustment_;

#ifdef ACL_DEVICE
  kernel::gpu_neighboringFlux deviceNfKrnlPrototype_;
  dynamicRupture::kernel::gpu_nodalFlux deviceDrKrnlPrototype_;
  device::DeviceInstance& device_ = device::DeviceInstance::getInstance();
#endif
};

} // namespace seissol::kernels::solver::linearck

#endif // SEISSOL_SRC_KERNELS_LINEARCK_NEIGHBOR_H_
