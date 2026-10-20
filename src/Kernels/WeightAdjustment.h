// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_WEIGHTADJUSTMENT_H_
#define SEISSOL_SRC_KERNELS_WEIGHTADJUSTMENT_H_

#include "Initializer/Typedefs.h"
#include "Kernels/StarOperands.h"
#include "Memory/Descriptor/LTS.h"
#include "Model/OperatorLayout.h"
#include "Monitoring/Metric.h"

namespace seissol::kernels {

/// Adds the inverse weight of the medium applied to the flux buffer of a weight-adjusted cell to
/// the field its faces lift into otherwise (see WeightAdjustable). A class of its own, so that a
/// build without the kernel never instantiates what names its operands.
template <bool Enabled, typename KernelT>
class WeightAdjustment {
  public:
  void setGlobalData(const CompoundGlobalData& global) { prototype_.bindGlobals(*global.onHost); }

  /// `bindTarget` hands the kernel that field: the degrees of freedom, or the extended field of a
  /// solver that keeps the memory variables apart.
  template <typename BindT>
  void apply(LTS::Ref& data, const BindT& bindTarget) const {
    KernelT krnl = prototype_;
    bindTarget(krnl);
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
  template <typename BindT>
  void apply(LTS::Ref& /*data*/, const BindT& /*bindTarget*/) const {}
  [[nodiscard]] static PerformanceEstimate metrics() { return {}; }
};

/// The one of this build.
using WeightAdjustmentT = WeightAdjustment<WeightAdjustable, kernel::weightAdjust>;

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_WEIGHTADJUSTMENT_H_
